# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
# -*- coding: utf-8 -*-

"""BAM file writing with filtered output generation.

Handles efficient BAM file writing with reference filtering, header generation,
and optional PMD tag addition for filtered alignment output.
"""

from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t, uint16_t
from libc.stdlib cimport malloc, free, calloc, realloc, qsort
from libc.string cimport memcpy, memset
from libc.stdint cimport intptr_t

# Import basic types from processor pxd
from bam_filter.processor cimport MemoryPool, samFile
from bam_filter.processor_mapping cimport ReferenceMapping, create_filtered_header_efficient
from bam_filter import logging as bf_logging

LOG_TAG = "BAM-WRITER"

# HTSlib declarations (full enough for writer)
cdef extern from "htslib/sam.h":
	ctypedef struct bam1_core_t:
		int32_t tid
		int32_t pos
		uint32_t bin
		uint8_t qual
		uint8_t l_qname
		uint16_t flag
		uint32_t n_cigar
		int32_t l_qseq
		int32_t mtid
		int32_t mpos
		int32_t isize

	ctypedef struct bam1_t:
		bam1_core_t core
		uint64_t id
		uint8_t* data
		int l_data
		uint32_t m_data

	ctypedef struct sam_hdr_t:
		int n_targets

	ctypedef struct hts_idx_t
	ctypedef struct hts_itr_t

	# we import samFile from bam_filter.processor to avoid duplicate typedefs

	ctypedef struct BGZF

	samFile* hts_open(const char* fn, const char* mode)
	int hts_close(samFile* fp)
	sam_hdr_t* sam_hdr_read(samFile* fp)
	void sam_hdr_destroy(sam_hdr_t* h)
	const char* sam_hdr_str(sam_hdr_t* header)
	const char* sam_hdr_tid2name(sam_hdr_t* header, int tid)
	int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid)
	sam_hdr_t* sam_hdr_parse(size_t l_text, const char* text)
	int32_t sam_hdr_nref(sam_hdr_t* header)
	int32_t bam_endpos(bam1_t* b) nogil

	int sam_read1(samFile* fp, sam_hdr_t* h, bam1_t* b) nogil
	int sam_write1(samFile* fp, const sam_hdr_t* h, const bam1_t* b) nogil
	int sam_hdr_write(samFile* fp, const sam_hdr_t* h)

	bam1_t* bam_init1() nogil
	void bam_destroy1(bam1_t* b) nogil

	hts_idx_t* sam_index_load(samFile* fp, const char* fn)
	void hts_idx_destroy(hts_idx_t* idx)
	hts_itr_t* sam_itr_queryi(const hts_idx_t* idx, int tid, int beg, int end)
	int sam_itr_next(samFile* fp, hts_itr_t* iter, bam1_t* b)
	void hts_itr_destroy(hts_itr_t* iter)
	int32_t sam_hdr_name2tid(sam_hdr_t* header, const char* name)
	int hts_set_threads(samFile* fp, int n)

from bam_filter.processor cimport bam1_t

# The actual implementation was moved here from processor.pyx. For brevity we re-implement only
# the large functions used by external callers: create_lookup_table,
# destroy_lookup_table, create_write_batch, destroy_write_batch, copy_bam_record,
# write_batch_to_bam, write_filtered_bam, write_bam_with_filtered_header.

cdef struct CompactAlignment:
	uint32_t original_position
	uint32_t pool_index
	uint32_t reference_id


cdef struct LookupTable:
	CompactAlignment* alignments
	uint64_t* ref_starts
	uint32_t* ref_counts
	uint64_t total_count
	uint32_t num_refs


cdef struct WriteBatch:
	bam1_t** records
	uint64_t* pool_indices
	uint32_t count
	uint32_t capacity
	uint32_t reference_id


cdef int _cmp_compact(const void* a, const void* b) noexcept nogil:
	"""Comparator for CompactAlignment by original_position (qsort-compatible)."""
	cdef CompactAlignment* A = <CompactAlignment*>a
	cdef CompactAlignment* B = <CompactAlignment*>b
	if A.original_position < B.original_position:
		return -1
	elif A.original_position > B.original_position:
		return 1
	else:
		return 0


cdef LookupTable* create_lookup_table(MemoryPool* pool, ReferenceMapping* mapping) except NULL nogil:
	"""Create lookup table mapping original BAM positions to filtered pool indices.

	Builds sorted index of alignments grouped by original reference ID for
	efficient BAM writing. Each reference's alignments are sorted by original
	position for sequential access during BAM iteration.

	Parameters
	----------
	pool : MemoryPool*
		Memory pool with filtered alignments
	mapping : ReferenceMapping*
		Mapping between original and retained reference IDs

	Returns
	-------
	LookupTable*
		Lookup table structure, or NULL on error
	"""
	cdef LookupTable* table = <LookupTable*>malloc(sizeof(LookupTable))
	if not table:
		return NULL

	table.total_count = pool.alignment_count
	table.num_refs = mapping.n_original_refs

	table.alignments = <CompactAlignment*>malloc(table.total_count * sizeof(CompactAlignment))
	table.ref_starts = <uint64_t*>malloc(table.num_refs * sizeof(uint64_t))
	table.ref_counts = <uint32_t*>calloc(table.num_refs, sizeof(uint32_t))

	if not table.alignments or not table.ref_starts or not table.ref_counts:
		if table.alignments: free(table.alignments)
		if table.ref_starts: free(table.ref_starts)
		if table.ref_counts: free(table.ref_counts)
		free(table)
		return NULL

	cdef uint64_t i
	cdef uint32_t compact_ref, original_ref

	for i in range(table.total_count):
		compact_ref = pool.alignments[i].reference_index
		if compact_ref < mapping.n_retained_refs:
			original_ref = mapping.new_to_old_tid[compact_ref]
			if original_ref < table.num_refs:
				table.ref_counts[original_ref] += 1
			else:
				destroy_lookup_table(table)
				return NULL
		else:
			destroy_lookup_table(table)
			return NULL

	table.ref_starts[0] = 0
	for i in range(1, table.num_refs):
		table.ref_starts[i] = table.ref_starts[i-1] + table.ref_counts[i-1]

	memset(table.ref_counts, 0, table.num_refs * sizeof(uint32_t))

	cdef uint64_t write_pos
	for i in range(table.total_count):
		compact_ref = pool.alignments[i].reference_index
		if compact_ref < mapping.n_retained_refs:
			original_ref = mapping.new_to_old_tid[compact_ref]
			if original_ref < table.num_refs:
				write_pos = table.ref_starts[original_ref] + table.ref_counts[original_ref]
				table.alignments[write_pos].original_position = pool.alignments[i].alignment_position
				table.alignments[write_pos].pool_index = i
				table.alignments[write_pos].reference_id = original_ref
				table.ref_counts[original_ref] += 1

	cdef uint64_t ref_start
	cdef uint32_t ref_count
	for i in range(table.num_refs):
		if table.ref_counts[i] > 1:
			ref_start = table.ref_starts[i]
			ref_count = table.ref_counts[i]
			qsort(<void*> &table.alignments[ref_start], ref_count, sizeof(CompactAlignment), _cmp_compact)

	return table


cdef void destroy_lookup_table(LookupTable* table) noexcept nogil:
	"""Free all memory associated with lookup table.

	Parameters
	----------
	table : LookupTable*
		Lookup table to destroy (safe to pass NULL)
	"""
	if not table:
		return
	if table.alignments: free(table.alignments)
	if table.ref_starts: free(table.ref_starts)
	if table.ref_counts: free(table.ref_counts)
	free(table)


cdef WriteBatch* create_write_batch(uint32_t capacity) except NULL nogil:
	"""Create write batch for buffered BAM output.

	Parameters
	----------
	capacity : uint32_t
		Maximum number of records in batch

	Returns
	-------
	WriteBatch*
		Batch structure with preallocated BAM records, or NULL on error
	"""
	cdef WriteBatch* batch = <WriteBatch*>malloc(sizeof(WriteBatch))
	if not batch:
		return NULL

	batch.capacity = capacity
	batch.count = 0
	batch.reference_id = 0

	batch.records = <bam1_t**>malloc(capacity * sizeof(bam1_t*))
	batch.pool_indices = <uint64_t*>malloc(capacity * sizeof(uint64_t))

	if not batch.records or not batch.pool_indices:
		if batch.records: free(batch.records)
		if batch.pool_indices: free(batch.pool_indices)
		free(batch)
		return NULL

	cdef uint32_t i, j
	for i in range(capacity):
		batch.records[i] = bam_init1()
		if not batch.records[i]:
			for j in range(i):
				if batch.records[j]:
					bam_destroy1(batch.records[j])
			free(batch.records)
			free(batch.pool_indices)
			free(batch)
			return NULL

	return batch


cdef void destroy_write_batch(WriteBatch* batch) noexcept nogil:
	"""Free all memory associated with write batch.

	Parameters
	----------
	batch : WriteBatch*
		Write batch to destroy (safe to pass NULL)
	"""
	if not batch:
		return
	if batch.records:
		for i in range(batch.capacity):
			if batch.records[i]:
				bam_destroy1(batch.records[i])
		free(batch.records)
	if batch.pool_indices:
		free(batch.pool_indices)
	free(batch)


cdef int copy_bam_record(bam1_t* src, bam1_t* dst, int extra_bytes) except -1 nogil:
	"""Copy BAM record with space for additional tags.

	Parameters
	----------
	src : bam1_t*
		Source BAM record
	dst : bam1_t*
		Destination BAM record (reallocated if needed)
	extra_bytes : int
		Additional bytes to allocate for tags

	Returns
	-------
	int
		0 on success, -1 on allocation failure
	"""
	cdef uint32_t needed_size = src.l_data + extra_bytes
	cdef uint32_t slack = (needed_size >> 3) + 64
	cdef uint32_t alloc_size = needed_size + slack

	if dst.m_data < needed_size:
		dst.data = <uint8_t*>realloc(dst.data, alloc_size)
		if not dst.data:
			return -1
		dst.m_data = alloc_size

	dst.core = src.core
	dst.id = src.id
	memcpy(dst.data, src.data, src.l_data)
	dst.l_data = src.l_data + extra_bytes

	return 0


cdef int write_batch_to_bam(samFile* out_bam, sam_hdr_t* header, WriteBatch* batch) except -1 nogil:
	"""Write all records in batch to BAM file.

	Parameters
	----------
	out_bam : samFile*
		Output BAM file handle
	header : sam_hdr_t*
		BAM header
	batch : WriteBatch*
		Write batch containing records

	Returns
	-------
	int
		0 on success, -1 on write error
	"""
	cdef uint32_t i
	cdef bam1_t** recs = batch.records
	cdef uint32_t n = batch.count
	for i in range(n):
		if sam_write1(out_bam, header, recs[i]) < 0:
			return -1
	return 0


cdef int write_filtered_bam(MemoryPool* pool,
							 const char* input_bam_path,
							 const char* output_bam_path,
							 sam_hdr_t* header,
							 ReferenceMapping* mapping,
							 int num_threads) except -1:
	"""Write filtered BAM file with updated reference IDs and additional tags.

	Main BAM writing function. Reads original BAM sequentially per reference,
	matches alignments via position counter to filtered pool, updates reference
	IDs, adds ZP/ZS/PM tags, and writes to output with filtered header.

	Parameters
	----------
	pool : MemoryPool*
		Memory pool with filtered alignments and precomputed ZP values
	input_bam_path : const char*
		Input BAM file path
	output_bam_path : const char*
		Output BAM file path
	header : sam_hdr_t*
		Original BAM header
	mapping : ReferenceMapping*
		Reference ID mapping between original and filtered
	num_threads : int
		Number of threads for BAM compression (capped at 4)

	Returns
	-------
	int
		0 on success, -1 on error

	Notes
	-----
	Adds three tags to filtered alignments:
	- ZP:f: Posterior probability from EM algorithm
	- ZS:f: Alignment score (log-likelihood)
	- PM:f: PMD score (if enabled)
	"""
	cdef samFile* in_bam = NULL
	cdef samFile* out_bam = NULL
	cdef hts_idx_t* bam_index = NULL
	cdef LookupTable* lookup = NULL
	cdef WriteBatch* batch = NULL
	cdef sam_hdr_t* filtered_header = NULL

	cdef uint8_t* zp_tag = <uint8_t*>malloc(7)
	cdef uint8_t* zs_tag = <uint8_t*>malloc(7)
	cdef uint8_t* pm_tag = <uint8_t*>malloc(7)
	cdef bint write_pmd_tags = pool.pmd_enabled_for_output
	cdef int tag_size = 14 if not write_pmd_tags else 21
	cdef uint8_t* tag_buffer = <uint8_t*>malloc(tag_size + 4)

	cdef uint64_t alignments_written = 0
	cdef uint64_t alignments_processed = 0
	cdef uint32_t batch_size = 8192
	cdef int cache_size = 256 * 1024 * 1024
	cdef int t

	cdef uint64_t ref_start, ref_count, search_idx
	cdef CompactAlignment* ref_alignments
	cdef hts_itr_t* iterator = NULL
	cdef bam1_t* in_record = NULL
	cdef bam1_t* out_record = NULL
	cdef uint64_t pool_idx
	cdef float zp_value, zs_value, pm_value
	cdef uint32_t ref_id
	cdef uint32_t alignment_position_counter

	# ASCII codes: 'Z' 90, 'P' 80, 'f' 102, 'S' 83, 'M' 77
	zp_tag[0] = 90; zp_tag[1] = 80; zp_tag[2] = 102
	zs_tag[0] = 90; zs_tag[1] = 83; zs_tag[2] = 102
	pm_tag[0] = 80; pm_tag[1] = 77; pm_tag[2] = 102

	verbosity = bf_logging.get_verbosity()

	bf_logging.log(LOG_TAG, "Starting filtered BAM writer (position counter matching enabled)")
	bf_logging.log(LOG_TAG, "Filtered alignments available: %d", pool.alignment_count)
	bf_logging.log(LOG_TAG, "PMD tag output: %s", "enabled" if write_pmd_tags else "disabled")

	if not pool.zp_values_computed or not pool.precomputed_zp_values:
		bf_logging.error("ZP values are not precomputed; aborting filtered BAM write")
		free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
		return -1

	try:
		# Create filtered header with only retained references
		filtered_header = create_filtered_header_efficient(header, mapping)
		if not filtered_header:
			bf_logging.error("Failed to create filtered BAM header")
			free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
			return -1

		bf_logging.log(LOG_TAG, "Created filtered header with %d references (was %d)",
			           filtered_header.n_targets, header.n_targets)

		lookup = create_lookup_table(pool, mapping)
		if not lookup:
			bf_logging.error("Failed to create lookup table for filtered BAM")
			return -1

		if verbosity >= 1:
			bf_logging.log(LOG_TAG, "Lookup table contains %lu entries", lookup.total_count)

		in_bam = hts_open(input_bam_path, "r")
		if not in_bam:
			bf_logging.error("Failed to open input BAM: %s", input_bam_path)
			return -1

		out_bam = hts_open(output_bam_path, "wb1")
		if not out_bam:
			bf_logging.error("Failed to open output BAM: %s", output_bam_path)
			return -1

		# Optionally set thread count (bound to 4)
		if num_threads > 1:
			if num_threads < 4:
				t = num_threads
			else:
				t = 4
			hts_set_threads(out_bam, t)

		# Write filtered header instead of original header
		if sam_hdr_write(out_bam, filtered_header) < 0:
			bf_logging.error("Failed to write filtered BAM header")
			return -1

		bam_index = sam_index_load(in_bam, input_bam_path)
		if not bam_index:
			bf_logging.error("Failed to load BAM index for %s", input_bam_path)
			return -1

		batch = create_write_batch(batch_size)
		if not batch:
			bf_logging.error("Failed to allocate write batch")
			return -1

		in_record = bam_init1()
		if not in_record:
			bf_logging.error("Failed to allocate input BAM record")
			return -1

		try:
			for ref_id in range(mapping.n_original_refs):
				ref_count = lookup.ref_counts[ref_id]
				if ref_count == 0:
					continue

				ref_start = lookup.ref_starts[ref_id]
				ref_alignments = &lookup.alignments[ref_start]
				search_idx = 0

				iterator = sam_itr_queryi(bam_index, ref_id, 0, 0x7fffffff)
				if not iterator:
					if verbosity >= 1:
						bf_logging.warn("No iterator available for reference %d", ref_id)
					continue

				batch.count = 0
				batch.reference_id = ref_id

				alignment_position_counter = 0

				while True:
					if sam_itr_next(in_bam, iterator, in_record) < 0:
						break

					alignments_processed += 1

					if (search_idx < ref_count and 
						ref_alignments[search_idx].original_position == alignment_position_counter):

						pool_idx = ref_alignments[search_idx].pool_index

						if pool_idx >= pool.alignment_count:
							bf_logging.error("Pool index %lu exceeds alignment count %lu", pool_idx, pool.alignment_count)
							hts_itr_destroy(iterator)
							free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
							return -1

						out_record = batch.records[batch.count]
						batch.pool_indices[batch.count] = pool_idx

						if copy_bam_record(in_record, out_record, tag_size) != 0:
							bf_logging.error("Failed to copy BAM record into filtered output")
							hts_itr_destroy(iterator)
							free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
							return -1

						out_record.core.tid = pool.alignments[pool_idx].reference_index

						zp_value = pool.precomputed_zp_values[pool_idx]
						zs_value = pool.alignments[pool_idx].alignment_score
						pm_value = pool.alignments[pool_idx].pmd_score

						if out_record.l_data < tag_size:
							bf_logging.error("Encountered BAM record too small for auxiliary tags")
							hts_itr_destroy(iterator)
							free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
							return -1

						memcpy(tag_buffer, zp_tag, 3)
						memcpy(tag_buffer + 3, &zp_value, 4)
						memcpy(tag_buffer + 7, zs_tag, 3)
						memcpy(tag_buffer + 10, &zs_value, 4)
						if write_pmd_tags:
							memcpy(tag_buffer + 14, pm_tag, 3)
							memcpy(tag_buffer + 17, &pm_value, 4)

						memcpy(out_record.data + (out_record.l_data - tag_size), tag_buffer, tag_size)

						batch.count += 1
						alignments_written += 1

						if batch.count >= batch_size:
							if write_batch_to_bam(out_bam, header, batch) != 0:
								bf_logging.error("Failed to flush batch to filtered BAM")
								hts_itr_destroy(iterator)
								free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
								return -1
							batch.count = 0

						search_idx += 1

					alignment_position_counter += 1

				if batch.count > 0:
					if write_batch_to_bam(out_bam, header, batch) != 0:
						bf_logging.error("Failed to flush final batch to filtered BAM")
						hts_itr_destroy(iterator)
						free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
						return -1
					batch.count = 0

				hts_itr_destroy(iterator)
				iterator = NULL

			bf_logging.log(LOG_TAG, "Wrote %lu of %lu alignments", alignments_written, alignments_processed)

			free(zp_tag); free(zs_tag); free(pm_tag); free(tag_buffer)
			return 0

		finally:
			if in_record:
				bam_destroy1(in_record)
			if iterator:
				hts_itr_destroy(iterator)

	finally:
		if filtered_header:
			sam_hdr_destroy(filtered_header)
		if batch:
			destroy_write_batch(batch)
		if lookup:
			destroy_lookup_table(lookup)
		if bam_index:
			hts_idx_destroy(bam_index)
		if in_bam:
			hts_close(in_bam)
		if out_bam:
			hts_close(out_bam)



cdef int write_bam_with_filtered_header(MemoryPool* pool,
										const char* input_bam_path,
										const char* output_bam_path,
										sam_hdr_t* original_header,
										ReferenceMapping* existing_mapping,
										int bam_write_threads) except -1:
	"""Write a filtered BAM using an updated header (compatibility wrapper).

	This thin wrapper exists for backward compatibility and simply delegates to
	:func:`write_filtered_bam`. It accepts the original header and a mapping
	between original and retained references and forwards the arguments.

	Parameters
	----------
	pool : MemoryPool*
		Memory pool containing filtered alignments and precomputed values.
	input_bam_path : const char*
		Path to the original input BAM file.
	output_bam_path : const char*
		Path where the filtered BAM should be written.
	original_header : sam_hdr_t*
		Original BAM header (used to construct the filtered header).
	existing_mapping : ReferenceMapping*
		Mapping between original and retained reference IDs.
	bam_write_threads : int
		Number of threads to use for BGZF compression (capped by caller).

	Returns
	-------
	int
		0 on success, -1 on error.
	"""
	return write_filtered_bam(pool, input_bam_path, output_bam_path, original_header, existing_mapping, bam_write_threads)
