import datatable as dt
import logging
import numpy as np
import pysam
import tqdm
import sys
import random
from bam_filter.utils import (
    create_empty_output_files,
    sort_keys_by_approx_weight,
    is_debug,
    get_arguments,
    check_tmp_dir_exists,
    handle_warning,
    create_output_files,
    allocate_threads,
)
from multiprocessing import Pool, Manager
from functools import partial
import gc
from collections import defaultdict
import os
import concurrent.futures
import math
import warnings
from bam_filter.sam_utils import check_bam_file
import shutil
import uuid
import psutil

log = logging.getLogger("my_logger")


def estimate_array_size(dtype, shape):
    """
    Estimate the size of a numpy array in bytes given its dtype and shape.
    """
    return np.dtype(dtype).itemsize * np.prod(shape)


def initialize_subject_weights(data, mmap_dir=None, max_memory=None):
    if data.shape[0] > 0:
        total_memory = max_memory if max_memory else psutil.virtual_memory().total
        max_source = np.max(data["source"])

        # Estimate sizes
        sum_weights_size = estimate_array_size(np.float64, (max_source + 1,))

        # Use memory-mapped arrays if estimated size exceeds available memory
        if sum_weights_size > total_memory * 0.8:
            sum_weights = np.memmap(
                os.path.join(mmap_dir, "sum_weights.mmap"),
                dtype="float64",
                mode="w+",
                shape=(max_source + 1,),
            )
        else:
            sum_weights = np.zeros(max_source + 1, dtype="float64")

        # Calculate the sum of weights for each unique source using np.add.at for efficiency
        np.add.at(sum_weights, data["source"], data["var"])

        # Calculate the normalized weights directly
        data["prob"] = data["var"] / sum_weights[data["source"]]

        del sum_weights

        return data
    else:
        return None


def resolve_multimaps(data, scale=0.9, iters=10, mmap_dir=None, max_memory=None):
    total_memory = max_memory if max_memory else psutil.virtual_memory().total
    current_iter = 0

    while True:
        progress_bar = tqdm.tqdm(
            total=9,
            desc=f"Iter {current_iter + 1}",
            unit=" step",
            disable=False,
            leave=False,
            ncols=80,
        )
        log.debug(f"::: Iter: {current_iter + 1} - Getting scores")
        # step 1
        progress_bar.update(1)
        n_alns = data.shape[0]
        log.debug(f"::: Iter: {current_iter + 1} - Total alignment: {n_alns:,}")

        log.debug(f"::: Iter: {current_iter + 1} - Calculating weights...")
        # step 2
        progress_bar.update(1)
        max_subject = np.max(data["subject"])
        subject_weights_size = estimate_array_size(
            np.float64, (np.int64(max_subject) + 1,)
        )

        if subject_weights_size > total_memory * 0.8:
            # Use memory-mapped arrays
            subject_weights = np.memmap(
                os.path.join(mmap_dir, f"subject_weights_{current_iter}.mmap"),
                dtype="float64",
                mode="w+",
                shape=(np.int64(max_subject) + 1,),
            )
        else:
            # Use in-memory arrays
            subject_weights = np.zeros(np.int64(max_subject) + 1, dtype="float64")

        np.add.at(subject_weights, data["subject"], data["prob"])
        data["s_W"] = subject_weights[data["subject"]] / data["slen"]
        del subject_weights

        log.debug(f"::: Iter: {current_iter + 1} - Calculating probabilities")
        # step 3
        progress_bar.update(1)
        new_prob = data["prob"] * data["s_W"]
        max_source = np.max(data["source"])
        prob_sum_array_size = estimate_array_size(
            np.float64, (np.int64(max_source) + 1,)
        )

        if prob_sum_array_size > total_memory * 0.8:
            # Use memory-mapped arrays
            prob_sum_array = np.memmap(
                os.path.join(mmap_dir, f"prob_sum_array_{current_iter}.mmap"),
                dtype="float64",
                mode="w+",
                shape=(np.int64(max_source) + 1,),
            )
        else:
            # Use in-memory arrays
            prob_sum_array = np.zeros(np.int64(max_source) + 1, dtype="float64")

        np.add.at(prob_sum_array, data["source"], new_prob)
        data["prob"] = new_prob / prob_sum_array[data["source"]]
        del new_prob, prob_sum_array

        log.debug("Calculating query counts")
        # step 4
        progress_bar.update(1)
        query_counts = np.bincount(data["source"])

        log.debug("Calculating query counts array")
        # step 5
        progress_bar.update(1)
        query_counts_array_size = estimate_array_size(
            np.int64, (np.int64(max_source) + 1,)
        )

        if query_counts_array_size > total_memory * 0.8:
            # Use memory-mapped arrays
            query_counts_array = np.memmap(
                os.path.join(mmap_dir, f"query_counts_array_{current_iter}.mmap"),
                dtype="int64",
                mode="w+",
                shape=(np.int64(max_source) + 1,),
            )
        else:
            # Use in-memory arrays
            query_counts_array = np.zeros(np.int64(max_source) + 1, dtype="int64")

        np.add.at(query_counts_array, data["source"], 1)

        log.debug(
            f"::: Iter: {current_iter + 1} - Calculating number of alignments per query"
        )
        # step 6
        progress_bar.update(1)
        data["n_aln"] = query_counts_array[data["source"]]

        unique_mask = data["n_aln"] == 1
        non_unique_mask = data["n_aln"] > 1

        if np.all(unique_mask):
            # step 7
            progress_bar.close()
            log.info("::: ::: No more multimapping reads. Early stopping.")
            return data

        log.debug("Calculating max_prob")
        max_prob_size = estimate_array_size(np.float64, (np.int64(max_source) + 1,))

        if max_prob_size > total_memory * 0.8:
            # Use memory-mapped arrays
            max_prob = np.memmap(
                os.path.join(mmap_dir, f"max_prob_{current_iter}.mmap"),
                dtype="float64",
                mode="w+",
                shape=(np.int64(max_source) + 1,),
            )
        else:
            # Use in-memory arrays
            max_prob = np.zeros(np.int64(max_source) + 1, dtype="float64")

        np.maximum.at(max_prob, data["source"], data["prob"])
        data["max_prob"] = max_prob[data["source"]] * scale
        del max_prob

        log.debug(
            f"::: Iter: {current_iter + 1} - Removing alignments with lower probability"
        )
        # step 8
        progress_bar.update(1)
        to_remove = np.sum(data["prob"] < data["max_prob"])

        filter_mask = data["prob"] >= data["max_prob"]
        final_mask = non_unique_mask & filter_mask

        current_iter += 1
        data["iter"][final_mask] = current_iter

        # Concatenate unique and filtered non-unique data
        data = np.concatenate([data[unique_mask], data[final_mask]])

        query_counts = np.bincount(data["source"])
        total_n_unique = np.sum(query_counts <= 1)

        keep_processing = to_remove != 0
        log.debug(f"::: Iter: {current_iter} - Removed {to_remove:,} alignments")
        log.debug(
            f"::: Iter: {current_iter} - Total mapping queries: {np.sum(unique_mask):,}"
        )
        log.debug(
            f"::: Iter: {current_iter} - New unique mapping queries: {total_n_unique:,}"
        )
        log.debug(f"::: Iter: {current_iter} - Alns left: {data.shape[0]:,}")
        # step 9
        progress_bar.update(1)
        progress_bar.close()
        log.info(
            f"::: Iter: {current_iter} - R: {to_remove:,} | U: {np.sum(unique_mask):,} | NU: {total_n_unique:,} | L: {data.shape[0]:,}"
        )
        log.debug(f"::: Iter: {current_iter} - done!")

        if iters > 0 and current_iter >= iters:
            log.info("::: ::: Reached maximum iterations. Stopping.")
            break
        elif to_remove == 0:
            log.info("::: ::: No more alignments to remove. Stopping.")
            break

    return data


def write_to_file(alns, out_bam_file, header=None):
    for aln in alns:
        out_bam_file.write(pysam.AlignedSegment.fromstring(aln, header))


def process_references_batch(references, entries, bam, refs_idx, threads=1):
    alns = []
    s_threads = threads
    with pysam.AlignmentFile(bam, "rb", threads=s_threads) as samfile:
        for reference in references:
            r_ids = entries[reference]
            for aln in samfile.fetch(
                reference=reference, multiple_iterators=False, until_eof=True
            ):
                if (aln.query_name, reference) in r_ids:
                    aln.reference_id = refs_idx[aln.reference_name]
                    alns.append(aln.to_string())

    return alns


def write_reassigned_bam(
    bam,
    ref_counts,
    out_files,
    threads,
    entries,
    sort_memory="1G",
    sort_by_name=False,
    min_read_ani=90,
    min_read_length=30,
    max_read_length=np.Inf,
    disable_sort=False,
):
    out_bam = out_files["bam_reassigned"]
    p_threads, s_threads = allocate_threads(threads, 1, 4)
    with pysam.AlignmentFile(bam, "rb", threads=s_threads) as samfile:
        references = list(entries.keys())
        refs_dict = {x: samfile.get_reference_length(x) for x in references}
        header = samfile.header

    log.info("::: Getting reference names and lengths...")
    (ref_names, ref_lengths) = zip(*refs_dict.items())
    refs_idx = {sys.intern(str(x)): i for i, x in enumerate(ref_names)}
    write_threads = s_threads

    new_header = header.to_dict()

    log.info("::: Creating new header...")
    ref_names_set = set(ref_names)
    new_header["SQ"] = [x for x in new_header["SQ"] if x["SN"] in ref_names_set]

    name_index = {name: idx for idx, name in enumerate(ref_names)}
    new_header["SQ"].sort(key=lambda x: name_index[x["SN"]])
    new_header["HD"]["SO"] = "unsorted"

    with pysam.AlignmentFile(
        out_files["bam_reassigned_tmp"],
        "wb",
        referencenames=list(ref_names),
        referencelengths=list(ref_lengths),
        threads=write_threads,
        header=new_header,
    ) as out_bam_file:

        num_cores = p_threads

        log.info("::: Creating reference chunks with uniform read amounts...")

        ref_chunks = sort_keys_by_approx_weight(
            input_dict=ref_counts,
            scale=1,
            num_cores=num_cores,
            verbose=False,
            max_entries_per_chunk=1_000_000,
        )

        log.info(f"::: Using {num_cores} processes to write {len(ref_chunks)} chunk(s)")

        with Manager() as manager:
            entries = manager.dict(dict(entries))

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=num_cores
            ) as executor:
                futures = []
                for batch_references in tqdm.tqdm(
                    ref_chunks,
                    total=len(ref_chunks),
                    desc="Submitted batches",
                    unit="batch",
                    leave=False,
                    ncols=80,
                    disable=is_debug(),
                ):
                    future = executor.submit(
                        process_references_batch,
                        batch_references,
                        entries,
                        bam,
                        refs_idx,
                        s_threads,
                    )
                    futures.append(future)

                log.info("::: Collecting batches...")

                completion_progress_bar = tqdm.tqdm(
                    total=len(futures),
                    desc="Completed",
                    unit="batch",
                    leave=False,
                    ncols=80,
                    disable=is_debug(),
                )
                completed_count = 0

                for completed_future in concurrent.futures.as_completed(futures):
                    alns = completed_future.result()
                    write_to_file(alns=alns, out_bam_file=out_bam_file, header=header)

                    completion_progress_bar.update(1)
                    completed_count += 1
                    completed_future.cancel()
                    gc.collect()

                completion_progress_bar.close()
    entries = None
    gc.collect()

    if not disable_sort:
        log.info("::: ::: Sorting BAM file...")
        w_threads = max(4, s_threads)
        if sort_by_name:
            log.info("::: ::: Sorting by name...")
            pysam.sort(
                "-n",
                "-@",
                str(w_threads),
                "-m",
                str(sort_memory),
                "-o",
                out_bam,
                out_files["bam_reassigned_tmp"],
            )
        else:
            pysam.sort(
                "-@",
                str(w_threads),
                "-m",
                str(sort_memory),
                "-o",
                out_bam,
                out_files["bam_reassigned_tmp"],
            )

        logging.info("BAM index not found. Indexing...")

        pysam.index(
            "-c",
            "-@",
            str(threads),
            out_bam,
        )

        os.remove(out_files["bam_reassigned_tmp"])
    else:
        logging.info("Skipping BAM file sorting...")
        shutil.move(out_files["bam_reassigned_tmp"], out_bam)


def calculate_alignment_score(
    num_matches,
    num_mismatches,
    num_gaps,
    gap_extensions,
    match_reward,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    precomputed_factor,  # This is lambda_value * match_reward / math.log(2)
    precomputed_log_K,  # This is math.log(K_value) / math.log(2)
):
    S = (
        (num_matches * match_reward)
        - (num_mismatches * mismatch_penalty)
        - (num_gaps * gap_open_penalty)
        - (gap_extensions * gap_extension_penalty)
    )

    bit_score = precomputed_factor * S - precomputed_log_K

    return bit_score


def get_bam_data(
    parms,
    ref_lengths=None,
    percid=90,
    min_read_length=30,
    max_read_length=np.Inf,
    threads=1,
    match_reward=1,
    mismatch_penalty=-1,
    gap_open_penalty=1,
    gap_extension_penalty=2,
    lambda_value=1.02,
    K_value=0.21,
    tmpdir=None,
):
    precomputed_factor = lambda_value * match_reward / math.log(2)
    precomputed_log_K = math.log(K_value) / math.log(2)

    bam, references = parms
    dt.options.progress.enabled = False
    dt.options.progress.clear_on_success = True

    results = []
    empty_df = 0

    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
        if ref_lengths is None:
            reference_lengths = {
                reference: np.int64(samfile.get_reference_length(reference))
                for reference in references
            }
        else:
            reference_lengths = {
                reference: np.int64(ref_lengths[reference]) for reference in references
            }

        for reference in references:
            reference_length = reference_lengths[reference]
            aln_data = []
            fetch = samfile.fetch(reference, multiple_iterators=False, until_eof=True)

            for aln in fetch:
                query_length = aln.query_length or aln.infer_query_length()
                if query_length >= min_read_length and query_length <= max_read_length:
                    num_mismatches = aln.get_tag("NM")
                    pident = (1 - (num_mismatches / query_length)) * 100
                    if pident >= percid:
                        num_matches = query_length - num_mismatches
                        num_gaps = aln.get_tag("XO") if aln.has_tag("XO") else 0
                        gap_extensions = aln.get_tag("XG") if aln.has_tag("XG") else 0

                        S = (
                            (num_matches * match_reward)
                            - (num_mismatches * mismatch_penalty)
                            - (num_gaps * gap_open_penalty)
                            - (gap_extensions * gap_extension_penalty)
                        )
                        bit_score = precomputed_factor * S - precomputed_log_K

                        aln_data.append(
                            (
                                aln.query_name,
                                aln.reference_name,
                                bit_score,
                                reference_length,
                            )
                        )

            if aln_data:
                aln_data_dt = dt.Frame(
                    aln_data, names=["queryId", "subjectId", "bitScore", "slen"]
                )
                aln_data_dt = aln_data_dt[
                    :1, :, dt.by(dt.f.queryId, dt.f.subjectId), dt.sort(-dt.f.bitScore)
                ]
                results.append(aln_data_dt)
            else:
                empty_df += 1

    if results:
        combined_results = dt.rbind(results)
        if tmpdir is not None:
            uuid_name = str(uuid.uuid4())
            jay_file = os.path.join(tmpdir, f"{uuid_name}.jay")
            combined_results.to_jay(jay_file)
            del combined_results
            return (jay_file, empty_df)
        else:
            return (combined_results, empty_df)
    else:
        return (None, empty_df)


def reassign_reads(
    bam,
    out_files,
    match_reward,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    lambda_value,
    K_value,
    reference_lengths=None,
    threads=1,
    min_read_count=1,
    min_read_ani=90,
    min_read_length=30,
    max_read_length=np.Inf,
    reassign_iters=25,
    reassign_scale=0.9,
    sort_memory="4G",
    disable_sort=False,
    tmp_dir=None,
    max_memory=None,
):

    p_threads, s_threads = allocate_threads(threads, 1, 4)
    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True
    if threads > 1:
        dt.options.nthreads = p_threads
    else:
        dt.options.nthreads = 1

    log.info("::: Loading BAM file")
    save = pysam.set_verbosity(0)

    log.info(f"::: IO Threads: {s_threads} | Processing Threads: {p_threads}")

    with pysam.AlignmentFile(bam, "rb", threads=s_threads) as samfile:
        references = samfile.references
        pysam.set_verbosity(save)

        total_refs = samfile.nreferences
        log.info(f"::: Found {total_refs:,} reference sequences")

        log.info(f"::: Removing references with less than {min_read_count} reads...")
        if reference_lengths is not None:
            ref_len_dt = dt.fread(reference_lengths)
            ref_len_dt.names = ["subjectId", "slen"]
            ref_len_dict = dict(
                zip(
                    ref_len_dt["subjectId"].to_list()[0],
                    ref_len_dt["slen"].to_list()[0],
                )
            )
            if not set(references).issubset(set(ref_len_dict.keys())):
                logging.error(
                    "The BAM file contains references not found in the reference lengths file"
                )
                sys.exit(1)
        else:
            ref_len_dict = None

        index_statistics = samfile.get_index_statistics()

    references_m = {
        chrom.contig: chrom.mapped
        for chrom in tqdm.tqdm(
            [chrom for chrom in index_statistics if chrom.mapped >= min_read_count],
            desc="Filtering references",
            total=len(index_statistics),
            unit="chrom",
            leave=False,
            ncols=80,
            unit_scale=True,
            unit_divisor=1000,
        )
    }
    del index_statistics
    n_alns = sum(references_m.values())
    log.info(f"::: Kept {n_alns:,} alignments")
    references = list(references_m.keys())

    if len(references) == 0:
        log.warning("::: No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    log.info(f"::: Keeping {len(references):,} references")

    log.info("::: Creating reference chunks with uniform read amounts...")
    ref_chunks = sort_keys_by_approx_weight(
        input_dict=references_m,
        scale=1,
        num_cores=threads,
        refinement_steps=10,
        verbose=False,
        max_entries_per_chunk=100_000_000,
    )

    log.info(f"::: ::: Created {len(ref_chunks):,} chunks")
    ref_chunks = random.sample(ref_chunks, len(ref_chunks))
    dt.options.progress.enabled = False
    dt.options.progress.clear_on_success = True
    dt.options.nthreads = 1
    del references_m
    gc.collect()

    parms = list(zip([bam] * len(ref_chunks), ref_chunks))

    log.info("::: Extracting reads from BAM file...")
    if is_debug():
        data = list(
            tqdm.tqdm(
                map(
                    partial(
                        get_bam_data,
                        ref_lengths=ref_len_dict,
                        percid=min_read_ani,
                        min_read_length=min_read_length,
                        max_read_length=max_read_length,
                        match_reward=match_reward,
                        mismatch_penalty=mismatch_penalty,
                        gap_open_penalty=gap_open_penalty,
                        gap_extension_penalty=gap_extension_penalty,
                        lambda_value=lambda_value,
                        K_value=K_value,
                        threads=s_threads,
                        tmpdir=out_files["tmp_dir"],
                    ),
                    parms,
                    chunksize=1,
                ),
                total=len(parms),
                leave=False,
                ncols=80,
                desc="Chunks processed",
            )
        )
    else:
        p = Pool(p_threads)
        data = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(
                        get_bam_data,
                        ref_lengths=ref_len_dict,
                        percid=min_read_ani,
                        min_read_length=min_read_length,
                        max_read_length=max_read_length,
                        match_reward=match_reward,
                        mismatch_penalty=mismatch_penalty,
                        gap_open_penalty=gap_open_penalty,
                        gap_extension_penalty=gap_extension_penalty,
                        lambda_value=lambda_value,
                        K_value=K_value,
                        threads=s_threads,
                        tmpdir=out_files["tmp_dir"],
                    ),
                    parms,
                    chunksize=1,
                ),
                total=len(parms),
                leave=False,
                ncols=80,
                desc="Chunks processed",
            )
        )

        p.close()
        p.join()

    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True
    if threads > 1:
        dt.options.nthreads = threads - 1
    else:
        dt.options.nthreads = 1

    log.info("::: Collecting results...")
    reads = set()
    refs = set()
    empty_df = 0

    for i in tqdm.tqdm(range(len(data)), total=len(data), leave=False, ncols=80):
        empty_df += data[i][1]
        df = dt.fread(data[i][0])
        data[i] = df
        query_ids = df[:, "queryId"].to_list()[0]
        subject_ids = df[:, "subjectId"].to_list()[0]

        reads.update(query_ids)
        refs.update(subject_ids)

        del df

    reads = list(reads)
    refs = list(refs)

    log.info(f"::: ::: Removed {empty_df:,} references without alignments")

    log.info("::: Indexing references...")
    refs = dt.Frame(list(set(refs)))
    refs.names = ["subjectId"]
    refs["sidx"] = dt.Frame(list(range(refs.shape[0])))
    refs.key = "subjectId"

    log.info("::: Indexing reads...")
    reads = dt.Frame(list(set(reads)))
    reads.names = ["queryId"]
    reads["qidx"] = dt.Frame([(i + refs.shape[0]) for i in range(reads.shape[0])])
    reads.key = "queryId"

    log.info("::: Allocating data...")
    total_rows = np.int64(sum(x.shape[0] for x in data))

    n_reads_0 = reads.shape[0]
    n_refs_0 = refs.shape[0]

    n_alns_0 = 0
    current_index = 0

    dtype = np.dtype(
        [
            ("source", "int64"),
            ("subject", "int64"),
            ("var", "float32"),
            ("slen", "int64"),
            ("s_W", "float32"),
            ("prob", "float32"),
            ("iter", "int32"),
            ("n_aln", "int64"),
            ("max_prob", "float32"),
        ]
    )

    total_memory = max_memory if max_memory else psutil.virtual_memory().total
    array_size = estimate_array_size(dtype, (total_rows,))

    if array_size > total_memory * 0.8:
        log.warning("::: Using memory-mapped arrays")
        # Use memory-mapped arrays
        m = np.memmap(
            os.path.join(tmp_dir.name, "m.mmap"),
            dtype=dtype,
            mode="w+",
            shape=(total_rows,),
        )
    else:
        # Use in-memory arrays
        m = np.zeros(total_rows, dtype=dtype)

    for i in tqdm.tqdm(
        range(len(data)),
        total=len(data),
        desc="Processing batches",
        unit="batch",
        disable=is_debug(),
        leave=False,
        ncols=80,
    ):
        x = data.pop(0)
        if x.shape[0] > 0:
            x = x[:, :, dt.join(reads)]
            x = x[:, :, dt.join(refs)]
            n_alns_0 += x.shape[0]

            x_processed = x[
                :, [dt.f.qidx, dt.f.sidx, dt.f.bitScore, dt.f.slen]
            ].to_numpy()
            num_rows = x_processed.shape[0]

            m["source"][current_index : current_index + num_rows] = x_processed[:, 0]
            m["subject"][current_index : current_index + num_rows] = x_processed[:, 1]
            m["var"][current_index : current_index + num_rows] = x_processed[:, 2]
            m["slen"][current_index : current_index + num_rows] = x_processed[:, 3]
            current_index += num_rows
            del x

    del data
    gc.collect()

    log.info(
        f"::: References: {n_refs_0:,} | Reads: {n_reads_0:,} | Alignments: {n_alns_0:,}"
    )
    log.info("::: Initializing data structures...")
    init_data = initialize_subject_weights(
        m, mmap_dir=tmp_dir.name, max_memory=total_memory
    )
    if reassign_iters > 0:
        log.info(f"::: Reassigning reads with {reassign_iters} iterations")
    else:
        log.info("::: Reassigning reads until convergence")
    no_multimaps = resolve_multimaps(
        init_data,
        iters=reassign_iters,
        scale=reassign_scale,
        mmap_dir=tmp_dir.name,
        max_memory=total_memory,
    )

    n_reads = len(list(set(no_multimaps["source"])))
    n_refs = len(list(set(no_multimaps["subject"])))
    n_alns = no_multimaps.shape[0]
    log.info(
        f"::: References: {n_refs:,} | Reads: {n_reads:,} | Alignments: {n_alns:,}"
    )

    log.info(
        f'::: Unique mapping reads: {no_multimaps[no_multimaps["n_aln"] == 1].shape[0]:,} | Multimapping reads: {len(np.unique(no_multimaps[no_multimaps["n_aln"] > 1]["source"])):,}'
    )

    log.info("::: Mapping back indices...")
    if threads > 1:
        dt.options.nthreads = threads - 1
    else:
        dt.options.nthreads = 1

    g = dt.Frame(no_multimaps["source"])
    g.names = ["qidx"]
    reads.key = "qidx"
    q = g[:, :, dt.join(reads)]

    g = dt.Frame(no_multimaps["subject"])
    g.names = ["sidx"]
    refs.key = "sidx"
    s = g[:, :, dt.join(refs)]

    log.info("::: Calculating reads per subject...")
    s_c = s[:, dt.count(dt.f.subjectId), dt.by(dt.f.subjectId)]
    s_c.names = ["subjectId", "counts"]
    references_m = dict()
    log.info(f"::: Removing references with less than {min_read_count:,}...")
    for i, k in zip(s_c[:, "subjectId"].to_list()[0], s_c[:, "counts"].to_list()[0]):
        if k >= min_read_count:
            references_m[i] = k
    log.info(f"::: ::: Keeping {len(references_m):,} references")
    s_c = None
    if len(references_m) == 0:
        log.warning("::: No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    log.info("::: Creating filtered set...")
    entries = defaultdict(set)
    q_query_ids = q[:, "queryId"].to_list()[0]
    s_subject_ids = s[:, "subjectId"].to_list()[0]

    for query_id, subject_id in zip(q_query_ids, s_subject_ids):
        if subject_id in references_m:
            entries[subject_id].add((query_id, subject_id))

    del m
    del init_data
    del no_multimaps
    del q
    del s
    del q_query_ids
    del s_subject_ids
    gc.collect()

    log.info("::: Writing to BAM file...")
    write_reassigned_bam(
        bam=bam,
        ref_counts=references_m,
        out_files=out_files,
        threads=threads,
        entries=entries,
        sort_memory=sort_memory,
        min_read_ani=min_read_ani,
        disable_sort=disable_sort,
    )


def reassign(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%H:%M:%S",
    )

    args = get_arguments()
    if args.max_read_length < args.min_read_length:
        logging.error("Maximum read length cannot be less than minimum read length")
        sys.exit(1)
    bam = args.bam
    tmp_dir = check_tmp_dir_exists(args.tmp_dir)
    log.info("Temporary directory: %s", tmp_dir.name)
    out_files = create_output_files(
        prefix=args.prefix,
        bam=args.bam,
        tmp_dir=tmp_dir,
        mode="reassign",
        bam_reassigned=args.bam_reassigned,
    )
    sorted_bam = bam.replace(".bam", ".bf-sorted.bam")
    bam = check_bam_file(
        bam=args.bam,
        threads=args.threads,
        reference_lengths=args.reference_lengths,
        sort_memory=args.sort_memory,
    )
    if bam is None:
        logging.warning("No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    if args.debug:
        warnings.showwarning = handle_warning
    else:
        warnings.filterwarnings("ignore")
    logging.info("Resolving multi-mapping reads...")
    reassign_reads(
        bam=bam,
        threads=args.threads,
        reference_lengths=args.reference_lengths,
        min_read_count=args.min_read_count,
        min_read_ani=args.min_read_ani,
        min_read_length=args.min_read_length,
        max_read_length=args.max_read_length,
        reassign_iters=args.reassign_iters,
        reassign_scale=args.reassign_scale,
        max_memory=args.max_memory,
        sort_memory=args.sort_memory,
        out_files=out_files,
        match_reward=args.match_reward,
        mismatch_penalty=args.mismatch_penalty,
        gap_open_penalty=args.gap_open_penalty,
        gap_extension_penalty=args.gap_extension_penalty,
        lambda_value=args.lambda_value,
        K_value=args.K_value,
        disable_sort=args.disable_sort,
        tmp_dir=tmp_dir,
    )
    if os.path.exists(sorted_bam):
        os.remove(sorted_bam)
    if os.path.exists(sorted_bam + ".bai"):
        os.remove(sorted_bam + ".bai")
    elif os.path.exists(sorted_bam + ".csi"):
        os.remove(sorted_bam + ".csi")
    log.info("Done!")
