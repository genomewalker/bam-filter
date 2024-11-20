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
from bam_filter.bam_utils import check_bam_file
import shutil
import uuid
import psutil
from memory_profiler import profile
from numba import jit, prange
import numba

log = logging.getLogger("my_logger")

# hide warnings from numba
logging.getLogger("numba").setLevel(logging.ERROR)


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


def configure_numba_threads(threads=None):
    """
    Configure Numba threading behavior.

    Args:
        threads: Number of threads to use. If None, uses CPU count.
    Returns:
        Original number of threads (for restoration if needed)
    """
    import numba
    from numba import config

    # Store original settings
    original_threads = numba.get_num_threads()

    if threads is None:
        threads = min(os.cpu_count(), 8)  # Limit to reasonable number

    # Set thread count for Numba
    numba.set_num_threads(threads)

    # Enable work stealing for better load balancing
    config.WORKQUEUE_THREAD_ALLOCATION = "steal"

    return original_threads


@jit(nopython=True, parallel=True)
def fast_e_step(prob_vector, subject_indices, slens, n_chunks):
    """E-step using chunked accumulation for deterministic results"""
    # Calculate chunk size
    array_len = len(prob_vector)
    chunk_size = (array_len + n_chunks - 1) // n_chunks  # Ceiling division
    n_chunks_actual = (array_len + chunk_size - 1) // chunk_size

    # Create per-chunk accumulators to ensure deterministic summation
    max_subject_idx = np.max(subject_indices) + 1
    chunk_accumulators = np.zeros((n_chunks_actual, max_subject_idx), dtype=np.float64)

    # First pass: Accumulate within chunks
    for chunk_idx in prange(n_chunks_actual):
        chunk_start = chunk_idx * chunk_size
        chunk_end = min(chunk_start + chunk_size, array_len)

        # Process chunk with local accumulator
        for i in range(chunk_start, chunk_end):
            subject_idx = subject_indices[i]
            chunk_accumulators[chunk_idx, subject_idx] += prob_vector[i]

    # Merge chunks deterministically
    subject_weights = np.zeros(max_subject_idx, dtype=np.float64)
    for subject_idx in range(max_subject_idx):
        for chunk_idx in range(n_chunks_actual):
            subject_weights[subject_idx] += chunk_accumulators[chunk_idx, subject_idx]

    # Calculate s_W deterministically
    s_w = np.empty_like(prob_vector, dtype=np.float64)
    for i in range(array_len):
        s_w[i] = subject_weights[subject_indices[i]] / slens[i]

    return s_w, subject_weights


@jit(nopython=True, parallel=True)
def fast_m_step(prob_vector, s_w, source_indices, n_chunks):
    """M-step using chunked accumulation for deterministic results"""
    # Calculate chunk size
    array_len = len(prob_vector)
    chunk_size = (array_len + n_chunks - 1) // n_chunks
    n_chunks_actual = (array_len + chunk_size - 1) // chunk_size

    # Initialize arrays
    new_prob = prob_vector * s_w
    max_source_idx = np.max(source_indices) + 1

    # Create per-chunk accumulators
    chunk_accumulators = np.zeros((n_chunks_actual, max_source_idx), dtype=np.float64)

    # First pass: Accumulate within chunks
    for chunk_idx in prange(n_chunks_actual):
        chunk_start = chunk_idx * chunk_size
        chunk_end = min(chunk_start + chunk_size, array_len)

        # Process chunk
        for i in range(chunk_start, chunk_end):
            source_idx = source_indices[i]
            chunk_accumulators[chunk_idx, source_idx] += new_prob[i]

    # Merge chunks deterministically
    source_sums = np.zeros(max_source_idx, dtype=np.float64)
    for source_idx in range(max_source_idx):
        for chunk_idx in range(n_chunks_actual):
            source_sums[source_idx] += chunk_accumulators[chunk_idx, source_idx]

    # Normalize probabilities deterministically
    for i in range(array_len):
        source_idx = source_indices[i]
        if source_sums[source_idx] > 0:
            new_prob[i] /= source_sums[source_idx]

    return new_prob


def parallel_em_step(
    prob_vector,
    subject_indices,
    source_indices,
    slens,
    n_chunks,
):
    """Parallel EM step with deterministic operations"""
    # Execute E and M steps
    s_w, subject_weights = fast_e_step(prob_vector, subject_indices, slens, n_chunks)

    new_prob = fast_m_step(prob_vector, s_w, source_indices, n_chunks)

    return new_prob


def squarem_resolve_multimaps(
    data,
    scale=0.9,
    iters=10,
    mmap_dir=None,
    max_memory=None,
    min_improvement=1e-4,
    threads=None,
):
    """SQUAREM-accelerated multimapping resolution."""
    # Configure threading
    original_threads = configure_numba_threads(threads)
    n_threads = numba.get_num_threads()

    try:
        total_memory = max_memory if max_memory else psutil.virtual_memory().total
        current_iter = 0

        # Pre-calculate constants and arrays
        max_subject = np.max(data["subject"])
        max_source = np.max(data["source"])
        source_indices = data["source"].astype(np.int64)
        subject_indices = data["subject"].astype(np.int64)
        slens = data["slen"].astype(np.float64)

        # Calculate optimal number of chunks
        n_chunks = min(32, len(data["prob"]) // 1000)
        if n_chunks < 1:
            n_chunks = 1

        def adaptive_squarem_step():
            """Adaptive SQUAREM step"""
            # First EM update
            r = parallel_em_step(
                data["prob"].astype(np.float64),
                subject_indices,
                source_indices,
                slens,
                n_chunks,
            )

            # Calculate first difference
            v = (
                parallel_em_step(
                    r,
                    subject_indices,
                    source_indices,
                    slens,
                    n_chunks,
                )
                - r
            )

            # Calculate second difference
            w = parallel_em_step(
                r + v,
                subject_indices,
                source_indices,
                slens,
                n_chunks,
            ) - (r + v)

            # Compute step size using stable computation
            v_norm = np.sqrt(np.sum(v * v))  # Explicit sum for stability
            w_v_norm = np.sqrt(np.sum((w - v) * (w - v)))

            if w_v_norm == 0:
                return r

            alpha = v_norm / w_v_norm
            alpha = min(max(alpha, 0.1), 4.0)  # Clip step size for stability

            # SQUAREM update with numerical stability checks
            new_prob = data["prob"] + 2 * alpha * v + (alpha * alpha) * (w - v)

            # Validate results
            if np.any(new_prob < 0) or np.any(np.isnan(new_prob)):
                return r

            return new_prob

        # Main iteration loop with convergence checking
        prev_likelihood = -np.inf

        while True:
            progress_bar = tqdm.tqdm(
                total=9,
                desc=f"Iter {current_iter + 1}",
                unit=" step",
                disable=False,
                leave=False,
                ncols=80,
            )

            # SQUAREM update
            new_prob = adaptive_squarem_step()
            data["prob"] = new_prob
            progress_bar.update(4)

            # Calculate convergence using stable log sum
            current_likelihood = np.sum(np.log(new_prob[new_prob > 0]))
            improvement = (
                (current_likelihood - prev_likelihood) / abs(prev_likelihood)
                if prev_likelihood != -np.inf
                else float("inf")
            )
            prev_likelihood = current_likelihood

            # Query counts
            query_counts = np.bincount(source_indices)
            data["n_aln"] = query_counts[source_indices]

            unique_mask = data["n_aln"] == 1
            non_unique_mask = data["n_aln"] > 1

            if np.all(unique_mask):
                progress_bar.close()
                log.info("::: ::: No more multimapping reads. Early stopping.")
                return data

            # Calculate maximum probabilities
            max_prob = np.zeros(max_source + 1, dtype=np.float64)
            np.maximum.at(max_prob, source_indices, data["prob"])
            data["max_prob"] = max_prob[source_indices] * scale

            # Filter alignments
            to_remove = np.sum(data["prob"] < data["max_prob"])
            filter_mask = data["prob"] >= data["max_prob"]
            final_mask = non_unique_mask & filter_mask

            current_iter += 1
            data["iter"][final_mask] = current_iter

            # Update data array with filtered results
            data = np.concatenate([data[unique_mask], data[final_mask]])

            # Update indices for next iteration
            source_indices = data["source"].astype(np.int64)
            subject_indices = data["subject"].astype(np.int64)
            slens = data["slen"].astype(np.float64)

            progress_bar.update(5)
            progress_bar.close()

            # Logging
            log.info(
                f"::: Iter: {current_iter} - R: {to_remove:,} | U: {np.sum(unique_mask):,} | "
                f"NU: {len(np.unique(data[data['n_aln'] > 1]['source'])):,} | L: {data.shape[0]:,} | "
                f"Improvement: {improvement:.6f}"
            )

            # Check stopping conditions
            if improvement < min_improvement:
                log.info("::: ::: Converged. Stopping.")
                break
            elif iters > 0 and current_iter >= iters:
                log.info("::: ::: Reached maximum iterations. Stopping.")
                break
            elif to_remove == 0:
                log.info("::: ::: No more alignments to remove. Stopping.")
                break

    finally:
        # Restore original Numba thread settings
        numba.set_num_threads(original_threads)

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
    max_read_length=np.inf,
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

        # get total entris in all chunks
        total_entries = sum(len(chunk) for chunk in ref_chunks)
        log.info(f"::: Total entries in all chunks: {total_entries:,}")

        with Manager() as manager:
            entries = manager.dict(dict(entries))

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=num_cores
            ) as executor:
                futures = []
                total_alignments = 0
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

                processed_futures = set()

                for completed_future in concurrent.futures.as_completed(futures):
                    try:
                        alns = completed_future.result()
                        total_alignments += len(alns)
                        write_to_file(
                            alns=alns, out_bam_file=out_bam_file, header=header
                        )
                        processed_futures.add(completed_future)
                    except Exception as e:
                        log.error(f"Error processing batch: {e}")

                    completion_progress_bar.update(1)
                    completed_count += 1
                    completed_future.cancel()

                completion_progress_bar.close()

                # Verify that all futures have been processed
                unprocessed_futures = set(futures) - processed_futures
                if unprocessed_futures:
                    log.error(f"{len(unprocessed_futures)} futures were not processed:")
                    for future in unprocessed_futures:
                        log.error(f"Unprocessed future: {future}")
                else:
                    log.info("All futures processed successfully.")

                log.info(f"Completed {completed_count} out of {len(futures)} futures.")

                log.info(f"Total alignments processed: {total_alignments:,}")

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


def calculate_global_alignment_score(
    num_matches,
    num_mismatches,
    num_gaps,
    gap_extensions,
    match_reward,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
):
    # Calculate the raw alignment score for global alignment
    S = (
        (num_matches * match_reward)
        - (num_mismatches * mismatch_penalty)
        - (num_gaps * gap_open_penalty)
        - (gap_extensions * gap_extension_penalty)
    )
    return S


def normalize_scores(scores):
    # Step 1: Ensure all scores are positive by shifting if necessary
    min_score = min(scores)
    if min_score < 0:
        scores = [
            score - min_score + 1 for score in scores
        ]  # Shift to make all scores positive

    # Step 2: Apply softmax to convert scores to probabilities
    exp_scores = np.exp(scores)
    probabilities = exp_scores / np.sum(exp_scores)

    return probabilities


def process_alignments(
    parms,
    ref_lengths=None,
    percid=90,
    min_read_length=30,
    max_read_length=np.inf,
    threads=1,
    match_reward=1,
    mismatch_penalty=-1,
    gap_open_penalty=1,
    gap_extension_penalty=2,
    tmpdir=None,
):
    """Process BAM alignments and calculate alignment scores."""
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
            alignment_info = []

            for aln in samfile.fetch(
                reference, multiple_iterators=False, until_eof=True
            ):
                query_length = aln.query_length or aln.infer_query_length()
                if min_read_length <= query_length <= max_read_length:
                    num_mismatches = aln.get_tag("NM")
                    pident = (1 - (num_mismatches / query_length)) * 100
                    if pident >= percid:
                        num_matches = query_length - num_mismatches
                        num_gaps = aln.get_tag("XO") if aln.has_tag("XO") else 0
                        gap_extensions = aln.get_tag("XG") if aln.has_tag("XG") else 0

                        score = (
                            (num_matches * match_reward)
                            - (num_mismatches * mismatch_penalty)
                            - (num_gaps * gap_open_penalty)
                            - (gap_extensions * gap_extension_penalty)
                        )

                        alignment_info.append(
                            (aln.query_name, aln.reference_name, score, query_length)
                        )

            if alignment_info:
                scores = np.array([info[2] for info in alignment_info])
                lengths = np.array([info[3] for info in alignment_info])
                normalized_scores = (scores - np.min(scores) + 1) / lengths

                aln_data = dt.Frame(
                    [
                        (info[0], info[1], normalized_scores[i], reference_length)
                        for i, info in enumerate(alignment_info)
                    ],
                    names=["queryId", "subjectId", "bitScore", "slen"],
                )

                aln_data = aln_data[
                    :1, :, dt.by(dt.f.queryId, dt.f.subjectId), dt.sort(-dt.f.bitScore)
                ]
                results.append(aln_data)
            else:
                empty_df += 1

    if results:
        combined_results = dt.rbind(results)
        if tmpdir is not None:
            jay_file = os.path.join(tmpdir, f"{uuid.uuid4()}.jay")
            combined_results.to_jay(jay_file)
            del combined_results
            return (jay_file, empty_df)
        return (combined_results, empty_df)
    return (None, empty_df)


def reassign_reads(
    bam,
    out_files,
    match_reward,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    reference_lengths=None,
    threads=1,
    min_read_count=1,
    min_read_ani=90,
    min_read_length=30,
    max_read_length=np.inf,
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
                        process_alignments,
                        ref_lengths=ref_len_dict,
                        percid=min_read_ani,
                        min_read_length=min_read_length,
                        max_read_length=max_read_length,
                        match_reward=match_reward,
                        mismatch_penalty=mismatch_penalty,
                        gap_open_penalty=gap_open_penalty,
                        gap_extension_penalty=gap_extension_penalty,
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
                        process_alignments,
                        ref_lengths=ref_len_dict,
                        percid=min_read_ani,
                        min_read_length=min_read_length,
                        max_read_length=max_read_length,
                        match_reward=match_reward,
                        mismatch_penalty=mismatch_penalty,
                        gap_open_penalty=gap_open_penalty,
                        gap_extension_penalty=gap_extension_penalty,
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
    refs["sidx"] = dt.Frame(list(range(refs.shape[0])), type=dt.int64)
    refs.key = "subjectId"

    log.info("::: Indexing reads...")
    reads = dt.Frame(list(set(reads)))
    reads.names = ["queryId"]
    reads["qidx"] = dt.Frame(
        [(i + refs.shape[0]) for i in range(reads.shape[0])], type=dt.int64
    )
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
    no_multimaps = squarem_resolve_multimaps(
        init_data,
        iters=reassign_iters,
        scale=reassign_scale,
        mmap_dir=tmp_dir.name,
        max_memory=total_memory,
        threads=threads,
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

    g = dt.Frame(no_multimaps["source"].tolist(), type=dt.int64)

    g.names = ["qidx"]
    reads.key = "qidx"
    q = g[:, :, dt.join(reads)]

    g = dt.Frame(no_multimaps["subject"].tolist(), type=dt.int64)
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
    # remove NA
    s_c = None
    if len(references_m) == 0:
        log.warning("::: No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    # count how many reads are assigned in total after filtering in references_m
    total_reads_refs = 0
    for k in references_m:
        total_reads_refs += references_m[k]

    log.info(f"::: Total refs/reads combination: {total_reads_refs:,}")
    log.info(f"::: Total references: {len(references_m):,}")

    log.info("::: Creating filtered set...")
    entries = defaultdict(set)
    q_query_ids = q[:, "queryId"].to_list()[0]
    s_subject_ids = s[:, "subjectId"].to_list()[0]

    for query_id, subject_id in zip(q_query_ids, s_subject_ids):
        if subject_id in references_m:
            entries[subject_id].add((query_id, subject_id))

    log.info(f"::: ::: Keeping {len(entries):,} references")

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
    bam = check_bam_file(
        bam=args.bam,
        threads=args.threads,
        reference_lengths=args.reference_lengths,
        sort_memory=args.sort_memory,
        sorted_bam=out_files["sorted_bam"],
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
        disable_sort=args.disable_sort,
        tmp_dir=tmp_dir,
    )

    log.info("Done!")
