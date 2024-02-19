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
)
from multiprocessing import Pool, Manager, cpu_count
from functools import partial
import numpy.lib.recfunctions as rf
import gc
from collections import defaultdict
import os
import concurrent.futures
import math
import warnings
from bam_filter.sam_utils import check_bam_file

# import cProfile as prof
# import pstats

log = logging.getLogger("my_logger")


def initialize_subject_weights(data):
    if data.shape[0] > 0:
        # Add a new column for inverse sequence length
        data["s_W"] = 1 / data["slen"]

        # Calculate the sum of weights for each unique source
        sum_weights = np.zeros(int(np.max(data["source"])) + 1)
        np.add.at(sum_weights, data["source"], data["var"])

        # Calculate the normalized weights based on the sum for each query
        query_sum_var = sum_weights[data["source"]]
        data["prob"] = data["var"] / query_sum_var

        return data
    else:
        return None


def resolve_multimaps(data, scale=0.9, iters=10):
    current_iter = 0
    while True:
        progress_bar = tqdm.tqdm(
            total=9,
            desc=f"Iter {current_iter + 1}",
            unit="step",
            disable=False,  # Replace with your logic or a boolean value
            leave=False,
            ncols=80,
        )
        log.debug(f"::: Iter: {current_iter + 1} - Getting scores")
        progress_bar.update(1)
        n_alns = data.shape[0]
        log.debug(f"::: Iter: {current_iter + 1} - Total alignment: {n_alns:,}")

        # Calculate the weights for each subject
        log.debug(f"::: Iter: {current_iter + 1} - Calculating weights...")
        progress_bar.update(1)
        subject_weights = np.zeros(int(np.max(data["subject"])) + 1)
        np.add.at(subject_weights, data["subject"], data["prob"])
        data["s_W"] = subject_weights[data["subject"]] / data["slen"]
        subject_weights = None

        log.debug(f"::: Iter: {current_iter + 1} - Calculating probabilities")
        progress_bar.update(1)
        # Calculate the alignment probabilities
        new_prob = data["prob"] * data["s_W"]
        log.debug("Calculating sum of probabilities")
        progress_bar.update(1)
        prob_sum = data["prob"] * data["s_W"]
        prob_sum_array = np.zeros(int(np.max(data["source"])) + 1)
        np.add.at(prob_sum_array, data["source"], prob_sum)
        prob_sum = None

        # data["prob_sum"] = prob_sum_array[data["source"]]
        data["prob"] = new_prob / prob_sum_array[data["source"]]
        prob_sum_array = None

        log.debug("Calculating query counts")
        progress_bar.update(1)
        # Calculate how many alignments are in each query
        # query_counts = np.zeros(int(np.max(data["source"])) + 1)
        # np.add.at(query_counts, data["source"], 1)

        query_counts = np.bincount(data["source"])

        log.debug("Calculating query counts array")
        progress_bar.update(1)
        # Use a separate array for query counts
        query_counts_array = np.zeros(int(np.max(data["source"])) + 1)
        np.add.at(
            query_counts_array,
            data["source"],
            query_counts[data["source"]],
        )

        log.debug(
            f"::: Iter: {current_iter + 1} - Calculating number of alignments per query"
        )
        progress_bar.update(1)
        data["n_aln"] = query_counts_array[data["source"]]

        log.debug("Calculating unique alignments")
        data["n_aln"] = query_counts_array[data["source"]]
        data_unique = data[data["n_aln"] == 1]
        n_unique = data_unique.shape[0]

        if n_unique == data.shape[0]:
            progress_bar.close()
            log.info("::: ::: No more multimapping reads. Early stopping.")
            return data
        data = data[(data["n_aln"] > 1) & (data["prob"] > 0)]

        # total_n_unique = np.sum(query_counts_array[data["source"]] <= 1)

        query_counts = None
        query_counts_array = None

        log.debug("Calculating max_prob")
        # Keep the ones that have a probability higher than the maximum scaled probability
        max_prob = np.zeros(int(np.max(data["source"])) + 1)
        np.maximum.at(max_prob, data["source"], data["prob"])

        data["max_prob"] = max_prob[data["source"]]
        data["max_prob"] = data["max_prob"] * scale
        # data["max_prob"] = max_prob[data["source"]]
        log.debug(
            f"::: Iter: {current_iter + 1} - Removing alignments with lower probability"
        )
        progress_bar.update(1)
        to_remove = np.sum(data["prob"] < data["max_prob"])

        data = data[data["prob"] >= data["max_prob"]]
        max_prob = None

        # Update the iteration count in the function call
        current_iter += 1
        data["iter"] = current_iter
        data_unique["iter"] = current_iter

        query_counts = np.bincount(data["source"])
        total_n_unique = np.sum(query_counts[data["source"]] <= 1)

        # data_unique["iter"] = current_iter

        # data = np.concatenate([data, data_unique])
        data = np.concatenate([data, data_unique])
        data_unique = None

        keep_processing = to_remove != 0
        log.debug(f"::: Iter: {current_iter} - Removed {to_remove:,} alignments")
        log.debug(f"::: Iter: {current_iter} - Total mapping queries: {n_unique:,}")
        log.debug(
            f"::: Iter: {current_iter} - New unique mapping queries: {total_n_unique:,}"
        )
        log.debug(f"::: Iter: {current_iter} - Alns left: {data.shape[0]:,}")
        progress_bar.update(1)
        progress_bar.close()
        log.info(
            f"::: Iter: {current_iter} - R: {to_remove:,} | U: {total_n_unique:,} | NU: {n_unique:,} | L: {data.shape[0]:,}"
        )
        log.debug(f"::: Iter: {current_iter} - done!")

        if iters > 0 and current_iter >= iters:
            log.info("::: ::: Reached maximum iterations. Stopping.")
            break
        elif not keep_processing:
            log.info("::: ::: No more alignments to remove. Stopping.")
            break
    return data


# def write_reassigned_bam(
#     bam, out_files, threads, entries, sort_memory="1G", min_read_ani=90
# ):
#     if out_files["bam_reassigned"] is not None:
#         out_bam = out_files["bam_reassigned"]
#     else:
#         out_bam = out_files["bam_reassigned_sorted"]

#     samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
#     references = list(entries.keys())
#     refs_dict = {x: samfile.get_reference_length(x) for x in list(entries.keys())}

#     (ref_names, ref_lengths) = zip(*refs_dict.items())

#     refs_idx = {sys.intern(str(x)): i for i, x in enumerate(ref_names)}
#     if threads > 4:
#         write_threads = 4
#     else:
#         write_threads = threads

#     out_bam_file = pysam.AlignmentFile(
#         out_files["bam_reassigned_tmp"],
#         "wb",
#         referencenames=list(ref_names),
#         referencelengths=list(ref_lengths),
#         threads=write_threads,
#     )

#     for reference in tqdm.tqdm(
#         references,
#         total=len(references),
#         leave=False,
#         ncols=80,
#         desc="References processed",
#     ):
#         r_ids = entries[reference]
#         for aln in samfile.fetch(
#             reference=reference, multiple_iterators=False, until_eof=True
#         ):
#             # ani_read = (1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100
#             if (aln.query_name, reference) in r_ids:
#                 aln.reference_id = refs_idx[aln.reference_name]
#                 out_bam_file.write(aln)
#     out_bam_file.close()


# def write_to_file(alns, out_bam_file):
#     for aln in tqdm.tqdm(alns, total=len(alns), leave=False, ncols=80, desc="Writing"):
#         out_bam_file.write(aln)


def write_to_file(alns, out_bam_file, header=None):
    for aln in alns:
        out_bam_file.write(pysam.AlignedSegment.fromstring(aln, header))


def process_references_batch(references, entries, bam, refs_idx, threads=4):
    alns = []
    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
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
):
    # if out_files["bam_reassigned"] is not None:
    #     out_bam = out_files["bam_reassigned"]
    # else:
    #     out_bam = out_files["bam_reassigned_sorted"]
    out_bam = out_files["bam_reassigned"]
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    references = list(entries.keys())
    refs_dict = {x: samfile.get_reference_length(x) for x in references}
    header = samfile.header
    samfile.close()
    (ref_names, ref_lengths) = zip(*refs_dict.items())

    refs_idx = {sys.intern(str(x)): i for i, x in enumerate(ref_names)}
    if threads > 4:
        write_threads = 4
    else:
        write_threads = threads

    out_bam_file = pysam.AlignmentFile(
        out_files["bam_reassigned_tmp"],
        "wb",
        referencenames=list(ref_names),
        referencelengths=list(ref_lengths),
        threads=write_threads,
    )

    num_cores = min(threads, cpu_count())
    # batch_size = len(references) // num_cores + 1  # Ensure non-zero batch size
    # batch_size = calc_chunksize(n_workers=num_cores, len_iterable=len(references))
    log.info("::: Creating reference chunks with uniform read amounts...")

    ref_chunks = sort_keys_by_approx_weight(
        input_dict=ref_counts, scale=1, num_cores=num_cores
    )
    num_cores = min(num_cores, len(ref_chunks))
    log.info(f"::: Using {num_cores} cores to write {len(ref_chunks)} chunk(s)")

    with Manager() as manager:
        # Use Manager to create a read-only proxy for the dictionary
        entries = manager.dict(dict(entries))

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
            # Use ProcessPoolExecutor to parallelize the processing of references in batches
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
                    process_references_batch, batch_references, entries, bam, refs_idx
                )
                futures.append(future)  # Store the future

            # Use a while loop to continuously check for completed futures
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

            # Use as_completed to iterate over completed futures as they become available
            for completed_future in concurrent.futures.as_completed(futures):
                alns = completed_future.result()
                write_to_file(alns=alns, out_bam_file=out_bam_file, header=header)

                # Update the progress bar for each completed write
                completion_progress_bar.update(1)
                completed_count += 1
                completed_future.cancel()  # Cancel the future to free memory
                gc.collect()  # Force garbage collection

            completion_progress_bar.close()
    out_bam_file.close()
    entries = None
    gc.collect()
    # prof.disable()
    # # print profiling output
    # stats = pstats.Stats(prof).strip_dirs().sort_stats("tottime")
    # stats.print_stats(5)  # top 10 rows
    log.info("::: ::: Sorting BAM file...")
    if sort_by_name:
        log.info("::: ::: Sorting by name...")
        pysam.sort(
            "-n",
            "-@",
            str(threads),
            "-m",
            str(sort_memory),
            "-o",
            out_bam,
            out_files["bam_reassigned_tmp"],
        )
    else:
        pysam.sort(
            "-@",
            str(threads),
            "-m",
            str(sort_memory),
            "-o",
            out_bam,
            out_files["bam_reassigned_tmp"],
        )

    logging.info("BAM index not found. Indexing...")
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(out_bam, "rb", threads=threads)
    chr_lengths = []
    for chrom in samfile.references:
        chr_lengths.append(samfile.get_reference_length(chrom))
    max_chr_length = np.max(chr_lengths)
    pysam.set_verbosity(save)
    samfile.close()

    if max_chr_length > 536870912:
        logging.info("A reference is longer than 2^29, indexing with csi")
        pysam.index(
            "-c",
            "-@",
            str(threads),
            out_bam,
        )
    else:
        pysam.index(
            "-@",
            str(threads),
            out_bam,
        )

    os.remove(out_files["bam_reassigned_tmp"])


def calculate_alignment_score(identity, read_length):
    return (identity / math.log(read_length)) * math.sqrt(read_length)


def get_bam_data(parms, ref_lengths=None, percid=90, min_read_length=30, threads=1):
    bam, references = parms
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)

    results = []
    reads = set()
    refs = set()
    empty_df = 0

    for reference in references:
        if ref_lengths is None:
            reference_length = int(samfile.get_reference_length(reference))
        else:
            reference_length = int(ref_lengths[reference])
        aln_data = []
        for aln in samfile.fetch(
            contig=reference, multiple_iterators=False, until_eof=True
        ):
            # AS_value = aln.get_tag("AS") if aln.has_tag("AS") else None
            # XS_value = aln.get_tag("XS") if aln.has_tag("XS") else None
            query_length = (
                aln.query_length if aln.query_length != 0 else aln.infer_query_length()
            )

            if query_length < min_read_length:
                continue

            pident = (1 - ((aln.get_tag("NM") / query_length))) * 100
            if pident >= percid:
                var = calculate_alignment_score(pident, query_length)
                aln_data.append(
                    (
                        aln.query_name,
                        aln.reference_name,
                        var,
                        reference_length,
                    )
                )
                reads.add(aln.query_name)
                refs.add(aln.reference_name)
        # remove duplicates
        # check if aln_data is not empty
        if len(aln_data) > 0:
            aln_data_dt = dt.Frame(aln_data)
            # "queryId", "subjectId", "bitScore", "slen"
            aln_data_dt.names = ["queryId", "subjectId", "var", "slen"]
            # remove duplicates and keep the ones with the largest mapq
            aln_data_dt = aln_data_dt[
                :1, :, dt.by(dt.f.queryId, dt.f.subjectId), dt.sort(-dt.f.var)
            ]
            # aln_data_dt["slen"] = dt.Frame([ref_lengths[x] for x in aln_data_dt[:,"subjectId"].to_list()[0]])
            results.append(aln_data_dt)
        else:
            results.append(dt.Frame())
            empty_df += 1
    samfile.close()
    return (dt.rbind(results), reads, refs, empty_df)


def reassign_reads(
    bam,
    out_files,
    reference_lengths=None,
    threads=1,
    min_read_count=1,
    min_read_ani=90,
    min_read_length=30,
    reassign_iters=25,
    reassign_scale=0.9,
    sort_memory="4G",
    sort_by_name=False,
):
    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True
    if threads > 1:
        dt.options.nthreads = threads - 1
    else:
        dt.options.nthreads = 1

    log.info("::: Loading BAM file")
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)

    references = samfile.references

    pysam.set_verbosity(save)

    if reference_lengths is not None:
        ref_len_dt = dt.fread(reference_lengths)
        ref_len_dt.names = ["subjectId", "slen"]
        # convert to dict
        ref_len_dict = dict(
            zip(ref_len_dt["subjectId"].to_list()[0], ref_len_dt["slen"].to_list()[0])
        )
        # check if the dataframe contains all the References in the BAM file
        if not set(references).issubset(set(ref_len_dict.keys())):
            logging.error(
                "The BAM file contains references not found in the reference lengths file"
            )
            sys.exit(1)
    else:
        ref_len_dict = None

    total_refs = samfile.nreferences
    log.info(f"::: Found {total_refs:,} reference sequences")
    # logging.info(f"Found {samfile.mapped:,} alignments")
    log.info(f"::: Removing references with less than {min_read_count} reads...")
    references_m = {
        chrom.contig: chrom.mapped
        for chrom in samfile.get_index_statistics()
        if chrom.mapped >= min_read_count
    }

    references = list(references_m.keys())

    if len(references) == 0:
        log.warning("::: No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    log.info(f"::: Keeping {len(references):,} references")

    log.info("::: Creating reference chunks with uniform read amounts...")
    # ify the number of chunks
    ref_chunks = sort_keys_by_approx_weight(
        input_dict=references_m, scale=1, num_cores=threads
    )
    log.info(f"::: Created {len(ref_chunks):,} chunks")
    ref_chunks = random.sample(ref_chunks, len(ref_chunks))

    dt.options.progress.enabled = False
    dt.options.progress.clear_on_success = True
    dt.options.nthreads = 1

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
                        threads=4,
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
        p = Pool(
            threads,
        )
        data = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(
                        get_bam_data,
                        ref_lengths=ref_len_dict,
                        percid=min_read_ani,
                        threads=4,
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
    reads = list()
    refs = list()
    empty_df = 0
    # new_data = list()
    # for i in tqdm.tqdm(range(len(data)), total=len(data), leave=False, ncols=80):
    #     empty_df += data[i][3]
    #     reads.extend(list(data[i][1]))
    #     refs.extend(list(data[i][2]))
    #     data[i] = data[i][0]
    new_data = list()
    for i in tqdm.tqdm(range(len(data)), total=len(data), leave=False, ncols=80):
        empty_df += data[i][3]
        reads.extend(list(data[i][1]))
        refs.extend(list(data[i][2]))

        # Check if the frame has more than 2 billion rows
        if data[i][0].nrows > 2e9:
            # Calculate the number of chunks needed
            num_chunks = (data[i][0].nrows // 1e9) + (data[i][0].nrows % 1e9 > 0)
            log.warning(f"Frame has more than 2 billion rows. Splitting into {num_chunks:,} chunks...")

            chunks = []

            # Create chunks of 1 billion rows each
            for chunk_idx in range(num_chunks):
                start = int(chunk_idx * 1e9)
                end = int(min((chunk_idx + 1) * 1e9, data[i][0].nrows))
                chunks.append(data[i][0][start:end, :])

            # Substitute data[i] with the first chunk
            data[i] = chunks[0]

            # Append the rest of the chunks to new_data
            for chunk in chunks[1:]:
                new_data.append(chunk)
        else:
            # If the frame is not larger than 2 billion rows, keep it as is
            data[i] = data[i][0]

    # Combine data with new_data
    data.extend(new_data)
    log.info(f"::: ::: Removed {empty_df:,} references without alignments")

    # data = dt.rbind([x for x in data])

    # log.info("::: Indexing references...")
    # refs = dt.Frame(list(set(refs)))
    # refs.names = ["subjectId"]
    # refs["sidx"] = dt.Frame(list(range(refs.shape[0])))
    # refs.key = "subjectId"

    # log.info("::: Indexing reads...")
    # reads = dt.Frame(list(set(reads)))
    # reads.names = ["queryId"]
    # reads["qidx"] = dt.Frame([(i + refs.shape[0]) for i in range(reads.shape[0])])
    # reads.key = "queryId"

    # log.info("::: Combining data...")
    # data = data[:, :, dt.join(reads)]
    # data = data[:, :, dt.join(refs)]

    # Initialize refs DataFrame
    log.info("::: Indexing references...")
    refs = dt.Frame(list(set(refs)))
    refs.names = ["subjectId"]
    refs["sidx"] = dt.Frame(list(range(refs.shape[0])))
    refs.key = "subjectId"

    # Initialize reads DataFrame
    log.info("::: Indexing reads...")
    reads = dt.Frame(list(set(reads)))
    reads.names = ["queryId"]
    reads["qidx"] = dt.Frame([(i + refs.shape[0]) for i in range(reads.shape[0])])
    reads.key = "queryId"
    n_alns_0 = 0
    # Loop through each DataFrame in the list and update it with the joined version
    log.info("::: Combining data...")
    for i, x in tqdm.tqdm(
        enumerate(data),
        total=len(data),
        desc="Processing batches",
        unit="batch",
        disable=is_debug(),
        leave=False,
        ncols=80,
    ):
        # Perform join with reads and then refs
        x = x[:, :, dt.join(reads)]
        x = x[:, :, dt.join(refs)]
        n_alns_0 += x.shape[0]
        del x["queryId"]
        del x["subjectId"]
        x = x[:, [dt.f.qidx, dt.f.sidx, dt.f.var, dt.f.slen]].to_numpy()
        # Substitute the original DataFrame with the joined version in the list
        data[i] = x

    # After the loop, use dt.rbind() to combine all the DataFrames in the list
    # data = dt.rbind([x for x in data])

    # del data["queryId"]
    # del data["subjectId"]
    # n_alns_0 = data.shape[0]
    n_reads_0 = reads.shape[0]
    n_refs_0 = refs.shape[0]
    log.info(
        f"::: References: {n_refs_0:,} | Reads: {n_reads_0:,} | Alignments: {n_alns_0:,}"
    )

    log.info("::: Allocating data structures...")
    # mat = data[:, [dt.f.qidx, dt.f.sidx, dt.f.var, dt.f.slen]].to_numpy()
    mat = np.vstack(data)
    data = None

    # Create a zeros array with the same number of rows as 'mat'
    zeros_array = np.zeros((mat.shape[0], 5))

    # Stack the zeros_array with the original 'mat'
    m = np.column_stack([mat, zeros_array])
    zeros_array = None

    dtype = np.dtype(
        [
            ("source", "int"),
            ("subject", "int"),
            ("var", "float"),
            ("slen", "int"),
            ("s_W", "float"),
            ("prob", "float"),
            ("iter", "int"),
            ("n_aln", "int"),
            ("max_prob", "float"),
        ]
    )

    # Convert the unstructured array to structured array
    m = rf.unstructured_to_structured(m, dtype)
    gc.collect()
    log.info("::: Initializing data structures...")
    init_data = initialize_subject_weights(m)
    if reassign_iters > 0:
        log.info(f"::: Reassigning reads with {reassign_iters} iterations")
    else:
        log.info("::: Reassigning reads until convergence")
    no_multimaps = resolve_multimaps(
        init_data, iters=reassign_iters, scale=reassign_scale
    )

    n_reads = len(list(set(no_multimaps["source"])))
    n_refs = len(list(set(no_multimaps["subject"])))
    n_alns = no_multimaps.shape[0]
    log.info(
        f"::: References: {n_refs:,} | Reads: {n_reads:,} | Alignments: {n_alns:,}"
    )
    log.info(
        f'::: Unique mapping reads: {no_multimaps[no_multimaps["n_aln"] == 1].shape[0]:,} | Multimapping reads: {no_multimaps[no_multimaps["n_aln"] > 1].shape[0]:,}'
    )

    # add this to the array
    # no_multimaps["n_aln"] = subject_counts_array[no_multimaps["subject"]]

    # log.info(f"::: Removing references with less than {min_read_count} reads...")
    # no_multimaps = no_multimaps[no_multimaps["n_aln"] >= min_read_count]
    # log.info(f"{no_multimaps.shape[0]:,} alignments left")
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
    # count how many alignments are in each subjectId
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
    # convert columns queryId from q and subjectId from s to a tuple
    log.info("::: Creating filtered set...")
    entries = defaultdict(set)
    q_query_ids = q[:, "queryId"].to_list()[0]
    s_subject_ids = s[:, "subjectId"].to_list()[0]

    for query_id, subject_id in zip(q_query_ids, s_subject_ids):
        if subject_id in references_m:
            entries[subject_id].add((query_id, subject_id))
    no_multimaps = None
    q = None
    s = None
    q_query_ids = None
    s_subject_ids = None
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
    )


def reassign(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%H:%M:%S",
    )

    args = get_arguments()
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
        reassign_iters=args.reassign_iters,
        reassign_scale=args.reassign_scale,
        sort_memory=args.sort_memory,
        sort_by_name=args.sort_by_name,
        out_files=out_files,
    )
    log.info("Done!")
