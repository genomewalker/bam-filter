import pysam
import os
import sys
import logging
import gc
import shutil
import tqdm
from multiprocessing import Process, Queue, Pool, Manager
from bam_filter.utils import (
    is_debug,
    create_empty_output_files,
    allocate_threads,
    sort_keys_by_approx_weight,
)
import concurrent.futures
import functools
import pandas as pd
import numpy as np
import datatable as dt

from bam_filter.stats_utils import initializer_lengths

log = logging.getLogger("my_logger")


def write_bam(bam, references, output_files, threads=1, sort_memory="1G"):
    logging.info("::: Writing temporary filtered BAM file... (be patient)")
    p_threads, s_threads = allocate_threads(threads, 2, 4)
    samfile = pysam.AlignmentFile(bam, "rb", threads=s_threads)
    header = samfile.header

    # Convert the dictionary to an array
    refs_dict = dict(zip(samfile.references, samfile.lengths))
    ref_names = [x for x in samfile.references if x in references]
    ref_lengths = [refs_dict[x] for x in ref_names]

    # Convert reference lengths to integers
    ref_lengths = [int(x) for x in ref_lengths]
    refs_idx = {sys.intern(str(x)): i for i, x in enumerate(ref_names)}

    write_threads = min(s_threads, 4)

    new_header = header.to_dict()
    new_header["SQ"] = [x for x in new_header["SQ"] if x["SN"] in ref_names]
    new_header["SQ"].sort(key=lambda x: ref_names.index(x["SN"]))
    # Change it to unsorted
    new_header["HD"]["SO"] = "unsorted"

    out_bam_file = pysam.AlignmentFile(
        output_files["bam_tmp"],
        "wb",
        referencenames=ref_names,
        referencelengths=ref_lengths,
        threads=write_threads,
        header=new_header,
    )

    logging.info(f"::: ::: Filtering {len(references):,} references sequentially...")
    for reference in tqdm.tqdm(
        references,
        total=len(references),
        leave=False,
        ncols=80,
        desc="References processed",
    ):
        for aln in samfile.fetch(
            reference=reference, multiple_iterators=False, until_eof=True
        ):
            aln.reference_id = refs_idx[aln.reference_name]
            out_bam_file.write(aln)
    out_bam_file.close()
    samfile.close()

    logging.info("::: ::: Sorting BAM file...")
    pysam.sort(
        "-@",
        str(s_threads),
        "-m",
        str(sort_memory),
        "-o",
        output_files["bam_tmp_sorted"],
        output_files["bam_tmp"],
    )

    logging.info("::: ::: BAM index not found. Indexing...")
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(
        output_files["bam_tmp_sorted"], "rb", threads=s_threads
    )
    pysam.set_verbosity(save)
    samfile.close()

    pysam.index(
        "-c",
        "-@",
        str(s_threads),
        output_files["bam_tmp_sorted"],
    )

    os.remove(output_files["bam_tmp"])
    return output_files["bam_tmp_sorted"]


def write_to_file(alns, out_bam_file, header):
    """Write alignments to BAM file with proper header"""
    for aln in alns:
        out_bam_file.write(pysam.AlignedSegment.fromstring(aln, header))


def process_references_batch(
    references, bam, refs_idx, min_read_ani, threads=1, batch_size=10000
):
    """Process a batch of references, returning multiple batches if needed"""
    all_batches = []
    current_batch = []
    total_size = 0
    max_batch_memory = 500 * 1024 * 1024  # 500MB batch size limit

    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
        for reference in references:
            for aln in samfile.fetch(reference=reference, multiple_iterators=False):
                ani_read = (1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100
                if ani_read >= min_read_ani:
                    aln.reference_id = refs_idx[aln.reference_name]
                    aln_str = aln.to_string()
                    current_batch.append(aln_str)
                    total_size += len(aln_str)

                    # If we hit our limits, add the current batch to all_batches and start a new one
                    if (
                        total_size >= max_batch_memory
                        or len(current_batch) >= batch_size
                    ):
                        all_batches.append(current_batch)
                        current_batch = []
                        total_size = 0

        # Don't forget the last batch if it exists
        if current_batch:
            all_batches.append(current_batch)

    return all_batches if all_batches else [[]]  # Return at least an empty batch


def create_output_bam(output_file, params, threads):
    """Create output BAM file with proper header configuration"""
    refs_idx, header_dict, _ = params
    ref_names = sorted(refs_idx.keys(), key=lambda x: refs_idx[x])
    ref_lengths = [entry["LN"] for entry in header_dict["SQ"]]

    return pysam.AlignmentFile(
        output_file,
        "wb",
        referencenames=ref_names,
        referencelengths=ref_lengths,
        header=header_dict,
        threads=threads,
    )


def check_bam_file(
    bam, threads=1, reference_lengths=None, sort_memory="1G", sorted_bam=None
):
    logging.info("Checking BAM file status")
    save = pysam.set_verbosity(0)

    def evaluate_bam(bam_filename, s_threads, sorted_bam=None):
        with pysam.AlignmentFile(bam_filename, "rb", threads=s_threads) as samfile:
            nreferences = samfile.nreferences
            log.info(f"::: Found {nreferences:,} reference sequences")

            if samfile.header["HD"]["SO"] != "coordinate":
                log.info("::: BAM file is not sorted by coordinates, sorting it...")
                if sorted_bam is None:
                    sorted_bam = bam_filename.replace(".bam", ".bf-sorted.bam")
                pysam.sort(
                    "-@",
                    str(s_threads),
                    "-m",
                    str(sort_memory),
                    "-o",
                    sorted_bam,
                    bam_filename,
                )
                bam_filename = sorted_bam
                pysam.index("-c", "-@", str(threads), bam_filename)
                return bam_filename, True

            if not samfile.has_index():
                logging.info("::: BAM index not found. Indexing...")
                pysam.index("-c", "-@", str(threads), bam_filename)

            logging.info("::: BAM file looks good.")

        gc.collect()
        return bam_filename, False

    try:
        s_threads = min(threads, 4)
        bam, reopened = evaluate_bam(bam, s_threads, sorted_bam)

        if reopened:
            bam, _ = evaluate_bam(bam, s_threads, sorted_bam)

        pysam.set_verbosity(save)
        return bam

    except ValueError as ve:
        if "file has no sequences defined (mode='rb')" in str(ve):
            return None
        else:
            log.error(f"Message:\n{ve}")
            exit(1)


def sort_and_index_bam(out_files, sort_memory, sort_by_name, s_threads, threads):
    """Sort and index the filtered BAM file"""
    try:
        logging.info("Sorting BAM file...")
        if sort_by_name:
            logging.info("Sorting BAM file by read name...")
            pysam.sort(
                "-n",
                "-@",
                str(s_threads),
                "-m",
                str(sort_memory),
                "-o",
                out_files["bam_filtered"],
                out_files["bam_filtered_tmp"],
            )
        else:
            logging.info("Sorting BAM file by coordinates...")
            pysam.sort(
                "-@",
                str(s_threads),
                "-m",
                str(sort_memory),
                "-o",
                out_files["bam_filtered"],
                out_files["bam_filtered_tmp"],
            )

            logging.info("Indexing BAM file...")
            pysam.index("-c", "-@", str(s_threads), out_files["bam_filtered"])

        os.remove(out_files["bam_filtered_tmp"])

    except Exception as e:
        logging.error(f"Error during sort/index: {str(e)}")
        raise


def process_bam(
    bam,
    threads=1,
    reference_lengths=None,
    min_read_ani=90.0,
    min_read_count=10,
    trim_ends=0,
    trim_min=10,
    trim_max=90,
    scale=1e6,
    plot=False,
    plots_dir="coverage-plots",
    read_length_freqs=False,
    output_files=None,
    low_memory=False,
    sort_memory="1G",
):
    """
    Processing function: calls pool of worker functions to extract metrics from a BAM file
    """
    from bam_filter.stats_utils import (
        get_bam_stats,
    )  # Move import here to avoid circular import

    logging.info("Loading BAM file")
    save = pysam.set_verbosity(0)

    p_threads, s_threads = allocate_threads(threads, 1, 4)

    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True
    if threads > 1:
        dt.options.nthreads = p_threads
    else:
        dt.options.nthreads = 1

    log.info(f"::: IO Threads: {s_threads} | Processing Threads: {p_threads}")

    with pysam.AlignmentFile(bam, "rb", threads=s_threads) as samfile:
        references = samfile.references
        pysam.set_verbosity(save)
        ref_lengths = None

        bam_reference_lengths = {
            reference: np.int64(samfile.get_reference_length(reference))
            for reference in references
        }

        if reference_lengths is not None:
            ref_len_dt = dt.fread(reference_lengths)
            ref_len_dt.names = ["subjectId", "slen"]
            ref_lengths = dict(
                zip(
                    ref_len_dt["subjectId"].to_list()[0],
                    map(int, ref_len_dt["slen"].to_list()[0]),
                )
            )
            del ref_len_dt
            if not set(references).issubset(set(ref_lengths.keys())):
                logging.error(
                    "The BAM file contains references not found in the reference lengths file"
                )
                sys.exit(1)
        else:
            ref_lengths = bam_reference_lengths

        total_refs = samfile.nreferences
        logging.info(f"Found {total_refs:,} reference sequences")

        if plot:
            if not os.path.exists(plots_dir):
                os.makedirs(plots_dir)
        logging.info(f"Removing references with less than {min_read_count} reads...")
        index_statistics = samfile.get_index_statistics()

        references_m = {
            chrom.contig: chrom.mapped
            for chrom in tqdm.tqdm(
                index_statistics,
                desc="Filtering references",
                unit="chrom",
                leave=False,
                ncols=80,
                unit_scale=True,
                unit_divisor=1000,
            )
            if chrom.mapped >= min_read_count
        }

        references = list(references_m.keys())

        if len(references) == 0:
            logging.warning(
                "No reference sequences with alignments found in the BAM file"
            )
            create_empty_output_files(output_files)
            sys.exit(0)

        ref_positions = {ref: samfile.get_tid(ref) for ref in references}

        logging.info(f"Keeping {len(references):,} references")

        if low_memory:
            logging.info("Low memory mode enabled, writing filtered BAM file to disk")
            samfile.close()
            bam = write_bam(
                bam=bam,
                references=references,
                output_files=output_files,
                sort_memory=sort_memory,
                threads=threads,
            )

    log.info("::: Creating reference chunks with uniform read amounts...")
    ref_chunks = sort_keys_by_approx_weight(
        input_dict=references_m,
        ref_positions=ref_positions,
        scale=1,
        num_cores=p_threads,
        verbose=False,
        max_entries_per_chunk=100_000_000,
    )
    log.info(f"::: Created {len(ref_chunks):,} chunks")

    params = zip([bam] * len(ref_chunks), ref_chunks)
    ref_lengths = {ref: ref_lengths[ref] for ref in references}

    bam_reference_lengths = {ref: bam_reference_lengths[ref] for ref in references}

    try:
        logging.info(f"Processing {len(ref_chunks):,} chunks")
        # logging.info("::: Distributing shared objects...")
        # manager = Manager()

        if is_debug():
            data = list(
                map(
                    functools.partial(
                        get_bam_stats,
                        # ref_lengths=filtered_ref_lengths,
                        # bam_reference_lengths=filtered_bam_ref_lengths,
                        min_read_ani=min_read_ani,
                        trim_ends=0,
                        trim_min=trim_min,
                        trim_max=trim_max,
                        scale=scale,
                        plot=plot,
                        plots_dir=plots_dir,
                        read_length_freqs=read_length_freqs,
                        threads=s_threads,
                    ),
                    params,
                )
            )
        else:
            with Pool(
                p_threads,
                initializer_lengths,
                (
                    ref_lengths,
                    bam_reference_lengths,
                ),
            ) as p:
                data = list(
                    tqdm.tqdm(
                        p.imap_unordered(
                            functools.partial(
                                get_bam_stats,
                                # ref_lengths=filtered_ref_lengths,
                                # bam_reference_lengths=filtered_bam_ref_lengths,
                                min_read_ani=min_read_ani,
                                trim_ends=0,
                                trim_min=trim_min,
                                trim_max=trim_max,
                                scale=scale,
                                plot=plot,
                                plots_dir=plots_dir,
                                read_length_freqs=read_length_freqs,
                                threads=s_threads,
                            ),
                            params,
                            chunksize=1,
                        ),
                        total=len(ref_chunks),
                        leave=False,
                        ncols=80,
                        desc="References processed",
                    )
                )

    except KeyboardInterrupt:
        logging.info("::: ::: User interrupted execution")
        p.terminate()
        p.join()
        sys.exit(1)

    return data


def filter_reference_BAM(
    bam,
    df,
    filter_conditions,
    threads,
    out_files,
    sort_memory="2G",
    min_read_ani=90,
    transform_cov_evenness=False,
    sort_by_name=False,
    disable_sort=False,
    write_bam=False,
    chunk_size=1_000_000,
):
    """Filter BAM files with optimized memory usage and parallel processing"""

    logging.info("Filtering stats...")

    dt.options.nthreads = threads
    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True

    # Store original columns order
    original_columns = df.columns.tolist()

    # Apply filtering conditions
    df_filtered = apply_filtering_conditions(
        df, filter_conditions, transform_cov_evenness
    )

    # Ensure filtered DataFrame has the same columns in the same order
    df_filtered = df_filtered[original_columns]

    # Write filtered stats to file with same format
    logging.info(f"Writing filtered statistics to {out_files['stats_filtered']}")

    dt.Frame(df_filtered).to_csv(out_files["stats_filtered"], sep="\t", header=True)

    if len(df_filtered.index) == 0:
        logging.info("No references meet the filter conditions.")
        create_empty_output_files(out_files)
        return

    if not write_bam:
        logging.info("Skipping BAM file creation...")
        return

    logging.info("Writing filtered BAM file...")

    # Initialize BAM processing
    references = df_filtered["reference"].tolist()
    p_threads, s_threads = allocate_threads(threads, 1, 4)

    # Setup BAM processing parameters
    processing_params = setup_bam_processing(
        bam, references, s_threads, chunk_size, out_files
    )

    if not processing_params:
        return

    refs_idx, header, ref_chunks = processing_params

    # Create output BAM file
    out_bam = create_output_bam(
        out_files["bam_filtered_tmp"], processing_params, s_threads
    )

    # Process in parallel with memory management
    process_bam_parallel(
        bam,
        ref_chunks,
        refs_idx,
        header,
        out_bam,
        p_threads,
        s_threads,
        min_read_ani,
        chunk_size,
    )

    # Sort and index if needed
    if not disable_sort:
        sort_and_index_bam(out_files, sort_memory, sort_by_name, s_threads, threads)
    else:
        shutil.move(out_files["bam_filtered_tmp"], out_files["bam_filtered"])


def apply_filtering_conditions(df, filter_conditions, transform_cov_evenness=False):
    """Apply filtering conditions to the dataframe using the original logic"""
    if "min_norm_entropy" in filter_conditions and "min_norm_gini" in filter_conditions:
        logging.info(
            f"::: min_read_count >= {filter_conditions['min_read_count']} "
            f"& min_read_length >= {filter_conditions['min_read_length']} "
            f"& max_read_length <= {filter_conditions['max_read_length']}"
            f"& min_avg_read_ani >= {filter_conditions['min_avg_read_ani']} "
            f"& min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} "
            f"& min_breadth >= {filter_conditions['min_breadth']} "
            f"& min_coverage_evenness >= {filter_conditions['min_coverage_evenness']} "
            f"& min_coeff_var =< {filter_conditions['min_coeff_var']} "
            f"& min_coverage_mean >= {filter_conditions['min_coverage_mean']} "
            f"& min_norm_entropy >= {filter_conditions['min_norm_entropy']} "
            f"& min_norm_gini <= {filter_conditions['min_norm_gini']}"
        )
        if transform_cov_evenness:
            df["cov_evenness_tmp"] = df["cov_evenness"]
            df["cov_evenness_tmp"] = np.where(
                np.rint(df.coverage_mean) < 1.0, 1.0, df.cov_evenness_tmp
            )
        else:
            df["cov_evenness_tmp"] = df["cov_evenness"]
        df_filtered = df.loc[
            (df["n_reads"] >= filter_conditions["min_read_count"])
            & (df["read_length_mean"] >= filter_conditions["min_read_length"])
            & (df["read_length_mean"] <= filter_conditions["max_read_length"])
            & (df["read_ani_mean"] >= filter_conditions["min_avg_read_ani"])
            & (
                df["breadth_exp_ratio"]
                >= filter_conditions["min_expected_breadth_ratio"]
            )
            & (df["breadth"] >= filter_conditions["min_breadth"])
            & (df["cov_evenness_tmp"] >= filter_conditions["min_coverage_evenness"])
            & (df["c_v"] <= filter_conditions["min_coeff_var"])
            & (df["coverage_mean"] >= filter_conditions["min_coverage_mean"])
            & (df["norm_entropy"] >= filter_conditions["min_norm_entropy"])
            & (df["norm_gini"] <= filter_conditions["min_norm_gini"])
        ]
    else:
        logging.info(
            f"::: min_read_count >= {filter_conditions['min_read_count']} "
            f"& min_read_length >= {filter_conditions['min_read_length']} "
            f"& max_read_length <= {filter_conditions['max_read_length']}"
            f"& min_avg_read_ani >= {filter_conditions['min_avg_read_ani']} "
            f"& min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} "
            f"& min_breadth >= {filter_conditions['min_breadth']} "
            f"& min_coverage_evenness >= {filter_conditions['min_coverage_evenness']} "
            f"& min_coeff_var <= {filter_conditions['min_coeff_var']} "
            f"& min_coverage_mean >= {filter_conditions['min_coverage_mean']}"
        )
        if transform_cov_evenness:
            df["cov_evenness_tmp"] = df["cov_evenness"]
            df["cov_evenness_tmp"] = np.where(
                np.rint(df.coverage_mean) < 1.0, 1.0, df.cov_evenness_tmp
            )
        else:
            df["cov_evenness_tmp"] = df["cov_evenness"]

        df_filtered = df.loc[
            (df["n_reads"] >= filter_conditions["min_read_count"])
            & (df["read_length_mean"] >= filter_conditions["min_read_length"])
            & (df["read_length_mean"] <= filter_conditions["max_read_length"])
            & (df["read_ani_mean"] >= filter_conditions["min_avg_read_ani"])
            & (
                df["breadth_exp_ratio"]
                >= filter_conditions["min_expected_breadth_ratio"]
            )
            & (df["breadth"] >= filter_conditions["min_breadth"])
            & (df["cov_evenness_tmp"] >= filter_conditions["min_coverage_evenness"])
            & (df["c_v"] <= filter_conditions["min_coeff_var"])
            & (df["coverage_mean"] >= filter_conditions["min_coverage_mean"])
        ]

    if "cov_evenness_tmp" in df_filtered.columns:
        del df_filtered["cov_evenness_tmp"]

    return df_filtered


def setup_bam_processing(bam, references, threads, chunk_size, out_files):
    """Setup BAM processing parameters"""
    try:
        log.info("::: Setting up BAM processing...")
        with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
            # Create set for fast lookups
            filtered_refs = set(references)  # O(1) lookups
            log.info(f"::: ::: Processing {len(filtered_refs):,} filtered references")

            # Get reference information - use set membership
            refs_dict = {}
            ref_dict_m = {}

            # Single pass through index stats instead of multiple lookups
            log.info("::: ::: Getting reference mapping statistics...")
            for chrom in samfile.get_index_statistics():
                if chrom.contig in filtered_refs:
                    refs_dict[chrom.contig] = samfile.get_reference_length(chrom.contig)
                    ref_dict_m[chrom.contig] = chrom.mapped

            log.info(f"::: ::: Processed {len(refs_dict):,} references")

            # Process header efficiently
            log.info("::: ::: Processing BAM header...")
            header_dict = {"HD": {"VN": "1.0", "SO": "unsorted"}}

            # Pre-sort references - this will be our final order
            ref_names = sorted(refs_dict.keys())

            # Build SQ entries directly in final order
            header_dict["SQ"] = [{"SN": ref, "LN": refs_dict[ref]} for ref in ref_names]

            log.info(
                f"::: ::: Created header with {len(header_dict['SQ']):,} references"
            )

            # Reference index matches our pre-sorted order
            refs_idx = {sys.intern(str(ref)): idx for idx, ref in enumerate(ref_names)}

            log.info("::: ::: Creating balanced reference chunks...")
            ref_chunks = sort_keys_by_approx_weight(
                input_dict=ref_dict_m,
                ref_positions=refs_idx,
                scale=1,
                num_cores=threads,
                verbose=False,
                max_entries_per_chunk=100_000_000,
            )
            log.info(f"::: ::: Created {len(ref_chunks):,} processing chunks")

            log.info("::: BAM setup completed successfully")
            return refs_idx, header_dict, ref_chunks

    except Exception as e:
        logging.error(f"Error setting up BAM processing: {str(e)}")
        return None


# def create_balanced_chunks(ref_dict, chunk_size, num_threads):
#     """Create balanced chunks based on reference sizes"""
#     total_size = sum(ref_dict.values())
#     chunk_target = max(chunk_size, total_size // (num_threads * 2))

#     chunks = []
#     current_chunk = []
#     current_size = 0

#     for ref, size in sorted(ref_dict.items(), key=lambda x: x[1], reverse=True):
#         if current_size + size > chunk_target and current_chunk:
#             chunks.append(current_chunk)
#             current_chunk = []
#             current_size = 0

#         current_chunk.append(ref)
#         current_size += size

#     if current_chunk:
#         chunks.append(current_chunk)

#     return chunks


def process_bam_parallel(
    bam,
    ref_chunks,
    refs_idx,
    header,
    out_bam,
    p_threads,
    s_threads,
    min_read_ani,
    chunk_size,
):
    """Process BAM file in parallel with memory management"""
    logging.info("::: Dumping alignments...")
    # Create proper pysam header
    with pysam.AlignmentFile(bam, "rb") as samfile:
        pysam_header = samfile.header

    # Create memory-managed queue
    result_queue = Queue(maxsize=p_threads * 2)

    # Start writer process
    writer = Process(
        target=bam_writer_process, args=(result_queue, out_bam, pysam_header)
    )
    writer.start()

    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=p_threads) as executor:
            futures = [
                executor.submit(
                    process_references_batch,
                    chunk,
                    bam,
                    refs_idx,
                    min_read_ani,
                    s_threads,
                    chunk_size,
                )
                for chunk in ref_chunks
            ]

            for future in concurrent.futures.as_completed(futures):
                try:
                    batch = future.result()
                    if batch:
                        result_queue.put(batch)
                except Exception as e:
                    logging.error(f"Error processing chunk: {str(e)}")
                finally:
                    gc.collect()

    finally:
        result_queue.put(None)
        writer.join()


def bam_writer_process(queue, out_bam_file, header):
    """Dedicated process for writing BAM output"""
    try:
        while True:
            batches = queue.get()
            if batches is None:
                break

            # Process all batches received
            for batch in batches:
                write_to_file(batch, out_bam_file, header)
                del batch
            gc.collect()
    except Exception as e:
        logging.error(f"Error in writer process: {str(e)}")
    finally:
        out_bam_file.close()
