import pysam
import taxopy as txp
from tqdm import tqdm
from multiprocessing import Pool, Manager
from functools import partial
import pandas as pd
import logging
import networkx as nx
import sys
from bam_filter.utils import (
    calc_chunksize,
    sort_keys_by_approx_weight,
    concat_df,
    is_debug,
    create_output_files,
    create_empty_output_files,
)
from bam_filter.sam_utils import check_bam_file
from collections import defaultdict
from functools import reduce
import operator
import random
import gzip
import datatable as dt

log = logging.getLogger("my_logger")

debug = is_debug()


def calculate_path_likelihood(path, graph):
    likelihood = 1.0
    for u, v in zip(path[:-1], path[1:]):
        likelihood *= graph[u][v]["cum_weight"]
    return likelihood


def find_most_likely_continuation_worker(partial_path_end, full_graph, index):
    # find where are we in the taxonomy tree, if we are in the root, return

    results = []
    try:
        descendants = list(nx.descendants(full_graph, partial_path_end))
    except nx.NetworkXError:
        descendants = []

    if not descendants:
        return (
            partial_path_end,
            {
                "reference": None,
                "best_path": None,
                "top_10_paths": None,
            },
        )

    if (
        len(nx.shortest_path(full_graph, source="root", target=partial_path_end))
        <= index + 1
    ):
        return (
            partial_path_end,
            {
                "reference": None,
                "best_path": None,
                "top_10_paths": None,
            },
        )

    for continuation in descendants:
        full_path = list(
            nx.all_simple_paths(
                full_graph, source=partial_path_end, target=continuation
            )
        )[0]
        likelihood = calculate_path_likelihood(full_path, full_graph)
        results.append((full_path, likelihood))

    results.sort(key=lambda x: x[1], reverse=True)

    return (
        partial_path_end,
        {
            "reference": results[0][0][-1] if results else None,
            "best_path": results[0] if results else None,
            "top_10_paths": results[:10] if results else None,
        },
    )


def find_most_likely_continuation(full_graph, leaves, index):
    result_dict = {}

    for partial_path_end in tqdm(leaves, leave=False, ncols=80):
        partial_path_end, result = find_most_likely_continuation_worker(
            partial_path_end, full_graph, index
        )
        result_dict[partial_path_end] = result

    return result_dict


def calculate_cumulative_weight_and_path(graph, node, lengths):
    # Base case: if the node is the root, return its cumulative weight and path
    if not list(graph.predecessors(node)):
        if node in lengths:
            return 0, 0, [node]
        else:
            # nei = list(graph.neighbors(node))[0]
            # cumulative_weight = graph[node][nei].get("weight", 0)
            # cumulative_norm_weight = graph[node][nei].get("norm_weight", 0)
            # return cumulative_weight, cumulative_norm_weight, [node]
            return 0, 0, [node]

    # Recursive case: calculate cumulative weight and path by summing the edge weight
    # and norm_weight with the cumulative weight and path of its parent(s)
    cumulative_weight = 0
    cumulative_norm_weight = 0
    cumulative_path = []

    for parent in graph.predecessors(node):
        edge_weight = graph[parent][node].get("weight", 0)
        norm_weight = graph[parent][node].get("norm_weight", 0)
        (
            parent_weight,
            parent_norm_weight,
            parent_path,
        ) = calculate_cumulative_weight_and_path(graph, parent, lengths)

        cumulative_weight += edge_weight + parent_weight
        cumulative_norm_weight += norm_weight + parent_norm_weight
        cumulative_path.extend(parent_path + [node])  # Fix: Use += to concatenate lists

    return cumulative_weight, cumulative_norm_weight, cumulative_path


def create_tax_graph_w(tax_path, weight, lengths, ref_stats=None, scale=1_000_000):
    root_row = pd.DataFrame({"source": ["root"], "target": ["root"], "weight": 0})
    res = list(zip(tax_path, tax_path[1:]))
    # get last element
    res = pd.DataFrame(res, columns=["source", "target"])
    res = pd.concat([root_row, res])
    res = res.drop_duplicates()
    res["weight"] = 0
    res["norm_weight"] = 0
    if ref_stats:
        target = res.iloc[-1, res.columns.get_loc("target")]
        if target in ref_stats:
            if weight <= ref_stats[target][0] and weight > ref_stats[target][1]:
                res.iloc[-1, res.columns.get_loc("weight")] = ref_stats[target][1]
                res.iloc[-1, res.columns.get_loc("norm_weight")] = round(
                    (ref_stats[target][1] / ref_stats[target][2]) * scale
                )
        else:
            res.iloc[-1, res.columns.get_loc("weight")] = weight
            res.iloc[-1, res.columns.get_loc("norm_weight")] = round(
                (weight / lengths[res.iloc[-1, res.columns.get_loc("target")]]) * scale
            )
    else:
        res.iloc[-1, res.columns.get_loc("weight")] = weight
        res.iloc[-1, res.columns.get_loc("norm_weight")] = round(
            (weight / lengths[res.iloc[-1, res.columns.get_loc("target")]]) * scale
        )
    return res


def create_lca_df(tax_path, weight):
    root_row = pd.DataFrame({"source": ["root"], "target": ["root"], "weight": 0})
    res = list(zip(tax_path, tax_path[1:]))
    # get last element
    res = pd.DataFrame(res, columns=["source", "target"])
    res = pd.concat([root_row, res])
    res = res.drop_duplicates()
    res["weight"] = 0
    # add 1 to the last row weight
    res.iloc[-1, res.columns.get_loc("weight")] = weight
    return res


def get_ref2read(params, dat, threads=1):
    bam, references = params
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    results = defaultdict(set)
    for reference in references:
        if reference not in dat:
            continue
        for aln in samfile.fetch(
            contig=reference, multiple_iterators=False, until_eof=True
        ):
            results[reference].add(aln.query_name)
    samfile.close()
    return results


def get_tax(ref, parms):
    taxdb = parms["taxdb"]
    acc2taxid = parms["acc2taxid"]
    custom = parms["custom"]
    missing = set()
    if ref in acc2taxid:
        taxid = acc2taxid[ref]
        # taxid = txp.taxid_from_name(ref, taxdb)[0]
        try:
            taxonomy_info = txp.Taxon(taxid, taxdb).rank_name_dictionary
        except txp.exceptions.TaxidError:
            log.debug(f"No taxid found for {ref}")
            taxonomy_info = None
            missing.add(ref)
            return taxonomy_info
        taxonomy_info["taxid"] = taxid
        taxonomy_info["ref"] = ref
        if custom:
            taxonomy_info["subspecies"] = f"S__{ref}"
    else:
        log.debug(f"No taxid found for {ref}")
        taxonomy_info = None
        missing.add(ref)
    return taxonomy_info


def get_taxonomy_info(refids, taxdb, acc2taxid, nprocs=1, custom=False):
    """Function to get the references taxonomic information for a given taxonomy id

    Args:
        taxids (list): A list of taxonomy ids
        taxdb (taxopy.TaxonomyDB): A taxopy DB

    Returns:
        dict: A list of taxonomy information
    """

    # acc2taxid_df = pd.read_csv(acc2taxid, sep="\t", index_col=None, engine="pyarrow")[
    #     ["accession", "taxid"]
    # ].rename(columns={"accession": "reference"}, inplace=False)
    # print("here0")
    # # Filter rows in refids from dataframe
    # acc2taxid_df = acc2taxid_df.loc[acc2taxid_df["reference"].isin(refids)]
    # print("here1")
    # acc2taxid_dict = acc2taxid_df.set_index("reference").T.to_dict("records")

    # parms = {"taxdb": taxdb, "acc2taxid": acc2taxid_dict[0]}
    # Initialize the progress bar without a total; it will update dynamically based on rows
    if custom:
        log.info("Using custom taxdump files...")
    log.info("Processing acc2taxid file...")

    # pbar = tqdm(
    #     desc="Processing Rows ",
    #     unit=" row",
    #     leave=False,
    #     ncols=100,
    # )
    # Initialize an empty DataFrame to hold the filtered results
    if custom:
        cols = ["accession", "taxid"]
        name = "accession"
    else:
        cols = ["accession.version", "taxid"]
        name = "accession.version"
    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True
    dt.options.nthreads = nprocs
    filtered_df = dt.fread(acc2taxid, verbose=False, nthreads=nprocs)
    filt = dt.Frame(refids, names=[name])
    filt["FOO"] = 1
    filt.key = name
    filtered_df = filtered_df[:, :, dt.join(filt)][~dt.isna(dt.f.FOO), :].to_pandas()
    filtered_df = filtered_df[cols]

    filtered_df.rename(columns={name: "reference"}, inplace=True)
    acc2taxid_dict = filtered_df.set_index("reference").to_dict()["taxid"]

    log.info(f"::: Read {filtered_df.shape[0]:,} rows from acc2taxid file...")

    parms = {"taxdb": taxdb, "acc2taxid": acc2taxid_dict, "custom": custom}
    func = partial(get_tax, parms=parms)
    if debug is True or len(refids) < 100000:
        taxonomy_info = list(map(func, refids))
    else:
        p = Pool(nprocs)
        c_size = calc_chunksize(nprocs, len(refids))
        taxonomy_info = list(
            tqdm(
                p.imap_unordered(func, refids, chunksize=c_size),
                total=len(refids),
                leave=False,
                ncols=100,
                desc="References processed",
            )
        )
        p.close()
        p.join()
    taxonomy_info = list(filter(None, taxonomy_info))
    exclude = ["taxid", "ref"]
    tax_ranks = []

    for k in taxonomy_info[0].keys():
        if k not in exclude:
            tax_ranks.append(k)

    taxonomy_info = {i["ref"]: i for i in taxonomy_info}
    return taxonomy_info, tax_ranks


def do_lca(args):
    bam = args.bam
    names = args.names
    nodes = args.nodes
    acc2taxid = args.acc2taxid
    sel_rank = args.rank_lca
    reference_lengths = args.reference_lengths
    threads = args.threads
    scale = args.scale

    out_files = create_output_files(
        bam=bam,
        prefix=args.prefix,
        lca_summary=args.lca_summary,
        tmp_dir=None,
        mode="lca",
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
        sys.exit(1)

    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    references = samfile.references
    references_m = {
        chrom.contig: chrom.mapped for chrom in samfile.get_index_statistics()
    }

    if reference_lengths is not None:
        ref_lengths = pd.read_csv(
            reference_lengths, sep="\t", index_col=False, names=["reference", "length"]
        )
        ref_lengths = dict(zip(ref_lengths["reference"], ref_lengths["length"]))
        # check if the dataframe contains all the References in the BAM file
        if not set(references).issubset(set(ref_lengths.keys())):
            logging.error(
                "The BAM file contains references not found in the reference lengths file"
            )
            sys.exit(1)
    else:
        ref_lengths = {x: samfile.get_reference_length(x) for x in references}

    log.info("Getting taxonomy information")
    taxdb = txp.TaxDb(
        nodes_dmp=nodes,
        names_dmp=names,
    )
    acc2taxid = acc2taxid
    taxonomy_info, tax_ranks = get_taxonomy_info(
        references, taxdb, acc2taxid, nprocs=threads, custom=args.custom
    )

    tax_ranks.reverse()
    index_r = tax_ranks.index(sel_rank)

    dat = {}
    for k, v in taxonomy_info.items():
        try:
            tax_list = ["root"]
            tax_list.extend([v[rank] for rank in tax_ranks])
            tax_list.append(k)
            dat[k] = tuple(tax_list)
        except KeyError:
            continue
    ref_stats = None
    if args.lca_stats:
        log.info("Loading reference stats...")
        ref_stats = pd.read_csv(args.lca_stats, sep="\t", index_col=False)
        ref_stats = ref_stats[
            ["reference", "n_reads", "n_reads_tad", "coverage_mean_trunc_len"]
        ]
        ref_stats = ref_stats[ref_stats["n_reads_tad"] > 0]
        if ref_stats.empty:
            del ref_stats
            ref_stats = None
        else:
            # convert to a named dictionary with reference as key and include n_reads and n_reads_tad
            ref_stats = dict(
                zip(
                    ref_stats["reference"],
                    zip(
                        ref_stats["n_reads"],
                        ref_stats["n_reads_tad"],
                        ref_stats["coverage_mean_trunc_len"],
                    ),
                )
            )
    log.info("Getting reference to read mapping")
    ref_chunks = sort_keys_by_approx_weight(
        references_m,
        scale=1,
        num_cores=threads,
        refinement_steps=100,
        max_entries_per_chunk=1_000_000,
    )

    ref_chunks = random.sample(ref_chunks, len(ref_chunks))
    dat_man = Manager().dict(dat)
    params = zip([bam] * len(ref_chunks), ref_chunks)
    p = Pool(
        threads,
    )
    data = list(
        tqdm(
            p.imap_unordered(
                partial(
                    get_ref2read,
                    dat=dat_man,
                    threads=1,
                ),
                params,
                chunksize=1,
            ),
            total=len(ref_chunks),
            leave=False,
            ncols=80,
            desc="Chunks processed",
        )
    )

    p.close()
    p.join()

    data = reduce(operator.ior, data, {})

    log.info("Getting reads taxonomic information")
    reads = defaultdict(set)
    for k, v in tqdm(data.items(), total=len(data), leave=False, ncols=80):
        tax = dat[k]
        for read in v:
            reads[read].add(tax)

    log.info("Creating taxonomic graph")
    root_row = pd.DataFrame({"source": ["root"], "target": ["root"], "weight": [0]})
    gs = []
    if len(dat) > 0:
        for t in tqdm(dat.values(), total=len(dat), leave=False, ncols=80):
            # remove last element
            # t = t[:-1]
            res = list(zip(t, t[1:]))
            res = pd.DataFrame(res, columns=["source", "target"])
            res = res.drop_duplicates()
            res["weight"] = 0
            res = pd.concat([root_row, res])
            gs.append(res)

    gs = concat_df(gs).drop_duplicates()
    # refine weights if the stats file has been loaded
    # if the weight in df is <= than the one in ref_s, then use the one in ref_stats

    G = nx.from_pandas_edgelist(gs, create_using=nx.DiGraph(), edge_attr=True)

    dat = None

    log.info("Calculating LCA...")
    cum_tax = defaultdict(int)
    unique = []

    discarded_lca = defaultdict(int)
    for_lca = defaultdict(int)

    # for k, v in tqdm(reads.items(), total=len(reads), leave=False, ncols=80):
    #     v = list(v)
    #     if len(v) > 1:
    #         tax_path = [edge for edge in v[0] if all(edge in t for t in v[1:])]
    #         if len(tax_path) <= index_r + 1:
    #             discarded_lca[tuple(tax_path)] += 1
    #             continue
    #         for_lca[tuple(tax_path)] += 1
    #     else:
    #         unique.append(k)
    #         cum_tax[v[0]] += 1

    for k, v in tqdm(reads.items(), total=len(reads), leave=False, ncols=80):
        v = list(v)
        if len(v) > 1:
            tax_path_set = set(v[0])
            for t in v[1:]:
                tax_path_set.intersection_update(t)

            tax_path = [edge for edge in v[0] if edge in tax_path_set]

            if len(tax_path) <= index_r + 1:
                discarded_lca[tuple(tax_path)] += 1
                continue

            for_lca[tuple(tax_path)] += 1
        else:
            if len(k) <= index_r + 1:
                discarded_lca[tuple(k)] += 1
                continue
            unique.append(k)
            cum_tax[v[0]] += 1

    log.info(
        f"Unique: {len(unique):,} | LCA: {sum(for_lca.values()):,} | Discarded: {sum(discarded_lca.values()):,}"
    )

    log.info("Creating unique mapping taxonomic graph")
    if ref_stats:
        log.info("::: Using TAD inferred reads from the stats file...")
    gs = []
    for k, v in tqdm(cum_tax.items(), total=len(cum_tax), leave=False, ncols=80):
        gs.append(create_tax_graph_w(list(k), v, ref_lengths, ref_stats, scale=scale))

    df = concat_df(gs)

    if len(for_lca) > 0:
        log.info("Creating taxonomic graph with LCA nodes")
        for_lca_df = []
        for k, v in tqdm(for_lca.items(), total=len(for_lca), leave=False, ncols=80):
            for_lca_df.append(create_lca_df(list(k), v))
        df_l = concat_df(for_lca_df)
        df_l["weight"] = 0
        df_l["norm_weight"] = 0
        df = concat_df([df, df_l])
    # group by source and target and sum weight
    df = df.groupby(["source", "target"]).sum().reset_index()

    G = nx.from_pandas_edgelist(df, create_using=nx.DiGraph(), edge_attr=True)
    G.remove_edges_from(nx.selfloop_edges(G))
    Gr = G.reverse()

    log.info("Calculating cumulative weights and paths on the taxonomic graph")
    cumulative_weights_and_paths = {}
    modified_graph = Gr.copy()  # Create a copy of the original graph
    for node in Gr.nodes:
        (
            cumulative_weight,
            cumulative_norm_weight,
            path,
        ) = calculate_cumulative_weight_and_path(modified_graph, node, ref_lengths)
        cumulative_weights_and_paths[node] = {
            "cumulative_weight": cumulative_weight,
            "cumulative_norm_weight": cumulative_norm_weight,
            "path": ";".join(nx.shortest_path(G, target=node)["root"]),
        }

        # Update the modified graph with inferred edge attributes after all weights have been calculated
        for parent in modified_graph.predecessors(node):
            modified_graph[parent][node]["cum_weight"] = (
                modified_graph[parent][node].get("weight", 0) + cumulative_weight
            )
            modified_graph[parent][node]["cum_norm_weight"] = (
                modified_graph[parent][node].get("norm_weight", 0)
                + cumulative_norm_weight
            )

    cumulative_weights_and_paths = dict(
        sorted(cumulative_weights_and_paths.items(), key=lambda item: item[1]["path"])
    )

    if len(for_lca) > 0:
        df1 = concat_df(for_lca_df)
        df1 = df1.groupby(["source", "target"]).sum().reset_index()
        G_lca = nx.from_pandas_edgelist(df1, create_using=nx.DiGraph(), edge_attr=True)
        G_lca.remove_edges_from(nx.selfloop_edges(G_lca))

        df1_d = dict(
            zip(df1[df1["weight"] > 0]["target"], df1[df1["weight"] > 0]["weight"])
        )
        leaves = list(df1_d.keys())

        log.info("Finding most likely reference for the LCA nodes")
        lca2ref = find_most_likely_continuation(
            full_graph=modified_graph.reverse(),
            leaves=leaves,
            index=index_r,
        )

        lca_dfs = []
        for k, v in lca2ref.items():
            tax_path = nx.shortest_path(G, target=k)["root"]

            root_row = pd.DataFrame(
                {"source": ["root"], "target": ["root"], "weight": 0}
            )
            res = list(zip(tax_path, tax_path[1:]))
            # get last element
            res = pd.DataFrame(res, columns=["source", "target"])
            res = pd.concat([root_row, res])
            res = res.drop_duplicates()
            res["weight"] = 0
            res["norm_weight"] = 0
            res.iloc[-1, res.columns.get_loc("weight")] = df1_d[k]
            if v["reference"] is None:
                res.iloc[-1, res.columns.get_loc("norm_weight")] = df1_d[k]
            else:
                res.iloc[-1, res.columns.get_loc("norm_weight")] = round(
                    scale * df1_d[k] / ref_lengths[v["reference"]]
                )
            lca_dfs.append(res)
        log.info("Adding LCA nodes to the taxonomic graph")
        df2 = concat_df(lca_dfs)
        df2 = concat_df([df2, df])
        df2 = df2.groupby(["source", "target"]).sum().reset_index()
        # %%
        G2 = nx.from_pandas_edgelist(df2, create_using=nx.DiGraph(), edge_attr=True)
        G2.remove_edges_from(nx.selfloop_edges(G2))
        Gr2 = G2.reverse()
        # %%
        log.info("Calculating cumulative weights and paths on the taxonomic graph")
        cumulative_weights_and_paths = {}
        modified_graph = Gr2.copy()  # Create a copy of the original graph

        for node in Gr2.nodes:
            (
                cumulative_weight,
                cumulative_norm_weight,
                path,
            ) = calculate_cumulative_weight_and_path(modified_graph, node, ref_lengths)

            cumulative_weights_and_paths[node] = {
                "cumulative_weight": cumulative_weight,
                "cumulative_norm_weight": cumulative_norm_weight,
                "path": ";".join(nx.shortest_path(G, target=node)["root"]),
            }

            # Update the modified graph with inferred edge attributes after all weights have been calculated
            for parent in modified_graph.predecessors(node):
                modified_graph[parent][node]["cum_weight"] = (
                    modified_graph[parent][node].get("weight", 0) + cumulative_weight
                )
                modified_graph[parent][node]["cum_norm_weight"] = (
                    modified_graph[parent][node].get("norm_weight", 0)
                    + cumulative_norm_weight
                )
        cumulative_weights_and_paths = dict(
            sorted(
                cumulative_weights_and_paths.items(), key=lambda item: item[1]["path"]
            )
        )

    log.info("Writing LCA results to file")
    nodes = list(
        set(
            [
                node
                for node, data in cumulative_weights_and_paths.items()
                if data["cumulative_weight"] > 0
            ]
        )
    )
    taxids = txp.taxid_from_name(nodes, taxdb)
    taxids = dict(zip(nodes, taxids))

    out_file = out_files["lca_summary"]
    # determine if the file is gzipped or not
    if out_file.endswith(".gz"):
        with gzip.open(out_file, "wt") as f:
            f.write("taxid\tname\trank\tn_reads\tabundance\ttax_path\n")
            for node, data in tqdm(
                cumulative_weights_and_paths.items(),
                total=len(cumulative_weights_and_paths),
                leave=False,
                ncols=80,
                desc="Writing results",
                unit="nodes",
            ):
                if data["cumulative_weight"] > 0:
                    taxid = taxids[node][0]
                    rank = taxdb.taxid2rank[taxid]
                    f.write(
                        f"{taxid}\t{node}\t{rank}\t{data['cumulative_weight']}\t{data['cumulative_norm_weight']}\t{data['path']}\n"
                    )
    else:
        with open(out_file, "wt") as f:
            f.write("taxid\tname\trank\tn_reads\tabundance\ttax_path\n")
            for node, data in tqdm(
                cumulative_weights_and_paths.items(),
                total=len(cumulative_weights_and_paths),
            ):
                if data["cumulative_weight"] > 0:
                    taxid = taxids[node][0]
                    rank = taxdb.taxid2rank[taxid]
                    f.write(
                        f"{taxid}\t{node}\t{rank}\t{data['cumulative_weight']}\t{data['cumulative_norm_weight']}\t{data['path']}\n"
                    )
