"""
CLI wrapper for the probabilistic taxonomic profiler.

This module provides the CLI interface for computing Bayesian hierarchical
profiles of ancient DNA samples, with dual-channel TAD (ancient/modern)
and credible intervals for P(ancient|taxon).

Uses DCMS (Damage-Calibrated Mixture-of-Sources) approach with belief propagation:
- Beta-Binomial damage model for P(ancient|ref)
- Logistic model for P(present|ref) from coverage features
- Belief propagation on taxonomy tree with latent epoch variables

Input modes:
1. Direct BAM processing (--bam): Computes all statistics de novo
2. Pre-computed stats (--stats-tsv): Uses stats from filter/reassign stage
"""

from pathlib import Path
from typing import Dict, List, Optional
import gzip

from bam_filter import logging as bf_logging

LOG_TAG = "PROB_PROFILE"


def _info(message: str, *args) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args) -> None:
    bf_logging.warn(message, *args)


def _parse_float(value: str, default: float = 0.0) -> float:
    """Parse a float from TSV, returning default for empty/invalid values."""
    if not value or value in ('', 'NA', 'nan', '-', 'None'):
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _parse_int(value: str, default: int = 0) -> int:
    """Parse an int from TSV, returning default for empty/invalid values."""
    if not value or value in ('', 'NA', 'nan', '-', 'None'):
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def _load_ref_stats_from_tsv(stats_path: str) -> Dict[str, dict]:
    """Load reference statistics from filter stage TSV output.

    Loads both basic abundance metrics and detection metrics (WCB, CPC, ORI,
    spatial entropy, gini) needed for computing P(present|taxon).
    """
    import gzip

    ref_stats = {}
    opener = gzip.open if stats_path.endswith('.gz') else open

    with opener(stats_path, 'rt') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue

            row = dict(zip(header, parts))
            ref_name = row.get('reference_name', row.get('reference', row.get('ref_name', '')))
            if not ref_name:
                continue

            ref_stats[ref_name] = {
                # Basic counts
                'n_reads': _parse_float(row.get('total_reads', row.get('n_reads', ''))),
                'n_alns': _parse_float(row.get('total_alignments', row.get('n_alns', row.get('n_reads', '')))),
                'ref_length': _parse_int(row.get('reference_length_bp', row.get('ref_length', row.get('reference_length', '')))),
                'em_posterior_sum': _parse_float(row.get('em_posterior_sum', row.get('total_reads', row.get('n_reads', '')))),

                # Abundance
                'tax_abund_tad': _parse_float(row.get('tax_abund_tad', '')),

                # Coverage metrics
                'breadth': _parse_float(row.get('breadth', '')),
                'mean_coverage': _parse_float(row.get('mean_coverage', row.get('coverage_mean', ''))),

                # Detection metrics for P(present)
                'weighted_contiguity_breadth': _parse_float(row.get('weighted_contiguity_breadth', '')),
                'norm_spatial_entropy': _parse_float(row.get('norm_spatial_entropy', '')),
                'norm_gini': _parse_float(row.get('norm_gini', '')),
                'overlap_redundancy_index': _parse_float(row.get('overlap_redundancy_index', '')),
                'dust_mean': _parse_float(row.get('dust_mean', '')),
                'ani_mean': _parse_float(row.get('read_ani_mean', row.get('ani_mean', ''))),

                # Damage model metrics
                'damage_amplitude': _parse_float(row.get('damage_amplitude', '')),
                'damage_log_bf': _parse_float(row.get('damage_log_bf', '')),
                'damage_baseline': _parse_float(row.get('damage_baseline', '')),

                # Graph-based multimapping metrics (if available)
                'unique_reads': _parse_float(row.get('unique_reads', '')),
                'shared_reads': _parse_float(row.get('shared_reads', '')),
                'connected_neighbors': _parse_float(row.get('connected_neighbors', '')),
                'neighbor_tax_entropy': _parse_float(row.get('neighbor_tax_entropy', '')),
            }

    return ref_stats


def _load_graph_stats(graph_tsv_path: str) -> Dict[str, dict]:
    """Load graph analysis statistics from separate graph TSV file.

    Graph stats include multimapping metrics:
    - unique_reads: Reads mapping only to this reference
    - shared_reads: Reads also mapping to other references
    - connected_neighbors: Number of references sharing reads
    - neighbor_tax_entropy: Taxonomic diversity of neighbors
    """
    import gzip

    graph_stats = {}
    opener = gzip.open if graph_tsv_path.endswith('.gz') else open

    with opener(graph_tsv_path, 'rt') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue

            row = dict(zip(header, parts))
            ref_name = row.get('reference_name', row.get('reference', ''))
            if not ref_name:
                continue

            total_reads = _parse_float(row.get('total_reads', '0'))
            unique_reads = _parse_float(row.get('unique_reads', '0'))
            shared_reads = _parse_float(row.get('shared_reads', '0'))

            exclusivity = unique_reads / total_reads if total_reads > 0 else 0.0

            graph_stats[ref_name] = {
                'total_reads': total_reads,
                'unique_reads': unique_reads,
                'shared_reads': shared_reads,
                'connected_neighbors': _parse_float(row.get('connected_neighbors', '0')),
                'neighbor_tax_entropy': _parse_float(row.get('neighbor_tax_entropy', '0')),
                'read_exclusivity': exclusivity,
                'gamma_ancient': _parse_float(row.get('gamma_ancient', '')),
                'damage_log_bf': _parse_float(row.get('damage_log_bf', '0')),
                'damage_amplitude': _parse_float(row.get('damage_amplitude', '0')),
                'damage_baseline': _parse_float(row.get('damage_baseline', '0')),
            }

    return graph_stats


def _compute_exclusivity_penalty(
    unique_reads: float,
    shared_reads: float,
    neighbors: float,
    eta: float = 0.05,
    alpha: float = 0.7,
    n0: float = 100.0,
    prior_strength: float = 1.0,
) -> tuple:
    """Compute penalty for low read exclusivity using Beta-Binomial model.

    Uses a Bayesian approach to handle uncertainty in exclusivity estimates:
    - Low read count: exclusivity estimate shrunk toward prior (gentle penalty)
    - High read count: exclusivity estimate dominated by data (strong penalty if low)

    Parameters
    ----------
    unique_reads : float
        Number of reads mapping only to this reference
    shared_reads : float
        Number of reads also mapping to other references
    neighbors : float
        Number of connected neighbors (references sharing reads)
    eta : float
        Floor penalty (minimum multiplier). Default 0.05 means max 95% reduction.
    alpha : float
        Steepness of exclusivity effect (0.5-1.0). Default 0.7.
    n0 : float
        Neighbor count for 50% of neighbor effect. Default 100.
    prior_strength : float
        Beta prior strength (higher = more shrinkage). Default 1.0.

    Returns
    -------
    tuple of (penalty, exclusivity_estimate)
        penalty: value in [eta, 1.0] that multiplies p_present
        exclusivity_estimate: Beta posterior mean of exclusivity
    """
    import math

    total = unique_reads + shared_reads
    if total <= 0:
        return (eta, 0.0)

    a0 = prior_strength
    b0 = prior_strength
    excl_post = (a0 + unique_reads) / (a0 + b0 + total)

    g_excl = eta + (1.0 - eta) * math.pow(excl_post, alpha)

    h_neigh = 1.0 / (1.0 + (neighbors / n0))

    penalty = g_excl * h_neigh

    penalty = max(eta * 0.5, min(1.0, penalty))

    return (penalty, excl_post)


def _apply_gmrf_smoothing(
    ref_posteriors: Dict[str, dict],
    graph_stats: Dict[str, dict],
    tau: float = 1.0,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> Dict[str, dict]:
    """Apply GMRF smoothing to reference-level P(ancient) estimates.

    Uses the read-sharing graph to smooth ancientness estimates:
    references that share many reads should have similar P(ancient) values.

    The GMRF model: P(eta) ~ exp(-tau/2 * sum_{r~s} w_rs * (eta_r - eta_s)^2)
    where w_rs = shared_reads_rs / sqrt(n_reads_r * n_reads_s)

    Parameters
    ----------
    ref_posteriors : dict
        Dictionary of {ref_name: {'p_ancient': float, 'log_bf': float, ...}}
    graph_stats : dict
        Graph statistics with neighbor information
    tau : float
        Smoothing precision (higher = more smoothing). Default 1.0.
    max_iter : int
        Maximum coordinate descent iterations
    tol : float
        Convergence tolerance

    Returns
    -------
    dict
        Updated ref_posteriors with smoothed p_ancient values

    Note
    ----
    Currently requires graph adjacency list to be exported separately.
    This is a placeholder until full graph export is implemented.
    The GMRF implementation is in probabilistic_profiler.pyx:apply_gmrf_smoothing().
    """
    import numpy as np

    if not ref_posteriors or not graph_stats:
        return ref_posteriors

    try:
        from bam_filter.probabilistic_profiler import apply_gmrf_smoothing, create_gmrf_hyperparams
    except ImportError:
        _warn("GMRF smoothing not available - module not compiled")
        return ref_posteriors

    ref_names = list(ref_posteriors.keys())
    ref_idx = {name: i for i, name in enumerate(ref_names)}
    n_refs = len(ref_names)

    log_bf = np.zeros(n_refs, dtype=np.float64)
    p_ancient_raw = np.zeros(n_refs, dtype=np.float64)
    n_reads = np.zeros(n_refs, dtype=np.int64)

    for i, name in enumerate(ref_names):
        post = ref_posteriors[name]
        log_bf[i] = post.get('log_bf', 0.0)
        p_ancient_raw[i] = post.get('p_ancient', 0.5)

        if name in graph_stats:
            g = graph_stats[name]
            n_reads[i] = int(g.get('total_reads', 0))
        else:
            n_reads[i] = 1

    neighbors = [[] for _ in range(n_refs)]
    weights = [[] for _ in range(n_refs)]

    _info("GMRF smoothing: adjacency list not available in graph_stats TSV")
    _info("  To enable GMRF smoothing, export graph edges with --export-graph-edges")
    _info("  Using unsmoothed P(ancient) values for now")

    return ref_posteriors


def _load_gamma_from_tsv(stats_path: str) -> Dict[str, float]:
    """Load gamma (P(ancient|ref)) values from filter stage TSV output."""
    import gzip

    gamma_dict = {}
    opener = gzip.open if stats_path.endswith('.gz') else open

    with opener(stats_path, 'rt') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue

            row = dict(zip(header, parts))
            ref_name = row.get('reference_name', row.get('reference', row.get('ref_name', '')))
            if not ref_name:
                continue

            gamma_val = row.get('gamma_ancient', row.get('gamma', ''))
            if gamma_val and gamma_val not in ('', '-1', '-1.0', 'NA', 'nan'):
                try:
                    gamma_dict[ref_name] = float(gamma_val)
                except ValueError:
                    gamma_dict[ref_name] = 0.5
            else:
                gamma_dict[ref_name] = 0.5

    return gamma_dict


def _load_accession_to_taxid(taxonomy_db_path: str, accession_filter: List[str] = None) -> Dict[str, int]:
    """Load accession to taxid mapping from taxonomy database.

    Parameters
    ----------
    taxonomy_db_path : str
        Path to taxonomy database directory
    accession_filter : list of str, optional
        If provided, only load mappings for these accessions (much faster)
    """
    from bam_filter.taxonomy_db import load_accession_map_from_file

    acc_map_path = Path(taxonomy_db_path) / "accession_map.parquet"
    if not acc_map_path.exists():
        raise FileNotFoundError(f"accession_map.parquet not found in {taxonomy_db_path}")

    acc_map = load_accession_map_from_file(str(acc_map_path), accession_filter=accession_filter)
    return acc_map


def _load_lca_per_read(lca_per_read_path: str) -> Dict[str, int]:
    """Load LCA per-read assignments (read_id -> taxid mapping)."""
    opener = gzip.open if lca_per_read_path.endswith('.gz') else open
    read_to_taxid = {}

    with opener(lca_per_read_path, 'rt') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue

            row = dict(zip(header, parts))
            read_id = row.get('read_id', row.get('read_name', ''))
            taxid_str = row.get('taxid', row.get('lca_taxid', ''))

            if read_id and taxid_str:
                try:
                    read_to_taxid[read_id] = int(taxid_str)
                except ValueError:
                    continue

    return read_to_taxid


def _compute_rank_authenticity_pvalues(
    ref_posteriors: Dict[str, dict],
    taxdb,
    min_taxa_per_rank: int = 10,
    verbose: bool = True,
) -> None:
    """
    Compute rank-specific authenticity p-values and update ref_posteriors in-place.

    This replaces the simple sigmoid-based P(present) with a more robust
    rank-specific empirical p-value. The approach:

    1. Group references by their taxon's rank
    2. For each rank with enough taxa, fit a normal to upper 50% of scores
       (assumed to be authentic taxa)
    3. Compute P(present) = CDF(score | authentic_distribution)

    This catches contaminants like E. coli that have positive authenticity
    scores but are in the lower tail for their rank.

    Parameters
    ----------
    ref_posteriors : dict
        Dict of ref_name -> {taxid, p_present, authenticity_score, ...}
        Modified in-place to update p_present values.
    taxdb : TaxonomyDatabase
        Taxonomy database for rank lookups
    min_taxa_per_rank : int
        Minimum taxa per rank to fit a distribution (default: 10)
    verbose : bool
        Print diagnostic info
    """
    import numpy as np
    from scipy import stats as scipy_stats

    # Group by rank and collect scores
    rank_data = {}  # rank -> [(ref_name, score)]

    for ref_name, post in ref_posteriors.items():
        taxid = post.get('taxid')
        score = post.get('authenticity_score')

        if taxid is None or score is None or np.isnan(score):
            continue

        rank = taxdb.get_rank(taxid) or 'unknown'
        if rank not in rank_data:
            rank_data[rank] = []
        rank_data[rank].append((ref_name, score))

    # Fit distribution per rank
    rank_params = {}  # rank -> (mu, std)
    for rank, items in rank_data.items():
        if len(items) < min_taxa_per_rank:
            continue

        scores = np.array([s for _, s in items])
        median_score = np.median(scores)
        upper_half = scores[scores >= median_score]

        mu_auth = np.mean(upper_half)
        std_auth = np.std(upper_half)

        if std_auth < 1e-6:
            std_auth = 0.1

        rank_params[rank] = (mu_auth, std_auth)

        if verbose:
            _info("Authenticity p-value [%s]: n=%d, mu=%.4f, std=%.4f",
                  rank, len(upper_half), mu_auth, std_auth)

    # Fallback: use genus params if available, or best available rank
    fallback_params = rank_params.get('species')
    if fallback_params is None:
        fallback_params = rank_params.get('genus')
    if fallback_params is None and rank_params:
        best_rank = max(rank_params.keys(), key=lambda r: len(rank_data.get(r, [])))
        fallback_params = rank_params[best_rank]
        if verbose:
            _info("Using %s as fallback for ranks without enough samples", best_rank)

    if not fallback_params:
        if verbose:
            _info("Not enough taxa to compute rank-specific p-values, keeping sigmoid P(present)")
        return

    # Update p_present for all refs using rank-specific or fallback params
    n_updated = 0
    for ref_name, post in ref_posteriors.items():
        taxid = post.get('taxid')
        score = post.get('authenticity_score')

        if taxid is None or score is None or np.isnan(score):
            continue

        rank = taxdb.get_rank(taxid) or 'unknown'
        mu, std = rank_params.get(rank, fallback_params)

        # P(present) = CDF(score | authentic_distribution)
        # High score -> high CDF -> high P(present) (authentic)
        # Low score -> low CDF -> low P(present) (contamination)
        p_present_new = scipy_stats.norm.cdf(score, loc=mu, scale=std)

        # Store both old and new values for diagnostics
        post['p_present_sigmoid'] = post['p_present']
        post['p_present'] = p_present_new
        post['authenticity_pvalue'] = p_present_new

        # Recompute TAD splits with new p_present
        tad = post.get('tad', 0.0)
        p_ancient = post.get('p_ancient', 0.5)
        post['tad_ancient'] = tad * p_ancient * p_present_new
        post['tad_modern'] = tad * (1.0 - p_ancient) * p_present_new

        n_updated += 1

    if verbose:
        _info("Updated P(present) using rank-specific p-values for %d references", n_updated)


def _compute_taxon_authenticity_pvalues(
    profile: Dict[int, dict],
    taxdb,
    min_taxa_per_rank: int = 10,
    verbose: bool = True,
) -> None:
    """
    Compute rank-specific authenticity p-values at the TAXON level.

    This is computed AFTER belief propagation, using the TAD-weighted
    average authenticity score aggregated from all descendant references.

    The p-value indicates how likely this taxon's authenticity score is
    compared to other taxa at the same rank (species, genus, etc.).

    Low p-value = lower authenticity than typical for this rank = likely contamination
    High p-value = higher authenticity than typical = likely authentic

    Parameters
    ----------
    profile : dict
        Dict of taxid -> {authenticity_score, ...}
        Modified in-place to add authenticity_pvalue.
    taxdb : TaxonomyDatabase
        Taxonomy database for rank lookups
    min_taxa_per_rank : int
        Minimum taxa per rank to fit a distribution (default: 10)
    verbose : bool
        Print diagnostic info
    """
    import numpy as np
    from scipy import stats as scipy_stats

    # Group by rank and collect scores
    rank_data = {}  # rank -> [(taxid, score)]

    for taxid, entry in profile.items():
        score = entry.get('authenticity_score', 0.0)

        if score == 0.0 or np.isnan(score):
            continue

        rank = taxdb.get_rank(taxid) or 'unknown'
        if rank not in rank_data:
            rank_data[rank] = []
        rank_data[rank].append((taxid, score))

    # Fit distribution per rank
    rank_params = {}  # rank -> (mu, std)
    for rank, items in rank_data.items():
        if len(items) < min_taxa_per_rank:
            continue

        scores = np.array([s for _, s in items])
        median_score = np.median(scores)
        upper_half = scores[scores >= median_score]

        mu_auth = np.mean(upper_half)
        std_auth = np.std(upper_half)

        if std_auth < 1e-6:
            std_auth = 0.1

        rank_params[rank] = (mu_auth, std_auth)

        if verbose:
            _info("Taxon authenticity [%s]: n=%d, mu=%.4f, std=%.4f",
                  rank, len(upper_half), mu_auth, std_auth)

    # Fallback: use species params if available, or best available rank
    fallback_params = rank_params.get('species')
    if fallback_params is None:
        fallback_params = rank_params.get('genus')
    if fallback_params is None and rank_params:
        best_rank = max(rank_params.keys(), key=lambda r: len(rank_data.get(r, [])))
        fallback_params = rank_params[best_rank]
        if verbose:
            _info("Using %s as fallback for ranks without enough samples", best_rank)

    if not fallback_params:
        if verbose:
            _info("Not enough taxa to compute authenticity p-values")
        for entry in profile.values():
            entry['authenticity_pvalue'] = 0.5
        return

    # Compute p-value for all taxa
    for taxid, entry in profile.items():
        score = entry.get('authenticity_score', 0.0)

        if score == 0.0 or np.isnan(score):
            entry['authenticity_pvalue'] = 0.5
            continue

        rank = taxdb.get_rank(taxid) or 'unknown'
        mu, std = rank_params.get(rank, fallback_params)

        # P(present) = CDF(score | authentic_distribution)
        # High score -> high CDF -> high p-value (authentic)
        # Low score -> low CDF -> low p-value (contamination)
        p_value = scipy_stats.norm.cdf(score, loc=mu, scale=std)
        entry['authenticity_pvalue'] = p_value


def _compute_dcms_stats_from_bam(
    bam_path: str,
    taxonomy_db_path: str,
    num_threads: int = 1,
    verbose: bool = True,
) -> str:
    """
    Compute per-reference statistics from BAM for DCMS profiling.

    Uses compute_bam_stats to get full coverage metrics including
    entropy, gini, authenticity_score, and damage parameters.
    """
    import tempfile
    from pathlib import Path
    from bam_filter.stats import compute_bam_stats

    _info("Computing per-reference statistics from BAM...")
    _info("BAM: %s", bam_path)

    # Create temp file for stats output
    prefix = Path(bam_path).stem
    stats_fd, stats_path = tempfile.mkstemp(prefix=f"{prefix}_dcms_", suffix=".tsv.gz")
    import os
    os.close(stats_fd)

    # Compute full reference statistics
    ret = compute_bam_stats(
        bam_file=bam_path,
        output=stats_path,
        num_threads=num_threads,
        verbose=verbose,
        min_read_count=1,
        min_read_length=30,
        max_read_length=10000,
        min_read_ani=0.0,
        damage_correction=True,
    )

    if ret != 0:
        raise RuntimeError(f"compute_bam_stats failed with code {ret}")

    _info("Per-reference statistics written to %s", stats_path)
    return stats_path


def _compute_stats_from_bam_capsule(
    bam_path: str,
    num_threads: int = 1,
    verbose: bool = True,
) -> dict:
    """
    Compute per-reference statistics directly from BAM file (capsule mode).

    Returns a PyCapsule containing RefStats* pointer for direct Cython access,
    avoiding Python dict conversion overhead.

    Returns
    -------
    dict with keys:
        'capsule': PyCapsule containing RefStats*
        'ref_names': list of reference names
        'n_refs': number of references
    """
    from bam_filter.stats import compute_bam_stats

    _info("Computing coverage and damage statistics from BAM (capsule mode)...")
    result = compute_bam_stats(
        bam_file=bam_path,
        num_threads=num_threads,
        verbose=verbose,
        min_read_count=1,
        min_read_length=30,
        max_read_length=10000,
        min_read_ani=0.0,
        damage_correction=True,
        return_stats='capsule',
    )

    if not isinstance(result, dict) or 'capsule' not in result:
        raise RuntimeError(f"compute_bam_stats capsule mode failed: {result}")

    _info("Computed stats for %d references (capsule mode)", result['n_refs'])
    return result


def _run_dcms_profile_from_capsule(
    stats_result: dict,
    taxonomy_db_path: str,
    output_path: str,
    hyperparams: Optional[dict] = None,
    num_threads: int = 4,
    verbose: bool = True,
) -> dict:
    """
    Run DCMS profiler directly from stats capsule (full Cython pipeline).

    This function uses a single Cython call that handles:
    - Posterior computation from RefStats
    - Accession-to-taxid mapping
    - Taxonomy tree construction
    - Belief propagation

    No Python intermediate data structures in the hot path.

    Parameters
    ----------
    stats_result : dict
        Result from _compute_stats_from_bam_capsule with keys:
        'capsule', 'ref_names', 'n_refs'
    taxonomy_db_path : str
        Path to taxonomy database directory
    output_path : str
        Output path for profile TSV (unused here, kept for API compatibility)
    hyperparams : dict, optional
        Override default hyperparameters
    num_threads : int
        Number of threads for parallel computation
    verbose : bool
        Print progress messages

    Returns
    -------
    dict
        Taxid-keyed profile with p_ancient, p_present, TAD values
    """
    from bam_filter.probabilistic_profiler import run_full_profile_from_capsule
    from bam_filter.taxonomy_db import TaxonomyDatabase, load_accession_map_from_file

    _info("Starting DCMS probabilistic profiler (full Cython pipeline)")

    capsule = stats_result['capsule']
    ref_names = stats_result['ref_names']
    n_refs = stats_result['n_refs']

    # Load taxonomy database (Cython object)
    _info("Loading taxonomy database...")
    taxdb = TaxonomyDatabase.from_parquet(taxonomy_db_path)

    # Load accession map (Cython object, filtered to only references in the stats)
    _info("Loading accession map (filtered to %d references)...", n_refs)
    acc_map_path = Path(taxonomy_db_path) / "accession_map.parquet"
    if not acc_map_path.exists():
        raise FileNotFoundError(f"accession_map.parquet not found in {taxonomy_db_path}")

    acc_map = load_accession_map_from_file(str(acc_map_path), accession_filter=ref_names)
    _info("Accession map loaded")

    # Run full Cython pipeline
    _info("Running full Cython pipeline...")
    profile = run_full_profile_from_capsule(
        stats_capsule=capsule,
        ref_names=ref_names,
        taxonomy_db=taxdb,
        acc_map=acc_map,
        num_threads=num_threads,
        hyperparams=hyperparams,
        verbose=verbose,
    )

    _info("Profile complete: %d taxa", len(profile))

    # Compute rank-specific authenticity p-values at taxon level
    _info("Computing rank-specific authenticity p-values at taxon level...")
    _compute_taxon_authenticity_pvalues(profile, taxdb, min_taxa_per_rank=10, verbose=verbose)

    return profile


def _format_dcms_profile_tsv(
    profile: dict,
    taxonomy,
    output_path: str,
    min_reads: float = 0.0,
    min_p_present: float = 0.0,
    verbose: bool = True,
) -> int:
    """Format DCMS profile results to TSV."""
    import gzip

    opener = gzip.open if output_path.endswith('.gz') else open

    header = [
        'taxid', 'parent_taxid', 'rank', 'name',
        'n_refs', 'tad_total', 'tad_ancient', 'tad_modern',
        'p_ancient', 'p_present',
        'authenticity_score', 'authenticity_pvalue',
        'read_exclusivity',
    ]

    n_written = 0
    with opener(output_path, 'wt') as f:
        f.write('\t'.join(header) + '\n')

        for taxid in sorted(profile.keys()):
            entry = profile[taxid]

            # Apply filters
            if entry['tad_total'] < min_reads:
                continue
            if entry['p_present'] < min_p_present:
                continue

            parent = taxonomy.get_parent(taxid)
            rank = taxonomy.get_rank(taxid)
            name = taxonomy.get_name(taxid)

            row = [
                str(taxid),
                str(parent) if parent else '0',
                rank if rank else 'no_rank',
                name if name else f'taxid_{taxid}',
                str(entry.get('n_refs', 0)),
                f"{entry['tad_total']:.4f}",
                f"{entry['tad_ancient']:.4f}",
                f"{entry['tad_modern']:.4f}",
                f"{entry['p_ancient']:.6f}",
                f"{entry['p_present']:.6f}",
                f"{entry.get('authenticity_score', 0.0):.4f}",
                f"{entry.get('authenticity_pvalue', 0.5):.6f}",
                f"{entry.get('read_exclusivity', 1.0):.4f}",
            ]
            f.write('\t'.join(row) + '\n')
            n_written += 1

    return n_written


def _compute_stats_from_bam(
    bam_path: str,
    lca_per_read_path: str,
    taxonomy_db_path: str,
    num_threads: int = 1,
    verbose: bool = False,
) -> Dict[int, dict]:
    """
    Compute per-taxon statistics directly from BAM and LCA assignments.

    Returns a dictionary mapping taxid -> stats dict with:
    - n_reads: number of reads assigned to this taxon
    - breadth, coverage, wcb, cpc, ori, entropy, gini metrics
    - damage_amplitude, damage_baseline from damage model fit
    - gamma_ancient estimated from damage pattern
    """
    from bam_filter.processor_lca_stats import (
        process_lca_stats_from_bam,
        configure_lca_stats_thresholds,
    )
    from bam_filter.taxonomy_db import TaxonomyDatabase

    _info("Computing statistics from BAM file...")
    _info("BAM: %s", bam_path)
    _info("LCA per-read: %s", lca_per_read_path)

    configure_lca_stats_thresholds(
        min_read_ani=0.0,
        min_read_length=0,
        max_read_length=(2**31) - 1,
        scale=1_000_000,
        trim_ends=0,
        trim_min=10,
        trim_max=90,
    )

    if verbose:
        _info("Loading taxonomy database...")
    taxdb = TaxonomyDatabase.from_parquet(taxonomy_db_path)

    taxon_stats = process_lca_stats_from_bam(
        bam_path.encode('utf-8'),
        lca_per_read_path.encode('utf-8'),
        taxdb,
        taxonomy_db_path.encode('utf-8'),
        num_threads,
        verbose,
    )

    _info("Computed stats for %d taxids", len(taxon_stats))
    return taxon_stats


def do_prob_profile(args):
    """
    Execute the profiler subcommand.

    Uses DCMS (Damage-Calibrated Mixture-of-Sources) approach with belief propagation.
    All statistics including gamma_ancient (P(ancient)) are computed fresh from the BAM
    using Cython-only computation for efficiency.

    Computes a probabilistic taxonomic profile with:
    - Dual-channel TAD (ancient/modern) weighted by EM posteriors
    - Beta-distributed P(ancient|taxon) with hierarchical shrinkage
    - Credible intervals for P(ancient|taxon)
    """
    from bam_filter.taxonomy_db import TaxonomyDatabase

    bam_path = getattr(args, "bam", None)
    if not bam_path:
        raise ValueError("--bam is required for the profiler command")

    taxonomy_db_path = Path(args.taxonomy_db)
    if not taxonomy_db_path.exists() or not taxonomy_db_path.is_dir():
        raise FileNotFoundError(f"Taxonomy database directory not found: {taxonomy_db_path}")

    num_threads = getattr(args, "threads", 1)
    verbose = not getattr(args, "quiet", False)

    # Compute stats from BAM using capsule mode (Cython-only)
    stats_result = _compute_stats_from_bam_capsule(
        bam_path=bam_path,
        num_threads=num_threads,
        verbose=verbose,
    )

    # Determine output path
    prefix_arg = getattr(args, "prefix", None)
    if prefix_arg:
        prefix = prefix_arg
    else:
        prefix = Path(bam_path).with_suffix("").name

    output_path = getattr(args, "output", None)
    if not output_path:
        output_path = f"{prefix}.profile.tsv"

    min_reads = float(getattr(args, "min_reads", 0.0))
    min_p_present = float(getattr(args, "min_p_present", 0.0))

    # Run DCMS profiler using Cython capsule (no intermediate TSV files)
    profile = _run_dcms_profile_from_capsule(
        stats_result=stats_result,
        taxonomy_db_path=str(taxonomy_db_path),
        output_path=output_path,
        num_threads=num_threads,
        verbose=verbose,
    )

    if not profile:
        _warn("No taxids in profile. Check input data.")
        return

    # Load taxonomy for output formatting
    taxdb = TaxonomyDatabase.from_parquet(str(taxonomy_db_path))

    # Write output
    _info("Writing profile to %s", output_path)
    n_written = _format_dcms_profile_tsv(
        profile=profile,
        taxonomy=taxdb,
        output_path=output_path,
        min_reads=min_reads,
        min_p_present=min_p_present,
        verbose=verbose,
    )

    _info("Wrote %d taxids to %s", n_written, output_path)

    # Summary
    if profile:
        n_ancient = sum(1 for e in profile.values() if e.get('p_ancient', 0) > 0.7)
        n_modern = sum(1 for e in profile.values() if e.get('p_ancient', 1) < 0.3)
        n_uncertain = len(profile) - n_ancient - n_modern
        n_confident = sum(1 for e in profile.values() if e.get('p_present', 0) > 0.7)

        print("")
        print("=" * 60)
        print("PROBABILISTIC PROFILE SUMMARY")
        print("=" * 60)
        print(f"Total taxids:                    {len(profile)}")
        print("")
        print("P(ancient) - Ancient DNA probability:")
        print(f"  Ancient (P > 0.7):             {n_ancient}")
        print(f"  Modern (P < 0.3):              {n_modern}")
        print(f"  Uncertain (0.3-0.7):           {n_uncertain}")
        print("")
        print("P(present) - Detection confidence:")
        print(f"  High confidence (P > 0.7):     {n_confident}")
        print("")
        print(f"Output:                          {output_path}")
        print("=" * 60)

    bf_logging.summary("Profile written to %s", Path(output_path).resolve())
