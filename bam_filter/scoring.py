import numpy as np
import logging
from typing import Tuple, List, Dict, Optional, Union, Any

log = logging.getLogger("my_logger")


def calculate_ani_from_tags(query_length: int, nm_value: int) -> float:
    """
    Calculate Average Nucleotide Identity (ANI) from alignment data.

    Args:
        query_length: Length of the query sequence
        nm_value: Number of mismatches (from NM tag)

    Returns:
        float: ANI value between 0 and 1
    """
    if query_length == 0:
        return 0.0

    # ANI = 1 - (mismatches / query_length)
    return 1.0 - (nm_value / query_length)


def calculate_alignment_score(
    num_matches: int,
    num_mismatches: int,
    num_gaps: int,
    gap_extensions: int,
    match_reward: int = 1,
    mismatch_penalty: int = -2,
    gap_open_penalty: int = 5,
    gap_extension_penalty: int = 2,
) -> int:
    """
    Calculate alignment score using the provided penalties and rewards.

    Args:
        num_matches: Number of matching bases
        num_mismatches: Number of mismatching bases
        num_gaps: Number of gap openings
        gap_extensions: Number of gap extensions
        match_reward: Score for a matching base (default: 1)
        mismatch_penalty: Penalty for a mismatching base (default: -2)
        gap_open_penalty: Penalty for opening a gap (default: 5)
        gap_extension_penalty: Penalty for extending a gap (default: 2)

    Returns:
        int: The calculated alignment score
    """
    return (
        (num_matches * match_reward)
        - (num_mismatches * mismatch_penalty)
        - (num_gaps * gap_open_penalty)
        - (gap_extensions * gap_extension_penalty)
    )


def calculate_shifted_scores(
    raw_scores: np.ndarray, aln_lengths: np.ndarray
) -> np.ndarray:
    """
    Calculate shifted and normalized scores.

    Args:
        raw_scores: Array of raw alignment scores
        aln_lengths: Array of alignment lengths

    Returns:
        np.ndarray: Array of shifted and normalized scores
    """
    if len(raw_scores) == 0 or len(aln_lengths) == 0:
        return np.array([])

    # Shift scores by subtracting the minimum score and adding 1
    # Then normalize by alignment length
    return (raw_scores - np.min(raw_scores) + 1) / aln_lengths


def process_alignment_batch_scores(
    alignment_info: List[Tuple[str, str, int, int, int]],
) -> List[Tuple[str, str, float, int]]:
    """
    Process a batch of alignments to calculate shifted scores.

    Args:
        alignment_info: List of tuples containing (query_name, reference_name, reference_length, score, query_alignment_length)

    Returns:
        List of tuples containing (query_name, reference_name, shifted_score, reference_length)
    """
    if not alignment_info:
        return []

    # Convert to numpy arrays for faster operations
    raw_scores = np.array([info[3] for info in alignment_info])
    aln_lengths = np.array([info[4] for info in alignment_info])

    # Calculate shifted and normalized scores
    shifted_scores = calculate_shifted_scores(raw_scores, aln_lengths)

    # Create final alignment data
    return [
        (info[0], info[1], score, info[2])
        for info, score in zip(alignment_info, shifted_scores)
    ]


def extract_alignment_features_from_tags(
    tag_map: Dict[str, Any], query_length: int
) -> Tuple[int, int, int, int, float]:
    """
    Extract alignment features from SAM/BAM tags.

    Args:
        tag_map: Dictionary of tag values
        query_length: Length of the query sequence

    Returns:
        Tuple containing:
            - num_mismatches: Number of mismatches
            - num_matches: Number of matches
            - num_gaps: Number of gap opens
            - gap_extensions: Number of gap extensions
            - ani: Average Nucleotide Identity
    """
    # Extract required tags with fallbacks
    num_mismatches = tag_map.get("NM", 0)
    num_gaps = tag_map.get("XO", 0)
    gap_extensions = tag_map.get("XG", 0)

    # Calculate derived values
    num_matches = query_length - num_mismatches
    ani = calculate_ani_from_tags(query_length, num_mismatches)

    return num_mismatches, num_matches, num_gaps, gap_extensions, ani
