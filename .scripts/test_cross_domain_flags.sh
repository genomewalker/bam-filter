#!/bin/bash
# Test script for the new cross-domain removal flags
#
# New flags (replacing old ones):
#   --remove-cross-domain-alignments  (replaces --remove-cross-domain-edges)
#   --remove-cross-domain-references  (replaces --taxonomy-strict-filter when used alone)
#   --remove-cross-domain-all         (shorthand for both - NEW combined mode)
#   --detect-misannotations           (replaces --flag-misannotations)
#
# Usage examples:

# Example 1: Remove only cross-domain alignments (preserves references)
# filterBAM reassign \
#     --bam input.bam \
#     -o output.bam \
#     --taxonomy-db /path/to/taxonomy \
#     --taxonomy-filter \
#     --clustering \
#     --remove-cross-domain-alignments \
#     --detect-misannotations

# Example 2: Remove cross-domain references (node-level removal)
# filterBAM reassign \
#     --bam input.bam \
#     -o output.bam \
#     --taxonomy-db /path/to/taxonomy \
#     --taxonomy-filter \
#     --clustering \
#     --remove-cross-domain-references \
#     --taxonomy-min-connections 5

# Example 3: Combined mode - most aggressive decontamination
# filterBAM reassign \
#     --bam input.bam \
#     -o output.bam \
#     --taxonomy-db /path/to/taxonomy \
#     --taxonomy-filter \
#     --clustering \
#     --remove-cross-domain-all \
#     --detect-misannotations

# Verify flags are available
echo "Testing new cross-domain removal flags..."
filterBAM reassign --help 2>&1 | grep -E "remove-cross-domain|detect-misannotation"

echo ""
echo "Flag verification complete."
