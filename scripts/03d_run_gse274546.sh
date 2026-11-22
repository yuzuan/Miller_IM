#!/usr/bin/env bash
# Purpose: Run GSE274546 raw-matrix quality control, then independently annotate the cells.
# Input: MILLER_IM_DATA_ROOT, or the repository data/ directory by default.
# Output: write/36_gse274546_raw_independent_reannotation/ and the downstream figure directories.
# Run: bash scripts/03d_run_gse274546.sh
# Activate the Python environment before running. The script uses python3 from PATH.

set -euo pipefail

# The R quality-control entry point takes the repository root.
root="$(cd "$(dirname "$0")/.." && pwd)"
Rscript "$root/scripts/03b_gse274546_quality_control.R" "$root"

# The reannotation reads the saved quality-control output.
python3 "$root/scripts/03c_gse274546_independent_reannotation.py"
