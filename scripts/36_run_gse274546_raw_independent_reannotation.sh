#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPO_ROOT="$(cd "$ROOT/.." && pwd)"
PYTHON="$REPO_ROOT/.venv/bin/python"

Rscript "$ROOT/scripts/36_gse274546_raw_qc_scdblfinder.R" "$ROOT"
"$PYTHON" "$ROOT/scripts/36_gse274546_raw_independent_reannotation.py"
"$PYTHON" "$ROOT/scripts/36_gse274546_paired_patient_harmony_clustering.py"
Rscript "$ROOT/scripts/36_gse274546_paired_cluster_da.R" "$ROOT"
Rscript "$ROOT/scripts/36_gse274546_author_label_concordance.R" "$ROOT"
"$PYTHON" "$ROOT/scripts/36_gse274546_annotation_figures.py"
