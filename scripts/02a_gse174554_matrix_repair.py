#!/usr/bin/env python3
"""Recover readable coordinate triplets from the malformed SF8963 MatrixMarket file.
Input: MILLER_IM_DATA_ROOT/GSE174554/GSE174554_RAW/GSM5319529_SF8963_matrix.mtx.gz.
Outputs under write/34_gse174554_raw_independent_discovery/01_qc_doublets/source_repairs:
    GSM5319529_SF8963_matrix.fixed.mtx.gz and GSM5319529_SF8963_matrix_repair_audit.json.
Run from the repository root with python scripts/02a_gse174554_matrix_repair.py.
The existing Python workflow uses SciPy sparse matrix parsing and mmwrite with the original repair rules.
"""

from __future__ import annotations

import gzip
import importlib
import json
import os
from pathlib import Path

from scipy.io import mmwrite


data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # Read the declared dimensions from gzip and skip blank lines or records with fewer than three numbers.
    root = Path(__file__).resolve().parents[1]
    data_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    source_path = data_root / "GSE174554/GSE174554_RAW/GSM5319529_SF8963_matrix.mtx.gz"
    matrix, declared_nnz, parsed_entries = data_functions.read_matrix_manual(source_path)
    n_rows, n_columns = matrix.shape

    # Write the repaired matrix to the existing directory with the original deterministic mtime=0 compression.
    output_dir = root / "write/34_gse174554_raw_independent_discovery/01_qc_doublets/source_repairs"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "GSM5319529_SF8963_matrix.fixed.mtx.gz"
    with gzip.GzipFile(filename=str(output_path), mode="wb", mtime=0) as handle:
        mmwrite(handle, matrix, field="integer", symmetry="general")

    # Keep the published repair record fields and notes, and calculate counts from the parsed matrix.
    audit = {
        "source": str(source_path),
        "output": str(output_path),
        "shape": [n_rows, n_columns],
        "declared_entries": declared_nnz,
        "parsed_entries": parsed_entries,
        "matrix_nnz_after_duplicate_sum": int(matrix.nnz),
        "repair_rule": "retain every parseable coordinate triplet; do not impute missing declared entries",
    }
    audit_path = output_dir / "GSM5319529_SF8963_matrix_repair_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    # Check that the two required repair outputs exist and are nonempty.
    for path in (output_path, audit_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")

    print(f"Parsed entries: {parsed_entries}; declared entries: {declared_nnz}; matrix nnz: {matrix.nnz}")
