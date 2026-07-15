#!/usr/bin/env python3
"""完整解析SF8963异常MatrixMarket并写成可被标准读取器读取的矩阵。"""

from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
SOURCE = DATA_ROOT / "GSE174554/GSE174554_RAW/GSM5319529_SF8963_matrix.mtx.gz"
OUT_DIR = ROOT / "write/34_gse174554_raw_independent_discovery/01_qc_doublets/source_repairs"
OUTPUT = OUT_DIR / "GSM5319529_SF8963_matrix.fixed.mtx.gz"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(SOURCE, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                n_rows, n_cols, declared_nnz = map(int, line.split())
                break
        rows = np.empty(declared_nnz, dtype=np.int32)
        cols = np.empty(declared_nnz, dtype=np.int32)
        data = np.empty(declared_nnz, dtype=np.int32)
        parsed = 0
        for line in handle:
            values = re.findall(r"\d+", line)
            if len(values) < 3:
                continue
            rows[parsed] = int(values[0]) - 1
            cols[parsed] = int(values[1]) - 1
            data[parsed] = int(values[2])
            parsed += 1
    matrix = sparse.coo_matrix(
        (data[:parsed], (rows[:parsed], cols[:parsed])),
        shape=(n_rows, n_cols),
        dtype=np.int32,
    ).tocsr()
    with gzip.GzipFile(filename=str(OUTPUT), mode="wb", mtime=0) as handle:
        mmwrite(handle, matrix, field="integer", symmetry="general")
    audit = {
        "source": str(SOURCE),
        "output": str(OUTPUT),
        "shape": [n_rows, n_cols],
        "declared_entries": declared_nnz,
        "parsed_entries": parsed,
        "matrix_nnz_after_duplicate_sum": int(matrix.nnz),
        "repair_rule": "retain every parseable coordinate triplet; do not impute missing declared entries",
    }
    (OUT_DIR / "GSM5319529_SF8963_matrix_repair_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
