#!/usr/bin/env python3
"""目的：完整解析 SF8963 的异常 MatrixMarket，保留每个能读取的坐标三元组。
输入：MILLER_IM_DATA_ROOT/GSE174554/GSE174554_RAW/GSM5319529_SF8963_matrix.mtx.gz。
输出：write/34_gse174554_raw_independent_discovery/01_qc_doublets/source_repairs/ 下的
    GSM5319529_SF8963_matrix.fixed.mtx.gz 和 GSM5319529_SF8963_matrix_repair_audit.json。
运行：在仓库根目录执行 python scripts/02a_GSE174554补充矩阵.py。
保留 Python：已有 scipy 稀疏矩阵解析及 mmwrite 路径是公开复现合同，本轮不更换语言或修复规则。
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

    # 从原始 gzip 文本读取声明的维度；空行及不足三个数字的记录沿用原规则跳过。
    root = Path(__file__).resolve().parents[1]
    data_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    source_path = data_root / "GSE174554/GSE174554_RAW/GSM5319529_SF8963_matrix.mtx.gz"
    matrix, declared_nnz, parsed_entries = data_functions.read_matrix_manual(source_path)
    n_rows, n_columns = matrix.shape

    # 修复矩阵写到既有目录；mtime=0 保留既有压缩设置，原始矩阵不覆盖。
    output_dir = root / "write/34_gse174554_raw_independent_discovery/01_qc_doublets/source_repairs"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "GSM5319529_SF8963_matrix.fixed.mtx.gz"
    with gzip.GzipFile(filename=str(output_path), mode="wb", mtime=0) as handle:
        mmwrite(handle, matrix, field="integer", symmetry="general")

    # 审计 JSON 为现有修复来源记录，字段及说明原文保留，数字取自本次解析对象。
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

    # 保存后核对两个必要产物非空；不重新解释或填补原始矩阵。
    for path in (output_path, audit_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")

    print(f"Parsed entries: {parsed_entries}; declared entries: {declared_nnz}; matrix nnz: {matrix.nnz}")
