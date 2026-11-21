#!/usr/bin/env python3
"""目的：从 GEO SOFT 元数据构建 GSE278456 髓系及 GSE276841 空间样本下载清单。
输入：--single-cell-soft、--spatial-soft；读取 SAMPLE、title、characteristics 和 supplementary_file。
输出：--output-dir 下的 GSE278456_sample_manifest.csv、GSE276841_sample_manifest.csv、
    eligible_download_manifest.csv 和 eligible_download_aria2.txt。
运行：python scripts/01b_GSE278456样本清单.py --single-cell-soft FILE --spatial-soft FILE --output-dir DIR。
保留 Python：已有 pandas、gzip 和 URL 解析构成公开复现路径，本轮保持筛选、字段及文件名。
"""

import argparse
import importlib
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # 命令行参数沿用既有接口；两个输入均为已下载 SOFT，本脚本只生成下载清单。
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-cell-soft", type=Path, required=True)
    parser.add_argument("--spatial-soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 输入 SOFT 路径，输出按原始顺序排列的样本字典列表；缺文件时保留读取错误。

    single_cell_samples = data_functions.parse_soft(args.single_cell_soft)
    single_cell_rows = []
    download_rows = []
    for sample in single_cell_samples:
        classification = str(sample.get("2021_classification", ""))
        is_idhwt_gbm = classification.lower() == "glioblastoma, idh-wt"
        myeloid_urls = [url for url in sample["supplementary_urls"] if "_Myl_filtered_counts.h5" in url]
        sample_row = {key: value for key, value in sample.items() if key != "supplementary_urls"}
        sample_row["is_idhwt_gbm"] = is_idhwt_gbm
        sample_row["myeloid_h5_url"] = myeloid_urls[0] if myeloid_urls else ""
        single_cell_rows.append(sample_row)

        # 仅严格 IDH-wt GBM 且有既定髓系 H5 的样本进入下载清单；保留每样本第一个 URL。
        if not (is_idhwt_gbm and myeloid_urls):
            continue

        download_rows.append(
            {
                "dataset": "GSE278456",
                "sample": sample["sample_title"],
                "geo_accession": sample["geo_accession"],
                "modality": "scRNA_myeloid",
                "url": myeloid_urls[0],
                "local_subdir": "source_metadata/GSE278456_myeloid_h5",
                "local_name": Path(urlparse(myeloid_urls[0]).path).name,
            }
        )

    # 空间 SOFT 不新增筛选，保留全部样本及各样本附件顺序。
    spatial_samples = data_functions.parse_soft(args.spatial_soft)
    spatial_rows = []
    for sample in spatial_samples:
        sample_row = {key: value for key, value in sample.items() if key != "supplementary_urls"}
        sample_row["n_supplementary_files"] = len(sample["supplementary_urls"])
        spatial_rows.append(sample_row)
        for url in sample["supplementary_urls"]:
            download_rows.append(
                {
                    "dataset": "GSE276841",
                    "sample": sample["sample_title"],
                    "geo_accession": sample["geo_accession"],
                    "modality": "Visium",
                    "url": url,
                    "local_subdir": "source_metadata/GSE276841_spatial",
                    "local_name": Path(urlparse(url).path).name,
                }
            )

    # 样本表一行一个 GEO 样本；下载总表一行一个附件，保持公开列名及输出顺序。
    single_cell_table = pd.DataFrame(single_cell_rows)
    spatial_table = pd.DataFrame(spatial_rows)
    download_table_total = pd.DataFrame(download_rows)
    single_cell_path = args.output_dir / "GSE278456_sample_manifest.csv"
    single_cell_table.to_csv(single_cell_path, index=False)
    spatial_path = args.output_dir / "GSE276841_sample_manifest.csv"
    spatial_table.to_csv(spatial_path, index=False)
    download_path = args.output_dir / "eligible_download_manifest.csv"
    download_table_total.to_csv(download_path, index=False)

    # aria2 文本每附件三行，与下载总表共用同一行顺序和目标相对目录。
    aria_lines = []
    for row in download_table_total.itertuples(index=False):
        aria_lines.extend([row.url, f"  dir={args.output_dir / row.local_subdir}", f"  out={row.local_name}"])

    aria_path = args.output_dir / "eligible_download_aria2.txt"
    aria_path.write_text("\n".join(aria_lines) + "\n", encoding="utf-8")

    # 回读清单核对行列数，检查 aria2 内容非空；任何损坏直接停止。
    for path, expected_table in (
        (single_cell_path, single_cell_table),
        (spatial_path, spatial_table),
        (download_path, download_table_total),
    ):
        saved_table = pd.read_csv(path)
        if saved_table.shape != expected_table.shape or saved_table.columns.tolist() != expected_table.columns.tolist():
            raise RuntimeError(f"Manifest dimensions or columns changed during save: {path}")

    if not aria_path.is_file() or aria_path.stat().st_size == 0:
        raise RuntimeError("Missing or empty aria2 manifest")

    print(f"GSE278456 samples: {len(single_cell_table)}")
    print(f"Strict IDH-wildtype GBM with myeloid H5: {int(single_cell_table.is_idhwt_gbm.sum())}")
    print(f"GSE276841 spatial samples: {len(spatial_table)}")
    print(f"Files selected for download: {len(download_table_total)}")
