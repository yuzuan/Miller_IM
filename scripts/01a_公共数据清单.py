#!/usr/bin/env python3
"""目的：构建GSE174554与GSE274546的样本、文库和合并计数矩阵清单。
输入：--project-root下的GEO样本注释、原始文库、肿瘤/正常细胞注释；--geo-cache-dir缓存文本；
    --gse274546-cache-dir标准MatrixMarket。输入字段与CLI接口沿用公开版本。
输出：--output-dir下gse274546_source_samples、gse274546_samples、gse174554_libraries、
    gse174554_samples、gse274546_standard_libraries、gbm_standard_libraries、dataset_summary的CSV及TSV。
运行：python scripts/01a_公共数据清单.py --project-root DIR --geo-cache-dir DIR --output-dir DIR。
保留Python：pandas元数据解析和本地gzip矩阵读取是既有复现路径，本轮不改数据字段或排序规则。
"""

from __future__ import annotations

import argparse
import gzip
import importlib
import re
from pathlib import Path

import pandas as pd

data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # 保留四个命令行参数及默认目录；GEO元数据只读本地缓存。
    parser = argparse.ArgumentParser(description="Build sample and library manifests for the GBM single-cell datasets.")
    parser.add_argument(
        "--project-root",
        type    = Path,
        default = Path(__file__).resolve().parents[1],
        help    = "Project root containing GSE174554 and GSE274546.",
    )
    parser.add_argument(
        "--geo-cache-dir",
        type    = Path,
        default = Path(__file__).resolve().parents[1] / "tmp_geo_gsm",
        help    = "Directory containing cached GEO sample text pages for GSE174554.",
    )
    parser.add_argument(
        "--output-dir",
        type    = Path,
        default = Path(__file__).resolve().parents[1] / "results" / "manifests",
        help    = "Directory where manifest TSV files will be written.",
    )
    parser.add_argument(
        "--gse274546-cache-dir",
        type    = Path,
        default = Path(__file__).resolve().parents[1] / "GSE274546",
        help    = "Directory containing standardized MatrixMarket files for GSE274546.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.geo_cache_dir.exists():
        raise FileNotFoundError(
            f"Missing GEO cache directory: {args.geo_cache_dir}. Fetch the GSM text pages before building manifests."
        )

    # 每行一个GSE274546样本；GEO时期标签、患者编号和源文件元数据按原规则提取。
    annotation_path = args.project_root / "GSE274546" / "GSE274546_GPL24676_sample_anno.csv"

    annotation_table = pd.read_csv(annotation_path)
    parsed_table = annotation_table["tissue:ch1"].str.extract(r"GBM single-nucleus,\s*([^,]+),\s*(.*)")
    annotation_table["patient_code"] = parsed_table[0]
    annotation_table["progression"] = parsed_table[1].map(data_functions.normalize_progression)
    annotation_table["tissue"] = annotation_table["progression"].map(data_functions.derive_tissue_from_progression)
    progression_order = {"Primary": 0, "1st Recurrent": 1, "Recurrent": 1, "2nd Recurrent": 2}
    annotation_table["progression_order"] = annotation_table["progression"].map(progression_order).astype("Int64")
    annotation_table["patient_numeric"] = annotation_table["title"].str.extract(r"Patient\s+(\d+)").astype("Int64")
    annotation_table["timepoint_numeric"] = annotation_table["title"].str.extract(r"Timepoint\s+(\d+)").astype("Int64")
    annotation_table["local_sample_id"] = (
        "P" + annotation_table["patient_numeric"].astype(str) + "T" + annotation_table["timepoint_numeric"].astype(str)
    )
    annotation_table["library_id"] = annotation_table["local_sample_id"]

    # 以文件名匹配源RDS文库；目录不存在时仍允许使用后续标准MatrixMarket输入。
    raw_dir = args.project_root / "GSE274546" / "GSE274546_RAW"
    file_rows = []
    if raw_dir.exists():
        for rds_path in sorted(raw_dir.glob("*.RDS.gz")):
            match = re.match(r"(?P<gsm>GSM\d+)_(?P<sample_id>[^_]+)_umi_counts\.RDS\.gz", rds_path.name)
            if not match:
                continue
            file_rows.append(
                {
                    "gsm_id": match.group("gsm"),
                    "local_sample_id": match.group("sample_id"),
                    "library_id": match.group("sample_id"),
                    "rds_gz_path": str(rds_path.resolve()),
                    "compressed_size_mb": round(rds_path.stat().st_size / 1024 / 1024, 3),
                }
            )
    file_table = pd.DataFrame(file_rows)

    # 按公开字段补齐队列和患者标识；缺RDS信息保留NA。
    gse274546_source_table_total = annotation_table.rename(columns={"geo_accession": "gsm_id", "title": "sample_title"})
    gse274546_source_table_total = gse274546_source_table_total.assign(
        dataset     = "GSE274546",
        technology  = "snRNA-seq",
        sample_id   = gse274546_source_table_total["gsm_id"],
        patient_id  = "GSE274546_" + gse274546_source_table_total["patient_code"].astype(str),
        rds_gz_path = pd.NA,
        compressed_size_mb = pd.NA,
    )
    if not file_table.empty:
        gse274546_source_table_total = gse274546_source_table_total.drop(
            columns = ["local_sample_id", "library_id", "rds_gz_path", "compressed_size_mb"]
        )
        gse274546_source_table_total = gse274546_source_table_total.merge(
            file_table,
            on  = "gsm_id",
            how = "left",
            validate = "1:1",
        )

    # 输出列序属于公开接口；正则提取的0/1列分别对应患者和时期捕获组。
    source_columns = [
        "dataset", "technology", "gsm_id", "sample_id", "library_id", "local_sample_id", "patient_id", "patient_code",
        "patient_numeric", "timepoint_numeric", "progression", "tissue", "progression_order", "sample_title",
        "rds_gz_path", "compressed_size_mb", "tissue:ch1",
    ]
    gse274546_source_table_total = gse274546_source_table_total[source_columns].sort_values(
        ["patient_numeric", "patient_code", "timepoint_numeric", "gsm_id"], kind="stable"
    )

    # GSE174554每行一个文库；按GSM连接已缓存GEO元数据，再汇总技术重复。
    sample_map_table = pd.read_csv(args.project_root / "GSE174554" / "GBM_GSM_samples.csv")
    sample_map_table["Sample_Name"] = sample_map_table["Sample_Name"].astype(str).str.strip()
    sample_map_table = sample_map_table.rename(columns={"GSM_ID": "gsm_id", "Sample_Name": "sample_name_csv"})

    geo_rows = [data_functions.parse_geo_sample_text(path) for path in sorted(args.geo_cache_dir.glob("GSM*.txt"))]
    geo_metadata_table = pd.DataFrame(geo_rows)

    # 肿瘤/正常注释每行一个barcode，按样本汇总；它不参与此处的髓系筛选。
    metadata_path = args.project_root / "GSE174554" / "GSE174554_Tumor_normal_metadata.txt.gz"
    tumor_normal_rows: list[dict[str, object]] = []
    with gzip.open(metadata_path, "rt") as handle:
        next(handle, None)
        for line in handle:
            line = line.strip()
            if not line:
                continue
            barcode, label = line.rsplit(" ", 1)
            sample_short = barcode.split("_", 1)[0]
            tumor_normal_rows.append({"sample_short": sample_short, "tumor_normal_label": label})

    if not tumor_normal_rows:
        tumor_normal_table = pd.DataFrame(columns=["sample_short", "tumor_cells_annotated", "normal_cells_annotated"])

    else:
        tumor_normal_table = (
            pd.DataFrame(tumor_normal_rows)
            .groupby(["sample_short", "tumor_normal_label"])
            .size()
            .unstack(fill_value=0)
            .rename(columns={"Tumor": "tumor_cells_annotated", "Normal": "normal_cells_annotated"})
            .reset_index()
        )
        for column in ["tumor_cells_annotated", "normal_cells_annotated"]:
            if column not in tumor_normal_table.columns:
                tumor_normal_table[column] = 0
    tumor_normal_table = tumor_normal_table[["sample_short", "tumor_cells_annotated", "normal_cells_annotated"]]

    # 扫描实际矩阵文件，每个技术重复保留一行和原matrix/features/barcodes路径。
    raw_dir = args.project_root / "GSE174554" / "GSE174554_RAW"
    library_rows: list[dict[str, object]] = []
    for matrix_path in sorted(raw_dir.glob("*_matrix.mtx.gz")):
        match = re.match(r"(?P<gsm>GSM\d+)_(?P<core>.+)_matrix\.mtx\.gz", matrix_path.name)
        if not match:
            continue
        gsm_id = match.group("gsm")
        library_core = match.group("core")
        if library_core.endswith("_batch2"):
            sample_short = library_core[: -len("_batch2")]
            tech_rep = "batch2"
        else:
            sample_short = library_core
            tech_rep = "main"

        n_genes, n_cells, nnz = data_functions.parse_matrix_header(matrix_path)
        base_prefix = matrix_path.name[: -len("matrix.mtx.gz")]
        features_path = raw_dir / f"{base_prefix}features.tsv.gz"
        barcodes_path = raw_dir / f"{base_prefix}barcodes.tsv.gz"
        library_rows.append(
            {
                "dataset": "GSE174554",
                "technology": "snRNA-seq",
                "gsm_id": gsm_id,
                "sample_id": gsm_id,
                "library_id": base_prefix[:-1],
                "sample_short": sample_short,
                "tech_rep": tech_rep,
                "matrix_path": str(matrix_path.resolve()),
                "features_path": str(features_path.resolve()),
                "barcodes_path": str(barcodes_path.resolve()),
                "n_genes": n_genes,
                "n_cells": n_cells,
                "nnz": nnz,
            }
        )

    # 文库逐次连接三个不同来源的元数据，保留原连接类型、校验和稳定排序。
    gse174554_library_table_total = pd.DataFrame(library_rows)
    gse174554_library_table_total = gse174554_library_table_total.merge(
        sample_map_table, on="gsm_id", how="left", validate="m:1"
    )
    gse174554_library_table_total = gse174554_library_table_total.merge(
        geo_metadata_table, on="gsm_id", how="left", validate="m:1"
    )
    gse174554_library_table_total = gse174554_library_table_total.merge(
        tumor_normal_table, on="sample_short", how="left"
    )
    gse174554_library_table_total["tissue"] = gse174554_library_table_total["progression"].map(
        data_functions.derive_tissue_from_progression
    )
    gse174554_library_table_total["progression_order"] = (
        gse174554_library_table_total["progression"].map(progression_order).astype("Int64")
    )
    gse174554_library_table_total["patient_id"] = gse174554_library_table_total["pair_id"].map(
        lambda value: f"GSE174554_pair{int(value)}" if pd.notna(value) else pd.NA
    )
    gse174554_library_table_total = gse174554_library_table_total.sort_values(
        ["pair_id", "progression_order", "gsm_id", "tech_rep", "library_id"],
        kind = "stable",
    )

    # 每个样本的元数据取第一个非缺失值；同一函数用于下列多个聚合列。

    gse174554_sample_table = (
        gse174554_library_table_total.groupby("gsm_id", dropna=False)
        .agg(
            dataset                = ("dataset", "first"),
            technology             = ("technology", "first"),
            sample_id              = ("sample_id", "first"),
            sample_short           = ("sample_short", "first"),
            sample_title           = ("sample_title", data_functions.first_nonmissing),
            sample_name_csv        = ("sample_name_csv", data_functions.first_nonmissing),
            diagnosis              = ("diagnosis", data_functions.first_nonmissing),
            progression            = ("progression", data_functions.first_nonmissing),
            tissue                 = ("tissue", data_functions.first_nonmissing),
            progression_order      = ("progression_order", data_functions.first_nonmissing),
            pair_label             = ("pair_label", data_functions.first_nonmissing),
            pair_id                = ("pair_id", data_functions.first_nonmissing),
            patient_id             = ("patient_id", data_functions.first_nonmissing),
            age                    = ("age", data_functions.first_nonmissing),
            gender                 = ("gender", data_functions.first_nonmissing),
            n_libraries            = ("library_id", "count"),
            total_cells            = ("n_cells", "sum"),
            total_nnz              = ("nnz", "sum"),
            tech_reps              = ("tech_rep", lambda values: ",".join(sorted(set(map(str, values))))),
            tumor_cells_annotated  = ("tumor_cells_annotated", "max"),
            normal_cells_annotated = ("normal_cells_annotated", "max"),
        )
        .reset_index()
        .sort_values(["pair_id", "progression_order", "gsm_id"], kind="stable")
    )

    # 标准矩阵每行一个文库；检查三个输入文件，并从MatrixMarket表头读取维度。
    standard_rows: list[dict[str, object]] = []
    for row in gse274546_source_table_total.to_dict("records"):
        sample_dir = args.gse274546_cache_dir / str(row["library_id"])
        matrix_path = sample_dir / "matrix.mtx.gz"
        features_path = sample_dir / "genes.tsv.gz"
        barcodes_path = sample_dir / "barcodes.tsv.gz"
        missing = [str(p) for p in [matrix_path, features_path, barcodes_path] if not p.exists()]
        if missing:
            raise FileNotFoundError(
                f"Missing standardized files for GSE274546. Run the RDS-to-sparse conversion first. Missing: {missing}"
            )
        n_genes, n_cells, nnz = data_functions.parse_matrix_header(matrix_path)
        has_rds_source = pd.notna(row["rds_gz_path"])
        source_format = "RDS.gz" if has_rds_source else "MatrixMarket"
        source_count_path = row["rds_gz_path"] if has_rds_source else str(matrix_path.resolve())
        standard_format = "matrix.mtx.gz + features/genes.tsv.gz + barcodes.tsv.gz"
        standard_rows.append(
            {
                "dataset": row["dataset"],
                "technology": row["technology"],
                "gsm_id": row["gsm_id"],
                "sample_id": row["sample_id"],
                "library_id": row["library_id"],
                "local_sample_id": row["local_sample_id"],
                "patient_id": row["patient_id"],
                "patient_code": row["patient_code"],
                "patient_numeric": row["patient_numeric"],
                "timepoint_numeric": row["timepoint_numeric"],
                "progression": row["progression"],
                "tissue": row["tissue"],
                "progression_order": row["progression_order"],
                "sample_title": row["sample_title"],
                "compressed_size_mb": row["compressed_size_mb"],
                "tissue:ch1": row["tissue:ch1"],
                "source_format": source_format,
                "standard_format": standard_format,
                "source_count_path": source_count_path,
                "matrix_path": str(matrix_path.resolve()),
                "features_path": str(features_path.resolve()),
                "barcodes_path": str(barcodes_path.resolve()),
                "n_genes": n_genes,
                "n_cells": n_cells,
                "nnz": nnz,
                "pair_id": pd.NA,
                "tech_rep": "main",
            }
        )
    gse274546_standard_table = pd.DataFrame(standard_rows).sort_values(
        ["patient_numeric", "patient_code", "timepoint_numeric", "gsm_id"], kind="stable"
    )

    # 两个队列按公开列序一次性拼接，总表稳定排序不变。
    gse174554_standard_table = gse174554_library_table_total.assign(
        source_format     = "10x MatrixMarket",
        standard_format   = standard_format,
        source_count_path = gse174554_library_table_total["matrix_path"],
        patient_code      = pd.NA,
        timepoint_numeric = pd.NA,
    )

    # 缺失字段按原合同补NA；tech_rep缺失时用main，不删除下游需要的列。
    for library_table in [gse274546_standard_table, gse174554_standard_table]:
        if "local_sample_id" not in library_table.columns:
            library_table["local_sample_id"] = pd.NA
        if "pair_id" not in library_table.columns:
            library_table["pair_id"] = pd.NA
        if "patient_numeric" not in library_table.columns:
            library_table["patient_numeric"] = pd.NA
        if "sample_title" not in library_table.columns:
            library_table["sample_title"] = pd.NA
        if "compressed_size_mb" not in library_table.columns:
            library_table["compressed_size_mb"] = pd.NA
        if "tech_rep" not in library_table.columns:
            library_table["tech_rep"] = "main"

    standard_columns = [
        "dataset", "technology", "gsm_id", "sample_id", "library_id", "local_sample_id", "patient_id", "patient_code",
        "patient_numeric", "pair_id", "timepoint_numeric", "progression", "tissue", "progression_order",
        "sample_title", "compressed_size_mb", "tech_rep", "source_format", "standard_format", "source_count_path",
        "matrix_path", "features_path", "barcodes_path", "n_genes", "n_cells", "nnz",
    ]
    standard_library_table_total = pd.concat(
        [
            gse274546_standard_table[standard_columns],
            gse174554_standard_table[standard_columns],
        ],
        ignore_index = True,
    )
    standard_library_table_total = standard_library_table_total.sort_values(
        ["dataset", "pair_id", "patient_code", "timepoint_numeric", "gsm_id", "tech_rep"], kind="stable"
    )

    # 汇总每个队列的样本数和时期分布；notes为已公开输出字段，保留既有文字。
    gse274546_note = (
        "Current folder uses standardized MatrixMarket sparse files."
        if "source_format" in gse274546_standard_table.columns
        and set(gse274546_standard_table["source_format"].dropna()) == {"MatrixMarket"}
        else "Standardized to MatrixMarket sparse format from RDS.gz source files."
    )
    batch2_samples = gse174554_library_table_total.loc[gse174554_library_table_total["tech_rep"].eq("batch2"), "gsm_id"]
    n_batch2_samples = int(batch2_samples.nunique())
    summary_rows = [
        {
            "dataset": "GSE274546",
            "entity_level": "library",
            "n_rows": int(len(gse274546_standard_table)),
            "n_unique_patients_or_pairs": int(gse274546_standard_table["patient_code"].nunique()),
            "n_primary": int((gse274546_standard_table["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse274546_standard_table["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": int((gse274546_standard_table["progression"] == "2nd Recurrent").sum()),
            "total_cells_local": int(gse274546_standard_table["n_cells"].sum()),
            "notes": gse274546_note,
        },
        {
            "dataset": "GSE174554",
            "entity_level": "library",
            "n_rows": int(len(gse174554_library_table_total)),
            "n_unique_patients_or_pairs": int(gse174554_library_table_total["pair_id"].dropna().nunique()),
            "n_primary": int((gse174554_library_table_total["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse174554_library_table_total["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": 0,
            "total_cells_local": int(gse174554_library_table_total["n_cells"].sum()),
            "notes": (
                f"{len(gse174554_library_table_total)} libraries map to "
                f"{len(gse174554_sample_table)} GSM samples; "
                f"{n_batch2_samples} "
                "GSMs contain an extra batch2 library."
            ),
        },
        {
            "dataset": "GSE174554",
            "entity_level": "sample",
            "n_rows": int(len(gse174554_sample_table)),
            "n_unique_patients_or_pairs": int(gse174554_sample_table["pair_id"].dropna().nunique()),
            "n_primary": int((gse174554_sample_table["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse174554_sample_table["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": 0,
            "total_cells_local": int(gse174554_sample_table["total_cells"].sum()),
            "notes": "Human GBM snRNA-seq GSM-level summary recovered from GEO sample text pages.",
        },
    ]
    summary_table = pd.DataFrame(summary_rows)

    # 七份既有清单各写TSV及CSV；回读维度和字段顺序，避免产物不完整。
    for manifest_table, stem in (
        (gse274546_source_table_total, "gse274546_source_samples"),
        (gse274546_standard_table, "gse274546_samples"),
        (gse174554_library_table_total, "gse174554_libraries"),
        (gse174554_sample_table, "gse174554_samples"),
        (gse274546_standard_table, "gse274546_standard_libraries"),
        (standard_library_table_total, "gbm_standard_libraries"),
        (summary_table, "dataset_summary"),
    ):
        for suffix, separator in (("tsv", "\t"), ("csv", ",")):
            output_path = args.output_dir / f"{stem}.{suffix}"
            manifest_table.to_csv(output_path, sep=separator, index=False)
            saved_table = pd.read_csv(output_path, sep=separator)
            if (
                saved_table.shape != manifest_table.shape
                or saved_table.columns.tolist() != manifest_table.columns.tolist()
            ):
                raise RuntimeError(f"Manifest dimensions or columns changed during save: {output_path}")

    print(
        f"GSE274546 libraries: {len(gse274546_standard_table)}; "
        f"GSE174554 libraries: {len(gse174554_library_table_total)}"
    )
    print(f"GSE174554 samples: {len(gse174554_sample_table)}; combined libraries: {len(standard_library_table_total)}")
