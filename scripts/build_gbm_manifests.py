#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


PROGRESSION_ORDER = {
    "Primary": 0,
    "1st Recurrent": 1,
    "Recurrent": 1,
    "2nd Recurrent": 2,
}
PAIR_RE = re.compile(r"#?\s*(\d+)")
STANDARD_FORMAT = "matrix.mtx.gz + features/genes.tsv.gz + barcodes.tsv.gz"


def normalize_progression(value: str | None) -> str | pd.NA:
    if value is None or pd.isna(value):
        return pd.NA
    value = str(value).strip()
    mapping = {
        "Primary": "Primary",
        "Recurrent": "1st Recurrent",
        "1st Recurrent": "1st Recurrent",
        "2nd Recurrent": "2nd Recurrent",
    }
    return mapping.get(value, value)


def derive_tissue_from_progression(value: str | None) -> str | pd.NA:
    progression = normalize_progression(value)
    if progression is pd.NA or pd.isna(progression):
        return pd.NA
    if progression == "Primary":
        return "Primary"
    if "Recurrent" in str(progression):
        return "Recurrent"
    return pd.NA


def parse_pair_id(value: str | None) -> int | pd.NA:
    if value is None or pd.isna(value):
        return pd.NA
    match = PAIR_RE.search(str(value))
    return int(match.group(1)) if match else pd.NA


def parse_gse274546(project_root: Path) -> pd.DataFrame:
    anno_path = project_root / "GSE274546" / "GSE274546_GPL24676_sample_anno.csv"
    raw_dir = project_root / "GSE274546" / "GSE274546_RAW"

    anno = pd.read_csv(anno_path)
    parsed = anno["tissue:ch1"].str.extract(r"GBM single-nucleus,\s*([^,]+),\s*(.*)")
    anno["patient_code"] = parsed[0]
    anno["progression"] = parsed[1].map(normalize_progression)
    anno["tissue"] = anno["progression"].map(derive_tissue_from_progression)
    anno["progression_order"] = anno["progression"].map(PROGRESSION_ORDER).astype("Int64")
    anno["patient_numeric"] = anno["title"].str.extract(r"Patient\s+(\d+)").astype("Int64")
    anno["timepoint_numeric"] = anno["title"].str.extract(r"Timepoint\s+(\d+)").astype("Int64")
    anno["local_sample_id"] = (
        "P"
        + anno["patient_numeric"].astype(str)
        + "T"
        + anno["timepoint_numeric"].astype(str)
    )
    anno["library_id"] = anno["local_sample_id"]

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
    files = pd.DataFrame(file_rows)

    df = anno.rename(columns={"geo_accession": "gsm_id", "title": "sample_title"}).assign(
        dataset="GSE274546",
        technology="snRNA-seq",
        sample_id=lambda x: x["gsm_id"],
        patient_id=lambda x: "GSE274546_" + x["patient_code"].astype(str),
        rds_gz_path=pd.NA,
        compressed_size_mb=pd.NA,
    )
    if not files.empty:
        df = df.drop(columns=["local_sample_id", "library_id", "rds_gz_path", "compressed_size_mb"]).merge(
            files,
            on="gsm_id",
            how="left",
            validate="1:1",
        )

    ordered_cols = [
        "dataset",
        "technology",
        "gsm_id",
        "sample_id",
        "library_id",
        "local_sample_id",
        "patient_id",
        "patient_code",
        "patient_numeric",
        "timepoint_numeric",
        "progression",
        "tissue",
        "progression_order",
        "sample_title",
        "rds_gz_path",
        "compressed_size_mb",
        "tissue:ch1",
    ]
    return df[ordered_cols].sort_values(
        ["patient_numeric", "patient_code", "timepoint_numeric", "gsm_id"],
        kind="stable",
    )


def parse_geo_sample_text(text_path: Path) -> dict[str, object]:
    meta: dict[str, object] = {
        "gsm_id": text_path.stem,
        "sample_title": pd.NA,
        "diagnosis": pd.NA,
        "progression": pd.NA,
        "tissue": pd.NA,
        "pair_label": pd.NA,
        "pair_id": pd.NA,
        "age": pd.NA,
        "gender": pd.NA,
        "source_name": pd.NA,
    }
    for raw_line in text_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if raw_line.startswith("!Sample_title = "):
            meta["sample_title"] = raw_line.split(" = ", 1)[1].strip()
        elif raw_line.startswith("!Sample_source_name_ch1 = "):
            meta["source_name"] = raw_line.split(" = ", 1)[1].strip()
        elif raw_line.startswith("!Sample_characteristics_ch1 = "):
            value = raw_line.split(" = ", 1)[1].strip()
            if ":" not in value:
                continue
            key, detail = value.split(":", 1)
            key = key.strip().lower()
            detail = detail.strip()
            if key == "progression":
                meta["progression"] = normalize_progression(detail)
                meta["tissue"] = derive_tissue_from_progression(detail)
            elif key == "diagnosis":
                meta["diagnosis"] = detail
            elif key == "pair#":
                meta["pair_label"] = detail
                meta["pair_id"] = parse_pair_id(detail)
            elif key == "age":
                meta["age"] = detail
            elif key == "gender":
                meta["gender"] = detail
    return meta


def parse_tumor_normal_counts(metadata_path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with gzip.open(metadata_path, "rt") as handle:
        next(handle, None)
        for line in handle:
            line = line.strip()
            if not line:
                continue
            barcode, label = line.rsplit(" ", 1)
            sample_short = barcode.split("_", 1)[0]
            rows.append({"sample_short": sample_short, "tumor_normal_label": label})

    if not rows:
        return pd.DataFrame(columns=["sample_short", "tumor_cells_annotated", "normal_cells_annotated"])

    counts = (
        pd.DataFrame(rows)
        .groupby(["sample_short", "tumor_normal_label"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={"Tumor": "tumor_cells_annotated", "Normal": "normal_cells_annotated"})
        .reset_index()
    )
    for column in ["tumor_cells_annotated", "normal_cells_annotated"]:
        if column not in counts.columns:
            counts[column] = 0
    return counts[["sample_short", "tumor_cells_annotated", "normal_cells_annotated"]]


def parse_matrix_header(matrix_path: Path) -> tuple[int, int, int]:
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_rows, n_cols, nnz = (int(x) for x in line.strip().split())
            return n_rows, n_cols, nnz
    raise ValueError(f"Could not parse Matrix Market header: {matrix_path}")


def build_gse274546_standard_libraries(gse274546: pd.DataFrame, cache_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in gse274546.to_dict("records"):
        sample_dir = cache_dir / str(row["library_id"])
        matrix_path = sample_dir / "matrix.mtx.gz"
        features_path = sample_dir / "genes.tsv.gz"
        barcodes_path = sample_dir / "barcodes.tsv.gz"
        missing = [str(p) for p in [matrix_path, features_path, barcodes_path] if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing standardized files for GSE274546. "
                f"Run the RDS-to-sparse conversion first. Missing: {missing}"
            )
        n_genes, n_cells, nnz = parse_matrix_header(matrix_path)
        has_rds_source = pd.notna(row["rds_gz_path"])
        source_format = "RDS.gz" if has_rds_source else "MatrixMarket"
        source_count_path = row["rds_gz_path"] if has_rds_source else str(matrix_path.resolve())
        rows.append(
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
                "standard_format": STANDARD_FORMAT,
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
    return pd.DataFrame(rows).sort_values(
        ["patient_numeric", "patient_code", "timepoint_numeric", "gsm_id"],
        kind="stable",
    )


def build_standard_library_manifest(
    gse274546_standard: pd.DataFrame,
    gse174554_libraries: pd.DataFrame,
) -> pd.DataFrame:
    gse174554_standard = gse174554_libraries.assign(
        source_format="10x MatrixMarket",
        standard_format=STANDARD_FORMAT,
        source_count_path=lambda x: x["matrix_path"],
        patient_code=pd.NA,
        timepoint_numeric=pd.NA,
    )

    ordered_cols = [
        "dataset",
        "technology",
        "gsm_id",
        "sample_id",
        "library_id",
        "local_sample_id",
        "patient_id",
        "patient_code",
        "patient_numeric",
        "pair_id",
        "timepoint_numeric",
        "progression",
        "tissue",
        "progression_order",
        "sample_title",
        "compressed_size_mb",
        "tech_rep",
        "source_format",
        "standard_format",
        "source_count_path",
        "matrix_path",
        "features_path",
        "barcodes_path",
        "n_genes",
        "n_cells",
        "nnz",
    ]

    for frame in [gse274546_standard, gse174554_standard]:
        if "local_sample_id" not in frame.columns:
            frame["local_sample_id"] = pd.NA
        if "pair_id" not in frame.columns:
            frame["pair_id"] = pd.NA
        if "patient_numeric" not in frame.columns:
            frame["patient_numeric"] = pd.NA
        if "sample_title" not in frame.columns:
            frame["sample_title"] = pd.NA
        if "compressed_size_mb" not in frame.columns:
            frame["compressed_size_mb"] = pd.NA
        if "tech_rep" not in frame.columns:
            frame["tech_rep"] = "main"

    combined = pd.concat(
        [
            gse274546_standard[ordered_cols],
            gse174554_standard[ordered_cols],
        ],
        ignore_index=True,
    )
    return combined.sort_values(
        ["dataset", "pair_id", "patient_code", "timepoint_numeric", "gsm_id", "tech_rep"],
        kind="stable",
    )


def write_table(df: pd.DataFrame, output_dir: Path, stem: str) -> None:
    df.to_csv(output_dir / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(output_dir / f"{stem}.csv", index=False)


def build_gse174554_manifests(project_root: Path, geo_cache_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_dir = project_root / "GSE174554" / "GSE174554_RAW"
    sample_map = pd.read_csv(project_root / "GSE174554" / "GBM_GSM_samples.csv")
    sample_map["Sample_Name"] = sample_map["Sample_Name"].astype(str).str.strip()
    sample_map = sample_map.rename(columns={"GSM_ID": "gsm_id", "Sample_Name": "sample_name_csv"})

    geo_rows = [parse_geo_sample_text(path) for path in sorted(geo_cache_dir.glob("GSM*.txt"))]
    geo_meta = pd.DataFrame(geo_rows)

    tumor_normal = parse_tumor_normal_counts(project_root / "GSE174554" / "GSE174554_Tumor_normal_metadata.txt.gz")

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

        n_genes, n_cells, nnz = parse_matrix_header(matrix_path)
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

    libraries = pd.DataFrame(library_rows)
    libraries = libraries.merge(sample_map, on="gsm_id", how="left", validate="m:1")
    libraries = libraries.merge(geo_meta, on="gsm_id", how="left", validate="m:1")
    libraries = libraries.merge(tumor_normal, on="sample_short", how="left")
    libraries["tissue"] = libraries["progression"].map(derive_tissue_from_progression)
    libraries["progression_order"] = libraries["progression"].map(PROGRESSION_ORDER).astype("Int64")
    libraries["patient_id"] = libraries["pair_id"].map(
        lambda x: f"GSE174554_pair{int(x)}" if pd.notna(x) else pd.NA
    )
    libraries = libraries.sort_values(
        ["pair_id", "progression_order", "gsm_id", "tech_rep", "library_id"],
        kind="stable",
    )

    first_value = lambda values: next((x for x in values if pd.notna(x)), pd.NA)

    sample_manifest = (
        libraries.groupby("gsm_id", dropna=False)
        .agg(
            dataset=("dataset", "first"),
            technology=("technology", "first"),
            sample_id=("sample_id", "first"),
            sample_short=("sample_short", "first"),
            sample_title=("sample_title", first_value),
            sample_name_csv=("sample_name_csv", first_value),
            diagnosis=("diagnosis", first_value),
            progression=("progression", first_value),
            tissue=("tissue", first_value),
            progression_order=("progression_order", first_value),
            pair_label=("pair_label", first_value),
            pair_id=("pair_id", first_value),
            patient_id=("patient_id", first_value),
            age=("age", first_value),
            gender=("gender", first_value),
            n_libraries=("library_id", "count"),
            total_cells=("n_cells", "sum"),
            total_nnz=("nnz", "sum"),
            tech_reps=("tech_rep", lambda x: ",".join(sorted(set(map(str, x))))),
            tumor_cells_annotated=("tumor_cells_annotated", "max"),
            normal_cells_annotated=("normal_cells_annotated", "max"),
        )
        .reset_index()
        .sort_values(["pair_id", "progression_order", "gsm_id"], kind="stable")
    )

    return libraries, sample_manifest


def build_summary(gse274546: pd.DataFrame, gse174554_libraries: pd.DataFrame, gse174554_samples: pd.DataFrame) -> pd.DataFrame:
    gse274546_note = (
        "Current folder uses standardized MatrixMarket sparse files."
        if "source_format" in gse274546.columns and set(gse274546["source_format"].dropna()) == {"MatrixMarket"}
        else "Standardized to MatrixMarket sparse format from RDS.gz source files."
    )
    rows = [
        {
            "dataset": "GSE274546",
            "entity_level": "library",
            "n_rows": int(len(gse274546)),
            "n_unique_patients_or_pairs": int(gse274546["patient_code"].nunique()),
            "n_primary": int((gse274546["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse274546["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": int((gse274546["progression"] == "2nd Recurrent").sum()),
            "total_cells_local": int(gse274546["n_cells"].sum()),
            "notes": gse274546_note,
        },
        {
            "dataset": "GSE174554",
            "entity_level": "library",
            "n_rows": int(len(gse174554_libraries)),
            "n_unique_patients_or_pairs": int(gse174554_libraries["pair_id"].dropna().nunique()),
            "n_primary": int((gse174554_libraries["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse174554_libraries["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": 0,
            "total_cells_local": int(gse174554_libraries["n_cells"].sum()),
            "notes": "91 libraries map to 81 GSM samples; 10 GSMs contain an extra batch2 library.",
        },
        {
            "dataset": "GSE174554",
            "entity_level": "sample",
            "n_rows": int(len(gse174554_samples)),
            "n_unique_patients_or_pairs": int(gse174554_samples["pair_id"].dropna().nunique()),
            "n_primary": int((gse174554_samples["progression"] == "Primary").sum()),
            "n_first_recurrent": int((gse174554_samples["progression"] == "1st Recurrent").sum()),
            "n_second_recurrent": 0,
            "total_cells_local": int(gse174554_samples["total_cells"].sum()),
            "notes": "Human GBM snRNA-seq GSM-level summary recovered from GEO sample text pages.",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build sample and library manifests for the GBM single-cell datasets.")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Project root containing GSE174554 and GSE274546.",
    )
    parser.add_argument(
        "--geo-cache-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tmp_geo_gsm",
        help="Directory containing cached GEO sample text pages for GSE174554.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "manifests",
        help="Directory where manifest TSV files will be written.",
    )
    parser.add_argument(
        "--gse274546-cache-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "GSE274546",
        help="Directory containing standardized MatrixMarket files for GSE274546.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.geo_cache_dir.exists():
        raise FileNotFoundError(
            f"Missing GEO cache directory: {args.geo_cache_dir}. Fetch the GSM text pages before building manifests."
        )

    gse274546_source = parse_gse274546(args.project_root)
    gse174554_libraries, gse174554_samples = build_gse174554_manifests(args.project_root, args.geo_cache_dir)
    gse274546_standard = build_gse274546_standard_libraries(gse274546_source, args.gse274546_cache_dir)
    standard_libraries = build_standard_library_manifest(gse274546_standard, gse174554_libraries)
    summary = build_summary(gse274546_standard, gse174554_libraries, gse174554_samples)

    write_table(gse274546_source, args.output_dir, "gse274546_source_samples")
    write_table(gse274546_standard, args.output_dir, "gse274546_samples")
    write_table(gse174554_libraries, args.output_dir, "gse174554_libraries")
    write_table(gse174554_samples, args.output_dir, "gse174554_samples")
    write_table(gse274546_standard, args.output_dir, "gse274546_standard_libraries")
    write_table(standard_libraries, args.output_dir, "gbm_standard_libraries")
    write_table(summary, args.output_dir, "dataset_summary")

    print("Wrote manifests:")
    print(f"  - {args.output_dir / 'gse274546_source_samples.csv'}")
    print(f"  - {args.output_dir / 'gse274546_samples.csv'}")
    print(f"  - {args.output_dir / 'gse174554_libraries.csv'}")
    print(f"  - {args.output_dir / 'gse174554_samples.csv'}")
    print(f"  - {args.output_dir / 'gse274546_standard_libraries.csv'}")
    print(f"  - {args.output_dir / 'gbm_standard_libraries.csv'}")
    print(f"  - {args.output_dir / 'dataset_summary.csv'}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
