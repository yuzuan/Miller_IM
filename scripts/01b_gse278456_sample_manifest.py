#!/usr/bin/env python3
"""Build GSE278456 myeloid and GSE276841 spatial download manifests from GEO SOFT metadata.
Inputs: --single-cell-soft and --spatial-soft (SAMPLE, title, characteristics and supplementary_file).
Outputs under --output-dir: GSE278456_sample_manifest.csv, GSE276841_sample_manifest.csv,
    eligible_download_manifest.csv and eligible_download_aria2.txt.
Run: python scripts/01b_gse278456_sample_manifest.py --single-cell-soft FILE --spatial-soft FILE --output-dir DIR.
The original Python workflow uses pandas, gzip and URL parsing with the existing filters and file names.
"""

import argparse
import importlib
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # Accept the original CLI arguments. Read both downloaded SOFT inputs and generate the manifests.
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-cell-soft", type=Path, required=True)
    parser.add_argument("--spatial-soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Read the SOFT file into sample dictionaries in source order. Let the reader report a missing file.

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

        # Keep the first eligible myeloid H5 URL for each sample meeting the strict IDH-wildtype GBM criteria.
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

    # Keep all spatial samples and list their attachments in the original order.
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

    # Sample rows represent GEO samples and download rows represent attachments. Keep the published order and fields.
    single_cell_table = pd.DataFrame(single_cell_rows)
    spatial_table = pd.DataFrame(spatial_rows)
    download_table_total = pd.DataFrame(download_rows)
    single_cell_path = args.output_dir / "GSE278456_sample_manifest.csv"
    single_cell_table.to_csv(single_cell_path, index=False)
    spatial_path = args.output_dir / "GSE276841_sample_manifest.csv"
    spatial_table.to_csv(spatial_path, index=False)
    download_path = args.output_dir / "eligible_download_manifest.csv"
    download_table_total.to_csv(download_path, index=False)

    # Write three aria2 lines per attachment in manifest order with the same relative destinations.
    aria_lines = []
    for row in download_table_total.itertuples(index=False):
        aria_lines.extend([row.url, f"  dir={args.output_dir / row.local_subdir}", f"  out={row.local_name}"])

    aria_path = args.output_dir / "eligible_download_aria2.txt"
    aria_path.write_text("\n".join(aria_lines) + "\n", encoding="utf-8")

    # Read back the manifest dimensions and require a nonempty aria2 file.
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
