#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


def parse_soft(path: Path) -> list[dict[str, object]]:
    opener = gzip.open if path.suffix == ".gz" else open
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"geo_accession": line.split("=", 1)[1].strip(), "supplementary_urls": []}
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["sample_title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, item = value.split(":", 1)
                    current[key.strip().lower().replace(" ", "_")] = item.strip()
            elif line.startswith("!Sample_supplementary_file"):
                url = line.split("=", 1)[1].strip().replace("ftp://", "https://")
                current["supplementary_urls"].append(url)
    if current is not None:
        samples.append(current)
    return samples


def file_name(url: str) -> str:
    return Path(urlparse(url).path).name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-cell-soft", type=Path, required=True)
    parser.add_argument("--spatial-soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sc_samples = parse_soft(args.single_cell_soft)
    sc_rows = []
    downloads = []
    for sample in sc_samples:
        classification = str(sample.get("2021_classification", ""))
        is_idhwt_gbm = classification.lower() == "glioblastoma, idh-wt"
        myeloid_urls = [url for url in sample["supplementary_urls"] if "_Myl_filtered_counts.h5" in url]
        row = {key: value for key, value in sample.items() if key != "supplementary_urls"}
        row["is_idhwt_gbm"] = is_idhwt_gbm
        row["myeloid_h5_url"] = myeloid_urls[0] if myeloid_urls else ""
        sc_rows.append(row)
        if is_idhwt_gbm and myeloid_urls:
            downloads.append(
                {
                    "dataset": "GSE278456",
                    "sample": sample["sample_title"],
                    "geo_accession": sample["geo_accession"],
                    "modality": "scRNA_myeloid",
                    "url": myeloid_urls[0],
                    "local_subdir": "source_metadata/GSE278456_myeloid_h5",
                    "local_name": file_name(myeloid_urls[0]),
                }
            )

    spatial_samples = parse_soft(args.spatial_soft)
    spatial_rows = []
    for sample in spatial_samples:
        row = {key: value for key, value in sample.items() if key != "supplementary_urls"}
        row["n_supplementary_files"] = len(sample["supplementary_urls"])
        spatial_rows.append(row)
        for url in sample["supplementary_urls"]:
            downloads.append(
                {
                    "dataset": "GSE276841",
                    "sample": sample["sample_title"],
                    "geo_accession": sample["geo_accession"],
                    "modality": "Visium",
                    "url": url,
                    "local_subdir": "source_metadata/GSE276841_spatial",
                    "local_name": file_name(url),
                }
            )

    sc_df = pd.DataFrame(sc_rows)
    spatial_df = pd.DataFrame(spatial_rows)
    download_df = pd.DataFrame(downloads)
    sc_df.to_csv(args.output_dir / "GSE278456_sample_manifest.csv", index=False)
    spatial_df.to_csv(args.output_dir / "GSE276841_sample_manifest.csv", index=False)
    download_df.to_csv(args.output_dir / "eligible_download_manifest.csv", index=False)
    aria_lines = []
    for row in download_df.itertuples(index=False):
        aria_lines.extend(
            [
                row.url,
                f"  dir={args.output_dir / row.local_subdir}",
                f"  out={row.local_name}",
            ]
        )
    (args.output_dir / "eligible_download_aria2.txt").write_text(
        "\n".join(aria_lines) + "\n", encoding="utf-8"
    )

    print(f"GSE278456 samples: {len(sc_df)}")
    print(f"Strict IDH-wildtype GBM with myeloid H5: {int(sc_df.is_idhwt_gbm.sum())}")
    print(f"GSE276841 spatial samples: {len(spatial_df)}")
    print(f"Files selected for download: {len(download_df)}")


if __name__ == "__main__":
    main()
