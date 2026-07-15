#!/usr/bin/env python3
"""Validate the frozen public Miller-IM code and figure source-data release."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / "source_data"
MANIFEST = ROOT / "manifests/source_data_sha256.tsv"
EXPECTED_GENES = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
TEXT_SUFFIXES = {".py", ".R", ".sh", ".md", ".tsv", ".cff", ".txt"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_files() -> list[Path]:
    return sorted(path for path in SOURCE_ROOT.rglob("*") if path.is_file())


def build_manifest() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "path": [path.relative_to(ROOT).as_posix() for path in source_files()],
            "bytes": [path.stat().st_size for path in source_files()],
            "sha256": [sha256(path) for path in source_files()],
        }
    )


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, compression="infer")


def validate() -> None:
    genes = pd.read_csv(ROOT / "config/miller_im_gene_set.tsv", sep="\t")
    check(genes["gene"].tolist() == EXPECTED_GENES, "The fixed 20-gene definition changed.")

    files = source_files()
    check(len(files) == 37, f"Expected 37 frozen Figure 1-5 source files, found {len(files)}.")
    for path in files:
        if path.name.endswith(".csv.gz"):
            with gzip.open(path, "rt") as handle:
                check(bool(handle.readline().strip()), f"Empty gzip CSV: {path}")
            read_csv(path)
        elif path.suffix == ".csv":
            read_csv(path)

    f1a = read_csv(SOURCE_ROOT / "Figure1/Figure1B_GSE174554_raw20_gsea_source.csv.gz").iloc[0]
    f1b = read_csv(SOURCE_ROOT / "Figure1/Figure1C_GSE274546_raw20_gsea_source.csv.gz").iloc[0]
    check(np.isclose(f1a["NES"], 2.3157685724298) and int(f1a["n_pairs"]) == 18, "Figure 1 GSE174554 mismatch.")
    check(np.isclose(f1b["NES"], 1.79932425427) and int(f1b["n_pairs"]) == 45, "Figure 1 GSE274546 mismatch.")

    f2 = read_csv(SOURCE_ROOT / "Figure2/Figure2_state_patient_waterfall_source.csv")
    mcg1 = f2.loc[f2["state"].eq("MCG1")].iloc[0]
    check(int(mcg1["n_positive"]) == 17 and int(mcg1["n_patients"]) == 18, "Figure 2 MCG1 mismatch.")

    f3 = read_csv(SOURCE_ROOT / "Figure3/Figure3A_geomx_paired_dumbbell_source.csv").iloc[0]
    check(int(f3["n_pairs"]) == 22 and int(f3["n_up"]) == 16, "Figure 3 GeoMx mismatch.")

    f4a = read_csv(SOURCE_ROOT / "Figure4/Figure4A_pdc_patient_delta_waterfall_source.csv").iloc[0]
    f4c = read_csv(SOURCE_ROOT / "Figure4/Figure4C_full_proteome_rank_landscape_source.csv").iloc[0]
    check(int(f4a["n_pairs"]) == 105 and int(f4a["n_recurrence_up"]) == 72, "Figure 4 paired score mismatch.")
    check(np.isclose(f4c["miller_im_nes"], 1.52795404604028), "Figure 4 NES mismatch.")
    check(np.isclose(f4c["miller_im_gsea_p"], 0.0490367775831874), "Figure 4 nominal P mismatch.")

    f5 = read_csv(SOURCE_ROOT / "Figure5/Figure5B_dual_cohort_raw20_enrichment_source.csv.gz")
    observed = f5.groupby("cohort", sort=True).first()
    check(np.isclose(observed.loc["GSE121810", "NES"], 1.89995642058625), "Figure 5 GSE121810 mismatch.")
    check(np.isclose(observed.loc["GSE154795", "NES"], 1.751, atol=0.01), "Figure 5 GSE154795 mismatch.")

    forbidden_paths = re.compile("/" + "Users/|/" + r"Volumes/|[A-Za-z]:\\")
    misleading_fdr = re.compile(r"(?:targeted|single-set)\s+(?:fgsea\s+)?FDR", re.IGNORECASE)
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts or path.suffix not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8")
        check(forbidden_paths.search(text) is None, f"Local absolute path in {path.relative_to(ROOT)}")
        check(misleading_fdr.search(text) is None, f"Misleading program-level multiplicity label in {path.relative_to(ROOT)}")

    current = build_manifest()
    check(MANIFEST.exists(), "Missing source-data SHA-256 manifest.")
    frozen = pd.read_csv(MANIFEST, sep="\t")
    check(current.equals(frozen), "Frozen source-data SHA-256 manifest does not match current files.")
    print(f"RELEASE_VALIDATION_PASS files={len(files)} figures=5")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    if args.write_manifest:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        build_manifest().to_csv(MANIFEST, sep="\t", index=False)
        print(f"WROTE_MANIFEST {MANIFEST}")
    validate()


if __name__ == "__main__":
    main()
