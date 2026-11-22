#!/usr/bin/env python3
"""Validate the frozen Miller-IM release tables.

Purpose: Check reported results, public data fields, readable files, and the release manifest.
Inputs: The fixed definition and IDH crosswalk in config/, Figures 1–5 in source_data/, and repository text.
Output: A validation summary. The --write-manifest option updates manifests/source_data_sha256.tsv.
Run: python scripts/10b_validate_release.py [--write-manifest]
Expected values and field names match the published release.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
sha256 = figure_style.sha256


def main() -> None:
    """Validate the release and regenerate the manifest only when explicitly requested."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()

    # Read back the requested manifest to check its rows and columns.
    if args.write_manifest:
        manifest_path = Path(__file__).resolve().parents[1] / "manifests/source_data_sha256.tsv"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_table = build_manifest()
        manifest_table.to_csv(manifest_path, sep="\t", index=False)
        written_manifest_table = pd.read_csv(manifest_path, sep="\t")
        check(manifest_table.equals(written_manifest_table), "Written source-data manifest differs from memory.")
        print(f"WROTE_MANIFEST {manifest_path}")

    validate()


def build_manifest() -> pd.DataFrame:
    """Return a path-sorted manifest of the frozen source file sizes and digests."""
    root = Path(__file__).resolve().parents[1]
    file_records = []

    # Keep the public path, bytes, and sha256 column order.
    for path in source_files():
        file_records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

    manifest_table = pd.DataFrame(file_records, columns=["path", "bytes", "sha256"])
    return manifest_table


def source_files() -> list[Path]:
    """Return all files in source_data/ in sorted order."""
    source_root = Path(__file__).resolve().parents[1] / "source_data"
    paths = []

    for path in source_root.rglob("*"):
        if path.is_file():
            paths.append(path)

    paths = sorted(paths)
    return paths


def validate() -> None:
    """Check the frozen release against the published validation targets."""

    # The configuration gives the author-ranked order of the 20 program genes.
    root = Path(__file__).resolve().parents[1]
    gene_table = pd.read_csv(root / "config/miller_im_gene_set.tsv", sep="\t")
    expected_genes = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3",
        "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]
    check(gene_table["gene"].tolist() == expected_genes, "The fixed 20-gene definition changed.")

    # The crosswalk records the 17 included pairs and the sample IDs for the one excluded pair.
    crosswalk_table = read_csv(root / "config/GSE174554_formal_pair_IDH_crosswalk.csv")
    required_crosswalk_columns = {
        "dataset", "pair_key", "primary_sample", "recurrent_sample", "primary_idh",
        "recurrent_idh", "include_idh_wildtype", "source_file_sha256",
    }
    check(
        required_crosswalk_columns.issubset(crosswalk_table.columns),
        "The GSE174554 IDH crosswalk is missing required columns.",
    )
    included = boolean_series(crosswalk_table["include_idh_wildtype"], "GSE174554 IDH crosswalk")
    excluded_table = crosswalk_table.loc[~included]
    check(len(crosswalk_table) == 18 and int(included.sum()) == 17, "GSE174554 IDH eligibility is not 17/18.")
    check(
        excluded_table["pair_key"].tolist() == ["GSE174554_pair33"],
        "The only excluded GSE174554 pair must be GSE174554_pair33.",
    )

    # Check the two sample identifiers and both IDH annotations for the excluded pair.
    pair33 = excluded_table.iloc[0]
    check(
        pair33["primary_sample"] == "SF8963"
        and pair33["recurrent_sample"] == "SF12165"
        and pair33["primary_idh"] == "IDH mutant"
        and pair33["recurrent_idh"] == "IDH mutant",
        "GSE174554_pair33 sample or IDH annotation changed.",
    )
    check(
        crosswalk_table["source_file_sha256"].nunique() == 1
        and crosswalk_table["source_file_sha256"].iloc[0]
        == "58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42",
        "The audited Wang supplementary-workbook hash changed.",
    )

    # Verify the 37 published source files, including complete reads of the compressed CSVs.
    files = source_files()
    check(len(files) == 37, f"Expected 37 frozen Figure 1-5 source files, found {len(files)}.")
    for path in files:
        if path.name.endswith(".csv.gz"):
            with gzip.open(path, "rt") as handle:
                check(bool(handle.readline().strip()), f"Empty gzip CSV: {path}")
            read_csv(path)
        elif path.suffix == ".csv":
            read_csv(path)

    # Figure 1A reports retained cells and pairs for each cohort.
    source_root = root / "source_data"
    design_table = read_csv(source_root / "Figure1/Figure1A_paired_cohort_design_source.csv")
    design_table = design_table.set_index("dataset")
    check(
        design_table["formal_pairs_threshold20"].astype(int).to_dict() == {"GSE174554": 17, "GSE274546": 45},
        "Figure 1 paired design mismatch.",
    )
    check(
        design_table["clean_myeloid"].astype(int).to_dict() == {"GSE174554": 12489, "GSE274546": 55568},
        "Figure 1 analyzed-cell counts mismatch.",
    )

    # Each curve row repeats the enrichment statistics. The first row supplies the comparison values.
    gse174554_table = read_csv(source_root / "Figure1/Figure1B_GSE174554_raw20_gsea_source.csv.gz")
    gse174554 = gse174554_table.iloc[0]
    check(
        np.isclose(gse174554["NES"], 2.4290789842056)
        and np.isclose(gse174554["nominal_P"], 4.00136234233155e-7)
        and int(gse174554["n_pairs"]) == 17,
        "Figure 1 GSE174554 mismatch.",
    )
    gse274546_table = read_csv(source_root / "Figure1/Figure1C_GSE274546_raw20_gsea_source.csv.gz")
    gse274546 = gse274546_table.iloc[0]
    check(
        np.isclose(gse274546["NES"], 1.79932425427) and int(gse274546["n_pairs"]) == 45,
        "Figure 1 GSE274546 mismatch.",
    )

    # Panels D and E share seven leading-edge genes. Panel E has 14 gene-by-cohort rows.
    rank_table = read_csv(source_root / "Figure1/Figure1D_shared_raw20_outliers_source.csv")
    shared = boolean_series(rank_table["shared_leading_edge"], "Figure 1D shared-leading-edge flag")
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}
    check(
        set(rank_table.loc[shared, "gene"]) == expected_shared_leading_edge,
        "Figure 1 shared leading-edge genes changed.",
    )
    effect_table = read_csv(source_root / "Figure1/Figure1E_shared_leading_edge_gene_effects_source.csv")
    check(
        len(effect_table) == 14
        and effect_table["gene"].nunique() == 7
        and set(effect_table["gene"]) == expected_shared_leading_edge
        and int((effect_table["FDR"] < 0.05).sum()) == 0,
        "Figure 1 gene-by-cohort audit mismatch.",
    )

    # The saved LOPO results record the omitted patient for each fold.
    lopo_table = read_csv(source_root / "Figure1/Figure1F_leave_one_patient_out_raw20_gsea.csv")
    lopo_table = lopo_table.loc[lopo_table["run_type"].eq("Leave-one-patient-out")].copy()
    check(
        lopo_table.groupby("dataset").size().to_dict() == {"GSE174554": 17, "GSE274546": 45},
        "Figure 1 LOPO fold counts mismatch.",
    )
    check(
        not lopo_table["omitted_patient"].eq("GSE174554_pair33").any(),
        "The excluded IDH-mutant pair appears in Figure 1 LOPO.",
    )

    # Use the frozen MCG1 patient counts to check the Figure 2 table.
    state_table = read_csv(source_root / "Figure2/Figure2_state_patient_waterfall_source.csv")
    mcg1 = state_table.loc[state_table["state"].eq("MCG1")].iloc[0]
    check(int(mcg1["n_positive"]) == 17 and int(mcg1["n_patients"]) == 18, "Figure 2 MCG1 mismatch.")

    # Check the Figure 3A paired total and the measured or clipped gene-by-patient entries in 3B.
    geomx_table = read_csv(source_root / "Figure3/Figure3A_geomx_paired_dumbbell_source.csv")
    geomx = geomx_table.iloc[0]
    check(int(geomx["n_pairs"]) == 22 and int(geomx["n_up"]) == 16, "Figure 3 GeoMx mismatch.")
    fingerprint_table = read_csv(source_root / "Figure3/Figure3B_shared_gene_patient_fingerprint_source.csv")
    measured = boolean_series(fingerprint_table["measured"], "Figure 3B measured flag")
    exceeds_color_limit = (
        fingerprint_table.loc[measured, "gene_delta"].abs() > fingerprint_table.loc[measured, "color_limit"]
    )
    clipped_count = int(exceeds_color_limit.sum())
    check(
        len(fingerprint_table) == 154 and int(measured.sum()) == 132 and clipped_count == 3,
        "Figure 3B 154/132/3 audit mismatch.",
    )

    # Check the frozen paired proteomic score and whole-proteome ranking statistics.
    protein_patient_table = read_csv(source_root / "Figure4/Figure4A_pdc_patient_delta_waterfall_source.csv")
    protein_patient = protein_patient_table.iloc[0]
    check(
        int(protein_patient["n_pairs"]) == 105 and int(protein_patient["n_recurrence_up"]) == 72,
        "Figure 4 paired score mismatch.",
    )
    protein_rank_table = read_csv(source_root / "Figure4/Figure4C_full_proteome_rank_landscape_source.csv")
    protein_rank = protein_rank_table.iloc[0]
    check(np.isclose(protein_rank["miller_im_nes"], 1.52795404604028), "Figure 4 NES mismatch.")
    check(np.isclose(protein_rank["miller_im_gsea_p"], 0.0490367775831874), "Figure 4 nominal P mismatch.")

    # The GSE154795 NES is checked against the reported value with an absolute tolerance of 0.01.
    enrichment_table = read_csv(source_root / "Figure5/Figure5B_dual_cohort_raw20_enrichment_source.csv.gz")
    cohort_table = enrichment_table.groupby("cohort", sort=True).first()
    check(np.isclose(cohort_table.loc["GSE121810", "NES"], 1.89995642058625), "Figure 5 GSE121810 mismatch.")
    check(
        np.isclose(cohort_table.loc["GSE154795", "NES"], 1.751, atol=0.01),
        "Figure 5 GSE154795 mismatch.",
    )

    # Split the path patterns so the source code itself can be included in the scan.
    forbidden_paths = re.compile("/" + "Users/|/" + r"Volumes/|[A-Za-z]:\\")
    misleading_fdr = re.compile(r"(?:targeted|single-set)\s+(?:fgsea\s+)?FDR", re.IGNORECASE)
    text_suffixes = {".py", ".R", ".sh", ".md", ".tsv", ".cff", ".txt"}
    for path in root.rglob("*"):
        if not path.is_file() or ".git" in path.parts or path.suffix not in text_suffixes:
            continue
        content = path.read_text(encoding="utf-8")
        check(forbidden_paths.search(content) is None, f"Local absolute path in {path.relative_to(root)}")
        check(
            misleading_fdr.search(content) is None,
            f"Misleading program-level multiplicity label in {path.relative_to(root)}",
        )

    # Also check headers and string columns in both plain and compressed CSVs.
    for path in tabular_files():
        character_table = read_csv(path)
        values = [str(column) for column in character_table.columns]
        for column in character_table.select_dtypes(include=["object", "string"]).columns:
            values.extend(character_table[column].dropna().astype(str).tolist())
        check(
            not any(forbidden_paths.search(value) for value in values),
            f"Local absolute path in tabular character field: {path.relative_to(root)}",
        )

    # Compare the current files to the frozen manifest by path, byte count, and digest.
    current_manifest_table = build_manifest()
    manifest_path = root / "manifests/source_data_sha256.tsv"
    check(manifest_path.exists(), "Missing source-data SHA-256 manifest.")
    frozen_manifest_table = pd.read_csv(manifest_path, sep="\t")
    check(
        current_manifest_table.equals(frozen_manifest_table),
        "Frozen source-data SHA-256 manifest does not match current files.",
    )
    print(f"RELEASE_VALIDATION_PASS files={len(files)} figures=5 idh_pairs=17+45")


def check(condition: bool, message: str) -> None:
    """Raise AssertionError when a release check fails."""
    if not condition:
        raise AssertionError(message)


def read_csv(path: Path) -> pd.DataFrame:
    """Read a plain or compressed CSV as a data frame."""
    source_table = pd.read_csv(path, compression="infer")
    return source_table


def boolean_series(series: pd.Series, label: str) -> pd.Series:
    """Return the parsed Boolean series and reject unexpected values."""
    normalized = series.astype(str).str.strip().str.lower()
    check(normalized.isin({"true", "false"}).all(), f"Invalid boolean values in {label}.")
    normalized = normalized.eq("true")
    return normalized


def tabular_files() -> list[Path]:
    """Return the repository CSV and gzip CSV paths outside .git in sorted order."""
    root = Path(__file__).resolve().parents[1]
    paths = []

    for path in root.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.name.endswith(".csv") or path.name.endswith(".csv.gz"):
            paths.append(path)

    paths = sorted(paths)
    return paths


if __name__ == "__main__":
    main()
