#!/usr/bin/env python3
# Purpose: reconstruct Miller raw20 state scores, patient effects and gene-by-state summaries from GSE278456 H5.
# Inputs: Step 30 sample and QC tables, plus the raw H5 files and author-marker workbook in the public data root.
# Outputs: GSE278456_*.csv and the compressed cell annotations in write/48_figure2_identity_state_rebuild/Figure2.
# Run: python scripts/06b_cell_identity_and_state.py from the repository root in the single-cell environment.
# Scanpy, AnnData and Matplotlib provide the released reconstruction and plotting implementation.
# Prepared input paths, output names and fields follow the original public analysis.

from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
sha256 = figure_style.sha256

def read_author_markers(
    marker_book: Path, raw20: list[str], myeloid_states: list[str]
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Read author markers, exclude raw20, and return marker sets and the audit. Raise for missing or sparse sets."""

    # Read worksheets in workbook order. Use the second header row when the first row lacks FDR.
    workbook = pd.ExcelFile(marker_book)
    markers: dict[str, list[str]] = {}
    rows: list[dict[str, object]] = []
    raw20_set = set(raw20)
    for sheet in workbook.sheet_names:
        data_table = pd.read_excel(marker_book, sheet_name=sheet)
        if "FDR" not in data_table.columns:
            data_table = pd.read_excel(marker_book, sheet_name=sheet, header=1)
        gene_col = data_table.columns[0]

        # Filter at FDR < 0.05 and Foldchange > 0, then take the first 40 in worksheet order before removing raw20.
        original = (
            data_table.loc[
                pd.to_numeric(data_table["FDR"], errors="coerce").lt(0.05)
                & pd.to_numeric(data_table["Foldchange"], errors="coerce").gt(0),
                gene_col,
            ]
            .astype(str)
            .tolist()[:40]
        )

        # Exclude raw20 genes from the selected 40. At least 30 markers must remain for state assignment.
        cleaned = [gene for gene in original if gene not in raw20_set]
        removed = [gene for gene in original if gene in raw20_set]
        if len(cleaned) < 30:
            raise RuntimeError(f"{sheet} has only {len(cleaned)} author markers after excluding raw20")
        markers[sheet] = cleaned
        rows.append(
            {
                "state": sheet,
                "n_original_top_markers": len(original),
                "n_after_raw20_exclusion": len(cleaned),
                "removed_raw20_genes": ";".join(removed),
                "markers_used": ";".join(cleaned),
            }
        )

    # Check all nine prespecified states. The audit lists removed genes and selected markers for each state.
    missing = sorted(set(myeloid_states) - set(markers))
    if missing:
        raise RuntimeError(f"The author-marker workbook lacks prespecified myeloid states: {missing}")
    marker_audit_table = pd.DataFrame(rows)
    result = (markers, marker_audit_table)
    return result

if __name__ == "__main__":

    # Use the original Scanpy settings and scoring implementation.
    sc.settings.verbosity = 0

    # Keep the Step 48 output paths for the primary tables, source tables and figures.
    root = Path(__file__).resolve().parents[1]

    write_root = root / "write/48_figure2_identity_state_rebuild/Figure2"

    write_root.mkdir(parents=True, exist_ok=True)

    source_out = write_root / "source_data"

    source_out.mkdir(parents=True, exist_ok=True)

    fig_root = root / "figures/48_figure2_identity_state_rebuild/Figure2"

    fig_root.mkdir(parents=True, exist_ok=True)

    # Read H5 files and the marker workbook under MILLER_IM_DATA_ROOT, which defaults to the project data directory.
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()

    marker_book = external_root / "write/30_dataset_rescue_search/source_metadata/NIHMS2115262-supplement-S2_S5_S6.xls"

    h5_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE278456_myeloid_h5"

    if not h5_root.exists() or not marker_book.exists():
        raise FileNotFoundError("The GSE278456 H5 files or the author-marker workbook are unavailable")

    # Keep the original 20 genes in Miller's published rank order. The symbols follow the source data.
    raw20 = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3",
        "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]

    # Use the original prespecified myeloid state set and order.
    myeloid_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]

    markers, marker_audit_table = read_author_markers(marker_book, raw20, myeloid_states)
    atlas_markers = markers

    step30 = root / "write/30_dataset_rescue_search"

    # Read the IDH-wild-type sample manifest with one row per GEO sample.
    atlas_manifest_table = pd.read_csv(step30 / "GSE278456_sample_manifest.csv")
    atlas_manifest_table["geo_accession"] = atlas_manifest_table["geo_accession"].astype(str)

    # Read the frozen Step 30 input and QC-pass cell counts for the sample-level checks below.
    atlas_existing_qc_table = pd.read_csv(step30 / "GSE278456_single_cell_qc.csv").set_index("geo_accession")
    atlas_paths = sorted(h5_root.glob("*.h5"))
    if len(atlas_paths) != 21:
        raise RuntimeError(f"Expected 21 frozen GSE278456 H5 files, found {len(atlas_paths)}")

    # Collect each sample's cell labels, state means, patient effects, gene expression, QC and input record.
    atlas_all_marker_genes = set().union(*atlas_markers.values())
    atlas_cell_frames: list[pd.DataFrame] = []
    atlas_state_score_frames: list[pd.DataFrame] = []
    atlas_effect_rows: list[dict[str, object]] = []
    atlas_gene_patient_rows: list[dict[str, object]] = []

    atlas_coverage_rows: list[dict[str, object]] = []
    atlas_locked_rows: list[dict[str, object]] = []
    atlas_qc_rows: list[dict[str, object]] = []

    # Process H5 files in filename order and match each to a unique IDH-wild-type manifest row.
    for atlas_path in atlas_paths:
        atlas_accession = atlas_path.name.split("_")[0]
        atlas_manifest_match_table = atlas_manifest_table.loc[atlas_manifest_table["geo_accession"].eq(atlas_accession)]
        if len(atlas_manifest_match_table) != 1:
            raise RuntimeError(f"{atlas_accession} does not match exactly one sample manifest row")
        atlas_manifest_row = atlas_manifest_match_table.iloc[0]
        if not bool(atlas_manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"Non-IDH-wild-type GBM in the frozen H5 set: {atlas_accession}")
        atlas_sample = str(atlas_manifest_row["sample_title"])

        # Read the 10x counts, make gene names unique and compute cell count and mitochondrial QC metrics.
        atlas_adata = sc.read_10x_h5(atlas_path)

        atlas_adata.var_names_make_unique()
        atlas_n_input = atlas_adata.n_obs
        atlas_n_genes = atlas_adata.n_vars
        atlas_adata.var["mt"] = atlas_adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(atlas_adata, qc_vars=["mt"], inplace=True, percent_top=None)

        # Apply the original QC bounds of 200–8,000 genes, at least 500 counts and at most 20% mitochondrial counts.
        atlas_keep = (
            atlas_adata.obs["n_genes_by_counts"].between(200, 8000)
            & atlas_adata.obs["total_counts"].ge(500)
            & atlas_adata.obs["pct_counts_mt"].le(20)
        )
        atlas_adata = atlas_adata[atlas_keep].copy()
        atlas_n_qc = atlas_adata.n_obs
        if atlas_accession not in atlas_existing_qc_table.index:
            raise RuntimeError(f"The frozen QC table lacks {atlas_accession}")
        atlas_old = atlas_existing_qc_table.loc[atlas_accession]

        if int(atlas_old["n_input"]) != atlas_n_input or int(atlas_old["n_qc_pass"]) != atlas_n_qc:
            raise RuntimeError(f"{atlas_accession} cell counts differ from the frozen QC counts")

        # Check coverage of the 19 measurable raw20 genes in each sample. AC253572.2 is the sole unmeasurable gene.
        atlas_present_raw20 = [atlas_gene for atlas_gene in raw20 if atlas_gene in atlas_adata.var_names]

        atlas_missing_raw20 = [atlas_gene for atlas_gene in raw20 if atlas_gene not in atlas_adata.var_names]
        if len(atlas_present_raw20) != 19 or atlas_missing_raw20 != ["AC253572.2"]:
            raise RuntimeError(
                f"{atlas_accession} has unexpected raw20 coverage: {len(atlas_present_raw20)}/20, "
                f"missing={atlas_missing_raw20}"
            )
        atlas_coverage_rows.append(
            {
                "geo_accession": atlas_accession,
                "sample": atlas_sample,
                "n_defined": len(raw20),
                "n_present": len(atlas_present_raw20),
                "present_genes": ";".join(atlas_present_raw20),
                "missing_genes": ";".join(atlas_missing_raw20),
            }
        )

        atlas_qc_rows.append(
            {
                "sample": atlas_sample,
                "geo_accession": atlas_accession,
                "n_input": atlas_n_input,
                "n_qc_pass": atlas_n_qc,
                "fraction_qc_pass": atlas_n_qc / atlas_n_input,
            }
        )
        atlas_locked_rows.append(
            {
                "geo_accession": atlas_accession,
                "sample": atlas_sample,
                "path": str(atlas_path),
                "sha256": sha256(atlas_path),
                "n_input_cells": atlas_n_input,
                "n_qc_pass_cells": atlas_n_qc,
                "n_genes": atlas_n_genes,
            }
        )

        # Keep the QC-passing raw counts for pseudobulk and prefix cell barcodes with their GSM accession.
        atlas_counts = atlas_adata.X.copy()
        atlas_raw20_indices = atlas_adata.var_names.get_indexer(atlas_present_raw20)
        atlas_adata.obs_names = [
            f"{atlas_accession}:{atlas_barcode}" for atlas_barcode in atlas_adata.obs_names.astype(str)
        ]

        # Normalize to 10,000 counts and apply log1p. Keep the original count matrix for pseudobulk.
        sc.pp.normalize_total(atlas_adata, target_sum=10000)
        sc.pp.log1p(atlas_adata)

        # Exclude raw20 and all author markers from the control pool. At least 5,000 neutral genes must remain.
        atlas_excluded_controls = set(raw20) | atlas_all_marker_genes
        atlas_neutral_control_pool = [
            atlas_gene for atlas_gene in atlas_adata.var_names if atlas_gene not in atlas_excluded_controls
        ]
        if len(atlas_neutral_control_pool) < 5000:
            raise RuntimeError(f"{atlas_accession} has fewer than 5,000 state-neutral control genes")

        # Score every author marker set. All sets compete for the cell label, including the non-myeloid sets.
        atlas_score_columns: list[str] = []
        for atlas_state, atlas_genes in atlas_markers.items():
            atlas_use = [atlas_gene for atlas_gene in atlas_genes if atlas_gene in atlas_adata.var_names]
            if len(atlas_use) < 30:
                raise RuntimeError(f"{atlas_accession} {atlas_state} has fewer than 30 measurable markers")
            atlas_column = f"state_score_{atlas_state}"

            # Keep 50 controls and 25 bins. Reset the original seed for each marker set to reproduce its score.
            sc.tl.score_genes(
                atlas_adata,
                gene_list    = atlas_use,
                gene_pool    = atlas_neutral_control_pool + atlas_use,
                ctrl_size    = 50,
                n_bins       = 25,
                score_name   = atlas_column,
                random_state = 20260713,
                ctrl_as_ref  = False,
                use_raw      = False,
            )

            atlas_score_columns.append(atlas_column)

        # Each row is a cell and each column is a marker-set score. Use the original argsort for the top two.
        atlas_score_matrix = atlas_adata.obs[atlas_score_columns].to_numpy(float)
        atlas_order = np.argsort(atlas_score_matrix, axis=1)
        atlas_state_names = np.array(list(atlas_markers), dtype=object)
        atlas_winner_index = atlas_order[:, -1]
        atlas_runner_index = atlas_order[:, -2]
        atlas_winners = atlas_state_names[atlas_winner_index]
        atlas_winning_score = atlas_score_matrix[np.arange(atlas_adata.n_obs), atlas_winner_index]
        atlas_runner_score = atlas_score_matrix[np.arange(atlas_adata.n_obs), atlas_runner_index]

        # Keep 50 controls, 25 bins and the original seed for the MG-Inflam20 score.
        sc.tl.score_genes(
            atlas_adata,
            gene_list    = atlas_present_raw20,
            gene_pool    = atlas_neutral_control_pool + atlas_present_raw20,
            ctrl_size    = 50,
            n_bins       = 25,
            score_name   = "raw20_score",
            random_state = 20260713,
            ctrl_as_ref  = False,
            use_raw      = False,
        )

        # Record all QC-passing cells, their winning states, runner-up margins and raw20 scores. Flag myeloid states.
        atlas_raw20_scores = atlas_adata.obs["raw20_score"].to_numpy(float)
        atlas_focus = np.isin(atlas_winners, myeloid_states)

        atlas_cell_frames.append(
            pd.DataFrame(
                {
                    "cell_id": atlas_adata.obs_names.astype(str),
                    "geo_accession": atlas_accession,
                    "sample": atlas_sample,
                    "author_marker_state": atlas_winners,
                    "state_score": atlas_winning_score,
                    "runner_up_state": atlas_state_names[atlas_runner_index],
                    "state_margin": atlas_winning_score - atlas_runner_score,
                    "raw20_score": atlas_raw20_scores,
                    "in_prespecified_myeloid_states": atlas_focus,
                }
            )
        )

        # Extract normalized raw20 expression. Densify the sparse branch while preserving its precision.
        atlas_normalized_raw20 = atlas_adata[:, atlas_present_raw20].X
        if sparse.issparse(atlas_normalized_raw20):
            atlas_normalized_raw20 = atlas_normalized_raw20.toarray()

        atlas_normalized_raw20 = np.asarray(atlas_normalized_raw20, dtype=np.float32)

        # Summarize expression for states with at least 30 cells, regardless of how many cells are in other states.
        for atlas_state in myeloid_states:
            atlas_target = atlas_focus & (atlas_winners == atlas_state)
            atlas_rest = atlas_focus & (atlas_winners != atlas_state)
            atlas_n_target = int(atlas_target.sum())
            atlas_n_rest = int(atlas_rest.sum())
            if atlas_n_target < 30:
                continue

            # Record the within-patient raw20 mean and each gene's mean expression and detection fraction.
            atlas_state_score_frames.append(
                pd.DataFrame(
                    {
                        "geo_accession": [atlas_accession],
                        "sample": [atlas_sample],
                        "state": [atlas_state],
                        "n_cells": [atlas_n_target],
                        "mean_raw20_score": [float(atlas_raw20_scores[atlas_target].mean())],
                    }
                )
            )

            atlas_state_expr = atlas_normalized_raw20[atlas_target]
            for atlas_gene_index, atlas_gene in enumerate(atlas_present_raw20):
                atlas_gene_patient_rows.append(
                    {
                        "geo_accession": atlas_accession,
                        "sample": atlas_sample,
                        "state": atlas_state,
                        "gene": atlas_gene,
                        "n_cells": atlas_n_target,
                        "mean_log_expression": float(atlas_state_expr[:, atlas_gene_index].mean()),
                        "pct_detected": float((atlas_state_expr[:, atlas_gene_index] > 0).mean() * 100),
                    }
                )

            # With 30 or more cells on each side, compute state-versus-rest score and raw-count pseudobulk contrasts.
            if atlas_n_rest < 30:
                continue

            atlas_delta = float(atlas_raw20_scores[atlas_target].mean() - atlas_raw20_scores[atlas_rest].mean())

            # Sum the state's counts across all genes for library size. Average log2(CPM + 1) over raw20.
            atlas_target_counts = atlas_counts[atlas_target]
            atlas_target_library = float(atlas_target_counts.sum())
            atlas_target_pb = np.nan
            if atlas_target_library > 0:
                atlas_target_gene_counts = np.asarray(atlas_target_counts[:, atlas_raw20_indices].sum(axis=0)).ravel()
                atlas_target_log_cpm = np.log2(atlas_target_gene_counts / atlas_target_library * 1_000_000 + 1)
                atlas_target_pb = float(atlas_target_log_cpm.mean())

            # Apply the same formula to other states. Zero-size libraries retain NaN for the alternative score.
            atlas_rest_counts = atlas_counts[atlas_rest]
            atlas_rest_library = float(atlas_rest_counts.sum())
            atlas_rest_pb = np.nan
            if atlas_rest_library > 0:
                atlas_rest_gene_counts = np.asarray(atlas_rest_counts[:, atlas_raw20_indices].sum(axis=0)).ravel()
                atlas_rest_log_cpm = np.log2(atlas_rest_gene_counts / atlas_rest_library * 1_000_000 + 1)
                atlas_rest_pb = float(atlas_rest_log_cpm.mean())
            atlas_effect_rows.append(
                {
                    "geo_accession": atlas_accession,
                    "sample": atlas_sample,
                    "state": atlas_state,
                    "n_target": atlas_n_target,
                    "n_rest": atlas_n_rest,
                    "raw20_delta": atlas_delta,
                    "raw20_pseudobulk_delta": atlas_target_pb - atlas_rest_pb,
                }
            )

        del atlas_adata, atlas_counts, atlas_normalized_raw20

    # Combine the sample-level records into input, quality control and gene coverage tables.
    atlas_locked_table = pd.DataFrame(atlas_locked_rows).sort_values("geo_accession")
    atlas_qc_table = pd.DataFrame(atlas_qc_rows).sort_values("geo_accession")
    atlas_coverage_table = pd.DataFrame(atlas_coverage_rows).sort_values("geo_accession")
    atlas_cells_table_total = pd.concat(atlas_cell_frames, ignore_index=True)
    atlas_state_scores_table_total = pd.concat(atlas_state_score_frames, ignore_index=True)
    atlas_effects_table = pd.DataFrame(atlas_effect_rows)
    atlas_gene_patient_table = pd.DataFrame(atlas_gene_patient_rows)

    # The frozen cohort has 21 samples, 131,780 input cells and 120,766 cells passing QC.
    if atlas_locked_table["geo_accession"].nunique() != 21 or atlas_locked_table["sample"].nunique() != 21:
        raise RuntimeError("The frozen GSE278456 manifest must contain 21 unique GEO accessions and patients")

    if int(atlas_qc_table["n_input"].sum()) != 131780 or int(atlas_qc_table["n_qc_pass"].sum()) != 120766:
        raise RuntimeError("The GSE278456 cell totals differ from the frozen Step 30 totals")
    if atlas_cells_table_total["cell_id"].duplicated().any():
        raise RuntimeError("GSE278456 cell IDs are not unique")

    # Use equal-weight patient contrasts for the original bootstrap interval and exact sign-flip P value.
    atlas_summary_rows = []
    for atlas_index, atlas_state in enumerate(myeloid_states):
        atlas_part_table = atlas_effects_table.loc[atlas_effects_table["state"].eq(atlas_state)].copy()
        atlas_values = atlas_part_table["raw20_delta"].to_numpy(float)

        def bootstrap_mean_ci(values: np.ndarray, seed: int) -> tuple[float, float, float]:
            """Return the mean and 95% interval from 20,000 seeded patient bootstraps. Empty input returns NaN."""

            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]

            if len(values) == 0:
                result = (np.nan, np.nan, np.nan)
                return result
            mean = float(values.mean())
            if len(values) == 1:
                result = (mean, mean, mean)
                return result

            # Use 20,000 seeded draws per state and take the 2.5th and 97.5th percentiles of the means.
            rng = np.random.default_rng(seed)
            draws = rng.choice(values, size=(20000, len(values)), replace=True).mean(axis=1)
            lo, hi = np.quantile(draws, [0.025, 0.975])

            result = (mean, float(lo), float(hi))
            return result

        atlas_mean, atlas_low, atlas_high = bootstrap_mean_ci(atlas_values, 20260713 + atlas_index)
        atlas_pb_values = atlas_part_table["raw20_pseudobulk_delta"].to_numpy(float)

        def exact_two_sided_signflip(values: np.ndarray) -> float:
            """Return exact two-sided sign-flip P with a fixed first sign and 1e-12 tolerance. Use NaN for n < 3."""

            values = np.asarray(values, dtype=float)

            values = values[np.isfinite(values)]
            n = len(values)
            if n < 3:
                result = np.nan
                return result
            observed = abs(values.mean())

            # Flipping all signs preserves the absolute statistic. Fix the first sign to halve enumeration.
            tail = values[1:]
            total = 1 << (n - 1)
            extreme = 0

            # Enumerate blocks of 65,536 sign combinations to limit memory while keeping the exact test.
            chunk_size = 65536
            bit_positions = np.arange(n - 1, dtype=np.uint64)
            for start in range(0, total, chunk_size):
                stop = min(start + chunk_size, total)
                codes = np.arange(start, stop, dtype=np.uint64)[:, None]
                signs = (((codes >> bit_positions) & 1).astype(np.float64) * 2.0) - 1.0
                means = (values[0] + signs @ tail) / n
                extreme += int(np.count_nonzero(np.abs(means) >= observed - 1e-12))
            result = extreme / total
            return result

        # Record each state's mean effect, positive patients, exact P and pseudobulk direction agreement.
        atlas_summary_rows.append(
            {
                "state": atlas_state,
                "n_patients": len(atlas_values),
                "mean_delta": atlas_mean,
                "ci95_low": atlas_low,
                "ci95_high": atlas_high,
                "median_delta": float(np.median(atlas_values)) if len(atlas_values) else np.nan,
                "n_positive": int((atlas_values > 0).sum()),
                "exact_two_sided_signflip_p": exact_two_sided_signflip(atlas_values),
                "mean_pseudobulk_delta": float(atlas_pb_values.mean()) if len(atlas_pb_values) else np.nan,
                "pseudobulk_n_positive": int((atlas_pb_values > 0).sum()),
                "score_pseudobulk_mean_direction_concordant": bool(
                    np.sign(atlas_mean) == np.sign(atlas_pb_values.mean())
                )
                if len(atlas_values)
                else False,
            }
        )

    # Apply BH correction across the nine prespecified state tests.
    atlas_effect_summary_table = pd.DataFrame(atlas_summary_rows)
    atlas_ok = atlas_effect_summary_table["exact_two_sided_signflip_p"].notna()
    atlas_effect_summary_table["fdr_9_states"] = np.nan
    atlas_effect_summary_table.loc[atlas_ok, "fdr_9_states"] = multipletests(
        atlas_effect_summary_table.loc[atlas_ok, "exact_two_sided_signflip_p"], method="fdr_bh"
    )[1]

    # Average expression within patients first, then compare state means within each gene with equal patient weights.
    atlas_gene_summary_table = (
        atlas_gene_patient_table.groupby(["state", "gene"], observed=True)
        .agg(
            n_patients                           = ("sample", "nunique"),
            patient_balanced_mean_log_expression = ("mean_log_expression", "mean"),
            patient_balanced_pct_detected        = ("pct_detected", "mean"),
        )
        .reset_index()
    )
    atlas_gene_summary_table["within_gene_z"] = np.nan

    # Standardize each gene across the nine state means with the population SD. Use zero for constant means.
    for _, group_table in atlas_gene_summary_table.groupby("gene", observed=True):
        group_values = group_table["patient_balanced_mean_log_expression"]
        group_sd = group_values.std(ddof=0)
        if group_sd > 0:
            group_z = (group_values - group_values.mean()) / group_sd
        else:
            group_z = 0
        atlas_gene_summary_table.loc[group_table.index, "within_gene_z"] = group_z

    # Keep the nine existing output suffixes. The cell annotation table uses gzip.
    results = {
        "locked_manifest": atlas_locked_table,
        "qc": atlas_qc_table,
        "coverage": atlas_coverage_table,
        "cells": atlas_cells_table_total,
        "patient_state_scores": atlas_state_scores_table_total,
        "patient_state_effects": atlas_effects_table,
        "state_effect_summary": atlas_effect_summary_table,
        "gene_patient": atlas_gene_patient_table,
        "gene_summary": atlas_gene_summary_table,
    }

    # Write the original result tables and check that each newly written file contains data.
    output_path = write_root / "GSE278456_author_marker_audit.csv"
    marker_audit_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
    for name, data_table in results.items():
        if name == "cells":

            # Write and check the compressed cell annotations.
            output_path = write_root / "GSE278456_cell_annotations_true_raw20.csv.gz"
            data_table.to_csv(output_path, index=False, compression={"method": "gzip", "mtime": 0})
            if not output_path.is_file() or output_path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty output: {output_path.name}")
        else:

            # Write and check each uncompressed CSV.
            output_path = write_root / f"GSE278456_{name}.csv"
            data_table.to_csv(output_path, index=False)
            if not output_path.is_file() or output_path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty output: {output_path.name}")
    print("State summaries:", len(atlas_effect_summary_table))
