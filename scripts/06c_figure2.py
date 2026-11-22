#!/usr/bin/env python3
# Purpose: draw the two Figure 2 panels from the frozen raw20 effects and recomputed author-marker scores.
# Inputs: Step 48 effects, gene summaries, marker audit, cell annotations and the original GSE278456 H5 files.
# Outputs: Source tables and panel manifest in the Step 50 write directory, with PDF and PNG panels in figures.
# Run: python scripts/06c_figure2.py from the repository root in the single-cell environment.
# The released workflow uses Scanpy and AnnData for the marker scores and Matplotlib for the panels.
# Prepared input paths, output names and fields follow the original public analysis.

from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import scanpy as sc
import seaborn as sns
import os
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

# The --plot-only option reads the frozen marker tables. The default branch replays the Step 48 assignments.
arguments_parser = argparse.ArgumentParser(description="Build the two Figure 2 panels")
arguments_parser.add_argument(
    "--plot-only",
    action = "store_true",
    help   = "Redraw the panels from the frozen source tables without reading the H5 files.",
)
args = arguments_parser.parse_args()

# Use the original font, stroke widths and white background for the PDF.
mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.2,
        "axes.labelsize": 7.6,
        "axes.titlesize": 8.2,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)

# Keep the original Scanpy settings.
sc.settings.verbosity = 0

# Create the source and figure directories at their original write and figures paths.
root = Path(__file__).resolve().parents[1]

write_root = root / "write/50_figure2_two_panel_final/Figure2"

source_out = write_root / "source_data"

source_out.mkdir(parents=True, exist_ok=True)

fig_root = root / "figures/50_figure2_two_panel_final/Figure2"

fig_root.mkdir(parents=True, exist_ok=True)

write_root.mkdir(parents=True, exist_ok=True)

step48 = root / "write/48_figure2_identity_state_rebuild/Figure2"

# Read the Step 48 patient contrasts, state and gene summaries, and author-marker audit.
effect_file = step48 / "GSE278456_patient_state_effects.csv"

inputs_effects_table = pd.read_csv(effect_file)

summary_file = step48 / "GSE278456_state_effect_summary.csv"

inputs_summary_table = pd.read_csv(summary_file)

gene_summary_file = step48 / "GSE278456_gene_summary.csv"

# Retain the original table columns, missing-value handling and row order.
inputs_gene_summary_table = pd.read_csv(gene_summary_file)

marker_audit_file = step48 / "GSE278456_author_marker_audit.csv"

inputs_marker_audit_table = pd.read_csv(marker_audit_file)

# The frozen inputs contain 108 patient-state effects, nine states and 171 gene-state combinations.
if inputs_effects_table.shape[0] != 108 or inputs_effects_table["state"].nunique() != 9:
    raise RuntimeError("The Step 48 patient-state effects table must contain 108 rows across nine states")

# Keep the nine-state display order. The workbook determines the order of competing marker sets.
state_order = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]

if set(inputs_summary_table["state"]) != set(state_order):
    raise RuntimeError("The Step 48 state list differs from the prespecified list")
if inputs_gene_summary_table.shape[0] != 171:
    raise RuntimeError("The Step 48 raw20 gene-by-state table must contain 171 rows")

# Keep the original 20 genes in Miller's published rank order. The symbols follow the source data.
raw20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]

# AC253572.2 is unmeasurable. Check the other 19 program genes across the nine prespecified states.
measurable_raw20 = [gene for gene in raw20 if gene != "AC253572.2"]

if set(inputs_gene_summary_table["gene"]) != set(measurable_raw20):
    raise RuntimeError("The Step 48 measurable raw20 genes differ from the frozen set")
if inputs_gene_summary_table.groupby("gene")["state"].nunique().ne(9).any():
    raise RuntimeError("Each measurable raw20 gene must cover all nine states")

# Collect the state-assignment markers and check that they exclude all measurable raw20 genes.
inputs_focus_audit_table = inputs_marker_audit_table.loc[inputs_marker_audit_table["state"].isin(state_order)].copy()
if inputs_focus_audit_table["state"].nunique() != 9:
    raise RuntimeError("The author-marker audit lacks one or more prespecified states")

inputs_used = set()
for inputs_value in inputs_focus_audit_table["markers_used"].fillna(""):
    inputs_used.update(inputs_gene for inputs_gene in str(inputs_value).split(";") if inputs_gene)
inputs_overlap = sorted(inputs_used & set(measurable_raw20))
if inputs_overlap:
    raise RuntimeError(f"State-assignment markers overlap the raw20 genes: {inputs_overlap}")
(effects_table, state_summary_table, gene_summary_table, _) = (
    inputs_effects_table,
    inputs_summary_table,
    inputs_gene_summary_table,
    inputs_marker_audit_table,
)

# Add state statistics and the pseudobulk direction flag to the patient contrast source table.
waterfall_effects_table = effects_table
waterfall_summary_table = state_summary_table

waterfall_summary_index_table = waterfall_summary_table.set_index("state")

# Rank patients within each prespecified state by their original contrast and retain the inference fields.
waterfall_rows: list[pd.DataFrame] = []
for waterfall_state in state_order:
    waterfall_data_table = (
        waterfall_effects_table.loc[waterfall_effects_table["state"].eq(waterfall_state)]
        .sort_values("raw20_delta")
        .copy()
    )
    waterfall_data_table["within_state_rank"] = np.arange(1, len(waterfall_data_table) + 1)
    waterfall_data_table["direction"] = np.where(waterfall_data_table["raw20_delta"] > 0, "Above rest", "Below rest")
    waterfall_data_table["n_patients"] = int(waterfall_summary_index_table.loc[waterfall_state, "n_patients"])

    waterfall_data_table["n_positive"] = int(waterfall_summary_index_table.loc[waterfall_state, "n_positive"])
    waterfall_data_table["mean_delta"] = float(waterfall_summary_index_table.loc[waterfall_state, "mean_delta"])
    waterfall_data_table["fdr_9_states"] = float(waterfall_summary_index_table.loc[waterfall_state, "fdr_9_states"])
    waterfall_data_table["score_pseudobulk_mean_direction_concordant"] = bool(
        waterfall_summary_index_table.loc[waterfall_state, "score_pseudobulk_mean_direction_concordant"]
    )
    waterfall_rows.append(waterfall_data_table)
waterfall_table_total = pd.concat(waterfall_rows, ignore_index=True)

# Read or compute the within-patient marker means and the 9 × 9 summary.
marker_patient_file = source_out / "GSE278456_patient_marker_set_scores_source.csv"
marker_summary_file = source_out / "GSE278456_marker_set_assignment_qc_source.csv"

# Plot-only mode reads the two marker-score tables. The default computes scores from the workbook and H5 files.
if args.plot_only:
    required = [marker_patient_file, marker_summary_file]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError(f"--plot-only requires these missing frozen source tables: {missing}")
    marker_patient_table = pd.read_csv(marker_patient_file)
    marker_summary_table = pd.read_csv(marker_summary_file)
else:

    # Import the Step 48 marker reader. Its main reconstruction is guarded by the script entry point.
    step48_script = root / "scripts/06b_cell_identity_and_state.py"
    loader_spec = importlib.util.spec_from_file_location("step48_fig2", step48_script)
    if loader_spec is None or loader_spec.loader is None:
        raise RuntimeError("Cannot load the Step 48 marker reader")
    loader_module = importlib.util.module_from_spec(loader_spec)
    loader_spec.loader.exec_module(loader_module)
    marker_step48 = loader_module
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    marker_book = external_root / "write/30_dataset_rescue_search/source_metadata/NIHMS2115262-supplement-S2_S5_S6.xls"
    marker_markers, _ = marker_step48.read_author_markers(marker_book, raw20, state_order)
    if len(marker_markers) != 17:
        raise RuntimeError(f"Expected 17 state-marker sets in the workbook, found {len(marker_markers)}")

    # Read the sample, QC and cell identifiers with their existing columns and missing-value handling.
    # Use the frozen Step 30 sample and QC tables and sort the 21 H5 files by filename.
    marker_manifest_table = pd.read_csv((root / "write/30_dataset_rescue_search") / "GSE278456_sample_manifest.csv")
    marker_manifest_table["geo_accession"] = marker_manifest_table["geo_accession"].astype(str)
    marker_existing_qc_table = pd.read_csv(
        (root / "write/30_dataset_rescue_search") / "GSE278456_single_cell_qc.csv"
    ).set_index("geo_accession")
    h5_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE278456_myeloid_h5"
    marker_paths = sorted(h5_root.glob("*.h5"))
    if len(marker_paths) != 21:
        raise RuntimeError(f"Expected 21 GSE278456 H5 files, found {len(marker_paths)}")

    # Match the frozen Step 48 assignments by cell_id so row order does not determine alignment.
    cell_file = step48 / "GSE278456_cell_annotations_true_raw20.csv.gz"
    marker_saved_cells_table = pd.read_csv(
        cell_file, usecols=["cell_id", "geo_accession", "author_marker_state", "state_score"]
    )
    if marker_saved_cells_table.shape[0] != 120766 or marker_saved_cells_table["cell_id"].duplicated().any():
        raise RuntimeError("The frozen Step 48 cell annotation table must contain 120,766 unique cells")
    marker_saved_cells_table = marker_saved_cells_table.set_index("cell_id")

    marker_all_marker_genes = set().union(*marker_markers.values())
    marker_patient_rows: list[dict[str, object]] = []
    marker_audit_rows: list[dict[str, object]] = []

    # Check IDH-wild-type eligibility and QC cell counts while preserving the order of competing marker sets.
    for marker_path in marker_paths:
        marker_accession = marker_path.name.split("_")[0]
        marker_manifest_match_table = marker_manifest_table.loc[
            marker_manifest_table["geo_accession"].eq(marker_accession)
        ]
        if len(marker_manifest_match_table) != 1:
            raise RuntimeError(f"{marker_accession} does not match exactly one manifest row")
        marker_manifest_row = marker_manifest_match_table.iloc[0]
        if not bool(marker_manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"Non-IDH-wild-type GBM in the frozen H5 set: {marker_accession}")
        marker_sample = str(marker_manifest_row["sample_title"])

        # Compute detected genes, total counts and mitochondrial percentage from the raw 10x counts.
        marker_adata = sc.read_10x_h5(marker_path)
        marker_adata.var_names_make_unique()
        marker_n_input = marker_adata.n_obs
        marker_adata.var["mt"] = marker_adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(marker_adata, qc_vars=["mt"], inplace=True, percent_top=None)

        # Apply the Step 48 QC bounds to reproduce the input cells for state assignment.
        marker_keep = (
            marker_adata.obs["n_genes_by_counts"].between(200, 8000)
            & marker_adata.obs["total_counts"].ge(500)
            & marker_adata.obs["pct_counts_mt"].le(20)
        )
        marker_adata = marker_adata[marker_keep].copy()
        marker_n_qc = marker_adata.n_obs
        marker_old = marker_existing_qc_table.loc[marker_accession]
        if int(marker_old["n_input"]) != marker_n_input or int(marker_old["n_qc_pass"]) != marker_n_qc:
            raise RuntimeError(f"{marker_accession} cell counts differ from the frozen QC counts")

        # Check that the measurable raw20 genes and their order match the frozen set.
        marker_present_raw20 = [marker_gene for marker_gene in raw20 if marker_gene in marker_adata.var_names]
        if marker_present_raw20 != measurable_raw20:
            raise RuntimeError(f"{marker_accession} raw20 coverage or gene order differs from the frozen set")
        marker_adata.obs_names = [
            f"{marker_accession}:{marker_barcode}" for marker_barcode in marker_adata.obs_names.astype(str)
        ]

        # Normalize to 10,000 counts and apply log1p. Keep the source counts intact.
        sc.pp.normalize_total(marker_adata, target_sum=10000)
        sc.pp.log1p(marker_adata)

        # Exclude raw20 and all author markers from the control pool. At least 5,000 neutral genes must remain.
        marker_excluded_controls = set(raw20) | marker_all_marker_genes
        marker_neutral_control_pool = [
            marker_gene for marker_gene in marker_adata.var_names if marker_gene not in marker_excluded_controls
        ]
        if len(marker_neutral_control_pool) < 5000:
            raise RuntimeError(f"{marker_accession} has fewer than 5,000 state-neutral control genes")

        # Score all 17 marker sets. Each set must have at least 30 measurable markers to compete.
        marker_score_columns: list[str] = []
        for marker_state, marker_genes in marker_markers.items():
            marker_use = [marker_gene for marker_gene in marker_genes if marker_gene in marker_adata.var_names]
            if len(marker_use) < 30:
                raise RuntimeError(f"{marker_accession} {marker_state} has fewer than 30 measurable markers")
            marker_column = f"state_score_{marker_state}"

            # Keep 50 controls and 25 bins. Reset the original seed for each marker set to reproduce its score.
            sc.tl.score_genes(
                marker_adata,
                gene_list    = marker_use,
                gene_pool    = marker_neutral_control_pool + marker_use,
                ctrl_size    = 50,
                n_bins       = 25,
                score_name   = marker_column,
                random_state = 20260713,
                ctrl_as_ref  = False,
                use_raw      = False,
            )

            marker_score_columns.append(marker_column)

        # Assign cells by the original argsort over all worksheet marker sets, then select the myeloid states.
        marker_score_matrix = marker_adata.obs[marker_score_columns].to_numpy(float)
        marker_order = np.argsort(marker_score_matrix, axis=1)
        marker_state_names = np.array(list(marker_markers), dtype=object)
        marker_winner_index = marker_order[:, -1]
        marker_winners = marker_state_names[marker_winner_index]
        marker_winning_score = marker_score_matrix[np.arange(marker_adata.n_obs), marker_winner_index]

        # Compare each cell's winning label and score with the frozen calls using the original 1e-10 tolerance.
        marker_frozen_table = marker_saved_cells_table.loc[marker_adata.obs_names]
        marker_winner_matches = marker_frozen_table["author_marker_state"].to_numpy(object) == marker_winners
        marker_score_difference = np.abs(marker_frozen_table["state_score"].to_numpy(float) - marker_winning_score)
        if not bool(marker_winner_matches.all()) or float(marker_score_difference.max()) > 1e-10:
            raise RuntimeError(f"{marker_accession} differs from the frozen Step 48 state assignments")
        marker_audit_rows.append(
            {
                "geo_accession": marker_accession,
                "sample": marker_sample,
                "n_input": marker_n_input,
                "n_qc_pass": marker_n_qc,
                "n_winner_exact_match": int(marker_winner_matches.sum()),
                "winner_match_fraction": float(marker_winner_matches.mean()),
                "max_abs_winner_score_difference": float(marker_score_difference.max()),
                "n_competing_marker_sets": len(marker_markers),
            }
        )

        # For each of the nine states, average all nine marker scores within patient when at least 30 cells qualify.
        marker_focus = np.isin(marker_winners, state_order)
        marker_focus_scores = {
            marker_state: marker_adata.obs[f"state_score_{marker_state}"].to_numpy(float)
            for marker_state in state_order
        }
        for marker_assigned_state in state_order:
            marker_target = marker_focus & (marker_winners == marker_assigned_state)
            marker_n_cells = int(marker_target.sum())
            if marker_n_cells < 30:
                continue

            for marker_tested_marker_set in state_order:
                marker_patient_rows.append(
                    {
                        "geo_accession": marker_accession,
                        "sample": marker_sample,
                        "assigned_state": marker_assigned_state,
                        "tested_marker_set": marker_tested_marker_set,
                        "n_cells": marker_n_cells,
                        "mean_marker_set_score": float(
                            marker_focus_scores[marker_tested_marker_set][marker_target].mean()
                        ),
                    }
                )

        del marker_adata, marker_score_matrix

    # Each row represents one patient, assigned state and tested marker set. Patients receive equal weight.
    marker_patient_scores_table = pd.DataFrame(marker_patient_rows)
    marker_audit_table = pd.DataFrame(marker_audit_rows).sort_values("geo_accession")
    if (
        int(marker_audit_table["n_qc_pass"].sum()) != 120766
        or not marker_audit_table["winner_match_fraction"].eq(1).all()
    ):
        raise RuntimeError("The 21 H5 inputs did not reproduce all frozen Step 48 state assignments")
    if marker_patient_scores_table.shape[0] != 972:
        raise RuntimeError(f"Expected 972 patient-marker rows, found {marker_patient_scores_table.shape[0]}")

    # Average patient means within each assigned state, then standardize across states for each marker set.
    marker_summary_table = (
        marker_patient_scores_table.groupby(["tested_marker_set", "assigned_state"], observed=True)
        .agg(n_patients=("sample", "nunique"), patient_balanced_mean_marker_set_score=("mean_marker_set_score", "mean"))
        .reset_index()
    )
    marker_summary_table["within_marker_set_z"] = np.nan

    # Use the population standard deviation for Z-scores. Assign zero for constant scores and keep row order.
    for _, group_table in marker_summary_table.groupby("tested_marker_set", observed=True):
        group_values = group_table["patient_balanced_mean_marker_set_score"]
        group_sd = group_values.std(ddof=0)
        if group_sd > 0:
            group_z = (group_values - group_values.mean()) / group_sd
        else:
            group_z = 0
        marker_summary_table.loc[group_table.index, "within_marker_set_z"] = group_z
    if marker_summary_table.shape[0] != 81:
        raise RuntimeError("The marker-set QC summary must cover a 9 × 9 state grid")
    (marker_patient_table, marker_summary_table, audit_table) = (
        marker_patient_scores_table,
        marker_summary_table,
        marker_audit_table,
    )

    # Write the source table and confirm the current output file is nonempty.
    output_path = marker_patient_file
    marker_patient_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

    # Write the source table and confirm the current output file is nonempty.
    output_path = marker_summary_file
    marker_summary_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

    # Save the sample-level assignment checks when recalculating the marker scores.
    output_path = source_out / "GSE278456_step48_assignment_reproduction_source.csv"
    audit_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

# Require 972 patient-state-marker rows and 81 summary rows in either execution mode.
if marker_patient_table.shape[0] != 972 or marker_summary_table.shape[0] != 81:
    raise RuntimeError("Unexpected row counts in the frozen marker-set source tables")
waterfall_plot_source_table = waterfall_table_total

waterfall_plot_stem = "Figure2_state_patient_waterfall"

waterfall_plot_source = source_out / f"{waterfall_plot_stem}_source.csv"

# Write the source table and confirm the current output file is nonempty.
output_path = waterfall_plot_source
waterfall_plot_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
waterfall_plot_summary_table = waterfall_plot_source_table.drop_duplicates("state").set_index("state")

mm = 1 / 25.4

# Draw the original 3 × 3 waterfall with shared scales and patient slots.
waterfall_plot_fig, waterfall_plot_axes = plt.subplots(3, 3, figsize=(183 * mm, 104 * mm), sharex=True, sharey=True)
waterfall_plot_y_min, waterfall_plot_y_max = -0.76, 0.56

# Use red and blue for positive and negative contrasts and the fixed state colors for titles and means.
blue = "#3C5488"

negative = blue

red = "#E64B35"

positive = red

salmon = "#F39B7F"

state_colors = {
    "MCG1": "#00A087",
    "MCG2": "#91D1C2",
    "MCG3": "#4DBBD5",
    "MCG4": "#8491B4",
    "MCG5": blue,
    "MAC1": salmon,
    "MAC2": red,
    "E-MDSC": "#B09C85",
    "M-MDSC": "#7E6148",
}

# Rank patients by the contrast within each state. Shading distinguishes the two sides of zero.
for waterfall_plot_index, (waterfall_plot_ax, waterfall_plot_state) in enumerate(
    zip(waterfall_plot_axes.flat, state_order, strict=True)
):
    waterfall_plot_data_table = waterfall_plot_source_table.loc[
        waterfall_plot_source_table["state"].eq(waterfall_plot_state)
    ]
    waterfall_plot_x = waterfall_plot_data_table["within_state_rank"].to_numpy(int)
    waterfall_plot_y = waterfall_plot_data_table["raw20_delta"].to_numpy(float)
    waterfall_plot_colors = np.where(waterfall_plot_y > 0, positive, negative)
    waterfall_plot_ax.axhspan(0, waterfall_plot_y_max, color="#FFF7F4", zorder=0)

    waterfall_plot_ax.axhspan(waterfall_plot_y_min, 0, color="#F4F6FA", zorder=0)
    black = "#222222"
    waterfall_plot_ax.axhline(0, color=black, linewidth=0.65, zorder=2)
    waterfall_plot_ax.bar(
        waterfall_plot_x,
        waterfall_plot_y,
        width     = 0.76,
        color     = waterfall_plot_colors,
        edgecolor = "none",
        alpha     = 0.88,
        zorder    = 3,
    )

    # Draw the patient-equal mean as a dashed line and flag states with fewer than six patients with a dagger.
    waterfall_plot_mean = float(waterfall_plot_summary_table.loc[waterfall_plot_state, "mean_delta"])
    waterfall_plot_ax.axhline(
        waterfall_plot_mean, color=state_colors[waterfall_plot_state], linewidth=1.05, linestyle=(0, (3, 2)), zorder=4
    )
    waterfall_plot_ax.set_xlim(0.25, 21.75)
    waterfall_plot_ax.set_ylim(waterfall_plot_y_min, waterfall_plot_y_max)
    waterfall_patient_count = int(waterfall_plot_summary_table.loc[waterfall_plot_state, "n_patients"])
    waterfall_small_sample_mark = "†" if waterfall_patient_count < 6 else ""
    waterfall_plot_state_title = f"{waterfall_plot_state}{waterfall_small_sample_mark}"
    waterfall_plot_ax.set_title(
        waterfall_plot_state_title, loc="left", color=state_colors[waterfall_plot_state], fontweight="bold", pad=7
    )

    waterfall_positive_count = int(waterfall_plot_summary_table.loc[waterfall_plot_state, "n_positive"])

    # Format the displayed FDR at the three original precisions. Inference uses the unrounded values.
    waterfall_fdr = float(waterfall_plot_summary_table.loc[waterfall_plot_state, "fdr_9_states"])
    if waterfall_fdr < 0.001:
        waterfall_fdr_label = f"{waterfall_fdr:.2e}"
    elif waterfall_fdr < 0.01:
        waterfall_fdr_label = f"{waterfall_fdr:.4f}"
    else:
        waterfall_fdr_label = f"{waterfall_fdr:.3f}"
    waterfall_plot_ax.text(
        1.0,
        1.025,
        (f"{waterfall_positive_count}/{waterfall_patient_count}  |  FDR {waterfall_fdr_label}"),
        transform = waterfall_plot_ax.transAxes,
        ha        = "right",
        va        = "bottom",
        fontsize  = 5.8,
        color     = black,
    )

    # Use dashed borders for MCG3 to mark the different directions of the score and pseudobulk contrasts.
    grey = "#868686"
    if waterfall_plot_state == "MCG3":
        for waterfall_plot_spine in waterfall_plot_ax.spines.values():
            waterfall_plot_spine.set_color(grey)
            waterfall_plot_spine.set_linestyle((0, (3, 2)))
        waterfall_plot_ax.text(
            0.98,
            0.06,
            "score / pseudobulk\ndirection differs",
            transform = waterfall_plot_ax.transAxes,
            ha        = "right",
            va        = "bottom",
            fontsize  = 5.3,
            color     = grey,
        )
    else:
        waterfall_plot_ax.spines[["top", "right"]].set_visible(False)

    # Show y-axis labels in the left column and patient ranks in the bottom row. Keep the same 21 slots.
    waterfall_plot_ax.spines[["left", "bottom"]].set_color(black)
    waterfall_plot_ax.tick_params(length=2.2, width=0.65)
    if waterfall_plot_index % 3 != 0:
        waterfall_plot_ax.tick_params(axis="y", left=False, labelleft=False)
    if waterfall_plot_index % 3 == 0:
        waterfall_plot_ax.set_ylabel("State − other myeloid\nMiller-IM score")
    if waterfall_plot_index >= 6:
        waterfall_plot_ax.set_xticks([1, 7, 14, 21])
        waterfall_plot_ax.set_xlabel("Patients ranked within state")
    else:
        waterfall_plot_ax.tick_params(axis="x", bottom=False, labelbottom=False)
waterfall_plot_fig.text(
    0.995,
    0.012,
    "† fewer than 6 patients; all panels use 21 fixed patient slots",
    ha       = "right",
    va       = "bottom",
    fontsize = 5.6,
    color    = grey,
)

# Use the original subplot spacing and tight PDF crop with 0.035 inches of padding.
waterfall_plot_fig.subplots_adjust(left=0.09, right=0.995, bottom=0.13, top=0.965, wspace=0.13, hspace=0.29)
waterfall_plot_pdf = fig_root / f"{waterfall_plot_stem}.pdf"
waterfall_plot_png = fig_root / f"{waterfall_plot_stem}.png"

# Draw the figure and save the PDF with the original crop, padding and metadata.
waterfall_plot_fig.canvas.draw()
waterfall_plot_fig.savefig(
    waterfall_plot_pdf,
    bbox_inches = "tight",
    pad_inches  = 0.035,
    metadata    = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None},
)

# Render the PDF to a 600 dpi PNG and verify the output through the shared function.
render_pdf_png(waterfall_plot_pdf, waterfall_plot_png, dpi=600)
plt.close(waterfall_plot_fig)
waterfall_record = {
    "panel": "2A",
    "source": str(waterfall_plot_source),
    "pdf": str(waterfall_plot_pdf),
    "png": str(waterfall_plot_png),
}

# Write separate source tables with the original held-out gene and state-assignment QC block names.
heatmap_gene_summary_table = gene_summary_table
heatmap_marker_summary_table = marker_summary_table

heatmap_state_summary_table = state_summary_table

heatmap_stem = "Figure2_raw20_marker_double_heatmap"
heatmap_gene_source_table = heatmap_gene_summary_table.copy()
heatmap_gene_source_table["block"] = "Held-out Miller-IM genes; not used for state assignment"
heatmap_marker_source_table = heatmap_marker_summary_table.copy()
heatmap_marker_source_table["block"] = "Author marker-set assignment scores; classification QC"
heatmap_source_top = source_out / f"{heatmap_stem}_raw20_source.csv"

# Write the source table and confirm the current output file is nonempty.
output_path = heatmap_source_top
heatmap_gene_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
heatmap_source_bottom = source_out / f"{heatmap_stem}_marker_qc_source.csv"

# Write the source table and confirm the current output file is nonempty.
output_path = heatmap_source_bottom
heatmap_marker_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# Keep the held-out gene expression above the marker-set quality control scores.
heatmap_top_table = heatmap_gene_summary_table.pivot(index="gene", columns="state", values="within_gene_z").reindex(
    index= measurable_raw20, columns=state_order
)
heatmap_bottom_table = heatmap_marker_summary_table.pivot(
    index= "tested_marker_set", columns="assigned_state", values="within_marker_set_z"
).reindex(index=state_order, columns=state_order)
if heatmap_top_table.isna().any().any() or heatmap_bottom_table.isna().any().any():
    raise RuntimeError("Missing cells in the two-block heatmap")

heatmap_n_by_state = heatmap_state_summary_table.set_index("state")["n_patients"].astype(int).to_dict()

# Label each state with its evaluable patient count. A dagger marks fewer than six patients.
heatmap_state_labels = [
    f"{heatmap_state}\n(n={heatmap_n_by_state[heatmap_state]}){'†' if heatmap_n_by_state[heatmap_state] < 6 else ''}"
    for heatmap_state in state_order
]

# Use four layout rows for the title, 19-gene heatmap, separator and nine-marker-set heatmap.
heatmap_fig = plt.figure(figsize=(183 * mm, 148 * mm))
heatmap_grid = heatmap_fig.add_gridspec(
    nrows         = 4,
    ncols         = 2,
    width_ratios  = [1, 0.034],
    height_ratios = [2.2, 19, 1.35, 5.3],
    left          = 0.16,
    right         = 0.94,
    bottom        = 0.10,
    top           = 0.98,
    wspace        = 0.16,
    hspace        = 0.05,
)

heatmap_title_top = heatmap_fig.add_subplot(heatmap_grid[0, 0])
heatmap_ax_top = heatmap_fig.add_subplot(heatmap_grid[1, 0])
heatmap_separator = heatmap_fig.add_subplot(heatmap_grid[2, 0])
heatmap_ax_bottom = heatmap_fig.add_subplot(heatmap_grid[3, 0])
heatmap_cax_top = heatmap_fig.add_subplot(heatmap_grid[1, 1])
heatmap_cax_bottom = heatmap_fig.add_subplot(heatmap_grid[3, 1])

heatmap_title_top.axis("off")
heatmap_title_top.text(
    0,
    0.93,
    "Held-out Miller-IM genes — not used for state assignment",
    ha         = "left",
    va         = "center",
    fontsize   = 7.4,
    fontweight = "bold",
    color      = black,
)

# Use the blue, white and salmon palette with fixed ±2.5 limits for within_gene_z.
raw20_cmap = LinearSegmentedColormap.from_list("raw20_npg", [blue, "#F7F7F7", salmon])

sns.heatmap(
    heatmap_top_table,
    ax          = heatmap_ax_top,
    cmap        = raw20_cmap,
    vmin        = -2.5,
    vmax        = 2.5,
    center      = 0,
    linewidths  = 0.45,
    linecolor   = "white",
    cbar        = True,
    cbar_ax     = heatmap_cax_top,
    xticklabels = heatmap_state_labels,
    yticklabels = measurable_raw20,
)

# Place state labels above and gene labels on the left. The separate colorbar shows within-gene Z-scores.
heatmap_ax_top.xaxis.tick_top()
heatmap_ax_top.tick_params(axis="x", rotation=0, length=0, pad=2)
heatmap_ax_top.tick_params(axis="y", rotation=0, length=0, pad=2)
heatmap_ax_top.set_xlabel("")
heatmap_ax_top.set_ylabel("")
heatmap_cax_top.set_ylabel("Within-gene z score", fontsize=6.3, labelpad=4)
heatmap_cax_top.tick_params(labelsize=5.7, length=2)

# Identify the lower block as assignment QC and display the previously verified zero gene overlap.
heatmap_separator.axis("off")
heatmap_separator.axhline(0.92, color=black, linewidth=0.8)

heatmap_separator.text(
    0,
    0.55,
    "Author marker-set assignment scores — used to define labels; classification QC",
    ha         = "left",
    va         = "center",
    fontsize   = 6.8,
    fontweight = "bold",
    color      = grey,
)

heatmap_separator.text(
    0, 0.08, "Direct marker–Miller-IM gene overlap after exclusion: 0", ha="left", va="center", fontsize=5.9, color=grey
)

# Use the original gray palette and the same ±2.5 limits for within-marker-set Z-scores.
qc_cmap = LinearSegmentedColormap.from_list("marker_qc", ["#ECECEC", "#A9B3BE", "#4B5563"])

sns.heatmap(
    heatmap_bottom_table,
    ax          = heatmap_ax_bottom,
    cmap        = qc_cmap,
    vmin        = -2.5,
    vmax        = 2.5,
    center      = 0,
    linewidths  = 0.45,
    linecolor   = "white",
    cbar        = True,
    cbar_ax     = heatmap_cax_bottom,
    xticklabels = False,
    yticklabels = [f"{heatmap_state} markers" for heatmap_state in state_order],
)

# Keep the marker labels and gray color scale. The footer states the purpose of the assignment scores.
heatmap_ax_bottom.tick_params(axis="y", rotation=0, length=0, pad=2, colors=grey)
heatmap_ax_bottom.tick_params(axis="x", bottom=False)
heatmap_ax_bottom.set_xlabel("")
heatmap_ax_bottom.set_ylabel("")
heatmap_cax_bottom.set_ylabel("Within-marker-set z score", fontsize=5.8, labelpad=4, color=grey)
heatmap_cax_bottom.tick_params(labelsize=5.5, length=2, colors=grey)
for heatmap_spine in heatmap_cax_bottom.spines.values():
    heatmap_spine.set_edgecolor(grey)

heatmap_fig.text(
    0.16,
    0.035,
    (
        "Columns are fixed a priori; † fewer than 6 patients. The lower block documents l"
        "abel construction and is not independent validation."
    ),
    ha       = "left",
    va       = "bottom",
    fontsize = 5.7,
    color    = grey,
)

heatmap_pdf = fig_root / f"{heatmap_stem}.pdf"
heatmap_png = fig_root / f"{heatmap_stem}.png"

# Draw the figure and save the PDF with the original crop, padding and metadata.
heatmap_fig.canvas.draw()
heatmap_fig.savefig(
    heatmap_pdf,
    bbox_inches = "tight",
    pad_inches  = 0.035,
    metadata    = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None},
)

# Render the PDF to a 600 dpi PNG and verify the output through the shared function.
render_pdf_png(heatmap_pdf, heatmap_png, dpi=600)
plt.close(heatmap_fig)
heatmap_record = {
    "panel": "2B",
    "source": f"{heatmap_source_top};{heatmap_source_bottom}",
    "pdf": str(heatmap_pdf),
    "png": str(heatmap_png),
}
records = [waterfall_record, heatmap_record]
preview_records = records

# Use Pillow to assemble the two PDF-derived panel images into the review preview.
preview_images = [Image.open(preview_record["png"]).convert("RGB") for preview_record in preview_records]
preview_target_width = 1800
preview_resized: list[Image.Image] = []
for preview_image in preview_images:
    preview_height = round(preview_image.height * preview_target_width / preview_image.width)
    preview_resized.append(preview_image.resize((preview_target_width, preview_height), Image.Resampling.LANCZOS))

preview_margin = 85
preview_gap = 40
preview_canvas = Image.new(
    "RGB",
    (
        preview_target_width + preview_margin * 2,
        sum(preview_image.height for preview_image in preview_resized) + preview_gap + preview_margin * 2,
    ),
    "white",
)

preview_draw = ImageDraw.Draw(preview_canvas)

# Keep the original preview font preference for Helvetica, then Arial.
preview_font_candidates = ["/System/Library/Fonts/Helvetica.ttc", "/Library/Fonts/Arial.ttf"]
preview_font_path = next(
    (Path(preview_path) for preview_path in preview_font_candidates if Path(preview_path).exists()), None
)
preview_font = ImageFont.truetype(str(preview_font_path), 54) if preview_font_path else ImageFont.load_default()
preview_y = preview_margin
for preview_label, preview_image in zip(["A", "B"], preview_resized, strict=True):
    preview_draw.text((22, preview_y + 2), preview_label, fill=black, font=preview_font)
    preview_canvas.paste(preview_image, (preview_margin, preview_y))
    preview_y += preview_image.height + preview_gap

# Save the color and grayscale previews from the two PDF-derived panels.
preview_preview = fig_root / "Figure2_two_panel_preview.png"
preview_canvas.save(preview_preview, dpi=(150, 150))
preview_grayscale = fig_root / "Figure2_two_panel_preview_grayscale.png"
preview_canvas.convert("L").convert("RGB").save(preview_grayscale, dpi=(150, 150))
(preview, grayscale) = (preview_preview, preview_grayscale)

# Keep the panel, source, PDF and PNG columns used by the figure assembler.
manifest = write_root / "Figure2_panel_manifest.csv"
panel_manifest_table = pd.DataFrame(records)

# Write the source table and confirm the current output file is nonempty.
output_path = manifest
panel_manifest_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
print("Figure2 panels:", len(records))

# Open each preview and check its dimensions.
for output_path in (preview, grayscale):
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
for preview_path in [preview, grayscale]:
    with Image.open(preview_path) as preview_image:
        if min(preview_image.size) <= 0:
            raise RuntimeError(f"Invalid preview dimensions: {preview_path.name}")
        preview_image.verify()
