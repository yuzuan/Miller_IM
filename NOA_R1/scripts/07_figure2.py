#!/usr/bin/env python3
# Purpose: draw patient-level MG-Inflam20 state differences and held-out gene and marker scores.
# Inputs: legacy Step 50 waterfall and two heatmap source CSVs. Recompute mode reads legacy Step 30 sample
# and QC tables, Step 48 effects, state and gene summaries, marker audit and cell calls, and the author XLS
# and 10x H5 files from NOA_PUBLIC_DATA_ROOT. The configured legacy root supplies the prepared output tables.
# Outputs: _panels/50_figure2_two_panel_final/Figure2/Figure2_{state_patient_waterfall,raw20_marker_double_heatmap}
# (.pdf/.png), Figure2_two_panel_preview.png and Figure2_two_panel_preview_grayscale.png. The matching
# _figure_build/50_figure2_two_panel_final/Figure2 directory holds Figure2_panel_manifest.csv and source_data/
# with Figure2_state_patient_waterfall_source.csv and Figure2_raw20_marker_double_heatmap_{raw20,marker_qc}_source.csv.
# Recompute mode also writes Figure2_legend.md, FINAL_RESULT.md and these source_data/ tables
# GSE278456_patient_marker_set_scores_source.csv, GSE278456_marker_set_assignment_qc_source.csv and
# GSE278456_step48_assignment_reproduction_source.csv.
# Run: python3 NOA_R1/scripts/07_figure2.py [--plot-only | --recompute-marker-qc]
# Keep Matplotlib and the original Scanpy scoring implementation to reproduce the released panels.

from __future__ import annotations

import argparse
import importlib
import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

figure_style = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
legacy_project_root = figure_style.legacy_root(revision_root)

# The exporter redirects these original output attributes to its scratch directory.
WRITE_ROOT = revision_root / "_figure_build/50_figure2_two_panel_final/Figure2"
FIG_ROOT = revision_root / "_panels/50_figure2_two_panel_final/Figure2"
SOURCE_OUT = WRITE_ROOT / "source_data"

# Use the October display order for both plots and the marker summary.
state_order = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]
raw20_genes = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
measurable_raw20_genes = [gene for gene in raw20_genes if gene != "AC253572.2"]


def configure_style() -> None:
    """Set the October Figure 2 typography and white background. Returns nothing."""
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


def save_panel(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    """Write the tightly cropped PDF master and derive the 600 dpi PNG."""
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    pdf = FIG_ROOT / f"{stem}.pdf"
    png = FIG_ROOT / f"{stem}.png"

    # The two panels use the original tight bounding box and 0.035 inch padding.
    fig.canvas.draw()
    fig.savefig(
        pdf,
        bbox_inches = "tight",
        pad_inches  = 0.035,
        metadata    = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None},
    )
    figure_style.render_pdf_png(pdf, png)
    plt.close(fig)
    output_paths = (pdf, png)
    return output_paths


def write_source(source_table: pd.DataFrame, filename: str) -> Path:
    """Write the original source-table schema and return its path."""
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    path = SOURCE_OUT / filename
    source_table.to_csv(path, index=False)
    return path


def validate_frozen_inputs(step48_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read Step 48 patient and gene effects and confirm marker exclusion."""
    effects_table = pd.read_csv(step48_root / "GSE278456_patient_state_effects.csv")
    summary_table = pd.read_csv(step48_root / "GSE278456_state_effect_summary.csv")
    gene_summary_table = pd.read_csv(step48_root / "GSE278456_gene_summary.csv")
    marker_audit_table = pd.read_csv(step48_root / "GSE278456_author_marker_audit.csv")

    # The frozen atlas has 108 eligible patient-state contrasts and a complete held-out gene grid.
    if len(effects_table) != 108 or effects_table["state"].nunique() != len(state_order):
        raise ValueError("Step 48 must contain 108 patient-state effects across the nine prespecified states")
    if set(summary_table["state"]) != set(state_order):
        raise ValueError("The Step 48 state list differs from the prespecified list")
    if len(gene_summary_table) != len(measurable_raw20_genes) * len(state_order):
        raise ValueError("The Step 48 held-out gene-by-state table is incomplete")
    if set(gene_summary_table["gene"]) != set(measurable_raw20_genes):
        raise ValueError("The Step 48 measurable MG-Inflam20 gene set differs from the fixed set")
    if gene_summary_table.groupby("gene")["state"].nunique().ne(len(state_order)).any():
        raise ValueError("Each measurable MG-Inflam20 gene must cover all nine states")

    # The nine displayed states use author markers with measurable MG-Inflam20 genes excluded.
    focus_audit_table = marker_audit_table.loc[marker_audit_table["state"].isin(state_order)]
    if focus_audit_table["state"].nunique() != len(state_order):
        raise ValueError("The author-marker audit is missing a prespecified state")
    used_markers = set()
    for value in focus_audit_table["markers_used"].fillna(""):
        used_markers.update(gene for gene in str(value).split(";") if gene)
    overlap = sorted(used_markers & set(measurable_raw20_genes))
    if overlap:
        raise ValueError(f"State-assignment markers overlap the held-out program: {overlap}")

    frozen_tables = (effects_table, summary_table, gene_summary_table)
    return frozen_tables


def recompute_marker_assignment_qc(step48_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Recompute patient marker scores and compare cell winners with the frozen Step 48 calls."""
    import scanpy as sc

    sc.settings.verbosity = 0

    # The prepared-data root can be separate from the Git checkout with the July marker reader.
    marker_script = legacy_project_root / "scripts/06b_cell_identity_and_state.py"
    if not marker_script.is_file():
        marker_script = revision_root.parent / "scripts/06b_cell_identity_and_state.py"
    if not marker_script.is_file():
        raise FileNotFoundError(f"Missing the July Figure 2 marker reader: {marker_script}")
    spec = importlib.util.spec_from_file_location("july_figure2_markers", marker_script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load the July Figure 2 marker reader: {marker_script}")
    marker_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(marker_module)

    # All 17 workbook states compete. The nine myeloid states are selected only after assignment.
    step30_relative = Path("write/30_dataset_rescue_search")
    public_metadata = figure_style.public_data_root(legacy_project_root) / step30_relative / "source_metadata"
    marker_book = public_metadata / "NIHMS2115262-supplement-S2_S5_S6.xls"
    h5_paths = sorted((public_metadata / "GSE278456_myeloid_h5").glob("*.h5"))
    if not marker_book.is_file():
        raise FileNotFoundError(f"Missing the GSE278456 author-marker workbook: {marker_book}")
    if len(h5_paths) != 21:
        raise ValueError(f"Expected 21 GSE278456 H5 files in {public_metadata}, found {len(h5_paths)}")
    markers, _ = marker_module.read_author_markers(marker_book, raw20_genes, state_order)
    if len(markers) != 17:
        raise ValueError(f"Expected 17 competing state-marker sets in the workbook, found {len(markers)}")

    # The manifest and QC table come from the prepared Step 30 analysis.
    step30_root = legacy_project_root / step30_relative
    manifest_table = pd.read_csv(step30_root / "GSE278456_sample_manifest.csv")
    manifest_table["geo_accession"] = manifest_table["geo_accession"].astype(str)
    qc_table = pd.read_csv(step30_root / "GSE278456_single_cell_qc.csv")
    qc_table = qc_table.set_index("geo_accession")

    # Preserve the original 120,766-cell assignment denominator and compare individual winners below.
    frozen_cells_table = pd.read_csv(
        step48_root / "GSE278456_cell_annotations_true_raw20.csv.gz",
        usecols=["cell_id", "geo_accession", "author_marker_state", "state_score"],
    )
    if len(frozen_cells_table) != 120766 or frozen_cells_table["cell_id"].duplicated().any():
        raise ValueError("Step 48 must contain the 120,766 unique frozen cell annotations")
    frozen_cells_table = frozen_cells_table.set_index("cell_id")

    all_marker_genes = set().union(*markers.values())
    patient_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    for h5_path in h5_paths:
        accession = h5_path.name.split("_")[0]
        sample_table = manifest_table.loc[manifest_table["geo_accession"].eq(accession)]
        if len(sample_table) != 1:
            raise ValueError(f"{accession} must match exactly one manifest row")
        manifest_row = sample_table.iloc[0]
        if not bool(manifest_row["is_idhwt_gbm"]):
            raise ValueError(f"Non-IDH-wild-type GBM in the frozen H5 input set: {accession}")
        sample = str(manifest_row["sample_title"])

        # Reapply the frozen bounds for detected genes, total counts and mitochondrial percentage.
        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        n_input = adata.n_obs
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].between(200, 8000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        adata = adata[keep].copy()
        n_qc = adata.n_obs
        expected_qc = qc_table.loc[accession]
        if int(expected_qc["n_input"]) != n_input or int(expected_qc["n_qc_pass"]) != n_qc:
            raise ValueError(f"{accession} differs from the frozen cell counts before or after QC")

        # Confirm the frozen gene coverage and use the original log-normalized counts.
        present_raw20 = [gene for gene in raw20_genes if gene in adata.var_names]
        if present_raw20 != measurable_raw20_genes:
            raise RuntimeError(f"{accession} differs from the frozen MG-Inflam20 gene coverage or order")
        adata.obs_names = [f"{accession}:{barcode}" for barcode in adata.obs_names.astype(str)]
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)

        # Exclude the held-out program and every author-marker gene from the neutral pool.
        excluded_controls = set(raw20_genes) | all_marker_genes
        neutral_control_pool = [gene for gene in adata.var_names if gene not in excluded_controls]
        if len(neutral_control_pool) < 5000:
            raise RuntimeError(f"{accession} has fewer than 5,000 state-neutral control genes")

        # Use the original control bins and seed for each competing marker set.
        score_columns: list[str] = []
        for state, genes in markers.items():
            use = [gene for gene in genes if gene in adata.var_names]
            if len(use) < 30:
                raise RuntimeError(f"{accession} {state} has fewer than 30 measurable markers")
            column = f"state_score_{state}"
            sc.tl.score_genes(
                adata,
                gene_list    = use,
                gene_pool    = neutral_control_pool + use,
                ctrl_size    = 50,
                n_bins       = 25,
                score_name   = column,
                random_state = 20260713,
                ctrl_as_ref  = False,
                use_raw      = False,
            )
            score_columns.append(column)

        # Preserve the original argsort winner selection, including its handling of ties.
        score_matrix = adata.obs[score_columns].to_numpy(float)
        order = np.argsort(score_matrix, axis=1)
        state_names = np.array(list(markers), dtype=object)
        winner_index = order[:, -1]
        winners = state_names[winner_index]
        winning_score = score_matrix[np.arange(adata.n_obs), winner_index]

        frozen_sample_table = frozen_cells_table.loc[adata.obs_names]
        winner_matches = frozen_sample_table["author_marker_state"].to_numpy(object) == winners
        score_difference = np.abs(frozen_sample_table["state_score"].to_numpy(float) - winning_score)
        if not bool(winner_matches.all()) or float(score_difference.max()) > 1e-10:
            raise RuntimeError(f"{accession} differs from the frozen Step 48 state assignment")
        audit_rows.append(
            {
                "geo_accession": accession,
                "sample": sample,
                "n_input": n_input,
                "n_qc_pass": n_qc,
                "n_winner_exact_match": int(winner_matches.sum()),
                "winner_match_fraction": float(winner_matches.mean()),
                "max_abs_winner_score_difference": float(score_difference.max()),
                "n_competing_marker_sets": len(markers),
            }
        )

        # Score each eligible patient-state once, retaining the 30-cell cutoff.
        focus = np.isin(winners, state_order)
        focus_scores = {state: adata.obs[f"state_score_{state}"].to_numpy(float) for state in state_order}
        for assigned_state in state_order:
            target = focus & (winners == assigned_state)
            n_cells = int(target.sum())
            if n_cells < 30:
                continue
            for tested_marker_set in state_order:
                patient_rows.append(
                    {
                        "geo_accession": accession,
                        "sample": sample,
                        "assigned_state": assigned_state,
                        "tested_marker_set": tested_marker_set,
                        "n_cells": n_cells,
                        "mean_marker_set_score": float(focus_scores[tested_marker_set][target].mean()),
                    }
                )
        del adata, score_matrix

    # Reconcile every cell with the Step 48 archive and use the original patient-state denominator.
    patient_scores_table = pd.DataFrame(patient_rows)
    audit_table = pd.DataFrame(audit_rows).sort_values("geo_accession")
    if int(audit_table["n_qc_pass"].sum()) != len(frozen_cells_table):
        raise ValueError("The H5 samples do not cover the entire frozen Step 48 cell table")
    if not audit_table["winner_match_fraction"].eq(1).all():
        raise ValueError("Some cell-level marker winners differ from the frozen Step 48 calls")
    if len(patient_scores_table) != 972:
        raise ValueError(f"Expected 972 frozen patient-state-marker rows, found {len(patient_scores_table)}")

    # Equal-weight patient means within each tested marker and assigned state form the classification QC block.
    marker_summary_table = (
        patient_scores_table.groupby(["tested_marker_set", "assigned_state"], observed=True)
        .agg(
            n_patients = ("sample", "nunique"),
            patient_balanced_mean_marker_set_score = ("mean_marker_set_score", "mean"),
        )
        .reset_index()
    )

    def standardize_marker_scores(scores: pd.Series) -> pd.Series | int:
        """Return within-marker-set population Z-scores, with zero for a constant row."""
        spread = scores.std(ddof=0)
        if spread > 0:
            standardized = (scores - scores.mean()) / spread
        else:
            standardized = 0
        return standardized

    marker_groups = marker_summary_table.groupby("tested_marker_set", observed=True)
    marker_means = marker_groups["patient_balanced_mean_marker_set_score"]
    marker_summary_table["within_marker_set_z"] = marker_means.transform(standardize_marker_scores)
    if len(marker_summary_table) != len(state_order) ** 2:
        raise ValueError("The marker-set QC summary must have one row per marker-set and state combination")

    marker_tables = (patient_scores_table, marker_summary_table, audit_table)
    return marker_tables


def prepare_waterfall_source(effects_table: pd.DataFrame, summary_table: pd.DataFrame) -> pd.DataFrame:
    """Return one row per patient-state contrast in the prespecified display order."""
    indexed_summary_table = summary_table.set_index("state")
    state_tables: list[pd.DataFrame] = []
    for state in state_order:
        state_table = effects_table.loc[effects_table["state"].eq(state)].sort_values("raw20_delta").copy()
        state_table["within_state_rank"] = np.arange(1, len(state_table) + 1)
        state_table["direction"] = np.where(state_table["raw20_delta"] > 0, "Above rest", "Below rest")
        state_table["n_patients"] = int(indexed_summary_table.loc[state, "n_patients"])
        state_table["n_positive"] = int(indexed_summary_table.loc[state, "n_positive"])
        state_table["mean_delta"] = float(indexed_summary_table.loc[state, "mean_delta"])
        state_table["fdr_9_states"] = float(indexed_summary_table.loc[state, "fdr_9_states"])
        state_table["score_pseudobulk_mean_direction_concordant"] = bool(
            indexed_summary_table.loc[state, "score_pseudobulk_mean_direction_concordant"]
        )
        state_tables.append(state_table)

    waterfall_table_total = pd.concat(state_tables, ignore_index=True)
    return waterfall_table_total


def make_patient_waterfall(waterfall_table: pd.DataFrame, *, final_layout: bool = False) -> dict[str, str]:
    """Write the patient waterfall table and plot one facet per prespecified state."""
    stem = "Figure2_state_patient_waterfall"
    source = write_source(waterfall_table, f"{stem}_source.csv")
    summary_table = waterfall_table.drop_duplicates("state").set_index("state")

    # The final print layout has more vertical space and keeps the same patient slots and score limits.
    panel_height = 89 if final_layout else 85
    fig, axes = plt.subplots(3, 3, figsize=(183 * (1 / 25.4), panel_height * (1 / 25.4)), sharex=True, sharey=True)
    y_min, y_max = -0.76, 0.56
    state_colors = {
        "MCG1": "#00A087",
        "MCG2": "#91D1C2",
        "MCG3": "#4DBBD5",
        "MCG4": "#8491B4",
        "MCG5": "#3C5488",
        "MAC1": "#F39B7F",
        "MAC2": "#E64B35",
        "E-MDSC": "#B09C85",
        "M-MDSC": "#7E6148",
    }

    def format_fdr(value: float) -> str:
        """Format a state FDR with the original precision used in the inset."""
        if value < 0.001:
            formatted = f"{value:.2e}"
        elif value < 0.01:
            formatted = f"{value:.4f}"
        else:
            formatted = f"{value:.3f}"
        return formatted

    # Bars are patient contrasts, with the dashed line showing the equal-weight patient mean.
    for index, (ax, state) in enumerate(zip(axes.flat, state_order, strict=True)):
        state_table = waterfall_table.loc[waterfall_table["state"].eq(state)]
        x = state_table["within_state_rank"].to_numpy(int)
        y = state_table["raw20_delta"].to_numpy(float)
        colors = np.where(y > 0, "#E64B35", "#3C5488")
        ax.axhspan(0, y_max, color="#FFF7F4", zorder=0)
        ax.axhspan(y_min, 0, color="#F4F6FA", zorder=0)
        ax.axhline(0, color="#222222", linewidth=0.65, zorder=2)
        ax.bar(x, y, width=0.76, color=colors, edgecolor="none", alpha=0.88, zorder=3)
        mean = float(summary_table.loc[state, "mean_delta"])
        ax.axhline(mean, color=state_colors[state], linewidth=1.05, linestyle=(0, (3, 2)), zorder=4)
        ax.set_xlim(0.25, 21.75)
        ax.set_ylim(y_min, y_max)

        # A dagger flags facets with fewer than six evaluable patients, as in the caption.
        n_patients = int(summary_table.loc[state, "n_patients"])
        n_positive = int(summary_table.loc[state, "n_positive"])
        fdr_label = format_fdr(float(summary_table.loc[state, "fdr_9_states"]))
        state_title = f"{state}{'†' if n_patients < 6 else ''}"
        ax.set_title(state_title, loc="left", color=state_colors[state], fontweight="bold", pad=7)
        ax.text(
            1.0,
            1.025,
            f"{n_positive}/{n_patients}  |  FDR {fdr_label}",
            transform = ax.transAxes,
            ha        = "right",
            va        = "bottom",
            fontsize  = 7.0,
            color     = "#222222",
        )
        if state == "MCG3":
            for spine in ax.spines.values():
                spine.set_color("#868686")
                spine.set_linestyle((0, (3, 2)))
        else:
            ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#222222")
        ax.tick_params(length=2.2, width=0.65)
        if index % 3 != 0:
            ax.tick_params(axis="y", left=False, labelleft=False)
        elif not final_layout:
            ax.set_ylabel("State − other myeloid\nMG-Inflam20 score")
        if index >= 6:
            ax.tick_params(axis="x", labelsize=7.0)
            ax.set_xticks([1, 7, 14, 21])
            ax.set_xlabel("Patients ranked within state" if index == 7 else "")
        else:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)

    # The October submission uses one vertical score label for all nine facets.
    if final_layout:
        fig.subplots_adjust(left=0.09, right=0.995, bottom=0.18, top=0.96, wspace=0.13, hspace=0.43)
        fig.text(
            0.025,
            0.55,
            "State − other myeloid MG-Inflam20 score",
            ha       = "center",
            va       = "center",
            rotation = 90,
            fontsize = 9.3,
            color    = "#222222",
        )
    else:
        fig.subplots_adjust(left=0.09, right=0.995, bottom=0.13, top=0.965, wspace=0.13, hspace=0.29)
    pdf, png = save_panel(fig, stem)
    panel_record = {"panel": "2A", "source": str(source), "pdf": str(pdf), "png": str(png)}
    return panel_record


def make_double_heatmap(
    gene_summary_table: pd.DataFrame,
    marker_summary_table: pd.DataFrame,
    state_summary_table: pd.DataFrame,
    *,
    final_layout: bool = False,
) -> dict[str, str]:
    """Write the held-out expression and marker QC source tables and draw the two heatmaps."""
    stem = "Figure2_raw20_marker_double_heatmap"
    gene_source_table = gene_summary_table.copy()
    gene_source_table["block"] = "Held-out MG-Inflam20 genes; not used for state assignment"
    marker_source_table = marker_summary_table.copy()
    marker_source_table["block"] = "Author marker-set assignment scores; classification QC"
    source_top = write_source(gene_source_table, f"{stem}_raw20_source.csv")
    source_bottom = write_source(marker_source_table, f"{stem}_marker_qc_source.csv")

    expression_table = gene_summary_table.pivot(index="gene", columns="state", values="within_gene_z").reindex(
        index=measurable_raw20_genes, columns=state_order
    )
    marker_z_table = marker_summary_table.pivot(
        index="tested_marker_set", columns="assigned_state", values="within_marker_set_z"
    ).reindex(index=state_order, columns=state_order)
    if expression_table.isna().any().any() or marker_z_table.isna().any().any():
        raise RuntimeError("Missing cells in the two-block heat map")

    n_by_state = state_summary_table.set_index("state")["n_patients"].astype(int).to_dict()
    state_labels = [f"{state}\n(n={n_by_state[state]}){'†' if n_by_state[state] < 6 else ''}" for state in state_order]

    # The final print layout reallocates space between the two blocks without resizing the page.
    fig = plt.figure(figsize=(183 * (1 / 25.4), 128 * (1 / 25.4)))
    grid = fig.add_gridspec(
        nrows         = 4,
        ncols         = 2,
        width_ratios  = [1, 0.034],
        height_ratios = [2.2, 16.5, 1.8, 7.05] if final_layout else [2.2, 19, 1.35, 5.3],
        left          = 0.16,
        right         = 0.94,
        bottom        = 0.10,
        top           = 0.98,
        wspace        = 0.16,
        hspace        = 0.05,
    )

    # Allocate the title, both blocks, the separator and their separate color scales.
    title_top = fig.add_subplot(grid[0, 0])
    ax_top = fig.add_subplot(grid[1, 0])
    separator = fig.add_subplot(grid[2, 0])
    ax_bottom = fig.add_subplot(grid[3, 0])
    cax_top = fig.add_subplot(grid[1, 1])
    cax_bottom = fig.add_subplot(grid[3, 1])

    title_top.axis("off")
    title_top.text(
        0,
        0.93,
        "MG-Inflam20 gene expression",
        ha         = "left",
        va         = "center",
        fontsize   = 7.4,
        fontweight = "bold",
        color      = "#222222",
    )
    sns.heatmap(
        expression_table,
        ax          = ax_top,
        cmap        = LinearSegmentedColormap.from_list("raw20_npg", ["#3C5488", "#F7F7F7", "#F39B7F"]),
        vmin        = -2.5,
        vmax        = 2.5,
        center      = 0,
        linewidths  = 0.45,
        linecolor   = "white",
        cbar        = True,
        cbar_ax     = cax_top,
        xticklabels = state_labels,
        yticklabels = measurable_raw20_genes,
    )
    ax_top.xaxis.tick_top()
    ax_top.tick_params(axis="x", rotation=0, length=0, pad=2)
    ax_top.tick_params(axis="y", rotation=0, length=0, pad=2)
    ax_top.set_xlabel("")
    ax_top.set_ylabel("")
    cax_top.set_ylabel("Within-gene\nZ-score", fontsize=9.3, labelpad=4, color="#222222")
    cax_top.tick_params(labelsize=6.2, length=2)

    # The lower block displays assignment QC with the original gray palette.
    separator.axis("off")
    separator.axhline(0.92, color="#222222", linewidth=0.8)
    separator.text(
        0,
        0.55,
        "State-marker scores",
        ha         = "left",
        va         = "center",
        fontsize   = 6.8,
        fontweight = "bold",
        color      = "#222222",
    )
    sns.heatmap(
        marker_z_table,
        ax          = ax_bottom,
        cmap        = LinearSegmentedColormap.from_list("marker_qc", ["#ECECEC", "#A9B3BE", "#4B5563"]),
        vmin        = -2.5,
        vmax        = 2.5,
        center      = 0,
        linewidths  = 0.45,
        linecolor   = "white",
        cbar        = True,
        cbar_ax     = cax_bottom,
        xticklabels = False,
        yticklabels = [f"{state} markers" for state in state_order],
    )
    ax_bottom.tick_params(axis="y", rotation=0, length=0, pad=2, colors="#222222", labelsize=7.0)
    ax_bottom.tick_params(axis="x", bottom=False)
    ax_bottom.set_xlabel("")
    ax_bottom.set_ylabel("")
    cax_bottom.set_ylabel("Within-marker-set\nZ-score", fontsize=9.3, labelpad=4, color="#222222")
    cax_bottom.tick_params(labelsize=7.2, length=2, colors="#222222")
    for spine in cax_bottom.spines.values():
        spine.set_edgecolor("#868686")

    pdf, png = save_panel(fig, stem)
    panel_record = {
        "panel": "2B",
        "source": f"{source_top};{source_bottom}",
        "pdf": str(pdf),
        "png": str(png),
    }
    return panel_record


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path]:
    """Build the two review images from the existing PDF-derived panel PNGs."""
    target_width = 1800
    resized: list[Image.Image] = []
    for record in records:
        with Image.open(record["png"]) as panel_image:
            image = panel_image.convert("RGB")
        height = round(image.height * target_width / image.width)
        resized.append(image.resize((target_width, height), Image.Resampling.LANCZOS))

    margin = 85
    gap = 40
    canvas = Image.new(
        "RGB",
        (target_width + margin * 2, sum(image.height for image in resized) + gap + margin * 2),
        "white",
    )
    draw = ImageDraw.Draw(canvas)

    # Keep the original figure-review font order and fallback on other platforms.
    font_candidates = [
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    font_path = next((Path(path) for path in font_candidates if Path(path).exists()), None)
    font = ImageFont.truetype(str(font_path), 54) if font_path else ImageFont.load_default()
    y = margin
    for label, image in zip(["A", "B"], resized, strict=True):
        draw.text((22, y + 2), label, fill="#222222", font=font)
        canvas.paste(image, (margin, y))
        y += image.height + gap
    preview = FIG_ROOT / "Figure2_two_panel_preview.png"
    canvas.save(preview, dpi=(150, 150))
    grayscale = FIG_ROOT / "Figure2_two_panel_preview_grayscale.png"
    canvas.convert("L").convert("RGB").save(grayscale, dpi=(150, 150))
    for image_path in (preview, grayscale):
        with Image.open(image_path) as image:
            if image.size != canvas.size:
                raise ValueError(f"Unexpected review image size: {image_path.name}")
            image.verify()
    previews = (preview, grayscale)
    return previews


def write_legend() -> Path:
    """Write the frozen October Figure 2 legend and return its path."""
    path = WRITE_ROOT / "Figure2_legend.md"
    path.write_text(
        "# Figure 2 | The MG-Inflam20 program preferentially maps to MCG1/MCG2-like states in an independent "
        "cross-sectional myeloid atlas\n"
        "\n"
        "**A.** Patient-level state-versus-rest waterfall plots in the 21 IDH-wild-type primary GBM samples "
        "of GSE278456. Each bar is one patient with at least 30 cells in the indicated author-marker-derived "
        "state and at least 30 cells in the remaining prespecified myeloid states. Bars show the difference "
        "in control-adjusted MG-Inflam20 score between that state and the remaining states; dashed lines show "
        "patient-equal means. Red and blue denote values above and below the remaining states, respectively. "
        "Numbers give patients above the remaining states / evaluable patients and Benjamini–Hochberg FDR "
        "across the nine state-wise exact two-sided sign-flip tests. All facets use the same 21 patient slots "
        "and y-axis. MCG3 is marked because its mean score direction differs from the raw-count pseudobulk "
        "sensitivity analysis.\n"
        "\n"
        "**B.** Two-block patient-balanced heatmap across the same fixed state order. The upper block shows "
        "relative expression of the 19 measurable MG-Inflam20 genes; values are patient-equal mean "
        "log-normalized expression standardized within each gene across the nine states. All measurable "
        "MG-Inflam20 genes were removed from the author marker sets before state assignment, so this block "
        "was held out from direct state-label construction. The lower, visually reduced block shows the nine "
        "complete author marker-set scores, averaged within patient and assigned state and then standardized "
        "within marker set. These scores were used to create the state labels; their diagonal structure is "
        "expected by construction and is shown only as classification quality control, not as independent "
        "validation. Direct overlap between the final label-defining marker sets and the measurable "
        "MG-Inflam20 genes was zero. State labels show the number of evaluable patients; dagger marks fewer "
        "than six patients.\n"
        "\n"
        "GSE278456 is cross-sectional and is used only to localize the recurrence-associated program across "
        "myeloid states; recurrence evidence is provided separately in Figure 1.\n",
        encoding="utf-8",
    )
    return path


def main(*, recompute_marker_qc: bool = False, final_layout: bool = False) -> None:
    """Redraw the frozen October tables, or recompute the marker QC on request. Returns nothing."""
    configure_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    WRITE_ROOT.mkdir(parents=True, exist_ok=True)

    # The final export reads the frozen Step 50 sources without opening any raw H5 files.
    if recompute_marker_qc:
        step48_root = legacy_project_root / "write/48_figure2_identity_state_rebuild/Figure2"
        effects_table, state_summary_table, gene_summary_table = validate_frozen_inputs(step48_root)
        waterfall_table = prepare_waterfall_source(effects_table, state_summary_table)
        marker_patient_table, marker_summary_table, audit_table = recompute_marker_assignment_qc(step48_root)
    else:
        frozen_root = legacy_project_root / "write/50_figure2_two_panel_final/Figure2/source_data"
        waterfall_table = pd.read_csv(frozen_root / "Figure2_state_patient_waterfall_source.csv")
        gene_summary_table = pd.read_csv(frozen_root / "Figure2_raw20_marker_double_heatmap_raw20_source.csv")
        marker_summary_table = pd.read_csv(frozen_root / "Figure2_raw20_marker_double_heatmap_marker_qc_source.csv")
        state_summary_table = waterfall_table.drop_duplicates("state")[["state", "n_patients"]]

    # The October tables contain 108 eligible contrasts and complete expression and marker grids.
    if len(waterfall_table) != 108:
        raise ValueError("The Figure 2 waterfall must contain 108 frozen patient-state contrasts")
    if len(gene_summary_table) != len(measurable_raw20_genes) * len(state_order):
        raise ValueError("The Figure 2 held-out gene table does not cover the frozen gene-by-state grid")
    if len(marker_summary_table) != len(state_order) ** 2:
        raise ValueError("The Figure 2 marker score table does not cover the frozen state-by-state grid")

    # Retain the original Step 50 source filenames when recalculating marker QC.
    if recompute_marker_qc:
        write_source(marker_patient_table, "GSE278456_patient_marker_set_scores_source.csv")
        write_source(marker_summary_table, "GSE278456_marker_set_assignment_qc_source.csv")
        write_source(audit_table, "GSE278456_step48_assignment_reproduction_source.csv")
        write_legend()

    records = [
        make_patient_waterfall(waterfall_table, final_layout=final_layout),
        make_double_heatmap(gene_summary_table, marker_summary_table, state_summary_table, final_layout=final_layout),
    ]
    preview, grayscale = make_preview(records)
    manifest_table = pd.DataFrame(records)
    manifest_path = WRITE_ROOT / "Figure2_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)
    if not manifest_path.stat().st_size:
        raise ValueError("The Figure 2 panel manifest is empty")

    # Summarize the recalculation using observed counts and the original report path.
    if recompute_marker_qc:
        matched_cells = int(audit_table["n_winner_exact_match"].sum())
        qc_cells = int(audit_table["n_qc_pass"].sum())
        max_score_difference = audit_table["max_abs_winner_score_difference"].max()
        report = (
            "# Figure 2 marker-score reconstruction\n\n"
            f"- Exact winner matches: {matched_cells:,}/{qc_cells:,}. "
            f"Maximum absolute score difference: {max_score_difference:.3g}.\n"
            f"- Source rows: {len(waterfall_table)} waterfall, {len(gene_summary_table)} held-out genes, "
            f"{len(marker_patient_table)} patient marker scores, and {len(marker_summary_table)} marker QC means.\n"
            f"- Review images: `{preview.name}` and `{grayscale.name}`.\n\n"
            f"{manifest_table.to_markdown(index=False)}\n"
        )
        (WRITE_ROOT / "FINAL_RESULT.md").write_text(report, encoding="utf-8")
    print(f"Figure 2 panels: {len(records)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw the October Figure 2 panels")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--plot-only",
        dest   = "recompute_marker_qc",
        action = "store_false",
        help   = "Draw from the frozen Step 50 source tables (default).",
    )
    mode.add_argument(
        "--recompute-marker-qc",
        action = "store_true",
        help   = "Recompute marker QC from the original GSE278456 inputs and draw the panels.",
    )
    parser.set_defaults(recompute_marker_qc=False)
    arguments = parser.parse_args()
    main(recompute_marker_qc=arguments.recompute_marker_qc)
