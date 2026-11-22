#!/usr/bin/env python3
# Purpose: draw the cohort workflow, clean-myeloid UMAPs, rank comparison and candidate state figures.
# Inputs: prepared Step 38 paired-cohort audit, Step 41 rank and patient score tables, Step 43 UMAP sources,
# and the Step 36 GSE274546 cluster annotation, all under NOA_LEGACY_ROOT. See each input read for the exact filename.
# Outputs: _panels/42_visual_story_rebuild/{Figure1,Figure2} contains candidate panel PDF/PNG pairs.
# The matching _figure_build/42_visual_story_rebuild/{Figure1,Figure2} directories contain the two
# Figure{1,2}_panel_manifest.csv files and source_data/{panel stem}_source.csv for each plotted panel.
# Run: python3 NOA_R1/scripts/06a_figure1_evidence_layout.py
# Keep Matplotlib to retain the workflow diagram, plot scales and the original annotation geometry.

from __future__ import annotations

import importlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr

figure_style = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
legacy_project_root = figure_style.legacy_root(revision_root)
mm = figure_style.mm

# The Figure 1 exporter redirects both output roots for the isolated build.
FIG_ROOT = revision_root / "_panels/42_visual_story_rebuild"
WRITE_ROOT = revision_root / "_figure_build/42_visual_story_rebuild"

# Reuse the same two colors, labels and RNG seed across the related candidate panels.
datasets = ["GSE174554", "GSE274546"]
dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}
text_color = "#252525"
neutral_color = "#717171"
grid_color = "#E6E6E6"
jitter_seed = 20260713


def configure_style() -> None:
    """Set the original candidate figure's font, stroke widths and colors. Returns nothing."""
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.8,
            "axes.linewidth": 0.65,
            "axes.edgecolor": text_color,
            "axes.labelcolor": text_color,
            "text.color": text_color,
            "xtick.color": text_color,
            "ytick.color": text_color,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def write_source(source_table: pd.DataFrame, figure_name: str, stem: str) -> Path:
    """Save the plotted source and return its original filename."""
    out = WRITE_ROOT / figure_name / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}_source.csv"
    source_table.to_csv(path, index=False)
    return path


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, object]:
    """Return the source and image file fields for the existing candidate manifest."""
    manifest_row = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source_csv": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": figure_style.sha256(source),
        "pdf_sha256": figure_style.sha256(pdf),
        "png_sha256": figure_style.sha256(png),
        "png_dpi": 600,
    }
    return manifest_row


def label_position(source_table: pd.DataFrame) -> tuple[float, float]:
    """Return the median coordinates within the central 60% of the cluster."""
    points = source_table[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
    center = np.nanmedian(points, axis=0)
    distance = np.sqrt(((points - center) ** 2).sum(axis=1))
    core = points[distance <= np.nanquantile(distance, 0.60)]
    target = np.nanmedian(core, axis=0) if len(core) else center
    coordinates = (float(target[0]), float(target[1]))
    return coordinates


def figure1_design(*, final_layout: bool = False) -> dict[str, object]:
    """Draw the paired-cohort workflow and return the Figure 1A manifest row."""
    stem = "Figure1A_compact_paired_design"
    pair_audit_path = (
        legacy_project_root
        / "write/38_independent_cohort_mg_inflammatory_recalculation"
        / "independent_input_pair_audit.csv"
    )
    audit_table = pd.read_csv(pair_audit_path)
    source_table = audit_table.loc[audit_table["dataset"].isin(datasets) & audit_table["threshold"].eq(20)].copy()
    if source_table.groupby("dataset").size().to_dict() != {
        "GSE174554": 1,
        "GSE274546": 1,
    }:
        raise ValueError("Step38 must contain one formal threshold-20 audit row per cohort.")

    # The 91 and 111 input units are the frozen design counts before reconstruction and pairing.
    source_table["input_libraries_or_matrices"] = source_table["dataset"].map({"GSE174554": 91, "GSE274546": 111})
    source_table["clean_myeloid"] = source_table["n_cells_primary"].astype(int) + source_table[
        "n_cells_recurrent"
    ].astype(int)
    source_table["formal_pairs_threshold20"] = source_table["n_pairs"].astype(int)
    observed = (
        source_table.set_index("dataset")[["formal_pairs_threshold20", "clean_myeloid"]].astype(int).to_dict("index")
    )
    expected = {
        "GSE174554": {"formal_pairs_threshold20": 17, "clean_myeloid": 12489},
        "GSE274546": {"formal_pairs_threshold20": 45, "clean_myeloid": 55568},
    }
    if observed != expected:
        raise ValueError(f"Unexpected IDH-restricted Figure1A analysis set: {observed}")
    source_table = source_table.sort_values("dataset").reset_index(drop=True)
    source_table["modality"] = "snRNA-seq"
    source_table["raw_input_unit"] = ["10x matrices", "10x libraries"]
    source_table["minimum_cells_per_condition"] = 20

    # Select the fields used to draw or interpret the diagram.
    source_table["processing_scope"] = "Each cohort reconstructed, annotated and tested independently"
    source_table["statistical_testing"] = (
        "Patient pseudobulk, paired edgeR and MG-Inflam20 GSEA performed separately within each cohort"
    )
    source_table["replication_assessment"] = "Cohort results compared only after independent testing"
    source_table["workflow"] = (
        "all-cell reconstruction -> pan-myeloid re-clustering -> marker and cluster-DEG review -> "
        "patient raw-count pseudobulk -> paired edgeR -> prespecified MG-Inflam20 program GSEA"
    )
    safe_source_columns = [
        "dataset",
        "modality",
        "raw_input_unit",
        "input_libraries_or_matrices",
        "n_pairs_before_idh",
        "n_pairs_excluded_idh",
        "excluded_pair_keys",
        "formal_pairs_threshold20",
        "n_cells_primary",
        "n_cells_recurrent",
        "clean_myeloid",
        "minimum_cells_per_condition",
        "processing_scope",
        "statistical_testing",
        "replication_assessment",
        "workflow",
    ]
    source_table = source_table.loc[:, safe_source_columns]
    source = write_source(source_table, "Figure1", stem)

    # The final export gives the two cohort tracks and the workflow boxes more vertical space.
    if final_layout:
        fig = plt.figure(figsize=(mm(183), mm(74)), facecolor="white")
    else:
        fig = plt.figure(figsize=(mm(183), mm(64)), facecolor="white")
    ax = fig.add_axes([0.025, 0.025, 0.95, 0.95])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def rounded_box(x0: float, x1: float, y0: float, y1: float, edge: str, fill: str | tuple, lw: float) -> None:
        """Draw a workflow box in axes coordinates. Returns nothing."""
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                boxstyle  = "round,pad=0.02,rounding_size=0.7",
                facecolor = fill,
                edgecolor = edge,
                linewidth = lw,
            )
        )

    # Both layouts share the two horizontal tracks. The final export lowers the pairing labels.
    cohort_boxes = [(4.0, 47.5), (52.5, 96.0)]
    cohort_rows = source_table.to_dict("records")
    title_y, subtitle_y = (95.0, 90.0) if final_layout else (94.5, 89.3)
    pair_bottom, pair_top, pair_text_y = (80.7, 85.7, 83.2) if final_layout else (81.7, 86.7, 84.2)
    for row, (x0, x1) in zip(cohort_rows, cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        rounded_box(x0, x1, 80.0, 97.0, color, pale, lw=1.0)
        ax.text(
            center,
            title_y,
            f"{row['dataset']}  |  snRNA-seq",
            fontsize = 8.3,
            weight   = "bold",
            color    = color,
            ha       = "center",
            va       = "top",
        )
        input_label = f"{int(row['input_libraries_or_matrices'])} {row['raw_input_unit']}"
        ax.text(center, subtitle_y, input_label, fontsize=6.1, color=neutral_color, ha="center", va="top")
        primary_x, recurrent_x = center - 8.0, center + 8.0
        rounded_box(primary_x - 4.9, primary_x + 4.9, pair_bottom, pair_top, color, "white", lw=0.9)
        rounded_box(recurrent_x - 5.5, recurrent_x + 5.5, pair_bottom, pair_top, color, color, lw=0.9)
        ax.text(primary_x, pair_text_y, "Primary", fontsize=6.3, color=text_color, ha="center", va="center")
        ax.text(recurrent_x, pair_text_y, "1st recurrence", fontsize=6.3, color="white", ha="center", va="center")
        ax.annotate(
            "",
            xy = (recurrent_x - 5.8, pair_text_y),
            xytext = (primary_x + 5.2, pair_text_y),
            arrowprops = dict(arrowstyle="-|>", color=color, lw=1.0),
        )

    # The boxes track cohort reconstruction, pan-myeloid reclustering and marker review.
    for center in [sum(box) / 2 for box in cohort_boxes]:
        ax.annotate(
            "",
            xy = (center, 77.0 if final_layout else 74.8),
            xytext = (center, 79.5),
            arrowprops = dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85),
        )
    if final_layout:
        reconstruction_y, reclustering_y, marker_review_y = (64.5, 76.5), (49.0, 61.0), (33.5, 45.5)
    else:
        reconstruction_y, reclustering_y, marker_review_y = (62.5, 75.5), (48.5, 61.5), (34.5, 47.5)
    workflow_boxes = [
        (
            8.0,
            92.0,
            *reconstruction_y,
            "Cohort-specific all-cell reconstruction",
            "Library QC · doublet removal · basic lineage annotation",
        ),
        (
            8.0,
            92.0,
            *reclustering_y,
            "Cohort-wise pan-myeloid re-clustering",
            "2,000 HVGs · PCA · neighbor graph · Leiden",
        ),
        (
            6.0,
            94.0,
            *marker_review_y,
            "Marker and cluster-DEG review",
            "Remove lymphoid/tumor contaminants and low-quality cells",
        ),
    ]
    for index, (x0, x1, y0, y1, title, subtitle) in enumerate(workflow_boxes):
        rounded_box(x0, x1, y0, y1, "#858585", "white", lw=0.75)
        ax.text(50, y1 - 2.2, title, fontsize=7.3, weight="bold", ha="center", va="top")
        ax.text(50, y0 + 2.2, subtitle, fontsize=6.1, color=neutral_color, ha="center", va="bottom")
        if index == len(workflow_boxes) - 1:
            continue

        # The print layout places arrow ends at the new box boundaries.
        next_y1 = workflow_boxes[index + 1][3]
        arrow_gap = 0.1 if final_layout else 0.5
        arrow_style = dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85)
        if final_layout:
            arrow_style.update(shrinkA=0, shrinkB=0)
        ax.annotate(
            "",
            xy = (50, next_y1 + arrow_gap),
            xytext = (50, y0 - arrow_gap),
            arrowprops = arrow_style,
        )

    # The cohort tracks join only after the paired pseudobulk and enrichment tests.
    analysis_y0, analysis_y1 = (16.0, 30.0) if final_layout else (14.5, 29.5)
    summary_title_y, cell_count_y, method_y = (27.5, 23.0, 18.7) if final_layout else (26.8, 22.0, 17.2)
    for row, (x0, x1) in zip(cohort_rows, cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        source_x = 35.0 if center < 50 else 65.0
        arc = "arc3,rad=0.08" if center < 50 else "arc3,rad=-0.08"
        ax.annotate(
            "",
            xy = (center, analysis_y1 + 0.6),
            xytext = (source_x, 32.5 if final_layout else 35.5),
            arrowprops = dict(arrowstyle="-|>", color=color, lw=0.9, connectionstyle=arc),
        )
        rounded_box(x0, x1, analysis_y0, analysis_y1, color, pale, lw=1.0)
        pairs_label = f"{row['dataset']}  |  {int(row['formal_pairs_threshold20'])} paired patients"
        ax.text(
            center, summary_title_y, pairs_label, fontsize=7.1, weight="bold", color=color, ha="center", va="center"
        )
        cell_count = f"{int(row['clean_myeloid']):,} analyzed myeloid cells"
        ax.text(center, cell_count_y, cell_count, fontsize=5.9, color=neutral_color, ha="center", va="center")
        ax.text(
            center,
            method_y,
            "Raw-count pseudobulk · paired edgeR · MG-Inflam20 camera and GSEA",
            fontsize = 5.9,
            color    = text_color,
            ha       = "center",
            va       = "center",
        )

    left_target, right_target = ((45.0, 12.4), (55.0, 12.4)) if final_layout else ((46.8, 10.8), (53.2, 10.8))
    ax.annotate(
        "",
        xy=left_target,
        xytext=(sum(cohort_boxes[0]) / 2, analysis_y0 - 0.5),
        arrowprops=dict(arrowstyle="-|>", color=dataset_colors["GSE174554"], lw=0.9, connectionstyle="arc3,rad=-0.08"),
    )
    ax.annotate(
        "",
        xy = right_target,
        xytext = (sum(cohort_boxes[1]) / 2, analysis_y0 - 0.5),
        arrowprops = dict(arrowstyle="-|>", color=dataset_colors["GSE274546"], lw=0.9, connectionstyle="arc3,rad=0.08"),
    )
    rounded_box(22.0, 78.0, 1.0, 11.8, "#159D88", mpl.colors.to_rgba("#159D88", 0.085), lw=1.0)
    ax.text(
        50, 8.8, "Paired recurrence analyses", fontsize=7.2, weight="bold", color=text_color, ha="center", va="center"
    )
    ax.text(
        50,
        3.8,
        "Correlation-aware camera and patient-level score changes",
        fontsize = 5.9,
        color    = neutral_color,
        ha       = "center",
        va       = "center",
    )

    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure1", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record("Figure1A", stem, "Independent paired-cohort myeloid workflow", source, pdf, png)
    return manifest_row


def figure1_clean_myeloid_umap(dataset: str, panel: str) -> dict[str, object]:
    """Draw annotated cells for one independently reconstructed cohort."""
    stem = f"{panel}_{dataset}_clean_myeloid_umap"
    atlas_root = legacy_project_root / "write/43_scgbm_style_clean_myeloid_reclustering"
    source_table = pd.read_csv(atlas_root / dataset / "clean_myeloid_umap_source.csv.gz")

    # Keep the reviewed state order, colors and cohort-specific annotation offsets.
    myeloid_state_order = [
        "Homeostatic microglia",
        "Transitional Mg/TAM",
        "Complement/HLA-II TAM",
        "Mg-inflammatory-high",
        "Inflammatory myeloid",
        "Scavenger/lipid TAM",
        "Hypoxia-glycolytic TAM",
        "IFN-responsive myeloid",
        "Heat-shock myeloid",
        "Cycling myeloid",
        "Monocyte-like",
        "cDC",
    ]
    myeloid_state_colors = {
        "Homeostatic microglia": "#159D88",
        "Transitional Mg/TAM": "#73B6B2",
        "Complement/HLA-II TAM": "#3B6FB6",
        "Mg-inflammatory-high": "#D8574F",
        "Inflammatory myeloid": "#D98270",
        "Scavenger/lipid TAM": "#D18A34",
        "Hypoxia-glycolytic TAM": "#7B61A8",
        "IFN-responsive myeloid": "#4C9A61",
        "Heat-shock myeloid": "#B06B88",
        "Cycling myeloid": "#E2B33B",
        "Monocyte-like": "#8C6D5A",
        "cDC": "#6B778D",
    }
    myeloid_label_offsets = {
        "GSE174554": {
            "Mg-inflammatory-high": (-0.035, 0.025),
            "Inflammatory myeloid": (-0.020, 0.000),
            "Transitional Mg/TAM": (0.060, 0.000),
            "IFN-responsive myeloid": (-0.020, -0.005),
            "Homeostatic microglia": (0.010, -0.015),
            "Complement/HLA-II TAM": (0.035, -0.005),
            "Scavenger/lipid TAM": (0.035, 0.000),
        },
        "GSE274546": {
            "Mg-inflammatory-high": (0.035, 0.000),
            "Inflammatory myeloid": (0.040, -0.010),
            "Transitional Mg/TAM": (-0.075, -0.010),
            "Complement/HLA-II TAM": (0.075, 0.000),
            "Scavenger/lipid TAM": (-0.040, -0.025),
            "Hypoxia-glycolytic TAM": (0.055, -0.010),
            "Heat-shock myeloid": (0.055, 0.010),
            "Monocyte-like": (0.000, -0.020),
            "cDC": (0.020, -0.010),
        },
    }
    unknown = sorted(set(source_table["curated_myeloid_subtype"]) - set(myeloid_state_order))
    if unknown:
        raise ValueError(f"Unexpected reviewed myeloid states for {dataset}: {unknown}")
    source = write_source(source_table, "Figure1", stem)

    fig = plt.figure(figsize=(mm(91), mm(88)), facecolor="white")
    ax = fig.add_axes([0.035, 0.205, 0.93, 0.64])

    # Plot the larger cell groups first and label the core of each state's largest cluster.
    counts = source_table["curated_myeloid_subtype"].value_counts()
    x_span = float(source_table["UMAP1"].max() - source_table["UMAP1"].min())
    y_span = float(source_table["UMAP2"].max() - source_table["UMAP2"].min())
    point_size = 0.62 if len(source_table) < 30_000 else 0.34
    for state in counts.sort_values(ascending=False).index:
        subset_table = source_table[source_table["curated_myeloid_subtype"].eq(state)]
        ax.scatter(
            subset_table["UMAP1"],
            subset_table["UMAP2"],
            s = point_size,
            color = myeloid_state_colors[state],
            alpha = 0.62,
            linewidths = 0,
            rasterized = True,
        )
    for state in myeloid_state_order:
        subset_table = source_table[source_table["curated_myeloid_subtype"].eq(state)]
        if subset_table.empty:
            continue
        largest_cluster = subset_table["leiden_myeloid"].astype(str).value_counts().index[0]
        label_cells_table = subset_table[subset_table["leiden_myeloid"].astype(str).eq(largest_cluster)]
        x, y = label_position(label_cells_table)
        dx, dy = myeloid_label_offsets.get(dataset, {}).get(state, (0.0, 0.0))
        text_x = x + dx * x_span
        text_y = y + dy * y_span
        short_label = {
            "Homeostatic microglia": "Homeostatic Mg",
            "Complement/HLA-II TAM": "Complement / HLA-II",
            "Mg-inflammatory-high": "Mg-inflammatory",
            "Inflammatory myeloid": "Inflammatory",
            "Scavenger/lipid TAM": "Scavenger / lipid",
            "Hypoxia-glycolytic TAM": "Hypoxia / glycolysis",
            "IFN-responsive myeloid": "IFN-responsive",
            "Heat-shock myeloid": "Heat-shock",
            "Cycling myeloid": "Cycling",
        }.get(state, state)
        if dx or dy:
            ax.plot(
                [x, text_x],
                [y, text_y],
                color  = myeloid_state_colors[state],
                lw     = 0.45,
                alpha  = 0.75,
                zorder = 4,
            )
        ax.text(
            text_x,
            text_y,
            short_label,
            fontsize = 5.2,
            weight   = "bold",
            color    = myeloid_state_colors[state],
            ha       = "center",
            va       = "center",
            path_effects = [pe.withStroke(linewidth=2.4, foreground="white")],
            zorder   = 5,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    handles = [
        mpl.lines.Line2D(
            [],
            [],
            marker     = "o",
            linestyle  = "none",
            markersize = 3.8,
            markerfacecolor = myeloid_state_colors[state],
            markeredgewidth = 0,
            label      = state,
        )
        for state in myeloid_state_order
        if state in counts.index
    ]
    fig.legend(
        handles        = handles,
        loc            = "lower center",
        bbox_to_anchor = (0.50, 0.025),
        ncol           = 3,
        frameon        = False,
        handletextpad  = 0.35,
        columnspacing  = 1.0,
        fontsize       = 5.25,
    )
    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(
        0.05,
        0.905,
        f"{len(source_table):,} clean myeloid cells  |  {source_table['leiden_myeloid'].nunique()} clusters",
        fontsize = 6.2,
        color    = neutral_color,
        va       = "top",
    )
    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure1", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record(
        panel, stem, f"{dataset} independently integrated and DEG-reviewed clean-myeloid UMAP", source, pdf, png
    )
    return manifest_row


def rank_density(rank_table: pd.DataFrame, bins: int = 75) -> np.ndarray:
    """Return log1p bin counts for the rank scatterplot background."""
    x = rank_table["rank_stat_GSE174554"].to_numpy(dtype=float)
    y = rank_table["rank_stat_GSE274546"].to_numpy(dtype=float)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)
    density = np.log1p(counts[xi, yi])
    return density


def figure1_rank_rank() -> dict[str, object]:
    """Draw the candidate paired-rank scatter and return the pane's manifest row."""
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    step41_root = legacy_project_root / "write/41_mg_inflammatory_sci_rebuild"
    curves_table = pd.read_csv(step41_root / "Figure1/source_data/raw20_gsea_curves.csv.gz")

    # Pair dataset-specific statistics and program membership by gene.
    cohort_tables = {}
    for dataset in datasets:
        cohort_table = curves_table.loc[curves_table["dataset"].eq(dataset), ["gene", "rank", "rank_stat", "hit"]]
        cohort_tables[dataset] = cohort_table.rename(
            columns={"rank": f"rank_{dataset}", "rank_stat": f"rank_stat_{dataset}", "hit": f"program_{dataset}"}
        )
    rank_table_total = cohort_tables["GSE174554"].merge(
        cohort_tables["GSE274546"], on="gene", how="inner", validate="one_to_one"
    )
    rank_table_total["miller_raw20"] = rank_table_total[["program_GSE174554", "program_GSE274546"]].any(axis=1)
    rank_table_total["density"] = rank_density(rank_table_total)
    program_table = rank_table_total[rank_table_total["miller_raw20"]].copy()
    program_table["label_score"] = (
        program_table["rank_stat_GSE174554"].abs() + program_table["rank_stat_GSE274546"].abs()
    )
    label_genes = set(program_table.nlargest(min(6, len(program_table)), "label_score")["gene"])
    rank_table_total["labelled"] = rank_table_total["gene"].isin(label_genes)
    source = write_source(rank_table_total, "Figure1", stem)

    fig = plt.figure(figsize=(mm(92), mm(83)), facecolor="white")
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
    background_table = rank_table_total[~rank_table_total["miller_raw20"]].sort_values("density")
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c     = background_table["density"],
        cmap  = "Greys",
        s     = 2.1,
        alpha = 0.70,
        linewidths = 0,
        rasterized = True,
    )
    ax.scatter(
        program_table["rank_stat_GSE174554"],
        program_table["rank_stat_GSE274546"],
        s         = 24,
        color     = "#D8574F",
        edgecolor = "white",
        linewidth = 0.65,
        zorder    = 4,
    )
    label_offsets = {
        "CCL3": (3, 7),
        "CCL4": (-26, 7),
        "CH25H": (5, 5),
        "SGK1": (-25, -10),
        "FOLR2": (5, 7),
        "BHLHE41": (5, -8),
    }
    for _, row in program_table[program_table["gene"].isin(label_genes)].iterrows():
        offset = label_offsets.get(row["gene"], (3, 3))
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext       = offset,
            textcoords   = "offset points",
            fontsize     = 5.7,
            color        = "#A63D38",
            path_effects = [pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", lw=0.55)
    ax.axvline(0, color="#A4A4A4", lw=0.55)
    limit = float(np.nanmax(np.abs(rank_table_total[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color=grid_color, lw=0.42, zorder=0)
    figure_style.clean_axis(ax)
    ax.tick_params(length=2.4)
    rho, _ = spearmanr(rank_table_total["rank_stat_GSE174554"], rank_table_total["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(rank_table_total):,}\nSpearman rho = {rho:.2f}",
        transform = ax.transAxes,
        va        = "top",
        fontsize  = 6.5,
        color     = neutral_color,
    )
    ax.text(
        0.98,
        0.03,
        f"MG-Inflam20: {len(program_table)} measured genes",
        transform = ax.transAxes,
        ha        = "right",
        va        = "bottom",
        fontsize  = 6.4,
        color     = "#A63D38",
    )
    fig.text(0.16, 0.965, "Cross-cohort all-gene concordance", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.16,
        0.915,
        "Transcriptome density with MG-Inflam20 genes highlighted",
        fontsize = 6.4,
        color    = neutral_color,
        va       = "top",
    )
    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure1", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record(
        "Figure1D", stem, "All-gene rank-stat density with MG-Inflam20 genes highlighted", source, pdf, png
    )
    return manifest_row


def figure2_raw20_umap(dataset: str, panel: str, vmin: float, vmax: float) -> dict[str, object]:
    """Draw an MG-Inflam20 UMAP using the shared display range."""
    stem = f"{panel}_{dataset}_raw20_terrain"
    step41_root = legacy_project_root / "write/41_mg_inflammatory_sci_rebuild"
    source_table = pd.read_csv(step41_root / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz")
    source = write_source(source_table, "Figure2", stem)
    ordered_table = source_table.sort_values("raw20_score")

    fig = plt.figure(figsize=(mm(88), mm(78)), facecolor="white")
    ax = fig.add_axes([0.045, 0.035, 0.80, 0.78])
    points = ax.scatter(
        ordered_table["UMAP1"],
        ordered_table["UMAP2"],
        c     = ordered_table["raw20_score"],
        s     = 0.55,
        cmap  = "magma",
        vmin  = vmin,
        vmax  = vmax,
        alpha = 0.94,
        linewidths = 0,
        rasterized = True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    cax = fig.add_axes([0.88, 0.20, 0.023, 0.48])
    cb = fig.colorbar(points, cax=cax)
    cb.ax.set_title("MG-Inflam20\nscore", fontsize=6.2, pad=4, loc="center")
    cb.ax.tick_params(labelsize=6.0, length=2)
    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(0.05, 0.905, "Continuous program terrain", fontsize=6.6, color=neutral_color, va="top")
    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure2", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record(panel, stem, f"{dataset} continuous MG-Inflam20 score terrain", source, pdf, png)
    return manifest_row


def ridge_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Return KDE heights with a unit peak, or a zero curve for constant or sparse data."""
    values = values[np.isfinite(values)]
    if len(values) < 2 or np.isclose(np.std(values), 0):
        density = np.zeros_like(grid)
    else:
        estimate = gaussian_kde(values)
        density = estimate(grid)
        peak = density.max()
        if peak > 0:
            density = density / peak
    return density


def figure2_ridgeline(dataset: str, panel: str, x_limits: tuple[float, float]) -> dict[str, object]:
    """Draw distributions of patient MG-Inflam20 scores by unsupervised cluster."""
    stem = f"{panel}_{dataset}_cluster_ridgeline"
    step41_root = legacy_project_root / "write/41_mg_inflammatory_sci_rebuild"
    interval_table = pd.read_csv(step41_root / "Figure2/source_data/fig2_patient_cluster_interval_source.csv")
    source_table = interval_table[
        interval_table["dataset"].eq(dataset) & interval_table["eligible"].astype(bool)
    ].copy()
    label_table = pd.read_csv(
        step41_root / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz",
        usecols=["paired_cluster", "cluster_label"],
    )
    label_table = label_table.dropna().drop_duplicates()
    if dataset == "GSE274546":

        # Use the Step 36 reviewed descriptive annotation for the GSE274546 cluster names.
        annotation_root = legacy_project_root / "write/36_gse274546_raw_independent_reannotation"
        annotation_path = (
            annotation_root / "04_paired_patient_harmony_clusters/paired_cluster_descriptive_annotation.csv"
        )
        annotation_table = pd.read_csv(annotation_path)
        label_table = annotation_table[["cluster", "descriptive_state"]].rename(
            columns={"cluster": "paired_cluster", "descriptive_state": "cluster_label"}
        )
        labels = []
        for cluster, label in label_table.itertuples(index=False, name=None):
            labels.append(f"C{int(cluster)} {label}")
        label_table["cluster_label"] = labels
    label_table["paired_cluster"] = pd.to_numeric(label_table["paired_cluster"]).astype(int).astype(str)
    label_map = label_table.set_index("paired_cluster")["cluster_label"].to_dict()
    source_table["paired_cluster"] = pd.to_numeric(source_table["paired_cluster"]).astype(int).astype(str)
    source_table["cluster_label_display"] = source_table["paired_cluster"].map(label_map)
    if source_table["cluster_label_display"].isna().any():
        raise ValueError(f"Missing full cluster labels for {dataset}")
    cluster_summary_table = (
        source_table.groupby(["paired_cluster", "cluster_label_display"], observed=True)["pooled_patient_raw20"]
        .agg(["median", "count"])
        .reset_index()
        .sort_values("median")
    )
    order = cluster_summary_table["paired_cluster"].astype(str).tolist()
    source_table["ridge_order"] = source_table["paired_cluster"].map({cluster: i for i, cluster in enumerate(order)})
    source = write_source(source_table, "Figure2", stem)

    fig = plt.figure(figsize=(mm(100), mm(91)), facecolor="white")
    ax = fig.add_axes([0.42, 0.13, 0.55, 0.70])
    grid = np.linspace(x_limits[0], x_limits[1], 300)
    color = dataset_colors[dataset]
    labels = []
    for i, cluster in enumerate(order):
        subset_table = source_table[source_table["paired_cluster"].eq(cluster)]
        values = subset_table["pooled_patient_raw20"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        baseline = float(i)
        ax.fill_between(grid, baseline, baseline + density * 0.78, color=color, alpha=0.36, linewidth=0)
        ax.plot(grid, baseline + density * 0.78, color=color, lw=0.9)
        jitter = np.random.default_rng(jitter_seed + i).uniform(-0.055, 0.055, len(values))
        ax.scatter(values, baseline - 0.10 + jitter, s=5.5, color=text_color, alpha=0.55, linewidths=0, zorder=3)
        ax.scatter(
            np.median(values),
            baseline - 0.10,
            marker    = "D",
            s         = 15,
            color     = color,
            edgecolor = "white",
            linewidth = 0.5,
            zorder    = 4,
        )
        labels.append(str(subset_table["cluster_label_display"].iloc[0]))
        ax.text(
            0.985,
            baseline + 0.49,
            f"n={len(values)}",
            transform = ax.get_yaxis_transform(),
            ha        = "right",
            va        = "center",
            fontsize  = 5.5,
            color     = neutral_color,
        )
    ax.set_yticks(range(len(order)), labels)
    ax.set_ylim(-0.35, len(order) - 0.02 + 0.82)
    ax.set_xlim(*x_limits)
    ax.set_xlabel("Patient-level pooled MG-Inflam20 score")
    ax.grid(axis="x", color=grid_color, lw=0.45)
    ax.tick_params(axis="y", length=0, labelsize=5.2, pad=2)
    figure_style.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(length=2.4)
    fig.text(0.08, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=color)
    fig.text(0.08, 0.913, "MG-Inflam20 score across unsupervised clusters", fontsize=6.5, color=neutral_color, va="top")
    fig.text(0.08, 0.87, "Each dot is one eligible patient-cluster unit", fontsize=6.2, color=neutral_color, va="top")
    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure2", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record(panel, stem, f"{dataset} patient-level cluster ridgeline", source, pdf, png)
    return manifest_row


def figure2_delta_raincloud() -> dict[str, object]:
    """Draw the original patient contrasts for the published atlas states."""
    stem = "Figure2E_GSE278456_patient_delta_raincloud"
    step41_root = legacy_project_root / "write/41_mg_inflammatory_sci_rebuild"
    source_table = pd.read_csv(step41_root / "Figure2/source_data/fig2_gse278456_all_state_patient_deltas.csv")

    # The candidate panel displays the original six named states and palette.
    state_order = ["MCG1", "MCG2", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]
    state_colors = {
        "MCG1": "#D8574F",
        "MCG2": "#159D88",
        "MAC1": "#277DA1",
        "MAC2": "#72B7B2",
        "E-MDSC": "#E09F3E",
        "M-MDSC": "#8C6D5A",
    }
    source_table = source_table[source_table["state"].isin(state_order)].copy()
    source_table["state_order"] = source_table["state"].map({state: i for i, state in enumerate(state_order)})
    source = write_source(source_table, "Figure2", stem)

    fig = plt.figure(figsize=(mm(118), mm(82)), facecolor="white")
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.67])
    x_min = min(-0.05, float(source_table["delta"].min()) - 0.06)
    x_max = float(source_table["delta"].max()) + 0.08
    grid = np.linspace(x_min, x_max, 320)
    for i, state in enumerate(state_order):
        values = source_table.loc[source_table["state"].eq(state), "delta"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        color = state_colors[state]
        ax.fill_between(grid, i, i + density * 0.34, color=color, alpha=0.40, linewidth=0)
        ax.plot(grid, i + density * 0.34, color=color, lw=0.9)
        jitter = np.random.default_rng(jitter_seed + 100 + i).uniform(-0.30, -0.08, len(values))
        ax.scatter(values, i + jitter, s=14, color=color, edgecolor="white", linewidth=0.35, alpha=0.80, zorder=3)
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        ax.hlines(i - 0.02, q1, q3, color=text_color, lw=2.0, zorder=4)
        ax.scatter(median, i - 0.02, s=20, marker="D", color=text_color, edgecolor="white", linewidth=0.45, zorder=5)
        ax.text(x_max, i, f"n={len(values)}", ha="right", va="center", fontsize=6.0, color=neutral_color)
    ax.axvline(0, color="#929292", lw=0.7, ls="--")
    ax.set_yticks(range(len(state_order)), state_order)
    ax.set_ylim(-0.45, len(state_order) - 0.45)
    ax.invert_yaxis()
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Within-patient MG-Inflam20 score delta (state mean - other states mean)")
    ax.grid(axis="x", color=grid_color, lw=0.45)
    ax.tick_params(axis="y", length=0)
    figure_style.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(length=2.4)
    fig.text(0.17, 0.965, "GSE278456 patient-level state contrasts", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.17,
        0.91,
        "Cloud: density | dots: patients | diamond and line: median and IQR",
        fontsize = 6.4,
        color    = neutral_color,
        va       = "top",
    )
    pdf, png = figure_style.save_panel(fig, FIG_ROOT / "Figure2", stem, "Step42 recurrent GBM visual story")
    manifest_row = panel_record(
        "Figure2E", stem, "GSE278456 patient delta raincloud across published states", source, pdf, png
    )
    return manifest_row


def write_manifest(figure_name: str, records: list[dict[str, object]]) -> None:
    """Write the per-figure source and panel manifest. Returns nothing."""
    out = WRITE_ROOT / figure_name
    out.mkdir(parents=True, exist_ok=True)
    manifest_table = pd.DataFrame(records)
    manifest_path = out / f"{figure_name}_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)
    if not manifest_path.stat().st_size:
        raise ValueError(f"Empty panel manifest for {figure_name}")


def main() -> None:
    """Draw the two candidate figures and write their source manifests. Returns nothing."""
    np.random.seed(jitter_seed)
    configure_style()
    for figure_name in ["Figure1", "Figure2"]:
        (FIG_ROOT / figure_name).mkdir(parents=True, exist_ok=True)
        (WRITE_ROOT / figure_name / "source_data").mkdir(parents=True, exist_ok=True)

    figure1_records = [
        figure1_design(),
        figure1_clean_myeloid_umap("GSE174554", "Figure1B"),
        figure1_clean_myeloid_umap("GSE274546", "Figure1C"),
        figure1_rank_rank(),
    ]
    write_manifest("Figure1", figure1_records)

    # The UMAPs share their score limits, and the ridgelines share the same patient score axis.
    step41_root = legacy_project_root / "write/41_mg_inflammatory_sci_rebuild"
    raw20_tables = {
        dataset: pd.read_csv(step41_root / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz")
        for dataset in datasets
    }
    all_scores = np.concatenate(
        [source_table["raw20_score"].to_numpy(dtype=float) for source_table in raw20_tables.values()]
    )
    vmin, vmax = np.nanpercentile(all_scores, [2, 98])
    cluster_source_table = pd.read_csv(step41_root / "Figure2/source_data/fig2_patient_cluster_interval_source.csv")
    eligible_scores = cluster_source_table.loc[
        cluster_source_table["eligible"].astype(bool), "pooled_patient_raw20"
    ].to_numpy(dtype=float)
    x_pad = max(0.03, (np.nanmax(eligible_scores) - np.nanmin(eligible_scores)) * 0.04)
    x_limits = (float(np.nanmin(eligible_scores) - x_pad), float(np.nanmax(eligible_scores) + x_pad))
    figure2_records = [
        figure2_raw20_umap("GSE174554", "Figure2A", float(vmin), float(vmax)),
        figure2_raw20_umap("GSE274546", "Figure2B", float(vmin), float(vmax)),
        figure2_ridgeline("GSE174554", "Figure2C", x_limits),
        figure2_ridgeline("GSE274546", "Figure2D", x_limits),
        figure2_delta_raincloud(),
    ]
    write_manifest("Figure2", figure2_records)
    print(f"Candidate panels: {len(figure1_records)} for Figure 1 and {len(figure2_records)} for Figure 2")


if __name__ == "__main__":
    main()
