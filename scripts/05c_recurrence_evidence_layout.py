#!/usr/bin/env python3
# Purpose: Draw the paired-cohort design and provide the plotting functions used for Figures 1 and 2.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05c_recurrence_evidence_layout.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
mm = figure_style.mm
sha256 = figure_style.sha256
save_panel = figure_style.save_panel
rank_density = figure_style.rank_density


def configure_style() -> None:
    """Set this step's Matplotlib fonts, colors and line widths. No arguments or return value."""
    text = "#252525"

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
            "axes.edgecolor": text,
            "axes.labelcolor": text,
            "text.color": text,
            "xtick.color": text,
            "ytick.color": text,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


root = Path(__file__).resolve().parents[1]

# Use the published output directory and basename.
write_root = root / "write/42_visual_story_rebuild"


def write_source(source_table: pd.DataFrame, figure_name: str, stem: str) -> Path:
    """Write the source table, check that the CSV is nonempty and return its path."""

    # Use the published output directory and basename.
    out = write_root / figure_name / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}_source.csv"
    csv_path = path
    source_table.to_csv(csv_path, index=False)

    # Use the existing input checks and branch conditions.
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    return path


def panel_record(
    panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path
) -> dict[str, object]:
    """Return the manifest record for a panel, its description and the source, PDF and PNG paths."""
    result = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source_csv": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": sha256(source),
        "pdf_sha256": sha256(pdf),
        "png_sha256": sha256(png),
        "png_dpi": 600,
    }
    return result


# Step 51 overrides this output root after importing this module.
fig_root = root / "figures/42_visual_story_rebuild"


def figure1_design() -> dict[str, object]:
    """Build the design from the Step 38 threshold-20 audit and return the ten-field panel record."""

    # Read the formal paired audit and apply the published threshold of 20 and cohort count checks.
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"

    audit_table = pd.read_csv(step38 / "independent_input_pair_audit.csv")
    datasets = ["GSE174554", "GSE274546"]

    source_table = audit_table.loc[
        audit_table["dataset"].isin(datasets) & audit_table["threshold"].eq(20)
    ].copy()
    if source_table.groupby("dataset").size().to_dict() != {"GSE174554": 1, "GSE274546": 1}:
        raise ValueError("Step38 must contain one formal threshold-20 audit row per cohort.")
    source_table["input_libraries_or_matrices"] = source_table["dataset"].map(
        {"GSE174554": 91, "GSE274546": 111}
    )

    # Organize the design table in the published column and cohort order.
    source_table["clean_myeloid"] = source_table["n_cells_primary"].astype(int) + source_table[
        "n_cells_recurrent"
    ].astype(int)
    source_table["formal_pairs_threshold20"] = source_table["n_pairs"].astype(int)
    observed = (
        source_table.set_index("dataset")[["formal_pairs_threshold20", "clean_myeloid"]]
        .astype(int)
        .to_dict("index")
    )
    expected = {
        "GSE174554": {"formal_pairs_threshold20": 17, "clean_myeloid": 12489},
        "GSE274546": {"formal_pairs_threshold20": 45, "clean_myeloid": 55568},
    }
    if observed != expected:
        raise ValueError(f"Unexpected IDH-restricted Figure1A analysis set: {observed}")

    # Organize the design table in the published column and cohort order.
    source_table = source_table.sort_values("dataset").reset_index(drop=True)
    source_table["modality"] = "snRNA-seq"
    source_table["raw_input_unit"] = ["10x matrices", "10x libraries"]
    source_table["minimum_cells_per_condition"] = 20
    source_table["processing_scope"] = "Each cohort reconstructed, annotated and tested independently"

    # Organize the design table in the published column and cohort order.
    source_table["statistical_testing"] = (
        "Patient pseudobulk, paired edgeR and Miller-IM GSEA performed separately within each cohort"
    )
    source_table["replication_assessment"] = "Cohort results compared only after independent testing"
    source_table["workflow"] = (
        "all-cell reconstruction -> pan-myeloid re-clustering -> marker an"
        "d cluster-DEG review -> patient raw-count pseudobulk -> paired ed"
        "geR -> prespecified Miller-IM program GSEA"
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

    stem = "Figure1A_compact_paired_design"
    source = write_source(source_table, "Figure1", stem)
    fig = plt.figure(figsize=(mm(183), mm(94)), facecolor="white")
    ax = fig.add_axes([0.025, 0.025, 0.95, 0.95])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.axis("off")

    def rounded_box(
        x0: float, x1: float, y0: float, y1: float, edge: str, fill: str, lw: float = 0.8
    ) -> None:
        """Draw the box using the supplied bounds, colors and line width. No return value."""
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                boxstyle="round,pad=0.02,rounding_size=0.7",
                facecolor=fill,
                edgecolor=edge,
                linewidth=lw,
            )
        )

    cohort_boxes = [(4.0, 47.5), (52.5, 96.0)]
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    neutral = "#717171"

    text = "#252525"

    for row, (x0, x1) in zip(source_table.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        rounded_box(x0, x1, 80.0, 97.0, color, pale, lw=1.0)
        ax.text(
            center,
            94.5,
            f"{row['dataset']}  |  snRNA-seq",
            fontsize=8.3,
            weight="bold",
            color=color,
            ha="center",
            va="top",
        )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(
            center,
            89.3,
            f"{int(row['input_libraries_or_matrices'])} {row['raw_input_unit']}",
            fontsize=6.1,
            color=neutral,
            ha="center",
            va="top",
        )
        primary_x, recurrent_x = (center - 8.0, center + 8.0)
        rounded_box(primary_x - 4.9, primary_x + 4.9, 81.7, 86.7, color, "white", lw=0.9)
        rounded_box(recurrent_x - 5.5, recurrent_x + 5.5, 81.7, 86.7, color, color, lw=0.9)
        ax.text(primary_x, 84.2, "Primary", fontsize=6.3, color=text, ha="center", va="center")

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(recurrent_x, 84.2, "1st recurrence", fontsize=6.3, color="white", ha="center", va="center")
        ax.annotate(
            "",
            xy=(recurrent_x - 5.8, 84.2),
            xytext=(primary_x + 5.2, 84.2),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0),
        )
    for center in [sum(box) / 2 for box in cohort_boxes]:
        ax.annotate(
            "",
            xy=(center, 74.8),
            xytext=(center, 79.5),
            arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85),
        )

    # Build the design panel from the audit rows, using the existing cell and patient-pair counts.
    workflow_boxes = [
        (
            8.0,
            92.0,
            64.0,
            74.0,
            "Independent all-cell reconstruction",
            "Library QC · doublet removal · basic lineage annotation",
        ),
        (
            8.0,
            92.0,
            50.0,
            60.0,
            "Cohort-wise pan-myeloid re-clustering",
            "2,000 HVGs · PCA · neighbor graph · Leiden",
        ),
        (
            6.0,
            94.0,
            36.0,
            46.0,
            "Marker and cluster-DEG review",
            "Remove lymphoid/tumour contaminants and low-quality cells",
        ),
    ]
    for index, (x0, x1, y0, y1, title, subtitle) in enumerate(workflow_boxes):
        rounded_box(x0, x1, y0, y1, "#858585", "white", lw=0.75)
        ax.text(50, y1 - 2.2, title, fontsize=7.3, weight="bold", ha="center", va="top")
        ax.text(50, y0 + 2.2, subtitle, fontsize=6.1, color=neutral, ha="center", va="bottom")
        if index < len(workflow_boxes) - 1:
            next_y1 = workflow_boxes[index + 1][3]
            ax.annotate(
                "",
                xy=(50, next_y1 + 0.5),
                xytext=(50, y0 - 0.5),
                arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85),
            )
    analysis_y0, analysis_y1 = (14.5, 29.5)
    for row, (x0, x1) in zip(source_table.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        source_x = 35.0 if center < 50 else 65.0
        ax.annotate(
            "",
            xy=(center, analysis_y1 + 0.6),
            xytext=(source_x, 35.5),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=0.9,
                connectionstyle="arc3,rad=0.08" if center < 50 else "arc3,rad=-0.08",
            ),
        )

        rounded_box(x0, x1, analysis_y0, analysis_y1, color, pale, lw=1.0)
        ax.text(
            center,
            26.8,
            f"{row['dataset']}  |  {int(row['formal_pairs_threshold20'])} paired patients",
            fontsize=7.1,
            weight="bold",
            color=color,
            ha="center",
            va="center",
        )
        ax.text(
            center,
            22.0,
            f"{int(row['clean_myeloid']):,} analyzed myeloid cells",
            fontsize=5.9,
            color=neutral,
            ha="center",
            va="center",
        )
        ax.text(
            center,
            17.2,
            "Raw-count pseudobulk · paired edgeR · Miller-IM GSEA",
            fontsize=5.9,
            color=text,
            ha="center",
            va="center",
        )
    ax.annotate(
        "",
        xy=(46.8, 10.8),
        xytext=(sum(cohort_boxes[0]) / 2, analysis_y0 - 0.5),
        arrowprops=dict(
            arrowstyle="-|>", color=dataset_colors["GSE174554"], lw=0.9, connectionstyle="arc3,rad=-0.08"
        ),
    )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.annotate(
        "",
        xy=(53.2, 10.8),
        xytext=(sum(cohort_boxes[1]) / 2, analysis_y0 - 0.5),
        arrowprops=dict(
            arrowstyle="-|>", color=dataset_colors["GSE274546"], lw=0.9, connectionstyle="arc3,rad=0.08"
        ),
    )
    rounded_box(22.0, 78.0, 1.0, 11.8, "#159D88", mpl.colors.to_rgba("#159D88", 0.085), lw=1.0)
    ax.text(
        50,
        8.8,
        "Cross-cohort replication assessment",
        fontsize=7.2,
        weight="bold",
        color=text,
        ha="center",
        va="center",
    )
    ax.text(
        50,
        3.8,
        "Compare effect direction and nominal P after independent cohort tests",
        fontsize=5.9,
        color=neutral,
        ha="center",
        va="center",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")

    # Return the result using the existing field names.
    result = panel_record("Figure1A", stem, "Independent paired-cohort myeloid workflow", source, pdf, png)
    return result


def label_position(coordinate_table: pd.DataFrame) -> tuple[float, float]:
    """Return the central label coordinates using the original quantile cutoff."""

    points = coordinate_table[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
    center = np.nanmedian(points, axis=0)
    distance = np.sqrt(((points - center) ** 2).sum(axis=1))
    core = points[distance <= np.nanquantile(distance, 0.6)]
    target = np.nanmedian(core, axis=0) if len(core) else center

    # Return the result using the existing field names.
    result = (float(target[0]), float(target[1]))
    return result


def figure1_clean_myeloid_umap(dataset: str, panel: str) -> dict[str, object]:
    """Save the specified cohort's myeloid coordinates and labels as a panel and return its manifest record."""

    # Use the panel's published stem for its source table and figure files.
    stem = f"{panel}_{dataset}_clean_myeloid_umap"
    myeloid_atlas_sources = {
        "GSE174554": root
        / "write/43_scgbm_style_clean_myeloid_reclustering/GSE174554/clean_myeloid_umap_source.csv.gz",
        "GSE274546": root
        / "write/43_scgbm_style_clean_myeloid_reclustering/GSE274546/clean_myeloid_umap_source.csv.gz",
    }

    umap_table = pd.read_csv(myeloid_atlas_sources[dataset])
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

    unknown = sorted(set(umap_table["curated_myeloid_subtype"]) - set(myeloid_state_order))
    if unknown:
        raise ValueError(f"Unexpected reviewed myeloid states for {dataset}: {unknown}")
    source = write_source(umap_table, "Figure1", stem)

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig = plt.figure(figsize=(mm(91), mm(88)), facecolor="white")
    ax = fig.add_axes([0.035, 0.205, 0.93, 0.64])
    counts = umap_table["curated_myeloid_subtype"].value_counts()
    x_span = float(umap_table["UMAP1"].max() - umap_table["UMAP1"].min())
    y_span = float(umap_table["UMAP2"].max() - umap_table["UMAP2"].min())

    point_size = 0.62 if len(umap_table) < 30000 else 0.34
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

    for state in counts.sort_values(ascending=False).index:
        current_table = umap_table[umap_table["curated_myeloid_subtype"].eq(state)]
        ax.scatter(
            current_table["UMAP1"],
            current_table["UMAP2"],
            s=point_size,
            color=myeloid_state_colors[state],
            alpha=0.62,
            linewidths=0,
            rasterized=True,
        )
    myeloid_label_offsets = {
        "GSE174554": {
            "Mg-inflammatory-high": (-0.035, 0.025),
            "Inflammatory myeloid": (-0.02, 0.0),
            "Transitional Mg/TAM": (0.06, 0.0),
            "IFN-responsive myeloid": (-0.02, -0.005),
            "Homeostatic microglia": (0.01, -0.015),
            "Complement/HLA-II TAM": (0.035, -0.005),
            "Scavenger/lipid TAM": (0.035, 0.0),
        },
        "GSE274546": {
            "Mg-inflammatory-high": (0.035, 0.0),
            "Inflammatory myeloid": (0.04, -0.01),
            "Transitional Mg/TAM": (-0.075, -0.01),
            "Complement/HLA-II TAM": (0.075, 0.0),
            "Scavenger/lipid TAM": (-0.04, -0.025),
            "Hypoxia-glycolytic TAM": (0.055, -0.01),
            "Heat-shock myeloid": (0.055, 0.01),
            "Monocyte-like": (0.0, -0.02),
            "cDC": (0.02, -0.01),
        },
    }

    for state in myeloid_state_order:
        current_table = umap_table[umap_table["curated_myeloid_subtype"].eq(state)]
        if current_table.empty:
            continue
        largest_cluster = current_table["leiden_myeloid"].astype(str).value_counts().index[0]
        label_cell_table = current_table[current_table["leiden_myeloid"].astype(str).eq(largest_cluster)]
        x, y = label_position(label_cell_table)

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
                [x, text_x], [y, text_y], color=myeloid_state_colors[state], lw=0.45, alpha=0.75, zorder=4
            )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(
            text_x,
            text_y,
            short_label,
            fontsize=5.2,
            weight="bold",
            color=myeloid_state_colors[state],
            ha="center",
            va="center",
            path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
            zorder=5,
        )
    ax.set_xticks([])
    ax.set_yticks([])

    # Process cohorts, genes or panels in the published order.
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    handles = [
        mpl.lines.Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=3.8,
            markerfacecolor=myeloid_state_colors[state],
            markeredgewidth=0,
            label=state,
        )
        for state in myeloid_state_order
        if state in counts.index
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=3,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.0,
        fontsize=5.25,
    )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(
        0.05,
        0.905,
        f"{len(umap_table):,} clean myeloid cells  |  {umap_table['leiden_myeloid'].nunique()} clusters",
        fontsize=6.2,
        color="#717171",
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")
    result = panel_record(
        panel,
        stem,
        f"{dataset} independently integrated and DEG-reviewed clean-myeloid UMAP",
        source,
        pdf,
        png,
    )
    return result


# The published panel uses a tick length of 2.4. The shared axis helper uses 2.5.


def clean_axis(ax: mpl.axes.Axes, *, left: bool = True, bottom: bool = True) -> None:
    """Set the given axes' border style in place without changing data. No return value."""

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(width=0.6, length=2.4)


def figure1_rank_rank() -> dict[str, object]:
    """Save merged Step 41 gene ranks and the density plot, then return the manifest record."""

    # Read the two cohort rank tables and merge one to one by gene.
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    curve_table = pd.read_csv(step41 / "Figure1/source_data/raw20_gsea_curves.csv.gz")
    left_table = curve_table[curve_table["dataset"].eq("GSE174554")][
        ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE174554", "rank_stat": "rank_stat_GSE174554", "hit": "program_GSE174554"}
    )
    right_table = curve_table[curve_table["dataset"].eq("GSE274546")][
        ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE274546", "rank_stat": "rank_stat_GSE274546", "hit": "program_GSE274546"}
    )
    rank_table_total = left_table.merge(right_table, on="gene", how="inner", validate="one_to_one")

    # Prepare the cohort gene ranks using the original order and merge keys.
    rank_table_total["miller_raw20"] = rank_table_total[["program_GSE174554", "program_GSE274546"]].any(
        axis=1
    )
    rank_table_total["density"] = rank_density(rank_table_total)
    program_table = rank_table_total[rank_table_total["miller_raw20"]].copy()
    program_table["label_score"] = (
        program_table["rank_stat_GSE174554"].abs() + program_table["rank_stat_GSE274546"].abs()
    )
    label_genes = set(program_table.nlargest(min(6, len(program_table)), "label_score")["gene"])

    # Prepare the cohort gene ranks using the original order and merge keys.
    rank_table_total["labelled"] = rank_table_total["gene"].isin(label_genes)
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    source = write_source(rank_table_total, "Figure1", stem)
    fig = plt.figure(figsize=(mm(92), mm(83)), facecolor="white")
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
    background_table = rank_table_total[~rank_table_total["miller_raw20"]].sort_values("density")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c=background_table["density"],
        cmap="Greys",
        s=2.1,
        alpha=0.7,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        program_table["rank_stat_GSE174554"],
        program_table["rank_stat_GSE274546"],
        s=24,
        color="#D8574F",
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
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
            xytext=offset,
            textcoords="offset points",
            fontsize=5.7,
            color="#A63D38",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", lw=0.55)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.axvline(0, color="#A4A4A4", lw=0.55)
    limit = (
        float(np.nanmax(np.abs(rank_table_total[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy())))
        * 1.04
    )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E6E6E6", lw=0.42, zorder=0)
    clean_axis(ax)
    rho, _ = spearmanr(rank_table_total["rank_stat_GSE174554"], rank_table_total["rank_stat_GSE274546"])

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    neutral = "#717171"

    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(rank_table_total):,}\nSpearman rho = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.5,
        color=neutral,
    )
    ax.text(
        0.98,
        0.03,
        f"Miller-IM: {len(program_table)} measured genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        color="#A63D38",
    )
    fig.text(0.16, 0.965, "Cross-cohort all-gene concordance", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.16,
        0.915,
        "Transcriptome density with Miller-IM genes highlighted",
        fontsize=6.4,
        color=neutral,
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")

    # Return the result using the existing field names.
    result = panel_record(
        "Figure1D", stem, "All-gene rank-stat density with Miller-IM genes highlighted", source, pdf, png
    )
    return result


def figure2_raw20_umap(dataset: str, panel: str, vmin: float, vmax: float) -> dict[str, object]:
    """Save the cohort raw20 UMAP and source table at the given color limits and return the manifest record."""

    # Use the panel's published stem for its source table and figure files.
    stem = f"{panel}_{dataset}_raw20_terrain"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    umap_table = pd.read_csv(step41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz")
    source = write_source(umap_table, "Figure2", stem)
    ordered_table = umap_table.sort_values("raw20_score")
    fig = plt.figure(figsize=(mm(88), mm(78)), facecolor="white")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax = fig.add_axes([0.045, 0.035, 0.8, 0.78])
    points = ax.scatter(
        ordered_table["UMAP1"],
        ordered_table["UMAP2"],
        c=ordered_table["raw20_score"],
        s=0.55,
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        alpha=0.94,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    cax = fig.add_axes([0.88, 0.2, 0.023, 0.48])
    cb = fig.colorbar(points, cax=cax)
    cb.ax.set_title("Miller-IM\nscore", fontsize=6.2, pad=4, loc="center")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    cb.ax.tick_params(labelsize=6.0, length=2)
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(0.05, 0.905, "Continuous program terrain", fontsize=6.6, color="#717171", va="top")
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(panel, stem, f"{dataset} continuous Miller-IM score terrain", source, pdf, png)

    return result


def ridge_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Return peak-scaled KDE values on the grid, or zeros for fewer than two values or zero variance."""

    values = values[np.isfinite(values)]
    if len(values) < 2 or np.isclose(np.std(values), 0):
        result = np.zeros_like(grid)
        return result
    density = gaussian_kde(values)(grid)
    peak = density.max()
    result = density / peak if peak > 0 else density

    return result


def figure2_ridgeline(dataset: str, panel: str, x_limits: tuple[float, float]) -> dict[str, object]:
    """Save the cohort score ridgeline and source table at the specified x limits and return the panel record."""

    # Use the panel's published stem for its source table and figure files.
    stem = f"{panel}_{dataset}_cluster_ridgeline"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    raw_table = pd.read_csv(step41 / "Figure2/source_data/fig2_patient_cluster_interval_source.csv")
    patient_table = raw_table[raw_table["dataset"].eq(dataset) & raw_table["eligible"].astype(bool)].copy()
    label_table = (
        pd.read_csv(
            step41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz",
            usecols=["paired_cluster", "cluster_label"],
        )
        .dropna()
        .drop_duplicates()
    )
    if dataset == "GSE274546":

        # Read the frozen inputs with their original row order, dataset names and column names.
        annotation_table = pd.read_csv(
            root
            / (
                "write/36_gse274546_raw_independent_reannotation/04_paired_patient"
                "_harmony_clusters/paired_cluster_descriptive_annotation.csv"
            )
        )
        label_table = annotation_table[["cluster", "descriptive_state"]].rename(
            columns={"cluster": "paired_cluster", "descriptive_state": "cluster_label"}
        )
        label_table["cluster_label"] = label_table.apply(
            lambda row: f"C{int(row['paired_cluster'])} {row['cluster_label']}", axis=1
        )

    # Apply the original labels and placement rules.
    label_table["paired_cluster"] = pd.to_numeric(label_table["paired_cluster"]).astype(int).astype(str)
    label_map = label_table.set_index("paired_cluster")["cluster_label"].to_dict()
    patient_table["paired_cluster"] = pd.to_numeric(patient_table["paired_cluster"]).astype(int).astype(str)
    patient_table["cluster_label_display"] = patient_table["paired_cluster"].map(label_map)
    if patient_table["cluster_label_display"].isna().any():
        raise ValueError(f"Missing full cluster labels for {dataset}")

    # Compute the cohort record count and extrema using the original output columns.
    summary_table = (
        patient_table.groupby(["paired_cluster", "cluster_label_display"], observed=True)[
            "pooled_patient_raw20"
        ]
        .agg(["median", "count"])
        .reset_index()
        .sort_values("median")
    )
    order = summary_table["paired_cluster"].astype(str).tolist()
    patient_table["ridge_order"] = patient_table["paired_cluster"].map(
        {cluster: i for i, cluster in enumerate(order)}
    )
    source = write_source(patient_table, "Figure2", stem)
    fig = plt.figure(figsize=(mm(100), mm(91)), facecolor="white")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax = fig.add_axes([0.42, 0.13, 0.55, 0.7])
    grid = np.linspace(x_limits[0], x_limits[1], 300)
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    color = dataset_colors[dataset]
    tick_labels = []
    neutral = "#717171"

    for i, cluster in enumerate(order):
        current_table = patient_table[patient_table["paired_cluster"].eq(cluster)]
        values = current_table["pooled_patient_raw20"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        baseline = float(i)
        ax.fill_between(grid, baseline, baseline + density * 0.78, color=color, alpha=0.36, linewidth=0)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.plot(grid, baseline + density * 0.78, color=color, lw=0.9)
        jitter = np.random.default_rng(20260713 + i).uniform(-0.055, 0.055, len(values))
        ax.scatter(values, baseline - 0.1 + jitter, s=5.5, color="#252525", alpha=0.55, linewidths=0, zorder=3)
        ax.scatter(
            np.median(values),
            baseline - 0.1,
            marker="D",
            s=15,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        tick_labels.append(str(current_table["cluster_label_display"].iloc[0]))

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(
            0.985,
            baseline + 0.49,
            f"n={len(values)}",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=5.5,
            color=neutral,
        )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_yticks(range(len(order)), tick_labels)
    ax.set_ylim(-0.35, len(order) - 0.02 + 0.82)
    ax.set_xlim(*x_limits)
    ax.set_xlabel("Patient-level pooled Miller-IM score")
    ax.grid(axis="x", color="#E6E6E6", lw=0.45)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.tick_params(axis="y", length=0, labelsize=5.2, pad=2)
    clean_axis(ax, left=False, bottom=True)
    fig.text(0.08, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=color)
    fig.text(
        0.08, 0.913, "Miller-IM score across unsupervised clusters", fontsize=6.5, color=neutral, va="top"
    )
    fig.text(
        0.08, 0.87, "Each dot is one eligible patient-cluster unit", fontsize=6.2, color=neutral, va="top"
    )

    # Save the PDF master at the published size and render its PNG from that PDF.
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(panel, stem, f"{dataset} patient-level cluster ridgeline", source, pdf, png)
    return result


def figure2_delta_raincloud() -> dict[str, object]:
    """Save the existing GSE278456 paired state changes as a raincloud panel and return the panel record."""

    # Use the panel's published stem for its source table and figure files.
    stem = "Figure2E_GSE278456_patient_delta_raincloud"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    delta_table = pd.read_csv(step41 / "Figure2/source_data/fig2_gse278456_all_state_patient_deltas.csv")
    state_order = ["MCG1", "MCG2", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]

    delta_table = delta_table[delta_table["state"].isin(state_order)].copy()
    delta_table["state_order"] = delta_table["state"].map({state: i for i, state in enumerate(state_order)})
    source = write_source(delta_table, "Figure2", stem)

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig = plt.figure(figsize=(mm(118), mm(82)), facecolor="white")
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.67])
    x_min = min(-0.05, float(delta_table["delta"].min()) - 0.06)
    x_max = float(delta_table["delta"].max()) + 0.08
    grid = np.linspace(x_min, x_max, 320)

    # Keep the seed at this point to reproduce the sampling or jitter sequence.
    neutral = "#717171"

    state_colors = {
        "MCG1": "#D8574F",
        "MCG2": "#159D88",
        "MAC1": "#277DA1",
        "MAC2": "#72B7B2",
        "E-MDSC": "#E09F3E",
        "M-MDSC": "#8C6D5A",
    }

    text = "#252525"

    for i, state in enumerate(state_order):
        values = delta_table.loc[delta_table["state"].eq(state), "delta"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        color = state_colors[state]
        ax.fill_between(grid, i, i + density * 0.34, color=color, alpha=0.4, linewidth=0)
        ax.plot(grid, i + density * 0.34, color=color, lw=0.9)

        # Keep the seed at this point to reproduce the sampling or jitter sequence.
        jitter = np.random.default_rng(20260713 + 100 + i).uniform(-0.3, -0.08, len(values))
        ax.scatter(
            values, i + jitter, s=14, color=color, edgecolor="white", linewidth=0.35, alpha=0.8, zorder=3
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.hlines(i - 0.02, q1, q3, color=text, lw=2.0, zorder=4)
        ax.scatter(
            median, i - 0.02, s=20, marker="D", color=text, edgecolor="white", linewidth=0.45, zorder=5
        )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(x_max, i, f"n={len(values)}", ha="right", va="center", fontsize=6.0, color=neutral)
    ax.axvline(0, color="#929292", lw=0.7, ls="--")
    ax.set_yticks(range(len(state_order)), state_order)
    ax.set_ylim(-0.45, len(state_order) - 0.45)
    ax.invert_yaxis()

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Within-patient Miller-IM score delta (state mean - other states mean)")
    ax.grid(axis="x", color="#E6E6E6", lw=0.45)
    ax.tick_params(axis="y", length=0)
    clean_axis(ax, left=False, bottom=True)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    fig.text(0.17, 0.965, "GSE278456 patient-level state contrasts", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.17,
        0.91,
        "Cloud: density | dots: patients | diamond and line: median and IQR",
        fontsize=6.4,
        color=neutral,
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(
        "Figure2E", stem, "GSE278456 patient delta raincloud across published states", source, pdf, png
    )
    return result


# Write the source table with its original columns, ordering and compression, then check that it is nonempty.
if __name__ == "__main__":

    # Keep the seed at this point to reproduce the sampling or jitter sequence.
    np.random.seed(20260713)
    configure_style()
    (fig_root / "Figure1").mkdir(parents=True, exist_ok=True)
    (write_root / "Figure1" / "source_data").mkdir(parents=True, exist_ok=True)
    panel_record_1 = figure1_design()

    figure1_records = [panel_record_1]
    figure_name = "Figure1"
    records = figure1_records
    out = write_root / figure_name
    out.mkdir(parents=True, exist_ok=True)

    # Build the table from the named records in their original order and with the same unit per row.
    output_table = pd.DataFrame(records)
    csv_path = out / f"{figure_name}_panel_manifest.csv"
    output_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
