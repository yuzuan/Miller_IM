#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr


ROOT = Path(__file__).resolve().parents[1]
STEP41 = ROOT / "write/41_mg_inflammatory_sci_rebuild"
STEP38 = ROOT / "write/38_independent_cohort_mg_inflammatory_recalculation"
FIG_ROOT = ROOT / "figures/42_visual_story_rebuild"
WRITE_ROOT = ROOT / "write/42_visual_story_rebuild"

DATASETS = ["GSE174554", "GSE274546"]
DATASET_COLORS = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}
MYELOID_ATLAS_SOURCES = {
    "GSE174554": ROOT / "write/43_scgbm_style_clean_myeloid_reclustering/GSE174554/clean_myeloid_umap_source.csv.gz",
    "GSE274546": ROOT / "write/43_scgbm_style_clean_myeloid_reclustering/GSE274546/clean_myeloid_umap_source.csv.gz",
}
MYELOID_STATE_ORDER = [
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
MYELOID_STATE_COLORS = {
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
MYELOID_LABEL_OFFSETS = {
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
IDENTITY_COLORS = {
    "Microglia": "#159D88",
    "Macrophage": "#277DA1",
    "Monocyte": "#E09F3E",
    "CDC": "#B565A7",
    "Neutrophil": "#6F6F6F",
}
STATE_ORDER = ["MCG1", "MCG2", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]
STATE_COLORS = {
    "MCG1": "#D8574F",
    "MCG2": "#159D88",
    "MAC1": "#277DA1",
    "MAC2": "#72B7B2",
    "E-MDSC": "#E09F3E",
    "M-MDSC": "#8C6D5A",
}
TEXT = "#252525"
NEUTRAL = "#717171"
GRID = "#E6E6E6"
SEED = 20260713
MM_PER_INCH = 25.4


def mm(value: float) -> float:
    return value / MM_PER_INCH


def configure_style() -> None:
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
            "axes.edgecolor": TEXT,
            "axes.labelcolor": TEXT,
            "text.color": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def clean_axis(ax: mpl.axes.Axes, *, left: bool = True, bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(width=0.6, length=2.4)


def save_panel(fig: mpl.figure.Figure, figure_name: str, stem: str) -> tuple[Path, Path]:
    out = FIG_ROOT / figure_name
    out.mkdir(parents=True, exist_ok=True)
    pdf = out / f"{stem}.pdf"
    png = out / f"{stem}.png"
    metadata = {"Creator": "Step42 recurrent GBM visual story", "CreationDate": None, "ModDate": None}
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(png, dpi=600, facecolor="white", edgecolor="none", metadata={"Software": "Step42 recurrent GBM visual story"})
    plt.close(fig)
    return pdf, png


def write_source(frame: pd.DataFrame, figure_name: str, stem: str) -> Path:
    out = WRITE_ROOT / figure_name / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}_source.csv"
    frame.to_csv(path, index=False)
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, object]:
    return {
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


def label_position(frame: pd.DataFrame) -> tuple[float, float]:
    points = frame[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
    center = np.nanmedian(points, axis=0)
    distance = np.sqrt(((points - center) ** 2).sum(axis=1))
    core = points[distance <= np.nanquantile(distance, 0.60)]
    target = np.nanmedian(core, axis=0) if len(core) else center
    return float(target[0]), float(target[1])


def figure1_design() -> dict[str, object]:
    stem = "Figure1A_compact_paired_design"
    audit = pd.read_csv(STEP38 / "independent_input_pair_audit.csv")
    source_frame = audit.loc[
        audit["dataset"].isin(DATASETS) & audit["threshold"].eq(20)
    ].copy()
    if source_frame.groupby("dataset").size().to_dict() != {
        "GSE174554": 1,
        "GSE274546": 1,
    }:
        raise ValueError("Step38 must contain one formal threshold-20 audit row per cohort.")
    source_frame["input_libraries_or_matrices"] = source_frame["dataset"].map(
        {"GSE174554": 91, "GSE274546": 111}
    )
    source_frame["clean_myeloid"] = (
        source_frame["n_cells_primary"].astype(int)
        + source_frame["n_cells_recurrent"].astype(int)
    )
    source_frame["formal_pairs_threshold20"] = source_frame["n_pairs"].astype(int)
    observed = (
        source_frame.set_index("dataset")[
            ["formal_pairs_threshold20", "clean_myeloid"]
        ]
        .astype(int)
        .to_dict("index")
    )
    expected = {
        "GSE174554": {"formal_pairs_threshold20": 17, "clean_myeloid": 12489},
        "GSE274546": {"formal_pairs_threshold20": 45, "clean_myeloid": 55568},
    }
    if observed != expected:
        raise ValueError(f"Unexpected IDH-restricted Figure1A analysis set: {observed}")
    source_frame = source_frame.sort_values("dataset").reset_index(drop=True)
    source_frame["modality"] = "snRNA-seq"
    source_frame["raw_input_unit"] = ["10x matrices", "10x libraries"]
    source_frame["minimum_cells_per_condition"] = 20
    source_frame["processing_scope"] = "Each cohort reconstructed, annotated and tested independently"
    source_frame["statistical_testing"] = "Patient pseudobulk, paired edgeR and Miller-IM GSEA performed separately within each cohort"
    source_frame["replication_assessment"] = "Cohort results compared only after independent testing"
    source_frame["workflow"] = (
        "all-cell reconstruction -> pan-myeloid re-clustering -> marker and cluster-DEG review -> "
        "patient raw-count pseudobulk -> paired edgeR -> prespecified Miller-IM program GSEA"
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
    source_frame = source_frame.loc[:, safe_source_columns]
    source = write_source(source_frame, "Figure1", stem)

    fig = plt.figure(figsize=(mm(183), mm(94)), facecolor="white")
    ax = fig.add_axes([0.025, 0.025, 0.95, 0.95])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def rounded_box(x0: float, x1: float, y0: float, y1: float, edge: str, fill: str, lw: float = 0.8) -> None:
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x0, y0), x1 - x0, y1 - y0,
                boxstyle="round,pad=0.02,rounding_size=0.7",
                facecolor=fill, edgecolor=edge, linewidth=lw,
            )
        )

    cohort_boxes = [(4.0, 47.5), (52.5, 96.0)]
    for row, (x0, x1) in zip(source_frame.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = DATASET_COLORS[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        rounded_box(x0, x1, 80.0, 97.0, color, pale, lw=1.0)
        ax.text(
            center, 94.5,
            f"{row['dataset']}  |  snRNA-seq",
            fontsize=8.3, weight="bold", color=color, ha="center", va="top",
        )
        ax.text(center, 89.3, f"{int(row['input_libraries_or_matrices'])} {row['raw_input_unit']}", fontsize=6.1, color=NEUTRAL, ha="center", va="top")
        primary_x, recurrent_x = center - 8.0, center + 8.0
        rounded_box(primary_x - 4.9, primary_x + 4.9, 81.7, 86.7, color, "white", lw=0.9)
        rounded_box(recurrent_x - 5.5, recurrent_x + 5.5, 81.7, 86.7, color, color, lw=0.9)
        ax.text(primary_x, 84.2, "Primary", fontsize=6.3, color=TEXT, ha="center", va="center")
        ax.text(recurrent_x, 84.2, "1st recurrence", fontsize=6.3, color="white", ha="center", va="center")
        ax.annotate("", xy=(recurrent_x - 5.8, 84.2), xytext=(primary_x + 5.2, 84.2), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0))

    for center in [sum(box) / 2 for box in cohort_boxes]:
        ax.annotate("", xy=(center, 74.8), xytext=(center, 79.5), arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85))

    workflow_boxes = [
        (8.0, 92.0, 64.0, 74.0, "Independent all-cell reconstruction", "Library QC · doublet removal · basic lineage annotation"),
        (8.0, 92.0, 50.0, 60.0, "Cohort-wise pan-myeloid re-clustering", "2,000 HVGs · PCA · neighbor graph · Leiden"),
        (6.0, 94.0, 36.0, 46.0, "Marker and cluster-DEG review", "Remove lymphoid/tumour contaminants and low-quality cells"),
    ]
    for index, (x0, x1, y0, y1, title, subtitle) in enumerate(workflow_boxes):
        rounded_box(x0, x1, y0, y1, "#858585", "white", lw=0.75)
        ax.text(50, y1 - 2.2, title, fontsize=7.3, weight="bold", ha="center", va="top")
        ax.text(50, y0 + 2.2, subtitle, fontsize=6.1, color=NEUTRAL, ha="center", va="bottom")
        if index < len(workflow_boxes) - 1:
            next_y1 = workflow_boxes[index + 1][3]
            ax.annotate("", xy=(50, next_y1 + 0.5), xytext=(50, y0 - 0.5), arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85))

    analysis_y0, analysis_y1 = 14.5, 29.5
    for row, (x0, x1) in zip(source_frame.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = DATASET_COLORS[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        source_x = 35.0 if center < 50 else 65.0
        ax.annotate(
            "", xy=(center, analysis_y1 + 0.6), xytext=(source_x, 35.5),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, connectionstyle="arc3,rad=0.08" if center < 50 else "arc3,rad=-0.08"),
        )
        rounded_box(x0, x1, analysis_y0, analysis_y1, color, pale, lw=1.0)
        ax.text(center, 26.8, f"{row['dataset']}  |  {int(row['formal_pairs_threshold20'])} paired patients", fontsize=7.1, weight="bold", color=color, ha="center", va="center")
        ax.text(center, 22.0, f"{int(row['clean_myeloid']):,} analyzed myeloid cells", fontsize=5.9, color=NEUTRAL, ha="center", va="center")
        ax.text(center, 17.2, "Raw-count pseudobulk · paired edgeR · Miller-IM GSEA", fontsize=5.9, color=TEXT, ha="center", va="center")

    ax.annotate("", xy=(46.8, 10.8), xytext=(sum(cohort_boxes[0]) / 2, analysis_y0 - 0.5), arrowprops=dict(arrowstyle="-|>", color=DATASET_COLORS["GSE174554"], lw=0.9, connectionstyle="arc3,rad=-0.08"))
    ax.annotate("", xy=(53.2, 10.8), xytext=(sum(cohort_boxes[1]) / 2, analysis_y0 - 0.5), arrowprops=dict(arrowstyle="-|>", color=DATASET_COLORS["GSE274546"], lw=0.9, connectionstyle="arc3,rad=0.08"))
    rounded_box(22.0, 78.0, 1.0, 11.8, "#159D88", mpl.colors.to_rgba("#159D88", 0.085), lw=1.0)
    ax.text(50, 8.8, "Cross-cohort replication assessment", fontsize=7.2, weight="bold", color=TEXT, ha="center", va="center")
    ax.text(50, 3.8, "Compare effect direction and nominal P after independent cohort tests", fontsize=5.9, color=NEUTRAL, ha="center", va="center")

    pdf, png = save_panel(fig, "Figure1", stem)
    return panel_record("Figure1A", stem, "Independent paired-cohort myeloid workflow", source, pdf, png)


def figure1_clean_myeloid_umap(dataset: str, panel: str) -> dict[str, object]:
    stem = f"{panel}_{dataset}_clean_myeloid_umap"
    frame = pd.read_csv(MYELOID_ATLAS_SOURCES[dataset])
    unknown = sorted(set(frame["curated_myeloid_subtype"]) - set(MYELOID_STATE_ORDER))
    if unknown:
        raise ValueError(f"Unexpected reviewed myeloid states for {dataset}: {unknown}")
    source = write_source(frame, "Figure1", stem)

    fig = plt.figure(figsize=(mm(91), mm(88)), facecolor="white")
    ax = fig.add_axes([0.035, 0.205, 0.93, 0.64])
    counts = frame["curated_myeloid_subtype"].value_counts()
    x_span = float(frame["UMAP1"].max() - frame["UMAP1"].min())
    y_span = float(frame["UMAP2"].max() - frame["UMAP2"].min())
    point_size = 0.62 if len(frame) < 30_000 else 0.34
    for state in counts.sort_values(ascending=False).index:
        current = frame[frame["curated_myeloid_subtype"].eq(state)]
        ax.scatter(
            current["UMAP1"], current["UMAP2"],
            s=point_size, color=MYELOID_STATE_COLORS[state], alpha=0.62,
            linewidths=0, rasterized=True,
        )
    for state in MYELOID_STATE_ORDER:
        current = frame[frame["curated_myeloid_subtype"].eq(state)]
        if current.empty:
            continue
        largest_cluster = current["leiden_myeloid"].astype(str).value_counts().index[0]
        label_cells = current[current["leiden_myeloid"].astype(str).eq(largest_cluster)]
        x, y = label_position(label_cells)
        dx, dy = MYELOID_LABEL_OFFSETS.get(dataset, {}).get(state, (0.0, 0.0))
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
                color=MYELOID_STATE_COLORS[state],
                lw=0.45,
                alpha=0.75,
                zorder=4,
            )
        ax.text(
            text_x, text_y, short_label, fontsize=5.2, weight="bold", color=MYELOID_STATE_COLORS[state],
            ha="center", va="center",
            path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
            zorder=5,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="none", markersize=3.8, markerfacecolor=MYELOID_STATE_COLORS[state], markeredgewidth=0, label=state)
        for state in MYELOID_STATE_ORDER if state in counts.index
    ]
    fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.025),
        ncol=3, frameon=False, handletextpad=0.35, columnspacing=1.0,
        fontsize=5.25,
    )
    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=DATASET_COLORS[dataset])
    fig.text(
        0.05,
        0.905,
        f"{len(frame):,} clean myeloid cells  |  {frame['leiden_myeloid'].nunique()} clusters",
        fontsize=6.2,
        color=NEUTRAL,
        va="top",
    )
    pdf, png = save_panel(fig, "Figure1", stem)
    return panel_record(panel, stem, f"{dataset} independently integrated and DEG-reviewed clean-myeloid UMAP", source, pdf, png)


def rank_density(frame: pd.DataFrame, bins: int = 75) -> np.ndarray:
    x = frame["rank_stat_GSE174554"].to_numpy(dtype=float)
    y = frame["rank_stat_GSE274546"].to_numpy(dtype=float)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)
    return np.log1p(counts[xi, yi])


def figure1_rank_rank() -> dict[str, object]:
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    curves = pd.read_csv(STEP41 / "Figure1/source_data/raw20_gsea_curves.csv.gz")
    left = curves[curves["dataset"].eq("GSE174554")][["gene", "rank", "rank_stat", "hit"]].rename(columns={"rank": "rank_GSE174554", "rank_stat": "rank_stat_GSE174554", "hit": "program_GSE174554"})
    right = curves[curves["dataset"].eq("GSE274546")][["gene", "rank", "rank_stat", "hit"]].rename(columns={"rank": "rank_GSE274546", "rank_stat": "rank_stat_GSE274546", "hit": "program_GSE274546"})
    frame = left.merge(right, on="gene", how="inner", validate="one_to_one")
    frame["miller_raw20"] = frame[["program_GSE174554", "program_GSE274546"]].any(axis=1)
    frame["density"] = rank_density(frame)
    program = frame[frame["miller_raw20"]].copy()
    program["label_score"] = program["rank_stat_GSE174554"].abs() + program["rank_stat_GSE274546"].abs()
    label_genes = set(program.nlargest(min(6, len(program)), "label_score")["gene"])
    frame["labelled"] = frame["gene"].isin(label_genes)
    source = write_source(frame, "Figure1", stem)

    fig = plt.figure(figsize=(mm(92), mm(83)), facecolor="white")
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
    background = frame[~frame["miller_raw20"]].sort_values("density")
    ax.scatter(background["rank_stat_GSE174554"], background["rank_stat_GSE274546"], c=background["density"], cmap="Greys", s=2.1, alpha=0.70, linewidths=0, rasterized=True)
    ax.scatter(program["rank_stat_GSE174554"], program["rank_stat_GSE274546"], s=24, color="#D8574F", edgecolor="white", linewidth=0.65, zorder=4)
    label_offsets = {
        "CCL3": (3, 7), "CCL4": (-26, 7), "CH25H": (5, 5), "SGK1": (-25, -10),
        "FOLR2": (5, 7), "BHLHE41": (5, -8),
    }
    for _, row in program[program["gene"].isin(label_genes)].iterrows():
        offset = label_offsets.get(row["gene"], (3, 3))
        ax.annotate(row["gene"], (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]), xytext=offset, textcoords="offset points", fontsize=5.7, color="#A63D38", path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
    ax.axhline(0, color="#A4A4A4", lw=0.55)
    ax.axvline(0, color="#A4A4A4", lw=0.55)
    limit = float(np.nanmax(np.abs(frame[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color=GRID, lw=0.42, zorder=0)
    clean_axis(ax)
    rho, _ = spearmanr(frame["rank_stat_GSE174554"], frame["rank_stat_GSE274546"])
    ax.text(0.02, 0.98, f"All shared genes: {len(frame):,}\nSpearman rho = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6.5, color=NEUTRAL)
    ax.text(0.98, 0.03, f"Miller-IM: {len(program)} measured genes", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4, color="#A63D38")
    fig.text(0.16, 0.965, "Cross-cohort all-gene concordance", fontsize=8.7, weight="bold", va="top")
    fig.text(0.16, 0.915, "Transcriptome density with Miller-IM genes highlighted", fontsize=6.4, color=NEUTRAL, va="top")
    pdf, png = save_panel(fig, "Figure1", stem)
    return panel_record("Figure1D", stem, "All-gene rank-stat density with Miller-IM genes highlighted", source, pdf, png)


def figure2_raw20_umap(dataset: str, panel: str, vmin: float, vmax: float) -> dict[str, object]:
    stem = f"{panel}_{dataset}_raw20_terrain"
    frame = pd.read_csv(STEP41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz")
    source = write_source(frame, "Figure2", stem)
    ordered = frame.sort_values("raw20_score")

    fig = plt.figure(figsize=(mm(88), mm(78)), facecolor="white")
    ax = fig.add_axes([0.045, 0.035, 0.80, 0.78])
    points = ax.scatter(ordered["UMAP1"], ordered["UMAP2"], c=ordered["raw20_score"], s=0.55, cmap="magma", vmin=vmin, vmax=vmax, alpha=0.94, linewidths=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    cax = fig.add_axes([0.88, 0.20, 0.023, 0.48])
    cb = fig.colorbar(points, cax=cax)
    cb.ax.set_title("Miller-IM\nscore", fontsize=6.2, pad=4, loc="center")
    cb.ax.tick_params(labelsize=6.0, length=2)
    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=DATASET_COLORS[dataset])
    fig.text(0.05, 0.905, "Continuous program terrain", fontsize=6.6, color=NEUTRAL, va="top")
    pdf, png = save_panel(fig, "Figure2", stem)
    return panel_record(panel, stem, f"{dataset} continuous Miller-IM score terrain", source, pdf, png)


def ridge_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    values = values[np.isfinite(values)]
    if len(values) < 2 or np.isclose(np.std(values), 0):
        return np.zeros_like(grid)
    density = gaussian_kde(values)(grid)
    peak = density.max()
    return density / peak if peak > 0 else density


def figure2_ridgeline(dataset: str, panel: str, x_limits: tuple[float, float]) -> dict[str, object]:
    stem = f"{panel}_{dataset}_cluster_ridgeline"
    raw = pd.read_csv(STEP41 / "Figure2/source_data/fig2_patient_cluster_interval_source.csv")
    frame = raw[raw["dataset"].eq(dataset) & raw["eligible"].astype(bool)].copy()
    labels = pd.read_csv(
        STEP41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz",
        usecols=["paired_cluster", "cluster_label"],
    ).dropna().drop_duplicates()
    if dataset == "GSE274546":
        annotation = pd.read_csv(
            ROOT / "write/36_gse274546_raw_independent_reannotation/04_paired_patient_harmony_clusters/paired_cluster_descriptive_annotation.csv"
        )
        labels = annotation[["cluster", "descriptive_state"]].rename(
            columns={"cluster": "paired_cluster", "descriptive_state": "cluster_label"}
        )
        labels["cluster_label"] = labels.apply(
            lambda row: f"C{int(row['paired_cluster'])} {row['cluster_label']}", axis=1
        )
    labels["paired_cluster"] = pd.to_numeric(labels["paired_cluster"]).astype(int).astype(str)
    label_map = labels.set_index("paired_cluster")["cluster_label"].to_dict()
    frame["paired_cluster"] = pd.to_numeric(frame["paired_cluster"]).astype(int).astype(str)
    frame["cluster_label_display"] = frame["paired_cluster"].map(label_map)
    if frame["cluster_label_display"].isna().any():
        raise ValueError(f"Missing full cluster labels for {dataset}")
    summary = frame.groupby(["paired_cluster", "cluster_label_display"], observed=True)["pooled_patient_raw20"].agg(["median", "count"]).reset_index().sort_values("median")
    order = summary["paired_cluster"].astype(str).tolist()
    frame["ridge_order"] = frame["paired_cluster"].map({cluster: i for i, cluster in enumerate(order)})
    source = write_source(frame, "Figure2", stem)

    fig = plt.figure(figsize=(mm(100), mm(91)), facecolor="white")
    ax = fig.add_axes([0.42, 0.13, 0.55, 0.70])
    grid = np.linspace(x_limits[0], x_limits[1], 300)
    color = DATASET_COLORS[dataset]
    labels = []
    for i, cluster in enumerate(order):
        current = frame[frame["paired_cluster"].eq(cluster)]
        values = current["pooled_patient_raw20"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        baseline = float(i)
        ax.fill_between(grid, baseline, baseline + density * 0.78, color=color, alpha=0.36, linewidth=0)
        ax.plot(grid, baseline + density * 0.78, color=color, lw=0.9)
        jitter = np.random.default_rng(SEED + i).uniform(-0.055, 0.055, len(values))
        ax.scatter(values, baseline - 0.10 + jitter, s=5.5, color=TEXT, alpha=0.55, linewidths=0, zorder=3)
        ax.scatter(np.median(values), baseline - 0.10, marker="D", s=15, color=color, edgecolor="white", linewidth=0.5, zorder=4)
        labels.append(str(current["cluster_label_display"].iloc[0]))
        ax.text(0.985, baseline + 0.49, f"n={len(values)}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.5, color=NEUTRAL)
    ax.set_yticks(range(len(order)), labels)
    ax.set_ylim(-0.35, len(order) - 0.02 + 0.82)
    ax.set_xlim(*x_limits)
    ax.set_xlabel("Patient-level pooled Miller-IM score")
    ax.grid(axis="x", color=GRID, lw=0.45)
    ax.tick_params(axis="y", length=0, labelsize=5.2, pad=2)
    clean_axis(ax, left=False, bottom=True)
    fig.text(0.08, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=color)
    fig.text(0.08, 0.913, "Miller-IM score across unsupervised clusters", fontsize=6.5, color=NEUTRAL, va="top")
    fig.text(0.08, 0.87, "Each dot is one eligible patient-cluster unit", fontsize=6.2, color=NEUTRAL, va="top")
    pdf, png = save_panel(fig, "Figure2", stem)
    return panel_record(panel, stem, f"{dataset} patient-level cluster ridgeline", source, pdf, png)


def figure2_delta_raincloud() -> dict[str, object]:
    stem = "Figure2E_GSE278456_patient_delta_raincloud"
    frame = pd.read_csv(STEP41 / "Figure2/source_data/fig2_gse278456_all_state_patient_deltas.csv")
    frame = frame[frame["state"].isin(STATE_ORDER)].copy()
    frame["state_order"] = frame["state"].map({state: i for i, state in enumerate(STATE_ORDER)})
    source = write_source(frame, "Figure2", stem)

    fig = plt.figure(figsize=(mm(118), mm(82)), facecolor="white")
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.67])
    x_min = min(-0.05, float(frame["delta"].min()) - 0.06)
    x_max = float(frame["delta"].max()) + 0.08
    grid = np.linspace(x_min, x_max, 320)
    for i, state in enumerate(STATE_ORDER):
        values = frame.loc[frame["state"].eq(state), "delta"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        color = STATE_COLORS[state]
        ax.fill_between(grid, i, i + density * 0.34, color=color, alpha=0.40, linewidth=0)
        ax.plot(grid, i + density * 0.34, color=color, lw=0.9)
        jitter = np.random.default_rng(SEED + 100 + i).uniform(-0.30, -0.08, len(values))
        ax.scatter(values, i + jitter, s=14, color=color, edgecolor="white", linewidth=0.35, alpha=0.80, zorder=3)
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        ax.hlines(i - 0.02, q1, q3, color=TEXT, lw=2.0, zorder=4)
        ax.scatter(median, i - 0.02, s=20, marker="D", color=TEXT, edgecolor="white", linewidth=0.45, zorder=5)
        ax.text(x_max, i, f"n={len(values)}", ha="right", va="center", fontsize=6.0, color=NEUTRAL)
    ax.axvline(0, color="#929292", lw=0.7, ls="--")
    ax.set_yticks(range(len(STATE_ORDER)), STATE_ORDER)
    ax.set_ylim(-0.45, len(STATE_ORDER) - 0.45)
    ax.invert_yaxis()
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Within-patient Miller-IM score delta (state mean - other states mean)")
    ax.grid(axis="x", color=GRID, lw=0.45)
    ax.tick_params(axis="y", length=0)
    clean_axis(ax, left=False, bottom=True)
    fig.text(0.17, 0.965, "GSE278456 patient-level state contrasts", fontsize=8.7, weight="bold", va="top")
    fig.text(0.17, 0.91, "Cloud: density | dots: patients | diamond and line: median and IQR", fontsize=6.4, color=NEUTRAL, va="top")
    pdf, png = save_panel(fig, "Figure2", stem)
    return panel_record("Figure2E", stem, "GSE278456 patient delta raincloud across published states", source, pdf, png)


def write_manifest(figure_name: str, records: list[dict[str, object]]) -> None:
    out = WRITE_ROOT / figure_name
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out / f"{figure_name}_panel_manifest.csv", index=False)


def main() -> None:
    np.random.seed(SEED)
    configure_style()
    (FIG_ROOT / "Figure1").mkdir(parents=True, exist_ok=True)
    (WRITE_ROOT / "Figure1" / "source_data").mkdir(parents=True, exist_ok=True)
    figure1_records = [figure1_design()]
    write_manifest("Figure1", figure1_records)
    print("STEP42_FIGURE1_DESIGN_COMPLETE panels=1")


if __name__ == "__main__":
    main()
