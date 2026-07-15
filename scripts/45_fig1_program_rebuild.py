#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw, ImageFont, ImageOps
from scipy.stats import spearmanr

from recurrent_figure_style import COLORS, apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP38 = ROOT / "write/38_independent_cohort_mg_inflammatory_recalculation"
STEP41_SOURCE = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure1/source_data"
STEP42_WRITE = ROOT / "write/42_visual_story_rebuild/Figure1"
STEP42_FIG = ROOT / "figures/42_visual_story_rebuild/Figure1"

WRITE_ROOT = ROOT / "write/45_figure1_program_rebuild/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/45_figure1_program_rebuild/Figure1"

DATASETS = ["GSE174554", "GSE274546"]
DATASET_COLORS = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}
PAIR_COUNTS = {"GSE174554": 18, "GSE274546": 45}
RAW20_NAME = "Miller_Microglial_Inflammatory_raw_top20"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def format_p(value: float) -> str:
    """Format nominal fgsea P values without implying multiplicity correction."""
    if value < 1e-4:
        return "<1e-4"
    if value < 0.01:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return f"{value:.3f}".rstrip("0").rstrip(".")


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    metadata = {
        "Creator": "Step45 recurrent GBM Figure 1 program rebuild",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step45 recurrent GBM Figure 1 program rebuild"},
    )
    plt.close(fig)
    return pdf, png


def copy_locked_panel() -> dict[str, str]:
    stem = "Figure1A_compact_paired_design"
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    source = SOURCE_OUT / f"{stem}_source.csv"
    png = FIG_OUT / f"{stem}.png"
    pdf = FIG_OUT / f"{stem}.pdf"
    shutil.copy2(STEP42_WRITE / "source_data" / f"{stem}_source.csv", source)
    shutil.copy2(STEP42_FIG / f"{stem}.png", png)
    shutil.copy2(STEP42_FIG / f"{stem}.pdf", pdf)
    return panel_record("Figure1A", stem, "Locked independent paired-cohort study design", source, pdf, png)


def load_gsea_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    curves = pd.read_csv(STEP41_SOURCE / "raw20_gsea_curves.csv.gz")
    statistics = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    statistics = statistics.loc[
        statistics["dataset"].isin(DATASETS)
        & statistics["threshold"].eq(20)
        & statistics["pathway"].eq(RAW20_NAME)
        & statistics["formal_testable"].astype(bool)
    ].copy()
    if len(statistics) != 2:
        raise ValueError(f"Expected two formal raw20 GSEA rows, found {len(statistics)}")
    return curves, statistics


def panel_gsea(
    dataset: str,
    panel: str,
    curves: pd.DataFrame,
    statistics: pd.DataFrame,
    running_limits: tuple[float, float],
    rank_limit: float,
) -> dict[str, str]:
    stem = f"{panel}_{dataset}_raw20_gsea"
    curve = curves.loc[curves["dataset"].eq(dataset)].sort_values("rank").copy()
    row = statistics.loc[statistics["dataset"].eq(dataset)].iloc[0]
    curve["NES"] = float(row["NES"])
    curve["nominal_P"] = float(row["pval"])
    curve["n_pairs"] = PAIR_COUNTS[dataset]
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    curve.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    fig = new_figure(86, 74)
    grid = fig.add_gridspec(
        3,
        1,
        left=0.16,
        right=0.96,
        bottom=0.17,
        top=0.86,
        height_ratios=[3.6, 0.45, 1.0],
        hspace=0.07,
    )
    enrichment_ax = fig.add_subplot(grid[0])
    hit_ax = fig.add_subplot(grid[1], sharex=enrichment_ax)
    rank_ax = fig.add_subplot(grid[2], sharex=enrichment_ax)

    ranks = curve["rank"].to_numpy(dtype=float)
    running = curve["running_enrichment"].to_numpy(dtype=float)
    rank_stat = curve["rank_stat"].to_numpy(dtype=float)
    color = DATASET_COLORS[dataset]

    enrichment_ax.plot(ranks, running, color=color, linewidth=1.55)
    enrichment_ax.fill_between(ranks, 0, running, color=color, alpha=0.14, linewidth=0)
    enrichment_ax.axhline(0, color="#999999", linewidth=0.55)
    enrichment_ax.set_ylim(*running_limits)
    enrichment_ax.set_ylabel("Running ES")
    enrichment_ax.set_xticks([])
    clean_axis(enrichment_ax, keep_left=True, keep_bottom=False)
    enrichment_ax.text(
        0.0,
        1.05,
        dataset,
        transform=enrichment_ax.transAxes,
        fontsize=8.7,
        fontweight="bold",
        color=color,
        ha="left",
        va="bottom",
    )
    enrichment_ax.text(
        1.0,
        1.05,
        f"NES {float(row['NES']):.2f}   P {format_p(float(row['pval']))}\n"
        f"{PAIR_COUNTS[dataset]} paired patients",
        transform=enrichment_ax.transAxes,
        fontsize=6.6,
        color="#555555",
        ha="right",
        va="bottom",
        linespacing=1.15,
    )

    hits = curve.loc[curve["hit"].astype(bool), "rank"].to_numpy(dtype=float)
    hit_ax.vlines(hits, 0, 1, color="#333333", linewidth=0.48)
    hit_ax.set_ylim(0, 1)
    hit_ax.axis("off")

    rank_ax.fill_between(ranks, 0, np.clip(rank_stat, 0, None), color=COLORS["recurrent"], alpha=0.78, linewidth=0)
    rank_ax.fill_between(ranks, 0, np.clip(rank_stat, None, 0), color=COLORS["primary"], alpha=0.78, linewidth=0)
    rank_ax.axhline(0, color="#777777", linewidth=0.45)
    rank_ax.set_ylim(-rank_limit, rank_limit)
    rank_ax.set_yticks([])
    rank_ax.set_xlabel("Ranked transcriptome")
    rank_ax.text(
        0.0,
        -0.62,
        "Recurrent-enriched",
        transform=rank_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        color=COLORS["recurrent"],
    )
    rank_ax.text(
        1.0,
        -0.62,
        "Primary-enriched",
        transform=rank_ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.5,
        color=COLORS["primary"],
    )
    clean_axis(rank_ax, keep_left=False, keep_bottom=True)
    rank_ax.spines["bottom"].set_color("#777777")
    pdf, png = save_panel(fig, stem)
    return panel_record(panel, stem, f"{dataset} paired pseudobulk raw20 GSEA", source, pdf, png)


def panel_rank_rank(statistics: pd.DataFrame) -> dict[str, str]:
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    source = SOURCE_OUT / f"{stem}_source.csv"
    frame = pd.read_csv(STEP42_WRITE / "source_data" / f"{stem}_source.csv")
    shared_genes = set(shared_leading_edge_table(statistics)["gene"])
    frame["shared_leading_edge"] = frame["gene"].isin(shared_genes)
    frame["labelled"] = frame["shared_leading_edge"]
    frame.to_csv(source, index=False)

    fig = new_figure(92, 83)
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
    background = frame.loc[~frame["miller_raw20"]].sort_values("density")
    raw20_other = frame.loc[frame["miller_raw20"] & ~frame["shared_leading_edge"]]
    shared = frame.loc[frame["shared_leading_edge"]]
    ax.scatter(
        background["rank_stat_GSE174554"],
        background["rank_stat_GSE274546"],
        c=background["density"],
        cmap="Greys",
        s=2.1,
        alpha=0.70,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        raw20_other["rank_stat_GSE174554"],
        raw20_other["rank_stat_GSE274546"],
        s=21,
        facecolor="white",
        edgecolor="#D8574F",
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared["rank_stat_GSE174554"],
        shared["rank_stat_GSE274546"],
        s=27,
        color="#D8574F",
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )
    label_offsets = {
        "CCL3": (4, 7),
        "CCL4": (-25, 7),
        "CH25H": (5, 5),
        "SGK1": (-23, -10),
        "FOLR2": (5, 7),
        "PDK4": (5, -10),
        "KLF6": (-18, -11),
        "SIGLEC8": (5, 4),
    }
    for _, row in shared.iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=label_offsets.get(row["gene"], (3, 3)),
            textcoords="offset points",
            fontsize=5.5,
            color="#A63D38",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = float(
        np.nanmax(np.abs(frame[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))
    ) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    clean_axis(ax)
    rho, _ = spearmanr(frame["rank_stat_GSE174554"], frame["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(frame):,}\nSpearman rho = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.5,
        color="#717171",
    )
    ax.text(
        0.98,
        0.03,
        f"8 shared leading-edge genes\n{len(raw20_other)} other measured raw20 genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.1,
        color="#A63D38",
        linespacing=1.15,
    )
    fig.text(0.16, 0.965, "Cross-cohort gene-rank comparison", fontsize=8.7, fontweight="bold", va="top")
    fig.text(
        0.16,
        0.915,
        "Transcriptome background with shared raw20 leading-edge genes highlighted",
        fontsize=6.2,
        color="#717171",
        va="top",
    )
    pdf, png = save_panel(fig, stem)
    return panel_record(
        "Figure1D",
        stem,
        "All-gene rank-stat comparison with eight shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )


def shared_leading_edge_table(statistics: pd.DataFrame) -> pd.DataFrame:
    leading_sets = {
        row["dataset"]: set(str(row["leadingEdge"]).split(";"))
        for _, row in statistics.iterrows()
    }
    shared = set.intersection(*(leading_sets[dataset] for dataset in DATASETS))
    if len(shared) != 8:
        raise ValueError(f"Expected eight shared leading-edge genes, found {sorted(shared)}")

    provenance = pd.read_csv(STEP38 / "fixed_program_provenance.csv")
    raw20_genes = provenance.loc[provenance["signature"].eq(RAW20_NAME), "genes"].iloc[0].split(";")
    order = [gene for gene in raw20_genes if gene in shared]

    direction = pd.read_csv(STEP38 / "independent_fixed_program_gene_direction.csv")
    direction = direction.loc[
        direction["dataset"].isin(DATASETS)
        & direction["threshold"].eq(20)
        & direction["signature"].eq(RAW20_NAME)
        & direction["gene"].isin(shared)
    ].copy()
    if len(direction) != 16:
        raise ValueError(f"Expected sixteen gene-by-cohort rows, found {len(direction)}")
    direction["shared_leading_edge"] = True
    direction["gene_order"] = direction["gene"].map({gene: idx for idx, gene in enumerate(order)})
    direction["dataset_order"] = direction["dataset"].map({dataset: idx for idx, dataset in enumerate(DATASETS)})
    direction = direction.sort_values(["gene_order", "dataset_order"]).reset_index(drop=True)
    return direction


def panel_shared_leading_edge_heatmap(statistics: pd.DataFrame) -> dict[str, str]:
    stem = "Figure1E_shared_leading_edge_heatmap"
    direction = shared_leading_edge_table(statistics)
    source = SOURCE_OUT / f"{stem}_source.csv"
    direction.to_csv(source, index=False)

    order = direction.sort_values("gene_order")["gene"].drop_duplicates().tolist()
    matrix = direction.pivot(index="gene", columns="dataset", values="logFC").reindex(index=order, columns=DATASETS)
    max_abs = float(np.nanmax(np.abs(matrix.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    fig = new_figure(86, 80)
    ax = fig.add_axes([0.24, 0.22, 0.55, 0.56])
    image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(
        range(len(DATASETS)),
        ["GSE174554\n18 pairs", "GSE274546\n45 pairs"],
        fontsize=6.8,
    )
    ax.set_yticks(range(len(order)), order, fontsize=7.1)
    ax.tick_params(length=0)
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(DATASETS):
            value = float(matrix.loc[gene, dataset])
            text_color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=text_color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(DATASETS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    color_ax = fig.add_axes([0.82, 0.30, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)
    colorbar.outline.set_linewidth(0.45)

    gsea_rows = statistics.set_index("dataset")
    fig.text(0.24, 0.95, "Shared leading-edge genes", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.24,
        0.90,
        "Intersection of the two formal raw20 GSEA leading edges",
        fontsize=6.2,
        color="#717171",
        va="top",
    )
    fig.text(
        0.24,
        0.145,
        f"Cohort GSEA: NES {gsea_rows.loc['GSE174554', 'NES']:.2f} / {gsea_rows.loc['GSE274546', 'NES']:.2f}\n"
        f"Nominal P {format_p(float(gsea_rows.loc['GSE174554', 'pval']))} / "
        f"{format_p(float(gsea_rows.loc['GSE274546', 'pval']))}",
        fontsize=6.1,
        color="#555555",
        va="top",
        linespacing=1.25,
    )
    fig.text(0.24, 0.065, "0/16 individual gene tests reached FDR < 0.05", fontsize=5.9, color="#777777", va="top")

    pdf, png = save_panel(fig, stem)
    return panel_record("Figure1E", stem, "Eight shared leading-edge genes across the two recurrence cohorts", source, pdf, png)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            try:
                return ImageFont.truetype(str(path), size=size)
            except OSError:
                continue
    return ImageFont.load_default()


def paste_contained(canvas: Image.Image, source: Path, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    width = right - left
    height = bottom - top
    with Image.open(source) as panel:
        panel = panel.convert("RGB")
        fitted = ImageOps.contain(panel, (width, height), Image.Resampling.LANCZOS)
        x = left + (width - fitted.width) // 2
        y = top + (height - fitted.height) // 2
        canvas.paste(fitted, (x, y))


def make_contact_sheet(records: list[dict[str, str]]) -> tuple[Path, Path, Path]:
    panel_paths = {record["panel"]: Path(record["png"]) for record in records}
    layout_rows = [
        {"panel": "Figure1A", "left": 270, "top": 130, "right": 3330, "bottom": 1700},
        {"panel": "Figure1B", "left": 120, "top": 1810, "right": 1740, "bottom": 3150},
        {"panel": "Figure1C", "left": 1860, "top": 1810, "right": 3480, "bottom": 3150},
        {"panel": "Figure1D", "left": 120, "top": 3260, "right": 1740, "bottom": 4830},
        {"panel": "Figure1E", "left": 1860, "top": 3260, "right": 3480, "bottom": 4830},
    ]
    canvas = Image.new("RGB", (3600, 4960), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = font(52, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        paste_contained(
            canvas,
            panel_paths[panel],
            (int(row["left"]), int(row["top"]), int(row["right"]), int(row["bottom"])),
        )
        drawer.text((int(row["left"]), int(row["top"]) - 56), panel[-1], fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1_candidate_contact_sheet.png"
    pdf = FIG_OUT / "Figure1_candidate_contact_sheet.pdf"
    canvas.save(png, format="PNG", dpi=(300, 300), optimize=True)
    contact_fig = plt.figure(figsize=(12, 4960 / 300), facecolor="white")
    contact_ax = contact_fig.add_axes([0, 0, 1, 1])
    contact_ax.imshow(np.asarray(canvas), interpolation="nearest")
    contact_ax.axis("off")
    contact_fig.savefig(
        pdf,
        dpi=300,
        facecolor="white",
        edgecolor="none",
        metadata={
            "Creator": "Step45 recurrent GBM Figure 1 program rebuild",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(contact_fig)
    layout = pd.DataFrame(layout_rows)
    source = SOURCE_OUT / "Figure1_candidate_contact_sheet_layout.csv"
    layout.to_csv(source, index=False)
    return pdf, png, source


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, str]:
    return {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": sha256(source),
        "pdf_sha256": sha256(pdf),
        "png_sha256": sha256(png),
    }


def main() -> None:
    np.random.seed(20260713)
    apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    curves, statistics = load_gsea_inputs()
    running_min = float(curves["running_enrichment"].min())
    running_max = float(curves["running_enrichment"].max())
    running_pad = 0.06 * (running_max - running_min)
    running_limits = (running_min - running_pad, running_max + running_pad)
    rank_limit = float(np.nanmax(np.abs(curves["rank_stat"].to_numpy()))) * 1.04

    records = [
        copy_locked_panel(),
        panel_gsea("GSE174554", "Figure1B", curves, statistics, running_limits, rank_limit),
        panel_gsea("GSE274546", "Figure1C", curves, statistics, running_limits, rank_limit),
        panel_rank_rank(statistics),
        panel_shared_leading_edge_heatmap(statistics),
    ]
    contact_pdf, contact_png, contact_source = make_contact_sheet(records)
    records.append(
        panel_record(
            "Figure1_preview",
            "Figure1_candidate_contact_sheet",
            "Review-only composite of Figure 1A-E",
            contact_source,
            contact_pdf,
            contact_png,
        )
    )
    manifest = pd.DataFrame(records)
    manifest.to_csv(WRITE_ROOT / "Figure1_panel_manifest.csv", index=False)
    print("STEP45_FIGURE1_COMPLETE panels=5 preview=1")


if __name__ == "__main__":
    main()
