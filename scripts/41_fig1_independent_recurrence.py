#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from recurrent_figure_style import COLORS, apply_publication_style, clean_axis, new_figure, save_figure


ROOT = Path(__file__).resolve().parents[1]
OLD_SOURCE = ROOT / "write/39_independent_miller_mg_inflammatory_figures/Figure1/source_data"
STEP38 = ROOT / "write/38_independent_cohort_mg_inflammatory_recalculation"
WRITE_ROOT = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/41_mg_inflammatory_sci_rebuild/Figure1/panel_library"

DATASETS = ["GSE174554", "GSE274546"]
DATASET_COLORS = {
    "GSE174554": COLORS["gse174554"],
    "GSE274546": COLORS["gse274546"],
}
IDENTITY_COLORS = {
    "Microglia": COLORS["microglia"],
    "Macrophage": COLORS["macrophage"],
    "Monocyte": COLORS["monocyte"],
    "CDC": COLORS["cdc"],
    "Neutrophil": COLORS["neutrophil"],
}
DISPLAY_IDENTITY = {"CDC": "cDC"}
MARKER_BLOCKS = {
    "Microglia": ["P2RY12", "TMEM119", "SALL1"],
    "Macrophage": ["GPNMB", "SPP1", "MSR1"],
    "Monocyte": ["VCAN", "FCN1", "S100A8"],
    "cDC": ["FCER1A", "CD1C", "CLEC10A"],
}
PAIR_COUNTS = {"GSE174554": 18, "GSE274546": 45}


def prepare_sources() -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    required = [
        "cohort_workflow.csv",
        "identity_marker_dotplot.csv",
        "raw20_gsea_curves.csv.gz",
        "umap_identity_GSE174554.csv.gz",
        "umap_identity_GSE274546.csv.gz",
    ]
    for name in required:
        src = OLD_SOURCE / name
        if not src.exists():
            raise FileNotFoundError(src)
        shutil.copy2(src, SOURCE_OUT / name)

    gsea = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    gsea = gsea[
        gsea["dataset"].isin(DATASETS)
        & gsea["threshold"].eq(20)
        & gsea["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")
    ].copy()
    if len(gsea) != 2:
        raise ValueError(f"Expected two primary GSEA rows, found {len(gsea)}")
    gsea.to_csv(SOURCE_OUT / "primary_gsea_statistics.csv", index=False)


def panel_study_design() -> None:
    workflow = pd.read_csv(SOURCE_OUT / "cohort_workflow.csv").set_index("dataset")
    fig = new_figure(183, 42)
    ax = fig.add_axes([0.015, 0.08, 0.97, 0.84])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.0, 0.94, "Independent paired-cohort design", fontsize=9, weight="bold", va="top")
    ax.text(
        1.0,
        0.94,
        "Primary vs first recurrence; patient-level inference",
        fontsize=7.2,
        color=COLORS["neutral"],
        va="top",
        ha="right",
    )

    x_positions = {"GSE174554": 0.0, "GSE274546": 0.38}
    widths = 0.34
    subtitles = {"GSE174554": "snRNA-seq", "GSE274546": "snRNA-seq"}
    for dataset in DATASETS:
        x = x_positions[dataset]
        color = DATASET_COLORS[dataset]
        row = workflow.loc[dataset]
        ax.add_patch(mpl.patches.Rectangle((x, 0.20), widths, 0.56, facecolor="#F7F7F7", edgecolor="none"))
        ax.add_patch(mpl.patches.Rectangle((x, 0.20), 0.012, 0.56, facecolor=color, edgecolor="none"))
        ax.text(x + 0.03, 0.69, dataset, fontsize=8.4, weight="bold", color=color, va="top")
        ax.text(x + widths - 0.02, 0.69, subtitles[dataset], fontsize=7, color=COLORS["neutral"], ha="right", va="top")
        metrics = [
            ("Input", f"{int(row.input_libraries_or_matrices):,} libraries/matrices"),
            ("Myeloid", f"{int(row.clean_myeloid):,} cells"),
            ("Paired", f"{int(row.formal_pairs_threshold20)} patients"),
        ]
        for i, (label, value) in enumerate(metrics):
            y = 0.55 - i * 0.14
            ax.text(x + 0.03, y, label.upper(), fontsize=6.2, color=COLORS["neutral"], va="center")
            ax.text(x + 0.12, y, value, fontsize=7.4, weight="bold", va="center")

    ax.annotate("", xy=(0.785, 0.48), xytext=(0.735, 0.48), arrowprops=dict(arrowstyle="-|>", lw=0.9, color=COLORS["neutral"]))
    ax.add_patch(mpl.patches.Rectangle((0.79, 0.20), 0.21, 0.56, facecolor="#F4F4F4", edgecolor="none"))
    ax.text(0.81, 0.67, "PAIRED ANALYSIS", fontsize=6.2, color=COLORS["neutral"], va="top")
    ax.text(0.81, 0.54, "Raw-count\npseudobulk", fontsize=8.0, weight="bold", va="top", linespacing=1.05)
    ax.text(0.81, 0.34, "edgeR ranking +\nfixed-program GSEA", fontsize=6.9, va="top", linespacing=1.08)

    save_figure(fig, FIG_OUT, "fig1a_independent_paired_design")


def label_position(points: np.ndarray) -> np.ndarray:
    center = np.nanmedian(points, axis=0)
    distances = np.sqrt(((points - center) ** 2).sum(axis=1))
    keep = points[distances <= np.nanquantile(distances, 0.65)]
    return np.nanmedian(keep, axis=0) if len(keep) else center


def panel_identity_umap(dataset: str, stem: str) -> None:
    data = pd.read_csv(SOURCE_OUT / f"umap_identity_{dataset}.csv.gz")
    fig = new_figure(88, 78)
    ax = fig.add_axes([0.04, 0.04, 0.92, 0.76])
    counts = data["myeloid_identity"].value_counts()
    for identity in counts.sort_values(ascending=False).index:
        current = data[data["myeloid_identity"].eq(identity)]
        ax.scatter(
            current["UMAP1"],
            current["UMAP2"],
            s=0.42,
            color=IDENTITY_COLORS.get(identity, "#BDBDBD"),
            alpha=0.56,
            linewidths=0,
            rasterized=True,
        )
    for identity in counts.index:
        current = data[data["myeloid_identity"].eq(identity)][["UMAP1", "UMAP2"]].to_numpy()
        x, y = label_position(current)
        label = DISPLAY_IDENTITY.get(identity, identity)
        ax.text(
            x,
            y,
            label,
            fontsize=7.2,
            weight="bold",
            ha="center",
            va="center",
            color=IDENTITY_COLORS.get(identity, COLORS["text"]),
            path_effects=[pe.withStroke(linewidth=2.6, foreground="white")],
        )
    fig.text(0.06, 0.95, dataset, ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(
        0.06,
        0.875,
        f"{len(data):,} independently annotated myeloid cells",
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["neutral"],
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.margins(0.015)
    ax.set_aspect("equal", adjustable="box")
    save_figure(fig, FIG_OUT, stem)


def panel_marker_dotplot() -> None:
    data = pd.read_csv(SOURCE_OUT / "identity_marker_dotplot.csv")
    genes = [gene for block in MARKER_BLOCKS.values() for gene in block]
    identities = ["Microglia", "Macrophage", "Monocyte", "CDC"]
    data = data[data["gene"].isin(genes) & data["identity"].isin(identities)].copy()
    x_map = {gene: i for i, gene in enumerate(genes)}
    y_map = {identity: i for i, identity in enumerate(identities)}

    fig = new_figure(183, 77)
    axes = [fig.add_axes([0.10, 0.20, 0.33, 0.56]), fig.add_axes([0.51, 0.20, 0.33, 0.56])]
    norm = mpl.colors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    cmap = mpl.colormaps["RdBu_r"]
    for ax, dataset in zip(axes, DATASETS):
        current = data[data["dataset"].eq(dataset)]
        x = current["gene"].map(x_map).to_numpy()
        y = current["identity"].map(y_map).to_numpy()
        size = 5 + np.clip(current["pct_positive"].to_numpy(), 0, 100) * 0.44
        ax.scatter(
            x,
            y,
            s=size,
            c=current["scaled_expression"],
            cmap=cmap,
            norm=norm,
            edgecolor="#555555",
            linewidth=0.22,
        )
        ax.set_xlim(-0.55, len(genes) - 0.45)
        ax.set_ylim(len(identities) - 0.45, -0.55)
        ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(identities)), [DISPLAY_IDENTITY.get(x, x) for x in identities])
        ax.set_title(dataset, color=DATASET_COLORS[dataset], pad=7)
        for boundary in [2.5, 5.5, 8.5]:
            ax.axvline(boundary, color="#CFCFCF", lw=0.7)
        ax.set_axisbelow(True)
        ax.grid(color=COLORS["grid"], lw=0.45)
        for spine in ax.spines.values():
            spine.set_visible(False)
    cax = fig.add_axes([0.88, 0.45, 0.014, 0.25])
    scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_map.set_array(np.array([-2.0, 0.0, 2.0]))
    cb = fig.colorbar(scalar_map, cax=cax)
    cb.set_label("Scaled mean expression", fontsize=6.8)
    cb.ax.tick_params(labelsize=6.2, length=2)
    lax = fig.add_axes([0.86, 0.19, 0.12, 0.19])
    lax.axis("off")
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.text(0, 0.98, "Cells expressing", fontsize=6.8, va="top")
    for i, pct in enumerate([25, 50, 75]):
        y = 0.68 - i * 0.27
        lax.scatter(0.13, y, s=5 + pct * 0.44, facecolor="#BDBDBD", edgecolor="#555555", linewidth=0.25)
        lax.text(0.34, y, f"{pct}%", fontsize=6.4, va="center")
    fig.text(0.10, 0.95, "Canonical markers support independent myeloid annotation", fontsize=8.5, weight="bold", va="top")
    save_figure(fig, FIG_OUT, "fig1d_cross_cohort_identity_markers")


def panel_gsea(dataset: str, stem: str) -> None:
    curve = pd.read_csv(SOURCE_OUT / "raw20_gsea_curves.csv.gz")
    curve = curve[curve["dataset"].eq(dataset)].sort_values("rank")
    stats = pd.read_csv(SOURCE_OUT / "primary_gsea_statistics.csv")
    row = stats[stats["dataset"].eq(dataset)].iloc[0]
    nominal_p = float(row["pval"])
    color = DATASET_COLORS[dataset]

    fig = new_figure(88, 76)
    gs = fig.add_gridspec(3, 1, left=0.16, right=0.96, bottom=0.15, top=0.79, height_ratios=[3.5, 0.48, 0.95], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    hit_ax = fig.add_subplot(gs[1], sharex=ax)
    rank_ax = fig.add_subplot(gs[2], sharex=ax)

    rank = curve["rank"].to_numpy()
    running = curve["running_enrichment"].to_numpy()
    ax.plot(rank, running, color=color, lw=1.45)
    ax.fill_between(rank, 0, running, color=color, alpha=0.12, linewidth=0)
    ax.axhline(0, color="#999999", lw=0.55)
    ax.set_ylabel("Running enrichment")
    ax.set_xticks([])
    clean_axis(ax, keep_left=True, keep_bottom=False)
    ax.text(0.0, 1.035, dataset, transform=ax.transAxes, fontsize=8.5, weight="bold", ha="left", va="bottom")
    ax.text(
        1.0,
        1.035,
        f"NES {row.NES:.2f}   nominal P {nominal_p:.2g}   n={PAIR_COUNTS[dataset]} pairs",
        transform=ax.transAxes,
        fontsize=7,
        ha="right",
        va="bottom",
    )

    hits = curve.loc[curve["hit"].astype(bool), "rank"].to_numpy()
    hit_ax.vlines(hits, 0, 1, color="#333333", lw=0.48)
    hit_ax.set_ylim(0, 1)
    hit_ax.axis("off")

    rank_stat = curve["rank_stat"].to_numpy()
    rank_ax.fill_between(rank, 0, np.clip(rank_stat, 0, None), color=COLORS["recurrent"], alpha=0.75, linewidth=0)
    rank_ax.fill_between(rank, 0, np.clip(rank_stat, None, 0), color=COLORS["primary"], alpha=0.75, linewidth=0)
    rank_ax.axhline(0, color="#777777", lw=0.45)
    rank_ax.set_yticks([])
    rank_ax.set_xlabel("Ranked genes")
    rank_ax.text(0.0, -0.58, "Recurrent-enriched", transform=rank_ax.transAxes, ha="left", va="top", fontsize=6.8, color=COLORS["recurrent"])
    rank_ax.text(1.0, -0.58, "Primary-enriched", transform=rank_ax.transAxes, ha="right", va="top", fontsize=6.8, color=COLORS["primary"])
    clean_axis(rank_ax, keep_left=False, keep_bottom=True)
    rank_ax.spines["bottom"].set_color("#777777")
    fig.text(0.16, 0.94, "Miller microglial-inflammatory program", fontsize=7.2, color=COLORS["neutral"], va="top")
    save_figure(fig, FIG_OUT, stem)


def write_manifest() -> None:
    rows = [
        ("Figure1A", "fig1a_independent_paired_design", "Study design and patient-level inference"),
        ("Figure1B", "fig1b_gse174554_identity_umap", "GSE174554 independent myeloid identities"),
        ("Figure1C", "fig1c_gse274546_identity_umap", "GSE274546 independent myeloid identities"),
        ("Figure1D", "fig1d_cross_cohort_identity_markers", "Canonical identity marker expression"),
        ("Figure1E", "fig1e_gse174554_raw20_gsea", "GSE174554 paired pseudobulk GSEA"),
        ("Figure1F", "fig1f_gse274546_raw20_gsea", "GSE274546 paired pseudobulk GSEA"),
    ]
    manifest = pd.DataFrame(rows, columns=["panel", "stem", "message"])
    manifest["png"] = manifest["stem"].map(lambda x: str(FIG_OUT / f"{x}.png"))
    manifest["pdf"] = manifest["stem"].map(lambda x: str(FIG_OUT / f"{x}.pdf"))
    manifest.to_csv(WRITE_ROOT / "Figure1_panel_manifest.csv", index=False)


def main() -> None:
    np.random.seed(20260713)
    apply_publication_style()
    prepare_sources()
    panel_study_design()
    panel_identity_umap("GSE174554", "fig1b_gse174554_identity_umap")
    panel_identity_umap("GSE274546", "fig1c_gse274546_identity_umap")
    panel_marker_dotplot()
    panel_gsea("GSE174554", "fig1e_gse174554_raw20_gsea")
    panel_gsea("GSE274546", "fig1f_gse274546_raw20_gsea")
    write_manifest()
    print("STEP41_FIGURE1_COMPLETE panels=6")


if __name__ == "__main__":
    main()
