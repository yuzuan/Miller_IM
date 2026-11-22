#!/usr/bin/env python3
# Purpose: draw Figure 1 program panels from paired-pseudobulk ranks, GSEA and shared leading edges.
# Inputs: legacy Step 38 fixed-program provenance, targeted GSEA, all-gene edgeR statistics and gene directions;
# the Step 42 panel A PDF and source; and step1_primary_recurrence/camera_results.csv (primary camera P values).
# Outputs: _panels/45_figure1_program_rebuild/Figure1 contains panel PDFs/PNGs with stems
# Figure1A_compact_paired_design, Figure1B_GSE174554_raw20_gsea, Figure1C_GSE274546_raw20_gsea,
# Figure1D_cross_cohort_all_gene_rank_density, Figure1E_shared_leading_edge_heatmap and Figure1_candidate_contact_sheet.
# The matching _figure_build directory contains Figure1_panel_manifest.csv and source_data/ with the panel sources
# (.csv.gz for B and C, .csv for A, D and E) and Figure1_candidate_contact_sheet_layout.csv.
# Run: python3 NOA_R1/scripts/06b_figure1_program_panels.py
# Keep Matplotlib and Pillow to retain the original graphics and panel placement.

from __future__ import annotations

import importlib
import shutil
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw
from scipy.stats import spearmanr

figure_style = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
legacy_project_root = figure_style.legacy_root(revision_root)

# The Figure 1 assembler redirects these original output attributes.
write_root = revision_root / "_figure_build/45_figure1_program_rebuild/Figure1"
SOURCE_OUT = write_root / "source_data"
FIG_OUT = revision_root / "_panels/45_figure1_program_rebuild/Figure1"

# Both ranking and leading-edge comparisons use the same independent cohorts and fixed signature.
datasets = ["GSE174554", "GSE274546"]
raw20_name = "Miller_Microglial_Inflammatory_raw_top20"


def format_p(value: float) -> str:
    """Format nominal fgsea P values without implying multiplicity correction."""
    if value < 1e-4:
        formatted = "<1e-4"
    elif value < 0.01:
        formatted = f"{value:.4f}".rstrip("0").rstrip(".")
    else:
        formatted = f"{value:.3f}".rstrip("0").rstrip(".")
    return formatted


def load_raw20_genes(step38_root: Path) -> list[str]:
    """Return the 20 genes in the frozen signature's published order."""
    provenance_table = pd.read_csv(step38_root / "fixed_program_provenance.csv")
    rows = provenance_table.loc[provenance_table["signature"].eq(raw20_name), "genes"]
    if len(rows) != 1:
        raise ValueError(f"Expected one {raw20_name} provenance row, found {len(rows)}")
    genes = str(rows.iloc[0]).split(";")
    if len(genes) != 20:
        raise ValueError(f"Expected 20 MG-Inflam20 genes, found {len(genes)}")
    return genes


def enrichment_curve(deg_table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Return running enrichment scores along the original signed square-root F rank."""
    rank_table = deg_table.copy()
    rank_table["rank_stat"] = np.sign(rank_table["logFC"]) * np.sqrt(np.maximum(rank_table["F"], 0))
    rank_table = (
        rank_table.loc[np.isfinite(rank_table["rank_stat"])]
        .sort_values("rank_stat", ascending=False)
        .reset_index(drop=True)
    )
    if rank_table["gene"].duplicated().any():
        duplicates = rank_table.loc[rank_table["gene"].duplicated(), "gene"].tolist()
        raise ValueError(f"Duplicate genes in Step38 ranked table: {duplicates[:10]}")
    hits = rank_table["gene"].isin(genes).to_numpy()
    weights = np.abs(rank_table["rank_stat"].to_numpy(dtype=float))
    hit_weights = np.where(hits, weights, 0.0)
    if hit_weights.sum() <= 0 or (~hits).sum() == 0:
        raise ValueError("Cannot construct enrichment curve from the Step38 ranked table")
    increments = np.where(
        hits,
        hit_weights / hit_weights.sum(),
        -1.0 / (~hits).sum(),
    )
    result_table = pd.DataFrame(
        {
            "rank": np.arange(1, len(rank_table) + 1),
            "gene": rank_table["gene"].to_numpy(),
            "rank_stat": rank_table["rank_stat"].to_numpy(),
            "hit": hits,
            "running_enrichment": np.cumsum(increments),
        }
    )
    return result_table


def copy_locked_panel() -> dict[str, str]:
    """Copy the legacy A source and PDF and return its panel manifest row."""
    stem = "Figure1A_compact_paired_design"
    source = SOURCE_OUT / f"{stem}_source.csv"
    png = FIG_OUT / f"{stem}.png"
    pdf = FIG_OUT / f"{stem}.pdf"
    step42_source = legacy_project_root / "write/42_visual_story_rebuild/Figure1/source_data"
    step42_figures = legacy_project_root / "figures/42_visual_story_rebuild/Figure1"
    shutil.copy2(step42_source / source.name, source)
    shutil.copy2(step42_figures / pdf.name, pdf)
    figure_style.render_pdf_png(pdf, png)
    manifest_row = panel_record("Figure1A", stem, "Locked independent paired-cohort study design", source, pdf, png)
    return manifest_row


def load_gsea_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the ranked curves and one formal GSEA result per paired cohort."""
    step38_root = legacy_project_root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    statistics_table = pd.read_csv(step38_root / "independent_fixed_program_targeted_gsea.csv")
    is_formal = statistics_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
    statistics_table = statistics_table.loc[
        statistics_table["dataset"].isin(datasets)
        & statistics_table["threshold"].eq(20)
        & statistics_table["pathway"].eq(raw20_name)
        & is_formal
    ].copy()
    if len(statistics_table) != 2:
        raise ValueError(f"Expected two formal raw20 GSEA rows, found {len(statistics_table)}")
    pair_counts = statistics_table.set_index("dataset")["n_pairs"].astype(int).to_dict()
    if pair_counts != {"GSE174554": 17, "GSE274546": 45}:
        raise ValueError(f"Unexpected formal pair counts in Step38 GSEA: {pair_counts}")

    deg_table = pd.read_csv(step38_root / "independent_paired_edger_all_genes.csv")
    genes = load_raw20_genes(step38_root)
    cohort_curve_tables = []
    for dataset in datasets:
        cohort_deg_table = deg_table.loc[deg_table["dataset"].eq(dataset) & deg_table["threshold"].eq(20)].copy()
        if "formal_testable" in cohort_deg_table.columns:
            current_is_formal = cohort_deg_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
            cohort_deg_table = cohort_deg_table.loc[current_is_formal].copy()
        if cohort_deg_table.empty:
            raise ValueError(f"No formal Step38 differential-expression rows for {dataset}")
        curve_table = enrichment_curve(cohort_deg_table, genes)
        curve_table["dataset"] = dataset
        cohort_curve_tables.append(curve_table)

    curves_table_total = pd.concat(cohort_curve_tables, ignore_index=True)
    outputs = (curves_table_total, statistics_table)
    return outputs


def panel_gsea(
    dataset: str,
    panel: str,
    curves_table_total: pd.DataFrame,
    statistics_table: pd.DataFrame,
    running_limits: tuple[float, float],
    rank_limit: float,
) -> dict[str, str]:
    """Draw one cohort's preranked GSEA and return its panel manifest row."""
    stem = f"{panel}_{dataset}_raw20_gsea"
    curve_table = curves_table_total.loc[curves_table_total["dataset"].eq(dataset)].sort_values("rank").copy()
    row = statistics_table.loc[statistics_table["dataset"].eq(dataset)].iloc[0]
    curve_table["NES"] = float(row["NES"])
    curve_table["nominal_P"] = float(row["pval"])
    curve_table["n_pairs"] = int(row["n_pairs"])
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    curve_table.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    fig = figure_style.new_figure(86, 53)
    grid = fig.add_gridspec(
        3,
        1,
        left   = 0.16,
        right  = 0.96,
        bottom = 0.17,
        top    = 0.80,
        height_ratios = [3.6, 0.45, 1.0],
        hspace = 0.07,
    )
    enrichment_ax = fig.add_subplot(grid[0])
    hit_ax = fig.add_subplot(grid[1], sharex=enrichment_ax)
    rank_ax = fig.add_subplot(grid[2], sharex=enrichment_ax)

    # These three tracks show the running score, program membership and signed rank statistic.
    ranks = curve_table["rank"].to_numpy(dtype=float)
    running = curve_table["running_enrichment"].to_numpy(dtype=float)
    rank_stat = curve_table["rank_stat"].to_numpy(dtype=float)
    color = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}[dataset]

    enrichment_ax.plot(ranks, running, color=color, linewidth=1.55)
    enrichment_ax.fill_between(ranks, 0, running, color=color, alpha=0.14, linewidth=0)
    enrichment_ax.axhline(0, color="#999999", linewidth=0.55)
    enrichment_ax.set_ylim(*running_limits)
    enrichment_ax.set_ylabel("Running ES")
    enrichment_ax.set_xticks([])
    figure_style.clean_axis(enrichment_ax, keep_left=True, keep_bottom=False)
    enrichment_ax.text(
        0.0,
        1.05,
        dataset,
        transform  = enrichment_ax.transAxes,
        fontsize   = 8.7,
        fontweight = "bold",
        color      = color,
        ha         = "left",
        va         = "bottom",
    )

    # Use the primary camera results and the original October display precision in each cohort.
    camera_table = pd.read_csv(revision_root / "step1_primary_recurrence/camera_results.csv")
    primary_result = camera_table["dataset"].eq(dataset) & camera_table["role"].eq("primary")
    camera_p = float(camera_table.loc[primary_result, "p_value"].item())
    camera_annotation = f"{camera_p:.2g}" if dataset == "GSE174554" else f"{camera_p:.3g}"
    enrichment_ax.text(
        1.0,
        1.05,
        f"NES {float(row['NES']):.2f}   P {format_p(float(row['pval']))}\ncamera P {camera_annotation}",
        transform   = enrichment_ax.transAxes,
        fontsize    = 6.6,
        color       = "#555555",
        ha          = "right",
        va          = "bottom",
        linespacing = 1.15,
    )

    hits = curve_table.loc[curve_table["hit"].astype(bool), "rank"].to_numpy(dtype=float)
    hit_ax.vlines(hits, 0, 1, color="#333333", linewidth=0.48)
    hit_ax.set_ylim(0, 1)
    hit_ax.axis("off")

    rank_ax.fill_between(
        ranks, 0, np.clip(rank_stat, 0, None), color=figure_style.colors["recurrent"], alpha=0.78, linewidth=0
    )
    rank_ax.fill_between(
        ranks, 0, np.clip(rank_stat, None, 0), color=figure_style.colors["primary"], alpha=0.78, linewidth=0
    )
    rank_ax.axhline(0, color="#777777", linewidth=0.45)
    rank_ax.set_ylim(-rank_limit, rank_limit)
    rank_ax.set_yticks([])
    rank_ax.set_xlabel("Ranked transcriptome")
    rank_ax.text(
        0.0,
        -0.62,
        "Recurrent-enriched",
        transform = rank_ax.transAxes,
        ha        = "left",
        va        = "top",
        fontsize  = 6.5,
        color     = figure_style.colors["recurrent"],
    )
    rank_ax.text(
        1.0,
        -0.62,
        "Primary-enriched",
        transform = rank_ax.transAxes,
        ha        = "right",
        va        = "top",
        fontsize  = 6.5,
        color     = figure_style.colors["primary"],
    )
    figure_style.clean_axis(rank_ax, keep_left=False, keep_bottom=True)
    rank_ax.spines["bottom"].set_color("#777777")
    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step45 recurrent GBM Figure 1 program rebuild")
    manifest_row = panel_record(panel, stem, f"{dataset} paired pseudobulk raw20 GSEA", source, pdf, png)
    return manifest_row


def panel_rank_rank(curves_table_total: pd.DataFrame, statistics_table: pd.DataFrame) -> dict[str, str]:
    """Draw all shared gene ranks and identify the shared leading-edge genes."""
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    source = SOURCE_OUT / f"{stem}_source.csv"
    gse174554_table = curves_table_total.loc[
        curves_table_total["dataset"].eq("GSE174554"),
        ["gene", "rank", "rank_stat", "hit"],
    ].rename(
        columns={
            "rank": "rank_GSE174554",
            "rank_stat": "rank_stat_GSE174554",
            "hit": "program_GSE174554",
        }
    )
    gse274546_table = curves_table_total.loc[
        curves_table_total["dataset"].eq("GSE274546"),
        ["gene", "rank", "rank_stat", "hit"],
    ].rename(
        columns={
            "rank": "rank_GSE274546",
            "rank_stat": "rank_stat_GSE274546",
            "hit": "program_GSE274546",
        }
    )
    rank_table_total = gse174554_table.merge(gse274546_table, on="gene", how="inner", validate="one_to_one")
    rank_table_total["miller_raw20"] = rank_table_total[["program_GSE174554", "program_GSE274546"]].any(axis=1)
    evidence_layout = importlib.import_module("06a_figure1_evidence_layout")
    rank_table_total["density"] = evidence_layout.rank_density(rank_table_total)
    direction_table = shared_leading_edge_table(statistics_table)
    shared_genes = set(direction_table["gene"])
    rank_table_total["shared_leading_edge"] = rank_table_total["gene"].isin(shared_genes)
    rank_table_total["labelled"] = rank_table_total["shared_leading_edge"]

    background_table = rank_table_total.loc[~rank_table_total["miller_raw20"]].sort_values("density")
    raw20_other_table = rank_table_total.loc[
        rank_table_total["miller_raw20"] & ~rank_table_total["shared_leading_edge"]
    ]
    shared_table = rank_table_total.loc[rank_table_total["shared_leading_edge"]]
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            "Expected 7 shared and 12 other measured MG-Inflam20 genes, "
            f"found {len(shared_table)} and {len(raw20_other_table)}"
        )
    rank_table_total.to_csv(source, index=False)

    fig = figure_style.new_figure(92, 83)
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
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
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s         = 21,
        facecolor = "white",
        edgecolor = "#D8574F",
        linewidth = 0.85,
        zorder    = 3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s         = 27,
        color     = "#D8574F",
        edgecolor = "white",
        linewidth = 0.65,
        zorder    = 4,
    )
    label_offsets = {
        "CCL4": (-25, 7),
        "CH25H": (5, 5),
        "SGK1": (-23, -10),
        "FOLR2": (5, 7),
        "PDK4": (5, -10),
        "KLF6": (-18, -11),
        "SIGLEC8": (5, 4),
    }
    for _, row in shared_table.iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext       = label_offsets.get(row["gene"], (3, 3)),
            textcoords   = "offset points",
            fontsize     = 5.5,
            color        = "#A63D38",
            path_effects = [pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = float(np.nanmax(np.abs(rank_table_total[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    figure_style.clean_axis(ax)
    rho, _ = spearmanr(rank_table_total["rank_stat_GSE174554"], rank_table_total["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(rank_table_total):,}\nSpearman rho = {rho:.2f}",
        transform = ax.transAxes,
        va        = "top",
        fontsize  = 6.5,
        color     = "#717171",
    )
    ax.text(
        0.98,
        0.03,
        f"{len(shared_table)} shared leading-edge genes\n{len(raw20_other_table)} other measured raw20 genes",
        transform   = ax.transAxes,
        ha          = "right",
        va          = "bottom",
        fontsize    = 6.1,
        color       = "#A63D38",
        linespacing = 1.15,
    )
    fig.text(0.16, 0.965, "Cross-cohort gene-rank comparison", fontsize=8.7, fontweight="bold", va="top")
    fig.text(
        0.16,
        0.915,
        "Transcriptome background with shared raw20 leading-edge genes highlighted",
        fontsize = 6.2,
        color    = "#717171",
        va       = "top",
    )
    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step45 recurrent GBM Figure 1 program rebuild")
    manifest_row = panel_record(
        "Figure1D",
        stem,
        "All-gene rank-stat comparison with seven shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )
    return manifest_row


def shared_leading_edge_table(statistics_table: pd.DataFrame) -> pd.DataFrame:
    """Return the gene-by-cohort effects for the intersection of the two GSEA leading edges."""
    leading_sets = {row["dataset"]: set(str(row["leadingEdge"]).split(";")) for _, row in statistics_table.iterrows()}
    shared = set.intersection(*(leading_sets[dataset] for dataset in datasets))
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}
    if shared != expected_shared_leading_edge:
        raise ValueError(
            "Unexpected shared leading-edge genes: "
            f"expected {sorted(expected_shared_leading_edge)}, found {sorted(shared)}"
        )

    # Keep the fixed signature's gene order and evaluate the unrounded single-gene FDR values.
    step38_root = legacy_project_root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    raw20_genes = load_raw20_genes(step38_root)
    order = [gene for gene in raw20_genes if gene in shared]

    direction_table = pd.read_csv(step38_root / "independent_fixed_program_gene_direction.csv")
    direction_table = direction_table.loc[
        direction_table["dataset"].isin(datasets)
        & direction_table["threshold"].eq(20)
        & direction_table["signature"].eq(raw20_name)
        & direction_table["gene"].isin(shared)
    ].copy()
    if len(direction_table) != len(shared) * len(datasets):
        raise ValueError(f"Expected one row per gene and cohort, found {len(direction_table)} rows")
    if set(direction_table["gene"]) != expected_shared_leading_edge:
        raise ValueError(f"Unexpected genes in shared direction table: {sorted(direction_table['gene'].unique())}")
    per_gene_counts = direction_table.groupby("gene")["dataset"].nunique().astype(int)
    if not per_gene_counts.eq(len(datasets)).all():
        raise ValueError(f"Each shared gene must be tested in both cohorts: {per_gene_counts.to_dict()}")
    if int((direction_table["FDR"] < 0.05).sum()) != 0:
        raise ValueError("Expected 0/14 shared gene-by-cohort tests at FDR < 0.05")
    direction_table["shared_leading_edge"] = True
    direction_table["gene_order"] = direction_table["gene"].map({gene: idx for idx, gene in enumerate(order)})
    direction_table["dataset_order"] = direction_table["dataset"].map(
        {dataset: idx for idx, dataset in enumerate(datasets)}
    )
    direction_table = direction_table.sort_values(["gene_order", "dataset_order"]).reset_index(drop=True)
    return direction_table


def panel_shared_leading_edge_heatmap(statistics_table: pd.DataFrame) -> dict[str, str]:
    """Draw the shared leading-edge effects on the original diverging scale."""
    stem = "Figure1E_shared_leading_edge_heatmap"
    direction_table = shared_leading_edge_table(statistics_table)
    source = SOURCE_OUT / f"{stem}_source.csv"
    direction_table.to_csv(source, index=False)

    order = direction_table.sort_values("gene_order")["gene"].drop_duplicates().tolist()
    effect_table = direction_table.pivot(index="gene", columns="dataset", values="logFC").reindex(
        index=order, columns=datasets
    )
    max_abs = float(np.nanmax(np.abs(effect_table.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    fig = figure_style.new_figure(86, 80)
    ax = fig.add_axes([0.24, 0.22, 0.55, 0.56])
    image = ax.imshow(effect_table.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest")
    indexed_gsea_table = statistics_table.set_index("dataset")
    ax.set_xticks(
        range(len(datasets)),
        [
            f"GSE174554\n{int(indexed_gsea_table.loc['GSE174554', 'n_pairs'])} pairs",
            f"GSE274546\n{int(indexed_gsea_table.loc['GSE274546', 'n_pairs'])} pairs",
        ],
        fontsize=6.8,
    )
    ax.set_yticks(range(len(order)), order, fontsize=7.1)
    ax.tick_params(length=0)
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(datasets):
            value = float(effect_table.loc[gene, dataset])
            text_color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=text_color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(datasets), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    color_ax = fig.add_axes([0.82, 0.30, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)
    colorbar.outline.set_linewidth(0.45)

    fig.text(0.24, 0.95, "Shared leading-edge genes", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.24,
        0.90,
        "Intersection of the two formal raw20 GSEA leading edges",
        fontsize = 6.2,
        color    = "#717171",
        va       = "top",
    )
    fig.text(
        0.24,
        0.145,
        f"Cohort GSEA: NES {indexed_gsea_table.loc['GSE174554', 'NES']:.2f} / "
        f"{indexed_gsea_table.loc['GSE274546', 'NES']:.2f}\n"
        f"Nominal P {format_p(float(indexed_gsea_table.loc['GSE174554', 'pval']))} / "
        f"{format_p(float(indexed_gsea_table.loc['GSE274546', 'pval']))}",
        fontsize    = 6.1,
        color       = "#555555",
        va          = "top",
        linespacing = 1.25,
    )
    n_significant = int((direction_table["FDR"] < 0.05).sum())
    fig.text(
        0.24,
        0.065,
        f"{n_significant}/{len(direction_table)} individual gene tests reached FDR < 0.05",
        fontsize = 5.9,
        color    = "#777777",
        va       = "top",
    )

    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step45 recurrent GBM Figure 1 program rebuild")
    manifest_row = panel_record(
        "Figure1E", stem, "Seven shared leading-edge genes across the two recurrence cohorts", source, pdf, png
    )
    return manifest_row


def make_contact_sheet(records: list[dict[str, str]]) -> tuple[Path, Path, Path]:
    """Assemble the original candidate layout and return its PDF, PNG and layout CSV."""
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
    label_font = figure_style.image_font(52, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        figure_style.paste_contained(
            canvas,
            panel_paths[panel],
            (int(row["left"]), int(row["top"]), int(row["right"]), int(row["bottom"])),
        )
        drawer.text((int(row["left"]), int(row["top"]) - 56), panel[-1], fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1_candidate_contact_sheet.png"
    pdf = FIG_OUT / "Figure1_candidate_contact_sheet.pdf"

    # Export the candidate preview at its original page size and 300 dpi.
    contact_fig = plt.figure(figsize=(12, 4960 / 300), facecolor="white")
    contact_ax = contact_fig.add_axes([0, 0, 1, 1])
    contact_ax.imshow(np.asarray(canvas), interpolation="nearest")
    contact_ax.axis("off")
    contact_fig.canvas.draw()
    contact_fig.savefig(
        pdf,
        dpi = 300,
        facecolor = "white",
        edgecolor = "none",
        metadata={
            "Creator": "Step45 recurrent GBM Figure 1 program rebuild",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(contact_fig)
    figure_style.render_pdf_png(pdf, png, dpi=300)

    layout_table = pd.DataFrame(layout_rows)
    source = SOURCE_OUT / "Figure1_candidate_contact_sheet_layout.csv"
    layout_table.to_csv(source, index=False)
    outputs = (pdf, png, source)
    return outputs


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, str]:
    """Return the existing Figure 1 candidate manifest fields and file digests."""
    manifest_row = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": figure_style.sha256(source),
        "pdf_sha256": figure_style.sha256(pdf),
        "png_sha256": figure_style.sha256(png),
    }
    return manifest_row


def main() -> None:
    """Build the Figure 1 candidates and write their manifest. Returns nothing."""
    np.random.seed(20260713)
    figure_style.apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    curves_table_total, statistics_table = load_gsea_inputs()
    running_min = float(curves_table_total["running_enrichment"].min())
    running_max = float(curves_table_total["running_enrichment"].max())
    running_pad = 0.06 * (running_max - running_min)
    running_limits = (running_min - running_pad, running_max + running_pad)
    rank_limit = float(np.nanmax(np.abs(curves_table_total["rank_stat"].to_numpy()))) * 1.04

    records = [
        copy_locked_panel(),
        panel_gsea("GSE174554", "Figure1B", curves_table_total, statistics_table, running_limits, rank_limit),
        panel_gsea("GSE274546", "Figure1C", curves_table_total, statistics_table, running_limits, rank_limit),
        panel_rank_rank(curves_table_total, statistics_table),
        panel_shared_leading_edge_heatmap(statistics_table),
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
    manifest_table = pd.DataFrame(records)
    manifest_table.to_csv(write_root / "Figure1_panel_manifest.csv", index=False)
    print(f"Figure 1 program panels: {int(manifest_table['panel'].ne('Figure1_preview').sum())}")


if __name__ == "__main__":
    main()
