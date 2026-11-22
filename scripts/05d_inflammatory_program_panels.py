#!/usr/bin/env python3
# Purpose: Rebuild Figure 1 GSEA, gene ranks and shared leading-edge effects from the formal Step 38 results.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05d_inflammatory_program_panels.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
import shutil
from pathlib import Path
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw
from scipy.stats import spearmanr
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
colors = figure_style.colors
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
render_pdf_png = figure_style.render_pdf_png
image_font = figure_style.image_font
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel
panel_record = figure_style.panel_record
enrichment_curve = figure_style.enrichment_curve
rank_density = figure_style.rank_density

if __name__ == "__main__":

    # Keep the seed at this point to reproduce the sampling or jitter sequence.
    np.random.seed(20260713)
    apply_publication_style()
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write/45_figure1_program_rebuild/Figure1"
    source_out = write_root / "source_data"

    # Use the published output directory and basename.
    source_out.mkdir(parents=True, exist_ok=True)
    fig_out = root / "figures/45_figure1_program_rebuild/Figure1"
    fig_out.mkdir(parents=True, exist_ok=True)
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    statistics_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")

    is_formal = statistics_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
    datasets = ["GSE174554", "GSE274546"]
    raw20_name = "Miller_Microglial_Inflammatory_raw_top20"
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
    deg_table = pd.read_csv(step38 / "independent_paired_edger_all_genes.csv")
    provenance_table = pd.read_csv(step38 / "fixed_program_provenance.csv")
    rows = provenance_table.loc[provenance_table["signature"].eq(raw20_name), "genes"]

    # Use the existing input checks and branch conditions.
    if len(rows) != 1:
        raise ValueError(f"Expected one {raw20_name} provenance row, found {len(rows)}")
    genes = str(rows.iloc[0]).split(";")
    if len(genes) != 20:
        raise ValueError(f"Expected 20 Miller-IM genes, found {len(genes)}")
    curves = []


    # Process cohorts, genes or panels in the published order.
    for dataset in datasets:
        current_table = deg_table.loc[deg_table["dataset"].eq(dataset) & deg_table["threshold"].eq(20)].copy()
        if "formal_testable" in current_table.columns:
            current_is_formal = (
                current_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
            )
            current_table = current_table.loc[current_is_formal].copy()
        if current_table.empty:
            raise ValueError(f"No formal Step38 differential-expression rows for {dataset}")
        curve_table = enrichment_curve(current_table, genes)
        curve_table["dataset"] = dataset

        curves.append(curve_table)
    curve_table_total = pd.concat(curves, ignore_index=True)
    running_min = float(curve_table_total["running_enrichment"].min())

    running_max = float(curve_table_total["running_enrichment"].max())
    running_pad = 0.06 * (running_max - running_min)
    running_limits = (running_min - running_pad, running_max + running_pad)
    rank_limit = float(np.nanmax(np.abs(curve_table_total["rank_stat"].to_numpy()))) * 1.04
    stem = "Figure1A_compact_paired_design"

    # Use the published output directory and basename.
    fig_out.mkdir(parents=True, exist_ok=True)
    source_out.mkdir(parents=True, exist_ok=True)
    source = source_out / f"{stem}_source.csv"
    png = fig_out / f"{stem}.png"
    pdf = fig_out / f"{stem}.pdf"

    # Use the published output directory and basename.
    step42_write = root / "write/42_visual_story_rebuild/Figure1"
    shutil.copy2(step42_write / "source_data" / f"{stem}_source.csv", source)
    step42_fig = root / "figures/42_visual_story_rebuild/Figure1"
    shutil.copy2(step42_fig / f"{stem}.pdf", pdf)
    render_pdf_png(pdf, png)

    # Record the panel source table, PDF and PNG in the original manifest field order.
    panel_record_1 = panel_record(
        "Figure1A", stem, "Locked independent paired-cohort study design", source, pdf, png
    )

    def format_p(value: float) -> str:
        """Return the display string for a raw P value using the original thresholds and precision."""

        # Use the existing input checks and branch conditions.
        if value < 0.0001:
            result = "<1e-4"
            return result
        if value < 0.01:
            result = f"{value:.4f}".rstrip("0").rstrip(".")
            return result
        result = f"{value:.3f}".rstrip("0").rstrip(".")
        return result

    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    def panel_gsea(
        dataset: str,
        panel: str,
        curve_table_total: pd.DataFrame,
        statistics_table: pd.DataFrame,
        running_limits: tuple[float, float],
        rank_limit: float,
    ) -> dict[str, str]:
        """Save the cohort three-track GSEA panel under the given stem and return its manifest record."""

        # Use the panel's published stem for its source table and figure files.
        stem = f"{panel}_{dataset}_raw20_gsea"
        curve_table = (
            curve_table_total.loc[curve_table_total["dataset"].eq(dataset)].sort_values("rank").copy()
        )
        row = statistics_table.loc[statistics_table["dataset"].eq(dataset)].iloc[0]
        curve_table["NES"] = float(row["NES"])
        curve_table["nominal_P"] = float(row["pval"])

        # Prepare the gene ranks and hits for the formal cohorts at the published threshold of 20 cells.
        curve_table["n_pairs"] = int(row["n_pairs"])
        source = source_out / f"{stem}_source.csv.gz"
        csv_path = source
        curve_table.to_csv(csv_path, index=False, compression={"method": "gzip", "mtime": 0})
        if not csv_path.is_file() or csv_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty output: {csv_path.name}")

        # Draw the panel with the published canvas dimensions, placements and styling.
        fig = new_figure(86, 74)
        grid = fig.add_gridspec(
            3, 1, left=0.16, right=0.96, bottom=0.17, top=0.86, height_ratios=[3.6, 0.45, 1.0], hspace=0.07
        )
        enrichment_ax = fig.add_subplot(grid[0])
        hit_ax = fig.add_subplot(grid[1], sharex=enrichment_ax)
        rank_ax = fig.add_subplot(grid[2], sharex=enrichment_ax)

        # Prepare the cohort gene ranks using the original order and merge keys.
        ranks = curve_table["rank"].to_numpy(dtype=float)
        running = curve_table["running_enrichment"].to_numpy(dtype=float)
        rank_stat = curve_table["rank_stat"].to_numpy(dtype=float)
        color = dataset_colors[dataset]
        enrichment_ax.plot(ranks, running, color=color, linewidth=1.55)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
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
            f"NES {float(row['NES']):.2f}   P {format_p(float(row['pval']))}\n{int(row['n_pairs'])} paired patients",
            transform=enrichment_ax.transAxes,
            fontsize=6.6,
            color="#555555",
            ha="right",
            va="bottom",
            linespacing=1.15,
        )
        hits = curve_table.loc[curve_table["hit"].astype(bool), "rank"].to_numpy(dtype=float)
        hit_ax.vlines(hits, 0, 1, color="#333333", linewidth=0.48)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        hit_ax.set_ylim(0, 1)
        hit_ax.axis("off")
        rank_ax.fill_between(
            ranks, 0, np.clip(rank_stat, 0, None), color=colors["recurrent"], alpha=0.78, linewidth=0
        )
        rank_ax.fill_between(
            ranks, 0, np.clip(rank_stat, None, 0), color=colors["primary"], alpha=0.78, linewidth=0
        )
        rank_ax.axhline(0, color="#777777", linewidth=0.45)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
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
            color=colors["recurrent"],
        )
        rank_ax.text(
            1.0,
            -0.62,
            "Primary-enriched",
            transform=rank_ax.transAxes,
            ha="right",
            va="top",
            fontsize=6.5,
            color=colors["primary"],
        )

        clean_axis(rank_ax, keep_left=False, keep_bottom=True)
        rank_ax.spines["bottom"].set_color("#777777")
        pdf, png = save_panel(fig, fig_out, stem, "Step45 recurrent GBM Figure 1 program rebuild")
        result = panel_record(panel, stem, f"{dataset} paired pseudobulk raw20 GSEA", source, pdf, png)
        return result

    # Record the panel source table, PDF and PNG in the original manifest field order.
    panel_record_2 = panel_gsea(
        "GSE174554", "Figure1B", curve_table_total, statistics_table, running_limits, rank_limit
    )
    panel_record_3 = panel_gsea(
        "GSE274546", "Figure1C", curve_table_total, statistics_table, running_limits, rank_limit
    )
    left_table = curve_table_total.loc[
        curve_table_total["dataset"].eq("GSE174554"), ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE174554", "rank_stat": "rank_stat_GSE174554", "hit": "program_GSE174554"}
    )
    right_table = curve_table_total.loc[
        curve_table_total["dataset"].eq("GSE274546"), ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE274546", "rank_stat": "rank_stat_GSE274546", "hit": "program_GSE274546"}
    )
    rank_table_total = left_table.merge(right_table, on="gene", how="inner", validate="one_to_one")
    rank_table_total["miller_raw20"] = rank_table_total[["program_GSE174554", "program_GSE274546"]].any(
        axis=1
    )


    # Prepare the cohort gene ranks using the original order and merge keys.
    rank_table_total["density"] = rank_density(rank_table_total)
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}

    def shared_leading_edge_table(statistics_table: pd.DataFrame) -> pd.DataFrame:
        """Use formal GSEA statistics and fixed gene definitions to return shared leading-edge effects by cohort."""

        leading_sets = {
            row["dataset"]: set(str(row["leadingEdge"]).split(";")) for _, row in statistics_table.iterrows()
        }
        shared = set.intersection(*(leading_sets[dataset] for dataset in datasets))
        if shared != expected_shared_leading_edge:
            raise ValueError(
                (
                    "Unexpected shared leading-edge genes: expected "
                    f"{sorted(expected_shared_leading_edge)}"
                    ", found "
                    f"{sorted(shared)}"
                )
            )
        provenance_table = pd.read_csv(step38 / "fixed_program_provenance.csv")
        raw20_genes = (
            provenance_table.loc[provenance_table["signature"].eq(raw20_name), "genes"].iloc[0].split(";")
        )

        order = [gene for gene in raw20_genes if gene in shared]
        direction_table = pd.read_csv(step38 / "independent_fixed_program_gene_direction.csv")
        direction_table = direction_table.loc[
            direction_table["dataset"].isin(datasets)
            & direction_table["threshold"].eq(20)
            & direction_table["signature"].eq(raw20_name)
            & direction_table["gene"].isin(shared)
        ].copy()
        if len(direction_table) != 14:
            raise ValueError(f"Expected fourteen gene-by-cohort rows, found {len(direction_table)}")
        if set(direction_table["gene"]) != expected_shared_leading_edge:
            raise ValueError(
                f"Unexpected genes in shared direction table: {sorted(direction_table['gene'].unique())}"
            )

        per_gene_counts = direction_table.groupby("gene")["dataset"].nunique().astype(int)
        if not per_gene_counts.eq(2).all():
            raise ValueError(f"Each shared gene must be tested in both cohorts: {per_gene_counts.to_dict()}")
        if int((direction_table["FDR"] < 0.05).sum()) != 0:
            raise ValueError("Expected 0/14 shared gene-by-cohort tests at FDR < 0.05")
        direction_table["shared_leading_edge"] = True
        direction_table["gene_order"] = direction_table["gene"].map(
            {gene: idx for idx, gene in enumerate(order)}
        )

        # Use the existing gene directions and effects from the fitted differential-expression models.
        direction_table["dataset_order"] = direction_table["dataset"].map(
            {dataset: idx for idx, dataset in enumerate(datasets)}
        )
        direction_table = direction_table.sort_values(["gene_order", "dataset_order"]).reset_index(drop=True)
        return direction_table

    shared_direction_table = shared_leading_edge_table(statistics_table)
    shared_genes = set(shared_direction_table["gene"])
    rank_table_total["shared_leading_edge"] = rank_table_total["gene"].isin(shared_genes)

    # Prepare the cohort gene ranks using the original order and merge keys.
    rank_table_total["labelled"] = rank_table_total["shared_leading_edge"]
    background_table = rank_table_total.loc[~rank_table_total["miller_raw20"]].sort_values("density")
    raw20_other_table = rank_table_total.loc[
        rank_table_total["miller_raw20"] & ~rank_table_total["shared_leading_edge"]
    ]
    shared_table = rank_table_total.loc[rank_table_total["shared_leading_edge"]]
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            (
                "Expected 7 shared and 12 other measured Miller-IM genes, found "
                f"{len(shared_table)}"
                " and "
                f"{len(raw20_other_table)}"
            )
        )

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    source = source_out / f"{stem}_source.csv"
    csv_path = source
    rank_table_total.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    fig = new_figure(92, 83)
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])

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
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s=21,
        facecolor="white",
        edgecolor="#D8574F",
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s=27,
        color="#D8574F",
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
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
            xytext=label_offsets.get(row["gene"], (3, 3)),
            textcoords="offset points",
            fontsize=5.5,
            color="#A63D38",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
        )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = (
        float(np.nanmax(np.abs(rank_table_total[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy())))
        * 1.04
    )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    clean_axis(ax)

    rho, _ = spearmanr(rank_table_total["rank_stat_GSE174554"], rank_table_total["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(rank_table_total):,}\nSpearman rho = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.5,
        color="#717171",
    )
    ax.text(
        0.98,
        0.03,
        f"{len(shared_table)} shared leading-edge genes\n{len(raw20_other_table)} other measured raw20 genes",
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

    # Save the PDF master at the published size and render its PNG from that PDF.
    pdf, png = save_panel(fig, fig_out, stem, "Step45 recurrent GBM Figure 1 program rebuild")
    panel_record_4 = panel_record(
        "Figure1D",
        stem,
        "All-gene rank-stat comparison with seven shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )
    stem = "Figure1E_shared_leading_edge_heatmap"
    direction_table = shared_leading_edge_table(statistics_table)

    # Use the published output directory and basename.
    source = source_out / f"{stem}_source.csv"
    csv_path = source
    direction_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    order = direction_table.sort_values("gene_order")["gene"].drop_duplicates().tolist()

    # Prepare the published gene effects at the original row granularity and with the original statistical fields.
    effect_table = direction_table.pivot(index="gene", columns="dataset", values="logFC").reindex(
        index=order, columns=datasets
    )
    max_abs = float(np.nanmax(np.abs(effect_table.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    fig = new_figure(86, 80)
    ax = fig.add_axes([0.24, 0.22, 0.55, 0.56])

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    image = ax.imshow(
        effect_table.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest"
    )
    gsea_table = statistics_table.set_index("dataset")
    ax.set_xticks(
        range(len(datasets)),
        [
            f"GSE174554\n{int(gsea_table.loc['GSE174554', 'n_pairs'])} pairs",
            f"GSE274546\n{int(gsea_table.loc['GSE274546', 'n_pairs'])} pairs",
        ],
        fontsize=6.8,
    )
    ax.set_yticks(range(len(order)), order, fontsize=7.1)
    ax.tick_params(length=0)

    # Process cohorts, genes or panels in the published order.
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(datasets):
            value = float(effect_table.loc[gene, dataset])
            text_color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(
                col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=text_color
            )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(datasets), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.tick_params(which="minor", bottom=False, left=False)
    color_ax = fig.add_axes([0.82, 0.3, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)

    colorbar.outline.set_linewidth(0.45)
    fig.text(0.24, 0.95, "Shared leading-edge genes", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.24,
        0.9,
        "Intersection of the two formal raw20 GSEA leading edges",
        fontsize=6.2,
        color="#717171",
        va="top",
    )
    fig.text(
        0.24,
        0.145,
        (
            "Cohort GSEA: NES "
            f"{gsea_table.loc['GSE174554', 'NES']:.2f}"
            " / "
            f"{gsea_table.loc['GSE274546', 'NES']:.2f}"
            "\nNominal P "
            f"{format_p(float(gsea_table.loc['GSE174554', 'pval']))}"
            " / "
            f"{format_p(float(gsea_table.loc['GSE274546', 'pval']))}"
        ),
        fontsize=6.1,
        color="#555555",
        va="top",
        linespacing=1.25,
    )
    fig.text(
        0.24, 0.065, "0/14 individual gene tests reached FDR < 0.05", fontsize=5.9, color="#777777", va="top"
    )

    # Save the PDF master at the published size and render its PNG from that PDF.
    pdf, png = save_panel(fig, fig_out, stem, "Step45 recurrent GBM Figure 1 program rebuild")
    panel_record_5 = panel_record(
        "Figure1E",
        stem,
        "Seven shared leading-edge genes across the two recurrence cohorts",
        source,
        pdf,
        png,
    )
    records = [panel_record_1, panel_record_2, panel_record_3, panel_record_4, panel_record_5]
    panel_paths = {record["panel"]: Path(record["png"]) for record in records}

    # Assemble the review preview on a Pillow canvas using the original panel rectangles.
    layout_rows = [
        {"panel": "Figure1A", "left": 270, "top": 130, "right": 3330, "bottom": 1700},
        {"panel": "Figure1B", "left": 120, "top": 1810, "right": 1740, "bottom": 3150},
        {"panel": "Figure1C", "left": 1860, "top": 1810, "right": 3480, "bottom": 3150},
        {"panel": "Figure1D", "left": 120, "top": 3260, "right": 1740, "bottom": 4830},
        {"panel": "Figure1E", "left": 1860, "top": 3260, "right": 3480, "bottom": 4830},
    ]
    canvas = Image.new("RGB", (3600, 4960), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = image_font(52, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        paste_contained(
            canvas,
            panel_paths[panel],
            (int(row["left"]), int(row["top"]), int(row["right"]), int(row["bottom"])),
        )
        drawer.text((int(row["left"]), int(row["top"]) - 56), panel[-1], fill="#111111", font=label_font)

    # Save this step's outputs to the published paths.
    png = fig_out / "Figure1_candidate_contact_sheet.png"
    pdf = fig_out / "Figure1_candidate_contact_sheet.pdf"
    canvas.save(png, format="PNG", dpi=(300, 300), optimize=True)
    contact_fig = plt.figure(figsize=(12, 4960 / 300), facecolor="white")
    contact_ax = contact_fig.add_axes([0, 0, 1, 1])

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    contact_ax.imshow(np.asarray(canvas), interpolation="nearest")
    contact_ax.axis("off")
    contact_fig.canvas.draw()
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

    # Check that both formats of the original-size Pillow preview exist and are nonempty.
    for preview_path in [pdf, png]:
        if not preview_path.is_file() or preview_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty preview: {preview_path.name}")
    layout_table = pd.DataFrame(layout_rows)

    # Use the published output directory and basename.
    source = source_out / "Figure1_candidate_contact_sheet_layout.csv"
    csv_path = source
    layout_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    (contact_pdf, contact_png, contact_source) = (pdf, png, source)

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
    csv_path = write_root / "Figure1_panel_manifest.csv"
    manifest_table.to_csv(csv_path, index=False)

    # Use the existing input checks and branch conditions.
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
