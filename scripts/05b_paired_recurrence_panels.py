#!/usr/bin/env python3
# Purpose: Replot cohort design, UMAP, marker dot plots and GSEA from Steps 38 and 39.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05b_paired_recurrence_panels.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
import shutil
from pathlib import Path
import matplotlib as mpl
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
colors = figure_style.colors
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
save_figure = figure_style.save_figure
enrichment_curve = figure_style.enrichment_curve

if __name__ == "__main__":

    # Set the original seed and style, then load the source tables for the panels.

    np.random.seed(20260713)
    apply_publication_style()
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write/41_mg_inflammatory_sci_rebuild/Figure1"

    # Use the published output directory and basename.
    source_out = write_root / "source_data"
    source_out.mkdir(parents=True, exist_ok=True)
    fig_out = root / "figures/41_mg_inflammatory_sci_rebuild/Figure1/panel_library"
    fig_out.mkdir(parents=True, exist_ok=True)
    required = [
        "identity_marker_dotplot.csv",
        "umap_identity_GSE174554.csv.gz",
        "umap_identity_GSE274546.csv.gz",
    ]

    # Use the published output directory and basename.
    old_source = root / "write/39_independent_miller_mg_inflammatory_figures/Figure1/source_data"
    for name in required:
        src = old_source / name
        if not src.exists():
            raise FileNotFoundError(src)
        shutil.copy2(src, source_out / name)
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    audit_table = pd.read_csv(step38 / "independent_input_pair_audit.csv")
    audit_required = {"dataset", "threshold", "n_pairs", "n_cells_primary", "n_cells_recurrent"}

    # Use the existing input checks and branch conditions.
    if not audit_required.issubset(audit_table.columns):
        raise ValueError(
            f"Missing Step38 pair-audit columns: {sorted(audit_required - set(audit_table.columns))}"
        )
    datasets = ["GSE174554", "GSE274546"]
    workflow_table = audit_table.loc[
        audit_table["dataset"].isin(datasets) & audit_table["threshold"].eq(20)
    ].copy()
    if "formal_testable" in workflow_table.columns:
        is_formal = workflow_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
        workflow_table = workflow_table.loc[is_formal].copy()
    if workflow_table.groupby("dataset").size().to_dict() != {dataset: 1 for dataset in datasets}:
        raise ValueError(
            (
                "Step38 pair audit must contain exactly one formal threshold-20 ro"
                "w per cohort: "
                f"{workflow_table.groupby('dataset').size().to_dict()}"
            )
        )

    # Build the design panel from the audit rows, using the existing cell and patient-pair counts.
    workflow_table["input_libraries_or_matrices"] = workflow_table["dataset"].map(
        {"GSE174554": 91, "GSE274546": 111}
    )
    workflow_table["clean_myeloid"] = (
        pd.to_numeric(workflow_table["n_cells_primary"], errors="raise")
        + pd.to_numeric(workflow_table["n_cells_recurrent"], errors="raise")
    ).astype(int)
    workflow_table["formal_pairs_threshold20"] = pd.to_numeric(
        workflow_table["n_pairs"], errors="raise"
    ).astype(int)
    workflow_table["cell_count_scope"] = "IDH-wildtype cells included in the formal threshold-20 analysis"
    output_table = workflow_table[
        [
            "dataset",
            "input_libraries_or_matrices",
            "clean_myeloid",
            "formal_pairs_threshold20",
            "n_cells_primary",
            "n_cells_recurrent",
            "cell_count_scope",
        ]
    ].sort_values("dataset")

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = source_out / "cohort_workflow.csv"
    output_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    gsea_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
    gsea_is_formal = gsea_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})

    # Keep the prespecified 20 genes in their original order.
    raw20_name = "Miller_Microglial_Inflammatory_raw_top20"
    gsea_table = gsea_table[
        gsea_table["dataset"].isin(datasets)
        & gsea_table["threshold"].eq(20)
        & gsea_table["pathway"].eq(raw20_name)
        & gsea_is_formal
    ].copy()
    if len(gsea_table) != 2:
        raise ValueError(f"Expected two primary GSEA rows, found {len(gsea_table)}")
    pair_counts = gsea_table.set_index("dataset")["n_pairs"].astype(int).to_dict()
    if pair_counts != {"GSE174554": 17, "GSE274546": 45}:
        raise ValueError(f"Unexpected formal pair counts in Step38 GSEA: {pair_counts}")

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = source_out / "primary_gsea_statistics.csv"
    gsea_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    provenance_table = pd.read_csv(step38 / "fixed_program_provenance.csv")
    raw20_rows = provenance_table.loc[provenance_table["signature"].eq(raw20_name), "genes"]

    # Use the existing input checks and branch conditions.
    if len(raw20_rows) != 1:
        raise ValueError(f"Expected one {raw20_name} provenance row, found {len(raw20_rows)}")
    raw20_genes = str(raw20_rows.iloc[0]).split(";")
    deg_table = pd.read_csv(step38 / "independent_paired_edger_all_genes.csv")
    curves = []


    # Process cohorts, genes or panels in the published order.
    for dataset in datasets:
        current_table = deg_table.loc[deg_table["dataset"].eq(dataset) & deg_table["threshold"].eq(20)].copy()
        if "formal_testable" in current_table.columns:
            is_formal = current_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
            current_table = current_table.loc[is_formal].copy()
        if current_table.empty:
            raise ValueError(f"No formal Step38 differential-expression rows for {dataset}")
        curve_table = enrichment_curve(current_table, raw20_genes)
        curve_table["dataset"] = dataset

        curves.append(curve_table)
    output_table_total = pd.concat(curves, ignore_index=True)
    csv_path = source_out / "raw20_gsea_curves.csv.gz"
    output_table_total.to_csv(csv_path, index=False, compression={"method": "gzip", "mtime": 0})
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Read the frozen inputs with their original row order, dataset names and column names.
    workflow_table = pd.read_csv(source_out / "cohort_workflow.csv")
    workflow_table = workflow_table.set_index("dataset")
    fig = new_figure(183, 42)
    ax = fig.add_axes([0.015, 0.08, 0.97, 0.84])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.axis("off")
    ax.text(0.0, 0.94, "Independent paired-cohort design", fontsize=9, weight="bold", va="top")
    ax.text(
        1.0,
        0.94,
        "Primary vs first recurrence; patient-level inference",
        fontsize=7.2,
        color=colors["neutral"],
        va="top",
        ha="right",
    )
    x_positions = {"GSE174554": 0.0, "GSE274546": 0.38}
    widths = 0.34

    subtitles = {"GSE174554": "snRNA-seq", "GSE274546": "snRNA-seq"}
    dataset_colors = {"GSE174554": colors["gse174554"], "GSE274546": colors["gse274546"]}
    for dataset in datasets:
        x = x_positions[dataset]
        color = dataset_colors[dataset]
        row = workflow_table.loc[dataset]
        ax.add_patch(mpl.patches.Rectangle((x, 0.2), widths, 0.56, facecolor="#F7F7F7", edgecolor="none"))
        ax.add_patch(mpl.patches.Rectangle((x, 0.2), 0.012, 0.56, facecolor=color, edgecolor="none"))

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(x + 0.03, 0.69, dataset, fontsize=8.4, weight="bold", color=color, va="top")
        ax.text(
            x + widths - 0.02,
            0.69,
            subtitles[dataset],
            fontsize=7,
            color=colors["neutral"],
            ha="right",
            va="top",
        )
        metrics = [
            ("Input", f"{int(row.input_libraries_or_matrices):,} libraries/matrices"),
            ("Analyzed", f"{int(row.clean_myeloid):,} myeloid cells"),
            ("Paired", f"{int(row.formal_pairs_threshold20)} patients"),
        ]
        for i, (label, value) in enumerate(metrics):
            y = 0.55 - i * 0.14
            ax.text(x + 0.03, y, label.upper(), fontsize=6.2, color=colors["neutral"], va="center")
            ax.text(x + 0.12, y, value, fontsize=7.4, weight="bold", va="center")
    ax.annotate(
        "",
        xy=(0.785, 0.48),
        xytext=(0.735, 0.48),
        arrowprops=dict(arrowstyle="-|>", lw=0.9, color=colors["neutral"]),
    )
    ax.add_patch(mpl.patches.Rectangle((0.79, 0.2), 0.21, 0.56, facecolor="#F4F4F4", edgecolor="none"))

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.text(0.81, 0.67, "PAIRED ANALYSIS", fontsize=6.2, color=colors["neutral"], va="top")
    ax.text(0.81, 0.54, "Raw-count\npseudobulk", fontsize=8.0, weight="bold", va="top", linespacing=1.05)
    ax.text(0.81, 0.34, "edgeR ranking +\nfixed-program GSEA", fontsize=6.9, va="top", linespacing=1.08)
    save_figure(fig, fig_out, "fig1a_independent_paired_design")
    identity_colors = {
        "Microglia": colors["microglia"],
        "Macrophage": colors["macrophage"],
        "Monocyte": colors["monocyte"],
        "CDC": colors["cdc"],
        "Neutrophil": colors["neutrophil"],
    }

    display_identity = {"CDC": "cDC"}

    def label_position(points: np.ndarray) -> np.ndarray:
        """Return the central label coordinates using the original quantile cutoff."""

        center = np.nanmedian(points, axis=0)
        distances = np.sqrt(((points - center) ** 2).sum(axis=1))
        keep = points[distances <= np.nanquantile(distances, 0.65)]
        result = np.nanmedian(keep, axis=0) if len(keep) else center
        return result

    def panel_identity_umap(dataset: str, stem: str) -> None:
        """Read the cohort identity coordinates and save the UMAP PDF and PNG. No return value."""

        # Read the frozen inputs with their original row order, dataset names and column names.
        umap_table = pd.read_csv(source_out / f"umap_identity_{dataset}.csv.gz")
        fig = new_figure(88, 78)
        ax = fig.add_axes([0.04, 0.04, 0.92, 0.76])
        counts = umap_table["myeloid_identity"].value_counts()
        for identity in counts.sort_values(ascending=False).index:
            current_table = umap_table[umap_table["myeloid_identity"].eq(identity)]
            ax.scatter(
                current_table["UMAP1"],
                current_table["UMAP2"],
                s=0.42,
                color=identity_colors.get(identity, "#BDBDBD"),
                alpha=0.56,
                linewidths=0,
                rasterized=True,
            )

        # Process cohorts, genes or panels in the published order.
        for identity in counts.index:

            # Apply the original labels and placement rules.
            label_points = umap_table[umap_table["myeloid_identity"].eq(identity)][
                ["UMAP1", "UMAP2"]
            ].to_numpy()
            (x, y) = label_position(label_points)
            label = display_identity.get(identity, identity)
            ax.text(
                x,
                y,
                label,
                fontsize=7.2,
                weight="bold",
                ha="center",
                va="center",
                color=identity_colors.get(identity, colors["text"]),
                path_effects=[pe.withStroke(linewidth=2.6, foreground="white")],
            )
        fig.text(0.06, 0.95, dataset, ha="left", va="top", fontsize=8.5, weight="bold")
        fig.text(
            0.06,
            0.875,
            f"{len(umap_table):,} independently annotated myeloid cells",
            ha="left",
            va="top",
            fontsize=6.8,
            color=colors["neutral"],
        )
        ax.set_xticks([])
        ax.set_yticks([])

        # Process cohorts, genes or panels in the published order.
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.margins(0.015)
        ax.set_aspect("equal", adjustable="box")
        save_figure(fig, fig_out, stem)

    panel_identity_umap("GSE174554", "fig1b_gse174554_identity_umap")
    panel_identity_umap("GSE274546", "fig1c_gse274546_identity_umap")
    marker_table = pd.read_csv(source_out / "identity_marker_dotplot.csv")

    # Summarize expression by cohort, cell identity and gene using the original markers and scaling.
    marker_blocks = {
        "Microglia": ["P2RY12", "TMEM119", "SALL1"],
        "Macrophage": ["GPNMB", "SPP1", "MSR1"],
        "Monocyte": ["VCAN", "FCN1", "S100A8"],
        "cDC": ["FCER1A", "CD1C", "CLEC10A"],
    }
    genes = [gene for block in marker_blocks.values() for gene in block]
    identities = ["Microglia", "Macrophage", "Monocyte", "CDC"]
    marker_table = marker_table[
        marker_table["gene"].isin(genes) & marker_table["identity"].isin(identities)
    ].copy()
    x_map = {gene: i for (i, gene) in enumerate(genes)}

    y_map = {identity: i for (i, identity) in enumerate(identities)}
    fig = new_figure(183, 77)
    axes = [fig.add_axes([0.1, 0.2, 0.33, 0.56]), fig.add_axes([0.51, 0.2, 0.33, 0.56])]
    norm = mpl.colors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    cmap = mpl.colormaps["RdBu_r"]

    # Process cohorts, genes or panels in the published order.
    for ax, dataset in zip(axes, datasets):
        current_table = marker_table[marker_table["dataset"].eq(dataset)]
        x = current_table["gene"].map(x_map).to_numpy()
        y = current_table["identity"].map(y_map).to_numpy()
        size = 5 + np.clip(current_table["pct_positive"].to_numpy(), 0, 100) * 0.44
        ax.scatter(
            x,
            y,
            s=size,
            c=current_table["scaled_expression"],
            cmap=cmap,
            norm=norm,
            edgecolor="#555555",
            linewidth=0.22,
        )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.set_xlim(-0.55, len(genes) - 0.45)
        ax.set_ylim(len(identities) - 0.45, -0.55)
        ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(identities)), [display_identity.get(x, x) for x in identities])
        ax.set_title(dataset, color=dataset_colors[dataset], pad=7)

        # Process cohorts, genes or panels in the published order.
        for boundary in [2.5, 5.5, 8.5]:
            ax.axvline(boundary, color="#CFCFCF", lw=0.7)
        ax.set_axisbelow(True)
        ax.grid(color=colors["grid"], lw=0.45)
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

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    lax.set_ylim(0, 1)
    lax.text(0, 0.98, "Cells expressing", fontsize=6.8, va="top")
    for i, pct in enumerate([25, 50, 75]):
        y = 0.68 - i * 0.27
        lax.scatter(0.13, y, s=5 + pct * 0.44, facecolor="#BDBDBD", edgecolor="#555555", linewidth=0.25)
        lax.text(0.34, y, f"{pct}%", fontsize=6.4, va="center")
    fig.text(
        0.1,
        0.95,
        "Canonical markers support independent myeloid annotation",
        fontsize=8.5,
        weight="bold",
        va="top",
    )
    save_figure(fig, fig_out, "fig1d_cross_cohort_identity_markers")

    # Read the frozen inputs with their original row order, dataset names and column names.
    def panel_gsea(dataset: str, stem: str) -> None:
        """Save the three-track GSEA panel using the cohort, output stem and frozen statistics. No return value."""

        # Read the frozen inputs with their original row order, dataset names and column names.
        curve_table = pd.read_csv(source_out / "raw20_gsea_curves.csv.gz")
        curve_table = curve_table[curve_table["dataset"].eq(dataset)].sort_values("rank")
        statistics_table = pd.read_csv(source_out / "primary_gsea_statistics.csv")
        row = statistics_table[statistics_table["dataset"].eq(dataset)].iloc[0]
        nominal_p = float(row["pval"])

        color = dataset_colors[dataset]
        fig = new_figure(88, 76)
        gs = fig.add_gridspec(
            3, 1, left=0.16, right=0.96, bottom=0.15, top=0.79, height_ratios=[3.5, 0.48, 0.95], hspace=0.08
        )
        ax = fig.add_subplot(gs[0])
        hit_ax = fig.add_subplot(gs[1], sharex=ax)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        rank_ax = fig.add_subplot(gs[2], sharex=ax)
        rank = curve_table["rank"].to_numpy()
        running = curve_table["running_enrichment"].to_numpy()
        ax.plot(rank, running, color=color, lw=1.45)
        ax.fill_between(rank, 0, running, color=color, alpha=0.12, linewidth=0)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.axhline(0, color="#999999", lw=0.55)
        ax.set_ylabel("Running enrichment")
        ax.set_xticks([])
        clean_axis(ax, keep_left=True, keep_bottom=False)
        ax.text(
            0.0, 1.035, dataset, transform=ax.transAxes, fontsize=8.5, weight="bold", ha="left", va="bottom"
        )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(
            1.0,
            1.035,
            f"NES {row.NES:.2f}   nominal P {nominal_p:.2g}   n={int(row.n_pairs)} pairs",
            transform=ax.transAxes,
            fontsize=7,
            ha="right",
            va="bottom",
        )
        hits = curve_table.loc[curve_table["hit"].astype(bool), "rank"].to_numpy()
        hit_ax.vlines(hits, 0, 1, color="#333333", lw=0.48)
        hit_ax.set_ylim(0, 1)
        hit_ax.axis("off")

        # Prepare the cohort gene ranks using the original order and merge keys.
        rank_stat = curve_table["rank_stat"].to_numpy()
        rank_ax.fill_between(
            rank, 0, np.clip(rank_stat, 0, None), color=colors["recurrent"], alpha=0.75, linewidth=0
        )
        rank_ax.fill_between(
            rank, 0, np.clip(rank_stat, None, 0), color=colors["primary"], alpha=0.75, linewidth=0
        )
        rank_ax.axhline(0, color="#777777", lw=0.45)
        rank_ax.set_yticks([])

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        rank_ax.set_xlabel("Ranked genes")
        rank_ax.text(
            0.0,
            -0.58,
            "Recurrent-enriched",
            transform=rank_ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.8,
            color=colors["recurrent"],
        )
        rank_ax.text(
            1.0,
            -0.58,
            "Primary-enriched",
            transform=rank_ax.transAxes,
            ha="right",
            va="top",
            fontsize=6.8,
            color=colors["primary"],
        )
        clean_axis(rank_ax, keep_left=False, keep_bottom=True)
        rank_ax.spines["bottom"].set_color("#777777")

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        fig.text(
            0.16,
            0.94,
            "Miller microglial-inflammatory program",
            fontsize=7.2,
            color=colors["neutral"],
            va="top",
        )
        save_figure(fig, fig_out, stem)

    panel_gsea("GSE174554", "fig1e_gse174554_raw20_gsea")
    panel_gsea("GSE274546", "fig1f_gse274546_raw20_gsea")
    rows = [
        ("Figure1A", "fig1a_independent_paired_design", "Study design and patient-level inference"),
        ("Figure1B", "fig1b_gse174554_identity_umap", "GSE174554 independent myeloid identities"),
        ("Figure1C", "fig1c_gse274546_identity_umap", "GSE274546 independent myeloid identities"),
        ("Figure1D", "fig1d_cross_cohort_identity_markers", "Canonical identity marker expression"),
        ("Figure1E", "fig1e_gse174554_raw20_gsea", "GSE174554 paired pseudobulk GSEA"),
        ("Figure1F", "fig1f_gse274546_raw20_gsea", "GSE274546 paired pseudobulk GSEA"),
    ]
    manifest_table = pd.DataFrame(rows, columns=["panel", "stem", "message"])

    manifest_table["png"] = manifest_table["stem"].map(lambda x: str(fig_out / f"{x}.png"))
    manifest_table["pdf"] = manifest_table["stem"].map(lambda x: str(fig_out / f"{x}.pdf"))
    csv_path = write_root / "Figure1_panel_manifest.csv"
    manifest_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
