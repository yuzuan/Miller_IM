#!/usr/bin/env python3
# Purpose: Build Figure 1E and 1F candidates and review previews from Step 45 ranks and existing LOPO results.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05f_shared_genes_and_robustness.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
from pathlib import Path
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy.stats import spearmanr
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
image_font = figure_style.image_font
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel
panel_record = figure_style.panel_record

if __name__ == "__main__":

    # Keep the seed at this point to reproduce the sampling or jitter sequence.
    np.random.seed(20260713)
    apply_publication_style()
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write/46_figure1_ef_candidates/Figure1"
    source_out = write_root / "source_data"

    # Use the published output directory and basename.
    source_out.mkdir(parents=True, exist_ok=True)
    fig_out = root / "figures/46_figure1_ef_candidates/Figure1"
    fig_out.mkdir(parents=True, exist_ok=True)
    stem = "Figure1E_cross_cohort_gene_rank_comparison"
    source = source_out / f"{stem}_source.csv"

    # Use the published output directory and basename.
    step45_source = root / "write/45_figure1_program_rebuild/Figure1/source_data"
    input_path = step45_source / "Figure1D_cross_cohort_all_gene_rank_density_source.csv"
    rank_table = pd.read_csv(input_path)
    required = {
        "gene",
        "rank_stat_GSE174554",
        "rank_stat_GSE274546",
        "miller_raw20",
        "density",
        "shared_leading_edge",
    }
    if not required.issubset(rank_table.columns):
        raise ValueError(f"Missing Figure1E source columns: {sorted(required - set(rank_table.columns))}")

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = source
    rank_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    fig = new_figure(89.5, 82)
    ax = fig.add_axes([0.17, 0.16, 0.73, 0.7])

    # Separate background genes from the highlighted Miller-IM genes and keep the original density order.
    background_table = rank_table.loc[~rank_table["miller_raw20"]].sort_values("density")
    raw20_other_table = rank_table.loc[rank_table["miller_raw20"] & ~rank_table["shared_leading_edge"]]
    shared_table = rank_table.loc[rank_table["shared_leading_edge"]]
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            (
                "Expected 7 shared and 12 other measured raw20 genes, found "
                f"{len(shared_table)}"
                " and "
                f"{len(raw20_other_table)}"
            )
        )
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}

    # Use the existing input checks and branch conditions.
    if set(shared_table["gene"]) != expected_shared_leading_edge:
        raise ValueError(
            (
                "Unexpected shared leading-edge genes: expected "
                f"{sorted(expected_shared_leading_edge)}"
                ", found "
                f"{sorted(shared_table['gene'])}"
            )
        )
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c=background_table["density"],
        cmap="Greys",
        s=2.0,
        alpha=0.68,
        linewidths=0,
        rasterized=True,
    )
    highlight = "#E64B35"
    ax.scatter(
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s=22,
        facecolor="white",
        edgecolor=highlight,
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s=28,
        color=highlight,
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )

    # Apply the original labels and placement rules.
    label_offsets = {"CCL4": (-25, 7), "CH25H": (5, 5), "SGK1": (-23, -10), "FOLR2": (5, 7), "PDK4": (5, -10)}
    for _, row in shared_table.loc[shared_table["gene"].isin(label_offsets)].iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=label_offsets[row["gene"]],
            textcoords="offset points",
            fontsize=5.7,
            color="#B53A2B",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = (
        float(np.nanmax(np.abs(rank_table[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    clean_axis(ax)
    rho, _ = spearmanr(rank_table["rank_stat_GSE174554"], rank_table["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"Shared genes = {len(rank_table):,}\nRank correlation $r_s$ = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.4,
        color="#686868",
    )
    ax.text(
        0.98,
        0.03,
        "7 shared leading-edge genes\n12 other measured raw20 genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.0,
        color="#B53A2B",
        linespacing=1.15,
    )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.text(
        0.98,
        0.98,
        "Positive = recurrent-enriched",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color="#777777",
    )
    fig.text(0.17, 0.96, "Cross-cohort gene-rank comparison", fontsize=8.6, fontweight="bold", va="top")
    pdf, png = save_panel(fig, fig_out, stem, "Step46 recurrent GBM Figure 1 E-F candidates")
    panel_record_1 = panel_record(
        "Figure1E",
        stem,
        "All-gene cross-cohort rank-stat comparison with shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )

    # Add the downstream files to the manifest in the published panel order.
    records = [panel_record_1]
    source = source_out / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    if not source.exists():
        raise FileNotFoundError(f"Run 05e_leave_one_patient_out.R first: {source}")
    fold_table = pd.read_csv(source)

    # Keep the prespecified 20 genes in their original order.
    if set(fold_table["pathway"]) != {"Miller_Microglial_Inflammatory_raw_top20"}:
        raise ValueError("Figure1F source must contain only the fixed Miller raw20 program")
    full_table = fold_table.loc[fold_table["run_type"].eq("Full")].set_index("dataset")
    lopo_table = fold_table.loc[fold_table["run_type"].eq("Leave-one-patient-out")].copy()
    expected_counts = {"GSE174554": 17, "GSE274546": 45}

    # Use the existing input checks and branch conditions.
    if lopo_table.groupby("dataset").size().to_dict() != expected_counts:
        raise ValueError(f"Unexpected LOPO counts: {lopo_table.groupby('dataset').size().to_dict()}")
    summary_rows = []
    datasets = ["GSE174554", "GSE274546"]
    for dataset in datasets:
        cohort_table = lopo_table.loc[lopo_table["dataset"].eq(dataset)]
        summary_rows.append(
            {
                "dataset": dataset,
                "n_omissions": len(cohort_table),
                "full_NES": float(full_table.loc[dataset, "NES"]),
                "full_nominal_P": float(full_table.loc[dataset, "nominal_P"]),
                "LOPO_NES_min": float(cohort_table["NES"].min()),
                "LOPO_NES_max": float(cohort_table["NES"].max()),
                "positive_NES_folds": int((cohort_table["NES"] > 0).sum()),
                "nominal_P_lt_0_05_folds": int((cohort_table["nominal_P"] < 0.05).sum()),
                "max_nominal_P": float(cohort_table["nominal_P"].max()),
            }
        )
    summary_table = pd.DataFrame(summary_rows)

    # Use the published output directory and basename.
    summary_path = source_out / "Figure1F_leave_one_patient_out_summary.csv"
    csv_path = summary_path
    summary_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    fig = new_figure(89.5, 82)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax = fig.add_axes([0.24, 0.21, 0.71, 0.65])
    y_positions = {"GSE174554": 1.0, "GSE274546": 0.0}
    rng = np.random.default_rng(20260713)
    all_values = []
    dataset_colors = {"GSE174554": "#3C5488", "GSE274546": "#00A087"}

    # Process cohorts, genes or panels in the published order.
    for dataset in datasets:

        # Select the current cohort's records in their original order and keep the statistical fields.
        cohort_table = (
            lopo_table.loc[lopo_table["dataset"].eq(dataset)].sort_values(["NES", "omitted_patient"]).copy()
        )
        y = y_positions[dataset]
        color = dataset_colors[dataset]
        all_values.extend(cohort_table["NES"].tolist())
        ax.hlines(
            y,
            cohort_table["NES"].min(),
            cohort_table["NES"].max(),
            color=color,
            linewidth=2.3,
            alpha=0.42,
            zorder=1,
        )

        jitter_grid = np.linspace(-0.13, 0.13, len(cohort_table))
        jitter = jitter_grid[rng.permutation(len(cohort_table))]
        significant = cohort_table["nominal_P"].to_numpy(dtype=float) < 0.05
        ax.scatter(
            cohort_table.loc[significant, "NES"],
            y + jitter[significant],
            s=17,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=0.92,
            zorder=3,
        )
        ax.scatter(
            cohort_table.loc[~significant, "NES"],
            y + jitter[~significant],
            s=18,
            facecolor="white",
            edgecolor=color,
            linewidth=0.85,
            alpha=1.0,
            zorder=3,
        )

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.scatter(
            [float(full_table.loc[dataset, "NES"])],
            [y],
            marker="D",
            s=48,
            color="#222222",
            edgecolor="white",
            linewidth=0.65,
            zorder=5,
        )
        row = summary_table.loc[summary_table["dataset"].eq(dataset)].iloc[0]
        ax.text(
            0.985,
            y + 0.23,
            (
                "NES "
                f"{row['LOPO_NES_min']:.2f}"
                "–"
                f"{row['LOPO_NES_max']:.2f}"
                "   positive "
                f"{int(row['positive_NES_folds'])}"
                "/"
                f"{int(row['n_omissions'])}"
                "\nP<0.05 "
                f"{int(row['nominal_P_lt_0_05_folds'])}"
                "/"
                f"{int(row['n_omissions'])}"
                "   max P "
                f"{row['max_nominal_P']:.3g}"
            ),
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=6.0,
            color="#555555",
            linespacing=1.15,
        )
    x_max = max(all_values + full_table["NES"].tolist())
    ax.set_xlim(-0.05, x_max + 0.18)
    ax.set_ylim(-0.48, 1.48)
    ax.axvline(0, color="#999999", linestyle=(0, (3, 2)), linewidth=0.7)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_yticks(
        [y_positions[dataset] for dataset in datasets],
        [f"{dataset}\n{expected_counts[dataset]} omissions" for dataset in datasets],
        fontsize=6.9,
    )
    for label, dataset in zip(ax.get_yticklabels(), datasets):
        label.set_color(dataset_colors[dataset])
        label.set_fontweight("bold")
    ax.set_xlabel("Leave-one-patient-out normalized enrichment score")
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.45)
    clean_axis(ax, keep_left=False, keep_bottom=True)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.tick_params(axis="y", length=0, pad=5)
    fig.text(0.16, 0.96, "Leave-one-patient-out GSEA stability", fontsize=8.6, fontweight="bold", va="top")
    fig.text(
        0.16,
        0.075,
        "Each circle: one omission   Filled: nominal P < 0.05   Diamond: full data",
        fontsize=5.8,
        color="#666666",
        va="bottom",
    )
    stem = "Figure1F_leave_one_patient_out_gsea_stability"
    pdf, png = save_panel(fig, fig_out, stem, "Step46 recurrent GBM Figure 1 E-F candidates")
    (panel_f, summary_path) = (
        panel_record(
            "Figure1F", stem, "Leave-one-patient-out paired pseudobulk raw20 GSEA stability", source, pdf, png
        ),
        summary_path,
    )

    records.append(panel_f)
    paths = {record["panel"]: Path(record["png"]) for record in records}
    layout_rows = [
        {"panel": "Figure1E", "left": 120, "top": 170, "right": 1750, "bottom": 1770},
        {"panel": "Figure1F", "left": 1850, "top": 170, "right": 3480, "bottom": 1770},
    ]
    canvas = Image.new("RGB", (3600, 1900), "white")

    # Add panel letters to the preview canvas.
    drawer = ImageDraw.Draw(canvas)
    label_font = image_font(54, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        paste_contained(canvas, paths[panel], (row["left"], row["top"], row["right"], row["bottom"]))
        drawer.text((int(row["left"]), 70), panel[-1], fill="#111111", font=label_font)
    png = fig_out / "Figure1EF_candidate_preview.png"
    pdf = fig_out / "Figure1EF_candidate_preview.pdf"

    # Use the published output directory and basename.
    source = source_out / "Figure1EF_candidate_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(300, 300), optimize=True)
    preview_fig = plt.figure(figsize=(12, 1900 / 300), facecolor="white")
    preview_ax = preview_fig.add_axes([0, 0, 1, 1])
    preview_ax.imshow(np.asarray(canvas), interpolation="nearest")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    preview_ax.axis("off")
    preview_fig.canvas.draw()
    preview_fig.savefig(
        pdf,
        dpi=300,
        facecolor="white",
        edgecolor="none",
        metadata={
            "Creator": "Step46 recurrent GBM Figure 1 E-F candidates",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(preview_fig)

    # Check that both formats of the original-size Pillow preview exist and are nonempty.
    for preview_path in [pdf, png]:
        if not preview_path.is_file() or preview_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty preview: {preview_path.name}")
    output_table = pd.DataFrame(layout_rows)
    csv_path = source

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    output_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    (preview_pdf, preview_png, preview_source) = (pdf, png, source)
    records.append(
        panel_record(
            "Figure1EF_preview",
            "Figure1EF_candidate_preview",
            "Review-only preview of Figure 1E-F candidates",
            preview_source,
            preview_pdf,
            preview_png,
        )
    )

    # Build the table from the named records in their original order and with the same unit per row.
    manifest_table = pd.DataFrame(records)
    manifest_table["lopo_summary"] = ""
    manifest_table.loc[manifest_table["panel"].eq("Figure1F"), "lopo_summary"] = str(summary_path)
    csv_path = write_root / "Figure1EF_panel_manifest.csv"
    manifest_table.to_csv(csv_path, index=False)

    # Use the existing input checks and branch conditions.
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
