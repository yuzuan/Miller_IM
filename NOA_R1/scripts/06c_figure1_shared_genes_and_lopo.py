#!/usr/bin/env python3
# Purpose: draw the Figure 1E rank comparison and Figure 1F leave-one-patient-out GSEA stability panels.
# Inputs: legacy Step 45 Figure1D_cross_cohort_all_gene_rank_density_source.csv (paired rank statistics,
# MG-Inflam20 membership and shared leading edges) and Step 46 Figure1F_leave_one_patient_out_raw20_gsea.csv
# (cohort, run type, omitted patient, NES and nominal P). Use an existing local fold table if already prepared.
# Outputs: figures/46_figure1_ef_candidates/Figure1/Figure1E_cross_cohort_gene_rank_comparison,
# Figure1F_leave_one_patient_out_gsea_stability and Figure1EF_candidate_preview (.pdf and .png);
# write/46_figure1_ef_candidates/Figure1/Figure1EF_panel_manifest.csv and source_data/ with the
# Figure 1E source table, Figure1F_leave_one_patient_out_summary.csv and Figure1EF_candidate_preview_layout.csv.
# Run: python3 NOA_R1/scripts/06c_figure1_shared_genes_and_lopo.py
# Keep Matplotlib and Pillow to retain the original placement of annotations and preview panels.

from __future__ import annotations

import importlib
import shutil
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw
from scipy.stats import spearmanr

figure_style = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
legacy_project_root = figure_style.legacy_root(revision_root)

# Output globals and new_figure can be redirected by the final figure export.
WRITE_ROOT = revision_root / "write/46_figure1_ef_candidates/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "figures/46_figure1_ef_candidates/Figure1"
new_figure = figure_style.new_figure


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, str]:
    """Return the original manifest fields for the source table and two image formats."""
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


def panel_e_rank_rank() -> dict[str, str]:
    """Plot the cross-cohort ranks and return the Figure 1E manifest row."""
    stem = "Figure1E_cross_cohort_gene_rank_comparison"
    source = SOURCE_OUT / f"{stem}_source.csv"
    step45_source = legacy_project_root / "write/45_figure1_program_rebuild/Figure1/source_data"
    input_path = step45_source / "Figure1D_cross_cohort_all_gene_rank_density_source.csv"
    source_table = pd.read_csv(input_path)
    required = {
        "gene",
        "rank_stat_GSE174554",
        "rank_stat_GSE274546",
        "miller_raw20",
        "density",
        "shared_leading_edge",
    }
    if not required.issubset(source_table.columns):
        raise ValueError(f"Missing Figure1E source columns: {sorted(required - set(source_table.columns))}")
    source_table.to_csv(source, index=False)

    fig = new_figure(89.5, 82)
    ax = fig.add_axes([0.17, 0.16, 0.73, 0.70])
    background_table = source_table.loc[~source_table["miller_raw20"]].sort_values("density")
    raw20_other_table = source_table.loc[source_table["miller_raw20"] & ~source_table["shared_leading_edge"]]
    shared_table = source_table.loc[source_table["shared_leading_edge"]]

    # Verify the displayed partition of the 19 measurable genes against the original leading edges.
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            "Expected 7 shared and 12 other measurable MG-Inflam20 genes, "
            f"found {len(shared_table)} and {len(raw20_other_table)}"
        )
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}
    if set(shared_table["gene"]) != expected_shared_leading_edge:
        raise ValueError(
            "Unexpected shared leading-edge genes: "
            f"expected {sorted(expected_shared_leading_edge)}, found {sorted(shared_table['gene'])}"
        )

    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c     = background_table["density"],
        cmap  = "Greys",
        s     = 2.0,
        alpha = 0.68,
        linewidths = 0,
        rasterized = True,
    )
    highlight = "#E64B35"
    ax.scatter(
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s         = 22,
        facecolor = "white",
        edgecolor = highlight,
        linewidth = 0.85,
        zorder    = 3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s         = 28,
        color     = highlight,
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
    }
    for _, row in shared_table.loc[shared_table["gene"].isin(label_offsets)].iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext       = label_offsets[row["gene"]],
            textcoords   = "offset points",
            fontsize     = 5.7,
            color        = "#B53A2B",
            path_effects = [pe.withStroke(linewidth=1.8, foreground="white")],
        )

    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = float(np.nanmax(np.abs(source_table[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    figure_style.clean_axis(ax)

    rho, _ = spearmanr(source_table["rank_stat_GSE174554"], source_table["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"Shared genes = {len(source_table):,}\nRank correlation $r_s$ = {rho:.2f}",
        transform = ax.transAxes,
        va        = "top",
        fontsize  = 6.4,
        color     = "#686868",
    )
    ax.text(
        0.98,
        0.03,
        f"{len(shared_table)} shared leading-edge genes\n{len(raw20_other_table)} other measured raw20 genes",
        transform   = ax.transAxes,
        ha          = "right",
        va          = "bottom",
        fontsize    = 6.0,
        color       = "#B53A2B",
        linespacing = 1.15,
    )
    ax.text(
        0.98,
        0.98,
        "Positive = recurrent-enriched",
        transform = ax.transAxes,
        ha        = "right",
        va        = "top",
        fontsize  = 5.8,
        color     = "#777777",
    )
    fig.text(0.17, 0.96, "Cross-cohort gene-rank comparison", fontsize=8.6, fontweight="bold", va="top")
    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step46 recurrent GBM Figure 1 E-F candidates")
    manifest_row = panel_record(
        "Figure1E",
        stem,
        "All-gene cross-cohort rank-stat comparison with shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )
    return manifest_row


def panel_f_lopo(*, final_layout: bool = False) -> tuple[dict[str, str], Path]:
    """Summarize the frozen GSEA folds and return the Figure 1F manifest row and summary path."""
    stem = "Figure1F_leave_one_patient_out_gsea_stability"
    source = SOURCE_OUT / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    if not source.exists():
        legacy_source = legacy_project_root / "write/46_figure1_ef_candidates/Figure1/source_data" / source.name
        if not legacy_source.is_file():
            raise FileNotFoundError(f"Prepare the fold table with scripts/05e_leave_one_patient_out.R: {legacy_source}")
        SOURCE_OUT.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_source, source)
    source_table = pd.read_csv(source)
    if set(source_table["pathway"]) != {"Miller_Microglial_Inflammatory_raw_top20"}:
        raise ValueError("Figure1F source must contain only the fixed Miller raw20 program")
    full_table = source_table.loc[source_table["run_type"].eq("Full")].set_index("dataset")
    lopo_table = source_table.loc[source_table["run_type"].eq("Leave-one-patient-out")].copy()

    # The omissions match the frozen 17- and 45-pair cohorts.
    expected_counts = {"GSE174554": 17, "GSE274546": 45}
    if lopo_table.groupby("dataset").size().to_dict() != expected_counts:
        raise ValueError(f"Unexpected LOPO counts: {lopo_table.groupby('dataset').size().to_dict()}")

    # Summarize the original fold tests at nominal P < 0.05 using unrounded values.
    datasets = ["GSE174554", "GSE274546"]
    summary_rows = []
    for dataset in datasets:
        folds_table = lopo_table.loc[lopo_table["dataset"].eq(dataset)]
        summary_rows.append(
            {
                "dataset": dataset,
                "n_omissions": len(folds_table),
                "full_NES": float(full_table.loc[dataset, "NES"]),
                "full_nominal_P": float(full_table.loc[dataset, "nominal_P"]),
                "LOPO_NES_min": float(folds_table["NES"].min()),
                "LOPO_NES_max": float(folds_table["NES"].max()),
                "positive_NES_folds": int((folds_table["NES"] > 0).sum()),
                "nominal_P_lt_0_05_folds": int((folds_table["nominal_P"] < 0.05).sum()),
                "max_nominal_P": float(folds_table["nominal_P"].max()),
            }
        )
    summary_table = pd.DataFrame(summary_rows)
    summary_path = SOURCE_OUT / "Figure1F_leave_one_patient_out_summary.csv"
    summary_table.to_csv(summary_path, index=False)

    # The October final export uses a wide panel and the candidate preview uses a single-column panel.
    if final_layout:
        fig = new_figure(170, 40)
    else:
        fig = new_figure(89.5, 82)
    ax = fig.add_axes([0.24, 0.40, 0.71, 0.42])
    y_positions = {"GSE174554": 1.0, "GSE274546": 0.0}
    dataset_colors = {"GSE174554": "#3C5488", "GSE274546": "#00A087"}

    # The frozen seed preserves the displayed position of each leave-one-out dot.
    rng = np.random.default_rng(20260713)
    all_values = []
    for dataset in datasets:
        folds_table = lopo_table.loc[lopo_table["dataset"].eq(dataset)].sort_values(["NES", "omitted_patient"]).copy()
        y = y_positions[dataset]
        color = dataset_colors[dataset]
        all_values.extend(folds_table["NES"].tolist())
        ax.hlines(
            y, folds_table["NES"].min(), folds_table["NES"].max(), color=color, linewidth=2.3, alpha=0.42, zorder=1
        )
        jitter_grid = np.linspace(-0.13, 0.13, len(folds_table))
        jitter = jitter_grid[rng.permutation(len(folds_table))]
        significant = folds_table["nominal_P"].to_numpy(dtype=float) < 0.05
        ax.scatter(
            folds_table.loc[significant, "NES"],
            y + jitter[significant],
            s         = 17,
            color     = color,
            edgecolor = "white",
            linewidth = 0.35,
            alpha     = 0.92,
            zorder    = 3,
        )
        ax.scatter(
            folds_table.loc[~significant, "NES"],
            y + jitter[~significant],
            s         = 18,
            facecolor = "white",
            edgecolor = color,
            linewidth = 0.85,
            alpha     = 1.0,
            zorder    = 3,
        )
        ax.scatter(
            [float(full_table.loc[dataset, "NES"])],
            [y],
            marker    = "D",
            s         = 48,
            color     = "#222222",
            edgecolor = "white",
            linewidth = 0.65,
            zorder    = 5,
        )

        row = summary_table.loc[summary_table["dataset"].eq(dataset)].iloc[0]
        ax.text(
            0.04,
            y,
            f"NES {row['LOPO_NES_min']:.2f}–{row['LOPO_NES_max']:.2f}\nmax P {row['max_nominal_P']:.3g}",
            transform   = ax.get_yaxis_transform(),
            ha          = "left",
            va          = "center",
            fontsize    = 6.0,
            color       = "#555555",
            linespacing = 1.15,
        )

    x_max = max(all_values + full_table["NES"].tolist())
    ax.set_xlim(-0.05, x_max + 0.18)
    ax.set_ylim(-0.48, 1.48)
    ax.axvline(0, color="#999999", linestyle=(0, (3, 2)), linewidth=0.7)
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
    figure_style.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=5)

    fig.text(0.16, 0.96, "Leave-one-patient-out GSEA stability", fontsize=8.6, fontweight="bold", va="top")
    fig.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", markersize=4, color="#555555", label="Leave-one-out"),
            Line2D([], [], marker="D", linestyle="none", markersize=4, color="#222222", label="Full cohort"),
        ],
        loc            = "lower left",
        bbox_to_anchor = (0.16, 0.015),
        ncol           = 2,
        frameon        = False,
        fontsize       = 5.8,
        handlelength   = 0.8,
        handletextpad  = 0.4,
        columnspacing  = 1.2,
        borderaxespad  = 0,
    )
    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step46 recurrent GBM Figure 1 E-F candidates")
    manifest_row = panel_record(
        "Figure1F",
        stem,
        "Leave-one-patient-out paired pseudobulk raw20 GSEA stability",
        source,
        pdf,
        png,
    )
    outputs = (manifest_row, summary_path)
    return outputs


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path]:
    """Assemble the candidate preview at the original size and return its PDF, PNG and layout paths."""
    paths = {record["panel"]: Path(record["png"]) for record in records}
    layout_rows = [
        {"panel": "Figure1E", "left": 120, "top": 170, "right": 1750, "bottom": 1770},
        {"panel": "Figure1F", "left": 1850, "top": 170, "right": 3480, "bottom": 1770},
    ]
    canvas = Image.new("RGB", (3600, 1900), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = figure_style.image_font(54, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        figure_style.paste_contained(canvas, paths[panel], (row["left"], row["top"], row["right"], row["bottom"]))
        drawer.text((int(row["left"]), 70), panel[-1], fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1EF_candidate_preview.png"
    pdf = FIG_OUT / "Figure1EF_candidate_preview.pdf"
    source = SOURCE_OUT / "Figure1EF_candidate_preview_layout.csv"

    # Save the candidate preview as a 300 dpi PDF master before rendering its PNG.
    preview_fig = plt.figure(figsize=(12, 1900 / 300), facecolor="white")
    preview_ax = preview_fig.add_axes([0, 0, 1, 1])
    preview_ax.imshow(np.asarray(canvas), interpolation="nearest")
    preview_ax.axis("off")
    preview_fig.canvas.draw()
    preview_fig.savefig(
        pdf,
        dpi = 300,
        facecolor = "white",
        edgecolor = "none",
        metadata={
            "Creator": "Step46 recurrent GBM Figure 1 E-F candidates",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(preview_fig)
    figure_style.render_pdf_png(pdf, png, dpi=300)

    layout_table = pd.DataFrame(layout_rows)
    layout_table.to_csv(source, index=False)
    outputs = (pdf, png, source)
    return outputs


def main() -> None:
    """Build the two Figure 1 candidates and write their original manifest. Returns nothing."""
    np.random.seed(20260713)
    figure_style.apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    records = [panel_e_rank_rank()]
    panel_f, summary_path = panel_f_lopo()
    records.append(panel_f)
    preview_pdf, preview_png, preview_source = make_preview(records)
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
    manifest_table = pd.DataFrame(records)
    manifest_table["lopo_summary"] = ""
    manifest_table.loc[manifest_table["panel"].eq("Figure1F"), "lopo_summary"] = str(summary_path)
    manifest_table.to_csv(WRITE_ROOT / "Figure1EF_panel_manifest.csv", index=False)
    print(f"Figure 1E–F panels: {int(manifest_table['panel'].isin(['Figure1E', 'Figure1F']).sum())}")


if __name__ == "__main__":
    main()
