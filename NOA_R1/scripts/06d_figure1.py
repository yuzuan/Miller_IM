#!/usr/bin/env python3
# Purpose: generate the six October Figure 1 panels and assemble the manuscript figure preview.
# Inputs: prepared Step 38, 45 and 46 tables under NOA_LEGACY_ROOT (cohort and threshold, GSEA, gene effects,
# leading edges and leave-one-patient-out folds), the Step 46 panel F PDF and the bundled Step 1 camera results.
# Outputs: _panels/51_figure1_final_reorganized/Figure1 contains the six A–F PDFs and PNGs, plus
# Figure1_final_preview.pdf, Figure1_final_preview.png and Figure1_final_preview_grayscale.png.
# _figure_build/51_figure1_final_reorganized/Figure1 contains Figure1_panel_manifest.csv, Figure1_legend.md,
# FINAL_RESULT.md, and source_data/ with the panel source tables, Figure1_raw20_definition_source.csv and
# Figure1_final_preview_layout.csv. Panel F retains its five existing leave-one-patient-out source CSVs.
# Run: python3 NOA_R1/scripts/06d_figure1.py
# Keep the Matplotlib/Pillow implementation and recorded October layouts for reproduction.

from __future__ import annotations

import importlib
import shutil
from pathlib import Path

import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw, ImageOps
from scipy.stats import spearmanr

figure_style = importlib.import_module("00a_functions")
apply_publication_style = figure_style.apply_publication_style
revision_root = Path(__file__).resolve().parents[1]
legacy_project_root = figure_style.legacy_root(revision_root)

# The final exporter redirects these paths and may use the original Step 46 source attribute.
WRITE_ROOT = revision_root / "_figure_build/51_figure1_final_reorganized/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "_panels/51_figure1_final_reorganized/Figure1"
STEP46_SOURCE = legacy_project_root / "write/46_figure1_ef_candidates/Figure1/source_data"

# The two paired cohorts and shared leading-edge genes are fixed for the October source tables.
datasets = ["GSE174554", "GSE274546"]
raw20_name = "Miller_Microglial_Inflammatory_raw_top20"
expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}


def formal_gsea_statistics() -> pd.DataFrame:
    """Return the two frozen GSEA results in cohort order after checking their pair counts."""
    step38_root = legacy_project_root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    gsea_table = pd.read_csv(step38_root / "independent_fixed_program_targeted_gsea.csv")
    is_formal = gsea_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
    gsea_table = gsea_table.loc[
        gsea_table["dataset"].isin(datasets)
        & gsea_table["threshold"].eq(20)
        & gsea_table["pathway"].eq(raw20_name)
        & is_formal
    ].copy()
    if gsea_table.groupby("dataset").size().to_dict() != {dataset: 1 for dataset in datasets}:
        raise ValueError(
            "Step38 must contain exactly one formal MG-Inflam20 GSEA row per cohort: "
            f"{gsea_table.groupby('dataset').size().to_dict()}"
        )
    pair_counts = gsea_table.set_index("dataset")["n_pairs"].astype(int).to_dict()
    if pair_counts != {"GSE174554": 17, "GSE274546": 45}:
        raise ValueError(f"Unexpected formal pair counts in Step38 GSEA: {pair_counts}")
    formal_table = gsea_table.set_index("dataset").loc[datasets]
    return formal_table


def legend_p(value: float) -> str:
    """Return the legend's formatted nominal GSEA P value."""
    if value < 1e-4:
        annotation = "nominal *P* < 1 × 10⁻⁴"
    else:
        formatted = f"{value:.5f}".rstrip("0").rstrip(".")
        annotation = f"nominal *P* = {formatted}"
    return annotation


def join_labels(labels: list[str]) -> str:
    """Return the existing comma-separated legend format for labels."""
    if not labels:
        annotation = "No genes"
    elif len(labels) == 1:
        annotation = labels[0]
    else:
        annotation = f"{', '.join(labels[:-1])}, and {labels[-1]}"
    return annotation


def flatten_png(source: Path, destination: Path) -> None:
    """Write a white-background copy of a PDF-derived panel raster. Returns nothing."""
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        white.convert("RGB").save(destination, format="PNG", dpi=(600, 600), optimize=True)


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
    provenance: str,
) -> dict[str, str]:
    """Return the original Figure 1 manifest fields after checking the output files."""
    for path in (source, pdf, png):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty Figure 1 output: {path.name}")
    manifest_row = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "provenance": provenance,
        "source_sha256": figure_style.sha256(source),
        "pdf_sha256": figure_style.sha256(pdf),
        "png_sha256": figure_style.sha256(png),
    }
    return manifest_row


def panel_a(*, final_layout: bool = False) -> dict[str, str]:
    """Redraw the design from the paired-cohort audit and return the panel manifest row."""
    design = importlib.import_module("06a_figure1_evidence_layout")
    design.FIG_ROOT = FIG_OUT.parent
    design.WRITE_ROOT = WRITE_ROOT.parent
    design.configure_style()
    generated = design.figure1_design(final_layout=final_layout)

    # Rename the generated source and PDF to the final panel's published stem.
    old_stem = "Figure1A_compact_paired_design"
    stem = "Figure1A_paired_cohort_design"
    old_source = SOURCE_OUT / f"{old_stem}_source.csv"
    old_pdf = FIG_OUT / f"{old_stem}.pdf"
    old_png = FIG_OUT / f"{old_stem}.png"
    source = SOURCE_OUT / f"{stem}_source.csv"
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    old_source.replace(source)
    old_pdf.replace(pdf)
    flatten_png(old_png, png)
    old_png.unlink()
    apply_publication_style()
    manifest_row = record(
        "A",
        stem,
        "Two independent paired recurrence cohorts and patient-level MG-Inflam20 analysis design",
        source,
        pdf,
        png,
        f"vector redraw from {generated['stem']}",
    )
    return manifest_row


def replot_gsea() -> list[dict[str, str]]:
    """Draw the two GSEA panels with the final figure's shared rank and enrichment scales."""
    programs = importlib.import_module("06b_figure1_program_panels")
    programs.SOURCE_OUT = SOURCE_OUT
    programs.FIG_OUT = FIG_OUT
    curves_table_total, statistics_table = programs.load_gsea_inputs()

    # The same scale includes ranks from both cohorts. The enrichment view is fixed for Figure 1 B and C.
    rank_limit = max(abs(curves_table_total.rank_stat.min()), abs(curves_table_total.rank_stat.max())) * 1.05
    records = []
    for dataset, panel in [("GSE174554", "Figure1B"), ("GSE274546", "Figure1C")]:
        row = programs.panel_gsea(dataset, panel, curves_table_total, statistics_table, (-0.15, 0.85), rank_limit)
        row["panel"] = panel[-1]
        records.append(row)
    return records


def panel_d() -> dict[str, str]:
    """Draw the paired rank comparison and return the panel manifest row."""
    stem = "Figure1D_shared_raw20_outliers"
    input_path = STEP46_SOURCE / "Figure1E_cross_cohort_gene_rank_comparison_source.csv"
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
        raise ValueError(f"Figure 1D is missing source columns: {sorted(required - set(rank_table.columns))}")

    shared_table = rank_table.loc[rank_table["shared_leading_edge"].astype(bool)].copy()
    raw20_other_table = rank_table.loc[
        rank_table["miller_raw20"].astype(bool) & ~rank_table["shared_leading_edge"].astype(bool)
    ].copy()
    background_table = rank_table.loc[~rank_table["miller_raw20"].astype(bool)].copy()
    if rank_table["gene"].duplicated().any() or background_table.empty:
        raise ValueError("Figure1D source must contain unique genes and a non-empty transcriptome background")
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            f"Unexpected Figure 1D row counts: all={len(rank_table)}, "
            f"shared={len(shared_table)}, other_raw20={len(raw20_other_table)}"
        )
    if set(shared_table["gene"]) != expected_shared_leading_edge:
        raise ValueError(
            "Unexpected Figure 1D shared leading-edge genes: "
            f"expected={sorted(expected_shared_leading_edge)}, found={sorted(shared_table['gene'])}"
        )

    # These five genes are the displayed labels, leaving the remaining measured genes as points.
    labelled_genes = {"CCL4", "CH25H", "FOLR2", "PDK4", "SGK1"}
    rank_table["display_role"] = "transcriptome_background"
    rank_table.loc[rank_table["miller_raw20"].astype(bool), "display_role"] = "other_measured_raw20"
    rank_table.loc[rank_table["shared_leading_edge"].astype(bool), "display_role"] = "shared_leading_edge"
    rank_table["display_label"] = rank_table["gene"].where(rank_table["gene"].isin(labelled_genes), "")
    source = SOURCE_OUT / f"{stem}_source.csv"
    rank_table.to_csv(source, index=False)

    fig = figure_style.new_figure(91, 70)
    ax = fig.add_axes([0.17, 0.23, 0.75, 0.57])
    background_table = background_table.sort_values("density")
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c      = background_table["density"],
        cmap   = "Greys",
        s      = 2.0,
        alpha  = 0.66,
        linewidths = 0,
        rasterized = True,
        zorder = 1,
    )
    highlight = "#E64B35"
    ax.scatter(
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s         = 23,
        facecolor = "white",
        edgecolor = highlight,
        linewidth = 0.85,
        zorder    = 3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s         = 29,
        color     = highlight,
        edgecolor = "white",
        linewidth = 0.65,
        zorder    = 4,
    )

    label_offsets = {
        "CCL4": (-43, 17),
        "CH25H": (10, 13),
        "SGK1": (-43, -17),
        "FOLR2": (15, 0),
        "PDK4": (15, -18),
    }
    for _, row in shared_table.loc[shared_table["gene"].isin(labelled_genes)].iterrows():
        label_position = (
            dict(
                xytext = (0.07, 0.76),
                textcoords = "axes fraction",
                arrowprops=dict(
                    arrowstyle = "-",
                    color      = "#999999",
                    lw         = 0.55,
                    connectionstyle = "angle,angleA=0,angleB=90,rad=0",
                    shrinkA    = 3,
                    shrinkB    = 5,
                ),
            )
            if row["gene"] == "SGK1"
            else dict(xytext=label_offsets[row["gene"]], textcoords="offset points")
        )
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            **label_position,
            fontsize = 5.8,
            color    = "#B53A2B",
            path_effects = [pe.withStroke(linewidth=1.8, foreground="white")],
            zorder   = 5,
        )

    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = float(np.nanmax(np.abs(rank_table[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    figure_style.clean_axis(ax)
    fig.text(
        0.17,
        0.965,
        "Shared MG-Inflam20 leading-edge genes",
        fontsize    = 8.3,
        fontweight  = "bold",
        va          = "top",
        linespacing = 1.05,
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                marker     = "o",
                linestyle  = "none",
                markersize = 4,
                markerfacecolor = highlight,
                markeredgecolor = highlight,
                label      = "Shared leading edge",
            ),
            Line2D(
                [],
                [],
                marker     = "o",
                linestyle  = "none",
                markersize = 4,
                markerfacecolor = "white",
                markeredgecolor = highlight,
                label      = "Other program genes",
            ),
        ],
        loc            = "lower left",
        bbox_to_anchor = (0.17, 0.007),
        ncol           = 2,
        frameon        = False,
        fontsize       = 5.8,
        handlelength   = 0.8,
        handletextpad  = 0.4,
        columnspacing  = 0.8,
        borderaxespad  = 0,
    )

    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step51 final recurrent GBM Figure 1")
    manifest_row = record(
        "D",
        stem,
        "Shared MG-Inflam20 leading-edge genes occupy the joint recurrent-enriched region",
        source,
        pdf,
        png,
        "redrawn from Step46 all-gene rank source; no naked genome-wide correlation annotation",
    )
    return manifest_row


def panel_e(*, final_layout: bool = False) -> dict[str, str]:
    """Draw the gene-level effects of the shared leading edge and return its manifest row."""
    stem = "Figure1E_shared_leading_edge_gene_effects"
    step45_source = legacy_project_root / "write/45_figure1_program_rebuild/Figure1/source_data"
    input_path = step45_source / "Figure1E_shared_leading_edge_heatmap_source.csv"
    direction_table = pd.read_csv(input_path)
    if len(direction_table) != 14 or direction_table["gene"].nunique() != 7:
        raise ValueError("Figure 1E must contain seven genes across two cohorts")
    if int((direction_table["FDR"] < 0.05).sum()) != 0:
        raise ValueError("Figure 1E single-gene FDR results differ from the frozen results")
    source = SOURCE_OUT / f"{stem}_source.csv"
    direction_table.to_csv(source, index=False)

    order = direction_table.sort_values("gene_order")["gene"].drop_duplicates().tolist()
    effect_table = direction_table.pivot(index="gene", columns="dataset", values="logFC").reindex(
        index=order, columns=datasets
    )
    max_abs = float(np.nanmax(np.abs(effect_table.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    fig = figure_style.new_figure(86, 70)
    ax = fig.add_axes([0.24, 0.27, 0.55, 0.53])
    image = ax.imshow(
        effect_table.to_numpy(),
        cmap   = "RdBu_r",
        norm   = norm,
        aspect = "auto",
        interpolation = "nearest",
    )
    ax.set_xticks(
        range(len(datasets)),
        [
            "GSE174554",
            "GSE274546",
        ],
        fontsize=6.8,
    )
    ax.set_yticks(range(len(order)), order, fontsize=7.1)
    ax.tick_params(length=0)
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(datasets):
            value = float(effect_table.loc[gene, dataset])
            color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(datasets), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    color_ax = fig.add_axes([0.82, 0.30, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    if final_layout:
        colorbar.set_label("Recurrent − primary log2 fold change", fontsize=6.5)
    else:
        colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)
    colorbar.outline.set_linewidth(0.45)
    fig.text(0.24, 0.95, "Shared leading-edge gene effects", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = figure_style.save_panel(fig, FIG_OUT, stem, "Step51 final recurrent GBM Figure 1")
    manifest_row = record(
        "E",
        stem,
        "Effect sizes of the seven shared leading-edge genes in both cohorts",
        source,
        pdf,
        png,
        "redrawn from Step45 gene-level effect source; redundant cohort GSEA text removed",
    )
    return manifest_row


def copy_panel_f(*, final_layout: bool = False) -> dict[str, str]:
    """Collect the frozen fold results and return the Figure 1F panel manifest row."""
    stem = "Figure1F_leave_one_patient_out_gsea_stability"
    source = SOURCE_OUT / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    shutil.copy2(STEP46_SOURCE / source.name, source)
    for filename in [
        "Figure1F_leave_one_patient_out_summary.csv",
        "Figure1F_full_data_reproduction_check.csv",
        "Figure1F_leave_one_patient_out_fold_audit.csv",
        "Figure1F_leave_one_patient_out_program_gsea.csv",
    ]:
        shutil.copy2(STEP46_SOURCE / filename, SOURCE_OUT / filename)

    # The final print export redraws F at 170 × 40 mm using the same frozen folds.
    if final_layout:
        lopo = importlib.import_module("06c_figure1_shared_genes_and_lopo")
        lopo.WRITE_ROOT = WRITE_ROOT
        lopo.SOURCE_OUT = SOURCE_OUT
        lopo.FIG_OUT = FIG_OUT
        row, _ = lopo.panel_f_lopo(final_layout=True)
        manifest_row = record(
            "F",
            row["stem"],
            "Leave-one-patient-out stability",
            source,
            Path(row["pdf"]),
            Path(row["png"]),
            "vector redraw from frozen data",
        )
        return manifest_row

    # The base layout uses the original panel F master and only regenerates its raster derivative.
    legacy_figures = legacy_project_root / "figures/46_figure1_ef_candidates/Figure1"
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    shutil.copy2(legacy_figures / pdf.name, pdf)
    figure_style.render_pdf_png(pdf, png)
    manifest_row = record(
        "F",
        stem,
        "Leave-one-patient-out paired pseudobulk MG-Inflam20 GSEA stability",
        source,
        pdf,
        png,
        "locked copy of Step46 Figure1F",
    )
    return manifest_row


def build_raw20_definition() -> Path:
    """Write the fixed program definition, cohort coverage and leading-edge membership."""
    step38_root = legacy_project_root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    provenance_table = pd.read_csv(step38_root / "fixed_program_provenance.csv")
    raw20_row = provenance_table.loc[provenance_table["signature"].eq(raw20_name)].iloc[0]
    genes = str(raw20_row["genes"]).split(";")
    if len(genes) != 20:
        raise ValueError(f"Expected 20 genes in the fixed MG-Inflam20 program, found {len(genes)}")

    direction_table = pd.read_csv(step38_root / "independent_fixed_program_gene_direction.csv")
    direction_table = direction_table.loc[
        direction_table["dataset"].isin(datasets)
        & direction_table["threshold"].eq(20)
        & direction_table["signature"].eq(raw20_name)
    ].copy()
    gsea_table = pd.read_csv(step38_root / "independent_fixed_program_targeted_gsea.csv")
    gsea_table = gsea_table.loc[
        gsea_table["dataset"].isin(datasets) & gsea_table["threshold"].eq(20) & gsea_table["pathway"].eq(raw20_name)
    ].copy()
    leading = {row["dataset"]: set(str(row["leadingEdge"]).split(";")) for _, row in gsea_table.iterrows()}
    shared = set.intersection(*(leading[dataset] for dataset in datasets))
    if shared != expected_shared_leading_edge:
        raise ValueError(
            "Unexpected shared leading-edge genes: "
            f"expected={sorted(expected_shared_leading_edge)}, found={sorted(shared)}"
        )

    rows: list[dict[str, object]] = []
    indexed_direction_table = direction_table.set_index(["dataset", "gene"])
    for order, gene in enumerate(genes, start=1):
        for dataset in datasets:
            row = indexed_direction_table.loc[(dataset, gene)]
            rows.append(
                {
                    "raw20_order": order,
                    "gene": gene,
                    "dataset": dataset,
                    "tested": bool(row["tested"]),
                    "logFC": row["logFC"],
                    "PValue": row["PValue"],
                    "FDR": row["FDR"],
                    "cohort_leading_edge": gene in leading[dataset],
                    "shared_leading_edge": gene in shared,
                }
            )
    definition_table = pd.DataFrame(rows)
    if len(definition_table) != len(genes) * len(datasets):
        raise ValueError("The program definition must contain every gene for both cohorts")
    tested_counts = definition_table.groupby("dataset")["tested"].sum().astype(int).to_dict()
    if tested_counts != {"GSE174554": 19, "GSE274546": 19}:
        raise ValueError(f"Unexpected measurable MG-Inflam20 gene counts: {tested_counts}")
    if set(definition_table.loc[~definition_table["tested"], "gene"]) != {"AC253572.2"}:
        raise ValueError("AC253572.2 must be the only unmeasurable MG-Inflam20 gene")
    if int(definition_table["shared_leading_edge"].sum()) != len(expected_shared_leading_edge) * len(datasets):
        raise ValueError("The definition must include all shared leading-edge genes in both cohorts")
    path = SOURCE_OUT / "Figure1_raw20_definition_source.csv"
    definition_table.to_csv(path, index=False)
    return path


def make_preview(records: list[dict[str, str]], *, final_layout: bool = False) -> tuple[Path, Path, Path, Path]:
    """Assemble the labeled Figure 1 preview and return the images and layout table."""
    by_panel = {row["panel"]: Path(row["png"]) for row in records}

    # The final print layout gives A more height while keeping the same 3600 × 4600 canvas.
    if final_layout:
        layout = [
            {"panel": "A", "left": 95, "top": 40, "right": 3505, "bottom": 1420},
            {"panel": "B", "left": 100, "top": 1470, "right": 1710, "bottom": 2430},
            {"panel": "C", "left": 1890, "top": 1470, "right": 3500, "bottom": 2430},
            {"panel": "D", "left": 55, "top": 2480, "right": 1745, "bottom": 3700},
            {"panel": "E", "left": 1835, "top": 2480, "right": 3555, "bottom": 3700},
            {"panel": "F", "left": 95, "top": 3760, "right": 3505, "bottom": 4560},
        ]
    else:
        layout = [
            {"panel": "A", "left": 95, "top": 40, "right": 3505, "bottom": 1240},
            {"panel": "B", "left": 100, "top": 1300, "right": 1710, "bottom": 2320},
            {"panel": "C", "left": 1890, "top": 1300, "right": 3500, "bottom": 2320},
            {"panel": "D", "left": 55, "top": 2380, "right": 1745, "bottom": 3680},
            {"panel": "E", "left": 1835, "top": 2380, "right": 3555, "bottom": 3680},
            {"panel": "F", "left": 95, "top": 3740, "right": 3505, "bottom": 4540},
        ]
    canvas = Image.new("RGB", (3600, 4600), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = figure_style.image_font(66, bold=True)
    for item in layout:
        panel = item["panel"]
        figure_style.paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], max(5, item["top"] - 60)), panel, fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1_final_preview.png"
    gray = FIG_OUT / "Figure1_final_preview_grayscale.png"
    pdf = FIG_OUT / "Figure1_final_preview.pdf"
    layout_path = SOURCE_OUT / "Figure1_final_preview_layout.csv"

    # Pillow writes the exact preview page dimensions. High-quality encoding preserves the panel raster.
    canvas.save(pdf, format="PDF", resolution=180.0, quality=100, subsampling=0)
    figure_style.render_pdf_png(pdf, png, dpi=180)
    with Image.open(png) as rendered:
        ImageOps.grayscale(rendered).save(gray, format="PNG", dpi=(180, 180), optimize=True)
    with Image.open(gray) as gray_image:
        gray_image.verify()

    layout_table = pd.DataFrame(layout)
    layout_table.to_csv(layout_path, index=False)
    preview_paths = (pdf, png, gray, layout_path)
    return preview_paths


def write_legend(raw20_path: Path) -> Path:
    """Write the original figure caption using the rendered source tables."""
    definition_table = pd.read_csv(raw20_path)
    raw20 = definition_table.sort_values("raw20_order")["gene"].drop_duplicates().tolist()
    raw20_text = ", ".join(raw20)
    gsea_table = formal_gsea_statistics()

    # Recalculate the descriptive correlation and confirm the source labels used in panel D.
    rank_table = pd.read_csv(SOURCE_OUT / "Figure1D_shared_raw20_outliers_source.csv")
    shared_count = int(rank_table["shared_leading_edge"].astype(bool).sum())
    other_count = int((rank_table["miller_raw20"].astype(bool) & ~rank_table["shared_leading_edge"].astype(bool)).sum())
    if (shared_count, other_count) != (7, 12):
        raise ValueError(
            f"Unexpected measured MG-Inflam20 partition in Figure1D: shared={shared_count}, other={other_count}"
        )
    rho, _ = spearmanr(
        rank_table["rank_stat_GSE174554"],
        rank_table["rank_stat_GSE274546"],
    )
    label_order = ["CCL4", "CH25H", "FOLR2", "PDK4", "SGK1"]
    labelled = [
        gene for gene in label_order if gene in set(rank_table.loc[rank_table["display_label"].ne(""), "display_label"])
    ]
    if labelled != label_order or "CCL3" in labelled:
        raise ValueError(f"Unexpected Figure1D display labels: {labelled}")

    gene_effect_table = pd.read_csv(SOURCE_OUT / "Figure1E_shared_leading_edge_gene_effects_source.csv")
    n_gene_tests = len(gene_effect_table)
    significant_gene_tests = int((gene_effect_table["FDR"] < 0.05).sum())
    if (gene_effect_table["gene"].nunique(), n_gene_tests, significant_gene_tests) != (7, 14, 0):
        raise ValueError(
            "Unexpected Figure1E test audit: "
            f"genes={gene_effect_table['gene'].nunique()}, tests={n_gene_tests}, FDR<0.05={significant_gene_tests}"
        )

    # Quantify the sensitivity from each patient omission using unrounded NES and P values.
    lopo_table = pd.read_csv(SOURCE_OUT / "Figure1F_leave_one_patient_out_raw20_gsea.csv")
    lopo_table = lopo_table.loc[lopo_table["run_type"].eq("Leave-one-patient-out")].copy()
    lopo_counts = lopo_table.groupby("dataset").size().astype(int).to_dict()
    expected_counts = gsea_table["n_pairs"].astype(int).to_dict()
    if lopo_counts != expected_counts:
        raise ValueError(f"Unexpected Figure1F LOPO counts: {lopo_counts}")
    positive_counts = (
        lopo_table.assign(positive=lopo_table["NES"] > 0).groupby("dataset")["positive"].sum().astype(int).to_dict()
    )
    significant_counts = (
        lopo_table.assign(significant=lopo_table["nominal_P"] < 0.05)
        .groupby("dataset")["significant"]
        .sum()
        .astype(int)
        .to_dict()
    )
    max_p_row = lopo_table.loc[lopo_table["nominal_P"].idxmax()]
    max_p = float(max_p_row["nominal_P"])
    max_p_dataset = str(max_p_row["dataset"])

    pair_174554 = int(gsea_table.loc["GSE174554", "n_pairs"])
    pair_274546 = int(gsea_table.loc["GSE274546", "n_pairs"])
    content = (
        "# Figure 1 legend\n"
        "\n"
        "**Figure 1 | The 20-gene Miller-derived inflammatory microglial program (MG-Inflam20 program) is "
        "enriched at recurrence in two independent paired GBM cohorts.**\n"
        "\n"
        "**A,** Study design. GSE174554 and GSE274546 were reconstructed independently, their pan-myeloid "
        "compartments were reclustered separately, and recurrence was tested using paired patient-level analyses"
        f" ({pair_174554} and {pair_274546} matched IDH-wild-type patients, respectively). **B-C,** Preranked "
        "GSEA of the fixed MG-Inflam20 program using paired raw-count pseudobulk differential-expression ranks. "
        "Positive normalized enrichment scores indicate recurrent enrichment. GSE174554: NES = "
        f"{float(gsea_table.loc['GSE174554', 'NES']):.2f}, {legend_p(float(gsea_table.loc['GSE174554', 'pval']))}"
        f"; GSE274546: NES = {float(gsea_table.loc['GSE274546', 'NES']):.2f}, "
        f"{legend_p(float(gsea_table.loc['GSE274546', 'pval']))}. **D,** Genome-wide paired rank statistics in "
        f"the two cohorts. Filled red symbols mark the {shared_count} measurable MG-Inflam20 genes present in "
        f"both cohort-specific GSEA leading edges; open red symbols mark the other {other_count} measurable "
        f"MG-Inflam20 genes. Only {join_labels(labelled)} are labelled to avoid over-annotation. Because "
        "filled-symbol membership is defined by the intersection of the two leading edges, their joint "
        "positive-quadrant position is descriptive and partly expected by construction; this panel is not an "
        f"independent single-gene replication test. Genome-wide rank correlation was low (Spearman rₛ = {rho:.2f}"
        ") and is provided here only as background context, not as the panel claim. **E,** "
        f"Recurrent-versus-primary log₂ fold changes of the {shared_count} shared leading-edge genes. "
        f"{significant_gene_tests}/{n_gene_tests} gene-by-cohort tests reached FDR < 0.05, supporting "
        "interpretation at the program rather than single-gene level. **F,** Leave-one-patient-out sensitivity "
        "analysis of the paired pseudobulk MG-Inflam20 GSEA. Normalized enrichment scores remained positive in "
        f"{positive_counts['GSE174554']}/{lopo_counts['GSE174554']} and {positive_counts['GSE274546']}/"
        f"{lopo_counts['GSE274546']} omissions; nominal *P* remained < 0.05 in {significant_counts['GSE174554']}/"
        f"{lopo_counts['GSE174554']} and {significant_counts['GSE274546']}/{lopo_counts['GSE274546']} omissions, "
        f"respectively. The maximum nominal *P* was {max_p:.4g} in {max_p_dataset}. Diamonds show full-data "
        "estimates; filled circles denote nominal *P* < 0.05.\n"
        "\n"
        "The fixed MG-Inflam20 gene set was defined a priori from Miller et al. Supplementary Table 2: "
        f"{raw20_text}. Nineteen of 20 genes were measurable in each cohort; AC253572.2 was not measured. "
        "Statistical unit: patient. Panel D is a descriptive cross-cohort background view; panels B, C, and F "
        "provide the formal program-level recurrence evidence.\n"
    )
    path = WRITE_ROOT / "Figure1_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest_path: Path, legend_path: Path, raw20_path: Path) -> None:
    """Write the original Figure 1 panel and file mapping. Returns nothing."""
    content = f"""# Figure 1 final reorganization

- A shows the design and B–C show independent MG-Inflam20 GSEA.
- D shows the all-gene background and shared leading-edge genes. E shows effects for the shared genes.
- F shows leave-one-patient-out stability.
- D labels CCL4, CH25H, FOLR2, PDK4, and SGK1. See the legend for the genome-wide correlation.
- The review preview adds the A–F labels to the six panel files.
- Manifest: `{manifest_path}`
- Legend: `{legend_path}`
- Internal raw20 version definition: `{raw20_path}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main(*, final_layout: bool = False) -> None:
    """Render panels A–F and write their source manifest and captions. Returns nothing."""
    apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    # Redraw A–E from the prepared cohort results. F uses the frozen leave-one-out data.
    records = [panel_a(final_layout=final_layout)]
    records.extend(replot_gsea())
    records.append(panel_d())
    records.append(panel_e(final_layout=final_layout))
    records.append(copy_panel_f(final_layout=final_layout))
    if [row["panel"] for row in records] != list("ABCDEF"):
        raise ValueError("The Figure 1 panels must appear in A–F order")
    if len({row["stem"] for row in records}) != 6:
        raise ValueError("The Figure 1 panel filenames must be unique")

    raw20_path = build_raw20_definition()
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records, final_layout=final_layout)
    manifest_path = WRITE_ROOT / "Figure1_panel_manifest.csv"
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)
    manifest_table["preview_layout"] = str(preview_layout)
    manifest_table.to_csv(manifest_path, index=False)
    legend_path = write_legend(raw20_path)
    write_result(manifest_path, legend_path, raw20_path)
    print(f"Figure 1 panels: {len(records)}")


if __name__ == "__main__":
    main()
