#!/usr/bin/env python3
# Purpose: Assemble Figure 1 panels A–F and generate the published legend and gene definition table.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05g_figure1.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
import importlib.util
import shutil
from pathlib import Path
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw, ImageOps
from scipy.stats import spearmanr
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
render_pdf_png = figure_style.render_pdf_png
sha256 = figure_style.sha256
image_font = figure_style.image_font
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel

if __name__ == "__main__":
    apply_publication_style()
    np.random.seed(20260714)
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write/51_figure1_final_reorganized/Figure1"
    source_out = write_root / "source_data"

    # Use the published output directory and basename.
    source_out.mkdir(parents=True, exist_ok=True)
    fig_out = root / "figures/51_figure1_final_reorganized/Figure1"
    fig_out.mkdir(parents=True, exist_ok=True)
    module_path = root / "scripts/05c_recurrence_evidence_layout.py"
    spec = importlib.util.spec_from_file_location("step42_visual_story", module_path)

    # Use the existing input checks and branch conditions.
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load the Figure 1A plot module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.fig_root = fig_out.parent
    module.write_root = write_root.parent

    module.configure_style()
    generated = module.figure1_design()
    old_stem = "Figure1A_compact_paired_design"
    stem = "Figure1A_paired_cohort_design"
    old_source = source_out / f"{old_stem}_source.csv"

    old_pdf = fig_out / f"{old_stem}.pdf"
    old_png = fig_out / f"{old_stem}.png"
    source = source_out / f"{stem}_source.csv"
    pdf = fig_out / f"{stem}.pdf"
    png = fig_out / f"{stem}.png"

    # Copy or rename the required files using the published input and output paths.
    old_source.replace(source)
    old_pdf.replace(pdf)
    render_pdf_png(pdf, png)
    old_png.unlink()
    apply_publication_style()

    # The helper below uses the existing arguments and return fields.
    def record(
        panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path, provenance: str
    ) -> dict[str, str]:
        """Return the existing manifest fields for the given panel name, files and provenance description."""
        result = {
            "panel": panel,
            "stem": stem,
            "message": message,
            "source": str(source),
            "pdf": str(pdf),
            "png": str(png),
            "provenance": provenance,
            "source_sha256": sha256(source),
            "pdf_sha256": sha256(pdf),
            "png_sha256": sha256(png),
        }
        return result

    panel_record_1 = record(
        "A",
        stem,
        "Two independent paired recurrence cohorts and patient-level Miller-IM analysis design",
        source,
        pdf,
        png,
        f"vector redraw from {generated['stem']}",
    )

    step45_fig = root / "figures/45_figure1_program_rebuild/Figure1"

    step45_source = root / "write/45_figure1_program_rebuild/Figure1/source_data"

    def copy_locked_panel(
        panel: str,
        old_stem: str,
        new_stem: str,
        message: str,
        figure_root: Path,
        source_root: Path,
        source_suffix: str,
    ) -> dict[str, str]:
        """Copy the source table and PDF, render the PNG, and return the manifest record under the new stem."""

        # Use the published output directory and basename.
        source_input = source_root / f"{old_stem}_source{source_suffix}"
        source = source_out / f"{new_stem}_source{source_suffix}"
        pdf = fig_out / f"{new_stem}.pdf"
        png = fig_out / f"{new_stem}.png"
        shutil.copy2(source_input, source)

        # Copy or rename the required files using the published input and output paths.
        shutil.copy2(figure_root / f"{old_stem}.pdf", pdf)
        render_pdf_png(pdf, png)
        result = record(panel, new_stem, message, source, pdf, png, f"locked copy of {old_stem}")
        return result


    # Use the published output directory and basename.
    panel_record_2 = copy_locked_panel(
        "B",
        "Figure1B_GSE174554_raw20_gsea",
        "Figure1B_GSE174554_raw20_gsea",
        "Formal paired raw-count pseudobulk Miller-IM GSEA in GSE174554",
        step45_fig,
        step45_source,
        ".csv.gz",
    )
    panel_record_3 = copy_locked_panel(
        "C",
        "Figure1C_GSE274546_raw20_gsea",
        "Figure1C_GSE274546_raw20_gsea",
        "Formal paired raw-count pseudobulk Miller-IM GSEA in GSE274546",
        step45_fig,
        step45_source,
        ".csv.gz",
    )
    step46_source = root / "write/46_figure1_ef_candidates/Figure1/source_data"

    input_path = step46_source / "Figure1E_cross_cohort_gene_rank_comparison_source.csv"
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
        raise ValueError(f"Figure1D is missing columns: {sorted(required - set(rank_table.columns))}")
    shared_table = rank_table.loc[rank_table["shared_leading_edge"].astype(bool)].copy()

    # Keep the prespecified 20 genes in their original order.
    raw20_other_table = rank_table.loc[
        rank_table["miller_raw20"].astype(bool) & ~rank_table["shared_leading_edge"].astype(bool)
    ].copy()
    background_table = rank_table.loc[~rank_table["miller_raw20"].astype(bool)].copy()
    if rank_table["gene"].duplicated().any() or background_table.empty:
        raise ValueError("Figure1D source must contain unique genes and a non-empty transcriptome background")
    if len(shared_table) != 7 or len(raw20_other_table) != 12:
        raise ValueError(
            f"Figure1D counts: all={len(rank_table)}, shared={len(shared_table)}, other_raw20={len(raw20_other_table)}"
        )
    expected_shared_leading_edge = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}

    # Use the existing input checks and branch conditions.
    if set(shared_table["gene"]) != expected_shared_leading_edge:
        raise ValueError(
            (
                "Figure1D shared leading-edge genes differ: expected="
                f"{sorted(expected_shared_leading_edge)}"
                ", found="
                f"{sorted(shared_table['gene'])}"
            )
        )
    labelled_genes = {"CCL4", "CH25H", "FOLR2", "PDK4", "SGK1"}
    rank_table["display_role"] = "transcriptome_background"
    rank_table.loc[rank_table["miller_raw20"].astype(bool), "display_role"] = "other_measured_raw20"
    rank_table.loc[rank_table["shared_leading_edge"].astype(bool), "display_role"] = "shared_leading_edge"

    # Prepare the cohort gene ranks using the original order and merge keys.
    rank_table["display_label"] = rank_table["gene"].where(rank_table["gene"].isin(labelled_genes), "")
    stem = "Figure1D_shared_raw20_outliers"
    source = source_out / f"{stem}_source.csv"
    csv_path = source
    rank_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig = new_figure(91, 84)
    ax = fig.add_axes([0.17, 0.23, 0.75, 0.57])
    background_table = background_table.sort_values("density")
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c=background_table["density"],
        cmap="Greys",
        s=2.0,
        alpha=0.66,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    highlight = "#E64B35"

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.scatter(
        raw20_other_table["rank_stat_GSE174554"],
        raw20_other_table["rank_stat_GSE274546"],
        s=23,
        facecolor="white",
        edgecolor=highlight,
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared_table["rank_stat_GSE174554"],
        shared_table["rank_stat_GSE274546"],
        s=29,
        color=highlight,
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )
    label_offsets = {"CCL4": (-25, 7), "CH25H": (5, 5), "SGK1": (-23, -10), "FOLR2": (5, 7), "PDK4": (5, -10)}
    for _, row in shared_table.loc[shared_table["gene"].isin(labelled_genes)].iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=label_offsets[row["gene"]],
            textcoords="offset points",
            fontsize=5.8,
            color="#B53A2B",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
            zorder=5,
        )
    ax.axhline(0, color="#A4A4A4", linewidth=0.55)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = (
        float(np.nanmax(np.abs(rank_table[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))) * 1.04
    )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    clean_axis(ax)
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

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    fig.text(
        0.17,
        0.965,
        "Shared Miller-IM leading-edge genes are\nrecurrent-enriched outliers",
        fontsize=8.3,
        fontweight="bold",
        va="top",
        linespacing=1.05,
    )
    fig.text(
        0.17,
        0.025,
        "Filled: shared GSEA leading edge   Open: other measured Miller-IM genes",
        fontsize=5.8,
        color="#666666",
        va="bottom",
    )
    pdf, png = save_panel(fig, fig_out, stem, "Step51 final recurrent GBM Figure 1")
    panel_record_4 = record(
        "D",
        stem,
        "Shared Miller-IM leading-edge genes occupy the joint recurrent-enriched region",
        source,
        pdf,
        png,
        "redrawn from Step46 all-gene rank source; no naked genome-wide correlation annotation",
    )

    # Use the panel's published stem for its source table and figure files.
    input_path = step45_source / "Figure1E_shared_leading_edge_heatmap_source.csv"
    direction_table = pd.read_csv(input_path)
    if len(direction_table) != 14 or direction_table["gene"].nunique() != 7:
        raise ValueError("Figure1E must contain 7 genes × 2 cohorts")
    if int((direction_table["FDR"] < 0.05).sum()) != 0:
        raise ValueError("Figure1E gene-level FDR results differ from the published analysis")

    # Use the published output directory and basename.
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    datasets = ["GSE174554", "GSE274546"]
    raw20_name = "Miller_Microglial_Inflammatory_raw_top20"

    def formal_gsea_statistics() -> pd.DataFrame:
        """Return Step 38 raw20 GSEA results at threshold 20 in the original cohort order."""

        # Read the frozen inputs with their original row order, dataset names and column names.
        gsea_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
        is_formal = gsea_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
        gsea_table = gsea_table.loc[
            gsea_table["dataset"].isin(datasets)
            & gsea_table["threshold"].eq(20)
            & gsea_table["pathway"].eq(raw20_name)
            & is_formal
        ].copy()
        if gsea_table.groupby("dataset").size().to_dict() != {dataset: 1 for dataset in datasets}:
            raise ValueError(
                (
                    "Step38 must contain exactly one formal Miller-IM GSEA row per coh"
                    "ort: "
                    f"{gsea_table.groupby('dataset').size().to_dict()}"
                )
            )
        pair_counts = gsea_table.set_index("dataset")["n_pairs"].astype(int).to_dict()

        # Use the existing input checks and branch conditions.
        if pair_counts != {"GSE174554": 17, "GSE274546": 45}:
            raise ValueError(f"Unexpected formal pair counts in Step38 GSEA: {pair_counts}")
        result_table = gsea_table.set_index("dataset").loc[datasets]
        return result_table

    gsea_table = formal_gsea_statistics()

    # Use the published output directory and basename.
    stem = "Figure1E_shared_leading_edge_gene_effects"
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
    ax = fig.add_axes([0.24, 0.19, 0.55, 0.6])

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    image = ax.imshow(
        effect_table.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest"
    )
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
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(datasets):
            value = float(effect_table.loc[gene, dataset])
            color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=color)

    # Process cohorts, genes or panels in the published order.
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(datasets), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    color_ax = fig.add_axes([0.82, 0.3, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)
    colorbar.outline.set_linewidth(0.45)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    fig.text(0.24, 0.95, "Shared leading-edge gene effects", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.24,
        0.085,
        "0/14 gene-by-cohort tests reached FDR < 0.05",
        fontsize=5.9,
        color="#777777",
        va="bottom",
    )
    pdf, png = save_panel(fig, fig_out, stem, "Step51 final recurrent GBM Figure 1")
    panel_record_5 = record(
        "E",
        stem,
        "Effect sizes of the seven shared leading-edge genes in both cohorts",
        source,
        pdf,
        png,
        "redrawn from Step45 gene-level effect source; redundant cohort GSEA text removed",
    )

    old_stem = "Figure1F_leave_one_patient_out_gsea_stability"
    new_stem = old_stem
    source_input = step46_source / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    source = source_out / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    shutil.copy2(source_input, source)

    # Process cohorts, genes or panels in the published order.
    for name in [
        "Figure1F_leave_one_patient_out_summary.csv",
        "Figure1F_full_data_reproduction_check.csv",
        "Figure1F_leave_one_patient_out_fold_audit.csv",
        "Figure1F_leave_one_patient_out_program_gsea.csv",
    ]:
        shutil.copy2(step46_source / name, source_out / name)
    pdf = fig_out / f"{new_stem}.pdf"
    png = fig_out / f"{new_stem}.png"
    step46_fig = root / "figures/46_figure1_ef_candidates/Figure1"
    shutil.copy2(step46_fig / f"{old_stem}.pdf", pdf)

    render_pdf_png(pdf, png)
    panel_record_6 = record(
        "F",
        new_stem,
        "Leave-one-patient-out paired pseudobulk Miller-IM GSEA stability",
        source,
        pdf,
        png,
        "locked copy of Step46 Figure1F",
    )
    records = [panel_record_1, panel_record_2, panel_record_3, panel_record_4, panel_record_5, panel_record_6]
    if [row["panel"] for row in records] != list("ABCDEF"):
        raise ValueError("Figure1 must have exactly one record for each panel A–F")

    # Use the existing input checks and branch conditions.
    if len({row["stem"] for row in records}) != 6:
        raise ValueError("Figure1 contains duplicate filenames")
    provenance_table = pd.read_csv(step38 / "fixed_program_provenance.csv")
    raw20_row = provenance_table.loc[provenance_table["signature"].eq(raw20_name)].iloc[0]
    genes = str(raw20_row["genes"]).split(";")
    if len(genes) != 20:
        raise ValueError(f"Expected 20 raw20 genes, found {len(genes)}")

    # Read the frozen inputs with their original row order, dataset names and column names.
    direction_table = pd.read_csv(step38 / "independent_fixed_program_gene_direction.csv")
    direction_table = direction_table.loc[
        direction_table["dataset"].isin(datasets)
        & direction_table["threshold"].eq(20)
        & direction_table["signature"].eq(raw20_name)
    ].copy()
    gsea_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
    gsea_table = gsea_table.loc[
        gsea_table["dataset"].isin(datasets)
        & gsea_table["threshold"].eq(20)
        & gsea_table["pathway"].eq(raw20_name)
    ].copy()
    leading = {row["dataset"]: set(str(row["leadingEdge"]).split(";")) for _, row in gsea_table.iterrows()}

    # Highlight the intersection of the two cohort leading edges and apply the original check.
    shared = set.intersection(*(leading[dataset] for dataset in datasets))
    if shared != expected_shared_leading_edge:
        raise ValueError(
            f"Shared leading-edge genes differ: expected={sorted(expected_shared_leading_edge)}, found={sorted(shared)}"
        )
    rows: list[dict[str, object]] = []
    direction_table = direction_table.set_index(["dataset", "gene"])
    for order, gene in enumerate(genes, start=1):
        for dataset in datasets:
            row = direction_table.loc[dataset, gene]
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

    # Build the table from the named records in their original order and with the same unit per row.
    definition_table = pd.DataFrame(rows)
    if len(definition_table) != 40:
        raise ValueError("The raw20 definition table must have 40 rows")
    tested_counts = definition_table.groupby("dataset")["tested"].sum().astype(int).to_dict()
    if tested_counts != {"GSE174554": 19, "GSE274546": 19}:
        raise ValueError(f"Unexpected number of measurable raw20 genes: {tested_counts}")
    if set(definition_table.loc[~definition_table["tested"], "gene"]) != {"AC253572.2"}:
        raise ValueError("AC253572.2 must be the only unmeasurable raw20 gene")

    # Use the existing input checks and branch conditions.
    if int(definition_table["shared_leading_edge"].sum()) != 14:
        raise ValueError("The shared leading edge must contain 7 genes × 2 cohorts")
    path = source_out / "Figure1_raw20_definition_source.csv"
    csv_path = path
    definition_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Keep the prespecified 20 genes in their original order.
    raw20_path = path
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout_rows = [
        {"panel": "A", "left": 180, "top": 120, "right": 3420, "bottom": 820},
        {"panel": "B", "left": 160, "top": 940, "right": 1700, "bottom": 1940},
        {"panel": "C", "left": 1900, "top": 940, "right": 3440, "bottom": 1940},
        {"panel": "D", "left": 80, "top": 2110, "right": 1160, "bottom": 3230},
        {"panel": "E", "left": 1260, "top": 2110, "right": 2340, "bottom": 3230},
        {"panel": "F", "left": 2440, "top": 2110, "right": 3520, "bottom": 3230},
    ]
    canvas = Image.new("RGB", (3600, 3320), "white")
    drawer = ImageDraw.Draw(canvas)

    # Apply the original labels and placement rules.
    label_font = image_font(58, bold=True)
    for item in layout_rows:
        panel = item["panel"]
        paste_contained(canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"]))
        drawer.text((item["left"], item["top"] - 60), panel, fill="#111111", font=label_font)
    png = fig_out / "Figure1_final_preview.png"
    gray = fig_out / "Figure1_final_preview_grayscale.png"
    pdf = fig_out / "Figure1_final_preview.pdf"

    # Use the published output directory and basename.
    layout_path = source_out / "Figure1_final_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(180, 180), optimize=True)
    ImageOps.grayscale(canvas).save(gray, format="PNG", dpi=(180, 180), optimize=True)
    canvas.save(pdf, format="PDF", resolution=180.0)

    # Check the color, grayscale and PDF review previews before recording their paths.
    for preview_path in [pdf, png, gray]:
        if not preview_path.is_file() or preview_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty preview: {preview_path.name}")
    output_table = pd.DataFrame(layout_rows)

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = layout_path
    output_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    (preview_pdf, preview_png, preview_gray, preview_layout) = (pdf, png, gray, layout_path)

    # Use the published output directory and basename.
    manifest_path = write_root / "Figure1_panel_manifest.csv"
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)

    manifest_table["preview_layout"] = str(preview_layout)
    csv_path = manifest_path
    manifest_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    raw20_table = pd.read_csv(raw20_path)
    raw20 = raw20_table.sort_values("raw20_order")["gene"].drop_duplicates().tolist()

    # Keep the prespecified 20 genes in their original order.
    raw20_text = ", ".join(raw20)
    gsea_table = formal_gsea_statistics()
    rank_table = pd.read_csv(source_out / "Figure1D_shared_raw20_outliers_source.csv")
    shared_count = int(rank_table["shared_leading_edge"].astype(bool).sum())
    other_count = int(
        (rank_table["miller_raw20"].astype(bool) & ~rank_table["shared_leading_edge"].astype(bool)).sum()
    )

    # Use the existing input checks and branch conditions.
    if (shared_count, other_count) != (7, 12):
        raise ValueError(
            f"Unexpected measured Miller-IM partition in Figure1D: shared={shared_count}, other={other_count}"
        )
    rho, _ = spearmanr(rank_table["rank_stat_GSE174554"], rank_table["rank_stat_GSE274546"])
    label_order = ["CCL4", "CH25H", "FOLR2", "PDK4", "SGK1"]
    labelled = [
        gene
        for gene in label_order
        if gene in set(rank_table.loc[rank_table["display_label"].ne(""), "display_label"])
    ]
    if labelled != label_order or "CCL3" in labelled:
        raise ValueError(f"Unexpected Figure1D display labels: {labelled}")

    # Read the frozen inputs with their original row order, dataset names and column names.
    effect_table = pd.read_csv(source_out / "Figure1E_shared_leading_edge_gene_effects_source.csv")
    n_gene_tests = len(effect_table)
    significant_gene_tests = int((effect_table["FDR"] < 0.05).sum())
    if (effect_table["gene"].nunique(), n_gene_tests, significant_gene_tests) != (7, 14, 0):
        raise ValueError(
            (
                "Unexpected Figure1E test audit: genes="
                f"{effect_table['gene'].nunique()}"
                ", tests="
                f"{n_gene_tests}"
                ", FDR<0.05="
                f"{significant_gene_tests}"
            )
        )
    lopo_table = pd.read_csv(source_out / "Figure1F_leave_one_patient_out_raw20_gsea.csv")

    # Summarize leave-one-patient-out results with the published denominator and nominal P threshold.
    lopo_table = lopo_table.loc[lopo_table["run_type"].eq("Leave-one-patient-out")].copy()
    lopo_counts = lopo_table.groupby("dataset").size().astype(int).to_dict()
    expected_counts = gsea_table["n_pairs"].astype(int).to_dict()
    if lopo_counts != expected_counts:
        raise ValueError(f"Unexpected Figure1F LOPO counts: {lopo_counts}")
    positive_counts = (
        lopo_table.assign(positive=lopo_table["NES"] > 0)
        .groupby("dataset")["positive"]
        .sum()
        .astype(int)
        .to_dict()
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

    # Format labels as in the published English text for zero, one or multiple terms.
    if not labelled:
        labelled_text = "No genes"
    elif len(labelled) == 1:
        labelled_text = labelled[0]
    else:
        labelled_text = f"{', '.join(labelled[:-1])}, and {labelled[-1]}"

    def legend_p(value: float) -> str:
        """Return the nominal P value formatted for the published legend."""

        # Use the existing input checks and branch conditions.
        if value < 0.0001:
            result = "nominal *P* < 1 × 10⁻⁴"
            return result
        formatted = f"{value:.5f}".rstrip("0").rstrip(".")
        result = f"nominal *P* = {formatted}"
        return result

    content = (
        "# Figure 1 legend\n\n**Figure 1 | The 20-gene Miller-derived inflam"
        "matory microglial program (Miller-IM program) is enriched at recu"
        "rrence in two independent paired GBM cohorts.**\n\n**A,** Study des"
        "ign. GSE174554 and GSE274546 were reconstructed independently, th"
        "eir pan-myeloid compartments were reclustered separately, and rec"
        "urrence was tested using paired patient-level analyses ("
        f"{pair_174554}"
        " and "
        f"{pair_274546}"
        " matched IDH-wild-type patients, respectively). **B-C,** Preranke"
        "d GSEA of the fixed Miller-IM program using paired raw-count pseu"
        "dobulk differential-expression ranks. Positive normalized enrichm"
        "ent scores indicate recurrent enrichment. GSE174554: NES = "
        f"{float(gsea_table.loc['GSE174554', 'NES']):.2f}"
        ", "
        f"{legend_p(float(gsea_table.loc['GSE174554', 'pval']))}"
        "; GSE274546: NES = "
        f"{float(gsea_table.loc['GSE274546', 'NES']):.2f}"
        ", "
        f"{legend_p(float(gsea_table.loc['GSE274546', 'pval']))}"
        ". **D,** Genome-wide paired rank statistics in the two cohorts. F"
        "illed red symbols mark the "
        f"{shared_count}"
        " measurable Miller-IM genes present in both cohort-specific GSEA "
        "leading edges; open red symbols mark the other "
        f"{other_count}"
        " measurable Miller-IM genes. Only "
        f"{labelled_text}"
        " are labelled to avoid over-annotation. Because filled-symbol mem"
        "bership is defined by the intersection of the two leading edges, "
        "their joint positive-quadrant position is descriptive and partly "
        "expected by construction; this panel is not an independent single"
        "-gene replication test. Genome-wide rank correlation was low (Spe"
        "arman rₛ = "
        f"{rho:.2f}"
        ") and is provided here only as background context, not as the pan"
        "el claim. **E,** Recurrent-versus-primary log₂ fold changes of th"
        "e "
        f"{shared_count}"
        " shared leading-edge genes. "
        f"{significant_gene_tests}"
        "/"
        f"{n_gene_tests}"
        " gene-by-cohort tests reached FDR < 0.05, supporting interpretati"
        "on at the program rather than single-gene level. **F,** Leave-one"
        "-patient-out sensitivity analysis of the paired pseudobulk Miller"
        "-IM GSEA. Normalized enrichment scores remained positive in "
        f"{positive_counts['GSE174554']}"
        "/"
        f"{lopo_counts['GSE174554']}"
        " and "
        f"{positive_counts['GSE274546']}"
        "/"
        f"{lopo_counts['GSE274546']}"
        " omissions; nominal *P* remained < 0.05 in "
        f"{significant_counts['GSE174554']}"
        "/"
        f"{lopo_counts['GSE174554']}"
        " and "
        f"{significant_counts['GSE274546']}"
        "/"
        f"{lopo_counts['GSE274546']}"
        " omissions, respectively. The maximum nominal *P* was "
        f"{max_p:.4g}"
        " in "
        f"{max_p_dataset}"
        ". Diamonds show full-data estimates; filled circles denote nomina"
        "l *P* < 0.05.\n\nThe fixed Miller-IM gene set was defined a priori "
        "from Miller et al. Supplementary Table 2: "
        f"{raw20_text}"
        ". Nineteen of 20 genes were measurable in each cohort; AC253572.2"
        " was not measured. Statistical unit: patient. Panel D is a descri"
        "ptive cross-cohort background view; panels B, C, and F provide th"
        "e formal program-level recurrence evidence.\n"
    )
    path = write_root / "Figure1_legend.md"

    path.write_text(content, encoding="utf-8")
    if path.read_text(encoding="utf-8") != content:
        raise ValueError("Figure1 legend content was not saved completely")
    legend_path = path
    content = (
        "# Figure 1 final reorganization\n\n- Final panel mapping: A design;"
        " B-C independent Miller-IM GSEA; D all-gene background with share"
        "d leading-edge Miller-IM genes; E shared leading-edge gene effect"
        "s; F leave-one-patient-out stability.\n- D contains no naked genom"
        "e-wide correlation annotation and labels only CCL4, CH25H, FOLR2,"
        " PDK4, and SGK1.\n- D and E are unique files; the historical Step4"
        "6 duplicate letter assignment is not carried forward.\n- All six i"
        "ndependent panels are letter-free. Letters occur only in the revi"
        "ew preview.\n- Manifest: `"
        f"{manifest_path}"
        "`\n- Legend: `"
        f"{legend_path}"
        "`\n- Internal raw20 version definition: `"
        f"{raw20_path}"
        "`\n"
    )
    final_result_path = write_root / "FINAL_RESULT.md"
    final_result_path.write_text(content, encoding="utf-8")
    if final_result_path.read_text(encoding="utf-8") != content:
        raise ValueError("Figure1 final-result note was not saved completely")
