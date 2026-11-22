#!/usr/bin/env python3
"""Rebuild the original Figure 3 panels from paired GeoMx and Visium data.

Inputs
    Under NOA_LEGACY_ROOT, use write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/:
        fig3a_geomx22_paired_change.csv, fig3b_geomx18_gene_forest.csv, fig3e_partial_rho_main.csv,
        fig3c_gbm030_spatial_raw20_mdsc_mes.csv and fig3d_gbm049_spatial_raw20_mdsc_mes.csv.
    Also use write/51_figure1_final_reorganized/Figure1/source_data/Figure1_raw20_definition_source.csv.
    Under NOA_PUBLIC_DATA_ROOT, use metadata.csv and normalized.txt in
    write/26_pure_bioinformatics_dataset_rescue/source_metadata/GeoMx_zenodo16839828 and images
    {GSM8506649_03422,GSM8506650_25526}_tissue_hires_image.png.gz in
    write/30_dataset_rescue_search/source_metadata/GSE276841_spatial.
Outputs
    Under NOA_LEGACY_ROOT/write/53_figure3_visual_variety_rebuild/Figure3, write Figure3_panel_manifest.csv,
    Figure3_legend.md, FINAL_RESULT.md and source_data/:
        Figure3A_geomx_paired_dumbbell_source.csv, Figure3B_shared_gene_patient_fingerprint_source.csv,
        Figure3C_spatial_score_atlas_source.csv.gz, Figure3D_partial_rho_effect_plot_source.csv,
        Figure3_visual_variety_preview_layout.csv.
    Under NOA_LEGACY_ROOT/figures/53_figure3_visual_variety_rebuild/Figure3, write PDF/PNG pairs for
    Figure3A_geomx_paired_dumbbell, Figure3B_shared_gene_patient_fingerprint, Figure3C_spatial_score_atlas,
    Figure3D_partial_rho_effect_plot and Figure3_visual_variety_preview.
    Also write Figure3_visual_variety_preview_grayscale.png.
Run
    python scripts/08a_figure3_spatial_panels.py from NOA_R1. The Figure 3 R1 renderer also imports this file.
"""

from __future__ import annotations

import gzip
import importlib
import io
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps


figure_tools = importlib.import_module("00a_functions")
apply_publication_style = figure_tools.apply_publication_style
clean_axis = figure_tools.clean_axis
new_figure = figure_tools.new_figure
project_root = figure_tools.legacy_root(Path(__file__).resolve().parents[1])

# The standalone renderer uses the Step 53 paths. The R1 exporter rebinds these attributes.
WRITE_ROOT = project_root / "write/53_figure3_visual_variety_rebuild/Figure3"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = project_root / "figures/53_figure3_visual_variety_rebuild/Figure3"

# Figure 3 uses these original colors across panels, previews, and the R1 renderer.
primary = "#3C5488"
recurrent = "#E64B35"
mdsc_color = "#159D88"
mes_color = "#7E6148"
text_color = "#222222"
grid_color = "#E7E7E7"


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    """Write the PDF master and its PNG. The Figure 3 R1 renderer replaces this binding."""
    output_paths = figure_tools.save_panel(fig, FIG_OUT, stem, "Step53 recurrent GBM Figure 3")
    return output_paths


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
) -> dict[str, str]:
    """Return the legacy manifest row for one panel and its source table."""
    for path in (source, pdf, png):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 3 panel output: {path.name}")
    panel_record = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": figure_tools.sha256(source),
        "pdf_sha256": figure_tools.sha256(pdf),
        "png_sha256": figure_tools.sha256(png),
    }
    return panel_record


def patient_table() -> pd.DataFrame:
    """Return one row per IDH-wild-type GeoMx pair in the original change order."""
    path = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3a_geomx22_paired_change.csv"
    paired_table = pd.read_csv(path)
    if len(paired_table) != 22 or int((paired_table["delta"] > 0).sum()) != 16:
        raise ValueError("The frozen GeoMx cohort must contain 22 pairs with 16 recurrent-higher scores")
    paired_table = paired_table.sort_values(["delta", "Patient_No"], ascending=[False, True]).reset_index(drop=True)
    paired_table["patient_label"] = [f"P{int(patient):02d}" for patient in paired_table["Patient_No"]]
    paired_table["patient_order"] = np.arange(len(paired_table))
    paired_table["direction"] = np.where(paired_table["delta"] > 0, "Recurrent higher", "Primary higher")
    return paired_table


def panel_a(paired_table: pd.DataFrame) -> dict[str, str]:
    """Write and plot the GeoMx patients grouped by trial exposure and ordered by change."""
    stem = "Figure3A_geomx_paired_dumbbell"

    # Separate exposure groups without changing the time-point color encoding.
    display_table = paired_table.copy()
    group_counts = display_table["trial_setting"].value_counts().to_dict()
    if group_counts != {"Nivolumab": 12, "Control": 10}:
        raise ValueError(f"Unexpected trial group counts for Figure 3A: {group_counts}")
    display_table["display_group"] = display_table["trial_setting"].map({"Control": 0, "Nivolumab": 1})
    display_table = display_table.sort_values(["display_group", "delta", "Patient_No"], ascending=[True, False, True])
    source = SOURCE_OUT / f"{stem}_source.csv"
    display_table.to_csv(source, index=False)

    # The fixed y positions preserve the gap and labels for the two frozen exposure groups.
    fig = new_figure(92, 92)
    ax = fig.add_axes([0.30, 0.20, 0.64, 0.69])
    y = np.arange(len(display_table), dtype=float) + display_table["display_group"].to_numpy()
    ax.axhspan(10.5, 22.5, color="#F3F4F6", zorder=0)
    ax.axhline(10.0, color="#D0D0D0", linewidth=0.6)
    for row_y, group in [(4.5, "Control"), (16.5, "Nivolumab")]:
        ax.text(
            -0.36,
            row_y,
            f"{group}\n(n = {group_counts[group]})",
            transform = ax.get_yaxis_transform(),
            rotation  = 90,
            ha        = "center",
            va        = "center",
            fontsize  = 7.4,
            color     = text_color,
        )

    # Connect primary and recurrent scores for each patient and mark each time point.
    for row_y, row in zip(y, display_table.itertuples()):
        ax.plot(
            [row.Primary, row.Recurrence],
            [row_y, row_y],
            color     = "#A6A6A6",
            linewidth = 1.15,
            alpha     = 0.72,
            zorder    = 1,
        )
    ax.scatter(
        display_table["Primary"],
        y,
        s         = 22,
        marker    = "o",
        color     = primary,
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Primary",
    )
    ax.scatter(
        display_table["Recurrence"],
        y,
        s         = 22,
        marker    = "s",
        color     = recurrent,
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Recurrence",
    )

    # The additional pooled-mean row uses the same time-point markers.
    mean_y = len(display_table) + 2.2
    primary_mean = float(display_table["Primary"].mean())
    recurrent_mean = float(display_table["Recurrence"].mean())
    ax.plot(
        [primary_mean, recurrent_mean],
        [mean_y, mean_y],
        color     = text_color,
        linewidth = 1.5,
        zorder    = 2,
    )
    ax.scatter(primary_mean, mean_y, s=34, marker="o", color=primary, edgecolor="white", linewidth=0.6, zorder=4)
    ax.scatter(recurrent_mean, mean_y, s=34, marker="s", color=recurrent, edgecolor="white", linewidth=0.6, zorder=4)

    ax.set_yticks(
        list(y) + [mean_y],
        display_table["patient_label"].tolist() + ["Pooled mean"],
        fontsize = 5.6,
    )
    ax.set_ylim(mean_y + 1.0, -1.0)
    values = display_table[["Primary", "Recurrence"]].to_numpy(dtype=float)
    margin = max(0.12, float(np.ptp(values)) * 0.08)
    ax.set_xlim(float(np.nanmin(values) - margin), float(np.nanmax(values) + margin))
    ax.set_xlabel("Patient-level MG-Inflam20 score")
    ax.grid(axis="x", color=grid_color, linewidth=0.5)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc            = "lower center",
        bbox_to_anchor = (0.58, 0.035),
        frameon        = False,
        fontsize       = 5.8,
        ncol           = 2,
        handletextpad  = 0.4,
        columnspacing  = 0.8,
    )

    fig.text(0.08, 0.985, "Paired IBA1+ tissue scores", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = record("A", stem, "Patient-level paired GeoMx MG-Inflam20 scores", source, pdf, png)
    return panel_record


def shared_genes() -> list[str]:
    """Return the prespecified shared leading-edge genes in the original Figure 1 order."""
    path = project_root / "write/51_figure1_final_reorganized/Figure1/source_data/Figure1_raw20_definition_source.csv"
    definition_table = pd.read_csv(path)
    shared_table = definition_table.loc[definition_table["shared_leading_edge"].astype(bool)]
    shared_table = shared_table.sort_values("raw20_order")
    genes = shared_table["gene"].drop_duplicates().tolist()
    if genes != ["PDK4", "SGK1", "CH25H", "SIGLEC8", "KLF6", "FOLR2", "CCL4"]:
        raise ValueError(f"Unexpected Figure 1 shared leading-edge gene order: {genes}")
    return genes


def panel_b(paired_table: pd.DataFrame) -> dict[str, str]:
    """Rebuild the GeoMx patient-by-gene change table and render panel B."""
    genes = shared_genes()
    geomx_root = figure_tools.public_data_root(project_root) / (
        "write/26_pure_bioinformatics_dataset_rescue/source_metadata/GeoMx_zenodo16839828"
    )

    # Keep quality-controlled IBA1-positive regions from primary or recurrent IDH-wild-type tumors.
    metadata_table = pd.read_csv(geomx_root / "metadata.csv")
    eligible_table = metadata_table.loc[
        metadata_table["iba1"].eq(True)
        & metadata_table["QCFlags"].isna()
        & metadata_table["IDH_status"].eq("IDH_WT")
        & metadata_table["tumor_setting"].isin(["Primary", "Recurrence"])
    ]

    # Average log2(normalized expression + 1) across regions within patient and time point.
    expression_table = pd.read_csv(geomx_root / "normalized.txt", sep="\t", index_col=0)
    measured_genes = [gene for gene in genes if gene in expression_table.index]
    missing_genes = [gene for gene in genes if gene not in expression_table.index]
    if measured_genes != genes[:-1] or missing_genes != ["CCL4"]:
        raise ValueError(f"Unexpected GeoMx coverage: measured={measured_genes}, missing={missing_genes}")
    regions = eligible_table["Unnamed: 0"].tolist()
    log_expression_table = np.log2(expression_table.loc[measured_genes, regions] + 1.0)
    region_table = eligible_table.set_index("Unnamed: 0").loc[log_expression_table.columns]
    region_expression_table = log_expression_table.T.assign(
        Patient_No    = region_table["Patient_No"].to_numpy(),
        tumor_setting = region_table["tumor_setting"].to_numpy(),
    )
    by_patient = region_expression_table.groupby(["Patient_No", "tumor_setting"], observed=True)
    patient_gene_table = by_patient[measured_genes].mean().reset_index()

    # Keep one row per patient and prespecified gene, including missing CCL4.
    patient_ids = paired_table["Patient_No"].astype(int).tolist()
    rows: list[dict[str, object]] = []
    for gene_order, gene in enumerate(genes):
        measured = gene in measured_genes
        gene_deltas = {}
        if measured:
            timepoint_table = patient_gene_table.pivot(index="Patient_No", columns="tumor_setting", values=gene)
            timepoint_table = timepoint_table.reindex(index=patient_ids, columns=["Primary", "Recurrence"])
            missing_pair = timepoint_table.isna().any(axis=1)
            if missing_pair.any():
                raise ValueError(f"GeoMx patients {missing_pair.index[missing_pair].tolist()} lack pairs for {gene}")
            gene_deltas = (timepoint_table["Recurrence"] - timepoint_table["Primary"]).to_dict()
        for patient_order, patient_id in enumerate(patient_ids):
            rows.append(
                {
                    "Patient_No": patient_id,
                    "patient_label": f"P{patient_id:02d}",
                    "patient_order": patient_order,
                    "gene": gene,
                    "gene_order": gene_order,
                    "measured": measured,
                    "gene_delta": gene_deltas.get(patient_id, np.nan),
                }
            )
    fingerprint_table = pd.DataFrame(rows)
    program_deltas = paired_table.set_index("Patient_No")["delta"].to_dict()
    fingerprint_table["raw20_delta"] = fingerprint_table["Patient_No"].map(program_deltas)

    # Keep the prior confidence limits and full-program FDR for each measured gene.
    summary_path = project_root / (
        "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3b_geomx18_gene_forest.csv"
    )
    gene_summary_table = pd.read_csv(summary_path)
    gene_summary_table = gene_summary_table.set_index("gene")
    fingerprint_table["mean_gene_delta"] = fingerprint_table["gene"].map(gene_summary_table["mean_log2_delta"])
    fingerprint_table["ci_low"] = fingerprint_table["gene"].map(gene_summary_table["ci_low"])
    fingerprint_table["ci_high"] = fingerprint_table["gene"].map(gene_summary_table["ci_high"])
    fingerprint_table["gene_fdr"] = fingerprint_table["gene"].map(gene_summary_table["fdr"])
    summary_fields = ["mean_gene_delta", "ci_low", "ci_high", "gene_fdr"]
    fingerprint_table.loc[~fingerprint_table["measured"], summary_fields] = np.nan

    # The symmetric 98th-percentile limit caps the heatmap colors only.
    measured_values = fingerprint_table.loc[fingerprint_table["measured"], "gene_delta"].abs().to_numpy(dtype=float)
    color_limit = float(np.quantile(measured_values, 0.98))
    fingerprint_table["color_limit"] = color_limit
    if len(fingerprint_table) != 154 or int(fingerprint_table["measured"].sum()) != 132:
        raise ValueError("Figure 3B needs 22 patients, seven prespecified genes and six measured genes")
    stem = "Figure3B_shared_gene_patient_fingerprint"
    source = SOURCE_OUT / f"{stem}_source.csv"
    fingerprint_table.to_csv(source, index=False)

    # Plot the prespecified gene rows and patient columns, with a separate color for unmeasured values.
    expression_table = fingerprint_table.pivot(index="gene", columns="patient_label", values="gene_delta")
    expression_table = expression_table.reindex(index=genes, columns=paired_table["patient_label"].tolist())
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#DADADA")
    masked = np.ma.masked_invalid(expression_table.to_numpy(dtype=float))

    fig = new_figure(116, 92)
    ax = fig.add_axes([0.12, 0.28, 0.66, 0.58])
    image = ax.imshow(masked, cmap=cmap, vmin=-color_limit, vmax=color_limit, aspect="auto")
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=6.4)
    ax.set_xticks(
        np.arange(len(paired_table)),
        paired_table["patient_label"],
        rotation = 90,
        fontsize = 4.8,
    )
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(paired_table), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.text(
        len(paired_table) / 2 - 0.5,
        len(genes) - 1,
        "not measured",
        ha       = "center",
        va       = "center",
        fontsize = 5.2,
        color    = "#777777",
    )

    color_ax = fig.add_axes([0.19, 0.16, 0.49, 0.022])
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.7, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
    colorbar.outline.set_visible(False)

    # Draw the archived paired mean changes and confidence intervals next to the corresponding gene rows.
    summary_table = fingerprint_table.drop_duplicates("gene").set_index("gene").reindex(genes)
    sx = fig.add_axes([0.81, 0.28, 0.16, 0.58])
    sx.axvline(0, color="#A0A0A0", linewidth=0.65)
    for index, gene in enumerate(genes):
        row = summary_table.loc[gene]
        if not bool(row["measured"]):
            sx.text(0.52, index, "NA", fontsize=5.0, color="#888888", va="center", ha="center")
            continue
        significant = float(row["gene_fdr"]) < 0.05
        sx.hlines(index, float(row["ci_low"]), float(row["ci_high"]), color=recurrent, linewidth=1.0)
        sx.scatter(
            float(row["mean_gene_delta"]),
            index,
            s         = 22,
            facecolor = recurrent if significant else "white",
            edgecolor = recurrent,
            linewidth = 0.8,
            zorder    = 3,
        )
    sx.set_xlim(-0.45, 0.92)
    sx.set_ylim(len(genes) - 0.5, -0.5)
    sx.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.8)
    sx.set_yticks([])
    sx.set_xlabel("Mean Δ", fontsize=5.5, labelpad=2)
    sx.grid(axis="x", color=grid_color, linewidth=0.45)
    clean_axis(sx, keep_left=False, keep_bottom=True)

    fig.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.81,
        0.16,
        "Filled: gene-level\nFDR < 0.05",
        fontsize    = 5.1,
        color       = "#707070",
        va          = "top",
        linespacing = 1.15,
    )
    pdf, png = save_panel(fig, stem)
    panel_record = record("B", stem, "Patient-by-gene GeoMx delta fingerprint", source, pdf, png)
    return panel_record


def soften_tissue(rgb: np.ndarray) -> np.ndarray:
    """Return a muted RGB background using the original Step 53 conversion."""
    gray = rgb.mean(axis=2, keepdims=True)
    muted = 0.55 * rgb + 0.45 * gray
    background = np.clip(0.20 + 0.80 * muted, 0.0, 1.0)
    return background


def panel_c() -> dict[str, str]:
    """Join the frozen spot data and draw each score on the corresponding H&E image."""
    sections = [
        ("GBM030", "GSM8506649_03422", "fig3c_gbm030_spatial_raw20_mdsc_mes.csv"),
        ("GBM049", "GSM8506650_25526", "fig3d_gbm049_spatial_raw20_mdsc_mes.csv"),
    ]
    source_dir = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data"
    section_tables = []
    for sample, _, source_name in sections:
        section_table = pd.read_csv(source_dir / source_name)
        if set(section_table["sample"]) != {sample}:
            raise ValueError(f"Unexpected sample identifiers in the spot table for {sample}")
        section_tables.append(section_table)
    spatial_table_total = pd.concat(section_tables, ignore_index=True)
    if len(spatial_table_total) != 7629:
        raise ValueError("The frozen Visium data must contain 7,629 displayed spots")

    # For a given score, use the pooled 2nd and 98th percentiles as display limits for both rows.
    score_columns = [
        ("raw20", "MG-Inflam20", "MG-Inflam20 score"),
        ("mdsc_combined", "MDSC-like", "MDSC-like score"),
        ("mes", "Mesenchymal/hypoxia\ncomposite", "Mesenchymal/hypoxia composite score"),
    ]
    scales = {}
    for variable, _, _ in score_columns:
        low, high = spatial_table_total[variable].quantile([0.02, 0.98]).tolist()
        scales[variable] = (float(low), float(high))
        spatial_table_total[f"{variable}_scale_low"] = float(low)
        spatial_table_total[f"{variable}_scale_high"] = float(high)
    stem = "Figure3C_spatial_score_atlas"
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    spatial_table_total.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    # Keep the original 2-by-3 grid that also determines the Figure 3C scale calibration.
    fig = new_figure(142, 96)
    grid = fig.add_gridspec(
        2,
        3,
        left   = 0.10,
        right  = 0.98,
        bottom = 0.18,
        top    = 0.88,
        wspace = 0.05,
        hspace = 0.08,
    )
    cmap = mpl.colormaps["viridis"]
    axes: dict[tuple[int, int], mpl.axes.Axes] = {}
    image_dir = figure_tools.public_data_root(project_root) / (
        "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
    )
    for row_index, (sample, prefix, _) in enumerate(sections):
        section_table = spatial_table_total.loc[spatial_table_total["sample"].eq(sample)]
        with gzip.open(image_dir / f"{prefix}_tissue_hires_image.png.gz", "rb") as handle:
            image_bytes = io.BytesIO(handle.read())
        with Image.open(image_bytes) as image:
            rgb = np.asarray(image.convert("RGB"), dtype=np.float32) / 255.0
        tissue = soften_tissue(rgb)

        # Crop the displayed spot extent with 10% padding, bounded by the source image.
        x_min, x_max = section_table["x_hires"].min(), section_table["x_hires"].max()
        y_min, y_max = section_table["y_hires"].min(), section_table["y_hires"].max()
        x_margin = (x_max - x_min) * 0.10
        y_margin = (y_max - y_min) * 0.10
        x0 = max(0.0, float(x_min - x_margin))
        x1 = min(float(tissue.shape[1]), float(x_max + x_margin))
        y0 = max(0.0, float(y_min - y_margin))
        y1 = min(float(tissue.shape[0]), float(y_max + y_margin))

        # Plot the same registered spot coordinates in each column.
        for col_index, (variable, title, _) in enumerate(score_columns):
            ax = fig.add_subplot(grid[row_index, col_index])
            axes[(row_index, col_index)] = ax
            low, high = scales[variable]
            norm = mcolors.Normalize(vmin=low, vmax=high, clip=True)
            ax.imshow(tissue, origin="upper", interpolation="nearest")
            ax.scatter(
                section_table["x_hires"],
                section_table["y_hires"],
                c          = section_table[variable],
                cmap       = cmap,
                norm       = norm,
                s          = 8.5,
                alpha      = 0.94,
                linewidths = 0,
                rasterized = True,
                zorder     = 3,
            )
            ax.set_xlim(x0, x1)
            ax.set_ylim(y1, y0)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row_index == 0:
                ax.set_title(title, fontsize=7.4, fontweight="bold", pad=4)
            if col_index == 0:
                ax.text(
                    -0.03,
                    0.98,
                    f"{sample}\n{len(section_table):,} spots",
                    transform   = ax.transAxes,
                    ha          = "right",
                    va          = "top",
                    fontsize    = 6.0,
                    color       = text_color,
                    linespacing = 1.15,
                )

    # Use a common color bar for each score and align it with the lower map.
    for col_index, (variable, _, colorbar_label) in enumerate(score_columns):
        low, high = scales[variable]
        position = axes[(1, col_index)].get_position()
        color_ax = fig.add_axes([position.x0 + 0.02, 0.105, position.width - 0.04, 0.018])
        scalar = mpl.cm.ScalarMappable(norm=mcolors.Normalize(vmin=low, vmax=high), cmap=cmap)
        colorbar = fig.colorbar(scalar, cax=color_ax, orientation="horizontal")
        colorbar.set_ticks([low, (low + high) / 2, high])
        colorbar.ax.tick_params(labelsize=5.6, length=1.5, pad=1)
        colorbar.outline.set_visible(False)
        colorbar.set_label(colorbar_label, fontsize=5.8, labelpad=1.5)

    fig.text(0.10, 0.955, "Spatial score maps in two primary GBMs", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = record(
        "C",
        stem,
        "Two-slice maps of MG-Inflam20, MDSC-like, and MES-like spatial scores",
        source,
        pdf,
        png,
    )
    return panel_record


def panel_d() -> dict[str, str]:
    """Write and plot the frozen partial correlations for the two section and score combinations."""
    stem = "Figure3D_partial_rho_effect_plot"
    path = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3e_partial_rho_main.csv"
    original_table = pd.read_csv(path)
    correlation_table = original_table[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
    correlation_table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
    if len(correlation_table) != 4 or not (correlation_table["partial_rho"] > 0).all():
        raise ValueError("Figure 3D needs four positive effects for two targets and two sections")
    source = SOURCE_OUT / f"{stem}_source.csv"
    correlation_table.to_csv(source, index=False)

    # Mark each sample separately using the same within-score y offsets as the original plot.
    fig = new_figure(58, 88)
    ax = fig.add_axes([0.28, 0.20, 0.68, 0.66])
    target_y = {"MDSC": 1.0, "MES": 0.0}
    sample_offset = {"GBM030": 0.09, "GBM049": -0.09}
    sample_marker = {"GBM030": "o", "GBM049": "s"}
    target_color = {"MDSC": mdsc_color, "MES": mes_color}
    for target in ["MDSC", "MES"]:
        target_table = correlation_table.loc[correlation_table["target"].eq(target)].set_index("sample")
        for sample in ["GBM030", "GBM049"]:
            effect = float(target_table.loc[sample, "partial_rho"])
            row_y = target_y[target] + sample_offset[sample]
            ax.scatter(
                effect,
                row_y,
                s         = 34,
                marker    = sample_marker[sample],
                facecolor = target_color[target],
                edgecolor = "white",
                linewidth = 0.6,
                zorder    = 3,
                label     = sample,
            )
            ax.text(effect + 0.012, row_y, f"{effect:.3f}", fontsize=5.7, va="center", color=text_color)

    # Retain the zero reference, 0–0.50 axis and sample legend.
    ax.axvline(0, color="#999999", linewidth=0.7)
    ax.set_xlim(0, 0.50)
    ax.set_ylim(-0.48, 1.48)
    ax.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.4)
    ax.set_xlabel("Partial Spearman r")
    ax.grid(axis="x", color=grid_color, linewidth=0.5)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles = []
    for sample, marker in sample_marker.items():
        handles.append(
            mpl.lines.Line2D(
                [], [], marker=marker, linestyle="none", color=text_color, markerfacecolor="white", label=sample
            )
        )
    ax.legend(handles=handles, frameon=False, fontsize=5.4, loc="lower right", handletextpad=0.3)
    fig.text(0.15, 0.955, "Myeloid-adjusted niche association", fontsize=8.2, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = record("D", stem, "Descriptive partial-rho effects in two Visium tumors", source, pdf, png)
    return panel_record


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    """Build the original color and grayscale previews and save their panel rectangles."""
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 100, "top": 235, "right": 1660, "bottom": 1360},
        {"panel": "B", "left": 1740, "top": 235, "right": 3700, "bottom": 1360},
        {"panel": "C", "left": 100, "top": 1630, "right": 2770, "bottom": 2820},
        {"panel": "D", "left": 2870, "top": 1630, "right": 3700, "bottom": 2820},
    ]
    canvas = Image.new("RGB", (3800, 2920), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = figure_tools.image_font(58, bold=True)
    section_font = figure_tools.image_font(36, bold=True)

    # Keep the section labels and placement used in the original composite.
    drawer.text((100, 55), "Independent recurrence support", fill=recurrent, font=section_font)
    drawer.line((100, 118, 3700, 118), fill="#E8B2AA", width=4)
    drawer.text(
        (100, 1450),
        "Inflammatory-microglial program in MDSC/MES-rich niches",
        fill = mdsc_color,
        font = section_font,
    )
    drawer.line((100, 1515, 3700, 1515), fill="#A8D9CF", width=4)

    for item in layout:
        panel = item["panel"]
        figure_tools.paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 62), panel, fill="#111111", font=panel_font)

    # Keep the two 180 dpi image previews and the composite PDF.
    png = FIG_OUT / "Figure3_visual_variety_preview.png"
    gray = FIG_OUT / "Figure3_visual_variety_preview_grayscale.png"
    pdf = FIG_OUT / "Figure3_visual_variety_preview.pdf"
    source = SOURCE_OUT / "Figure3_visual_variety_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(180, 180), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(180, 180), optimize=True)
    canvas.save(pdf, format="PDF", resolution=180.0)
    layout_table = pd.DataFrame(layout)
    layout_table.to_csv(source, index=False)
    output_paths = (pdf, png, gray, source)
    return output_paths


def main() -> None:
    """Rebuild the four Step 53 panels, compose the previews, and write the manifest."""
    apply_publication_style()
    np.random.seed(20260714)

    # Confirm the downloaded image cache is available for this raw-panel entry point.
    image_dir = figure_tools.public_data_root(project_root) / (
        "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
    )
    if not image_dir.is_dir():
        raise FileNotFoundError(f"Missing GSE276841 spatial images at {image_dir}")
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    # Keep the original panel order and the shared GeoMx patient list for panels A and B.
    paired_table = patient_table()
    records = [panel_a(paired_table), panel_b(paired_table), panel_c(), panel_d()]
    if [row["panel"] for row in records] != list("ABCD"):
        raise ValueError("Figure 3 requires panels A through D in order")

    # Write one manifest row per panel, including paths to the shared composite.
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)
    manifest_table["preview_layout"] = str(preview_layout)
    manifest_path = WRITE_ROOT / "Figure3_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)

    # Write the archived Step 53 figure legend with the original wording and statistics.
    legend_path = WRITE_ROOT / "Figure3_legend.md"
    legend_text = (
        "# Figure 3 legend\n\n"
        "**Figure 3 | Independent myeloid-enriched tissue validation and spatial niche context of the "
        "MG-Inflam20 program.**\n\n"
        "**A,** Patient-level paired MG-Inflam20 scores in the independent Artzi et al. GeoMx cohort. "
        "Quality-controlled IBA1+ areas of illumination were averaged within each patient and time point before "
        "calculating the MG-Inflam20 score from the 18 measurable genes. Lines connect primary and recurrent specimens "
        "from the same patient; circle and square markers denote primary and recurrent samples, respectively. "
        "Sixteen of 22 paired IDH-wild-type patients were recurrent-higher; the mean paired change was +0.273 "
        "(95% CI, +0.094 to +0.452; BH-adjusted sign-flip FDR across three prespecified GeoMx entry definitions = "
        "0.00721). **B,** Patient-level recurrent-minus-primary expression changes for the seven shared leading-edge "
        "genes defined in Figure 1, displayed in the same patient order as panel A. Six genes were measured; CCL4 was "
        "not measured. Heatmap colors were clipped symmetrically at the 98th percentile of the absolute patient-gene "
        "changes; 3 of 132 measured cells exceeded the display limits. The right-hand summary shows mean paired "
        "changes and 95% CIs; filled symbols denote gene-level FDR < 0.05 after correction across all 18 measurable "
        "MG-Inflam20 genes. This panel is a gene-level fingerprint and does not replace the MG-Inflam20 score based "
        "on 18 measurable genes used in panel A. **C,** H&E-backed Visium maps from two independent untreated primary "
        "IDH-wild-type GBMs. Rows denote tumors and columns show the MG-Inflam20, combined MDSC-like, and MES-like "
        "scores. Nineteen of 20 MG-Inflam20 genes were measurable; AC253572.2 was not measured. For each score, the "
        "color scale is shared across the two tumors and clipped at the pooled 2nd and 98th percentiles for display. "
        "These maps describe multicellular neighborhood context and do not assign MDSC identity to MG-Inflam20-high "
        "spots. **D,** Partial Spearman correlations between the MG-Inflam20 score and combined MDSC-like or MES-like "
        "scores after within-tumor rank residualization for total myeloid score. The two tumors are shown as "
        "independent points and are not connected. Only effect sizes are shown because the analysis contains two "
        "biological tumors and the spot-level tests did not correct for spatial autocorrelation.\n\n"
        "Panel A provides independent tissue-level recurrence support with patient as the statistical unit; panel B "
        "decomposes this effect across prespecified leading-edge genes. Figure 2 localized the MG-Inflam20 program "
        "preferentially to MCG1/MCG2-like inflammatory-microglial states rather than MDSC-like transcriptional states. "
        "Here, panels C-D show that MG-Inflam20-high Visium spots have positive myeloid-adjusted associations with "
        "MDSC/MES-like spatial scores. Cell-state identity and multicellular spatial niche are therefore distinct "
        "analytical levels: an inflammatory-microglial program can occur within an MDSC/MES-rich neighborhood. "
        "Panels C-D are cross-sectional and do not test recurrence. Treatment-stratified GeoMx analyses, the complete "
        "18-gene forest, and E-/M-MDSC sensitivity analyses are retained as supplementary analyses.\n"
    )
    legend_path.write_text(legend_text, encoding="utf-8")

    # Retain the existing build summary with the panel outputs and source-table manifest.
    result_path = WRITE_ROOT / "FINAL_RESULT.md"
    result_text = (
        "# Figure 3 visual-variety rebuild\n\n"
        f"- A. Paired GeoMx scores for {len(paired_table)} patients.\n"
        "- B. Patient-by-shared-leading-edge gene changes in the original order.\n"
        "- C. H&E-backed spatial maps of inflammatory, MDSC-like, and MES-like scores.\n"
        "- D. Descriptive myeloid-adjusted spatial effects.\n"
        f"- Manifest `{manifest_path}`\n"
        f"- Legend `{legend_path}`\n"
    )
    result_path.write_text(result_text, encoding="utf-8")
    for path in (preview_pdf, preview_png, preview_gray, preview_layout, manifest_path, legend_path, result_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 3 output: {path.name}")
    print(f"Figure 3 panels: {len(records)}")


if __name__ == "__main__":
    main()
