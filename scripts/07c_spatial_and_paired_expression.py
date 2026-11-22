#!/usr/bin/env python3
# Purpose: Draw the Figure 3 panels from the strict GeoMx pairs and descriptive Visium data.
# Inputs: Step 41 Figure3/source_data patient, gene, and spatial tables and the Step 51 shared leading-edge gene list.
# Under MILLER_IM_DATA_ROOT, use GeoMx metadata.csv and normalized.txt and the two compressed Visium H&E images.
# Outputs: source_data, Figure3_panel_manifest.csv, Figure3_legend.md, and FINAL_RESULT.md under
# write/53_figure3_visual_variety_rebuild/Figure3, with panels A–D and Figure3_visual_variety_preview under figures.
# Run python scripts/07c_spatial_and_paired_expression.py from the repository root.
# The Python workflow uses pandas, NumPy, Matplotlib, and Pillow for the released calculations and figures.
# Step 55 imports the patient table and panels A and C. The main workflow runs when invoked as a script.

# Each panel stem below produces a PDF, PNG, and source_data/<stem>_source.csv.
#   Figure3A_geomx_paired_dumbbell
#   Figure3B_shared_gene_patient_fingerprint
#   Figure3C_spatial_score_atlas
#   Figure3D_partial_rho_effect_plot
# Panel C uses .csv.gz. Preview filenames use the _grayscale.png and _layout.csv suffixes for those assets.

from __future__ import annotations

import gzip
import io
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
image_font = figure_style.image_font
new_figure = figure_style.new_figure
panel_record = figure_style.panel_record
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel

# A direct run initializes its directories and styles. Step 55 supplies its own when importing this file.
if __name__ == "__main__":
    apply_publication_style()
    np.random.seed(20260714)
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write/53_figure3_visual_variety_rebuild/Figure3"
    source_out = write_root / "source_data"
    source_out.mkdir(parents=True, exist_ok=True)
    figure_out = root / "figures/53_figure3_visual_variety_rebuild/Figure3"
    figure_out.mkdir(parents=True, exist_ok=True)

    # Check the Visium image directory from the public data setting before drawing.
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    spatial_image_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
    if not spatial_image_root.exists():
        raise FileNotFoundError(f"Visium tissue image directory is unavailable: {spatial_image_root}")


# Both Figure 3 builds use one row per strict GeoMx patient pair, in the same order.
def patient_table() -> pd.DataFrame:
    """Return the ordered GeoMx table and require the frozen 22 pairs and 16 positive changes."""
    root = Path(__file__).resolve().parents[1]
    patient_path = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3a_geomx22_paired_change.csv"
    patients_table = pd.read_csv(patient_path)
    if len(patients_table) != 22 or int((patients_table["delta"] > 0).sum()) != 16:
        raise ValueError("The strict GeoMx cohort requires 22 pairs with 16 recurrent-higher scores")

    # Keep the original patient ordering, labels, direction, and eligibility.
    patients_table = patients_table.sort_values(["delta", "Patient_No"], ascending=[False, True]).reset_index(drop=True)
    patients_table["patient_label"] = patients_table["Patient_No"].map(lambda value: f"P{int(value):02d}")
    patients_table["patient_order"] = np.arange(len(patients_table))
    patients_table["direction"] = np.where(patients_table["delta"] > 0, "Recurrent higher", "Primary higher")
    return patients_table


if __name__ == "__main__":
    patients_table = patient_table()


# Both Figure 3 builds draw panel A using the same paired patient table.
def panel_a(patients_table: pd.DataFrame, source_out: Path, figure_out: Path, creator: str) -> dict[str, str]:
    """Write panel A for the ordered patient data and return its manifest record."""
    stem = "Figure3A_geomx_paired_dumbbell"

    # Write the panel source table with the existing columns in their original order.
    source_path = source_out / f"{stem}_source.csv"
    patients_table.to_csv(source_path, index=False)

    # Use the original figure size in millimeters and the same PDF axes positions.
    figure = new_figure(92, 92)
    axis = figure.add_axes([0.22, 0.20, 0.72, 0.65])
    patient_positions = np.arange(len(patients_table), dtype=float)
    for y_position, row in zip(patient_positions, patients_table.itertuples()):
        axis.plot(
            [row.Primary, row.Recurrence],
            [y_position, y_position],
            color     = "#A6A6A6",
            linewidth = 1.15,
            alpha     = 0.72,
            zorder    = 1,
        )
    axis.scatter(
        patients_table["Primary"],
        patient_positions,
        s         = 22,
        marker    = "o",
        color     = "#3C5488",
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Primary",
    )
    axis.scatter(
        patients_table["Recurrence"],
        patient_positions,
        s         = 22,
        marker    = "s",
        color     = "#E64B35",
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Recurrence",
    )

    # Draw the pooled mean row separately using the original patient scores.
    mean_y = len(patients_table) + 0.8
    primary_mean = float(patients_table["Primary"].mean())
    recurrent_mean = float(patients_table["Recurrence"].mean())
    axis.plot(
        [primary_mean, recurrent_mean],
        [mean_y, mean_y],
        color     = "#222222",
        linewidth = 1.5,
        zorder    = 2,
    )
    axis.scatter(primary_mean, mean_y, s=34, marker="o", color="#3C5488", edgecolor="white", linewidth=0.6, zorder=4)
    axis.scatter(recurrent_mean, mean_y, s=34, marker="s", color="#E64B35", edgecolor="white", linewidth=0.6, zorder=4)

    axis.set_yticks(
        list(patient_positions) + [mean_y],
        patients_table["patient_label"].tolist() + ["Mean"],
        fontsize = 5.6,
    )

    # Use the original plot limits and data order.
    axis.set_ylim(mean_y + 1.0, -1.0)
    values = patients_table[["Primary", "Recurrence"]].to_numpy(dtype=float)
    margin = max(0.12, float(np.ptp(values)) * 0.08)

    # Keep the horizontal scale and labels for the original comparison.
    axis.set_xlim(float(np.nanmin(values) - margin), float(np.nanmax(values) + margin))
    axis.set_xlabel("Patient-level Miller-IM score")
    axis.grid(axis="x", color="#E7E7E7", linewidth=0.5)
    clean_axis(axis, keep_left=False, keep_bottom=True)
    axis.tick_params(axis="y", length=0, pad=3)
    handles, labels = axis.get_legend_handles_labels()
    figure.legend(
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

    figure.text(0.12, 0.955, "IBA1+ tissue recurrence support", fontsize=8.5, fontweight="bold", va="top")

    # Save the panel as PDF and render its 600 dpi PNG from that file.
    pdf_path, png_path = save_panel(figure, figure_out, stem, creator)
    panel_result = panel_record(
        "A", stem, "Patient-level paired GeoMx Miller-IM scores", source_path, pdf_path, png_path
    )
    return panel_result


if __name__ == "__main__":
    records = [panel_a(patients_table, source_out, figure_out, "Step53 recurrent GBM Figure 3")]


# Both Figure 3 builds read the shared leading-edge genes in the original Figure 1 order.
def shared_genes() -> list[str]:
    """Return the Figure 1 shared genes in order and require the frozen seven-gene definition."""
    root = Path(__file__).resolve().parents[1]
    definition_path = (
        root / "write/51_figure1_final_reorganized/Figure1/source_data/Figure1_raw20_definition_source.csv"
    )
    definition_table = pd.read_csv(definition_path)
    genes = (
        definition_table.loc[definition_table["shared_leading_edge"].astype(bool)]
        .sort_values("raw20_order")["gene"]
        .drop_duplicates()
        .tolist()
    )
    expected = ["PDK4", "SGK1", "CH25H", "SIGLEC8", "KLF6", "FOLR2", "CCL4"]
    if genes != expected:
        raise ValueError(f"Unexpected Figure 1 shared leading-edge genes: {genes}")
    return genes


if __name__ == "__main__":
    genes = shared_genes()


# Use the original ROI eligibility and log2(x + 1) values to compute paired changes for both builds.
def strict_gene_delta(patients_table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Return paired patient-by-gene changes and require the original coverage and complete pairs."""
    root = Path(__file__).resolve().parents[1]
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    geomx_root = external_root / "write/26_pure_bioinformatics_dataset_rescue/source_metadata/GeoMx_zenodo16839828"
    metadata_path = geomx_root / "metadata.csv"
    expression_path = geomx_root / "normalized.txt"
    if not metadata_path.exists() or not expression_path.exists():
        raise FileNotFoundError("The normalized GeoMx expression matrix is unavailable")
    metadata_table = pd.read_csv(metadata_path)
    eligible_roi_table = metadata_table.loc[
        metadata_table["iba1"].eq(True)
        & metadata_table["QCFlags"].isna()
        & metadata_table["IDH_status"].eq("IDH_WT")
        & metadata_table["tumor_setting"].isin(["Primary", "Recurrence"])
    ].copy()

    # The normalized matrix has one row per gene and one column per ROI. Apply log2(x + 1).
    expression_table = pd.read_csv(expression_path, sep="\t", index_col=0)
    measured = [gene for gene in genes if gene in expression_table.index]
    missing = [gene for gene in genes if gene not in expression_table.index]
    if measured != genes[:-1] or missing != ["CCL4"]:
        raise ValueError(f"Unexpected GeoMx shared leading-edge coverage: measured={measured}, missing={missing}")
    ids = eligible_roi_table["Unnamed: 0"].tolist()
    log_expression_table = np.log2(expression_table.loc[measured, ids] + 1.0)

    # Align the area metadata to the expression columns, then average within patient and time point.
    roi_index_table = eligible_roi_table.set_index("Unnamed: 0")
    roi_patient_ids = roi_index_table.loc[log_expression_table.columns, "Patient_No"].to_numpy()
    roi_settings = roi_index_table.loc[log_expression_table.columns, "tumor_setting"].to_numpy()
    patient_gene_table = log_expression_table.T.assign(Patient_No=roi_patient_ids, tumor_setting=roi_settings)
    patient_gene_table = patient_gene_table.groupby(["Patient_No", "tumor_setting"], observed=True)[measured].mean()
    patient_gene_table = patient_gene_table.reset_index()

    # Keep one row per patient and gene. Unmeasured CCL4 retains missing values.
    patient_ids = patients_table["Patient_No"].astype(int).tolist()
    rows: list[dict[str, object]] = []
    for gene_order, gene in enumerate(genes):
        if gene not in measured:
            for patient_order, patient_id in enumerate(patient_ids):
                rows.append(
                    {
                        "Patient_No": patient_id,
                        "patient_label": f"P{patient_id:02d}",
                        "patient_order": patient_order,
                        "gene": gene,
                        "gene_order": gene_order,
                        "measured": False,
                        "gene_delta": np.nan,
                    }
                )
            continue
        paired_gene_table = patient_gene_table.pivot(index="Patient_No", columns="tumor_setting", values=gene)
        for patient_order, patient_id in enumerate(patient_ids):
            if (
                patient_id not in paired_gene_table.index
                or pd.isna(paired_gene_table.loc[patient_id, ["Primary", "Recurrence"]]).any()
            ):
                raise ValueError(f"Missing paired {gene} values for strict GeoMx patient {patient_id}")
            rows.append(
                {
                    "Patient_No": patient_id,
                    "patient_label": f"P{patient_id:02d}",
                    "patient_order": patient_order,
                    "gene": gene,
                    "gene_order": gene_order,
                    "measured": True,
                    "gene_delta": float(
                        paired_gene_table.loc[patient_id, "Recurrence"] - paired_gene_table.loc[patient_id, "Primary"]
                    ),
                }
            )
    gene_delta_table = pd.DataFrame(rows)
    delta_by_patient = patients_table.set_index("Patient_No")["delta"].to_dict()
    gene_delta_table["raw20_delta"] = gene_delta_table["Patient_No"].map(delta_by_patient)

    # Annotate from the original 18-gene table of effects, P values, and intervals.
    summary_path = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3b_geomx18_gene_forest.csv"
    gene_summary_table = pd.read_csv(summary_path)
    gene_summary_table = gene_summary_table.set_index("gene")
    gene_delta_table["mean_gene_delta"] = gene_delta_table["gene"].map(gene_summary_table["mean_log2_delta"])
    gene_delta_table["ci_low"] = gene_delta_table["gene"].map(gene_summary_table["ci_low"])
    gene_delta_table["ci_high"] = gene_delta_table["gene"].map(gene_summary_table["ci_high"])
    gene_delta_table["gene_fdr"] = gene_delta_table["gene"].map(gene_summary_table["fdr"])
    gene_delta_table.loc[~gene_delta_table["measured"], ["mean_gene_delta", "ci_low", "ci_high", "gene_fdr"]] = np.nan

    # Use the 98th percentile of absolute changes to clip the heatmap color range.
    measured_values = gene_delta_table.loc[gene_delta_table["measured"], "gene_delta"].abs().to_numpy(dtype=float)
    color_limit = float(np.quantile(measured_values, 0.98))
    gene_delta_table["color_limit"] = color_limit
    return gene_delta_table


if __name__ == "__main__":
    gene_delta_table = strict_gene_delta(patients_table, genes)


# Draw panel B with all seven genes, gray for the unmeasured gene, and the original 98th-percentile color limit.
if __name__ == "__main__":
    if len(gene_delta_table) != 154 or int(gene_delta_table["measured"].sum()) != 132:
        raise ValueError("Figure 3B requires 22 patients and seven genes, of which six are measurable")
    stem = "Figure3B_shared_gene_patient_fingerprint"

    # Write the panel source table with the existing columns in their original order.
    source_path = source_out / f"{stem}_source.csv"
    gene_delta_table.to_csv(source_path, index=False)

    delta_pivot_table = gene_delta_table.pivot(index="gene", columns="patient_label", values="gene_delta").reindex(
        index=genes, columns=patients_table["patient_label"].tolist()
    )
    color_limit = float(gene_delta_table["color_limit"].iloc[0])
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#DADADA")
    masked = np.ma.masked_invalid(delta_pivot_table.to_numpy(dtype=float))

    # Use the original figure size in millimeters and the same PDF axes positions.
    figure = new_figure(116, 92)
    axis = figure.add_axes([0.12, 0.28, 0.66, 0.58])
    image = axis.imshow(masked, cmap=cmap, vmin=-color_limit, vmax=color_limit, aspect="auto")
    axis.set_yticks(np.arange(len(genes)), genes, fontsize=6.4)
    axis.set_xticks(
        np.arange(len(patients_table)),
        patients_table["patient_label"],
        rotation = 90,
        fontsize = 4.8,
    )
    axis.tick_params(length=0, pad=2)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_xticks(np.arange(-0.5, len(patients_table), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.35)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.text(
        len(patients_table) / 2 - 0.5,
        len(genes) - 1,
        "not measured",
        ha       = "center",
        va       = "center",
        fontsize = 5.2,
        color    = "#777777",
    )

    # The color bar shows the display scale for the gene changes.
    color_axis = figure.add_axes([0.19, 0.16, 0.49, 0.022])
    colorbar = figure.colorbar(image, cax=color_axis, orientation="horizontal")
    colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.7, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
    colorbar.outline.set_visible(False)

    # Use the original effect estimates, confidence intervals, and FDRs for the side plot.
    gene_plot_summary_table = gene_delta_table.drop_duplicates("gene").set_index("gene").reindex(genes)
    summary_axis = figure.add_axes([0.81, 0.28, 0.16, 0.58])
    summary_axis.axvline(0, color="#A0A0A0", linewidth=0.65)
    for index, gene in enumerate(genes):
        row = gene_plot_summary_table.loc[gene]
        if not bool(row["measured"]):
            summary_axis.text(0.52, index, "NA", fontsize=5.0, color="#888888", va="center", ha="center")
            continue
        significant = float(row["gene_fdr"]) < 0.05
        summary_axis.hlines(index, float(row["ci_low"]), float(row["ci_high"]), color="#E64B35", linewidth=1.0)
        summary_axis.scatter(
            float(row["mean_gene_delta"]),
            index,
            s         = 22,
            facecolor = "#E64B35" if significant else "white",
            edgecolor = "#E64B35",
            linewidth = 0.8,
            zorder    = 3,
        )

    # Keep the zero reference, axis limits, and marker definitions.
    summary_axis.set_xlim(-0.45, 0.92)
    summary_axis.set_ylim(len(genes) - 0.5, -0.5)
    summary_axis.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.8)
    summary_axis.set_yticks([])
    summary_axis.set_xlabel("Mean Δ", fontsize=5.5, labelpad=2)
    summary_axis.grid(axis="x", color="#E7E7E7", linewidth=0.45)
    clean_axis(summary_axis, keep_left=False, keep_bottom=True)

    figure.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=8.5, fontweight="bold", va="top")
    figure.text(
        0.81,
        0.16,
        "Filled: gene-level\nFDR < 0.05",
        fontsize    = 5.1,
        color       = "#707070",
        va          = "top",
        linespacing = 1.15,
    )

    # Save the panel as PDF and render its 600 dpi PNG from that file.
    pdf_path, png_path = save_panel(figure, figure_out, stem, "Step53 recurrent GBM Figure 3")
    records.append(panel_record("B", stem, "Patient-by-gene GeoMx delta fingerprint", source_path, pdf_path, png_path))


# Both builds share panel C. The grayscale setting selects the original Step 53 or Step 55 conversion.
def panel_c(source_out: Path, figure_out: Path, creator: str, grayscale: bool) -> dict[str, str]:
    """Write panel C and return its manifest record for the requested directory, PDF metadata, and tissue rendering
    mode."""
    root = Path(__file__).resolve().parents[1]
    step41_source = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data"
    gbm030_table = pd.read_csv(step41_source / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv")
    gbm049_table = pd.read_csv(step41_source / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv")
    spatial_tables = {"GBM030": gbm030_table, "GBM049": gbm049_table}
    for sample, spatial_table in spatial_tables.items():
        if set(spatial_table["sample"]) != {sample}:
            raise ValueError(f"Unexpected sample name in the Visium spot table for {sample}")
    spatial_table_total = pd.concat(spatial_tables.values(), ignore_index=True)
    if len(spatial_table_total) != 7629:
        raise ValueError("The Visium table must contain 7,629 displayed spots")

    # Use the 2nd and 98th percentiles of the pooled spot data to set a shared color scale per score.
    for variable in ["raw20", "mdsc_combined", "mes"]:
        low, high = spatial_table_total[variable].quantile([0.02, 0.98]).tolist()
        spatial_table_total[f"{variable}_scale_low"] = float(low)
        spatial_table_total[f"{variable}_scale_high"] = float(high)
    spatial_tables = {}
    for sample in ["GBM030", "GBM049"]:
        spatial_tables[sample] = spatial_table_total.loc[spatial_table_total["sample"].eq(sample)].copy()
    stem = "Figure3C_spatial_score_atlas"

    # Write the panel source table with the existing columns in their original order.
    source_path = source_out / f"{stem}_source.csv.gz"
    spatial_table_total.to_csv(source_path, index=False, compression={"method": "gzip", "mtime": 0})

    score_columns = [
        ("raw20", "Miller-IM", "Miller-IM score"),
        ("mdsc_combined", "MDSC-like", "MDSC-like score"),
        ("mes", "MES-like", "MES-like score"),
    ]

    # Read each score's pooled lower and upper color limits for use in both sections.
    scales = {
        variable: (
            float(spatial_table_total[f"{variable}_scale_low"].iloc[0]),
            float(spatial_table_total[f"{variable}_scale_high"].iloc[0]),
        )
        for variable, _, _ in score_columns
    }

    # Use the original figure size in millimeters and the same PDF axes positions.
    figure = new_figure(142, 96)
    grid = figure.add_gridspec(
        2,
        3,
        left   = 0.10,
        right  = 0.98,
        bottom = 0.18,
        top    = 0.88,
        wspace = 0.05,
        hspace = 0.08,
    )

    # Use the appropriate archived tissue display conversion for each build.
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    spatial_image_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
    spatial_prefix = {"GBM030": "GSM8506649_03422", "GBM049": "GSM8506650_25526"}
    images = {}
    for sample in spatial_tables:
        prefix = spatial_prefix[sample]
        image_path = spatial_image_root / f"{prefix}_tissue_hires_image.png.gz"
        with gzip.open(image_path, "rb") as handle:
            tissue_image = Image.open(io.BytesIO(handle.read())).convert("RGB")
        rgb = np.asarray(tissue_image, dtype=np.float32) / 255.0
        if grayscale:
            luminance = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
            display = 0.68 + 0.32 * np.clip(luminance, 0.0, 1.0)
            images[sample] = np.repeat(display[:, :, None], 3, axis=2)
        else:
            gray_channel = rgb.mean(axis=2, keepdims=True)
            muted = 0.55 * rgb + 0.45 * gray_channel
            images[sample] = np.clip(0.20 + 0.80 * muted, 0.0, 1.0)

    # Prepare the tissue backgrounds before drawing the sections and scores in their original order.
    cmap = mpl.colormaps["viridis"]
    axes: dict[tuple[int, int], mpl.axes.Axes] = {}
    for row_index, sample in enumerate(["GBM030", "GBM049"]):
        spatial_table = spatial_tables[sample]
        tissue = images[sample]

        # Crop to the spot extent plus its original 10% padding.
        x_min, x_max = spatial_table["x_hires"].min(), spatial_table["x_hires"].max()
        y_min, y_max = spatial_table["y_hires"].min(), spatial_table["y_hires"].max()
        x_margin = (x_max - x_min) * 0.10
        y_margin = (y_max - y_min) * 0.10
        x0 = max(0.0, float(x_min - x_margin))
        x1 = min(float(tissue.shape[1]), float(x_max + x_margin))
        y0 = max(0.0, float(y_min - y_margin))
        y1 = min(float(tissue.shape[0]), float(y_max + y_margin))
        for col_index, (variable, title, _) in enumerate(score_columns):
            axis = figure.add_subplot(grid[row_index, col_index])
            axes[(row_index, col_index)] = axis
            low, high = scales[variable]
            norm = mcolors.Normalize(vmin=low, vmax=high, clip=True)
            axis.imshow(tissue, origin="upper", interpolation="nearest")
            axis.scatter(
                spatial_table["x_hires"],
                spatial_table["y_hires"],
                c          = spatial_table[variable],
                cmap       = cmap,
                norm       = norm,
                s          = 8.5,
                alpha      = 0.94,
                linewidths = 0,
                rasterized = True,
                zorder     = 3,
            )

            # Keep the spatial coordinates and axis labels in their original orientation.
            axis.set_xlim(x0, x1)

            # Set the map bounds to the original tissue crop.
            axis.set_ylim(y1, y0)

            # Hide the spatial ticks and retain the tissue extent and orientation.
            axis.set_xticks([])
            axis.set_yticks([])
            for spine in axis.spines.values():
                spine.set_visible(False)
            if row_index == 0:
                axis.set_title(title, fontsize=7.4, fontweight="bold", pad=4)
            if col_index == 0:
                axis.text(
                    -0.03,
                    0.98,
                    f"{sample}\n{len(spatial_table):,} spots",
                    transform   = axis.transAxes,
                    ha          = "right",
                    va          = "top",
                    fontsize    = 6.0,
                    color       = "#222222",
                    linespacing = 1.15,
                )

    for col_index, (variable, _, colorbar_label) in enumerate(score_columns):
        low, high = scales[variable]
        position = axes[(1, col_index)].get_position()

        # The color bar shows the display scale for the spatial scores.
        color_axis = figure.add_axes([position.x0 + 0.02, 0.105, position.width - 0.04, 0.018])
        scalar = mpl.cm.ScalarMappable(norm=mcolors.Normalize(vmin=low, vmax=high), cmap=cmap)
        colorbar = figure.colorbar(scalar, cax=color_axis, orientation="horizontal")
        colorbar.set_ticks([low, (low + high) / 2, high])
        colorbar.ax.tick_params(labelsize=5.6, length=1.5, pad=1)
        colorbar.outline.set_visible(False)
        colorbar.set_label(colorbar_label, fontsize=5.8, labelpad=1.5)

    figure.text(0.10, 0.955, "Spatial score maps in two primary GBMs", fontsize=8.5, fontweight="bold", va="top")

    # Save the panel as PDF and render its 600 dpi PNG from that file.
    pdf_path, png_path = save_panel(figure, figure_out, stem, creator)
    panel_result = panel_record(
        "C",
        stem,
        "Two-slice maps of Miller-IM, MDSC-like, and MES-like spatial scores",
        source_path,
        pdf_path,
        png_path,
    )
    return panel_result


if __name__ == "__main__":
    records.append(panel_c(source_out, figure_out, "Step53 recurrent GBM Figure 3", False))


# Draw panel D from the four archived partial correlations.
if __name__ == "__main__":
    partial_path = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3e_partial_rho_main.csv"
    partial_source_table = pd.read_csv(partial_path)
    partial_table = partial_source_table[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
    partial_table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
    if len(partial_table) != 4 or not (partial_table["partial_rho"] > 0).all():
        raise ValueError("Figure 3D requires four positive effects for the MDSC and MES scores in two sections")
    stem = "Figure3D_partial_rho_effect_plot"

    # Write the panel source table with the existing columns in their original order.
    source_path = source_out / f"{stem}_source.csv"
    partial_table.to_csv(source_path, index=False)

    # Use the original figure size in millimeters and the same PDF axes positions.
    figure = new_figure(58, 88)
    axis = figure.add_axes([0.28, 0.20, 0.68, 0.66])
    target_y = {"MDSC": 1.0, "MES": 0.0}
    sample_offset = {"GBM030": 0.09, "GBM049": -0.09}
    sample_marker = {"GBM030": "o", "GBM049": "s"}
    target_color = {"MDSC": "#159D88", "MES": "#7E6148"}
    for target in ["MDSC", "MES"]:
        target_table = partial_table.loc[partial_table["target"].eq(target)].set_index("sample")
        x_values = [float(target_table.loc[sample, "partial_rho"]) for sample in ["GBM030", "GBM049"]]
        y_values = [target_y[target] + sample_offset[sample] for sample in ["GBM030", "GBM049"]]
        for sample, x_position, y_position in zip(["GBM030", "GBM049"], x_values, y_values):
            axis.scatter(
                x_position,
                y_position,
                s         = 34,
                marker    = sample_marker[sample],
                facecolor = target_color[target],
                edgecolor = "white",
                linewidth = 0.6,
                zorder    = 3,
                label     = sample,
            )
            axis.text(x_position + 0.012, y_position, f"{x_position:.3f}", fontsize=5.7, va="center", color="#222222")
    axis.axvline(0, color="#999999", linewidth=0.7)

    # Keep the horizontal scale and labels for the original comparison.
    axis.set_xlim(0, 0.50)

    # Use the original plot limits and data order.
    axis.set_ylim(-0.48, 1.48)
    axis.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.4)
    axis.set_xlabel("Partial Spearman r")
    axis.grid(axis="x", color="#E7E7E7", linewidth=0.5)
    clean_axis(axis, keep_left=False, keep_bottom=True)
    axis.tick_params(axis="y", length=0, pad=3)

    # Keep the original colors, order, and symbols in the section legend.
    handles = [
        mpl.lines.Line2D(
            [], [], marker="o", linestyle="none", color="#222222", markerfacecolor="white", label="GBM030"
        ),
        mpl.lines.Line2D(
            [], [], marker="s", linestyle="none", color="#222222", markerfacecolor="white", label="GBM049"
        ),
    ]
    axis.legend(handles=handles, frameon=False, fontsize=5.4, loc="lower right", handletextpad=0.3)
    figure.text(0.15, 0.955, "Myeloid-adjusted niche association", fontsize=8.2, fontweight="bold", va="top")

    # Save the panel as PDF and render its 600 dpi PNG from that file.
    pdf_path, png_path = save_panel(figure, figure_out, stem, "Step53 recurrent GBM Figure 3")
    records.append(
        panel_record("D", stem, "Descriptive partial-rho effects in two Visium tumors", source_path, pdf_path, png_path)
    )


# Keep the Pillow composite layout and 180 dpi resolution, using the shared font and placement helpers.
if __name__ == "__main__":
    canvas = Image.new("RGB", (3800, 2920), "white")
    drawer = ImageDraw.Draw(canvas)
    section_font = image_font(36, bold=True)

    drawer.text((100, 55), "Independent recurrence support", fill="#E64B35", font=section_font)
    drawer.line((100, 118, 3700, 118), fill="#E8B2AA", width=4)
    drawer.text(
        (100, 1450),
        "Inflammatory-microglial program in MDSC/MES-rich niches",
        fill = "#159D88",
        font = section_font,
    )
    drawer.line((100, 1515, 3700, 1515), fill="#A8D9CF", width=4)

    if [row["panel"] for row in records] != list("ABCD"):
        raise ValueError("Figure 3 requires panels A through D in order")

    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 100, "top": 235, "right": 1660, "bottom": 1360},
        {"panel": "B", "left": 1740, "top": 235, "right": 3700, "bottom": 1360},
        {"panel": "C", "left": 100, "top": 1630, "right": 2770, "bottom": 2820},
        {"panel": "D", "left": 2870, "top": 1630, "right": 3700, "bottom": 2820},
    ]
    panel_font = image_font(58, bold=True)
    for item in layout:
        panel = item["panel"]
        paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 62), panel, fill="#111111", font=panel_font)

    # Save the Pillow preview with the original pixel size and resolution for layout review.
    png_path = figure_out / "Figure3_visual_variety_preview.png"
    canvas.save(png_path, format="PNG", dpi=(180, 180), optimize=True)
    preview_gray_path = figure_out / "Figure3_visual_variety_preview_grayscale.png"
    ImageOps.grayscale(canvas).convert("RGB").save(preview_gray_path, format="PNG", dpi=(180, 180), optimize=True)
    pdf_path = figure_out / "Figure3_visual_variety_preview.pdf"
    canvas.save(pdf_path, format="PDF", resolution=180.0)
    layout_table = pd.DataFrame(layout)
    source_path = source_out / "Figure3_visual_variety_preview_layout.csv"
    layout_table.to_csv(source_path, index=False)


if __name__ == "__main__":

    # Write one row per panel with the existing manifest fields, digests, and preview references.
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(pdf_path)
    manifest_table["preview_png"] = str(png_path)
    manifest_table["preview_grayscale"] = str(preview_gray_path)
    manifest_table["preview_layout"] = str(source_path)
    manifest_path = write_root / "Figure3_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)

    # Keep the archived figure legend and standalone summary text.
    content = (
        "# Figure 3 legend\n\n**Figure 3 | Independent myeloid-enriched tissue validation and spatial "
        "niche context of the Miller-IM program.**\n\n**A,** Patient-level paired Miller-IM scores in the "
        "independent Artzi et al. GeoMx cohort. Quality-controlled IBA1+ areas of illumination were "
        "averaged within each patient and time point before calculating the Miller-IM score from the 18 "
        "measurable genes. Lines connect primary and recurrent specimens from the same patient; circle "
        "and square markers denote primary and recurrent samples, respectively. Sixteen of 22 paired "
        "IDH-wild-type patients were recurrent-higher; the mean paired change was +0.273 (95% CI, +0.094"
        " to +0.452; BH-adjusted sign-flip FDR across three prespecified GeoMx entry definitions = "
        "0.00721). **B,** Patient-level recurrent-minus-primary expression changes for the seven shared "
        "leading-edge genes defined in Figure 1, displayed in the same patient order as panel A. Six "
        "genes were measured; CCL4 was not measured. Heatmap colors were clipped symmetrically at the "
        "98th percentile of the absolute patient-gene changes; 3 of 132 measured cells exceeded the "
        "display limits. The right-hand summary shows mean paired changes and 95% CIs; filled symbols "
        "denote gene-level FDR < 0.05 after correction across all 18 measurable Miller-IM genes. This "
        "panel is a gene-level fingerprint and does not replace the Miller-IM score based on 18 "
        "measurable genes used in panel A. **C,** H&E-backed Visium maps from two independent untreated "
        "primary IDH-wild-type GBMs. Rows denote tumors and columns show the Miller-IM, combined "
        "MDSC-like, and MES-like scores. Nineteen of 20 Miller-IM genes were measurable; AC253572.2 was "
        "not measured. For each score, the color scale is shared across the two tumors and clipped at "
        "the pooled 2nd and 98th percentiles for display. These maps describe multicellular neighborhood"
        " context and do not assign MDSC identity to Miller-IM-high spots. **D,** Partial Spearman "
        "correlations between the Miller-IM score and combined MDSC-like or MES-like scores after "
        "within-tumor rank residualization for total myeloid score. The two tumors are shown as "
        "independent points and are not connected. Only effect sizes are shown because the analysis "
        "contains two biological tumors and the spot-level tests did not correct for spatial "
        "autocorrelation.\n\nPanel A provides independent tissue-level recurrence support with patient as "
        "the statistical unit; panel B decomposes this effect across prespecified leading-edge genes. "
        "Figure 2 localized the Miller-IM program preferentially to MCG1/MCG2-like "
        "inflammatory-microglial states rather than MDSC-like transcriptional states. Here, panels C-D "
        "show that Miller-IM-high Visium spots have positive myeloid-adjusted associations with "
        "MDSC/MES-like spatial scores. Cell-state identity and multicellular spatial niche are therefore"
        " distinct analytical levels: an inflammatory-microglial program can occur within an "
        "MDSC/MES-rich neighborhood. Panels C-D are cross-sectional and do not test recurrence. "
        "Treatment-stratified GeoMx analyses, the complete 18-gene forest, and E-/M-MDSC sensitivity "
        "analyses are retained as supplementary analyses.\n"
    )
    legend_path = write_root / "Figure3_legend.md"
    legend_path.write_text(content, encoding="utf-8")

    content = (
        "# Figure 3 visual-variety rebuild\n\n- A: 22-patient paired GeoMx dumbbell.\n- B: "
        "patient-by-shared-leading-edge gene delta fingerprint, with the same patient order as A.\n- C: "
        "two tumors by three single-score H&E-backed spatial maps.\n- D: four descriptive "
        "myeloid-adjusted spatial effects.\n- Protein and anti-PD-1 results remain outside Figure 3.\n- "
        "Manifest: `"
        f"{manifest_path}"
        "`\n- Legend: `"
        f"{legend_path}"
        "`\n"
    )
    (write_root / "FINAL_RESULT.md").write_text(content, encoding="utf-8")

    # Verify the written files are nonempty and the preview dimensions are correct.
    output_paths = [
        manifest_path,
        legend_path,
        write_root / "FINAL_RESULT.md",
        pdf_path,
        png_path,
        preview_gray_path,
        source_path,
    ]
    for panel_result in records:
        for field in ["source", "pdf", "png"]:
            output_paths.append(Path(panel_result[field]))
    for output_path in output_paths:

        # Check the actual files produced by this build for existence and nonzero size.
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty output: {output_path.name}")
    with Image.open(png_path) as preview_image:
        if preview_image.size != canvas.size:
            raise ValueError("Preview dimensions differ from the original canvas")
        preview_image.verify()

    # Read the written manifest and verify its shape against the panel records.
    manifest_check_table = pd.read_csv(manifest_path)
    if manifest_check_table.shape != manifest_table.shape:
        raise ValueError("Panel manifest dimensions changed during writing")
    print(f"Figure3 panels={len(records)} preview=1")
