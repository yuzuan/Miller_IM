#!/usr/bin/env python3
# Purpose: Supply the shared GeoMx and Visium functions used by 07d_figure3.py.
# Inputs: Step 41 GeoMx and Visium source tables and the Step 51 shared leading-edge genes.
# Use GeoMx metadata and normalized expression and the registered Visium images under MILLER_IM_DATA_ROOT.
# Outputs: Figure 3A and 3C panel files and source tables in the directories passed by 07d.
# Run: python scripts/07d_figure3.py from the repository root.

from __future__ import annotations

import gzip
import io
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from PIL import Image

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
panel_record = figure_style.panel_record
save_panel = figure_style.save_panel

# Keep one row per strict GeoMx patient pair in the frozen order.
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


# Draw panel A using the frozen paired patient table.
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


# Read the shared leading-edge genes in the original Figure 1 order.
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


# Use the original ROI eligibility and log2(x + 1) values to compute paired changes.
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


def panel_c(source_out: Path, figure_out: Path, creator: str) -> dict[str, str]:
    """Write panel C over grayscale tissue images and return its manifest record."""
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

    # Apply the same low-contrast grayscale conversion to both sections.
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
        luminance = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
        display = 0.68 + 0.32 * np.clip(luminance, 0.0, 1.0)
        images[sample] = np.repeat(display[:, :, None], 3, axis=2)

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
