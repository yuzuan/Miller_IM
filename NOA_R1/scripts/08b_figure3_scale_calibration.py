"""Calculate Figure 3C scale bars from the registered Visium array.

Inputs
    The frozen Figure3C_spatial_score_atlas_source.csv.gz passed by the caller.
    figure3_scale_calibration/sources/{GSM8506649_03422,GSM8506650_25526}_tissue_positions.csv.gz
    and the matching *_scalefactors_json.json.gz files for full-resolution registration.
Output
    figure3_scale_calibration/Figure3C_scale_calibration.csv, one row per tissue section.
Run
    The Figure 3 step of scripts/09b_export_final_figures.py calls calibrate(frozen_source).
"""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def calibrate(frozen_source):
    """Return the calibration table after checking registration and writing the CSV."""
    displayed_spots_table = pd.read_csv(frozen_source)
    calibration_dir = Path(__file__).resolve().parents[1] / "figure3_scale_calibration"
    sample_prefixes = {"GBM030": "GSM8506649_03422", "GBM049": "GSM8506650_25526"}
    calibration_rows = []

    # Use the original 142 by 96 mm canvas at 600 dpi and the original GridSpec spacing.
    canvas_width = 142 / 25.4 * 600
    canvas_height = 96 / 25.4 * 600
    box_width = (0.98 - 0.10) * canvas_width / (3 + 2 * 0.05)
    box_height = (0.88 - 0.18) * canvas_height / (2 + 0.08)

    for sample, prefix in sample_prefixes.items():
        positions_path = calibration_dir / "sources" / f"{prefix}_tissue_positions.csv.gz"
        position_table = pd.read_csv(positions_path)
        scale_path = calibration_dir / "sources" / f"{prefix}_scalefactors_json.json.gz"
        with gzip.open(scale_path, "rt") as handle:
            scalefactors = json.load(handle)

        # The complete CytAssist array contains 4,992 distinct positions in each source file.
        if len(position_table) != 4992 or position_table["barcode"].duplicated().any():
            raise ValueError(f"{sample} needs 4,992 unique Visium array positions")
        section_table = displayed_spots_table.loc[displayed_spots_table["sample"].eq(sample)].copy()
        if section_table.empty:
            raise ValueError(f"{sample} is missing from the displayed spot table")

        # Adjacent spot centers in a row are 100 micrometers apart and differ by two array columns.
        # Staggered rows differ in height by sqrt(3) * 50 micrometers.
        array_x = position_table["array_col"] * 50
        array_y = position_table["array_row"] * np.sqrt(3) * 50
        array_microns = np.column_stack([array_x, array_y, np.ones(len(position_table))])
        fullres_pixels = position_table[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        coefficients, _, _, _ = np.linalg.lstsq(array_microns, fullres_pixels, rcond=None)
        fitted_pixels = (
            array_microns[:, 0, None] * coefficients[0]
            + array_microns[:, 1, None] * coefficients[1]
            + coefficients[2]
        )
        residuals = np.linalg.norm(fitted_pixels - fullres_pixels, axis=1)
        max_residual = residuals.max()
        if not max_residual < 1.0:
            raise ValueError(f"{sample} full-resolution registration residual is {max_residual} pixels")

        # Every plotted spot must match the original array registration on its barcode.
        registered_table = section_table.merge(
            position_table, on="barcode", suffixes=("", "_raw"), validate="one_to_one"
        )
        if len(registered_table) != len(section_table):
            raise ValueError(f"{sample} displayed barcodes do not all match the array")
        for column in ["array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]:
            if not np.array_equal(registered_table[column], registered_table[f"{column}_raw"]):
                raise ValueError(f"{sample} differs in registered coordinates for {column}")

        hires_scale = float(scalefactors["tissue_hires_scalef"])
        if not np.allclose(section_table["tissue_hires_scalef"], hires_scale, rtol=0, atol=1e-12):
            raise ValueError(f"{sample} differs in high-resolution scale factor")
        for display_column, fullres_column in [("x_hires", "pxl_col_in_fullres"), ("y_hires", "pxl_row_in_fullres")]:
            expected_pixels = section_table[fullres_column] * hires_scale
            if not np.allclose(section_table[display_column], expected_pixels, rtol=0, atol=1e-9):
                raise ValueError(f"{sample} differs in scaled coordinates for {display_column}")

        # Preserve slight anisotropy by using the horizontal vector of the inverse affine transform.
        inverse = np.linalg.inv(coefficients[:2] * hires_scale)
        microns_per_hires_x_pixel = float(np.linalg.norm(inverse[0]))
        x_min, x_max = section_table["x_hires"].min(), section_table["x_hires"].max()
        y_min, y_max = section_table["y_hires"].min(), section_table["y_hires"].max()
        x_margin = (x_max - x_min) * 0.10
        y_margin = (y_max - y_min) * 0.10

        # Match the maps' 10% crop padding and aspect-preserving fit inside their original axes.
        x0_hires = max(0.0, x_min - x_margin)
        x1_hires = min(float(section_table["image_width_px"].iloc[0]), x_max + x_margin)
        y0_hires = max(0.0, y_min - y_margin)
        y1_hires = min(float(section_table["image_height_px"].iloc[0]), y_max + y_margin)
        canvas_per_hires = min(box_width / (x1_hires - x0_hires), box_height / (y1_hires - y0_hires))
        bar_um = 1000.0
        calibration_rows.append(
            {
                "sample": sample,
                "original_positions": len(position_table),
                "displayed_spots": len(section_table),
                "tissue_hires_scalef": hires_scale,
                "max_registration_residual_fullres_px": float(max_residual),
                "microns_per_horizontal_hires_pixel": microns_per_hires_x_pixel,
                "canvas_pixels_per_hires_pixel": float(canvas_per_hires),
                "bar_length_um": bar_um,
                "bar_length_canvas_px": bar_um / microns_per_hires_x_pixel * canvas_per_hires,
                "x0_hires": float(x0_hires),
                "x1_hires": float(x1_hires),
                "y0_hires": float(y0_hires),
                "y1_hires": float(y1_hires),
                "affine_fullres": json.dumps(coefficients.tolist()),
                "calibration": (
                    "registered 100-micrometre Visium CytAssist array pitch; horizontal inverse-affine metric"
                ),
            }
        )

    calibration_table = pd.DataFrame(calibration_rows)
    output_path = calibration_dir / "Figure3C_scale_calibration.csv"
    calibration_table.to_csv(output_path, index=False)
    written_table = pd.read_csv(output_path)
    if written_table.shape != calibration_table.shape or not written_table.columns.equals(calibration_table.columns):
        raise ValueError("The written Figure 3C calibration has unexpected rows or columns")
    return calibration_table
