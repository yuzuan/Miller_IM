#!/usr/bin/env python3
"""Assemble the five October NOA R1 main figures at the final print size.

Inputs: Panel and source-data inputs of 06d, 07, 08c, 08d and 08e.
    Figure 3C uses the archived Step 55 tissue raster and spot table under NOA_LEGACY_ROOT
    and the registered-array source files under figure3_scale_calibration/sources.
    Comparison sheets use subssion/Figures/QC/Upload_Ready_Previews/Figure{1–5}.png under NOA_LEGACY_ROOT.
Outputs: figures_R1/Figure{1–5}_R1.{png,tiff}, figures_R1/_compare_old_vs_new_Fig{1–5}.png,
    panel builds under figures_R1/_font_build/fig{1–5}, and
    figure3_scale_calibration/Figure3C_scale_calibration.csv when building Figure 3.
Run: python NOA_R1/scripts/09b_export_final_figures.py [--figures 1 2 3 4 5] [--skip-comparisons]
"""

import argparse
import importlib
import shutil
from functools import partial
from pathlib import Path

import matplotlib

# Run as a headless export using the regular Matplotlib cache.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.text import Text
from PIL import Image, ImageDraw, ImageOps


figure_tools = importlib.import_module("00a_functions")


def draw_tissue_panel(figure_module, project_root):
    """Draw the October vector labels and calibrated scale bars on the archived Figure 3C map."""
    stem = "Figure3C_spatial_score_atlas"
    panel_dir = Path("55_figure3_visual_polish/Figure3")
    image_path = project_root / "figures" / panel_dir / f"{stem}.png"
    if not image_path.is_file():
        raise FileNotFoundError(f"Figure 3C requires the archived tissue raster at {image_path}")

    # Keep the original tissue pixels and replace only the surrounding text bands.
    with Image.open(image_path) as original:
        image = original.convert("RGB")
    if image.size != (3354, 2267):
        raise ValueError(f"Figure 3C tissue raster must be 3354 × 2267 pixels, found {image.size}")
    drawer = ImageDraw.Draw(image)
    drawer.rectangle((0, 0, 3354, 263), fill="white")
    drawer.rectangle((0, 264, 400, 1925), fill="white")
    drawer.rectangle((0, 2035, 3354, 2267), fill="white")

    # Use the registered-array fit and archived color limits for the labels.
    source_path = figure_module.SOURCE_OUT / f"{stem}_source.csv.gz"
    original_source = project_root / "write" / panel_dir / "source_data" / source_path.name
    shutil.copy2(original_source, source_path)
    calibration = importlib.import_module("08b_figure3_scale_calibration")
    calibration_table = calibration.calibrate(source_path)
    calibration_table = calibration_table.set_index("sample")
    spot_table = pd.read_csv(source_path)

    # Draw in tissue-image pixels so all label positions match the archived raster.
    figure = plt.figure(figsize=(3354 / 600, 2267 / 600), dpi=600, facecolor="white")
    axis = figure.add_axes([0, 0, 1, 1])
    axis.imshow(image, extent=(0, 3354, 2267, 0), interpolation="nearest")
    axis.set_xlim(0, 3354)
    axis.set_ylim(2267, 0)
    axis.axis("off")

    def label(x, y, text, size=8, align="center", weight="normal") -> None:
        """Place one vector label on the tissue canvas. No return value."""
        axis.text(
            x, y, text,
            fontsize    = size,
            color       = "#222222",
            ha          = align,
            va          = "center",
            fontweight  = weight,
            linespacing = 1.05,
        )

    label(330, 80, "Spatial score maps in two primary GBMs", 10, "left", "bold")
    for x, text in ((800, "MG-Inflam20"), (1800, "MDSC-like"), (2800, "Mesenchymal/hypoxia\ncomposite")):
        label(x, 205, text, weight="bold")

    # Place the two section counts at their original positions above the 1 mm scale bars.
    for sample, label_y, bar_y in (("GBM030", 350, 610), ("GBM049", 1140, 1450)):
        spots = int(calibration_table.loc[sample, "displayed_spots"])
        label(195, label_y, f"{sample}\n{spots:,} spots")
        half_length = float(calibration_table.loc[sample, "bar_length_canvas_px"]) / 2
        ends = [195 - half_length, 195 + half_length]
        axis.plot(ends, [bar_y, bar_y], color="#222222", linewidth=1.2, solid_capstyle="butt")
        for x in ends:
            axis.plot([x, x], [bar_y - 9, bar_y + 9], color="#222222", linewidth=0.8)
        label(195, bar_y + 65, "1 mm")

    # The original color-bar text uses the lower, middle and upper display limits.
    for x, column, digits, title in (
        (810, "raw20", 4, "MG-Inflam20 score"),
        (1810, "mdsc_combined", 3, "MDSC-like score"),
        (2810, "mes", 3, "Mesenchymal/hypoxia\ncomposite score"),
    ):
        limit_columns = [f"{column}_scale_low", f"{column}_scale_high"]
        color_limit_table = spot_table.loc[:, limit_columns].drop_duplicates()
        if len(color_limit_table) != 1:
            raise ValueError(f"Figure 3C requires shared display limits for {column}")
        lower, upper = color_limit_table.iloc[0]

        for shift, value in zip((-325, 0, 325), (lower, (lower + upper) / 2, upper)):
            number = f"{value:.{digits}f}".replace("-", "−")
            label(x + shift, 2095, number)
        label(x, 2205, title)

    pdf_path, png_path = figure_module.save_panel(figure, stem)
    panel_record = figure_module.base.record(
        "C", stem, "Frozen tissue map; registered-array 1 mm scale bars and vector labels",
        source_path, pdf_path, png_path,
    )
    return panel_record


def build_figure_preview(number, output_dir, project_root) -> Path:
    """Run the selected numbered panel module and return its assembled preview path."""
    figure_modules = {
        1: ("06d_figure1", "Figure1_final_preview.png"),
        2: ("07_figure2", "Figure2_two_panel_preview.png"),
        3: ("08c_figure3", "Figure3_visual_polish_preview.png"),
        4: ("08d_figure4", "Figure4_proteome_preview.png"),
        5: ("08e_figure5", "Figure5_pd1_preview.png"),
    }
    module_name, preview_name = figure_modules[number]
    figure_module = importlib.import_module(module_name)

    # Keep the panel modules' output bindings for the existing isolated build tree.
    target = output_dir / "_font_build" / f"fig{number}"
    figure_dir = target / f"Figure{number}"
    figure_module.WRITE_ROOT = target / "data" / f"Figure{number}"
    figure_module.SOURCE_OUT = figure_module.WRITE_ROOT / "source_data"
    if number == 2:
        figure_module.FIG_ROOT = figure_dir
    else:
        figure_module.FIG_OUT = figure_dir
    figure_dir.mkdir(parents=True, exist_ok=True)
    figure_module.SOURCE_OUT.mkdir(parents=True, exist_ok=True)

    # The base module also writes panel A into the redirected Figure 3 directories.
    if number == 3:
        figure_module.base.WRITE_ROOT = figure_module.WRITE_ROOT
        figure_module.base.SOURCE_OUT = figure_module.SOURCE_OUT
        figure_module.base.FIG_OUT = figure_dir
        figure_module.base.save_panel = figure_module.save_panel
        figure_module.panel_c_from_original = partial(draw_tissue_panel, figure_module, project_root)

    # Figures 1 and 2 expose the October geometry as an explicit layout option.
    if number in (1, 2):
        figure_module.main(final_layout=True)
    else:
        figure_module.main()

    preview_path = figure_dir / preview_name
    return preview_path


def export_figure(number, preview_path, output_dir, comparison_path) -> None:
    """Resize the assembled panel preview to 170 mm and write the submission raster formats."""
    width = round(figure_tools.mm(170) * 600)
    with Image.open(preview_path) as preview:
        image = preview.convert("RGB")
    height = round(image.height * width / image.width)
    image = image.resize((width, height), Image.Resampling.LANCZOS)
    if height > round(figure_tools.mm(220) * 600):
        raise ValueError(f"Figure {number} exceeds the 220 mm height limit: {image.size}")

    # The data panels already have PDF masters. These raster formats resize the assembled preview.
    stem = f"Figure{number}_R1"
    png_path = output_dir / f"{stem}.png"
    tiff_path = output_dir / f"{stem}.tiff"
    image.save(png_path, dpi=(600, 600), optimize=True)
    image.save(tiff_path, dpi=(600, 600), compression="tiff_lzw")

    # Retain the comparison sheet's original image sizes, labels and positions.
    if comparison_path is not None:
        with Image.open(comparison_path) as original:
            original_image = original.convert("RGB")
        before = ImageOps.contain(original_image, (1800, 2350), Image.Resampling.LANCZOS)
        after = ImageOps.contain(image, (1800, 2350), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (3720, max(before.height, after.height) + 130), "white")
        drawer = ImageDraw.Draw(canvas)
        font = figure_tools.image_font(48, bold=True)
        drawer.text((30, 22), "Original submission", fill="#222222", font=font)
        drawer.text((1890, 22), "R1 larger type", fill="#222222", font=font)
        canvas.paste(before, ((1800 - before.width) // 2, 105))
        canvas.paste(after, (1860 + (1800 - after.width) // 2, 105))
        preview_comparison_path = output_dir / f"_compare_old_vs_new_Fig{number}.png"
        canvas.save(preview_comparison_path, optimize=True)
        with Image.open(preview_comparison_path) as saved:
            if saved.size != canvas.size:
                raise ValueError(f"Unexpected Figure {number} comparison dimensions")
            saved.verify()

    # Verify both submission rasters before reporting the actual dimensions.
    for output_path in (png_path, tiff_path):
        with Image.open(output_path) as saved:
            if saved.size != image.size or saved.mode != "RGB":
                raise ValueError(f"Unexpected output size or mode in {output_path.name}")
            saved.verify()
    print(f"Figure{number}: {width}x{height} pixels; {width / 600 * 25.4:.2f}x{height / 600 * 25.4:.2f} mm")


def main() -> None:
    """Build the selected figures and restore the original save method on exit."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", nargs="+", type=int, choices=range(1, 6), default=list(range(1, 6)))
    parser.add_argument(
        "--skip-comparisons",
        action = "store_true",
        help   = "Build final figures without the archived original-submission comparison previews.",
    )
    args = parser.parse_args()

    # The actual panel builds write only under the revision root.
    revision_dir = Path(__file__).resolve().parents[1]
    project_root = figure_tools.legacy_root(revision_dir)
    output_dir = revision_dir / "figures_R1"

    # Check the original review previews before starting any requested comparison export.
    comparison_paths = {}
    if not args.skip_comparisons:
        original_dir = project_root / "subssion/Figures/QC/Upload_Ready_Previews"
        for number in args.figures:
            comparison_path = original_dir / f"Figure{number}.png"
            if not comparison_path.is_file():
                raise FileNotFoundError(f"Missing original preview at {comparison_path}. Use --skip-comparisons.")
            comparison_paths[number] = comparison_path

    # Use the original print-size factors only for the matching panel stems.
    panel_scales = {
        "Figure1D": 1.18, "Figure1E": 1.15, "Figure1F": 1.06,
        "Figure2_state": 1.00, "Figure2_raw20": 1.00,
        "Figure3A": 1.22, "Figure3B": 1.40, "Figure3D": 1.36,
        "Figure4A": 1.00, "Figure4B": 1.23, "Figure4C": 1.00,
        "Figure5A": 1.10, "Figure5B": 1.38, "Figure5C": 1.08,
    }
    original_savefig = Figure.savefig

    def save_with_print_labels(figure, fname, *args, **kwargs):
        """Save a named panel using the print layout and Matplotlib's fname interface."""
        stem = Path(fname).stem
        factor = 1.0
        for prefix, scale in panel_scales.items():
            if stem.startswith(prefix):
                factor = scale
                break

        # Apply the final label if the panel contains the older model title.
        if stem.startswith("Figure5C"):
            labels = [
                item for item in figure.findobj(match=Text) if item.get_text() == "Myeloid-proxy-adjusted model"
            ]
            for item in labels:
                item.set_text("Adjusted model")

        # Draw once to materialize tick labels. Scale text below 7.8 pt only once per figure.
        if factor != 1 and not getattr(figure, "_r1b_scaled", False):
            figure.canvas.draw()
            labels = [
                item for item in figure.findobj(match=Text)
                if item.get_text().strip() and item.get_fontsize() < 7.8
            ]
            for item in labels:
                item.set_fontsize(item.get_fontsize() * factor)
            figure._r1b_scaled = True

        # Keep the original extra line spacing on multiline labels in Figures 1 and 2.
        if stem.startswith(("Figure1", "Figure2")):
            labels = [item for item in figure.findobj(match=Text) if "\n" in item.get_text()]
            for item in labels:
                item.set_linespacing(1.65)

        result = original_savefig(figure, fname, *args, **kwargs)
        return result

    Figure.savefig = save_with_print_labels
    try:
        for number in args.figures:
            preview_path = build_figure_preview(number, output_dir, project_root)
            export_figure(number, preview_path, output_dir, comparison_paths.get(number))
    finally:
        Figure.savefig = original_savefig


if __name__ == "__main__":
    main()
