#!/usr/bin/env python3
# Purpose: Draw six measurable genes, a fixed grayscale H&E background, and four Figure 3 spatial effects.
# Inputs: Paired GeoMx data, Step 41 spatial tables, and the two Visium images loaded by 07c.
# Outputs: source tables and Figure3_panel_manifest.csv under
# write/55_figure3_visual_polish/Figure3, with panels A–D and Figure3_visual_polish_preview under figures.
# Run python scripts/07d_figure3.py from the repository root.
# The Python workflow uses pandas, NumPy, Matplotlib, and Pillow for the released calculations and figures.
# Keep the public paths, fields, and data names. Render panel PNGs from PDF and compose previews with Pillow.

# Each panel stem below produces a PDF, PNG, and source_data/<stem>_source.csv.
#   Figure3A_geomx_paired_dumbbell
#   Figure3B_shared_gene_patient_fingerprint
#   Figure3C_spatial_score_atlas
#   Figure3D_partial_rho_effect_strip
# Panel C uses .csv.gz. Preview filenames use the _grayscale.png and _layout.csv suffixes for those assets.

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib as mpl
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

# Load the shared GeoMx and Visium functions used by this figure.
root = Path(__file__).resolve().parents[1]
base_script = root / "scripts/07c_spatial_and_paired_expression.py"
spec = importlib.util.spec_from_file_location("step53_figure3", base_script)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load the shared Figure 3 functions: {base_script}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

# Initialize in the original order with the existing global seed. Pass output directories to the imported functions.
apply_publication_style()
np.random.seed(20260714)
write_root = root / "write/55_figure3_visual_polish/Figure3"
source_out = write_root / "source_data"
source_out.mkdir(parents=True, exist_ok=True)
figure_out = root / "figures/55_figure3_visual_polish/Figure3"
figure_out.mkdir(parents=True, exist_ok=True)
external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
spatial_image_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
if not spatial_image_root.exists():
    raise FileNotFoundError(f"Visium tissue image directory is unavailable: {spatial_image_root}")

# Use the shared paired-patient plot, supplying the output directories and PDF metadata.
patients_table = base.patient_table()
records = [base.panel_a(patients_table, source_out, figure_out, "Step55 recurrent GBM Figure 3 visual polish")]

# Keep seven prespecified genes in the source table and six measurable genes in panel B, with the same order, FDR and
# scale.

genes = base.shared_genes()
gene_delta_table = base.strict_gene_delta(patients_table, genes)
if len(gene_delta_table) != 154 or int(gene_delta_table["measured"].sum()) != 132:
    raise ValueError("Figure 3B requires the source rows for all seven prespecified genes, with six measurable genes")
if gene_delta_table.loc[~gene_delta_table["measured"], "gene"].drop_duplicates().tolist() != ["CCL4"]:
    raise ValueError("CCL4 must be the only unmeasured gene in Figure 3B")
stem = "Figure3B_shared_gene_patient_fingerprint"

# Write the panel source table with the existing columns in their original order.
source_path = source_out / f"{stem}_source.csv"
gene_delta_table.to_csv(source_path, index=False)


# Select the six measured genes for display and require three changes beyond the 98th-percentile color limit.
measured_genes = gene_delta_table.loc[gene_delta_table["measured"], "gene"].drop_duplicates().tolist()
if measured_genes != genes[:-1]:
    raise ValueError(f"Unexpected order of measured genes: {measured_genes}")
delta_pivot_table = (
    gene_delta_table.loc[gene_delta_table["measured"]]
    .pivot(index="gene", columns="patient_label", values="gene_delta")
    .reindex(index=measured_genes, columns=patients_table["patient_label"].tolist())
)
color_limit = float(gene_delta_table["color_limit"].iloc[0])
clipped = int((gene_delta_table.loc[gene_delta_table["measured"], "gene_delta"].abs() > color_limit).sum())
if clipped != 3:
    raise ValueError(f"Unexpected count beyond the 98th-percentile color limit in Figure 3B: {clipped}")

# Use the original figure size in millimeters and the same PDF axes positions.
figure = new_figure(116, 86)
axis = figure.add_axes([0.12, 0.29, 0.66, 0.56])
image = axis.imshow(
    delta_pivot_table.to_numpy(dtype=float),
    cmap   = mpl.colormaps["RdBu_r"],
    vmin   = -color_limit,
    vmax   = color_limit,
    aspect = "auto",
)
axis.set_yticks(np.arange(len(measured_genes)), measured_genes, fontsize=6.6)
axis.set_xticks(
    np.arange(len(patients_table)),
    patients_table["patient_label"],
    rotation = 90,
    fontsize = 4.9,
)
axis.tick_params(length=0, pad=2)
for spine in axis.spines.values():
    spine.set_visible(False)
axis.set_xticks(np.arange(-0.5, len(patients_table), 1), minor=True)
axis.set_yticks(np.arange(-0.5, len(measured_genes), 1), minor=True)
axis.grid(which="minor", color="white", linewidth=0.35)
axis.tick_params(which="minor", bottom=False, left=False)

# The color bar shows the display scale for the gene changes.
color_axis = figure.add_axes([0.19, 0.16, 0.49, 0.023])
colorbar = figure.colorbar(image, cax=color_axis, orientation="horizontal")
colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.8, labelpad=1.5)
colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
colorbar.outline.set_visible(False)

# Align the original summary statistics to the order of the measurable genes.
gene_summary_table = (
    gene_delta_table.loc[gene_delta_table["measured"]].drop_duplicates("gene").set_index("gene").reindex(measured_genes)
)
summary_axis = figure.add_axes([0.81, 0.29, 0.16, 0.56])
summary_axis.axvline(0, color="#A0A0A0", linewidth=0.65)
for index, gene in enumerate(measured_genes):
    row = gene_summary_table.loc[gene]
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
summary_axis.set_ylim(len(measured_genes) - 0.5, -0.5)
summary_axis.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.9)
summary_axis.set_yticks([])
summary_axis.set_xlabel("Mean Δ", fontsize=5.6, labelpad=2)
summary_axis.grid(axis="x", color="#E7E7E7", linewidth=0.45)
clean_axis(summary_axis, keep_left=False, keep_bottom=True)

figure.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=9.0, fontweight="bold", va="top")

# Save the panel as PDF and render its 600 dpi PNG from that file.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step55 recurrent GBM Figure 3 visual polish")
records.append(
    panel_record(
        "B",
        stem,
        "Six measurable genes shown from seven prespecified shared leading-edge genes",
        source_path,
        pdf_path,
        png_path,
    )
)


# In the shared panel C function, select the original Step 55 grayscale conversion and use the three existing score
# scales.
records.append(base.panel_c(source_out, figure_out, "Step55 recurrent GBM Figure 3 visual polish"))

# Draw the four existing partial correlations against zero on the original 0–0.5 scale.

partial_path = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3e_partial_rho_main.csv"
partial_source_table = pd.read_csv(partial_path)
partial_table = partial_source_table[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
partial_table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
if len(partial_table) != 4 or not (partial_table["partial_rho"] > 0).all():
    raise ValueError("Figure 3D requires four positive effects for the MDSC and MES scores in two sections")
stem = "Figure3D_partial_rho_effect_strip"

# Write the panel source table with the existing columns in their original order.
source_path = source_out / f"{stem}_source.csv"
partial_table.to_csv(source_path, index=False)

# Use the original figure size in millimeters and the same PDF axes positions.
figure = new_figure(58, 68)
axis = figure.add_axes([0.31, 0.25, 0.64, 0.55])
target_y = {"MDSC": 1.0, "MES": 0.0}
sample_offset = {"GBM030": 0.105, "GBM049": -0.105}
sample_marker = {"GBM030": "o", "GBM049": "s"}
target_color = {"MDSC": "#159D88", "MES": "#7E6148"}

for target in ["MDSC", "MES"]:
    target_table = partial_table.loc[partial_table["target"].eq(target)].set_index("sample")
    for sample in ["GBM030", "GBM049"]:
        x_position = float(target_table.loc[sample, "partial_rho"])
        y_position = target_y[target] + sample_offset[sample]
        axis.hlines(y_position, 0, x_position, color=target_color[target], linewidth=1.15, alpha=0.55, zorder=1)
        axis.scatter(
            x_position,
            y_position,
            s         = 34,
            marker    = sample_marker[sample],
            facecolor = target_color[target],
            edgecolor = "white",
            linewidth = 0.6,
            zorder    = 3,
        )
        axis.text(x_position + 0.014, y_position, f"{x_position:.3f}", fontsize=5.9, va="center", color="#222222")

axis.axvline(0, color="#8E8E8E", linewidth=0.75)

# Keep the horizontal scale and labels for the original comparison.
axis.set_xlim(0, 0.50)

# Use the original plot limits and data order.
axis.set_ylim(-0.46, 1.46)
axis.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.5)
axis.set_xticks([0, 0.2, 0.4])
axis.set_xlabel(r"Partial Spearman $r_s$", fontsize=6.2)
axis.grid(axis="x", color="#E7E7E7", linewidth=0.5)
clean_axis(axis, keep_left=False, keep_bottom=True)
axis.tick_params(axis="y", length=0, pad=3)

# Keep the original legend order and colors and match the symbols in the panel.
handles = [
    mpl.lines.Line2D([], [], marker="o", linestyle="none", color="#222222", markerfacecolor="white", label="GBM030"),
    mpl.lines.Line2D([], [], marker="s", linestyle="none", color="#222222", markerfacecolor="white", label="GBM049"),
]
axis.legend(
    handles        = handles,
    frameon        = False,
    fontsize       = 5.5,
    loc            = "upper center",
    bbox_to_anchor = (0.52, -0.30),
    ncol           = 2,
    handletextpad  = 0.3,
    columnspacing  = 0.8,
)
figure.text(0.12, 0.955, "Miller-IM score vs niche", fontsize=9.2, fontweight="bold", va="top")

# Save the panel as PDF and render its 600 dpi PNG from that file.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step55 recurrent GBM Figure 3 visual polish")
records.append(panel_record("D", stem, "Compact descriptive partial-rho effect strip", source_path, pdf_path, png_path))


# Use the 4320 × 4370 pixel preview and the existing panel rectangles.

canvas = Image.new("RGB", (4320, 4370), "white")
drawer = ImageDraw.Draw(canvas)
section_font = image_font(62, bold=True)

drawer.text((115, 55), "Independent recurrence support", fill="#E64B35", font=section_font)
drawer.line((115, 145, 4205, 145), fill="#E8B2AA", width=5)
drawer.text((115, 2190), "Miller-IM program in MDSC/MES-rich niches", fill="#159D88", font=section_font)
drawer.line((115, 2285, 4205, 2285), fill="#A8D9CF", width=5)

if [row["panel"] for row in records] != list("ABCD"):
    raise ValueError("Figure 3 requires panels A through D in order")

by_panel = {row["panel"]: Path(row["png"]) for row in records}
layout = [
    {"panel": "A", "left": 115, "top": 270, "right": 1885, "bottom": 2110},
    {"panel": "B", "left": 1975, "top": 270, "right": 4205, "bottom": 2110},
    {"panel": "C", "left": 115, "top": 2450, "right": 2465, "bottom": 4040},
    {"panel": "D", "left": 2610, "top": 2450, "right": 4205, "bottom": 4235},
]
panel_font = image_font(78, bold=True)
for item in layout:
    panel = item["panel"]
    paste_contained(canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"]))
    drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

# Save the Pillow preview with the original pixel size and resolution for layout review.
png_path = figure_out / "Figure3_visual_polish_preview.png"
canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
preview_gray_path = figure_out / "Figure3_visual_polish_preview_grayscale.png"
ImageOps.grayscale(canvas).convert("RGB").save(preview_gray_path, format="PNG", dpi=(600, 600), optimize=True)
pdf_path = figure_out / "Figure3_visual_polish_preview.pdf"
canvas.save(pdf_path, format="PDF", resolution=600.0, creationDate="", modDate="")
layout_table = pd.DataFrame(layout)
source_path = source_out / "Figure3_visual_polish_preview_layout.csv"
layout_table.to_csv(source_path, index=False)


# Write one manifest row per panel and retain the four preview references.
manifest_table = pd.DataFrame(records)
manifest_table["preview_pdf"] = str(pdf_path)
manifest_table["preview_png"] = str(png_path)
manifest_table["preview_grayscale"] = str(preview_gray_path)
manifest_table["preview_layout"] = str(source_path)
manifest_path = write_root / "Figure3_panel_manifest.csv"
manifest_table.to_csv(manifest_path, index=False)

# Verify the output files, preview dimensions, and manifest shape. Stop if an output is missing or empty.
output_paths = [
    manifest_path,
    pdf_path,
    png_path,
    preview_gray_path,
    source_path,
]
for panel_result in records:
    for field in ["source", "pdf", "png"]:
        output_paths.append(Path(panel_result[field]))
for output_path in output_paths:

    # Check that each file from the current build is present and nonempty.
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
