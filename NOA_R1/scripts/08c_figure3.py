#!/usr/bin/env python3
"""Render Figure 3 with paired GeoMx observations and patient-matched spatial maps.

Inputs
    Under NOA_LEGACY_ROOT, use the paired scores, shared-gene definition and partial rho table
    read by scripts/08a_figure3_spatial_panels.py, plus the archived Step 55 assets:
    write/55_figure3_visual_polish/Figure3/source_data/Figure3B_shared_gene_patient_fingerprint_source.csv
    write/55_figure3_visual_polish/Figure3/source_data/Figure3C_spatial_score_atlas_source.csv.gz
    figures/55_figure3_visual_polish/Figure3/Figure3C_spatial_score_atlas.png
Outputs
    Under _figure_build/55_figure3_visual_polish/Figure3, write Figure3_panel_manifest.csv,
    Figure3_legend.md, FINAL_RESULT.md and source_data/:
        Figure3A_geomx_paired_dumbbell_source.csv, Figure3B_shared_gene_patient_fingerprint_source.csv,
        Figure3C_spatial_score_atlas_source.csv.gz, Figure3D_partial_rho_effect_strip_source.csv,
        Figure3_visual_polish_preview_layout.csv.
    Under _panels/55_figure3_visual_polish/Figure3, write 600 dpi PNGs and PDF masters for
    Figure3A_geomx_paired_dumbbell, Figure3B_shared_gene_patient_fingerprint, Figure3C_spatial_score_atlas,
    Figure3D_partial_rho_effect_strip and Figure3_visual_polish_preview.
    Also write Figure3_visual_polish_preview_grayscale.png.
Run
    python scripts/08c_figure3.py from NOA_R1, or use scripts/09b_export_final_figures.py
    for the final Figure 3C with calibrated scale bars and vector annotations.
"""

from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps


figure_tools = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
project_root = figure_tools.legacy_root(revision_root)

# These output bindings are overridden by the final exporter.
WRITE_ROOT = revision_root / "_figure_build/55_figure3_visual_polish/Figure3"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "_panels/55_figure3_visual_polish/Figure3"


def load_base():
    """Load a local instance of 08a and set the output directories used by its panel A."""
    base_path = revision_root / "scripts/08a_figure3_spatial_panels.py"
    spec = importlib.util.spec_from_file_location("step53_figure3", base_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load the Figure 3 panel script at {base_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.WRITE_ROOT = WRITE_ROOT
    module.SOURCE_OUT = SOURCE_OUT
    module.FIG_OUT = FIG_OUT
    return module


base = load_base()


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    """Save a PDF master and 600 dpi PNG at the location selected by the exporter."""
    output_paths = figure_tools.save_panel(fig, FIG_OUT, stem, "Step55 recurrent GBM Figure 3 visual polish")
    return output_paths


base.save_panel = save_panel


def fixed_grayscale_tissue(rgb: np.ndarray) -> np.ndarray:
    """Return low-contrast grayscale RGB using a fixed conversion for both sections."""
    luminance = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
    display = 0.68 + 0.32 * np.clip(luminance, 0.0, 1.0)
    grayscale = np.repeat(display[:, :, None], 3, axis=2)
    return grayscale


base.soften_tissue = fixed_grayscale_tissue


def panel_b(patient_table: pd.DataFrame) -> dict[str, str]:
    """Write the frozen patient-by-gene table and render the measurable genes in panel B."""
    stem = "Figure3B_shared_gene_patient_fingerprint"
    archived_source = project_root / "write/55_figure3_visual_polish/Figure3/source_data" / f"{stem}_source.csv"
    fingerprint_table = pd.read_csv(archived_source)
    if len(fingerprint_table) != 154 or int(fingerprint_table["measured"].sum()) != 132:
        raise ValueError("Figure 3B needs all seven prespecified genes and six measured genes")
    if fingerprint_table.loc[~fingerprint_table["measured"], "gene"].drop_duplicates().tolist() != ["CCL4"]:
        raise ValueError("CCL4 must be the only unmeasured gene in Figure 3B")
    source = SOURCE_OUT / f"{stem}_source.csv"
    fingerprint_table.to_csv(source, index=False)

    # Keep one heatmap row per measured gene and the same patient ordering as the archived plot.
    genes = base.shared_genes()
    measured_table = fingerprint_table.loc[fingerprint_table["measured"]]
    measured_genes = measured_table["gene"].drop_duplicates().tolist()
    if measured_genes != genes[:-1]:
        raise ValueError(f"Unexpected order of measured genes: {measured_genes}")
    expression_table = measured_table.pivot(index="gene", columns="patient_label", values="gene_delta")
    expression_table = expression_table.reindex(index=measured_genes, columns=patient_table["patient_label"])
    color_limit = float(fingerprint_table["color_limit"].iloc[0])
    clipped = int((measured_table["gene_delta"].abs() > color_limit).sum())
    if clipped != 3:
        raise ValueError(f"Figure 3B needs three values beyond the frozen 98th-percentile color limit, got {clipped}")

    # Keep the shared color bar and each patient's recurrent-minus-primary log2 change.
    fig = base.new_figure(116, 86)
    ax = fig.add_axes([0.14, 0.29, 0.64, 0.56])
    image = ax.imshow(
        expression_table.to_numpy(dtype=float),
        cmap   = mpl.colormaps["RdBu_r"],
        vmin   = -color_limit,
        vmax   = color_limit,
        aspect = "auto",
    )
    ax.set_yticks(np.arange(len(measured_genes)), measured_genes, fontsize=6.6)
    ax.set_xticks(
        np.arange(len(patient_table)),
        patient_table["patient_label"],
        rotation = 90,
        fontsize = 4.9,
    )
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(patient_table), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(measured_genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Label the symmetric color scale under the display.
    color_ax = fig.add_axes([0.19, 0.16, 0.49, 0.023])
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.8, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
    colorbar.outline.set_visible(False)

    # The right-side strip uses the existing confidence limits and marks FDR below 0.05.
    summary_table = measured_table.drop_duplicates("gene").set_index("gene").reindex(measured_genes)
    sx = fig.add_axes([0.81, 0.29, 0.16, 0.56])
    sx.axvline(0, color="#A0A0A0", linewidth=0.65)
    for index, gene in enumerate(measured_genes):
        row = summary_table.loc[gene]
        significant = float(row["gene_fdr"]) < 0.05
        sx.hlines(index, float(row["ci_low"]), float(row["ci_high"]), color=base.recurrent, linewidth=1.0)
        sx.scatter(
            float(row["mean_gene_delta"]),
            index,
            s         = 22,
            facecolor = base.recurrent if significant else "white",
            edgecolor = base.recurrent,
            linewidth = 0.8,
            zorder    = 3,
        )
    sx.set_xlim(-0.45, 0.92)
    sx.set_ylim(len(measured_genes) - 0.5, -0.5)
    sx.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.9)
    sx.set_yticks([])
    sx.set_xlabel("Mean Δ", fontsize=5.6, labelpad=2)
    sx.grid(axis="x", color=base.grid_color, linewidth=0.45)
    base.clean_axis(sx, keep_left=False, keep_bottom=True)

    fig.text(0.10, 0.955, "Shared leading-edge gene changes", fontsize=9.0, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = base.record(
        "B", stem, "Six measurable genes shown from seven prespecified shared leading-edge genes", source, pdf, png
    )
    return panel_record


def panel_c_from_original() -> dict[str, str]:
    """Reuse the frozen H&E map and source table. The final exporter replaces this callback."""
    stem = "Figure3C_spatial_score_atlas"
    original_source = project_root / "write/55_figure3_visual_polish/Figure3/source_data" / f"{stem}_source.csv.gz"
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    shutil.copy2(original_source, source)

    original = project_root / "figures/55_figure3_visual_polish/Figure3" / f"{stem}.png"
    with Image.open(original) as old:
        canvas = old.convert("RGB")

    # The label replacements are confined to the archived white text bands.
    draw = ImageDraw.Draw(canvas)
    font = figure_tools.image_font(42, bold=True)
    label_font = figure_tools.image_font(39)
    draw.rectangle((580, 175, 1010, 263), fill="white")
    draw.text((795, 195), "MG-Inflam20", fill="#222222", font=font, anchor="mm")
    draw.rectangle((2600, 175, 3075, 264), fill="white")
    draw.text((2838, 187), "Mesenchymal/hypoxia", fill="#222222", font=font, anchor="mm")
    draw.text((2838, 234), "composite", fill="#222222", font=font, anchor="mm")
    draw.rectangle((575, 2095, 1050, 2198), fill="white")
    draw.text((812, 2150), "MG-Inflam20 score", fill="#222222", font=label_font, anchor="mm")
    draw.rectangle((2510, 2095, 3150, 2198), fill="white")
    draw.text(
        (2830, 2150),
        "Mesenchymal/hypoxia composite score",
        fill   = "#222222",
        font   = figure_tools.image_font(30),
        anchor = "mm",
    )

    # Embed the original raster at native scale in a 600 dpi PDF canvas.
    fig = plt.figure(figsize=(canvas.width / 600, canvas.height / 600), dpi=600, facecolor="white")
    axis = fig.add_axes([0, 0, 1, 1])
    axis.imshow(canvas, interpolation="none")
    axis.axis("off")
    pdf, png = save_panel(fig, stem)
    panel_record = base.record("C", stem, "Frozen tissue-backed maps with R1 labels", source, pdf, png)
    return panel_record


def panel_d() -> dict[str, str]:
    """Write and plot the four archived myeloid-adjusted partial correlations."""
    stem = "Figure3D_partial_rho_effect_strip"
    rho_path = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3e_partial_rho_main.csv"
    original_table = pd.read_csv(rho_path)
    association_table = original_table[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
    association_table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
    if len(association_table) != 4 or not (association_table["partial_rho"] > 0).all():
        raise ValueError("Figure 3D needs four positive effects for two targets and two sections")
    source = SOURCE_OUT / f"{stem}_source.csv"
    association_table.to_csv(source, index=False)

    # Target colors and sample markers match the frozen effect strip.
    fig = base.new_figure(58, 68)
    ax = fig.add_axes([0.40, 0.25, 0.55, 0.55])
    target_y = {"MDSC": 1.0, "MES": 0.0}
    sample_offset = {"GBM030": 0.105, "GBM049": -0.105}
    sample_marker = {"GBM030": "o", "GBM049": "s"}
    target_color = {"MDSC": base.mdsc_color, "MES": base.mes_color}

    for target in ["MDSC", "MES"]:
        target_table = association_table.loc[association_table["target"].eq(target)].set_index("sample")
        for sample in ["GBM030", "GBM049"]:
            effect = float(target_table.loc[sample, "partial_rho"])
            row_y = target_y[target] + sample_offset[sample]
            ax.hlines(row_y, 0, effect, color=target_color[target], linewidth=1.15, alpha=0.55, zorder=1)
            ax.scatter(
                effect,
                row_y,
                s         = 34,
                marker    = sample_marker[sample],
                facecolor = target_color[target],
                edgecolor = "white",
                linewidth = 0.6,
                zorder    = 3,
            )
            ax.text(effect + 0.014, row_y, f"{effect:.3f}", fontsize=5.9, va="center", color=base.text_color)

    # Keep the zero reference, the 0–0.50 axis and the two section labels.
    ax.axvline(0, color="#8E8E8E", linewidth=0.75)
    ax.set_xlim(0, 0.50)
    ax.set_ylim(-0.46, 1.46)
    ax.set_yticks([1, 0], ["MDSC-like", "Mesenchymal/\nhypoxia\ncomposite"], fontsize=6.5)
    ax.set_xticks([0, 0.2, 0.4])
    ax.set_xlabel(r"Partial Spearman $r_s$", fontsize=6.2)
    ax.grid(axis="x", color=base.grid_color, linewidth=0.5)
    base.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles = []
    for sample, marker in sample_marker.items():
        handles.append(
            mpl.lines.Line2D(
                [], [], marker=marker, linestyle="none", color=base.text_color, markerfacecolor="white", label=sample
            )
        )
    ax.legend(
        handles        = handles,
        frameon        = False,
        fontsize       = 5.5,
        loc            = "upper center",
        bbox_to_anchor = (0.52, -0.30),
        ncol           = 2,
        handletextpad  = 0.3,
        columnspacing  = 0.8,
    )
    fig.text(0.12, 0.955, "MG-Inflam20 score vs niche", fontsize=9.2, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = base.record("D", stem, "Compact descriptive partial-rho effect strip", source, pdf, png)
    return panel_record


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    """Compose the existing image previews from the saved panel PNGs."""
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 115, "top": 270, "right": 1885, "bottom": 2110},
        {"panel": "B", "left": 1975, "top": 270, "right": 4205, "bottom": 2110},
        {"panel": "C", "left": 115, "top": 2450, "right": 2650, "bottom": 4350},
        {"panel": "D", "left": 2750, "top": 2450, "right": 4205, "bottom": 4350},
    ]
    canvas = Image.new("RGB", (4320, 4450), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = figure_tools.image_font(78, bold=True)
    section_font = figure_tools.image_font(62, bold=True)

    # Preserve the section titles and panel rectangles used in the author preview.
    drawer.text((115, 55), "Paired myeloid-enriched tissue measurements", fill=base.recurrent, font=section_font)
    drawer.line((115, 145, 4205, 145), fill="#E8B2AA", width=5)
    drawer.text((115, 2190), "MG-Inflam20 program in spatial niches", fill=base.mdsc_color, font=section_font)
    drawer.line((115, 2285, 4205, 2285), fill="#A8D9CF", width=5)

    for item in layout:
        panel = item["panel"]
        figure_tools.paste_contained(
            canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"])
        )
        drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

    # Retain the color and grayscale image previews alongside the composite PDF.
    png = FIG_OUT / "Figure3_visual_polish_preview.png"
    gray = FIG_OUT / "Figure3_visual_polish_preview_grayscale.png"
    pdf = FIG_OUT / "Figure3_visual_polish_preview.pdf"
    source = SOURCE_OUT / "Figure3_visual_polish_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(600, 600), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(pdf, format="PDF", resolution=600.0, creationDate="", modDate="")
    layout_table = pd.DataFrame(layout)
    layout_table.to_csv(source, index=False)
    output_paths = (pdf, png, gray, source)
    return output_paths


def main() -> None:
    """Render A–D and write the existing source tables, previews, and manifest."""
    base.apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    # Use the prepared GeoMx table for panels A and B and the two section records for C and D.
    patient_table = base.patient_table()
    records = [base.panel_a(patient_table), panel_b(patient_table), panel_c_from_original(), panel_d()]
    if [row["panel"] for row in records] != list("ABCD"):
        raise ValueError("Figure 3 requires panels A through D in order")

    # Keep the manifest schema used by the exporter and the existing preview layout.
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)
    manifest_table["preview_layout"] = str(preview_layout)
    manifest_path = WRITE_ROOT / "Figure3_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)

    # Retain the frozen Step 55 standalone figure legend.
    legend_path = WRITE_ROOT / "Figure3_legend.md"
    legend_text = (
        "# Figure 3 legend\n\n"
        "**Figure 3 | Paired myeloid-enriched tissue support and patient-matched spatial context of the "
        "MG-Inflam20 program.**\n\n"
        "**A,** Paired IBA1-positive GeoMx MG-Inflam20 scores in 22 IDH-wild-type patients; 16/22 were higher at "
        "recurrence (mean change 0.273, 95% CI 0.094–0.452; adjusted sign-flip FDR 0.00721). **B,** Changes for "
        "seven shared leading-edge genes; six were measurable and *CCL4* was absent. **C,** Visium maps of "
        "MG-Inflam20, MDSC-like, and MES-like scores in two untreated primary IDH-wild-type tumors. "
        "Both cases, GBM030 and GBM049, "
        "were also represented in the single-cell atlas used for Figure 2. **D,** Myeloid-adjusted partial Spearman "
        "correlations were positive for MDSC-like (*r*<sub>s</sub> 0.399/0.355) and MES-like scores (0.280/0.270). "
        "Panels C,D provide patient-matched, cross-modal context rather than independent patient replication and are "
        "descriptive because spatial dependence was not modeled. State identity and multicellular niche context are "
        "distinct analytical levels; these maps do not assign MDSC identity to cells expressing the "
        "MG-Inflam20 program.\n"
    )
    legend_path.write_text(legend_text, encoding="utf-8")

    # Retain the build summary and report the panel count from the manifest records.
    result_path = WRITE_ROOT / "FINAL_RESULT.md"
    result_text = (
        "# Figure 3 visual polish\n\n"
        f"- A. Paired GeoMx measurements for {len(patient_table)} patients.\n"
        "- B. Measurable genes are shown. The source table retains every prespecified gene and its measurement "
        "status.\n"
        "- C. Both H&E backgrounds use the same grayscale conversion, with a common scale per score.\n"
        "- D. Partial correlations shown as effect sizes with a zero baseline and a 0–0.5 horizontal axis.\n"
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
