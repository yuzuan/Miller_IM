#!/usr/bin/env python3
"""Figure 3视觉精修：删去缺失空行、统一灰度组织底图、压缩空间效应条。"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts/53_fig3_visual_variety_rebuild.py"
WRITE_ROOT = ROOT / "write/55_figure3_visual_polish/Figure3"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/55_figure3_visual_polish/Figure3"


def load_base():
    spec = importlib.util.spec_from_file_location("step53_figure3", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法载入Step53脚本: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.WRITE_ROOT = WRITE_ROOT
    module.SOURCE_OUT = SOURCE_OUT
    module.FIG_OUT = FIG_OUT
    return module


base = load_base()


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    metadata = {
        "Creator": "Step55 recurrent GBM Figure 3 visual polish",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step55 recurrent GBM Figure 3 visual polish"},
    )
    base.plt.close(fig)
    return pdf, png


base.save_panel = save_panel


def fixed_grayscale_tissue(rgb: np.ndarray) -> np.ndarray:
    """对所有切片应用同一固定灰度转换，不做逐切片强度归一化。"""
    luminance = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
    display = 0.68 + 0.32 * np.clip(luminance, 0.0, 1.0)
    return np.repeat(display[:, :, None], 3, axis=2)


base.soften_tissue = fixed_grayscale_tissue


def panel_b(patients: pd.DataFrame) -> dict[str, str]:
    stem = "Figure3B_shared_gene_patient_fingerprint"
    genes = base.shared_genes()
    table = base.strict_gene_delta(patients, genes)
    if len(table) != 176 or int(table["measured"].sum()) != 154:
        raise ValueError("B面板必须保留8个预定义基因的完整审计源表，且7个基因可测")
    if table.loc[~table["measured"], "gene"].drop_duplicates().tolist() != ["CCL4"]:
        raise ValueError("B面板唯一缺失基因必须为CCL4")
    source = SOURCE_OUT / f"{stem}_source.csv"
    table.to_csv(source, index=False)

    measured_genes = table.loc[table["measured"], "gene"].drop_duplicates().tolist()
    if measured_genes != genes[:-1]:
        raise ValueError(f"可测基因顺序改变: {measured_genes}")
    matrix = (
        table.loc[table["measured"]]
        .pivot(index="gene", columns="patient_label", values="gene_delta")
        .reindex(index=measured_genes, columns=patients["patient_label"].tolist())
    )
    color_limit = float(table["color_limit"].iloc[0])
    clipped = int((table.loc[table["measured"], "gene_delta"].abs() > color_limit).sum())
    if clipped != 4:
        raise ValueError(f"B面板98%色限截色数量改变: {clipped}")

    fig = base.new_figure(116, 86)
    ax = fig.add_axes([0.12, 0.29, 0.66, 0.56])
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        cmap=mpl.colormaps["RdBu_r"],
        vmin=-color_limit,
        vmax=color_limit,
        aspect="auto",
    )
    ax.set_yticks(np.arange(len(measured_genes)), measured_genes, fontsize=6.6)
    ax.set_xticks(
        np.arange(len(patients)),
        patients["patient_label"],
        rotation=90,
        fontsize=4.9,
    )
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(patients), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(measured_genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    color_ax = fig.add_axes([0.19, 0.16, 0.49, 0.023])
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.8, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
    colorbar.outline.set_visible(False)

    summary = table.loc[table["measured"]].drop_duplicates("gene").set_index("gene").reindex(measured_genes)
    sx = fig.add_axes([0.81, 0.29, 0.16, 0.56])
    sx.axvline(0, color="#A0A0A0", linewidth=0.65)
    for index, gene in enumerate(measured_genes):
        row = summary.loc[gene]
        significant = float(row["gene_fdr"]) < 0.05
        sx.hlines(index, float(row["ci_low"]), float(row["ci_high"]), color=base.RECURRENT, linewidth=1.0)
        sx.scatter(
            float(row["mean_gene_delta"]),
            index,
            s=22,
            facecolor=base.RECURRENT if significant else "white",
            edgecolor=base.RECURRENT,
            linewidth=0.8,
            zorder=3,
        )
    sx.set_xlim(-0.45, 0.92)
    sx.set_ylim(len(measured_genes) - 0.5, -0.5)
    sx.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.9)
    sx.set_yticks([])
    sx.set_xlabel("Mean Δ", fontsize=5.6, labelpad=2)
    sx.grid(axis="x", color=base.GRID, linewidth=0.45)
    base.clean_axis(sx, keep_left=False, keep_bottom=True)

    fig.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=9.0, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return base.record("B", stem, "Seven measurable genes shown from eight prespecified shared leading-edge genes", source, pdf, png)


def panel_d() -> dict[str, str]:
    stem = "Figure3D_partial_rho_effect_strip"
    original = pd.read_csv(base.STEP41 / "fig3e_partial_rho_main.csv")
    table = original[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
    table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
    if len(table) != 4 or not (table["partial_rho"] > 0).all():
        raise ValueError("D面板必须包含两切片×MDSC/MES四个正向效应")
    source = SOURCE_OUT / f"{stem}_source.csv"
    table.to_csv(source, index=False)

    fig = base.new_figure(58, 68)
    ax = fig.add_axes([0.31, 0.25, 0.64, 0.55])
    target_y = {"MDSC": 1.0, "MES": 0.0}
    sample_offset = {"GBM030": 0.105, "GBM049": -0.105}
    sample_marker = {"GBM030": "o", "GBM049": "s"}
    target_color = {"MDSC": base.MDSC, "MES": base.MES}

    for target in ["MDSC", "MES"]:
        current = table.loc[table["target"].eq(target)].set_index("sample")
        for sample in ["GBM030", "GBM049"]:
            xx = float(current.loc[sample, "partial_rho"])
            yy = target_y[target] + sample_offset[sample]
            ax.hlines(yy, 0, xx, color=target_color[target], linewidth=1.15, alpha=0.55, zorder=1)
            ax.scatter(
                xx,
                yy,
                s=34,
                marker=sample_marker[sample],
                facecolor=target_color[target],
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            ax.text(xx + 0.014, yy, f"{xx:.3f}", fontsize=5.9, va="center", color=base.TEXT)

    ax.axvline(0, color="#8E8E8E", linewidth=0.75)
    ax.set_xlim(0, 0.50)
    ax.set_ylim(-0.46, 1.46)
    ax.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.5)
    ax.set_xticks([0, 0.2, 0.4])
    ax.set_xlabel(r"Partial Spearman $r_s$", fontsize=6.2)
    ax.grid(axis="x", color=base.GRID, linewidth=0.5)
    base.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="none", color=base.TEXT, markerfacecolor="white", label="GBM030"),
        mpl.lines.Line2D([], [], marker="s", linestyle="none", color=base.TEXT, markerfacecolor="white", label="GBM049"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=5.5,
        loc="upper center",
        bbox_to_anchor=(0.52, -0.30),
        ncol=2,
        handletextpad=0.3,
        columnspacing=0.8,
    )
    fig.text(0.12, 0.955, "Miller-IM score vs niche", fontsize=9.2, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return base.record("D", stem, "Compact descriptive partial-rho effect strip", source, pdf, png)


def image_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            try:
                return ImageFont.truetype(candidate, size=size)
            except OSError:
                pass
    return ImageFont.load_default()


def paste_contained(canvas: Image.Image, source: Path, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        fitted = ImageOps.contain(white.convert("RGB"), (right - left, bottom - top), Image.Resampling.LANCZOS)
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 115, "top": 270, "right": 1885, "bottom": 2110},
        {"panel": "B", "left": 1975, "top": 270, "right": 4205, "bottom": 2110},
        {"panel": "C", "left": 115, "top": 2450, "right": 2465, "bottom": 4040},
        {"panel": "D", "left": 2610, "top": 2450, "right": 4205, "bottom": 4235},
    ]
    canvas = Image.new("RGB", (4320, 4370), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = image_font(78, bold=True)
    section_font = image_font(62, bold=True)

    drawer.text((115, 55), "Independent recurrence validation", fill=base.RECURRENT, font=section_font)
    drawer.line((115, 145, 4205, 145), fill="#E8B2AA", width=5)
    drawer.text((115, 2190), "Miller-IM program in MDSC/MES-rich niches", fill=base.MDSC, font=section_font)
    drawer.line((115, 2285, 4205, 2285), fill="#A8D9CF", width=5)

    for item in layout:
        panel = item["panel"]
        paste_contained(canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"]))
        drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

    png = FIG_OUT / "Figure3_visual_polish_preview.png"
    gray = FIG_OUT / "Figure3_visual_polish_preview_grayscale.png"
    pdf = FIG_OUT / "Figure3_visual_polish_preview.pdf"
    source = SOURCE_OUT / "Figure3_visual_polish_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(600, 600), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(pdf, format="PDF", resolution=600.0, creationDate="", modDate="")
    pd.DataFrame(layout).to_csv(source, index=False)
    return pdf, png, gray, source


def write_legend() -> Path:
    content = """# Figure 3 legend

**Figure 3 | Paired myeloid-enriched tissue support and patient-matched spatial context of the Miller-IM program.**

**A,** Paired IBA1-positive GeoMx Miller-IM scores in 22 IDH-wild-type patients; 16/22 were higher at recurrence (mean change 0.273, 95% CI 0.094–0.452; adjusted sign-flip FDR 0.00721). **B,** Changes for eight shared leading-edge genes; seven were measurable and *CCL4* was absent. **C,** Visium maps of Miller-IM, MDSC-like, and MES-like scores in two untreated primary IDH-wild-type tumors. Both cases, GBM030 and GBM049, were also represented in the single-cell atlas used for Figure 2. **D,** Myeloid-adjusted partial Spearman correlations were positive for MDSC-like (*r*<sub>s</sub> 0.399/0.355) and MES-like scores (0.280/0.270). Panels C,D provide patient-matched, cross-modal context rather than independent patient replication and are descriptive because spatial dependence was not modeled. State identity and multicellular niche context are distinct analytical levels; these maps do not assign MDSC identity to cells expressing the Miller-IM program.
"""
    path = WRITE_ROOT / "Figure3_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest: Path, legend: Path) -> None:
    content = f"""# Figure 3 visual polish

- A：保留22对患者GeoMx paired dumbbell。
- B：图面只展示7个可测基因；8个预定义基因和CCL4缺失信息仍完整保留在源表与图注。
- C：两张真实H&E统一采用同一固定灰度低对比显示，三类score及各自跨切片共享色标不变。
- D：压缩为保留0基线和0–0.5横轴的四效应lollipop strip，不新增亚型、p值、FDR或CI。
- Manifest：`{manifest}`
- Legend：`{legend}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main() -> None:
    base.apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    if not base.SPATIAL_IMAGE_ROOT.exists():
        raise FileNotFoundError(f"空间组织图目录未挂载: {base.SPATIAL_IMAGE_ROOT}")

    patients = base.patient_table()
    records = [base.panel_a(patients), panel_b(patients), base.panel_c(), panel_d()]
    if [row["panel"] for row in records] != list("ABCD"):
        raise ValueError("Figure3面板必须唯一锁定为A-D")
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_path = WRITE_ROOT / "Figure3_panel_manifest.csv"
    manifest = pd.DataFrame(records)
    manifest["preview_pdf"] = str(preview_pdf)
    manifest["preview_png"] = str(preview_png)
    manifest["preview_grayscale"] = str(preview_gray)
    manifest["preview_layout"] = str(preview_layout)
    manifest.to_csv(manifest_path, index=False)
    legend_path = write_legend()
    write_result(manifest_path, legend_path)
    print("STEP55_FIGURE3_POLISH_COMPLETE panels=4 preview=1 mapping=A-D")


if __name__ == "__main__":
    main()
