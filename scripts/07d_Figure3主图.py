#!/usr/bin/env python3
# 目的：沿用Figure3科学口径，保留六可测基因图面、固定灰度H&E和四个空间效应。
# 输入：Step53复用核心函数读取的GeoMx、Step41空间源表与两张Visium组织图。
# 输出：write/55_figure3_visual_polish/Figure3/source_data、Figure3_panel_manifest.csv、
#       Figure3_legend.md、FINAL_RESULT.md；figures同模块的A-D面板及Figure3_visual_polish_preview。
# 运行：在仓库根目录执行python scripts/07d_Figure3主图.py。
# 保留Python及pandas/NumPy/Matplotlib/Pillow，维持现有复现路径、公式和PDF画面。
# 公开路径、字段与数据文件名不变；PNG从PDF渲染，纯预览PIL画布保留。

# 独立面板主名如下，各写同名.pdf/.png及source_data下的主名_source.csv：
#   Figure3A_geomx_paired_dumbbell
#   Figure3B_shared_gene_patient_fingerprint
#   Figure3C_spatial_score_atlas
#   Figure3D_partial_rho_effect_strip
# C源表使用.csv.gz；预览同主名另存_grayscale.png和source_data下_layout.csv。

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

# 动态导入只加载Step53可复用定义；基础文件的主流程受__main__保护。
root = Path(__file__).resolve().parents[1]
base_script = root / "scripts/07c_空间与配对表达图.py"
spec = importlib.util.spec_from_file_location("step53_figure3", base_script)
if spec is None or spec.loader is None:
    raise RuntimeError(f"无法载入Step53脚本: {base_script}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

# 原初始化顺序与全局seed保留；输出目录不再通过修改被导入模块的全局量传递。
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
    raise FileNotFoundError(f"空间组织图目录未挂载: {spatial_image_root}")

# A直接复用患者配对图，输出路径和PDF来源文字作为显式参数传入。
patients_table = base.patient_table()
records = [base.panel_a(patients_table, source_out, figure_out, "Step55 recurrent GBM Figure 3 visual polish")]

# B源表仍保留七个预定义基因，仅图面显示六个可测基因；不改变排序、FDR或色限。

genes = base.shared_genes()
gene_delta_table = base.strict_gene_delta(patients_table, genes)
if len(gene_delta_table) != 154 or int(gene_delta_table["measured"].sum()) != 132:
    raise ValueError("B面板必须保留7个预定义基因的完整审计源表，且6个基因可测")
if gene_delta_table.loc[~gene_delta_table["measured"], "gene"].drop_duplicates().tolist() != ["CCL4"]:
    raise ValueError("B面板唯一缺失基因必须为CCL4")
stem = "Figure3B_shared_gene_patient_fingerprint"

# 写出与原入口同字段、同顺序的面板源表。
source_path = source_out / f"{stem}_source.csv"
gene_delta_table.to_csv(source_path, index=False)


# 原七基因源表保留，仅选六个可测基因进入画面；98%截色数量仍核对为3。
measured_genes = gene_delta_table.loc[gene_delta_table["measured"], "gene"].drop_duplicates().tolist()
if measured_genes != genes[:-1]:
    raise ValueError(f"可测基因顺序改变: {measured_genes}")
delta_pivot_table = (
    gene_delta_table.loc[gene_delta_table["measured"]]
    .pivot(index="gene", columns="patient_label", values="gene_delta")
    .reindex(index=measured_genes, columns=patients_table["patient_label"].tolist())
)
color_limit = float(gene_delta_table["color_limit"].iloc[0])
clipped = int((gene_delta_table.loc[gene_delta_table["measured"], "gene_delta"].abs() > color_limit).sum())
if clipped != 3:
    raise ValueError(f"B面板98%色限截色数量改变: {clipped}")

# 使用原毫米尺寸和轴位置构图，确保PDF布局保持原样。
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

# 色条仅解释原显示色限，不参与统计量计算。
color_axis = figure.add_axes([0.19, 0.16, 0.49, 0.023])
colorbar = figure.colorbar(image, cax=color_axis, orientation="horizontal")
colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.8, labelpad=1.5)
colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
colorbar.outline.set_visible(False)

# 右侧汇总按可测基因顺序排列，统计值仍来自原源表。
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

# 汇总图保留原零线、范围和标记含义。
summary_axis.set_xlim(-0.45, 0.92)
summary_axis.set_ylim(len(measured_genes) - 0.5, -0.5)
summary_axis.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.9)
summary_axis.set_yticks([])
summary_axis.set_xlabel("Mean Δ", fontsize=5.6, labelpad=2)
summary_axis.grid(axis="x", color="#E7E7E7", linewidth=0.45)
clean_axis(summary_axis, keep_left=False, keep_bottom=True)

figure.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=9.0, fontweight="bold", va="top")

# 保存原PDF图面，再由同一PDF生成600 dpi PNG。
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


# C在共享核心函数中选择原Step55固定灰度公式，三项score和色标计算均未改变。
records.append(base.panel_c(source_out, figure_out, "Step55 recurrent GBM Figure 3 visual polish", True))

# D只读取四个已有partial-rho，保留原0基线及0至0.5横轴。

partial_path = root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3e_partial_rho_main.csv"
partial_source_table = pd.read_csv(partial_path)
partial_table = partial_source_table[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
partial_table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
if len(partial_table) != 4 or not (partial_table["partial_rho"] > 0).all():
    raise ValueError("D面板必须包含两切片×MDSC/MES四个正向效应")
stem = "Figure3D_partial_rho_effect_strip"

# 写出与原入口同字段、同顺序的面板源表。
source_path = source_out / f"{stem}_source.csv"
partial_table.to_csv(source_path, index=False)

# 使用原毫米尺寸和轴位置构图，确保PDF布局保持原样。
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

# 横轴及标签沿用原图定义，不改变比较方向。
axis.set_xlim(0, 0.50)

# 坐标范围与原图一致，不裁剪或重排分析数据。
axis.set_ylim(-0.46, 1.46)
axis.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.5)
axis.set_xticks([0, 0.2, 0.4])
axis.set_xlabel(r"Partial Spearman $r_s$", fontsize=6.2)
axis.grid(axis="x", color="#E7E7E7", linewidth=0.5)
clean_axis(axis, keep_left=False, keep_bottom=True)
axis.tick_params(axis="y", length=0, pad=3)

# 图例标记与实际面板一致，保留原颜色和顺序。
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

# 保存原PDF图面，再由同一PDF生成600 dpi PNG。
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step55 recurrent GBM Figure 3 visual polish")
records.append(panel_record("D", stem, "Compact descriptive partial-rho effect strip", source_path, pdf_path, png_path))


# 组合预览沿用原4320×4370画布及位置；不调整独立面板比例或数据图形。

canvas = Image.new("RGB", (4320, 4370), "white")
drawer = ImageDraw.Draw(canvas)
section_font = image_font(62, bold=True)

drawer.text((115, 55), "Independent recurrence support", fill="#E64B35", font=section_font)
drawer.line((115, 145, 4205, 145), fill="#E8B2AA", width=5)
drawer.text((115, 2190), "Miller-IM program in MDSC/MES-rich niches", fill="#159D88", font=section_font)
drawer.line((115, 2285, 4205, 2285), fill="#A8D9CF", width=5)

if [row["panel"] for row in records] != list("ABCD"):
    raise ValueError("Figure3面板必须唯一锁定为A-D")

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

# PIL预览保留原像素画布和分辨率，供组合排版查看。
png_path = figure_out / "Figure3_visual_polish_preview.png"
canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
preview_gray_path = figure_out / "Figure3_visual_polish_preview_grayscale.png"
ImageOps.grayscale(canvas).convert("RGB").save(preview_gray_path, format="PNG", dpi=(600, 600), optimize=True)
pdf_path = figure_out / "Figure3_visual_polish_preview.pdf"
canvas.save(pdf_path, format="PDF", resolution=600.0, creationDate="", modDate="")
layout_table = pd.DataFrame(layout)
source_path = source_out / "Figure3_visual_polish_preview_layout.csv"
layout_table.to_csv(source_path, index=False)


# 一行一个面板，保留原manifest字段及四种预览引用。
manifest_table = pd.DataFrame(records)
manifest_table["preview_pdf"] = str(pdf_path)
manifest_table["preview_png"] = str(png_path)
manifest_table["preview_grayscale"] = str(preview_gray_path)
manifest_table["preview_layout"] = str(source_path)
manifest_path = write_root / "Figure3_panel_manifest.csv"
manifest_table.to_csv(manifest_path, index=False)

# 图注与说明文字逐字保留，仅改变源码断行，不修改稿件。

content = (
    "# Figure 3 legend\n\n**Figure 3 | Paired myeloid-enriched tissue support and patient-matched "
    "spatial context of the Miller-IM program.**\n\n**A,** Paired IBA1-positive GeoMx Miller-IM scores"
    " in 22 IDH-wild-type patients; 16/22 were higher at recurrence (mean change 0.273, 95% CI "
    "0.094–0.452; adjusted sign-flip FDR 0.00721). **B,** Changes for seven shared leading-edge "
    "genes; six were measurable and *CCL4* was absent. **C,** Visium maps of Miller-IM, MDSC-like, "
    "and MES-like scores in two untreated primary IDH-wild-type tumors. Both cases, GBM030 and "
    "GBM049, were also represented in the single-cell atlas used for Figure 2. **D,** "
    "Myeloid-adjusted partial Spearman correlations were positive for MDSC-like (*r*<sub>s</sub> "
    "0.399/0.355) and MES-like scores (0.280/0.270). Panels C,D provide patient-matched, cross-modal"
    " context rather than independent patient replication and are descriptive because spatial "
    "dependence was not modeled. State identity and multicellular niche context are distinct "
    "analytical levels; these maps do not assign MDSC identity to cells expressing the Miller-IM "
    "program.\n"
)
legend_path = write_root / "Figure3_legend.md"
legend_path.write_text(content, encoding="utf-8")


content = (
    "# Figure 3 visual polish\n\n- A：保留22对患者GeoMx paired dumbbell。\n- "
    "B：图面只展示6个可测基因；7个预定义基因和CCL4缺失信息仍完整保留在源表与图注。\n- C：两张真实H&E统一采用同一固定灰度低对比显示，三类score及各自跨切片共享色标不变。\n- "
    "D：压缩为保留0基线和0–0.5横轴的四效应lollipop strip，不新增亚型、p值、FDR或CI。\n- Manifest：`"
    f"{manifest_path}"
    "`\n- Legend：`"
    f"{legend_path}"
    "`\n"
)
(write_root / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


# 核对实际文件、预览尺寸与清单行列数；缺失或空产物直接停止。
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

    # 存在性与体积检查针对本次输出，不接受空文件。
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {output_path.name}")
with Image.open(png_path) as preview_image:
    if preview_image.size != canvas.size:
        raise ValueError("Preview dimensions differ from the original canvas")
    preview_image.verify()

# 回读实际清单，核对行列数与本次面板记录一致。
manifest_check_table = pd.read_csv(manifest_path)
if manifest_check_table.shape != manifest_table.shape:
    raise ValueError("Panel manifest dimensions changed during writing")
print(f"Figure3 panels={len(records)} preview=1")
