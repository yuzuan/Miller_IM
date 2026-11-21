#!/usr/bin/env python3
# 目的：使用前序空间源表及组织图像重排Figure3与SupplementaryFigure5。
# 输入：write/39_independent_miller_mg_inflammatory_figures内Figure3和SupplementaryFigure6源表；组织图像和缩放因子。
# 输出：write/41_mg_inflammatory_sci_rebuild内Figure3/SupplementaryFigure5源表及面板清单；
# figures/41_mg_inflammatory_sci_rebuild内对应panel_library的PDF与PNG。
# 运行：在项目根使用既有sc环境执行 python scripts/07b_空间关联图.py。
# 保留Python：Scanpy、AnnData及Matplotlib是已公开复现路径；改成R超出本轮纯风格整理。
# 存量source_data/config、数据文件名、输出字段、write及figures路径保留原公开合同。

from __future__ import annotations
import gzip
import io
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
colors = figure_style.colors
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
new_figure = figure_style.new_figure
save_figure = figure_style.save_figure

# 保留原全局随机种子及设置位置，不改变后续随机序列。
np.random.seed(20260713)
apply_publication_style()

# 创建Figure3与SupplementaryFigure5的源表和panel_library输出目录。
root = Path(__file__).resolve().parents[1]

fig3_out = root / "figures/41_mg_inflammatory_sci_rebuild/Figure3/panel_library"

write_fig3 = root / "write/41_mg_inflammatory_sci_rebuild/Figure3"

fig3_source_out = write_fig3 / "source_data"

s5_out = root / "figures/41_mg_inflammatory_sci_rebuild/SupplementaryFigure5/panel_library"

write_s5 = root / "write/41_mg_inflammatory_sci_rebuild/SupplementaryFigure5"

s5_source_out = write_s5 / "source_data"

for directories_path in [fig3_source_out, s5_source_out, fig3_out, s5_out]:
    directories_path.mkdir(parents=True, exist_ok=True)

# 读取前序患者配对表和入口统计，只保留strict_pass_idhwt用于主配对面板。
old_fig3_source = root / "write/39_independent_miller_mg_inflammatory_figures/Figure3/source_csv"

paired_data_table = pd.read_csv(old_fig3_source / "geomx_raw20_paired_deltas.csv")
paired_summary_table = pd.read_csv(old_fig3_source / "geomx_raw20_entry_level_summary.csv").set_index("entry")
paired_data_table = paired_data_table.loc[paired_data_table["entry"].eq("strict_pass_idhwt")].copy()

def mean_ci(values: pd.Series | np.ndarray) -> tuple[float, float, float]:
    """输入差值向量；返回原均值和1.96倍标准误区间，空输入返回NaN，单值区间等于该值。"""

    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        result = (np.nan, np.nan, np.nan)
        return result
    mean = float(array.mean())
    if len(array) == 1:
        result = (mean, mean, mean)
        return result
    sem = float(array.std(ddof=1) / np.sqrt(len(array)))

    delta = 1.96 * sem
    result = (mean, mean - delta, mean + delta)
    return result

# 主配对源表补入患者差值均值、95%区间、人数、上升人数及原入口FDR。
paired_mean_delta, paired_ci_low, paired_ci_high = mean_ci(paired_data_table["delta"])
paired_data_table["mean_delta"] = paired_mean_delta
paired_data_table["ci_low"] = paired_ci_low
paired_data_table["ci_high"] = paired_ci_high
paired_data_table["n_pairs"] = len(paired_data_table)

paired_data_table["n_up"] = int((paired_data_table["delta"] > 0).sum())
paired_data_table["fdr"] = float(paired_summary_table.loc["strict_pass_idhwt", "fdr"])
paired_out_path = fig3_source_out / "fig3a_geomx22_paired_change.csv"

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = paired_out_path
paired_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
fig3a_table = paired_data_table

# 读取前序基因效应，按均值差降序和基因名升序写出森林图源表。
genes_data_table = pd.read_csv(old_fig3_source / "fig3c_geomx_gene_forest.csv")
genes_data_table = genes_data_table.sort_values(["mean_log2_delta", "gene"], ascending=[False, True]).reset_index(
    drop = True
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_source_out / "fig3b_geomx18_gene_forest.csv"
genes_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
fig3b_table = genes_data_table

# 按GBM030、GBM049顺序读取空间点并补充组织图缩放信息。
spatial_all_spots_table_total = pd.read_csv(old_fig3_source / "gse276841_raw20_spot_scores.csv.gz")
spatial_spatial_tables: dict[str, pd.DataFrame] = {}

data_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()

spatial_image_root = data_root / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"

spatial_prefix = {"GBM030": "GSM8506649_03422", "GBM049": "GSM8506650_25526"}

for spatial_sample in ["GBM030", "GBM049"]:

    def load_grayscale_image(sample: str) -> tuple[np.ndarray, float, float]:
        """输入原切片名称；返回灰度图、组织图缩放因子及缩放后的点直径；缺输入由读盘报错。"""

        # 读取gzip组织图并转灰度，再读取缩放因子以计算高分辨率点直径。
        prefix = spatial_prefix[sample]
        with gzip.open(spatial_image_root / f"{prefix}_tissue_hires_image.png.gz", "rb") as handle:
            image = Image.open(io.BytesIO(handle.read())).convert("L")
        with gzip.open(spatial_image_root / f"{prefix}_scalefactors_json.json.gz", "rt") as handle:
            scale_meta = json.load(handle)
        gray = np.asarray(image, dtype=np.float32)
        scale = float(scale_meta["tissue_hires_scalef"])
        spot_diameter = float(scale_meta["spot_diameter_fullres"]) * scale
        result = (gray, scale, spot_diameter)
        return result

    spatial_gray, spatial_scale, spatial_spot_diameter = load_grayscale_image(spatial_sample)

    # 按切片提取空间点，并用原缩放因子转换为组织图像像素坐标。
    spatial_current_table = spatial_all_spots_table_total.loc[
        spatial_all_spots_table_total["sample"].eq(spatial_sample)
    ].copy()
    spatial_current_table["x_hires"] = spatial_current_table["pxl_col_in_fullres"] * spatial_scale
    spatial_current_table["y_hires"] = spatial_current_table["pxl_row_in_fullres"] * spatial_scale
    spatial_current_table["image_width_px"] = spatial_gray.shape[1]
    spatial_current_table["image_height_px"] = spatial_gray.shape[0]
    spatial_current_table["tissue_hires_scalef"] = spatial_scale
    spatial_current_table["spot_diameter_hires_px"] = spatial_spot_diameter

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = (
        fig3_source_out
        / f"{('fig3c' if spatial_sample == 'GBM030' else 'fig3d')}_{spatial_sample.lower()}_spatial_raw20_mdsc_mes.csv"
    )
    spatial_current_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
    spatial_spatial_tables[spatial_sample] = spatial_current_table

spatial_tables = spatial_spatial_tables

# 主相关面板保留MDSC和MES，Fisher区间使用既有相关系数及有效点数。
rho_data_table = pd.read_csv(old_fig3_source / "fig3f_partial_rho.csv")
rho_data_table = rho_data_table.loc[rho_data_table["target"].isin(["MDSC", "MES"])].copy()

def fisher_ci(rho: float, n_obs: int) -> tuple[float, float]:
    """输入相关系数和点数；返回原Fisher变换95%区间，无效边界返回NaN。"""

    if not np.isfinite(rho) or n_obs <= 3 or abs(rho) >= 1:
        result = (np.nan, np.nan)
        return result
    z_value = np.arctanh(rho)
    sem = 1 / np.sqrt(n_obs - 3)
    low, high = z_value - 1.96 * sem, z_value + 1.96 * sem
    result = (float(np.tanh(low)), float(np.tanh(high)))
    return result

# 按输入行顺序计算区间并写回原索引，避免切片/目标对应关系改变。
ci_rows = []

for row in rho_data_table.itertuples(index=False):
    ci_low, ci_high = fisher_ci(float(row.partial_rho), int(row.n_spots))
    ci_rows.append((ci_low, ci_high))
rho_ci_bounds_table = pd.DataFrame(ci_rows, index=rho_data_table.index)
rho_data_table["ci_low"] = rho_ci_bounds_table[0].to_numpy()
rho_data_table["ci_high"] = rho_ci_bounds_table[1].to_numpy()
rho_data_table["panel_label"] = rho_data_table["sample_label"] + " · " + rho_data_table["target"]

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_source_out / "fig3e_partial_rho_main.csv"
rho_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
fig3e_table = rho_data_table

# 敏感性源表复用前序配对差值、入口FDR及各分层FDR。
sensitivity_paired_table = pd.read_csv(old_fig3_source / "geomx_raw20_paired_deltas.csv")
sensitivity_entry_summary_table = pd.read_csv(old_fig3_source / "geomx_raw20_entry_level_summary.csv").set_index(
    "entry"
)

old_s6_source = root / "write/39_independent_miller_mg_inflammatory_figures/SupplementaryFigure6/source_csv"

sensitivity_subgroup_summary_table = pd.read_csv(old_s6_source / "geomx_raw20_subgroup_summary.csv")
sensitivity_subgroup_summary_table = sensitivity_subgroup_summary_table.set_index(["entry", "subgroup"])

sensitivity_rows: list[dict[str, object]] = []

def add_row(
    section: str,
    sensitivity_order: int,
    sensitivity_label: str,
    deltas: pd.Series,
    fdr: float,
    n_up: int | None = None,
) -> None:
    """输入分层名称、患者差值、FDR及可选上升人数；向既有敏感性记录列表追加一行；无返回值。"""

    sensitivity_effect, sensitivity_ci_low, sensitivity_ci_high = mean_ci(deltas.to_numpy(dtype=float))
    sensitivity_rows.append(
        {
            "section": section,
            "section_order": sensitivity_order,
            "label": sensitivity_label,
            "effect": sensitivity_effect,
            "ci_low": sensitivity_ci_low,
            "ci_high": sensitivity_ci_high,
            "n_pairs": int(deltas.notna().sum()),
            "n_up": int((deltas > 0).sum()) if n_up is None else int(n_up),
            "fdr": float(fdr),
        }
    )

# 先写严格IDH-wt主分析，再依次加入入口、IDH、治疗和交互分层。
sensitivity_main = sensitivity_paired_table.loc[sensitivity_paired_table["entry"].eq("strict_pass_idhwt"), "delta"]
add_row(
    "Main analysis",
    1,
    "IBA1+ strict + IDH-wt",
    sensitivity_main,
    float(sensitivity_entry_summary_table.loc["strict_pass_idhwt", "fdr"]),
    int(sensitivity_entry_summary_table.loc["strict_pass_idhwt", "n_positive"]),
)

# 入口敏感性分别放宽AOI筛选和纳入全部IDH，沿用各入口的配对差值。
for sensitivity_order, sensitivity_entry, sensitivity_label in [
    (2, "all_aoi_idhwt", "IBA1+ all AOI + IDH-wt"),
    (2, "strict_pass_all_idh", "IBA1+ strict + all IDH"),
]:
    sensitivity_subset = sensitivity_paired_table.loc[sensitivity_paired_table["entry"].eq(sensitivity_entry), "delta"]
    add_row(
        section           = "Entry sensitivity",
        sensitivity_order = sensitivity_order,
        sensitivity_label = sensitivity_label,
        deltas            = sensitivity_subset,
        fdr               = float(sensitivity_entry_summary_table.loc[sensitivity_entry, "fdr"]),
    )

# IDH敏感性在strict_pass_all_idh入口内分别计算野生型和突变型。
for sensitivity_subgroup, sensitivity_label in [
    ("IDH_WT", "All-IDH entry: IDH-wt"),
    ("IDH_mutated", "All-IDH entry: IDH-mutant"),
]:
    sensitivity_subset = sensitivity_paired_table.loc[
        sensitivity_paired_table["entry"].eq("strict_pass_all_idh")
        & sensitivity_paired_table["IDH_status"].eq("IDH_WT" if sensitivity_subgroup == "IDH_WT" else "IDH_mutated"),
        "delta",
    ]
    add_row(
        section = "IDH sensitivity",
        sensitivity_order = 3,
        sensitivity_label = sensitivity_label,
        deltas = sensitivity_subset,
        fdr = float(sensitivity_subgroup_summary_table.loc[("strict_pass_all_idh", sensitivity_subgroup), "fdr"]),
    )

# 治疗分层仅使用strict_pass_idhwt，分别汇总Control和Nivolumab。
for sensitivity_subgroup, sensitivity_label in [("Control", "Control"), ("Nivolumab", "Nivolumab")]:
    sensitivity_subset = sensitivity_paired_table.loc[
        sensitivity_paired_table["entry"].eq("strict_pass_idhwt")
        & sensitivity_paired_table["trial_setting"].eq(sensitivity_subgroup),
        "delta",
    ]
    add_row(
        section = "Treatment strata",
        sensitivity_order = 4,
        sensitivity_label = sensitivity_label,
        deltas = sensitivity_subset,
        fdr = float(sensitivity_subgroup_summary_table.loc[("strict_pass_idhwt", sensitivity_subgroup), "fdr"]),
    )

# 交互效应使用两治疗组配对差值的均值差，区间沿用独立组方差公式。
sensitivity_control = sensitivity_paired_table.loc[
    sensitivity_paired_table["entry"].eq("strict_pass_idhwt") & sensitivity_paired_table["trial_setting"].eq("Control"),
    "delta",
]
sensitivity_nivolumab = sensitivity_paired_table.loc[
    sensitivity_paired_table["entry"].eq("strict_pass_idhwt")
    & sensitivity_paired_table["trial_setting"].eq("Nivolumab"),
    "delta",
]

sensitivity_control_array = sensitivity_control.to_numpy(dtype=float)
sensitivity_nivolumab_array = sensitivity_nivolumab.to_numpy(dtype=float)
sensitivity_effect = float(sensitivity_nivolumab_array.mean() - sensitivity_control_array.mean())
if len(sensitivity_control_array) <= 1 or len(sensitivity_nivolumab_array) <= 1:
    sensitivity_ci_low = np.nan
    sensitivity_ci_high = np.nan
else:
    sensitivity_se = np.sqrt(
        sensitivity_control_array.var(ddof=1) / len(sensitivity_control_array)
        + sensitivity_nivolumab_array.var(ddof=1) / len(sensitivity_nivolumab_array)
    )
    sensitivity_ci_delta = 1.96 * sensitivity_se
    sensitivity_ci_low = float(sensitivity_effect - sensitivity_ci_delta)
    sensitivity_ci_high = float(sensitivity_effect + sensitivity_ci_delta)

sensitivity_rows.append(
    {
        "section": "Interaction",
        "section_order": 5,
        "label": "Nivolumab minus control",
        "effect": sensitivity_effect,
        "ci_low": sensitivity_ci_low,
        "ci_high": sensitivity_ci_high,
        "n_pairs": int(len(sensitivity_control) + len(sensitivity_nivolumab)),
        "n_up": np.nan,
        "fdr": float(
            sensitivity_subgroup_summary_table.loc[("strict_pass_idhwt", "Interaction_Nivolumab_minus_Control"), "fdr"]
        ),
    }
)

# 保留分层追加顺序作为森林图row_order，写出原敏感性源表。
sensitivity_forest_table = pd.DataFrame(sensitivity_rows)
sensitivity_forest_table["row_order"] = range(len(sensitivity_forest_table))

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s5_source_out / "s5a_geomx_sensitivity_forest.csv"
sensitivity_forest_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
s5a_table = sensitivity_forest_table

# 对四种空间目标逐切片补充Fisher区间；此处保留前序完整目标集合。
slice_data_table = pd.read_csv(old_s6_source / "gse276841_raw20_slice_sensitivity.csv")
ci_rows = []
for row in slice_data_table.itertuples(index=False):
    ci_low, ci_high = fisher_ci(float(row.partial_rho), int(row.n_spots))

    ci_rows.append((ci_low, ci_high))
slice_ci_bounds_table = pd.DataFrame(ci_rows, index=slice_data_table.index)
slice_data_table["ci_low"] = slice_ci_bounds_table[0].to_numpy()
slice_data_table["ci_high"] = slice_ci_bounds_table[1].to_numpy()
slice_data_table["panel_label"] = slice_data_table["sample"] + " · " + slice_data_table["target"]

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s5_source_out / "s5b_slice_sensitivity_forest.csv"
slice_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
s5b_table = slice_data_table

# 88毫米配对面板用灰线连接患者，原发和复发端点分别着色。
paired_plot_data_table = fig3a_table

paired_plot_fig = new_figure(88, 88)
paired_plot_ax = paired_plot_fig.add_axes([0.16, 0.15, 0.78, 0.70])

paired_plot_x_positions = np.array([0.0, 1.0])
paired_plot_primary = paired_plot_data_table["Primary"].to_numpy(dtype=float)
paired_plot_recurrent = paired_plot_data_table["Recurrence"].to_numpy(dtype=float)

for paired_plot_left, paired_plot_right in zip(paired_plot_primary, paired_plot_recurrent):
    paired_plot_ax.plot(
        paired_plot_x_positions, [paired_plot_left, paired_plot_right], color="#C4C4C4", linewidth=0.8, zorder=1
    )
paired_plot_ax.scatter(
    np.full(len(paired_plot_primary), paired_plot_x_positions[0]),
    paired_plot_primary,
    s          = 16,
    color      = colors["primary"],
    linewidths = 0,
    zorder     = 3,
)

paired_plot_ax.scatter(
    np.full(len(paired_plot_recurrent), paired_plot_x_positions[1]),
    paired_plot_recurrent,
    s          = 16,
    color      = colors["recurrent"],
    linewidths = 0,
    zorder     = 3,
)

# 标题与底部区间文字读取源表统计值，人数和上升数由当前患者差值计算。
paired_plot_mean_delta = float(paired_plot_data_table["mean_delta"].iloc[0])

paired_plot_ci_low = float(paired_plot_data_table["ci_low"].iloc[0])
paired_plot_ci_high = float(paired_plot_data_table["ci_high"].iloc[0])
paired_plot_n_pairs = int(paired_plot_data_table["n_pairs"].iloc[0])
paired_plot_n_up = int(paired_plot_data_table["n_up"].iloc[0])
paired_plot_fdr = float(paired_plot_data_table["fdr"].iloc[0])

paired_plot_ax.set_xticks(paired_plot_x_positions, ["Primary", "Recurrent"])
paired_plot_ax.set_ylabel("Raw20 score")
paired_plot_ax.set_xlim(-0.3, 1.3)
clean_axis(paired_plot_ax)

paired_plot_ax.grid(axis="y", color=colors["grid"], linewidth=0.6)
paired_plot_fig.text(0.16, 0.96, "GeoMx IBA1+ paired change", ha="left", va="top", fontsize=8.5, weight="bold")
paired_plot_fig.text(
    0.16,
    0.90,
    (
        f"{paired_plot_n_pairs}"
        " IDH-wt pairs | "
        f"{paired_plot_n_up}"
        "/"
        f"{paired_plot_n_pairs}"
        " recurrent-higher | mean delta "
        f"{paired_plot_mean_delta:+.3f}"
    ),
    ha       = "left",
    va       = "top",
    fontsize = 7.0,
    color    = colors["neutral"],
)

paired_plot_fig.text(
    0.16,
    0.05,
    f"95% CI {paired_plot_ci_low:+.3f} to {paired_plot_ci_high:+.3f} | FDR {paired_plot_fdr:.3f}",
    ha       = "left",
    va       = "bottom",
    fontsize = 7.0,
    color    = colors["text"],
)

# 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
save_figure(paired_plot_fig, fig3_out, "fig3a_geomx22_paired_change")

# 基因森林图保留原120×112毫米画布，点与区间按差值正负着色。
genes_plot_data_table = fig3b_table

# 面板标题的18个可测基因为原冻结文案，绘图前核对当前基因数。
if genes_plot_data_table["gene"].nunique() != 18:
    raise RuntimeError("Figure3 frozen 18-gene label no longer matches inputs")

genes_plot_fig = new_figure(120, 112)
genes_plot_ax = genes_plot_fig.add_axes([0.23, 0.10, 0.71, 0.82])
genes_plot_y_values = np.arange(len(genes_plot_data_table))
genes_plot_colors = np.where(
    genes_plot_data_table["mean_log2_delta"].to_numpy(dtype=float) >= 0, colors["recurrent"], colors["primary"]
)
genes_plot_ax.hlines(
    genes_plot_y_values,
    genes_plot_data_table["ci_low"],
    genes_plot_data_table["ci_high"],
    color     = genes_plot_colors,
    linewidth = 1.05,
)

genes_plot_ax.scatter(
    genes_plot_data_table["mean_log2_delta"], genes_plot_y_values, s=20, color=genes_plot_colors, zorder=3
)
genes_plot_ax.axvline(0, color=colors["neutral"], linewidth=0.7)
genes_plot_ax.set_yticks(genes_plot_y_values, genes_plot_data_table["gene"])
genes_plot_ax.invert_yaxis()
genes_plot_ax.set_xlabel("Mean paired log2 delta")
clean_axis(genes_plot_ax)
genes_plot_ax.grid(axis="x", color=colors["grid"], linewidth=0.6)
genes_plot_x_span = float(genes_plot_data_table["ci_high"].max() - genes_plot_data_table["ci_low"].min())

# 根据区间跨度留出右侧人数/FDR文字区，保持原比例边距。
genes_plot_x_text = float(genes_plot_data_table["ci_high"].max()) + genes_plot_x_span * 0.08
genes_plot_ax.set_xlim(
    float(genes_plot_data_table["ci_low"].min()) - genes_plot_x_span * 0.08,
    float(genes_plot_data_table["ci_high"].max()) + genes_plot_x_span * 0.46,
)
for genes_plot_yy, genes_plot_n_pairs, genes_plot_fdr in zip(
    genes_plot_y_values, genes_plot_data_table["n_pairs"], genes_plot_data_table["fdr"]
):
    genes_plot_ax.text(
        genes_plot_x_text,
        genes_plot_yy,
        f"n={int(genes_plot_n_pairs)}  FDR={genes_plot_fdr:.3f}",
        va       = "center",
        ha       = "left",
        fontsize = 6.4,
        color    = colors["text"],
    )

genes_plot_fig.text(0.23, 0.97, "Measurable raw20 genes in GeoMx", ha="left", va="top", fontsize=8.5, weight="bold")
genes_plot_fig.text(
    0.23,
    0.93,
    "18 genes detectable in strict IBA1+ IDH-wt pairs",
    ha       = "left",
    va       = "top",
    fontsize = 7.0,
    color    = colors["neutral"],
)

# 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
save_figure(genes_plot_fig, fig3_out, "fig3b_geomx18_gene_forest")

for sample in ["GBM030", "GBM049"]:
    stem = "fig3c_gbm030_spatial_raw20_mdsc_mes" if sample == "GBM030" else "fig3d_gbm049_spatial_raw20_mdsc_mes"

    spatial_plot_sample = sample
    spatial_plot_data_table = spatial_tables[sample]
    spatial_plot_stem = stem

    spatial_plot_gray, _, _ = load_grayscale_image(spatial_plot_sample)

    # 组织灰度按原2%/98%分位裁切，再映射至0.12至0.90，保持背景对比度。
    spatial_plot_gray_low = float(np.nanpercentile(spatial_plot_gray, 2))
    spatial_plot_gray_high = float(np.nanpercentile(spatial_plot_gray, 98))
    spatial_plot_gray_scaled = np.clip(
        (spatial_plot_gray - spatial_plot_gray_low) / max(spatial_plot_gray_high - spatial_plot_gray_low, 1e-6),
        0,
        1,
    )
    spatial_plot_gray_img = 0.12 + 0.78 * spatial_plot_gray_scaled

    # 读回两切片源表，逐指标从合并点表计算共享2%和98%分位色限。
    spatial_plot_all_spots_table_total = pd.concat(
        [
            pd.read_csv(fig3_source_out / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv"),
            pd.read_csv(fig3_source_out / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv"),
        ],
        ignore_index = True,
    )

    # 两张组织图均显示raw20、MDSC、MES三指标，色带统一为viridis。
    metric_meta = {
        "raw20": {"title": "Mg-inflammatory raw20", "cmap": "viridis"},
        "mdsc_combined": {"title": "MDSC score", "cmap": "viridis"},
        "mes": {"title": "MES score", "cmap": "viridis"},
    }

    spatial_plot_limits = {
        spatial_plot_metric: np.nanquantile(
            spatial_plot_all_spots_table_total[spatial_plot_metric].to_numpy(dtype=float), [0.02, 0.98]
        )
        for spatial_plot_metric in metric_meta
    }

    # 183×70毫米画布横排三个组织图和底部水平色条。
    spatial_plot_fig = new_figure(183, 70)
    spatial_plot_lefts = [0.04, 0.36, 0.68]
    spatial_plot_axes = [
        spatial_plot_fig.add_axes([spatial_plot_left, 0.19, 0.23, 0.63]) for spatial_plot_left in spatial_plot_lefts
    ]

    spatial_plot_cbar_axes = [
        spatial_plot_fig.add_axes([spatial_plot_left, 0.08, 0.23, 0.03]) for spatial_plot_left in spatial_plot_lefts
    ]

    for spatial_plot_ax, spatial_plot_cax, (spatial_plot_metric, spatial_plot_meta) in zip(
        spatial_plot_axes, spatial_plot_cbar_axes, metric_meta.items()
    ):
        spatial_plot_ax.imshow(spatial_plot_gray_img, cmap="gray", origin="upper", interpolation="nearest")
        spatial_plot_mask = spatial_plot_gray_img < 0.96
        spatial_plot_ax.contour(
            spatial_plot_mask.astype(float), levels=[0.5], colors=["#8B8B8B"], linewidths=0.55, origin="upper"
        )

        # 灰色大点形成底层轮廓，上层较小彩点显示当前指标分数。
        spatial_plot_ax.scatter(
            spatial_plot_data_table["x_hires"],
            spatial_plot_data_table["y_hires"],
            s          = 11,
            color      = "#D1D1D1",
            alpha      = 0.85,
            linewidths = 0,
            rasterized = True,
            zorder     = 2,
        )

        spatial_plot_points = spatial_plot_ax.scatter(
            spatial_plot_data_table["x_hires"],
            spatial_plot_data_table["y_hires"],
            s          = 7.5,
            c          = spatial_plot_data_table[spatial_plot_metric],
            cmap       = spatial_plot_meta["cmap"],
            vmin       = float(spatial_plot_limits[spatial_plot_metric][0]),
            vmax       = float(spatial_plot_limits[spatial_plot_metric][1]),
            linewidths = 0,
            rasterized = True,
            zorder     = 3,
        )

        # 坐标范围严格采用组织图宽高，保留左上原点并隐藏坐标轴。
        spatial_plot_ax.set_title(spatial_plot_meta["title"], pad=3)
        spatial_plot_ax.set_xticks([])
        spatial_plot_ax.set_yticks([])
        for spatial_plot_spine in spatial_plot_ax.spines.values():
            spatial_plot_spine.set_visible(False)
        spatial_plot_ax.set_xlim(0, spatial_plot_gray.shape[1])
        spatial_plot_ax.set_ylim(spatial_plot_gray.shape[0], 0)
        spatial_plot_colorbar = spatial_plot_fig.colorbar(
            spatial_plot_points, cax=spatial_plot_cax, orientation="horizontal"
        )

        spatial_plot_colorbar.outline.set_visible(False)
        spatial_plot_colorbar.ax.tick_params(labelsize=6.2, length=1.8, pad=1)

    spatial_plot_fig.text(0.04, 0.97, spatial_plot_sample, ha="left", va="top", fontsize=8.5, weight="bold")
    spatial_plot_fig.text(
        0.96,
        0.97,
        f"{len(spatial_plot_data_table):,} in-tissue Visium spots | shared scale across slices",
        ha       = "right",
        va       = "top",
        fontsize = 7.0,
        color    = colors["neutral"],
    )

    # 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
    save_figure(spatial_plot_fig, fig3_out, spatial_plot_stem)

# 主相关森林图按MDSC两切片、MES两切片固定顺序展示点估计与Fisher区间。
rho_plot_data_table = fig3e_table

rho_plot_ordered_table = (
    rho_plot_data_table.set_index("panel_label")
    .loc[["GBM030 · MDSC", "GBM049 · MDSC", "GBM030 · MES", "GBM049 · MES"]]
    .reset_index()
)
rho_plot_fig = new_figure(88, 82)
rho_plot_ax = rho_plot_fig.add_axes([0.28, 0.12, 0.66, 0.68])

rho_plot_y_values = np.arange(len(rho_plot_ordered_table))

partial_colors = {"MDSC": "#00A087", "MES": "#7E6148", "E-MDSC": "#4DBBD5", "M-MDSC": "#3C5488"}

for rho_plot_yy, rho_plot_row in zip(rho_plot_y_values, rho_plot_ordered_table.itertuples(index=False)):
    rho_plot_color = partial_colors[rho_plot_row.target]
    rho_plot_ax.hlines(rho_plot_yy, rho_plot_row.ci_low, rho_plot_row.ci_high, color=rho_plot_color, linewidth=1.2)
    rho_plot_ax.scatter(rho_plot_row.partial_rho, rho_plot_yy, s=28, color=rho_plot_color, zorder=3)
    rho_plot_ax.text(
        rho_plot_row.partial_rho + 0.018,
        rho_plot_yy,
        f"n={int(rho_plot_row.n_spots)}",
        va       = "center",
        ha       = "left",
        fontsize = 6.4,
        color    = colors["text"],
    )

# 保留零线和原横轴留白；标题注明切片层面的描述性相关解释。
rho_plot_ax.axvline(0, color=colors["neutral"], linewidth=0.7)
rho_plot_ax.set_yticks(rho_plot_y_values, rho_plot_ordered_table["panel_label"])
rho_plot_ax.invert_yaxis()
rho_plot_ax.set_xlabel("Partial Spearman rho")
clean_axis(rho_plot_ax)
rho_plot_ax.grid(axis="x", color=colors["grid"], linewidth=0.6)
rho_plot_ax.set_xlim(0.0, max(0.64, float(rho_plot_ordered_table["ci_high"].max()) + 0.13))
rho_plot_fig.text(0.28, 0.97, "Slice-level spatial effect sizes", ha="left", va="top", fontsize=8.5, weight="bold")
rho_plot_fig.text(
    0.28,
    0.92,
    "Descriptive only; adjusted for total myeloid score\nNo smoothing; no spatial-autocorrelation correction",
    ha          = "left",
    va          = "top",
    fontsize    = 6.6,
    linespacing = 1.15,
    color       = colors["neutral"],
)

# 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
save_figure(rho_plot_fig, fig3_out, "fig3e_partial_rho_effects")

# 敏感性森林图分左右两轴放分层标签及效应区间，各分节之间留0.55行。
sensitivity_plot_data_table = s5a_table

sensitivity_forest_data_table = sensitivity_plot_data_table

sensitivity_forest_sections = []
sensitivity_forest_y_cursor = 0.0
sensitivity_forest_row_positions: list[float] = []
for sensitivity_forest_section in sensitivity_forest_data_table["section"].drop_duplicates():
    sensitivity_forest_sections.append((sensitivity_forest_section, sensitivity_forest_y_cursor - 0.55))
    sensitivity_forest_current_table = sensitivity_forest_data_table.loc[
        sensitivity_forest_data_table["section"].eq(sensitivity_forest_section)
    ]
    for _ in range(len(sensitivity_forest_current_table)):
        sensitivity_forest_row_positions.append(sensitivity_forest_y_cursor)
        sensitivity_forest_y_cursor += 1.0
    sensitivity_forest_y_cursor += 0.55
sensitivity_forest_data_table = sensitivity_forest_data_table.copy()
sensitivity_forest_data_table["y"] = sensitivity_forest_row_positions

sensitivity_forest_fig = new_figure(120, 140)
sensitivity_forest_label_ax = sensitivity_forest_fig.add_axes([0.04, 0.13, 0.40, 0.76])
sensitivity_forest_forest_ax = sensitivity_forest_fig.add_axes([0.47, 0.13, 0.49, 0.76])
sensitivity_forest_label_ax.axis("off")

# 横轴优先覆盖全部有限区间；无有限区间时沿用原效应范围回退。
sensitivity_forest_x_min = float(np.nanmin(sensitivity_forest_data_table["ci_low"].to_numpy(dtype=float)))
sensitivity_forest_x_max = float(np.nanmax(sensitivity_forest_data_table["ci_high"].to_numpy(dtype=float)))
if not np.isfinite(sensitivity_forest_x_min):
    sensitivity_forest_x_min = float(np.nanmin(sensitivity_forest_data_table["effect"].to_numpy(dtype=float)))
if not np.isfinite(sensitivity_forest_x_max):
    sensitivity_forest_x_max = float(np.nanmax(sensitivity_forest_data_table["effect"].to_numpy(dtype=float)))
sensitivity_forest_span = max(sensitivity_forest_x_max - sensitivity_forest_x_min, 0.16)
sensitivity_forest_forest_ax.set_xlim(
    sensitivity_forest_x_min - sensitivity_forest_span * 0.16, sensitivity_forest_x_max + sensitivity_forest_span * 0.58
)

for sensitivity_forest_section, sensitivity_forest_section_y in sensitivity_forest_sections:
    sensitivity_forest_label_ax.text(
        0.0,
        sensitivity_forest_section_y,
        sensitivity_forest_section,
        ha       = "left",
        va       = "center",
        fontsize = 6.5,
        color    = colors["neutral"],
        weight   = "bold",
    )

# 逐行绘制有限区间、点估计、分层名称，并显示人数、可用上升数及FDR。
for sensitivity_forest_row in sensitivity_forest_data_table.itertuples(index=False):
    sensitivity_forest_color = colors["recurrent"] if sensitivity_forest_row.effect >= 0 else colors["primary"]
    if np.isfinite(sensitivity_forest_row.ci_low) and np.isfinite(sensitivity_forest_row.ci_high):
        sensitivity_forest_forest_ax.hlines(
            sensitivity_forest_row.y,
            sensitivity_forest_row.ci_low,
            sensitivity_forest_row.ci_high,
            color     = sensitivity_forest_color,
            linewidth = 1.1,
        )

    sensitivity_forest_forest_ax.scatter(
        sensitivity_forest_row.effect, sensitivity_forest_row.y, s=24, color=sensitivity_forest_color, zorder=3
    )
    sensitivity_forest_label_ax.text(
        0.0,
        sensitivity_forest_row.y,
        sensitivity_forest_row.label,
        ha       = "left",
        va       = "center",
        fontsize = 7.0,
        color    = colors["text"],
    )

    sensitivity_forest_side_bits = [f"n={int(sensitivity_forest_row.n_pairs)}"]
    if pd.notna(sensitivity_forest_row.n_up):
        sensitivity_forest_side_bits.append(f"{int(sensitivity_forest_row.n_up)} up")
    sensitivity_forest_side_bits.append(f"FDR={sensitivity_forest_row.fdr:.3f}")
    sensitivity_forest_forest_ax.text(
        sensitivity_forest_row.effect + sensitivity_forest_span * 0.08,
        sensitivity_forest_row.y,
        "  ".join(sensitivity_forest_side_bits),
        ha       = "left",
        va       = "center",
        fontsize = 6.2,
        color    = colors["text"],
    )

# 统一左右轴纵向范围与倒序，保留同一效应尺度和零线。
sensitivity_forest_forest_ax.axvline(0, color=colors["neutral"], linewidth=0.7)
sensitivity_forest_forest_ax.set_yticks([])
sensitivity_forest_forest_ax.set_xlabel("Mean paired raw20 delta")
clean_axis(sensitivity_forest_forest_ax, keep_left=False, keep_bottom=True)
sensitivity_forest_forest_ax.grid(axis="x", color=colors["grid"], linewidth=0.6)
sensitivity_forest_forest_ax.invert_yaxis()
sensitivity_forest_label_ax.set_ylim(sensitivity_forest_forest_ax.get_ylim())
sensitivity_forest_fig.text(0.04, 0.97, "GeoMx sensitivity analyses", ha="left", va="top", fontsize=8.5, weight="bold")
sensitivity_forest_fig.text(
    0.04,
    0.93,
    "Entry, IDH, treatment and interaction on the same effect scale",
    ha       = "left",
    va       = "top",
    fontsize = 7.0,
    color    = colors["neutral"],
)

# 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
save_figure(sensitivity_forest_fig, s5_out, "s5a_geomx_sensitivity_forest")

# 切片敏感性先分MDSC家族与MES两节，再保留每节输入行顺序。
slice_plot_data_table = s5b_table

slice_plot_order = []
for slice_plot_section, slice_plot_targets in [("MDSC family", ["MDSC", "E-MDSC", "M-MDSC"]), ("MES", ["MES"])]:
    slice_plot_current_table = slice_plot_data_table.loc[
        slice_plot_data_table["target"].isin(slice_plot_targets)
    ].copy()
    slice_plot_current_table["section"] = slice_plot_section
    slice_plot_current_table["label"] = slice_plot_current_table["sample"] + " · " + slice_plot_current_table["target"]

    slice_plot_order.append(slice_plot_current_table)
slice_plot_forest_table_total = pd.concat(slice_plot_order, ignore_index=True)
slice_plot_forest_table_total["effect"] = slice_plot_forest_table_total["partial_rho"]
slice_plot_forest_table_total = slice_plot_forest_table_total[
    ["section", "label", "effect", "ci_low", "ci_high", "n_spots", "fdr", "target"]
].rename(columns={"n_spots": "n_pairs"})
slice_plot_forest_table_total["n_up"] = np.nan
slice_plot_forest_table_total["section_order"] = slice_plot_forest_table_total["section"].map(
    {"MDSC family": 1, "MES": 2}
)
slice_plot_forest_table_total["row_order"] = range(len(slice_plot_forest_table_total))

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s5_source_out / "s5b_slice_sensitivity_forest.csv"
slice_plot_forest_table_total.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 以分节间距构造绘图纵坐标；覆盖同名源表时保留原下游字段合同。
slice_plot_sections = []
slice_plot_y_cursor = 0.0
slice_plot_row_positions: list[float] = []
for slice_plot_section in slice_plot_forest_table_total["section"].drop_duplicates():
    slice_plot_sections.append((slice_plot_section, slice_plot_y_cursor - 0.35))

    slice_plot_current_table = slice_plot_forest_table_total.loc[
        slice_plot_forest_table_total["section"].eq(slice_plot_section)
    ]
    for _ in range(len(slice_plot_current_table)):
        slice_plot_row_positions.append(slice_plot_y_cursor)
        slice_plot_y_cursor += 1.0
    slice_plot_y_cursor += 0.55
slice_plot_forest_table_total = slice_plot_forest_table_total.copy()
slice_plot_forest_table_total["y"] = slice_plot_row_positions

slice_plot_fig = new_figure(120, 118)
slice_plot_label_ax = slice_plot_fig.add_axes([0.04, 0.10, 0.36, 0.70])

slice_plot_forest_ax = slice_plot_fig.add_axes([0.44, 0.10, 0.52, 0.70])
slice_plot_label_ax.axis("off")
for slice_plot_section, slice_plot_section_y in slice_plot_sections:
    slice_plot_label_ax.text(
        0.0,
        slice_plot_section_y,
        slice_plot_section,
        ha       = "left",
        va       = "center",
        fontsize = 6.5,
        color    = colors["neutral"],
        weight   = "bold",
    )

# 切片森林图横轴覆盖所有Fisher区间，右侧留出有效点数标注。
slice_plot_x_min = float(np.nanmin(slice_plot_forest_table_total["ci_low"]))
slice_plot_x_max = float(np.nanmax(slice_plot_forest_table_total["ci_high"]))
slice_plot_span = max(slice_plot_x_max - slice_plot_x_min, 0.18)
slice_plot_forest_ax.set_xlim(slice_plot_x_min - slice_plot_span * 0.16, slice_plot_x_max + slice_plot_span * 0.50)
for slice_plot_row in slice_plot_forest_table_total.itertuples(index=False):
    slice_plot_color = partial_colors[slice_plot_row.target]
    slice_plot_forest_ax.hlines(
        slice_plot_row.y, slice_plot_row.ci_low, slice_plot_row.ci_high, color=slice_plot_color, linewidth=1.1
    )
    slice_plot_forest_ax.scatter(slice_plot_row.effect, slice_plot_row.y, s=24, color=slice_plot_color, zorder=3)

    slice_plot_label_ax.text(
        0.0, slice_plot_row.y, slice_plot_row.label, ha="left", va="center", fontsize=7.0, color=colors["text"]
    )
    slice_plot_forest_ax.text(
        slice_plot_row.effect + slice_plot_span * 0.08,
        slice_plot_row.y,
        f"n={int(slice_plot_row.n_pairs)}",
        ha       = "left",
        va       = "center",
        fontsize = 6.3,
        color    = colors["text"],
    )

# 以零线、纵向对齐及原描述性说明完成切片敏感性面板。
slice_plot_forest_ax.axvline(0, color=colors["neutral"], linewidth=0.7)
slice_plot_forest_ax.set_yticks([])
slice_plot_forest_ax.set_xlabel("Partial Spearman rho")
clean_axis(slice_plot_forest_ax, keep_left=False, keep_bottom=True)
slice_plot_forest_ax.grid(axis="x", color=colors["grid"], linewidth=0.6)
slice_plot_forest_ax.invert_yaxis()
slice_plot_label_ax.set_ylim(slice_plot_forest_ax.get_ylim())
slice_plot_fig.text(0.04, 0.97, "Spatial slice sensitivity", ha="left", va="top", fontsize=8.5, weight="bold")
slice_plot_fig.text(
    0.04,
    0.91,
    "Descriptive partial rho only; each row is one slice-by-target partial correlation; no smoothing",
    ha       = "left",
    va       = "top",
    fontsize = 7.0,
    color    = colors["neutral"],
)

# 沿用Step41共享保存接口；PDF参数不变，PNG由该PDF渲染。
save_figure(slice_plot_fig, s5_out, "s5b_slice_sensitivity_forest")

# Figure3清单逐项连接五个面板、源CSV及PDF/PNG路径；字段维持原六列。
manifest_fig3_manifest_table = pd.DataFrame(
    [
        (
            "Figure3A",
            "fig3a_geomx22_paired_change",
            "GeoMx paired raw20 change",
            fig3_source_out / "fig3a_geomx22_paired_change.csv",
        ),
        (
            "Figure3B",
            "fig3b_geomx18_gene_forest",
            "Measurable gene forest",
            fig3_source_out / "fig3b_geomx18_gene_forest.csv",
        ),
        (
            "Figure3C",
            "fig3c_gbm030_spatial_raw20_mdsc_mes",
            "GBM030 spatial maps",
            fig3_source_out / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv",
        ),
        (
            "Figure3D",
            "fig3d_gbm049_spatial_raw20_mdsc_mes",
            "GBM049 spatial maps",
            fig3_source_out / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv",
        ),
        (
            "Figure3E",
            "fig3e_partial_rho_effects",
            "Main slice partial-rho effect plot",
            fig3_source_out / "fig3e_partial_rho_main.csv",
        ),
    ],
    columns = ["panel", "stem", "message", "source_csv"],
)

manifest_fig3_manifest_table["png"] = manifest_fig3_manifest_table["stem"].map(
    lambda stem: str(fig3_out / f"{stem}.png")
)
manifest_fig3_manifest_table["pdf"] = manifest_fig3_manifest_table["stem"].map(
    lambda stem: str(fig3_out / f"{stem}.pdf")
)
manifest_fig3_manifest_table["source_csv"] = manifest_fig3_manifest_table["source_csv"].map(str)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = write_fig3 / "Figure3_panel_manifest.csv"
manifest_fig3_manifest_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# SupplementaryFigure5清单连接两张敏感性面板及其源CSV和图像。
manifest_s5_manifest_table = pd.DataFrame(
    [
        (
            "SupplementaryFigure5A",
            "s5a_geomx_sensitivity_forest",
            "GeoMx sensitivity forest",
            s5_source_out / "s5a_geomx_sensitivity_forest.csv",
        ),
        (
            "SupplementaryFigure5B",
            "s5b_slice_sensitivity_forest",
            "Spatial slice sensitivity forest",
            s5_source_out / "s5b_slice_sensitivity_forest.csv",
        ),
    ],
    columns = ["panel", "stem", "message", "source_csv"],
)

manifest_s5_manifest_table["png"] = manifest_s5_manifest_table["stem"].map(lambda stem: str(s5_out / f"{stem}.png"))
manifest_s5_manifest_table["pdf"] = manifest_s5_manifest_table["stem"].map(lambda stem: str(s5_out / f"{stem}.pdf"))
manifest_s5_manifest_table["source_csv"] = manifest_s5_manifest_table["source_csv"].map(str)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = write_s5 / "SupplementaryFigure5_panel_manifest.csv"
manifest_s5_manifest_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
print(
    "Figure3 panels:", len(manifest_fig3_manifest_table),
    "SupplementaryFigure5 panels:", len(manifest_s5_manifest_table)
)
