#!/usr/bin/env python3
# 目的：用既有raw20状态效应及作者marker重算结果生成Figure2两面板。
# 输入：write/48_figure2_identity_state_rebuild/Figure2/GSE278456_*.csv、细胞注释及原始H5。
# 输出：write/50_figure2_two_panel_final/Figure2内源表、清单、图注与结果说明；同名figures目录内PDF、PNG及预览。
# 运行：在项目根使用既有sc环境执行 python scripts/06c_Figure2主图.py。
# 保留Python：Scanpy、AnnData及Matplotlib是已公开复现路径；改成R超出本轮纯风格整理。
# 存量source_data/config、数据文件名、输出字段、write及figures路径保留原公开合同。

from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import scanpy as sc
import seaborn as sns
import os
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png
sha256 = figure_style.sha256

# --plot-only读取本步骤已冻结的marker源表；默认分支重读H5核对Step48赋值。
arguments_parser = argparse.ArgumentParser(description="生成 Figure2 两面板正式候选版")
arguments_parser.add_argument("--plot-only", action="store_true", help="只从本 Step 已写出的冻结源数据重画，不重读 H5")
args = arguments_parser.parse_args()

# 沿用原图的字体、线宽和背景配置，以保持PDF画面。
mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.2,
        "axes.labelsize": 7.6,
        "axes.titlesize": 8.2,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)

# 沿用原Scanpy设置；本轮不更换库或评分实现。
sc.settings.verbosity = 0

# 创建Figure2源表和面板输出目录，保留公开版本的write/figures目录名。
root = Path(__file__).resolve().parents[1]

write_root = root / "write/50_figure2_two_panel_final/Figure2"

source_out = write_root / "source_data"

source_out.mkdir(parents=True, exist_ok=True)

fig_root = root / "figures/50_figure2_two_panel_final/Figure2"

fig_root.mkdir(parents=True, exist_ok=True)

write_root.mkdir(parents=True, exist_ok=True)

step48 = root / "write/48_figure2_identity_state_rebuild/Figure2"

# 读取患者状态差值、九状态统计、基因表达及marker审查四张Step48冻结表。
effect_file = step48 / "GSE278456_patient_state_effects.csv"

inputs_effects_table = pd.read_csv(effect_file)

summary_file = step48 / "GSE278456_state_effect_summary.csv"

inputs_summary_table = pd.read_csv(summary_file)

gene_summary_file = step48 / "GSE278456_gene_summary.csv"

# 读取该步骤原有输入表，保留列名、缺失值处理和读入顺序。
inputs_gene_summary_table = pd.read_csv(gene_summary_file)

marker_audit_file = step48 / "GSE278456_author_marker_audit.csv"

inputs_marker_audit_table = pd.read_csv(marker_audit_file)

# 108个患者状态效应、9个状态和171个基因状态组合是既有冻结尺寸。
if inputs_effects_table.shape[0] != 108 or inputs_effects_table["state"].nunique() != 9:
    raise RuntimeError("Step48 患者状态效应表不再是 108 行 / 9 状态")

# 展示顺序沿用既有九状态；仅影响排版，不重排marker竞争顺序。
state_order = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]

if set(inputs_summary_table["state"]) != set(state_order):
    raise RuntimeError("Step48 状态全集发生变化")
if inputs_gene_summary_table.shape[0] != 171:
    raise RuntimeError("Step48 raw20 基因状态表不再是 171 行")

# Miller原始排序前20基因固定不变；基因大写是外部数据标识。
raw20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]

# AC253572.2为Step48冻结的唯一不可测raw20项，核对剩余19基因及九状态覆盖。
measurable_raw20 = [gene for gene in raw20 if gene != "AC253572.2"]

if set(inputs_gene_summary_table["gene"]) != set(measurable_raw20):
    raise RuntimeError("Step48 可测 raw20 基因集发生变化")
if inputs_gene_summary_table.groupby("gene")["state"].nunique().ne(9).any():
    raise RuntimeError("raw20 基因未完整覆盖 9 状态")

# 收集九状态实际使用的marker，检查其与19个可测raw20基因完全分离。
inputs_focus_audit_table = inputs_marker_audit_table.loc[inputs_marker_audit_table["state"].isin(state_order)].copy()
if inputs_focus_audit_table["state"].nunique() != 9:
    raise RuntimeError("marker audit 不完整")

inputs_used = set()
for inputs_value in inputs_focus_audit_table["markers_used"].fillna(""):
    inputs_used.update(inputs_gene for inputs_gene in str(inputs_value).split(";") if inputs_gene)
inputs_overlap = sorted(inputs_used & set(measurable_raw20))
if inputs_overlap:
    raise RuntimeError(f"赋值后 marker 与 raw20 仍重叠: {inputs_overlap}")
(effects_table, state_summary_table, gene_summary_table, _) = (
    inputs_effects_table,
    inputs_summary_table,
    inputs_gene_summary_table,
    inputs_marker_audit_table,
)

# 瀑布图源表以患者状态差值为主体，补入状态统计与pseudobulk方向标记。
waterfall_effects_table = effects_table
waterfall_summary_table = state_summary_table

waterfall_summary_index_table = waterfall_summary_table.set_index("state")

# 按固定状态顺序逐组整理患者排序，原始差值和显著性字段不变。
waterfall_rows: list[pd.DataFrame] = []
for waterfall_state in state_order:
    waterfall_data_table = (
        waterfall_effects_table.loc[waterfall_effects_table["state"].eq(waterfall_state)]
        .sort_values("raw20_delta")
        .copy()
    )
    waterfall_data_table["within_state_rank"] = np.arange(1, len(waterfall_data_table) + 1)
    waterfall_data_table["direction"] = np.where(waterfall_data_table["raw20_delta"] > 0, "Above rest", "Below rest")
    waterfall_data_table["n_patients"] = int(waterfall_summary_index_table.loc[waterfall_state, "n_patients"])

    waterfall_data_table["n_positive"] = int(waterfall_summary_index_table.loc[waterfall_state, "n_positive"])
    waterfall_data_table["mean_delta"] = float(waterfall_summary_index_table.loc[waterfall_state, "mean_delta"])
    waterfall_data_table["fdr_9_states"] = float(waterfall_summary_index_table.loc[waterfall_state, "fdr_9_states"])
    waterfall_data_table["score_pseudobulk_mean_direction_concordant"] = bool(
        waterfall_summary_index_table.loc[waterfall_state, "score_pseudobulk_mean_direction_concordant"]
    )
    waterfall_rows.append(waterfall_data_table)
waterfall_table_total = pd.concat(waterfall_rows, ignore_index=True)

# 三张marker源表分别保存患者均值、9×9汇总和逐样本赋值复核结果。
marker_patient_file = source_out / "GSE278456_patient_marker_set_scores_source.csv"

marker_summary_file = source_out / "GSE278456_marker_set_assignment_qc_source.csv"
marker_audit_file = source_out / "GSE278456_step48_assignment_reproduction_source.csv"

# 绘图分支要求三张冻结表齐全；默认分支从工作簿与H5重新构造同一赋值。
if args.plot_only:
    required = [marker_patient_file, marker_summary_file, marker_audit_file]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError(f"--plot-only 缺少冻结源数据: {missing}")
    marker_patient_table = pd.read_csv(marker_patient_file)
    marker_summary_table = pd.read_csv(marker_summary_file)
    audit_table = pd.read_csv(marker_audit_file)
else:

    # 只复用Step48的marker读取函数；该模块的主分析受__main__入口保护。
    step48_script = root / "scripts/06b_细胞身份与状态图.py"
    loader_spec = importlib.util.spec_from_file_location("step48_fig2", step48_script)
    if loader_spec is None or loader_spec.loader is None:
        raise RuntimeError("Step48 模块无法加载")
    loader_module = importlib.util.module_from_spec(loader_spec)
    loader_spec.loader.exec_module(loader_module)
    marker_step48 = loader_module
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    marker_book = external_root / "write/30_dataset_rescue_search/source_metadata/NIHMS2115262-supplement-S2_S5_S6.xls"
    marker_markers, _ = marker_step48.read_author_markers(marker_book, raw20, state_order)
    if len(marker_markers) != 17:
        raise RuntimeError(f"工作簿应提供 17 套状态 marker，实际 {len(marker_markers)}")

    # 读取该步骤原有输入表，保留列名、缺失值处理和读入顺序。
    # 加载样本清单和旧QC计数，再按文件名排序锁定21份H5。
    marker_manifest_table = pd.read_csv((root / "write/30_dataset_rescue_search") / "GSE278456_sample_manifest.csv")
    marker_manifest_table["geo_accession"] = marker_manifest_table["geo_accession"].astype(str)
    marker_existing_qc_table = pd.read_csv(
        (root / "write/30_dataset_rescue_search") / "GSE278456_single_cell_qc.csv"
    ).set_index("geo_accession")
    h5_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE278456_myeloid_h5"
    marker_paths = sorted(h5_root.glob("*.h5"))
    if len(marker_paths) != 21:
        raise RuntimeError(f"GSE278456 H5 数量应为 21，实际 {len(marker_paths)}")

    # 读取Step48冻结细胞标签；用cell_id对齐，不按行号匹配。
    cell_file = step48 / "GSE278456_cell_annotations_true_raw20.csv.gz"
    marker_saved_cells_table = pd.read_csv(
        cell_file, usecols=["cell_id", "geo_accession", "author_marker_state", "state_score"]
    )
    if marker_saved_cells_table.shape[0] != 120766 or marker_saved_cells_table["cell_id"].duplicated().any():
        raise RuntimeError("Step48 冻结 cell annotation 不再是 120,766 个唯一细胞")
    marker_saved_cells_table = marker_saved_cells_table.set_index("cell_id")

    marker_all_marker_genes = set().union(*marker_markers.values())
    marker_patient_rows: list[dict[str, object]] = []
    marker_audit_rows: list[dict[str, object]] = []

    # 逐H5核对IDH-wt入口及QC细胞数；每样本保留原工作簿marker竞争顺序。
    for marker_path in marker_paths:
        marker_accession = marker_path.name.split("_")[0]
        marker_manifest_match_table = marker_manifest_table.loc[
            marker_manifest_table["geo_accession"].eq(marker_accession)
        ]
        if len(marker_manifest_match_table) != 1:
            raise RuntimeError(f"{marker_accession} 无法唯一匹配 manifest")
        marker_manifest_row = marker_manifest_match_table.iloc[0]
        if not bool(marker_manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"锁定 H5 出现非 IDH-wt GBM: {marker_accession}")
        marker_sample = str(marker_manifest_row["sample_title"])

        # 从10x原始计数计算基因数、总计数和线粒体比例，为同一QC门槛提供输入。
        marker_adata = sc.read_10x_h5(marker_path)
        marker_adata.var_names_make_unique()
        marker_n_input = marker_adata.n_obs
        marker_adata.var["mt"] = marker_adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(marker_adata, qc_vars=["mt"], inplace=True, percent_top=None)

        # 与Step48使用相同的原QC阈值，避免状态重算改变输入细胞。
        marker_keep = (
            marker_adata.obs["n_genes_by_counts"].between(200, 8000)
            & marker_adata.obs["total_counts"].ge(500)
            & marker_adata.obs["pct_counts_mt"].le(20)
        )
        marker_adata = marker_adata[marker_keep].copy()
        marker_n_qc = marker_adata.n_obs
        marker_old = marker_existing_qc_table.loc[marker_accession]
        if int(marker_old["n_input"]) != marker_n_input or int(marker_old["n_qc_pass"]) != marker_n_qc:
            raise RuntimeError(f"{marker_accession} QC 细胞数未复现")

        # 检查当前H5的可测raw20集合与冻结表完全同序。
        marker_present_raw20 = [marker_gene for marker_gene in raw20 if marker_gene in marker_adata.var_names]
        if marker_present_raw20 != measurable_raw20:
            raise RuntimeError(f"{marker_accession} raw20 覆盖或顺序异常")
        marker_adata.obs_names = [
            f"{marker_accession}:{marker_barcode}" for marker_barcode in marker_adata.obs_names.astype(str)
        ]

        # 按原10,000计数尺度归一化，随后log1p；不改变原始计数副本。
        sc.pp.normalize_total(marker_adata, target_sum=10000)
        sc.pp.log1p(marker_adata)

        # 从对照池排除raw20和全部作者marker，保留原至少5000基因门槛。
        marker_excluded_controls = set(raw20) | marker_all_marker_genes
        marker_neutral_control_pool = [
            marker_gene for marker_gene in marker_adata.var_names if marker_gene not in marker_excluded_controls
        ]
        if len(marker_neutral_control_pool) < 5000:
            raise RuntimeError(f"{marker_accession} state-neutral control pool 过小")

        # 对工作簿17套marker分别评分，每套至少30个可测基因后才进入竞争。
        marker_score_columns: list[str] = []
        for marker_state, marker_genes in marker_markers.items():
            marker_use = [marker_gene for marker_gene in marker_genes if marker_gene in marker_adata.var_names]
            if len(marker_use) < 30:
                raise RuntimeError(f"{marker_accession} {marker_state} 可测 marker 少于 30")
            marker_column = f"state_score_{marker_state}"

            # 保留ctrl_size=50、n_bins=25和原种子；每套评分重置种子是原方法，移出会改变结果。
            sc.tl.score_genes(
                marker_adata,
                gene_list    = marker_use,
                gene_pool    = marker_neutral_control_pool + marker_use,
                ctrl_size    = 50,
                n_bins       = 25,
                score_name   = marker_column,
                random_state = 20260713,
                ctrl_as_ref  = False,
                use_raw      = False,
            )

            marker_score_columns.append(marker_column)

        # 在工作簿全部marker集之间按原argsort竞争，再截取预设髓系状态。
        marker_score_matrix = marker_adata.obs[marker_score_columns].to_numpy(float)
        marker_order = np.argsort(marker_score_matrix, axis=1)
        marker_state_names = np.array(list(marker_markers), dtype=object)
        marker_winner_index = marker_order[:, -1]
        marker_winners = marker_state_names[marker_winner_index]
        marker_winning_score = marker_score_matrix[np.arange(marker_adata.n_obs), marker_winner_index]

        # 比较同一细胞的winner标签及分数；保留原1e-10精确复现容差。
        marker_frozen_table = marker_saved_cells_table.loc[marker_adata.obs_names]
        marker_winner_matches = marker_frozen_table["author_marker_state"].to_numpy(object) == marker_winners
        marker_score_difference = np.abs(marker_frozen_table["state_score"].to_numpy(float) - marker_winning_score)
        if not bool(marker_winner_matches.all()) or float(marker_score_difference.max()) > 1e-10:
            raise RuntimeError(f"{marker_accession} 未精确复现 Step48 状态赋值")
        marker_audit_rows.append(
            {
                "geo_accession": marker_accession,
                "sample": marker_sample,
                "n_input": marker_n_input,
                "n_qc_pass": marker_n_qc,
                "n_winner_exact_match": int(marker_winner_matches.sum()),
                "winner_match_fraction": float(marker_winner_matches.mean()),
                "max_abs_winner_score_difference": float(marker_score_difference.max()),
                "n_competing_marker_sets": len(marker_markers),
            }
        )

        # 只对九个髓系获胜状态汇总；每个患者状态至少30个细胞才计算九套marker均值。
        marker_focus = np.isin(marker_winners, state_order)
        marker_focus_scores = {
            marker_state: marker_adata.obs[f"state_score_{marker_state}"].to_numpy(float)
            for marker_state in state_order
        }
        for marker_assigned_state in state_order:
            marker_target = marker_focus & (marker_winners == marker_assigned_state)
            marker_n_cells = int(marker_target.sum())
            if marker_n_cells < 30:
                continue

            for marker_tested_marker_set in state_order:
                marker_patient_rows.append(
                    {
                        "geo_accession": marker_accession,
                        "sample": marker_sample,
                        "assigned_state": marker_assigned_state,
                        "tested_marker_set": marker_tested_marker_set,
                        "n_cells": marker_n_cells,
                        "mean_marker_set_score": float(
                            marker_focus_scores[marker_tested_marker_set][marker_target].mean()
                        ),
                    }
                )

        del marker_adata, marker_score_matrix

    # 每行是患者、被赋予状态和被测试marker集的组合，后续仍按患者等权。
    marker_patient_scores_table = pd.DataFrame(marker_patient_rows)
    marker_audit_table = pd.DataFrame(marker_audit_rows).sort_values("geo_accession")
    if (
        int(marker_audit_table["n_qc_pass"].sum()) != 120766
        or not marker_audit_table["winner_match_fraction"].eq(1).all()
    ):
        raise RuntimeError("21 个 H5 的 Step48 赋值复现未全部通过")
    if marker_patient_scores_table.shape[0] != 972:
        raise RuntimeError(f"患者等权 marker score 长表应为 972 行，实际 {marker_patient_scores_table.shape[0]}")

    # 先对同一标签状态的患者均值等权平均，再在每套marker内部跨九状态标准化。
    marker_summary_table = (
        marker_patient_scores_table.groupby(["tested_marker_set", "assigned_state"], observed=True)
        .agg(n_patients=("sample", "nunique"), patient_balanced_mean_marker_set_score=("mean_marker_set_score", "mean"))
        .reset_index()
    )
    marker_summary_table["within_marker_set_z"] = np.nan

    # 使用总体标准差ddof=0计算组内z分数；零方差时填0，并按原索引写回。
    for _, group_table in marker_summary_table.groupby("tested_marker_set", observed=True):
        group_values = group_table["patient_balanced_mean_marker_set_score"]
        group_sd = group_values.std(ddof=0)
        if group_sd > 0:
            group_z = (group_values - group_values.mean()) / group_sd
        else:
            group_z = 0
        marker_summary_table.loc[group_table.index, "within_marker_set_z"] = group_z
    if marker_summary_table.shape[0] != 81:
        raise RuntimeError("marker-set QC 汇总表不是 9×9")
    (marker_patient_table, marker_summary_table, audit_table) = (
        marker_patient_scores_table,
        marker_summary_table,
        marker_audit_table,
    )

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = marker_patient_file
    marker_patient_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = marker_summary_file
    marker_summary_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = marker_audit_file
    audit_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 无论重算还是只绘图，患者marker长表须为972行，汇总须为81行。
if marker_patient_table.shape[0] != 972 or marker_summary_table.shape[0] != 81:
    raise RuntimeError("冻结 marker-set 源数据尺寸异常")
waterfall_plot_source_table = waterfall_table_total

waterfall_plot_stem = "Figure2_state_patient_waterfall"

waterfall_plot_source = source_out / f"{waterfall_plot_stem}_source.csv"

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = waterfall_plot_source
waterfall_plot_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
waterfall_plot_summary_table = waterfall_plot_source_table.drop_duplicates("state").set_index("state")

mm = 1 / 25.4

# 绘制原3×3患者waterfall，统一画布和患者槽位。
waterfall_plot_fig, waterfall_plot_axes = plt.subplots(3, 3, figsize=(183 * mm, 104 * mm), sharex=True, sharey=True)
waterfall_plot_y_min, waterfall_plot_y_max = -0.76, 0.56

# 瀑布图正负差值分别使用红蓝，状态标题和均值虚线使用固定九状态配色。
blue = "#3C5488"

negative = blue

red = "#E64B35"

positive = red

salmon = "#F39B7F"

state_colors = {
    "MCG1": "#00A087",
    "MCG2": "#91D1C2",
    "MCG3": "#4DBBD5",
    "MCG4": "#8491B4",
    "MCG5": blue,
    "MAC1": salmon,
    "MAC2": red,
    "E-MDSC": "#B09C85",
    "M-MDSC": "#7E6148",
}

# 每格按该状态内差值排序绘制患者柱；背景阴影区分高于和低于其余髓系。
for waterfall_plot_index, (waterfall_plot_ax, waterfall_plot_state) in enumerate(
    zip(waterfall_plot_axes.flat, state_order, strict=True)
):
    waterfall_plot_data_table = waterfall_plot_source_table.loc[
        waterfall_plot_source_table["state"].eq(waterfall_plot_state)
    ]
    waterfall_plot_x = waterfall_plot_data_table["within_state_rank"].to_numpy(int)
    waterfall_plot_y = waterfall_plot_data_table["raw20_delta"].to_numpy(float)
    waterfall_plot_colors = np.where(waterfall_plot_y > 0, positive, negative)
    waterfall_plot_ax.axhspan(0, waterfall_plot_y_max, color="#FFF7F4", zorder=0)

    waterfall_plot_ax.axhspan(waterfall_plot_y_min, 0, color="#F4F6FA", zorder=0)
    black = "#222222"
    waterfall_plot_ax.axhline(0, color=black, linewidth=0.65, zorder=2)
    waterfall_plot_ax.bar(
        waterfall_plot_x,
        waterfall_plot_y,
        width     = 0.76,
        color     = waterfall_plot_colors,
        edgecolor = "none",
        alpha     = 0.88,
        zorder    = 3,
    )

    # 虚线给出患者等权均值；不足6位患者的状态标题附加dagger标记。
    waterfall_plot_mean = float(waterfall_plot_summary_table.loc[waterfall_plot_state, "mean_delta"])
    waterfall_plot_ax.axhline(
        waterfall_plot_mean, color=state_colors[waterfall_plot_state], linewidth=1.05, linestyle=(0, (3, 2)), zorder=4
    )
    waterfall_plot_ax.set_xlim(0.25, 21.75)
    waterfall_plot_ax.set_ylim(waterfall_plot_y_min, waterfall_plot_y_max)
    waterfall_patient_count = int(waterfall_plot_summary_table.loc[waterfall_plot_state, "n_patients"])
    waterfall_small_sample_mark = "†" if waterfall_patient_count < 6 else ""
    waterfall_plot_state_title = f"{waterfall_plot_state}{waterfall_small_sample_mark}"
    waterfall_plot_ax.set_title(
        waterfall_plot_state_title, loc="left", color=state_colors[waterfall_plot_state], fontweight="bold", pad=7
    )

    waterfall_positive_count = int(waterfall_plot_summary_table.loc[waterfall_plot_state, "n_positive"])

    # FDR只按原三段精度格式化；显示规则不参与显著性判断。
    waterfall_fdr = float(waterfall_plot_summary_table.loc[waterfall_plot_state, "fdr_9_states"])
    if waterfall_fdr < 0.001:
        waterfall_fdr_label = f"{waterfall_fdr:.2e}"
    elif waterfall_fdr < 0.01:
        waterfall_fdr_label = f"{waterfall_fdr:.4f}"
    else:
        waterfall_fdr_label = f"{waterfall_fdr:.3f}"
    waterfall_plot_ax.text(
        1.0,
        1.025,
        (f"{waterfall_positive_count}/{waterfall_patient_count}  |  FDR {waterfall_fdr_label}"),
        transform = waterfall_plot_ax.transAxes,
        ha        = "right",
        va        = "bottom",
        fontsize  = 5.8,
        color     = black,
    )

    # MCG3保留虚线边框与score/pseudobulk方向不一致提示。
    grey = "#868686"
    if waterfall_plot_state == "MCG3":
        for waterfall_plot_spine in waterfall_plot_ax.spines.values():
            waterfall_plot_spine.set_color(grey)
            waterfall_plot_spine.set_linestyle((0, (3, 2)))
        waterfall_plot_ax.text(
            0.98,
            0.06,
            "score / pseudobulk\ndirection differs",
            transform = waterfall_plot_ax.transAxes,
            ha        = "right",
            va        = "bottom",
            fontsize  = 5.3,
            color     = grey,
        )
    else:
        waterfall_plot_ax.spines[["top", "right"]].set_visible(False)

    # 3×3面板只在左列显示纵轴、底行显示患者排序横轴，统一21个患者槽位。
    waterfall_plot_ax.spines[["left", "bottom"]].set_color(black)
    waterfall_plot_ax.tick_params(length=2.2, width=0.65)
    if waterfall_plot_index % 3 != 0:
        waterfall_plot_ax.tick_params(axis="y", left=False, labelleft=False)
    if waterfall_plot_index % 3 == 0:
        waterfall_plot_ax.set_ylabel("State − other myeloid\nMiller-IM score")
    if waterfall_plot_index >= 6:
        waterfall_plot_ax.set_xticks([1, 7, 14, 21])
        waterfall_plot_ax.set_xlabel("Patients ranked within state")
    else:
        waterfall_plot_ax.tick_params(axis="x", bottom=False, labelbottom=False)
waterfall_plot_fig.text(
    0.995,
    0.012,
    "† fewer than 6 patients; all panels use 21 fixed patient slots",
    ha       = "right",
    va       = "bottom",
    fontsize = 5.6,
    color    = grey,
)

# 固定3×3面板留白后，保存原tight裁边及0.035英寸边距的PDF。
waterfall_plot_fig.subplots_adjust(left=0.09, right=0.995, bottom=0.13, top=0.965, wspace=0.13, hspace=0.29)
waterfall_plot_pdf = fig_root / f"{waterfall_plot_stem}.pdf"
waterfall_plot_png = fig_root / f"{waterfall_plot_stem}.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
waterfall_plot_fig.canvas.draw()
waterfall_plot_fig.savefig(
    waterfall_plot_pdf,
    bbox_inches = "tight",
    pad_inches  = 0.035,
    metadata    = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None},
)

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(waterfall_plot_pdf, waterfall_plot_png, dpi=600)
plt.close(waterfall_plot_fig)
waterfall_record = {
    "panel": "2A",
    "source": str(waterfall_plot_source),
    "pdf": str(waterfall_plot_pdf),
    "png": str(waterfall_plot_png),
}

# 热图上块源表标注held-out基因，下块标注标签构建QC，分别保留源CSV。
heatmap_gene_summary_table = gene_summary_table
heatmap_marker_summary_table = marker_summary_table

heatmap_state_summary_table = state_summary_table

heatmap_stem = "Figure2_raw20_marker_double_heatmap"
heatmap_gene_source_table = heatmap_gene_summary_table.copy()
heatmap_gene_source_table["block"] = "Held-out Miller-IM genes; not used for state assignment"
heatmap_marker_source_table = heatmap_marker_summary_table.copy()
heatmap_marker_source_table["block"] = "Author marker-set assignment scores; classification QC"
heatmap_source_top = source_out / f"{heatmap_stem}_raw20_source.csv"

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = heatmap_source_top
heatmap_gene_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
heatmap_source_bottom = source_out / f"{heatmap_stem}_marker_qc_source.csv"

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = heatmap_source_bottom
heatmap_marker_source_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 上块保留held-out基因表达，下块保留作者marker标签构建QC；不交换两种证据。
heatmap_top_table = heatmap_gene_summary_table.pivot(index="gene", columns="state", values="within_gene_z").reindex(
    index= measurable_raw20, columns=state_order
)
heatmap_bottom_table = heatmap_marker_summary_table.pivot(
    index= "tested_marker_set", columns="assigned_state", values="within_marker_set_z"
).reindex(index=state_order, columns=state_order)
if heatmap_top_table.isna().any().any() or heatmap_bottom_table.isna().any().any():
    raise RuntimeError("双块热图存在缺失单元格")

heatmap_n_by_state = heatmap_state_summary_table.set_index("state")["n_patients"].astype(int).to_dict()

# 列标题显示各状态可评估患者数；少于6位患者继续标dagger。
heatmap_state_labels = [
    f"{heatmap_state}\n(n={heatmap_n_by_state[heatmap_state]}){'†' if heatmap_n_by_state[heatmap_state] < 6 else ''}"
    for heatmap_state in state_order
]

# 四行布局依次放上块标题、19基因热图、证据分隔说明与9套marker热图。
heatmap_fig = plt.figure(figsize=(183 * mm, 148 * mm))
heatmap_grid = heatmap_fig.add_gridspec(
    nrows         = 4,
    ncols         = 2,
    width_ratios  = [1, 0.034],
    height_ratios = [2.2, 19, 1.35, 5.3],
    left          = 0.16,
    right         = 0.94,
    bottom        = 0.10,
    top           = 0.98,
    wspace        = 0.16,
    hspace        = 0.05,
)

heatmap_title_top = heatmap_fig.add_subplot(heatmap_grid[0, 0])
heatmap_ax_top = heatmap_fig.add_subplot(heatmap_grid[1, 0])
heatmap_separator = heatmap_fig.add_subplot(heatmap_grid[2, 0])
heatmap_ax_bottom = heatmap_fig.add_subplot(heatmap_grid[3, 0])
heatmap_cax_top = heatmap_fig.add_subplot(heatmap_grid[1, 1])
heatmap_cax_bottom = heatmap_fig.add_subplot(heatmap_grid[3, 1])

heatmap_title_top.axis("off")
heatmap_title_top.text(
    0,
    0.93,
    "Held-out Miller-IM genes — not used for state assignment",
    ha         = "left",
    va         = "center",
    fontsize   = 7.4,
    fontweight = "bold",
    color      = black,
)

# 上块使用蓝白鲑红色带，色限固定为±2.5，展示within_gene_z。
raw20_cmap = LinearSegmentedColormap.from_list("raw20_npg", [blue, "#F7F7F7", salmon])

sns.heatmap(
    heatmap_top_table,
    ax          = heatmap_ax_top,
    cmap        = raw20_cmap,
    vmin        = -2.5,
    vmax        = 2.5,
    center      = 0,
    linewidths  = 0.45,
    linecolor   = "white",
    cbar        = True,
    cbar_ax     = heatmap_cax_top,
    xticklabels = heatmap_state_labels,
    yticklabels = measurable_raw20,
)

# 上块状态标签置顶，基因标签置左，独立色条说明基因内标准化。
heatmap_ax_top.xaxis.tick_top()
heatmap_ax_top.tick_params(axis="x", rotation=0, length=0, pad=2)
heatmap_ax_top.tick_params(axis="y", rotation=0, length=0, pad=2)
heatmap_ax_top.set_xlabel("")
heatmap_ax_top.set_ylabel("")
heatmap_cax_top.set_ylabel("Within-gene z score", fontsize=6.3, labelpad=4)
heatmap_cax_top.tick_params(labelsize=5.7, length=2)

# 分隔带说明下块属于标签构建QC，并显示已在输入处核对的零基因重叠。
heatmap_separator.axis("off")
heatmap_separator.axhline(0.92, color=black, linewidth=0.8)

heatmap_separator.text(
    0,
    0.55,
    "Author marker-set assignment scores — used to define labels; classification QC",
    ha         = "left",
    va         = "center",
    fontsize   = 6.8,
    fontweight = "bold",
    color      = grey,
)

heatmap_separator.text(
    0, 0.08, "Direct marker–Miller-IM gene overlap after exclusion: 0", ha="left", va="center", fontsize=5.9, color=grey
)

# 下块使用灰色色带和同样±2.5色限，显示marker集内标准化分数。
qc_cmap = LinearSegmentedColormap.from_list("marker_qc", ["#ECECEC", "#A9B3BE", "#4B5563"])

sns.heatmap(
    heatmap_bottom_table,
    ax          = heatmap_ax_bottom,
    cmap        = qc_cmap,
    vmin        = -2.5,
    vmax        = 2.5,
    center      = 0,
    linewidths  = 0.45,
    linecolor   = "white",
    cbar        = True,
    cbar_ax     = heatmap_cax_bottom,
    xticklabels = False,
    yticklabels = [f"{heatmap_state} markers" for heatmap_state in state_order],
)

# 下块仅保留marker行标签与灰色色条，底部说明其不属于独立验证。
heatmap_ax_bottom.tick_params(axis="y", rotation=0, length=0, pad=2, colors=grey)
heatmap_ax_bottom.tick_params(axis="x", bottom=False)
heatmap_ax_bottom.set_xlabel("")
heatmap_ax_bottom.set_ylabel("")
heatmap_cax_bottom.set_ylabel("Within-marker-set z score", fontsize=5.8, labelpad=4, color=grey)
heatmap_cax_bottom.tick_params(labelsize=5.5, length=2, colors=grey)
for heatmap_spine in heatmap_cax_bottom.spines.values():
    heatmap_spine.set_edgecolor(grey)

heatmap_fig.text(
    0.16,
    0.035,
    (
        "Columns are fixed a priori; † fewer than 6 patients. The lower block documents l"
        "abel construction and is not independent validation."
    ),
    ha       = "left",
    va       = "bottom",
    fontsize = 5.7,
    color    = grey,
)

heatmap_pdf = fig_root / f"{heatmap_stem}.pdf"
heatmap_png = fig_root / f"{heatmap_stem}.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
heatmap_fig.canvas.draw()
heatmap_fig.savefig(
    heatmap_pdf,
    bbox_inches = "tight",
    pad_inches  = 0.035,
    metadata    = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None},
)

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(heatmap_pdf, heatmap_png, dpi=600)
plt.close(heatmap_fig)
heatmap_record = {
    "panel": "2B",
    "source": f"{heatmap_source_top};{heatmap_source_bottom}",
    "pdf": str(heatmap_pdf),
    "png": str(heatmap_png),
}
records = [waterfall_record, heatmap_record]
preview_records = records

# PIL仅用于两张已导出面板的预览拼版；正式面板PNG来自PDF。
preview_images = [Image.open(preview_record["png"]).convert("RGB") for preview_record in preview_records]
preview_target_width = 1800
preview_resized: list[Image.Image] = []
for preview_image in preview_images:
    preview_height = round(preview_image.height * preview_target_width / preview_image.width)
    preview_resized.append(preview_image.resize((preview_target_width, preview_height), Image.Resampling.LANCZOS))

preview_margin = 85
preview_gap = 40
preview_canvas = Image.new(
    "RGB",
    (
        preview_target_width + preview_margin * 2,
        sum(preview_image.height for preview_image in preview_resized) + preview_gap + preview_margin * 2,
    ),
    "white",
)

preview_draw = ImageDraw.Draw(preview_canvas)

# 预览沿用原Helvetica、Arial候选顺序；与共享image_font的顺序不同，故保留。
preview_font_candidates = ["/System/Library/Fonts/Helvetica.ttc", "/Library/Fonts/Arial.ttf"]
preview_font_path = next(
    (Path(preview_path) for preview_path in preview_font_candidates if Path(preview_path).exists()), None
)
preview_font = ImageFont.truetype(str(preview_font_path), 54) if preview_font_path else ImageFont.load_default()
preview_y = preview_margin
for preview_label, preview_image in zip(["A", "B"], preview_resized, strict=True):
    preview_draw.text((22, preview_y + 2), preview_label, fill=black, font=preview_font)
    preview_canvas.paste(preview_image, (preview_margin, preview_y))
    preview_y += preview_image.height + preview_gap

# 保存彩色和灰度PIL拼版供预览；两张正式面板的PNG均由各自PDF生成。
preview_preview = fig_root / "Figure2_two_panel_preview.png"
preview_canvas.save(preview_preview, dpi=(150, 150))
preview_grayscale = fig_root / "Figure2_two_panel_preview_grayscale.png"
preview_canvas.convert("L").convert("RGB").save(preview_grayscale, dpi=(150, 150))
(preview, grayscale) = (preview_preview, preview_grayscale)

# 保留原图注文本；这里只改变源码换行，不修改科学表述。
legend_path = write_root / "Figure2_legend.md"
legend_path.write_text(
    (
        "# Figure 2 | The Miller-IM program preferentially maps to MCG1/MCG2-like states "
        "in an independent cross-sectional myeloid atlas\n\n**A.** Patient-level state-vers"
        "us-rest waterfall plots in the 21 IDH-wild-type primary GBM samples of GSE278456"
        ". Each bar is one patient with at least 30 cells in the indicated author-marker-"
        "derived state and at least 30 cells in the remaining prespecified myeloid states"
        ". Bars show the difference in control-adjusted Miller-IM score between that stat"
        "e and the remaining states; dashed lines show patient-equal means. Red and blue "
        "denote values above and below the remaining states, respectively. Numbers give p"
        "atients above the remaining states / evaluable patients and Benjamini–Hochberg F"
        "DR across the nine state-wise exact two-sided sign-flip tests. All facets use th"
        "e same 21 patient slots and y-axis. MCG3 is marked because its mean score direct"
        "ion differs from the raw-count pseudobulk sensitivity analysis.\n\n**B.** Two-bloc"
        "k patient-balanced heatmap across the same fixed state order. The upper block sh"
        "ows relative expression of the 19 measurable Miller-IM genes; values are patient"
        "-equal mean log-normalized expression standardized within each gene across the n"
        "ine states. All measurable Miller-IM genes were removed from the author marker s"
        "ets before state assignment, so this block was held out from direct state-label "
        "construction. The lower, visually reduced block shows the nine complete author m"
        "arker-set scores, averaged within patient and assigned state and then standardiz"
        "ed within marker set. These scores were used to create the state labels; their d"
        "iagonal structure is expected by construction and is shown only as classificatio"
        "n quality control, not as independent validation. Direct overlap between the fin"
        "al label-defining marker sets and the measurable Miller-IM genes was zero. State"
        " labels show the number of evaluable patients; dagger marks fewer than six patie"
        "nts.\n\nGSE278456 is cross-sectional and is used only to localize the recurrence-a"
        "ssociated program across myeloid states; recurrence evidence is provided separat"
        "ely in Figure 1.\n"
    ),
    encoding= "utf-8",
)

# 面板清单保留panel/source/pdf/png四字段，供报告与下游合版定位。
manifest = write_root / "Figure2_panel_manifest.csv"
panel_manifest_table = pd.DataFrame(records)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = manifest
panel_manifest_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
report_marker_patient_table = marker_patient_table
report_marker_summary_table = marker_summary_table
report_audit_table = audit_table
report_panel_table = pd.DataFrame(records)
report_preview = preview
report_grayscale = grayscale

# 用原结果字段生成既有说明文件；文案和输出名维持公开合同。
report_path = write_root / "FINAL_RESULT.md"
report_path.write_text(
    (
        "# Step50 Figure2 两面板候选版\n\n- 2A：9 状态患者级 waterfall，固定 21 个患者槽位和统一 y 轴。\n- 2B 上块：19×9"
        " held-out Miller-IM genes 患者等权 `within_gene_z`。\n- 2B 下块：9×9 完整 author marker-set"
        " score QC，先在 17 套 marker 间完成赋值，再截取 9 个髓系状态。\n- H5 复核："
        f"{int(report_audit_table['n_qc_pass'].sum()):,}"
        "/"
        f"{int(report_audit_table['n_qc_pass'].sum()):,}"
        " 个 QC 细胞的 winner label 与 Step48 精确一致；最大 winner-score 差 "
        f"{report_audit_table['max_abs_winner_score_difference'].max():.3g}"
        "。\n- 源数据：waterfall 108 行；Miller-IM 上块 171 行；marker 患者长表 "
        f"{report_marker_patient_table.shape[0]}"
        " 行；marker QC 汇总 "
        f"{report_marker_summary_table.shape[0]}"
        " 行。\n- 预览：`"
        f"{report_preview}"
        "`；灰度预览：`"
        f"{report_grayscale}"
        "`。\n- 边界：GSE278456 只做横断面状态定位；下块只是标签构建 QC，不是独立验证。\n\n## 面板\n\n"
        f"{report_panel_table.to_markdown(index=False)}"
        "\n"
    ),
    encoding= "utf-8",
)

# 原输出摘要清单仍保留给既有复现检查使用，不新增其他审计产物。
inventory = []
for path in sorted(list(fig_root.glob("*.pdf")) + list(fig_root.glob("*.png"))):
    inventory.append({"path": str(path), "sha256": sha256(path)})
inventory_table = pd.DataFrame(inventory)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = write_root / "Figure2_output_hashes.csv"
inventory_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")
print("Figure2 panels:", len(records))

# 核对本次图注、结果说明和两张PIL预览存在；预览还须可解码且尺寸为正。
for output_path in [legend_path, report_path, preview, grayscale]:
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
for preview_path in [preview, grayscale]:
    with Image.open(preview_path) as preview_image:
        if min(preview_image.size) <= 0:
            raise RuntimeError(f"Invalid preview dimensions: {preview_path.name}")
        preview_image.verify()
