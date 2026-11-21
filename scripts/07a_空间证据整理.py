#!/usr/bin/env python3
# 目的：按原公式重算GeoMx配对raw20和Visium空间相关，生成Figure3与SupplementaryFigure6。
# 输入：外部数据根的GeoMx metadata.csv/normalized.txt、GSE276841 H5/位置表及作者marker工作簿；Step38 GSEA表。
# 输出：write/39_independent_miller_mg_inflammatory_figures内Figure3/SupplementaryFigure6/source_csv及汇总；
# figures/39_independent_miller_mg_inflammatory_plottie内同名面板PDF与其PNG。
# 运行：在项目根使用既有sc环境执行 python scripts/07a_空间证据整理.py。
# 保留Python：Scanpy、AnnData及Matplotlib是已公开复现路径；改成R超出本轮纯风格整理。
# 存量source_data/config、数据文件名、输出字段、write及figures路径保留原公开合同。

from __future__ import annotations
import itertools
import os
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

# 沿用原图的字体、线宽和背景配置，以保持PDF画面。
mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#222222",
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.transparent": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
)

sns.set_style("white")

# 创建Figure3、SupplementaryFigure6及各自源表目录，保留原公开路径。
root = Path(__file__).resolve().parents[1]

fig_base = root / "figures" / "39_independent_miller_mg_inflammatory_plottie"

fig3_fig = fig_base / "Figure3"

write_base = root / "write" / "39_independent_miller_mg_inflammatory_figures"

fig3_write = write_base / "Figure3"

s6_fig = fig_base / "SupplementaryFigure6"

s6_write = write_base / "SupplementaryFigure6"

for directories_path in [fig3_write, s6_write, fig3_fig, s6_fig]:
    directories_path.mkdir(parents=True, exist_ok=True)
for directories_path in [fig3_write / "source_csv", s6_write / "source_csv"]:
    directories_path.mkdir(parents=True, exist_ok=True)

# 外部数据根读取作者marker工作簿、GeoMx表达和Visium原始输入。
ext_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()

marker_xls = (
    ext_root / "write" / "30_dataset_rescue_search" / "source_metadata" / "NIHMS2115262-supplement-S2_S5_S6.xls"
)

author_workbook = pd.ExcelFile(marker_xls)
author_markers: dict[str, list[str]] = {}

# Miller原始排序前20基因固定不变；基因大写是外部数据标识。
raw20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]

# 先筛FDR<0.05且Foldchange>0，再去除raw20，最后取前40个；此顺序与Step48不同。
author_excluded = set(raw20)
for author_sheet in author_workbook.sheet_names:
    author_data_table = pd.read_excel(marker_xls, sheet_name=author_sheet)
    if "FDR" not in author_data_table.columns:
        author_data_table = pd.read_excel(marker_xls, sheet_name=author_sheet, header=1)
    author_gene_col = author_data_table.columns[0]
    author_chosen = author_data_table.loc[
        author_data_table["FDR"].lt(0.05) & author_data_table["Foldchange"].gt(0), author_gene_col
    ].astype(str)
    author_markers[author_sheet] = [author_gene for author_gene in author_chosen if author_gene not in author_excluded][
        :40
    ]

markers = author_markers

geomx_source = ext_root / "write" / "26_pure_bioinformatics_dataset_rescue" / "source_metadata" / "GeoMx_zenodo16839828"

# 读取GeoMx区域元数据及表达矩阵；区域标识必须完全匹配。
geomx_metadata_table = pd.read_csv(geomx_source / "metadata.csv")
geomx_expression_table = pd.read_csv(geomx_source / "normalized.txt", sep="\t", index_col=0)
if set(geomx_expression_table.columns) != set(geomx_metadata_table["Unnamed: 0"]):
    raise ValueError("GeoMx matrix and metadata do not match.")

geomx_raw20_present = [geomx_gene for geomx_gene in raw20 if geomx_gene in geomx_expression_table.index]
geomx_raw20_missing = [geomx_gene for geomx_gene in raw20 if geomx_gene not in geomx_expression_table.index]
geomx_coverage_table = pd.DataFrame(
    {
        "signature": ["Miller_raw20"],
        "n_defined": [len(raw20)],
        "n_present": [len(geomx_raw20_present)],
        "coverage": [len(geomx_raw20_present) / len(raw20)],
        "present_genes": [";".join(geomx_raw20_present)],
        "missing_genes": [";".join(geomx_raw20_missing)],
    }
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "geomx_raw20_coverage.csv"
geomx_coverage_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

geomx_metadata_table = geomx_metadata_table.copy()

# 保留严格IBA1+、QC通过、IDH_WT入口及两套原敏感性入口。
geomx_metadata_table["entry_strict_pass_idhwt"] = (
    geomx_metadata_table["iba1"].eq(True)
    & geomx_metadata_table["QCFlags"].isna()
    & geomx_metadata_table["IDH_status"].eq("IDH_WT")
)
geomx_metadata_table["entry_strict_pass_all_idh"] = (
    geomx_metadata_table["iba1"].eq(True) & geomx_metadata_table["QCFlags"].isna()
)
geomx_metadata_table["entry_all_aoi_idhwt"] = geomx_metadata_table["iba1"].eq(True) & geomx_metadata_table[
    "IDH_status"
].eq("IDH_WT")

geomx_pair_stats = []
geomx_paired_scores = []
geomx_gene_stats = []
geomx_patient_score_tables = []

geomx_subgroup_stats = []
geomx_patient_map_rows = []

# 按严格IDH-wt、严格全IDH、全部IBA1+ IDH-wt顺序处理三种既定入口。
geomx_entries = {
    "strict_pass_idhwt": "entry_strict_pass_idhwt",
    "strict_pass_all_idh": "entry_strict_pass_all_idh",
    "all_aoi_idhwt": "entry_all_aoi_idhwt",
}

# 每入口仅取Primary和Recurrence区域，对表达作log2(x+1)后汇总患者时间点。
for geomx_entry_name, geomx_column in geomx_entries.items():
    geomx_current_table = geomx_metadata_table.loc[geomx_metadata_table[geomx_column]].copy()
    geomx_current_table = geomx_current_table.loc[
        geomx_current_table["tumor_setting"].isin(["Primary", "Recurrence"])
    ].copy()
    geomx_ids = geomx_current_table["Unnamed: 0"].tolist()
    geomx_log_expression_table = np.log2(geomx_expression_table.loc[:, geomx_ids] + 1.0)

    # 先按患者、治疗组、时间点及IDH状态平均区域表达，再计算原逐基因z分数。
    geomx_patient_gene_table = (
        geomx_log_expression_table.T.assign(
            Patient_No    = geomx_current_table.set_index("Unnamed: 0")
            .loc[geomx_log_expression_table.columns, "Patient_No"]
            .to_numpy(),
            trial_setting = geomx_current_table.set_index("Unnamed: 0")
            .loc[geomx_log_expression_table.columns, "trial_setting"]
            .to_numpy(),
            tumor_setting = geomx_current_table.set_index("Unnamed: 0")
            .loc[geomx_log_expression_table.columns, "tumor_setting"]
            .to_numpy(),
            IDH_status    = geomx_current_table.set_index("Unnamed: 0")
            .loc[geomx_log_expression_table.columns, "IDH_status"]
            .to_numpy(),
        )
        .groupby(["Patient_No", "trial_setting", "tumor_setting", "IDH_status"], observed=True)[geomx_raw20_present]
        .mean()
        .reset_index()
    )

    # 在当前入口的患者时间点间按基因计算总体z分数，再对可测raw20基因取均值。
    geomx_z_table = geomx_patient_gene_table[geomx_raw20_present].copy()
    geomx_z_table = (geomx_z_table - geomx_z_table.mean(axis=0)) / geomx_z_table.std(axis=0, ddof=0).replace(0, np.nan)
    geomx_patient_gene_table["raw20_score"] = geomx_z_table.mean(axis=1)
    geomx_patient_score_tables.append(geomx_patient_gene_table.copy())

    # 按同一患者和治疗/IDH分层配对Primary与Recurrence；缺任一端点的行仍排除。
    geomx_paired_table = (
        geomx_patient_gene_table.pivot_table(
            index   = ["Patient_No", "trial_setting", "IDH_status"],
            columns = "tumor_setting",
            values  = "raw20_score",
            aggfunc = "first",
        )
        .dropna(subset=["Primary", "Recurrence"])
        .reset_index()
    )

    geomx_paired_table["delta"] = geomx_paired_table["Recurrence"] - geomx_paired_table["Primary"]
    geomx_paired_table["entry"] = geomx_entry_name
    geomx_paired_scores.append(geomx_paired_table.copy())

    def sign_flip_p(delta: np.ndarray, seed: int = 20260713) -> float:
        """输入患者差值与原种子；返回原精确或300,000次抽样符号翻转P值，保留加一校正。"""

        delta = np.asarray(delta, dtype=float)
        delta = delta[np.isfinite(delta)]
        n = len(delta)

        # 无有限差值返回NaN；20对以内穷举符号，大样本保留300000次随机符号翻转。
        if n == 0:
            result = np.nan
            return result
        observed = abs(delta.mean())
        if n <= 20:
            signs = np.asarray(list(itertools.product([-1.0, 1.0], repeat=n)))
            null = np.abs((signs * delta).mean(axis=1))
        else:
            rng = np.random.default_rng(seed)
            signs = rng.choice(np.array([-1.0, 1.0]), size=(300_000, n))

            null = np.abs((signs * delta).mean(axis=1))
        result = float((np.count_nonzero(null >= observed - 1e-15) + 1) / (len(null) + 1))
        return result

    geomx_pair_stats.append(
        {
            "entry": geomx_entry_name,
            "n_pairs": int(len(geomx_paired_table)),
            "n_positive": int((geomx_paired_table["delta"] > 0).sum()),
            "positive_fraction": float((geomx_paired_table["delta"] > 0).mean()) if len(geomx_paired_table) else np.nan,
            "mean_delta": float(geomx_paired_table["delta"].mean()) if len(geomx_paired_table) else np.nan,
            "median_delta": float(geomx_paired_table["delta"].median()) if len(geomx_paired_table) else np.nan,
            "sign_flip_p": sign_flip_p(geomx_paired_table["delta"].to_numpy()),
        }
    )

    geomx_patient_map_rows.append(geomx_paired_table.copy())

    # 分别保留总体、Control、Nivolumab及两类IDH分层；空分层统计仍为NaN。
    for geomx_subgroup, geomx_data_table in {
        "All": geomx_paired_table,
        "Control": geomx_paired_table.loc[geomx_paired_table["trial_setting"].eq("Control")],
        "Nivolumab": geomx_paired_table.loc[geomx_paired_table["trial_setting"].eq("Nivolumab")],
        "IDH_WT": geomx_paired_table.loc[geomx_paired_table["IDH_status"].eq("IDH_WT")],
        "IDH_mutated": geomx_paired_table.loc[geomx_paired_table["IDH_status"].eq("IDH_mutated")],
    }.items():
        geomx_subgroup_stats.append(
            {
                "entry": geomx_entry_name,
                "subgroup": geomx_subgroup,
                "n_pairs": int(len(geomx_data_table)),
                "mean_delta": float(geomx_data_table["delta"].mean()) if len(geomx_data_table) else np.nan,
                "median_delta": float(geomx_data_table["delta"].median()) if len(geomx_data_table) else np.nan,
                "positive_fraction": float((geomx_data_table["delta"] > 0).mean()) if len(geomx_data_table) else np.nan,
                "sign_flip_p": sign_flip_p(geomx_data_table["delta"].to_numpy()) if len(geomx_data_table) else np.nan,
            }
        )

    # 治疗交互用Nivolumab减Control的患者配对差值均值，检验采用标签置换。
    geomx_control = geomx_paired_table.loc[geomx_paired_table["trial_setting"].eq("Control"), "delta"].to_numpy()
    geomx_nivolumab = geomx_paired_table.loc[geomx_paired_table["trial_setting"].eq("Nivolumab"), "delta"].to_numpy()

    def label_permutation_p(control: np.ndarray, treated: np.ndarray, seed: int = 20260713) -> float:
        """输入Control和Nivolumab差值及原种子；返回300,000次标签置换P值，保留原调用顺序。"""

        control = np.asarray(control, dtype=float)
        treated = np.asarray(treated, dtype=float)
        observed = abs(treated.mean() - control.mean())
        joined = np.concatenate([control, treated])

        # 每次置换合并后的差值并保持原Control人数；原种子和300000次循环不变。
        n_control = len(control)
        rng = np.random.default_rng(seed)
        null = np.empty(300_000, dtype=float)
        for idx in range(len(null)):
            perm = rng.permutation(joined)
            null[idx] = abs(perm[n_control:].mean() - perm[:n_control].mean())
        result = float((np.count_nonzero(null >= observed - 1e-15) + 1) / (len(null) + 1))
        return result

    geomx_subgroup_stats.append(
        {
            "entry": geomx_entry_name,
            "subgroup": "Interaction_Nivolumab_minus_Control",
            "n_pairs": int(len(geomx_paired_table)),
            "mean_delta": float(geomx_nivolumab.mean() - geomx_control.mean()),
            "median_delta": np.nan,
            "positive_fraction": np.nan,
            "sign_flip_p": label_permutation_p(geomx_control, geomx_nivolumab),
        }
    )

    # 基因层面沿用相同患者时间点配对；区间和检验均作用于log2表达差值。
    geomx_gene_patient_table = geomx_patient_gene_table[
        ["Patient_No", "trial_setting", "tumor_setting", "IDH_status"] + geomx_raw20_present
    ].copy()
    for geomx_gene in geomx_raw20_present:

        # 每次仅处理一个可测基因；沿用原配对顺序、均值区间和符号翻转检验。
        geomx_gene_paired_table = (
            geomx_gene_patient_table.pivot_table(
                index   = ["Patient_No", "trial_setting", "IDH_status"],
                columns = "tumor_setting",
                values  = geomx_gene,
                aggfunc = "first",
            )
            .dropna(subset=["Primary", "Recurrence"])
            .reset_index()
        )

        geomx_gene_paired_table["delta"] = geomx_gene_paired_table["Recurrence"] - geomx_gene_paired_table["Primary"]

        def mean_ci(values: np.ndarray) -> tuple[float, float, float]:
            """输入差值向量；返回原均值和1.96倍标准误区间，空输入返回NaN，单值区间等于该值。"""

            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]
            if len(values) == 0:
                result = (np.nan, np.nan, np.nan)
                return result

            mean = float(values.mean())
            if len(values) == 1:
                result = (mean, mean, mean)
                return result
            sem = float(stats.sem(values))
            delta = 1.96 * sem
            result = (mean, mean - delta, mean + delta)
            return result

        geomx_mean, geomx_ci_low, geomx_ci_high = mean_ci(geomx_gene_paired_table["delta"].to_numpy())

        geomx_gene_stats.append(
            {
                "entry": geomx_entry_name,
                "gene": geomx_gene,
                "n_pairs": int(len(geomx_gene_paired_table)),
                "mean_log2_delta": geomx_mean,
                "ci_low": geomx_ci_low,
                "ci_high": geomx_ci_high,
                "median_log2_delta": float(geomx_gene_paired_table["delta"].median()),
                "positive_fraction": float((geomx_gene_paired_table["delta"] > 0).mean()),
                "sign_flip_p": sign_flip_p(geomx_gene_paired_table["delta"].to_numpy()),
            }
        )

# 构建入口、分层和基因三类统计表，分别保留原BH校正范围。
geomx_pair_stats_table = pd.DataFrame(geomx_pair_stats)

def bh(series: pd.Series) -> np.ndarray:
    """输入含缺失的P值序列；返回原顺序的BH校正值，缺失保持NaN。"""

    out = np.full(len(series), np.nan)
    ok = series.notna().to_numpy()
    if ok.any():
        out[ok] = multipletests(series.loc[ok], method="fdr_bh")[1]
    return out

# 入口表整体BH；分层及基因表均在各入口内部BH，最后拼接患者记录。
geomx_pair_stats_table["fdr"] = bh(geomx_pair_stats_table["sign_flip_p"])
geomx_subgroup_table = pd.DataFrame(geomx_subgroup_stats)
geomx_subgroup_table["fdr"] = geomx_subgroup_table.groupby("entry")["sign_flip_p"].transform(bh)
geomx_gene_stats_table = pd.DataFrame(geomx_gene_stats)
geomx_gene_stats_table["fdr"] = geomx_gene_stats_table.groupby("entry")["sign_flip_p"].transform(bh)
geomx_patient_scores_table_total = pd.concat(geomx_patient_score_tables, ignore_index=True)
geomx_paired_scores_table_total = pd.concat(geomx_paired_scores, ignore_index=True)
geomx_patient_map_table_total = pd.concat(geomx_patient_map_rows, ignore_index=True)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "geomx_raw20_entry_level_summary.csv"
geomx_pair_stats_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "geomx_raw20_subgroup_summary.csv"
geomx_subgroup_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "geomx_raw20_gene_level_summary.csv"
geomx_gene_stats_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "geomx_raw20_patient_timepoint_scores.csv"
geomx_patient_scores_table_total.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "geomx_raw20_paired_deltas.csv"
geomx_paired_scores_table_total.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "geomx_raw20_patient_map.csv"
geomx_patient_map_table_total.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

geomx = {
    "coverage": geomx_coverage_table,
    "entry_summary": geomx_pair_stats_table,
    "subgroup_summary": geomx_subgroup_table,
    "gene_summary": geomx_gene_stats_table,
    "patient_scores": geomx_patient_scores_table_total,
    "paired_scores": geomx_paired_scores_table_total,
}

# 空间评分复用作者marker集合，并去除MDSC及MES集合中与raw20直接重叠的基因。
visium_markers = markers

gse276841_manifest = root / "write" / "30_dataset_rescue_search" / "GSE276841_sample_manifest.csv"

visium_manifest_table = pd.read_csv(gse276841_manifest).set_index("geo_accession")
visium_spot_rows = []
visium_partial_rows = []
visium_detected_counts: dict[str, int] = {}

visium_mmdsc_genes = [visium_gene for visium_gene in visium_markers["M-MDSC"] if visium_gene not in raw20]
visium_emdsc_genes = [visium_gene for visium_gene in visium_markers["E-MDSC"] if visium_gene not in raw20]

# MES与总髓系基因集保持原固定定义；基因大写为外部标识。
mes = [
    "CHI3L1", "CD44", "SERPINE1", "VEGFA", "ADM", "CA9", "BNIP3", "NDRG1", "SLC2A1", "HK2", "LDHA", "ENO1",
    "VIM", "LGALS3", "ANXA1", "TGFBI",
]

visium_mes_genes = [visium_gene for visium_gene in mes if visium_gene not in raw20]

spatial_source = ext_root / "write" / "30_dataset_rescue_search" / "source_metadata" / "GSE276841_spatial"

# 按文件名排序读取Visium矩阵，并记录每切片raw20可测基因数。
for visium_matrix_path in sorted(spatial_source.glob("*_filtered_feature_bc_matrix.h5")):
    visium_accession = visium_matrix_path.name.split("_")[0]

    visium_sample = visium_manifest_table.loc[visium_accession, "sample_title"]
    visium_prefix = visium_matrix_path.name.replace("_filtered_feature_bc_matrix.h5", "")

    visium_adata = sc.read_10x_h5(visium_matrix_path)
    visium_adata.var_names_make_unique()
    visium_detected_counts[visium_sample] = int(sum(visium_gene in visium_adata.var_names for visium_gene in raw20))

    # 按原10,000计数尺度归一化，随后log1p；不改变原始计数副本。
    sc.pp.normalize_total(visium_adata, target_sum=10_000)
    sc.pp.log1p(visium_adata)

    visium_scores_table = pd.DataFrame(index=visium_adata.obs_names)

    def control_gene_score(adata, genes: list[str], score_name: str) -> np.ndarray:
        """输入AnnData、基因及临时列名；返回Scanpy对照校正分数，并删除临时obs列；无可测基因返回NaN。"""

        use = [gene for gene in genes if gene in adata.var_names]
        if not use:
            result = np.full(adata.n_obs, np.nan)
            return result

        # 保留ctrl_size=50、n_bins=25和原种子；每套评分重置种子是原方法，移出会改变结果。
        sc.tl.score_genes(
            adata, gene_list=use, score_name=score_name, ctrl_size=50, n_bins=25, random_state=20260713, use_raw=False
        )
        values = adata.obs[score_name].to_numpy()
        del adata.obs[score_name]
        return values

    visium_scores_table["raw20"] = control_gene_score(visium_adata, raw20, "raw20_score")
    myeloid = ["PTPRC", "TYROBP", "AIF1", "LST1", "FCER1G", "CTSS", "CSF1R"]
    visium_scores_table["myeloid"] = control_gene_score(visium_adata, myeloid, "myeloid_score")
    visium_scores_table["mmdsc"] = control_gene_score(visium_adata, visium_mmdsc_genes, "mmdsc_score")

    # 两类MDSC分数按行平均形成组合分数；MES仍单独进行对照基因评分。
    visium_scores_table["emdsc"] = control_gene_score(visium_adata, visium_emdsc_genes, "emdsc_score")
    visium_scores_table["mdsc_combined"] = visium_scores_table[["mmdsc", "emdsc"]].mean(axis=1)
    visium_scores_table["mes"] = control_gene_score(visium_adata, visium_mes_genes, "mes_score")
    visium_scores_table["sample"] = visium_sample
    visium_scores_table["geo_accession"] = visium_accession

    # 按barcode加入组织位置，只保留in_tissue=1；纵坐标仍按原定义取负值。
    visium_positions_table = pd.read_csv(spatial_source / f"{visium_prefix}_tissue_positions.csv.gz").set_index(
        "barcode"
    )

    visium_scores_table = visium_scores_table.join(visium_positions_table, how="left")

    visium_scores_table = visium_scores_table.loc[visium_scores_table["in_tissue"].eq(1)].copy()
    visium_scores_table["pixel_y_plot"] = -visium_scores_table["pxl_row_in_fullres"]
    visium_spot_rows.append(visium_scores_table.reset_index(names="barcode"))

    # 逐切片计算MDSC、MES及两类MDSC与raw20的髓系校正秩相关。
    for visium_target_key, visium_target_label in [
        ("mdsc_combined", "MDSC"),
        ("mes", "MES"),
        ("emdsc", "E-MDSC"),
        ("mmdsc", "M-MDSC"),
    ]:

        def partial_spearman(x: pd.Series, y: pd.Series, covariate: pd.Series) -> tuple[float, float, int]:
            """输入两个分数及髓系协变量；返回秩残差Pearson相关、P值和有效点数；保留原少于20点门槛。"""

            # 仅保留三列共同非缺失点；少于20点或协变量秩零方差时返回不可估计。
            data_table_total = pd.concat([x, y, covariate], axis=1).dropna()
            if len(data_table_total) < 20:
                result = (np.nan, np.nan, len(data_table_total))
                return result
            ranked = data_table_total.rank(method="average").to_numpy(dtype=np.float64, copy=True)
            x_rank = ranked[:, 0] - ranked[:, 0].mean()
            y_rank = ranked[:, 1] - ranked[:, 1].mean()
            cov_rank = ranked[:, 2] - ranked[:, 2].mean()
            denominator = float(np.dot(cov_rank, cov_rank))
            if denominator <= 0:
                result = (np.nan, np.nan, len(data_table_total))
                return result

            # 分别去除髓系秩的线性投影，再对两组秩残差计算Pearson相关。
            rx = x_rank - cov_rank * (float(np.dot(cov_rank, x_rank)) / denominator)
            ry = y_rank - cov_rank * (float(np.dot(cov_rank, y_rank)) / denominator)
            rho, p_value = stats.pearsonr(rx, ry)
            result = (float(rho), float(p_value), len(data_table_total))
            return result

        visium_rho, visium_p_value, visium_n_spots = partial_spearman(
            visium_scores_table["raw20"], visium_scores_table[visium_target_key], visium_scores_table["myeloid"]
        )

        visium_partial_rows.append(
            {
                "sample": visium_sample,
                "geo_accession": visium_accession,
                "target": visium_target_label,
                "n_spots": visium_n_spots,
                "partial_rho": visium_rho,
                "p_value": visium_p_value,
            }
        )

# 合并切片点表及相关结果；BH校正覆盖原全部切片与四种目标组合。
visium_spot_table_total = pd.concat(visium_spot_rows, ignore_index=True)
visium_partial_table = pd.DataFrame(visium_partial_rows)
visium_partial_table["fdr"] = bh(visium_partial_table["p_value"])

visium_coverage_table = pd.DataFrame(
    [
        {
            "signature": "raw20",
            "n_defined": len(raw20),
            "n_detected_gbm030": visium_detected_counts.get("GBM030", np.nan),
            "n_detected_gbm049": visium_detected_counts.get("GBM049", np.nan),
        }
    ]
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "gse276841_raw20_spot_scores.csv.gz"
visium_spot_table_total.to_csv(output_path, index=False, compression={"method": "gzip", "mtime": 0})
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "gse276841_raw20_partial_rho.csv"
visium_partial_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "gse276841_raw20_slice_sensitivity.csv"
visium_partial_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "gse276841_raw20_coverage.csv"
visium_coverage_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

spatial = {"spots": visium_spot_table_total, "partial": visium_partial_table, "coverage": visium_coverage_table}

step38_gsea = (
    root
    / "write"
    / "38_independent_cohort_mg_inflammatory_recalculation"
    / "independent_fixed_program_targeted_gsea.csv"
)

# 读取该步骤原有输入表，保留列名、缺失值处理和读入顺序。
# 从Step38仅提取threshold=20的Miller原始top20程序及其NES、nominal P。
discovery_data_table = pd.read_csv(step38_gsea)
discovery_data_table = discovery_data_table.loc[
    discovery_data_table["threshold"].eq(20)
    & discovery_data_table["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")
].copy()
discovery_data_table["label"] = discovery_data_table["dataset"].map(
    {"GSE174554": "GSE174554 sc/snRNA", "GSE274546": "GSE274546 scRNA"}
)
discovery_data_table["nominal_p"] = discovery_data_table["pval"]
step38_table = discovery_data_table

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "step38_raw20_main_gsea.csv"
step38_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 流程图连接两个发现队列、GeoMx配对证据和两张Visium切片，保存对应源表。
flow_step38_table = step38_table
flow_geomx = geomx
flow_spatial = spatial

flow_geomx_main = flow_geomx["entry_summary"].loc[flow_geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]

flow_partial_main_table = flow_spatial["partial"].loc[flow_spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()

# 图中22对患者和2张切片沿用既有冻结文案，先与当前输入汇总核对。
if int(flow_geomx_main["n_pairs"]) != 22 or visium_spot_table_total["sample"].nunique() != 2:
    raise RuntimeError("Figure3 frozen 22-pair / 2-slice labels no longer match inputs")

flow_flow_source_table_total = pd.concat(
    [
        flow_step38_table[["dataset", "label", "n_pairs", "NES", "nominal_p"]].rename(
            columns = {"nominal_p": "step_p_or_fdr"}
        ),
        pd.DataFrame(
            [
                {
                    "dataset": "Artzi2025_GeoMx",
                    "label": "Artzi IBA1+ GeoMx",
                    "n_pairs": int(flow_geomx_main["n_pairs"]),
                    "NES": np.nan,
                    "step_p_or_fdr": flow_geomx_main["fdr"],
                },
                {
                    "dataset": "GSE276841",
                    "label": "GSE276841 Visium",
                    "n_pairs": 2,
                    "NES": np.nan,
                    "step_p_or_fdr": float(flow_partial_main_table["fdr"].min()),
                },
            ]
        ),
    ],
    ignore_index = True,
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "figure3_flowchart_summary.csv"
flow_flow_source_table_total.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

# 固定五个流程框位置；框中文字从发现队列和当前空间汇总表读取。
flow_fig, flow_ax = plt.subplots(figsize=(7.8, 3.6))
flow_ax.axis("off")
flow_boxes = {
    "discovery": (0.03, 0.22, 0.26, 0.56),
    "g174": (0.38, 0.56, 0.24, 0.24),
    "g274": (0.38, 0.20, 0.24, 0.24),
    "geomx": (0.71, 0.56, 0.24, 0.24),
    "visium": (0.71, 0.20, 0.24, 0.24),
}

# 各证据框和后续面板使用原NPG配色，键名仅表示颜色角色。
npg = {
    "red": "#E64B35",
    "blue": "#4DBBD5",
    "teal": "#00A087",
    "navy": "#3C5488",
    "salmon": "#F39B7F",
    "grayblue": "#8491B4",
    "mint": "#91D1C2",
    "brown": "#7E6148",
    "gray": "#6E7781",
}

def draw_box(key: str, lines: list[str], color: str) -> None:
    """输入原流程框位置名、文字及颜色；向当前Figure3流程图添加图形和文字；无返回值。"""

    flow_x, flow_y, flow_w, flow_h = flow_boxes[key]
    flow_rect = mpl.patches.FancyBboxPatch(
        (flow_x, flow_y),
        flow_w,
        flow_h,
        boxstyle  = "round,pad=0.012,rounding_size=0.02",
        linewidth = 1.0,
        edgecolor = color,
        facecolor = "white",
    )

    flow_ax.add_patch(flow_rect)
    flow_ax.text(
        flow_x + 0.02, flow_y + flow_h - 0.05, "\n".join(lines), ha="left", va="top", color="#222222", fontsize=8
    )

draw_box(
    "discovery", ["Independent raw20 anchor", "Miller raw top20 only", "No project-curated 16 replacement"], npg["navy"]
)
flow_g174 = flow_step38_table.loc[flow_step38_table["dataset"].eq("GSE174554")].iloc[0]
flow_g274 = flow_step38_table.loc[flow_step38_table["dataset"].eq("GSE274546")].iloc[0]
draw_box(
    "g174",
    [
        flow_g174["label"],
        f"{int(flow_g174['n_pairs'])} pairs",
        f"NES={flow_g174['NES']:.3f}",
        f"nominal P={flow_g174['nominal_p']:.2g}",
    ],
    npg["blue"],
)

draw_box(
    "g274",
    [
        flow_g274["label"],
        f"{int(flow_g274['n_pairs'])} pairs",
        f"NES={flow_g274['NES']:.3f}",
        f"nominal P={flow_g274['nominal_p']:.2g}",
    ],
    npg["teal"],
)

draw_box(
    "geomx",
    [
        "Artzi IBA1+ GeoMx",
        "22 IDH-wt pairs",
        f"mean delta={flow_geomx_main['mean_delta']:.3f}",
        f"FDR={flow_geomx_main['fdr']:.2g}",
    ],
    npg["red"],
)

draw_box(
    "visium",
    [
        "GSE276841 Visium",
        "2 untreated IDH-wt slices",
        f"MDSC min FDR={flow_partial_main_table.loc[flow_partial_main_table['target'].eq('MDSC'), 'fdr'].min():.2g}",
        f"MES min FDR={flow_partial_main_table.loc[flow_partial_main_table['target'].eq('MES'), 'fdr'].min():.2g}",
    ],
    npg["salmon"],
)

# 箭头维持发现队列通向GeoMx与Visium的原布局。
flow_arrows = [
    ((0.29, 0.50), (0.38, 0.68)),
    ((0.29, 0.50), (0.38, 0.32)),
    ((0.62, 0.68), (0.71, 0.68)),
    ((0.62, 0.32), (0.71, 0.32)),
]
for flow_start, flow_end in flow_arrows:
    flow_ax.annotate("", xy=flow_end, xytext=flow_start, arrowprops=dict(arrowstyle="->", lw=1.0, color="#555555"))
pdf_path = fig3_fig / "fig3a_flowchart.pdf"
png_path = fig3_fig / "fig3a_flowchart.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
flow_fig.canvas.draw()
flow_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(flow_fig)

# 主配对面板仅用strict_pass_idhwt患者，每条线连同一患者两个时间点。
paired_plot_geomx = geomx

paired_plot_data_table = (
    paired_plot_geomx["paired_scores"].loc[paired_plot_geomx["paired_scores"]["entry"].eq("strict_pass_idhwt")].copy()
)
paired_plot_summary_row = (
    paired_plot_geomx["entry_summary"].loc[paired_plot_geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "fig3b_geomx_paired_change.csv"
paired_plot_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

paired_plot_fig, paired_plot_ax = plt.subplots(figsize=(3.0, 3.8))
for _, paired_plot_row in paired_plot_data_table.iterrows():
    paired_plot_color = npg["red"] if paired_plot_row["trial_setting"] == "Nivolumab" else npg["blue"]
    paired_plot_ax.plot(
        [0, 1],
        [paired_plot_row["Primary"], paired_plot_row["Recurrence"]],
        color     = paired_plot_color,
        alpha     = 0.45,
        linewidth = 0.9,
    )

    paired_plot_ax.scatter(
        [0, 1],
        [paired_plot_row["Primary"], paired_plot_row["Recurrence"]],
        color      = paired_plot_color,
        s          = 14,
        zorder     = 3,
        linewidths = 0,
    )

# 黑线显示两时间点患者等权均值；标题读取入口均值差和FDR。
paired_plot_means = paired_plot_data_table[["Primary", "Recurrence"]].mean(axis=0)
paired_plot_ax.plot([0, 1], paired_plot_means, color="#111111", linewidth=2.0, marker="o", markersize=4.2, zorder=4)

paired_plot_ax.set_xticks([0, 1], ["Primary", "Recurrence"])
paired_plot_ax.set_ylabel("Raw20 score")
paired_plot_ax.text(
    0.5,
    1.01,
    f"22 pairs | mean delta={paired_plot_summary_row['mean_delta']:.3f} | FDR={paired_plot_summary_row['fdr']:.2g}",
    transform = paired_plot_ax.transAxes,
    ha        = "center",
    va        = "bottom",
    fontsize  = 7.2,
)

paired_plot_ax.spines[["top", "right"]].set_visible(False)
paired_plot_ax.grid(axis="y", color="#E8E8E8", linewidth=0.6)
pdf_path = fig3_fig / "fig3b_geomx_paired_change.pdf"
png_path = fig3_fig / "fig3b_geomx_paired_change.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
paired_plot_fig.canvas.draw()
paired_plot_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(paired_plot_fig)

# 基因森林图按平均log2差值降序、基因名升序排列，右侧显示原基因FDR。
genes_geomx = geomx

genes_data_table = genes_geomx["gene_summary"].loc[genes_geomx["gene_summary"]["entry"].eq("strict_pass_idhwt")].copy()
genes_data_table = genes_data_table.sort_values(["mean_log2_delta", "gene"], ascending=[False, True]).reset_index(
    drop = True
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "fig3c_geomx_gene_forest.csv"
genes_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

genes_fig, genes_ax = plt.subplots(figsize=(4.2, 4.8))
genes_y = np.arange(len(genes_data_table))
genes_colors = np.where(genes_data_table["mean_log2_delta"] >= 0, npg["red"], npg["navy"])
genes_ax.hlines(genes_y, genes_data_table["ci_low"], genes_data_table["ci_high"], color=genes_colors, linewidth=1.2)
genes_ax.scatter(genes_data_table["mean_log2_delta"], genes_y, color=genes_colors, s=24, zorder=3)

# 零线区分上升和下降；横轴边界根据置信区间及右侧文字留白计算。
genes_ax.axvline(0, color="#666666", linewidth=0.8)
genes_x_text = float(genes_data_table["ci_high"].max()) + 0.06
genes_x_min = float(genes_data_table["ci_low"].min()) - 0.08
genes_x_max = genes_x_text + 0.26
genes_ax.set_yticks(genes_y, genes_data_table["gene"])
genes_ax.invert_yaxis()
genes_ax.set_xlabel("Mean paired log2 delta")
genes_ax.set_xlim(genes_x_min, genes_x_max)
genes_ax.spines[["top", "right"]].set_visible(False)
genes_ax.grid(axis="x", color="#ECECEC", linewidth=0.6)

for genes_yy, genes_fdr in zip(genes_y, genes_data_table["fdr"]):
    genes_ax.text(genes_x_text, genes_yy, f"FDR={genes_fdr:.2g}", va="center", ha="left", fontsize=6.8, color="#333333")
pdf_path = fig3_fig / "fig3c_geomx_measurable_gene_forest.pdf"
png_path = fig3_fig / "fig3c_geomx_measurable_gene_forest.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
genes_fig.canvas.draw()
genes_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(genes_fig)

# 空间点图跨切片共享各指标2%至98%分位色限，逐切片绘制raw20、MDSC和MES。
maps_spatial = spatial

maps_spots_table = maps_spatial["spots"].copy()
maps_targets = [("raw20", "raw20"), ("mdsc_combined", "MDSC"), ("mes", "MES")]
maps_limits = {maps_key: maps_spots_table[maps_key].quantile([0.02, 0.98]).to_numpy() for maps_key, _ in maps_targets}

for maps_sample in sorted(maps_spots_table["sample"].unique()):
    maps_sample_table = maps_spots_table.loc[maps_spots_table["sample"].eq(maps_sample)].copy()

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = fig3_write / "source_csv" / f"{maps_sample}_spatial_maps.csv.gz"
    maps_sample_table.to_csv(output_path, index=False, compression={"method": "gzip", "mtime": 0})
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
    maps_fig, maps_axes = plt.subplots(1, 3, figsize=(6.8, 2.25))
    for maps_ax, (maps_key, maps_label) in zip(maps_axes, maps_targets):
        maps_low, maps_high = maps_limits[maps_key]
        maps_sca = maps_ax.scatter(
            maps_sample_table["pxl_col_in_fullres"],
            maps_sample_table["pixel_y_plot"],
            c          = maps_sample_table[maps_key],
            s          = 4,
            cmap       = "RdYlBu_r",
            vmin       = maps_low,
            vmax       = maps_high,
            linewidths = 0,
        )

        maps_ax.set_xticks([])
        maps_ax.set_yticks([])
        maps_ax.set_aspect("equal")
        maps_ax.set_title(maps_label, fontsize=8.5, pad=2)
        for maps_spine in maps_ax.spines.values():
            maps_spine.set_visible(False)
        maps_cbar = maps_fig.colorbar(maps_sca, ax=maps_ax, fraction=0.045, pad=0.02)
        maps_cbar.ax.tick_params(labelsize=6.5, length=2)
    maps_pdf_stem = (
        f"fig3d_{maps_sample}_spatial_maps" if maps_sample == "GBM030" else f"fig3e_{maps_sample}_spatial_maps"
    )

    pdf_path = fig3_fig / (f"{maps_pdf_stem}.pdf")
    maps_png_stem = (
        f"fig3d_{maps_sample}_spatial_maps" if maps_sample == "GBM030" else f"fig3e_{maps_sample}_spatial_maps"
    )
    png_path = fig3_fig / (f"{maps_png_stem}.png")

    # 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
    maps_fig.canvas.draw()
    maps_fig.savefig(pdf_path, bbox_inches="tight")

    # 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
    render_pdf_png(pdf_path, png_path, dpi=600)
    plt.close(maps_fig)

# 主相关图仅展示MDSC和MES，并以切片为纵轴、目标为纵向偏移。
rho_spatial = spatial

rho_data_table = rho_spatial["partial"].loc[rho_spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()
rho_data_table["sample_label"] = rho_data_table["sample"].map({"GBM030": "GBM030", "GBM049": "GBM049"})

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = fig3_write / "source_csv" / "fig3f_partial_rho.csv"
rho_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

rho_fig, rho_ax = plt.subplots(figsize=(4.2, 2.6))
rho_y_positions = {"GBM030": 1, "GBM049": 0}
rho_offset = {"MDSC": 0.10, "MES": -0.10}
rho_color = {"MDSC": npg["teal"], "MES": npg["brown"]}
for _, rho_row in rho_data_table.iterrows():
    rho_yy = rho_y_positions[rho_row["sample_label"]] + rho_offset[rho_row["target"]]
    rho_ax.scatter(
        rho_row["partial_rho"],
        rho_yy,
        s         = 44,
        color     = rho_color[rho_row["target"]],
        edgecolor = "white",
        linewidth = 0.6,
        zorder    = 3,
    )

    rho_ax.text(
        rho_row["partial_rho"] + 0.015, rho_yy, f"FDR={rho_row['fdr']:.2g}", va="center", fontsize=6.8, color="#333333"
    )

# 相关图保留零线、原横轴范围及目标颜色图例，FDR标注来自当前汇总表。
rho_ax.axvline(0, color="#666666", linewidth=0.8)
rho_ax.set_yticks([0, 1], ["GBM049", "GBM030"])
rho_ax.set_xlabel("Partial Spearman rho\n(adjusted for total myeloid score)")
rho_ax.set_xlim(0, max(0.62, rho_data_table["partial_rho"].max() + 0.12))
rho_ax.spines[["top", "right", "left"]].set_visible(False)
rho_ax.tick_params(axis="y", length=0)
rho_ax.grid(axis="x", color="#ECECEC", linewidth=0.6)

rho_handles = [
    mpl.lines.Line2D(
        [0],
        [0],
        marker          = "o",
        color           = "w",
        markerfacecolor = rho_color[rho_label],
        markeredgecolor = "white",
        markersize      = 7,
        label           = rho_label,
    )
    for rho_label in ["MDSC", "MES"]
]

rho_ax.legend(handles=rho_handles, frameon=False, loc="lower right")
pdf_path = fig3_fig / "fig3f_partial_rho.pdf"
png_path = fig3_fig / "fig3f_partial_rho.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
rho_fig.canvas.draw()
rho_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(rho_fig)

# 入口敏感性面板按均值差排序，对比严格筛选、放宽AOI及全IDH三种入口。
entry_geomx = geomx

entry_data_table = entry_geomx["entry_summary"].copy()
entry_data_table["label"] = entry_data_table["entry"].map(
    {
        "strict_pass_idhwt": "IBA1+ strict + IDH-wt",
        "all_aoi_idhwt": "IBA1+ all AOI + IDH-wt",
        "strict_pass_all_idh": "IBA1+ strict + all IDH",
    }
)

entry_data_table = entry_data_table.sort_values("mean_delta", ascending=False)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "s6a_entry_sensitivity.csv"
entry_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

entry_fig, entry_ax = plt.subplots(figsize=(4.0, 2.8))
entry_y = np.arange(len(entry_data_table))
entry_ax.scatter(
    entry_data_table["mean_delta"],
    entry_y,
    s         = 54,
    color     = [npg["red"], npg["salmon"], npg["navy"]][: len(entry_data_table)],
    edgecolor = "white",
    linewidth = 0.7,
    zorder    = 3,
)

entry_ax.axvline(0, color="#666666", linewidth=0.8)
entry_ax.set_yticks(entry_y, entry_data_table["label"])
entry_ax.invert_yaxis()
entry_ax.set_xlabel("Mean paired raw20 delta")
entry_ax.spines[["top", "right"]].set_visible(False)
entry_ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
for entry_x, entry_yy, entry_n_pairs, entry_fdr in zip(
    entry_data_table["mean_delta"], entry_y, entry_data_table["n_pairs"], entry_data_table["fdr"]
):
    entry_ax.text(
        entry_x + 0.012, entry_yy, f"n={entry_n_pairs}, FDR={entry_fdr:.2g}", va="center", fontsize=7, color="#333333"
    )

pdf_path = s6_fig / "s6a_entry_sensitivity.pdf"
png_path = s6_fig / "s6a_entry_sensitivity.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
entry_fig.canvas.draw()
entry_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(entry_fig)

# IDH敏感性保留主分析及全IDH入口中的野生型、突变型两个分层。
idh_geomx = geomx

idh_data_table = idh_geomx["subgroup_summary"].copy()
idh_data_table = idh_data_table.loc[
    (idh_data_table["entry"].eq("strict_pass_all_idh") & idh_data_table["subgroup"].isin(["IDH_WT", "IDH_mutated"]))
    | (idh_data_table["entry"].eq("strict_pass_idhwt") & idh_data_table["subgroup"].eq("All"))
].copy()
idh_labels = {
    ("strict_pass_idhwt", "All"): "Main analysis IDH-wt",
    ("strict_pass_all_idh", "IDH_WT"): "All-IDH entry: IDH-wt only",
    ("strict_pass_all_idh", "IDH_mutated"): "All-IDH entry: IDH-mutant only",
}

idh_display_labels = []
for row in idh_data_table.itertuples(index=False):
    idh_display_labels.append(idh_labels[(row.entry, row.subgroup)])
idh_data_table["label"] = idh_display_labels

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "s6b_idh_sensitivity.csv"
idh_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

idh_fig, idh_ax = plt.subplots(figsize=(4.2, 2.8))
idh_y = np.arange(len(idh_data_table))
idh_colors = [npg["red"], npg["blue"], npg["grayblue"]]
idh_ax.scatter(
    idh_data_table["mean_delta"],
    idh_y,
    s         = 54,
    color     = idh_colors[: len(idh_data_table)],
    edgecolor = "white",
    linewidth = 0.7,
    zorder    = 3,
)

idh_ax.axvline(0, color="#666666", linewidth=0.8)
idh_ax.set_yticks(idh_y, idh_data_table["label"])
idh_ax.invert_yaxis()
idh_ax.set_xlabel("Mean paired raw20 delta")
idh_ax.spines[["top", "right"]].set_visible(False)
idh_ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
for idh_x, idh_yy, idh_n_pairs, idh_fdr in zip(
    idh_data_table["mean_delta"], idh_y, idh_data_table["n_pairs"], idh_data_table["fdr"]
):
    idh_label = f"n={idh_n_pairs}"

    if pd.notna(idh_fdr):
        idh_label += f", FDR={idh_fdr:.2g}"
    idh_ax.text(idh_x + 0.012, idh_yy, idh_label, va="center", fontsize=7, color="#333333")
pdf_path = s6_fig / "s6b_idh_sensitivity.pdf"
png_path = s6_fig / "s6b_idh_sensitivity.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
idh_fig.canvas.draw()
idh_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(idh_fig)

# 治疗分层面板分别连接Control和Nivolumab患者的Primary与Recurrence分数。
treatment_geomx = geomx

treatment_data_table = (
    treatment_geomx["paired_scores"].loc[treatment_geomx["paired_scores"]["entry"].eq("strict_pass_idhwt")].copy()
)
treatment_summary_table = (
    treatment_geomx["subgroup_summary"]
    .loc[
        treatment_geomx["subgroup_summary"]["entry"].eq("strict_pass_idhwt")
        & treatment_geomx["subgroup_summary"]["subgroup"].isin(["Control", "Nivolumab"])
    ]
    .copy()
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "s6c_treatment_stratified_paired_change.csv"
treatment_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

treatment_fig, treatment_axes = plt.subplots(1, 2, figsize=(5.1, 3.2), sharey=True)
for treatment_ax, treatment_subgroup, treatment_color in zip(
    treatment_axes, ["Control", "Nivolumab"], [npg["blue"], npg["red"]]
):
    treatment_part_table = treatment_data_table.loc[treatment_data_table["trial_setting"].eq(treatment_subgroup)].copy()

    for _, treatment_row in treatment_part_table.iterrows():
        treatment_ax.plot(
            [0, 1],
            [treatment_row["Primary"], treatment_row["Recurrence"]],
            color     = treatment_color,
            alpha     = 0.45,
            linewidth = 0.9,
        )
        treatment_ax.scatter(
            [0, 1],
            [treatment_row["Primary"], treatment_row["Recurrence"]],
            color      = treatment_color,
            s          = 14,
            linewidths = 0,
            zorder     = 3,
        )

    # 每个治疗格叠加组内患者均值连线，并标注该治疗分层的配对人数及FDR。
    treatment_means = treatment_part_table[["Primary", "Recurrence"]].mean(axis=0)
    treatment_ax.plot([0, 1], treatment_means, color="#111111", linewidth=2.0, marker="o", markersize=4.0, zorder=4)
    treatment_row = treatment_summary_table.loc[treatment_summary_table["subgroup"].eq(treatment_subgroup)].iloc[0]
    treatment_ax.text(
        0.5,
        1.01,
        f"n={int(treatment_row['n_pairs'])} | FDR={treatment_row['fdr']:.2g}",
        transform = treatment_ax.transAxes,
        ha        = "center",
        va        = "bottom",
        fontsize  = 7,
    )

    treatment_ax.set_title(treatment_subgroup, fontsize=8.5, pad=3)
    treatment_ax.set_xticks([0, 1], ["Primary", "Recurrence"])
    treatment_ax.spines[["top", "right"]].set_visible(False)
    treatment_ax.grid(axis="y", color="#ECECEC", linewidth=0.6)
treatment_axes[0].set_ylabel("Raw20 score")
pdf_path = s6_fig / "s6c_treatment_stratified.pdf"
png_path = s6_fig / "s6c_treatment_stratified.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
treatment_fig.canvas.draw()
treatment_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(treatment_fig)

# 治疗交互面板展示各入口的Nivolumab减Control差值及原标签置换P值。
interaction_geomx = geomx

interaction_data_table = (
    interaction_geomx["subgroup_summary"]
    .loc[interaction_geomx["subgroup_summary"]["subgroup"].eq("Interaction_Nivolumab_minus_Control")]
    .copy()
)
interaction_data_table["label"] = interaction_data_table["entry"].map(
    {
        "strict_pass_idhwt": "IBA1+ strict + IDH-wt",
        "all_aoi_idhwt": "IBA1+ all AOI + IDH-wt",
        "strict_pass_all_idh": "IBA1+ strict + all IDH",
    }
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "s6d_interaction.csv"
interaction_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

interaction_fig, interaction_ax = plt.subplots(figsize=(4.0, 2.8))
interaction_y = np.arange(len(interaction_data_table))
interaction_ax.scatter(
    interaction_data_table["mean_delta"],
    interaction_y,
    s         = 54,
    color     = npg["brown"],
    edgecolor = "white",
    linewidth = 0.7,
    zorder    = 3,
)

interaction_ax.axvline(0, color="#666666", linewidth=0.8)
interaction_ax.set_yticks(interaction_y, interaction_data_table["label"])
interaction_ax.invert_yaxis()
interaction_ax.set_xlabel("Nivolumab minus control\nmean paired delta")
interaction_ax.spines[["top", "right"]].set_visible(False)
interaction_ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
for interaction_x, interaction_yy, interaction_p_value in zip(
    interaction_data_table["mean_delta"], interaction_y, interaction_data_table["sign_flip_p"]
):
    interaction_ax.text(
        interaction_x + 0.012, interaction_yy, f"P={interaction_p_value:.3g}", va="center", fontsize=7, color="#333333"
    )

pdf_path = s6_fig / "s6d_interaction.pdf"
png_path = s6_fig / "s6d_interaction.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
interaction_fig.canvas.draw()
interaction_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(interaction_fig)

# 逐切片敏感性热图按MDSC、E-MDSC、M-MDSC、MES固定顺序展示相关和FDR。
slice_spatial = spatial

slice_data_table = slice_spatial["partial"].copy()
slice_data_table = slice_data_table.loc[slice_data_table["target"].isin(["MDSC", "E-MDSC", "M-MDSC", "MES"])].copy()
slice_pivot_table = slice_data_table.pivot(index="target", columns="sample", values="partial_rho").reindex(
    ["MDSC", "E-MDSC", "M-MDSC", "MES"]
)

# 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
output_path = s6_write / "source_csv" / "s6e_slice_sensitivity.csv"
slice_data_table.to_csv(output_path, index=False)
if not output_path.is_file() or output_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {output_path.name}")

slice_fig, slice_ax = plt.subplots(figsize=(3.2, 2.9))
slice_data = slice_pivot_table.to_numpy(dtype=float)
slice_im = slice_ax.imshow(slice_data, cmap="RdYlBu_r", vmin=0, vmax=max(0.55, np.nanmax(slice_data)))
slice_ax.set_xticks(range(slice_pivot_table.shape[1]), slice_pivot_table.columns)
slice_ax.set_yticks(range(slice_pivot_table.shape[0]), slice_pivot_table.index)

slice_ax.tick_params(length=0)
for slice_i in range(slice_pivot_table.shape[0]):
    for slice_j in range(slice_pivot_table.shape[1]):
        slice_value = slice_data[slice_i, slice_j]
        slice_fdr = slice_data_table.loc[
            (slice_data_table["target"].eq(slice_pivot_table.index[slice_i]))
            & (slice_data_table["sample"].eq(slice_pivot_table.columns[slice_j])),
            "fdr",
        ].iloc[0]
        slice_ax.text(
            slice_j,
            slice_i,
            f"{slice_value:.2f}\nFDR={slice_fdr:.1g}",
            ha       = "center",
            va       = "center",
            fontsize = 6.5,
            color    = "white" if slice_value > 0.33 else "#222222",
        )

for slice_spine in slice_ax.spines.values():
    slice_spine.set_visible(False)
slice_cbar = slice_fig.colorbar(slice_im, ax=slice_ax, fraction=0.05, pad=0.04)
slice_cbar.set_label("Partial rho", fontsize=7)
slice_cbar.ax.tick_params(labelsize=6.5)
pdf_path = s6_fig / "s6e_slice_sensitivity.pdf"
png_path = s6_fig / "s6e_slice_sensitivity.png"

# 先实际绘制画布，再用原裁边、留白和元数据保存PDF母版。
slice_fig.canvas.draw()
slice_fig.savefig(pdf_path, bbox_inches="tight")

# 从刚保存的PDF渲染600 dpi PNG，文件与图像尺寸检查由共享工具完成。
render_pdf_png(pdf_path, png_path, dpi=600)
plt.close(slice_fig)

# 汇总文字继续读取发现队列NES/P值、GeoMx变化和两切片相关范围；正文文案不改。
report_step38_table = step38_table
report_geomx = geomx
report_spatial = spatial

report_geomx_main = (
    report_geomx["entry_summary"].loc[report_geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]
)

report_partial_main_table = (
    report_spatial["partial"].loc[report_spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()
)
report_gse174554_nes = report_step38_table.loc[report_step38_table["dataset"].eq("GSE174554"), "NES"].iloc[0]
report_gse174554_p = report_step38_table.loc[report_step38_table["dataset"].eq("GSE174554"), "nominal_p"].iloc[0]
report_gse274546_nes = report_step38_table.loc[report_step38_table["dataset"].eq("GSE274546"), "NES"].iloc[0]
report_gse274546_p = report_step38_table.loc[report_step38_table["dataset"].eq("GSE274546"), "nominal_p"].iloc[0]
report_mdsc_rho_min = report_partial_main_table.loc[report_partial_main_table["target"].eq("MDSC"), "partial_rho"].min()
report_mdsc_rho_max = report_partial_main_table.loc[report_partial_main_table["target"].eq("MDSC"), "partial_rho"].max()
report_mes_rho_min = report_partial_main_table.loc[report_partial_main_table["target"].eq("MES"), "partial_rho"].min()

report_mes_rho_max = report_partial_main_table.loc[report_partial_main_table["target"].eq("MES"), "partial_rho"].max()
report_lines = [
    "# Figure3 / Supplementary Figure6 raw20 summary",
    "",
    "## 主结果",
    "",
    (
        "- Step38 raw20 主程序在 GSE174554 为 NES="
        f"{report_gse174554_nes:.3f}"
        ", nominal P="
        f"{report_gse174554_p:.2g}"
        "; 在 GSE274546 为 NES="
        f"{report_gse274546_nes:.3f}"
        ", nominal P="
        f"{report_gse274546_p:.2g}"
        "."
    ),
    (
        "- Artzi IBA1+ GeoMx 严格 IDH-wt 22 对原发-复发患者平均 raw20 变化 "
        f"{report_geomx_main['mean_delta']:.3f}"
        ", FDR="
        f"{report_geomx_main['fdr']:.2g}"
        ", 正向 "
        f"{int(report_geomx_main['n_positive'])}"
        "/"
        f"{int(report_geomx_main['n_pairs'])}"
        "."
    ),
    (
        "- GSE276841 两张切片在扣除总髓系后，raw20 与 MDSC 的 partial rho 范围 "
        f"{report_mdsc_rho_min:.3f}"
        "-"
        f"{report_mdsc_rho_max:.3f}"
        "; 与 MES 的 partial rho 范围 "
        f"{report_mes_rho_min:.3f}"
        "-"
        f"{report_mes_rho_max:.3f}"
        "."
    ),
    "",
    "## 补图边界",
    "",
    "- GeoMx 入口、IDH 和治疗分层均未改变方向；治疗 interaction 仍不显著。",
    "- GSE276841 的逐切片敏感性继续支持 raw20 靠近 MDSC / MES 生态位。",
]

report_path = write_base / "39_raw20_figure_summary.md"
report_path.write_text("\n".join(report_lines), encoding="utf-8")
if not report_path.is_file() or report_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty output: {report_path.name}")
print("GeoMx pairs:", int(report_geomx_main["n_pairs"]), "Visium slices:", visium_spot_table_total["sample"].nunique())
