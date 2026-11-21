#!/usr/bin/env python3
# 目的：从原始GSE278456 H5重建Miller raw20状态分数、患者效应与基因状态表。
# 输入：write/30_dataset_rescue_search内样本清单、QC表及外部数据根内原始H5和作者marker工作簿。
# 输出：write/48_figure2_identity_state_rebuild/Figure2/GSE278456_*.csv及细胞注释.csv.gz。
# 运行：在项目根使用既有sc环境执行 python scripts/06b_细胞身份与状态图.py。
# 保留Python：Scanpy、AnnData及Matplotlib是已公开复现路径；改成R超出本轮纯风格整理。
# 存量source_data/config、数据文件名、输出字段、write及figures路径保留原公开合同。

from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
sha256 = figure_style.sha256

def read_author_markers(
    marker_book: Path, raw20: list[str], myeloid_states: list[str]
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """输入作者marker工作簿；返回排除raw20重叠后的marker字典及逐状态审查表，缺状态或marker不足时报错。"""

    # 按工作簿sheet顺序读入状态marker；首行不含FDR时使用第二行表头。
    workbook = pd.ExcelFile(marker_book)
    markers: dict[str, list[str]] = {}
    rows: list[dict[str, object]] = []
    raw20_set = set(raw20)
    for sheet in workbook.sheet_names:
        data_table = pd.read_excel(marker_book, sheet_name=sheet)
        if "FDR" not in data_table.columns:
            data_table = pd.read_excel(marker_book, sheet_name=sheet, header=1)
        gene_col = data_table.columns[0]

        # 先筛FDR<0.05且Foldchange>0，再取前40个；之后才去除raw20重叠。
        original = (
            data_table.loc[
                pd.to_numeric(data_table["FDR"], errors="coerce").lt(0.05)
                & pd.to_numeric(data_table["Foldchange"], errors="coerce").gt(0),
                gene_col,
            ]
            .astype(str)
            .tolist()[:40]
        )

        # 从入选40个marker中删除raw20；剩余至少30个才允许用于状态赋值。
        cleaned = [gene for gene in original if gene not in raw20_set]
        removed = [gene for gene in original if gene in raw20_set]
        if len(cleaned) < 30:
            raise RuntimeError(f"{sheet} 删除raw20后只剩{len(cleaned)}个作者marker")
        markers[sheet] = cleaned
        rows.append(
            {
                "state": sheet,
                "n_original_top_markers": len(original),
                "n_after_raw20_exclusion": len(cleaned),
                "removed_raw20_genes": ";".join(removed),
                "markers_used": ";".join(cleaned),
            }
        )

    # 确认预设九个髓系状态均存在；审查表记录每套marker的删除项与最终集合。
    missing = sorted(set(myeloid_states) - set(markers))
    if missing:
        raise RuntimeError(f"作者marker工作簿缺少预设髓系状态: {missing}")
    marker_audit_table = pd.DataFrame(rows)
    result = (markers, marker_audit_table)
    return result

if __name__ == "__main__":

    # 沿用原Scanpy设置；本轮不更换库或评分实现。
    sc.settings.verbosity = 0

    # 输出沿用Step48目录；同时保留既有source_data和figures目录入口。
    root = Path(__file__).resolve().parents[1]

    write_root = root / "write/48_figure2_identity_state_rebuild/Figure2"

    write_root.mkdir(parents=True, exist_ok=True)

    source_out = write_root / "source_data"

    source_out.mkdir(parents=True, exist_ok=True)

    fig_root = root / "figures/48_figure2_identity_state_rebuild/Figure2"

    fig_root.mkdir(parents=True, exist_ok=True)

    # 原始H5和marker工作簿由MILLER_IM_DATA_ROOT定位；未设置时使用项目data目录。
    external_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()

    marker_book = external_root / "write/30_dataset_rescue_search/source_metadata/NIHMS2115262-supplement-S2_S5_S6.xls"

    h5_root = external_root / "write/30_dataset_rescue_search/source_metadata/GSE278456_myeloid_h5"

    if not h5_root.exists() or not marker_book.exists():
        raise FileNotFoundError("GSE278456原始H5或作者marker工作簿不可用")

    # Miller原始排序前20基因固定不变；基因大写是外部数据标识。
    raw20 = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3",
        "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]

    # 预设髓系状态及其原顺序保持不变；不由本轮结果重选状态。
    myeloid_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]

    markers, marker_audit_table = read_author_markers(marker_book, raw20, myeloid_states)
    atlas_markers = markers

    step30 = root / "write/30_dataset_rescue_search"

    # 读取原样本清单；每行是一份GEO样本，沿用IDH-wt入口。
    atlas_manifest_table = pd.read_csv(step30 / "GSE278456_sample_manifest.csv")
    atlas_manifest_table["geo_accession"] = atlas_manifest_table["geo_accession"].astype(str)

    # 读取Step30已冻结的输入细胞数和QC通过数；下面逐样本核对。
    atlas_existing_qc_table = pd.read_csv(step30 / "GSE278456_single_cell_qc.csv").set_index("geo_accession")
    atlas_paths = sorted(h5_root.glob("*.h5"))
    if len(atlas_paths) != 21:
        raise RuntimeError(f"GSE278456锁定H5数量应为21，实际为{len(atlas_paths)}")

    # 按样本分别收集细胞标签、状态均值、患者效应、基因表达、QC及输入锁定记录。
    atlas_all_marker_genes = set().union(*atlas_markers.values())
    atlas_cell_frames: list[pd.DataFrame] = []
    atlas_state_score_frames: list[pd.DataFrame] = []
    atlas_effect_rows: list[dict[str, object]] = []
    atlas_gene_patient_rows: list[dict[str, object]] = []

    atlas_coverage_rows: list[dict[str, object]] = []
    atlas_locked_rows: list[dict[str, object]] = []
    atlas_qc_rows: list[dict[str, object]] = []

    # 按H5文件名排序处理；每份输入必须唯一匹配IDH-wt样本清单。
    for atlas_path in atlas_paths:
        atlas_accession = atlas_path.name.split("_")[0]
        atlas_manifest_match_table = atlas_manifest_table.loc[atlas_manifest_table["geo_accession"].eq(atlas_accession)]
        if len(atlas_manifest_match_table) != 1:
            raise RuntimeError(f"{atlas_accession}无法唯一匹配样本manifest")
        atlas_manifest_row = atlas_manifest_match_table.iloc[0]
        if not bool(atlas_manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"锁定H5中出现非IDH-wt GBM: {atlas_accession}")
        atlas_sample = str(atlas_manifest_row["sample_title"])

        # 读取10x计数，唯一化基因名，并计算线粒体比例及细胞计数QC指标。
        atlas_adata = sc.read_10x_h5(atlas_path)

        atlas_adata.var_names_make_unique()
        atlas_n_input = atlas_adata.n_obs
        atlas_n_genes = atlas_adata.n_vars
        atlas_adata.var["mt"] = atlas_adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(atlas_adata, qc_vars=["mt"], inplace=True, percent_top=None)

        # 保留原QC阈值：200至8000个基因、至少500计数、线粒体比例不超过20%。
        atlas_keep = (
            atlas_adata.obs["n_genes_by_counts"].between(200, 8000)
            & atlas_adata.obs["total_counts"].ge(500)
            & atlas_adata.obs["pct_counts_mt"].le(20)
        )
        atlas_adata = atlas_adata[atlas_keep].copy()
        atlas_n_qc = atlas_adata.n_obs
        if atlas_accession not in atlas_existing_qc_table.index:
            raise RuntimeError(f"旧QC表缺少{atlas_accession}")
        atlas_old = atlas_existing_qc_table.loc[atlas_accession]

        if int(atlas_old["n_input"]) != atlas_n_input or int(atlas_old["n_qc_pass"]) != atlas_n_qc:
            raise RuntimeError(f"{atlas_accession} QC细胞数未复现")

        # 逐样本核对19个可测raw20基因；既有唯一缺失项为AC253572.2。
        atlas_present_raw20 = [atlas_gene for atlas_gene in raw20 if atlas_gene in atlas_adata.var_names]

        atlas_missing_raw20 = [atlas_gene for atlas_gene in raw20 if atlas_gene not in atlas_adata.var_names]
        if len(atlas_present_raw20) != 19 or atlas_missing_raw20 != ["AC253572.2"]:
            raise RuntimeError(
                f"{atlas_accession} raw20覆盖异常: {len(atlas_present_raw20)}/20, missing={atlas_missing_raw20}"
            )
        atlas_coverage_rows.append(
            {
                "geo_accession": atlas_accession,
                "sample": atlas_sample,
                "n_defined": len(raw20),
                "n_present": len(atlas_present_raw20),
                "present_genes": ";".join(atlas_present_raw20),
                "missing_genes": ";".join(atlas_missing_raw20),
            }
        )

        atlas_qc_rows.append(
            {
                "sample": atlas_sample,
                "geo_accession": atlas_accession,
                "n_input": atlas_n_input,
                "n_qc_pass": atlas_n_qc,
                "fraction_qc_pass": atlas_n_qc / atlas_n_input,
            }
        )
        atlas_locked_rows.append(
            {
                "geo_accession": atlas_accession,
                "sample": atlas_sample,
                "path": str(atlas_path),
                "sha256": sha256(atlas_path),
                "n_input_cells": atlas_n_input,
                "n_qc_pass_cells": atlas_n_qc,
                "n_genes": atlas_n_genes,
            }
        )

        # 保留QC后原始计数供pseudobulk使用；以GSM前缀避免跨样本barcode重名。
        atlas_counts = atlas_adata.X.copy()
        atlas_raw20_indices = atlas_adata.var_names.get_indexer(atlas_present_raw20)
        atlas_adata.obs_names = [
            f"{atlas_accession}:{atlas_barcode}" for atlas_barcode in atlas_adata.obs_names.astype(str)
        ]

        # 按原10,000计数尺度归一化，随后log1p；不改变原始计数副本。
        sc.pp.normalize_total(atlas_adata, target_sum=10000)
        sc.pp.log1p(atlas_adata)

        # 对照池同时排除raw20及全部作者marker，保持至少5000个中性候选基因。
        atlas_excluded_controls = set(raw20) | atlas_all_marker_genes
        atlas_neutral_control_pool = [
            atlas_gene for atlas_gene in atlas_adata.var_names if atlas_gene not in atlas_excluded_controls
        ]
        if len(atlas_neutral_control_pool) < 5000:
            raise RuntimeError(f"{atlas_accession} state-neutral control pool过小")

        # 逐套计算所有作者marker分数；竞争集合不局限于九个髓系状态。
        atlas_score_columns: list[str] = []
        for atlas_state, atlas_genes in atlas_markers.items():
            atlas_use = [atlas_gene for atlas_gene in atlas_genes if atlas_gene in atlas_adata.var_names]
            if len(atlas_use) < 30:
                raise RuntimeError(f"{atlas_accession} {atlas_state}可测marker少于30")
            atlas_column = f"state_score_{atlas_state}"

            # 保留ctrl_size=50、n_bins=25和原种子；每套评分重置种子是原方法，移出会改变结果。
            sc.tl.score_genes(
                atlas_adata,
                gene_list    = atlas_use,
                gene_pool    = atlas_neutral_control_pool + atlas_use,
                ctrl_size    = 50,
                n_bins       = 25,
                score_name   = atlas_column,
                random_state = 20260713,
                ctrl_as_ref  = False,
                use_raw      = False,
            )

            atlas_score_columns.append(atlas_column)

        # 每行是一个细胞、每列是一套作者marker分数；按原argsort取第一和第二名。
        atlas_score_matrix = atlas_adata.obs[atlas_score_columns].to_numpy(float)
        atlas_order = np.argsort(atlas_score_matrix, axis=1)
        atlas_state_names = np.array(list(atlas_markers), dtype=object)
        atlas_winner_index = atlas_order[:, -1]
        atlas_runner_index = atlas_order[:, -2]
        atlas_winners = atlas_state_names[atlas_winner_index]
        atlas_winning_score = atlas_score_matrix[np.arange(atlas_adata.n_obs), atlas_winner_index]
        atlas_runner_score = atlas_score_matrix[np.arange(atlas_adata.n_obs), atlas_runner_index]

        # 保留ctrl_size=50、n_bins=25和原种子；每套评分重置种子是原方法，移出会改变结果。
        sc.tl.score_genes(
            atlas_adata,
            gene_list    = atlas_present_raw20,
            gene_pool    = atlas_neutral_control_pool + atlas_present_raw20,
            ctrl_size    = 50,
            n_bins       = 25,
            score_name   = "raw20_score",
            random_state = 20260713,
            ctrl_as_ref  = False,
            use_raw      = False,
        )

        # 保存全部QC细胞的获胜状态、runner-up差值和raw20分数，并标记预设髓系细胞。
        atlas_raw20_scores = atlas_adata.obs["raw20_score"].to_numpy(float)
        atlas_focus = np.isin(atlas_winners, myeloid_states)

        atlas_cell_frames.append(
            pd.DataFrame(
                {
                    "cell_id": atlas_adata.obs_names.astype(str),
                    "geo_accession": atlas_accession,
                    "sample": atlas_sample,
                    "author_marker_state": atlas_winners,
                    "state_score": atlas_winning_score,
                    "runner_up_state": atlas_state_names[atlas_runner_index],
                    "state_margin": atlas_winning_score - atlas_runner_score,
                    "raw20_score": atlas_raw20_scores,
                    "in_prespecified_myeloid_states": atlas_focus,
                }
            )
        )

        # 提取raw20标准化表达；只在原稀疏矩阵分支转稠密，不改变精度。
        atlas_normalized_raw20 = atlas_adata[:, atlas_present_raw20].X
        if sparse.issparse(atlas_normalized_raw20):
            atlas_normalized_raw20 = atlas_normalized_raw20.toarray()

        atlas_normalized_raw20 = np.asarray(atlas_normalized_raw20, dtype=np.float32)

        # 每状态至少30个细胞才汇总表达；其余髓系细胞不足30时仍保留表达汇总。
        for atlas_state in myeloid_states:
            atlas_target = atlas_focus & (atlas_winners == atlas_state)
            atlas_rest = atlas_focus & (atlas_winners != atlas_state)
            atlas_n_target = int(atlas_target.sum())
            atlas_n_rest = int(atlas_rest.sum())
            if atlas_n_target < 30:
                continue

            # 先记录患者内该状态的平均raw20分数，再逐基因记录均值和检出比例。
            atlas_state_score_frames.append(
                pd.DataFrame(
                    {
                        "geo_accession": [atlas_accession],
                        "sample": [atlas_sample],
                        "state": [atlas_state],
                        "n_cells": [atlas_n_target],
                        "mean_raw20_score": [float(atlas_raw20_scores[atlas_target].mean())],
                    }
                )
            )

            atlas_state_expr = atlas_normalized_raw20[atlas_target]
            for atlas_gene_index, atlas_gene in enumerate(atlas_present_raw20):
                atlas_gene_patient_rows.append(
                    {
                        "geo_accession": atlas_accession,
                        "sample": atlas_sample,
                        "state": atlas_state,
                        "gene": atlas_gene,
                        "n_cells": atlas_n_target,
                        "mean_log_expression": float(atlas_state_expr[:, atlas_gene_index].mean()),
                        "pct_detected": float((atlas_state_expr[:, atlas_gene_index] > 0).mean() * 100),
                    }
                )

            # 两侧均至少30个细胞时计算状态减其余髓系差值，以及原始计数pseudobulk差值。
            if atlas_n_rest < 30:
                continue

            atlas_delta = float(atlas_raw20_scores[atlas_target].mean() - atlas_raw20_scores[atlas_rest].mean())

            # 目标状态先汇总全部基因的文库总量，再对raw20计算log2(CPM+1)均值。
            atlas_target_counts = atlas_counts[atlas_target]
            atlas_target_library = float(atlas_target_counts.sum())
            atlas_target_pb = np.nan
            if atlas_target_library > 0:
                atlas_target_gene_counts = np.asarray(atlas_target_counts[:, atlas_raw20_indices].sum(axis=0)).ravel()
                atlas_target_log_cpm = np.log2(atlas_target_gene_counts / atlas_target_library * 1_000_000 + 1)
                atlas_target_pb = float(atlas_target_log_cpm.mean())

            # 其余髓系使用相同公式；零文库仍为NaN，不纳入替代分数。
            atlas_rest_counts = atlas_counts[atlas_rest]
            atlas_rest_library = float(atlas_rest_counts.sum())
            atlas_rest_pb = np.nan
            if atlas_rest_library > 0:
                atlas_rest_gene_counts = np.asarray(atlas_rest_counts[:, atlas_raw20_indices].sum(axis=0)).ravel()
                atlas_rest_log_cpm = np.log2(atlas_rest_gene_counts / atlas_rest_library * 1_000_000 + 1)
                atlas_rest_pb = float(atlas_rest_log_cpm.mean())
            atlas_effect_rows.append(
                {
                    "geo_accession": atlas_accession,
                    "sample": atlas_sample,
                    "state": atlas_state,
                    "n_target": atlas_n_target,
                    "n_rest": atlas_n_rest,
                    "raw20_delta": atlas_delta,
                    "raw20_pseudobulk_delta": atlas_target_pb - atlas_rest_pb,
                }
            )

        del atlas_adata, atlas_counts, atlas_normalized_raw20

    # 把逐样本记录拼成冻结输入、QC和覆盖表；一行仍对应原样本。
    atlas_locked_table = pd.DataFrame(atlas_locked_rows).sort_values("geo_accession")
    atlas_qc_table = pd.DataFrame(atlas_qc_rows).sort_values("geo_accession")
    atlas_coverage_table = pd.DataFrame(atlas_coverage_rows).sort_values("geo_accession")
    atlas_cells_table_total = pd.concat(atlas_cell_frames, ignore_index=True)
    atlas_state_scores_table_total = pd.concat(atlas_state_score_frames, ignore_index=True)
    atlas_effects_table = pd.DataFrame(atlas_effect_rows)
    atlas_gene_patient_table = pd.DataFrame(atlas_gene_patient_rows)

    # 21份样本、131780个输入细胞及120766个QC细胞均为既有冻结对账值。
    if atlas_locked_table["geo_accession"].nunique() != 21 or atlas_locked_table["sample"].nunique() != 21:
        raise RuntimeError("GSE278456锁定manifest不是21个唯一GSM/患者")

    if int(atlas_qc_table["n_input"].sum()) != 131780 or int(atlas_qc_table["n_qc_pass"].sum()) != 120766:
        raise RuntimeError("GSE278456总细胞数未复现Step30")
    if atlas_cells_table_total["cell_id"].duplicated().any():
        raise RuntimeError("GSE278456 cell ID不唯一")

    # 每状态用患者差值计算原bootstrap区间和精确符号翻转P值，患者等权。
    atlas_summary_rows = []
    for atlas_index, atlas_state in enumerate(myeloid_states):
        atlas_part_table = atlas_effects_table.loc[atlas_effects_table["state"].eq(atlas_state)].copy()
        atlas_values = atlas_part_table["raw20_delta"].to_numpy(float)

        def bootstrap_mean_ci(values: np.ndarray, seed: int) -> tuple[float, float, float]:
            """输入患者差值及原种子；返回均值和20,000次有放回抽样95%区间，空输入返回NaN。"""

            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]

            if len(values) == 0:
                result = (np.nan, np.nan, np.nan)
                return result
            mean = float(values.mean())
            if len(values) == 1:
                result = (mean, mean, mean)
                return result

            # 每状态使用原种子独立抽样20000次；区间取均值分布的2.5%和97.5%分位。
            rng = np.random.default_rng(seed)
            draws = rng.choice(values, size=(20000, len(values)), replace=True).mean(axis=1)
            lo, hi = np.quantile(draws, [0.025, 0.975])

            result = (mean, float(lo), float(hi))
            return result

        atlas_mean, atlas_low, atlas_high = bootstrap_mean_ci(atlas_values, 20260713 + atlas_index)
        atlas_pb_values = atlas_part_table["raw20_pseudobulk_delta"].to_numpy(float)

        def exact_two_sided_signflip(values: np.ndarray) -> float:
            """输入患者差值向量；返回精确双侧符号翻转P值。保留首项固定及1e-12容差，少于3项返回NaN。"""

            values = np.asarray(values, dtype=float)

            values = values[np.isfinite(values)]
            n = len(values)
            if n < 3:
                result = np.nan
                return result
            observed = abs(values.mean())

            # 双侧绝对值分布具有整体符号对称性，固定首项为正可精确减半枚举空间。
            tail = values[1:]
            total = 1 << (n - 1)
            extreme = 0

            # 以65536个符号组合为块枚举；分块只限制内存，不抽样或重排患者差值。
            chunk_size = 65536
            bit_positions = np.arange(n - 1, dtype=np.uint64)
            for start in range(0, total, chunk_size):
                stop = min(start + chunk_size, total)
                codes = np.arange(start, stop, dtype=np.uint64)[:, None]
                signs = (((codes >> bit_positions) & 1).astype(np.float64) * 2.0) - 1.0
                means = (values[0] + signs @ tail) / n
                extreme += int(np.count_nonzero(np.abs(means) >= observed - 1e-12))
            result = extreme / total
            return result

        # 汇总患者均值、正向人数、精确P值及pseudobulk方向一致性，保留原输出字段。
        atlas_summary_rows.append(
            {
                "state": atlas_state,
                "n_patients": len(atlas_values),
                "mean_delta": atlas_mean,
                "ci95_low": atlas_low,
                "ci95_high": atlas_high,
                "median_delta": float(np.median(atlas_values)) if len(atlas_values) else np.nan,
                "n_positive": int((atlas_values > 0).sum()),
                "exact_two_sided_signflip_p": exact_two_sided_signflip(atlas_values),
                "mean_pseudobulk_delta": float(atlas_pb_values.mean()) if len(atlas_pb_values) else np.nan,
                "pseudobulk_n_positive": int((atlas_pb_values > 0).sum()),
                "score_pseudobulk_mean_direction_concordant": bool(
                    np.sign(atlas_mean) == np.sign(atlas_pb_values.mean())
                )
                if len(atlas_values)
                else False,
            }
        )

    # 每行是一种预设状态；BH校正范围仍为这九个状态。
    atlas_effect_summary_table = pd.DataFrame(atlas_summary_rows)
    atlas_ok = atlas_effect_summary_table["exact_two_sided_signflip_p"].notna()
    atlas_effect_summary_table["fdr_9_states"] = np.nan
    atlas_effect_summary_table.loc[atlas_ok, "fdr_9_states"] = multipletests(
        atlas_effect_summary_table.loc[atlas_ok, "exact_two_sided_signflip_p"], method="fdr_bh"
    )[1]

    # 先在患者内汇总，再在同一基因的状态间标准化，保留患者等权。
    atlas_gene_summary_table = (
        atlas_gene_patient_table.groupby(["state", "gene"], observed=True)
        .agg(
            n_patients                           = ("sample", "nunique"),
            patient_balanced_mean_log_expression = ("mean_log_expression", "mean"),
            patient_balanced_pct_detected        = ("pct_detected", "mean"),
        )
        .reset_index()
    )
    atlas_gene_summary_table["within_gene_z"] = np.nan

    # 在同一基因的九状态均值间计算总体标准差z分数；零方差按原规则记0。
    for _, group_table in atlas_gene_summary_table.groupby("gene", observed=True):
        group_values = group_table["patient_balanced_mean_log_expression"]
        group_sd = group_values.std(ddof=0)
        if group_sd > 0:
            group_z = (group_values - group_values.mean()) / group_sd
        else:
            group_z = 0
        atlas_gene_summary_table.loc[group_table.index, "within_gene_z"] = group_z

    # 将九张结果表映射到既有文件后缀；cells保留gzip格式，其余写CSV。
    results = {
        "locked_manifest": atlas_locked_table,
        "qc": atlas_qc_table,
        "coverage": atlas_coverage_table,
        "cells": atlas_cells_table_total,
        "patient_state_scores": atlas_state_scores_table_total,
        "patient_state_effects": atlas_effects_table,
        "state_effect_summary": atlas_effect_summary_table,
        "gene_patient": atlas_gene_patient_table,
        "gene_summary": atlas_gene_summary_table,
    }

    # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
    output_path = write_root / "GSE278456_author_marker_audit.csv"
    marker_audit_table.to_csv(output_path, index=False)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output_path.name}")
    for name, data_table in results.items():
        if name == "cells":

            # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
            output_path = write_root / "GSE278456_cell_annotations_true_raw20.csv.gz"
            data_table.to_csv(output_path, index=False, compression={"method": "gzip", "mtime": 0})
            if not output_path.is_file() or output_path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty output: {output_path.name}")
        else:

            # 写出原合同表格并核对本次文件存在且非空；不以历史文件替代当前输出。
            output_path = write_root / f"GSE278456_{name}.csv"
            data_table.to_csv(output_path, index=False)
            if not output_path.is_file() or output_path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty output: {output_path.name}")
    print("State summaries:", len(atlas_effect_summary_table))
