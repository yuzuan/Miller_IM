#!/usr/bin/env python3
"""目的：以GSE278456单细胞和GSE276841空间数据验证既定髓系程序的状态分布及空间相关。
输入：--root/write/30_dataset_rescue_search/source_metadata下的髓系/空间H5、位置文件和
    NIHMS2115262-supplement-S2_S5_S6.xls；同write模块的两个GEO sample_manifest.csv。
输出：同write模块GSE278456_author_marker_sets_used、single_cell_qc、program_gene_coverage、
    author_like_state_counts、patient_state_program_scores、patient_level_state_tests的CSV，
    GSE278456_cell_annotations_scores.csv.gz、GSE276841_spatial_partial_correlations.csv、
    GSE276841_spot_scores.csv.gz和GSE278456_GSE276841_validation_summary.md；
    同figures模块GSE278456_patient_level_Mg_program_by_state及GSE276841_spatial_program_maps的PDF/PNG。
运行：python scripts/06a_GSE278456髓系与空间验证.py --root DIR。
保留Python：scanpy/anndata/seaborn及既有统计实现构成公开复现路径；不改QC、marker、阈值或随机顺序。
"""

from __future__ import annotations

import argparse
import importlib
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # CLI、输出位置及全局种子沿用原入口，数据源为本地已下载文件。
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    write = args.root / "write/30_dataset_rescue_search"
    figures = args.root / "figures/30_dataset_rescue_search"
    source = write / "source_metadata"
    write.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 1
    np.random.seed(20260712)

    # 作者marker排除MG_FULL成员，保留FDR<0.05、Foldchange>0及原工作表前40顺序。
    mg_full = [
        "PDK4", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CH25H", "JUN", "SIGLEC8",
        "KLF6", "FOLR2", "AC253572.2", "NLRP3",
    ]
    marker_path = source / "NIHMS2115262-supplement-S2_S5_S6.xls"
    workbook = pd.ExcelFile(marker_path)
    excluded = set(mg_full)
    markers = {}
    for sheet in workbook.sheet_names:
        marker_source_table = pd.read_excel(marker_path, sheet_name=sheet)
        if "FDR" not in marker_source_table.columns:
            marker_source_table = pd.read_excel(marker_path, sheet_name=sheet, header=1)
        gene_column = marker_source_table.columns[0]
        chosen = marker_source_table.loc[
            marker_source_table["FDR"].lt(0.05) & marker_source_table["Foldchange"].gt(0), gene_column
        ].astype(str)
        markers[sheet] = [gene for gene in chosen if gene not in excluded][:40]

    # 实际使用的作者marker按原工作表顺序写出，每状态一行。
    marker_rows = [
        {"state": state, "n_genes": len(genes), "genes": ";".join(genes)} for state, genes in markers.items()
    ]
    author_marker_table = pd.DataFrame(marker_rows)
    author_marker_table.to_csv(write / "GSE278456_author_marker_sets_used.csv", index=False)

    # 每H5为一个样本；QC仍为200–8000基因、>=500 counts、线粒体比例<=20%。
    manifest_table = pd.read_csv(write / "GSE278456_sample_manifest.csv")
    manifest_table = manifest_table.set_index("geo_accession")
    objects = []
    qc_rows = []
    for path in sorted((source / "GSE278456_myeloid_h5").glob("*.h5")):
        accession = re.match(r"(GSM\d+)", path.name).group(1)
        sample = manifest_table.loc[accession, "sample_title"]
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        adata.obs_names = [f"{sample}:{barcode}" for barcode in adata.obs_names]
        adata.obs["sample"] = sample
        adata.obs["geo_accession"] = accession

        # 原10x基因名以MT-识别线粒体；QC在未标准化counts上计算。
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].ge(200)
            & adata.obs["n_genes_by_counts"].le(8_000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        qc_rows.append(
            {
                "sample": sample,
                "geo_accession": accession,
                "n_input": adata.n_obs,
                "n_qc_pass": int(keep.sum()),
                "fraction_qc_pass": float(keep.mean()),
            }
        )
        objects.append(adata[keep].copy())
    combined = ad.concat(objects, join="inner", merge="same")
    combined.uns["qc_summary"] = pd.DataFrame(qc_rows).to_dict("list")
    adata = combined

    # 程序评分及患者内检验按原调用顺序；符号翻转每次使用原种子及200000次抽样。
    qc_table = pd.DataFrame(adata.uns.pop("qc_summary"))
    qc_table.to_csv(write / "GSE278456_single_cell_qc.csv", index=False)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    # 每状态对照基因评分按原顺序计算，最大score决定author_like_state。
    for state, genes in markers.items():
        adata.obs[f"author_{state}"] = data_functions.control_gene_score(adata, genes, f"tmp_author_{state}")
    score_columns = [f"author_{state}" for state in markers]
    score_matrix = adata.obs[score_columns].to_numpy()
    state_names = np.array(list(markers))
    adata.obs["author_like_state"] = state_names[np.nanargmax(score_matrix, axis=1)]
    adata.obs["Mg_full"] = data_functions.control_gene_score(adata, mg_full, "tmp_Mg_full")
    mg_leading6 = ["BHLHE41", "CH25H", "FOLR2", "JUN", "SGK1", "SIGLEC8"]
    adata.obs["Mg_leading6"] = data_functions.control_gene_score(adata, mg_leading6, "tmp_Mg_leading6")
    myeloid_genes = ["PTPRC", "TYROBP", "AIF1", "LST1", "FCER1G", "CTSS", "CSF1R"]
    adata.obs["myeloid_abundance"] = data_functions.control_gene_score(adata, myeloid_genes, "tmp_myeloid")

    # 程序基因覆盖从当前AnnData变量名计算，不把定义基因数当成可测基因数。
    coverage = []
    for name, genes in {"Mg_full": mg_full, "Mg_leading6": mg_leading6, "myeloid_abundance": myeloid_genes}.items():
        use = [gene for gene in genes if gene in adata.var_names]
        coverage.append(
            {"program": name, "n_defined": len(genes), "n_detected": len(use), "detected_genes": ";".join(use)}
        )
    program_coverage_table = pd.DataFrame(coverage)
    program_coverage_table.to_csv(write / "GSE278456_program_gene_coverage.csv", index=False)

    # 每样本每状态计数，再连接该样本总数计算fraction，分母范围保持全部输入髓系细胞。
    state_count_table_total = (
        adata.obs.groupby(["sample", "author_like_state"], observed=True).size().rename("n_cells").reset_index()
    )
    sample_count_table = adata.obs.groupby("sample", observed=True).size().rename("n_total").reset_index()
    state_count_table_total = state_count_table_total.merge(sample_count_table, on="sample")
    state_count_table_total["fraction"] = state_count_table_total["n_cells"] / state_count_table_total["n_total"]
    state_count_table_total.to_csv(write / "GSE278456_author_like_state_counts.csv", index=False)

    # 患者状态表限原九个髓系状态，保留原均值与每组细胞数。
    myeloid_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]
    state_cell_table = adata.obs.loc[
        adata.obs["author_like_state"].isin(myeloid_states), ["sample", "author_like_state", "Mg_full", "Mg_leading6"]
    ]
    patient_state_table = (
        state_cell_table.groupby(["sample", "author_like_state"], observed=True)
        .agg(n_cells=("Mg_full", "size"), Mg_full=("Mg_full", "mean"), Mg_leading6=("Mg_leading6", "mean"))
        .reset_index()
    )
    patient_state_table.to_csv(write / "GSE278456_patient_state_program_scores.csv", index=False)

    test_rows = []

    # 输入样本细胞表、状态和程序列，返回原groupby顺序的患者内差值；不满足30/30细胞时跳过。
    def collect_state_deltas(cell_table: pd.DataFrame, state: str, program: str) -> list[tuple[str, float]]:
        """每个sample以目标状态均值减其余髓系状态均值，返回(sample,delta)列表。"""
        deltas = []
        for sample, sample_table in cell_table.groupby("sample", observed=True):
            target = sample_table.loc[sample_table["author_like_state"].eq(state), program]
            other = sample_table.loc[~sample_table["author_like_state"].eq(state), program]
            if len(target) >= 30 and len(other) >= 30:
                deltas.append((sample, target.mean() - other.mean()))

        return deltas

    for program in ["Mg_full", "Mg_leading6"]:
        for state in myeloid_states:
            deltas = collect_state_deltas(state_cell_table, state, program)
            values = np.array([value for _, value in deltas])
            test_rows.append(
                {
                    "program": program,
                    "state": state,
                    "n_patients": len(values),
                    "mean_within_patient_delta": float(values.mean()) if len(values) else np.nan,
                    "median_within_patient_delta": float(np.median(values)) if len(values) else np.nan,
                    "n_positive": int((values > 0).sum()),
                    "sign_flip_p": data_functions.sign_flip_p(values),
                    "wilcoxon_p": stats.wilcoxon(values).pvalue if len(values) >= 5 and np.any(values) else np.nan,
                }
            )
    test_table = pd.DataFrame(test_rows)
    test_table["fdr"] = test_table.groupby("program")["sign_flip_p"].transform(data_functions.bh)
    test_table.to_csv(write / "GSE278456_patient_level_state_tests.csv", index=False)

    # 图的状态顺序来自患者均值；box/strip绘图顺序及抖动随机状态不变。
    state_order = patient_state_table.groupby("author_like_state")["Mg_full"].mean().sort_values(ascending=False).index
    plt.figure(figsize=(9, 5.2))
    sns.boxplot(
        data = patient_state_table, x="author_like_state", y="Mg_full", order=state_order, color="#d9e4dc", fliersize=0
    )
    sns.stripplot(
        data  = patient_state_table,
        x     = "author_like_state",
        y     = "Mg_full",
        order = state_order,
        color = "#1f4e5f",
        size  = 3,
        alpha = 0.65,
    )
    plt.xlabel("")
    plt.ylabel("Mg-inflammatory program score")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.gcf().canvas.draw()
    pdf_path = figures / "GSE278456_patient_level_Mg_program_by_state.pdf"
    plt.savefig(pdf_path)
    render_pdf_png(pdf_path, figures / "GSE278456_patient_level_Mg_program_by_state.png", dpi=220)
    plt.close()

    # 输出每细胞的公开注释与score字段，压缩mtime继续固定为0。
    obs_columns = ["sample", "geo_accession", "author_like_state", "Mg_full", "Mg_leading6", "myeloid_abundance"]
    adata.obs[obs_columns].to_csv(
        write / "GSE278456_cell_annotations_scores.csv.gz",
        compression = {"method": "gzip", "mtime": 0},
    )

    # 每个空间样本独立标准化和评分，限in_tissue==1；偏Spearman以髓系信号的秩作协变量。
    spatial_manifest_table = pd.read_csv(write / "GSE276841_sample_manifest.csv")
    spatial_manifest_table = spatial_manifest_table.set_index("geo_accession")
    correlation_rows = []
    spot_tables = []
    map_frames = []
    for matrix_path in sorted((source / "GSE276841_spatial").glob("*_filtered_feature_bc_matrix.h5")):
        accession = re.match(r"(GSM\d+)", matrix_path.name).group(1)
        sample = spatial_manifest_table.loc[accession, "sample_title"]
        prefix = matrix_path.name.replace("_filtered_feature_bc_matrix.h5", "")
        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=10_000)
        sc.pp.log1p(adata)

        # 每spot一行score，所有对照基因评分沿用同一固定种子。
        spatial_score_table_total = pd.DataFrame(index=adata.obs_names)
        spatial_score_table_total["Mg_full"] = data_functions.control_gene_score(adata, mg_full, "tmp_Mg_full")
        spatial_score_table_total["Mg_leading6"] = data_functions.control_gene_score(
            adata, mg_leading6, "tmp_Mg_leading6"
        )
        spatial_score_table_total["myeloid_abundance"] = data_functions.control_gene_score(
            adata, myeloid_genes, "tmp_myeloid"
        )
        spatial_score_table_total["E_MDSC"] = data_functions.control_gene_score(adata, markers["E-MDSC"], "tmp_E_MDSC")
        spatial_score_table_total["M_MDSC"] = data_functions.control_gene_score(adata, markers["M-MDSC"], "tmp_M_MDSC")
        t4_mes_metabolic = [
            "CHI3L1", "CD44", "SERPINE1", "VEGFA", "ADM", "CA9", "BNIP3", "NDRG1", "SLC2A1", "HK2", "LDHA", "ENO1",
            "VIM", "LGALS3", "ANXA1", "TGFBI",
        ]
        spatial_score_table_total["T4_MES_metabolic"] = data_functions.control_gene_score(
            adata, t4_mes_metabolic, "tmp_T4_MES"
        )
        spatial_score_table_total["sample"] = sample
        spatial_score_table_total["geo_accession"] = accession

        positions_path = source / "GSE276841_spatial" / f"{prefix}_tissue_positions.csv.gz"
        position_table = pd.read_csv(positions_path)
        position_table = position_table.set_index("barcode")
        spatial_score_table_total = spatial_score_table_total.join(position_table, how="left")
        spatial_score_table_total = spatial_score_table_total.loc[spatial_score_table_total["in_tissue"].eq(1)].copy()
        spot_tables.append(spatial_score_table_total.reset_index(names="barcode"))


        # 每样本分别检验三个原target；统计量由共享函数按原秩残差公式计算。
        for target in ["E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            rho, p_value, n_spots = data_functions.partial_spearman(
                spatial_score_table_total["Mg_full"],
                spatial_score_table_total[target],
                spatial_score_table_total["myeloid_abundance"],
            )
            correlation_rows.append(
                {
                    "sample": sample,
                    "geo_accession": accession,
                    "target": target,
                    "n_spots": n_spots,
                    "partial_spearman_rho": rho,
                    "p_value": p_value,
                }
            )

        for program in ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            map_table = spatial_score_table_total[["pxl_col_in_fullres", "pxl_row_in_fullres", program]].copy()
            map_table["program"] = program
            map_table["value"] = map_table.pop(program)
            map_table["sample"] = sample
            map_frames.append(map_table)

    spatial_test_table = pd.DataFrame(correlation_rows)
    spatial_test_table["fdr"] = data_functions.bh(spatial_test_table["p_value"])
    spatial_test_table.to_csv(write / "GSE276841_spatial_partial_correlations.csv", index=False)
    spot_score_table_total = pd.concat(spot_tables, ignore_index=True)
    spot_score_table_total.to_csv(
        write / "GSE276841_spot_scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )

    # 拼接所有sample/program坐标；每小图按原2%–98%分位裁定色阶。
    map_table_total = pd.concat(map_frames, ignore_index=True)
    samples = map_table_total["sample"].drop_duplicates().tolist()
    programs = ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]
    fig, axes = plt.subplots(len(samples), len(programs), figsize=(13, 6.4), squeeze=False)
    for row, sample in enumerate(samples):
        for column_index, program in enumerate(programs):
            ax = axes[row, column_index]
            map_table = map_table_total.loc[
                map_table_total["sample"].eq(sample) & map_table_total["program"].eq(program)
            ]
            low, high = map_table["value"].quantile([0.02, 0.98])
            points = ax.scatter(
                map_table["pxl_col_in_fullres"],
                -map_table["pxl_row_in_fullres"],
                c    = map_table["value"],
                s    = 4,
                cmap = "magma",
                vmin = low,
                vmax = high,
                linewidths = 0,
            )
            ax.set_title(f"{sample} | {program}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02)
    plt.tight_layout()
    fig.canvas.draw()
    pdf_path = figures / "GSE276841_spatial_program_maps.pdf"
    plt.savefig(pdf_path)
    render_pdf_png(pdf_path, figures / "GSE276841_spatial_program_maps.png", dpi=220)
    plt.close()

    # 摘要直接读取本次两张统计表，不手填排名或数值；原结论文字保留。
    saved_test_table = pd.read_csv(write / "GSE278456_patient_level_state_tests.csv")
    saved_spatial_table = pd.read_csv(write / "GSE276841_spatial_partial_correlations.csv")
    top_test_table = saved_test_table.loc[saved_test_table["program"].eq("Mg_full")].sort_values("fdr").head(5)
    lines = [
        "# GSE278456 / GSE276841 independent validation",
        "",
        "## Single-cell patient-level state mapping",
        "",
        top_test_table.to_markdown(index=False),
        "",
        "## Spatial partial correlations adjusted for total myeloid signal",
        "",
        saved_spatial_table.to_markdown(index=False),
        "",
        (
            "Interpretation must remain claim-specific: this cohort is untreated primary IDH-wildtype GBM, "
            "so it tests phenotype and spatial ecology, not recurrence direction."
        ),
    ]
    (write / "GSE278456_GSE276841_validation_summary.md").write_text("\n".join(lines), encoding="utf-8")

    # 核对既有表、图和摘要的文件存在性及非空；不另建审核产物。
    for output_name in (
        "GSE278456_author_marker_sets_used.csv", "GSE278456_single_cell_qc.csv", "GSE278456_program_gene_coverage.csv",
        "GSE278456_author_like_state_counts.csv", "GSE278456_patient_state_program_scores.csv",
        "GSE278456_patient_level_state_tests.csv", "GSE278456_cell_annotations_scores.csv.gz",
        "GSE276841_spatial_partial_correlations.csv", "GSE276841_spot_scores.csv.gz",
        "GSE278456_GSE276841_validation_summary.md",
    ):
        output_path = write / output_name
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty validation output: {output_path}")

    print(f"Single-cell samples: {qc_table['sample'].nunique()}; patient-state tests: {len(test_table)}")
    print(f"Spatial samples: {spatial_manifest_table.shape[0]}; partial correlations: {len(spatial_test_table)}")
