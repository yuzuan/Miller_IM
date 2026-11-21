#!/usr/bin/env python3
"""目的：从GSE274546原始snRNA矩阵及singlet分类，重建全细胞图谱并独立筛选clean myeloid。
输入：MILLER_IM_DATA_ROOT/GSE274546标准矩阵及results/manifests/gse274546_standard_libraries.csv；
    write/36_gse274546_raw_independent_reannotation/01_qc_doublets/scDblFinder_cell_classifications.csv.gz。
输出：同write模块02_full_cell的training_sketch_*、全细胞注释/metadata/marker/pseudobulk和embedding；
    03_myeloid的clean h5ad、metadata、cluster复核、Harmony比较和三组patient_condition*_pseudobulk；
    paired_patient_threshold_audit.csv、模块根目录raw_independent_reannotation_audit.json；
    figures/36_gse274546_raw_independent_reannotation/01_fullcell_umap_lineage、02_myeloid_umap_cluster、
    03_myeloid_umap_identity的PDF及由对应PDF转换的PNG。既有字段和文件名保持不变。
运行：python scripts/03c_GSE274546独立重注释.py。
保留Python：scanpy/anndata/harmonypy及sklearn SVD是原复现路径，本轮不改抽样、阈值或分析口径。
"""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path

# 线程设置必须先于数值库导入，保持原运行约束及随机计算路径。
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import anndata as ad
import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # 多段调用共用原种子；不移入或移出原抽样循环，不调整抽样次数和顺序。
    seed = 274546
    np.random.seed(seed)
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    root = Path(__file__).resolve().parents[1]
    source_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    manifest_path = source_root / "results" / "manifests" / "gse274546_standard_libraries.csv"

    manifest_table = pd.read_csv(manifest_path)
    manifest_table = manifest_table.loc[manifest_table["progression"].isin(["Primary", "1st Recurrent"])].copy()
    manifest_table["pair_id"] = manifest_table["patient_id"]
    manifest_table["condition"] = np.where(manifest_table["progression"].eq("Primary"), "Primary", "Recurrent")
    manifest_table = manifest_table.sort_values(["patient_id", "progression_order", "library_id"]).reset_index(
        drop = True
    )

    # 111文库/59患者是公开manifest冻结检查，本轮保留，不能放宽。
    expected_library_count = 111
    expected_patient_count = 59
    if (
        len(manifest_table) != expected_library_count
        or manifest_table["patient_id"].nunique() != expected_patient_count
    ):
        raise ValueError(f"Expected {expected_library_count} libraries from {expected_patient_count} patients")

    # singlet分类每行一个capture/barcode；先保留singlet，再检查重复标识。
    base = root / "write" / "36_gse274546_raw_independent_reannotation"
    classification_path = base / "01_qc_doublets" / "scDblFinder_cell_classifications.csv.gz"
    classification_table = pd.read_csv(classification_path, dtype=str)
    classification_table = classification_table.loc[classification_table["scDblFinder_class"].eq("singlet")].copy()
    if classification_table.duplicated(["capture_id", "barcode"]).any():
        raise ValueError("Duplicated singlet classifications")

    # 每文库最多500个singlet，保留rng.choice及排序顺序；在原计数上选2000 HVG。
    raw_root = source_root / "GSE274546"
    objects = []
    reference_genes: list[str] | None = None
    rng = np.random.default_rng(seed)
    sketch_rows = []
    for row in manifest_table.itertuples(index=False):
        matrix, genes, barcodes = data_functions.read_gse274546_library(row.library_id, raw_root)
        matrix, cell_metadata_table = data_functions.select_singlets(
            row.library_id, matrix, barcodes, classification_table
        )
        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")

        # 每文库独立抽取至多500个细胞，先抽样后排序，与原随机序列一致。
        n_keep = min(500, matrix.shape[0])
        selected = np.sort(rng.choice(matrix.shape[0], size=n_keep, replace=False))
        observation_table = cell_metadata_table.iloc[selected].copy()
        observation_table.index = [f"{row.library_id}:{barcode}" for barcode in observation_table["barcode"]]
        obj = ad.AnnData(X=matrix[selected], obs=observation_table, var=pd.DataFrame(index=reference_genes))
        objects.append(obj)
        sketch_rows.append({"library_id": row.library_id, "singlets": matrix.shape[0], "sketch_cells": n_keep})
    if reference_genes is None:
        raise RuntimeError("No libraries were loaded")

    # 拼接sketch后仅在原计数选择HVG；SVD输入保持1e4标准化后的log表达。
    sketch = ad.concat(objects, join="inner", merge="same", index_unique=None)
    sketch.X = sketch.X.tocsr().astype(np.int32)
    sketch.obs_names_make_unique()
    full_dir = base / "02_full_cell"
    full_dir.mkdir(parents=True, exist_ok=True)
    sketch_table = pd.DataFrame(sketch_rows)
    sketch_table.to_csv(full_dir / "training_sketch_by_library.csv", index=False)

    # 固定2000 HVG用于训练和髓系重聚类；保留原seurat_v3与batch_key=None。
    n_hvg = 2000
    sc.pp.highly_variable_genes(
        sketch,
        n_top_genes  = n_hvg,
        flavor       = "seurat_v3",
        batch_key    = None,
        check_values = True,
    )
    hvg_genes = sketch.var_names[sketch.var["highly_variable"]].astype(str).tolist()
    if len(hvg_genes) != n_hvg:
        raise RuntimeError(f"Expected {n_hvg} HVGs, found {len(hvg_genes)}")
    hvg_table = pd.DataFrame({"gene": hvg_genes})
    hvg_table.to_csv(full_dir / "training_sketch_hvg2000.csv", index=False)
    sc.pp.normalize_total(sketch, target_sum=1e4)
    sc.pp.log1p(sketch)
    hvg_matrix = sketch[:, hvg_genes].X.tocsr()
    full_n_pcs = 40
    svd = TruncatedSVD(n_components=full_n_pcs, random_state=seed)
    sketch.obsm["X_pca"] = svd.fit_transform(hvg_matrix).astype(np.float32)
    if not np.isfinite(svd.components_).all() or not np.isfinite(sketch.obsm["X_pca"]).all():
        raise RuntimeError("Training-sketch SVD produced non-finite values")

    # SVD组件和训练embedding是公开中间产物，保留原数组类型及文件名。
    np.savez_compressed(
        full_dir / "training_sketch_svd.npz",
        components = svd.components_,
        explained_variance_ratio = svd.explained_variance_ratio_,
    )
    embedding = ad.AnnData(X=sparse.csr_matrix((sketch.n_obs, 1)), obs=sketch.obs.copy())
    embedding.obsm["X_pca"] = sketch.obsm["X_pca"].copy()
    embedding.write_h5ad(full_dir / "training_sketch_embedding.h5ad", compression="gzip")

    # 固定谱系marker只用于谱系入口；字典顺序和所有基因保留。
    lineage_markers = {
        "Oligodendrocyte": ["MAG", "CLDN11", "APLP1", "TMEM144", "CNDP1", "EDIL3", "LARP6", "MOG", "AMER2", "TUBB4A"],
        "MES-like": ["VIM", "CD44", "CHI3L1", "HILPDA", "DDIT3", "ENO2"],
        "AC-like": ["GFAP", "S100B", "HOPX", "SLC1A3", "MLC1"],
        "OPC-like": ["PLP1", "OLIG1", "OMG", "TNR", "ALCAM"],
        "NPC-like": ["SOX4", "SOX11", "DCX", "CD24", "STMN1", "STMN2"],
        "NK/T": ["CD3D", "CD3E", "TRBC2", "CD3G", "CD2", "TRAC", "IL7R", "TRBC1", "LCK", "SKAP1", "CD48"],
        "Prol.NK/T": ["MKI67", "TOP2A", "CD3D", "CD3E", "TRBC2", "CD3G"],
        "Pericyte": ["COL1A2", "COL6A2", "PDGFRB", "COL3A1", "EDNRA", "LUM", "COL1A1", "PLAC9", "FRZB"],
        "Endothelial": ["ESAM", "FLT1", "RAMP2", "VWF", "EGFL7", "SLC9A3R2", "CAVIN2", "ADGRL4", "ABLIM1", "ABCB1"],
        "Pan-myeloid": ["PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "SPI1", "CSF1R", "LILRB1", "CTSS", "CYBA"],
        "Monocyte": ["FCN1", "VCAN", "S100A8", "S100A9"],
        "Microglia": ["TMEM119", "CX3CR1", "TREM2", "GPNMB", "CTSD", "CD68"],
        "Macrophage": ["SPP1", "APOE", "GPNMB", "CTSD", "CD68", "C1QA", "C1QB"],
        "CDC": ["FCER1A", "HLA-DPB1", "CD83", "AREG"],
        "Neutrophil": ["S100A8", "FCGR3B", "IFITM2"],
        "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD37"],
        "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
    }
    present = set(reference_genes)
    coverage_rows = []
    for name, marker_genes in lineage_markers.items():
        found = [gene for gene in marker_genes if gene in present]
        if len(found) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")
        coverage_rows.append(
            {
                "marker_set": name,
                "n_defined": len(marker_genes),
                "n_present": len(found),
                "genes_present": ";".join(found),
            }
        )
    coverage_table = pd.DataFrame(coverage_rows)
    coverage_table.to_csv(full_dir / "scgbm_lineage_marker_coverage.csv", index=False)

    gene_position = {gene: index for index, gene in enumerate(reference_genes)}
    hvg_indices = np.fromiter((gene_position[gene] for gene in hvg_genes), dtype=np.int64)
    metadata_blocks = []
    pca_blocks = []

    # 将全部singlet投影到同一sketch SVD；module分数使用全基因总计数，不以子集总计数替代。
    def module_scores(
        matrix: sparse.csr_matrix,
        totals: np.ndarray,
        genes: list[str],
        marker_sets: dict[str, list[str]],
    ) -> pd.DataFrame:
        """输入细胞×基因raw counts、全基因总计数及marker，返回每细胞一行的marker均值表。"""
        position = {gene: index for index, gene in enumerate(genes)}

        # marker并集保持字典遍历和基因首次出现顺序，平均分母不变。
        union_markers = []
        for marker_genes in marker_sets.values():
            union_markers.extend([gene for gene in marker_genes if gene in position])

        union_markers = list(dict.fromkeys(union_markers))
        union_indices = [position[gene] for gene in union_markers]
        normalized_union = data_functions.normalize_sparse(matrix[:, union_indices], totals)
        union_position = {gene: index for index, gene in enumerate(union_markers)}
        values = {}
        for name, marker_genes in marker_sets.items():
            indices = [union_position[gene] for gene in marker_genes if gene in union_position]
            values[f"raw_score_{name}"] = (
                np.asarray(normalized_union[:, indices].mean(axis=1)).ravel().astype(np.float32)
            )
        score_table = pd.DataFrame(values)
        return score_table

    for row in manifest_table.itertuples(index=False):
        matrix, genes, barcodes = data_functions.read_gse274546_library(row.library_id, raw_root)
        if genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")
        matrix, cell_metadata_table = data_functions.select_singlets(
            row.library_id, matrix, barcodes, classification_table
        )
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
        hvg = data_functions.normalize_sparse(matrix[:, hvg_indices], totals)
        pca_blocks.append(svd.transform(hvg).astype(np.float32))
        score_table = module_scores(matrix, totals, genes, lineage_markers)
        library_metadata_table = cell_metadata_table.reset_index(drop=True).copy()
        library_metadata_table.index = [f"{row.library_id}:{barcode}" for barcode in library_metadata_table["barcode"]]
        library_metadata_table.index.name = "cell_id"
        metadata_blocks.append(
            pd.concat([library_metadata_table, score_table.set_axis(library_metadata_table.index)], axis=1)
        )
    metadata_table_total = pd.concat(metadata_blocks, axis=0)
    pca = np.vstack(pca_blocks)
    if (
        len(metadata_table_total) != pca.shape[0]
        or metadata_table_total.index.duplicated().any()
        or not np.isfinite(pca).all()
    ):
        raise RuntimeError("Full-cell projection metadata is not aligned")

    # 每个raw_score在全体细胞内转z-score；保留float32、总体std及1e-8分母。
    raw_columns = [column for column in metadata_table_total if column.startswith("raw_score_")]
    for column in raw_columns:
        values = metadata_table_total[column].to_numpy(np.float32)
        scale = float(values.std())
        metadata_table_total[column.replace("raw_score_", "z_score_")] = (values - values.mean()) / (scale + 1e-8)
    del sketch

    # 全细胞主图按原capture Harmony、40PC、15邻居和三个分辨率依次计算。
    harmony = hm.run_harmony(
        pca,
        metadata_table_total,
        "capture_id",
        max_iter_harmony = 20,
        random_state = seed,
        verbose = False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == pca.shape:
        corrected = corrected.T
    elif corrected.shape != pca.shape:
        raise RuntimeError(f"Unexpected Harmony shape {corrected.shape}")
    embedding = ad.AnnData(X=sparse.csr_matrix((len(metadata_table_total), 1)), obs=metadata_table_total.copy())
    embedding.obsm["X_pca"] = pca
    embedding.obsm["X_pca_harmony"] = corrected.astype(np.float32)

    # 全细胞使用40PC和15邻居，Leiden主分辨率0.5，其余两个结果仍保留。
    n_neighbors = 15
    sc.pp.neighbors(
        embedding,
        n_neighbors  = n_neighbors,
        n_pcs        = full_n_pcs,
        use_rep      = "X_pca_harmony",
        random_state = seed,
    )
    sc.tl.umap(embedding, random_state=seed)
    for resolution in (0.5, 0.8, 1.0):
        sc.tl.leiden(
            embedding,
            resolution   = resolution,
            key_added    = f"leiden_full_{resolution}",
            random_state = seed,
            flavor       = "igraph",
            n_iterations = 2,
            directed     = False,
        )
    embedding.obs["leiden_full"] = embedding.obs["leiden_full_0.5"].astype("category")

    # 当前调用未传marker表；先按谱系z-score注释，Pan-myeloid入口保持mean z_score>0.5。
    cluster_key = "leiden_full"
    z_columns = [f"z_score_{name}" for name in lineage_markers]
    cluster_score_table = embedding.obs.groupby(cluster_key, observed=True)[z_columns].mean()
    cell_winner = embedding.obs[z_columns].idxmax(axis=1).str.replace("z_score_", "", regex=False)
    annotation_rows = []
    mapping = {}
    entry_mapping = {}
    for cluster in sorted(cluster_score_table.index.astype(str), key=int):
        row = cluster_score_table.loc[cluster]
        label = row.idxmax().replace("z_score_", "")
        mask = embedding.obs[cluster_key].astype(str).eq(cluster)
        support = float(cell_winner.loc[mask].eq(label).mean())
        mapping[cluster] = label
        is_pan_myeloid_candidate = bool(row["z_score_Pan-myeloid"] > 0.5)
        entry_mapping[cluster] = is_pan_myeloid_candidate
        annotation_rows.append(
            {
                "cluster": cluster,
                "lineage_label": label,
                "support_fraction": support,
                "n_cells": int(mask.sum()),
                "top20_markers": "",
                "is_pan_myeloid_candidate": is_pan_myeloid_candidate,
                "lineage_score_margin": float(np.sort(row.to_numpy())[-1] - np.sort(row.to_numpy())[-2]),
                **{column: float(row[column]) for column in z_columns},
            }
        )
    embedding.obs["lineage_label"] = embedding.obs[cluster_key].astype(str).map(mapping).astype("category")
    embedding.obs["is_pan_myeloid_candidate"] = embedding.obs[cluster_key].astype(str).map(entry_mapping).astype(bool)
    annotation_table = pd.DataFrame(annotation_rows)

    # 按原顺序逐文库汇总cluster counts/检出率，并仅保留已标记的pan-myeloid细胞。
    cluster_names = sorted(embedding.obs["leiden_full"].astype(str).unique(), key=int)
    cluster_position = {cluster: index for index, cluster in enumerate(cluster_names)}
    sums = np.zeros((len(cluster_names), len(reference_genes)), dtype=np.int64)
    detected = np.zeros_like(sums)
    cell_counts = np.zeros(len(cluster_names), dtype=np.int64)
    myeloid_blocks = []
    for row in manifest_table.itertuples(index=False):
        matrix, genes, barcodes = data_functions.read_gse274546_library(row.library_id, raw_root)
        if genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")
        matrix, cell_metadata_table = data_functions.select_singlets(
            row.library_id, matrix, barcodes, classification_table
        )
        cell_ids = [f"{row.library_id}:{barcode}" for barcode in cell_metadata_table["barcode"]]
        labels = embedding.obs.loc[cell_ids, "leiden_full"].astype(str).map(cluster_position).to_numpy(np.int32)
        indicator = sparse.csr_matrix(
            (np.ones(len(labels), dtype=np.int8), (labels, np.arange(len(labels)))),
            shape = (len(cluster_names), len(labels)),
        )
        sums += np.asarray((indicator @ matrix).todense(), dtype=np.int64)
        detected += np.asarray((indicator @ (matrix > 0)).todense(), dtype=np.int64)
        cell_counts += np.bincount(labels, minlength=len(cluster_names))
        keep = embedding.obs.loc[cell_ids, "is_pan_myeloid_candidate"].astype(bool).to_numpy()
        if keep.any():
            observation_table = embedding.obs.loc[pd.Index(cell_ids)[keep]].copy()
            myeloid_blocks.append(
                ad.AnnData(X=matrix[keep], obs=observation_table, var=pd.DataFrame(index=reference_genes))
            )

    # sums/detected为cluster×gene；按行总计数算CPM，rest为其余cluster的原始总和。
    totals = sums.sum(axis=1)
    cpm = np.divide(sums, totals[:, None], out=np.zeros_like(sums, dtype=np.float64), where=totals[:, None] > 0) * 1e6
    total_sums = sums.sum(axis=0)
    total_detected = detected.sum(axis=0)
    marker_rows = []
    for index, cluster in enumerate(cluster_names):
        rest_sums = total_sums - sums[index]
        rest_total = rest_sums.sum()
        rest_cpm = rest_sums / rest_total * 1e6 if rest_total > 0 else np.zeros_like(rest_sums, dtype=float)
        logfc = np.log2(cpm[index] + 1) - np.log2(rest_cpm + 1)
        pct_in = detected[index] / max(cell_counts[index], 1)
        rest_n = max(int(cell_counts.sum() - cell_counts[index]), 1)
        pct_rest = (total_detected - detected[index]) / rest_n
        rank_score = logfc * np.sqrt(np.clip(pct_in, 0, 1))
        order = np.argsort(-rank_score)[:100]

        # 排序指标仍为log2FC×sqrt(pct_in)，输出前100，不新增P值检验。
        for rank, gene_index in enumerate(order, start=1):
            marker_rows.append(
                {
                    "cluster": cluster,
                    "rank": rank,
                    "gene": reference_genes[gene_index],
                    "log2FC_cluster_vs_rest": logfc[gene_index],
                    "pct_in": pct_in[gene_index],
                    "pct_rest": pct_rest[gene_index],
                    "rank_score": rank_score[gene_index],
                }
            )
    if not myeloid_blocks:
        raise RuntimeError("No pan-myeloid cells were selected")
    myeloid = ad.concat(myeloid_blocks, join="inner", merge="same", index_unique=None)
    myeloid.X = myeloid.X.tocsr().astype(np.int32)
    if myeloid.obs_names.duplicated().any():
        raise RuntimeError("Duplicated pan-myeloid cell IDs")
    myeloid.layers["counts"] = myeloid.X.copy()
    full_marker_table = pd.DataFrame(marker_rows)

    # 导出原全细胞表和embedding；UMAP布局及点大小保持，PNG从同一PDF转换。
    full_marker_table.to_csv(
        full_dir / "full_cell_cluster_top100_aggregate_markers.csv.gz",
        index = False,
        compression = {"method": "gzip", "mtime": 0},
    )

    # 从实际聚合marker排序提取前20，补入注释表；cluster计数仍按原顺序转置。
    top20 = (
        full_marker_table.sort_values(["cluster", "rank"])
        .groupby("cluster", observed=True)["gene"]
        .apply(lambda values: ";".join(values.head(20)))
    )
    annotation_table["top20_markers"] = annotation_table["cluster"].map(top20)
    annotation_table.to_csv(full_dir / "full_cell_cluster_annotation.csv", index=False)
    cluster_count_table = pd.DataFrame(
        sums.T,
        index   = reference_genes,
        columns = [f"cluster_{x}" for x in sorted(embedding.obs["leiden_full"].astype(str).unique(), key=int)],
    )
    cluster_count_table.index.name = "gene"
    cluster_count_table.to_csv(
        full_dir / "full_cell_cluster_pseudobulk_counts.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    cluster_metadata_table = pd.DataFrame(
        {"cluster": sorted(embedding.obs["leiden_full"].astype(str).unique(), key=int), "n_cells": cell_counts}
    )
    cluster_metadata_table.to_csv(full_dir / "full_cell_cluster_pseudobulk_metadata.csv", index=False)

    # 全细胞元数据一行一个cell；原AnnData存盘后再绘制UMAP。
    full_cell_table = embedding.obs[
        [
            "barcode", "lineage_label", "is_pan_myeloid_candidate", "leiden_full", "condition", "capture_id",
            "patient_id",
        ]
    ].copy()
    full_cell_table["UMAP1"] = embedding.obsm["X_umap"][:, 0]
    full_cell_table["UMAP2"] = embedding.obsm["X_umap"][:, 1]
    full_cell_table.to_csv(full_dir / "full_cell_metadata.csv.gz", compression={"method": "gzip", "mtime": 0})
    embedding.write_h5ad(full_dir / "GSE274546_raw_independent_fullcell_embedding.h5ad", compression="gzip")
    figure_dir = root / "figures" / "36_gse274546_raw_independent_reannotation"
    figure_dir.mkdir(parents=True, exist_ok=True)

    # 同一UMAP图形逻辑复用三次；输入为细胞坐标表，返回值为空，仅保存图。
    def save_umap(cell_table: pd.DataFrame, label: str, path: Path, point_size: float) -> None:
        """输入含UMAP1/UMAP2和标签的细胞表；保存PDF及派生PNG，无返回值。"""
        labels = sorted(cell_table[label].astype(str).unique())
        colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(labels), 2)))
        fig, ax = plt.subplots(figsize=(6.0, 4.8))
        for value, color in zip(labels, colors):
            selected_table = cell_table.loc[cell_table[label].astype(str).eq(value)]
            ax.scatter(
                selected_table["UMAP1"],
                selected_table["UMAP2"],
                s          = point_size,
                color      = color,
                linewidths = 0,
                rasterized = True,
                label      = value,
            )

        # 轴、图例及点大小沿用原版；先draw，再保存PDF并转换400 dpi PNG。
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=4, fontsize=6)
        fig.canvas.draw()
        pdf_path = path.with_suffix(".pdf")
        fig.savefig(pdf_path, bbox_inches="tight")
        render_pdf_png(pdf_path, path.with_suffix(".png"), dpi=400)
        plt.close(fig)

    save_umap(full_cell_table, "lineage_label", figure_dir / "01_fullcell_umap_lineage", 0.18)

    # 髓系主聚类及Harmony敏感性顺序保留；污染top20>=3、线粒体/核糖体>=6仍为剔除界限。
    pan_myeloid_cells = myeloid.n_obs
    sc.pp.normalize_total(myeloid, target_sum=1e4)
    sc.pp.log1p(myeloid)
    sc.pp.highly_variable_genes(myeloid, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    hvg = myeloid[:, myeloid.var["highly_variable"].fillna(False)].copy()
    sc.pp.scale(hvg, max_value=10)
    myeloid_n_pcs = 30
    sc.tl.pca(hvg, n_comps=myeloid_n_pcs, svd_solver="arpack", random_state=seed)
    myeloid.obsm["X_pca"] = np.asarray(hvg.obsm["X_pca"], dtype=np.float32)
    sc.pp.neighbors(myeloid, n_neighbors=n_neighbors, n_pcs=myeloid_n_pcs, use_rep="X_pca", random_state=seed)
    sc.tl.umap(myeloid, random_state=seed)
    for resolution in (0.4, 0.6, 0.8):
        sc.tl.leiden(
            myeloid,
            resolution   = resolution,
            key_added    = f"leiden_myeloid_{resolution}",
            random_state = seed,
            flavor       = "igraph",
            n_iterations = 2,
            directed     = False,
        )
    myeloid.obs["leiden_myeloid"] = myeloid.obs["leiden_myeloid_0.6"].astype("category")

    harmony = hm.run_harmony(
        np.asarray(myeloid.obsm["X_pca"]),
        myeloid.obs,
        "capture_id",
        max_iter_harmony = 20,
        random_state = seed,
        verbose = False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == myeloid.obsm["X_pca"].shape:
        corrected = corrected.T
    elif corrected.shape != myeloid.obsm["X_pca"].shape:
        raise RuntimeError(f"Unexpected myeloid Harmony shape {corrected.shape}")
    myeloid.obsm["X_pca_harmony_sensitivity"] = corrected.astype(np.float32)
    sc.pp.neighbors(
        myeloid,
        n_neighbors  = n_neighbors,
        n_pcs        = myeloid_n_pcs,
        use_rep      = "X_pca_harmony_sensitivity",
        random_state = seed,
        key_added    = "harmony_neighbors",
    )
    sc.tl.umap(
        myeloid,
        random_state  = seed,
        neighbors_key = "harmony_neighbors",
        key_added     = "X_umap_harmony_sensitivity",
    )
    sc.tl.leiden(
        myeloid,
        resolution    = 0.6,
        key_added     = "leiden_myeloid_harmony_sensitivity",
        random_state  = seed,
        flavor        = "igraph",
        n_iterations  = 2,
        directed      = False,
        neighbors_key = "harmony_neighbors",
    )

    # 主聚类与Harmony敏感性聚类的ARI/NMI和交叉表继续供下游使用。
    myeloid_dir = base / "03_myeloid"
    myeloid_dir.mkdir(parents=True, exist_ok=True)
    harmony_agreement_table = pd.DataFrame(
        {
            "metric": ["adjusted_rand_index", "normalized_mutual_information"],
            "value": [
                adjusted_rand_score(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]),
                normalized_mutual_info_score(
                    myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]
                ),
            ],
        }
    )
    harmony_agreement_table.to_csv(myeloid_dir / "main_vs_harmony_cluster_agreement.csv", index=False)
    harmony_crosstab_table = pd.crosstab(
        myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]
    )
    harmony_crosstab_table.to_csv(myeloid_dir / "main_vs_harmony_cluster_crosstab.csv")

    # identity及污染marker评分每集合沿用seed=274546；集合顺序不变。
    myeloid_identity_labels = ("Monocyte", "Microglia", "Macrophage", "CDC", "Neutrophil")
    identity_sets = {name: lineage_markers[name] for name in myeloid_identity_labels}
    identity_coverage_table = data_functions.score_marker_sets(
        myeloid, identity_sets, "score_identity_", random_state=seed
    )
    identity_coverage_table.to_csv(myeloid_dir / "myeloid_identity_marker_coverage.csv", index=False)
    contamination_markers = {
        "B_cell": [
            "MS4A1", "CD79A", "CD79B", "BANK1", "CD37", "AFF3", "BLK", "BCL11A", "EBF1", "POU2AF1", "PAX5", "IKZF3",
            "RIPOR2",
        ],
        "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
        "T_NK_cell": [
            "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD2", "LCK", "NKG7", "GNLY", "KLRD1", "SKAP1", "FYN",
            "SLFN12L", "CD96", "THEMIS", "CD247", "BCL11B", "TOX", "ETS1", "ITK",
        ],
        "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5", "RAMP2", "FLT1", "ESAM"],
        "Neural_glial": [
            "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "TUBB3", "RBFOX3", "NLGN1", "LRP1B",
            "LSAMP", "PTPRZ1", "NPAS3", "DPP6", "NOVA1", "MAGI2", "NRCAM", "PCDH9", "SNTG1", "RBFOX1",
        ],
        "Stromal": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "AEBP1", "CCDC80", "COL14A1"],
        "Erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1", "GYPA"],
    }
    contamination_coverage_table = data_functions.score_marker_sets(
        myeloid, contamination_markers, "score_contam_", random_state=seed
    )
    contamination_coverage_table.to_csv(myeloid_dir / "myeloid_contamination_marker_coverage.csv", index=False)

    # t-test_overestim_var排序后使用原top20污染门槛，输出top100仍保留。
    myeloid_marker_table = data_functions.rank_clusters(myeloid, "leiden_myeloid", "rank_genes_myeloid")
    myeloid_marker_table.loc[myeloid_marker_table["rank"] <= 100].to_csv(
        myeloid_dir / "myeloid_cluster_top100_DEG.csv.gz",
        index = False,
        compression = {"method": "gzip", "mtime": 0},
    )
    top20 = myeloid_marker_table.loc[myeloid_marker_table["rank"] <= 20].groupby("cluster")["gene"].apply(list)
    identity_columns = [f"score_identity_{name}" for name in myeloid_identity_labels]
    contamination_columns = [f"score_contam_{name}" for name in contamination_markers]
    cluster_score_table = myeloid.obs.groupby("leiden_myeloid", observed=True)[
        identity_columns + contamination_columns
    ].mean()
    review_rows = []
    clean_clusters = []
    for cluster in cluster_score_table.index.astype(str):
        row = cluster_score_table.loc[cluster]
        identity = row[identity_columns].idxmax().replace("score_identity_", "")
        contamination = row[contamination_columns].idxmax().replace("score_contam_", "")
        genes = top20.get(cluster, [])
        hit_map = {
            name: [gene for gene in genes if gene in marker_set] for name, marker_set in contamination_markers.items()
        }
        winner = max(hit_map, key=lambda name: len(hit_map[name]))
        hits = len(hit_map[winner])
        mitochondrial_hits = sum(gene.startswith("MT-") for gene in genes)
        ribosomal_hits = sum(gene.startswith(("RPL", "RPS")) for gene in genes)
        clean = hits < 3 and mitochondrial_hits < 6 and ribosomal_hits < 6
        if clean:
            clean_clusters.append(cluster)
        review_rows.append(
            {
                "cluster": cluster,
                "n_cells": int(myeloid.obs["leiden_myeloid"].astype(str).eq(cluster).sum()),
                "identity": identity,
                "top_contamination": contamination,
                "top20_contamination_hits": hits,
                "top20_contamination_genes": ";".join(hit_map[winner]),
                "is_clean_myeloid": clean,
                "mitochondrial_top20_hits": mitochondrial_hits,
                "ribosomal_top20_hits": ribosomal_hits,
                "top20_markers": ";".join(genes),
            }
        )

    # 按数值cluster排序写复核表，再从is_clean_myeloid选择clean对象。
    review_table = pd.DataFrame(review_rows).sort_values("cluster", key=lambda x: x.astype(int))
    myeloid.obs["is_clean_myeloid"] = myeloid.obs["leiden_myeloid"].astype(str).isin(clean_clusters)
    myeloid.obs["myeloid_identity"] = (
        myeloid.obs["leiden_myeloid"].astype(str).map(review_table.set_index("cluster")["identity"])
    )
    review_table.to_csv(myeloid_dir / "myeloid_cluster_review.csv", index=False)
    clean = myeloid[myeloid.obs["is_clean_myeloid"].to_numpy()].copy()

    # clean myeloid导出及三种pseudobulk分组沿用原顺序；末尾汇总真实人数。
    clean.write_h5ad(myeloid_dir / "GSE274546_raw_independent_clean_myeloid.h5ad", compression="gzip")
    data_functions.write_pseudobulk(
        clean, ["patient_id", "pair_id", "condition"], "patient_condition_pseudobulk", myeloid_dir
    )
    data_functions.write_pseudobulk(
        clean,
        ["patient_id", "pair_id", "condition", "leiden_myeloid"],
        "patient_condition_cluster_pseudobulk",
        myeloid_dir,
    )
    data_functions.write_pseudobulk(
        clean,
        ["patient_id", "pair_id", "condition", "leiden_myeloid_harmony_sensitivity"],
        "patient_condition_harmony_cluster_pseudobulk",
        myeloid_dir,
    )
    clean_cell_table = clean.obs[
        ["barcode", "leiden_myeloid", "myeloid_identity", "condition", "capture_id", "patient_id"]
    ].copy()
    clean_cell_table["UMAP1"] = clean.obsm["X_umap"][:, 0]
    clean_cell_table["UMAP2"] = clean.obsm["X_umap"][:, 1]
    clean_cell_table.to_csv(
        myeloid_dir / "clean_myeloid_cell_metadata.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    save_umap(clean_cell_table, "leiden_myeloid", figure_dir / "02_myeloid_umap_cluster", 0.6)
    save_umap(clean_cell_table, "myeloid_identity", figure_dir / "03_myeloid_umap_identity", 0.6)

    # 每患者每时期计数，20/50/100/200四个原门槛只做原资格人数汇总。
    pair_count_table = clean.obs.groupby(["patient_id", "condition"], observed=True).size().unstack(fill_value=0)
    threshold_rows = []
    for threshold in (20, 50, 100, 200):
        eligible = pair_count_table.index[
            (pair_count_table.get("Primary", 0) >= threshold) & (pair_count_table.get("Recurrent", 0) >= threshold)
        ]
        threshold_rows.append({"min_cells_per_endpoint": threshold, "n_eligible_pairs": len(eligible)})
    threshold_table = pd.DataFrame(threshold_rows)
    threshold_table.to_csv(myeloid_dir / "paired_patient_threshold_audit.csv", index=False)

    # 52为原manifest声明的配对数；从输入重新计数写出，不能由印象填入审计字段。
    endpoint_table = pd.crosstab(manifest_table["patient_id"], manifest_table["condition"])
    n_paired_patients = int(((endpoint_table.get("Primary", 0) > 0) & (endpoint_table.get("Recurrent", 0) > 0)).sum())
    saved_sketch_table = pd.read_csv(full_dir / "training_sketch_by_library.csv")
    audit = {
        "dataset": "GSE274546",
        "input_libraries_primary_first_recurrent": len(manifest_table),
        "input_patients": int(manifest_table["patient_id"].nunique()),
        "paired_patients_before_cell_threshold": n_paired_patients,
        "singlet_cells": int(len(metadata_table_total)),
        "full_cell_training_sketch_cells": int(sum(saved_sketch_table["sketch_cells"])),
        "full_cell_hvg_method": "pooled seurat_v3 on a library-balanced sketch",
        "full_cell_hvg_batch_attempt": (
            "capture_id batch_key failed with a near-singular LOESS fit; no condition or patient input used"
        ),
        "full_cell_projection": (
            "all singlets projected to sketch-trained SVD; all-cell capture Harmony, kNN and Leiden"
        ),
        "full_cell_marker_method": "streamed raw-count cluster aggregation across all singlets",
        "pan_myeloid_entry": "full-cell clusters with mean z_score_Pan-myeloid > 0.5",
        "full_cell_clusters": int(embedding.obs["leiden_full"].nunique()),
        "pan_myeloid_cells": int(pan_myeloid_cells),
        "clean_myeloid_cells": int(clean.n_obs),
        "myeloid_clusters": int(clean.obs["leiden_myeloid"].nunique()),
        "activity_marker_sets_loaded_for_clustering": [],
        "external_validation_inputs_loaded": [],
    }
    (base / "raw_independent_reannotation_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"Singlets: {len(metadata_table_total)}; pan-myeloid: {pan_myeloid_cells}; clean myeloid: {clean.n_obs}")

    # 检查全部表格、矩阵和审计产物非空；不重新估计或替换任何结果。
    for path in (
        full_dir / "training_sketch_by_library.csv",
        full_dir / "training_sketch_hvg2000.csv",
        full_dir / "training_sketch_svd.npz",
        full_dir / "training_sketch_embedding.h5ad",
        full_dir / "GSE274546_raw_independent_fullcell_embedding.h5ad",
        full_dir / "scgbm_lineage_marker_coverage.csv",
        full_dir / "full_cell_cluster_top100_aggregate_markers.csv.gz",
        full_dir / "full_cell_cluster_annotation.csv",
        full_dir / "full_cell_cluster_pseudobulk_counts.csv.gz",
        full_dir / "full_cell_cluster_pseudobulk_metadata.csv",
        full_dir / "full_cell_metadata.csv.gz",
        myeloid_dir / "main_vs_harmony_cluster_agreement.csv",
        myeloid_dir / "main_vs_harmony_cluster_crosstab.csv",
        myeloid_dir / "myeloid_identity_marker_coverage.csv",
        myeloid_dir / "myeloid_contamination_marker_coverage.csv",
        myeloid_dir / "myeloid_cluster_top100_DEG.csv.gz",
        myeloid_dir / "myeloid_cluster_review.csv",
        myeloid_dir / "GSE274546_raw_independent_clean_myeloid.h5ad",
        myeloid_dir / "clean_myeloid_cell_metadata.csv.gz",
        myeloid_dir / "paired_patient_threshold_audit.csv",
        base / "raw_independent_reannotation_audit.json",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty reconstruction output: {path}")
