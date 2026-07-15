#!/usr/bin/env python3
"""GSE174554原始计数矩阵的独立全细胞注释与pan-myeloid重聚类。"""

from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import anndata as ad
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
RAW_ROOT = SOURCE_ROOT / "GSE174554" / "GSE174554_RAW"
MANIFEST = SOURCE_ROOT / "results" / "manifests" / "gse174554_libraries.csv"
ANALYSIS_ROOT = ROOT / "write" / "34_gse174554_raw_independent_discovery"
QC_DIR = ANALYSIS_ROOT / "01_qc_doublets"
FULL_DIR = ANALYSIS_ROOT / "02_full_cell"
MYELOID_DIR = ANALYSIS_ROOT / "03_myeloid"
FIGURE_DIR = ROOT / "figures" / "34_gse174554_raw_independent_discovery"

RESEQUENCING_OVERLAP_MIN = 0.80
INDEPENDENT_OVERLAP_MAX = 0.05
FORCE_MANUAL_MATRIX_LIBRARY = "GSM5319529_SF8963"
MYELOID_LABELS = ("Monocyte", "Microglia", "Macrophage", "CDC", "Neutrophil")

# 与scGBM正式全细胞粗注释一致，只用于谱系入口，不用于复发程序选基因。
SCGBM_LINEAGE_MARKERS = {
    "Oligodendrocyte": ["MAG", "CLDN11", "APLP1", "TMEM144", "CNDP1", "EDIL3", "LARP6", "MOG", "AMER2", "TUBB4A"],
    "MES-like": ["VIM", "CD44", "CHI3L1", "HILPDA", "DDIT3", "ENO2"],
    "AC-like": ["GFAP", "S100B", "HOPX", "SLC1A3", "MLC1"],
    "OPC-like": ["PLP1", "OLIG1", "OMG", "TNR", "ALCAM"],
    "NPC-like": ["SOX4", "SOX11", "DCX", "CD24", "STMN1", "STMN2"],
    "NK/T": ["CD3D", "CD3E", "TRBC2", "CD3G", "CD2", "TRAC", "IL7R", "TRBC1", "LCK", "SKAP1", "CD48"],
    "Prol.NK/T": ["MKI67", "TOP2A", "CD3D", "CD3E", "TRBC2", "CD3G"],
    "Pericyte": ["COL1A2", "COL6A2", "PDGFRB", "COL3A1", "EDNRA", "LUM", "COL1A1", "PLAC9", "FRZB"],
    "Endothelial": ["ESAM", "FLT1", "RAMP2", "VWF", "EGFL7", "SLC9A3R2", "CAVIN2", "ADGRL4", "ABLIM1", "ABCB1"],
    "Monocyte": ["FCN1", "VCAN", "S100A8", "S100A9"],
    "Microglia": ["TMEM119", "CX3CR1", "TREM2", "GPNMB", "CTSD", "CD68"],
    "Macrophage": ["SPP1", "APOE", "GPNMB", "CTSD", "CD68", "C1QA", "C1QB"],
    "CDC": ["FCER1A", "HLA-DPB1", "CD83", "AREG"],
    "Neutrophil": ["S100A8", "FCGR3B", "IFITM2"],
}

CONTAMINATION_MARKERS = {
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD37"],
    "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
    "T_NK_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "NKG7", "GNLY", "KLRD1"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5", "RAMP2", "FLT1", "ESAM"],
    "Neural_glial": ["GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "TUBB3", "RBFOX3", "NLGN1", "LRP1B", "LSAMP", "PTPRZ1", "NPAS3", "DPP6", "NOVA1", "MAGI2", "NRCAM", "PCDH9", "SNTG1", "RBFOX1"],
    "Stromal": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "AEBP1", "CCDC80", "COL14A1"],
    "Erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1", "GYPA"],
}


def make_unique(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    result = []
    for name in names:
        if name not in seen:
            seen[name] = 0
            result.append(name)
        else:
            seen[name] += 1
            result.append(f"{name}-{seen[name]}")
    return result


def read_matrix_manual(path: Path) -> sparse.csr_matrix:
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                n_rows, n_cols, declared_nnz = map(int, line.split())
                break
        rows = np.empty(declared_nnz, dtype=np.int32)
        cols = np.empty(declared_nnz, dtype=np.int32)
        data = np.empty(declared_nnz, dtype=np.int32)
        cursor = 0
        for line in handle:
            values = re.findall(r"\d+", line)
            if len(values) < 3:
                continue
            rows[cursor] = int(values[0]) - 1
            cols[cursor] = int(values[1]) - 1
            data[cursor] = int(values[2])
            cursor += 1
    return sparse.coo_matrix(
        (data[:cursor], (rows[:cursor], cols[:cursor])), shape=(n_rows, n_cols), dtype=np.int32
    ).tocsr()


def source_path(value: str) -> Path:
    path = RAW_ROOT / Path(value).name
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_library(row: pd.Series) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    matrix_path = source_path(row["matrix_path"])
    features_path = source_path(row["features_path"])
    barcodes_path = source_path(row["barcodes_path"])
    features = pd.read_csv(features_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    if str(row["library_id"]) == FORCE_MANUAL_MATRIX_LIBRARY:
        matrix = read_matrix_manual(matrix_path)
    else:
        with gzip.open(matrix_path, "rb") as handle:
            matrix = mmread(handle).tocsr().astype(np.int32)
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(f"Matrix shape mismatch for {row['library_id']}: {matrix.shape}")
    return matrix.T.tocsr(), features, barcodes


def merge_resequencing(
    main: tuple[sparse.csr_matrix, list[str], list[str]],
    batch2: tuple[sparse.csr_matrix, list[str], list[str]],
) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    main_matrix, features, main_barcodes = main
    batch_matrix, batch_features, batch_barcodes = batch2
    if features != batch_features:
        raise ValueError("Feature order differs within resequenced sample")
    union = list(dict.fromkeys(main_barcodes + batch_barcodes))
    position = {barcode: index for index, barcode in enumerate(union)}
    source = sparse.vstack([main_matrix, batch_matrix], format="csr")
    union_index = np.fromiter((position[x] for x in main_barcodes + batch_barcodes), dtype=np.int64)
    mapper = sparse.csr_matrix(
        (np.ones(len(union_index), dtype=np.int8), (union_index, np.arange(len(union_index)))),
        shape=(len(union), len(union_index)),
    )
    return (mapper @ source).tocsr().astype(np.int32), features, union


def build_capture_objects() -> list[ad.AnnData]:
    manifest = pd.read_csv(MANIFEST)
    if len(manifest) != 91 or manifest["gsm_id"].nunique() != 81:
        raise ValueError("GSE174554 manifest does not contain 91 matrices and 81 GSMs")
    classifications = pd.read_csv(QC_DIR / "scDblFinder_cell_classifications.csv.gz", dtype=str)
    classifications = classifications.loc[classifications["scDblFinder_class"].eq("singlet")].copy()
    captures: list[ad.AnnData] = []
    reference_features: list[str] | None = None

    def append_capture(capture_id: str, row: pd.Series, payload: tuple[sparse.csr_matrix, list[str], list[str]]) -> None:
        nonlocal reference_features
        matrix, features, barcodes = payload
        if reference_features is None:
            reference_features = features
        elif features != reference_features:
            raise ValueError(f"Feature order differs for {capture_id}")
        singlets = classifications.loc[classifications["capture_id"].eq(capture_id), "barcode"].tolist()
        if not singlets:
            return
        barcode_position = {barcode: index for index, barcode in enumerate(barcodes)}
        missing = [barcode for barcode in singlets if barcode not in barcode_position]
        if missing:
            raise ValueError(f"{capture_id} has {len(missing)} classification barcodes absent from raw matrix")
        positions = np.fromiter((barcode_position[x] for x in singlets), dtype=np.int64)
        selected = matrix[positions].tocsr()
        obs = pd.DataFrame(index=[f"{capture_id}:{barcode}" for barcode in singlets])
        obs["capture_id"] = capture_id
        for column in ["gsm_id", "sample_id", "sample_short", "patient_id", "pair_id", "progression", "tech_rep"]:
            obs[column] = row[column]
        obs["condition"] = "Primary" if row["progression"] == "Primary" else "Recurrent"
        obs["barcode"] = singlets
        var = pd.DataFrame(index=pd.Index(make_unique(features), name="gene"))
        obj = ad.AnnData(X=selected, obs=obs, var=var)
        captures.append(obj)

    for gsm_id, rows in manifest.groupby("gsm_id", sort=True):
        if len(rows) == 1:
            row = rows.iloc[0]
            append_capture(str(row["library_id"]), row, read_library(row))
            continue
        main_row = rows.loc[rows["tech_rep"].eq("main")].iloc[0]
        batch_row = rows.loc[rows["tech_rep"].eq("batch2")].iloc[0]
        main = read_library(main_row)
        batch = read_library(batch_row)
        overlap_fraction = len(set(main[2]).intersection(batch[2])) / min(len(main[2]), len(batch[2]))
        if overlap_fraction >= RESEQUENCING_OVERLAP_MIN:
            append_capture(f"{gsm_id}__merged_resequencing", main_row, merge_resequencing(main, batch))
        elif overlap_fraction <= INDEPENDENT_OVERLAP_MAX:
            append_capture(str(main_row["library_id"]), main_row, main)
            append_capture(str(batch_row["library_id"]), batch_row, batch)
        else:
            raise RuntimeError(f"Ambiguous barcode overlap for {gsm_id}: {overlap_fraction:.3f}")
    return captures


def score_marker_sets(adata: ad.AnnData, marker_sets: dict[str, list[str]], prefix: str) -> pd.DataFrame:
    rows = []
    for name, genes in marker_sets.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if len(present) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")
        sc.tl.score_genes(adata, present, score_name=f"{prefix}{name}", use_raw=False, random_state=0)
        rows.append({"marker_set": name, "n_defined": len(genes), "n_present": len(present), "genes_present": ";".join(present)})
    return pd.DataFrame(rows)


def cluster_annotation(adata: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    score_columns = [f"score_lineage_{name}" for name in SCGBM_LINEAGE_MARKERS]
    scores = adata.obs[score_columns]
    values = np.sort(scores.to_numpy(float), axis=1)
    margin = values[:, -1] - values[:, -2]
    predicted = scores.idxmax(axis=1).str.replace("score_lineage_", "", regex=False)
    predicted.loc[margin < 0.05] = "Ambiguous"
    adata.obs["predicted_lineage"] = predicted.astype("category")
    rows = []
    mapping = {}
    for cluster in sorted(adata.obs[cluster_key].astype(str).unique(), key=int):
        mask = adata.obs[cluster_key].astype(str).eq(cluster)
        votes = predicted.loc[mask]
        non_ambiguous = votes.loc[~votes.eq("Ambiguous")]
        if len(non_ambiguous):
            counts = non_ambiguous.value_counts()
            label = str(counts.index[0])
            support = float(counts.iloc[0] / mask.sum())
        else:
            label, support = "Unknown", 0.0
        mapping[cluster] = label
        rows.append({"cluster": cluster, "lineage_label": label, "support_fraction": support, "n_cells": int(mask.sum()), "n_ambiguous": int(votes.eq("Ambiguous").sum())})
    adata.obs["lineage_label"] = adata.obs[cluster_key].astype(str).map(mapping).astype("category")
    return pd.DataFrame(rows)


def rank_clusters(adata: ad.AnnData, cluster_key: str, key: str) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="t-test_overestim_var", use_raw=False, pts=True, key_added=key)
    result = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    result = result.rename(columns={"group": "cluster", "names": "gene"})
    result["cluster"] = result["cluster"].astype(str)
    result["rank"] = result.groupby("cluster", observed=True).cumcount() + 1
    return result


def write_pseudobulk(adata: ad.AnnData, group_columns: list[str], prefix: str) -> None:
    obs = adata.obs[group_columns].copy()
    obs["row"] = np.arange(adata.n_obs)
    matrix = adata.layers["counts"].tocsr()
    vectors = []
    metadata = []
    for index, (keys, part) in enumerate(obs.groupby(group_columns, observed=True, sort=True)):
        if not isinstance(keys, tuple):
            keys = (keys,)
        sample_id = f"PB{index + 1:04d}"
        vector = np.asarray(matrix[part["row"].to_numpy()].sum(axis=0)).ravel().astype(np.int64)
        vectors.append(pd.Series(vector, name=sample_id))
        row = dict(zip(group_columns, keys))
        row.update(sample_id=sample_id, n_cells=len(part), library_size=int(vector.sum()))
        metadata.append(row)
    counts = pd.concat(vectors, axis=1)
    counts.insert(0, "gene", adata.var_names.astype(str))
    counts.to_csv(MYELOID_DIR / f"{prefix}_counts.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    pd.DataFrame(metadata).to_csv(MYELOID_DIR / f"{prefix}_metadata.csv", index=False)


def main() -> None:
    for path in (FULL_DIR, MYELOID_DIR, FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)
    captures = build_capture_objects()
    adata = ad.concat(captures, join="inner", merge="same", index_unique=None)
    adata.X = adata.X.tocsr().astype(np.int32)
    adata.layers["counts"] = adata.X.copy()
    adata.obs_names_make_unique()
    if adata.n_obs == 0 or adata.obs_names.duplicated().any():
        raise ValueError("No singlets or duplicated cell ids after reconstruction")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    hvg = adata[:, adata.var["highly_variable"].fillna(False)].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=40, svd_solver="arpack", random_state=0)
    adata.obsm["X_pca"] = np.asarray(hvg.obsm["X_pca"], dtype=np.float32)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep="X_pca", random_state=0)
    sc.tl.umap(adata, random_state=0)
    for resolution in (0.5, 0.8, 1.0):
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}", random_state=0, flavor="igraph", n_iterations=2, directed=False)
    adata.obs["leiden"] = adata.obs["leiden_0.5"].astype("category")

    coverage = score_marker_sets(adata, SCGBM_LINEAGE_MARKERS, "score_lineage_")
    annotation = cluster_annotation(adata, "leiden")
    markers = rank_clusters(adata, "leiden", "rank_genes_full")
    coverage.to_csv(FULL_DIR / "scgbm_lineage_marker_coverage.csv", index=False)
    annotation.to_csv(FULL_DIR / "full_cell_cluster_annotation.csv", index=False)
    markers.loc[markers["rank"] <= 50].to_csv(FULL_DIR / "full_cell_cluster_top50_DEG.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    adata.write_h5ad(FULL_DIR / "GSE174554_raw_independent_fullcell_annotated.h5ad", compression="gzip")

    myeloid = adata[adata.obs["lineage_label"].astype(str).isin(MYELOID_LABELS)].copy()
    sc.pp.highly_variable_genes(myeloid, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    myeloid_hvg = myeloid[:, myeloid.var["highly_variable"].fillna(False)].copy()
    sc.pp.scale(myeloid_hvg, max_value=10)
    sc.tl.pca(myeloid_hvg, n_comps=30, svd_solver="arpack", random_state=0)
    myeloid.obsm["X_pca"] = np.asarray(myeloid_hvg.obsm["X_pca"], dtype=np.float32)
    sc.pp.neighbors(myeloid, n_neighbors=15, n_pcs=30, use_rep="X_pca", random_state=0)
    sc.tl.umap(myeloid, random_state=0)
    for resolution in (0.4, 0.6, 0.8):
        sc.tl.leiden(myeloid, resolution=resolution, key_added=f"leiden_myeloid_{resolution}", random_state=0, flavor="igraph", n_iterations=2, directed=False)
    myeloid.obs["leiden_myeloid"] = myeloid.obs["leiden_myeloid_0.6"].astype("category")

    harmony = hm.run_harmony(
        np.asarray(myeloid.obsm["X_pca"]), myeloid.obs, "capture_id",
        max_iter_harmony=20, random_state=0, verbose=False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == myeloid.obsm["X_pca"].shape:
        corrected = corrected.T
    elif corrected.shape != myeloid.obsm["X_pca"].shape:
        raise ValueError(f"Unexpected Harmony shape: {corrected.shape}")
    myeloid.obsm["X_pca_harmony_sensitivity"] = corrected.astype(np.float32)
    sc.pp.neighbors(
        myeloid, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony_sensitivity",
        random_state=0, key_added="harmony_neighbors",
    )
    sc.tl.umap(myeloid, random_state=0, neighbors_key="harmony_neighbors", key_added="X_umap_harmony_sensitivity")
    sc.tl.leiden(
        myeloid, resolution=0.6, key_added="leiden_myeloid_harmony_sensitivity",
        random_state=0, flavor="igraph", n_iterations=2, directed=False,
        neighbors_key="harmony_neighbors",
    )
    pd.DataFrame({
        "metric": ["adjusted_rand_index", "normalized_mutual_information"],
        "value": [
            adjusted_rand_score(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]),
            normalized_mutual_info_score(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]),
        ],
    }).to_csv(MYELOID_DIR / "main_vs_harmony_cluster_agreement.csv", index=False)
    pd.crosstab(
        myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]
    ).to_csv(MYELOID_DIR / "main_vs_harmony_cluster_crosstab.csv")
    score_marker_sets(myeloid, {name: SCGBM_LINEAGE_MARKERS[name] for name in MYELOID_LABELS}, "score_identity_").to_csv(
        MYELOID_DIR / "myeloid_identity_marker_coverage.csv", index=False
    )
    score_marker_sets(myeloid, CONTAMINATION_MARKERS, "score_contam_").to_csv(
        MYELOID_DIR / "myeloid_contamination_marker_coverage.csv", index=False
    )
    myeloid_markers = rank_clusters(myeloid, "leiden_myeloid", "rank_genes_myeloid")
    top20 = myeloid_markers.loc[myeloid_markers["rank"] <= 20].groupby("cluster")["gene"].apply(list)
    identity_columns = [f"score_identity_{name}" for name in MYELOID_LABELS]
    contamination_columns = [f"score_contam_{name}" for name in CONTAMINATION_MARKERS]
    cluster_scores = myeloid.obs.groupby("leiden_myeloid", observed=True)[identity_columns + contamination_columns].mean()
    review_rows = []
    clean_clusters = []
    for cluster in cluster_scores.index.astype(str):
        row = cluster_scores.loc[cluster]
        identity = row[identity_columns].idxmax().replace("score_identity_", "")
        contamination = row[contamination_columns].idxmax().replace("score_contam_", "")
        genes = top20.get(cluster, [])
        hit_map = {name: [gene for gene in genes if gene in markers_] for name, markers_ in CONTAMINATION_MARKERS.items()}
        winner = max(hit_map, key=lambda name: len(hit_map[name]))
        hits = len(hit_map[winner])
        mitochondrial_hits = sum(gene.startswith("MT-") for gene in genes)
        ribosomal_hits = sum(gene.startswith(("RPL", "RPS")) for gene in genes)
        clean = hits < 3 and mitochondrial_hits < 6 and ribosomal_hits < 6
        if clean:
            clean_clusters.append(cluster)
        review_rows.append({
            "cluster": cluster, "n_cells": int((myeloid.obs["leiden_myeloid"].astype(str) == cluster).sum()),
            "identity": identity, "top_contamination": contamination, "top20_contamination_hits": hits,
            "top20_contamination_genes": ";".join(hit_map[winner]), "is_clean_myeloid": clean,
            "mitochondrial_top20_hits": mitochondrial_hits, "ribosomal_top20_hits": ribosomal_hits,
            "top20_markers": ";".join(genes),
        })
    review = pd.DataFrame(review_rows).sort_values("cluster", key=lambda x: x.astype(int))
    myeloid.obs["is_clean_myeloid"] = myeloid.obs["leiden_myeloid"].astype(str).isin(clean_clusters)
    myeloid.obs["myeloid_identity"] = myeloid.obs["leiden_myeloid"].astype(str).map(review.set_index("cluster")["identity"])
    review.to_csv(MYELOID_DIR / "myeloid_cluster_review.csv", index=False)
    myeloid_markers.loc[myeloid_markers["rank"] <= 100].to_csv(
        MYELOID_DIR / "myeloid_cluster_top100_DEG.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    clean = myeloid[myeloid.obs["is_clean_myeloid"].to_numpy()].copy()
    clean.write_h5ad(MYELOID_DIR / "GSE174554_raw_independent_clean_myeloid.h5ad", compression="gzip")
    write_pseudobulk(clean, ["patient_id", "pair_id", "condition"], "patient_condition_pseudobulk")
    write_pseudobulk(clean, ["patient_id", "pair_id", "condition", "leiden_myeloid"], "patient_condition_cluster_pseudobulk")
    write_pseudobulk(
        clean, ["patient_id", "pair_id", "condition", "leiden_myeloid_harmony_sensitivity"],
        "patient_condition_harmony_cluster_pseudobulk",
    )

    cell_metadata = clean.obs.reset_index(names="cell_id")
    cell_metadata["UMAP1"] = clean.obsm["X_umap"][:, 0]
    cell_metadata["UMAP2"] = clean.obsm["X_umap"][:, 1]
    cell_metadata.to_csv(MYELOID_DIR / "clean_myeloid_cell_metadata.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    audit = {
        "dataset": "GSE174554",
        "raw_matrices": 91,
        "gsm": 81,
        "singlet_cells": int(adata.n_obs),
        "full_cell_clusters": int(adata.obs["leiden"].nunique()),
        "pan_myeloid_cells": int(myeloid.n_obs),
        "clean_myeloid_cells": int(clean.n_obs),
        "myeloid_clusters": int(myeloid.obs["leiden_myeloid"].nunique()),
        "batch_correction_primary_clustering": "none",
        "harmony_sensitivity_batch_key": "capture_id",
        "hvg_method": "seurat_v3 on pooled GSE174554 singlets; no condition or patient input",
        "activity_marker_sets_loaded": [],
        "external_validation_inputs_loaded": [],
    }
    (ANALYSIS_ROOT / "raw_independent_reconstruction_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
