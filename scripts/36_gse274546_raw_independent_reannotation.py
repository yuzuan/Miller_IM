#!/usr/bin/env python3
"""GSE274546原始snRNA矩阵的独立全细胞重建和髓系重聚类。"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

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
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
RAW_ROOT = SOURCE_ROOT / "GSE274546"
MANIFEST = SOURCE_ROOT / "results" / "manifests" / "gse274546_standard_libraries.csv"
BASE = ROOT / "write" / "36_gse274546_raw_independent_reannotation"
QC_DIR = BASE / "01_qc_doublets"
FULL_DIR = BASE / "02_full_cell"
MYELOID_DIR = BASE / "03_myeloid"
FIGURE_DIR = ROOT / "figures" / "36_gse274546_raw_independent_reannotation"

SEED = 274546
SKETCH_PER_LIBRARY = 500
N_HVG = 2000
FULL_N_PCS = 40
MYELOID_N_PCS = 30
N_NEIGHBORS = 15
FULL_RESOLUTIONS = (0.5, 0.8, 1.0)
MYELOID_RESOLUTIONS = (0.4, 0.6, 0.8)
FULL_PRIMARY_RESOLUTION = 0.5
MYELOID_PRIMARY_RESOLUTION = 0.6

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
    "Pan-myeloid": ["PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "SPI1", "CSF1R", "LILRB1", "CTSS", "CYBA"],
    "Monocyte": ["FCN1", "VCAN", "S100A8", "S100A9"],
    "Microglia": ["TMEM119", "CX3CR1", "TREM2", "GPNMB", "CTSD", "CD68"],
    "Macrophage": ["SPP1", "APOE", "GPNMB", "CTSD", "CD68", "C1QA", "C1QB"],
    "CDC": ["FCER1A", "HLA-DPB1", "CD83", "AREG"],
    "Neutrophil": ["S100A8", "FCGR3B", "IFITM2"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD37"],
    "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
}

MYELOID_IDENTITY_LABELS = ("Monocyte", "Microglia", "Macrophage", "CDC", "Neutrophil")

CONTAMINATION_MARKERS = {
    "B_cell": [
        "MS4A1", "CD79A", "CD79B", "BANK1", "CD37", "AFF3", "BLK",
        "BCL11A", "EBF1", "POU2AF1", "PAX5", "IKZF3", "RIPOR2",
    ],
    "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
    "T_NK_cell": [
        "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD2", "LCK",
        "NKG7", "GNLY", "KLRD1", "SKAP1", "FYN", "SLFN12L", "CD96",
        "THEMIS", "CD247", "BCL11B", "TOX", "ETS1", "ITK",
    ],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5", "RAMP2", "FLT1", "ESAM"],
    "Neural_glial": [
        "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2",
        "TUBB3", "RBFOX3", "NLGN1", "LRP1B", "LSAMP", "PTPRZ1", "NPAS3", "DPP6",
        "NOVA1", "MAGI2", "NRCAM", "PCDH9", "SNTG1", "RBFOX1",
    ],
    "Stromal": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "AEBP1", "CCDC80", "COL14A1"],
    "Erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1", "GYPA"],
}


def make_unique(names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    result = []
    for name in names:
        if name not in counts:
            counts[name] = 0
            result.append(name)
        else:
            counts[name] += 1
            result.append(f"{name}-{counts[name]}")
    return result


def source_paths(library_id: str) -> tuple[Path, Path, Path]:
    library_dir = RAW_ROOT / library_id
    paths = (
        library_dir / "matrix.mtx.gz",
        library_dir / "genes.tsv.gz",
        library_dir / "barcodes.tsv.gz",
    )
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("; ".join(map(str, missing)))
    return paths


def read_library(library_id: str) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    matrix_path, genes_path, barcodes_path = source_paths(library_id)
    genes = pd.read_csv(genes_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr().astype(np.int32)
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(f"Matrix dimension mismatch for {library_id}: {matrix.shape}")
    return matrix.T.tocsr(), make_unique(genes), barcodes


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(MANIFEST)
    manifest = manifest.loc[manifest["progression"].isin(["Primary", "1st Recurrent"])].copy()
    manifest["pair_id"] = manifest["patient_id"]
    manifest["condition"] = np.where(manifest["progression"].eq("Primary"), "Primary", "Recurrent")
    manifest = manifest.sort_values(["patient_id", "progression_order", "library_id"]).reset_index(drop=True)
    if len(manifest) != 111 or manifest["patient_id"].nunique() != 59:
        raise ValueError("Expected 111 libraries from 59 patients")
    classifications = pd.read_csv(QC_DIR / "scDblFinder_cell_classifications.csv.gz", dtype=str)
    classifications = classifications.loc[classifications["scDblFinder_class"].eq("singlet")].copy()
    if classifications.duplicated(["capture_id", "barcode"]).any():
        raise ValueError("Duplicated singlet classifications")
    return manifest, classifications


def select_singlets(
    library_id: str,
    matrix: sparse.csr_matrix,
    barcodes: list[str],
    classifications: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    part = classifications.loc[classifications["capture_id"].eq(library_id)].copy()
    positions = {barcode: index for index, barcode in enumerate(barcodes)}
    missing = part.loc[~part["barcode"].isin(positions), "barcode"]
    if len(missing):
        raise ValueError(f"{library_id} has {len(missing)} classified barcodes absent from raw matrix")
    selected = np.fromiter((positions[x] for x in part["barcode"]), dtype=np.int64)
    return matrix[selected].tocsr(), part.reset_index(drop=True)


def normalize_sparse(matrix: sparse.csr_matrix, totals: np.ndarray) -> sparse.csr_matrix:
    factors = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    result = matrix.astype(np.float32).multiply(factors[:, None]).tocsr()
    np.log1p(result.data, out=result.data)
    return result


def marker_coverage(genes: list[str], marker_sets: dict[str, list[str]]) -> pd.DataFrame:
    present = set(genes)
    rows = []
    for name, markers in marker_sets.items():
        found = [gene for gene in markers if gene in present]
        if len(found) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")
        rows.append({"marker_set": name, "n_defined": len(markers), "n_present": len(found), "genes_present": ";".join(found)})
    return pd.DataFrame(rows)


def module_scores(
    matrix: sparse.csr_matrix,
    totals: np.ndarray,
    genes: list[str],
    marker_sets: dict[str, list[str]],
) -> pd.DataFrame:
    position = {gene: index for index, gene in enumerate(genes)}
    union_markers = list(dict.fromkeys(
        gene for markers in marker_sets.values() for gene in markers if gene in position
    ))
    union_indices = [position[gene] for gene in union_markers]
    normalized_union = normalize_sparse(matrix[:, union_indices], totals)
    union_position = {gene: index for index, gene in enumerate(union_markers)}
    values = {}
    for name, markers in marker_sets.items():
        indices = [union_position[gene] for gene in markers if gene in union_position]
        values[f"raw_score_{name}"] = np.asarray(
            normalized_union[:, indices].mean(axis=1)
        ).ravel().astype(np.float32)
    return pd.DataFrame(values)


def build_training_sketch(
    manifest: pd.DataFrame,
    classifications: pd.DataFrame,
) -> tuple[list[str], list[str], TruncatedSVD, ad.AnnData]:
    objects = []
    reference_genes: list[str] | None = None
    rng = np.random.default_rng(SEED)
    sketch_rows = []
    for row in manifest.itertuples(index=False):
        matrix, genes, barcodes = read_library(row.library_id)
        matrix, cell_meta = select_singlets(row.library_id, matrix, barcodes, classifications)
        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")
        n_keep = min(SKETCH_PER_LIBRARY, matrix.shape[0])
        selected = np.sort(rng.choice(matrix.shape[0], size=n_keep, replace=False))
        obs = cell_meta.iloc[selected].copy()
        obs.index = [f"{row.library_id}:{barcode}" for barcode in obs["barcode"]]
        obj = ad.AnnData(X=matrix[selected], obs=obs, var=pd.DataFrame(index=reference_genes))
        objects.append(obj)
        sketch_rows.append({"library_id": row.library_id, "singlets": matrix.shape[0], "sketch_cells": n_keep})
    if reference_genes is None:
        raise RuntimeError("No libraries were loaded")
    sketch = ad.concat(objects, join="inner", merge="same", index_unique=None)
    sketch.X = sketch.X.tocsr().astype(np.int32)
    sketch.obs_names_make_unique()
    pd.DataFrame(sketch_rows).to_csv(FULL_DIR / "training_sketch_by_library.csv", index=False)
    sc.pp.highly_variable_genes(
        sketch,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        batch_key=None,
        check_values=True,
    )
    hvg_genes = sketch.var_names[sketch.var["highly_variable"]].astype(str).tolist()
    if len(hvg_genes) != N_HVG:
        raise RuntimeError(f"Expected {N_HVG} HVGs, found {len(hvg_genes)}")
    pd.DataFrame({"gene": hvg_genes}).to_csv(FULL_DIR / "training_sketch_hvg2000.csv", index=False)
    sc.pp.normalize_total(sketch, target_sum=1e4)
    sc.pp.log1p(sketch)
    hvg_matrix = sketch[:, hvg_genes].X.tocsr()
    svd = TruncatedSVD(n_components=FULL_N_PCS, random_state=SEED)
    sketch.obsm["X_pca"] = svd.fit_transform(hvg_matrix).astype(np.float32)
    if not np.isfinite(svd.components_).all() or not np.isfinite(sketch.obsm["X_pca"]).all():
        raise RuntimeError("Training-sketch SVD produced non-finite values")
    np.savez_compressed(
        FULL_DIR / "training_sketch_svd.npz",
        components=svd.components_,
        explained_variance_ratio=svd.explained_variance_ratio_,
    )
    embedding = ad.AnnData(X=sparse.csr_matrix((sketch.n_obs, 1)), obs=sketch.obs.copy())
    embedding.obsm["X_pca"] = sketch.obsm["X_pca"].copy()
    embedding.write_h5ad(FULL_DIR / "training_sketch_embedding.h5ad", compression="gzip")
    return reference_genes, hvg_genes, svd, sketch


def project_all_cells(
    manifest: pd.DataFrame,
    classifications: pd.DataFrame,
    reference_genes: list[str],
    hvg_genes: list[str],
    svd: TruncatedSVD,
) -> tuple[pd.DataFrame, np.ndarray]:
    gene_position = {gene: index for index, gene in enumerate(reference_genes)}
    hvg_indices = np.fromiter((gene_position[gene] for gene in hvg_genes), dtype=np.int64)
    metadata_blocks = []
    pca_blocks = []
    for row in manifest.itertuples(index=False):
        matrix, genes, barcodes = read_library(row.library_id)
        if genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")
        matrix, cell_meta = select_singlets(row.library_id, matrix, barcodes, classifications)
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
        hvg = normalize_sparse(matrix[:, hvg_indices], totals)
        pca_blocks.append(svd.transform(hvg).astype(np.float32))
        scores = module_scores(matrix, totals, genes, SCGBM_LINEAGE_MARKERS)
        meta = cell_meta.reset_index(drop=True).copy()
        meta.index = [f"{row.library_id}:{barcode}" for barcode in meta["barcode"]]
        meta.index.name = "cell_id"
        metadata_blocks.append(pd.concat([meta, scores.set_axis(meta.index)], axis=1))
    metadata = pd.concat(metadata_blocks, axis=0)
    pca = np.vstack(pca_blocks)
    if len(metadata) != pca.shape[0] or metadata.index.duplicated().any() or not np.isfinite(pca).all():
        raise RuntimeError("Full-cell projection metadata is not aligned")
    raw_columns = [column for column in metadata if column.startswith("raw_score_")]
    for column in raw_columns:
        values = metadata[column].to_numpy(np.float32)
        scale = float(values.std())
        metadata[column.replace("raw_score_", "z_score_")] = (values - values.mean()) / (scale + 1e-8)
    return metadata, pca


def run_fullcell_clustering(metadata: pd.DataFrame, pca: np.ndarray) -> ad.AnnData:
    harmony = hm.run_harmony(
        pca,
        metadata,
        "capture_id",
        max_iter_harmony=20,
        random_state=SEED,
        verbose=False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == pca.shape:
        corrected = corrected.T
    elif corrected.shape != pca.shape:
        raise RuntimeError(f"Unexpected Harmony shape {corrected.shape}")
    embedding = ad.AnnData(X=sparse.csr_matrix((len(metadata), 1)), obs=metadata.copy())
    embedding.obsm["X_pca"] = pca
    embedding.obsm["X_pca_harmony"] = corrected.astype(np.float32)
    sc.pp.neighbors(
        embedding,
        n_neighbors=N_NEIGHBORS,
        n_pcs=FULL_N_PCS,
        use_rep="X_pca_harmony",
        random_state=SEED,
    )
    sc.tl.umap(embedding, random_state=SEED)
    for resolution in FULL_RESOLUTIONS:
        sc.tl.leiden(
            embedding,
            resolution=resolution,
            key_added=f"leiden_full_{resolution}",
            random_state=SEED,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    embedding.obs["leiden_full"] = embedding.obs[f"leiden_full_{FULL_PRIMARY_RESOLUTION}"].astype("category")
    return embedding


def aggregate_clusters(
    embedding: ad.AnnData,
    manifest: pd.DataFrame,
    classifications: pd.DataFrame,
    reference_genes: list[str],
    cluster_key: str,
    extract_myeloid: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, ad.AnnData | None]:
    cluster_names = sorted(embedding.obs[cluster_key].astype(str).unique(), key=int)
    cluster_position = {cluster: index for index, cluster in enumerate(cluster_names)}
    sums = np.zeros((len(cluster_names), len(reference_genes)), dtype=np.int64)
    detected = np.zeros_like(sums)
    cell_counts = np.zeros(len(cluster_names), dtype=np.int64)
    myeloid_blocks = []
    for row in manifest.itertuples(index=False):
        matrix, genes, barcodes = read_library(row.library_id)
        if genes != reference_genes:
            raise ValueError(f"Feature order differs for {row.library_id}")
        matrix, cell_meta = select_singlets(row.library_id, matrix, barcodes, classifications)
        cell_ids = [f"{row.library_id}:{barcode}" for barcode in cell_meta["barcode"]]
        labels = embedding.obs.loc[cell_ids, cluster_key].astype(str).map(cluster_position).to_numpy(np.int32)
        indicator = sparse.csr_matrix(
            (np.ones(len(labels), dtype=np.int8), (labels, np.arange(len(labels)))),
            shape=(len(cluster_names), len(labels)),
        )
        sums += np.asarray((indicator @ matrix).todense(), dtype=np.int64)
        detected += np.asarray((indicator @ (matrix > 0)).todense(), dtype=np.int64)
        cell_counts += np.bincount(labels, minlength=len(cluster_names))
        if extract_myeloid:
            keep = embedding.obs.loc[cell_ids, "is_pan_myeloid_candidate"].astype(bool).to_numpy()
            if keep.any():
                obs = embedding.obs.loc[pd.Index(cell_ids)[keep]].copy()
                myeloid_blocks.append(
                    ad.AnnData(X=matrix[keep], obs=obs, var=pd.DataFrame(index=reference_genes))
                )
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
        for rank, gene_index in enumerate(order, start=1):
            marker_rows.append({
                "cluster": cluster,
                "rank": rank,
                "gene": reference_genes[gene_index],
                "log2FC_cluster_vs_rest": logfc[gene_index],
                "pct_in": pct_in[gene_index],
                "pct_rest": pct_rest[gene_index],
                "rank_score": rank_score[gene_index],
            })
    myeloid = None
    if extract_myeloid:
        if not myeloid_blocks:
            raise RuntimeError("No pan-myeloid cells were selected")
        myeloid = ad.concat(myeloid_blocks, join="inner", merge="same", index_unique=None)
        myeloid.X = myeloid.X.tocsr().astype(np.int32)
        if myeloid.obs_names.duplicated().any():
            raise RuntimeError("Duplicated pan-myeloid cell IDs")
        myeloid.layers["counts"] = myeloid.X.copy()
    return sums, detected, cell_counts, pd.DataFrame(marker_rows), myeloid


def annotate_full_clusters(embedding: ad.AnnData, marker_table: pd.DataFrame | None = None) -> pd.DataFrame:
    cluster_key = "leiden_full"
    z_columns = [f"z_score_{name}" for name in SCGBM_LINEAGE_MARKERS]
    cluster_scores = embedding.obs.groupby(cluster_key, observed=True)[z_columns].mean()
    cell_winner = embedding.obs[z_columns].idxmax(axis=1).str.replace("z_score_", "", regex=False)
    rows = []
    mapping = {}
    entry_mapping = {}
    for cluster in sorted(cluster_scores.index.astype(str), key=int):
        row = cluster_scores.loc[cluster]
        label = row.idxmax().replace("z_score_", "")
        mask = embedding.obs[cluster_key].astype(str).eq(cluster)
        support = float(cell_winner.loc[mask].eq(label).mean())
        if marker_table is None:
            top_markers = ""
        else:
            top = marker_table.loc[marker_table["cluster"].eq(cluster)].sort_values("rank").head(20)
            top_markers = ";".join(top["gene"])
        mapping[cluster] = label
        is_pan_myeloid_candidate = bool(row["z_score_Pan-myeloid"] > 0.5)
        entry_mapping[cluster] = is_pan_myeloid_candidate
        rows.append({
            "cluster": cluster,
            "lineage_label": label,
            "support_fraction": support,
            "n_cells": int(mask.sum()),
            "top20_markers": top_markers,
            "is_pan_myeloid_candidate": is_pan_myeloid_candidate,
            "lineage_score_margin": float(np.sort(row.to_numpy())[-1] - np.sort(row.to_numpy())[-2]),
            **{column: float(row[column]) for column in z_columns},
        })
    embedding.obs["lineage_label"] = embedding.obs[cluster_key].astype(str).map(mapping).astype("category")
    embedding.obs["is_pan_myeloid_candidate"] = embedding.obs[cluster_key].astype(str).map(entry_mapping).astype(bool)
    return pd.DataFrame(rows)


def score_scanpy_sets(adata: ad.AnnData, marker_sets: dict[str, list[str]], prefix: str) -> pd.DataFrame:
    rows = []
    for name, markers in marker_sets.items():
        present = [gene for gene in markers if gene in adata.var_names]
        if len(present) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")
        sc.tl.score_genes(adata, present, score_name=f"{prefix}{name}", use_raw=False, random_state=SEED)
        rows.append({"marker_set": name, "n_defined": len(markers), "n_present": len(present), "genes_present": ";".join(present)})
    return pd.DataFrame(rows)


def rank_clusters(adata: ad.AnnData, cluster_key: str, key: str) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="t-test_overestim_var",
        use_raw=False,
        pts=True,
        key_added=key,
    )
    result = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    result = result.rename(columns={"group": "cluster", "names": "gene"})
    result["cluster"] = result["cluster"].astype(str)
    result["rank"] = result.groupby("cluster", observed=True).cumcount() + 1
    return result


def reprocess_myeloid(myeloid: ad.AnnData) -> tuple[ad.AnnData, pd.DataFrame]:
    sc.pp.normalize_total(myeloid, target_sum=1e4)
    sc.pp.log1p(myeloid)
    sc.pp.highly_variable_genes(myeloid, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
    hvg = myeloid[:, myeloid.var["highly_variable"].fillna(False)].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=MYELOID_N_PCS, svd_solver="arpack", random_state=SEED)
    myeloid.obsm["X_pca"] = np.asarray(hvg.obsm["X_pca"], dtype=np.float32)
    sc.pp.neighbors(myeloid, n_neighbors=N_NEIGHBORS, n_pcs=MYELOID_N_PCS, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(myeloid, random_state=SEED)
    for resolution in MYELOID_RESOLUTIONS:
        sc.tl.leiden(
            myeloid,
            resolution=resolution,
            key_added=f"leiden_myeloid_{resolution}",
            random_state=SEED,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    myeloid.obs["leiden_myeloid"] = myeloid.obs[f"leiden_myeloid_{MYELOID_PRIMARY_RESOLUTION}"].astype("category")

    harmony = hm.run_harmony(
        np.asarray(myeloid.obsm["X_pca"]),
        myeloid.obs,
        "capture_id",
        max_iter_harmony=20,
        random_state=SEED,
        verbose=False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == myeloid.obsm["X_pca"].shape:
        corrected = corrected.T
    elif corrected.shape != myeloid.obsm["X_pca"].shape:
        raise RuntimeError(f"Unexpected myeloid Harmony shape {corrected.shape}")
    myeloid.obsm["X_pca_harmony_sensitivity"] = corrected.astype(np.float32)
    sc.pp.neighbors(
        myeloid,
        n_neighbors=N_NEIGHBORS,
        n_pcs=MYELOID_N_PCS,
        use_rep="X_pca_harmony_sensitivity",
        random_state=SEED,
        key_added="harmony_neighbors",
    )
    sc.tl.umap(
        myeloid,
        random_state=SEED,
        neighbors_key="harmony_neighbors",
        key_added="X_umap_harmony_sensitivity",
    )
    sc.tl.leiden(
        myeloid,
        resolution=MYELOID_PRIMARY_RESOLUTION,
        key_added="leiden_myeloid_harmony_sensitivity",
        random_state=SEED,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        neighbors_key="harmony_neighbors",
    )
    pd.DataFrame({
        "metric": ["adjusted_rand_index", "normalized_mutual_information"],
        "value": [
            adjusted_rand_score(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]),
            normalized_mutual_info_score(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]),
        ],
    }).to_csv(MYELOID_DIR / "main_vs_harmony_cluster_agreement.csv", index=False)
    pd.crosstab(myeloid.obs["leiden_myeloid"], myeloid.obs["leiden_myeloid_harmony_sensitivity"]).to_csv(
        MYELOID_DIR / "main_vs_harmony_cluster_crosstab.csv"
    )

    identity_sets = {name: SCGBM_LINEAGE_MARKERS[name] for name in MYELOID_IDENTITY_LABELS}
    score_scanpy_sets(myeloid, identity_sets, "score_identity_").to_csv(MYELOID_DIR / "myeloid_identity_marker_coverage.csv", index=False)
    score_scanpy_sets(myeloid, CONTAMINATION_MARKERS, "score_contam_").to_csv(MYELOID_DIR / "myeloid_contamination_marker_coverage.csv", index=False)
    markers = rank_clusters(myeloid, "leiden_myeloid", "rank_genes_myeloid")
    markers.loc[markers["rank"] <= 100].to_csv(
        MYELOID_DIR / "myeloid_cluster_top100_DEG.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    top20 = markers.loc[markers["rank"] <= 20].groupby("cluster")["gene"].apply(list)
    identity_columns = [f"score_identity_{name}" for name in MYELOID_IDENTITY_LABELS]
    contamination_columns = [f"score_contam_{name}" for name in CONTAMINATION_MARKERS]
    cluster_scores = myeloid.obs.groupby("leiden_myeloid", observed=True)[identity_columns + contamination_columns].mean()
    review_rows = []
    clean_clusters = []
    for cluster in cluster_scores.index.astype(str):
        row = cluster_scores.loc[cluster]
        identity = row[identity_columns].idxmax().replace("score_identity_", "")
        contamination = row[contamination_columns].idxmax().replace("score_contam_", "")
        genes = top20.get(cluster, [])
        hit_map = {name: [gene for gene in genes if gene in marker_set] for name, marker_set in CONTAMINATION_MARKERS.items()}
        winner = max(hit_map, key=lambda name: len(hit_map[name]))
        hits = len(hit_map[winner])
        mitochondrial_hits = sum(gene.startswith("MT-") for gene in genes)
        ribosomal_hits = sum(gene.startswith(("RPL", "RPS")) for gene in genes)
        clean = hits < 3 and mitochondrial_hits < 6 and ribosomal_hits < 6
        if clean:
            clean_clusters.append(cluster)
        review_rows.append({
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
        })
    review = pd.DataFrame(review_rows).sort_values("cluster", key=lambda x: x.astype(int))
    myeloid.obs["is_clean_myeloid"] = myeloid.obs["leiden_myeloid"].astype(str).isin(clean_clusters)
    myeloid.obs["myeloid_identity"] = myeloid.obs["leiden_myeloid"].astype(str).map(review.set_index("cluster")["identity"])
    review.to_csv(MYELOID_DIR / "myeloid_cluster_review.csv", index=False)
    return myeloid[myeloid.obs["is_clean_myeloid"].to_numpy()].copy(), review


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


def save_umap(frame: pd.DataFrame, label: str, path: Path, point_size: float) -> None:
    labels = sorted(frame[label].astype(str).unique())
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(labels), 2)))
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    for value, color in zip(labels, colors):
        part = frame.loc[frame[label].astype(str).eq(value)]
        ax.scatter(part["UMAP1"], part["UMAP2"], s=point_size, color=color, linewidths=0, rasterized=True, label=value)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=4, fontsize=6)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    np.random.seed(SEED)
    for path in (FULL_DIR, MYELOID_DIR, FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    manifest, classifications = load_inputs()
    reference_genes, hvg_genes, svd, sketch = build_training_sketch(manifest, classifications)
    marker_coverage(reference_genes, SCGBM_LINEAGE_MARKERS).to_csv(FULL_DIR / "scgbm_lineage_marker_coverage.csv", index=False)
    metadata, pca = project_all_cells(manifest, classifications, reference_genes, hvg_genes, svd)
    del sketch
    embedding = run_fullcell_clustering(metadata, pca)
    annotation = annotate_full_clusters(embedding)
    sums, detected, cell_counts, full_markers, myeloid = aggregate_clusters(
        embedding, manifest, classifications, reference_genes, "leiden_full", extract_myeloid=True
    )
    if myeloid is None:
        raise RuntimeError("Pan-myeloid extraction did not return an AnnData object")
    full_markers.to_csv(FULL_DIR / "full_cell_cluster_top100_aggregate_markers.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    top20 = (
        full_markers.sort_values(["cluster", "rank"])
        .groupby("cluster", observed=True)["gene"]
        .apply(lambda values: ";".join(values.head(20)))
    )
    annotation["top20_markers"] = annotation["cluster"].map(top20)
    annotation.to_csv(FULL_DIR / "full_cell_cluster_annotation.csv", index=False)
    cluster_counts = pd.DataFrame(sums.T, index=reference_genes, columns=[f"cluster_{x}" for x in sorted(embedding.obs["leiden_full"].astype(str).unique(), key=int)])
    cluster_counts.index.name = "gene"
    cluster_counts.to_csv(FULL_DIR / "full_cell_cluster_pseudobulk_counts.csv.gz", compression={"method": "gzip", "mtime": 0})
    pd.DataFrame({"cluster": sorted(embedding.obs["leiden_full"].astype(str).unique(), key=int), "n_cells": cell_counts}).to_csv(
        FULL_DIR / "full_cell_cluster_pseudobulk_metadata.csv", index=False
    )
    full_frame = embedding.obs[[
        "barcode", "lineage_label", "is_pan_myeloid_candidate", "leiden_full",
        "condition", "capture_id", "patient_id",
    ]].copy()
    full_frame["UMAP1"] = embedding.obsm["X_umap"][:, 0]
    full_frame["UMAP2"] = embedding.obsm["X_umap"][:, 1]
    full_frame.to_csv(FULL_DIR / "full_cell_metadata.csv.gz", compression={"method": "gzip", "mtime": 0})
    embedding.write_h5ad(FULL_DIR / "GSE274546_raw_independent_fullcell_embedding.h5ad", compression="gzip")
    save_umap(full_frame, "lineage_label", FIGURE_DIR / "01_fullcell_umap_lineage", 0.18)

    pan_myeloid_cells = myeloid.n_obs
    clean, review = reprocess_myeloid(myeloid)
    clean.write_h5ad(MYELOID_DIR / "GSE274546_raw_independent_clean_myeloid.h5ad", compression="gzip")
    write_pseudobulk(clean, ["patient_id", "pair_id", "condition"], "patient_condition_pseudobulk")
    write_pseudobulk(clean, ["patient_id", "pair_id", "condition", "leiden_myeloid"], "patient_condition_cluster_pseudobulk")
    write_pseudobulk(
        clean,
        ["patient_id", "pair_id", "condition", "leiden_myeloid_harmony_sensitivity"],
        "patient_condition_harmony_cluster_pseudobulk",
    )
    clean_frame = clean.obs[["barcode", "leiden_myeloid", "myeloid_identity", "condition", "capture_id", "patient_id"]].copy()
    clean_frame["UMAP1"] = clean.obsm["X_umap"][:, 0]
    clean_frame["UMAP2"] = clean.obsm["X_umap"][:, 1]
    clean_frame.to_csv(MYELOID_DIR / "clean_myeloid_cell_metadata.csv.gz", compression={"method": "gzip", "mtime": 0})
    save_umap(clean_frame, "leiden_myeloid", FIGURE_DIR / "02_myeloid_umap_cluster", 0.6)
    save_umap(clean_frame, "myeloid_identity", FIGURE_DIR / "03_myeloid_umap_identity", 0.6)

    pair_counts = clean.obs.groupby(["patient_id", "condition"], observed=True).size().unstack(fill_value=0)
    threshold_rows = []
    for threshold in (20, 50, 100, 200):
        eligible = pair_counts.index[(pair_counts.get("Primary", 0) >= threshold) & (pair_counts.get("Recurrent", 0) >= threshold)]
        threshold_rows.append({"min_cells_per_endpoint": threshold, "n_eligible_pairs": len(eligible)})
    pd.DataFrame(threshold_rows).to_csv(MYELOID_DIR / "paired_patient_threshold_audit.csv", index=False)

    audit = {
        "dataset": "GSE274546",
        "input_libraries_primary_first_recurrent": len(manifest),
        "input_patients": int(manifest["patient_id"].nunique()),
        "paired_patients_before_cell_threshold": 52,
        "singlet_cells": int(len(metadata)),
        "full_cell_training_sketch_cells": int(sum(pd.read_csv(FULL_DIR / "training_sketch_by_library.csv")["sketch_cells"])),
        "full_cell_hvg_method": "pooled seurat_v3 on a library-balanced sketch",
        "full_cell_hvg_batch_attempt": "capture_id batch_key failed with a near-singular LOESS fit; no condition or patient input used",
        "full_cell_projection": "all singlets projected to sketch-trained SVD; all-cell capture Harmony, kNN and Leiden",
        "full_cell_marker_method": "streamed raw-count cluster aggregation across all singlets",
        "pan_myeloid_entry": "full-cell clusters with mean z_score_Pan-myeloid > 0.5",
        "full_cell_clusters": int(embedding.obs["leiden_full"].nunique()),
        "pan_myeloid_cells": int(pan_myeloid_cells),
        "clean_myeloid_cells": int(clean.n_obs),
        "myeloid_clusters": int(clean.obs["leiden_myeloid"].nunique()),
        "activity_marker_sets_loaded_for_clustering": [],
        "external_validation_inputs_loaded": [],
    }
    (BASE / "raw_independent_reannotation_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
