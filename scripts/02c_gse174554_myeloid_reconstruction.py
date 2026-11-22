#!/usr/bin/env python3
"""Rebuild all-cell and clean-myeloid objects from GSE174554 raw counts and scDblFinder singlets.
Inputs under MILLER_IM_DATA_ROOT: GSE174554/GSE174554_RAW and results/manifests/gse174554_libraries.csv.
    write/34_gse174554_raw_independent_discovery/01_qc_doublets/scDblFinder_cell_classifications.csv.gz.
Outputs under the same write module: 02_full_cell coverage, annotations, top 50 DE genes and all-cell H5AD;
    03_myeloid cluster review, top 100 DE genes, Harmony comparison, clean H5AD, metadata and three pseudobulks;
    raw_independent_reconstruction_audit.json at the module root. The published file names are retained.
Run with python scripts/02c_gse174554_myeloid_reconstruction.py.
The existing Python workflow uses Scanpy, AnnData, harmonypy and SciPy with the original QC, markers and seeds.
"""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path

import anndata as ad
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # Use published paths, sorted GSMs and library order to keep the original resequencing merge.
    root = Path(__file__).resolve().parents[1]
    source_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
    manifest_path = source_root / "results" / "manifests" / "gse174554_libraries.csv"
    manifest_table = pd.read_csv(manifest_path)

    # The frozen public manifest requires 91 matrices and 81 GSM accessions.
    expected_library_count = 91
    expected_sample_count = 81
    if len(manifest_table) != expected_library_count or manifest_table["gsm_id"].nunique() != expected_sample_count:
        raise ValueError(
            f"GSE174554 manifest does not contain {expected_library_count} matrices and {expected_sample_count} GSMs"
        )

    analysis_root = root / "write" / "34_gse174554_raw_independent_discovery"
    classification_path = analysis_root / "01_qc_doublets" / "scDblFinder_cell_classifications.csv.gz"
    classification_table = pd.read_csv(classification_path, dtype=str)
    classification_table = classification_table.loc[classification_table["scDblFinder_class"].eq("singlet")].copy()
    captures = []
    reference_features = None
    raw_root = source_root / "GSE174554" / "GSE174554_RAW"

    # Read a manifest row into cell-by-gene CSR counts, genes and barcodes. Stop on missing files or size mismatch.

    for gsm_id, library_table in manifest_table.groupby("gsm_id", sort=True):
        capture_payloads = []
        if len(library_table) == 1:
            row = library_table.iloc[0]
            capture_payloads.append((str(row["library_id"]), row, data_functions.read_gse174554_library(row, raw_root)))
        else:
            main_row = library_table.loc[library_table["tech_rep"].eq("main")].iloc[0]
            batch_row = library_table.loc[library_table["tech_rep"].eq("batch2")].iloc[0]
            main_payload = data_functions.read_gse174554_library(main_row, raw_root)
            batch_payload = data_functions.read_gse174554_library(batch_row, raw_root)
            main_matrix, main_features, main_barcodes = main_payload
            batch_matrix, batch_features, batch_barcodes = batch_payload
            overlap_fraction = len(set(main_barcodes).intersection(batch_barcodes)) / min(
                len(main_barcodes), len(batch_barcodes)
            )

            # Overlap >=0.80 denotes resequencing and <=0.05 denotes separate captures. Stop between those limits.
            if overlap_fraction >= 0.80 and main_features != batch_features:
                raise ValueError("Feature order differs within resequenced sample")

            if overlap_fraction >= 0.80:
                union_barcodes = list(dict.fromkeys(main_barcodes + batch_barcodes))
                positions = {barcode: index for index, barcode in enumerate(union_barcodes)}
                source_matrix = sparse.vstack([main_matrix, batch_matrix], format="csr")
                union_indices = np.fromiter(
                    (positions[barcode] for barcode in main_barcodes + batch_barcodes), dtype=np.int64
                )
                mapper = sparse.csr_matrix(
                    (np.ones(len(union_indices), dtype=np.int8), (union_indices, np.arange(len(union_indices)))),
                    shape = (len(union_barcodes), len(union_indices)),
                )
                merged_matrix = (mapper @ source_matrix).tocsr().astype(np.int32)
                capture_payloads.append(
                    (f"{gsm_id}__merged_resequencing", main_row, (merged_matrix, main_features, union_barcodes))
                )
            elif overlap_fraction <= 0.05:
                capture_payloads.append((str(main_row["library_id"]), main_row, main_payload))
                capture_payloads.append((str(batch_row["library_id"]), batch_row, batch_payload))
            else:
                raise RuntimeError(f"Ambiguous barcode overlap for {gsm_id}: {overlap_fraction:.3f}")

        # Keep only the raw matrices for the current GSM and collect singlets in the original order.
        for capture_id, row, payload in capture_payloads:
            matrix, features, barcodes = payload
            if reference_features is None:
                reference_features = features
            elif features != reference_features:
                raise ValueError(f"Feature order differs for {capture_id}")

            singlets = classification_table.loc[classification_table["capture_id"].eq(capture_id), "barcode"].tolist()
            if not singlets:
                continue

            barcode_positions = {barcode: index for index, barcode in enumerate(barcodes)}
            missing_barcodes = [barcode for barcode in singlets if barcode not in barcode_positions]
            if missing_barcodes:
                raise ValueError(
                    f"{capture_id} has {len(missing_barcodes)} classification barcodes absent from raw matrix"
                )

            selected_indices = np.fromiter((barcode_positions[barcode] for barcode in singlets), dtype=np.int64)
            selected_matrix = matrix[selected_indices].tocsr()
            observation_table = pd.DataFrame(index=[f"{capture_id}:{barcode}" for barcode in singlets])
            observation_table["capture_id"] = capture_id
            for column in ["gsm_id", "sample_id", "sample_short", "patient_id", "pair_id", "progression", "tech_rep"]:
                observation_table[column] = row[column]

            observation_table["condition"] = "Primary" if row["progression"] == "Primary" else "Recurrent"
            observation_table["barcode"] = singlets
            feature_table = pd.DataFrame(index=pd.Index(data_functions.make_unique(features), name="gene"))
            captures.append(ad.AnnData(X=selected_matrix, obs=observation_table, var=feature_table))

    # Keep the original random-call order through all-cell processing, clustering and myeloid review.
    adata = ad.concat(captures, join="inner", merge="same", index_unique=None)
    adata.X = adata.X.tocsr().astype(np.int32)
    adata.layers["counts"] = adata.X.copy()
    adata.obs_names_make_unique()
    if adata.n_obs == 0 or adata.obs_names.duplicated().any():
        raise ValueError("No singlets or duplicated cell ids after reconstruction")

    # Keep raw counts and use the frozen 1e4 normalization, 2,000 HVGs, scaling limit of 10 and 40 PCs.
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
        sc.tl.leiden(
            adata,
            resolution   = resolution,
            key_added    = f"leiden_{resolution}",
            random_state = 0,
            flavor       = "igraph",
            n_iterations = 2,
            directed     = False,
        )
    adata.obs["leiden"] = adata.obs["leiden_0.5"].astype("category")

    # Use the original lineage markers and published case-sensitive score suffixes.
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
        "Monocyte": ["FCN1", "VCAN", "S100A8", "S100A9"],
        "Microglia": ["TMEM119", "CX3CR1", "TREM2", "GPNMB", "CTSD", "CD68"],
        "Macrophage": ["SPP1", "APOE", "GPNMB", "CTSD", "CD68", "C1QA", "C1QB"],
        "CDC": ["FCER1A", "HLA-DPB1", "CD83", "AREG"],
        "Neutrophil": ["S100A8", "FCGR3B", "IFITM2"],
    }

    coverage_table = data_functions.score_marker_sets(adata, lineage_markers, "score_lineage_", random_state=0)

    # Mark close score winners Ambiguous using the 0.05 margin. The cluster vote uses all cells as its denominator.
    score_columns = [f"score_lineage_{name}" for name in lineage_markers]
    score_table = adata.obs[score_columns]
    values = np.sort(score_table.to_numpy(float), axis=1)
    margin = values[:, -1] - values[:, -2]
    predicted = score_table.idxmax(axis=1).str.replace("score_lineage_", "", regex=False)
    predicted.loc[margin < 0.05] = "Ambiguous"
    adata.obs["predicted_lineage"] = predicted.astype("category")

    # Count non-Ambiguous votes within each cluster. With no valid vote, record Unknown and zero support.
    annotation_rows = []
    mapping = {}
    for cluster in sorted(adata.obs["leiden"].astype(str).unique(), key=int):
        mask = adata.obs["leiden"].astype(str).eq(cluster)
        votes = predicted.loc[mask]
        non_ambiguous = votes.loc[~votes.eq("Ambiguous")]
        if len(non_ambiguous):
            counts = non_ambiguous.value_counts()
            label = str(counts.index[0])
            support = float(counts.iloc[0] / mask.sum())
        else:
            label, support = "Unknown", 0.0
        mapping[cluster] = label
        annotation_rows.append(
            {
                "cluster": cluster,
                "lineage_label": label,
                "support_fraction": support,
                "n_cells": int(mask.sum()),
                "n_ambiguous": int(votes.eq("Ambiguous").sum()),
            }
        )
    adata.obs["lineage_label"] = adata.obs["leiden"].astype(str).map(mapping).astype("category")
    annotation_table = pd.DataFrame(annotation_rows)

    # Keep t-test_overestim_var for DE and export the top 50 plus the complete all-cell AnnData.
    marker_table = data_functions.rank_clusters(adata, "leiden", "rank_genes_full")
    full_dir = analysis_root / "02_full_cell"
    full_dir.mkdir(parents=True, exist_ok=True)
    coverage_table.to_csv(full_dir / "scgbm_lineage_marker_coverage.csv", index=False)
    annotation_table.to_csv(full_dir / "full_cell_cluster_annotation.csv", index=False)
    marker_table.loc[marker_table["rank"] <= 50].to_csv(
        full_dir / "full_cell_cluster_top50_DEG.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    adata.write_h5ad(full_dir / "GSE174554_raw_independent_fullcell_annotated.h5ad", compression="gzip")

    # Select the five original myeloid lineages and recluster with 2,000 HVGs and 30 PCs.
    myeloid_labels = ("Monocyte", "Microglia", "Macrophage", "CDC", "Neutrophil")
    myeloid = adata[adata.obs["lineage_label"].astype(str).isin(myeloid_labels)].copy()
    sc.pp.highly_variable_genes(myeloid, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    myeloid_hvg = myeloid[:, myeloid.var["highly_variable"].fillna(False)].copy()
    sc.pp.scale(myeloid_hvg, max_value=10)
    sc.tl.pca(myeloid_hvg, n_comps=30, svd_solver="arpack", random_state=0)
    myeloid.obsm["X_pca"] = np.asarray(myeloid_hvg.obsm["X_pca"], dtype=np.float32)
    sc.pp.neighbors(myeloid, n_neighbors=15, n_pcs=30, use_rep="X_pca", random_state=0)
    sc.tl.umap(myeloid, random_state=0)
    for resolution in (0.4, 0.6, 0.8):
        sc.tl.leiden(
            myeloid,
            resolution   = resolution,
            key_added    = f"leiden_myeloid_{resolution}",
            random_state = 0,
            flavor       = "igraph",
            n_iterations = 2,
            directed     = False,
        )
    myeloid.obs["leiden_myeloid"] = myeloid.obs["leiden_myeloid_0.6"].astype("category")

    # Keep the Harmony sensitivity clustering with 20 iterations, capture_id batches and seed 0.
    harmony = hm.run_harmony(
        np.asarray(myeloid.obsm["X_pca"]),
        myeloid.obs,
        "capture_id",
        max_iter_harmony = 20,
        random_state = 0,
        verbose = False,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.T.shape == myeloid.obsm["X_pca"].shape:
        corrected = corrected.T
    elif corrected.shape != myeloid.obsm["X_pca"].shape:
        raise ValueError(f"Unexpected Harmony shape: {corrected.shape}")
    myeloid.obsm["X_pca_harmony_sensitivity"] = corrected.astype(np.float32)
    sc.pp.neighbors(
        myeloid,
        n_neighbors  = 15,
        n_pcs        = 30,
        use_rep      = "X_pca_harmony_sensitivity",
        random_state = 0,
        key_added    = "harmony_neighbors",
    )
    sc.tl.umap(myeloid, random_state=0, neighbors_key="harmony_neighbors", key_added="X_umap_harmony_sensitivity")
    sc.tl.leiden(
        myeloid,
        resolution    = 0.6,
        key_added     = "leiden_myeloid_harmony_sensitivity",
        random_state  = 0,
        flavor        = "igraph",
        n_iterations  = 2,
        directed      = False,
        neighbors_key = "harmony_neighbors",
    )

    # Save the original ARI, NMI and cross-tabulation for downstream tables and figures.
    myeloid_dir = analysis_root / "03_myeloid"
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

    # Keep the original identity and contamination markers, score prefixes and seed 0 for each set.
    identity_coverage_table = data_functions.score_marker_sets(
        myeloid, {name: lineage_markers[name] for name in myeloid_labels}, "score_identity_", random_state=0
    )
    identity_coverage_table.to_csv(myeloid_dir / "myeloid_identity_marker_coverage.csv", index=False)
    contamination_markers = {
        "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD37"],
        "Plasma_cell": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "SDC1"],
        "T_NK_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "NKG7", "GNLY", "KLRD1"],
        "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5", "RAMP2", "FLT1", "ESAM"],
        "Neural_glial": [
            "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "TUBB3", "RBFOX3", "NLGN1", "LRP1B",
            "LSAMP", "PTPRZ1", "NPAS3", "DPP6", "NOVA1", "MAGI2", "NRCAM", "PCDH9", "SNTG1", "RBFOX1",
        ],
        "Stromal": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "AEBP1", "CCDC80", "COL14A1"],
        "Erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1", "GYPA"],
    }
    contamination_coverage_table = data_functions.score_marker_sets(
        myeloid, contamination_markers, "score_contam_", random_state=0
    )
    contamination_coverage_table.to_csv(myeloid_dir / "myeloid_contamination_marker_coverage.csv", index=False)

    # Check the top 20 markers. Exclude at least 3 of one contaminating type, 6 mitochondrial or 6 ribosomal genes.
    myeloid_marker_table = data_functions.rank_clusters(myeloid, "leiden_myeloid", "rank_genes_myeloid")
    top20 = myeloid_marker_table.loc[myeloid_marker_table["rank"] <= 20].groupby("cluster")["gene"].apply(list)
    identity_columns = [f"score_identity_{name}" for name in myeloid_labels]
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
            name: [gene for gene in genes if gene in markers_] for name, markers_ in contamination_markers.items()
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
                "n_cells": int((myeloid.obs["leiden_myeloid"].astype(str) == cluster).sum()),
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

    # Write the review in numeric cluster order, then select clean myeloid cells using the frozen labels.
    review_table = pd.DataFrame(review_rows).sort_values("cluster", key=lambda x: x.astype(int))
    myeloid.obs["is_clean_myeloid"] = myeloid.obs["leiden_myeloid"].astype(str).isin(clean_clusters)
    myeloid.obs["myeloid_identity"] = (
        myeloid.obs["leiden_myeloid"].astype(str).map(review_table.set_index("cluster")["identity"])
    )
    review_table.to_csv(myeloid_dir / "myeloid_cluster_review.csv", index=False)
    myeloid_marker_table.loc[myeloid_marker_table["rank"] <= 100].to_csv(
        myeloid_dir / "myeloid_cluster_top100_DEG.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    clean = myeloid[myeloid.obs["is_clean_myeloid"].to_numpy()].copy()
    clean.write_h5ad(myeloid_dir / "GSE174554_raw_independent_clean_myeloid.h5ad", compression="gzip")
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

    # Write one metadata row per cell with primary UMAP coordinates and count the actual objects for the audit.
    cell_metadata_table = clean.obs.reset_index(names="cell_id")
    cell_metadata_table["UMAP1"] = clean.obsm["X_umap"][:, 0]
    cell_metadata_table["UMAP2"] = clean.obsm["X_umap"][:, 1]
    cell_metadata_table.to_csv(
        myeloid_dir / "clean_myeloid_cell_metadata.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    audit = {
        "dataset": "GSE174554",
        "raw_matrices": len(manifest_table),
        "gsm": int(manifest_table["gsm_id"].nunique()),
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
    (analysis_root / "raw_independent_reconstruction_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(f"Singlets: {adata.n_obs}; pan-myeloid: {myeloid.n_obs}; clean myeloid: {clean.n_obs}")

    # Require all expected downstream files to exist and be nonempty.
    figure_dir = root / "figures" / "34_gse174554_raw_independent_discovery"
    figure_dir.mkdir(parents=True, exist_ok=True)
    for path in (
        full_dir / "scgbm_lineage_marker_coverage.csv",
        full_dir / "full_cell_cluster_annotation.csv",
        full_dir / "full_cell_cluster_top50_DEG.csv.gz",
        full_dir / "GSE174554_raw_independent_fullcell_annotated.h5ad",
        myeloid_dir / "main_vs_harmony_cluster_agreement.csv",
        myeloid_dir / "main_vs_harmony_cluster_crosstab.csv",
        myeloid_dir / "myeloid_identity_marker_coverage.csv",
        myeloid_dir / "myeloid_contamination_marker_coverage.csv",
        myeloid_dir / "myeloid_cluster_top100_DEG.csv.gz",
        myeloid_dir / "GSE174554_raw_independent_clean_myeloid.h5ad",
        myeloid_dir / "clean_myeloid_cell_metadata.csv.gz",
        myeloid_dir / "myeloid_cluster_review.csv",
        analysis_root / "raw_independent_reconstruction_audit.json",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty reconstruction output: {path}")
