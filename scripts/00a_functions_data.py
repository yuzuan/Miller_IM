#!/usr/bin/env python3
# ======================================================================
# 00a_functions_data.py — Miller-IM shared data functions
# Public data reconstruction and sample manifests for paired glioblastomas.
#
# Purpose
# Share file readers, marker scoring, and pseudobulk aggregation across the data workflows.
# The analysis scripts show how each helper is used.
# Use go-to-definition in the editor to read the implementation from its call site.
# Data and analysis settings are passed by the callers.

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    import anndata as ad
    import numpy as np
    from scipy import sparse


# Append serial suffixes -1, -2, and so on to duplicated names.
def make_unique(names: list[str]) -> list[str]:
    """Return unique names in the same order as the input."""
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


# Read repaired Matrix Market counts with the published record parser and sum repeated coordinates.
def read_matrix_manual(path: Path) -> tuple[sparse.csr_matrix, int, int]:
    """Return the int32 gene-by-cell CSR matrix and the declared and parsed entry counts."""
    import numpy as np
    from scipy import sparse

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                n_rows, n_columns, declared_nnz = map(int, line.split())
                break

        # Allocate from the header and retain only the records that were actually read.
        row_indices = np.empty(declared_nnz, dtype=np.int32)
        column_indices = np.empty(declared_nnz, dtype=np.int32)
        count_values = np.empty(declared_nnz, dtype=np.int32)
        cursor = 0
        for line in handle:
            values = re.findall(r"\d+", line)
            if len(values) < 3:
                continue

            row_indices[cursor] = int(values[0]) - 1
            column_indices[cursor] = int(values[1]) - 1
            count_values[cursor] = int(values[2])
            cursor += 1

    matrix = sparse.coo_matrix(
        (count_values[:cursor], (row_indices[:cursor], column_indices[:cursor])),
        shape = (n_rows, n_columns),
        dtype = np.int32,
    )
    matrix = matrix.tocsr()
    result = matrix, declared_nnz, cursor
    return result


# Score the supplied marker sets in order and record their measured gene coverage.
def score_marker_sets(
    adata: ad.AnnData,
    marker_sets: dict[str, list[str]],
    prefix: str,
    random_state: int,
) -> pd.DataFrame:
    """Save scores to AnnData and return one coverage row per marker set."""
    import scanpy as sc

    coverage_rows = []
    for name, genes in marker_sets.items():
        present_genes = [gene for gene in genes if gene in adata.var_names]
        if len(present_genes) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")

        # Repeat the caller's seed for each set to match the published control-gene draws.
        sc.tl.score_genes(adata, present_genes, score_name=f"{prefix}{name}", use_raw=False, random_state=random_state)
        coverage_rows.append(
            {
                "marker_set": name,
                "n_defined": len(genes),
                "n_present": len(present_genes),
                "genes_present": ";".join(present_genes),
            }
        )

    coverage_table = pd.DataFrame(coverage_rows)
    return coverage_table


# Rank the detected genes for each cluster using the published Scanpy test.
def rank_clusters(adata: ad.AnnData, cluster_key: str, key: str) -> pd.DataFrame:
    """Save the ranking and return cluster, gene, and rank from t-test_overestim_var."""
    import scanpy as sc

    sc.tl.rank_genes_groups(
        adata, groupby=cluster_key, method="t-test_overestim_var", use_raw=False, pts=True, key_added=key
    )
    marker_table = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    marker_table = marker_table.rename(columns={"group": "cluster", "names": "gene"})
    marker_table["cluster"] = marker_table["cluster"].astype(str)
    marker_table["rank"] = marker_table.groupby("cluster", observed=True).cumcount() + 1
    return marker_table


# Aggregate the raw counts by the supplied fields and write the count and sample tables.
def write_pseudobulk(adata: ad.AnnData, group_columns: list[str], prefix: str, output_dir: Path) -> None:
    """Write a gene-by-sample count matrix and matching sample metadata in group order."""
    import numpy as np

    observation_table = adata.obs[group_columns].copy()
    observation_table["row"] = np.arange(adata.n_obs)
    count_matrix = adata.layers["counts"].tocsr()
    vectors = []
    metadata_rows = []
    for index, (group_values, group_table) in enumerate(
        observation_table.groupby(group_columns, observed=True, sort=True)
    ):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)

        # Sum the cells within each gene and keep the original integer dtype and sample numbering.
        sample_id = f"PB{index + 1:04d}"
        vector = np.asarray(count_matrix[group_table["row"].to_numpy()].sum(axis=0)).ravel().astype(np.int64)
        vectors.append(pd.Series(vector, name=sample_id))
        metadata_row = dict(zip(group_columns, group_values))
        metadata_row.update(sample_id=sample_id, n_cells=len(group_table), library_size=int(vector.sum()))
        metadata_rows.append(metadata_row)

    # Stack the sample vectors in gene order and create the output directory before saving.
    count_table_total = pd.concat(vectors, axis=1)
    count_table_total.insert(0, "gene", adata.var_names.astype(str))
    metadata_table = pd.DataFrame(metadata_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    count_path = output_dir / f"{prefix}_counts.csv.gz"
    count_table_total.to_csv(count_path, index=False, compression={"method": "gzip", "mtime": 0})
    metadata_path = output_dir / f"{prefix}_metadata.csv"
    metadata_table.to_csv(metadata_path, index=False)

    # Verify that the pseudobulk files for the downstream analyses have matching sample dimensions.
    if count_table_total.shape != (adata.n_vars, len(metadata_table) + 1):
        raise RuntimeError("Pseudobulk counts and metadata dimensions differ")

    for path in (count_path, metadata_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty pseudobulk output: {path}")


# Calculate the original float32 library-size normalization and sparse log1p transform.
def normalize_sparse(matrix: sparse.csr_matrix, totals: np.ndarray) -> sparse.csr_matrix:
    """Normalize to 10000 per cell using the supplied totals from the full count matrix."""
    import numpy as np

    factors = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    normalized_matrix = matrix.astype(np.float32).multiply(factors[:, None]).tocsr()
    np.log1p(normalized_matrix.data, out=normalized_matrix.data)
    return normalized_matrix


# Read the dimensions and declared nonzero entries from a gzip Matrix Market header.
def parse_matrix_header(matrix_path: Path) -> tuple[int, int, int]:
    """Return row count, column count, and declared entries from the first non-comment line."""
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue

            n_rows, n_columns, nnz = (int(value) for value in line.strip().split())
            dimensions = n_rows, n_columns, nnz
            return dimensions

    raise ValueError(f"Could not parse Matrix Market header: {matrix_path}")


# Standardize known recurrence-stage labels from the GEO metadata.
def normalize_progression(value: str | None) -> object:
    """Return the published stage label after trimming surrounding spaces."""
    if value is None or pd.isna(value):
        result = pd.NA
        return result

    value = str(value).strip()
    progression_labels = {
        "Primary": "Primary",
        "Recurrent": "1st Recurrent",
        "1st Recurrent": "1st Recurrent",
        "2nd Recurrent": "2nd Recurrent",
    }
    result = progression_labels.get(value, value)
    return result


# Use the original stage substrings to classify the sample time point.
def derive_tissue_from_progression(value: str | None) -> object:
    """Return Primary, Recurrent, or pd.NA for unavailable or unknown stages."""
    progression = normalize_progression(value)
    if progression is pd.NA or pd.isna(progression):
        result = pd.NA
    elif progression == "Primary":
        result = "Primary"
    elif "Recurrent" in str(progression):
        result = "Recurrent"
    else:
        result = pd.NA

    return result


# Extract the first pair number from the original GEO pairing field.
def parse_pair_id(value: str | None) -> object:
    """Return the first matching integer, or pd.NA when the field has no match."""
    if value is None or pd.isna(value):
        result = pd.NA
        return result

    match = re.search(r"#?\s*(\d+)", str(value))
    result = int(match.group(1)) if match else pd.NA
    return result


# Read the ten original sample fields from the local cached GSM text.
def parse_geo_sample_text(text_path: Path) -> dict[str, object]:
    """Return the parsed GEO metadata, using the last value when a field is repeated."""
    metadata = {
        "gsm_id": text_path.stem,
        "sample_title": pd.NA,
        "diagnosis": pd.NA,
        "progression": pd.NA,
        "tissue": pd.NA,
        "pair_label": pd.NA,
        "pair_id": pd.NA,
        "age": pd.NA,
        "gender": pd.NA,
        "source_name": pd.NA,
    }
    for raw_line in text_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if raw_line.startswith("!Sample_title = "):
            metadata["sample_title"] = raw_line.split(" = ", 1)[1].strip()
            continue

        if raw_line.startswith("!Sample_source_name_ch1 = "):
            metadata["source_name"] = raw_line.split(" = ", 1)[1].strip()
            continue

        # Parse study characteristics using their original key and value fields.
        if not raw_line.startswith("!Sample_characteristics_ch1 = "):
            continue

        value = raw_line.split(" = ", 1)[1].strip()
        if ":" not in value:
            continue

        field, detail = value.split(":", 1)
        field = field.strip().lower()
        detail = detail.strip()
        if field == "progression":
            metadata["progression"] = normalize_progression(detail)
            metadata["tissue"] = derive_tissue_from_progression(detail)
        elif field == "diagnosis":
            metadata["diagnosis"] = detail
        elif field == "pair#":
            metadata["pair_label"] = detail
            metadata["pair_id"] = parse_pair_id(detail)
        elif field == "age":
            metadata["age"] = detail
        elif field == "gender":
            metadata["gender"] = detail

    return metadata


# Read standard Matrix Market counts and identifiers for one GSE274546 library.
def read_gse274546_library(library_id: str, raw_root: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Return cell-by-gene CSR counts, gene names, and barcodes after checking dimensions."""
    import numpy as np
    from scipy.io import mmread

    library_dir = raw_root / library_id
    matrix_path = library_dir / "matrix.mtx.gz"
    genes_path = library_dir / "genes.tsv.gz"
    barcodes_path = library_dir / "barcodes.tsv.gz"
    missing_paths = [path for path in (matrix_path, genes_path, barcodes_path) if not path.exists()]
    if missing_paths:
        raise FileNotFoundError("; ".join(map(str, missing_paths)))

    # Take the first column of each headerless TSV in the original file order.
    gene_table = pd.read_csv(genes_path, sep="\t", header=None)
    genes = gene_table.iloc[:, 0].astype(str).tolist()
    barcode_table = pd.read_csv(barcodes_path, sep="\t", header=None)
    barcodes = barcode_table.iloc[:, 0].astype(str).tolist()
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr().astype(np.int32)

    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(f"Matrix dimension mismatch for {library_id}: {matrix.shape}")

    result = matrix.T.tocsr(), make_unique(genes), barcodes
    return result


# Select the singlet cells in the order in which they appear in the classification table.
def select_singlets(
    library_id: str,
    matrix: sparse.csr_matrix,
    barcodes: list[str],
    classification_table: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    """Return the filtered counts and aligned annotations, checking barcode membership."""
    import numpy as np

    selected_table = classification_table.loc[classification_table["capture_id"].eq(library_id)].copy()
    barcode_positions = {barcode: index for index, barcode in enumerate(barcodes)}
    missing_barcodes = selected_table.loc[~selected_table["barcode"].isin(barcode_positions), "barcode"]
    if len(missing_barcodes):
        raise ValueError(f"{library_id} has {len(missing_barcodes)} classified barcodes absent from raw matrix")

    selected_indices = np.fromiter(
        (barcode_positions[barcode] for barcode in selected_table["barcode"]), dtype=np.int64
    )
    result = matrix[selected_indices].tocsr(), selected_table.reset_index(drop=True)
    return result


# Adjust finite P values within the supplied family and preserve missing entries.
def bh(values: pd.Series) -> np.ndarray:
    """Return the BH-adjusted P values aligned to the input series."""
    import numpy as np
    from statsmodels.stats.multitest import multipletests

    adjusted_values = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        adjusted_values[valid] = multipletests(values.loc[valid], method="fdr_bh")[1]

    return adjusted_values


# Compute a paired two-sided sign-flip P value from finite patient differences.
def sign_flip_p(delta: np.ndarray, seed: int = 20260712) -> float:
    """Use 200000 sign permutations and the original plus-one correction."""
    import numpy as np

    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    if len(delta) < 3:
        result = np.nan
        return result

    observed = abs(delta.mean())
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(200_000, len(delta)))
    null = np.abs((signs * delta).mean(axis=1))
    result = float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))
    return result


# Score detected genes in AnnData with the published Scanpy control-gene settings.
def control_gene_score(adata: ad.AnnData, genes: list[str], score_name: str) -> np.ndarray:
    """Return the per-cell scores and remove the temporary score column."""
    import numpy as np
    import scanpy as sc

    present_genes = [gene for gene in genes if gene in adata.var_names]
    if not present_genes:
        scores = np.full(adata.n_obs, np.nan)
        return scores

    # Reset the original seed so each call reproduces its published control genes.
    sc.tl.score_genes(
        adata,
        gene_list    = present_genes,
        score_name   = score_name,
        ctrl_size    = 50,
        n_bins       = 25,
        random_state = 20260712,
        use_raw      = False,
    )
    scores = adata.obs.pop(score_name).to_numpy()
    return scores


# Rank-transform the three indexed series and correlate their covariate-adjusted residuals.
def partial_spearman(x_values: pd.Series, y_values: pd.Series, covariate: pd.Series) -> tuple[float, float, int]:
    """Return the residual correlation, P value, and the number of complete rows."""
    import numpy as np
    from scipy import stats

    complete_table_total = pd.concat([x_values, y_values, covariate], axis=1).dropna()
    if len(complete_table_total) < 20:
        result = np.nan, np.nan, len(complete_table_total)
        return result

    ranks = complete_table_total.rank(method="average").to_numpy(float)
    design = np.column_stack([np.ones(len(ranks)), ranks[:, 2]])
    x_residual = ranks[:, 0] - design @ np.linalg.lstsq(design, ranks[:, 0], rcond=None)[0]
    y_residual = ranks[:, 1] - design @ np.linalg.lstsq(design, ranks[:, 1], rcond=None)[0]
    rho, p_value = stats.pearsonr(x_residual, y_residual)
    result = float(rho), float(p_value), len(complete_table_total)
    return result


# Read the local SOFT file and preserve its sample order and public download links.
def parse_soft(path: Path) -> list[dict[str, object]]:
    """Parse the sample records and return their selected metadata fields."""
    opener = gzip.open if path.suffix == ".gz" else open
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                samples.extend([current] if current is not None else [])
                current = {"geo_accession": line.split("=", 1)[1].strip(), "supplementary_urls": []}
                continue

            if current is None:
                continue

            # Parse colon-separated characteristics and convert FTP links to HTTPS.
            if line.startswith("!Sample_title = "):
                current["sample_title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = ") and ":" in line.split("=", 1)[1]:
                key, item = line.split("=", 1)[1].strip().split(":", 1)
                current[key.strip().lower().replace(" ", "_")] = item.strip()
            elif line.startswith("!Sample_supplementary_file"):
                url = line.split("=", 1)[1].strip().replace("ftp://", "https://")
                current["supplementary_urls"].append(url)

    if current is not None:
        samples.append(current)

    return samples


# Load raw counts and feature identifiers for one GSE174554 library.
def read_gse174554_library(row: pd.Series, raw_root: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Use the repaired parser for SF8963 and scipy.mmread for standard libraries."""
    import numpy as np
    from scipy.io import mmread

    matrix_path = raw_root / Path(row["matrix_path"]).name
    features_path = raw_root / Path(row["features_path"]).name
    barcodes_path = raw_root / Path(row["barcodes_path"]).name
    for path in (matrix_path, features_path, barcodes_path):
        if not path.exists():
            raise FileNotFoundError(path)

    # Read the first column from each headerless gene and barcode TSV.
    feature_table = pd.read_csv(features_path, sep="\t", header=None)
    features = feature_table.iloc[:, 0].astype(str).tolist()
    barcode_table = pd.read_csv(barcodes_path, sep="\t", header=None)
    barcodes = barcode_table.iloc[:, 0].astype(str).tolist()
    if str(row["library_id"]) == "GSM5319529_SF8963":
        matrix, _, _ = read_matrix_manual(matrix_path)
    else:
        with gzip.open(matrix_path, "rb") as handle:
            matrix = mmread(handle).tocsr().astype(np.int32)

    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(f"Matrix shape mismatch for {row['library_id']}: {matrix.shape}")

    result = matrix.T.tocsr(), features, barcodes
    return result


# Resolve repeated metadata within a group using the first available entry.
def first_nonmissing(values):
    """Return the first nonmissing value, or pd.NA for an entirely missing group."""
    result = next((value for value in values if pd.notna(value)), pd.NA)
    return result
