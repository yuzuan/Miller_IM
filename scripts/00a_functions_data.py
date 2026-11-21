#!/usr/bin/env python3
# ======================================================================
# 00a_functions_data.py —— Miller_IM 公共数据函数库
# Miller_IM：公开数据重建与样本清单。
#
# 目的：
# 集中多个脚本共用的数据读取、marker评分及pseudobulk汇总，保持原始算法和随机数设置。
# 本函数脚本不用单独浏览学习，在其它分析脚本调用到它时再学。
# 在编辑器调用处跳转到定义，可直接查看对应函数；函数查找依靠跳转，不靠文件头清单。
# 输入由各调用脚本提供；本模块导入时不读取数据，不运行分析，也不直接加载绘图库。

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


# 输入名称列表，返回同顺序唯一名称；重复项按原计数追加-1、-2，不改用其他去重规则。
def make_unique(names: list[str]) -> list[str]:
    """返回追加重复计数后缀的名称列表；输入和输出顺序一一对应。"""
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


# 输入异常MatrixMarket gzip路径，返回基因×细胞CSR及声明/解析条目数；保留原逐行正则和重复坐标求和规则。
def read_matrix_manual(path: Path) -> tuple[sparse.csr_matrix, int, int]:
    """返回(int32 CSR,声明条目数,解析条目数)；不填补缺失声明条目。"""
    import numpy as np
    from scipy import sparse

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                n_rows, n_columns, declared_nnz = map(int, line.split())
                break

        # 预分配长度来自文件表头；有效记录不足声明条目时，仅截取实际读到部分。
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


# 输入AnnData、保序marker字典、列名前缀和原调用种子，返回每marker集合一行的覆盖表。
def score_marker_sets(
    adata: ad.AnnData,
    marker_sets: dict[str, list[str]],
    prefix: str,
    random_state: int,
) -> pd.DataFrame:
    """原位写入score列并返回覆盖表；任一集合少于两个可测基因时停止。"""
    import scanpy as sc

    coverage_rows = []
    for name, genes in marker_sets.items():
        present_genes = [gene for gene in genes if gene in adata.var_names]
        if len(present_genes) < 2:
            raise ValueError(f"Marker set {name} has fewer than two genes")

        # 每集合沿用调用者原种子；移出循环或更换种子会改变对照基因抽样，不能作风格调整。
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


# 输入AnnData及聚类/结果字段名，返回每个cluster-gene一行的原Scanpy差异基因排序表。
def rank_clusters(adata: ad.AnnData, cluster_key: str, key: str) -> pd.DataFrame:
    """原位保存rank结果，返回包含cluster、gene、rank的表；保持t-test_overestim_var。"""
    import scanpy as sc

    sc.tl.rank_genes_groups(
        adata, groupby=cluster_key, method="t-test_overestim_var", use_raw=False, pts=True, key_added=key
    )
    marker_table = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    marker_table = marker_table.rename(columns={"group": "cluster", "names": "gene"})
    marker_table["cluster"] = marker_table["cluster"].astype(str)
    marker_table["rank"] = marker_table.groupby("cluster", observed=True).cumcount() + 1
    return marker_table


# 输入AnnData、分组列、文件前缀和输出目录；写计数矩阵及样本表，无返回值。
def write_pseudobulk(adata: ad.AnnData, group_columns: list[str], prefix: str, output_dir: Path) -> None:
    """按分组顺序汇总counts；输出基因×样本矩阵与每样本一行元数据，无返回值。"""
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

        # counts每行一个细胞，按axis=0汇总每个基因；PB编号、int64和压缩设置均沿用原实现。
        sample_id = f"PB{index + 1:04d}"
        vector = np.asarray(count_matrix[group_table["row"].to_numpy()].sum(axis=0)).ravel().astype(np.int64)
        vectors.append(pd.Series(vector, name=sample_id))
        metadata_row = dict(zip(group_columns, group_values))
        metadata_row.update(sample_id=sample_id, n_cells=len(group_table), library_size=int(vector.sum()))
        metadata_rows.append(metadata_row)

    # 一次拼接各样本向量，保留基因顺序；输出目录在写盘前创建。
    count_table_total = pd.concat(vectors, axis=1)
    count_table_total.insert(0, "gene", adata.var_names.astype(str))
    metadata_table = pd.DataFrame(metadata_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    count_path = output_dir / f"{prefix}_counts.csv.gz"
    count_table_total.to_csv(count_path, index=False, compression={"method": "gzip", "mtime": 0})
    metadata_path = output_dir / f"{prefix}_metadata.csv"
    metadata_table.to_csv(metadata_path, index=False)

    # 两个文件均为下游统计模型输入；核对非空和样本维度，不新增审核产物。
    if count_table_total.shape != (adata.n_vars, len(metadata_table) + 1):
        raise RuntimeError("Pseudobulk counts and metadata dimensions differ")

    for path in (count_path, metadata_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty pseudobulk output: {path}")


# 输入细胞×基因稀疏计数和每细胞总计数，返回按原float32顺序计算的log1p标准化矩阵。
def normalize_sparse(matrix: sparse.csr_matrix, totals: np.ndarray) -> sparse.csr_matrix:
    """按1e4总计数标准化并log1p；输入totals来自完整基因矩阵，不能用子矩阵重新计算。"""
    import numpy as np

    factors = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    normalized_matrix = matrix.astype(np.float32).multiply(factors[:, None]).tocsr()
    np.log1p(normalized_matrix.data, out=normalized_matrix.data)
    return normalized_matrix


# 输入gzip MatrixMarket路径，返回表头的行数、列数及非零条目数；找不到表头则停止。
def parse_matrix_header(matrix_path: Path) -> tuple[int, int, int]:
    """读取首个非%行，返回(n_rows,n_columns,nnz)，不读取或改写矩阵数值。"""
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue

            n_rows, n_columns, nnz = (int(value) for value in line.strip().split())
            dimensions = n_rows, n_columns, nnz
            return dimensions

    raise ValueError(f"Could not parse Matrix Market header: {matrix_path}")


# 输入GEO时期原值，返回清除首尾空格及统一Recurrent后的标签；缺失保持pd.NA。
def normalize_progression(value: str | None) -> object:
    """返回原公开时期标签映射；未识别的非缺失值不猜测、不改写。"""
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


# 输入原时期标签，返回Primary/Recurrent或pd.NA，保留原字符串包含判断。
def derive_tissue_from_progression(value: str | None) -> object:
    """由统一时期产生组织类别；未知或缺失输入返回pd.NA。"""
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


# 输入GEO配对说明，返回首个匹配整数；缺失或未匹配保持pd.NA。
def parse_pair_id(value: str | None) -> object:
    """按原正则提取pair编号，保留首个数字匹配规则。"""
    if value is None or pd.isna(value):
        result = pd.NA
        return result

    match = re.search(r"#?\s*(\d+)", str(value))
    result = int(match.group(1)) if match else pd.NA
    return result


# 输入缓存GSM文本路径，返回原十个字段的样本字典；只读取本地文件，不联网。
def parse_geo_sample_text(text_path: Path) -> dict[str, object]:
    """逐行提取GEO样本元数据；未出现字段保持pd.NA，重复字段仍由最后一行覆盖。"""
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

        # 只有characteristics行进入字段解析，提前跳过可避免五层嵌套。
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


# 输入GSE274546文库ID及原始目录，返回细胞×基因CSR、唯一基因名和barcode列表。
def read_gse274546_library(library_id: str, raw_root: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """读取标准MatrixMarket文库；缺文件或矩阵与无表头TSV维度不符时停止。"""
    import numpy as np
    from scipy.io import mmread

    library_dir = raw_root / library_id
    matrix_path = library_dir / "matrix.mtx.gz"
    genes_path = library_dir / "genes.tsv.gz"
    barcodes_path = library_dir / "barcodes.tsv.gz"
    missing_paths = [path for path in (matrix_path, genes_path, barcodes_path) if not path.exists()]
    if missing_paths:
        raise FileNotFoundError("; ".join(map(str, missing_paths)))

    # 无表头TSV使用首列；先保存整表再提取向量，保留原基因和barcode顺序。
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


# 输入文库ID、细胞×基因矩阵、barcode顺序和singlet分类表，返回筛选矩阵及对齐后的分类表。
def select_singlets(
    library_id: str,
    matrix: sparse.csr_matrix,
    barcodes: list[str],
    classification_table: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    """保持分类表行序筛选singlet；分类barcode不在原始矩阵时停止。"""
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


# 输入一列P值，返回同长度BH校正数组；缺失值位置保持NaN。
def bh(values: pd.Series) -> np.ndarray:
    """仅对非缺失P值执行fdr_bh，返回与输入位置一致的数组。"""
    import numpy as np
    from statsmodels.stats.multitest import multipletests

    adjusted_values = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        adjusted_values[valid] = multipletests(values.loc[valid], method="fdr_bh")[1]

    return adjusted_values


# 输入患者内差值及原种子，返回双侧符号翻转P值；少于3个有限差值时返回NaN。
def sign_flip_p(delta: np.ndarray, seed: int = 20260712) -> float:
    """沿用200000次随机符号翻转及加一校正，不改变抽样矩阵或计算顺序。"""
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


# 输入AnnData、基因列表及暂存score列名，返回每细胞一值的原Scanpy对照基因评分。
def control_gene_score(adata: ad.AnnData, genes: list[str], score_name: str) -> np.ndarray:
    """返回score数组并移除暂存列；ctrl_size=50、n_bins=25、seed=20260712均为原固定设置。"""
    import numpy as np
    import scanpy as sc

    present_genes = [gene for gene in genes if gene in adata.var_names]
    if not present_genes:
        scores = np.full(adata.n_obs, np.nan)
        return scores

    # 每次调用重设相同种子是原设计；改成共享随机流会改变对照基因选择。
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


# 输入三个同索引Series，返回调整协变量后的rank残差相关系数、P值和完整观测数。
def partial_spearman(x_values: pd.Series, y_values: pd.Series, covariate: pd.Series) -> tuple[float, float, int]:
    """至少20个完整观测；秩、最小二乘残差、Pearson顺序及公式与原实现相同。"""
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


# 输入本地SOFT路径，返回保序样本字典列表；保持原GEO字段解析和FTP转换规则。
def parse_soft(path: Path) -> list[dict[str, object]]:
    """解析 SAMPLE 区块；输入本地 SOFT，返回样本列表，不执行下载。"""
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

            # characteristics 的冒号分隔规则和 FTP 转 HTTPS 规则均保持原样。
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


# 输入GSE174554文库行和原始目录，返回细胞×基因CSR、基因与barcode列表；维度不符时停止。
def read_gse174554_library(row: pd.Series, raw_root: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """读取一个文库；SF8963沿用逐行解析，其他文库沿用scipy.mmread。"""
    import numpy as np
    from scipy.io import mmread

    matrix_path = raw_root / Path(row["matrix_path"]).name
    features_path = raw_root / Path(row["features_path"]).name
    barcodes_path = raw_root / Path(row["barcodes_path"]).name
    for path in (matrix_path, features_path, barcodes_path):
        if not path.exists():
            raise FileNotFoundError(path)

    # 两份TSV原本无表头，首列分别是基因及barcode；读取后再明确取列。
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


# 输入分组Series，返回首个非缺失元数据值；供多个样本汇总字段共用。
def first_nonmissing(values):
    """输入分组Series，返回首个非缺失值；全缺失时返回pd.NA。"""
    result = next((value for value in values if pd.notna(value)), pd.NA)
    return result
