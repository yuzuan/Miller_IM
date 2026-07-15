#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


MG_FULL = [
    "PDK4", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C",
    "GSTM3", "CH25H", "JUN", "SIGLEC8", "KLF6", "FOLR2", "AC253572.2", "NLRP3",
]
MG_LEADING6 = ["BHLHE41", "CH25H", "FOLR2", "JUN", "SGK1", "SIGLEC8"]
MYELOID = ["PTPRC", "TYROBP", "AIF1", "LST1", "FCER1G", "CTSS", "CSF1R"]
T4_MES_METABOLIC = [
    "CHI3L1", "CD44", "SERPINE1", "VEGFA", "ADM", "CA9", "BNIP3", "NDRG1",
    "SLC2A1", "HK2", "LDHA", "ENO1", "VIM", "LGALS3", "ANXA1", "TGFBI",
]
MYELOID_STATES = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]


def bh(values: pd.Series) -> np.ndarray:
    out = np.full(len(values), np.nan)
    ok = values.notna().to_numpy()
    if ok.any():
        out[ok] = multipletests(values.loc[ok], method="fdr_bh")[1]
    return out


def sign_flip_p(delta: np.ndarray, seed: int = 20260712) -> float:
    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    if len(delta) < 3:
        return np.nan
    observed = abs(delta.mean())
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(200_000, len(delta)))
    null = np.abs((signs * delta).mean(axis=1))
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def mean_expression(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    use = [gene for gene in genes if gene in adata.var_names]
    if not use:
        return np.full(adata.n_obs, np.nan)
    values = adata[:, use].X.mean(axis=1)
    return np.asarray(values).ravel()


def control_gene_score(adata: ad.AnnData, genes: list[str], score_name: str) -> np.ndarray:
    use = [gene for gene in genes if gene in adata.var_names]
    if not use:
        return np.full(adata.n_obs, np.nan)
    sc.tl.score_genes(
        adata,
        gene_list=use,
        score_name=score_name,
        ctrl_size=50,
        n_bins=25,
        random_state=20260712,
        use_raw=False,
    )
    return adata.obs.pop(score_name).to_numpy()


def read_markers(path: Path, top_n: int = 40) -> dict[str, list[str]]:
    workbook = pd.ExcelFile(path)
    excluded = set(MG_FULL)
    markers = {}
    for sheet in workbook.sheet_names:
        table = pd.read_excel(path, sheet_name=sheet)
        if "FDR" not in table.columns:
            table = pd.read_excel(path, sheet_name=sheet, header=1)
        gene_col = table.columns[0]
        chosen = table.loc[
            table["FDR"].lt(0.05) & table["Foldchange"].gt(0), gene_col
        ].astype(str)
        markers[sheet] = [gene for gene in chosen if gene not in excluded][:top_n]
    return markers


def read_single_cell(source: Path, manifest_path: Path) -> ad.AnnData:
    manifest = pd.read_csv(manifest_path).set_index("geo_accession")
    objects = []
    qc_rows = []
    for path in sorted(source.glob("*.h5")):
        accession = re.match(r"(GSM\d+)", path.name).group(1)
        sample = manifest.loc[accession, "sample_title"]
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        adata.obs_names = [f"{sample}:{barcode}" for barcode in adata.obs_names]
        adata.obs["sample"] = sample
        adata.obs["geo_accession"] = accession
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].ge(200)
            & adata.obs["n_genes_by_counts"].le(8_000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        qc_rows.append({
            "sample": sample,
            "geo_accession": accession,
            "n_input": adata.n_obs,
            "n_qc_pass": int(keep.sum()),
            "fraction_qc_pass": float(keep.mean()),
        })
        objects.append(adata[keep].copy())
    combined = ad.concat(objects, join="inner", merge="same")
    combined.uns["qc_summary"] = pd.DataFrame(qc_rows).to_dict("list")
    return combined


def annotate_and_test(adata: ad.AnnData, markers: dict[str, list[str]], write: Path, figures: Path) -> None:
    qc = pd.DataFrame(adata.uns.pop("qc_summary"))
    qc.to_csv(write / "GSE278456_single_cell_qc.csv", index=False)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    for state, genes in markers.items():
        adata.obs[f"author_{state}"] = control_gene_score(adata, genes, f"tmp_author_{state}")
    score_cols = [f"author_{state}" for state in markers]
    score_matrix = adata.obs[score_cols].to_numpy()
    state_names = np.array(list(markers))
    adata.obs["author_like_state"] = state_names[np.nanargmax(score_matrix, axis=1)]
    adata.obs["Mg_full"] = control_gene_score(adata, MG_FULL, "tmp_Mg_full")
    adata.obs["Mg_leading6"] = control_gene_score(adata, MG_LEADING6, "tmp_Mg_leading6")
    adata.obs["myeloid_abundance"] = control_gene_score(adata, MYELOID, "tmp_myeloid")

    coverage = []
    for name, genes in {"Mg_full": MG_FULL, "Mg_leading6": MG_LEADING6, "myeloid_abundance": MYELOID}.items():
        use = [gene for gene in genes if gene in adata.var_names]
        coverage.append({"program": name, "n_defined": len(genes), "n_detected": len(use), "detected_genes": ";".join(use)})
    pd.DataFrame(coverage).to_csv(write / "GSE278456_program_gene_coverage.csv", index=False)

    counts = adata.obs.groupby(["sample", "author_like_state"], observed=True).size().rename("n_cells").reset_index()
    totals = adata.obs.groupby("sample", observed=True).size().rename("n_total").reset_index()
    counts = counts.merge(totals, on="sample")
    counts["fraction"] = counts["n_cells"] / counts["n_total"]
    counts.to_csv(write / "GSE278456_author_like_state_counts.csv", index=False)

    long = adata.obs.loc[adata.obs["author_like_state"].isin(MYELOID_STATES), ["sample", "author_like_state", "Mg_full", "Mg_leading6"]]
    patient_state = long.groupby(["sample", "author_like_state"], observed=True).agg(
        n_cells=("Mg_full", "size"), Mg_full=("Mg_full", "mean"), Mg_leading6=("Mg_leading6", "mean")
    ).reset_index()
    patient_state.to_csv(write / "GSE278456_patient_state_program_scores.csv", index=False)

    tests = []
    for program in ["Mg_full", "Mg_leading6"]:
        for state in MYELOID_STATES:
            deltas = []
            for sample, frame in long.groupby("sample", observed=True):
                target = frame.loc[frame["author_like_state"].eq(state), program]
                other = frame.loc[~frame["author_like_state"].eq(state), program]
                if len(target) >= 30 and len(other) >= 30:
                    deltas.append((sample, target.mean() - other.mean()))
            values = np.array([value for _, value in deltas])
            tests.append({
                "program": program,
                "state": state,
                "n_patients": len(values),
                "mean_within_patient_delta": float(values.mean()) if len(values) else np.nan,
                "median_within_patient_delta": float(np.median(values)) if len(values) else np.nan,
                "n_positive": int((values > 0).sum()),
                "sign_flip_p": sign_flip_p(values),
                "wilcoxon_p": stats.wilcoxon(values).pvalue if len(values) >= 5 and np.any(values) else np.nan,
            })
    tests = pd.DataFrame(tests)
    tests["fdr"] = tests.groupby("program")["sign_flip_p"].transform(lambda x: bh(x))
    tests.to_csv(write / "GSE278456_patient_level_state_tests.csv", index=False)

    state_order = patient_state.groupby("author_like_state")["Mg_full"].mean().sort_values(ascending=False).index
    plt.figure(figsize=(9, 5.2))
    sns.boxplot(data=patient_state, x="author_like_state", y="Mg_full", order=state_order, color="#d9e4dc", fliersize=0)
    sns.stripplot(data=patient_state, x="author_like_state", y="Mg_full", order=state_order, color="#1f4e5f", size=3, alpha=0.65)
    plt.xlabel("")
    plt.ylabel("Mg-inflammatory program score")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.savefig(figures / "GSE278456_patient_level_Mg_program_by_state.pdf")
    plt.savefig(figures / "GSE278456_patient_level_Mg_program_by_state.png", dpi=220)
    plt.close()

    obs_columns = ["sample", "geo_accession", "author_like_state", "Mg_full", "Mg_leading6", "myeloid_abundance"]
    adata.obs[obs_columns].to_csv(
        write / "GSE278456_cell_annotations_scores.csv.gz",
        compression={"method": "gzip", "mtime": 0},
    )


def partial_spearman(x: pd.Series, y: pd.Series, covariate: pd.Series) -> tuple[float, float, int]:
    table = pd.concat([x, y, covariate], axis=1).dropna()
    if len(table) < 20:
        return np.nan, np.nan, len(table)
    ranks = table.rank(method="average").to_numpy(float)
    design = np.column_stack([np.ones(len(ranks)), ranks[:, 2]])
    rx = ranks[:, 0] - design @ np.linalg.lstsq(design, ranks[:, 0], rcond=None)[0]
    ry = ranks[:, 1] - design @ np.linalg.lstsq(design, ranks[:, 1], rcond=None)[0]
    rho, p = stats.pearsonr(rx, ry)
    return float(rho), float(p), len(table)


def spatial_validation(source: Path, manifest_path: Path, markers: dict[str, list[str]], write: Path, figures: Path) -> None:
    manifest = pd.read_csv(manifest_path).set_index("geo_accession")
    rows = []
    spot_tables = []
    map_frames = []
    for matrix_path in sorted(source.glob("*_filtered_feature_bc_matrix.h5")):
        accession = re.match(r"(GSM\d+)", matrix_path.name).group(1)
        sample = manifest.loc[accession, "sample_title"]
        prefix = matrix_path.name.replace("_filtered_feature_bc_matrix.h5", "")
        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=10_000)
        sc.pp.log1p(adata)
        scores = pd.DataFrame(index=adata.obs_names)
        scores["Mg_full"] = control_gene_score(adata, MG_FULL, "tmp_Mg_full")
        scores["Mg_leading6"] = control_gene_score(adata, MG_LEADING6, "tmp_Mg_leading6")
        scores["myeloid_abundance"] = control_gene_score(adata, MYELOID, "tmp_myeloid")
        scores["E_MDSC"] = control_gene_score(adata, markers["E-MDSC"], "tmp_E_MDSC")
        scores["M_MDSC"] = control_gene_score(adata, markers["M-MDSC"], "tmp_M_MDSC")
        scores["T4_MES_metabolic"] = control_gene_score(adata, T4_MES_METABOLIC, "tmp_T4_MES")
        scores["sample"] = sample
        scores["geo_accession"] = accession

        positions_path = source / f"{prefix}_tissue_positions.csv.gz"
        positions = pd.read_csv(positions_path).set_index("barcode")
        scores = scores.join(positions, how="left")
        scores = scores.loc[scores["in_tissue"].eq(1)].copy()
        spot_tables.append(scores.reset_index(names="barcode"))

        for target in ["E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            rho, p, n = partial_spearman(scores["Mg_full"], scores[target], scores["myeloid_abundance"])
            rows.append({"sample": sample, "geo_accession": accession, "target": target, "n_spots": n, "partial_spearman_rho": rho, "p_value": p})

        for program in ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            frame = scores[["pxl_col_in_fullres", "pxl_row_in_fullres", program]].copy()
            frame["program"] = program
            frame["value"] = frame.pop(program)
            frame["sample"] = sample
            map_frames.append(frame)

    tests = pd.DataFrame(rows)
    tests["fdr"] = bh(tests["p_value"])
    tests.to_csv(write / "GSE276841_spatial_partial_correlations.csv", index=False)
    pd.concat(spot_tables, ignore_index=True).to_csv(
        write / "GSE276841_spot_scores.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    maps = pd.concat(map_frames, ignore_index=True)
    samples = maps["sample"].drop_duplicates().tolist()
    programs = ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]
    fig, axes = plt.subplots(len(samples), len(programs), figsize=(13, 6.4), squeeze=False)
    for row, sample in enumerate(samples):
        for col, program in enumerate(programs):
            ax = axes[row, col]
            frame = maps.loc[maps["sample"].eq(sample) & maps["program"].eq(program)]
            low, high = frame["value"].quantile([0.02, 0.98])
            points = ax.scatter(
                frame["pxl_col_in_fullres"], -frame["pxl_row_in_fullres"], c=frame["value"],
                s=4, cmap="magma", vmin=low, vmax=high, linewidths=0,
            )
            ax.set_title(f"{sample} | {program}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02)
    plt.tight_layout()
    plt.savefig(figures / "GSE276841_spatial_program_maps.pdf")
    plt.savefig(figures / "GSE276841_spatial_program_maps.png", dpi=220)
    plt.close()


def write_summary(write: Path) -> None:
    tests = pd.read_csv(write / "GSE278456_patient_level_state_tests.csv")
    spatial = pd.read_csv(write / "GSE276841_spatial_partial_correlations.csv")
    top = tests.loc[tests["program"].eq("Mg_full")].sort_values("fdr").head(5)
    lines = [
        "# GSE278456 / GSE276841 independent validation",
        "",
        "## Single-cell patient-level state mapping",
        "",
        top.to_markdown(index=False),
        "",
        "## Spatial partial correlations adjusted for total myeloid signal",
        "",
        spatial.to_markdown(index=False),
        "",
        "Interpretation must remain claim-specific: this cohort is untreated primary IDH-wildtype GBM, so it tests phenotype and spatial ecology, not recurrence direction.",
    ]
    (write / "GSE278456_GSE276841_validation_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
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

    markers = read_markers(source / "NIHMS2115262-supplement-S2_S5_S6.xls")
    marker_rows = [{"state": state, "n_genes": len(genes), "genes": ";".join(genes)} for state, genes in markers.items()]
    pd.DataFrame(marker_rows).to_csv(write / "GSE278456_author_marker_sets_used.csv", index=False)

    adata = read_single_cell(source / "GSE278456_myeloid_h5", write / "GSE278456_sample_manifest.csv")
    annotate_and_test(adata, markers, write, figures)
    spatial_validation(source / "GSE276841_spatial", write / "GSE276841_sample_manifest.csv", markers, write, figures)
    write_summary(write)


if __name__ == "__main__":
    main()
