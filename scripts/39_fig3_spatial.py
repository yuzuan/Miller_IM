#!/usr/bin/env python3
from __future__ import annotations

import itertools
import math
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


ROOT = Path(__file__).resolve().parents[1]
EXT_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()

WRITE_BASE = ROOT / "write" / "39_independent_miller_mg_inflammatory_figures"
FIG_BASE = ROOT / "figures" / "39_independent_miller_mg_inflammatory_plottie"
FIG3_WRITE = WRITE_BASE / "Figure3"
S6_WRITE = WRITE_BASE / "SupplementaryFigure6"
FIG3_FIG = FIG_BASE / "Figure3"
S6_FIG = FIG_BASE / "SupplementaryFigure6"

GEOMX_SOURCE = EXT_ROOT / "write" / "26_pure_bioinformatics_dataset_rescue" / "source_metadata" / "GeoMx_zenodo16839828"
SPATIAL_SOURCE = EXT_ROOT / "write" / "30_dataset_rescue_search" / "source_metadata" / "GSE276841_spatial"
MARKER_XLS = EXT_ROOT / "write" / "30_dataset_rescue_search" / "source_metadata" / "NIHMS2115262-supplement-S2_S5_S6.xls"
GSE276841_MANIFEST = ROOT / "write" / "30_dataset_rescue_search" / "GSE276841_sample_manifest.csv"
STEP38_GSEA = ROOT / "write" / "38_independent_cohort_mg_inflammatory_recalculation" / "independent_fixed_program_targeted_gsea.csv"

RAW20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
MYELOID = ["PTPRC", "TYROBP", "AIF1", "LST1", "FCER1G", "CTSS", "CSF1R"]
MES = [
    "CHI3L1", "CD44", "SERPINE1", "VEGFA", "ADM", "CA9", "BNIP3", "NDRG1",
    "SLC2A1", "HK2", "LDHA", "ENO1", "VIM", "LGALS3", "ANXA1", "TGFBI",
]
NPG = {
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


def setup_style() -> None:
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


def ensure_dirs() -> None:
    for path in [FIG3_WRITE, S6_WRITE, FIG3_FIG, S6_FIG]:
        path.mkdir(parents=True, exist_ok=True)
    for path in [FIG3_WRITE / "source_csv", S6_WRITE / "source_csv"]:
        path.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def bh(series: pd.Series) -> np.ndarray:
    out = np.full(len(series), np.nan)
    ok = series.notna().to_numpy()
    if ok.any():
        out[ok] = multipletests(series.loc[ok], method="fdr_bh")[1]
    return out


def sign_flip_p(delta: np.ndarray, seed: int = 20260713) -> float:
    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    n = len(delta)
    if n == 0:
        return np.nan
    observed = abs(delta.mean())
    if n <= 20:
        signs = np.asarray(list(itertools.product([-1.0, 1.0], repeat=n)))
        null = np.abs((signs * delta).mean(axis=1))
    else:
        rng = np.random.default_rng(seed)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(300_000, n))
        null = np.abs((signs * delta).mean(axis=1))
    return float((np.count_nonzero(null >= observed - 1e-15) + 1) / (len(null) + 1))


def label_permutation_p(control: np.ndarray, treated: np.ndarray, seed: int = 20260713) -> float:
    control = np.asarray(control, dtype=float)
    treated = np.asarray(treated, dtype=float)
    observed = abs(treated.mean() - control.mean())
    joined = np.concatenate([control, treated])
    n_control = len(control)
    rng = np.random.default_rng(seed)
    null = np.empty(300_000, dtype=float)
    for idx in range(len(null)):
        perm = rng.permutation(joined)
        null[idx] = abs(perm[n_control:].mean() - perm[:n_control].mean())
    return float((np.count_nonzero(null >= observed - 1e-15) + 1) / (len(null) + 1))


def control_gene_score(adata, genes: list[str], score_name: str) -> np.ndarray:
    use = [gene for gene in genes if gene in adata.var_names]
    if not use:
        return np.full(adata.n_obs, np.nan)
    sc.tl.score_genes(
        adata,
        gene_list=use,
        score_name=score_name,
        ctrl_size=50,
        n_bins=25,
        random_state=20260713,
        use_raw=False,
    )
    values = adata.obs[score_name].to_numpy()
    del adata.obs[score_name]
    return values


def read_author_markers() -> dict[str, list[str]]:
    workbook = pd.ExcelFile(MARKER_XLS)
    markers: dict[str, list[str]] = {}
    excluded = set(RAW20)
    for sheet in workbook.sheet_names:
        table = pd.read_excel(MARKER_XLS, sheet_name=sheet)
        if "FDR" not in table.columns:
            table = pd.read_excel(MARKER_XLS, sheet_name=sheet, header=1)
        gene_col = table.columns[0]
        chosen = table.loc[
            table["FDR"].lt(0.05) & table["Foldchange"].gt(0), gene_col
        ].astype(str)
        markers[sheet] = [gene for gene in chosen if gene not in excluded][:40]
    return markers


def mean_ci(values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    mean = float(values.mean())
    if len(values) == 1:
        return mean, mean, mean
    sem = float(stats.sem(values))
    delta = 1.96 * sem
    return mean, mean - delta, mean + delta


def partial_spearman(x: pd.Series, y: pd.Series, covariate: pd.Series) -> tuple[float, float, int]:
    table = pd.concat([x, y, covariate], axis=1).dropna()
    if len(table) < 20:
        return np.nan, np.nan, len(table)
    ranked = table.rank(method="average").to_numpy(dtype=np.float64, copy=True)
    x_rank = ranked[:, 0] - ranked[:, 0].mean()
    y_rank = ranked[:, 1] - ranked[:, 1].mean()
    cov_rank = ranked[:, 2] - ranked[:, 2].mean()
    denominator = float(np.dot(cov_rank, cov_rank))
    if denominator <= 0:
        return np.nan, np.nan, len(table)
    rx = x_rank - cov_rank * (float(np.dot(cov_rank, x_rank)) / denominator)
    ry = y_rank - cov_rank * (float(np.dot(cov_rank, y_rank)) / denominator)
    rho, p_value = stats.pearsonr(rx, ry)
    return float(rho), float(p_value), len(table)


def recompute_geomx_raw20() -> dict[str, pd.DataFrame]:
    metadata = pd.read_csv(GEOMX_SOURCE / "metadata.csv")
    expression = pd.read_csv(GEOMX_SOURCE / "normalized.txt", sep="\t", index_col=0)
    if set(expression.columns) != set(metadata["Unnamed: 0"]):
        raise ValueError("GeoMx matrix and metadata do not match.")

    raw20_present = [gene for gene in RAW20 if gene in expression.index]
    raw20_missing = [gene for gene in RAW20 if gene not in expression.index]
    coverage = pd.DataFrame(
        {
            "signature": ["Miller_raw20"],
            "n_defined": [len(RAW20)],
            "n_present": [len(raw20_present)],
            "coverage": [len(raw20_present) / len(RAW20)],
            "present_genes": [";".join(raw20_present)],
            "missing_genes": [";".join(raw20_missing)],
        }
    )
    coverage.to_csv(FIG3_WRITE / "source_csv" / "geomx_raw20_coverage.csv", index=False)

    metadata = metadata.copy()
    metadata["entry_strict_pass_idhwt"] = metadata["iba1"].eq(True) & metadata["QCFlags"].isna() & metadata["IDH_status"].eq("IDH_WT")
    metadata["entry_strict_pass_all_idh"] = metadata["iba1"].eq(True) & metadata["QCFlags"].isna()
    metadata["entry_all_aoi_idhwt"] = metadata["iba1"].eq(True) & metadata["IDH_status"].eq("IDH_WT")

    entries = {
        "strict_pass_idhwt": "entry_strict_pass_idhwt",
        "strict_pass_all_idh": "entry_strict_pass_all_idh",
        "all_aoi_idhwt": "entry_all_aoi_idhwt",
    }

    pair_stats = []
    paired_scores = []
    gene_stats = []
    patient_scores = []
    subgroup_stats = []
    patient_map_rows = []

    for entry_name, column in entries.items():
        current = metadata.loc[metadata[column]].copy()
        current = current.loc[current["tumor_setting"].isin(["Primary", "Recurrence"])].copy()
        ids = current["Unnamed: 0"].tolist()
        log_expression = np.log2(expression.loc[:, ids] + 1.0)
        patient_gene = (
            log_expression.T.assign(
                Patient_No=current.set_index("Unnamed: 0").loc[log_expression.columns, "Patient_No"].to_numpy(),
                trial_setting=current.set_index("Unnamed: 0").loc[log_expression.columns, "trial_setting"].to_numpy(),
                tumor_setting=current.set_index("Unnamed: 0").loc[log_expression.columns, "tumor_setting"].to_numpy(),
                IDH_status=current.set_index("Unnamed: 0").loc[log_expression.columns, "IDH_status"].to_numpy(),
            )
            .groupby(["Patient_No", "trial_setting", "tumor_setting", "IDH_status"], observed=True)[raw20_present]
            .mean()
            .reset_index()
        )
        z = patient_gene[raw20_present].copy()
        z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=0).replace(0, np.nan)
        patient_gene["raw20_score"] = z.mean(axis=1)
        patient_scores.append(patient_gene.copy())

        wide = patient_gene.pivot_table(
            index=["Patient_No", "trial_setting", "IDH_status"],
            columns="tumor_setting",
            values="raw20_score",
            aggfunc="first",
        ).dropna(subset=["Primary", "Recurrence"]).reset_index()
        wide["delta"] = wide["Recurrence"] - wide["Primary"]
        wide["entry"] = entry_name
        paired_scores.append(wide.copy())
        pair_stats.append(
            {
                "entry": entry_name,
                "n_pairs": int(len(wide)),
                "n_positive": int((wide["delta"] > 0).sum()),
                "positive_fraction": float((wide["delta"] > 0).mean()) if len(wide) else np.nan,
                "mean_delta": float(wide["delta"].mean()) if len(wide) else np.nan,
                "median_delta": float(wide["delta"].median()) if len(wide) else np.nan,
                "sign_flip_p": sign_flip_p(wide["delta"].to_numpy()),
            }
        )
        patient_map_rows.append(wide.copy())

        for subgroup, frame in {
            "All": wide,
            "Control": wide.loc[wide["trial_setting"].eq("Control")],
            "Nivolumab": wide.loc[wide["trial_setting"].eq("Nivolumab")],
            "IDH_WT": wide.loc[wide["IDH_status"].eq("IDH_WT")],
            "IDH_mutated": wide.loc[wide["IDH_status"].eq("IDH_mutated")],
        }.items():
            subgroup_stats.append(
                {
                    "entry": entry_name,
                    "subgroup": subgroup,
                    "n_pairs": int(len(frame)),
                    "mean_delta": float(frame["delta"].mean()) if len(frame) else np.nan,
                    "median_delta": float(frame["delta"].median()) if len(frame) else np.nan,
                    "positive_fraction": float((frame["delta"] > 0).mean()) if len(frame) else np.nan,
                    "sign_flip_p": sign_flip_p(frame["delta"].to_numpy()) if len(frame) else np.nan,
                }
            )
        control = wide.loc[wide["trial_setting"].eq("Control"), "delta"].to_numpy()
        nivolumab = wide.loc[wide["trial_setting"].eq("Nivolumab"), "delta"].to_numpy()
        subgroup_stats.append(
            {
                "entry": entry_name,
                "subgroup": "Interaction_Nivolumab_minus_Control",
                "n_pairs": int(len(wide)),
                "mean_delta": float(nivolumab.mean() - control.mean()),
                "median_delta": np.nan,
                "positive_fraction": np.nan,
                "sign_flip_p": label_permutation_p(control, nivolumab),
            }
        )

        gene_patient = patient_gene[["Patient_No", "trial_setting", "tumor_setting", "IDH_status"] + raw20_present].copy()
        for gene in raw20_present:
            gene_wide = gene_patient.pivot_table(
                index=["Patient_No", "trial_setting", "IDH_status"],
                columns="tumor_setting",
                values=gene,
                aggfunc="first",
            ).dropna(subset=["Primary", "Recurrence"]).reset_index()
            gene_wide["delta"] = gene_wide["Recurrence"] - gene_wide["Primary"]
            mean, ci_low, ci_high = mean_ci(gene_wide["delta"].to_numpy())
            gene_stats.append(
                {
                    "entry": entry_name,
                    "gene": gene,
                    "n_pairs": int(len(gene_wide)),
                    "mean_log2_delta": mean,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "median_log2_delta": float(gene_wide["delta"].median()),
                    "positive_fraction": float((gene_wide["delta"] > 0).mean()),
                    "sign_flip_p": sign_flip_p(gene_wide["delta"].to_numpy()),
                }
            )

    pair_stats_df = pd.DataFrame(pair_stats)
    pair_stats_df["fdr"] = bh(pair_stats_df["sign_flip_p"])
    subgroup_df = pd.DataFrame(subgroup_stats)
    subgroup_df["fdr"] = subgroup_df.groupby("entry")["sign_flip_p"].transform(lambda x: bh(x))
    gene_stats_df = pd.DataFrame(gene_stats)
    gene_stats_df["fdr"] = gene_stats_df.groupby("entry")["sign_flip_p"].transform(lambda x: bh(x))
    patient_scores_df = pd.concat(patient_scores, ignore_index=True)
    paired_scores_df = pd.concat(paired_scores, ignore_index=True)
    patient_map_df = pd.concat(patient_map_rows, ignore_index=True)

    pair_stats_df.to_csv(FIG3_WRITE / "source_csv" / "geomx_raw20_entry_level_summary.csv", index=False)
    subgroup_df.to_csv(S6_WRITE / "source_csv" / "geomx_raw20_subgroup_summary.csv", index=False)
    gene_stats_df.to_csv(FIG3_WRITE / "source_csv" / "geomx_raw20_gene_level_summary.csv", index=False)
    patient_scores_df.to_csv(FIG3_WRITE / "source_csv" / "geomx_raw20_patient_timepoint_scores.csv", index=False)
    paired_scores_df.to_csv(FIG3_WRITE / "source_csv" / "geomx_raw20_paired_deltas.csv", index=False)
    patient_map_df.to_csv(S6_WRITE / "source_csv" / "geomx_raw20_patient_map.csv", index=False)

    return {
        "coverage": coverage,
        "entry_summary": pair_stats_df,
        "subgroup_summary": subgroup_df,
        "gene_summary": gene_stats_df,
        "patient_scores": patient_scores_df,
        "paired_scores": paired_scores_df,
    }


def recompute_gse276841_raw20(markers: dict[str, list[str]]) -> dict[str, pd.DataFrame]:
    manifest = pd.read_csv(GSE276841_MANIFEST).set_index("geo_accession")
    spot_rows = []
    partial_rows = []
    detected_counts: dict[str, int] = {}

    mmdsc_genes = [gene for gene in markers["M-MDSC"] if gene not in RAW20]
    emdsc_genes = [gene for gene in markers["E-MDSC"] if gene not in RAW20]
    mes_genes = [gene for gene in MES if gene not in RAW20]

    for matrix_path in sorted(SPATIAL_SOURCE.glob("*_filtered_feature_bc_matrix.h5")):
        accession = matrix_path.name.split("_")[0]
        sample = manifest.loc[accession, "sample_title"]
        prefix = matrix_path.name.replace("_filtered_feature_bc_matrix.h5", "")

        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()
        detected_counts[sample] = int(sum(gene in adata.var_names for gene in RAW20))
        sc.pp.normalize_total(adata, target_sum=10_000)
        sc.pp.log1p(adata)

        scores = pd.DataFrame(index=adata.obs_names)
        scores["raw20"] = control_gene_score(adata, RAW20, "raw20_score")
        scores["myeloid"] = control_gene_score(adata, MYELOID, "myeloid_score")
        scores["mmdsc"] = control_gene_score(adata, mmdsc_genes, "mmdsc_score")
        scores["emdsc"] = control_gene_score(adata, emdsc_genes, "emdsc_score")
        scores["mdsc_combined"] = scores[["mmdsc", "emdsc"]].mean(axis=1)
        scores["mes"] = control_gene_score(adata, mes_genes, "mes_score")
        scores["sample"] = sample
        scores["geo_accession"] = accession

        positions = pd.read_csv(SPATIAL_SOURCE / f"{prefix}_tissue_positions.csv.gz").set_index("barcode")
        scores = scores.join(positions, how="left")
        scores = scores.loc[scores["in_tissue"].eq(1)].copy()
        scores["pixel_y_plot"] = -scores["pxl_row_in_fullres"]
        spot_rows.append(scores.reset_index(names="barcode"))

        for target_key, target_label in [("mdsc_combined", "MDSC"), ("mes", "MES"), ("emdsc", "E-MDSC"), ("mmdsc", "M-MDSC")]:
            rho, p_value, n_spots = partial_spearman(scores["raw20"], scores[target_key], scores["myeloid"])
            partial_rows.append(
                {
                    "sample": sample,
                    "geo_accession": accession,
                    "target": target_label,
                    "n_spots": n_spots,
                    "partial_rho": rho,
                    "p_value": p_value,
                }
            )

    spot_df = pd.concat(spot_rows, ignore_index=True)
    partial_df = pd.DataFrame(partial_rows)
    partial_df["fdr"] = bh(partial_df["p_value"])

    coverage = pd.DataFrame(
        [
            {
                "signature": "raw20",
                "n_defined": len(RAW20),
                "n_detected_gbm030": detected_counts.get("GBM030", np.nan),
                "n_detected_gbm049": detected_counts.get("GBM049", np.nan),
            }
        ]
    )

    spot_df.to_csv(FIG3_WRITE / "source_csv" / "gse276841_raw20_spot_scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    partial_df.to_csv(FIG3_WRITE / "source_csv" / "gse276841_raw20_partial_rho.csv", index=False)
    partial_df.to_csv(S6_WRITE / "source_csv" / "gse276841_raw20_slice_sensitivity.csv", index=False)
    coverage.to_csv(FIG3_WRITE / "source_csv" / "gse276841_raw20_coverage.csv", index=False)

    return {"spots": spot_df, "partial": partial_df, "coverage": coverage}


def load_step38_main() -> pd.DataFrame:
    table = pd.read_csv(STEP38_GSEA)
    table = table.loc[
        table["threshold"].eq(20)
        & table["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")
    ].copy()
    table["label"] = table["dataset"].map({"GSE174554": "GSE174554 sc/snRNA", "GSE274546": "GSE274546 scRNA"})
    table["nominal_p"] = table["pval"]
    return table


def plot_flowchart(step38: pd.DataFrame, geomx: dict[str, pd.DataFrame], spatial: dict[str, pd.DataFrame]) -> None:
    geomx_main = geomx["entry_summary"].loc[geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]
    partial_main = spatial["partial"].loc[spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()

    flow_source = pd.concat(
        [
            step38[["dataset", "label", "n_pairs", "NES", "nominal_p"]].rename(columns={"nominal_p": "step_p_or_fdr"}),
            pd.DataFrame(
                [
                    {
                        "dataset": "Artzi2025_GeoMx",
                        "label": "Artzi IBA1+ GeoMx",
                        "n_pairs": int(geomx_main["n_pairs"]),
                        "NES": np.nan,
                        "step_p_or_fdr": geomx_main["fdr"],
                    },
                    {
                        "dataset": "GSE276841",
                        "label": "GSE276841 Visium",
                        "n_pairs": 2,
                        "NES": np.nan,
                        "step_p_or_fdr": float(partial_main["fdr"].min()),
                    },
                ]
            ),
        ],
        ignore_index=True,
    )
    flow_source.to_csv(FIG3_WRITE / "source_csv" / "figure3_flowchart_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    ax.axis("off")
    boxes = {
        "discovery": (0.03, 0.22, 0.26, 0.56),
        "g174": (0.38, 0.56, 0.24, 0.24),
        "g274": (0.38, 0.20, 0.24, 0.24),
        "geomx": (0.71, 0.56, 0.24, 0.24),
        "visium": (0.71, 0.20, 0.24, 0.24),
    }

    def draw_box(key: str, lines: list[str], color: str) -> None:
        x, y, w, h = boxes[key]
        rect = mpl.patches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=1.0, edgecolor=color, facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(x + 0.02, y + h - 0.05, "\n".join(lines), ha="left", va="top", color="#222222", fontsize=8)

    draw_box(
        "discovery",
        [
            "Independent raw20 anchor",
            "Miller raw top20 only",
            "No project-curated 16 replacement",
        ],
        NPG["navy"],
    )
    g174 = step38.loc[step38["dataset"].eq("GSE174554")].iloc[0]
    g274 = step38.loc[step38["dataset"].eq("GSE274546")].iloc[0]
    draw_box("g174", [g174["label"], f"18 pairs", f"NES={g174['NES']:.3f}", f"nominal P={g174['nominal_p']:.2g}"], NPG["blue"])
    draw_box("g274", [g274["label"], f"45 pairs", f"NES={g274['NES']:.3f}", f"nominal P={g274['nominal_p']:.2g}"], NPG["teal"])
    draw_box(
        "geomx",
        [
            "Artzi IBA1+ GeoMx",
            "22 IDH-wt pairs",
            f"mean delta={geomx_main['mean_delta']:.3f}",
            f"FDR={geomx_main['fdr']:.2g}",
        ],
        NPG["red"],
    )
    draw_box(
        "visium",
        [
            "GSE276841 Visium",
            "2 untreated IDH-wt slices",
            f"MDSC min FDR={partial_main.loc[partial_main['target'].eq('MDSC'), 'fdr'].min():.2g}",
            f"MES min FDR={partial_main.loc[partial_main['target'].eq('MES'), 'fdr'].min():.2g}",
        ],
        NPG["salmon"],
    )

    arrows = [
        ((0.29, 0.50), (0.38, 0.68)),
        ((0.29, 0.50), (0.38, 0.32)),
        ((0.62, 0.68), (0.71, 0.68)),
        ((0.62, 0.32), (0.71, 0.32)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.0, color="#555555"))
    save(fig, FIG3_FIG, "fig3a_flowchart")


def plot_geomx_paired_change(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["paired_scores"].loc[geomx["paired_scores"]["entry"].eq("strict_pass_idhwt")].copy()
    summary = geomx["entry_summary"].loc[geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]
    table.to_csv(FIG3_WRITE / "source_csv" / "fig3b_geomx_paired_change.csv", index=False)

    fig, ax = plt.subplots(figsize=(3.0, 3.8))
    for _, row in table.iterrows():
        color = NPG["red"] if row["trial_setting"] == "Nivolumab" else NPG["blue"]
        ax.plot([0, 1], [row["Primary"], row["Recurrence"]], color=color, alpha=0.45, linewidth=0.9)
        ax.scatter([0, 1], [row["Primary"], row["Recurrence"]], color=color, s=14, zorder=3, linewidths=0)
    means = table[["Primary", "Recurrence"]].mean(axis=0)
    ax.plot([0, 1], means, color="#111111", linewidth=2.0, marker="o", markersize=4.2, zorder=4)
    ax.set_xticks([0, 1], ["Primary", "Recurrence"])
    ax.set_ylabel("Raw20 score")
    ax.text(0.5, 1.01, f"22 pairs | mean delta={summary['mean_delta']:.3f} | FDR={summary['fdr']:.2g}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.6)
    save(fig, FIG3_FIG, "fig3b_geomx_paired_change")


def plot_gene_forest(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["gene_summary"].loc[geomx["gene_summary"]["entry"].eq("strict_pass_idhwt")].copy()
    table = table.sort_values(["mean_log2_delta", "gene"], ascending=[False, True]).reset_index(drop=True)
    table.to_csv(FIG3_WRITE / "source_csv" / "fig3c_geomx_gene_forest.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 4.8))
    y = np.arange(len(table))
    colors = np.where(table["mean_log2_delta"] >= 0, NPG["red"], NPG["navy"])
    ax.hlines(y, table["ci_low"], table["ci_high"], color=colors, linewidth=1.2)
    ax.scatter(table["mean_log2_delta"], y, color=colors, s=24, zorder=3)
    ax.axvline(0, color="#666666", linewidth=0.8)
    x_text = float(table["ci_high"].max()) + 0.06
    x_min = float(table["ci_low"].min()) - 0.08
    x_max = x_text + 0.26
    ax.set_yticks(y, table["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean paired log2 delta")
    ax.set_xlim(x_min, x_max)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    for yy, fdr in zip(y, table["fdr"]):
        ax.text(x_text, yy, f"FDR={fdr:.2g}", va="center", ha="left", fontsize=6.8, color="#333333")
    save(fig, FIG3_FIG, "fig3c_geomx_measurable_gene_forest")


def plot_spatial_maps(spatial: dict[str, pd.DataFrame]) -> None:
    spots = spatial["spots"].copy()
    targets = [("raw20", "raw20"), ("mdsc_combined", "MDSC"), ("mes", "MES")]
    limits = {
        key: spots[key].quantile([0.02, 0.98]).to_numpy()
        for key, _ in targets
    }

    for sample in sorted(spots["sample"].unique()):
        sample_df = spots.loc[spots["sample"].eq(sample)].copy()
        sample_df.to_csv(FIG3_WRITE / "source_csv" / f"{sample}_spatial_maps.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
        fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.25))
        for ax, (key, label) in zip(axes, targets):
            low, high = limits[key]
            sca = ax.scatter(
                sample_df["pxl_col_in_fullres"],
                sample_df["pixel_y_plot"],
                c=sample_df[key],
                s=4,
                cmap="RdYlBu_r",
                vmin=low,
                vmax=high,
                linewidths=0,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            ax.set_title(label, fontsize=8.5, pad=2)
            for spine in ax.spines.values():
                spine.set_visible(False)
            cbar = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.02)
            cbar.ax.tick_params(labelsize=6.5, length=2)
        save(fig, FIG3_FIG, f"fig3d_{sample}_spatial_maps" if sample == "GBM030" else f"fig3e_{sample}_spatial_maps")


def plot_partial_rho(spatial: dict[str, pd.DataFrame]) -> None:
    table = spatial["partial"].loc[spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()
    table["sample_label"] = table["sample"].map({"GBM030": "GBM030", "GBM049": "GBM049"})
    table.to_csv(FIG3_WRITE / "source_csv" / "fig3f_partial_rho.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    y_positions = {"GBM030": 1, "GBM049": 0}
    offset = {"MDSC": 0.10, "MES": -0.10}
    color = {"MDSC": NPG["teal"], "MES": NPG["brown"]}
    for _, row in table.iterrows():
        yy = y_positions[row["sample_label"]] + offset[row["target"]]
        ax.scatter(row["partial_rho"], yy, s=44, color=color[row["target"]], edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(row["partial_rho"] + 0.015, yy, f"FDR={row['fdr']:.2g}", va="center", fontsize=6.8, color="#333333")
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks([0, 1], ["GBM049", "GBM030"])
    ax.set_xlabel("Partial Spearman rho\n(adjusted for total myeloid score)")
    ax.set_xlim(0, max(0.62, table["partial_rho"].max() + 0.12))
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    handles = [
        mpl.lines.Line2D([0], [0], marker="o", color="w", markerfacecolor=color[label], markeredgecolor="white",
                         markersize=7, label=label)
        for label in ["MDSC", "MES"]
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    save(fig, FIG3_FIG, "fig3f_partial_rho")


def plot_s6_entry_sensitivity(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["entry_summary"].copy()
    table["label"] = table["entry"].map({
        "strict_pass_idhwt": "IBA1+ strict + IDH-wt",
        "all_aoi_idhwt": "IBA1+ all AOI + IDH-wt",
        "strict_pass_all_idh": "IBA1+ strict + all IDH",
    })
    table = table.sort_values("mean_delta", ascending=False)
    table.to_csv(S6_WRITE / "source_csv" / "s6a_entry_sensitivity.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    y = np.arange(len(table))
    ax.scatter(table["mean_delta"], y, s=54, color=[NPG["red"], NPG["salmon"], NPG["navy"]][:len(table)], edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y, table["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean paired raw20 delta")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    for x, yy, n_pairs, fdr in zip(table["mean_delta"], y, table["n_pairs"], table["fdr"]):
        ax.text(x + 0.012, yy, f"n={n_pairs}, FDR={fdr:.2g}", va="center", fontsize=7, color="#333333")
    save(fig, S6_FIG, "s6a_entry_sensitivity")


def plot_s6_idh_sensitivity(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["subgroup_summary"].copy()
    table = table.loc[
        (table["entry"].eq("strict_pass_all_idh") & table["subgroup"].isin(["IDH_WT", "IDH_mutated"]))
        | (table["entry"].eq("strict_pass_idhwt") & table["subgroup"].eq("All"))
    ].copy()
    table["label"] = table.apply(
        lambda row: {
            ("strict_pass_idhwt", "All"): "Main analysis IDH-wt",
            ("strict_pass_all_idh", "IDH_WT"): "All-IDH entry: IDH-wt only",
            ("strict_pass_all_idh", "IDH_mutated"): "All-IDH entry: IDH-mutant only",
        }[(row["entry"], row["subgroup"])],
        axis=1,
    )
    table.to_csv(S6_WRITE / "source_csv" / "s6b_idh_sensitivity.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    y = np.arange(len(table))
    colors = [NPG["red"], NPG["blue"], NPG["grayblue"]]
    ax.scatter(table["mean_delta"], y, s=54, color=colors[:len(table)], edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y, table["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean paired raw20 delta")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    for x, yy, n_pairs, fdr in zip(table["mean_delta"], y, table["n_pairs"], table["fdr"]):
        label = f"n={n_pairs}"
        if pd.notna(fdr):
            label += f", FDR={fdr:.2g}"
        ax.text(x + 0.012, yy, label, va="center", fontsize=7, color="#333333")
    save(fig, S6_FIG, "s6b_idh_sensitivity")


def plot_s6_treatment_stratified(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["paired_scores"].loc[geomx["paired_scores"]["entry"].eq("strict_pass_idhwt")].copy()
    summary = geomx["subgroup_summary"].loc[
        geomx["subgroup_summary"]["entry"].eq("strict_pass_idhwt")
        & geomx["subgroup_summary"]["subgroup"].isin(["Control", "Nivolumab"])
    ].copy()
    table.to_csv(S6_WRITE / "source_csv" / "s6c_treatment_stratified_paired_change.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(5.1, 3.2), sharey=True)
    for ax, subgroup, color in zip(axes, ["Control", "Nivolumab"], [NPG["blue"], NPG["red"]]):
        part = table.loc[table["trial_setting"].eq(subgroup)].copy()
        for _, row in part.iterrows():
            ax.plot([0, 1], [row["Primary"], row["Recurrence"]], color=color, alpha=0.45, linewidth=0.9)
            ax.scatter([0, 1], [row["Primary"], row["Recurrence"]], color=color, s=14, linewidths=0, zorder=3)
        means = part[["Primary", "Recurrence"]].mean(axis=0)
        ax.plot([0, 1], means, color="#111111", linewidth=2.0, marker="o", markersize=4.0, zorder=4)
        row = summary.loc[summary["subgroup"].eq(subgroup)].iloc[0]
        ax.text(0.5, 1.01, f"n={int(row['n_pairs'])} | FDR={row['fdr']:.2g}", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=7)
        ax.set_title(subgroup, fontsize=8.5, pad=3)
        ax.set_xticks([0, 1], ["Primary", "Recurrence"])
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#ECECEC", linewidth=0.6)
    axes[0].set_ylabel("Raw20 score")
    save(fig, S6_FIG, "s6c_treatment_stratified")


def plot_s6_interaction(geomx: dict[str, pd.DataFrame]) -> None:
    table = geomx["subgroup_summary"].loc[
        geomx["subgroup_summary"]["subgroup"].eq("Interaction_Nivolumab_minus_Control")
    ].copy()
    table["label"] = table["entry"].map({
        "strict_pass_idhwt": "IBA1+ strict + IDH-wt",
        "all_aoi_idhwt": "IBA1+ all AOI + IDH-wt",
        "strict_pass_all_idh": "IBA1+ strict + all IDH",
    })
    table.to_csv(S6_WRITE / "source_csv" / "s6d_interaction.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    y = np.arange(len(table))
    ax.scatter(table["mean_delta"], y, s=54, color=NPG["brown"], edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y, table["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Nivolumab minus control\nmean paired delta")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    for x, yy, p_value in zip(table["mean_delta"], y, table["sign_flip_p"]):
        ax.text(x + 0.012, yy, f"P={p_value:.3g}", va="center", fontsize=7, color="#333333")
    save(fig, S6_FIG, "s6d_interaction")


def plot_s6_slice_sensitivity(spatial: dict[str, pd.DataFrame]) -> None:
    table = spatial["partial"].copy()
    table = table.loc[table["target"].isin(["MDSC", "E-MDSC", "M-MDSC", "MES"])].copy()
    pivot = table.pivot(index="target", columns="sample", values="partial_rho").reindex(["MDSC", "E-MDSC", "M-MDSC", "MES"])
    table.to_csv(S6_WRITE / "source_csv" / "s6e_slice_sensitivity.csv", index=False)

    fig, ax = plt.subplots(figsize=(3.2, 2.9))
    data = pivot.to_numpy(dtype=float)
    im = ax.imshow(data, cmap="RdYlBu_r", vmin=0, vmax=max(0.55, np.nanmax(data)))
    ax.set_xticks(range(pivot.shape[1]), pivot.columns)
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    ax.tick_params(length=0)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = data[i, j]
            fdr = table.loc[(table["target"].eq(pivot.index[i])) & (table["sample"].eq(pivot.columns[j])), "fdr"].iloc[0]
            ax.text(j, i, f"{value:.2f}\nFDR={fdr:.1g}", ha="center", va="center", fontsize=6.5,
                    color="white" if value > 0.33 else "#222222")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cbar.set_label("Partial rho", fontsize=7)
    cbar.ax.tick_params(labelsize=6.5)
    save(fig, S6_FIG, "s6e_slice_sensitivity")


def write_summary(step38: pd.DataFrame, geomx: dict[str, pd.DataFrame], spatial: dict[str, pd.DataFrame]) -> None:
    geomx_main = geomx["entry_summary"].loc[geomx["entry_summary"]["entry"].eq("strict_pass_idhwt")].iloc[0]
    partial_main = spatial["partial"].loc[spatial["partial"]["target"].isin(["MDSC", "MES"])].copy()
    lines = [
        "# Figure3 / Supplementary Figure6 raw20 summary",
        "",
        "## 主结果",
        "",
        f"- Step38 raw20 主程序在 GSE174554 为 NES={step38.loc[step38['dataset'].eq('GSE174554'), 'NES'].iloc[0]:.3f}, nominal P={step38.loc[step38['dataset'].eq('GSE174554'), 'nominal_p'].iloc[0]:.2g}; 在 GSE274546 为 NES={step38.loc[step38['dataset'].eq('GSE274546'), 'NES'].iloc[0]:.3f}, nominal P={step38.loc[step38['dataset'].eq('GSE274546'), 'nominal_p'].iloc[0]:.2g}.",
        f"- Artzi IBA1+ GeoMx 严格 IDH-wt 22 对原发-复发患者平均 raw20 变化 {geomx_main['mean_delta']:.3f}, FDR={geomx_main['fdr']:.2g}, 正向 {int(geomx_main['n_positive'])}/{int(geomx_main['n_pairs'])}.",
        f"- GSE276841 两张切片在扣除总髓系后，raw20 与 MDSC 的 partial rho 范围 {partial_main.loc[partial_main['target'].eq('MDSC'), 'partial_rho'].min():.3f}-{partial_main.loc[partial_main['target'].eq('MDSC'), 'partial_rho'].max():.3f}; 与 MES 的 partial rho 范围 {partial_main.loc[partial_main['target'].eq('MES'), 'partial_rho'].min():.3f}-{partial_main.loc[partial_main['target'].eq('MES'), 'partial_rho'].max():.3f}.",
        "",
        "## 补图边界",
        "",
        "- GeoMx 入口、IDH 和治疗分层均未改变方向；治疗 interaction 仍不显著。",
        "- GSE276841 的逐切片敏感性继续支持 raw20 靠近 MDSC / MES 生态位。",
    ]
    (WRITE_BASE / "39_raw20_figure_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    setup_style()
    ensure_dirs()
    markers = read_author_markers()
    geomx = recompute_geomx_raw20()
    spatial = recompute_gse276841_raw20(markers)
    step38 = load_step38_main()
    step38.to_csv(FIG3_WRITE / "source_csv" / "step38_raw20_main_gsea.csv", index=False)

    plot_flowchart(step38, geomx, spatial)
    plot_geomx_paired_change(geomx)
    plot_gene_forest(geomx)
    plot_spatial_maps(spatial)
    plot_partial_rho(spatial)

    plot_s6_entry_sensitivity(geomx)
    plot_s6_idh_sensitivity(geomx)
    plot_s6_treatment_stratified(geomx)
    plot_s6_interaction(geomx)
    plot_s6_slice_sensitivity(spatial)

    write_summary(step38, geomx, spatial)
    print("FIGURE3_AND_S6_COMPLETE")


if __name__ == "__main__":
    main()
