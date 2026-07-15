#!/usr/bin/env python3
"""重建正式 Figure 2：身份可检验性与 GSE278456 真正 raw20 状态定位。"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from PIL import Image, ImageDraw, ImageFont
from scipy import sparse
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
H5_ROOT = EXTERNAL_ROOT / "write/30_dataset_rescue_search/source_metadata/GSE278456_myeloid_h5"
MARKER_BOOK = EXTERNAL_ROOT / "write/30_dataset_rescue_search/source_metadata/NIHMS2115262-supplement-S2_S5_S6.xls"
STEP30 = ROOT / "write/30_dataset_rescue_search"
STEP39_IDENTITY = ROOT / "write/39_independent_miller_mg_inflammatory_figures/SupplementaryFigure4/source_data"

WRITE_ROOT = ROOT / "write/48_figure2_identity_state_rebuild/Figure2"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_ROOT = ROOT / "figures/48_figure2_identity_state_rebuild/Figure2"

RAW20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
MYELOID_STATES = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]
IDENTITY_ORDER = ["Macrophage", "Microglia", "Monocyte", "CDC", "Neutrophil"]

BLUE = "#3C5488"
RED = "#E64B35"
TEAL = "#00A087"
CYAN = "#4DBBD5"
SALMON = "#F39B7F"
GREY_BLUE = "#8491B4"
MINT = "#91D1C2"
BROWN = "#7E6148"
BLACK = "#222222"
GREY = "#8A8A8A"
LIGHT_GREY = "#D9D9D9"
GRID = "#E8E8E8"
DATASET_COLORS = {"GSE174554": BLUE, "GSE274546": RED}
STATE_COLORS = {
    "MCG1": TEAL,
    "MCG2": MINT,
    "MCG3": CYAN,
    "MCG4": GREY_BLUE,
    "MCG5": BLUE,
    "MAC1": SALMON,
    "MAC2": RED,
    "M-MDSC": BROWN,
    "E-MDSC": "#B09C85",
}
MM = 1 / 25.4


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.2,
            "axes.labelsize": 7.6,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 6.2,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    sc.settings.verbosity = 0


def clean_axes(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(BLACK)
    ax.tick_params(colors=BLACK, length=2.4, width=0.7)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6, zorder=0)


def save_panel(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    pdf = FIG_ROOT / f"{stem}.pdf"
    png = FIG_ROOT / f"{stem}.png"
    fig.savefig(
        pdf,
        bbox_inches="tight",
        pad_inches=0.035,
        metadata={"Creator": "Step48 Figure2", "CreationDate": None, "ModDate": None},
    )
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    return pdf, png


def write_source(frame: pd.DataFrame, stem: str) -> Path:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    path = SOURCE_OUT / f"{stem}_source.csv"
    frame.to_csv(path, index=False)
    return path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_author_markers() -> tuple[dict[str, list[str]], pd.DataFrame]:
    workbook = pd.ExcelFile(MARKER_BOOK)
    markers: dict[str, list[str]] = {}
    rows: list[dict[str, object]] = []
    raw20_set = set(RAW20)
    for sheet in workbook.sheet_names:
        table = pd.read_excel(MARKER_BOOK, sheet_name=sheet)
        if "FDR" not in table.columns:
            table = pd.read_excel(MARKER_BOOK, sheet_name=sheet, header=1)
        gene_col = table.columns[0]
        original = table.loc[
            pd.to_numeric(table["FDR"], errors="coerce").lt(0.05)
            & pd.to_numeric(table["Foldchange"], errors="coerce").gt(0),
            gene_col,
        ].astype(str).tolist()[:40]
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
    missing = sorted(set(MYELOID_STATES) - set(markers))
    if missing:
        raise RuntimeError(f"作者marker工作簿缺少预设髓系状态: {missing}")
    return markers, pd.DataFrame(rows)


def exact_two_sided_signflip(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = len(values)
    if n < 3:
        return np.nan
    observed = abs(values.mean())
    # 双侧绝对值分布具有整体符号对称性，固定首项为正可精确减半枚举空间。
    tail = values[1:]
    total = 1 << (n - 1)
    extreme = 0
    chunk_size = 65536
    bit_positions = np.arange(n - 1, dtype=np.uint64)
    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        codes = np.arange(start, stop, dtype=np.uint64)[:, None]
        signs = (((codes >> bit_positions) & 1).astype(np.float64) * 2.0) - 1.0
        means = (values[0] + signs @ tail) / n
        extreme += int(np.count_nonzero(np.abs(means) >= observed - 1e-12))
    return extreme / total


def bootstrap_mean_ci(values: np.ndarray, seed: int) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    mean = float(values.mean())
    if len(values) == 1:
        return mean, mean, mean
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(20000, len(values)), replace=True).mean(axis=1)
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return mean, float(lo), float(hi)


def dense_vector(matrix: object) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def aggregate_raw20_pseudobulk_score(
    counts: sparse.spmatrix | np.ndarray,
    raw20_indices: np.ndarray,
    mask: np.ndarray,
) -> float:
    subset = counts[mask]
    library = float(subset.sum())
    if library <= 0:
        return np.nan
    gene_counts = np.asarray(subset[:, raw20_indices].sum(axis=0)).ravel()
    log_cpm = np.log2(gene_counts / library * 1_000_000 + 1)
    return float(log_cpm.mean())


def process_gse278456(
    markers: dict[str, list[str]],
) -> dict[str, pd.DataFrame]:
    manifest_all = pd.read_csv(STEP30 / "GSE278456_sample_manifest.csv")
    manifest_all["geo_accession"] = manifest_all["geo_accession"].astype(str)
    existing_qc = pd.read_csv(STEP30 / "GSE278456_single_cell_qc.csv").set_index("geo_accession")
    paths = sorted(H5_ROOT.glob("*.h5"))
    if len(paths) != 21:
        raise RuntimeError(f"GSE278456锁定H5数量应为21，实际为{len(paths)}")

    all_marker_genes = set().union(*markers.values())
    cell_frames: list[pd.DataFrame] = []
    state_score_frames: list[pd.DataFrame] = []
    effect_rows: list[dict[str, object]] = []
    gene_patient_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    locked_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []

    for file_index, path in enumerate(paths, start=1):
        accession = path.name.split("_")[0]
        manifest_row = manifest_all.loc[manifest_all["geo_accession"].eq(accession)]
        if len(manifest_row) != 1:
            raise RuntimeError(f"{accession}无法唯一匹配样本manifest")
        manifest_row = manifest_row.iloc[0]
        if not bool(manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"锁定H5中出现非IDH-wt GBM: {accession}")
        sample = str(manifest_row["sample_title"])
        print(f"[{file_index:02d}/21] {accession} {sample}", flush=True)

        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        n_input = adata.n_obs
        n_genes = adata.n_vars
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].between(200, 8000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        adata = adata[keep].copy()
        n_qc = adata.n_obs
        if accession not in existing_qc.index:
            raise RuntimeError(f"旧QC表缺少{accession}")
        old = existing_qc.loc[accession]
        if int(old["n_input"]) != n_input or int(old["n_qc_pass"]) != n_qc:
            raise RuntimeError(f"{accession} QC细胞数未复现")

        present_raw20 = [gene for gene in RAW20 if gene in adata.var_names]
        missing_raw20 = [gene for gene in RAW20 if gene not in adata.var_names]
        if len(present_raw20) != 19 or missing_raw20 != ["AC253572.2"]:
            raise RuntimeError(f"{accession} raw20覆盖异常: {len(present_raw20)}/20, missing={missing_raw20}")
        coverage_rows.append(
            {
                "geo_accession": accession,
                "sample": sample,
                "n_defined": len(RAW20),
                "n_present": len(present_raw20),
                "present_genes": ";".join(present_raw20),
                "missing_genes": ";".join(missing_raw20),
            }
        )
        qc_rows.append(
            {
                "sample": sample,
                "geo_accession": accession,
                "n_input": n_input,
                "n_qc_pass": n_qc,
                "fraction_qc_pass": n_qc / n_input,
            }
        )
        locked_rows.append(
            {
                "geo_accession": accession,
                "sample": sample,
                "path": str(path),
                "sha256": file_sha256(path),
                "n_input_cells": n_input,
                "n_qc_pass_cells": n_qc,
                "n_genes": n_genes,
            }
        )

        counts = adata.X.copy()
        raw20_indices = adata.var_names.get_indexer(present_raw20)
        adata.obs_names = [f"{accession}:{barcode}" for barcode in adata.obs_names.astype(str)]
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)

        excluded_controls = set(RAW20) | all_marker_genes
        neutral_control_pool = [gene for gene in adata.var_names if gene not in excluded_controls]
        if len(neutral_control_pool) < 5000:
            raise RuntimeError(f"{accession} state-neutral control pool过小")

        score_columns: list[str] = []
        for state, genes in markers.items():
            use = [gene for gene in genes if gene in adata.var_names]
            if len(use) < 30:
                raise RuntimeError(f"{accession} {state}可测marker少于30")
            column = f"state_score_{state}"
            sc.tl.score_genes(
                adata,
                gene_list=use,
                gene_pool=neutral_control_pool + use,
                ctrl_size=50,
                n_bins=25,
                score_name=column,
                random_state=20260713,
                ctrl_as_ref=False,
                use_raw=False,
            )
            score_columns.append(column)

        score_matrix = adata.obs[score_columns].to_numpy(float)
        order = np.argsort(score_matrix, axis=1)
        state_names = np.array(list(markers), dtype=object)
        winner_index = order[:, -1]
        runner_index = order[:, -2]
        winners = state_names[winner_index]
        winning_score = score_matrix[np.arange(adata.n_obs), winner_index]
        runner_score = score_matrix[np.arange(adata.n_obs), runner_index]

        sc.tl.score_genes(
            adata,
            gene_list=present_raw20,
            gene_pool=neutral_control_pool + present_raw20,
            ctrl_size=50,
            n_bins=25,
            score_name="raw20_score",
            random_state=20260713,
            ctrl_as_ref=False,
            use_raw=False,
        )
        raw20_scores = adata.obs["raw20_score"].to_numpy(float)
        focus = np.isin(winners, MYELOID_STATES)

        cell_frames.append(
            pd.DataFrame(
                {
                    "cell_id": adata.obs_names.astype(str),
                    "geo_accession": accession,
                    "sample": sample,
                    "author_marker_state": winners,
                    "state_score": winning_score,
                    "runner_up_state": state_names[runner_index],
                    "state_margin": winning_score - runner_score,
                    "raw20_score": raw20_scores,
                    "in_prespecified_myeloid_states": focus,
                }
            )
        )

        normalized_raw20 = adata[:, present_raw20].X
        if sparse.issparse(normalized_raw20):
            normalized_raw20 = normalized_raw20.toarray()
        normalized_raw20 = np.asarray(normalized_raw20, dtype=np.float32)

        for state in MYELOID_STATES:
            target = focus & (winners == state)
            rest = focus & (winners != state)
            n_target = int(target.sum())
            n_rest = int(rest.sum())
            if n_target >= 30:
                state_score_frames.append(
                    pd.DataFrame(
                        {
                            "geo_accession": [accession],
                            "sample": [sample],
                            "state": [state],
                            "n_cells": [n_target],
                            "mean_raw20_score": [float(raw20_scores[target].mean())],
                        }
                    )
                )
                state_expr = normalized_raw20[target]
                for gene_index, gene in enumerate(present_raw20):
                    gene_patient_rows.append(
                        {
                            "geo_accession": accession,
                            "sample": sample,
                            "state": state,
                            "gene": gene,
                            "n_cells": n_target,
                            "mean_log_expression": float(state_expr[:, gene_index].mean()),
                            "pct_detected": float((state_expr[:, gene_index] > 0).mean() * 100),
                        }
                    )
            if n_target >= 30 and n_rest >= 30:
                delta = float(raw20_scores[target].mean() - raw20_scores[rest].mean())
                target_pb = aggregate_raw20_pseudobulk_score(counts, raw20_indices, target)
                rest_pb = aggregate_raw20_pseudobulk_score(counts, raw20_indices, rest)
                effect_rows.append(
                    {
                        "geo_accession": accession,
                        "sample": sample,
                        "state": state,
                        "n_target": n_target,
                        "n_rest": n_rest,
                        "raw20_delta": delta,
                        "raw20_pseudobulk_delta": target_pb - rest_pb,
                    }
                )

        del adata, counts, normalized_raw20

    locked = pd.DataFrame(locked_rows).sort_values("geo_accession")
    qc = pd.DataFrame(qc_rows).sort_values("geo_accession")
    coverage = pd.DataFrame(coverage_rows).sort_values("geo_accession")
    cells = pd.concat(cell_frames, ignore_index=True)
    state_scores = pd.concat(state_score_frames, ignore_index=True)
    effects = pd.DataFrame(effect_rows)
    gene_patient = pd.DataFrame(gene_patient_rows)

    if locked["geo_accession"].nunique() != 21 or locked["sample"].nunique() != 21:
        raise RuntimeError("GSE278456锁定manifest不是21个唯一GSM/患者")
    if int(qc["n_input"].sum()) != 131780 or int(qc["n_qc_pass"].sum()) != 120766:
        raise RuntimeError("GSE278456总细胞数未复现Step30")
    if cells["cell_id"].duplicated().any():
        raise RuntimeError("GSE278456 cell ID不唯一")

    summary_rows = []
    for index, state in enumerate(MYELOID_STATES):
        part = effects.loc[effects["state"].eq(state)].copy()
        values = part["raw20_delta"].to_numpy(float)
        mean, low, high = bootstrap_mean_ci(values, 20260713 + index)
        pb_values = part["raw20_pseudobulk_delta"].to_numpy(float)
        summary_rows.append(
            {
                "state": state,
                "n_patients": len(values),
                "mean_delta": mean,
                "ci95_low": low,
                "ci95_high": high,
                "median_delta": float(np.median(values)) if len(values) else np.nan,
                "n_positive": int((values > 0).sum()),
                "exact_two_sided_signflip_p": exact_two_sided_signflip(values),
                "mean_pseudobulk_delta": float(pb_values.mean()) if len(pb_values) else np.nan,
                "pseudobulk_n_positive": int((pb_values > 0).sum()),
                "score_pseudobulk_mean_direction_concordant": bool(np.sign(mean) == np.sign(pb_values.mean())) if len(values) else False,
            }
        )
    effect_summary = pd.DataFrame(summary_rows)
    ok = effect_summary["exact_two_sided_signflip_p"].notna()
    effect_summary["fdr_9_states"] = np.nan
    effect_summary.loc[ok, "fdr_9_states"] = multipletests(
        effect_summary.loc[ok, "exact_two_sided_signflip_p"], method="fdr_bh"
    )[1]
    gene_summary = (
        gene_patient.groupby(["state", "gene"], observed=True)
        .agg(
            n_patients=("sample", "nunique"),
            patient_balanced_mean_log_expression=("mean_log_expression", "mean"),
            patient_balanced_pct_detected=("pct_detected", "mean"),
        )
        .reset_index()
    )
    gene_summary["within_gene_z"] = gene_summary.groupby("gene", observed=True)[
        "patient_balanced_mean_log_expression"
    ].transform(lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0)

    return {
        "locked_manifest": locked,
        "qc": qc,
        "coverage": coverage,
        "cells": cells,
        "patient_state_scores": state_scores,
        "patient_state_effects": effects,
        "state_effect_summary": effect_summary,
        "gene_patient": gene_patient,
        "gene_summary": gene_summary,
    }


def load_identity_sources() -> dict[str, pd.DataFrame]:
    pair_counts = pd.read_csv(STEP39_IDENTITY / "identity_pair_testability.csv")
    summary = pd.read_csv(STEP39_IDENTITY / "identity_raw20_patient_delta_summary.csv")
    deltas = pd.read_csv(STEP39_IDENTITY / "identity_raw20_patient_deltas.csv")
    threshold20 = pair_counts.loc[pair_counts["threshold"].eq(20)].copy()
    merged = threshold20.merge(summary, on=["dataset", "identity"], how="left", validate="one_to_one")
    if int(merged.loc[merged["dataset"].eq("GSE174554") & merged["identity"].eq("Macrophage"), "n_pairs_x"].iloc[0]) != 18:
        raise RuntimeError("GSE174554 Macrophage配对数异常")
    if int(merged.loc[merged["dataset"].eq("GSE274546") & merged["identity"].eq("Macrophage"), "n_pairs_x"].iloc[0]) != 39:
        raise RuntimeError("GSE274546 Macrophage配对数异常")
    merged = merged.rename(columns={"n_pairs_x": "n_pairs"}).drop(columns=["n_pairs_y"])
    return {"testability": merged, "deltas": deltas}


def panel_a_identity_testability(frame: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2A_identity_raw20_testability"
    source = write_source(frame, stem)
    datasets = ["GSE174554", "GSE274546"]
    fig, ax = plt.subplots(figsize=(86 * MM, 60 * MM))
    norm = mpl.colors.TwoSlopeNorm(vmin=-0.10, vcenter=0, vmax=0.10)
    cmap = mpl.colormaps["RdBu_r"]
    for yi, identity in enumerate(IDENTITY_ORDER):
        for xi, dataset in enumerate(datasets):
            row = frame.loc[frame["dataset"].eq(dataset) & frame["identity"].eq(identity)].iloc[0]
            n_pairs = int(row["n_pairs"])
            if n_pairs == 0:
                ax.scatter(xi, yi, marker="x", s=22, color="#BDBDBD", linewidth=0.8)
            else:
                value = float(row["mean_delta"]) if np.isfinite(row["mean_delta"]) else 0
                ax.scatter(
                    xi,
                    yi,
                    s=42 + n_pairs * 7.5,
                    color=cmap(norm(np.clip(value, -0.10, 0.10))),
                    edgecolor="white",
                    linewidth=0.7,
                )
            ax.text(xi, yi, str(n_pairs), ha="center", va="center", fontsize=6.0, color="white" if n_pairs >= 10 else BLACK)
    ax.set_xticks(range(2), datasets)
    ax.set_yticks(range(len(IDENTITY_ORDER)), IDENTITY_ORDER)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(len(IDENTITY_ORDER) - 0.5, -0.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.06)
    cbar.set_label("Mean recurrent - primary raw20 score", fontsize=6.4)
    cbar.ax.tick_params(labelsize=5.8, length=2)
    fig.tight_layout(pad=0.45)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2A", "description": "Identity-resolved testability at 20 cells per endpoint", "source": str(source), "pdf": str(pdf), "png": str(png)}


def panel_b_macrophage_deltas(deltas: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2B_macrophage_patient_raw20_deltas"
    frame = deltas.loc[deltas["identity"].eq("Macrophage")].copy()
    source = write_source(frame, stem)
    datasets = ["GSE174554", "GSE274546"]
    fig, ax = plt.subplots(figsize=(86 * MM, 64 * MM))
    rng = np.random.default_rng(20260713)
    for index, dataset in enumerate(datasets):
        values = frame.loc[frame["dataset"].eq(dataset), "delta"].to_numpy(float)
        color = DATASET_COLORS[dataset]
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        low, high = np.min(values), np.max(values)
        ax.add_patch(mpl.patches.Rectangle((index - 0.24, q1), 0.48, q3 - q1, facecolor=color, alpha=0.18, edgecolor=color, linewidth=0.9))
        ax.plot([index - 0.24, index + 0.24], [median, median], color=color, linewidth=1.2)
        ax.plot([index, index], [low, q1], color=color, linewidth=0.7)
        ax.plot([index, index], [q3, high], color=color, linewidth=0.7)
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(index + jitter, values, s=13, color=color, alpha=0.75, edgecolor="white", linewidth=0.3, zorder=3)
        mean, ci_low, ci_high = bootstrap_mean_ci(values, 20260713 + index)
        ax.vlines(index + 0.30, ci_low, ci_high, color=BLACK, linewidth=1.1, zorder=4)
        ax.scatter(index + 0.30, mean, marker="D", s=18, color=BLACK, edgecolor="white", linewidth=0.35, zorder=5)
        ax.text(
            index,
            max(values) + 0.035,
            f"{int((values > 0).sum())}/{len(values)} higher\nmean {mean:+.3f}",
            ha="center",
            va="bottom",
            fontsize=6.1,
            color=color,
        )
    ax.axhline(0, color="#888888", linewidth=0.7)
    ax.set_xticks(range(2), datasets)
    ax.set_ylabel("Patient raw20 score delta")
    clean_axes(ax, grid_axis="y")
    ax.set_axisbelow(True)
    ymax = frame["delta"].max() + 0.16
    ymin = frame["delta"].min() - 0.06
    ax.set_ylim(ymin, ymax)
    fig.tight_layout(pad=0.55)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2B", "description": "Descriptive raw20 change within Macrophage-labelled compartments", "source": str(source), "pdf": str(pdf), "png": str(png)}


def panel_c_state_scores(scores: pd.DataFrame, effect_summary: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2C_GSE278456_patient_state_raw20_scores"
    order = effect_summary.sort_values("mean_delta", ascending=True)["state"].tolist()
    frame = scores.loc[scores["state"].isin(order)].copy()
    frame["state"] = pd.Categorical(frame["state"], categories=order, ordered=True)
    source = write_source(frame, stem)
    fig, ax = plt.subplots(figsize=(92 * MM, 72 * MM))
    for yi, state in enumerate(order):
        values = frame.loc[frame["state"].eq(state), "mean_raw20_score"].to_numpy(float)
        if not len(values):
            continue
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.add_patch(mpl.patches.Rectangle((q1, yi - 0.22), q3 - q1, 0.44, facecolor=STATE_COLORS[state], alpha=0.20, edgecolor=STATE_COLORS[state], linewidth=0.8))
        ax.plot([median, median], [yi - 0.22, yi + 0.22], color=STATE_COLORS[state], linewidth=1.2)
        rng = np.random.default_rng(20260713 + yi)
        jitter = rng.uniform(-0.15, 0.15, size=len(values))
        ax.scatter(values, yi + jitter, s=12, color=STATE_COLORS[state], alpha=0.72, edgecolor="white", linewidth=0.25)
        ax.text(max(values) + 0.02, yi, f"n={len(values)}", va="center", fontsize=5.8, color=GREY)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("Patient-level Miller raw20 score")
    clean_axes(ax, grid_axis="x")
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.55)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2C", "description": "Patient-balanced raw20 scores across author-marker-derived states", "source": str(source), "pdf": str(pdf), "png": str(png)}


def panel_d_state_effects(summary: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2D_GSE278456_state_vs_rest_effects"
    frame = summary.sort_values("mean_delta", ascending=True).reset_index(drop=True)
    source = write_source(frame, stem)
    fig, ax = plt.subplots(figsize=(102 * MM, 72 * MM))
    for index, row in frame.iterrows():
        significant = bool(row["fdr_9_states"] < 0.05)
        color = STATE_COLORS[row["state"]] if significant else "#AFAFAF"
        ax.hlines(index, row["ci95_low"], row["ci95_high"], color=color, linewidth=1.25)
        ax.scatter(row["mean_delta"], index, s=30, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        fdr = row["fdr_9_states"]
        label = f"{int(row['n_positive'])}/{int(row['n_patients'])}  FDR {fdr:.3g}"
        ax.text(max(frame["ci95_high"].max(), 0) + 0.025, index, label, va="center", fontsize=5.8, color=BLACK)
    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.set_yticks(range(len(frame)), frame["state"])
    ax.set_xlabel("Within-patient raw20 delta (state - other myeloid states)")
    xmin = min(float(frame["ci95_low"].min()), 0) - 0.03
    xmax = max(float(frame["ci95_high"].max()), 0) + 0.20
    ax.set_xlim(xmin, xmax)
    clean_axes(ax, grid_axis="x")
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.55)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2D", "description": "GSE278456 patient-level state-versus-rest raw20 effects", "source": str(source), "pdf": str(pdf), "png": str(png)}


def panel_e_patient_heatmap(effects: pd.DataFrame, summary: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2E_GSE278456_patient_state_effect_heatmap"
    state_order = summary.sort_values("mean_delta", ascending=False)["state"].tolist()
    matrix = effects.pivot(index="sample", columns="state", values="raw20_delta").reindex(columns=state_order)
    patient_order = matrix.mean(axis=1, skipna=True).sort_values(ascending=False).index
    matrix = matrix.reindex(patient_order)
    source_frame = matrix.reset_index().melt(id_vars="sample", var_name="state", value_name="raw20_delta")
    source = write_source(source_frame, stem)

    fig, ax = plt.subplots(figsize=(102 * MM, 78 * MM))
    values = matrix.to_numpy(float)
    vmax = float(np.nanquantile(np.abs(values), 0.98))
    image = ax.imshow(np.ma.masked_invalid(values), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(state_order)), state_order, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.tick_params(length=0)
    ax.set_facecolor("#EFEFEF")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("State - rest raw20 delta", fontsize=6.4)
    cbar.ax.tick_params(labelsize=5.7, length=2)
    fig.tight_layout(pad=0.45)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2E", "description": "Patient-level stability of GSE278456 state effects", "source": str(source), "pdf": str(pdf), "png": str(png)}


def panel_f_gene_state_bubble(gene_summary: pd.DataFrame, effect_summary: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2F_GSE278456_raw20_gene_state_matrix"
    state_order = effect_summary.sort_values("mean_delta", ascending=False)["state"].tolist()
    frame = gene_summary.copy()
    frame["state"] = pd.Categorical(frame["state"], categories=state_order, ordered=True)
    frame["gene"] = pd.Categorical(frame["gene"], categories=[gene for gene in RAW20 if gene != "AC253572.2"][::-1], ordered=True)
    source = write_source(frame, stem)

    fig, ax = plt.subplots(figsize=(132 * MM, 98 * MM))
    x_lookup = {state: index for index, state in enumerate(state_order)}
    gene_order = [gene for gene in RAW20 if gene != "AC253572.2"][::-1]
    y_lookup = {gene: index for index, gene in enumerate(gene_order)}
    x = frame["state"].astype(str).map(x_lookup).to_numpy(float)
    y = frame["gene"].astype(str).map(y_lookup).to_numpy(float)
    sizes = 6 + frame["patient_balanced_pct_detected"].to_numpy(float) * 0.55
    points = ax.scatter(
        x,
        y,
        s=sizes,
        c=frame["within_gene_z"],
        cmap="RdBu_r",
        norm=mpl.colors.Normalize(vmin=-1.8, vmax=1.8),
        edgecolor="white",
        linewidth=0.25,
    )
    ax.set_xticks(range(len(state_order)), state_order, rotation=45, ha="right")
    ax.set_yticks(range(len(gene_order)), gene_order)
    ax.set_xlim(-0.6, len(state_order) - 0.4)
    ax.set_ylim(-0.6, len(gene_order) - 0.4)
    ax.tick_params(length=0)
    ax.grid(color="#F0F0F0", linewidth=0.5)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(points, ax=ax, fraction=0.028, pad=0.025)
    cbar.set_label("Within-gene z score", fontsize=6.4)
    cbar.ax.tick_params(labelsize=5.7, length=2)
    handles = [
        ax.scatter([], [], s=6 + pct * 0.55, color="#8A8A8A", edgecolor="white", linewidth=0.25, label=f"{pct}%")
        for pct in [10, 30, 50]
    ]
    ax.legend(
        handles=handles,
        title="Detected",
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=3,
        borderaxespad=0,
        handletextpad=0.4,
        columnspacing=0.8,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1), pad=0.5)
    pdf, png = save_panel(fig, stem)
    return {"panel": "Figure2F", "description": "Patient-balanced expression of all 19 measurable raw20 genes", "source": str(source), "pdf": str(pdf), "png": str(png)}


def make_preview(records: list[dict[str, str]]) -> Path:
    paths = {record["panel"]: Path(record["png"]) for record in records}
    rows = [["Figure2A", "Figure2B"], ["Figure2C", "Figure2D"], ["Figure2E", "Figure2F"]]
    cell_w, cell_h = 1500, 1120
    margin, gutter = 85, 50
    canvas = Image.new("RGB", (margin * 2 + 2 * cell_w + gutter, margin * 2 + 3 * cell_h + 2 * gutter), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("Arial Bold.ttf", 52)
    except OSError:
        font = ImageFont.load_default()
    for row_index, panels in enumerate(rows):
        for column_index, panel in enumerate(panels):
            image = Image.open(paths[panel]).convert("RGB")
            image.thumbnail((cell_w - 45, cell_h - 45), Image.Resampling.LANCZOS)
            x0 = margin + column_index * (cell_w + gutter)
            y0 = margin + row_index * (cell_h + gutter)
            x = x0 + (cell_w - image.width) // 2
            y = y0 + (cell_h - image.height) // 2
            canvas.paste(image, (x, y))
            draw.text((x0 + 3, y0 + 2), panel[-1], fill=BLACK, font=font)
    path = FIG_ROOT / "Figure2_formal_candidate_preview.png"
    canvas.save(path, dpi=(180, 180))
    return path


def write_legend(effect_summary: pd.DataFrame) -> Path:
    ordered = effect_summary.sort_values("mean_delta", ascending=False)
    top = ordered.iloc[0]
    text = f"""# Figure 2 legend

**Working title:** The recurrence-associated Miller program forms a continuous myeloid functional axis and preferentially maps to defined author-marker-derived states.

- **A** Number of patients testable within each myeloid identity at a fixed minimum of 20 cells per endpoint. Circle color shows the descriptive mean recurrent-minus-primary raw20 score; `0` indicates insufficient paired cells, not absence of the identity.
- **B** Patient-level raw20 score changes within Macrophage-labelled compartments. Each point is one paired patient; diamonds and vertical lines show the mean and patient-bootstrap 95% confidence interval. These are secondary cell-score summaries, whereas Figure 1 provides the formal raw-count pseudobulk GSEA evidence.
- **C** GSE278456 patient-balanced Miller raw20 scores across all nine prespecified myeloid states. States were assigned from published author markers after removing every measurable raw20 gene from the state-marker sets.
- **D** Within-patient state-versus-rest raw20 effects in GSE278456. Both the target state and the remaining eight myeloid states required at least 30 cells. Points show patient-equal means and bootstrap 95% confidence intervals; FDR is BH-adjusted across all nine states using exact two-sided sign-flip tests.
- **E** Patient-level state-versus-rest raw20 effects. Grey cells denote an insufficient target or rest cell count, not a zero effect.
- **F** Expression of all 19 measurable Miller raw20 genes across the nine states. Color shows the patient-balanced within-gene z score and size shows the patient-balanced detected-cell percentage; `AC253572.2` was absent from all 21 matrices and was not imputed.

GSE278456 included 21 primary IDH-wildtype GBM myeloid-enriched samples. The top state was `{top['state']}` (mean within-patient delta {top['mean_delta']:+.3f}; {int(top['n_positive'])}/{int(top['n_patients'])} positive; FDR={top['fdr_9_states']:.3g}). GSE278456 provides cross-sectional state localization, not recurrence validation. The data do not support a shared pure discrete Mg-inflammatory cluster or a Microglia-specific recurrence effect.
"""
    path = WRITE_ROOT / "Figure2_legend.md"
    path.write_text(text, encoding="utf-8")
    return path


def write_summary(results: dict[str, pd.DataFrame], identity: dict[str, pd.DataFrame]) -> Path:
    state = results["state_effect_summary"].sort_values("mean_delta", ascending=False)
    macrophage = identity["deltas"].loc[identity["deltas"]["identity"].eq("Macrophage")]
    lines = [
        "# Step48 Figure2 result",
        "",
        "- 正式Figure2不使用UMAP，也不再把连续程序包装成跨队列共同离散簇。",
        "- 当前独立重建入口在20细胞/端点门槛下，Macrophage可检验18/39对，Microglia仅0/1对；因此不能主张Microglia特异复发效应。",
        f"- Macrophage患者层raw20细胞分数：GSE174554 {int((macrophage.loc[macrophage['dataset'].eq('GSE174554'), 'delta'] > 0).sum())}/18上升；GSE274546 {int((macrophage.loc[macrophage['dataset'].eq('GSE274546'), 'delta'] > 0).sum())}/39上升。",
        "- GSE278456从21个原始H5重算真正Miller raw20，21/21均覆盖19/20，只缺AC253572.2；CCL3与CCL4均进入。",
        "- author-marker-derived state赋值前删除全部raw20基因，状态效应只在其余八个预设髓系状态作rest。",
        "",
        "## GSE278456 state effects",
        "",
        state.to_markdown(index=False),
    ]
    path = WRITE_ROOT / "FINAL_RESULT.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> None:
    configure_style()
    WRITE_ROOT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_ROOT.mkdir(parents=True, exist_ok=True)

    if not H5_ROOT.exists() or not MARKER_BOOK.exists():
        raise FileNotFoundError("GSE278456原始H5或作者marker工作簿不可用")

    markers, marker_audit = read_author_markers()
    results = process_gse278456(markers)

    marker_audit.to_csv(WRITE_ROOT / "GSE278456_author_marker_audit.csv", index=False)
    for name, frame in results.items():
        if name == "cells":
            frame.to_csv(
                WRITE_ROOT / "GSE278456_cell_annotations_true_raw20.csv.gz",
                index=False,
                compression={"method": "gzip", "mtime": 0},
            )
        else:
            frame.to_csv(WRITE_ROOT / f"GSE278456_{name}.csv", index=False)
    print("STEP48_GSE278456_STATE_RECALCULATION_COMPLETE")
    print(f"State tables: {WRITE_ROOT}")


if __name__ == "__main__":
    main()
