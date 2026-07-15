#!/usr/bin/env python3
"""生成 Figure 2 两面板正式候选版。

2A：9 个 author-marker-derived 髓系状态的患者级 Miller-IM score state-vs-rest waterfall。
2B：19 个 held-out Miller-IM genes 与 9 套 author marker-set score 的双块热图。

下块 marker-set score 必须从 21 个原始 H5 重算：状态赋值时依然让工作簿全部
17 套 marker score 参与竞争，然后才截取 9 个预设髓系状态做展示。
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import scanpy as sc
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
STEP48_SCRIPT = ROOT / "scripts/48_fig2_identity_state_rebuild.py"
STEP48 = ROOT / "write/48_figure2_identity_state_rebuild/Figure2"
WRITE_ROOT = ROOT / "write/50_figure2_two_panel_final/Figure2"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_ROOT = ROOT / "figures/50_figure2_two_panel_final/Figure2"

EFFECT_FILE = STEP48 / "GSE278456_patient_state_effects.csv"
SUMMARY_FILE = STEP48 / "GSE278456_state_effect_summary.csv"
GENE_SUMMARY_FILE = STEP48 / "GSE278456_gene_summary.csv"
CELL_FILE = STEP48 / "GSE278456_cell_annotations_true_raw20.csv.gz"
MARKER_AUDIT_FILE = STEP48 / "GSE278456_author_marker_audit.csv"

RAW20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
MEASURABLE_RAW20 = [gene for gene in RAW20 if gene != "AC253572.2"]
STATE_ORDER = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]

BLUE = "#3C5488"
RED = "#E64B35"
TEAL = "#00A087"
CYAN = "#4DBBD5"
SALMON = "#F39B7F"
GREY_BLUE = "#8491B4"
MINT = "#91D1C2"
BROWN = "#7E6148"
BLACK = "#222222"
GREY = "#868686"
LIGHT_GREY = "#D8D8D8"
VERY_LIGHT = "#F5F5F5"
POSITIVE = RED
NEGATIVE = BLUE
STATE_COLORS = {
    "MCG1": TEAL,
    "MCG2": MINT,
    "MCG3": CYAN,
    "MCG4": GREY_BLUE,
    "MCG5": BLUE,
    "MAC1": SALMON,
    "MAC2": RED,
    "E-MDSC": "#B09C85",
    "M-MDSC": BROWN,
}
RAW20_CMAP = LinearSegmentedColormap.from_list("raw20_npg", [BLUE, "#F7F7F7", SALMON])
QC_CMAP = LinearSegmentedColormap.from_list("marker_qc", ["#ECECEC", "#A9B3BE", "#4B5563"])
MM = 1 / 25.4


def load_step48_module():
    spec = importlib.util.spec_from_file_location("step48_fig2", STEP48_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Step48 模块无法加载")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.2,
            "axes.labelsize": 7.6,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.3,
            "legend.fontsize": 6.2,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    sc.settings.verbosity = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="生成 Figure2 两面板正式候选版")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="只从本 Step 已写出的冻结源数据重画，不重读 H5",
    )
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_panel(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    pdf = FIG_ROOT / f"{stem}.pdf"
    png = FIG_ROOT / f"{stem}.png"
    metadata = {"Creator": "Step50 Figure2", "CreationDate": None, "ModDate": None}
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.035, metadata=metadata)
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    return pdf, png


def write_source(frame: pd.DataFrame, filename: str) -> Path:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    path = SOURCE_OUT / filename
    frame.to_csv(path, index=False)
    return path


def format_fdr(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}"
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.3f}"


def validate_frozen_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    effects = pd.read_csv(EFFECT_FILE)
    summary = pd.read_csv(SUMMARY_FILE)
    gene_summary = pd.read_csv(GENE_SUMMARY_FILE)
    marker_audit = pd.read_csv(MARKER_AUDIT_FILE)

    if effects.shape[0] != 108 or effects["state"].nunique() != 9:
        raise RuntimeError("Step48 患者状态效应表不再是 108 行 / 9 状态")
    if set(summary["state"]) != set(STATE_ORDER):
        raise RuntimeError("Step48 状态全集发生变化")
    if gene_summary.shape[0] != 171:
        raise RuntimeError("Step48 raw20 基因状态表不再是 171 行")
    if set(gene_summary["gene"]) != set(MEASURABLE_RAW20):
        raise RuntimeError("Step48 可测 raw20 基因集发生变化")
    if gene_summary.groupby("gene")["state"].nunique().ne(9).any():
        raise RuntimeError("raw20 基因未完整覆盖 9 状态")

    focus_audit = marker_audit.loc[marker_audit["state"].isin(STATE_ORDER)].copy()
    if focus_audit["state"].nunique() != 9:
        raise RuntimeError("marker audit 不完整")
    used = set()
    for value in focus_audit["markers_used"].fillna(""):
        used.update(gene for gene in str(value).split(";") if gene)
    overlap = sorted(used & set(MEASURABLE_RAW20))
    if overlap:
        raise RuntimeError(f"赋值后 marker 与 raw20 仍重叠: {overlap}")
    return effects, summary, gene_summary, marker_audit


def recompute_marker_assignment_qc() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    step48 = load_step48_module()
    markers, _ = step48.read_author_markers()
    if len(markers) != 17:
        raise RuntimeError(f"工作簿应提供 17 套状态 marker，实际 {len(markers)}")

    manifest_all = pd.read_csv(step48.STEP30 / "GSE278456_sample_manifest.csv")
    manifest_all["geo_accession"] = manifest_all["geo_accession"].astype(str)
    existing_qc = pd.read_csv(step48.STEP30 / "GSE278456_single_cell_qc.csv").set_index("geo_accession")
    paths = sorted(step48.H5_ROOT.glob("*.h5"))
    if len(paths) != 21:
        raise RuntimeError(f"GSE278456 H5 数量应为 21，实际 {len(paths)}")

    saved_cells = pd.read_csv(
        CELL_FILE,
        usecols=["cell_id", "geo_accession", "author_marker_state", "state_score"],
    )
    if saved_cells.shape[0] != 120766 or saved_cells["cell_id"].duplicated().any():
        raise RuntimeError("Step48 冻结 cell annotation 不再是 120,766 个唯一细胞")
    saved_cells = saved_cells.set_index("cell_id")

    all_marker_genes = set().union(*markers.values())
    patient_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    for file_index, path in enumerate(paths, start=1):
        accession = path.name.split("_")[0]
        manifest_row = manifest_all.loc[manifest_all["geo_accession"].eq(accession)]
        if len(manifest_row) != 1:
            raise RuntimeError(f"{accession} 无法唯一匹配 manifest")
        manifest_row = manifest_row.iloc[0]
        if not bool(manifest_row["is_idhwt_gbm"]):
            raise RuntimeError(f"锁定 H5 出现非 IDH-wt GBM: {accession}")
        sample = str(manifest_row["sample_title"])
        print(f"[marker QC {file_index:02d}/21] {accession} {sample}", flush=True)

        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        n_input = adata.n_obs
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].between(200, 8000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        adata = adata[keep].copy()
        n_qc = adata.n_obs
        old = existing_qc.loc[accession]
        if int(old["n_input"]) != n_input or int(old["n_qc_pass"]) != n_qc:
            raise RuntimeError(f"{accession} QC 细胞数未复现")

        present_raw20 = [gene for gene in RAW20 if gene in adata.var_names]
        if present_raw20 != MEASURABLE_RAW20:
            raise RuntimeError(f"{accession} raw20 覆盖或顺序异常")
        adata.obs_names = [f"{accession}:{barcode}" for barcode in adata.obs_names.astype(str)]
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)

        excluded_controls = set(RAW20) | all_marker_genes
        neutral_control_pool = [gene for gene in adata.var_names if gene not in excluded_controls]
        if len(neutral_control_pool) < 5000:
            raise RuntimeError(f"{accession} state-neutral control pool 过小")

        score_columns: list[str] = []
        for state, genes in markers.items():
            use = [gene for gene in genes if gene in adata.var_names]
            if len(use) < 30:
                raise RuntimeError(f"{accession} {state} 可测 marker 少于 30")
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
        winners = state_names[winner_index]
        winning_score = score_matrix[np.arange(adata.n_obs), winner_index]

        frozen = saved_cells.loc[adata.obs_names]
        winner_matches = frozen["author_marker_state"].to_numpy(object) == winners
        score_difference = np.abs(frozen["state_score"].to_numpy(float) - winning_score)
        if not bool(winner_matches.all()) or float(score_difference.max()) > 1e-10:
            raise RuntimeError(f"{accession} 未精确复现 Step48 状态赋值")
        audit_rows.append(
            {
                "geo_accession": accession,
                "sample": sample,
                "n_input": n_input,
                "n_qc_pass": n_qc,
                "n_winner_exact_match": int(winner_matches.sum()),
                "winner_match_fraction": float(winner_matches.mean()),
                "max_abs_winner_score_difference": float(score_difference.max()),
                "n_competing_marker_sets": len(markers),
            }
        )

        focus = np.isin(winners, STATE_ORDER)
        focus_scores = {
            state: adata.obs[f"state_score_{state}"].to_numpy(float)
            for state in STATE_ORDER
        }
        for assigned_state in STATE_ORDER:
            target = focus & (winners == assigned_state)
            n_cells = int(target.sum())
            if n_cells < 30:
                continue
            for tested_marker_set in STATE_ORDER:
                patient_rows.append(
                    {
                        "geo_accession": accession,
                        "sample": sample,
                        "assigned_state": assigned_state,
                        "tested_marker_set": tested_marker_set,
                        "n_cells": n_cells,
                        "mean_marker_set_score": float(focus_scores[tested_marker_set][target].mean()),
                    }
                )
        del adata, score_matrix

    patient_scores = pd.DataFrame(patient_rows)
    audit = pd.DataFrame(audit_rows).sort_values("geo_accession")
    if int(audit["n_qc_pass"].sum()) != 120766 or not audit["winner_match_fraction"].eq(1).all():
        raise RuntimeError("21 个 H5 的 Step48 赋值复现未全部通过")
    if patient_scores.shape[0] != 972:
        raise RuntimeError(f"患者等权 marker score 长表应为 972 行，实际 {patient_scores.shape[0]}")

    summary = (
        patient_scores.groupby(["tested_marker_set", "assigned_state"], observed=True)
        .agg(
            n_patients=("sample", "nunique"),
            patient_balanced_mean_marker_set_score=("mean_marker_set_score", "mean"),
        )
        .reset_index()
    )
    summary["within_marker_set_z"] = summary.groupby("tested_marker_set", observed=True)[
        "patient_balanced_mean_marker_set_score"
    ].transform(lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0)
    if summary.shape[0] != 81:
        raise RuntimeError("marker-set QC 汇总表不是 9×9")
    return patient_scores, summary, audit


def prepare_waterfall_source(effects: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    summary_index = summary.set_index("state")
    rows: list[pd.DataFrame] = []
    for state in STATE_ORDER:
        frame = effects.loc[effects["state"].eq(state)].sort_values("raw20_delta").copy()
        frame["within_state_rank"] = np.arange(1, len(frame) + 1)
        frame["direction"] = np.where(frame["raw20_delta"] > 0, "Above rest", "Below rest")
        frame["n_patients"] = int(summary_index.loc[state, "n_patients"])
        frame["n_positive"] = int(summary_index.loc[state, "n_positive"])
        frame["mean_delta"] = float(summary_index.loc[state, "mean_delta"])
        frame["fdr_9_states"] = float(summary_index.loc[state, "fdr_9_states"])
        frame["score_pseudobulk_mean_direction_concordant"] = bool(
            summary_index.loc[state, "score_pseudobulk_mean_direction_concordant"]
        )
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def make_patient_waterfall(source_frame: pd.DataFrame) -> dict[str, str]:
    stem = "Figure2_state_patient_waterfall"
    source = write_source(source_frame, f"{stem}_source.csv")
    summary = source_frame.drop_duplicates("state").set_index("state")

    fig, axes = plt.subplots(3, 3, figsize=(183 * MM, 104 * MM), sharex=True, sharey=True)
    y_min, y_max = -0.76, 0.56
    for index, (ax, state) in enumerate(zip(axes.flat, STATE_ORDER, strict=True)):
        frame = source_frame.loc[source_frame["state"].eq(state)]
        x = frame["within_state_rank"].to_numpy(int)
        y = frame["raw20_delta"].to_numpy(float)
        colors = np.where(y > 0, POSITIVE, NEGATIVE)
        ax.axhspan(0, y_max, color="#FFF7F4", zorder=0)
        ax.axhspan(y_min, 0, color="#F4F6FA", zorder=0)
        ax.axhline(0, color=BLACK, linewidth=0.65, zorder=2)
        ax.bar(x, y, width=0.76, color=colors, edgecolor="none", alpha=0.88, zorder=3)
        mean = float(summary.loc[state, "mean_delta"])
        ax.axhline(mean, color=STATE_COLORS[state], linewidth=1.05, linestyle=(0, (3, 2)), zorder=4)
        ax.set_xlim(0.25, 21.75)
        ax.set_ylim(y_min, y_max)
        state_title = f"{state}{'†' if int(summary.loc[state, 'n_patients']) < 6 else ''}"
        ax.set_title(state_title, loc="left", color=STATE_COLORS[state], fontweight="bold", pad=7)
        ax.text(
            1.0,
            1.025,
            f"{int(summary.loc[state, 'n_positive'])}/{int(summary.loc[state, 'n_patients'])}  |  FDR {format_fdr(float(summary.loc[state, 'fdr_9_states']))}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.8,
            color=BLACK,
        )
        if state == "MCG3":
            for spine in ax.spines.values():
                spine.set_color(GREY)
                spine.set_linestyle((0, (3, 2)))
            ax.text(
                0.98,
                0.06,
                "score / pseudobulk\ndirection differs",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=5.3,
                color=GREY,
            )
        else:
            ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(BLACK)
        ax.tick_params(length=2.2, width=0.65)
        if index % 3 != 0:
            ax.tick_params(axis="y", left=False, labelleft=False)
        if index % 3 == 0:
            ax.set_ylabel("State − other myeloid\nMiller-IM score")
        if index >= 6:
            ax.set_xticks([1, 7, 14, 21])
            ax.set_xlabel("Patients ranked within state")
        else:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
    fig.text(
        0.995,
        0.012,
        "† fewer than 6 patients; all panels use 21 fixed patient slots",
        ha="right",
        va="bottom",
        fontsize=5.6,
        color=GREY,
    )
    fig.subplots_adjust(left=0.09, right=0.995, bottom=0.13, top=0.965, wspace=0.13, hspace=0.29)
    pdf, png = save_panel(fig, stem)
    return {"panel": "2A", "source": str(source), "pdf": str(pdf), "png": str(png)}


def make_double_heatmap(
    gene_summary: pd.DataFrame,
    marker_summary: pd.DataFrame,
    state_summary: pd.DataFrame,
) -> dict[str, str]:
    stem = "Figure2_raw20_marker_double_heatmap"
    gene_source = gene_summary.copy()
    gene_source["block"] = "Held-out Miller-IM genes; not used for state assignment"
    marker_source = marker_summary.copy()
    marker_source["block"] = "Author marker-set assignment scores; classification QC"
    source_top = write_source(gene_source, f"{stem}_raw20_source.csv")
    source_bottom = write_source(marker_source, f"{stem}_marker_qc_source.csv")

    top = (
        gene_summary.pivot(index="gene", columns="state", values="within_gene_z")
        .reindex(index=MEASURABLE_RAW20, columns=STATE_ORDER)
    )
    bottom = (
        marker_summary.pivot(index="tested_marker_set", columns="assigned_state", values="within_marker_set_z")
        .reindex(index=STATE_ORDER, columns=STATE_ORDER)
    )
    if top.isna().any().any() or bottom.isna().any().any():
        raise RuntimeError("双块热图存在缺失单元格")

    n_by_state = state_summary.set_index("state")["n_patients"].astype(int).to_dict()
    state_labels = [
        f"{state}\n(n={n_by_state[state]}){'†' if n_by_state[state] < 6 else ''}"
        for state in STATE_ORDER
    ]

    fig = plt.figure(figsize=(183 * MM, 148 * MM))
    grid = fig.add_gridspec(
        nrows=4,
        ncols=2,
        width_ratios=[1, 0.034],
        height_ratios=[2.2, 19, 1.35, 5.3],
        left=0.16,
        right=0.94,
        bottom=0.10,
        top=0.98,
        wspace=0.16,
        hspace=0.05,
    )
    title_top = fig.add_subplot(grid[0, 0])
    ax_top = fig.add_subplot(grid[1, 0])
    separator = fig.add_subplot(grid[2, 0])
    ax_bottom = fig.add_subplot(grid[3, 0])
    cax_top = fig.add_subplot(grid[1, 1])
    cax_bottom = fig.add_subplot(grid[3, 1])

    title_top.axis("off")
    title_top.text(
        0,
        0.93,
        "Held-out Miller-IM genes — not used for state assignment",
        ha="left",
        va="center",
        fontsize=7.4,
        fontweight="bold",
        color=BLACK,
    )
    sns.heatmap(
        top,
        ax=ax_top,
        cmap=RAW20_CMAP,
        vmin=-2.5,
        vmax=2.5,
        center=0,
        linewidths=0.45,
        linecolor="white",
        cbar=True,
        cbar_ax=cax_top,
        xticklabels=state_labels,
        yticklabels=MEASURABLE_RAW20,
    )
    ax_top.xaxis.tick_top()
    ax_top.tick_params(axis="x", rotation=0, length=0, pad=2)
    ax_top.tick_params(axis="y", rotation=0, length=0, pad=2)
    ax_top.set_xlabel("")
    ax_top.set_ylabel("")
    cax_top.set_ylabel("Within-gene z score", fontsize=6.3, labelpad=4)
    cax_top.tick_params(labelsize=5.7, length=2)

    separator.axis("off")
    separator.axhline(0.92, color=BLACK, linewidth=0.8)
    separator.text(
        0,
        0.55,
        "Author marker-set assignment scores — used to define labels; classification QC",
        ha="left",
        va="center",
        fontsize=6.8,
        fontweight="bold",
        color=GREY,
    )
    separator.text(
        0,
        0.08,
        "Direct marker–Miller-IM gene overlap after exclusion: 0",
        ha="left",
        va="center",
        fontsize=5.9,
        color=GREY,
    )
    sns.heatmap(
        bottom,
        ax=ax_bottom,
        cmap=QC_CMAP,
        vmin=-2.5,
        vmax=2.5,
        center=0,
        linewidths=0.45,
        linecolor="white",
        cbar=True,
        cbar_ax=cax_bottom,
        xticklabels=False,
        yticklabels=[f"{state} markers" for state in STATE_ORDER],
    )
    ax_bottom.tick_params(axis="y", rotation=0, length=0, pad=2, colors=GREY)
    ax_bottom.tick_params(axis="x", bottom=False)
    ax_bottom.set_xlabel("")
    ax_bottom.set_ylabel("")
    cax_bottom.set_ylabel("Within-marker-set z score", fontsize=5.8, labelpad=4, color=GREY)
    cax_bottom.tick_params(labelsize=5.5, length=2, colors=GREY)
    for spine in cax_bottom.spines.values():
        spine.set_edgecolor(GREY)

    fig.text(
        0.16,
        0.035,
        "Columns are fixed a priori; † fewer than 6 patients. The lower block documents label construction and is not independent validation.",
        ha="left",
        va="bottom",
        fontsize=5.7,
        color=GREY,
    )
    pdf, png = save_panel(fig, stem)
    return {
        "panel": "2B",
        "source": f"{source_top};{source_bottom}",
        "pdf": str(pdf),
        "png": str(png),
    }


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path]:
    images = [Image.open(record["png"]).convert("RGB") for record in records]
    target_width = 1800
    resized: list[Image.Image] = []
    for image in images:
        height = round(image.height * target_width / image.width)
        resized.append(image.resize((target_width, height), Image.Resampling.LANCZOS))

    margin = 85
    gap = 40
    canvas = Image.new(
        "RGB",
        (target_width + margin * 2, sum(image.height for image in resized) + gap + margin * 2),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    font_candidates = [
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    font_path = next((Path(path) for path in font_candidates if Path(path).exists()), None)
    font = ImageFont.truetype(str(font_path), 54) if font_path else ImageFont.load_default()
    y = margin
    for label, image in zip(["A", "B"], resized, strict=True):
        draw.text((22, y + 2), label, fill=BLACK, font=font)
        canvas.paste(image, (margin, y))
        y += image.height + gap
    preview = FIG_ROOT / "Figure2_two_panel_preview.png"
    canvas.save(preview, dpi=(150, 150))
    grayscale = FIG_ROOT / "Figure2_two_panel_preview_grayscale.png"
    canvas.convert("L").convert("RGB").save(grayscale, dpi=(150, 150))
    return preview, grayscale


def write_legend() -> Path:
    path = WRITE_ROOT / "Figure2_legend.md"
    path.write_text(
        """# Figure 2 | The Miller-IM program preferentially maps to MCG1/MCG2-like states in an independent cross-sectional myeloid atlas

**A.** Patient-level state-versus-rest waterfall plots in the 21 IDH-wild-type primary GBM samples of GSE278456. Each bar is one patient with at least 30 cells in the indicated author-marker-derived state and at least 30 cells in the remaining prespecified myeloid states. Bars show the difference in control-adjusted Miller-IM score between that state and the remaining states; dashed lines show patient-equal means. Red and blue denote values above and below the remaining states, respectively. Numbers give patients above the remaining states / evaluable patients and Benjamini–Hochberg FDR across the nine state-wise exact two-sided sign-flip tests. All facets use the same 21 patient slots and y-axis. MCG3 is marked because its mean score direction differs from the raw-count pseudobulk sensitivity analysis.

**B.** Two-block patient-balanced heatmap across the same fixed state order. The upper block shows relative expression of the 19 measurable Miller-IM genes; values are patient-equal mean log-normalized expression standardized within each gene across the nine states. All measurable Miller-IM genes were removed from the author marker sets before state assignment, so this block was held out from direct state-label construction. The lower, visually reduced block shows the nine complete author marker-set scores, averaged within patient and assigned state and then standardized within marker set. These scores were used to create the state labels; their diagonal structure is expected by construction and is shown only as classification quality control, not as independent validation. Direct overlap between the final label-defining marker sets and the measurable Miller-IM genes was zero. State labels show the number of evaluable patients; dagger marks fewer than six patients.

GSE278456 is cross-sectional and is used only to localize the recurrence-associated program across myeloid states; recurrence evidence is provided separately in Figure 1.
""",
        encoding="utf-8",
    )
    return path


def write_summary(
    marker_patient: pd.DataFrame,
    marker_summary: pd.DataFrame,
    audit: pd.DataFrame,
    records: list[dict[str, str]],
    preview: Path,
    grayscale: Path,
) -> Path:
    summary = WRITE_ROOT / "FINAL_RESULT.md"
    summary.write_text(
        f"""# Step50 Figure2 两面板候选版

- 2A：9 状态患者级 waterfall，固定 21 个患者槽位和统一 y 轴。
- 2B 上块：19×9 held-out Miller-IM genes 患者等权 `within_gene_z`。
- 2B 下块：9×9 完整 author marker-set score QC，先在 17 套 marker 间完成赋值，再截取 9 个髓系状态。
- H5 复核：{int(audit['n_qc_pass'].sum()):,}/{int(audit['n_qc_pass'].sum()):,} 个 QC 细胞的 winner label 与 Step48 精确一致；最大 winner-score 差 {audit['max_abs_winner_score_difference'].max():.3g}。
- 源数据：waterfall 108 行；Miller-IM 上块 171 行；marker 患者长表 {marker_patient.shape[0]} 行；marker QC 汇总 {marker_summary.shape[0]} 行。
- 预览：`{preview}`；灰度预览：`{grayscale}`。
- 边界：GSE278456 只做横断面状态定位；下块只是标签构建 QC，不是独立验证。

## 面板

{pd.DataFrame(records).to_markdown(index=False)}
""",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    args = parse_args()
    configure_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    WRITE_ROOT.mkdir(parents=True, exist_ok=True)

    effects, state_summary, gene_summary, _ = validate_frozen_inputs()
    waterfall = prepare_waterfall_source(effects, state_summary)

    marker_patient_file = SOURCE_OUT / "GSE278456_patient_marker_set_scores_source.csv"
    marker_summary_file = SOURCE_OUT / "GSE278456_marker_set_assignment_qc_source.csv"
    marker_audit_file = SOURCE_OUT / "GSE278456_step48_assignment_reproduction_source.csv"
    if args.plot_only:
        required = [marker_patient_file, marker_summary_file, marker_audit_file]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise RuntimeError(f"--plot-only 缺少冻结源数据: {missing}")
        marker_patient = pd.read_csv(marker_patient_file)
        marker_summary = pd.read_csv(marker_summary_file)
        audit = pd.read_csv(marker_audit_file)
    else:
        marker_patient, marker_summary, audit = recompute_marker_assignment_qc()
        marker_patient.to_csv(marker_patient_file, index=False)
        marker_summary.to_csv(marker_summary_file, index=False)
        audit.to_csv(marker_audit_file, index=False)

    if marker_patient.shape[0] != 972 or marker_summary.shape[0] != 81:
        raise RuntimeError("冻结 marker-set 源数据尺寸异常")
    records = [
        make_patient_waterfall(waterfall),
        make_double_heatmap(gene_summary, marker_summary, state_summary),
    ]
    preview, grayscale = make_preview(records)
    legend = write_legend()
    manifest = WRITE_ROOT / "Figure2_panel_manifest.csv"
    pd.DataFrame(records).to_csv(manifest, index=False)
    summary = write_summary(marker_patient, marker_summary, audit, records, preview, grayscale)

    inventory = []
    for path in sorted(list(FIG_ROOT.glob("*.pdf")) + list(FIG_ROOT.glob("*.png"))):
        inventory.append({"path": str(path), "sha256": file_sha256(path)})
    pd.DataFrame(inventory).to_csv(WRITE_ROOT / "Figure2_output_hashes.csv", index=False)
    print(f"Wrote {preview}")
    print(f"Wrote {legend}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
