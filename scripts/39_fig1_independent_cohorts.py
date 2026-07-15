#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
WRITE_ROOT = ROOT / "write" / "39_independent_miller_mg_inflammatory_figures"
FIG_ROOT = ROOT / "figures" / "39_independent_miller_mg_inflammatory_plottie"
FIG1_DATA = WRITE_ROOT / "Figure1" / "source_data"
FIG1_OUT = FIG_ROOT / "Figure1" / "panel_library"
S1_DATA = WRITE_ROOT / "SupplementaryFigure1" / "source_data"
S1_OUT = FIG_ROOT / "SupplementaryFigure1" / "panel_library"
S2_DATA = WRITE_ROOT / "SupplementaryFigure2" / "source_data"
S2_OUT = FIG_ROOT / "SupplementaryFigure2" / "panel_library"

H5ADS = {
    "GSE174554": ROOT / "write/34_gse174554_raw_independent_discovery/03_myeloid/GSE174554_raw_independent_clean_myeloid.h5ad",
    "GSE274546": ROOT / "write/36_gse274546_raw_independent_reannotation/03_myeloid/GSE274546_raw_independent_clean_myeloid.h5ad",
}
STEP_DIRS = {
    "GSE174554": ROOT / "write/34_gse174554_raw_independent_discovery",
    "GSE274546": ROOT / "write/36_gse274546_raw_independent_reannotation",
}
STEP38 = ROOT / "write/38_independent_cohort_mg_inflammatory_recalculation"

RAW20 = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]
MARKERS = {
    "Microglia": ["P2RY12", "TMEM119", "SALL1", "CX3CR1"],
    "Macrophage": ["GPNMB", "SPP1", "MRC1", "MSR1"],
    "Monocyte": ["VCAN", "FCN1", "S100A8", "S100A9"],
    "CDC": ["FCER1A", "CD1C", "CLEC10A", "HLA-DPB1"],
}
DATASET_COLORS = {"GSE174554": "#3C5488", "GSE274546": "#E64B35"}
IDENTITY_COLORS = {
    "Macrophage": "#3C5488",
    "Microglia": "#00A087",
    "Monocyte": "#F39B7F",
    "CDC": "#8491B4",
    "Neutrophil": "#7E6148",
}
PRIMARY = "#4DBBD5"
RECURRENT = "#E64B35"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })
    sns.set_style("white")


def save_panel(fig: plt.Figure, out_dir: Path, name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def dense_mean_and_pct(x) -> tuple[np.ndarray, np.ndarray]:
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel()
        pct = np.asarray((x > 0).mean(axis=0)).ravel() * 100
    else:
        mean = np.asarray(x).mean(axis=0)
        pct = (np.asarray(x) > 0).mean(axis=0) * 100
    return mean, pct


def panel_workflow() -> None:
    fig, ax = plt.subplots(figsize=(8.1, 3.1))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    def box(x, y, w, h, text, edge, fill="white", lw=1.2, weight="normal"):
        rect = mpl.patches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
            facecolor=fill, edgecolor=edge, linewidth=lw,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8, weight=weight)

    box(0.3, 4.8, 5.1, 1.0, "GSE174554  |  91 matrices  |  snRNA-seq\nPrimary -> 1st recurrent", DATASET_COLORS["GSE174554"], fill="#F3F6FA", weight="bold")
    box(6.6, 4.8, 5.1, 1.0, "GSE274546  |  111 libraries  |  scRNA-seq\nPrimary -> 1st recurrent", DATASET_COLORS["GSE274546"], fill="#FDF2EF", weight="bold")
    for x in (2.85, 9.15):
        ax.annotate("", xy=(x, 4.15), xytext=(x, 4.75), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0))
    box(0.65, 3.15, 4.4, 0.82, "Independent QC, doublet removal,\nHVG/PCA/Harmony/Leiden", "#666666")
    box(6.95, 3.15, 4.4, 0.82, "Independent QC, doublet removal,\nHVG/PCA/Harmony/Leiden", "#666666")
    for x in (2.85, 9.15):
        ax.annotate("", xy=(x, 2.48), xytext=(x, 3.12), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0))
    box(0.65, 1.45, 4.4, 0.85, "23,344 clean myeloid\n18 paired patients (>=20 cells/side)", DATASET_COLORS["GSE174554"], fill="#F3F6FA")
    box(6.95, 1.45, 4.4, 0.85, "68,109 clean myeloid\n45 paired patients (>=20 cells/side)", DATASET_COLORS["GSE274546"], fill="#FDF2EF")
    for x in (2.85, 9.15):
        ax.annotate("", xy=(6.0, 0.75), xytext=(x, 1.4), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0))
    box(3.75, 0.05, 4.5, 0.68, "Paired raw-count pseudobulk\nedgeR + fixed-program GSEA", "#00A087", fill="#F0FAF7", weight="bold")
    ax.text(6, -0.18, "Patient is the statistical replicate", ha="center", va="top", color="#555555", fontsize=7)
    source = pd.DataFrame({
        "dataset": ["GSE174554", "GSE274546"],
        "input_libraries_or_matrices": [91, 111],
        "clean_myeloid": [23344, 68109],
        "formal_pairs_threshold20": [18, 45],
    })
    source.to_csv(FIG1_DATA / "cohort_workflow.csv", index=False)
    save_panel(fig, FIG1_OUT, "cohort_independent_processing")


def panel_umap(dataset: str, adata: ad.AnnData) -> None:
    coords = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs["myeloid_identity"].astype(str).to_numpy()
    order = np.argsort(labels)
    coords = coords[order]
    labels = labels[order]
    fig, ax = plt.subplots(figsize=(4.0, 3.8))
    for identity in sorted(np.unique(labels)):
        mask = labels == identity
        ax.scatter(coords[mask, 0], coords[mask, 1], s=1.0, alpha=0.7, linewidths=0,
                   color=IDENTITY_COLORS.get(identity, "#BDBDBD"), rasterized=True, label=identity)
        center = np.median(coords[mask], axis=0)
        ax.text(center[0], center[1], identity, ha="center", va="center", fontsize=7,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.8))
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), markerscale=4, handletextpad=0.2)
    pd.DataFrame({
        "cell_id": adata.obs_names.to_numpy(),
        "UMAP1": np.asarray(adata.obsm["X_umap"])[:, 0],
        "UMAP2": np.asarray(adata.obsm["X_umap"])[:, 1],
        "myeloid_identity": adata.obs["myeloid_identity"].astype(str).to_numpy(),
    }).to_csv(FIG1_DATA / f"umap_identity_{dataset}.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    save_panel(fig, FIG1_OUT, f"UMAP_identity_{dataset}")


def build_dotplot(adatas: dict[str, ad.AnnData]) -> pd.DataFrame:
    genes = [g for group in MARKERS.values() for g in group]
    rows = []
    for dataset, adata in adatas.items():
        for identity in sorted(adata.obs["myeloid_identity"].astype(str).unique()):
            mask = adata.obs["myeloid_identity"].astype(str).eq(identity).to_numpy()
            present = [g for g in genes if g in adata.var_names]
            mean, pct = dense_mean_and_pct(adata[mask, present].X)
            for gene, value, fraction in zip(present, mean, pct):
                rows.append({"dataset": dataset, "identity": identity, "gene": gene,
                             "mean_expression": value, "pct_positive": fraction})
    out = pd.DataFrame(rows)
    out["scaled_expression"] = out.groupby(["dataset", "gene"])["mean_expression"].transform(
        lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1)
    )
    out.to_csv(FIG1_DATA / "identity_marker_dotplot.csv", index=False)
    return out


def panel_dotplot(dot: pd.DataFrame) -> None:
    genes = [g for group in MARKERS.values() for g in group]
    identities = [x for x in ["Microglia", "Macrophage", "Monocyte", "CDC"] if x in set(dot["identity"])]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), sharey=True, gridspec_kw={"wspace": 0.06})
    for ax, dataset in zip(axes, ["GSE174554", "GSE274546"]):
        d = dot[dot["dataset"].eq(dataset)].copy()
        d["x"] = d["gene"].map({g: i for i, g in enumerate(genes)})
        d["y"] = d["identity"].map({g: i for i, g in enumerate(identities)})
        ax.scatter(d["x"], d["y"], s=8 + d["pct_positive"] * 0.8,
                   c=d["scaled_expression"], cmap="coolwarm", vmin=-2, vmax=2,
                   edgecolors="none")
        ax.set_xlim(-0.6, len(genes) - 0.4)
        ax.set_ylim(len(identities) - 0.4, -0.6)
        ax.set_xticks(range(len(genes)), genes, rotation=55, ha="right")
        ax.set_yticks(range(len(identities)), identities)
        ax.grid(color="#ECECEC", lw=0.5)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        ax.text(0.5, 1.03, dataset, transform=ax.transAxes, ha="center", va="bottom", weight="bold")
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-2, 2), cmap="coolwarm")
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Scaled expression")
    save_panel(fig, FIG1_OUT, "identity_marker_dotplot")


def enrichment_curve(deg: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    d = deg.copy()
    d["rank_stat"] = np.sign(d["logFC"]) * np.sqrt(np.maximum(d["F"], 0))
    d = d[np.isfinite(d["rank_stat"])].sort_values("rank_stat", ascending=False).reset_index(drop=True)
    hits = d["gene"].isin(genes).to_numpy()
    weights = np.abs(d["rank_stat"].to_numpy())
    hit_weights = np.where(hits, weights, 0)
    hit_norm = hit_weights.sum()
    miss_norm = (~hits).sum()
    increments = np.where(hits, hit_weights / hit_norm, -1 / miss_norm)
    running = np.cumsum(increments)
    return pd.DataFrame({
        "rank": np.arange(1, len(d) + 1), "gene": d["gene"], "rank_stat": d["rank_stat"],
        "hit": hits, "running_enrichment": running,
    })


def panel_gsea_curves() -> None:
    deg = pd.read_csv(STEP38 / "independent_paired_edger_all_genes.csv")
    gsea = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 5.6), sharex=False)
    sources = []
    for ax, dataset in zip(axes, ["GSE174554", "GSE274546"]):
        current = deg[(deg["dataset"].eq(dataset)) & (deg["threshold"].eq(20))]
        curve = enrichment_curve(current, RAW20)
        curve["dataset"] = dataset
        sources.append(curve)
        color = DATASET_COLORS[dataset]
        ax.plot(curve["rank"], curve["running_enrichment"], color=color, lw=1.6)
        ax.fill_between(curve["rank"], 0, curve["running_enrichment"], color=color, alpha=0.10)
        hit_ranks = curve.loc[curve["hit"], "rank"]
        ax.vlines(hit_ranks, ymin=-0.025, ymax=0.025, color="#555555", lw=0.5)
        row = gsea[(gsea["dataset"].eq(dataset)) & (gsea["threshold"].eq(20)) &
                   (gsea["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20"))].iloc[0]
        ax.text(0.02, 0.90, f"{dataset}   NES {row.NES:.2f}   nominal P {row.pval:.3g}",
                transform=ax.transAxes, ha="left", va="top", fontsize=8, weight="bold")
        ax.axhline(0, color="#999999", lw=0.6)
        ax.set_ylabel("Running enrichment")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Ranked genes")
    pd.concat(sources, ignore_index=True).to_csv(FIG1_DATA / "raw20_gsea_curves.csv.gz", index=False,
                                                 compression={"method": "gzip", "mtime": 0})
    save_panel(fig, FIG1_OUT, "Miller_raw20_GSEA_curves")


def panel_threshold_forest() -> None:
    gsea = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    d = gsea[gsea["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()
    d.to_csv(FIG1_DATA / "raw20_threshold_forest.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    ymap = {"GSE174554": 1, "GSE274546": 0}
    offsets = {20: 0.10, 50: -0.10}
    for _, row in d.iterrows():
        y = ymap[row.dataset] + offsets[int(row.threshold)]
        formal = bool(row.formal_testable)
        ax.scatter(row.NES, y, s=40, color=DATASET_COLORS[row.dataset] if formal else "white",
                   edgecolor=DATASET_COLORS[row.dataset], linewidth=1.1,
                   marker="o" if row.threshold == 20 else "s", zorder=3)
        ax.text(row.NES + 0.025, y, f"{row.NES:.2f}", va="center", fontsize=7)
    ax.axvline(0, color="#999999", lw=0.7)
    ax.set_yticks([0, 1], ["GSE274546", "GSE174554"])
    ax.set_xlabel("Normalized enrichment score")
    ax.set_ylim(-0.45, 1.45)
    ax.spines[["top", "right", "left"]].set_visible(False)
    handles = [
        mpl.lines.Line2D([], [], marker="o", color="#555555", linestyle="", label=">=20 cells/side"),
        mpl.lines.Line2D([], [], marker="s", color="#555555", linestyle="", label=">=50 cells/side"),
        mpl.lines.Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#555555", linestyle="", label="<10 pairs: descriptive"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    save_panel(fig, FIG1_OUT, "raw20_threshold_stability")


def supplementary_qc(adatas: dict[str, ad.AnnData]) -> None:
    summary_rows, capture_rows = [], []
    for dataset, base in STEP_DIRS.items():
        capture = pd.read_csv(base / "01_qc_doublets/capture_qc_doublet_summary.csv")
        capture["dataset"] = dataset
        capture_rows.append(capture)
        summary_rows.append({
            "dataset": dataset,
            "raw": int(capture["n_raw"].sum()),
            "qc": int(capture["n_qc"].sum()),
            "singlet": int(capture.loc[capture["status"].eq("ok"), "n_singlet"].sum()),
            "clean_myeloid": int(adatas[dataset].n_obs),
        })
    summary = pd.DataFrame(summary_rows)
    captures = pd.concat(capture_rows, ignore_index=True)
    summary.to_csv(S1_DATA / "qc_stage_counts.csv", index=False)
    captures.to_csv(S1_DATA / "capture_qc_doublets.csv", index=False)

    long = summary.melt(id_vars="dataset", var_name="stage", value_name="cells")
    fig, ax = plt.subplots(figsize=(4.4, 3.1))
    sns.barplot(data=long, x="stage", y="cells", hue="dataset", palette=DATASET_COLORS, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Cells")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_panel(fig, S1_OUT, "qc_stage_counts")

    fig, ax = plt.subplots(figsize=(4.3, 3.0))
    sns.boxplot(data=captures[captures["status"].eq("ok")], x="dataset", y="doublet_rate",
                hue="dataset", palette=DATASET_COLORS, legend=False, width=0.5, fliersize=0, ax=ax)
    sns.stripplot(data=captures[captures["status"].eq("ok")], x="dataset", y="doublet_rate",
                  color="#333333", size=1.8, alpha=0.45, jitter=0.18, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Doublet rate")
    ax.spines[["top", "right"]].set_visible(False)
    save_panel(fig, S1_OUT, "capture_doublet_rates")

    meta_rows = []
    threshold_rows = []
    identity_rows = []
    review_rows = []
    for dataset, base in STEP_DIRS.items():
        meta = pd.read_csv(base / "03_myeloid/patient_condition_pseudobulk_metadata.csv")
        meta["dataset"] = dataset
        meta_rows.append(meta)
        for threshold in [20, 50, 100, 200]:
            eligible = meta[meta["n_cells"].ge(threshold)]
            counts = eligible.groupby("patient_id")["condition"].nunique()
            threshold_rows.append({"dataset": dataset, "threshold": threshold, "n_pairs": int(counts.ge(2).sum())})
        identities = adatas[dataset].obs["myeloid_identity"].astype(str).value_counts()
        for identity, n in identities.items():
            identity_rows.append({"dataset": dataset, "identity": identity, "n_cells": int(n)})
        review = pd.read_csv(base / "03_myeloid/myeloid_cluster_review.csv")
        review["dataset"] = dataset
        review_rows.append(review)
    meta = pd.concat(meta_rows, ignore_index=True)
    thresholds = pd.DataFrame(threshold_rows)
    identities = pd.DataFrame(identity_rows)
    reviews = pd.concat(review_rows, ignore_index=True)
    meta.to_csv(S1_DATA / "patient_condition_cell_counts.csv", index=False)
    thresholds.to_csv(S1_DATA / "pair_count_thresholds.csv", index=False)
    identities.to_csv(S1_DATA / "identity_cell_counts.csv", index=False)
    reviews.to_csv(S1_DATA / "cluster_contamination_review.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.5, 3.1))
    sns.boxplot(data=meta, x="condition", y="n_cells", hue="dataset", palette=DATASET_COLORS,
                fliersize=0, width=0.65, ax=ax)
    sns.stripplot(data=meta, x="condition", y="n_cells", hue="dataset", dodge=True,
                  palette=DATASET_COLORS, size=1.5, alpha=0.35, ax=ax, legend=False)
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.set_ylabel("Clean myeloid cells per patient endpoint")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_panel(fig, S1_OUT, "patient_endpoint_cell_counts")

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for dataset, d in thresholds.groupby("dataset"):
        ax.plot(d["threshold"], d["n_pairs"], marker="o", lw=1.4, color=DATASET_COLORS[dataset], label=dataset)
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200], [20, 50, 100, 200])
    ax.set_xlabel("Minimum cells per endpoint")
    ax.set_ylabel("Eligible patient pairs")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_panel(fig, S1_OUT, "pair_count_thresholds")

    fig, ax = plt.subplots(figsize=(4.5, 3.1))
    pivot = identities.pivot(index="dataset", columns="identity", values="n_cells").fillna(0)
    pivot = pivot[[c for c in IDENTITY_COLORS if c in pivot.columns]]
    bottom = np.zeros(len(pivot))
    for identity in pivot.columns:
        values = pivot[identity].to_numpy()
        ax.bar(pivot.index, values, bottom=bottom, color=IDENTITY_COLORS[identity], label=identity)
        bottom += values
    ax.set_ylabel("Clean myeloid cells")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_panel(fig, S1_OUT, "identity_cell_counts")

    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    plot = reviews.copy()
    plot["cluster_label"] = plot["dataset"].str.replace("GSE", "") + ":" + plot["cluster"].astype(str)
    colors = plot["is_clean_myeloid"].map({True: "#00A087", False: "#BDBDBD"})
    ax.scatter(plot["top20_contamination_hits"], plot["ribosomal_top20_hits"],
               s=10 + np.sqrt(plot["n_cells"]) * 2, c=colors, alpha=0.75, edgecolors="white", linewidth=0.3)
    ax.axvline(3, color="#888888", ls="--", lw=0.7)
    ax.axhline(6, color="#888888", ls="--", lw=0.7)
    ax.set_xlabel("Contamination markers among top 20")
    ax.set_ylabel("Ribosomal markers among top 20")
    ax.spines[["top", "right"]].set_visible(False)
    save_panel(fig, S1_OUT, "cluster_contamination_audit")


def supplementary_sensitivity() -> None:
    targeted = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    scores = pd.read_csv(STEP38 / "independent_fixed_program_paired_scores.csv")
    genes = pd.read_csv(STEP38 / "independent_fixed_program_gene_direction.csv")
    panel = pd.read_csv(STEP38 / "independent_old_panel_gsea_recomputed.csv")
    raw = targeted[targeted["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()
    raw_scores = scores[scores["signature"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()
    raw_genes = genes[genes["signature"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()
    raw.to_csv(S2_DATA / "raw20_gsea.csv", index=False)
    raw_scores.to_csv(S2_DATA / "raw20_paired_scores.csv", index=False)
    raw_genes.to_csv(S2_DATA / "raw20_gene_directions.csv", index=False)
    panel.to_csv(S2_DATA / "full_panel_gsea.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.8, 3.1))
    sns.pointplot(data=raw, x="threshold", y="NES", hue="dataset", palette=DATASET_COLORS,
                  dodge=0.12, markers="o", linestyles="-", errorbar=None, ax=ax)
    ax.set_xlabel("Minimum cells per endpoint")
    ax.set_ylabel("NES")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_panel(fig, S2_OUT, "raw20_threshold_sensitivity")

    merged = raw_scores.merge(raw, left_on=["dataset", "threshold", "signature"],
                              right_on=["dataset", "threshold", "pathway"], how="inner")
    merged = merged[merged["threshold"].eq(20)]
    fig, ax = plt.subplots(figsize=(4.3, 3.2))
    for _, row in merged.iterrows():
        ax.scatter(row.mean_delta, row.NES, color=DATASET_COLORS[row.dataset], s=50,
                   marker="o" if "top20" in row.signature else "s")
        ax.text(row.mean_delta, row.NES + 0.025, row.dataset.replace("GSE", ""), ha="center", fontsize=6.5)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.set_xlabel("Mean patient score change")
    ax.set_ylabel("GSEA NES")
    ax.spines[["top", "right"]].set_visible(False)
    save_panel(fig, S2_OUT, "average_score_vs_gsea")

    gd = raw_genes[raw_genes["threshold"].eq(20)]
    mat = gd.pivot(index="gene", columns="dataset", values="logFC").reindex(RAW20)
    fig, ax = plt.subplots(figsize=(3.8, 5.4))
    sns.heatmap(mat, cmap="coolwarm", center=0, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Paired logFC", "shrink": 0.55}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    save_panel(fig, S2_OUT, "raw20_gene_logfc")

    hp = panel[panel["threshold"].eq(20)].copy()
    keep = hp.groupby("pathway")["pval"].min().sort_values().head(14).index
    hp = hp[hp["pathway"].isin(keep)]
    mat = hp.pivot(index="pathway", columns="dataset", values="NES").reindex(keep)
    fig, ax = plt.subplots(figsize=(4.2, 5.0))
    sns.heatmap(mat, cmap="coolwarm", center=0, linewidths=0.5, linecolor="white",
                cbar_kws={"label": "NES", "shrink": 0.55}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    save_panel(fig, S2_OUT, "independent_program_gsea_landscape")


def main() -> None:
    np.random.seed(20260713)
    set_style()
    for path in [FIG1_DATA, FIG1_OUT, S1_DATA, S1_OUT, S2_DATA, S2_OUT]:
        path.mkdir(parents=True, exist_ok=True)
    adatas = {dataset: ad.read_h5ad(path) for dataset, path in H5ADS.items()}
    panel_workflow()
    for dataset, adata in adatas.items():
        panel_umap(dataset, adata)
    panel_dotplot(build_dotplot(adatas))
    panel_gsea_curves()
    panel_threshold_forest()
    supplementary_qc(adatas)
    supplementary_sensitivity()
    print("STEP39_FIG1_S1_S2_COMPLETE")


if __name__ == "__main__":
    main()
