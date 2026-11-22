#!/usr/bin/env python3
# Purpose: Draw cohort evidence, quality-control and sensitivity panels from two AnnData cohorts and Step 38.
# Inputs: The source tables named below, with published fields such as dataset, gene, logFC, F and NES.
# Outputs: Source tables and manifests under write, plus PDF and PNG panels under figures.
# Run: python scripts/05a_cohort_evidence.py from the project environment.
# Python preserves the existing AnnData, NumPy, pandas and Matplotlib reproduction workflow.
# Keep the original thresholds, gene sets, seeds and test order. Render the 600 dpi PNGs from PDF.
# Dataset names, third-party API parameters and output column names are public interfaces.

from __future__ import annotations
from pathlib import Path
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

if __name__ == "__main__":

    # Keep the seed at this point to reproduce the sampling or jitter sequence.
    np.random.seed(20260713)
    mpl.rcParams.update(
        {
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
        }
    )
    sns.set_style("white")
    root = Path(__file__).resolve().parents[1]
    write_root = root / "write" / "39_independent_miller_mg_inflammatory_figures"

    # Use the published output directory and basename.
    fig_root = root / "figures" / "39_independent_miller_mg_inflammatory_plottie"
    fig1_data = write_root / "Figure1" / "source_data"
    fig1_out = fig_root / "Figure1" / "panel_library"
    s1_data = write_root / "SupplementaryFigure1" / "source_data"
    s1_out = fig_root / "SupplementaryFigure1" / "panel_library"

    # Use the published output directory and basename.
    s2_data = write_root / "SupplementaryFigure2" / "source_data"
    s2_out = fig_root / "SupplementaryFigure2" / "panel_library"
    for path in [fig1_data, fig1_out, s1_data, s1_out, s2_data, s2_out]:
        path.mkdir(parents=True, exist_ok=True)
    h5ads = {
        "GSE174554": root
        / "write/34_gse174554_raw_independent_discovery/03_myeloid/GSE174554_raw_independent_clean_myeloid.h5ad",
        "GSE274546": root
        / "write/36_gse274546_raw_independent_reannotation/03_myeloid/GSE274546_raw_independent_clean_myeloid.h5ad",
    }
    adatas = {dataset: ad.read_h5ad(path) for dataset, path in h5ads.items()}

    # Use the published output directory and basename.
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"
    audit_path = step38 / "independent_input_pair_audit.csv"
    audit_table = pd.read_csv(audit_path)
    required = {"dataset", "threshold", "n_pairs", "n_cells_primary", "n_cells_recurrent"}
    if not required.issubset(audit_table.columns):
        raise ValueError(f"Missing Step38 pair-audit columns: {sorted(required - set(audit_table.columns))}")

    datasets = ["GSE174554", "GSE274546"]
    formal_table = audit_table.loc[
        audit_table["dataset"].isin(datasets) & audit_table["threshold"].eq(20)
    ].copy()
    if "formal_testable" in formal_table.columns:
        is_formal = formal_table["formal_testable"].astype(str).str.lower().isin({"true", "t", "1"})
        formal_table = formal_table.loc[is_formal].copy()
    if formal_table.groupby("dataset").size().to_dict() != {dataset: 1 for dataset in datasets}:
        raise ValueError(
            (
                "Step38 pair audit must contain exactly one formal threshold-20 ro"
                "w per cohort: "
                f"{formal_table.groupby('dataset').size().to_dict()}"
            )
        )
    formal_table["clean_myeloid"] = (
        pd.to_numeric(formal_table["n_cells_primary"], errors="raise")
        + pd.to_numeric(formal_table["n_cells_recurrent"], errors="raise")
    ).astype(int)

    # Keep the formal analysis rows in the existing audit and check the unit of analysis for each cohort.
    formal_table["formal_pairs_threshold20"] = pd.to_numeric(formal_table["n_pairs"], errors="raise").astype(
        int
    )
    formal_table["input_libraries_or_matrices"] = formal_table["dataset"].map(
        {"GSE174554": 91, "GSE274546": 111}
    )
    formal_table["cell_count_scope"] = "IDH-wildtype cells included in the formal threshold-20 analysis"
    source_table = (
        formal_table[
            [
                "dataset",
                "input_libraries_or_matrices",
                "clean_myeloid",
                "formal_pairs_threshold20",
                "n_cells_primary",
                "n_cells_recurrent",
                "cell_count_scope",
            ]
        ]
        .sort_values("dataset")
        .reset_index(drop=True)
    )

    # Build the design panel from the audit rows, using the existing cell and patient-pair counts.
    workflow_table = source_table.set_index("dataset")
    fig, ax = plt.subplots(figsize=(8.1, 3.1))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    # The helper below uses the existing arguments and return fields.
    dataset_colors = {"GSE174554": "#3C5488", "GSE274546": "#E64B35"}

    def box(x, y, w, h, text, edge, fill="white", lw=1.2, weight="normal"):
        """Draw a labeled box at the requested position and size on the workflow axes. No return value."""

        rect = mpl.patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.04",
            facecolor=fill,
            edgecolor=edge,
            linewidth=lw,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8, weight=weight)

    box(
        0.3,
        4.8,
        5.1,
        1.0,
        "GSE174554  |  91 matrices  |  snRNA-seq\nPrimary -> 1st recurrent",
        dataset_colors["GSE174554"],
        fill="#F3F6FA",
        weight="bold",
    )
    box(
        6.6,
        4.8,
        5.1,
        1.0,
        "GSE274546  |  111 libraries  |  scRNA-seq\nPrimary -> 1st recurrent",
        dataset_colors["GSE274546"],
        fill="#FDF2EF",
        weight="bold",
    )
    for x in (2.85, 9.15):
        ax.annotate(
            "", xy=(x, 4.15), xytext=(x, 4.75), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0)
        )

    box(0.65, 3.15, 4.4, 0.82, "Independent QC, doublet removal,\nHVG/PCA/Harmony/Leiden", "#666666")
    box(6.95, 3.15, 4.4, 0.82, "Independent QC, doublet removal,\nHVG/PCA/Harmony/Leiden", "#666666")
    for x in (2.85, 9.15):
        ax.annotate(
            "", xy=(x, 2.48), xytext=(x, 3.12), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0)
        )
    gse174554 = workflow_table.loc["GSE174554"]
    gse274546 = workflow_table.loc["GSE274546"]

    box(
        0.65,
        1.45,
        4.4,
        0.85,
        (
            f"{int(gse174554.clean_myeloid):,}"
            " analyzed myeloid cells\n"
            f"{int(gse174554.formal_pairs_threshold20)}"
            " paired patients (>=20 cells/side)"
        ),
        dataset_colors["GSE174554"],
        fill="#F3F6FA",
    )
    box(
        6.95,
        1.45,
        4.4,
        0.85,
        (
            f"{int(gse274546.clean_myeloid):,}"
            " analyzed myeloid cells\n"
            f"{int(gse274546.formal_pairs_threshold20)}"
            " paired patients (>=20 cells/side)"
        ),
        dataset_colors["GSE274546"],
        fill="#FDF2EF",
    )
    for x in (2.85, 9.15):
        ax.annotate(
            "", xy=(6.0, 0.75), xytext=(x, 1.4), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0)
        )
    box(
        3.75,
        0.05,
        4.5,
        0.68,
        "Paired raw-count pseudobulk\nedgeR + fixed-program GSEA",
        "#00A087",
        fill="#F0FAF7",
        weight="bold",
    )
    ax.text(
        6, -0.18, "Patient is the statistical replicate", ha="center", va="top", color="#555555", fontsize=7
    )

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = fig1_data / "cohort_workflow.csv"
    source_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    def save_tight_panel(fig: plt.Figure, out_dir: Path, name: str) -> None:
        """Save a tightly cropped PDF under the given stem, then render its PNG. No return value."""

        # Use the published output directory and basename.
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.canvas.draw()
        fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
        render_pdf_png(out_dir / f"{name}.pdf", out_dir / f"{name}.png")
        plt.close(fig)

    save_tight_panel(fig, fig1_out, "cohort_independent_processing")

    identity_colors = {
        "Macrophage": "#3C5488",
        "Microglia": "#00A087",
        "Monocyte": "#F39B7F",
        "CDC": "#8491B4",
        "Neutrophil": "#7E6148",
    }
    for dataset, adata in adatas.items():
        coords = np.asarray(adata.obsm["X_umap"])
        labels = adata.obs["myeloid_identity"].astype(str).to_numpy()
        order = np.argsort(labels)
        coords = coords[order]
        labels = labels[order]

        # Draw the panel with the published canvas dimensions, placements and styling.
        fig, ax = plt.subplots(figsize=(4.0, 3.8))
        for identity in sorted(np.unique(labels)):
            mask = labels == identity
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=1.0,
                alpha=0.7,
                linewidths=0,
                color=identity_colors.get(identity, "#BDBDBD"),
                rasterized=True,
                label=identity,
            )
            center = np.median(coords[mask], axis=0)
            ax.text(
                center[0],
                center[1],
                identity,
                ha="center",
                va="center",
                fontsize=7,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.8),
            )
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        ax.legend(
            frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), markerscale=4, handletextpad=0.2
        )
        output_table = pd.DataFrame(
            {
                "cell_id": adata.obs_names.to_numpy(),
                "UMAP1": np.asarray(adata.obsm["X_umap"])[:, 0],
                "UMAP2": np.asarray(adata.obsm["X_umap"])[:, 1],
                "myeloid_identity": adata.obs["myeloid_identity"].astype(str).to_numpy(),
            }
        )
        csv_path = fig1_data / f"umap_identity_{dataset}.csv.gz"

        # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
        output_table.to_csv(csv_path, index=False, compression={"method": "gzip", "mtime": 0})
        if not csv_path.is_file() or csv_path.stat().st_size == 0:
            raise ValueError(f"Missing or empty output: {csv_path.name}")
        save_tight_panel(fig, fig1_out, f"UMAP_identity_{dataset}")
    markers = {
        "Microglia": ["P2RY12", "TMEM119", "SALL1", "CX3CR1"],
        "Macrophage": ["GPNMB", "SPP1", "MRC1", "MSR1"],
        "Monocyte": ["VCAN", "FCN1", "S100A8", "S100A9"],
        "CDC": ["FCER1A", "CD1C", "CLEC10A", "HLA-DPB1"],
    }
    marker_genes = [g for group in markers.values() for g in group]
    rows = []

    # The helper below uses the existing arguments and return fields.
    def dense_mean_and_pct(x) -> tuple[np.ndarray, np.ndarray]:
        """Return gene means and percentages of positive cells from a sparse or dense cell-by-gene matrix."""

        # Use the existing input checks and branch conditions.
        if sparse.issparse(x):
            mean = np.asarray(x.mean(axis=0)).ravel()
            pct = np.asarray((x > 0).mean(axis=0)).ravel() * 100
        else:
            mean = np.asarray(x).mean(axis=0)
            pct = (np.asarray(x) > 0).mean(axis=0) * 100
        result = (mean, pct)
        return result

    for dataset, adata in adatas.items():
        for identity in sorted(adata.obs["myeloid_identity"].astype(str).unique()):
            mask = adata.obs["myeloid_identity"].astype(str).eq(identity).to_numpy()
            present = [g for g in marker_genes if g in adata.var_names]
            mean, pct = dense_mean_and_pct(adata[mask, present].X)
            for gene, value, fraction in zip(present, mean, pct):
                rows.append(
                    {
                        "dataset": dataset,
                        "identity": identity,
                        "gene": gene,
                        "mean_expression": value,
                        "pct_positive": fraction,
                    }
                )
    marker_table = pd.DataFrame(rows)

    # Store group results in a Series of the same dtype to preserve the floating-point values written to CSV.
    scaled_expression = marker_table["mean_expression"].copy()

    # Standardize each cohort and gene across cell identities using ddof=0 and a denominator of 1 for zero variance.
    for _, expression in marker_table.groupby(["dataset", "gene"])["mean_expression"]:
        expression_sd = expression.std(ddof=0)
        if expression_sd > 0:
            denominator = expression_sd
        else:
            denominator = 1
        scaled_expression.loc[expression.index] = (
            expression - expression.mean()
        ) / denominator

    # Assign the completed column by the original row index, preserving its dtype for export and plotting.
    marker_table["scaled_expression"] = scaled_expression

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = fig1_data / "identity_marker_dotplot.csv"
    marker_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    dot_table = marker_table
    genes = [g for group in markers.values() for g in group]

    identities = [
        x for x in ["Microglia", "Macrophage", "Monocyte", "CDC"] if x in set(dot_table["identity"])
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), sharey=True, gridspec_kw={"wspace": 0.06})
    for ax, dataset in zip(axes, ["GSE174554", "GSE274546"]):
        dataset_table = dot_table[dot_table["dataset"].eq(dataset)].copy()
        dataset_table["x"] = dataset_table["gene"].map({g: i for i, g in enumerate(genes)})
        dataset_table["y"] = dataset_table["identity"].map({g: i for i, g in enumerate(identities)})
        ax.scatter(
            dataset_table["x"],
            dataset_table["y"],
            s    = 8 + dataset_table["pct_positive"] * 0.8,
            c    = dataset_table["scaled_expression"],
            cmap = "coolwarm",
            vmin = -2,
            vmax = 2,
            edgecolors = "none",
        )
        ax.set_xlim(-0.6, len(genes) - 0.4)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.set_ylim(len(identities) - 0.4, -0.6)
        ax.set_xticks(range(len(genes)), genes, rotation=55, ha="right")
        ax.set_yticks(range(len(identities)), identities)
        ax.grid(color="#ECECEC", lw=0.5)
        ax.set_axisbelow(True)

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        ax.text(0.5, 1.03, dataset, transform=ax.transAxes, ha="center", va="bottom", weight="bold")
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-2, 2), cmap="coolwarm")
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)

    cbar.set_label("Scaled expression")
    save_tight_panel(fig, fig1_out, "identity_marker_dotplot")
    deg_table = pd.read_csv(step38 / "independent_paired_edger_all_genes.csv")
    gsea_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 5.6), sharex=False)

    sources = []

    raw20 = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3",
        "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]

    def enrichment_curve_legacy(deg_table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
        """Return the running enrichment table for a differential expression table and gene set using Step 39 rules."""

        # Prepare the cohort gene ranks using the original order and merge keys.
        ranked_table = deg_table.copy()
        ranked_table["rank_stat"] = np.sign(ranked_table["logFC"]) * np.sqrt(np.maximum(ranked_table["F"], 0))
        ranked_table = (
            ranked_table[np.isfinite(ranked_table["rank_stat"])]
            .sort_values("rank_stat", ascending=False)
            .reset_index(drop=True)
        )
        hits = ranked_table["gene"].isin(genes).to_numpy()
        weights = np.abs(ranked_table["rank_stat"].to_numpy())

        hit_weights = np.where(hits, weights, 0)
        hit_norm = hit_weights.sum()
        miss_norm = (~hits).sum()
        increments = np.where(hits, hit_weights / hit_norm, -1 / miss_norm)
        running = np.cumsum(increments)

        # Build the table from the named records in their original order and with the same unit per row.
        result_table = pd.DataFrame(
            {
                "rank": np.arange(1, len(ranked_table) + 1),
                "gene": ranked_table["gene"],
                "rank_stat": ranked_table["rank_stat"],
                "hit": hits,
                "running_enrichment": running,
            }
        )
        return result_table

    for ax, dataset in zip(axes, ["GSE174554", "GSE274546"]):
        current_table = deg_table[deg_table["dataset"].eq(dataset) & deg_table["threshold"].eq(20)]
        curve_table = enrichment_curve_legacy(current_table, raw20)
        curve_table["dataset"] = dataset
        sources.append(curve_table)
        color = dataset_colors[dataset]

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.plot(curve_table["rank"], curve_table["running_enrichment"], color=color, lw=1.6)
        ax.fill_between(curve_table["rank"], 0, curve_table["running_enrichment"], color=color, alpha=0.1)
        hit_ranks = curve_table.loc[curve_table["hit"], "rank"]
        ax.vlines(hit_ranks, ymin=-0.025, ymax=0.025, color="#555555", lw=0.5)
        row = gsea_table[
            gsea_table["dataset"].eq(dataset)
            & gsea_table["threshold"].eq(20)
            & gsea_table["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")
        ].iloc[0]

        # Set the published coordinates, marks, labels and axis ranges for this panel.
        ax.text(
            0.02,
            0.9,
            f"{dataset}   NES {row.NES:.2f}   nominal P {row.pval:.3g}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            weight="bold",
        )
        ax.axhline(0, color="#999999", lw=0.6)
        ax.set_ylabel("Running enrichment")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Ranked genes")

    # Build the table from the named records in their original order and with the same unit per row.
    output_table_total = pd.concat(sources, ignore_index=True)
    csv_path = fig1_data / "raw20_gsea_curves.csv.gz"
    output_table_total.to_csv(csv_path, index=False, compression={"method": "gzip", "mtime": 0})
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    save_tight_panel(fig, fig1_out, "Miller_raw20_GSEA_curves")

    # Read the frozen inputs with their original row order, dataset names and column names.
    gsea_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
    threshold_table = gsea_table[gsea_table["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()
    csv_path = fig1_data / "raw20_threshold_forest.csv"
    threshold_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    ymap = {"GSE174554": 1, "GSE274546": 0}
    offsets = {20: 0.1, 50: -0.1}
    for _, row in threshold_table.iterrows():
        y = ymap[row.dataset] + offsets[int(row.threshold)]
        is_formal = bool(row.formal_testable)
        ax.scatter(
            row.NES,
            y,
            s=40,
            color=dataset_colors[row.dataset] if is_formal else "white",
            edgecolor=dataset_colors[row.dataset],
            linewidth=1.1,
            marker="o" if row.threshold == 20 else "s",
            zorder=3,
        )
        ax.text(row.NES + 0.025, y, f"{row.NES:.2f}", va="center", fontsize=7)
    ax.axvline(0, color="#999999", lw=0.7)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_yticks([0, 1], ["GSE274546", "GSE174554"])
    ax.set_xlabel("Normalized enrichment score")
    ax.set_ylim(-0.45, 1.45)
    ax.spines[["top", "right", "left"]].set_visible(False)
    handles = [
        mpl.lines.Line2D([], [], marker="o", color="#555555", linestyle="", label=">=20 cells/side"),
        mpl.lines.Line2D([], [], marker="s", color="#555555", linestyle="", label=">=50 cells/side"),
        mpl.lines.Line2D(
            [],
            [],
            marker="o",
            markerfacecolor="white",
            markeredgecolor="#555555",
            linestyle="",
            label="<10 pairs: descriptive",
        ),
    ]

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.legend(handles=handles, frameon=False, loc="lower right")
    save_tight_panel(fig, fig1_out, "raw20_threshold_stability")
    summary_rows, capture_rows = ([], [])
    step_dirs = {
        "GSE174554": root / "write/34_gse174554_raw_independent_discovery",
        "GSE274546": root / "write/36_gse274546_raw_independent_reannotation",
    }
    for dataset, base in step_dirs.items():

        # Read the frozen inputs with their original row order, dataset names and column names.
        capture_table = pd.read_csv(base / "01_qc_doublets/capture_qc_doublet_summary.csv")
        capture_table["dataset"] = dataset
        capture_rows.append(capture_table)
        summary_rows.append(
            {
                "dataset": dataset,
                "raw": int(capture_table["n_raw"].sum()),
                "qc": int(capture_table["n_qc"].sum()),
                "singlet": int(capture_table.loc[capture_table["status"].eq("ok"), "n_singlet"].sum()),
                "clean_myeloid": int(adatas[dataset].n_obs),
            }
        )

    # Build the table from the named records in their original order and with the same unit per row.
    summary_table = pd.DataFrame(summary_rows)
    capture_table_total = pd.concat(capture_rows, ignore_index=True)
    csv_path = s1_data / "qc_stage_counts.csv"
    summary_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = s1_data / "capture_qc_doublets.csv"
    capture_table_total.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    stage_table = summary_table.melt(id_vars="dataset", var_name="stage", value_name="cells")
    fig, ax = plt.subplots(figsize=(4.4, 3.1))

    sns.barplot(data=stage_table, x="stage", y="cells", hue="dataset", palette=dataset_colors, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Cells")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.legend(frameon=False)
    save_tight_panel(fig, s1_out, "qc_stage_counts")
    fig, ax = plt.subplots(figsize=(4.3, 3.0))
    sns.boxplot(
        data=capture_table_total[capture_table_total["status"].eq("ok")],
        x="dataset",
        y="doublet_rate",
        hue="dataset",
        palette=dataset_colors,
        legend=False,
        width=0.5,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=capture_table_total[capture_table_total["status"].eq("ok")],
        x="dataset",
        y="doublet_rate",
        color="#333333",
        size=1.8,
        alpha=0.45,
        jitter=0.18,
        ax=ax,
    )

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlabel("")
    ax.set_ylabel("Doublet rate")
    ax.spines[["top", "right"]].set_visible(False)
    save_tight_panel(fig, s1_out, "capture_doublet_rates")
    meta_rows = []

    threshold_rows = []
    identity_rows = []
    review_rows = []
    for dataset, base in step_dirs.items():

        # Read the frozen inputs with their original row order, dataset names and column names.
        metadata_table = pd.read_csv(base / "03_myeloid/patient_condition_pseudobulk_metadata.csv")
        metadata_table["dataset"] = dataset
        meta_rows.append(metadata_table)
        for threshold in [20, 50, 100, 200]:
            eligible_table = metadata_table[metadata_table["n_cells"].ge(threshold)]
            counts = eligible_table.groupby("patient_id")["condition"].nunique()
            threshold_rows.append(
                {"dataset": dataset, "threshold": threshold, "n_pairs": int(counts.ge(2).sum())}
            )
        identities = adatas[dataset].obs["myeloid_identity"].astype(str).value_counts()

        # Process cohorts, genes or panels in the published order.
        for identity, n in identities.items():
            identity_rows.append({"dataset": dataset, "identity": identity, "n_cells": int(n)})
        review_table = pd.read_csv(base / "03_myeloid/myeloid_cluster_review.csv")
        review_table["dataset"] = dataset
        review_rows.append(review_table)
    metadata_table_total = pd.concat(meta_rows, ignore_index=True)

    # Build the table from the named records in their original order and with the same unit per row.
    threshold_table = pd.DataFrame(threshold_rows)
    identity_table = pd.DataFrame(identity_rows)
    review_table_total = pd.concat(review_rows, ignore_index=True)
    csv_path = s1_data / "patient_condition_cell_counts.csv"
    metadata_table_total.to_csv(csv_path, index=False)

    # Use the existing input checks and branch conditions.
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s1_data / "pair_count_thresholds.csv"
    threshold_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s1_data / "identity_cell_counts.csv"

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    identity_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s1_data / "cluster_contamination_review.csv"
    review_table_total.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig, ax = plt.subplots(figsize=(4.5, 3.1))
    sns.boxplot(
        data=metadata_table_total,
        x="condition",
        y="n_cells",
        hue="dataset",
        palette=dataset_colors,
        fliersize=0,
        width=0.65,
        ax=ax,
    )
    sns.stripplot(
        data=metadata_table_total,
        x="condition",
        y="n_cells",
        hue="dataset",
        dodge=True,
        palette=dataset_colors,
        size=1.5,
        alpha=0.35,
        ax=ax,
        legend=False,
    )
    ax.set_yscale("log")
    ax.set_xlabel("")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_ylabel("Clean myeloid cells per patient endpoint")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_tight_panel(fig, s1_out, "patient_endpoint_cell_counts")
    fig, ax = plt.subplots(figsize=(4.4, 3.0))

    # Process cohorts, genes or panels in the published order.
    for dataset, dataset_table in threshold_table.groupby("dataset"):
        ax.plot(
            dataset_table["threshold"],
            dataset_table["n_pairs"],
            marker="o",
            lw=1.4,
            color=dataset_colors[dataset],
            label=dataset,
        )
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200], [20, 50, 100, 200])
    ax.set_xlabel("Minimum cells per endpoint")
    ax.set_ylabel("Eligible patient pairs")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_tight_panel(fig, s1_out, "pair_count_thresholds")
    fig, ax = plt.subplots(figsize=(4.5, 3.1))
    identity_pivot_table = identity_table.pivot(index="dataset", columns="identity", values="n_cells").fillna(
        0
    )

    identity_pivot_table = identity_pivot_table[
        [c for c in identity_colors if c in identity_pivot_table.columns]
    ]
    bottom = np.zeros(len(identity_pivot_table))
    for identity in identity_pivot_table.columns:
        values = identity_pivot_table[identity].to_numpy()
        ax.bar(
            identity_pivot_table.index, values, bottom=bottom, color=identity_colors[identity], label=identity
        )
        bottom += values
    ax.set_ylabel("Clean myeloid cells")
    ax.spines[["top", "right"]].set_visible(False)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_tight_panel(fig, s1_out, "identity_cell_counts")
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    review_plot_table = review_table_total.copy()
    review_plot_table["cluster_label"] = (
        review_plot_table["dataset"].str.replace("GSE", "") + ":" + review_plot_table["cluster"].astype(str)
    )

    colors = review_plot_table["is_clean_myeloid"].map({True: "#00A087", False: "#BDBDBD"})
    ax.scatter(
        review_plot_table["top20_contamination_hits"],
        review_plot_table["ribosomal_top20_hits"],
        s=10 + np.sqrt(review_plot_table["n_cells"]) * 2,
        c=colors,
        alpha=0.75,
        edgecolors="white",
        linewidth=0.3,
    )
    ax.axvline(3, color="#888888", ls="--", lw=0.7)
    ax.axhline(6, color="#888888", ls="--", lw=0.7)
    ax.set_xlabel("Contamination markers among top 20")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_ylabel("Ribosomal markers among top 20")
    ax.spines[["top", "right"]].set_visible(False)
    save_tight_panel(fig, s1_out, "cluster_contamination_audit")
    targeted_table = pd.read_csv(step38 / "independent_fixed_program_targeted_gsea.csv")
    score_table = pd.read_csv(step38 / "independent_fixed_program_paired_scores.csv")

    # Read the frozen inputs with their original row order, dataset names and column names.
    gene_table = pd.read_csv(step38 / "independent_fixed_program_gene_direction.csv")
    panel_table = pd.read_csv(step38 / "independent_old_panel_gsea_recomputed.csv")
    raw_table = targeted_table[
        targeted_table["pathway"].eq("Miller_Microglial_Inflammatory_raw_top20")
    ].copy()
    raw_score_table = score_table[
        score_table["signature"].eq("Miller_Microglial_Inflammatory_raw_top20")
    ].copy()
    raw_gene_table = gene_table[gene_table["signature"].eq("Miller_Microglial_Inflammatory_raw_top20")].copy()

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    csv_path = s2_data / "raw20_gsea.csv"
    raw_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s2_data / "raw20_paired_scores.csv"
    raw_score_table.to_csv(csv_path, index=False)

    # Use the existing input checks and branch conditions.
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s2_data / "raw20_gene_directions.csv"
    raw_gene_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    csv_path = s2_data / "full_panel_gsea.csv"

    # Write the source table with its original columns, ordering and compression, then check that it is nonempty.
    panel_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    fig, ax = plt.subplots(figsize=(4.8, 3.1))
    sns.pointplot(
        data=raw_table,
        x="threshold",
        y="NES",
        hue="dataset",
        palette=dataset_colors,
        dodge=0.12,
        markers="o",
        linestyles="-",
        errorbar=None,
        ax=ax,
    )
    ax.set_xlabel("Minimum cells per endpoint")

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_ylabel("NES")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    save_tight_panel(fig, s2_out, "raw20_threshold_sensitivity")
    score_gsea_table_total = raw_score_table.merge(
        raw_table,
        left_on=["dataset", "threshold", "signature"],
        right_on=["dataset", "threshold", "pathway"],
        how="inner",
    )

    score_gsea_table_total = score_gsea_table_total[score_gsea_table_total["threshold"].eq(20)]
    fig, ax = plt.subplots(figsize=(4.3, 3.2))
    for _, row in score_gsea_table_total.iterrows():
        ax.scatter(
            row.mean_delta,
            row.NES,
            color=dataset_colors[row.dataset],
            s=50,
            marker="o" if "top20" in row.signature else "s",
        )
        ax.text(row.mean_delta, row.NES + 0.025, row.dataset.replace("GSE", ""), ha="center", fontsize=6.5)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.axhline(0, color="#999999", lw=0.6)

    # Set the published coordinates, marks, labels and axis ranges for this panel.
    ax.set_xlabel("Mean patient score change")
    ax.set_ylabel("GSEA NES")
    ax.spines[["top", "right"]].set_visible(False)
    save_tight_panel(fig, s2_out, "average_score_vs_gsea")
    gene_direction_table = raw_gene_table[raw_gene_table["threshold"].eq(20)]

    # Prepare the published gene effects at the original row granularity and with the original statistical fields.
    effect_table = gene_direction_table.pivot(index="gene", columns="dataset", values="logFC").reindex(raw20)
    fig, ax = plt.subplots(figsize=(3.8, 5.4))
    sns.heatmap(
        effect_table,
        cmap="coolwarm",
        center=0,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Paired logFC", "shrink": 0.55},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Save the PDF master at the published size and render its PNG from that PDF.
    save_tight_panel(fig, s2_out, "raw20_gene_logfc")
    selected_pathway_table = panel_table[panel_table["threshold"].eq(20)].copy()
    keep = selected_pathway_table.groupby("pathway")["pval"].min().sort_values().head(14).index
    selected_pathway_table = selected_pathway_table[selected_pathway_table["pathway"].isin(keep)]
    effect_table = selected_pathway_table.pivot(index="pathway", columns="dataset", values="NES").reindex(
        keep
    )

    # Draw the panel with the published canvas dimensions, placements and styling.
    fig, ax = plt.subplots(figsize=(4.2, 5.0))
    sns.heatmap(
        effect_table,
        cmap="coolwarm",
        center=0,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "NES", "shrink": 0.55},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    save_tight_panel(fig, s2_out, "independent_program_gsea_landscape")
