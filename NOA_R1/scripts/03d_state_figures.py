#!/usr/bin/env python3
"""Draw the three archived state plots with the existing Matplotlib layout.

Inputs in step3_state_analysis: expression_pair_coverage.csv (n_pairs), state_abundance_statistics.csv
(paired differences and intervals), expression_model_audit.csv (fitted units), target_gene_DE.csv (logFC and BH).
Input: step2_gene_sets.csv or step2_extended_specificity/locked_gene_sets_used.csv (ordered target_top20 genes).
Outputs in step3_state_analysis: state_pair_coverage, state_abundance_changes, state_target_gene_changes,
each as a PDF master and its PNG preview.
Run: python NOA_R1/scripts/03d_state_figures.py
"""

from __future__ import annotations

from importlib import import_module
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize


def main():
    """Render the state coverage, abundance, and differential expression figures. Returns nothing."""
    figure_style = import_module("00a_functions")
    figure_style.apply_publication_style()

    # Preserve the frozen state ordering on all three plots, with the two broad states first.
    state_order = [
        "MCG-like", "MAC_MDSC-like", "MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"
    ]
    state_labels = {state: f"{state}-like" for state in state_order}
    state_labels["MCG-like"] = "MCG-like"
    state_labels["MAC_MDSC-like"] = "MAC/MDSC-like"
    display_labels = [state_labels[state] for state in state_order]
    datasets = ("GSE174554", "GSE274546")

    # Plot one coverage count for each state, cohort, and prespecified per-endpoint minimum.
    revision = Path(__file__).resolve().parents[1]
    output_dir = revision / "step3_state_analysis"
    coverage_table = pd.read_csv(output_dir / "expression_pair_coverage.csv")
    coverage_table = coverage_table.set_index(["dataset", "state", "threshold"])
    cohort_thresholds = list(product(datasets, (20, 10)))
    coverage_matrix = np.empty((len(state_order), len(cohort_thresholds)), dtype=np.int64)
    for state_index, state in enumerate(state_order):
        for column_index, (dataset, threshold) in enumerate(cohort_thresholds):
            coverage_matrix[state_index, column_index] = coverage_table.loc[(dataset, state, threshold), "n_pairs"]

    # Use the original palette range, cell labels and six-pair fitting threshold.
    fig, axis = plt.subplots(figsize=(figure_style.mm(155), figure_style.mm(125)))
    image = axis.imshow(coverage_matrix, cmap="Blues", vmin=0, vmax=45, aspect="auto")
    for state_index, row in enumerate(coverage_matrix):
        for column_index, pair_count in enumerate(row):
            axis.text(
                column_index, state_index, str(pair_count),
                ha         = "center",
                va         = "center",
                color      = "white" if pair_count > 25 else "#222222",
                fontweight = "bold" if pair_count >= 6 else "normal"
            )

    axis.set_yticks(range(len(state_order)), display_labels)
    column_labels = [f"{dataset}\n≥{threshold} cells per endpoint" for dataset, threshold in cohort_thresholds]
    axis.set_xticks(range(len(cohort_thresholds)), column_labels)
    axis.axvline(1.5, color="white", linewidth=2)
    axis.axhline(1.5, color="white", linewidth=2)
    axis.set_title("Complete patient pairs available for expression analysis", pad=10)
    fig.colorbar(image, ax=axis, pad=.025, label="Paired patients")

    fig.text(
        .53, .018, "Bold: ≥6 pairs eligible for a count model. No missing-state expression is imputed.",
        ha       = "center",
        fontsize = 6.8
    )
    fig.tight_layout(rect=(0, .045, 1, 1))
    figure_style.save_figure(fig, output_dir, "state_pair_coverage")

    # Use pooled interval limits so viewers can compare the two cohort panels on one scale.
    abundance_table = pd.read_csv(output_dir / "state_abundance_statistics.csv")
    fig, axes = plt.subplots(1, 2, figsize=(figure_style.mm(200), figure_style.mm(135)), sharey=True, sharex=True)
    lower = float(abundance_table.ci_low_pp.min())
    upper = float(abundance_table.ci_high_pp.max())
    span = upper - lower
    colors = figure_style.colors
    for axis, dataset in zip(axes, datasets):
        cohort_table = abundance_table[abundance_table.dataset == dataset].set_index("state").loc[state_order]
        y_positions = np.arange(len(state_order))
        means = cohort_table.mean_delta_pp.to_numpy()
        errors = np.vstack([means - cohort_table.ci_low_pp, cohort_table.ci_high_pp - means])

        axis.axvline(0, color=colors["neutral"], linestyle="--", linewidth=.75)
        axis.axhline(1.5, color=colors["grid"], linewidth=.7)
        axis.errorbar(
            means, y_positions,
            xerr       = errors,
            fmt        = "o",
            color      = colors[dataset.lower()],
            markersize = 3.8,
            capsize    = 2,
            linewidth  = 1
        )
        axis.set_yticks(y_positions, display_labels)
        axis.set_ylim(len(state_order) - .35, -.65)
        axis.set_xlim(lower - .10 * span, upper + .10 * span)
        axis.set_xlabel("Recurrent − primary\n(relative percentage points)")
        axis.set_title(f"{dataset} | {int(cohort_table.n_pairs.iloc[0])} paired patients", loc="left", pad=10)
        figure_style.clean_axis(axis)

    fig.text(
        .52, .045, "Points and bars: mean paired change and 95% paired-t CI (unadjusted).",
        ha       = "center",
        fontsize = 7
    )
    fig.text(
        .52, .015,
        "Fine-state inference: sign-flip P with BY9. Broad groups are complementary, not independent tests.",
        ha       = "center",
        fontsize = 6.8
    )
    fig.tight_layout(rect=(0, .09, 1, 1), w_pad=2)
    figure_style.save_figure(fig, output_dir, "state_abundance_changes")

    # Keep all fitted 20-cell models and the frozen target order, excluding the absent raw gene.
    audit_table = pd.read_csv(output_dir / "expression_model_audit.csv")
    audit_table = audit_table[audit_table.threshold == 20].copy()
    audit_table["state_order"] = audit_table.state.map({state: index for index, state in enumerate(state_order)})
    audit_table = audit_table.sort_values(["dataset", "state_order"])
    definitions_path = revision / "step2_gene_sets.csv"
    if not definitions_path.is_file():
        definitions_path = revision / "step2_extended_specificity/locked_gene_sets_used.csv"
    definitions_table = pd.read_csv(definitions_path)
    plotted_genes = (definitions_table.set_id == "target_top20") & definitions_table.gene.ne("AC253572.2")
    genes = definitions_table.loc[plotted_genes, "gene"].tolist()

    # Missing or filtered genes remain grey. Mark BH <0.05 in the frozen cohort-wide fitted family.
    de_table = pd.read_csv(output_dir / "target_gene_DE.csv")
    de_table = de_table[de_table.threshold == 20].set_index(["dataset", "state", "gene"])
    fold_change_matrix = np.full((len(genes), len(audit_table)), np.nan)
    significant = np.zeros(fold_change_matrix.shape, dtype=bool)
    for column_index, model in enumerate(audit_table.itertuples()):
        for gene_index, gene in enumerate(genes):
            index = (model.dataset, model.state, gene)
            if index not in de_table.index:
                continue
            result = de_table.loc[index]
            fold_change_matrix[gene_index, column_index] = result.logFC
            significant[gene_index, column_index] = result.bh_all_fitted_state_gene_tests < .05

    # Keep the original symmetric log-fold-change limits and separate the two annotation lines.
    fig, axis = plt.subplots(figsize=(figure_style.mm(190), figure_style.mm(145)))
    palette = plt.get_cmap("RdBu_r").copy()
    palette.set_bad("#DDDDDD")
    image = axis.imshow(fold_change_matrix, cmap=palette, norm=Normalize(-2.5, 2.5), aspect="auto")
    for gene_index, column_index in zip(*np.where(significant)):
        axis.plot(column_index, gene_index, "o", markersize=2, color="black")

    axis.set_yticks(range(len(genes)), genes, fontstyle="italic")
    model_labels = [
        f"{model.dataset.removeprefix('GSE')}\n{state_labels[model.state]}\nn={model.n_pairs}"
        for model in audit_table.itertuples()
    ]
    axis.set_xticks(range(len(audit_table)), model_labels, rotation=35, ha="right", fontsize=6.6)
    axis.set_title("Within-group target-gene changes | ≥20 cells per endpoint", pad=10)
    fig.colorbar(image, ax=axis, pad=.025, extend="both", label="edgeR log2 fold change (recurrent − primary)")

    fig.text(
        .53, .043,
        "Grey: gene did not pass condition-based expression filtering. AC253572.2 absent in both raw objects.",
        ha       = "center",
        fontsize = 6.7
    )
    fig.text(
        .53, .016,
        "Dot: BH<0.05 across all fitted state × gene tests within a cohort (20-cell threshold). "
        "All eligible models shown.",
        ha       = "center",
        fontsize = 6.5
    )
    fig.tight_layout(rect=(0, .075, 1, 1))
    figure_style.save_figure(fig, output_dir, "state_target_gene_changes")


if __name__ == "__main__":
    main()
