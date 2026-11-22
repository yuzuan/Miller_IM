#!/usr/bin/env python3
"""Plot the historical five-set program comparison and empirical null distributions.

Inputs: step2_gene_sets.csv or step2_extended_specificity/locked_gene_sets_used.csv (set_id, display_name);
    step2_gene_set_specificity/camera_results.csv, random_set_summary.csv and <dataset>_<scheme>_random_sets.csv.
Outputs: step2_gene_set_specificity/program_comparison.pdf and .png, random_gene_set_nulls.pdf and .png.
Run: python NOA_R1/scripts/02g_program_specificity_figure.py after 02c, with Poppler on PATH.
"""

import csv
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# Load the fixed plot order from the five-set definitions, or the corresponding archived rows.
revision = Path(__file__).resolve().parents[1]
definition_path = revision / "step2_gene_sets.csv"
if not definition_path.is_file():
    definition_path = revision / "step2_extended_specificity" / "locked_gene_sets_used.csv"

with definition_path.open(encoding="utf-8", newline="") as handle:
    definition_table = list(csv.DictReader(handle))

display_names = {}
for row in definition_table:
    if row["set_id"] != "Miller_systemic_top20":
        display_names[row["set_id"]] = row["display_name"]
set_ids = list(display_names)

output_dir = revision / "step2_gene_set_specificity"
with (output_dir / "camera_results.csv").open(encoding="utf-8", newline="") as handle:
    camera_table = list(csv.DictReader(handle))

# Plot each set's estimated-correlation p value and its original five-set BH result.
figure_tools = importlib.import_module("00a_functions")
figure_tools.apply_publication_style()
datasets = ("GSE174554", "GSE274546")
figure, axes = plt.subplots(1, 2, figsize=(figure_tools.mm(205), figure_tools.mm(110)), sharey=True)
for axis, dataset in zip(axes, datasets, strict=True):
    rows_by_set = {}
    for row in camera_table:
        if row["dataset"] == dataset and row["correlation_mode"] == "estimated":
            rows_by_set[row["set_id"]] = row

    values = np.array([-np.log10(float(rows_by_set[set_id]["p_value"])) for set_id in set_ids])
    color = figure_tools.colors[dataset.lower()]
    axis.barh(np.arange(len(set_ids)), values, height=0.58, color=color, alpha=0.83)
    axis.axvline(-np.log10(.05), color=figure_tools.colors["neutral"], linestyle="--", linewidth=.8)
    for y, set_id in enumerate(set_ids):
        row = rows_by_set[set_id]
        label = (
            f"{row['n_genes']} genes; {row['direction']}; P={float(row['p_value']):.3g}; "
            f"q={float(row['bh_five_sets_within_cohort']):.3g}"
        )
        axis.text(0.02, y + .35, label, fontsize=6.8, va="center", color=figure_tools.colors["text"])

    # The shared axis names the sets. The annotation in each panel gives that cohort's coverage.
    axis.set_yticks(np.arange(len(set_ids)), [display_names[set_id] for set_id in set_ids])
    axis.set_ylim(len(set_ids) - .3, -.7)
    axis.set_xlim(0, max(3.0, max(values) + .35))
    axis.set_xlabel("−log10(P), estimated-correlation CAMERA")
    axis.set_title(dataset, loc="left", pad=9)
    figure_tools.clean_axis(axis)

figure.suptitle("Prespecified program-width and external-program comparisons", fontsize=9, y=.99)
figure.text(
    .5,
    .015,
    "Paired design; dashed line P=0.05. q: BH across five sets within each cohort. All sets shown.",
    ha = "center",
    fontsize = 7,
)
figure.tight_layout(rect=(0, .07, 1, .95), w_pad=2)
figure_tools.save_panel(figure, output_dir, "program_comparison", "NOA R1 step2 program comparison")

# Compare stored null ES draws against each cohort's fixed target under both matching schemes.
with (output_dir / "random_set_summary.csv").open(encoding="utf-8", newline="") as handle:
    summary_table = list(csv.DictReader(handle))
summary_by_scheme = {(row["dataset"], row["scheme"]): row for row in summary_table}
figure, axes = plt.subplots(2, 2, figsize=(figure_tools.mm(190), figure_tools.mm(135)), sharex=True)
for i, dataset in enumerate(datasets):
    for j, scheme in enumerate(("size_only", "size_and_expression")):
        axis = axes[i, j]
        draw_path = output_dir / f"{dataset}_{scheme}_random_sets.csv"
        with draw_path.open(encoding="utf-8", newline="") as handle:
            draw_table = list(csv.DictReader(handle))
        values = np.array([float(row["ES"]) for row in draw_table])
        summary = summary_by_scheme[dataset, scheme]
        n_draws = int(summary["n_random_sets"])
        if len(values) != n_draws:
            raise ValueError(f"Draw count differs from the summary for {dataset}, {scheme}")

        # The original 60 bins span the full ES range. Label the observed ES and empirical upper tail.
        axis.hist(
            values,
            bins      = np.linspace(-1, 1, 61),
            color     = figure_tools.colors[dataset.lower()],
            alpha     = .75,
            edgecolor = "white",
            linewidth = .15,
        )
        axis.axvline(float(summary["observed_ES"]), color=figure_tools.colors["recurrent"], linewidth=1.5)
        label = (
            f"{dataset}\nObserved ES={float(summary['observed_ES']):.3f}\n"
            f"{summary['positive_extreme_count']}/{n_draws:,} ≥ observed\n"
            f"Empirical P={float(summary['empirical_positive_p']):.3g}"
        )
        axis.text(.03, .95, label, transform=axis.transAxes, va="top", fontsize=7)

        # Add shared axes and the original scheme labels.
        if i == 0:
            axis.set_title("Size-matched" if j == 0 else "Size + expression matched", pad=8)
        if i == 1:
            axis.set_xlabel("Weighted GSEA enrichment score (ES)")
        if j == 0:
            axis.set_ylabel("Random gene sets")
        axis.set_xlim(-1, 1)
        axis.set_ylim(0, axis.get_ylim()[1] * 1.15)
        figure_tools.clean_axis(axis)

# Both cohorts use the same measured target size. Read it from the summary for the caption.
set_sizes = {int(row["set_size"]) for row in summary_table}
if len(set_sizes) != 1:
    raise ValueError("Target set size differs across random-set comparisons")
set_size = set_sizes.pop()
figure.text(
    .5,
    .012,
    f"Red line: fixed top20 ({set_size} measurable). "
    "Gene-membership resampling does not preserve intergene correlation.",
    ha = "center",
    fontsize = 7,
)
figure.tight_layout(rect=(0, .04, 1, 1), h_pad=1.8, w_pad=2)
figure_tools.save_panel(figure, output_dir, "random_gene_set_nulls", "NOA R1 step2 random gene-set nulls")
