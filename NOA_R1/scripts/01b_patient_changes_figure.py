#!/usr/bin/env python3
"""Plot the cohort-specific paired recurrence score changes.

Inputs: step1_primary_recurrence/patient_paired_changes.csv (dataset, delta) and
    patient_score_effects.csv (dataset, n_pairs, mean_delta, ci_low, ci_high, t_p_value).
Outputs: step1_primary_recurrence/patient_paired_changes.pdf and .png.
Run: python NOA_R1/scripts/01b_patient_changes_figure.py after 01a, with Poppler on PATH.
"""

import csv
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# One input table has one row per patient. The other has one paired effect estimate per cohort.
revision = Path(__file__).resolve().parents[1]
output_dir = revision / "step1_primary_recurrence"
with (output_dir / "patient_paired_changes.csv").open(encoding="utf-8", newline="") as handle:
    pair_table = list(csv.DictReader(handle))
with (output_dir / "patient_score_effects.csv").open(encoding="utf-8", newline="") as handle:
    effect_table = list(csv.DictReader(handle))
effect_by_dataset = {row["dataset"]: row for row in effect_table}

# Preserve the original single-panel coordinates, cohort colors, and 200 x 100 mm layout.
figure_tools = importlib.import_module("00a_functions")
figure_tools.apply_publication_style()
figure, axes = plt.subplots(1, 2, figsize=(figure_tools.mm(200), figure_tools.mm(100)))
for axis, dataset in zip(axes, ("GSE174554", "GSE274546"), strict=True):
    deltas = np.array(sorted(float(row["delta"]) for row in pair_table if row["dataset"] == dataset))
    effect = effect_by_dataset[dataset]
    mean = float(effect["mean_delta"])
    low = float(effect["ci_low"])
    high = float(effect["ci_high"])
    n_pairs = int(effect["n_pairs"])
    if len(deltas) != n_pairs or not np.isfinite(deltas).all():
        raise ValueError(f"Invalid patient changes for {dataset}")

    # Draw the ordered patient changes and the paired t interval at the original right-hand summary position.
    color = figure_tools.colors[dataset.lower()]
    axis.axhline(0, color=figure_tools.colors["neutral"], linestyle="--", linewidth=0.7, zorder=1)
    axis.scatter(np.arange(1, n_pairs + 1), deltas, color=color, s=15, edgecolors="white", linewidth=0.3, zorder=3)
    summary_x = n_pairs + max(3, n_pairs * 0.12)
    axis.axvline(n_pairs + max(1.5, n_pairs * 0.06), color=figure_tools.colors["light_neutral"], linewidth=0.6)
    axis.errorbar(
        summary_x,
        mean,
        yerr       = [[mean - low], [high - mean]],
        fmt        = "D",
        color      = "black",
        capsize    = 3,
        markersize = 4,
        linewidth  = 1.1,
        zorder     = 4,
    )

    # Use the saved effect and p value for the panel annotation.
    axis.set_xticks([1, n_pairs, summary_x], ["1", str(n_pairs), "Mean\n95% CI"])
    axis.set_xlim(0, summary_x + max(2, n_pairs * 0.05))
    axis.set_xlabel("Patients ordered by score change")
    axis.set_ylabel("Recurrent − primary\n(cohort-standardized program score)")
    axis.set_title(f"{dataset}  |  {n_pairs} paired patients", loc="left", pad=9)
    effect_label = f"Mean {mean:.3f} (95% CI {low:.3f} to {high:.3f})\nPaired t P = {float(effect['t_p_value']):.3g}"
    axis.text(0.02, 0.97, effect_label, transform=axis.transAxes, va="top", fontsize=7)
    axis.margins(y=0.22)
    figure_tools.clean_axis(axis)

# Export the PDF master and derive the original 600 dpi PNG with the shared renderer.
figure.tight_layout(pad=1.4, w_pad=2.6)
figure_tools.save_panel(figure, output_dir, "patient_paired_changes", "NOA R1 step1 paired score changes")
