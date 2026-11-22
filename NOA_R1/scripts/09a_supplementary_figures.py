#!/usr/bin/env python3
"""Rebuild Supplementary Figures S1–S3 using the frozen October NOA R1 layout.

Inputs: Paired changes and effects from step1_primary_recurrence, random-set draws
    and random_summary.csv from step2_extended_specificity, camera_seven_set_family.csv
    from step2_hallmark_completion, and state abundance and camera results from step3_state_analysis.
Outputs: supplementary_R1/Supplementary_Figure_S{1,2,3}.{png,tiff}, matching _text.txt
    files and Supplementary_Figure_S2.pdf. S1 and S3 use temporary PDF masters.
Run: python NOA_R1/scripts/09a_supplementary_figures.py
"""

import importlib
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib

# Use the headless backend and let Matplotlib select its cache directory.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text
from PIL import Image


figure_tools = importlib.import_module("00a_functions")


def write_supplementary_figure(figure, output_dir, stem, *, fixed_canvas=False) -> None:
    """Render the figure's PDF as its original PNG and TIFF exports. No return value."""
    text_objects = [item for item in figure.findobj(match=Text) if item.get_text().strip()]
    small_labels = [item.get_text() for item in text_objects if item.get_fontsize() < 8.2]
    if small_labels:
        raise ValueError(f"{stem} contains text below 8.2 points: {small_labels[:10]}")

    # Capture the original label list before drawing again to measure raster bounds.
    labels = [item.get_text() for item in text_objects]
    label_path = output_dir / f"{stem}_text.txt"
    label_path.write_text("\n".join(labels) + "\n", encoding="utf-8")

    # S1 and S3 retain their original tight crop and padding. S2 uses its full print canvas.
    raster_sizes = {}
    crop_bounds = None
    for dpi in (220, 600):
        figure.set_dpi(dpi)
        figure.canvas.draw()
        width, height = figure.get_size_inches()
        if not fixed_canvas:
            renderer = figure.canvas.get_renderer()
            crop_bounds = figure.get_tightbbox(renderer).padded(0.12)
            width, height = crop_bounds.size
        raster_sizes[dpi] = (int(width * dpi), int(height * dpi))

    # Use one PDF master for both raster sizes. Keep a permanent PDF only for S2.
    with TemporaryDirectory() as temporary_dir:
        scratch = Path(temporary_dir)
        pdf_dir = output_dir if fixed_canvas else scratch
        pdf_path = pdf_dir / f"{stem}.pdf"
        figure.savefig(pdf_path, bbox_inches=crop_bounds)

        for dpi, suffix in ((220, ".png"), (600, ".tiff")):
            rendered_path = scratch / f"{dpi}.png"
            figure_tools.render_pdf_png(pdf_path, rendered_path, dpi=dpi)

            # Preserve the published pixel rounding and the white background of both exports.
            mode = "RGBA" if suffix == ".png" else "RGB"
            with Image.open(rendered_path) as rendered:
                raster = Image.new(mode, raster_sizes[dpi], "white")
                raster.paste(rendered.convert(mode), (0, 0))

            output_path = output_dir / f"{stem}{suffix}"
            save_options = {"compression": "tiff_lzw"} if suffix == ".tiff" else {}
            raster.save(output_path, dpi=(dpi, dpi), **save_options)

            with Image.open(output_path) as saved:
                valid_raster = saved.size == raster_sizes[dpi] and saved.mode == mode
                saved.verify()
            if not valid_raster:
                raise ValueError(f"Unexpected raster format or dimensions for {output_path.name}")

    if not label_path.is_file() or label_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty label file for {stem}")
    plt.close(figure)


def main() -> None:
    """Export the three figures in supplement order. No return value."""
    revision_dir = Path(__file__).resolve().parents[1]
    output_dir = revision_dir / "supplementary_R1"
    output_dir.mkdir(parents=True, exist_ok=True)

    # The supplement uses its own typography and colors for the October print layout.
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.2,
            "axes.titlesize": 8.5,
            "axes.labelsize": 8.2,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "figure.titlesize": 11,
            "savefig.facecolor": "white",
            "savefig.bbox": None,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    blue = "#3C5488"
    red = "#E64B35"
    datasets = ("GSE174554", "GSE274546")

    # S1. Order paired patients by recurrent-minus-primary change within each cohort.
    recurrence_dir = revision_dir / "step1_primary_recurrence"
    patient_table = pd.read_csv(recurrence_dir / "patient_paired_changes.csv")
    effect_table = pd.read_csv(recurrence_dir / "patient_score_effects.csv")
    figure, axes = plt.subplots(2, 1, figsize=(6.69, 7.3), sharey=True)
    for panel_index, dataset in enumerate(datasets):
        axis = axes[panel_index]
        cohort_table = patient_table.loc[patient_table.dataset == dataset].sort_values("delta")
        effect = effect_table.loc[effect_table.dataset == dataset].iloc[0]
        colors = np.where(cohort_table.delta > 0, red, blue)
        axis.bar(np.arange(len(cohort_table)), cohort_table.delta, color=colors, width=0.9)

        # Read the paired mean, confidence limits and counts from the frozen summary table.
        mean_label = f"Mean {effect.mean_delta:.3f} [95% CI {effect.ci_low:.3f}, {effect.ci_high:.3f}]"
        axis.axhline(effect.mean_delta, color="black", lw=1.1, label=mean_label)
        axis.axhline(0, color="#888", lw=0.7)
        axis.legend(frameon=False, fontsize=8.2)
        axis.set(xlabel="Paired patient (ordered by change)", ylabel="Recurrent − primary score")
        title = f"{'AB'[panel_index]}  {dataset}: {int(effect.positive_pairs)}/{int(effect.n_pairs)} increased"
        axis.set_title(title, loc="left", fontweight="bold", pad=7)

    figure.tight_layout(rect=[0, 0, 1, 0.98])
    write_supplementary_figure(figure, output_dir, "Supplementary_Figure_S1")

    # S2A–B. Draw both random-set distributions on the fixed 170 × 218 mm page.
    figure, axes = plt.subplots(
        3, 1,
        figsize = (figure_tools.mm(170), figure_tools.mm(218)),
        gridspec_kw = {"height_ratios": [1, 1, 1.08]},
    )
    random_dir = revision_dir / "step2_extended_specificity"
    random_summary_table = pd.read_csv(random_dir / "random_summary.csv")
    for panel_index, dataset in enumerate(datasets):
        axis = axes[panel_index]
        for scheme, color in (("size_only", blue), ("size_and_expression", red)):
            draw_table = pd.read_csv(random_dir / f"{dataset}_{scheme}_random_sets.csv")

            # Keep the original fallback to a column name that contains "es".
            if "ES" in draw_table.columns:
                score_column = "ES"
            elif "es" in draw_table.columns:
                score_column = "es"
            else:
                score_column = next(column for column in draw_table.columns if "es" in column.lower())
            axis.hist(
                draw_table[score_column],
                bins     = 45,
                histtype = "step",
                lw       = 1.2,
                color    = color,
                label    = scheme.replace("_", " "),
            )

        observed_es = random_summary_table.loc[random_summary_table.dataset == dataset, "observed_es"].iloc[0]
        axis.axvline(observed_es, color="black", lw=1, label=f"Observed ES={observed_es:.3f}")
        axis.legend(frameon=False, fontsize=8.2)
        axis.set(xlabel="Random-set weighted ES", ylabel="Sets")
        axis.set_title(f"{'AB'[panel_index]}  {dataset}", loc="left", fontweight="bold", pad=7)

    # S2C. Use the estimated-correlation camera rows in their original gene-set order.
    camera_table = pd.read_csv(revision_dir / "step2_hallmark_completion/camera_seven_set_family.csv")
    estimated_table = camera_table.loc[camera_table["mode"] == "estimated"].copy()
    set_ids = estimated_table.set_id.drop_duplicates().tolist()
    set_positions = np.arange(len(set_ids))
    hallmark_labels = []
    axis = axes[2]
    for cohort_index, (dataset, color) in enumerate(zip(datasets, (blue, red))):
        cohort_table = estimated_table.loc[estimated_table.dataset == dataset].set_index("set_id").loc[set_ids]
        direction = np.where(cohort_table.direction.astype(str).str.lower().eq("up"), 1, -1)
        height = direction * (-np.log10(cohort_table.p_value.to_numpy(float)))
        axis.bar(set_positions + (cohort_index - 0.5) * 0.36, height, width=0.36, label=dataset, color=color)

        hallmark = cohort_table.loc["HALLMARK_INFLAMMATORY_RESPONSE"]
        hallmark_labels.append(f"{dataset} {hallmark.direction} (P={hallmark.p_value:.3f})")

    # Dashed lines mark the original nominal P = 0.05 guide for each direction.
    threshold = -np.log10(0.05)
    axis.axhline(0, color="#333333", lw=0.8)
    axis.axhline(threshold, color="#777777", ls="--", lw=0.8)
    axis.axhline(-threshold, color="#777777", ls="--", lw=0.8)

    short_labels = {
        "target_top20": "Top 20",
        "target_top50": "Top 50",
        "target_top100": "Top 100",
        "Abdelfattah_iMic_top50": "i-Mic",
        "Sankowski_C1": "C1",
        "Miller_systemic_top20": "Systemic",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Hallmark",
    }
    axis.set_xticks(set_positions, [short_labels.get(set_id, set_id) for set_id in set_ids], rotation=25, ha="right")
    axis.set_ylabel("Signed −log10(camera P)\nUp > 0; Down < 0")
    axis.legend(frameon=False, ncol=2, fontsize=8.2)
    axis.text(
        0.02, 0.07, "Hallmark: " + "; ".join(hallmark_labels),
        transform = axis.transAxes,
        fontsize  = 8.2,
        va        = "bottom",
        color     = "#333333",
    )
    axis.set_title("C  Seven-set estimated-correlation camera comparison", loc="left", fontweight="bold", pad=7)
    figure.tight_layout(pad=0.8, h_pad=1.1)
    figure.subplots_adjust(left=0.14)
    write_supplementary_figure(figure, output_dir, "Supplementary_Figure_S2", fixed_canvas=True)

    # S3A. Plot the frozen paired changes and CIs for broad myeloid groups.
    figure, axes = plt.subplots(2, 1, figsize=(6.69, 8.3), gridspec_kw={"height_ratios": [0.8, 2.5]})
    state_dir = revision_dir / "step3_state_analysis"
    abundance_table = pd.read_csv(state_dir / "state_abundance_statistics.csv")
    broad_table = abundance_table.loc[abundance_table.level == "broad"].sort_values(["dataset", "state"])
    state_labels = broad_table.state.str.replace("MAC_MDSC-like", "MAC/MDSC-like", regex=False)
    state_labels = broad_table.dataset + " / " + state_labels
    y_positions = np.arange(len(broad_table))[::-1]
    ci_lower = broad_table.mean_delta_pp - broad_table.ci_low_pp
    ci_upper = broad_table.ci_high_pp - broad_table.mean_delta_pp
    axis = axes[0]
    axis.errorbar(
        broad_table.mean_delta_pp, y_positions,
        xerr    = [ci_lower, ci_upper],
        fmt     = "o",
        ms      = 4,
        capsize = 2,
        color   = blue,
        lw      = 1.2,
    )
    axis.axvline(0, color="#999999", lw=0.7)
    axis.set_yticks(y_positions, state_labels.astype(str))
    axis.grid(axis="x", alpha=0.15)
    axis.set_xlabel("Recurrent − primary percentage points")
    axis.set_title("A  Paired broad myeloid-group fractions", loc="left", fontweight="bold", pad=7)

    # S3B. Show coverage at the formal 20-cell threshold, colored by test eligibility.
    camera_table = pd.read_csv(state_dir / "state_program_camera.csv")
    covered = (camera_table.threshold == 20) & (camera_table.correlation_mode == "estimated")
    coverage_table = camera_table.loc[covered].copy()
    state_labels = coverage_table.state.str.replace("MAC_MDSC-like", "MAC/MDSC-like", regex=False)
    coverage_table["key"] = coverage_table.dataset + " / " + state_labels
    coverage_table = coverage_table.sort_values(["dataset", "level", "state"])
    y_positions = np.arange(len(coverage_table))[::-1]
    axis = axes[1]
    axis.barh(y_positions, coverage_table.n_pairs, color=np.where(coverage_table.status == "analyzed", red, blue))
    axis.set_yticks(y_positions, coverage_table.key, fontsize=8.2)
    axis.set_xlabel("Eligible paired patients (≥20 cells each time point)")
    axis.set_title(
        "B  State-paired coverage; red = program test eligible (≥16/20 genes)",
        loc = "left",
        fontweight = "bold",
        pad = 7,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.98])
    write_supplementary_figure(figure, output_dir, "Supplementary_Figure_S3")


if __name__ == "__main__":
    main()
