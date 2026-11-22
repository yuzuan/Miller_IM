#!/usr/bin/env python3
"""Render Figure 5 from the archived exposure timelines, ranks, and patient models.

Inputs
    Under NOA_LEGACY_ROOT, use write/41_mg_inflammatory_sci_rebuild/Figure5/source_data/
    figure5_treatment_timeline.csv and figure5d_patient_score_samples.csv. In the parent Figure5 folder, use
    GSE121810_raw20_score_models.csv and GSE154795_raw20_score_models_auto_strict_idhwt.csv.
    Use write/56_figure5_pd1_rebuild/Figure5/frozen_inputs/GSE154795_myeloid_sample_audit.csv
    for the GSE154795 covariates, and the following R1 tables for the rank curves and annotations:
    step5_rank_unification_v2/ranked_genes.csv, step5_rank_unification_v2/sqrtF_fgsea.csv,
    and step2_pd1_camera/camera_results.csv.
Outputs
    Under _figure_build/56_figure5_pd1_rebuild/Figure5, write Figure5_panel_manifest.csv,
    Figure5_legend.md, FINAL_RESULT.md and source_data/:
        Figure5A_anti_pd1_exposure_timeline_source.csv, Figure5B_dual_cohort_raw20_enrichment_source.csv.gz,
        Figure5C_patient_scores_and_adjusted_effects_source.csv, Figure5_pd1_preview_layout.csv.
    Under _panels/56_figure5_pd1_rebuild/Figure5, write PDF masters and 600 dpi PNGs for all three panels,
    plus Figure5_pd1_preview.pdf, Figure5_pd1_preview.png and Figure5_pd1_preview_grayscale.png.
Run
    python scripts/08e_figure5.py from NOA_R1, or use scripts/09b_export_final_figures.py for the final figure.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps


figure_tools = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
project_root = figure_tools.legacy_root(revision_root)

# The final exporter rebinds these directories for its isolated build.
WRITE_ROOT = revision_root / "_figure_build/56_figure5_pd1_rebuild/Figure5"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "_panels/56_figure5_pd1_rebuild/Figure5"

# The original exposure color is shared across the timeline, enrichment and adjusted effects.
exposed = "#E64B35"
text_color = "#222222"
grid_color = "#E8E8E8"


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
) -> dict[str, str]:
    """Return the original per-panel manifest fields from the three written files."""
    for path in (source, pdf, png):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 5 panel output: {path.name}")
    panel_record = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": figure_tools.sha256(source),
        "pdf_sha256": figure_tools.sha256(pdf),
        "png_sha256": figure_tools.sha256(png),
    }
    return panel_record


def panel_a() -> dict[str, str]:
    """Write and draw the treatment events and resections for the four comparison groups."""
    timeline_path = project_root / (
        "write/41_mg_inflammatory_sci_rebuild/Figure5/source_data/figure5_treatment_timeline.csv"
    )
    timeline_table = pd.read_csv(timeline_path)
    if len(timeline_table) != 7:
        raise ValueError("Figure 5A needs seven timeline events")
    expected = {
        "GSE121810 neo-adjuvant (n=14)",
        "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)",
        "GSE154795 untreated recurrent (n=8)",
    }
    if set(timeline_table["lane_label"]) != expected:
        raise ValueError("Unexpected cohort labels or sample sizes for the Figure 5A timeline")

    # Apply the October display labels before writing the same four groups to the panel source table.
    timeline_table["lane_label"] = timeline_table["lane_label"].replace(
        {
            "GSE121810 neo-adjuvant (n=14)":
                "GSE121810 neoadjuvant (n=14)",
            "GSE154795 untreated recurrent (n=8)":
                "GSE154795 PD-1-unexposed recurrent (n=8)",
        }
    )
    stem = "Figure5A_anti_pd1_exposure_timeline"
    source = SOURCE_OUT / f"{stem}_source.csv"
    timeline_table.to_csv(source, index=False)

    # Use one lane per frozen comparison group and a common resection x coordinate.
    lane_order = [
        "GSE121810 neoadjuvant (n=14)",
        "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)",
        "GSE154795 PD-1-unexposed recurrent (n=8)",
    ]
    y_lookup = {label: len(lane_order) - 1 - index for index, label in enumerate(lane_order)}

    fig = figure_tools.new_figure(178, 50)
    ax = fig.add_axes([0.31, 0.20, 0.66, 0.60])
    for label in lane_order:
        row_y = y_lookup[label]
        ax.hlines(row_y, 0, 2, color="#B9B9B9", linewidth=1.4, zorder=1)
        ax.scatter(1, row_y, marker="D", s=32, color=text_color, edgecolor="white", linewidth=0.5, zorder=4)

    # Exposure starts before resection for step 1 and after resection for the other exposure event.
    exposure_table = timeline_table.loc[timeline_table["event_group"].eq("Exposure")]
    for row in exposure_table.itertuples():
        row_y = y_lookup[row.lane_label]
        if row.step == 1:
            exposure_x, arrow_start, arrow_end = 0.24, 0.31, 0.94
        else:
            exposure_x, arrow_start, arrow_end = 1.76, 1.06, 1.69
        ax.scatter(exposure_x, row_y, marker="o", s=42, color=exposed, edgecolor="white", linewidth=0.55, zorder=5)
        ax.annotate(
            "",
            xy     = (arrow_end, row_y),
            xytext = (arrow_start, row_y),
            arrowprops = {"arrowstyle": "-|>", "color": exposed, "linewidth": 1.3},
            zorder = 3,
        )

    # Keep the event labels at their frozen timeline coordinates.
    ax.text(0.24, 3.32, "Pembrolizumab", ha="center", va="bottom", fontsize=6.1)
    ax.text(1.76, 2.32, "Starts after surgery", ha="center", va="bottom", fontsize=6.1)
    ax.text(0.24, 1.32, "Anti-PD-1", ha="center", va="bottom", fontsize=6.1)
    ax.text(1.00, 2.65, "Resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, 1.65, "Resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, 0.65, "Recurrent resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, -0.35, "Recurrent resection", ha="center", va="top", fontsize=5.7)

    ax.set_yticks([3, 2, 1, 0], lane_order, fontsize=6.3)
    ax.set_xticks([0.24, 1.0, 1.76], ["Before surgery", "Resection", "After surgery"], fontsize=6.1)
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.65, 3.55)
    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(
        0.04,
        0.955,
        "Anti-PD-1 treatment and tissue sampling",
        fontsize   = 9.0,
        fontweight = "bold",
        va         = "top",
    )
    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step56 recurrent GBM Figure 5")
    panel_record = record("A", stem, "Treatment exposure timing in the two Figure5 cohorts", source, pdf, png)
    return panel_record


def sci_label(value: float) -> str:
    """Format the plotted p-value using the original 0.001 display threshold."""
    if value >= 0.001:
        label = f"{value:.3f}"
        return label
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / (10 ** exponent)
    label = rf"${coefficient:.1f}\times10^{{{exponent}}}$"
    return label


def panel_b() -> dict[str, str]:
    """Draw signed square-root F enrichment curves with the archived fgsea and camera results."""
    ranked_table = pd.read_csv(revision_root / "step5_rank_unification_v2/ranked_genes.csv")
    fgsea_table = pd.read_csv(revision_root / "step5_rank_unification_v2/sqrtF_fgsea.csv")
    fgsea_table = fgsea_table.set_index("dataset")
    camera_table = pd.read_csv(revision_root / "step2_pd1_camera/camera_results.csv")

    # Use the prespecified 20 genes and preserve the observed coverage of 19 genes in both cohorts.
    miller_im_genes = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
        "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]
    curve_tables = {}
    cohort_counts = {}
    for cohort in ("GSE121810", "GSE154795"):
        rank_columns = ["gene", "logFC", "original_F", "rank_sqrtF", "target"]
        cohort_rank_table = ranked_table.loc[ranked_table["dataset"].eq(cohort), rank_columns].copy()
        cohort_rank_table = cohort_rank_table.rename(columns={"original_F": "F"})
        statistics = fgsea_table.loc[cohort]
        n_present = int(cohort_rank_table["target"].sum())
        if len(cohort_rank_table) != int(statistics["n_tested"]) or n_present != int(statistics["size"]):
            raise ValueError(f"{cohort} ranked-gene counts differ from the frozen fgsea results")

        # Order ties by gene name and use absolute signed square-root F weights for hits.
        curve_table = cohort_rank_table[["gene", "logFC", "F"]].dropna().copy()
        curve_table["rank_stat"] = np.sign(curve_table["logFC"]) * np.sqrt(np.maximum(curve_table["F"], 0))
        curve_table = curve_table.sort_values(["rank_stat", "gene"], ascending=[False, True]).reset_index(drop=True)
        curve_table["rank"] = np.arange(1, len(curve_table) + 1)
        curve_table["rank_percentile"] = 100 * (curve_table["rank"] - 1) / max(1, len(curve_table) - 1)
        curve_table["hit"] = curve_table["gene"].isin(miller_im_genes)
        if n_present != 19 or int(curve_table["hit"].sum()) != n_present:
            raise ValueError(f"{cohort} must contain 19 of the 20 MG-Inflam20 genes")
        weights = curve_table["rank_stat"].abs().to_numpy(dtype=float)
        hits = curve_table["hit"].to_numpy(dtype=bool)
        hit_weight_sum = float(weights[hits].sum())
        increments = np.where(hits, weights / hit_weight_sum, -1.0 / int((~hits).sum()))
        curve_table["running_enrichment"] = np.cumsum(increments)

        # The curve is descriptive. Keep the archived fgsea and both camera test statistics.
        test_table = camera_table.loc[camera_table["dataset"].eq(cohort)].set_index("mode")
        curve_table["cohort"] = cohort
        curve_table["NES"] = float(statistics["NES"])
        curve_table["fgsea_nominal_p"] = float(statistics["p_value"])
        curve_table["camera_direction"] = "Up"
        curve_table["camera_nominal_p"] = float(test_table.loc["estimated", "p_value"])
        curve_table["camera_fixed_p"] = float(test_table.loc["fixed_0.01", "p_value"])
        curve_table["n_present"] = n_present
        curve_table["n_defined"] = len(miller_im_genes)
        curve_tables[cohort] = curve_table
        cohort_counts[cohort] = n_present

    stem = "Figure5B_dual_cohort_raw20_enrichment"
    curve_table_total = pd.concat(curve_tables.values(), ignore_index=True)
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    curve_table_total.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    # Keep the GSE121810 curve on the left and the GSE154795 curve on the right.
    fig = figure_tools.new_figure(178, 90)
    comparisons = [
        ("GSE121810", 0.08, "Pre-op exposed vs post-op only"),
        ("GSE154795", 0.56, "PD-1 exposed vs PD-1-unexposed recurrent"),
    ]
    for cohort, left, comparison in comparisons:
        curve_table = curve_tables[cohort]
        ax = fig.add_axes([left, 0.35, 0.40, 0.44])
        ax.axhline(0, color="#A5A5A5", linewidth=0.6)
        ax.fill_between(
            curve_table["rank_percentile"],
            0,
            curve_table["running_enrichment"],
            color     = "#F8DEDA",
            alpha     = 0.92,
            linewidth = 0,
            zorder    = 1,
        )
        ax.plot(
            curve_table["rank_percentile"],
            curve_table["running_enrichment"],
            color     = exposed,
            linewidth = 1.2,
            zorder    = 3,
        )
        hit_percentiles = curve_table.loc[curve_table["hit"], "rank_percentile"]
        ax.vlines(hit_percentiles, -0.105, -0.035, color=text_color, linewidth=0.65, zorder=4)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.13, 0.96)
        ax.set_xticks([0, 50, 100])
        ax.set_xlabel("Transcriptome rank percentile")
        if cohort == "GSE121810":
            ax.set_ylabel("Running enrichment score")
        else:
            ax.set_yticklabels([])
        ax.grid(axis="y", color=grid_color, linewidth=0.45)
        figure_tools.clean_axis(ax, keep_left=True, keep_bottom=True)
        ax.set_title(f"{cohort}\n{comparison}", fontsize=7.0, fontweight="bold", loc="left", pad=5)

        # The archived tests are constant within a cohort's rank table.
        statistics = curve_table.iloc[0]
        stat_text = (
            f"NES {statistics['NES']:.2f}; P {sci_label(float(statistics['fgsea_nominal_p']))}\n"
            f"camera (estimated) P {float(statistics['camera_nominal_p']):.3f}\n"
            f"camera (fixed 0.01) P {sci_label(float(statistics['camera_fixed_p']))}"
        )
        fig.text(
            left + 0.40,
            0.17,
            stat_text,
            ha          = "right",
            va          = "top",
            fontsize    = 5.7,
            color       = text_color,
            linespacing = 1.25,
        )

    fig.text(
        0.52, 0.225, f"{cohort_counts['GSE121810']}/{len(miller_im_genes)} genes measured in each cohort",
        ha       = "center",
        va       = "top",
        fontsize = 5.7,
        color    = text_color,
    )
    fig.text(
        0.04,
        0.965,
        "Preranked enrichment and correlation-aware tests",
        fontsize   = 9.0,
        fontweight = "bold",
        va         = "top",
    )
    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step56 recurrent GBM Figure 5")
    panel_record = record("B", stem, "Two-cohort MG-Inflam20 rank-enrichment evidence", source, pdf, png)
    return panel_record


def load_patient_evidence() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return patient scores and frozen adjusted and unadjusted effects after checking their cohort labels."""
    source_dir = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure5"
    score_table = pd.read_csv(source_dir / "source_data/figure5d_patient_score_samples.csv")
    gse121810_model_table = pd.read_csv(source_dir / "GSE121810_raw20_score_models.csv")
    gse154795_model_table = pd.read_csv(source_dir / "GSE154795_raw20_score_models_auto_strict_idhwt.csv")
    model_table_total = pd.concat([gse121810_model_table, gse154795_model_table], ignore_index=True, sort=False)
    if len(score_table) != 50 or set(score_table["cohort"]) != {"GSE121810", "GSE154795"}:
        raise ValueError("Figure 5C needs the 50 frozen patient scores from both cohorts")

    # Match the frozen GSE154795 sample metadata to the patient-score rows.
    audit_path = project_root / (
        "write/56_figure5_pd1_rebuild/Figure5/frozen_inputs/GSE154795_myeloid_sample_audit.csv"
    )
    covariate_table = pd.read_csv(audit_path)
    covariate_table = covariate_table.rename(columns={"ID": "sample", "condition": "audit_condition"})
    audit_columns = [
        "sample", "audit_condition", "IDH_status", "MGMT", "n_myeloid_cells", "GSM",
        "sequencing_method", "gel_bead_version", "acquisition_wave",
    ]
    if not set(audit_columns).issubset(covariate_table.columns):
        raise ValueError("The GSE154795 sample audit is missing model covariates")
    score_table = score_table.merge(
        covariate_table[audit_columns],
        on  = "sample",
        how = "left",
        validate = "many_to_one",
    )
    gse154795_table = score_table.loc[score_table["cohort"].eq("GSE154795")]
    covariates = ["IDH_status", "MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave"]
    if len(gse154795_table) != 21 or gse154795_table[covariates].isna().any().any():
        raise ValueError("All 21 GSE154795 patients must match the model covariates")
    if not gse154795_table["IDH_status"].eq("WT").all():
        raise ValueError("All GSE154795 patients in the frozen cohort must have IDH-wild-type tumors")
    if not gse154795_table["condition"].eq(gse154795_table["audit_condition"]).all():
        raise ValueError("GSE154795 patient groups differ between the score and covariate tables")

    # Retain the archived raw and adjusted models, which define each cohort's two-model FDR family.
    adjusted_table = model_table_total.loc[model_table_total["metric"].eq("raw20_score_myeloid_adjusted")].copy()
    raw_model_table = model_table_total.loc[model_table_total["metric"].eq("raw20_score")].copy()
    if len(adjusted_table) != 2 or len(raw_model_table) != 2:
        raise ValueError("Figure 5C needs two archived adjusted models and two raw-score models")
    expected = {
        "GSE121810": (
            0.416017028595708, 0.137781145186939, 0.694252912004477,
            0.00492100468020317, 0.00984200936040635,
        ),
        "GSE154795": (
            0.427369843346401, -0.222316569830202, 1.077056256523,
            0.180119103326531, 0.180119103326531,
        ),
    }
    for row in adjusted_table.itertuples():
        values = expected[row.cohort]
        observed = (row.effect, row.conf_low, row.conf_high, row.p_value, row.fdr)
        if not all(np.isclose(x, y, atol=1e-12) for x, y in zip(observed, values)):
            raise ValueError(f"{row.cohort} patient model differs from the frozen estimates")
    patient_evidence = (score_table, adjusted_table, raw_model_table)
    return patient_evidence


def draw_box_swarm(
    ax: mpl.axes.Axes,
    values_by_group: list[np.ndarray],
    labels: list[str],
    colors: list[str],
    seed: int,
) -> None:
    """Draw group distributions at fixed positions using the cohort seed. Returns nothing."""
    rng = np.random.default_rng(seed)
    for index, (values, color) in enumerate(zip(values_by_group, colors)):
        values = np.asarray(values, dtype=float)
        jitter = rng.uniform(-0.10, 0.10, len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s         = 18,
            color     = color,
            alpha     = 0.78,
            edgecolor = "white",
            linewidth = 0.4,
            zorder    = 3,
        )

        # Bound the 1.5-IQR whiskers by the observed values, matching the original display.
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        whisker_low = max(values.min(), q1 - 1.5 * (q3 - q1))
        whisker_high = min(values.max(), q3 + 1.5 * (q3 - q1))
        ax.vlines(index, whisker_low, whisker_high, color=text_color, linewidth=0.8, zorder=4)
        rectangle = mpl.patches.Rectangle(
            (index - 0.16, q1),
            0.32,
            q3 - q1,
            facecolor = "white",
            edgecolor = text_color,
            linewidth = 0.8,
            zorder    = 5,
        )
        ax.add_patch(rectangle)
        ax.hlines(median, index - 0.16, index + 0.16, color=text_color, linewidth=1.0, zorder=6)
    ax.axhline(0, color="#8A8A8A", linewidth=0.65, linestyle=(0, (3, 2)), zorder=1)
    ax.set_xticks(range(len(labels)), labels, fontsize=5.8)
    ax.set_xlim(-0.45, len(labels) - 0.55)
    ax.set_ylim(-1.65, 1.45)
    ax.grid(axis="y", color=grid_color, linewidth=0.45)
    figure_tools.clean_axis(ax, keep_left=True, keep_bottom=True)
    ax.tick_params(axis="x", length=0, pad=3)


def draw_adjusted_effect(ax: mpl.axes.Axes, row: pd.Series, show_xlabel: bool) -> None:
    """Draw one archived adjusted estimate and interval. Returns nothing."""
    ax.axvline(0, color="#8A8A8A", linewidth=0.7)
    ax.hlines(0, float(row["conf_low"]), float(row["conf_high"]), color=text_color, linewidth=1.1, zorder=2)
    significant = float(row["fdr"]) < 0.05
    ax.scatter(
        float(row["effect"]),
        0,
        marker    = "D",
        s         = 38,
        facecolor = exposed if significant else "white",
        edgecolor = exposed,
        linewidth = 0.9,
        zorder    = 4,
    )
    ax.text(
        0.98,
        0.96,
        f"FDR {float(row['fdr']):.3f}",
        transform  = ax.transAxes,
        ha         = "right",
        va         = "top",
        fontsize   = 6.1,
        color      = exposed if significant else text_color,
        fontweight = "bold" if significant else "normal",
    )
    ax.text(
        0.98,
        0.79,
        f"effect {float(row['effect']):+.3f}",
        transform = ax.transAxes,
        ha        = "right",
        va        = "top",
        fontsize  = 5.6,
        color     = "#777777",
    )
    ax.set_xlim(-0.40, 1.25)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0])
    if show_xlabel:
        ax.set_xlabel("Adjusted exposure-associated MG-Inflam20 score shift")
    else:
        ax.set_xticklabels([])
    ax.grid(axis="x", color=grid_color, linewidth=0.45)
    figure_tools.clean_axis(ax, keep_left=False, keep_bottom=True)


def panel_c() -> dict[str, str]:
    """Plot observed patient scores alongside the original adjusted effects for both cohorts."""
    score_table, model_table, raw_model_table = load_patient_evidence()

    # Record the original cohort-specific formula and the family used for the archived FDR.
    score_table = score_table.assign(
        analysis_group = score_table["group"].fillna(score_table["condition"]),
        model_formula  = score_table["cohort"].map(
            {
                "GSE121810": "raw20_score ~ myeloid_abundance + group",
                "GSE154795": (
                    "raw20_score ~ myeloid_abundance + MGMT + sequencing_method + "
                    "gel_bead_version + acquisition_wave + condition"
                ),
            }
        ),
        fdr_family     = "BH across raw and myeloid-adjusted patient-score models within cohort (n=2)",
    )
    source_table = score_table.merge(
        model_table[
            [
                "cohort", "effect", "conf_low", "conf_high", "p_value", "fdr",
                "n_positive", "n_negative",
            ]
        ],
        on  = "cohort",
        how = "left",
        validate = "many_to_one",
    )
    source_table = source_table.merge(
        raw_model_table[["cohort", "p_value", "fdr"]].rename(
            columns = {"p_value": "raw_model_p_value", "fdr": "raw_model_fdr"}
        ),
        on  = "cohort",
        how = "left",
        validate = "many_to_one",
    )
    stem = "Figure5C_patient_scores_and_adjusted_effects"
    source = SOURCE_OUT / f"{stem}_source.csv"
    source_table.to_csv(source, index=False)

    # Place the patient distributions to the left of the adjusted effects on each row.
    fig = figure_tools.new_figure(178, 78)
    fig.text(
        0.04,
        0.965,
        "Patient scores and adjusted associations",
        fontsize   = 9.0,
        fontweight = "bold",
        va         = "top",
    )

    configurations = [
        {
            "cohort": "GSE121810",
            "bottom": 0.53,
            "groups": ["Adj", "Neo"],
            "labels": ["Post-op only\n(n=15)", "Pre-op exposed\n(n=14)"],
            "title": "GSE121810 bulk tumors",
            "seed": 121810,
        },
        {
            "cohort": "GSE154795",
            "bottom": 0.13,
            "groups": ["GBM.rec", "GBM.PD1"],
            "labels": ["PD-1-unexposed recurrent\n(n=8)", "PD-1 exposed\n(n=13)"],
            "title": "GSE154795 IDH-wildtype myeloid pseudobulks",
            "seed": 154795,
        },
    ]
    model_index_table = model_table.set_index("cohort")
    for index, config in enumerate(configurations):
        cohort = config["cohort"]
        cohort_score_table = score_table.loc[score_table["cohort"].eq(cohort)].copy()
        group_column = "group" if cohort == "GSE121810" else "condition"
        values = [
            cohort_score_table.loc[cohort_score_table[group_column].eq(group), "raw20_score"].to_numpy(dtype=float)
            for group in config["groups"]
        ]
        expected_n = [15, 14] if cohort == "GSE121810" else [8, 13]
        if [len(item) for item in values] != expected_n:
            raise ValueError(f"Unexpected patient counts by exposure group for {cohort}")
        dist_ax = fig.add_axes([0.10, config["bottom"], 0.47, 0.27])
        draw_box_swarm(dist_ax, values, config["labels"], ["#3C5488", exposed], config["seed"])
        dist_ax.set_title(config["title"], fontsize=6.8, fontweight="bold", loc="left", pad=3)
        dist_ax.set_ylabel("MG-Inflam20 score")

        effect_ax = fig.add_axes([0.64, config["bottom"] + 0.01, 0.31, 0.24])
        draw_adjusted_effect(effect_ax, model_index_table.loc[cohort], show_xlabel=index == 1)

    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step56 recurrent GBM Figure 5")
    panel_record = record("C", stem, "MG-Inflam20 scores with cohort-specific adjusted model effects", source, pdf, png)
    return panel_record


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    """Create the original 600 dpi image previews and composite PDF."""
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 115, "top": 155, "right": 4205, "bottom": 1260},
        {"panel": "B", "left": 115, "top": 1375, "right": 4205, "bottom": 3415},
        {"panel": "C", "left": 115, "top": 3530, "right": 4205, "bottom": 5310},
    ]
    canvas = Image.new("RGB", (4320, 5425), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = figure_tools.image_font(76, bold=True)
    for item in layout:
        panel = item["panel"]
        figure_tools.paste_contained(
            canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"])
        )
        drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

    png = FIG_OUT / "Figure5_pd1_preview.png"
    gray = FIG_OUT / "Figure5_pd1_preview_grayscale.png"
    pdf = FIG_OUT / "Figure5_pd1_preview.pdf"
    source = SOURCE_OUT / "Figure5_pd1_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(600, 600), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(pdf, format="PDF", resolution=600.0)
    layout_table = pd.DataFrame(layout)
    layout_table.to_csv(source, index=False)
    output_paths = (pdf, png, gray, source)
    return output_paths


def main() -> None:
    """Render the three Figure 5 panels and write their preview and source manifest."""
    figure_tools.apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    records = [panel_a(), panel_b(), panel_c()]
    if [row["panel"] for row in records] != list("ABC"):
        raise ValueError("Figure 5 requires panels A through C in order")

    # Write one manifest row per panel with the existing preview and layout fields.
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)
    manifest_table["preview_layout"] = str(preview_layout)
    manifest_path = WRITE_ROOT / "Figure5_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)

    # Keep the October Figure 5 legend so the figure can build from the code archive alone.
    legend_path = WRITE_ROOT / "Figure5_legend.md"
    legend_text = (
        "# Figure 5. Exploratory anti-PD-1 exposure associations\n\n"
        "A, GSE121810 compared preoperative/adjuvant-only pembrolizumab (14/15 patients); GSE154795 compared "
        "PD-1-exposed/unexposed myeloid pseudobulks (13/8). B, Preranked enrichment curves annotated with "
        "estimated-correlation camera (primary) and fixed-correlation camera (sensitivity) results. "
        "C, Patient-score distributions and adjusted exposure-associated differences. Points represent patients; "
        "boxes show medians and interquartile ranges, with whiskers bounded by 1.5 interquartile ranges and "
        "observed values. Diamonds and bars show adjusted estimates and 95% CIs (filled, FDR < 0.05). Cohorts "
        "may overlap. Abbreviations: CI, confidence interval; FDR, false discovery rate; IDH, isocitrate "
        "dehydrogenase; MG-Inflam20, 20-gene representation of the microglial inflammatory program; "
        "NES, normalized enrichment score; PD-1, programmed cell death protein 1.\n\n"
        "Alt text: Sampling timelines, nonsignificant correlation-aware tests, and adjusted score differences "
        "significant only in the first cohort.\n"
    )
    legend_path.write_text(legend_text, encoding="utf-8")

    # Keep the existing build summary and verify the written composite and its references.
    result_path = WRITE_ROOT / "FINAL_RESULT.md"
    result_text = (
        "# Figure 5 anti-PD-1 exposure rebuild\n\n"
        "- A. Anti-PD-1 exposure timing at sampling in the two cohorts.\n"
        "- B. MG-Inflam20 enrichment curves with archived fgsea and camera results.\n"
        "- C. Observed patient scores and cohort-specific adjusted associations.\n"
        f"- Manifest `{manifest_path}`\n"
        f"- Legend `{legend_path}`\n"
    )
    result_path.write_text(result_text, encoding="utf-8")
    for path in (preview_pdf, preview_png, preview_gray, preview_layout, manifest_path, legend_path, result_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 5 output: {path.name}")
    print(f"Figure 5 panels: {len(records)}")


if __name__ == "__main__":
    main()
