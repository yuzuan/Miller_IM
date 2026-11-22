#!/usr/bin/env python3
# Purpose  Draw the Figure 5 exposure timelines, rank enrichment, scores and adjusted effects.
# Inputs   EdgeR, fgsea, camera, coverage and model tables under write/41_mg_inflammatory_sci_rebuild/Figure5,
#          and GSE154795_myeloid_sample_audit.csv under Step 56 frozen_inputs.
# Outputs  Source tables, Figure5_panel_manifest.csv, Figure5_legend.md and FINAL_RESULT.md under
#          write/56_figure5_pd1_rebuild/Figure5. The figures tree contains the three panels and previews.
# Run      python scripts/09b_figure5.py from the repository root.
# Use the existing Python libraries to preserve ordering, random draws and the canvas.
# Keep frozen fields, formulas and results. Render PNG panels from PDFs and use PIL for previews.

# Each panel below has a PDF, PNG and matching source CSV in source_data.
#   Figure5A_anti_pd1_exposure_timeline
#   Figure5B_dual_cohort_raw20_enrichment
#   Figure5C_patient_scores_and_adjusted_effects
# Panel B uses CSV.GZ. The preview has a grayscale PNG and a layout CSV with the same stem.

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
image_font = figure_style.image_font
new_figure = figure_style.new_figure
panel_record = figure_style.panel_record
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel

# Set the original global seed before analysis. Keep the separate default_rng jitter seeds for both cohorts.
apply_publication_style()
np.random.seed(20260714)
root = Path(__file__).resolve().parents[1]
write_root = root / "write/56_figure5_pd1_rebuild/Figure5"
source_out = write_root / "source_data"
source_out.mkdir(parents=True, exist_ok=True)
figure_out = root / "figures/56_figure5_pd1_rebuild/Figure5"
figure_out.mkdir(parents=True, exist_ok=True)

# Each panel A row is a treatment event. Check seven events, four sequences and frozen cohort counts.
step41_root = root / "write/41_mg_inflammatory_sci_rebuild/Figure5"
step41_source = step41_root / "source_data"

timeline_table = pd.read_csv(step41_source / "figure5_treatment_timeline.csv")
if len(timeline_table) != 7:
    raise ValueError("Figure 5A requires four treatment sequences and seven events.")
expected = {
    "GSE121810 neo-adjuvant (n=14)",
    "GSE121810 adjuvant-only (n=15)",
    "GSE154795 PD-1 exposed (n=13)",
    "GSE154795 untreated recurrent (n=8)",
}
if set(timeline_table["lane_label"]) != expected:
    raise ValueError("The cohort labels or sample sizes in Figure 5A have changed.")
timeline_table["lane_label"] = timeline_table["lane_label"].replace(
    {
        "GSE154795 untreated recurrent (n=8)": "GSE154795 PD-1-unexposed recurrent (n=8)",
    }
)


records = []

stem = "Figure5A_anti_pd1_exposure_timeline"

# Write the panel source table in the original field and row order.
source_path = source_out / f"{stem}_source.csv"
timeline_table.to_csv(source_path, index=False)

lane_order = [
    "GSE121810 neo-adjuvant (n=14)", "GSE121810 adjuvant-only (n=15)", "GSE154795 PD-1 exposed (n=13)",
    "GSE154795 PD-1-unexposed recurrent (n=8)"
]
y_lookup = {label: 3 - index for index, label in enumerate(lane_order)}

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(178, 50)
axis = figure.add_axes([0.31, 0.20, 0.66, 0.60])
for label in lane_order:
    y_position = y_lookup[label]
    axis.hlines(y_position, 0, 2, color="#B9B9B9", linewidth=1.4, zorder=1)
    axis.scatter(1, y_position, marker="D", s=32, color="#222222", edgecolor="white", linewidth=0.5, zorder=4)


# Draw exposure arrows from the original events in the direction of surgery and treatment.
exposure_table = timeline_table.loc[timeline_table["event_group"].eq("Exposure")]
for row in exposure_table.itertuples():
    y_position = y_lookup[row.lane_label]
    if row.step == 1:
        exposure_x, arrow_start, arrow_end = 0.24, 0.31, 0.94
    else:
        exposure_x, arrow_start, arrow_end = 1.76, 1.06, 1.69
    axis.scatter(exposure_x, y_position, marker="o", s=42, color="#E64B35", edgecolor="white", linewidth=0.55, zorder=5)
    axis.annotate(
        "",
        xy     = (arrow_end, y_position),
        xytext = (arrow_start, y_position),
        arrowprops = {"arrowstyle": "-|>", "color": "#E64B35", "linewidth": 1.3},
        zorder = 3,
    )


# Keep the original exposure timing text and comparison groups.
axis.text(0.24, 3.32, "Pembrolizumab", ha="center", va="bottom", fontsize=6.1)
axis.text(1.76, 2.32, "Starts after surgery", ha="center", va="bottom", fontsize=6.1)
axis.text(0.24, 1.32, "Anti-PD-1", ha="center", va="bottom", fontsize=6.1)
axis.text(1.00, 2.65, "Resection", ha="center", va="top", fontsize=5.7)
axis.text(1.00, 1.65, "Resection", ha="center", va="top", fontsize=5.7)
axis.text(1.00, 0.65, "Recurrent resection", ha="center", va="top", fontsize=5.7)
axis.text(1.00, -0.35, "Recurrent resection", ha="center", va="top", fontsize=5.7)

axis.set_yticks([3, 2, 1, 0], lane_order, fontsize=6.3)
axis.set_xticks([0.24, 1.0, 1.76], ["Before surgery", "Resection", "After surgery"], fontsize=6.1)

# Use the original axis and labels for the comparison.
axis.set_xlim(0, 2)

# Keep the original axis limits and data order.
axis.set_ylim(-0.65, 3.55)
axis.tick_params(axis="both", length=0, pad=3)
for spine in axis.spines.values():
    spine.set_visible(False)
figure.text(
    0.04,
    0.955,
    "Anti-PD-1 exposure precedes tissue collection in the analyzed exposed groups",
    fontsize   = 9.0,
    fontweight = "bold",
    va         = "top",
)

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step56 recurrent GBM Figure 5")

records.append(
    panel_record("A", stem, "Treatment exposure timing in the two Figure5 cohorts", source_path, pdf_path, png_path)
)


# Read both edgeR tables and rank by F multiplied by the log-fold-change sign for panel B.


gse121810_differential_table = pd.read_csv(step41_root / "GSE121810_raw20_edger.csv")
gse154795_differential_table = pd.read_csv(step41_root / "GSE154795_raw20_edger_auto_strict_idhwt.csv")
differential_tables = {
    "GSE121810": gse121810_differential_table,
    "GSE154795": gse154795_differential_table,
}


# Combine the original enrichment, competitive-test and coverage summaries with one row per cohort.
enrichment_records = []
for cohort, suffix in [
    ("GSE121810", ""),
    ("GSE154795", "_auto_strict_idhwt"),
]:
    fgsea_table = pd.read_csv(step41_root / f"{cohort}_raw20_fgsea{suffix}.csv")
    fgsea_record = fgsea_table.iloc[0]
    camera_table = pd.read_csv(step41_root / f"{cohort}_raw20_camera{suffix}.csv")
    camera_record = camera_table.iloc[0]
    coverage_table = pd.read_csv(step41_root / f"{cohort}_raw20_coverage{suffix}.csv")
    coverage_record = coverage_table.iloc[0]
    enrichment_records.append(
        {
            "cohort": cohort,
            "NES": float(fgsea_record["NES"]),
            "fgsea_nominal_p": float(fgsea_record["pval"]),
            "camera_direction": str(camera_record["Direction"]),
            "camera_nominal_p": float(camera_record["PValue"]),
            "n_defined": int(coverage_record["n_defined"]),
            "n_present": int(coverage_record["n_present"]),
            "leading_edge": str(fgsea_record["leadingEdge"]),
        }
    )
enrichment_statistics_table = pd.DataFrame(enrichment_records)
expected = {
    "GSE121810": (1.89995642058625, 6.3505755649683e-05, 2.81295477233785e-05),
    "GSE154795": (1.75124003154227, 0.000290544571159221, 0.015367499497791),
}
for row in enrichment_statistics_table.itertuples():
    values = expected[row.cohort]
    if not (
        np.isclose(row.NES, values[0], atol=1e-12)
        and np.isclose(row.fgsea_nominal_p, values[1], atol=1e-15)
        and np.isclose(row.camera_nominal_p, values[2], atol=1e-15)
        and row.n_present == 19
        and row.n_defined == 20
    ):
        raise ValueError(f"The Step 41 program-level statistics for {row.cohort} have changed.")


# Retain the frozen program order, weighting, denominators and cumulative rank sequence.
miller_im_genes = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3"
]
enrichment_index_table = enrichment_statistics_table.set_index("cohort")
curve_tables = {}
for cohort, differential_table in differential_tables.items():
    enrichment_record = enrichment_index_table.loc[cohort]

    required = {"gene", "logFC", "F"}
    if not required.issubset(differential_table.columns):
        raise ValueError(f"{cohort} edgeR table lacks columns {sorted(required - set(differential_table.columns))}.")
    curve_table = differential_table[["gene", "logFC", "F"]].dropna().copy()
    curve_table["rank_stat"] = curve_table["F"] * np.sign(curve_table["logFC"])
    curve_table = curve_table.sort_values(["rank_stat", "gene"], ascending=[False, True]).reset_index(drop=True)
    curve_table["rank"] = np.arange(1, len(curve_table) + 1)
    curve_table["rank_percentile"] = 100 * (curve_table["rank"] - 1) / max(1, len(curve_table) - 1)
    curve_table["hit"] = curve_table["gene"].isin(miller_im_genes)
    if int(curve_table["hit"].sum()) != 19:
        raise ValueError(f"Expected coverage of 19/20 Miller-IM genes in {cohort}.")

    # Normalize hit weights using absolute rank statistics. Use the original miss step and accumulation order.
    weights = curve_table["rank_stat"].abs().to_numpy(dtype=float)
    hits    = curve_table["hit"].to_numpy(dtype=bool)
    hit_weight_sum = float(weights[hits].sum())
    increments = np.where(hits, weights / hit_weight_sum, -1.0 / int((~hits).sum()))
    curve_table["running_enrichment"] = np.cumsum(increments)

    # Keep source statistics and coverage fields for review of panel B.
    curve_table["cohort"] = cohort
    curve_table["NES"] = float(enrichment_record["NES"])
    curve_table["fgsea_nominal_p"] = float(enrichment_record["fgsea_nominal_p"])
    curve_table["camera_direction"] = str(enrichment_record["camera_direction"])
    curve_table["camera_nominal_p"] = float(enrichment_record["camera_nominal_p"])
    curve_table["n_present"] = int(enrichment_record["n_present"])
    curve_table["n_defined"] = int(enrichment_record["n_defined"])
    curve_tables[cohort] = curve_table

curve_table_total = pd.concat(curve_tables.values(), ignore_index=True)
stem = "Figure5B_dual_cohort_raw20_enrichment"

# Write the panel source table in the original field and row order.
source_path = source_out / f"{stem}_source.csv.gz"
curve_table_total.to_csv(source_path, index=False, compression={"method": "gzip", "mtime": 0})

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(178, 78)
positions = {
    "GSE121810": 0.08,
    "GSE154795": 0.56,
}
for cohort in ["GSE121810", "GSE154795"]:
    left = positions[cohort]
    curve_table = curve_tables[cohort]
    enrichment_record = enrichment_index_table.loc[cohort]
    axis = figure.add_axes([left, 0.28, 0.40, 0.56])
    axis.axhline(0, color="#A5A5A5", linewidth=0.6)
    axis.fill_between(
        curve_table["rank_percentile"],
        0,
        curve_table["running_enrichment"],
        color     = "#F8DEDA",
        alpha     = 0.92,
        linewidth = 0,
        zorder    = 1,
    )
    axis.plot(
        curve_table["rank_percentile"],
        curve_table["running_enrichment"],
        color     = "#E64B35",
        linewidth = 1.2,
        zorder    = 3,
    )

    # Mark the positions of the fixed program genes in the ranking.
    hit_table = curve_table.loc[curve_table["hit"]]
    axis.vlines(hit_table["rank_percentile"], -0.105, -0.035, color="#222222", linewidth=0.65, zorder=4)

    # Use the original axis and labels for the comparison.
    axis.set_xlim(0, 100)

    # Keep the original axis limits and data order.
    axis.set_ylim(-0.13, 0.96)
    axis.set_xticks([0, 50, 100])
    axis.set_xlabel("Transcriptome rank percentile")
    if cohort == "GSE121810":
        axis.set_ylabel("Running enrichment score")
    else:
        axis.set_yticklabels([])
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.45)
    clean_axis(axis, keep_left=True, keep_bottom=True)
    comparison = {
        "GSE121810": "Pre-op exposed vs post-op only",
        "GSE154795": "PD-1 exposed vs PD-1-unexposed recurrent",
    }
    axis.set_title(f"{cohort}\n{comparison[cohort]}", fontsize=7.0, fontweight="bold", loc="left", pad=5)

    # Format the two p values in scientific notation after all decisions.
    p_labels = []
    for p_value in [float(enrichment_record["fgsea_nominal_p"]), float(enrichment_record["camera_nominal_p"])]:
        if p_value >= 0.001:
            p_label = f"{p_value:.3f}"
        else:
            exponent = int(np.floor(np.log10(p_value)))
            coefficient = p_value / (10**exponent)
            p_label = rf"${coefficient:.1f}\times10^{{{exponent}}}$"
        p_labels.append(p_label)

    stat_text = (
        "19/20 genes measured\n"
        f"NES {enrichment_record['NES']:.2f}; P {p_labels[0]}\n"
        f"camera {enrichment_record['camera_direction']}; P {p_labels[1]}"
    )
    axis.text(
        0.98,
        0.96,
        stat_text,
        transform   = axis.transAxes,
        ha          = "right",
        va          = "top",
        fontsize    = 5.7,
        color       = "#222222",
        linespacing = 1.25,
    )

figure.text(
    0.04,
    0.965,
    "The Miller-IM program is enriched toward anti-PD-1 exposure in both cohorts",
    fontsize   = 9.0,
    fontweight = "bold",
    va         = "top",
)

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step56 recurrent GBM Figure 5")

records.append(
    panel_record("B", stem, "Two-cohort Miller-IM rank-enrichment evidence", source_path, pdf_path, png_path)
)


# Panel C uses one row per patient and combines the saved model and FDR results from the two cohorts.

patient_table_total = pd.read_csv(step41_source / "figure5d_patient_score_samples.csv")
gse121810_model_table = pd.read_csv(step41_root / "GSE121810_raw20_score_models.csv")
gse154795_model_table = pd.read_csv(step41_root / "GSE154795_raw20_score_models_auto_strict_idhwt.csv")
model_table_total = pd.concat([gse121810_model_table, gse154795_model_table], ignore_index=True, sort=False)
if len(patient_table_total) != 50 or set(patient_table_total["cohort"]) != {"GSE121810", "GSE154795"}:
    raise ValueError("Figure 5C patient-score source data must contain 29+21 samples.")

# Match covariates from the frozen audit table by the original fields and many-to-one join.
covariate_path = write_root / "frozen_inputs/GSE154795_myeloid_sample_audit.csv"
covariate_table = pd.read_csv(covariate_path)
covariate_table = covariate_table.rename(columns={"ID": "sample", "condition": "audit_condition"})
audit_columns = [
    "sample", "audit_condition", "IDH_status", "MGMT", "n_myeloid_cells", "GSM", "sequencing_method",
    "gel_bead_version", "acquisition_wave"
]
if not set(audit_columns).issubset(covariate_table.columns):
    raise ValueError("The GSE154795 patient covariate audit table is missing model columns.")
patient_table_total = patient_table_total.merge(
    covariate_table[audit_columns],
    on  = "sample",
    how = "left",
    validate = "many_to_one",
)
gse154795_table = patient_table_total.loc[patient_table_total["cohort"].eq("GSE154795")]
covariates = ["IDH_status", "MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave"]
if len(gse154795_table) != 21 or gse154795_table[covariates].isna().any().any():
    raise ValueError("All 21 GSE154795 patients must match the model covariates.")
if not gse154795_table["IDH_status"].eq("WT").all():
    raise ValueError("All patients in the formal GSE154795 analysis must be IDH-wild-type.")
if not gse154795_table["condition"].eq(gse154795_table["audit_condition"]).all():
    raise ValueError("GSE154795 patient groups differ from the covariate audit table.")


# Select each model's metric and retain the two raw and myeloid-adjusted FDR families.
adjusted_model_table = model_table_total.loc[model_table_total["metric"].eq("raw20_score_myeloid_adjusted")].copy()
raw_model_table = model_table_total.loc[model_table_total["metric"].eq("raw20_score")].copy()
if len(adjusted_model_table) != 2:
    raise ValueError("Figure 5C requires one myeloid-adjusted model for each of the two cohorts.")
if len(raw_model_table) != 2:
    raise ValueError("The Figure 5C FDR audit requires one raw-score model for each of the two cohorts.")
expected = {
    "GSE121810": (
        0.416017028595708, 0.137781145186939, 0.694252912004477, 0.00492100468020317, 0.00984200936040635
    ),
    "GSE154795": (
        0.427369843346401,
        -0.222316569830202,
        1.077056256523,
        0.180119103326531,
        0.180119103326531,
    ),
}
for row in adjusted_model_table.itertuples():
    values = expected[row.cohort]
    observed = (row.effect, row.conf_low, row.conf_high, row.p_value, row.fdr)
    if not all(np.isclose(x, y, atol=1e-12) for x, y in zip(observed, values)):
        raise ValueError(f"The patient-level model results for {row.cohort} have changed.")

model_index_table = adjusted_model_table.set_index("cohort")

# Keep the original model formulas, group fields and FDR family description in the patient table.
patient_table_total = patient_table_total.assign(
    analysis_group = patient_table_total["group"].fillna(patient_table_total["condition"]),
    model_formula  = patient_table_total["cohort"].map(
        {
            "GSE121810": "raw20_score ~ myeloid_abundance + group",
            "GSE154795": (
                "raw20_score ~ myeloid_abundance + MGMT + sequencing_method + "
                "gel_bead_version + acquisition_wave + condition"
            ),
        }
    ),
    fdr_family = "BH across raw and myeloid-adjusted patient-score models within cohort (n=2)",
)
patient_source_table_total = patient_table_total.merge(
    adjusted_model_table[
        [
            "cohort", "effect", "conf_low", "conf_high", "p_value", "fdr", "n_positive", "n_negative"
        ]
    ],
    on  = "cohort",
    how = "left",
    validate = "many_to_one",
)
patient_source_table_total = patient_source_table_total.merge(
    raw_model_table[["cohort", "p_value", "fdr"]].rename(
        columns = {"p_value": "raw_model_p_value", "fdr": "raw_model_fdr"}
    ),
    on  = "cohort",
    how = "left",
    validate = "many_to_one",
)
stem = "Figure5C_patient_scores_and_adjusted_effects"

# Write the panel source table in the original field and row order.
source_path = source_out / f"{stem}_source.csv"
patient_source_table_total.to_csv(source_path, index=False)

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(178, 78)
figure.text(
    0.04,
    0.965,
    "Adjusted Miller-IM score support is significant only in GSE121810",
    fontsize   = 9.0,
    fontweight = "bold",
    va         = "top",
)
figure.text(0.69, 0.885, "Myeloid-proxy-adjusted model", fontsize=6.4, fontweight="bold", ha="center")


# Keep the original group sequence, labels and seeds in both cohorts.
cohort_configurations = [
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
        "title": "GSE154795 strict IDH-wt myeloid pseudobulk",
        "seed": 154795,
    },
]
for index, cohort_config in enumerate(cohort_configurations):
    cohort = cohort_config["cohort"]
    cohort_patient_table = patient_table_total.loc[patient_table_total["cohort"].eq(cohort)].copy()
    group_column = "group" if cohort == "GSE121810" else "condition"
    group_values = [
        cohort_patient_table.loc[cohort_patient_table[group_column].eq(group), "raw20_score"].to_numpy(dtype=float)
        for group in cohort_config["groups"]
    ]
    expected_n = [15, 14] if cohort == "GSE121810" else [8, 13]
    if [len(item) for item in group_values] != expected_n:
        raise ValueError(f"The patient counts by group for {cohort} have changed.")
    distribution_axis = figure.add_axes([0.10, cohort_config["bottom"], 0.47, 0.27])

    # Restart the original seed within each cohort to preserve point jitter.
    rng = np.random.default_rng(cohort_config["seed"])
    for group_index, (group_scores, group_color) in enumerate(zip(group_values, ["#3C5488", "#E64B35"])):
        group_scores = np.asarray(group_scores, dtype=float)
        jitter = rng.uniform(-0.10, 0.10, len(group_scores))
        distribution_axis.scatter(
            np.full(len(group_scores), group_index) + jitter,
            group_scores,
            s         = 18,
            color     = group_color,
            alpha     = 0.78,
            edgecolor = "white",
            linewidth = 0.4,
            zorder    = 3,
        )

        # Use the original 25th, 50th and 75th percentiles and 1.5 IQR whiskers.
        q1, median, q3 = np.quantile(group_scores, [0.25, 0.50, 0.75])
        whisker_low = max(group_scores.min(), q1 - 1.5 * (q3 - q1))
        whisker_high = min(group_scores.max(), q3 + 1.5 * (q3 - q1))
        distribution_axis.vlines(group_index, whisker_low, whisker_high, color="#222222", linewidth=0.8, zorder=4)
        rectangle = mpl.patches.Rectangle(
            (group_index - 0.16, q1),
            0.32,
            q3 - q1,
            facecolor = "white",
            edgecolor = "#222222",
            linewidth = 0.8,
            zorder    = 5,
        )
        distribution_axis.add_patch(rectangle)
        distribution_axis.hlines(
            median, group_index - 0.16, group_index + 0.16, color="#222222", linewidth=1.0, zorder=6
        )

    # Draw each cohort on its original scale with the zero reference line.
    distribution_axis.axhline(0, color="#8A8A8A", linewidth=0.65, linestyle=(0, (3, 2)), zorder=1)
    distribution_axis.set_xticks(range(len(cohort_config["labels"])), cohort_config["labels"], fontsize=5.8)
    distribution_axis.set_xlim(-0.45, len(cohort_config["labels"]) - 0.55)
    distribution_axis.set_ylim(-1.65, 1.45)
    distribution_axis.grid(axis="y", color="#E8E8E8", linewidth=0.45)
    clean_axis(distribution_axis, keep_left=True, keep_bottom=True)
    distribution_axis.tick_params(axis="x", length=0, pad=3)
    distribution_axis.set_title(cohort_config["title"], fontsize=6.8, fontweight="bold", loc="left", pad=3)
    distribution_axis.set_ylabel("Miller-IM score")

    # Use frozen estimates and 95% intervals. Determine fill from full-precision FDR before formatting labels.
    effect_axis = figure.add_axes([0.64, cohort_config["bottom"] + 0.01, 0.31, 0.24])
    model_record = model_index_table.loc[cohort]
    effect_axis.axvline(0, color="#8A8A8A", linewidth=0.7)
    effect_axis.hlines(
        0, float(model_record["conf_low"]), float(model_record["conf_high"]), color="#222222", linewidth=1.1, zorder=2
    )
    significant = float(model_record["fdr"]) < 0.05
    effect_axis.scatter(
        float(model_record["effect"]),
        0,
        marker    = "D",
        s         = 38,
        facecolor = "#E64B35" if significant else "white",
        edgecolor = "#E64B35",
        linewidth = 0.9,
        zorder    = 4,
    )
    effect_axis.text(
        0.98,
        0.76,
        f"FDR {float(model_record['fdr']):.3f}",
        transform  = effect_axis.transAxes,
        ha         = "right",
        va         = "top",
        fontsize   = 6.1,
        color      = "#E64B35" if significant else "#222222",
        fontweight = "bold" if significant else "normal",
    )
    effect_axis.text(
        0.98,
        0.58,
        f"effect {float(model_record['effect']):+.3f}",
        transform = effect_axis.transAxes,
        ha        = "right",
        va        = "top",
        fontsize  = 5.6,
        color     = "#777777",
    )
    effect_axis.set_xlim(-0.40, 1.25)
    effect_axis.set_ylim(-0.55, 0.55)
    effect_axis.set_yticks([])
    effect_axis.set_xticks([0, 0.5, 1.0])
    if index == 1:
        effect_axis.set_xlabel("Adjusted exposure-associated Miller-IM score shift")
    else:
        effect_axis.set_xticklabels([])
    effect_axis.grid(axis="x", color="#E8E8E8", linewidth=0.45)
    clean_axis(effect_axis, keep_left=False, keep_bottom=True)

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step56 recurrent GBM Figure 5")

records.append(
    panel_record(
        "C", stem, "Miller-IM scores with cohort-specific adjusted model effects", source_path, pdf_path, png_path
    )
)


# Keep the 4320 by 5165 PIL canvas, original panel placement, 600 dpi and grayscale preview.

canvas = Image.new("RGB", (4320, 5165), "white")
drawer = ImageDraw.Draw(canvas)
if [row["panel"] for row in records] != list("ABC"):
    raise ValueError("Figure 5 must contain panels A, B and C exactly once.")

by_panel = {row["panel"]: Path(row["png"]) for row in records}
layout = [
    {"panel": "A", "left": 115, "top": 155, "right": 4205, "bottom": 1260},
    {"panel": "B", "left": 115, "top": 1375, "right": 4205, "bottom": 3155},
    {"panel": "C", "left": 115, "top": 3270, "right": 4205, "bottom": 5050},
]
panel_font = image_font(76, bold=True)
for item in layout:
    panel = item["panel"]
    paste_contained(canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"]))
    drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

# Retain the preview's pixel dimensions and resolution for composite layout review.
png_path = figure_out / "Figure5_pd1_preview.png"
canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
preview_gray_path = figure_out / "Figure5_pd1_preview_grayscale.png"
ImageOps.grayscale(canvas).convert("RGB").save(preview_gray_path, format="PNG", dpi=(600, 600), optimize=True)
pdf_path = figure_out / "Figure5_pd1_preview.pdf"
canvas.save(pdf_path, format="PDF", resolution=600.0)
layout_table = pd.DataFrame(layout)
source_path = source_out / "Figure5_pd1_preview_layout.csv"
layout_table.to_csv(source_path, index=False)


# Each manifest row records one panel with the original preview fields.
manifest_table = pd.DataFrame(records)
manifest_table["preview_pdf"] = str(pdf_path)
manifest_table["preview_png"] = str(png_path)
manifest_table["preview_grayscale"] = str(preview_gray_path)
manifest_table["preview_layout"] = str(source_path)
manifest_path = write_root / "Figure5_panel_manifest.csv"
manifest_table.to_csv(manifest_path, index=False)

# Retain the exact Figure 5 legend.

content = (
    "# Figure 5 legend\n\n**Figure 5 | Anti-PD-1 exposure is associated with a coordinated shift of "
    "the Miller-IM program.**\n\n**A,** Sampling and exposure designs. In GSE121810, tumors from 14 "
    "patients were resected after preoperative pembrolizumab exposure, whereas tumors from 15 "
    "adjuvant-only patients were collected before pembrolizumab began. In GSE154795, strict "
    "IDH-wild-type patient-level myeloid pseudobulks compared 13 anti-PD-1-exposed recurrent tumors "
    "with eight recurrent tumors without neoadjuvant anti-PD-1 exposure. The diagram distinguishes "
    "exposure at resection from treatment initiated only after tissue collection. **B,** Rank-based "
    "enrichment of the prespecified Miller-IM program. Nineteen of 20 genes were measurable in each "
    "cohort; AC253572.2 was not measured. In GSE121810, the Miller-IM program was enriched toward "
    "the preoperative-exposure end of the transcriptome-wide ranking (NES = 1.900; nominal fgsea *P*"
    " = 6.35 × 10^-5), with a complementary competitive test agreeing in direction (camera, Up; "
    "nominal *P* = 2.81 × 10^-5). In the strict IDH-wild-type GSE154795 myeloid pseudobulk analysis,"
    " the program was likewise enriched toward anti-PD-1 exposure (NES = 1.751; nominal fgsea *P* = "
    "2.91 × 10^-4; camera, Up; nominal *P* = 0.0154). Curves show the weighted running enrichment "
    "score; ticks mark Miller-IM genes in the ranked transcriptome. **C,** Left, unadjusted "
    "patient-level Miller-IM scores, displayed as individual patients with box summaries. Right, "
    "model-based exposure coefficients and 95% CIs after adjustment for the myeloid-expression "
    "proxy; the GSE154795 model additionally retained MGMT and technical covariates. GSE121810 "
    "showed an adjusted effect of +0.416 (95% CI, +0.138 to +0.694; FDR = 0.00984), whereas "
    "GSE154795 showed a positive but nonsignificant adjusted effect of +0.427 (95% CI, -0.222 to "
    "+1.077; FDR = 0.180).\n\nThese cohorts test an exposure-associated coordinated program shift, not"
    " clinical response prediction. GSE121810 is a whole-tumor bulk RNA assay and cannot localize "
    "the signal to myeloid cells; GSE154795 provides patient-level myeloid pseudobulk support but "
    "its patient-level score did not reach FDR < 0.05. The main conclusion therefore concerns "
    "coordinated rank-based program enrichment under anti-PD-1 exposure. It does not establish "
    "drug-induced causality, treatment benefit, or validity across all IDH groups; the GSE154795 "
    "all-IDH sensitivity analysis was not significant. Entry sensitivity and historical "
    "response-labeled analyses remain supplementary and are not presented as Miller-IM response "
    "validation.\n"
)
legend_path = write_root / "Figure5_legend.md"
legend_path.write_text(content, encoding="utf-8")


content = (
    "# Figure 5 anti-PD-1 exposure rebuild\n\n"
    "- A. Timing of anti-PD-1 exposure at tissue collection in both cohorts.\n"
    "- B. Miller-IM program rank enrichment and agreement with camera in both cohorts.\n"
    "- C. Raw patient scores alongside myeloid-proxy-adjusted model effects. The GSE154795 patient-score result "
    "remains nonsignificant (FDR=0.180).\n"
    "- The layout omits volcano plots, the four-cohort leading-edge UpSet plot and response-prediction cards.\n"
    "- Main claim. Anti-PD-1 exposure-associated coordinated shift. "
    "Claims of induction or response prediction fall outside this scope.\n"
    "- Manifest `"
    f"{manifest_path}"
    "`\n- Legend `"
    f"{legend_path}"
    "`\n"
)
(write_root / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


# Check the actual output files and preview canvas dimensions.
output_paths = [
    manifest_path,
    legend_path,
    write_root / "FINAL_RESULT.md",
    pdf_path,
    png_path,
    preview_gray_path,
    source_path,
]
for panel_result in records:
    for field in ["source", "pdf", "png"]:
        output_paths.append(Path(panel_result[field]))
for output_path in output_paths:

    # Check that each current output exists and has content.
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {output_path.name}")
with Image.open(png_path) as preview_image:
    if preview_image.size != canvas.size:
        raise ValueError("Preview dimensions differ from the original canvas")
    preview_image.verify()

# Read the saved manifest and compare its dimensions to the current panel records.
manifest_check_table = pd.read_csv(manifest_path)
if manifest_check_table.shape != manifest_table.shape:
    raise ValueError("Panel manifest dimensions changed during writing")
print(f"Figure5 panels={len(records)} preview=1")
