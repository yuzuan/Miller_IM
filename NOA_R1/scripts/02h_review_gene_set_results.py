#!/usr/bin/env python3
"""Independently check the archived gene-set and PD-1 source-table calculations.

Inputs: provenance/step1_input_manifest.csv, the Miller source workbook under sources/, the 01a,
    02c, 02d and 02f result tables, archived patient_paired_deltas.csv, and the legacy Figure 5 CAMERA references.
Outputs: step2_independent_review/validation_checks.csv and independent_review.md.
Run: python NOA_R1/scripts/02h_review_gene_set_results.py after 02c, 02d and 02f.
"""

import csv
import importlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats


# Git and the archive extraction do not retain the original file time sequence.
revision = Path(__file__).resolve().parents[1]
check_table = []


def check(area, label, passed, detail=""):
    """Append one validation row. A missing result means original evidence is required. Returns nothing."""
    row = {
        "area": area,
        "check": label,
        "passed": None if passed is None else bool(passed),
        "detail": str(detail),
    }
    check_table.append(row)


check(
    "timing",
    "plans_before_new_results",
    None,
    "Requires the original timestamp-preserving plans and result files. "
    "See step2_analysis_plan.md for the lock record.",
)

# Verify the archived systemic set against the ranked Miller source.
extended_dir = revision / "step2_extended_specificity"
definition_table = pd.read_csv(extended_dir / "locked_gene_sets_used.csv")
workbook_path = revision / "sources" / "Miller_2025_Supplementary_Table_2.xlsx"
workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
try:
    workbook_rows = list(workbook["Myeloid Top 100 genes per prog."].values)
finally:
    workbook.close()
systemic_column = workbook_rows[0].index("Systemic Inflammatory")
original_systemic = [row[systemic_column] for row in workbook_rows[1:21]]
actual_systemic = definition_table.loc[definition_table["set_id"] == "Miller_systemic_top20", "gene"].tolist()
check("source", "systemic_complete_original_top20", actual_systemic == original_systemic)
expected_sets = {
    "target_top20", "target_top50", "target_top100", "Abdelfattah_iMic_top50", "Sankowski_C1", "Miller_systemic_top20",
}
check("source", "historical_source_sets_preserved", set(definition_table["set_id"]) == expected_sets)
status_table = pd.read_csv(extended_dir / "subitem_status.csv")
check(
    "source",
    "Hallmark_stopped_without_download",
    len(status_table) == 1 and status_table["status"].iat[0] == "stopped" and
        "cache absent" in status_table["reason"].iat[0],
)

# The manifest supplies each cohort's formal pair count and expected filtered gene universe.
manifest_table = pd.read_csv(revision / "provenance" / "step1_input_manifest.csv")
coverage_table = pd.read_csv(extended_dir / "gene_coverage.csv")
camera_table = pd.read_csv(extended_dir / "camera_results.csv")
effect_table = pd.read_csv(extended_dir / "patient_effects.csv")
score_table = pd.read_csv(extended_dir / "patient_scores.csv")
random_table = pd.read_csv(extended_dir / "random_summary.csv")
pair_delta_table = pd.read_csv(extended_dir / "patient_paired_deltas.csv")
primary_camera_table = pd.read_csv(revision / "step1_primary_recurrence" / "camera_results.csv")
expected_deltas = int(manifest_table["expected_pairs"].sum()) * len(expected_sets)
expected_effects = len(manifest_table) * len(expected_sets)
check("structure", "patient_delta_rows_complete", len(pair_delta_table) == expected_deltas)
check(
    "structure",
    "complete_camera_rows",
    len(camera_table) == 2 * expected_effects and set(camera_table["mode"]) == {"estimated", "fixed_0.01"},
)
check(
    "structure",
    "complete_patient_effects",
    len(effect_table) == expected_effects and len(score_table) == 2 * expected_deltas,
)


def bh(p_values):
    """Return Benjamini-Hochberg adjusted p values in their input order."""
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    scaled = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum.accumulate(scaled[::-1])[::-1]
    result = np.minimum(adjusted, 1)
    return result


for cohort in manifest_table.itertuples(index=False):
    dataset = cohort.dataset
    n_pairs = cohort.expected_pairs
    background_size = cohort.expected_tested_genes
    cohort_camera_table = camera_table[camera_table["dataset"] == dataset]
    cohort_effect_table = effect_table[effect_table["dataset"] == dataset]
    n_camera_modes = cohort_camera_table.groupby("mode").size()
    check(dataset, "six_sets_both_camera_modes", len(cohort_camera_table) == 2 * len(expected_sets) and
        all(n_camera_modes == len(expected_sets)))
    fixed_camera_table = cohort_camera_table[cohort_camera_table["mode"] == "fixed_0.01"]
    check(dataset, "same_pair_counts_and_background", set(cohort_camera_table["n_pairs"]) == {n_pairs} and
        set(fixed_camera_table["camera_reference_df"]) == {background_size - 2})

    # Compute BH within each mode and compare the target result with the corresponding single-set 01a p value.
    for mode in ("estimated", "fixed_0.01"):
        mode_table = cohort_camera_table[cohort_camera_table["mode"] == mode]
        adjusted = bh(mode_table["p_value"])
        check(dataset, f"BH_{mode}", np.allclose(mode_table["bh_within_cohort_mode"], adjusted, rtol=0, atol=1e-13))
        primary_rows = (primary_camera_table["dataset"] == dataset) & (primary_camera_table["correlation_mode"] == mode)
        original_p = primary_camera_table.loc[primary_rows, "p_value"].item()
        current_p = mode_table.loc[mode_table["set_id"] == "target_top20", "p_value"].item()
        check(dataset, f"Step1_target_camera_{mode}", abs(current_p - original_p) < 1e-12)

    # Independently pair the patient scores, then recalculate each t test and 95% interval.
    cohort_score_table = score_table[score_table["dataset"] == dataset]
    for set_id, group_table in cohort_score_table.groupby("set_id"):
        primary_table = group_table[group_table["condition"] == "Primary"].set_index("patient_id")
        recurrent_table = group_table[group_table["condition"] == "Recurrent"].set_index("patient_id")
        check(dataset, f"{set_id}_pairing", len(primary_table) == len(recurrent_table) == n_pairs and
            set(primary_table.index) == set(recurrent_table.index))
        delta = (recurrent_table.loc[primary_table.index, "score"] - primary_table["score"]).to_numpy()
        saved_rows = (pair_delta_table["dataset"] == dataset) & (pair_delta_table["set_id"] == set_id)
        archived_table = pair_delta_table[saved_rows].set_index("patient_id")
        archived_delta = archived_table.loc[primary_table.index, "delta_recurrent_minus_primary"]
        check(dataset, f"{set_id}_saved_patient_deltas", np.allclose(archived_delta, delta, rtol=0, atol=1e-12))
        set_effect_table = cohort_effect_table[cohort_effect_table["set_id"] == set_id]
        paired_test = stats.ttest_1samp(delta, 0)
        standard_error = delta.std(ddof=1) / math.sqrt(n_pairs)
        ci_low, ci_high = stats.t.interval(.95, n_pairs - 1, loc=delta.mean(), scale=standard_error)
        errors = (
            abs(delta.mean() - set_effect_table["mean_delta"].item()),
            abs(ci_low - set_effect_table["ci_low"].item()),
            abs(ci_high - set_effect_table["ci_high"].item()),
            abs(paired_test.pvalue - set_effect_table["paired_t_p"].item()),
        )
        check(dataset, f"{set_id}_patient_effect", max(errors) < 1e-11)

    # Check the eligible gene universe and the declared target exclusion before testing saved draws.
    target_rows = (coverage_table["dataset"] == dataset) & (coverage_table["set_id"] == "target_top20")
    target_table = coverage_table[target_rows]
    check(dataset, "target_nineteen_coverage", len(target_table) == 20 and target_table["tested"].sum() == 19)
    background_table = pd.read_csv(extended_dir / f"{dataset}_background.csv").set_index("gene")
    eligible = set(background_table.index[background_table["random_eligible"]])
    target = set(target_table["gene"])
    check(dataset, "random_background_size_and_exclusion", len(background_table) == background_size and
        len(eligible) == background_size - 19 and not eligible.intersection(target))
    target_bins = Counter(background_table.loc[list(target.intersection(background_table.index)), "expression_decile"])
    cohort_random_table = random_table[random_table["dataset"] == dataset]

    # Each scheme uses 10,000 draws of 19 distinct eligible genes. Match expression bins in the matched scheme.
    for summary in cohort_random_table.itertuples(index=False):
        draw_table = pd.read_csv(extended_dir / f"{dataset}_{summary.scheme}_random_sets.csv")
        check(dataset, f"{summary.scheme}_B", len(draw_table) == 10000)
        valid = True
        matched = True
        for membership in draw_table["members"]:
            genes = membership.split(";")
            valid = valid and len(genes) == 19 and len(set(genes)) == 19 and set(genes) <= eligible
            if summary.scheme == "size_and_expression":
                matched = matched and Counter(background_table.loc[genes, "expression_decile"]) == target_bins
        check(dataset, f"{summary.scheme}_members_valid", valid)
        if summary.scheme == "size_and_expression":
            check(dataset, f"{summary.scheme}_decile_match", matched)

        # Recalculate both tails from the stored ES and the original plus-one correction.
        upper_count = int(np.sum(draw_table["ES"] >= summary.observed_es))
        absolute_count = int(np.sum(abs(draw_table["ES"]) >= abs(summary.observed_es)))
        check(
            dataset,
            f"{summary.scheme}_empirical_formula",
            upper_count == summary.n_es_ge and absolute_count == summary.n_abs_ge and
                abs((upper_count + 1) / 10001 - summary.empirical_p) < 1e-12 and
                abs((absolute_count + 1) / 10001 - summary.abs_empirical_p) < 1e-12,
        )

check(
    "random",
    "four_test_BH",
    np.allclose(random_table["bh_four_tests"], bh(random_table["empirical_p"]), atol=1e-12) and
        np.allclose(random_table["abs_bh_four_tests"], bh(random_table["abs_empirical_p"]), atol=1e-12),
)


# Use the saved five-set ranking to independently recompute weighted ES for the three original draw indices.
def independent_es(order, rank, genes):
    """Calculate standard weighted GSEA ES for one stored set and a fixed signed ranking."""
    selected = set(genes)
    weight_total = sum(abs(rank[gene]) for gene in selected)
    hit = 0.0
    miss = 0.0
    maximum = -float("inf")
    minimum = float("inf")
    for gene in order:
        if gene in selected:
            hit += abs(rank[gene]) / weight_total
        else:
            miss += 1 / (len(order) - len(selected))
        running = hit - miss
        maximum = max(maximum, running)
        minimum = min(minimum, running)
    result = maximum if maximum >= -minimum else minimum
    return result


for dataset in manifest_table["dataset"]:
    rank_table = pd.read_csv(revision / "step2_gene_set_specificity" / f"{dataset}_ranking.csv")
    order = rank_table["gene"].tolist()
    rank = dict(zip(order, rank_table["statistic"], strict=True))
    for scheme in ("size_only", "size_and_expression"):
        draw_table = pd.read_csv(extended_dir / f"{dataset}_{scheme}_random_sets.csv").set_index("draw")
        for draw in (1, 5000, 10000):
            stored_genes = draw_table.at[draw, "members"].split(";")
            es = independent_es(order, rank, stored_genes)
            check(dataset, f"{scheme}_independent_ES_draw_{draw}", abs(es - draw_table.at[draw, "ES"]) < 1e-9)

# Check the original fixed PD-1 P values and formal sample/design alignment.
pd1_dir = revision / "step2_pd1_camera"
pd1_path = pd1_dir / "camera_results.csv"
if pd1_path.is_file():
    pd1_table = pd.read_csv(pd1_path)
    revision_tools = importlib.import_module("00a_functions")
    project = revision_tools.legacy_root(revision)
    reference_dir = project / "write" / "41_mg_inflammatory_sci_rebuild" / "Figure5"
    check("PD1", "four_modes", len(pd1_table) == 4 and set(pd1_table["mode"]) == {"estimated", "fixed_0.01"})
    references = (
        ("GSE121810", "GSE121810_raw20_camera.csv", 2.81295477233785e-5),
        ("GSE154795", "GSE154795_raw20_camera_auto_strict_idhwt.csv", .015367499497791),
    )
    for dataset, name, expected_p in references:
        reference_table = pd.read_csv(reference_dir / name)
        source_p = reference_table.loc[reference_table["pathway"] == "Miller_raw20", "PValue"].item()
        cohort_pd1_table = pd1_table[pd1_table["dataset"] == dataset]
        fixed_table = cohort_pd1_table[cohort_pd1_table["mode"] == "fixed_0.01"]
        fixed_p = fixed_table["p_value"].item()
        check("PD1", f"{dataset}_fixed_reproduces_original",
            abs(fixed_p - source_p) <= 1e-10 and abs(fixed_p - expected_p) <= 1e-10)
        check("PD1", f"{dataset}_same_19_genes", set(cohort_pd1_table["n_genes"]) == {19})

        design_table = pd.read_csv(pd1_dir / f"{dataset}_design.csv")
        sample_table = pd.read_csv(pd1_dir / "sample_audit.csv")
        cohort_sample_table = sample_table[sample_table["dataset"] == dataset]
        check(
            "PD1",
            f"{dataset}_sample_design_alignment",
            len(design_table) == len(cohort_sample_table) == fixed_table["n_samples"].item() and
                set(design_table["sample"]) == set(cohort_sample_table["sample"]),
        )
else:
    check("PD1", "camera_not_yet_available", False, "PD-1 run incomplete")

# Save the original validation columns. The chronology row is blank pending original timestamp-preserving files.
review_dir = revision / "step2_independent_review"
review_dir.mkdir(parents=True, exist_ok=True)
validation_path = review_dir / "validation_checks.csv"
fields = ["area", "check", "passed", "detail"]
with validation_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(check_table)

with validation_path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    written_fields = reader.fieldnames
    written_count = sum(1 for row in reader)
if written_fields != fields or written_count != len(check_table):
    raise ValueError("Written validation table differs from the expected header or row count")

passed = sum(row["passed"] is True for row in check_table)
failed = sum(row["passed"] is False for row in check_table)
not_assessed = len(check_table) - passed - failed
summary = (
    "# Step 2 independent review\n\n"
    "The checks use the published Miller workbook, archived source tables, and legacy CAMERA reference results.\n\n"
    f"{passed} passed, {failed} failed, and {not_assessed} not assessed. See validation_checks.csv.\n\n"
    "Chronology requires the original timestamp-preserving plans and outputs. "
    "The documented protocol lock is in step2_analysis_plan.md.\n"
)
summary_path = review_dir / "independent_review.md"
summary_path.write_text(summary, encoding="utf-8")
if summary_path.read_text(encoding="utf-8") != summary:
    raise ValueError("Written summary differs from the review results")
print(f"Independent review: {passed} passed, {failed} failed, {not_assessed} not assessed.")
