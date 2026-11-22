#!/usr/bin/env python3
"""Build editable Supplementary Tables S1–S8 from frozen October NOA R1 CSVs.

Inputs by table (paths relative to NOA_R1):
    S1: step5_closeout_v2/{overlap_source_evidence,overlap_28_pair_reaudit}.csv.
    S2: step1_primary_recurrence/{camera_results,patient_score_effects}.csv and
        step5_filter_consistency/{camera_results,patient_effects}.csv.
    S3: step2_hallmark_completion/{camera_seven_set_family,hallmark_patient_effects}.csv and
        step2_extended_specificity/{patient_effects,random_summary}.csv.
    S4: step3_state_analysis/{state_abundance_statistics,state_program_camera,state_score_effects}.csv.
    S5: step4_steroids_geomx/camera_results.csv and clinical_case_normalized/
        {steroid_summary_all_statuses,steroid_missing_source_audit,GSE274546_center_summary}.csv.
    S6: step4_steroids_geomx/{GeoMx_subgroup_interaction,GeoMx_Welch_statistic_df}.csv.
    S7: step5_protein_common/effects.csv.
    S8: step5_closeout_v2/{model_results_with_original_raw,cook_threshold_patients}.csv,
        step5_pd1_stability_final/{leave_one_out,vif,covariate_distributions}.csv,
        step2_pd1_camera/camera_results.csv and
        step5_rank_unification_v2/{sqrtF_fgsea,historical_F_fgsea,negative_F_genes}.csv.
Output: supplementary_R1/Supplementary_Tables_S1-S8.docx.
Run: python NOA_R1/scripts/09c_supplementary_tables.py
"""

from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as etree
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd


word_namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
word_tag = f"{{{word_namespace}}}"
etree.register_namespace("w", word_namespace)


def add_element(parent, name, attributes=None):
    """Append a WordprocessingML element and return it with qualified attributes."""
    qualified_attributes = {}
    for key, value in (attributes or {}).items():
        qualified_attributes[word_tag + key] = str(value)
    element = etree.SubElement(parent, word_tag + name, qualified_attributes)
    return element


def append_paragraph(parent, text, bold=False, size=16) -> None:
    """Append a paragraph using the original half-point sizes and spacing. No return value."""
    paragraph = add_element(parent, "p")
    properties = add_element(paragraph, "pPr")
    add_element(properties, "spacing", {"before": "90", "after": "50"})

    run = add_element(paragraph, "r")
    run_properties = add_element(run, "rPr")
    add_element(run_properties, "sz", {"val": size})
    if bold:
        add_element(run_properties, "b")
    text_element = add_element(run, "t")
    text_element.text = str(text)


def append_section(parent, number, title, note) -> None:
    """Add a numbered supplementary-table title and its frozen caption. No return value."""
    append_paragraph(parent, f"Supplementary Table S{number}. {title}", True, 21)
    append_paragraph(parent, note, size=15)


def append_table(parent, display_labels, expected_shapes, column_names, rows) -> None:
    """Add a table with a repeating header and check the fields in each row. No return value."""
    table_element = add_element(parent, "tbl")
    properties = add_element(table_element, "tblPr")
    add_element(properties, "tblW", {"w": 0, "type": "auto"})
    borders = add_element(properties, "tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        add_element(borders, edge, {"val": "single", "sz": "4", "color": "C5CBD0"})

    # Keep the original automatic table and cell widths for the landscape page.
    for row_index, values in enumerate([column_names, *rows]):
        if len(values) != len(column_names):
            raise ValueError(f"Table {len(expected_shapes) + 1}, row {row_index + 1} has an invalid field count")
        row_element = add_element(table_element, "tr")
        if row_index == 0:
            row_properties = add_element(row_element, "trPr")
            add_element(row_properties, "tblHeader", {"val": 1})

        for value in values:
            cell = add_element(row_element, "tc")
            cell_properties = add_element(cell, "tcPr")
            add_element(cell_properties, "tcW", {"w": 0, "type": "auto"})
            paragraph = add_element(cell, "p")
            paragraph_properties = add_element(paragraph, "pPr")
            add_element(paragraph_properties, "spacing", {"after": 0})

            # Preserve blank cells and the original translation of display-only labels.
            run = add_element(paragraph, "r")
            run_properties = add_element(run, "rPr")
            add_element(run_properties, "sz", {"val": 15})
            if row_index == 0:
                add_element(run_properties, "b")
            text = add_element(run, "t")
            text.text = "" if pd.isna(value) else display_labels.get(str(value), str(value))

    expected_shapes.append((len(rows) + 1, len(column_names)))
    append_paragraph(parent, "")


def format_number(value, digits=3) -> str:
    """Format an output cell with the original dash, precision and scientific notation."""
    if pd.isna(value):
        label = "—"
    elif isinstance(value, str):
        label = value
    elif value != 0 and abs(value) < 0.0001:
        label = f"{value:.2e}"
    else:
        label = f"{value:.{digits}f}"
    return label


def format_ci(row) -> str:
    """Return the confidence limits at the original three-decimal precision."""
    label = f"{format_number(row.ci_low)} to {format_number(row.ci_high)}"
    return label


def main() -> None:
    """Build the eight supplements in order and verify the generated document."""
    revision_dir = Path(__file__).resolve().parents[1]
    document = etree.Element(word_tag + "document")
    body = add_element(document, "body")
    write_paragraph = partial(append_paragraph, body)
    write_section = partial(append_section, body)

    # Display labels belong to the Word output. The CSV field names and values stay intact.
    display_labels = {
        "GeoMx_Zenodo16839828": "Artzi GeoMx",
        "target_top20": "MG-Inflam20 (top 20)",
        "target_top50": "Microglial inflammatory (top 50)",
        "target_top100": "Microglial inflammatory (top 100)",
        "Abdelfattah_iMic_top50": "i-Mic (top 50)",
        "Sankowski_C1": "Sankowski C1",
        "Miller_systemic_top20": "Systemic inflammatory (top 20)",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Hallmark inflammatory response",
        "GOBP_RESPONSE_TO_CORTICOSTEROID": "Corticosteroid response",
        "target_minus_SGK1_PDK4": "MG-Inflam20 minus SGK1/PDK4",
        "fixed_0.01": "Fixed 0.01", "estimated": "Estimated",
        "size_only": "Size matched", "size_and_expression": "Size and expression matched",
        "MAC_MDSC-like": "MAC/MDSC-like",
        "insufficient_paired_patients": "Not tested: insufficient pairs",
        "program_coverage_below_16_of_20": "Not tested: fewer than 16 genes",
        "analyzed": "Tested", "fixed_proteins_all105": "Fixed proteins, all patients",
        "original11_complete_patients": "All 11 proteins, complete cases",
        "unique_peptides": "Unique peptides", "all_peptides": "All peptides",
        "unadjusted": "Unadjusted", "original_raw": "Original covariates",
        "proxy_only": "Proxy only", "full": "All covariates plus proxy",
        "Control_paired_delta": "Paired change, controls",
        "Nivolumab_paired_delta": "Paired change, nivolumab",
        "Nivolumab_minus_Control_interaction": "Difference in change, nivolumab minus controls",
        "Adj": "Adjuvant-only pembrolizumab", "Neo": "Preoperative pembrolizumab",
        "GBM.rec": "PD-1 unexposed", "GBM.PD1": "PD-1 exposed",
        "myeloid_abundance": "Myeloid-expression proxy", "n_myeloid_cells": "Myeloid cell count",
        "sequencing_method": "Sequencing method", "gel_bead_version": "Gel-bead version",
        "acquisition_wave": "Acquisition wave", "MGMT": "MGMT promoter methylation",
        "Meth": "Methylated", "Un": "Unmethylated",
        "No->No": "No → No", "No->Yes": "No → Yes", "Yes->No": "Yes → No", "Yes->Yes": "Yes → Yes",
        "Missing": "Incomplete steroid records", "Recorded": "Complete steroid records", "All": "All patients",
    }
    expected_shapes = []
    write_table = partial(append_table, body, display_labels, expected_shapes)

    # S1. Each source supplies one provenance row. The second table lists all archived comparisons.
    closeout_dir = revision_dir / "step5_closeout_v2"
    site_table = pd.read_csv(closeout_dir / "overlap_source_evidence.csv")
    overlap_table = pd.read_csv(closeout_dir / "overlap_28_pair_reaudit.csv")
    if len(site_table) != 8 or len(overlap_table) != 28:
        raise ValueError("S1 requires the eight sources and 28 comparisons described in its frozen caption")
    write_section(
        1, "Source provenance and possible patient overlap",
        "Eight public sources, 28 pairwise comparisons. Categories describe source evidence; no cross-institution "
        "identity linkage was available. The GSE121810–GSE274546 shared-site inference uses a later publication "
        "from the same trial; this does not establish the sites of the original analyzed samples.",
    )
    rows = []
    for row in site_table.itertuples():
        provenance = str(row.sites).replace("|", ", ")
        if row.dataset == "GSE121810":
            provenance = "multicenter trial (NCT02852655); later-stage site list includes UCLA, UCSF and MD Anderson"
        rows.append((row.dataset, provenance, row.paper))
    write_table(["Dataset", "Collection provenance", "Source publication"], rows)

    # Recognize the three original archive codes as well as the English published codes.
    archive_categories = {
        "\u5df2\u77e5\u91cd\u53e0": "known overlap",
        "\u6765\u6e90\u76f8\u5173\u3001\u53ef\u80fd\u91cd\u53e0": "possible overlap due to shared sources",
        "\u6709\u8bc1\u636e\u652f\u6301\u65e0\u91cd\u53e0": "evidence supports no overlap",
    }
    category_labels = {
        "known overlap": "Known overlap",
        "possible overlap due to shared sources": "Related sources; overlap possible",
        "evidence supports no overlap": "Different institutions; no evidence of overlap",
    }
    rows = []
    for row in overlap_table.itertuples():
        category = archive_categories.get(row.classification, row.classification)
        if row.pair_number == 14:
            evidence = "Shared patient IDs GBM030 and GBM049 in the Jackson datasets."
        elif row.pair_number == 12:
            evidence = "CARE MD Anderson contribution; MD Anderson appears in a later site list for NCT02852655."
        elif row.pair_number in (17, 18, 21, 22):
            evidence = "Jackson collection at Johns Hopkins; comparison collection at another institution."
        elif category == "evidence supports no overlap":
            evidence = "Collection institutions differ in the source publications."
        elif row.pair_number == 28:
            evidence = "UCLA was a trial site; Lee patients were treated at UCLA, with no ID mapping."
        elif row.pair_number == 11:
            evidence = "Both publications list SNU and Pitié-Salpêtrière as collection sites."
        else:
            evidence = "A contributing site is shared; no patient-ID mapping available."
        rows.append((row.pair_number, row.source_a, row.source_b, category_labels[category], evidence))
    write_table(["No.", "Source A", "Source B", "Classification", "Evidence summary"], rows)
    write_paragraph(
        "Abbreviations: CARE, GBM-CARE; SNU, Seoul National University; UCLA, University of California Los Angeles; "
        "UCSF, University of California San Francisco.",
        size=15,
    )

    # S2. Export the paired camera tests and score effects at the primary and sensitivity filters.
    write_section(
        2, "Paired recurrence testing and expression-filter sensitivity",
        "Competitive camera used a paired voom design. The primary analysis estimated inter-gene correlation; "
        "a fixed correlation of 0.01 was evaluated in a sensitivity analysis. Score changes are patient-level "
        "paired means with t-based 95% confidence intervals (CIs).",
    )
    recurrence_dir = revision_dir / "step1_primary_recurrence"
    camera_table = pd.read_csv(recurrence_dir / "camera_results.csv")
    rows = []
    for row in camera_table.itertuples():
        background = f"{row.n_genes}/{row.tested_universe}"
        values = (
            row.dataset, row.correlation_mode, row.n_pairs, background, row.direction,
            format_number(row.p_value, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Mode", "Pairs", "Genes / background", "Direction", "camera P"], rows)

    effect_table = pd.read_csv(recurrence_dir / "patient_score_effects.csv")
    rows = []
    for row in effect_table.itertuples():
        pairs = f"{row.positive_pairs}/{row.n_pairs}"
        paired_p = format_number(row.t_p_value, 5)
        rows.append((row.dataset, pairs, format_number(row.mean_delta), format_ci(row), paired_p))
    write_table(["Cohort", "Pairs higher / total", "Mean change", "95% CI", "Paired t P"], rows)

    filter_dir = revision_dir / "step5_filter_consistency"
    filter_camera_table = pd.read_csv(filter_dir / "camera_results.csv")
    rows = []
    for _, row in filter_camera_table.iterrows():
        mode = "Condition / " + str(row["mode"])
        background = f"{row.n_genes}/{row.background}"
        rows.append((row.dataset, mode, row.n_pairs, background, format_number(row.p_value, 5)))
    write_table(["Cohort", "Filter / mode", "Pairs", "Genes / background", "camera P"], rows)

    filter_effect_table = pd.read_csv(filter_dir / "patient_effects.csv")
    rows = []
    for row in filter_effect_table.itertuples():
        values = (
            row.dataset, row.n_pairs, row.n_genes, format_number(row.mean_delta), format_ci(row),
            format_number(row.paired_t_p, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Pairs", "Genes", "Mean change", "95% CI", "t P"], rows)

    # S3. Keep camera, score effects and random-set tests in the original table order.
    write_section(
        3, "Gene-set specificity and random-set comparisons",
        "Seven sets were tested per cohort. Benjamini–Hochberg q values cover the seven-set family separately "
        "within each cohort and camera mode. Random sets exclude the 20 target genes and do not preserve "
        "gene–gene correlation.",
    )
    hallmark_dir = revision_dir / "step2_hallmark_completion"
    seven_set_table = pd.read_csv(hallmark_dir / "camera_seven_set_family.csv")
    rows = []
    for _, row in seven_set_table.iterrows():
        values = (
            row.dataset, row.set_id, row["mode"], row.n_pairs, row.n_genes, row.direction,
            format_number(row.p_value, 5), format_number(row.bh_seven, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Gene set", "Mode", "Pairs", "Genes", "Direction", "camera P", "BH q"], rows)

    specificity_dir = revision_dir / "step2_extended_specificity"
    effect_table = pd.read_csv(specificity_dir / "patient_effects.csv")
    rows = []
    for row in effect_table.itertuples():
        pairs = f"{row.positive_pairs}/{row.n_pairs}"
        values = (
            row.dataset, row.set_id, pairs, format_number(row.mean_delta), format_ci(row),
            format_number(row.paired_t_p, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Set", "Pairs higher / total", "Mean change", "95% CI", "Paired t P"], rows)

    hallmark_table = pd.read_csv(hallmark_dir / "hallmark_patient_effects.csv")
    rows = [(row.dataset, format_number(row.mean_delta), format_ci(row)) for row in hallmark_table.itertuples()]
    write_table(["Cohort", "Hallmark mean change", "95% CI"], rows)

    random_table = pd.read_csv(specificity_dir / "random_summary.csv")
    rows = []
    for row in random_table.itertuples():
        values = (
            row.dataset, row.scheme, row.B, row.n_es_ge, format_number(row.empirical_p, 5),
            format_number(row.bh_four_tests, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Matching", "Sets", "Sets ≥ observed ES", "Empirical P", "BH q"], rows)
    write_paragraph(
        "Abbreviations: BH, Benjamini–Hochberg; CI, confidence interval; ES, enrichment score; "
        "i-Mic, inflammatory microglia.",
        size=15,
    )

    # S4. Abundance records are one per cohort and state. The following tables use the stated cell thresholds.
    write_section(
        4, "Paired myeloid-state abundance, program coverage, and camera",
        "Abundance is the fraction of captured myeloid cells; changes are recurrent minus primary percentage "
        "points. Within-state models retain patient blocking but filter by condition, requiring ≥20 cells at each "
        "endpoint and ≥6 pairs; program testing additionally requires ≥16/20 filtered genes. Program P values "
        "use Benjamini–Yekutieli correction across 11 tests per cohort, assigning unavailable tests P=1 for "
        "correction but reporting them as not tested. Fine-state abundance uses nine-test correction; complementary "
        "broad fractions represent one hypothesis.",
    )
    state_dir = revision_dir / "step3_state_analysis"
    abundance_table = pd.read_csv(state_dir / "state_abundance_statistics.csv")
    rows = []
    for row in abundance_table.itertuples():
        interval = f"{format_number(row.ci_low_pp)} to {format_number(row.ci_high_pp)}"
        values = (
            row.dataset, row.state, row.n_pairs, format_number(row.mean_delta_pp), interval,
            format_number(row.signflip_p, 5), format_number(row.primary_adjusted_p, 5),
        )
        rows.append(values)
    write_table(["Cohort", "State", "Pairs", "Change (pp)", "95% CI (pp)", "Sign-flip P", "Adjusted P"], rows)

    camera_table = pd.read_csv(state_dir / "state_program_camera.csv")
    column_names = ["Cohort", "State", "Eligible pairs", "Filtered genes", "Program status", "camera P", "BY q"]
    primary_table = camera_table.loc[(camera_table.threshold == 20) & (camera_table.correlation_mode == "estimated")]
    rows = []
    for row in primary_table.itertuples():
        values = (
            row.dataset, row.state, row.n_pairs, format_number(row.n_genes, 0), row.status,
            format_number(row.p_value, 5), format_number(row.by_eleven_primary, 5),
        )
        rows.append(values)
    write_table(column_names, rows)
    write_paragraph(
        "Abbreviations: BY, Benjamini–Yekutieli; MCG, microglia; MAC, macrophage; MDSC, myeloid-derived suppressor "
        "cell; pp, percentage points. State labels follow Jackson et al. A dash for P or q indicates that the "
        "prespecified coverage criteria prevented testing.",
        size=15,
    )

    # The 10-cell sensitivity lists the full 11-group family for each of the two cohorts in both modes.
    for mode, label in (
        ("estimated", "Estimated inter-gene correlation"), ("fixed_0.01", "Fixed inter-gene correlation of 0.01")
    ):
        write_paragraph("Ten-cell threshold sensitivity: " + label, True)
        sensitivity_table = camera_table.loc[(camera_table.threshold == 10) & (camera_table.correlation_mode == mode)]
        if len(sensitivity_table) != 22:
            raise ValueError(f"S4 is missing one or more of the 22 prespecified 10-cell records for {mode}")
        rows = []
        for row in sensitivity_table.itertuples():
            values = (
                row.dataset, row.state, row.n_pairs, format_number(row.n_genes, 0), row.status,
                format_number(row.p_value, 5), format_number(row.by_eleven_primary, 5),
            )
            rows.append(values)
        write_table(column_names, rows)
    write_paragraph(
        "The 10-cell sensitivity retained the same minimum of six paired patients and 16 filtered target genes. "
        "Benjamini–Yekutieli adjustment used 11 prespecified groupings within each cohort and correlation setting; "
        "unavailable tests were assigned P = 1 only for correction and are displayed as not tested. A dash indicates "
        "an unavailable value.",
        size=15,
    )

    # Keep the descriptive within-state score records separate from the competitive tests.
    write_paragraph("Descriptive within-state RNA-score changes", True)
    score_table = pd.read_csv(state_dir / "state_score_effects.csv")
    if len(score_table) != 19:
        raise ValueError("S4 requires the original 19 descriptive within-state score records")
    rows = []
    for row in score_table.itertuples():
        pairs = f"{row.positive_pairs}/{row.n_pairs}"
        values = (
            row.dataset, row.state, row.threshold, pairs, row.n_score_genes, format_number(row.mean_delta),
            format_ci(row), format_number(row.t_p_value, 5),
        )
        rows.append(values)
    column_names = [
        "Cohort", "State", "Cells per endpoint", "Pairs higher / total", "Genes",
        "Mean change", "95% CI", "Nominal paired t P",
    ]
    write_table(column_names, rows)
    write_paragraph(
        "These descriptive scores average 19 target-gene logCPM values from TMM-normalized, unfiltered state "
        "matrices (prior count = 1), without gene-wise z scoring. Changes are recurrent minus primary, with "
        "paired t-based 95% confidence intervals and nominal P values. They differ from the standardized pan-myeloid "
        "scores and are not competitive camera tests. Score availability does not imply that the filtered-gene "
        "coverage criterion for camera was met.",
        size=15,
    )

    # S5. The source table supplies one row per steroid trajectory or completeness group.
    steroid_dir = revision_dir / "step4_steroids_geomx"
    clinical_dir = steroid_dir / "clinical_case_normalized"
    trajectory_table = pd.read_csv(clinical_dir / "steroid_summary_all_statuses.csv")
    all_pairs = int(trajectory_table.loc[trajectory_table.group == "All", "n"].item())
    write_section(
        5, "Steroid trajectories and corticosteroid-related sensitivity analyses",
        "GSE274546 steroid status was classified as Yes or No at primary and recurrent sampling. Arrows indicate "
        "the transition from primary to recurrent sampling. The table includes four observed trajectories, "
        f"patients with incomplete records, patients with complete records, and all {all_pairs} patients. Group "
        "effects are descriptive paired score changes with t-based 95% CIs.",
    )
    rows = []
    for row in trajectory_table.itertuples():
        pairs = f"{row.positive}/{row.n}"
        rows.append((row.group, pairs, format_number(row.mean_delta), format_ci(row), format_number(row.paired_t_p, 5)))
    write_table(["Steroid trajectory / records", "Pairs higher / total", "Mean change", "95% CI", "Paired t P"], rows)

    camera_table = pd.read_csv(steroid_dir / "camera_results.csv")
    rows = []
    for _, row in camera_table.iterrows():
        rows.append((row.dataset, row.set_id, row["mode"], row.n_pairs, row.n_genes, format_number(row.p_value, 5)))
    write_table(["Cohort", "Set", "Mode", "Pairs", "Genes", "camera P"], rows)
    write_paragraph(
        "The gene sets include GOBP_RESPONSE_TO_CORTICOSTEROID and MG-Inflam20 without SGK1/PDK4. Camera with "
        "estimated inter-gene correlation is the primary analysis; a fixed correlation of 0.01 is evaluated in "
        "a sensitivity analysis.",
        size=15,
    )

    # The frozen source audit gives collection site and processing lab for each incomplete record.
    write_paragraph("Source information for patients with incomplete steroid records", True)
    audit_table = pd.read_csv(clinical_dir / "steroid_missing_source_audit.csv")
    rows = [
        (row.analysis_patient_id, row.source_patient_id, row.tissue_source, row.processing_laboratory)
        for row in audit_table.itertuples()
    ]
    write_table(["Analysis patient ID", "Source patient ID", "Tissue source", "Processing laboratory"], rows)

    center_counts = audit_table.tissue_source.value_counts().to_dict()
    if center_counts != {"Pitié-Salpêtrière Hospital": 6, "Duke University": 1, "Tokyo University": 1}:
        raise ValueError("S5 source counts do not match the frozen statement about incomplete records")
    center_table = pd.read_csv(clinical_dir / "GSE274546_center_summary.csv")
    center_table = center_table.set_index("center_group")
    paris = center_table.loc["Pitié-Salpêtrière Hospital"]
    other_centers = center_table.loc["Other centers"]
    write_paragraph(
        "Patient and tumor identifiers were matched to Supplementary Table 1 of the GBM-CARE source study. Tissue "
        "source and processing laboratory agreed between primary and recurrent samples for each patient. Six patients "
        "were from Pitié-Salpêtrière Hospital, one from Duke University, and one from Tokyo University. "
        f"Pitié-Salpêtrière Hospital contributed {int(paris.n_pairs)} of the {all_pairs} pairs; steroid records were "
        f"incomplete for {int(paris.n_incomplete_records)}/{int(paris.n_pairs)} patients at this center and "
        f"{int(other_centers.n_incomplete_records)}/{int(other_centers.n_pairs)} at the other centers. Mean paired "
        f"score changes were {paris.mean_paired_change:.2f} at this center and "
        f"{other_centers.mean_paired_change:.2f} at the other centers. These center summaries are descriptive. "
        "Lab1 and Lab2 denote processing laboratories, not sequencing batches. Sequencing-run or batch identifiers "
        "were not identified in the clinical and cell metadata examined.",
        size=15,
    )

    # S6. Preserve the three GeoMx contrast rows and the separate Welch statistic and degrees of freedom.
    write_section(
        6, "Paired GeoMx changes by nivolumab exposure",
        "Patient-level MG-Inflam20 scores in paired IBA1-positive regions. Within-group P values are paired t tests; "
        "the between-group difference uses Welch’s two-sample t test.",
    )
    geomx_table = pd.read_csv(steroid_dir / "GeoMx_subgroup_interaction.csv")
    rows = []
    for row in geomx_table.itertuples():
        pairs = "—" if pd.isna(row.positive) else f"{int(row.positive)}/{int(row.n)}"
        rows.append((row.contrast, pairs, format_number(row.estimate), format_ci(row), format_number(row.p_value, 5)))
    write_table(["Contrast", "Pairs higher / total", "Effect", "95% CI", "P"], rows)

    welch_table = pd.read_csv(steroid_dir / "GeoMx_Welch_statistic_df.csv")
    welch = welch_table.iloc[0]
    write_paragraph(
        f"Welch between-group contrast: t={welch.t:.3f}; df={welch.df:.2f}; n={welch.n_nivolumab} nivolumab and "
        f"{welch.n_control} controls. CI, confidence interval.",
        size=15,
    )

    # S7. The number of sign flips is a fixed setting of the October protein analysis.
    write_section(
        7, "Fixed protein-set and complete-case sensitivity analyses",
        "P values are two-sided paired sign-flip estimates based on 100,000 random configurations. The minimum "
        "attainable Monte Carlo P is 1.0×10⁻⁵.",
    )
    protein_table = pd.read_csv(revision_dir / "step5_protein_common/effects.csv")
    rows = []
    for _, row in protein_table.iterrows():
        pairs = f"{row.positive}/{row.n_pairs}"
        values = (
            row.analysis, row["mode"], pairs, row.n_proteins, format_number(row.mean_delta),
            format_ci(row), format_number(row.sign_flip_p, 5),
        )
        rows.append(values)
    column_names = [
        "Peptide matrix", "Analysis", "Pairs higher / total", "Proteins", "Mean change", "95% CI", "Sign-flip P"
    ]
    write_table(column_names, rows)

    fixed_unique = (protein_table.analysis == "unique_peptides") & (protein_table["mode"] == "fixed_proteins_all105")
    proteins = protein_table.loc[fixed_unique, "proteins"].item().split(";")
    if len(proteins) != 6:
        raise ValueError("S7 requires the October six-protein unique-peptide set")
    write_paragraph(f"Fixed six-protein unique-peptide set: {', '.join(proteins)}. CI, confidence interval.", size=15)

    # S8. Export the frozen model coefficients, influence diagnostics, covariates and enrichment tests.
    write_section(
        8, "Exploratory anti-PD-1 model and enrichment sensitivity analyses",
        "Patient-score models show exposure coefficients and t-based 95% CIs. In GSE121810, the original-covariate "
        "and unadjusted models coincide because additional baseline covariates were unavailable. In GSE154795, the "
        "original-covariate model includes MGMT promoter methylation status and technical covariates. The "
        "myeloid-expression proxy is defined in Supplementary Methods. Cell counts were unavailable for the bulk "
        "RNA-sequencing cohort GSE121810; n = 0 in these rows denotes no available observations, not an absence "
        "of myeloid cells.",
    )
    model_table = pd.read_csv(closeout_dir / "model_results_with_original_raw.csv")
    rows = []
    for row in model_table.itertuples():
        model_size = f"{row.parameters}/{row.residual_df}"
        values = (
            row.dataset, row.model, row.n, model_size, format_number(row.effect),
            format_ci(row), format_number(row.p_value, 5),
        )
        rows.append(values)
    write_table(["Cohort", "Model", "n", "Parameters / df", "Effect", "95% CI", "P"], rows)

    stability_dir = revision_dir / "step5_pd1_stability_final"
    omission_table = pd.read_csv(stability_dir / "leave_one_out.csv")
    rows = []
    for cohort, cohort_table in omission_table.groupby("dataset"):
        minimum = format_number(cohort_table.left_out_effect.min())
        maximum = format_number(cohort_table.left_out_effect.max())
        rows.append((cohort, f"{minimum} to {maximum}", len(cohort_table)))
    write_table(["Cohort", "Leave-one-out effects: minimum to maximum", "Omissions"], rows)

    cook_table = pd.read_csv(closeout_dir / "cook_threshold_patients.csv")
    rows = []
    for cohort, cohort_table in cook_table.groupby("dataset"):
        above_threshold = cohort_table.above_4_over_n
        patient_ids = ", ".join(cohort_table.loc[above_threshold, "sample"].astype(str)) or "None"
        rows.append((cohort, patient_ids, int(above_threshold.sum())))
    write_table(["Cohort", "Patients above Cook 4/n", "Count"], rows)

    vif_table = pd.read_csv(stability_dir / "vif.csv")
    rows = [(cohort, format_number(group_table.vif.max())) for cohort, group_table in vif_table.groupby("dataset")]
    write_table(["Cohort", "Largest VIF"], rows)

    covariate_table = pd.read_csv(stability_dir / "covariate_distributions.csv")
    rows = []
    for _, row in covariate_table.iterrows():
        average = row["median"] if pd.isna(row["mean"]) else row["mean"]
        rows.append((row.dataset, row.group, row.variable, row.level, row.n, format_number(average)))
    write_table(["Cohort", "Exposure group", "Covariate", "Level", "n", "Mean / median"], rows)

    camera_table = pd.read_csv(revision_dir / "step2_pd1_camera/camera_results.csv")
    rows = []
    for _, row in camera_table.iterrows():
        rows.append((row.dataset, row["mode"], row.n_samples, row.n_genes, format_number(row.p_value, 6)))
    write_table(["Cohort", "Camera mode", "Samples", "Genes", "P"], rows)

    # Report recorded fgsea statistics and the affected negative F values without re-ranking any genes.
    rank_dir = revision_dir / "step5_rank_unification_v2"
    for label, filename in (
        ("sign(logFC)×sqrt(F)", "sqrtF_fgsea.csv"), ("Historical sign(logFC)×F", "historical_F_fgsea.csv")
    ):
        enrichment_table = pd.read_csv(rank_dir / filename)
        rows = []
        for row in enrichment_table.itertuples():
            rows.append((row.dataset, label, row.size, format_number(row.NES), format_number(row.p_value, 7)))
        write_table(["Cohort", "Ranking", "Genes", "NES", "fgsea P"], rows)

    negative_table = pd.read_csv(rank_dir / "negative_F_genes.csv")
    negative_table = negative_table.loc[negative_table.dataset == "GSE154795"]
    minimum_f = f"{negative_table.F.min():.1f}".replace("-", "−")
    write_paragraph(
        "For the square-root ranking, negative F statistics with P=1 were set to zero "
        f"({len(negative_table)} genes in GSE154795; minimum F approximately {minimum_f}). The difference from "
        "historical results therefore includes both the square root and this explicit handling, not the square "
        "root alone.",
        size=15,
    )
    write_paragraph(
        "Abbreviations: CI, confidence interval; df, residual degrees of freedom; NES, normalized enrichment score; "
        "VIF, variance inflation factor. The historical F ranking is a sensitivity comparison.",
        size=15,
    )

    # Retain the landscape page, margins and original three-part Word package.
    page_size = {"w": 16838, "h": 11906, "orient": "landscape"}
    page_margins = {"top": 720, "right": 600, "bottom": 720, "left": 600, "header": 360, "footer": 360, "gutter": 0}
    section = add_element(body, "sectPr")
    add_element(section, "pgSz", page_size)
    add_element(section, "pgMar", page_margins)
    document_xml = etree.tostring(document, encoding="utf-8", xml_declaration=True)

    content_types = (
        b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        b'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        b'<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        b'<Default Extension="xml" ContentType="application/xml"/>'
        b'<Override PartName="/word/document.xml" '
        b'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        b'</Types>'
    )
    relationships = (
        b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        b'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        b'<Relationship Id="rId1" '
        b'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        b'Target="word/document.xml"/>'
        b'</Relationships>'
    )
    output_path = revision_dir / "supplementary_R1/Supplementary_Tables_S1-S8.docx"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check the saved package before replacing the final document in the output directory.
    with TemporaryDirectory(dir=output_path.parent) as temporary_dir:
        staged_path = Path(temporary_dir) / output_path.name
        with ZipFile(staged_path, "w", ZIP_DEFLATED) as package:
            package.writestr("[Content_Types].xml", content_types)
            package.writestr("_rels/.rels", relationships)
            package.writestr("word/document.xml", document_xml)

        if not staged_path.is_file() or staged_path.stat().st_size == 0:
            raise ValueError("Supplementary table export is missing or empty")
        with ZipFile(staged_path) as package:
            if package.testzip() is not None:
                raise ValueError("Supplementary table export contains an unreadable Word package part")
            saved_document = etree.fromstring(package.read("word/document.xml"))

        saved_size = saved_document.find(f".//{word_tag}pgSz")
        saved_margins = saved_document.find(f".//{word_tag}pgMar")
        if saved_size is None or saved_margins is None:
            raise ValueError("Supplementary table export is missing page settings")
        if not all(saved_size.get(word_tag + key) == str(value) for key, value in page_size.items()):
            raise ValueError("Supplementary table export has the wrong page size")
        if not all(saved_margins.get(word_tag + key) == str(value) for key, value in page_margins.items()):
            raise ValueError("Supplementary table export has the wrong page margins")

        # Each output table must contain its header and every source-derived row and cell.
        saved_tables = saved_document.findall(f".//{word_tag}tbl")
        if len(saved_tables) != len(expected_shapes):
            raise ValueError("Supplementary table export does not contain all tables")
        for saved_table, (row_count, column_count) in zip(saved_tables, expected_shapes):
            saved_rows = saved_table.findall(word_tag + "tr")
            valid_columns = all(len(row.findall(word_tag + "tc")) == column_count for row in saved_rows)
            if len(saved_rows) != row_count or not valid_columns:
                raise ValueError("Supplementary table export has an incorrect number of rows or cells")

        staged_path.replace(output_path)
    print(f"Supplementary Word tables: {len(expected_shapes)}")


if __name__ == "__main__":
    main()
