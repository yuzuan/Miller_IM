#!/usr/bin/env python3
"""Describe the completeness of paired steroid records by tissue source.

Inputs: step4_steroids_geomx/patient_deltas.csv and
step4_steroids_geomx/clinical_case_normalized/GSE274546_steroid_patient_status_normalized.csv.
Use the original GBM-CARE Supplementary Table 1 workbook (Patient ID, Tumor ID, Tissue source).
Outputs in step4_steroids_geomx/clinical_case_normalized: GSE274546_center_patient_mapping.csv and
GSE274546_center_summary.csv.
Run with the original Python openpyxl reader:
python NOA_R1/scripts/04c_steroid_center_summary.py --root NOA_R1 --clinical-workbook FILE
"""

import argparse
import csv
import math
from pathlib import Path
from statistics import mean

import openpyxl


def summarize(root, workbook):
    """Validate the existing clinical mapping and return patient and center rows after writing both CSVs."""
    step_dir = root / "step4_steroids_geomx"
    output_dir = step_dir / "clinical_case_normalized"
    with (output_dir / "GSE274546_steroid_patient_status_normalized.csv").open(encoding="utf-8") as handle:
        patient_rows = list(csv.DictReader(handle))

    # Use the formal target_top20 differences as the independent score reference.
    with (step_dir / "patient_deltas.csv").open(encoding="utf-8") as handle:
        reference_rows = csv.DictReader(handle)
        deltas = {
            row["patient_id"]: float(row["delta"])
            for row in reference_rows
            if row["dataset"] == "GSE274546" and row["set_id"] == "target_top20"
        }

    # Read the source workbook in place. The original analysis used its active worksheet.
    source_workbook = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    try:
        workbook_rows = list(source_workbook.active.values)
    finally:
        source_workbook.close()

    columns = [str(value).strip() for value in workbook_rows[0]]
    clinical_rows = []
    for row in workbook_rows[1:]:
        if len(row) != len(columns):
            raise ValueError("Clinical workbook row width differs from its header")
        clinical_rows.append(dict(zip(columns, row)))

    # Check the frozen cohort size and retain one row per formal pair.
    patient_ids = {row["patient_id"] for row in patient_rows}
    if not len(patient_rows) == len(patient_ids) == len(deltas) == 45:
        raise ValueError("The center summary requires the 45 formal GBM-CARE pairs")
    mapping_rows = []
    for row in patient_rows:
        patient_id = row["patient_id"]
        prefix = row["tumor_prefix"]
        tumor_ids = (prefix + "_0", prefix + "_1")
        matched = [
            source for source in clinical_rows
            if source["Patient ID"] == row["clinical_patient_id"] and source["Tumor ID"] in tumor_ids
        ]
        if len(matched) != 2 or len({source["Tissue source"] for source in matched}) != 1:
            raise ValueError(f"Expected one tissue source across both tumors for {patient_id}")

        delta = float(row["delta"])
        if not math.isclose(delta, deltas[patient_id], abs_tol=1e-12, rel_tol=0):
            raise ValueError(f"Paired score change differs from the frozen reference for {patient_id}")
        mapping_row = dict(
            analysis_patient_id   = patient_id,
            source_patient_id     = row["clinical_patient_id"],
            tumor_prefix          = prefix,
            tissue_source         = matched[0]["Tissue source"],
            steroid_record_status = "Incomplete" if row["group"] == "Missing" else "Complete",
            primary_steroid       = row["primary_steroid"],
            recurrent_steroid     = row["recurrent_steroid"],
            paired_score_change   = delta
        )
        mapping_rows.append(mapping_row)

    # Preserve the original source label and collapse every other site into the existing Other centers group.
    index_center = "Pitié-Salpêtrière Hospital"
    groups = {
        index_center: [row for row in mapping_rows if row["tissue_source"] == index_center],
        "Other centers": [row for row in mapping_rows if row["tissue_source"] != index_center],
    }
    summary_rows = []
    for name, rows in groups.items():
        complete = [row for row in rows if row["steroid_record_status"] == "Complete"]
        incomplete = [row for row in rows if row["steroid_record_status"] == "Incomplete"]
        summary_row = dict(
            center_group         = name,
            n_pairs              = len(rows),
            n_incomplete_records = len(incomplete),
            n_complete_records   = len(complete),
            mean_paired_change   = mean(row["paired_score_change"] for row in rows),
            mean_change_complete_records   = mean(row["paired_score_change"] for row in complete),
            mean_change_incomplete_records = mean(row["paired_score_change"] for row in incomplete)
        )
        summary_rows.append(summary_row)

    # Export genuine one-patient and one-center rows with their original column order.
    outputs = (
        ("GSE274546_center_patient_mapping.csv", mapping_rows),
        ("GSE274546_center_summary.csv", summary_rows),
    )
    for name, rows in outputs:
        with (output_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f"Verified {len(mapping_rows)} paired source mappings and {len(deltas)} score changes.")
    result = (mapping_rows, summary_rows)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical-workbook", type=Path, required=True, help="GBM-CARE Supplementary Table 1 workbook")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    summarize(args.root, args.clinical_workbook)
