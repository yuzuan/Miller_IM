#!/usr/bin/env python3
"""Verify the Miller systemic program and extend the five historical sets.

Inputs: sources/Miller_2025_Supplementary_Table_2.xlsx (Systemic Inflammatory column),
    step2_gene_sets.csv and step2_extended_specificity/locked_gene_sets_used.csv.
Output: step2_locked_gene_sets.csv, in the original six-set order.
Run: python NOA_R1/scripts/02b_prepare_prespecified_sets.py after 02a.
"""

import csv
from pathlib import Path

import openpyxl


# Retain the five definitions from 02a, including their original ranks.
revision = Path(__file__).resolve().parents[1]
with (revision / "step2_gene_sets.csv").open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    definition_table = list(reader)
if not fieldnames:
    raise ValueError("Missing historical five-set table header")

# Compare the prepared base sets with the archived rows, then check the extra systemic set from the workbook.
archive_path = revision / "step2_extended_specificity" / "locked_gene_sets_used.csv"
with archive_path.open(encoding="utf-8", newline="") as handle:
    archive_table = list(csv.DictReader(handle))
archived_base_table = [row for row in archive_table if row["set_id"] != "Miller_systemic_top20"]
if definition_table != archived_base_table:
    raise ValueError("Prepared five-set definitions differ from the archived membership and metadata")

source_path = revision / "sources" / "Miller_2025_Supplementary_Table_2.xlsx"
workbook = openpyxl.load_workbook(source_path, read_only=True, data_only=True)
try:
    workbook_rows = list(workbook["Myeloid Top 100 genes per prog."].values)
finally:
    workbook.close()
systemic_column = workbook_rows[0].index("Systemic Inflammatory")
systemic_genes = [row[systemic_column] for row in workbook_rows[1:21]]

# The comparator is the 20 author-ranked genes as recorded in the six-set archive.
systemic_table = [row for row in archive_table if row["set_id"] == "Miller_systemic_top20"]
archived_members = [(int(row["original_rank"]), row["gene"]) for row in systemic_table]
if archived_members != list(enumerate(systemic_genes, 1)) or len(set(systemic_genes)) != 20:
    raise ValueError("Workbook systemic gene order differs from the archived comparator")
if not all(isinstance(gene, str) and gene for gene in systemic_genes):
    raise ValueError("Invalid gene symbol in the Miller systemic comparator")

# Write base rows followed by the original systemic comparator rows and source metadata.
definition_table_total = definition_table + systemic_table
output_path = revision / "step2_locked_gene_sets.csv"
with output_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(definition_table_total)
with output_path.open(encoding="utf-8", newline="") as handle:
    written_table = list(csv.DictReader(handle))
if written_table != definition_table_total:
    raise ValueError("Written six-set definitions differ from the verified rows")
set_count = len({row["set_id"] for row in definition_table_total})
print(f"Validated {set_count} gene sets and {len(definition_table_total)} definition rows.")
