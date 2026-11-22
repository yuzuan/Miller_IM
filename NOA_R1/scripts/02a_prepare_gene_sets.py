#!/usr/bin/env python3
"""Verify the historical five gene sets against the original published workbooks.

Inputs: sources/Miller_2025_Supplementary_Table_2.xlsx, the Abdelfattah Supplementary Data 9 workbook,
    sources/Sankowski_2019_Supplementary_Tables.xlsx and step2_extended_specificity/locked_gene_sets_used.csv.
    Use NOA_ABDEL_WORKBOOK for the downloaded Abdelfattah workbook. The original Miller workbook is also
    supported in literature_markers/ under NOA_LEGACY_ROOT or its parent.
Outputs: step2_gene_sets.csv and the original Miller and Abdelfattah workbook copies under sources/.
Run: python NOA_R1/scripts/02a_prepare_gene_sets.py
"""

import csv
import importlib
import os
import shutil
from pathlib import Path

import openpyxl
import pandas as pd


# The published order and comparator extraction rules are checked against the archived gene definitions.
revision = Path(__file__).resolve().parents[1]
source_dir = revision / "sources"
miller_copy = source_dir / "Miller_2025_Supplementary_Table_2.xlsx"
miller_path = miller_copy
if not miller_path.is_file():
    revision_tools = importlib.import_module("00a_functions")
    project = revision_tools.legacy_root(revision)
    workbook_name = "Miller_2025_Nature_myeloid_glioma_Supplementary_Table_2.xlsx"
    miller_path = project / "literature_markers" / workbook_name
    if not miller_path.is_file():
        miller_path = project.parent / "literature_markers" / workbook_name

# The Miller source is a ranked workbook. Read the fixed inflammatory column in author order.
workbook = openpyxl.load_workbook(miller_path, read_only=True, data_only=True)
try:
    workbook_rows = list(workbook["Myeloid Top 100 genes per prog."].values)
finally:
    workbook.close()
inflammatory_column = workbook_rows[0].index("Microglial Inflammatory")
target_genes = [row[inflammatory_column] for row in workbook_rows[1:101]]

# The original Abdelfattah workbook names the i-Mic cluster MC01 and has its headers on the second row.
abdel_copy = source_dir / "Abdelfattah_2022_Supplementary_Data_9.xlsx"
abdel_path = Path(os.environ.get("NOA_ABDEL_WORKBOOK") or abdel_copy).expanduser()
abdel_table = pd.read_excel(abdel_path, sheet_name="Sheet1", header=1)
imic_table = abdel_table.loc[abdel_table["cluster"] == "MC01"]
if len(imic_table) != 85:
    raise ValueError("Abdelfattah MC01 source row count differs from the original workbook")
ranked_imic_table = imic_table.sort_values("avg_log2FC", ascending=False, kind="stable")
imic_genes = ranked_imic_table["gene"].head(50).tolist()
if imic_genes != imic_table["gene"].head(50).tolist():
    raise ValueError("Abdelfattah author row order and fold-change rank disagree")

# Keep every original Cluster 1 row from Sankowski, including its compound feature.
sankowski_path = source_dir / "Sankowski_2019_Supplementary_Tables.xlsx"
sankowski_table = pd.read_excel(sankowski_path, sheet_name="Suppl. table 3")
sankowski_genes = sankowski_table.loc[sankowski_table["Cluster"] == 1, "GENEID"].tolist()
if "CCL3L3|CCL3" not in sankowski_genes:
    raise ValueError("Sankowski Cluster 1 is missing the original compound feature")

# Match the exact archived symbol and rank sequence. These CSV rows replace the private ranked catalog check.
set_genes = {
    "target_top20": target_genes[:20],
    "target_top50": target_genes[:50],
    "target_top100": target_genes,
    "Abdelfattah_iMic_top50": imic_genes,
    "Sankowski_C1": sankowski_genes,
}
archive_path = revision / "step2_extended_specificity" / "locked_gene_sets_used.csv"
with archive_path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    definition_table = [row for row in reader if row["set_id"] in set_genes]

if not fieldnames or list(dict.fromkeys(row["set_id"] for row in definition_table)) != list(set_genes):
    raise ValueError("Archived five-set order does not match the published definitions")
for set_id, genes in set_genes.items():
    member_table = [row for row in definition_table if row["set_id"] == set_id]
    archived_members = [(int(row["original_rank"]), row["gene"]) for row in member_table]
    if archived_members != list(enumerate(genes, 1)) or len(genes) != len(set(genes)):
        raise ValueError(f"Workbook genes or order differ from the archived set {set_id}")
    if not all(isinstance(gene, str) and gene and gene.strip() == gene for gene in genes):
        raise ValueError(f"Invalid source gene symbols in {set_id}")

# Keep the existing filenames and portable source fields stored with the archived definitions.
source_dir.mkdir(parents=True, exist_ok=True)
for source_path, destination in ((miller_path, miller_copy), (abdel_path, abdel_copy)):
    if source_path.resolve() != destination.resolve():
        shutil.copyfile(source_path, destination)
        if destination.stat().st_size != source_path.stat().st_size:
            raise ValueError(f"Incomplete workbook copy: {destination.name}")

output_path = revision / "step2_gene_sets.csv"
with output_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(definition_table)
with output_path.open(encoding="utf-8", newline="") as handle:
    written_table = list(csv.DictReader(handle))
if written_table != definition_table:
    raise ValueError("Written five-set definitions differ from the verified rows")
print(f"Validated {len(set_genes)} gene sets and {len(definition_table)} definition rows.")
