#!/usr/bin/env python3
"""Validate the published revision scripts and derived source tables.

Inputs: provenance/step1_input_manifest.csv, the archived CAMERA and figure source tables, and release files.
Output: a validation summary and an optional report at the requested JSON path.
Run: python NOA_R1/scripts/11_validate_revision_package.py [directory] [--json report.json]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument("--json", type=Path, help="Write a machine-readable validation report.")
args = parser.parse_args()
package_root = args.directory.resolve(strict=True)
if not package_root.is_dir():
    raise NotADirectoryError(package_root)

# The archived provenance manifest also defines the expected pair counts.
manifest_path = package_root / "provenance" / "step1_input_manifest.csv"
with manifest_path.open(encoding="utf-8", newline="") as handle:
    manifest_table = list(csv.DictReader(handle))
expected_counts = {}
for row in manifest_table:
    expected_counts[row["dataset"]] = int(row["expected_pairs"])
if len(expected_counts) != len(manifest_table) or len(expected_counts) != 2:
    raise ValueError("The provenance manifest must contain two distinct recurrence cohorts.")

camera_path = package_root / "step1_primary_recurrence" / "camera_results.csv"
with camera_path.open(encoding="utf-8", newline="") as handle:
    camera_table = list(csv.DictReader(handle))
observed_counts = {}
for row in camera_table:
    if row["role"] == "primary":
        if row["dataset"] in observed_counts:
            raise ValueError(f"Duplicate primary CAMERA result for {row['dataset']}.")
        observed_counts[row["dataset"]] = int(row["n_pairs"])
if observed_counts != expected_counts:
    raise ValueError("The primary camera table does not match the provenance pair counts.")

# Confirm that each main figure has the archived panel manifest and source tables.
figure_count = 0
panel_count = 0
for figure_id in range(1, 6):
    figure_dir = package_root / "figures_R1" / "_font_build" / f"fig{figure_id}"
    data_dir = figure_dir / "data" / f"Figure{figure_id}"
    figure_manifest = data_dir / f"Figure{figure_id}_panel_manifest.csv"
    source_paths = list((data_dir / "source_data").glob("*.csv*"))
    if not figure_manifest.is_file() or not source_paths:
        raise FileNotFoundError(f"Missing manifest or source tables for Figure {figure_id}.")

    with figure_manifest.open(encoding="utf-8", newline="") as handle:
        panel_table = list(csv.DictReader(handle))
    if not panel_table:
        raise ValueError(f"Empty panel manifest for Figure {figure_id}.")

    for row in panel_table:
        for source in row["source"].split(";"):
            source_path = (package_root / source.strip()).resolve()
            if not source_path.is_relative_to(package_root) or not source_path.is_file():
                raise FileNotFoundError(f"Missing package source for Figure {figure_id} panel {row['panel']}.")
        panel_count += 1
    figure_count += 1

# Check authored filenames and UTF-8 source text. Scientific symbols are permitted.
han_text = re.compile(r"[\u3400-\u9fff\uf900-\ufaff\U00020000-\U000323af]")
local_path = re.compile(
    r"(?<![\w/])/(?:Users|Volumes|home|(?:private/)?tmp)/|"
    r"(?<![\w/])[A-Za-z]:[/\\](?!Windows[/\\]Fonts[/\\])|"
    r"\\\\[\w.-]+\\\w",
    re.IGNORECASE
)
text_suffixes = {".csv", ".tsv", ".json", ".txt", ".md", ".r", ".py", ".sh"}
cache_directories = {".git", ".ruff_cache", "__pycache__"}
report_path = args.json.resolve() if args.json is not None else None
file_count = 0
table_count = 0
issues = []
for path in sorted(package_root.rglob("*")):
    relative_path = path.relative_to(package_root)
    if not path.is_file() or not cache_directories.isdisjoint(relative_path.parts) or path == report_path:
        continue
    file_count += 1
    if han_text.search(relative_path.as_posix()):
        issues.append(f"Non-English filename: {relative_path}")

    suffix = path.suffix.lower()
    if suffix == ".gz":
        suffix = Path(path.stem).suffix.lower()
    if suffix not in text_suffixes and path.name != ".gitignore":
        continue

    if path.suffix == ".gz":
        with gzip.open(path, encoding="utf-8", mode="rt") as handle:
            source_text = handle.read()
    else:
        source_text = path.read_text(encoding="utf-8")
    if han_text.search(source_text):
        issues.append(f"Non-English text: {relative_path}")
    if local_path.search(source_text):
        issues.append(f"Local absolute path: {relative_path}")

    if suffix in {".csv", ".tsv"}:
        delimiter = "," if suffix == ".csv" else "\t"
        reader = csv.reader(io.StringIO(source_text, newline=""), delimiter=delimiter, strict=True)
        header = next(reader, [])
        if not header:
            issues.append(f"Missing table header: {relative_path}")
        table_count += 1

report = {
    "files": file_count,
    "tables": table_count,
    "figures": figure_count,
    "panels": panel_count,
    "paired_cohorts": observed_counts,
    "issues": issues,
}
if report_path is not None:
    report_text = json.dumps(report, indent=2) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    if report_path.read_text(encoding="utf-8") != report_text:
        raise ValueError("Written validation report differs from memory.")

print(
    f"R1 validation: {file_count} files, {table_count} tables, "
    f"{figure_count} figures, {panel_count} panels, {len(issues)} issues"
)
if issues:
    for issue in issues:
        print(issue)
    raise SystemExit(1)
