#!/usr/bin/env python3
"""Check the initial inventory and its documented Jackson cross-modality overlap.

Inputs are step5_overlap_audit/source_provenance.csv and pairwise_overlap_classification.csv.
Output is step5_overlap_audit/independent_checks.csv.
Run python3 NOA_R1/scripts/10b_check_initial_overlap_inventory.py <scratch_revision_dir> after 10a.
"""

import csv
import sys
from itertools import combinations
from pathlib import Path


# Read the source and pair records from the initial inventory.
if len(sys.argv) < 2:
    raise SystemExit("Provide the revision directory.")
output_dir = Path(sys.argv[1]).resolve() / "step5_overlap_audit"
with (output_dir / "source_provenance.csv").open(encoding="utf-8", newline="") as handle:
    sources_table = list(csv.DictReader(handle))

with (output_dir / "pairwise_overlap_classification.csv").open(encoding="utf-8", newline="") as handle:
    pairs_table = list(csv.DictReader(handle))

# Show exhaustive coverage and the single pair with linked patient IDs.
source_ids = [row["dataset"] for row in sources_table]
expected_pairs = set(combinations(source_ids, 2))
recorded_pairs = {(row["source_a"], row["source_b"]) for row in pairs_table}
known_table = [row for row in pairs_table if row["classification"] == "known overlap"]
known_jackson_pair = False
if len(known_table) == 1:
    known_record = known_table[0]
    known_sources = {known_record["source_a"], known_record["source_b"]}
    known_jackson_pair = known_sources == {"GSE278456", "GSE276841"} and all(
        patient_id in known_record["evidence"] for patient_id in ("GBM030", "GBM049")
    )

checks = [
    ("independent_eight_sources", len(sources_table) == 8 and len(set(source_ids)) == 8),
    ("independent_exhaustive_28_pairs", len(pairs_table) == 28 and expected_pairs == recorded_pairs),
    ("independent_GBM030_GBM049_cross_modal", known_jackson_pair),
    (
        "independent_no_unsupported_nonoverlap",
        not any(row["classification"] == "evidence supports no overlap" for row in pairs_table),
    ),
    (
        "independent_27_uncertain",
        sum(row["classification"] == "insufficient information to determine overlap" for row in pairs_table) == 27,
    ),
]

# Record the checks in their original order and retain the exclusive-write contract.
checks_path = output_dir / "independent_checks.csv"
with checks_path.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["check", "passed"])
    writer.writerows(checks)
if checks_path.stat().st_size == 0:
    raise OSError("Empty initial overlap checks report.")

failed = sum(not passed for _, passed in checks)
print("OVERLAP_INDEPENDENT", len(checks), "failed", failed)
if failed:
    raise SystemExit(1)
