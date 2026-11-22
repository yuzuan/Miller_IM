#!/usr/bin/env python3
"""Write the initial eight-source overlap inventory.

Inputs are the recorded papers, accession records and patient identifiers below.
Outputs are step5_overlap_audit/source_provenance.csv, pairwise_overlap_classification.csv and checks.csv.
Run python3 NOA_R1/scripts/10a_initial_overlap_inventory.py <scratch_revision_dir>.
"""

import csv
import sys
from itertools import combinations
from pathlib import Path


if len(sys.argv) < 2:
    raise SystemExit("Provide the revision directory.")
output_dir = Path(sys.argv[1]).resolve() / "step5_overlap_audit"
output_dir.mkdir(exist_ok=False)

# Retain the original source order and report sampling years only when documented as sampling years.
sources_table = [
    dict(
        dataset="GSE174554",
        original_paper="Wang et al., Nature Cancer 2022; doi:10.1038/s43018-022-00475-x",
        contribution="UCSF Neurosurgery Tissue Bank",
        sampling_years="Individual dates unavailable. The paper describes decades of biobanking.",
        patient_ids="Sample IDs include PB0001. Revision pairing keys use GSE174554_pair*.",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554",
    ),
    dict(
        dataset="GSE274546",
        original_paper="Spitzer et al., Nature Genetics 2025; doi:10.1038/s41588-025-02168-4",
        contribution="CARE: MD Anderson, Duke, Tokyo, Pitié-Salpêtrière, St Michael’s, SNU, NORLUX",
        sampling_years="Sampling years for the whole cohort were unavailable from the paper.",
        patient_ids="Deidentified IDs include DU01 and FR01. Prefixes alone do not establish collection sites.",
        evidence_url="https://www.nature.com/articles/s41588-025-02168-4",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546",
    ),
    dict(
        dataset="GSE278456",
        original_paper="Jackson et al., Distinct myeloid-derived suppressor cell populations in human glioblastoma; "
        "PMID:39818911",
        contribution="Johns Hopkins Hospital. GEO also lists a University of Pennsylvania submission contact.",
        sampling_years="Individual sampling years unavailable.",
        patient_ids="GBM006, GBM030, GBM049 and other IDs",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12836367/",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278456",
    ),
    dict(
        dataset="GSE276841",
        original_paper="Jackson et al., Distinct myeloid-derived suppressor cell populations in human glioblastoma; "
        "PMID:39818911",
        contribution="Johns Hopkins Hospital",
        sampling_years="Individual sampling years unavailable.",
        patient_ids="Spatial samples from GBM030 and GBM049.",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12836367/",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276841",
    ),
    dict(
        dataset="GeoMx_Zenodo16839828",
        original_paper="Artzi et al., Acta Neuropathologica 2025; doi:10.1007/s00401-025-02937-9",
        contribution="Copenhagen University Hospital Rigshospitalet",
        sampling_years="Diagnoses were in 2015–2021. Dates of the two samples per patient are unavailable.",
        patient_ids="Numeric Patient_No identifiers. Revision figure IDs include P01.",
        evidence_url="https://link.springer.com/article/10.1007/s00401-025-02937-9",
        metadata_url="https://zenodo.org/records/16839828",
    ),
    dict(
        dataset="PDC000514",
        original_paper="Integrated Proteogenomic Characterization of Glioblastoma Evolution; PMCID:PMC10939876",
        contribution="NCC, SNU, SMC and CNU in Korea and Pitié-Salpêtrière in France",
        sampling_years="Individual sampling years unavailable.",
        patient_ids="IDs include KNCC_GBM0001 with T1/T2 suffixes.",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10939876/",
        metadata_url="https://proteomic.datacommons.cancer.gov/pdc/study/PDC000514",
    ),
    dict(
        dataset="GSE121810",
        original_paper="Cloughesy et al., Nature Medicine 2019; doi:10.1038/s41591-018-0337-7",
        contribution="UCLA clinical trial and tumor samples",
        sampling_years="Individual sampling years unavailable.",
        patient_ids="IDs include Pt3_A.",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC6408961/",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810",
    ),
    dict(
        dataset="GSE154795",
        original_paper="Lee et al., Nature Communications 2021; PMCID:PMC8626557",
        contribution="UCLA clinical study and tumor samples",
        sampling_years="Individual sampling years unavailable.",
        patient_ids="IDs include LB3018T.",
        evidence_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC8626557/",
        metadata_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154795",
    ),
]
source_path = output_dir / "source_provenance.csv"
with source_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(sources_table[0]))
    writer.writeheader()
    writer.writerows(sources_table)

# The historical review classified only the documented Jackson cross-modality pair as known overlap.
known_overlap = "known overlap"
uncertain_overlap = "insufficient information to determine overlap"
pairs_table = []
for source_a, source_b in combinations(sources_table, 2):
    pair = {source_a["dataset"], source_b["dataset"]}
    if pair == {"GSE278456", "GSE276841"}:
        classification = known_overlap
        evidence = "Both GEO records list GBM030 and GBM049. One paper reports the single-cell and spatial experiments."
        url_a = source_a["metadata_url"]
        url_b = source_b["metadata_url"]
        limitation = ""
    else:
        classification = uncertain_overlap
        evidence = "Reviewed the institutions and deidentified IDs. Patient-level mappings and full dates were absent."
        url_a = source_a["evidence_url"]
        url_b = source_b["evidence_url"]
        if pair == {"GSE121810", "GSE154795"}:
            limitation = "Both studies include UCLA. Shared patients remain possible."
        else:
            limitation = "Different identifiers or institutions alone cannot establish independent patients."
    pairs_table.append(
        dict(
            source_a=source_a["dataset"],
            source_b=source_b["dataset"],
            classification=classification,
            evidence=evidence,
            evidence_url_a=url_a,
            evidence_url_b=url_b,
            limitation=limitation,
        )
    )
pair_path = output_dir / "pairwise_overlap_classification.csv"
with pair_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(pairs_table[0]))
    writer.writeheader()
    writer.writerows(pairs_table)

# Check the full set of eight sources and its 28 unique comparisons.
unique_pairs = {(row["source_a"], row["source_b"]) for row in pairs_table}
checks = [
    ("eight_sources", len(sources_table) == 8),
    ("28_pairs", len(pairs_table) == 28),
    ("unique_pairs", len(unique_pairs) == 28),
    ("known_JHU_cross_modal_overlap", sum(row["classification"] == known_overlap for row in pairs_table) == 1),
    (
        "all_sources_have_original_urls",
        all(row["evidence_url"].startswith("https://") and row["metadata_url"].startswith("https://")
            for row in sources_table),
    ),
    (
        "no_year_inferred_from_publication",
        all("2025" not in row["sampling_years"] and "2022" not in row["sampling_years"] for row in sources_table),
    ),
]
checks_path = output_dir / "checks.csv"
with checks_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["check", "passed"])
    writer.writerows(checks)

for output_path in (source_path, pair_path, checks_path):
    if output_path.stat().st_size == 0:
        raise OSError(f"Empty provenance report {output_path.name}")
failed = sum(not passed for _, passed in checks)
print("OVERLAP_AUDIT", len(checks), "failed", failed)
if failed:
    raise SystemExit(1)
