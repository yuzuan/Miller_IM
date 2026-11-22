#!/usr/bin/env python3
"""Write the October patient-overlap reassessment from the recorded source evidence.

Inputs are the eight source publications, accession records and pair-specific findings embedded below.
Outputs beside the scripts directory are step5_closeout_v2/overlap_source_evidence.csv,
overlap_28_pair_reaudit.csv and overlap_access_attempts.csv. Existing files are protected.
Run python3 NOA_R1/scripts/05f_overlap_reassessment.py from a scratch copy of the revision.
"""

import csv
from collections import Counter
from itertools import combinations
from pathlib import Path


# The recorded sources retain their accession order and the original citation locators.
source_table = {
    "GSE174554": dict(
        paper="Wang et al., Nature Cancer 2022",
        url="https://www.nature.com/articles/s43018-022-00475-x",
        locator="Methods > Tumor tissue acquisition; all specimens UCSF Neurosurgery Tissue Bank, "
        "decades of biobanking",
        sites="UCSF",
        dates="Individual sampling years unavailable. The publication describes decades of biobanking.",
        ids="Deidentified GEO PB and Pt IDs. No cross-cohort mapping.",
        metadata="https://www.ncbi.nlm.nih.gov/bioproject/730532",
        metadata_locator="BioProject PRJNA730532 > GEO GSE174554; UCSF Diaz Lab",
    ),
    "GSE274546": dict(
        paper="Spitzer et al., Nature Genetics 2025; companion Nomura et al. generated GEO",
        url="https://www.nature.com/articles/s41588-025-02168-4",
        locator="Methods > Human participants and ethical approval, 7 sites; "
        "Results > Longitudinal profiling, CARE 59; Results > independent cohort ref.14",
        sites="MD Anderson|Duke|Tokyo University|Pitié-Salpêtrière|St Michael's Toronto|SNU|NORLUX/CHL",
        dates="Individual sampling years unavailable. GEO submission dates refer to the deposit.",
        ids="Deidentified IDs include DU01. GSM8452584 is DU01 T1.",
        metadata="https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8452584",
        metadata_locator="GSM8452584 > Source name / Characteristics: DU01 Primary",
    ),
    "GSE278456": dict(
        paper="Jackson et al., Science 2025",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12836367/",
        locator="Methods > Human samples: prospectively enrolled at Johns Hopkins Hospital; "
        "GEO indexed record GSE278456 samples GBM030/GBM049",
        sites="Johns Hopkins Hospital",
        dates="Individual sampling years unavailable. GEO submission dates refer to the deposit.",
        ids="GBM030, GBM049 and other IDs",
        metadata="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278456",
        metadata_locator="GEO Samples > GSM8546885 GBM030; GSM8546893 GBM049",
    ),
    "GSE276841": dict(
        paper="Jackson et al., Science 2025 (same paper)",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12836367/",
        locator="Methods > Human samples; GEO GSE276841 Summary: two spatial patients, untreated before resection",
        sites="Johns Hopkins Hospital",
        dates="Individual sampling years unavailable.",
        ids="GBM030, GBM049",
        metadata="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276841",
        metadata_locator="GEO Samples > GSM8506649 GBM030; GSM8506650 GBM049",
    ),
    "GeoMx_Zenodo16839828": dict(
        paper="Artzi et al., Acta Neuropathologica 2025",
        url="https://link.springer.com/article/10.1007/s00401-025-02937-9",
        locator="Materials and methods > Patient selection: CA209-9UP/NCT03890952 plus 10 local controls; all treated "
        "Rigshospitalet",
        sites="Copenhagen University Hospital Rigshospitalet",
        dates="Diagnoses in 2015–2021. Dates of the two samples per patient are unavailable.",
        ids="Deidentified publication IDs and Zenodo Figure 3A mapping IDs.",
        metadata="https://zenodo.org/records/16839828",
        metadata_locator="Data record for article; no cross-study patient key",
    ),
    "PDC000514": dict(
        paper="Kim et al., Cancer Cell 2024",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10939876/",
        locator="Key resources table > Human GBM tumor tissue samples and blood samples; Fig.1A 123 longitudinal pairs",
        sites="NCC Korea|SNU|Samsung Medical Center|Chonnam National University|Pitié-Salpêtrière",
        dates="Individual sampling years unavailable in the paper or PDC JSON. Table S1 unavailable at review.",
        ids="PDC case_submitter_id and KNCC_GBM IDs. No crosswalk to CARE.",
        metadata="https://proteomic.datacommons.cancer.gov/pdc/study/PDC000514",
        metadata_locator="PDC study page only JS shell; existing public clinical/biospecimen JSON lacks sampled year",
    ),
    "GSE121810": dict(
        paper="Cloughesy et al., Nature Medicine 2019; NCT02852655",
        url="https://escholarship.org/content/qt7gf847t7/qt7gf847t7_noSplash_e1f45a73ebcbea28921e3801ba049f36.pdf",
        locator="Results > Trial patient characteristics, PDF p.3: 35 randomized Oct 2016–Sep 2017 at 7 institutions; "
        "Fig.2 caption Pt A/B",
        sites="multi-site US; NCT02852655 later stage Methods lists DFCI, Huntsman, MD Anderson, MGH, MSK, UCLA, UCSF",
        dates="Trial enrollment October 2016 to September 2017. Individual sampling dates unavailable.",
        ids="Pt#_A/B (trial treatment arm)",
        metadata="https://www.nature.com/articles/s41467-024-54326-7",
        metadata_locator="Follow-on trial article Methods > Study design and patients: same NCT, 7 IRB sites; "
        "stage 2 roster may not prove stage 1 per-patient site",
    ),
    "GSE154795": dict(
        paper="Lee et al., Nature Communications 2021",
        url="https://www.nature.com/articles/s41467-021-26940-2",
        locator="Results > Patient characteristics: 69 all UCLA, 20 neo-aPD1; "
        "Methods > Human samples: recurrent patients off-label/off-trial pembrolizumab",
        sites="UCLA surgical cohort",
        dates="Individual surgery and enrollment dates unavailable. GEO was submitted in 2020.",
        ids="LBxxxxT",
        metadata="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154795",
        metadata_locator="GEO Samples > LB patient identifiers and anti-pd1-therapy; no NCT or Pt crosswalk",
    ),
}

# Use the released English category labels. Collection provenance supports the stated inference.
known_overlap = "known overlap"
possible_overlap = "possible overlap due to shared sources"
separate_sources = "evidence supports no overlap"
pair_evidence_table = {
    ("GSE174554", "GSE274546"): (
        separate_sources,
        "The 59 CARE patients came from seven sites outside UCSF. CARE used Wang 2022 as external validation.",
        "CARE Results independent cohort ref.14; Methods seven institutions; Wang Methods UCSF",
        "Individual sampling years are unavailable. The conclusion concerns the GSE274546 core cohort.",
    ),
    ("GSE174554", "GSE278456"): (
        separate_sources,
        "Wang used UCSF surgical tissue. Jackson prospectively enrolled at Johns Hopkins for GSE278456.",
        "Wang Methods tissue acquisition; Jackson Methods enrollment",
        "The collection sites differ. Individual identities across sites were not compared.",
    ),
    ("GSE174554", "GSE276841"): (
        separate_sources,
        "Wang used UCSF tissue. Jackson collected the GBM030 and GBM049 spatial samples at Johns Hopkins.",
        "Wang Methods tissue acquisition; Jackson Methods; GSE276841 GEO Samples",
        "Individual identities across sites were not compared.",
    ),
    ("GSE174554", "GeoMx_Zenodo16839828"): (
        separate_sources,
        "Wang used UCSF specimens. Artzi's NCT03890952 and local control tissues came from Rigshospitalet.",
        "Wang Methods tissue acquisition; Artzi Methods Patient selection",
        "2015–2021 are diagnosis years. Individual sampling years are unavailable.",
    ),
    ("GSE174554", "PDC000514"): (
        separate_sources,
        "Wang used UCSF specimens. Kim's five contributing hospitals were in Korea and Paris.",
        "Wang Methods tissue acquisition; Kim Key resources table",
        "The inference uses collection sites. Full dates are unavailable from the PDC portal.",
    ),
    ("GSE174554", "GSE121810"): (
        possible_overlap,
        "UCSF contributed to NCT02852655 and supplied Wang's tissue bank. Trial dates fall within its decades.",
        "Wang Methods tissue acquisition; Cloughesy PDF Results trial dates; UCLA NCT page Locations; same-NCT stage2 "
        "Methods",
        "Only site and time overlap is plausible. Trial Pt IDs lack a crosswalk to PB IDs.",
    ),
    ("GSE174554", "GSE154795"): (
        separate_sources,
        "Wang used UCSF surgical specimens. Lee's 69 patients had surgery at UCLA, also listed in GEO.",
        "Wang Methods tissue acquisition; Lee Results patient characteristics",
        "The inference uses collection sites. Individual identities across sites were not compared.",
    ),
    ("GSE274546", "GSE278456"): (
        separate_sources,
        "The seven CARE sites exclude Johns Hopkins, where Jackson enrolled the prospective cohort.",
        "CARE Methods seven institutions; Jackson Methods enrollment",
        "Inference based on collection sites.",
    ),
    ("GSE274546", "GSE276841"): (
        separate_sources,
        "The seven CARE sites exclude Johns Hopkins, where Jackson collected spatial specimens GBM030 and GBM049.",
        "CARE Methods seven institutions; Jackson Methods and GSE276841 GEO Samples",
        "Inference based on collection sites.",
    ),
    ("GSE274546", "GeoMx_Zenodo16839828"): (
        separate_sources,
        "CARE's site list excludes Rigshospitalet, where both Artzi groups were treated.",
        "CARE Methods seven institutions; Artzi Methods Patient selection",
        "Inference based on collection sites.",
    ),
    ("GSE274546", "PDC000514"): (
        possible_overlap,
        "Both reports list SNU and Pitié-Salpêtrière. CARE reports 20 and 28 tumors, respectively.",
        "CARE Methods seven institutions; Kim Key resources table; CARE Supplementary Table 1 link",
        "A common institution supports possible overlap. Individual identifiers and sampling years were unavailable.",
    ),
    ("GSE274546", "GSE121810"): (
        possible_overlap,
        "CARE includes MD Anderson. Later NCT02852655 Methods list that site. "
        "The original trial enrolled in 2016–2017.",
        "CARE Methods seven institutions; Cloughesy PDF Results dates; same-NCT stage2 Methods sites",
        "The later site list does not locate the original RNA patients. No Pt-to-CARE mapping was available.",
    ),
    ("GSE274546", "GSE154795"): (
        separate_sources,
        "CARE's seven collecting institutions exclude UCLA, where Lee's 69 patients had surgery.",
        "CARE Methods seven institutions; Lee Results patient characteristics",
        "Inference based on collection sites.",
    ),
    ("GSE278456", "GSE276841"): (
        known_overlap,
        "The two Jackson datasets list GBM030 and GBM049 in both single-cell and spatial records.",
        "GSE278456 GEO Samples GSM8546885/GSM8546893; GSE276841 GEO Samples GSM8506649/GSM8506650",
        "At least two documented patients are shared across modalities.",
    ),
    ("GSE278456", "GeoMx_Zenodo16839828"): (
        separate_sources,
        "Jackson enrolled untreated primary cases at Johns Hopkins. Artzi enrolled at Rigshospitalet.",
        "Jackson Methods; GSE278456 GEO Summary; Artzi Methods",
        "Inference based on collection sites and enrollment definitions.",
    ),
    ("GSE278456", "PDC000514"): (
        separate_sources,
        "Jackson collected untreated specimens at Johns Hopkins. Kim's five sites were in Korea and Paris.",
        "Jackson Methods; Kim Key resources table",
        "Inference based on collection sites.",
    ),
    ("GSE278456", "GSE121810"): (
        separate_sources,
        "Jackson's cases at Johns Hopkins had no preoperative chemoradiation. The trial enrolled treated recurrences.",
        "Jackson Methods/GEO Summary; Cloughesy PDF Results/Methods; same-NCT stage2 Methods sites",
        "The later site list lacks Johns Hopkins and cannot resolve individual stage 1 cases. Treatment stages differ.",
    ),
    ("GSE278456", "GSE154795"): (
        separate_sources,
        "Jackson's untreated cases were at Johns Hopkins. Lee's recurrent 21-patient comparison was from UCLA.",
        "Jackson Methods and GSE278456 GEO Summary; Lee Results/Methods",
        "Inference based on collection sites and case definitions.",
    ),
    ("GSE276841", "GeoMx_Zenodo16839828"): (
        separate_sources,
        "Jackson's untreated GBM030/049 spatial specimens were from Johns Hopkins. Artzi used Rigshospitalet tissue.",
        "GSE276841 GEO Summary/Samples; Artzi Methods Patient selection",
        "Inference based on collection sites.",
    ),
    ("GSE276841", "PDC000514"): (
        separate_sources,
        "Jackson's two spatial patients came from Johns Hopkins, which is absent from Kim's five sites.",
        "GSE276841 GEO Summary; Kim Key resources table",
        "Inference based on collection sites.",
    ),
    ("GSE276841", "GSE121810"): (
        separate_sources,
        "Jackson's spatial samples were treatment-naive. The trial enrolled treated recurrences at other listed sites.",
        "GSE276841 GEO Summary; Cloughesy PDF Results/Methods; same-NCT stage2 Methods sites",
        "Inference based on sample collection sites and treatment stages.",
    ),
    ("GSE276841", "GSE154795"): (
        separate_sources,
        "Jackson's two spatial specimens were collected before treatment at Johns Hopkins. Lee had UCLA recurrences.",
        "GSE276841 GEO Summary; Lee Results/Methods",
        "Inference based on collection sites.",
    ),
    ("GeoMx_Zenodo16839828", "PDC000514"): (
        separate_sources,
        "Artzi's trial and controls came from Rigshospitalet. Kim listed hospitals in Korea and Paris.",
        "Artzi Methods Patient selection; Kim Key resources table",
        "Inference based on collection sites.",
    ),
    ("GeoMx_Zenodo16839828", "GSE121810"): (
        separate_sources,
        "Artzi used Copenhagen trial NCT03890952 with nivolumab and local controls. "
        "Cloughesy ran a US pembrolizumab trial.",
        "Artzi Methods Patient selection; Cloughesy PDF Results trial; same-NCT stage2 Methods sites",
        "Inference uses sites and regions in addition to the separate trial identifiers.",
    ),
    ("GeoMx_Zenodo16839828", "GSE154795"): (
        separate_sources,
        "Artzi's patients were at Rigshospitalet and received nivolumab. Lee enrolled at UCLA and used pembrolizumab.",
        "Artzi Methods Patient selection; Lee Results/Methods",
        "Inference based on collection sites.",
    ),
    ("PDC000514", "GSE121810"): (
        separate_sources,
        "Kim's sites were in Korea and France. Later NCT02852655 Methods describe seven US sites.",
        "Kim Key resources table; Cloughesy PDF Results; same-NCT stage2 Methods sites",
        "The later roster does not locate individual original RNA cases. Inference uses the collection regions.",
    ),
    ("PDC000514", "GSE154795"): (
        separate_sources,
        "Kim's five sites exclude UCLA, where all Lee patients had surgery.",
        "Kim Key resources table; Lee Results patient characteristics",
        "Inference based on collection sites.",
    ),
    ("GSE121810", "GSE154795"): (
        possible_overlap,
        "NCT02852655 enrolled in 2016–2017 across seven sites, including UCLA. "
        "Lee cites the trial and used off-label, off-trial pembrolizumab. "
        "Pt-to-LB mappings, individual surgery dates and prior enrollment were unavailable.",
        "Cloughesy PDF Results dates/Fig2 Pt A-B; Lee Methods Human samples, Results 13 vs8; GSE154795 GEO Samples; "
        "Lu 2021 follow-on article Results",
        "The 13 exposed patients in Lee received off-trial treatment. Zero crosswalk matches were confirmable "
        "and up to 13 could overlap. Independence of the 21-patient cohort remains unresolved.",
    ),
}
expected_pairs = list(combinations(source_table, 2))
if len(source_table) != 8 or set(pair_evidence_table) != set(expected_pairs):
    raise ValueError("The recorded evidence must cover every pair of the eight sources.")

# Record one row per source, then one row per unordered pair in the original accession order.
sources_table = [{"dataset": dataset, **source} for dataset, source in source_table.items()]
pairs_table = []
for pair_number, pair in enumerate(expected_pairs, start=1):
    source_a, source_b = pair
    classification, finding, locator, limitation = pair_evidence_table[pair]
    references = []
    if "GSE121810" in pair:
        references.extend(
            (
                "https://www.nature.com/articles/s41467-024-54326-7 [Methods > Study design and patients, "
                "NCT02852655 seven IRB sites; Results expansion dates 2020–2021 (not original GSE121810 sampling)]",
                "https://ucla.clinicaltrials.researcherprofiles.org/trial/NCT02852655 [Locations; Start Date]",
            )
        )
    if pair == ("GSE121810", "GSE154795"):
        references.append(
            "https://www.nature.com/articles/s41467-021-24293-4.pdf "
            "[Results p.2 follow-on off-label 13; separate from original trial]"
        )
    if pair == ("GSE174554", "GSE274546"):
        references.append(f"{source_table['GSE274546']['url']} [Results ref.14=Wang 2022; external validation]")
    inference = (
        "no" if classification == known_overlap else "yes (source or timing evidence; no patient identity match)"
    )
    pairs_table.append(
        dict(
            pair_number=pair_number,
            source_a=source_a,
            source_b=source_b,
            classification=classification,
            observed_evidence=finding,
            evidence_url_a=source_table[source_a]["url"],
            evidence_locator_a=source_table[source_a]["locator"],
            evidence_url_b=source_table[source_b]["url"],
            evidence_locator_b=source_table[source_b]["locator"],
            pair_locator=locator,
            additional_primary_evidence=" | ".join(references),
            inference=inference,
            limitation=limitation,
        )
    )

# These are the original access attempts. The script records the review and does not reopen URLs.
attempts = [
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6408961/", "Cloughesy full text", "reCAPTCHA",
        "Read the archived 35-page UCSF eScholarship PDF.",
    ),
    (
        "https://www.nature.com/articles/s41591-018-0337-7.pdf", "Cloughesy publisher PDF",
        "redirect to subscription HTML", "Read the UCSF archive PDF.",
    ),
    (
        "https://clinicaltrials.gov/study/NCT02852655", "trial registry details",
        "JavaScript shell, no study fields in web view",
        "Read the UCLA trial site and original and follow-on trial papers.",
    ),
    (
        "https://clinicaltrials.gov/api/v2/studies/NCT02852655", "trial registry API", "web tool Internal Error",
        "Read the UCLA trial site and original and follow-on trial papers.",
    ),
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10939876/", "Kim PDC paper full text", "reCAPTCHA",
        "Indexed Key resources and Figure 1A were available. ScienceDirect returned HTTP 403.",
    ),
    (
        "https://www.sciencedirect.com/science/article/pii/S1535610823004439", "Kim publisher full text", "HTTP 403",
        "Used indexed paper content. Table S1 was unavailable.",
    ),
    (
        "https://proteomic.datacommons.cancer.gov/pdc/study/PDC000514", "PDC metadata web page", "only JS shell",
        "Checked public clinical and biospecimen JSON. Individual sampling years were unavailable.",
    ),
    (
        "https://www.science.org/doi/10.1126/science.abm5214", "Jackson publisher full text", "HTTP 403",
        "Used indexed PMC Methods and GEO records.",
    ),
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC12836367/", "Jackson PMC full text", "reCAPTCHA",
        "Used indexed Methods and GEO records.",
    ),
    (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278456", "Jackson scRNA GEO direct", "reCAPTCHA",
        "Used the indexed GBM030 and GBM049 sample records.",
    ),
    (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276841", "Jackson spatial GEO direct", "reCAPTCHA",
        "Used the indexed GBM030 and GBM049 sample records.",
    ),
    (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554", "Wang GEO direct", "reCAPTCHA",
        "Read the Nature Cancer Methods.",
    ),
    (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546", "CARE GEO direct", "reCAPTCHA",
        "Read Nature Genetics Methods and indexed GSM records.",
    ),
]
attempt_fields = ("url", "object", "failure", "fallback")
access_table = []
for attempt in attempts:
    if len(attempt) != len(attempt_fields):
        raise ValueError("Each access attempt must have the four recorded fields.")
    access_table.append(dict(zip(attempt_fields, attempt)))

# Create the reports exclusively and preserve the original column order.
output_dir = Path(__file__).resolve().parents[1] / "step5_closeout_v2"
output_dir.mkdir(parents=True, exist_ok=True)
for filename, report_table in (
    ("overlap_source_evidence.csv", sources_table),
    ("overlap_28_pair_reaudit.csv", pairs_table),
    ("overlap_access_attempts.csv", access_table),
):
    output_path = output_dir / filename
    fieldnames = list(report_table[0])
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_table)
    if output_path.stat().st_size == 0:
        raise OSError(f"Empty provenance report {filename}")

classification_counts = Counter(row["classification"] for row in pairs_table)
print("Sources", len(sources_table), "pairs", len(pairs_table), "classifications", dict(classification_counts))
