# Prepare inputs for NOA R1

Use this workflow to prepare the source data from public deposits. To render the five main figures from the supplied inputs, follow the [NOA R1 figure build](../NOA_R1/README.md#rebuild-figures-and-tables). Run preparation scripts from the repository root in script-number order.

| Stage | Scripts | Output used downstream |
|---|---|---|
| Public-data manifests | `01a_public_data_manifests.py`, `01b_gse278456_sample_manifest.py` | Library/sample manifests |
| GSE174554 reconstruction | `02a_gse174554_matrix_repair.py`, `02b_gse174554_quality_control.R`, `02c_gse174554_myeloid_reconstruction.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| GSE274546 reconstruction | `03a_gse274546_matrix_conversion.R`, `03b_gse274546_quality_control.R`, `03c_gse274546_independent_reannotation.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| Paired recurrence tests | `04_paired_recurrence_statistics.R` | IDH-restricted paired edgeR statistics, nominal single-set fgsea results, gene directions, and GSE174554 exclusion audit |
| Figure 1 | `05a_cohort_evidence.py`, `05b_paired_recurrence_panels.py`, `05c_recurrence_evidence_layout.py`, `05d_inflammatory_program_panels.py`, `05e_leave_one_patient_out.R`, `05f_shared_genes_and_robustness.py`, `05g_figure1.py` | Six-panel recurrence figure and LOPO audit |
| Figure 2 | `06a_gse278456_myeloid_spatial_validation.py`, `06b_cell_identity_and_state.py`, `06c_figure2.py` | State localization and marker-leakage control |
| Figure 3 | `07a_spatial_evidence.py`, `07b_spatial_associations.py`, `07d_figure3.py`, which imports `07c_spatial_and_paired_expression.py` | Paired GeoMx and linked Visium analyses |
| Figure 4 | `08a_paired_proteome_statistics.R`, `08b_figure4.py` | Paired protein score, protein effects, nominal competitive enrichment |
| Figure 5 | `09a_pd1_association_analysis.R`, `09b_figure5.py` | Anti-PD-1 gene-set and adjusted patient-score analyses |

After activating the recorded Python environment, `03d_run_gse274546.sh` runs the GSE274546 QC and reannotation stages in sequence. Matrix conversion must already be complete.

The main recurrence analysis is prespecified at at least 20 retained myeloid cells per patient and time point. The 50-cell analysis is a sensitivity analysis. GSE278456 and the two linked GSE276841 Visium cases originate from the same study and are not treated as independent patient replication.

## IDH restriction and analysis order

`config/GSE174554_formal_pair_IDH_crosswalk.csv` is a sample-level reconciliation of the 18 GSE174554 pairs that otherwise satisfy the formal 20-cell endpoint threshold against Wang *et al.* Supplementary Table 1. It retains 17 IDH-wild-type pairs and excludes `GSE174554_pair33` (`SF8963`/`SF12165`), for which both specimens are annotated as IDH-mutant. Module 04 applies this eligibility rule before library normalization, paired edgeR fitting, score calculation, and fgsea. Script 05e reads the same crosswalk independently before constructing the full-data and leave-one-patient-out jobs; the excluded pair therefore cannot appear as a fold.

The formal recurrence sets are:

| Dataset | Retained pairs | Analyzed myeloid cells | IDH basis |
|---|---:|---:|---|
| GSE174554 | 17 | 12,489 | Sample-level crosswalk to Wang Supplementary Table 1 |
| GSE274546 | 45 | 55,568 | Original-study WHO 2021 IDH-wild-type diagnosis and bulk-DNA confirmation |

Continue with the [NOA R1 workflow](../NOA_R1/ARCHITECTURE.md) after preparing the inputs for the analyses you plan to run.
