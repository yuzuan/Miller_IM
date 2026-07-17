# Analysis workflow

Run commands from the repository root. The scripts preserve their original numbered prefixes so that provenance can be traced to the frozen figure source tables.

| Stage | Scripts | Output used downstream |
|---|---|---|
| Public-data manifests | `build_gbm_manifests.py`, `30_prepare_gse278456_manifest.py` | Library/sample manifests |
| GSE174554 reconstruction | `34_prepare_sf8963_matrix.py`, `34_gse174554_raw_qc_scdblfinder.R`, `34_gse174554_raw_fullcell_myeloid.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| GSE274546 reconstruction | `convert_gse274546_rds_to_sparse.R`, `36_gse274546_raw_qc_scdblfinder.R`, `36_gse274546_raw_independent_reannotation.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| Paired recurrence tests | `38_independent_cohort_mg_inflammatory_recalculation.R` | IDH-restricted paired edgeR statistics, nominal single-set fgsea results, gene directions, and GSE174554 exclusion audit |
| Figure 1 | `39_fig1_independent_cohorts.py`, `41_fig1_independent_recurrence.py`, `42_fig1_fig2_visual_story.py`, `45_fig1_program_rebuild.py`, `46_fig1_lopo_gsea.R`, `46_fig1_ef_candidates.py`, `51_fig1_final_reorganized.py` | Six-panel recurrence figure and LOPO audit |
| Figure 2 | `30_gse278456_myeloid_spatial_validation.py`, `48_fig2_identity_state_rebuild.py`, `50_fig2_two_panel_final.py` | State localization and marker-leakage control |
| Figure 3 | `39_fig3_spatial.py`, `41_fig3_spatial.py`, `53_fig3_visual_variety_rebuild.py`, `55_fig3_visual_polish.py` | Paired GeoMx and linked Visium analyses |
| Figure 4 | `41_fig4_proteome.R`, `54_fig4_proteome_rebuild.py` | Paired protein score, protein effects, nominal competitive enrichment |
| Figure 5 | `41_fig5_pd1.R`, `56_fig5_pd1_rebuild.py` | Anti-PD-1 gene-set and adjusted patient-score analyses |
| Final checks | `57_main_figure_consistency_audit.py`, `validate_release.py` | Figure-value audit and frozen-release audit |

The main recurrence analysis is prespecified at at least 20 retained myeloid cells per patient and time point. The 50-cell analysis is a sensitivity analysis. GSE278456 and the two linked GSE276841 Visium cases originate from the same study and are not treated as independent patient replication.

## IDH restriction and analysis order

`config/GSE174554_formal_pair_IDH_crosswalk.csv` is a sample-level reconciliation of the 18 GSE174554 pairs that otherwise satisfy the formal 20-cell endpoint threshold against Wang *et al.* Supplementary Table 1. It retains 17 IDH-wild-type pairs and excludes `GSE174554_pair33` (`SF8963`/`SF12165`), for which both specimens are annotated as IDH-mutant. Step 38 applies this eligibility rule before library normalization, paired edgeR fitting, score calculation, and fgsea. Step 46 reads the same crosswalk independently before constructing the full-data and leave-one-patient-out jobs; the excluded pair therefore cannot appear as a fold.

The formal recurrence sets are:

| Dataset | Retained pairs | Analyzed myeloid cells | IDH basis |
|---|---:|---:|---|
| GSE174554 | 17 | 12,489 | Sample-level crosswalk to Wang Supplementary Table 1 |
| GSE274546 | 45 | 55,568 | Original-study WHO 2021 IDH-wild-type diagnosis and bulk-DNA confirmation |

The resulting GSE174554 Miller-IM enrichment is NES 2.4290789842 with nominal *P* 4.001362 × 10⁻⁷. The Figure 1 leading-edge intersection contains seven genes and its gene-level effect panel contains 14 cohort-specific tests, with 0/14 at FDR < 0.05.
