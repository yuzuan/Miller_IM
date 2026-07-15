# Analysis workflow

Run commands from the repository root. The scripts preserve their original numbered prefixes so that provenance can be traced to the frozen figure source tables.

| Stage | Scripts | Output used downstream |
|---|---|---|
| Public-data manifests | `build_gbm_manifests.py`, `30_prepare_gse278456_manifest.py` | Library/sample manifests |
| GSE174554 reconstruction | `34_prepare_sf8963_matrix.py`, `34_gse174554_raw_qc_scdblfinder.R`, `34_gse174554_raw_fullcell_myeloid.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| GSE274546 reconstruction | `convert_gse274546_rds_to_sparse.R`, `36_gse274546_raw_qc_scdblfinder.R`, `36_gse274546_raw_independent_reannotation.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| Paired recurrence tests | `38_independent_cohort_mg_inflammatory_recalculation.R` | Paired edgeR statistics, nominal single-set fgsea results, gene directions |
| Figure 1 | `39_fig1_independent_cohorts.py`, `41_fig1_independent_recurrence.py`, `42_fig1_fig2_visual_story.py`, `45_fig1_program_rebuild.py`, `46_fig1_lopo_gsea.R`, `46_fig1_ef_candidates.py`, `51_fig1_final_reorganized.py` | Six-panel recurrence figure and LOPO audit |
| Figure 2 | `30_gse278456_myeloid_spatial_validation.py`, `48_fig2_identity_state_rebuild.py`, `50_fig2_two_panel_final.py` | State localization and marker-leakage control |
| Figure 3 | `39_fig3_spatial.py`, `41_fig3_spatial.py`, `53_fig3_visual_variety_rebuild.py`, `55_fig3_visual_polish.py` | Paired GeoMx and linked Visium analyses |
| Figure 4 | `41_fig4_proteome.R`, `54_fig4_proteome_rebuild.py` | Paired protein score, protein effects, nominal competitive enrichment |
| Figure 5 | `41_fig5_pd1.R`, `56_fig5_pd1_rebuild.py` | Anti-PD-1 gene-set and adjusted patient-score analyses |
| Final checks | `57_main_figure_consistency_audit.py`, `validate_release.py` | Figure-value audit and frozen-release audit |

The main recurrence analysis is prespecified at at least 20 retained myeloid cells per patient and time point. The 50-cell analysis is a sensitivity analysis. GSE278456 and the two linked GSE276841 Visium cases originate from the same study and are not treated as independent patient replication.

