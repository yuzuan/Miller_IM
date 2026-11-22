# NOA R1 workflow

Use the [rebuild commands](README.md#rebuild-figures-and-tables) to render Figures 1–5 from `inputs/`. The input tree retains the paths expected by the existing panel scripts. Figure outputs go to `figures_R1/`. For raw-data reconstruction, follow the separate [preparation workflow](../docs/analysis_workflow.md).

R analysis scripts take the revision directory as their argument. Python figure entry points resolve it from their script location. Set `NOA_LEGACY_ROOT` to the supplied `NOA_R1/inputs/` tree for frozen upstream tables. The [README](README.md#refit-the-statistical-analyses) shows how to use the supplied pseudobulk counts and metadata for the primary fit.

| Workflow | Entry points in `scripts/` | Derived results |
|---|---|---|
| Paired recurrence | `01a_primary_recurrence_camera.R` | Primary CAMERA, paired changes, and voom models for 02d, 02e, 04a and 05b. Source data for Figure 1, Figure S1 and Table S2. |
| Program comparators | `02d_extended_specificity.R`, then `02e_hallmark_comparators.R`. Use 02a and 02b to verify the fixed gene sets against the source workbooks. | Fixed and random comparator analyses for Figure S2 and Table S3. |
| PD-1 exposure CAMERA | `02f_pd1_camera.R` | CAMERA results and the tested samples and genes for 05c, 05d, Figure 5 and Table S8. |
| Myeloid states | `03a_aggregate_state_counts.py`, then `03b_state_expression.R` and `03c_state_abundance.py` | Paired state-expression and abundance analyses for Figure S3 and Table S4. The distinct GSE278456 inputs for main Figure 2 come from the upstream preparation workflow. |
| Steroids and GeoMx | `04a_steroid_and_geomx_models.R`, `04b_normalize_steroid_metadata.R`, `04c_steroid_center_summary.py`, `04d_complete_geomx_and_steroid_reports.R` | Paired steroid and GeoMx results for Tables S5–S6. Pass the GBM-CARE workbook to 04c with `--clinical-workbook`. |
| Protein, filtering, PD-1 models and overlap | `05a_fixed_protein_subset.R`, `05b_filter_consistency.R`, `05c_pd1_model_stability.R`, `05d_unify_pd1_ranks.R`, `05e_pd1_covariate_models.R`, `05f_overlap_reassessment.py` | Protein subsets, expression-filter sensitivity, PD-1 models and rank comparisons, and dataset overlap for Figures 4–5 and Tables S1, S2, S7 and S8. Run 05c before 05e. |
| Figure and table exports | `09a_supplementary_figures.py`, `09b_export_final_figures.py`, `09c_supplementary_tables.py` | Main and supplementary figures, and editable supplementary tables. Main figure layout and panel modules are in 06–08. |

The state-expression fit (03b) reads the archived count matrices, metadata and coverage tables already in `step3_state_analysis/`. To redo the aggregation (03a), use the frozen cell assignments and markers under `NOA_LEGACY_ROOT/write/71_recurrence_state_decomposition/` and the two clean-myeloid H5AD files produced by repository scripts 02c and 03c. Place those H5AD files and the paired pseudobulk files below `NOA_PUBLIC_DATA_ROOT` at the paths in the manifests.

The archived clinical source table `step4_steroids_geomx/clinical_case_normalized/steroid_missing_source_audit.csv` supplies the tissue-source and processing-laboratory records for Table S5.

The figure scripts use the palette and export helpers in `00a_functions.py`. The headers of the analysis and plotting scripts list the input files and output locations.
