# NOA R1 workflow

The revision reads prepared public-data inputs from the original pipeline and writes its outputs within the revision directory. R analysis scripts take the revision directory as their argument. Python figure entry points resolve it from their script location. Run the commands from the repository root.

| Workflow | Entry points in `scripts/` | Derived results |
|---|---|---|
| Paired recurrence | `01a_primary_recurrence_camera.R`, then `01b_patient_changes_figure.py` | Figure 1, Table S2 and Figure S1. The paired voom models feed steps 02 and 04. |
| Comparator and exposure analyses | `02a_prepare_gene_sets.py`, `02b_prepare_prespecified_sets.py`, `02c_gene_set_specificity.R`, `02d_extended_specificity.R`, `02e_hallmark_comparators.R`, `02f_pd1_camera.R`, `02g_program_specificity_figure.py` | Fixed and random comparator sets, anti-PD-1 CAMERA, Table S3 and Figure S2. The bundled locked definitions can be used directly. |
| Myeloid states | `03a_aggregate_state_counts.py`, `03b_state_expression.R`, `03c_state_abundance.py`, `03d_state_figures.py` | Figure 2, Table S4 and Figure S3. |
| Steroids and GeoMx | `04a_steroid_and_geomx_models.R`, `04b_normalize_steroid_metadata.R`, `04c_steroid_center_summary.py`, `04d_complete_geomx_and_steroid_reports.R` | Figure 3 and Tables S5–S6. |
| Protein, filtering, exposure models and provenance | `05a_fixed_protein_subset.R` through `05g_check_overlap_reassessment.py` | Figures 4–5 and Tables S1, S2, S7 and S8. |
| Figure and table exports | `09a_supplementary_figures.py`, `09b_export_final_figures.py`, `09c_supplementary_tables.py` | Main and supplementary figures, and editable supplementary tables. Main figure layout and panel modules are in 06–08. |

Figure 3C reads the checked-in Visium calibration inputs. The figure scripts use the palette and export helpers in `00a_functions.py`. The published panel layouts and source tables retain their original names and column definitions.

`02c_gene_set_specificity.R` reproduces the earlier five-set analysis used by the comparison figure and its independent review. Run it after `01a` and before `02g` or `02h`. The `10*` scripts preserve earlier analysis variants. The `02h_*`, `04e_*`, `04f_*`, and `05g_*`–`05j_*` scripts verify specific outputs or source provenance. Some historical audits use the original archived intermediate files or source materials. See the inputs in the respective script headers.

Use `05i_check_protein_and_pd1_ranks.R` to verify the 105 complete proteomic pairs. In the historical `05j_check_filtering_and_protein_models.R` workflow, the protein section records an error when the original patient list includes unmatched T1 samples.
