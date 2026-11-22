# Miller-IM NOA R1 revision

Final English analysis scripts and derived source tables for the October 2026 first revision of the paired glioblastoma study.

## Contents

- `scripts/` contains the data processing, statistical analyses, and figure and table exports.
- `step1_*` through `step5_*` contain the final result tables and inputs used by downstream analyses and exports.
- `figures_R1/_font_build/fig*/data/Figure*/` contains the Figure 1–5 panel manifests and source tables.
- `figure3_scale_calibration/` contains the public Visium coordinates and the Figure 3C calibration results.
- [The input manifest](provenance/step1_input_manifest.csv) records the paired pseudobulk inputs.
- [The fixed gene sets](step2_locked_gene_sets.csv) define the target and comparator programs.

## Rebuild figures and tables

Install the Python and Poppler dependencies in [environment/](../environment/). Run the commands from the repository root in a working copy. The supplementary figures and tables read the supplied result tables.

```sh
python NOA_R1/scripts/09a_supplementary_figures.py
python NOA_R1/scripts/09c_supplementary_tables.py
```

For the five main figures, first prepare the upstream inputs listed in the [preparation workflow](../docs/analysis_workflow.md). Set `NOA_LEGACY_ROOT` if the prepared outputs are stored outside the repository root. Figure 3C uses the original 3354 × 2267 Step 55 tissue raster, its spot source table, and the [Visium calibration inputs](figure3_scale_calibration/README.md).

```sh
python NOA_R1/scripts/09b_export_final_figures.py
```

## Refit the statistical analyses

1. Install the R dependencies in [environment/](../environment/), along with `jsonlite` and `msigdbr` 25.1.1. The GO comparison uses the official MSigDB 2025.1 human cache. The Hallmark workflow can also use the supplied membership and source-version tables.
2. Download the inputs in the [dataset inventory](../config/datasets.tsv) and follow the [external data layout](../docs/data_layout.md). Set `NOA_PUBLIC_DATA_ROOT` to the public-data directory, with the prepared pseudobulk counts and metadata at the paths in [the input manifest](provenance/step1_input_manifest.csv).
3. Follow the [analysis order and input requirements](ARCHITECTURE.md) in a working copy. Move existing outputs aside before running scripts that require a fresh output location. Use the supplied fixed gene sets for ordinary runs. To verify them against the source workbooks, use 02a and 02b with the workbooks named in the [fixed gene-set table](step2_locked_gene_sets.csv) and set `NOA_ABDEL_WORKBOOK` to Abdelfattah Supplementary Data 9.

```sh
Rscript NOA_R1/scripts/01a_primary_recurrence_camera.R NOA_R1
```

Raw data are available from GEO (GSE174554, GSE274546, GSE278456, GSE276841, GSE121810, GSE154795), [Zenodo](https://doi.org/10.5281/zenodo.16839828), and PDC (PDC000514). The gene-set supplements are from [Miller et al.](https://doi.org/10.1038/s41586-025-08633-8), [Abdelfattah et al.](https://doi.org/10.1038/s41467-022-28372-y), and [Sankowski et al.](https://doi.org/10.1038/s41593-019-0532-y). Cite the study and the source datasets when using this material.
