# Miller-IM NOA R1 revision

Final English analysis scripts and derived source tables for the October 2026 first revision of the paired glioblastoma study.

## Contents

- `scripts/` contains the analysis, figure, supplementary table, and verification workflows.
- `step1_*` through `step5_*` contain the archived revision results and supplementary figure data.
- `figures_R1/_font_build/fig*/data/Figure*/` contains the Figure 1–5 panel manifests and source tables.
- `figure3_scale_calibration/` contains the public Visium coordinates and the Figure 3C calibration results.
- `provenance/` records the paired pseudobulk inputs and the [script names in the October archive](provenance/script_names.tsv).

## Validate the archived tables

From the repository root, run:

```sh
python NOA_R1/scripts/11_validate_revision_package.py
```

The validator uses the Python 3.10 standard library. It checks the cohort counts in the primary CAMERA tables, the five figure source directories, and the readable release files.

## Reproduce the analyses

1. Install the baseline Python, R, and Poppler dependencies in [environment/](../environment/). The revision additionally uses the R packages `jsonlite`, `msigdbr` 25.1.1, and `dqrng`. Hallmark and GO analyses use the official MSigDB 2025.1 human data cache.
2. Obtain the public data and run the preparation scripts. Set `NOA_LEGACY_ROOT` if the prepared outputs are stored outside the repository root. [Dataset inventory](../config/datasets.tsv) · [Data layout](../docs/data_layout.md) · [Preparation workflow](../docs/analysis_workflow.md)
3. Set `NOA_PUBLIC_DATA_ROOT` to the external public-data directory. Make the two prepared pseudobulk directories available under that root at the paths in [the input manifest](provenance/step1_input_manifest.csv). Set `NOA_ABDEL_WORKBOOK` to the Abdelfattah Supplementary Data 9 workbook when rebuilding the gene-set definitions. Place the Miller and Sankowski workbooks under `NOA_R1/sources/` using the names in `step2_locked_gene_sets.csv`.
4. Follow [the revision workflow](ARCHITECTURE.md) in a working copy. Analysis scripts that create a fresh output directory require its name to be free in that copy. The archived tables can be used as reference or as inputs to downstream scripts.

```sh
Rscript NOA_R1/scripts/01a_primary_recurrence_camera.R NOA_R1
python NOA_R1/scripts/09b_export_final_figures.py --figures 1 2 3 4 5 --skip-comparisons
```

Figure 3C reads the original pipeline's Step 55 tissue raster and spot source table. The supplied calibration files reproduce the registered 1-mm scale bars.

Raw data are available from GEO (GSE174554, GSE274546, GSE278456, GSE276841, GSE121810, GSE154795), [Zenodo](https://doi.org/10.5281/zenodo.16839828), and PDC (PDC000514). The gene-set supplements are from [Miller et al.](https://doi.org/10.1038/s41586-025-08633-8), [Abdelfattah et al.](https://doi.org/10.1038/s41467-022-28372-y), and [Sankowski et al.](https://doi.org/10.1038/s41593-019-0532-y). Cite the study and the source datasets when using this material.
