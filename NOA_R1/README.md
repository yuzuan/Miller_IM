# Miller-IM NOA R1 revision

Analysis scripts, derived data and figure inputs for the October 2026 first revision of the paired glioblastoma study.

## Contents

- `scripts/` contains the data processing, statistical analyses, and figure and table exports.
- `inputs/` contains the frozen inputs for Figures 1–5 and the prepared paired recurrence inputs.
- `step1_*` through `step5_*` contain the final result tables and inputs used by downstream analyses and exports.
- `figures_R1/_font_build/fig*/data/Figure*/` contains the Figure 1–5 panel manifests and source tables.
- `figure3_scale_calibration/` contains the public Visium coordinates and the Figure 3C calibration results.
- [The input manifest](provenance/step1_input_manifest.csv) records the paired pseudobulk inputs.
- [The fixed gene sets](step2_locked_gene_sets.csv) define the target and comparator programs.

## Rebuild figures and tables

Download and extract the repository ZIP, or clone the repository. Run the commands below from the repository root in a working copy. Use Python 3.10 on macOS or Linux, with Poppler's `pdftoppm` on PATH and Arial installed for the figure text. The [environment files](../environment/) record the dependencies.

```sh
python3.10 -m venv .venv
. .venv/bin/activate
python -m pip --isolated install --no-cache-dir --index-url https://pypi.org/simple -r environment/python-requirements.txt
```

The following commands use the supplied inputs to render all five main figures, the supplementary figures, and the supplementary tables. The main figure command rasterizes the panel PDFs with Poppler and writes 600-dpi PNG and TIFF files under `NOA_R1/figures_R1/`. Add `--figures 3` to run only Figure 3.

```sh
NOA_LEGACY_ROOT="$PWD/NOA_R1/inputs" python NOA_R1/scripts/09b_export_final_figures.py
python NOA_R1/scripts/09a_supplementary_figures.py
python NOA_R1/scripts/09c_supplementary_tables.py
```

Figure 3C uses the [archived 3354 × 2267 tissue map](inputs/figures/55_figure3_visual_polish/Figure3/Figure3C_spatial_score_atlas.png) as its background. The export uses the matching spot source table and [Visium registration files](figure3_scale_calibration/README.md) to place the labels and 1-mm scale bars.

## Refit the statistical analyses

Install the R dependencies in [environment/](../environment/) and `jsonlite`. The primary recurrence fit uses the prepared patient pseudobulk counts and metadata listed in the [input manifest](provenance/step1_input_manifest.csv). The state-expression fit uses the archived count matrices and matching metadata under `step3_state_analysis/pseudobulk/`. Run both in a working copy.

```sh
export NOA_LEGACY_ROOT="$PWD/NOA_R1/inputs"
export NOA_PUBLIC_DATA_ROOT="$PWD/NOA_R1/inputs"
Rscript NOA_R1/scripts/01a_primary_recurrence_camera.R NOA_R1
Rscript NOA_R1/scripts/03b_state_expression.R NOA_R1
```

For further refits and raw-data reconstruction, see the [analysis order](ARCHITECTURE.md), [dataset inventory](../config/datasets.tsv), and [external data layout](../docs/data_layout.md). Put any required downloads in the separate public-data directory or add them to the input tree in your working copy. If you use a separate public-data directory, place the prepared counts and metadata at the same manifest paths and point `NOA_PUBLIC_DATA_ROOT` at that directory.

Use the supplied [fixed gene sets](step2_locked_gene_sets.csv) for ordinary runs. To check them against the original workbooks, run 02a and 02b with the workbooks named in that table and set `NOA_ABDEL_WORKBOOK` to Abdelfattah Supplementary Data 9. The GO comparison uses `msigdbr` 25.1.1 and the official MSigDB 2025.1 human cache. The Hallmark workflow can use the supplied membership and source-version tables.

Raw data are available from GEO (GSE174554, GSE274546, GSE278456, GSE276841, GSE121810, GSE154795), [Zenodo](https://doi.org/10.5281/zenodo.16839828), and PDC (PDC000514). Cite the study and the source datasets when using this material.
