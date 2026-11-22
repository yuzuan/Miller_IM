# Data preparation for NOA R1

The repository-level scripts generate the inputs used by NOA R1. Script numbers indicate execution order. Output directory names match the paths read by the revision scripts.

| Module | Role | Main inputs | Main consumers |
|---|---|---|---|
| 01 | Public-data manifests | Downloaded study metadata | Modules 02, 03 and 06 |
| 02 | GSE174554 reconstruction | Raw matrices and library manifest | Module 04 |
| 03 | GSE274546 reconstruction | Raw matrices and library manifest | Module 04 |
| 04 | Paired recurrence statistics | Patient pseudobulks, fixed genes and eligibility metadata | Module 05 and R1 recurrence analyses |
| 05 | Recurrence figure inputs | Recurrence tests and leave-one-patient-out results | Module 07 and R1 Figure 1 |
| 06 | State figure inputs | GSE278456 expression and state assignments | R1 Figure 2 |
| 07 | Paired spatial figure inputs | GeoMx, Visium and shared genes from Figure 1 | R1 GeoMx analysis and Figure 3 |
| 08 | Paired proteomic figure inputs | Public paired proteomic matrices | R1 protein analysis and Figure 4 |
| 09 | Anti-PD-1 figure inputs | Public anti-PD-1 expression and clinical metadata | R1 exposure analyses and Figure 5 |

Shared plotting functions are in `scripts/00a_functions.py`. They retain the original fonts, dimensions and palette, save PDF masters, and render PNGs through Poppler. Python loads this numbered module with `importlib.import_module`. Data reconstruction helpers are in `scripts/00a_functions_data.py`, separate from plotting dependencies. R scripts source `scripts/00a_functions.R` for shared input, enrichment and export functions.

Python runs the AnnData, Scanpy, SciPy and Matplotlib workflows. R runs the edgeR, limma, fgsea and scDblFinder workflows.

See the [execution sequence](docs/analysis_workflow.md) and [public data layout](docs/data_layout.md). The scripts write intermediate tables to `write/` and rendered outputs to `figures/`.

See the [NOA R1 workflow](NOA_R1/ARCHITECTURE.md) for the revision analyses and figure exports.
