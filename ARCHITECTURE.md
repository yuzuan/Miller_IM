# Repository architecture

The patient-level analysis and the public figure tables use the same scientific definitions as the previous release. Script numbers indicate execution order; output directories retain their published names for compatibility.

| Module | Role | Main inputs | Main consumers |
|---|---|---|---|
| 01 | Public-data manifests | Downloaded study metadata | Modules 02, 03 and 06 |
| 02 | GSE174554 reconstruction | Raw matrices and library manifest | Module 04 |
| 03 | GSE274546 reconstruction | Raw matrices and library manifest | Module 04 |
| 04 | Paired recurrence statistics | Patient pseudobulks, fixed genes and eligibility metadata | Module 05 and the comparison panels in 09 |
| 05 | Figure 1 | Recurrence tests and leave-one-patient-out results | Module 07 and final audit |
| 06 | Figure 2 | GSE278456 expression and state assignments | Final audit |
| 07 | Figure 3 | GeoMx, Visium and shared genes from Figure 1 | Final audit |
| 08 | Figure 4 | Public paired proteomic matrices | Final audit |
| 09 | Figure 5 | Public anti-PD-1 expression and clinical metadata | Final audit |
| 10 | Validation | Generated figures or frozen public source tables | Repository users |

Shared plotting functions are in `scripts/00a_functions.py`. They retain the original fonts, dimensions and palette, save PDF masters, and render PNGs through Poppler. Python loads this numbered module with `importlib.import_module`. Data reconstruction helpers are in `scripts/00a_functions_data.py`, separate from plotting dependencies. R scripts source `scripts/00a_functions.R` for shared input, enrichment and export functions.

Python remains the implementation for the existing AnnData, Scanpy, SciPy and Matplotlib workflow. R remains the implementation for edgeR, limma, fgsea and scDblFinder. The style update does not migrate statistical methods between languages or change seeds, cohort rules, data fields, or frozen source tables.

The script-name map is maintained in `docs/script_names.tsv`. The execution sequence is documented in `docs/analysis_workflow.md`; downloaded inputs are described in `docs/data_layout.md`. Generated `write/` and `figures/` files stay outside version control.

The October R1 workflow, its figure source data, and its dependency order are documented in [NOA_R1](NOA_R1/ARCHITECTURE.md).
