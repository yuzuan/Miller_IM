# External data layout

Set `MILLER_IM_DATA_ROOT` to the directory containing downloaded public data. The scripts expect the following paths below that root. Source filenames are retained from the public deposits.

```text
$MILLER_IM_DATA_ROOT/
├── GSE174554/
│   └── GSE174554_RAW/
├── GSE274546/
│   └── <library_id>/{matrix.mtx.gz,genes.tsv.gz,barcodes.tsv.gz}
├── results/manifests/
│   ├── gse174554_libraries.csv
│   └── gse274546_standard_libraries.csv
└── write/
    ├── 26_pure_bioinformatics_dataset_rescue/source_metadata/
    │   ├── GeoMx_zenodo16839828/
    │   └── GSE154795_*
    ├── 28_supportive_dataset_search/source_metadata/GSE121810_*.xlsx
    ├── 30_dataset_rescue_search/source_metadata/
    │   ├── GSE278456_myeloid_h5/
    │   ├── GSE276841_spatial/
    │   └── NIHMS2115262-supplement-S2_S5_S6.xls
    └── 31_proteomic_validation/source_metadata/
        ├── PDC000514/
        ├── MSV000087947/
        └── Zenodo7646550/
```

The manifest-building scripts create the two library manifests from downloaded GEO files. No raw participant identifiers are used; the workflow retains only deidentified study identifiers supplied by the public deposits.

## Repository-resident eligibility metadata

The IDH eligibility table is versioned with the code and is not expected under `MILLER_IM_DATA_ROOT`:

```text
config/
├── miller_im_gene_set.tsv
└── GSE174554_formal_pair_IDH_crosswalk.csv
```

The crosswalk contains public deidentified sample identifiers, the IDH annotation reconciled from Wang *et al.* Supplementary Table 1, the source URL and workbook SHA-256, and the formal inclusion flag. Step 38 and Step 46 both read this same file from `config/`. The only excluded formal GSE174554 pair is `GSE174554_pair33` (`SF8963`/`SF12165`); the remaining 17 pairs are IDH-wild-type.

Generated pseudobulks remain in the repository-local `write/34_gse174554_raw_independent_discovery/03_myeloid/` and `write/36_gse274546_raw_independent_reannotation/03_myeloid/` paths. `MILLER_IM_DATA_ROOT` is used only for downloaded public inputs and large external data required by reconstruction or Figure 3; in particular, `53_fig3_visual_variety_rebuild.py` never requires a machine-specific absolute path.
