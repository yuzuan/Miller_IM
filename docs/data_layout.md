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

