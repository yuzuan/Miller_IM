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

The manifest-building scripts create the two library manifests from downloaded GEO files and the supplied sample annotation tables. The workflow retains the deidentified study identifiers supplied by the public deposits.

## Paired-cohort library manifests

Use the [dataset inventory](../config/datasets.tsv) to obtain the two GEO series. Extract the GSE174554 RAW archive into `$MILLER_IM_DATA_ROOT/GSE174554/GSE174554_RAW/` and place its `GSE174554_Tumor_normal_metadata.txt.gz` supplement in `$MILLER_IM_DATA_ROOT/GSE174554/`. Place the GSE274546 matrices under `$MILLER_IM_DATA_ROOT/GSE274546/<library_id>/`. Convert each source RDS with `scripts/03a_gse274546_matrix_conversion.R`. The supplied GSE274546 sample annotation lists the patient and time point used to name each library `P<patient>T<timepoint>`.

Save the [GSE174554 family SOFT file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174554/soft/GSE174554_family.soft.gz) in the GSE174554 directory. Copy the supplied sample annotations, then split the SOFT sample records into the local cache and build the manifests. Run from the repository root.

```sh
mkdir -p "$MILLER_IM_DATA_ROOT/tmp_geo_gsm"
cp NOA_R1/inputs/GSE174554/GBM_GSM_samples.csv "$MILLER_IM_DATA_ROOT/GSE174554/"
cp NOA_R1/inputs/GSE274546/GSE274546_GPL24676_sample_anno.csv "$MILLER_IM_DATA_ROOT/GSE274546/"

# Cache the public sample fields under their GEO accession for the manifest builder.
gzip -dc "$MILLER_IM_DATA_ROOT/GSE174554/GSE174554_family.soft.gz" |
    awk -v cache="$MILLER_IM_DATA_ROOT/tmp_geo_gsm" '
        /^\^SAMPLE = GSM[0-9]+$/ {
            if (path != "") close(path)
            path = cache "/" $3 ".txt"
        }

        /^!Sample_/ && path != "" {
            print > path
        }
    '

python scripts/01a_public_data_manifests.py \
    --project-root "$MILLER_IM_DATA_ROOT" \
    --geo-cache-dir "$MILLER_IM_DATA_ROOT/tmp_geo_gsm" \
    --gse274546-cache-dir "$MILLER_IM_DATA_ROOT/GSE274546" \
    --output-dir "$MILLER_IM_DATA_ROOT/results/manifests"
```

## Repository-resident eligibility metadata

The IDH eligibility table is versioned with the code and is not expected under `MILLER_IM_DATA_ROOT`:

```text
config/
├── miller_im_gene_set.tsv
└── GSE174554_formal_pair_IDH_crosswalk.csv
```

The crosswalk contains public deidentified sample identifiers, the IDH annotation reconciled from Wang *et al.* Supplementary Table 1, the source URL and workbook SHA-256, and the formal inclusion flag. The paired recurrence tests (module 04) and leave-one-patient-out analysis (script 05e) both read this file from `config/`. The only excluded formal GSE174554 pair is `GSE174554_pair33` (`SF8963`/`SF12165`); the remaining 17 pairs are IDH-wild-type.

Generated pseudobulks remain in the repository-local `write/34_gse174554_raw_independent_discovery/03_myeloid/` and `write/36_gse274546_raw_independent_reannotation/03_myeloid/` paths. `MILLER_IM_DATA_ROOT` is used only for downloaded public inputs and large external data required by reconstruction or Figure 3; in particular, `07c_spatial_and_paired_expression.py` never requires a machine-specific absolute path.

For refits from the archived pseudobulks use the [NOA R1 instructions](../NOA_R1/README.md#refit-the-statistical-analyses). To reaggregate states from cells, use the clean-myeloid H5AD files from the same two `03_myeloid/` directories. The [state input manifest](../NOA_R1/inputs/write/71_recurrence_state_decomposition/input_manifest.csv) lists the expected filenames under `NOA_PUBLIC_DATA_ROOT`.
