# Miller-IM recurrence analysis

Code and figure source data for **An inflammatory myeloid program is enriched at recurrence in paired IDH-wildtype glioblastomas**.

The study tests a fixed 20-gene representation of the inflammatory microglial program reported by Miller *et al.* in paired primary and recurrent IDH-wild-type glioblastoma datasets. The patient is the inferential unit. Single-set enrichment results are reported as nominal *P* values; FDR is used only for explicitly defined multi-hypothesis families.

## Repository contents

| Path | Contents |
|---|---|
| `config/` | Fixed Miller-IM gene set, public-dataset inventory, and audited GSE174554 IDH crosswalk |
| `scripts/` | Raw-data reconstruction, statistical analysis, figure generation, and validation |
| `source_data/Figure1`–`Figure5` | Frozen source tables underlying the current main figures |
| `environment/` | Recorded Python and R package versions |
| `manifests/` | SHA-256 manifest for the frozen source tables |
| `docs/` | Data layout and analysis workflow |

Raw public datasets and large intermediate objects are not redistributed. The bundled figure source tables contain only measurements derived from public, deidentified datasets and public study identifiers.

## Quick validation

Create a Python environment using the recorded versions, then run:

```bash
python scripts/validate_release.py
```

This verifies the fixed gene-set definition and the 18-row GSE174554 IDH crosswalk, reads every bundled CSV/CSV.GZ file, checks headline values used in Figures 1–5, scans text and tabular character fields for local absolute paths and misleading program-level multiplicity labels, and verifies the SHA-256 manifest.

## Re-running from public data

1. Download the source datasets listed in `config/datasets.tsv`.
2. Arrange them as described in `docs/data_layout.md`.
3. Set a single external-data location:

```bash
export MILLER_IM_DATA_ROOT=/absolute/path/to/miller_im_public_data
```

4. Run the cohort reconstruction and downstream scripts in the order shown in `docs/analysis_workflow.md`.

Scripts write generated intermediates to `write/` and figures to `figures/`; both directories are intentionally excluded from version control. Random seeds are fixed in the analysis scripts.

## IDH eligibility

The formal recurrence analysis is restricted to IDH-wild-type glioblastoma. GSE174554 sample identifiers were reconciled to Wang *et al.* Supplementary Table 1 in `config/GSE174554_formal_pair_IDH_crosswalk.csv`. The crosswalk contains the 18 pairs that otherwise passed the 20-cell endpoint threshold; `GSE174554_pair33` (`SF8963`/`SF12165`) is annotated as IDH-mutant at both time points and is excluded before normalization, differential expression, enrichment, and leave-one-patient-out analysis. The retained analysis therefore contains 17 GSE174554 pairs and 12,489 analyzed myeloid cells. The 45 analyzed GSE274546 pairs are IDH-wild-type by the original study's WHO 2021 diagnosis and bulk-DNA confirmation.

After this restriction, the GSE174554 Miller-IM result is NES 2.4290789842 with nominal *P* 4.001362 × 10⁻⁷. The two cohort-specific leading edges share seven genes, none of the 14 gene-by-cohort tests reaches FDR < 0.05, and the GSE174554 leave-one-patient-out analysis contains 17 folds without `GSE174554_pair33`.

## Main statistical conventions

- IDH eligibility is applied before formal GSE174554 recurrence testing; the exclusion is not performed as a post hoc figure filter.
- Paired recurrence analyses use patient-level raw-count pseudobulks and an edgeR design of `~ patient + condition`.
- Miller-IM enrichment uses the fixed gene set in `config/miller_im_gene_set.tsv` and a ranking statistic of `sign(logFC) * sqrt(F)` for the principal recurrence analysis.
- Single-set fgsea and camera results are nominal *P* values; no FDR is claimed for these single-hypothesis tests.
- Cells, regions, spatial spots, and proteins are not treated as independent patient replicates.

## Data availability

Transcriptomic datasets are available through GEO under GSE174554, GSE274546, GSE278456, GSE276841, GSE121810, and GSE154795. Paired GeoMx data are available from Zenodo ([10.5281/zenodo.16839828](https://doi.org/10.5281/zenodo.16839828)); proteomic data are available through the Proteomic Data Commons under PDC000514.

## Citation

Please cite the associated manuscript and the original source studies. Repository citation metadata are provided in `CITATION.cff`.

No open-source license has been selected for this repository; absent a license, reuse beyond viewing and citation requires permission from the copyright holders.
