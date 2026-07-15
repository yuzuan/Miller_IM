# Miller-IM recurrence analysis

Code and figure source data for **An inflammatory myeloid program is enriched at recurrence in paired IDH-wildtype glioblastomas**.

The study tests a fixed 20-gene representation of the inflammatory microglial program reported by Miller *et al.* in paired primary and recurrent glioblastoma datasets. The patient is the inferential unit. Single-set enrichment results are reported as nominal *P* values; FDR is used only for explicitly defined multi-hypothesis families.

## Repository contents

| Path | Contents |
|---|---|
| `config/` | Fixed Miller-IM gene set and public-dataset inventory |
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

This verifies the fixed gene-set definition, reads every bundled CSV/CSV.GZ file, checks headline values used in Figures 1–5, scans for local absolute paths and misleading program-level multiplicity labels, and verifies the SHA-256 manifest.

## Re-running from public data

1. Download the source datasets listed in `config/datasets.tsv`.
2. Arrange them as described in `docs/data_layout.md`.
3. Set a single external-data location:

```bash
export MILLER_IM_DATA_ROOT=/absolute/path/to/miller_im_public_data
```

4. Run the cohort reconstruction and downstream scripts in the order shown in `docs/analysis_workflow.md`.

Scripts write generated intermediates to `write/` and figures to `figures/`; both directories are intentionally excluded from version control. Random seeds are fixed in the analysis scripts.

## Main statistical conventions

- Paired recurrence analyses use patient-level raw-count pseudobulks and an edgeR design of `~ patient + condition`.
- Miller-IM enrichment uses the fixed gene set in `config/miller_im_gene_set.tsv` and a ranking statistic of `sign(logFC) * sqrt(F)` for the principal recurrence analysis.
- Single-set fgsea and camera results are nominal *P* values; no FDR is claimed for these single-hypothesis tests.
- Cells, regions, spatial spots, and proteins are not treated as independent patient replicates.

## Data availability

Transcriptomic datasets are available through GEO under GSE174554, GSE274546, GSE278456, GSE276841, GSE121810, and GSE154795. Paired GeoMx data are available from Zenodo ([10.5281/zenodo.16839828](https://doi.org/10.5281/zenodo.16839828)); proteomic data are available through the Proteomic Data Commons under PDC000514.

## Citation

Please cite the associated manuscript and the original source studies. Repository citation metadata are provided in `CITATION.cff`.

No open-source license has been selected for this repository; absent a license, reuse beyond viewing and citation requires permission from the copyright holders.
