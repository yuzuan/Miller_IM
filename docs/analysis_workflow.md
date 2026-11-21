# Analysis workflow

Run commands from the repository root. Scripts are numbered consecutively by data flow. `script_names.tsv` maps the previous names to the current entry points; existing output directory names are retained for compatibility with the frozen tables.

| Stage | Scripts | Output used downstream |
|---|---|---|
| Public-data manifests | `01a_公共数据清单.py`, `01b_GSE278456样本清单.py` | Library/sample manifests |
| GSE174554 reconstruction | `02a_GSE174554补充矩阵.py`, `02b_GSE174554质控.R`, `02c_GSE174554髓系重建.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| GSE274546 reconstruction | `03a_GSE274546矩阵转换.R`, `03b_GSE274546质控.R`, `03c_GSE274546独立重注释.py` | Quality-filtered pan-myeloid raw-count pseudobulk |
| Paired recurrence tests | `04_配对复发统计.R` | IDH-restricted paired edgeR statistics, nominal single-set fgsea results, gene directions, and GSE174554 exclusion audit |
| Figure 1 | `05a_队列证据图.py`, `05b_配对复发图.py`, `05c_复发证据布局.py`, `05d_炎症程序图.py`, `05e_留一患者分析.R`, `05f_共享基因与稳健性图.py`, `05g_Figure1主图.py` | Six-panel recurrence figure and LOPO audit |
| Figure 2 | `06a_GSE278456髓系与空间验证.py`, `06b_细胞身份与状态图.py`, `06c_Figure2主图.py` | State localization and marker-leakage control |
| Figure 3 | `07a_空间证据整理.py`, `07b_空间关联图.py`, `07c_空间与配对表达图.py`, `07d_Figure3主图.py` | Paired GeoMx and linked Visium analyses |
| Figure 4 | `08a_配对蛋白统计.R`, `08b_Figure4主图.py` | Paired protein score, protein effects, nominal competitive enrichment |
| Figure 5 | `09a_PD1治疗相关分析.R`, `09b_Figure5主图.py` | Anti-PD-1 gene-set and adjusted patient-score analyses |
| Final checks | `10a_主图结果核验.py`, `10b_公开数据核验.py` | Figure-value audit and frozen-release audit |

After activating the recorded Python environment, `03d_GSE274546运行入口.sh` runs the GSE274546 QC and reannotation stages in sequence. Matrix conversion must already be complete.

The main recurrence analysis is prespecified at at least 20 retained myeloid cells per patient and time point. The 50-cell analysis is a sensitivity analysis. GSE278456 and the two linked GSE276841 Visium cases originate from the same study and are not treated as independent patient replication.

## IDH restriction and analysis order

`config/GSE174554_formal_pair_IDH_crosswalk.csv` is a sample-level reconciliation of the 18 GSE174554 pairs that otherwise satisfy the formal 20-cell endpoint threshold against Wang *et al.* Supplementary Table 1. It retains 17 IDH-wild-type pairs and excludes `GSE174554_pair33` (`SF8963`/`SF12165`), for which both specimens are annotated as IDH-mutant. Module 04 applies this eligibility rule before library normalization, paired edgeR fitting, score calculation, and fgsea. Script 05e reads the same crosswalk independently before constructing the full-data and leave-one-patient-out jobs; the excluded pair therefore cannot appear as a fold.

The formal recurrence sets are:

| Dataset | Retained pairs | Analyzed myeloid cells | IDH basis |
|---|---:|---:|---|
| GSE174554 | 17 | 12,489 | Sample-level crosswalk to Wang Supplementary Table 1 |
| GSE274546 | 45 | 55,568 | Original-study WHO 2021 IDH-wild-type diagnosis and bulk-DNA confirmation |

The resulting GSE174554 Miller-IM enrichment is NES 2.4290789842 with nominal *P* 4.001362 × 10⁻⁷. The Figure 1 leading-edge intersection contains seven genes and its gene-level effect panel contains 14 cohort-specific tests, with 0/14 at FDR < 0.05.
