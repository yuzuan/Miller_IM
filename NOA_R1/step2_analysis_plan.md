# Step 2 analysis specification

The October revision protocol was locked on 2026-10-19 23:45:58 +0800 after the earlier five-set analysis. This summary records the settings used by the revision scripts.

Use the Step 1 paired voom fits for the same 17 and 45 patients and the same filtered gene backgrounds. For each fixed set report coverage and CAMERA with estimated within-set correlation and with correlation fixed at 0.01. Apply BH within each cohort and correlation mode for the fixed comparison family. Patient scores use gene-wise standardized TMM log-CPM and are summarized by paired differences.

Use the author-ranked 20-, 50-, and 100-gene Miller inflammatory definitions, the 20-gene Miller systemic inflammatory program, and the published Abdelfattah and Sankowski comparators in the locked table. The Hallmark workflow uses the official human inflammatory response set. Draw 10,000 size-matched and 10,000 size-and-expression-matched sets per cohort from its tested genes after excluding all 20 target genes. Use the locked cohort- and scheme-specific seeds and compare each random set's weighted GSEA enrichment score with the fixed target score.

In each anti-PD-1 cohort, reproduce CAMERA with the original 0.01 correlation on the original formal samples and design, then rerun with the estimated correlation.
