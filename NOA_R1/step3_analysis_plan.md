# Step 3 analysis specification

The state analysis settings were fixed on 2026-10-19. Preserve the nine existing marker-defined myeloid states, the two broader state groups, and the 17 and 45 formal patient pairs.

Aggregate integer counts by patient, time point and state. For expression, require at least six complete pairs, each with at least 20 cells per time point in the primary analysis and 10 in the sensitivity analysis. Use condition-based filtering, TMM and the paired edgeR model, then apply CAMERA on the corresponding voom fit when at least 16 target genes remain. Apply BY to the 11 state hypotheses per cohort, threshold and correlation mode, with BH also recorded.

For abundance, divide each state count by that sample's total number of analyzed myeloid cells. Compare within-patient changes using two-sided paired sign flips, exhaustive for 17 pairs and 100,000 draws for 45 pairs. Keep the coarse-group contrast separate from the nine-state correction family.
