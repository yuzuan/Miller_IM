# Step 5 analysis specification

The sensitivity analyses were locked on 2026-10-20 00:55:04 +0800. Fit the specified unadjusted, myeloid-proxy-adjusted and full anti-PD-1 score models on each formal cohort's complete analysis set. Record coefficients, confidence intervals, influence diagnostics and leave-one-patient-out results. Repeat the recurrence CAMERA analysis with condition-based gene filtering and retain the original filter as the reference.

Use the original 105 PDC000514 patient pairs and the unique-peptide measurements for the principal fixed-protein analysis. Use the common measured subset when at least five proteins are available in all pairs. Separately use the original 11 measured proteins with their complete-case patient set when at least ten pairs are available. Process the all-peptide data as a technical sensitivity analysis.

The rank audit amendment was locked on 2026-10-20 00:57:50 +0800. Preserve the original F statistic in the audit and use `sign(logFC) * sqrt(pmax(F, 0))` for the common rank definition. Break ties by gene name.

## Step 5 closeout addendum

Locked on 2026-10-20 01:17:18 +0800. Check the 28 pairwise provenance relationships among the eight study sources. Reproduce the original anti-PD-1 raw models for the same formal patients and verify the influence calculations independently.

The independent paired protein sign-flip check uses R/dqrng PCG64 as specified in the technical amendment of 2026-10-20 01:20:38 +0800. The check uses seed 20261020, 100,000 draws, two-sided exceedances and the plus-one correction.
