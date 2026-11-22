# Step 4 analysis specification

The protocol used by the steroid and GeoMx analyses was locked on 2026-10-20 00:31:53 +0800.

Compare the original twenty-gene program, the fixed definition excluding SGK1 and PDK4, and the MSigDB 2025.1 human corticosteroid-response set on each cohort's Step 1 paired voom model. Use estimated inter-gene correlation for the primary CAMERA results, fixed 0.01 for sensitivity, and BH across the three sets per cohort and mode.

Summarize score changes by the two recorded preoperative steroid statuses in GSE274546. For the 22 GeoMx patient pairs, calculate the within-group effects for the recorded Nivolumab and Control groups and use a Welch test for the difference in patient-level changes.

The clinical field normalization amendment was locked on 2026-10-20 00:35:46 +0800. Trim whitespace and normalize case before mapping the affirmative and negative steroid fields to Yes and No. Other values are missing.
