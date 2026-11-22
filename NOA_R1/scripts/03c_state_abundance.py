#!/usr/bin/env python3
"""Estimate paired changes in captured myeloid-state proportions with the original NumPy sign-flip stream.

Input: step3_state_analysis/patient_state_cell_counts.csv (patient, condition, state, fraction).
Outputs in step3_state_analysis: state_abundance_statistics.csv, patient_state_abundance_changes.csv,
abundance_checks.csv, abundance_parameters.json.
Run: python NOA_R1/scripts/03c_state_abundance.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def main():
    """Read the state fractions and write the abundance results. Returns nothing."""
    output_dir = Path(__file__).resolve().parents[1] / "step3_state_analysis"
    counts_table = pd.read_csv(output_dir / "patient_state_cell_counts.csv")

    # Retain the nine fine states, one tested broad fraction, and the original cohort order.
    fine_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]
    tested_states = fine_states + ["MCG-like"]
    fine_state_count = len(fine_states)
    expected_pairs = {"GSE174554": 17, "GSE274546": 45}
    summary_rows = []
    change_rows = []
    check_rows = []

    def check(dataset, name, passed):
        """Record one frozen abundance check and stop on failure. Returns nothing."""
        check_rows.append(dict(dataset=dataset, check=name, passed=bool(passed)))
        pd.DataFrame(check_rows).to_csv(output_dir / "abundance_checks.csv", index=False)
        if not passed:
            raise ValueError(f"{dataset}: {name}")

    # Each matrix row is one formal patient, with explicit zeros for uncaptured states.
    for dataset_index, (dataset, expected_pair_count) in enumerate(expected_pairs.items()):
        cohort_table = counts_table[counts_table.dataset == dataset]
        patients = sorted(cohort_table.patient_id.unique())
        pair_count = len(patients)
        fractions_by_condition = {}
        for condition in ("Primary", "Recurrent"):
            selected_rows = (cohort_table.condition == condition) & cohort_table.state.isin(tested_states)
            selected_table = cohort_table[selected_rows]
            wide_table = selected_table.pivot(index="patient_id", columns="state", values="fraction")
            fraction_matrix = wide_table.loc[patients, tested_states].to_numpy()
            complete = (
                fraction_matrix.shape == (pair_count, len(tested_states))
                and np.isfinite(fraction_matrix).all()
                and (fraction_matrix >= 0).all()
                and (fraction_matrix <= 1).all()
            )
            check(dataset, f"{condition}_complete_proportions", complete)

            # Check the all-myeloid denominator and the published MCG definition (MCG1 to MCG5).
            fine_sum = fraction_matrix[:, :fine_state_count].sum(axis=1)
            check(dataset, f"{condition}_fine_proportions_sum_one", np.allclose(fine_sum, 1, atol=1e-12))
            mcg_sum = fraction_matrix[:, :5].sum(axis=1)
            broad_fraction = fraction_matrix[:, fine_state_count]
            check(dataset, f"{condition}_broad_matches_fine_sum", np.allclose(mcg_sum, broad_fraction, atol=1e-12))
            fractions_by_condition[condition] = fraction_matrix

        # The same patient signs are applied across states and the test uses absolute sums.
        primary_fractions = fractions_by_condition["Primary"]
        recurrent_fractions = fractions_by_condition["Recurrent"]
        delta = 100 * (recurrent_fractions - primary_fractions)
        observed = np.abs(delta.sum(axis=0))
        tolerance = 1e-12 * np.maximum(1, observed)
        extreme = np.zeros(len(tested_states), dtype=np.int64)

        # Enumerate up to 18 pairs, otherwise use 100,000 original seeded draws in 4,096-row batches.
        exact = pair_count <= 18
        trials = 2**pair_count if exact else 100000
        seed = 202610191 + dataset_index
        rng = np.random.default_rng(seed)
        for start in range(0, trials, 4096):
            batch_size = min(4096, trials - start)
            if exact:
                integers = np.arange(start, start + batch_size, dtype=np.uint64)
                patient_bits = np.arange(pair_count, dtype=np.uint64)
                signs = ((integers[:, None] >> patient_bits) & 1).astype(np.int8) * 2 - 1
            else:
                signs = rng.integers(0, 2, size=(batch_size, pair_count), dtype=np.int8) * 2 - 1

            permuted = (signs[:, :, None] * delta[None, :, :]).sum(axis=1)
            extreme += (np.abs(permuted) >= observed[None, :] - tolerance[None, :]).sum(axis=0)

        # Apply the Monte Carlo plus-one correction and adjust only the nine fine-state tests.
        p_values = extreme / trials if exact else (extreme + 1) / (trials + 1)
        valid_p = np.isfinite(p_values).all() and (p_values >= 0).all() and (p_values <= 1).all()
        check(dataset, "signflip_probabilities_valid", valid_p)
        by_adjusted = multipletests(p_values[:fine_state_count], method="fdr_by")[1]
        bh_adjusted = multipletests(p_values[:fine_state_count], method="fdr_bh")[1]

        # An unadjusted paired-t interval describes the patient-level percentage-point changes.
        for state_index, state in enumerate(tested_states):
            values = delta[:, state_index]
            mean_delta = float(values.mean())
            sd_delta = float(values.std(ddof=1))
            standard_error = sd_delta / np.sqrt(pair_count)
            half_width = stats.t.ppf(.975, pair_count - 1) * standard_error
            if standard_error > 0:
                t_p = float(2 * stats.t.sf(abs(mean_delta / standard_error), pair_count - 1))
            else:
                t_p = 1.0 if mean_delta == 0 else np.nan

            fine = state_index < fine_state_count
            level = "fine" if fine else "broad"
            summary_row = dict(
                dataset                 = dataset,
                state                   = state,
                level                   = level,
                n_pairs                 = pair_count,
                mean_delta_pp           = mean_delta,
                median_delta_pp         = float(np.median(values)),
                sd_delta_pp             = sd_delta,
                ci_low_pp               = mean_delta - half_width,
                ci_high_pp              = mean_delta + half_width,
                t_p_value               = t_p,
                positive_pairs          = int((values > 0).sum()),
                negative_pairs          = int((values < 0).sum()),
                zero_pairs              = int((values == 0).sum()),
                signflip_p              = float(p_values[state_index]),
                signflip_extreme        = int(extreme[state_index]),
                signflip_trials         = trials,
                signflip_method         = "exact_signflip" if exact else "monte_carlo_signflip",
                seed                    = seed,
                primary_adjusted_p      = float(by_adjusted[state_index]) if fine else float(p_values[state_index]),
                bh_sensitivity          = float(bh_adjusted[state_index]) if fine else float(p_values[state_index]),
                family_size             = fine_state_count if fine else 1,
                adjustment              = "BY9" if fine else "single_prespecified_broad_proportion",
                primary_mean_fraction   = float(primary_fractions[:, state_index].mean()),
                recurrent_mean_fraction = float(recurrent_fractions[:, state_index].mean()),
                interpretation          = "relative fraction among captured quality-filtered myeloid cells"
            )
            summary_rows.append(summary_row)

            # Keep one row per cohort, state, and patient, including zero changes.
            for patient_index, patient in enumerate(patients):
                change_row = dict(
                    dataset            = dataset,
                    state              = state,
                    level              = level,
                    patient_id         = patient,
                    primary_fraction   = primary_fractions[patient_index, state_index],
                    recurrent_fraction = recurrent_fractions[patient_index, state_index],
                    delta_pp           = values[patient_index]
                )
                change_rows.append(change_row)

        # Display the complement of the tested MCG-like fraction with its reflected change and interval.
        complement_row = summary_rows[-1].copy()
        complement_row.update(
            state                   = "MAC_MDSC-like",
            mean_delta_pp           = -complement_row["mean_delta_pp"],
            median_delta_pp         = -complement_row["median_delta_pp"],
            ci_low_pp               = -complement_row["ci_high_pp"],
            ci_high_pp              = -complement_row["ci_low_pp"],
            positive_pairs          = complement_row["negative_pairs"],
            negative_pairs          = complement_row["positive_pairs"],
            primary_mean_fraction   = 1 - complement_row["primary_mean_fraction"],
            recurrent_mean_fraction = 1 - complement_row["recurrent_mean_fraction"],
            adjustment              = "complement_of_same_broad_hypothesis_not_independent"
        )
        summary_rows.append(complement_row)
        for patient_index, patient in enumerate(patients):
            change_row = dict(
                dataset            = dataset,
                state              = "MAC_MDSC-like",
                level              = "broad",
                patient_id         = patient,
                primary_fraction   = 1 - primary_fractions[patient_index, fine_state_count],
                recurrent_fraction = 1 - recurrent_fractions[patient_index, fine_state_count],
                delta_pp           = -delta[patient_index, fine_state_count]
            )
            change_rows.append(change_row)

        check(dataset, "all_patients_retained", pair_count == expected_pair_count)
        print(f"{dataset}: {pair_count} pairs, {trials} sign patterns", flush=True)

    # Check both cohorts and all eleven displayed groups before exporting the two tables.
    summary_table_total = pd.DataFrame(summary_rows)
    change_table_total = pd.DataFrame(change_rows)
    display_state_count = len(tested_states) + 1
    complete = (
        len(summary_table_total) == len(expected_pairs) * display_state_count
        and len(change_table_total) == sum(expected_pairs.values()) * display_state_count
    )
    check("Both", "two_cohorts_eleven_display_groups", complete)
    summary_table_total.to_csv(output_dir / "state_abundance_statistics.csv", index=False)
    change_table_total.to_csv(output_dir / "patient_state_abundance_changes.csv", index=False)

    # Preserve the metadata fields and values paired with the original abundance results.
    parameters = dict(
        date            = "2026-10-19",
        denominator     = "all captured quality-filtered myeloid cells per patient/timepoint",
        exact_limit     = 18,
        monte_carlo     = 100000,
        seeds           = [202610191, 202610192],
        null_assumption = "patient difference signs exchangeable under null; not randomized timepoints",
        statistic       = "absolute sum of paired percentage-point differences",
        ci              = "paired t 95%, all groups, unadjusted descriptive intervals",
        family          = "nine fine states BY9 per cohort; broad MCG proportion one hypothesis; complementary group "
            "not independent"
    )
    (output_dir / "abundance_parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    print(f"Abundance checks passed: {len(check_rows)}", flush=True)


if __name__ == "__main__":
    main()
