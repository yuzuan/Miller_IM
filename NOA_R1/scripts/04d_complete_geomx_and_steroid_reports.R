#!/usr/bin/env Rscript
# Purpose: Complete the frozen clinical strata summary and add the GeoMx Welch statistic and degrees of freedom.
# Inputs: step4_steroids_geomx/clinical_case_normalized/GSE274546_steroid_patient_status_normalized.csv,
# step4_steroids_geomx/{GeoMx_formal_patient_source_copy,GeoMx_subgroup_interaction}.csv.
# Outputs in step4_steroids_geomx: GeoMx_Welch_statistic_df.csv and report_completion_checks.csv.
# Outputs in step4_steroids_geomx/clinical_case_normalized: steroid_summary_all_statuses.csv
# and steroid_missing_eight.csv.
# Run: Rscript NOA_R1/scripts/04d_complete_geomx_and_steroid_reports.R NOA_R1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04d_complete_geomx_and_steroid_reports.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
step_dir <- file.path(revision, "step4_steroids_geomx")
clinical_dir <- file.path(step_dir, "clinical_case_normalized")
output_paths <- c(
    strata = file.path(clinical_dir, "steroid_summary_all_statuses.csv"),
    missing = file.path(clinical_dir, "steroid_missing_eight.csv"),
    welch   = file.path(step_dir, "GeoMx_Welch_statistic_df.csv"),
    checks  = file.path(step_dir, "report_completion_checks.csv")
)
if (any(file.exists(output_paths))) stop("Report completion output already exists")

# Use the normalized 45-pair table. Preserve blank and NA input codes as missing.
status_table <- read.csv(
    file.path(clinical_dir, "GSE274546_steroid_patient_status_normalized.csv"),
    na.strings = c("", "NA")
)
check_rows <- list()

# Append one original report check and stop on failure. Returns nothing.
check <- function(name, value) {
    passed <- isTRUE(value)
    check_rows[[length(check_rows) + 1L]] <<- data.frame(check = name, passed = passed)
    if (!passed) stop(name)
}

check("formal_45_unique", nrow(status_table) == 45L && !anyDuplicated(status_table$patient_id))

# The four transitions, missing records, recorded records, and all pairs retain their original definitions.
group_tables <- list()
for (group in c("No->No", "No->Yes", "Yes->No", "Yes->Yes", "Missing", "Recorded", "All")) {
    values <- status_table$delta[status_table$group == group]
    if (group == "Recorded") values <- status_table$delta[status_table$group != "Missing"]
    if (group == "All") values <- status_table$delta
    paired_test <- t.test(values)
    group_tables[[group]] <- data.frame(
        group      = group,
        n          = length(values),
        positive   = sum(values > 0),
        mean_delta = mean(values),
        ci_low     = unname(paired_test$conf.int[1]),
        ci_high    = unname(paired_test$conf.int[2]),
        paired_t_p = paired_test$p.value,
        sum_delta  = sum(values)
    )
}

# These values are the fixed October reporting checkpoints for the original cohort.
recorded_table <- group_tables[["Recorded"]]
missing_table <- group_tables[["Missing"]]
all_table <- group_tables[["All"]]
check(
    "recorded_37",
    recorded_table$n == 37L && recorded_table$positive == 19L && abs(recorded_table$mean_delta - .05138585) < 1e-7
)
check(
    "missing_8",
    missing_table$n == 8L && missing_table$positive == 7L && abs(missing_table$mean_delta - .7045843) < 1e-7
)
check("missing_ci", abs(missing_table$ci_low - .1235453) < 1e-7 && abs(missing_table$ci_high - 1.285623) < 1e-6)
check("all_45", all_table$n == 45L && all_table$positive == 26L && abs(all_table$mean_delta - .167510016) < 1e-8)
strata_table_total <- do.call(rbind, group_tables)
subgroup_sums <- strata_table_total$sum_delta[strata_table_total$group %in% c("Recorded", "Missing")]
check("partition_sum", abs(sum(subgroup_sums) - all_table$sum_delta) < 1e-12)
missing_fields <- c("patient_id", "primary_raw", "recurrent_raw", "delta")
missing_patient_table <- status_table[status_table$group == "Missing", missing_fields, drop = FALSE]
check("missing_eight_unique", nrow(missing_patient_table) == 8L && !anyDuplicated(missing_patient_table$patient_id))

# Reproduce the Welch contrast from the formal Figure 3A patient differences.
geomx_table <- read.csv(file.path(step_dir, "GeoMx_formal_patient_source_copy.csv"))
nivolumab <- geomx_table$delta[geomx_table$trial_setting == "Nivolumab"]
control <- geomx_table$delta[geomx_table$trial_setting == "Control"]
formal <- nrow(geomx_table) == 22L && length(nivolumab) == 12L && length(control) == 10L &&
    all(geomx_table$entry == "strict_pass_idhwt")
check("geomx_formal_22", formal)
welch_test <- t.test(nivolumab, control, var.equal = FALSE)
nivolumab_term <- var(nivolumab) / length(nivolumab)
control_term <- var(control) / length(control)
variance <- nivolumab_term + control_term
estimate <- mean(nivolumab) - mean(control)
manual_t <- estimate / sqrt(variance)
manual_df <- variance^2 / (nivolumab_term^2 / (length(nivolumab) - 1) + control_term^2 / (length(control) - 1))

# Check the independently calculated t and df and the original effect, interval and P value.
check(
    "Welch_independent_t_df",
    abs(manual_t - unname(welch_test$statistic)) < 1e-12 && abs(manual_df - unname(welch_test$parameter)) < 1e-12
)
check("Welch_reference_t_df", abs(manual_t - 1.0541) < .0001 && abs(manual_df - 19.7739) < .0001)
reference_table <- read.csv(file.path(step_dir, "GeoMx_subgroup_interaction.csv"))
reference_table <- reference_table[reference_table$contrast == "Nivolumab_minus_Control_interaction", , drop = FALSE]
reference_interval <- c(reference_table$ci_low, reference_table$ci_high)
reproduces <- nrow(reference_table) == 1L && abs(reference_table$estimate - estimate) < 1e-12 &&
    max(abs(reference_interval - welch_test$conf.int)) < 1e-12 &&
    abs(reference_table$p_value - welch_test$p.value) < 1e-12
check("Welch_old_effect_ci_p", reproduces)

# Retain the precise interval and P from t.test and the independently checked Welch statistics.
welch_table <- data.frame(
    contrast    = "Nivolumab_minus_Control",
    n_nivolumab = length(nivolumab),
    n_control   = length(control),
    estimate    = estimate,
    se          = sqrt(variance),
    t           = manual_t,
    df          = manual_df,
    ci_low      = welch_test$conf.int[1],
    ci_high     = welch_test$conf.int[2],
    p_value     = welch_test$p.value,
    source      = "formal Figure3A patient source; no upstream score reconstruction"
)
write.csv(strata_table_total, output_paths[["strata"]], row.names = FALSE, na = "")
write.csv(missing_patient_table, output_paths[["missing"]], row.names = FALSE, na = "")
write.csv(welch_table, output_paths[["welch"]], row.names = FALSE)
check_table_total <- do.call(rbind, check_rows)
write.csv(check_table_total, output_paths[["checks"]], row.names = FALSE)
cat(nrow(check_table_total), " passing report completion checks\n", sep = "")
