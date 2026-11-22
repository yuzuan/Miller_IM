#!/usr/bin/env Rscript
# Purpose: Independently verify the corrected steroid assignments and paired changes from the source records.
# Inputs: step4_steroids_geomx/sample_scores.csv and clinical_case_normalized/
# {GSE274546_steroid_patient_status_normalized,GSE274546_steroid_group_summary_normalized}.csv;
# NOA_PUBLIC_DATA_ROOT/reference_metadata/GBM_CARE_WT/spitzer_supptable1.xlsx.
# Output: step4_steroids_geomx/clinical_case_normalized/independent_checks_v2.csv.
# Run: Rscript NOA_R1/scripts/04e_check_steroid_metadata.R NOA_R1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04e_check_steroid_metadata.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
step_dir <- file.path(revision, "step4_steroids_geomx")
clinical_dir <- file.path(step_dir, "clinical_case_normalized")
output_path <- file.path(clinical_dir, "independent_checks_v2.csv")
if (file.exists(output_path)) stop("Independent clinical output already exists")

# Read the source workbook without using the normalized status to choose workbook rows.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)
workbook_path <- file.path(public_root, "reference_metadata", "GBM_CARE_WT", "spitzer_supptable1.xlsx")
clinical_table <- as.data.frame(readxl::read_excel(workbook_path))
clinical_table$prefix <- sub("_[012]$", "", clinical_table$`Tumor ID`)
clinical_table <- clinical_table[clinical_table$`Time point` %in% c("T1", "T2"), , drop = FALSE]

# Pair original sample scores by patient and retain the normalized rows for comparison.
score_table <- read.csv(file.path(step_dir, "sample_scores.csv"))
target_rows <- score_table$dataset == "GSE274546" & score_table$set_id == "target_top20"
score_table <- score_table[target_rows, , drop = FALSE]
primary_table <- score_table[score_table$condition == "Primary", , drop = FALSE]
recurrent_table <- score_table[score_table$condition == "Recurrent", , drop = FALSE]
recurrent_table <- recurrent_table[match(primary_table$patient_id, recurrent_table$patient_id), , drop = FALSE]
fixed_table <- read.csv(file.path(clinical_dir, "GSE274546_steroid_patient_status_normalized.csv"))
check_rows <- list()

# Keep one row per audit check for the final independent review table. Returns nothing.
add_check <- function(name, value) {
    check_rows[[length(check_rows) + 1L]] <<- data.frame(check = name, passed = isTRUE(value))
}

# Reapply the case rule to source T1 and T2 and compare with the two recorded endpoints.
for (i in seq_len(nrow(primary_table))) {
    patient <- primary_table$patient_id[i]
    prefix <- sub("^GSE274546_", "", patient)
    selected_table <- clinical_table[clinical_table$prefix == prefix, , drop = FALSE]
    selected_table <- selected_table[match(c("T1", "T2"), selected_table$`Time point`), , drop = FALSE]
    canonical <- tolower(trimws(as.character(selected_table$`Steroid treatment before the surgery`)))
    expected <- ifelse(canonical == "yes", "Yes", ifelse(canonical == "no", "No", NA_character_))
    patient_table <- fixed_table[fixed_table$patient_id == patient, , drop = FALSE]
    observed <- as.character(c(patient_table$primary_steroid, patient_table$recurrent_steroid))
    observed[trimws(observed) == ""] <- NA_character_
    expected_delta <- recurrent_table$score[i] - primary_table$score[i]
    matches <- nrow(patient_table) == 1L && identical(unname(observed), unname(as.character(expected))) &&
        abs(patient_table$delta - expected_delta) < 1e-10
    add_check(paste("clinical_from_raw_and_sample_scores", patient), matches)
}

# Confirm the original complete count and recompute each recorded transition's t interval.
add_check("37_complete", sum(fixed_table$group != "Missing") == 37L)
summary_table <- read.csv(file.path(clinical_dir, "GSE274546_steroid_group_summary_normalized.csv"))
for (i in seq_len(nrow(summary_table))) {
    group_row_table <- summary_table[i, , drop = FALSE]
    values <- fixed_table$delta[fixed_table$group == group_row_table$group]
    paired_test <- t.test(values)
    interval <- c(group_row_table$ci_low, group_row_table$ci_high)
    matches <- length(values) == group_row_table$n &&
        sum(values > 0) == group_row_table$positive &&
        abs(mean(values) - group_row_table$mean_delta) < 1e-12 &&
        max(abs(paired_test$conf.int - interval)) < 1e-12
    add_check(paste("independent_group", group_row_table$group), matches)
}

check_table_total <- do.call(rbind, check_rows)
write.csv(check_table_total, output_path, row.names = FALSE)
failed <- sum(!check_table_total$passed)
cat(nrow(check_table_total) - failed, " passing independent clinical checks, ", failed, " failed\n", sep = "")
if (failed > 0L) quit(status = 1)
