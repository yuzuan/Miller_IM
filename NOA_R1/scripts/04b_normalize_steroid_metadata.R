#!/usr/bin/env Rscript
# Purpose: Apply the recorded capitalization correction to the frozen GBM-CARE steroid fields.
# Inputs: step4_steroids_geomx/{patient_deltas,GSE274546_steroid_patient_status}.csv;
# NOA_PUBLIC_DATA_ROOT/reference_metadata/GBM_CARE_WT/spitzer_supptable1.xlsx (Tumor ID, Time point,
# Patient ID, Steroid treatment before the surgery).
# Outputs in step4_steroids_geomx/clinical_case_normalized: GSE274546_steroid_patient_status_normalized.csv,
# GSE274546_steroid_group_summary_normalized.csv, changed_classification.csv and validation_checks.csv.
# Run: Rscript NOA_R1/scripts/04b_normalize_steroid_metadata.R NOA_R1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04b_normalize_steroid_metadata.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
step_dir <- file.path(revision, "step4_steroids_geomx")
output_dir <- file.path(step_dir, "clinical_case_normalized")
if (dir.exists(output_dir)) stop("Normalized clinical output directory already exists")
dir.create(output_dir)

delta_table <- read.csv(file.path(step_dir, "patient_deltas.csv"))
target_rows <- delta_table$dataset == "GSE274546" & delta_table$set_id == "target_top20"
delta_table <- delta_table[target_rows, , drop = FALSE]

# Load the source workbook locally and use the original prefixes for the formal patient-time pairs.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)
workbook_path <- file.path(public_root, "reference_metadata", "GBM_CARE_WT", "spitzer_supptable1.xlsx")
clinical_table <- as.data.frame(readxl::read_excel(workbook_path))
clinical_table$prefix <- sub("_[012]$", "", clinical_table$`Tumor ID`)
clinical_table$time <- trimws(as.character(clinical_table$`Time point`))
clinical_table <- clinical_table[clinical_table$time %in% c("T1", "T2"), , drop = FALSE]

check_rows <- list()

# Save the original correction check rows on the first failure. Returns nothing.
check <- function(name, value) {
    passed <- isTRUE(value)
    check_rows[[length(check_rows) + 1L]] <<- data.frame(check = name, passed = passed)
    if (!passed) {
        check_table_total <- do.call(rbind, check_rows)
        write.csv(check_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
        stop(name)
    }
}

check("clinical_prefix_time_unique", !anyDuplicated(paste(clinical_table$prefix, clinical_table$time)))

# Compare the trimmed lowercase text and restore the original Yes and No display values.
canonical <- tolower(trimws(as.character(clinical_table$`Steroid treatment before the surgery`)))
clinical_table$steroid <- ifelse(canonical == "yes", "Yes", ifelse(canonical == "no", "No", NA_character_))
primary_table <- clinical_table[clinical_table$time == "T1", , drop = FALSE]
recurrent_table <- clinical_table[clinical_table$time == "T2", , drop = FALSE]
prefixes <- sub("^GSE274546_", "", delta_table$patient_id)
primary_index <- match(prefixes, primary_table$prefix)
recurrent_index <- match(prefixes, recurrent_table$prefix)
all_mapped <- nrow(delta_table) == 45L && all(!is.na(primary_index)) && all(!is.na(recurrent_index)) &&
    all(primary_table$`Patient ID`[primary_index] == recurrent_table$`Patient ID`[recurrent_index])
check("formal_45_all_mapped", all_mapped)

# Keep the raw source cells and the unchanged paired scores alongside the normalized status.
status_table <- data.frame(
    patient_id          = delta_table$patient_id,
    clinical_patient_id = primary_table$`Patient ID`[primary_index],
    tumor_prefix        = prefixes,
    primary_raw         = primary_table$`Steroid treatment before the surgery`[primary_index],
    recurrent_raw       = recurrent_table$`Steroid treatment before the surgery`[recurrent_index],
    primary_steroid     = primary_table$steroid[primary_index],
    recurrent_steroid   = recurrent_table$steroid[recurrent_index],
    primary_score       = delta_table$primary_score,
    recurrent_score     = delta_table$recurrent_score,
    delta               = delta_table$delta
)
unknown <- is.na(status_table$primary_steroid) | is.na(status_table$recurrent_steroid)
transition <- paste(status_table$primary_steroid, status_table$recurrent_steroid, sep = "->")
status_table$group <- ifelse(unknown, "Missing", transition)
write.csv(
    status_table, file.path(output_dir, "GSE274546_steroid_patient_status_normalized.csv"),
    row.names = FALSE,
    na = ""
)

# Describe each recorded transition using its original within-patient differences.
group_tables <- list()
for (group in c("No->No", "No->Yes", "Yes->No", "Yes->Yes")) {
    values <- status_table$delta[status_table$group == group]
    interval <- if (length(values) > 1L) t.test(values)$conf.int else c(NA_real_, NA_real_)
    group_tables[[group]] <- data.frame(
        group      = group,
        n          = length(values),
        positive   = sum(values > 0),
        mean_delta = mean(values),
        ci_low     = interval[1],
        ci_high    = interval[2]
    )
}
group_table_total <- do.call(rbind, group_tables)
write.csv(
    group_table_total, file.path(output_dir, "GSE274546_steroid_group_summary_normalized.csv"),
    row.names = FALSE,
    na = ""
)

# Confirm the original amendment changed only the documented patient and reproduced the frozen partition.
check("complete_37", sum(status_table$group != "Missing") == 37L && sum(status_table$group == "Missing") == 8L)
baseline_table <- read.csv(file.path(step_dir, "GSE274546_steroid_patient_status.csv"))
baseline_groups <- baseline_table$group[match(status_table$patient_id, baseline_table$patient_id)]
changed <- status_table$patient_id[status_table$group != baseline_groups]
check("only_SM07_changed", identical(changed, "GSE274546_SM07"))
check("SM07_Yes_to_Yes", status_table$group[status_table$patient_id == "GSE274546_SM07"] == "Yes->Yes")
changed_table <- data.frame(
    patient_id      = changed,
    old_group       = baseline_table$group[match(changed, baseline_table$patient_id)],
    corrected_group = status_table$group[match(changed, status_table$patient_id)],
    reason          = "lowercase yes normalized to Yes"
)
write.csv(changed_table, file.path(output_dir, "changed_classification.csv"), row.names = FALSE, na = "")
check_table_total <- do.call(rbind, check_rows)
write.csv(check_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
cat(nrow(check_table_total), " passing clinical correction checks\n", sep = "")
