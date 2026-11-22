#!/usr/bin/env Rscript
# Purpose: Repeat paired edgeR/GSEA with one patient omitted at a time and check full data against Step 38.
# Inputs: Both 03_myeloid pseudobulk count/metadata sets, the config IDH crosswalk and the Step 38 gene set.
# Outputs: Figure1F_leave_one_patient_out_program_gsea.csv,
# Figure1F_leave_one_patient_out_raw20_gsea.csv, Figure1F_full_data_reproduction_check.csv,
# Figure1F_leave_one_patient_out_fold_audit.csv, in the requested output directory.
# Run: Rscript scripts/05e_leave_one_patient_out.R [output_dir] [step38_dir].
# The published analysis uses edgeR, fgsea, seed 42 for each GSEA, mclapply and SerialParam.
# Keep the original mclapply and mc.set.seed settings to reproduce the published random execution.

required_packages <- c("edgeR", "fgsea", "BiocParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required Bioconductor packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
    library(edgeR)
    library(fgsea)
    library(BiocParallel)
    library(parallel)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(dirname(script_path))
source(file.path(project_dir, "scripts", "00a_functions.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: 05e_leave_one_patient_out.R [output_dir] [step38_dir]")
output_dir <- if (length(args) >= 1L) {
    normalizePath(args[[1]], mustWork = FALSE)
} else {
    file.path(project_dir, "write", "46_figure1_ef_candidates", "Figure1", "source_data")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Check the frozen IDH pairing against the published eligibility rule, which retains 17 of 18 pairs.
idh_crosswalk_file <- file.path(
    project_dir, "config", "GSE174554_formal_pair_IDH_crosswalk.csv"
)

if (!file.exists(idh_crosswalk_file)) stop("Missing GSE174554 IDH crosswalk: ", idh_crosswalk_file)
idh_crosswalk_table <- read.csv(
    idh_crosswalk_file, stringsAsFactors = FALSE, check.names = FALSE
)
required_idh_columns <- c(
    "dataset", "pair_key", "primary_sample", "recurrent_sample",
    "primary_idh", "recurrent_idh", "include_idh_wildtype",
    "source_citation", "source_file_sha256"
)
if (!all(required_idh_columns %in% colnames(idh_crosswalk_table))) {
    stop("GSE174554 IDH crosswalk is missing required columns.")
}
idh_flag <- toupper(trimws(as.character(idh_crosswalk_table$include_idh_wildtype)))
if (any(!idh_flag %in% c("TRUE", "FALSE"))) stop("GSE174554 IDH eligibility contains values other than TRUE/FALSE.")
idh_crosswalk_table$include_idh_wildtype <- idh_flag == "TRUE"
if (
    nrow(idh_crosswalk_table) != 18L ||
        any(idh_crosswalk_table$dataset != "GSE174554") ||
        anyDuplicated(idh_crosswalk_table$pair_key) ||
        sum(idh_crosswalk_table$include_idh_wildtype) != 17L ||
        !identical(
            idh_crosswalk_table$pair_key[!idh_crosswalk_table$include_idh_wildtype],
            "GSE174554_pair33"
        )
        ||
        length(unique(idh_crosswalk_table$source_file_sha256)) != 1L ||
        unique(idh_crosswalk_table$source_file_sha256) !=
            "58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42"
) {
    stop("GSE174554 IDH crosswalk does not reproduce the audited 17/18 eligibility.")
}

# Read the Step 38 program in its original gene order.
step38_dir <- if (length(args) >= 2L) {
    normalizePath(args[[2]], mustWork = TRUE)
} else {
    file.path(project_dir, "write", "38_independent_cohort_mg_inflammatory_recalculation")
}
raw20_name <- "Miller_Microglial_Inflammatory_raw_top20"
provenance_path <- file.path(step38_dir, "fixed_program_provenance.csv")
provenance_table <- read.csv(provenance_path, stringsAsFactors = FALSE, check.names = FALSE)
raw20_table <- provenance_table[provenance_table$signature == raw20_name, , drop = FALSE]
if (nrow(raw20_table) != 1L) stop("The fixed Miller-IM program definition is missing or duplicated.")
fixed_sets <- setNames(list(strsplit(raw20_table$genes[1], ";", fixed = TRUE)[[1]]), raw20_name)

# Index the input directories by the published cohort names for the loader below.
independent_dirs <- c(
    GSE174554 = file.path(project_dir, "write", "34_gse174554_raw_independent_discovery", "03_myeloid"),
    GSE274546 = file.path(project_dir, "write", "36_gse274546_raw_independent_reannotation", "03_myeloid")
)

# Read counts and aligned metadata for one cohort. Stop if files, columns or sample ordering are invalid.
load_cohort <- function(dataset) {
    input_dir <- independent_dirs[[dataset]]
    count_path <- file.path(input_dir, "patient_condition_pseudobulk_counts.csv.gz")
    meta_path <- file.path(input_dir, "patient_condition_pseudobulk_metadata.csv")
    if (!file.exists(count_path) || !file.exists(meta_path)) stop("Missing pseudobulk input for ", dataset)
    counts <- read_count_matrix(count_path)
    metadata_table <- read.csv(meta_path, stringsAsFactors = FALSE)
    required_meta <- c("sample_id", "patient_id", "condition", "n_cells")
    if (!all(required_meta %in% colnames(metadata_table))) stop("Missing metadata columns for ", dataset)
    if (!identical(colnames(counts), metadata_table$sample_id)) stop(dataset, " counts and metadata are not aligned.")
    metadata_table <- data.frame(
        unit_id          = metadata_table$sample_id,
        pair_key         = metadata_table$patient_id,
        condition        = metadata_table$condition,
        n_myeloid_cells  = metadata_table$n_cells,
        stringsAsFactors = FALSE
    )
    result <- list(counts = counts, meta = metadata_table)
    return(result)
}

cohorts <- vector("list", length(independent_dirs))
names(cohorts) <- names(independent_dirs)
for (dataset in names(independent_dirs)) cohorts[[dataset]] <- load_cohort(dataset)
threshold <- 20L

# Filter metadata at the original cell minimum and return one primary and one recurrent sample per eligible patient.
eligible_pairs <- function(metadata_table, dataset) {
    keep <- !is.na(metadata_table$pair_key) & nzchar(metadata_table$pair_key) &
        metadata_table$n_myeloid_cells >= threshold & metadata_table$condition %in% c("Primary", "Recurrent")
    metadata_table <- metadata_table[keep, , drop = FALSE]
    pair_count_matrix <- table(metadata_table$pair_key, metadata_table$condition)
    complete <- rownames(pair_count_matrix)[pair_count_matrix[, "Primary"] == 1 & pair_count_matrix[, "Recurrent"] == 1]
    metadata_table <- metadata_table[metadata_table$pair_key %in% complete, , drop = FALSE]
    if (dataset != "GSE174554") return(metadata_table)
    missing_pairs <- setdiff(unique(metadata_table$pair_key), idh_crosswalk_table$pair_key)
    if (length(missing_pairs)) {
        stop(
            "GSE174554 eligible pairs are absent from the IDH crosswalk: ",
            paste(missing_pairs, collapse = ";")
        )
    }
    retained_pairs <- idh_crosswalk_table$pair_key[idh_crosswalk_table$include_idh_wildtype]
    metadata_table <- metadata_table[metadata_table$pair_key %in% retained_pairs, , drop = FALSE]
    return(metadata_table)
}

pair_counts_before <- setNames(integer(length(cohorts)), names(cohorts))
for (dataset in names(cohorts)) {
    eligible_table <- eligible_pairs(cohorts[[dataset]]$meta, dataset)
    pair_counts_before[dataset] <- length(unique(eligible_table$pair_key))
}
pair_counts_match <- identical(
    as.integer(pair_counts_before[c("GSE174554", "GSE274546")]),
    c(17L, 45L)
)
if (!pair_counts_match) stop("Formal IDH-restricted analyses must contain 17 and 45 pairs.")

# Run the full model first, then omit patients individually in the original patient order.
jobs <- list()
index <- 1L
for (dataset in names(independent_dirs)) {
    eligible_table <- eligible_pairs(cohorts[[dataset]]$meta, dataset)
    patients <- sort(unique(eligible_table$pair_key))
    jobs[[index]] <- list(dataset = dataset, omitted_patient = NA_character_)
    index <- index + 1L
    for (patient in patients) {
        jobs[[index]] <- list(dataset = dataset, omitted_patient = patient)
        index <- index + 1L
    }
}

# Fit paired edgeR/GSEA for one cohort with an optional omitted patient and return the result table.
fit_one <- function(cohort, dataset, omitted_patient = NA_character_) {
    metadata_table <- eligible_pairs(cohort$meta, dataset)
    if (!is.na(omitted_patient)) {
        metadata_table <- metadata_table[metadata_table$pair_key != omitted_patient, , drop = FALSE]
    }
    counts <- cohort$counts[, metadata_table$unit_id, drop = FALSE]
    metadata_table <- metadata_table[match(colnames(counts), metadata_table$unit_id), , drop = FALSE]
    metadata_table$condition <- factor(metadata_table$condition, levels = c("Primary", "Recurrent"))
    metadata_table$pair_key <- factor(metadata_table$pair_key)
    n_pairs <- nlevels(metadata_table$pair_key)

    # Use the published paired design and robust edgeR, with ranks defined as sign(logFC) × sqrt(F).
    design_matrix <- model.matrix(~ pair_key + condition, data = metadata_table)
    dge <- DGEList(counts = counts)
    keep <- filterByExpr(dge, design = design_matrix)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design_matrix, robust = TRUE)
    fit_model <- glmQLFit(dge, design_matrix, robust = TRUE)
    coefficient <- which(colnames(design_matrix) == "conditionRecurrent")
    if (length(coefficient) != 1) stop("Cannot find recurrent coefficient for ", dataset)
    qlf_model <- glmQLFTest(fit_model, coef = coefficient)
    deg_table <- topTags(qlf_model, n = Inf, sort.by = "PValue")$table
    rank_stat <- sign(deg_table$logFC) * sqrt(pmax(deg_table$F, 0))
    names(rank_stat) <- rownames(deg_table)
    rank_stat <- sort(rank_stat[is.finite(rank_stat)], decreasing = TRUE)

    # The shared enrichment function sets seed 42 here and retains the published result columns.
    result_table <- run_fixed_gsea(fixed_sets, rank_stat)
    result_table$nominal_P <- result_table$pval
    result_table$dataset <- dataset
    result_table$run_type <- if (is.na(omitted_patient)) "Full" else "Leave-one-patient-out"
    result_table$omitted_patient <- if (is.na(omitted_patient)) "None" else omitted_patient
    result_table$n_pairs <- n_pairs
    result_table$n_tested_genes <- length(rank_stat)
    return(result_table)
}

# Run a single leave-one-patient-out fit and return its table for the original mclapply call.
run_job <- function(job) {
    result_table <- fit_one(cohorts[[job$dataset]], job$dataset, job$omitted_patient)
    return(result_table)
}

workers <- if (.Platform$OS.type == "windows") 1L else min(9L, max(1L, detectCores() - 1L))
all_rows <- mclapply(
    jobs,
    run_job,
    mc.cores       = workers,
    mc.preschedule = FALSE,
    mc.set.seed    = FALSE
)

# Order results by cohort, full versus left-out fit and patient, with the published column order.
result_table_total <- do.call(rbind, all_rows)
column_order <- c(
    "dataset", "run_type", "omitted_patient", "n_pairs", "n_tested_genes", "pathway",
    "size", "ES", "NES", "nominal_P", "log2err", "leadingEdge"
)
result_table_total <- result_table_total[, column_order]
result_order <- order(
    match(result_table_total$dataset, names(independent_dirs)),
    match(result_table_total$run_type, c("Full", "Leave-one-patient-out")),
    result_table_total$omitted_patient,
    result_table_total$pathway
)
result_table_total <- result_table_total[result_order, ]

raw20_result_table <- result_table_total[result_table_total$pathway == raw20_name, , drop = FALSE]
expected_rows <- length(cohorts) + sum(pair_counts_before)
if (nrow(raw20_result_table) != expected_rows) {
    stop("Expected ", expected_rows, " raw20 rows, found ", nrow(raw20_result_table))
}

input_gene_counts <- setNames(integer(length(cohorts)), names(cohorts))
for (dataset in names(cohorts)) input_gene_counts[dataset] <- nrow(cohorts[[dataset]]$counts)
fold_columns <- c(
    "dataset", "run_type", "omitted_patient", "n_pairs", "n_tested_genes"
)
fold_audit_table <- unique(raw20_result_table[, fold_columns])
fold_audit_table$n_pairs_before <- unname(pair_counts_before[fold_audit_table$dataset])
fold_audit_table$n_pairs_after <- fold_audit_table$n_pairs
fold_audit_table$n_units <- 2L * fold_audit_table$n_pairs
fold_audit_table$n_genes_input <- unname(input_gene_counts[fold_audit_table$dataset])
fold_audit_table$design_rank <- fold_audit_table$n_pairs + 1L
fold_audit_table$residual_df <- fold_audit_table$n_pairs - 1L
fold_audit_table$coefficient <- "conditionRecurrent"
fold_audit_table$seed <- 42L
fold_audit_table$threshold_min_cells_per_endpoint <- threshold
fold_audit_table$workers <- workers
fold_columns <- c(
    "dataset", "run_type", "omitted_patient", "n_pairs_before", "n_pairs_after", "n_units",
    "n_genes_input", "n_tested_genes", "design_rank", "residual_df", "coefficient", "seed",
    "threshold_min_cells_per_endpoint", "workers"
)
fold_audit_table <- fold_audit_table[, fold_columns]

# Check that full-data results reproduce Step 38 within the original 1e-10 tolerance.
official_table <- read.csv(
    file.path(step38_dir, "independent_fixed_program_targeted_gsea.csv"),
    stringsAsFactors = FALSE
)
official_table <- official_table[
    official_table$dataset %in% names(independent_dirs) & official_table$threshold == threshold &
        official_table$pathway == raw20_name & official_table$formal_testable,
    ,
    drop = FALSE
]
full_raw20_table <- raw20_result_table[raw20_result_table$run_type == "Full", , drop = FALSE]
comparison_table_total <- merge(
    full_raw20_table[, c("dataset", "NES", "nominal_P")],
    official_table[, c("dataset", "NES", "pval")],
    by = "dataset",
    suffixes = c("_rerun", "_step38")
)
reproduction_differences <- c(
    abs(comparison_table_total$NES_rerun - comparison_table_total$NES_step38),
    abs(comparison_table_total$nominal_P - comparison_table_total$pval)
)
if (nrow(comparison_table_total) != 2L || any(reproduction_differences > 1e-10)) {
    stop("Full-data GSEA does not reproduce Step38.")
}

lopo_counts <- table(
    raw20_result_table$dataset[raw20_result_table$run_type == "Leave-one-patient-out"]
)
if (
    unname(lopo_counts[["GSE174554"]]) != 17L ||
        unname(lopo_counts[["GSE274546"]]) != 45L ||
        any(raw20_result_table$omitted_patient == "GSE174554_pair33")
) {
    stop("IDH-restricted LOPO folds do not reproduce the expected 17/45 design.")
}

# Write the four tables consumed downstream using the published filenames and columns.
write.csv(
    result_table_total,
    file.path(output_dir, "Figure1F_leave_one_patient_out_program_gsea.csv"),
    row.names = FALSE
)
write.csv(
    raw20_result_table,
    file.path(output_dir, "Figure1F_leave_one_patient_out_raw20_gsea.csv"),
    row.names = FALSE
)
write.csv(
    comparison_table_total,
    file.path(output_dir, "Figure1F_full_data_reproduction_check.csv"),
    row.names = FALSE
)
write.csv(
    fold_audit_table,
    file.path(output_dir, "Figure1F_leave_one_patient_out_fold_audit.csv"),
    row.names = FALSE
)

cat(
    "STEP46_LOPO_COMPLETE raw20_rows=", nrow(raw20_result_table),
    " lopo_runs=", sum(raw20_result_table$run_type == "Leave-one-patient-out"),
    "\n",
    sep = ""
)

# Check all LOPO files for existence and size. Row counts were checked against the patient tasks above.
output_names <- c(
    "Figure1F_leave_one_patient_out_program_gsea.csv", "Figure1F_leave_one_patient_out_raw20_gsea.csv",
    "Figure1F_full_data_reproduction_check.csv", "Figure1F_leave_one_patient_out_fold_audit.csv"
)
output_paths <- file.path(output_dir, output_names)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
