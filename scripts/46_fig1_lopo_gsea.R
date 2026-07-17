#!/usr/bin/env Rscript

required_packages <- c("edgeR", "fgsea", "BiocParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required Bioconductor packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(BiocParallel)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: 46_fig1_lopo_gsea.R [output_dir] [step38_dir]")
}
output_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = FALSE)
} else {
  file.path(project_dir, "write", "46_figure1_ef_candidates", "Figure1", "source_data")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

independent_dirs <- c(
  GSE174554 = file.path(project_dir, "write", "34_gse174554_raw_independent_discovery", "03_myeloid"),
  GSE274546 = file.path(project_dir, "write", "36_gse274546_raw_independent_reannotation", "03_myeloid")
)
step38_dir <- if (length(args) >= 2L) {
  normalizePath(args[[2]], mustWork = TRUE)
} else {
  file.path(project_dir, "write", "38_independent_cohort_mg_inflammatory_recalculation")
}
idh_crosswalk_file <- file.path(
  project_dir, "config", "GSE174554_formal_pair_IDH_crosswalk.csv"
)
threshold <- 20L
raw20_name <- "Miller_Microglial_Inflammatory_raw_top20"

if (!file.exists(idh_crosswalk_file)) {
  stop("Missing GSE174554 IDH crosswalk: ", idh_crosswalk_file)
}
idh_crosswalk <- read.csv(
  idh_crosswalk_file, stringsAsFactors = FALSE, check.names = FALSE
)
required_idh_columns <- c(
  "dataset", "pair_key", "primary_sample", "recurrent_sample",
  "primary_idh", "recurrent_idh", "include_idh_wildtype",
  "source_citation", "source_file_sha256"
)
if (!all(required_idh_columns %in% names(idh_crosswalk))) {
  stop("GSE174554 IDH crosswalk is missing required columns.")
}
idh_flag <- toupper(trimws(as.character(idh_crosswalk$include_idh_wildtype)))
if (any(!idh_flag %in% c("TRUE", "FALSE"))) {
  stop("GSE174554 IDH eligibility contains values other than TRUE/FALSE.")
}
idh_crosswalk$include_idh_wildtype <- idh_flag == "TRUE"
expected_supplement_sha256 <- "58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42"
if (
  nrow(idh_crosswalk) != 18L ||
    any(idh_crosswalk$dataset != "GSE174554") ||
    anyDuplicated(idh_crosswalk$pair_key) ||
    sum(idh_crosswalk$include_idh_wildtype) != 17L ||
    !identical(
      idh_crosswalk$pair_key[!idh_crosswalk$include_idh_wildtype],
      "GSE174554_pair33"
    ) ||
    length(unique(idh_crosswalk$source_file_sha256)) != 1L ||
    unique(idh_crosswalk$source_file_sha256) != expected_supplement_sha256
) {
  stop("GSE174554 IDH crosswalk does not reproduce the audited 17/18 eligibility.")
}

read_counts <- function(path) {
  x <- read.csv(gzfile(path), row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  if (any(!is.finite(x)) || any(abs(x - round(x)) > 1e-8)) stop("Non-integer counts in ", path)
  storage.mode(x) <- "integer"
  if (anyDuplicated(rownames(x))) stop("Duplicated genes in ", path)
  if (any(x < 0L)) stop("Negative counts in ", path)
  x
}

load_cohort <- function(dataset) {
  input_dir <- independent_dirs[[dataset]]
  count_path <- file.path(input_dir, "patient_condition_pseudobulk_counts.csv.gz")
  meta_path <- file.path(input_dir, "patient_condition_pseudobulk_metadata.csv")
  if (!file.exists(count_path) || !file.exists(meta_path)) stop("Missing pseudobulk input for ", dataset)
  counts <- read_counts(count_path)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  required_meta <- c("sample_id", "patient_id", "condition", "n_cells")
  if (!all(required_meta %in% names(meta))) stop("Missing metadata columns for ", dataset)
  if (!identical(colnames(counts), meta$sample_id)) stop(dataset, " counts and metadata are not aligned.")
  meta <- data.frame(
    unit_id = meta$sample_id,
    pair_key = meta$patient_id,
    condition = meta$condition,
    n_myeloid_cells = meta$n_cells,
    stringsAsFactors = FALSE
  )
  list(counts = counts, meta = meta)
}

eligible_pairs <- function(meta, dataset) {
  keep <- !is.na(meta$pair_key) & nzchar(meta$pair_key) &
    meta$n_myeloid_cells >= threshold & meta$condition %in% c("Primary", "Recurrent")
  meta <- meta[keep, , drop = FALSE]
  tab <- table(meta$pair_key, meta$condition)
  complete <- rownames(tab)[tab[, "Primary"] == 1 & tab[, "Recurrent"] == 1]
  meta <- meta[meta$pair_key %in% complete, , drop = FALSE]
  if (dataset != "GSE174554") return(meta)
  missing_pairs <- setdiff(unique(meta$pair_key), idh_crosswalk$pair_key)
  if (length(missing_pairs)) {
    stop(
      "GSE174554 eligible pairs are absent from the IDH crosswalk: ",
      paste(missing_pairs, collapse = ";")
    )
  }
  retained_pairs <- idh_crosswalk$pair_key[idh_crosswalk$include_idh_wildtype]
  meta[meta$pair_key %in% retained_pairs, , drop = FALSE]
}

provenance_path <- file.path(step38_dir, "fixed_program_provenance.csv")
provenance <- read.csv(provenance_path, stringsAsFactors = FALSE, check.names = FALSE)
raw20_row <- provenance[provenance$signature == raw20_name, , drop = FALSE]
if (nrow(raw20_row) != 1L) stop("The fixed Miller-IM program definition is missing or duplicated.")
fixed_sets <- setNames(list(strsplit(raw20_row$genes[[1]], ";", fixed = TRUE)[[1]]), raw20_name)

run_fgsea <- function(rank_stat) {
  tested <- lapply(fixed_sets, function(genes) intersect(unique(genes), names(rank_stat)))
  tested <- tested[lengths(tested) >= 5]
  set.seed(42)
  result <- fgseaMultilevel(
    pathways = tested,
    stats = rank_stat,
    minSize = 5,
    maxSize = 500,
    eps = 0,
    BPPARAM = SerialParam(progressbar = FALSE)
  )
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  result$nominal_P <- result$pval
  result
}

fit_one <- function(cohort, dataset, omitted_patient = NA_character_) {
  meta <- eligible_pairs(cohort$meta, dataset)
  if (!is.na(omitted_patient)) meta <- meta[meta$pair_key != omitted_patient, , drop = FALSE]
  counts <- cohort$counts[, meta$unit_id, drop = FALSE]
  meta <- meta[match(colnames(counts), meta$unit_id), , drop = FALSE]
  meta$condition <- factor(meta$condition, levels = c("Primary", "Recurrent"))
  meta$pair_key <- factor(meta$pair_key)
  n_pairs <- nlevels(meta$pair_key)
  design <- model.matrix(~ pair_key + condition, data = meta)
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  coefficient <- which(colnames(design) == "conditionRecurrent")
  if (length(coefficient) != 1) stop("Cannot find recurrent coefficient for ", dataset)
  qlf <- glmQLFTest(fit, coef = coefficient)
  deg <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  rank_stat <- sign(deg$logFC) * sqrt(pmax(deg$F, 0))
  names(rank_stat) <- rownames(deg)
  rank_stat <- sort(rank_stat[is.finite(rank_stat)], decreasing = TRUE)
  result <- run_fgsea(rank_stat)
  result$dataset <- dataset
  result$run_type <- if (is.na(omitted_patient)) "Full" else "Leave-one-patient-out"
  result$omitted_patient <- if (is.na(omitted_patient)) "None" else omitted_patient
  result$n_pairs <- n_pairs
  result$n_tested_genes <- length(rank_stat)
  result
}

cohorts <- setNames(lapply(names(independent_dirs), load_cohort), names(independent_dirs))
pair_counts_before <- vapply(
  names(cohorts),
  function(dataset) {
    length(unique(eligible_pairs(cohorts[[dataset]]$meta, dataset)$pair_key))
  },
  integer(1)
)
if (!identical(
  as.integer(pair_counts_before[c("GSE174554", "GSE274546")]),
  c(17L, 45L)
)) {
  stop("Formal IDH-restricted analyses must contain 17 and 45 pairs.")
}
jobs <- list()
index <- 1L
for (dataset in names(independent_dirs)) {
  eligible <- eligible_pairs(cohorts[[dataset]]$meta, dataset)
  patients <- sort(unique(eligible$pair_key))
  jobs[[index]] <- list(dataset = dataset, omitted_patient = NA_character_)
  index <- index + 1L
  for (patient in patients) {
    jobs[[index]] <- list(dataset = dataset, omitted_patient = patient)
    index <- index + 1L
  }
}

run_job <- function(job) {
  fit_one(cohorts[[job$dataset]], job$dataset, job$omitted_patient)
}

workers <- if (.Platform$OS.type == "windows") 1L else min(9L, max(1L, parallel::detectCores() - 1L))
all_rows <- parallel::mclapply(
  jobs,
  run_job,
  mc.cores = workers,
  mc.preschedule = FALSE,
  mc.set.seed = FALSE
)

all_results <- do.call(rbind, all_rows)
column_order <- c(
  "dataset", "run_type", "omitted_patient", "n_pairs", "n_tested_genes", "pathway",
  "size", "ES", "NES", "nominal_P", "log2err", "leadingEdge"
)
all_results <- all_results[, column_order]
all_results <- all_results[order(
  match(all_results$dataset, names(independent_dirs)),
  match(all_results$run_type, c("Full", "Leave-one-patient-out")),
  all_results$omitted_patient,
  all_results$pathway
), ]

raw20_results <- all_results[all_results$pathway == raw20_name, , drop = FALSE]
expected_rows <- 2L + 17L + 45L
if (nrow(raw20_results) != expected_rows) stop("Expected ", expected_rows, " raw20 rows, found ", nrow(raw20_results))

input_gene_counts <- vapply(cohorts, function(cohort) nrow(cohort$counts), integer(1))
fold_audit <- unique(raw20_results[, c(
  "dataset", "run_type", "omitted_patient", "n_pairs", "n_tested_genes"
)])
fold_audit$n_pairs_before <- unname(pair_counts_before[fold_audit$dataset])
fold_audit$n_pairs_after <- fold_audit$n_pairs
fold_audit$n_units <- 2L * fold_audit$n_pairs
fold_audit$n_genes_input <- unname(input_gene_counts[fold_audit$dataset])
fold_audit$design_rank <- fold_audit$n_pairs + 1L
fold_audit$residual_df <- fold_audit$n_pairs - 1L
fold_audit$coefficient <- "conditionRecurrent"
fold_audit$seed <- 42L
fold_audit$threshold_min_cells_per_endpoint <- threshold
fold_audit$workers <- workers
fold_audit <- fold_audit[, c(
  "dataset", "run_type", "omitted_patient", "n_pairs_before", "n_pairs_after", "n_units",
  "n_genes_input", "n_tested_genes", "design_rank", "residual_df", "coefficient", "seed",
  "threshold_min_cells_per_endpoint", "workers"
)]

official <- read.csv(file.path(step38_dir, "independent_fixed_program_targeted_gsea.csv"), stringsAsFactors = FALSE)
official <- official[
  official$dataset %in% names(independent_dirs) & official$threshold == threshold &
    official$pathway == raw20_name & official$formal_testable,
  ,
  drop = FALSE
]
full_raw20 <- raw20_results[raw20_results$run_type == "Full", , drop = FALSE]
comparison <- merge(
  full_raw20[, c("dataset", "NES", "nominal_P")],
  official[, c("dataset", "NES", "pval")],
  by = "dataset",
  suffixes = c("_rerun", "_step38")
)
reproduction_differences <- c(
  abs(comparison$NES_rerun - comparison$NES_step38),
  abs(comparison$nominal_P - comparison$pval)
)
if (nrow(comparison) != 2L || any(reproduction_differences > 1e-10)) {
  stop("Full-data GSEA does not reproduce Step38.")
}

lopo_counts <- table(
  raw20_results$dataset[raw20_results$run_type == "Leave-one-patient-out"]
)
if (
  unname(lopo_counts[["GSE174554"]]) != 17L ||
    unname(lopo_counts[["GSE274546"]]) != 45L ||
    any(raw20_results$omitted_patient == "GSE174554_pair33")
) {
  stop("IDH-restricted LOPO folds do not reproduce the expected 17/45 design.")
}

write.csv(
  all_results,
  file.path(output_dir, "Figure1F_leave_one_patient_out_program_gsea.csv"),
  row.names = FALSE
)
write.csv(
  raw20_results,
  file.path(output_dir, "Figure1F_leave_one_patient_out_raw20_gsea.csv"),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(output_dir, "Figure1F_full_data_reproduction_check.csv"),
  row.names = FALSE
)
write.csv(
  fold_audit,
  file.path(output_dir, "Figure1F_leave_one_patient_out_fold_audit.csv"),
  row.names = FALSE
)

cat(
  "STEP46_LOPO_COMPLETE raw20_rows=", nrow(raw20_results),
  " lopo_runs=", sum(raw20_results$run_type == "Leave-one-patient-out"),
  "\n",
  sep = ""
)
