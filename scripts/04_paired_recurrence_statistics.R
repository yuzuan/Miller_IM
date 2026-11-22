#!/usr/bin/env Rscript
# Purpose: Reproduce the paired fixed-program scores, edgeR and GSEA in the independent cohorts.
# Inputs: write/34_gse174554_raw_independent_discovery and write/36_gse274546_raw_independent_reannotation
# 03_myeloid/patient_condition_pseudobulk_* files, config/miller_im_gene_set.tsv and the IDH crosswalk.
# Outputs in the requested directory: fixed_program_provenance.csv, GSE174554_IDH_filter_audit.csv,
# and the independent model and score tables (independent_*.csv).
# Run with Rscript scripts/04_paired_recurrence_statistics.R <output_dir>.
# Use the original R/edgeR/fgsea workflow, 20/50-cell thresholds, seeds and serial BiocParallel.

suppressPackageStartupMessages({
    library(edgeR)
    library(fgsea)
    library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 04_paired_recurrence_statistics.R <output_dir>")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(dirname(script_path))
source(file.path(project_dir, "scripts", "00a_functions.R"))
output_dir <- normalizePath(args[[1]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

independent_dirs <- c(
    GSE174554 = file.path(project_dir, "write", "34_gse174554_raw_independent_discovery", "03_myeloid"),
    GSE274546 = file.path(project_dir, "write", "36_gse274546_raw_independent_reannotation", "03_myeloid")
)
miller_file <- file.path(project_dir, "config", "miller_im_gene_set.tsv")
idh_crosswalk_file <- file.path(
    project_dir, "config", "GSE174554_formal_pair_IDH_crosswalk.csv"
)

required_files <- c(miller_file, idh_crosswalk_file)
for (input_dir in independent_dirs) {
    required_files <- c(
        required_files,
        file.path(input_dir, "patient_condition_pseudobulk_counts.csv.gz"),
        file.path(input_dir, "patient_condition_pseudobulk_metadata.csv")
    )
}
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = "; "))

idh_crosswalk_table <- read.csv(
    idh_crosswalk_file, stringsAsFactors = FALSE, check.names = FALSE
)
required_idh_columns <- c(
    "dataset", "pair_key", "geo_pair_number", "primary_sample", "recurrent_sample",
    "primary_idh", "recurrent_idh", "include_idh_wildtype",
    "wang_supp_pair_number", "pair_stage_source", "idh_source", "source_citation",
    "source_url", "source_file_sha256", "stage_mapping_note"
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
) {
    stop("GSE174554 IDH crosswalk does not encode the audited 17/18 eligibility.")
}
if (
    any(!nzchar(idh_crosswalk_table$primary_sample)) ||
        any(!nzchar(idh_crosswalk_table$recurrent_sample)) ||
        any(!idh_crosswalk_table$primary_idh %in% c("IDH wildtype", "IDH mutant")) ||
        any(!idh_crosswalk_table$recurrent_idh %in% c("IDH wildtype", "IDH mutant"))
) {
    stop("GSE174554 IDH crosswalk contains incomplete or unexpected sample annotations.")
}
expected_supplement_sha256 <- "58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42"
if (
    length(unique(idh_crosswalk_table$source_file_sha256)) != 1L ||
        unique(idh_crosswalk_table$source_file_sha256) != expected_supplement_sha256
) {
    stop("GSE174554 IDH crosswalk does not retain the audited supplementary-workbook hash.")
}

miller_table <- read.delim(miller_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(colnames(miller_table), c("rank", "gene")) || nrow(miller_table) != 20L) {
    stop("config/miller_im_gene_set.tsv must contain exactly rank and gene columns with 20 rows.")
}
miller_table <- miller_table[order(miller_table$rank), , drop = FALSE]
raw20 <- miller_table$gene

expected_raw20 <- c(
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)
if (!identical(raw20, expected_raw20)) stop("Miller raw top20 changed from the audited definition.")

fixed_sets <- list(
    Miller_Microglial_Inflammatory_raw_top20 = raw20
)

provenance_table <- data.frame(
    signature = names(fixed_sets),
    source    = "First 20 genes of the Miller et al. Nature 2025 Supplementary Table 2 cNMF ranking",
    n_genes   = lengths(fixed_sets),
    genes     = vapply(fixed_sets, paste, collapse = ";", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
)
write.csv(provenance_table, file.path(output_dir, "fixed_program_provenance.csv"), row.names = FALSE)

# Read a cohort's count matrix and metadata, align them and stop on missing columns or invalid counts.
load_independent_cohort <- function(dataset) {
    input_dir <- independent_dirs[[dataset]]
    count_path <- file.path(input_dir, "patient_condition_pseudobulk_counts.csv.gz")
    meta_path <- file.path(input_dir, "patient_condition_pseudobulk_metadata.csv")
    count_matrix <- read_count_matrix(count_path)
    sample_table <- read.csv(meta_path, stringsAsFactors = FALSE)
    required_meta <- c("sample_id", "patient_id", "condition", "n_cells")
    if (!all(required_meta %in% colnames(sample_table))) stop("Missing metadata columns for ", dataset)
    if (!identical(colnames(count_matrix), sample_table$sample_id)) {
        stop(dataset, " counts and metadata are not aligned.")
    }
    if (anyDuplicated(sample_table$sample_id)) stop("Duplicated pseudobulk sample IDs for ", dataset)
    if (!all(sample_table$condition %in% c("Primary", "Recurrent"))) stop("Unexpected condition for ", dataset)
    if (anyNA(sample_table$patient_id) || any(!nzchar(sample_table$patient_id))) {
        keep <- !is.na(sample_table$patient_id) & nzchar(sample_table$patient_id)
        count_matrix <- count_matrix[, keep, drop = FALSE]
        sample_table <- sample_table[keep, , drop = FALSE]
    }
    source_label <- if (dataset == "GSE174554") {
        "Step34 independently reconstructed and annotated clean myeloid"
    } else {
        "Step36 independently reconstructed and annotated clean myeloid"
    }
    metadata_table <- data.frame(
        unit_id          = sample_table$sample_id,
        dataset          = dataset,
        pair_key         = sample_table$patient_id,
        condition        = sample_table$condition,
        n_myeloid_cells  = sample_table$n_cells,
        source           = source_label,
        count_path       = normalizePath(count_path),
        metadata_path    = normalizePath(meta_path),
        stringsAsFactors = FALSE
    )
    result <- list(counts = count_matrix, meta = metadata_table)
    return(result)
}

cohorts <- vector("list", length(independent_dirs))
names(cohorts) <- names(independent_dirs)
for (dataset in names(independent_dirs)) cohorts[[dataset]] <- load_independent_cohort(dataset)

# From metadata and the cell minimum, return eligible patients with one primary and one recurrent sample.
eligible_pairs <- function(metadata_table, threshold) {
    sample_table <- metadata_table[
        metadata_table$n_myeloid_cells >= threshold &
            metadata_table$condition %in% c("Primary", "Recurrent"),
        ,
        drop = FALSE
    ]
    pair_count_matrix <- table(sample_table$pair_key, sample_table$condition)
    complete <- rownames(pair_count_matrix)[pair_count_matrix[, "Primary"] == 1 & pair_count_matrix[, "Recurrent"] == 1]
    eligible_table <- sample_table[sample_table$pair_key %in% complete, , drop = FALSE]
    return(eligible_table)
}

# From metadata and the cell minimum, return each pair's original IDH eligibility. Stop on an unmatched pair.
build_idh_filter_audit <- function(metadata_table, threshold) {
    eligible_before_table <- eligible_pairs(metadata_table, threshold)
    pair_keys <- sort(unique(eligible_before_table$pair_key))
    missing_pairs <- setdiff(pair_keys, idh_crosswalk_table$pair_key)
    if (length(missing_pairs)) stop("Cannot audit GSE174554 IDH status for: ", paste(missing_pairs, collapse = ";"))
    audit_table <- idh_crosswalk_table[match(pair_keys, idh_crosswalk_table$pair_key), , drop = FALSE]
    primary_table <- eligible_before_table[eligible_before_table$condition == "Primary", , drop = FALSE]
    recurrent_table <- eligible_before_table[eligible_before_table$condition == "Recurrent", , drop = FALSE]
    audit_table$threshold <- threshold
    audit_table$primary_n_myeloid_cells <- primary_table$n_myeloid_cells[
        match(audit_table$pair_key, primary_table$pair_key)
    ]
    audit_table$recurrent_n_myeloid_cells <- recurrent_table$n_myeloid_cells[
        match(audit_table$pair_key, recurrent_table$pair_key)
    ]
    audit_table$eligible_before_idh_filter <- TRUE
    audit_table$retained_after_idh_filter <- audit_table$include_idh_wildtype
    audit_table$exclusion_reason <- ifelse(
        audit_table$retained_after_idh_filter,
        "",
        "Primary and recurrent samples are IDH mutant in Wang Supplementary Table 1"
    )
    return(audit_table)
}

idh_filter_audit_table <- build_idh_filter_audit(cohorts[["GSE174554"]]$meta, 20L)
if (
    nrow(idh_filter_audit_table) != 18L ||
        sum(idh_filter_audit_table$retained_after_idh_filter) != 17L ||
        idh_filter_audit_table$primary_n_myeloid_cells[
            idh_filter_audit_table$pair_key == "GSE174554_pair33"
        ] != 151L ||
        idh_filter_audit_table$recurrent_n_myeloid_cells[
            idh_filter_audit_table$pair_key == "GSE174554_pair33"
        ] != 477L
) {
    stop("GSE174554 formal-pair IDH audit does not reproduce the expected 17/18 filter.")
}
pair_audit <- list()
score_summary <- list()
score_units <- list()
deg_results <- list()
targeted_results <- list()
gene_results <- list()

# Select eligible samples from the cohort metadata using the cell minimum and the original GSE174554 IDH crosswalk.
eligible_idhwt_pairs <- function(metadata_table, dataset, threshold) {
    eligible_table <- eligible_pairs(metadata_table, threshold)
    if (dataset != "GSE174554") return(eligible_table)
    missing_pairs <- setdiff(unique(eligible_table$pair_key), idh_crosswalk_table$pair_key)
    if (length(missing_pairs)) {
        stop(
            "GSE174554 eligible pairs are absent from the IDH crosswalk: ",
            paste(missing_pairs, collapse = ";")
        )
    }
    retained_pairs <- idh_crosswalk_table$pair_key[idh_crosswalk_table$include_idh_wildtype]
    eligible_table <- eligible_table[eligible_table$pair_key %in% retained_pairs, , drop = FALSE]
    return(eligible_table)
}

# Use paired differences and the original seed. Enumerate all signs through 20 pairs and return P, count and flag.
sign_flip_test <- function(delta, seed) {
    delta <- delta[is.finite(delta)]
    n <- length(delta)
    observed <- abs(mean(delta))
    if (!n) {
        result <- c(p_value = NA_real_, n_permutations = 0, exact = FALSE)
        return(result)
    }
    if (n <= 20) {
        total <- 2^n
        extreme <- 0
        chunk <- 20000
        powers <- 2^(seq_len(n) - 1)
        for (start in seq(0, total - 1, by = chunk)) {
            ids <- start:min(start + chunk - 1, total - 1)
            bits <- outer(ids, powers, bitwAnd) != 0
            permuted <- rowMeans((2 * bits - 1) * matrix(delta, nrow = length(ids), ncol = n, byrow = TRUE))
            extreme <- extreme + sum(abs(permuted) >= observed - 1e-12)
        }
        result <- c(p_value = extreme / total, n_permutations = total, exact = TRUE)
        return(result)
    }
    set.seed(seed)
    total <- 100000
    extreme <- 0
    chunk <- 5000
    for (start in seq(1, total, by = chunk)) {
        current <- min(chunk, total - start + 1)
        signs <- matrix(sample(c(-1, 1), current * n, replace = TRUE), nrow = current)
        permuted <- rowMeans(signs * matrix(delta, nrow = current, ncol = n, byrow = TRUE))
        extreme <- extreme + sum(abs(permuted) >= observed - 1e-12)
    }
    result <- c(p_value = (extreme + 1) / (total + 1), n_permutations = total, exact = FALSE)
    return(result)
}

# With a cohort and cell minimum, return the paired summary and sample scores using frozen scaling and sign flips.
score_fixed_sets <- function(cohort, dataset, threshold) {
    eligible_table <- eligible_idhwt_pairs(cohort$meta, dataset, threshold)
    eligible_counts <- cohort$counts[, eligible_table$unit_id, drop = FALSE]
    eligible_table <- eligible_table[match(colnames(eligible_counts), eligible_table$unit_id), , drop = FALSE]
    eligible_dge <- DGEList(counts = eligible_counts)
    eligible_dge <- calcNormFactors(eligible_dge)
    log_cpm_matrix <- cpm(eligible_dge, log = TRUE, prior.count = 1)
    rows <- list()
    unit_rows <- list()
    for (set_name in names(fixed_sets)) {
        genes <- fixed_sets[[set_name]]
        present <- intersect(genes, rownames(log_cpm_matrix))
        z_matrix <- t(scale(t(log_cpm_matrix[present, , drop = FALSE])))
        z_matrix[!is.finite(z_matrix)] <- 0
        scores <- colMeans(z_matrix)
        unit_table <- eligible_table
        unit_table$threshold <- threshold
        unit_table$signature <- set_name
        unit_table$score <- scores[unit_table$unit_id]
        unit_order <- order(
            unit_table$pair_key,
            factor(unit_table$condition, levels = c("Primary", "Recurrent"))
        )
        unit_table <- unit_table[unit_order, ]
        wide_table <- reshape(
            unit_table[,
                c("pair_key", "condition", "score")],
            idvar     = "pair_key",
            timevar   = "condition",
            direction = "wide"
        )
        delta <- wide_table$score.Recurrent - wide_table$score.Primary
        test <- sign_flip_test(delta, 20260713 + threshold + match(dataset, names(independent_dirs)) * 100)
        rows[[set_name]] <- data.frame(
            dataset           = dataset,
            threshold         = threshold,
            formal_testable   = nrow(wide_table) >= 10,
            signature         = set_name,
            n_pairs           = nrow(wide_table),
            n_input_genes     = length(genes),
            n_present_genes   = length(present),
            coverage          = length(present) / length(genes),
            mean_delta        = mean(delta),
            median_delta      = median(delta),
            positive_pairs    = sum(delta > 0),
            positive_fraction = mean(delta > 0),
            p_value           = unname(test["p_value"]),
            sign_flip_exact   = as.logical(test["exact"]),
            n_permutations    = unname(test["n_permutations"]),
            stringsAsFactors  = FALSE
        )
        unit_table$delta_recurrent_minus_primary <- delta[match(unit_table$pair_key, wide_table$pair_key)]
        unit_rows[[set_name]] <- unit_table
    }
    summary_table_total <- do.call(rbind, rows)
    unit_table_total <- do.call(rbind, unit_rows)
    result <- list(summary = summary_table_total, units = unit_table_total)
    return(result)
}

# With a cohort and cell minimum, return the edgeR, fixed-program GSEA and gene-direction results.
run_pseudobulk <- function(cohort, dataset, threshold) {
    sample_table <- eligible_idhwt_pairs(cohort$meta, dataset, threshold)
    count_matrix <- cohort$counts[, sample_table$unit_id, drop = FALSE]
    sample_table <- sample_table[match(colnames(count_matrix), sample_table$unit_id), , drop = FALSE]
    sample_table$condition <- factor(sample_table$condition, levels = c("Primary", "Recurrent"))
    sample_table$pair_key <- factor(sample_table$pair_key)
    n_pairs <- length(unique(sample_table$pair_key))
    if (n_pairs < 5) stop(dataset, " has fewer than five pairs at threshold ", threshold)
    design_matrix <- model.matrix(~ pair_key + condition, data = sample_table)
    dge <- DGEList(counts = count_matrix)
    keep <- filterByExpr(dge, design = design_matrix)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design_matrix, robust = TRUE)
    fit_model <- glmQLFit(dge, design_matrix, robust = TRUE)
    coef_index <- which(colnames(design_matrix) == "conditionRecurrent")
    if (length(coef_index) != 1) stop("Cannot find recurrent coefficient for ", dataset)
    qlf_model <- glmQLFTest(fit_model, coef = coef_index)
    deg_table <- topTags(qlf_model, n = Inf, sort.by = "PValue")$table
    deg_table$gene <- rownames(deg_table)
    deg_table <- deg_table[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
    deg_table$dataset <- dataset
    deg_table$threshold <- threshold
    deg_table$formal_testable <- n_pairs >= 10
    rank_stat <- sign(deg_table$logFC) * sqrt(pmax(deg_table$F, 0))
    names(rank_stat) <- deg_table$gene
    rank_stat <- sort(rank_stat[is.finite(rank_stat)], decreasing = TRUE)

    targeted_table <- run_fixed_gsea(fixed_sets, rank_stat)
    targeted_table$padj <- NULL
    targeted_table$dataset <- dataset
    targeted_table$threshold <- threshold
    targeted_table$n_pairs <- n_pairs
    targeted_table$formal_testable <- n_pairs >= 10

    gene_rows <- vector("list", length(fixed_sets))
    names(gene_rows) <- names(fixed_sets)
    for (set_name in names(fixed_sets)) {
        genes <- fixed_sets[[set_name]]
        gene_table <- deg_table[match(genes, deg_table$gene), , drop = FALSE]
        gene_rows[[set_name]] <- data.frame(
            dataset   = dataset,
            threshold = threshold,
            signature = set_name,
            gene      = genes,
            tested    = genes %in% deg_table$gene,
            logFC     = gene_table$logFC,
            PValue    = gene_table$PValue,
            FDR       = gene_table$FDR,
            stringsAsFactors = FALSE
        )
    }
    gene_table_total <- do.call(rbind, gene_rows)
    result <- list(deg = deg_table, targeted = targeted_table, genes = gene_table_total)
    return(result)
}

for (dataset in names(cohorts)) {
    cohort <- cohorts[[dataset]]
    for (threshold in c(20, 50)) {
        eligible_before_idh_table <- eligible_pairs(cohort$meta, threshold)
        eligible_table <- eligible_idhwt_pairs(cohort$meta, dataset, threshold)
        n_pairs <- length(unique(eligible_table$pair_key))
        excluded_pairs <- setdiff(
            unique(eligible_before_idh_table$pair_key), unique(eligible_table$pair_key)
        )
        pair_audit[[paste(dataset, threshold)]] <- data.frame(
            dataset              = dataset,
            threshold            = threshold,
            n_pairs_before_idh   = length(unique(eligible_before_idh_table$pair_key)),
            n_pairs_excluded_idh = length(excluded_pairs),
            excluded_pair_keys   = paste(excluded_pairs, collapse = ";"),
            n_pairs              = n_pairs,
            n_units              = nrow(eligible_table),
            n_cells_primary      = sum(eligible_table$n_myeloid_cells[eligible_table$condition == "Primary"]),
            n_cells_recurrent    = sum(eligible_table$n_myeloid_cells[eligible_table$condition == "Recurrent"]),
            formal_testable      = n_pairs >= 10,
            source               = unique(cohort$meta$source),
            count_path           = unique(cohort$meta$count_path),
            metadata_path        = unique(cohort$meta$metadata_path),
            stringsAsFactors     = FALSE
        )
        scores <- score_fixed_sets(cohort, dataset, threshold)
        score_summary[[paste(dataset, threshold)]] <- scores$summary
        score_units[[paste(dataset, threshold)]] <- scores$units
        model <- run_pseudobulk(cohort, dataset, threshold)
        deg_results[[paste(dataset, threshold)]] <- model$deg
        targeted_results[[paste(dataset, threshold)]] <- model$targeted
        gene_results[[paste(dataset, threshold)]] <- model$genes
    }
}

pair_audit_table_total <- do.call(rbind, pair_audit)
score_summary_table_total <- do.call(rbind, score_summary)
score_units_table_total <- do.call(rbind, score_units)
deg_results_table_total <- do.call(rbind, deg_results)
targeted_results_table_total <- do.call(rbind, targeted_results)
gene_results_table_total <- do.call(rbind, gene_results)

expected_pair_table <- data.frame(
    dataset           = c("GSE174554", "GSE174554", "GSE274546", "GSE274546"),
    threshold         = c(20, 50, 20, 50),
    n_pairs           = c(17, 7, 45, 42),
    n_units           = c(34, 14, 90, 84),
    n_cells_primary   = c(4764, 3521, 28494, 27793),
    n_cells_recurrent = c(7725, 4875, 27074, 26834),
    stringsAsFactors  = FALSE
)
observed_pair_table_total <- merge(
    expected_pair_table,
    pair_audit_table_total,
    by       = c("dataset", "threshold"),
    all.x    = TRUE,
    sort     = FALSE,
    suffixes = c("_expected", "_observed")
)
audit_columns <- c(
    "n_pairs", "n_units", "n_cells_primary", "n_cells_recurrent"
)
for (column in audit_columns) {
    expected_column <- paste0(column, "_expected")
    observed_column <- paste0(column, "_observed")
    if (any(observed_pair_table_total[, expected_column] != observed_pair_table_total[, observed_column])) {
        stop("IDH-restricted pair audit failed for ", column, ".")
    }
}
main_count_table <- pair_audit_table_total[pair_audit_table_total$threshold == 20, c("dataset", "n_pairs")]
main_counts <- setNames(main_count_table$n_pairs, main_count_table$dataset)
if (!identical(as.integer(main_counts[c("GSE174554", "GSE274546")]), c(17L, 45L))) {
    stop("Formal IDH-restricted analyses must contain 17 and 45 pairs.")
}

write.csv(
    idh_filter_audit_table,
    file.path(output_dir, "GSE174554_IDH_filter_audit.csv"),
    row.names = FALSE
)
write.csv(pair_audit_table_total, file.path(output_dir, "independent_input_pair_audit.csv"), row.names = FALSE)
write.csv(
    score_summary_table_total,
    file.path(output_dir, "independent_fixed_program_paired_scores.csv"),
    row.names = FALSE
)
write.csv(
    score_units_table_total,
    file.path(output_dir, "independent_fixed_program_scores_by_unit.csv"),
    row.names = FALSE
)
write.csv(deg_results_table_total, file.path(output_dir, "independent_paired_edger_all_genes.csv"), row.names = FALSE)
write.csv(
    targeted_results_table_total,
    file.path(output_dir, "independent_fixed_program_targeted_gsea.csv"),
    row.names = FALSE
)
write.csv(
    gene_results_table_total,
    file.path(output_dir, "independent_fixed_program_gene_direction.csv"),
    row.names = FALSE
)

main_targeted_table <- targeted_results_table_total[targeted_results_table_total$threshold == 20, ]
cat("cohorts=", nrow(main_targeted_table), " paired_units=", sum(pair_audit_table_total$n_units), "\n", sep = "")

# Verify the expected output files from this run.
output_names <- c(
    "fixed_program_provenance.csv", "GSE174554_IDH_filter_audit.csv", "independent_input_pair_audit.csv",
    "independent_fixed_program_paired_scores.csv", "independent_fixed_program_scores_by_unit.csv",
    "independent_paired_edger_all_genes.csv", "independent_fixed_program_targeted_gsea.csv",
    "independent_fixed_program_gene_direction.csv"
)
output_paths <- file.path(output_dir, output_names)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
stopifnot(nrow(pair_audit_table_total) == length(pair_audit))
