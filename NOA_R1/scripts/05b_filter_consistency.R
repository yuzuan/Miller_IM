#!/usr/bin/env Rscript

# Purpose  Refit recurrence camera tests after filtering by condition. Compare with paired-design filtering.
# Inputs   step5_analysis_plan.md and provenance/step1_input_manifest.csv (or step1_input_manifest.csv).
#          The manifest lists each dataset's counts_path, expected_pairs and expected_tested_genes.
#          step1_primary_recurrence/{GSE174554,GSE274546}_paired_voom.rds (metadata, design, voom and signature),
#          camera_results.csv and sample_audit.csv. The manifest points to the prepared pseudobulk counts.
# Outputs  step5_filter_consistency/{camera_results.csv,sample_scores.csv,patient_effects.csv,
#          background_coverage.csv,gene_background_members.csv,checks.csv}.
# Run      Rscript NOA_R1/scripts/05b_filter_consistency.R <scratch_revision_dir>

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_filter_consistency")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir)

# Record the frozen comparability checks before reporting the condition-filtered results. No return value.
check_tables <- list()
record_check <- function(check_name, passed) {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(check = check_name, passed = isTRUE(passed))

    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
        stop(check_name)
    }
}

plan <- readLines(file.path(revision, "step5_analysis_plan.md"), warn = FALSE)
record_check("plan_frozen", any(grepl("2026-10-20 00:55:04", plan, fixed = TRUE)))

# The public manifest keeps the same input basenames and records a portable counts-path prefix.
manifest_path <- file.path(revision, "provenance", "step1_input_manifest.csv")
if (!file.exists(manifest_path)) manifest_path <- file.path(revision, "step1_input_manifest.csv")
manifest_table <- read.csv(manifest_path)
primary_dir <- file.path(revision, "step1_primary_recurrence")
legacy_camera_table <- read.csv(file.path(primary_dir, "camera_results.csv"))
legacy_scores_table <- read.csv(file.path(primary_dir, "sample_audit.csv"))

camera_tables <- list()
score_tables <- list()
effect_tables <- list()
coverage_tables <- list()
universe_tables <- list()

for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, ]
    dataset <- entry_table$dataset
    paired_input <- readRDS(file.path(primary_dir, paste0(dataset, "_paired_voom.rds")))
    metadata_table <- paired_input$metadata

    # Substitute the public-data root only for the placeholder used in the published manifest.
    counts_path <- entry_table$counts_path
    if (startsWith(counts_path, "<public_data_root>/")) {
        project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
        public_data_root <- Sys.getenv(
            "NOA_PUBLIC_DATA_ROOT",
            unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
        )
        counts_path <- file.path(public_data_root, sub("^<public_data_root>/", "", counts_path))
    }
    counts_table <- read.csv(gzfile(counts_path), row.names = 1, check.names = FALSE)
    sample_ids <- as.character(metadata_table$sample_id)
    counts_matrix <- as.matrix(counts_table[, sample_ids, drop = FALSE])
    record_check(
        paste(dataset, "formal_samples"),
        ncol(counts_matrix) == 2 * entry_table$expected_pairs &&
            identical(colnames(counts_matrix), sample_ids) && all(metadata_table$n_cells >= 20L)
    )

    # Hold the original patient-blocked design fixed while changing only the filtering criterion.
    metadata_table$condition <- factor(metadata_table$condition, levels = c("Primary", "Recurrent"))
    metadata_table$pair_key <- factor(metadata_table$patient_id)
    design_matrix <- model.matrix(~ pair_key + condition, metadata_table)
    contrast <- match("conditionRecurrent", colnames(design_matrix))
    record_check(
        paste(dataset, "design_same"),
        identical(unname(design_matrix), unname(paired_input$design)) && length(contrast) == 1L
    )
    all_counts <- DGEList(counts = counts_matrix)
    keep_design <- filterByExpr(all_counts, design = design_matrix)
    keep_condition <- filterByExpr(all_counts, group = metadata_table$condition)
    record_check(
        paste(dataset, "design_background_reproduced"),
        setequal(rownames(counts_matrix)[keep_design], rownames(paired_input$voom)) &&
            sum(keep_design) == entry_table$expected_tested_genes
    )
    condition_genes <- rownames(counts_matrix)[keep_condition]
    target_genes <- intersect(paired_input$signature, condition_genes)
    record_check(paste(dataset, "condition_target_minimum"), length(target_genes) >= 5L)

    # Store the design and condition backgrounds at the gene and dataset levels.
    universe_tables[[i]] <- data.frame(
        dataset          = dataset,
        gene             = rownames(counts_matrix),
        design_tested    = keep_design,
        condition_tested = keep_condition,
        target           = rownames(counts_matrix) %in% paired_input$signature
    )
    coverage_tables[[i]] <- data.frame(
        dataset              = dataset,
        design_background    = sum(keep_design),
        condition_background = sum(keep_condition),
        intersection         = sum(keep_design & keep_condition),
        design_only          = sum(keep_design & !keep_condition),
        condition_only       = sum(keep_condition & !keep_design),
        design_target        = paste(paired_input$present, collapse = ";"),
        condition_target     = paste(target_genes, collapse = ";"),
        condition_target_n   = length(target_genes)
    )

    # Apply the frozen TMM and voom settings. Run estimated correlation, then fixed 0.01.
    condition_counts <- calcNormFactors(all_counts[keep_condition, , keep.lib.sizes = FALSE], method = "TMM")
    voom_data <- voom(condition_counts, design_matrix, normalize.method = "none", adaptive.span = TRUE, plot = FALSE)
    for (mode in c("estimated", "fixed_0.01")) {
        correlation <- if (mode == "estimated") NA_real_ else .01
        camera_table <- camera(
            voom_data, list(target = match(target_genes, rownames(voom_data))),
            design         = design_matrix,
            contrast       = contrast,
            inter.gene.cor = correlation,
            allow.neg.cor  = FALSE,
            use.ranks      = FALSE,
            directional    = TRUE,
            trend.var      = FALSE,
            sort           = FALSE
        )
        legacy_row_table <- legacy_camera_table[
            legacy_camera_table$dataset == dataset & legacy_camera_table$correlation_mode == mode,
        ]
        record_check(paste(dataset, mode, "one_old_row"), nrow(legacy_row_table) == 1L)
        camera_tables[[length(camera_tables) + 1L]] <- data.frame(
            dataset          = dataset,
            filter           = "condition",
            mode             = mode,
            n_pairs          = entry_table$expected_pairs,
            background       = nrow(voom_data),
            n_genes          = camera_table$NGenes[1],
            direction        = camera_table$Direction[1],
            correlation      = if (mode == "estimated") camera_table$Correlation[1] else .01,
            p_value          = camera_table$PValue[1],
            design_p_value   = legacy_row_table$p_value,
            design_direction = legacy_row_table$direction
        )
    }

    # Score the filtered signature from TMM-normalized full-library counts as in the primary analysis.
    log_cpm <- cpm(calcNormFactors(all_counts, method = "TMM"), log = TRUE, prior.count = 1)
    target_log_cpm <- log_cpm[target_genes, , drop = FALSE]
    record_check(paste(dataset, "score_sd"), all(apply(target_log_cpm, 1, sd) > 0))
    score_table <- data.frame(
        dataset    = dataset,
        sample_id  = sample_ids,
        patient_id = as.character(metadata_table$patient_id),
        condition  = as.character(metadata_table$condition),
        score      = unname(colMeans(t(scale(t(target_log_cpm)))))
    )
    score_tables[[i]] <- score_table
    if (identical(target_genes, paired_input$present)) {
        reference_table <- legacy_scores_table[legacy_scores_table$dataset == dataset, ]
        reference_order <- match(score_table$sample_id, reference_table$sample_id)
        differences <- score_table$score - reference_table$score[reference_order]
        record_check(paste(dataset, "identical_members_score_reproduces"), max(abs(differences)) < 1e-10)
    }

    # Each patient contributes one recurrent-minus-primary score difference.
    primary_table <- score_table[score_table$condition == "Primary", ]
    recurrent_table <- score_table[score_table$condition == "Recurrent", ]
    recurrent_order <- match(primary_table$patient_id, recurrent_table$patient_id)
    delta <- recurrent_table$score[recurrent_order] - primary_table$score
    paired_test <- t.test(delta)
    effect_tables[[i]] <- data.frame(
        dataset    = dataset,
        n_pairs    = length(delta),
        n_genes    = length(target_genes),
        positive   = sum(delta > 0),
        mean_delta = mean(delta),
        ci_low     = paired_test$conf.int[1],
        ci_high    = paired_test$conf.int[2],
        paired_t_p = paired_test$p.value
    )
    rm(counts_table, counts_matrix, all_counts, condition_counts, voom_data, log_cpm, target_log_cpm)
    gc(verbose = FALSE)
}

# Assemble results in manifest order, with estimated camera before fixed correlation.
camera_table_total <- do.call(rbind, camera_tables)
score_table_total <- do.call(rbind, score_tables)
effect_table_total <- do.call(rbind, effect_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
universe_table_total <- do.call(rbind, universe_tables)
checks_table_total <- do.call(rbind, check_tables)
write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "sample_scores.csv"), row.names = FALSE, na = "")
write.csv(effect_table_total, file.path(output_dir, "patient_effects.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "background_coverage.csv"), row.names = FALSE, na = "")
write.csv(universe_table_total, file.path(output_dir, "gene_background_members.csv"), row.names = FALSE, na = "")
write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")

cat("FILTER_CONSISTENCY ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
