#!/usr/bin/env Rscript
# Purpose: Fit the frozen paired state-expression models and descriptive within-state target scores.
# Inputs: step2_locked_gene_sets.csv, or step2_gene_sets.csv if rebuilt from the source workbooks.
# Read step3_state_analysis/{pseudobulk_manifest,expression_pair_coverage}.csv and the listed matrices/metadata.
# Outputs in step3_state_analysis: state_program_camera.csv, state_all_gene_DE.csv.gz, target_gene_DE.csv,
# state_score_effects.csv, patient_state_score_changes.csv, state_program_gene_coverage.csv,
# expression_model_audit.csv, expression_sample_audit.csv, expression_checks.csv,
# expression_parameters.json and expression_sessionInfo.txt.
# Outputs in step3_state_analysis/models: <model>_design.csv and, for a fresh fit, <model>.rds and
# <model>_{BCV,voom}.pdf. An explicit resume reads and verifies existing model RDS files.
# Run: Rscript NOA_R1/scripts/03b_state_expression.R NOA_R1 [--resume-cache]

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) stop("Usage: 03b_state_expression.R <revision_dir> [--resume-cache]")
resume_cache <- length(args) == 2L
if (resume_cache && args[[2]] != "--resume-cache") stop("Unknown second argument")
revision <- normalizePath(args[[1]], mustWork = TRUE)

# Use the written count manifests and the frozen author-ranked gene order.
output_dir <- file.path(revision, "step3_state_analysis")
manifest_table <- read.csv(file.path(output_dir, "pseudobulk_manifest.csv"))
coverage_table <- read.csv(file.path(output_dir, "expression_pair_coverage.csv"))
definitions_path <- file.path(revision, "step2_gene_sets.csv")
if (!file.exists(definitions_path)) definitions_path <- file.path(revision, "step2_locked_gene_sets.csv")
definitions_table <- read.csv(definitions_path)
target <- definitions_table$gene[definitions_table$set_id == "target_top20"]
model_dir <- file.path(output_dir, "models")
dir.create(model_dir, showWarnings = FALSE)

program_rows <- list()
gene_tables <- list()
model_rows <- list()
score_rows <- list()
delta_tables <- list()
coverage_tables <- list()
sample_tables <- list()
check_rows <- list()

# Record a frozen alignment/model check immediately so a failing check is also retained. Returns nothing.
check <- function(dataset, label, value, detail = "") {
    passed <- isTRUE(value)
    check_rows[[length(check_rows) + 1L]] <<- data.frame(
        dataset = dataset,
        check   = label,
        passed  = passed,
        detail  = detail
    )
    checks_table <- do.call(rbind, check_rows)
    write.csv(checks_table, file.path(output_dir, "expression_checks.csv"), row.names = FALSE)
    if (!passed) stop(dataset, ". ", label, ". ", detail)
}

# The archived manifest uses revision-relative paths. Original local manifests used absolute paths.
for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, , drop = FALSE]
    dataset <- entry_table$dataset
    state <- entry_table$state
    state_id <- paste(dataset, state, sep = "_")
    input_paths <- c(counts = entry_table$counts_path, metadata = entry_table$metadata_path)
    relative <- !grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", input_paths)
    input_paths[relative] <- file.path(revision, input_paths[relative])
    count_table <- read.csv(gzfile(input_paths[["counts"]]), row.names = 1, check.names = FALSE)
    raw_matrix <- as.matrix(count_table)
    metadata_table <- read.csv(input_paths[["metadata"]], stringsAsFactors = FALSE)

    aligned <- identical(colnames(raw_matrix), metadata_table$sample_id) &&
        !anyDuplicated(rownames(raw_matrix)) &&
        all(is.finite(raw_matrix)) &&
        all(raw_matrix >= 0) &&
        all(raw_matrix == round(raw_matrix)) &&
        all(colSums(raw_matrix) == metadata_table$library_size)
    check(dataset, paste0(state, "_raw_counts_alignment"), aligned)
    raw_target <- intersect(target, rownames(raw_matrix))
    check(dataset, paste0(state, "_raw_target_coverage"), identical(raw_target, target[target != "AC253572.2"]))

    # Both endpoints need at least 20 cells in the primary analysis or 10 in the sensitivity analysis.
    for (threshold in c(20L, 10L)) {
        selected <- coverage_table$dataset == dataset &
            coverage_table$state == state & coverage_table$threshold == threshold
        coverage_row_table <- coverage_table[selected, , drop = FALSE]
        patients <- character()
        if (!is.na(coverage_row_table$eligible_patients) && nzchar(coverage_row_table$eligible_patients)) {
            patients <- strsplit(coverage_row_table$eligible_patients, ";", fixed = TRUE)[[1]]
        }
        pair_count <- length(patients)
        check(
            dataset, paste0(state, "_", threshold, "_coverage"),
            nrow(coverage_row_table) == 1L && pair_count == coverage_row_table$n_pairs
        )

        # The prespecified minimum is six complete pairs. Retain unfitted rows in the CAMERA family.
        if (pair_count < 6L) {
            for (mode in c("estimated", "fixed_0.01")) {
                program_rows[[length(program_rows) + 1L]] <- data.frame(
                    dataset   = dataset,
                    state     = state,
                    level     = entry_table$level,
                    threshold = threshold,
                    n_pairs   = pair_count,
                    n_genes   = NA_integer_,
                    tested_universe  = NA_integer_,
                    status           = "insufficient_paired_patients",
                    correlation_mode = mode,
                    correlation      = NA_real_,
                    effective_vif    = NA_real_,
                    direction        = NA_character_,
                    residual_df      = NA_integer_,
                    camera_reference_df = NA_integer_,
                    p_value          = NA_real_
                )
            }
            next
        }

        # Fit patients as blocks and keep Primary as the reference endpoint.
        model_table <- metadata_table[metadata_table$patient_id %in% patients, , drop = FALSE]
        count_matrix <- raw_matrix[, model_table$sample_id, drop = FALSE]
        model_table$patient <- factor(model_table$patient_id)
        model_table$condition <- factor(model_table$condition, levels = c("Primary", "Recurrent"))
        design <- model.matrix(~ patient + condition, data = model_table)
        coefficient <- which(colnames(design) == "conditionRecurrent")
        design_rank <- qr(design)$rank
        residual_df <- nrow(model_table) - design_rank
        model_id <- paste(state_id, paste0("min", threshold), sep = "_")

        complete_pairs <- nrow(model_table) == 2L * pair_count &&
            all(table(model_table$patient, model_table$condition) == 1L) &&
            all(model_table$n_cells >= threshold) && design_rank == ncol(design) &&
            residual_df == pair_count - 1L && length(coefficient) == 1L
        check(dataset, paste0(model_id, "_complete_paired_design"), complete_pairs)

        # Use the frozen condition-based filter and TMM factors before robust edgeR fitting.
        counts <- DGEList(counts = count_matrix)
        keep <- filterByExpr(counts, group = model_table$condition)
        counts <- counts[keep, , keep.lib.sizes = FALSE]
        counts <- calcNormFactors(counts, method = "TMM")
        valid_filter <- nrow(counts) > 100L && all(counts$samples$lib.size > 0) &&
            all(is.finite(counts$samples$norm.factors)) && all(counts$samples$norm.factors > 0)
        check(dataset, paste0(model_id, "_filtered_matrix_valid"), valid_filter)

        # An explicit resume must reproduce the current design, filtered counts and normalization.
        cache_path <- file.path(model_dir, paste0(model_id, ".rds"))
        if (resume_cache) {
            if (!file.exists(cache_path)) stop("Explicit cache resume requested but model absent. ", model_id)
            cached <- readRDS(cache_path)
            same_design <- identical(cached$design, design) &&
                identical(as.character(cached$metadata$sample_id), model_table$sample_id) &&
                identical(cached$coefficient, coefficient)
            check(dataset, paste0(model_id, "_cache_design_and_patients"), same_design)
            same_counts <- identical(rownames(cached$y), rownames(counts)) &&
                identical(colnames(cached$y), colnames(counts)) &&
                all(cached$y$counts == counts$counts) &&
                max(abs(cached$y$samples$norm.factors - counts$samples$norm.factors)) < 1e-12 &&
                all(cached$y$samples$lib.size == counts$samples$lib.size)
            check(dataset, paste0(model_id, "_cache_counts_filter_TMM"), same_counts)
            counts <- cached$y
            test <- cached$edgeR_test
        } else {
            counts <- estimateDisp(counts, design, robust = TRUE)
            fit <- glmQLFit(counts, design, robust = TRUE, legacy = FALSE)
            test <- glmQLFTest(fit, coef = coefficient)
        }

        # Preserve the package's gene order and the original within-model BH column.
        de_table <- topTags(test, n = Inf, sort.by = "none")$table
        de_table$gene <- rownames(de_table)
        de_table$dataset <- dataset
        de_table$state <- state
        de_table$level <- entry_table$level
        de_table$threshold <- threshold
        de_table$n_pairs <- pair_count
        de_table$model_id <- model_id
        colnames(de_table)[colnames(de_table) == "FDR"] <- "bh_within_model_genome"
        valid_de <- nrow(de_table) == nrow(counts) && all(is.finite(de_table$PValue)) &&
            all(de_table$PValue >= 0 & de_table$PValue <= 1) && all(is.finite(de_table$logFC))
        check(dataset, paste0(model_id, "_edgeR_finite_results"), valid_de)
        gene_tables[[model_id]] <- de_table

        # plotBCV and voom draw directly on a base PDF device. Save the voom plot data in the model.
        if (resume_cache) {
            weighted <- cached$voom
        } else {
            bcv_path <- file.path(model_dir, paste0(model_id, "_BCV.pdf"))
            pdf(bcv_path, width = 6, height = 4)
                plotBCV(counts)
            dev.off()

            voom_path <- file.path(model_dir, paste0(model_id, "_voom.pdf"))
            pdf(voom_path, width = 6, height = 4)
                weighted <- voom(
                    counts, design,
                    normalize.method = "none",
                    adaptive.span = TRUE,
                    plot          = TRUE,
                    save.plot     = TRUE
                )
            dev.off()
            if (!all(file.info(c(bcv_path, voom_path))$size > 0, na.rm = FALSE)) {
                stop("Missing or empty expression diagnostic for ", model_id)
            }
        }
        valid_weights <- all(is.finite(weighted$E)) && all(is.finite(weighted$weights)) && all(weighted$weights > 0)
        check(dataset, paste0(model_id, "_voom_weights_valid"), valid_weights)

        # CAMERA uses all tested genes as the comparison universe and requires at least 16 target genes.
        present <- intersect(target, rownames(counts))
        coverage_tables[[model_id]] <- data.frame(
            dataset     = dataset,
            state       = state,
            threshold   = threshold,
            gene        = target,
            raw_present = target %in% rownames(count_matrix),
            tested      = target %in% rownames(counts)
        )
        status <- if (length(present) >= 16L) "analyzed" else "program_coverage_below_16_of_20"
        for (mode in c("estimated", "fixed_0.01")) {
            p_value <- NA_real_
            correlation <- NA_real_
            effective_vif <- NA_real_
            direction <- NA_character_
            rho <- if (mode == "estimated") NA_real_ else .01
            if (status == "analyzed") {
                target_index <- list(target = match(present, rownames(weighted)))
                camera_table <- camera(
                    weighted,
                    index = target_index,
                    design   = design,
                    contrast = coefficient,
                    inter.gene.cor = rho,
                    allow.neg.cor  = FALSE,
                    use.ranks      = FALSE,
                    directional    = TRUE,
                    trend.var      = FALSE,
                    sort           = FALSE
                )
                p_value <- camera_table$PValue
                direction <- camera_table$Direction
                correlations <- list(estimated = camera_table$Correlation, fixed_0.01 = rho)
                correlation <- correlations[[mode]]
                effective_vif <- max(1, 1 + (length(present) - 1L) * correlation)
                finite_result <- is.finite(p_value) && p_value >= 0 && p_value <= 1 && is.finite(correlation)
                check(dataset, paste0(model_id, "_camera_", mode), finite_result)
            }

            reference_df <- nrow(counts) - 2L
            if (mode == "estimated") reference_df <- min(residual_df, reference_df)
            program_rows[[length(program_rows) + 1L]] <- data.frame(
                dataset      = dataset,
                state        = state,
                level        = entry_table$level,
                threshold    = threshold,
                n_pairs      = pair_count,
                n_genes      = length(present),
                tested_universe  = nrow(counts),
                status           = status,
                correlation_mode = mode,
                correlation      = correlation,
                effective_vif    = effective_vif,
                direction        = direction,
                residual_df      = residual_df,
                camera_reference_df = reference_df,
                p_value          = p_value
            )
        }

        # The descriptive score is the mean of 19 measurable genes on unfiltered TMM logCPM.
        score_counts <- DGEList(counts = count_matrix)
        score_counts <- calcNormFactors(score_counts, method = "TMM")
        raw_present <- intersect(target, rownames(count_matrix))
        log_cpm <- cpm(score_counts, log = TRUE, prior.count = 1)
        scores <- colMeans(log_cpm[raw_present, , drop = FALSE])
        delta_table <- data.frame(
            dataset    = dataset,
            state      = state,
            level      = entry_table$level,
            threshold  = threshold,
            patient_id = patients
        )
        primary <- model_table$condition == "Primary"
        recurrent <- model_table$condition == "Recurrent"
        delta_table$primary_score <- scores[primary][match(patients, model_table$patient_id[primary])]
        delta_table$recurrent_score <- scores[recurrent][match(patients, model_table$patient_id[recurrent])]
        delta_table$delta <- delta_table$recurrent_score - delta_table$primary_score
        check(
            dataset, paste0(model_id, "_score_finite"),
            length(scores) == nrow(model_table) && all(is.finite(delta_table$delta))
        )

        if (sd(delta_table$delta) > 0) {
            paired_test <- t.test(delta_table$delta)
            interval <- unname(paired_test$conf.int)
            t_p <- paired_test$p.value
        } else {
            interval <- c(NA_real_, NA_real_)
            t_p <- NA_real_
        }
        score_rows[[model_id]] <- data.frame(
            dataset       = dataset,
            state         = state,
            level         = entry_table$level,
            threshold     = threshold,
            n_pairs       = pair_count,
            n_score_genes = length(raw_present),
            mean_delta    = mean(delta_table$delta),
            ci_low        = interval[1],
            ci_high       = interval[2],
            t_p_value     = t_p,
            positive_pairs   = sum(delta_table$delta > 0),
            score_definition = paste0(
                "mean19 TMM logCPM prior.count1 on unfiltered state matrix",
                "; not pan-myeloid z-score"
            ),
            program_status   = status
        )
        delta_tables[[model_id]] <- delta_table

        # Keep one sample audit row per fitted patient/endpoint and one summary per fitted model.
        model_table$included_in_model <- model_id
        model_table$tmm_factor_DE <- counts$samples$norm.factors
        model_table$filtered_library_size <- counts$samples$lib.size
        model_table$effective_library_size <- counts$samples$lib.size * counts$samples$norm.factors
        model_table$tmm_factor_score <- score_counts$samples$norm.factors
        audit_columns <- c(
            "dataset", "state", "patient_id", "condition", "sample_id", "n_cells", "library_size",
            "included_in_model", "tmm_factor_DE", "filtered_library_size", "effective_library_size", "tmm_factor_score"
        )
        sample_tables[[model_id]] <- model_table[, audit_columns, drop = FALSE]
        model_rows[[model_id]] <- data.frame(
            dataset       = dataset,
            state         = state,
            level         = entry_table$level,
            threshold     = threshold,
            n_pairs       = pair_count,
            primary_cells = sum(model_table$n_cells[primary]),
            recurrent_cells = sum(model_table$n_cells[recurrent]),
            input_genes     = nrow(count_matrix),
            tested_genes    = nrow(counts),
            n_program_genes = length(present),
            residual_df     = residual_df,
            design_rank     = design_rank,
            voom_span       = weighted$span,
            program_status  = status,
            filter_definition = "filterByExpr(group=condition)"
        )

        # Retain the design CSV and the public fields used by explicit cache recovery.
        design_table <- data.frame(sample_id = model_table$sample_id, design, check.names = FALSE)
        write.csv(design_table, file.path(model_dir, paste0(model_id, "_design.csv")), row.names = FALSE)
        if (!resume_cache) {
            cache <- list(
                y        = counts,
                voom     = weighted,
                metadata = model_table,
                design   = design,
                coefficient = coefficient,
                present     = present,
                target      = target,
                edgeR_test  = test
            )
            saveRDS(cache, cache_path)
        }
    }
}

# BY and BH operate on each complete eleven-state family. Missing tests enter as P=1 for adjustment.
program_table_total <- do.call(rbind, program_rows)
program_table_total$by_eleven_primary <- NA_real_
program_table_total$bh_eleven_sensitivity <- NA_real_
program_table_total$family_size <- NA_integer_
for (dataset in unique(program_table_total$dataset)) {
    for (threshold in c(20L, 10L)) {
        for (mode in c("estimated", "fixed_0.01")) {
            selected <- program_table_total$dataset == dataset &
                program_table_total$threshold == threshold & program_table_total$correlation_mode == mode
            indices <- which(selected)
            check(dataset, paste0("eleven_hypothesis_family_", threshold, "_", mode), length(indices) == 11L)
            p_values <- program_table_total$p_value[indices]
            missing <- is.na(p_values)
            p_values[missing] <- 1
            by <- p.adjust(p_values, "BY")
            bh <- p.adjust(p_values, "BH")
            by[missing] <- NA_real_
            bh[missing] <- NA_real_
            program_table_total$by_eleven_primary[indices] <- by
            program_table_total$bh_eleven_sensitivity[indices] <- bh
            program_table_total$family_size[indices] <- 11L
        }
    }
}

# Separately adjust across all fitted state-by-gene hypotheses within each cohort and threshold.
gene_table_total <- do.call(rbind, gene_tables)
gene_table_total$bh_all_fitted_state_gene_tests <- NA_real_
gene_table_total$pooled_family_size <- NA_integer_
for (dataset in unique(gene_table_total$dataset)) {
    for (threshold in c(20L, 10L)) {
        indices <- which(gene_table_total$dataset == dataset & gene_table_total$threshold == threshold)
        gene_table_total$bh_all_fitted_state_gene_tests[indices] <- p.adjust(gene_table_total$PValue[indices], "BH")
        gene_table_total$pooled_family_size[indices] <- length(indices)
    }
}

# Combine the original model, sample and gene rows once, then retain the published file names and schemas.
model_table_total <- do.call(rbind, model_rows)
score_table_total <- do.call(rbind, score_rows)
delta_table_total <- do.call(rbind, delta_tables)
sample_table_total <- do.call(rbind, sample_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
target_de_table <- gene_table_total[gene_table_total$gene %in% target, , drop = FALSE]

write.csv(program_table_total, file.path(output_dir, "state_program_camera.csv"), row.names = FALSE, na = "")
gene_connection <- gzfile(file.path(output_dir, "state_all_gene_DE.csv.gz"), "wt")
write.csv(gene_table_total, gene_connection, row.names = FALSE, na = "")
close(gene_connection)
write.csv(model_table_total, file.path(output_dir, "expression_model_audit.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "state_score_effects.csv"), row.names = FALSE, na = "")
write.csv(delta_table_total, file.path(output_dir, "patient_state_score_changes.csv"), row.names = FALSE, na = "")
write.csv(sample_table_total, file.path(output_dir, "expression_sample_audit.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "state_program_gene_coverage.csv"), row.names = FALSE, na = "")
write.csv(target_de_table, file.path(output_dir, "target_gene_DE.csv"), row.names = FALSE, na = "")

# These existing provenance files record the exact package versions and analysis options.
capture.output(sessionInfo(), file = file.path(output_dir, "expression_sessionInfo.txt"))
package_versions <- list()
for (package in c("edgeR", "limma", "statmod")) package_versions[[package]] <- as.character(packageVersion(package))
parameters <- list(
    date = "2026-10-19",
    recovered_from_explicit_cache = resume_cache,
    thresholds           = c(20, 10),
    minimum_pairs        = 6,
    minimum_program_genes = 16,
    filter               = "filterByExpr(group=condition)",
    design               = "~patient+condition",
    normalization        = "TMM",
    robust_dispersion    = TRUE,
    robust_QL            = TRUE,
    legacy               = FALSE,
    camera               = list(primary = "estimated correlation NA", sensitivity = .01, two_sided = TRUE),
    program_adjustment   = "BY11 primary; BH11 sensitivity; missing tests P=1 for correction, results NA",
    gene_adjustment      = "BH within model genome; separately BH across fitted gene-state tests per cohort/threshold",
    score                = paste0(
        "19 raw measurable genes; TMM on unfiltered model matrix; ",
        "mean logCPM prior.count1; paired t CI/P"
    ),
    packages             = package_versions
)
write_json(parameters, file.path(output_dir, "expression_parameters.json"), auto_unbox = TRUE, pretty = TRUE)
cat(nrow(model_table_total), " fitted models and ", length(check_rows), " passing expression checks\n", sep = "")
