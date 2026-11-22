#!/usr/bin/env Rscript
# Purpose: Reproduce the six-set paired CAMERA, program scores, and target-excluding random gene-set tests.
# Inputs: provenance/step1_input_manifest.csv (counts_path, expected_pairs and gene totals);
#   step2_locked_gene_sets.csv or the archived step2_extended_specificity/locked_gene_sets_used.csv (set_id, gene);
#   step1_primary_recurrence/<dataset>_paired_voom.rds, camera_results.csv, sample_audit.csv, and legacy Step 38
#   independent_paired_edger_all_genes.csv (gene, logFC, F, PValue) and independent_fixed_program_targeted_gsea.csv.
# Outputs under step2_extended_specificity/: subitem_status.csv, locked_gene_sets_used.csv, camera_results.csv,
#   gene_coverage.csv, patient_scores.csv, patient_effects.csv, random_summary.csv, validation_checks.csv,
#   sessionInfo.txt, and <dataset>_background.csv, _expression_bins.csv and _<scheme>_random_sets.csv.
# Run: Rscript NOA_R1/scripts/02d_extended_specificity.R NOA_R1 after 01a.
# NOA_LEGACY_ROOT supplies Step 38, and NOA_PUBLIC_DATA_ROOT supplies the manifest's pseudobulk inputs.

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(fgsea)
})

# Read the frozen six-set definitions in their original row order.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 02d_extended_specificity.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- file.path(revision, "step2_extended_specificity")
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
if (!file.exists(definition_path)) definition_path <- file.path(output_dir, "locked_gene_sets_used.csv")
definition_table <- read.csv(definition_path, stringsAsFactors = FALSE)
set_ids <- unique(definition_table$set_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Store the reconstruction checks in the original validation CSV, including the first failed check.
checks <- list()

check <- function(dataset, label, condition, detail = "") {
    passed <- isTRUE(condition)
    checks[[length(checks) + 1L]] <<- data.frame(dataset = dataset, check = label, passed = passed, detail = detail)
    validation_table_total <- do.call(rbind, checks)
    write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
    if (!passed) stop(dataset, ": ", label, " ", detail)
}

check("both", "six_source_verified_sets", length(set_ids) == 6L && "Miller_systemic_top20" %in% set_ids)
sets <- setNames(vector("list", length(set_ids)), set_ids)
valid_sets <- TRUE
for (set_id in set_ids) {
    genes <- definition_table$gene[definition_table$set_id == set_id]
    sets[[set_id]] <- genes
    valid_sets <- valid_sets && !anyDuplicated(genes) && !anyNA(genes) && all(nzchar(genes))
}
check("both", "all_sets_unique", valid_sets)

# Record the archived six-set status. Hallmark was completed in the separate seven-set analysis (02e).
status_table <- data.frame(
    subitem = "MSigDB_HALLMARK_INFLAMMATORY_RESPONSE",
    status  = "stopped",
    reason  = "msigdbr human 2025.1 cache absent; no download allowed"
)
write.csv(status_table, file.path(output_dir, "subitem_status.csv"), row.names = FALSE, na = "")
write.csv(definition_table, file.path(output_dir, "locked_gene_sets_used.csv"), row.names = FALSE, na = "")

# The original tested genes and ES come from Step 38 at the 20-cell threshold.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
reference_dir <- file.path(project, "write", "38_independent_cohort_mg_inflammatory_recalculation")
reference_de_table <- read.csv(file.path(reference_dir, "independent_paired_edger_all_genes.csv"))
reference_de_table <- reference_de_table[reference_de_table$threshold == 20L, ]
reference_es_table <- read.csv(file.path(reference_dir, "independent_fixed_program_targeted_gsea.csv"))

# The paired model and target scores are from 01a. Resolve only the manifest path used to open counts.
primary_dir <- file.path(revision, "step1_primary_recurrence")
primary_camera_table <- read.csv(file.path(primary_dir, "camera_results.csv"))
primary_score_table <- read.csv(file.path(primary_dir, "sample_audit.csv"))
manifest_table <- read.csv(file.path(revision, "provenance", "step1_input_manifest.csv"), stringsAsFactors = FALSE)
public_root <- Sys.getenv("NOA_PUBLIC_DATA_ROOT")
if (!nzchar(public_root)) public_root <- Sys.getenv("MILLER_IM_DATA_ROOT")
if (!nzchar(public_root)) public_root <- file.path(project, "data")
public_root <- normalizePath(path.expand(public_root), winslash = "/", mustWork = FALSE)

camera_tables <- list()
coverage_tables <- list()
score_tables <- list()
effect_tables <- list()
random_tables <- list()
for (cohort_index in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[cohort_index, ]
    dataset <- entry_table$dataset
    model <- readRDS(file.path(primary_dir, paste0(dataset, "_paired_voom.rds")))
    genes <- rownames(model$voom)
    sample_table <- model$metadata
    present <- lapply(sets, intersect, y = genes)
    check(
        dataset,
        "formal_pairs_and_background",
        length(unique(sample_table$patient_id)) == entry_table$expected_pairs &&
            length(genes) == entry_table$expected_tested_genes && nrow(sample_table) == 2 * entry_table$expected_pairs
    )
    check(
        dataset,
        "target_nineteen_same_order",
        identical(model$signature, sets$target_top20) && length(present$target_top20) == 19L &&
            identical(present$target_top20, model$present)
    )
    check(dataset, "all_sets_minimum_coverage", all(lengths(present) >= 5L))

    # Retain one coverage row per gene per set and annotate membership in the cohort-specific tested universe.
    for (set_id in set_ids) {
        coverage_table <- definition_table[definition_table$set_id == set_id, , drop = FALSE]
        coverage_table$dataset <- dataset
        coverage_table$tested <- coverage_table$gene %in% genes
        coverage_tables[[length(coverage_tables) + 1L]] <- coverage_table
    }

    # Keep the six comparisons in the same cohort and correlation mode for BH correction.
    indices <- lapply(present, match, table = genes)
    residual_df <- nrow(model$design) - qr(model$design)$rank
    for (mode in c("estimated", "fixed_0.01")) {
        rho <- if (mode == "estimated") NA_real_ else 0.01
        camera_fit_table <- camera(
            model$voom,
            indices,
            design         = model$design,
            contrast       = model$coefficient,
            inter.gene.cor = rho,
            allow.neg.cor  = FALSE,
            use.ranks      = FALSE,
            directional    = TRUE,
            trend.var      = FALSE,
            sort           = FALSE
        )
        camera_table <- data.frame(
            dataset               = dataset,
            set_id                = rownames(camera_fit_table),
            mode                  = mode,
            n_pairs               = entry_table$expected_pairs,
            n_genes               = camera_fit_table$NGenes,
            direction             = camera_fit_table$Direction,
            correlation           = if (mode == "estimated") camera_fit_table$Correlation else rho,
            residual_df           = residual_df,
            camera_reference_df   = if (mode == "estimated") residual_df else nrow(model$voom) - 2L,
            p_value               = camera_fit_table$PValue,
            bh_within_cohort_mode = p.adjust(camera_fit_table$PValue, "BH")
        )
        camera_table$effective_vif <- pmax(1, 1 + (camera_table$n_genes - 1) * camera_table$correlation)

        # The target result must match the single-set 01a p value in both correlation modes.
        check(dataset, paste0("camera_", mode), nrow(camera_table) == length(set_ids) &&
            all(is.finite(camera_table$p_value)) && all(camera_table$p_value >= 0 & camera_table$p_value <= 1))
        primary_table <- primary_camera_table[primary_camera_table$dataset == dataset &
            primary_camera_table$correlation_mode == mode, ]
        target_p <- camera_table$p_value[camera_table$set_id == "target_top20"]
        check(dataset, paste0("top20_reproduces_step1_", mode), nrow(primary_table) == 1L &&
            abs(target_p - primary_table$p_value) < 1e-12)
        camera_tables[[length(camera_tables) + 1L]] <- camera_table
    }

    # Scores use unfiltered counts and TMM logCPM with prior count 1, as in 01a.
    counts_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$counts_path)
    if (!file.exists(counts_path)) stop("Missing pseudobulk counts: ", counts_path)
    counts_table <- read.csv(gzfile(counts_path), row.names = 1L, check.names = FALSE)
    counts_matrix <- as.matrix(counts_table)
    counts_matrix <- counts_matrix[, as.character(sample_table$sample_id), drop = FALSE]
    rm(counts_table)
    check(dataset, "raw_input_sample_match", identical(colnames(counts_matrix), as.character(sample_table$sample_id)) &&
        nrow(counts_matrix) == entry_table$expected_input_genes)
    normalized_counts <- calcNormFactors(DGEList(counts = counts_matrix), method = "TMM")
    log_cpm_matrix <- cpm(normalized_counts, log = TRUE, prior.count = 1)

    for (set_id in set_ids) {
        set_genes <- present[[set_id]]
        set_matrix <- log_cpm_matrix[set_genes, , drop = FALSE]
        gene_sd <- apply(set_matrix, 1, sd)
        check(dataset, paste0(set_id, "_score_valid"), all(is.finite(set_matrix)) &&
            all(is.finite(gene_sd)) && all(gene_sd > 0))
        scores <- colMeans(t(scale(t(set_matrix))))
        score_table <- data.frame(
            dataset    = dataset,
            set_id     = set_id,
            sample_id  = as.character(sample_table$sample_id),
            patient_id = as.character(sample_table$patient_id),
            condition  = as.character(sample_table$condition),
            n_genes    = length(set_genes),
            score      = unname(scores)
        )
        if (set_id == "target_top20") {
            primary_table <- primary_score_table[primary_score_table$dataset == dataset, ]
            primary_scores <- primary_table$score[match(score_table$sample_id, primary_table$sample_id)]
            check(dataset, "top20_score_reproduces_step1", max(abs(score_table$score - primary_scores)) < 1e-10)
        }
        score_tables[[length(score_tables) + 1L]] <- score_table

        # Align recurrent and primary samples by patient before calculating the 95% t interval on the deltas.
        primary_table <- score_table[score_table$condition == "Primary", ]
        recurrent_table <- score_table[score_table$condition == "Recurrent", ]
        check(
            dataset,
            paste0(set_id, "_paired_scores"),
            setequal(primary_table$patient_id, recurrent_table$patient_id) &&
                nrow(primary_table) == entry_table$expected_pairs
        )
        recurrent_rows <- match(primary_table$patient_id, recurrent_table$patient_id)
        delta <- recurrent_table$score[recurrent_rows] - primary_table$score
        paired_test <- t.test(delta, mu = 0, conf.level = .95)
        effect_tables[[length(effect_tables) + 1L]] <- data.frame(
            dataset    = dataset,
            set_id     = set_id,
            n_pairs    = length(delta),
            n_genes    = length(set_genes),
            positive_pairs = sum(delta > 0),
            mean_delta = mean(delta),
            ci_low     = paired_test$conf.int[1],
            ci_high    = paired_test$conf.int[2],
            paired_t_p = paired_test$p.value
        )
    }
    rm(counts_matrix, log_cpm_matrix, normalized_counts)
    gc(verbose = FALSE)

    # The original signed square-root edgeR F statistic clips only numerical negative roundoff.
    cohort_de_table <- reference_de_table[reference_de_table$dataset == dataset, ]
    check(dataset, "same_rank_universe", setequal(cohort_de_table$gene, genes) && !anyDuplicated(cohort_de_table$gene))
    negative <- cohort_de_table$F < 0
    check(dataset, "negative_F_roundoff_only", all(is.finite(cohort_de_table$F)) &&
        all(cohort_de_table$F[negative] > -1e-6) && all(cohort_de_table$PValue[negative] == 1))
    rank_statistic <- setNames(sign(cohort_de_table$logFC) * sqrt(pmax(cohort_de_table$F, 0)), cohort_de_table$gene)
    rank_statistic <- rank_statistic[order(-rank_statistic, names(rank_statistic))]
    target_ranks <- sort(match(present$target_top20, names(rank_statistic)))
    observed_es <- calcGseaStat(rank_statistic, target_ranks, gseaParam = 1, scoreType = "std")
    cohort_es_table <- reference_es_table[reference_es_table$dataset == dataset & reference_es_table$threshold == 20L, ]
    check(dataset, "target_ES_reproduces_original", nrow(cohort_es_table) == 1L &&
        abs(observed_es - cohort_es_table$ES) < 1e-12)

    # Exclude all fixed target members and use the filtered expression rank to create ten matching bins.
    excluded_genes <- sets$target_top20
    pool <- setdiff(genes, excluded_genes)
    check(dataset, "random_pool_excludes_all_twenty", length(pool) == length(genes) - 19L &&
        !any(pool %in% excluded_genes))
    abundance <- aveLogCPM(model$y)
    abundance_order <- order(abundance, genes)
    bins <- integer(length(genes))
    bins[abundance_order] <- pmin(10L, ceiling(seq_along(genes) * 10L / length(genes)))
    names(bins) <- genes
    target_counts <- tabulate(bins[present$target_top20], nbins = 10L)
    eligible_counts <- tabulate(bins[pool], nbins = 10L)
    pools <- split(pool, bins[pool])
    check(dataset, "expression_pool_sufficient", all(eligible_counts >= target_counts),
        paste(target_counts, eligible_counts, collapse = ";"))

    background_table <- data.frame(
        gene              = genes,
        average_logCPM    = abundance,
        expression_decile = bins,
        random_eligible   = genes %in% pool,
        target            = genes %in% present$target_top20
    )
    bins_table <- data.frame(decile = seq_len(10), target_count = target_counts, eligible_count = eligible_counts)
    write.csv(background_table, file.path(output_dir, paste0(dataset, "_background.csv")), row.names = FALSE, na = "")
    write.csv(bins_table, file.path(output_dir, paste0(dataset, "_expression_bins.csv")), row.names = FALSE, na = "")

    # The four original seeds belong to cohort x scheme. Setting each seed once preserves the draw order.
    schemes <- c("size_only", "size_and_expression")
    for (scheme_index in seq_along(schemes)) {
        scheme <- schemes[scheme_index]
        seed <- 20261100L + 10L * cohort_index + scheme_index
        set.seed(seed)
        n_draws <- 10000L
        members <- character(n_draws)
        es <- numeric(n_draws)
        overlap <- integer(n_draws)
        valid <- TRUE
        matching_bins <- if (scheme == "size_and_expression") which(target_counts > 0L) else integer()

        # Sample 19 genes without replacement, traversing nonempty bins in ascending order for expression matching.
        for (draw in seq_len(n_draws)) {
            selected <- character()
            for (bin in matching_bins) {
                selected <- c(selected, sample(pools[[as.character(bin)]], target_counts[bin], replace = FALSE))
            }
            if (scheme == "size_only") selected <- sample(pool, 19L, replace = FALSE)
            matched <- scheme == "size_only" || identical(tabulate(bins[selected], nbins = 10L), target_counts)
            valid <- valid && length(selected) == 19L && !anyDuplicated(selected) &&
                !any(selected %in% excluded_genes) && all(selected %in% pool) && matched
            members[draw] <- paste(sort(selected), collapse = ";")
            overlap[draw] <- sum(selected %in% present$target_top20)
            selected_ranks <- sort(match(selected, names(rank_statistic)))
            es[draw] <- calcGseaStat(rank_statistic, selected_ranks, gseaParam = 1, scoreType = "std")
        }

        # Store the upper and absolute tails using the original plus-one correction.
        check(dataset, paste0("random_valid_", scheme), valid && all(is.finite(es)) && all(overlap == 0L))
        upper_count <- sum(es >= observed_es)
        absolute_count <- sum(abs(es) >= abs(observed_es))
        random_tables[[length(random_tables) + 1L]] <- data.frame(
            dataset     = dataset,
            scheme      = scheme,
            seed        = seed,
            B           = n_draws,
            set_size    = 19L,
            observed_es = observed_es,
            n_es_ge     = upper_count,
            empirical_p = (1 + upper_count) / (n_draws + 1),
            n_abs_ge    = absolute_count,
            abs_empirical_p = (1 + absolute_count) / (n_draws + 1)
        )
        draw_table <- data.frame(draw = seq_len(n_draws), members = members, ES = es, overlap_with_target = overlap)
        draw_path <- file.path(output_dir, paste0(dataset, "_", scheme, "_random_sets.csv"))
        write.csv(draw_table, draw_path, row.names = FALSE, na = "")
    }
}

# Combine results in original cohort, mode and set order. Adjust each random-test tail across the four tests.
camera_table_total <- do.call(rbind, camera_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
score_table_total <- do.call(rbind, score_tables)
effect_table_total <- do.call(rbind, effect_tables)
random_table_total <- do.call(rbind, random_tables)
random_table_total$bh_four_tests <- p.adjust(random_table_total$empirical_p, "BH")
random_table_total$abs_bh_four_tests <- p.adjust(random_table_total$abs_empirical_p, "BH")

write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "gene_coverage.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "patient_scores.csv"), row.names = FALSE, na = "")
write.csv(effect_table_total, file.path(output_dir, "patient_effects.csv"), row.names = FALSE, na = "")
write.csv(random_table_total, file.path(output_dir, "random_summary.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Validation checks passed: ", length(checks), "\n", sep = "")
