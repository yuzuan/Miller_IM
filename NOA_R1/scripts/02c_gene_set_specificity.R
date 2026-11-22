#!/usr/bin/env Rscript
# Purpose: Reproduce the historical five-set CAMERA, preranked fgsea and random-set comparisons.
# Inputs: step2_gene_sets.csv (set_id, gene) or the five matching sets from
#   step2_extended_specificity/locked_gene_sets_used.csv; step1_primary_recurrence/<dataset>_paired_voom.rds and
#   camera_results.csv; legacy Step 38 independent_paired_edger_all_genes.csv and
#   independent_fixed_program_targeted_gsea.csv.
# Outputs under step2_gene_set_specificity/: camera_results.csv, fgsea_results.csv, gene_coverage.csv,
#   set_overlap.csv, random_set_summary.csv, validation_checks.csv, parameters.json, sessionInfo.txt, and
#   <dataset>_negative_F_roundoff_audit.csv, _ranking.csv, _random_background.csv, _expression_matching.csv
#   and _<scheme>_random_sets.csv.
# Run: Rscript NOA_R1/scripts/02c_gene_set_specificity.R NOA_R1 after 01a. NOA_LEGACY_ROOT supplies Step 38.

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(fgsea)
    library(jsonlite)
})

# Prefer the original five-set table. The frozen extended definitions include those exact rows in the same order.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 02c_gene_set_specificity.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- file.path(revision, "step2_gene_set_specificity")
definition_path <- file.path(revision, "step2_gene_sets.csv")
if (file.exists(definition_path)) {
    definition_table <- read.csv(definition_path, stringsAsFactors = FALSE)
} else {
    definition_path <- file.path(revision, "step2_extended_specificity", "locked_gene_sets_used.csv")
    definition_table <- read.csv(definition_path, stringsAsFactors = FALSE)
    definition_table <- definition_table[definition_table$set_id != "Miller_systemic_top20", , drop = FALSE]
}
set_ids <- unique(definition_table$set_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Report each reconstruction check in the existing validation CSV before proceeding.
checks <- list()

check <- function(dataset, label, condition, detail = "") {
    passed <- isTRUE(condition)
    checks[[length(checks) + 1L]] <<- data.frame(dataset = dataset, check = label, passed = passed, detail = detail)
    validation_table_total <- do.call(rbind, checks)
    write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
    if (!passed) stop(dataset, ": ", label, "; ", detail)
}

check(
    "Both",
    "five_fixed_gene_sets",
    length(set_ids) == 5L && all(c("target_top20", "target_top50", "target_top100") %in% set_ids)
)
sets <- setNames(vector("list", length(set_ids)), set_ids)
valid_sets <- TRUE
for (set_id in set_ids) {
    genes <- definition_table$gene[definition_table$set_id == set_id]
    sets[[set_id]] <- genes
    valid_sets <- valid_sets && !anyNA(genes) && all(nzchar(genes)) && !anyDuplicated(genes)
}
check("Both", "unique_genes_in_each_set", valid_sets)
check(
    "Both",
    "nested_target_widths",
    identical(sets$target_top20, head(sets$target_top100, 20)) &&
        identical(sets$target_top50, head(sets$target_top100, 50)) && length(sets$target_top100) == 100
)

# Step 38 supplies the 20-cell tested universe and the target ES for the same formal cohort.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
reference_dir <- file.path(project, "write", "38_independent_cohort_mg_inflammatory_recalculation")
reference_de_table <- read.csv(file.path(reference_dir, "independent_paired_edger_all_genes.csv"))
reference_de_table <- reference_de_table[reference_de_table$threshold == 20L, ]
reference_es_table <- read.csv(file.path(reference_dir, "independent_fixed_program_targeted_gsea.csv"))
primary_dir <- file.path(revision, "step1_primary_recurrence")
primary_camera_table <- read.csv(file.path(primary_dir, "camera_results.csv"))

camera_tables <- list()
fgsea_tables <- list()
coverage_tables <- list()
overlap_tables <- list()
random_tables <- list()
datasets <- c("GSE174554", "GSE274546")
for (cohort_index in seq_along(datasets)) {
    dataset <- datasets[[cohort_index]]
    model <- readRDS(file.path(primary_dir, paste0(dataset, "_paired_voom.rds")))
    genes <- rownames(model$voom)
    cohort_de_table <- reference_de_table[reference_de_table$dataset == dataset, ]
    check(
        dataset,
        "same_tested_universe",
        !anyDuplicated(cohort_de_table$gene) && setequal(genes, cohort_de_table$gene)
    )
    check(dataset, "same_fixed_top20", identical(model$signature, sets$target_top20))
    cohort_de_table <- cohort_de_table[match(genes, cohort_de_table$gene), ]

    # Match the original signed square-root F ranking by clipping numerical negative roundoff to zero.
    negative <- cohort_de_table$F < 0
    check(dataset, "negative_F_only_numerical_roundoff", all(is.finite(cohort_de_table$F)) &&
        all(cohort_de_table$F[negative] > -1e-6) && all(cohort_de_table$PValue[negative] == 1))
    roundoff_table <- cohort_de_table[negative, c("gene", "F", "PValue", "logFC")]
    roundoff_path <- file.path(output_dir, paste0(dataset, "_negative_F_roundoff_audit.csv"))
    write.csv(roundoff_table, roundoff_path, row.names = FALSE, na = "")

    rank_statistic <- setNames(sign(cohort_de_table$logFC) * sqrt(pmax(cohort_de_table$F, 0)), cohort_de_table$gene)
    check(dataset, "finite_signed_sqrtF", all(is.finite(rank_statistic)) && !anyDuplicated(names(rank_statistic)))

    # Break tied statistics by gene name, preserving the original deterministic ranking.
    rank_statistic <- rank_statistic[order(-rank_statistic, names(rank_statistic))]
    rank_table <- data.frame(
        rank = seq_along(rank_statistic),
        gene = names(rank_statistic),
        statistic = unname(rank_statistic)
    )
    write.csv(rank_table, file.path(output_dir, paste0(dataset, "_ranking.csv")), row.names = FALSE, na = "")

    # The original tests accept sets with 5 to 500 measured genes in the filtered cohort universe.
    present <- lapply(sets, intersect, y = genes)
    check(dataset, "minimum_coverage", all(lengths(present) >= 5L) && all(lengths(present) <= 500L))
    for (set_id in set_ids) {
        definition_subset_table <- definition_table[definition_table$set_id == set_id, ]
        coverage_tables[[length(coverage_tables) + 1L]] <- data.frame(
            dataset = dataset,
            definition_subset_table,
            tested = definition_subset_table$gene %in% genes
        )
        overlap <- intersect(present$target_top20, present[[set_id]])
        overlap_tables[[length(overlap_tables) + 1L]] <- data.frame(
            dataset       = dataset,
            set_id        = set_id,
            defined_genes = length(sets[[set_id]]),
            tested_genes  = length(present[[set_id]]),
            overlap_with_target20 = length(overlap),
            overlap_genes = paste(overlap, collapse = ";")
        )
    }

    # Keep the model and CAMERA settings from 01a. Correct over five sets within each cohort and mode.
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
        correlation <- if (mode == "estimated") camera_fit_table$Correlation else rep(rho, nrow(camera_fit_table))
        camera_table <- data.frame(
            dataset             = dataset,
            set_id              = rownames(camera_fit_table),
            correlation_mode    = mode,
            n_genes             = camera_fit_table$NGenes,
            direction           = camera_fit_table$Direction,
            correlation         = correlation,
            effective_vif       = pmax(1, 1 + (camera_fit_table$NGenes - 1) * correlation),
            residual_df         = residual_df,
            camera_reference_df = if (mode == "estimated") residual_df else nrow(model$voom) - 2L,
            p_value             = camera_fit_table$PValue,
            bh_five_sets_within_cohort = p.adjust(camera_fit_table$PValue, "BH")
        )
        check(dataset, paste0("camera_five_sets_", mode), nrow(camera_table) == 5L &&
            all(is.finite(camera_table$p_value)) && all(camera_table$p_value >= 0 & camera_table$p_value <= 1))
        primary_table <- primary_camera_table[primary_camera_table$dataset == dataset &
            primary_camera_table$correlation_mode == mode, ]
        target_p <- camera_table$p_value[camera_table$set_id == "target_top20"]
        check(dataset, paste0("top20_camera_reproduces_step1_", mode), length(target_p) == 1L &&
            abs(target_p - primary_table$p_value) < 1e-12)
        camera_tables[[length(camera_tables) + 1L]] <- camera_table
    }

    # Retain the original serial fgsea settings and separate seed for each cohort.
    fgsea_seed <- 20261019L + cohort_index
    set.seed(fgsea_seed)
    fgsea_table <- fgseaMultilevel(
        present,
        rank_statistic,
        minSize     = 5L,
        maxSize     = 500L,
        eps         = 0,
        gseaParam   = 1,
        scoreType   = "std",
        BPPARAM     = BiocParallel::SerialParam(),
        sampleSize  = 101L,
        nPermSimple = 10000L
    )
    fgsea_table <- as.data.frame(fgsea_table)
    check(dataset, "fgsea_five_valid_results", nrow(fgsea_table) == 5L &&
        all(is.finite(fgsea_table$pval)) && all(is.finite(fgsea_table$NES)) && all(is.finite(fgsea_table$log2err)))
    adjustment_error <- max(abs(fgsea_table$padj - p.adjust(fgsea_table$pval, "BH")))
    check(dataset, "fgsea_BH_matches_declared_family", adjustment_error < 1e-14)
    fgsea_table$leadingEdge <- vapply(fgsea_table$leadingEdge, paste, character(1), collapse = ";")
    colnames(fgsea_table)[colnames(fgsea_table) == "pathway"] <- "set_id"
    colnames(fgsea_table)[colnames(fgsea_table) == "padj"] <- "bh_five_sets_within_cohort"
    fgsea_table$dataset <- dataset
    fgsea_table$seed <- fgsea_seed
    fgsea_table$ranking <- "sign(logFC)*sqrt(pmax(edgeR_QL_F,0)); original Step38 definition"
    fgsea_tables[[dataset]] <- fgsea_table

    # Check the target ES against both the fitted fgsea result and the prior Step 38 result.
    target_ranks <- sort(match(present$target_top20, names(rank_statistic)))
    observed_es <- calcGseaStat(rank_statistic, target_ranks, gseaParam = 1, scoreType = "std")
    fgsea_es <- fgsea_table$ES[fgsea_table$set_id == "target_top20"]
    check(dataset, "observed_ES_matches_fgsea", abs(observed_es - fgsea_es) < 1e-12)
    cohort_es_table <- reference_es_table[reference_es_table$dataset == dataset & reference_es_table$threshold == 20L, ]
    check(dataset, "top20_ES_reproduces_original", nrow(cohort_es_table) == 1L &&
        abs(cohort_es_table$ES - observed_es) < 1e-12)

    # Match on filtered mean logCPM rank-deciles, with gene names breaking abundance ties.
    abundance <- aveLogCPM(model$y)
    abundance_order <- order(abundance, genes)
    bins <- integer(length(genes))
    bins[abundance_order] <- pmin(10L, ceiling(seq_along(genes) * 10L / length(genes)))
    names(bins) <- genes
    target_counts <- tabulate(bins[present$target_top20], nbins = 10L)
    pools <- split(genes, bins)
    background_table <- data.frame(
        gene              = genes,
        average_logCPM    = abundance,
        expression_decile = bins,
        target_top20      = genes %in% present$target_top20
    )
    matching_table <- data.frame(
        expression_decile = seq_len(10),
        target_count      = target_counts,
        background_count  = tabulate(bins, nbins = 10)
    )
    background_path <- file.path(output_dir, paste0(dataset, "_random_background.csv"))
    matching_path <- file.path(output_dir, paste0(dataset, "_expression_matching.csv"))
    write.csv(background_table, background_path, row.names = FALSE, na = "")
    write.csv(matching_table, matching_path, row.names = FALSE, na = "")

    # The five-set analysis samples from the full tested universe, including target genes.
    n_draws <- 10000L
    schemes <- c("size_only", "size_and_expression")
    for (scheme_index in seq_along(schemes)) {
        scheme <- schemes[[scheme_index]]
        seed <- 20261100L + 10L * cohort_index + scheme_index
        set.seed(seed)
        memberships <- character(n_draws)
        es <- numeric(n_draws)
        overlap_counts <- integer(n_draws)
        size_ok <- TRUE
        unique_ok <- TRUE
        background_ok <- TRUE
        matched_ok <- TRUE
        matching_bins <- if (scheme == "size_and_expression") which(target_counts > 0L) else integer()

        # Sample by ascending bin when matching mean expression. The seed and sequence of draws are unchanged.
        for (draw in seq_len(n_draws)) {
            selected <- character()
            for (bin in matching_bins) {
                selected <- c(selected, sample(pools[[as.character(bin)]], target_counts[bin], replace = FALSE))
            }
            if (scheme == "size_only") selected <- sample(genes, length(present$target_top20), replace = FALSE)
            size_ok <- size_ok && length(selected) == length(present$target_top20)
            unique_ok <- unique_ok && !anyDuplicated(selected)
            background_ok <- background_ok && all(selected %in% genes)
            matched_ok <- matched_ok && (scheme == "size_only" ||
                identical(tabulate(bins[selected], nbins = 10L), target_counts))
            memberships[[draw]] <- paste(sort(selected), collapse = ";")
            overlap_counts[[draw]] <- length(intersect(selected, present$target_top20))
            selected_ranks <- sort(match(selected, names(rank_statistic)))
            es[[draw]] <- calcGseaStat(rank_statistic, selected_ranks, gseaParam = 1, scoreType = "std")
        }

        # Save the original empirical p values, with a plus-one correction for both tails, and all draws.
        check(dataset, paste0(scheme, "_all_replicates_valid"), size_ok && unique_ok && background_ok &&
            matched_ok && all(is.finite(es)))
        upper_count <- sum(es >= observed_es)
        absolute_count <- sum(abs(es) >= abs(observed_es))
        random_tables[[length(random_tables) + 1L]] <- data.frame(
            dataset                = dataset,
            scheme                 = scheme,
            n_random_sets          = n_draws,
            set_size               = length(present$target_top20),
            observed_ES            = observed_es,
            null_mean_ES           = mean(es),
            null_median_ES         = median(es),
            null_q025              = unname(quantile(es, .025)),
            null_q975              = unname(quantile(es, .975)),
            positive_extreme_count = upper_count,
            empirical_positive_p   = (1 + upper_count) / (n_draws + 1),
            absolute_extreme_count = absolute_count,
            empirical_absolute_p   = (1 + absolute_count) / (n_draws + 1),
            seed                   = seed
        )
        draw_table <- data.frame(
            replicate = seq_len(n_draws),
            ES    = es,
            overlap_with_target20 = overlap_counts,
            genes = memberships
        )
        draw_path <- file.path(output_dir, paste0(dataset, "_", scheme, "_random_sets.csv"))
        write.csv(draw_table, draw_path, row.names = FALSE, na = "")
    }
}

# Adjust each random-set tail across its original four cohort-scheme tests and save the result tables.
camera_table_total <- do.call(rbind, camera_tables)
fgsea_table_total <- do.call(rbind, fgsea_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
overlap_table_total <- do.call(rbind, overlap_tables)
random_table_total <- do.call(rbind, random_tables)
random_table_total$bh_four_random_tests <- p.adjust(random_table_total$empirical_positive_p, "BH")
random_table_total$bh_four_absolute_tests <- p.adjust(random_table_total$empirical_absolute_p, "BH")
write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(fgsea_table_total, file.path(output_dir, "fgsea_results.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "gene_coverage.csv"), row.names = FALSE, na = "")
write.csv(overlap_table_total, file.path(output_dir, "set_overlap.csv"), row.names = FALSE, na = "")
write.csv(random_table_total, file.path(output_dir, "random_set_summary.csv"), row.names = FALSE, na = "")

# Record the original plan reference, thresholds, seeds, and software versions for traceability.
package_versions <- list()
for (package in c("edgeR", "limma", "fgsea", "BiocParallel")) {
    package_versions[[package]] <- as.character(packageVersion(package))
}
parameters <- list(
    date   = "2026-10-19",
    plan   = "step2_analysis_plan.md",
    sets   = set_ids,
    comparisons = "five sets per cohort per test method/settings; BH=5",
    camera = list(
        primary_correlation = "estimated (NA)",
        sensitivity_correlation = 0.01,
        paired_design       = "~ pair_key + condition",
        observation_weights = TRUE,
        allow_neg_cor       = FALSE,
        use_ranks           = FALSE,
        directional         = TRUE,
        two_sided           = TRUE
    ),
    fgsea = list(
        rank_statistic = "sign(logFC)*sqrt(pmax(edgeR_QL_F,0)); original Step38 definition",
        gseaParam      = 1,
        scoreType      = "std",
        eps            = 0,
        sampleSize     = 101,
        nPermSimple    = 10000,
        parallel       = "SerialParam",
        seeds          = 20261019L + seq_along(datasets)
    ),
    random = list(
        B                      = n_draws,
        statistic              = "same weighted GSEA ES, std, exponent=1",
        null                   = "gene membership resampling, not sample/patient permutation",
        primary_empirical_test = "positive/upper-tail: (1+sum(nullES>=observedES))/(B+1)",
        sensitivity_empirical_test = "absolute ES upper tail; not identical to fgsea sign-conditional P",
        multiple_testing       = "BH across two cohorts x two schemes, 4 tests per tail",
        background             = "cohort-specific filterByExpr universe, includes target members",
        expression_matching    = "average logCPM rank-deciles, tie order by gene, exact target decile counts"
    ),
    package_versions = package_versions
)
write_json(parameters, file.path(output_dir, "parameters.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Validation checks passed: ", length(checks), "\n", sep = "")
