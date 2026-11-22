#!/usr/bin/env Rscript

# Purpose  Independently verify full-model diagnostics and the fixed-protein sign-flip p values.
# Inputs   step5_pd1_stability_final/{formal_score_samples.csv,influence.csv,vif.csv},
#          and step5_protein_common/{patient_scores.csv,effects.csv}.
#          Use sample, cohort, exposure and covariates for models, and use analysis, mode and delta for sign flips.
# Outputs  step5_closeout_v2/{independent_cook_dfbeta.csv,independent_vif.csv,independent_diagnostics_checks.csv,
#          independent_signflip.csv,independent_signflip_checks.csv,independent_signflip_engine.csv}.
# Run      Rscript NOA_R1/scripts/05h_check_pd1_covariate_models.R <scratch_revision_dir>
# Use a scratch closeout directory without the independent diagnostic outputs.

suppressPackageStartupMessages(library(dqrng))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_closeout_v2")
if (!dir.exists(output_dir)) stop("Run the covariate closeout in this revision directory first.")

# Read formal patient scores and prior diagnostics. The matrix calculations use only the patient scores.
stability_dir <- file.path(revision, "step5_pd1_stability_final")
formal_table <- read.csv(file.path(stability_dir, "formal_score_samples.csv"))
prior_influence_table <- read.csv(file.path(stability_dir, "influence.csv"))
prior_vif_table <- read.csv(file.path(stability_dir, "vif.csv"))
influence_tables <- list()
vif_tables <- list()

# Write each review table exclusively. Missing values use empty CSV fields. No return value.
write_new_table <- function(input_table, filename) {
    output_path <- file.path(output_dir, filename)
    if (file.exists(output_path)) stop("Output file already exists: ", output_path)
    write.csv(input_table, output_path, row.names = FALSE, na = "")
}

# Preserve the first failed check in the diagnostic file whenever it does not already exist. No return value.
check_tables <- list()
record_check <- function(check_name, passed, detail = "") {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(
        check  = check_name,
        passed = isTRUE(passed),
        detail = as.character(detail)
    )
    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write_new_table(checks_table_total, "independent_diagnostics_checks.csv")
        stop(check_name)
    }
}

for (dataset in c("GSE121810", "GSE154795")) {
    cohort_table <- formal_table[formal_table$cohort == dataset, ]
    if (dataset == "GSE121810") {
        cohort_table$group <- factor(cohort_table$group, levels = c("Adj", "Neo"))
        formula <- raw20_score ~ myeloid_abundance + group
        term <- "groupNeo"
    } else {
        cohort_table$condition <- factor(cohort_table$condition, levels = c("GBM.rec", "GBM.PD1"))
        for (field in c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave")) {
            cohort_table[, field] <- factor(cohort_table[, field])
        }
        formula <- raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version +
            acquisition_wave + condition
        term <- "conditionGBM.PD1"
    }

    # Compute coefficients, residuals and leverage directly from the same full-model design.
    design_matrix <- model.matrix(formula, data = cohort_table)
    response <- cohort_table$raw20_score
    n_samples <- nrow(design_matrix)
    n_parameters <- ncol(design_matrix)
    inverse_crossproduct <- solve(crossprod(design_matrix))
    coefficients <- as.vector(inverse_crossproduct %*% crossprod(design_matrix, response))
    residuals <- as.vector(response - design_matrix %*% coefficients)
    leverage <- rowSums((design_matrix %*% inverse_crossproduct) * design_matrix)
    residual_variance <- sum(residuals^2) / (n_samples - n_parameters)
    cook <- residuals^2 / (n_parameters * residual_variance) * leverage / (1 - leverage)^2

    # DFBETA uses the same offset from the full-model coefficient for each omitted sample.
    exposure_index <- match(term, colnames(design_matrix))
    exposure_dfbeta <- numeric(n_samples)
    for (i in seq_len(n_samples)) {
        coefficient_shift <- as.vector(inverse_crossproduct %*% design_matrix[i, ])
        coefficient_shift <- coefficient_shift * residuals[i] / (1 - leverage[i])
        exposure_dfbeta[i] <- coefficient_shift[exposure_index]
    }

    reference_table <- prior_influence_table[prior_influence_table$dataset == dataset, ]
    sample_order <- match(cohort_table$sample, reference_table$sample)
    cook_difference <- abs(cook - reference_table$cooks_distance[sample_order])
    dfbeta_difference <- abs(exposure_dfbeta - reference_table$exposure_dfbeta[sample_order])
    max_cook_difference <- max(cook_difference)
    max_dfbeta_difference <- max(dfbeta_difference)
    record_check(
        paste(dataset, "independent_cook"),
        length(sample_order) == n_samples && all(!is.na(sample_order)) && max_cook_difference < 1e-8,
        max_cook_difference
    )
    record_check(paste(dataset, "independent_exposure_dfbeta"), max_dfbeta_difference < 1e-8, max_dfbeta_difference)
    influence_tables[[dataset]] <- data.frame(
        dataset               = dataset,
        sample                = cohort_table$sample,
        cook_independent      = cook,
        cook_old              = reference_table$cooks_distance[sample_order],
        cook_abs_difference   = cook_difference,
        dfbeta_independent    = exposure_dfbeta,
        dfbeta_old            = reference_table$exposure_dfbeta[sample_order],
        dfbeta_abs_difference = dfbeta_difference,
        above_4_over_n        = cook > 4 / n_samples
    )

    # Center the non-intercept columns to calculate VIF from the inverse cross-product.
    predictor_columns <- setdiff(colnames(design_matrix), "(Intercept)")
    centered_matrix <- scale(design_matrix[, predictor_columns, drop = FALSE], center = TRUE, scale = FALSE)
    crossproduct <- crossprod(centered_matrix)
    vif <- diag(solve(crossproduct)) * diag(crossproduct)
    names(vif) <- colnames(centered_matrix)
    reference_vif_table <- prior_vif_table[prior_vif_table$dataset == dataset, ]
    vif_order <- match(names(vif), reference_vif_table$column)
    vif_difference <- abs(vif - reference_vif_table$vif[vif_order])
    max_vif_difference <- max(vif_difference)
    record_check(
        paste(dataset, "independent_vif"),
        all(!is.na(vif_order)) && max_vif_difference < 1e-8,
        max_vif_difference
    )
    vif_tables[[dataset]] <- data.frame(
        dataset = dataset,
        column  = names(vif),
        vif_independent = as.vector(vif),
        vif_old = reference_vif_table$vif[vif_order],
        absolute_difference = vif_difference
    )
}

# Write patient and coefficient diagnostics before the separate protein calculations.
influence_table_total <- do.call(rbind, influence_tables)
vif_table_total <- do.call(rbind, vif_tables)
checks_table_total <- do.call(rbind, check_tables)
write_new_table(influence_table_total, "independent_cook_dfbeta.csv")
write_new_table(vif_table_total, "independent_vif.csv")
write_new_table(checks_table_total, "independent_diagnostics_checks.csv")
cat("INDEPENDENT_DIAGNOSTICS ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")

# Use the saved paired differences and compare the independent PCG64 stream with the original sign flips.
check_tables <- list()
result_tables <- list()
protein_dir <- file.path(revision, "step5_protein_common")
score_table <- read.csv(file.path(protein_dir, "patient_scores.csv"))
effects_table <- read.csv(file.path(protein_dir, "effects.csv"))
record_check("four_existing_protein_effects", nrow(effects_table) == 4L)

for (i in seq_len(nrow(effects_table))) {
    effect_table <- effects_table[i, ]
    score_rows <- score_table$analysis == effect_table$analysis & score_table$mode == effect_table$mode
    delta <- score_table$delta[score_rows]
    n_pairs <- length(delta)
    record_check(
        paste(effect_table$analysis, effect_table$mode, "source_deltas"),
        n_pairs == effect_table$n_pairs && abs(mean(delta) - effect_table$mean_delta) < 1e-12
    )

    # Restart the approved PCG64 stream for each mode. Moving the seed changes the p values.
    dqRNGkind("pcg64")
    dqset.seed(20261020)
    simulations <- 100000L
    exceedances <- 0L
    observed <- abs(mean(delta))
    for (start in seq.int(1L, simulations, by = 2000L)) {
        batch_size <- min(2000L, simulations - start + 1L)
        signs <- matrix(dqrrademacher(batch_size * n_pairs), nrow = batch_size, ncol = n_pairs)
        null_mean <- as.vector(signs %*% delta) / n_pairs
        exceedances <- exceedances + sum(abs(null_mean) >= observed - 1e-12)
    }

    # The accepted difference is three Monte Carlo standard errors plus one simulation increment.
    p_previous <- effect_table$sign_flip_p
    p_independent <- (exceedances + 1) / (simulations + 1)
    standard_error <- sqrt(p_previous * (1 - p_previous) / simulations)
    tolerance <- 3 * standard_error + 1 / (simulations + 1)
    difference <- abs(p_independent - p_previous)
    result_tables[[i]] <- data.frame(
        analysis            = effect_table$analysis,
        mode                = effect_table$mode,
        n_pairs             = n_pairs,
        seed                = 20261020,
        engine              = "dqrng pcg64",
        simulations         = simulations,
        exceedances         = exceedances,
        p_independent       = p_independent,
        p_previous          = p_previous,
        absolute_difference = difference,
        mc_standard_error   = standard_error,
        prefixed_tolerance  = tolerance,
        within_mc_tolerance = difference <= tolerance,
        previous_method     = if (n_pairs <= 20L) "exact" else "base_R_sample_Monte_Carlo"
    )
    record_check(
        paste(effect_table$analysis, effect_table$mode, "mc_tolerance"),
        difference <= tolerance,
        paste("new", p_independent, "old", p_previous, "tolerance", tolerance)
    )
}

# Retain the engine and its version with the Monte Carlo results and checks.
result_table_total <- do.call(rbind, result_tables)
checks_table_total <- do.call(rbind, check_tables)
engine_table <- data.frame(
    package = "dqrng",
    version = as.character(packageVersion("dqrng")),
    engine  = "pcg64",
    seed    = 20261020,
    simulations = 100000
)
write_new_table(result_table_total, "independent_signflip.csv")
write_new_table(checks_table_total, "independent_signflip_checks.csv")
write_new_table(engine_table, "independent_signflip_engine.csv")
cat("INDEPENDENT_SIGNFLIP ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
