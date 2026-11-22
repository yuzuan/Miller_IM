#!/usr/bin/env Rscript

# Purpose  Compare unadjusted and myeloid-adjusted PD-1 associations and assess patient influence.
# Inputs   step2_pd1_camera/{sample_audit.csv,GSE154795_formal_samples.csv}.
#          NOA_LEGACY_ROOT/write/56_figure5_pd1_rebuild/Figure5/source_data/
#          Figure5C_patient_scores_and_adjusted_effects_source.csv (sample, cohort, raw20_score and covariates).
#          NOA_LEGACY_ROOT/write/41_mg_inflammatory_sci_rebuild/Figure5/:
#          GSE121810_raw20_score_models.csv and GSE154795_raw20_score_models_auto_strict_idhwt.csv.
# Outputs  step5_pd1_stability_final/{model_results.csv,covariate_distributions.csv,leave_one_out.csv,
#          influence.csv,vif.csv,formal_score_samples.csv,checks.csv}.
# Run      Rscript NOA_R1/scripts/05c_pd1_model_stability.R <scratch_revision_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_pd1_stability_final")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir)

# Record membership and model-comparability checks in the original order. No return value.
check_tables <- list()
record_check <- function(check_name, passed) {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(check = check_name, passed = isTRUE(passed))

    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
        stop(check_name)
    }
}

# Match each patient-level score to the formal PD-1 camera cohort before fitting.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
source_path <- file.path(
    project, "write", "56_figure5_pd1_rebuild", "Figure5", "source_data",
    "Figure5C_patient_scores_and_adjusted_effects_source.csv"
)
source_table <- read.csv(source_path, stringsAsFactors = FALSE)
camera_dir <- file.path(revision, "step2_pd1_camera")
audit_table <- read.csv(file.path(camera_dir, "sample_audit.csv"))
formal_table <- read.csv(file.path(camera_dir, "GSE154795_formal_samples.csv"))
samples_gse121810 <- source_table$sample[source_table$cohort == "GSE121810"]
samples_gse154795 <- source_table$sample[source_table$cohort == "GSE154795"]
record_check(
    "formal_121810",
    length(samples_gse121810) == 29L &&
        setequal(samples_gse121810, audit_table$sample[audit_table$dataset == "GSE121810"])
)
record_check("formal_154795", length(samples_gse154795) == 21L && setequal(samples_gse154795, formal_table$ID))

model_tables <- list()
covariate_tables <- list()
leave_one_out_tables <- list()
influence_tables <- list()
vif_tables <- list()
sample_tables <- list()

for (dataset in c("GSE121810", "GSE154795")) {

    # Match the original reference levels and preserve each formula's term order.
    cohort_table <- source_table[source_table$cohort == dataset, ]
    if (dataset == "GSE154795") {
        cohort_table$condition <- factor(cohort_table$condition, levels = c("GBM.rec", "GBM.PD1"))
        categorical_fields <- c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave")
        for (field in categorical_fields) cohort_table[, field] <- factor(cohort_table[, field])
        formulas <- list(
            unadjusted = raw20_score ~ condition,
            proxy_only = raw20_score ~ myeloid_abundance + condition,
            full       = raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version +
                acquisition_wave + condition
        )
        exposure_groups <- cohort_table$condition
    } else {
        cohort_table$group <- factor(cohort_table$group, levels = c("Adj", "Neo"))
        categorical_fields <- character()
        formulas <- list(
            unadjusted = raw20_score ~ group,
            proxy_only = raw20_score ~ myeloid_abundance + group,
            full       = raw20_score ~ myeloid_abundance + group
        )
        exposure_groups <- cohort_table$group
    }
    expected_samples <- if (dataset == "GSE154795") 21L else 29L
    record_check(
        paste(dataset, "scores_complete"),
        nrow(cohort_table) == expected_samples &&
            all(is.finite(cohort_table$raw20_score)) && all(is.finite(cohort_table$myeloid_abundance))
    )
    available_fields <- intersect(
        c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave", "n_myeloid_cells"),
        colnames(cohort_table)
    )
    sample_fields <- c("sample", "cohort", "group", "condition", "raw20_score", "myeloid_abundance", available_fields)
    sample_tables[[dataset]] <- cohort_table[, sample_fields, drop = FALSE]

    # Summarize each categorical level by exposure group, retaining explicit missing counts.
    for (field in categorical_fields) {
        for (group in levels(exposure_groups)) {
            values <- cohort_table[exposure_groups == group, field]
            for (level in levels(cohort_table[, field])) {
                covariate_tables[[length(covariate_tables) + 1L]] <- data.frame(
                    dataset  = dataset,
                    group    = group,
                    variable = field,
                    level    = level,
                    n        = sum(values == level, na.rm = TRUE),
                    missing  = sum(is.na(values)),
                    mean     = NA_real_,
                    sd       = NA_real_,
                    median   = NA_real_,
                    q1       = NA_real_,
                    q3       = NA_real_
                )
            }
        }
    }

    # Continuous summaries use the same group order and the original quantile defaults.
    for (field in c("myeloid_abundance", "n_myeloid_cells")) {
        if (!field %in% colnames(cohort_table)) next

        for (group in levels(exposure_groups)) {
            values <- cohort_table[exposure_groups == group, field]
            covariate_tables[[length(covariate_tables) + 1L]] <- data.frame(
                dataset  = dataset,
                group    = group,
                variable = field,
                level    = "continuous",
                n        = sum(!is.na(values)),
                missing  = sum(is.na(values)),
                mean     = mean(values, na.rm = TRUE),
                sd       = sd(values, na.rm = TRUE),
                median   = median(values, na.rm = TRUE),
                q1       = unname(quantile(values, .25, na.rm = TRUE)),
                q3       = unname(quantile(values, .75, na.rm = TRUE))
            )
        }
    }

    # Each row is an exposure coefficient from one formal patient-level model.
    for (label in names(formulas)) {
        fit <- lm(formulas[[label]], data = cohort_table)
        design_matrix <- model.matrix(fit)
        term <- grep("^(condition|group)", names(coef(fit)), value = TRUE)
        record_check(
            paste(dataset, label, "unique_exposure"),
            length(term) == 1L && qr(design_matrix)$rank == ncol(design_matrix)
        )
        coefficient <- coef(summary(fit))[term, ]
        interval <- confint(fit, term)
        model_table <- data.frame(
            dataset     = dataset,
            model       = label,
            formula     = paste(deparse(formulas[[label]]), collapse = ""),
            n           = nobs(fit),
            parameters  = length(coef(fit)),
            residual_df = df.residual(fit),
            contrast    = term,
            effect      = unname(coefficient["Estimate"]),
            ci_low      = unname(interval[1]),
            ci_high     = unname(interval[2]),
            t           = unname(coefficient["t value"]),
            p_value     = unname(coefficient["Pr(>|t|)"])
        )
        model_tables[[length(model_tables) + 1L]] <- model_table
        if (label != "full") next

        # Leave out one patient at a time using the full model. Keep a row for unestimable fits.
        for (i in seq_len(nrow(cohort_table))) {
            subset_table <- cohort_table[-i, , drop = FALSE]
            subset_fit <- lm(formulas[[label]], data = subset_table)
            subset_design <- model.matrix(subset_fit)
            subset_coefficients <- coef(subset_fit)
            estimate <- NA_real_
            if (term %in% names(subset_coefficients) && is.finite(subset_coefficients[term]) &&
                qr(subset_design)$rank == ncol(subset_design)) {
                estimate <- subset_coefficients[term]
            }
            leave_one_out_tables[[length(leave_one_out_tables) + 1L]] <- data.frame(
                dataset = dataset,
                sample  = cohort_table$sample[i],
                left_out_effect = estimate,
                reason  = if (is.na(estimate)) "rank_or_coefficient_unestimable" else ""
            )
        }

        # Cook distances and exposure DFBETA refer to the complete full-model fit.
        influence_tables[[length(influence_tables) + 1L]] <- data.frame(
            dataset         = dataset,
            sample          = cohort_table$sample,
            cooks_distance  = unname(cooks.distance(fit)),
            cooks_4_over_n  = 4 / nrow(cohort_table),
            exposure_dfbeta = unname(dfbeta(fit)[, term])
        )
        for (column in setdiff(colnames(design_matrix), "(Intercept)")) {
            design_column <- design_matrix[, column]
            other_columns <- setdiff(colnames(design_matrix), column)
            auxiliary_fit <- lm.fit(design_matrix[, other_columns, drop = FALSE], design_column)
            r_squared <- 1 - sum(auxiliary_fit$residuals^2) / sum((design_column - mean(design_column))^2)
            vif <- if (is.finite(r_squared) && r_squared < 1) 1 / (1 - r_squared) else NA_real_
            vif_tables[[length(vif_tables) + 1L]] <- data.frame(
                dataset   = dataset,
                column    = column,
                vif       = vif,
                r_squared = r_squared
            )
        }
    }
    if (dataset == "GSE154795") {
        record_check("full_7_parameters_14_df", model_table$parameters == 7L && model_table$residual_df == 14L)
    }

    # Verify the model shared with the July analysis before retaining additional sensitivity estimates.
    legacy_figure_dir <- file.path(project, "write", "41_mg_inflammatory_sci_rebuild", "Figure5")
    model_filename <- if (dataset == "GSE154795") {
        "GSE154795_raw20_score_models_auto_strict_idhwt.csv"
    } else {
        "GSE121810_raw20_score_models.csv"
    }
    legacy_table <- read.csv(file.path(legacy_figure_dir, model_filename))
    legacy_raw_table <- legacy_table[legacy_table$metric == "raw20_score", ]
    legacy_adjusted_table <- legacy_table[legacy_table$metric == "raw20_score_myeloid_adjusted", ]
    model_table_total <- do.call(rbind, model_tables)
    cohort_model_table <- model_table_total[model_table_total$dataset == dataset, ]
    unadjusted_table <- cohort_model_table[cohort_model_table$model == "unadjusted", ]
    full_table <- cohort_model_table[cohort_model_table$model == "full", ]
    record_check(
        paste(dataset, "original_adjusted_score_model_reproduced"),
        abs(full_table$effect - legacy_adjusted_table$effect) < 1e-12 &&
            abs(full_table$p_value - legacy_adjusted_table$p_value) < 1e-12
    )
    if (dataset == "GSE121810") {
        record_check(
            "GSE121810_original_unadjusted_reproduced",
            abs(unadjusted_table$effect - legacy_raw_table$effect) < 1e-12
        )
    }
}

# Preserve the original model order, patient order and output columns.
covariate_table_total <- do.call(rbind, covariate_tables)
leave_one_out_table_total <- do.call(rbind, leave_one_out_tables)
influence_table_total <- do.call(rbind, influence_tables)
vif_table_total <- do.call(rbind, vif_tables)
sample_table_total <- do.call(rbind, sample_tables)
checks_table_total <- do.call(rbind, check_tables)
write.csv(model_table_total, file.path(output_dir, "model_results.csv"), row.names = FALSE, na = "")
write.csv(covariate_table_total, file.path(output_dir, "covariate_distributions.csv"), row.names = FALSE, na = "")
write.csv(leave_one_out_table_total, file.path(output_dir, "leave_one_out.csv"), row.names = FALSE, na = "")
write.csv(influence_table_total, file.path(output_dir, "influence.csv"), row.names = FALSE, na = "")
write.csv(vif_table_total, file.path(output_dir, "vif.csv"), row.names = FALSE, na = "")
write.csv(sample_table_total, file.path(output_dir, "formal_score_samples.csv"), row.names = FALSE, na = "")
write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")

cat("PD1_STABILITY ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
