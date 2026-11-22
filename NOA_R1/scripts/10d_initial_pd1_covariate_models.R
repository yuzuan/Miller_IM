#!/usr/bin/env Rscript

# Purpose  Reproduce the July raw-score PD-1 models alongside the revision sensitivity models.
# Inputs   step5_analysis_plan.md and step5_pd1_stability_final/{formal_score_samples.csv,
#          model_results.csv,influence.csv} (cohort, exposure, raw20_score and model covariates).
#          NOA_LEGACY_ROOT/write/41_mg_inflammatory_sci_rebuild/Figure5/:
#          GSE121810_raw20_score_models.csv and GSE154795_raw20_score_models_auto_strict_idhwt.csv
#          (metric, effect, conf_low, conf_high and p_value).
# Outputs  step5_closeout/{original_raw_reproduced.csv,model_results_with_original_raw.csv,
#          cook_threshold_patients.csv,raw_checks.csv}.
# Run      Rscript NOA_R1/scripts/10d_initial_pd1_covariate_models.R <scratch_revision_dir>
# This historical entry keeps its original output directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_closeout")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir)

# Record checks in the original order and stop if reproduction differs from the July values. No return value.
check_tables <- list()
record_check <- function(check_name, passed) {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(check = check_name, passed = isTRUE(passed))

    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write.csv(checks_table_total, file.path(output_dir, "raw_checks.csv"), row.names = FALSE, na = "")
        stop(check_name)
    }
}

# Check the heading of the frozen closeout amendment in the English protocol.
plan <- readLines(file.path(revision, "step5_analysis_plan.md"), warn = FALSE, encoding = "UTF-8")
record_check("closeout_plan_frozen", any(grepl("Step 5 closeout addendum", plan, fixed = TRUE)))

stability_dir <- file.path(revision, "step5_pd1_stability_final")
formal_table <- read.csv(file.path(stability_dir, "formal_score_samples.csv"))
prior_table <- read.csv(file.path(stability_dir, "model_results.csv"))

# GSE121810 compares Neo with Adj. GSE154795 uses the original clinical and technical adjustment.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
legacy_figure_dir <- file.path(project, "write", "41_mg_inflammatory_sci_rebuild", "Figure5")
raw_model_tables <- list()
for (dataset in c("GSE121810", "GSE154795")) {
    sample_table <- formal_table[formal_table$cohort == dataset, ]
    if (dataset == "GSE121810") {
        sample_table$group <- factor(sample_table$group, levels = c("Adj", "Neo"))
        formula <- raw20_score ~ group
        term <- "groupNeo"
        covariates <- "group"
        legacy_model_table <- read.csv(file.path(legacy_figure_dir, "GSE121810_raw20_score_models.csv"))
    } else {
        sample_table$condition <- factor(sample_table$condition, levels = c("GBM.rec", "GBM.PD1"))
        for (field in c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave")) {
            sample_table[, field] <- factor(sample_table[, field])
        }
        formula <- raw20_score ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition
        term <- "conditionGBM.PD1"
        covariates <- "MGMT;sequencing_method;gel_bead_version;acquisition_wave;condition"
        legacy_model_path <- file.path(legacy_figure_dir, "GSE154795_raw20_score_models_auto_strict_idhwt.csv")
        legacy_model_table <- read.csv(legacy_model_path)
    }

    # Keep the full precision coefficients, interval and p value for the archived comparison.
    fit <- lm(formula, data = sample_table)
    coefficient <- coef(summary(fit))[term, ]
    interval <- confint(fit, term)
    effect <- unname(coefficient["Estimate"])
    p_value <- unname(coefficient["Pr(>|t|)"])
    legacy_raw_table <- legacy_model_table[legacy_model_table$metric == "raw20_score", ]
    record_check(paste(dataset, "original_raw_one_row"), nrow(legacy_raw_table) == 1L)
    record_check(paste(dataset, "original_raw_effect_reproduced"), abs(effect - legacy_raw_table$effect) < 1e-12)

    # Compare each endpoint as a scalar so base R can take the absolute differences.
    differences <- c(
        interval[1] - legacy_raw_table$conf_low, interval[2] - legacy_raw_table$conf_high,
        p_value - legacy_raw_table$p_value
    )
    record_check(paste(dataset, "original_raw_ci_p_reproduced"), max(abs(differences)) < 1e-12)

    # The review response records these rounded GSE154795 values.
    if (dataset == "GSE154795") {
        record_check("reviewer_raw_0.41943_p0.17015", abs(effect - .41943) < .00001 && abs(p_value - .17015) < .00001)
    }
    raw_model_tables[[dataset]] <- data.frame(
        dataset     = dataset,
        model       = "original_raw",
        covariates  = covariates,
        n           = nobs(fit),
        parameters  = length(coef(fit)),
        residual_df = df.residual(fit),
        contrast    = term,
        effect      = effect,
        ci_low      = unname(interval[1]),
        ci_high     = unname(interval[2]),
        t           = unname(coefficient["t value"]),
        p_value     = p_value
    )
}

# The six archived rows are ordered by cohort and then unadjusted, proxy-only and full models.
prior_table$covariates <- c(
    "group", "myeloid_abundance;group", "myeloid_abundance;group",
    "condition", "myeloid_abundance;condition",
    "myeloid_abundance;MGMT;sequencing_method;gel_bead_version;acquisition_wave;condition"
)
raw_model_table_total <- do.call(rbind, raw_model_tables)
prior_table <- prior_table[, colnames(raw_model_table_total)]
model_table_total <- rbind(prior_table, raw_model_table_total)
model_order <- match(model_table_total$model, c("unadjusted", "original_raw", "proxy_only", "full"))
model_table_total <- model_table_total[order(model_table_total$dataset, model_order), ]
record_check("eight_parallel_rows", nrow(model_table_total) == 8L && all(table(model_table_total$dataset) == 4L))
write.csv(raw_model_table_total, file.path(output_dir, "original_raw_reproduced.csv"), row.names = FALSE, na = "")
write.csv(model_table_total, file.path(output_dir, "model_results_with_original_raw.csv"), row.names = FALSE, na = "")

# Flag Cook distances against 4 / n within the original 29 and 21 analytical samples.
influence_table <- read.csv(file.path(stability_dir, "influence.csv"))
cook_tables <- list()
for (dataset in c("GSE121810", "GSE154795")) {
    cohort_table <- influence_table[influence_table$dataset == dataset, ]
    record_check(paste(dataset, "cook_threshold"), nrow(cohort_table) == if (dataset == "GSE121810") 29L else 21L)
    cook_tables[[dataset]] <- data.frame(
        dataset         = dataset,
        sample          = cohort_table$sample,
        cooks_distance  = cohort_table$cooks_distance,
        threshold       = 4 / nrow(cohort_table),
        above_4_over_n  = cohort_table$cooks_distance > 4 / nrow(cohort_table),
        exposure_dfbeta = cohort_table$exposure_dfbeta
    )
}

cook_table_total <- do.call(rbind, cook_tables)
checks_table_total <- do.call(rbind, check_tables)
write.csv(cook_table_total, file.path(output_dir, "cook_threshold_patients.csv"), row.names = FALSE, na = "")
write.csv(checks_table_total, file.path(output_dir, "raw_checks.csv"), row.names = FALSE, na = "")

cat("RAW_CLOSEOUT ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
