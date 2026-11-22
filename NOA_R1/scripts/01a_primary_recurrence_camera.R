#!/usr/bin/env Rscript
# Purpose: Rebuild the frozen paired recurrence CAMERA and patient score analyses.
# Inputs: provenance/step1_input_manifest.csv (cohort, matrix paths and expected totals), counts and sample metadata;
#   legacy Figure 1 fold audit (omitted_patient), raw20 definition (raw20_order, gene), and IDH crosswalk (pair_key);
#   legacy Step 38 scores by unit, paired score summary and edgeR tested genes at the original threshold.
# Outputs under step1_primary_recurrence/: camera_results.csv, patient_score_effects.csv, sample_audit.csv,
#   patient_paired_changes.csv, input_audit.csv, gene_coverage.csv, validation_checks.csv, parameters.json,
#   sessionInfo.txt, and <dataset>_tested_genes.csv, _design_matrix.csv, _paired_voom.rds and _voom_trend.pdf.
# Run: Rscript NOA_R1/scripts/01a_primary_recurrence_camera.R NOA_R1
# NOA_LEGACY_ROOT supplies the prior analyses. NOA_PUBLIC_DATA_ROOT supplies the manifest's public-data tree.

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(jsonlite)
})

# The revision directory contains the manifest and receives the results.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 01a_primary_recurrence_camera.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)

manifest_table <- read.csv(file.path(revision, "provenance", "step1_input_manifest.csv"), check.names = FALSE)

# The fold audit supplies formal patient membership, and the definition supplies the ranked fixed program.
figure_source <- file.path(project, "write", "51_figure1_final_reorganized", "Figure1", "source_data")
formal_table <- read.csv(file.path(figure_source, "Figure1F_leave_one_patient_out_fold_audit.csv"))
definition_table <- read.csv(file.path(figure_source, "Figure1_raw20_definition_source.csv"))
signature_table <- unique(definition_table[, c("raw20_order", "gene")])
signature_table <- signature_table[order(signature_table$raw20_order), ]
signature <- as.character(signature_table$gene)

# Step 38 reference tables use the original 20-cell threshold.
source_dir <- file.path(project, "write", "38_independent_cohort_mg_inflammatory_recalculation")
reference_unit_table <- read.csv(file.path(source_dir, "independent_fixed_program_scores_by_unit.csv"))
reference_unit_table <- reference_unit_table[reference_unit_table$threshold == 20L, ]
reference_summary_table <- read.csv(file.path(source_dir, "independent_fixed_program_paired_scores.csv"))
reference_summary_table <- reference_summary_table[reference_summary_table$threshold == 20L, ]
reference_gene_table <- read.csv(file.path(source_dir, "independent_paired_edger_all_genes.csv"))
reference_gene_table <- reference_gene_table[reference_gene_table$threshold == 20L, ]

# Prefer the versioned public crosswalk. The original prepared tree keeps the same file under clinical_metadata.
crosswalk_path <- file.path(project, "config", "GSE174554_formal_pair_IDH_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
    crosswalk_path <- file.path(project, "data", "clinical_metadata", "GSE174554_formal_pair_IDH_crosswalk.csv")
}
idh_table <- read.csv(crosswalk_path)

output_dir <- file.path(revision, "step1_primary_recurrence")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Record each matching check and save the failure to the same validation table before stopping.
checks <- list()

check <- function(dataset, label, condition, detail = "") {
    passed <- isTRUE(condition)
    checks[[length(checks) + 1L]] <<- data.frame(dataset = dataset, check = label, passed = passed, detail = detail)

    if (!passed) {
        validation_table_total <- do.call(rbind, checks)
        write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
        stop(dataset, ": ", label, "; ", detail)
    }
}

check(
    "Both",
    "fixed_twenty_gene_order",
    identical(signature_table$raw20_order, 1:20) && length(unique(signature)) == 20L
)
check("Both", "manifest_two_cohorts", identical(manifest_table$dataset, c("GSE174554", "GSE274546")))

# Use configured roots to open the two manifest input columns. Keep their original tokens in the audit table.
public_root <- Sys.getenv("NOA_PUBLIC_DATA_ROOT")
if (!nzchar(public_root)) public_root <- Sys.getenv("MILLER_IM_DATA_ROOT")
if (!nzchar(public_root)) public_root <- file.path(project, "data")
public_root <- normalizePath(path.expand(public_root), winslash = "/", mustWork = FALSE)

camera_tables <- list()
score_tables <- list()
sample_tables <- list()
pair_tables <- list()
input_tables <- list()
coverage_tables <- list()
model_files <- character()

for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, ]
    dataset <- entry_table$dataset
    counts_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$counts_path)
    metadata_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$metadata_path)
    input_paths <- c(counts_path, metadata_path)
    missing_paths <- input_paths[!file.exists(input_paths)]
    check(dataset, "inputs_exist", length(missing_paths) == 0L, paste(missing_paths, collapse = ", "))

    counts_table <- read.csv(gzfile(counts_path), row.names = 1L, check.names = FALSE)
    counts_matrix <- as.matrix(counts_table)
    sample_table <- read.csv(metadata_path, stringsAsFactors = FALSE)
    rm(counts_table)

    # Rows are genes and columns are pseudobulk samples. Sample metadata must already match column order.
    check(dataset, "raw_counts_valid", all(is.finite(counts_matrix)) && all(counts_matrix >= 0) &&
        all(counts_matrix == round(counts_matrix)))
    check(
        dataset,
        "unique_gene_and_sample_ids",
        !anyDuplicated(rownames(counts_matrix)) && !anyDuplicated(colnames(counts_matrix)) &&
            !anyDuplicated(sample_table$sample_id)
    )
    check(dataset, "counts_metadata_order", identical(colnames(counts_matrix), sample_table$sample_id))
    check(dataset, "input_gene_count", nrow(counts_matrix) == entry_table$expected_input_genes)

    # Keep the original formal patients and their two time points at the frozen cell threshold.
    in_fold_audit <- formal_table$dataset == dataset & formal_table$run_type == "Leave-one-patient-out"
    patients <- sort(unique(formal_table$omitted_patient[in_fold_audit]))
    check(dataset, "formal_patient_count", length(patients) == entry_table$expected_pairs)
    sample_table <- sample_table[sample_table$patient_id %in% patients, , drop = FALSE]
    counts_matrix <- counts_matrix[, sample_table$sample_id, drop = FALSE]
    check(
        dataset,
        "one_sample_per_patient_time",
        nrow(sample_table) == 2L * entry_table$expected_pairs &&
            all(table(sample_table$patient_id, sample_table$condition) == 1L)
    )
    check(
        dataset,
        "time_labels_and_cell_threshold",
        setequal(sample_table$condition, c("Primary", "Recurrent")) && all(sample_table$n_cells >= 20L)
    )
    primary_cells <- sum(sample_table$n_cells[sample_table$condition == "Primary"])
    recurrent_cells <- sum(sample_table$n_cells[sample_table$condition == "Recurrent"])
    check(dataset, "formal_cell_counts", primary_cells == entry_table$expected_primary_cells &&
        recurrent_cells == entry_table$expected_recurrent_cells)

    cohort_unit_table <- reference_unit_table[reference_unit_table$dataset == dataset, ]
    check(dataset, "same_formal_units_as_step38", setequal(sample_table$sample_id, cohort_unit_table$unit_id) &&
        setequal(patients, cohort_unit_table$pair_key))
    if (dataset == "GSE174554") {
        wildtype_rows <- toupper(as.character(idh_table$include_idh_wildtype)) == "TRUE"
        check(
            dataset,
            "idh_wildtype_patients",
            setequal(patients, idh_table$pair_key[wildtype_rows]) && !"GSE174554_pair33" %in% patients
        )
    }

    # AC253572.2 is the fixed program member absent from both source matrices.
    expected_present <- signature[signature != "AC253572.2"]
    check(
        dataset,
        "nineteen_genes_in_raw_counts",
        identical(intersect(signature, rownames(counts_matrix)), expected_present)
    )

    # Pair fixed effects and Primary as the reference define the recurrent-minus-primary contrast.
    sample_table$pair_key <- factor(sample_table$patient_id)
    sample_table$condition <- factor(sample_table$condition, levels = c("Primary", "Recurrent"))
    design <- model.matrix(~ pair_key + condition, data = sample_table)
    coefficient <- which(colnames(design) == "conditionRecurrent")
    design_rank <- qr(design)$rank
    residual_df <- nrow(design) - design_rank
    check(dataset, "paired_full_rank_design", design_rank == ncol(design) && length(coefficient) == 1L &&
        residual_df == entry_table$expected_pairs - 1L)

    # CAMERA uses the original filterByExpr universe, filtered TMM and voom observation weights.
    filtered_counts <- DGEList(counts = counts_matrix)
    keep <- filterByExpr(filtered_counts, design = design)
    filtered_counts <- filtered_counts[keep, , keep.lib.sizes = FALSE]
    filtered_counts <- calcNormFactors(filtered_counts, method = "TMM")
    reference_genes <- reference_gene_table$gene[reference_gene_table$dataset == dataset]
    check(dataset, "original_tested_gene_universe", nrow(filtered_counts) == entry_table$expected_tested_genes &&
        setequal(rownames(filtered_counts), reference_genes))
    present <- intersect(signature, rownames(filtered_counts))
    check(dataset, "same_nineteen_genes_after_filter", identical(present, expected_present))

    # voom draws its base diagnostic while choosing the adaptive span. Save the original 6 x 4 inch PDF.
    voom_path <- file.path(output_dir, paste0(dataset, "_voom_trend.pdf"))
    pdf(voom_path, width = 6, height = 4)
        voom_fit <- voom(
            filtered_counts,
            design           = design,
            normalize.method = "none",
            adaptive.span    = TRUE,
            plot             = TRUE,
            save.plot        = TRUE
        )
    dev.off()
    if (!file.exists(voom_path) || file.size(voom_path) == 0) stop("Missing voom diagnostic for ", dataset)
    check(dataset, "finite_voom_expression_and_weights", all(is.finite(voom_fit$E)) &&
        all(is.finite(voom_fit$weights)) && all(voom_fit$weights > 0))
    check(dataset, "actual_voom_span_recorded", length(voom_fit$span) == 1L && is.finite(voom_fit$span) &&
        voom_fit$span > 0 && voom_fit$span <= 1)

    # Run the estimated-correlation primary test and the 0.01 sensitivity test on the same contrast.
    indices <- list(inflammatory_program_top20 = match(present, rownames(voom_fit)))
    for (mode in c("estimated", "fixed_0.01")) {
        rho <- if (mode == "estimated") NA_real_ else 0.01
        camera_table <- camera(
            voom_fit,
            index          = indices,
            design         = design,
            contrast       = coefficient,
            inter.gene.cor = rho,
            allow.neg.cor  = FALSE,
            use.ranks      = FALSE,
            directional    = TRUE,
            trend.var      = FALSE,
            sort           = FALSE
        )
        observed_rho <- if (mode == "estimated") camera_table$Correlation else rho
        check(
            dataset,
            paste0("camera_valid_", mode),
            nrow(camera_table) == 1L && camera_table$NGenes == 19L && is.finite(camera_table$PValue) &&
                camera_table$PValue >= 0 && camera_table$PValue <= 1 && is.finite(observed_rho)
        )

        # Estimated correlation uses the residual DF cap. Fixed correlation uses the tested-gene reference DF.
        reference_df <- if (mode == "estimated") min(residual_df, nrow(voom_fit) - 2L) else nrow(voom_fit) - 2L
        camera_tables[[length(camera_tables) + 1L]] <- data.frame(
            dataset               = dataset,
            role                  = if (mode == "estimated") "primary" else "sensitivity",
            correlation_mode      = mode,
            n_pairs               = entry_table$expected_pairs,
            n_genes               = camera_table$NGenes,
            tested_universe       = nrow(voom_fit),
            residual_df           = residual_df,
            camera_reference_df   = reference_df,
            correlation           = observed_rho,
            effective_correlation = max(0, observed_rho),
            effective_vif         = max(1, 1 + (length(present) - 1L) * observed_rho),
            direction             = camera_table$Direction,
            p_value               = camera_table$PValue,
            test                  = "two-sided competitive camera",
            contrast              = "Recurrent minus Primary"
        )
    }

    # Patient scores use TMM on all genes and prior count 1, then average gene-wise z scores within each cohort.
    score_counts <- calcNormFactors(DGEList(counts = counts_matrix), method = "TMM")
    log_cpm_matrix <- cpm(score_counts, log = TRUE, prior.count = 1)
    score_gene_matrix <- log_cpm_matrix[present, , drop = FALSE]
    gene_sd <- apply(score_gene_matrix, 1L, sd)
    check(dataset, "score_gene_sd_positive", all(is.finite(gene_sd)) && all(gene_sd > 0))
    standardized_matrix <- t(scale(t(score_gene_matrix)))
    scores <- colMeans(standardized_matrix)
    check(dataset, "scores_finite", all(is.finite(standardized_matrix)) && all(is.finite(scores)))
    reference_scores <- cohort_unit_table$score[match(sample_table$sample_id, cohort_unit_table$unit_id)]
    score_error <- max(abs(scores - reference_scores))
    check(
        dataset,
        "score_reproduces_step38",
        is.finite(score_error) && score_error < 1e-10,
        paste("maximum absolute error", signif(score_error, 4))
    )

    # Pair scores by patient. Keep the original order and the within-patient subtraction.
    pair_table <- data.frame(dataset = dataset, patient_id = patients)
    primary <- sample_table$condition == "Primary"
    recurrent <- sample_table$condition == "Recurrent"
    primary_rows <- match(patients, as.character(sample_table$patient_id[primary]))
    recurrent_rows <- match(patients, as.character(sample_table$patient_id[recurrent]))
    pair_table$primary_score <- scores[primary][primary_rows]
    pair_table$recurrent_score <- scores[recurrent][recurrent_rows]
    pair_table$delta <- pair_table$recurrent_score - pair_table$primary_score
    pair_table$primary_n_cells <- sample_table$n_cells[primary][primary_rows]
    pair_table$recurrent_n_cells <- sample_table$n_cells[recurrent][recurrent_rows]
    reference_deltas <- cohort_unit_table$delta_recurrent_minus_primary[match(patients, cohort_unit_table$pair_key)]
    check(dataset, "patient_deltas_reproduce_step38", all(is.finite(pair_table$delta)) &&
        max(abs(pair_table$delta - reference_deltas)) < 1e-10)

    # The paired t interval is a one-sample interval on the deltas. Retain the original sign-flip result fields.
    paired_test <- t.test(pair_table$delta, mu = 0, alternative = "two.sided", conf.level = 0.95)
    cohort_summary_table <- reference_summary_table[reference_summary_table$dataset == dataset, ]
    check(dataset, "effect_reproduces_step38", nrow(cohort_summary_table) == 1L &&
        abs(mean(pair_table$delta) - cohort_summary_table$mean_delta) < 1e-10 &&
        sum(pair_table$delta > 0) == cohort_summary_table$positive_pairs)
    score_tables[[dataset]] <- data.frame(
        dataset        = dataset,
        n_pairs        = nrow(pair_table),
        n_genes        = length(present),
        mean_delta     = mean(pair_table$delta),
        sd_delta       = sd(pair_table$delta),
        se_delta       = sd(pair_table$delta) / sqrt(nrow(pair_table)),
        ci_low         = paired_test$conf.int[[1]],
        ci_high        = paired_test$conf.int[[2]],
        ci_method      = "two-sided paired t, df=n_pairs-1",
        t_statistic    = unname(paired_test$statistic),
        residual_df    = unname(paired_test$parameter),
        t_p_value      = paired_test$p.value,
        original_sign_flip_p = cohort_summary_table$p_value,
        original_sign_flip_exact = cohort_summary_table$sign_flip_exact,
        positive_pairs = sum(pair_table$delta > 0),
        median_delta   = median(pair_table$delta),
        score_unit     = "mean gene-wise z-score of TMM logCPM; cohort-specific scaling"
    )

    # Save one audit row per sample, pair, gene in the fixed program, and cohort, respectively.
    sample_tables[[dataset]] <- data.frame(
        dataset                  = dataset,
        sample_id                = sample_table$sample_id,
        patient_id               = as.character(sample_table$patient_id),
        condition                = as.character(sample_table$condition),
        n_cells                  = sample_table$n_cells,
        raw_library_size         = colSums(counts_matrix),
        filtered_library_size    = filtered_counts$samples$lib.size,
        tmm_factor_camera        = filtered_counts$samples$norm.factors,
        effective_library_camera = filtered_counts$samples$lib.size * filtered_counts$samples$norm.factors,
        tmm_factor_score         = score_counts$samples$norm.factors,
        score                    = scores
    )
    pair_tables[[dataset]] <- pair_table
    coverage_tables[[dataset]] <- data.frame(
        dataset = dataset,
        order   = seq_along(signature),
        gene    = signature,
        present_in_raw = signature %in% rownames(counts_matrix),
        included_in_camera = signature %in% rownames(filtered_counts)
    )
    input_tables[[dataset]] <- data.frame(
        dataset            = dataset,
        n_pairs            = nrow(pair_table),
        n_samples          = ncol(counts_matrix),
        primary_cells      = sum(pair_table$primary_n_cells),
        recurrent_cells    = sum(pair_table$recurrent_n_cells),
        input_genes        = nrow(counts_matrix),
        tested_genes       = nrow(filtered_counts),
        design_rank        = design_rank,
        residual_df        = residual_df,
        contrast           = "conditionRecurrent",
        selected_voom_span = voom_fit$span,
        maximum_score_reproduction_error = score_error,
        counts_path        = entry_table$counts_path,
        metadata_path      = entry_table$metadata_path
    )

    # The comparator scripts read the paired model RDS with its original named components.
    tested_gene_table <- data.frame(gene = rownames(filtered_counts))
    design_table <- data.frame(sample_id = sample_table$sample_id, design, check.names = FALSE)
    tested_gene_path <- file.path(output_dir, paste0(dataset, "_tested_genes.csv"))
    write.csv(tested_gene_table, tested_gene_path, row.names = FALSE, na = "")
    write.csv(design_table, file.path(output_dir, paste0(dataset, "_design_matrix.csv")), row.names = FALSE, na = "")
    model <- list(
        y           = filtered_counts,
        voom        = voom_fit,
        metadata    = sample_table,
        design      = design,
        coefficient = coefficient,
        signature   = signature,
        present     = present
    )
    model_file <- paste0("step1_primary_recurrence/", dataset, "_paired_voom.rds")
    saveRDS(model, file.path(revision, model_file))
    model_files <- c(model_files, model_file)
}

# Combine cohort results in manifest order and retain the existing schemas and paths.
camera_table_total <- do.call(rbind, camera_tables)
score_table_total <- do.call(rbind, score_tables)
sample_table_total <- do.call(rbind, sample_tables)
pair_table_total <- do.call(rbind, pair_tables)
input_table_total <- do.call(rbind, input_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
validation_table_total <- do.call(rbind, checks)
write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "patient_score_effects.csv"), row.names = FALSE, na = "")
write.csv(sample_table_total, file.path(output_dir, "sample_audit.csv"), row.names = FALSE, na = "")
write.csv(pair_table_total, file.path(output_dir, "patient_paired_changes.csv"), row.names = FALSE, na = "")
write.csv(input_table_total, file.path(output_dir, "input_audit.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "gene_coverage.csv"), row.names = FALSE, na = "")
write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")

# Keep the declared methods and record the versions used to reconstruct the cached models.
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
package_versions <- list()
for (package in c("edgeR", "limma", "statmod")) package_versions[[package]] <- as.character(packageVersion(package))
parameters <- list(
    analysis_date  = "2026-10-19",
    design         = "~ pair_key + condition",
    reference      = "Primary",
    gene_set       = "original fixed top20; 19 measurable, AC253572.2 absent",
    primary_camera = list(
        inter_gene_cor = "estimated (NA)",
        allow_neg_cor  = FALSE,
        use_ranks      = FALSE,
        directional    = TRUE,
        two_sided      = TRUE
    ),
    sensitivity_camera = list(
        inter_gene_cor = 0.01,
        allow_neg_cor  = FALSE,
        use_ranks      = FALSE,
        directional    = TRUE,
        two_sided      = TRUE
    ),
    voom             = list(normalize_method = "none", adaptive_span = TRUE),
    patient_score    = "TMM on unfiltered formal matrix; logCPM prior.count=1; per-gene cohort z-score; mean19",
    patient_interval = "paired t 95% CI; conditional on the cohort-defined standardized scores",
    packages         = package_versions,
    cached_models    = model_files
)
write_json(parameters, file.path(output_dir, "parameters.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Validation checks passed: ", nrow(validation_table_total), "\n", sep = "")
