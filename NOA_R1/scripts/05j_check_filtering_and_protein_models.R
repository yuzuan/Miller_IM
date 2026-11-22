#!/usr/bin/env Rscript

# Purpose  Independently refit the PD-1, filtering, protein and rank analyses from their source inputs.
# Inputs   NOA_PUBLIC_DATA_ROOT/write/31_proteomic_validation/source_metadata/PDC000514/
#          KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv (Gene and paired Log/Unshared Log columns).
#          step5_protein_common/{protein_fixed_members.csv,patient_scores.csv,effects.csv}.
#          step5_rank_unification_v2/{sqrtF_fgsea.csv,ranked_genes.csv} and target_top20 in
#          step2_locked_gene_sets.csv or step2_extended_specificity/locked_gene_sets_used.csv.
#          NOA_LEGACY_ROOT/write/41_mg_inflammatory_sci_rebuild/Figure5/:
#          GSE121810_raw20_edger.csv and GSE154795_raw20_edger_auto_strict_idhwt.csv (gene, logFC and F).
#          Legacy Figure5C_patient_scores_and_adjusted_effects_source.csv (write/56_figure5_pd1_rebuild).
#          step5_pd1_stability_final/{model_results.csv,leave_one_out.csv}, step5_filter_consistency/camera_results.csv,
#          provenance/step1_input_manifest.csv (or step1_input_manifest.csv), its pseudobulk counts, and
#          step1_primary_recurrence/{GSE174554,GSE274546}_paired_voom.rds (metadata and signature).
# Outputs  independent_checks.csv in step5_pd1_stability_final, step5_filter_consistency,
#          step5_protein_common and step5_rank_unification_v2.
# Run      Rscript NOA_R1/scripts/05j_check_filtering_and_protein_models.R <scratch_revision_dir>

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(data.table)
    library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)

# A failed independent check stops before writing the report for its analysis. No return value.
check_tables <- list()
record_check <- function(check_name, passed, detail = "") {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(
        check  = check_name,
        passed = isTRUE(passed),
        detail = as.character(detail)
    )
    if (!isTRUE(passed)) stop(check_name)
}

# Refit the archived models from July patient scores using the original factor levels.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
source_path <- file.path(
    project, "write", "56_figure5_pd1_rebuild", "Figure5", "source_data",
    "Figure5C_patient_scores_and_adjusted_effects_source.csv"
)
source_table <- read.csv(source_path)
stability_dir <- file.path(revision, "step5_pd1_stability_final")
model_reference_table <- read.csv(file.path(stability_dir, "model_results.csv"))
loo_reference_table <- read.csv(file.path(stability_dir, "leave_one_out.csv"))
for (dataset in c("GSE121810", "GSE154795")) {
    cohort_table <- source_table[source_table$cohort == dataset, ]
    if (dataset == "GSE121810") {
        cohort_table$group <- factor(cohort_table$group, levels = c("Adj", "Neo"))
        formulas <- list(
            unadjusted = raw20_score ~ group,
            proxy_only = raw20_score ~ myeloid_abundance + group,
            full       = raw20_score ~ myeloid_abundance + group
        )
    } else {
        cohort_table$condition <- factor(cohort_table$condition, levels = c("GBM.rec", "GBM.PD1"))
        for (field in c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave")) {
            cohort_table[, field] <- factor(cohort_table[, field])
        }
        formulas <- list(
            unadjusted = raw20_score ~ condition,
            proxy_only = raw20_score ~ myeloid_abundance + condition,
            full       = raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version +
                acquisition_wave + condition
        )
    }

    # The exposure term is last in the frozen formula. Compare its estimate, interval and p value.
    for (label in names(formulas)) {
        fit <- lm(formulas[[label]], data = cohort_table)
        term <- tail(names(coef(fit)), 1L)
        reference_table <- model_reference_table[
            model_reference_table$dataset == dataset & model_reference_table$model == label,
        ]
        interval <- confint(fit, term)
        effect <- unname(coef(fit)[term])
        p_value <- coef(summary(fit))[term, "Pr(>|t|)"]
        record_check(
            paste(dataset, label, "independent_lm"),
            nrow(reference_table) == 1L && nobs(fit) == reference_table$n &&
                df.residual(fit) == reference_table$residual_df && abs(effect - reference_table$effect) < 1e-10 &&
                abs(interval[1] - reference_table$ci_low) < 1e-10 &&
                abs(interval[2] - reference_table$ci_high) < 1e-10 &&
                abs(p_value - reference_table$p_value) < 1e-10
        )
    }

    # Leave each patient out in the original sample order and compare the exposure coefficients.
    cohort_loo_table <- loo_reference_table[loo_reference_table$dataset == dataset, ]
    estimates <- numeric(nrow(cohort_table))
    for (i in seq_len(nrow(cohort_table))) {
        subset_table <- cohort_table[-i, , drop = FALSE]
        fit <- lm(formulas$full, data = subset_table)
        estimates[i] <- tail(coef(fit), 1L)
    }
    sample_order <- match(cohort_table$sample, cohort_loo_table$sample)
    difference <- estimates - cohort_loo_table$left_out_effect[sample_order]
    record_check(
        paste(dataset, "independent_leave_one_out"),
        all(is.finite(estimates)) && max(abs(difference)) < 1e-10
    )
}
checks_table_total <- do.call(rbind, check_tables)
write.csv(checks_table_total, file.path(stability_dir, "independent_checks.csv"), row.names = FALSE)
cat("MODEL_INDEPENDENT ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")

# Independently refit the condition-filtered models from paired Step 1 pseudobulk counts.
check_tables <- list()
manifest_path <- file.path(revision, "provenance", "step1_input_manifest.csv")
if (!file.exists(manifest_path)) manifest_path <- file.path(revision, "step1_input_manifest.csv")
manifest_table <- read.csv(manifest_path)
filter_dir <- file.path(revision, "step5_filter_consistency")
camera_reference_table <- read.csv(file.path(filter_dir, "camera_results.csv"))

for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, ]
    dataset <- entry_table$dataset
    paired_input <- readRDS(file.path(revision, "step1_primary_recurrence", paste0(dataset, "_paired_voom.rds")))
    metadata_table <- paired_input$metadata
    metadata_table$condition <- factor(metadata_table$condition, levels = c("Primary", "Recurrent"))
    metadata_table$pair_key <- factor(metadata_table$patient_id)

    # Substitute the public-data root for the placeholder used in the published manifest.
    counts_path <- entry_table$counts_path
    if (startsWith(counts_path, "<public_data_root>/")) {
        public_data_root <- Sys.getenv(
            "NOA_PUBLIC_DATA_ROOT",
            unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
        )
        counts_path <- file.path(public_data_root, sub("^<public_data_root>/", "", counts_path))
    }
    counts_table <- read.csv(gzfile(counts_path), row.names = 1, check.names = FALSE)
    counts_matrix <- as.matrix(counts_table[, as.character(metadata_table$sample_id), drop = FALSE])
    design_matrix <- model.matrix(~ pair_key + condition, data = metadata_table)
    all_counts <- DGEList(counts = counts_matrix)
    keep <- filterByExpr(all_counts, group = metadata_table$condition)
    condition_counts <- calcNormFactors(all_counts[keep, , keep.lib.sizes = FALSE], method = "TMM")
    voom_data <- voom(condition_counts, design_matrix, normalize.method = "none", adaptive.span = TRUE, plot = FALSE)

    # Re-evaluate the same target and recurrence coefficient under both frozen correlation settings.
    target_index <- which(rownames(voom_data) %in% paired_input$signature)
    contrast <- which(colnames(design_matrix) == "conditionRecurrent")
    for (mode in c("estimated", "fixed_0.01")) {
        correlation <- if (mode == "estimated") NA_real_ else .01
        camera_table <- camera(
            voom_data, list(program = target_index),
            design         = design_matrix,
            contrast       = contrast,
            inter.gene.cor = correlation,
            allow.neg.cor  = FALSE,
            use.ranks      = FALSE,
            directional    = TRUE,
            trend.var      = FALSE,
            sort           = FALSE
        )
        reference_table <- camera_reference_table[
            camera_reference_table$dataset == dataset & camera_reference_table$mode == mode,
        ]
        record_check(
            paste(dataset, mode, "independent_camera"),
            nrow(reference_table) == 1L && nrow(voom_data) == reference_table$background &&
                camera_table$NGenes[1] == reference_table$n_genes &&
                camera_table$Direction[1] == reference_table$direction &&
                abs(camera_table$PValue[1] - reference_table$p_value) < 1e-10
        )
    }
}
checks_table_total <- do.call(rbind, check_tables)
write.csv(checks_table_total, file.path(filter_dir, "independent_checks.csv"), row.names = FALSE)
cat("FILTER_INDEPENDENT ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")

# Read the same raw matrix and the frozen tables for comparison. Keep the independent averaging step.
check_tables <- list()
public_data_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)
protein_path <- file.path(
    public_data_root, "write", "31_proteomic_validation", "source_metadata", "PDC000514",
    "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
)
source_table <- as.data.frame(fread(protein_path, check.names = FALSE))
is_protein <- !source_table$Gene %in% c("Mean", "Median", "StdDev") &
    !is.na(source_table$Gene) & source_table$Gene != ""
protein_table <- source_table[is_protein, , drop = FALSE]
protein_dir <- file.path(revision, "step5_protein_common")
member_table <- read.csv(file.path(protein_dir, "protein_fixed_members.csv"))
effect_table <- read.csv(file.path(protein_dir, "effects.csv"))
patient_score_table <- read.csv(file.path(protein_dir, "patient_scores.csv"))

for (analysis in c("unique_peptides", "all_peptides")) {
    prefix <- if (analysis == "unique_peptides") "Unshared Log" else "Log"
    columns <- grep(paste0("^", prefix, " KNCC_GBM[0-9]+_T[12]:"), colnames(protein_table), value = TRUE)
    sample_ids <- sub(paste0("^", prefix, " (KNCC_GBM[0-9]+_T[12]):.*$"), "\\1", columns)
    raw_matrix <- as.matrix(protein_table[, columns, drop = FALSE])
    storage.mode(raw_matrix) <- "double"
    unique_samples <- sort(unique(sample_ids))
    expression_matrix <- matrix(
        NA_real_, nrow = nrow(raw_matrix), ncol = length(unique_samples),
        dimnames = list(protein_table$Gene, unique_samples)
    )

    # Each row is a protein. Average only finite measurements within each sample.
    for (sample_id in unique_samples) {
        sample_matrix <- raw_matrix[, sample_ids == sample_id, drop = FALSE]
        for (i in seq_len(nrow(sample_matrix))) {
            measurements <- sample_matrix[i, ]
            finite_measurements <- measurements[is.finite(measurements)]
            if (!length(finite_measurements)) next

            expression_matrix[i, sample_id] <- mean(finite_measurements)
        }
    }

    # Compare coverage for the archived protein members under the historical eligibility rules.
    genes <- member_table$gene[member_table$analysis == analysis]
    patients <- sort(unique(sub("_T[12]$", "", unique_samples)))
    primary_ids <- paste0(patients, "_T1")
    recurrent_ids <- paste0(patients, "_T2")

    # This combined review retains the original union of observed T1/T2 patient IDs.
    if (!all(c(primary_ids, recurrent_ids) %in% colnames(expression_matrix))) {
        stop("The historical patient union requires matching T1 and T2 columns for every included patient.")
    }
    joint <- is.finite(expression_matrix[genes, primary_ids, drop = FALSE]) &
        is.finite(expression_matrix[genes, recurrent_ids, drop = FALSE])
    fixed_genes <- genes[rowSums(joint) == length(patients)]
    complete_patients <- patients[colSums(joint) == length(genes)]
    reference_fixed <- member_table$gene[member_table$analysis == analysis & member_table$all_105]
    expected_complete <- if (analysis == "unique_peptides") 22L else 29L
    record_check(
        paste(analysis, "independent_coverage"),
        length(patients) == 105L && setequal(fixed_genes, reference_fixed) &&
            length(complete_patients) == expected_complete
    )

    # The original independent check standardizes before the fixed-protein or complete-case restriction.
    scaled_matrix <- t(scale(t(expression_matrix[genes, , drop = FALSE])))
    for (mode in c("fixed_proteins_all105", "original11_complete_patients")) {
        selected_genes <- if (mode == "fixed_proteins_all105") fixed_genes else genes
        selected_patients <- if (mode == "fixed_proteins_all105") patients else complete_patients
        delta <- numeric(length(selected_patients))
        for (i in seq_along(selected_patients)) {
            patient_id <- selected_patients[i]
            primary <- mean(scaled_matrix[selected_genes, paste0(patient_id, "_T1")])
            recurrent <- mean(scaled_matrix[selected_genes, paste0(patient_id, "_T2")])
            delta[i] <- recurrent - primary
        }

        # One saved effect row and the matching patient differences identify this analysis and mode.
        reference_table <- effect_table[effect_table$analysis == analysis & effect_table$mode == mode, ]
        reference_scores_table <- patient_score_table[
            patient_score_table$analysis == analysis & patient_score_table$mode == mode,
        ]
        patient_order <- match(selected_patients, reference_scores_table$patient)
        score_difference <- delta - reference_scores_table$delta[patient_order]
        record_check(
            paste(analysis, mode, "independent_scores_effect"),
            length(selected_patients) == reference_table$n_pairs &&
                max(abs(score_difference)) < 1e-10 && abs(mean(delta) - reference_table$mean_delta) < 1e-10 &&
                abs(t.test(delta)$p.value - reference_table$paired_t_p) < 1e-10 &&
                sum(delta > 0) == reference_table$positive
        )
    }
}
checks_table_total <- do.call(rbind, check_tables)
write.csv(checks_table_total, file.path(protein_dir, "independent_checks.csv"), row.names = FALSE)
cat("PROTEIN_INDEPENDENT ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")

# Read the final rank archive and build a new ranked vector from the original edgeR statistics.
check_tables <- list()
rank_dir <- file.path(revision, "step5_rank_unification_v2")
enrichment_reference_table <- read.csv(file.path(rank_dir, "sqrtF_fgsea.csv"))
rank_reference_table <- read.csv(file.path(rank_dir, "ranked_genes.csv"))
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
if (!file.exists(definition_path)) {
    definition_path <- file.path(revision, "step2_extended_specificity", "locked_gene_sets_used.csv")
}
definition_table <- read.csv(definition_path)
target_genes <- definition_table$gene[definition_table$set_id == "target_top20"]
legacy_figure_dir <- file.path(project, "write", "41_mg_inflammatory_sci_rebuild", "Figure5")
for (dataset in c("GSE121810", "GSE154795")) {
    filename <- if (dataset == "GSE121810") {
        "GSE121810_raw20_edger.csv"
    } else {
        "GSE154795_raw20_edger_auto_strict_idhwt.csv"
    }
    de_table <- read.csv(file.path(legacy_figure_dir, filename))
    rank_score <- sign(de_table$logFC) * sqrt(pmax(de_table$F, 0))
    reference_table <- rank_reference_table[rank_reference_table$dataset == dataset, ]
    enrichment_table <- enrichment_reference_table[enrichment_reference_table$dataset == dataset, ]
    reference_order <- match(de_table$gene, reference_table$gene)
    score_difference <- rank_score - reference_table$rank_sqrtF[reference_order]
    record_check(
        paste(dataset, "independent_rank"),
        setequal(de_table$gene, reference_table$gene) && max(abs(score_difference)) < 1e-12 &&
            sum(de_table$F < 0) == enrichment_table$negative_F_clamped
    )

    # The two fixed cohort seeds preserve the original independent serial fgsea runs.
    rank_order <- order(-rank_score, de_table$gene)
    ranked_genes <- setNames(rank_score[rank_order], de_table$gene[rank_order])
    set.seed(if (dataset == "GSE121810") 20261021L else 20261022L)
    result_table <- fgseaMultilevel(
        list(program = target_genes),
        stats   = ranked_genes,
        minSize = 5,
        maxSize = 500,
        eps     = 0,
        BPPARAM = BiocParallel::SerialParam()
    )
    record_check(
        paste(dataset, "independent_fgsea"),
        nrow(result_table) == 1L && result_table$size == 19L &&
            sign(result_table$NES) == sign(enrichment_table$NES) && result_table$pval < .001,
        paste("NES", result_table$NES, "P", result_table$pval)
    )
}
checks_table_total <- do.call(rbind, check_tables)
write.csv(checks_table_total, file.path(rank_dir, "independent_checks.csv"), row.names = FALSE)
cat("RANK_INDEPENDENT ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
