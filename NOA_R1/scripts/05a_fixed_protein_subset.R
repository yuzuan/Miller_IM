#!/usr/bin/env Rscript

# Purpose  Recalculate paired PDC000514 scores for the same proteins in each patient.
# Inputs   target_top20 from step2_locked_gene_sets.csv (set_id and gene).
#          NOA_PUBLIC_DATA_ROOT/write/31_proteomic_validation/source_metadata/PDC000514/
#          KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv (Gene and paired Log/Unshared Log columns).
#          NOA_LEGACY_ROOT/write/54_figure4_proteome_rebuild/Figure4/source_data/
#          Figure4A_pdc_patient_delta_waterfall_source.csv (analysis, patient, Primary, Recurrent).
# Outputs  step5_protein_common/{protein_patient_coverage.csv,protein_fixed_members.csv,
#          analysis_feasibility.csv,patient_scores.csv,effects.csv,checks.csv,sessionInfo.txt}.
# Run      Rscript NOA_R1/scripts/05a_fixed_protein_subset.R <scratch_revision_dir>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_protein_common")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir)

# Preserve the check record and stop at the first failure of the frozen protocol. No return value.
check_tables <- list()
record_check <- function(check_name, passed) {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(check = check_name, passed = isTRUE(passed))

    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
        stop(check_name)
    }
}

# Keep the original protein symbols and measurement columns. The TSV also includes summary rows.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_data_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)
source_path <- file.path(
    public_data_root, "write", "31_proteomic_validation", "source_metadata", "PDC000514",
    "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
)
record_check("original_source_exists", file.exists(source_path))
source_table <- as.data.frame(fread(source_path, check.names = FALSE))
is_protein <- !source_table$Gene %in% c("Mean", "Median", "StdDev") &
    !is.na(source_table$Gene) & source_table$Gene != ""
protein_table <- source_table[is_protein, , drop = FALSE]

# The packaged gene-set table retains the original target genes in author-ranked order.
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
definition_table <- read.csv(definition_path)
target_genes <- definition_table$gene[definition_table$set_id == "target_top20"]
record_check("target_20", length(target_genes) == 20L)

coverage_tables <- list()
member_tables <- list()
score_tables <- list()
effect_tables <- list()
feasibility_tables <- list()

for (analysis in c("unique_peptides", "all_peptides")) {

    # Match the original T1/T2 columns and order the complete patient pairs as in Figure 4.
    prefix <- if (analysis == "unique_peptides") "Unshared Log" else "Log"
    columns <- grep(paste0("^", prefix, " KNCC_GBM[0-9]+_T[12]:"), colnames(protein_table), value = TRUE)
    sample_ids <- sub(paste0("^", prefix, " (KNCC_GBM[0-9]+_T[12]):.*$"), "\\1", columns)
    sample_map_table <- data.frame(
        sample    = sample_ids,
        patient   = sub("_T[12]$", "", sample_ids),
        condition = ifelse(grepl("_T2$", sample_ids), "Recurrent", "Primary")
    )
    sample_map_table <- unique(sample_map_table)
    pair_counts <- table(sample_map_table$patient, sample_map_table$condition)
    patients <- rownames(pair_counts)[pair_counts[, "Primary"] == 1L & pair_counts[, "Recurrent"] == 1L]
    sample_map_table <- sample_map_table[sample_map_table$patient %in% patients, , drop = FALSE]
    sample_map_table <- sample_map_table[order(sample_map_table$patient, sample_map_table$condition), ]
    record_check(paste(analysis, "formal_105_pairs"), length(patients) == 105L && nrow(sample_map_table) == 210L)

    # Each row is a protein. Average duplicate measurements for each retained sample.
    raw_matrix <- as.matrix(protein_table[, columns, drop = FALSE])
    storage.mode(raw_matrix) <- "double"
    expression_matrix <- matrix(
        NA_real_, nrow = nrow(raw_matrix), ncol = nrow(sample_map_table),
        dimnames = list(protein_table$Gene, sample_map_table$sample)
    )
    for (sample_id in sample_map_table$sample) {
        sample_matrix <- raw_matrix[, sample_ids == sample_id, drop = FALSE]
        expression_matrix[, sample_id] <- rowMeans(sample_matrix, na.rm = TRUE)
    }
    expression_matrix[!is.finite(expression_matrix)] <- NA_real_
    present_genes <- intersect(target_genes, rownames(expression_matrix))
    record_check(paste(analysis, "eleven_original"), length(present_genes) == 11L)

    # Standardize each of the 11 measured proteins across all 210 samples before selecting subsets.
    scaled_matrix <- t(scale(t(expression_matrix[present_genes, , drop = FALSE])))
    scaled_matrix[!is.finite(scaled_matrix)] <- NA_real_

    if (analysis == "unique_peptides") {

        # The original per-patient score uses any three or more proteins measured on both sides.
        original_score_tables <- vector("list", length(patients))
        for (i in seq_along(patients)) {
            patient_id <- patients[i]
            primary <- scaled_matrix[, paste0(patient_id, "_T1")]
            recurrent <- scaled_matrix[, paste0(patient_id, "_T2")]
            joint <- is.finite(primary) & is.finite(recurrent)
            if (sum(joint) < 3L) next

            original_score_tables[[i]] <- data.frame(
                patient   = patient_id,
                Primary   = mean(primary[joint]),
                Recurrent = mean(recurrent[joint]),
                n_common  = sum(joint)
            )
        }
        original_score_table_total <- do.call(rbind, original_score_tables)
        legacy_score_path <- file.path(
            project, "write", "54_figure4_proteome_rebuild", "Figure4", "source_data",
            "Figure4A_pdc_patient_delta_waterfall_source.csv"
        )
        legacy_score_table <- read.csv(legacy_score_path)
        legacy_score_table <- legacy_score_table[legacy_score_table$analysis == "unique_peptides", ]
        score_order <- match(original_score_table_total$patient, legacy_score_table$patient)
        primary_difference <- original_score_table_total$Primary - legacy_score_table$Primary[score_order]
        recurrent_difference <- original_score_table_total$Recurrent - legacy_score_table$Recurrent[score_order]
        record_check(
            "original_105_score_reproduced",
            nrow(original_score_table_total) == 105L && all(!is.na(score_order)) &&
                max(abs(primary_difference)) < 1e-10 && max(abs(recurrent_difference)) < 1e-10
        )
    }

    # A finite pair requires both samples. Rows give protein coverage and columns give patient coverage.
    coverage_matrix <- is.finite(expression_matrix[present_genes, paste0(patients, "_T1"), drop = FALSE]) &
        is.finite(expression_matrix[present_genes, paste0(patients, "_T2"), drop = FALSE])
    colnames(coverage_matrix) <- patients
    finite_pairs <- rowSums(coverage_matrix)
    fixed_genes <- present_genes[finite_pairs == length(patients)]
    complete_patients <- patients[colSums(coverage_matrix) == length(present_genes)]
    coverage_tables[[analysis]] <- data.frame(
        analysis    = analysis,
        gene        = rep(present_genes, each = length(patients)),
        patient     = rep(patients, times = length(present_genes)),
        both_finite = as.vector(t(coverage_matrix))
    )
    member_tables[[analysis]] <- data.frame(
        analysis = analysis,
        gene     = present_genes,
        finite_pairs = finite_pairs,
        all_105  = present_genes %in% fixed_genes
    )
    feasibility_tables[[analysis]] <- data.frame(
        analysis           = analysis,
        original_proteins  = length(present_genes),
        fixed_all_105      = length(fixed_genes),
        complete_11_pairs  = length(complete_patients),
        primary_feasible   = length(fixed_genes) >= 5L,
        secondary_feasible = length(complete_patients) >= 10L
    )

    # Estimate one effect per eligibility rule. The protocol requires five proteins or ten complete pairs.
    for (mode in c("fixed_proteins_all105", "original11_complete_patients")) {
        if (mode == "fixed_proteins_all105" && length(fixed_genes) < 5L) next
        if (mode == "original11_complete_patients" && length(complete_patients) < 10L) next

        selected_genes <- if (mode == "fixed_proteins_all105") fixed_genes else present_genes
        selected_patients <- if (mode == "fixed_proteins_all105") patients else complete_patients
        patient_score_tables <- vector("list", length(selected_patients))

        # Each row stores both scores and the recurrent-minus-primary difference for one patient.
        for (i in seq_along(selected_patients)) {
            patient_id <- selected_patients[i]
            primary <- scaled_matrix[selected_genes, paste0(patient_id, "_T1")]
            recurrent <- scaled_matrix[selected_genes, paste0(patient_id, "_T2")]
            record_check(paste(analysis, mode, patient_id, "complete"), all(is.finite(primary) & is.finite(recurrent)))
            primary_score <- mean(primary)
            recurrent_score <- mean(recurrent)
            patient_score_tables[[i]] <- data.frame(
                analysis   = analysis,
                mode       = mode,
                patient    = patient_id,
                n_proteins = length(selected_genes),
                Primary    = primary_score,
                Recurrent  = recurrent_score,
                delta      = recurrent_score - primary_score
            )
        }
        patient_score_table_total <- do.call(rbind, patient_score_tables)
        score_tables[[length(score_tables) + 1L]] <- patient_score_table_total

        # Apply a paired t test and the original two-sided sign-flip test to the same patient differences.
        delta <- patient_score_table_total$delta
        paired_test <- t.test(delta)
        n_pairs <- length(delta)
        observed <- abs(mean(delta))
        if (n_pairs <= 20L) {
            signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n_pairs)))
            null_means <- as.vector(signs %*% delta) / n_pairs
            sign_flip_p <- mean(abs(null_means) >= observed - 1e-12)
        } else {

            # Each eligibility rule restarts its analysis-specific stream. Moving the seed changes the p values.
            set.seed(if (analysis == "unique_peptides") 20260713L else 20260714L)
            simulations <- 100000L
            exceedances <- 0L
            for (start in seq.int(1L, simulations, by = 5000L)) {
                batch_size <- min(5000L, simulations - start + 1L)
                signs <- matrix(sample(c(-1, 1), batch_size * n_pairs, replace = TRUE), nrow = batch_size)
                differences <- matrix(delta, nrow = batch_size, ncol = n_pairs, byrow = TRUE)
                exceedances <- exceedances + sum(abs(rowMeans(signs * differences)) >= observed - 1e-12)
            }
            sign_flip_p <- (exceedances + 1) / (simulations + 1)
        }
        effect_tables[[length(effect_tables) + 1L]] <- data.frame(
            analysis    = analysis,
            mode        = mode,
            n_pairs     = n_pairs,
            n_proteins  = length(selected_genes),
            proteins    = paste(selected_genes, collapse = ";"),
            positive    = sum(delta > 0),
            mean_delta  = mean(delta),
            ci_low      = paired_test$conf.int[1],
            ci_high     = paired_test$conf.int[2],
            paired_t_p  = paired_test$p.value,
            sign_flip_p = sign_flip_p
        )
    }
    rm(raw_matrix, expression_matrix, scaled_matrix, sample_matrix)
    gc(verbose = FALSE)
}

# Concatenate in analysis order, then eligibility order, with the archived column names.
coverage_table_total <- do.call(rbind, coverage_tables)
member_table_total <- do.call(rbind, member_tables)
feasibility_table_total <- do.call(rbind, unname(feasibility_tables))
write.csv(coverage_table_total, file.path(output_dir, "protein_patient_coverage.csv"), row.names = FALSE, na = "")
write.csv(member_table_total, file.path(output_dir, "protein_fixed_members.csv"), row.names = FALSE, na = "")
write.csv(feasibility_table_total, file.path(output_dir, "analysis_feasibility.csv"), row.names = FALSE, na = "")

# The frozen feasibility rules can leave the score and effect outputs empty.
if (length(score_tables)) {
    score_table_total <- do.call(rbind, score_tables)
    write.csv(score_table_total, file.path(output_dir, "patient_scores.csv"), row.names = FALSE, na = "")
}
if (length(effect_tables)) {
    effect_table_total <- do.call(rbind, effect_tables)
    write.csv(effect_table_total, file.path(output_dir, "effects.csv"), row.names = FALSE, na = "")
}
checks_table_total <- do.call(rbind, check_tables)
write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

cat("PROTEIN_COMMON ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
