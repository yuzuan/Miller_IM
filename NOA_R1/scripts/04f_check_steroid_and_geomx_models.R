#!/usr/bin/env Rscript
# Purpose: Independently reproduce the original Hallmark, steroid, clinical crosswalk and GeoMx checks.
# Inputs: step2_locked_gene_sets.csv, provenance/step1_input_manifest.csv (or a prepared root manifest),
# step1_primary_recurrence/<dataset>_paired_voom.rds and the pan-myeloid counts listed in the manifest.
# Read step2_hallmark_completion/{camera_seven_set_family,hallmark_patient_scores,hallmark_original_members}.csv
# and step4_steroids_geomx/{camera_results,sample_scores,corticosteroid_original_members,
# GSE274546_steroid_patient_status,GeoMx_subgroup_interaction}.csv.
# External inputs: msigdbr 25.1.1 human 2025.1 cache, the locally held GBM-CARE workbook under
# NOA_PUBLIC_DATA_ROOT/reference_metadata/GBM_CARE_WT/spitzer_supptable1.xlsx, and NOA_LEGACY_ROOT/
# write/55_figure3_visual_polish/Figure3/source_data/Figure3A_geomx_paired_dumbbell_source.csv.
# Outputs: independent_checks.csv in step2_hallmark_completion and step4_steroids_geomx.
# Run: Rscript NOA_R1/scripts/04f_check_steroid_and_geomx_models.R NOA_R1

suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
    library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04f_check_steroid_and_geomx_models.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
hallmark_dir <- file.path(revision, "step2_hallmark_completion")
steroid_dir <- file.path(revision, "step4_steroids_geomx")
hallmark_output <- file.path(hallmark_dir, "independent_checks.csv")
steroid_output <- file.path(steroid_dir, "independent_checks.csv")
if (any(file.exists(c(hallmark_output, steroid_output)))) stop("Independent review output already exists")

# Confirm the frozen package has its local human cache before making either MSigDB call.
cache_path <- file.path(tools::R_user_dir("msigdbr", "cache"), "msigdb.2025.1.Hs.rds")
if (!file.exists(cache_path)) stop("Missing frozen MSigDB cache. ", cache_path)
if (as.character(packageVersion("msigdbr")) != "25.1.1") stop("The frozen local cache requires msigdbr 25.1.1")
hallmark_msigdb_table <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))
steroid_msigdb_table <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
hallmark_rows <- hallmark_msigdb_table$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE"
hallmark_genes <- unique(hallmark_msigdb_table$gene_symbol[hallmark_rows])
steroid_rows <- steroid_msigdb_table$gs_name == "GOBP_RESPONSE_TO_CORTICOSTEROID"
steroid_genes <- unique(steroid_msigdb_table$gene_symbol[steroid_rows])

# Compare package membership with each original source table in its recorded order.
hallmark_members_table <- read.csv(file.path(hallmark_dir, "hallmark_original_members.csv"))
steroid_members_table <- read.csv(file.path(steroid_dir, "corticosteroid_original_members.csv"))
checks <- list(hallmark = list(), steroid = list())

# Collect each independent check in order, including failures. Returns nothing.
add_check <- function(section, name, value, detail = "") {
    check_table <- data.frame(check = name, passed = isTRUE(value), detail = detail)
    rows <- checks[[section]]
    rows[[length(rows) + 1L]] <- check_table
    checks[[section]] <<- rows
}

add_check(
    "hallmark", "source_full_200", length(hallmark_genes) == 200L &&
        identical(hallmark_genes, hallmark_members_table$gene)
)
add_check(
    "steroid", "source_full_153", length(steroid_genes) == 153L &&
        identical(steroid_genes, steroid_members_table$gene)
)

# Read reference model outputs and the three locked sets independently of the main model scripts.
hallmark_camera_table <- read.csv(file.path(hallmark_dir, "camera_seven_set_family.csv"))
steroid_camera_table <- read.csv(file.path(steroid_dir, "camera_results.csv"))
hallmark_score_table <- read.csv(file.path(hallmark_dir, "hallmark_patient_scores.csv"))
steroid_score_table <- read.csv(file.path(steroid_dir, "sample_scores.csv"))
definitions_table <- read.csv(file.path(revision, "step2_locked_gene_sets.csv"))
target <- definitions_table$gene[definitions_table$set_id == "target_top20"]
sets <- list(
    HALLMARK_INFLAMMATORY_RESPONSE  = hallmark_genes,
    target_minus_SGK1_PDK4 = setdiff(target, c("SGK1", "PDK4")),
    GOBP_RESPONSE_TO_CORTICOSTEROID = steroid_genes
)

manifest_path <- file.path(revision, "step1_input_manifest.csv")
if (!file.exists(manifest_path)) manifest_path <- file.path(revision, "provenance", "step1_input_manifest.csv")
manifest_table <- read.csv(manifest_path)
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)
for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, , drop = FALSE]
    dataset <- entry_table$dataset
    model <- readRDS(file.path(revision, "step1_primary_recurrence", paste0(dataset, "_paired_voom.rds")))

    # Call CAMERA for a single set at a time on the same design and fitted gene universe.
    for (set_id in names(sets)) {
        hallmark <- set_id == "HALLMARK_INFLAMMATORY_RESPONSE"
        section <- if (hallmark) "hallmark" else "steroid"
        source_table <- if (hallmark) hallmark_camera_table else steroid_camera_table
        indices <- which(rownames(model$voom) %in% sets[[set_id]])
        for (mode in c("estimated", "fixed_0.01")) {
            rho <- if (mode == "estimated") NA_real_ else .01
            camera_table <- camera(
                model$voom, list(independent = indices),
                design   = model$design,
                contrast = model$coefficient,
                inter.gene.cor = rho,
                allow.neg.cor  = FALSE,
                use.ranks      = FALSE,
                directional    = TRUE,
                trend.var      = FALSE,
                sort           = FALSE
            )
            selected <- source_table$dataset == dataset & source_table$set_id == set_id & source_table$mode == mode
            reference_table <- source_table[selected, , drop = FALSE]
            result_p <- camera_table$PValue
            matches <- nrow(reference_table) == 1L && camera_table$NGenes == reference_table$n_genes &&
                identical(as.character(camera_table$Direction), reference_table$direction) &&
                abs(result_p - reference_table$p_value) < 1e-12
            check_name <- paste(dataset, set_id, mode, "independent_camera", sep = "_")
            add_check(section, check_name, matches, paste("independent P", format(result_p, digits = 14)))
        }
    }

    # Independently calculate per-sample means of the standardized target-set logCPM.
    counts_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$counts_path)
    count_table <- read.csv(gzfile(counts_path), row.names = 1L, check.names = FALSE)
    sample_ids <- as.character(model$metadata$sample_id)
    count_matrix <- as.matrix(count_table)[, sample_ids, drop = FALSE]
    counts <- DGEList(counts = count_matrix)
    counts <- calcNormFactors(counts, method = "TMM")
    log_cpm <- cpm(counts, log = TRUE, prior.count = 1)
    for (set_id in names(sets)) {
        tested_genes <- intersect(sets[[set_id]], rownames(model$voom))
        values <- log_cpm[tested_genes, , drop = FALSE]
        centered <- sweep(values, 1, rowMeans(values), "-")
        scaled <- sweep(centered, 1, apply(values, 1, sd), "/")
        scores <- colMeans(scaled)

        hallmark <- set_id == "HALLMARK_INFLAMMATORY_RESPONSE"
        section <- if (hallmark) "hallmark" else "steroid"
        source_table <- if (hallmark) hallmark_score_table else steroid_score_table
        selected <- source_table$dataset == dataset & source_table$set_id == set_id
        reference_table <- source_table[selected, , drop = FALSE]
        indices <- match(sample_ids, reference_table$sample_id)
        matches <- length(indices) == ncol(values) && all(!is.na(indices)) &&
            max(abs(scores - reference_table$score[indices])) < 1e-10
        add_check(section, paste(dataset, set_id, "independent_score", sep = "_"), matches)
    }
    rm(count_table, count_matrix, counts, log_cpm)
    gc(verbose = FALSE)
}

# Verify the original seven-set and three-set BH families separately for each cohort and correlation mode.
add_check(
    "hallmark", "seven_family_has_28_rows", nrow(hallmark_camera_table) == 28L &&
        all(table(hallmark_camera_table$dataset, hallmark_camera_table$mode) == 7L)
)
for (dataset in unique(hallmark_camera_table$dataset)) {
    for (mode in unique(hallmark_camera_table$mode)) {
        selected <- hallmark_camera_table$dataset == dataset & hallmark_camera_table$mode == mode
        family_table <- hallmark_camera_table[selected, , drop = FALSE]
        error <- max(abs(p.adjust(family_table$p_value, "BH") - family_table$bh_seven))
        add_check("hallmark", paste(dataset, mode, "BH7", sep = "_"), error < 1e-12)
    }
}

add_check(
    "steroid", "step4_family_has_12_rows", nrow(steroid_camera_table) == 12L &&
        all(table(steroid_camera_table$dataset, steroid_camera_table$mode) == 3L)
)
for (dataset in unique(steroid_camera_table$dataset)) {
    for (mode in unique(steroid_camera_table$mode)) {
        selected <- steroid_camera_table$dataset == dataset & steroid_camera_table$mode == mode
        family_table <- steroid_camera_table[selected, , drop = FALSE]
        error <- max(abs(p.adjust(family_table$p_value, "BH") - family_table$bh_three))
        add_check("steroid", paste(dataset, mode, "BH3", sep = "_"), error < 1e-12)
    }
}

# Verify the baseline crosswalk against the source workbook, before the separate case correction.
workbook_path <- file.path(public_root, "reference_metadata", "GBM_CARE_WT", "spitzer_supptable1.xlsx")
clinical_table <- as.data.frame(readxl::read_excel(workbook_path))
clinical_table$prefix <- sub("_[012]$", "", clinical_table$`Tumor ID`)
status_table <- read.csv(file.path(steroid_dir, "GSE274546_steroid_patient_status.csv"))
for (i in seq_len(nrow(status_table))) {
    prefix <- status_table$tumor_prefix[i]
    selected <- clinical_table$prefix == prefix & clinical_table$`Time point` %in% c("T1", "T2")
    patient_table <- clinical_table[selected, , drop = FALSE]
    matches <- nrow(patient_table) == 2L && setequal(patient_table$`Time point`, c("T1", "T2")) &&
        identical(
            as.character(patient_table$`Patient ID`[1]),
            as.character(status_table$clinical_patient_id[i])
        )
    add_check("steroid", paste("clinical", status_table$patient_id[i]), matches)
}
complete <- nrow(status_table) == 45L && sum(status_table$group != "Missing") == 36L
add_check("steroid", "clinical_45_complete_36", complete)

# Read the original Figure 3A source and reproduce the paired and between-trial calculations.
geomx_path <- file.path(
    project, "write", "55_figure3_visual_polish", "Figure3", "source_data", "Figure3A_geomx_paired_dumbbell_source.csv"
)
geomx_table <- read.csv(geomx_path)
geomx_stats_table <- read.csv(file.path(steroid_dir, "GeoMx_subgroup_interaction.csv"))
formal <- nrow(geomx_table) == 22L && !anyDuplicated(geomx_table$Patient_No) &&
    max(abs(geomx_table$delta - geomx_table$Recurrence + geomx_table$Primary)) < 1e-12
add_check("steroid", "geomx_22_unique_and_delta", formal)
for (group in c("Control", "Nivolumab")) {
    differences <- geomx_table$delta[geomx_table$trial_setting == group]
    paired_test <- t.test(differences)
    selected <- geomx_stats_table$contrast == paste0(group, "_paired_delta")
    reference_table <- geomx_stats_table[selected, , drop = FALSE]
    interval <- c(reference_table$ci_low, reference_table$ci_high)
    matches <- nrow(reference_table) == 1L && reference_table$n == length(differences) &&
        reference_table$positive == sum(differences > 0) &&
        abs(reference_table$estimate - mean(differences)) < 1e-12 &&
        abs(reference_table$p_value - paired_test$p.value) < 1e-12 &&
        max(abs(interval - paired_test$conf.int)) < 1e-12
    add_check("steroid", paste("geomx", group), matches)
}

nivolumab <- geomx_table$delta[geomx_table$trial_setting == "Nivolumab"]
control <- geomx_table$delta[geomx_table$trial_setting == "Control"]
welch_test <- t.test(nivolumab, control, var.equal = FALSE)
selected <- geomx_stats_table$contrast == "Nivolumab_minus_Control_interaction"
reference_table <- geomx_stats_table[selected, , drop = FALSE]
interval <- c(reference_table$ci_low, reference_table$ci_high)
matches <- nrow(reference_table) == 1L && abs(reference_table$estimate - (mean(nivolumab) - mean(control))) < 1e-12 &&
    abs(reference_table$p_value - welch_test$p.value) < 1e-12 &&
    max(abs(interval - welch_test$conf.int)) < 1e-12
add_check("steroid", "geomx_Welch_interaction", matches)

# Keep both archived audit tables and return a nonzero status for any failed check.
hallmark_check_table_total <- do.call(rbind, checks$hallmark)
steroid_check_table_total <- do.call(rbind, checks$steroid)
write.csv(hallmark_check_table_total, hallmark_output, row.names = FALSE)
write.csv(steroid_check_table_total, steroid_output, row.names = FALSE)
failed_hallmark <- sum(!hallmark_check_table_total$passed)
failed_steroid <- sum(!steroid_check_table_total$passed)
cat(
    "Hallmark checks passed: ", nrow(hallmark_check_table_total) - failed_hallmark,
    ", failed: ", failed_hallmark, "\n",
    sep = ""
)
cat(
    "Steroid and GeoMx checks passed: ", nrow(steroid_check_table_total) - failed_steroid,
    ", failed: ", failed_steroid, "\n",
    sep = ""
)
if (failed_hallmark + failed_steroid > 0L) quit(status = 1)
