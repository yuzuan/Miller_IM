#!/usr/bin/env Rscript
# Purpose: Add the frozen human inflammatory Hallmark to the six-set paired CAMERA family.
# Inputs: step2_extended_specificity/{locked_gene_sets_used,camera_results}.csv, step2_locked_gene_sets.csv if present,
#   step1_primary_recurrence/<dataset>_paired_voom.rds, and the counts in provenance/step1_input_manifest.csv.
#   Use archived hallmark_original_members.csv and source_versions.csv, or the verified local MSigDB 2025.1 cache.
# Outputs under step2_hallmark_completion/: hallmark_original_members.csv, hallmark_cohort_members.csv,
#   camera_seven_set_family.csv, hallmark_patient_scores.csv, hallmark_patient_effects.csv, source_versions.csv,
#   validation_checks.csv and sessionInfo.txt.
# Run: Rscript NOA_R1/scripts/02e_hallmark_comparators.R NOA_R1 after 01a and 02d, in a scratch copy for reruns.

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

# Resolve the revision directory and the original Hallmark source exported in the release package.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 02e_hallmark_comparators.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- file.path(revision, "step2_hallmark_completion")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
members_path <- file.path(output_dir, "hallmark_original_members.csv")
versions_path <- file.path(output_dir, "source_versions.csv")

# Store cohort, source, and test checks in the existing validation table. Returns nothing.
checks <- list()

check <- function(label, condition, detail = "") {
    passed <- isTRUE(condition)
    checks[[length(checks) + 1L]] <<- data.frame(check = label, passed = passed, detail = detail)

    if (!passed) {
        validation_table_total <- do.call(rbind, checks)
        write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
        stop(label, " ", detail)
    }
}

# The archived gene order and source version are sufficient for testing the fixed set.
hallmark_set_id <- "HALLMARK_INFLAMMATORY_RESPONSE"
expected_cache_digest <- "d898c396b109d61ef9b20b0632f0bfa7"
if (file.exists(members_path) && file.exists(versions_path)) {
    hallmark_table <- read.csv(members_path, stringsAsFactors = FALSE)
    version_table <- read.csv(versions_path, stringsAsFactors = FALSE)
    check(
        "archived_source_version",
        nrow(version_table) == 1L && identical(version_table$db_version, "2025.1.Hs") &&
            identical(version_table$md5, expected_cache_digest) &&
            identical(hallmark_table$rank, seq_len(nrow(hallmark_table))) &&
            all(hallmark_table$gs_name == hallmark_set_id)
    )
    hallmark_genes <- hallmark_table$gene
    db_version <- unique(hallmark_table$db_version)
} else {

    # To start without the archive, resolve and verify the local cache through the R tools namespace.
    cache_path <- file.path(tools::R_user_dir("msigdbr", "cache"), "msigdb.2025.1.Hs.rds")
    cache_digest <- unname(tools::md5sum(cache_path))
    check(
        "cache_md5",
        file.exists(cache_path) && identical(cache_digest, expected_cache_digest),
        "A verified local MSigDB cache or both archived Hallmark source tables are required"
    )
    package_version <- as.character(packageVersion("msigdbr"))
    if (package_version != "25.1.1") stop("The local MSigDB cache requires msigdbr 25.1.1")
    msigdb_table <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    db_version <- unique(msigdb_table$db_version)
    hallmark_genes <- unique(msigdb_table$gene_symbol[msigdb_table$gs_name == hallmark_set_id])
    hallmark_table <- data.frame(
        rank       = seq_along(hallmark_genes),
        gene       = hallmark_genes,
        gs_name    = hallmark_set_id,
        db_version = db_version
    )
    version_table <- data.frame(
        msigdbr_version = package_version,
        db_version = db_version,
        cache = "<public_data_root>/msigdb.2025.1.Hs.rds",
        md5   = cache_digest,
        note  = "official cached human data; no package installation"
    )
}
check("db_version", identical(db_version, "2025.1.Hs"))
check("hallmark_200", length(hallmark_genes) == 200L && !anyDuplicated(hallmark_genes) &&
    !anyNA(hallmark_genes) && all(nzchar(hallmark_genes)))
write.csv(hallmark_table, members_path, row.names = FALSE, na = "")

# Preserve the six historical set positions. Add the original Hallmark as the last set in the family.
extended_dir <- file.path(revision, "step2_extended_specificity")
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
if (!file.exists(definition_path)) definition_path <- file.path(extended_dir, "locked_gene_sets_used.csv")
definition_table <- read.csv(definition_path, stringsAsFactors = FALSE)
set_ids <- unique(definition_table$set_id)
check("six_historical_sets", length(set_ids) == 6L)
sets <- setNames(vector("list", length(set_ids)), set_ids)
for (set_id in set_ids) sets[[set_id]] <- definition_table$gene[definition_table$set_id == set_id]
sets[[hallmark_set_id]] <- hallmark_genes

# The prior six-set p values and paired design identify the analysis being extended.
extended_camera_table <- read.csv(file.path(extended_dir, "camera_results.csv"), stringsAsFactors = FALSE)
manifest_table <- read.csv(file.path(revision, "provenance", "step1_input_manifest.csv"), stringsAsFactors = FALSE)
public_root <- Sys.getenv("NOA_PUBLIC_DATA_ROOT")
if (!nzchar(public_root)) public_root <- Sys.getenv("MILLER_IM_DATA_ROOT")
if (!nzchar(public_root)) {
    project <- Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision))
    public_root <- file.path(project, "data")
}
public_root <- normalizePath(path.expand(public_root), winslash = "/", mustWork = FALSE)

camera_tables <- list()
coverage_tables <- list()
score_tables <- list()
effect_tables <- list()
for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, ]
    dataset <- entry_table$dataset
    model_path <- file.path(revision, "step1_primary_recurrence", paste0(dataset, "_paired_voom.rds"))
    model <- readRDS(model_path)
    genes <- rownames(model$voom)
    sample_table <- model$metadata
    check(
        paste0(dataset, "_same_formal_input"),
        nrow(sample_table) == 2 * entry_table$expected_pairs &&
            length(unique(sample_table$patient_id)) == entry_table$expected_pairs &&
            length(genes) == entry_table$expected_tested_genes && identical(model$signature, sets$target_top20)
    )

    # Report all original Hallmark members and their overlap with the measured target.
    present <- lapply(sets, intersect, y = genes)
    hallmark_present <- present[[hallmark_set_id]]
    check(paste0(dataset, "_hallmark_coverage"), length(hallmark_present) >= 5L)
    coverage_tables[[i]] <- data.frame(
        dataset          = dataset,
        rank             = seq_along(hallmark_genes),
        gene             = hallmark_genes,
        tested           = hallmark_genes %in% genes,
        overlap_target20 = hallmark_genes %in% sets$target_top20,
        overlap_target19 = hallmark_genes %in% model$present
    )

    # Match the six-set CAMERA arguments and adjust BH over seven p values within each cohort and mode.
    indices <- lapply(present, match, table = genes)
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
            dataset     = dataset,
            set_id      = rownames(camera_fit_table),
            mode        = mode,
            n_pairs     = entry_table$expected_pairs,
            n_genes     = camera_fit_table$NGenes,
            direction   = camera_fit_table$Direction,
            correlation = if (mode == "estimated") camera_fit_table$Correlation else rho,
            p_value     = camera_fit_table$PValue,
            bh_seven    = p.adjust(camera_fit_table$PValue, "BH")
        )
        historical_table <- extended_camera_table[extended_camera_table$dataset == dataset &
            extended_camera_table$mode == mode, ]
        camera_table$historical_bh_six <- historical_table$bh_within_cohort_mode[
            match(camera_table$set_id, historical_table$set_id)
        ]
        historical_rows <- match(historical_table$set_id, camera_table$set_id)
        check(
            paste0(dataset, "_", mode, "_old_six_p_reproduced"),
            nrow(historical_table) == 6L &&
                max(abs(camera_table$p_value[historical_rows] - historical_table$p_value)) < 1e-12
        )
        camera_tables[[length(camera_tables) + 1L]] <- camera_table
    }

    # Hallmark scores follow the unfiltered TMM, prior count 1 and cohort gene-standardization from 02d.
    counts_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$counts_path)
    if (!file.exists(counts_path)) stop("Missing pseudobulk counts: ", counts_path)
    counts_table <- read.csv(gzfile(counts_path), row.names = 1L, check.names = FALSE)
    counts_matrix <- as.matrix(counts_table)
    counts_matrix <- counts_matrix[, as.character(sample_table$sample_id), drop = FALSE]
    rm(counts_table)
    check(
        paste0(dataset, "_raw_sample_match"),
        identical(colnames(counts_matrix), as.character(sample_table$sample_id)) &&
            nrow(counts_matrix) == entry_table$expected_input_genes
    )
    normalized_counts <- calcNormFactors(DGEList(counts = counts_matrix), method = "TMM")
    log_cpm_matrix <- cpm(normalized_counts, log = TRUE, prior.count = 1)
    score_matrix <- log_cpm_matrix[hallmark_present, , drop = FALSE]
    check(paste0(dataset, "_score_sd"), all(apply(score_matrix, 1, sd) > 0))
    scores <- colMeans(t(scale(t(score_matrix))))
    score_table <- data.frame(
        dataset    = dataset,
        set_id     = hallmark_set_id,
        sample_id  = as.character(sample_table$sample_id),
        patient_id = as.character(sample_table$patient_id),
        condition  = as.character(sample_table$condition),
        n_genes    = nrow(score_matrix),
        score      = unname(scores)
    )
    score_tables[[i]] <- score_table

    # Compute the paired t interval on aligned recurrent-minus-primary differences.
    primary_table <- score_table[score_table$condition == "Primary", ]
    recurrent_table <- score_table[score_table$condition == "Recurrent", ]
    check(paste0(dataset, "_score_pairing"), setequal(primary_table$patient_id, recurrent_table$patient_id) &&
        nrow(primary_table) == entry_table$expected_pairs)
    recurrent_rows <- match(primary_table$patient_id, recurrent_table$patient_id)
    delta <- recurrent_table$score[recurrent_rows] - primary_table$score
    paired_test <- t.test(delta)
    effect_tables[[i]] <- data.frame(
        dataset    = dataset,
        set_id     = hallmark_set_id,
        n_pairs    = length(delta),
        n_genes    = nrow(score_matrix),
        positive_pairs = sum(delta > 0),
        mean_delta = mean(delta),
        ci_low     = paired_test$conf.int[1],
        ci_high    = paired_test$conf.int[2],
        paired_t_p = paired_test$p.value
    )
    rm(counts_matrix, normalized_counts, log_cpm_matrix, score_matrix)
    gc(verbose = FALSE)
}

# Keep the archived filenames and inherited source metadata. The archive check identifies source reuse.
camera_table_total <- do.call(rbind, camera_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
score_table_total <- do.call(rbind, score_tables)
effect_table_total <- do.call(rbind, effect_tables)
validation_table_total <- do.call(rbind, checks)
write.csv(camera_table_total, file.path(output_dir, "camera_seven_set_family.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "hallmark_cohort_members.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "hallmark_patient_scores.csv"), row.names = FALSE, na = "")
write.csv(effect_table_total, file.path(output_dir, "hallmark_patient_effects.csv"), row.names = FALSE, na = "")
write.csv(version_table, versions_path, row.names = FALSE, na = "")
write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Validation checks passed: ", nrow(validation_table_total), "\n", sep = "")
