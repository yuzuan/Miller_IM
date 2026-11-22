#!/usr/bin/env Rscript
# Purpose: Reproduce the three steroid-related gene-set comparisons and the formal GeoMx subgroup analysis.
# Inputs: step4_analysis_plan.md, step2_locked_gene_sets.csv, provenance/step1_input_manifest.csv
# (or a prepared root manifest), and step1_primary_recurrence/{<dataset>_paired_voom.rds,camera_results.csv,
# sample_audit.csv}. The manifest lists pan-myeloid count files under NOA_PUBLIC_DATA_ROOT.
# External inputs: msigdbr 25.1.1 human cache msigdb.2025.1.Hs.rds; NOA_PUBLIC_DATA_ROOT/
# reference_metadata/GBM_CARE_WT/spitzer_supptable1.xlsx; and NOA_LEGACY_ROOT/
# write/55_figure3_visual_polish/Figure3/source_data/Figure3A_geomx_paired_dumbbell_source.csv.
# Outputs in a fresh step4_steroids_geomx: corticosteroid_original_members.csv, locked_members_all_sets.csv,
# camera_results.csv, gene_coverage.csv, sample_scores.csv, patient_deltas.csv, patient_effects.csv,
# source_versions.csv, GSE274546_steroid_patient_status.csv, GSE274546_steroid_group_summary.csv,
# GeoMx_formal_patient_source_copy.csv, GeoMx_subgroup_interaction.csv, validation_checks.csv, sessionInfo.txt.
# Run: Rscript NOA_R1/scripts/04a_steroid_and_geomx_models.R NOA_R1

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(msigdbr)
    library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04a_steroid_and_geomx_models.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- file.path(revision, "step4_steroids_geomx")
if (dir.exists(output_dir)) stop("Steroid and GeoMx output directory already exists")
dir.create(output_dir)

# Check the recorded analysis time before using the frozen cached gene membership.
plan <- readLines(file.path(revision, "step4_analysis_plan.md"), warn = FALSE)
check_rows <- list()

# Retain the original audit row and stop with the accumulated checks on failure. Returns nothing.
check <- function(name, ok, detail = "") {
    passed <- isTRUE(ok)
    check_rows[[length(check_rows) + 1L]] <<- data.frame(check = name, passed = passed, detail = detail)
    if (!passed) {
        check_table_total <- do.call(rbind, check_rows)
        write.csv(check_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
        stop(name, " ", detail)
    }
}

check("plan_locked", any(grepl("2026-10-20 00:31:53", plan, fixed = TRUE)))
cache_path <- file.path(tools::R_user_dir("msigdbr", "cache"), "msigdb.2025.1.Hs.rds")
cache_md5 <- if (file.exists(cache_path)) unname(tools::md5sum(cache_path)) else NA_character_
check("cache_md5", identical(cache_md5, "d898c396b109d61ef9b20b0632f0bfa7"))
msigdbr_version <- as.character(packageVersion("msigdbr"))
if (msigdbr_version != "25.1.1") stop("The frozen local cache requires msigdbr 25.1.1")
msigdb_table <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
db_version <- unique(msigdb_table$db_version)
check("source_version", identical(db_version, "2025.1.Hs") && msigdbr_version == "25.1.1")
check("first_choice_absent", !any(msigdb_table$gs_name == "GOBP_RESPONSE_TO_GLUCOCORTICOID"))

# The written plan uses the corticosteroid-response set when the first-choice name is absent.
steroid_genes <- unique(msigdb_table$gene_symbol[msigdb_table$gs_name == "GOBP_RESPONSE_TO_CORTICOSTEROID"])
check("locked_set_153", length(steroid_genes) == 153L && !anyDuplicated(steroid_genes))
steroid_members_table <- data.frame(
    rank       = seq_along(steroid_genes),
    gene       = steroid_genes,
    gs_name    = "GOBP_RESPONSE_TO_CORTICOSTEROID",
    db_version = db_version
)
write.csv(
    steroid_members_table, file.path(output_dir, "corticosteroid_original_members.csv"),
    row.names = FALSE,
    na = ""
)

definitions_table <- read.csv(file.path(revision, "step2_locked_gene_sets.csv"), stringsAsFactors = FALSE)
target <- definitions_table$gene[definitions_table$set_id == "target_top20"]
check("target_20", length(target) == 20L && all(c("SGK1", "PDK4") %in% target))
minus <- setdiff(target, c("SGK1", "PDK4"))
check("target_minus2_18", length(minus) == 18L)
sets <- list(
    target_top20 = target,
    target_minus_SGK1_PDK4 = minus,
    GOBP_RESPONSE_TO_CORTICOSTEROID = steroid_genes
)
members_table_total <- data.frame(set_id = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE))
write.csv(members_table_total, file.path(output_dir, "locked_members_all_sets.csv"), row.names = FALSE, na = "")

# Reuse the Step 1 fitted universe, weighted expression and within-patient design for each cohort.
manifest_path <- file.path(revision, "step1_input_manifest.csv")
if (!file.exists(manifest_path)) manifest_path <- file.path(revision, "provenance", "step1_input_manifest.csv")
manifest_table <- read.csv(manifest_path, stringsAsFactors = FALSE)
model_dir <- file.path(revision, "step1_primary_recurrence")
reference_camera_table <- read.csv(file.path(model_dir, "camera_results.csv"), stringsAsFactors = FALSE)
reference_scores_table <- read.csv(file.path(model_dir, "sample_audit.csv"), stringsAsFactors = FALSE)
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_root <- Sys.getenv(
    "NOA_PUBLIC_DATA_ROOT",
    unset = Sys.getenv("MILLER_IM_DATA_ROOT", unset = file.path(project, "data"))
)

camera_tables <- list()
coverage_tables <- list()
sample_tables <- list()
delta_tables <- list()
effect_tables <- list()
for (i in seq_len(nrow(manifest_table))) {
    entry_table <- manifest_table[i, , drop = FALSE]
    dataset <- entry_table$dataset
    model <- readRDS(file.path(model_dir, paste0(dataset, "_paired_voom.rds")))
    metadata_table <- model$metadata
    genes <- rownames(model$voom)
    expected_pairs <- entry_table$expected_pairs
    formal_inputs <- nrow(metadata_table) == 2 * expected_pairs &&
        length(unique(metadata_table$patient_id)) == expected_pairs &&
        length(genes) == entry_table$expected_tested_genes && identical(model$signature, target)
    check(paste0(dataset, "_formal_inputs"), formal_inputs)

    # Record which members are tested in the frozen genome-wide background for this cohort.
    present <- list()
    for (set_id in names(sets)) {
        members <- sets[[set_id]]
        present[[set_id]] <- intersect(members, genes)
        coverage_tables[[length(coverage_tables) + 1L]] <- data.frame(
            dataset = dataset,
            set_id  = set_id,
            rank    = seq_along(members),
            gene    = members,
            tested  = members %in% genes,
            overlap_target20 = members %in% target,
            overlap_target19 = members %in% model$present
        )
    }
    target_coverage <- length(present$target_top20) == 19L && length(present$target_minus_SGK1_PDK4) == 17L
    check(paste0(dataset, "_target_coverage"), target_coverage)
    check(paste0(dataset, "_steroid_coverage"), length(present$GOBP_RESPONSE_TO_CORTICOSTEROID) >= 5L)

    # Fit the same three-set family with estimated or fixed 0.01 inter-gene correlation.
    indices <- lapply(present, match, table = genes)
    for (mode in c("estimated", "fixed_0.01")) {
        rho <- if (mode == "estimated") NA_real_ else .01
        fit_table <- camera(
            model$voom, indices,
            design   = model$design,
            contrast = model$coefficient,
            inter.gene.cor = rho,
            allow.neg.cor  = FALSE,
            use.ranks      = FALSE,
            directional    = TRUE,
            trend.var      = FALSE,
            sort           = FALSE
        )
        camera_table <- data.frame(
            dataset     = dataset,
            set_id      = rownames(fit_table),
            mode        = mode,
            n_pairs     = expected_pairs,
            n_genes     = fit_table$NGenes,
            direction   = fit_table$Direction,
            correlation = if (mode == "estimated") fit_table$Correlation else rho,
            p_value     = fit_table$PValue,
            bh_three    = p.adjust(fit_table$PValue, "BH")
        )
        reference_rows <- reference_camera_table$dataset == dataset &
            reference_camera_table$correlation_mode == mode
        reference_table <- reference_camera_table[reference_rows, , drop = FALSE]
        target_p <- camera_table$p_value[camera_table$set_id == "target_top20"]
        reproduces <- nrow(reference_table) == 1L && abs(target_p - reference_table$p_value) < 1e-12
        check(paste0(dataset, "_", mode, "_full_target_reproduces"), reproduces)
        camera_tables[[length(camera_tables) + 1L]] <- camera_table
    }

    # Normalize the unfiltered formal counts before z-scoring each set's tested genes across samples.
    counts_path <- sub("^<public_data_root>/", paste0(public_root, "/"), entry_table$counts_path)
    raw_table <- read.csv(gzfile(counts_path), row.names = 1L, check.names = FALSE)
    sample_ids <- as.character(metadata_table$sample_id)
    raw_matrix <- as.matrix(raw_table)[, sample_ids, drop = FALSE]
    correct_samples <- identical(colnames(raw_matrix), sample_ids) &&
        nrow(raw_matrix) == entry_table$expected_input_genes
    check(paste0(dataset, "_raw_samples"), correct_samples)
    counts <- DGEList(counts = raw_matrix)
    counts <- calcNormFactors(counts, method = "TMM")
    log_cpm <- cpm(counts, log = TRUE, prior.count = 1)

    for (set_id in names(sets)) {
        values <- log_cpm[present[[set_id]], , drop = FALSE]
        check(paste0(dataset, "_", set_id, "_score_sd"), all(apply(values, 1, sd) > 0))
        scores <- colMeans(t(scale(t(values))))
        score_table <- data.frame(
            dataset    = dataset,
            set_id     = set_id,
            sample_id  = sample_ids,
            patient_id = as.character(metadata_table$patient_id),
            condition  = as.character(metadata_table$condition),
            n_genes    = nrow(values),
            score      = unname(scores)
        )
        if (set_id == "target_top20") {
            reference_table <- reference_scores_table[reference_scores_table$dataset == dataset, , drop = FALSE]
            matched_scores <- reference_table$score[match(score_table$sample_id, reference_table$sample_id)]
            check(paste0(dataset, "_full_score_reproduces"), max(abs(score_table$score - matched_scores)) < 1e-10)
        }
        sample_tables[[length(sample_tables) + 1L]] <- score_table

        # Compute one recurrence-minus-primary contrast for each formal patient.
        primary_table <- score_table[score_table$condition == "Primary", , drop = FALSE]
        recurrent_table <- score_table[score_table$condition == "Recurrent", , drop = FALSE]
        paired <- setequal(primary_table$patient_id, recurrent_table$patient_id) &&
            nrow(primary_table) == expected_pairs
        check(paste0(dataset, "_", set_id, "_pairs"), paired)
        recurrent_scores <- recurrent_table$score[match(primary_table$patient_id, recurrent_table$patient_id)]
        differences <- recurrent_scores - primary_table$score
        delta_tables[[length(delta_tables) + 1L]] <- data.frame(
            dataset         = dataset,
            set_id          = set_id,
            patient_id      = primary_table$patient_id,
            primary_score   = primary_table$score,
            recurrent_score = recurrent_scores,
            delta           = differences
        )
        paired_test <- t.test(differences)
        effect_tables[[length(effect_tables) + 1L]] <- data.frame(
            dataset     = dataset,
            set_id      = set_id,
            n_pairs     = length(differences),
            n_genes     = nrow(values),
            positive_pairs = sum(differences > 0),
            mean_delta  = mean(differences),
            ci_low      = paired_test$conf.int[1],
            ci_high     = paired_test$conf.int[2],
            paired_t_p  = paired_test$p.value
        )
    }
    rm(raw_table, raw_matrix, counts, log_cpm)
    gc(verbose = FALSE)
}

# Save the same cohort, sample and patient tables plus the original cache provenance.
camera_table_total <- do.call(rbind, camera_tables)
coverage_table_total <- do.call(rbind, coverage_tables)
score_table_total <- do.call(rbind, sample_tables)
delta_table_total <- do.call(rbind, delta_tables)
effect_table_total <- do.call(rbind, effect_tables)
write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(coverage_table_total, file.path(output_dir, "gene_coverage.csv"), row.names = FALSE, na = "")
write.csv(score_table_total, file.path(output_dir, "sample_scores.csv"), row.names = FALSE, na = "")
write.csv(delta_table_total, file.path(output_dir, "patient_deltas.csv"), row.names = FALSE, na = "")
write.csv(effect_table_total, file.path(output_dir, "patient_effects.csv"), row.names = FALSE, na = "")
source_table <- data.frame(
    msigdbr_version = msigdbr_version,
    db_version = db_version,
    cache      = paste0("<public_data_root>/", basename(cache_path)),
    md5        = cache_md5
)
write.csv(source_table, file.path(output_dir, "source_versions.csv"), row.names = FALSE, na = "")

# Match the exact Tumor ID prefix and both time points in the locally held GBM-CARE workbook.
workbook_path <- file.path(public_root, "reference_metadata", "GBM_CARE_WT", "spitzer_supptable1.xlsx")
clinical_table <- as.data.frame(read_excel(workbook_path))
clinical_table$prefix <- sub("_[012]$", "", clinical_table$`Tumor ID`)
clinical_table$time <- trimws(as.character(clinical_table$`Time point`))
clinical_table$steroid <- trimws(as.character(clinical_table$`Steroid treatment before the surgery`))
clinical_table$steroid[!clinical_table$steroid %in% c("Yes", "No")] <- NA_character_
clinical_table <- clinical_table[clinical_table$time %in% c("T1", "T2"), , drop = FALSE]
unique_time <- !anyDuplicated(paste(clinical_table$prefix, clinical_table$time)) &&
    all(c("T1", "T2") %in% clinical_table$time)
check("clinical_prefix_time_unique", unique_time)

# Keep the original case-sensitive extraction for the baseline audited by the later correction.
target_rows <- delta_table_total$dataset == "GSE274546" & delta_table_total$set_id == "target_top20"
target_delta_table <- delta_table_total[target_rows, , drop = FALSE]
prefixes <- sub("^GSE274546_", "", target_delta_table$patient_id)
check("clinical_45_crosswalk", length(prefixes) == 45L && all(prefixes %in% clinical_table$prefix))
primary_table <- clinical_table[clinical_table$time == "T1", , drop = FALSE]
recurrent_table <- clinical_table[clinical_table$time == "T2", , drop = FALSE]
primary_index <- match(prefixes, primary_table$prefix)
recurrent_index <- match(prefixes, recurrent_table$prefix)
matched <- all(!is.na(primary_index)) && all(!is.na(recurrent_index)) &&
    all(primary_table$`Patient ID`[primary_index] == recurrent_table$`Patient ID`[recurrent_index])
check("clinical_both_timepoints", matched)

status_table <- data.frame(
    patient_id          = target_delta_table$patient_id,
    clinical_patient_id = primary_table$`Patient ID`[primary_index],
    tumor_prefix        = prefixes,
    primary_steroid     = primary_table$steroid[primary_index],
    recurrent_steroid   = recurrent_table$steroid[recurrent_index],
    primary_raw         = primary_table$`Steroid treatment before the surgery`[primary_index],
    recurrent_raw       = recurrent_table$`Steroid treatment before the surgery`[recurrent_index],
    primary_score       = target_delta_table$primary_score,
    recurrent_score     = target_delta_table$recurrent_score,
    delta               = target_delta_table$delta
)
unknown <- is.na(status_table$primary_steroid) | is.na(status_table$recurrent_steroid)
transition <- paste(status_table$primary_steroid, status_table$recurrent_steroid, sep = "->")
status_table$group <- ifelse(unknown, "Missing", transition)
write.csv(
    status_table, file.path(output_dir, "GSE274546_steroid_patient_status.csv"),
    row.names = FALSE,
    na = ""
)

# Describe the recorded transitions. A singleton has no estimated interval.
group_tables <- list()
for (group in c("No->No", "No->Yes", "Yes->No", "Yes->Yes")) {
    values <- status_table$delta[status_table$group == group]
    interval <- if (length(values) > 1L) t.test(values)$conf.int else c(NA_real_, NA_real_)
    group_tables[[group]] <- data.frame(
        group      = group,
        n          = length(values),
        positive   = sum(values > 0),
        mean_delta = mean(values),
        ci_low     = interval[1],
        ci_high    = interval[2]
    )
}
group_table_total <- do.call(rbind, group_tables)
write.csv(
    group_table_total, file.path(output_dir, "GSE274546_steroid_group_summary.csv"),
    row.names = FALSE,
    na = ""
)
check("clinical_status_partition", sum(table(status_table$group)) == 45L)

# Use one published pair per strictly IDH-wildtype GeoMx patient and the stated trial labels.
geomx_path <- file.path(
    project, "write", "55_figure3_visual_polish", "Figure3", "source_data", "Figure3A_geomx_paired_dumbbell_source.csv"
)
geomx_table <- read.csv(geomx_path)
formal_geomx <- nrow(geomx_table) == 22L && anyDuplicated(geomx_table$Patient_No) == 0L &&
    all(geomx_table$entry == "strict_pass_idhwt") && all(geomx_table$trial_setting %in% c("Control", "Nivolumab"))
check("geomx_formal_22", formal_geomx)
difference_error <- max(abs(geomx_table$delta - (geomx_table$Recurrence - geomx_table$Primary)))
known_pairs <- unique(geomx_table$n_pairs)
check("geomx_delta_exact", difference_error < 1e-12 && length(known_pairs) == 1L && known_pairs == 22L)
write.csv(geomx_table, file.path(output_dir, "GeoMx_formal_patient_source_copy.csv"), row.names = FALSE, na = "")

# Calculate within-subgroup paired effects and the prespecified between-subgroup Welch contrast.
geomx_tables <- list()
for (group in c("Control", "Nivolumab")) {
    differences <- geomx_table$delta[geomx_table$trial_setting == group]
    paired_test <- t.test(differences)
    geomx_tables[[length(geomx_tables) + 1L]] <- data.frame(
        contrast = paste0(group, "_paired_delta"),
        n        = length(differences),
        positive = sum(differences > 0),
        estimate = mean(differences),
        ci_low   = paired_test$conf.int[1],
        ci_high  = paired_test$conf.int[2],
        p_value  = paired_test$p.value
    )
}
nivolumab <- geomx_table$delta[geomx_table$trial_setting == "Nivolumab"]
control <- geomx_table$delta[geomx_table$trial_setting == "Control"]
welch_test <- t.test(nivolumab, control, var.equal = FALSE)
geomx_tables[[length(geomx_tables) + 1L]] <- data.frame(
    contrast = "Nivolumab_minus_Control_interaction",
    n        = 22L,
    positive = NA_integer_,
    estimate = unname(diff(rev(welch_test$estimate))),
    ci_low   = welch_test$conf.int[1],
    ci_high  = welch_test$conf.int[2],
    p_value  = welch_test$p.value
)
interaction_table_total <- do.call(rbind, geomx_tables)
interaction_table_total$bh_three <- p.adjust(interaction_table_total$p_value, "BH")
write.csv(interaction_table_total, file.path(output_dir, "GeoMx_subgroup_interaction.csv"), row.names = FALSE, na = "")
check("geomx_group_counts", length(control) == 10L && length(nivolumab) == 12L)
check("geomx_published_overall", abs(mean(geomx_table$delta) - unique(geomx_table$mean_delta)) < 1e-12)

check_table_total <- do.call(rbind, check_rows)
write.csv(check_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat(nrow(check_table_total), " passing steroid and GeoMx checks\n", sep = "")
