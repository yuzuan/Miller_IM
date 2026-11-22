#!/usr/bin/env Rscript
# Purpose  Rebuild Figure 4 and Supplementary Figure 6 from the original paired proteome models.
# Inputs   PDC000514, MSV000087947 and Zenodo7646550 under write/31_proteomic_validation/source_metadata,
#          and the public Step 38 and 39 tables. MILLER_IM_DATA_ROOT can point to the input tree.
# Outputs  Source tables, manifests and notes under write/41_mg_inflammatory_sci_rebuild/{Figure4,SupplementaryFigure6}.
#          The corresponding figures directory contains the PDFs and PNGs under panel_library.
# Run      Rscript scripts/08a_paired_proteome_statistics.R with R, limma and fgsea.
# Preserve data.table grouping, joins and sorting to keep input order, formulas, thresholds and seeds stable.

suppressPackageStartupMessages({
    library(data.table)
    library(fgsea)
    library(BiocParallel)
    library(ggplot2)
    library(limma)
    library(patchwork)
    library(readxl)
    library(grid)
})

set.seed(20260713)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(script_path)
project_dir <- dirname(project_dir)
source(file.path(project_dir, "scripts", "00a_functions.R"))
data_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(project_dir, "data"))

step38_dir <- file.path(project_dir, "write", "38_independent_cohort_mg_inflammatory_recalculation")
step38_gsea_path <- file.path(step38_dir, "independent_fixed_program_targeted_gsea.csv")
step38_gene_path <- file.path(step38_dir, "independent_fixed_program_gene_direction.csv")
geomx_summary_path <- file.path(
    project_dir, "write", "39_independent_miller_mg_inflammatory_figures",
    "Figure3", "source_csv", "geomx_raw20_entry_level_summary.csv"
)
geomx_coverage_path <- file.path(
    project_dir, "write", "39_independent_miller_mg_inflammatory_figures",
    "Figure3", "source_csv", "geomx_raw20_coverage.csv"
)
required_local <- c(step38_gsea_path, step38_gene_path, geomx_summary_path, geomx_coverage_path)
missing_local <- required_local[!file.exists(required_local)]
if (length(missing_local)) stop("Missing local inputs: ", paste(missing_local, collapse = "; "))

pdc_source_dir <- locate_required(
    c(
        file.path(project_dir, "write", "31_proteomic_validation", "source_metadata", "PDC000514"),
        file.path(data_root, "write", "31_proteomic_validation", "source_metadata", "PDC000514")
    ),
    "PDC000514 source directory"
)

# Read the PDC matrix, remove summary rows and analyze the two peptide modes in their original order.
matrix_path <- file.path(pdc_source_dir, "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv")
source_table <- fread(matrix_path, check.names = FALSE)
source_table <- source_table[!Gene %in% c("Mean", "Median", "StdDev") & !is.na(Gene) & Gene != ""]

# Inputs are paired differences and the seed. Return the original two-sided sign-flip p value for their mean.
# The PDC, MSV and Zenodo analyses call this test with their original enumeration and sampling settings.
paired_sign_flip <- function(delta, seed) {
    delta <- delta[is.finite(delta)]
    n <- length(delta)
    if (n < 3L) {
        result <- NA_real_
        return(result)
    }
    observed <- abs(mean(delta))
    if (n <= 20L) {
        signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
        null <- abs(as.vector(signs %*% delta) / n)
        result <- mean(null >= observed - 1e-12)
        return(result)
    }
    set.seed(seed)
    total <- 100000L
    chunk <- 5000L
    extreme <- 0L
    for (start in seq.int(1L, total, by = chunk)) {
        current <- min(chunk, total - start + 1L)
        signs <- matrix(sample(c(-1, 1), current * n, replace = TRUE), nrow = current)
        null <- abs(rowMeans(signs * matrix(delta, nrow = current, ncol = n, byrow = TRUE)))
        extreme <- extreme + sum(null >= observed - 1e-12)
    }
    result <- (extreme + 1) / (total + 1)
    return(result)
}

# Inputs are ranking statistics and a gene set. Return the weighted running enrichment curve if coverage is adequate.
# Both PDC peptide modes use the same curve through analyse_matrix.
running_enrichment_curve <- function(stats, geneset) {
    stats <- stats[is.finite(stats)]
    stats <- sort(stats, decreasing = TRUE)
    geneset <- intersect(geneset, names(stats))
    hits <- names(stats) %in% geneset
    n_hits <- sum(hits)
    n_total <- length(stats)
    if (n_hits < 3L || n_hits == n_total) stop("Cannot build enrichment curve with current gene coverage.")
    hit_weight <- abs(stats[hits])
    hit_scale <- sum(hit_weight)
    phit <- ifelse(hits, abs(stats) / hit_scale, 0)
    pmiss <- ifelse(!hits, 1 / (n_total - n_hits), 0)
    result_table <- data.table(
        rank = seq_along(stats),
        gene = names(stats),
        rank_stat = unname(stats),
        hit  = hits,
        running_enrichment = cumsum(phit - pmiss)
    )
    return(result_table)
}

# Map measurement column names and a peptide prefix to sample IDs, patient IDs and conditions.
parse_samples <- function(columns, prefix) {
    selected <- grep(paste0("^", prefix, " KNCC_GBM[0-9]+_T[12]:"), columns, value = TRUE)
    sample_table <- data.table(
        column = selected,
        sample = sub(paste0("^", prefix, " (KNCC_GBM[0-9]+_T[12]):.*$"), "\\1", selected)
    )[, `:=`(
        patient   = sub("_T[12]$", "", sample),
        condition = fifelse(grepl("_T2$", sample), "Recurrent", "Primary")
    )]
    return(sample_table)
}

raw20 <- c(
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)

# Use the peptide prefix, analysis label and seed to return the PDC model, scores and enrichment curve.
analyse_matrix <- function(prefix, analysis, seed) {
    source_map_table <- parse_samples(colnames(source_table), prefix)
    sample_map_table_total <- unique(source_map_table[, .(sample, patient, condition)])
    sample_map_table_total <- merge(
        sample_map_table_total,
        source_map_table[, .(n_source_columns = .N), by = sample],
        by = "sample", sort = FALSE
    )
    pair_count_table <- dcast(sample_map_table_total, patient ~ condition, value.var = "sample", fun.aggregate = length)
    paired_patients <- pair_count_table[Primary == 1L & Recurrent == 1L, patient]
    sample_map_table_total <- sample_map_table_total[patient %in% paired_patients]
    setorder(sample_map_table_total, patient, condition)

    source_map_table[, analysis := analysis]
    source_map_table[, included_pair := patient %in% paired_patients]

    # Average duplicate columns for each sample, giving a matrix of proteins by patient-condition samples.
    source_columns <- source_map_table$column
    source_matrix <- as.matrix(source_table[, ..source_columns])
    storage.mode(source_matrix) <- "double"
    rownames(source_matrix) <- source_table$Gene
    colnames(source_matrix) <- source_map_table$sample
    average_matrix <- matrix(NA_real_, nrow(source_matrix), nrow(sample_map_table_total))
    for (sample_index in seq_len(nrow(sample_map_table_total))) {
        sample_id <- sample_map_table_total$sample[sample_index]
        average_matrix[, sample_index] <- rowMeans(
            source_matrix[, colnames(source_matrix) == sample_id, drop = FALSE], na.rm = TRUE
        )
    }
    average_matrix[!is.finite(average_matrix)] <- NA_real_
    rownames(average_matrix) <- source_table$Gene
    colnames(average_matrix) <- sample_map_table_total$sample

    present <- intersect(raw20, rownames(average_matrix))
    coverage_table <- data.table(
        analysis       = analysis,
        n_defined      = length(raw20),
        n_detected     = length(present),
        coverage       = length(present) / length(raw20),
        detected_genes = paste(present, collapse = ";"),
        missing_genes  = paste(setdiff(raw20, present), collapse = ";")
    )

    # Fit robust paired limma with patients as blocks and the last condition coefficient as the contrast.
    design_matrix <- model.matrix(
        ~ factor(sample_map_table_total$patient) +
            factor(sample_map_table_total$condition, levels = c("Primary", "Recurrent"))
    )
    coef_index <- ncol(design_matrix)
    fit_model <- eBayes(lmFit(average_matrix, design_matrix), robust = TRUE)
    gene_result_table <- as.data.table(
        topTable(fit_model, coef = coef_index, number = Inf, sort.by = "none"),
        keep.rownames = "gene"
    )
    gene_result_table[, analysis := analysis]
    gene_result_table[, measured_raw20 := gene %in% present]
    gene_result_table[, significant := adj.P.Val < 0.05]

    # Run fgsea on moderated t statistics with 3–500 genes, eps = 0 and SerialParam.
    ranks <- gene_result_table$t
    names(ranks) <- gene_result_table$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    fgsea_table <- as.data.table(
        fgseaMultilevel(
            pathways = list(Miller_raw20 = present),
            stats    = ranks,
            minSize  = 3,
            maxSize  = 500,
            eps      = 0,
            BPPARAM  = SerialParam(progressbar = FALSE)
        )
    )
    fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    fgsea_table[, analysis := analysis]
    gsea_curve_table <- running_enrichment_curve(ranks, present)
    gsea_curve_table[, analysis := analysis]

    # Score each pair from protein z scores measured at both time points. At least three shared proteins are required.
    z_matrix <- t(scale(t(average_matrix[present, , drop = FALSE])))
    z_matrix[!is.finite(z_matrix)] <- NA_real_
    paired_rows <- vector("list", length(unique(sample_map_table_total$patient)))
    patient_ids <- unique(sample_map_table_total$patient)
    for (patient_index in seq_along(patient_ids)) {
        patient_id <- patient_ids[patient_index]
        primary_sample <- sample_map_table_total[patient == patient_id & condition == "Primary", sample]
        recurrent_sample <- sample_map_table_total[patient == patient_id & condition == "Recurrent", sample]
        primary_values <- z_matrix[, primary_sample, drop = TRUE]
        recurrent_values <- z_matrix[, recurrent_sample, drop = TRUE]
        common <- is.finite(primary_values) & is.finite(recurrent_values)
        if (sum(common) < 3L) next
        paired_rows[[patient_index]] <- data.table(
            analysis       = analysis,
            patient        = patient_id,
            Primary        = mean(primary_values[common]),
            Recurrent      = mean(recurrent_values[common]),
            n_common_raw20 = sum(common),
            common_genes   = paste(present[common], collapse = ";")
        )
    }
    paired_score_table_total <- rbindlist(paired_rows)
    paired_score_table_total[, delta := Recurrent - Primary]
    summary_ci <- ci_mean(paired_score_table_total$delta)
    score_test_table <- data.table(
        analysis        = analysis,
        n_pairs         = nrow(paired_score_table_total),
        n_recurrence_up = sum(paired_score_table_total$delta > 0),
        mean_delta      = summary_ci[["mean"]],
        ci_low          = summary_ci[["low"]],
        ci_high         = summary_ci[["high"]],
        median_delta    = median(paired_score_table_total$delta),
        sign_flip_p     = paired_sign_flip(paired_score_table_total$delta, seed = seed),
        paired_t_p      = t.test(paired_score_table_total$delta)$p.value,
        wilcoxon_p      = wilcox.test(paired_score_table_total$delta, exact = FALSE)$p.value,
        n_defined       = coverage_table$n_defined,
        n_detected      = coverage_table$n_detected,
        coverage        = coverage_table$coverage
    )
    score_test_table[, gsea_NES := fgsea_table$NES]
    score_test_table[, gsea_p := fgsea_table$pval]
    score_test_table[, leading_edge := fgsea_table$leadingEdge]

    sample_flow_table <- data.table(
        analysis              = analysis,
        total_patients        = uniqueN(source_map_table$patient),
        primary_patients      = uniqueN(source_map_table[condition == "Primary", patient]),
        recurrent_patients    = uniqueN(source_map_table[condition == "Recurrent", patient]),
        paired_patients       = uniqueN(sample_map_table_total$patient),
        primary_only_patients = uniqueN(source_map_table[condition == "Primary" & !included_pair, patient]),
        duplicated_samples    = sum(source_map_table[, .N, by = sample]$N > 1L),
        extra_aliquot_columns = sum(pmax(source_map_table[, .N, by = sample]$N - 1L, 0L))
    )

    result <- list(
        source_map    = source_map_table,
        sample_map    = sample_map_table_total,
        coverage      = coverage_table,
        gene_results  = gene_result_table,
        fgsea         = fgsea_table,
        gsea_curve    = gsea_curve_table,
        paired_scores = paired_score_table_total,
        score_tests   = score_test_table,
        sample_flow   = sample_flow_table
    )
    return(result)
}

unique_res <- analyse_matrix("Unshared Log", "unique_peptides", 20260713L)
all_res <- analyse_matrix("Log", "all_peptides", 20260714L)
pdc <- list(unique = unique_res, all = all_res)

# Analyze the canonical MSV region first, then the average of replicate regions from the same input files.
msv_source_dir <- locate_required(
    c(
        file.path(project_dir, "write", "31_proteomic_validation", "source_metadata", "MSV000087947"),
        file.path(data_root, "write", "31_proteomic_validation", "source_metadata", "MSV000087947")
    ),
    "MSV000087947 source directory"
)
matrix_path <- file.path(msv_source_dir, "401_2022_2506_MOESM3_ESM.xlsx")
source_table <- as.data.table(read_excel(matrix_path, sheet = "Imputed values", skip = 3))
source_table[, protein_row := .I]
sample_columns <- grep("^GBM[0-9]+-(PM2?|RM2?)$", colnames(source_table), value = TRUE)
gene_rows <- vector("list", nrow(source_table))
for (index in seq_len(nrow(source_table))) {
    gene_names <- trimws(unlist(strsplit(as.character(source_table$Gene.names[index]), ";", fixed = TRUE)))
    gene_rows[[index]] <- data.table(protein_row = source_table$protein_row[index], gene = gene_names)
}
gene_map_table_total <- rbindlist(gene_rows)
gene_map_table_total <- gene_map_table_total[!is.na(gene) & gene != ""]
expanded_table_total <- merge(
    gene_map_table_total,
    source_table[,
        c("protein_row", sample_columns),
        with = FALSE],
    by = "protein_row"
)
gene_table <- expanded_table_total[, lapply(.SD, median, na.rm = TRUE), by = gene, .SDcols = sample_columns]
expression_matrix_all <- as.matrix(gene_table[, ..sample_columns])
storage.mode(expression_matrix_all) <- "double"
expression_matrix_all[!is.finite(expression_matrix_all)] <- NA_real_
rownames(expression_matrix_all) <- gene_table$gene

# Return the expression matrix and sample table for the canonical or averaged-replicate MSV regions.
prepare <- function(mode) {
    if (mode == "canonical_regions") {
        selected <- grep("^GBM[0-9]+-(PM|RM)$", colnames(expression_matrix_all), value = TRUE)
        expression_matrix <- expression_matrix_all[, selected, drop = FALSE]
        sample_map_table_total <- data.table(sample = selected)
        sample_map_table_total[, `:=`(
            patient   = sub("-(PM|RM)$", "", sample),
            condition = fifelse(grepl("-RM$", sample), "Recurrent", "Primary"),
            n_regions = 1L
        )]
    } else {
        source_sample_table <- data.table(source_sample = colnames(expression_matrix_all))
        source_sample_table[, `:=`(
            patient   = sub("-(PM2?|RM2?)$", "", source_sample),
            condition = fifelse(grepl("-RM2?$", source_sample), "Recurrent", "Primary")
        )]
        sample_map_table_total <- unique(source_sample_table[, .(patient, condition)])
        sample_map_table_total <- merge(
            sample_map_table_total, source_sample_table[, .(n_regions = .N), by = .(patient, condition)],
            by = c("patient", "condition"), sort = FALSE
        )
        setorder(sample_map_table_total, patient, condition)
        sample_map_table_total[, sample := paste(patient, condition, sep = "_")]
        expression_matrix <- matrix(NA_real_, nrow(expression_matrix_all), nrow(sample_map_table_total))
        for (i in seq_len(nrow(sample_map_table_total))) {
            selected_columns <- source_sample_table[
                patient == sample_map_table_total$patient[i] &
                    condition == sample_map_table_total$condition[i], source_sample
            ]
            expression_matrix[, i] <- rowMeans(expression_matrix_all[, selected_columns, drop = FALSE], na.rm = TRUE)
        }
        rownames(expression_matrix) <- rownames(expression_matrix_all)
        colnames(expression_matrix) <- sample_map_table_total$sample
    }
    setorder(sample_map_table_total, patient, condition)
    result <- list(expr = expression_matrix, sample_map = sample_map_table_total)
    return(result)
}

# Use the original mode and seed to return the paired MSV result tables for that region definition.
analyse <- function(mode, seed) {
    prepared <- prepare(mode)
    expression_matrix <- prepared$expr
    sample_map_table_total <- prepared$sample_map
    present <- intersect(raw20, rownames(expression_matrix))
    design_matrix <- model.matrix(
        ~ factor(sample_map_table_total$patient) +
            factor(sample_map_table_total$condition, levels = c("Primary", "Recurrent"))
    )
    coef_index <- ncol(design_matrix)
    fit_model <- eBayes(lmFit(expression_matrix, design_matrix), robust = TRUE)
    gene_result_table <- as.data.table(
        topTable(fit_model, coef = coef_index, number = Inf, sort.by = "none"),
        keep.rownames = "gene"
    )

    # Rank by moderated t from the paired MSV model and test the measurable raw20 proteins by GSEA.
    ranks <- gene_result_table$t
    names(ranks) <- gene_result_table$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    fgsea_table <- as.data.table(
        fgseaMultilevel(
            pathways = list(Miller_raw20 = present),
            stats    = ranks,
            minSize  = 3,
            maxSize  = 500,
            eps      = 0,
            BPPARAM  = SerialParam(progressbar = FALSE)
        )
    )
    fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]

    # Standardize each protein across MSV samples. Average the raw20 proteins measured at both times for each pair.
    z_matrix <- t(scale(t(expression_matrix[present, , drop = FALSE])))
    z_matrix[!is.finite(z_matrix)] <- NA_real_
    paired_rows <- vector("list", length(unique(sample_map_table_total$patient)))
    patient_ids <- unique(sample_map_table_total$patient)
    for (patient_index in seq_along(patient_ids)) {
        patient_id <- patient_ids[patient_index]
        primary_sample <- sample_map_table_total[patient == patient_id & condition == "Primary", sample]
        recurrent_sample <- sample_map_table_total[patient == patient_id & condition == "Recurrent", sample]
        primary_values <- z_matrix[, primary_sample, drop = TRUE]
        recurrent_values <- z_matrix[, recurrent_sample, drop = TRUE]
        common <- is.finite(primary_values) & is.finite(recurrent_values)
        if (sum(common) < 3L) next
        paired_rows[[patient_index]] <- data.table(
            analysis  = mode,
            patient   = patient_id,
            Primary   = mean(primary_values[common]),
            Recurrent = mean(recurrent_values[common]),
            n_common_raw20 = sum(common)
        )
    }
    paired_score_table_total <- rbindlist(paired_rows)
    paired_score_table_total[, delta := Recurrent - Primary]
    ci <- ci_mean(paired_score_table_total$delta)
    result_table <- data.table(
        analysis        = mode,
        n_pairs         = nrow(paired_score_table_total),
        n_recurrence_up = sum(paired_score_table_total$delta > 0),
        mean_delta      = ci[["mean"]],
        ci_low          = ci[["low"]],
        ci_high         = ci[["high"]],
        median_delta    = median(paired_score_table_total$delta),
        sign_flip_p     = paired_sign_flip(paired_score_table_total$delta, seed = seed),
        gsea_NES        = fgsea_table$NES,
        gsea_p          = fgsea_table$pval,
        n_detected      = length(present),
        coverage        = length(present) / length(raw20),
        detected_genes  = paste(present, collapse = ";"),
        paired_scores   = list(paired_score_table_total)
    )
    return(result_table)
}

canonical_table <- analyse("canonical_regions", 20260715L)
averaged_table <- analyse("averaged_repeated_regions", 20260716L)
msv_table_total <- rbindlist(list(canonical_table, averaged_table), fill = TRUE)

# Use the published processed Zenodo matrix and limma results with the original handling of duplicate genes.
zenodo_source_root <- locate_required(
    c(
        file.path(
            project_dir, "write", "31_proteomic_validation", "source_metadata", "Zenodo7646550",
            "extracted", "MiguelCos-gbm_manuscript_data_analysis-3c4349d"
        ),
        file.path(
            data_root, "write", "31_proteomic_validation", "source_metadata", "Zenodo7646550",
            "extracted", "MiguelCos-gbm_manuscript_data_analysis-3c4349d"
        )
    ),
    "Zenodo7646550 source directory"
)

matrix_path <- file.path(
    zenodo_source_root, "data", "specific_search_fragpipe17",
    "specific_no_ptms_2", "tmt-report", "abundance_protein_MD.tsv"
)
protein_table <- fread(matrix_path, check.names = FALSE)
sample_columns <- grep("^[0-9]+-(prim|rec)$", colnames(protein_table), value = TRUE)
protein_table <- protein_table[!is.na(Gene) & Gene != ""]
gene_table <- protein_table[, lapply(.SD, median, na.rm = TRUE), by = Gene, .SDcols = sample_columns]
expression_matrix <- as.matrix(gene_table[, ..sample_columns])
storage.mode(expression_matrix) <- "double"
expression_matrix[!is.finite(expression_matrix)] <- NA_real_
rownames(expression_matrix) <- gene_table$Gene
present <- intersect(raw20, rownames(expression_matrix))

# Standardize within the Zenodo matrix and average raw20 genes measurable at both time points for each patient.
z_matrix <- t(scale(t(expression_matrix[present, , drop = FALSE])))
z_matrix[!is.finite(z_matrix)] <- NA_real_
paired_rows <- vector("list", length(sample_columns) / 2L)
for (i in seq_len(length(sample_columns) / 2L)) {
    patient_id <- as.character(i)
    primary_sample <- paste0(patient_id, "-prim")
    recurrent_sample <- paste0(patient_id, "-rec")
    primary_values <- z_matrix[, primary_sample, drop = TRUE]
    recurrent_values <- z_matrix[, recurrent_sample, drop = TRUE]
    common <- is.finite(primary_values) & is.finite(recurrent_values)
    if (sum(common) < 3L) next
    paired_rows[[i]] <- data.table(
        patient   = patient_id,
        Primary   = mean(primary_values[common]),
        Recurrent = mean(recurrent_values[common]),
        n_common_raw20 = sum(common)
    )
}
paired_score_table_total <- rbindlist(paired_rows)
paired_score_table_total[, delta := Recurrent - Primary]
ci <- ci_mean(paired_score_table_total$delta)

# Within each gene, sort by absolute t and retain the first row. Run GSEA on the fixed gene set.
limma_path <- file.path(zenodo_source_root, "results", "limma_proteomics_results.tsv")
limma_result_table <- fread(limma_path)
limma_result_table <- limma_result_table[!is.na(Gene) & Gene != ""]
limma_result_table[, abs_t := abs(t)]
setorder(limma_result_table, Gene, -abs_t)
limma_result_table <- limma_result_table[!duplicated(Gene)]
ranks <- limma_result_table$t
names(ranks) <- limma_result_table$Gene
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
fgsea_table <- as.data.table(
    fgseaMultilevel(
        pathways = list(Miller_raw20 = present),
        stats    = ranks,
        minSize  = 3,
        maxSize  = 500,
        eps      = 0,
        BPPARAM  = SerialParam(progressbar = FALSE)
    )
)
fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]

zenodo_table <- data.table(
    analysis        = "single_matrix",
    n_pairs         = nrow(paired_score_table_total),
    n_recurrence_up = sum(paired_score_table_total$delta > 0),
    mean_delta      = ci[["mean"]],
    ci_low          = ci[["low"]],
    ci_high         = ci[["high"]],
    median_delta    = median(paired_score_table_total$delta),
    sign_flip_p     = paired_sign_flip(paired_score_table_total$delta, seed = 20260717L),
    gsea_NES        = fgsea_table$NES,
    gsea_p          = fgsea_table$pval,
    n_detected      = length(present),
    coverage        = length(present) / length(raw20),
    detected_genes  = paste(present, collapse = ";"),
    paired_scores   = list(paired_score_table_total)
)

pdc_unique_gene_table <- copy(pdc$unique$gene_results)[gene %in% raw20]
pdc_unique_gene_table[, `:=`(
    ci_low   = logFC - 1.96 * abs(logFC / t),
    ci_high  = logFC + 1.96 * abs(logFC / t),
    analysis = "unique_peptides"
)]



figure4_delta_table <- copy(pdc$unique$paired_scores)
figure4_delta_table[, patient := factor(patient, levels = patient[order(delta)])]
delta_test_table <- pdc$unique$score_tests

# Build the left half violin with nrd0 bandwidth, 256 density points, center 1 and width 0.25.
violin_values <- figure4_delta_table$delta
violin_values <- violin_values[is.finite(violin_values)]
violin_density <- density(violin_values, n = 256, bw = "nrd0")
violin_scaled <- if (max(violin_density$y) > 0) {
    violin_density$y / max(violin_density$y) * 0.25
} else {
    rep(0, length(violin_density$y))
}
violin_side <- 1 - violin_scaled
delta_violin_table <- data.table(
    x = c(rep(1, length(violin_density$x)), rev(violin_side)),
    y = c(violin_density$x, rev(violin_density$x))
)

# Keep the original point jitter seed, locations and order of random calls.
set.seed(20260713)
figure4_delta_table[, jitter_x := 1.08 + runif(.N, -0.06, 0.06)]
figure4_source_dir <- file.path(
    project_dir, "write", "41_mg_inflammatory_sci_rebuild", "Figure4", "source_data"
)
dir.create(figure4_source_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(figure4_delta_table, file.path(figure4_source_dir, "Figure4A_patient_delta_source.csv"))
npg <- c(
    red       = "#E64B35",
    blue      = "#4DBBD5",
    teal      = "#00A087",
    navy      = "#3C5488",
    orange    = "#F39B7F",
    purple    = "#8491B4",
    green     = "#91D1C2",
    darkred   = "#DC0000",
    grey      = "#7A7A7A",
    lightgrey = "#D9D9D9"
)

theme_main <- theme_classic(base_size = 10)
theme_main <- theme_main + theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    legend.title     = element_text(color = "black"),
    legend.text      = element_text(color = "black"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text       = element_text(face = "bold"),
    plot.margin      = margin(6, 6, 6, 6)
)


# In panel A, use the calculated paired differences, mean interval and sign-flip p value.
p4a <- ggplot()
p4a <- p4a + geom_polygon(
    data = delta_violin_table,
    aes(x = x, y = y),
    fill  = "#C8D8DC",
    alpha = 0.72,
    color = "#76A5AF",
    linewidth = 0.4
)
p4a <- p4a + geom_boxplot(
    data  = figure4_delta_table, aes(x = 1, y = delta, group = 1),
    width = 0.10, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35
)
p4a <- p4a + geom_point(
    data = figure4_delta_table, aes(x = jitter_x, y = delta),
    size = 1.15, shape = 21, stroke = 0.2, fill = npg["red"], color = "white", alpha = 0.78
)
p4a <- p4a + geom_segment(
    aes(x = 0.80, xend = 0.80, y = delta_test_table$ci_low, yend = delta_test_table$ci_high),
    linewidth = 1.0, color = npg["red"]
)
p4a <- p4a + geom_point(
    aes(x = 0.80, y = delta_test_table$mean_delta),
    size   = 2.9,
    shape  = 21,
    fill   = npg["red"],
    color  = "black",
    stroke = 0.25
)
p4a <- p4a + geom_hline(yintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4)
p4a <- p4a + annotate(
    "text",
    x     = 1.23,
    y     = max(figure4_delta_table$delta) + 0.06,
    hjust = 1,
    vjust = 1,
    size  = 2.95,
    label = sprintf(
        "mean %.3f\n95%% CI %.3f to %.3f\nn = %d\n%d/%d up\nP %.2g",
        delta_test_table$mean_delta, delta_test_table$ci_low, delta_test_table$ci_high,
        delta_test_table$n_pairs, delta_test_table$n_recurrence_up,
        delta_test_table$n_pairs, delta_test_table$sign_flip_p
    )
)
p4a <- p4a + coord_cartesian(
    xlim = c(0.56, 1.24),
    ylim = c(min(figure4_delta_table$delta) - 0.12, max(figure4_delta_table$delta) + 0.12),
    clip = "off"
)
p4a <- p4a + labs(
    x        = NULL,
    y        = "Recurrent - primary raw20 score",
    title    = "PDC unique peptides",
    subtitle = "Patient-level raw20 delta distribution"
)
p4a <- p4a + theme_main
p4a <- p4a + theme(
    plot.title      = element_text(size = 8.4, face = "bold"),
    plot.subtitle   = element_text(size = 7.0, face = "plain"),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "none"
)
figure4_plot_dir <- file.path(
    project_dir, "figures", "41_mg_inflammatory_sci_rebuild", "Figure4", "panel_library"
)
save_plot_pair(p4a, "fig4a_pdc_patient_delta", figure4_plot_dir, width = 3.25, height = 4.1)
figure4_panel_1_table <- panel_manifest(
    "Figure4A",
    "fig4a_pdc_patient_delta",
    "PDC unique-peptide patient delta distribution with mean and 95% CI.",
    figure4_plot_dir
)

# In panel B, show measurable raw20 proteins with the original moderated log fold changes and 95% intervals.
figure4_forest_table <- copy(pdc_unique_gene_table[measured_raw20 == TRUE])
setorder(figure4_forest_table, logFC)
figure4_forest_table[, gene := factor(gene, levels = gene)]
figure4_forest_table[, significant_label := factor(as.character(significant), levels = c("FALSE", "TRUE"))]
fwrite(figure4_forest_table, file.path(figure4_source_dir, "Figure4B_measurable_protein_forest_source.csv"))
p4b <- ggplot(figure4_forest_table, aes(x = logFC, y = gene))
p4b <- p4b + geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4)
p4b <- p4b + geom_segment(aes(x = ci_low, xend = ci_high, yend = gene), linewidth = 0.85, color = npg["grey"])
p4b <- p4b + geom_point(aes(fill = significant_label), shape = 21, size = 2.8, color = "black", stroke = 0.25)
p4b <- p4b + scale_fill_manual(
    values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"),
    labels = c(`TRUE` = "Yes", `FALSE` = "No")
)
p4b <- p4b + labs(
    x        = "Proteome logFC (recurrent vs primary)",
    y        = NULL,
    fill     = "FDR < 0.05",
    title    = "Measured raw20 proteins in PDC unique peptides",
    subtitle = "Moderated logFC with 95% CI"
)
p4b <- p4b + theme_main
p4b <- p4b + theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 8.2, face = "bold"),
    plot.subtitle   = element_text(size = 7.0, face = "plain")
)
save_plot_pair(p4b, "fig4b_measurable_protein_forest", figure4_plot_dir, width = 4.9, height = 4.3)
figure4_panel_2_table <- panel_manifest(
    "Figure4B",
    "fig4b_measurable_protein_forest",
    "Measured raw20 proteins in PDC unique-peptide analysis with 95% CI.",
    figure4_plot_dir
)

# In panel C, retain the ranked curve, hit positions and moderated t track in their original height ratios.
figure4_gsea_curve_table <- copy(pdc$unique$gsea_curve)
leading_edge_genes <- unlist(strsplit(pdc$unique$score_tests$leading_edge, ";", fixed = TRUE))
figure4_gsea_curve_table[, leading_edge := gene %in% leading_edge_genes]
fwrite(figure4_gsea_curve_table, file.path(figure4_source_dir, "Figure4C_protein_ranked_gsea_curve_source.csv"))
p4c_top <- ggplot(figure4_gsea_curve_table, aes(x = rank, y = running_enrichment))
p4c_top <- p4c_top + geom_line(color = npg["navy"], linewidth = 1.05)
p4c_top <- p4c_top + geom_area(fill = npg["navy"], alpha = 0.12)
p4c_top <- p4c_top + geom_hline(yintercept = 0, color = npg["grey"], linewidth = 0.4)
p4c_top <- p4c_top + labs(y = "Running enrichment", x = NULL)
p4c_top <- p4c_top + theme_main
p4c_top <- p4c_top + theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
)

# Mark raw20 proteins at their positions in the full proteome ranking.
p4c_mid <- ggplot(figure4_gsea_curve_table[hit == TRUE], aes(x = rank, xend = rank, y = 0, yend = 1))
p4c_mid <- p4c_mid + geom_segment(color = "black", linewidth = 0.28)
p4c_mid <- p4c_mid + coord_cartesian(ylim = c(0, 1), expand = FALSE)
p4c_mid <- p4c_mid + theme_void()

# Plot moderated t on the same ranking axis with the original colors for positive and negative values.
p4c_bottom <- ggplot(figure4_gsea_curve_table, aes(x = rank, y = rank_stat))
p4c_bottom <- p4c_bottom + geom_area(data = figure4_gsea_curve_table[rank_stat >= 0], fill = npg["red"], alpha = 0.75)
p4c_bottom <- p4c_bottom + geom_area(data = figure4_gsea_curve_table[rank_stat < 0], fill = npg["blue"], alpha = 0.75)
p4c_bottom <- p4c_bottom + geom_hline(yintercept = 0, color = npg["grey"], linewidth = 0.35)
p4c_bottom <- p4c_bottom + labs(x = "Ranked proteins", y = "Moderated t")
p4c_bottom <- p4c_bottom + theme_main
p4c_bottom <- p4c_bottom + theme(axis.title.y = element_text(size = 7))
p4c <- (p4c_top / p4c_mid / p4c_bottom)
p4c <- p4c + plot_layout(heights = c(3.2, 0.45, 1.0))
p4c <- p4c + plot_annotation(
    title    = "PDC unique peptides",
    subtitle = sprintf(
        "NES %.2f, nominal P %.2g, %d/%d proteins measured",
        pdc$unique$score_tests$gsea_NES, pdc$unique$score_tests$gsea_p,
        pdc$unique$score_tests$n_detected, length(raw20)
    )
)
p4c <- p4c & theme(
    plot.title    = element_text(size = 8.4, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.0, face = "plain", hjust = 0)
)
save_plot_pair(p4c, "fig4c_protein_ranked_gsea", figure4_plot_dir, width = 4.9, height = 4.2)
figure4_panel_3_table <- panel_manifest(
    "Figure4C",
    "fig4c_protein_ranked_gsea",
    "Complete ranked-list GSEA curve for PDC unique-peptide proteome.",
    figure4_plot_dir
)

# In panel D, align GeoMx and proteome estimates for shared genes and keep each study's effect scale.
geomx_gene_path <- file.path(
    project_dir, "write", "39_independent_miller_mg_inflammatory_figures", "Figure3", "source_csv",
    "geomx_raw20_gene_level_summary.csv"
)
geomx_gene_table <- fread(geomx_gene_path)

# Check the frozen title counts of 22 and 105 pairs against the current GeoMx and PDC inputs.
geomx_pair_counts <- unique(geomx_gene_table[entry == "strict_pass_idhwt", n_pairs])
stopifnot(identical(as.integer(geomx_pair_counts), 22L))
stopifnot(identical(as.integer(pdc$unique$score_tests$n_pairs), 105L))

geomx_dual_table <- geomx_gene_table[
    entry == "strict_pass_idhwt",
    .(gene, geomx_mean = mean_log2_delta, geomx_ci_low = ci_low, geomx_ci_high = ci_high, geomx_fdr = fdr)
]
pdc_dual_table <- pdc_unique_gene_table[
    measured_raw20 == TRUE,
    .(gene, proteome_mean = logFC, proteome_ci_low = ci_low, proteome_ci_high = ci_high, proteome_fdr = adj.P.Val)
]
dual_modal_table_total <- merge(geomx_dual_table, pdc_dual_table, by = "gene")
setorder(dual_modal_table_total, proteome_mean)
dual_modal_table_total[, gene := factor(gene, levels = gene)]
dual_modal_table_total[, `:=`(
    geomx_sig_label    = factor(as.character(geomx_fdr < 0.05), levels = c("FALSE", "TRUE")),
    proteome_sig_label = factor(as.character(proteome_fdr < 0.05), levels = c("FALSE", "TRUE"))
)]
fwrite(dual_modal_table_total, file.path(figure4_source_dir, "Figure4D_cross_modal_dual_scale_source.csv"))
p4d_left <- ggplot(dual_modal_table_total, aes(x = geomx_mean, y = gene))
p4d_left <- p4d_left + geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.35)
p4d_left <- p4d_left + geom_segment(
    aes(x = geomx_ci_low, xend = geomx_ci_high, yend = gene),
    linewidth = 0.8,
    color = npg["grey"]
)
p4d_left <- p4d_left + geom_point(aes(fill = geomx_sig_label), shape = 21, size = 2.7, color = "black", stroke = 0.25)
p4d_left <- p4d_left + scale_fill_manual(values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"), guide = "none")
p4d_left <- p4d_left + labs(x = "GeoMx mean log2 delta", y = NULL, title = "IBA1+ GeoMx (22 pairs)")
p4d_left <- p4d_left + theme_main
p4d_left <- p4d_left + theme(plot.title = element_text(size = 8.0, face = "bold", color = npg["navy"]))
p4d_right <- ggplot(dual_modal_table_total, aes(x = proteome_mean, y = gene))
p4d_right <- p4d_right + geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.35)
p4d_right <- p4d_right + geom_segment(
    aes(x = proteome_ci_low, xend = proteome_ci_high, yend = gene),
    linewidth = 0.8,
    color = npg["grey"]
)
p4d_right <- p4d_right + geom_point(
    aes(fill = proteome_sig_label),
    shape  = 21,
    size   = 2.7,
    color  = "black",
    stroke = 0.25
)
p4d_right <- p4d_right + scale_fill_manual(values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"), guide = "none")
p4d_right <- p4d_right + labs(x = "Proteome logFC", y = NULL, title = "PDC unique peptides (105 pairs)")
p4d_right <- p4d_right + theme_main
p4d_right <- p4d_right + theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title   = element_text(size = 8.0, face = "bold", color = npg["darkred"])
)
p4d <- p4d_left
p4d <- p4d + p4d_right
p4d <- p4d + plot_layout(widths = c(1, 1))
p4d <- p4d + plot_annotation(
    title    = "Same measurable genes across GeoMx and proteome",
    subtitle = "Each column keeps its own effect scale"
)
p4d <- p4d & theme(
    plot.title    = element_text(size = 8.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.0, face = "plain", hjust = 0)
)
save_plot_pair(p4d, "fig4d_cross_modal_dual_scale", figure4_plot_dir, width = 7.3, height = 4.3)
figure4_panel_4_table <- panel_manifest(
    "Figure4D",
    "fig4d_cross_modal_dual_scale",
    "Same-gene GeoMx and proteome effects shown in separate columns with separate scales.",
    figure4_plot_dir
)

# Write the source notes, panel manifest and analysis summary.
figure4_summary_table <- data.table(
    input = c(
        "Step38 paired scRNA pseudobulk GSEA",
        "PDC000514 raw TMT matrix",
        "Artzi IBA1+ GeoMx raw20 gene summary",
        "MSV000087947 paired proteome sensitivity",
        "Zenodo7646550 paired proteome sensitivity"
    ),
    location = c(
        step38_dir,
        pdc_source_dir,
        file.path(
            project_dir,
            "write",
            "39_independent_miller_mg_inflammatory_figures",
            "Figure3",
            "source_csv",
            "geomx_raw20_gene_level_summary.csv"
        ),
        msv_source_dir,
        zenodo_source_root
    )
)
figure4_write_dir <- dirname(figure4_source_dir)
writeLines(c(
    "# Figure4 summary",
    "",
    sprintf(
        "- Main protein result uses only Miller raw20 and PDC unique peptides: %d/%d proteins measured.",
        pdc$unique$score_tests$n_detected,
        length(raw20)
    ),
    sprintf(
        "- Patient delta: mean %.3f, 95%% CI %.3f to %.3f, %d/%d recurrence-up, P %.2g.",
        delta_test_table$mean_delta,
        delta_test_table$ci_low,
        delta_test_table$ci_high,
        delta_test_table$n_recurrence_up,
        delta_test_table$n_pairs,
        delta_test_table$sign_flip_p
    ),
    sprintf(
        "- Protein GSEA: NES %.3f, nominal P %.3g.",
        pdc$unique$score_tests$gsea_NES,
        pdc$unique$score_tests$gsea_p
    ),
    sprintf(
        "- Cross-modal dual-column panel is limited to %d genes measurable in both GeoMx and PDC.",
        nrow(dual_modal_table_total)
    )
), file.path(figure4_write_dir, "Figure4_summary.md"), useBytes = TRUE)
figure4_manifest_table_total <- rbindlist(list(
    figure4_panel_1_table, figure4_panel_2_table, figure4_panel_3_table, figure4_panel_4_table
))
fwrite(figure4_manifest_table_total, file.path(figure4_write_dir, "Figure4_panel_manifest.csv"))
fwrite(figure4_summary_table, file.path(figure4_write_dir, "Figure4_inputs.csv"))

# In supplementary panel A, count inclusions, exclusions and merged repeated samples from the current sample map.
supp6_flow_table_total <- rbindlist(list(
    data.table(step = "Patients in matrix", value = pdc$unique$sample_flow$total_patients),
    data.table(step = "Primary available", value = pdc$unique$sample_flow$primary_patients),
    data.table(step = "Recurrent available", value = pdc$unique$sample_flow$recurrent_patients),
    data.table(step = "Paired patients", value = pdc$unique$sample_flow$paired_patients),
    data.table(step = "Primary-only excluded", value = pdc$unique$sample_flow$primary_only_patients),
    data.table(step = "Samples merged", value = pdc$unique$sample_flow$duplicated_samples),
    data.table(step = "Extra aliquot columns", value = pdc$unique$sample_flow$extra_aliquot_columns)
))
supp6_source_dir <- file.path(
    project_dir, "write", "41_mg_inflammatory_sci_rebuild", "SupplementaryFigure6", "source_data"
)
dir.create(supp6_source_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(supp6_flow_table_total, file.path(supp6_source_dir, "SupplementaryFigure6A_sample_flow_source.csv"))
supp6_box_table <- data.table(
    label = c(
        sprintf("All\n%d", pdc$unique$sample_flow$total_patients),
        sprintf("Primary\n%d", pdc$unique$sample_flow$primary_patients),
        sprintf("Recurrent\n%d", pdc$unique$sample_flow$recurrent_patients),
        sprintf("Paired\n%d", pdc$unique$sample_flow$paired_patients),
        sprintf("Excluded\n%d", pdc$unique$sample_flow$primary_only_patients),
        sprintf("Merged\n%d", pdc$unique$sample_flow$duplicated_samples)
    ),
    xmin = c(0.0, 1.8, 3.6, 5.4, 3.6, 5.4),
    xmax = c(1.2, 3.0, 4.8, 6.6, 4.8, 6.6),
    ymin = c(0.15, 0.15, 0.15, 0.15, -1.15, -1.15),
    ymax = c(1.15, 1.15, 1.15, 1.15, -0.15, -0.15),
    fill = c(npg["blue"], npg["teal"], npg["teal"], npg["red"], npg["grey"], npg["orange"])
)
p6a <- ggplot()
p6a <- p6a + geom_rect(
    data = supp6_box_table,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = supp6_box_table$fill, color = "white", linewidth = 0.6
)
p6a <- p6a + geom_text(
    data = supp6_box_table,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    color = "white", size = 3.2, fontface = "bold", lineheight = 0.95
)
p6a <- p6a + annotate(
    "segment",
    x     = 1.2,
    xend  = 1.8,
    y     = 0.65,
    yend  = 0.65,
    arrow = arrow(length = unit(0.10, "inches"))
)
p6a <- p6a + annotate(
    "segment",
    x     = 3.0,
    xend  = 3.6,
    y     = 0.65,
    yend  = 0.65,
    arrow = arrow(length = unit(0.10, "inches"))
)
p6a <- p6a + annotate(
    "segment",
    x     = 4.8,
    xend  = 5.4,
    y     = 0.65,
    yend  = 0.65,
    arrow = arrow(length = unit(0.10, "inches"))
)
p6a <- p6a + annotate(
    "segment",
    x     = 4.2,
    xend  = 4.2,
    y     = 0.15,
    yend  = -0.15,
    arrow = arrow(length = unit(0.10, "inches"))
)
p6a <- p6a + annotate(
    "segment",
    x     = 6.0,
    xend  = 6.0,
    y     = 0.15,
    yend  = -0.15,
    arrow = arrow(length = unit(0.10, "inches"))
)
p6a <- p6a + coord_cartesian(xlim = c(-0.3, 6.75), ylim = c(-1.28, 1.28), clip = "off")
p6a <- p6a + theme_void()
supp6_plot_dir <- file.path(
    project_dir, "figures", "41_mg_inflammatory_sci_rebuild", "SupplementaryFigure6", "panel_library"
)
save_plot_pair(p6a, "supp6a_sample_flow", supp6_plot_dir, width = 6.8, height = 2.5)
supp6_panel_1_table <- panel_manifest(
    "SupplementaryFigure6A",
    "supp6a_sample_flow",
    "PDC000514 sample flow and exclusions.",
    supp6_plot_dir
)

# In supplementary panel B, report raw20 coverage by cohort and peptide mode against the fixed 20 genes.
supp6_coverage_table_total <- rbindlist(list(
    data.table(
        dataset  = "PDC unique",
        pairs    = pdc$unique$score_tests$n_pairs,
        detected = pdc$unique$coverage$n_detected,
        total    = 20L,
        coverage = pdc$unique$coverage$coverage
    ),
    data.table(
        dataset  = "PDC all",
        pairs    = pdc$all$score_tests$n_pairs,
        detected = pdc$all$coverage$n_detected,
        total    = 20L,
        coverage = pdc$all$coverage$coverage
    ),
    data.table(
        dataset = "MSV canonical",
        pairs   = msv_table_total[analysis == "canonical_regions",
            n_pairs],
        detected = msv_table_total[analysis == "canonical_regions",
            n_detected],
        total    = 20L,
        coverage = msv_table_total[analysis == "canonical_regions",
            coverage]
    ),
    data.table(
        dataset = "MSV averaged",
        pairs   = msv_table_total[analysis == "averaged_repeated_regions",
            n_pairs],
        detected = msv_table_total[analysis == "averaged_repeated_regions",
            n_detected],
        total    = 20L,
        coverage = msv_table_total[analysis == "averaged_repeated_regions",
            coverage]
    ),
    data.table(
        dataset  = "Zenodo",
        pairs    = zenodo_table$n_pairs,
        detected = zenodo_table$n_detected,
        total    = 20L,
        coverage = zenodo_table$coverage
    )
))
supp6_coverage_table_total[, dataset := factor(dataset, levels = rev(dataset))]
fwrite(supp6_coverage_table_total, file.path(supp6_source_dir, "SupplementaryFigure6B_coverage_source.csv"))
p6b <- ggplot(supp6_coverage_table_total, aes(x = coverage, y = dataset, color = dataset))
p6b <- p6b + geom_segment(aes(x = 0, xend = coverage, yend = dataset), linewidth = 1.0, color = npg["lightgrey"])
p6b <- p6b + geom_point(size = 3.2)
p6b <- p6b + geom_text(
    aes(label = sprintf("%d/20  |  n=%d", detected, pairs)),
    nudge_x = 0.04,
    hjust   = 0,
    size    = 2.75,
    color   = "black"
)
p6b <- p6b + scale_color_manual(values = c(
    "PDC unique" = unname(npg["red"]),
    "PDC all" = unname(npg["teal"]),
    "MSV canonical" = unname(npg["blue"]),
    "MSV averaged" = unname(npg["orange"]),
    "Zenodo" = unname(npg["navy"])
))
p6b <- p6b + coord_cartesian(xlim = c(0, 0.92), clip = "off")
p6b <- p6b + labs(x = "Raw20 protein coverage", y = NULL, color = NULL)
p6b <- p6b + theme_main
p6b <- p6b + theme(legend.position = "none")
save_plot_pair(p6b, "supp6b_dataset_coverage", supp6_plot_dir, width = 5.2, height = 3.0)
supp6_panel_2_table <- panel_manifest(
    "SupplementaryFigure6B",
    "supp6b_dataset_coverage",
    "Raw20 protein coverage across proteome datasets.",
    supp6_plot_dir
)

# In supplementary panel C, show proteins measured in each pair with the original jitter and random draw order.
supp6_count_table_total <- rbindlist(list(
    copy(pdc$unique$paired_scores)[, analysis := "Unique peptides"],
    copy(pdc$all$paired_scores)[, analysis := "All peptides"]
))
fwrite(
    supp6_count_table_total,
    file.path(supp6_source_dir, "SupplementaryFigure6C_measurable_protein_counts_source.csv")
)
p6c <- ggplot(supp6_count_table_total, aes(x = analysis, y = n_common_raw20, fill = analysis, color = analysis))
p6c <- p6c + geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35)
p6c <- p6c + geom_jitter(
    width  = 0.08,
    height = 0.05,
    size   = 0.95,
    alpha  = 0.72,
    shape  = 21,
    stroke = 0.15,
    fill   = "white"
)
p6c <- p6c + scale_fill_manual(values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"])))
p6c <- p6c + scale_color_manual(
    values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"]))
)
p6c <- p6c + labs(x = NULL, y = "Detectable raw20 proteins per pair")
p6c <- p6c + theme_main
p6c <- p6c + theme(legend.position = "none")
save_plot_pair(p6c, "supp6c_measurable_protein_counts", supp6_plot_dir, width = 3.9, height = 4.0)
supp6_panel_3_table <- panel_manifest(
    "SupplementaryFigure6C",
    "supp6c_measurable_protein_counts",
    "Detectable raw20 protein counts per paired patient in PDC.",
    supp6_plot_dir
)

# In supplementary panel D, report the score changes and confidence intervals on the score scale.
supp6_score_table_total <- rbindlist(list(
    data.table(
        analysis   = "Unique peptides",
        mean_delta = pdc$unique$score_tests$mean_delta,
        ci_low     = pdc$unique$score_tests$ci_low,
        ci_high    = pdc$unique$score_tests$ci_high,
        fdr        = pdc$unique$score_tests$sign_flip_p,
        n_pairs    = pdc$unique$score_tests$n_pairs,
        n_recurrence_up = pdc$unique$score_tests$n_recurrence_up
    ),
    data.table(
        analysis   = "All peptides",
        mean_delta = pdc$all$score_tests$mean_delta,
        ci_low     = pdc$all$score_tests$ci_low,
        ci_high    = pdc$all$score_tests$ci_high,
        fdr        = pdc$all$score_tests$sign_flip_p,
        n_pairs    = pdc$all$score_tests$n_pairs,
        n_recurrence_up = pdc$all$score_tests$n_recurrence_up
    )
))
supp6_score_table_total[, analysis := factor(analysis, levels = rev(c("Unique peptides", "All peptides")))]
fwrite(
    supp6_score_table_total,
    file.path(supp6_source_dir, "SupplementaryFigure6D_peptide_score_sensitivity_source.csv")
)
p6d <- ggplot(supp6_score_table_total, aes(x = mean_delta, y = analysis, color = analysis))
p6d <- p6d + geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4)
p6d <- p6d + geom_segment(aes(x = ci_low, xend = ci_high, yend = analysis), linewidth = 0.95)
p6d <- p6d + geom_point(size = 3.1)
p6d <- p6d + geom_text(
    aes(x = pmax(mean_delta, ci_high) + 0.06, label = sprintf("%d/%d up\nP %.2g", n_recurrence_up, n_pairs, fdr)),
    hjust = 0, size = 2.7, color = "black"
)
p6d <- p6d + scale_color_manual(
    values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"]))
)
p6d <- p6d + coord_cartesian(xlim = c(-0.02, max(supp6_score_table_total$ci_high) + 0.25), clip = "off")
p6d <- p6d + labs(x = "Mean paired raw20 score delta", y = NULL, color = NULL)
p6d <- p6d + theme_main
p6d <- p6d + theme(legend.position = "none")
save_plot_pair(p6d, "supp6d_peptide_score_sensitivity", supp6_plot_dir, width = 4.8, height = 2.5)
supp6_panel_4_table <- panel_manifest(
    "SupplementaryFigure6D",
    "supp6d_peptide_score_sensitivity",
    "Score-scale sensitivity for unique versus all peptides.",
    supp6_plot_dir
)

# In supplementary panel E, report NES and nominal p values calculated for the five original modes.
supp6_gsea_table_total <- rbindlist(list(
    data.table(
        dataset  = "PDC unique",
        NES      = pdc$unique$score_tests$gsea_NES,
        p        = pdc$unique$score_tests$gsea_p,
        detected = pdc$unique$coverage$n_detected,
        pairs    = pdc$unique$score_tests$n_pairs
    ),
    data.table(
        dataset  = "PDC all",
        NES      = pdc$all$score_tests$gsea_NES,
        p        = pdc$all$score_tests$gsea_p,
        detected = pdc$all$coverage$n_detected,
        pairs    = pdc$all$score_tests$n_pairs
    ),
    data.table(
        dataset = "MSV canonical",
        NES = msv_table_total[analysis == "canonical_regions",
            gsea_NES],
        p = msv_table_total[analysis == "canonical_regions",
            gsea_p],
        detected = msv_table_total[analysis == "canonical_regions",
            n_detected],
        pairs = msv_table_total[analysis == "canonical_regions",
            n_pairs]
    ),
    data.table(
        dataset = "MSV averaged",
        NES = msv_table_total[analysis == "averaged_repeated_regions",
            gsea_NES],
        p = msv_table_total[analysis == "averaged_repeated_regions",
            gsea_p],
        detected = msv_table_total[analysis == "averaged_repeated_regions",
            n_detected],
        pairs = msv_table_total[analysis == "averaged_repeated_regions",
            n_pairs]
    ),
    data.table(
        dataset  = "Zenodo",
        NES      = zenodo_table$gsea_NES,
        p        = zenodo_table$gsea_p,
        detected = zenodo_table$n_detected,
        pairs    = zenodo_table$n_pairs
    )
))
supp6_gsea_table_total[, dataset := factor(dataset, levels = rev(dataset))]
fwrite(supp6_gsea_table_total, file.path(supp6_source_dir, "SupplementaryFigure6E_gsea_sensitivity_source.csv"))
p6e <- ggplot(supp6_gsea_table_total, aes(x = NES, y = dataset, color = dataset))
p6e <- p6e + geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4)
p6e <- p6e + geom_point(size = 3.1)
p6e <- p6e + geom_text(
    aes(label = sprintf("nominal P %.2g  |  %d/20  |  n=%d", p, detected, pairs)),
    nudge_x = 0.09, hjust = 0, size = 2.7, color = "black"
)
p6e <- p6e + scale_color_manual(values = c(
    "PDC unique" = unname(npg["red"]),
    "PDC all" = unname(npg["teal"]),
    "MSV canonical" = unname(npg["blue"]),
    "MSV averaged" = unname(npg["orange"]),
    "Zenodo" = unname(npg["navy"])
))
p6e <- p6e + coord_cartesian(
    xlim = c(min(supp6_gsea_table_total$NES) - 0.12, max(supp6_gsea_table_total$NES) + 1.0),
    clip = "off"
)
p6e <- p6e + labs(x = "GSEA NES", y = NULL, color = NULL)
p6e <- p6e + theme_main
p6e <- p6e + theme(legend.position = "none")
save_plot_pair(p6e, "supp6e_gsea_sensitivity", supp6_plot_dir, width = 5.5, height = 3.1)
supp6_panel_5_table <- panel_manifest(
    "SupplementaryFigure6E",
    "supp6e_gsea_sensitivity",
    "GSEA-scale sensitivity across proteome datasets and peptide handling.",
    supp6_plot_dir
)

supp6_input_table <- data.table(
    input = c(
        "PDC000514 raw TMT matrix",
        "MSV000087947 processed paired proteome matrix",
        "Zenodo7646550 processed paired proteome matrix"
    ),
    location = c(pdc_source_dir, msv_source_dir, zenodo_source_root)
)
supp6_write_dir <- dirname(supp6_source_dir)
writeLines(c(
    "# Supplementary Figure 6 summary",
    "",
    sprintf(
        "- PDC unique/all both measure %d/20 raw20 proteins; main inference remains unique-peptide PDC.",
        pdc$unique$coverage$n_detected
    ),
    sprintf(
        "- MSV sensitivity is limited to %d/20 proteins and Zenodo to %d/20 proteins.",
        msv_table_total[analysis == "canonical_regions",
            n_detected],
        zenodo_table$n_detected
    ),
    sprintf("- GSEA sensitivity is shown on its own NES axis and is not mixed with score delta.")
), file.path(supp6_write_dir, "SupplementaryFigure6_summary.md"), useBytes = TRUE)
supp6_manifest_table_total <- rbindlist(list(
    supp6_panel_1_table, supp6_panel_2_table, supp6_panel_3_table, supp6_panel_4_table, supp6_panel_5_table
))
fwrite(supp6_manifest_table_total, file.path(supp6_write_dir, "SupplementaryFigure6_panel_manifest.csv"))
fwrite(supp6_input_table, file.path(supp6_write_dir, "SupplementaryFigure6_inputs.csv"))

# Verify that the main and supplementary manifests refer to nonempty figures.
plot_paths <- c(
    figure4_manifest_table_total$pdf, figure4_manifest_table_total$png,
    supp6_manifest_table_total$pdf, supp6_manifest_table_total$png
)
stopifnot(all(file.exists(plot_paths)), all(file.info(plot_paths)$size > 0))
cat("panels=", nrow(figure4_manifest_table_total) + nrow(supp6_manifest_table_total), "\n", sep = "")
