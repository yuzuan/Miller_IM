#!/usr/bin/env Rscript
# Purpose: Reproduce the fixed and estimated correlation CAMERA results for the two formal PD-1 contrasts.
# Inputs: GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx (Genes and sample counts);
#   GSE154795_GBM.AllCell.Integrated.Scaled.ClusterRes.0.1.rds.gz (RNA counts and metadata);
#   GSE154795_family.soft.gz (sample technical covariates); step2_locked_gene_sets.csv (set_id, gene);
#   legacy Figure 5 raw20 CAMERA results and frozen GSE154795_myeloid_sample_audit.csv
#   (ID, condition, IDH_status, n_myeloid_cells).
# Outputs under step2_pd1_camera/: camera_results.csv, sample_audit.csv, validation_checks.csv, sessionInfo.txt,
#   <dataset>_tested_genes.csv, <dataset>_design.csv, GSE154795_cluster_audit.csv and GSE154795_formal_samples.csv.
# Run: Rscript NOA_R1/scripts/02f_pd1_camera.R NOA_R1
# NOA_PUBLIC_DATA_ROOT holds the source workbook, RDS and SOFT file. NOA_LEGACY_ROOT holds the Figure 5 checks.

suppressPackageStartupMessages({
    library(Matrix)
    library(edgeR)
    library(limma)
})

# Read the released fixed definitions. The same target was used for both exposure contrasts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 02f_pd1_camera.R <revision_dir>")
revision <- normalizePath(args[[1]], mustWork = TRUE)
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
definition_table <- read.csv(definition_path)
signature <- definition_table$gene[definition_table$set_id == "target_top20"]
stopifnot(length(signature) == 20L)

# Resolve the original public input paths and the prior analyses used to check sample membership and fixed p values.
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
public_root <- Sys.getenv("NOA_PUBLIC_DATA_ROOT")
if (!nzchar(public_root)) public_root <- Sys.getenv("MILLER_IM_DATA_ROOT")
if (!nzchar(public_root)) public_root <- file.path(project, "data")
public_root <- normalizePath(path.expand(public_root), winslash = "/", mustWork = FALSE)

gse121810_path <- file.path(
    public_root, "write", "28_supportive_dataset_search", "source_metadata",
    "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx"
)
gse154795_dir <- file.path(public_root, "write", "26_pure_bioinformatics_dataset_rescue", "source_metadata")
gse154795_path <- file.path(gse154795_dir, "GSE154795_GBM.AllCell.Integrated.Scaled.ClusterRes.0.1.rds.gz")
soft_path <- file.path(gse154795_dir, "GSE154795_family.soft.gz")
input_paths <- c(gse121810_path, gse154795_path, soft_path)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0L) stop("Missing PD-1 source files: ", paste(missing_paths, collapse = ", "))

output_dir <- file.path(revision, "step2_pd1_camera")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
reference_dir <- file.path(project, "write", "41_mg_inflammatory_sci_rebuild", "Figure5")

# Retain the existing validation output and record the first failed check before stopping.
checks <- list()

check <- function(dataset, label, condition, detail = "") {
    passed <- isTRUE(condition)
    checks[[length(checks) + 1L]] <<- data.frame(dataset = dataset, check = label, passed = passed, detail = detail)
    validation_table_total <- do.call(rbind, checks)
    write.csv(validation_table_total, file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
    if (!passed) stop(dataset, ": ", label, " ", detail)
}

# GSE121810 workbook rows are gene symbols. Sum duplicate symbols in the original order.
counts_table <- as.data.frame(readxl::read_excel(gse121810_path), check.names = FALSE)
rownames(counts_table) <- counts_table$Genes
counts_table$Genes <- NULL
counts_matrix <- as.matrix(counts_table)
counts_matrix <- rowsum(counts_matrix, group = rownames(counts_table), reorder = FALSE)
storage.mode(counts_matrix) <- "numeric"

# The workbook suffix _A marks neoadjuvant exposure. Retain the original 5/15 filter and unpaired group design.
group <- factor(ifelse(grepl("_A$", colnames(counts_matrix)), "Neo", "Adj"), levels = c("Adj", "Neo"))
design <- model.matrix(~ group)
filtered_counts <- DGEList(counts_matrix)
keep <- filterByExpr(filtered_counts, design = design, min.count = 5, min.total.count = 15)
filtered_counts <- calcNormFactors(filtered_counts[keep, , keep.lib.sizes = FALSE])
voom_fit <- voom(filtered_counts, design, plot = FALSE)
indices <- ids2indices(list(Miller_raw20 = signature), rownames(voom_fit))
check(
    "GSE121810",
    "full_rank_and_19_genes",
    qr(design)$rank == ncol(design) && length(indices$Miller_raw20) == 19L
)

# First match the original fixed-correlation result, retaining its exact comparison and rounding checks.
fixed_table <- camera(voom_fit, indices, design, contrast = "groupNeo", inter.gene.cor = 0.01)
reference_table <- read.csv(file.path(reference_dir, "GSE121810_raw20_camera.csv"))
gse121810_reference_p <- reference_table$PValue[reference_table$pathway == "Miller_raw20"]
check(
    "GSE121810",
    "original_fixed_P_reproduced",
    abs(fixed_table$PValue - gse121810_reference_p) <= 1e-10 &&
        round(fixed_table$PValue, 7) == round(2.81295477233785e-5, 7),
    paste(fixed_table$PValue, gse121810_reference_p)
)
fixed_tables <- list(GSE121810 = fixed_table)
models <- list(GSE121810 = list(voom = voom_fit, design = design, contrast = "groupNeo"))
sample_tables <- list()
sample_tables[["GSE121810"]] <- data.frame(
    dataset = "GSE121810",
    sample  = colnames(voom_fit),
    group   = as.character(group),
    n_myeloid_cells = NA_integer_
)

tested_gene_table <- data.frame(gene = rownames(voom_fit))
design_table <- data.frame(sample = colnames(voom_fit), design, check.names = FALSE)
write.csv(tested_gene_table, file.path(output_dir, "GSE121810_tested_genes.csv"), row.names = FALSE, na = "")
write.csv(design_table, file.path(output_dir, "GSE121810_design.csv"), row.names = FALSE, na = "")
rm(counts_table, counts_matrix, filtered_counts, voom_fit, design, group, fixed_table)
invisible(gc(verbose = FALSE))

# The source RDS may be gzip-wrapped. Keep the original order of the three read methods and close connections.
read_pd1_rds <- function(path) {
    result <- tryCatch(readRDS(path), error = identity)
    if (!inherits(result, "error")) return(result)

    connection <- gzfile(path, "rb")
    result <- tryCatch(readRDS(connection), error = identity, finally = close(connection))
    if (!inherits(result, "error")) return(result)

    connection <- gzcon(gzfile(path, "rb"))
    result <- tryCatch(readRDS(connection), error = identity, finally = close(connection))
    if (inherits(result, "error")) stop("Cannot read RDS: ", conditionMessage(result))
    return(result)
}

seurat_object <- read_pd1_rds(gse154795_path)
counts_matrix <- seurat_object@assays$RNA@counts
cell_table <- seurat_object@meta.data
rm(seurat_object)
invisible(gc(verbose = FALSE))
check(
    "GSE154795",
    "count_meta_alignment",
    ncol(counts_matrix) == nrow(cell_table) && identical(colnames(counts_matrix), rownames(cell_table))
)
cell_table$cluster <- as.character(cell_table$seurat_clusters)

# Group cell counts by cluster and compute the original six lineage scores from cluster logCPM.
lineages <- list(
    Myeloid      = c("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "CSF1R", "CD68", "C1QA", "C1QB", "C1QC"),
    T_NK         = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1"),
    B_cell       = c("MS4A1", "CD79A", "CD79B", "CD37", "CD22", "BANK1"),
    Plasma       = c("MZB1", "JCHAIN", "SDC1", "IGHG1", "IGHG3"),
    Neural_glial = c("GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "OLIG1", "OLIG2"),
    Erythroid    = c("HBA1", "HBA2", "HBB", "ALAS2", "GYPA")
)
cluster_factor <- factor(cell_table$cluster, levels = sort(unique(cell_table$cluster)))
cluster_design <- sparse.model.matrix(~ 0 + cluster_factor)
colnames(cluster_design) <- levels(cluster_factor)
cluster_counts <- counts_matrix %*% cluster_design
cluster_log_cpm <- cpm(cluster_counts, log = TRUE, prior.count = 1)
lineage_scores <- matrix(
    NA_real_,
    nrow = ncol(cluster_log_cpm),
    ncol = length(lineages),
    dimnames = list(colnames(cluster_log_cpm), names(lineages))
)
for (lineage in names(lineages)) {
    present <- intersect(lineages[[lineage]], rownames(cluster_log_cpm))
    lineage_scores[, lineage] <- colMeans(cluster_log_cpm[present, , drop = FALSE])
}

# A cluster needs at least 100 cells for the standardized lineage comparison and strict myeloid call.
cluster_sizes <- as.integer(table(cluster_factor)[rownames(lineage_scores)])
substantial <- cluster_sizes >= 100
standardized_scores <- scale(lineage_scores[substantial, , drop = FALSE])
best_lineage <- colnames(lineage_scores)[max.col(lineage_scores, ties.method = "first")]
best_lineage[substantial] <- colnames(standardized_scores)[max.col(standardized_scores, ties.method = "first")]
cluster_table <- data.frame(
    cluster         = rownames(lineage_scores),
    n_cells         = cluster_sizes,
    best_lineage    = best_lineage,
    is_auto_myeloid = best_lineage == "Myeloid" & cluster_sizes >= 100
)
selected_clusters <- cluster_table$cluster[cluster_table$is_auto_myeloid]
write.csv(cluster_table, file.path(output_dir, "GSE154795_cluster_audit.csv"), row.names = FALSE, na = "")
check("GSE154795", "strict_clusters_match_frozen", setequal(selected_clusters, c("0", "8")),
    paste(selected_clusters, collapse = ","))

# Aggregate the selected myeloid cells by source sample ID in encounter order.
myeloid_cells <- cell_table$cluster %in% selected_clusters
myeloid_cell_table <- cell_table[myeloid_cells, , drop = FALSE]
myeloid_counts <- counts_matrix[, myeloid_cells, drop = FALSE]
rm(counts_matrix, cell_table, cluster_counts, cluster_log_cpm, lineage_scores, standardized_scores, cluster_design)
invisible(gc(verbose = FALSE))
sample_factor <- factor(myeloid_cell_table$ID, levels = unique(myeloid_cell_table$ID))
sample_design <- sparse.model.matrix(~ 0 + sample_factor)
colnames(sample_design) <- levels(sample_factor)
pseudobulk_matrix <- myeloid_counts %*% sample_design
sample_cell_counts <- as.integer(table(sample_factor)[colnames(pseudobulk_matrix)])
rm(myeloid_counts, sample_design)
invisible(gc(verbose = FALSE))

# Sample labels come from the source RDS. Reorder after the SOFT join to keep matrix columns aligned.
sample_table <- unique(myeloid_cell_table[, c("ID", "condition", "IDH_status", "MGMT")])
rownames(sample_table) <- sample_table$ID
sample_table <- sample_table[colnames(pseudobulk_matrix), , drop = FALSE]
sample_table$n_myeloid_cells <- sample_cell_counts

# Each SOFT sample record contributes the first source occurrence of each technical covariate.
soft_lines <- readLines(gzfile(soft_path), warn = FALSE)
sample_starts <- grep("^\\^SAMPLE = ", soft_lines)
sample_ends <- c(sample_starts[-1] - 1, length(soft_lines))
field_patterns <- c(
    ID                = "^!Sample_title = ",
    GSM               = "^\\^SAMPLE = ",
    sequencing_method = "^!Sample_characteristics_ch1 = sequencing_method: ",
    gel_bead_version  = "^!Sample_characteristics_ch1 = gel_bead_version: "
)
technical_tables <- vector("list", length(sample_starts))
for (i in seq_along(sample_starts)) {
    block <- soft_lines[sample_starts[i]:sample_ends[i]]
    annotations <- setNames(rep(NA_character_, length(field_patterns)), names(field_patterns))
    for (field in names(field_patterns)) {
        pattern <- field_patterns[[field]]
        matches <- grep(pattern, block, value = TRUE)
        if (length(matches) > 0L) annotations[[field]] <- sub(pattern, "", matches[1])
    }
    gsm <- annotations[["GSM"]]
    technical_tables[[i]] <- data.frame(
        ID                = annotations[["ID"]],
        GSM               = gsm,
        sequencing_method = annotations[["sequencing_method"]],
        gel_bead_version  = annotations[["gel_bead_version"]],
        acquisition_wave  = ifelse(grepl("^GSM467", gsm), "early", "late"),
        stringsAsFactors  = FALSE
    )
}
technical_table_total <- unique(do.call(rbind, technical_tables))
sample_table_total <- merge(sample_table, technical_table_total, by = "ID", all.x = TRUE, sort = FALSE)
rownames(sample_table_total) <- sample_table_total$ID
sample_table_total <- sample_table_total[colnames(pseudobulk_matrix), , drop = FALSE]

# Keep the frozen untreated recurrent and PD-1 arms, wildtype IDH, and the 20-myeloid-cell sample threshold.
formal <- sample_table_total$condition %in% c("GBM.rec", "GBM.PD1") &
    sample_table_total$n_myeloid_cells >= 20 & sample_table_total$IDH_status == "WT"
formal_counts <- pseudobulk_matrix[, formal, drop = FALSE]
formal_table <- sample_table_total[formal, , drop = FALSE]
formal_table$condition <- factor(formal_table$condition, levels = c("GBM.rec", "GBM.PD1"))
frozen_path <- file.path(
    project, "write", "56_figure5_pd1_rebuild", "Figure5", "frozen_inputs", "GSE154795_myeloid_sample_audit.csv"
)
frozen_table <- read.csv(frozen_path)
frozen_rows <- frozen_table$condition %in% c("GBM.rec", "GBM.PD1") &
    frozen_table$n_myeloid_cells >= 20 & frozen_table$IDH_status == "WT"
frozen_table <- frozen_table[frozen_rows, ]
frozen_match <- match(formal_table$ID, frozen_table$ID)
check(
    "GSE154795",
    "formal_samples_match_frozen",
    setequal(formal_table$ID, frozen_table$ID) &&
        all(formal_table$n_myeloid_cells == frozen_table$n_myeloid_cells[frozen_match]),
    paste(nrow(formal_table), nrow(frozen_table))
)
write.csv(formal_table, file.path(output_dir, "GSE154795_formal_samples.csv"), row.names = FALSE, na = "")

# Preserve the original 5/10 gene filter, TMM and MGMT/technical-adjusted condition contrast.
filtered_counts <- calcNormFactors(DGEList(formal_counts))
keep <- filterByExpr(filtered_counts, group = formal_table$condition, min.count = 5, min.total.count = 10)
filtered_counts <- calcNormFactors(filtered_counts[keep, , keep.lib.sizes = FALSE])
for (field in c("MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave", "IDH_status")) {
    formal_table[, field] <- factor(formal_table[, field])
}
design <- model.matrix(~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition, formal_table)
check("GSE154795", "full_rank", qr(design)$rank == ncol(design))
coefficient <- grep("^condition", colnames(design))
check("GSE154795", "unique_condition_coef", length(coefficient) == 1L)
contrast <- colnames(design)[coefficient]
voom_fit <- voom(filtered_counts, design, plot = FALSE)
indices <- ids2indices(list(Miller_raw20 = signature), rownames(voom_fit))
check("GSE154795", "nineteen_genes", length(indices$Miller_raw20) == 19L)

# Reproduce the original fixed-correlation p value before adding the estimated-correlation sensitivity comparison.
fixed_table <- camera(voom_fit, indices, design, contrast = contrast, inter.gene.cor = 0.01)
reference_table <- read.csv(file.path(reference_dir, "GSE154795_raw20_camera_auto_strict_idhwt.csv"))
gse154795_reference_p <- reference_table$PValue[reference_table$pathway == "Miller_raw20"]
check(
    "GSE154795",
    "original_fixed_P_reproduced",
    abs(fixed_table$PValue - gse154795_reference_p) <= 1e-10 &&
        round(fixed_table$PValue, 4) == round(0.015367499497791, 4),
    paste(fixed_table$PValue, gse154795_reference_p)
)
fixed_tables[["GSE154795"]] <- fixed_table
models[["GSE154795"]] <- list(voom = voom_fit, design = design, contrast = contrast)
sample_tables[["GSE154795"]] <- data.frame(
    dataset = "GSE154795",
    sample  = colnames(voom_fit),
    group   = as.character(formal_table$condition),
    n_myeloid_cells = formal_table$n_myeloid_cells
)

tested_gene_table <- data.frame(gene = rownames(voom_fit))
design_table <- data.frame(sample = colnames(voom_fit), design, check.names = FALSE)
write.csv(tested_gene_table, file.path(output_dir, "GSE154795_tested_genes.csv"), row.names = FALSE, na = "")
write.csv(design_table, file.path(output_dir, "GSE154795_design.csv"), row.names = FALSE, na = "")
rm(pseudobulk_matrix, myeloid_cell_table, sample_table, filtered_counts, formal_counts, formal_table, soft_lines)
invisible(gc(verbose = FALSE))

# Both fixed p values must agree to 1e-10 before the estimated-correlation tests.
check(
    "both",
    "fixed_original_Ps_passed",
    abs(fixed_tables$GSE121810$PValue - gse121810_reference_p) <= 1e-10 &&
        abs(fixed_tables$GSE154795$PValue - gse154795_reference_p) <= 1e-10
)
camera_tables <- list()
for (dataset in c("GSE121810", "GSE154795")) {
    model <- models[[dataset]]
    design_rank <- qr(model$design)$rank
    residual_df <- nrow(model$design) - design_rank
    for (mode in c("estimated", "fixed_0.01")) {
        if (mode == "fixed_0.01") {
            camera_table <- fixed_tables[[dataset]]
            rho <- 0.01
        } else {
            indices <- ids2indices(list(Miller_raw20 = signature), rownames(model$voom))
            camera_table <- camera(
                model$voom,
                indices,
                model$design,
                contrast       = model$contrast,
                inter.gene.cor = NA,
                allow.neg.cor  = FALSE,
                use.ranks      = FALSE,
                directional    = TRUE,
                trend.var      = FALSE
            )
            rho <- camera_table$Correlation
        }

        # The reference DF differs between the estimated and fixed modes, as in the paired CAMERA scripts.
        camera_tables[[length(camera_tables) + 1L]] <- data.frame(
            dataset       = dataset,
            mode          = mode,
            contrast      = model$contrast,
            n_samples     = ncol(model$voom),
            n_genes       = camera_table$NGenes,
            tested_genes  = nrow(model$voom),
            design_rank   = design_rank,
            residual_df   = residual_df,
            camera_reference_df = if (mode == "estimated") residual_df else nrow(model$voom) - 2L,
            correlation   = rho,
            effective_vif = max(1, 1 + (camera_table$NGenes - 1) * rho),
            direction     = camera_table$Direction,
            p_value       = camera_table$PValue
        )
    }
}

# Write one row per cohort and mode, and one audit row per formal sample.
camera_table_total <- do.call(rbind, camera_tables)
check("both", "all_four_P_finite", nrow(camera_table_total) == 4L && all(is.finite(camera_table_total$p_value)))
sample_table_total <- do.call(rbind, sample_tables)
write.csv(camera_table_total, file.path(output_dir, "camera_results.csv"), row.names = FALSE, na = "")
write.csv(sample_table_total, file.path(output_dir, "sample_audit.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Validation checks passed: ", length(checks), "\n", sep = "")
