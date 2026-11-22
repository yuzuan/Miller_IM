#!/usr/bin/env Rscript
# Purpose: Apply the frozen GSE274546 per-library QC thresholds and classify doublets with scDblFinder.
# Inputs under MILLER_IM_DATA_ROOT: results/manifests/gse274546_standard_libraries.csv and the GSE274546 cache.
# Outputs under write/36_gse274546_raw_independent_reannotation/01_qc_doublets:
# scDblFinder_cell_classifications.csv.gz, capture_qc_doublet_summary.csv,
# patient_condition_qc_summary.csv and qc_run_summary.csv with their published column names.
# Run with Rscript scripts/03b_gse274546_quality_control.R [project_root].
# Retain the R/Bioconductor workflow, per-library seed placement and SerialParam for the original classifications.

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(SingleCellExperiment)
    library(scDblFinder)
    library(BiocParallel)
})

# Set the original global seed. Each library model below gets a separate seed.
set.seed(274546)
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1L) args[1L] else ".", mustWork = TRUE)
source_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root, "data"))
manifest_path <- file.path(source_root, "results", "manifests", "gse274546_standard_libraries.csv")
manifest_table <- as.data.frame(fread(manifest_path))
manifest_table <- manifest_table[manifest_table$progression %in% c("Primary", "1st Recurrent"), ]

# The frozen manifest specifies 111 libraries and 59 patients.
if (nrow(manifest_table) != 111L || uniqueN(manifest_table$patient_id) != 59L) {
    stop("Expected 111 Primary/1st Recurrent libraries from 59 patients.")
}
manifest_table$pair_id <- manifest_table$patient_id
manifest_table$condition <- fifelse(manifest_table$progression == "Primary", "Primary", "Recurrent")
manifest_order <- order(
    manifest_table$patient_id, manifest_table$progression_order, manifest_table$library_id,
    na.last = FALSE, method = "radix"
)
manifest_table <- manifest_table[manifest_order, ]

# Keep one result per library and process libraries sequentially in the original random-call order.
classification_rows <- vector("list", nrow(manifest_table))
summary_rows <- vector("list", nrow(manifest_table))
raw_root <- file.path(source_root, "GSE274546")
for (index in seq_len(nrow(manifest_table))) {
    library_table <- manifest_table[index, , drop = FALSE]
    library_id <- as.character(library_table$library_id)
    library_dir <- file.path(raw_root, library_id)
    source_paths <- file.path(library_dir, c("matrix.mtx.gz", "genes.tsv.gz", "barcodes.tsv.gz"))
    if (any(!file.exists(source_paths))) stop("Missing source files for ", library_id)

    # The headerless gene and barcode TSVs identify the rows and columns of the gene-by-cell matrix.
    gene_path <- file.path(library_dir, "genes.tsv.gz")
    gene_table <- as.data.frame(fread(cmd = paste("gzip -cd", shQuote(gene_path)), header = FALSE))
    genes <- as.character(gene_table[, 1L])
    barcode_path <- file.path(library_dir, "barcodes.tsv.gz")
    barcode_table <- as.data.frame(fread(cmd = paste("gzip -cd", shQuote(barcode_path)), header = FALSE))
    barcodes <- as.character(barcode_table[, 1L])
    matrix_path <- file.path(library_dir, "matrix.mtx.gz")
    connection <- gzfile(matrix_path, open = "rt")
    count_matrix <- readMM(connection)
    close(connection)
    count_matrix <- as(count_matrix, "CsparseMatrix")
    if (!identical(dim(count_matrix), c(length(genes), length(barcodes)))) {
        stop("Matrix dimension mismatch for ", library_id)
    }
    rownames(count_matrix) <- make.unique(genes)
    colnames(count_matrix) <- barcodes

    # Require at least 200 detected genes and at most 3% mitochondrial counts with the original zero-total rule.
    symbols <- sub("\\.[0-9]+$", "", rownames(count_matrix))
    is_mt <- startsWith(toupper(symbols), "MT-")
    n_genes <- colSums(count_matrix > 0)
    total_counts <- colSums(count_matrix)
    mt_counts <- rep(0, ncol(count_matrix))
    if (any(is_mt)) mt_counts <- colSums(count_matrix[is_mt, , drop = FALSE])
    pct_mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)
    qc_keep <- n_genes >= 200L & pct_mt <= 3
    n_qc <- sum(qc_keep)

    # The minimum is 50 passing cells per capture. Keep only the summary for excluded libraries.
    if (n_qc < 50L) {
        summary_rows[[index]] <- data.frame(
            capture_id   = library_id,
            library_id   = library_id,
            gsm_id       = library_table$gsm_id,
            patient_id   = library_table$patient_id,
            pair_id      = library_table$pair_id,
            progression  = library_table$progression,
            condition    = library_table$condition,
            n_raw        = ncol(count_matrix),
            n_qc         = n_qc,
            n_doublet    = NA_integer_,
            n_singlet    = 0L,
            doublet_rate = NA_real_,
            status       = "excluded_too_few_qc_cells"
        )
        rm(count_matrix)
        gc(verbose = FALSE)
        next
    }

    # The per-library seed reset and serial BiocParallel preserve the original doublet classifications.
    qc_matrix <- count_matrix[, qc_keep, drop = FALSE]
    sce <- SingleCellExperiment(list(counts = qc_matrix))
    set.seed(274546 + sum(utf8ToInt(library_id)))
    sce <- scDblFinder(sce, BPPARAM = SerialParam(progressbar = FALSE))
    cell_class <- as.character(colData(sce)$scDblFinder.class)
    score <- as.numeric(colData(sce)$scDblFinder.score)
    if (anyNA(cell_class) || !all(cell_class %in% c("singlet", "doublet"))) {
        stop("scDblFinder returned invalid classes for ", library_id)
    }

    # Each row is one passing cell, with the published scDblFinder fields and types.
    classification_rows[[index]] <- data.frame(
        capture_id        = library_id,
        library_id        = library_id,
        gsm_id            = as.character(library_table$gsm_id),
        patient_id        = as.character(library_table$patient_id),
        pair_id           = as.character(library_table$pair_id),
        progression       = as.character(library_table$progression),
        condition         = as.character(library_table$condition),
        barcode           = colnames(qc_matrix),
        n_genes_by_counts = as.integer(n_genes[qc_keep]),
        total_counts      = as.integer(total_counts[qc_keep]),
        pct_counts_mt     = as.numeric(pct_mt[qc_keep]),
        scDblFinder_score = score,
        scDblFinder_class = cell_class
    )
    summary_rows[[index]] <- data.frame(
        capture_id   = library_id,
        library_id   = library_id,
        gsm_id       = library_table$gsm_id,
        patient_id   = library_table$patient_id,
        pair_id      = library_table$pair_id,
        progression  = library_table$progression,
        condition    = library_table$condition,
        n_raw        = ncol(count_matrix),
        n_qc         = n_qc,
        n_doublet    = sum(cell_class == "doublet"),
        n_singlet    = sum(cell_class == "singlet"),
        doublet_rate = mean(cell_class == "doublet"),
        status       = "ok"
    )
    rm(count_matrix, qc_matrix, sce)
    gc(verbose = FALSE)
}

# Check capture-barcode uniqueness after combining the results and retain the existing output names.
classification_table_total <- as.data.frame(rbindlist(classification_rows, use.names = TRUE, fill = TRUE))
summary_table_total <- as.data.frame(rbindlist(summary_rows, use.names = TRUE, fill = TRUE))
cell_ids <- paste(classification_table_total$capture_id, classification_table_total$barcode, sep = ":")
if (anyDuplicated(cell_ids)) stop("Duplicated capture-barcode classifications.")
out_dir <- file.path(root, "write", "36_gse274546_raw_independent_reannotation", "01_qc_doublets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(classification_table_total, file.path(out_dir, "scDblFinder_cell_classifications.csv.gz"), compress = "gzip")
fwrite(summary_table_total, file.path(out_dir, "capture_qc_doublet_summary.csv"))

# Each row is one patient-condition unit. Sum successful libraries in their order of first appearance.
success_table <- summary_table_total[summary_table_total$status == "ok", , drop = FALSE]
pair_table <- unique(success_table[, c("patient_id", "condition"), drop = FALSE])
pair_table$raw_cells <- 0L
pair_table$qc_cells <- 0L
pair_table$singlets <- 0L
for (index in seq_len(nrow(pair_table))) {
    selected <- success_table$patient_id == pair_table$patient_id[index] &
        success_table$condition == pair_table$condition[index]
    pair_table$raw_cells[index] <- sum(success_table$n_raw[selected])
    pair_table$qc_cells[index] <- sum(success_table$n_qc[selected])
    pair_table$singlets[index] <- sum(success_table$n_singlet[selected])
}
fwrite(pair_table, file.path(out_dir, "patient_condition_qc_summary.csv"))

# Compute the run summary from the library results, with NA doublet counts for excluded libraries.
run_summary_table <- data.frame(
    metric = c(
        "input_libraries", "input_patients", "raw_cells", "qc_cells", "doublets",
        "singlets", "successful_libraries", "excluded_libraries"
    ),
    value = c(
        nrow(manifest_table), uniqueN(manifest_table$patient_id), sum(summary_table_total$n_raw),
        sum(summary_table_total$n_qc), sum(summary_table_total$n_doublet, na.rm = TRUE),
        sum(summary_table_total$n_singlet), sum(summary_table_total$status == "ok"),
        sum(summary_table_total$status != "ok")
    )
)
fwrite(run_summary_table, file.path(out_dir, "qc_run_summary.csv"))

# Require nonempty outputs and one classification row per QC-passing cell from the successful libraries.
output_paths <- file.path(
    out_dir,
    c(
        "scDblFinder_cell_classifications.csv.gz", "capture_qc_doublet_summary.csv",
        "patient_condition_qc_summary.csv", "qc_run_summary.csv"
    )
)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
stopifnot(nrow(classification_table_total) == sum(success_table$n_qc))
cat("libraries=", nrow(manifest_table), " singlets=", sum(summary_table_total$n_singlet), "\n", sep = "")
