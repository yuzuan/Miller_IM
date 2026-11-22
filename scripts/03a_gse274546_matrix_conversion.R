#!/usr/bin/env Rscript
# Purpose: Convert the GSE274546 gene-by-cell count RDS to the existing sparse cache with the original rounding.
# Input: an RDS passed on the command line, with gene row names and cell barcode column names.
# Outputs in the requested directory: matrix.mtx.gz, genes.tsv.gz and barcodes.tsv.gz.
# Run with Rscript scripts/03a_gse274546_matrix_conversion.R <input.RDS> <output_dir>.
# Keep the R/Matrix conversion and round-then-integer order used by the downstream Python workflow.

suppressPackageStartupMessages(library(Matrix))

# Require the two command-line paths so callers can supply the source and output locations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 03a_gse274546_matrix_conversion.R <input.RDS> <output_dir>")
input_rds <- args[1]
count_matrix <- readRDS(input_rds)
if (!is.matrix(count_matrix) && !inherits(count_matrix, "Matrix")) {
    stop("Unsupported RDS object class: ", paste(class(count_matrix), collapse = ", "))
}

# Select the conversion by input matrix type while retaining gene and cell order.
if (is.matrix(count_matrix)) {
    sparse_matrix <- Matrix(count_matrix, sparse = TRUE)
} else if (inherits(count_matrix, "dgCMatrix")) {
    sparse_matrix <- count_matrix
} else {
    sparse_matrix <- as(count_matrix, "dgCMatrix")
}
gene_names <- rownames(count_matrix)
cell_names <- colnames(count_matrix)
if (is.null(gene_names) || is.null(cell_names)) {
    stop("Input matrix must have both rownames (genes) and colnames (cells).")
}

# Each Matrix summary row is one sparse entry. Apply the original rounding and integer conversion.
triplet_table <- as.data.frame(summary(sparse_matrix))
export_table <- data.frame(i = triplet_table$i, j = triplet_table$j, x = as.integer(round(triplet_table$x)))
stopifnot(nrow(export_table) == nrow(triplet_table), ncol(export_table) == 3L)
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
matrix_path <- file.path(output_dir, "matrix.mtx.gz")
matrix_connection <- gzfile(matrix_path, open = "wt")
writeLines("%%MatrixMarket matrix coordinate integer general", matrix_connection)
writeLines("% generated from GSE274546 dense RDS counts", matrix_connection)
writeLines(paste(nrow(sparse_matrix), ncol(sparse_matrix), nrow(triplet_table)), matrix_connection)
write.table(export_table, matrix_connection, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
close(matrix_connection)

# Write one gene or barcode per line. Their counts must match the corresponding matrix dimensions.
stopifnot(length(gene_names) == nrow(sparse_matrix), length(cell_names) == ncol(sparse_matrix))
genes_path <- file.path(output_dir, "genes.tsv.gz")
genes_connection <- gzfile(genes_path, open = "wt")
writeLines(gene_names, genes_connection)
close(genes_connection)

barcodes_path <- file.path(output_dir, "barcodes.tsv.gz")
barcodes_connection <- gzfile(barcodes_path, open = "wt")
writeLines(cell_names, barcodes_connection)
close(barcodes_connection)

# Check that the converted output files exist and are nonempty.
output_paths <- c(matrix_path, genes_path, barcodes_path)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
cat("genes=", length(gene_names), " cells=", length(cell_names), " nonzero=", nrow(triplet_table), "\n", sep = "")
