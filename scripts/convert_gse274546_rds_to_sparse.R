#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: convert_gse274546_rds_to_sparse.R <input.RDS> <output_dir>")
}

input_rds <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!is.matrix(obj) && !inherits(obj, "Matrix")) {
  stop(sprintf("Unsupported RDS object class: %s", paste(class(obj), collapse = ", ")))
}

if (is.matrix(obj)) {
  sparse_obj <- Matrix(obj, sparse = TRUE)
} else if (inherits(obj, "dgCMatrix")) {
  sparse_obj <- obj
} else {
  sparse_obj <- as(obj, "dgCMatrix")
}

gene_names <- rownames(obj)
cell_names <- colnames(obj)
if (is.null(gene_names) || is.null(cell_names)) {
  stop("Input matrix must have both rownames (genes) and colnames (cells).")
}

triplets <- summary(sparse_obj)
matrix_path <- file.path(output_dir, "matrix.mtx.gz")
genes_path <- file.path(output_dir, "genes.tsv.gz")
barcodes_path <- file.path(output_dir, "barcodes.tsv.gz")

matrix_con <- gzfile(matrix_path, open = "wt")
writeLines("%%MatrixMarket matrix coordinate integer general", matrix_con)
writeLines("% generated from GSE274546 dense RDS counts", matrix_con)
writeLines(sprintf("%d %d %d", nrow(sparse_obj), ncol(sparse_obj), nrow(triplets)), matrix_con)
utils::write.table(
  data.frame(i = triplets$i, j = triplets$j, x = as.integer(round(triplets$x))),
  file = matrix_con,
  quote = FALSE,
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
)
close(matrix_con)

genes_con <- gzfile(genes_path, open = "wt")
writeLines(gene_names, genes_con)
close(genes_con)

barcodes_con <- gzfile(barcodes_path, open = "wt")
writeLines(cell_names, barcodes_con)
close(barcodes_con)

cat(sprintf("Wrote sparse cache to %s\n", output_dir))
