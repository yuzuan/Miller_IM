#!/usr/bin/env Rscript
# 目的：把 GSE274546 的基因×细胞 RDS 计数转为既有稀疏缓存，不改变整数化规则。
# 输入：命令行指定的 RDS；矩阵行名为基因、列名为细胞条码。
# 输出：指定目录中的 matrix.mtx.gz、genes.tsv.gz、barcodes.tsv.gz；文件名是公开接口。
# 运行：Rscript scripts/03a_GSE274546矩阵转换.R <input.RDS> <output_dir>。
# 保留 R/Matrix 实现和 round 后转整数的顺序，避免改变已有 Python 下游的复现入口。

suppressPackageStartupMessages(library(Matrix))

# 检查两个命令行参数；输入路径由调用者提供，不写入本机绝对路径。
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 03a_GSE274546矩阵转换.R <input.RDS> <output_dir>")
input_rds <- args[1]
count_matrix <- readRDS(input_rds)
if (!is.matrix(count_matrix) && !inherits(count_matrix, "Matrix")) {
    stop("Unsupported RDS object class: ", paste(class(count_matrix), collapse = ", "))
}

# 原矩阵类型决定转换方式；稀疏矩阵的行列顺序保持不变。
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

# Matrix summary 的每行是一个非零元素；沿用原有 round→integer，不更改计数口径。
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

# 基因和条码分别逐行写出；长度与矩阵对应维度一致才能供下游读取。
stopifnot(length(gene_names) == nrow(sparse_matrix), length(cell_names) == ncol(sparse_matrix))
genes_path <- file.path(output_dir, "genes.tsv.gz")
genes_connection <- gzfile(genes_path, open = "wt")
writeLines(gene_names, genes_connection)
close(genes_connection)

barcodes_path <- file.path(output_dir, "barcodes.tsv.gz")
barcodes_connection <- gzfile(barcodes_path, open = "wt")
writeLines(cell_names, barcodes_connection)
close(barcodes_connection)

# 文件必须实际存在且非空；本核对不读取或修改原始 RDS。
output_paths <- c(matrix_path, genes_path, barcodes_path)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
cat("genes=", length(gene_names), " cells=", length(cell_names), " nonzero=", nrow(triplet_table), "\n", sep = "")
