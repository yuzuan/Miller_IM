#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

set.seed(274546)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1L) args[[1L]] else ".", mustWork = TRUE)
source_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root, "data"))
manifest_path <- file.path(source_root, "results", "manifests", "gse274546_standard_libraries.csv")
raw_root <- file.path(source_root, "GSE274546")
out_dir <- file.path(root, "write", "36_gse274546_raw_independent_reannotation", "01_qc_doublets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_GENES <- 200L
MAX_MT <- 3
MIN_QC_CELLS_PER_CAPTURE <- 50L

manifest <- fread(manifest_path)
manifest <- manifest[progression %in% c("Primary", "1st Recurrent")]
if (nrow(manifest) != 111L || uniqueN(manifest$patient_id) != 59L) {
  stop("Expected 111 Primary/1st Recurrent libraries from 59 patients.")
}
manifest[, pair_id := patient_id]
manifest[, condition := fifelse(progression == "Primary", "Primary", "Recurrent")]
setorder(manifest, patient_id, progression_order, library_id)

source_paths <- function(library_id) {
  library_dir <- file.path(raw_root, library_id)
  paths <- c(
    matrix = file.path(library_dir, "matrix.mtx.gz"),
    genes = file.path(library_dir, "genes.tsv.gz"),
    barcodes = file.path(library_dir, "barcodes.tsv.gz")
  )
  if (any(!file.exists(paths))) stop("Missing source files for ", library_id)
  paths
}

read_library <- function(library_id) {
  paths <- source_paths(library_id)
  genes <- as.character(fread(cmd = paste("gzip -cd", shQuote(paths[["genes"]])), header = FALSE)[[1L]])
  barcodes <- as.character(fread(cmd = paste("gzip -cd", shQuote(paths[["barcodes"]])), header = FALSE)[[1L]])
  con <- gzfile(paths[["matrix"]], open = "rt")
  counts <- readMM(con)
  close(con)
  counts <- as(counts, "CsparseMatrix")
  if (!identical(dim(counts), c(length(genes), length(barcodes)))) {
    stop("Matrix dimension mismatch for ", library_id)
  }
  rownames(counts) <- make.unique(genes)
  colnames(counts) <- barcodes
  counts
}

classification_rows <- vector("list", nrow(manifest))
summary_rows <- vector("list", nrow(manifest))

for (index in seq_len(nrow(manifest))) {
  row <- manifest[index]
  library_id <- as.character(row$library_id)
  counts <- read_library(library_id)
  symbols <- sub("\\.[0-9]+$", "", rownames(counts))
  mt <- startsWith(toupper(symbols), "MT-")
  n_genes <- Matrix::colSums(counts > 0)
  total_counts <- Matrix::colSums(counts)
  mt_counts <- if (any(mt)) Matrix::colSums(counts[mt, , drop = FALSE]) else rep(0, ncol(counts))
  pct_mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)
  qc_keep <- n_genes >= MIN_GENES & pct_mt <= MAX_MT
  n_qc <- sum(qc_keep)

  if (n_qc < MIN_QC_CELLS_PER_CAPTURE) {
    summary_rows[[index]] <- data.table(
      capture_id = library_id, library_id = library_id, gsm_id = row$gsm_id,
      patient_id = row$patient_id, pair_id = row$pair_id, progression = row$progression,
      condition = row$condition, n_raw = ncol(counts), n_qc = n_qc,
      n_doublet = NA_integer_, n_singlet = 0L, doublet_rate = NA_real_,
      status = "excluded_too_few_qc_cells"
    )
    rm(counts); gc(verbose = FALSE)
    next
  }

  qc_counts <- counts[, qc_keep, drop = FALSE]
  sce <- SingleCellExperiment(list(counts = qc_counts))
  set.seed(274546 + sum(utf8ToInt(library_id)))
  sce <- scDblFinder(sce, BPPARAM = SerialParam(progressbar = FALSE))
  cls <- as.character(colData(sce)$scDblFinder.class)
  score <- as.numeric(colData(sce)$scDblFinder.score)
  if (anyNA(cls) || !all(cls %in% c("singlet", "doublet"))) {
    stop("scDblFinder returned invalid classes for ", library_id)
  }

  classification_rows[[index]] <- data.table(
    capture_id = library_id, library_id = library_id, gsm_id = as.character(row$gsm_id),
    patient_id = as.character(row$patient_id), pair_id = as.character(row$pair_id),
    progression = as.character(row$progression), condition = as.character(row$condition),
    barcode = colnames(qc_counts), n_genes_by_counts = as.integer(n_genes[qc_keep]),
    total_counts = as.integer(total_counts[qc_keep]), pct_counts_mt = as.numeric(pct_mt[qc_keep]),
    scDblFinder_score = score, scDblFinder_class = cls
  )
  summary_rows[[index]] <- data.table(
    capture_id = library_id, library_id = library_id, gsm_id = row$gsm_id,
    patient_id = row$patient_id, pair_id = row$pair_id, progression = row$progression,
    condition = row$condition, n_raw = ncol(counts), n_qc = n_qc,
    n_doublet = sum(cls == "doublet"), n_singlet = sum(cls == "singlet"),
    doublet_rate = mean(cls == "doublet"), status = "ok"
  )
  message(library_id, ": raw=", ncol(counts), " qc=", n_qc, " singlet=", sum(cls == "singlet"))
  rm(counts, qc_counts, sce); gc(verbose = FALSE)
}

classification <- rbindlist(classification_rows, use.names = TRUE, fill = TRUE)
summary <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
if (anyDuplicated(classification[, paste(capture_id, barcode, sep = ":")])) {
  stop("Duplicated capture-barcode classifications.")
}
fwrite(classification, file.path(out_dir, "scDblFinder_cell_classifications.csv.gz"), compress = "gzip")
fwrite(summary, file.path(out_dir, "capture_qc_doublet_summary.csv"))

pair_table <- summary[status == "ok", .(
  raw_cells = sum(n_raw), qc_cells = sum(n_qc), singlets = sum(n_singlet)
), by = .(patient_id, condition)]
fwrite(pair_table, file.path(out_dir, "patient_condition_qc_summary.csv"))

run_summary <- data.table(
  metric = c(
    "input_libraries", "input_patients", "raw_cells", "qc_cells", "doublets",
    "singlets", "successful_libraries", "excluded_libraries"
  ),
  value = c(
    nrow(manifest), uniqueN(manifest$patient_id), sum(summary$n_raw), sum(summary$n_qc),
    sum(summary$n_doublet, na.rm = TRUE), sum(summary$n_singlet),
    sum(summary$status == "ok"), sum(summary$status != "ok")
  )
)
fwrite(run_summary, file.path(out_dir, "qc_run_summary.csv"))
cat(sprintf("GSE274546_QC_COMPLETE libraries=%d singlets=%d\n", nrow(manifest), sum(summary$n_singlet)))
