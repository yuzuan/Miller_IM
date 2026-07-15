#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

set.seed(174554)

root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else ".", mustWork = TRUE)
source_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root, "data"))
manifest_path <- file.path(source_root, "results", "manifests", "gse174554_libraries.csv")
raw_root <- file.path(source_root, "GSE174554", "GSE174554_RAW")
out_dir <- file.path(root, "write", "34_gse174554_raw_independent_discovery", "01_qc_doublets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_GENES <- 200L
MAX_MT <- 3
MIN_QC_CELLS_PER_CAPTURE <- 50L
RESEQUENCING_OVERLAP_MIN <- 0.80
INDEPENDENT_OVERLAP_MAX <- 0.05
FORCE_MANUAL_MATRIX_LIBRARY <- "GSM5319529_SF8963"
FIXED_MANUAL_MATRIX <- file.path(out_dir, "source_repairs", "GSM5319529_SF8963_matrix.fixed.mtx.gz")

manifest <- fread(manifest_path)
stopifnot(nrow(manifest) == 91L, uniqueN(manifest$gsm_id) == 81L)
manifest <- manifest[progression %in% c("Primary", "1st Recurrent")]
for (column in c("matrix_path", "features_path", "barcodes_path")) {
  manifest[, (column) := file.path(raw_root, basename(get(column)))]
  missing_source <- manifest[!file.exists(get(column)), get(column)]
  if (length(missing_source)) stop("Missing relocated source files: ", paste(missing_source, collapse = ", "))
}
setorder(manifest, gsm_id, tech_rep)

read_features <- function(path) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(path)), header = FALSE)
  if (ncol(tab) == 1L) {
    symbols <- as.character(tab[[1L]])
    return(data.table(gene_id = symbols, gene_symbol = symbols))
  }
  data.table(gene_id = as.character(tab[[1L]]), gene_symbol = as.character(tab[[2L]]))
}

read_barcodes <- function(path) {
  as.character(fread(cmd = paste("gzip -cd", shQuote(path)), header = FALSE)[[1L]])
}

read_library <- function(row) {
  features <- read_features(row$features_path)
  barcodes <- read_barcodes(row$barcodes_path)
  if (as.character(row$library_id) == FORCE_MANUAL_MATRIX_LIBRARY) {
    if (!file.exists(FIXED_MANUAL_MATRIX)) stop("Missing repaired SF8963 matrix: ", FIXED_MANUAL_MATRIX)
    con <- gzfile(FIXED_MANUAL_MATRIX, open = "rt")
    counts <- readMM(con)
    close(con)
    counts <- as(counts, "CsparseMatrix")
  } else {
    con <- gzfile(row$matrix_path, open = "rt")
    counts <- readMM(con)
    close(con)
    counts <- as(counts, "CsparseMatrix")
  }
  if (!identical(dim(counts), c(nrow(features), length(barcodes)))) {
    stop("Matrix dimension mismatch: ", row$library_id)
  }
  rownames(counts) <- make.unique(features$gene_symbol)
  colnames(counts) <- barcodes
  list(counts = counts, features = features, barcodes = barcodes)
}

merge_resequencing <- function(main, batch2, capture_id) {
  if (!identical(main$features, batch2$features)) stop("Feature order differs for ", capture_id)
  union_barcodes <- union(main$barcodes, batch2$barcodes)
  joined <- cbind(main$counts, batch2$counts)
  source_barcodes <- c(main$barcodes, batch2$barcodes)
  map <- sparseMatrix(
    i = seq_along(source_barcodes),
    j = match(source_barcodes, union_barcodes),
    x = 1,
    dims = c(length(source_barcodes), length(union_barcodes))
  )
  merged <- as(joined %*% map, "dgCMatrix")
  colnames(merged) <- union_barcodes
  rownames(merged) <- rownames(main$counts)
  list(counts = merged, features = main$features, barcodes = union_barcodes)
}

classification_rows <- list()
summary_rows <- list()
process_capture <- function(capture_id, row, capture_data, merge_mode) {
  counts <- capture_data$counts
  symbols <- sub("\\.[0-9]+$", "", rownames(counts))
  mt <- startsWith(toupper(symbols), "MT-")
  n_genes <- Matrix::colSums(counts > 0)
  total_counts <- Matrix::colSums(counts)
  mt_counts <- if (any(mt)) Matrix::colSums(counts[mt, , drop = FALSE]) else rep(0, ncol(counts))
  pct_mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)
  qc_keep <- n_genes >= MIN_GENES & pct_mt <= MAX_MT
  n_qc <- sum(qc_keep)
  if (n_qc < MIN_QC_CELLS_PER_CAPTURE) {
    summary <- data.table(
      capture_id = capture_id, gsm_id = row$gsm_id, patient_id = row$patient_id,
      pair_id = row$pair_id, progression = row$progression, merge_mode = merge_mode,
      n_raw = ncol(counts), n_qc = n_qc, n_doublet = NA_integer_, n_singlet = 0L,
      doublet_rate = NA_real_, status = "excluded_too_few_qc_cells"
    )
    return(list(classification = NULL, summary = summary))
  }
  qc_counts <- counts[, qc_keep, drop = FALSE]
  sce <- SingleCellExperiment(list(counts = qc_counts))
  set.seed(174554 + sum(utf8ToInt(capture_id)))
  sce <- scDblFinder(sce, BPPARAM = SerialParam(progressbar = FALSE))
  cls <- as.character(colData(sce)$scDblFinder.class)
  score <- as.numeric(colData(sce)$scDblFinder.score)
  if (anyNA(cls) || !all(cls %in% c("singlet", "doublet"))) stop("scDblFinder returned invalid classes for ", capture_id)
  classification <- data.table(
    capture_id = capture_id,
    gsm_id = as.character(row$gsm_id),
    patient_id = as.character(row$patient_id),
    pair_id = as.character(row$pair_id),
    progression = as.character(row$progression),
    barcode = colnames(qc_counts),
    n_genes_by_counts = as.integer(n_genes[qc_keep]),
    total_counts = as.integer(total_counts[qc_keep]),
    pct_counts_mt = as.numeric(pct_mt[qc_keep]),
    scDblFinder_score = score,
    scDblFinder_class = cls
  )
  summary <- data.table(
    capture_id = capture_id, gsm_id = row$gsm_id, patient_id = row$patient_id,
    pair_id = row$pair_id, progression = row$progression, merge_mode = merge_mode,
    n_raw = ncol(counts), n_qc = n_qc, n_doublet = sum(cls == "doublet"),
    n_singlet = sum(cls == "singlet"), doublet_rate = mean(cls == "doublet"), status = "ok"
  )
  message(capture_id, ": raw=", summary$n_raw, " qc=", n_qc, " singlet=", summary$n_singlet)
  list(classification = classification, summary = summary)
}

overlap_rows <- list()
constructed_captures <- 0L
for (gsm in unique(manifest$gsm_id)) {
  rows <- manifest[gsm_id == gsm]
  if (nrow(rows) == 1L) {
    capture_id <- as.character(rows$library_id[[1L]])
    result <- process_capture(capture_id, rows[1L], read_library(rows[1L]), "single_library")
    classification_rows[[capture_id]] <- result$classification
    summary_rows[[capture_id]] <- result$summary
    constructed_captures <- constructed_captures + 1L
    rm(result); gc(verbose = FALSE)
    next
  }
  if (nrow(rows) != 2L || !setequal(rows$tech_rep, c("main", "batch2"))) stop("Unexpected technical replicate structure for ", gsm)
  main_row <- rows[tech_rep == "main"][1L]
  batch_row <- rows[tech_rep == "batch2"][1L]
  main <- read_library(main_row)
  batch2 <- read_library(batch_row)
  if (!identical(main$features, batch2$features)) stop("Feature tables differ within ", gsm)
  overlap <- intersect(main$barcodes, batch2$barcodes)
  overlap_fraction_min <- length(overlap) / min(length(main$barcodes), length(batch2$barcodes))
  if (overlap_fraction_min >= RESEQUENCING_OVERLAP_MIN) {
    capture_id <- paste0(gsm, "__merged_resequencing")
    merged <- merge_resequencing(main, batch2, capture_id)
    result <- process_capture(capture_id, main_row, merged, "merged_resequencing")
    classification_rows[[capture_id]] <- result$classification
    summary_rows[[capture_id]] <- result$summary
    constructed_captures <- constructed_captures + 1L
    decision <- "merge_umi_by_barcode"
    rm(merged, result)
  } else if (overlap_fraction_min <= INDEPENDENT_OVERLAP_MAX) {
    for (entry in list(list(row = main_row, data = main), list(row = batch_row, data = batch2))) {
      capture_id <- as.character(entry$row$library_id)
      result <- process_capture(capture_id, entry$row, entry$data, "independent_capture")
      classification_rows[[capture_id]] <- result$classification
      summary_rows[[capture_id]] <- result$summary
      constructed_captures <- constructed_captures + 1L
      rm(result)
    }
    decision <- "retain_as_independent_captures"
  } else {
    stop("Ambiguous barcode overlap for ", gsm, ": ", signif(overlap_fraction_min, 4))
  }
  overlap_rows[[gsm]] <- data.table(
    gsm_id = gsm, sample_short = as.character(main_row$sample_short),
    n_main = length(main$barcodes), n_batch2 = length(batch2$barcodes),
    n_overlap = length(overlap), overlap_fraction_min = overlap_fraction_min, decision = decision
  )
  rm(main, batch2); gc(verbose = FALSE)
}

overlap_audit <- rbindlist(overlap_rows, fill = TRUE)
fwrite(overlap_audit, file.path(out_dir, "batch2_barcode_overlap_audit.csv"))

classification <- rbindlist(classification_rows, use.names = TRUE)
summary <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
fwrite(classification, file.path(out_dir, "scDblFinder_cell_classifications.csv.gz"), compress = "gzip")
fwrite(summary, file.path(out_dir, "capture_qc_doublet_summary.csv"))

run_summary <- data.table(
  metric = c("input_matrices", "input_gsm", "constructed_captures", "raw_cells_after_resequencing_merge", "qc_cells", "singlets"),
  value = c(nrow(manifest), uniqueN(manifest$gsm_id), constructed_captures, sum(summary$n_raw), sum(summary$n_qc), sum(summary$n_singlet))
)
fwrite(run_summary, file.path(out_dir, "qc_run_summary.csv"))
cat(sprintf("GSE174554_QC_COMPLETE captures=%d singlets=%d\n", constructed_captures, sum(summary$n_singlet)))
