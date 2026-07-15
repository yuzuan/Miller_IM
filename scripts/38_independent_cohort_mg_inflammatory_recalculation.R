#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: 38_independent_cohort_mg_inflammatory_recalculation.R <output_dir>")
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(dirname(script_path))
output_dir <- normalizePath(args[[1]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

independent_dirs <- c(
  GSE174554 = file.path(project_dir, "write", "34_gse174554_raw_independent_discovery", "03_myeloid"),
  GSE274546 = file.path(project_dir, "write", "36_gse274546_raw_independent_reannotation", "03_myeloid")
)
miller_file <- file.path(project_dir, "config", "miller_im_gene_set.tsv")

required_files <- c(
  unlist(lapply(independent_dirs, function(x) c(
    file.path(x, "patient_condition_pseudobulk_counts.csv.gz"),
    file.path(x, "patient_condition_pseudobulk_metadata.csv")
  ))),
  miller_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = "; "))

read_counts <- function(path) {
  x <- read.csv(gzfile(path), row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  if (any(!is.finite(x)) || any(abs(x - round(x)) > 1e-8)) stop("Non-integer counts in ", path)
  storage.mode(x) <- "integer"
  if (anyDuplicated(rownames(x))) stop("Duplicated genes in ", path)
  if (any(x < 0L)) stop("Negative counts in ", path)
  x
}

load_independent_cohort <- function(dataset) {
  input_dir <- independent_dirs[[dataset]]
  count_path <- file.path(input_dir, "patient_condition_pseudobulk_counts.csv.gz")
  meta_path <- file.path(input_dir, "patient_condition_pseudobulk_metadata.csv")
  x <- read_counts(count_path)
  m <- read.csv(meta_path, stringsAsFactors = FALSE)
  required_meta <- c("sample_id", "patient_id", "condition", "n_cells")
  if (!all(required_meta %in% names(m))) stop("Missing metadata columns for ", dataset)
  if (!identical(colnames(x), m$sample_id)) stop(dataset, " counts and metadata are not aligned.")
  if (anyDuplicated(m$sample_id)) stop("Duplicated pseudobulk sample IDs for ", dataset)
  if (!all(m$condition %in% c("Primary", "Recurrent"))) stop("Unexpected condition for ", dataset)
  if (anyNA(m$patient_id) || any(!nzchar(m$patient_id))) {
    keep <- !is.na(m$patient_id) & nzchar(m$patient_id)
    x <- x[, keep, drop = FALSE]
    m <- m[keep, , drop = FALSE]
  }
  source_label <- if (dataset == "GSE174554") {
    "Step34 independently reconstructed and annotated clean myeloid"
  } else {
    "Step36 independently reconstructed and annotated clean myeloid"
  }
  meta <- data.frame(
    unit_id = m$sample_id,
    dataset = dataset,
    pair_key = m$patient_id,
    condition = m$condition,
    n_myeloid_cells = m$n_cells,
    source = source_label,
    count_path = normalizePath(count_path),
    metadata_path = normalizePath(meta_path),
    stringsAsFactors = FALSE
  )
  list(counts = x, meta = meta)
}

miller <- read.delim(miller_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(miller), c("rank", "gene")) || nrow(miller) != 20L) {
  stop("config/miller_im_gene_set.tsv must contain exactly rank and gene columns with 20 rows.")
}
miller <- miller[order(miller$rank), , drop = FALSE]
raw20 <- miller$gene

expected_raw20 <- c(
  "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
  "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
  "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)
if (!identical(raw20, expected_raw20)) stop("Miller raw top20 changed from the audited definition.")

fixed_sets <- list(
  Miller_Microglial_Inflammatory_raw_top20 = raw20
)

provenance <- data.frame(
  signature = names(fixed_sets),
  source = "First 20 genes of the Miller et al. Nature 2025 Supplementary Table 2 cNMF ranking",
  n_genes = lengths(fixed_sets),
  genes = vapply(fixed_sets, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(output_dir, "fixed_program_provenance.csv"), row.names = FALSE)

eligible_pairs <- function(meta, threshold) {
  m <- meta[meta$n_myeloid_cells >= threshold & meta$condition %in% c("Primary", "Recurrent"), , drop = FALSE]
  tab <- table(m$pair_key, m$condition)
  complete <- rownames(tab)[tab[, "Primary"] == 1 & tab[, "Recurrent"] == 1]
  m[m$pair_key %in% complete, , drop = FALSE]
}

sign_flip_test <- function(delta, seed) {
  delta <- delta[is.finite(delta)]
  n <- length(delta)
  observed <- abs(mean(delta))
  if (!n) return(c(p_value = NA_real_, n_permutations = 0, exact = FALSE))
  if (n <= 20) {
    total <- 2^n
    extreme <- 0
    chunk <- 20000
    powers <- 2^(seq_len(n) - 1)
    for (start in seq(0, total - 1, by = chunk)) {
      ids <- start:min(start + chunk - 1, total - 1)
      bits <- outer(ids, powers, function(a, b) bitwAnd(a, b) != 0)
      permuted <- rowMeans((2 * bits - 1) * matrix(delta, nrow = length(ids), ncol = n, byrow = TRUE))
      extreme <- extreme + sum(abs(permuted) >= observed - 1e-12)
    }
    return(c(p_value = extreme / total, n_permutations = total, exact = TRUE))
  }
  set.seed(seed)
  total <- 100000
  extreme <- 0
  chunk <- 5000
  for (start in seq(1, total, by = chunk)) {
    current <- min(chunk, total - start + 1)
    signs <- matrix(sample(c(-1, 1), current * n, replace = TRUE), nrow = current)
    permuted <- rowMeans(signs * matrix(delta, nrow = current, ncol = n, byrow = TRUE))
    extreme <- extreme + sum(abs(permuted) >= observed - 1e-12)
  }
  c(p_value = (extreme + 1) / (total + 1), n_permutations = total, exact = FALSE)
}

score_fixed_sets <- function(cohort, dataset, threshold) {
  y_all <- DGEList(counts = cohort$counts)
  y_all <- calcNormFactors(y_all)
  log_cpm <- cpm(y_all, log = TRUE, prior.count = 1)
  eligible <- eligible_pairs(cohort$meta, threshold)
  rows <- list()
  unit_rows <- list()
  for (set_name in names(fixed_sets)) {
    genes <- fixed_sets[[set_name]]
    present <- intersect(genes, rownames(log_cpm))
    z <- t(scale(t(log_cpm[present, , drop = FALSE])))
    z[!is.finite(z)] <- 0
    scores <- colMeans(z)
    u <- eligible
    u$threshold <- threshold
    u$signature <- set_name
    u$score <- scores[u$unit_id]
    u <- u[order(u$pair_key, factor(u$condition, levels = c("Primary", "Recurrent"))), ]
    wide <- reshape(u[, c("pair_key", "condition", "score")], idvar = "pair_key", timevar = "condition", direction = "wide")
    delta <- wide$score.Recurrent - wide$score.Primary
    test <- sign_flip_test(delta, 20260713 + threshold + match(dataset, names(independent_dirs)) * 100)
    rows[[set_name]] <- data.frame(
      dataset = dataset,
      threshold = threshold,
      formal_testable = nrow(wide) >= 10,
      signature = set_name,
      n_pairs = nrow(wide),
      n_input_genes = length(genes),
      n_present_genes = length(present),
      coverage = length(present) / length(genes),
      mean_delta = mean(delta),
      median_delta = median(delta),
      positive_pairs = sum(delta > 0),
      positive_fraction = mean(delta > 0),
      p_value = unname(test["p_value"]),
      sign_flip_exact = as.logical(test["exact"]),
      n_permutations = unname(test["n_permutations"]),
      stringsAsFactors = FALSE
    )
    u$delta_recurrent_minus_primary <- delta[match(u$pair_key, wide$pair_key)]
    unit_rows[[set_name]] <- u
  }
  list(summary = do.call(rbind, rows), units = do.call(rbind, unit_rows))
}

run_fgsea <- function(pathways, rank_stat) {
  tested <- lapply(pathways, function(genes) intersect(unique(genes), names(rank_stat)))
  tested <- tested[lengths(tested) >= 5]
  set.seed(42)
  out <- fgseaMultilevel(
    pathways = tested,
    stats = rank_stat,
    minSize = 5,
    maxSize = 500,
    eps = 0,
    BPPARAM = SerialParam(progressbar = FALSE)
  )
  out <- as.data.frame(out)
  out$leadingEdge <- vapply(out$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  out
}

run_pseudobulk <- function(cohort, dataset, threshold) {
  m <- eligible_pairs(cohort$meta, threshold)
  x <- cohort$counts[, m$unit_id, drop = FALSE]
  m <- m[match(colnames(x), m$unit_id), , drop = FALSE]
  m$condition <- factor(m$condition, levels = c("Primary", "Recurrent"))
  m$pair_key <- factor(m$pair_key)
  n_pairs <- length(unique(m$pair_key))
  if (n_pairs < 5) stop(dataset, " has fewer than five pairs at threshold ", threshold)
  design <- model.matrix(~ pair_key + condition, data = m)
  y <- DGEList(counts = x)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  coef_index <- which(colnames(design) == "conditionRecurrent")
  if (length(coef_index) != 1) stop("Cannot find recurrent coefficient for ", dataset)
  qlf <- glmQLFTest(fit, coef = coef_index)
  deg <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  deg$gene <- rownames(deg)
  deg <- deg[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  deg$dataset <- dataset
  deg$threshold <- threshold
  deg$formal_testable <- n_pairs >= 10
  rank_stat <- sign(deg$logFC) * sqrt(pmax(deg$F, 0))
  names(rank_stat) <- deg$gene
  rank_stat <- sort(rank_stat[is.finite(rank_stat)], decreasing = TRUE)

  targeted <- run_fgsea(fixed_sets, rank_stat)
  targeted$padj <- NULL
  targeted$dataset <- dataset
  targeted$threshold <- threshold
  targeted$n_pairs <- n_pairs
  targeted$formal_testable <- n_pairs >= 10

  gene_rows <- lapply(names(fixed_sets), function(set_name) {
    genes <- fixed_sets[[set_name]]
    part <- deg[match(genes, deg$gene), , drop = FALSE]
    data.frame(
      dataset = dataset,
      threshold = threshold,
      signature = set_name,
      gene = genes,
      tested = genes %in% deg$gene,
      logFC = part$logFC,
      PValue = part$PValue,
      FDR = part$FDR,
      stringsAsFactors = FALSE
    )
  })
  list(deg = deg, targeted = targeted, genes = do.call(rbind, gene_rows))
}

cohorts <- lapply(names(independent_dirs), load_independent_cohort)
names(cohorts) <- names(independent_dirs)
thresholds <- c(20, 50)
pair_audit <- list()
score_summary <- list()
score_units <- list()
deg_results <- list()
targeted_results <- list()
gene_results <- list()

for (dataset in names(cohorts)) {
  cohort <- cohorts[[dataset]]
  for (threshold in thresholds) {
    eligible <- eligible_pairs(cohort$meta, threshold)
    n_pairs <- length(unique(eligible$pair_key))
    pair_audit[[paste(dataset, threshold)]] <- data.frame(
      dataset = dataset,
      threshold = threshold,
      n_pairs = n_pairs,
      n_units = nrow(eligible),
      n_cells_primary = sum(eligible$n_myeloid_cells[eligible$condition == "Primary"]),
      n_cells_recurrent = sum(eligible$n_myeloid_cells[eligible$condition == "Recurrent"]),
      formal_testable = n_pairs >= 10,
      source = unique(cohort$meta$source),
      count_path = unique(cohort$meta$count_path),
      metadata_path = unique(cohort$meta$metadata_path),
      stringsAsFactors = FALSE
    )
    scores <- score_fixed_sets(cohort, dataset, threshold)
    score_summary[[paste(dataset, threshold)]] <- scores$summary
    score_units[[paste(dataset, threshold)]] <- scores$units
    model <- run_pseudobulk(cohort, dataset, threshold)
    deg_results[[paste(dataset, threshold)]] <- model$deg
    targeted_results[[paste(dataset, threshold)]] <- model$targeted
    gene_results[[paste(dataset, threshold)]] <- model$genes
  }
}

pair_audit <- do.call(rbind, pair_audit)
score_summary <- do.call(rbind, score_summary)
score_units <- do.call(rbind, score_units)
deg_results <- do.call(rbind, deg_results)
targeted_results <- do.call(rbind, targeted_results)
gene_results <- do.call(rbind, gene_results)

write.csv(pair_audit, file.path(output_dir, "independent_input_pair_audit.csv"), row.names = FALSE)
write.csv(score_summary, file.path(output_dir, "independent_fixed_program_paired_scores.csv"), row.names = FALSE)
write.csv(score_units, file.path(output_dir, "independent_fixed_program_scores_by_unit.csv"), row.names = FALSE)
write.csv(deg_results, file.path(output_dir, "independent_paired_edger_all_genes.csv"), row.names = FALSE)
write.csv(targeted_results, file.path(output_dir, "independent_fixed_program_targeted_gsea.csv"), row.names = FALSE)
write.csv(gene_results, file.path(output_dir, "independent_fixed_program_gene_direction.csv"), row.names = FALSE)

main_targeted <- targeted_results[targeted_results$threshold == 20, ]
main_raw <- main_targeted[main_targeted$pathway == "Miller_Microglial_Inflammatory_raw_top20", ]
strict_raw <- nrow(main_raw) == 2 && all(main_raw$NES > 0) && all(main_raw$pval < 0.05)

summary_lines <- c(
  "# Independently reconstructed cohort Mg-inflammatory recalculation",
  "",
  "## Input rule",
  "",
  "- GSE174554 uses only the Step34 independently reconstructed and annotated clean-myeloid raw-count pseudobulk.",
  "- GSE274546 uses only the Step36 independently reconstructed and annotated clean-myeloid raw-count pseudobulk.",
  "- No merged GSE174554/GSE274546 object contributes cells, annotations or counts to these tests.",
  "",
  "## Main paired result at >=20 clean myeloid cells per endpoint",
  "",
  "| dataset | signature | pairs | NES | nominal p |",
  "|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(main_targeted))) {
  row <- main_targeted[i, ]
  summary_lines <- c(summary_lines, sprintf(
    "| %s | %s | %d | %.3f | %.4g |",
    row$dataset, row$pathway, row$n_pairs, row$NES, row$pval
  ))
}
summary_lines <- c(
  summary_lines,
  "",
  "## Decision",
  "",
  paste0("- The prespecified Miller top20 program is positive with nominal P<0.05 in both independent cohorts: ", strict_raw, "."),
  "- Exact NES, pair counts and leading-edge genes are allowed to change because the independently reconstructed clean-myeloid entries are different.",
  "- A positive answer means the biological conclusion is reproduced, not that every numerical value is identical.",
  "",
  "## Statistical method",
  "",
  "- edgeR paired raw-count pseudobulk model: ~ patient + condition.",
  "- GSEA ranking: sign(logFC) x sqrt(F); fgseaMultilevel seed 42; serial execution.",
  "- The >=20-cell threshold is the formal main analysis. The >=50-cell GSE174554 result has fewer than 10 pairs and is descriptive."
)
writeLines(summary_lines, file.path(output_dir, "FINAL_INDEPENDENT_RECALCULATION.md"), useBytes = TRUE)

cat("STEP38_INDEPENDENT_MG_INFLAMMATORY_RECALCULATION_COMPLETE\n")
print(pair_audit)
print(main_targeted[, c("dataset", "pathway", "n_pairs", "NES", "pval", "size")])
