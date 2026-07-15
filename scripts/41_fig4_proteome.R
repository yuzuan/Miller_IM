#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
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
if (length(script_arg) != 1) {
  stop("Cannot locate the current script.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(script_path)
project_dir <- dirname(project_dir)
data_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(project_dir, "data"))

write_root <- file.path(project_dir, "write", "41_mg_inflammatory_sci_rebuild")
figure_root <- file.path(project_dir, "figures", "41_mg_inflammatory_sci_rebuild")
figure4_write_dir <- file.path(write_root, "Figure4")
supp6_write_dir <- file.path(write_root, "SupplementaryFigure6")
figure4_plot_dir <- file.path(figure_root, "Figure4", "panel_library")
supp6_plot_dir <- file.path(figure_root, "SupplementaryFigure6", "panel_library")
figure4_source_dir <- file.path(figure4_write_dir, "source_data")
supp6_source_dir <- file.path(supp6_write_dir, "source_data")
for (path in c(
  write_root, figure_root, figure4_write_dir, supp6_write_dir,
  figure4_plot_dir, supp6_plot_dir, figure4_source_dir, supp6_source_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

npg <- c(
  red = "#E64B35",
  blue = "#4DBBD5",
  teal = "#00A087",
  navy = "#3C5488",
  orange = "#F39B7F",
  purple = "#8491B4",
  green = "#91D1C2",
  darkred = "#DC0000",
  grey = "#7A7A7A",
  lightgrey = "#D9D9D9"
)

theme_main <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

save_plot_pair <- function(plot, stem, dir_path, width, height) {
  pdf_path <- file.path(dir_path, paste0(stem, ".pdf"))
  ggsave(
    filename = pdf_path,
    plot = plot, width = width, height = height, device = cairo_pdf
  )
  ggsave(
    filename = file.path(dir_path, paste0(stem, ".png")),
    plot = plot, width = width, height = height, dpi = 600, bg = "white"
  )
  normalise_pdf_metadata(pdf_path)
}

normalise_pdf_metadata <- function(pdf_path) {
  tmp_path <- paste0(pdf_path, ".tmp")
  py_script <- tempfile(pattern = "normalise_pdf_", fileext = ".py")
  writeLines(c(
    "from pathlib import Path",
    "from pypdf import PdfReader, PdfWriter",
    "import sys",
    "",
    "src = Path(sys.argv[1])",
    "dst = Path(sys.argv[2])",
    "reader = PdfReader(str(src))",
    "writer = PdfWriter(clone_from=reader)",
    "writer.add_metadata({",
    "    '/Producer': 'Step41 recurrent GBM figure rebuild',",
    "    '/Creator': 'Codex',",
    "    '/CreationDate': \"D:20260713000000+00'00\",",
    "    '/ModDate': \"D:20260713000000+00'00\",",
    "})",
    "with dst.open('wb') as handle:",
    "    writer.write(handle)"
  ), py_script, useBytes = TRUE)
  on.exit(unlink(py_script), add = TRUE)
  status <- system2("python3", c(py_script, pdf_path, tmp_path), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(tmp_path)) {
    stop("Failed to normalise PDF metadata for ", pdf_path)
  }
  if (file.exists(pdf_path)) {
    unlink(pdf_path)
  }
  ok <- file.rename(tmp_path, pdf_path)
  if (!ok) {
    stop("Failed to replace PDF after metadata normalisation: ", pdf_path)
  }
}

locate_required <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop(label, " not found. Checked: ", paste(candidates, collapse = "; "))
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

paired_sign_flip <- function(delta, seed) {
  delta <- delta[is.finite(delta)]
  n <- length(delta)
  if (n < 3L) {
    return(NA_real_)
  }
  observed <- abs(mean(delta))
  if (n <= 20L) {
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
    null <- abs(as.vector(signs %*% delta) / n)
    return(mean(null >= observed - 1e-12))
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
  (extreme + 1) / (total + 1)
}

ci_mean <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) {
    return(c(mean = NA_real_, low = NA_real_, high = NA_real_))
  }
  se <- stats::sd(x) / sqrt(n)
  m <- mean(x)
  c(mean = m, low = m - 1.96 * se, high = m + 1.96 * se)
}

running_enrichment_curve <- function(stats, geneset) {
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  geneset <- intersect(geneset, names(stats))
  hits <- names(stats) %in% geneset
  n_hits <- sum(hits)
  n_total <- length(stats)
  if (n_hits < 3L || n_hits == n_total) {
    stop("Cannot build enrichment curve with current gene coverage.")
  }
  hit_weight <- abs(stats[hits])
  hit_scale <- sum(hit_weight)
  phit <- ifelse(hits, abs(stats) / hit_scale, 0)
  pmiss <- ifelse(!hits, 1 / (n_total - n_hits), 0)
  data.table(
    rank = seq_along(stats),
    gene = names(stats),
    rank_stat = unname(stats),
    hit = hits,
    running_enrichment = cumsum(phit - pmiss)
  )
}

half_violin_polygon <- function(values, center = 1, width = 0.26, side = "left") {
  values <- values[is.finite(values)]
  dens <- density(values, n = 256, bw = "nrd0")
  scaled <- if (max(dens$y) > 0) dens$y / max(dens$y) * width else rep(0, length(dens$y))
  x_side <- if (side == "left") center - scaled else center + scaled
  data.table(
    x = c(rep(center, length(dens$x)), rev(x_side)),
    y = c(dens$x, rev(dens$x))
  )
}

panel_manifest <- function(panel, stem, description, out_dir) {
  data.table(
    panel = panel,
    stem = stem,
    description = description,
    pdf = file.path(out_dir, paste0(stem, ".pdf")),
    png = file.path(out_dir, paste0(stem, ".png"))
  )
}

raw20 <- c(
  "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
  "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
  "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)

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
if (length(missing_local)) {
  stop("Missing local inputs: ", paste(missing_local, collapse = "; "))
}

pdc_source_dir <- locate_required(
  c(
    file.path(project_dir, "write", "31_proteomic_validation", "source_metadata", "PDC000514"),
    file.path(data_root, "write", "31_proteomic_validation", "source_metadata", "PDC000514")
  ),
  "PDC000514 source directory"
)
msv_source_dir <- locate_required(
  c(
    file.path(project_dir, "write", "31_proteomic_validation", "source_metadata", "MSV000087947"),
    file.path(data_root, "write", "31_proteomic_validation", "source_metadata", "MSV000087947")
  ),
  "MSV000087947 source directory"
)
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

analyse_pdc <- function() {
  matrix_path <- file.path(pdc_source_dir, "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv")
  raw <- fread(matrix_path, check.names = FALSE)
  raw <- raw[!Gene %in% c("Mean", "Median", "StdDev") & !is.na(Gene) & Gene != ""]

  parse_samples <- function(columns, prefix) {
    selected <- grep(paste0("^", prefix, " KNCC_GBM[0-9]+_T[12]:"), columns, value = TRUE)
    data.table(
      column = selected,
      sample = sub(paste0("^", prefix, " (KNCC_GBM[0-9]+_T[12]):.*$"), "\\1", selected)
    )[, `:=`(
      patient = sub("_T[12]$", "", sample),
      condition = fifelse(grepl("_T2$", sample), "Recurrent", "Primary")
    )]
  }

  analyse_matrix <- function(prefix, analysis, seed) {
    source_map <- parse_samples(names(raw), prefix)
    sample_map <- unique(source_map[, .(sample, patient, condition)])
    sample_map <- merge(
      sample_map,
      source_map[, .(n_source_columns = .N), by = sample],
      by = "sample", sort = FALSE
    )
    pair_counts <- dcast(sample_map, patient ~ condition, value.var = "sample", fun.aggregate = length)
    paired_patients <- pair_counts[Primary == 1L & Recurrent == 1L, patient]
    sample_map <- sample_map[patient %in% paired_patients]
    setorder(sample_map, patient, condition)

    source_map[, analysis := analysis]
    source_map[, included_pair := patient %in% paired_patients]

    source_columns <- source_map$column
    expr_source <- as.matrix(raw[, ..source_columns])
    storage.mode(expr_source) <- "double"
    rownames(expr_source) <- raw$Gene
    colnames(expr_source) <- source_map$sample
    expr_avg <- vapply(sample_map$sample, function(sample_id) {
      rowMeans(expr_source[, colnames(expr_source) == sample_id, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(expr_source)))
    expr_avg[!is.finite(expr_avg)] <- NA_real_
    rownames(expr_avg) <- raw$Gene
    colnames(expr_avg) <- sample_map$sample

    present <- intersect(raw20, rownames(expr_avg))
    coverage <- data.table(
      analysis = analysis,
      n_defined = length(raw20),
      n_detected = length(present),
      coverage = length(present) / length(raw20),
      detected_genes = paste(present, collapse = ";"),
      missing_genes = paste(setdiff(raw20, present), collapse = ";")
    )

    design <- model.matrix(~ factor(sample_map$patient) + factor(sample_map$condition, levels = c("Primary", "Recurrent")))
    coef_index <- ncol(design)
    fit <- eBayes(lmFit(expr_avg, design), robust = TRUE)
    gene_results <- as.data.table(
      topTable(fit, coef = coef_index, number = Inf, sort.by = "none"),
      keep.rownames = "gene"
    )
    gene_results[, analysis := analysis]
    gene_results[, measured_raw20 := gene %in% present]
    gene_results[, significant := adj.P.Val < 0.05]

    ranks <- gene_results$t
    names(ranks) <- gene_results$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    fgsea_out <- as.data.table(
      fgseaMultilevel(
        pathways = list(Miller_raw20 = present),
        stats = ranks,
        minSize = 3,
        maxSize = 500,
        eps = 0,
        BPPARAM = SerialParam(progressbar = FALSE)
      )
    )
    fgsea_out[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    fgsea_out[, analysis := analysis]
    gsea_curve <- running_enrichment_curve(ranks, present)
    gsea_curve[, analysis := analysis]

    z_expr <- t(scale(t(expr_avg[present, , drop = FALSE])))
    z_expr[!is.finite(z_expr)] <- NA_real_
    paired_scores <- rbindlist(lapply(unique(sample_map$patient), function(patient_id) {
      primary_sample <- sample_map[patient == patient_id & condition == "Primary", sample]
      recurrent_sample <- sample_map[patient == patient_id & condition == "Recurrent", sample]
      primary_values <- z_expr[, primary_sample, drop = TRUE]
      recurrent_values <- z_expr[, recurrent_sample, drop = TRUE]
      common <- is.finite(primary_values) & is.finite(recurrent_values)
      if (sum(common) < 3L) {
        return(NULL)
      }
      data.table(
        analysis = analysis,
        patient = patient_id,
        Primary = mean(primary_values[common]),
        Recurrent = mean(recurrent_values[common]),
        n_common_raw20 = sum(common),
        common_genes = paste(present[common], collapse = ";")
      )
    }))
    paired_scores[, delta := Recurrent - Primary]
    summary_ci <- ci_mean(paired_scores$delta)
    score_tests <- data.table(
      analysis = analysis,
      n_pairs = nrow(paired_scores),
      n_recurrence_up = sum(paired_scores$delta > 0),
      mean_delta = summary_ci[["mean"]],
      ci_low = summary_ci[["low"]],
      ci_high = summary_ci[["high"]],
      median_delta = median(paired_scores$delta),
      sign_flip_p = paired_sign_flip(paired_scores$delta, seed = seed),
      paired_t_p = t.test(paired_scores$delta)$p.value,
      wilcoxon_p = wilcox.test(paired_scores$delta, exact = FALSE)$p.value,
      n_defined = coverage$n_defined,
      n_detected = coverage$n_detected,
      coverage = coverage$coverage
    )
    score_tests[, gsea_NES := fgsea_out$NES]
    score_tests[, gsea_p := fgsea_out$pval]
    score_tests[, leading_edge := fgsea_out$leadingEdge]

    sample_flow <- data.table(
      analysis = analysis,
      total_patients = uniqueN(source_map$patient),
      primary_patients = uniqueN(source_map[condition == "Primary", patient]),
      recurrent_patients = uniqueN(source_map[condition == "Recurrent", patient]),
      paired_patients = uniqueN(sample_map$patient),
      primary_only_patients = uniqueN(source_map[condition == "Primary" & !included_pair, patient]),
      duplicated_samples = sum(source_map[, .N, by = sample]$N > 1L),
      extra_aliquot_columns = sum(pmax(source_map[, .N, by = sample]$N - 1L, 0L))
    )

    list(
      source_map = source_map,
      sample_map = sample_map,
      coverage = coverage,
      gene_results = gene_results,
      fgsea = fgsea_out,
      gsea_curve = gsea_curve,
      paired_scores = paired_scores,
      score_tests = score_tests,
      sample_flow = sample_flow
    )
  }

  unique_res <- analyse_matrix("Unshared Log", "unique_peptides", 20260713L)
  all_res <- analyse_matrix("Log", "all_peptides", 20260714L)
  list(unique = unique_res, all = all_res)
}

analyse_msv <- function() {
  matrix_path <- file.path(msv_source_dir, "401_2022_2506_MOESM3_ESM.xlsx")
  raw <- as.data.table(read_excel(matrix_path, sheet = "Imputed values", skip = 3))
  raw[, protein_row := .I]
  sample_columns <- grep("^GBM[0-9]+-(PM2?|RM2?)$", names(raw), value = TRUE)
  gene_map <- rbindlist(Map(function(index, genes) {
    data.table(protein_row = index, gene = trimws(unlist(strsplit(as.character(genes), ";", fixed = TRUE))))
  }, raw$protein_row, raw$`Gene.names`))
  gene_map <- gene_map[!is.na(gene) & gene != ""]
  expanded <- merge(gene_map, raw[, c("protein_row", sample_columns), with = FALSE], by = "protein_row")
  gene_matrix <- expanded[, lapply(.SD, median, na.rm = TRUE), by = gene, .SDcols = sample_columns]
  expr_all <- as.matrix(gene_matrix[, ..sample_columns])
  storage.mode(expr_all) <- "double"
  expr_all[!is.finite(expr_all)] <- NA_real_
  rownames(expr_all) <- gene_matrix$gene

  prepare <- function(mode) {
    if (mode == "canonical_regions") {
      selected <- grep("^GBM[0-9]+-(PM|RM)$", colnames(expr_all), value = TRUE)
      expr <- expr_all[, selected, drop = FALSE]
      sample_map <- data.table(sample = selected)
      sample_map[, `:=`(
        patient = sub("-(PM|RM)$", "", sample),
        condition = fifelse(grepl("-RM$", sample), "Recurrent", "Primary"),
        n_regions = 1L
      )]
    } else {
      src <- data.table(source_sample = colnames(expr_all))
      src[, `:=`(
        patient = sub("-(PM2?|RM2?)$", "", source_sample),
        condition = fifelse(grepl("-RM2?$", source_sample), "Recurrent", "Primary")
      )]
      sample_map <- unique(src[, .(patient, condition)])
      sample_map <- merge(
        sample_map, src[, .(n_regions = .N), by = .(patient, condition)],
        by = c("patient", "condition"), sort = FALSE
      )
      setorder(sample_map, patient, condition)
      sample_map[, sample := paste(patient, condition, sep = "_")]
      expr <- vapply(seq_len(nrow(sample_map)), function(i) {
        cols <- src[patient == sample_map$patient[i] & condition == sample_map$condition[i], source_sample]
        rowMeans(expr_all[, cols, drop = FALSE], na.rm = TRUE)
      }, numeric(nrow(expr_all)))
      rownames(expr) <- rownames(expr_all)
      colnames(expr) <- sample_map$sample
    }
    setorder(sample_map, patient, condition)
    list(expr = expr, sample_map = sample_map)
  }

  analyse <- function(mode, seed) {
    prepared <- prepare(mode)
    expr <- prepared$expr
    sample_map <- prepared$sample_map
    present <- intersect(raw20, rownames(expr))
    design <- model.matrix(~ factor(sample_map$patient) + factor(sample_map$condition, levels = c("Primary", "Recurrent")))
    coef_index <- ncol(design)
    fit <- eBayes(lmFit(expr, design), robust = TRUE)
    gene_results <- as.data.table(
      topTable(fit, coef = coef_index, number = Inf, sort.by = "none"),
      keep.rownames = "gene"
    )
    ranks <- gene_results$t
    names(ranks) <- gene_results$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    fgsea_out <- as.data.table(
      fgseaMultilevel(
        pathways = list(Miller_raw20 = present),
        stats = ranks,
        minSize = 3,
        maxSize = 500,
        eps = 0,
        BPPARAM = SerialParam(progressbar = FALSE)
      )
    )
    fgsea_out[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    z_expr <- t(scale(t(expr[present, , drop = FALSE])))
    z_expr[!is.finite(z_expr)] <- NA_real_
    paired_scores <- rbindlist(lapply(unique(sample_map$patient), function(patient_id) {
      primary_sample <- sample_map[patient == patient_id & condition == "Primary", sample]
      recurrent_sample <- sample_map[patient == patient_id & condition == "Recurrent", sample]
      primary_values <- z_expr[, primary_sample, drop = TRUE]
      recurrent_values <- z_expr[, recurrent_sample, drop = TRUE]
      common <- is.finite(primary_values) & is.finite(recurrent_values)
      if (sum(common) < 3L) {
        return(NULL)
      }
      data.table(
        analysis = mode,
        patient = patient_id,
        Primary = mean(primary_values[common]),
        Recurrent = mean(recurrent_values[common]),
        n_common_raw20 = sum(common)
      )
    }))
    paired_scores[, delta := Recurrent - Primary]
    ci <- ci_mean(paired_scores$delta)
    data.table(
      analysis = mode,
      n_pairs = nrow(paired_scores),
      n_recurrence_up = sum(paired_scores$delta > 0),
      mean_delta = ci[["mean"]],
      ci_low = ci[["low"]],
      ci_high = ci[["high"]],
      median_delta = median(paired_scores$delta),
      sign_flip_p = paired_sign_flip(paired_scores$delta, seed = seed),
      gsea_NES = fgsea_out$NES,
      gsea_p = fgsea_out$pval,
      n_detected = length(present),
      coverage = length(present) / length(raw20),
      detected_genes = paste(present, collapse = ";"),
      paired_scores = list(paired_scores)
    )
  }

  rbindlist(list(
    analyse("canonical_regions", 20260715L),
    analyse("averaged_repeated_regions", 20260716L)
  ), fill = TRUE)
}

analyse_zenodo <- function() {
  matrix_path <- file.path(
    zenodo_source_root, "data", "specific_search_fragpipe17",
    "specific_no_ptms_2", "tmt-report", "abundance_protein_MD.tsv"
  )
  limma_path <- file.path(zenodo_source_root, "results", "limma_proteomics_results.tsv")
  protein <- fread(matrix_path, check.names = FALSE)
  sample_columns <- grep("^[0-9]+-(prim|rec)$", names(protein), value = TRUE)
  protein <- protein[!is.na(Gene) & Gene != ""]
  gene_matrix <- protein[, lapply(.SD, function(x) median(x, na.rm = TRUE)), by = Gene, .SDcols = sample_columns]
  expr <- as.matrix(gene_matrix[, ..sample_columns])
  storage.mode(expr) <- "double"
  expr[!is.finite(expr)] <- NA_real_
  rownames(expr) <- gene_matrix$Gene
  present <- intersect(raw20, rownames(expr))

  z_expr <- t(scale(t(expr[present, , drop = FALSE])))
  z_expr[!is.finite(z_expr)] <- NA_real_
  paired_scores <- rbindlist(lapply(seq_len(length(sample_columns) / 2L), function(i) {
    patient_id <- as.character(i)
    primary_sample <- paste0(patient_id, "-prim")
    recurrent_sample <- paste0(patient_id, "-rec")
    primary_values <- z_expr[, primary_sample, drop = TRUE]
    recurrent_values <- z_expr[, recurrent_sample, drop = TRUE]
    common <- is.finite(primary_values) & is.finite(recurrent_values)
    if (sum(common) < 3L) {
      return(NULL)
    }
    data.table(
      patient = patient_id,
      Primary = mean(primary_values[common]),
      Recurrent = mean(recurrent_values[common]),
      n_common_raw20 = sum(common)
    )
  }))
  paired_scores[, delta := Recurrent - Primary]
  ci <- ci_mean(paired_scores$delta)

  limma_results <- fread(limma_path)
  limma_results <- limma_results[!is.na(Gene) & Gene != ""]
  limma_results[, abs_t := abs(t)]
  setorder(limma_results, Gene, -abs_t)
  limma_results <- limma_results[!duplicated(Gene)]
  ranks <- limma_results$t
  names(ranks) <- limma_results$Gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  fgsea_out <- as.data.table(
    fgseaMultilevel(
      pathways = list(Miller_raw20 = present),
      stats = ranks,
      minSize = 3,
      maxSize = 500,
      eps = 0,
      BPPARAM = SerialParam(progressbar = FALSE)
    )
  )
  fgsea_out[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]

  data.table(
    analysis = "single_matrix",
    n_pairs = nrow(paired_scores),
    n_recurrence_up = sum(paired_scores$delta > 0),
    mean_delta = ci[["mean"]],
    ci_low = ci[["low"]],
    ci_high = ci[["high"]],
    median_delta = median(paired_scores$delta),
    sign_flip_p = paired_sign_flip(paired_scores$delta, seed = 20260717L),
    gsea_NES = fgsea_out$NES,
    gsea_p = fgsea_out$pval,
    n_detected = length(present),
    coverage = length(present) / length(raw20),
    detected_genes = paste(present, collapse = ";"),
    paired_scores = list(paired_scores)
  )
}

pdc <- analyse_pdc()
msv <- analyse_msv()
zenodo <- analyse_zenodo()

step38_gsea <- fread(step38_gsea_path)
step38_gsea <- step38_gsea[
  pathway == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20
]
step38_gene <- fread(step38_gene_path)
step38_gene <- step38_gene[
  signature == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20 & tested == TRUE
]
geomx_summary <- fread(geomx_summary_path)
geomx_coverage <- fread(geomx_coverage_path)

pdc_gene_unique <- copy(pdc$unique$gene_results)[gene %in% raw20]
pdc_gene_all <- copy(pdc$all$gene_results)[gene %in% raw20]
pdc_gene_unique[, `:=`(
  ci_low = logFC - 1.96 * abs(logFC / t),
  ci_high = logFC + 1.96 * abs(logFC / t),
  analysis = "unique_peptides"
)]
pdc_gene_all[, `:=`(
  ci_low = logFC - 1.96 * abs(logFC / t),
  ci_high = logFC + 1.96 * abs(logFC / t),
  analysis = "all_peptides"
)]

figure4_panel_manifest <- rbindlist(list())
supp6_panel_manifest <- rbindlist(list())

figure4_delta <- copy(pdc$unique$paired_scores)
figure4_delta[, patient := factor(patient, levels = patient[order(delta)])]
delta_test <- pdc$unique$score_tests
delta_violin <- half_violin_polygon(figure4_delta$delta, center = 1, width = 0.25, side = "left")
set.seed(20260713)
figure4_delta[, jitter_x := 1.08 + runif(.N, -0.06, 0.06)]
fwrite(figure4_delta, file.path(figure4_source_dir, "Figure4A_patient_delta_source.csv"))
p4a <- ggplot() +
  geom_polygon(data = delta_violin, aes(x = x, y = y), fill = "#C8D8DC", alpha = 0.72, color = "#76A5AF", linewidth = 0.4) +
  geom_boxplot(
    data = figure4_delta, aes(x = 1, y = delta, group = 1),
    width = 0.10, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35
  ) +
  geom_point(
    data = figure4_delta, aes(x = jitter_x, y = delta),
    size = 1.15, shape = 21, stroke = 0.2, fill = npg["red"], color = "white", alpha = 0.78
  ) +
  geom_segment(
    aes(x = 0.80, xend = 0.80, y = delta_test$ci_low, yend = delta_test$ci_high),
    linewidth = 1.0, color = npg["red"]
  ) +
  geom_point(aes(x = 0.80, y = delta_test$mean_delta), size = 2.9, shape = 21, fill = npg["red"], color = "black", stroke = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4) +
  annotate(
    "text",
    x = 1.23,
    y = max(figure4_delta$delta) + 0.06,
    hjust = 1,
    vjust = 1,
    size = 2.95,
    label = sprintf(
      "mean %.3f\n95%% CI %.3f to %.3f\nn = %d\n%d/%d up\nP %.2g",
      delta_test$mean_delta, delta_test$ci_low, delta_test$ci_high,
      delta_test$n_pairs, delta_test$n_recurrence_up, delta_test$n_pairs, delta_test$sign_flip_p
    )
  ) +
  coord_cartesian(xlim = c(0.56, 1.24), ylim = c(min(figure4_delta$delta) - 0.12, max(figure4_delta$delta) + 0.12), clip = "off") +
  labs(
    x = NULL,
    y = "Recurrent - primary raw20 score",
    title = "PDC unique peptides",
    subtitle = "Patient-level raw20 delta distribution"
  ) +
  theme_main +
  theme(
    plot.title = element_text(size = 8.4, face = "bold"),
    plot.subtitle = element_text(size = 7.0, face = "plain"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
save_plot_pair(p4a, "fig4a_pdc_patient_delta", figure4_plot_dir, width = 3.25, height = 4.1)
figure4_panel_manifest <- rbindlist(list(
  figure4_panel_manifest,
  panel_manifest("Figure4A", "fig4a_pdc_patient_delta", "PDC unique-peptide patient delta distribution with mean and 95% CI.", figure4_plot_dir)
))

figure4_forest <- copy(pdc_gene_unique[measured_raw20 == TRUE])
setorder(figure4_forest, logFC)
figure4_forest[, gene := factor(gene, levels = gene)]
figure4_forest[, significant_label := factor(as.character(significant), levels = c("FALSE", "TRUE"))]
fwrite(figure4_forest, file.path(figure4_source_dir, "Figure4B_measurable_protein_forest_source.csv"))
p4b <- ggplot(figure4_forest, aes(x = logFC, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = gene), linewidth = 0.85, color = npg["grey"]) +
  geom_point(aes(fill = significant_label), shape = 21, size = 2.8, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"), labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
  labs(
    x = "Proteome logFC (recurrent vs primary)",
    y = NULL,
    fill = "FDR < 0.05",
    title = "Measured raw20 proteins in PDC unique peptides",
    subtitle = "Moderated logFC with 95% CI"
  ) +
  theme_main +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 8.2, face = "bold"),
    plot.subtitle = element_text(size = 7.0, face = "plain")
  )
save_plot_pair(p4b, "fig4b_measurable_protein_forest", figure4_plot_dir, width = 4.9, height = 4.3)
figure4_panel_manifest <- rbindlist(list(
  figure4_panel_manifest,
  panel_manifest("Figure4B", "fig4b_measurable_protein_forest", "Measured raw20 proteins in PDC unique-peptide analysis with 95% CI.", figure4_plot_dir)
))

figure4_gsea_curve <- copy(pdc$unique$gsea_curve)
figure4_gsea_curve[, leading_edge := gene %in% unlist(strsplit(pdc$unique$score_tests$leading_edge, ";", fixed = TRUE))]
fwrite(figure4_gsea_curve, file.path(figure4_source_dir, "Figure4C_protein_ranked_gsea_curve_source.csv"))
p4c_top <- ggplot(figure4_gsea_curve, aes(x = rank, y = running_enrichment)) +
  geom_line(color = npg["navy"], linewidth = 1.05) +
  geom_area(fill = npg["navy"], alpha = 0.12) +
  geom_hline(yintercept = 0, color = npg["grey"], linewidth = 0.4) +
  labs(y = "Running enrichment", x = NULL) +
  theme_main +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p4c_mid <- ggplot(figure4_gsea_curve[hit == TRUE], aes(x = rank, xend = rank, y = 0, yend = 1)) +
  geom_segment(color = "black", linewidth = 0.28) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_void()
p4c_bottom <- ggplot(figure4_gsea_curve, aes(x = rank, y = rank_stat)) +
  geom_area(data = figure4_gsea_curve[rank_stat >= 0], fill = npg["red"], alpha = 0.75) +
  geom_area(data = figure4_gsea_curve[rank_stat < 0], fill = npg["blue"], alpha = 0.75) +
  geom_hline(yintercept = 0, color = npg["grey"], linewidth = 0.35) +
  labs(x = "Ranked proteins", y = "Moderated t") +
  theme_main +
  theme(axis.title.y = element_text(size = 7))
p4c <- (p4c_top / p4c_mid / p4c_bottom) +
  plot_layout(heights = c(3.2, 0.45, 1.0)) +
  plot_annotation(
    title = "PDC unique peptides",
    subtitle = sprintf(
      "NES %.2f, nominal P %.2g, %d/%d proteins measured",
      pdc$unique$score_tests$gsea_NES, pdc$unique$score_tests$gsea_p,
      pdc$unique$score_tests$n_detected, length(raw20)
    )
  ) &
  theme(
    plot.title = element_text(size = 8.4, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.0, face = "plain", hjust = 0)
  )
save_plot_pair(p4c, "fig4c_protein_ranked_gsea", figure4_plot_dir, width = 4.9, height = 4.2)
figure4_panel_manifest <- rbindlist(list(
  figure4_panel_manifest,
  panel_manifest("Figure4C", "fig4c_protein_ranked_gsea", "Complete ranked-list GSEA curve for PDC unique-peptide proteome.", figure4_plot_dir)
))

geomx_gene <- fread(file.path(project_dir, "write", "39_independent_miller_mg_inflammatory_figures", "Figure3", "source_csv", "geomx_raw20_gene_level_summary.csv"))
geomx_dual <- geomx_gene[
  entry == "strict_pass_idhwt",
  .(gene, geomx_mean = mean_log2_delta, geomx_ci_low = ci_low, geomx_ci_high = ci_high, geomx_fdr = fdr)
]
pdc_dual <- pdc_gene_unique[
  measured_raw20 == TRUE,
  .(gene, proteome_mean = logFC, proteome_ci_low = ci_low, proteome_ci_high = ci_high, proteome_fdr = adj.P.Val)
]
dual_modal <- merge(geomx_dual, pdc_dual, by = "gene")
setorder(dual_modal, proteome_mean)
dual_modal[, gene := factor(gene, levels = gene)]
dual_modal[, `:=`(
  geomx_sig_label = factor(as.character(geomx_fdr < 0.05), levels = c("FALSE", "TRUE")),
  proteome_sig_label = factor(as.character(proteome_fdr < 0.05), levels = c("FALSE", "TRUE"))
)]
fwrite(dual_modal, file.path(figure4_source_dir, "Figure4D_cross_modal_dual_scale_source.csv"))
p4d_left <- ggplot(dual_modal, aes(x = geomx_mean, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.35) +
  geom_segment(aes(x = geomx_ci_low, xend = geomx_ci_high, yend = gene), linewidth = 0.8, color = npg["grey"]) +
  geom_point(aes(fill = geomx_sig_label), shape = 21, size = 2.7, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"), guide = "none") +
  labs(x = "GeoMx mean log2 delta", y = NULL, title = "IBA1+ GeoMx (22 pairs)") +
  theme_main +
  theme(plot.title = element_text(size = 8.0, face = "bold", color = npg["navy"]))
p4d_right <- ggplot(dual_modal, aes(x = proteome_mean, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.35) +
  geom_segment(aes(x = proteome_ci_low, xend = proteome_ci_high, yend = gene), linewidth = 0.8, color = npg["grey"]) +
  geom_point(aes(fill = proteome_sig_label), shape = 21, size = 2.7, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c(`TRUE` = unname(npg["red"]), `FALSE` = "white"), guide = "none") +
  labs(x = "Proteome logFC", y = NULL, title = "PDC unique peptides (105 pairs)") +
  theme_main +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 8.0, face = "bold", color = npg["darkred"])
  )
p4d <- p4d_left + p4d_right + plot_layout(widths = c(1, 1))
p4d <- p4d + plot_annotation(
  title = "Same measurable genes across GeoMx and proteome",
  subtitle = "Each column keeps its own effect scale"
) & theme(
  plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
  plot.subtitle = element_text(size = 7.0, face = "plain", hjust = 0)
)
save_plot_pair(p4d, "fig4d_cross_modal_dual_scale", figure4_plot_dir, width = 7.3, height = 4.3)
figure4_panel_manifest <- rbindlist(list(
  figure4_panel_manifest,
  panel_manifest("Figure4D", "fig4d_cross_modal_dual_scale", "Same-gene GeoMx and proteome effects shown in separate columns with separate scales.", figure4_plot_dir)
))

figure4_summary <- data.table(
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
    file.path(project_dir, "write", "39_independent_miller_mg_inflammatory_figures", "Figure3", "source_csv", "geomx_raw20_gene_level_summary.csv"),
    msv_source_dir,
    zenodo_source_root
  )
)
writeLines(c(
  "# Figure4 summary",
  "",
  sprintf("- Main protein result uses only Miller raw20 and PDC unique peptides: %d/%d proteins measured.", pdc$unique$score_tests$n_detected, length(raw20)),
  sprintf("- Patient delta: mean %.3f, 95%% CI %.3f to %.3f, %d/%d recurrence-up, P %.2g.", delta_test$mean_delta, delta_test$ci_low, delta_test$ci_high, delta_test$n_recurrence_up, delta_test$n_pairs, delta_test$sign_flip_p),
  sprintf("- Protein GSEA: NES %.3f, nominal P %.3g.", pdc$unique$score_tests$gsea_NES, pdc$unique$score_tests$gsea_p),
  sprintf("- Cross-modal dual-column panel is limited to %d genes measurable in both GeoMx and PDC.", nrow(dual_modal))
), file.path(figure4_write_dir, "Figure4_summary.md"), useBytes = TRUE)
fwrite(figure4_panel_manifest, file.path(figure4_write_dir, "Figure4_panel_manifest.csv"))
fwrite(figure4_summary, file.path(figure4_write_dir, "Figure4_inputs.csv"))

supp6_flow <- rbindlist(list(
  data.table(step = "Patients in matrix", value = pdc$unique$sample_flow$total_patients),
  data.table(step = "Primary available", value = pdc$unique$sample_flow$primary_patients),
  data.table(step = "Recurrent available", value = pdc$unique$sample_flow$recurrent_patients),
  data.table(step = "Paired patients", value = pdc$unique$sample_flow$paired_patients),
  data.table(step = "Primary-only excluded", value = pdc$unique$sample_flow$primary_only_patients),
  data.table(step = "Samples merged", value = pdc$unique$sample_flow$duplicated_samples),
  data.table(step = "Extra aliquot columns", value = pdc$unique$sample_flow$extra_aliquot_columns)
))
fwrite(supp6_flow, file.path(supp6_source_dir, "SupplementaryFigure6A_sample_flow_source.csv"))
supp6_boxes <- data.table(
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
p6a <- ggplot() +
  geom_rect(
    data = supp6_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = supp6_boxes$fill, color = "white", linewidth = 0.6
  ) +
  geom_text(
    data = supp6_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    color = "white", size = 3.2, fontface = "bold", lineheight = 0.95
  ) +
  annotate("segment", x = 1.2, xend = 1.8, y = 0.65, yend = 0.65, arrow = arrow(length = unit(0.10, "inches"))) +
  annotate("segment", x = 3.0, xend = 3.6, y = 0.65, yend = 0.65, arrow = arrow(length = unit(0.10, "inches"))) +
  annotate("segment", x = 4.8, xend = 5.4, y = 0.65, yend = 0.65, arrow = arrow(length = unit(0.10, "inches"))) +
  annotate("segment", x = 4.2, xend = 4.2, y = 0.15, yend = -0.15, arrow = arrow(length = unit(0.10, "inches"))) +
  annotate("segment", x = 6.0, xend = 6.0, y = 0.15, yend = -0.15, arrow = arrow(length = unit(0.10, "inches"))) +
  coord_cartesian(xlim = c(-0.3, 6.75), ylim = c(-1.28, 1.28), clip = "off") +
  theme_void()
save_plot_pair(p6a, "supp6a_sample_flow", supp6_plot_dir, width = 6.8, height = 2.5)
supp6_panel_manifest <- rbindlist(list(
  supp6_panel_manifest,
  panel_manifest("SupplementaryFigure6A", "supp6a_sample_flow", "PDC000514 sample flow and exclusions.", supp6_plot_dir)
))

supp6_coverage <- rbindlist(list(
  data.table(dataset = "PDC unique", pairs = pdc$unique$score_tests$n_pairs, detected = pdc$unique$coverage$n_detected, total = 20L, coverage = pdc$unique$coverage$coverage),
  data.table(dataset = "PDC all", pairs = pdc$all$score_tests$n_pairs, detected = pdc$all$coverage$n_detected, total = 20L, coverage = pdc$all$coverage$coverage),
  data.table(dataset = "MSV canonical", pairs = msv[analysis == "canonical_regions", n_pairs], detected = msv[analysis == "canonical_regions", n_detected], total = 20L, coverage = msv[analysis == "canonical_regions", coverage]),
  data.table(dataset = "MSV averaged", pairs = msv[analysis == "averaged_repeated_regions", n_pairs], detected = msv[analysis == "averaged_repeated_regions", n_detected], total = 20L, coverage = msv[analysis == "averaged_repeated_regions", coverage]),
  data.table(dataset = "Zenodo", pairs = zenodo$n_pairs, detected = zenodo$n_detected, total = 20L, coverage = zenodo$coverage)
))
supp6_coverage[, dataset := factor(dataset, levels = rev(dataset))]
fwrite(supp6_coverage, file.path(supp6_source_dir, "SupplementaryFigure6B_coverage_source.csv"))
p6b <- ggplot(supp6_coverage, aes(x = coverage, y = dataset, color = dataset)) +
  geom_segment(aes(x = 0, xend = coverage, yend = dataset), linewidth = 1.0, color = npg["lightgrey"]) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%d/20  |  n=%d", detected, pairs)), nudge_x = 0.04, hjust = 0, size = 2.75, color = "black") +
  scale_color_manual(values = c(
    "PDC unique" = unname(npg["red"]),
    "PDC all" = unname(npg["teal"]),
    "MSV canonical" = unname(npg["blue"]),
    "MSV averaged" = unname(npg["orange"]),
    "Zenodo" = unname(npg["navy"])
  )) +
  coord_cartesian(xlim = c(0, 0.92), clip = "off") +
  labs(x = "Raw20 protein coverage", y = NULL, color = NULL) +
  theme_main +
  theme(legend.position = "none")
save_plot_pair(p6b, "supp6b_dataset_coverage", supp6_plot_dir, width = 5.2, height = 3.0)
supp6_panel_manifest <- rbindlist(list(
  supp6_panel_manifest,
  panel_manifest("SupplementaryFigure6B", "supp6b_dataset_coverage", "Raw20 protein coverage across proteome datasets.", supp6_plot_dir)
))

supp6_counts <- rbindlist(list(
  copy(pdc$unique$paired_scores)[, analysis := "Unique peptides"],
  copy(pdc$all$paired_scores)[, analysis := "All peptides"]
))
fwrite(supp6_counts, file.path(supp6_source_dir, "SupplementaryFigure6C_measurable_protein_counts_source.csv"))
p6c <- ggplot(supp6_counts, aes(x = analysis, y = n_common_raw20, fill = analysis, color = analysis)) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35) +
  geom_jitter(width = 0.08, height = 0.05, size = 0.95, alpha = 0.72, shape = 21, stroke = 0.15, fill = "white") +
  scale_fill_manual(values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"]))) +
  scale_color_manual(values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"]))) +
  labs(x = NULL, y = "Detectable raw20 proteins per pair") +
  theme_main +
  theme(legend.position = "none")
save_plot_pair(p6c, "supp6c_measurable_protein_counts", supp6_plot_dir, width = 3.9, height = 4.0)
supp6_panel_manifest <- rbindlist(list(
  supp6_panel_manifest,
  panel_manifest("SupplementaryFigure6C", "supp6c_measurable_protein_counts", "Detectable raw20 protein counts per paired patient in PDC.", supp6_plot_dir)
))

supp6_score_sensitivity <- rbindlist(list(
  data.table(
    analysis = "Unique peptides",
    mean_delta = pdc$unique$score_tests$mean_delta,
    ci_low = pdc$unique$score_tests$ci_low,
    ci_high = pdc$unique$score_tests$ci_high,
    fdr = pdc$unique$score_tests$sign_flip_p,
    n_pairs = pdc$unique$score_tests$n_pairs,
    n_recurrence_up = pdc$unique$score_tests$n_recurrence_up
  ),
  data.table(
    analysis = "All peptides",
    mean_delta = pdc$all$score_tests$mean_delta,
    ci_low = pdc$all$score_tests$ci_low,
    ci_high = pdc$all$score_tests$ci_high,
    fdr = pdc$all$score_tests$sign_flip_p,
    n_pairs = pdc$all$score_tests$n_pairs,
    n_recurrence_up = pdc$all$score_tests$n_recurrence_up
  )
))
supp6_score_sensitivity[, analysis := factor(analysis, levels = rev(c("Unique peptides", "All peptides")))]
fwrite(supp6_score_sensitivity, file.path(supp6_source_dir, "SupplementaryFigure6D_peptide_score_sensitivity_source.csv"))
p6d <- ggplot(supp6_score_sensitivity, aes(x = mean_delta, y = analysis, color = analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = analysis), linewidth = 0.95) +
  geom_point(size = 3.1) +
  geom_text(
    aes(x = pmax(mean_delta, ci_high) + 0.06, label = sprintf("%d/%d up\nP %.2g", n_recurrence_up, n_pairs, fdr)),
    hjust = 0, size = 2.7, color = "black"
  ) +
  scale_color_manual(values = c("Unique peptides" = unname(npg["red"]), "All peptides" = unname(npg["teal"]))) +
  coord_cartesian(xlim = c(-0.02, max(supp6_score_sensitivity$ci_high) + 0.25), clip = "off") +
  labs(x = "Mean paired raw20 score delta", y = NULL, color = NULL) +
  theme_main +
  theme(legend.position = "none")
save_plot_pair(p6d, "supp6d_peptide_score_sensitivity", supp6_plot_dir, width = 4.8, height = 2.5)
supp6_panel_manifest <- rbindlist(list(
  supp6_panel_manifest,
  panel_manifest("SupplementaryFigure6D", "supp6d_peptide_score_sensitivity", "Score-scale sensitivity for unique versus all peptides.", supp6_plot_dir)
))

supp6_gsea_sensitivity <- rbindlist(list(
  data.table(dataset = "PDC unique", NES = pdc$unique$score_tests$gsea_NES, p = pdc$unique$score_tests$gsea_p, detected = pdc$unique$coverage$n_detected, pairs = pdc$unique$score_tests$n_pairs),
  data.table(dataset = "PDC all", NES = pdc$all$score_tests$gsea_NES, p = pdc$all$score_tests$gsea_p, detected = pdc$all$coverage$n_detected, pairs = pdc$all$score_tests$n_pairs),
  data.table(dataset = "MSV canonical", NES = msv[analysis == "canonical_regions", gsea_NES], p = msv[analysis == "canonical_regions", gsea_p], detected = msv[analysis == "canonical_regions", n_detected], pairs = msv[analysis == "canonical_regions", n_pairs]),
  data.table(dataset = "MSV averaged", NES = msv[analysis == "averaged_repeated_regions", gsea_NES], p = msv[analysis == "averaged_repeated_regions", gsea_p], detected = msv[analysis == "averaged_repeated_regions", n_detected], pairs = msv[analysis == "averaged_repeated_regions", n_pairs]),
  data.table(dataset = "Zenodo", NES = zenodo$gsea_NES, p = zenodo$gsea_p, detected = zenodo$n_detected, pairs = zenodo$n_pairs)
))
supp6_gsea_sensitivity[, dataset := factor(dataset, levels = rev(dataset))]
fwrite(supp6_gsea_sensitivity, file.path(supp6_source_dir, "SupplementaryFigure6E_gsea_sensitivity_source.csv"))
p6e <- ggplot(supp6_gsea_sensitivity, aes(x = NES, y = dataset, color = dataset)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = npg["grey"], linewidth = 0.4) +
  geom_point(size = 3.1) +
  geom_text(
    aes(label = sprintf("nominal P %.2g  |  %d/20  |  n=%d", p, detected, pairs)),
    nudge_x = 0.09, hjust = 0, size = 2.7, color = "black"
  ) +
  scale_color_manual(values = c(
    "PDC unique" = unname(npg["red"]),
    "PDC all" = unname(npg["teal"]),
    "MSV canonical" = unname(npg["blue"]),
    "MSV averaged" = unname(npg["orange"]),
    "Zenodo" = unname(npg["navy"])
  )) +
  coord_cartesian(xlim = c(min(supp6_gsea_sensitivity$NES) - 0.12, max(supp6_gsea_sensitivity$NES) + 1.0), clip = "off") +
  labs(x = "GSEA NES", y = NULL, color = NULL) +
  theme_main +
  theme(legend.position = "none")
save_plot_pair(p6e, "supp6e_gsea_sensitivity", supp6_plot_dir, width = 5.5, height = 3.1)
supp6_panel_manifest <- rbindlist(list(
  supp6_panel_manifest,
  panel_manifest("SupplementaryFigure6E", "supp6e_gsea_sensitivity", "GSEA-scale sensitivity across proteome datasets and peptide handling.", supp6_plot_dir)
))

supp6_inputs <- data.table(
  input = c(
    "PDC000514 raw TMT matrix",
    "MSV000087947 processed paired proteome matrix",
    "Zenodo7646550 processed paired proteome matrix"
  ),
  location = c(pdc_source_dir, msv_source_dir, zenodo_source_root)
)
writeLines(c(
  "# Supplementary Figure 6 summary",
  "",
  sprintf("- PDC unique/all both measure %d/20 raw20 proteins; main inference remains unique-peptide PDC.", pdc$unique$coverage$n_detected),
  sprintf("- MSV sensitivity is limited to %d/20 proteins and Zenodo to %d/20 proteins.", msv[analysis == "canonical_regions", n_detected], zenodo$n_detected),
  sprintf("- GSEA sensitivity is shown on its own NES axis and is not mixed with score delta.")
), file.path(supp6_write_dir, "SupplementaryFigure6_summary.md"), useBytes = TRUE)
fwrite(supp6_panel_manifest, file.path(supp6_write_dir, "SupplementaryFigure6_panel_manifest.csv"))
fwrite(supp6_inputs, file.path(supp6_write_dir, "SupplementaryFigure6_inputs.csv"))
