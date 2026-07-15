#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(readxl)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(ggplot2)
  library(methods)
  library(patchwork)
  library(cowplot)
  library(gridExtra)
  library(grid)
})

set.seed(20260713)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else dirname(dirname(script_path))

write_root <- file.path(root_dir, "write", "41_mg_inflammatory_sci_rebuild")
figure_root <- file.path(root_dir, "figures", "41_mg_inflammatory_sci_rebuild")
figure5_write_dir <- file.path(write_root, "Figure5")
s7_write_dir <- file.path(write_root, "SupplementaryFigure7")
figure5_dir <- file.path(figure_root, "Figure5", "panel_library")
s7_dir <- file.path(figure_root, "SupplementaryFigure7", "panel_library")
figure5_source_dir <- file.path(figure5_write_dir, "source_data")
s7_source_dir <- file.path(s7_write_dir, "source_data")

dirs_to_make <- c(
  write_root, figure_root,
  figure5_write_dir, s7_write_dir,
  figure5_dir, s7_dir,
  figure5_source_dir, s7_source_dir
)
invisible(lapply(dirs_to_make, dir.create, recursive = TRUE, showWarnings = FALSE))

external_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root_dir, "data"))
gse121810_source <- file.path(external_root, "write", "28_supportive_dataset_search", "source_metadata", "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx")
gse154795_rds <- file.path(external_root, "write", "26_pure_bioinformatics_dataset_rescue", "source_metadata", "GSE154795_GBM.AllCell.Integrated.Scaled.ClusterRes.0.1.rds.gz")
gse154795_soft <- file.path(external_root, "write", "26_pure_bioinformatics_dataset_rescue", "source_metadata", "GSE154795_family.soft.gz")

step38_dir <- file.path(root_dir, "write", "38_independent_cohort_mg_inflammatory_recalculation")
step38_gsea_file <- file.path(step38_dir, "independent_fixed_program_targeted_gsea.csv")
step38_gene_dir_file <- file.path(step38_dir, "independent_fixed_program_gene_direction.csv")

pal <- c(
  primary = "#4C78A8",
  recurrent = "#E45756",
  neutral = "#777777",
  light = "#D9D9D9",
  grid = "#E8E8E8",
  text = "#222222",
  teal = "#00A087",
  blue = "#3C5488",
  mint = "#91D1C2",
  grayblue = "#8491B4"
)

raw20 <- c(
  "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
  "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
  "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)

myeloid_genes <- c(
  "TREM2", "TMEM119", "P2RY12", "CX3CR1", "AIF1", "TYROBP", "FCER1G", "CSF1R",
  "CD68", "C1QA", "C1QB", "C1QC", "APOC1", "GPNMB", "LGALS3", "FCGR1A"
)

theme_panel <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(color = pal["text"], face = "bold", size = base_size + 0.7, hjust = 0),
      plot.subtitle = element_text(color = pal["neutral"], size = base_size - 0.5, hjust = 0),
      axis.title = element_text(color = pal["text"]),
      axis.text = element_text(color = pal["text"]),
      axis.line = element_line(linewidth = 0.35, color = pal["text"]),
      axis.ticks = element_line(linewidth = 0.35, color = pal["text"]),
      legend.title = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = pal["text"], face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

save_panel <- function(plot_obj, stem, out_dir, width_mm, height_mm) {
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::pdf(
    file = pdf_file,
    width = width_in,
    height = height_in,
    onefile = FALSE,
    paper = "special",
    bg = "white",
    useDingbats = FALSE,
    compress = TRUE,
    timestamp = FALSE,
    producer = "R grDevices",
    author = "Codex"
  )
  print(plot_obj)
  dev.off()
  ggsave(png_file, plot_obj, width = width_in, height = height_in, units = "in", bg = "white", dpi = 600)
}

write_csv <- function(x, path_stem, out_dir) {
  fwrite(as.data.table(x), file.path(out_dir, paste0(path_stem, ".csv")))
}

write_csv_gz <- function(x, path_stem, out_dir) {
  fwrite(as.data.table(x), file.path(out_dir, paste0(path_stem, ".csv.gz")))
}

write_markdown_table <- function(df, file_path, title_text, intro_lines = character()) {
  display_df <- as.data.frame(df, stringsAsFactors = FALSE)
  header <- paste(names(display_df), collapse = " | ")
  separator <- paste(rep("---", ncol(display_df)), collapse = " | ")
  body <- apply(display_df, 1, function(row) paste(as.character(row), collapse = " | "))
  lines <- c(
    paste0("# ", title_text),
    "",
    intro_lines,
    "",
    paste0("| ", header, " |"),
    paste0("| ", separator, " |"),
    paste0("| ", body, " |")
  )
  writeLines(lines, file_path, useBytes = TRUE)
}

read_rds_auto <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e1) {
      tryCatch(
        readRDS(gzfile(path, "rb")),
        error = function(e2) {
          tryCatch(
            readRDS(gzcon(gzfile(path, "rb"))),
            error = function(e3) stop("Cannot read RDS file: ", path, "\n", conditionMessage(e3))
          )
        }
      )
    }
  )
}

safe_zscore_rows <- function(x) {
  m <- rowMeans(x, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- NA_real_
  z <- sweep(x, 1, m, "-")
  z <- sweep(z, 1, s, "/")
  z[!is.finite(z)] <- 0
  z
}

build_p_label <- function(p) {
  if (!is.finite(p) || is.na(p)) return("NA")
  if (p < 1e-3) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

extract_model_row <- function(model, coef_name, cohort, metric, metric_label, cohort_label, comparator, n_positive, n_negative, extra = list()) {
  coef_tab <- summary(model)$coefficients
  ci <- suppressMessages(confint(model, parm = coef_name, level = 0.95))
  row <- data.frame(
    cohort = cohort,
    cohort_label = cohort_label,
    metric = metric,
    metric_label = metric_label,
    comparator = comparator,
    n_positive = n_positive,
    n_negative = n_negative,
    effect = unname(coef(model)[coef_name]),
    conf_low = unname(ci[1]),
    conf_high = unname(ci[2]),
    p_value = unname(coef_tab[coef_name, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
  if (length(extra)) {
    for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  }
  row
}

compute_gsea_curve <- function(ranks, genes, cohort_label, contrast_label) {
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  gene_names <- names(ranks)
  present <- intersect(genes, gene_names)
  hits <- gene_names %in% present
  hit_weight <- abs(ranks)
  nr <- sum(hit_weight[hits])
  if (nr <= 0 || sum(hits) == 0) stop("No valid hits for GSEA curve: ", cohort_label)
  miss_penalty <- 1 / sum(!hits)
  increments <- ifelse(hits, hit_weight / nr, -miss_penalty)
  running <- cumsum(increments)
  data.frame(
    cohort = cohort_label,
    contrast = contrast_label,
    rank = seq_along(ranks),
    gene = gene_names,
    rank_stat = unname(ranks),
    running_enrichment = running,
    hit = hits,
    hit_gene = ifelse(hits, gene_names, NA_character_),
    stringsAsFactors = FALSE
  )
}

build_gsea_panel <- function(curve_df, stats_row, stem_label, entry_label) {
  running_df <- as.data.table(copy(curve_df))
  hit_df <- running_df[hit == TRUE, .(rank)]
  rank_df <- data.frame(rank = running_df$rank, rank_stat = running_df$rank_stat)
  stats_lines <- c(
    stem_label,
    entry_label,
    sprintf("%d/20 genes measured", stats_row$n_present),
    sprintf("%s", stats_row$comparison_label),
    sprintf("NES %.2f   P %s", stats_row$NES, build_p_label(stats_row$fgsea_p)),
    sprintf("camera %s   P %s", stats_row$camera_direction, build_p_label(stats_row$camera_p)),
    sprintf("n = %d vs %d", stats_row$n_positive, stats_row$n_negative)
  )

  top_plot <- ggplot(running_df, aes(x = rank, y = running_enrichment)) +
    geom_hline(yintercept = 0, color = pal["light"], linewidth = 0.35) +
    geom_line(color = pal["recurrent"], linewidth = 0.75) +
    geom_area(fill = pal["recurrent"], alpha = 0.12) +
    labs(x = NULL, y = "Running enrichment") +
    theme_panel(7.2) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  hit_plot <- ggplot(hit_df, aes(x = rank, y = 1)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), color = pal["text"], linewidth = 0.22) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_void(base_family = "sans")

  rank_plot <- ggplot(rank_df, aes(x = rank, y = rank_stat)) +
    geom_hline(yintercept = 0, color = pal["light"], linewidth = 0.3) +
    geom_area(data = rank_df[rank_df$rank_stat >= 0, ], fill = pal["recurrent"], alpha = 0.85) +
    geom_area(data = rank_df[rank_df$rank_stat < 0, ], fill = pal["primary"], alpha = 0.85) +
    labs(x = "Ranked genes", y = NULL) +
    theme_panel(7.2) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  left <- top_plot / hit_plot / rank_plot + plot_layout(heights = c(5.0, 0.65, 1.4))

  text_df <- data.frame(
    x = 0,
    y = rev(seq_along(stats_lines)),
    label = stats_lines,
    face = c("bold", rep("plain", length(stats_lines) - 1)),
    stringsAsFactors = FALSE
  )
  right <- ggplot(text_df, aes(x = x, y = y, label = label)) +
    geom_text(data = text_df[1, , drop = FALSE], aes(fontface = face), hjust = 0, vjust = 1, size = 3.1) +
    geom_text(data = text_df[-1, , drop = FALSE], aes(fontface = face), hjust = 0, vjust = 1, size = 2.75) +
    xlim(0, 1) +
    ylim(0.4, length(stats_lines) + 0.6) +
    theme_void(base_family = "sans")

  (left | right) + plot_layout(widths = c(4.0, 2.8))
}

build_swimlane_plot <- function(timeline_df, title_text = NULL, subtitle_text = NULL) {
  ggplot(timeline_df, aes(x = step, y = lane)) +
    geom_segment(aes(x = 1, xend = 3, y = lane, yend = lane), color = pal["light"], linewidth = 1.7, lineend = "round") +
    geom_segment(
      data = timeline_df[timeline_df$exposed_preop == "Yes", ],
      aes(x = 1, xend = 2, y = lane, yend = lane),
      color = pal["recurrent"],
      linewidth = 1.7,
      lineend = "round"
    ) +
    geom_segment(
      data = timeline_df[timeline_df$starts_postop == "Yes", ],
      aes(x = 2, xend = 3, y = lane, yend = lane),
      color = pal["primary"],
      linewidth = 1.7,
      lineend = "round"
    ) +
    geom_point(aes(fill = event_group, shape = event), size = 2.7, color = pal["text"], stroke = 0.3) +
    geom_text(aes(x = step + label_nudge, y = lane + label_offset, label = label, hjust = label_hjust), size = 2.55, lineheight = 0.9) +
    scale_shape_manual(values = c(Exposure = 21, Surgery = 24)) +
    scale_fill_manual(values = c("Exposure" = pal["recurrent"], "Surgery" = "white")) +
    scale_y_continuous(NULL, breaks = timeline_df$lane, labels = timeline_df$lane_label) +
    scale_x_continuous(
      NULL,
      breaks = 1:3,
      labels = c("Pre-op", "Surgery", "Post-op"),
      expand = expansion(mult = c(0.04, 0.10))
    ) +
    labs(title = title_text, subtitle = subtitle_text) +
    theme_panel(7.3) +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
      plot.margin = margin(8, 12, 6, 6)
    )
}

build_forest_plot <- function(df, x_label, point_colors, facet_var = NULL, title_text = NULL, subtitle_text = NULL) {
  plot_df <- copy(df)
  plot_df$row_id <- seq_len(nrow(plot_df))
  plot_df$label_y <- plot_df$row_label
  x_rng <- range(c(plot_df$conf_low, plot_df$conf_high), finite = TRUE)
  span <- diff(x_rng)
  if (!is.finite(span) || span <= 0) span <- 0.2
  show_n <- any(nchar(plot_df$n_label) > 0)
  text_x1 <- x_rng[2] + span * 0.16
  text_x2 <- if (show_n) x_rng[2] + span * 0.46 else x_rng[2] + span * 0.24
  xmax <- x_rng[2] + span * 0.95
  p <- ggplot(plot_df, aes(x = effect, y = factor(row_label, levels = rev(unique(row_label))))) +
    geom_vline(xintercept = 0, color = pal["light"], linewidth = 0.35) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high, color = metric_label), height = 0.18, linewidth = 0.55) +
    geom_point(aes(color = metric_label, shape = metric_label), size = 2.1, stroke = 0.25) +
    geom_text(aes(x = text_x2, label = fdr_label), size = 2.5, hjust = 0) +
    scale_color_manual(values = setNames(unname(point_colors), names(point_colors))) +
    scale_shape_manual(values = c("Raw score" = 16, "Myeloid-adjusted" = 17)) +
    scale_x_continuous(x_label, expand = expansion(mult = c(0.08, 0.18)), limits = c(x_rng[1] - span * 0.08, xmax)) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, title = title_text, subtitle = subtitle_text) +
    theme_panel(7.2) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
      plot.margin = margin(12, 42, 6, 6)
    )
  if (show_n) {
    p <- p +
      geom_text(aes(x = text_x1, label = n_label), size = 2.5, hjust = 0) +
      annotate("text", x = text_x1, y = length(unique(plot_df$row_label)) + 0.24, label = "n", hjust = 0, size = 2.5, fontface = "bold")
  }
  p <- p + annotate("text", x = text_x2, y = length(unique(plot_df$row_label)) + 0.24, label = "FDR", hjust = 0, size = 2.5, fontface = "bold")
  if (!is.null(facet_var)) {
    p <- p + facet_grid(as.formula(paste(facet_var, "~ .")), scales = "free_y", space = "free_y")
  }
  p
}

build_nes_plot <- function(df, x_label, title_text = NULL, subtitle_text = NULL) {
  plot_df <- copy(df)
  plot_df$row_label <- factor(plot_df$row_label, levels = rev(plot_df$row_label))
  x_rng <- range(plot_df$NES, finite = TRUE)
  span <- diff(x_rng)
  if (!is.finite(span) || span <= 0) span <- 0.4
  text_x <- x_rng[2] + span * 0.35
  xmax <- x_rng[2] + span * 0.95
  ggplot(plot_df, aes(x = NES, y = row_label)) +
    geom_vline(xintercept = 0, color = pal["light"], linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = NES, y = row_label, yend = row_label), color = pal["grayblue"], linewidth = 0.55) +
    geom_point(size = 2.2, color = pal["recurrent"]) +
    geom_text(aes(x = text_x, label = p_label), hjust = 0, size = 2.55) +
    annotate("text", x = text_x, y = length(unique(plot_df$row_label)) + 0.36, label = "P", hjust = 0, size = 2.45, fontface = "bold") +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.82))) +
    scale_x_continuous(x_label, expand = expansion(mult = c(0.08, 0.18)), limits = c(min(0, x_rng[1] - span * 0.08), xmax)) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, title = title_text, subtitle = subtitle_text) +
    theme_panel(7.2) +
    theme(
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
      plot.margin = margin(12, 36, 6, 6)
    )
}

build_algorithm_table <- function(df) {
  display_df <- as.data.frame(df[, .(
    Cohort = cohort_label,
    Analysis = analysis_label,
    `n exp / ctrl` = n_label,
    fgsea = fgsea_label,
    camera = camera_label
  )])
  tbl <- tableGrob(
    display_df,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontsize = 7, fontfamily = "sans", col = pal["text"]),
        bg_params = list(fill = c(rep("white", nrow(display_df))), col = "grey85", lwd = 0.4)
      ),
      colhead = list(
        fg_params = list(fontsize = 7.2, fontface = "bold", fontfamily = "sans", col = pal["text"]),
        bg_params = list(fill = "#F5F5F5", col = "grey80", lwd = 0.5)
      )
    )
  )
  ggdraw() +
    draw_label("Algorithm summary kept off shared axes", x = 0, y = 1, hjust = 0, vjust = 1, fontfamily = "sans", fontface = "bold", size = 8.5) +
    draw_grob(tbl, x = 0, y = 0.03, width = 1, height = 0.86)
}

build_leading_edge_heatmap <- function(df) {
  plot_df <- copy(df)
  plot_df$cohort_label <- factor(
    plot_df$cohort_label,
    levels = c("GSE174554", "GSE274546", "GSE121810", "GSE154795"),
    labels = c(
      "GSE174554\nRecurrent logFC",
      "GSE274546\nRecurrent logFC",
      "GSE121810\nExposure logFC",
      "GSE154795\nExposure logFC"
    )
  )
  plot_df$gene <- factor(plot_df$gene, levels = rev(unique(plot_df$gene)))
  leading_df <- plot_df[leading_edge == TRUE, ]
  ggplot(plot_df, aes(x = cohort_label, y = gene)) +
    geom_tile(aes(fill = logFC), color = "white", linewidth = 0.32) +
    geom_point(
      data = leading_df,
      aes(shape = "Leading-edge member"),
      size = 1.7,
      stroke = 0.15,
      color = "black",
      inherit.aes = TRUE
    ) +
    scale_fill_gradient2(
      low = pal["primary"],
      mid = "white",
      high = pal["recurrent"],
      midpoint = 0,
      name = "logFC"
    ) +
    scale_shape_manual(values = c("Leading-edge member" = 16), name = NULL) +
    guides(
      fill = guide_colorbar(order = 1, title.position = "top"),
      shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.2))
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Leading-edge genes keep the same direction across recurrence and exposure cohorts",
      subtitle = "First two columns are paired recurrence logFC; last two columns are anti-PD-1 exposure logFC"
    ) +
    theme_panel(7) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 8.0, hjust = 0),
      plot.subtitle = element_text(size = 6.6, color = pal["neutral"], hjust = 0)
    )
}

analyze_gse121810 <- function(source_file) {
  x <- as.data.frame(read_excel(source_file), check.names = FALSE)
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x <- rowsum(as.matrix(x), group = rownames(x), reorder = FALSE)
  storage.mode(x) <- "numeric"
  group <- factor(ifelse(grepl("_A$", colnames(x)), "Neo", "Adj"), levels = c("Adj", "Neo"))
  design <- model.matrix(~ group)
  y <- DGEList(x)
  keep <- filterByExpr(y, design = design, min.count = 5, min.total.count = 15)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupNeo")
  de <- topTags(qlf, n = Inf, sort.by = "none")$table
  de$gene <- rownames(de)
  de$FDR <- p.adjust(de$PValue, method = "BH")

  ranks <- de$F
  names(ranks) <- de$gene
  ranks <- sort(ranks * sign(de$logFC), decreasing = TRUE)
  fg <- as.data.table(fgseaMultilevel(
    pathways = list(Miller_raw20 = raw20),
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 0,
    nproc = 1
  ))
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  fg[, cohort := "GSE121810"]
  fg[, contrast := "Neoadjuvant_minus_AdjuvantOnly"]

  v <- voom(y, design, plot = FALSE)
  raw20_indices <- ids2indices(list(Miller_raw20 = raw20), rownames(v))
  camera_dt <- as.data.table(camera(v, raw20_indices, design, contrast = "groupNeo"), keep.rownames = "pathway")
  if (!"FDR" %in% names(camera_dt)) camera_dt[, FDR := p.adjust(PValue, method = "BH")]
  camera_dt[, cohort := "GSE121810"]
  camera_dt[, contrast := "Neoadjuvant_minus_AdjuvantOnly"]

  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  z <- safe_zscore_rows(logcpm)
  raw20_present <- intersect(raw20, rownames(z))
  myeloid_present <- intersect(myeloid_genes, rownames(z))
  scores <- data.frame(
    sample = colnames(z),
    cohort = "GSE121810",
    group = group,
    raw20_score = colMeans(z[raw20_present, , drop = FALSE]),
    myeloid_abundance = colMeans(z[myeloid_present, , drop = FALSE]),
    stringsAsFactors = FALSE
  )

  raw_model <- lm(raw20_score ~ group, data = scores)
  adj_model <- lm(raw20_score ~ myeloid_abundance + group, data = scores)
  score_models <- rbindlist(list(
    extract_model_row(raw_model, "groupNeo", "GSE121810", "raw20_score", "Raw score", "GSE121810", "Neo_minus_Adj", sum(scores$group == "Neo"), sum(scores$group == "Adj")),
    extract_model_row(adj_model, "groupNeo", "GSE121810", "raw20_score_myeloid_adjusted", "Myeloid-adjusted", "GSE121810", "Neo_minus_Adj", sum(scores$group == "Neo"), sum(scores$group == "Adj"))
  ))
  score_models[, fdr := p.adjust(p_value, method = "BH")]
  score_models[, n_label := sprintf("%d / %d", n_positive, n_negative)]
  score_models[, fdr_label := vapply(fdr, build_p_label, character(1))]
  coverage <- data.frame(
    cohort = "GSE121810",
    signature = "Miller_raw20",
    n_defined = length(raw20),
    n_present = length(raw20_present),
    coverage = length(raw20_present) / length(raw20),
    present_genes = paste(raw20_present, collapse = ";"),
    missing_genes = paste(setdiff(raw20, raw20_present), collapse = ";"),
    stringsAsFactors = FALSE
  )
  list(
    de = as.data.table(de),
    fgsea = fg,
    camera = camera_dt,
    scores = as.data.table(scores),
    score_models = score_models,
    coverage = coverage,
    ranks = ranks
  )
}

build_gse154795_audit <- function(obj) {
  counts <- obj@assays$RNA@counts
  meta <- obj@meta.data
  stopifnot(ncol(counts) == nrow(meta))
  meta$cell_id <- rownames(meta)
  meta$cluster <- as.character(meta$seurat_clusters)
  lineages <- list(
    Myeloid = c("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "CSF1R", "CD68", "C1QA", "C1QB", "C1QC"),
    T_NK = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "CD37", "CD22", "BANK1"),
    Plasma = c("MZB1", "JCHAIN", "SDC1", "IGHG1", "IGHG3"),
    Neural_glial = c("GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "OLIG1", "OLIG2"),
    Erythroid = c("HBA1", "HBA2", "HBB", "ALAS2", "GYPA")
  )
  cluster_factor <- factor(meta$cluster, levels = sort(unique(meta$cluster)))
  cluster_design <- sparse.model.matrix(~ 0 + cluster_factor)
  colnames(cluster_design) <- levels(cluster_factor)
  cluster_counts <- counts %*% cluster_design
  cluster_cpm <- cpm(cluster_counts, log = TRUE, prior.count = 1)
  lineage_score <- matrix(
    NA_real_,
    nrow = ncol(cluster_cpm),
    ncol = length(lineages),
    dimnames = list(colnames(cluster_cpm), names(lineages))
  )
  for (name in names(lineages)) {
    present <- intersect(lineages[[name]], rownames(cluster_cpm))
    lineage_score[, name] <- colMeans(cluster_cpm[present, , drop = FALSE])
  }
  substantial <- as.integer(table(cluster_factor)[rownames(lineage_score)]) >= 100
  lineage_z <- scale(lineage_score[substantial, , drop = FALSE])
  best <- colnames(lineage_score)[max.col(lineage_score, ties.method = "first")]
  best[substantial] <- colnames(lineage_z)[max.col(lineage_z, ties.method = "first")]
  audit <- data.frame(
    cluster = rownames(lineage_score),
    n_cells = as.integer(table(cluster_factor)[rownames(lineage_score)]),
    best_lineage = best,
    lineage_score,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  audit$is_auto_myeloid <- audit$best_lineage == "Myeloid" & audit$n_cells >= 100
  list(counts = counts, meta = meta, audit = audit)
}

get_gse154795_cluster_selection <- function(audit, entry) {
  if (entry == "reviewed_scgbm") return(c("0", "1", "2", "5", "8"))
  if (entry == "strict") return(audit$cluster[audit$is_auto_myeloid])
  if (entry == "erythroid_ambient_inclusive") {
    extra <- audit$cluster[audit$best_lineage == "Erythroid" & audit$Myeloid >= 9.5 & audit$n_cells >= 100]
    return(unique(c(audit$cluster[audit$is_auto_myeloid], extra)))
  }
  stop("Unknown GSE154795 entry: ", entry)
}

extract_gse154795_tech <- function(soft_file) {
  soft <- readLines(gzfile(soft_file), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", soft)
  ends <- c(starts[-1] - 1, length(soft))
  rows <- lapply(seq_along(starts), function(i) {
    block <- soft[starts[i]:ends[i]]
    get_one <- function(pattern) {
      hit <- grep(pattern, block, value = TRUE)
      if (length(hit) == 0) return(NA_character_)
      sub(pattern, "", hit[1])
    }
    gsm <- get_one("^\\^SAMPLE = ")
    title <- get_one("^!Sample_title = ")
    data.frame(
      ID = title,
      GSM = gsm,
      sequencing_method = get_one("^!Sample_characteristics_ch1 = sequencing_method: "),
      gel_bead_version = get_one("^!Sample_characteristics_ch1 = gel_bead_version: "),
      acquisition_wave = ifelse(grepl("^GSM467", gsm), "early", "late"),
      stringsAsFactors = FALSE
    )
  })
  unique(rbindlist(rows, fill = TRUE))
}

analyze_gse154795_variant <- function(prepped, tech, entry, idh_entry) {
  meta <- prepped$meta
  counts <- prepped$counts
  audit <- prepped$audit
  selected_clusters <- get_gse154795_cluster_selection(audit, entry)
  keep <- meta$cluster %in% selected_clusters
  myeloid_meta <- meta[keep, , drop = FALSE]
  myeloid_counts <- counts[, keep, drop = FALSE]

  sample_factor <- factor(myeloid_meta$ID, levels = unique(myeloid_meta$ID))
  sample_design <- sparse.model.matrix(~ 0 + sample_factor)
  colnames(sample_design) <- levels(sample_factor)
  pseudobulk <- myeloid_counts %*% sample_design
  sample_cells <- as.integer(table(sample_factor)[colnames(pseudobulk)])

  sample_meta <- unique(myeloid_meta[, c("ID", "condition", "IDH_status", "MGMT")])
  rownames(sample_meta) <- sample_meta$ID
  sample_meta <- sample_meta[colnames(pseudobulk), , drop = FALSE]
  sample_meta$n_myeloid_cells <- sample_cells
  sample_meta <- merge(sample_meta, tech, by = "ID", all.x = TRUE, sort = FALSE)
  rownames(sample_meta) <- sample_meta$ID
  sample_meta <- sample_meta[colnames(pseudobulk), , drop = FALSE]

  formal <- sample_meta$condition %in% c("GBM.rec", "GBM.PD1") & sample_meta$n_myeloid_cells >= 20
  if (idh_entry == "strict_idhwt") formal <- formal & sample_meta$IDH_status == "WT"
  formal_counts <- pseudobulk[, formal, drop = FALSE]
  formal_meta <- sample_meta[formal, , drop = FALSE]
  formal_meta$condition <- factor(formal_meta$condition, levels = c("GBM.rec", "GBM.PD1"))

  y_all <- calcNormFactors(DGEList(formal_counts))
  y <- y_all
  keep_gene <- filterByExpr(y, group = formal_meta$condition, min.count = 5, min.total.count = 10)
  y <- calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
  formal_meta$MGMT <- factor(formal_meta$MGMT)
  formal_meta$sequencing_method <- factor(formal_meta$sequencing_method)
  formal_meta$gel_bead_version <- factor(formal_meta$gel_bead_version)
  formal_meta$acquisition_wave <- factor(formal_meta$acquisition_wave)
  formal_meta$IDH_status <- factor(formal_meta$IDH_status)
  if (idh_entry == "strict_idhwt") {
    design <- model.matrix(~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition, formal_meta)
  } else {
    design <- model.matrix(~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + IDH_status + condition, formal_meta)
  }
  if (qr(design)$rank != ncol(design)) stop("GSE154795 design is not full rank for ", entry, " / ", idh_entry)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  coef_index <- grep("^condition", colnames(design))
  if (length(coef_index) != 1) stop("Cannot find unique treatment coefficient for ", entry, " / ", idh_entry)
  qlf <- glmQLFTest(fit, coef = coef_index)
  de <- topTags(qlf, n = Inf, sort.by = "none")$table
  de$gene <- rownames(de)
  de$FDR <- p.adjust(de$PValue, method = "BH")

  ranks <- de$F
  names(ranks) <- de$gene
  ranks <- sort(ranks * sign(de$logFC), decreasing = TRUE)
  fg <- as.data.table(fgseaMultilevel(
    pathways = list(Miller_raw20 = raw20),
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 0,
    nproc = 1
  ))
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  fg[, cohort := "GSE154795"]
  fg[, entry := entry]
  fg[, idh_entry := idh_entry]
  fg[, contrast := "PD1_minus_UntreatedRecurrent"]

  v <- voom(y, design, plot = FALSE)
  raw20_indices <- ids2indices(list(Miller_raw20 = raw20), rownames(v))
  camera_dt <- as.data.table(camera(v, raw20_indices, design, contrast = colnames(design)[coef_index]), keep.rownames = "pathway")
  if (!"FDR" %in% names(camera_dt)) camera_dt[, FDR := p.adjust(PValue, method = "BH")]
  camera_dt[, cohort := "GSE154795"]
  camera_dt[, entry := entry]
  camera_dt[, idh_entry := idh_entry]
  camera_dt[, contrast := "PD1_minus_UntreatedRecurrent"]

  logcpm <- cpm(y_all, log = TRUE, prior.count = 1)
  z <- safe_zscore_rows(logcpm)
  raw20_present <- intersect(raw20, rownames(z))
  myeloid_present <- intersect(myeloid_genes, rownames(z))
  scores <- data.frame(
    sample = colnames(z),
    cohort = "GSE154795",
    entry = entry,
    idh_entry = idh_entry,
    condition = formal_meta$condition,
    n_myeloid_cells = formal_meta$n_myeloid_cells,
    raw20_score = colMeans(z[raw20_present, , drop = FALSE]),
    myeloid_abundance = colMeans(z[myeloid_present, , drop = FALSE]),
    MGMT = formal_meta$MGMT,
    sequencing_method = formal_meta$sequencing_method,
    gel_bead_version = formal_meta$gel_bead_version,
    acquisition_wave = formal_meta$acquisition_wave,
    IDH_status = formal_meta$IDH_status,
    stringsAsFactors = FALSE
  )
  raw_formula <- if (idh_entry == "all_idh") {
    raw20_score ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + IDH_status + condition
  } else {
    raw20_score ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition
  }
  adj_formula <- if (idh_entry == "all_idh") {
    raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version + acquisition_wave + IDH_status + condition
  } else {
    raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition
  }
  raw_model <- lm(raw_formula, data = scores)
  adj_model <- lm(adj_formula, data = scores)
  coef_name <- grep("^condition", names(coef(raw_model)), value = TRUE)

  score_models <- rbindlist(list(
    extract_model_row(raw_model, coef_name, "GSE154795", "raw20_score", "Raw score", "GSE154795", "PD1_minus_recurrent", sum(scores$condition == "GBM.PD1"), sum(scores$condition == "GBM.rec"), extra = list(entry = entry, idh_entry = idh_entry)),
    extract_model_row(adj_model, coef_name, "GSE154795", "raw20_score_myeloid_adjusted", "Myeloid-adjusted", "GSE154795", "PD1_minus_recurrent", sum(scores$condition == "GBM.PD1"), sum(scores$condition == "GBM.rec"), extra = list(entry = entry, idh_entry = idh_entry))
  ))
  score_models[, fdr := p.adjust(p_value, method = "BH")]
  score_models[, n_label := sprintf("%d / %d", n_positive, n_negative)]
  score_models[, fdr_label := vapply(fdr, build_p_label, character(1))]

  coverage <- data.frame(
    cohort = "GSE154795",
    entry = entry,
    idh_entry = idh_entry,
    signature = "Miller_raw20",
    n_defined = length(raw20),
    n_present = length(raw20_present),
    coverage = length(raw20_present) / length(raw20),
    present_genes = paste(raw20_present, collapse = ";"),
    missing_genes = paste(setdiff(raw20, raw20_present), collapse = ";"),
    stringsAsFactors = FALSE
  )
  list(
    audit = as.data.table(audit),
    sample_meta = as.data.table(sample_meta),
    formal_meta = as.data.table(formal_meta),
    de = as.data.table(de),
    fgsea = fg,
    camera = camera_dt,
    scores = as.data.table(scores),
    score_models = score_models,
    coverage = coverage,
    ranks = ranks
  )
}

message("Re-running Figure5 treatment datasets with Miller raw20...")
gse121810 <- analyze_gse121810(gse121810_source)
gse154795_obj <- read_rds_auto(gse154795_rds)
gse154795_prep <- build_gse154795_audit(gse154795_obj)
rm(gse154795_obj)
gc()
gse154795_tech <- extract_gse154795_tech(gse154795_soft)

gse154795_main <- analyze_gse154795_variant(gse154795_prep, gse154795_tech, entry = "strict", idh_entry = "strict_idhwt")
gse154795_reviewed <- analyze_gse154795_variant(gse154795_prep, gse154795_tech, entry = "reviewed_scgbm", idh_entry = "strict_idhwt")
gse154795_inclusive <- analyze_gse154795_variant(gse154795_prep, gse154795_tech, entry = "erythroid_ambient_inclusive", idh_entry = "strict_idhwt")
gse154795_all_idh <- analyze_gse154795_variant(gse154795_prep, gse154795_tech, entry = "strict", idh_entry = "all_idh")

step38_gsea <- fread(step38_gsea_file)
step38_gene_dir <- fread(step38_gene_dir_file)
step38_gsea <- step38_gsea[pathway == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20]
step38_gene_dir <- step38_gene_dir[signature == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20]

gse121810_leading <- unlist(strsplit(gse121810$fgsea$leadingEdge[1], ";"))
gse154795_leading <- unlist(strsplit(gse154795_main$fgsea$leadingEdge[1], ";"))
gse174554_leading <- unlist(strsplit(step38_gsea[dataset == "GSE174554"]$leadingEdge[1], ";"))
gse274546_leading <- unlist(strsplit(step38_gsea[dataset == "GSE274546"]$leadingEdge[1], ";"))

timeline_df <- rbindlist(list(
  data.frame(
    cohort = "GSE121810",
    lane = 4,
    lane_label = sprintf("GSE121810 neo-adjuvant (n=%d)", gse121810$score_models[1, n_positive]),
    step = c(1, 2),
    event = c("Exposure", "Surgery"),
    event_group = c("Exposure", "Surgery"),
    label = c("Pembrolizumab", "Resection"),
    label_offset = c(0.22, -0.22),
    exposed_preop = "Yes",
    starts_postop = "No",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE121810",
    lane = 3,
    lane_label = sprintf("GSE121810 adjuvant-only (n=%d)", gse121810$score_models[1, n_negative]),
    step = c(2, 3),
    event = c("Surgery", "Exposure"),
    event_group = c("Surgery", "Exposure"),
    label = c("Resection", "Starts after surgery"),
    label_offset = c(-0.22, 0.22),
    exposed_preop = "No",
    starts_postop = "Yes",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE154795",
    lane = 2,
    lane_label = sprintf("GSE154795 PD-1 exposed (n=%d)", gse154795_main$score_models[1, n_positive]),
    step = c(1, 2),
    event = c("Exposure", "Surgery"),
    event_group = c("Exposure", "Surgery"),
    label = c("Anti-PD-1", "Recurrent resection"),
    label_offset = c(0.22, -0.22),
    exposed_preop = "Yes",
    starts_postop = "No",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE154795",
    lane = 1,
    lane_label = sprintf("GSE154795 untreated recurrent (n=%d)", gse154795_main$score_models[1, n_negative]),
    step = 2,
    event = "Surgery",
    event_group = "Surgery",
    label = "Recurrent resection",
    label_offset = -0.22,
    exposed_preop = "No",
    starts_postop = "No",
    stringsAsFactors = FALSE
  )
))
timeline_df[, label_hjust := fifelse(step == 1, 0, fifelse(step == 3, 1, 0.5))]
timeline_df[, label_nudge := fifelse(step == 1, 0.06, fifelse(step == 3, -0.03, 0))]
write_csv(timeline_df, "figure5_treatment_timeline", figure5_source_dir)
save_panel(
  build_swimlane_plot(
    timeline_df,
    title_text = "Anti-PD-1 exposure timing differs across the two recurrent GBM cohorts",
    subtitle_text = "GSE121810 contrasts pre-op exposure with post-op-only start; GSE154795 compares exposed with untreated recurrent tumors"
  ),
  "fig5a_treatment_timeline",
  figure5_dir,
  183,
  58
)

gse121810_curve <- compute_gsea_curve(gse121810$ranks, raw20, "GSE121810", "Neo minus adjuvant-only")
gse154795_curve <- compute_gsea_curve(gse154795_main$ranks, raw20, "GSE154795", "PD1 minus untreated recurrent")
write_csv_gz(gse121810_curve, "figure5b_gse121810_gsea_curve", figure5_source_dir)
write_csv_gz(gse154795_curve, "figure5c_gse154795_gsea_curve", figure5_source_dir)

gsea_stats_df <- rbindlist(list(
  data.frame(
    cohort = "GSE121810",
    comparison_label = "Pre-op exposed\nvs post-op only",
    NES = gse121810$fgsea$NES[1],
    fgsea_p = gse121810$fgsea$pval[1],
    camera_direction = gse121810$camera$Direction[1],
    camera_p = gse121810$camera$PValue[1],
    n_positive = gse121810$score_models[1, n_positive],
    n_negative = gse121810$score_models[1, n_negative],
    n_present = gse121810$coverage$n_present[1],
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE154795",
    comparison_label = "PD-1 exposed\nvs untreated recurrent",
    NES = gse154795_main$fgsea$NES[1],
    fgsea_p = gse154795_main$fgsea$pval[1],
    camera_direction = gse154795_main$camera$Direction[1],
    camera_p = gse154795_main$camera$PValue[1],
    n_positive = gse154795_main$score_models[1, n_positive],
    n_negative = gse154795_main$score_models[1, n_negative],
    n_present = gse154795_main$coverage$n_present[1],
    stringsAsFactors = FALSE
  )
))
write_csv(gsea_stats_df, "figure5_gsea_side_stats", figure5_source_dir)

save_panel(
  build_gsea_panel(gse121810_curve, gsea_stats_df[cohort == "GSE121810", ], "GSE121810", "Randomized timing cohort"),
  "fig5b_gse121810_raw20_gsea",
  figure5_dir,
  135,
  82
)
save_panel(
  build_gsea_panel(gse154795_curve, gsea_stats_df[cohort == "GSE154795", ], "GSE154795", "Strict IDH-wt auto-myeloid entry"),
  "fig5c_gse154795_raw20_gsea",
  figure5_dir,
  135,
  82
)

main_effect_df <- rbindlist(list(
  gse121810$score_models[, .(cohort, cohort_label, metric, metric_label, effect, conf_low, conf_high, fdr, n_positive, n_negative)],
  gse154795_main$score_models[, .(cohort, cohort_label, metric, metric_label, effect, conf_low, conf_high, fdr, n_positive, n_negative)]
))
main_effect_df[, row_label := c(
  "GSE121810 raw score", "GSE121810 myeloid-adjusted",
  "GSE154795 raw score", "GSE154795 myeloid-adjusted"
)]
main_effect_df[, n_label := sprintf("%d / %d", n_positive, n_negative)]
main_effect_df[, fdr_label := vapply(fdr, build_p_label, character(1))]
write_csv(main_effect_df, "figure5d_patient_score_models", figure5_source_dir)

main_score_samples <- rbindlist(list(
  gse121810$scores[, .(sample, cohort, group, raw20_score, myeloid_abundance)],
  gse154795_main$scores[, .(sample, cohort, condition, raw20_score, myeloid_abundance, entry, idh_entry)]
), fill = TRUE)
write_csv(main_score_samples, "figure5d_patient_score_samples", figure5_source_dir)

save_panel(
  build_forest_plot(
    main_effect_df,
    x_label = "Treatment-associated score shift",
    point_colors = c("Raw score" = pal["primary"], "Myeloid-adjusted" = pal["recurrent"]),
    title_text = "Patient-level raw and myeloid-adjusted effects",
    subtitle_text = "Linear models with 95% CI; GSE154795 raw score keeps technical covariates"
  ),
  "fig5d_patient_score_forest",
  figure5_dir,
  140,
  88
)

recurrence_gene_df <- step38_gene_dir[, .(dataset, gene, logFC)]
gse121810_gene_df <- gse121810$de[, .(dataset = "GSE121810", gene, logFC)]
gse154795_gene_df <- gse154795_main$de[, .(dataset = "GSE154795", gene, logFC)]
matrix_df <- rbindlist(list(
  recurrence_gene_df[dataset %in% c("GSE174554", "GSE274546")],
  gse121810_gene_df,
  gse154795_gene_df
), fill = TRUE)
matrix_df <- matrix_df[gene %in% raw20]
matrix_df[, leading_edge := FALSE]
matrix_df[dataset == "GSE174554" & gene %in% gse174554_leading, leading_edge := TRUE]
matrix_df[dataset == "GSE274546" & gene %in% gse274546_leading, leading_edge := TRUE]
matrix_df[dataset == "GSE121810" & gene %in% gse121810_leading, leading_edge := TRUE]
matrix_df[dataset == "GSE154795" & gene %in% gse154795_leading, leading_edge := TRUE]
gene_order <- matrix_df[, .(leading_count = sum(leading_edge), mean_logFC = mean(logFC, na.rm = TRUE)), by = gene][order(-leading_count, -mean_logFC)]$gene
matrix_df[, gene := factor(gene, levels = gene_order)]
matrix_df[, cohort_label := factor(dataset, levels = c("GSE174554", "GSE274546", "GSE121810", "GSE154795"))]
write_csv(matrix_df, "figure5e_four_cohort_leading_edge_logfc", figure5_source_dir)

save_panel(
  build_leading_edge_heatmap(matrix_df),
  "fig5e_four_cohort_leading_edge_heatmap",
  figure5_dir,
  118,
  115
)

entry_fgsea_df <- rbindlist(list(
  data.frame(row_label = "Auto strict IDH-wt", NES = gse154795_main$fgsea$NES[1], p = gse154795_main$fgsea$pval[1]),
  data.frame(row_label = "Reviewed IDH-wt", NES = gse154795_reviewed$fgsea$NES[1], p = gse154795_reviewed$fgsea$pval[1]),
  data.frame(row_label = "Inclusive IDH-wt", NES = gse154795_inclusive$fgsea$NES[1], p = gse154795_inclusive$fgsea$pval[1])
))
entry_fgsea_df[, p_label := vapply(p, build_p_label, character(1))]
write_csv(entry_fgsea_df, "supp7a_entry_fgsea_sensitivity", s7_source_dir)
save_panel(
  build_nes_plot(
    entry_fgsea_df,
    x_label = "fgsea NES",
    title_text = "GSE154795 entry sensitivity",
    subtitle_text = "Only fgsea shown here; camera stays in the separate algorithm table"
  ),
  "supp7a_entry_fgsea_sensitivity",
  s7_dir,
  125,
  72
)

idh_fgsea_df <- rbindlist(list(
  data.frame(row_label = "Strict IDH-wt", NES = gse154795_main$fgsea$NES[1], p = gse154795_main$fgsea$pval[1]),
  data.frame(row_label = "All IDH", NES = gse154795_all_idh$fgsea$NES[1], p = gse154795_all_idh$fgsea$pval[1])
))
idh_fgsea_df[, p_label := vapply(p, build_p_label, character(1))]
write_csv(idh_fgsea_df, "supp7b_idh_fgsea_sensitivity", s7_source_dir)
save_panel(
  build_nes_plot(
    idh_fgsea_df,
    x_label = "fgsea NES",
    title_text = "GSE154795 IDH sensitivity",
    subtitle_text = "All-IDH inclusion weakens the coordinated shift"
  ),
  "supp7b_idh_fgsea_sensitivity",
  s7_dir,
  125,
  62
)

score_sensitivity_df <- rbindlist(list(
  gse154795_main$score_models[, .(scenario = "Auto strict IDH-wt", metric_label, effect, conf_low, conf_high, fdr)],
  gse154795_reviewed$score_models[, .(scenario = "Reviewed IDH-wt", metric_label, effect, conf_low, conf_high, fdr)],
  gse154795_inclusive$score_models[, .(scenario = "Inclusive IDH-wt", metric_label, effect, conf_low, conf_high, fdr)],
  gse154795_all_idh$score_models[, .(scenario = "All IDH", metric_label, effect, conf_low, conf_high, fdr)]
))
score_sensitivity_df[, row_label := scenario]
score_sensitivity_df[, n_label := ""]
score_sensitivity_df[, fdr_label := vapply(fdr, build_p_label, character(1))]
write_csv(score_sensitivity_df, "supp7c_score_sensitivity_models", s7_source_dir)
save_panel(
  build_forest_plot(
    score_sensitivity_df,
    x_label = "Treatment-associated score shift",
    point_colors = c("Raw score" = pal["primary"], "Myeloid-adjusted" = pal["recurrent"]),
    facet_var = "metric_label",
    title_text = "GSE154795 score sensitivity",
    subtitle_text = "Raw and myeloid-adjusted effects stay on their own score scale"
  ),
  "supp7c_score_sensitivity",
  s7_dir,
  145,
  106
)

algorithm_table_df <- rbindlist(list(
  data.frame(
    cohort_label = "GSE121810",
    analysis_label = "Main exposure comparison",
    n_label = sprintf("%d / %d", gse121810$score_models[1, n_positive], gse121810$score_models[1, n_negative]),
    fgsea_label = sprintf("NES %.2f; P %s", gse121810$fgsea$NES[1], build_p_label(gse121810$fgsea$pval[1])),
    camera_label = sprintf("%s; P %s", gse121810$camera$Direction[1], build_p_label(gse121810$camera$PValue[1]))
  ),
  data.frame(
    cohort_label = "GSE154795",
    analysis_label = "Auto strict IDH-wt",
    n_label = sprintf("%d / %d", gse154795_main$score_models[1, n_positive], gse154795_main$score_models[1, n_negative]),
    fgsea_label = sprintf("NES %.2f; P %s", gse154795_main$fgsea$NES[1], build_p_label(gse154795_main$fgsea$pval[1])),
    camera_label = sprintf("%s; P %s", gse154795_main$camera$Direction[1], build_p_label(gse154795_main$camera$PValue[1]))
  ),
  data.frame(
    cohort_label = "GSE154795",
    analysis_label = "Reviewed IDH-wt",
    n_label = sprintf("%d / %d", gse154795_reviewed$score_models[1, n_positive], gse154795_reviewed$score_models[1, n_negative]),
    fgsea_label = sprintf("NES %.2f; P %s", gse154795_reviewed$fgsea$NES[1], build_p_label(gse154795_reviewed$fgsea$pval[1])),
    camera_label = sprintf("%s; P %s", gse154795_reviewed$camera$Direction[1], build_p_label(gse154795_reviewed$camera$PValue[1]))
  ),
  data.frame(
    cohort_label = "GSE154795",
    analysis_label = "Inclusive IDH-wt",
    n_label = sprintf("%d / %d", gse154795_inclusive$score_models[1, n_positive], gse154795_inclusive$score_models[1, n_negative]),
    fgsea_label = sprintf("NES %.2f; P %s", gse154795_inclusive$fgsea$NES[1], build_p_label(gse154795_inclusive$fgsea$pval[1])),
    camera_label = sprintf("%s; P %s", gse154795_inclusive$camera$Direction[1], build_p_label(gse154795_inclusive$camera$PValue[1]))
  ),
  data.frame(
    cohort_label = "GSE154795",
    analysis_label = "All IDH",
    n_label = sprintf("%d / %d", gse154795_all_idh$score_models[1, n_positive], gse154795_all_idh$score_models[1, n_negative]),
    fgsea_label = sprintf("NES %.2f; P %s", gse154795_all_idh$fgsea$NES[1], build_p_label(gse154795_all_idh$fgsea$pval[1])),
    camera_label = sprintf("%s; P %s", gse154795_all_idh$camera$Direction[1], build_p_label(gse154795_all_idh$camera$PValue[1]))
  )
))
write_csv(algorithm_table_df, "supp7d_algorithm_summary_table", s7_source_dir)
write_markdown_table(
  algorithm_table_df[, .(
    Cohort = cohort_label,
    Analysis = analysis_label,
    `n exp / ctrl` = n_label,
    fgsea = fgsea_label,
    camera = camera_label
  )],
  file.path(s7_write_dir, "SupplementaryFigure7D_algorithm_summary.md"),
  title_text = "Supplementary Figure 7D algorithm consistency table",
  intro_lines = c(
    "This algorithm summary is kept as a formal supplementary table and is not exported as a figure panel.",
    "It reports fgsea and camera side by side without forcing them onto a shared axis."
  )
)
unlink(file.path(s7_dir, c("supp7d_algorithm_summary_table.pdf", "supp7d_algorithm_summary_table.png")))

write_csv(gse121810$de, "GSE121810_raw20_edger", figure5_write_dir)
write_csv(gse121810$fgsea, "GSE121810_raw20_fgsea", figure5_write_dir)
write_csv(gse121810$camera, "GSE121810_raw20_camera", figure5_write_dir)
write_csv(gse121810$score_models, "GSE121810_raw20_score_models", figure5_write_dir)
write_csv(gse121810$coverage, "GSE121810_raw20_coverage", figure5_write_dir)

write_csv(gse154795_main$de, "GSE154795_raw20_edger_auto_strict_idhwt", figure5_write_dir)
write_csv(gse154795_main$fgsea, "GSE154795_raw20_fgsea_auto_strict_idhwt", figure5_write_dir)
write_csv(gse154795_main$camera, "GSE154795_raw20_camera_auto_strict_idhwt", figure5_write_dir)
write_csv(gse154795_main$score_models, "GSE154795_raw20_score_models_auto_strict_idhwt", figure5_write_dir)
write_csv(gse154795_main$coverage, "GSE154795_raw20_coverage_auto_strict_idhwt", figure5_write_dir)

write_csv(gse154795_reviewed$fgsea, "GSE154795_raw20_fgsea_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_reviewed$camera, "GSE154795_raw20_camera_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_reviewed$score_models, "GSE154795_raw20_score_models_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$fgsea, "GSE154795_raw20_fgsea_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$camera, "GSE154795_raw20_camera_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$score_models, "GSE154795_raw20_score_models_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_all_idh$fgsea, "GSE154795_raw20_fgsea_auto_all_idh", s7_write_dir)
write_csv(gse154795_all_idh$camera, "GSE154795_raw20_camera_auto_all_idh", s7_write_dir)
write_csv(gse154795_all_idh$score_models, "GSE154795_raw20_score_models_auto_all_idh", s7_write_dir)

figure5_manifest <- data.frame(
  panel = c("Figure5A", "Figure5B", "Figure5C", "Figure5D", "Figure5E"),
  stem = c(
    "fig5a_treatment_timeline",
    "fig5b_gse121810_raw20_gsea",
    "fig5c_gse154795_raw20_gsea",
    "fig5d_patient_score_forest",
    "fig5e_four_cohort_leading_edge_heatmap"
  ),
  description = c(
    "Two anti-PD-1 cohort timelines",
    "GSE121810 full Miller raw20 GSEA",
    "GSE154795 full Miller raw20 GSEA",
    "Patient-level raw and myeloid-adjusted forest",
    "Four-cohort leading-edge logFC heatmap"
  ),
  stringsAsFactors = FALSE
)
figure5_manifest$png <- file.path(figure5_dir, paste0(figure5_manifest$stem, ".png"))
figure5_manifest$pdf <- file.path(figure5_dir, paste0(figure5_manifest$stem, ".pdf"))
write_csv(figure5_manifest, "Figure5_panel_manifest", figure5_write_dir)

s7_manifest <- data.frame(
  panel = c("SupplementaryFigure7A", "SupplementaryFigure7B", "SupplementaryFigure7C"),
  stem = c(
    "supp7a_entry_fgsea_sensitivity",
    "supp7b_idh_fgsea_sensitivity",
    "supp7c_score_sensitivity"
  ),
  description = c(
    "GSE154795 entry sensitivity on fgsea NES scale",
    "GSE154795 IDH sensitivity on fgsea NES scale",
    "GSE154795 raw and myeloid-adjusted score sensitivity forest"
  ),
  stringsAsFactors = FALSE
)
s7_manifest$png <- file.path(s7_dir, paste0(s7_manifest$stem, ".png"))
s7_manifest$pdf <- file.path(s7_dir, paste0(s7_manifest$stem, ".pdf"))
write_csv(s7_manifest, "SupplementaryFigure7_panel_manifest", s7_write_dir)

summary_lines <- c(
  "# Step41 Figure5 / SupplementaryFigure7 rebuild",
  "",
  "- Figure5 keeps only the anti-PD-1 exposure-associated coordinated shift story.",
  "- Main panels use Miller raw20 only; no project-curated 16-gene replacement.",
  sprintf("- GSE121810 fgsea NES %.3f, nominal P %s; camera %s, nominal P %s.", gse121810$fgsea$NES[1], build_p_label(gse121810$fgsea$pval[1]), gse121810$camera$Direction[1], build_p_label(gse121810$camera$PValue[1])),
  sprintf("- GSE154795 strict IDH-wt fgsea NES %.3f, nominal P %s; camera %s, nominal P %s.", gse154795_main$fgsea$NES[1], build_p_label(gse154795_main$fgsea$pval[1]), gse154795_main$camera$Direction[1], build_p_label(gse154795_main$camera$PValue[1])),
  sprintf("- Main forest uses model-based 95%% CI; GSE121810 myeloid-adjusted FDR %s, GSE154795 myeloid-adjusted FDR %s.", build_p_label(gse121810$score_models[metric == "raw20_score_myeloid_adjusted", fdr][1]), build_p_label(gse154795_main$score_models[metric == "raw20_score_myeloid_adjusted", fdr][1])),
  "- Historical legacy16 response tables are excluded from the current rebuild and kept only in the obsolete archive.",
  "- Supplementary Figure 7D is kept as CSV plus markdown table only; no PNG or PDF is exported."
)
writeLines(summary_lines, file.path(write_root, "Figure5_SupplementaryFigure7_summary.md"))

cat("STEP41_FIG5_PD1_COMPLETE panels=8\n")
