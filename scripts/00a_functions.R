# ======================================================================
# 00a_functions.R — Miller-IM shared R functions
# Paired glioblastoma recurrence and external cohort analyses.
#
# Purpose
# Share count input, fixed-program enrichment, and figure exports between the R workflows.
# Preserve test settings, output fields, and random number order across callers.
# The analysis scripts show which functions are used at each step.
# Use go-to-definition in the editor to read a function from its call site.
# The signatures document the parameters passed by each workflow.

# Read a gzip CSV of gene counts with gene identifiers in the first column.
# Use the original count tolerance of 1e-8 and reject nonfinite or duplicate entries.
read_count_matrix <- function(path) {
    count_table <- read.csv(gzfile(path), row.names = 1, check.names = FALSE)
    count_matrix <- as.matrix(count_table)
    if (any(!is.finite(count_matrix)) || any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
        stop("Non-integer counts in ", path)
    }
    storage.mode(count_matrix) <- "integer"
    if (anyDuplicated(rownames(count_matrix))) stop("Duplicated genes in ", path)
    if (any(count_matrix < 0L)) stop("Negative counts in ", path)
    return(count_matrix)
}

# Run fgsea from a named gene set list and a named ranked vector. Join leading-edge genes with semicolons.
# The recurrence and leave-one-out analyses use seed 42, sizes 5–500, eps = 0, and SerialParam.
# The caller loads fgsea and BiocParallel before sourcing these functions.
run_fixed_gsea <- function(pathways, rank_stat) {
    tested <- vector("list", length(pathways))
    names(tested) <- names(pathways)
    for (set_name in names(pathways)) tested[[set_name]] <- intersect(unique(pathways[[set_name]]), names(rank_stat))
    tested <- tested[lengths(tested) >= 5]

    set.seed(42)
    result_table <- fgseaMultilevel(
        pathways = tested,
        stats    = rank_stat,
        minSize  = 5,
        maxSize  = 500,
        eps      = 0,
        BPPARAM  = SerialParam(progressbar = FALSE)
    )
    result_table <- as.data.frame(result_table)
    result_table$leadingEdge <- vapply(result_table$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    return(result_table)
}

# Convert a saved one-page PDF to a PNG at the requested resolution and verify it is nonempty.
# The proteome and PD-1 plots use pdftools to render the saved single-page PDF.
render_pdf_png <- function(pdf_path, png_path, dpi = 600) {
    if (!file.exists(pdf_path)) stop("Missing PDF: ", pdf_path)
    dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
    pdftools::pdf_convert(pdf_path, format = "png", pages = 1, filenames = png_path, dpi = dpi, verbose = FALSE)
    output_paths <- c(pdf_path, png_path)
    stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
}

# Save a plot as a PDF and render its PNG using the supplied dimensions in inches.
# Save the first draw's random state for the PDF. A second preview preserves the next plot's jitter sequence.
save_plot_pair <- function(plot, stem, dir_path, width, height) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    plot_seed <- .Random.seed
    pdf(NULL)
    print(plot)
    assign(".Random.seed", plot_seed, envir = .GlobalEnv)
    pdf_path <- file.path(dir_path, paste0(stem, ".pdf"))
    ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf)

    # The original PNG render drew the plot a second time. Consume the same draws here.
    print(plot)
    dev.off()
    normalise_pdf_metadata(pdf_path)
    png_path <- file.path(dir_path, paste0(stem, ".png"))
    render_pdf_png(pdf_path, png_path, dpi = 600)
}

# Replace the PDF metadata with the original pypdf settings.
# Quote the Python argument so a PDF path containing spaces is passed as a single value.
normalise_pdf_metadata <- function(pdf_path) {
    py_script <- tempfile(pattern = "normalise_pdf_", fileext = ".py")
    python_lines <- c(
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
    )
    writeLines(python_lines, py_script, useBytes = TRUE)
    on.exit(unlink(py_script), add = TRUE)
    tmp_path <- paste0(pdf_path, ".tmp")
    status <- system2("python3", shQuote(c(py_script, pdf_path, tmp_path)), stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L) || !file.exists(tmp_path)) stop("Failed to normalise PDF metadata for ", pdf_path)
    if (file.exists(pdf_path)) unlink(pdf_path)
    ok <- file.rename(tmp_path, pdf_path)
    if (!ok) stop("Failed to replace PDF after metadata normalisation: ", pdf_path)
}

# Return the first existing candidate as an absolute path, or raise an error naming the input.
locate_required <- function(candidates, label) {
    existing <- candidates[file.exists(candidates)]
    if (!length(existing)) stop(label, " not found. Checked: ", paste(candidates, collapse = "; "))
    result <- normalizePath(existing[1], mustWork = TRUE)
    return(result)
}

# Return the mean and the original interval of 1.96 standard errors. Use NA when no values are valid.
ci_mean <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (!n) {
        result <- c(mean = NA_real_, low = NA_real_, high = NA_real_)
        return(result)
    }
    se <- stats::sd(x) / sqrt(n)
    m <- mean(x)
    result <- c(mean = m, low = m - 1.96 * se, high = m + 1.96 * se)
    return(result)
}

# Return the original five manifest columns for a single panel.
panel_manifest <- function(panel, stem, description, out_dir) {
    result_table <- data.table(
        panel = panel,
        stem  = stem,
        description = description,
        pdf   = file.path(out_dir, paste0(stem, ".pdf")),
        png   = file.path(out_dir, paste0(stem, ".png"))
    )
    return(result_table)
}

# Save a plot in millimeters using the original PDF device settings and render the PNG.
# Use the first draw's random state for the PDF and consume the original PNG draw in the second preview.
save_panel <- function(plot_obj, stem, out_dir, width_mm, height_mm) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    plot_seed <- .Random.seed
    pdf(NULL)
    print(plot_obj)
    assign(".Random.seed", plot_seed, envir = .GlobalEnv)
    pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
    pdf(
        file        = pdf_file,
        width       = width_mm / 25.4,
        height      = height_mm / 25.4,
        onefile     = FALSE,
        paper       = "special",
        bg          = "white",
        useDingbats = FALSE,
        compress    = TRUE,
        timestamp   = FALSE,
        producer    = "R grDevices",
        author      = "Codex"
    )
        print(plot_obj)
    dev.off()
    print(plot_obj)
    dev.off()
    png_file <- file.path(out_dir, paste0(stem, ".png"))
    render_pdf_png(pdf_file, png_file, dpi = 600)
}

# Write the table with the published fwrite CSV format and verify the file is nonempty.
write_csv <- function(input_table, path_stem, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(out_dir, paste0(path_stem, ".csv"))
    fwrite(as.data.table(input_table), output_path)
    stopifnot(file.exists(output_path), file.info(output_path)$size > 0)
}

# Write the gzip CSV with the published fwrite settings and verify the file is nonempty.
write_csv_gz <- function(input_table, path_stem, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(out_dir, paste0(path_stem, ".csv.gz"))
    fwrite(as.data.table(input_table), output_path)
    stopifnot(file.exists(output_path), file.info(output_path)$size > 0)
}

# Read the cached object as plain, gzip, or double-gzip RDS, in that order.
# Return the error condition to the caller so the next reader can be attempted.
read_rds_auto <- function(path) {
    result <- tryCatch(readRDS(path), error = identity)
    if (!inherits(result, "error")) return(result)
    result <- tryCatch(readRDS(gzfile(path, "rb")), error = identity)
    if (!inherits(result, "error")) return(result)
    result <- tryCatch(readRDS(gzcon(gzfile(path, "rb"))), error = identity)
    if (inherits(result, "error")) stop("Cannot read RDS file: ", path, "\n", conditionMessage(result))
    return(result)
}

# Compute gene-wise z scores, then set nonfinite standardized values to zero.
safe_zscore_rows <- function(x) {
    m <- rowMeans(x, na.rm = TRUE)
    s <- apply(x, 1, sd, na.rm = TRUE)
    s[!is.finite(s) | s == 0] <- NA_real_
    z <- sweep(x, 1, m, "-")
    z <- sweep(z, 1, s, "/")
    z[!is.finite(z)] <- 0
    return(z)
}

# Format finite P values at the original threshold and precision. Display missing values as NA.
build_p_label <- function(p) {
    if (!is.finite(p) || is.na(p)) {
        result <- "NA"
    } else if (p < 1e-3) {
        result <- formatC(p, format = "e", digits = 2)
    } else {
        result <- formatC(p, format = "f", digits = 3)
    }
    return(result)
}

# Return one model summary row, including the specified effect, its 95% CI, and the P value.
extract_model_row <- function(
    model,
    coef_name,
    cohort,
    metric,
    metric_label,
    cohort_label,
    comparator,
    n_positive,
    n_negative,
    extra = list()
) {
    coefficient_matrix <- summary(model)$coefficients
    ci <- suppressMessages(confint(model, parm = coef_name, level = 0.95))
    model_table <- data.frame(
        cohort       = cohort,
        cohort_label = cohort_label,
        metric       = metric,
        metric_label = metric_label,
        comparator   = comparator,
        n_positive   = n_positive,
        n_negative   = n_negative,
        effect       = unname(coef(model)[coef_name]),
        conf_low     = unname(ci[1]),
        conf_high    = unname(ci[2]),
        p_value      = unname(coefficient_matrix[coef_name, "Pr(>|t|)"]),
        stringsAsFactors = FALSE
    )
    if (length(extra)) for (column in names(extra)) model_table[, column] <- extra[[column]]
    return(model_table)
}

# Create a ggplot theme from the requested type size and the caller's color palette.
theme_panel <- function(base_size = 7.5) {
    result_theme <- theme_classic(base_size = base_size, base_family = "sans")
    result_theme <- result_theme + theme(
        plot.title       = element_text(color = pal["text"], face = "bold", size = base_size + 0.7, hjust = 0),
        plot.subtitle    = element_text(color = pal["neutral"], size = base_size - 0.5, hjust = 0),
        axis.title       = element_text(color = pal["text"]),
        axis.text        = element_text(color = pal["text"]),
        axis.line        = element_line(linewidth = 0.35, color = pal["text"]),
        axis.ticks       = element_line(linewidth = 0.35, color = pal["text"]),
        legend.title     = element_blank(),
        legend.key       = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(color = pal["text"], face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    )
    return(result_theme)
}
