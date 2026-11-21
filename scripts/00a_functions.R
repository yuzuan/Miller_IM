# ======================================================================
# 00a_functions.R —— Miller_IM 公共函数库
# Miller_IM：复发胶质母细胞瘤公开分析代码
#
# 目的：
# 统一多个 R 脚本共用的计数输入检查、固定程序富集和图件导出。
# 只归并已有实现；原分析阈值、输出字段与随机数顺序保持不变。
# 本函数脚本不用单独浏览学习，在其它分析脚本调用到它时再学。
# 用 RStudio 打开本脚本，在调用处 Ctrl + 左键（Mac 为 Cmd + 左键）点函数名可跳转到定义。
# 函数定义的查找靠上面的跳转，不靠文件头清单。

# 输入：基因×样本的 gzip CSV，首列为行名；输出：保留行列顺序的整数矩阵。
# 非有限、非整数、负计数或重复基因时停止；1e-8 是原有整数容差。
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

# 输入：命名基因集列表和命名排序向量；输出：fgsea 结果 data.frame，leadingEdge 为分号文本。
# 供独立队列与 LOPO 共用；原种子 42、5–500 基因边界、eps=0 和 SerialParam 均保留。
# 调用脚本已加载 fgsea 和 BiocParallel；这里不更换并行后端或额外消耗随机数。
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

# 输入：已保存的单页 PDF、PNG 目标及分辨率；输出：无返回值，写 PNG 并核对非空。
# pdftools 是明确的公开依赖；两个绘图入口均从各自 PDF 母版转换，不重画 PNG。
render_pdf_png <- function(pdf_path, png_path, dpi = 600) {
    if (!file.exists(pdf_path)) stop("Missing PDF: ", pdf_path)
    dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
    pdftools::pdf_convert(pdf_path, format = "png", pages = 1, filenames = png_path, dpi = dpi, verbose = FALSE)
    output_paths <- c(pdf_path, png_path)
    stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
}

# 输入：图对象、文件名、输出目录及英寸尺寸；输出：无返回值，保存 PDF 与派生 PNG。
# 图对象在预览前保留随机状态；PDF 使用原首轮绘图状态，第二次预览替代旧 PNG 的随机绘图消耗。
save_plot_pair <- function(plot, stem, dir_path, width, height) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    plot_seed <- .Random.seed
    pdf(NULL)
    print(plot)
    assign(".Random.seed", plot_seed, envir = .GlobalEnv)
    pdf_path <- file.path(dir_path, paste0(stem, ".pdf"))
    ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf)

    # 旧代码在写 PNG 时还会绘制一次；保留该次随机消耗，避免影响后面的 jitter。
    print(plot)
    dev.off()
    normalise_pdf_metadata(pdf_path)
    png_path <- file.path(dir_path, paste0(stem, ".png"))
    render_pdf_png(pdf_path, png_path, dpi = 600)
}

# 输入：PDF 文件；输出：无返回值，沿用原 pypdf 元数据规则替换文件；失败立即停止。
# Python 是原公开实现的既有依赖；shQuote 保证含空格的项目路径作为完整参数传递。
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

# 输入：按优先级排列的候选路径及标签；输出：第一个存在路径的绝对形式；无输入时报错。
locate_required <- function(candidates, label) {
    existing <- candidates[file.exists(candidates)]
    if (!length(existing)) stop(label, " not found. Checked: ", paste(candidates, collapse = "; "))
    result <- normalizePath(existing[1], mustWork = TRUE)
    return(result)
}

# 输入：一组数值；输出：均值及原定义 1.96 标准误区间的命名向量，无有效值时为 NA。
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

# 输入：面板标识、文件名、说明和目录；输出：保留原五列的单行清单表。
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

# 输入：图对象、文件名、目录及毫米尺寸；输出：无返回值，保留原 PDF 设备参数并派生 PNG。
# 预览保持 PDF 原随机状态；第二次绘图消耗替代旧版直接输出 PNG 的绘图，不更改后续序列。
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

# 输入：结果表、文件名和目录；输出：无返回值，按既有 fwrite 格式写 CSV 并核对非空。
write_csv <- function(input_table, path_stem, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(out_dir, paste0(path_stem, ".csv"))
    fwrite(as.data.table(input_table), output_path)
    stopifnot(file.exists(output_path), file.info(output_path)$size > 0)
}

# 输入：结果表、文件名和目录；输出：无返回值，沿用 fwrite 的 gzip 写出并核对非空。
write_csv_gz <- function(input_table, path_stem, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(out_dir, paste0(path_stem, ".csv.gz"))
    fwrite(as.data.table(input_table), output_path)
    stopifnot(file.exists(output_path), file.info(output_path)$size > 0)
}

# 输入：表、Markdown 路径、标题和说明；输出：无返回值，写原列序的 Markdown 表。
write_markdown_table <- function(input_table, file_path, title_text, intro_lines = character()) {
    display_table <- as.data.frame(input_table, stringsAsFactors = FALSE)
    header <- paste(colnames(display_table), collapse = " | ")
    separator <- paste(rep("---", ncol(display_table)), collapse = " | ")
    display_matrix <- as.matrix(display_table)
    body <- character(nrow(display_table))
    for (index in seq_len(nrow(display_table))) body[index] <- paste(display_matrix[index, ], collapse = " | ")
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
    stopifnot(file.exists(file_path), file.info(file_path)$size > 0)
}

# 输入：原始 RDS 路径；输出：R 对象。按原顺序尝试普通、gzip、双 gzip，三次均失败才停止。
# error 处理直接返回条件对象，由主流程判断，不再嵌套三层匿名错误处理。
read_rds_auto <- function(path) {
    result <- tryCatch(readRDS(path), error = identity)
    if (!inherits(result, "error")) return(result)
    result <- tryCatch(readRDS(gzfile(path, "rb")), error = identity)
    if (!inherits(result, "error")) return(result)
    result <- tryCatch(readRDS(gzcon(gzfile(path, "rb"))), error = identity)
    if (inherits(result, "error")) stop("Cannot read RDS file: ", path, "\n", conditionMessage(result))
    return(result)
}

# 输入：基因×样本表达矩阵；输出：逐基因 z 分数矩阵，保留原无效标准差置 NA/最终非有限值置 0。
safe_zscore_rows <- function(x) {
    m <- rowMeans(x, na.rm = TRUE)
    s <- apply(x, 1, sd, na.rm = TRUE)
    s[!is.finite(s) | s == 0] <- NA_real_
    z <- sweep(x, 1, m, "-")
    z <- sweep(z, 1, s, "/")
    z[!is.finite(z)] <- 0
    return(z)
}

# 输入：P 值；输出：原阈值和小数位格式的字符串，无有效 P 值时返回 NA 文本。
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

# 输入：线性模型、系数名称及队列说明；输出：原效应、95% CI、P 值和附加字段的单行表。
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

# 输入：基础字号；输出：沿用调用脚本 pal 配色的 ggplot 主题；不修改会话选项。
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
