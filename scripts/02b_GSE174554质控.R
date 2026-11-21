#!/usr/bin/env Rscript
# 目的：按公开规则合并重测文库，并在每个 capture 内完成质控与 scDblFinder 分类。
# 输入：MILLER_IM_DATA_ROOT 下 results/manifests/gse174554_libraries.csv 和 GSE174554/GSE174554_RAW。
# 输出：write/34_gse174554_raw_independent_discovery/01_qc_doublets 下的分类、capture 与运行汇总 CSV。
# 输出还包括 batch2_barcode_overlap_audit.csv；所有文件名、列名及修复矩阵入口均为存量公开接口。
# 运行：Rscript scripts/02b_GSE174554质控.R [project_root]。
# 保留 R/Bioconductor 实现；每个 capture 的种子位置及 SerialParam 不改，以维持原随机序列。

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(SingleCellExperiment)
    library(scDblFinder)
    library(BiocParallel)
})

# 全局种子保留在原入口位置；后面每个 capture 的种子仍在 scDblFinder 前单独设置。
set.seed(174554)
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
source_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root, "data"))
manifest_path <- file.path(source_root, "results", "manifests", "gse174554_libraries.csv")
manifest_table <- as.data.frame(fread(manifest_path))

# 91 个矩阵和 81 个 GSM 是冻结输入清单的约束，保留旧检查而非放宽队列。
stopifnot(nrow(manifest_table) == 91L, uniqueN(manifest_table$gsm_id) == 81L)
manifest_table <- manifest_table[manifest_table$progression %in% c("Primary", "1st Recurrent"), ]
raw_root <- file.path(source_root, "GSE174554", "GSE174554_RAW")
for (column in c("matrix_path", "features_path", "barcodes_path")) {
    manifest_table[, column] <- file.path(raw_root, basename(manifest_table[, column]))
    missing_source <- manifest_table[!file.exists(manifest_table[, column]), column]
    if (length(missing_source)) stop("Missing relocated source files: ", paste(missing_source, collapse = ", "))
}
manifest_order <- order(manifest_table$gsm_id, manifest_table$tech_rep, na.last = FALSE, method = "radix")
manifest_table <- manifest_table[manifest_order, ]

# 每份列表保存一个 capture；输出目录也提供旧版人工修复矩阵的相对入口。
out_dir <- file.path(root, "write", "34_gse174554_raw_independent_discovery", "01_qc_doublets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
classification_rows <- list()
summary_rows <- list()
overlap_rows <- list()
constructed_captures <- 0L

# 输入：capture 标识、单行文库资料、已合并计数和合并方式；输出：分类表及单行汇总表。
# 低于原定 50 个合格细胞时只保留排除汇总；无效 scDblFinder 分类直接报错。
process_capture <- function(capture_id, library_table, capture_data, merge_mode) {
    count_matrix <- capture_data$counts
    symbols <- sub("\\.[0-9]+$", "", rownames(count_matrix))
    is_mt <- startsWith(toupper(symbols), "MT-")
    n_genes <- colSums(count_matrix > 0)
    total_counts <- colSums(count_matrix)
    mt_counts <- rep(0, ncol(count_matrix))
    if (any(is_mt)) mt_counts <- colSums(count_matrix[is_mt, , drop = FALSE])
    pct_mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)

    # 原质控阈值为至少 200 个基因、线粒体比例至多 3%；不改变边界方向。
    qc_keep <- n_genes >= 200L & pct_mt <= 3
    n_qc <- sum(qc_keep)
    if (n_qc < 50L) {
        summary_table <- data.frame(
            capture_id   = capture_id,
            gsm_id       = library_table$gsm_id,
            patient_id   = library_table$patient_id,
            pair_id      = library_table$pair_id,
            progression  = library_table$progression,
            merge_mode   = merge_mode,
            n_raw        = ncol(count_matrix),
            n_qc         = n_qc,
            n_doublet    = NA_integer_,
            n_singlet    = 0L,
            doublet_rate = NA_real_,
            status       = "excluded_too_few_qc_cells"
        )
        result <- list(classification = NULL, summary = summary_table)
        return(result)
    }

    # capture 内重设种子是原实现的一部分；移出会改变双细胞结果。SerialParam 保持兼容。
    qc_matrix <- count_matrix[, qc_keep, drop = FALSE]
    sce <- SingleCellExperiment(list(counts = qc_matrix))
    set.seed(174554 + sum(utf8ToInt(capture_id)))
    sce <- scDblFinder(sce, BPPARAM = SerialParam(progressbar = FALSE))
    cell_class <- as.character(colData(sce)$scDblFinder.class)
    score <- as.numeric(colData(sce)$scDblFinder.score)
    if (anyNA(cell_class) || !all(cell_class %in% c("singlet", "doublet"))) {
        stop("scDblFinder returned invalid classes for ", capture_id)
    }

    # 一行对应一个合格细胞；scDblFinder_score/class 是下游既有字段，保留大小写。
    classification_table <- data.frame(
        capture_id        = capture_id,
        gsm_id            = as.character(library_table$gsm_id),
        patient_id        = as.character(library_table$patient_id),
        pair_id           = as.character(library_table$pair_id),
        progression       = as.character(library_table$progression),
        barcode           = colnames(qc_matrix),
        n_genes_by_counts = as.integer(n_genes[qc_keep]),
        total_counts      = as.integer(total_counts[qc_keep]),
        pct_counts_mt     = as.numeric(pct_mt[qc_keep]),
        scDblFinder_score = score,
        scDblFinder_class = cell_class
    )
    summary_table <- data.frame(
        capture_id   = capture_id,
        gsm_id       = library_table$gsm_id,
        patient_id   = library_table$patient_id,
        pair_id      = library_table$pair_id,
        progression  = library_table$progression,
        merge_mode   = merge_mode,
        n_raw        = ncol(count_matrix),
        n_qc         = n_qc,
        n_doublet    = sum(cell_class == "doublet"),
        n_singlet    = sum(cell_class == "singlet"),
        doublet_rate = mean(cell_class == "doublet"),
        status       = "ok"
    )
    result <- list(classification = classification_table, summary = summary_table)
    return(result)
}
for (gsm in unique(manifest_table$gsm_id)) {
    library_table <- manifest_table[manifest_table$gsm_id == gsm, , drop = FALSE]

    # 输入：一行文库资料；输出：计数矩阵、基因表、条码组成的列表；尺寸不符时停止。
    # SF8963 必须使用既有修复文件，其余文库沿用清单路径，不重新修补原始矩阵。
    read_library <- function(source_table) {
        feature_table <- fread(cmd = paste("gzip -cd", shQuote(source_table$features_path)), header = FALSE)
        feature_table <- as.data.frame(feature_table)
        gene_symbols <- as.character(feature_table[, 1L])
        if (ncol(feature_table) > 1L) gene_symbols <- as.character(feature_table[, 2L])
        feature_table <- data.frame(gene_id = as.character(feature_table[, 1L]), gene_symbol = gene_symbols)
        barcode_table <- fread(cmd = paste("gzip -cd", shQuote(source_table$barcodes_path)), header = FALSE)
        barcode_table <- as.data.frame(barcode_table)
        barcodes <- as.character(barcode_table[, 1L])

        # features/barcodes 原始文件没有列名，所以这里保留按位置取列。
        matrix_path <- source_table$matrix_path
        if (as.character(source_table$library_id) == "GSM5319529_SF8963") {
            matrix_path <- file.path(out_dir, "source_repairs", "GSM5319529_SF8963_matrix.fixed.mtx.gz")
            if (!file.exists(matrix_path)) stop("Missing repaired SF8963 matrix: ", matrix_path)
        }
        connection <- gzfile(matrix_path, open = "rt")
        count_matrix <- readMM(connection)
        close(connection)
        count_matrix <- as(count_matrix, "CsparseMatrix")
        if (!identical(dim(count_matrix), c(nrow(feature_table), length(barcodes)))) {
            stop("Matrix dimension mismatch: ", source_table$library_id)
        }
        rownames(count_matrix) <- make.unique(feature_table$gene_symbol)
        colnames(count_matrix) <- barcodes
        result <- list(counts = count_matrix, features = feature_table, barcodes = barcodes)
        return(result)
    }
    if (nrow(library_table) == 1L) {
        capture_id <- as.character(library_table$library_id[1L])
        capture_data <- read_library(library_table[1L, , drop = FALSE])
        result <- process_capture(capture_id, library_table[1L, , drop = FALSE], capture_data, "single_library")
        classification_rows[[capture_id]] <- result$classification
        summary_rows[[capture_id]] <- result$summary
        constructed_captures <- constructed_captures + 1L
        rm(result, capture_data)
        gc(verbose = FALSE)
        next
    }

    # 技术重复只允许 main/batch2 两份；读取及后续 scDblFinder 顺序与原实现相同。
    if (nrow(library_table) != 2L || !setequal(library_table$tech_rep, c("main", "batch2"))) {
        stop("Unexpected technical replicate structure for ", gsm)
    }
    main_table <- library_table[library_table$tech_rep == "main", , drop = FALSE]
    batch_table <- library_table[library_table$tech_rep == "batch2", , drop = FALSE]
    main <- read_library(main_table[1L, , drop = FALSE])
    batch2 <- read_library(batch_table[1L, , drop = FALSE])
    if (!identical(main$features, batch2$features)) stop("Feature tables differ within ", gsm)
    overlap <- intersect(main$barcodes, batch2$barcodes)
    overlap_fraction_min <- length(overlap) / min(length(main$barcodes), length(batch2$barcodes))

    # 重叠比例至少 0.80 时逐条码合并 UMI；至多 0.05 时按独立 capture 分别处理。
    if (overlap_fraction_min >= 0.80) {
        capture_id <- paste0(gsm, "__merged_resequencing")
        union_barcodes <- union(main$barcodes, batch2$barcodes)
        joined_matrix <- cbind(main$counts, batch2$counts)
        source_barcodes <- c(main$barcodes, batch2$barcodes)
        map_matrix <- sparseMatrix(
            i    = seq_along(source_barcodes),
            j    = match(source_barcodes, union_barcodes),
            x    = 1,
            dims = c(length(source_barcodes), length(union_barcodes))
        )
        merged_matrix <- as(joined_matrix %*% map_matrix, "dgCMatrix")
        colnames(merged_matrix) <- union_barcodes
        rownames(merged_matrix) <- rownames(main$counts)
        merged <- list(counts = merged_matrix, features = main$features, barcodes = union_barcodes)

        # 合并只改变 capture 的输入组织，不改变双细胞模型和种子。
        result <- process_capture(capture_id, main_table, merged, "merged_resequencing")
        classification_rows[[capture_id]] <- result$classification
        summary_rows[[capture_id]] <- result$summary
        constructed_captures <- constructed_captures + 1L
        decision <- "merge_umi_by_barcode"
        rm(merged, merged_matrix, joined_matrix, map_matrix, result)
    } else if (overlap_fraction_min <= 0.05) {
        entries <- list(list(row = main_table, data = main), list(row = batch_table, data = batch2))
        for (entry in entries) {
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

    # 一行对应一个技术重复 GSM，供既有下游核对合并依据。
    overlap_rows[[gsm]] <- data.frame(
        gsm_id       = gsm,
        sample_short = as.character(main_table$sample_short),
        n_main       = length(main$barcodes),
        n_batch2     = length(batch2$barcodes),
        n_overlap    = length(overlap),
        overlap_fraction_min = overlap_fraction_min,
        decision     = decision
    )
    rm(main, batch2)
    gc(verbose = FALSE)
}

# 按 capture 处理顺序合并；rbindlist 保留原有缺失列填充和 gzip 导出行为。
overlap_table_total <- as.data.frame(rbindlist(overlap_rows, fill = TRUE))
fwrite(overlap_table_total, file.path(out_dir, "batch2_barcode_overlap_audit.csv"))
classification_table_total <- as.data.frame(rbindlist(classification_rows, use.names = TRUE))
summary_table_total <- as.data.frame(rbindlist(summary_rows, use.names = TRUE, fill = TRUE))
fwrite(classification_table_total, file.path(out_dir, "scDblFinder_cell_classifications.csv.gz"), compress = "gzip")
fwrite(summary_table_total, file.path(out_dir, "capture_qc_doublet_summary.csv"))

# 原汇总字段逐项从本次输入和 capture 汇总计算，不写死实际细胞数。
run_summary_table <- data.frame(
    metric = c(
        "input_matrices", "input_gsm", "constructed_captures", "raw_cells_after_resequencing_merge",
        "qc_cells", "singlets"
    ),
    value = c(
        nrow(manifest_table), uniqueN(manifest_table$gsm_id), constructed_captures,
        sum(summary_table_total$n_raw), sum(summary_table_total$n_qc), sum(summary_table_total$n_singlet)
    )
)
fwrite(run_summary_table, file.path(out_dir, "qc_run_summary.csv"))

# 输出文件须非空，分类行数须对应实际进入双细胞检测的细胞数。
output_paths <- file.path(
    out_dir,
    c(
        "batch2_barcode_overlap_audit.csv", "scDblFinder_cell_classifications.csv.gz",
        "capture_qc_doublet_summary.csv", "qc_run_summary.csv"
    )
)
stopifnot(all(file.exists(output_paths)), all(file.info(output_paths)$size > 0))
stopifnot(nrow(classification_table_total) == sum(summary_table_total$n_qc[summary_table_total$status == "ok"]))
cat("captures=", constructed_captures, " singlets=", sum(summary_table_total$n_singlet), "\n", sep = "")
