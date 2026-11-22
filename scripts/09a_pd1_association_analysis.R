#!/usr/bin/env Rscript
# Purpose  Rebuild Figure 5 and Supplementary Figure 7 using the original PD-1 exposure models.
# Inputs   GSE121810 workbook and GSE154795 RDS and SOFT under MILLER_IM_DATA_ROOT, plus the write/38_* tables.
# Outputs  Source tables and manifests under write/41_mg_inflammatory_sci_rebuild/{Figure5,SupplementaryFigure7}.
#          The figures tree contains panel_library PDF and PNG files. The algorithm comparison is a CSV table.
# Run      Rscript scripts/09a_pd1_association_analysis.R [project_root] with R, edgeR, limma and fgsea.
# Keep the original fields, data.table groups, model designs, seeds and fit order. Render PNGs from PDFs.

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
})

set.seed(20260713)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot locate the current script.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else dirname(dirname(script_path))

source(file.path(root_dir, "scripts", "00a_functions.R"))

external_root <- Sys.getenv("MILLER_IM_DATA_ROOT", file.path(root_dir, "data"))
gse121810_source <- file.path(
    external_root,
    "write",
    "28_supportive_dataset_search",
    "source_metadata",
    "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx"
)

# Use the original GSE121810 workbook to produce edgeR, fgsea, camera, score-model and coverage results.

count_table <- as.data.frame(read_excel(gse121810_source), check.names = FALSE)
rownames(count_table) <- count_table[, 1L]
count_table[, 1L] <- NULL
count_matrix <- rowsum(as.matrix(count_table), group = rownames(count_table), reorder = FALSE)
storage.mode(count_matrix) <- "numeric"

# Determine neoadjuvant versus adjuvant exposure from the sample suffix. Retain filterByExpr thresholds of 5 and 15.
group <- factor(ifelse(grepl("_A$", colnames(count_matrix)), "Neo", "Adj"), levels = c("Adj", "Neo"))
design_matrix <- model.matrix(~group)
y <- DGEList(count_matrix)
keep <- filterByExpr(y, design = design_matrix, min.count = 5, min.total.count = 15)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
y <- estimateDisp(y, design_matrix, robust = TRUE)
fit_model <- glmQLFit(y, design_matrix, robust = TRUE)
qlf <- glmQLFTest(fit_model, coef = "groupNeo")
deg_table <- topTags(qlf, n = Inf, sort.by = "none")$table
deg_table$gene <- rownames(deg_table)
deg_table$FDR <- p.adjust(deg_table$PValue, method = "BH")

# Rank by sign(logFC) times F. Run the original single-process fgsea and camera tests.
ranks <- deg_table$F
names(ranks) <- deg_table$gene
ranks <- sort(ranks * sign(deg_table$logFC), decreasing = TRUE)

raw20 <- c(
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3"
)
fgsea_table <- as.data.table(fgseaMultilevel(
    pathways = list(Miller_raw20 = raw20),
    stats    = ranks,
    minSize  = 5,
    maxSize  = 500,
    eps      = 0,
    nproc    = 1
))
fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
fgsea_table[, cohort := "GSE121810"]
fgsea_table[, contrast := "Neoadjuvant_minus_AdjuvantOnly"]

v <- voom(y, design_matrix, plot = FALSE)
raw20_indices <- ids2indices(list(Miller_raw20 = raw20), rownames(v))
camera_table <- as.data.table(
    camera(v, raw20_indices, design_matrix, contrast = "groupNeo"),
    keep.rownames = "pathway"
)
if (!"FDR" %in% colnames(camera_table)) camera_table[, FDR := p.adjust(PValue, method = "BH")]
camera_table[, cohort := "GSE121810"]
camera_table[, contrast := "Neoadjuvant_minus_AdjuvantOnly"]

# Standardize log CPM by gene, average scores over the original gene set and fit both patient models.
myeloid_genes <- c(
    "TREM2", "TMEM119", "P2RY12", "CX3CR1", "AIF1", "TYROBP", "FCER1G", "CSF1R",
    "CD68", "C1QA", "C1QB", "C1QC", "APOC1", "GPNMB", "LGALS3", "FCGR1A"
)
logcpm <- cpm(y, log = TRUE, prior.count = 1)
z <- safe_zscore_rows(logcpm)
raw20_present <- intersect(raw20, rownames(z))
myeloid_present <- intersect(myeloid_genes, rownames(z))
score_table <- data.frame(
    sample            = colnames(z),
    cohort            = "GSE121810",
    group             = group,
    raw20_score       = colMeans(z[raw20_present, , drop = FALSE]),
    myeloid_abundance = colMeans(z[myeloid_present, , drop = FALSE]),
    stringsAsFactors  = FALSE
)

raw_model <- lm(raw20_score ~ group, data = score_table)
adj_model <- lm(raw20_score ~ myeloid_abundance + group, data = score_table)
score_model_table_total <- rbindlist(list(
    extract_model_row(
        raw_model,
        "groupNeo",
        "GSE121810",
        "raw20_score",
        "Raw score",
        "GSE121810",
        "Neo_minus_Adj",
        sum(score_table$group == "Neo"),
        sum(score_table$group == "Adj")
    ),
    extract_model_row(
        adj_model,
        "groupNeo",
        "GSE121810",
        "raw20_score_myeloid_adjusted",
        "Myeloid-adjusted",
        "GSE121810",
        "Neo_minus_Adj",
        sum(score_table$group == "Neo"),
        sum(score_table$group == "Adj")
    )
))
score_model_table_total[, fdr := p.adjust(p_value, method = "BH")]
score_model_table_total[, n_label := sprintf("%d / %d", n_positive, n_negative)]
score_model_table_total[, fdr_label := vapply(fdr, build_p_label, character(1))]
coverage_table <- data.frame(
    cohort           = "GSE121810",
    signature        = "Miller_raw20",
    n_defined        = length(raw20),
    n_present        = length(raw20_present),
    coverage         = length(raw20_present) / length(raw20),
    present_genes    = paste(raw20_present, collapse = ";"),
    missing_genes    = paste(setdiff(raw20, raw20_present), collapse = ";"),
    stringsAsFactors = FALSE
)
gse121810 <- list(
    de       = as.data.table(deg_table),
    fgsea    = fgsea_table,
    camera   = camera_table,
    scores   = as.data.table(score_table),
    score_models = score_model_table_total,
    coverage = coverage_table,
    ranks    = ranks
)
gse154795_rds <- file.path(
    external_root,
    "write",
    "26_pure_bioinformatics_dataset_rescue",
    "source_metadata",
    "GSE154795_GBM.AllCell.Integrated.Scaled.ClusterRes.0.1.rds.gz"
)
gse154795_obj <- read_rds_auto(gse154795_rds)

# Use the original Seurat object to derive counts, metadata and lineage annotations at the frozen thresholds.

counts <- gse154795_obj@assays$RNA@counts
metadata_table <- gse154795_obj@meta.data
stopifnot(ncol(counts) == nrow(metadata_table))
metadata_table$cell_id <- rownames(metadata_table)
metadata_table$cluster <- as.character(metadata_table$seurat_clusters)
lineages <- list(
    Myeloid      = c("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "CSF1R", "CD68", "C1QA", "C1QB", "C1QC"),
    T_NK         = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1"),
    B_cell       = c("MS4A1", "CD79A", "CD79B", "CD37", "CD22", "BANK1"),
    Plasma       = c("MZB1", "JCHAIN", "SDC1", "IGHG1", "IGHG3"),
    Neural_glial = c("GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "OLIG1", "OLIG2"),
    Erythroid    = c("HBA1", "HBA2", "HBB", "ALAS2", "GYPA")
)

# Summarize counts per cluster. Use standardized lineage scores for clusters with at least 100 cells.
cluster_factor <- factor(metadata_table$cluster, levels = sort(unique(metadata_table$cluster)))
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

# Keep the 100-cell threshold, first-tie max.col behavior and automatic myeloid definition.
substantial <- as.integer(table(cluster_factor)[rownames(lineage_score)]) >= 100
lineage_z <- scale(lineage_score[substantial, , drop = FALSE])
best <- colnames(lineage_score)[max.col(lineage_score, ties.method = "first")]
best[substantial] <- colnames(lineage_z)[max.col(lineage_z, ties.method = "first")]
audit_table <- data.frame(
    cluster = rownames(lineage_score),
    n_cells = as.integer(table(cluster_factor)[rownames(lineage_score)]),
    best_lineage = best,
    lineage_score,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
audit_table$is_auto_myeloid <- audit_table$best_lineage == "Myeloid" & audit_table$n_cells >= 100
gse154795_prep <- list(counts = counts, meta = metadata_table, audit = audit_table)
rm(gse154795_obj)
gc()

gse154795_soft <- file.path(
    external_root,
    "write",
    "26_pure_bioinformatics_dataset_rescue",
    "source_metadata",
    "GSE154795_family.soft.gz"
)

# Read GSE154795 SOFT and keep the first matched technical annotation for each sample.

soft <- readLines(gzfile(gse154795_soft), warn = FALSE)
starts <- grep("^\\^SAMPLE = ", soft)
ends <- c(starts[-1] - 1, length(soft))
rows <- vector("list", length(starts))
for (i in seq_along(starts)) {
    block <- soft[starts[i]:ends[i]]
    patterns <- c(
        "^\\^SAMPLE = ", "^!Sample_title = ",
        "^!Sample_characteristics_ch1 = sequencing_method: ",
        "^!Sample_characteristics_ch1 = gel_bead_version: "
    )
    values <- rep(NA_character_, length(patterns))
    for (index in seq_along(patterns)) {
        hits <- grep(patterns[index], block, value = TRUE)
        if (length(hits)) values[index] <- sub(patterns[index], "", hits[1])
    }
    rows[[i]] <- data.frame(
        ID                = values[2],
        GSM               = values[1],
        sequencing_method = values[3],
        gel_bead_version  = values[4],
        acquisition_wave  = ifelse(grepl("^GSM467", values[1]), "early", "late"),
        stringsAsFactors  = FALSE
    )
}
gse154795_technical_table_total <- unique(rbindlist(rows, fill = TRUE))

# Return the cluster IDs for the requested lineage entry. Stop for an unknown definition.
get_gse154795_cluster_selection <- function(audit_table, entry) {
    if (entry == "reviewed_scgbm") {
        result <- c("0", "1", "2", "5", "8")
        return(result)
    }
    if (entry == "strict") {
        result <- audit_table$cluster[audit_table$is_auto_myeloid]
        return(result)
    }
    if (entry == "erythroid_ambient_inclusive") {
        selected <- audit_table$best_lineage == "Erythroid" & audit_table$Myeloid >= 9.5 & audit_table$n_cells >= 100
        extra <- audit_table$cluster[selected]
        result <- unique(c(audit_table$cluster[audit_table$is_auto_myeloid], extra))
        return(result)
    }
    stop("Unknown GSE154795 entry: ", entry)
}

# Use counts, metadata, technical annotations and lineage and IDH definitions to return model and score results.
analyze_gse154795_variant <- function(prepped, technical_table, entry, idh_entry) {
    metadata_table <- prepped$meta
    counts <- prepped$counts
    audit_table <- prepped$audit
    selected_clusters <- get_gse154795_cluster_selection(audit_table, entry)
    keep <- metadata_table$cluster %in% selected_clusters
    myeloid_metadata_table <- metadata_table[keep, , drop = FALSE]
    myeloid_counts <- counts[, keep, drop = FALSE]

    # Sum myeloid counts per sample and match technical covariates to the corresponding sample.
    sample_factor <- factor(myeloid_metadata_table$ID, levels = unique(myeloid_metadata_table$ID))
    sample_design <- sparse.model.matrix(~ 0 + sample_factor)
    colnames(sample_design) <- levels(sample_factor)
    pseudobulk <- myeloid_counts %*% sample_design
    sample_cells <- as.integer(table(sample_factor)[colnames(pseudobulk)])

    sample_metadata_table_total <- unique(myeloid_metadata_table[, c("ID", "condition", "IDH_status", "MGMT")])
    rownames(sample_metadata_table_total) <- sample_metadata_table_total$ID
    sample_metadata_table_total <- sample_metadata_table_total[colnames(pseudobulk), , drop = FALSE]
    sample_metadata_table_total$n_myeloid_cells <- sample_cells
    sample_metadata_table_total <- merge(
        sample_metadata_table_total,
        technical_table,
        by    = "ID",
        all.x = TRUE,
        sort  = FALSE
    )
    rownames(sample_metadata_table_total) <- sample_metadata_table_total$ID
    sample_metadata_table_total <- sample_metadata_table_total[colnames(pseudobulk), , drop = FALSE]

    # Keep recurrent unexposed and PD-1 groups and at least 20 cells. Apply the IDH rule for the active analysis.
    formal <- sample_metadata_table_total$condition %in% c(
        "GBM.rec",
        "GBM.PD1"
    ) & sample_metadata_table_total$n_myeloid_cells >= 20
    if (idh_entry == "strict_idhwt") formal <- formal & sample_metadata_table_total$IDH_status == "WT"
    formal_counts <- pseudobulk[, formal, drop = FALSE]
    formal_metadata_table <- sample_metadata_table_total[formal, , drop = FALSE]
    formal_metadata_table$condition <- factor(formal_metadata_table$condition, levels = c("GBM.rec", "GBM.PD1"))

    # Keep filtering thresholds of 5 and 10 and use the original sources for raw and adjusted score standardization.
    y_all <- calcNormFactors(DGEList(formal_counts))
    y <- y_all
    keep_gene <- filterByExpr(y, group = formal_metadata_table$condition, min.count = 5, min.total.count = 10)
    y <- calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
    formal_metadata_table$MGMT <- factor(formal_metadata_table$MGMT)
    formal_metadata_table$sequencing_method <- factor(formal_metadata_table$sequencing_method)
    formal_metadata_table$gel_bead_version <- factor(formal_metadata_table$gel_bead_version)
    formal_metadata_table$acquisition_wave <- factor(formal_metadata_table$acquisition_wave)
    formal_metadata_table$IDH_status <- factor(formal_metadata_table$IDH_status)
    if (idh_entry == "strict_idhwt") {
        design_matrix <- model.matrix(
            ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition,
            formal_metadata_table
        )
    } else {
        design_matrix <- model.matrix(
            ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + IDH_status + condition,
            formal_metadata_table
        )
    }
    if (qr(design_matrix)$rank != ncol(design_matrix)) {
        stop("GSE154795 design is not full rank for ", entry, " / ", idh_entry)
    }

    # Fit robust edgeR with the original adjustment, treatment contrast and BH correction family.
    y <- estimateDisp(y, design_matrix, robust = TRUE)
    fit_model <- glmQLFit(y, design_matrix, robust = TRUE)
    coef_index <- grep("^condition", colnames(design_matrix))
    if (length(coef_index) != 1) stop("Cannot find unique treatment coefficient for ", entry, " / ", idh_entry)
    qlf <- glmQLFTest(fit_model, coef = coef_index)
    deg_table <- topTags(qlf, n = Inf, sort.by = "none")$table
    deg_table$gene <- rownames(deg_table)
    deg_table$FDR <- p.adjust(deg_table$PValue, method = "BH")

    # Rank by F multiplied by the log-fold-change sign and test raw20 enrichment with the original settings.
    ranks <- deg_table$F
    names(ranks) <- deg_table$gene
    ranks <- sort(ranks * sign(deg_table$logFC), decreasing = TRUE)
    fgsea_table <- as.data.table(fgseaMultilevel(
        pathways = list(Miller_raw20 = raw20),
        stats    = ranks,
        minSize  = 5,
        maxSize  = 500,
        eps      = 0,
        nproc    = 1
    ))
    fgsea_table[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    fgsea_table[, cohort := "GSE154795"]
    fgsea_table[, entry := entry]
    fgsea_table[, idh_entry := idh_entry]
    fgsea_table[, contrast := "PD1_minus_UntreatedRecurrent"]

    # Use voom on the same filtered matrix and run camera on the same treatment coefficient.
    v <- voom(y, design_matrix, plot = FALSE)
    raw20_indices <- ids2indices(list(Miller_raw20 = raw20), rownames(v))
    camera_table <- as.data.table(
        camera(v, raw20_indices, design_matrix, contrast = colnames(design_matrix)[coef_index]),
        keep.rownames = "pathway"
    )
    if (!"FDR" %in% colnames(camera_table)) camera_table[, FDR := p.adjust(PValue, method = "BH")]
    camera_table[, cohort := "GSE154795"]
    camera_table[, entry := entry]
    camera_table[, idh_entry := idh_entry]
    camera_table[, contrast := "PD1_minus_UntreatedRecurrent"]

    # Standardize log CPM in the unfiltered y_all and average the raw20 and myeloid gene sets separately.
    logcpm <- cpm(y_all, log = TRUE, prior.count = 1)
    z <- safe_zscore_rows(logcpm)
    raw20_present <- intersect(raw20, rownames(z))
    myeloid_present <- intersect(myeloid_genes, rownames(z))
    score_table <- data.frame(
        sample            = colnames(z),
        cohort            = "GSE154795",
        entry             = entry,
        idh_entry         = idh_entry,
        condition         = formal_metadata_table$condition,
        n_myeloid_cells   = formal_metadata_table$n_myeloid_cells,
        raw20_score       = colMeans(z[raw20_present, , drop = FALSE]),
        myeloid_abundance = colMeans(z[myeloid_present, , drop = FALSE]),
        MGMT              = formal_metadata_table$MGMT,
        sequencing_method = formal_metadata_table$sequencing_method,
        gel_bead_version  = formal_metadata_table$gel_bead_version,
        acquisition_wave  = formal_metadata_table$acquisition_wave,
        IDH_status        = formal_metadata_table$IDH_status,
        stringsAsFactors  = FALSE
    )
    raw_formula <- if (idh_entry == "all_idh") {
        raw20_score ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + IDH_status + condition
    } else {
        raw20_score ~ MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition
    }
    adj_formula <- if (idh_entry == "all_idh") {
        raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version + acquisition_wave +
            IDH_status + condition
    } else {
        raw20_score ~ myeloid_abundance + MGMT + sequencing_method + gel_bead_version + acquisition_wave + condition
    }
    raw_model <- lm(raw_formula, data = score_table)
    adj_model <- lm(adj_formula, data = score_table)
    coef_name <- grep("^condition", names(coef(raw_model)), value = TRUE)

    # Extract effects and 95% intervals from both patient models. Apply BH within the two results for this entry.
    score_model_table_total <- rbindlist(list(
        extract_model_row(
            raw_model,
            coef_name,
            "GSE154795",
            "raw20_score",
            "Raw score",
            "GSE154795",
            "PD1_minus_recurrent",
            sum(score_table$condition == "GBM.PD1"),
            sum(score_table$condition == "GBM.rec"),
            extra = list(entry = entry, idh_entry = idh_entry)
        ),
        extract_model_row(
            adj_model,
            coef_name,
            "GSE154795",
            "raw20_score_myeloid_adjusted",
            "Myeloid-adjusted",
            "GSE154795",
            "PD1_minus_recurrent",
            sum(score_table$condition == "GBM.PD1"),
            sum(score_table$condition == "GBM.rec"),
            extra = list(entry = entry, idh_entry = idh_entry)
        )
    ))
    score_model_table_total[, fdr := p.adjust(p_value, method = "BH")]
    score_model_table_total[, n_label := sprintf("%d / %d", n_positive, n_negative)]
    score_model_table_total[, fdr_label := vapply(fdr, build_p_label, character(1))]

    coverage_table <- data.frame(
        cohort           = "GSE154795",
        entry            = entry,
        idh_entry        = idh_entry,
        signature        = "Miller_raw20",
        n_defined        = length(raw20),
        n_present        = length(raw20_present),
        coverage         = length(raw20_present) / length(raw20),
        present_genes    = paste(raw20_present, collapse = ";"),
        missing_genes    = paste(setdiff(raw20, raw20_present), collapse = ";"),
        stringsAsFactors = FALSE
    )
    result <- list(
        audit        = as.data.table(audit_table),
        sample_meta  = as.data.table(sample_metadata_table_total),
        formal_meta  = as.data.table(formal_metadata_table),
        de           = as.data.table(deg_table),
        fgsea        = fgsea_table,
        camera       = camera_table,
        scores       = as.data.table(score_table),
        score_models = score_model_table_total,
        coverage     = coverage_table,
        ranks        = ranks
    )
    return(result)
}

gse154795_main <- analyze_gse154795_variant(
    gse154795_prep,
    gse154795_technical_table_total,
    entry = "strict",
    idh_entry = "strict_idhwt"
)
gse154795_reviewed <- analyze_gse154795_variant(
    gse154795_prep,
    gse154795_technical_table_total,
    entry = "reviewed_scgbm",
    idh_entry = "strict_idhwt"
)
gse154795_inclusive <- analyze_gse154795_variant(
    gse154795_prep,
    gse154795_technical_table_total,
    entry = "erythroid_ambient_inclusive",
    idh_entry = "strict_idhwt"
)
gse154795_all_idh <- analyze_gse154795_variant(
    gse154795_prep,
    gse154795_technical_table_total,
    entry = "strict",
    idh_entry = "all_idh"
)

step38_dir <- file.path(root_dir, "write", "38_independent_cohort_mg_inflammatory_recalculation")
step38_gsea_file <- file.path(step38_dir, "independent_fixed_program_targeted_gsea.csv")
step38_gene_dir_file <- file.path(step38_dir, "independent_fixed_program_gene_direction.csv")

step38_gsea_table <- fread(step38_gsea_file)
step38_gene_table <- fread(step38_gene_dir_file)
step38_gsea_table <- step38_gsea_table[pathway == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20]
step38_gene_table <- step38_gene_table[signature == "Miller_Microglial_Inflammatory_raw_top20" & threshold == 20]

gse121810_leading <- unlist(strsplit(gse121810$fgsea$leadingEdge[1], ";"))
gse154795_leading <- unlist(strsplit(gse154795_main$fgsea$leadingEdge[1], ";"))
gse174554_leading <- unlist(strsplit(step38_gsea_table[dataset == "GSE174554"]$leadingEdge[1], ";"))
gse274546_leading <- unlist(strsplit(step38_gsea_table[dataset == "GSE274546"]$leadingEdge[1], ";"))

# Draw four timelines from the cohort exposure definitions and read sample sizes from the patient models.
timeline_table_total <- rbindlist(list(
    data.frame(
        cohort           = "GSE121810",
        lane             = 4,
        lane_label       = sprintf("GSE121810 neo-adjuvant (n=%d)", gse121810$score_models[1, n_positive]),
        step             = c(1, 2),
        event            = c("Exposure", "Surgery"),
        event_group      = c("Exposure", "Surgery"),
        label            = c("Pembrolizumab", "Resection"),
        label_offset     = c(0.22, -0.22),
        exposed_preop    = "Yes",
        starts_postop    = "No",
        stringsAsFactors = FALSE
    ),
    data.frame(
        cohort           = "GSE121810",
        lane             = 3,
        lane_label       = sprintf("GSE121810 adjuvant-only (n=%d)", gse121810$score_models[1, n_negative]),
        step             = c(2, 3),
        event            = c("Surgery", "Exposure"),
        event_group      = c("Surgery", "Exposure"),
        label            = c("Resection", "Starts after surgery"),
        label_offset     = c(-0.22, 0.22),
        exposed_preop    = "No",
        starts_postop    = "Yes",
        stringsAsFactors = FALSE
    ),
    data.frame(
        cohort           = "GSE154795",
        lane             = 2,
        lane_label       = sprintf("GSE154795 PD-1 exposed (n=%d)", gse154795_main$score_models[1, n_positive]),
        step             = c(1, 2),
        event            = c("Exposure", "Surgery"),
        event_group      = c("Exposure", "Surgery"),
        label            = c("Anti-PD-1", "Recurrent resection"),
        label_offset     = c(0.22, -0.22),
        exposed_preop    = "Yes",
        starts_postop    = "No",
        stringsAsFactors = FALSE
    ),
    data.frame(
        cohort           = "GSE154795",
        lane             = 1,
        lane_label       = sprintf("GSE154795 untreated recurrent (n=%d)", gse154795_main$score_models[1, n_negative]),
        step             = 2,
        event            = "Surgery",
        event_group      = "Surgery",
        label            = "Recurrent resection",
        label_offset     = -0.22,
        exposed_preop    = "No",
        starts_postop    = "No",
        stringsAsFactors = FALSE
    )
))
timeline_table_total[, label_hjust := fifelse(step == 1, 0, fifelse(step == 3, 1, 0.5))]
timeline_table_total[, label_nudge := fifelse(step == 1, 0.06, fifelse(step == 3, -0.03, 0))]
figure5_source_dir <- file.path(
    root_dir, "write", "41_mg_inflammatory_sci_rebuild", "Figure5", "source_data"
)
write_csv(timeline_table_total, "figure5_treatment_timeline", figure5_source_dir)
pal <- c(
    primary   = "#4C78A8",
    recurrent = "#E45756",
    neutral   = "#777777",
    light     = "#D9D9D9",
    grid      = "#E8E8E8",
    text      = "#222222",
    teal      = "#00A087",
    blue      = "#3C5488",
    mint      = "#91D1C2",
    grayblue  = "#8491B4"
)


# Draw the original exposure and surgery timeline given the event table and title.
timeline_plot <- ggplot(timeline_table_total, aes(x = step, y = lane))
timeline_plot <- timeline_plot + geom_segment(
    aes(x = 1, xend = 3, y = lane, yend = lane),
    color     = pal["light"],
    linewidth = 1.7,
    lineend   = "round"
)
timeline_plot <- timeline_plot + geom_segment(
    data = timeline_table_total[timeline_table_total$exposed_preop == "Yes", ],
    aes(x = 1, xend = 2, y = lane, yend = lane),
    color     = pal["recurrent"],
    linewidth = 1.7,
    lineend   = "round"
)
timeline_plot <- timeline_plot + geom_segment(
    data = timeline_table_total[timeline_table_total$starts_postop == "Yes", ],
    aes(x = 2, xend = 3, y = lane, yend = lane),
    color     = pal["primary"],
    linewidth = 1.7,
    lineend   = "round"
)
timeline_plot <- timeline_plot + geom_point(
    aes(fill = event_group, shape = event),
    size   = 2.7,
    color  = pal["text"],
    stroke = 0.3
)
timeline_plot <- timeline_plot + geom_text(
    aes(x = step + label_nudge, y = lane + label_offset, label = label, hjust = label_hjust),
    size = 2.55,
    lineheight = 0.9
)
timeline_plot <- timeline_plot + scale_shape_manual(values = c(Exposure = 21, Surgery = 24))
timeline_plot <- timeline_plot + scale_fill_manual(values = c("Exposure" = pal["recurrent"], "Surgery" = "white"))
timeline_plot <- timeline_plot + scale_y_continuous(
    NULL,
    breaks = timeline_table_total$lane,
    labels = timeline_table_total$lane_label
)
timeline_plot <- timeline_plot + scale_x_continuous(
    NULL,
    breaks = 1:3,
    labels = c("Pre-op", "Surgery", "Post-op"),
    expand = expansion(mult = c(0.04, 0.10))
)
timeline_plot <- timeline_plot + labs(
    title = "Anti-PD-1 exposure timing differs across the two recurrent GBM cohorts",
    subtitle = paste0(
        "GSE121810 contrasts pre-op exposure with post-op-only start; ",
        "GSE154795 compares exposed with untreated recurrent tumors"
    )
)
timeline_plot <- timeline_plot + theme_panel(7.3)
timeline_plot <- timeline_plot + theme(
    legend.position = "none",
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    plot.title      = element_text(face = "bold", size = 8.2, hjust = 0),
    plot.subtitle   = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
    plot.margin     = margin(8, 12, 6, 6)
)

figure5_dir <- file.path(
    root_dir, "figures", "41_mg_inflammatory_sci_rebuild", "Figure5", "panel_library"
)
save_panel(
    timeline_plot,
    "fig5a_treatment_timeline",
    figure5_dir,
    183,
    58
)

# Use the ranking vector, gene set, cohort and contrast labels to return the GSEA curve.
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
    result_table <- data.frame(
        cohort             = cohort_label,
        contrast           = contrast_label,
        rank               = seq_along(ranks),
        gene               = gene_names,
        rank_stat          = unname(ranks),
        running_enrichment = running,
        hit                = hits,
        hit_gene           = ifelse(hits, gene_names, NA_character_),
        stringsAsFactors   = FALSE
    )
    return(result_table)
}

gse121810_curve_table <- compute_gsea_curve(gse121810$ranks, raw20, "GSE121810", "Neo minus adjuvant-only")
gse154795_curve_table <- compute_gsea_curve(gse154795_main$ranks, raw20, "GSE154795", "PD1 minus untreated recurrent")
write_csv_gz(gse121810_curve_table, "figure5b_gse121810_gsea_curve", figure5_source_dir)
write_csv_gz(gse154795_curve_table, "figure5c_gse154795_gsea_curve", figure5_source_dir)

gsea_stats_table_total <- rbindlist(list(
    data.frame(
        cohort           = "GSE121810",
        comparison_label = "Pre-op exposed\nvs post-op only",
        NES              = gse121810$fgsea$NES[1],
        fgsea_p          = gse121810$fgsea$pval[1],
        camera_direction = gse121810$camera$Direction[1],
        camera_p         = gse121810$camera$PValue[1],
        n_positive       = gse121810$score_models[1, n_positive],
        n_negative       = gse121810$score_models[1, n_negative],
        n_present        = gse121810$coverage$n_present[1],
        stringsAsFactors = FALSE
    ),
    data.frame(
        cohort           = "GSE154795",
        comparison_label = "PD-1 exposed\nvs untreated recurrent",
        NES              = gse154795_main$fgsea$NES[1],
        fgsea_p          = gse154795_main$fgsea$pval[1],
        camera_direction = gse154795_main$camera$Direction[1],
        camera_p         = gse154795_main$camera$PValue[1],
        n_positive       = gse154795_main$score_models[1, n_positive],
        n_negative       = gse154795_main$score_models[1, n_negative],
        n_present        = gse154795_main$coverage$n_present[1],
        stringsAsFactors = FALSE
    )
))
write_csv(gsea_stats_table_total, "figure5_gsea_side_stats", figure5_source_dir)

# Combine the curve, statistics and labels into the original plot layout without writing files.
build_gsea_panel <- function(curve_table, stats_table, stem_label, entry_label) {
    running_table <- as.data.table(copy(curve_table))
    hit_table <- running_table[hit == TRUE, .(rank)]
    rank_table <- data.frame(rank = running_table$rank, rank_stat = running_table$rank_stat)
    stats_lines <- c(
        stem_label,
        entry_label,
        sprintf("%d/20 genes measured", stats_table$n_present),
        sprintf("%s", stats_table$comparison_label),
        sprintf("NES %.2f   P %s", stats_table$NES, build_p_label(stats_table$fgsea_p)),
        sprintf("camera %s   P %s", stats_table$camera_direction, build_p_label(stats_table$camera_p)),
        sprintf("n = %d vs %d", stats_table$n_positive, stats_table$n_negative)
    )

    # Plot the running enrichment curve and its area across the full ranking. The zero line marks direction.
    top_plot <- ggplot(running_table, aes(x = rank, y = running_enrichment))
    top_plot <- top_plot + geom_hline(yintercept = 0, color = pal["light"], linewidth = 0.35)
    top_plot <- top_plot + geom_line(color = pal["recurrent"], linewidth = 0.75)
    top_plot <- top_plot + geom_area(fill = pal["recurrent"], alpha = 0.12)
    top_plot <- top_plot + labs(x = NULL, y = "Running enrichment")
    top_plot <- top_plot + theme_panel(7.2)
    top_plot <- top_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    # Mark each raw20 hit at its position in the gene ranking.
    hit_plot <- ggplot(hit_table, aes(x = rank, y = 1))
    hit_plot <- hit_plot + geom_segment(aes(xend = rank, y = 0, yend = 1), color = pal["text"], linewidth = 0.22)
    hit_plot <- hit_plot + scale_y_continuous(NULL, breaks = NULL)
    hit_plot <- hit_plot + labs(x = NULL, y = NULL)
    hit_plot <- hit_plot + theme_void(base_family = "sans")

    # Plot the signed ranking statistic and use the original colors for the two treatment-associated directions.
    rank_plot <- ggplot(rank_table, aes(x = rank, y = rank_stat))
    rank_plot <- rank_plot + geom_hline(yintercept = 0, color = pal["light"], linewidth = 0.3)
    rank_plot <- rank_plot + geom_area(
        data = rank_table[rank_table$rank_stat >= 0,
        ],
        fill  = pal["recurrent"],
        alpha = 0.85
    )
    rank_plot <- rank_plot + geom_area(
        data = rank_table[rank_table$rank_stat < 0,
        ],
        fill  = pal["primary"],
        alpha = 0.85
    )
    rank_plot <- rank_plot + labs(x = "Ranked genes", y = NULL)
    rank_plot <- rank_plot + theme_panel(7.2)
    rank_plot <- rank_plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    left <- top_plot / hit_plot / rank_plot
    left <- left + plot_layout(heights = c(5.0, 0.65, 1.4))

    # Annotate coverage, NES, both tests and group sample sizes from the calculation.
    text_table <- data.frame(
        x     = 0,
        y     = rev(seq_along(stats_lines)),
        label = stats_lines,
        face  = c("bold", rep("plain", length(stats_lines) - 1)),
        stringsAsFactors = FALSE
    )
    right <- ggplot(text_table, aes(x = x, y = y, label = label))
    right <- right + geom_text(
        data = text_table[1,
            ,
            drop = FALSE],
        aes(fontface = face),
        hjust = 0,
        vjust = 1,
        size  = 3.1
    )
    right <- right + geom_text(
        data = text_table[-1,
            ,
            drop = FALSE],
        aes(fontface = face),
        hjust = 0,
        vjust = 1,
        size  = 2.75
    )
    right <- right + xlim(0, 1)
    right <- right + ylim(0.4, length(stats_lines) + 0.6)
    right <- right + theme_void(base_family = "sans")

    result_plot <- (left | right)
    result_plot <- result_plot + plot_layout(widths = c(4.0, 2.8))
    return(result_plot)
}

save_panel(
    build_gsea_panel(
        gse121810_curve_table,
        gsea_stats_table_total[cohort == "GSE121810",
        ],
        "GSE121810",
        "Randomized timing cohort"
    ),
    "fig5b_gse121810_raw20_gsea",
    figure5_dir,
    135,
    82
)
save_panel(
    build_gsea_panel(
        gse154795_curve_table,
        gsea_stats_table_total[cohort == "GSE154795",
        ],
        "GSE154795",
        "Strict IDH-wt auto-myeloid entry"
    ),
    "fig5c_gse154795_raw20_gsea",
    figure5_dir,
    135,
    82
)

# Assemble the raw and myeloid-adjusted model effects for the two formal cohorts used in the main figure.
main_effect_table_total <- rbindlist(list(
    gse121810$score_models[, .(
        cohort,
        cohort_label,
        metric,
        metric_label,
        effect,
        conf_low,
        conf_high,
        fdr,
        n_positive,
        n_negative
    )],
    gse154795_main$score_models[, .(
        cohort,
        cohort_label,
        metric,
        metric_label,
        effect,
        conf_low,
        conf_high,
        fdr,
        n_positive,
        n_negative
    )]
))
main_effect_table_total[, row_label := c(
    "GSE121810 raw score", "GSE121810 myeloid-adjusted",
    "GSE154795 raw score", "GSE154795 myeloid-adjusted"
)]
main_effect_table_total[, n_label := sprintf("%d / %d", n_positive, n_negative)]
main_effect_table_total[, fdr_label := vapply(fdr, build_p_label, character(1))]
write_csv(main_effect_table_total, "figure5d_patient_score_models", figure5_source_dir)

main_score_table_total <- rbindlist(list(
    gse121810$scores[, .(sample, cohort, group, raw20_score, myeloid_abundance)],
    gse154795_main$scores[, .(sample, cohort, condition, raw20_score, myeloid_abundance, entry, idh_entry)]
), fill = TRUE)
write_csv(main_score_table_total, "figure5d_patient_score_samples", figure5_source_dir)

# Plot model effects and intervals from the effect table, axis labels, colors and facet settings.
build_forest_plot <- function(
    input_table,
    x_label,
    point_colors,
    facet_var     = NULL,
    title_text    = NULL,
    subtitle_text = NULL
) {
    plot_table <- copy(input_table)
    x_rng <- range(c(plot_table$conf_low, plot_table$conf_high), finite = TRUE)
    span <- diff(x_rng)
    if (!is.finite(span) || span <= 0) span <- 0.2
    show_n <- any(nchar(plot_table$n_label) > 0)
    text_x1 <- x_rng[2] + span * 0.16
    text_x2 <- if (show_n) x_rng[2] + span * 0.46 else x_rng[2] + span * 0.24
    xmax <- x_rng[2] + span * 0.95

    # Show the original estimate and 95% interval per model. Point shapes distinguish raw and adjusted scores.
    p <- ggplot(plot_table, aes(x = effect, y = factor(row_label, levels = rev(unique(row_label)))))
    p <- p + geom_vline(xintercept = 0, color = pal["light"], linewidth = 0.35)
    p <- p + geom_errorbarh(
        aes(xmin = conf_low, xmax = conf_high, color = metric_label),
        height    = 0.18,
        linewidth = 0.55
    )
    p <- p + geom_point(aes(color = metric_label, shape = metric_label), size = 2.1, stroke = 0.25)
    p <- p + geom_text(aes(x = text_x2, label = fdr_label), size = 2.5, hjust = 0)
    p <- p + scale_color_manual(values = setNames(unname(point_colors), names(point_colors)))
    p <- p + scale_shape_manual(values = c("Raw score" = 16, "Myeloid-adjusted" = 17))
    p <- p + scale_x_continuous(
        x_label,
        expand = expansion(mult = c(0.08, 0.18)),
        limits = c(x_rng[1] - span * 0.08, xmax)
    )
    p <- p + coord_cartesian(clip = "off")
    p <- p + labs(y = NULL, title = title_text, subtitle = subtitle_text)
    p <- p + theme_panel(7.2)
    p <- p + theme(
        legend.position = "bottom",
        plot.title      = element_text(face = "bold", size = 8.2, hjust = 0),
        plot.subtitle   = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
        plot.margin     = margin(12, 42, 6, 6)
    )

    # Include sample sizes where available, always show FDR and use the requested facets.
    if (show_n) {
        p <- p + geom_text(aes(x = text_x1, label = n_label), size = 2.5, hjust = 0)
        p <- p + annotate(
            "text",
            x        = text_x1,
            y        = length(unique(plot_table$row_label)) + 0.24,
            label    = "n",
            hjust    = 0,
            size     = 2.5,
            fontface = "bold"
        )
    }
    p <- p + annotate(
        "text",
        x        = text_x2,
        y        = length(unique(plot_table$row_label)) + 0.24,
        label    = "FDR",
        hjust    = 0,
        size     = 2.5,
        fontface = "bold"
    )
    if (!is.null(facet_var)) {
        p <- p + facet_grid(as.formula(paste(facet_var, "~ .")), scales = "free_y", space = "free_y")
    }
    return(p)
}

save_panel(
    build_forest_plot(
        main_effect_table_total,
        x_label       = "Treatment-associated score shift",
        point_colors  = c("Raw score" = pal["primary"], "Myeloid-adjusted" = pal["recurrent"]),
        title_text    = "Patient-level raw and myeloid-adjusted effects",
        subtitle_text = "Linear models with 95% CI; GSE154795 raw score keeps technical covariates"
    ),
    "fig5d_patient_score_forest",
    figure5_dir,
    140,
    88
)

# Retain each cohort's log fold change scale and order genes by leading-edge count and mean log fold change.
recurrence_gene_table <- step38_gene_table[, .(dataset, gene, logFC)]
gse121810_gene_table <- gse121810$de[, .(dataset = "GSE121810", gene, logFC)]
gse154795_gene_table <- gse154795_main$de[, .(dataset = "GSE154795", gene, logFC)]
matrix_table_total <- rbindlist(list(
    recurrence_gene_table[dataset %in% c("GSE174554", "GSE274546")],
    gse121810_gene_table,
    gse154795_gene_table
), fill = TRUE)
matrix_table_total <- matrix_table_total[gene %in% raw20]
matrix_table_total[, leading_edge := FALSE]
matrix_table_total[dataset == "GSE174554" & gene %in% gse174554_leading, leading_edge := TRUE]
matrix_table_total[dataset == "GSE274546" & gene %in% gse274546_leading, leading_edge := TRUE]
matrix_table_total[dataset == "GSE121810" & gene %in% gse121810_leading, leading_edge := TRUE]
matrix_table_total[dataset == "GSE154795" & gene %in% gse154795_leading, leading_edge := TRUE]
gene_order <- matrix_table_total[, .(
    leading_count = sum(leading_edge),
    mean_logFC    = mean(logFC, na.rm = TRUE)
), by = gene][order(-leading_count, -mean_logFC)]$gene
matrix_table_total[, gene := factor(gene, levels = gene_order)]
matrix_table_total[, cohort_label := factor(dataset, levels = c("GSE174554", "GSE274546", "GSE121810", "GSE154795"))]
write_csv(matrix_table_total, "figure5e_four_cohort_leading_edge_logfc", figure5_source_dir)

# Plot the four-cohort gene effects with the original log fold change scale and leading-edge dots.
plot_table <- copy(matrix_table_total)
plot_table$cohort_label <- factor(
    plot_table$cohort_label,
    levels = c("GSE174554", "GSE274546", "GSE121810", "GSE154795"),
    labels = c(
        "GSE174554\nRecurrent logFC",
        "GSE274546\nRecurrent logFC",
        "GSE121810\nExposure logFC",
        "GSE154795\nExposure logFC"
    )
)
plot_table$gene <- factor(plot_table$gene, levels = rev(unique(plot_table$gene)))
leading_table <- plot_table[leading_edge == TRUE, ]
heatmap_plot <- ggplot(plot_table, aes(x = cohort_label, y = gene))
heatmap_plot <- heatmap_plot + geom_tile(aes(fill = logFC), color = "white", linewidth = 0.32)
heatmap_plot <- heatmap_plot + geom_point(
    data = leading_table,
    aes(shape = "Leading-edge member"),
    size   = 1.7,
    stroke = 0.15,
    color  = "black",
    inherit.aes = TRUE
)
heatmap_plot <- heatmap_plot + scale_fill_gradient2(
    low  = pal["primary"],
    mid  = "white",
    high = pal["recurrent"],
    midpoint = 0,
    name = "logFC"
)
heatmap_plot <- heatmap_plot + scale_shape_manual(values = c("Leading-edge member" = 16), name = NULL)
heatmap_plot <- heatmap_plot + guides(
    fill  = guide_colorbar(order = 1, title.position = "top"),
    shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.2))
)
heatmap_plot <- heatmap_plot + labs(
    x        = NULL,
    y        = NULL,
    title    = "Leading-edge genes keep the same direction across recurrence and exposure cohorts",
    subtitle = "First two columns are paired recurrence logFC; last two columns are anti-PD-1 exposure logFC"
)
heatmap_plot <- heatmap_plot + theme_panel(7)
heatmap_plot <- heatmap_plot + theme(
    axis.text.x     = element_text(angle = 25, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.box      = "vertical",
    plot.title      = element_text(face = "bold", size = 8.0, hjust = 0),
    plot.subtitle   = element_text(size = 6.6, color = pal["neutral"], hjust = 0)
)

save_panel(
    heatmap_plot,
    "fig5e_four_cohort_leading_edge_heatmap",
    figure5_dir,
    118,
    115
)

# Generate the original supplementary panels for myeloid definition, IDH definition and score effects.
entry_fgsea_table_total <- rbindlist(list(
    data.frame(row_label = "Auto strict IDH-wt", NES = gse154795_main$fgsea$NES[1], p = gse154795_main$fgsea$pval[1]),
    data.frame(
        row_label = "Reviewed IDH-wt",
        NES = gse154795_reviewed$fgsea$NES[1],
        p   = gse154795_reviewed$fgsea$pval[1]
    ),
    data.frame(
        row_label = "Inclusive IDH-wt",
        NES = gse154795_inclusive$fgsea$NES[1],
        p   = gse154795_inclusive$fgsea$pval[1]
    )
))
entry_fgsea_table_total[, p_label := vapply(p, build_p_label, character(1))]
s7_source_dir <- file.path(
    root_dir, "write", "41_mg_inflammatory_sci_rebuild", "SupplementaryFigure7", "source_data"
)
write_csv(entry_fgsea_table_total, "supp7a_entry_fgsea_sensitivity", s7_source_dir)

# Plot the NES sensitivity table with the original axes and annotation locations.
build_nes_plot <- function(input_table, x_label, title_text = NULL, subtitle_text = NULL) {
    plot_table <- copy(input_table)
    plot_table$row_label <- factor(plot_table$row_label, levels = rev(plot_table$row_label))
    x_rng <- range(plot_table$NES, finite = TRUE)
    span <- diff(x_rng)
    if (!is.finite(span) || span <= 0) span <- 0.4
    text_x <- x_rng[2] + span * 0.35
    xmax <- x_rng[2] + span * 0.95
    result_plot <- ggplot(plot_table, aes(x = NES, y = row_label))
    result_plot <- result_plot + geom_vline(xintercept = 0, color = pal["light"], linewidth = 0.35)
    result_plot <- result_plot + geom_segment(
        aes(x = 0, xend = NES, y = row_label, yend = row_label),
        color = pal["grayblue"],
        linewidth = 0.55
    )
    result_plot <- result_plot + geom_point(size = 2.2, color = pal["recurrent"])
    result_plot <- result_plot + geom_text(aes(x = text_x, label = p_label), hjust = 0, size = 2.55)
    result_plot <- result_plot + annotate(
        "text",
        x        = text_x,
        y        = length(unique(plot_table$row_label)) + 0.36,
        label    = "P",
        hjust    = 0,
        size     = 2.45,
        fontface = "bold"
    )
    result_plot <- result_plot + scale_y_discrete(expand = expansion(add = c(0.45, 0.82)))
    result_plot <- result_plot + scale_x_continuous(
        x_label,
        expand = expansion(mult = c(0.08, 0.18)),
        limits = c(min(0, x_rng[1] - span * 0.08), xmax)
    )
    result_plot <- result_plot + coord_cartesian(clip = "off")
    result_plot <- result_plot + labs(y = NULL, title = title_text, subtitle = subtitle_text)
    result_plot <- result_plot + theme_panel(7.2)
    result_plot <- result_plot + theme(
        plot.title    = element_text(face = "bold", size = 8.2, hjust = 0),
        plot.subtitle = element_text(size = 6.8, color = pal["neutral"], hjust = 0),
        plot.margin   = margin(12, 36, 6, 6)
    )
    return(result_plot)
}

s7_dir <- file.path(
    root_dir, "figures", "41_mg_inflammatory_sci_rebuild", "SupplementaryFigure7", "panel_library"
)
save_panel(
    build_nes_plot(
        entry_fgsea_table_total,
        x_label       = "fgsea NES",
        title_text    = "GSE154795 entry sensitivity",
        subtitle_text = "Only fgsea shown here; camera stays in the separate algorithm table"
    ),
    "supp7a_entry_fgsea_sensitivity",
    s7_dir,
    125,
    72
)

idh_fgsea_table_total <- rbindlist(list(
    data.frame(row_label = "Strict IDH-wt", NES = gse154795_main$fgsea$NES[1], p = gse154795_main$fgsea$pval[1]),
    data.frame(row_label = "All IDH", NES = gse154795_all_idh$fgsea$NES[1], p = gse154795_all_idh$fgsea$pval[1])
))
idh_fgsea_table_total[, p_label := vapply(p, build_p_label, character(1))]
write_csv(idh_fgsea_table_total, "supp7b_idh_fgsea_sensitivity", s7_source_dir)
save_panel(
    build_nes_plot(
        idh_fgsea_table_total,
        x_label       = "fgsea NES",
        title_text    = "GSE154795 IDH sensitivity",
        subtitle_text = "All-IDH inclusion weakens the coordinated shift"
    ),
    "supp7b_idh_fgsea_sensitivity",
    s7_dir,
    125,
    62
)

score_sensitivity_table_total <- rbindlist(list(
    gse154795_main$score_models[, .(scenario = "Auto strict IDH-wt", metric_label, effect, conf_low, conf_high, fdr)],
    gse154795_reviewed$score_models[, .(scenario = "Reviewed IDH-wt", metric_label, effect, conf_low, conf_high, fdr)],
    gse154795_inclusive$score_models[, .(
        scenario = "Inclusive IDH-wt",
        metric_label,
        effect,
        conf_low,
        conf_high,
        fdr
    )],
    gse154795_all_idh$score_models[, .(scenario = "All IDH", metric_label, effect, conf_low, conf_high, fdr)]
))
score_sensitivity_table_total[, row_label := scenario]
score_sensitivity_table_total[, n_label := ""]
score_sensitivity_table_total[, fdr_label := vapply(fdr, build_p_label, character(1))]
write_csv(score_sensitivity_table_total, "supp7c_score_sensitivity_models", s7_source_dir)
save_panel(
    build_forest_plot(
        score_sensitivity_table_total,
        x_label       = "Treatment-associated score shift",
        point_colors  = c("Raw score" = pal["primary"], "Myeloid-adjusted" = pal["recurrent"]),
        facet_var     = "metric_label",
        title_text    = "GSE154795 score sensitivity",
        subtitle_text = "Raw and myeloid-adjusted effects stay on their own score scale"
    ),
    "supp7c_score_sensitivity",
    s7_dir,
    145,
    106
)

# Write the original fgsea and camera supplementary tables for downstream use.
algorithm_table_total <- rbindlist(list(
    data.frame(
        cohort_label = "GSE121810",
        analysis_label = "Main exposure comparison",
        n_label = sprintf("%d / %d", gse121810$score_models[1, n_positive], gse121810$score_models[1, n_negative]),
        fgsea_label = sprintf("NES %.2f; P %s", gse121810$fgsea$NES[1], build_p_label(gse121810$fgsea$pval[1])),
        camera_label = sprintf("%s; P %s", gse121810$camera$Direction[1], build_p_label(gse121810$camera$PValue[1]))
    ),
    data.frame(
        cohort_label   = "GSE154795",
        analysis_label = "Auto strict IDH-wt",
        n_label        = sprintf(
            "%d / %d",
            gse154795_main$score_models[1,
                n_positive],
            gse154795_main$score_models[1,
                n_negative]
        ),
        fgsea_label = sprintf(
            "NES %.2f; P %s",
            gse154795_main$fgsea$NES[1],
            build_p_label(gse154795_main$fgsea$pval[1])
        ),
        camera_label = sprintf(
            "%s; P %s",
            gse154795_main$camera$Direction[1],
            build_p_label(gse154795_main$camera$PValue[1])
        )
    ),
    data.frame(
        cohort_label   = "GSE154795",
        analysis_label = "Reviewed IDH-wt",
        n_label        = sprintf(
            "%d / %d",
            gse154795_reviewed$score_models[1,
                n_positive],
            gse154795_reviewed$score_models[1,
                n_negative]
        ),
        fgsea_label = sprintf(
            "NES %.2f; P %s",
            gse154795_reviewed$fgsea$NES[1],
            build_p_label(gse154795_reviewed$fgsea$pval[1])
        ),
        camera_label = sprintf(
            "%s; P %s",
            gse154795_reviewed$camera$Direction[1],
            build_p_label(gse154795_reviewed$camera$PValue[1])
        )
    ),
    data.frame(
        cohort_label   = "GSE154795",
        analysis_label = "Inclusive IDH-wt",
        n_label        = sprintf(
            "%d / %d",
            gse154795_inclusive$score_models[1,
                n_positive],
            gse154795_inclusive$score_models[1,
                n_negative]
        ),
        fgsea_label = sprintf(
            "NES %.2f; P %s",
            gse154795_inclusive$fgsea$NES[1],
            build_p_label(gse154795_inclusive$fgsea$pval[1])
        ),
        camera_label = sprintf(
            "%s; P %s",
            gse154795_inclusive$camera$Direction[1],
            build_p_label(gse154795_inclusive$camera$PValue[1])
        )
    ),
    data.frame(
        cohort_label   = "GSE154795",
        analysis_label = "All IDH",
        n_label        = sprintf(
            "%d / %d",
            gse154795_all_idh$score_models[1,
                n_positive],
            gse154795_all_idh$score_models[1,
                n_negative]
        ),
        fgsea_label = sprintf(
            "NES %.2f; P %s",
            gse154795_all_idh$fgsea$NES[1],
            build_p_label(gse154795_all_idh$fgsea$pval[1])
        ),
        camera_label = sprintf(
            "%s; P %s",
            gse154795_all_idh$camera$Direction[1],
            build_p_label(gse154795_all_idh$camera$PValue[1])
        )
    )
))
write_csv(algorithm_table_total, "supp7d_algorithm_summary_table", s7_source_dir)

figure5_write_dir <- dirname(figure5_source_dir)
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

s7_write_dir <- dirname(s7_source_dir)
write_csv(gse154795_reviewed$fgsea, "GSE154795_raw20_fgsea_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_reviewed$camera, "GSE154795_raw20_camera_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_reviewed$score_models, "GSE154795_raw20_score_models_reviewed_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$fgsea, "GSE154795_raw20_fgsea_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$camera, "GSE154795_raw20_camera_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_inclusive$score_models, "GSE154795_raw20_score_models_inclusive_idhwt", s7_write_dir)
write_csv(gse154795_all_idh$fgsea, "GSE154795_raw20_fgsea_auto_all_idh", s7_write_dir)
write_csv(gse154795_all_idh$camera, "GSE154795_raw20_camera_auto_all_idh", s7_write_dir)
write_csv(gse154795_all_idh$score_models, "GSE154795_raw20_score_models_auto_all_idh", s7_write_dir)

# Record the existing filenames in panel save order and check each PDF and PNG.
figure5_manifest_table <- data.frame(
    panel = c("Figure5A", "Figure5B", "Figure5C", "Figure5D", "Figure5E"),
    stem  = c(
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
figure5_manifest_table$png <- file.path(figure5_dir, paste0(figure5_manifest_table$stem, ".png"))
figure5_manifest_table$pdf <- file.path(figure5_dir, paste0(figure5_manifest_table$stem, ".pdf"))
write_csv(figure5_manifest_table, "Figure5_panel_manifest", figure5_write_dir)

s7_manifest_table <- data.frame(
    panel = c("SupplementaryFigure7A", "SupplementaryFigure7B", "SupplementaryFigure7C"),
    stem  = c(
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
s7_manifest_table$png <- file.path(s7_dir, paste0(s7_manifest_table$stem, ".png"))
s7_manifest_table$pdf <- file.path(s7_dir, paste0(s7_manifest_table$stem, ".pdf"))
write_csv(s7_manifest_table, "SupplementaryFigure7_panel_manifest", s7_write_dir)

cat("panels=", nrow(figure5_manifest_table) + nrow(s7_manifest_table), "\n", sep = "")

# Verify that the manifests refer to nonempty figure files.
plot_paths <- c(figure5_manifest_table$pdf, figure5_manifest_table$png, s7_manifest_table$pdf, s7_manifest_table$png)
stopifnot(all(file.exists(plot_paths)), all(file.info(plot_paths)$size > 0))
