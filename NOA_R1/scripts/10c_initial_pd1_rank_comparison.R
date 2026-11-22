#!/usr/bin/env Rscript

# Purpose  Compare the frozen historical F-based enrichment with a common signed-square-root F ranking.
# Inputs   step5_analysis_plan.md, target_top20 in step2_locked_gene_sets.csv or
#          step2_extended_specificity/locked_gene_sets_used.csv (set_id and gene), and
#          step2_pd1_camera/{GSE121810,GSE154795}_tested_genes.csv (gene).
#          NOA_LEGACY_ROOT/write/41_mg_inflammatory_sci_rebuild/Figure5/:
#          GSE121810_raw20_{edger,fgsea}.csv and GSE154795_raw20_{edger,fgsea}_auto_strict_idhwt.csv.
#          DE fields are gene, F, logFC and PValue. Historical enrichment supplies size, ES, NES and pval.
# Outputs  step5_rank_unification/{sqrtF_fgsea.csv,historical_F_fgsea.csv,ranked_genes.csv,
#          negative_F_genes.csv,checks.csv,sessionInfo.txt}.
# Run      Rscript NOA_R1/scripts/10c_initial_pd1_rank_comparison.R <scratch_revision_dir>
# The historical entry retains fgsea's one-process dispatch and its original output directory.

suppressPackageStartupMessages({
    library(fgsea)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the revision directory.")
revision <- normalizePath(args[1], mustWork = TRUE)
output_dir <- file.path(revision, "step5_rank_unification")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir)

# Keep the frozen check labels, including the protocol amendment for negative F values. No return value.
check_tables <- list()
record_check <- function(check_name, passed) {
    check_tables[[length(check_tables) + 1L]] <<- data.frame(check = check_name, passed = isTRUE(passed))

    if (!isTRUE(passed)) {
        checks_table_total <- do.call(rbind, check_tables)
        write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
        stop(check_name)
    }
}

plan <- readLines(file.path(revision, "step5_analysis_plan.md"), warn = FALSE)
record_check(
    "plan_and_F_addendum",
    any(grepl("2026-10-20 00:55:04", plan, fixed = TRUE)) &&
        any(grepl("2026-10-20 00:57:50", plan, fixed = TRUE))
)

# The packaged gene-set table retains the original target genes in author-ranked order.
definition_path <- file.path(revision, "step2_locked_gene_sets.csv")
if (!file.exists(definition_path)) {
    definition_path <- file.path(revision, "step2_extended_specificity", "locked_gene_sets_used.csv")
}
definition_table <- read.csv(definition_path)
target_genes <- definition_table$gene[definition_table$set_id == "target_top20"]
record_check("target_20", length(target_genes) == 20L)

# Retain source-file names and the original dataset order for concatenation.
result_tables <- list()
historical_tables <- list()
negative_tables <- list()
rank_tables <- list()
project <- normalizePath(Sys.getenv("NOA_LEGACY_ROOT", unset = dirname(revision)), mustWork = TRUE)
legacy_figure_dir <- file.path(project, "write", "41_mg_inflammatory_sci_rebuild", "Figure5")

for (dataset in c("GSE121810", "GSE154795")) {
    if (dataset == "GSE121810") {
        de_table <- read.csv(file.path(legacy_figure_dir, "GSE121810_raw20_edger.csv"))
        historical_table <- read.csv(file.path(legacy_figure_dir, "GSE121810_raw20_fgsea.csv"))
    } else {
        de_table <- read.csv(file.path(legacy_figure_dir, "GSE154795_raw20_edger_auto_strict_idhwt.csv"))
        historical_table <- read.csv(file.path(legacy_figure_dir, "GSE154795_raw20_fgsea_auto_strict_idhwt.csv"))
    }

    # Match the tested genes from camera and confirm the 19 measured target genes.
    background_table <- read.csv(file.path(revision, "step2_pd1_camera", paste0(dataset, "_tested_genes.csv")))
    background <- background_table$gene
    record_check(
        paste(dataset, "same_background"),
        nrow(de_table) == length(background) && setequal(de_table$gene, background) && !anyDuplicated(de_table$gene)
    )
    record_check(
        paste(dataset, "valid_F_and_target"),
        all(is.finite(de_table$F)) && all(is.finite(de_table$logFC)) && sum(target_genes %in% de_table$gene) == 19L
    )

    # Keep the original F in its own column. Clamp negative values only for the new ranking.
    negative_table <- de_table[de_table$F < 0, c("gene", "logFC", "F", "PValue")]
    if (nrow(negative_table)) {
        negative_table$dataset <- dataset
        negative_tables[[dataset]] <- negative_table
    }
    rank_table <- data.frame(
        dataset    = dataset,
        gene       = de_table$gene,
        logFC      = de_table$logFC,
        original_F = de_table$F,
        rank_sqrtF = sign(de_table$logFC) * sqrt(pmax(de_table$F, 0)),
        negative_F = de_table$F < 0,
        target     = de_table$gene %in% target_genes
    )
    rank_table <- rank_table[order(-rank_table$rank_sqrtF, rank_table$gene), ]
    rank_tables[[dataset]] <- rank_table
    ranked_genes <- setNames(rank_table$rank_sqrtF, rank_table$gene)

    # Use the fixed cohort-specific streams and keep ties in the prescribed gene order.
    seed <- if (dataset == "GSE121810") 20261021L else 20261022L
    set.seed(seed)
    fgsea_table <- fgseaMultilevel(
        pathways = list(Miller_raw20 = target_genes),
        stats    = ranked_genes,
        minSize  = 5,
        maxSize  = 500,
        eps      = 0,
        nproc    = 1
    )
    setDF(fgsea_table)
    record_check(
        paste(dataset, "one_fgsea"),
        nrow(fgsea_table) == 1L && is.finite(fgsea_table$pval) && fgsea_table$size == 19L
    )

    # Store one row per cohort and ranking method. ES, NES and F retain their frozen names.
    result_tables[[dataset]] <- data.frame(
        dataset  = dataset,
        method   = "signed_sqrt_F",
        seed     = seed,
        n_tested = nrow(rank_table),
        negative_F_clamped = nrow(negative_table),
        size     = fgsea_table$size,
        ES       = fgsea_table$ES,
        NES      = fgsea_table$NES,
        p_value  = fgsea_table$pval,
        log2err  = fgsea_table$log2err,
        leading_edge = paste(unlist(fgsea_table$leadingEdge, use.names = FALSE), collapse = ";")
    )
    historical_tables[[dataset]] <- data.frame(
        dataset  = dataset,
        method   = "historical_signed_F",
        seed     = NA_integer_,
        n_tested = nrow(de_table),
        negative_F_clamped = NA_integer_,
        size     = historical_table$size,
        ES       = historical_table$ES,
        NES      = historical_table$NES,
        p_value  = historical_table$pval,
        log2err  = historical_table$log2err,
        leading_edge = historical_table$leadingEdge
    )
}

# Save in the original file order. The negative-F list is empty when neither cohort has a negative F.
result_table_total <- do.call(rbind, result_tables)
historical_table_total <- do.call(rbind, historical_tables)
rank_table_total <- do.call(rbind, rank_tables)
negative_table_total <- do.call(rbind, negative_tables)
checks_table_total <- do.call(rbind, check_tables)
write.csv(result_table_total, file.path(output_dir, "sqrtF_fgsea.csv"), row.names = FALSE, na = "")
write.csv(historical_table_total, file.path(output_dir, "historical_F_fgsea.csv"), row.names = FALSE, na = "")
write.csv(rank_table_total, file.path(output_dir, "ranked_genes.csv"), row.names = FALSE, na = "")
write.csv(negative_table_total, file.path(output_dir, "negative_F_genes.csv"), row.names = FALSE, na = "")
write.csv(checks_table_total, file.path(output_dir, "checks.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

cat("RANK_UNIFICATION ", nrow(checks_table_total), " failed ", sum(!checks_table_total$passed), "\n", sep = "")
