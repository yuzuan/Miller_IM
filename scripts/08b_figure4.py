#!/usr/bin/env python3
# Purpose  Draw Figure 4 scores, protein changes and rankings from the original PDC000514 pairs.
# Inputs   Step 41 Figure4/source_data patient differences, protein statistics and ranks, Step 39 peptide sensitivity,
#          and PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv under MILLER_IM_DATA_ROOT.
# Outputs  Source tables and Figure4_panel_manifest.csv under
#          write/54_figure4_proteome_rebuild/Figure4. The figures tree contains three panels and the proteome preview.
# Run      python scripts/08b_figure4.py from the repository root.
# Keep the existing Python plotting and data libraries to preserve the formulas and PDF layout.
# The source_data fields, filenames and write/figures directories are public interfaces.

# Each panel below has a PDF, PNG and matching source CSV in source_data.
#   Figure4A_pdc_patient_delta_waterfall
#   Figure4B_raw20_protein_patient_distributions
#   Figure4C_full_proteome_rank_landscape
# The preview also has a grayscale PNG and a layout CSV using the same filename stem.

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
apply_publication_style = figure_style.apply_publication_style
clean_axis = figure_style.clean_axis
image_font = figure_style.image_font
new_figure = figure_style.new_figure
panel_record = figure_style.panel_record
paste_contained = figure_style.paste_contained
save_panel = figure_style.save_panel

# Apply the publication style and create the output directories before reading existing source data.
apply_publication_style()
root = Path(__file__).resolve().parents[1]
write_root = root / "write/54_figure4_proteome_rebuild/Figure4"
source_out = write_root / "source_data"
source_out.mkdir(parents=True, exist_ok=True)
figure_out = root / "figures/54_figure4_proteome_rebuild/Figure4"
figure_out.mkdir(parents=True, exist_ok=True)

# Each row is one paired patient. Check 105 pairs, 72 increases and the mean against the frozen source.
step41_source = root / "write/41_mg_inflammatory_sci_rebuild/Figure4/source_data"

patient_path = step41_source / "Figure4A_patient_delta_source.csv"
program_delta_table = pd.read_csv(patient_path)
if len(program_delta_table) != 105:
    raise ValueError(f"The primary PDC analysis requires 105 patient pairs. Found {len(program_delta_table)}.")
if int((program_delta_table["delta"] > 0).sum()) != 72:
    raise ValueError("The primary PDC analysis requires higher scores at recurrence in 72/105 pairs.")
if not np.isclose(program_delta_table["delta"].mean(), 0.252534746498684, atol=1e-12):
    raise ValueError("The mean patient-level PDC delta has changed.")

# Each peptide summary row is one analysis. Use the unique-peptide primary analysis.
sensitivity_path = root / (
    "write/39_independent_miller_mg_inflammatory_figures/Figure4/source_data/Figure4D_peptide_sensitivity_source.csv"
)
sensitivity_table = pd.read_csv(sensitivity_path)
sensitivity_table = sensitivity_table.loc[sensitivity_table["analysis"].eq("Unique peptides")]
if len(sensitivity_table) != 1:
    raise ValueError("The PDC unique-peptide statistics must have exactly one summary row.")
row = sensitivity_table.iloc[0]
if not (
    int(row["n_pairs"]) == 105
    and int(row["n_recurrence_up"]) == 72
    and np.isclose(float(row["mean_delta"]), program_delta_table["delta"].mean(), atol=1e-12)
):
    raise ValueError("The PDC patient-level summary statistics differ from the patient-level source table.")

# Copy the summary fields directly to the patient table, preserving the test and FDR definitions.
program_delta_table["sign_flip_p"] = float(row["sign_flip_fdr"])
program_delta_table["mean_delta"] = float(row["mean_delta"])
program_delta_table["ci_low"] = float(row["ci_low"])
program_delta_table["ci_high"] = float(row["ci_high"])
program_delta_table["n_recurrence_up"] = int(row["n_recurrence_up"])
program_delta_table["n_pairs"] = int(row["n_pairs"])
program_delta_table = program_delta_table.sort_values(["delta", "patient"], ascending=[True, True]).reset_index(
    drop = True
)
program_delta_table["patient_rank"] = np.arange(1, len(program_delta_table) + 1)
program_delta_table["direction"] = np.where(program_delta_table["delta"] > 0, "Recurrent higher", "Primary higher")


# Each row is one measured protein. Use the original proteome-wide FDR threshold of 0.05 and significant set.

statistics_path = step41_source / "Figure4B_measurable_protein_forest_source.csv"
gene_statistics_table = pd.read_csv(statistics_path)
gene_statistics_table["significant"] = gene_statistics_table["adj.P.Val"] < 0.05
if len(gene_statistics_table) != 11 or int(gene_statistics_table["significant"].sum()) != 5:
    raise ValueError("Expected 11 measurable PDC Miller-IM proteins, of which 5 pass the proteome-wide FDR threshold.")
expected = {"RHOB", "GSTM3", "FOLR2", "PDK4", "CTTNBP2"}
observed = set(gene_statistics_table.loc[gene_statistics_table["significant"], "gene"])
if observed != expected:
    raise ValueError(f"The set of significant PDC Miller-IM proteins has changed. Found {sorted(observed)}.")


# Retain the first record per Gene and use the original Unshared Log columns and T1/T2 pairing.
data_root = Path(os.environ.get("MILLER_IM_DATA_ROOT", root / "data")).expanduser().resolve()
pdc_matrix_path = data_root / (
    "write/31_proteomic_validation/source_metadata/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
)

if not pdc_matrix_path.exists():
    raise FileNotFoundError(f"The raw PDC matrix is unavailable at {pdc_matrix_path}.")

header_table = pd.read_csv(pdc_matrix_path, sep="\t", nrows=0)
usecols = ["Gene"] + [column for column in header_table.columns if column.startswith("Unshared Log KNCC_GBM")]
protein_table = pd.read_csv(pdc_matrix_path, sep="\t", usecols=usecols)

# Use the 20 fixed genes in their original program order.
miller_im_genes = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
    "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3"
]
protein_table = protein_table.loc[
    protein_table["Gene"].notna()
    & ~protein_table["Gene"].isin(["Mean", "Median", "StdDev"])
    & protein_table["Gene"].isin(miller_im_genes)
].drop_duplicates("Gene", keep="first")
protein_table = protein_table.set_index("Gene")

sample_pattern = re.compile(r"^Unshared Log (KNCC_GBM[0-9]+_T[12]):")
source_columns: dict[str, list[str]] = {}
for column in protein_table.columns:
    match = sample_pattern.match(column)
    if match:
        source_columns.setdefault(match.group(1), []).append(column)

# Average source columns per protein within each sample, giving one column per sample.
sample_expression = {}
for sample, columns in source_columns.items():
    sample_expression[sample] = protein_table[columns].mean(axis=1)
expression_table = pd.DataFrame(sample_expression)


# Use the original T1/T2 sample labels to identify strict patient pairs.
sample_rows = []
for sample in expression_table.columns:
    patient = re.sub(r"_T[12]$", "", sample)
    condition = "Recurrent" if sample.endswith("_T2") else "Primary"
    sample_rows.append((sample, patient, condition))
sample_table = pd.DataFrame(sample_rows, columns=["sample", "patient", "condition"])
sample_count_table = sample_table.pivot_table(
    index="patient", columns="condition", values="sample", aggfunc="nunique", fill_value=0
)
paired_patients = sample_count_table.index[
    (sample_count_table.get("Primary", 0) == 1) & (sample_count_table.get("Recurrent", 0) == 1)
]
if len(paired_patients) != 105:
    raise ValueError(f"Expected 105 patient pairs in the raw PDC matrix. Found {len(paired_patients)}.")


# Retain finite values and the original order for each patient-protein difference.
delta_records: list[dict[str, object]] = []
for patient in sorted(paired_patients):
    primary_sample = f"{patient}_T1"
    recurrent_sample = f"{patient}_T2"
    for gene in expression_table.index:
        primary_value = expression_table.at[gene, primary_sample]
        recurrent_value = expression_table.at[gene, recurrent_sample]
        if np.isfinite(primary_value) and np.isfinite(recurrent_value):
            delta_records.append(
                {
                    "patient": patient,
                    "gene": gene,
                    "Primary": primary_value,
                    "Recurrent": recurrent_value,
                    "delta_log2_abundance": recurrent_value - primary_value,
                }
            )
gene_delta_table_total = pd.DataFrame(delta_records)

present_genes = [gene for gene in miller_im_genes if gene in expression_table.index]

# Order each patient's T1 before T2 and keep the paired patient sequence for standardization.
paired_samples = []
for patient in paired_patients:
    paired_samples.append(f"{patient}_T1")
    paired_samples.append(f"{patient}_T2")

# Standardize each protein across all paired samples with sample standard deviation ddof = 1.
standardized_expression_table = (
    expression_table.loc[present_genes, paired_samples]
    .sub(expression_table.loc[present_genes, paired_samples].mean(axis=1), axis=0)
    .div(expression_table.loc[present_genes, paired_samples].std(axis=1, ddof=1), axis=0)
)

# Average proteins measured on both sides of each pair. The original score requires three or more shared proteins.
reconstructed_scores = []
for patient in sorted(paired_patients):
    primary = standardized_expression_table[f"{patient}_T1"]
    recurrent = standardized_expression_table[f"{patient}_T2"]
    common = primary.notna() & recurrent.notna()
    if int(common.sum()) < 3:
        continue
    reconstructed_scores.append(
        {
            "patient": patient,
            "reconstructed_delta": recurrent[common].mean() - primary[common].mean(),
            "reconstructed_n_common": int(common.sum()),
        }
    )
reconstructed_table = pd.DataFrame(reconstructed_scores)

# Match rebuilt patient differences and counts of shared proteins to Step 41.
score_comparison_table_total = program_delta_table[["patient", "delta", "n_common_raw20"]].merge(
    reconstructed_table, on="patient", how="outer", validate="one_to_one"
)
if score_comparison_table_total.isna().any().any():
    raise ValueError("Reconstructed PDC patient scores contain missing values.")
if np.max(np.abs(score_comparison_table_total["delta"] - score_comparison_table_total["reconstructed_delta"])) > 1e-10:
    raise ValueError("Patient-level PDC deltas reconstructed from the raw matrix differ from Step 41.")
if not np.array_equal(
    score_comparison_table_total["n_common_raw20"], score_comparison_table_total["reconstructed_n_common"]
):
    raise ValueError("The number of shared measurable PDC proteins per pair differs from Step 41.")


# Combine pair counts, increases and limma statistics in the original order for each protein.
counts_per_gene = gene_delta_table_total.groupby("gene").size().rename("n_pairs_gene")
up_per_gene = (
    gene_delta_table_total.assign(up=gene_delta_table_total["delta_log2_abundance"] > 0)
    .groupby("gene")["up"]
    .sum()
    .rename("n_recurrence_up_gene")
)
gene_delta_table_total = gene_delta_table_total.merge(counts_per_gene, on="gene").merge(up_per_gene, on="gene")
gene_delta_table_total = gene_delta_table_total.merge(
    gene_statistics_table[["gene", "logFC", "ci_low", "ci_high", "P.Value", "adj.P.Val", "significant"]],
    on  = "gene",
    how = "left",
    validate = "many_to_one",
)
observed_means = gene_delta_table_total.groupby("gene")["delta_log2_abundance"].mean()
expected_means = gene_statistics_table.set_index("gene")["logFC"]
common = observed_means.index.intersection(expected_means.index)
if np.max(np.abs(observed_means.loc[common] - expected_means.loc[common])) > 1e-10:
    raise ValueError("The mean raw patient-level delta differs from the paired limma logFC.")


# Use the same patient table and frozen 1.96 multiplier for the panel A mean and 95% interval.
records = []

stem = "Figure4A_pdc_patient_delta_waterfall"

# Write the panel source table in the original field and row order.
source_path = source_out / f"{stem}_source.csv"
program_delta_table.to_csv(source_path, index=False)

mean_delta = float(program_delta_table["delta"].mean())
standard_error = float(program_delta_table["delta"].std(ddof=1) / np.sqrt(len(program_delta_table)))
ci_low, ci_high = mean_delta - 1.96 * standard_error, mean_delta + 1.96 * standard_error

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(88, 67)
axis = figure.add_axes([0.13, 0.21, 0.82, 0.63])
colors = np.where(program_delta_table["delta"] > 0, "#E64B35", "#3C5488")
axis.bar(
    program_delta_table["patient_rank"],
    program_delta_table["delta"],
    width     = 0.88,
    color     = colors,
    edgecolor = "none",
    alpha     = 0.88,
    zorder    = 2,
)
axis.axhline(0, color="#222222", linewidth=0.72, zorder=3)
mean_x = 108.2
axis.plot([mean_x, mean_x], [ci_low, ci_high], color="#222222", linewidth=1.0, zorder=4)
axis.scatter([mean_x], [mean_delta], marker="D", s=28, color="white", edgecolor="#222222", linewidth=0.8, zorder=5)
axis.text(mean_x, ci_high + 0.08, "mean", ha="center", va="bottom", fontsize=5.8)
axis.text(
    3,
    max(program_delta_table["delta"]) * 0.92,
    "72/105 recurrent-higher",
    ha         = "left",
    va         = "top",
    fontsize   = 6.2,
    color      = "#E64B35",
    fontweight = "bold",
)

# Use the original axis and labels for the comparison.
axis.set_xlim(0, 111)
bound = max(abs(program_delta_table["delta"].min()), abs(program_delta_table["delta"].max())) * 1.08

# Keep the original axis limits and data order.
axis.set_ylim(-bound, bound)
axis.set_xticks([])
axis.set_xlabel("Paired patients ranked by change")
axis.set_ylabel("Recurrent − primary Miller-IM protein score")
clean_axis(axis, keep_left=True, keep_bottom=True)
axis.tick_params(axis="x", length=0)
figure.text(0.10, 0.955, "Miller-IM protein score rises at recurrence", fontsize=8.3, fontweight="bold", va="top")

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step54 recurrent GBM Figure 4")
records.append(
    panel_record("A", stem, "105 paired-patient Miller-IM protein-score changes", source_path, pdf_path, png_path)
)


# Panel B uses one row per patient-protein change and the original jitter sequence.

# Keep the patient-protein order and generate jitter once after sorting with the fixed seed.
protein_plot_table_total = gene_delta_table_total.copy()
order = gene_statistics_table.sort_values("logFC", ascending=False)["gene"].tolist()
y_lookup = {gene: len(order) - 1 - index for index, gene in enumerate(order)}
rng = np.random.default_rng(20260714)
protein_plot_table_total["y_base"] = protein_plot_table_total["gene"].map(y_lookup).astype(float)
protein_plot_table_total["y_jitter"] = protein_plot_table_total["y_base"] + rng.uniform(
    -0.20, 0.20, len(protein_plot_table_total)
)

# Write the panel source table in the original field and row order.
stem = "Figure4B_raw20_protein_patient_distributions"
source_path = source_out / f"{stem}_source.csv"
protein_plot_table_total.to_csv(source_path, index=False)

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(178, 83)
axis = figure.add_axes([0.115, 0.16, 0.855, 0.73])
for index, gene in enumerate(order):
    y_position = y_lookup[gene]
    if index % 2 == 0:
        axis.axhspan(y_position - 0.46, y_position + 0.46, color="#F7F7F7", zorder=0)
positive = protein_plot_table_total["delta_log2_abundance"] > 0
axis.scatter(
    protein_plot_table_total.loc[~positive, "delta_log2_abundance"],
    protein_plot_table_total.loc[~positive, "y_jitter"],
    s         = 7,
    color     = "#3C5488",
    alpha     = 0.35,
    edgecolor = "none",
    zorder    = 2,
)
axis.scatter(
    protein_plot_table_total.loc[positive, "delta_log2_abundance"],
    protein_plot_table_total.loc[positive, "y_jitter"],
    s         = 7,
    color     = "#E64B35",
    alpha     = 0.35,
    edgecolor = "none",
    zorder    = 2,
)

# Read the protein estimates, intervals and significance symbols from the frozen statistics.
protein_summary_table = gene_statistics_table.set_index("gene").loc[order].reset_index()
protein_summary_table["y"] = protein_summary_table["gene"].map(y_lookup)
axis.hlines(
    protein_summary_table["y"],
    protein_summary_table["ci_low"],
    protein_summary_table["ci_high"],
    color     = "#222222",
    linewidth = 1.1,
    zorder    = 4,
)
axis.scatter(
    protein_summary_table["logFC"],
    protein_summary_table["y"],
    marker    = "D",
    s         = 34,
    facecolor = np.where(protein_summary_table["significant"], "#E64B35", "white"),
    edgecolor = "#222222",
    linewidth = 0.75,
    zorder    = 5,
)
axis.axvline(0, color="#222222", linewidth=0.7, zorder=1)

# Use the original axis and labels for the comparison.
axis.set_xlim(-2.75, 2.75)

# Keep the original axis limits and data order.
axis.set_ylim(-0.6, len(order) - 0.4)
labels = []
for gene in reversed(order):
    n_pairs = int(protein_plot_table_total.loc[protein_plot_table_total["gene"] == gene, "n_pairs_gene"].iloc[0])
    labels.append(f"{gene}  ({n_pairs})")
axis.set_yticks(np.arange(len(order)), labels)
axis.set_xlabel("Paired protein change (recurrent − primary)")
axis.set_ylabel(None)
axis.grid(axis="x", color="#E7E7E7", linewidth=0.45, zorder=0)
clean_axis(axis, keep_left=False, keep_bottom=True)
axis.tick_params(axis="y", length=0, pad=3)

# Use legend symbols, colors and order that match the panel.
handles = [
    mpl.lines.Line2D(
        [],
        [],
        marker          = "D",
        linestyle       = "none",
        markerfacecolor = "#E64B35",
        markeredgecolor = "#222222",
        markersize      = 4.7,
        label           = "Full-proteome FDR < 0.05",
    ),
    mpl.lines.Line2D(
        [],
        [],
        marker          = "D",
        linestyle       = "none",
        markerfacecolor = "white",
        markeredgecolor = "#222222",
        markersize      = 4.7,
        label           = "FDR ≥ 0.05",
    ),
]
axis.legend(handles=handles, frameon=False, loc="lower right", fontsize=5.7, handletextpad=0.4)
figure.text(0.075, 0.96, "Five Miller-IM proteins increase significantly", fontsize=8.3, fontweight="bold", va="top")

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step54 recurrent GBM Figure 4")
records.append(
    panel_record(
        "B", stem, "Patient-level paired changes for 11 measurable Miller-IM proteins", source_path, pdf_path, png_path
    )
)


# Panel C uses the full proteome ranking, five label locations and the frozen enrichment summary.

curve_path = step41_source / "Figure4C_protein_ranked_gsea_curve_source.csv"
rank_curve_table = pd.read_csv(curve_path)
rank_plot_table_total = rank_curve_table.merge(
    gene_statistics_table[["gene", "logFC", "adj.P.Val", "significant"]],
    on  = "gene",
    how = "left",
    validate = "one_to_one",
)
rank_plot_table_total["measured_miller_im"] = rank_plot_table_total["logFC"].notna()
rank_sensitivity_table = pd.read_csv(sensitivity_path)
rank_sensitivity_table = rank_sensitivity_table.loc[rank_sensitivity_table["analysis"].eq("Unique peptides")]
if len(rank_sensitivity_table) != 1:
    raise ValueError("The PDC unique-peptide GSEA summary must have exactly one row.")
rank_plot_table_total["miller_im_nes"] = float(rank_sensitivity_table.iloc[0]["gsea_NES"])
rank_plot_table_total["miller_im_gsea_p"] = float(rank_sensitivity_table.iloc[0]["gsea_p"])
stem = "Figure4C_full_proteome_rank_landscape"

# Write the panel source table in the original field and row order.
source_path = source_out / f"{stem}_source.csv"
rank_plot_table_total.to_csv(source_path, index=False)

# Preserve the original panel dimensions and axes positions in millimeters.
figure = new_figure(88, 67)
axis = figure.add_axes([0.13, 0.21, 0.83, 0.63])
positive = rank_curve_table["rank_stat"] >= 0
axis.fill_between(
    rank_curve_table.loc[positive, "rank"],
    0,
    rank_curve_table.loc[positive, "rank_stat"],
    color     = "#F8DEDA",
    linewidth = 0,
    zorder    = 1,
)
axis.fill_between(
    rank_curve_table.loc[~positive, "rank"],
    0,
    rank_curve_table.loc[~positive, "rank_stat"],
    color     = "#DCE5F1",
    linewidth = 0,
    zorder    = 1,
)
axis.plot(rank_curve_table["rank"], rank_curve_table["rank_stat"], color="#9A9A9A", linewidth=0.55, zorder=2)
measured_protein_table = rank_plot_table_total.loc[rank_plot_table_total["measured_miller_im"]].copy()
significant = measured_protein_table["significant"].eq(True)
axis.scatter(
    measured_protein_table.loc[~significant, "rank"],
    measured_protein_table.loc[~significant, "rank_stat"],
    s         = 28,
    facecolor = "white",
    edgecolor = "#222222",
    linewidth = 0.7,
    zorder    = 5,
)
axis.scatter(
    measured_protein_table.loc[significant, "rank"],
    measured_protein_table.loc[significant, "rank_stat"],
    s         = 32,
    facecolor = "#E64B35",
    edgecolor = "#222222",
    linewidth = 0.7,
    zorder    = 6,
)

# Place the five protein labels at their original coordinates.
label_positions = {
    "RHOB": (700, 8.0, "left"),
    "GSTM3": (1350, 5.9, "left"),
    "FOLR2": (1850, 6.9, "left"),
    "PDK4": (2450, 3.8, "left"),
    "P2RY13": (7000, -3.8, "center"),
}
for row in measured_protein_table.itertuples():
    if row.gene not in label_positions:
        continue
    label_x, label_y, ha = label_positions[row.gene]
    axis.annotate(
        row.gene,
        (row.rank, row.rank_stat),
        xytext     = (label_x, label_y),
        textcoords = "data",
        ha         = ha,
        va         = "center",
        fontsize   = 5.9,
        color      = "#222222",
        arrowprops = {"arrowstyle": "-", "color": "#8A8A8A", "linewidth": 0.45},
    )
axis.axhline(0, color="#222222", linewidth=0.72, zorder=3)
axis.text(
    0.01,
    0.96,
    "recurrent-enriched",
    transform  = axis.transAxes,
    color      = "#E64B35",
    fontsize   = 6.2,
    fontweight = "bold",
    ha         = "left",
    va         = "top",
)
axis.text(
    0.99,
    0.05,
    "primary-enriched",
    transform  = axis.transAxes,
    color      = "#3C5488",
    fontsize   = 6.2,
    fontweight = "bold",
    ha         = "right",
    va         = "bottom",
)
axis.text(
    0.99, 0.96, "11/20 proteins measured", transform=axis.transAxes, color="#7A7A7A", fontsize=5.8, ha="right", va="top"
)

# Use the original axis and labels for the comparison.
axis.set_xlim(rank_curve_table["rank"].min(), rank_curve_table["rank"].max())
axis.set_xlabel("All proteins ranked by paired differential statistic")
axis.set_ylabel("Moderated t statistic")
axis.grid(axis="y", color="#E7E7E7", linewidth=0.45, zorder=0)
clean_axis(axis, keep_left=True, keep_bottom=True)
figure.text(
    0.10, 0.955, "Miller-IM proteins favor the recurrent proteome rank", fontsize=8.3, fontweight="bold", va="top"
)

# Save the original PDF layout and render the 600 dpi PNG from it.
pdf_path, png_path = save_panel(figure, figure_out, stem, "Step54 recurrent GBM Figure 4")
records.append(
    panel_record(
        "C", stem, "All-protein rank landscape with 11 measurable Miller-IM proteins", source_path, pdf_path, png_path
    )
)


# Keep the original PIL canvas, arrangement and 180 dpi setting for the preview.

canvas = Image.new("RGB", (3800, 3000), "white")
drawer = ImageDraw.Draw(canvas)
if [row["panel"] for row in records] != list("ABC"):
    raise ValueError("Figure 4 must contain panels A, B and C exactly once.")

by_panel = {row["panel"]: Path(row["png"]) for row in records}
layout = [
    {"panel": "A", "left": 100, "top": 120, "right": 1870, "bottom": 1420},
    {"panel": "C", "left": 1930, "top": 120, "right": 3700, "bottom": 1420},
    {"panel": "B", "left": 100, "top": 1540, "right": 3700, "bottom": 2920},
]
panel_font = image_font(58, bold=True)
for item in layout:
    panel = item["panel"]
    paste_contained(
        canvas,
        by_panel[panel],
        (item["left"], item["top"], item["right"], item["bottom"]),
    )
    drawer.text((item["left"], item["top"] - 60), panel, fill="#111111", font=panel_font)

# Keep the preview's pixel dimensions and resolution for composite layout review.
png_path = figure_out / "Figure4_proteome_preview.png"
canvas.save(png_path, format="PNG", dpi=(180, 180), optimize=True)
preview_gray_path = figure_out / "Figure4_proteome_preview_grayscale.png"
ImageOps.grayscale(canvas).convert("RGB").save(preview_gray_path, format="PNG", dpi=(180, 180), optimize=True)
preview = plt.figure(figsize=(3800 / 180, 3000 / 180), dpi=180, facecolor="white")
preview_ax = preview.add_axes([0, 0, 1, 1])
preview_ax.imshow(canvas)
preview_ax.axis("off")
pdf_path = figure_out / "Figure4_proteome_preview.pdf"
preview.canvas.draw()
preview.savefig(
    pdf_path,
    dpi       = 180,
    facecolor = "white",
    edgecolor = "none",
    metadata  = {"Creator": "Step54 recurrent GBM Figure 4", "CreationDate": None, "ModDate": None},
)
plt.close(preview)
layout_table = pd.DataFrame(layout)
source_path = source_out / "Figure4_proteome_preview_layout.csv"
layout_table.to_csv(source_path, index=False)


# Each manifest row records one panel with the existing public fields and preview references.
manifest_table = pd.DataFrame(records)
manifest_table["preview_pdf"] = str(pdf_path)
manifest_table["preview_png"] = str(png_path)
manifest_table["preview_grayscale"] = str(preview_gray_path)
manifest_table["preview_layout"] = str(source_path)
manifest_path = write_root / "Figure4_panel_manifest.csv"
manifest_table.to_csv(manifest_path, index=False)

# Check the files and preview dimensions produced by this run. Raise on missing or empty output.
output_paths = [
    manifest_path,
    pdf_path,
    png_path,
    preview_gray_path,
    source_path,
]
for panel_result in records:
    for field in ["source", "pdf", "png"]:
        output_paths.append(Path(panel_result[field]))
for output_path in output_paths:

    # Check that each current output exists and has content.
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {output_path.name}")
with Image.open(png_path) as preview_image:
    if preview_image.size != canvas.size:
        raise ValueError("Preview dimensions differ from the original canvas")
    preview_image.verify()

# Read the saved manifest and compare its dimensions to the current panel records.
manifest_check_table = pd.read_csv(manifest_path)
if manifest_check_table.shape != manifest_table.shape:
    raise ValueError("Panel manifest dimensions changed during writing")
print(f"Figure4 panels={len(records)} preview=1")
