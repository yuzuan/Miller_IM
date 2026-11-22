#!/usr/bin/env python3
"""Render the paired PDC000514 proteome panels for Figure 4.

Inputs
    Under NOA_LEGACY_ROOT, use write/41_mg_inflammatory_sci_rebuild/Figure4/source_data/:
        Figure4A_patient_delta_source.csv, Figure4B_measurable_protein_forest_source.csv,
        Figure4C_protein_ranked_gsea_curve_source.csv.
    Read the unique-peptide row from
    write/39_independent_miller_mg_inflammatory_figures/Figure4/source_data/Figure4D_peptide_sensitivity_source.csv
    and the archived per-patient protein changes from
    write/54_figure4_proteome_rebuild/Figure4/source_data/Figure4B_raw20_protein_patient_distributions_source.csv.
    The separate raw reconstruction uses NOA_PUBLIC_DATA_ROOT/write/31_proteomic_validation/source_metadata/
    PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv.
Outputs
    Under _figure_build/54_figure4_proteome_rebuild/Figure4, write Figure4_panel_manifest.csv,
    Figure4_legend.md, FINAL_RESULT.md and source_data/:
        Figure4A_pdc_patient_delta_waterfall_source.csv, Figure4B_raw20_protein_patient_distributions_source.csv,
        Figure4C_full_proteome_rank_landscape_source.csv, Figure4_proteome_preview_layout.csv.
    Under _panels/54_figure4_proteome_rebuild/Figure4, write PDF/PNG pairs for the three data panels,
    plus Figure4_proteome_preview.pdf, Figure4_proteome_preview.png and Figure4_proteome_preview_grayscale.png.
Run
    python scripts/08d_figure4.py from NOA_R1, or use scripts/09b_export_final_figures.py for the final figure.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps


figure_tools = importlib.import_module("00a_functions")
revision_root = Path(__file__).resolve().parents[1]
project_root = figure_tools.legacy_root(revision_root)

# The final exporter rebinds these directories for its isolated build.
WRITE_ROOT = revision_root / "_figure_build/54_figure4_proteome_rebuild/Figure4"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "_panels/54_figure4_proteome_rebuild/Figure4"

# Reuse the original across-panel colors for signed effects, text, and grid lines.
primary = "#3C5488"
recurrent = "#E64B35"
text_color = "#222222"
grid_color = "#E7E7E7"


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
) -> dict[str, str]:
    """Return the original per-panel manifest fields from the three written files."""
    for path in (source, pdf, png):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 4 panel output: {path.name}")
    panel_record = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": figure_tools.sha256(source),
        "pdf_sha256": figure_tools.sha256(pdf),
        "png_sha256": figure_tools.sha256(png),
    }
    return panel_record


def load_program_deltas() -> pd.DataFrame:
    """Return the frozen unique-peptide program changes and statistics in waterfall order."""
    path = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure4/source_data/Figure4A_patient_delta_source.csv"
    paired_table = pd.read_csv(path)
    if len(paired_table) != 105:
        raise ValueError(f"PDC000514 needs 105 paired patients, got {len(paired_table)}")
    if int((paired_table["delta"] > 0).sum()) != 72:
        raise ValueError("PDC000514 needs 72 positive paired changes in the frozen cohort")
    if not np.isclose(paired_table["delta"].mean(), 0.252534746498684, atol=1e-12):
        raise ValueError("The PDC000514 mean paired change differs from the frozen result")

    # Check the archived gene-set summary against the paired changes before attaching its statistics.
    summary_path = project_root / (
        "write/39_independent_miller_mg_inflammatory_figures/Figure4/source_data/"
        "Figure4D_peptide_sensitivity_source.csv"
    )
    summary_table = pd.read_csv(summary_path)
    unique_peptide_table = summary_table.loc[summary_table["analysis"].eq("Unique peptides")]
    if len(unique_peptide_table) != 1:
        raise ValueError("PDC000514 needs exactly one unique-peptide summary row")
    row = unique_peptide_table.iloc[0]
    if not (
        int(row["n_pairs"]) == 105
        and int(row["n_recurrence_up"]) == 72
        and np.isclose(float(row["mean_delta"]), paired_table["delta"].mean(), atol=1e-12)
    ):
        raise ValueError("The unique-peptide summary and paired PDC000514 scores disagree")
    paired_table["sign_flip_p"] = float(row["sign_flip_fdr"])
    paired_table["mean_delta"] = float(row["mean_delta"])
    paired_table["ci_low"] = float(row["ci_low"])
    paired_table["ci_high"] = float(row["ci_high"])
    paired_table["n_recurrence_up"] = int(row["n_recurrence_up"])
    paired_table["n_pairs"] = int(row["n_pairs"])

    # Rank patients by recurrent-minus-primary change, breaking ties by patient identifier.
    paired_table = paired_table.sort_values(["delta", "patient"], ascending=[True, True]).reset_index(drop=True)
    paired_table["patient_rank"] = np.arange(1, len(paired_table) + 1)
    paired_table["direction"] = np.where(paired_table["delta"] > 0, "Recurrent higher", "Primary higher")
    return paired_table


def load_gene_statistics() -> pd.DataFrame:
    """Return the paired limma effects with their full-proteome FDR status."""
    path = project_root / (
        "write/41_mg_inflammatory_sci_rebuild/Figure4/source_data/Figure4B_measurable_protein_forest_source.csv"
    )
    statistics_table = pd.read_csv(path)
    statistics_table["significant"] = statistics_table["adj.P.Val"] < 0.05
    if len(statistics_table) != 11 or int(statistics_table["significant"].sum()) != 5:
        raise ValueError("Figure 4 needs 11 measurable proteins and five with full-proteome FDR below 0.05")
    expected = {"RHOB", "GSTM3", "FOLR2", "PDK4", "CTTNBP2"}
    observed = set(statistics_table.loc[statistics_table["significant"], "gene"])
    if observed != expected:
        raise ValueError(f"Unexpected significant MG-Inflam20 proteins: {sorted(observed)}")
    return statistics_table


def rebuild_patient_gene_deltas(
    program_delta_table: pd.DataFrame,
    gene_statistics_table: pd.DataFrame,
) -> pd.DataFrame:
    """Return raw-matrix paired protein changes after matching the prepared effects."""
    matrix_path = figure_tools.public_data_root(project_root) / (
        "write/31_proteomic_validation/source_metadata/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
    )
    header_table = pd.read_csv(matrix_path, sep="\t", nrows=0)
    use_columns = ["Gene"] + [column for column in header_table.columns if column.startswith("Unshared Log KNCC_GBM")]
    raw_table = pd.read_csv(matrix_path, sep="\t", usecols=use_columns)

    # Keep the frozen MG-Inflam20 definition and the first row of each measured protein.
    miller_im_genes = [
        "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CCL3", "CH25H",
        "P2RY12", "JUN", "SIGLEC8", "KLF6", "FOLR2", "CCL4", "AC253572.2", "NLRP3",
    ]
    raw_table = raw_table.loc[
        raw_table["Gene"].notna()
        & ~raw_table["Gene"].isin(["Mean", "Median", "StdDev"])
        & raw_table["Gene"].isin(miller_im_genes)
    ]
    raw_table = raw_table.drop_duplicates("Gene", keep="first")
    raw_table = raw_table.set_index("Gene")

    # Average the unique-peptide abundance columns that belong to the same sample.
    sample_pattern = re.compile(r"^Unshared Log (KNCC_GBM[0-9]+_T[12]):")
    source_columns: dict[str, list[str]] = {}
    for column in raw_table.columns:
        match = sample_pattern.match(column)
        if match:
            source_columns.setdefault(match.group(1), []).append(column)
    sample_abundances = {}
    for sample, columns in source_columns.items():
        sample_abundances[sample] = raw_table[columns].mean(axis=1)
    abundance_table = pd.DataFrame(sample_abundances)

    # Retain patients with one primary and one recurrent sample.
    sample_rows = []
    for sample in abundance_table.columns:
        patient = re.sub(r"_T[12]$", "", sample)
        condition = "Recurrent" if sample.endswith("_T2") else "Primary"
        sample_rows.append((sample, patient, condition))
    sample_table = pd.DataFrame(sample_rows, columns=["sample", "patient", "condition"])
    count_table = sample_table.pivot_table(
        index      = "patient",
        columns    = "condition",
        values     = "sample",
        aggfunc    = "nunique",
        fill_value = 0,
    )
    paired = count_table.index[(count_table.get("Primary", 0) == 1) & (count_table.get("Recurrent", 0) == 1)]
    if len(paired) != 105:
        raise ValueError(f"The PDC000514 matrix must contain 105 paired patients, got {len(paired)}")

    # A protein contributes a change only when both values in a pair are finite.
    records: list[dict[str, object]] = []
    for patient in sorted(paired):
        primary_sample = f"{patient}_T1"
        recurrent_sample = f"{patient}_T2"
        for gene in abundance_table.index:
            primary_value = abundance_table.at[gene, primary_sample]
            recurrent_value = abundance_table.at[gene, recurrent_sample]
            if np.isfinite(primary_value) and np.isfinite(recurrent_value):
                records.append(
                    {
                        "patient": patient,
                        "gene": gene,
                        "Primary": primary_value,
                        "Recurrent": recurrent_value,
                        "delta_log2_abundance": recurrent_value - primary_value,
                    }
                )
    delta_table = pd.DataFrame(records)

    # Standardize each protein across paired specimens and score pairs with at least three common proteins.
    present = [gene for gene in miller_im_genes if gene in abundance_table.index]
    paired_samples = [item for patient in paired for item in (f"{patient}_T1", f"{patient}_T2")]
    paired_abundance_table = abundance_table.loc[present, paired_samples]
    centered_table = paired_abundance_table.sub(paired_abundance_table.mean(axis=1), axis=0)
    z_expression_table = centered_table.div(paired_abundance_table.std(axis=1, ddof=1), axis=0)
    reconstructed_scores = []
    for patient in sorted(paired):
        primary_values = z_expression_table[f"{patient}_T1"]
        recurrent_values = z_expression_table[f"{patient}_T2"]
        common = primary_values.notna() & recurrent_values.notna()
        if int(common.sum()) < 3:
            continue
        reconstructed_scores.append(
            {
                "patient": patient,
                "reconstructed_delta": recurrent_values[common].mean() - primary_values[common].mean(),
                "reconstructed_n_common": int(common.sum()),
            }
        )
    reconstructed_table = pd.DataFrame(reconstructed_scores)
    matched_table = program_delta_table[["patient", "delta", "n_common_raw20"]].merge(
        reconstructed_table, on="patient", how="outer", validate="one_to_one"
    )
    if matched_table.isna().any().any():
        raise ValueError("A PDC000514 patient is missing a reconstructed or archived score")
    if np.max(np.abs(matched_table["delta"] - matched_table["reconstructed_delta"])) > 1e-10:
        raise ValueError("Raw-matrix reconstruction differs from the archived PDC000514 patient changes")
    if not np.array_equal(matched_table["n_common_raw20"], matched_table["reconstructed_n_common"]):
        raise ValueError("Counts of common proteins differ from the archived PDC000514 patient data")

    # Add the paired counts and the archived full-proteome statistics to each protein record.
    counts_per_gene = delta_table.groupby("gene").size().rename("n_pairs_gene")
    signed_delta_table = delta_table.assign(up=delta_table["delta_log2_abundance"] > 0)
    up_per_gene = signed_delta_table.groupby("gene")["up"].sum().rename("n_recurrence_up_gene")
    delta_table = delta_table.merge(counts_per_gene, on="gene").merge(up_per_gene, on="gene")
    delta_table = delta_table.merge(
        gene_statistics_table[["gene", "logFC", "ci_low", "ci_high", "P.Value", "adj.P.Val", "significant"]],
        on  = "gene",
        how = "left",
        validate = "many_to_one",
    )
    observed_means = delta_table.groupby("gene")["delta_log2_abundance"].mean()
    expected_means = gene_statistics_table.set_index("gene")["logFC"]
    common = observed_means.index.intersection(expected_means.index)
    if np.max(np.abs(observed_means.loc[common] - expected_means.loc[common])) > 1e-10:
        raise ValueError("Raw paired protein changes disagree with the archived limma log fold changes")
    return delta_table


def panel_a(program_delta_table: pd.DataFrame) -> dict[str, str]:
    """Plot paired score changes in rank order and the pooled mean with its original interval."""
    stem = "Figure4A_pdc_patient_delta_waterfall"
    source = SOURCE_OUT / f"{stem}_source.csv"
    program_delta_table.to_csv(source, index=False)

    # The pooled marker uses the original normal-approximation interval for the paired mean.
    changes = program_delta_table["delta"]
    mean_delta = float(changes.mean())
    standard_error = float(changes.std(ddof=1) / np.sqrt(len(program_delta_table)))
    ci_low, ci_high = mean_delta - 1.96 * standard_error, mean_delta + 1.96 * standard_error

    fig = figure_tools.new_figure(88, 67)
    ax = fig.add_axes([0.13, 0.21, 0.82, 0.63])
    colors = np.where(changes > 0, recurrent, primary)
    ax.bar(
        program_delta_table["patient_rank"],
        changes,
        width     = 0.88,
        color     = colors,
        edgecolor = "none",
        alpha     = 0.88,
        zorder    = 2,
    )
    ax.axhline(0, color=text_color, linewidth=0.72, zorder=3)

    # The mean diamond is offset beyond the final patient bar.
    mean_x = 108.2
    ax.plot([mean_x, mean_x], [ci_low, ci_high], color=text_color, linewidth=1.0, zorder=4)
    ax.scatter([mean_x], [mean_delta], marker="D", s=28, color="white", edgecolor=text_color, linewidth=0.8, zorder=5)
    ax.text(mean_x, ci_high + 0.08, "mean", ha="center", va="bottom", fontsize=6.1)
    positive_count = int((changes > 0).sum())
    ax.text(
        3,
        max(changes) * 0.92,
        f"{positive_count}/{len(program_delta_table)} recurrent-higher",
        ha         = "left",
        va         = "top",
        fontsize   = 6.2,
        color      = recurrent,
        fontweight = "bold",
    )
    ax.set_xlim(0, 111)
    bound = max(abs(changes.min()), abs(changes.max())) * 1.08
    ax.set_ylim(-bound, bound)
    ax.set_xticks([])
    ax.set_xlabel("Paired patients ranked by change")
    ax.set_ylabel("Recurrent − primary MG-Inflam20 protein score")
    figure_tools.clean_axis(ax, keep_left=True, keep_bottom=True)
    ax.tick_params(axis="x", length=0)
    fig.text(0.10, 0.955, "MG-Inflam20 protein score rises at recurrence", fontsize=8.3, fontweight="bold", va="top")
    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step54 recurrent GBM Figure 4")
    panel_record = record(
        "A", stem, f"{len(program_delta_table)} paired-patient MG-Inflam20 protein-score changes", source, pdf, png
    )
    return panel_record


def panel_b(gene_delta_table: pd.DataFrame, gene_statistics_table: pd.DataFrame) -> dict[str, str]:
    """Plot the patient change distributions with paired limma estimates and full-proteome FDR."""
    stem = "Figure4B_raw20_protein_patient_distributions"
    source = SOURCE_OUT / f"{stem}_source.csv"
    plot_table = gene_delta_table.copy()
    order = gene_statistics_table.sort_values("logFC", ascending=False)["gene"].tolist()
    y_lookup = {gene: len(order) - 1 - index for index, gene in enumerate(order)}

    # The fixed seed and original row order preserve the stored jittered y coordinates.
    rng = np.random.default_rng(20260714)
    plot_table["y_base"] = plot_table["gene"].map(y_lookup).astype(float)
    plot_table["y_jitter"] = plot_table["y_base"] + rng.uniform(-0.20, 0.20, len(plot_table))
    plot_table.to_csv(source, index=False)

    # Alternate row shading and separate positive and nonpositive patient-level changes by color.
    fig = figure_tools.new_figure(178, 83)
    ax = fig.add_axes([0.16, 0.16, 0.80, 0.69])
    for index, gene in enumerate(order):
        row_y = y_lookup[gene]
        if index % 2 == 0:
            ax.axhspan(row_y - 0.46, row_y + 0.46, color="#F7F7F7", zorder=0)
    positive = plot_table["delta_log2_abundance"] > 0
    ax.scatter(
        plot_table.loc[~positive, "delta_log2_abundance"],
        plot_table.loc[~positive, "y_jitter"],
        s         = 7,
        color     = primary,
        alpha     = 0.35,
        edgecolor = "none",
        zorder    = 2,
    )
    ax.scatter(
        plot_table.loc[positive, "delta_log2_abundance"],
        plot_table.loc[positive, "y_jitter"],
        s         = 7,
        color     = recurrent,
        alpha     = 0.35,
        edgecolor = "none",
        zorder    = 2,
    )

    # Mark the limma estimates with archived confidence limits and FDR status.
    forest_table = gene_statistics_table.set_index("gene").loc[order].reset_index()
    forest_table["y"] = forest_table["gene"].map(y_lookup)
    ax.hlines(
        forest_table["y"], forest_table["ci_low"], forest_table["ci_high"], color=text_color, linewidth=1.1, zorder=4
    )
    ax.scatter(
        forest_table["logFC"],
        forest_table["y"],
        marker    = "D",
        s         = 34,
        facecolor = np.where(forest_table["significant"], recurrent, "white"),
        edgecolor = text_color,
        linewidth = 0.75,
        zorder    = 5,
    )
    ax.axvline(0, color=text_color, linewidth=0.7, zorder=1)
    ax.set_xlim(-2.75, 2.75)
    ax.set_ylim(-0.6, len(order) - 0.4)
    labels = []
    for gene in reversed(order):
        n_pairs = int(plot_table.loc[plot_table["gene"] == gene, "n_pairs_gene"].iloc[0])
        labels.append(f"{gene}  ({n_pairs})")
    ax.set_yticks(np.arange(len(order)), labels)
    ax.set_xlabel("Paired protein change (recurrent − primary)")
    ax.set_ylabel(None)
    ax.grid(axis="x", color=grid_color, linewidth=0.45, zorder=0)
    figure_tools.clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)

    # The legend and heading describe the five proteins verified in load_gene_statistics.
    handles = []
    for fill, label in [(recurrent, "Full-proteome FDR < 0.05"), ("white", "FDR ≥ 0.05")]:
        handles.append(
            mpl.lines.Line2D(
                [], [],
                marker          = "D",
                linestyle       = "none",
                markerfacecolor = fill,
                markeredgecolor = text_color,
                markersize      = 4.7,
                label           = label,
            )
        )
    ax.legend(
        handles        = handles,
        frameon        = False,
        loc            = "lower right",
        bbox_to_anchor = (1, 1.02),
        borderaxespad  = 0,
        ncol           = 2,
        fontsize       = 5.7,
        handletextpad  = 0.4,
    )
    fig.text(0.075, 0.96, "Five MG-Inflam20 proteins increase significantly", fontsize=8.3, fontweight="bold", va="top")
    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step54 recurrent GBM Figure 4")
    panel_record = record(
        "B", stem, f"Patient-level paired changes for {len(order)} measurable MG-Inflam20 proteins", source, pdf, png
    )
    return panel_record


def panel_c(gene_statistics_table: pd.DataFrame) -> dict[str, str]:
    """Plot the full proteome rank and annotate the frozen measured gene set."""
    stem = "Figure4C_full_proteome_rank_landscape"
    curve_path = project_root / (
        "write/41_mg_inflammatory_sci_rebuild/Figure4/source_data/Figure4C_protein_ranked_gsea_curve_source.csv"
    )
    rank_table = pd.read_csv(curve_path)
    plot_table = rank_table.merge(
        gene_statistics_table[["gene", "logFC", "adj.P.Val", "significant"]],
        on  = "gene",
        how = "left",
        validate = "one_to_one",
    )
    plot_table["measured_miller_im"] = plot_table["logFC"].notna()

    # Record the enrichment statistics for the same unique-peptide proteome rank.
    summary_path = project_root / (
        "write/39_independent_miller_mg_inflammatory_figures/Figure4/source_data/"
        "Figure4D_peptide_sensitivity_source.csv"
    )
    summary_table = pd.read_csv(summary_path)
    summary_table = summary_table.loc[summary_table["analysis"].eq("Unique peptides")]
    if len(summary_table) != 1:
        raise ValueError("PDC000514 needs exactly one unique-peptide GSEA summary row")
    plot_table["miller_im_nes"] = float(summary_table.iloc[0]["gsea_NES"])
    plot_table["miller_im_gsea_p"] = float(summary_table.iloc[0]["gsea_p"])
    source = SOURCE_OUT / f"{stem}_source.csv"
    plot_table.to_csv(source, index=False)

    # Shade positive and negative ranks using the original recurrent and primary tints.
    fig = figure_tools.new_figure(88, 67)
    ax = fig.add_axes([0.20, 0.21, 0.76, 0.63])
    positive = rank_table["rank_stat"] >= 0
    ax.fill_between(
        rank_table.loc[positive, "rank"],
        0,
        rank_table.loc[positive, "rank_stat"],
        color     = "#F8DEDA",
        linewidth = 0,
        zorder    = 1,
    )
    ax.fill_between(
        rank_table.loc[~positive, "rank"],
        0,
        rank_table.loc[~positive, "rank_stat"],
        color     = "#DCE5F1",
        linewidth = 0,
        zorder    = 1,
    )
    ax.plot(rank_table["rank"], rank_table["rank_stat"], color="#9A9A9A", linewidth=0.55, zorder=2)

    # Highlight all measured MG-Inflam20 proteins and fill those that pass the full-proteome FDR threshold.
    measured_table = plot_table.loc[plot_table["measured_miller_im"]].copy()
    significant = measured_table["significant"].eq(True)
    ax.scatter(
        measured_table.loc[~significant, "rank"],
        measured_table.loc[~significant, "rank_stat"],
        s         = 28,
        facecolor = "white",
        edgecolor = text_color,
        linewidth = 0.7,
        zorder    = 5,
    )
    ax.scatter(
        measured_table.loc[significant, "rank"],
        measured_table.loc[significant, "rank_stat"],
        s         = 32,
        facecolor = recurrent,
        edgecolor = text_color,
        linewidth = 0.7,
        zorder    = 6,
    )

    # Use the original offsets for the five selected rank labels.
    label_positions = {
        "RHOB": (650, 7.8, "left"),
        "GSTM3": (1300, 5.1, "left"),
        "FOLR2": (1750, 7.1, "left"),
        "PDK4": (2450, 3.8, "left"),
        "P2RY13": (7000, -3.8, "center"),
    }
    for row in measured_table.itertuples():
        if row.gene not in label_positions:
            continue
        label_x, label_y, ha = label_positions[row.gene]
        ax.annotate(
            row.gene,
            (row.rank, row.rank_stat),
            xytext     = (label_x, label_y),
            textcoords = "data",
            ha         = ha,
            va         = "center",
            fontsize   = 5.9,
            color      = text_color,
            arrowprops = {"arrowstyle": "-", "color": "#8A8A8A", "linewidth": 0.45},
        )

    # Label the sign of the moderated t-statistic and the coverage of the prespecified 20-gene program.
    ax.axhline(0, color=text_color, linewidth=0.72, zorder=3)
    ax.text(
        0.01, 0.96, "recurrent-enriched",
        transform  = ax.transAxes,
        color      = recurrent,
        fontsize   = 6.2,
        fontweight = "bold",
        ha         = "left",
        va         = "top",
    )
    ax.text(
        0.99, 0.05, "primary-enriched",
        transform  = ax.transAxes,
        color      = primary,
        fontsize   = 6.2,
        fontweight = "bold",
        ha         = "right",
        va         = "bottom",
    )
    ax.text(
        0.99, 0.96, f"{len(measured_table)}/20 proteins measured",
        transform = ax.transAxes,
        color     = "#7A7A7A",
        fontsize  = 6.1,
        ha        = "right",
        va        = "top",
    )
    ax.set_xlim(rank_table["rank"].min(), rank_table["rank"].max())
    ax.set_xlabel("All proteins ranked by paired differential statistic")
    ax.set_ylabel("Moderated t-statistic", fontsize=10)
    ax.grid(axis="y", color=grid_color, linewidth=0.45, zorder=0)
    figure_tools.clean_axis(ax, keep_left=True, keep_bottom=True)
    fig.text(
        0.10, 0.955, "MG-Inflam20 proteins favor the recurrent proteome rank",
        fontsize   = 8.3,
        fontweight = "bold",
        va         = "top",
    )
    pdf, png = figure_tools.save_panel(fig, FIG_OUT, stem, "Step54 recurrent GBM Figure 4")
    panel_record = record(
        "C", stem, f"All-protein rank landscape with {len(measured_table)} measurable MG-Inflam20 proteins",
        source, pdf, png,
    )
    return panel_record


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    """Place the three panel images into the original color and grayscale composites."""
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 100, "top": 120, "right": 1870, "bottom": 1420},
        {"panel": "C", "left": 1930, "top": 120, "right": 3700, "bottom": 1420},
        {"panel": "B", "left": 100, "top": 1540, "right": 3700, "bottom": 2920},
    ]
    canvas = Image.new("RGB", (3800, 3000), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = figure_tools.image_font(66, bold=True)
    for item in layout:
        panel = item["panel"]
        figure_tools.paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 60), panel, fill="#111111", font=panel_font)

    # The composite's archived PNG and PDF resolutions are 180 dpi.
    png = FIG_OUT / "Figure4_proteome_preview.png"
    gray = FIG_OUT / "Figure4_proteome_preview_grayscale.png"
    pdf = FIG_OUT / "Figure4_proteome_preview.pdf"
    source = SOURCE_OUT / "Figure4_proteome_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(180, 180), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(180, 180), optimize=True)
    preview = plt.figure(figsize=(3800 / 180, 3000 / 180), dpi=180, facecolor="white")
    preview_ax = preview.add_axes([0, 0, 1, 1])
    preview_ax.imshow(canvas)
    preview_ax.axis("off")
    preview.canvas.draw()
    preview.savefig(
        pdf,
        dpi       = 180,
        facecolor = "white",
        edgecolor = "none",
        metadata  = {"Creator": "Step54 recurrent GBM Figure 4", "CreationDate": None, "ModDate": None},
    )
    plt.close(preview)
    layout_table = pd.DataFrame(layout)
    layout_table.to_csv(source, index=False)
    output_paths = (pdf, png, gray, source)
    return output_paths


def main() -> None:
    """Render all three panels using the archived source tables and write the preview manifest."""
    figure_tools.apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    # The standard figure uses the archived per-patient protein changes, as in the October build.
    program_delta_table = load_program_deltas()
    gene_statistics_table = load_gene_statistics()
    gene_delta_path = project_root / (
        "write/54_figure4_proteome_rebuild/Figure4/source_data/Figure4B_raw20_protein_patient_distributions_source.csv"
    )
    gene_delta_table = pd.read_csv(gene_delta_path)
    records = [
        panel_a(program_delta_table),
        panel_b(gene_delta_table, gene_statistics_table),
        panel_c(gene_statistics_table),
    ]
    if [row["panel"] for row in records] != list("ABC"):
        raise ValueError("Figure 4 requires panels A through C in order")

    # Preserve the manifest fields and keep a source layout table for the composite.
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_table = pd.DataFrame(records)
    manifest_table["preview_pdf"] = str(preview_pdf)
    manifest_table["preview_png"] = str(preview_png)
    manifest_table["preview_grayscale"] = str(preview_gray)
    manifest_table["preview_layout"] = str(preview_layout)
    manifest_path = WRITE_ROOT / "Figure4_panel_manifest.csv"
    manifest_table.to_csv(manifest_path, index=False)

    # Write the archived Figure 4 legend used by the standalone renderer.
    legend_path = WRITE_ROOT / "Figure4_legend.md"
    legend_text = (
        "# Figure 4 legend\n\n"
        "**Figure 4 | Paired tumor proteomes support elevation of the measurable MG-Inflam20 protein program at "
        "recurrence.**\n\n"
        "**A,** Patient-level recurrent-minus-primary changes in the prespecified MG-Inflam20 protein score in "
        "PDC000514. Unique-peptide protein abundances were used for the primary analysis. Each protein was "
        "standardized across the paired-sample matrix, and each patient score was calculated from proteins measured "
        "in both specimens of that pair; at least three common proteins were required. Eleven of 20 MG-Inflam20 "
        "proteins were measurable globally, while individual pairs contained 7-11 common proteins. Seventy-two of "
        "105 paired patients were recurrent-higher; the mean paired change was +0.253 (95% CI, +0.141 to +0.364; "
        "prespecified paired sign-flip P = 4.0 × 10^-5). The diamond and interval show the mean and 95% CI. "
        "**B,** Patient-level recurrent-minus-primary changes for the 11 measurable MG-Inflam20 proteins. Small "
        "points denote individual patient pairs; parentheses after gene labels show the number of pairs with both "
        "measurements. Diamonds and horizontal intervals show paired limma effect estimates and approximate 95% CIs. "
        "Filled diamonds denote full-proteome Benjamini-Hochberg FDR < 0.05. RHOB, GSTM3, FOLR2, PDK4, and CTTNBP2 "
        "met this threshold. **C,** All 11,320 quantified proteins ranked by the paired moderated t statistic. "
        "The 11 measurable MG-Inflam20 proteins are overlaid; filled points denote full-proteome FDR < 0.05. "
        "The prespecified MG-Inflam20 gene set showed positive single-set competitive enrichment in the unique-peptide "
        "analysis (NES = 1.53; nominal P = 0.049), whereas the all-peptide sensitivity analysis was nonsignificant "
        "and is retained in the supplement.\n\n"
        "PDC000514 is a whole-tumor TMT proteomic dataset rather than a myeloid-enriched assay. These results "
        "therefore provide orthogonal support for recurrence-associated elevation of the measurable protein program "
        "but do not identify the cellular source, prove an MCG1/MCG2 identity, or establish that every MG-Inflam20 "
        "gene is translated and detected. Clinical metadata available to this analysis did not contain an IDH field, "
        "so the protein cohort is not labeled as IDH-wild-type.\n"
    )
    legend_path.write_text(legend_text, encoding="utf-8")

    # Keep the existing panel summary and verify the written composite and manifest.
    result_path = WRITE_ROOT / "FINAL_RESULT.md"
    result_text = (
        "# Figure 4 proteome rebuild\n\n"
        f"- A. MG-Inflam20 protein score changes for {len(program_delta_table)} patient pairs.\n"
        f"- B. Patient-level changes for {len(gene_statistics_table)} measurable MG-Inflam20 proteins.\n"
        "- C. Measurable MG-Inflam20 proteins placed in the complete paired proteome rank.\n"
        f"- Manifest `{manifest_path}`\n"
        f"- Legend `{legend_path}`\n"
    )
    result_path.write_text(result_text, encoding="utf-8")
    for path in (preview_pdf, preview_png, preview_gray, preview_layout, manifest_path, legend_path, result_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 4 output: {path.name}")
    print(f"Figure 4 panels: {len(records)}")


if __name__ == "__main__":
    main()
