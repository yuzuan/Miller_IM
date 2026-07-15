#!/usr/bin/env python3
"""Figure 4：PDC000514 配对蛋白组对 Miller-IM program 的正交验证。"""

from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

from recurrent_figure_style import apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP41 = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure4/source_data"
STEP39_SENSITIVITY = (
    ROOT
    / "write/39_independent_miller_mg_inflammatory_figures/Figure4/source_data/"
    "Figure4D_peptide_sensitivity_source.csv"
)
DATA_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
PDC_MATRIX = DATA_ROOT / (
    "write/31_proteomic_validation/source_metadata/PDC000514/"
    "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
)

WRITE_ROOT = ROOT / "write/54_figure4_proteome_rebuild/Figure4"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/54_figure4_proteome_rebuild/Figure4"

MILLER_IM_GENES = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]

PRIMARY = "#3C5488"
RECURRENT = "#E64B35"
TEXT = "#222222"
NEUTRAL = "#7A7A7A"
LIGHT = "#E7E7E7"
PALE_RED = "#F8DEDA"
PALE_BLUE = "#DCE5F1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    metadata = {
        "Creator": "Step54 recurrent GBM Figure 4",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step54 recurrent GBM Figure 4"},
    )
    plt.close(fig)
    return pdf, png


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
) -> dict[str, str]:
    return {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": sha256(source),
        "pdf_sha256": sha256(pdf),
        "png_sha256": sha256(png),
    }


def load_program_deltas() -> pd.DataFrame:
    path = STEP41 / "Figure4A_patient_delta_source.csv"
    table = pd.read_csv(path)
    if len(table) != 105:
        raise ValueError(f"PDC主分析必须为105对，当前为{len(table)}")
    if int((table["delta"] > 0).sum()) != 72:
        raise ValueError("PDC主分析必须为72/105复发升高")
    if not np.isclose(table["delta"].mean(), 0.252534746498684, atol=1e-12):
        raise ValueError("PDC患者层平均delta已改变")
    summary = pd.read_csv(STEP39_SENSITIVITY)
    summary = summary.loc[summary["analysis"].eq("Unique peptides")]
    if len(summary) != 1:
        raise ValueError("PDC unique-peptide统计摘要必须唯一")
    row = summary.iloc[0]
    if not (
        int(row["n_pairs"]) == 105
        and int(row["n_recurrence_up"]) == 72
        and np.isclose(float(row["mean_delta"]), table["delta"].mean(), atol=1e-12)
    ):
        raise ValueError("PDC患者层统计摘要与逐患者源表不一致")
    table["sign_flip_p"] = float(row["sign_flip_fdr"])
    table["mean_delta"] = float(row["mean_delta"])
    table["ci_low"] = float(row["ci_low"])
    table["ci_high"] = float(row["ci_high"])
    table["n_recurrence_up"] = int(row["n_recurrence_up"])
    table["n_pairs"] = int(row["n_pairs"])
    table = table.sort_values(["delta", "patient"], ascending=[True, True]).reset_index(drop=True)
    table["patient_rank"] = np.arange(1, len(table) + 1)
    table["direction"] = np.where(table["delta"] > 0, "Recurrent higher", "Primary higher")
    return table


def load_gene_statistics() -> pd.DataFrame:
    path = STEP41 / "Figure4B_measurable_protein_forest_source.csv"
    table = pd.read_csv(path)
    table["significant"] = table["adj.P.Val"] < 0.05
    if len(table) != 11 or int(table["significant"].sum()) != 5:
        raise ValueError("PDC Miller-IM蛋白必须为11个可测、5个全蛋白组FDR显著")
    expected = {"RHOB", "GSTM3", "FOLR2", "PDK4", "CTTNBP2"}
    observed = set(table.loc[table["significant"], "gene"])
    if observed != expected:
        raise ValueError(f"PDC显著Miller-IM蛋白集合已改变: {sorted(observed)}")
    return table


def rebuild_patient_gene_deltas(
    program_deltas: pd.DataFrame,
    gene_stats: pd.DataFrame,
) -> pd.DataFrame:
    if not PDC_MATRIX.exists():
        raise FileNotFoundError(f"PDC原始矩阵未挂载: {PDC_MATRIX}")

    header = pd.read_csv(PDC_MATRIX, sep="\t", nrows=0)
    usecols = ["Gene"] + [
        column for column in header.columns if column.startswith("Unshared Log KNCC_GBM")
    ]
    raw = pd.read_csv(PDC_MATRIX, sep="\t", usecols=usecols)
    raw = raw.loc[
        raw["Gene"].notna()
        & ~raw["Gene"].isin(["Mean", "Median", "StdDev"])
        & raw["Gene"].isin(MILLER_IM_GENES)
    ].drop_duplicates("Gene", keep="first")
    raw = raw.set_index("Gene")

    sample_pattern = re.compile(r"^Unshared Log (KNCC_GBM[0-9]+_T[12]):")
    source_columns: dict[str, list[str]] = {}
    for column in raw.columns:
        match = sample_pattern.match(column)
        if match:
            source_columns.setdefault(match.group(1), []).append(column)
    expression = pd.DataFrame(
        {sample: raw[columns].mean(axis=1) for sample, columns in source_columns.items()}
    )

    sample_rows = []
    for sample in expression.columns:
        patient = re.sub(r"_T[12]$", "", sample)
        condition = "Recurrent" if sample.endswith("_T2") else "Primary"
        sample_rows.append((sample, patient, condition))
    sample_map = pd.DataFrame(sample_rows, columns=["sample", "patient", "condition"])
    counts = sample_map.pivot_table(
        index="patient", columns="condition", values="sample", aggfunc="nunique", fill_value=0
    )
    paired = counts.index[(counts.get("Primary", 0) == 1) & (counts.get("Recurrent", 0) == 1)]
    if len(paired) != 105:
        raise ValueError(f"PDC原始矩阵必须形成105对，当前为{len(paired)}")

    records: list[dict[str, object]] = []
    for patient in sorted(paired):
        primary_sample = f"{patient}_T1"
        recurrent_sample = f"{patient}_T2"
        for gene in expression.index:
            primary_value = expression.at[gene, primary_sample]
            recurrent_value = expression.at[gene, recurrent_sample]
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
    deltas = pd.DataFrame(records)

    present = [gene for gene in MILLER_IM_GENES if gene in expression.index]
    paired_samples = [item for patient in paired for item in (f"{patient}_T1", f"{patient}_T2")]
    z_expression = expression.loc[present, paired_samples].sub(
        expression.loc[present, paired_samples].mean(axis=1), axis=0
    ).div(expression.loc[present, paired_samples].std(axis=1, ddof=1), axis=0)
    reconstructed_scores = []
    for patient in sorted(paired):
        primary = z_expression[f"{patient}_T1"]
        recurrent = z_expression[f"{patient}_T2"]
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
    reconstructed = pd.DataFrame(reconstructed_scores)
    audit = program_deltas[["patient", "delta", "n_common_raw20"]].merge(
        reconstructed, on="patient", how="outer", validate="one_to_one"
    )
    if audit.isna().any().any():
        raise ValueError("PDC患者层分数重建出现缺失")
    if np.max(np.abs(audit["delta"] - audit["reconstructed_delta"])) > 1e-10:
        raise ValueError("PDC原始矩阵重建与Step41患者delta不一致")
    if not np.array_equal(audit["n_common_raw20"], audit["reconstructed_n_common"]):
        raise ValueError("PDC每对共同可测蛋白数与Step41不一致")

    counts_per_gene = deltas.groupby("gene").size().rename("n_pairs_gene")
    up_per_gene = (
        deltas.assign(up=deltas["delta_log2_abundance"] > 0)
        .groupby("gene")["up"]
        .sum()
        .rename("n_recurrence_up_gene")
    )
    deltas = deltas.merge(counts_per_gene, on="gene").merge(up_per_gene, on="gene")
    deltas = deltas.merge(
        gene_stats[
            ["gene", "logFC", "ci_low", "ci_high", "P.Value", "adj.P.Val", "significant"]
        ],
        on="gene",
        how="left",
        validate="many_to_one",
    )
    observed_means = deltas.groupby("gene")["delta_log2_abundance"].mean()
    expected_means = gene_stats.set_index("gene")["logFC"]
    common = observed_means.index.intersection(expected_means.index)
    if np.max(np.abs(observed_means.loc[common] - expected_means.loc[common])) > 1e-10:
        raise ValueError("患者原始delta均值与paired limma logFC不一致")
    return deltas


def panel_a(program_deltas: pd.DataFrame) -> dict[str, str]:
    stem = "Figure4A_pdc_patient_delta_waterfall"
    source = SOURCE_OUT / f"{stem}_source.csv"
    program_deltas.to_csv(source, index=False)

    mean_delta = float(program_deltas["delta"].mean())
    se = float(program_deltas["delta"].std(ddof=1) / np.sqrt(len(program_deltas)))
    ci_low, ci_high = mean_delta - 1.96 * se, mean_delta + 1.96 * se

    fig = new_figure(88, 67)
    ax = fig.add_axes([0.13, 0.21, 0.82, 0.63])
    colors = np.where(program_deltas["delta"] > 0, RECURRENT, PRIMARY)
    ax.bar(
        program_deltas["patient_rank"],
        program_deltas["delta"],
        width=0.88,
        color=colors,
        edgecolor="none",
        alpha=0.88,
        zorder=2,
    )
    ax.axhline(0, color=TEXT, linewidth=0.72, zorder=3)
    mean_x = 108.2
    ax.plot([mean_x, mean_x], [ci_low, ci_high], color=TEXT, linewidth=1.0, zorder=4)
    ax.scatter([mean_x], [mean_delta], marker="D", s=28, color="white", edgecolor=TEXT, linewidth=0.8, zorder=5)
    ax.text(mean_x, ci_high + 0.08, "mean", ha="center", va="bottom", fontsize=5.8)
    ax.text(
        3,
        max(program_deltas["delta"]) * 0.92,
        "72/105 recurrent-higher",
        ha="left",
        va="top",
        fontsize=6.2,
        color=RECURRENT,
        fontweight="bold",
    )
    ax.set_xlim(0, 111)
    bound = max(abs(program_deltas["delta"].min()), abs(program_deltas["delta"].max())) * 1.08
    ax.set_ylim(-bound, bound)
    ax.set_xticks([])
    ax.set_xlabel("Paired patients ranked by change")
    ax.set_ylabel("Recurrent − primary Miller-IM protein score")
    clean_axis(ax, keep_left=True, keep_bottom=True)
    ax.tick_params(axis="x", length=0)
    fig.text(0.10, 0.955, "Miller-IM protein score rises at recurrence", fontsize=8.3, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record("A", stem, "105 paired-patient Miller-IM protein-score changes", source, pdf, png)


def panel_b(gene_deltas: pd.DataFrame, gene_stats: pd.DataFrame) -> dict[str, str]:
    stem = "Figure4B_raw20_protein_patient_distributions"
    source = SOURCE_OUT / f"{stem}_source.csv"
    plot_data = gene_deltas.copy()
    order = gene_stats.sort_values("logFC", ascending=False)["gene"].tolist()
    y_lookup = {gene: len(order) - 1 - index for index, gene in enumerate(order)}
    rng = np.random.default_rng(20260714)
    plot_data["y_base"] = plot_data["gene"].map(y_lookup).astype(float)
    plot_data["y_jitter"] = plot_data["y_base"] + rng.uniform(-0.20, 0.20, len(plot_data))
    plot_data.to_csv(source, index=False)

    fig = new_figure(178, 83)
    ax = fig.add_axes([0.115, 0.16, 0.855, 0.73])
    for index, gene in enumerate(order):
        yy = y_lookup[gene]
        if index % 2 == 0:
            ax.axhspan(yy - 0.46, yy + 0.46, color="#F7F7F7", zorder=0)
    positive = plot_data["delta_log2_abundance"] > 0
    ax.scatter(
        plot_data.loc[~positive, "delta_log2_abundance"],
        plot_data.loc[~positive, "y_jitter"],
        s=7,
        color=PRIMARY,
        alpha=0.35,
        edgecolor="none",
        zorder=2,
    )
    ax.scatter(
        plot_data.loc[positive, "delta_log2_abundance"],
        plot_data.loc[positive, "y_jitter"],
        s=7,
        color=RECURRENT,
        alpha=0.35,
        edgecolor="none",
        zorder=2,
    )
    stats = gene_stats.set_index("gene").loc[order].reset_index()
    stats["y"] = stats["gene"].map(y_lookup)
    ax.hlines(stats["y"], stats["ci_low"], stats["ci_high"], color=TEXT, linewidth=1.1, zorder=4)
    ax.scatter(
        stats["logFC"],
        stats["y"],
        marker="D",
        s=34,
        facecolor=np.where(stats["significant"], RECURRENT, "white"),
        edgecolor=TEXT,
        linewidth=0.75,
        zorder=5,
    )
    ax.axvline(0, color=TEXT, linewidth=0.7, zorder=1)
    ax.set_xlim(-2.75, 2.75)
    ax.set_ylim(-0.6, len(order) - 0.4)
    labels = []
    for gene in reversed(order):
        n_pairs = int(plot_data.loc[plot_data["gene"] == gene, "n_pairs_gene"].iloc[0])
        labels.append(f"{gene}  ({n_pairs})")
    ax.set_yticks(np.arange(len(order)), labels)
    ax.set_xlabel("Paired protein change (recurrent − primary)")
    ax.set_ylabel(None)
    ax.grid(axis="x", color=LIGHT, linewidth=0.45, zorder=0)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles = [
        mpl.lines.Line2D([], [], marker="D", linestyle="none", markerfacecolor=RECURRENT, markeredgecolor=TEXT, markersize=4.7, label="Full-proteome FDR < 0.05"),
        mpl.lines.Line2D([], [], marker="D", linestyle="none", markerfacecolor="white", markeredgecolor=TEXT, markersize=4.7, label="FDR ≥ 0.05"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", fontsize=5.7, handletextpad=0.4)
    fig.text(0.075, 0.96, "Five Miller-IM proteins increase significantly", fontsize=8.3, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record("B", stem, "Patient-level paired changes for 11 measurable Miller-IM proteins", source, pdf, png)


def panel_c(gene_stats: pd.DataFrame) -> dict[str, str]:
    stem = "Figure4C_full_proteome_rank_landscape"
    curve_path = STEP41 / "Figure4C_protein_ranked_gsea_curve_source.csv"
    curve = pd.read_csv(curve_path)
    plot_data = curve.merge(
        gene_stats[["gene", "logFC", "adj.P.Val", "significant"]],
        on="gene",
        how="left",
        validate="one_to_one",
    )
    plot_data["measured_miller_im"] = plot_data["logFC"].notna()
    summary = pd.read_csv(STEP39_SENSITIVITY)
    summary = summary.loc[summary["analysis"].eq("Unique peptides")]
    if len(summary) != 1:
        raise ValueError("PDC unique-peptide GSEA摘要必须唯一")
    plot_data["miller_im_nes"] = float(summary.iloc[0]["gsea_NES"])
    plot_data["miller_im_gsea_p"] = float(summary.iloc[0]["gsea_p"])
    source = SOURCE_OUT / f"{stem}_source.csv"
    plot_data.to_csv(source, index=False)

    fig = new_figure(88, 67)
    ax = fig.add_axes([0.13, 0.21, 0.83, 0.63])
    positive = curve["rank_stat"] >= 0
    ax.fill_between(
        curve.loc[positive, "rank"],
        0,
        curve.loc[positive, "rank_stat"],
        color=PALE_RED,
        linewidth=0,
        zorder=1,
    )
    ax.fill_between(
        curve.loc[~positive, "rank"],
        0,
        curve.loc[~positive, "rank_stat"],
        color=PALE_BLUE,
        linewidth=0,
        zorder=1,
    )
    ax.plot(curve["rank"], curve["rank_stat"], color="#9A9A9A", linewidth=0.55, zorder=2)
    measured = plot_data.loc[plot_data["measured_miller_im"]].copy()
    significant = measured["significant"].eq(True)
    ax.scatter(
        measured.loc[~significant, "rank"],
        measured.loc[~significant, "rank_stat"],
        s=28,
        facecolor="white",
        edgecolor=TEXT,
        linewidth=0.7,
        zorder=5,
    )
    ax.scatter(
        measured.loc[significant, "rank"],
        measured.loc[significant, "rank_stat"],
        s=32,
        facecolor=RECURRENT,
        edgecolor=TEXT,
        linewidth=0.7,
        zorder=6,
    )
    label_positions = {
        "RHOB": (700, 8.0, "left"),
        "GSTM3": (1350, 5.9, "left"),
        "FOLR2": (1850, 6.9, "left"),
        "PDK4": (2450, 3.8, "left"),
        "P2RY13": (7000, -3.8, "center"),
    }
    for row in measured.itertuples():
        if row.gene not in label_positions:
            continue
        label_x, label_y, ha = label_positions[row.gene]
        ax.annotate(
            row.gene,
            (row.rank, row.rank_stat),
            xytext=(label_x, label_y),
            textcoords="data",
            ha=ha,
            va="center",
            fontsize=5.9,
            color=TEXT,
            arrowprops={"arrowstyle": "-", "color": "#8A8A8A", "linewidth": 0.45},
        )
    ax.axhline(0, color=TEXT, linewidth=0.72, zorder=3)
    ax.text(0.01, 0.96, "recurrent-enriched", transform=ax.transAxes, color=RECURRENT, fontsize=6.2, fontweight="bold", ha="left", va="top")
    ax.text(0.99, 0.05, "primary-enriched", transform=ax.transAxes, color=PRIMARY, fontsize=6.2, fontweight="bold", ha="right", va="bottom")
    ax.text(0.99, 0.96, "11/20 proteins measured", transform=ax.transAxes, color=NEUTRAL, fontsize=5.8, ha="right", va="top")
    ax.set_xlim(curve["rank"].min(), curve["rank"].max())
    ax.set_xlabel("All proteins ranked by paired differential statistic")
    ax.set_ylabel("Moderated t statistic")
    ax.grid(axis="y", color=LIGHT, linewidth=0.45, zorder=0)
    clean_axis(ax, keep_left=True, keep_bottom=True)
    fig.text(0.10, 0.955, "Miller-IM proteins favor the recurrent proteome rank", fontsize=8.3, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record("C", stem, "All-protein rank landscape with 11 measurable Miller-IM proteins", source, pdf, png)


def image_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            try:
                return ImageFont.truetype(candidate, size=size)
            except OSError:
                pass
    return ImageFont.load_default()


def paste_contained(canvas: Image.Image, source: Path, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        fitted = ImageOps.contain(
            white.convert("RGB"),
            (right - left, bottom - top),
            Image.Resampling.LANCZOS,
        )
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 100, "top": 120, "right": 1870, "bottom": 1420},
        {"panel": "C", "left": 1930, "top": 120, "right": 3700, "bottom": 1420},
        {"panel": "B", "left": 100, "top": 1540, "right": 3700, "bottom": 2920},
    ]
    canvas = Image.new("RGB", (3800, 3000), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = image_font(58, bold=True)
    for item in layout:
        panel = item["panel"]
        paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 60), panel, fill="#111111", font=panel_font)
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
    preview.savefig(
        pdf,
        dpi=180,
        facecolor="white",
        edgecolor="none",
        metadata={"Creator": "Step54 recurrent GBM Figure 4", "CreationDate": None, "ModDate": None},
    )
    plt.close(preview)
    pd.DataFrame(layout).to_csv(source, index=False)
    return pdf, png, gray, source


def write_legend() -> Path:
    content = """# Figure 4 legend

**Figure 4 | Paired tumor proteomes support elevation of the measurable Miller-IM protein program at recurrence.**

**A,** Patient-level recurrent-minus-primary changes in the prespecified Miller-IM protein score in PDC000514. Unique-peptide protein abundances were used for the primary analysis. Each protein was standardized across the paired-sample matrix, and each patient score was calculated from proteins measured in both specimens of that pair; at least three common proteins were required. Eleven of 20 Miller-IM proteins were measurable globally, while individual pairs contained 7-11 common proteins. Seventy-two of 105 paired patients were recurrent-higher; the mean paired change was +0.253 (95% CI, +0.141 to +0.364; prespecified paired sign-flip P = 4.0 × 10^-5). The diamond and interval show the mean and 95% CI. **B,** Patient-level recurrent-minus-primary changes for the 11 measurable Miller-IM proteins. Small points denote individual patient pairs; parentheses after gene labels show the number of pairs with both measurements. Diamonds and horizontal intervals show paired limma effect estimates and approximate 95% CIs. Filled diamonds denote full-proteome Benjamini-Hochberg FDR < 0.05. RHOB, GSTM3, FOLR2, PDK4, and CTTNBP2 met this threshold. **C,** All 11,320 quantified proteins ranked by the paired moderated t statistic. The 11 measurable Miller-IM proteins are overlaid; filled points denote full-proteome FDR < 0.05. The prespecified Miller-IM gene set showed positive single-set competitive enrichment in the unique-peptide analysis (NES = 1.53; nominal P = 0.049), whereas the all-peptide sensitivity analysis was nonsignificant and is retained in the supplement.

PDC000514 is a whole-tumor TMT proteomic dataset rather than a myeloid-enriched assay. These results therefore provide orthogonal support for recurrence-associated elevation of the measurable protein program but do not identify the cellular source, prove an MCG1/MCG2 identity, or establish that every Miller-IM gene is translated and detected. Clinical metadata available to this analysis did not contain an IDH field, so the protein cohort is not labeled as IDH-wild-type.
"""
    path = WRITE_ROOT / "Figure4_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest: Path, legend: Path) -> None:
    content = f"""# Figure 4 proteome rebuild

- A: 105 paired-patient Miller-IM protein-score changes.
- B: patient-level distributions for all 11 measurable Miller-IM proteins.
- C: the 11 measurable proteins in the complete paired proteome rank.
- GeoMx-proteome correlation and the standard GSEA curve are not used as main panels.
- Manifest: `{manifest}`
- Legend: `{legend}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main() -> None:
    apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    program_deltas = load_program_deltas()
    gene_stats = load_gene_statistics()
    gene_deltas = rebuild_patient_gene_deltas(program_deltas, gene_stats)
    records = [panel_a(program_deltas), panel_b(gene_deltas, gene_stats), panel_c(gene_stats)]
    if [row["panel"] for row in records] != list("ABC"):
        raise ValueError("Figure4面板必须唯一锁定为A-C")
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_path = WRITE_ROOT / "Figure4_panel_manifest.csv"
    manifest = pd.DataFrame(records)
    manifest["preview_pdf"] = str(preview_pdf)
    manifest["preview_png"] = str(preview_png)
    manifest["preview_grayscale"] = str(preview_gray)
    manifest["preview_layout"] = str(preview_layout)
    manifest.to_csv(manifest_path, index=False)
    legend_path = write_legend()
    write_result(manifest_path, legend_path)
    print("STEP54_FIGURE4_COMPLETE panels=3 preview=1 mapping=A-C")


if __name__ == "__main__":
    main()
