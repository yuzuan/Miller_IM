#!/usr/bin/env python3
"""Figure 3 去单调重画：配对患者、基因指纹、单变量空间图谱与空间效应。"""

from __future__ import annotations

import gzip
import hashlib
import io
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

from recurrent_figure_style import apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP39 = ROOT / "write/39_independent_miller_mg_inflammatory_figures/Figure3/source_csv"
STEP41 = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data"
STEP51 = ROOT / "write/51_figure1_final_reorganized/Figure1/source_data"

EXTERNAL_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
GEOMX_ROOT = (
    EXTERNAL_ROOT
    / "write/26_pure_bioinformatics_dataset_rescue/source_metadata/GeoMx_zenodo16839828"
)
SPATIAL_IMAGE_ROOT = (
    EXTERNAL_ROOT
    / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"
)

WRITE_ROOT = ROOT / "write/53_figure3_visual_variety_rebuild/Figure3"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/53_figure3_visual_variety_rebuild/Figure3"

SPATIAL_PREFIX = {
    "GBM030": "GSM8506649_03422",
    "GBM049": "GSM8506650_25526",
}
SPATIAL_SOURCE = {
    "GBM030": STEP41 / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv",
    "GBM049": STEP41 / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv",
}

PRIMARY = "#3C5488"
RECURRENT = "#E64B35"
MDSC = "#159D88"
MES = "#7E6148"
TEXT = "#222222"
NEUTRAL = "#707070"
LIGHT = "#D8D8D8"
GRID = "#E7E7E7"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    metadata = {
        "Creator": "Step53 recurrent GBM Figure 3",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step53 recurrent GBM Figure 3"},
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


def patient_table() -> pd.DataFrame:
    table = pd.read_csv(STEP41 / "fig3a_geomx22_paired_change.csv")
    if len(table) != 22 or int((table["delta"] > 0).sum()) != 16:
        raise ValueError("GeoMx严格口径必须为22对且16/22复发升高")
    table = table.sort_values(["delta", "Patient_No"], ascending=[False, True]).reset_index(drop=True)
    table["patient_label"] = table["Patient_No"].map(lambda value: f"P{int(value):02d}")
    table["patient_order"] = np.arange(len(table))
    table["direction"] = np.where(table["delta"] > 0, "Recurrent higher", "Primary higher")
    return table


def shared_genes() -> list[str]:
    table = pd.read_csv(STEP51 / "Figure1_raw20_definition_source.csv")
    genes = (
        table.loc[table["shared_leading_edge"].astype(bool)]
        .sort_values("raw20_order")["gene"]
        .drop_duplicates()
        .tolist()
    )
    expected = ["PDK4", "SGK1", "CCL3", "CH25H", "SIGLEC8", "KLF6", "FOLR2", "CCL4"]
    if genes != expected:
        raise ValueError(f"Figure1共享leading-edge定义改变: {genes}")
    return genes


def panel_a(patients: pd.DataFrame) -> dict[str, str]:
    stem = "Figure3A_geomx_paired_dumbbell"
    source = SOURCE_OUT / f"{stem}_source.csv"
    patients.to_csv(source, index=False)

    fig = new_figure(92, 92)
    ax = fig.add_axes([0.22, 0.20, 0.72, 0.65])
    y = np.arange(len(patients), dtype=float)
    for yy, row in zip(y, patients.itertuples()):
        ax.plot(
            [row.Primary, row.Recurrence],
            [yy, yy],
            color="#A6A6A6",
            linewidth=1.15,
            alpha=0.72,
            zorder=1,
        )
    ax.scatter(
        patients["Primary"],
        y,
        s=22,
        marker="o",
        color=PRIMARY,
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
        label="Primary",
    )
    ax.scatter(
        patients["Recurrence"],
        y,
        s=22,
        marker="s",
        color=RECURRENT,
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
        label="Recurrence",
    )

    mean_y = len(patients) + 0.8
    primary_mean = float(patients["Primary"].mean())
    recurrent_mean = float(patients["Recurrence"].mean())
    ax.plot(
        [primary_mean, recurrent_mean],
        [mean_y, mean_y],
        color=TEXT,
        linewidth=1.5,
        zorder=2,
    )
    ax.scatter(primary_mean, mean_y, s=34, marker="o", color=PRIMARY, edgecolor="white", linewidth=0.6, zorder=4)
    ax.scatter(recurrent_mean, mean_y, s=34, marker="s", color=RECURRENT, edgecolor="white", linewidth=0.6, zorder=4)

    ax.set_yticks(
        list(y) + [mean_y],
        patients["patient_label"].tolist() + ["Mean"],
        fontsize=5.6,
    )
    ax.set_ylim(mean_y + 1.0, -1.0)
    values = patients[["Primary", "Recurrence"]].to_numpy(dtype=float)
    margin = max(0.12, float(np.ptp(values)) * 0.08)
    ax.set_xlim(float(np.nanmin(values) - margin), float(np.nanmax(values) + margin))
    ax.set_xlabel("Patient-level Miller-IM score")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.58, 0.035),
        frameon=False,
        fontsize=5.8,
        ncol=2,
        handletextpad=0.4,
        columnspacing=0.8,
    )

    fig.text(0.12, 0.955, "IBA1+ tissue recurrence validation", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record("A", stem, "Patient-level paired GeoMx Miller-IM scores", source, pdf, png)


def strict_gene_delta(patients: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    metadata_path = GEOMX_ROOT / "metadata.csv"
    expression_path = GEOMX_ROOT / "normalized.txt"
    if not metadata_path.exists() or not expression_path.exists():
        raise FileNotFoundError("GeoMx原始归一化矩阵未挂载")
    metadata = pd.read_csv(metadata_path)
    current = metadata.loc[
        metadata["iba1"].eq(True)
        & metadata["QCFlags"].isna()
        & metadata["IDH_status"].eq("IDH_WT")
        & metadata["tumor_setting"].isin(["Primary", "Recurrence"])
    ].copy()
    expression = pd.read_csv(expression_path, sep="\t", index_col=0)
    measured = [gene for gene in genes if gene in expression.index]
    missing = [gene for gene in genes if gene not in expression.index]
    if measured != genes[:-1] or missing != ["CCL4"]:
        raise ValueError(f"GeoMx shared-leading-edge覆盖改变: measured={measured}, missing={missing}")
    ids = current["Unnamed: 0"].tolist()
    log_expression = np.log2(expression.loc[measured, ids] + 1.0)
    patient_gene = (
        log_expression.T.assign(
            Patient_No=current.set_index("Unnamed: 0").loc[log_expression.columns, "Patient_No"].to_numpy(),
            tumor_setting=current.set_index("Unnamed: 0").loc[log_expression.columns, "tumor_setting"].to_numpy(),
        )
        .groupby(["Patient_No", "tumor_setting"], observed=True)[measured]
        .mean()
        .reset_index()
    )
    patient_ids = patients["Patient_No"].astype(int).tolist()
    rows: list[dict[str, object]] = []
    for gene_order, gene in enumerate(genes):
        if gene not in measured:
            for patient_order, patient_id in enumerate(patient_ids):
                rows.append(
                    {
                        "Patient_No": patient_id,
                        "patient_label": f"P{patient_id:02d}",
                        "patient_order": patient_order,
                        "gene": gene,
                        "gene_order": gene_order,
                        "measured": False,
                        "gene_delta": np.nan,
                    }
                )
            continue
        wide = patient_gene.pivot(index="Patient_No", columns="tumor_setting", values=gene)
        for patient_order, patient_id in enumerate(patient_ids):
            if patient_id not in wide.index or pd.isna(wide.loc[patient_id, ["Primary", "Recurrence"]]).any():
                raise ValueError(f"严格GeoMx患者缺少{gene}配对值: {patient_id}")
            rows.append(
                {
                    "Patient_No": patient_id,
                    "patient_label": f"P{patient_id:02d}",
                    "patient_order": patient_order,
                    "gene": gene,
                    "gene_order": gene_order,
                    "measured": True,
                    "gene_delta": float(wide.loc[patient_id, "Recurrence"] - wide.loc[patient_id, "Primary"]),
                }
            )
    table = pd.DataFrame(rows)
    program = patients.set_index("Patient_No")["delta"].to_dict()
    table["raw20_delta"] = table["Patient_No"].map(program)
    summary = pd.read_csv(STEP41 / "fig3b_geomx18_gene_forest.csv").set_index("gene")
    table["mean_gene_delta"] = table["gene"].map(summary["mean_log2_delta"])
    table["ci_low"] = table["gene"].map(summary["ci_low"])
    table["ci_high"] = table["gene"].map(summary["ci_high"])
    table["gene_fdr"] = table["gene"].map(summary["fdr"])
    table.loc[~table["measured"], ["mean_gene_delta", "ci_low", "ci_high", "gene_fdr"]] = np.nan
    measured_values = table.loc[table["measured"], "gene_delta"].abs().to_numpy(dtype=float)
    color_limit = float(np.quantile(measured_values, 0.98))
    table["color_limit"] = color_limit
    return table


def panel_b(patients: pd.DataFrame) -> dict[str, str]:
    stem = "Figure3B_shared_gene_patient_fingerprint"
    genes = shared_genes()
    table = strict_gene_delta(patients, genes)
    if len(table) != 176 or int(table["measured"].sum()) != 154:
        raise ValueError("B面板必须为22患者×8基因，且7基因可测")
    source = SOURCE_OUT / f"{stem}_source.csv"
    table.to_csv(source, index=False)

    matrix = (
        table.pivot(index="gene", columns="patient_label", values="gene_delta")
        .reindex(index=genes, columns=patients["patient_label"].tolist())
    )
    color_limit = float(table["color_limit"].iloc[0])
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#DADADA")
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))

    fig = new_figure(116, 92)
    ax = fig.add_axes([0.12, 0.28, 0.66, 0.58])
    image = ax.imshow(masked, cmap=cmap, vmin=-color_limit, vmax=color_limit, aspect="auto")
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=6.4)
    ax.set_xticks(
        np.arange(len(patients)),
        patients["patient_label"],
        rotation=90,
        fontsize=4.8,
    )
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(patients), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.text(
        len(patients) / 2 - 0.5,
        len(genes) - 1,
        "not measured",
        ha="center",
        va="center",
        fontsize=5.2,
        color="#777777",
    )

    color_ax = fig.add_axes([0.19, 0.16, 0.49, 0.022])
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    colorbar.set_label("Recurrent − primary log2 expression", fontsize=5.7, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=5.2, length=1.5, pad=1)
    colorbar.outline.set_visible(False)

    summary = table.drop_duplicates("gene").set_index("gene").reindex(genes)
    sx = fig.add_axes([0.81, 0.28, 0.16, 0.58])
    sx.axvline(0, color="#A0A0A0", linewidth=0.65)
    for index, gene in enumerate(genes):
        row = summary.loc[gene]
        if not bool(row["measured"]):
            sx.text(0.52, index, "NA", fontsize=5.0, color="#888888", va="center", ha="center")
            continue
        significant = float(row["gene_fdr"]) < 0.05
        sx.hlines(index, float(row["ci_low"]), float(row["ci_high"]), color=RECURRENT, linewidth=1.0)
        sx.scatter(
            float(row["mean_gene_delta"]),
            index,
            s=22,
            facecolor=RECURRENT if significant else "white",
            edgecolor=RECURRENT,
            linewidth=0.8,
            zorder=3,
        )
    sx.set_xlim(-0.45, 0.92)
    sx.set_ylim(len(genes) - 0.5, -0.5)
    sx.set_xticks([0, 0.5], ["0", "+0.5"], fontsize=4.8)
    sx.set_yticks([])
    sx.set_xlabel("Mean Δ", fontsize=5.5, labelpad=2)
    sx.grid(axis="x", color=GRID, linewidth=0.45)
    clean_axis(sx, keep_left=False, keep_bottom=True)

    fig.text(0.10, 0.955, "Shared leading-edge gene fingerprint", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.81,
        0.16,
        "Filled: gene-level\nFDR < 0.05",
        fontsize=5.1,
        color=NEUTRAL,
        va="top",
        linespacing=1.15,
    )
    pdf, png = save_panel(fig, stem)
    return record("B", stem, "Patient-by-gene GeoMx delta fingerprint", source, pdf, png)


def load_rgb_image(sample: str) -> np.ndarray:
    prefix = SPATIAL_PREFIX[sample]
    path = SPATIAL_IMAGE_ROOT / f"{prefix}_tissue_hires_image.png.gz"
    with gzip.open(path, "rb") as handle:
        image = Image.open(io.BytesIO(handle.read())).convert("RGB")
    return np.asarray(image, dtype=np.float32) / 255.0


def soften_tissue(rgb: np.ndarray) -> np.ndarray:
    gray = rgb.mean(axis=2, keepdims=True)
    muted = 0.55 * rgb + 0.45 * gray
    return np.clip(0.20 + 0.80 * muted, 0.0, 1.0)


def crop_limits(table: pd.DataFrame, image: np.ndarray) -> tuple[float, float, float, float]:
    x_min, x_max = table["x_hires"].min(), table["x_hires"].max()
    y_min, y_max = table["y_hires"].min(), table["y_hires"].max()
    x_margin = (x_max - x_min) * 0.10
    y_margin = (y_max - y_min) * 0.10
    return (
        max(0.0, float(x_min - x_margin)),
        min(float(image.shape[1]), float(x_max + x_margin)),
        max(0.0, float(y_min - y_margin)),
        min(float(image.shape[0]), float(y_max + y_margin)),
    )


def spatial_tables() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    tables = {sample: pd.read_csv(path) for sample, path in SPATIAL_SOURCE.items()}
    for sample, table in tables.items():
        if set(table["sample"]) != {sample}:
            raise ValueError(f"空间表样本名不符: {sample}")
    combined = pd.concat(tables.values(), ignore_index=True)
    if len(combined) != 7629:
        raise ValueError("Visium主表必须包含7,629个spots")
    for variable in ["raw20", "mdsc_combined", "mes"]:
        low, high = combined[variable].quantile([0.02, 0.98]).tolist()
        combined[f"{variable}_scale_low"] = float(low)
        combined[f"{variable}_scale_high"] = float(high)
    split = {
        sample: combined.loc[combined["sample"].eq(sample)].copy()
        for sample in ["GBM030", "GBM049"]
    }
    return split, combined


def panel_c() -> dict[str, str]:
    stem = "Figure3C_spatial_score_atlas"
    tables, combined = spatial_tables()
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    combined.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    columns = [
        ("raw20", "Miller-IM", "Miller-IM score"),
        ("mdsc_combined", "MDSC-like", "MDSC-like score"),
        ("mes", "MES-like", "MES-like score"),
    ]
    scales = {
        variable: (
            float(combined[f"{variable}_scale_low"].iloc[0]),
            float(combined[f"{variable}_scale_high"].iloc[0]),
        )
        for variable, _, _ in columns
    }

    fig = new_figure(142, 96)
    grid = fig.add_gridspec(
        2,
        3,
        left=0.10,
        right=0.98,
        bottom=0.18,
        top=0.88,
        wspace=0.05,
        hspace=0.08,
    )
    cmap = mpl.colormaps["viridis"]
    axes: dict[tuple[int, int], mpl.axes.Axes] = {}
    images = {sample: soften_tissue(load_rgb_image(sample)) for sample in tables}
    for row_index, sample in enumerate(["GBM030", "GBM049"]):
        table = tables[sample]
        tissue = images[sample]
        x0, x1, y0, y1 = crop_limits(table, tissue)
        for col_index, (variable, title, _) in enumerate(columns):
            ax = fig.add_subplot(grid[row_index, col_index])
            axes[(row_index, col_index)] = ax
            low, high = scales[variable]
            norm = mcolors.Normalize(vmin=low, vmax=high, clip=True)
            ax.imshow(tissue, origin="upper", interpolation="nearest")
            ax.scatter(
                table["x_hires"],
                table["y_hires"],
                c=table[variable],
                cmap=cmap,
                norm=norm,
                s=8.5,
                alpha=0.94,
                linewidths=0,
                rasterized=True,
                zorder=3,
            )
            ax.set_xlim(x0, x1)
            ax.set_ylim(y1, y0)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row_index == 0:
                ax.set_title(title, fontsize=7.4, fontweight="bold", pad=4)
            if col_index == 0:
                ax.text(
                    -0.03,
                    0.98,
                    f"{sample}\n{len(table):,} spots",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=6.0,
                    color=TEXT,
                    linespacing=1.15,
                )

    for col_index, (variable, _, colorbar_label) in enumerate(columns):
        low, high = scales[variable]
        position = axes[(1, col_index)].get_position()
        color_ax = fig.add_axes([position.x0 + 0.02, 0.105, position.width - 0.04, 0.018])
        scalar = mpl.cm.ScalarMappable(norm=mcolors.Normalize(vmin=low, vmax=high), cmap=cmap)
        colorbar = fig.colorbar(scalar, cax=color_ax, orientation="horizontal")
        colorbar.set_ticks([low, (low + high) / 2, high])
        colorbar.ax.tick_params(labelsize=5.6, length=1.5, pad=1)
        colorbar.outline.set_visible(False)
        colorbar.set_label(colorbar_label, fontsize=5.8, labelpad=1.5)

    fig.text(0.10, 0.955, "Spatial score maps in two primary GBMs", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record(
        "C",
        stem,
        "Two-slice maps of Miller-IM, MDSC-like, and MES-like spatial scores",
        source,
        pdf,
        png,
    )


def panel_d() -> dict[str, str]:
    stem = "Figure3D_partial_rho_effect_plot"
    original = pd.read_csv(STEP41 / "fig3e_partial_rho_main.csv")
    table = original[["sample", "geo_accession", "target", "n_spots", "partial_rho"]].copy()
    table["inference_boundary"] = "descriptive_only_no_spatial_autocorrelation_correction"
    if len(table) != 4 or not (table["partial_rho"] > 0).all():
        raise ValueError("D面板必须包含两切片×MDSC/MES四个正向效应")
    source = SOURCE_OUT / f"{stem}_source.csv"
    table.to_csv(source, index=False)

    fig = new_figure(58, 88)
    ax = fig.add_axes([0.28, 0.20, 0.68, 0.66])
    target_y = {"MDSC": 1.0, "MES": 0.0}
    sample_offset = {"GBM030": 0.09, "GBM049": -0.09}
    sample_marker = {"GBM030": "o", "GBM049": "s"}
    target_color = {"MDSC": MDSC, "MES": MES}
    for target in ["MDSC", "MES"]:
        current = table.loc[table["target"].eq(target)].set_index("sample")
        x_values = [float(current.loc[sample, "partial_rho"]) for sample in ["GBM030", "GBM049"]]
        y_values = [target_y[target] + sample_offset[sample] for sample in ["GBM030", "GBM049"]]
        for sample, xx, yy in zip(["GBM030", "GBM049"], x_values, y_values):
            ax.scatter(
                xx,
                yy,
                s=34,
                marker=sample_marker[sample],
                facecolor=target_color[target],
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
                label=sample,
            )
            ax.text(xx + 0.012, yy, f"{xx:.3f}", fontsize=5.7, va="center", color=TEXT)
    ax.axvline(0, color="#999999", linewidth=0.7)
    ax.set_xlim(0, 0.50)
    ax.set_ylim(-0.48, 1.48)
    ax.set_yticks([1, 0], ["MDSC-like", "MES-like"], fontsize=6.4)
    ax.set_xlabel("Partial Spearman r")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="none", color=TEXT, markerfacecolor="white", label="GBM030"),
        mpl.lines.Line2D([], [], marker="s", linestyle="none", color=TEXT, markerfacecolor="white", label="GBM049"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=5.4, loc="lower right", handletextpad=0.3)
    fig.text(0.15, 0.955, "Myeloid-adjusted niche association", fontsize=8.2, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return record("D", stem, "Descriptive partial-rho effects in two Visium tumors", source, pdf, png)


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
        {"panel": "A", "left": 100, "top": 235, "right": 1660, "bottom": 1360},
        {"panel": "B", "left": 1740, "top": 235, "right": 3700, "bottom": 1360},
        {"panel": "C", "left": 100, "top": 1630, "right": 2770, "bottom": 2820},
        {"panel": "D", "left": 2870, "top": 1630, "right": 3700, "bottom": 2820},
    ]
    canvas = Image.new("RGB", (3800, 2920), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = image_font(58, bold=True)
    section_font = image_font(36, bold=True)

    drawer.text((100, 55), "Independent recurrence validation", fill=RECURRENT, font=section_font)
    drawer.line((100, 118, 3700, 118), fill="#E8B2AA", width=4)
    drawer.text(
        (100, 1450),
        "Inflammatory-microglial program in MDSC/MES-rich niches",
        fill=MDSC,
        font=section_font,
    )
    drawer.line((100, 1515, 3700, 1515), fill="#A8D9CF", width=4)

    for item in layout:
        panel = item["panel"]
        paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 62), panel, fill="#111111", font=panel_font)
    png = FIG_OUT / "Figure3_visual_variety_preview.png"
    gray = FIG_OUT / "Figure3_visual_variety_preview_grayscale.png"
    pdf = FIG_OUT / "Figure3_visual_variety_preview.pdf"
    source = SOURCE_OUT / "Figure3_visual_variety_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(180, 180), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(180, 180), optimize=True)
    canvas.save(pdf, format="PDF", resolution=180.0)
    pd.DataFrame(layout).to_csv(source, index=False)
    return pdf, png, gray, source


def write_legend() -> Path:
    content = """# Figure 3 legend

**Figure 3 | Independent myeloid-enriched tissue validation and spatial niche context of the Miller-IM program.**

**A,** Patient-level paired Miller-IM scores in the independent Artzi et al. GeoMx cohort. Quality-controlled IBA1+ areas of illumination were averaged within each patient and time point before calculating the Miller-IM score from the 18 measurable genes. Lines connect primary and recurrent specimens from the same patient; circle and square markers denote primary and recurrent samples, respectively. Sixteen of 22 paired IDH-wild-type patients were recurrent-higher; the mean paired change was +0.273 (95% CI, +0.094 to +0.452; BH-adjusted sign-flip FDR across three prespecified GeoMx entry definitions = 0.00721). **B,** Patient-level recurrent-minus-primary expression changes for the eight shared leading-edge genes defined in Figure 1, displayed in the same patient order as panel A. Seven genes were measured; CCL4 was not measured. Heatmap colors were clipped symmetrically at the 98th percentile of the absolute patient-gene changes; 4 of 154 measured cells exceeded the display limits. The right-hand summary shows mean paired changes and 95% CIs; filled symbols denote gene-level FDR < 0.05 after correction across all 18 measurable Miller-IM genes. This panel is a gene-level fingerprint and does not replace the Miller-IM score based on 18 measurable genes used in panel A. **C,** H&E-backed Visium maps from two independent untreated primary IDH-wild-type GBMs. Rows denote tumors and columns show the Miller-IM, combined MDSC-like, and MES-like scores. Nineteen of 20 Miller-IM genes were measurable; AC253572.2 was not measured. For each score, the color scale is shared across the two tumors and clipped at the pooled 2nd and 98th percentiles for display. These maps describe multicellular neighborhood context and do not assign MDSC identity to Miller-IM-high spots. **D,** Partial Spearman correlations between the Miller-IM score and combined MDSC-like or MES-like scores after within-tumor rank residualization for total myeloid score. The two tumors are shown as independent points and are not connected. Only effect sizes are shown because the analysis contains two biological tumors and the spot-level tests did not correct for spatial autocorrelation.

Panel A provides independent tissue-level recurrence validation with patient as the statistical unit; panel B decomposes this effect across prespecified leading-edge genes. Figure 2 localized the Miller-IM program preferentially to MCG1/MCG2-like inflammatory-microglial states rather than MDSC-like transcriptional states. Here, panels C-D show that Miller-IM-high Visium spots have positive myeloid-adjusted associations with MDSC/MES-like spatial scores. Cell-state identity and multicellular spatial niche are therefore distinct analytical levels: an inflammatory-microglial program can occur within an MDSC/MES-rich neighborhood. Panels C-D are cross-sectional and do not test recurrence. Treatment-stratified GeoMx analyses, the complete 18-gene forest, and E-/M-MDSC sensitivity analyses are retained as supplementary analyses.
"""
    path = WRITE_ROOT / "Figure3_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest: Path, legend: Path) -> None:
    content = f"""# Figure 3 visual-variety rebuild

- A: 22-patient paired GeoMx dumbbell.
- B: patient-by-shared-leading-edge gene delta fingerprint, with the same patient order as A.
- C: two tumors by three single-score H&E-backed spatial maps.
- D: four descriptive myeloid-adjusted spatial effects.
- Protein and anti-PD-1 results remain outside Figure 3.
- Manifest: `{manifest}`
- Legend: `{legend}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main() -> None:
    apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    if not SPATIAL_IMAGE_ROOT.exists():
        raise FileNotFoundError(f"空间组织图目录未挂载: {SPATIAL_IMAGE_ROOT}")

    patients = patient_table()
    records = [panel_a(patients), panel_b(patients), panel_c(), panel_d()]
    if [row["panel"] for row in records] != list("ABCD"):
        raise ValueError("Figure3面板必须唯一锁定为A-D")
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_path = WRITE_ROOT / "Figure3_panel_manifest.csv"
    manifest = pd.DataFrame(records)
    manifest["preview_pdf"] = str(preview_pdf)
    manifest["preview_png"] = str(preview_png)
    manifest["preview_grayscale"] = str(preview_gray)
    manifest["preview_layout"] = str(preview_layout)
    manifest.to_csv(manifest_path, index=False)
    legend_path = write_legend()
    write_result(manifest_path, legend_path)
    print("STEP53_FIGURE3_COMPLETE panels=4 preview=1 mapping=A-D")


if __name__ == "__main__":
    main()
