#!/usr/bin/env python3
"""锁定 Figure 1 的 A-F 映射并生成最终独立面板库。"""

from __future__ import annotations

import hashlib
import importlib.util
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw, ImageFont, ImageOps

from recurrent_figure_style import apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP38 = ROOT / "write/38_independent_cohort_mg_inflammatory_recalculation"
STEP45_FIG = ROOT / "figures/45_figure1_program_rebuild/Figure1"
STEP45_SOURCE = ROOT / "write/45_figure1_program_rebuild/Figure1/source_data"
STEP46_FIG = ROOT / "figures/46_figure1_ef_candidates/Figure1"
STEP46_SOURCE = ROOT / "write/46_figure1_ef_candidates/Figure1/source_data"

WRITE_ROOT = ROOT / "write/51_figure1_final_reorganized/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/51_figure1_final_reorganized/Figure1"

DATASETS = ["GSE174554", "GSE274546"]
PAIR_COUNTS = {"GSE174554": 18, "GSE274546": 45}
HIGHLIGHT = "#E64B35"
RAW20_NAME = "Miller_Microglial_Inflammatory_raw_top20"


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
        "Creator": "Step51 final recurrent GBM Figure 1",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step51 final recurrent GBM Figure 1"},
    )
    plt.close(fig)
    return pdf, png


def flatten_png(source: Path, destination: Path) -> None:
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        white.convert("RGB").save(destination, format="PNG", dpi=(600, 600), optimize=True)


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
    provenance: str,
) -> dict[str, str]:
    return {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "provenance": provenance,
        "source_sha256": sha256(source),
        "pdf_sha256": sha256(pdf),
        "png_sha256": sha256(png),
    }


def copy_locked_panel(
    panel: str,
    old_stem: str,
    new_stem: str,
    message: str,
    figure_root: Path,
    source_root: Path,
    source_suffix: str,
) -> dict[str, str]:
    source_input = source_root / f"{old_stem}_source{source_suffix}"
    source = SOURCE_OUT / f"{new_stem}_source{source_suffix}"
    pdf = FIG_OUT / f"{new_stem}.pdf"
    png = FIG_OUT / f"{new_stem}.png"
    shutil.copy2(source_input, source)
    shutil.copy2(figure_root / f"{old_stem}.pdf", pdf)
    flatten_png(figure_root / f"{old_stem}.png", png)
    return record(
        panel,
        new_stem,
        message,
        source,
        pdf,
        png,
        f"locked copy of {old_stem}",
    )


def panel_a() -> dict[str, str]:
    """从原始队列流程源表重绘A，避免在旧PDF上覆盖术语。"""
    module_path = ROOT / "scripts/42_fig1_fig2_visual_story.py"
    spec = importlib.util.spec_from_file_location("step42_visual_story", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法载入Figure1A绘图模块: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIG_ROOT = FIG_OUT.parent
    module.WRITE_ROOT = WRITE_ROOT.parent
    module.configure_style()
    generated = module.figure1_design()

    old_stem = "Figure1A_compact_paired_design"
    stem = "Figure1A_paired_cohort_design"
    old_source = SOURCE_OUT / f"{old_stem}_source.csv"
    old_pdf = FIG_OUT / f"{old_stem}.pdf"
    old_png = FIG_OUT / f"{old_stem}.png"
    source = SOURCE_OUT / f"{stem}_source.csv"
    pdf = FIG_OUT / f"{stem}.pdf"
    png = FIG_OUT / f"{stem}.png"
    old_source.replace(source)
    old_pdf.replace(pdf)
    flatten_png(old_png, png)
    old_png.unlink()
    apply_publication_style()
    return record(
        "A",
        stem,
        "Two independent paired recurrence cohorts and patient-level Miller-IM analysis design",
        source,
        pdf,
        png,
        f"vector redraw from {generated['stem']}",
    )


def panel_d() -> dict[str, str]:
    stem = "Figure1D_shared_raw20_outliers"
    input_path = STEP46_SOURCE / "Figure1E_cross_cohort_gene_rank_comparison_source.csv"
    frame = pd.read_csv(input_path)
    required = {
        "gene",
        "rank_stat_GSE174554",
        "rank_stat_GSE274546",
        "miller_raw20",
        "density",
        "shared_leading_edge",
    }
    if not required.issubset(frame.columns):
        raise ValueError(f"Figure1D 缺少列: {sorted(required - set(frame.columns))}")

    shared = frame.loc[frame["shared_leading_edge"].astype(bool)].copy()
    raw20_other = frame.loc[
        frame["miller_raw20"].astype(bool) & ~frame["shared_leading_edge"].astype(bool)
    ].copy()
    background = frame.loc[~frame["miller_raw20"].astype(bool)].copy()
    if len(frame) != 9432 or len(shared) != 8 or len(raw20_other) != 11:
        raise ValueError(
            f"Figure1D 数量不符: all={len(frame)}, shared={len(shared)}, other_raw20={len(raw20_other)}"
        )

    labelled_genes = {"CCL3", "CCL4", "CH25H", "FOLR2", "SGK1"}
    frame["display_role"] = "transcriptome_background"
    frame.loc[frame["miller_raw20"].astype(bool), "display_role"] = "other_measured_raw20"
    frame.loc[frame["shared_leading_edge"].astype(bool), "display_role"] = "shared_leading_edge"
    frame["display_label"] = frame["gene"].where(frame["gene"].isin(labelled_genes), "")
    source = SOURCE_OUT / f"{stem}_source.csv"
    frame.to_csv(source, index=False)

    fig = new_figure(91, 84)
    ax = fig.add_axes([0.17, 0.23, 0.75, 0.57])
    background = background.sort_values("density")
    ax.scatter(
        background["rank_stat_GSE174554"],
        background["rank_stat_GSE274546"],
        c=background["density"],
        cmap="Greys",
        s=2.0,
        alpha=0.66,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        raw20_other["rank_stat_GSE174554"],
        raw20_other["rank_stat_GSE274546"],
        s=23,
        facecolor="white",
        edgecolor=HIGHLIGHT,
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared["rank_stat_GSE174554"],
        shared["rank_stat_GSE274546"],
        s=29,
        color=HIGHLIGHT,
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )

    label_offsets = {
        "CCL3": (4, 7),
        "CCL4": (-25, 7),
        "CH25H": (5, 5),
        "SGK1": (-23, -10),
        "FOLR2": (5, 7),
    }
    for _, row in shared.loc[shared["gene"].isin(labelled_genes)].iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=label_offsets[row["gene"]],
            textcoords="offset points",
            fontsize=5.8,
            color="#B53A2B",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
            zorder=5,
        )

    ax.axhline(0, color="#A4A4A4", linewidth=0.55)
    ax.axvline(0, color="#A4A4A4", linewidth=0.55)
    limit = float(
        np.nanmax(np.abs(frame[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy()))
    ) * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E7E7E7", linewidth=0.42, zorder=0)
    clean_axis(ax)
    ax.text(
        0.98,
        0.98,
        "Positive = recurrent-enriched",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color="#777777",
    )
    fig.text(
        0.17,
        0.965,
        "Shared Miller-IM leading-edge genes are\nrecurrent-enriched outliers",
        fontsize=8.3,
        fontweight="bold",
        va="top",
        linespacing=1.05,
    )
    fig.text(
        0.17,
        0.025,
        "Filled: shared GSEA leading edge   Open: other measured Miller-IM genes",
        fontsize=5.8,
        color="#666666",
        va="bottom",
    )

    pdf, png = save_panel(fig, stem)
    return record(
        "D",
        stem,
        "Shared Miller-IM leading-edge genes occupy the joint recurrent-enriched region",
        source,
        pdf,
        png,
        "redrawn from Step46 all-gene rank source; no naked genome-wide correlation annotation",
    )


def panel_e() -> dict[str, str]:
    stem = "Figure1E_shared_leading_edge_gene_effects"
    input_path = STEP45_SOURCE / "Figure1E_shared_leading_edge_heatmap_source.csv"
    direction = pd.read_csv(input_path)
    if len(direction) != 16 or direction["gene"].nunique() != 8:
        raise ValueError("Figure1E 必须是 8 个基因 × 2 个队列")
    if int((direction["FDR"] < 0.05).sum()) != 0:
        raise ValueError("Figure1E 的正式单基因 FDR 结果与既有结论不符")
    source = SOURCE_OUT / f"{stem}_source.csv"
    direction.to_csv(source, index=False)

    order = (
        direction.sort_values("gene_order")["gene"].drop_duplicates().tolist()
    )
    matrix = (
        direction.pivot(index="gene", columns="dataset", values="logFC")
        .reindex(index=order, columns=DATASETS)
    )
    max_abs = float(np.nanmax(np.abs(matrix.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    fig = new_figure(86, 80)
    ax = fig.add_axes([0.24, 0.19, 0.55, 0.60])
    image = ax.imshow(
        matrix.to_numpy(),
        cmap="RdBu_r",
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_xticks(
        range(len(DATASETS)),
        [f"GSE174554\n{PAIR_COUNTS['GSE174554']} pairs", f"GSE274546\n{PAIR_COUNTS['GSE274546']} pairs"],
        fontsize=6.8,
    )
    ax.set_yticks(range(len(order)), order, fontsize=7.1)
    ax.tick_params(length=0)
    for row_index, gene in enumerate(order):
        for col_index, dataset in enumerate(DATASETS):
            value = float(matrix.loc[gene, dataset])
            color = "white" if abs(value) > 0.9 else "#222222"
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=6.4, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(DATASETS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    color_ax = fig.add_axes([0.82, 0.30, 0.025, 0.38])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Recurrent − primary\nlog2 fold change", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6.2, length=2)
    colorbar.outline.set_linewidth(0.45)
    fig.text(0.24, 0.95, "Shared leading-edge gene effects", fontsize=8.5, fontweight="bold", va="top")
    fig.text(
        0.24,
        0.085,
        "0/16 gene-by-cohort tests reached FDR < 0.05",
        fontsize=5.9,
        color="#777777",
        va="bottom",
    )
    pdf, png = save_panel(fig, stem)
    return record(
        "E",
        stem,
        "Effect sizes of the eight shared leading-edge genes in both cohorts",
        source,
        pdf,
        png,
        "redrawn from Step45 gene-level effect source; redundant cohort GSEA text removed",
    )


def copy_panel_f() -> dict[str, str]:
    old_stem = "Figure1F_leave_one_patient_out_gsea_stability"
    new_stem = old_stem
    source_input = STEP46_SOURCE / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    source = SOURCE_OUT / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    shutil.copy2(source_input, source)
    for name in [
        "Figure1F_leave_one_patient_out_summary.csv",
        "Figure1F_full_data_reproduction_check.csv",
        "Figure1F_leave_one_patient_out_fold_audit.csv",
        "Figure1F_leave_one_patient_out_program_gsea.csv",
    ]:
        shutil.copy2(STEP46_SOURCE / name, SOURCE_OUT / name)
    pdf = FIG_OUT / f"{new_stem}.pdf"
    png = FIG_OUT / f"{new_stem}.png"
    shutil.copy2(STEP46_FIG / f"{old_stem}.pdf", pdf)
    flatten_png(STEP46_FIG / f"{old_stem}.png", png)
    return record(
        "F",
        new_stem,
        "Leave-one-patient-out paired pseudobulk Miller-IM GSEA stability",
        source,
        pdf,
        png,
        "locked copy of Step46 Figure1F",
    )


def build_raw20_definition() -> Path:
    provenance = pd.read_csv(STEP38 / "fixed_program_provenance.csv")
    raw20_row = provenance.loc[provenance["signature"].eq(RAW20_NAME)].iloc[0]
    genes = str(raw20_row["genes"]).split(";")
    if len(genes) != 20:
        raise ValueError(f"raw20 定义不是 20 个基因: {len(genes)}")

    direction = pd.read_csv(STEP38 / "independent_fixed_program_gene_direction.csv")
    direction = direction.loc[
        direction["dataset"].isin(DATASETS)
        & direction["threshold"].eq(20)
        & direction["signature"].eq(RAW20_NAME)
    ].copy()
    gsea = pd.read_csv(STEP38 / "independent_fixed_program_targeted_gsea.csv")
    gsea = gsea.loc[
        gsea["dataset"].isin(DATASETS)
        & gsea["threshold"].eq(20)
        & gsea["pathway"].eq(RAW20_NAME)
    ].copy()
    leading = {
        row["dataset"]: set(str(row["leadingEdge"]).split(";"))
        for _, row in gsea.iterrows()
    }
    shared = set.intersection(*(leading[dataset] for dataset in DATASETS))

    rows: list[dict[str, object]] = []
    direction_index = direction.set_index(["dataset", "gene"])
    for order, gene in enumerate(genes, start=1):
        for dataset in DATASETS:
            row = direction_index.loc[(dataset, gene)]
            rows.append(
                {
                    "raw20_order": order,
                    "gene": gene,
                    "dataset": dataset,
                    "tested": bool(row["tested"]),
                    "logFC": row["logFC"],
                    "PValue": row["PValue"],
                    "FDR": row["FDR"],
                    "cohort_leading_edge": gene in leading[dataset],
                    "shared_leading_edge": gene in shared,
                }
            )
    table = pd.DataFrame(rows)
    if len(table) != 40:
        raise ValueError("raw20 定义表必须有 40 行")
    tested_counts = table.groupby("dataset")["tested"].sum().astype(int).to_dict()
    if tested_counts != {"GSE174554": 19, "GSE274546": 19}:
        raise ValueError(f"raw20 可测数量不符: {tested_counts}")
    if set(table.loc[~table["tested"], "gene"]) != {"AC253572.2"}:
        raise ValueError("raw20 唯一缺失基因应为 AC253572.2")
    if int(table["shared_leading_edge"].sum()) != 16:
        raise ValueError("共享 leading-edge 应为 8 基因 × 2 队列")
    path = SOURCE_OUT / "Figure1_raw20_definition_source.csv"
    table.to_csv(path, index=False)
    return path


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
        {"panel": "A", "left": 180, "top": 120, "right": 3420, "bottom": 820},
        {"panel": "B", "left": 160, "top": 940, "right": 1700, "bottom": 1940},
        {"panel": "C", "left": 1900, "top": 940, "right": 3440, "bottom": 1940},
        {"panel": "D", "left": 80, "top": 2110, "right": 1160, "bottom": 3230},
        {"panel": "E", "left": 1260, "top": 2110, "right": 2340, "bottom": 3230},
        {"panel": "F", "left": 2440, "top": 2110, "right": 3520, "bottom": 3230},
    ]
    canvas = Image.new("RGB", (3600, 3320), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = image_font(58, bold=True)
    for item in layout:
        panel = item["panel"]
        paste_contained(
            canvas,
            by_panel[panel],
            (item["left"], item["top"], item["right"], item["bottom"]),
        )
        drawer.text((item["left"], item["top"] - 60), panel, fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1_final_preview.png"
    gray = FIG_OUT / "Figure1_final_preview_grayscale.png"
    pdf = FIG_OUT / "Figure1_final_preview.pdf"
    layout_path = SOURCE_OUT / "Figure1_final_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(180, 180), optimize=True)
    ImageOps.grayscale(canvas).save(gray, format="PNG", dpi=(180, 180), optimize=True)
    canvas.save(pdf, format="PDF", resolution=180.0)
    pd.DataFrame(layout).to_csv(layout_path, index=False)
    return pdf, png, gray, layout_path


def write_legend(raw20_path: Path) -> Path:
    raw20 = pd.read_csv(raw20_path).sort_values("raw20_order")["gene"].drop_duplicates().tolist()
    raw20_text = ", ".join(raw20)
    content = f"""# Figure 1 legend

**Figure 1 | The 20-gene Miller-derived inflammatory microglial program (Miller-IM program) is enriched at recurrence in two independent paired GBM cohorts.**

**A,** Study design. GSE174554 and GSE274546 were reconstructed independently, pan-myeloid cells were reclustered with patient-level Harmony correction, and recurrence was tested using paired patient-level analyses (18 and 45 matched patients, respectively). **B-C,** Preranked GSEA of the fixed Miller-IM program using paired raw-count pseudobulk differential-expression ranks. Positive normalized enrichment scores indicate recurrent enrichment. GSE174554: NES = 2.32, nominal *P* < 1 × 10⁻⁴; GSE274546: NES = 1.80, nominal *P* = 0.00398. **D,** Genome-wide paired rank statistics in the two cohorts. Filled red symbols mark the eight measurable Miller-IM genes present in both cohort-specific GSEA leading edges; open red symbols mark the other measurable Miller-IM genes. Only CCL3, CCL4, CH25H, FOLR2, and SGK1 are labelled to avoid over-annotation. Because filled-symbol membership is defined by the intersection of the two leading edges, their joint positive-quadrant position is descriptive and partly expected by construction; this panel is not an independent single-gene replication test. Genome-wide rank correlation was low (Spearman rₛ = 0.16) and is provided here only as background context, not as the panel claim. **E,** Recurrent-versus-primary log₂ fold changes of the eight shared leading-edge genes. None of the 16 gene-by-cohort tests reached FDR < 0.05, supporting interpretation at the program rather than single-gene level. **F,** Leave-one-patient-out sensitivity analysis of the paired pseudobulk Miller-IM GSEA. Normalized enrichment scores remained positive and nominal *P* remained < 0.05 for all 18/18 and 45/45 omissions; the maximum nominal *P* was 0.0311 in GSE274546. Diamonds show full-data estimates; filled circles denote nominal *P* < 0.05.

The fixed Miller-IM gene set was defined a priori from Miller et al. Supplementary Table 2: {raw20_text}. Nineteen of 20 genes were measurable in each cohort; AC253572.2 was not measured. Statistical unit: patient. Panel D is a descriptive cross-cohort background view; panels B, C, and F provide the formal program-level recurrence evidence.
"""
    path = WRITE_ROOT / "Figure1_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest_path: Path, legend_path: Path, raw20_path: Path) -> None:
    content = f"""# Figure 1 final reorganization

- Final panel mapping: A design; B-C independent Miller-IM GSEA; D all-gene background with shared leading-edge Miller-IM genes; E shared leading-edge gene effects; F leave-one-patient-out stability.
- D contains no naked genome-wide correlation annotation and labels only CCL3, CCL4, CH25H, FOLR2, and SGK1.
- D and E are unique files; the historical Step46 duplicate letter assignment is not carried forward.
- All six independent panels are letter-free. Letters occur only in the review preview.
- Manifest: `{manifest_path}`
- Legend: `{legend_path}`
- Internal raw20 version definition: `{raw20_path}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main() -> None:
    apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    records = [
        panel_a(),
        copy_locked_panel(
            "B",
            "Figure1B_GSE174554_raw20_gsea",
            "Figure1B_GSE174554_raw20_gsea",
            "Formal paired raw-count pseudobulk Miller-IM GSEA in GSE174554",
            STEP45_FIG,
            STEP45_SOURCE,
            ".csv.gz",
        ),
        copy_locked_panel(
            "C",
            "Figure1C_GSE274546_raw20_gsea",
            "Figure1C_GSE274546_raw20_gsea",
            "Formal paired raw-count pseudobulk Miller-IM GSEA in GSE274546",
            STEP45_FIG,
            STEP45_SOURCE,
            ".csv.gz",
        ),
        panel_d(),
        panel_e(),
        copy_panel_f(),
    ]
    if [row["panel"] for row in records] != list("ABCDEF"):
        raise ValueError("Figure1 面板映射必须唯一锁定为 A-F")
    if len({row["stem"] for row in records}) != 6:
        raise ValueError("Figure1 存在重复文件名")

    raw20_path = build_raw20_definition()
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_path = WRITE_ROOT / "Figure1_panel_manifest.csv"
    manifest = pd.DataFrame(records)
    manifest["preview_pdf"] = str(preview_pdf)
    manifest["preview_png"] = str(preview_png)
    manifest["preview_grayscale"] = str(preview_gray)
    manifest["preview_layout"] = str(preview_layout)
    manifest.to_csv(manifest_path, index=False)
    legend_path = write_legend(raw20_path)
    write_result(manifest_path, legend_path, raw20_path)
    print("STEP51_FIGURE1_COMPLETE panels=6 preview=1 mapping=A-F")


if __name__ == "__main__":
    main()
