#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps
from scipy.stats import spearmanr

from recurrent_figure_style import apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP45_SOURCE = ROOT / "write/45_figure1_program_rebuild/Figure1/source_data"
WRITE_ROOT = ROOT / "write/46_figure1_ef_candidates/Figure1"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/46_figure1_ef_candidates/Figure1"

DATASETS = ["GSE174554", "GSE274546"]
DATASET_COLORS = {"GSE174554": "#3C5488", "GSE274546": "#00A087"}
HIGHLIGHT = "#E64B35"
RAW20_NAME = "Miller_Microglial_Inflammatory_raw_top20"
EXPECTED_SHARED_LEADING_EDGE = {"CCL4", "CH25H", "FOLR2", "KLF6", "PDK4", "SGK1", "SIGLEC8"}


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
        "Creator": "Step46 recurrent GBM Figure 1 E-F candidates",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step46 recurrent GBM Figure 1 E-F candidates"},
    )
    plt.close(fig)
    return pdf, png


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, str]:
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


def panel_e_rank_rank() -> dict[str, str]:
    stem = "Figure1E_cross_cohort_gene_rank_comparison"
    source = SOURCE_OUT / f"{stem}_source.csv"
    input_path = STEP45_SOURCE / "Figure1D_cross_cohort_all_gene_rank_density_source.csv"
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
        raise ValueError(f"Missing Figure1E source columns: {sorted(required - set(frame.columns))}")
    frame.to_csv(source, index=False)

    fig = new_figure(89.5, 82)
    ax = fig.add_axes([0.17, 0.16, 0.73, 0.70])
    background = frame.loc[~frame["miller_raw20"]].sort_values("density")
    raw20_other = frame.loc[frame["miller_raw20"] & ~frame["shared_leading_edge"]]
    shared = frame.loc[frame["shared_leading_edge"]]
    if len(shared) != 7 or len(raw20_other) != 12:
        raise ValueError(f"Expected 7 shared and 12 other measured raw20 genes, found {len(shared)} and {len(raw20_other)}")
    if set(shared["gene"]) != EXPECTED_SHARED_LEADING_EDGE:
        raise ValueError(
            "Unexpected shared leading-edge genes: "
            f"expected {sorted(EXPECTED_SHARED_LEADING_EDGE)}, found {sorted(shared['gene'])}"
        )

    ax.scatter(
        background["rank_stat_GSE174554"],
        background["rank_stat_GSE274546"],
        c=background["density"],
        cmap="Greys",
        s=2.0,
        alpha=0.68,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        raw20_other["rank_stat_GSE174554"],
        raw20_other["rank_stat_GSE274546"],
        s=22,
        facecolor="white",
        edgecolor=HIGHLIGHT,
        linewidth=0.85,
        zorder=3,
    )
    ax.scatter(
        shared["rank_stat_GSE174554"],
        shared["rank_stat_GSE274546"],
        s=28,
        color=HIGHLIGHT,
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )

    label_offsets = {
        "CCL4": (-25, 7),
        "CH25H": (5, 5),
        "SGK1": (-23, -10),
        "FOLR2": (5, 7),
        "PDK4": (5, -10),
    }
    for _, row in shared.loc[shared["gene"].isin(label_offsets)].iterrows():
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=label_offsets[row["gene"]],
            textcoords="offset points",
            fontsize=5.7,
            color="#B53A2B",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
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

    rho, _ = spearmanr(frame["rank_stat_GSE174554"], frame["rank_stat_GSE274546"])
    ax.text(
        0.02,
        0.98,
        f"Shared genes = {len(frame):,}\nRank correlation $r_s$ = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.4,
        color="#686868",
    )
    ax.text(
        0.98,
        0.03,
        "7 shared leading-edge genes\n12 other measured raw20 genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.0,
        color="#B53A2B",
        linespacing=1.15,
    )
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
    fig.text(0.17, 0.96, "Cross-cohort gene-rank comparison", fontsize=8.6, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    return panel_record(
        "Figure1E",
        stem,
        "All-gene cross-cohort rank-stat comparison with shared leading-edge genes highlighted",
        source,
        pdf,
        png,
    )


def panel_f_lopo() -> tuple[dict[str, str], Path]:
    stem = "Figure1F_leave_one_patient_out_gsea_stability"
    source = SOURCE_OUT / "Figure1F_leave_one_patient_out_raw20_gsea.csv"
    if not source.exists():
        raise FileNotFoundError(f"Run 46_fig1_lopo_gsea.R first: {source}")
    frame = pd.read_csv(source)
    if set(frame["pathway"]) != {RAW20_NAME}:
        raise ValueError("Figure1F source must contain only the fixed Miller raw20 program")
    full = frame.loc[frame["run_type"].eq("Full")].set_index("dataset")
    lopo = frame.loc[frame["run_type"].eq("Leave-one-patient-out")].copy()
    expected_counts = {"GSE174554": 17, "GSE274546": 45}
    if lopo.groupby("dataset").size().to_dict() != expected_counts:
        raise ValueError(f"Unexpected LOPO counts: {lopo.groupby('dataset').size().to_dict()}")

    summary_rows = []
    for dataset in DATASETS:
        part = lopo.loc[lopo["dataset"].eq(dataset)]
        summary_rows.append(
            {
                "dataset": dataset,
                "n_omissions": len(part),
                "full_NES": float(full.loc[dataset, "NES"]),
                "full_nominal_P": float(full.loc[dataset, "nominal_P"]),
                "LOPO_NES_min": float(part["NES"].min()),
                "LOPO_NES_max": float(part["NES"].max()),
                "positive_NES_folds": int((part["NES"] > 0).sum()),
                "nominal_P_lt_0_05_folds": int((part["nominal_P"] < 0.05).sum()),
                "max_nominal_P": float(part["nominal_P"].max()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary_path = SOURCE_OUT / "Figure1F_leave_one_patient_out_summary.csv"
    summary.to_csv(summary_path, index=False)

    fig = new_figure(89.5, 82)
    ax = fig.add_axes([0.24, 0.21, 0.71, 0.65])
    y_positions = {"GSE174554": 1.0, "GSE274546": 0.0}
    rng = np.random.default_rng(20260713)
    all_values = []
    for dataset in DATASETS:
        part = lopo.loc[lopo["dataset"].eq(dataset)].sort_values(["NES", "omitted_patient"]).copy()
        y = y_positions[dataset]
        color = DATASET_COLORS[dataset]
        all_values.extend(part["NES"].tolist())
        ax.hlines(y, part["NES"].min(), part["NES"].max(), color=color, linewidth=2.3, alpha=0.42, zorder=1)
        jitter_grid = np.linspace(-0.13, 0.13, len(part))
        jitter = jitter_grid[rng.permutation(len(part))]
        significant = part["nominal_P"].to_numpy(dtype=float) < 0.05
        ax.scatter(
            part.loc[significant, "NES"],
            y + jitter[significant],
            s=17,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=0.92,
            zorder=3,
        )
        ax.scatter(
            part.loc[~significant, "NES"],
            y + jitter[~significant],
            s=18,
            facecolor="white",
            edgecolor=color,
            linewidth=0.85,
            alpha=1.0,
            zorder=3,
        )
        ax.scatter(
            [float(full.loc[dataset, "NES"])],
            [y],
            marker="D",
            s=48,
            color="#222222",
            edgecolor="white",
            linewidth=0.65,
            zorder=5,
        )

        row = summary.loc[summary["dataset"].eq(dataset)].iloc[0]
        ax.text(
            0.985,
            y + 0.23,
            f"NES {row['LOPO_NES_min']:.2f}–{row['LOPO_NES_max']:.2f}   "
            f"positive {int(row['positive_NES_folds'])}/{int(row['n_omissions'])}\n"
            f"P<0.05 {int(row['nominal_P_lt_0_05_folds'])}/{int(row['n_omissions'])}   "
            f"max P {row['max_nominal_P']:.3g}",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=6.0,
            color="#555555",
            linespacing=1.15,
        )

    x_max = max(all_values + full["NES"].tolist())
    ax.set_xlim(-0.05, x_max + 0.18)
    ax.set_ylim(-0.48, 1.48)
    ax.axvline(0, color="#999999", linestyle=(0, (3, 2)), linewidth=0.7)
    ax.set_yticks(
        [y_positions[dataset] for dataset in DATASETS],
        [f"{dataset}\n{expected_counts[dataset]} omissions" for dataset in DATASETS],
        fontsize=6.9,
    )
    for label, dataset in zip(ax.get_yticklabels(), DATASETS):
        label.set_color(DATASET_COLORS[dataset])
        label.set_fontweight("bold")
    ax.set_xlabel("Leave-one-patient-out normalized enrichment score")
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.45)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=5)

    fig.text(0.16, 0.96, "Leave-one-patient-out GSEA stability", fontsize=8.6, fontweight="bold", va="top")
    fig.text(
        0.16,
        0.075,
        "Each circle: one omission   Filled: nominal P < 0.05   Diamond: full data",
        fontsize=5.8,
        color="#666666",
        va="bottom",
    )
    pdf, png = save_panel(fig, stem)
    return (
        panel_record(
            "Figure1F",
            stem,
            "Leave-one-patient-out paired pseudobulk raw20 GSEA stability",
            source,
            pdf,
            png,
        ),
        summary_path,
    )


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            try:
                return ImageFont.truetype(str(path), size=size)
            except OSError:
                continue
    return ImageFont.load_default()


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path]:
    paths = {record["panel"]: Path(record["png"]) for record in records}
    layout_rows = [
        {"panel": "Figure1E", "left": 120, "top": 170, "right": 1750, "bottom": 1770},
        {"panel": "Figure1F", "left": 1850, "top": 170, "right": 3480, "bottom": 1770},
    ]
    canvas = Image.new("RGB", (3600, 1900), "white")
    drawer = ImageDraw.Draw(canvas)
    label_font = font(54, bold=True)
    for row in layout_rows:
        panel = row["panel"]
        with Image.open(paths[panel]) as image:
            image = image.convert("RGB")
            width = int(row["right"] - row["left"])
            height = int(row["bottom"] - row["top"])
            fitted = ImageOps.contain(image, (width, height), Image.Resampling.LANCZOS)
            x = int(row["left"] + (width - fitted.width) / 2)
            y = int(row["top"] + (height - fitted.height) / 2)
            canvas.paste(fitted, (x, y))
        drawer.text((int(row["left"]), 70), panel[-1], fill="#111111", font=label_font)

    png = FIG_OUT / "Figure1EF_candidate_preview.png"
    pdf = FIG_OUT / "Figure1EF_candidate_preview.pdf"
    source = SOURCE_OUT / "Figure1EF_candidate_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(300, 300), optimize=True)
    preview_fig = plt.figure(figsize=(12, 1900 / 300), facecolor="white")
    preview_ax = preview_fig.add_axes([0, 0, 1, 1])
    preview_ax.imshow(np.asarray(canvas), interpolation="nearest")
    preview_ax.axis("off")
    preview_fig.savefig(
        pdf,
        dpi=300,
        facecolor="white",
        edgecolor="none",
        metadata={
            "Creator": "Step46 recurrent GBM Figure 1 E-F candidates",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(preview_fig)
    pd.DataFrame(layout_rows).to_csv(source, index=False)
    return pdf, png, source


def main() -> None:
    np.random.seed(20260713)
    apply_publication_style()
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    records = [panel_e_rank_rank()]
    panel_f, summary_path = panel_f_lopo()
    records.append(panel_f)
    preview_pdf, preview_png, preview_source = make_preview(records)
    records.append(
        panel_record(
            "Figure1EF_preview",
            "Figure1EF_candidate_preview",
            "Review-only preview of Figure 1E-F candidates",
            preview_source,
            preview_pdf,
            preview_png,
        )
    )
    manifest = pd.DataFrame(records)
    manifest["lopo_summary"] = ""
    manifest.loc[manifest["panel"].eq("Figure1F"), "lopo_summary"] = str(summary_path)
    manifest.to_csv(WRITE_ROOT / "Figure1EF_panel_manifest.csv", index=False)
    print("STEP46_FIGURE1_EF_COMPLETE panels=2 preview=1")


if __name__ == "__main__":
    main()
