#!/usr/bin/env python3
# ======================================================================
# 00a_functions.py — Miller-IM shared figure functions
# Public paired recurrence analysis in glioblastoma.
#
# Purpose
# Share the plotting palette, PDF and PNG export, and figure-preview layout.
# Use the published Matplotlib settings and panel geometry.
# The analysis scripts provide examples at each function call.
# Use go-to-definition in the editor to follow a call to its implementation.
# Function signatures describe the inputs used by the figure scripts.

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageFont, ImageOps


def mm(value: float) -> float:
    """Convert millimeters to inches using the exact factor 25.4."""
    inches = value / 25.4
    return inches


# These colors match the labels in the published figure panels.
colors = {
    "primary": "#4C78A8",
    "recurrent": "#E45756",
    "neutral": "#777777",
    "light_neutral": "#D9D9D9",
    "grid": "#E8E8E8",
    "text": "#222222",
    "microglia": "#009E73",
    "macrophage": "#0072B2",
    "monocyte": "#E69F00",
    "cdc": "#CC79A7",
    "neutrophil": "#7A7A7A",
    "gse174554": "#3C5488",
    "gse274546": "#00A087",
}


def apply_publication_style() -> None:
    """Apply the published font, line width, and background settings."""
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.65,
            "axes.edgecolor": colors["text"],
            "axes.labelcolor": colors["text"],
            "text.color": colors["text"],
            "xtick.color": colors["text"],
            "ytick.color": colors["text"],
            "lines.linewidth": 1.15,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "figure.dpi": 120,
            "savefig.dpi": 600,
        }
    )


def new_figure(width_mm: float, height_mm: float, **kwargs):
    """Return a white Matplotlib figure with dimensions specified in millimeters."""
    figure = plt.figure(figsize=(mm(width_mm), mm(height_mm)), facecolor="white", **kwargs)
    return figure


def clean_axis(ax, *, keep_left: bool = True, keep_bottom: bool = True) -> None:
    """Apply the requested spine visibility and axis styling."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)
    ax.tick_params(width=0.6, length=2.5)


def render_pdf_png(pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """Render a single-page PDF as a PNG and verify that the image can be opened."""
    if not pdf_path.is_file() or pdf_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PDF: {pdf_path.name}")

    # Poppler reads the saved PDF at the requested resolution, normally 600 dpi.
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["pdftoppm", "-singlefile", "-png", "-r", str(dpi), str(pdf_path), str(png_path.with_suffix(""))],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    # Open the PNG to confirm that the renderer produced a readable image.
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PNG: {png_path.name}")
    with Image.open(png_path) as image:
        if min(image.size) <= 0:
            raise ValueError(f"Invalid PNG dimensions: {png_path.name}")
        image.verify()


def save_panel(fig, out_dir: Path, stem: str, creator: str) -> tuple[Path, Path]:
    """Save a figure with the original creator metadata and return the PDF and PNG paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{stem}.pdf"

    # Draw the figure before saving its PDF with the original size and background.
    fig.canvas.draw()
    fig.savefig(
        pdf_path,
        facecolor="white",
        edgecolor="none",
        metadata={"Creator": creator, "CreationDate": None, "ModDate": None},
    )

    # Render the PNG from the saved PDF so both exports use the same drawing.
    png_path = out_dir / f"{stem}.png"
    render_pdf_png(pdf_path, png_path)
    plt.close(fig)
    output_paths = (pdf_path, png_path)
    return output_paths


def save_figure(fig, out_dir: Path, stem: str) -> None:
    """Save the PDF and PNG with the original Step 41 creator metadata."""
    save_panel(fig, out_dir, stem, "Step41 recurrent GBM figure rebuild")


def sha256(path: Path) -> str:
    """Return the file digest used in the published source manifest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    checksum = digest.hexdigest()
    return checksum


def panel_record(panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path) -> dict[str, str]:
    """Return the nine manifest fields for a panel and its three output paths."""
    record = {
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
    return record


def image_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return the first available preview font from the original candidate order."""
    font_name = "Arial Bold.ttf" if bold else "Arial.ttf"
    candidates = [Path("/System/Library/Fonts/Supplemental") / font_name, Path("/System/Library/Fonts/Helvetica.ttc")]

    # Try the configured Arial and Helvetica locations before the Pillow default.
    for candidate in candidates:
        if not candidate.exists():
            continue
        try:
            font = ImageFont.truetype(str(candidate), size=size)
            return font
        except OSError:
            continue
    font = ImageFont.load_default()
    return font


def paste_contained(canvas: Image.Image, source: Path, box: tuple[int, int, int, int]) -> None:
    """Scale a panel proportionally and center it in the requested preview rectangle."""
    left, top, right, bottom = box
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        fitted = ImageOps.contain(white.convert("RGB"), (right - left, bottom - top), Image.Resampling.LANCZOS)

    # Divide spare width and height evenly around the centered panel.
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def enrichment_curve(deg_table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Build the ranked gene table with hit flags and running enrichment scores."""

    # Keep the original sort order and join keys when comparing the cohort gene ranks.
    ranked_table = deg_table.copy()
    ranked_table["rank_stat"] = np.sign(ranked_table["logFC"]) * np.sqrt(np.maximum(ranked_table["F"], 0))
    ranked_table = (
        ranked_table.loc[np.isfinite(ranked_table["rank_stat"])]
        .sort_values("rank_stat", ascending=False)
        .reset_index(drop=True)
    )
    if ranked_table["gene"].duplicated().any():
        duplicates = ranked_table.loc[ranked_table["gene"].duplicated(), "gene"].tolist()
        raise ValueError(f"Duplicate genes in Step38 ranked table: {duplicates[:10]}")
    hits = ranked_table["gene"].isin(genes).to_numpy()

    weights = np.abs(ranked_table["rank_stat"].to_numpy(dtype=float))
    hit_weights = np.where(hits, weights, 0.0)
    if hit_weights.sum() <= 0 or (~hits).sum() == 0:
        raise ValueError("Cannot construct enrichment curve from the Step38 ranked table")
    increments = np.where(hits, hit_weights / hit_weights.sum(), -1.0 / (~hits).sum())
    result_table = pd.DataFrame(
        {
            "rank": np.arange(1, len(ranked_table) + 1),
            "gene": ranked_table["gene"].to_numpy(),
            "rank_stat": ranked_table["rank_stat"].to_numpy(),
            "hit": hits,
            "running_enrichment": np.cumsum(increments),
        }
    )

    return result_table


def rank_density(rank_table: pd.DataFrame, bins: int = 75) -> np.ndarray:
    """Return the log1p bin count for each row in the two-cohort rank table."""

    x = rank_table["rank_stat_GSE174554"].to_numpy(dtype=float)
    y = rank_table["rank_stat_GSE274546"].to_numpy(dtype=float)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)

    # Index the two-dimensional counts at each row's assigned bin.
    result = np.log1p(counts[xi, yi])
    return result
