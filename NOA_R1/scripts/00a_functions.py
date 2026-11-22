#!/usr/bin/env python3
"""Shared paths, figure styling, and output utilities for the NOA R1 scripts."""

from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image, ImageFont, ImageOps


def legacy_root(revision_root: Path) -> Path:
    """Return the upstream analysis root, which defaults to the repository parent."""
    configured = os.environ.get("NOA_LEGACY_ROOT")
    root = Path(configured).expanduser() if configured else revision_root.parent
    root = root.resolve()
    return root


def public_data_root(project_root: Path) -> Path:
    """Return the configured root of downloaded public source data."""
    configured = os.environ.get("NOA_PUBLIC_DATA_ROOT") or os.environ.get("MILLER_IM_DATA_ROOT")
    root = Path(configured).expanduser() if configured else project_root / "data"
    root = root.resolve()
    return root


def mm(value: float) -> float:
    """Convert a measurement in millimeters to inches."""
    inches = value / 25.4
    return inches


# Keep the palette used for the original paired recurrence figures.
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
    """Set the existing display font, line widths, and white background. Returns nothing."""
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
    """Return a white Figure with the specified dimensions in millimeters."""
    figure = plt.figure(figsize=(mm(width_mm), mm(height_mm)), facecolor="white", **kwargs)
    return figure


def clean_axis(ax, *, keep_left: bool = True, keep_bottom: bool = True) -> None:
    """Set border and tick appearance for one Axes. Returns nothing."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)
    ax.tick_params(width=0.6, length=2.5)


def render_pdf_png(pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """Rasterize a single-page PDF and check the resulting PNG. Returns nothing."""
    if not pdf_path.is_file() or pdf_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PDF: {pdf_path.name}")

    # Use Poppler so the raster output keeps the PDF master layout.
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["pdftoppm", "-singlefile", "-png", "-r", str(dpi), str(pdf_path), str(png_path.with_suffix(""))],
        check  = True,
        stdout = subprocess.DEVNULL,
    )

    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PNG: {png_path.name}")
    with Image.open(png_path) as image:
        if min(image.size) <= 0:
            raise ValueError(f"Invalid PNG dimensions: {png_path.name}")
        image.verify()


def save_panel(fig, out_dir: Path, stem: str, creator: str) -> tuple[Path, Path]:
    """Write a Figure as PDF and PNG and return the two output paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{stem}.pdf"
    fig.canvas.draw()
    fig.savefig(
        pdf_path,
        facecolor = "white",
        edgecolor = "none",
        metadata  = {"Creator": creator, "CreationDate": None, "ModDate": None},
    )

    png_path = out_dir / f"{stem}.png"
    render_pdf_png(pdf_path, png_path)
    plt.close(fig)
    output_paths = (pdf_path, png_path)
    return output_paths


def save_figure(fig, out_dir: Path, stem: str) -> None:
    """Save a panel with the shared Step 41 metadata. Returns nothing."""
    save_panel(fig, out_dir, stem, "Step41 recurrent GBM figure rebuild")


def sha256(path: Path) -> str:
    """Return the file digest for the existing panel manifest schema."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    checksum = digest.hexdigest()
    return checksum


def image_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return the preview font using the original Arial and Helvetica fallback order."""
    font_name = "Arial Bold.ttf" if bold else "Arial.ttf"
    candidates = [Path("/System/Library/Fonts/Supplemental") / font_name, Path("/System/Library/Fonts/Helvetica.ttc")]
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
    """Fit a panel inside its preview rectangle and center it. Returns nothing."""
    left, top, right, bottom = box
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        fitted = ImageOps.contain(white.convert("RGB"), (right - left, bottom - top), Image.Resampling.LANCZOS)

    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))
