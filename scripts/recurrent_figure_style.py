#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


MM_PER_INCH = 25.4

COLORS = {
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


def mm(value: float) -> float:
    return value / MM_PER_INCH


def apply_publication_style() -> None:
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
            "axes.edgecolor": COLORS["text"],
            "axes.labelcolor": COLORS["text"],
            "text.color": COLORS["text"],
            "xtick.color": COLORS["text"],
            "ytick.color": COLORS["text"],
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
    return plt.figure(figsize=(mm(width_mm), mm(height_mm)), facecolor="white", **kwargs)


def clean_axis(ax, *, keep_left: bool = True, keep_bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)
    ax.tick_params(width=0.6, length=2.5)


def save_figure(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = out_dir / f"{stem}.pdf"
    png = out_dir / f"{stem}.png"
    fig.savefig(
        pdf,
        facecolor="white",
        edgecolor="none",
        metadata={"Creator": "Step41 recurrent GBM figure rebuild", "CreationDate": None, "ModDate": None},
    )
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step41 recurrent GBM figure rebuild"},
    )
    plt.close(fig)

