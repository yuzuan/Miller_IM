#!/usr/bin/env python3
"""Supply the common Figure 3 panel functions and style for 08c_figure3.py.

Inputs: The prepared GeoMx paired scores and the shared leading-edge genes from Figure 1 under NOA_LEGACY_ROOT.
Outputs: The panel A figure and source table in the directories set by the Figure 3 renderer.
Run: python NOA_R1/scripts/09b_export_final_figures.py --figures 3
"""

from __future__ import annotations

import importlib
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd


figure_tools = importlib.import_module("00a_functions")
apply_publication_style = figure_tools.apply_publication_style
clean_axis = figure_tools.clean_axis
new_figure = figure_tools.new_figure
revision_root = Path(__file__).resolve().parents[1]
project_root = figure_tools.legacy_root(revision_root)

# The final exporter redirects these Figure 3 output paths.
WRITE_ROOT = revision_root / "_figure_build/55_figure3_visual_polish/Figure3"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = revision_root / "_panels/55_figure3_visual_polish/Figure3"

# Figure 3 uses these original colors across panels, previews, and the R1 renderer.
primary = "#3C5488"
recurrent = "#E64B35"
mdsc_color = "#159D88"
mes_color = "#7E6148"
text_color = "#222222"
grid_color = "#E7E7E7"


def save_panel(fig: mpl.figure.Figure, stem: str) -> tuple[Path, Path]:
    """Write the PDF master and its PNG. The Figure 3 R1 renderer replaces this binding."""
    output_paths = figure_tools.save_panel(fig, FIG_OUT, stem, "Step55 recurrent GBM Figure 3 visual polish")
    return output_paths


def record(
    panel: str,
    stem: str,
    message: str,
    source: Path,
    pdf: Path,
    png: Path,
) -> dict[str, str]:
    """Return the legacy manifest row for one panel and its source table."""
    for path in (source, pdf, png):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty Figure 3 panel output: {path.name}")
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


def patient_table() -> pd.DataFrame:
    """Return one row per IDH-wild-type GeoMx pair in the original change order."""
    path = project_root / "write/41_mg_inflammatory_sci_rebuild/Figure3/source_data/fig3a_geomx22_paired_change.csv"
    paired_table = pd.read_csv(path)
    if len(paired_table) != 22 or int((paired_table["delta"] > 0).sum()) != 16:
        raise ValueError("The frozen GeoMx cohort must contain 22 pairs with 16 recurrent-higher scores")
    paired_table = paired_table.sort_values(["delta", "Patient_No"], ascending=[False, True]).reset_index(drop=True)
    paired_table["patient_label"] = [f"P{int(patient):02d}" for patient in paired_table["Patient_No"]]
    paired_table["patient_order"] = np.arange(len(paired_table))
    paired_table["direction"] = np.where(paired_table["delta"] > 0, "Recurrent higher", "Primary higher")
    return paired_table


def panel_a(paired_table: pd.DataFrame) -> dict[str, str]:
    """Write and plot the GeoMx patients grouped by trial exposure and ordered by change."""
    stem = "Figure3A_geomx_paired_dumbbell"

    # Separate exposure groups without changing the time-point color encoding.
    display_table = paired_table.copy()
    group_counts = display_table["trial_setting"].value_counts().to_dict()
    if group_counts != {"Nivolumab": 12, "Control": 10}:
        raise ValueError(f"Unexpected trial group counts for Figure 3A: {group_counts}")
    display_table["display_group"] = display_table["trial_setting"].map({"Control": 0, "Nivolumab": 1})
    display_table = display_table.sort_values(["display_group", "delta", "Patient_No"], ascending=[True, False, True])
    source = SOURCE_OUT / f"{stem}_source.csv"
    display_table.to_csv(source, index=False)

    # The fixed y positions preserve the gap and labels for the two frozen exposure groups.
    fig = new_figure(92, 92)
    ax = fig.add_axes([0.30, 0.20, 0.64, 0.69])
    y = np.arange(len(display_table), dtype=float) + display_table["display_group"].to_numpy()
    ax.axhspan(10.5, 22.5, color="#F3F4F6", zorder=0)
    ax.axhline(10.0, color="#D0D0D0", linewidth=0.6)
    for row_y, group in [(4.5, "Control"), (16.5, "Nivolumab")]:
        ax.text(
            -0.36,
            row_y,
            f"{group}\n(n = {group_counts[group]})",
            transform = ax.get_yaxis_transform(),
            rotation  = 90,
            ha        = "center",
            va        = "center",
            fontsize  = 7.4,
            color     = text_color,
        )

    # Connect primary and recurrent scores for each patient and mark each time point.
    for row_y, row in zip(y, display_table.itertuples()):
        ax.plot(
            [row.Primary, row.Recurrence],
            [row_y, row_y],
            color     = "#A6A6A6",
            linewidth = 1.15,
            alpha     = 0.72,
            zorder    = 1,
        )
    ax.scatter(
        display_table["Primary"],
        y,
        s         = 22,
        marker    = "o",
        color     = primary,
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Primary",
    )
    ax.scatter(
        display_table["Recurrence"],
        y,
        s         = 22,
        marker    = "s",
        color     = recurrent,
        edgecolor = "white",
        linewidth = 0.55,
        zorder    = 3,
        label     = "Recurrence",
    )

    # The additional pooled-mean row uses the same time-point markers.
    mean_y = len(display_table) + 2.2
    primary_mean = float(display_table["Primary"].mean())
    recurrent_mean = float(display_table["Recurrence"].mean())
    ax.plot(
        [primary_mean, recurrent_mean],
        [mean_y, mean_y],
        color     = text_color,
        linewidth = 1.5,
        zorder    = 2,
    )
    ax.scatter(primary_mean, mean_y, s=34, marker="o", color=primary, edgecolor="white", linewidth=0.6, zorder=4)
    ax.scatter(recurrent_mean, mean_y, s=34, marker="s", color=recurrent, edgecolor="white", linewidth=0.6, zorder=4)

    ax.set_yticks(
        list(y) + [mean_y],
        display_table["patient_label"].tolist() + ["Pooled mean"],
        fontsize = 5.6,
    )
    ax.set_ylim(mean_y + 1.0, -1.0)
    values = display_table[["Primary", "Recurrence"]].to_numpy(dtype=float)
    margin = max(0.12, float(np.ptp(values)) * 0.08)
    ax.set_xlim(float(np.nanmin(values) - margin), float(np.nanmax(values) + margin))
    ax.set_xlabel("Patient-level MG-Inflam20 score")
    ax.grid(axis="x", color=grid_color, linewidth=0.5)
    clean_axis(ax, keep_left=False, keep_bottom=True)
    ax.tick_params(axis="y", length=0, pad=3)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc            = "lower center",
        bbox_to_anchor = (0.58, 0.035),
        frameon        = False,
        fontsize       = 5.8,
        ncol           = 2,
        handletextpad  = 0.4,
        columnspacing  = 0.8,
    )

    fig.text(0.08, 0.985, "Paired IBA1+ tissue scores", fontsize=8.5, fontweight="bold", va="top")
    pdf, png = save_panel(fig, stem)
    panel_record = record("A", stem, "Patient-level paired GeoMx MG-Inflam20 scores", source, pdf, png)
    return panel_record


def shared_genes() -> list[str]:
    """Return the prespecified shared leading-edge genes in the original Figure 1 order."""
    path = project_root / "write/51_figure1_final_reorganized/Figure1/source_data/Figure1_raw20_definition_source.csv"
    definition_table = pd.read_csv(path)
    shared_table = definition_table.loc[definition_table["shared_leading_edge"].astype(bool)]
    shared_table = shared_table.sort_values("raw20_order")
    genes = shared_table["gene"].drop_duplicates().tolist()
    if genes != ["PDK4", "SGK1", "CH25H", "SIGLEC8", "KLF6", "FOLR2", "CCL4"]:
        raise ValueError(f"Unexpected Figure 1 shared leading-edge gene order: {genes}")
    return genes
