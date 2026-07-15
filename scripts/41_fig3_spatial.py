#!/usr/bin/env python3
from __future__ import annotations

import gzip
import io
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from recurrent_figure_style import COLORS, apply_publication_style, clean_axis, new_figure, save_figure


ROOT = Path(__file__).resolve().parents[1]
OLD_FIG3_SOURCE = ROOT / "write/39_independent_miller_mg_inflammatory_figures/Figure3/source_csv"
OLD_S6_SOURCE = ROOT / "write/39_independent_miller_mg_inflammatory_figures/SupplementaryFigure6/source_csv"
DATA_ROOT = Path(os.environ.get("MILLER_IM_DATA_ROOT", ROOT / "data")).expanduser().resolve()
SPATIAL_IMAGE_ROOT = DATA_ROOT / "write/30_dataset_rescue_search/source_metadata/GSE276841_spatial"

WRITE_FIG3 = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure3"
WRITE_S5 = ROOT / "write/41_mg_inflammatory_sci_rebuild/SupplementaryFigure5"
FIG3_SOURCE_OUT = WRITE_FIG3 / "source_data"
S5_SOURCE_OUT = WRITE_S5 / "source_data"
FIG3_OUT = ROOT / "figures/41_mg_inflammatory_sci_rebuild/Figure3/panel_library"
S5_OUT = ROOT / "figures/41_mg_inflammatory_sci_rebuild/SupplementaryFigure5/panel_library"

SPATIAL_PREFIX = {
    "GBM030": "GSM8506649_03422",
    "GBM049": "GSM8506650_25526",
}
METRIC_META = {
    "raw20": {"title": "Mg-inflammatory raw20", "cmap": "viridis"},
    "mdsc_combined": {"title": "MDSC score", "cmap": "viridis"},
    "mes": {"title": "MES score", "cmap": "viridis"},
}
PARTIAL_COLORS = {
    "MDSC": "#00A087",
    "MES": "#7E6148",
    "E-MDSC": "#4DBBD5",
    "M-MDSC": "#3C5488",
}


def ensure_dirs() -> None:
    for path in [FIG3_SOURCE_OUT, S5_SOURCE_OUT, FIG3_OUT, S5_OUT]:
        path.mkdir(parents=True, exist_ok=True)


def mean_ci(values: pd.Series | np.ndarray) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        return np.nan, np.nan, np.nan
    mean = float(array.mean())
    if len(array) == 1:
        return mean, mean, mean
    sem = float(array.std(ddof=1) / np.sqrt(len(array)))
    delta = 1.96 * sem
    return mean, mean - delta, mean + delta


def fisher_ci(rho: float, n_obs: int) -> tuple[float, float]:
    if not np.isfinite(rho) or n_obs <= 3 or abs(rho) >= 1:
        return np.nan, np.nan
    z_value = np.arctanh(rho)
    sem = 1 / np.sqrt(n_obs - 3)
    low, high = z_value - 1.96 * sem, z_value + 1.96 * sem
    return float(np.tanh(low)), float(np.tanh(high))


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def load_grayscale_image(sample: str) -> tuple[np.ndarray, float, float]:
    prefix = SPATIAL_PREFIX[sample]
    with gzip.open(SPATIAL_IMAGE_ROOT / f"{prefix}_tissue_hires_image.png.gz", "rb") as handle:
        image = Image.open(io.BytesIO(handle.read())).convert("L")
    with gzip.open(SPATIAL_IMAGE_ROOT / f"{prefix}_scalefactors_json.json.gz", "rt") as handle:
        scale_meta = json.load(handle)
    gray = np.asarray(image, dtype=np.float32)
    scale = float(scale_meta["tissue_hires_scalef"])
    spot_diameter = float(scale_meta["spot_diameter_fullres"]) * scale
    return gray, scale, spot_diameter


def build_fig3a_table() -> pd.DataFrame:
    table = load_table(OLD_FIG3_SOURCE / "geomx_raw20_paired_deltas.csv")
    summary = load_table(OLD_FIG3_SOURCE / "geomx_raw20_entry_level_summary.csv").set_index("entry")
    table = table.loc[table["entry"].eq("strict_pass_idhwt")].copy()
    mean_delta, ci_low, ci_high = mean_ci(table["delta"])
    table["mean_delta"] = mean_delta
    table["ci_low"] = ci_low
    table["ci_high"] = ci_high
    table["n_pairs"] = len(table)
    table["n_up"] = int((table["delta"] > 0).sum())
    table["fdr"] = float(summary.loc["strict_pass_idhwt", "fdr"])
    out_path = FIG3_SOURCE_OUT / "fig3a_geomx22_paired_change.csv"
    table.to_csv(out_path, index=False)
    return table


def build_fig3b_table() -> pd.DataFrame:
    table = load_table(OLD_FIG3_SOURCE / "fig3c_geomx_gene_forest.csv")
    table = table.sort_values(["mean_log2_delta", "gene"], ascending=[False, True]).reset_index(drop=True)
    table.to_csv(FIG3_SOURCE_OUT / "fig3b_geomx18_gene_forest.csv", index=False)
    return table


def build_spatial_tables() -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, float]]]:
    all_spots = load_table(OLD_FIG3_SOURCE / "gse276841_raw20_spot_scores.csv.gz")
    spatial_tables: dict[str, pd.DataFrame] = {}
    spatial_meta: dict[str, dict[str, float]] = {}
    for sample in ["GBM030", "GBM049"]:
        gray, scale, spot_diameter = load_grayscale_image(sample)
        current = all_spots.loc[all_spots["sample"].eq(sample)].copy()
        current["x_hires"] = current["pxl_col_in_fullres"] * scale
        current["y_hires"] = current["pxl_row_in_fullres"] * scale
        current["image_width_px"] = gray.shape[1]
        current["image_height_px"] = gray.shape[0]
        current["tissue_hires_scalef"] = scale
        current["spot_diameter_hires_px"] = spot_diameter
        current.to_csv(FIG3_SOURCE_OUT / f"{'fig3c' if sample == 'GBM030' else 'fig3d'}_{sample.lower()}_spatial_raw20_mdsc_mes.csv", index=False)
        spatial_tables[sample] = current
        spatial_meta[sample] = {
            "width": float(gray.shape[1]),
            "height": float(gray.shape[0]),
            "scale": scale,
            "spot_diameter": spot_diameter,
        }
    return spatial_tables, spatial_meta


def build_fig3e_table() -> pd.DataFrame:
    table = load_table(OLD_FIG3_SOURCE / "fig3f_partial_rho.csv")
    table = table.loc[table["target"].isin(["MDSC", "MES"])].copy()
    ci_bounds = table.apply(lambda row: fisher_ci(float(row["partial_rho"]), int(row["n_spots"])), axis=1, result_type="expand")
    table["ci_low"] = ci_bounds[0].to_numpy()
    table["ci_high"] = ci_bounds[1].to_numpy()
    table["panel_label"] = table["sample_label"] + " · " + table["target"]
    table.to_csv(FIG3_SOURCE_OUT / "fig3e_partial_rho_main.csv", index=False)
    return table


def _entry_ci(deltas: pd.Series) -> tuple[float, float, float]:
    return mean_ci(deltas.to_numpy(dtype=float))


def _interaction_ci(control: pd.Series, nivolumab: pd.Series) -> tuple[float, float, float]:
    control_arr = control.to_numpy(dtype=float)
    nivo_arr = nivolumab.to_numpy(dtype=float)
    effect = float(nivo_arr.mean() - control_arr.mean())
    if len(control_arr) <= 1 or len(nivo_arr) <= 1:
        return effect, np.nan, np.nan
    se = np.sqrt(control_arr.var(ddof=1) / len(control_arr) + nivo_arr.var(ddof=1) / len(nivo_arr))
    delta = 1.96 * se
    return effect, float(effect - delta), float(effect + delta)


def build_s5a_table() -> pd.DataFrame:
    paired = load_table(OLD_FIG3_SOURCE / "geomx_raw20_paired_deltas.csv")
    entry_summary = load_table(OLD_FIG3_SOURCE / "geomx_raw20_entry_level_summary.csv").set_index("entry")
    subgroup_summary = load_table(OLD_S6_SOURCE / "geomx_raw20_subgroup_summary.csv")
    subgroup_summary = subgroup_summary.set_index(["entry", "subgroup"])

    rows: list[dict[str, object]] = []

    def add_row(section: str, order: int, label: str, deltas: pd.Series, fdr: float, n_up: int | None = None) -> None:
        effect, ci_low, ci_high = _entry_ci(deltas)
        rows.append(
            {
                "section": section,
                "section_order": order,
                "label": label,
                "effect": effect,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "n_pairs": int(deltas.notna().sum()),
                "n_up": int((deltas > 0).sum()) if n_up is None else int(n_up),
                "fdr": float(fdr),
            }
        )

    main = paired.loc[paired["entry"].eq("strict_pass_idhwt"), "delta"]
    add_row(
        "Main analysis",
        1,
        "IBA1+ strict + IDH-wt",
        main,
        float(entry_summary.loc["strict_pass_idhwt", "fdr"]),
        int(entry_summary.loc["strict_pass_idhwt", "n_positive"]),
    )

    for order, entry, label in [
        (2, "all_aoi_idhwt", "IBA1+ all AOI + IDH-wt"),
        (2, "strict_pass_all_idh", "IBA1+ strict + all IDH"),
    ]:
        subset = paired.loc[paired["entry"].eq(entry), "delta"]
        add_row(section="Entry sensitivity", order=order, label=label, deltas=subset, fdr=float(entry_summary.loc[entry, "fdr"]))

    for subgroup, label in [
        ("IDH_WT", "All-IDH entry: IDH-wt"),
        ("IDH_mutated", "All-IDH entry: IDH-mutant"),
    ]:
        subset = paired.loc[
            paired["entry"].eq("strict_pass_all_idh")
            & paired["IDH_status"].eq("IDH_WT" if subgroup == "IDH_WT" else "IDH_mutated"),
            "delta",
        ]
        add_row(
            section="IDH sensitivity",
            order=3,
            label=label,
            deltas=subset,
            fdr=float(subgroup_summary.loc[("strict_pass_all_idh", subgroup), "fdr"]),
        )

    for subgroup, label in [("Control", "Control"), ("Nivolumab", "Nivolumab")]:
        subset = paired.loc[
            paired["entry"].eq("strict_pass_idhwt") & paired["trial_setting"].eq(subgroup),
            "delta",
        ]
        add_row(
            section="Treatment strata",
            order=4,
            label=label,
            deltas=subset,
            fdr=float(subgroup_summary.loc[("strict_pass_idhwt", subgroup), "fdr"]),
        )

    control = paired.loc[
        paired["entry"].eq("strict_pass_idhwt") & paired["trial_setting"].eq("Control"),
        "delta",
    ]
    nivolumab = paired.loc[
        paired["entry"].eq("strict_pass_idhwt") & paired["trial_setting"].eq("Nivolumab"),
        "delta",
    ]
    effect, ci_low, ci_high = _interaction_ci(control, nivolumab)
    rows.append(
        {
            "section": "Interaction",
            "section_order": 5,
            "label": "Nivolumab minus control",
            "effect": effect,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_pairs": int(len(control) + len(nivolumab)),
            "n_up": np.nan,
            "fdr": float(subgroup_summary.loc[("strict_pass_idhwt", "Interaction_Nivolumab_minus_Control"), "fdr"]),
        }
    )

    forest = pd.DataFrame(rows)
    forest["row_order"] = range(len(forest))
    forest.to_csv(S5_SOURCE_OUT / "s5a_geomx_sensitivity_forest.csv", index=False)
    return forest


def build_s5b_table() -> pd.DataFrame:
    table = load_table(OLD_S6_SOURCE / "gse276841_raw20_slice_sensitivity.csv")
    ci_bounds = table.apply(lambda row: fisher_ci(float(row["partial_rho"]), int(row["n_spots"])), axis=1, result_type="expand")
    table["ci_low"] = ci_bounds[0].to_numpy()
    table["ci_high"] = ci_bounds[1].to_numpy()
    table["panel_label"] = table["sample"] + " · " + table["target"]
    table.to_csv(S5_SOURCE_OUT / "s5b_slice_sensitivity_forest.csv", index=False)
    return table


def plot_fig3a(table: pd.DataFrame) -> None:
    fig = new_figure(88, 88)
    ax = fig.add_axes([0.16, 0.15, 0.78, 0.70])
    x_positions = np.array([0.0, 1.0])
    primary = table["Primary"].to_numpy(dtype=float)
    recurrent = table["Recurrence"].to_numpy(dtype=float)

    for left, right in zip(primary, recurrent):
        ax.plot(x_positions, [left, right], color="#C4C4C4", linewidth=0.8, zorder=1)
    ax.scatter(np.full(len(primary), x_positions[0]), primary, s=16, color=COLORS["primary"], linewidths=0, zorder=3)
    ax.scatter(np.full(len(recurrent), x_positions[1]), recurrent, s=16, color=COLORS["recurrent"], linewidths=0, zorder=3)

    mean_delta = float(table["mean_delta"].iloc[0])
    ci_low = float(table["ci_low"].iloc[0])
    ci_high = float(table["ci_high"].iloc[0])
    n_pairs = int(table["n_pairs"].iloc[0])
    n_up = int(table["n_up"].iloc[0])
    fdr = float(table["fdr"].iloc[0])

    ax.set_xticks(x_positions, ["Primary", "Recurrent"])
    ax.set_ylabel("Raw20 score")
    ax.set_xlim(-0.3, 1.3)
    clean_axis(ax)
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.6)
    fig.text(0.16, 0.96, "GeoMx IBA1+ paired change", ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(
        0.16,
        0.90,
        f"22 IDH-wt pairs | 16/22 recurrent-higher | mean delta {mean_delta:+.3f}",
        ha="left",
        va="top",
        fontsize=7.0,
        color=COLORS["neutral"],
    )
    fig.text(
        0.16,
        0.05,
        f"95% CI {ci_low:+.3f} to {ci_high:+.3f} | FDR {fdr:.3f}",
        ha="left",
        va="bottom",
        fontsize=7.0,
        color=COLORS["text"],
    )
    save_figure(fig, FIG3_OUT, "fig3a_geomx22_paired_change")


def plot_fig3b(table: pd.DataFrame) -> None:
    fig = new_figure(120, 112)
    ax = fig.add_axes([0.23, 0.10, 0.71, 0.82])
    y_values = np.arange(len(table))
    colors = np.where(table["mean_log2_delta"].to_numpy(dtype=float) >= 0, COLORS["recurrent"], COLORS["primary"])
    ax.hlines(y_values, table["ci_low"], table["ci_high"], color=colors, linewidth=1.05)
    ax.scatter(table["mean_log2_delta"], y_values, s=20, color=colors, zorder=3)
    ax.axvline(0, color=COLORS["neutral"], linewidth=0.7)
    ax.set_yticks(y_values, table["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean paired log2 delta")
    clean_axis(ax)
    ax.grid(axis="x", color=COLORS["grid"], linewidth=0.6)
    x_span = float(table["ci_high"].max() - table["ci_low"].min())
    x_text = float(table["ci_high"].max()) + x_span * 0.08
    ax.set_xlim(float(table["ci_low"].min()) - x_span * 0.08, float(table["ci_high"].max()) + x_span * 0.46)
    for yy, n_pairs, fdr in zip(y_values, table["n_pairs"], table["fdr"]):
        ax.text(x_text, yy, f"n={int(n_pairs)}  FDR={fdr:.3f}", va="center", ha="left", fontsize=6.4, color=COLORS["text"])
    fig.text(0.23, 0.97, "Measurable raw20 genes in GeoMx", ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(0.23, 0.93, "18 genes detectable in strict IBA1+ IDH-wt pairs", ha="left", va="top", fontsize=7.0, color=COLORS["neutral"])
    save_figure(fig, FIG3_OUT, "fig3b_geomx18_gene_forest")


def _format_grayscale(gray: np.ndarray) -> np.ndarray:
    low = float(np.nanpercentile(gray, 2))
    high = float(np.nanpercentile(gray, 98))
    scaled = np.clip((gray - low) / max(high - low, 1e-6), 0, 1)
    return 0.12 + 0.78 * scaled


def plot_spatial_panel(sample: str, table: pd.DataFrame, stem: str) -> None:
    gray, scale, _ = load_grayscale_image(sample)
    gray_img = _format_grayscale(gray)

    all_spots = pd.concat(
        [load_table(FIG3_SOURCE_OUT / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv"), load_table(FIG3_SOURCE_OUT / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv")],
        ignore_index=True,
    )
    limits = {
        metric: np.nanquantile(all_spots[metric].to_numpy(dtype=float), [0.02, 0.98])
        for metric in METRIC_META
    }

    fig = new_figure(183, 70)
    lefts = [0.04, 0.36, 0.68]
    axes = [fig.add_axes([left, 0.19, 0.23, 0.63]) for left in lefts]
    cbar_axes = [fig.add_axes([left, 0.08, 0.23, 0.03]) for left in lefts]

    for ax, cax, (metric, meta) in zip(axes, cbar_axes, METRIC_META.items()):
        ax.imshow(gray_img, cmap="gray", origin="upper", interpolation="nearest")
        mask = gray_img < 0.96
        ax.contour(mask.astype(float), levels=[0.5], colors=["#8B8B8B"], linewidths=0.55, origin="upper")
        ax.scatter(table["x_hires"], table["y_hires"], s=11, color="#D1D1D1", alpha=0.85, linewidths=0, rasterized=True, zorder=2)
        points = ax.scatter(
            table["x_hires"],
            table["y_hires"],
            s=7.5,
            c=table[metric],
            cmap=meta["cmap"],
            vmin=float(limits[metric][0]),
            vmax=float(limits[metric][1]),
            linewidths=0,
            rasterized=True,
            zorder=3,
        )
        ax.set_title(meta["title"], pad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlim(0, gray.shape[1])
        ax.set_ylim(gray.shape[0], 0)
        colorbar = fig.colorbar(points, cax=cax, orientation="horizontal")
        colorbar.outline.set_visible(False)
        colorbar.ax.tick_params(labelsize=6.2, length=1.8, pad=1)

    fig.text(0.04, 0.97, sample, ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(
        0.96,
        0.97,
        f"{len(table):,} in-tissue Visium spots | shared scale across slices",
        ha="right",
        va="top",
        fontsize=7.0,
        color=COLORS["neutral"],
    )
    save_figure(fig, FIG3_OUT, stem)


def plot_fig3e(table: pd.DataFrame) -> None:
    ordered = table.set_index("panel_label").loc[
        ["GBM030 · MDSC", "GBM049 · MDSC", "GBM030 · MES", "GBM049 · MES"]
    ].reset_index()
    fig = new_figure(88, 82)
    ax = fig.add_axes([0.28, 0.12, 0.66, 0.68])
    y_values = np.arange(len(ordered))
    for yy, row in zip(y_values, ordered.itertuples(index=False)):
        color = PARTIAL_COLORS[row.target]
        ax.hlines(yy, row.ci_low, row.ci_high, color=color, linewidth=1.2)
        ax.scatter(row.partial_rho, yy, s=28, color=color, zorder=3)
        ax.text(row.partial_rho + 0.018, yy, f"n={int(row.n_spots)}", va="center", ha="left", fontsize=6.4, color=COLORS["text"])
    ax.axvline(0, color=COLORS["neutral"], linewidth=0.7)
    ax.set_yticks(y_values, ordered["panel_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Partial Spearman rho")
    clean_axis(ax)
    ax.grid(axis="x", color=COLORS["grid"], linewidth=0.6)
    ax.set_xlim(0.0, max(0.64, float(ordered["ci_high"].max()) + 0.13))
    fig.text(0.28, 0.97, "Slice-level spatial effect sizes", ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(
        0.28,
        0.92,
        "Descriptive only; adjusted for total myeloid score\nNo smoothing; no spatial-autocorrelation correction",
        ha="left",
        va="top",
        fontsize=6.6,
        linespacing=1.15,
        color=COLORS["neutral"],
    )
    save_figure(fig, FIG3_OUT, "fig3e_partial_rho_effects")


def _plot_grouped_forest(
    table: pd.DataFrame,
    stem: str,
    out_dir: Path,
    title: str,
    subtitle: str,
    x_label: str,
    height_mm: float,
) -> None:
    sections = []
    y_cursor = 0.0
    row_positions: list[float] = []
    for section in table["section"].drop_duplicates():
        sections.append((section, y_cursor - 0.55))
        current = table.loc[table["section"].eq(section)]
        for _ in range(len(current)):
            row_positions.append(y_cursor)
            y_cursor += 1.0
        y_cursor += 0.55
    table = table.copy()
    table["y"] = row_positions

    fig = new_figure(120, height_mm)
    label_ax = fig.add_axes([0.04, 0.13, 0.40, 0.76])
    forest_ax = fig.add_axes([0.47, 0.13, 0.49, 0.76])
    label_ax.axis("off")

    x_min = float(np.nanmin(table["ci_low"].to_numpy(dtype=float)))
    x_max = float(np.nanmax(table["ci_high"].to_numpy(dtype=float)))
    if not np.isfinite(x_min):
        x_min = float(np.nanmin(table["effect"].to_numpy(dtype=float)))
    if not np.isfinite(x_max):
        x_max = float(np.nanmax(table["effect"].to_numpy(dtype=float)))
    span = max(x_max - x_min, 0.16)
    forest_ax.set_xlim(x_min - span * 0.16, x_max + span * 0.58)

    for section, section_y in sections:
        label_ax.text(0.0, section_y, section, ha="left", va="center", fontsize=6.5, color=COLORS["neutral"], weight="bold")

    for row in table.itertuples(index=False):
        color = COLORS["recurrent"] if row.effect >= 0 else COLORS["primary"]
        if np.isfinite(row.ci_low) and np.isfinite(row.ci_high):
            forest_ax.hlines(row.y, row.ci_low, row.ci_high, color=color, linewidth=1.1)
        forest_ax.scatter(row.effect, row.y, s=24, color=color, zorder=3)
        label_ax.text(0.0, row.y, row.label, ha="left", va="center", fontsize=7.0, color=COLORS["text"])
        side_bits = [f"n={int(row.n_pairs)}"]
        if pd.notna(row.n_up):
            side_bits.append(f"{int(row.n_up)} up")
        side_bits.append(f"FDR={row.fdr:.3f}")
        forest_ax.text(row.effect + span * 0.08, row.y, "  ".join(side_bits), ha="left", va="center", fontsize=6.2, color=COLORS["text"])

    forest_ax.axvline(0, color=COLORS["neutral"], linewidth=0.7)
    forest_ax.set_yticks([])
    forest_ax.set_xlabel(x_label)
    clean_axis(forest_ax, keep_left=False, keep_bottom=True)
    forest_ax.grid(axis="x", color=COLORS["grid"], linewidth=0.6)
    forest_ax.invert_yaxis()
    label_ax.set_ylim(forest_ax.get_ylim())
    fig.text(0.04, 0.97, title, ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(0.04, 0.93, subtitle, ha="left", va="top", fontsize=7.0, color=COLORS["neutral"])
    save_figure(fig, out_dir, stem)


def plot_s5a(table: pd.DataFrame) -> None:
    _plot_grouped_forest(
        table=table,
        stem="s5a_geomx_sensitivity_forest",
        out_dir=S5_OUT,
        title="GeoMx sensitivity analyses",
        subtitle="Entry, IDH, treatment and interaction on the same effect scale",
        x_label="Mean paired raw20 delta",
        height_mm=140,
    )


def plot_s5b(table: pd.DataFrame) -> None:
    order = []
    for section, targets in [("MDSC family", ["MDSC", "E-MDSC", "M-MDSC"]), ("MES", ["MES"])]:
        current = table.loc[table["target"].isin(targets)].copy()
        current["section"] = section
        current["label"] = current["sample"] + " · " + current["target"]
        order.append(current)
    forest = pd.concat(order, ignore_index=True)
    forest["effect"] = forest["partial_rho"]
    forest = forest[["section", "label", "effect", "ci_low", "ci_high", "n_spots", "fdr", "target"]].rename(columns={"n_spots": "n_pairs"})
    forest["n_up"] = np.nan
    forest["section_order"] = forest["section"].map({"MDSC family": 1, "MES": 2})
    forest["row_order"] = range(len(forest))
    forest.to_csv(S5_SOURCE_OUT / "s5b_slice_sensitivity_forest.csv", index=False)

    sections = []
    y_cursor = 0.0
    row_positions: list[float] = []
    for section in forest["section"].drop_duplicates():
        sections.append((section, y_cursor - 0.35))
        current = forest.loc[forest["section"].eq(section)]
        for _ in range(len(current)):
            row_positions.append(y_cursor)
            y_cursor += 1.0
        y_cursor += 0.55
    forest = forest.copy()
    forest["y"] = row_positions

    fig = new_figure(120, 118)
    label_ax = fig.add_axes([0.04, 0.10, 0.36, 0.70])
    forest_ax = fig.add_axes([0.44, 0.10, 0.52, 0.70])
    label_ax.axis("off")
    for section, section_y in sections:
        label_ax.text(0.0, section_y, section, ha="left", va="center", fontsize=6.5, color=COLORS["neutral"], weight="bold")

    x_min = float(np.nanmin(forest["ci_low"]))
    x_max = float(np.nanmax(forest["ci_high"]))
    span = max(x_max - x_min, 0.18)
    forest_ax.set_xlim(x_min - span * 0.16, x_max + span * 0.50)
    for row in forest.itertuples(index=False):
        color = PARTIAL_COLORS[row.target]
        forest_ax.hlines(row.y, row.ci_low, row.ci_high, color=color, linewidth=1.1)
        forest_ax.scatter(row.effect, row.y, s=24, color=color, zorder=3)
        label_ax.text(0.0, row.y, row.label, ha="left", va="center", fontsize=7.0, color=COLORS["text"])
        forest_ax.text(row.effect + span * 0.08, row.y, f"n={int(row.n_pairs)}", ha="left", va="center", fontsize=6.3, color=COLORS["text"])
    forest_ax.axvline(0, color=COLORS["neutral"], linewidth=0.7)
    forest_ax.set_yticks([])
    forest_ax.set_xlabel("Partial Spearman rho")
    clean_axis(forest_ax, keep_left=False, keep_bottom=True)
    forest_ax.grid(axis="x", color=COLORS["grid"], linewidth=0.6)
    forest_ax.invert_yaxis()
    label_ax.set_ylim(forest_ax.get_ylim())
    fig.text(0.04, 0.97, "Spatial slice sensitivity", ha="left", va="top", fontsize=8.5, weight="bold")
    fig.text(
        0.04,
        0.91,
        "Descriptive partial rho only; each row is one slice-by-target partial correlation; no smoothing",
        ha="left",
        va="top",
        fontsize=7.0,
        color=COLORS["neutral"],
    )
    save_figure(fig, S5_OUT, "s5b_slice_sensitivity_forest")


def write_manifests() -> None:
    fig3_manifest = pd.DataFrame(
        [
            ("Figure3A", "fig3a_geomx22_paired_change", "GeoMx paired raw20 change", FIG3_SOURCE_OUT / "fig3a_geomx22_paired_change.csv"),
            ("Figure3B", "fig3b_geomx18_gene_forest", "Measurable gene forest", FIG3_SOURCE_OUT / "fig3b_geomx18_gene_forest.csv"),
            ("Figure3C", "fig3c_gbm030_spatial_raw20_mdsc_mes", "GBM030 spatial maps", FIG3_SOURCE_OUT / "fig3c_gbm030_spatial_raw20_mdsc_mes.csv"),
            ("Figure3D", "fig3d_gbm049_spatial_raw20_mdsc_mes", "GBM049 spatial maps", FIG3_SOURCE_OUT / "fig3d_gbm049_spatial_raw20_mdsc_mes.csv"),
            ("Figure3E", "fig3e_partial_rho_effects", "Main slice partial-rho effect plot", FIG3_SOURCE_OUT / "fig3e_partial_rho_main.csv"),
        ],
        columns=["panel", "stem", "message", "source_csv"],
    )
    fig3_manifest["png"] = fig3_manifest["stem"].map(lambda stem: str(FIG3_OUT / f"{stem}.png"))
    fig3_manifest["pdf"] = fig3_manifest["stem"].map(lambda stem: str(FIG3_OUT / f"{stem}.pdf"))
    fig3_manifest["source_csv"] = fig3_manifest["source_csv"].map(str)
    fig3_manifest.to_csv(WRITE_FIG3 / "Figure3_panel_manifest.csv", index=False)

    s5_manifest = pd.DataFrame(
        [
            ("SupplementaryFigure5A", "s5a_geomx_sensitivity_forest", "GeoMx sensitivity forest", S5_SOURCE_OUT / "s5a_geomx_sensitivity_forest.csv"),
            ("SupplementaryFigure5B", "s5b_slice_sensitivity_forest", "Spatial slice sensitivity forest", S5_SOURCE_OUT / "s5b_slice_sensitivity_forest.csv"),
        ],
        columns=["panel", "stem", "message", "source_csv"],
    )
    s5_manifest["png"] = s5_manifest["stem"].map(lambda stem: str(S5_OUT / f"{stem}.png"))
    s5_manifest["pdf"] = s5_manifest["stem"].map(lambda stem: str(S5_OUT / f"{stem}.pdf"))
    s5_manifest["source_csv"] = s5_manifest["source_csv"].map(str)
    s5_manifest.to_csv(WRITE_S5 / "SupplementaryFigure5_panel_manifest.csv", index=False)


def main() -> None:
    np.random.seed(20260713)
    apply_publication_style()
    ensure_dirs()

    fig3a = build_fig3a_table()
    fig3b = build_fig3b_table()
    spatial_tables, _ = build_spatial_tables()
    fig3e = build_fig3e_table()
    s5a = build_s5a_table()
    s5b = build_s5b_table()

    plot_fig3a(fig3a)
    plot_fig3b(fig3b)
    plot_spatial_panel("GBM030", spatial_tables["GBM030"], "fig3c_gbm030_spatial_raw20_mdsc_mes")
    plot_spatial_panel("GBM049", spatial_tables["GBM049"], "fig3d_gbm049_spatial_raw20_mdsc_mes")
    plot_fig3e(fig3e)
    plot_s5a(s5a)
    plot_s5b(s5b)
    write_manifests()
    print("STEP41_FIGURE3_S5_COMPLETE figure3_panels=5 s5_panels=2")


if __name__ == "__main__":
    main()
