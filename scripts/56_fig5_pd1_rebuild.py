#!/usr/bin/env python3
"""Figure 5：两套anti-PD-1暴露队列的Miller-IM program关联。"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

from recurrent_figure_style import apply_publication_style, clean_axis, new_figure


ROOT = Path(__file__).resolve().parents[1]
STEP41 = ROOT / "write/41_mg_inflammatory_sci_rebuild/Figure5"
STEP41_SOURCE = STEP41 / "source_data"
GSE154795_AUDIT = (
    ROOT
    / "write/56_figure5_pd1_rebuild"
    / "Figure5"
    / "frozen_inputs"
    / "GSE154795_myeloid_sample_audit.csv"
)
WRITE_ROOT = ROOT / "write/56_figure5_pd1_rebuild/Figure5"
SOURCE_OUT = WRITE_ROOT / "source_data"
FIG_OUT = ROOT / "figures/56_figure5_pd1_rebuild/Figure5"

MILLER_IM_GENES = [
    "PDK4", "P2RY13", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1",
    "ITM2C", "GSTM3", "CCL3", "CH25H", "P2RY12", "JUN", "SIGLEC8", "KLF6",
    "FOLR2", "CCL4", "AC253572.2", "NLRP3",
]

PRIMARY = "#3C5488"
EXPOSED = "#E64B35"
TEAL = "#00A087"
TEXT = "#222222"
NEUTRAL = "#777777"
LIGHT = "#E8E8E8"
PALE_RED = "#F8DEDA"


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
        "Creator": "Step56 recurrent GBM Figure 5",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(pdf, facecolor="white", edgecolor="none", metadata=metadata)
    fig.savefig(
        png,
        dpi=600,
        facecolor="white",
        edgecolor="none",
        metadata={"Software": "Step56 recurrent GBM Figure 5"},
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


def load_timeline() -> pd.DataFrame:
    table = pd.read_csv(STEP41_SOURCE / "figure5_treatment_timeline.csv")
    if len(table) != 7:
        raise ValueError("Figure5A必须包含四条治疗路径、七个事件")
    expected = {
        "GSE121810 neo-adjuvant (n=14)",
        "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)",
        "GSE154795 untreated recurrent (n=8)",
    }
    if set(table["lane_label"]) != expected:
        raise ValueError("Figure5A队列标签或样本量已改变")
    table["lane_label"] = table["lane_label"].replace(
        {
            "GSE154795 untreated recurrent (n=8)":
                "GSE154795 PD-1-unexposed recurrent (n=8)",
        }
    )
    return table


def panel_a() -> dict[str, str]:
    stem = "Figure5A_anti_pd1_exposure_timeline"
    table = load_timeline()
    source = SOURCE_OUT / f"{stem}_source.csv"
    table.to_csv(source, index=False)

    lane_order = [
        "GSE121810 neo-adjuvant (n=14)",
        "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)",
        "GSE154795 PD-1-unexposed recurrent (n=8)",
    ]
    y_lookup = {label: 3 - index for index, label in enumerate(lane_order)}

    fig = new_figure(178, 50)
    ax = fig.add_axes([0.31, 0.20, 0.66, 0.60])
    for label in lane_order:
        yy = y_lookup[label]
        ax.hlines(yy, 0, 2, color="#B9B9B9", linewidth=1.4, zorder=1)
        ax.scatter(1, yy, marker="D", s=32, color=TEXT, edgecolor="white", linewidth=0.5, zorder=4)

    exposed_rows = table.loc[table["event_group"].eq("Exposure")]
    for row in exposed_rows.itertuples():
        yy = y_lookup[row.lane_label]
        if row.step == 1:
            exposure_x, arrow_start, arrow_end = 0.24, 0.31, 0.94
        else:
            exposure_x, arrow_start, arrow_end = 1.76, 1.06, 1.69
        ax.scatter(exposure_x, yy, marker="o", s=42, color=EXPOSED, edgecolor="white", linewidth=0.55, zorder=5)
        ax.annotate(
            "",
            xy=(arrow_end, yy),
            xytext=(arrow_start, yy),
            arrowprops={"arrowstyle": "-|>", "color": EXPOSED, "linewidth": 1.3},
            zorder=3,
        )

    ax.text(0.24, 3.32, "Pembrolizumab", ha="center", va="bottom", fontsize=6.1)
    ax.text(1.76, 2.32, "Starts after surgery", ha="center", va="bottom", fontsize=6.1)
    ax.text(0.24, 1.32, "Anti-PD-1", ha="center", va="bottom", fontsize=6.1)
    ax.text(1.00, 2.65, "Resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, 1.65, "Resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, 0.65, "Recurrent resection", ha="center", va="top", fontsize=5.7)
    ax.text(1.00, -0.35, "Recurrent resection", ha="center", va="top", fontsize=5.7)

    ax.set_yticks([3, 2, 1, 0], lane_order, fontsize=6.3)
    ax.set_xticks([0.24, 1.0, 1.76], ["Before surgery", "Resection", "After surgery"], fontsize=6.1)
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.65, 3.55)
    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(
        0.04,
        0.955,
        "Anti-PD-1 exposure precedes tissue collection in the analyzed exposed groups",
        fontsize=9.0,
        fontweight="bold",
        va="top",
    )
    pdf, png = save_panel(fig, stem)
    return record("A", stem, "Treatment exposure timing in the two Figure5 cohorts", source, pdf, png)


def load_gsea_inputs() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    paths = {
        "GSE121810": STEP41 / "GSE121810_raw20_edger.csv",
        "GSE154795": STEP41 / "GSE154795_raw20_edger_auto_strict_idhwt.csv",
    }
    de_tables = {cohort: pd.read_csv(path) for cohort, path in paths.items()}

    stats_rows = []
    for cohort, suffix in [
        ("GSE121810", ""),
        ("GSE154795", "_auto_strict_idhwt"),
    ]:
        fgsea = pd.read_csv(STEP41 / f"{cohort}_raw20_fgsea{suffix}.csv").iloc[0]
        camera = pd.read_csv(STEP41 / f"{cohort}_raw20_camera{suffix}.csv").iloc[0]
        coverage = pd.read_csv(STEP41 / f"{cohort}_raw20_coverage{suffix}.csv").iloc[0]
        stats_rows.append(
            {
                "cohort": cohort,
                "NES": float(fgsea["NES"]),
                "fgsea_nominal_p": float(fgsea["pval"]),
                "camera_direction": str(camera["Direction"]),
                "camera_nominal_p": float(camera["PValue"]),
                "n_defined": int(coverage["n_defined"]),
                "n_present": int(coverage["n_present"]),
                "leading_edge": str(fgsea["leadingEdge"]),
            }
        )
    stats = pd.DataFrame(stats_rows)
    expected = {
        "GSE121810": (1.89995642058625, 6.3505755649683e-05, 2.81295477233785e-05),
        "GSE154795": (1.75124003154227, 0.000290544571159221, 0.015367499497791),
    }
    for row in stats.itertuples():
        values = expected[row.cohort]
        if not (
            np.isclose(row.NES, values[0], atol=1e-12)
            and np.isclose(row.fgsea_nominal_p, values[1], atol=1e-15)
            and np.isclose(row.camera_nominal_p, values[2], atol=1e-15)
            and row.n_present == 19
            and row.n_defined == 20
        ):
            raise ValueError(f"{row.cohort} Step41程序级统计已改变")
    return de_tables, stats


def compute_gsea_curve(de: pd.DataFrame, cohort: str, stats: pd.Series) -> pd.DataFrame:
    required = {"gene", "logFC", "F"}
    if not required.issubset(de.columns):
        raise ValueError(f"{cohort} edgeR表缺少列: {sorted(required - set(de.columns))}")
    table = de[["gene", "logFC", "F"]].dropna().copy()
    table["rank_stat"] = table["F"] * np.sign(table["logFC"])
    table = table.sort_values(["rank_stat", "gene"], ascending=[False, True]).reset_index(drop=True)
    table["rank"] = np.arange(1, len(table) + 1)
    table["rank_percentile"] = 100 * (table["rank"] - 1) / max(1, len(table) - 1)
    table["hit"] = table["gene"].isin(MILLER_IM_GENES)
    if int(table["hit"].sum()) != 19:
        raise ValueError(f"{cohort}必须覆盖19/20 Miller-IM基因")
    weights = table["rank_stat"].abs().to_numpy(dtype=float)
    hits = table["hit"].to_numpy(dtype=bool)
    nr = float(weights[hits].sum())
    increments = np.where(hits, weights / nr, -1.0 / int((~hits).sum()))
    table["running_enrichment"] = np.cumsum(increments)
    table["cohort"] = cohort
    table["NES"] = float(stats["NES"])
    table["fgsea_nominal_p"] = float(stats["fgsea_nominal_p"])
    table["camera_direction"] = str(stats["camera_direction"])
    table["camera_nominal_p"] = float(stats["camera_nominal_p"])
    table["n_present"] = int(stats["n_present"])
    table["n_defined"] = int(stats["n_defined"])
    return table


def sci_label(value: float, digits: int = 1) -> str:
    if value >= 0.001:
        return f"{value:.3f}"
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / (10 ** exponent)
    return rf"${coefficient:.{digits}f}\times10^{{{exponent}}}$"


def panel_b() -> dict[str, str]:
    stem = "Figure5B_dual_cohort_raw20_enrichment"
    de_tables, stats = load_gsea_inputs()
    stats_indexed = stats.set_index("cohort")
    curves = {
        cohort: compute_gsea_curve(de, cohort, stats_indexed.loc[cohort])
        for cohort, de in de_tables.items()
    }
    source_table = pd.concat(curves.values(), ignore_index=True)
    source = SOURCE_OUT / f"{stem}_source.csv.gz"
    source_table.to_csv(source, index=False, compression={"method": "gzip", "mtime": 0})

    fig = new_figure(178, 78)
    positions = {
        "GSE121810": (0.08, 0.54),
        "GSE154795": (0.56, 0.54),
    }
    comparison = {
        "GSE121810": "Pre-op exposed vs post-op only",
        "GSE154795": "PD-1 exposed vs PD-1-unexposed recurrent",
    }
    main_axes = {}
    for cohort in ["GSE121810", "GSE154795"]:
        left, _ = positions[cohort]
        curve = curves[cohort]
        stat = stats_indexed.loc[cohort]
        ax = fig.add_axes([left, 0.28, 0.40, 0.56])
        main_axes[cohort] = ax
        ax.axhline(0, color="#A5A5A5", linewidth=0.6)
        ax.fill_between(
            curve["rank_percentile"],
            0,
            curve["running_enrichment"],
            color=PALE_RED,
            alpha=0.92,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            curve["rank_percentile"],
            curve["running_enrichment"],
            color=EXPOSED,
            linewidth=1.2,
            zorder=3,
        )
        hits = curve.loc[curve["hit"]]
        rug_y0, rug_y1 = -0.105, -0.035
        ax.vlines(hits["rank_percentile"], rug_y0, rug_y1, color=TEXT, linewidth=0.65, zorder=4)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.13, 0.96)
        ax.set_xticks([0, 50, 100])
        ax.set_xlabel("Transcriptome rank percentile")
        if cohort == "GSE121810":
            ax.set_ylabel("Running enrichment score")
        else:
            ax.set_yticklabels([])
        ax.grid(axis="y", color=LIGHT, linewidth=0.45)
        clean_axis(ax, keep_left=True, keep_bottom=True)
        ax.set_title(f"{cohort}\n{comparison[cohort]}", fontsize=7.0, fontweight="bold", loc="left", pad=5)
        stat_text = (
            f"19/20 genes measured\n"
            f"NES {stat['NES']:.2f}; P {sci_label(float(stat['fgsea_nominal_p']))}\n"
            f"camera {stat['camera_direction']}; P {sci_label(float(stat['camera_nominal_p']))}"
        )
        ax.text(
            0.98,
            0.96,
            stat_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=5.7,
            color=TEXT,
            linespacing=1.25,
        )

    fig.text(
        0.04,
        0.965,
        "The Miller-IM program is enriched toward anti-PD-1 exposure in both cohorts",
        fontsize=9.0,
        fontweight="bold",
        va="top",
    )
    pdf, png = save_panel(fig, stem)
    return record("B", stem, "Two-cohort Miller-IM rank-enrichment evidence", source, pdf, png)


def load_patient_evidence() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples = pd.read_csv(STEP41_SOURCE / "figure5d_patient_score_samples.csv")
    models = pd.concat(
        [
            pd.read_csv(STEP41 / "GSE121810_raw20_score_models.csv"),
            pd.read_csv(STEP41 / "GSE154795_raw20_score_models_auto_strict_idhwt.csv"),
        ],
        ignore_index=True,
        sort=False,
    )
    if len(samples) != 50 or set(samples["cohort"]) != {"GSE121810", "GSE154795"}:
        raise ValueError("Figure5C患者分数源表必须为29+21个样本")

    audit = pd.read_csv(GSE154795_AUDIT).rename(
        columns={"ID": "sample", "condition": "audit_condition"}
    )
    audit_columns = [
        "sample", "audit_condition", "IDH_status", "MGMT", "n_myeloid_cells", "GSM",
        "sequencing_method", "gel_bead_version", "acquisition_wave",
    ]
    if not set(audit_columns).issubset(audit.columns):
        raise ValueError("GSE154795患者协变量审计表缺少模型所需列")
    samples = samples.merge(
        audit[audit_columns],
        on="sample",
        how="left",
        validate="many_to_one",
    )
    gse154795 = samples.loc[samples["cohort"].eq("GSE154795")]
    covariates = ["IDH_status", "MGMT", "sequencing_method", "gel_bead_version", "acquisition_wave"]
    if len(gse154795) != 21 or gse154795[covariates].isna().any().any():
        raise ValueError("GSE154795的21名患者必须完整匹配模型协变量")
    if not gse154795["IDH_status"].eq("WT").all():
        raise ValueError("GSE154795正式入口必须全部为IDH-wild-type")
    if not gse154795["condition"].eq(gse154795["audit_condition"]).all():
        raise ValueError("GSE154795患者分组与协变量审计表不一致")

    adjusted = models.loc[models["metric"].eq("raw20_score_myeloid_adjusted")].copy()
    raw_models = models.loc[models["metric"].eq("raw20_score")].copy()
    if len(adjusted) != 2:
        raise ValueError("Figure5C必须有两队列各一个myeloid-adjusted模型")
    if len(raw_models) != 2:
        raise ValueError("Figure5C必须有两队列各一个raw-score模型用于FDR审计")
    expected = {
        "GSE121810": (
            0.416017028595708, 0.137781145186939, 0.694252912004477,
            0.00492100468020317, 0.00984200936040635,
        ),
        "GSE154795": (
            0.427369843346401, -0.222316569830202, 1.077056256523,
            0.180119103326531, 0.180119103326531,
        ),
    }
    for row in adjusted.itertuples():
        values = expected[row.cohort]
        observed = (row.effect, row.conf_low, row.conf_high, row.p_value, row.fdr)
        if not all(np.isclose(x, y, atol=1e-12) for x, y in zip(observed, values)):
            raise ValueError(f"{row.cohort}患者模型结果已改变")
    return samples, adjusted, raw_models


def draw_box_swarm(
    ax: mpl.axes.Axes,
    values_by_group: list[np.ndarray],
    labels: list[str],
    colors: list[str],
    seed: int,
) -> None:
    rng = np.random.default_rng(seed)
    for index, (values, color) in enumerate(zip(values_by_group, colors)):
        values = np.asarray(values, dtype=float)
        jitter = rng.uniform(-0.10, 0.10, len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=18,
            color=color,
            alpha=0.78,
            edgecolor="white",
            linewidth=0.4,
            zorder=3,
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        whisker_low = max(values.min(), q1 - 1.5 * (q3 - q1))
        whisker_high = min(values.max(), q3 + 1.5 * (q3 - q1))
        ax.vlines(index, whisker_low, whisker_high, color=TEXT, linewidth=0.8, zorder=4)
        rectangle = mpl.patches.Rectangle(
            (index - 0.16, q1),
            0.32,
            q3 - q1,
            facecolor="white",
            edgecolor=TEXT,
            linewidth=0.8,
            zorder=5,
        )
        ax.add_patch(rectangle)
        ax.hlines(median, index - 0.16, index + 0.16, color=TEXT, linewidth=1.0, zorder=6)
    ax.axhline(0, color="#8A8A8A", linewidth=0.65, linestyle=(0, (3, 2)), zorder=1)
    ax.set_xticks(range(len(labels)), labels, fontsize=5.8)
    ax.set_xlim(-0.45, len(labels) - 0.55)
    ax.set_ylim(-1.65, 1.45)
    ax.grid(axis="y", color=LIGHT, linewidth=0.45)
    clean_axis(ax, keep_left=True, keep_bottom=True)
    ax.tick_params(axis="x", length=0, pad=3)


def draw_adjusted_effect(ax: mpl.axes.Axes, row: pd.Series, color: str, show_xlabel: bool) -> None:
    ax.axvline(0, color="#8A8A8A", linewidth=0.7)
    ax.hlines(0, float(row["conf_low"]), float(row["conf_high"]), color=TEXT, linewidth=1.1, zorder=2)
    significant = float(row["fdr"]) < 0.05
    ax.scatter(
        float(row["effect"]),
        0,
        marker="D",
        s=38,
        facecolor=color if significant else "white",
        edgecolor=color,
        linewidth=0.9,
        zorder=4,
    )
    ax.text(
        0.98,
        0.76,
        f"FDR {float(row['fdr']):.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.1,
        color=color if significant else TEXT,
        fontweight="bold" if significant else "normal",
    )
    ax.text(
        0.98,
        0.58,
        f"effect {float(row['effect']):+.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.6,
        color=NEUTRAL,
    )
    ax.set_xlim(-0.40, 1.25)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0])
    if show_xlabel:
        ax.set_xlabel("Adjusted exposure-associated Miller-IM score shift")
    else:
        ax.set_xticklabels([])
    ax.grid(axis="x", color=LIGHT, linewidth=0.45)
    clean_axis(ax, keep_left=False, keep_bottom=True)


def panel_c() -> dict[str, str]:
    stem = "Figure5C_patient_scores_and_adjusted_effects"
    samples, models, raw_models = load_patient_evidence()
    model_index = models.set_index("cohort")
    samples = samples.assign(
        analysis_group=samples["group"].fillna(samples["condition"]),
        model_formula=samples["cohort"].map(
            {
                "GSE121810": "raw20_score ~ myeloid_abundance + group",
                "GSE154795": (
                    "raw20_score ~ myeloid_abundance + MGMT + sequencing_method + "
                    "gel_bead_version + acquisition_wave + condition"
                ),
            }
        ),
        fdr_family="BH across raw and myeloid-adjusted patient-score models within cohort (n=2)",
    )
    source_table = samples.merge(
        models[
            [
                "cohort", "effect", "conf_low", "conf_high", "p_value", "fdr",
                "n_positive", "n_negative",
            ]
        ],
        on="cohort",
        how="left",
        validate="many_to_one",
    ).merge(
        raw_models[["cohort", "p_value", "fdr"]].rename(
            columns={"p_value": "raw_model_p_value", "fdr": "raw_model_fdr"}
        ),
        on="cohort",
        how="left",
        validate="many_to_one",
    )
    source = SOURCE_OUT / f"{stem}_source.csv"
    source_table.to_csv(source, index=False)

    fig = new_figure(178, 78)
    fig.text(
        0.04,
        0.965,
        "Adjusted Miller-IM score support is significant only in GSE121810",
        fontsize=9.0,
        fontweight="bold",
        va="top",
    )
    fig.text(0.69, 0.885, "Myeloid-proxy-adjusted model", fontsize=6.4, fontweight="bold", ha="center")

    configurations = [
        {
            "cohort": "GSE121810",
            "bottom": 0.53,
            "groups": ["Adj", "Neo"],
            "labels": ["Post-op only\n(n=15)", "Pre-op exposed\n(n=14)"],
            "title": "GSE121810 bulk tumors",
            "seed": 121810,
        },
        {
            "cohort": "GSE154795",
            "bottom": 0.13,
            "groups": ["GBM.rec", "GBM.PD1"],
            "labels": ["PD-1-unexposed recurrent\n(n=8)", "PD-1 exposed\n(n=13)"],
            "title": "GSE154795 strict IDH-wt myeloid pseudobulk",
            "seed": 154795,
        },
    ]
    for index, config in enumerate(configurations):
        cohort = config["cohort"]
        cohort_samples = samples.loc[samples["cohort"].eq(cohort)].copy()
        group_column = "group" if cohort == "GSE121810" else "condition"
        values = [
            cohort_samples.loc[cohort_samples[group_column].eq(group), "raw20_score"].to_numpy(dtype=float)
            for group in config["groups"]
        ]
        expected_n = [15, 14] if cohort == "GSE121810" else [8, 13]
        if [len(item) for item in values] != expected_n:
            raise ValueError(f"{cohort}患者分组数量改变")
        dist_ax = fig.add_axes([0.10, config["bottom"], 0.47, 0.27])
        draw_box_swarm(dist_ax, values, config["labels"], [PRIMARY, EXPOSED], config["seed"])
        dist_ax.set_title(config["title"], fontsize=6.8, fontweight="bold", loc="left", pad=3)
        if index == 0:
            dist_ax.set_ylabel("Miller-IM score")
        else:
            dist_ax.set_ylabel("Miller-IM score")

        effect_ax = fig.add_axes([0.64, config["bottom"] + 0.01, 0.31, 0.24])
        draw_adjusted_effect(effect_ax, model_index.loc[cohort], EXPOSED, show_xlabel=index == 1)

    pdf, png = save_panel(fig, stem)
    return record("C", stem, "Miller-IM scores with cohort-specific adjusted model effects", source, pdf, png)


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
        fitted = ImageOps.contain(white.convert("RGB"), (right - left, bottom - top), Image.Resampling.LANCZOS)
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def make_preview(records: list[dict[str, str]]) -> tuple[Path, Path, Path, Path]:
    by_panel = {row["panel"]: Path(row["png"]) for row in records}
    layout = [
        {"panel": "A", "left": 115, "top": 155, "right": 4205, "bottom": 1260},
        {"panel": "B", "left": 115, "top": 1375, "right": 4205, "bottom": 3155},
        {"panel": "C", "left": 115, "top": 3270, "right": 4205, "bottom": 5050},
    ]
    canvas = Image.new("RGB", (4320, 5165), "white")
    drawer = ImageDraw.Draw(canvas)
    panel_font = image_font(76, bold=True)
    for item in layout:
        panel = item["panel"]
        paste_contained(canvas, by_panel[panel], (item["left"], item["top"], item["right"], item["bottom"]))
        drawer.text((item["left"], item["top"] - 85), panel, fill="#111111", font=panel_font)

    png = FIG_OUT / "Figure5_pd1_preview.png"
    gray = FIG_OUT / "Figure5_pd1_preview_grayscale.png"
    pdf = FIG_OUT / "Figure5_pd1_preview.pdf"
    source = SOURCE_OUT / "Figure5_pd1_preview_layout.csv"
    canvas.save(png, format="PNG", dpi=(600, 600), optimize=True)
    ImageOps.grayscale(canvas).convert("RGB").save(gray, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(pdf, format="PDF", resolution=600.0)
    pd.DataFrame(layout).to_csv(source, index=False)
    return pdf, png, gray, source


def write_legend() -> Path:
    content = """# Figure 5 legend

**Figure 5 | Anti-PD-1 exposure is associated with a coordinated shift of the Miller-IM program.**

**A,** Sampling and exposure designs. In GSE121810, tumors from 14 patients were resected after preoperative pembrolizumab exposure, whereas tumors from 15 adjuvant-only patients were collected before pembrolizumab began. In GSE154795, strict IDH-wild-type patient-level myeloid pseudobulks compared 13 anti-PD-1-exposed recurrent tumors with eight recurrent tumors without neoadjuvant anti-PD-1 exposure. The diagram distinguishes exposure at resection from treatment initiated only after tissue collection. **B,** Rank-based enrichment of the prespecified Miller-IM program. Nineteen of 20 genes were measurable in each cohort; AC253572.2 was not measured. In GSE121810, the Miller-IM program was enriched toward the preoperative-exposure end of the transcriptome-wide ranking (NES = 1.900; nominal fgsea *P* = 6.35 × 10^-5), with a complementary competitive test agreeing in direction (camera, Up; nominal *P* = 2.81 × 10^-5). In the strict IDH-wild-type GSE154795 myeloid pseudobulk analysis, the program was likewise enriched toward anti-PD-1 exposure (NES = 1.751; nominal fgsea *P* = 2.91 × 10^-4; camera, Up; nominal *P* = 0.0154). Curves show the weighted running enrichment score; ticks mark Miller-IM genes in the ranked transcriptome. **C,** Left, unadjusted patient-level Miller-IM scores, displayed as individual patients with box summaries. Right, model-based exposure coefficients and 95% CIs after adjustment for the myeloid-expression proxy; the GSE154795 model additionally retained MGMT and technical covariates. GSE121810 showed an adjusted effect of +0.416 (95% CI, +0.138 to +0.694; FDR = 0.00984), whereas GSE154795 showed a positive but nonsignificant adjusted effect of +0.427 (95% CI, -0.222 to +1.077; FDR = 0.180).

These cohorts test an exposure-associated coordinated program shift, not clinical response prediction. GSE121810 is a whole-tumor bulk RNA assay and cannot localize the signal to myeloid cells; GSE154795 provides patient-level myeloid pseudobulk support but its patient-level score did not reach FDR < 0.05. The main conclusion therefore concerns coordinated rank-based program enrichment under anti-PD-1 exposure. It does not establish drug-induced causality, treatment benefit, or validity across all IDH groups; the GSE154795 all-IDH sensitivity analysis was not significant. Entry sensitivity and historical response-labeled analyses remain supplementary and are not presented as Miller-IM response validation.
"""
    path = WRITE_ROOT / "Figure5_legend.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_result(manifest: Path, legend: Path) -> None:
    content = f"""# Figure 5 anti-PD-1 exposure rebuild

- A：两队列取材时的真实anti-PD-1暴露时序。
- B：两队列Miller-IM program的rank enrichment与camera一致性。
- C：患者原始分数与myeloid-proxy-adjusted模型效应并列，保留GSE154795患者分数FDR=0.180的阴性边界。
- 不纳入volcano、四队列leading-edge UpSet或疗效预测卡片。
- Main claim：anti-PD-1 exposure-associated coordinated shift；不写induction或response prediction。
- Manifest：`{manifest}`
- Legend：`{legend}`
"""
    (WRITE_ROOT / "FINAL_RESULT.md").write_text(content, encoding="utf-8")


def main() -> None:
    apply_publication_style()
    np.random.seed(20260714)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    records = [panel_a(), panel_b(), panel_c()]
    if [row["panel"] for row in records] != list("ABC"):
        raise ValueError("Figure5面板必须唯一锁定为A-C")
    preview_pdf, preview_png, preview_gray, preview_layout = make_preview(records)
    manifest_path = WRITE_ROOT / "Figure5_panel_manifest.csv"
    manifest = pd.DataFrame(records)
    manifest["preview_pdf"] = str(preview_pdf)
    manifest["preview_png"] = str(preview_png)
    manifest["preview_grayscale"] = str(preview_gray)
    manifest["preview_layout"] = str(preview_layout)
    manifest.to_csv(manifest_path, index=False)
    legend = write_legend()
    write_result(manifest_path, legend)
    print("STEP56_FIGURE5_COMPLETE panels=3 preview=1 mapping=A-C")


if __name__ == "__main__":
    main()
