#!/usr/bin/env python3
# 目的：绘制独立配对队列设计图；保留供下游调用的图一与图二绘图函数。
# 输入：下列读取位置所列write源表；dataset、gene、logFC、F、NES等公开字段保持原名。
# 输出：本步骤write源数据、清单及figures PDF/PNG；原文件名与图面参数保持不变。
# 运行：在项目环境执行 python scripts/05c_复发证据布局.py。
# 保留Python：既有AnnData/NumPy/Pandas/Matplotlib复现路径不在本轮改为R。
# 阈值、基因集、种子与检验顺序均沿用；600 dpi PNG由对应PDF母版渲染。
# 数据集名、第三方API参数和输出字段是公开接口，不按本轮变量规范改写。

from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr
from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
mm = figure_style.mm
sha256 = figure_style.sha256
save_panel = figure_style.save_panel
rank_density = figure_style.rank_density


def configure_style() -> None:
    """无输入；设置本步骤既有Matplotlib字体、颜色及线宽；无返回值。"""
    text = "#252525"

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.8,
            "axes.linewidth": 0.65,
            "axes.edgecolor": text,
            "axes.labelcolor": text,
            "text.color": text,
            "xtick.color": text,
            "ytick.color": text,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


root = Path(__file__).resolve().parents[1]

# 在使用处确定公开输出路径，沿用原目录和文件合同。
write_root = root / "write/42_visual_story_rebuild"


def write_source(source_table: pd.DataFrame, figure_name: str, stem: str) -> Path:
    """输入源数据表、图编号与文件主名；保存原字段CSV并核对非空，返回输出路径。"""

    # 在使用处确定公开输出路径，沿用原目录和文件合同。
    out = write_root / figure_name / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}_source.csv"
    csv_path = path
    source_table.to_csv(csv_path, index=False)

    # 按既定条件选择处理分支；原输入断言和异常保持不变。
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
    return path


def panel_record(
    panel: str, stem: str, message: str, source: Path, pdf: Path, png: Path
) -> dict[str, object]:
    """输入面板标识、说明和源表/PDF/PNG路径；返回含既有字段及文件摘要的清单记录。"""
    result = {
        "panel": panel,
        "stem": stem,
        "message": message,
        "source_csv": str(source),
        "pdf": str(pdf),
        "png": str(png),
        "source_sha256": sha256(source),
        "pdf_sha256": sha256(pdf),
        "png_sha256": sha256(png),
        "png_dpi": 600,
    }
    return result


# 此输出根由51动态导入后覆写，必须保留模块级接口。
fig_root = root / "figures/42_visual_story_rebuild"


def figure1_design() -> dict[str, object]:
    """读取Step38阈值20队列审计表；保存设计源表和PDF/PNG，返回本模块原十字段清单记录。"""

    # 读取正式配对审计表；阈值20与原队列计数断言保持不变。
    step38 = root / "write/38_independent_cohort_mg_inflammatory_recalculation"

    audit_table = pd.read_csv(step38 / "independent_input_pair_audit.csv")
    datasets = ["GSE174554", "GSE274546"]

    source_table = audit_table.loc[
        audit_table["dataset"].isin(datasets) & audit_table["threshold"].eq(20)
    ].copy()
    if source_table.groupby("dataset").size().to_dict() != {"GSE174554": 1, "GSE274546": 1}:
        raise ValueError("Step38 must contain one formal threshold-20 audit row per cohort.")
    source_table["input_libraries_or_matrices"] = source_table["dataset"].map(
        {"GSE174554": 91, "GSE274546": 111}
    )

    # 按公开列顺序组织设计源表，维持队列顺序。
    source_table["clean_myeloid"] = source_table["n_cells_primary"].astype(int) + source_table[
        "n_cells_recurrent"
    ].astype(int)
    source_table["formal_pairs_threshold20"] = source_table["n_pairs"].astype(int)
    observed = (
        source_table.set_index("dataset")[["formal_pairs_threshold20", "clean_myeloid"]]
        .astype(int)
        .to_dict("index")
    )
    expected = {
        "GSE174554": {"formal_pairs_threshold20": 17, "clean_myeloid": 12489},
        "GSE274546": {"formal_pairs_threshold20": 45, "clean_myeloid": 55568},
    }
    if observed != expected:
        raise ValueError(f"Unexpected IDH-restricted Figure1A analysis set: {observed}")

    # 按公开列顺序组织设计源表，维持队列顺序。
    source_table = source_table.sort_values("dataset").reset_index(drop=True)
    source_table["modality"] = "snRNA-seq"
    source_table["raw_input_unit"] = ["10x matrices", "10x libraries"]
    source_table["minimum_cells_per_condition"] = 20
    source_table["processing_scope"] = "Each cohort reconstructed, annotated and tested independently"

    # 按公开列顺序组织设计源表，维持队列顺序。
    source_table["statistical_testing"] = (
        "Patient pseudobulk, paired edgeR and Miller-IM GSEA performed separately within each cohort"
    )
    source_table["replication_assessment"] = "Cohort results compared only after independent testing"
    source_table["workflow"] = (
        "all-cell reconstruction -> pan-myeloid re-clustering -> marker an"
        "d cluster-DEG review -> patient raw-count pseudobulk -> paired ed"
        "geR -> prespecified Miller-IM program GSEA"
    )
    safe_source_columns = [
        "dataset",
        "modality",
        "raw_input_unit",
        "input_libraries_or_matrices",
        "n_pairs_before_idh",
        "n_pairs_excluded_idh",
        "excluded_pair_keys",
        "formal_pairs_threshold20",
        "n_cells_primary",
        "n_cells_recurrent",
        "clean_myeloid",
        "minimum_cells_per_condition",
        "processing_scope",
        "statistical_testing",
        "replication_assessment",
        "workflow",
    ]
    source_table = source_table.loc[:, safe_source_columns]

    stem = "Figure1A_compact_paired_design"
    source = write_source(source_table, "Figure1", stem)
    fig = plt.figure(figsize=(mm(183), mm(94)), facecolor="white")
    ax = fig.add_axes([0.025, 0.025, 0.95, 0.95])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.axis("off")

    def rounded_box(
        x0: float, x1: float, y0: float, y1: float, edge: str, fill: str, lw: float = 0.8
    ) -> None:
        """输入：x0、x1、y0、y1、edge、fill、lw。输出为文件或原位绘图设置；无返回值。"""
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                boxstyle="round,pad=0.02,rounding_size=0.7",
                facecolor=fill,
                edgecolor=edge,
                linewidth=lw,
            )
        )

    cohort_boxes = [(4.0, 47.5), (52.5, 96.0)]
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    neutral = "#717171"

    text = "#252525"

    for row, (x0, x1) in zip(source_table.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        rounded_box(x0, x1, 80.0, 97.0, color, pale, lw=1.0)
        ax.text(
            center,
            94.5,
            f"{row['dataset']}  |  snRNA-seq",
            fontsize=8.3,
            weight="bold",
            color=color,
            ha="center",
            va="top",
        )

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.text(
            center,
            89.3,
            f"{int(row['input_libraries_or_matrices'])} {row['raw_input_unit']}",
            fontsize=6.1,
            color=neutral,
            ha="center",
            va="top",
        )
        primary_x, recurrent_x = (center - 8.0, center + 8.0)
        rounded_box(primary_x - 4.9, primary_x + 4.9, 81.7, 86.7, color, "white", lw=0.9)
        rounded_box(recurrent_x - 5.5, recurrent_x + 5.5, 81.7, 86.7, color, color, lw=0.9)
        ax.text(primary_x, 84.2, "Primary", fontsize=6.3, color=text, ha="center", va="center")

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.text(recurrent_x, 84.2, "1st recurrence", fontsize=6.3, color="white", ha="center", va="center")
        ax.annotate(
            "",
            xy=(recurrent_x - 5.8, 84.2),
            xytext=(primary_x + 5.2, 84.2),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0),
        )
    for center in [sum(box) / 2 for box in cohort_boxes]:
        ax.annotate(
            "",
            xy=(center, 74.8),
            xytext=(center, 79.5),
            arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85),
        )

    # 由正式审计行生成设计图输入；细胞和配对计数仍取既定源表。
    workflow_boxes = [
        (
            8.0,
            92.0,
            64.0,
            74.0,
            "Independent all-cell reconstruction",
            "Library QC · doublet removal · basic lineage annotation",
        ),
        (
            8.0,
            92.0,
            50.0,
            60.0,
            "Cohort-wise pan-myeloid re-clustering",
            "2,000 HVGs · PCA · neighbor graph · Leiden",
        ),
        (
            6.0,
            94.0,
            36.0,
            46.0,
            "Marker and cluster-DEG review",
            "Remove lymphoid/tumour contaminants and low-quality cells",
        ),
    ]
    for index, (x0, x1, y0, y1, title, subtitle) in enumerate(workflow_boxes):
        rounded_box(x0, x1, y0, y1, "#858585", "white", lw=0.75)
        ax.text(50, y1 - 2.2, title, fontsize=7.3, weight="bold", ha="center", va="top")
        ax.text(50, y0 + 2.2, subtitle, fontsize=6.1, color=neutral, ha="center", va="bottom")
        if index < len(workflow_boxes) - 1:
            next_y1 = workflow_boxes[index + 1][3]
            ax.annotate(
                "",
                xy=(50, next_y1 + 0.5),
                xytext=(50, y0 - 0.5),
                arrowprops=dict(arrowstyle="-|>", color="#8C8C8C", lw=0.85),
            )
    analysis_y0, analysis_y1 = (14.5, 29.5)
    for row, (x0, x1) in zip(source_table.reset_index(drop=True).to_dict("records"), cohort_boxes):
        color = dataset_colors[row["dataset"]]
        pale = mpl.colors.to_rgba(color, 0.085)
        center = (x0 + x1) / 2
        source_x = 35.0 if center < 50 else 65.0
        ax.annotate(
            "",
            xy=(center, analysis_y1 + 0.6),
            xytext=(source_x, 35.5),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=0.9,
                connectionstyle="arc3,rad=0.08" if center < 50 else "arc3,rad=-0.08",
            ),
        )

        rounded_box(x0, x1, analysis_y0, analysis_y1, color, pale, lw=1.0)
        ax.text(
            center,
            26.8,
            f"{row['dataset']}  |  {int(row['formal_pairs_threshold20'])} paired patients",
            fontsize=7.1,
            weight="bold",
            color=color,
            ha="center",
            va="center",
        )
        ax.text(
            center,
            22.0,
            f"{int(row['clean_myeloid']):,} analyzed myeloid cells",
            fontsize=5.9,
            color=neutral,
            ha="center",
            va="center",
        )
        ax.text(
            center,
            17.2,
            "Raw-count pseudobulk · paired edgeR · Miller-IM GSEA",
            fontsize=5.9,
            color=text,
            ha="center",
            va="center",
        )
    ax.annotate(
        "",
        xy=(46.8, 10.8),
        xytext=(sum(cohort_boxes[0]) / 2, analysis_y0 - 0.5),
        arrowprops=dict(
            arrowstyle="-|>", color=dataset_colors["GSE174554"], lw=0.9, connectionstyle="arc3,rad=-0.08"
        ),
    )

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.annotate(
        "",
        xy=(53.2, 10.8),
        xytext=(sum(cohort_boxes[1]) / 2, analysis_y0 - 0.5),
        arrowprops=dict(
            arrowstyle="-|>", color=dataset_colors["GSE274546"], lw=0.9, connectionstyle="arc3,rad=0.08"
        ),
    )
    rounded_box(22.0, 78.0, 1.0, 11.8, "#159D88", mpl.colors.to_rgba("#159D88", 0.085), lw=1.0)
    ax.text(
        50,
        8.8,
        "Cross-cohort replication assessment",
        fontsize=7.2,
        weight="bold",
        color=text,
        ha="center",
        va="center",
    )
    ax.text(
        50,
        3.8,
        "Compare effect direction and nominal P after independent cohort tests",
        fontsize=5.9,
        color=neutral,
        ha="center",
        va="center",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")

    # 形成已计算结果的具名返回对象。
    result = panel_record("Figure1A", stem, "Independent paired-cohort myeloid workflow", source, pdf, png)
    return result


def label_position(coordinate_table: pd.DataFrame) -> tuple[float, float]:
    """输入二维坐标；返回保留原分位截取规则的中心标签位置。"""

    points = coordinate_table[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
    center = np.nanmedian(points, axis=0)
    distance = np.sqrt(((points - center) ** 2).sum(axis=1))
    core = points[distance <= np.nanquantile(distance, 0.6)]
    target = np.nanmedian(core, axis=0) if len(core) else center

    # 形成已计算结果的具名返回对象。
    result = (float(target[0]), float(target[1]))
    return result


def figure1_clean_myeloid_umap(dataset: str, panel: str) -> dict[str, object]:
    """输入队列及面板标识；读取既有髓系坐标和标签，保存源表及图，返回清单记录。"""

    # 确定当前面板的公开文件主名，后续源表与图使用同一主名。
    stem = f"{panel}_{dataset}_clean_myeloid_umap"
    myeloid_atlas_sources = {
        "GSE174554": root
        / "write/43_scgbm_style_clean_myeloid_reclustering/GSE174554/clean_myeloid_umap_source.csv.gz",
        "GSE274546": root
        / "write/43_scgbm_style_clean_myeloid_reclustering/GSE274546/clean_myeloid_umap_source.csv.gz",
    }

    umap_table = pd.read_csv(myeloid_atlas_sources[dataset])
    myeloid_state_order = [
        "Homeostatic microglia",
        "Transitional Mg/TAM",
        "Complement/HLA-II TAM",
        "Mg-inflammatory-high",
        "Inflammatory myeloid",
        "Scavenger/lipid TAM",
        "Hypoxia-glycolytic TAM",
        "IFN-responsive myeloid",
        "Heat-shock myeloid",
        "Cycling myeloid",
        "Monocyte-like",
        "cDC",
    ]

    unknown = sorted(set(umap_table["curated_myeloid_subtype"]) - set(myeloid_state_order))
    if unknown:
        raise ValueError(f"Unexpected reviewed myeloid states for {dataset}: {unknown}")
    source = write_source(umap_table, "Figure1", stem)

    # 按原画布尺寸绘制当前面板，以下位置和样式保持原图参数。
    fig = plt.figure(figsize=(mm(91), mm(88)), facecolor="white")
    ax = fig.add_axes([0.035, 0.205, 0.93, 0.64])
    counts = umap_table["curated_myeloid_subtype"].value_counts()
    x_span = float(umap_table["UMAP1"].max() - umap_table["UMAP1"].min())
    y_span = float(umap_table["UMAP2"].max() - umap_table["UMAP2"].min())

    point_size = 0.62 if len(umap_table) < 30000 else 0.34
    myeloid_state_colors = {
        "Homeostatic microglia": "#159D88",
        "Transitional Mg/TAM": "#73B6B2",
        "Complement/HLA-II TAM": "#3B6FB6",
        "Mg-inflammatory-high": "#D8574F",
        "Inflammatory myeloid": "#D98270",
        "Scavenger/lipid TAM": "#D18A34",
        "Hypoxia-glycolytic TAM": "#7B61A8",
        "IFN-responsive myeloid": "#4C9A61",
        "Heat-shock myeloid": "#B06B88",
        "Cycling myeloid": "#E2B33B",
        "Monocyte-like": "#8C6D5A",
        "cDC": "#6B778D",
    }

    for state in counts.sort_values(ascending=False).index:
        current_table = umap_table[umap_table["curated_myeloid_subtype"].eq(state)]
        ax.scatter(
            current_table["UMAP1"],
            current_table["UMAP2"],
            s=point_size,
            color=myeloid_state_colors[state],
            alpha=0.62,
            linewidths=0,
            rasterized=True,
        )
    myeloid_label_offsets = {
        "GSE174554": {
            "Mg-inflammatory-high": (-0.035, 0.025),
            "Inflammatory myeloid": (-0.02, 0.0),
            "Transitional Mg/TAM": (0.06, 0.0),
            "IFN-responsive myeloid": (-0.02, -0.005),
            "Homeostatic microglia": (0.01, -0.015),
            "Complement/HLA-II TAM": (0.035, -0.005),
            "Scavenger/lipid TAM": (0.035, 0.0),
        },
        "GSE274546": {
            "Mg-inflammatory-high": (0.035, 0.0),
            "Inflammatory myeloid": (0.04, -0.01),
            "Transitional Mg/TAM": (-0.075, -0.01),
            "Complement/HLA-II TAM": (0.075, 0.0),
            "Scavenger/lipid TAM": (-0.04, -0.025),
            "Hypoxia-glycolytic TAM": (0.055, -0.01),
            "Heat-shock myeloid": (0.055, 0.01),
            "Monocyte-like": (0.0, -0.02),
            "cDC": (0.02, -0.01),
        },
    }

    for state in myeloid_state_order:
        current_table = umap_table[umap_table["curated_myeloid_subtype"].eq(state)]
        if current_table.empty:
            continue
        largest_cluster = current_table["leiden_myeloid"].astype(str).value_counts().index[0]
        label_cell_table = current_table[current_table["leiden_myeloid"].astype(str).eq(largest_cluster)]
        x, y = label_position(label_cell_table)

        dx, dy = myeloid_label_offsets.get(dataset, {}).get(state, (0.0, 0.0))
        text_x = x + dx * x_span
        text_y = y + dy * y_span
        short_label = {
            "Homeostatic microglia": "Homeostatic Mg",
            "Complement/HLA-II TAM": "Complement / HLA-II",
            "Mg-inflammatory-high": "Mg-inflammatory",
            "Inflammatory myeloid": "Inflammatory",
            "Scavenger/lipid TAM": "Scavenger / lipid",
            "Hypoxia-glycolytic TAM": "Hypoxia / glycolysis",
            "IFN-responsive myeloid": "IFN-responsive",
            "Heat-shock myeloid": "Heat-shock",
            "Cycling myeloid": "Cycling",
        }.get(state, state)
        if dx or dy:
            ax.plot(
                [x, text_x], [y, text_y], color=myeloid_state_colors[state], lw=0.45, alpha=0.75, zorder=4
            )

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.text(
            text_x,
            text_y,
            short_label,
            fontsize=5.2,
            weight="bold",
            color=myeloid_state_colors[state],
            ha="center",
            va="center",
            path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
            zorder=5,
        )
    ax.set_xticks([])
    ax.set_yticks([])

    # 按原队列、基因或面板顺序逐项处理，保持循环内计算次序。
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    handles = [
        mpl.lines.Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=3.8,
            markerfacecolor=myeloid_state_colors[state],
            markeredgewidth=0,
            label=state,
        )
        for state in myeloid_state_order
        if state in counts.index
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=3,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.0,
        fontsize=5.25,
    )

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(
        0.05,
        0.905,
        f"{len(umap_table):,} clean myeloid cells  |  {umap_table['leiden_myeloid'].nunique()} clusters",
        fontsize=6.2,
        color="#717171",
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")
    result = panel_record(
        panel,
        stem,
        f"{dataset} independently integrated and DEG-reviewed clean-myeloid UMAP",
        source,
        pdf,
        png,
    )
    return result


# 原模块轴线长度为2.4，保留与公共工具2.5不同的既有图形参数。


def clean_axis(ax: mpl.axes.Axes, *, left: bool = True, bottom: bool = True) -> None:
    """输入Axes及边框开关；原位设置坐标轴外观，不更改数据；无返回值。"""

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(width=0.6, length=2.4)


def figure1_rank_rank() -> dict[str, object]:
    """读取Step41逐基因排名曲线；保存合并排名源表及密度图，返回清单记录。"""

    # 读取既有两队列逐基因排名，按gene一对一合并。
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    curve_table = pd.read_csv(step41 / "Figure1/source_data/raw20_gsea_curves.csv.gz")
    left_table = curve_table[curve_table["dataset"].eq("GSE174554")][
        ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE174554", "rank_stat": "rank_stat_GSE174554", "hit": "program_GSE174554"}
    )
    right_table = curve_table[curve_table["dataset"].eq("GSE274546")][
        ["gene", "rank", "rank_stat", "hit"]
    ].rename(
        columns={"rank": "rank_GSE274546", "rank_stat": "rank_stat_GSE274546", "hit": "program_GSE274546"}
    )
    rank_table_total = left_table.merge(right_table, on="gene", how="inner", validate="one_to_one")

    # 构建或筛选两队列基因排名，维持原排序和合并键。
    rank_table_total["miller_raw20"] = rank_table_total[["program_GSE174554", "program_GSE274546"]].any(
        axis=1
    )
    rank_table_total["density"] = rank_density(rank_table_total)
    program_table = rank_table_total[rank_table_total["miller_raw20"]].copy()
    program_table["label_score"] = (
        program_table["rank_stat_GSE174554"].abs() + program_table["rank_stat_GSE274546"].abs()
    )
    label_genes = set(program_table.nlargest(min(6, len(program_table)), "label_score")["gene"])

    # 构建或筛选两队列基因排名，维持原排序和合并键。
    rank_table_total["labelled"] = rank_table_total["gene"].isin(label_genes)
    stem = "Figure1D_cross_cohort_all_gene_rank_density"
    source = write_source(rank_table_total, "Figure1", stem)
    fig = plt.figure(figsize=(mm(92), mm(83)), facecolor="white")
    ax = fig.add_axes([0.16, 0.14, 0.76, 0.69])
    background_table = rank_table_total[~rank_table_total["miller_raw20"]].sort_values("density")

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.scatter(
        background_table["rank_stat_GSE174554"],
        background_table["rank_stat_GSE274546"],
        c=background_table["density"],
        cmap="Greys",
        s=2.1,
        alpha=0.7,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        program_table["rank_stat_GSE174554"],
        program_table["rank_stat_GSE274546"],
        s=24,
        color="#D8574F",
        edgecolor="white",
        linewidth=0.65,
        zorder=4,
    )
    label_offsets = {
        "CCL3": (3, 7),
        "CCL4": (-26, 7),
        "CH25H": (5, 5),
        "SGK1": (-25, -10),
        "FOLR2": (5, 7),
        "BHLHE41": (5, -8),
    }
    for _, row in program_table[program_table["gene"].isin(label_genes)].iterrows():
        offset = label_offsets.get(row["gene"], (3, 3))
        ax.annotate(
            row["gene"],
            (row["rank_stat_GSE174554"], row["rank_stat_GSE274546"]),
            xytext=offset,
            textcoords="offset points",
            fontsize=5.7,
            color="#A63D38",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
        )
    ax.axhline(0, color="#A4A4A4", lw=0.55)

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.axvline(0, color="#A4A4A4", lw=0.55)
    limit = (
        float(np.nanmax(np.abs(rank_table_total[["rank_stat_GSE174554", "rank_stat_GSE274546"]].to_numpy())))
        * 1.04
    )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.set_xlabel("GSE174554 paired rank statistic")
    ax.set_ylabel("GSE274546 paired rank statistic")
    ax.grid(color="#E6E6E6", lw=0.42, zorder=0)
    clean_axis(ax)
    rho, _ = spearmanr(rank_table_total["rank_stat_GSE174554"], rank_table_total["rank_stat_GSE274546"])

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    neutral = "#717171"

    ax.text(
        0.02,
        0.98,
        f"All shared genes: {len(rank_table_total):,}\nSpearman rho = {rho:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.5,
        color=neutral,
    )
    ax.text(
        0.98,
        0.03,
        f"Miller-IM: {len(program_table)} measured genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        color="#A63D38",
    )
    fig.text(0.16, 0.965, "Cross-cohort all-gene concordance", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.16,
        0.915,
        "Transcriptome density with Miller-IM genes highlighted",
        fontsize=6.4,
        color=neutral,
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure1", stem, "Step42 recurrent GBM visual story")

    # 形成已计算结果的具名返回对象。
    result = panel_record(
        "Figure1D", stem, "All-gene rank-stat density with Miller-IM genes highlighted", source, pdf, png
    )
    return result


def figure2_raw20_umap(dataset: str, panel: str, vmin: float, vmax: float) -> dict[str, object]:
    """输入队列、面板标识和色标界限；保存既有raw20连续分值UMAP及源表，返回清单记录。"""

    # 确定当前面板的公开文件主名，后续源表与图使用同一主名。
    stem = f"{panel}_{dataset}_raw20_terrain"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    umap_table = pd.read_csv(step41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz")
    source = write_source(umap_table, "Figure2", stem)
    ordered_table = umap_table.sort_values("raw20_score")
    fig = plt.figure(figsize=(mm(88), mm(78)), facecolor="white")

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax = fig.add_axes([0.045, 0.035, 0.8, 0.78])
    points = ax.scatter(
        ordered_table["UMAP1"],
        ordered_table["UMAP2"],
        c=ordered_table["raw20_score"],
        s=0.55,
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        alpha=0.94,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.012)
    cax = fig.add_axes([0.88, 0.2, 0.023, 0.48])
    cb = fig.colorbar(points, cax=cax)
    cb.ax.set_title("Miller-IM\nscore", fontsize=6.2, pad=4, loc="center")

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    cb.ax.tick_params(labelsize=6.0, length=2)
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    fig.text(0.05, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=dataset_colors[dataset])
    fig.text(0.05, 0.905, "Continuous program terrain", fontsize=6.6, color="#717171", va="top")
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(panel, stem, f"{dataset} continuous Miller-IM score terrain", source, pdf, png)

    return result


def ridge_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """输入数值向量和密度网格；返回峰值归一化KDE，少于两值或零方差时返回零向量。"""

    values = values[np.isfinite(values)]
    if len(values) < 2 or np.isclose(np.std(values), 0):
        result = np.zeros_like(grid)
        return result
    density = gaussian_kde(values)(grid)
    peak = density.max()
    result = density / peak if peak > 0 else density

    return result


def figure2_ridgeline(dataset: str, panel: str, x_limits: tuple[float, float]) -> dict[str, object]:
    """输入队列、面板标识和横轴界限；保存患者簇分值山脊图及源表，返回清单记录。"""

    # 确定当前面板的公开文件主名，后续源表与图使用同一主名。
    stem = f"{panel}_{dataset}_cluster_ridgeline"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    raw_table = pd.read_csv(step41 / "Figure2/source_data/fig2_patient_cluster_interval_source.csv")
    patient_table = raw_table[raw_table["dataset"].eq(dataset) & raw_table["eligible"].astype(bool)].copy()
    label_table = (
        pd.read_csv(
            step41 / f"Figure2/source_data/fig2_{dataset}_raw20_umap_source.csv.gz",
            usecols=["paired_cluster", "cluster_label"],
        )
        .dropna()
        .drop_duplicates()
    )
    if dataset == "GSE274546":

        # 读取此步骤的既定输入，保留原行序、数据集名称及公开字段。
        annotation_table = pd.read_csv(
            root
            / (
                "write/36_gse274546_raw_independent_reannotation/04_paired_patient"
                "_harmony_clusters/paired_cluster_descriptive_annotation.csv"
            )
        )
        label_table = annotation_table[["cluster", "descriptive_state"]].rename(
            columns={"cluster": "paired_cluster", "descriptive_state": "cluster_label"}
        )
        label_table["cluster_label"] = label_table.apply(
            lambda row: f"C{int(row['paired_cluster'])} {row['cluster_label']}", axis=1
        )

    # 沿用既定标签及位置规则，避免改变原图的标注内容。
    label_table["paired_cluster"] = pd.to_numeric(label_table["paired_cluster"]).astype(int).astype(str)
    label_map = label_table.set_index("paired_cluster")["cluster_label"].to_dict()
    patient_table["paired_cluster"] = pd.to_numeric(patient_table["paired_cluster"]).astype(int).astype(str)
    patient_table["cluster_label_display"] = patient_table["paired_cluster"].map(label_map)
    if patient_table["cluster_label_display"].isna().any():
        raise ValueError(f"Missing full cluster labels for {dataset}")

    # 汇总当前队列的实际记录数及极值，维持原输出列。
    summary_table = (
        patient_table.groupby(["paired_cluster", "cluster_label_display"], observed=True)[
            "pooled_patient_raw20"
        ]
        .agg(["median", "count"])
        .reset_index()
        .sort_values("median")
    )
    order = summary_table["paired_cluster"].astype(str).tolist()
    patient_table["ridge_order"] = patient_table["paired_cluster"].map(
        {cluster: i for i, cluster in enumerate(order)}
    )
    source = write_source(patient_table, "Figure2", stem)
    fig = plt.figure(figsize=(mm(100), mm(91)), facecolor="white")

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax = fig.add_axes([0.42, 0.13, 0.55, 0.7])
    grid = np.linspace(x_limits[0], x_limits[1], 300)
    dataset_colors = {"GSE174554": "#3B5B92", "GSE274546": "#159D88"}

    color = dataset_colors[dataset]
    tick_labels = []
    neutral = "#717171"

    for i, cluster in enumerate(order):
        current_table = patient_table[patient_table["paired_cluster"].eq(cluster)]
        values = current_table["pooled_patient_raw20"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        baseline = float(i)
        ax.fill_between(grid, baseline, baseline + density * 0.78, color=color, alpha=0.36, linewidth=0)

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.plot(grid, baseline + density * 0.78, color=color, lw=0.9)
        jitter = np.random.default_rng(20260713 + i).uniform(-0.055, 0.055, len(values))
        ax.scatter(values, baseline - 0.1 + jitter, s=5.5, color="#252525", alpha=0.55, linewidths=0, zorder=3)
        ax.scatter(
            np.median(values),
            baseline - 0.1,
            marker="D",
            s=15,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        tick_labels.append(str(current_table["cluster_label_display"].iloc[0]))

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.text(
            0.985,
            baseline + 0.49,
            f"n={len(values)}",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=5.5,
            color=neutral,
        )

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.set_yticks(range(len(order)), tick_labels)
    ax.set_ylim(-0.35, len(order) - 0.02 + 0.82)
    ax.set_xlim(*x_limits)
    ax.set_xlabel("Patient-level pooled Miller-IM score")
    ax.grid(axis="x", color="#E6E6E6", lw=0.45)

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.tick_params(axis="y", length=0, labelsize=5.2, pad=2)
    clean_axis(ax, left=False, bottom=True)
    fig.text(0.08, 0.965, dataset, fontsize=8.7, weight="bold", va="top", color=color)
    fig.text(
        0.08, 0.913, "Miller-IM score across unsupervised clusters", fontsize=6.5, color=neutral, va="top"
    )
    fig.text(
        0.08, 0.87, "Each dot is one eligible patient-cluster unit", fontsize=6.2, color=neutral, va="top"
    )

    # 保存原尺寸PDF母版；数据图的PNG仅由对应PDF派生。
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(panel, stem, f"{dataset} patient-level cluster ridgeline", source, pdf, png)
    return result


def figure2_delta_raincloud() -> dict[str, object]:
    """读取既有GSE278456患者内状态差值；保存原状态顺序的雨云图及源表，返回清单记录。"""

    # 确定当前面板的公开文件主名，后续源表与图使用同一主名。
    stem = "Figure2E_GSE278456_patient_delta_raincloud"
    step41 = root / "write/41_mg_inflammatory_sci_rebuild"

    delta_table = pd.read_csv(step41 / "Figure2/source_data/fig2_gse278456_all_state_patient_deltas.csv")
    state_order = ["MCG1", "MCG2", "MAC1", "MAC2", "E-MDSC", "M-MDSC"]

    delta_table = delta_table[delta_table["state"].isin(state_order)].copy()
    delta_table["state_order"] = delta_table["state"].map({state: i for i, state in enumerate(state_order)})
    source = write_source(delta_table, "Figure2", stem)

    # 按原画布尺寸绘制当前面板，以下位置和样式保持原图参数。
    fig = plt.figure(figsize=(mm(118), mm(82)), facecolor="white")
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.67])
    x_min = min(-0.05, float(delta_table["delta"].min()) - 0.06)
    x_max = float(delta_table["delta"].max()) + 0.08
    grid = np.linspace(x_min, x_max, 320)

    # 沿用原随机种子及调用位置，不调整抽样或抖动顺序。
    neutral = "#717171"

    state_colors = {
        "MCG1": "#D8574F",
        "MCG2": "#159D88",
        "MAC1": "#277DA1",
        "MAC2": "#72B7B2",
        "E-MDSC": "#E09F3E",
        "M-MDSC": "#8C6D5A",
    }

    text = "#252525"

    for i, state in enumerate(state_order):
        values = delta_table.loc[delta_table["state"].eq(state), "delta"].to_numpy(dtype=float)
        density = ridge_density(values, grid)
        color = state_colors[state]
        ax.fill_between(grid, i, i + density * 0.34, color=color, alpha=0.4, linewidth=0)
        ax.plot(grid, i + density * 0.34, color=color, lw=0.9)

        # 沿用原随机种子及调用位置，不调整抽样或抖动顺序。
        jitter = np.random.default_rng(20260713 + 100 + i).uniform(-0.3, -0.08, len(values))
        ax.scatter(
            values, i + jitter, s=14, color=color, edgecolor="white", linewidth=0.35, alpha=0.8, zorder=3
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.hlines(i - 0.02, q1, q3, color=text, lw=2.0, zorder=4)
        ax.scatter(
            median, i - 0.02, s=20, marker="D", color=text, edgecolor="white", linewidth=0.45, zorder=5
        )

        # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
        ax.text(x_max, i, f"n={len(values)}", ha="right", va="center", fontsize=6.0, color=neutral)
    ax.axvline(0, color="#929292", lw=0.7, ls="--")
    ax.set_yticks(range(len(state_order)), state_order)
    ax.set_ylim(-0.45, len(state_order) - 0.45)
    ax.invert_yaxis()

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Within-patient Miller-IM score delta (state mean - other states mean)")
    ax.grid(axis="x", color="#E6E6E6", lw=0.45)
    ax.tick_params(axis="y", length=0)
    clean_axis(ax, left=False, bottom=True)

    # 设置该面板的坐标、标记和文字，保留原图面表达与显示范围。
    fig.text(0.17, 0.965, "GSE278456 patient-level state contrasts", fontsize=8.7, weight="bold", va="top")
    fig.text(
        0.17,
        0.91,
        "Cloud: density | dots: patients | diamond and line: median and IQR",
        fontsize=6.4,
        color=neutral,
        va="top",
    )
    pdf, png = save_panel(fig, fig_root / "Figure2", stem, "Step42 recurrent GBM visual story")
    result = panel_record(
        "Figure2E", stem, "GSE278456 patient delta raincloud across published states", source, pdf, png
    )
    return result


# 写出当前源表并检查非空文件，保持原列名、排序和压缩参数。
if __name__ == "__main__":

    # 沿用原随机种子及调用位置，不调整抽样或抖动顺序。
    np.random.seed(20260713)
    configure_style()
    (fig_root / "Figure1").mkdir(parents=True, exist_ok=True)
    (write_root / "Figure1" / "source_data").mkdir(parents=True, exist_ok=True)
    panel_record_1 = figure1_design()

    figure1_records = [panel_record_1]
    figure_name = "Figure1"
    records = figure1_records
    out = write_root / figure_name
    out.mkdir(parents=True, exist_ok=True)

    # 将当前具名记录整理成表；保持原记录顺序及每行的统计单位。
    output_table = pd.DataFrame(records)
    csv_path = out / f"{figure_name}_panel_manifest.csv"
    output_table.to_csv(csv_path, index=False)
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty output: {csv_path.name}")
