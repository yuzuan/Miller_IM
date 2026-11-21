#!/usr/bin/env python3
# ======================================================================
# 00a_functions.py —— Miller-IM公共函数库
# 项目：公开配对复发GBM分析。
#
# 目的：
# 统一多个脚本使用的配色、PDF保存、PNG渲染与预览排版。
# 保留原Matplotlib实现，避免语言替换改变既有图件。
# 本函数脚本不用单独浏览学习，在其它分析脚本调用到它时再学。
# 编辑器中点击调用处函数名，可以跳到对应定义。
# 函数定义的查找靠编辑器跳转，不靠文件头清单。

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
    """输入毫米数，返回英寸数；25.4为毫米到英寸的精确换算。"""
    inches = value / 25.4
    return inches


# 颜色由多个独立绘图入口读取；名称是原图中的组别，不重新分组。
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
    """无输入；设置原公开图的字体、线宽和背景；无返回值。"""
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
    """输入画布毫米尺寸及Matplotlib参数，返回白底Figure。"""
    figure = plt.figure(figsize=(mm(width_mm), mm(height_mm)), facecolor="white", **kwargs)
    return figure


def clean_axis(ax, *, keep_left: bool = True, keep_bottom: bool = True) -> None:
    """输入Axes及边框开关，修改原图轴线样式；无返回值。"""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)
    ax.tick_params(width=0.6, length=2.5)


def render_pdf_png(pdf_path: Path, png_path: Path, dpi: int = 600) -> None:
    """将单页PDF母版渲染为PNG；缺文件、渲染失败或空图片时停止；无返回值。"""
    if not pdf_path.is_file() or pdf_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PDF: {pdf_path.name}")

    # Poppler从已保存PDF读取图形；600 dpi沿用既有期刊图分辨率。
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["pdftoppm", "-singlefile", "-png", "-r", str(dpi), str(pdf_path), str(png_path.with_suffix(""))],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    # 检查实际产物，防止仅有文件名而没有可读图片。
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty PNG: {png_path.name}")
    with Image.open(png_path) as image:
        if min(image.size) <= 0:
            raise ValueError(f"Invalid PNG dimensions: {png_path.name}")
        image.verify()


def save_panel(fig, out_dir: Path, stem: str, creator: str) -> tuple[Path, Path]:
    """输入Figure、目录、主名及原作者标记，返回已核验的PDF和PNG路径。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{stem}.pdf"

    # 先绘制画布再保存PDF，沿用原尺寸、白底与PDF元数据。
    fig.canvas.draw()
    fig.savefig(
        pdf_path,
        facecolor="white",
        edgecolor="none",
        metadata={"Creator": creator, "CreationDate": None, "ModDate": None},
    )

    # PNG只从刚生成的母版转换，两个格式使用同一图形。
    png_path = out_dir / f"{stem}.png"
    render_pdf_png(pdf_path, png_path)
    plt.close(fig)
    output_paths = (pdf_path, png_path)
    return output_paths


def save_figure(fig, out_dir: Path, stem: str) -> None:
    """兼容原Step41公共调用；保留原Creator及无返回值接口。"""
    save_panel(fig, out_dir, stem, "Step41 recurrent GBM figure rebuild")


def sha256(path: Path) -> str:
    """输入已有文件路径，返回原发布清单所用的摘要字符串。"""
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
    """输入单个面板及三份产物路径，返回原发布清单的九个字段。"""
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
    """按原候选顺序读取预览字体，返回Pillow字体；沿用原缺字体处理。"""
    font_name = "Arial Bold.ttf" if bold else "Arial.ttf"
    candidates = [Path("/System/Library/Fonts/Supplemental") / font_name, Path("/System/Library/Fonts/Helvetica.ttc")]

    # 字体只是预览依赖；依次尝试原Arial、Helvetica及Pillow默认字体。
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
    """将已有面板等比居中放入预览矩形，不改变矩形布局；无返回值。"""
    left, top, right, bottom = box
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        fitted = ImageOps.contain(white.convert("RGB"), (right - left, bottom - top), Image.Resampling.LANCZOS)

    # 留白在左右和上下各平分，沿用原预览定位。
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def enrichment_curve(deg_table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """输入含gene/logFC/F的差异表与固定基因列表；返回逐基因排名、命中及累积富集表。"""

    # 构建或筛选两队列基因排名，维持原排序和合并键。
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
    """输入两队列排名统计表及分箱数；返回每行所在二维箱计数的log1p值。"""

    x = rank_table["rank_stat_GSE174554"].to_numpy(dtype=float)
    y = rank_table["rank_stat_GSE274546"].to_numpy(dtype=float)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    xi = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    yi = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)

    # 形成已计算结果的具名返回对象。
    result = np.log1p(counts[xi, yi])
    return result
