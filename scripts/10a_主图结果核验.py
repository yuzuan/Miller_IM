#!/usr/bin/env python3
"""Figure 1–5 正式产物的一致性与命名核验，保留原 Python API。

目的：逐条核对既有图源、统计量、命名、文件摘要与输出格式；不重做分析或出图。
输入：write/51_figure1_final_reorganized/、50_figure2_two_panel_final/、
    55_figure3_visual_polish/、54_figure4_proteome_rebuild/、56_figure5_pd1_rebuild/。
    各 Figure 目录提供 panel_manifest.csv、legend.md 和 source_data/ 中的原命名文件。
输出：write/57_main_figure_consistency_audit/main_figure_consistency_checks.csv
    和同目录 MAIN_FIGURE_AUDIT.md；任一既有核验失败时返回非零状态。
运行：python scripts/10a_主图结果核验.py。
固定数值均为既有发布结果的核验目标；字段、路径、容差及报表字段保持兼容。
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
sha256 = figure_style.sha256


def main() -> None:
    """按交付文件、Figure 1–5、汇总报表顺序核验；失败抛错，无返回值。"""
    audit_delivery()
    audit_figure1()
    audit_figure2()
    audit_figure3()
    audit_figure4()
    audit_figure5()
    write_report()


def audit_delivery() -> None:
    """核对 manifest 指向的图源、图件与图注，追加核验记录；无返回值。"""

    # 每个条目对应正式交付目录和 panel 顺序；旧路径是现有图表的读取接口。
    root = Path(__file__).resolve().parents[1]
    figures = {
        "Figure1": {
            "manifest": root / "write/51_figure1_final_reorganized/Figure1/Figure1_panel_manifest.csv",
            "legend": root / "write/51_figure1_final_reorganized/Figure1/Figure1_legend.md",
            "panels": list("ABCDEF"),
        },
        "Figure2": {
            "manifest": root / "write/50_figure2_two_panel_final/Figure2/Figure2_panel_manifest.csv",
            "legend": root / "write/50_figure2_two_panel_final/Figure2/Figure2_legend.md",
            "panels": ["2A", "2B"],
        },
        "Figure3": {
            "manifest": root / "write/55_figure3_visual_polish/Figure3/Figure3_panel_manifest.csv",
            "legend": root / "write/55_figure3_visual_polish/Figure3/Figure3_legend.md",
            "panels": list("ABCD"),
        },
        "Figure4": {
            "manifest": root / "write/54_figure4_proteome_rebuild/Figure4/Figure4_panel_manifest.csv",
            "legend": root / "write/54_figure4_proteome_rebuild/Figure4/Figure4_legend.md",
            "panels": list("ABC"),
        },
        "Figure5": {
            "manifest": root / "write/56_figure5_pd1_rebuild/Figure5/Figure5_panel_manifest.csv",
            "legend": root / "write/56_figure5_pd1_rebuild/Figure5/Figure5_legend.md",
            "panels": list("ABC"),
        },
    }
    forbidden = re.compile(r"raw20|Miller\s+raw|Mg-inflammatory|Mg-IM|pharmacodynamic", re.IGNORECASE)
    for figure, spec in figures.items():
        manifest_table = pd.read_csv(spec["manifest"])
        observed_panels = manifest_table["panel"].astype(str).tolist()
        add(figure, "panel mapping", observed_panels == spec["panels"], str(observed_panels))

        # 一行对应一个 panel；多来源仍按原分号分隔，不改变相对/绝对路径解释。
        for row in manifest_table.to_dict("records"):
            panel = str(row["panel"])
            sources = [Path(item) for item in str(row["source"]).split(";")]
            pdf = Path(row["pdf"])
            png = Path(row["png"])
            files_exist = all(path.exists() for path in [*sources, pdf, png])
            add(figure, f"{panel} files exist", files_exist, "source/PDF/PNG")
            if "source_sha256" in row and pd.notna(row["source_sha256"]):
                source_matches = len(sources) == 1 and sha256(sources[0]) == row["source_sha256"]
                add(figure, f"{panel} source hash", source_matches, sources[0].name)
                add(figure, f"{panel} PDF hash", sha256(pdf) == row["pdf_sha256"], pdf.name)
                add(figure, f"{panel} PNG hash", sha256(png) == row["png_sha256"], png.name)

            # DPI 两轴均须与 600 相差小于 1；字体检查沿用 pdffonts 的 Type 3 判据。
            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0, 0))
            passed_dpi = len(dpi) >= 2 and all(abs(float(value) - 600) < 1 for value in dpi[:2])
            add(figure, f"{panel} PNG 600 dpi", passed_dpi, str(dpi))
            font_result = subprocess.run(["pdffonts", str(pdf)], check=True, capture_output=True, text=True)
            add(figure, f"{panel} PDF no Type3", "Type 3" not in font_result.stdout, pdf.name)

            # 只扫描可见 PDF 文字；内部文件名中的 raw20 仍属兼容标识。
            content = pdf_text(pdf)
            match = forbidden.search(content)
            detail = "no legacy public name" if match is None else match.group(0)
            add(figure, f"{panel} visible nomenclature", match is None, detail)

        # 同一命名规则覆盖图注与 manifest 的 message 字段。
        legend = Path(spec["legend"]).read_text(encoding="utf-8")
        match = forbidden.search(legend)
        detail = "Miller-IM public terminology" if match is None else match.group(0)
        add(figure, "legend nomenclature", match is None, detail)
        misleading_fdr = re.search(r"targeted\s+(?:fgsea\s+)?FDR", legend, re.IGNORECASE)
        detail = (
            "nominal P used for single-set enrichment"
            if misleading_fdr is None else misleading_fdr.group(0)
        )
        add(figure, "legend program-level P labeling", misleading_fdr is None, detail)
        if "message" in manifest_table.columns:
            messages = "\n".join(manifest_table["message"].fillna("").astype(str))
            match = forbidden.search(messages)
            detail = "Miller-IM public terminology" if match is None else match.group(0)
            add(figure, "manifest message nomenclature", match is None, detail)

    # Figure 1 图注必须包含原来的首次完整定义。
    figure1_legend = Path(figures["Figure1"]["legend"]).read_text(encoding="utf-8")
    definition = "20-gene Miller-derived inflammatory microglial program (Miller-IM program)"
    add("Cross-figure", "first formal definition", definition in figure1_legend, definition)


# 多个 audit_* API 共用同一记录列表；保留原来累加后统一写报表的行为。
check_records: list[dict[str, str]] = []


def add(section: str, check: str, passed: bool, detail: str) -> None:
    """按原四字段追加一项核验记录；无返回值。"""
    check_records.append(
        {
            "section": section,
            "check": check,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        }
    )


def pdf_text(path: Path) -> str:
    """调用 pdftotext 返回全部文字；外部命令失败时直接抛错。"""
    result = subprocess.run(["pdftotext", str(path), "-"], check=True, capture_output=True, text=True)
    content = result.stdout
    return content


def audit_figure1() -> None:
    """核对 Figure 1 六个 panel 的既有源表，追加结果；无返回值。"""

    # A 每行一个队列；保留输入规模、分析细胞、合格配对及本地路径检查。
    base = Path(__file__).resolve().parents[1] / "write/51_figure1_final_reorganized/Figure1/source_data"
    design_table = pd.read_csv(base / "Figure1A_paired_cohort_design_source.csv")
    design_table = design_table.set_index("dataset")
    design_text = design_table.reset_index().astype(str).to_csv(index=False)
    local_roots = ("/" + "Users/", "/" + "Volumes/")
    add(
        "Figure1", "A source contains no local paths",
        all(root not in design_text for root in local_roots), "public-safe fields only",
    )
    input_counts = design_table["input_libraries_or_matrices"].to_dict()
    add("Figure1", "A cohort inputs", input_counts == {"GSE174554": 91, "GSE274546": 111}, str(input_counts))
    cell_counts = design_table["clean_myeloid"].to_dict()
    add(
        "Figure1", "A analyzed myeloid cells",
        cell_counts == {"GSE174554": 12489, "GSE274546": 55568}, str(cell_counts),
    )
    pair_counts = design_table["formal_pairs_threshold20"].to_dict()
    add("Figure1", "A paired patients", pair_counts == {"GSE174554": 17, "GSE274546": 45}, str(pair_counts))

    # B/C 的每行曲线点重复携带统计量；仍用首行和原有绝对容差。
    expected_gsea = {
        "GSE174554": (2.4290789842056, 4.00136234233155e-7, 17),
        "GSE274546": (1.79932425427, 0.00397813985660199, 45),
    }
    for panel, cohort in [("B", "GSE174554"), ("C", "GSE274546")]:
        gsea_table = pd.read_csv(base / f"Figure1{panel}_{cohort}_raw20_gsea_source.csv.gz")
        observed = (gsea_table["NES"].iat[0], gsea_table["nominal_P"].iat[0], int(gsea_table["n_pairs"].iat[0]))
        target = expected_gsea[cohort]
        passed = close(observed[0], target[0]) and close(observed[1], target[1]) and observed[2] == target[2]
        detail = f"NES={observed[0]:.6f}; nominal P={observed[1]:.6g}; n={observed[2]}"
        add("Figure1", f"{panel} GSEA statistics", passed, detail)

    # D 一行一个全基因排序记录；冻结 19/7/5 数量和五个显示标签。
    rank_table = pd.read_csv(base / "Figure1D_shared_raw20_outliers_source.csv")
    counts = (
        len(rank_table), int(rank_table["miller_raw20"].sum()), int(rank_table["shared_leading_edge"].sum()),
        int(rank_table["display_label"].fillna("").ne("").sum()),
    )
    observed_labels = set(rank_table.loc[rank_table["display_label"].fillna("").ne(""), "display_label"])
    rank_passed = (
        len(rank_table) > 9000
        and rank_table["gene"].nunique() == len(rank_table)
        and counts[1:] == (19, 7, 5)
        and observed_labels == {"CCL4", "CH25H", "FOLR2", "PDK4", "SGK1"}
    )
    add("Figure1", "D all-gene/shared-gene counts", rank_passed, f"{counts}; labels={sorted(observed_labels)}")

    # E 的行粒度为基因×队列；先用未取位的 FDR 计数，再格式化报表。
    effect_table = pd.read_csv(base / "Figure1E_shared_leading_edge_gene_effects_source.csv")
    significant_count = int((effect_table["FDR"] < 0.05).sum())
    passed = len(effect_table) == 14 and effect_table["gene"].nunique() == 7 and significant_count == 0
    detail = f"rows={len(effect_table)}; genes={effect_table['gene'].nunique()}; FDR<0.05={significant_count}"
    add("Figure1", "E gene-by-cohort tests", passed, detail)

    # F 为每队列的 LOPO 汇总；保留每项折数与最大 nominal P 的原容差。
    lopo_table = pd.read_csv(base / "Figure1F_leave_one_patient_out_summary.csv")
    lopo_table = lopo_table.set_index("dataset")
    passed = (
        int(lopo_table.loc["GSE174554", "positive_NES_folds"]) == 17
        and int(lopo_table.loc["GSE274546", "positive_NES_folds"]) == 45
        and int(lopo_table.loc["GSE174554", "nominal_P_lt_0_05_folds"]) == 17
        and int(lopo_table.loc["GSE274546", "nominal_P_lt_0_05_folds"]) == 45
        and close(lopo_table.loc["GSE274546", "max_nominal_P"], 0.0311, atol=5e-5)
    )
    add("Figure1", "F LOPO stability", passed, "positive 17/17, 45/45; nominal P<0.05 17/17, 45/45; max 0.0311")


def close(value: float, target: float, atol: float = 1e-9) -> bool:
    """返回原有数值比较结果；默认 atol=1e-9、rtol=0，显式参数仍可覆盖 atol。"""
    passed = bool(np.isclose(float(value), float(target), atol=atol, rtol=0))
    return passed


def audit_figure2() -> None:
    """核对 Figure 2 状态、均值/FDR、marker QC 和重现人数；无返回值。"""

    # A 每行一个状态中的患者，组内重复携带状态汇总；先核对九组患者计数。
    base = Path(__file__).resolve().parents[1] / "write/50_figure2_two_panel_final/Figure2/source_data"
    waterfall_table = pd.read_csv(base / "Figure2_state_patient_waterfall_source.csv")
    expected = {
        "MCG1": (18, 17), "MCG2": (20, 16), "MCG3": (5, 1),
        "MCG4": (12, 1), "MCG5": (10, 3), "MAC1": (3, 0),
        "MAC2": (12, 5), "E-MDSC": (17, 0), "M-MDSC": (11, 0),
    }
    observed = {}
    for state, state_table in waterfall_table.groupby("state", sort=False):
        observed[state] = (int(state_table["n_patients"].iat[0]), int(state_table["n_positive"].iat[0]))
    add("Figure2", "A state patient counts", len(waterfall_table) == 108 and observed == expected, str(observed))

    # 同一状态还须满足原均值和九状态 FDR；集合不一致时保留原短路失败行为。
    expected_statistics = {
        "MCG1": (0.2785499722190481, 6.866455078125e-05),
        "MCG2": (0.0977333217439854, 0.0040615081787109),
        "MCG3": (-0.0847913201303971, 0.1607142857142857),
        "MCG4": (-0.1940185334308926, 0.002197265625),
        "MCG5": (-0.1221195292707085, 0.10546875),
        "MAC1": (-0.2835956082055430, 0.28125),
        "MAC2": (-0.0329282403218937, 0.4990234375),
        "E-MDSC": (-0.3026058457897609, 6.866455078125e-05),
        "M-MDSC": (-0.3471333914445482, 0.002197265625),
    }
    observed_statistics = {}
    for state, state_table in waterfall_table.groupby("state", sort=False):
        observed_statistics[state] = (
            float(state_table["mean_delta"].iat[0]), float(state_table["fdr_9_states"].iat[0]),
        )

    passed_statistics = set(observed_statistics) == set(expected_statistics)
    if passed_statistics:
        for state, target in expected_statistics.items():
            matches = (
                close(observed_statistics[state][0], target[0]) and close(observed_statistics[state][1], target[1])
            )
            if not matches:
                passed_statistics = False
                break

    # 取位只用于报表显示，不影响上面的核验。
    rounded_statistics = {}
    for state, value in observed_statistics.items():
        rounded_statistics[state] = (round(value[0], 6), round(value[1], 6))
    add("Figure2", "A state means and FDRs", passed_statistics, str(rounded_statistics))
    discordant = set(
        waterfall_table.loc[~waterfall_table["score_pseudobulk_mean_direction_concordant"].astype(bool), "state"]
    )
    add("Figure2", "A score-pseudobulk boundary", discordant == {"MCG3"}, str(discordant))

    # B 一行一个基因×状态；原始 marker 的状态峰值用于分类 QC。
    gene_table = pd.read_csv(base / "Figure2_raw20_marker_double_heatmap_raw20_source.csv")
    marker_table = pd.read_csv(base / "Figure2_raw20_marker_double_heatmap_marker_qc_source.csv")
    passed = len(gene_table) == 171 and gene_table["gene"].nunique() == 19 and gene_table["state"].nunique() == 9
    detail = f"{len(gene_table)} rows; {gene_table['gene'].nunique()} genes; {gene_table['state'].nunique()} states"
    add("Figure2", "B held-out gene block", passed, detail)
    marker_peak_table = marker_table.pivot(
        index="tested_marker_set", columns="assigned_state", values="within_marker_set_z"
    )
    passed = len(marker_table) == 81 and all(marker_peak_table.idxmax(axis=1) == marker_peak_table.index)
    add("Figure2", "B classification QC diagonal", passed, "9/9 marker sets peak in matched state")

    # 汇总每份 assignment 重现记录，保留 QC 数与精确匹配数同时等于 120766 的判据。
    reproduction_table = pd.read_csv(base / "GSE278456_step48_assignment_reproduction_source.csv")
    total_qc = int(reproduction_table["n_qc_pass"].sum())
    total_match = int(reproduction_table["n_winner_exact_match"].sum())
    add("Figure2", "state-label reproduction", total_qc == total_match == 120766, f"{total_match}/{total_qc}")


def audit_figure3() -> None:
    """核对 Figure 3 患者差值、基因、空间相关和预览布局；无返回值。"""

    # A 一行一个配对患者；均值、CI 和 FDR 都按既有结果核对。
    base = Path(__file__).resolve().parents[1] / "write/55_figure3_visual_polish/Figure3/source_data"
    patient_table = pd.read_csv(base / "Figure3A_geomx_paired_dumbbell_source.csv")
    positive_count = int((patient_table["delta"] > 0).sum())
    mean_delta = patient_table["delta"].mean()
    passed = len(patient_table) == 22 and positive_count == 16 and close(mean_delta, 0.2732193449779852)
    detail = f"n={len(patient_table)}; up={positive_count}; mean={mean_delta:.6f}"
    add("Figure3", "A paired GeoMx score", passed, detail)
    passed = (
        close(patient_table["ci_low"].iat[0], 0.0942090397046188)
        and close(patient_table["ci_high"].iat[0], 0.4522296502513515)
        and close(patient_table["fdr"].iat[0], 0.0072066426445245)
    )
    detail = (
        f"CI={patient_table['ci_low'].iat[0]:.6f},{patient_table['ci_high'].iat[0]:.6f}; "
        f"FDR={patient_table['fdr'].iat[0]:.6g}"
    )
    add("Figure3", "A GeoMx CI/FDR", passed, detail)

    # B 一行一个基因×患者；测量缺口、色阶截断数和显著基因集合均保持原判据。
    fingerprint_table = pd.read_csv(base / "Figure3B_shared_gene_patient_fingerprint_source.csv")
    measured = fingerprint_table["measured"].astype(bool)
    exceeds_color_limit = (
        fingerprint_table.loc[measured, "gene_delta"].abs() > fingerprint_table.loc[measured, "color_limit"]
    )
    clipped_count = int(exceeds_color_limit.sum())
    significant = set(fingerprint_table.loc[measured & (fingerprint_table["gene_fdr"] < 0.05), "gene"])
    passed = (
        len(fingerprint_table) == 154
        and int(measured.sum()) == 132
        and set(fingerprint_table.loc[~measured, "gene"]) == {"CCL4"}
        and clipped_count == 3
        and significant == {"CH25H", "FOLR2"}
    )
    detail = f"154 rows; measured={measured.sum()}; clipped={clipped_count}; significant={sorted(significant)}"
    add("Figure3", "B gene fingerprint", passed, detail)

    # C 一行一个空间 spot；D 一行一个样本×目标相关系数。
    spatial_table = pd.read_csv(base / "Figure3C_spatial_score_atlas_source.csv.gz")
    counts = spatial_table["sample"].value_counts().to_dict()
    passed = len(spatial_table) == 7629 and counts == {"GBM030": 4501, "GBM049": 3128}
    add("Figure3", "C spatial spots", passed, str(counts))
    rho_table = pd.read_csv(base / "Figure3D_partial_rho_effect_strip_source.csv")
    observed = {}
    for row in rho_table.itertuples():
        observed[(row.sample, row.target)] = row.partial_rho

    target = {
        ("GBM030", "MDSC"): 0.3992505547833288,
        ("GBM049", "MDSC"): 0.3547268989165864,
        ("GBM030", "MES"): 0.2799413751534758,
        ("GBM049", "MES"): 0.2704516311329281,
    }
    passed = set(observed) == set(target)
    if passed:
        for key, value in target.items():
            if not close(observed[key], value):
                passed = False
                break

    # 展示取三位小数，判断仍使用未取位的相关系数。
    rounded_correlations = {}
    for key, value in observed.items():
        rounded_correlations[key] = round(value, 3)
    add("Figure3", "D partial correlations", passed, str(rounded_correlations))

    # 预览布局是既有图件的交付检查；left/top/right/bottom 顺序不变。
    layout_table = pd.read_csv(base / "Figure3_visual_polish_preview_layout.csv")
    layout_table = layout_table.set_index("panel")
    expected_layout = {"C": (115, 2450, 2465, 4040), "D": (2610, 2450, 4205, 4235)}
    observed_layout = {}
    for panel in expected_layout:
        coordinates = []
        for column in ["left", "top", "right", "bottom"]:
            coordinates.append(int(layout_table.loc[panel, column]))
        observed_layout[panel] = tuple(coordinates)
    add("Figure3", "C-D preview layout", observed_layout == expected_layout, str(observed_layout))


def audit_figure4() -> None:
    """核对 Figure 4 配对蛋白评分、显著蛋白和全蛋白排序；无返回值。"""

    # A 一行一个患者；正向数、均值、必需字段及 sign-flip P 均按冻结结果核对。
    base = Path(__file__).resolve().parents[1] / "write/54_figure4_proteome_rebuild/Figure4/source_data"
    patient_table = pd.read_csv(base / "Figure4A_pdc_patient_delta_waterfall_source.csv")
    positive_count = int((patient_table["delta"] > 0).sum())
    mean_delta = patient_table["delta"].mean()
    passed = len(patient_table) == 105 and positive_count == 72 and close(mean_delta, 0.252534746498684)
    detail = f"n={len(patient_table)}; up={positive_count}; mean={mean_delta:.6f}"
    add("Figure4", "A paired protein score", passed, detail)
    passed = (
        {"sign_flip_p", "mean_delta", "ci_low", "ci_high", "n_recurrence_up", "n_pairs"}.issubset(patient_table.columns)
        and close(patient_table["sign_flip_p"].iat[0], 3.9999600004e-05)
    )
    add("Figure4", "A source completeness", passed, f"P={patient_table['sign_flip_p'].iat[0]:.6g}")

    # B 一行一个蛋白×患者，保留可测蛋白数、显著集合以及 PDK4 的独立样本量/FDR。
    protein_table = pd.read_csv(base / "Figure4B_raw20_protein_patient_distributions_source.csv")
    significant = set(protein_table.loc[protein_table["significant"].astype(bool), "gene"])
    pdk4_count = int(protein_table.loc[protein_table["gene"].eq("PDK4"), "n_pairs_gene"].iat[0])
    passed = (
        protein_table["gene"].nunique() == 11
        and significant == {"RHOB", "GSTM3", "FOLR2", "PDK4", "CTTNBP2"}
        and pdk4_count == 30
    )
    detail = f"11 proteins; significant={sorted(significant)}; PDK4 n={pdk4_count}"
    add("Figure4", "B measurable proteins", passed, detail)
    pdk4_fdr = float(protein_table.loc[protein_table["gene"].eq("PDK4"), "adj.P.Val"].iat[0])
    detail = f"FDR={pdk4_fdr:.15f}; rounded={pdk4_fdr:.4f}"
    add("Figure4", "B PDK4 exact FDR", close(pdk4_fdr, 0.0356470353680158), detail)

    # C 一行一个排序蛋白；总数、Miller-IM 测量数/显著数和富集统计全部保留。
    rank_table = pd.read_csv(base / "Figure4C_full_proteome_rank_landscape_source.csv")
    significant = rank_table["significant"].eq(True)
    passed = (
        len(rank_table) == 11320
        and int(rank_table["measured_miller_im"].sum()) == 11
        and int((rank_table["measured_miller_im"] & significant).sum()) == 5
    )
    detail = f"{len(rank_table)} proteins; measured={rank_table['measured_miller_im'].sum()}"
    add("Figure4", "C proteome rank", passed, detail)
    passed = (
        close(rank_table["miller_im_nes"].iat[0], 1.52795404604028)
        and close(rank_table["miller_im_gsea_p"].iat[0], 0.0490367775831874)
    )
    detail = f"NES={rank_table['miller_im_nes'].iat[0]:.6f}; nominal P={rank_table['miller_im_gsea_p'].iat[0]:.6f}"
    add("Figure4", "C source completeness", passed, detail)


def audit_figure5() -> None:
    """核对 Figure 5 暴露分组、两队列富集及调整效应；无返回值。"""

    # A 使用已交付时间轴的标签，不改动治疗分组或人数。
    base = Path(__file__).resolve().parents[1] / "write/56_figure5_pd1_rebuild/Figure5/source_data"
    timeline_table = pd.read_csv(base / "Figure5A_anti_pd1_exposure_timeline_source.csv")
    labels = set(timeline_table["lane_label"])
    expected_labels = {
        "GSE121810 neo-adjuvant (n=14)", "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)", "GSE154795 PD-1-unexposed recurrent (n=8)",
    }
    add("Figure5", "A exposure groups", labels == expected_labels, str(sorted(labels)))

    # B 每队列曲线重复携带 NES/fgsea P/camera P；按原顺序比较三项数值。
    enrichment_table = pd.read_csv(base / "Figure5B_dual_cohort_raw20_enrichment_source.csv.gz")
    expected = {
        "GSE121810": (1.89995642058625, 6.3505755649683e-05, 2.81295477233785e-05),
        "GSE154795": (1.75124003154227, 2.905445711592e-04, 0.015367499497791),
    }
    for cohort, target in expected.items():
        cohort_table = enrichment_table.loc[enrichment_table["cohort"].eq(cohort)]
        observed = (
            cohort_table["NES"].iat[0], cohort_table["fgsea_nominal_p"].iat[0], cohort_table["camera_nominal_p"].iat[0],
        )
        passed = int(cohort_table["hit"].sum()) == 19 and all(close(a, b) for a, b in zip(observed, target))
        detail = f"19/20; NES={observed[0]:.6f}; fgsea P={observed[1]:.6g}; camera P={observed[2]:.6g}"
        add("Figure5", f"B {cohort} enrichment", passed, detail)

    # C 一行一个患者；患者分组计数和每队列调整效应/FDR 分别核对。
    score_table = pd.read_csv(base / "Figure5C_patient_scores_and_adjusted_effects_source.csv")
    counts = score_table.groupby(["cohort", "analysis_group"]).size().to_dict()
    passed = len(score_table) == 50 and counts == {
        ("GSE121810", "Adj"): 15, ("GSE121810", "Neo"): 14,
        ("GSE154795", "GBM.PD1"): 13, ("GSE154795", "GBM.rec"): 8,
    }
    add("Figure5", "C patient groups", passed, str(counts))
    model_table = score_table.drop_duplicates("cohort").set_index("cohort")
    passed = (
        close(model_table.loc["GSE121810", "effect"], 0.416017028595708)
        and close(model_table.loc["GSE121810", "fdr"], 0.0098420093604063)
        and close(model_table.loc["GSE154795", "effect"], 0.427369843346401)
        and close(model_table.loc["GSE154795", "fdr"], 0.180119103326531)
    )
    add("Figure5", "C adjusted effects", passed, "GSE121810 +0.416/FDR 0.00984; GSE154795 +0.427/FDR 0.180")


def write_report() -> None:
    """写出既有核验 CSV/Markdown；失败仍写报表后退出 1，无返回值。"""

    # 汇总表一行一项核验，保留四个字段及原输出路径供发布审查读取。
    out = Path(__file__).resolve().parents[1] / "write/57_main_figure_consistency_audit"
    out.mkdir(parents=True, exist_ok=True)
    audit_table_total = pd.DataFrame(check_records)
    csv_path = out / "main_figure_consistency_checks.csv"
    audit_table_total.to_csv(csv_path, index=False)
    written_table = pd.read_csv(csv_path)
    if csv_path.stat().st_size == 0 or written_table.shape != audit_table_total.shape:
        raise AssertionError("Written main-figure audit CSV is empty or has an unexpected shape.")
    if written_table.columns.tolist() != audit_table_total.columns.tolist():
        raise AssertionError("Written main-figure audit CSV columns changed.")

    # 计数从记录现算；正文和分组汇总均使用同一张表。
    passed = int(audit_table_total["status"].eq("PASS").sum())
    failed = int(audit_table_total["status"].eq("FAIL").sum())
    lines = [
        "# Main Figure 1-5 consistency audit",
        "",
        f"- Result: {'PASS' if failed == 0 else 'FAIL'}",
        f"- Checks: {passed} passed, {failed} failed",
        "- Public naming: Miller-derived inflammatory microglial program (Miller-IM program); "
        "scores are Miller-IM scores; raw20 is retained only in internal method identifiers.",
        "- Formal delivery: independent vector PDFs; preview PDFs are review-only composites.",
        "",
        "| Figure | Passed | Failed |",
        "|---|---:|---:|",
    ]
    for section, section_table in audit_table_total.groupby("section", sort=False):
        section_passed = int(section_table["status"].eq("PASS").sum())
        section_failed = int(section_table["status"].eq("FAIL").sum())
        lines.append(f"| {section} | {section_passed} | {section_failed} |")

    # 原失败明细和退出状态保留；Markdown 回读只核对写盘完整性。
    if failed:
        lines.extend(["", "## Failed checks", ""])
        for row in audit_table_total.loc[audit_table_total["status"].eq("FAIL")].itertuples():
            lines.append(f"- {row.section} / {row.check}: {row.detail}")
    report_path = out / "MAIN_FIGURE_AUDIT.md"
    report_content = "\n".join(lines) + "\n"
    report_path.write_text(report_content, encoding="utf-8")
    if report_path.stat().st_size == 0 or report_path.read_text(encoding="utf-8") != report_content:
        raise AssertionError("Written main-figure audit report differs from memory.")

    print(f"MAIN_FIGURE_AUDIT {'PASS' if failed == 0 else 'FAIL'} passed={passed} failed={failed}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
