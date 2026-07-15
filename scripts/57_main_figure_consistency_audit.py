#!/usr/bin/env python3
"""Figure 1-5 正式产物的一致性与命名验收。"""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "write/57_main_figure_consistency_audit"

FIGURES = {
    "Figure1": {
        "manifest": ROOT / "write/51_figure1_final_reorganized/Figure1/Figure1_panel_manifest.csv",
        "legend": ROOT / "write/51_figure1_final_reorganized/Figure1/Figure1_legend.md",
        "panels": list("ABCDEF"),
    },
    "Figure2": {
        "manifest": ROOT / "write/50_figure2_two_panel_final/Figure2/Figure2_panel_manifest.csv",
        "legend": ROOT / "write/50_figure2_two_panel_final/Figure2/Figure2_legend.md",
        "panels": ["2A", "2B"],
    },
    "Figure3": {
        "manifest": ROOT / "write/55_figure3_visual_polish/Figure3/Figure3_panel_manifest.csv",
        "legend": ROOT / "write/55_figure3_visual_polish/Figure3/Figure3_legend.md",
        "panels": list("ABCD"),
    },
    "Figure4": {
        "manifest": ROOT / "write/54_figure4_proteome_rebuild/Figure4/Figure4_panel_manifest.csv",
        "legend": ROOT / "write/54_figure4_proteome_rebuild/Figure4/Figure4_legend.md",
        "panels": list("ABC"),
    },
    "Figure5": {
        "manifest": ROOT / "write/56_figure5_pd1_rebuild/Figure5/Figure5_panel_manifest.csv",
        "legend": ROOT / "write/56_figure5_pd1_rebuild/Figure5/Figure5_legend.md",
        "panels": list("ABC"),
    },
}

CHECKS: list[dict[str, str]] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add(section: str, check: str, passed: bool, detail: str) -> None:
    CHECKS.append(
        {
            "section": section,
            "check": check,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        }
    )


def close(value: float, target: float, atol: float = 1e-9) -> bool:
    return bool(np.isclose(float(value), float(target), atol=atol, rtol=0))


def pdf_text(path: Path) -> str:
    result = subprocess.run(
        ["pdftotext", str(path), "-"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def audit_delivery() -> None:
    forbidden = re.compile(r"raw20|Miller\s+raw|Mg-inflammatory|Mg-IM|pharmacodynamic", re.IGNORECASE)
    for figure, spec in FIGURES.items():
        manifest = pd.read_csv(spec["manifest"])
        observed_panels = manifest["panel"].astype(str).tolist()
        add(figure, "panel mapping", observed_panels == spec["panels"], str(observed_panels))

        for row in manifest.to_dict("records"):
            panel = str(row["panel"])
            sources = [Path(item) for item in str(row["source"]).split(";")]
            pdf = Path(row["pdf"])
            png = Path(row["png"])
            add(figure, f"{panel} files exist", all(path.exists() for path in [*sources, pdf, png]), "source/PDF/PNG")
            if "source_sha256" in row and pd.notna(row["source_sha256"]):
                add(figure, f"{panel} source hash", len(sources) == 1 and sha256(sources[0]) == row["source_sha256"], sources[0].name)
                add(figure, f"{panel} PDF hash", sha256(pdf) == row["pdf_sha256"], pdf.name)
                add(figure, f"{panel} PNG hash", sha256(png) == row["png_sha256"], png.name)

            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0, 0))
                passed_dpi = len(dpi) >= 2 and all(abs(float(value) - 600) < 1 for value in dpi[:2])
            add(figure, f"{panel} PNG 600 dpi", passed_dpi, str(dpi))

            fonts = subprocess.run(["pdffonts", str(pdf)], check=True, capture_output=True, text=True).stdout
            add(figure, f"{panel} PDF no Type3", "Type 3" not in fonts, pdf.name)

            text = pdf_text(pdf)
            match = forbidden.search(text)
            add(figure, f"{panel} visible nomenclature", match is None, "no legacy public name" if match is None else match.group(0))

        legend = Path(spec["legend"]).read_text(encoding="utf-8")
        match = forbidden.search(legend)
        add(figure, "legend nomenclature", match is None, "Miller-IM public terminology" if match is None else match.group(0))
        misleading_fdr = re.search(r"targeted\s+(?:fgsea\s+)?FDR", legend, re.IGNORECASE)
        add(figure, "legend program-level P labeling", misleading_fdr is None, "nominal P used for single-set enrichment" if misleading_fdr is None else misleading_fdr.group(0))
        if "message" in manifest.columns:
            messages = "\n".join(manifest["message"].fillna("").astype(str))
            match = forbidden.search(messages)
            add(figure, "manifest message nomenclature", match is None, "Miller-IM public terminology" if match is None else match.group(0))

    figure1_legend = Path(FIGURES["Figure1"]["legend"]).read_text(encoding="utf-8")
    definition = "20-gene Miller-derived inflammatory microglial program (Miller-IM program)"
    add("Cross-figure", "first formal definition", definition in figure1_legend, definition)


def audit_figure1() -> None:
    base = ROOT / "write/51_figure1_final_reorganized/Figure1/source_data"
    design = pd.read_csv(base / "Figure1A_paired_cohort_design_source.csv").set_index("dataset")
    add("Figure1", "A cohort inputs", design["input_libraries_or_matrices"].to_dict() == {"GSE174554": 91, "GSE274546": 111}, str(design["input_libraries_or_matrices"].to_dict()))
    add("Figure1", "A clean myeloid cells", design["clean_myeloid"].to_dict() == {"GSE174554": 23344, "GSE274546": 68109}, str(design["clean_myeloid"].to_dict()))
    add("Figure1", "A paired patients", design["formal_pairs_threshold20"].to_dict() == {"GSE174554": 18, "GSE274546": 45}, str(design["formal_pairs_threshold20"].to_dict()))

    expected_gsea = {
        "GSE174554": (2.3157685724298, 4.13267973361141e-7, 18),
        "GSE274546": (1.79932425427, 0.00397813985660199, 45),
    }
    for panel, cohort in [("B", "GSE174554"), ("C", "GSE274546")]:
        table = pd.read_csv(base / f"Figure1{panel}_{cohort}_raw20_gsea_source.csv.gz")
        observed = (table["NES"].iat[0], table["nominal_P"].iat[0], int(table["n_pairs"].iat[0]))
        target = expected_gsea[cohort]
        passed = close(observed[0], target[0]) and close(observed[1], target[1]) and observed[2] == target[2]
        add("Figure1", f"{panel} GSEA statistics", passed, f"NES={observed[0]:.6f}; nominal P={observed[1]:.6g}; n={observed[2]}")

    rank = pd.read_csv(base / "Figure1D_shared_raw20_outliers_source.csv")
    counts = (len(rank), int(rank["miller_raw20"].sum()), int(rank["shared_leading_edge"].sum()), int(rank["display_label"].fillna("").ne("").sum()))
    add("Figure1", "D all-gene/shared-gene counts", counts == (9432, 19, 8, 5), str(counts))
    effects = pd.read_csv(base / "Figure1E_shared_leading_edge_gene_effects_source.csv")
    add("Figure1", "E gene-by-cohort tests", len(effects) == 16 and effects["gene"].nunique() == 8 and int((effects["FDR"] < 0.05).sum()) == 0, f"rows={len(effects)}; genes={effects['gene'].nunique()}; FDR<0.05={(effects['FDR'] < 0.05).sum()}")
    lopo = pd.read_csv(base / "Figure1F_leave_one_patient_out_summary.csv").set_index("dataset")
    passed = (
        int(lopo.loc["GSE174554", "positive_NES_folds"]) == 18
        and int(lopo.loc["GSE274546", "positive_NES_folds"]) == 45
        and int(lopo.loc["GSE174554", "nominal_P_lt_0_05_folds"]) == 18
        and int(lopo.loc["GSE274546", "nominal_P_lt_0_05_folds"]) == 45
        and close(lopo.loc["GSE274546", "max_nominal_P"], 0.0311, atol=5e-5)
    )
    add("Figure1", "F LOPO stability", passed, "positive 18/18, 45/45; nominal P<0.05 18/18, 45/45; max 0.0311")


def audit_figure2() -> None:
    base = ROOT / "write/50_figure2_two_panel_final/Figure2/source_data"
    waterfall = pd.read_csv(base / "Figure2_state_patient_waterfall_source.csv")
    expected = {
        "MCG1": (18, 17), "MCG2": (20, 16), "MCG3": (5, 1),
        "MCG4": (12, 1), "MCG5": (10, 3), "MAC1": (3, 0),
        "MAC2": (12, 5), "E-MDSC": (17, 0), "M-MDSC": (11, 0),
    }
    observed = {
        state: (int(group["n_patients"].iat[0]), int(group["n_positive"].iat[0]))
        for state, group in waterfall.groupby("state", sort=False)
    }
    add("Figure2", "A state patient counts", len(waterfall) == 108 and observed == expected, str(observed))
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
    observed_statistics = {
        state: (float(group["mean_delta"].iat[0]), float(group["fdr_9_states"].iat[0]))
        for state, group in waterfall.groupby("state", sort=False)
    }
    passed_statistics = set(observed_statistics) == set(expected_statistics) and all(
        close(observed_statistics[state][0], target[0]) and close(observed_statistics[state][1], target[1])
        for state, target in expected_statistics.items()
    )
    add(
        "Figure2",
        "A state means and FDRs",
        passed_statistics,
        str({state: (round(value[0], 6), round(value[1], 6)) for state, value in observed_statistics.items()}),
    )
    discordant = set(waterfall.loc[~waterfall["score_pseudobulk_mean_direction_concordant"].astype(bool), "state"])
    add("Figure2", "A score-pseudobulk boundary", discordant == {"MCG3"}, str(discordant))

    genes = pd.read_csv(base / "Figure2_raw20_marker_double_heatmap_raw20_source.csv")
    marker = pd.read_csv(base / "Figure2_raw20_marker_double_heatmap_marker_qc_source.csv")
    add("Figure2", "B held-out gene block", len(genes) == 171 and genes["gene"].nunique() == 19 and genes["state"].nunique() == 9, f"{len(genes)} rows; {genes['gene'].nunique()} genes; {genes['state'].nunique()} states")
    pivot = marker.pivot(index="tested_marker_set", columns="assigned_state", values="within_marker_set_z")
    add("Figure2", "B classification QC diagonal", len(marker) == 81 and all(pivot.idxmax(axis=1) == pivot.index), "9/9 marker sets peak in matched state")
    reproduction = pd.read_csv(base / "GSE278456_step48_assignment_reproduction_source.csv")
    total_qc = int(reproduction["n_qc_pass"].sum())
    total_match = int(reproduction["n_winner_exact_match"].sum())
    add("Figure2", "state-label reproduction", total_qc == total_match == 120766, f"{total_match}/{total_qc}")


def audit_figure3() -> None:
    base = ROOT / "write/55_figure3_visual_polish/Figure3/source_data"
    patients = pd.read_csv(base / "Figure3A_geomx_paired_dumbbell_source.csv")
    passed = len(patients) == 22 and int((patients["delta"] > 0).sum()) == 16 and close(patients["delta"].mean(), 0.2732193449779852)
    add("Figure3", "A paired GeoMx score", passed, f"n={len(patients)}; up={(patients['delta'] > 0).sum()}; mean={patients['delta'].mean():.6f}")
    add("Figure3", "A GeoMx CI/FDR", close(patients["ci_low"].iat[0], 0.0942090397046188) and close(patients["ci_high"].iat[0], 0.4522296502513515) and close(patients["fdr"].iat[0], 0.0072066426445245), f"CI={patients['ci_low'].iat[0]:.6f},{patients['ci_high'].iat[0]:.6f}; FDR={patients['fdr'].iat[0]:.6g}")

    fingerprint = pd.read_csv(base / "Figure3B_shared_gene_patient_fingerprint_source.csv")
    measured = fingerprint["measured"].astype(bool)
    clipped = int((fingerprint.loc[measured, "gene_delta"].abs() > fingerprint.loc[measured, "color_limit"]).sum())
    significant = set(fingerprint.loc[measured & (fingerprint["gene_fdr"] < 0.05), "gene"])
    add("Figure3", "B gene fingerprint", len(fingerprint) == 176 and int(measured.sum()) == 154 and set(fingerprint.loc[~measured, "gene"]) == {"CCL4"} and clipped == 4 and significant == {"CH25H", "FOLR2"}, f"176 rows; measured={measured.sum()}; clipped={clipped}; significant={sorted(significant)}")

    spatial = pd.read_csv(base / "Figure3C_spatial_score_atlas_source.csv.gz")
    counts = spatial["sample"].value_counts().to_dict()
    add("Figure3", "C spatial spots", len(spatial) == 7629 and counts == {"GBM030": 4501, "GBM049": 3128}, str(counts))
    rho = pd.read_csv(base / "Figure3D_partial_rho_effect_strip_source.csv")
    observed = {(row.sample, row.target): row.partial_rho for row in rho.itertuples()}
    target = {
        ("GBM030", "MDSC"): 0.3992505547833288,
        ("GBM049", "MDSC"): 0.3547268989165864,
        ("GBM030", "MES"): 0.2799413751534758,
        ("GBM049", "MES"): 0.2704516311329281,
    }
    add("Figure3", "D partial correlations", set(observed) == set(target) and all(close(observed[key], value) for key, value in target.items()), str({key: round(value, 3) for key, value in observed.items()}))
    layout = pd.read_csv(base / "Figure3_visual_polish_preview_layout.csv").set_index("panel")
    expected_layout = {
        "C": (115, 2450, 2465, 4040),
        "D": (2610, 2450, 4205, 4235),
    }
    observed_layout = {
        panel: tuple(int(layout.loc[panel, column]) for column in ["left", "top", "right", "bottom"])
        for panel in expected_layout
    }
    add("Figure3", "C-D preview layout", observed_layout == expected_layout, str(observed_layout))


def audit_figure4() -> None:
    base = ROOT / "write/54_figure4_proteome_rebuild/Figure4/source_data"
    patients = pd.read_csv(base / "Figure4A_pdc_patient_delta_waterfall_source.csv")
    passed = len(patients) == 105 and int((patients["delta"] > 0).sum()) == 72 and close(patients["delta"].mean(), 0.252534746498684)
    add("Figure4", "A paired protein score", passed, f"n={len(patients)}; up={(patients['delta'] > 0).sum()}; mean={patients['delta'].mean():.6f}")
    add("Figure4", "A source completeness", {"sign_flip_p", "mean_delta", "ci_low", "ci_high", "n_recurrence_up", "n_pairs"}.issubset(patients.columns) and close(patients["sign_flip_p"].iat[0], 3.9999600004e-05), f"P={patients['sign_flip_p'].iat[0]:.6g}")

    proteins = pd.read_csv(base / "Figure4B_raw20_protein_patient_distributions_source.csv")
    significant = set(proteins.loc[proteins["significant"].astype(bool), "gene"])
    expected = {"RHOB", "GSTM3", "FOLR2", "PDK4", "CTTNBP2"}
    pdk4_n = int(proteins.loc[proteins["gene"].eq("PDK4"), "n_pairs_gene"].iat[0])
    add("Figure4", "B measurable proteins", proteins["gene"].nunique() == 11 and significant == expected and pdk4_n == 30, f"11 proteins; significant={sorted(significant)}; PDK4 n={pdk4_n}")
    pdk4_fdr = float(proteins.loc[proteins["gene"].eq("PDK4"), "adj.P.Val"].iat[0])
    add("Figure4", "B PDK4 exact FDR", close(pdk4_fdr, 0.0356470353680158), f"FDR={pdk4_fdr:.15f}; rounded={pdk4_fdr:.4f}")

    rank = pd.read_csv(base / "Figure4C_full_proteome_rank_landscape_source.csv")
    significant = rank["significant"].eq(True)
    passed = len(rank) == 11320 and int(rank["measured_miller_im"].sum()) == 11 and int((rank["measured_miller_im"] & significant).sum()) == 5
    add("Figure4", "C proteome rank", passed, f"{len(rank)} proteins; measured={rank['measured_miller_im'].sum()}")
    add("Figure4", "C source completeness", close(rank["miller_im_nes"].iat[0], 1.52795404604028) and close(rank["miller_im_gsea_p"].iat[0], 0.0490367775831874), f"NES={rank['miller_im_nes'].iat[0]:.6f}; nominal P={rank['miller_im_gsea_p'].iat[0]:.6f}")


def audit_figure5() -> None:
    base = ROOT / "write/56_figure5_pd1_rebuild/Figure5/source_data"
    timeline = pd.read_csv(base / "Figure5A_anti_pd1_exposure_timeline_source.csv")
    labels = set(timeline["lane_label"])
    expected_labels = {
        "GSE121810 neo-adjuvant (n=14)", "GSE121810 adjuvant-only (n=15)",
        "GSE154795 PD-1 exposed (n=13)", "GSE154795 PD-1-unexposed recurrent (n=8)",
    }
    add("Figure5", "A exposure groups", labels == expected_labels, str(sorted(labels)))

    enrichment = pd.read_csv(base / "Figure5B_dual_cohort_raw20_enrichment_source.csv.gz")
    expected = {
        "GSE121810": (1.89995642058625, 6.3505755649683e-05, 2.81295477233785e-05),
        "GSE154795": (1.75124003154227, 2.905445711592e-04, 0.015367499497791),
    }
    for cohort, target in expected.items():
        table = enrichment.loc[enrichment["cohort"].eq(cohort)]
        observed = (table["NES"].iat[0], table["fgsea_nominal_p"].iat[0], table["camera_nominal_p"].iat[0])
        passed = int(table["hit"].sum()) == 19 and all(close(a, b) for a, b in zip(observed, target))
        add("Figure5", f"B {cohort} enrichment", passed, f"19/20; NES={observed[0]:.6f}; fgsea P={observed[1]:.6g}; camera P={observed[2]:.6g}")

    scores = pd.read_csv(base / "Figure5C_patient_scores_and_adjusted_effects_source.csv")
    counts = scores.groupby(["cohort", "analysis_group"]).size().to_dict()
    add("Figure5", "C patient groups", len(scores) == 50 and counts == {("GSE121810", "Adj"): 15, ("GSE121810", "Neo"): 14, ("GSE154795", "GBM.PD1"): 13, ("GSE154795", "GBM.rec"): 8}, str(counts))
    model = scores.drop_duplicates("cohort").set_index("cohort")
    passed = (
        close(model.loc["GSE121810", "effect"], 0.416017028595708)
        and close(model.loc["GSE121810", "fdr"], 0.0098420093604063)
        and close(model.loc["GSE154795", "effect"], 0.427369843346401)
        and close(model.loc["GSE154795", "fdr"], 0.180119103326531)
    )
    add("Figure5", "C adjusted effects", passed, "GSE121810 +0.416/FDR 0.00984; GSE154795 +0.427/FDR 0.180")


def write_report() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(CHECKS)
    csv_path = OUT / "main_figure_consistency_checks.csv"
    frame.to_csv(csv_path, index=False)
    passed = int(frame["status"].eq("PASS").sum())
    failed = int(frame["status"].eq("FAIL").sum())
    lines = [
        "# Main Figure 1-5 consistency audit",
        "",
        f"- Result: {'PASS' if failed == 0 else 'FAIL'}",
        f"- Checks: {passed} passed, {failed} failed",
        "- Public naming: Miller-derived inflammatory microglial program (Miller-IM program); scores are Miller-IM scores; raw20 is retained only in internal method identifiers.",
        "- Formal delivery: independent vector PDFs; preview PDFs are review-only composites.",
        "",
        "| Figure | Passed | Failed |",
        "|---|---:|---:|",
    ]
    for section, group in frame.groupby("section", sort=False):
        lines.append(f"| {section} | {int(group['status'].eq('PASS').sum())} | {int(group['status'].eq('FAIL').sum())} |")
    if failed:
        lines.extend(["", "## Failed checks", ""])
        for row in frame.loc[frame["status"].eq("FAIL")].itertuples():
            lines.append(f"- {row.section} / {row.check}: {row.detail}")
    (OUT / "MAIN_FIGURE_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"MAIN_FIGURE_AUDIT {'PASS' if failed == 0 else 'FAIL'} passed={passed} failed={failed}")
    if failed:
        raise SystemExit(1)


def main() -> None:
    audit_delivery()
    audit_figure1()
    audit_figure2()
    audit_figure3()
    audit_figure4()
    audit_figure5()
    write_report()


if __name__ == "__main__":
    main()
