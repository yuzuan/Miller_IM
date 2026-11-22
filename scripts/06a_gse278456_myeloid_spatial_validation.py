#!/usr/bin/env python3
"""Purpose: localize the predefined myeloid program in GSE278456 and examine its spatial associations in GSE276841.
Inputs: myeloid and spatial H5 files, tissue positions and NIHMS2115262-supplement-S2_S5_S6.xls
    under DIR/write/30_dataset_rescue_search/source_metadata, plus the two GEO sample manifests in that module.
Outputs: GSE278456_author_marker_sets_used, single_cell_qc, program_gene_coverage, author_like_state_counts,
    patient_state_program_scores and patient_level_state_tests CSVs in the Step 30 write directory,
    plus GSE278456_cell_annotations_scores.csv.gz, GSE276841_spatial_partial_correlations.csv,
    and GSE276841_spot_scores.csv.gz.
    The matching figures directory contains GSE278456_patient_level_Mg_program_by_state and
    GSE276841_spatial_program_maps as PDF and PNG. Run: python scripts/06a_gse278456_myeloid_spatial_validation.py
    --root DIR. The released single-cell and plotting workflow uses Scanpy, AnnData and Seaborn.
"""

from __future__ import annotations

import argparse
import importlib
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

from importlib import import_module as load_figure_style

figure_style = load_figure_style("00a_functions")
render_pdf_png = figure_style.render_pdf_png

data_functions = importlib.import_module("00a_functions_data")


if __name__ == "__main__":

    # Use the original command-line options, output paths and seed with the locally downloaded inputs.
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    write = args.root / "write/30_dataset_rescue_search"
    figures = args.root / "figures/30_dataset_rescue_search"
    source = write / "source_metadata"
    write.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 1
    np.random.seed(20260712)

    # Exclude MG_FULL. Keep the first 40 positive markers at FDR < 0.05 in worksheet order.
    mg_full = [
        "PDK4", "USP53", "KLF2", "RHOB", "BHLHE41", "CTTNBP2", "SGK1", "ITM2C", "GSTM3", "CH25H", "JUN", "SIGLEC8",
        "KLF6", "FOLR2", "AC253572.2", "NLRP3",
    ]
    marker_path = source / "NIHMS2115262-supplement-S2_S5_S6.xls"
    workbook = pd.ExcelFile(marker_path)
    excluded = set(mg_full)
    markers = {}
    for sheet in workbook.sheet_names:
        marker_source_table = pd.read_excel(marker_path, sheet_name=sheet)
        if "FDR" not in marker_source_table.columns:
            marker_source_table = pd.read_excel(marker_path, sheet_name=sheet, header=1)
        gene_column = marker_source_table.columns[0]
        chosen = marker_source_table.loc[
            marker_source_table["FDR"].lt(0.05) & marker_source_table["Foldchange"].gt(0), gene_column
        ].astype(str)
        markers[sheet] = [gene for gene in chosen if gene not in excluded][:40]

    # Write the selected author markers in workbook order, with one row per state.
    marker_rows = [
        {"state": state, "n_genes": len(genes), "genes": ";".join(genes)} for state, genes in markers.items()
    ]
    author_marker_table = pd.DataFrame(marker_rows)
    author_marker_table.to_csv(write / "GSE278456_author_marker_sets_used.csv", index=False)

    # Apply QC per H5 sample using 200–8,000 genes, at least 500 counts and at most 20% mitochondrial counts.
    manifest_table = pd.read_csv(write / "GSE278456_sample_manifest.csv")
    manifest_table = manifest_table.set_index("geo_accession")
    objects = []
    qc_rows = []
    for path in sorted((source / "GSE278456_myeloid_h5").glob("*.h5")):
        accession = re.match(r"(GSM\d+)", path.name).group(1)
        sample = manifest_table.loc[accession, "sample_title"]
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        adata.obs_names = [f"{sample}:{barcode}" for barcode in adata.obs_names]
        adata.obs["sample"] = sample
        adata.obs["geo_accession"] = accession

        # Use the MT- prefix to identify mitochondrial genes and calculate QC from raw 10x counts.
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            adata.obs["n_genes_by_counts"].ge(200)
            & adata.obs["n_genes_by_counts"].le(8_000)
            & adata.obs["total_counts"].ge(500)
            & adata.obs["pct_counts_mt"].le(20)
        )
        qc_rows.append(
            {
                "sample": sample,
                "geo_accession": accession,
                "n_input": adata.n_obs,
                "n_qc_pass": int(keep.sum()),
                "fraction_qc_pass": float(keep.mean()),
            }
        )
        objects.append(adata[keep].copy())
    combined = ad.concat(objects, join="inner", merge="same")
    combined.uns["qc_summary"] = pd.DataFrame(qc_rows).to_dict("list")
    adata = combined

    # Preserve the order of program scoring and within-patient tests. Sign flips use 200,000 draws and the fixed seed.
    qc_table = pd.DataFrame(adata.uns.pop("qc_summary"))
    qc_table.to_csv(write / "GSE278456_single_cell_qc.csv", index=False)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    # Score states in workbook order. The highest score determines author_like_state.
    for state, genes in markers.items():
        adata.obs[f"author_{state}"] = data_functions.control_gene_score(adata, genes, f"tmp_author_{state}")
    score_columns = [f"author_{state}" for state in markers]
    score_matrix = adata.obs[score_columns].to_numpy()
    state_names = np.array(list(markers))
    adata.obs["author_like_state"] = state_names[np.nanargmax(score_matrix, axis=1)]
    adata.obs["Mg_full"] = data_functions.control_gene_score(adata, mg_full, "tmp_Mg_full")
    mg_leading6 = ["BHLHE41", "CH25H", "FOLR2", "JUN", "SGK1", "SIGLEC8"]
    adata.obs["Mg_leading6"] = data_functions.control_gene_score(adata, mg_leading6, "tmp_Mg_leading6")
    myeloid_genes = ["PTPRC", "TYROBP", "AIF1", "LST1", "FCER1G", "CTSS", "CSF1R"]
    adata.obs["myeloid_abundance"] = data_functions.control_gene_score(adata, myeloid_genes, "tmp_myeloid")

    # Count program genes present among the AnnData features for each sample.
    coverage = []
    for name, genes in {"Mg_full": mg_full, "Mg_leading6": mg_leading6, "myeloid_abundance": myeloid_genes}.items():
        use = [gene for gene in genes if gene in adata.var_names]
        coverage.append(
            {"program": name, "n_defined": len(genes), "n_detected": len(use), "detected_genes": ";".join(use)}
        )
    program_coverage_table = pd.DataFrame(coverage)
    program_coverage_table.to_csv(write / "GSE278456_program_gene_coverage.csv", index=False)

    # Divide each sample's state count by that sample's total number of input myeloid cells.
    state_count_table_total = (
        adata.obs.groupby(["sample", "author_like_state"], observed=True).size().rename("n_cells").reset_index()
    )
    sample_count_table = adata.obs.groupby("sample", observed=True).size().rename("n_total").reset_index()
    state_count_table_total = state_count_table_total.merge(sample_count_table, on="sample")
    state_count_table_total["fraction"] = state_count_table_total["n_cells"] / state_count_table_total["n_total"]
    state_count_table_total.to_csv(write / "GSE278456_author_like_state_counts.csv", index=False)

    # Summarize scores and cell counts within patient for the nine prespecified myeloid states.
    myeloid_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]
    state_cell_table = adata.obs.loc[
        adata.obs["author_like_state"].isin(myeloid_states), ["sample", "author_like_state", "Mg_full", "Mg_leading6"]
    ]
    patient_state_table = (
        state_cell_table.groupby(["sample", "author_like_state"], observed=True)
        .agg(n_cells=("Mg_full", "size"), Mg_full=("Mg_full", "mean"), Mg_leading6=("Mg_leading6", "mean"))
        .reset_index()
    )
    patient_state_table.to_csv(write / "GSE278456_patient_state_program_scores.csv", index=False)

    test_rows = []

    # Calculate within-patient contrasts in group order for states with at least 30 cells on each side.
    def collect_state_deltas(cell_table: pd.DataFrame, state: str, program: str) -> list[tuple[str, float]]:
        """Return pairs of sample and contrast, using the target-state mean minus the mean of other states."""
        deltas = []
        for sample, sample_table in cell_table.groupby("sample", observed=True):
            target = sample_table.loc[sample_table["author_like_state"].eq(state), program]
            other = sample_table.loc[~sample_table["author_like_state"].eq(state), program]
            if len(target) >= 30 and len(other) >= 30:
                deltas.append((sample, target.mean() - other.mean()))

        return deltas

    for program in ["Mg_full", "Mg_leading6"]:
        for state in myeloid_states:
            deltas = collect_state_deltas(state_cell_table, state, program)
            values = np.array([value for _, value in deltas])
            test_rows.append(
                {
                    "program": program,
                    "state": state,
                    "n_patients": len(values),
                    "mean_within_patient_delta": float(values.mean()) if len(values) else np.nan,
                    "median_within_patient_delta": float(np.median(values)) if len(values) else np.nan,
                    "n_positive": int((values > 0).sum()),
                    "sign_flip_p": data_functions.sign_flip_p(values),
                    "wilcoxon_p": stats.wilcoxon(values).pvalue if len(values) >= 5 and np.any(values) else np.nan,
                }
            )
    test_table = pd.DataFrame(test_rows)
    test_table["fdr"] = test_table.groupby("program")["sign_flip_p"].transform(data_functions.bh)
    test_table.to_csv(write / "GSE278456_patient_level_state_tests.csv", index=False)

    # Order states by their patient means and retain the plotting sequence and jitter settings.
    state_order = patient_state_table.groupby("author_like_state")["Mg_full"].mean().sort_values(ascending=False).index
    plt.figure(figsize=(9, 5.2))
    sns.boxplot(
        data = patient_state_table, x="author_like_state", y="Mg_full", order=state_order, color="#d9e4dc", fliersize=0
    )
    sns.stripplot(
        data  = patient_state_table,
        x     = "author_like_state",
        y     = "Mg_full",
        order = state_order,
        color = "#1f4e5f",
        size  = 3,
        alpha = 0.65,
    )
    plt.xlabel("")
    plt.ylabel("Mg-inflammatory program score")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.gcf().canvas.draw()
    pdf_path = figures / "GSE278456_patient_level_Mg_program_by_state.pdf"
    plt.savefig(pdf_path)
    render_pdf_png(pdf_path, figures / "GSE278456_patient_level_Mg_program_by_state.png", dpi=220)
    plt.close()

    # Write the public cell annotations and scores with the gzip modification time fixed at zero.
    obs_columns = ["sample", "geo_accession", "author_like_state", "Mg_full", "Mg_leading6", "myeloid_abundance"]
    adata.obs[obs_columns].to_csv(
        write / "GSE278456_cell_annotations_scores.csv.gz",
        compression = {"method": "gzip", "mtime": 0},
    )

    # Normalize and score each spatial sample independently in tissue. Adjust rank correlations for myeloid signal.
    spatial_manifest_table = pd.read_csv(write / "GSE276841_sample_manifest.csv")
    spatial_manifest_table = spatial_manifest_table.set_index("geo_accession")
    correlation_rows = []
    spot_tables = []
    map_frames = []
    for matrix_path in sorted((source / "GSE276841_spatial").glob("*_filtered_feature_bc_matrix.h5")):
        accession = re.match(r"(GSM\d+)", matrix_path.name).group(1)
        sample = spatial_manifest_table.loc[accession, "sample_title"]
        prefix = matrix_path.name.replace("_filtered_feature_bc_matrix.h5", "")
        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=10_000)
        sc.pp.log1p(adata)

        # Record one score row per spot using the original seed for each control-gene score.
        spatial_score_table_total = pd.DataFrame(index=adata.obs_names)
        spatial_score_table_total["Mg_full"] = data_functions.control_gene_score(adata, mg_full, "tmp_Mg_full")
        spatial_score_table_total["Mg_leading6"] = data_functions.control_gene_score(
            adata, mg_leading6, "tmp_Mg_leading6"
        )
        spatial_score_table_total["myeloid_abundance"] = data_functions.control_gene_score(
            adata, myeloid_genes, "tmp_myeloid"
        )
        spatial_score_table_total["E_MDSC"] = data_functions.control_gene_score(adata, markers["E-MDSC"], "tmp_E_MDSC")
        spatial_score_table_total["M_MDSC"] = data_functions.control_gene_score(adata, markers["M-MDSC"], "tmp_M_MDSC")
        t4_mes_metabolic = [
            "CHI3L1", "CD44", "SERPINE1", "VEGFA", "ADM", "CA9", "BNIP3", "NDRG1", "SLC2A1", "HK2", "LDHA", "ENO1",
            "VIM", "LGALS3", "ANXA1", "TGFBI",
        ]
        spatial_score_table_total["T4_MES_metabolic"] = data_functions.control_gene_score(
            adata, t4_mes_metabolic, "tmp_T4_MES"
        )
        spatial_score_table_total["sample"] = sample
        spatial_score_table_total["geo_accession"] = accession

        positions_path = source / "GSE276841_spatial" / f"{prefix}_tissue_positions.csv.gz"
        position_table = pd.read_csv(positions_path)
        position_table = position_table.set_index("barcode")
        spatial_score_table_total = spatial_score_table_total.join(position_table, how="left")
        spatial_score_table_total = spatial_score_table_total.loc[spatial_score_table_total["in_tissue"].eq(1)].copy()
        spot_tables.append(spatial_score_table_total.reset_index(names="barcode"))


        # Test the three original targets per sample using the shared rank-residual correlation.
        for target in ["E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            rho, p_value, n_spots = data_functions.partial_spearman(
                spatial_score_table_total["Mg_full"],
                spatial_score_table_total[target],
                spatial_score_table_total["myeloid_abundance"],
            )
            correlation_rows.append(
                {
                    "sample": sample,
                    "geo_accession": accession,
                    "target": target,
                    "n_spots": n_spots,
                    "partial_spearman_rho": rho,
                    "p_value": p_value,
                }
            )

        for program in ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]:
            map_table = spatial_score_table_total[["pxl_col_in_fullres", "pxl_row_in_fullres", program]].copy()
            map_table["program"] = program
            map_table["value"] = map_table.pop(program)
            map_table["sample"] = sample
            map_frames.append(map_table)

    spatial_test_table = pd.DataFrame(correlation_rows)
    spatial_test_table["fdr"] = data_functions.bh(spatial_test_table["p_value"])
    spatial_test_table.to_csv(write / "GSE276841_spatial_partial_correlations.csv", index=False)
    spot_score_table_total = pd.concat(spot_tables, ignore_index=True)
    spot_score_table_total.to_csv(
        write / "GSE276841_spot_scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )

    # Combine sample and program coordinates. Use the 2nd to 98th percentiles for each map's color scale.
    map_table_total = pd.concat(map_frames, ignore_index=True)
    samples = map_table_total["sample"].drop_duplicates().tolist()
    programs = ["Mg_full", "E_MDSC", "M_MDSC", "T4_MES_metabolic"]
    fig, axes = plt.subplots(len(samples), len(programs), figsize=(13, 6.4), squeeze=False)
    for row, sample in enumerate(samples):
        for column_index, program in enumerate(programs):
            ax = axes[row, column_index]
            map_table = map_table_total.loc[
                map_table_total["sample"].eq(sample) & map_table_total["program"].eq(program)
            ]
            low, high = map_table["value"].quantile([0.02, 0.98])
            points = ax.scatter(
                map_table["pxl_col_in_fullres"],
                -map_table["pxl_row_in_fullres"],
                c    = map_table["value"],
                s    = 4,
                cmap = "magma",
                vmin = low,
                vmax = high,
                linewidths = 0,
            )
            ax.set_title(f"{sample} | {program}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02)
    plt.tight_layout()
    fig.canvas.draw()
    pdf_path = figures / "GSE276841_spatial_program_maps.pdf"
    plt.savefig(pdf_path)
    render_pdf_png(pdf_path, figures / "GSE276841_spatial_program_maps.png", dpi=220)
    plt.close()

    # Check that the exported tables exist and contain data.
    for output_name in (
        "GSE278456_author_marker_sets_used.csv", "GSE278456_single_cell_qc.csv", "GSE278456_program_gene_coverage.csv",
        "GSE278456_author_like_state_counts.csv", "GSE278456_patient_state_program_scores.csv",
        "GSE278456_patient_level_state_tests.csv", "GSE278456_cell_annotations_scores.csv.gz",
        "GSE276841_spatial_partial_correlations.csv", "GSE276841_spot_scores.csv.gz",
    ):
        output_path = write / output_name
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty validation output: {output_path}")

    print(f"Single-cell samples: {qc_table['sample'].nunique()}; patient-state tests: {len(test_table)}")
    print(f"Spatial samples: {spatial_manifest_table.shape[0]}; partial correlations: {len(spatial_test_table)}")
