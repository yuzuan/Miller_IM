#!/usr/bin/env python3
"""Sum the frozen state assignments from the original H5AD CSR integer counts using AnnData.

Inputs: provenance/step1_input_manifest.csv (cohorts and pan-myeloid references).
A locally prepared step1_input_manifest.csv can supply the same reference fields.
Inputs under NOA_LEGACY_ROOT: write/71_recurrence_state_decomposition/input_manifest.csv,
author_marker_provenance.csv, and <dataset>_cell_scores.csv.gz. The manifest identifies the formal patient list,
raw20 definition, and H5AD files. Use NOA_PUBLIC_DATA_ROOT for H5AD and pan-myeloid reference counts and metadata.
Outputs in step3_state_analysis: pseudobulk/<dataset>_<state>_{counts.csv.gz,metadata.csv},
<dataset>_cell_assignment_audit.csv.gz, aggregation_checks.csv, frozen_marker_provenance.csv, input_audit.csv,
patient_state_cell_counts.csv, expression_pair_coverage.csv, pseudobulk_manifest.csv, aggregation_parameters.json.
Run: python NOA_R1/scripts/03a_aggregate_state_counts.py
"""

from __future__ import annotations

import gc
import json
from importlib import import_module
from importlib.metadata import version
from itertools import product
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse


def main():
    """Aggregate and export raw state counts and their coverage. Returns nothing."""
    revision = Path(__file__).resolve().parents[1]

    # Use the frozen cohort manifest and the paths of the original state assignment inputs.
    manifest_path = revision / "step1_input_manifest.csv"
    if not manifest_path.is_file():
        manifest_path = revision / "provenance/step1_input_manifest.csv"
    formal_table = pd.read_csv(manifest_path)
    functions = import_module("00a_functions")
    project = functions.legacy_root(revision)
    state_source_dir = project / "write/71_recurrence_state_decomposition"
    input_table = pd.read_csv(state_source_dir / "input_manifest.csv")
    input_paths = input_table.set_index("role").path

    def write_input_path(reference, root):
        """Rebase a frozen write-tree reference on the current root and return its path."""
        _, separator, suffix = reference.rpartition("/write/")
        relative_path = Path("write") / suffix if separator else Path(reference)
        path = root / relative_path
        return path

    # The fold audit defines the formal patients. The raw20 list checks marker independence.
    folds_table = pd.read_csv(write_input_path(input_paths["formal_patients"], project))
    markers_table = pd.read_csv(state_source_dir / "author_marker_provenance.csv")
    definition_table = pd.read_csv(write_input_path(input_paths["signature"], project))
    ordered_table = definition_table.sort_values("raw20_order").drop_duplicates("gene")
    raw20 = set(ordered_table.gene)

    output_dir = revision / "step3_state_analysis"
    pseudobulk_dir = output_dir / "pseudobulk"
    pseudobulk_dir.mkdir(parents=True, exist_ok=True)
    fine_states = ["MCG1", "MCG2", "MCG3", "MCG4", "MCG5", "MAC1", "MAC2", "M-MDSC", "E-MDSC"]
    broad_states = ["MCG-like", "MAC_MDSC-like"]
    check_rows = []

    def check(dataset, label, passed, detail=""):
        """Record a frozen input or alignment check and stop on failure. Returns nothing."""
        check_rows.append(dict(dataset=dataset, check=label, passed=bool(passed), detail=detail))
        pd.DataFrame(check_rows).to_csv(output_dir / "aggregation_checks.csv", index=False)
        if not passed:
            raise ValueError(f"{dataset}: {label}: {detail}")

    for state in fine_states:
        markers = markers_table.loc[markers_table.state == state, "markers_used"].iloc[0].split(";")
        check("Both", f"{state}_markers_exclude_target20", not set(markers) & raw20)
    markers_table.to_csv(output_dir / "frozen_marker_provenance.csv", index=False)

    # Aggregate each cohort in the original manifest order using the frozen per-cell labels.
    count_rows = []
    coverage_rows = []
    manifest_rows = []
    input_rows = []
    public_root = functions.public_data_root(project)
    for cohort in formal_table.itertuples():
        dataset = cohort.dataset
        cohort_folds = (folds_table.dataset == dataset) & (folds_table.run_type == "Leave-one-patient-out")
        patients = sorted(folds_table.loc[cohort_folds, "omitted_patient"].unique())
        check(dataset, "formal_patient_number", len(patients) == cohort.expected_pairs)
        labels_table = pd.read_csv(
            state_source_dir / f"{dataset}_cell_scores.csv.gz",
            usecols=["cell_id", "dataset", "patient_id", "condition", "state", "assignment_score", "assignment_margin"]
        )
        valid_labels = (
            labels_table.cell_id.is_unique
            and set(labels_table.patient_id) == set(patients)
            and set(labels_table.state) <= set(fine_states)
            and set(labels_table.dataset) == {dataset}
        )
        check(dataset, "frozen_labels_unique_and_formal", valid_labels)

        # Align labels to H5AD row order so that counts and cell-level assignments share the same index.
        raw_path = write_input_path(input_paths[f"raw_counts_{dataset}"], public_root)
        with h5py.File(raw_path, "r") as handle:
            obs_table = read_elem(handle["obs"])
            var_table = read_elem(handle["var"])
            check(dataset, "unique_h5ad_ids", obs_table.index.is_unique and var_table.index.is_unique)
            check(dataset, "cell_id_coverage", set(labels_table.cell_id) <= set(obs_table.index))
            keep = obs_table.index.isin(labels_table.cell_id)
            selected_table = obs_table.loc[keep]
            labels_table = labels_table.set_index("cell_id").loc[selected_table.index].reset_index(names="cell_id")
            aligned = (
                np.array_equal(selected_table.patient_id.astype(str), labels_table.patient_id)
                and np.array_equal(selected_table.condition.astype(str), labels_table.condition)
            )
            check(dataset, "cell_patient_time_alignment", aligned)

            # The frozen manifest states each formal endpoint total and the size of the original gene universe.
            primary_cells = int((labels_table.condition == "Primary").sum())
            recurrent_cells = int((labels_table.condition == "Recurrent").sum())
            correct_totals = (
                int(keep.sum()) == cohort.expected_primary_cells + cohort.expected_recurrent_cells
                and primary_cells == cohort.expected_primary_cells
                and recurrent_cells == cohort.expected_recurrent_cells
            )
            check(dataset, "formal_cells", correct_totals)
            check(dataset, "formal_gene_count", len(var_table) == cohort.expected_input_genes)
            raw_counts = sparse_dataset(handle["layers/counts"])
            encoding = str(handle["layers/counts"].attrs["encoding-type"])
            integer_counts = np.issubdtype(handle["layers/counts/data"].dtype, np.integer)
            check(dataset, "csr_raw_integer_layer", encoding == "csr_matrix" and integer_counts)

            # One row per observed patient, endpoint and fine state, with genes in H5AD column order.
            unit_index = pd.MultiIndex.from_frame(labels_table[["patient_id", "condition", "state"]])
            codes, unique_units = pd.factorize(unit_index, sort=True)
            units_table = unique_units.to_frame(index=False)
            units_table.columns = ["patient_id", "condition", "state"]
            units_table["n_cells"] = np.bincount(codes, minlength=len(units_table))
            pseudobulk = np.zeros((len(units_table), len(var_table)), dtype=np.int64)
            row_codes = np.full(len(obs_table), -1, dtype=np.int64)
            row_codes[keep] = codes

            # Use the same 2,048-row CSR blocks, visiting only blocks with formal cells.
            nonnegative = True
            starts = np.unique(np.flatnonzero(keep) // 2048) * 2048
            for start in starts:
                stop = min(start + 2048, len(obs_table))
                block_codes = row_codes[start:stop]
                valid = block_codes >= 0
                block = raw_counts[start:stop, :][valid].tocsr().astype(np.int64)
                nonnegative = nonnegative and bool((block.data >= 0).all())
                groups = block_codes[valid]
                active_groups = np.unique(groups)
                group_index = np.searchsorted(active_groups, groups)
                weights = np.ones(len(groups), dtype=np.int64)
                incidence = sparse.csr_matrix(
                    (weights, (group_index, np.arange(len(groups)))),
                    shape=(len(active_groups), len(groups))
                )
                pseudobulk[active_groups] += (incidence @ block).toarray()

            check(dataset, "all_raw_counts_nonnegative", nonnegative)
            genes = var_table.index.astype(str)

        # Check totals across all fine states against the pan-myeloid matrix, gene by gene.
        valid_counts = bool((pseudobulk >= 0).all()) and int(units_table.n_cells.sum()) == len(labels_table)
        check(dataset, "aggregated_counts_valid", valid_counts)
        units_table["library_size"] = pseudobulk.sum(axis=1)
        check(dataset, "observed_units_positive_libraries", bool((units_table.library_size > 0).all()))
        reference_metadata_path = public_root / cohort.metadata_path.removeprefix("<public_data_root>/")
        reference_metadata_table = pd.read_csv(reference_metadata_path)
        reference_metadata_table = reference_metadata_table[reference_metadata_table.patient_id.isin(patients)].copy()
        reference_counts_path = public_root / cohort.counts_path.removeprefix("<public_data_root>/")
        reference_counts_table = pd.read_csv(reference_counts_path, index_col=0)
        reference_counts_table = reference_counts_table.loc[genes, reference_metadata_table.sample_id]
        valid_reference = (
            reference_counts_table.index.is_unique
            and reference_counts_table.columns.is_unique
            and len(reference_counts_table) == len(genes)
            and len(reference_metadata_table) == 2 * cohort.expected_pairs
        )
        check(dataset, "reference_gene_sample_ids", valid_reference)

        max_error = 0
        for reference in reference_metadata_table.itertuples():
            patient_match = units_table.patient_id == reference.patient_id
            endpoint_match = units_table.condition == reference.condition
            unit_rows = patient_match & endpoint_match
            total_counts = pseudobulk[unit_rows.to_numpy()].sum(axis=0)
            reference_counts = reference_counts_table[reference.sample_id].to_numpy(dtype=np.int64)
            max_error = max(max_error, int(np.max(np.abs(total_counts - reference_counts))))
        check(
            dataset, "whole_genome_pan_myeloid_reproduction", max_error == 0,
            f"maximum integer difference={max_error}"
        )
        del reference_counts_table

        # Pool broad groups in reference sample order, retaining only observed patient-endpoint rows.
        broad_rows = {state: [] for state in broad_states}
        broad_counts = {state: [] for state in broad_states}
        members = {"MCG-like": fine_states[:5], "MAC_MDSC-like": fine_states[5:]}
        for state, reference in product(broad_states, reference_metadata_table.itertuples()):
            state_match = units_table.state.isin(members[state])
            patient_match = units_table.patient_id == reference.patient_id
            endpoint_match = units_table.condition == reference.condition
            select = state_match & patient_match & endpoint_match
            cell_count = int(units_table.loc[select, "n_cells"].sum())
            if cell_count == 0:
                continue

            values = pseudobulk[select.to_numpy()].sum(axis=0)
            metadata_row = dict(
                patient_id   = reference.patient_id,
                condition    = reference.condition,
                n_cells      = cell_count,
                library_size = int(values.sum())
            )
            broad_rows[state].append(metadata_row)
            broad_counts[state].append(values)

        # Fine groups use observed state rows. Broad groups use the pooled MCG1–5 or MAC/MDSC counts.
        for state in fine_states + broad_states:
            if state in fine_states:
                level = "fine"
                select = units_table.state == state
                state_columns = ["patient_id", "condition", "n_cells", "library_size"]
                metadata_table = units_table.loc[select, state_columns].copy().reset_index(drop=True)
                counts_matrix = pseudobulk[select.to_numpy()]
            else:
                level = "broad"
                metadata_table = pd.DataFrame(broad_rows[state])
                counts_matrix = np.asarray(broad_counts[state], dtype=np.int64)

            # Export gene-by-sample counts and metadata in the same sample order.
            sample_ids = [f"{dataset}_{state}_U{i:03d}" for i in range(1, len(metadata_table) + 1)]
            metadata_table.insert(0, "sample_id", sample_ids)
            metadata_table.insert(1, "dataset", dataset)
            metadata_table.insert(2, "state", state)
            aligned = len(metadata_table) == len(counts_matrix)
            aligned = aligned and bool(np.array_equal(metadata_table.library_size, counts_matrix.sum(axis=1)))
            check(dataset, f"{state}_counts_metadata_alignment", aligned)
            counts_path = pseudobulk_dir / f"{dataset}_{state}_counts.csv.gz"
            metadata_path = pseudobulk_dir / f"{dataset}_{state}_metadata.csv"
            counts_table = pd.DataFrame(counts_matrix.T, index=genes, columns=metadata_table.sample_id)
            counts_table.rename_axis("gene").to_csv(counts_path, compression={"method": "gzip", "mtime": 0})
            metadata_table.to_csv(metadata_path, index=False)
            if min(counts_path.stat().st_size, metadata_path.stat().st_size) <= 0:
                raise ValueError(f"Empty pseudobulk output for {dataset} {state}")

            manifest_row = dict(
                dataset       = dataset,
                state         = state,
                level         = level,
                counts_path   = counts_path.relative_to(revision).as_posix(),
                metadata_path = metadata_path.relative_to(revision).as_posix()
            )
            manifest_rows.append(manifest_row)

            # Include observed zeros only in the abundance table. Expression still requires captured cells.
            for patient, condition in product(patients, ("Primary", "Recurrent")):
                select = (metadata_table.patient_id == patient) & (metadata_table.condition == condition)
                endpoint_counts = metadata_table.loc[select, "n_cells"]
                cell_count = 0 if endpoint_counts.empty else int(endpoint_counts.iloc[0])
                patient_match = reference_metadata_table.patient_id == patient
                endpoint_match = reference_metadata_table.condition == condition
                denominator = int(reference_metadata_table.loc[patient_match & endpoint_match, "n_cells"].iloc[0])
                count_row = dict(
                    dataset            = dataset,
                    patient_id         = patient,
                    condition          = condition,
                    state              = state,
                    level              = level,
                    n_cells            = cell_count,
                    total_myeloid_cells = denominator,
                    fraction           = cell_count / denominator
                )
                count_rows.append(count_row)

            # Both endpoints must meet the 20-cell primary or 10-cell sensitivity threshold.
            pairs_table = metadata_table.pivot(index="patient_id", columns="condition", values="n_cells")
            pairs_table = pairs_table.reindex(patients).fillna(0)
            for threshold in (20, 10):
                eligible = (pairs_table.Primary >= threshold) & (pairs_table.Recurrent >= threshold)
                available = pairs_table.index[eligible].tolist()
                coverage_row = dict(
                    dataset                  = dataset,
                    state                    = state,
                    level                    = level,
                    threshold                = threshold,
                    n_pairs                  = len(available),
                    formal_pairs             = len(patients),
                    eligible_patients        = ";".join(available),
                    model_eligible           = len(available) >= 6,
                    included_primary_cells   = int(pairs_table.loc[available, "Primary"].sum()),
                    included_recurrent_cells = int(pairs_table.loc[available, "Recurrent"].sum())
                )
                coverage_rows.append(coverage_row)

        # Retain assignment margin and score diagnostics without refitting or relabeling cells.
        labels_table["broad_state"] = np.where(labels_table.state.str.startswith("MCG"), "MCG-like", "MAC_MDSC-like")
        labels_table.to_csv(
            output_dir / f"{dataset}_cell_assignment_audit.csv.gz",
            index = False,
            compression = {"method": "gzip", "mtime": 0}
        )
        input_row = dict(
            dataset                = dataset,
            n_cells                = len(labels_table),
            n_pairs                = len(patients),
            n_genes                = len(genes),
            n_observed_state_units = len(units_table),
            low_margin_fraction    = float((labels_table.assignment_margin < .05).mean()),
            nonpositive_winner      = int((labels_table.assignment_score <= 0).sum()),
            pan_count_max_error    = max_error,
            raw_h5ad               = "<public_data_root>/" + raw_path.relative_to(public_root).as_posix()
        )
        input_rows.append(input_row)
        print(f"{dataset}: {len(patients)} pairs, {len(labels_table)} cells, maximum count difference {max_error}")
        del pseudobulk, labels_table
        gc.collect()

    # The state fractions sum to one within each patient, endpoint and classification level.
    counts_table_total = pd.DataFrame(count_rows)
    fraction_sums = counts_table_total.groupby(["dataset", "patient_id", "condition", "level"]).fraction.sum()
    check("Both", "relative_fractions_sum_one", bool(np.allclose(fraction_sums, 1, atol=1e-12)))
    counts_table_total.to_csv(output_dir / "patient_state_cell_counts.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(output_dir / "expression_pair_coverage.csv", index=False)
    pd.DataFrame(manifest_rows).to_csv(output_dir / "pseudobulk_manifest.csv", index=False)
    pd.DataFrame(input_rows).to_csv(output_dir / "input_audit.csv", index=False)

    # Export the frozen aggregation settings alongside the versions used to read and sum the counts.
    packages = {package: version(package) for package in ["anndata", "h5py", "numpy", "scipy", "pandas"]}
    parameters = dict(
        source     = "frozen Step71 assignments; raw counts layer",
        chunk_size = 2048,
        sum_dtype  = "int64",
        thresholds = [20, 10],
        min_pairs  = 6,
        packages   = packages
    )
    parameter_text = json.dumps(parameters, ensure_ascii=False, indent=2)
    (output_dir / "aggregation_parameters.json").write_text(parameter_text, encoding="utf-8")
    print(f"Aggregation checks passed: {len(check_rows)}", flush=True)


if __name__ == "__main__":
    main()
