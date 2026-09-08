import csv
import gzip
import os
from collections import Counter
from pathlib import Path

import pytest

from tools.split_ortholog_files import split_ortholog_files

EXPECTED_RESULTS = Path(__file__).parent / "data" / "proteome_expected_results"


def prepare_compressed_results(expected_dir, results_dir):
    orthologues_dir = results_dir / "Orthologues"
    orthologues_dir.mkdir(parents=True)

    expected_rows = {}
    for expected_file in sorted(expected_dir.glob("*.tsv")):
        with expected_file.open(newline="") as infile:
            rows = list(csv.reader(infile, delimiter="\t"))

        source_species = rows[0][2]
        expected_rows[source_species] = rows

        compressed_file = orthologues_dir / f"{source_species}.tsv.gz"
        with gzip.open(compressed_file, "wt", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)

    return orthologues_dir, expected_rows


def read_tsv(path):
    with path.open(newline="") as infile:
        return list(csv.reader(infile, delimiter="\t"))


@pytest.mark.order(6)
def test_split_default_run_matches_famsa_pof_run(tmp_path):
    default_results = os.environ.get("ORTHOFINDER_TEST_RESULTS_default")
    famsa_results = os.environ.get("ORTHOFINDER_TEST_RESULTS_famsa")

    assert default_results and famsa_results, (
        "This comparison requires the default and famsa cases from "
        "test_orthofinder_runs.py to run first"
    )

    split_dir = tmp_path / "Orthologues"
    output_count, row_count, output_dir = split_ortholog_files(
        default_results, output_dir=split_dir, processes=2
    )

    expected_dir = Path(famsa_results) / "Orthologues"
    expected_files = sorted(expected_dir.glob("Orthologues_*/*.tsv"))
    actual_files = sorted(split_dir.glob("Orthologues_*/*.tsv"))
    expected_relative_paths = [
        path.relative_to(expected_dir) for path in expected_files
    ]
    actual_relative_paths = [path.relative_to(split_dir) for path in actual_files]

    assert expected_files, f"No pairwise orthologue files found in {expected_dir}"
    assert actual_relative_paths == expected_relative_paths
    assert output_count == len(expected_files)
    assert output_dir == str(split_dir)

    expected_row_count = 0
    for relative_path in expected_relative_paths:
        expected_rows = read_tsv(expected_dir / relative_path)
        actual_rows = read_tsv(split_dir / relative_path)
        assert actual_rows[0] == expected_rows[0], str(relative_path)
        assert Counter(map(tuple, actual_rows[1:])) == Counter(
            map(tuple, expected_rows[1:])
        ), str(relative_path)
        expected_row_count += len(expected_rows) - 1

    assert row_count == expected_row_count


def test_existing_outputs_are_skipped_unless_forced(tmp_path):
    expected_dir = EXPECTED_RESULTS / "Core" / "orthologues"
    results_dir = tmp_path / "Results_Core"
    prepare_compressed_results(expected_dir, results_dir)

    first_counts = split_ortholog_files(results_dir, processes=1)[:2]
    second_counts = split_ortholog_files(results_dir, processes=1)[:2]
    forced_counts = split_ortholog_files(results_dir, processes=1, force=True)[:2]

    assert first_counts[0] > 0
    assert first_counts[1] > 0
    assert second_counts == (0, 0)
    assert forced_counts == first_counts


def test_uncompressed_input_is_preferred_when_both_exist(tmp_path):
    expected_dir = EXPECTED_RESULTS / "Core" / "orthologues"
    results_dir = tmp_path / "Results_Core"
    orthologues_dir, _ = prepare_compressed_results(expected_dir, results_dir)

    compressed_file = next(orthologues_dir.glob("*.tsv.gz"))
    uncompressed_file = compressed_file.with_suffix("")
    with gzip.open(compressed_file, "rt") as infile:
        uncompressed_file.write_text(infile.read())

    output_count, _, _ = split_ortholog_files(
        results_dir, output_dir=tmp_path / "PairwiseOrthologues", processes=1
    )

    assert output_count > 0
