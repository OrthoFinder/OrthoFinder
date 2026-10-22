import csv
import gzip
from collections import Counter
from pathlib import Path

import pytest

from tools.split_ortholog_files import split_ortholog_files

EXPECTED_RESULTS = Path(__file__).parent / "data" / "proteome_expected_results"


def prepare_compressed_results(expected_dir, results_dir):
    orthologs_dir = results_dir / "Orthologs"
    orthologs_dir.mkdir(parents=True)

    expected_rows = {}
    for expected_file in sorted(expected_dir.glob("*.tsv")):
        with expected_file.open(newline="") as infile:
            rows = list(csv.reader(infile, delimiter="\t"))

        source_species = rows[0][2]
        expected_rows[source_species] = rows

        compressed_file = orthologs_dir / f"{source_species}.tsv.gz"
        with gzip.open(compressed_file, "wt", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)

    return orthologs_dir, expected_rows


def read_tsv(path):
    with path.open(newline="") as infile:
        return list(csv.reader(infile, delimiter="\t"))


@pytest.mark.order(6)
def test_split_default_run_matches_famsa_pof_run(tmp_path, orthofinder_results):
    # famsa is the default run with -pof: its pairwise files are those to expect
    default_results = orthofinder_results("default")
    famsa_results = orthofinder_results("famsa")

    split_dir = tmp_path / "Orthologs"
    output_count, row_count, output_dir = split_ortholog_files(
        default_results, output_dir=split_dir, processes=2
    )

    expected_dir = Path(famsa_results) / "Orthologs"
    expected_files = sorted(expected_dir.glob("Orthologs_*/*.tsv"))
    actual_files = sorted(split_dir.glob("Orthologs_*/*.tsv"))
    expected_relative_paths = [
        path.relative_to(expected_dir) for path in expected_files
    ]
    actual_relative_paths = [path.relative_to(split_dir) for path in actual_files]

    assert expected_files, f"No pairwise ortholog files found in {expected_dir}"
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
    expected_dir = EXPECTED_RESULTS / "Core" / "orthologs"
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
    expected_dir = EXPECTED_RESULTS / "Core" / "orthologs"
    results_dir = tmp_path / "Results_Core"
    orthologs_dir, _ = prepare_compressed_results(expected_dir, results_dir)

    compressed_file = next(orthologs_dir.glob("*.tsv.gz"))
    uncompressed_file = compressed_file.with_suffix("")
    with gzip.open(compressed_file, "rt") as infile:
        uncompressed_file.write_text(infile.read())

    output_count, _, _ = split_ortholog_files(
        results_dir, output_dir=tmp_path / "PairwiseOrthologs", processes=1
    )

    assert output_count > 0


def _write_species_files(d, n_species):
    for i in range(n_species):
        rows = ["Orthogroup\tSpecies\tsp%d\tOrthologs\n" % i]
        for j in range(n_species):
            if j != i:
                rows.append("OG%07d\tsp%d\tg%d_%d\th%d_%d\n" % (j, j, i, j, j, i))
        (d / ("sp%d.tsv" % i)).write_text("".join(rows))


def _read_tree(d):
    return {
        str(p.relative_to(d)): p.read_text()
        for p in sorted(d.glob("Orthologs_*/*.tsv"))
    }


def test_package_split_whatever_the_buffer(tmp_path):
    # -pof: each species' file read once, its rows held in buffers and appended in
    # pieces (ortholog_writers.TextBuffers): the same files with a buffer of a few rows
    # (written again and again) as with one that holds them all
    from orthofinder.orthology.ortholog_writers import split_ortholog_files as split_pof

    dirs = {name: tmp_path / name for name in ("default", "tiny", "one_row", "parallel")}
    for d in dirs.values():
        d.mkdir()
        _write_species_files(d, 7)
    split_pof(str(dirs["default"]))
    split_pof(str(dirs["tiny"]), buffer_bytes=30)
    split_pof(str(dirs["one_row"]), buffer_bytes=1)
    split_pof(str(dirs["parallel"]), buffer_bytes=30, n_processes=3)   # (species split at the same time)

    files = _read_tree(dirs["default"])
    assert len(files) == 7 * 6
    assert files == _read_tree(dirs["tiny"]) == _read_tree(dirs["one_row"]) == _read_tree(dirs["parallel"])
    assert files["Orthologs_sp2/sp2__v__sp5.tsv"] == (
        "Orthogroup\tsp2\tsp5\nOG0000005\tg2_5\th5_2\n"
    )


@pytest.mark.parametrize("in_memory", [True, False])
def test_tool_split_in_memory_and_in_batches(tmp_path, monkeypatch, in_memory):
    import tools.split_ortholog_files as tool

    d = tmp_path / "Orthologs"
    d.mkdir()
    _write_species_files(d, 7)
    if not in_memory:
        monkeypatch.setattr(tool, "IN_MEMORY_BYTES_PER_PROCESS", 0)
        monkeypatch.setattr(tool, "max_open_per_process", lambda n_processes: 2)
    output_count, row_count, _ = tool.split_ortholog_files(str(d), processes=1)

    assert (output_count, row_count) == (7 * 6, 7 * 6)
    assert (d / "Orthologs_sp2" / "sp2__v__sp5.tsv").read_text() == (
        "Orthogroup\tsp2\tsp5\nOG0000005\tg2_5\th5_2\n"
    )
    assert (d / "Orthologs_sp2" / "sp2__v__sp0.tsv").read_text() == (
        "Orthogroup\tsp2\tsp0\nOG0000000\tg2_0\th0_2\n"
    )


def test_tool_split_selected_species(tmp_path):
    import tools.split_ortholog_files as tool

    d = tmp_path / "Orthologs"
    d.mkdir()
    _write_species_files(d, 5)
    output_count, _, _ = tool.split_ortholog_files(
        str(d), processes=1, species=["sp1", "sp3"], partners=["sp0", "sp3"]
    )

    assert output_count == 3
    assert sorted(str(p.relative_to(d)) for p in d.glob("Orthologs_*/*.tsv")) == [
        "Orthologs_sp1/sp1__v__sp0.tsv",
        "Orthologs_sp1/sp1__v__sp3.tsv",
        "Orthologs_sp3/sp3__v__sp0.tsv",
    ]
    with pytest.raises(ValueError, match="--species: no ortholog file for spX"):
        tool.split_ortholog_files(str(d), processes=1, species=["spX"])


def test_package_split_many_rows_per_file_and_compressed(tmp_path):
    # rows of a pairwise file in the order of the species' file, whatever the buffer;
    # and the compressed files (gzip members) the same text
    from orthofinder.orthology.ortholog_writers import split_ortholog_files as split_pof
    import random
    rng = random.Random(4)
    for name, kwargs in (("a", {}), ("b", {"buffer_bytes": 200}), ("gz", {"buffer_bytes": 200, "q_compress": True})):
        d = tmp_path / name
        d.mkdir()
        for i in range(4):
            rows = ["Orthogroup\tSpecies\tsp%d\tOrthologs\n" % i]
            rnd = random.Random(i)
            for k in range(300):
                j = rnd.choice([x for x in range(4) if x != i])
                rows.append("OG%07d\tsp%d\tg%d_%d\th%d_%d\n" % (k, j, i, k, j, k))
            (d / ("sp%d.tsv" % i)).write_text("".join(rows))
        split_pof(str(d), **kwargs)
    a, b = _read_tree(tmp_path / "a"), _read_tree(tmp_path / "b")
    assert len(a) == 12 and a == b
    gz = {str(p.relative_to(tmp_path / "gz"))[:-3]: gzip.open(p, "rt").read()
          for p in sorted((tmp_path / "gz").glob("Orthologs_*/*.tsv.gz"))}
    assert gz == a
