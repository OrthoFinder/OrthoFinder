"""The data files of the pages to view the results (view_pages): small and quick for thousands of species."""
import json
import os
import re
import shutil
import subprocess

import pytest

from orthofinder.orthology import stats_page
from orthofinder.utils import view_pages


def data_file(path):
    m = re.match(r'ofData\("([^"]+)",(.*)\);\n$', path.read_text(), re.S)
    assert m and m.group(1) == path.stem
    return json.loads(m.group(2))


def test_sparse_cells():
    # the non-empty cells and their columns, without splitting at every tab
    assert list(view_pages.sparse_cells("OG1\t\tb1, b2\t\t\tc1\n")) == [(0, "OG1"), (2, "b1, b2"), (5, "c1")]
    assert list(view_pages.sparse_cells("\tx\t\r\n")) == [(1, "x")]
    assert list(view_pages.sparse_cells("\t\t\n")) == []


def test_compact_ids():
    ids = [f"OG{i:07d}" for i in range(12)]
    assert view_pages.compact_ids(ids) == {"prefix": "OG", "width": 7, "n": 12}
    assert view_pages.compact_ids(["N3.HOG0000000", "N3.HOG0000001"]) == {"prefix": "N3.HOG", "width": 7, "n": 2}
    # not in order, not from 0, not all alike: the list
    assert view_pages.compact_ids(["OG0000000", "OG0000002"]) == ["OG0000000", "OG0000002"]
    assert view_pages.compact_ids(["OG0000001", "OG0000002"]) == ["OG0000001", "OG0000002"]
    assert view_pages.compact_ids(["OG0000000", "OG000001"]) == ["OG0000000", "OG000001"]
    assert view_pages.compact_ids([]) == []


def test_write_chunks(tmp_path):
    rows = [["a"] * 5, ["b"] * 5, ["c"] * 30, ["d"], ["e"]]
    starts = view_pages.write_chunks(str(tmp_path), "T", rows, len, limit=12)
    # rows of about limit in weight each (weight + 1 a row); a row heavier than limit on its own
    assert starts == [0, 2, 3]
    assert [data_file(tmp_path / f"T_{c}.js") for c in range(3)] == [rows[:2], rows[2:3], rows[3:]]


def test_write_gene_index(tmp_path):
    genes = [("b2", 1), ("A1", 0), ("a2", 0), ("c1", 2), ("B1", 1)]
    index = view_pages.write_gene_index(str(tmp_path), genes, shard=2)
    # sorted ignoring case, in files of shard genes; the first gene of each, in lower case
    assert index == ["a1", "b1", "c1"]
    assert [data_file(tmp_path / f"Genes_{i}.js") for i in range(3)] == [
        [["A1", 0], ["a2", 0]], [["B1", 1], ["b2", 1]], [["c1", 2]]]


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to run the page's script")
def test_find_genes(tmp_path):
    # the page's search: the genes whose name starts with it, across the files of the index, with the
    # same name in two species (in two files) found in both
    genes = [("g1", 0), ("g10", 1), ("g11", 2), ("g2", 3), ("h1", 4), ("x", 5), ("x", 6), ("x", 7), ("y", 8)]
    index = view_pages.write_gene_index(str(tmp_path), genes, shard=2)
    script = tmp_path / "find.js"
    script.write_text(f"""
const fs = require("fs"), window = globalThis;
const D = {{data: {json.dumps(str(tmp_path))}}};
{view_pages.DATA_JS}
loadData = name => {{ eval(fs.readFileSync(`${{D.data}}/${{name}}.js`, "utf8")); return Promise.resolve(dataCache[name]); }};
(async () => {{
  const index = {json.dumps(index)}, out = {{}};
  for (const q of ["g1", "G", "g2", "h", "x", "y", "z", "a"]) out[q] = (await findGenes(index, q)).map(e => e[1]);
  out.limit = (await findGenes(index, "g", 2)).map(e => e[1]);
  console.log(JSON.stringify(out));
}})();
""")
    result = subprocess.run(["node", str(script)], capture_output=True, text=True, check=True)
    assert json.loads(result.stdout) == {"g1": [0, 1, 2], "G": [0, 1, 2, 3], "g2": [3], "h": [4], "x": [5, 6, 7],
                                         "y": [8], "z": [], "a": [], "limit": [0, 1]}


def test_gene_counts_sparse(tmp_path):
    # only the cells other than 0 (most, for many species), not the total
    fn = tmp_path / "Orthogroups.GeneCount.tsv"
    fn.write_text("Orthogroup\tSp_A\tSp_B\tSp_C\tTotal\nOG0000000\t0\t12\t1\t13\nOG0000001\t10\t0\t0\t10\n")
    assert stats_page.read_gene_counts(str(fn)) == (["Sp_A", "Sp_B", "Sp_C"],
                                                    [["OG0000000", [1, 12, 2, 1]], ["OG0000001", [0, 10]]])


def test_shared_orthogroups():
    counts = [["OG0", [0, 2, 1, 1]], ["OG1", [0, 1, 1, 1]], ["OG2", [0, 3]], ["OG3", [2, 1]]]
    assert stats_page.shared_orthogroups(counts, 3) == [[3, 2, 0], [2, 2, 0], [0, 0, 1]]


def test_chunks_from_an_iterator_as_from_a_list(tmp_path):
    import random
    rng = random.Random(3)
    rows = [["x"] * rng.randint(0, 30) for _ in range(200)]
    (tmp_path / "a").mkdir(), (tmp_path / "b").mkdir()
    starts_a = view_pages.write_chunks(str(tmp_path / "a"), "T", rows, len, limit=50)
    starts_b = view_pages.write_chunks(str(tmp_path / "b"), "T", iter(rows), len, limit=50)
    assert starts_a == starts_b
    assert sorted(os.listdir(tmp_path / "a")) == sorted(os.listdir(tmp_path / "b"))
    for fn in os.listdir(tmp_path / "a"):
        assert (tmp_path / "a" / fn).read_text() == (tmp_path / "b" / fn).read_text()


def test_gene_index_sorted_in_runs_as_in_memory(tmp_path):
    # merged runs give what sorted() gives, equal names (ignoring case) in the order given
    import random
    rng = random.Random(8)
    genes = [(rng.choice(["g", "G", "a", "Ab", "ab", "z"]) + str(rng.randint(0, 30)), i) for i in range(5000)]
    (tmp_path / "mem").mkdir(), (tmp_path / "runs").mkdir()
    index_mem = view_pages.write_gene_index(str(tmp_path / "mem"), genes, shard=300)
    index_runs = view_pages.write_gene_index(str(tmp_path / "runs"), iter(genes), shard=300, run=137)
    assert index_mem == index_runs
    assert os.listdir(tmp_path / "runs") and not [f for f in os.listdir(tmp_path / "runs") if not f.endswith(".js")]
    for fn in os.listdir(tmp_path / "mem"):
        assert (tmp_path / "mem" / fn).read_text() == (tmp_path / "runs" / fn).read_text()


def test_command_threads_whatever_the_environment(monkeypatch):
    import subprocess
    from orthofinder.utils import command_runner, parallel_task_manager
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        monkeypatch.setitem(parallel_task_manager.my_env, var, "32")
    seen = {}

    def run(command, env, stdout=None, stderr=None):
        seen.update({k: env[k] for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")})
        return 0, b"", b""
    monkeypatch.setattr(parallel_task_manager, "run_monitored_command", run)
    command_runner.run_command("prog -t METHODTHREAD", "1")
    assert seen == {"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    command_runner.run_command("prog -t METHODTHREAD", "2")
    assert seen["OMP_NUM_THREADS"] == "2" and seen["MKL_NUM_THREADS"] == "1"
