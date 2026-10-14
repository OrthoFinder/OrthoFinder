"""The data files of the pages to view the results (view_pages): small and quick for thousands of species."""
import json
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
