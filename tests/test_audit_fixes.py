"""
Regression tests of the package review of 2026-10-16 (audit/package_review_2026-10-16.md):
each asserts the behaviour wanted, where the review's probes reproduced the defect.
"""
import queue
import threading
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import ete4
import numpy as np
import pytest

from orthofinder.gene_tree import msa_trees, trim
from orthofinder.hog_orthogroups import convert, update_trees
from orthofinder.orthology import stats, tree_analyser
from orthofinder.phylo import tree
from orthofinder.species_tree import concatenated_alignment as ca
from orthofinder.species_tree import stag
from orthofinder.utils import file_handler

ROOT = Path(__file__).resolve().parents[1]


def _process_hog_tree(t, rows, names):
    incoming, outgoing, reports = queue.Queue(), queue.Queue(), queue.Queue()
    incoming.put(("OG0000000", t, {}))
    incoming.put(None)
    update_trees.process_task(incoming, outgoing, {"OG0000000": rows}, {}, names, threading.Event(), reports)
    return outgoing.get_nowait()


def test_hog_pruning_keeps_the_distances_between_genes():
    t = ete4.Tree("((0_0:1,1_0:2)n1:3,2_0:4)n0;", parser=1)
    row = {"HOG": "H0", "OG": "OG0000000", "Gene Tree Parent Clade": "n0", "A": "0_0", "C": "2_0"}
    msg = _process_hog_tree(t, [row], ["A", "C"])
    out = ete4.Tree(msg[2][0][1], parser=5)
    assert out.get_distance("0_0", "2_0") == t.get_distance("0_0", "2_0") == 8.0


def test_first_pass_resolved_trees_have_the_gene_ids(tmp_path):
    # with -X two species can have genes of the same name: the update of the
    # orthogroups reads these trees, so they have the IDs, not the names
    ta = object.__new__(tree_analyser.TreeAnalyser)
    ta.recon_trees_renamed_dir = str(tmp_path) + "/"
    ta.spec_seq_dict = {"0_0": "shared", "1_0": "shared", "2_0": "b"}
    t = tree.Tree("((0_0:1,1_0:1)n1:1,2_0:1)n0;", format=1)
    ta.write_hog_tree = True            # the first pass
    ta.write_resolved_tree(t.copy(), 3)
    assert set(tree.Tree(str(tmp_path / "OG0000003.txt"), format=1).get_leaf_names()) == {"0_0", "1_0", "2_0"}
    ta.write_hog_tree = False           # the second pass: the results, with the names
    ta.write_resolved_tree(t.copy(), 4)
    assert sorted(tree.Tree(str(tmp_path / "OG0000004.txt"), format=1).get_leaf_names()) == ["b", "shared", "shared"]


@pytest.mark.parametrize("program", ["raxml", "raxml-ng", "fasttree"])
def test_tree_program_does_not_remove_small_orthogroups(tmp_path, program):
    fn = tmp_path / "groups.txt"
    fn.write_text("0_0, 1_0, 2_0, 3_0\n0_1, 1_1\n2_1\n")
    obj = convert.HogOrthogroupSet.__new__(convert.HogOrthogroupSet)
    obj.ogs_all, obj.ogs_all_fn, obj.tree_program, obj.min_seq = None, str(fn), program, 4
    assert [len(og) for og in obj.get_ogs_all()] == [4, 2, 1]
    assert [len(og) for og in obj.all_ogs()] == [4, 2, 1]


def test_raxml_trees_only_of_four_or_more_sequences():
    class Caller:
        def get_tree_commands(self, program, alignments, outfns, ids, n_seqs, method_threads=None):
            return ids
    ogs = [[0] * 5, [0] * 4, [0] * 3, [0] * 2]
    for program, expected in (("raxml", ["OG0000000", "OG0000001"]), ("fasttree", ["OG0000000", "OG0000001", "OG0000002"])):
        gen = msa_trees.TreesForOrthogroups(Caller(), "famsa", program)
        gen.get_tree_filename = lambda i: "t%d" % i
        ids, iogs = gen.get_tree_commands(["a0", "a1", "a2", "a3"], [0, 1, 2, 3], ogs)
        assert ids == expected


def test_malformed_alignment_is_an_error_not_a_successful_exit(tmp_path):
    fn = tmp_path / "bad.fa"
    fn.write_text(">0_0\nAA\n>1_0\nA\n")
    with pytest.raises(ValueError, match="Malformed alignment"):
        trim.MSA(str(fn))
    with pytest.raises(FileNotFoundError):
        trim.main(str(tmp_path / "missing.fa"), str(tmp_path / "out.fa"))


def test_trimming_uses_the_threshold_on_the_limit():
    gaps = np.array([0, 5, 8, 9])
    f, keep = trim.get_satifactory_f(gaps, 10 - gaps, 10, .6, 2, .5)
    assert f == pytest.approx(0.3) and len(keep[0]) == 2


def test_stag_node_with_one_child():
    d = stag.get_distances_fast(tree.Tree("((0_0:1):2,1_0:3);"), 2, {"0_0": 0, "1_0": 1})
    assert d[0, 1] == d[1, 0] == 6.0


def test_species_tree_orthogroups_of_one_species_column():
    assert ca.determine_orthogroups_for_species_tree(np.array([[4.]]), [0]) == ([], 1.0)


def test_species_tree_orthogroups_are_the_ids_given():
    selected, fraction = ca.determine_orthogroups_for_species_tree(np.ones((1001, 4)), list(range(7, 1008)))
    assert (selected[0], selected[-1], fraction) == (7, 1007, 1.0)


def test_occupancy_counts_as_the_test_per_fraction():
    rng = np.random.default_rng(1)
    for _ in range(50):
        m = rng.choice([0, 1, 2, 3], (int(rng.integers(1, 200)), int(rng.integers(1, 40))),
                       p=rng.dirichlet(np.ones(4))).astype(float)
        n = m.shape[1]
        expected = [float(len(ca.single_copy_with_probability_test(f - 1e-5, m))) for f in [1.0 - i / n for i in range(n)]]
        assert ca.get_orthogroup_occupancy_info(m)[1] == expected


class _Rows(list):
    def writerow(self, row):
        self.append(list(row))


def test_size_tables_count_each_species_genes_per_orthogroup():
    rng = np.random.default_rng(4)
    i_species = [0, 2, 5, 7]
    ogs = [[[int(rng.choice(i_species + [9])), k] for k in range(int(rng.integers(2, 25)))] for _ in range(80)]
    all_genes = Counter(g[0] for og in ogs for g in og)
    all_genes.update({i: 1 for i in i_species})
    summary, per_species = _Rows(), _Rows()
    stats.stats_size_table(summary, per_species, ogs, all_genes, i_species, [set(g[0] for g in og) for og in ogs])
    # the number of orthogroups with each number of a species' genes, as counted gene by gene
    expected = [Counter(sum(1 for g in og if g[0] == i_sp) for og in ogs) for i_sp in i_species]
    bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 51, 101, 151, 201, 501, 1001, 9e99]
    table = [r[1:] for r in per_species if r and r[0] != "Number of genes per-species in orthogroup"][:len(bins) - 1]
    assert table == [[sum(n for size, n in c.items() if start <= size < end) for c in expected]
                     for start, end in zip(bins, bins[1:])]


def test_statistics_with_every_gene_unassigned(tmp_path):
    d = str(tmp_path) + "/"
    handler = SimpleNamespace(get_ogs_stats_results_directory=lambda: d,
                              get_orthogroup_results_fn_base=lambda: d + "Orthogroups",
                              get_results_seqs_dir_single_copy=lambda: d)
    with patch.object(file_handler, "handler", handler):
        stats.stats([{"0_0"}, {"1_0"}], {0: "A", 1: "B"}, [0, 1], 0, None, {})
    rows = dict(line.rstrip("\n").split("\t")[:2] for line in open(d + "Statistics_Overall.tsv") if "\t" in line)
    assert rows["Number of orthogroups"] == "0" and rows["G50 (assigned genes)"] == "NA"


def test_species_fasta_from_the_first_directory(tmp_path):
    newer, older = tmp_path / "newer", tmp_path / "older"
    for d, seq in ((newer, "NEW"), (older, "OLD")):
        d.mkdir()
        (d / "Species0.fa").write_text(">0_0\n%s\n" % seq)
    writer = msa_trees.FastaWriter([str(newer) + "/", str(older) + "/"], [0])
    assert writer.seq_lists["0_0"].strip() == "NEW"


def test_packaging_has_the_bundled_mafft_and_no_run_results():
    assert '"orthofinder.bin.mafft.libexec" = ["*"]' in (ROOT / "pyproject.toml").read_text()
    assert "prune ExampleData/OrthoFinder" in (ROOT / "MANIFEST.in").read_text()


@pytest.mark.parametrize("ogs, species_b_genes", [([], "0"), ([{"0_0", "0_1"}, {"0_2"}], "0")])
def test_statistics_with_no_genes_or_a_species_without_genes(tmp_path, ogs, species_b_genes):
    # (OrthoFinder does not take an empty FASTA file, but stats may be given such data)
    d = str(tmp_path) + "/"
    handler = SimpleNamespace(get_ogs_stats_results_directory=lambda: d,
                              get_orthogroup_results_fn_base=lambda: d + "Orthogroups",
                              get_results_seqs_dir_single_copy=lambda: d)
    with patch.object(file_handler, "handler", handler):
        stats.stats(ogs, {0: "A", 1: "B"}, [0, 1], 0, None, {})
    per_species = dict((line.split("\t")[0], line.rstrip("\n").split("\t")[1:])
                       for line in open(d + "Statistics_PerSpecies.tsv") if "\t" in line)
    assert per_species["Number of genes"][1] == species_b_genes
    assert per_species["Percentage of genes in orthogroups"][1] == "0.0"


def test_stride_species_tree_that_cannot_be_read_is_an_error(tmp_path):
    # (it printed an error and went on, then failed with UnboundLocalError)
    from orthofinder.species_tree import stride
    fn = tmp_path / "SpeciesTree.txt"
    fn.write_text("((A,B;")
    with pytest.raises(RuntimeError, match="cannot be read"):
        stride.get_root(str(fn), str(tmp_path) + "/", None, 1)


def _old_supported_distances(st, clusters):
    top = min([st.get_distance(n, topology_only=True) for n in st.traverse("preorder")
               if frozenset(n.get_leaf_names()) in clusters])
    length = min([st.get_distance(n) for n in st.traverse("preorder")
                  if frozenset(n.get_leaf_names()) in clusters and st.get_distance(n, topology_only=True) == top])
    return top, length


def test_stride_root_distances_as_get_distance():
    import random
    from orthofinder.species_tree import stride
    rng = random.Random(3)
    for case in range(100):
        nodes = [tree.TreeNode(name="S%d" % i, dist=rng.choice([rng.random(), 0.1, 0.0])) for i in range(rng.randint(3, 40))]
        while len(nodes) > 1:
            a = nodes.pop(0 if case % 4 == 0 else rng.randrange(len(nodes)))   # (some caterpillars)
            b = nodes.pop(rng.randrange(len(nodes)))
            p = tree.TreeNode(dist=rng.choice([rng.random(), 0.0]))
            p.add_child(a)
            p.add_child(b)
            nodes.append(p)
        t = nodes[0]
        clades = [frozenset(n.get_leaf_names()) for n in t.traverse()]
        clusters = set(rng.sample(clades, rng.randint(1, len(clades))))
        for r in rng.sample(clades[1:], 3):
            st = stride.root_at_clade(t, r)
            supported = stride._supported_clade_distances(st, clusters)
            top = min(d for d, _ in supported)
            assert repr((top, min(l for d, l in supported if d == top))) == repr(_old_supported_distances(st, clusters))


def test_csr_from_sorted_as_csr_matrix():
    from scipy import sparse
    from orthofinder.utils import matrices
    rng = np.random.default_rng(1)
    for _ in range(300):
        n, m = rng.integers(1, 100, 2)
        rows, cols = rng.integers(0, n, rng.integers(0, 200)), None
        cols = rng.integers(0, m, rows.size)
        if rng.random() < 0.6:     # (sorted; else given to csr_matrix)
            key = np.unique(rows * m + cols)
            rows, cols = key // m, key % m
        data = rng.random(rows.size)
        a = sparse.csr_matrix((data, (rows, cols)), shape=(n, m))
        b = matrices.csr_from_sorted(data, rows, cols, (n, m))
        for x, y in ((a.data, b.data), (a.indices, b.indices), (a.indptr, b.indptr)):
            assert x.dtype == y.dtype and np.array_equal(x, y)
