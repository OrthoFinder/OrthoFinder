"""--tree-view: a page that draws the trees (in the browser) and each gene tree's alignment."""
import json
import re
import shutil
import subprocess

import pytest
from ete4 import Tree

from orthofinder.phylo import tree_pages


def test_leaf_species():
    species = sorted(["Sp_A", "Sp_A_long", "Sp_B"], key=len, reverse=True)
    genes = {"g7": "Sp_B"}
    assert tree_pages.leaf_species("Sp_A_long_g1", species, genes) == ("Sp_A_long", "g1")
    assert tree_pages.leaf_species("Sp_A_g2", species, genes) == ("Sp_A", "g2")
    assert tree_pages.leaf_species("g7", species, genes) == ("Sp_B", "g7")  # -X: no species in the names
    assert tree_pages.leaf_species("unknown", species, genes) == (None, "unknown")


def test_flat_tree():
    nodes, leaves = tree_pages.flat_tree(Tree("((a:0.1,b:0.2)n1:0.3,c:0.4)n0;", parser=1))
    # preorder, a parent before its children: [parent, branch length, leaf number or node name]
    assert nodes == [[-1, 0, "n0"], [0, 0.3, "n1"], [1, 0.1, 0], [1, 0.2, 1], [0, 0.4, 2]]
    assert leaves == ["a", "b", "c"]


def test_flat_tree_deep():
    newick = "(" * 3000 + "a" + ",b)" * 3000 + ";"  # deeper than Python's recursion limit
    nodes, leaves = tree_pages.flat_tree(Tree(newick))
    assert len(leaves) == 3001 and len(nodes) == 6001
    json.dumps(nodes)  # a flat list: no recursion to write it


def make_results(d):
    (d / "Orthogroups").mkdir(parents=True)
    (d / "Orthogroups" / "Orthogroups.tsv").write_text(
        "Orthogroup\tSp_A\tSp_B\nOG0000000\ta1, a2\tb1\nOG0000001\ta3\tb2, b3\n")
    (d / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text("Orthogroup\tSp_A\tSp_B\nOG0000002\ta4\t\n")
    (d / "Species_Tree").mkdir()
    (d / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt").write_text("(Sp_A:0.5,Sp_B:0.5)N0;")
    (d / "Resolved_Gene_Trees").mkdir()
    (d / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text(
        "OG0000000: ((Sp_A_a1:0.1,Sp_A_a2:0.2)n1:0.3,Sp_B_b1:0.4)n0;\n"
        "OG0000001: (Sp_A_a3:0.1,(Sp_B_b2:0.2,Sp_B_b3:0.2)n1:0.1)n0;\n"
        "OG0000009: ((not a tree;\n")
    (d / "Gene_Duplication_Events").mkdir()
    (d / "Gene_Duplication_Events" / "Duplications.tsv").write_text(
        "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n"
        "OG0000000\tSp_A\tn1\t1.0\tTerminal\tSp_A_a1\tSp_A_a2\n"
        "OG0000000\tN0\tn0\t0.4\tNon-Terminal\tSp_A_a1, Sp_A_a2\tSp_B_b1\n")
    (d / "Phylogenetic_Hierarchical_Orthogroups").mkdir()
    (d / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv").write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\n"
        "N0.HOG0000000\tOG0000000\tn0\ta1, a2\tb1\n"
        "N0.HOG0000001\tOG0000001\tn0\ta3\tb2, b3\n")
    (d / "MultipleSequenceAlignments").mkdir()  # OG0000001 has none
    (d / "MultipleSequenceAlignments" / "OG0000000.fa").write_text(
        ">Sp_B_b1\nMK-L\nV\n>Sp_A_a2\nMKAL-\n>Sp_A_a1\nM--LV\n")


def page_data(d):
    page = (d / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").read_text()
    return json.loads(re.search(r"const D = (.*);\n", page).group(1))


def tree_data(d, og):
    text = (d / "Resolved_Gene_Trees" / "Tree_Data" / f"{og}.js").read_text()
    assert text.startswith("treeLoaded(") and text.endswith(");\n")
    return json.loads(text[len("treeLoaded("):-len(");\n")])


def test_write_tree_pages(tmp_path, capsys):
    make_results(tmp_path)
    assert tree_pages.write_tree_pages(str(tmp_path)) == 2
    assert "1 gene trees could not be read" in capsys.readouterr().out
    # each gene tree; the explorer's HOGs of each species tree node and genes of each species
    assert sorted(p.name for p in (tmp_path / "Resolved_Gene_Trees" / "Tree_Data").iterdir()) == [
        "All_0.js", "All_1.js", "All_2.js", "All_Gained.js", "Genes_0.js", "HOGs_N0.js", "HOGs_N0_0.js",
        "OG0000000.js", "OG0000001.js", "Species_0.js", "Species_1.js"]
    colour = dict(zip(["Sp_A", "Sp_B"], tree_pages.PALETTE))

    data = page_data(tmp_path)
    # the trees in columns: their genes, species, length, duplications, genes per species (for at most 20
    # species) and alignment ([columns, % gaps, shortest, median and longest sequence, % identical columns])
    # MK-LV, MKAL-, M--LV: 4 of 15 positions are gaps; columns 1 (M) and 4 (L) are identical in all;
    # branch lengths 0.1 + 0.2 + 0.3 + 0.4; the duplication with support 0.4 is not counted as supported
    assert data["trees"] == {"ids": {"prefix": "OG", "width": 7, "n": 2}, "genes": [3, 3], "n_species": [2, 2],
                             "length": [1.0, 0.6], "dups": [2, 0], "dups05": [1, 0],
                             "species": [[[0, 2], [1, 1]], [[0, 1], [1, 2]]], "alignment": [[5, 26.7, 3, 4, 4, 40.0], 0]}
    # the genes of the trees, sorted, each with its tree's number: to find a gene by the start of its name
    assert data["genes"] == ["a1"]
    assert data_file(tmp_path, "Genes_0") == [["a1", 0], ["a2", 0], ["a3", 1], ["b1", 0], ["b2", 1], ["b3", 1]]
    assert data["species"] == [["Sp_A", colour["Sp_A"]], ["Sp_B", colour["Sp_B"]]]
    assert data["data"] == "Tree_Data"
    assert data["species_tree"] == {"genes": [["Sp_A", 0], ["Sp_B", 1]],
                                    "nodes": [[-1, 0, "N0"], [0, 0.5, 0], [0, 0.5, 1]]}


def test_tree_data(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    og0 = tree_data(tmp_path, "OG0000000")
    # the genes in the tree's order, each with its label and species number; the sequences
    # in the same order (not the file's); the duplications on their gene tree nodes
    assert og0 == {"og": "OG0000000",
                   "genes": [["Sp_A: a1", 0], ["Sp_A: a2", 0], ["Sp_B: b1", 1]],
                   "nodes": [[-1, 0, "n0"], [0, 0.3, "n1"], [1, 0.1, 0], [1, 0.2, 1], [0, 0.4, 2]],
                   "dups": [["n1", 1.0, "Sp_A"], ["n0", 0.4, "N0"]],
                   "seqs": ["M--LV", "MKAL-", "MK-LV"]}
    og1 = tree_data(tmp_path, "OG0000001")
    assert "seqs" not in og1 and og1["dups"] == []


def test_alignment_sequences_not_in_the_tree(tmp_path):
    make_results(tmp_path)
    (tmp_path / "MultipleSequenceAlignments" / "OG0000000.fa").write_text(">Sp_A_a1\nMK\n>Sp_B_extra\nMV\n")
    tree_pages.write_tree_pages(str(tmp_path))
    og0 = tree_data(tmp_path, "OG0000000")
    # a gene of the tree without a sequence has None; a sequence not in the tree comes after the tree's genes
    assert og0["genes"] == [["Sp_A: a1", 0], ["Sp_A: a2", 0], ["Sp_B: b1", 1], ["Sp_B: extra", 1]]
    assert og0["seqs"] == ["MK", None, None, "MV"]


def test_read_alignment(tmp_path):
    fn = tmp_path / "a.fa"
    fn.write_text(">x y\nAC-\nGT\n\n>z\nAC\n")
    assert tree_pages.read_alignment(str(fn)) == {"x y": "AC-GT", "z": "AC"}


def test_no_gene_trees_nothing_written(tmp_path):
    (tmp_path / "Orthogroups").mkdir()  # e.g. -og: no gene trees
    assert tree_pages.write_tree_pages(str(tmp_path)) == 0
    assert [p.name for p in tmp_path.iterdir()] == ["Orthogroups"]


def test_empty_gene_trees_file_nothing_written(tmp_path):
    make_results(tmp_path)
    (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("")  # as after -og
    assert tree_pages.write_tree_pages(str(tmp_path)) == 0
    assert sorted(p.name for p in (tmp_path / "Resolved_Gene_Trees").iterdir()) == ["Resolved_Gene_Trees.txt"]


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the page's script")
def test_page_script_is_valid_javascript(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    js = tmp_path / "page.js"
    js.write_text((tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").read_text()
                  .split("<script>")[-1].split("</script>")[0])
    result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_needs_no_pyqt6(tmp_path, monkeypatch):
    make_results(tmp_path)
    monkeypatch.setitem(__import__("sys").modules, "PyQt6", None)  # as if PyQt6 were not installed
    monkeypatch.setitem(__import__("sys").modules, "ete4.treeview", None)
    assert tree_pages.write_tree_pages(str(tmp_path)) == 2


def test_command(tmp_path, capsys):
    make_results(tmp_path)
    assert tree_pages.main([str(tmp_path)]) == 0
    assert "Gene tree page (2 trees)" in capsys.readouterr().out
    assert (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").exists()


def test_links_to_the_other_pages(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    data = page_data(tmp_path)
    assert data["stats_page"] is None and data["groups_page"] is None  # no other pages
    assert data["nav"] == [["Gene trees", "../Resolved_Gene_Trees/Resolved_Gene_Trees.html", 1, "trees"]]
    # as a run with --view: the other pages are written after this one
    tree_pages.write_tree_pages(str(tmp_path), {"trees", "stats", "groups"})
    data = page_data(tmp_path)
    assert data["stats_page"] == "../Comparative_Genomics_Statistics/Statistics.html"
    assert data["groups_page"] == "../Orthogroups/Groups.html"
    assert [name for name, _, _, _ in data["nav"]] == ["Gene trees", "Statistics", "Orthogroups & HOGs"]
    # the command, for results with a statistics page
    (tmp_path / "Comparative_Genomics_Statistics").mkdir()
    (tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html").write_text("")
    tree_pages.write_tree_pages(str(tmp_path))
    data = page_data(tmp_path)
    assert data["stats_page"] == "../Comparative_Genomics_Statistics/Statistics.html" and data["groups_page"] is None


def data_file(d, name):
    """The value of a data file (view_pages.write_data): ofData(name, value);"""
    path = d / "Resolved_Gene_Trees" / "Tree_Data" / f"{name}.js"
    m = re.match(r'ofData\("([^"]+)",(.*)\);\n$', path.read_text(), re.S)
    assert m and m.group(1) == name
    return json.loads(m.group(2))


def test_explorer_data(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    data = page_data(tmp_path)
    assert data["hog_nodes"] == [["N0", 2]]
    assert data["species_genes"] == [4, 3]  # assigned and unassigned
    # the HOGs of a node in columns, and each HOG's genes by species number in chunks
    assert data_file(tmp_path, "HOGs_N0") == {"ids": {"prefix": "N0.HOG", "width": 7, "n": 2},
                                              "ogs": ["OG0000000", "OG0000001"], "clades": ["n0", "n0"],
                                              "genes": [3, 3], "species": [2, 2], "chunks": [0]}
    assert data_file(tmp_path, "HOGs_N0_0") == [[[0, ["a1", "a2"]], [1, ["b1"]]], [[0, ["a3"]], [1, ["b2", "b3"]]]]
    assert data_file(tmp_path, "Species_0") == {
        "orthogroups": [["OG0000000", ["a1", "a2"]], ["OG0000001", ["a3"]]], "unassigned": ["a4"]}
    assert data_file(tmp_path, "Species_1") == {
        "orthogroups": [["OG0000000", ["b1"]], ["OG0000001", ["b2", "b3"]]], "unassigned": []}


def test_explorer_without_hogs(tmp_path):
    make_results(tmp_path)
    shutil.rmtree(tmp_path / "Phylogenetic_Hierarchical_Orthogroups")  # e.g. -og
    tree_pages.write_tree_pages(str(tmp_path))
    assert page_data(tmp_path)["hog_nodes"] == []


def test_reconcile_duplication_and_losses():
    # A duplication at the root: one copy goes on to A and B (lost in C), the other to A and C (lost in B)
    s_nodes, s_names = tree_pages.flat_tree(Tree("((A:1,B:1)N1:1,C:1)N0;", parser=1))
    number = {"A": 0, "B": 1, "C": 2}
    smap = tree_pages.SpeciesTreeMap(s_nodes, [number[n] for n in s_names])
    g_nodes, g_names = tree_pages.flat_tree(Tree("((A_a1,B_b1)n1,(A_a2,C_c1)n2)n0;", parser=1))
    copies, losses, root = tree_pages.reconcile(g_nodes, [number[n[0]] for n in g_names], smap)
    node = {(s_names[x] if isinstance(x, int) else x): k for k, (_, _, x) in enumerate(s_nodes)}
    assert root == node["N0"]
    # only the nodes with any
    assert {name: copies[k] for name, k in node.items() if k in copies} == {"N0": 1, "N1": 2, "A": 2, "B": 1, "C": 1}
    assert {name: losses[k] for name, k in node.items() if k in losses} == {"B": 1, "C": 1}


def test_all_orthogroups(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    # species tree nodes: N0, Sp_A, Sp_B; both gene trees are a speciation at N0 (duplications within species)
    assert page_data(tmp_path)["all_orthogroups"] == {
        "n_trees": 2, "n_orthogroups": 2, "copies": [2, 2, 2], "losses": [0, 0, 0], "gains": [2, 0, 0],
        "dups": [[1, 0], [1, 1], [0, 0]]}  # from Duplications.tsv: n0 at N0 (support 0.4), n1 at Sp_A
    # the orthogroups gained at each node, and those with losses and duplications at each
    assert data_file(tmp_path, "All_Gained") == [["OG0000000", "OG0000001"], [], []]
    assert data_file(tmp_path, "All_0") == {"lost": [], "duplicated": [["OG0000000", 1]]}
    assert data_file(tmp_path, "All_1") == {"lost": [], "duplicated": [["OG0000000", 1]]}
    assert data_file(tmp_path, "All_2") == {"lost": [], "duplicated": []}
