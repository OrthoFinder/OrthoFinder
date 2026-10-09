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
    assert sorted(p.name for p in (tmp_path / "Resolved_Gene_Trees" / "Tree_Data").iterdir()) == \
        ["OG0000000.js", "OG0000001.js"]
    colour = dict(zip(["Sp_A", "Sp_B"], tree_pages.PALETTE))

    data = page_data(tmp_path)
    assert [t[:2] for t in data["trees"]] == [["OG0000000", ["a1", "a2", "b1"]], ["OG0000001", ["a3", "b2", "b3"]]]
    assert data["species"] == [["Sp_A", colour["Sp_A"]], ["Sp_B", colour["Sp_B"]]]
    assert data["data"] == "Tree_Data"
    assert data["species_tree"] == {"genes": [["Sp_A", 0], ["Sp_B", 1]],
                                    "nodes": [[-1, 0, "N0"], [0, 0.5, 0], [0, 0.5, 1]]}
    # MK-LV, MKAL-, M--LV: 4 of 15 positions are gaps; columns 1 (M) and 4 (L) are identical in all;
    # branch lengths 0.1 + 0.2 + 0.3 + 0.4; the duplication with support 0.4 is not counted as supported
    assert data["trees"][0][2] == {"species": [[0, 2], [1, 1]], "length": 1.0, "duplications": [2, 1],
                                   "alignment": {"columns": 5, "gaps": 26.7, "lengths": [3, 4, 4], "identical": 40.0}}
    assert data["trees"][1][2] == {"species": [[0, 1], [1, 2]], "length": 0.6, "duplications": [0, 0]}


def test_tree_data(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    og0 = tree_data(tmp_path, "OG0000000")
    # the genes in the tree's order, each with its label and species number; the sequences
    # in the same order (not the file's); the duplications on their gene tree nodes
    assert og0 == {"og": "OG0000000",
                   "genes": [["Sp_A: a1", 0], ["Sp_A: a2", 0], ["Sp_B: b1", 1]],
                   "nodes": [[-1, 0, "n0"], [0, 0.3, "n1"], [1, 0.1, 0], [1, 0.2, 1], [0, 0.4, 2]],
                   "dups": [["n1", 1.0], ["n0", 0.4]],
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
                  .split("<script>")[1].split("</script>")[0])
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


def test_link_to_the_stats_page(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    assert page_data(tmp_path)["stats_page"] is None  # no statistics page (--stats-view)
    tree_pages.write_tree_pages(str(tmp_path), link_stats_page=True)  # as a run with --stats-view: written after this page
    assert page_data(tmp_path)["stats_page"] == "../Comparative_Genomics_Statistics/Statistics.html"
    (tmp_path / "Comparative_Genomics_Statistics").mkdir()
    (tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html").write_text("")
    tree_pages.write_tree_pages(str(tmp_path))  # the command, for results with a statistics page
    assert page_data(tmp_path)["stats_page"] == "../Comparative_Genomics_Statistics/Statistics.html"
