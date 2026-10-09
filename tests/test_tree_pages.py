"""--tree-view: pictures of the trees (ete4, which needs PyQt6) and a page for them."""
import json
import re
import shutil
import subprocess

import pytest

from orthofinder.phylo import tree_pages


def test_leaf_species():
    species = sorted(["Sp_A", "Sp_A_long", "Sp_B"], key=len, reverse=True)
    genes = {"g7": "Sp_B"}
    assert tree_pages.leaf_species("Sp_A_long_g1", species, genes) == ("Sp_A_long", "g1")
    assert tree_pages.leaf_species("Sp_A_g2", species, genes) == ("Sp_A", "g2")
    assert tree_pages.leaf_species("g7", species, genes) == ("Sp_B", "g7")  # -X: no species in the names
    assert tree_pages.leaf_species("unknown", species, genes) == (None, "unknown")


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


def treeview_available():
    try:
        tree_pages.load_treeview()
        return True
    except tree_pages.TreeViewUnavailable:
        return False


needs_treeview = pytest.mark.skipif(not treeview_available(), reason="needs PyQt6 for ete4 to draw trees")


@needs_treeview
def test_write_tree_pages(tmp_path, capsys):
    make_results(tmp_path)
    assert tree_pages.write_tree_pages(str(tmp_path)) == 2
    assert "1 gene trees could not be drawn" in capsys.readouterr().out
    colour = dict(zip(["Sp_A", "Sp_B"], tree_pages.PALETTE))

    svg = (tmp_path / "Resolved_Gene_Trees" / "Tree_Images" / "OG0000000.svg").read_text()
    for gene in ("a1", "a2", "b1"):
        assert gene in svg
    assert "Sp_A_a1" not in svg  # the gene name, the species is its colour
    assert colour["Sp_A"].lower() in svg.lower() and colour["Sp_B"].lower() in svg.lower()
    assert "n1" in svg  # the gene tree node labels
    assert not (tmp_path / "Resolved_Gene_Trees" / "Tree_Images" / "OG0000009.svg").exists()
    species_svg = (tmp_path / "Species_Tree" / "SpeciesTree_rooted_node_labels.svg").read_text()
    assert "Sp_A" in species_svg and "N0" in species_svg

    page = (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").read_text()
    data = json.loads(re.search(r"const D = (.*);\n", page).group(1))
    assert data["trees"] == [["OG0000000", ["a1", "a2", "b1"]], ["OG0000001", ["a3", "b2", "b3"]]]
    assert data["species"] == [["Sp_A", colour["Sp_A"]], ["Sp_B", colour["Sp_B"]]]
    assert data["species_tree"] == "../Species_Tree/SpeciesTree_rooted_node_labels.svg"


@needs_treeview
def test_no_gene_trees_nothing_written(tmp_path):
    (tmp_path / "Orthogroups").mkdir()  # e.g. -og: no gene trees
    assert tree_pages.write_tree_pages(str(tmp_path)) == 0
    assert [p.name for p in tmp_path.iterdir()] == ["Orthogroups"]


@needs_treeview
@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the page's script")
def test_page_script_is_valid_javascript(tmp_path):
    make_results(tmp_path)
    tree_pages.write_tree_pages(str(tmp_path))
    js = tmp_path / "page.js"
    js.write_text((tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").read_text()
                  .split("<script>")[1].split("</script>")[0])
    result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_without_pyqt6(tmp_path, monkeypatch):
    make_results(tmp_path)
    monkeypatch.setitem(__import__("sys").modules, "ete4.treeview", None)  # as if PyQt6 were not installed
    with pytest.raises(tree_pages.TreeViewUnavailable):
        tree_pages.write_tree_pages(str(tmp_path))
    assert not (tmp_path / "Resolved_Gene_Trees" / "Tree_Images").exists()


@needs_treeview
def test_empty_gene_trees_file_nothing_written(tmp_path):
    make_results(tmp_path)
    (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("")  # as after -og
    assert tree_pages.write_tree_pages(str(tmp_path)) == 0
    assert sorted(p.name for p in (tmp_path / "Resolved_Gene_Trees").iterdir()) == ["Resolved_Gene_Trees.txt"]
    assert not (tmp_path / "Species_Tree" / "SpeciesTree_rooted_node_labels.svg").exists()


def test_missing_pyqt6_message_says_why_and_how():
    msg = tree_pages.missing_pyqt6_message("/data/Results_Oct09")
    assert "Qt graphics library (PyQt6)" in msg and "only this option needs it" in msg
    assert "pip install PyQt6" in msg and 'pip install "orthofinder[trees]"' in msg
    assert "python -m orthofinder.phylo.tree_pages /data/Results_Oct09" in msg


def test_command_without_pyqt6(tmp_path, monkeypatch, capsys):
    make_results(tmp_path)
    monkeypatch.setitem(__import__("sys").modules, "ete4.treeview", None)
    assert tree_pages.main([str(tmp_path)]) == 1
    assert "pip install PyQt6" in capsys.readouterr().out


@needs_treeview
def test_command(tmp_path, capsys):
    make_results(tmp_path)
    assert tree_pages.main([str(tmp_path)]) == 0
    assert "Gene tree pictures (2)" in capsys.readouterr().out
    assert (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").exists()
