"""--groups-view: a page of the genes of the orthogroups and HOGs, a table loaded per tab."""
import gzip
import json
import re
import shutil
import subprocess

import pytest

from orthofinder.orthology import groups_page


def make_results(d):
    og = d / "Orthogroups"
    og.mkdir(parents=True)
    (og / "Orthogroups.tsv").write_text(
        "Orthogroup\tSp_A\tSp_B\tSp_C\n"
        "OG0000000\ta0\tb0, b1\tc0\n"
        "OG0000001\ta1\tb2\t\n")
    (og / "Orthogroups_UnassignedGenes.tsv").write_text("Orthogroup\tSp_A\tSp_B\tSp_C\nOG0000002\t\t\tc1\n")
    hogs = d / "Phylogenetic_Hierarchical_Orthogroups"
    hogs.mkdir()
    (hogs / "N1.tsv").write_text("HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\tSp_C\n"
                                 "N1.HOG0000000\tOG0000000\tn0\ta0\tb0, b1\t\n"
                                 "N1.HOG0000001\tOG0000001\tn2\ta1\tb2\t\n")
    (hogs / "N0.tsv").write_text("HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\tSp_C\n"
                                 "N0.HOG0000000\tOG0000000\tn0\ta0\tb0, b1\tc0\n")
    orth = d / "Orthologs"
    orth.mkdir()
    for name, text in {
        "Sp_A": "Orthogroup\tSpecies\tSp_A\tOrthologs\nOG0000000\tSp_B\ta0\tb0, b1\nOG0000000\tSp_C\ta0\tc0\n",
        "Sp_B": "Orthogroup\tSpecies\tSp_B\tOrthologs\nOG0000000\tSp_A\tb0, b1\ta0\n",
        "Sp_C": "Orthogroup\tSpecies\tSp_C\tOrthologs\nOG0000000\tSp_A\tc0\ta0\n",
    }.items():
        with gzip.open(orth / f"{name}.tsv.gz", "wt") as f:
            f.write(text)
    xen = d / "Putative_Xenologs"
    xen.mkdir()
    (xen / "Sp_A.tsv").write_text("Orthogroup\tSp_A\tOther\nOG0000001\ta1\tc1\n")  # the other species: from its genes
    (xen / "Sp_B.tsv").write_text("Orthogroup\tSp_B\tOther\n")
    (xen / "Sp_C.tsv").write_text("Orthogroup\tSp_C\tOther\n")


def page_data(d):
    page = (d / "Orthogroups" / "Groups.html").read_text()
    assert page.count("</script>") == 2  # the theme's and the page's: a gene name cannot end the page's script
    return json.loads(re.search(r"const D = (.*);\n", page).group(1))


def data_file(path):
    """The value of a page's data file (view_pages.write_data): ofData(name, value);"""
    text = path.read_text()
    m = re.match(r'ofData\("([^"]+)",(.*)\);\n$', text, re.S)
    assert m and m.group(1) == path.stem
    return json.loads(m.group(2))


def table(d, key):
    return data_file(d / "Orthogroups" / "Groups_Data" / f"{key}.js")


def test_orthogroups_and_hogs(tmp_path):
    make_results(tmp_path)
    assert groups_page.write_groups_page(str(tmp_path)) == str(tmp_path / "Orthogroups" / "Groups.html")
    data = page_data(tmp_path)
    assert [t["id"] for t in data["tabs"]] == ["orthogroups", "hogs"]
    assert [sp for sp, _ in data["species"]] == ["Sp_A", "Sp_B", "Sp_C"]
    # a table: its rows' IDs, genes and species, the species of its columns, and its genes in chunks
    assert table(tmp_path, "Orthogroups") == {"ids": {"prefix": "OG", "width": 7, "n": 2}, "cols": [0, 1, 2],
                                              "tree": [0, 0], "genes": [4, 2], "species": [3, 2], "chunks": [0]}
    assert table(tmp_path, "Orthogroups_0") == [[[0, ["a0"]], [1, ["b0", "b1"]], [2, ["c0"]]], [[0, ["a1"]], [1, ["b2"]]]]
    # the genes, sorted, each with its orthogroup's row: to find a gene by the start of its name
    assert data["genes"] == ["a0"]
    assert table(tmp_path, "Genes_0") == [["a0", 0], ["a1", 1], ["b0", 0], ["b1", 0], ["b2", 1], ["c0", 0]]
    # the HOGs: a table for each species tree node, in node order, with the species of its clade
    hogs = data["tabs"][1]["views"]
    assert [(v["name"], v["key"], v["rows"]) for v in hogs] == [
        ("N0: Sp_A, Sp_B, Sp_C", "HOGs_N0", 1), ("N1: Sp_A, Sp_B", "HOGs_N1", 2)]
    assert table(tmp_path, "HOGs_N1") == {"ids": {"prefix": "N1.HOG", "width": 7, "n": 2}, "ogs": [0, 1],
                                          "clades": ["n0", "n2"], "cols": [0, 1], "genes": [3, 2], "species": [2, 2],
                                          "chunks": [0]}
    assert table(tmp_path, "HOGs_N1_0")[0] == [[0, ["a0"]], [1, ["b0", "b1"]]]


def test_links(tmp_path):
    make_results(tmp_path)
    (tmp_path / "Resolved_Gene_Trees").mkdir()
    (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("OG0000000: (a:1,b:1);\n")
    groups_page.write_groups_page(str(tmp_path), written={"trees", "stats", "groups"})  # as a run with --view
    data = page_data(tmp_path)
    assert data["tree_page"] == "../Resolved_Gene_Trees/Resolved_Gene_Trees.html"
    assert table(tmp_path, "Orthogroups")["tree"] == [1, 0]  # the orthogroups with a tree
    assert data["stats_page"] == "../Comparative_Genomics_Statistics/Statistics.html"
    assert data["nav"] == [["Gene trees", "../Resolved_Gene_Trees/Resolved_Gene_Trees.html", 0, "trees"],
                           ["Statistics", "../Comparative_Genomics_Statistics/Statistics.html", 0, "stats"],
                           ["Orthogroups & HOGs", "../Orthogroups/Groups.html", 1, "groups"]]


def test_no_orthogroups(tmp_path):
    assert groups_page.write_groups_page(str(tmp_path)) is None
    assert list(tmp_path.iterdir()) == []


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the page's script")
def test_page_script_is_valid_javascript(tmp_path):
    make_results(tmp_path)
    groups_page.write_groups_page(str(tmp_path), written={"trees", "stats", "groups"})
    js = tmp_path / "page.js"
    js.write_text((tmp_path / "Orthogroups" / "Groups.html").read_text().split("<script>")[-1].split("</script>")[0])
    result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_command(tmp_path, capsys):
    make_results(tmp_path)
    assert groups_page.main([str(tmp_path)]) == 0
    assert "Orthogroups and HOGs page:" in capsys.readouterr().out
