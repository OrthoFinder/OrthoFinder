"""--stats-view: a page of the statistics of a run, with plots (drawn in the browser)."""
import json
import os
import re
import shutil
import subprocess

import pytest

from orthofinder.orthology import stats_page


def make_results(d):
    (d / "Orthogroups").mkdir(parents=True)
    (d / "Orthogroups" / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\tSp_A\tSp_B\tTotal\n"
        "OG0000000\t2\t1\t3\n"
        "OG0000001\t1\t1\t2\n"
        "OG0000002\t3\t0\t3\n")
    stats = d / "Comparative_Genomics_Statistics"
    stats.mkdir()
    (stats / "Statistics_Overall.tsv").write_text(
        "Number of species\t2\nNumber of genes\t10\nPercentage of genes in orthogroups\t80.0\nDate\t2026-10-09\n\n"
        "Average number of genes per-species in orthogroup\tNumber of orthogroups\n'1\t3\n")
    (stats / "Statistics_PerSpecies.tsv").write_text(
        "\tSp_A\tSp_B\nNumber of genes in orthogroups\t6\t2\nNumber of unassigned genes\t1\t1\n\n"
        "Number of genes per-species in orthogroup\tNumber of orthogroups\tNumber of orthogroups\n'1\t1\t2\n")
    (stats / "Duplications_per_Orthogroup.tsv").write_text(
        "Orthogroup\tDuplications (all)\tDuplications (50% support)\nOG0000000\t2\t1\n")
    (stats / "Duplications_per_Species_Tree_Node.tsv").write_text(
        "Species Tree Node\tDuplications (all)\tDuplications (50% support)\nSp_A\t2\t1\nN0\t0\t0\n")
    (stats / "OrthologsStats_one-to-one.tsv").write_text("\tSp_A\tSp_B\nSp_A\t0.0\t1.0\nSp_B\t1.0\t0.0\n")
    (d / "Species_Tree").mkdir()
    (d / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt").write_text("(Sp_A:0.5,Sp_B:0.5)N0;")
    hogs = d / "Phylogenetic_Hierarchical_Orthogroups"
    hogs.mkdir()
    (hogs / "N0.tsv").write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\n"
        "N0.HOG0000000\tOG0000000\tn0\ta1, a2\tb1\n"
        "N0.HOG0000001\tOG0000002\tn0\ta3, a4, a5\t\n")


def page_data(d):
    page = (d / "Comparative_Genomics_Statistics" / "Statistics.html").read_text()
    return json.loads(re.search(r"const D = (.*);\n", page).group(1))


def test_write_stats_page(tmp_path):
    make_results(tmp_path)
    page = stats_page.write_stats_page(str(tmp_path))
    assert page == str(tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html")
    data = page_data(tmp_path)
    assert [sp for sp, _ in data["species"]] == ["Sp_A", "Sp_B"]
    # only the first block of the statistics files; numbers as numbers
    assert data["overall"] == [["Number of species", 2], ["Number of genes", 10],
                               ["Percentage of genes in orthogroups", 80.0], ["Date", "2026-10-09"]]
    assert data["per_species"] == [["Number of genes in orthogroups", [6, 2]], ["Number of unassigned genes", [1, 1]]]
    # [orthogroup, genes of each species, duplications, with support >= 0.5, has a tree page]
    assert data["orthogroups"] == [["OG0000000", [2, 1], 2, 1, 0], ["OG0000001", [1, 1], 0, 0, 0],
                                   ["OG0000002", [3, 0], 0, 0, 0]]
    assert dict(data["og_sizes"])["2"] == 1 and dict(data["og_sizes"])["3"] == 2
    assert data["matrices"] == {"Orthogroups shared": [[3, 2], [2, 2]], "Orthologs: one-to-one": [[0, 1], [1, 0]]}
    assert data["node_duplications"] == {"Sp_A": [2, 1], "N0": [0, 0]}
    assert data["species_tree"]["genes"] == [["Sp_A", 0], ["Sp_B", 1]]
    assert data["hogs"] == [{"node": "N0", "hogs": 2, "genes": 6, "multi_species": 1,
                             "sizes": [[label, {"3": 2}.get(label, 0)] for label, _ in data["hogs"][0]["sizes"]]}]
    assert "tree_page" not in data


def test_hog_data(tmp_path):
    make_results(tmp_path)
    stats_page.write_stats_page(str(tmp_path))
    text = (tmp_path / "Comparative_Genomics_Statistics" / "Statistics_Data" / "HOGs_N0.js").read_text()
    assert text.startswith("hogsLoaded(") and text.endswith(");\n")
    assert json.loads(text[len("hogsLoaded("):-len(");\n")]) == {
        "node": "N0", "hogs": [["N0.HOG0000000", "OG0000000", "n0", [2, 1]], ["N0.HOG0000001", "OG0000002", "n0", [3, 0]]]}


def test_links_to_the_tree_page(tmp_path):
    make_results(tmp_path)
    (tmp_path / "Resolved_Gene_Trees").mkdir()
    (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").write_text("")
    (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("OG0000000: (a:1,b:1);\n")
    stats_page.write_stats_page(str(tmp_path))
    data = page_data(tmp_path)
    assert data["tree_page"] == "../Resolved_Gene_Trees/Resolved_Gene_Trees.html"
    assert [og[4] for og in data["orthogroups"]] == [1, 0, 0]  # only the orthogroups with a tree


def test_only_gene_counts(tmp_path):
    # e.g. an early stop: no statistics files, species tree or HOGs; the page is still written
    (tmp_path / "Orthogroups").mkdir()
    (tmp_path / "Orthogroups" / "Orthogroups.GeneCount.tsv").write_text("Orthogroup\tSp_A\tSp_B\tTotal\nOG0000000\t1\t1\t2\n")
    assert stats_page.write_stats_page(str(tmp_path))
    data = page_data(tmp_path)
    assert data["hogs"] == [] and "overall" not in data and "species_tree" not in data


def test_no_orthogroups(tmp_path):
    assert stats_page.write_stats_page(str(tmp_path)) is None
    assert list(tmp_path.iterdir()) == []


def test_number():
    assert stats_page.number("81.0") == 81.0 and stats_page.number("2,733") == 2733
    assert stats_page.number("'1") == 1 and stats_page.number("Orthogroups.tsv") == "Orthogroups.tsv"


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the page's script")
def test_page_script_is_valid_javascript(tmp_path):
    make_results(tmp_path)
    stats_page.write_stats_page(str(tmp_path))
    js = tmp_path / "page.js"
    js.write_text((tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html").read_text()
                  .split("<script>")[1].split("</script>")[0])
    result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_command(tmp_path, capsys):
    make_results(tmp_path)
    assert stats_page.main([str(tmp_path)]) == 0
    assert "Statistics page:" in capsys.readouterr().out


@pytest.mark.skipif(shutil.which("orthofinder") is None, reason="needs the orthofinder command")
def test_options_in_the_help():  # the pages and --view
    result = subprocess.run(["orthofinder", "-h"], capture_output=True, text=True, check=False,
                            env={**os.environ, "COLUMNS": "200"})
    for option in ("--tree-view", "--stats-view", "--groups-view", "--view", "--orthoxml"):
        assert option in result.stdout
