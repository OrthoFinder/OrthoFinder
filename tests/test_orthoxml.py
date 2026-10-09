"""-x: the orthogroups, HOGs and orthologs as OrthoXML and as web pages."""
import gzip
import json
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import pytest

from orthofinder.utils import orthoxml

NS = {"o": orthoxml.NAMESPACE}

SPECIES = [("Species_A", "Species_A.faa"), ("Species B", "Species B.fa")]
N_SEQS = {0: 3, 1: 2, 5: 9}  # species 5 is not in the analysis
ID_DICT = {"0_0": "a0", "0_1": "a1&x", "0_2": "a2</script>", "1_0": "b0", "1_1": "b1", "5_0": "unused"}
OGS = [{"0_0", "1_1"}, {"0_1"}, {"0_2", "1_0", "0_1"}]  # OG0000001 has one gene: not a group


def page_data(path):
    page = path.read_text()
    assert page.count("</script>") == 1  # a gene name cannot end the page's script
    return json.loads(re.search(r"const D = (.*);\n", page).group(1))


def parse(path):
    root = ET.parse(path).getroot()
    assert root.tag == f"{{{orthoxml.NAMESPACE}}}orthoXML"
    genes = {g.get("id"): (s.get("name"), g.get("geneId"))
             for s in root.findall("o:species", NS) for g in s.iterfind("o:database/o:genes/o:gene", NS)}
    groups = {g.get("id"): ({p.get("name"): p.get("value") for p in g.findall("o:property", NS)},
                            sorted(genes[r.get("id")] for r in g.findall("o:geneRef", NS)))
              for g in root.findall("o:groups/o:orthologGroup", NS)}
    return root, genes, groups


def test_species_for_analysis(tmp_path):
    ids = tmp_path / "SpeciesIDs.txt"
    ids.write_text("0: Species_A.faa\n#1: Excluded.fa\n2: Species B.fa\n\n")
    assert orthoxml.species_for_analysis(str(ids), [2, 0]) == [
        ("Species B", "Species B.fa"), ("Species_A", "Species_A.faa")]


def test_orthogroups(tmp_path):
    base = str(tmp_path / "Orthogroups")
    orthoxml.write_results(SPECIES, OGS, N_SEQS, ID_DICT, [0, 1], base)

    root, genes, groups = parse(base + ".orthoxml")
    assert root.get("origin") == "OrthoFinder" and root.get("version") == "0.3"
    species = root.findall("o:species", NS)
    assert [(s.get("name"), s.get("NCBITaxId")) for s in species] == [("Species_A", "0"), ("Species B", "0")]
    databases = [s.find("o:database", NS) for s in species]
    assert [(d.get("name"), d.get("version")) for d in databases] == [("Species_A.faa", "NA"), ("Species B.fa", "NA")]
    # every gene of the species used, once, with a unique id
    assert {i: g for i, (_, g) in genes.items()} == {"0": "a0", "1": "a1&x", "2": "a2</script>", "3": "b0", "4": "b1"}
    # the orthogroups, named as in Orthogroups.tsv; single-gene orthogroups are not groups
    assert {k: sorted(g for _, g in v) for k, (_, v) in groups.items()} == {
        "OG0000000": ["a0", "b1"], "OG0000002": ["a1&x", "a2</script>", "b0"]}

    data = page_data(tmp_path / "Orthogroups.html")
    view = data["views"][0]
    assert view["columns"] == ["@Orthogroup", "#Genes", "*Species_A", "*Species B"]
    assert view["rows"] == [["OG0000000", 2, ["a0"], ["b1"]], ["OG0000002", 3, ["a1&x", "a2</script>"], ["b0"]]]
    assert data["info"]["rows"] == [["Species_A", "Species_A.faa", 3, "100.0%"], ["Species B", "Species B.fa", 2, "100.0%"]]


def make_results(tmp_path, gz=True):
    """Orthogroups.orthoxml and small HOG and ortholog files, as OrthoFinder writes them."""
    species = [("Sp_A", "Sp_A.fa"), ("Sp_B", "Sp_B.fa"), ("Sp_C", "Sp_C.fa")]
    ids = {"0_0": "a0", "0_1": "a1", "1_0": "b0", "1_1": "b1", "1_2": "b2", "2_0": "c0"}
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir()
    orthoxml.write_results(species, [{"0_0", "1_0", "2_0"}, {"0_1", "1_1"}], {0: 2, 1: 3, 2: 1}, ids, [0, 1, 2],
                           str(og_dir / "Orthogroups"))
    hog_dir = tmp_path / "Phylogenetic_Hierarchical_Orthogroups"
    hog_dir.mkdir()
    head = "HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\tSp_C\n"
    (hog_dir / "N1.tsv").write_text(head + "N1.HOG0000000\tOG0000000\tn0\ta0\tb0\t\n"
                                           "N1.HOG0000001\tOG0000001\tn3\ta1\tb1\t\n"
                                           "N1.HOG0000002\tOG0000002\tn5\t\tb2\t\n")  # a single gene
    (hog_dir / "N2.tsv").write_text(head)  # a node without HOGs
    orth_dir = tmp_path / "Orthologues"
    orth_dir.mkdir()
    files = {  # each relationship is in both species' files
        "Sp_A": "Orthogroup\tSpecies\tSp_A\tOrthologs\nOG0000000\tSp_B\ta0\tb0\nOG0000000\tSp_C\ta0\tc0\n"
                "OG0000001\tSp_B\ta1\tb1\n",
        "Sp_B": "Orthogroup\tSpecies\tSp_B\tOrthologs\nOG0000000\tSp_A\tb0\ta0\nOG0000001\tSp_A\tb1\ta1\n"
                "OG0000000\tSp_C\tb0, b1\tc0\n",
        "Sp_C": "Orthogroup\tSpecies\tSp_C\tOrthologs\nOG0000000\tSp_A\tc0\ta0\nOG0000000\tSp_B\tc0\tb0, b1\n",
    }
    for name, text in files.items():
        if gz:
            with gzip.open(orth_dir / f"{name}.tsv.gz", "wt") as f:
                f.write(text)
        else:
            (orth_dir / f"{name}.tsv").write_text(text)
    return og_dir / "Orthogroups.orthoxml", hog_dir, orth_dir


@pytest.mark.parametrize("gz", [True, False])
def test_hogs_and_orthologs(tmp_path, gz):
    og_xml, hog_dir, orth_dir = make_results(tmp_path, gz)
    _, og_genes, _ = parse(og_xml)
    written = orthoxml.write_hogs_and_orthologs(str(og_xml), str(hog_dir), str(orth_dir))
    assert sorted(p.rsplit("/", 1)[1] for p in written) == [
        "HOGs.html", "N1.orthoxml", "Orthologues.html", "Orthologues.orthoxml"]  # no N2: it has no HOGs

    # HOGs: the groups of N1.tsv (OrthoXML groups have two or more genes), the species of its
    # clade, gene ids as in Orthogroups.orthoxml
    root, genes, groups = parse(hog_dir / "N1.orthoxml")
    assert [s.get("name") for s in root.findall("o:species", NS)] == ["Sp_A", "Sp_B"]
    assert all(og_genes[i] == g for i, g in genes.items())
    assert groups == {
        "N1.HOG0000000": ({"orthogroup": "OG0000000", "gene tree parent clade": "n0"}, [("Sp_A", "a0"), ("Sp_B", "b0")]),
        "N1.HOG0000001": ({"orthogroup": "OG0000001", "gene tree parent clade": "n3"}, [("Sp_A", "a1"), ("Sp_B", "b1")])}
    view = page_data(hog_dir / "HOGs.html")["views"][0]
    assert view["name"] == "N1: Sp_A, Sp_B"
    assert view["rows"][0] == ["N1.HOG0000000", "OG0000000", "n0", 2, ["a0"], ["b0"]]
    assert view["rows"][2] == ["N1.HOG0000002", "OG0000002", "n5", 1, [], ["b2"]]  # the page has every HOG

    # orthologs: each relationship once, by species pair
    root, genes, groups = parse(orth_dir / "Orthologues.orthoxml")
    assert all(og_genes[i] == g for i, g in genes.items())
    assert sorted((p["orthogroup"], tuple(g)) for p, g in groups.values()) == [
        ("OG0000000", (("Sp_A", "a0"), ("Sp_B", "b0"))),
        ("OG0000000", (("Sp_A", "a0"), ("Sp_C", "c0"))),
        ("OG0000000", (("Sp_B", "b0"), ("Sp_B", "b1"), ("Sp_C", "c0"))),
        ("OG0000001", (("Sp_A", "a1"), ("Sp_B", "b1")))]
    views = page_data(orth_dir / "Orthologues.html")["views"]
    assert [v["name"] for v in views] == ["Sp_A - Sp_B", "Sp_A - Sp_C", "Sp_B - Sp_C"]
    assert views[2]["rows"] == [["OG0000000", "many-to-one", ["b0", "b1"], ["c0"]]]


def test_nothing_to_write(tmp_path):
    og_xml, _, _ = make_results(tmp_path)
    # e.g. after -og, which removes the HOG and ortholog directories
    assert orthoxml.write_hogs_and_orthologs(str(og_xml), str(tmp_path / "none"), str(tmp_path / "none2")) == []


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the pages' script")
def test_pages_script_is_valid_javascript(tmp_path):
    og_xml, hog_dir, orth_dir = make_results(tmp_path)
    orthoxml.write_hogs_and_orthologs(str(og_xml), str(hog_dir), str(orth_dir))
    for page in (tmp_path / "Orthogroups" / "Orthogroups.html", hog_dir / "HOGs.html", orth_dir / "Orthologues.html"):
        js = tmp_path / "page.js"
        js.write_text(page.read_text().split("<script>")[1].split("</script>")[0])
        result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
        assert result.returncode == 0, (page, result.stderr)
