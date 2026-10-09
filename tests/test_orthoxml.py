"""-x: the orthogroups as OrthoXML (Orthogroups.orthoxml) and as a web page (Orthogroups.html)."""
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


def test_species_for_analysis(tmp_path):
    ids = tmp_path / "SpeciesIDs.txt"
    ids.write_text("0: Species_A.faa\n#1: Excluded.fa\n2: Species B.fa\n\n")
    assert orthoxml.species_for_analysis(str(ids), [2, 0]) == [
        ("Species B", "Species B.fa"), ("Species_A", "Species_A.faa")]


def test_write_results(tmp_path):
    base = str(tmp_path / "Orthogroups")
    orthoxml.write_results(SPECIES, OGS, N_SEQS, ID_DICT, [0, 1], base)

    # the OrthoXML
    root = ET.parse(base + ".orthoxml").getroot()
    assert root.tag == f"{{{orthoxml.NAMESPACE}}}orthoXML"
    assert root.get("origin") == "OrthoFinder" and root.get("version") == "0.3"
    species = root.findall("o:species", NS)
    assert [(s.get("name"), s.get("NCBITaxId")) for s in species] == [("Species_A", "0"), ("Species B", "0")]
    databases = [s.find("o:database", NS) for s in species]
    assert [(d.get("name"), d.get("version")) for d in databases] == [("Species_A.faa", "NA"), ("Species B.fa", "NA")]
    # every gene of the species used, once, with a unique id
    genes = {g.get("id"): g.get("geneId") for s in species for g in s.iterfind("o:database/o:genes/o:gene", NS)}
    assert genes == {"0": "a0", "1": "a1&x", "2": "a2</script>", "3": "b0", "4": "b1"}
    # the orthogroups, named as in Orthogroups.tsv; single-gene orthogroups are not groups
    groups = {g.get("id"): sorted(genes[r.get("id")] for r in g.findall("o:geneRef", NS))
              for g in root.findall("o:groups/o:orthologGroup", NS)}
    assert groups == {"OG0000000": ["a0", "b1"], "OG0000002": ["a1&x", "a2</script>", "b0"]}

    # the web page: the same orthogroups, genes by species
    page = (tmp_path / "Orthogroups.html").read_text()
    assert page.count("</script>") == 1  # a gene name cannot end the page's script
    data = json.loads(re.search(r"const D = (.*);\n", page).group(1))
    assert [s["name"] for s in data["species"]] == ["Species_A", "Species B"]
    assert [(s["n_genes"], s["n_assigned"]) for s in data["species"]] == [(3, 3), (2, 2)]
    assert data["groups"] == [["OG0000000", [["a0"], ["b1"]]],
                              ["OG0000002", [["a1&x", "a2</script>"], ["b0"]]]]


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the page's script")
def test_html_script_is_valid_javascript(tmp_path):
    base = str(tmp_path / "Orthogroups")
    orthoxml.write_results(SPECIES, OGS, N_SEQS, ID_DICT, [0, 1], base)
    script = (tmp_path / "Orthogroups.html").read_text().split("<script>")[1].split("</script>")[0]
    js = tmp_path / "page.js"
    js.write_text(script)
    result = subprocess.run(["node", "--check", str(js)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
