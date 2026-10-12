"""--orthoxml (-x): the orthogroups and HOGs as OrthoXML."""
import xml.etree.ElementTree as ET

from orthofinder.utils import orthoxml

NS = {"o": orthoxml.NAMESPACE}

SPECIES = [("Species_A", "Species_A.faa"), ("Species B", "Species B.fa")]
N_SEQS = {0: 3, 1: 2, 5: 9}  # species 5 is not in the analysis
ID_DICT = {"0_0": "a0", "0_1": "a1&x", "0_2": "a2</script>", "1_0": "b0", "1_1": "b1", "5_0": "unused"}
OGS = [{"0_0", "1_1"}, {"0_1"}, {"0_2", "1_0", "0_1"}]  # OG0000001 has one gene: not a group


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
    assert sorted(p.name for p in tmp_path.iterdir()) == ["Orthogroups.orthoxml"]  # only the OrthoXML


def make_results(tmp_path):
    """Orthogroups.orthoxml and small HOG files, as OrthoFinder writes them."""
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
    return og_dir / "Orthogroups.orthoxml", hog_dir


def test_hogs(tmp_path):
    og_xml, hog_dir = make_results(tmp_path)
    _, og_genes, _ = parse(og_xml)
    written = orthoxml.write_hog_files(str(og_xml), str(hog_dir))
    assert [p.rsplit("/", 1)[1] for p in written] == ["N1.orthoxml"]  # no N2: it has no HOGs

    # HOGs: the groups of N1.tsv (OrthoXML groups have two or more genes), the species of its
    # clade, gene ids as in Orthogroups.orthoxml
    root, genes, groups = parse(hog_dir / "N1.orthoxml")
    assert [s.get("name") for s in root.findall("o:species", NS)] == ["Sp_A", "Sp_B"]
    assert all(og_genes[i] == g for i, g in genes.items())
    assert groups == {
        "N1.HOG0000000": ({"orthogroup": "OG0000000", "gene tree parent clade": "n0"}, [("Sp_A", "a0"), ("Sp_B", "b0")]),
        "N1.HOG0000001": ({"orthogroup": "OG0000001", "gene tree parent clade": "n3"}, [("Sp_A", "a1"), ("Sp_B", "b1")])}


def test_nothing_to_write(tmp_path):
    og_xml, _ = make_results(tmp_path)
    # e.g. after -og, which removes the HOG directory
    assert orthoxml.write_hog_files(str(og_xml), str(tmp_path / "none")) == []
