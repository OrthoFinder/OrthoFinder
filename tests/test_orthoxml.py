"""--orthoxml (-x): the orthogroups and HOGs as OrthoXML."""
import xml.etree.ElementTree as ET

from orthofinder.utils import orthoxml

NS = {"o": orthoxml.NAMESPACE}


def parse(path):
    root = ET.parse(path).getroot()
    assert root.tag == f"{{{orthoxml.NAMESPACE}}}orthoXML"
    genes = {g.get("id"): (s.get("name"), g.get("geneId"))
             for s in root.findall("o:species", NS) for g in s.iterfind("o:database/o:genes/o:gene", NS)}
    return root, genes


def write_table(path, species, rows, first="OG", start=0):
    """An Orthogroups.tsv-like table: a row for each list of genes, a gene in the column of its species."""
    lines = ["Orthogroup\t" + "\t".join(name for name, _ in species)]
    for i, genes in enumerate(rows, start):
        lines.append(f"{first}{i:07d}\t" + "\t".join(
            ", ".join(g for g in genes if g.startswith(prefix)) for _, prefix in species))
    path.write_text("\n".join(lines) + "\n")


def test_species_for_analysis(tmp_path):
    ids = tmp_path / "SpeciesIDs.txt"
    ids.write_text("0: Species_A.faa\n#1: Excluded.fa\n2: Species B.fa\n\n", encoding="utf-8")
    # (the names of the other results files: gene_ids.clean_species_name)
    assert orthoxml.species_for_analysis(str(ids), [2, 0]) == [
        ("Species_B", "Species B.fa"), ("Species_A", "Species_A.faa")]


def test_orthogroups(tmp_path):
    """Without the HOGs (e.g. -og): a group for each orthogroup of Orthogroups.tsv."""
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir()
    columns = [("Species_A", "a"), ("Species B", "b")]
    write_table(og_dir / "Orthogroups.tsv", columns, [["a0", "b1"], ["a2</script>", "b0", "a1&x"], ["a3"]])
    write_table(og_dir / "Orthogroups_UnassignedGenes.tsv", columns, [["b2"]], start=3)
    species = [("Species_A", "Species_A.faa"), ("Species B", "Species B.fa")]

    out = orthoxml.write_results(str(tmp_path), species)
    assert out == str(og_dir / "Orthogroups.orthoxml")
    root, genes = parse(out)
    assert root.get("origin") == "OrthoFinder" and root.get("version") == "0.3"
    assert root.find("o:taxonomy", NS) is None
    species_el = root.findall("o:species", NS)
    assert [(s.get("name"), s.get("NCBITaxId")) for s in species_el] == [("Species_A", "0"), ("Species B", "0")]
    databases = [s.find("o:database", NS) for s in species_el]
    assert [(d.get("name"), d.get("version")) for d in databases] == [("Species_A.faa", "NA"), ("Species B.fa", "NA")]
    # every gene, also those in no orthogroup, once, with a unique id
    assert sorted(genes.values()) == sorted([("Species_A", g) for g in ("a0", "a2</script>", "a1&x", "a3")]
                                            + [("Species B", g) for g in ("b0", "b1", "b2")])
    # the orthogroups, named as in Orthogroups.tsv; single-gene orthogroups are not groups
    groups = {g.get("id"): sorted(genes[r.get("id")][1] for r in g.findall("o:geneRef", NS))
              for g in root.findall("o:groups/o:orthologGroup", NS)}
    assert groups == {"OG0000000": ["a0", "b1"], "OG0000001": ["a1&x", "a2</script>", "b0"]}


def test_no_orthogroups(tmp_path):
    """No Orthogroups.tsv in the results (e.g. -ogt, -ost): nothing is written."""
    assert orthoxml.write_results(str(tmp_path), []) is None
    assert list(tmp_path.iterdir()) == []


def hog_results(tmp_path):
    """Results with a species tree ((Sp_A,Sp_B)N1,Sp_C)N0: Orthogroups.tsv and N1.tsv."""
    columns = [("Sp_A", "a"), ("Sp_B", "b"), ("Sp_C", "c")]
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir()
    write_table(og_dir / "Orthogroups.tsv", columns,
                [["a0", "a1", "b0", "b1", "c0"],   # a duplication on the branch to N1
                 ["a2", "b2", "c1", "c2"],         # in-paralogs in Sp_C
                 ["a3", "b3"],                     # genes of N1 only
                 ["a4", "c3"],                     # no HOG at N1 (e.g. no gene tree)
                 ["c4", "c5"]])                    # genes of Sp_C only, in-paralogs
    write_table(og_dir / "Orthogroups_UnassignedGenes.tsv", columns, [["b4"]], start=5)
    (tmp_path / "Species_Tree").mkdir()
    (tmp_path / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt").write_text("((Sp_A:1,Sp_B:1)N1:1,Sp_C:1)N0;\n")
    hog_dir = tmp_path / "Phylogenetic_Hierarchical_Orthogroups"
    hog_dir.mkdir()
    (hog_dir / "N1.tsv").write_text("HOG\tOG\tGene Tree Parent Clade\tSp_A\tSp_B\tSp_C\n"
                                     "N1.HOG0000000\tOG0000000\tn1\ta0\tb0\t\n"
                                     "N1.HOG0000001\tOG0000000\tn2\ta1\tb1\t\n"
                                     "N1.HOG0000002\tOG0000001\tn1\ta2\tb2\t\n"
                                     "N1.HOG0000003\tOG0000002\tn0\ta3\tb3\t\n")
    return [(name, name + ".fa") for name, _ in columns]


def shape(el, genes, taxa):
    """A group as (tag, id, taxon name, {score: value}, {property: value}, [members]), a gene as its name."""
    tag = el.tag.split("}")[1]
    if tag == "geneRef":
        return genes[el.get("id")][1]
    members = [shape(m, genes, taxa) for m in el if m.tag.split("}")[1] in ("geneRef", "orthologGroup", "paralogGroup")]
    scores = {s.get("id"): float(s.get("value")) for s in el.findall("o:score", NS)}
    properties = {p.get("name"): p.get("value") for p in el.findall("o:property", NS)}
    return tag, el.get("id"), taxa.get(el.get("taxonId")), scores, properties, members


def test_with_hogs(tmp_path):
    species = hog_results(tmp_path)
    out = orthoxml.write_results(str(tmp_path), species)
    assert out == str(tmp_path / "Orthogroups" / "Orthogroups.orthoxml")
    assert sorted(p.name for p in (tmp_path / "Orthogroups").iterdir()) == [
        "Orthogroups.orthoxml", "Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"]

    root, genes = parse(out)
    assert root.get("version") == "0.5"
    assert ("Sp_B", "b4") in genes.values()   # every gene, also those in no orthogroup
    # the species tree as the taxonomy, the species with their taxon
    taxa = {t.get("id"): t.get("name") for t in root.iter(f"{{{orthoxml.NAMESPACE}}}taxon")}
    assert sorted(taxa.values()) == ["N0", "N1", "Sp_A", "Sp_B", "Sp_C"]
    assert [taxa[s.get("taxonId")] for s in root.findall("o:species", NS)] == ["Sp_A", "Sp_B", "Sp_C"]
    assert [d.get("id") for d in root.findall("o:scores/o:scoreDef", NS)] == ["CompletenessScore"]

    groups = [shape(g, genes, taxa) for g in root.find("o:groups", NS)]
    complete = {"CompletenessScore": 1.0}
    n1 = lambda hog, clade, *members: ("orthologGroup", hog, "N1", complete,
                                       {"TaxRange": "N1", "gene tree parent clade": clade}, list(members))
    n0 = {"TaxRange": "N0"}
    paralogs = lambda *members: ("paralogGroup", None, None, {}, {}, list(members))
    assert groups == [
        ("orthologGroup", "OG0000000", "N0", complete, n0, [
            paralogs(n1("N1.HOG0000000", "n1", "a0", "b0"), n1("N1.HOG0000001", "n2", "a1", "b1")),
            "c0"]),
        ("orthologGroup", "OG0000001", "N0", complete, n0, [n1("N1.HOG0000002", "n1", "a2", "b2"), paralogs("c1", "c2")]),
        # a single member: the HOG of N1, with the orthogroup's id and its own as "HOG"
        ("orthologGroup", "OG0000002", "N1", complete,
         {"TaxRange": "N1", "gene tree parent clade": "n0", "HOG": "N1.HOG0000003"}, ["a3", "b3"]),
        # not in a HOG of N1: in the orthogroup's group directly; 2 of the 3 species
        ("orthologGroup", "OG0000003", "N0", {"CompletenessScore": 0.6667}, n0, ["c3", "a4"]),
        # a root HOG is an orthologGroup, also of a single member (as in OMA and orthoxml-tools)
        ("orthologGroup", "OG0000004", "N0", {"CompletenessScore": 0.3333}, n0, [paralogs("c4", "c5")]),
    ]


def test_without_species_tree(tmp_path):
    """HOG files but no species tree: the orthogroups only."""
    species = hog_results(tmp_path)
    (tmp_path / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt").unlink()
    root, _ = parse(orthoxml.write_results(str(tmp_path), species))
    assert root.get("version") == "0.3" and root.find("o:taxonomy", NS) is None
    assert [g.get("id") for g in root.findall("o:groups/o:orthologGroup", NS)] == [
        "OG0000000", "OG0000001", "OG0000002", "OG0000003", "OG0000004"]
