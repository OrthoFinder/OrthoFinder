"""OrthoXML output (-x): the species file, the XML written, and the -x argument."""
import xml.etree.ElementTree as ET

import pytest

from orthofinder.run import process_args
from orthofinder.utils import orthoxml

NS = {"o": orthoxml.NAMESPACE}

SPECIES_FILE = (
    "Mycoplasma_agalactiae.faa\tMycoplasma agalactiae\t2110\tEnsembl\tEnsembl 46\n"
    "\n"
    "Mycoplasma_genitalium.fa\tMycoplasma genitalium\t2097\tUniProt\t2024_01\n"
)


def write(tmp_path, text, name="species.tsv"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# --- the species file ---------------------------------------------------------

def test_read_species_file(tmp_path):
    species = orthoxml.read_species_file(write(tmp_path, SPECIES_FILE))
    # keyed by the species name OrthoFinder uses: the FASTA filename without its extension
    assert species == {
        "Mycoplasma_agalactiae": ("Mycoplasma agalactiae", "2110", "Ensembl", "Ensembl 46"),
        "Mycoplasma_genitalium": ("Mycoplasma genitalium", "2097", "UniProt", "2024_01"),
    }


@pytest.mark.parametrize("text, message", [
    ("a.fa\tA\t1\tdb\n", "Line 1 "),                                   # 4 fields
    ("a.fa\tA\t1\tdb\tv\textra\n", "Line 1 "),                         # 6 fields
    ("a.fa\tA\t1\tdb\tv\nb.fa\tB\t\tdb\tv\n", "Line 2 "),              # an empty field
    ("a.fa\tA\tHomo sapiens\tdb\tv\n", "should be a number"),          # taxonomy ID not a number
    ("a.fa\tA\t1\tdb\tv\na.faa\tA2\t2\tdb\tv\n", "more than one line"),  # same species twice
    ("\n\n", "is empty"),
])
def test_read_species_file_errors(tmp_path, text, message):
    with pytest.raises(orthoxml.SpeciesFileError, match=message):
        orthoxml.read_species_file(write(tmp_path, text))


def test_read_species_file_missing(tmp_path):
    with pytest.raises(orthoxml.SpeciesFileError, match="Cannot read"):
        orthoxml.read_species_file(str(tmp_path / "no_such_file.tsv"))


def test_species_for_analysis(tmp_path):
    fn = write(tmp_path, SPECIES_FILE)
    names = {0: "Mycoplasma_agalactiae", 1: "Mycoplasma_genitalium", 2: "Other"}
    # in the order of the species used; species not in the analysis are not needed
    assert orthoxml.species_for_analysis(fn, names, [1, 0]) == [
        ("Mycoplasma genitalium", "2097", "UniProt", "2024_01"),
        ("Mycoplasma agalactiae", "2110", "Ensembl", "Ensembl 46"),
    ]
    with pytest.raises(orthoxml.SpeciesFileError, match="no line for:\n  Other"):
        orthoxml.species_for_analysis(fn, names, [0, 1, 2])


# --- the XML --------------------------------------------------------------------

def test_write_orthoxml(tmp_path):
    species_xml = [("Species A", "1", "db & co", "v<1>"), ("Species B", "2", "db", 'v"2"')]
    n_seqs = {0: 3, 1: 2, 5: 9}                  # species 5 is not in the analysis
    id_dict = {"0_0": "a0", "0_1": "a1&x", "0_2": "a<2>", "1_0": "b0", "1_1": "b1", "5_0": "unused"}
    ogs = [{"0_0", "1_1"}, {"0_1"}, {"0_2", "1_0", "0_1"}]
    fn = str(tmp_path / "Orthogroups.orthoxml")
    orthoxml.write_orthoxml(species_xml, ogs, n_seqs, id_dict, [0, 1], fn)

    root = ET.parse(fn).getroot()
    assert root.tag == f"{{{orthoxml.NAMESPACE}}}orthoXML"
    assert root.get("origin") == "OrthoFinder" and root.get("version") == "0.3"

    species = root.findall("o:species", NS)
    assert [(s.get("name"), s.get("NCBITaxId")) for s in species] == [("Species A", "1"), ("Species B", "2")]
    databases = [s.find("o:database", NS) for s in species]
    assert [(d.get("name"), d.get("version")) for d in databases] == [("db & co", "v<1>"), ("db", 'v"2"')]

    # every gene of the species used, once, with a unique id
    genes = {g.get("id"): g.get("geneId") for s in species for g in s.iterfind("o:database/o:genes/o:gene", NS)}
    assert genes == {"0": "a0", "1": "a1&x", "2": "a<2>", "3": "b0", "4": "b1"}

    # the orthogroups, named as in Orthogroups.tsv; single-gene orthogroups are not groups
    groups = {g.get("id"): sorted(genes[r.get("id")] for r in g.findall("o:geneRef", NS))
              for g in root.findall("o:groups/o:orthologGroup", NS)}
    assert groups == {"OG0000000": ["a0", "b1"], "OG0000002": ["a1&x", "a<2>", "b0"]}


# --- the -x argument -----------------------------------------------------------

def test_x_argument_bad_file_stops_at_once(tmp_path, capsys):
    bad = write(tmp_path, "a.fa\tA\tnot_a_number\tdb\tv\n")
    with pytest.raises(SystemExit):
        process_args.process_args(["-f", str(tmp_path), "-x", bad])
    assert "the NCBI taxonomy ID should be a number" in capsys.readouterr().out


def test_x_argument_missing_file_name(capsys):
    with pytest.raises(SystemExit):
        process_args.process_args(["-x"])
    assert "Missing option for command line argument -x" in capsys.readouterr().out
