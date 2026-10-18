"""
Gene and species names from FASTA headers and file names of all sorts: each is
cleaned once (gene_ids.clean_name) and written the same in every results file.
"""
import pytest

from orthofinder.orthology import groups_page
from orthofinder.utils import fasta_processor, gene_ids, orthoxml


@pytest.mark.parametrize("header, name", [
    ("sp|P12345|GENE_HUMAN Protein X OS=Homo sapiens", "sp|P12345|GENE_HUMAN_Protein_X_OS_Homo_sapiens"),
    ('"q1"', "_q1_"),
    ('g"x', "g_x"),
    ("it's", "it_s"),
    ("name with spaces; semi", "name_with_spaces__semi"),
    ("gene:1", "gene_1"), ("gene,1", "gene_1"), ("x[1]", "x_1_"), ("(a)=b", "_a__b"),
    ("tab\there", "tab_here"), ("nul\x00x", "nul_x"), ("nbsp x", "nbsp_x"), ("ls x", "ls_x"),
    ("α-tubulin µ", "α-tubulin_µ"), ("🧬gene", "🧬gene"),
    ("<script>alert(1)</script>", "<script>alert_1_</script>"), ("a&b", "a&b"),
    ("back\\slash", "back\\slash"), ("#hash", "#hash"),
    ("  trailing   ", "trailing"),
])
def test_clean_name(header, name):
    assert gene_ids.clean_name(header) == name


def test_clean_species_name():
    # also a file name (Orthologs/<species>.tsv): nor the characters Windows does not allow
    assert gene_ids.clean_species_name('Myco "agalactiae".faa') == "Myco__agalactiae_"
    assert gene_ids.clean_species_name("Myco gallisepticum µ.fa") == "Myco_gallisepticum_µ"
    assert gene_ids.clean_species_name("a*b?c<d>e|f\\g.fasta") == "a_b_c_d_e_f_g"
    assert gene_ids.clean_species_name("Homo.sapiens.v2.fa") == "Homo.sapiens.v2"


def test_gene_names_unique():
    names = gene_ids.GeneNames()
    assert [names.add(h) for h in ["gene:1", "gene,1", "gene 1", "gene_1_2", "x"]] == \
        ["gene_1", "gene_1_2", "gene_1_3", "gene_1_2_2", "x"]


def test_extractors_agree(tmp_path):
    ids = tmp_path / "SequenceIDs.txt"
    ids.write_text("0_0: g\"x1 first protein\n0_1: α:2 second\n1_0: g\"x1 other species\n", encoding="utf-8")
    first = gene_ids.FirstWordExtractor(str(ids)).get_id_to_name_dict()
    full = gene_ids.FullAccession(str(ids)).get_id_to_name_dict()
    assert first == {"0_0": "g_x1", "0_1": "α_2", "1_0": "g_x1"}
    assert full == {"0_0": "g_x1_first_protein", "0_1": "α_2_second", "1_0": "g_x1_other_species"}
    ids.write_text("0_0: x[1] a\n0_1: x_1_ b\n", encoding="utf-8")
    with pytest.raises(RuntimeError) as e:   # the same first word once cleaned: FullAccession is used
        gene_ids.FirstWordExtractor(str(ids))
    assert not str(e.value).startswith("ERROR")


def test_species_names(tmp_path):
    ids = tmp_path / "SpeciesIDs.txt"
    ids.write_text('0: Myco "agalactiae".faa\n#1: Excluded x.fa\n2: B: colon.fa\n', encoding="utf-8")
    assert gene_ids.species_names(str(ids)) == {"0": "Myco__agalactiae_", "1": "Excluded_x", "2": "B__colon"}
    assert gene_ids.species_names(str(ids), q_unused=False) == {"0": "Myco__agalactiae_", "2": "B__colon"}


@pytest.mark.parametrize("data, headers", [
    ("﻿>a1 x\nMK\n>a2\nLL\n".encode("utf-8"), ["a1 x", "a2"]),       # byte order mark
    (b">a1\r\nMK\r\n>a2\r\nLL\r\n", ["a1", "a2"]),                           # Windows line endings
    (b">a1\rMK\r>a2\rLL\r", ["a1", "a2"]),                                    # old Mac line endings
    (">α µ\nMK\n".encode("utf-8"), ["α µ"]),
])
def test_input_fasta_forms(tmp_path, data, headers):
    fn = tmp_path / "in.fa"
    fn.write_bytes(data)
    with fasta_processor.open_input_fasta(str(fn)) as f:
        lines = list(f)
    assert [fasta_processor.input_header(l[1:].rstrip()) for l in lines if l.startswith(">")] == \
        [(h, True) for h in headers]
    assert all(not l.startswith(">") and l.rstrip("\n") in ("MK", "LL") for l in lines if not l.startswith(">"))


def test_input_fasta_not_utf8(tmp_path):
    # e.g. saved by a Windows program as cp1252: read, as cp1252, and reported
    fn = tmp_path / "in.fa"
    fn.write_bytes(">µ-gene caf\xe9\nMK\n".encode("latin-1"))
    with fasta_processor.open_input_fasta(str(fn)) as f:
        header = next(f)[1:].rstrip()
    assert fasta_processor.input_header(header) == ("µ-gene café", False)


def test_read_groups(tmp_path):
    names = [gene_ids.clean_name(h) for h in ['"q1"', "it's", "α µ", "a&b", "<b>"]]
    og = tmp_path / "Orthogroups.tsv"
    og.write_text("Orthogroup\tMyco__agalactiae_\tB\nOG0000000\t%s\tb1\n" % ", ".join(names), encoding="utf-8")
    hogs = tmp_path / "N1.tsv"
    hogs.write_text("HOG\tOG\tGene Tree Parent Clade\tMyco__agalactiae_\tB\n"
                    "N1.HOG0000000\tOG0000000\tn1\t%s\tb1\n" % ", ".join(names), encoding="utf-8")
    index = {"Myco__agalactiae_": 0, "B": 1}
    header, rows = groups_page.read_groups(str(og), index)
    assert header[1] == "Myco__agalactiae_" and rows[0][1] == [(0, names), (1, ["b1"])]
    _, rows = groups_page.read_groups(str(hogs), index, first=3)
    assert rows[0][1] == [(0, names), (1, ["b1"])]


def test_orthoxml_reads_tables(tmp_path):
    names = [gene_ids.clean_name(h) for h in ['"q1"', "it's", "α µ", "a&b"]]
    hogs = tmp_path / "N1.tsv"
    hogs.write_text("HOG\tOG\tGene Tree Parent Clade\tA\nN1.HOG0000000\tOG0000000\tn1\t%s\n" % ", ".join(names),
                    encoding="utf-8")
    ids = {"A": {name: i for i, name in enumerate(names)}}
    assert orthoxml.read_table(str(hogs), 3, ids) == [("N1.HOG0000000", "n1", {"A": [0, 1, 2, 3]})]


def test_user_species_tree_names(tmp_path):
    # A user's species tree may name the species as in their FASTA file names
    from orthofinder.run import previous_run
    from orthofinder.gene_tree import infer_trees
    fn = tmp_path / "SpeciesTree.txt"
    fn.write_text("(('Homo sapiens',Myco_gallisepticum_µ),'it''s');\n", encoding="utf-8-sig")
    species = {"0": "Homo_sapiens", "1": "Myco_gallisepticum_µ", "2": "it_s"}
    previous_run.check_user_species_tree(str(fn), list(species.values()))
    out = tmp_path / "ids.txt"
    infer_trees.convert_user_species_tree(str(fn), species, str(out))
    from orthofinder.phylo import tree
    assert sorted(tree.Tree(str(out)).get_leaf_names()) == ["0", "1", "2"]


def test_names_cache(tmp_path):
    # Parsed once and shared while used; never kept by the cache itself; a changed file is parsed again
    import gc
    ids = tmp_path / "SequenceIDs.txt"
    ids.write_text("0_0: a x\n0_1: b y\n", encoding="utf-8")
    a = gene_ids.FirstWordExtractor(str(ids)).get_id_to_name_dict()
    assert gene_ids.FirstWordExtractor(str(ids)).get_id_to_name_dict() is a
    del a
    gc.collect()
    assert not [k for k in gene_ids._names_cache.keys() if k[1] == str(ids.resolve())]
    # first words now the same (and the size changed: SequenceIDs.txt is only appended to,
    # and some filesystems keep the time of a change only to the second)
    ids.write_text("0_0: a x\n0_1: a yy\n", encoding="utf-8")
    for _ in range(2):    # (the error is remembered, and raised each time)
        with pytest.raises(RuntimeError) as e:
            gene_ids.FirstWordExtractor(str(ids))
        assert not str(e.value).startswith("ERROR")
    assert gene_ids.FullAccession(str(ids)).get_id_to_name_dict() == {"0_0": "a_x", "0_1": "a_yy"}


def test_files_of_older_versions(tmp_path):
    # Older versions wrote SpeciesIDs.txt etc. in the system's encoding (cp1252 on
    # Windows): each line is read as UTF-8 or else cp1252, also in one file (--assign
    # appends UTF-8 lines to it), and a byte order mark is dropped
    from orthofinder.utils import file_io
    ids = tmp_path / "SpeciesIDs.txt"
    ids.write_bytes("﻿0: Myco µ.fa\n".encode("utf-8") + "1: café µ.fa\n".encode("cp1252") + "2: α.fa\n".encode("utf-8"))
    with file_io.open_text_any(str(ids)) as f:
        assert list(f) == ["0: Myco µ.fa\n", "1: café µ.fa\n", "2: α.fa\n"]
    assert gene_ids.species_names(str(ids)) == {"0": "Myco_µ", "1": "café_µ", "2": "α"}
    tree = tmp_path / "SpeciesTree.txt"
    tree.write_bytes("('Myco µ','café µ');\n".encode("cp1252"))     # e.g. a user's tree saved on Windows
    from orthofinder.phylo import tree as tree_lib
    assert sorted(gene_ids.user_tree_species(n) for n in tree_lib.Tree(str(tree)).get_leaf_names()) == \
        ["Myco_µ", "café_µ"]


def test_species_prefixes_where_a_file_gives_no_species():
    # trees, alignments, FASTA files, Orthogroups.txt, Duplications.tsv: <species>_<gene>
    prefixes = gene_ids.species_prefixes({0: "Homo.sapiens", 1: "Mus_musculus"})
    assert prefixes == {"0": "Homo_sapiens_", "1": "Mus_musculus_"}
    assert gene_ids.with_species({"0_3": "g1", "1_0": "g1"}, prefixes) == {"0_3": "Homo_sapiens_g1", "1_0": "Mus_musculus_g1"}


def test_xenologue_other_genes_have_their_species():
    # Putative_Xenologs/<species>.tsv: its species' genes alone (the column is the
    # species), the "Other" genes (of any species) as <species>_<gene>
    from collections import defaultdict
    from orthofinder.orthology import tree_analyser

    def genes(sp, *g):
        d = defaultdict(list)
        d[sp] = list(g)
        return d
    names = {"0_1": "a1", "1_2": "b2", "1_3": "b3"}
    other = gene_ids.with_species(names, gene_ids.species_prefixes({0: "A", 1: "B"}))
    olog, sus = ["", ""], ["", ""]
    # the genes of species 0 are suspect: xenologues of those of species 1
    orthologs = [(defaultdict(list), genes("1", "2"), genes("0", "1"), defaultdict(list))]
    tree_analyser.get_lines_for_ortholog_files([(5, orthologs)], {"0": "A", "1": "B"}, [0, 1], names, True,
                                               olog, sus, other_names=other)
    assert sus == ["OG0000005\ta1\tB_b2\n", "OG0000005\tb2\tA_a1\n"]


def test_gene_names_unique_as_the_smallest_free_number():
    # the same names as checking every "_<n>" from 2 (also with names that are a suffix already)
    import random

    class Reference(object):
        def __init__(self):
            self.used = set()

        def add(self, accession):
            name = gene_ids.clean_name(accession)
            if name in self.used:
                n = 2
                while "%s_%d" % (name, n) in self.used:
                    n += 1
                name = "%s_%d" % (name, n)
            self.used.add(name)
            return name
    rng = random.Random(6)
    for _ in range(50):
        headers = [rng.choice(["g", "g:1", "g,1", "g_1", "g_2", "g_1_2", "x[1]", "x_1_", "g 3"]) for _ in range(200)]
        names, reference = gene_ids.GeneNames(), Reference()
        assert [names.add(h) for h in headers] == [reference.add(h) for h in headers]
    many = gene_ids.GeneNames()
    assert [many.add("gene") for _ in range(20000)][-1] == "gene_20000"    # (linear: was 2 x 10^8 checks)


def _old_full_accession_names(fn):
    names, names_in_species = {}, {}
    for line in open(fn, encoding="utf-8"):
        line = line.rstrip("\r\n")
        if not line.strip():
            continue
        id, accession = line.split(": ", 1)
        id = id.replace("#", "").strip()
        if id in names:
            raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
        names_in_species.setdefault(id.split("_")[0], gene_ids.GeneNames())
        names[id] = names_in_species[id.split("_")[0]].add(accession)
    return names


def _old_first_word_names(fn):
    names, accs = {}, []
    for line in open(fn, encoding="utf-8"):
        if not line.strip():
            continue
        id, rest = line.split(": ", 1)
        acc, i_sp = gene_ids.clean_name(rest.split(None, 1)[0]), int(id.split("_")[0])
        while len(accs) < i_sp + 1:
            accs.append(set())
        if acc in accs[i_sp]:
            raise RuntimeError("A duplicate accession was found using just first part: % s" % acc)
        accs[i_sp].add(acc)
        if id in names:
            raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
        names[id] = acc
    return names


def test_gene_name_tables_as_the_dicts_they_replace(tmp_path):
    # the names, their order (SequenceIDs.txt's: unassigned genes are numbered in it) and
    # the errors (a duplicate ID; first words not unique: FullAccession then) of the dicts
    import random
    rng = random.Random(11)
    fn = tmp_path / "SequenceIDs.txt"

    def outcome(parse):
        try:
            names = parse(str(fn))
            return "ok", list(names.items()), len(names)
        except RuntimeError as e:
            return "error", str(e)
    for _ in range(400):
        lines = ["%d_%d: %s %s\n" % (sp, q, rng.choice(["g1", "x:1", "x,1", "a b", "名", "g%d" % q, "u%d_%d" % (sp, q)]),
                                     rng.choice(["", "desc"]))
                 for sp in range(rng.randint(1, 4)) for q in range(rng.randint(0, 8))]
        if lines and rng.random() < 0.2:
            lines.append(rng.choice(lines))
        if rng.random() < 0.4:
            rng.shuffle(lines)
        fn.write_text("".join(lines), encoding="utf-8")
        assert outcome(gene_ids._full_accession_names) == outcome(_old_full_accession_names)
        assert outcome(gene_ids._first_word_names) == outcome(_old_first_word_names)


def test_gene_name_table_lookups_and_species_view():
    pairs = [("0_0", "a"), ("0_2", "名"), ("1_0", "b"), ("1_1", "")]
    table = gene_ids.GeneNameTable.build(pairs)
    assert dict(table) == dict(pairs) and list(table) == [k for k, _ in pairs]
    for key in ("0_1", "01_0", "0_00", "1_0_0", "+1_0", "x", "2_0"):
        assert key not in table and table.get(key) is None
    view = gene_ids.SpeciesNames(table, [1], {"0": "A_", "1": "B_"})
    assert dict(view) == {"1_0": "B_b", "1_1": "B_"} and "0_0" not in view
    view.update({"1": "SpeciesB"})          # (as the callers add the species' names)
    assert view["1"] == "SpeciesB" and len(view) == 3
