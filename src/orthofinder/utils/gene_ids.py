"""Reading gene identifiers and names from sequence accessions."""
import os
import re
import weakref
from collections.abc import Mapping

import numpy as np

from . import file_io

# The characters a name cannot hold, each replaced by "_": whitespace and control
# characters (a tab or line break would split a cell or line of the tab-separated
# files, a space a Newick label), those with a meaning in Newick (: ; ( ) , [ ] = and
# the quote '), and the quote " (spreadsheets and csv readers take a cell starting
# with it as quoted). A name is then written as it is in every file (tables, trees,
# FASTA, pages) on every OS, and no file needs quoting.
_UNSAFE_IN_NAME = re.compile(r"""[\s\x00-\x1f\x7f-\x9f:;(),\[\]='"]""")
# A species name is also a file name (e.g. Orthologs/<species>.tsv): nor those Windows
# does not allow in file names
_UNSAFE_IN_FILE_NAME = re.compile(r"""[*?<>|\\/]""")


def clean_name(text):
    """A FASTA header (or its first word) as a gene name: see _UNSAFE_IN_NAME."""
    return _UNSAFE_IN_NAME.sub("_", text.strip())


def clean_species_label(label):
    """A species' name as it is written in the results (clean_name), also e.g. a leaf of a user's species tree."""
    return _UNSAFE_IN_FILE_NAME.sub("_", clean_name(label))


def user_tree_species(label):
    """
    A leaf of a user's species tree as a species' name (clean_species_label): written
    as in the results, as in the FASTA file name, or as a quoted Newick label ('Homo
    sapiens', with '' for a quote).
    """
    if len(label) > 1 and label[0] == label[-1] == "'":
        label = label[1:-1].replace("''", "'")
    return clean_species_label(label)


def clean_species_name(fasta_filename):
    """The name of a species: its FASTA file name without the extension (clean_species_label)."""
    return clean_species_label(fasta_filename.strip().rsplit(".", 1)[0])


def species_prefixes(species_names_dict):
    """
    {species ID (str): the start of its genes' names in the files that do not give
    the species of each gene} (unless -X): its name, "." and " " as "_", then "_".
    A file that gives the species (a column per species, e.g. Orthogroups.tsv)
    has the gene names alone; the others (trees, alignments, FASTA files,
    Orthogroups.txt, Duplications.tsv, the other genes of Putative_Xenologs) have
    <species>_<gene>.
    """
    return {str(k): v.replace(".", "_").replace(" ", "_") + "_" for k, v in species_names_dict.items()}


def with_species(names, prefixes):
    """{sequence ID: name} as {sequence ID: <species>_<name>} (species_prefixes)."""
    return {seq_id: prefixes[seq_id.split("_", 1)[0]] + name for seq_id, name in names.items()}


def species_names(species_ids_fn, q_unused=True):
    """
    {species ID (str): name} of SpeciesIDs.txt ("<ID>: <FASTA filename>" lines);
    q_unused: also the species not used in this analysis ("#<ID>: ..." lines).
    """
    names = {}
    with file_io.open_text_any(species_ids_fn) as infile:
        for line in infile:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if not q_unused:
                    continue
                line = line[1:]
            if not line.strip():
                continue
            i_species, fasta_filename = line.split(": ", 1)
            names[i_species.strip()] = clean_species_name(fasta_filename)
    return names


def get_id_pair_from_string(line):
    return list(map(int, line.split("_")))


class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different
    ways."""

    def get_id_to_name_dict(self):
        raise NotImplementedError("Should not be implemented")


class GeneNames(object):
    """
    The gene names of one species (clean_name), kept unique: two headers that
    differ only in the characters a name cannot hold ("gene:1" and "gene,1",
    or "x[1]" and "x_1_") would give the same name, and the genes could not be
    told apart in the results: the later one gets a suffix ("gene_1_2"), the
    next free "_<n>" in the order of the input file.
    """

    def __init__(self):
        self.used = set()
        # name -> the next n to try for its "_<n>": the numbers below it are all used
        # (names are never removed), so each copy of a name is not checked from 2 again
        # (many copies of one name: n^2 checks)
        self.next_suffix = {}

    def add(self, accession):
        """The unique name of the next gene with this header."""
        name = clean_name(accession)
        if name in self.used:
            n = self.next_suffix.get(name, 2)
            while "%s_%d" % (name, n) in self.used:
                n += 1
            self.next_suffix[name] = n + 1
            name = "%s_%d" % (name, n)
        self.used.add(name)
        return name


def _canonical_id(key):
    """(species, sequence) of a key that is a gene's ID "<species>_<sequence>" exactly, else None."""
    sp, sep, seq = key.partition("_")
    if not sep or not sp.isdigit() or not seq.isdigit() or \
            (len(sp) > 1 and sp[0] == "0") or (len(seq) > 1 and seq[0] == "0"):
        return None
    return int(sp), int(seq)


class GeneNameTable(Mapping):
    """
    {sequence ID ("<species>_<sequence>"): gene name}, as a dict gives them, kept
    compactly: the names of each species as one UTF-8 buffer and the offsets of its
    names in it (numpy), not a pair of strings a gene (~145 bytes a gene in a dict:
    GBs with 10s of millions of genes, ~25 here). Above all, a forked process (the
    tree analysis) reads a name without changing anything it inherited: a lookup in
    a dict changes the reference counts of the strings it finds, so the process
    gets its own copy of most of the dict's memory, one per process.
    In the order of the IDs file (as a dict filled from it); other keys (and names
    set, e.g. the species: update) are kept in a dict, after the genes.
    """
    __slots__ = ("_buffers", "_offsets", "_present", "_order", "_n", "_extra", "__weakref__")

    def __init__(self):
        self._buffers = {}      # species -> bytes of its names
        self._offsets = {}      # species -> numpy offsets (n + 1) of the names in the buffer
        self._present = {}      # species -> numpy bool (the sequences with a name), if not all
        self._order = None      # (species, sequence) arrays of the IDs in their order, if not species-by-species, sequence order
        self._n = 0
        self._extra = {}

    @classmethod
    def build(cls, entries):
        """The table of (ID, name) pairs in their order (each ID once), in any order."""
        builder = _TableBuilder(natural_only=False)
        for key, name in entries:
            builder.add(key, name)
        return builder.finish()

    def _name(self, sp, seq):
        offsets = self._offsets.get(sp)
        if offsets is None or seq >= len(offsets) - 1:
            return None
        present = self._present.get(sp)
        if present is not None and not present[seq]:
            return None
        return self._buffers[sp][offsets[seq]:offsets[seq + 1]].decode("utf-8")

    def __getitem__(self, key):
        if self._extra:
            value = self._extra.get(key, self)
            if value is not self:
                return value
        numbers = _canonical_id(key)
        name = None if numbers is None else self._name(*numbers)
        if name is None:
            raise KeyError(key)
        return name

    def __contains__(self, key):
        if key in self._extra:
            return True
        numbers = _canonical_id(key) if isinstance(key, str) else None
        return numbers is not None and self._name(*numbers) is not None

    def _gene_ids(self):
        if self._order is not None:
            for sp, seq in zip(self._order[0].tolist(), self._order[1].tolist()):
                yield sp, seq
            return
        for sp, offsets in self._offsets.items():
            present = self._present.get(sp)
            for seq in range(len(offsets) - 1):
                if present is None or present[seq]:
                    yield sp, seq

    def __iter__(self):
        # (a gene's ID where it is, also if its name was set again, as in a dict; then the other keys)
        for sp, seq in self._gene_ids():
            yield "%d_%d" % (sp, seq)
        for key in self._extra:
            numbers = _canonical_id(key)
            if numbers is None or self._name(*numbers) is None:
                yield key

    def __len__(self):
        overlap = sum(1 for key in self._extra if (lambda n: n is not None and self._name(*n) is not None)(_canonical_id(key)))
        return self._n + len(self._extra) - overlap

    def items(self):
        return ((key, self[key]) for key in self)

    def __setitem__(self, key, name):
        self._extra[key] = name

    def update(self, other=(), **kwargs):
        self._extra.update(other, **kwargs)

    def __eq__(self, other):
        return Mapping.__eq__(self, other)

    __hash__ = None


class SpeciesNames(Mapping):
    """
    The names of a GeneNameTable of some species only, as <prefix of its species><name>
    (species_prefixes) or the names alone (prefixes None): computed when read, not a
    second table of every gene. Other keys set (update) are kept in a dict.
    """
    __slots__ = ("_table", "_species", "_prefixes", "_extra")

    def __init__(self, table, species, prefixes=None):
        self._table = table
        self._species = set(map(str, species))
        self._prefixes = prefixes
        self._extra = {}

    def __getitem__(self, key):
        if self._extra:
            value = self._extra.get(key, self)
            if value is not self:
                return value
        sp = key.partition("_")[0]
        if sp not in self._species:
            raise KeyError(key)
        name = self._table[key]
        return name if self._prefixes is None else self._prefixes[sp] + name

    def __contains__(self, key):
        if key in self._extra:
            return True
        return isinstance(key, str) and key.partition("_")[0] in self._species and key in self._table

    def __iter__(self):
        own = set()
        for key in self._table:
            if key.partition("_")[0] in self._species:
                own.add(key) if key in self._extra else None
                yield key
        for key in self._extra:
            if key not in own:
                yield key

    def __len__(self):
        return sum(1 for _ in self)

    def items(self):
        return ((key, self[key]) for key in self)

    def __setitem__(self, key, name):
        self._extra[key] = name

    def update(self, other=(), **kwargs):
        self._extra.update(other, **kwargs)

    def __eq__(self, other):
        return Mapping.__eq__(self, other)

    __hash__ = None


# The names of SequenceIDs.txt are parsed once and shared, for as long as anything
# uses them: the parts of a run that need them at the same time (OrthogroupSet,
# the orthogroup files, ...) get the same GeneNameTable, not one each. The entries
# are weak, so the cache never keeps the names in memory itself. The table is
# shared: callers must not change it.
_names_cache = weakref.WeakValueDictionary()   # (kind, path, mtime, size) -> GeneNameTable
_names_errors = {}                             # (kind, path, mtime, size) -> the RuntimeError's args


def _cached_names(kind, ids_filename, parse):
    st = os.stat(ids_filename)
    key = (kind, os.path.realpath(ids_filename), st.st_mtime_ns, st.st_size)
    if key in _names_errors:
        raise RuntimeError(*_names_errors[key])
    names = _names_cache.get(key)
    if names is None:
        try:
            names = parse(ids_filename)
        except RuntimeError as error:
            # (e.g. the first words are not unique: FullAccession is used instead)
            _names_errors[key] = error.args
            raise
        _names_cache[key] = names
    return names


class _NotInOrder(Exception):
    """The IDs are not each species' together in increasing order (_TableBuilder)."""


class _TableBuilder(object):
    """
    A GeneNameTable, from the (ID, name) pairs of an IDs file in its order. With
    natural_only (the order OrthoFinder writes: each species' IDs together, its
    sequences in increasing order) each species is finished (encoded, its strings
    let go) as its IDs end, so only one species' names are strings at a time, not
    every gene's (GBs, kept by the process, which Python does not give back);
    another order raises _NotInOrder (the file is then read again, natural_only
    False: all species kept until the end, and the order of the IDs too).
    first_word_names: a species' names must be unique (else a RuntimeError, not
    starting "ERROR", FirstWordExtractor: FullAccession is used instead).
    """

    def __init__(self, natural_only, first_word_names=False):
        self.table = GeneNameTable()
        self.natural_only = natural_only
        self.first_word_names = first_word_names
        self.names = {}             # species -> its names by sequence (None: none), until finished
        self.finished = set()
        self.current, self.last_seq = None, -1
        self.order_sp, self.order_seq = [], []
        self.names_used = {}        # (first_word_names) species -> its names so far (natural_only: the current)

    def add(self, key, name):
        if self.first_word_names:
            # (checked as each ID comes, before the ID itself, as FirstWordExtractor always did)
            sp_text = key.split("_")[0]
            if self.natural_only and sp_text not in self.names_used:
                self.names_used.clear()
            used = self.names_used.setdefault(sp_text, set())
            if name in used:
                raise RuntimeError("A duplicate accession was found using just first part: % s" % name)
            used.add(name)
        numbers = _canonical_id(key)
        if numbers is None:
            if key in self.table._extra:
                raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % key)
            self.table._extra[key] = name
            return
        sp, seq = numbers
        names = self.names.get(sp)
        if names is not None and seq < len(names) and names[seq] is not None:
            raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % key)
        if self.natural_only:
            if sp != self.current:
                if sp in self.finished:
                    raise _NotInOrder()
                if self.current is not None:
                    self._finish(self.current)
                self.current, self.last_seq = sp, -1
            elif seq < self.last_seq:
                raise _NotInOrder()
        else:
            self.order_sp.append(sp)
            self.order_seq.append(seq)
        if names is None:
            names = self.names[sp] = []
        if seq >= len(names):
            names.extend([None] * (seq + 1 - len(names)))
        names[seq] = name
        self.last_seq = seq
        self.table._n += 1

    def _finish(self, sp):
        names = self.names.pop(sp)
        self.finished.add(sp)
        encoded = [b"" if n is None else n.encode("utf-8") for n in names]
        offsets = np.zeros(len(encoded) + 1, dtype=np.int64)
        np.cumsum(np.fromiter((len(e) for e in encoded), dtype=np.int64, count=len(encoded)), out=offsets[1:])
        table = self.table
        table._offsets[sp] = offsets.astype(np.int32) if offsets[-1] < 2 ** 31 else offsets
        table._buffers[sp] = b"".join(encoded)
        if any(n is None for n in names):
            table._present[sp] = np.fromiter((n is not None for n in names), dtype=bool, count=len(names))

    def finish(self):
        for sp in list(self.names):
            self._finish(sp)
        if not self.natural_only:
            order = (self.order_sp, self.order_seq)
            natural, seen, last = True, set(), (None, -1)
            for sp, seq in zip(*order):
                if sp != last[0]:
                    if sp in seen:
                        natural = False
                        break
                    seen.add(sp)
                elif seq < last[1]:
                    natural = False
                    break
                last = (sp, seq)
            if not natural:
                self.table._order = (np.array(order[0], dtype=np.int32), np.array(order[1], dtype=np.int32))
            else:
                # (the species in the order they came: a dict keeps its keys in that order)
                self.table._offsets = {sp: self.table._offsets[sp] for sp in dict.fromkeys(order[0])}
        else:
            # (the species finished in the order they came already)
            pass
        return self.table


def _full_accession_names(ids_filename):
    try:
        return _full_accession_table(ids_filename, natural_only=True)
    except _NotInOrder:
        return _full_accession_table(ids_filename, natural_only=False)


def _full_accession_table(ids_filename, natural_only):
    builder = _TableBuilder(natural_only)
    names_in_species = {}       # (natural_only: of the current species only)
    current = None
    with file_io.open_text_any(ids_filename) as ids_file:
        for line in ids_file:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            id, accession = line.split(": ", 1)
            id = id.replace("#", "")
            id = id.strip()
            # Names are unique within a species (see GeneNames)
            i_sp = id.split("_")[0]
            if natural_only and i_sp != current:
                names_in_species.pop(current, None)
                current = i_sp
            if i_sp not in names_in_species:
                names_in_species[i_sp] = GeneNames()
            builder.add(id, names_in_species[i_sp].add(accession))     # (a duplicate ID: a RuntimeError "ERROR...")
    return builder.finish()


def _first_word_names(ids_filename):
    try:
        return _first_word_table(ids_filename, natural_only=True)
    except _NotInOrder:
        return _first_word_table(ids_filename, natural_only=False)


def _first_word_table(ids_filename, natural_only):
    builder = _TableBuilder(natural_only, first_word_names=True)
    with file_io.open_text_any(ids_filename) as ids_file:
        for line in ids_file:
            if not line.strip():
                continue
            id, rest = line.split(": ", 1)
            builder.add(id, clean_name(rest.split(None, 1)[0]))
    return builder.finish()


class FullAccession(IDExtractor):
    """The genes named by their whole FASTA header (GeneNames)."""

    def __init__(self, ids_filename):
        self.id_to_name_dict = _cached_names("full", ids_filename, _full_accession_names)

    def get_id_to_name_dict(self):
        return self.id_to_name_dict


class FirstWordExtractor(IDExtractor):
    """
    The genes named by the first word of their FASTA header (clean_name). If two
    genes of a species would have the same name, a RuntimeError (not starting
    "ERROR"): the callers then use FullAccession.
    """

    def __init__(self, ids_filename):
        self.id_to_name_dict = _cached_names("first", ids_filename, _first_word_names)

    def get_id_to_name_dict(self):
        return self.id_to_name_dict
