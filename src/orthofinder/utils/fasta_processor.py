# -*- coding: utf-8 -*-
import sys
import gzip
import os.path
import glob
from collections import Counter
from rich.markup import escape
try:
    from rich import print
except ImportError:
    ...
from . import util, file_handler as files
from .printing import printer
from . import species_info
from . import gene_ids, file_io

# count = 0 

# Input FASTA files are those with these extensions.
FASTA_EXTENSIONS = {"fa", "faa", "fasta", "fas", "pep", "fna"}

# Telling nucleotide from protein sequences. Protein uses ~20 letters, DNA
# almost only A, C, G, T/U and N, so:
#   - a letter that is never a nucleotide (IUPAC) code means protein: E F I L
#     P Q (amino acids), J O Z (rare amino acids), X (unknown amino acid), *;
#   - otherwise the sequences are nucleotide if at least NUCLEOTIDE_FRACTION
#     of the letters are A, C, G, T, U or N, and protein if not (a protein
#     can lack all the letters above, e.g. a short or low-complexity one).
#     In protein files 22-37% of the letters are A/C/G/T/N (whole files;
#     measured on Mycoplasma, Chlamydomonas and the test proteomes), in DNA
#     ~100%: 0.75 leaves a wide margin, and allows DNA with many ambiguity
#     codes (R, Y, K, ...).
# Up to SEQUENCE_LINES_TO_CHECK lines / SEQUENCE_LETTERS_TO_CHECK letters of a
# file are looked at.
PROTEIN_ONLY_LETTERS = frozenset("EFILPQJOZX*")
NUCLEOTIDE_LETTERS = frozenset("ACGTUN")
NUCLEOTIDE_FRACTION = 0.75
SEQUENCE_LINES_TO_CHECK = 1000
SEQUENCE_LETTERS_TO_CHECK = 100000


def open_input_fasta(fn):
    """
    An input FASTA file, for reading as text in any of the forms it comes in: a
    byte order mark at the start (e.g. from Windows Notepad) is dropped, any line
    ending (LF, CRLF or CR) is read as LF, and bytes that are not UTF-8
    are kept (surrogateescape: see input_header), so that reading never fails.
    """
    return open(fn, encoding="utf-8-sig", errors="surrogateescape")


def input_header(text):
    """
    A header line of open_input_fasta, and whether it was UTF-8: if not (e.g. a
    file from a Windows program, in cp1252), it is read as cp1252.
    """
    if text.isascii():      # (nearly every header: no copy made)
        return text, True
    try:
        text.encode("utf-8")
        return text, True
    except UnicodeEncodeError:
        return text.encode("utf-8", "surrogateescape").decode("cp1252", errors="replace"), False


def is_input_fasta_name(f):
    """Whether a file name is that of an input FASTA file (by its extension; not macOS "._" files)."""
    return len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in FASTA_EXTENSIONS and not f.startswith("._")


class SequenceTypeSample(object):
    """
    The letters of the first sequence lines of a FASTA file, enough to tell
    nucleotide from protein sequences (see above). add() each sequence line
    as the file is read; it ignores lines once enough have been seen.
    """

    def __init__(self):
        self.counts = Counter()
        self.n_lines = 0
        self.n_letters = 0

    @property
    def full(self):
        return self.n_lines >= SEQUENCE_LINES_TO_CHECK or self.n_letters >= SEQUENCE_LETTERS_TO_CHECK

    def add(self, line):
        if self.full:
            return
        line = line.strip().upper().replace("-", "").replace(".", "")
        if line:
            self.counts.update(line)
            self.n_letters += len(line)
            self.n_lines += 1

    def type(self):
        """"protein", "dna" or "empty"."""
        if self.n_letters == 0:
            return "empty"
        if any(self.counts[c] for c in PROTEIN_ONLY_LETTERS):
            return "protein"
        n_nucleotide = sum(self.counts[c] for c in NUCLEOTIDE_LETTERS)
        return "dna" if n_nucleotide >= NUCLEOTIDE_FRACTION * self.n_letters else "protein"


def sequence_type(fn):
    """"protein", "dna" or "empty", from the first sequences of a FASTA file."""
    sample = SequenceTypeSample()
    # (only the sequences are looked at: a header that is not UTF-8 does not matter)
    with open(fn, encoding="utf-8", errors="replace") as infile:
        for line in infile:
            if not line.startswith(">"):
                sample.add(line)
                if sample.full:
                    break
    return sample.type()


class FastaWriter(object):
    def __init__(self, source_fasta_filename, q_use_only_second_part=False, q_glob=False, q_first_word=False, q_whole_line=False):
        self.seq_lists = dict()
        q_first = True
        accession = ""
        sequence = ""
        files = glob.glob(source_fasta_filename) if q_glob else [source_fasta_filename]
        for fn in files:
            with (gzip.open(fn, 'rt', encoding="utf-8") if fn.endswith('.gz') else open(fn, 'r', encoding="utf-8")) as fasta_file:
                for line in fasta_file:
                    if line[0] == ">":
                        # deal with old sequence
                        if not q_first:
                            self.seq_lists[accession] = sequence
                            sequence = ""
                        q_first = False
                        # get id for new sequence
                        if q_whole_line:
                            accession = line[1:].rstrip()
                        else:
                            accession = line[1:].rstrip().split()[0]
                        if q_use_only_second_part:
                            try:
                                accession = accession.split("|")[1]
                            except IndexError:     # (no "|": the header is reported, the accession kept)
                                sys.stderr.write(line + "\n")
                        elif q_first_word:
                            accession = accession.split()[0]
                    else:
                        sequence += line
            try:
                self.seq_lists[accession] = sequence
            except Exception:
                sys.stderr.write(accession + "\n")
    
    def write_seqs_to_fasta(self, seqs, out_filename):
        with open(out_filename, 'w', encoding="utf-8", newline="\n") as out_file:
            for seq in seqs:
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % seq)
                    out_file.write(self.seq_lists[seq])
                else:
                    sys.stderr.write("ERROR: %s not found\n" % seq)
    
    def write_seqs_to_fasta_with_new_accessions(self, seqs, out_filename, id_dict):
        with open(out_filename, 'w', encoding="utf-8", newline="\n") as out_file:
            for seq in seqs:
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % id_dict[seq])
                    out_file.write(self.seq_lists[seq])
                else:
                    sys.stderr.write(seq + "\n")


def processes_new_fasta(
        fasta_dir, 
        q_dna, 
        species_info_obj_prev = None, 
        species_to_use_prev_names=(),
        species_id_fn="",
        sequence_id_fn="",
    ):
    """
    Process fasta files and return a Directory object with all paths completed.
    """

    fasta_extensions = FASTA_EXTENSIONS
    # Check files present
    q_ok = True
    if not os.path.exists(fasta_dir):
        print("\nDirectory does not exist: %s" % fasta_dir)
        util.fail()
    files_in_directory = sorted([f for f in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir,f))])
    original_fasta_filenames = []
    excluded_files = []

    for f in files_in_directory:
        if is_input_fasta_name(f):
            original_fasta_filenames.append(f)
        else:
            excluded_files.append(f)

    if len(excluded_files) != 0:
        print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
        for f in excluded_files:
            print(f)
        print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fasta_extensions)))
    
    species_to_use_prev_names = set(species_to_use_prev_names)
    if len(original_fasta_filenames) + len(species_to_use_prev_names) < 2:
        print("ERROR: At least two species are required")
        util.fail()

    if any([fn in species_to_use_prev_names for fn in original_fasta_filenames]):
        print("ERROR: Attempted to add a second copy of a previously included species:")
        for fn in original_fasta_filenames:
            if fn in species_to_use_prev_names: print(fn)
        print("")
        util.fail()

    # A species is named by its file name without the extension (as in the
    # species tree and the results): "A.fa" and "A.fna" would both be "A".
    files_of_name = {}
    for fn in sorted(species_to_use_prev_names) + sorted(original_fasta_filenames):
        files_of_name.setdefault(gene_ids.clean_species_name(fn), []).append(fn)
    same_name = {name: fns for name, fns in files_of_name.items() if len(fns) > 1}
    if same_name:
        print("ERROR: Species must have different names. These files would give species of the same name "
              "(the file name without its extension, with the characters a name cannot hold replaced by _):")
        for name, fns in sorted(same_name.items()):
            print("  %s: %s" % (name, ", ".join(fns)))
        if species_to_use_prev_names:
            print("(including the species already in the analysis)")
        print("")
        util.fail()

    species_to_use_prev_names = sorted([*species_to_use_prev_names])
    original_fasta_filenames = sorted([*original_fasta_filenames])

    if len(original_fasta_filenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fasta_dir)
        util.fail()

    if species_info_obj_prev == None:
        species_info_obj = species_info.SpeciesInfo()
    else:
        species_info_obj = species_info_obj_prev

    if not species_id_fn:
        if files.handler.wd_current:
            species_id_fn = os.path.join(files.handler.wd_current, "SpeciesIDs.txt")
        else:
            species_id_fn = files.handler.get_species_ids_fn()

    if not sequence_id_fn:
        if files.handler.wd_current:
            sequence_id_fn = os.path.join(files.handler.wd_current, "SequenceIDs.txt")
        else:
            sequence_id_fn = files.handler.get_sequence_ids_fn()

    i_seq = 0
    i_species = 0
    # If it's a previous analysis:
    if len(species_to_use_prev_names) != 0:
        with file_io.open_text_any(species_id_fn) as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        i_species = int(line.split(":")[0]) + 1
    species_info_obj.i_first_new_species = i_species
    new_species_ids = []
    duplicated = False

    species_seen_dict = {}
    renamed = []      # (file, header, name): headers that clean to a name already used

    not_utf8 = []     # the files whose headers are not UTF-8 (read as cp1252: input_header)

    with open(sequence_id_fn, 'a', encoding="utf-8", newline="\n") as ids_file, \
            open(species_id_fn, 'a', encoding="utf-8", newline="\n") as species_file:
        for fasta_filename in original_fasta_filenames:
            new_species_ids.append(i_species)
            # (UTF-8: a byte of a sequence line that is not, which open_input_fasta keeps,
            # is written as U+FFFD, so that the file is read like any other)
            output_fasta = open(files.handler.get_species_fasta_fn(i_species, q_for_creation=True), 'w',
                                encoding="utf-8", errors="replace", newline="\n")
            fasta_filename = fasta_filename.rstrip()
            # (a file name need not be UTF-8 either: input_header)
            species_file.write("%d: %s\n" % (i_species, input_header(fasta_filename)[0]))
            base_filename, extension = os.path.splitext(fasta_filename)

            species_seen_dict[fasta_filename] = {}
            gene_names = gene_ids.GeneNames()     # the names the genes get (gene_ids.FullAccession)
            # -d declares nucleotide input; without it the start of each file
            # is checked as it is read (it costs a few ms per proteome)
            sequence_sample = None if q_dna else SequenceTypeSample()

            q_utf8 = True
            with open_input_fasta(fasta_dir + os.sep + fasta_filename) as fasta_file:
                for i_line, line in enumerate(fasta_file):
                    
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        new_id = "%d_%d" % (i_species, i_seq)
                        assign = line[1:].rstrip()
                        if not assign.isascii():    # (nearly every header is: no call)
                            assign, q_header_utf8 = input_header(assign)
                            q_utf8 = q_utf8 and q_header_utf8
                        if len(assign) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fasta_dir + os.sep + fasta_filename, i_line+1))
                            util.fail()
                        if assign in species_seen_dict[fasta_filename]:
                            species_seen_dict[fasta_filename][assign] += 1
                            # if species_seen_dict[fastaFilename][acc] < 10:
                            #     printer.print(f"ERROR: Duplicated gene names found in '{fastaFilename}' - <{acc}>", style="error")
                            duplicated = True
                        else:
                            species_seen_dict[fasta_filename][assign] = 0
                            name = gene_names.add(assign)
                            if name != gene_ids.clean_name(assign):
                                renamed.append((fasta_filename, assign, name))
                        # acc = f"{acc}_{seen[acc]}"
                        ids_file.write("%s: %s\n" % (new_id, assign))
                        output_fasta.write(">%s\n" % new_id)    
                        i_seq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if sequence_sample is not None:
                            sequence_sample.add(line)
                        output_fasta.write(line)
                output_fasta.write("\n")
            output_fasta.close()
            if not q_utf8:
                not_utf8.append(fasta_filename)
            # Nucleotide input needs -d: the search programs (DIAMOND by
            # default) would otherwise compare it as protein, without error
            if sequence_sample is not None and sequence_sample.type() == "dna":
                q_ok = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fasta_filename)
            if i_seq == 0:
                # (the search programs would fail on it, e.g. "Error detecting input file")
                q_ok = False
                print("ERROR: %s contains no sequences" % fasta_filename)
            i_species += 1
            i_seq = 0
        if not q_ok:
            util.fail()

    if duplicated:
        print()
        printer.print("ERROR: Duplicated gene names found.", style="error")
        for filename in species_seen_dict:
            for assign, acc_count in species_seen_dict[filename].items():
                if acc_count > 0:
                    printer.print(f"{acc_count+1}: {escape(assign)} - [orange3]{escape(filename)}[/orange3]")
        printer.print("Please check the input file and make sure the gene names are unique.\n", style="error")
        util.fail()

    if not_utf8:
        printer.print("\nWARNING: The headers of these files are not UTF-8 text: they have been read as "
                      "Windows-1252 (cp1252) text. Please check the gene names (SequenceIDs.txt), or save the "
                      "files as UTF-8:", style="warning")
        for fasta_filename in not_utf8:
            printer.print(f"    [orange3]{escape(fasta_filename)}[/orange3]")
        print()

    if renamed:
        printer.print("\nWARNING: Some gene names would be the same once the characters a name cannot "
                      "hold (spaces, tabs, : ; ( ) , [ ] = ' \") are replaced by _. These genes have been "
                      "given distinct names:", style="warning")
        for fasta_filename, assign, name in renamed:
            printer.print(f"    {escape(assign)} -> {escape(name)}  [orange3]{escape(fasta_filename)}[/orange3]")
        print()

    if len(original_fasta_filenames) > 0: output_fasta.close()
    species_info_obj.species_to_use = species_info_obj.species_to_use + new_species_ids
    species_info_obj.n_sp_all = max(species_info_obj.species_to_use) + 1      # will be one of the new species
    
    return species_info_obj

