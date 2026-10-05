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
from . import util, files 
from ..utils.util import printer

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
    with open(fn) as infile:
        for line in infile:
            if not line.startswith(">"):
                sample.add(line)
                if sample.full:
                    break
    return sample.type()


class FastaWriter(object):
    def __init__(self, sourceFastaFilename, qUseOnlySecondPart=False, qGlob=False, qFirstWord=False, qWholeLine=False):
        self.SeqLists = dict()
        qFirst = True
        accession = ""
        sequence = ""
        files = glob.glob(sourceFastaFilename) if qGlob else [sourceFastaFilename]
        for fn in files:
            with (gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn, 'r')) as fastaFile:
                for line in fastaFile:
                    if line[0] == ">":
                        # deal with old sequence
                        if not qFirst:
                            self.SeqLists[accession] = sequence
                            sequence = ""
                        qFirst = False
                        # get id for new sequence
                        if qWholeLine:
                            accession = line[1:].rstrip()
                        else:
                            accession = line[1:].rstrip().split()[0]
                        if qUseOnlySecondPart:
                            try:
                                accession = accession.split("|")[1]
                            except:
                                sys.stderr(line + "\n")
                        elif qFirstWord:
                            accession = accession.split()[0]
                    else:
                        sequence += line
            try:
                self.SeqLists[accession] = sequence
            except:
                sys.stderr(accession + "\n")
    
    def WriteSeqsToFasta(self, seqs, outFilename):
        with open(outFilename, 'w') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % seq)
                    outFile.write(self.SeqLists[seq])
                else:
                    sys.stderr.write("ERROR: %s not found\n" % seq)
    
    def AppendToStockholm(self, msa_id, outFilename, seqs=None):
        """
        Append the msa to a stockholm format file. File must already be aligned
        Args:
            msa_id - the ID to use for the accession line, "#=GF AC" 
            outFilename - the file to write to
            seqs - Only write selected sequences o/w write all
        """
        l = [len(s) for s in self.SeqLists.values()]
        if min(l) != max(l):
            raise Exception("Sequences must be aligned")
        with open(outFilename, 'a') as outFile:
            outFile.write("# STOCKHOLM 1.0\n#=GF AC %s\n" % msa_id)
            seqs = self.SeqLists.keys() if seqs is None else seqs
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write("%s " % seq)
                    outFile.write(self.SeqLists[seq].replace("\n", "") + "\n")
                else:
                    sys.stderr("ERROR: %s not found\n" % seq)
            outFile.write("//\n")
            
    def Print(self, seqs):
        if type(seqs[0]) is not str:
            seqs = [s.ToString() for s in seqs]
        for seq in seqs:
            if seq in self.SeqLists:
                sys.stdout.write(">%s\n" % seq)
                sys.stdout.write(self.SeqLists[seq])
            else:
                sys.stderr.write("ERROR: %s not found\n" % seq)

                    
    def WriteSeqsToFasta_withNewAccessions(self, seqs, outFilename, idDict):
        with open(outFilename, 'w') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % idDict[seq])
                    outFile.write(self.SeqLists[seq])
                else:
                    sys.stderr.write(seq + "\n")


def ProcessesNewFasta(
        fastaDir, 
        q_dna, 
        speciesInfoObj_prev = None, 
        speciesToUse_prev_names=[],
        species_id_fn="",
        sequence_id_fn="",
    ):
    """
    Process fasta files and return a Directory object with all paths completed.
    """

    fastaExtensions = FASTA_EXTENSIONS
    # Check files present
    qOk = True
    if not os.path.exists(fastaDir):
        print("\nDirectory does not exist: %s" % fastaDir)
        util.Fail()
    files_in_directory = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
    originalFastaFilenames = []
    excludedFiles = []

    for f in files_in_directory:
        if is_input_fasta_name(f):
            originalFastaFilenames.append(f)
        else:
            excludedFiles.append(f)

    if len(excludedFiles) != 0:
        print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
        for f in excludedFiles:
            print(f)
        print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fastaExtensions)))
    
    speciesToUse_prev_names = set(speciesToUse_prev_names)
    if len(originalFastaFilenames) + len(speciesToUse_prev_names) < 2:
        print("ERROR: At least two species are required")
        util.Fail()

    if any([fn in speciesToUse_prev_names for fn in originalFastaFilenames]):
        print("ERROR: Attempted to add a second copy of a previously included species:")
        for fn in originalFastaFilenames:
            if fn in speciesToUse_prev_names: print(fn)
        print("")
        util.Fail()

    # A species is named by its file name without the extension (as in the
    # species tree and the results): "A.fa" and "A.fna" would both be "A".
    files_of_name = {}
    for fn in sorted(speciesToUse_prev_names) + sorted(originalFastaFilenames):
        files_of_name.setdefault(fn.rsplit(".", 1)[0], []).append(fn)
    same_name = {name: fns for name, fns in files_of_name.items() if len(fns) > 1}
    if same_name:
        print("ERROR: Species must have different names. These files would give species of the same name "
              "(the file name without its extension):")
        for name, fns in sorted(same_name.items()):
            print("  %s: %s" % (name, ", ".join(fns)))
        if speciesToUse_prev_names:
            print("(including the species already in the analysis)")
        print("")
        util.Fail()

    speciesToUse_prev_names = sorted([*speciesToUse_prev_names])
    originalFastaFilenames = sorted([*originalFastaFilenames])

    if len(originalFastaFilenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fastaDir)
        util.Fail()

    if speciesInfoObj_prev == None:
        speciesInfoObj = util.SpeciesInfo()
    else:
        speciesInfoObj = speciesInfoObj_prev

    if not species_id_fn:
        if files.FileHandler.wd_current:
            species_id_fn = os.path.join(files.FileHandler.wd_current, "SpeciesIDs.txt")
        else:
            species_id_fn = files.FileHandler.GetSpeciesIDsFN()

    if not sequence_id_fn:
        if files.FileHandler.wd_current:
            sequence_id_fn = os.path.join(files.FileHandler.wd_current, "SequenceIDs.txt")
        else:
            sequence_id_fn = files.FileHandler.GetSequenceIDsFN()

    iSeq = 0
    iSpecies = 0
    # If it's a previous analysis:
    if len(speciesToUse_prev_names) != 0:
        with open(species_id_fn, 'r') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    speciesInfoObj.iFirstNewSpecies = iSpecies
    newSpeciesIDs = []
    duplicated = False

    species_seen_dict = {}
    renamed = []      # (file, header, name): headers that clean to a name already used

    with open(sequence_id_fn, 'a') as idsFile, open(species_id_fn, 'a') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)

            species_seen_dict[fastaFilename] = {}
            gene_names = util.GeneNames()     # the names the genes get (util.FullAccession)
            # -d declares nucleotide input; without it the start of each file
            # is checked as it is read (it costs a few ms per proteome)
            sequence_sample = None if q_dna else SequenceTypeSample()

            with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                for iLine, line in enumerate(fastaFile):
                    
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        acc = line[1:].rstrip()
                        if len(acc) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                            util.Fail()
                        if acc in species_seen_dict[fastaFilename]:
                            species_seen_dict[fastaFilename][acc] += 1
                            # if species_seen_dict[fastaFilename][acc] < 10:
                            #     printer.print(f"ERROR: Duplicated gene names found in '{fastaFilename}' - <{acc}>", style="error")
                            duplicated = True
                        else:
                            species_seen_dict[fastaFilename][acc] = 0
                            name = gene_names.add(acc)
                            if name != util.CleanAccession(acc):
                                renamed.append((fastaFilename, acc, name))
                        # acc = f"{acc}_{seen[acc]}"
                        idsFile.write("%s: %s\n" % (newID, acc))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if sequence_sample is not None:
                            sequence_sample.add(line)
                        outputFasta.write(line)
                outputFasta.write("\n")
            outputFasta.close()
            # Nucleotide input needs -d: the search programs (DIAMOND by
            # default) would otherwise compare it as protein, without error
            if sequence_sample is not None and sequence_sample.type() == "dna":
                qOk = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
            iSpecies += 1
            iSeq = 0
        if not qOk:
            util.Fail()

    if duplicated:
        print()
        printer.print("ERROR: Duplicated gene names found.", style="error")
        for filename in species_seen_dict:
            for acc, acc_count in species_seen_dict[filename].items():
                if acc_count > 0:
                    printer.print(f"{acc_count+1}: {acc} - [orange3]{filename}[/orange3]")
        printer.print("Please check the input file and make sure the gene names are unique.\n", style="error")
        util.Fail()

    if renamed:
        printer.print("\nWARNING: Some gene names would be the same once the characters that output "
                      "files cannot hold (: , ( ) ; [ ] =) are replaced by _. These genes have been "
                      "given distinct names:", style="warning")
        for fastaFilename, acc, name in renamed:
            printer.print(f"    {escape(acc)} -> {escape(name)}  [orange3]{escape(fastaFilename)}[/orange3]")
        print()

    if len(originalFastaFilenames) > 0: outputFasta.close()
    speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
    speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    
    return speciesInfoObj

