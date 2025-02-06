# -*- coding: utf-8 -*-
import sys
import gzip
import string
import os.path
import glob
from . import util, files 


# count = 0 

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


def ProcessesNewFasta(fastaDir, q_dna, speciesInfoObj_prev = None, speciesToUse_prev_names=[]):
    """
    Process fasta files and return a Directory object with all paths completed.
    """
    fastaExtensions = {"fa", "faa", "fasta", "fas", "pep"}
    # Check files present
    qOk = True
    if not os.path.exists(fastaDir):
        print("\nDirectory does not exist: %s" % fastaDir)
        util.Fail()
    files_in_directory = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
    originalFastaFilenames = []
    excludedFiles = []

    for f in files_in_directory:
        if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions and not f.startswith("._"):
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

    if len(originalFastaFilenames) == 0:
        print("\nNo fasta files found in supplied directory: %s" % fastaDir)
        util.Fail()

    if speciesInfoObj_prev == None:
        # Then this is a new, clean analysis 
        speciesInfoObj = util.SpeciesInfo()
    else:
        speciesInfoObj = speciesInfoObj_prev

    iSeq = 0
    iSpecies = 0
    # If it's a previous analysis:
    if len(speciesToUse_prev_names) != 0:
        with open(files.FileHandler.GetSpeciesIDsFN(), 'r') as infile:
            for line in infile: pass
        if line.startswith("#"): line = line[1:]
        iSpecies = int(line.split(":")[0]) + 1
    speciesInfoObj.iFirstNewSpecies = iSpecies
    newSpeciesIDs = []

    with open(files.FileHandler.GetSequenceIDsFN(), 'a') as idsFile, open(files.FileHandler.GetSpeciesIDsFN(), 'a') as speciesFile:
        for fastaFilename in originalFastaFilenames:
            newSpeciesIDs.append(iSpecies)
            outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
            fastaFilename = fastaFilename.rstrip()
            speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
            baseFilename, extension = os.path.splitext(fastaFilename)
            mLinesToCheck = 100
            qHasAA = False
            with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                for iLine, line in enumerate(fastaFile):
                    if line.isspace(): continue
                    if len(line) > 0 and line[0] == ">":
                        newID = "%d_%d" % (iSpecies, iSeq)
                        acc = line[1:].rstrip()
                        if len(acc) == 0:
                            print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                            util.Fail()
                        idsFile.write("%s: %s\n" % (newID, acc))
                        outputFasta.write(">%s\n" % newID)    
                        iSeq += 1
                    else:
                        line = line.upper()    # allow lowercase letters in sequences
                        if not qHasAA and (iLine < mLinesToCheck):
#                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                            qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                        outputFasta.write(line)
                outputFasta.write("\n")
            if (not qHasAA) and (not q_dna):
                qOk = False
                print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
            iSpecies += 1
            iSeq = 0
            outputFasta.close()
        if not qOk:
            util.Fail()

    if len(originalFastaFilenames) > 0: outputFasta.close()
    speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
    speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    
    return speciesInfoObj