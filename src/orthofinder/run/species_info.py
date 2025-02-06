from ..utils import util, files
import os
import csv

try:
    from rich import print
except ImportError:
    ...


def IDsFileOK(filename):
    """
    It is best to detect any issues with input files at start, perform all required checks here
    """
    with open(filename, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if len(line) == 0:
                continue
            tokens = line.split(": ", 1)
            if len(tokens) != 2 or len(tokens[1]) == 0:
                return False, line
    return True, None


# 0
def ProcessPreviousFiles(workingDir_list, qDoubleBlast, check_blast=True):
    """Checks for:
    workingDir should be the WorkingDirectory containing Blast*.txt files

    SpeciesIDs.txt
    Species*.fa
    Blast*txt
    SequenceIDs.txt

    Checks which species should be included

    """
    # check BLAST results directory exists
    if not os.path.exists(workingDir_list[0]):
        err_text = (
            "ERROR: Previous/Pre-calculated BLAST results directory does not exist: %s\n"
            % workingDir_list[0]
        )
        files.FileHandler.LogFailAndExit(err_text)

    speciesInfo = util.SpeciesInfo()
    if not os.path.exists(files.FileHandler.GetSpeciesIDsFN()):
        err_text = (
            "ERROR: %s file must be provided if using previously calculated BLAST results"
            % files.FileHandler.GetSpeciesIDsFN()
        )
        files.FileHandler.LogFailAndExit(err_text)
    file_ok, err_line = IDsFileOK(files.FileHandler.GetSpeciesIDsFN())
    if not file_ok:
        files.FileHandler.LogFailAndExit(
            "ERROR: %s file contains a blank accession. Line:\n %s"
            % (files.FileHandler.GetSpeciesIDsFN(), err_line)
        )
    speciesInfo.speciesToUse, speciesInfo.nSpAll, speciesToUse_names = (
        util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN())
    )

    # check fasta files are present
    previousFastaFiles = files.FileHandler.GetSortedSpeciesFastaFiles()
    if len(previousFastaFiles) == 0:
        err_text = (
            "ERROR: No processed fasta files in the supplied previous working directories:\n"
            + "\n".join(workingDir_list)
            + "\n"
        )
        files.FileHandler.LogFailAndExit(err_text)
    tokens = previousFastaFiles[-1][:-3].split("Species")
    lastFastaNumberString = tokens[-1]
    iLastFasta = 0
    nFasta = len(previousFastaFiles)
    try:
        iLastFasta = int(lastFastaNumberString)
    except:
        files.FileHandler.LogFailAndExit(
            "ERROR: Filenames for processed fasta files are incorrect: %s\n"
            % previousFastaFiles[-1]
        )
    if nFasta != iLastFasta + 1:
        files.FileHandler.LogFailAndExit(
            "ERROR: Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n"
            % (lastFastaNumberString, len(previousFastaFiles))
        )

    # check BLAST files
    if check_blast:
        blast_fns_triangular = [
            files.FileHandler.GetBlastResultsFN(
                iSpecies, jSpecies, raise_exception=False
            )
            for iSpecies in speciesInfo.speciesToUse
            for jSpecies in speciesInfo.speciesToUse
            if jSpecies >= iSpecies
        ]
        have_triangular = [
            (os.path.exists(fn) or os.path.exists(fn + ".gz"))
            for fn in blast_fns_triangular
        ]
        for qHave, fn in zip(have_triangular, blast_fns_triangular):
            if not qHave:
                print("Required BLAST results file not found: %s" % fn)

        if qDoubleBlast:
            blast_fns_remainder = [
                files.FileHandler.GetBlastResultsFN(
                    iSpecies, jSpecies, raise_exception=False
                )
                for iSpecies in speciesInfo.speciesToUse
                for jSpecies in speciesInfo.speciesToUse
                if jSpecies < iSpecies
            ]
            have_remainder = [
                (os.path.exists(fn) or os.path.exists(fn + ".gz"))
                for fn in blast_fns_remainder
            ]
            if not (all(have_triangular) and all(have_remainder)):
                for qHave, fn in zip(have_remainder, blast_fns_remainder):
                    if not qHave:
                        print("Required BLAST results file not found: %s" % fn)
                if not all(have_triangular):
                    files.FileHandler.LogFailAndExit()
                else:
                    # would be able to do it using just one-way blast
                    files.FileHandler.LogFailAndExit(
                        "ERROR: Required BLAST results files are present for using the one-way sequence search option (default) but not the double BLAST search ('-d' option)"
                    )
        else:
            if not all(have_triangular):
                files.FileHandler.LogFailAndExit()

    # check SequenceIDs.txt and SpeciesIDs.txt files are present
    if not os.path.exists(files.FileHandler.GetSequenceIDsFN()):
        files.FileHandler.LogFailAndExit(
            "ERROR: %s file must be provided if using previous calculated BLAST results"
            % files.FileHandler.GetSequenceIDsFN()
        )

    file_ok, err_line = IDsFileOK(files.FileHandler.GetSequenceIDsFN())
    if not file_ok:
        files.FileHandler.LogFailAndExit(
            "ERROR: %s file contains a blank accession. Line:\n %s"
            % (files.FileHandler.GetSequenceIDsFN(), err_line)
        )
    return speciesInfo, speciesToUse_names


def GetXMLSpeciesInfo(seqsInfoObj, options):
    # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
    util.PrintUnderline("Reading species information file")
    # do this now so that we can alert user to any errors prior to running the algorithm
    speciesXML = [[] for i_ in seqsInfoObj.speciesToUse]
    speciesNamesDict = SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
    speciesRevDict = {v: k for k, v in speciesNamesDict.items()}
    userFastaFilenames = [
        os.path.split(speciesNamesDict[i])[1] for i in seqsInfoObj.speciesToUse
    ]
    with open(options.speciesXMLInfoFN, "r") as speciesInfoFile:
        reader = csv.reader(speciesInfoFile, delimiter="\t")
        for iLine, line in enumerate(reader):
            if len(line) != 5:
                # allow for an extra empty line at the end
                if len(line) == 0 and iLine == len(userFastaFilenames):
                    continue
                print("ERROR")
                print(
                    "Species information file %s line %d is incorrectly formatted."
                    % (options.speciesXMLInfoFN, iLine + 1)
                )
                print("File should be contain one line per species")
                print("Each line should contain 5 tab-delimited fields:")
                print(
                    "  fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseFastaFilename"
                )
                print("See README file for more information.")
                util.Fail()
            (
                fastaFilename,
                speciesName,
                NCBITaxID,
                sourceDatabaseName,
                databaseVersionFastaFile,
            ) = line
            try:
                iSpecies = speciesRevDict[os.path.splitext(fastaFilename)[0]]
            except KeyError:
                print(
                    "Skipping %s from line %d as it is not being used in this analysis"
                    % (fastaFilename, iLine + 1)
                )
                continue
            speciesXML[seqsInfoObj.speciesToUse.index(iSpecies)] = line
    # check information has been provided for all species
    speciesMissing = False
    for iPos, iSpecies in enumerate(seqsInfoObj.speciesToUse):
        if speciesXML[iPos] == []:
            if not speciesMissing:
                print("ERROR")
                print(
                    "Species information file %s does not contain information for all species."
                    % options.speciesXMLInfoFN
                )
                print("Information is missing for:")
                speciesMissing = True
            print(speciesNamesDict[iSpecies])
    if speciesMissing:
        util.Fail()
    return speciesXML


def SpeciesNameDict(speciesIDsFN):
    # type: (str) -> Dict[int, str]
    speciesNamesDict = dict()
    with open(speciesIDsFN, "r") as speciesNamesFile:
        for line in speciesNamesFile:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            if not line:
                continue
            short, full = line.split(": ")
            speciesNamesDict[int(short)] = full.rsplit(".", 1)[0]
    return speciesNamesDict


def SpeciesSequenceNameDict(sequenceIDsFN):
    try:
        seqIDsEx = util.FirstWordExtractor(sequenceIDsFN)
    except RuntimeError as error:
        print(str(error))
        if str(error).startswith("ERROR"): 
            files.FileHandler.LogFailAndExit()
        else:
            print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup")
            print("more concisely but these were not unique. The full accession line will be used instead.\n")
            seqIDsEx = util.FullAccession(sequenceIDsFN)
    return seqIDsEx.GetIDToNameDict()