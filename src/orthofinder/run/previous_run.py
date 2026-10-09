from ..utils import util, file_handler, file_io
import os
from collections import Counter

try:
    from rich import print
except ImportError:
    ...

from ..phylo import tree

def ids_file_ok(filename):
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


# # 0
# def process_previous_files(workingDir_list, qDoubleBlast, check_blast=True):
#     """Checks for:
#     workingDir should be the WorkingDirectory containing Blast*.txt files

#     SpeciesIDs.txt
#     Species*.fa
#     Blast*txt
#     SequenceIDs.txt

#     Checks which species should be included

#     """
#     # check BLAST results directory exists
#     if not os.path.exists(workingDir_list[0]):
#         err_text = (
#             "ERROR: Previous/Pre-calculated BLAST results directory does not exist: %s\n"
#             % workingDir_list[0]
#         )
#         files.FileHandler.log_fail_and_exit(err_text)

#     speciesInfo = util.SpeciesInfo()
#     if not os.path.exists(files.FileHandler.get_species_ids_fn()):
#         err_text = (
#             "ERROR: %s file must be provided if using previously calculated BLAST results"
#             % files.FileHandler.get_species_ids_fn()
#         )
#         files.FileHandler.log_fail_and_exit(err_text)
#     file_ok, err_line = ids_file_ok(files.FileHandler.get_species_ids_fn())
#     if not file_ok:
#         files.FileHandler.log_fail_and_exit(
#             "ERROR: %s file contains a blank accession. Line:\n %s"
#             % (files.FileHandler.get_species_ids_fn(), err_line)
#         )
#     print(f"[DEBUG] species_file_name: {files.FileHandler.get_species_ids_fn()}")
#     speciesInfo.speciesToUse, speciesInfo.nSpAll, speciesToUse_names = (
#         util.get_species_to_use(files.FileHandler.get_species_ids_fn())
#     )

#     # check fasta files are present
#     previousFastaFiles = files.FileHandler.get_sorted_species_fasta_files()
#     if len(previousFastaFiles) == 0:
#         err_text = (
#             "ERROR: No processed fasta files in the supplied previous working directories:\n"
#             + "\n".join(workingDir_list)
#             + "\n"
#         )
#         files.FileHandler.log_fail_and_exit(err_text)
#     tokens = previousFastaFiles[-1][:-3].split("Species")
#     lastFastaNumberString = tokens[-1]
#     iLastFasta = 0
#     nFasta = len(previousFastaFiles)
#     try:
#         iLastFasta = int(lastFastaNumberString)
#     except:
#         files.FileHandler.log_fail_and_exit(
#             "ERROR: Filenames for processed fasta files are incorrect: %s\n"
#             % previousFastaFiles[-1]
#         )
#     if nFasta != iLastFasta + 1:
#         files.FileHandler.log_fail_and_exit(
#             "ERROR: Not all expected fasta files are present. Index of last fasta file is %s but found %d fasta files.\n"
#             % (lastFastaNumberString, len(previousFastaFiles))
#         )

#     # check BLAST files
#     if check_blast:
#         blast_fns_triangular = [
#             files.FileHandler.get_blast_results_fn(
#                 iSpecies, jSpecies, raise_exception=False
#             )
#             for iSpecies in speciesInfo.speciesToUse
#             for jSpecies in speciesInfo.speciesToUse
#             if jSpecies >= iSpecies
#         ]
#         have_triangular = [
#             (os.path.exists(fn) or os.path.exists(fn + ".gz"))
#             for fn in blast_fns_triangular
#         ]
#         for qHave, fn in zip(have_triangular, blast_fns_triangular):
#             if not qHave:
#                 print("Required BLAST results file not found: %s" % fn)

#         if qDoubleBlast:
#             blast_fns_remainder = [
#                 files.FileHandler.get_blast_results_fn(
#                     iSpecies, jSpecies, raise_exception=False
#                 )
#                 for iSpecies in speciesInfo.speciesToUse
#                 for jSpecies in speciesInfo.speciesToUse
#                 if jSpecies < iSpecies
#             ]
#             have_remainder = [
#                 (os.path.exists(fn) or os.path.exists(fn + ".gz"))
#                 for fn in blast_fns_remainder
#             ]
#             if not (all(have_triangular) and all(have_remainder)):
#                 for qHave, fn in zip(have_remainder, blast_fns_remainder):
#                     if not qHave:
#                         print("Required BLAST results file not found: %s" % fn)
#                 if not all(have_triangular):
#                     files.FileHandler.log_fail_and_exit()
#                 else:
#                     # would be able to do it using just one-way blast
#                     files.FileHandler.log_fail_and_exit(
#                         "ERROR: Required BLAST results files are present for using the one-way sequence search option (default) but not the double BLAST search ('-d' option)"
#                     )
#         else:
#             if not all(have_triangular):
#                 files.FileHandler.log_fail_and_exit()

#     # check SequenceIDs.txt and SpeciesIDs.txt files are present
#     if not os.path.exists(files.FileHandler.get_sequence_ids_fn()):
#         files.FileHandler.log_fail_and_exit(
#             "ERROR: %s file must be provided if using previous calculated BLAST results"
#             % files.FileHandler.get_sequence_ids_fn()
#         )

#     file_ok, err_line = ids_file_ok(files.FileHandler.get_sequence_ids_fn())
#     if not file_ok:
#         files.FileHandler.log_fail_and_exit(
#             "ERROR: %s file contains a blank accession. Line:\n %s"
#             % (files.FileHandler.get_sequence_ids_fn(), err_line)
#         )
#     return speciesInfo, speciesToUse_names

def get_blast_results_to_create(species_to_use, q_double_blast, only=None):
    """
    The required search results that are missing or incomplete (their search
    was interrupted), as (filename without .gz, incomplete file or None).
    only: if given, just these filenames (without .gz) are checked.

    No file is read: checkpoint.txt in the results directory tells whether
    the searches were interrupted, and if so which ones finished, with the
    size of their results (file_io.search_checkpoint).
    """
    required = [
        file_handler.handler.get_blast_results_fn(
            i_species, j_species, raise_exception=False
        )
        for i_species in species_to_use
        for j_species in species_to_use
        if q_double_blast or j_species >= i_species
    ]
    if only is not None:
        only = set(only)
        required = [fn for fn in required if fn in only]

    return results_to_create(required)


def results_to_create(required):
    """
    Of the search results files required (names without .gz), those missing
    or incomplete, as (filename, incomplete file or None): see
    get_blast_results_to_create.
    """
    checkpoints = {}
    def complete(path):
        d, name = os.path.split(path)
        if d not in checkpoints:
            checkpoints[d] = file_io.search_checkpoint(d)
        interrupted, done = checkpoints[d]
        return not interrupted or done.get(name) == os.path.getsize(path)

    to_create = []
    for fn in required:
        path = fn + ".gz" if os.path.exists(fn + ".gz") else fn if os.path.exists(fn) else None
        if path is None:
            to_create.append((fn, None))
        elif not complete(path):
            to_create.append((fn, path))
    return to_create


def process_previous_files(working_dir_list, q_double_blast, check_blast=True):
    if not working_dir_list:
        file_handler.handler.log_fail_and_exit(
            "ERROR: No previous WorkingDirectory paths supplied to ProcessPreviousFiles"
        )

    primary_wd = working_dir_list[0]

    if not os.path.exists(primary_wd):
        err_text = (
            "ERROR: Previous/Pre-calculated BLAST results directory does not exist: %s\n"
            % primary_wd
        )
        file_handler.handler.log_fail_and_exit(err_text)

    # ---- SpeciesIDs.txt: use the primary working directory explicitly ----
    species_ids_fn = os.path.join(primary_wd, "SpeciesIDs.txt")
    if not os.path.exists(species_ids_fn):
        err_text = (
            "ERROR: %s file must be provided if using previously calculated BLAST results"
            % species_ids_fn
        )
        file_handler.handler.log_fail_and_exit(err_text)

    file_ok, err_line = ids_file_ok(species_ids_fn)
    if not file_ok:
        file_handler.handler.log_fail_and_exit(
            "ERROR: %s file contains a blank accession. Line:\n %s"
            % (species_ids_fn, err_line)
        )

    species_info = util.SpeciesInfo()
    species_info.species_to_use, species_info.n_sp_all, species_to_use_names = (
        util.get_species_to_use(species_ids_fn)
    )

    # ---- Check fasta files are present (still uses FileHandler for the list) ----
    previous_fasta_files = file_handler.handler.get_sorted_species_fasta_files()
    if len(previous_fasta_files) == 0:
        err_text = (
            "ERROR: No processed fasta files in the supplied previous working directories:\n"
            + "\n".join(working_dir_list)
            + "\n"
        )
        file_handler.handler.log_fail_and_exit(err_text)

    tokens = previous_fasta_files[-1][:-3].split("Species")
    last_fasta_number_string = tokens[-1]
    n_fasta = len(previous_fasta_files)
    try:
        i_last_fasta = int(last_fasta_number_string)
    except Exception:
        file_handler.handler.log_fail_and_exit(
            "ERROR: Filenames for processed fasta files are incorrect: %s\n"
            % previous_fasta_files[-1]
        )
    if n_fasta != i_last_fasta + 1:
        file_handler.handler.log_fail_and_exit(
            "ERROR: Not all expected fasta files are present. Index of last fasta file is %s "
            "but found %d fasta files.\n"
            % (last_fasta_number_string, len(previous_fasta_files))
        )

    # ---- Check BLAST files as before ----
    if check_blast:
        blast_fns_triangular = [
            file_handler.handler.get_blast_results_fn(
                i_species, j_species, raise_exception=False
            )
            for i_species in species_info.species_to_use
            for j_species in species_info.species_to_use
            if j_species >= i_species
        ]
        have_triangular = [
            (os.path.exists(fn) or os.path.exists(fn + ".gz"))
            for fn in blast_fns_triangular
        ]
        for q_have, fn in zip(have_triangular, blast_fns_triangular):
            if not q_have:
                print("Required BLAST results file not found: %s" % fn)

        if q_double_blast:
            blast_fns_remainder = [
                file_handler.handler.get_blast_results_fn(
                    i_species, j_species, raise_exception=False
                )
                for i_species in species_info.species_to_use
                for j_species in species_info.species_to_use
                if j_species < i_species
            ]
            have_remainder = [
                (os.path.exists(fn) or os.path.exists(fn + ".gz"))
                for fn in blast_fns_remainder
            ]
            if not (all(have_triangular) and all(have_remainder)):
                for q_have, fn in zip(have_remainder, blast_fns_remainder):
                    if not q_have:
                        print("Required BLAST results file not found: %s" % fn)
                if not all(have_triangular):
                    file_handler.handler.log_fail_and_exit()
                else:
                    file_handler.handler.log_fail_and_exit(
                        "ERROR: Required BLAST results files are present for using the one-way "
                        "sequence search option (default) but not the double BLAST search ('-d' option)"
                    )
        else:
            if not all(have_triangular):
                file_handler.handler.log_fail_and_exit()

    # ---- SequenceIDs.txt: again, use the primary working directory explicitly ----
    sequence_ids_fn = os.path.join(primary_wd, "SequenceIDs.txt")
    if not os.path.exists(sequence_ids_fn):
        file_handler.handler.log_fail_and_exit(
            "ERROR: %s file must be provided if using previous calculated BLAST results"
            % sequence_ids_fn
        )

    file_ok, err_line = ids_file_ok(sequence_ids_fn)
    if not file_ok:
        file_handler.handler.log_fail_and_exit(
            "ERROR: %s file contains a blank accession. Line:\n %s"
            % (sequence_ids_fn, err_line)
        )

    return species_info, species_to_use_names


def check_user_species_tree(species_tree_fn, exp_species):
    # File exists
    if not os.path.exists(species_tree_fn):
        print(("Species tree file does not exist: %s" % species_tree_fn))
        util.fail()
    # Species in tree are unique
    try:
        t = tree.Tree(species_tree_fn, format=1)
    except Exception as e:
        print("\nERROR: Incorrectly formated user-supplied species tree")
        print(str(e))
        util.fail()
    act_species = (t.get_leaf_names())
    c = Counter(act_species)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names in species tree are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.fail()
    # All required species are present
    act_species = set(act_species)
    ok = True
    for sp in exp_species:
        if sp not in act_species:
            print(("ERROR: '%s' is missing from species tree" % sp))
            ok = False
    # expected species are unique
    c = Counter(exp_species)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.fail()
    exp_species = set(exp_species)
    for sp in act_species:
        if sp not in exp_species:
            print(("ERROR: Additional species '%s' in species tree" % sp))
            ok = False
    if not ok: util.fail()
    # Tree is rooted
    if len(t.get_children()) != 2:
        print("ERROR: Species tree is not rooted")
        util.fail()




# def GetXMLSpeciesInfo(seqsInfoObj, options):
#     # speciesInfo:  name, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
#     util.print_underline("Reading species information file")
#     # do this now so that we can alert user to any errors prior to running the algorithm
#     speciesXML = [[] for i_ in seqsInfoObj.speciesToUse]
#     speciesNamesDict = species_name_dict(files.FileHandler.get_species_ids_fn())
#     speciesRevDict = {v: k for k, v in speciesNamesDict.items()}
#     userFastaFilenames = [
#         os.path.split(speciesNamesDict[i])[1] for i in seqsInfoObj.speciesToUse
#     ]
#     with open(options.speciesXMLInfoFN, "r") as speciesInfoFile:
#         reader = csv.reader(speciesInfoFile, delimiter="\t")
#         for iLine, line in enumerate(reader):
#             if len(line) != 5:
#                 # allow for an extra empty line at the end
#                 if len(line) == 0 and iLine == len(userFastaFilenames):
#                     continue
#                 print("ERROR")
#                 print(
#                     "Species information file %s line %d is incorrectly formatted."
#                     % (options.speciesXMLInfoFN, iLine + 1)
#                 )
#                 print("File should be contain one line per species")
#                 print("Each line should contain 5 tab-delimited fields:")
#                 print(
#                     "  fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseFastaFilename"
#                 )
#                 print("See README file for more information.")
#                 util.Fail()
#             (
#                 fastaFilename,
#                 speciesName,
#                 NCBITaxID,
#                 sourceDatabaseName,
#                 databaseVersionFastaFile,
#             ) = line
#             try:
#                 iSpecies = speciesRevDict[os.path.splitext(fastaFilename)[0]]
#             except KeyError:
#                 print(
#                     "Skipping %s from line %d as it is not being used in this analysis"
#                     % (fastaFilename, iLine + 1)
#                 )
#                 continue
#             speciesXML[seqsInfoObj.speciesToUse.index(iSpecies)] = line
#     # check information has been provided for all species
#     speciesMissing = False
#     for iPos, iSpecies in enumerate(seqsInfoObj.speciesToUse):
#         if speciesXML[iPos] == []:
#             if not speciesMissing:
#                 print("ERROR")
#                 print(
#                     "Species information file %s does not contain information for all species."
#                     % options.speciesXMLInfoFN
#                 )
#                 print("Information is missing for:")
#                 speciesMissing = True
#             print(speciesNamesDict[iSpecies])
#     if speciesMissing:
#         util.Fail()
#     return speciesXML


def species_name_dict(species_ids_fn):
    # type: (str) -> Dict[int, str]
    species_names_dict = dict()
    with open(species_ids_fn, "r") as species_names_file:
        for line in species_names_file:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            if not line:
                continue
            short, full = line.split(": ")
            species_names_dict[int(short)] = full.rsplit(".", 1)[0]
    return species_names_dict


# def SpeciesSequenceNameDict(sequenceIDsFN):
#     try:
#         seqIDsEx = util.FirstWordExtractor(sequenceIDsFN)
#     except RuntimeError as error:
#         print(str(error))
#         if str(error).startswith("ERROR"): 
#             files.FileHandler.log_fail_and_exit()
#         else:
#             print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup")
#             print("more concisely but these were not unique. The full accession line will be used instead.\n")
#             seqIDsEx = util.FullAccession(sequenceIDsFN)
#     return seqIDsEx.get_id_to_name_dict()
