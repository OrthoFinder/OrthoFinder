from ..utils import util, file_handler, file_io
import os
from collections import Counter

try:
    from rich import print
except ImportError:
    ...

from ..phylo import tree
from ..utils import species_info as utils_species_info
import glob
from ..utils import run_state
from ..utils.file_handler import ASSIGN_FIRST_OG_LOG_KEY, EARLIER_WD_KEYS, PROGRAM_LOG_KEYS, SPECIES_TREE_LOG_KEY, handler

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

#     speciesInfo = species_info.SpeciesInfo()
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
#         species_info.get_species_to_use(files.FileHandler.get_species_ids_fn())
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

    species_info = utils_species_info.SpeciesInfo()
    species_info.species_to_use, species_info.n_sp_all, species_to_use_names = (
        utils_species_info.get_species_to_use(species_ids_fn)
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
#         seqIDsEx = gene_ids.FirstWordExtractor(sequenceIDsFN)
#     except RuntimeError as error:
#         print(str(error))
#         if str(error).startswith("ERROR"): 
#             files.FileHandler.log_fail_and_exit()
#         else:
#             print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup")
#             print("more concisely but these were not unique. The full accession line will be used instead.\n")
#             seqIDsEx = gene_ids.FullAccession(sequenceIDsFN)
#     return seqIDsEx.get_id_to_name_dict()


def state_from_log(continuation_dir):
    """
    What a restart continues with, from the Log.txt of the results directory
    it is given (or of the one holding that WorkingDirectory), in the form of
    a run_state state, for results directories without run_state.json (older
    versions, or the file removed): the programs used, the species tree
    (FN_SpeciesTree) and, for --assign, the first orthogroup whose gene tree
    was not inferred to place the new species. The last record of each counts.
    Read before the restart writes its own Log.txt. {} if there is no Log.txt.
    """
    d = os.path.abspath(continuation_dir)
    log_fn = next((os.path.join(r, "Log.txt") for r in (d, os.path.dirname(d))
                   if os.path.exists(os.path.join(r, "Log.txt"))), None)
    if log_fn is None:
        return {}
    state = {}
    with open(log_fn, errors="replace") as infile:
        for line in infile:
            line = line.rstrip("\r\n")
            for key, name in PROGRAM_LOG_KEYS:
                if line.startswith(key):
                    state.setdefault("settings", {})[name] = line[len(key):].strip()
            if line.startswith(SPECIES_TREE_LOG_KEY):
                state.setdefault("files", {})["species_tree"] = PreviousFilesLocatorNew._in_this_results_dir(
                    line[len(SPECIES_TREE_LOG_KEY):], log_fn, 2)
            if line.startswith(ASSIGN_FIRST_OG_LOG_KEY):
                state.setdefault("assign", {})["i_og_restart"] = int(line[len(ASSIGN_FIRST_OG_LOG_KEY):])
    return state


class Unprocessable(Exception):
    pass


class PreviousFilesLocator(object):
    def __init__(self):
        self.wd_base_prev = []
        self.clusters_filename_pairs = None
        self.wd_trees = None
        self.home_for_results = None
        self.species_tree_rooted_ids_fn = None
        self.species_tree_unrooted_ids_fn = None
        self.species_ids_lines = None
        self.q_assign = False       # an --assign analysis (Log.txt lists WorkingDirectory_Core)
                
    def get_home_for_results(self):
        return self.home_for_results

    def get_start_from_blast(self):
        return self.wd_base_prev

    def get_start_from_ogs(self):
        return self.wd_base_prev, self.clusters_filename_pairs

    def get_start_from_trees(self):
        return self.wd_base_prev, self.clusters_filename_pairs, self.wd_trees, self.species_tree_unrooted_ids_fn, self.species_tree_rooted_ids_fn


class PreviousFilesLocatorNew(PreviousFilesLocator):
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        if not continuation_dir.endswith("/"): continuation_dir += "/"
        # self.home_for_results = continuationDir + "../"
        self.home_for_results = os.path.abspath(continuation_dir + "..") + os.sep
        if (options.q_start_from_fasta and not options.q_start_from_blast):
            # there are no files to find
            return
        if not self._is_new_dir_structure(continuation_dir): raise Unprocessable("Input directory structure is not processable as new structure")
        self._process_log(os.path.join(continuation_dir, "Log.txt"))
    
    def _is_new_dir_structure(self, input_dir):
        return os.path.exists(os.path.join(input_dir, "Log.txt"))
    
    def _process_log(self, log_fn):
        """
        Get all relevant data from log file. 
        Checks the paths saved do exist still
        Should work with relevant paths to allow directory to move
        Other methods can then check that the data required for a particular run is available
        """
        wd_base_anchor = None
        earlier_wds = []     # WorkingDirectory_Core/_Previous lines after the last WorkingDirectory_Base
        with open(log_fn, 'r') as infile:
            for line in infile:
                for key in EARLIER_WD_KEYS:
                    if line.startswith(key + ": "):
                        earlier_wds.append(line.rstrip()[len(key) + 2:])
                        self.q_assign = key == "WorkingDirectory_Core"

                if line.startswith("Species used:"):
                    self.species_ids_lines = ""
                    line = next(infile)
                    while line.rstrip() != "":
                        self.species_ids_lines += line
                        line = next(infile)
                wd_base_str = "WorkingDirectory_Base: "
                wd_trees_str = "WorkingDirectory_Trees: "
                clusters_str = "FN_Orthogroups: "
                if line.startswith(wd_base_str): 
                    wd_base_anchor = self._in_this_results_dir(line.rstrip()[len(wd_base_str):], log_fn, 1)
                    if not os.path.exists(wd_base_anchor):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(wd_base_anchor[:-1])
                        path, d_res = os.path.split(path)
                        # wd_base_anchor = os.path.split(logFN)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        wd_base_anchor = os.path.abspath(os.path.split(log_fn)[0] + ("/../%s/%s/" % (d_res, d_wd)))
                        if not os.path.exists(wd_base_anchor):
                            print("ERROR: Missing directory: %s" % wd_base_anchor)
                            util.fail()
                    earlier_wds = []
                    self.wd_base_prev = [wd_base_anchor]
                    self.wd_trees = wd_base_anchor
                if line.startswith(clusters_str): 
                    clusters_fn_full_path = line.rstrip()[len(clusters_str):]
                    self.clusters_filename_pairs = self._in_this_results_dir(clusters_fn_full_path, log_fn, 2)
                    if not os.path.exists(self.clusters_filename_pairs):
                        # try to see if it's a relative directory to current one
                        path, clusters_fn = os.path.split(self.clusters_filename_pairs)
                        path, d_wd = os.path.split(path)
                        path, d_res = os.path.split(path)
                        # self.clustersFilename_pairs = os.path.split(logFN)[0] + ("/../%s/%s/%s" %(d_res, d_wd, clusters_fn))
                        self.clusters_filename_pairs = os.path.abspath(os.path.split(log_fn)[0] + ("/../%s/%s/%s" %(d_res, d_wd, clusters_fn)))
                        if not os.path.exists(self.clusters_filename_pairs):
                            print("ERROR: Missing orthogroups file: %s or %s" % (self.clusters_filename_pairs, clusters_fn_full_path))
                            util.fail()
#                    self._get_ogs_file(wd_ogs_path)
                if line.startswith(wd_trees_str): 
                    self.wd_trees = self._in_this_results_dir(line.rstrip()[len(wd_trees_str):], log_fn, 1)
                    if not os.path.exists(self.wd_trees):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(self.wd_trees[:-1])
                        path, d_res = os.path.split(path)
                        self.wd_trees = os.path.split(log_fn)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        self.wd_trees = os.path.abspath(self.wd_trees)
                        if not os.path.exists(self.wd_trees):
                            print("ERROR: Missing directory: %s" % self.wd_trees)
                            util.fail()
                    self.species_tree_rooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_rooted_ids.txt") 
                    self.species_tree_unrooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_unrooted_ids.txt")

        if wd_base_anchor is not None:
            # The earlier working directories the analysis reads from: listed in
            # Log.txt (this version), or followed from previous_wd.txt (older
            # versions, which wrote one in each working directory).
            if earlier_wds:
                for d in earlier_wds:
                    if not os.path.isdir(d):
                        print("ERROR: Missing directory: %s" % d)
                        util.fail()
                self.wd_base_prev = [wd_base_anchor] + earlier_wds
            else:
                self.wd_base_prev = self.get_wd_base_chain(wd_base_anchor)

    @staticmethod
    def _in_this_results_dir(recorded, log_fn, n_parts):
        """
        A path recorded in Log.txt, in the results directory being restarted
        if it has it there: the last n_parts of the recorded path (e.g.
        "WorkingDirectory/") under the directory of logFN. Paths are recorded
        absolute, so for a results directory that was copied or moved they
        name the original location, which must not be read or written then
        (whatever the copy is called).
        A working directory is used here only if it is a complete one (it has
        the run's SpeciesIDs.txt): older versions restarted into a new results
        directory whose WorkingDirectory held only the new files, and whose
        Log.txt names the results directory it continued from.
        Otherwise the recorded path (checked by the caller).
        """
        trailing = recorded.endswith(os.sep)
        parts = recorded.rstrip(os.sep).split(os.sep)[-n_parts:]
        local = os.path.join(os.path.dirname(os.path.abspath(log_fn)), *parts)
        if n_parts == 1 and not os.path.exists(os.path.join(local, "SpeciesIDs.txt")):
            return recorded
        if os.path.exists(local):
            return local + (os.sep if trailing else "")
        return recorded

    @staticmethod           
    def get_wd_base_chain(wd_base_anchor):
        chain = [wd_base_anchor]
        while os.path.exists(chain[-1] + "previous_wd.txt"):
            with open(chain[-1] + "previous_wd.txt", 'r') as infile:
                wd = infile.readline().rstrip()
                if not os.path.exists(wd):
                    # try to see if it's a relative directory to current one
                    path, d_wd = os.path.split(wd[:-1])
                    path, d_res = os.path.split(path)
                    wd = wd_base_anchor + ("/../../%s/%s/" % (d_res, d_wd))
                    wd = os.path.abspath(wd)
                    if wd[-1] != os.sep:
                        wd += os.sep
                chain.append(wd)
        return chain


class PreviousFilesLocatorState(PreviousFilesLocator):
    """
    The files of a results directory from its run_state.json (see run_state):
    the working directories it reads from (base), its orthogroups file, its
    trees directory. Raises Unprocessable if the directory has no state
    (e.g. made by an older version): Log.txt is then read instead.
    """
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        results_dir, state = run_state.find(continuation_dir)
        if not state.get("dirs", {}).get("base"):
            raise Unprocessable("No run_state.json in %s" % continuation_dir)
        self.results_dir = results_dir
        self.home_for_results = os.path.abspath(os.path.join(results_dir, "..")) + os.sep
        if options.q_start_from_fasta and not options.q_start_from_blast:
            return
        dirs = run_state.dirs(results_dir, state)
        for d in dirs["base"]:
            if not os.path.isdir(d):
                print("ERROR: A working directory of this analysis is missing: %s" % d)
                print("(recorded in %s)" % run_state.state_path(results_dir))
                util.fail()
        self.wd_base_prev = dirs["base"]
        self.wd_trees = dirs.get("trees") or self.wd_base_prev[0]
        self.clusters_filename_pairs = run_state.file(results_dir, state, "orthogroups")
        self.species_ids_lines = state.get("species", {}).get("ids_lines")
        self.species_tree_rooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_rooted_ids.txt")
        self.species_tree_unrooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_unrooted_ids.txt")


class PreviousFilesLocatorOld(PreviousFilesLocator):
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        if not continuation_dir.endswith("/"): continuation_dir += "/"
        self.home_for_results = continuation_dir + "OrthoFinder/"
        if options.q_start_from_groups or options.q_start_from_trees:
            # User can specify it using clusters_id_pairs file, process this first to get the workingDirectory
            ogs_dir = os.path.abspath(continuation_dir + "../" if options.q_start_from_trees else continuation_dir)

            self.wd_base_prev, self.orthofinder_results_dir, self.clusters_filename_pairs = self._get_ogs_file(ogs_dir)
            if options.q_start_from_trees:
                self._find_from_trees(continuation_dir, options.species_tree_fn, unrooted=True)

            elif options.q_start_from_species_trees:
                self._find_from_trees(continuation_dir, options.species_tree_fn)
        elif options.q_start_from_blast:
            if self._is_working_directory(continuation_dir): 
                self.wd_base_prev = continuation_dir
            elif self._is_working_directory(os.path.join(continuation_dir, "WorkingDirectory") + os.sep):
                self.wd_base_prev = os.path.join(continuation_dir, "WorkingDirectory") + os.sep
            else:
                self.wd_base_prev = continuation_dir   # nothing much to do, set this as the one to try and fail later
        if self.wd_base_prev:
            self.wd_base_prev = [self.wd_base_prev]
        else:
            self.wd_base_prev = []
                
    def _get_ogs_file(self, user_arg):
        """returns the WorkingDirectory, ResultsDirectory and clusters_id_pairs filename"""
        q_specified_results_file = False
        if user_arg == None:
            print("ERROR: orthofinder_results_directory has not been specified")
            util.fail()
        if os.path.isfile(user_arg):
            fn = os.path.split(user_arg)[1]
            if ("clusters_OrthoFinder_" not in fn) or ("txt_id_pairs.txt" not in fn):
                print("ERROR:\n    %s\nis neither a directory or a clusters_OrthoFinder_*.txt_id_pairs.txt file." % user_arg)
                util.fail()
            q_specified_results_file = True
            # user has specified specific results file
        elif user_arg[-1] != os.path.sep: 
            user_arg += os.path.sep
        
        # find required files
        if q_specified_results_file:
            orthofinder_working_dir = os.path.split(user_arg)[0] + os.sep
            if not self._is_working_directory(orthofinder_working_dir):
                print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s" % orthofinder_working_dir)
                util.fail()
        else:
            orthofinder_working_dir = os.path.split(user_arg)[0] if q_specified_results_file else user_arg
            if not self._is_working_directory(orthofinder_working_dir):
                orthofinder_working_dir = user_arg + "WorkingDirectory" + os.sep   
                if not self._is_working_directory(orthofinder_working_dir):
                    print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s\nor\n   %s\n" % (user_arg, orthofinder_working_dir))
                    util.fail()
                
        if q_specified_results_file:
            print("\nUsing orthogroups in file:\n    %s" % user_arg)
            return orthofinder_working_dir, orthofinder_working_dir, user_arg
        else:     
            # identify orthogroups file
            clusters_files = glob.glob(orthofinder_working_dir + "clusters_OrthoFinder_*.txt_id_pairs.txt")
            orthogroup_files = glob.glob(orthofinder_working_dir + "OrthologousGroups*.txt") + glob.glob(orthofinder_working_dir + "Orthogroups*.txt")
            if orthofinder_working_dir != user_arg:
                orthogroup_files += glob.glob(user_arg + "OrthologousGroups*.txt")
                orthogroup_files += glob.glob(user_arg + "Orthogroups*.txt")
            # User may have specified a WorkingDirectory and results could be in directory above
            if len(orthogroup_files) < len(clusters_files):
                # orthogroupFiles += glob.glob(userArg + ".." + os.sep + "OrthologousGroups*.txt")
                # orthogroupFiles += glob.glob(userArg + ".." + os.sep + "Orthogroups*.txt")
                og_path = os.path.abspath(user_arg + "..") 
                orthogroup_files += glob.glob(os.path.join(og_path, "OrthologousGroups*.txt"))
                orthogroup_files += glob.glob(os.path.join(og_path,"Orthogroups*.txt"))
            clusters_files = sorted(clusters_files)
            orthogroup_files = sorted(orthogroup_files)
            if len(clusters_files) > 1 or len(orthogroup_files) > 1:
                print("ERROR: Results from multiple OrthoFinder runs found\n")
                print("Tab-delimiter Orthogroups*.txt/OrthologousGroups*.txt files:")
                for fn in orthogroup_files:
                    print("    " + fn)
                print("With corresponding cluster files:")
                for fn in clusters_files:
                    print("    " + fn)
                print("\nPlease run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line")
                util.fail()        
                
            if len(clusters_files) != 1 or len(orthogroup_files) != 1:
                print("ERROR: Results not found in <orthofinder_results_directory> or <orthofinder_results_directory>/WorkingDirectory")
                print("\nCould not find:\n    Orthogroups*.txt/OrthologousGroups*.txt\nor\n    clusters_OrthoFinder_*.txt_id_pairs.txt")
                util.fail()
                
            print("\nUsing orthogroups in file:\n    %s" % orthogroup_files[0])
            print("and corresponding clusters file:\n    %s" % clusters_files[0])
            return orthofinder_working_dir, user_arg, clusters_files[0]
            
    def _is_working_directory(self, orthofinder_working_dir):
        ok = True
        ok = ok and len(glob.glob(orthofinder_working_dir + "clusters_OrthoFinder_*.txt_id_pairs.txt")) > 0
        ok = ok and len(glob.glob(orthofinder_working_dir + "Species*.fa")) > 0
        return ok
        
    def _find_from_trees(self, orthologs_dir, user_species_tree, unrooted=False):
        """
        if userSpeciesTree == None: Use existing tree
        """
        print("\nFind from trees:")
        print((orthologs_dir, user_species_tree))
        self.wd_trees = os.path.join(orthologs_dir, "WorkingDirectory") + os.sep
        # Find species tree
        if user_species_tree == None:
            possibilities = ["SpeciesTree_ids_0_rooted.txt", "SpeciesTree_ids_1_rooted.txt", "SpeciesTree_user_ids.txt", "SpeciesTree_unrooted_0_rooted.txt", "STAG_SpeciesTree_ids_0_rooted.txt"] # etc (only need to determine if unique)
            n_trees = 0
            for p in possibilities:
                for d in [self.wd_trees, self.wd_trees + "Trees_ids/"]:
                    fn = d + p
                    if os.path.exists(fn): 
                        n_trees += 1
                        species_tree_fn = fn
            if n_trees == 0:
                print("\nERROR: There is a problem with the specified directory. The rooted species tree %s or %s is not present." % (possibilities[0], possibilities[2]))
                print("Please rectify the problem or alternatively use the -s option to specify the species tree to use.\n")
                util.fail()
            if n_trees > 1:
                print("\nERROR: There is more than one rooted species tree in the specified directory structure. Please use the -s option to specify which species tree should be used\n")
                util.fail()
            if unrooted and "unrooted" in fn:
                self.species_tree_unrooted_ids_fn = species_tree_fn
            else:
                self.species_tree_rooted_ids_fn = species_tree_fn
        else:
            if not os.path.exists(user_species_tree):
                print("\nERROR: %s does not exist\n" % user_species_tree)
                util.fail()
            if unrooted:
                self.species_tree_unrooted_ids_fn = species_tree_fn
            else:
                self.species_tree_rooted_ids_fn = species_tree_fn


def initialise_file_handler(
        options, 
        fasta_dir=None, 
        continuation_dir=None, 
        results_dir_non_default=None, 
        pickle_dir_non_default=None,
        working_dir="",
    ):
    """
    Creates a file handler object which will determine the location of all the files:
    Results will be under the user specified directory of the default results location. Defaults:
        - New, from start: 
                FastaDir/OrthoFinder/Results_Date
            or
                resultsDir_nonDefault/Results_Date
        - New, continuation: Existing_OrthoFinder_Dir/Results_Date
        - Old, continuation: 
                ContinuationDir/OrthoFinder/Results_Date
            or 
                resultsDir_nonDefault/Results_Date
        
    
    Implementation
    1. Working out if an old directory structure is being used
    2. Construct and appropriate PreviousFilesLocator if necessary - this locates all required files
    3. Pass this to FileHandler - this creates the directory structure required for this run
    4. if error: print and exit
    5. Return FileHandler
    
    Tasks:
    - Switch this round, I can tell if it's and old or new directory right from the start - read log and check info present,
    perhaps just psss it to the new file handler and let it decide if everything is there
    """
    # 1 & 2
    # If starting from scratch, no need for a PreviousFileLocator
    if fasta_dir is not None and fasta_dir[-1] != os.sep:
        fasta_dir += os.sep
    if options.q_start_from_fasta and not options.q_start_from_blast:
        pfl = None
        base_dir = results_dir_non_default if results_dir_non_default != None else fasta_dir + "OrthoFinder/"
    else:
        try:
            try:
                # The locations recorded in the results directory (run_state.json)
                pfl = PreviousFilesLocatorState(options, continuation_dir)
            except Unprocessable:
                # Otherwise from its Log.txt (results of older versions)
                pfl = PreviousFilesLocatorNew(options, continuation_dir)
            # don't create any new directory, it already exists
            base_dir = pfl.get_home_for_results()

        except Unprocessable:
            pfl = PreviousFilesLocatorOld(options, continuation_dir)
            if not pfl.wd_base_prev:
                print("ERROR: Couldn't process as an OrthoFinder results directory: %s" % continuation_dir)
                util.fail()
            base_dir = results_dir_non_default if results_dir_non_default is not None else pfl.get_home_for_results()
    if not os.path.exists(base_dir): os.mkdir(base_dir)
    if base_dir[-1] != os.sep:
        base_dir += os.sep
    # 3 
    # RefactorDS - this might be suitable as a constructor now
    # base_dir - should now exist
    
    """The previous file locator should decide where the output directory should be rooted? 
    Rules:
    - If starting from Fasta then Fasta/OrthoFinder/Results_Date
          - Or, SpecifiedDirectory/Results_Date if user specified
    - If starting from a previous new-structure directory (TopLevel/Results_X) then TopLevel/Results_Date
    - If starting from a previous old-structure directory then, as high up as we can go and still be in the directory structure:
        - Fasta/Results_OldDate/OrthoFinder/Results_Date
    """
    if not working_dir:
        handler.create_output_directories(options, pfl, base_dir, fasta_dir)
    else:
        if working_dir[-1] != os.sep:
            working_dir += os.sep
        handler.wd_base = [working_dir]
        handler.wd_current = working_dir
        handler.rd1 = os.path.dirname(working_dir)
        handler.wd_trees = working_dir
