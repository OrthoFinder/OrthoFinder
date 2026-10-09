import os
import re
from ..utils import util, file_handler, program_caller
from ..utils.printing import printer
from ..orthology import workflow
from .. import __version__, g_mcl_inflation, n_threads_default, __location__
from . import helpinfo, previous_run
import shutil
import traceback
from typing import Optional
import time
from rich.console import Console
from datetime import datetime

try:
    from rich import print
except ImportError:
    ...
from ..utils import printing

configfile_location = os.path.join(__location__, "run")

# # Default DIAMOND custom scoring matricies and their corresponding gapopen and gapextend values
# diamond_sm_options = {
#     "BLOSUM45": [{2: 14}, {3: (10, 13), 2: (12, 16), 1: (16, 19)}],
#     "BLOSUM50": [{2: 13}, {3: (9, 13), 2: (12, 16), 1: (15, 19)}],
#     "BLOSUM62": [{1: 11}, {2: (6, 11), 1: (9, 13)}],
#     "BLOSUM80": [{1: 10}, {2: [6, 7, 8, 9, 13, 25], 1: (9, 11)}],
#     "BLOSUM90": [{1: 10}, {2: (6, 11), 1: (9, 11)}],
#     "PAM250": [{2: 14}, {3: (11, 15), 2: (13, 17), 1: (17, 21)}],
#     "PAM70": [{1: 10}, {2: (6, 8), 1: (9, 11)}],
#     "PAM30": [{1: 9}, {2: (5, 7), 1: (8, 10)}],
# }

# diamond_sm_options_table = """
# The following matrices are supported by DIAMOND, with the default being BLOSUM62.

# +----------+---------------------------------+-------------------------+
# |  Matrix  |         Supported values        |  Default gap penalties  |
# +----------+---------------------------------+-------------------------+
# | BLOSUM45 | (10-13)/3; (12-16)/2; (16-19)/1 |           14/2          |
# | BLOSUM50 | (9-13)/3; (12-16)/2; (15-19)/1  |           13/2          |
# | BLOSUM62 | (6-11)/2; (9-13)/1              |           11/1          |
# | BLOSUM80 | (6-9)/2; 13/2; 25/2; (9-11)/1   |           10/1          |
# | BLOSUM90 | (6-9)/2; (9-11)/1               |           10/1          |
# | PAM250   | (11-15)/3; (13-17)/2; (17-21)/1 |           14/2          |
# | PAM70    | (6-8)/2; (9-11)/1               |           10/1          |
# | PAM30    | (5-7)/2; (8-10)/1               |           9/1           |
# +----------+---------------------------------+-------------------------+
# NOTE: int (gap open) / int (gap extend)

# EXPLANATION:
# The default gap open and gap extend penalties for BLOSUM62 are 11 and 1, respectively. Apart from the default gap penalties, BLOSUM62 also support gap extend to be either 2 or 1, with different gap open values allowed. For example, when gap extend is 2, gap open can be chosen between 6 and 11.

# USAGE 1:
# Available scoring matrices supported by DIAMOND:

# 1. orthofinder -f ExampleData
# When no scoring matrix is provided, OrthoFinder will run the default BLOSUM62 with gap penalties 11/1.

# 2. orthofinder -f ExampleData --matrix BLOSUM45
# When no gap penalties are specified, OrthoFinder will use the default gap penalties.

# 3. orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2 --gapopen 10
# When specifying the gap penalties, gap extend penalty must be defined before the gap open penalty

# 4. orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2
# If the gap open penalty is unavailable, OrthoFinder will use the largest gap open in allowed range defined by the provided gap extend penalty.

# USAGE 2:
# Custom scoring matrices:

# 1. orthofinder -f ExampleData -S diamond_custom --custom-matrix scoring-matrix-file.txt --gapopen 10 --gapextend 2
# When use a custom scoring matrix, the search program needs be changed to the custom version.
# In the meantime, the gap penalties must be provided. There are no default gap penalties available in OrthoFinder for a custom scoring matrix.
# The order of the gap extend and gap open penalties are unimportant in this case. All entries in a custom scoring matrix must be integers.

# ATTENTION:
# By default, the output directory for each run is saved inside the data folder under the OrthoFinder directory with the name convetion in the following format. Naming format: Results + abbreviated month name + date + (number of runs, the first one is unshown). Example: Results_Jan08_2

# To distinguish the results obtained via different scoring matrices, attaching `-efn` to the command can enable OrthoFinder to save the output in a folder with the corresponding scoring matrix name, gap penalties as well as the used search program, MSA method and tree method appending to the convetional name.
# Command: orthofinder -f ExampleData --matrix BLOSUM80 -efn
# Example: Results_Feb28_5_BLOSUM80-10-1-diamond-mafft-fasttree
# (By default, orthofinder uses DIAMOND as the search program, MAFFT as the MSA method, FASTTREE as the tree method)
# """


# Control
class Options(object):  #
    def __init__(self):
        self.n_blast = n_threads_default
        self.n_process_alg = None
        self.q_fast_add = False  # Add species in near-linear time
        self.q_start_from_blast = False  # remove, just store BLAST to do
        self.q_start_from_fasta = False  # local to argument checking
        self.q_start_from_groups = False
        self.q_start_from_trees = False
        self.q_start_from_species_trees = False
        self.q_stop_after_prepare = False
        self.q_stop_after_groups = False
        self.q_stop_after_mcl_groups = False
        self.q_stop_after_seqs = False
        self.q_stop_after_alignments = False
        self.q_stop_after_trees = False
        self.q_stop_after_species_trees = False
        self.q_msa_trees = True  # Updated default
        self.q_add_species_to_ids = True
        self.q_trim = True
        self.gathering_version = (1, 0)  # < 3 is the original method
        self.search_program = "diamond"
        self.search_program_given = False   # -S on the command line
        self.search_program_from_run = False   # restored from the run being continued
        self.msa_program = "famsa"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.name = None  # name to identify this set of results
        self.q_double_blast = True
        self.q_split_para_clades = False
        self.q_phyldog = False
        self.orthoxml = False  # -x: also write the orthogroups as OrthoXML and as a web page
        self.tree_view = False  # --tree-view (not in the help): pictures of the trees and a page for them
        self.species_tree_fn = None
        self.mcl_inflation = g_mcl_inflation
        self.dna = False
        self.fewer_open_files = (
            True  # By default only open O(n) orthologs files at a time
        )
        self.save_space = True  # On complete, have only one orthologs file per species
        self.v2_scores = False
        self.root_from_previous = False
        self.score_matrix = "BLOSUM62"  # None
        self.gapopen = "11"  # None
        self.gapextend = "1"  # None
        self.extended_filename = False
        self.method_threads = "1"
        self.rm_gene_trees = True
        self.rm_resolved_gene_trees = True
        self.cmd_order = "descending"  # "ascending"
        self.threshold = None
        self.method_threads_large = None
        self.method_threads_small = None
        self.fix_files = True
        self.config = None
        self.min_seq = 4
        self.msa_min_seq = 10
        self.astral = False
        self.species_tree_program = "stag"  # astral
        self.dynamic_threads = False
        self.n_skip = 50
        # self.save_blast_commands = False
        self.rm_legacy = True
        self.print_info = True
        self.fd_limit = None
        self.use_base = False

    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)


# def period_of_day(hour):
#     if hour < 12:
#         return "in the morning"
#     elif hour < 18:
#         return "in the afternoon"
#     else:
#         return "in the evening"

# def ordinal(n):
#     if 11 <= n % 100 <= 13:
#         suffix = "th"
#     else:
#         suffix = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
#     return str(n) + suffix


def get_directory_argument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        util.fail()
    directory = os.path.abspath(args.pop(0))
    if not os.path.isfile(directory) and directory[-1] != os.sep:
        directory += os.sep
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.fail()
    return directory


def get_file_argument(arg: str) -> str:
    file_path = os.path.abspath(arg)
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        print("Directory points to the file doesn't exist: %s" % directory)
        util.fail()
    if not os.path.isfile(file_path):
        print("Specified file doesn't exist: %s" % file_path)
        util.fail()
    return file_path


# def GetScoreMatrix(matrixid: str):
#     if matrixid.upper() in diamond_sm_options:
#         return matrixid.upper()
#     else:
#         if os.path.isfile(matrixid):
#             return matrixid
#         else:
#             raise Exception("The custom scoring matrix file doesn't exist!")


# def GetGapExtend(
#         matrixid: str,
#         gapextend: Optional[str] = None,
#         scoring_matrix_info = "--matrix"
#     ):

#     if gapextend is None and scoring_matrix_info == "--matrix":
#         if matrixid.upper() in diamond_sm_options:
#             return str([*diamond_sm_options[matrixid][0].keys()][0])
#         else:
#             print(f"{matrixid} is not allowed by DIAMOND")
#             print(f"Available scoring matrices allowed by DIAMOND: {[*diamond_sm_options.keys()]}")
#             util.Fail()

#     if gapextend is None and scoring_matrix_info == "--custom-matrix":
#         print("--gapextend must be provided when using '--custom-matrix'")
#         util.Fail()

#     if len(gapextend) != 0:
#         try:
#             gapextend_penalty = abs(int(gapextend))
#         except ValueError as e:
#             print(f"Unrecognisable --gapextend value {gapextend}!")
#             print(f"The gapextend penalty needs to be positive integers!")
#             util.Fail()

#         if matrixid.upper() in diamond_sm_options:
#             allowed_gapextend = [*diamond_sm_options[matrixid][1].keys()]
#             if gapextend_penalty in allowed_gapextend:
#                 return str(gapextend_penalty)
#             else:
#                 raise Exception(
#                     f"User defined --gapextend penalty is not allowed by DIAMOND. Acceptable gapextend is {allowed_gapextend}"
#                 )
#         else:
#             return str(gapextend_penalty)
#     else:
#         print(f"Unrecognisable --gapextend value {gapextend}!")
#         print(f"The gapextend penalty needs to be positive integers!")
#         util.Fail()


# def GetGapOpen(
#         matrixid: str,
#         gapopen: Optional[str] = None,
#         gapextend: Optional[str] = None,
#         scoring_matrix_info = "--matrix",
#     ):

#     if gapopen is None and gapextend is None and scoring_matrix_info == "--matrix":
#         if matrixid.upper() in diamond_sm_options:
#             return str([*diamond_sm_options[matrixid][0].values()][0])
#         else:
#             print(f"{matrixid} is not allowed by DIAMOND")
#             print(f"Available scoring matrices allowed by DIAMOND: {[*diamond_sm_options.keys()]}")
#             util.Fail()

#     if gapopen is None and gapextend is not None and scoring_matrix_info == "--matrix":
#         if matrixid.upper() in diamond_sm_options:
#             gapextend = int(GetGapExtend(matrixid, gapextend))
#             gapopen_range = diamond_sm_options[matrixid][1][gapextend]
#             if isinstance(gapopen_range, tuple):
#                 return str(gapopen_range[1])
#             elif isinstance(gapopen_range, list):
#                 return str(max(gapopen_range))

#     if gapopen is None and gapextend is not None and scoring_matrix_info == "--custom-matrix":
#         print("--gapopen must be provided when using '--custom-matrix'")
#         util.Fail()

#     if len(gapopen) != 0:
#         try:
#             gapopen_penalty = abs(int(gapopen))
#         except ValueError as e:
#             print(f"Unrecognisable --gapopen value {gapopen}!")
#             print(f"The gapopen penalty needs to be positive integers!")
#             util.Fail()

#         if matrixid.upper() in diamond_sm_options:
#             allowed_gapextend = [*diamond_sm_options[matrixid][1].keys()]
#             gapextend = int(GetGapExtend(matrixid, gapextend))
#             if gapextend in allowed_gapextend:
#                 gapopen_range = diamond_sm_options[matrixid][1][gapextend]
#                 if isinstance(gapopen_range, tuple):
#                     if (
#                         gapopen_penalty >= gapopen_range[0]
#                         and gapopen_penalty <= gapopen_range[1]
#                     ):
#                         return str(gapopen_penalty)
#                     else:
#                         raise Exception(
#                             f"User defined --gapopen penalty is not allowed by DIAMOND. Acceptable gapopen when gapextend={gapextend} is between {gapopen_range}"
#                         )

#                 elif isinstance(gapopen_range, list):
#                     if gapopen_penalty in gapopen_range:
#                         return str(gapopen_penalty)
#                     else:
#                         raise Exception(
#                             f"User defined --gapopen penalty is not allowed by DIAMOND. When gapextend={gapextend} the acceptable gapopen penalties are {gapopen_range}"
#                         )

#         else:
#             return str(gapopen_penalty)

#     else:
#         print(f"Unrecognisable --gapopen value {gapextend}!")
#         print(f"The gapeopen penalty needs to be positive integers!")
#         util.Fail()


def process_args(args):
    """
    Workflow
    | 1. Fasta Files | 2.  Prepare files    | 3.   Blast    | 4. Orthogroups    | 5.   Gene Trees     | 6.   Reconciliations/Orthologs   |

    Options
    Start from:
    -f: 1,2,..,6    (start from fasta files, --fasta)
    -b: 4,5,6       (start from blast results, --blast)
    -fg: 5,6         (start from orthogroups/do ortholog workflow, --from-groups)
    -ft: 6           (start from gene tree/do reconciliation, --from-trees)
    Stop at:
    -op: 2           (only prepare, --only-prepare)
    -og: 4           (orthogroups, --only-groups)
    """

    config_file = os.path.join(configfile_location, "config.json")
    prog_caller = program_caller.ProgramCaller(
        config_file if os.path.exists(config_file) else None
    )
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        prog_caller.add(pc_user)

    if "--config" in args:
        config_index = args.index("--config")
        config_file = args[config_index + 1]
        config_path = get_file_argument(config_file)
        if os.path.exists(config_path):
            pc_user = program_caller.ProgramCaller(config_path)
            prog_caller.add(pc_user)
        args.remove("--config")
        args.remove(config_file)

    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.print_help(prog_caller)
        util.success()

    if args[0] == "-v" or args[0] == "--version":
        printer.print(
            f"[dark_goldenrod]OrthoFinder[/dark_goldenrod]:v[deep_sky_blue2]{__version__}[/deep_sky_blue2]"
        )
        util.success()

    # if args[0] == "-sm" or args[0] == "--scoring-matrix":
    #     printer.print(diamond_sm_options_table)
    #     util.Success()

    options = Options()
    fasta_dir = None
    continuation_dir = None
    results_dir_non_default = None
    pickle_dir_non_default = None
    q_selected_msa_options = False
    q_selected_tree_options = False
    q_selected_search_option = False
    user_specified_m = False
    scoring_matrix_info = "--matrix"

    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologsDir 
    -ft: store orthologsDir 
    -x: orthoxml
    """

    while len(args) > 0:
        arg = args.pop(0)

        if arg == "-f" or arg == "--fasta":
            if options.q_start_from_fasta:
                print("Repeated argument: -f/--fasta\n")
                util.fail()
            options.q_start_from_fasta = True
            fasta_dir = get_directory_argument(arg, args)

        elif arg == "-b" or arg == "--blast":
            if options.q_start_from_blast:
                print("Repeated argument: -b/--blast\n")
                util.fail()
            options.q_start_from_blast = True
            options.use_base = True
            continuation_dir = get_directory_argument(arg, args)

        elif arg == "--assign":
            options.q_fast_add = True
            fasta_dir = get_directory_argument(arg, args)

        elif arg == "--core":
            options.q_fast_add = True
            continuation_dir = get_directory_argument(arg, args)

        elif arg == "-fg" or arg == "--from-groups":
            if options.q_start_from_groups:
                print("Repeated argument: -fg/--from-groups\n")
                util.fail()
            options.q_start_from_groups = True
            options.use_base = True
            continuation_dir = get_directory_argument(arg, args)

        elif arg == "-fgt" or arg == "--from-trees":
            if options.q_start_from_trees:
                print("Repeated argument: -ft/--from-trees\n")
                util.fail()
            options.q_start_from_trees = True
            options.use_base = True
            continuation_dir = get_directory_argument(arg, args)

        elif arg == "-fst" or arg == "--from-species-trees":
            if options.q_start_from_species_trees:
                print("Repeated argument: -fst/--from-species-trees\n")
                util.fail()
            options.q_start_from_species_trees = True
            options.use_base = True
            continuation_dir = get_directory_argument(arg, args)

        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            try:
                options.n_blast = int(arg)
                if options.n_blast > n_threads_default:
                    print(
                        "WARNING: The specified number of BLAST threads exceed the number of cores: %s"
                        % str(n_threads_default)
                    )
                    print(
                        "The number of BLAST threads is reset to %s"
                        % str(n_threads_default)
                    )
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.fail()

        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            try:
                options.n_process_alg = int(arg)
                if options.n_process_alg > n_threads_default:
                    print(
                        "WARNING: The specified number of OrthoFinder threads exceed the number of cores: %s"
                        % str(n_threads_default)
                    )
                    print(
                        "The number of OrthoFinder threads is reset to %s"
                        % str(n_threads_default)
                    )
                    options.n_process_alg = n_threads_default
            except:
                print(
                    "Incorrect argument for number of OrthoFinder threads: %s\n" % arg
                )
                util.fail()

        elif arg == "-mt" or arg == "--method-threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            try:
                options.method_threads = str(arg)
                if int(options.method_threads) > n_threads_default:
                    print(
                        "WARNING: The specified number of threads for external commands exceed the number of cores: %d"
                        % n_threads_default
                    )
                    print(
                        "The number of threads for external commands is reset to %d"
                        % n_threads_default
                    )
                    options.method_threads = str(n_threads_default)
                if int(options.method_threads) < 1:
                    raise ValueError
            except:
                print("Incorrect argument for number of method threads: %s\n" % arg)
                util.fail()
            if int(options.method_threads) > 1:
                # Each alignment/tree command of a program that takes a thread
                # count (FAMSA, IQ-TREE, RAxML-NG) then reserves this many
                # cores: fewer run at once, which is slower for many small ones.
                printer.print(
                    "WARNING: -mt %s: each alignment and tree command that takes a thread count "
                    "(e.g. FAMSA, IQ-TREE, RAxML-NG) will use %s cores, so fewer commands run at "
                    "the same time. Most gene trees are small and run fastest on 1 core each: "
                    "for many species/orthogroups -mt 1 (the default) is usually faster.\n"
                    % (options.method_threads, options.method_threads), style="warning")

        # elif arg == "-mtl" or arg == "--method-threads-large":
        #     if len(args) == 0:
        #         print("Missing option for command line argument %s\n" % arg)
        #         util.Fail()
        #     arg = args.pop(0)
        #     try:
        #         options.method_threads_large = str(arg)
        #         if int(options.method_threads_large) > nThreadsDefault:
        #             print(
        #                 "WARNING: The specified number of threads for external commands to process large file exceed the number of cores: %d"
        #                 % nThreadsDefault
        #             )
        #             print(
        #                 "The number of threads for external commands to process large file is reset to %d"
        #                 % nThreadsDefault
        #             )
        #             options.method_threads_large = str(nThreadsDefault)

        #     except:
        #         print(
        #             "Incorrect argument for number of method threads for the large files: %s\n"
        #             % arg
        #         )
        #         util.Fail()

        # elif arg == "-mts" or arg == "--method-threads-small":
        #     if len(args) == 0:
        #         print("Missing option for command line argument %s\n" % arg)
        #         util.Fail()
        #     arg = args.pop(0)
        #     try:
        #         options.method_threads_small = str(arg)
        #         if int(options.method_threads_small) > nThreadsDefault:
        #             print(
        #                 "WARNING: The specified number of threads for external commands to process small file exceed the number of cores: %"
        #                 % nThreadsDefault
        #             )
        #             print(
        #                 "The number of threads for external commands to process small file is reset to %d"
        #                 % nThreadsDefault
        #             )
        #             options.method_threads_small = str(nThreadsDefault)
        #     except:
        #         print(
        #             "Incorrect argument for number of method threads for the small files: %s\n"
        #             % arg
        #         )
        #         util.Fail()

        # elif arg == "--order":
        #     if len(args) == 0:
        #         print("Missing option for command line argument %s\n" % arg)
        #         util.Fail()
        #     arg = args.pop(0)
        #     try:
        #         options.cmd_order = str(arg)
        #     except:
        #         print(
        #             "Incorrect argument for command order: %s. It must be either descending or ascending.\n"
        #             % arg
        #         )
        #         util.Fail()

        # elif arg == "--threshold":
        #     if len(args) == 0:
        #         print("Missing option for command line argument %s\n" % arg)
        #         util.Fail()
        #     arg = args.pop(0)
        #     try:
        #         options.threshold = int(arg)
        #     except:
        #         print(
        #             "Incorrect argument for threshold: %s. Values must be between 0 and 100 inclusive.\n"
        #             % arg
        #         )
        #         util.Fail()

        elif arg == "-ms" or arg == "--min-seq":
            try:
                arg = int(args.pop(0))
                options.min_seq = arg if arg >= 4 else 4
            except:
                print(
                    f"Incorrect argument {arg} for the minmum number of sequence. Values must be an integer equal to or greater than 4."
                )
                util.fail()
        elif arg == "--msa-min-seq":
            try:
                arg = int(args.pop(0))
                options.msa_min_seq = arg if arg >= 0 else options.msa_min_seq
            except:
                print(
                    f"Incorrect argument {arg} for the minmum number of sequence MSA analysis. Values must be an integer equal to or greater than 4."
                )
                util.fail()

        elif arg == "-fd":
            try:
                arg = args.pop(0)
                soft, hard = re.split(r"[,\s;]+", arg)
                soft = int(soft)
                hard = int(hard)
                if soft > hard:
                    temp = soft
                    soft, hard = hard, temp
                options.fd_limit = (soft, hard) if soft > 0 and hard > 0 else None
            except:
                print(f"Incorrect argument {arg} for the minmum number of sequence.")
                util.fail()

        # elif arg == "--friendly":

        #     console = Console()
        #     now = datetime.now()
        #     time_str = now.strftime("%H:%M")
        #     period = period_of_day(now.hour)
        #     formatted_date = f"{now.strftime('%A')} the {ordinal(now.day)} {now.strftime('%B, %Y')}"
        #     console.print(f"Hey Yo! How's going Bro? It's {time_str} {period} {formatted_date}", style="orange3")
        #     console.print("You know what you are using right now??", style="orange3")
        #     console.print("The best software on this planet!!! :muscle::sunglasses:", style="orange3")

        # elif arg == "-ub" or arg == "--use-base":
        #     options.use_base = True

        elif arg == "-rmgt" or arg == "--rm-gene-trees":
            options.rm_gene_trees = False

        elif arg == "-rmrgt" or arg == "--rm-resolved-gene-trees":
            options.rm_resolved_gene_trees = False

        elif arg == "--no-print-info":
            options.print_info = False

        elif arg == "-ST" or arg == "--species-tree-program":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            if "astral" in arg.lower():
                options.astral = True

        # elif arg == "--save-blast-commands":
        #     options.save_blast_commands = True

        elif arg == "-1":
            options.q_double_blast = False

        elif arg == "-d" or arg == "--dna":
            options.dna = True

        elif arg == "-X":
            options.q_add_species_to_ids = False

        elif arg == "-y":
            options.q_split_para_clades = True

        elif arg == "-z":
            options.q_trim = False

        elif arg == "-nk":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()

            options.n_skip = int(args.pop(0))

        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            try:
                options.mcl_inflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s\n" % arg)
                util.fail()

        # elif arg == "-c1":
        #     print("\nThe option 'c1' has been renamed '--c-homologs'")
        #     util.Fail()

        # elif arg == "--c-homologs":
        #     options.gathering_version = (3, 2)

        elif arg == "-pof" or arg == "--pairwise-ortholog-files":
            options.save_space = False

        elif arg == "--no-fix-files":
            options.fix_files = False

        elif arg == "-rmlg":
            # keep the legacy files (they are removed by default), incl. the
            # HOGs with gene names from before the orthogroups are updated
            options.rm_legacy = False

        elif arg == "-x" or arg == "--orthoxml":
            options.orthoxml = True

        elif arg == "--tree-view":  # not in the help: needs PyQt6 for ete4 to draw the trees
            options.tree_view = True

        elif arg == "-n" or arg == "--name":
            if options.name:
                print("Repeated argument: -n/--name")
                util.fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            options.name = args.pop(0)
            while options.name.endswith("/"):
                options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]):
                print("Invalid symbol for command line argument %s\n" % arg)
                util.fail()

        elif arg == "-o" or arg == "--output":
            if results_dir_non_default != None:
                print("Repeated argument: -o/--output")
                util.fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            results_dir_non_default = args.pop(0)

            while results_dir_non_default.endswith("/"):
                results_dir_non_default = results_dir_non_default[:-1]

            results_dir_non_default += "/"
            # Tests put multiple named Results_* directories under one shared
            # output base, including runs whose FASTA input is ExampleData.
            # Keep the normal CLI safeguard, but allow that shared base in the
            # suite's explicitly isolated test mode.
            if os.path.exists(results_dir_non_default) and os.getenv("ORTHOFINDER_TEST_ISOLATE") != "1":
                print(
                    "ERROR: non-default output directory already exists: %s\n"
                    % results_dir_non_default
                )
                util.fail()

            if " " in results_dir_non_default:
                print(
                    "ERROR: non-default output directory cannot include spaces: %s\n"
                    % results_dir_non_default
                )
                util.fail()
            check_dir_name = results_dir_non_default

            while check_dir_name.endswith("/"):
                check_dir_name = check_dir_name[:-1]

            path, new_dir = os.path.split(check_dir_name)
            if path != "" and not os.path.exists(path):
                print(
                    "ERROR: location '%s' for results directory '%s' does not exist.\n"
                    % (path, new_dir)
                )
                util.fail()

        elif arg == "-s" or arg == "--speciestree":
            if options.species_tree_fn:
                print("Repeated argument: -s/--speciestree")
                util.fail()

            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()

            options.species_tree_fn = args.pop(0)

        elif arg == "--scores-v2":
            options.v2_scores = True

        elif arg == "-S" or arg == "--search":
            choices = prog_caller.list_search_methods()
            if "blast" not in choices:
                choices = ["blast"] + choices   # another name for "blastp"
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()

            arg = args.pop(0)
            # "-S blast" becomes blastn (with -d) or blastp (resolve_search_program)
            if arg in choices:
                options.search_program = arg
                options.search_program_given = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.fail()

        elif arg == "-M" or arg == "--method":
            arg_m_or_msa = arg
            user_specified_m = True
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()

            arg = args.pop(0)
            if arg == "msa":
                options.q_msa_trees = True
            # elif arg == "phyldog":
            #     options.qPhyldog = True
            #     options.recon_method = "phyldog"
            #     options.qMSATrees = False
            elif arg == "dendroblast":
                options.q_msa_trees = False
            else:
                print("Invalid argument for option %s: %s" % (arg_m_or_msa, arg))
                print("Valid options are 'dendroblast' and 'msa'\n")
                util.fail()

        elif arg == "-A" or arg == "--msa-program":
            choices = ["mafft"] + prog_caller.list_msa_methods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            if arg in choices:
                options.msa_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.fail()

        elif arg == "-T" or arg == "--tree-program":
            choices = ["fasttree"] + prog_caller.list_tree_methods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.fail()
            arg = args.pop(0)
            if arg in choices:
                options.tree_program = arg
                q_selected_tree_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.fail()

        # elif arg == "-R" or arg == "--recon_method":
        #     choices = ["of_recon", "dlcpar", "dlcpar_convergedsearch", "only_overlap"]
        #     switch_used = arg
        #     if len(args) == 0:
        #         print("Missing option for command line argument %s\n" % arg)
        #         util.Fail()
        #     arg = args.pop(0)
        #     if arg in choices:
        #         options.recon_method = arg
        #     else:
        #         print("Invalid argument for option %s: %s" % (switch_used, arg))
        #         print("Valid options are: {%s}\n" % (", ".join(choices)))
        #         util.Fail()

        # elif arg == "-p":
        #     pickleDir_nonDefault = get_directory_argument(arg, args)

        elif arg == "-op" or arg == "--only-prepare":
            options.q_stop_after_prepare = True

        elif arg == "-og" or arg == "--only-groups":
            options.q_stop_after_mcl_groups = True


        elif arg == "-ogh" or arg == "--only-hierarchical-groups":
            options.q_stop_after_groups = True

        # elif arg == "-os" or arg == "--only-seqs":
        #     options.qStopAfterSeqs = False

        # elif arg == "-oa" or arg == "--only-alignments":
        #     options.qStopAfterAlignments = False

        elif arg == "-ogt" or arg == "--only-trees":
            options.q_stop_after_trees = True

        elif arg == "-ost" or arg == "--only-species-trees":
            options.q_stop_after_species_trees = True

        elif arg == "-h" or arg == "--help":
            helpinfo.print_help(prog_caller)
            util.success()
        ### ----------------------------------------------------------
        # elif arg == "--matrix" or arg == "--custom-matrix":
        #     options.score_matrix = GetScoreMatrix(args.pop(0))
        #     if arg == "--custom-matrix":
        #         scoring_matrix_info = "--custom-matrix"

        # elif arg == "-ge" or arg == "--gapextend":
        #     options.gapextend = args.pop(0)

        # elif arg == "-go" or arg == "--gapopen":
        #     options.gapopen =  args.pop(0)

        ### ----------------------------------------------------------

        # elif arg == "-ge" or arg == "--gapextend":
        #     options.gapextend = GetGapExtend(options.score_matrix, args.pop(0))

        # elif arg == "-go" or arg == "--gapopen":
        #     if options.score_matrix in diamond_sm_options and options.gapextend is None:
        #         raise Exception("The gapopen penalty cannot be define before gapextend")
        #     options.gapopen = GetGapOpen(
        #         options.score_matrix, args.pop(0), options.gapextend
        #     )

        elif arg == "-efn" or arg == "--extended-filename":
            options.extended_filename = True

        # elif arg == "--dythreads":
        #     options.dynamic_threads = True

        else:
            print("Unrecognised argument: %s\n" % arg)
            util.fail()

    # if "diamond" in options.search_program and not options.score_matrix:
    #     options.score_matrix = "BLOSUM62"

    # if options.score_matrix:
    #     if not options.gapextend and not options.gapopen:
    #         options.gapextend = GetGapExtend(options.score_matrix)
    #         options.gapopen = GetGapOpen(options.score_matrix)

    #     elif not options.gapopen and options.gapextend:
    #         options.gapopen = GetGapOpen(
    #             options.score_matrix, gapextend=options.gapextend
    #         )

    # set a default for number of algorithm threads
    if options.n_process_alg is None:
        options.n_process_alg = min(16, max(1, int(options.n_blast / 8)))

    if options.n_blast < 1:
        print(
            "ERROR: Number of '-t' threads cannot be fewer than 1, got %d"
            % options.n_blast
        )
        util.fail()

    if options.n_process_alg < 1:
        print(
            "ERROR: Number of '-a' threads cannot be fewer than 1, got %d"
            % options.n_process_alg
        )
        util.fail()

    # check argument combinations
    if not (
        options.q_start_from_fasta
        or options.q_start_from_blast
        or options.q_start_from_groups
        or options.q_start_from_trees
        or options.q_fast_add
        or options.q_start_from_species_trees
    ):
        print(
            "ERROR: Please specify the input directory for OrthoFinder using one of the options: '-f', '-b', '-fg' or '-ft', '--assign'."
        )
        util.fail()

    if options.q_fast_add:
        if (
            options.q_start_from_fasta
            or options.q_start_from_blast
            or options.q_start_from_groups
            or options.q_start_from_trees
            or options.q_start_from_species_trees
        ):
            print(
                "ERROR: Incompatible options used with --assign, cannot accept: '-f', '-b', '-fg', '-ft', '-fst'"
            )
            util.fail()
        if fasta_dir is None:
            print(
                "ERROR: '--assign' option also requires '--core' directory to be specified"
            )
            util.fail()
        if continuation_dir is None:
            print(
                "ERROR: '--core' option also requires '--assign' directory to be specified"
            )
            util.fail()
        if not options.q_msa_trees:
            print(
                "ERROR: --assign requires MSA trees, option '-M dendroblast' is invalid"
            )
            util.fail()

    if options.q_start_from_fasta and (
        options.q_start_from_trees
        # or options.qStartFromGroups
        or options.q_start_from_species_trees
    ):
        err_msg = "ERROR: Incompatible arguments, -f (start from fasta files) and"
        if options.q_start_from_groups:
            err_msg += " -fg (start from orthogroups)"
        elif options.q_start_from_trees:
            err_msg += " -ft (start from gene trees)"
        else:
            err_msg += "-fst (start from species trees)"
        print(err_msg)

        util.fail()

    if options.q_start_from_blast and (
        options.q_start_from_trees
        or options.q_start_from_groups
        or options.q_start_from_species_trees
    ):
        err_msg = "ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and"

        if options.q_start_from_groups:
            err_msg += " -fg (start from orthogroups)"
        elif options.q_start_from_groups:
            err_msg += " -ft (start from gene trees)"
        else:
            err_msg += "-fst (start from species trees)"
        print(err_msg)
        util.fail()

    if options.q_start_from_trees and options.q_start_from_groups:
        print(
            "ERROR: Incompatible arguments, -fg (start from orthogroups) and -ft (start from gene trees)"
        )
        util.fail()

    if options.q_start_from_species_trees and options.q_start_from_groups:
        print(
            "ERROR: Incompatible arguments, -fg (start from orthogroups) and -fst (start from species gene trees)"
        )
        util.fail()

    if options.q_start_from_species_trees and options.q_start_from_trees:
        print(
            "ERROR: Incompatible arguments, -ft (start from gene trees) and -fst (start from species gene trees)"
        )
        util.fail()

    if options.q_stop_after_seqs and (not options.q_msa_trees):
        print(
            "ERROR: Argument '-os' (stop after sequences) also requires option '-M msa'"
        )
        util.fail()

    if options.q_stop_after_alignments and (not options.q_msa_trees):
        print(
            "ERROR: Argument '-oa' (stop after alignments) also requires option '-M msa'"
        )
        util.fail()

    if (q_selected_msa_options or q_selected_tree_options) and (
        not options.q_msa_trees and not options.q_phyldog
    ):
        print(
            "ERROR: Argument '-A' or '-T' (multiple sequence alignment/tree inference program) also requires option '-M msa'"
        )
        util.fail()

    # if options.qPhyldog and (not options.speciesTreeFN):
    #     print("ERROR: Phyldog currently needs a species tree to be provided")
    #     util.Fail()

    if results_dir_non_default != None and (
        (not options.q_start_from_fasta) or options.q_start_from_blast
    ):
        print(
            "ERROR: Incompatible arguments, -o (non-default output directory) can only be used with a new OrthoFinder run using option '-f'"
        )
        util.fail()

    if options.search_program not in prog_caller.list_search_methods() + ["blast"]:
        print(
            "ERROR: Search program (%s) not configured in config.json file"
            % options.search_program
        )
        util.fail()

    if options.print_info:
        print()
        util.print_time(
            f"Starting [dark_goldenrod]OrthoFinder[/dark_goldenrod] v[deep_sky_blue2]{__version__}[/deep_sky_blue2]"
        )
        print(
            "%d thread(s) for highly parallel tasks (BLAST searches etc.)"
            % options.n_blast
        )
        print(
            "%d thread(s) for [dark_goldenrod]OrthoFinder[/dark_goldenrod] algorithm\n"
            % options.n_process_alg
        )

    if options.q_fast_add and not q_selected_msa_options:
        if options.print_info:
            print("INFO: For --assign defaulting to 'famsa' to reduce RAM usage\n")
        # options.msa_program = "mafft_memsave"
        options.msa_program = "famsa"

    if options.q_fast_add and not q_selected_tree_options:
        if options.print_info:
            print(
                "INFO: For --assign defaulting to 'FastTree -fastest' to reduce RAM usage\n"
            )
        options.tree_program = "fasttree_fastest"
        # options.tree_program = "veryfasttree"
    ### ---------------------------------------------------------
    # if options.search_program == "diamond":
    #     # check gapextend
    #     options.gapextend = GetGapExtend(
    #         options.score_matrix,
    #         options.gapextend,
    #         scoring_matrix_info
    #     )

    #     if options.score_matrix in diamond_sm_options and options.gapextend is None:
    #         print("The gapopen penalty cannot be define before gapextend")
    #         util.Fail()

    #     # check gapopen
    #     options.gapopen = GetGapOpen(
    #         options.score_matrix,
    #         options.gapopen,
    #         options.gapextend,
    #         scoring_matrix_info
    #     )
    ### ---------------------------------------------------------
    if results_dir_non_default is not None:
        results_dir_non_default = os.path.abspath(results_dir_non_default) + os.path.sep

    return (
        prog_caller,
        options,
        fasta_dir,
        continuation_dir,
        results_dir_non_default,
        pickle_dir_non_default,
        user_specified_m,
    )


def delete_directory_tree(d):
    if os.path.exists(d):
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)


def check_options(options, species_to_use):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.species_tree_fn:
        exp_species = list(
            previous_run.species_name_dict(file_handler.handler.get_species_ids_fn()).values()
        )
        previous_run.check_user_species_tree(options.species_tree_fn, exp_species)

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any(
        (
            options.q_stop_after_prepare,
            options.q_stop_after_groups,
            options.q_stop_after_mcl_groups,
            options.q_stop_after_seqs,
            options.q_stop_after_alignments,
            options.q_stop_after_trees,
            options.q_stop_after_species_trees,
        )
    )
    if q_do_orthologs:
        n_sp = len(species_to_use)
        wd = file_handler.handler.get_working_directory_write()
        wd_files_test = wd + "Files_test/"
        fh = []
        try:
            if not os.path.exists(wd_files_test):
                os.mkdir(wd_files_test)
            for i_sp in range(n_sp):
                di = wd_files_test + "Sp%d/" % i_sp
                if not os.path.exists(di):
                    os.mkdir(di)
                for j_sp in range(
                    1
                ):  # We only create a linear number of ortholog files now
                    fnij = di + "Sp%d.txt" % j_sp
                    fh.append(open(fnij, "w"))
            # create a few extra files to be safe
            for i_extra in range(n_extra):
                fh.append(open(wd_files_test + "Extra%d.txt" % i_extra, "w"))
            # close the files again and delete
            for fhh in fh:
                fhh.close()
            delete_directory_tree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                printing.number_open_files_exception_advice(len(species_to_use), False)
                for fhh in fh:
                    fhh.close()
                delete_directory_tree(wd_files_test)
                util.fail()
            else:
                for fhh in fh:
                    fhh.close()
                delete_directory_tree(wd_files_test)
                print(
                    "ERROR: Attempted to open required files for OrthoFinder run but an unexpected error occurred. \n\nStacktrace:"
                )
                raise
    return options


def resolve_search_program(options, prog_caller):
    """
    The search programs whose commands depend on whether the input is
    nucleotide (-d) or protein:
      - "-S blast" is BLAST+ blastn with -d and blastp otherwise, as
        configured for "blastn"/"blastp" in config.json;
      - MMseqs2 (any configured method that runs it) needs --search-type 3
        for nucleotide sequences, which is then added to its commands with
        a k-mer length (program_caller.add_nucleotide_search_type).
    DIAMOND compares protein sequences only: with -d the default search is
    BLAST+ blastn instead, and DIAMOND given with -S is an error (a restart
    keeps the program of the run it continues, which ran its searches).
    The input is not read here: nucleotide input must be given with -d
    (without it, fasta_processor.processes_new_fasta stops with an error).
    """
    if options.dna and prog_caller.uses_diamond(options.search_program) and not options.search_program_from_run:
        if options.search_program_given:
            print('ERROR: The search program "%s" (DIAMOND) compares protein sequences only, and the '
                  'input is nucleotide (-d). Use "-S blast" (BLAST+ blastn) or "-S mmseqs" (MMseqs2), '
                  'or leave out -S to use blastn.' % options.search_program)
            util.fail()
        print('Search program: nucleotide input (-d) is searched with BLAST+ blastn by default, as '
              'DIAMOND compares protein sequences only ("-S mmseqs" is faster on large inputs)\n')
        if "blastn" in prog_caller.list_search_methods():
            options.search_program = "blastn"
            return options
        options.search_program = "blast"
    q_blast = options.search_program == "blast" and "blast" not in prog_caller.list_search_methods()
    if q_blast:
        options.search_program = "blastn" if options.dna else "blastp"
        print('Search program: "-S blast" uses BLAST+ %s (%s sequences), as configured for "%s" in config.json\n'
              % (options.search_program, "nucleotide" if options.dna else "protein", options.search_program))
    elif options.dna and prog_caller.uses_mmseqs(options.search_program):
        print('Search program: "%s" (MMseqs2) compares the nucleotide sequences as nucleotides: '
              '"--search-type 3 -k %d" is added to its search commands, and no index is created '
              '(see the mmseqs comment in config.json)\n'
              % (options.search_program, program_caller.MMSEQS_NUCLEOTIDE_KMER))
    return options
