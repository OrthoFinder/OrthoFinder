import os
from ..utils import util, files, program_caller
from ..utils.util import printer
from ..comparative_genomics import orthologues
from .. import __version__, g_mclInflation, nThreadsDefault, __location__
from . import helpinfo, species_info
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

configfile_location = os.path.join(__location__, "run")

# Default DIAMOND custom scoring matricies and their corresponding gapopen and gapextend values
diamond_sm_options = {
    "BLOSUM45": [{2: 14}, {3: (10, 13), 2: (12, 16), 1: (16, 19)}],
    "BLOSUM50": [{2: 13}, {3: (9, 13), 2: (12, 16), 1: (15, 19)}],
    "BLOSUM62": [{1: 11}, {2: (6, 11), 1: (9, 13)}],
    "BLOSUM80": [{1: 10}, {2: [6, 7, 8, 9, 13, 25], 1: (9, 11)}],
    "BLOSUM90": [{1: 10}, {2: (6, 11), 1: (9, 11)}],
    "PAM250": [{2: 14}, {3: (11, 15), 2: (13, 17), 1: (17, 21)}],
    "PAM70": [{1: 10}, {2: (6, 8), 1: (9, 11)}],
    "PAM30": [{1: 9}, {2: (5, 7), 1: (8, 10)}],
}

diamond_sm_options_table = """
The following matrices are supported by DIAMOND, with the default being BLOSUM62.

+----------+---------------------------------+-------------------------+
|  Matrix  |         Supported values        |  Default gap penalties  |
+----------+---------------------------------+-------------------------+
| BLOSUM45 | (10-13)/3; (12-16)/2; (16-19)/1 |           14/2          |
| BLOSUM50 | (9-13)/3; (12-16)/2; (15-19)/1  |           13/2          |
| BLOSUM62 | (6-11)/2; (9-13)/1              |           11/1          |
| BLOSUM80 | (6-9)/2; 13/2; 25/2; (9-11)/1   |           10/1          |
| BLOSUM90 | (6-9)/2; (9-11)/1               |           10/1          |
| PAM250   | (11-15)/3; (13-17)/2; (17-21)/1 |           14/2          |
| PAM70    | (6-8)/2; (9-11)/1               |           10/1          |
| PAM30    | (5-7)/2; (8-10)/1               |           9/1           |
+----------+---------------------------------+-------------------------+
NOTE: int (gap open) / int (gap extend)

EXPLANATION:
The default gap open and gap extend penalties for BLOSUM62 are 11 and 1, respectively. Apart from the default gap penalties, BLOSUM62 also support gap extend to be either 2 or 1, with different gap open values allowed. For example, when gap extend is 2, gap open can be chosen between 6 and 11.

USAGE 1:
Available scoring matrices supported by DIAMOND:

1. orthofinder -f ExampleData
When no scoring matrix is provided, OrthoFinder will run the default BLOSUM62 with gap penalties 11/1.

2. orthofinder -f ExampleData --matrix BLOSUM45
When no gap penalties are specified, OrthoFinder will use the default gap penalties.

3. orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2 --gapopen 10
When specifying the gap penalties, gap extend penalty must be defined before the gap open penalty

4. orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2
If the gap open penalty is unavailable, OrthoFinder will use the largest gap open in allowed range defined by the provided gap extend penalty. 

USAGE 2:
Custom scoring matrices:

1. orthofinder -f ExampleData -S diamond_custom --custom-matrix scoring-matrix-file.txt --gapopen 10 --gapextend 2
When use a custom scoring matrix, the search program needs be changed to the custom version.
In the meantime, the gap penalties must be provided. There are no default gap penalties available in OrthoFinder for a custom scoring matrix.
The order of the gap extend and gap open penalties are unimportant in this case. All entries in a custom scoring matrix must be integers.

ATTENTION:
By default, the output directory for each run is saved inside the data folder under the OrthoFinder directory with the name convetion in the following format. Naming format: Results + abbreviated month name + date + (number of runs, the first one is unshown). Example: Results_Jan08_2

To distinguish the results obtained via different scoring matrices, attaching `-efn` to the command can enable OrthoFinder to save the output in a folder with the corresponding scoring matrix name, gap penalties as well as the used search program, MSA method and tree method appending to the convetional name.
Command: orthofinder -f ExampleData --matrix BLOSUM80 -efn
Example: Results_Feb28_5_BLOSUM80-10-1-diamond-mafft-fasttree
(By default, orthofinder uses DIAMOND as the search program, MAFFT as the MSA method, FASTTREE as the tree method)
"""

# Control
class Options(object):  #
    def __init__(self):
        self.nBlast = nThreadsDefault
        self.nProcessAlg = None
        self.qFastAdd = False  # Add species in near-linear time
        self.qStartFromBlast = False  # remove, just store BLAST to do
        self.qStartFromFasta = False  # local to argument checking
        self.qStartFromGroups = False
        self.qStartFromTrees = False
        self.qStopAfterPrepare = False
        self.qStopAfterGroups = False
        self.qStopAfterSeqs = False
        self.qStopAfterAlignments = False
        self.qStopAfterTrees = False
        self.qMSATrees = True  # Updated default
        self.qAddSpeciesToIDs = True
        self.qTrim = True
        self.gathering_version = (1, 0)  # < 3 is the original method
        self.search_program = "diamond"
        self.msa_program = "famsa"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.name = None  # name to identify this set of results
        self.qDoubleBlast = True
        self.qSplitParaClades = False
        self.qPhyldog = False
        self.speciesXMLInfoFN = None
        self.speciesTreeFN = None
        self.mclInflation = g_mclInflation
        self.dna = False
        self.fewer_open_files = (
            True  # By default only open O(n) orthologs files at a time
        )
        self.save_space = False  # On complete, have only one orthologs file per species
        self.v2_scores = False
        self.root_from_previous = False
        self.score_matrix = None
        self.gapopen = None
        self.gapextend = None  #
        self.extended_filename = False
        self.method_threads = "1"
        self.rm_gene_trees = True
        self.rm_resolved_gene_trees = True
        self.cmd_order = "descending"  # "ascending"
        self.threshold = None
        self.method_threads_large = None
        self.method_threads_small = None
        self.old_version = False
        self.fix_files = True
        self.config = None
        self.min_seq = 4

    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)


def period_of_day(hour):
    if hour < 12:
        return "in the morning"
    elif hour < 18:
        return "in the afternoon"
    else:
        return "in the evening"

def ordinal(n):
    if 11 <= n % 100 <= 13:
        suffix = "th"
    else:
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return str(n) + suffix

def GetDirectoryArgument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        util.Fail()
    directory = os.path.abspath(args.pop(0))
    if not os.path.isfile(directory) and directory[-1] != os.sep:
        directory += os.sep
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
    return directory

def GetFileArgument(arg: str) -> str:
    file_path = os.path.abspath(arg)
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        print("Directory points to the file doesn't exist: %s" % directory)
        util.Fail()
    if not os.path.isfile(file_path):
        print("Specified file doesn't exist: %s" % file_path)
        util.Fail()
    return file_path


def GetScoreMatrix(matrixid: str):
    if matrixid.upper() in diamond_sm_options:
        return matrixid.upper()
    else:
        if os.path.isfile(matrixid):
            return matrixid
        else:
            raise Exception("The custom scoring matrix file doesn't exist!")


def GetGapExtend(
        matrixid: str, 
        gapextend: Optional[str] = None,
        scoring_matrix_info = "--matrix"
    ):

    if gapextend is None and scoring_matrix_info == "--matrix":
        if matrixid.upper() in diamond_sm_options:
            return str([*diamond_sm_options[matrixid][0].keys()][0])
        else:
            print(f"{matrixid} is not allowed by DIAMOND")
            print(f"Available scoring matrices allowed by DIAMOND: {[*diamond_sm_options.keys()]}")
            util.Fail()
        
    if gapextend is None and scoring_matrix_info == "--custom-matrix":
        print("--gapextend must be provided when using '--custom-matrix'")
        util.Fail()
    
    if len(gapextend) != 0:
        try:
            gapextend_penalty = abs(int(gapextend))
        except ValueError as e:
            print(f"Unrecognisable --gapextend value {gapextend}!")
            print(f"The gapextend penalty needs to be positive integers!")
            util.Fail()

        if matrixid.upper() in diamond_sm_options:
            allowed_gapextend = [*diamond_sm_options[matrixid][1].keys()]
            if gapextend_penalty in allowed_gapextend:
                return str(gapextend_penalty)
            else:
                raise Exception(
                    f"User defined --gapextend penalty is not allowed by DIAMOND. Acceptable gapextend is {allowed_gapextend}"
                )
        else:
            return str(gapextend_penalty)
    else:
        print(f"Unrecognisable --gapextend value {gapextend}!")
        print(f"The gapextend penalty needs to be positive integers!")
        util.Fail()


def GetGapOpen(
        matrixid: str, 
        gapopen: Optional[str] = None, 
        gapextend: Optional[str] = None,
        scoring_matrix_info = "--matrix",
    ):

    if gapopen is None and gapextend is None and scoring_matrix_info == "--matrix":
        if matrixid.upper() in diamond_sm_options:
            return str([*diamond_sm_options[matrixid][0].values()][0])
        else:
            print(f"{matrixid} is not allowed by DIAMOND")
            print(f"Available scoring matrices allowed by DIAMOND: {[*diamond_sm_options.keys()]}")
            util.Fail()

    if gapopen is None and gapextend is not None and scoring_matrix_info == "--matrix":
        if matrixid.upper() in diamond_sm_options:
            gapextend = int(GetGapExtend(matrixid, gapextend))
            gapopen_range = diamond_sm_options[matrixid][1][gapextend]
            if isinstance(gapopen_range, tuple):
                return str(gapopen_range[1])
            elif isinstance(gapopen_range, list):
                return str(max(gapopen_range))
            
    if gapopen is None and gapextend is not None and scoring_matrix_info == "--custom-matrix":
        print("--gapopen must be provided when using '--custom-matrix'")
        util.Fail()

    if len(gapopen) != 0:
        try:
            gapopen_penalty = abs(int(gapopen))
        except ValueError as e:
            print(f"Unrecognisable --gapopen value {gapopen}!")
            print(f"The gapopen penalty needs to be positive integers!")
            util.Fail()

        if matrixid.upper() in diamond_sm_options:
            allowed_gapextend = [*diamond_sm_options[matrixid][1].keys()]
            gapextend = int(GetGapExtend(matrixid, gapextend))
            if gapextend in allowed_gapextend:
                gapopen_range = diamond_sm_options[matrixid][1][gapextend]
                if isinstance(gapopen_range, tuple):
                    if (
                        gapopen_penalty >= gapopen_range[0]
                        and gapopen_penalty <= gapopen_range[1]
                    ):
                        return str(gapopen_penalty)
                    else:
                        raise Exception(
                            f"User defined --gapopen penalty is not allowed by DIAMOND. Acceptable gapopen when gapextend={gapextend} is between {gapopen_range}"
                        )

                elif isinstance(gapopen_range, list):
                    if gapopen_penalty in gapopen_range:
                        return str(gapopen_penalty)
                    else:
                        raise Exception(
                            f"User defined --gapopen penalty is not allowed by DIAMOND. When gapextend={gapextend} the acceptable gapopen penalties are {gapopen_range}"
                        )

        else:
            return str(gapopen_penalty)
        
    else:
        print(f"Unrecognisable --gapopen value {gapextend}!")
        print(f"The gapeopen penalty needs to be positive integers!")
        util.Fail()

def GetProgramCaller():
    config_file = os.path.join(configfile_location, "config.json")
    pc = program_caller.ProgramCaller(
        config_file if os.path.exists(config_file) else None
    )
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        pc.Add(pc_user)
    return pc


def ProcessArgs(args):
    """
    Workflow
    | 1. Fasta Files | 2.  Prepare files    | 3.   Blast    | 4. Orthogroups    | 5.   Gene Trees     | 6.   Reconciliations/Orthologues   |

    Options
    Start from:
    -f: 1,2,..,6    (start from fasta files, --fasta)
    -b: 4,5,6       (start from blast results, --blast)
    -fg: 5,6         (start from orthogroups/do orthologue workflow, --from-groups)
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
        prog_caller.Add(pc_user)

    if "--config" in args:
        config_index = args.index("--config")
        config_file = args[config_index + 1]
        config_path = GetFileArgument(config_file)
        if os.path.exists(config_path):
            pc_user = program_caller.ProgramCaller(config_path)
            prog_caller.Add(pc_user)
        args.remove("--config")
        args.remove(config_file)
        
    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.PrintHelp(prog_caller)
        util.Success()

    if args[0] == "-v" or args[0] == "--version":
        printer.print(f"[dark_goldenrod]OrthoFinder[/dark_goldenrod]:v[deep_sky_blue2]{__version__}[/deep_sky_blue2]")
        util.Success()

    if args[0] == "-sm" or args[0] == "--scoring-matrix":
        printer.print(diamond_sm_options_table)
        util.Success()
            

    options = Options()
    fastaDir = None
    continuationDir = None
    resultsDir_nonDefault = None
    pickleDir_nonDefault = None
    q_selected_msa_options = False
    q_selected_tree_options = False
    q_selected_search_option = False
    user_specified_M = False
    scoring_matrix_info = "--matrix"

    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """
    

    while len(args) > 0:
        arg = args.pop(0)

        if arg == "-f" or arg == "--fasta":
            if options.qStartFromFasta:
                print("Repeated argument: -f/--fasta\n")
                util.Fail()
            options.qStartFromFasta = True
            fastaDir = GetDirectoryArgument(arg, args)

        elif arg == "-b" or arg == "--blast":
            if options.qStartFromBlast:
                print("Repeated argument: -b/--blast\n")
                util.Fail()
            options.qStartFromBlast = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "--assign":
            options.qFastAdd = True
            fastaDir = GetDirectoryArgument(arg, args)

        elif arg == "--core":
            options.qFastAdd = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-fg" or arg == "--from-groups":
            if options.qStartFromGroups:
                print("Repeated argument: -fg/--from-groups\n")
                util.Fail()
            options.qStartFromGroups = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-ft" or arg == "--from-trees":
            if options.qStartFromTrees:
                print("Repeated argument: -ft/--from-trees\n")
                util.Fail()
            options.qStartFromTrees = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nBlast = int(arg)
                if options.nBlast > nThreadsDefault:
                    print(
                        "WARNING: The specified number of BLAST threads exceed the number of cores: %s"
                        % str(nThreadsDefault)
                    )
                    print(
                        "The number of BLAST threads is reset to %s"
                        % str(nThreadsDefault)
                    )
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()

        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
                if options.nProcessAlg > nThreadsDefault:
                    print(
                        "WARNING: The specified number of OrthoFinder threads exceed the number of cores: %s"
                        % str(nThreadsDefault)
                    )
                    print(
                        "The number of OrthoFinder threads is reset to %s"
                        % str(nThreadsDefault)
                    )
                    options.nProcessAlg = nThreadsDefault
            except:
                print(
                    "Incorrect argument for number of OrthoFinder threads: %s\n" % arg
                )
                util.Fail()

        elif arg == "-mt" or arg == "--method-threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.method_threads = str(arg)
                if int(options.method_threads) > nThreadsDefault:
                    print(
                        "WARNING: The specified number of threads for external commands exceed the number of cores: %d"
                        % nThreadsDefault
                    )
                    print(
                        "The number of threads for external commands is reset to %d"
                        % nThreadsDefault
                    )
                    options.method_threads = str(nThreadsDefault)
            except:
                print("Incorrect argument for number of method threads: %s\n" % arg)
                util.Fail()

        elif arg == "-mtl" or arg == "--method-threads-large":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.method_threads_large = str(arg)
                if int(options.method_threads_large) > nThreadsDefault:
                    print(
                        "WARNING: The specified number of threads for external commands to process large file exceed the number of cores: %d"
                        % nThreadsDefault
                    )
                    print(
                        "The number of threads for external commands to process large file is reset to %d"
                        % nThreadsDefault
                    )
                    options.method_threads_large = str(nThreadsDefault)

            except:
                print(
                    "Incorrect argument for number of method threads for the large files: %s\n"
                    % arg
                )
                util.Fail()

        elif arg == "-mts" or arg == "--method-threads-small":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.method_threads_small = str(arg)
                if int(options.method_threads_small) > nThreadsDefault:
                    print(
                        "WARNING: The specified number of threads for external commands to process small file exceed the number of cores: %"
                        % nThreadsDefault
                    )
                    print(
                        "The number of threads for external commands to process small file is reset to %d"
                        % nThreadsDefault
                    )
                    options.method_threads_small = str(nThreadsDefault)
            except:
                print(
                    "Incorrect argument for number of method threads for the small files: %s\n"
                    % arg
                )
                util.Fail()

        elif arg == "--order":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.cmd_order = str(arg)
            except:
                print(
                    "Incorrect argument for command order: %s. It must be either descending or ascending.\n"
                    % arg
                )
                util.Fail()

        elif arg == "--threshold":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.threshold = int(arg)
            except:
                print(
                    "Incorrect argument for threshold: %s. Values must be between 0 and 100 inclusive.\n"
                    % arg
                )
                util.Fail()

        elif arg == "-ms" or arg == "--min-seq":
            try:
                arg = int(args.pop(0))
                options.min_seq = arg if arg >= 4 else 4 
            except:
                print(f"Incorrect argument {arg} for the minmum number of sequence. Values must be an integer equal to or greater than 4.")
                util.Fail()
        
        elif arg == "--friendly":
            
            console = Console()
            now = datetime.now()
            time_str = now.strftime("%H:%M")
            period = period_of_day(now.hour)
            formatted_date = f"{now.strftime('%A')} the {ordinal(now.day)} {now.strftime('%B, %Y')}"
            console.print(f"Hey Yo! How's going Bro? It's {time_str} {period} {formatted_date}", style="orange3")
            console.print("You know what you are using right now??", style="orange3")
            console.print("The best software on this planet!!! :muscle::sunglasses:", style="orange3")
           
        elif arg == "-rmgt" or arg == "--rm-gene-trees":
            options.rm_gene_trees = False

        elif arg == "-rmrgt" or arg == "--rm-resolved-gene-trees":
            options.rm_resolved_gene_trees = False

        elif arg == "--old-version":
            options.old_version = True

        elif arg == "-1":
            options.qDoubleBlast = False

        elif arg == "-d" or arg == "--dna":
            options.dna = True

        elif arg == "-X":
            options.qAddSpeciesToIDs = False

        elif arg == "-y":
            options.qSplitParaClades = True

        elif arg == "-z":
            options.qTrim = False

        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.mclInflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s\n" % arg)
                util.Fail()

        elif arg == "-c1":
            print("\nThe option 'c1' has been renamed '--c-homologs'")
            util.Fail()

        elif arg == "--c-homologs":
            options.gathering_version = (3, 2)

        elif arg == "--save-space":
            options.save_space = True

        elif arg == "--no-fix-files":
            options.fix_files = False

        elif arg == "-x" or arg == "--orthoxml":
            if options.speciesXMLInfoFN:
                print("Repeated argument: -x/--orthoxml")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesXMLInfoFN = args.pop(0)

        elif arg == "-n" or arg == "--name":
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.name = args.pop(0)
            while options.name.endswith("/"):
                options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]):
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()

        elif arg == "-o" or arg == "--output":
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)

            while resultsDir_nonDefault.endswith("/"):
                resultsDir_nonDefault = resultsDir_nonDefault[:-1]

            resultsDir_nonDefault += "/"
            if os.path.exists(resultsDir_nonDefault):
                print(
                    "ERROR: non-default output directory already exists: %s\n"
                    % resultsDir_nonDefault
                )
                util.Fail()

            if " " in resultsDir_nonDefault:
                print(
                    "ERROR: non-default output directory cannot include spaces: %s\n"
                    % resultsDir_nonDefault
                )
                util.Fail()
            checkDirName = resultsDir_nonDefault

            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]

            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print(
                    "ERROR: location '%s' for results directory '%s' does not exist.\n"
                    % (path, newDir)
                )
                util.Fail()

        elif arg == "-s" or arg == "--speciestree":
            if options.speciesTreeFN:
                print("Repeated argument: -s/--speciestree")
                util.Fail()

            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.speciesTreeFN = args.pop(0)

        elif arg == "--scores-v2":
            options.v2_scores = True

        elif arg == "-S" or arg == "--search":
            choices = ["blast"] + prog_caller.ListSearchMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0)
            if arg in choices:
                options.search_program = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-M" or arg == "--method":
            arg_M_or_msa = arg
            user_specified_M = True
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0)
            if arg == "msa":
                options.qMSATrees = True
            elif arg == "phyldog":
                options.qPhyldog = True
                options.recon_method = "phyldog"
                options.qMSATrees = False
            elif arg == "dendroblast":
                options.qMSATrees = False
            else:
                print("Invalid argument for option %s: %s" % (arg_M_or_msa, arg))
                print("Valid options are 'dendroblast' and 'msa'\n")
                util.Fail()

        elif arg == "-A" or arg == "--msa_program":
            choices = ["mafft"] + prog_caller.ListMSAMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.msa_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-T" or arg == "--tree_program":
            choices = ["fasttree"] + prog_caller.ListTreeMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.tree_program = arg
                q_selected_tree_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-R" or arg == "--recon_method":
            choices = ["of_recon", "dlcpar", "dlcpar_convergedsearch", "only_overlap"]
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.recon_method = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-p":
            pickleDir_nonDefault = GetDirectoryArgument(arg, args)

        elif arg == "-op" or arg == "--only-prepare":
            options.qStopAfterPrepare = True

        elif arg == "-og" or arg == "--only-groups":
            options.qStopAfterGroups = False

        elif arg == "-os" or arg == "--only-seqs":
            options.qStopAfterSeqs = False

        elif arg == "-oa" or arg == "--only-alignments":
            options.qStopAfterAlignments = False

        elif arg == "-ot" or arg == "--only-trees":
            options.qStopAfterTrees = False

        elif arg == "-h" or arg == "--help":
            helpinfo.PrintHelp(prog_caller)
            util.Success()

        elif arg == "--matrix" or arg == "--custom-matrix":
            options.score_matrix = GetScoreMatrix(args.pop(0))
            if arg == "--custom-matrix":
                scoring_matrix_info = "--custom-matrix"

        elif arg == "-ge" or arg == "--gapextend":
            options.gapextend = args.pop(0)

        elif arg == "-go" or arg == "--gapopen":
            options.gapopen =  args.pop(0)
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

        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()

    if "diamond" in options.search_program and not options.score_matrix:
        options.score_matrix = "BLOSUM62"

    if options.score_matrix:
        if not options.gapextend and not options.gapopen:
            options.gapextend = GetGapExtend(options.score_matrix)
            options.gapopen = GetGapOpen(options.score_matrix)

        elif not options.gapopen and options.gapextend:
            options.gapopen = GetGapOpen(
                options.score_matrix, gapextend=options.gapextend
            )

    # set a default for number of algorithm threads
    if options.nProcessAlg is None:
        options.nProcessAlg = min(16, max(1, int(options.nBlast / 8)))

    if options.nBlast < 1:
        print(
            "ERROR: Number of '-t' threads cannot be fewer than 1, got %d"
            % options.nBlast
        )
        util.Fail()

    if options.nProcessAlg < 1:
        print(
            "ERROR: Number of '-a' threads cannot be fewer than 1, got %d"
            % options.nProcessAlg
        )
        util.Fail()

    # check argument combinations
    if not (
        options.qStartFromFasta
        or options.qStartFromBlast
        or options.qStartFromGroups
        or options.qStartFromTrees
        or options.qFastAdd
    ):
        print(
            "ERROR: Please specify the input directory for OrthoFinder using one of the options: '-f', '-b', '-fg' or '-ft', '--assign'."
        )
        util.Fail()

    if options.qFastAdd:
        if (
            options.qStartFromFasta
            or options.qStartFromBlast
            or options.qStartFromGroups
            or options.qStartFromTrees
        ):
            print(
                "ERROR: Incompatible options used with --assign, cannot accept: '-f', '-b', '-fg' or '-ft'"
            )
            util.Fail()
        if fastaDir is None:
            print(
                "ERROR: '--assign' option also requires '--core' directory to be specified"
            )
            util.Fail()
        if continuationDir is None:
            print(
                "ERROR: '--core' option also requires '--assign' directory to be specified"
            )
            util.Fail()
        if not options.qMSATrees:
            print(
                "ERROR: --assign requires MSA trees, option '-M dendroblast' is invalid"
            )
            util.Fail()

    if options.qStartFromFasta and (
        options.qStartFromTrees or options.qStartFromGroups
    ):
        print(
            "ERROR: Incompatible arguments, -f (start from fasta files) and"
            + (
                " -fg (start from orthogroups)"
                if options.qStartFromGroups
                else " -ft (start from trees)"
            )
        )
        util.Fail()

    if options.qStartFromBlast and (
        options.qStartFromTrees or options.qStartFromGroups
    ):
        print(
            "ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and"
            + (
                " -fg (start from orthogroups)"
                if options.qStartFromGroups
                else " -ft (start from trees)"
            )
        )
        util.Fail()

    if options.qStartFromTrees and options.qStartFromGroups:
        print(
            "ERROR: Incompatible arguments, -fg (start from orthogroups) and -ft (start from trees)"
        )
        util.Fail()

    if options.qStopAfterSeqs and (not options.qMSATrees):
        print(
            "ERROR: Argument '-os' (stop after sequences) also requires option '-M msa'"
        )
        util.Fail()

    if options.qStopAfterAlignments and (not options.qMSATrees):
        print(
            "ERROR: Argument '-oa' (stop after alignments) also requires option '-M msa'"
        )
        util.Fail()

    if (q_selected_msa_options or q_selected_tree_options) and (
        not options.qMSATrees and not options.qPhyldog
    ):
        print(
            "ERROR: Argument '-A' or '-T' (multiple sequence alignment/tree inference program) also requires option '-M msa'"
        )
        util.Fail()

    if options.qPhyldog and (not options.speciesTreeFN):
        print("ERROR: Phyldog currently needs a species tree to be provided")
        util.Fail()

    if resultsDir_nonDefault != None and (
        (not options.qStartFromFasta) or options.qStartFromBlast
    ):
        print(
            "ERROR: Incompatible arguments, -o (non-default output directory) can only be used with a new OrthoFinder run using option '-f'"
        )
        util.Fail()

    if options.search_program not in (prog_caller.ListSearchMethods() + ["blast"]):
        print(
            "ERROR: Search program (%s) not configured in config.json file"
            % options.search_program
        )
        util.Fail()

    print()
    util.PrintTime(f"Starting [dark_goldenrod]OrthoFinder[/dark_goldenrod] v[deep_sky_blue2]{__version__}[/deep_sky_blue2]")
    print(
        "%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.nBlast
    )
    print("%d thread(s) for [dark_goldenrod]OrthoFinder[/dark_goldenrod] algorithm\n" % options.nProcessAlg)

    if options.qFastAdd and not q_selected_msa_options:
        print("INFO: For --assign defaulting to 'famsa' to reduce RAM usage\n")
        # options.msa_program = "mafft_memsave"
        options.msa_program = "famsa"

    if options.qFastAdd and not q_selected_tree_options:
        print( "INFO: For --assign defaulting to 'FastTree -fastest' to reduce RAM usage\n")
        options.tree_program = "fasttree_fastest"
        # options.tree_program = "veryfasttree"
    
    if options.search_program == "diamond":
        # check gapextend
        options.gapextend = GetGapExtend(
            options.score_matrix, 
            options.gapextend, 
            scoring_matrix_info
        )

        if options.score_matrix in diamond_sm_options and options.gapextend is None:
            print("The gapopen penalty cannot be define before gapextend")
            util.Fail()

        # check gapopen
        options.gapopen = GetGapOpen(
            options.score_matrix, 
            options.gapopen, 
            options.gapextend,
            scoring_matrix_info
        )
   
    if resultsDir_nonDefault is not None:
        resultsDir_nonDefault = os.path.abspath(resultsDir_nonDefault) + os.path.sep

    return (
        prog_caller, 
        options,
        fastaDir,
        continuationDir,
        resultsDir_nonDefault,
        pickleDir_nonDefault,
        user_specified_M,
    )


def DeleteDirectoryTree(d):
    if os.path.exists(d):
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)


def CheckOptions(options, speciesToUse):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = list(
            species_info.SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN()).values()
        )
        orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any(
        (
            options.qStopAfterPrepare,
            options.qStopAfterGroups,
            options.qStopAfterSeqs,
            options.qStopAfterAlignments,
            options.qStopAfterTrees,
        )
    )
    if q_do_orthologs:
        n_sp = len(speciesToUse)
        wd = files.FileHandler.GetWorkingDirectory_Write()
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
            DeleteDirectoryTree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                util.number_open_files_exception_advice(len(speciesToUse), False)
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                util.Fail()
            else:
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                print(
                    "ERROR: Attempted to open required files for OrthoFinder run but an unexpected error occurred. \n\nStacktrace:"
                )
                raise
    return options
