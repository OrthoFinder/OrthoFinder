"""
The command-line options of OrthoFinder: process_args reads them into an
Options object, checks they can be used together and returns them with the
directories they name.

The steps of an analysis, and the options that start or stop it at a step:
    | 1. FASTA files | 2. Prepare files | 3. Search | 4. Orthogroups | 5. Gene trees | 6. Orthologs |
    Start from:
    -f:   1-6   (FASTA files, --fasta)
    -b:   4-6   (the search results of an earlier analysis, --blast; with -f: add species to it)
    -fg:  5-6   (its orthogroups, --from-groups)
    -fgt: 6     (its gene trees, --from-trees)
    -fst:       (its species tree, --from-species-trees)
    --assign <dir> --core <dir>: assign the species in <dir> to the orthogroups of a finished analysis
    Stop after:
    -op:  2     (--only-prepare)
    -og:  4     (orthogroups, --only-groups; -ogh: hierarchical orthogroups)
    -ogt: 5     (gene trees, --only-trees; -ost: species tree)
"""
import os
import re
import shutil
import time

from .. import __version__, g_mcl_inflation, n_cores, n_threads_default, __location__
from ..search import diamond_scoring
from ..utils import util, file_handler, program_caller, printing
from ..utils.printing import printer
from . import cli_options, helpinfo, previous_run

try:
    from rich import print
except ImportError:
    ...

configfile_location = os.path.join(__location__, "run")


class Options(object):
    def __init__(self):
        # where the analysis starts and stops
        self.q_start_from_fasta = False
        self.q_start_from_blast = False
        self.q_start_from_groups = False
        self.q_start_from_trees = False
        self.q_start_from_species_trees = False
        self.q_fast_add = False  # --assign: add species to a finished analysis
        self.use_base = False  # continues an earlier analysis (-b, -fg, -fgt, -fst)
        self.q_stop_after_prepare = False
        self.q_stop_after_groups = False
        self.q_stop_after_mcl_groups = False
        self.q_stop_after_seqs = False
        self.q_stop_after_alignments = False
        self.q_stop_after_trees = False
        self.q_stop_after_species_trees = False
        # threads
        self.n_blast = n_threads_default
        self.n_process_alg = None
        self.method_threads = "1"
        self.method_threads_large = None
        self.method_threads_small = None
        self.dynamic_threads = False
        self.cmd_order = "descending"
        self.threshold = None
        # programs and methods
        self.search_program = "diamond"
        self.search_program_given = False   # -S on the command line
        self.search_program_from_run = False   # restored from the run being continued
        self.msa_program = "famsa"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.q_msa_trees = True
        self.q_phyldog = False
        self.species_tree_program = None  # -ST: from the gene trees (None: the default for -M)
        self.n_skip = 50
        self.dna = False
        self.q_double_blast = True
        self.v2_scores = False
        self.score_matrix = "BLOSUM62"  # DIAMOND's --matrix, or a file (--custom-matrix)
        self.gapopen = "11"
        self.gapextend = "1"
        self.score_matrix_given = False  # --matrix, --custom-matrix, --gapopen or --gapextend given
        self.gathering_version = (1, 0)  # < 3 is the original method
        self.mcl_inflation = g_mcl_inflation
        self.min_seq = 4
        self.msa_min_seq = 10
        self.q_trim = True
        self.q_add_species_to_ids = True
        self.q_split_para_clades = False
        self.species_tree_fn = None
        # output
        self.name = None  # name to identify this set of results
        self.extended_filename = False
        self.fewer_open_files = True  # By default only open O(n) orthologs files at a time
        self.save_space = True  # On complete, have only one orthologs file per species
        self.fix_files = True
        self.rm_gene_trees = True
        self.rm_resolved_gene_trees = True
        self.rm_legacy = True
        self.orthoxml = False  # --orthoxml (-x): also write the orthogroups and HOGs as OrthoXML
        self.species_tree_alignment = False  # --concat-alignment: with -ST, also write the concatenated alignment
        self.mcl_orthogroups = False  # -mog: also keep the MCL orthogroups (MCL_Orthogroups.*)
        # pages to look at the results in a web browser (--view: all three)
        self.tree_view = False  # --tree-view: the trees
        self.stats_view = False  # --stats-view: the statistics, with plots
        self.groups_view = False  # --groups-view: the genes of the orthogroups and HOGs
        self.print_info = True
        self.fd_limit = None


class ArgumentError(Exception):
    """A command line that cannot be run: the message says why (and what to do instead)."""


# What each option that starts an analysis does (in messages)
START_DOES = {
    "-f": "starts a new analysis from FASTA files",
    "-b": "continues an earlier analysis from its search results",
    "-fg": "continues an earlier analysis from its orthogroups",
    "-fgt": "continues an earlier analysis from its gene trees",
    "-fst": "continues an earlier analysis from its species tree",
}

# The start options in the order of the steps, as (attribute, short name)
START_ORDER = (
    ("q_start_from_fasta", "-f"),
    ("q_start_from_blast", "-b"),
    ("q_start_from_groups", "-fg"),
    ("q_start_from_trees", "-fgt"),
    ("q_start_from_species_trees", "-fst"),
)
# The pairs of them that can be used together: -f with -b adds species to an
# earlier analysis
START_PAIRS_ALLOWED = {("-f", "-b")}

ASSIGN_USAGE = "orthofinder --assign <new species dir> --core <results dir>"

# Other names for species tree programs (-ST)
SPECIES_TREE_ALIASES = {"astral": "astral-pro"}
# The species tree program --assign uses, unless -ST gives another
ASSIGN_SPECIES_TREE_PROGRAM = "astral-pro"


def assign_species_tree_program(options):
    """The species tree program an --assign analysis uses: the one given with -ST, or ASTRAL-Pro."""
    return options.species_tree_program or ASSIGN_SPECIES_TREE_PROGRAM


def get_directory_argument(arg, args):
    """The directory given after option arg (an absolute path ending with a separator)."""
    directory = os.path.abspath(_value(arg, args))
    if not os.path.isfile(directory) and directory[-1] != os.sep:
        directory += os.sep
    if not os.path.exists(directory):
        raise ArgumentError("Specified directory doesn't exist: %s" % directory)
    return directory


def get_file_argument(arg: str) -> str:
    """The absolute path of an existing file."""
    file_path = os.path.abspath(arg)
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        raise ArgumentError("Directory points to the file doesn't exist: %s" % directory)
    if not os.path.isfile(file_path):
        raise ArgumentError("Specified file doesn't exist: %s" % file_path)
    return file_path


def _value(arg, args):
    """The value given after option arg."""
    if len(args) == 0:
        raise ArgumentError("Missing option for command line argument %s" % arg)
    return args.pop(0)


def _int_value(arg, args, what):
    value = _value(arg, args)
    try:
        return int(value)
    except ValueError:
        raise ArgumentError("Incorrect argument for %s: %s" % (what, value))


def _threads_value(arg, args, what):
    """A number of threads, at most the number of cores."""
    n = _int_value(arg, args, "number of %s" % what)
    if n > n_cores:
        print("WARNING: The specified number of %s exceed the number of cores: %d" % (what, n_cores))
        print("The number of %s is reset to %d" % (what, n_cores))
        n = n_cores
    return n


def _choice(arg, args, choices):
    value = _value(arg, args)
    if value not in choices:
        raise ArgumentError("Invalid argument for option %s: %s\nValid options are: {%s}"
                            % (arg, value, ", ".join(choices)))
    return value


def program_caller_with_configs(args):
    """
    The configured programs: config.json, ~/config_orthofinder_user.json and
    those of --config <file> (whose options are removed from args).
    """
    config_file = os.path.join(configfile_location, "config.json")
    prog_caller = program_caller.ProgramCaller(config_file if os.path.exists(config_file) else None)
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        prog_caller.add(program_caller.ProgramCaller(config_file_user))
    while "--config" in args:
        i = args.index("--config")
        if i + 1 == len(args):
            raise ArgumentError("Missing option for command line argument --config")
        config_path = get_file_argument(args[i + 1])
        prog_caller.add(program_caller.ProgramCaller(config_path))
        del args[i:i + 2]
    return prog_caller


def process_args(args):
    """
    Read the command line (args, a list, which is emptied). Exits for -h and
    -v, and with an error message for a command line that cannot be run.
    Returns (prog_caller, options, fasta_dir, continuation_dir,
    results_dir_non_default, pickle_dir_non_default, user_specified_m).
    """
    try:
        prog_caller = program_caller_with_configs(args)
        if len(args) == 0 or args[0] == "help":
            helpinfo.print_help(prog_caller)
            util.success()
        parsed = _read_options(args, prog_caller)
        _check_options(prog_caller, *parsed)
    except ArgumentError as e:
        print("ERROR: %s\n" % e)
        util.fail()
    options, fasta_dir, continuation_dir, results_dir_non_default, q_selected_msa, q_selected_tree, user_specified_m = parsed

    if options.print_info:
        print()
        util.print_time(f"Starting [dark_goldenrod]OrthoFinder[/dark_goldenrod] v[deep_sky_blue2]{__version__}[/deep_sky_blue2]")
        print("%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.n_blast)
        print("%d thread(s) for [dark_goldenrod]OrthoFinder[/dark_goldenrod] algorithm\n" % options.n_process_alg)

    if options.q_fast_add and not q_selected_msa:
        if options.print_info:
            print("INFO: For --assign defaulting to 'famsa' to reduce RAM usage\n")
        options.msa_program = "famsa"
    if options.q_fast_add and not q_selected_tree:
        if options.print_info:
            print("INFO: For --assign defaulting to 'FastTree -fastest' to reduce RAM usage\n")
        options.tree_program = "fasttree_fastest"

    if results_dir_non_default is not None:
        results_dir_non_default = os.path.abspath(results_dir_non_default) + os.path.sep
    pickle_dir_non_default = None
    return (
        prog_caller,
        options,
        fasta_dir,
        continuation_dir,
        results_dir_non_default,
        pickle_dir_non_default,
        user_specified_m,
    )


class _Reading(object):
    """What has been read of the command line so far."""

    def __init__(self, prog_caller):
        self.prog_caller = prog_caller
        self.options = Options()
        self.fasta_dir = None
        self.continuation_dir = None
        self.results_dir_non_default = None
        self.q_selected_msa = False
        self.q_selected_tree = False
        self.user_specified_m = False
        self.given = set()   # the options that can only be given once, by their first name
        self.scoring = {}    # --matrix, --custom-matrix, --gapopen, --gapextend: their values


def _names(name):
    """All the names of an option, for messages: e.g. -f/--fasta."""
    return "/".join(cli_options.OPTIONS_BY_NAME[name]["names"])


def _once(r, option):
    """An option that can only be given once."""
    if option["names"][0] in r.given:
        raise ArgumentError("Repeated argument: %s" % "/".join(option["names"]))
    r.given.add(option["names"][0])


def _read_options(args, prog_caller):
    """The options of the command line (cli_options), before checking how they are combined."""
    r = _Reading(prog_caller)
    while len(args) > 0:
        arg = args.pop(0)
        option = cli_options.OPTIONS_BY_NAME.get(arg)
        if option is None:
            raise ArgumentError("Unrecognised argument: %s" % arg)
        if "set" in option:
            for attribute, value in option["set"].items():
                setattr(r.options, attribute, value)
        else:
            READERS[option["read"]](r, option, arg, args)

    options = r.options
    if r.scoring:
        try:
            options.score_matrix, options.gapopen, options.gapextend = diamond_scoring.scoring(
                r.scoring.get("--matrix"), r.scoring.get("--custom-matrix"),
                r.scoring.get("--gapopen"), r.scoring.get("--gapextend"))
        except ValueError as e:
            raise ArgumentError(str(e))
        options.score_matrix_given = True

    # the default number of algorithm threads
    if options.n_process_alg is None:
        options.n_process_alg = min(16, max(1, int(options.n_blast / 8)))
    return (options, r.fasta_dir, r.continuation_dir, r.results_dir_non_default, r.q_selected_msa,
            r.q_selected_tree, r.user_specified_m)


# ------------------------------------------- the readers of the options ----
# Each reads the value of an option (cli_options: its "read"): reader(r, option, arg, args)

def _read_start(r, option, arg, args):
    _once(r, option)
    setattr(r.options, option["field"], True)
    if option["field"] == "q_start_from_fasta":
        r.fasta_dir = get_directory_argument(arg, args)
    else:
        r.options.use_base = True
        r.continuation_dir = get_directory_argument(arg, args)


def _read_assign(r, option, arg, args):
    _once(r, option)
    r.options.q_fast_add = True
    r.fasta_dir = get_directory_argument(arg, args)


def _read_core(r, option, arg, args):
    _once(r, option)
    r.options.q_fast_add = True
    r.continuation_dir = get_directory_argument(arg, args)


def _read_threads(r, option, arg, args):
    setattr(r.options, option["field"], _threads_value(arg, args, option["what"]))


def _read_method_threads(r, option, arg, args):
    n = _threads_value(arg, args, "threads for external commands")
    if n < 1:
        raise ArgumentError("Incorrect argument for number of method threads: %d" % n)
    r.options.method_threads = str(n)
    if n > 1:
        # Each alignment/tree command of a program that takes a thread
        # count (FAMSA, IQ-TREE, RAxML-NG) then reserves this many
        # cores: fewer run at once, which is slower for many small ones.
        printer.print(
            "WARNING: -mt %d: each alignment and tree command that takes a thread count "
            "(e.g. FAMSA, IQ-TREE, RAxML-NG) will use %d cores, so fewer commands run at "
            "the same time. Most gene trees are small and run fastest on 1 core each: "
            "for many species/orthogroups -mt 1 (the default) is usually faster.\n"
            % (n, n), style="warning")


def _read_min_seq(r, option, arg, args):
    n = _int_value(arg, args, "the minimum number of sequences (an integer, 4 or more)")
    if n < 4:
        print("WARNING: %s %d: the minimum number of sequences is 4, which is used\n" % (arg, n))
    r.options.min_seq = max(n, 4)


def _read_msa_min_seq(r, option, arg, args):
    what = "the minimum number of sequences for trimming an alignment (an integer, 0 or more)"
    n = _int_value(arg, args, what)
    if n < 0:
        raise ArgumentError("Incorrect argument for %s: %d" % (what, n))
    r.options.msa_min_seq = n


def _read_fd_limit(r, option, arg, args):
    value = _value(arg, args)
    try:
        soft, hard = sorted(int(v) for v in re.split(r"[,\s;]+", value))
    except ValueError:
        raise ArgumentError("Incorrect argument for -fd (the soft and hard limits on the number of open "
                            "files, e.g. 10000,20000): %s" % value)
    r.options.fd_limit = (soft, hard) if soft > 0 and hard > 0 else None


def _read_int(r, option, arg, args):
    setattr(r.options, option["field"], _int_value(arg, args, option["what"]))


def _read_float(r, option, arg, args):
    value = _value(arg, args)
    try:
        setattr(r.options, option["field"], float(value))
    except ValueError:
        raise ArgumentError("Incorrect argument for %s: %s" % (option["what"], value))


def _read_species_tree_program(r, option, arg, args):
    if args and args[0].lower() in SPECIES_TREE_ALIASES:
        args[0] = SPECIES_TREE_ALIASES[args[0].lower()]
    r.options.species_tree_program = _choice(arg, args, r.prog_caller.list_species_tree_methods())


def _read_name(r, option, arg, args):
    _once(r, option)
    r.options.name = _value(arg, args).rstrip("/")
    if " " in r.options.name or "/" in r.options.name:
        raise ArgumentError("Invalid symbol for command line argument %s: %s" % (arg, r.options.name))


def _read_output(r, option, arg, args):
    _once(r, option)
    r.results_dir_non_default = _output_directory(_value(arg, args))


def _read_species_tree(r, option, arg, args):
    _once(r, option)
    r.options.species_tree_fn = _value(arg, args)


def _read_search(r, option, arg, args):
    choices = r.prog_caller.list_search_methods()
    if "blast" not in choices:
        choices = ["blast"] + choices   # another name for "blastp"
    # "-S blast" becomes blastn (with -d) or blastp (resolve_search_program)
    r.options.search_program = _choice(arg, args, choices)
    r.options.search_program_given = True


def _read_method(r, option, arg, args):
    r.user_specified_m = True
    value = _value(arg, args)
    if value == "msa":
        r.options.q_msa_trees = True
    elif value == "dendroblast":
        r.options.q_msa_trees = False
    else:
        raise ArgumentError("Invalid argument for option %s: %s\nValid options are 'dendroblast' and 'msa'"
                            % (arg, value))


def _read_msa(r, option, arg, args):
    r.options.msa_program = _choice(arg, args, ["mafft"] + r.prog_caller.list_msa_methods())
    r.q_selected_msa = True


def _read_tree(r, option, arg, args):
    r.options.tree_program = _choice(arg, args, ["fasttree"] + r.prog_caller.list_tree_methods())
    r.q_selected_tree = True


def _read_scoring(r, option, arg, args):
    # by the long name: --matrix, --custom-matrix, --gapopen, --gapextend
    name = next(n for n in option["names"] if n.startswith("--"))
    if name in r.scoring:
        raise ArgumentError("Repeated argument: %s" % name)
    if name == "--matrix":
        r.scoring[name] = _value(arg, args)
    elif name == "--custom-matrix":
        r.scoring[name] = get_file_argument(_value(arg, args))
    else:
        r.scoring[name] = _int_value(arg, args, "%s (an integer)" % name)


def _read_help(r, option, arg, args):
    helpinfo.print_help(r.prog_caller)
    util.success()


def _read_version(r, option, arg, args):
    printer.print(f"[dark_goldenrod]OrthoFinder[/dark_goldenrod]:v[deep_sky_blue2]{__version__}[/deep_sky_blue2]")
    util.success()


READERS = {
    "start": _read_start,
    "assign": _read_assign,
    "core": _read_core,
    "threads": _read_threads,
    "method_threads": _read_method_threads,
    "min_seq": _read_min_seq,
    "msa_min_seq": _read_msa_min_seq,
    "fd_limit": _read_fd_limit,
    "int": _read_int,
    "float": _read_float,
    "species_tree_program": _read_species_tree_program,
    "name": _read_name,
    "output": _read_output,
    "species_tree": _read_species_tree,
    "search": _read_search,
    "method": _read_method,
    "msa": _read_msa,
    "tree": _read_tree,
    "scoring": _read_scoring,
    "help": _read_help,
    "version": _read_version,
}


def _output_directory(value):
    """The results directory given with -o (it must not exist yet), ending with '/'."""
    results_dir = value.rstrip("/") + "/"
    # Tests put multiple named Results_* directories under one shared
    # output base, including runs whose FASTA input is ExampleData.
    # Keep the normal CLI safeguard, but allow that shared base in the
    # suite's explicitly isolated test mode.
    if os.path.exists(results_dir) and os.getenv("ORTHOFINDER_TEST_ISOLATE") != "1":
        raise ArgumentError("non-default output directory already exists: %s" % results_dir)
    if " " in results_dir:
        raise ArgumentError("non-default output directory cannot include spaces: %s" % results_dir)
    path, new_dir = os.path.split(results_dir.rstrip("/"))
    if path != "" and not os.path.exists(path):
        raise ArgumentError("location '%s' for results directory '%s' does not exist." % (path, new_dir))
    return results_dir


def _check_options(prog_caller, options, fasta_dir, continuation_dir, results_dir_non_default,
                   q_selected_msa, q_selected_tree, user_specified_m):
    """Raise an ArgumentError if the options cannot be used together."""
    if options.n_blast < 1:
        raise ArgumentError("Number of '-t' threads cannot be fewer than 1, got %d" % options.n_blast)
    if options.n_process_alg < 1:
        raise ArgumentError("Number of '-a' threads cannot be fewer than 1, got %d" % options.n_process_alg)

    starts = [name for attribute, name in START_ORDER if getattr(options, attribute)]
    if not starts and not options.q_fast_add:
        raise ArgumentError("Please specify the input directory for OrthoFinder using one of the options: "
                            "'-f', '-b', '-fg', '-fgt', '-fst' or '--assign' (with '--core').")

    if options.q_fast_add:
        if starts:
            message = "--assign cannot be used with %s.\n--assign adds the species in a directory of FASTA files to a finished analysis:\n  %s" % (
                " or ".join(starts), ASSIGN_USAGE)
            if "-f" in starts:
                message += "\nTo analyse all the species together, start a new analysis with -f alone."
            raise ArgumentError(message)
        if continuation_dir is None:
            raise ArgumentError("--assign also needs --core: the results directory of the analysis to add the "
                                "species to\n  %s" % ASSIGN_USAGE)
        if fasta_dir is None:
            raise ArgumentError("--core also needs --assign: the directory of FASTA files of the species to "
                                "add\n  %s" % ASSIGN_USAGE)
        if not options.q_msa_trees:
            raise ArgumentError("--assign requires MSA trees, option '-M dendroblast' is invalid")

    for i, first in enumerate(starts):
        for second in starts[i + 1:]:
            if (first, second) not in START_PAIRS_ALLOWED:
                first_names, first_does = _names(first), START_DOES[first]
                second_names, second_does = _names(second), START_DOES[second]
                raise ArgumentError("Incompatible arguments, %s and %s cannot be used together:\n"
                                    "  %s %s\n  %s %s\nUse one of them."
                                    % (first, second, first_names, first_does, second_names, second_does))

    if (q_selected_msa or q_selected_tree) and not options.q_msa_trees and not options.q_phyldog:
        raise ArgumentError("Argument '-A' or '-T' (multiple sequence alignment/tree inference program) "
                            "also requires option '-M msa'")
    if options.species_tree_program is not None and not options.q_msa_trees:
        raise ArgumentError("Argument '-ST' (species tree program) also requires option '-M msa'")

    if results_dir_non_default is not None and (not options.q_start_from_fasta or options.q_start_from_blast):
        raise ArgumentError("Incompatible arguments, -o (non-default output directory) can only be used with "
                            "a new OrthoFinder run using option '-f'")

    if options.search_program not in prog_caller.list_search_methods() + ["blast"]:
        raise ArgumentError("Search program (%s) not configured in config.json file" % options.search_program)


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


def resolve_scoring(options, prog_caller):
    """
    The scoring matrix with the search program (once a restart has restored
    those of the run it continues): --matrix etc. are for DIAMOND only, and a
    custom matrix is used by the search method of config.json for one
    (diamond -> diamond_custom etc.).
    """
    if options.score_matrix_given and not prog_caller.uses_diamond(options.search_program):
        print('ERROR: --matrix, --custom-matrix, --gapopen and --gapextend are for DIAMOND searches, and the search '
              'program is "%s"%s.' % (options.search_program, " (nucleotide input, -d)" if options.dna else ""))
        util.fail()
    custom = diamond_scoring.is_custom_matrix(options.score_matrix)
    if custom:
        try:
            _, warnings = diamond_scoring.read_custom_matrix(options.score_matrix)
        except ValueError as e:
            print("ERROR: %s\n" % e)
            util.fail()
        for warning in warnings:
            print("WARNING: %s\n" % warning)
    method = options.search_program
    if custom and method in diamond_scoring.CUSTOM_MATRIX_METHODS:
        method = diamond_scoring.CUSTOM_MATRIX_METHODS[method]
    elif not custom and method in diamond_scoring.CUSTOM_MATRIX_METHODS.values():
        method = {v: k for k, v in diamond_scoring.CUSTOM_MATRIX_METHODS.items()}[method]
    if method != options.search_program:
        if method not in prog_caller.list_search_methods():
            print('ERROR: The search method "%s" for %s is not configured in config.json'
                  % (method, "a custom scoring matrix" if custom else "the scoring matrix " + options.score_matrix))
            util.fail()
        print('Search program: "%s", as the scoring matrix is %s\n'
              % (method, "a custom one (--custom-matrix)" if custom else options.score_matrix))
        options.search_program = method


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
    Then the scoring matrix with it (resolve_scoring).
    """
    _resolve_search_method(options, prog_caller)
    resolve_scoring(options, prog_caller)
    return options


def _resolve_search_method(options, prog_caller):
    """The search program for the input (see resolve_search_program)."""
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
