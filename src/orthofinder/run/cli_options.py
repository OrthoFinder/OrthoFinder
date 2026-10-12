"""
The command-line options of OrthoFinder, in groups: process_args reads them,
helpinfo prints them (a section for each group, in this order).

A group: {"title": its section in the help, "options": [option, ...]}. An option:
    "names":   its names, e.g. ("-t", "--threads"); the first is shown in the help
    "help":    its text in the help (rich markup); None: not in the help
    "value":   what it takes, shown in the help ("int", "txt", "dir", "file"); absent: a flag
    "set":     (a flag) {Options attribute: the value it sets}
    "read":    (with a value) the reader of process_args that reads it (READERS)
    "field":   the Options attribute a reader sets, "what": what it is, in messages
    "choices": a function(prog_caller): the values listed under it in the help
--config <file> is read before the others (process_args.program_caller_with_configs).
"""
from .. import g_mcl_inflation, n_threads_default
from ..search import diamond_scoring

OF = "[dark_goldenrod]OrthoFinder[/dark_goldenrod]"


def _dir(text):
    return text.replace("<dir>", "<[bright_magenta]dir[/bright_magenta]>")


INPUT_OPTIONS = {
    "title": "INPUT OPTIONS",
    "options": [
        {"names": ("-f", "--fasta"), "value": "dir", "read": "start", "field": "q_start_from_fasta",
         "help": _dir(f"Run a full {OF} analysis on [red1]FASTA[/red1] format proteomes in <dir>")},
        {"names": ("--assign",), "value": "dir", "read": "assign",
         "help": _dir(f"Assign the species in <dir> to the orthogroups of a finished {OF} analysis, "
                      "requires [dark_cyan]--core[/dark_cyan]")},
        {"names": ("--core",), "value": "dir", "read": "core",
         "help": _dir("The results directory of the finished analysis, for [dark_cyan]--assign[/dark_cyan]")},
        {"names": ("-d", "--dna"), "set": {"dna": True},
         "help": "Input is DNA sequences (searched with BLAST+ blastn by default, or -S mmseqs; DIAMOND cannot be used)."},
        {"names": ("-s", "--speciestree"), "value": "file", "read": "species_tree",
         "help": "User-specified rooted species tree"},
    ],
}

METHOD_OPTIONS = {
    "title": "METHOD OPTIONS",
    "options": [
        {"names": ("-M", "--method"), "value": "txt", "read": "method",
         "help": 'Method for gene tree inference. Options [dark_cyan]"dendroblast"[/dark_cyan] & '
                 '[dark_cyan]"msa"[/dark_cyan] [Default = [dark_cyan]msa[/dark_cyan]]'},
        {"names": ("-S", "--search"), "value": "txt", "read": "search",
         "help": "Sequence search program [Default = [dark_cyan]diamond[/dark_cyan]]",
         "choices": lambda pc: [m for m in pc.list_search_methods()
                                if m not in diamond_scoring.CUSTOM_MATRIX_METHODS.values()]},
        {"names": ("-A", "--msa-program"), "value": "txt", "read": "msa",
         "help": 'MSA program, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]famsa[/dark_cyan]]',
         "choices": lambda pc: [m for m in pc.list_msa_methods() if m != "mafft_memsave"]},
        {"names": ("-T", "--tree-program"), "value": "txt", "read": "tree",
         "help": 'Tree inference method, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]FastTree[/dark_cyan]]',
         "choices": lambda pc: [m for m in pc.list_tree_methods() if m != "raxml-ng"]},
        {"names": ("-ST", "--species-tree-program"), "value": "txt", "read": "species_tree_program",
         "help": 'Species tree program, from the gene trees, requires [dark_cyan]"-M msa"[/dark_cyan] '
                 "[Default = a tree of the concatenated alignment of single-copy orthogroups]",
         "choices": lambda pc: pc.list_species_tree_methods()},
        {"names": ("-I", "--inflation"), "value": "int", "read": "float", "field": "mcl_inflation",
         "what": "MCL inflation parameter",
         "help": f"MCL inflation parameter [Default = [deep_sky_blue2]{g_mcl_inflation:0.1f}[/deep_sky_blue2]]"},
        {"names": ("-ms", "--min-seq"), "value": "int", "read": "min_seq", "help": None},
        {"names": ("--msa-min-seq",), "value": "int", "read": "msa_min_seq", "help": None},
        {"names": ("-nk",), "value": "int", "read": "int", "field": "n_skip", "what": "-nk", "help": None},
        {"names": ("--scores-v2",), "set": {"v2_scores": True}, "help": None},
    ],
}

PARALLEL_OPTIONS = {
    "title": "PARALLEL PROCESSING OPTIONS",
    "options": [
        {"names": ("-t", "--threads"), "value": "int", "read": "threads", "field": "n_blast", "what": "BLAST threads",
         "help": f"Number of parallel sequence search threads [Default = [deep_sky_blue2]{n_threads_default}[/deep_sky_blue2]]"},
        {"names": ("-a", "--algthreads"), "value": "int", "read": "threads", "field": "n_process_alg",
         "what": "OrthoFinder threads", "help": "Number of parallel analysis threads"},
        {"names": ("-mt", "--method-threads"), "value": "int", "read": "method_threads", "help": None},
        {"names": ("-fd",), "value": "txt", "read": "fd_limit", "help": None},
    ],
}

OUTPUT_OPTIONS = {
    "title": "OUTPUT OPTIONS",
    "options": [
        {"names": ("-n", "--name"), "value": "txt", "read": "name", "help": "Name to append to the results directory"},
        {"names": ("-o", "--output"), "value": "txt", "read": "output", "help": "Non-default results directory"},
        {"names": ("-X",), "set": {"q_add_species_to_ids": False}, "help": "Don't add species names to sequence IDs"},
        {"names": ("-pof", "--pairwise-ortholog-files"), "set": {"save_space": False},
         "help": "Save the pairwise ortholog files"},
        {"names": ("-efn", "--extended-filename"), "set": {"extended_filename": True}, "help": None},
        {"names": ("--no-fix-files",), "set": {"fix_files": False}, "help": None},
        {"names": ("-rmgt", "--rm-gene-trees"), "set": {"rm_gene_trees": False}, "help": None},
        {"names": ("-rmrgt", "--rm-resolved-gene-trees"), "set": {"rm_resolved_gene_trees": False}, "help": None},
        # keep the legacy files (they are removed by default), incl. the HOGs
        # with gene names from before the orthogroups are updated
        {"names": ("-rmlg",), "set": {"rm_legacy": False}, "help": None},
        {"names": ("--no-print-info",), "set": {"print_info": False}, "help": None},
    ],
}

ADDITIONAL_OUTPUT_OPTIONS = {
    "title": "ADDITIONAL OUTPUT OPTIONS",
    "options": [
        {"names": ("--orthoxml", "-x"), "set": {"orthoxml": True},
         "help": "The orthogroups and HOGs in OrthoXML format, for other tools: Orthogroups.orthoxml and "
                 "Phylogenetic_Hierarchical_Orthogroups/N*.orthoxml"},
        {"names": ("--concat-alignment",), "set": {"species_tree_alignment": True},
         "help": "With [dark_cyan]-ST[/dark_cyan], also write the concatenated alignment of the orthogroups "
                 "chosen for a species tree (as a run without [dark_cyan]-ST[/dark_cyan] does): "
                 "WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa and "
                 "Species_Tree/Orthogroups_for_concatenated_alignment.txt. For the results of a finished run: "
                 "[dark_cyan]species_tree_alignment[/dark_cyan] <[bright_magenta]dir[/bright_magenta]>"},
        {"names": ("--tree-view",), "set": {"tree_view": True},
         "help": "A page to view the trees in a web browser: each resolved gene tree with its "
                 "alignment (clades can be collapsed, a tree saved as SVG), a species tree explorer "
                 "(the HOGs at each node, an orthogroup's duplications and losses, each species' "
                 "genes) and all genes on the species tree: Resolved_Gene_Trees/Resolved_Gene_Trees.html"},
        {"names": ("--stats-view",), "set": {"stats_view": True},
         "help": "A page of the run's statistics with interactive plots (overall, species, "
                 "orthogroups, duplications, HOGs): Comparative_Genomics_Statistics/Statistics.html"},
        {"names": ("--groups-view",), "set": {"groups_view": True},
         "help": "A page of the genes of each orthogroup and HOG, searchable: Orthogroups/Groups.html"},
        {"names": ("--view",), "set": {"tree_view": True, "stats_view": True, "groups_view": True},
         "help": "All three pages ([dark_cyan]--tree-view[/dark_cyan], [dark_cyan]--stats-view[/dark_cyan], "
                 "[dark_cyan]--groups-view[/dark_cyan]), linked to each other, and OrthoFinder_Results.html to "
                 "open them, a tab for each (also with two of them). With many species they can be large: each "
                 "can also be written on its own. For the results of a finished run: "
                 "[dark_cyan]view_results[/dark_cyan] <[bright_magenta]dir[/bright_magenta]>"},
    ],
}

STOPPING_OPTIONS = {
    "title": "WORKFLOW STOPPING OPTIONS",
    "options": [
        {"names": ("-op", "--only-prepare"), "set": {"q_stop_after_prepare": True},
         "help": "Stop after preparing input files for BLAST"},
        {"names": ("-og", "--only-groups"), "set": {"q_stop_after_mcl_groups": True}, "help": "Stop after MCL Orthogroups"},
        {"names": ("-ogh", "--only-hierarchical-groups"), "set": {"q_stop_after_groups": True}, "help": None},
        {"names": ("-ogt", "--only-trees"), "set": {"q_stop_after_trees": True}, "help": "Stop after gene trees"},
        {"names": ("-ost", "--only-species-trees"), "set": {"q_stop_after_species_trees": True},
         "help": "Stop after species trees"},
    ],
}

RESTART_OPTIONS = {
    "title": "WORKFLOW RESTART COMMANDS",
    "options": [
        {"names": ("-b", "--blast"), "value": "dir", "read": "start", "field": "q_start_from_blast",
         "help": _dir(f"Restart {OF} from pre-computed sequence search results in <dir>")},
        {"names": ("-fg", "--from-groups"), "value": "dir", "read": "start", "field": "q_start_from_groups",
         "help": _dir(f"Restart {OF} from MCL Orthogroups results in <dir>")},
        {"names": ("-fgt", "--from-trees"), "value": "dir", "read": "start", "field": "q_start_from_trees",
         "help": _dir(f"Restart {OF} from gene trees results in <dir>")},
        {"names": ("-fst", "--from-species-trees"), "value": "dir", "read": "start", "field": "q_start_from_species_trees",
         "help": _dir(f"Restart {OF} from species tree results in <dir>")},
    ],
}

OTHER_OPTIONS = {
    "title": "OTHER OPTIONS",
    "options": [
        {"names": ("-y",), "set": {"q_split_para_clades": True},
         "help": "Split paralogous clades below root of a HOG into separate HOGs"},
        {"names": ("-z",), "set": {"q_trim": False},
         "help": "Don't trim MSAs (columns>=[deep_sky_blue2]90[/deep_sky_blue2]% gap, min. alignment length "
                 "[deep_sky_blue2]500[/deep_sky_blue2])"},
        {"names": ("-1",), "set": {"q_double_blast": False}, "help": None},
        {"names": ("--matrix",), "value": "txt", "read": "scoring",
         "help": "Scoring matrix of the DIAMOND searches [Default = [dark_cyan]BLOSUM62[/dark_cyan]]",
         "choices": lambda pc: list(diamond_scoring.DIAMOND_MATRICES)},
        {"names": ("--custom-matrix",), "value": "file", "read": "scoring",
         "help": "A scoring matrix file for the DIAMOND searches, requires [dark_cyan]--gapopen[/dark_cyan] and "
                 "[dark_cyan]--gapextend[/dark_cyan]\n"
                 "In NCBI's format (e.g. the files at https://ftp.ncbi.nlm.nih.gov/blast/matrices/): a line of "
                 "one-letter amino acid codes, then for each code, in that order, a line of the code and its integer "
                 "scores. The 20 standard amino acids are needed; lines starting with # are comments."},
        {"names": ("-go", "--gapopen"), "value": "int", "read": "scoring",
         "help": "Gap open penalty of the DIAMOND searches, [dark_cyan]--gapopen[/dark_cyan] "
                 "[Default = the matrix's, e.g. [deep_sky_blue2]11[/deep_sky_blue2] for BLOSUM62]"},
        {"names": ("-ge", "--gapextend"), "value": "int", "read": "scoring",
         "help": "Gap extend penalty of the DIAMOND searches, [dark_cyan]--gapextend[/dark_cyan] "
                 "[Default = the matrix's, e.g. [deep_sky_blue2]1[/deep_sky_blue2] for BLOSUM62]"},
        {"names": ("-h", "--help"), "read": "help", "help": "Print this help text"},
    ],
}

VERSION_OPTIONS = {
    "title": "VERSION",
    "options": [
        {"names": ("-v", "--version"), "read": "version", "help": "Show the current version number"},
    ],
}

# The groups, in the order of the help
OPTION_GROUPS = (
    INPUT_OPTIONS,
    METHOD_OPTIONS,
    PARALLEL_OPTIONS,
    OUTPUT_OPTIONS,
    STOPPING_OPTIONS,
    RESTART_OPTIONS,
    OTHER_OPTIONS,
    ADDITIONAL_OUTPUT_OPTIONS,
    VERSION_OPTIONS,
)

# Each option by each of its names
OPTIONS_BY_NAME = {name: option for group in OPTION_GROUPS for option in group["options"] for name in option["names"]}


def help_label(option):
    """An option as the help shows it, e.g. "-t <int>"."""
    if "value" not in option:
        return option["names"][0]
    return "%s <[bright_magenta]%s[/bright_magenta]>" % (option["names"][0], option["value"])
