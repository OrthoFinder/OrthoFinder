"""
The command-line options of OrthoFinder, in groups: process_args reads them,
helpinfo prints them (a section for each group, in this order).

A group: {"title": its section in the help, "options": [option, ...]}, with "hidden": True
not in the help (DEVELOPMENT_OPTIONS). An option:
    "names":   its names, e.g. ("-t", "--threads"); the first is shown in the help
    "help":    its text in the help (rich markup); None: not in the help (development-only
               options go in DEVELOPMENT_OPTIONS, hidden, with a help text for developers)
    "value":   what it takes, shown in the help ("int", "num", "txt", "dir", "file"); absent: a flag
    "set":     (a flag) {Options attribute: the value it sets}
    "read":    (with a value) the reader of process_args that reads it (READERS)
    "field":   the Options attribute a reader sets, "what": what it is, in messages
    "choices": a function(prog_caller): the values listed under it in the help
    "preview": True for an option of a feature not yet released (PREVIEW_OPTIONS)
--config <file> is read before the others (process_args.program_caller_with_configs).
"""
from .. import g_mcl_inflation, n_cores, n_threads_default
from ..search import diamond_scoring

OF = "[dark_goldenrod]OrthoFinder[/dark_goldenrod]"

# The options marked "preview" (the pages to view the results, the scoring of the DIAMOND
# searches): True on the development branch (dev_master), where they are in the help and
# can be used; False on main, the released OrthoFinder, where they are neither (the
# script that makes main from dev_master sets it, and removes their documentation).
PREVIEW_OPTIONS = True


def _dir(text):
    return text.replace("<dir>", "<[bright_magenta]dir[/bright_magenta]>")


def _search_choices(prog_caller):
    """The -S choices: the search methods configured (not those for a custom matrix), with
    "blast" (BLAST+ blastp or blastn: process_args.resolve_search_program) after DIAMOND's."""
    methods = [m for m in prog_caller.list_search_methods() if m not in diamond_scoring.CUSTOM_MATRIX_METHODS.values()]
    if "blast" not in methods and ("blastp" in methods or "blastn" in methods):
        i = max((k + 1 for k, m in enumerate(methods) if m.startswith("diamond")), default=0)
        methods.insert(i, "blast")
    return methods


def _path(path):
    """A path in the results."""
    return "[medium_purple]%s[/medium_purple]" % path


def _link(url):
    """A web address, in the colour of the links of the citations, and clickable."""
    return "[dark_cyan][link=%s]%s[/link][/dark_cyan]" % (url, url)


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
         "help": "Input is DNA sequences (searched with [dark_cyan]BLAST+ blastn[/dark_cyan] by default, or [dark_cyan]-S mmseqs[/dark_cyan]; [dark_cyan]DIAMOND[/dark_cyan] cannot be used)."},
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
         "help": "Sequence search program [Default = [dark_cyan]diamond[/dark_cyan]; "
                 "[dark_cyan]blastn[/dark_cyan] with [dark_cyan]-d[/dark_cyan]]. "
                 "[dark_cyan]blast[/dark_cyan]: BLAST+ for either input, [dark_cyan]blastp[/dark_cyan] "
                 "for protein sequences, [dark_cyan]blastn[/dark_cyan] with [dark_cyan]-d[/dark_cyan]",
         "choices": lambda pc: _search_choices(pc)},
        {"names": ("-A", "--msa-program"), "value": "txt", "read": "msa",
         "help": 'MSA program, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]famsa[/dark_cyan]]',
         "choices": lambda pc: [m for m in pc.list_msa_methods() if m != "mafft_memsave"]},
        {"names": ("-T", "--tree-program"), "value": "txt", "read": "tree",
         "help": 'Tree inference method, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]FastTree[/dark_cyan]]',
         "choices": lambda pc: [m for m in pc.list_tree_methods() if m != "raxml-ng"]},
        {"names": ("-ST", "--species-tree-program"), "value": "txt", "read": "species_tree_program",
         "help": 'Species tree program, from the gene trees, requires [dark_cyan]"-M msa"[/dark_cyan] '
                 "[Default = a tree of the concatenated alignment of single-copy orthogroups; "
                 "with --assign, [dark_cyan]astral-pro[/dark_cyan]]",
         "choices": lambda pc: pc.list_species_tree_methods()},
        {"names": ("-I", "--inflation"), "value": "num", "read": "float", "field": "mcl_inflation",
         "what": "MCL inflation parameter",
         "help": f"[dark_cyan]MCL[/dark_cyan] inflation parameter, usually [deep_sky_blue2]1.2[/deep_sky_blue2]-[deep_sky_blue2]5[/deep_sky_blue2]; larger values give smaller orthogroups "
                 f"[Default = [deep_sky_blue2]{g_mcl_inflation:0.1f}[/deep_sky_blue2]]"},
    ],
}

PARALLEL_OPTIONS = {
    "title": "PARALLEL PROCESSING OPTIONS",
    "options": [
        {"names": ("-t", "--threads"), "value": "int", "read": "threads", "field": "n_blast", "what": "BLAST threads",
         "help": f"Number of parallel sequence search threads [Default = half the CPU cores available, "
                 f"here [deep_sky_blue2]{n_threads_default}[/deep_sky_blue2] of [deep_sky_blue2]{n_cores}[/deep_sky_blue2]]"},
        {"names": ("-a", "--algthreads"), "value": "int", "read": "threads", "field": "n_process_alg",
         "what": "OrthoFinder threads",
         "help": "Number of parallel analysis threads [Default = a quarter of [dark_cyan]-t[/dark_cyan], "
                 "at least [deep_sky_blue2]1[/deep_sky_blue2], "
                 f"here [deep_sky_blue2]{max(1, n_threads_default // 4)}[/deep_sky_blue2]]"},
    ],
}

OUTPUT_OPTIONS = {
    "title": "OUTPUT OPTIONS",
    "options": [
        {"names": ("-n", "--name"), "value": "txt", "read": "name", "help": "Name to append to the results directory"},
        {"names": ("-o", "--output"), "value": "txt", "read": "output", "help": "Non-default results directory"},
        {"names": ("-X",), "set": {"q_add_species_to_ids": False}, "help": "Gene names without their species prefix in gene trees, alignments and Duplications.tsv"},
        {"names": ("-pof", "--pairwise-ortholog-files"), "set": {"save_space": False},
         "help": "Save the pairwise ortholog files. With many species (a file for each pair: n x (n-1) files), "
                 "it is better to run without it and then make only the files needed with the "
                 "[dark_cyan]split_ortholog_files[/dark_cyan] tool ([dark_cyan]tools/split_ortholog_files.py[/dark_cyan] "
                 "in the source code)"},
        {"names": ("--orthoxml", "-x"), "set": {"orthoxml": True},
         "help": "The orthogroups in OrthoXML format, for other tools, with their HOGs nested in them: " +
                 _path("Orthogroups/Orthogroups.orthoxml")},
        {"names": ("-mog", "--mcl-orthogroups"), "set": {"mcl_orthogroups": True},
         "help": "Also keep the [dark_cyan]MCL[/dark_cyan] orthogroups (before the HOGs): " +
                 _path("Orthogroups/MCL_Orthogroups.tsv") + " etc. (with [dark_cyan]--orthoxml[/dark_cyan] also " +
                 _path("MCL_Orthogroups.orthoxml") + "). They are always written with a stop option "
                 "([dark_cyan]-og[/dark_cyan], [dark_cyan]-ogt[/dark_cyan], [dark_cyan]-ost[/dark_cyan]) "
                 "or if the run fails"},
    ],
}

ADDITIONAL_OUTPUT_OPTIONS = {
    "title": "ADDITIONAL OUTPUT OPTIONS",
    "options": [
        {"names": ("--tree-view",), "preview": True, "set": {"tree_view": True},
         "help": "A page to view the trees in a web browser: each resolved gene tree with its "
                 "alignment (clades can be collapsed, a tree saved as SVG), a species tree explorer "
                 "(the HOGs at each node, an orthogroup's duplications and losses, each species' "
                 "genes) and all genes on the species tree: " + _path("Resolved_Gene_Trees/Resolved_Gene_Trees.html")},
        {"names": ("--stats-view",), "preview": True, "set": {"stats_view": True},
         "help": "A page of the run's statistics with interactive plots (overall, species, "
                 "orthogroups, duplications, HOGs): " + _path("Comparative_Genomics_Statistics/Statistics.html")},
        {"names": ("--groups-view",), "preview": True, "set": {"groups_view": True},
         "help": "A page of the genes of each orthogroup and HOG, searchable: " + _path("Orthogroups/Groups.html")},
        {"names": ("--view",), "preview": True, "set": {"tree_view": True, "stats_view": True, "groups_view": True},
         "help": "All three pages ([dark_cyan]--tree-view[/dark_cyan], [dark_cyan]--stats-view[/dark_cyan], "
                 "[dark_cyan]--groups-view[/dark_cyan]), linked to each other, and " + _path("OrthoFinder_Results.html") + " to "
                 "open them, a tab for each (also with two of them). With many species they can be large: each "
                 "can also be written on its own. For the results of a finished run: "
                 "[dark_cyan]view_results[/dark_cyan] <[bright_magenta]dir[/bright_magenta]>"},
    ],
}

STOPPING_OPTIONS = {
    "title": "WORKFLOW STOPPING OPTIONS",
    "options": [
        {"names": ("-op", "--only-prepare"), "set": {"q_stop_after_prepare": True},
         "help": "Stop after preparing input files for [dark_cyan]BLAST[/dark_cyan]"},
        {"names": ("-og", "--only-groups"), "set": {"q_stop_after_mcl_groups": True}, "help": "Stop after the [dark_cyan]MCL[/dark_cyan] orthogroups (no gene trees, species tree or HOGs)"},
        {"names": ("-ogt", "--only-trees"), "set": {"q_stop_after_trees": True}, "help": "Stop after gene trees"},
        {"names": ("-ost", "--only-species-trees"), "set": {"q_stop_after_species_trees": True},
         "help": "Stop after the rooted species tree (rooted by [dark_cyan]STRIDE[/dark_cyan] from the gene trees; no HOGs or orthologs)"},
    ],
}

RESTART_OPTIONS = {
    "title": "WORKFLOW RESTART COMMANDS",
    "options": [
        {"names": ("-b", "--blast"), "value": "dir", "read": "start", "field": "q_start_from_blast",
         "help": _dir(f"Restart {OF} from pre-computed sequence search results in <dir>")},
        {"names": ("-fg", "--from-groups"), "value": "dir", "read": "start", "field": "q_start_from_groups",
         "help": _dir(f"Restart {OF} from [dark_cyan]MCL[/dark_cyan] orthogroups results in <dir>")},
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
        {"names": ("--matrix",), "preview": True, "value": "txt", "read": "scoring",
         "help": "Scoring matrix of the [dark_cyan]DIAMOND[/dark_cyan] searches [Default = [dark_cyan]BLOSUM62[/dark_cyan]]",
         "choices": lambda pc: list(diamond_scoring.DIAMOND_MATRICES)},
        {"names": ("--custom-matrix",), "preview": True, "value": "file", "read": "scoring",
         "help": "A scoring matrix file for the [dark_cyan]DIAMOND[/dark_cyan] searches, requires [dark_cyan]--gapopen[/dark_cyan] and "
                 "[dark_cyan]--gapextend[/dark_cyan]\n"
                 "In NCBI's format (e.g. the files at " + _link("https://ftp.ncbi.nlm.nih.gov/blast/matrices/") + "): a line of "
                 "one-letter amino acid codes, then for each code, in that order, a line of the code and its integer "
                 "scores. The 20 standard amino acids are needed; lines starting with # are comments."},
        {"names": ("-go", "--gapopen"), "preview": True, "value": "int", "read": "scoring",
         "help": "Gap open penalty of the [dark_cyan]DIAMOND[/dark_cyan] searches, [dark_cyan]--gapopen[/dark_cyan] "
                 "[Default = the matrix's, e.g. [deep_sky_blue2]11[/deep_sky_blue2] for BLOSUM62]"},
        {"names": ("-ge", "--gapextend"), "preview": True, "value": "int", "read": "scoring",
         "help": "Gap extend penalty of the [dark_cyan]DIAMOND[/dark_cyan] searches, [dark_cyan]--gapextend[/dark_cyan] "
                 "[Default = the matrix's, e.g. [deep_sky_blue2]1[/deep_sky_blue2] for BLOSUM62]"},
        {"names": ("-h", "--help"), "read": "help", "help": "Print this help text"},
    ],
}

DEVELOPMENT_OPTIONS = {
    "title": "DEVELOPMENT OPTIONS",
    "hidden": True,   # not in the help (-h)
    "options": [
        {"names": ("-ms", "--min-seq"), "value": "int", "read": "min_seq",
         "help": "Minimum number of genes in an orthogroup for a gene tree "
                 "[Default = [deep_sky_blue2]4[/deep_sky_blue2]]"},
        {"names": ("--msa-min-seq",), "value": "int", "read": "msa_min_seq",
         "help": "Minimum number of sequences in an alignment for it to be trimmed "
                 "[Default = [deep_sky_blue2]10[/deep_sky_blue2]]"},
        {"names": ("-nk",), "value": "int", "read": "int", "field": "n_skip", "what": "-nk",
         "help": "With [dark_cyan]-ST[/dark_cyan], leave out the gene trees of the largest <int> orthogroups "
                 "[Default = [deep_sky_blue2]50[/deep_sky_blue2]]"},
        {"names": ("--scores-v2",), "set": {"v2_scores": True},
         "help": "Normalise the sequence search scores as normalised bit scores (v2)"},
        {"names": ("-1",), "set": {"q_double_blast": False},
         "help": "Search each pair of species in one direction only"},
        {"names": ("-mt", "--method-threads"), "value": "int", "read": "method_threads",
         "help": "Threads for each alignment and tree command that takes a thread count "
                 "[Default = [deep_sky_blue2]1[/deep_sky_blue2]]"},
        {"names": ("-fd",), "value": "txt", "read": "fd_limit",
         "help": "The soft and hard limits on the number of open files, e.g. "
                 "[deep_sky_blue2]10000,20000[/deep_sky_blue2]"},
        {"names": ("-efn", "--extended-filename"), "set": {"extended_filename": True},
         "help": "Add the scoring and the search, MSA and tree programs to the results directory name"},
        {"names": ("--keep-gene-trees",), "set": {"keep_gene_trees": True},
         "help": "Keep the gene trees in the working directory (removed by default)"},
        {"names": ("--keep-resolved-gene-trees",), "set": {"keep_resolved_gene_trees": True},
         "help": "Keep the resolved gene trees as separate files in the working directory (removed by default)"},
        {"names": ("--keep-legacy",), "set": {"keep_legacy": True},
         "help": "Write the legacy files in WorkingDirectory/Legacy (not written by default): the HOGs from "
                 "before the orthogroups are updated and the orthogroups as sequence IDs; and keep the "
                 "working files the update uses"},
        {"names": ("--no-print-info",), "set": {"print_info": False},
         "help": "Don't print the progress messages of the ortholog inference"},
        {"names": ("--concat-alignment",), "set": {"species_tree_alignment": True},
         "help": "With [dark_cyan]-ST[/dark_cyan], also write the concatenated alignment of the orthogroups "
                 "for the species tree: " + _path("WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa")},
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
    DEVELOPMENT_OPTIONS,
    VERSION_OPTIONS,
)

# Each option by each of its names
OPTIONS_BY_NAME = {name: option for group in OPTION_GROUPS for option in group["options"] for name in option["names"]}


def is_available(option):
    """Whether an option can be used and is in the help: not a preview one, unless PREVIEW_OPTIONS."""
    return PREVIEW_OPTIONS or not option.get("preview")


def help_label(option):
    """An option as the help shows it, e.g. "-t <int>"."""
    if "value" not in option:
        return option["names"][0]
    return "%s <[bright_magenta]%s[/bright_magenta]>" % (option["names"][0], option["value"])
