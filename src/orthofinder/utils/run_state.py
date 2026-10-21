"""
The state of a results directory: where the files of its analysis are, and
what a restart needs to continue it exactly as it would have gone on. Kept in
<results directory>/WorkingDirectory/run_state.json, written as the run goes
on, and read by the restart options (-b, -fg, -fgt, -fst). Without it (e.g.
results directories of older versions) they read Log.txt, which lists the
working directories too, and for older versions previous_wd.txt, which those
wrote in each working directory.

Log.txt is the record of a run for people; checkpoint.txt the record of its
progress (stages, finished searches); run_state.json the locations and
settings a restart needs.

  run         "analysis": "core" (-f) or "assign" (--assign); "invocations": the
              OrthoFinder commands that made and continued it, each
              {"number", "command", "directory" (where it was run: its
              paths are relative to it), "started", "version"} (results of older
              versions: "commands", the command lines only).
  commands    one example of each external command the run used, for
              debugging: {"<program> [<sub-command>]": {"command" (as it was
              run: the program's full path, the threads filled in), "program",
              "stage", "invocation" (its number)}}, the first of each in the
              latest invocation that ran it (a dependency test's only if the
              run ran none)
              (note_command; "(species tree)": the tree of the concatenated
              alignment).
  dirs        "working": this results directory's WorkingDirectory (the
              current dir, where the run writes);
              "base": the working directories read from, this one first
              (species and sequence IDs, species FASTA, search results), then
              those it extends (also in Log.txt: WorkingDirectory_Core
              or WorkingDirectory_Previous);
              "trees": the working directory with the gene and species trees;
              --assign: "core_results", "core_working": the results and
              working directory of the core analysis (-f) the new species were
              assigned to. The base can differ from the core (it starts with
              this --assign working directory).
  files       "orthogroups": the orthogroups (clusters) file;
              "species_tree": the species tree the run uses, a copy in its
              working directory (one given with -s, or for --assign the one
              inferred to find the clades of new species).
  assign      "new_species": the directory of the new species' FASTA files;
              "first_new_species": the ID of the first new species;
              "i_og_restart": the first orthogroup whose gene tree is still to
              be inferred (the others were inferred to find the clades).
  programs    {program: full path} of the programs the run uses (the search,
              MCL, MSA, tree and species tree programs, as the dependency
              check ran them; parallel_task_manager.program_dirs: where they
              are looked for). A restart's check updates them.
  python      the Python the run uses: "executable", "prefix", "version", "orthofinder"
              (its version and path) and "packages": {name: {"version", "path"}}
              of the Python packages it needs (null if one is not installed).
  system      the computer the run is on (utils.system_info: OS, glibc, CPU and
              the cores available to the run, memory and a container's or job's
              limit, open-file limits, container, cluster job), for debugging;
              null where a value cannot be read. A restart's check updates it.
  settings    the settings a restart continues with (SETTINGS), unless given
              again on its command line. A custom scoring matrix
              (--custom-matrix) is the copy in the working directory that
              the searches used, so a restart needs neither the file given
              nor --custom-matrix again.

Paths inside the results directory are kept relative to it, so a results
directory that is copied or moved keeps working; others are absolute.
"""
import contextlib
import datetime
import json
import os
import re
import shutil
import threading

from ..search import diamond_scoring

VERSION = 1
STATE_FN = "run_state.json"
WORKING_DIR_NAME = "WorkingDirectory"
SPECIES_TREE_FN = "run_state_species_tree.txt"


# ------------------------------------------------------------------ paths ----

def _results_dir(results_dir):
    return os.path.join(os.path.abspath(results_dir), "")


def state_path(results_dir):
    return os.path.join(_results_dir(results_dir), WORKING_DIR_NAME, STATE_FN)


def _record(results_dir, path):
    """A path as kept in the state: relative if inside results_dir (keeping a trailing separator)."""
    if path is None:
        return None
    trailing = path.endswith(os.sep)
    full = os.path.abspath(path)
    rd = _results_dir(results_dir)
    if (full + os.sep).startswith(rd):
        rel = os.path.relpath(full, rd)
        rel = "" if rel == "." else rel
    else:
        rel = full
    return rel + os.sep if trailing and rel else rel


def _resolve(results_dir, recorded):
    """The absolute path of a path kept in the state."""
    if recorded is None:
        return None
    trailing = recorded.endswith(os.sep) or recorded == ""
    full = recorded if os.path.isabs(recorded) else os.path.join(_results_dir(results_dir), recorded)
    full = os.path.abspath(full)
    return full + os.sep if trailing else full


# ------------------------------------------------------------ read/write ----

def find(continuation_dir):
    """
    The results directory of a directory given to a restart option (the
    results directory, or its WorkingDirectory) and its state: (results_dir,
    state), or (None, {}) if it has no run_state.json.
    """
    d = os.path.abspath(continuation_dir)
    for results_dir in (d, os.path.dirname(d)):
        if os.path.exists(state_path(results_dir)):
            return _results_dir(results_dir), read(results_dir)
    return None, {}


def read(results_dir):
    """The state of a results directory ({} if none)."""
    fn = state_path(results_dir)
    if not os.path.exists(fn):
        return {}
    with open(fn, encoding="utf-8") as infile:
        return json.load(infile)


def _write(results_dir, state):
    fn = state_path(results_dir)
    tmp = fn + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as outfile:
        json.dump(state, outfile, indent=2, sort_keys=True)
    os.replace(tmp, fn)        # never a half-written state


# One update at a time: the commands are noted from the threads, and the forked
# processes, that run them. (A lock of the working directory itself, not a file.)
_update_lock = threading.RLock()
_lock_depth = [0]


@contextlib.contextmanager
def _locked(results_dir):
    with _update_lock:
        fd = None
        if _lock_depth[0] == 0:
            try:
                import fcntl
                fd = os.open(os.path.join(_results_dir(results_dir), WORKING_DIR_NAME), os.O_RDONLY)
                fcntl.flock(fd, fcntl.LOCK_EX)
            except (ImportError, OSError):     # (Windows; or not lockable: the threads only)
                if fd is not None:
                    os.close(fd)
                fd = None
        _lock_depth[0] += 1
        try:
            yield
        finally:
            _lock_depth[0] -= 1
            if fd is not None:
                os.close(fd)                    # (releases the lock)


def update(results_dir, section, **values):
    """Set values in a section of the state (paths are given as they are: see record_dirs etc.)."""
    if not os.path.isdir(os.path.join(_results_dir(results_dir), WORKING_DIR_NAME)):
        return
    with _locked(results_dir):
        state = read(results_dir)
        state["version"] = VERSION
        state.setdefault(section, {}).update(values)
        _write(results_dir, state)


# ------------------------------------------------------------- the parts ----

def record_run(results_dir, analysis, command):
    """
    The analysis a results directory holds (set by the run that made it) and the
    OrthoFinder commands run on it (the latest 20): each with the directory it was
    run in, when it started and the version that ran it.
    """
    from .. import __version__
    if not os.path.isdir(os.path.join(_results_dir(results_dir), WORKING_DIR_NAME)):
        return
    with _locked(results_dir):
        state = read(results_dir)
        run = state.get("run", {})
        invocations = run.get("invocations", [])
        _invocation[0] = invocations[-1].get("number", len(invocations)) + 1 if invocations else 1
        invocations.append({"number": _invocation[0], "command": command, "directory": os.getcwd(),
                            "started": datetime.datetime.now().isoformat(timespec="seconds"),
                            "version": __version__})
        update(results_dir, "run", analysis=run.get("analysis", analysis), invocations=invocations[-20:])


# ------------------------------------------------- the external commands ----
# The stage of the run (the workflow's steps, utils.logging.RunLogger.stage and
# run.context.RunContext.step), for the commands noted while it runs
_stage = ["Initialisation"]


@contextlib.contextmanager
def stage(name):
    """The commands run in this block are noted as of the stage name."""
    previous = _stage[0]
    _stage[0] = name
    try:
        yield
    finally:
        _stage[0] = previous


_noted = set()            # the keys noted in this invocation (by this process)
_tested = set()           # those noted by a dependency test only (a run's command replaces it)
DEPENDENCY_CHECK = "Dependency check"
_invocation = [0]         # its number (record_run)
_noted_lock = threading.Lock()
_SUB_COMMAND = re.compile(r"[A-Za-z][A-Za-z0-9_-]*$")
_ASSIGNMENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*=")
# (the concatenated alignment of the species tree: file_handler.get_species_tree_concat_align_fn)
_SPECIES_TREE_ALIGNMENT = "SpeciesTreeAlignment.fa"


def command_key(command):
    """
    "<program> [<sub-command>]" of an external command: its program, and
    the word after it if a sub-command (diamond makedb / blastp, mmseqs createdb /
    search: not an option or a file); "(species tree)" if it reads the concatenated
    alignment of the species tree. None if it runs no program.
    """
    first = re.split(r"\n|;|&&|\|\||\|", command, maxsplit=1)[0].split()
    while first and _ASSIGNMENT.match(first[0]):
        first.pop(0)
    if not first:
        return None
    key = os.path.basename(first[0])
    if len(first) > 1 and _SUB_COMMAND.match(first[1]):
        key += " " + first[1]
    if _SPECIES_TREE_ALIGNMENT in command:
        key += " (species tree)"
    return key


def note_command(command):
    """
    Keep one example of each kind of external command (command_key) in the run's
    state ("commands"): the first the invocation runs (by any of its processes),
    as it is run (the program's full path, utils.programs; the threads filled in),
    with its stage. A dependency test's command is kept only until the run's
    first command of its kind, and one of an earlier invocation is replaced.
    Never fails the run.
    """
    try:
        key = command_key(command)
        test = _stage[0] == DEPENDENCY_CHECK
        if key is None or key in _noted or (test and key in _tested):
            return
        with _noted_lock:
            if key in _noted or (test and key in _tested):
                return
            (_tested if test else _noted).add(key)
        from ..utils import file_handler, programs
        results_dir = getattr(file_handler.handler, "rd1", None)
        if not results_dir:
            return
        program = programs.command_programs(command)
        path = None
        if program:
            # (a dependency test, before a copy is chosen, runs the first on PATH)
            from . import parallel_task_manager
            path = programs.chosen().get(program[0]) or shutil.which(
                program[0], path=parallel_task_manager.my_env.get("PATH"))
        if not os.path.isdir(os.path.join(_results_dir(results_dir), WORKING_DIR_NAME)):
            return
        with _locked(results_dir):
            noted = read(results_dir).get("commands", {}).get(key)
            if (noted is not None and noted.get("invocation") == _invocation[0]
                    and (test or noted.get("stage") != DEPENDENCY_CHECK)):
                return             # (noted by another process)
            update(results_dir, "commands", **{key: {"command": programs.with_paths(command), "program": path,
                                                     "stage": _stage[0], "invocation": _invocation[0]}})
    except Exception:
        pass


def record_dirs(results_dir, working=None, base=None, trees=None, core_results=None, core_working=None):
    """Directories of the analysis (None: unchanged)."""
    values = {}
    if working is not None:
        values["working"] = _record(results_dir, working)
    if base is not None:
        values["base"] = [_record(results_dir, d) for d in base]
    if trees is not None:
        values["trees"] = _record(results_dir, trees)
    if core_results is not None:
        values["core_results"] = _record(results_dir, core_results)
    if core_working is not None:
        values["core_working"] = _record(results_dir, core_working)
    if values:
        update(results_dir, "dirs", **values)


def record_file(results_dir, name, path):
    update(results_dir, "files", **{name: _record(results_dir, path)})


def dirs(results_dir, state):
    """The directories in a state, as absolute paths ending in the separator (base: a list)."""
    def as_dir(recorded):
        d = _resolve(results_dir, recorded)
        return None if d is None else os.path.join(d, "")
    out = {}
    for name, value in state.get("dirs", {}).items():
        out[name] = [as_dir(d) for d in value] if isinstance(value, list) else as_dir(value)
    return out


def file(results_dir, state, name):
    """A file in a state, as an absolute path (None if not recorded or missing)."""
    fn = _resolve(results_dir, state.get("files", {}).get(name))
    return fn if fn is not None and os.path.exists(fn) else None


def species_tree_copy_fn(results_dir):
    """The copy of the species tree the run uses, in its working directory (keep_species_tree)."""
    return os.path.join(_results_dir(results_dir), WORKING_DIR_NAME, SPECIES_TREE_FN)


def keep_species_tree(results_dir, species_tree_fn):
    """Keep a copy of the species tree the run uses, in its working directory (a restart uses it)."""
    copy_fn = species_tree_copy_fn(results_dir)
    if os.path.abspath(species_tree_fn) != os.path.abspath(copy_fn):
        shutil.copyfile(species_tree_fn, copy_fn)
    record_file(results_dir, "species_tree", copy_fn)
    return copy_fn


def record_assign(results_dir, core_results, new_species, first_new_species):
    update(results_dir, "assign", new_species=_record(results_dir, new_species),
           first_new_species=first_new_species)
    record_dirs(results_dir, core_results=core_results)


def assign_info(results_dir, state):
    """The --assign part of a state, with absolute paths ({} if not an --assign analysis)."""
    info = dict(state.get("assign", {}))
    if "new_species" in info:
        info["new_species"] = _resolve(results_dir, info["new_species"])
    return info


# The scoring of DIAMOND searches is one setting: if any of these options is
# given again, none of the matrix and gap penalties is restored.
SCORING_FLAGS = ("--matrix", "--custom-matrix", "-go", "--gapopen", "-ge", "--gapextend")
SCORING_SETTINGS = ("score_matrix", "gapopen", "gapextend")

# The settings a restart continues with (a run's, unless given again): name
# in options, and the command-line options that set it.
SETTINGS = (
    ("search_program", ("-S", "--search")),
    ("score_matrix", SCORING_FLAGS),
    ("gapopen", SCORING_FLAGS),
    ("gapextend", SCORING_FLAGS),
    ("dna", ("-d", "--dna")),
    ("q_double_blast", ("-1",)),
    ("v2_scores", ("--scores-v2",)),
    ("mcl_inflation", ("-I", "--inflation")),
    ("q_msa_trees", ("-M", "--method")),
    ("msa_program", ("-A", "--msa-program")),
    ("tree_program", ("-T", "--tree-program")),
    ("min_seq", ("-ms", "--min-seq")),
    ("msa_min_seq", ("--msa-min-seq",)),
    ("q_trim", ("-z",)),
    ("species_tree_program", ("-ST", "--species-tree-program")),
    ("n_skip", ("-nk",)),
    ("q_split_para_clades", ("-y",)),
    ("q_add_species_to_ids", ("-X",)),
)


def _is_path(name, value):
    """Whether a setting is a path: a custom scoring matrix (a file, not one of DIAMOND's)."""
    return name == "score_matrix" and isinstance(value, str) and diamond_scoring.is_custom_matrix(value)


def save_settings(results_dir, options):
    """Keep the settings of this run (see SETTINGS)."""
    settings = {name: getattr(options, name) for name, _ in SETTINGS}
    settings = {name: _record(results_dir, v) if _is_path(name, v) else v for name, v in settings.items()}
    update(results_dir, "settings", **settings)


def _saved_settings(state, results_dir):
    """The settings kept in a state, with paths absolute (results_dir: where the state is; None: as kept)."""
    saved = dict(state.get("settings", {}))
    if results_dir is not None:
        saved = {name: _resolve(results_dir, v) if _is_path(name, v) else v for name, v in saved.items()}
    return saved


def restore_settings(state, options, args, names=None, results_dir=None):
    """
    Set the options to those of the run being continued, except those given
    on this command line (args), or only those in names. results_dir: the
    results directory of the state (find), for its paths. Returns
    [(name, value)] of those changed.
    """
    given = set(args)
    saved = _saved_settings(state, results_dir)
    if "astral" in saved and "species_tree_program" not in saved:
        # runs before -ST named a configured program: ASTRAL-Pro or the default
        saved["species_tree_program"] = "astral-pro" if saved["astral"] else None
    changed = []
    for name, flags in SETTINGS:
        if names is not None and name not in names:
            continue
        if name in saved and not given.intersection(flags) and getattr(options, name) != saved[name]:
            setattr(options, name, saved[name])
            changed.append((name, saved[name]))
    return changed


def scoring_differences(state, options, results_dir=None):
    """
    The scoring settings given that differ from those of the run being
    continued, as [(name, its value, this run's value)] ([] if it has none).
    """
    saved = _saved_settings(state, results_dir)
    return [(name, saved[name], getattr(options, name)) for name in SCORING_SETTINGS
            if name in saved and str(saved[name]) != str(getattr(options, name))]
