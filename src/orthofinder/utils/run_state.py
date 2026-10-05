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

  run         "analysis": "core" (-f) or "assign" (--assign); the commands
              that made and continued it.
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
  settings    the settings a restart continues with (SETTINGS), unless given
              again on its command line.

Paths inside the results directory are kept relative to it, so a results
directory that is copied or moved keeps working; others are absolute.
"""
import json
import os
import shutil

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
    with open(fn) as infile:
        return json.load(infile)


def _write(results_dir, state):
    fn = state_path(results_dir)
    tmp = fn + ".tmp"
    with open(tmp, "w") as outfile:
        json.dump(state, outfile, indent=2, sort_keys=True)
    os.replace(tmp, fn)        # never a half-written state


def update(results_dir, section, **values):
    """Set values in a section of the state (paths are given as they are: see record_dirs etc.)."""
    if not os.path.isdir(os.path.join(_results_dir(results_dir), WORKING_DIR_NAME)):
        return
    state = read(results_dir)
    state["version"] = VERSION
    state.setdefault(section, {}).update(values)
    _write(results_dir, state)


# ------------------------------------------------------------- the parts ----

def record_run(results_dir, analysis, command):
    """The analysis a results directory holds (set by the run that made it) and the commands run on it."""
    state = read(results_dir)
    run = state.get("run", {})
    commands = run.get("commands", [])
    commands.append(command)
    update(results_dir, "run", analysis=run.get("analysis", analysis), commands=commands[-20:])


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
    """The directories in a state, as absolute paths (base: a list)."""
    out = {}
    for name, value in state.get("dirs", {}).items():
        out[name] = ([_resolve(results_dir, d) for d in value] if isinstance(value, list)
                     else _resolve(results_dir, value))
    return out


def file(results_dir, state, name):
    """A file in a state, as an absolute path (None if not recorded or missing)."""
    fn = _resolve(results_dir, state.get("files", {}).get(name))
    return fn if fn is not None and os.path.exists(fn) else None


def keep_species_tree(results_dir, species_tree_fn):
    """Keep a copy of the species tree the run uses, in its working directory (a restart uses it)."""
    copy_fn = os.path.join(_results_dir(results_dir), WORKING_DIR_NAME, SPECIES_TREE_FN)
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


# The settings a restart continues with (a run's, unless given again): name
# in options, and the command-line options that set it.
SETTINGS = (
    ("search_program", ("-S", "--search")),
    ("dna", ("-d", "--dna")),
    ("qDoubleBlast", ("-1",)),
    ("v2_scores", ("--scores-v2",)),
    ("mclInflation", ("-I", "--inflation")),
    ("qMSATrees", ("-M", "--method")),
    ("msa_program", ("-A", "--msa-program")),
    ("tree_program", ("-T", "--tree-program")),
    ("min_seq", ("-ms", "--min-seq")),
    ("msa_min_seq", ("--msa-min-seq",)),
    ("qTrim", ("-z",)),
    ("astral", ("-ST", "--species-tree-program")),
    ("n_skip", ("-nk",)),
    ("qSplitParaClades", ("-y",)),
    ("qAddSpeciesToIDs", ("-X",)),
)


def save_settings(results_dir, options):
    """Keep the settings of this run (see SETTINGS)."""
    update(results_dir, "settings", **{name: getattr(options, name) for name, _ in SETTINGS})


def restore_settings(state, options, args):
    """
    Set the options to those of the run being continued, except those given
    on this command line (args). Returns [(name, value)] of those changed.
    """
    given = set(args)
    saved = state.get("settings", {})
    changed = []
    for name, flags in SETTINGS:
        if name in saved and not given.intersection(flags) and getattr(options, name) != saved[name]:
            setattr(options, name, saved[name])
            changed.append((name, saved[name]))
    return changed
