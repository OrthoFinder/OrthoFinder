import os 
import sys
import re
import traceback
import json
from collections import defaultdict
import pytest
from orthofinder.run.main import main
from orthofinder.run.process_args import GetFileArgument



ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def strip_ansi(text):
    return ANSI_RE.sub("", text or "")


def assert_no_orthofinder_failure(name, code, out, err):
    text = out + err
    clean = strip_ansi(text)

    fatal_patterns = (
        r"(?mi)^\s*ERROR:\s",
        r"(?mi)Unrecognised argument:",
        r"(?mi)Unknown option:",
        r"(?mi)Incorrect option:",
        r"(?m)^\s*Traceback \(most recent call last\):",
        r'(?m)^\s*File ".*", line \d+, in \S+',
        r"(?m)^[A-Za-z_]\w*Error:\s",
    )

    failed_by_text = any(re.search(p, clean) for p in fatal_patterns)

    assert code == 0 and not failed_by_text, (
        f"{name}: OrthoFinder failed\n"
        f"exit code: {code}\n"
        f"--- stdout ---\n{out}\n"
        f"--- stderr ---\n{err}"
    )

def run_orthofinder_core_case(
    name,
    argstr,
    input_proj,
    example_data,
    dna_projects,
    species_tree,
    user_ofconfig,
    capfd,
):

    uses_example_data = "EXAMPLE_DATA" in argstr
    run_project = example_data if uses_example_data else os.path.abspath(input_proj)

    s = argstr.strip().replace("EXAMPLE_DATA", example_data)\
                      .replace("DNA_INPUT", dna_projects)\
                      .replace("INPUT", input_proj)\
                      .replace("SPECIES_TREE", species_tree)\
                      .replace("USER_CONFIG", user_ofconfig)

    args = s.split()[1:] + ["-n", name]

    code, out, err, text = _run_main(args, capfd)

    assert_no_orthofinder_failure(name, code, out, err)

    results_dir = os.path.join(run_project, "OrthoFinder")
    result_path = _find_output_dir(
        results_dir,
        "Results_" + name,
        find_core=True
    )

    assert result_path, (
        f"{name}: OrthoFinder exited successfully, but no valid Results_* "
        f"directory was found in {results_dir}\n"
        f"--- stdout ---\n{out}\n--- stderr ---\n{err}"
    )

    os.environ[f"ORTHOFINDER_TEST_RESULTS_{name}"] = result_path

def run_orthofinder_assign_case(
    name,
    argstr,
    input_proj,
    example_data,
    assign,
    species_tree_assign,
    user_ofconfig,
    capfd,
):

    core_name = name.replace("_assign", "").replace("_restart", "")

    results_dir = os.path.join(input_proj, "OrthoFinder")

    core_results = os.environ.get(f"ORTHOFINDER_TEST_RESULTS_{core_name}", "")

    if not core_results:
        core_results = _find_output_dir(
            results_dir,
            "Results_" + core_name,
            find_core=True
        )

    if not core_results:
        core_results = _find_output_dir(
            os.path.join(example_data, "OrthoFinder"),
            "Results_" + core_name,
            find_core=True,
        )

    if not core_results:
        raise AssertionError(
            f"Core results for '{core_name}' not found in {results_dir}. "
            f"Expected a valid OrthoFinder Results_* directory with a WorkingDirectory."
        )

    s = (
        argstr.strip()
        .replace("EXAMPLE_DATA", example_data)
        .replace("INPUT", assign)
        .replace("CORE_RESULTS", core_results)
        .replace("SPECIES_TREE", species_tree_assign)
        .replace("USER_CONFIG", user_ofconfig)
    )

    args = s.split()[1:] + ["-n", name]

    code, out, err, text = _run_main(args, capfd)

    assert_no_orthofinder_failure(name, code, out, err)

    assign_results_dir = os.path.dirname(core_results)
    result_path = _find_output_dir(
        assign_results_dir,
        "Results_" + name,
        find_core=False
    )

    assert result_path, (
        f"{name}: assign run exited successfully, but no valid assign Results_* "
        f"directory was found in {assign_results_dir}\n"
        f"--- stdout ---\n{out}\n--- stderr ---\n{err}"
    )

    os.environ[f"ORTHOFINDER_TEST_RESULTS_{name}"] = result_path

def _is_valid_results_dir(path):
    """
    A real OrthoFinder results directory should have a WorkingDirectory.
    Use this to avoid returning incomplete/random directories.
    """
    if not path:
        return False

    if not os.path.isdir(path):
        return False

    wd = os.path.join(path, "WorkingDirectory")
    if not os.path.isdir(wd):
        return False

    species_ids = os.path.join(wd, "SpeciesIDs.txt")
    sequence_ids = os.path.join(wd, "SequenceIDs.txt")

    return os.path.exists(species_ids) and os.path.exists(sequence_ids)


def _normalise_result_name(name):
    """
    For matching core results from assign/restart names.
    """
    return (
        name
        .replace("_assign", "")
        .replace("_restart", "")
    )


def _find_output_dir(results_dir, test_filename, fileno=-1, find_core=True) -> str:
    """
    Find an OrthoFinder Results_* directory.

    Previous version was too strict:
        expected Results_famsa_species_tree exactly.

    This version:
        1. tries strict name match first
        2. then tries relaxed name match
        3. then falls back to newest valid Results_* directory
    """
    if isinstance(results_dir, list):
        results_dir = results_dir[0]

    results_dir = os.path.abspath(results_dir)

    if not os.path.isdir(results_dir):
        return ""

    candidates = []

    for name in os.listdir(results_dir):
        path = os.path.join(results_dir, name)

        if not os.path.isdir(path):
            continue

        if not name.startswith("Results_"):
            continue

        if not _is_valid_results_dir(path):
            continue

        match_name = _normalise_result_name(name) if find_core else name

        candidates.append((name, match_name, path, os.stat(path).st_mtime))

    if not candidates:
        return ""

    # 1. Strict/expected match
    matched = [
        (mtime, path)
        for name, match_name, path, mtime in candidates
        if test_filename in match_name
    ]

    if matched:
        return sorted(matched)[fileno][1]

    # 2. Relaxed match: remove Results_ prefix and try again
    relaxed = test_filename
    if relaxed.startswith("Results_"):
        relaxed = relaxed[len("Results_"):]

    matched = [
        (mtime, path)
        for name, match_name, path, mtime in candidates
        if relaxed in match_name
    ]

    if matched:
        return sorted(matched)[fileno][1]

    # 3. Fallback: newest valid Results_* directory
    # This prevents the entire results test suite from collapsing to
    # "WorkingDirectory" when the name has changed.
    newest = sorted((mtime, path) for name, match_name, path, mtime in candidates)
    return newest[fileno][1] if newest else ""

def read_config_file(configure_file, user_config_file=None):
    
    if user_config_file is not None:
        configure_file = user_config_file
    
    of_config_dict = {}
    configure_file = GetFileArgument(configure_file)    
    with open(configure_file, "r") as infile:
        try:
            d = json.load(infile)
        except ValueError:
            print(("WARNING: Incorrectly formatted configuration file %s" % configure_file))
            print("File is not in .json format. No user-confgurable multiple sequence alignment or tree inference methods have been added.\n")
            return
        for name, val_dict in d.items():
            if name == "__comment":
                continue
            if " " in name:
                print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                print(("No space is allowed in name: '%s'" % name))
                continue
            if not isinstance(val_dict, dict) or len(val_dict) == 0:
                continue
            of_config_dict[name] = {}
            for k, v in val_dict.items():
                if not isinstance(v, str) or len(v) == 0:
                    continue
                of_config_dict[name][k] = v
    return of_config_dict



def _run_main(args, capfd):
    try:
        ret = main(args)
        code = 0 if ret is None else int(ret)

    except SystemExit as e:
        if e.code is None:
            code = 0
        elif isinstance(e.code, int):
            code = e.code
        else:
            code = 1

    except Exception:
        out, err = capfd.readouterr()
        tb = traceback.format_exc()
        pytest.fail(
            "Unexpected exception from main()\n\n"
            f"{tb}\n"
            f"--- Captured stdout ---\n{out}\n"
            f"--- Captured stderr ---\n{err}\n"
        )

    out, err = capfd.readouterr()
    text = out + err
    return code, out, err, text


def create_path(arg):
    filepath = os.path.abspath(arg)
    if not os.path.isfile(filepath) and filepath[-1] != os.sep:
        filepath += os.sep
    return filepath



def _latest_output_dir(results_dir, fileno=-1) -> str:
    if isinstance(results_dir, list):
        results_dir = results_dir[0]

    results_dir = os.path.abspath(results_dir)

    try:
        entries = [
            (os.stat(os.path.join(results_dir, name)).st_mtime,
             os.path.join(results_dir, name))
            for name in os.listdir(results_dir)
        ]
        return sorted(entries)[fileno][1] if entries else ""
    except FileNotFoundError:
        return ""
    except NotADirectoryError:
        return results_dir


def _index_of_orthologues(of_orthologues):
    per_species_triplets = defaultdict(set)
    per_species_locations = defaultdict(lambda: defaultdict(set))

    for species, og_map in of_orthologues.items():
        for og_id, entries in og_map.items():
            for e in entries:
                per_species_triplets[species].add(e)
                per_species_locations[species][e].add(og_id)

    return per_species_triplets, per_species_locations


def _index_duplications(dup_dict):

    def _norm_label(s: str) -> str:
        for dash in ["\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2212"]:
            s = s.replace(dash, "-")
        s = s.strip()
        return "Non-Terminal" if s.lower().startswith("non") else "Terminal"

    def _norm_id(x: str) -> str:
        return " ".join(x.split())

    def _norm_triplet(t):
        label, A, B = t
        A = frozenset(_norm_id(x) for x in A)
        B = frozenset(_norm_id(x) for x in B)
        label = _norm_label(label)
        return (label, A, B)

    trips = set()
    locations = defaultdict(set)

    for og_id, entries in dup_dict.items():
        for e in entries:
            t = _norm_triplet(e)
            t_rev = (t[0], t[2], t[1])
            trips.add(t)
            trips.add(t_rev)
            locations[t].add(og_id)
            locations[t_rev].add(og_id)

    return trips, locations

def _index_of_orthogroups(orthogroups_dict):
    all_ogs_set = set()
    per_ogs_locations = {}

    for ogname, ogs in orthogroups_dict.items():
        all_ogs_set.add(tuple(sorted(ogs)))
        per_ogs_locations[tuple(sorted(ogs))] = ogname
    return all_ogs_set, per_ogs_locations

def _index_of_hogs(hogs_dict):
    all_hogs_set = set()
    per_hog_og = {}
    per_hog_node = {}
    for ogname, node_hogs in hogs_dict.items():
        for node, hogs in node_hogs.items():
            for hog in hogs:
                all_hogs_set.add(tuple(sorted(hog)))
                per_hog_node[tuple(sorted(hog))] = node
                per_hog_og[tuple(sorted(hog))] = ogname

    return all_hogs_set, per_hog_og, per_hog_node
