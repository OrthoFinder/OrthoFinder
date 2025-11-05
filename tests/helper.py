import os 
import sys
import re
import traceback
import json
from collections import defaultdict
import pytest
from orthofinder.run.__main__ import main
from orthofinder.run.process_args import GetFileArgument


def run_orthofinder_core_case(name, argstr, input_proj, dna_projects, species_tree, user_ofconfig, capfd):
    s = argstr.strip().replace("INPUT", input_proj)\
                      .replace("DNA_INPUT", dna_projects)\
                      .replace("SPECIES_TREE", species_tree)\
                      .replace("USER_CONFIG", user_ofconfig)
                      
    args = s.split()[1:] + ["-n", name]

    code, out, err, text = _run_main(args, capfd)
    assert code == 0, f"{name}: exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    fatal_patterns = (
        r"(?mi)^\s*ERROR:\s",
        r"(?m)^\s*Traceback \(most recent call last\):",
        r'(?m)^\s*File ".*", line \d+, in \S+',
        r"(?m)^[A-Za-z_]\w*Error:\s",
    )
    assert not any(re.search(p, text) for p in fatal_patterns), text


def run_orthofinder_assign_case(name, argstr, input_proj, assign, species_tree_assign, user_ofconfig, capfd):
    
    results_dir = os.path.join(input_proj, "OrthoFinder")
    projects_results = _find_output_dir(results_dir, "Results_" + name.rsplit("_", 1)[0])
    
    s = argstr.strip().replace("INPUT", assign)\
                      .replace("CORE_RESULTS", projects_results)\
                      .replace("SPECIES_TREE", species_tree_assign)\
                      .replace("USER_CONFIG", user_ofconfig)
                      
    args = s.split()[1:] + ["-n", name]

    code, out, err, text = _run_main(args, capfd)
    assert code == 0, f"{name}: exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    fatal_patterns = (
        r"(?mi)^\s*ERROR:\s",
        r"(?m)^\s*Traceback \(most recent call last\):",
        r'(?m)^\s*File ".*", line \d+, in \S+',
        r"(?m)^[A-Za-z_]\w*Error:\s",
    )
    assert not any(re.search(p, text) for p in fatal_patterns), text



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
    """
    Call main(args) in-process and return (exit_code, out, err, text).
    - Success: exit_code == 0
    - If main() returns None, we treat it as 0
    - If main() calls sys.exit(None), treat as 0
    """
    try:
        ret = main(args)
        code = 0 if ret is None else int(ret)
    except SystemExit as e:
        # sys.exit(0) or sys.exit(None) are both "success"
        code = e.code if isinstance(e.code, int) else 0
    except Exception:
        # Unexpected crash: fail with traceback + whatever was printed so far
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

def _find_output_dir(results_dir, test_filename, fileno=-1) -> str:
    if isinstance(results_dir, list):
        results_dir = results_dir[0]

    results_dir = os.path.abspath(results_dir)

    try:
        available_paths = [
            os.path.join(results_dir, name)
            for name in os.listdir(results_dir)
            if test_filename == name.split(".")[0].rsplit("_", 1)[0]
        ]
                
        entries = [
            (os.stat(path).st_mtime, path)
            for path in available_paths
            if path is not None 
        ]
        return sorted(entries)[fileno][1] if entries else ""
    except FileNotFoundError:
        return ""
    except NotADirectoryError:
        return results_dir

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
            trips.add(t)
            locations[t].add(og_id)

    return trips, locations

def _index_of_orthogroups(orthogroups_dict):
    all_ogs_set = set()
    per_ogs_locations = {}

    for ogname, ogs in orthogroups_dict.items():
        all_ogs_set.add(ogs)
        per_ogs_locations[ogs] = ogname
    return all_ogs_set, per_ogs_locations

def _index_of_hogs(hogs_dict):
    all_hogs_set = set()
    per_hog_og = {}
    per_hog_node = {}
    for ogname, node_hogs in hogs_dict.items():
        for node, hogs in node_hogs.items():
            for hog in hogs:
                all_hogs_set.add(hog)
                per_hog_node[hog] = node
                per_hog_og[hog] = ogname

    return all_hogs_set, per_hog_og, per_hog_node
