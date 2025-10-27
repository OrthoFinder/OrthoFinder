import os 
import sys
import re
import traceback
from collections import defaultdict
import pytest
from orthofinder.run.__main__ import main



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

# def get_dir_path(arg):
#     directory = os.path.abspath(arg)
#     if not os.path.isfile(directory) and directory[-1] != os.sep:
#         directory += os.sep
#     if not os.path.exists(directory):
#         print("Specified directory doesn't exist: %s" % directory)
#         sys.exit(1)
#     return directory

# def get_file_path(arg):
#     file_path = os.path.abspath(arg)
#     directory = os.path.dirname(file_path)
#     if not os.path.exists(directory):
#         print("Directory points to the file doesn't exist: %s" % directory)
#         sys.exit(1)
#     if not os.path.isfile(file_path):
#         print("Specified file doesn't exist: %s" % file_path)
#         sys.exit(1)
#     return file_path
