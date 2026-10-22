#!/usr/bin/env python3

"""
Write the concatenated alignment of the orthogroups chosen for a species tree,
from the results of a finished run (the same files as a run that infers the
species tree from it, or a run with -ST --concat-alignment):

    WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa
    Species_Tree/Orthogroups_for_concatenated_alignment.txt

The rows are named by species ID (0, 1, ...), as in the run, or with
--species-names by species name. INPUT is an OrthoFinder results directory;
it needs the alignments in its working directory (WorkingDirectory/Alignments_ids).

Examples
--------
species_tree_alignment Results_*/

species_tree_alignment Results_*/ -o species_tree_alignment/ --species-names
"""

import argparse
import os
import sys
import types

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)
src_dir = os.path.join(repo_root, "src")

if os.path.isdir(os.path.join(src_dir, "orthofinder")) and src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from orthofinder.mcl_orthogroups import mcl
from orthofinder.mcl_orthogroups.orthogroups_set import Seq
from orthofinder.run import previous_run
from orthofinder.species_tree import concatenated_alignment
from orthofinder.utils import gene_ids, run_state

ALIGNMENT_FN = "SpeciesTreeAlignment.fa"
ORTHOGROUPS_FN = "Orthogroups_for_concatenated_alignment.txt"


def run_files(results_dir):
    """
    (orthogroups file, the working directories to look in, settings) of a
    results directory, from its run_state.json or else its Log.txt (as a
    restart finds them).
    """
    found_dir, state = run_state.find(results_dir)
    if state:
        if state.get("run", {}).get("analysis") == "assign":
            raise ValueError("this is an --assign run: use the results of the core run (-f)")
        dirs = run_state.dirs(found_dir, state)
        working_dirs = [d for d in [dirs.get("trees"), dirs.get("working")] + dirs.get("base", []) if d is not None]
        orthogroups_fn = run_state.file(found_dir, state, "orthogroups")
        if orthogroups_fn is None:
            raise ValueError("the orthogroups file recorded in run_state.json is missing")
        return orthogroups_fn, working_dirs, state.get("settings", {})
    try:
        locator = previous_run.PreviousFilesLocatorNew(
            types.SimpleNamespace(q_start_from_fasta=False, q_start_from_blast=False), results_dir)
    except previous_run.Unprocessable:
        raise ValueError("no run_state.json or Log.txt to find the run's files in")
    if locator.q_assign:
        raise ValueError("this is an --assign run: use the results of the core run (-f)")
    wd_base, orthogroups_fn, wd_trees, _, _ = locator.get_start_from_trees()
    if orthogroups_fn is None or not os.path.exists(orthogroups_fn):
        raise ValueError("the orthogroups file in Log.txt is missing: %s" % orthogroups_fn)
    settings = previous_run.state_from_log(results_dir).get("settings", {})
    return orthogroups_fn, [d for d in [wd_trees] + list(wd_base) if d is not None], settings


def first_with(working_dirs, path):
    """The first working directory that has path in it, or None."""
    return next((d for d in working_dirs if os.path.exists(os.path.join(d, path))), None)


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="species_tree_alignment",
        description="Write the concatenated alignment of the orthogroups chosen for a species tree, "
                    "from the results of a finished run (as a run without -ST does).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  species_tree_alignment Results_*/\n"
               "  species_tree_alignment Results_*/ -o species_tree_alignment/ --species-names")
    parser.add_argument("results_dir", metavar="INPUT", help="an OrthoFinder results directory (Results_<name>)")
    parser.add_argument("-o", "--output-dir", metavar="DIR",
                        help="write both files to DIR (by default where the run writes them: "
                             "WorkingDirectory/Alignments_ids/%s and Species_Tree/%s)" % (ALIGNMENT_FN, ORTHOGROUPS_FN))
    parser.add_argument("--species-names", action="store_true",
                        help="name the rows by species name (by default by species ID, as in the run)")
    a = parser.parse_args(args)
    results_dir = os.path.join(os.path.normpath(a.results_dir), "")
    if not os.path.isdir(os.path.join(results_dir, "Orthogroups")):
        parser.error("not an OrthoFinder results directory (no Orthogroups directory): %s" % a.results_dir)
    try:
        orthogroups_fn, working_dirs, settings = run_files(results_dir)
    except ValueError as e:
        parser.error(str(e))

    # The orthogroups as the run has them (OrthoGroupsSet.get_ogs_all)
    min_seq = int(settings.get("min_seq", 4))
    ogs = [[Seq(g) for g in og] for og in mcl.get_predicted_ogs(orthogroups_fn)]
    if settings.get("tree_program") in ["raxml", "raxml-ng"]:
        ogs = [og for og in ogs if len(og) >= min_seq]
    occupancy, iogs4 = concatenated_alignment.orthogroup_occupancy(ogs, min_seq)
    if occupancy.n_species < 4:
        print("ERROR: a species tree is only inferred for four or more species", file=sys.stderr)
        return 1
    i_ogs, f_single_copy = concatenated_alignment.determine_orthogroups_for_species_tree(occupancy, iogs4)

    align_dir = first_with(working_dirs, "Alignments_ids")
    if align_dir is None:
        print("ERROR: no alignments (Alignments_ids) in the working directory: %s" % working_dirs[0], file=sys.stderr)
        return 1
    align_dir = os.path.join(align_dir, "Alignments_ids", "")

    def alignment_fn(i_og):
        return align_dir + "OG%07d.fa" % i_og

    missing = [i_og for i_og in i_ogs if not os.path.exists(alignment_fn(i_og))]
    if missing:
        print("ERROR: %d of the %d alignments are not in %s, e.g. %s" % (
            len(missing), len(i_ogs), align_dir, os.path.basename(alignment_fn(missing[0]))), file=sys.stderr)
        return 1

    species_names = None
    if a.species_names:
        ids_dir = first_with(working_dirs, "SpeciesIDs.txt")
        if ids_dir is None:
            print("ERROR: no SpeciesIDs.txt in the working directory: %s" % working_dirs[0], file=sys.stderr)
            return 1
        species_names = gene_ids.species_names(os.path.join(ids_dir, "SpeciesIDs.txt"))

    if a.output_dir:
        os.makedirs(a.output_dir, exist_ok=True)
        out_alignment = os.path.join(a.output_dir, ALIGNMENT_FN)
        out_orthogroups = os.path.join(a.output_dir, ORTHOGROUPS_FN)
    else:
        out_alignment = align_dir + ALIGNMENT_FN
        os.makedirs(os.path.join(results_dir, "Species_Tree"), exist_ok=True)
        out_orthogroups = os.path.join(results_dir, "Species_Tree", ORTHOGROUPS_FN)

    print("Species tree: Using %d orthogroups with minimum of %0.1f%% of species having single-copy genes in any orthogroup"
          % (len(i_ogs), 100. * f_single_copy))
    concatenated_alignment.create_concatenated_alignment(
        i_ogs, ogs, alignment_fn, out_alignment, f_single_copy,
        q_progress=True, species_names=species_names)
    with open(out_orthogroups, "w", encoding="utf-8", newline="\n") as outfile:
        for i_og in i_ogs:
            outfile.write("OG%07d\n" % i_og)
    print("Concatenated alignment: %s" % out_alignment)
    print("Orthogroups used: %s" % out_orthogroups)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
