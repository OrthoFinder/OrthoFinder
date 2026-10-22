#!/usr/bin/env python3

"""
Write the sequences of the MCL orthogroups (before the HOGs) of an OrthoFinder
run, a FASTA file for each, from its working directory:

    MCL_Orthogroup_Sequences/OG0000000.fa, ...

The orthogroups are those of Orthogroups/MCL_Orthogroups.tsv (written e.g. by
-og), numbered as there, with the gene names as in the input FASTA files.
INPUT is an OrthoFinder results directory; it needs its working directory
(WorkingDirectory: the species' FASTA files and the MCL clusters). The
sequences of the final orthogroups (from the HOGs) are in
Orthogroup_Sequences/ after a full run.

Examples
--------
mcl_orthogroup_sequences Results_*/

mcl_orthogroup_sequences Results_*/ -o mcl_sequences/
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

from orthofinder.gene_tree.msa_trees import FastaWriter
from orthofinder.mcl_orthogroups import mcl
from orthofinder.run import previous_run
from orthofinder.utils import gene_ids, run_state

OUT_DIR = "MCL_Orthogroup_Sequences"


def run_files(results_dir):
    """
    (MCL clusters file, the working directories to look in) of a results
    directory, from its run_state.json or else its Log.txt (as a restart finds them).
    """
    found_dir, state = run_state.find(results_dir)
    if state:
        if state.get("run", {}).get("analysis") == "assign":
            raise ValueError("this is an --assign run: use the results of the core run (-f)")
        dirs = run_state.dirs(found_dir, state)
        working_dirs = [d for d in [dirs.get("working")] + dirs.get("base", []) if d is not None]
        clusters_fn = run_state.file(found_dir, state, "orthogroups")
        if clusters_fn is None:
            raise ValueError("the orthogroups file recorded in run_state.json is missing")
        return clusters_fn, working_dirs
    try:
        locator = previous_run.PreviousFilesLocatorNew(
            types.SimpleNamespace(q_start_from_fasta=False, q_start_from_blast=False), results_dir)
    except previous_run.Unprocessable:
        raise ValueError("no run_state.json or Log.txt to find the run's files in")
    if locator.q_assign:
        raise ValueError("this is an --assign run: use the results of the core run (-f)")
    wd_base, clusters_fn = locator.get_start_from_ogs()
    if clusters_fn is None or not os.path.exists(clusters_fn):
        raise ValueError("the orthogroups file in Log.txt is missing: %s" % clusters_fn)
    return clusters_fn, list(wd_base)


def first_with(working_dirs, path):
    """The first working directory that has path in it, or None."""
    return next((d for d in working_dirs if os.path.exists(os.path.join(d, path))), None)


def gene_names(sequence_ids_fn):
    """{sequence ID: gene name}, as the run names the genes in its results."""
    try:
        return mcl.id_full_dict([sequence_ids_fn], func=gene_ids.FirstWordExtractor)
    except RuntimeError:   # the first words are not unique: the full accessions
        return mcl.id_full_dict([sequence_ids_fn], func=gene_ids.FullAccession)


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="mcl_orthogroup_sequences",
        description="Write the sequences of the MCL orthogroups (before the HOGs) of an OrthoFinder run, "
                    "a FASTA file for each, from its working directory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  mcl_orthogroup_sequences Results_*/\n"
               "  mcl_orthogroup_sequences Results_*/ -o mcl_sequences/")
    parser.add_argument("results_dir", metavar="INPUT", help="an OrthoFinder results directory (Results_<name>)")
    parser.add_argument("-o", "--output-dir", metavar="DIR",
                        help="write the FASTA files to DIR (by default %s/ in the results directory)" % OUT_DIR)
    a = parser.parse_args(args)
    results_dir = os.path.join(os.path.normpath(a.results_dir), "")
    if not os.path.isdir(results_dir):
        parser.error("no such directory: %s" % a.results_dir)
    try:
        clusters_fn, working_dirs = run_files(results_dir)
    except ValueError as e:
        parser.error(str(e))

    ids_dir = first_with(working_dirs, "SequenceIDs.txt")
    if ids_dir is None or first_with(working_dirs, "SpeciesIDs.txt") is None:
        print("ERROR: no SequenceIDs.txt and SpeciesIDs.txt in the working directory: %s" % working_dirs[0],
              file=sys.stderr)
        return 1
    species_to_use = sorted(previous_run.species_name_dict(
        os.path.join(first_with(working_dirs, "SpeciesIDs.txt"), "SpeciesIDs.txt")))
    fasta_dirs = [os.path.join(d, "") for d in working_dirs]
    missing = [i for i in species_to_use if first_with(fasta_dirs, "Species%d.fa" % i) is None]
    if missing:
        print("ERROR: the FASTA files of species %s are not in the working directory: %s"
              % (", ".join(map(str, missing)), working_dirs[0]), file=sys.stderr)
        return 1

    ogs = mcl.get_predicted_ogs(clusters_fn)
    names = gene_names(os.path.join(ids_dir, "SequenceIDs.txt"))
    fasta_writer = FastaWriter(fasta_dirs, species_to_use)
    out_dir = a.output_dir or os.path.join(results_dir, OUT_DIR)
    os.makedirs(out_dir, exist_ok=True)
    for i_og, og in enumerate(ogs):
        fasta_writer.write_seqs_to_fasta_with_new_accessions(
            list(og), os.path.join(out_dir, "OG%07d.fa" % i_og), names)
    print("Wrote the sequences of %d MCL orthogroups to %s" % (len(ogs), out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
