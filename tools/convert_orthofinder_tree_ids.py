#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert OrthoFinder tree sequence IDs to gene accessions.

INPUT
  TreeInput may be:
    - a single tree file
    - a directory containing tree files

REQUIRED FILES
  SequenceIDs is normally found in:

      Results_<Date>/WorkingDirectory/SequenceIDs.txt

  SpeciesIDs is optional and normally found in:

      Results_<Date>/WorkingDirectory/SpeciesIDs.txt

OUTPUT
  For each input tree, a new tree is written with "_accessions" added
  before the original file extension.

Examples
--------
convert_orthofinder_tree_ids Gene_Trees/OG0000000_tree.txt SequenceIDs.txt

convert_orthofinder_tree_ids Gene_Trees/ SequenceIDs.txt SpeciesIDs.txt
"""

import os
import sys
import glob
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)
src_dir = os.path.join(repo_root, "src")

if os.path.isdir(os.path.join(src_dir, "orthofinder")) and src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from orthofinder.tools import tree
from orthofinder.utils import util



def invert_mapping(ids_map):
    acc_to_ids = {}
    duplicated = set()

    for seq_id, accession in ids_map.items():
        if accession in acc_to_ids:
            duplicated.add(accession)
        else:
            acc_to_ids[accession] = seq_id

    for accession in duplicated:
        del acc_to_ids[accession]

    return acc_to_ids


def convert_resolved_gene_trees_file(idsMap, treeFilename, outputFilename):
    with open(treeFilename, "r") as infile, open(outputFilename, "w") as outfile:
        for line in infile:
            line = line.strip()

            if not line:
                continue

            if ":" not in line:
                outfile.write(line + "\n")
                continue

            og_name, newick = line.split(":", 1)
            newick = newick.strip()

            tmp_tree = tree.Tree(newick, format=1)

            for node in tmp_tree.get_leaves():
                if node.name in idsMap:
                    node.name = idsMap[node.name]

            outfile.write("%s: %s\n" % (og_name, tmp_tree.write(format=5).strip()))

def ReplaceFileWithNewIDs(idsMap, treeFilename, newTreeFilename):
    with open(treeFilename, "r") as infile:
        first_line = infile.readline().strip()

    if first_line.startswith("OG") and ":" in first_line:
        convert_resolved_gene_trees_file(idsMap, treeFilename, newTreeFilename)
        return

    qHaveSupport = False
    qHaveInternalNames = False

    try:
        t = tree.Tree(treeFilename, format=2)
        qHaveSupport = True
    except Exception:
        try:
            t = tree.Tree(treeFilename, format=3)
            qHaveInternalNames = True
        except Exception:
            t = tree.Tree(treeFilename)

    for node in t.get_leaves():
        if node.name in idsMap:
            node.name = idsMap[node.name]

    if qHaveSupport:
        t.write(outfile=newTreeFilename)
    elif qHaveInternalNames:
        t.write(outfile=newTreeFilename, format=3)
    else:
        t.write(outfile=newTreeFilename, format=5)


def GetSpeciesSequenceIDsDict(sequenceIDsFilename, speciesIDsFN=None):
    try:
        extract = util.FirstWordExtractor(sequenceIDsFilename)
    except RuntimeError as error:
        msg = str(error)
        print(msg)

        if msg.startswith("ERROR"):
            util.Fail()

        print(
            "Tried to use only the first part of the accession in order to list "
            "the sequences in each orthogroup more concisely, but these were not "
            "unique. The full accession line will be used instead.\n"
        )

        extract = util.FullAccession(sequenceIDsFilename)

    idsDict = extract.GetIDToNameDict()

    if speciesIDsFN is not None:
        speciesDict = util.FullAccession(speciesIDsFN).GetIDToNameDict()
        speciesDict = {
            k: v.rsplit(".", 1)[0].replace(".", "_").replace(" ", "_")
            for k, v in speciesDict.items()
        }
        idsDict = {
            seqID: speciesDict[seqID.split("_")[0]] + "_" + name
            for seqID, name in idsDict.items()
        }

    return idsDict


def get_tree_files(tree_input):
    if os.path.isfile(tree_input):
        return [tree_input]

    if os.path.isdir(tree_input):
        return sorted(
            f for f in glob.glob(os.path.join(tree_input, "*"))
            if os.path.isfile(f)
        )

    print("ERROR: tree input does not exist: %s" % tree_input)
    util.Fail()


def convert_tree_ids(tree_input, sequence_ids, species_ids=None, direction="id-to-accession"):
    idsDict = GetSpeciesSequenceIDsDict(sequence_ids, species_ids)

    if direction == "accession-to-id":
        idsDict = invert_mapping(idsDict)

    filesToDo = get_tree_files(tree_input)

    for treeFilename in filesToDo:
        pathfilename, ext = os.path.splitext(treeFilename)

        if direction == "id-to-accession":
            suffix = "_accessions"
        else:
            suffix = "_ids"

        newFilename = pathfilename + suffix + ext

        sys.stdout.write(newFilename)

        try:
            ReplaceFileWithNewIDs(idsDict, treeFilename, newFilename)
        except Exception as error:
            sys.stdout.write(" - skipped: %s" % error)

        print("")


def main(args=None):
    with util.Finalise():
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

        parser.add_argument(
            "tree_input",
            help="Input tree file or directory containing tree files",
        )

        parser.add_argument(
            "sequence_ids",
            help="SequenceIDs.txt file from the OrthoFinder WorkingDirectory",
        )

        parser.add_argument(
            "species_ids",
            nargs="?",
            default=None,
            help="Optional SpeciesIDs.txt file from the OrthoFinder WorkingDirectory",
        )

        parser.add_argument(
            "-d",
            "--direction",
            choices=["id-to-accession", "accession-to-id"],
            default="id-to-accession",
            help="Conversion direction [default: id-to-accession]",
        )

        parsed = parser.parse_args(args)

        convert_tree_ids(
            tree_input=parsed.tree_input,
            sequence_ids=parsed.sequence_ids,
            species_ids=parsed.species_ids,
            direction=parsed.direction,
        )


if __name__ == "__main__":
    main()