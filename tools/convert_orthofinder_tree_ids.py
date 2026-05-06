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

try:
    from orthofinder.tools import tree
    from orthofinder.utils import util
except ImportError:
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_dir = os.path.join(repo_root, "src")
    sys.path.insert(0, src_dir)

    from orthofinder.tools import tree
    from orthofinder.utils import util


def ReplaceFileWithNewIDs(idsMap, treeFilename, newTreeFilename):
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


def convert_tree_ids(tree_input, sequence_ids, species_ids=None):
    idsDict = GetSpeciesSequenceIDsDict(sequence_ids, species_ids)
    filesToDo = get_tree_files(tree_input)

    for treeFilename in filesToDo:
        pathfilename, ext = os.path.splitext(treeFilename)
        newFilename = pathfilename + "_accessions" + ext

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

        parsed = parser.parse_args(args)

        convert_tree_ids(
            tree_input=parsed.tree_input,
            sequence_ids=parsed.sequence_ids,
            species_ids=parsed.species_ids,
        )


if __name__ == "__main__":
    main()