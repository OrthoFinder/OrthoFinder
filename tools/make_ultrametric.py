#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Modify branch lengths on rooted tree files so that they are ultrametric.

INPUT
  The input may be:
    - a single rooted Newick tree file
    - a directory containing rooted tree files

OUTPUT
  For a single input tree, the default output is:

      INPUT.ultrametric.tre

  For a directory input, ultrametric trees are written to:

      INPUT/ultrametric_trees/

Examples
--------
make_ultrametric species_tree.tre

make_ultrametric species_tree.tre --root-age 100

make_ultrametric species_tree.tre -o species_tree.ultrametric.tre

make_ultrametric Gene_Trees/

make_ultrametric Gene_Trees/ -o Ultrametric_Gene_Trees/
"""

import os
import sys
import argparse
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)
src_dir = os.path.join(repo_root, "src")

if os.path.isdir(os.path.join(src_dir, "orthofinder")) and src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from orthofinder.tools import tree
from orthofinder.utils import util


TREE_EXTENSIONS = {
    ".tre",
    ".tree",
    ".nwk",
    ".newick",
    ".txt",
}


def AveDist(node):
    return np.average([node.get_distance(l) for l in node.get_leaf_names()])


def Fail():
    print("ERROR: An error occurred, please review error messages for more information.")
    sys.exit(1)


def CheckTree(t):
    if len(t.get_children()) != 2:
        print("Input tree must be rooted")
        Fail()


def is_tree_file(path):
    return os.path.isfile(path) and os.path.splitext(path)[1].lower() in TREE_EXTENSIONS

def get_tree_files(input_path):
    if os.path.isfile(input_path):
        return [input_path]

    if os.path.isdir(input_path):
        tree_files = []

        for root, _, files in os.walk(input_path):
            for filename in sorted(files):
                path = os.path.join(root, filename)

                if is_tree_file(path):
                    tree_files.append(path)

        return sorted(tree_files)

    print("ERROR: input path does not exist: %s" % input_path)
    Fail()


def process_input(input_path, root_age=None, output=None):
    input_path = os.path.abspath(input_path)

    if os.path.isfile(input_path):
        make_ultrametric(
            tree_fn=input_path,
            root_age=root_age,
            output_file=output,
        )
        return

    tree_files = get_tree_files(input_path)

    if not tree_files:
        print("ERROR: no tree files found in directory: %s" % input_path)
        Fail()

    if output is None:
        output_dir = os.path.join(input_path, "ultrametric_trees")
    else:
        output_dir = output

    os.makedirs(output_dir, exist_ok=True)

    for tree_fn in tree_files:
        rel_path = os.path.relpath(tree_fn, input_path)
        rel_root, _ = os.path.splitext(rel_path)

        outfn = os.path.join(output_dir, rel_root + ".ultrametric.tre")
        os.makedirs(os.path.dirname(outfn), exist_ok=True)

        print("\nProcessing tree: %s" % tree_fn)
        make_ultrametric(
            tree_fn=tree_fn,
            root_age=root_age,
            output_file=outfn,
        )

def make_ultrametric(tree_fn, root_age=None, output_file=None):
    if not os.path.exists(tree_fn):
        print("Input tree file does not exist: %s" % tree_fn)
        Fail()

    t = tree.Tree(tree_fn, format=1)
    CheckTree(t)

    d = AveDist(t)
    print("Average distance from root to leaves: %f" % d)

    for n in t.traverse("preorder"):
        if n.is_root():
            n.dist = 0
            continue

        # Work downwards, setting the branch distances from the top down.
        x = t.get_distance(n) - n.dist
        y = n.dist

        print("\nTaxa:")
        print(", ".join(n.get_leaf_names()))

        if n.is_leaf():
            z = 0.0
        else:
            z = AveDist(n)

        print("Distance of parent node from root: %f" % x)
        print("Current branch length: %f" % y)
        print("Average distance to leaves: %f" % z)

        if (y + z) == 0.0:
            n.dist = 0
        else:
            f = (d - x) / (y + z)
            n.dist = f * n.dist

        print("Branch length for ultrametric tree: %f" % n.dist)

    if root_age is not None:
        x = root_age / d
        print(
            "\nRescaling branch lengths by factor of %0.2f so that root age is %f"
            % (x, root_age)
        )

        for n in t.traverse():
            if n.is_root():
                continue
            n.dist = x * n.dist

    if output_file is None:
        output_file = tree_fn + ".ultrametric.tre"

    t.write(outfile=output_file, format=5)
    print("\nUltrametric tree written to: %s\n" % output_file)


def main(args=None):
    with util.Finalise():
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

        parser.add_argument(
            "input",
            help="Input rooted Newick tree file or directory containing tree files",
        )
        parser.add_argument(
            "-r",
            "--root-age",
            dest="root_age",
            type=float,
            default=None,
            help="Rescale branch lengths so that the root age equals this value",
        )

        parser.add_argument(
            "-o",
            "--output",
            default=None,
            help=(
                "Output filename for single-tree input, or output directory for "
                "directory input [default: INPUT.ultrametric.tre or INPUT/ultrametric_trees]"
            ),
        )

        parsed = parser.parse_args(args)

        process_input(
            input_path=parsed.input,
            root_age=parsed.root_age,
            output=parsed.output,
        )

if __name__ == "__main__":
    main()