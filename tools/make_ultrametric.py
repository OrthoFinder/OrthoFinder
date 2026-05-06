#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Modify branch lengths on a rooted tree so that it is ultrametric.

INPUT
  The input must be a rooted Newick tree.

OUTPUT
  By default, the output file is written as:

      INPUT.ultrametric.tre

Examples
--------
make_ultrametric species_tree.tre

make_ultrametric species_tree.tre --root-age 100

make_ultrametric species_tree.tre -o species_tree.ultrametric.tre
"""

import os
import sys
import argparse
import numpy as np

try:
    from orthofinder.tools import tree
    from orthofinder.utils import util
except ImportError:
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_dir = os.path.join(repo_root, "src")
    sys.path.insert(0, src_dir)

    from orthofinder.tools import tree
    from orthofinder.utils import util


def AveDist(node):
    return np.average([node.get_distance(l) for l in node.get_leaf_names()])


def Fail():
    print("ERROR: An error occurred, please review error messages for more information.")
    sys.exit(1)


def CheckTree(t):
    if len(t.get_children()) != 2:
        print("Input tree must be rooted")
        Fail()


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
            "tree_fn",
            help="Input rooted Newick tree file",
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
            help="Output tree filename [default: INPUT.ultrametric.tre]",
        )

        parsed = parser.parse_args(args)

        make_ultrametric(
            tree_fn=parsed.tree_fn,
            root_age=parsed.root_age,
            output_file=parsed.output,
        )


if __name__ == "__main__":
    main()