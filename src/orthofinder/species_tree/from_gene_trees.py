"""
The species tree from the gene trees, by a species tree program (-ST, e.g.
ASTRAL-Pro), configured in config.json with program_type "species_tree".
"""
import os
import glob
import re

from ..phylo import tree
from ..utils import parallel_task_manager


def create_input_file(d, output_fn, n_skip=50):
    """
    The gene trees for the species tree program, one per line, with species
    IDs as leaf names. The trees of the first n_skip orthogroups (the largest,
    -nk) are left out, unless that would leave none (an analysis with few
    orthogroups): then all are used. Raises RuntimeError if there are no gene
    trees (ASTRAL-Pro would crash on an empty input).
    """
    # in orthogroup order (glob's order varies), so the input is the same every run
    trees = sorted((int(os.path.basename(fn).split(".")[0][2:]), fn) for fn in glob.glob(d + "/*txt"))
    if not trees:
        raise RuntimeError("No gene trees to infer the species tree from (in %s)" % d)
    kept = [fn for iog, fn in trees if iog >= n_skip]
    if not kept:
        print("Only %d gene tree(s), all among the first %d orthogroups (-nk): all are used for the "
              "species tree" % (len(trees), n_skip))
        kept = [fn for _, fn in trees]
    with open(output_fn, 'w', encoding="utf-8", newline="\n") as outfile:
        for fn in kept:
            t = tree.Tree(fn,format=1)
            for n in t:
                n.name = n.name.split("_")[0]
            outfile.write(t.write(format=1) + "\n")


def infer_species_tree(prog_caller, program, trees_dir, input_fn, species_tree_fn, threads, n_skip=50):
    """
    Infer the species tree (unrooted, species IDs) from the gene trees in
    trees_dir with a species tree program, writing its input to input_fn.
    """
    create_input_file(trees_dir, input_fn, n_skip=n_skip)
    cmd, actual_target = prog_caller.get_species_tree_command(program, input_fn, species_tree_fn, threads)
    # Such programs can write their progress log to stderr (e.g. ASTRAL-Pro):
    # only shown if it fails.
    parallel_task_manager.run_command(cmd, q_print_stderr=False, raise_on_error=True)
    if actual_target is not None:
        actual, target = actual_target
        if os.path.exists(actual):
            os.rename(actual, target)
    if os.path.exists(species_tree_fn):
        remove_internal_labels(species_tree_fn)


def remove_internal_labels(fn):
    """
    The species tree without the labels of its internal nodes, which programs
    write e.g. as support values (ASTER >= 1.20: "(1:0.4,2:0.5)1:0.48"): read
    as node names, one can be the same as a species ID.
    """
    with open(fn, encoding="utf-8") as infile:
        text = infile.read()
    with open(fn, "w", encoding="utf-8", newline="\n") as outfile:
        outfile.write(re.sub(r"\)[^:,;()\[\]]+", ")", text))
