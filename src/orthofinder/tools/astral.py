import os
import sys
import glob

from . import tree


def create_input_file(d, output_fn, n_skip=50):
    """
    The gene trees for ASTRAL-Pro, with species IDs as leaf names. The trees
    of the first n_skip orthogroups (the largest, -nk) are left out, unless
    that would leave none (an analysis with few orthogroups): then all are
    used. Raises RuntimeError if there are no gene trees (ASTRAL-Pro would
    crash on an empty input).
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
    with open(output_fn, 'w') as outfile:
        for fn in kept:
            t = tree.Tree(fn,format=1)
            for n in t:
                n.name = n.name.split("_")[0]
            outfile.write(t.write(format=1) + "\n")

def get_astral_command(astral_input, species_tree, threads):
    return " ".join(["astral-pro", "-i", astral_input, "-o", species_tree, "-t", str(threads)])


# if __name__ == "__main__":
#     d = sys.argv[1]
#     output_fn = "astral_input.nwk"
#     create_input_file(d, output_fn)
