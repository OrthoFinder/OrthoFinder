"""
Draw the tree figures of the tutorial PDF from an OrthoFinder results directory.

    QT_QPA_PLATFORM=offscreen python draw_trees.py RESULTS_DIR OG OUT_DIR [ASSIGN_RESULTS_DIR]
"""
import os
import sys

from ete4 import Tree
from ete4.treeview import TextFace, TreeStyle


def style(scale=120):
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_scale = False
    ts.scale = scale
    ts.margin_left = ts.margin_right = ts.margin_top = ts.margin_bottom = 10
    return ts


def label_internal(ts, colour="#1f5fa8"):
    def layout(node):
        if not node.is_leaf and node.name:
            node.add_face(TextFace(" " + node.name, fsize=11, fgcolor=colour), 0, position="branch-top")
    ts.layout_fn = layout
    return ts


def read_tree(fn):
    with open(fn, encoding="utf-8") as infile:
        return Tree(infile.read().strip(), parser=1)


def gene_tree(results_dir, og):
    with open(os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt"), encoding="utf-8") as infile:
        for line in infile:
            name, newick = line.split(": ", 1)
            if name == og:
                return Tree(newick.strip(), parser=1)
    raise SystemExit("%s not found" % og)


def main():
    results_dir, og, out_dir = sys.argv[1:4]
    assign_dir = sys.argv[4] if len(sys.argv) > 4 else None
    os.makedirs(out_dir, exist_ok=True)
    read_tree(os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")).render(
        os.path.join(out_dir, "species_tree.png"), w=900, units="px", tree_style=label_internal(style()))
    gene_tree(results_dir, og).render(
        os.path.join(out_dir, "gene_tree.png"), w=1100, units="px", tree_style=label_internal(style(60), "#b0451f"))
    read_tree(os.path.join(results_dir, "Gene_Duplication_Events", "SpeciesTree_Gene_Duplications_0.5_Support.txt")).render(
        os.path.join(out_dir, "duplication_tree.png"), w=900, units="px", tree_style=label_internal(style()))
    if assign_dir:
        read_tree(os.path.join(assign_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")).render(
            os.path.join(out_dir, "assign_species_tree.png"), w=900, units="px", tree_style=label_internal(style()))


if __name__ == "__main__":
    main()
