"""Renaming the taxa of trees."""

from . import tree


def have_support_values(species_tree_fn_ids):
    q_have_support = False
    try:
        tree.Tree(species_tree_fn_ids, format=2)
        q_have_support = True
    except:
        pass
    return q_have_support


def rename_tree_taxa(
    tree_fn_or_tree,
    new_tree_filename,
    ids_map,
    q_support,
    q_fix_negatives=False,
    in_format=None,
    label=None,
    q_via_copy=False,
):
    """
    qViaCopy - create a copy of the tree and edit this copy. I.e. don't make changes to the original
    """
    if label != None:
        q_support = False
    q_have_support = False
    try:
        if type(tree_fn_or_tree) is tree.TreeNode:
            if q_via_copy:
                t = tree_fn_or_tree.copy("newick")
            else:
                t = tree_fn_or_tree
        else:
            q_have_support = False
            if in_format == None:
                try:
                    t = tree.Tree(tree_fn_or_tree, format=2)
                    q_have_support = True
                except:
                    t = tree.Tree(tree_fn_or_tree)
            else:
                t = tree.Tree(tree_fn_or_tree, format=in_format)
        for node in t.get_leaves():
            node.name = ids_map[node.name]
        if q_fix_negatives:
            tree_length = sum([n.dist for n in t.traverse() if n != t])
            sliver = tree_length * 1e-6
        i_node = 1
        for n in t.traverse():
            if q_fix_negatives and n.dist < 0.0:
                n.dist = sliver
            if label != None:
                if (not n.is_leaf()) and (not n.is_root()):
                    n.name = label + ("%d" % i_node)
                    i_node += 1
        if label != None:
            with open(new_tree_filename, "w") as outfile:
                outfile.write(
                    t.write(format=3)[:-1] + label + "0;"
                )  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        elif t.name == "N0" or t.name == "n0":
            with open(new_tree_filename, "w") as outfile:
                outfile.write(
                    t.write(format=3)[:-1] + t.name + ";"
                )  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        else:
            if q_support or q_have_support:
                t.write(outfile=new_tree_filename, format=2)
            else:
                t.write(outfile=new_tree_filename, format=5)
    except:
        pass
