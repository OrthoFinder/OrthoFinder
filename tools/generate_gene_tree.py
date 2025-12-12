import os
import argparse

from ete4 import Tree


def get_random_gene_tree(target_dir):
    Resolved_genes_tree_path = "/".join([target_dir,"Resolved_Gene_Trees","Resolved_Gene_Trees.txt"])
    with open(Resolved_genes_tree_path) as gene_tree_file:
        for line in gene_tree_file:
            OG_Tree_Pair = line.rstrip().split(": ")
            ete3_tree = Tree(OG_Tree_Pair[1], parser=1)

            if len(list(ete3_tree.leaves())) < 12 and 4 < len(list(ete3_tree.leaves())):
                #new_ete3_tree = reformat_ete3_tree(ete3_tree)
                new_ete3_tree = ete3_tree
                orthogroup = OG_Tree_Pair[0]
                break
    return orthogroup, str(new_ete3_tree)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()
    target_dir = args.path

    if not os.path.exists(target_dir):
        print("The target directory doesn't exist")
        raise SystemExit(1)
    
    if not os.path.isdir(target_dir):
        print("The input path is not a directory")
        raise SystemExit(1)
    
    gene_tree_orthogroup, new_ete3_tree = get_random_gene_tree(target_dir)
    print(new_ete3_tree)

    
    
    
