import os
import argparse

from ete4 import Tree
from tabulate import tabulate


def gene_duplications(target_dir):
    gene_duplications = "/".join([target_dir,"Gene_Duplication_Events"])
    gene_duplication_tree = Tree(open("/".join([gene_duplications,"SpeciesTree_Gene_Duplications_0.5_Support.txt"])).read(),parser=1)

    ## Example NODE for now N2
    for node in gene_duplication_tree.traverse():
        if node.name.startswith("N2"):
            node_example = node.name.split("_")[0]
            number_of_duplications_for_node = node.name.split("_")[1]
    
    
    #print(gene_duplication_tree)
    first_duplication_lines = []
    gene_duplication_table = "/".join([gene_duplications,"Duplications.tsv"])
    with open(gene_duplication_table) as duplication_table:
        for lines in range(0,6):
            line = next(duplication_table)
            first_duplication_lines.append(line.rstrip())
    
    gap = 0
    for col_name in first_duplication_lines[0].split("\t"):
        if gap < len(col_name):
            gap = len(col_name)
    reformatted_table = []
    
    for line in first_duplication_lines:
        new_line = []
        for element in line.split("\t"):
            if len(element) < gap + 3:
                new_line.append(element)
            else:
                new_line.append(element[:gap + 3])
        reformatted_table.append(new_line)
    return tabulate(reformatted_table[1:], reformatted_table[0], tablefmt="github"),\
        gene_duplication_tree.to_str(props=['name'], compact=True),\
        node_example,\
        number_of_duplications_for_node
    




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
    
    duplication_table, duplication_tree, node_example, number_of_duplications_for_node = gene_duplications(target_dir)
    
    print(duplication_tree)