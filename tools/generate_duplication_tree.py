import os
import argparse

from ete4 import Tree
from ete4 import Tree
from ete4.treeview import TreeStyle, TextFace

def save_tree_image(tree, output_path):
    def layout(node):
        if node.is_leaf == False:
            node.add_face(TextFace(node.name), 0, position='branch-right')

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False   
    ts.show_branch_support = False   
    ts.show_scale = False   
    ts.layout_fn = layout
    ts.scale = 120    

    tree.render(output_path, w=800, units="px", tree_style=ts)

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
    return  gene_duplication_tree

    
        # gene_duplication_tree.to_str(props=['name'], compact=True),\
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
    
    duplication_tree = gene_duplications(target_dir)
    
    save_tree_image(duplication_tree, r"/tmp/duplication_tree.png")