import os
import argparse

from ete4 import Tree
from ete4.treeview import TreeStyle

def save_tree_image(tree, output_path):
    # def layout(node):
    #     if node.is_leaf == False:
    #         node.add_face(TextFace(node.name), 0, position='branch-right')

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False   
    ts.show_branch_support = False   
    ts.show_scale = False   
    # ts.layout_fn = layout
    ts.scale = 120    

    tree.render(output_path, w=800, units="px", tree_style=ts)

def Species_Tree_output(target_dir):
    
    species_tree = "/".join([target_dir,"Species_Tree","SpeciesTree_rooted.txt"]) 
    Species_Tree_ETE = Tree(open(species_tree).read())
    return Species_Tree_ETE


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
    

    species_tree = Species_Tree_output(target_dir)

    save_tree_image(species_tree, r"/tmp/species_tree.png")

    
    
    
