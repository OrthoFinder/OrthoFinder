import os
import argparse

from ete4 import Tree

def Species_Tree_output(target_dir):
    
    species_tree = "/".join([target_dir,"Species_Tree","SpeciesTree_rooted.txt"]) 
    Species_Tree_ETE = Tree(open(species_tree).read())
    return str(Species_Tree_ETE)


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
    print(species_tree)

    
    
    
