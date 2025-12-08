from directory_tree import DisplayTree
import os
import argparse


def get_path(Path):
	
    customPath: str = Path
    stringRepresentation: str = DisplayTree(customPath, stringRep=True, showHidden=True)
    altered_files = ""
    Last_Blast = ""
    Last_Species = ""

    HOGs = 0

    for file in stringRepresentation.split("\n"):
        if "OG0" not in file and "Blast" not in file and "Species" not in file and "__v__" not in file:
            altered_files += (file) + "\n"


        if "\u2514" in file and "OG" in file:
            new_files_start = "".join(file.split("OG")[:-1]).replace("\u2514","\u251C") + "OG" + "0"*7 + file.split("OG")[-1][7:]
            new_files_mid = "".join(file.split("OG")[:-1]).replace("\u2514","\u251C") + "OG" + "0"*3 + "####" + file.split("OG")[-1][7:]
            altered_files += (new_files_start) + "\n"
            altered_files += (new_files_mid) + "\n"
            altered_files += (file) + "\n"

        if "Blast" in file:
            if "Blast0_0.txt" in file:
                altered_files += (file) + "\n"
            elif "Blast0_1.txt" in file:
                altered_files += (file.replace("0_1","#_#")) + "\n"

        if "Species" in file:
            if "Species0.fa" in file:
                altered_files += (file) + "\n"
            if "Species1.fa" in file:
                altered_files += (file).replace("1","#") + "\n"

    Orthologs = (altered_files.split("\n").index("\u251C" + "\u2500" + "\u2500" + " Orthologues/"))
    HOGS = altered_files.split("\n").index("\u251C" + "\u2500" + "\u2500" + " Phylogenetic_Hierarchical_Orthogroups/")
    Misplaced = altered_files.split("\n").index("\u251C" + "\u2500" + "\u2500" + " Phylogenetically_Misplaced_Genes/")
    Xenologs = altered_files.split("\n").index("\u251C" + "\u2500" + "\u2500" + " Putative_Xenologs/")
    Gene_Trees = altered_files.split("\n").index("\u251C" + "\u2500" + "\u2500" + " Resolved_Gene_Trees/")
    altered_files_list = altered_files.split("\n")

    stripped_altered_files = altered_files_list[:Orthologs] 
    stripped_altered_files = stripped_altered_files + altered_files_list[Orthologs:Orthologs+1] 
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Orthologs+2].split(" ")[:-1]) + " Species_1.tsv"]
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Orthologs+2].split(" ")[:-1]) + " Species_#.tsv"]
    stripped_altered_files = stripped_altered_files + altered_files_list[HOGS:HOGS+1] 
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[HOGS+2].split(" ")[:-1]) + " N0.tsv"]
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[HOGS+2].split(" ")[:-1]) + " N#.tsv"]
    stripped_altered_files = stripped_altered_files + altered_files_list[Misplaced :Misplaced +1] 
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Misplaced +2].split(" ")[:-1]) + " Species_1.tsv"]
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Misplaced +2].split(" ")[:-1]) + " Species_#.tsv"]
    stripped_altered_files = stripped_altered_files + altered_files_list[Xenologs :Xenologs +1] 
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Xenologs +2].split(" ")[:-1]) + " Species_1.tsv"]
    stripped_altered_files = stripped_altered_files + [" ".join(altered_files_list[Xenologs +2].split(" ")[:-1]) + " Species_#.tsv"]
    stripped_altered_files = stripped_altered_files + altered_files_list[Gene_Trees:] 

    stripped_altered_files_str = "\n".join(stripped_altered_files)

    return stripped_altered_files_str


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()
    target_dir = args.path
    print(target_dir)
    if not os.path.exists(target_dir):
        print("The target directory doesn't exist")
        raise SystemExit(1)
    
    if not os.path.isdir(target_dir):
        print("The input path is not a directory")
        raise SystemExit(1)
    
    print(get_path(target_dir))

