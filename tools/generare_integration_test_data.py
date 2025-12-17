# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 13:16:52 2025

@author: Biol0216

1)
CREATE 1000x random amino acid sequences
2) 
Given a tree run alisim to generate sequnce - reformat
3)
Move Alisim alignment into proteome file (remove gaps etc...)



"""



def Generate_Random_Amino_Acids(Output_File_Name, Seq_to_generate, gene_length_range):
    # Output_File_Name = str : output file name
    # Seq_to_generate = int : number of sequences to generate
    # gene_length_range = list : upper and low bound for sequence length
    AA = "ARNDBCEQZGHILKMFPSTWYV"
    ## gene name convention = "Gene_" + int
    if Output_File_Name in os.listdir():
        os.remove(Output_File_Name)
    
    
    for i in range(0,Seq_to_generate):
        length = random.randint(gene_length_range[0], gene_length_range[1])
        seq = ""
        for I in range(0,length):
            seq += random.choice(AA)
        
        with open(Output_File_Name,"a") as output_file:
            output_file.write(">Gene_" + str(i) + "\n" + seq + "\n")
        
    
def Generate_Proteomes(Tree_Folder,Number_of_Sequences,Root_file,Proteome_Folder):
    # Folder containing newick files with trees to test = tree leaves should be in the format = TEST_PROTEOME_GENECODE
    # Number_of_Sequences = number of sequences in root fasta
    # Root_file = fasta path containing the sequences
    # Proteome Folder = proteome fold path
    ## Protomes will the in the Format A.fasta = >TEST_A_## 
    ## The second element after "_" will the the proteome file file...
    
    if Proteome_Folder in os.listdir():
        for file in os.listdir(Proteome_Folder):
          os.remove(Proteome_Folder + "/" + file)          
    
    starting_genes  = list(range(0, Number_of_Sequences))
    for tree_file in os.listdir(Tree_Folder):
        starting_gene = "Gene_" + str(random.choice(starting_genes))
        starting_genes.remove(int(starting_gene.replace("Gene_","")))
        
        ############# Make tree ultrametric and write to temp file.
        new_tree_length = 0.2
        from ete3 import Tree
        print(tree_file)
        Input_tree = Tree(Tree_Folder + "/" + tree_file)
        Input_tree.convert_to_ultrametric(tree_length=new_tree_length, strategy='balanced')
        new_tree_path = (Tree_Folder + "/" + tree_file).replace(".nwk","_ultrametric.nwk")
        Input_tree.write(format=1, outfile=new_tree_path)

        command = "iqtree --alisim %s -t %s --root-seq %s,%s -af fasta -m JTT+G4" % ("Temp_File",new_tree_path,Root_file,starting_gene)
        ## should replace with subprocess for release - to avoid output spam...
        os.system(command)
 
        ###### sequence file   : 
        sequence_fasta = open("Temp_File.fa").read().split(">")[1:]
        for seq in sequence_fasta:
            proteome = seq.split("_")[1]
            with open(Proteome_Folder + "/" + proteome + ".fa","a") as proteome_file:
                proteome_file.write(">" + seq)
        
        os.remove(new_tree_path + ".log")
        os.remove("Temp_File.fa")        
        os.remove(new_tree_path)
        ##### 
        
        
        
        
        
    #     command = "diamond blastp -d %s -q %s -o %s --outfmt 5 --quiet --evalue 0.05 -k 0 --max-hsps 0 --threads 1" % (Input,Sequence,blast_results)

    #iqtree --alisim TEST -t Gene_Tree.nwk --root-seq STARTING_SEQUENCE.fa,START_SEQUENCE -af fasta -m JTT+G4
    
def Generate_random_orthogroups(Species_Tree,Trees_to_generate,number_of_single_copies,max_size):
    
    # Species tree to create single copies from
    # trees_to_generate = directory to add the trees too
    # number_of_single_copies : how many single copy genes to make
    # Makes many single copy genes trees for OrthoFinder to increase gene counts
    from ete3 import Tree
    input_tree = open(Species_Tree).read()
    Species = list(input_tree.replace(")","").replace("(","").replace(",","").replace(";",""))
    pos = 0
    for i in range(0,number_of_single_copies):
        new_tree = input_tree
        t = Tree(new_tree)
        ## must contain more than 4 species to make a tree.....
        new_number_of_leaves = random.choice(list(range(4,max_size)))   
        rand_tree  = Tree()
        rand_tree.populate(new_number_of_leaves)        
        for node in rand_tree:
            if node.is_leaf():
                pos = pos + 1
                Species_node = random.choice(Species)
                new_leaf_name = "OG_" + Species_node + "_" + str(pos)    
                node.name = new_leaf_name
        #t.write(format=1, outfile="new_tree.nw")
        path =  Trees_to_generate + "/" + str(i) + "_og.nwk"
        rand_tree.write(format=1, outfile=path)

    

if __name__ == "__main__":           
    import sys
    import random
    import os
    Trees_to_generate = "Tree_Tests"
    Species_Tree = "species_tree_assign.txt"
    #make_tree(Species_Tree)
    #### Generates a file oontaining random sequences strings
    Generate_Random_Amino_Acids("Random_Sequence.fasta",1000,[75,500])  
    
    #### generates many random orthogroups 
    #Generate_random_orthogroups(Species_Tree,Trees_to_generate,2,20)
    
    ##### creates proteomes from a set of input trees
    Generate_Proteomes(Trees_to_generate,1000,"Random_Sequence.fasta","Proteomes")
    
    











