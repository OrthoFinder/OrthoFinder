import os
import traceback
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import ete3


def index_files(id_dir, extension=".fa"):
    file_index = {}
    if id_dir is None:
        return file_index
    for entry in os.scandir(id_dir):
        if entry.is_file() and entry.name.endswith(extension):
            key = entry.name[:-len(extension)]
            file_index[key] = entry.path
    return file_index

# def read_fasta(file_path):
#     genes_dict = {}
#     with open(file_path, 'r') as fastaFile:
#         accession = None
#         sequence = []
#         for line in fastaFile:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if accession is not None:
#                     genes_dict[accession] = "".join(sequence)
#                 accession = line[1:]
#                 sequence = []
#             else:
#                 sequence.append(line)
#         if accession is not None:
#             genes_dict[accession] = "".join(sequence)
#     return genes_dict


def read_fasta(file_path):
    genes_dict = {}
    qFirst = True
    accession = ""
    sequence = ""
    with open(file_path, 'r') as fastaFile:
        for line in fastaFile:
            if line[0] == ">":
                # deal with old sequence
                if not qFirst:
                    genes_dict[accession] = sequence
                    sequence = ""
                qFirst = False
                # get id for new sequence
                accession = line[1:].rstrip()
            else:
                sequence += line
        genes_dict[accession] = sequence
    return genes_dict

def process_unique_og(
        unique_og, 
        hog_entries, 
        species_names, 
        name_dict,
        tree_file_index, 
        fasta_file_index,
    ):

    results = []
    try:
        gene_tree = None
        if unique_og in tree_file_index:
            with open(tree_file_index[unique_og], "r") as file:
                tree_data = file.read().strip()
                gene_tree = ete3.Tree(tree_data, quoted_node_names=True, format=1)
        else:
            print(f"WARNING: Tree file not found for {unique_og}")
            return results 
        
        gene_dict = None
        if unique_og in fasta_file_index:
            gene_dict = read_fasta(fasta_file_index[unique_og])
        else:
            print(f"WARNING: FASTA file not found for {unique_og}")
        
        for row in hog_entries:
            hog_name = name_dict.get(row["HOG"], row["HOG"])
            parent_node = row["Gene Tree Parent Clade"]

            if parent_node == "n0":
                subtree = gene_tree.copy()
            else:
                subtree_nodes = gene_tree.get_tree_root().search_nodes(name=parent_node)
                if not subtree_nodes:
                    print(f"WARNING: Parent node '{parent_node}' not found for {unique_og}")
                    continue
                subtree = subtree_nodes[0]
            
            current_leaves = [leaf.name for leaf in subtree.get_leaves()]
            
            expected_leaves = []
            for col in species_names:
                if row.get(col):
                    expected_leaves.extend([x.strip() for x in row[col].split(',')])
            
            if set(current_leaves) != set(expected_leaves):
                try:
                    subtree.prune(expected_leaves)
                except Exception as e:
                    print(f"WARNING: Could not prune subtree for {hog_name}: {e}")
            
            pruned_alignments = None
            if gene_dict is not None:
                pruned_alignments = {
                    gene: gene_dict[gene]
                    for gene in expected_leaves
                    if gene in gene_dict
                }
            results.append((hog_name, subtree, pruned_alignments))
    except Exception as e:
        print(f"ERROR processing {unique_og}: {e}")
        print(traceback.format_exc())
    return results


def write_tree(resolved_trees_working_dir, hog_name, subtree):
    try:
        tree_path = os.path.join(resolved_trees_working_dir, hog_name + ".txt")
        subtree.write(outfile=tree_path)
    except Exception as e:
        print(f"ERROR writing tree {hog_name}: {e}")

def write_fasta(align_id_dir, hog_name, sequences):
    try:
        fasta_path = os.path.join(align_id_dir, hog_name + ".fa")
        sorted_seqs = sorted(
            sequences.keys(), 
            key=lambda x: list(map(int, x.split("_"))) if "_" in x else x
        )
        with open(fasta_path, 'w') as outFile:
            for gene in sorted_seqs:
                outFile.write(f">{gene}\n")
                outFile.write(sequences[gene])
    except Exception as e:
        print(f"ERROR writing FASTA for {hog_name}: {e}")

def post_ogs_processing(
    unique_ogs,
    resolved_trees_working_dir, 
    resolved_trees_working_dir2, 
    hog_n0_over4genes, 
    name_dict, 
    species_names, 
    nprocess,
    align_id_dir=None,
    align_id_dir2=None,
):

    tree_file_index = index_files(resolved_trees_working_dir2, ".txt")
    fasta_file_index = index_files(align_id_dir2, ".fa") if align_id_dir2 is not None else {}

    hog_index = {
        unique_og: 
        [row for row in hog_n0_over4genes 
        if unique_og in row["OG"]] 
        for unique_og in unique_ogs
    }
    
    results_by_hog = []

    with ProcessPoolExecutor(max_workers=nprocess) as executor:
        futures = {}
        for unique_og in unique_ogs:
            hog_entries = hog_index.get(unique_og, [])
            if not hog_entries:
                continue
            future = executor.submit(
                process_unique_og,
                unique_og,
                hog_entries,
                species_names,
                name_dict,
                tree_file_index,
                fasta_file_index
            )
            futures[future] = unique_og
        
        for future in as_completed(futures):
            unique_og = futures[future]
            try:
                res = future.result()
                results_by_hog.extend(res)
            except Exception as e:
                print(f"ERROR in processing {unique_og}: {e}")
    
    for hog_name, subtree, pruned_alignments in results_by_hog:
        write_tree(resolved_trees_working_dir, hog_name, subtree)
        if align_id_dir is not None and pruned_alignments is not None:
            write_fasta(align_id_dir, hog_name, pruned_alignments)