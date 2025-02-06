import os
import traceback
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, as_completed
import ete3

def extract_gene_name(leaf_name, species_names):
    # Find the matching species name at the start of leaf_name
    matching_species = next((species for species in species_names if leaf_name.startswith(species + "_")), None)
    # If a matching species is found, extract the gene name
    if matching_species:
        # Remove the species prefix and underscore to get the gene name
        gene_name = leaf_name[len(matching_species) + 1:]  # +1 to account for the underscore
        return gene_name
    else:
        return None  # or handle as appropriate if no match is found
    
def get_file_path(unique_og, id_dir, extension=".fa"):
    for entry in os.scandir(id_dir):
        if entry.is_file() and entry.name == unique_og + extension:
            return entry.path

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


# Function to extract gene names from species_gene format
def extract_gene_name(leaf_name, species_names):
    # Find the matching species name at the start of leaf_name
    matching_species = next((species for species in species_names if leaf_name.startswith(species + "_")), None)
    # If a matching species is found, extract the gene name
    if matching_species:
        # Remove the species prefix and underscore to get the gene name
        gene_name = leaf_name[len(matching_species) + 1:]  # +1 to account for the underscore
        return gene_name
    else:
        return None  # or handle as appropriate if no match is found


def read_tree_and_fasta(unique_og, resolved_trees_working_dir, id_dir, read_queue, tree_extension=".txt", fasta_extension=".fa"):
    try:
        tree_id_path = get_file_path(unique_og, resolved_trees_working_dir, tree_extension)
        gene_tree = None
        if tree_id_path:
            with open(tree_id_path, "r") as file:
                gene_tree = ete3.Tree(file.read().strip(), quoted_node_names=True, format=1)
        else:
            print(f"WARNING: Tree file not found for {unique_og}")
        
        gene_dict = None
        if id_dir is not None:
            fasta_id_path = get_file_path(unique_og, id_dir, fasta_extension)
            if fasta_id_path:
                gene_dict = read_fasta(fasta_id_path)
            else:
                print(f"WARNING: FASTA file not found for {unique_og}")

        read_queue.put((unique_og, gene_tree, gene_dict))

    except Exception as e:
        print(f"ERROR reading data for {unique_og}: {e}")


def process_tree_id(hog_n0_over4genes, name_dict, species_names, read_queue, process_queue):
    while True:
        task = read_queue.get()
        if task is None: 
            process_queue.put(None)
            break

        unique_og, gene_tree, gene_dict = task
        try:
            mini_hog = [row for row in hog_n0_over4genes if unique_og in row["OG"]]
            results = []

            for row in mini_hog:
                hog_name = name_dict.get(row["HOG"], row["HOG"])
                parent_node = row["Gene Tree Parent Clade"]

                if parent_node == "n0":
                    subtree = gene_tree.copy()
                else:
                    subtree = gene_tree.get_tree_root().search_nodes(name=parent_node)
                    if not subtree:
                        print(f"WARNING: Parent node '{parent_node}' not found for {unique_og}")
                        continue
                    subtree = subtree[0]

                current_leaves = [leaf.name for leaf in subtree.get_leaves()]
                expected_leaves = ', '.join([row[col] for col in species_names if row[col]]).split(', ')

                if set(current_leaves) != set(expected_leaves):
                    subtree.prune(expected_leaves)

                pruned_alignments = None    
                if gene_dict is not None:
                    pruned_alignments = {
                            gene: gene_dict[gene] 
                            for gene in expected_leaves 
                            if gene in gene_dict
                        }

                results.append((hog_name, subtree, pruned_alignments))

            process_queue.put(results)
        except Exception as e:
            print(f"ERROR processing tree for {unique_og}: {e}")
            print(traceback.format_exc(e))

def write_tree_id(resolved_trees_working_dir, hog_name, subtree):
    try:
        tree_id_path = os.path.join(resolved_trees_working_dir, hog_name + ".txt")
        subtree.write(outfile=tree_id_path)

    except Exception as e:
        print(f"ERROR writing tree {hog_name}: {e}")

def write_to_fasta(id_output_dir, hog_name, sequences):
    pruned_seq_id_path = os.path.join(id_output_dir, hog_name + ".fa")
    sorted_seqs = sorted([*sequences.keys()], key=lambda x: list(map(int, x.split("_"))))
    with open(pruned_seq_id_path, 'w') as outFile:
        for gene in sorted_seqs:
            outFile.write(">%s\n" % gene)
            outFile.write(sequences[gene])

def write_to_files(
        process_queue,
        resolved_trees_working_dir, 
        align_id_dir, 
    ):
    while True:
        task = process_queue.get()
        if task is None: 
            break
        
        for hog_name, subtree, pruned_alignments in task:
            try:
                write_tree_id(resolved_trees_working_dir, hog_name, subtree)

                if align_id_dir is not None:
                    write_to_fasta(align_id_dir, hog_name, pruned_alignments)

            except Exception as e:
                print(f"ERROR writing tree or sequences for {hog_name}: {e}")


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

    read_queue = mp.Queue(maxsize=2 * nprocess)  
    process_queue = mp.Queue(maxsize=nprocess)  

    def start_reading():
        with ThreadPoolExecutor(max_workers=2) as reader_executor:
            reader_futures = [
                reader_executor.submit(
                    read_tree_and_fasta, 
                    unique_og, 
                    resolved_trees_working_dir2,
                    align_id_dir2, 
                    read_queue
                )
                for unique_og in unique_ogs
            ]

            for future in as_completed(reader_futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"ERROR in reading task: {e}")
        for _ in range(nprocess):
            read_queue.put(None)

    def start_processing():
        processor_workers = []
        for _ in range(nprocess):
            p = mp.Process(
                target=process_tree_id, 
                args=(hog_n0_over4genes, name_dict, species_names, read_queue, process_queue)
            )
            p.start()
            processor_workers.append(p)
        for p in processor_workers:
            p.join()
        process_queue.put(None)

    def start_writing():
        with ThreadPoolExecutor(max_workers=1) as writer_executor:
            writer_future = writer_executor.submit(
                write_to_files, 
                process_queue,
                resolved_trees_working_dir,
                align_id_dir, 
            )
            try:
                writer_future.result()
            except Exception as e:
                print(f"ERROR in writing task: {e}")
                print(traceback.format_exc(e))

    with ThreadPoolExecutor(max_workers=3) as executor:
        executor.submit(start_reading)
        executor.submit(start_processing)
        executor.submit(start_writing)
