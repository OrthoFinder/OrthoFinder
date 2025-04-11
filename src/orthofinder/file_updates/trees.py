import os
import io 
import numpy as np
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
import ete3


def write_tree(hog_name, subtree, resolved_trees_id_dir):
    try:
        tree_id_file = os.path.join(resolved_trees_id_dir, hog_name + ".txt")
        tree_string = subtree.write(outfile=None, format=5)
        with open(tree_id_file, "wb") as f:
            f.write(tree_string.encode("utf-8"))
    except Exception as e:
        print(f"ERROR writing tree {hog_name}: {e}")


def read_fasta(file_path):
    genes_dict = {}
    qFirst = True
    accession = ""
    sequence = ""
    try:
        with open(file_path, 'r') as fastaFile:
            for line in fastaFile:
                if line[0] == ">":
                    if not qFirst:
                        genes_dict[accession] = sequence
                        sequence = ""
                    qFirst = False
                    accession = line[1:].rstrip()
                else:
                    sequence += line
            genes_dict[accession] = sequence
    except Exception as e:
        print(f"ERROR reading FASTA file {file_path}: {e}")
        raise
    return genes_dict


def write_fasta(align_dir, hog_name, sequences, idDict):
    try:
        fasta_path = os.path.join(align_dir, hog_name + ".fa")
        sorted_seqs = sorted(
            sequences.keys(),
            key=lambda x: list(map(int, x.split("_"))) if "_" in x else x
        )
        buffer = io.StringIO()
        for gene in sorted_seqs:
            gene_name = idDict.get(gene)
            buffer.write(f">{gene_name}\n")
            buffer.write(sequences[gene])
        with open(fasta_path, 'w', buffering=1024 * 1024) as outFile:
            outFile.write(buffer.getvalue())
    except Exception as e:
        print(f"ERROR writing FASTA for {hog_name}: {e}")


def read_files(unique_og, spec_seq_id_dict, tree_file_index, fasta_file_index):
    gene_tree = None
    gene_dict = None

    if unique_og in tree_file_index:
        try:
            with open(tree_file_index[unique_og], "r") as file:
                tree_data = file.read().strip()
                gene_tree = ete3.Tree(tree_data, quoted_node_names=True, format=1)
                for leaf in gene_tree.iter_leaves():
                    original = leaf.name
                    if original is None or not original.strip():
                        print(f"Warning: Null or empty leaf name in tree {unique_og}")
                        continue
                    if original not in spec_seq_id_dict:
                        print(f"Warning: Leaf name '{original}' not found in mapping dictionary for {unique_og}")
                    leaf.name = spec_seq_id_dict.get(original, original)
        except Exception as e:
            print(f"ERROR reading tree for {unique_og}: {e}")
            raise
    else:
        print(f"WARNING: Tree file not found for {unique_og}")

    if unique_og in fasta_file_index:
        try:
            gene_dict = read_fasta(fasta_file_index[unique_og])
        except Exception as e:
            print(f"ERROR reading FASTA for {unique_og}: {e}")
            raise
    else:
        print(f"WARNING: FASTA file not found for {unique_og}")

    return (unique_og, gene_tree, gene_dict)


def process_task(read_queue, process_queue, hog_index, name_dict, species_names, stop_event):
    try:
        while not stop_event.is_set():
            try:
                task = read_queue.get(timeout=1)
            except Exception:
                continue

            if task is None:
                break

            unique_og, gene_tree, gene_dict = task
            hog_entries = hog_index.get(unique_og, [])
            results = []

            if gene_tree is None:
                print(f"Skipping {unique_og}: No gene tree loaded")
                continue

            for row in hog_entries:
                hog_name = name_dict.get(row.get("HOG"), row.get("HOG"))
                parent_node = row.get("Gene Tree Parent Clade")
                if not parent_node:
                    print(f"WARNING: Missing parent node for HOG {hog_name} in {unique_og}")
                    continue
                if parent_node == "n0":
                    subtree = gene_tree.copy()
                else:
                    subtree_nodes = gene_tree.search_nodes(name=parent_node)
                    if not subtree_nodes:
                        print(f"WARNING: Parent node '{parent_node}' not found in gene tree for {unique_og}, HOG {hog_name}")
                        continue
                    subtree = subtree_nodes[0]
                current_leaves = [leaf.name for leaf in subtree.get_leaves() if leaf.name]
                expected_leaves = [x.strip() for col in species_names if row.get(col) for x in row[col].split(',')]
                if not expected_leaves:
                    print(f"WARNING: No expected leaves found for {hog_name} in {unique_og}")
                    continue
                missing_in_tree = [name for name in expected_leaves if name not in current_leaves]
                if missing_in_tree:
                    print(f"INFO: Some expected leaves not found in subtree for {hog_name} ({unique_og}): {missing_in_tree}")
                try:
                    valid_leaves = [leaf for leaf in expected_leaves if leaf in current_leaves]
                    subtree.prune(valid_leaves)
                except Exception as e:
                    print(f"WARNING: Could not prune subtree for {hog_name} in {unique_og}: {e}")
                    raise

                pruned_alignments = None
                if gene_dict is not None:
                    pruned_alignments = {gene: gene_dict[gene] for gene in expected_leaves if gene in gene_dict}
                results.append((hog_name, subtree, pruned_alignments))

            process_queue.put(results)
    except Exception as e:
        print(f"ERROR in process_task: {e}")
        os._exit(1) 
    return


def writer_task(process_queue, min_seq, idDict, resolved_trees_id_dir, align_dir, stop_event):
    try:
        while not stop_event.is_set():
            try:
                task = process_queue.get(timeout=1)
            except Exception:
                continue

            if task is None:
                break

            for hog_name, subtree, pruned_alignments in task:
                if pruned_alignments is None:
                    continue
                if len(pruned_alignments) >= min_seq:
                    write_tree(hog_name, subtree, resolved_trees_id_dir)
                if align_dir is not None and pruned_alignments is not None:
                    write_fasta(align_dir, hog_name, pruned_alignments, idDict)
    except Exception as e:
        print(f"ERROR in writer_task: {e}")
        os._exit(1)
    return


def threaded_reader(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index, n_threads=4, stop_event=None):
    try:
        def worker(unique_og):
            if stop_event is not None and stop_event.is_set():
                return
            task = read_files(unique_og, spec_seq_id_dict, tree_file_index, fasta_file_index)
            read_queue.put(task)
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            executor.map(worker, unique_ogs)
    except Exception as e:
        print(f"ERROR in threaded_reader: {e}")
        os._exit(1)
    return


def post_ogs_processing(
    unique_ogs,  
    resolved_trees_id_dir,
    hog_index, 
    name_dict, 
    idDict,
    spec_seq_id_dict,
    species_names, 
    nprocess,
    tree_file_index,
    fasta_file_index,            
    align_dir=None,   
    min_seq=4       
):
    if nprocess >= 128:
        n_reader_threads = min(max(nprocess // 4, 2), 32)
        n_processor_processes = max(nprocess * 3 // 4, 1)
        n_writer_processes = min(max(nprocess // 4, 4), 32)
    else:
        n_reader_threads = max(min(nprocess // 2, 16), 4)
        n_processor_processes = max(nprocess // 2, 1)
        n_writer_processes = max(1, min(int(np.ceil(np.abs(nprocess // 2 - 1))), max(4, nprocess // 4)))
   
    process_queue = mp.Queue()
    read_queue = mp.Queue()
    stop_event = mp.Event()

    file_reader = mp.Process(
        target=threaded_reader,
        args=(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index, n_reader_threads, stop_event)
    )
    file_reader.start()

    file_processors = []
    for _ in range(n_processor_processes):
        p = mp.Process(
            target=process_task,
            args=(read_queue, process_queue, hog_index, name_dict, species_names, stop_event)
        )
        p.start()
        file_processors.append(p)

    writer_processes = []
    for _ in range(n_writer_processes):
        writer_process = mp.Process(
            target=writer_task,
            args=(process_queue, min_seq, idDict, resolved_trees_id_dir, align_dir, stop_event)
        )
        writer_process.start()
        writer_processes.append(writer_process)

    all_processes = [file_reader] + file_processors + writer_processes

    try:
        file_reader.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt detected during file reading. Initiating shutdown.")
        stop_event.set()
        file_reader.terminate()

    if file_reader.exitcode is not None and file_reader.exitcode != 0:
        print(f"Reader process failed with exit code {file_reader.exitcode}. Initiating shutdown.")
        stop_event.set()

    for _ in range(n_processor_processes):
        read_queue.put(None)
    for p in file_processors:
        p.join()

    for _ in range(n_writer_processes):
        process_queue.put(None)
    for p in writer_processes:
        p.join()

    for proc in all_processes:
        if proc.exitcode is not None and proc.exitcode != 0:
            print(f"Process {proc.pid} terminated with exit code {proc.exitcode}.")
