import os
import io 
import numpy as np
import multiprocessing as mp
import queue
from concurrent.futures import ThreadPoolExecutor
import ete4
import tempfile
import traceback
import time
from ..utils import util


def write_tree(hog_name, newick_string, resolved_trees_id_dir):

    try:

        if resolved_trees_id_dir is None:
            raise ValueError("resolved_trees_id_dir is None")

        os.makedirs(resolved_trees_id_dir, exist_ok=True)
        tree_id_file = os.path.join(resolved_trees_id_dir, f"{hog_name}.txt")

        if newick_string is None:
            raise ValueError(f"Empty or None tree string for '{hog_name}'")

        if isinstance(newick_string, str):
            data = newick_string.encode("utf-8")
        elif isinstance(newick_string, bytes):
            data = newick_string
        else:
            raise TypeError(f"Unexpected type for newick_string: {type(newick_string)}")

        with tempfile.NamedTemporaryFile("wb", dir=resolved_trees_id_dir, delete=False) as tmp:
            tmp.write(data)
            tmp.flush()
            os.fsync(tmp.fileno())
            tmp_name = tmp.name

        os.replace(tmp_name, tree_id_file)
        return True

    except Exception as e:
        print(f"ERROR writing tree '{hog_name}': {e}\n{traceback.format_exc()}")
        return False


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


def read_files(unique_og, spec_seq_id_dict, tree_file_index, fasta_file_index, exist_msa=True):
   
    gene_tree = read_tree_file(unique_og, tree_file_index, spec_seq_id_dict)
    gene_dict = read_fasta_file(unique_og, fasta_file_index, exist_msa=exist_msa)
    return (unique_og, gene_tree, gene_dict)

def check_path(s):
    s = s.strip().strip('"').strip("'")
    
    return (
        os.path.exists(s) or
        '/' in s or '\\' in s or
        os.path.dirname(s) != ''
    )

def update_leaves(unique_og, gene_tree, spec_seq_id_dict=None):
    for leaf in gene_tree.leaves(): #.iter_leaves():
        original = leaf.name
        if original is None or not original.strip():
            print(f"Warning: Null or empty leaf name in tree {unique_og}")
            continue
        if spec_seq_id_dict is not None and original not in spec_seq_id_dict:
            print(f"Warning: Leaf name '{original}' not found in mapping dictionary for {unique_og}")
        if spec_seq_id_dict is not None:
            leaf.name = spec_seq_id_dict.get(original, original)
    return gene_tree

def read_tree_file(unique_og, tree_file_index, spec_seq_id_dict=None):
    gene_tree = None
    if unique_og in tree_file_index:
        try:
            if check_path(tree_file_index[unique_og]):
                with open(tree_file_index[unique_og], "r") as file:
                    tree_data = file.read().strip()
                    gene_tree = ete4.Tree(tree_data, parser=1) #quoted_node_names=True, format=1
            else:
                gene_tree = ete4.Tree(tree_file_index[unique_og], parser=1)
            gene_tree = update_leaves(unique_og, gene_tree, spec_seq_id_dict=spec_seq_id_dict)
        except Exception as e:
            print(f"ERROR reading tree for {unique_og}: {e}")
            raise
    else:
        print(f"WARNING: Tree file not found for {unique_og}")
    return gene_tree

def read_fasta_file(unique_og, fasta_file_index, exist_msa=True):
    gene_dict = {}
    if unique_og in fasta_file_index:
        try:
            gene_dict = read_fasta(fasta_file_index[unique_og])
        except Exception as e:
            print(f"ERROR reading FASTA for {unique_og}: {e}")
            raise
    else:
        if exist_msa:
            print(f"WARNING: FASTA file not found for {unique_og}")
    return gene_dict

def process_task(
        read_queue, 
        process_queue, 
        hog_index, 
        name_dict, 
        species_names, 
        stop_event,
        strict_prune_fail=False
    ):
    try:
        while not stop_event.is_set():
            try:
                task = read_queue.get(timeout=1)
            except queue.Empty:
                continue

            if task is None:
                break

            unique_og, gene_tree, gene_dict = task
            hog_entries = hog_index.get(unique_og, [])
            results = []

            if gene_tree is None:
                process_queue.put(("skip", unique_og, "no_gene_tree"))
                continue

            if not hog_entries:
                process_queue.put(("skip", unique_og, "no_hog_entries"))
                continue

            hog_entries = sorted(hog_entries, key=lambda r: str(r.get("HOG", "")))

            for row in hog_entries:
                hog_name = name_dict.get(row.get("HOG"), row.get("HOG"))
                parent_node = row.get("Gene Tree Parent Clade")

                if not parent_node:
                    continue

                if parent_node == "n0":
                    subtree = gene_tree.copy()
                else:
                    subtree_nodes = list(gene_tree.search_nodes(name=parent_node))
                    if not subtree_nodes:
                        continue
                    subtree = subtree_nodes[0].copy()

                current_leaves = [leaf.name for leaf in subtree.leaves() if leaf.name]

                expected_leaves = []
                for col in species_names:
                    v = row.get(col)
                    if not v:
                        continue
                    for x in str(v).split(","):
                        x = x.strip()
                        if x:
                            expected_leaves.append(x)

                expected_leaves = sorted(set(expected_leaves))
                if not expected_leaves:
                    continue
                
                valid_leaves = [leaf for leaf in expected_leaves if leaf in current_leaves]
                if len(valid_leaves) < 2:
                    continue
                valid_leaves = sorted(set(valid_leaves))
                try:
                    subtree.prune(valid_leaves)
                except Exception:
                    if strict_prune_fail:
                        process_queue.put(("error", unique_og, f"prune_failed:{hog_name}", traceback.format_exc()))
                        stop_event.set()
                        break
                    else:
                        continue
               
                pruned_alignments = None
                if gene_dict:
                    pruned_alignments = {g: gene_dict[g] for g in valid_leaves if g in gene_dict}
                try:
                    newick = subtree.write(outfile=None, parser=5)
                except Exception:
                    if strict_prune_fail:
                        process_queue.put(("error", unique_og, f"prune_failed:{hog_name}"))
                        stop_event.set()
                        break
                    else:
                        continue
                results.append((hog_name, newick, pruned_alignments))

            if not results:
                process_queue.put(("skip", unique_og, "no_outputs_from_hog_entries"))
            else:
                process_queue.put(("og", unique_og, results))

            if stop_event.is_set():
                break
          
    except Exception:
        process_queue.put(("error", None, "process_task_crash"))
        stop_event.set()
        raise


def writer_task(
        process_queue, 
        min_seq, 
        idDict,
        resolved_trees_id_dir, 
        align_dir, 
        stop_event, 
        exist_msa=True
    ):

    try:
        while not stop_event.is_set():
            try:
                msg = process_queue.get(timeout=1.0)
            except queue.Empty:
                continue

            if msg is None:
                break

            kind = msg[0]

            if kind == "skip":
                continue
            elif kind == "error":
                stop_event.set()
                continue
            elif kind != "og":
                continue

            _, unique_og, results = msg

            for out_name, newick_string, pruned_alignments in results:
                try:
                    ok = write_tree(out_name, newick_string, resolved_trees_id_dir)
                    if not ok:
                        raise RuntimeError("write_tree returned False")
                except Exception:
                    stop_event.set()
                    break

                if exist_msa and pruned_alignments is not None and len(pruned_alignments) >= min_seq:
                    if align_dir is not None:
                        write_fasta(align_dir, out_name, pruned_alignments, idDict)

    except Exception as e:
        stop_event.set()
        raise


def threaded_reader(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index, n_threads=4, stop_event=None, exist_msa=True):
    try:
        def worker(unique_og):
            if stop_event is not None and stop_event.is_set():
                return
            task = read_files(unique_og, spec_seq_id_dict, tree_file_index, fasta_file_index, exist_msa=exist_msa)
            read_queue.put(task)
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            executor.map(worker, unique_ogs)
    except Exception:
        if stop_event is not None:
            stop_event.set()
        raise

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
    min_seq=4,
    exist_msa=True,
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
    report_queue = mp.Queue()

    # Start reader
    file_reader = mp.Process(
        target=threaded_reader,
        args=(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index,
              n_reader_threads, stop_event, exist_msa)
    )
    file_reader.start()

    # Start processors
    file_processors = []
    for _ in range(n_processor_processes):
        p = mp.Process(
            target=process_task,
            args=(read_queue, process_queue, hog_index, name_dict, species_names,
                  stop_event)
        )
        p.start()
        file_processors.append(p)

    # Start writers
    writer_processes = []
    for _ in range(n_writer_processes):
        w = mp.Process(
            target=writer_task,
            args=(process_queue, min_seq, idDict, resolved_trees_id_dir,
                  align_dir, stop_event, exist_msa)
        )
        w.start()
        writer_processes.append(w)

    all_processes = [file_reader] + file_processors + writer_processes

    # Join reader
    try:
        file_reader.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt detected during file reading. Initiating shutdown.", flush=True)
        stop_event.set()
        file_reader.terminate()
        file_reader.join()

    if file_reader.exitcode not in (0, None):
        print(f"Reader process failed with exit code {file_reader.exitcode}. Initiating shutdown.", flush=True)
        stop_event.set()

    # Stop processors
    for _ in range(n_processor_processes):
        read_queue.put(None)
    for p in file_processors:
        p.join()

    # Stop writers
    for _ in range(n_writer_processes):
        process_queue.put(None)
    for w in writer_processes:
        w.join()

    for proc in all_processes:
        if proc.exitcode not in (0, None):
            print(f"ERROR: process {proc.pid} terminated with exit code {proc.exitcode}.", flush=True)
            util.Fail()

