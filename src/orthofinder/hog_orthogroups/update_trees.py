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
from ..utils import util, parallel_task_manager


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
    sequence = []
    try:
        with open(file_path, 'r') as fastaFile:
            for line in fastaFile:
                if line[0] == ">":
                    if not qFirst:
                        genes_dict[accession] = "".join(sequence)
                        sequence = []
                    qFirst = False
                    accession = line[1:].rstrip()
                else:
                    sequence.append(line)
            genes_dict[accession] = "".join(sequence)
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
        raise


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

class TreeLeafNames(object):
    """
    Sequence IDs for the gene names on the leaves of the resolved gene trees.

    A name of one gene maps to its ID. A name shared by several genes (the
    same name after ID extraction, or the same once Newick-unsafe characters
    are replaced by "_") is resolved per tree: its leaf gets the gene of that
    name that is in the tree's orthogroup (from the converted N0 HOGs, which
    hold the exact IDs), each gene once. Behaves as a dict of the
    unambiguous names otherwise (in, get).
    """

    def __init__(self, name_to_id, ambiguous, og_ids):
        self.name_to_id = name_to_id     # name -> ID, names of one gene
        self.ambiguous = ambiguous       # name -> [IDs], names of several genes
        self.og_ids = og_ids             # OG -> set of the IDs in it

    def __contains__(self, name):
        return name in self.name_to_id or name in self.ambiguous

    def get(self, name, default=None):
        return self.name_to_id.get(name, default)

    def leaf_ids(self, unique_og, names):
        """The ID of each leaf name of one orthogroup's tree (None if unknown)."""
        in_og = self.og_ids.get(unique_og, set())
        used = set()
        ids = []
        for name in names:
            seq_id = self.name_to_id.get(name)
            if seq_id is None and name in self.ambiguous:
                seq_id = next((i for i in self.ambiguous[name] if i in in_og and i not in used), None)
            if seq_id is not None:
                used.add(seq_id)
            ids.append(seq_id)
        return ids


def update_leaves(unique_og, gene_tree, spec_seq_id_dict=None):
    leaves = [leaf for leaf in gene_tree.leaves() if leaf.name is not None and leaf.name.strip()]
    if len(leaves) != len(list(gene_tree.leaves())):
        print(f"Warning: Null or empty leaf name in tree {unique_og}")
    if spec_seq_id_dict is None:
        return gene_tree
    names = [leaf.name for leaf in leaves]
    if isinstance(spec_seq_id_dict, TreeLeafNames):
        ids = spec_seq_id_dict.leaf_ids(unique_og, names)
    else:
        ids = [spec_seq_id_dict.get(name) for name in names]
    for leaf, name, seq_id in zip(leaves, names, ids):
        if seq_id is None:
            print(f"Warning: Leaf name '{name}' not found in mapping dictionary for {unique_og}")
        else:
            leaf.name = seq_id
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

def _put(q, item, stop_event, timeout=1.0):
    """Put with back-pressure, giving up (False) once stop_event is set."""
    while not stop_event.is_set():
        try:
            q.put(item, timeout=timeout)
            return True
        except queue.Full:
            continue
    return False


def _report_error(report_queue, stop_event, where):
    stop_event.set()
    try:
        report_queue.put(("error", "%s:\n%s" % (where, traceback.format_exc())))
    except Exception:
        pass


def process_task(
        read_queue, 
        process_queue, 
        hog_index, 
        name_dict, 
        species_names, 
        stop_event,
        report_queue,
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
                _put(process_queue, ("skip", unique_og, "no_gene_tree"), stop_event)
                continue

            if not hog_entries:
                _put(process_queue, ("skip", unique_og, "no_hog_entries"), stop_event)
                continue

            hog_entries = sorted(hog_entries, key=lambda r: str(r.get("HOG", "")))

            # One pass over the tree instead of a search per HOG row. setdefault
            # keeps the first node in traversal order, as search_nodes did.
            nodes_by_name = {}
            for node in gene_tree.traverse():
                if node.name:
                    nodes_by_name.setdefault(node.name, node)

            for row in hog_entries:
                hog_name = name_dict.get(row.get("HOG"), row.get("HOG"))
                parent_node = row.get("Gene Tree Parent Clade")

                if not parent_node:
                    continue

                if parent_node == "n0":
                    subtree = gene_tree.copy()
                else:
                    node = nodes_by_name.get(parent_node)
                    if node is None:
                        continue
                    subtree = node.copy()

                current_leaves = {leaf.name for leaf in subtree.leaves() if leaf.name}

                expected_leaves = set()
                for col in species_names:
                    v = row.get(col)
                    if not v:
                        continue
                    for x in str(v).split(","):
                        x = x.strip()
                        if x:
                            expected_leaves.add(x)

                if not expected_leaves:
                    continue

                valid_leaves = sorted(expected_leaves & current_leaves)
                if len(valid_leaves) < 2:
                    continue
                try:
                    subtree.prune(valid_leaves)
                except Exception:
                    if strict_prune_fail:
                        raise
                    continue
               
                pruned_alignments = None
                if gene_dict:
                    pruned_alignments = {g: gene_dict[g] for g in valid_leaves if g in gene_dict}
                try:
                    newick = subtree.write(outfile=None, parser=5)
                except Exception:
                    if strict_prune_fail:
                        raise
                    continue
                results.append((hog_name, newick, pruned_alignments))

            if not results:
                msg = ("skip", unique_og, "no_outputs_from_hog_entries")
            else:
                msg = ("og", unique_og, results)
            if not _put(process_queue, msg, stop_event):
                break
          
    except Exception:
        _report_error(report_queue, stop_event, "tree processing")
        raise
    finally:
        if stop_event.is_set():
            # Do not block at exit on data nobody will read.
            process_queue.cancel_join_thread()


def writer_task(
        process_queue, 
        min_seq, 
        idDict,
        resolved_trees_id_dir, 
        align_dir, 
        stop_event, 
        report_queue,
        exist_msa=True
    ):
    n_handled = 0
    try:
        while not stop_event.is_set():
            try:
                msg = process_queue.get(timeout=1.0)
            except queue.Empty:
                continue

            if msg is None:
                report_queue.put(("done", n_handled))
                return

            kind = msg[0]
            if kind == "og":
                _, unique_og, results = msg
                for out_name, newick_string, pruned_alignments in results:
                    if not write_tree(out_name, newick_string, resolved_trees_id_dir):
                        raise RuntimeError("Failed to write tree %s" % out_name)

                    if exist_msa and pruned_alignments is not None and len(pruned_alignments) >= min_seq:
                        if align_dir is not None:
                            write_fasta(align_dir, out_name, pruned_alignments, idDict)
            elif kind != "skip":
                raise TypeError("Unexpected message: %r" % (msg,))
            n_handled += 1

    except Exception:
        _report_error(report_queue, stop_event, "writing trees/alignments")
        raise


def threaded_reader(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index, n_threads, stop_event, report_queue, exist_msa=True):
    try:
        def worker(unique_og):
            if stop_event.is_set():
                return
            task = read_files(unique_og, spec_seq_id_dict, tree_file_index, fasta_file_index, exist_msa=exist_msa)
            _put(read_queue, task, stop_event)
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            # Consume the results so an exception in any read is raised here.
            for _ in executor.map(worker, unique_ogs):
                pass
    except Exception:
        _report_error(report_queue, stop_event, "reading trees/alignments")
        raise
    finally:
        if stop_event.is_set():
            read_queue.cancel_join_thread()


def _put_sentinels(q, n_left):
    """Put up to n_left None sentinels without blocking; return how many remain."""
    while n_left:
        try:
            q.put_nowait(None)
        except queue.Full:
            break
        n_left -= 1
    return n_left


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

    # Bounded queues: the reader cannot run ahead and hold every tree and
    # alignment in memory at once.
    read_queue = mp.Queue(maxsize=max(4 * n_processor_processes, 16))
    process_queue = mp.Queue(maxsize=max(4 * n_writer_processes, 16))
    report_queue = mp.Queue()
    stop_event = mp.Event()

    file_reader = mp.Process(
        target=threaded_reader,
        args=(read_queue, unique_ogs, spec_seq_id_dict, tree_file_index, fasta_file_index,
              n_reader_threads, stop_event, report_queue, exist_msa)
    )
    file_processors = [
        mp.Process(
            target=process_task,
            args=(read_queue, process_queue, hog_index, name_dict, species_names,
                  stop_event, report_queue)
        )
        for _ in range(n_processor_processes)
    ]
    writer_processes = [
        mp.Process(
            target=writer_task,
            args=(process_queue, min_seq, idDict, resolved_trees_id_dir,
                  align_dir, stop_event, report_queue, exist_msa)
        )
        for _ in range(n_writer_processes)
    ]
    all_processes = [file_reader] + file_processors + writer_processes
    for proc in all_processes:
        proc.start()

    n_expected = len(unique_ogs)
    n_handled = 0
    writers_done = 0
    processor_sentinels = None   # None: reader not finished yet
    writer_sentinels = None      # None: processors not finished yet
    failure = None

    try:
        while failure is None:
            try:
                while True:
                    msg = report_queue.get(timeout=0.2)
                    if msg[0] == "error":
                        failure = msg[1]
                        break
                    if msg[0] == "done":
                        n_handled += msg[1]
                        writers_done += 1
            except queue.Empty:
                pass
            if failure is not None:
                break

            for proc in all_processes:
                if proc.exitcode not in (None, 0):
                    failure = "process %d terminated with exit code %d" % (proc.pid, proc.exitcode)
                    break
            if failure is not None:
                break

            # Each stage is told to stop only after the previous one finished cleanly.
            if processor_sentinels is None and file_reader.exitcode == 0:
                processor_sentinels = n_processor_processes
            if processor_sentinels:
                processor_sentinels = _put_sentinels(read_queue, processor_sentinels)

            if (
                writer_sentinels is None
                and processor_sentinels == 0
                and all(p.exitcode == 0 for p in file_processors)
            ):
                writer_sentinels = n_writer_processes
            if writer_sentinels:
                writer_sentinels = _put_sentinels(process_queue, writer_sentinels)

            if writers_done == n_writer_processes:
                break

        if failure is None and n_handled != n_expected:
            failure = "only %d of %d orthogroups were processed" % (n_handled, n_expected)

    except BaseException:
        failure = failure or "interrupted"
        raise

    finally:
        if failure is not None:
            stop_event.set()
            parallel_task_manager.TerminateProcesses(all_processes)
            for q in (read_queue, process_queue, report_queue):
                q.cancel_join_thread()
        for proc in all_processes:
            proc.join()

    if failure is not None:
        # Raised (not Fail()) so main() records the error in checkpoint.txt.
        raise parallel_task_manager.WorkerError(
            "Updating gene trees and alignments failed: %s" % failure
        )
