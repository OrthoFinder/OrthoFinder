import time
import numpy as np
from collections import defaultdict
import itertools
import multiprocessing as mp
import warnings
try: 
    import queue
except ImportError:
    import Queue as queue
try:
    from rich import print
except ImportError:
    ...

from ..phylo import tree
from ..species_tree import stag
from ..utils import util, file_handler, parallel_task_manager, program_caller
from ..utils import blast_file_processor as blast_file_processor


class DendroBLASTTrees(object):
    def __init__(self, og_set, n_processes_alg, n_process_std, q_double_blast):
        self.og_set = og_set
        self.n_processes = n_processes_alg
        self.n_process_std = n_process_std
        self.q_double_blast = q_double_blast
        # Check files exist
    
    def tree_filename_ids(self, iog):
        return file_handler.handler.get_ogs_tree_fn(iog)
        
    def get_og_matrices_full_parallel(self):
        """
        read the blast files as well, remove need for intermediate pickle and unpickle
        ogMatrices contains matrix M for each OG where:
            Mij = 0.5*max(Bij, Bmin_i)/Bmax_i
        """
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            ogs_all = self.og_set.get_ogs_all()
            # DendroBLAST cannot be applied for --fast-add since we don't have hits between all genes, therefore can
            # assume that orthogroups are ordered
            lengths = [len(og) for og in ogs_all]
            for iog in range(1, len(ogs_all)):
                assert lengths[iog] <= lengths[iog-1], iog
            try:
                iog4 = next(i for i, og in enumerate(ogs_all) if len(og) < 4)
            except:
                iog4 = len(ogs_all)
            ogs = ogs_all[:iog4]
            ogs_per_species = [[[(g, i) for i, g in enumerate(og) if g.i_sp == i_sp] for i_sp in self.og_set.seqs_info.species_to_use] for og in ogs]
            n_genes = [len(og) for og in ogs]
            n_seqs = self.og_set.seqs_info.n_seqs_per_species
            og_matrices = [[mp.Array('d', n, lock=False) for _ in range(n)] for n in n_genes]
            blast_dir_list = file_handler.handler.get_blast_results_dir()
            cmd_queue = mp.Queue()
            for ii_sp, sp1 in enumerate(self.og_set.seqs_info.species_to_use):
                cmd_queue.put((ii_sp, sp1, n_seqs[sp1]))
            worker_status_queue = mp.Queue()
            # Should use PTM?
            running_processes = [mp.Process(target=worker_og_matrices_read_blast_and_update_distances,
                                           args=(cmd_queue, worker_status_queue, i_worker, og_matrices, n_genes,
                                                 self.og_set.seqs_info, blast_dir_list, ogs_per_species, self.q_double_blast))
                                for i_worker in range(self.n_processes)]
            for proc in running_processes:
                proc.start()
            rota = [None for i_worker in range(self.n_processes)]
            unfinished = []
            while True:
                # get process alive/dead
                time.sleep(1)
                alive = [proc.is_alive() for proc in running_processes]
                # read latest updates from queue, update rota
                try:
                    while True:
                        status, i_worker, i_task = worker_status_queue.get(True, 0.1)
                        if status == "start":
                            rota[i_worker] = i_task
                        elif status == "finish":
                            rota[i_worker] = None
                        elif status == "empty":
                            rota[i_worker] = "empty"
                except queue.Empty:
                    pass
                # A worker that died without finishing its task: stop the
                # others and raise, so main() records the error in checkpoint.txt.
                for proc, al, r in zip(running_processes, alive, rota):
                    if (not al) and (r != "empty"):
                        parallel_task_manager.terminate_processes(running_processes)
                        cmd_queue.cancel_join_thread()
                        raise parallel_task_manager.WorkerError(
                            "DendroBLAST worker (PID %s) exited with status %s "
                            "before finishing species %s.\n%s"
                            % (proc.pid, proc.exitcode, r, get_ram_error_text())
                        )
                if not any(alive):
                    break
#                print("WARNING: Computer ran out of RAM and killed OrthoFinder processes")
#                print("OrthoFinder will attempt to run these processes once more. If it is")
#                print("unsuccessful again then it will have to exit. Consider using")
#                print("the option '-a 1' or running on a machine with more RAM")
            #ogMatrices = [np.matrix(m) for m in ogMatrices]
            return ogs, og_matrices      
                   
    def complete_og_matrices(self, ogs, og_matrices):
        new_matrices = []
        for iog, (og, m) in enumerate(zip(ogs, og_matrices)):
            # dendroblast scores
            n = m.shape[0]
            m2 = np.zeros(m.shape)
            max_og = -9e99
            for i in range(n):
                for j in range(i):
                    m2[i, j] = -np.log(m[i,j] + m[j,i])  
                    m2[j, i] = m2[i, j]  
                    max_og = max(max_og, m2[i,j])
            new_matrices.append(m2)
        return new_matrices
        
    def complete_and_write_og_matrices(self, ogs, og_matrices):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        for iog, (og, m) in enumerate(zip(ogs, og_matrices)):
            # dendroblast scores
            n = len(m)
            max_og = -9e99
            # Careful not to over-write a value and then attempt to try to use the old value
            for i in range(n):
                for j in range(i):
                    m[i][j] = -np.log(m[i][j] + m[j][i])  
                    m[j][i] = m[i][j]  
                    max_og = max(max_og, m[i][j])
            self.write_phylip_matrix(m, [g.to_string() for g in og], file_handler.handler.get_ogs_dist_mat_fn(iog), max_og)
        return og_matrices
    
    @staticmethod
    def write_phylip_matrix(m, names, out_fn, max_og):
        """
        m - list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        max_og = 1.1*max_og
        sliver = 1e-6
        with open(out_fn, 'w') as outfile:
            n = len(m)
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(names[i] + " ")
                # values could be -inf, these are the most distantly related so replace with max_og
                V = [0. + (0. if i==j else m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")
    
    def species_tree_distances(self, ogs, og_matrices, method = 0):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        sp_pairs = list(itertools.combinations(self.og_set.seqs_info.species_to_use, 2))
        D = [[] for _ in sp_pairs]
        if method == 0:
            """ closest distance for each species pair in each orthogroup"""
            for og, m in zip(ogs, og_matrices):
                sp_dict = defaultdict(list)
                for i, g in enumerate(og):
                    sp_dict[g.i_sp].append(i)
                for (sp1, sp2), d_list in zip(sp_pairs, D):
                    distances = [m[i][j] for i in sp_dict[sp1] for j in sp_dict[sp2]]
                    if len(distances) > 0: d_list.append(min(distances))
#                    d_list.append(min(distances) if len(distances) > 0 else None)
        return D, sp_pairs
    
    def prepare_species_tree_command(self, D, sp_pairs, q_put_in_working_dir=False):
        n = len(self.og_set.seqs_info.species_to_use)
        M = np.zeros((n, n))
        for (sp1, sp2), d in zip(sp_pairs, D):
            sp1 = self.og_set.seqs_info.species_to_use.index(sp1)
            sp2 = self.og_set.seqs_info.species_to_use.index(sp2)
            x = np.median(d)
            M[sp1, sp2] = x
            M[sp2, sp1] = x
        species_matrix_fn = file_handler.handler.get_species_tree_matrix_fn(q_put_in_working_dir)  
        sliver = 1e-6
        with open(species_matrix_fn, 'w') as outfile:
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(str(self.og_set.seqs_info.species_to_use[i]) + " ")
                V = [(0. + M[i,j]) for j in range(n)]  # hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")       
        tree_fn = file_handler.handler.get_species_tree_unrooted_fn()
        cmd = " ".join(["fastme", "-i", species_matrix_fn, "-o", tree_fn, "-N", "-w", "O"] + (["-s"] if n < 1000 else []))
        return cmd, tree_fn
                
    def prepare_gene_tree_command(self):
        cmds = []
        ogs = self.og_set.get_ogs_all()
        for iog in range(len(ogs)):
            n_taxa = len(ogs[iog])
            if n_taxa < 4:
                continue
            cmds.append([" ".join(["fastme", "-i", file_handler.handler.get_ogs_dist_mat_fn(iog), "-o",
                                   file_handler.handler.get_ogs_tree_fn(iog), "-N", "-w", "O"] + (["-s"] if n_taxa < 1000 else []))])
        return cmds

    @staticmethod    
    def enough_ogs_for_stag(ogs, species_to_use):
        n_sp = len(species_to_use)
        n_sp_per_og = [len(set([g.i_sp for g in og])) for og in ogs]
        return (n_sp_per_og.count(n_sp) >= 100)
    
    def run_analysis(self, q_species_tree=True):
        """
        Args:
            qSpeciesTree - Bool: infer the species tree
        """
        util.print_underline("Calculating gene distances")
        ogs, og_matrices_partial = self.get_og_matrices_full_parallel()
        og_matrices = self.complete_and_write_og_matrices(ogs, og_matrices_partial)
        del og_matrices_partial
        util.print_time("Done")
        cmds_trees = self.prepare_gene_tree_command()
        q_less_than_four_species = len(self.og_set.seqs_info.species_to_use) < 4
        if not q_species_tree:
            q_stag = False
        elif q_less_than_four_species:
            q_stag = False
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
            write_species_tree_ids_two_three(self.og_set.seqs_info.species_to_use, sp_tree_fn_ids)
        else:
            q_stag = self.enough_ogs_for_stag(ogs, self.og_set.seqs_info.species_to_use)
            if not q_stag:
                print("Using fallback species tree inference method")
                D, sp_pairs = self.species_tree_distances(ogs, og_matrices)
                cmd_sp_tree, sp_tree_fn_ids = self.prepare_species_tree_command(D, sp_pairs)
                cmds_trees = [[cmd_sp_tree]] + cmds_trees
        del og_matrices
        util.print_underline("Inferring gene and species trees" if q_species_tree else "Inferring gene trees")
        program_caller.run_parallel_commands(self.n_process_std, cmds_trees, q_list_of_list=True)
        if q_stag:
            # Trees must have been completed
            print("")
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
            stag.run_for_ortho_finder(file_handler.handler.get_ogs_tree_dir(), file_handler.handler.get_working_directory_write(),
                                    self.og_set.seqs_info.species_to_use, sp_tree_fn_ids)
        if q_species_tree:
            util.rename_tree_taxa(sp_tree_fn_ids, file_handler.handler.get_species_tree_unrooted_fn(True), self.og_set.species_dict(),
                                q_support=False, q_fix_negatives=True)
            return sp_tree_fn_ids, q_stag
        else:      
            return None, q_stag
            
    def species_tree_only(self):
        q_less_than_four_species = len(self.og_set.seqs_info.species_to_use) < 4
        if q_less_than_four_species:
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
            write_species_tree_ids_two_three(self.og_set.seqs_info.species_to_use, sp_tree_fn_ids)
        else:
            ogs, og_matrices_partial = self.get_og_matrices_full_parallel()
            og_matrices = self.complete_og_matrices(ogs, og_matrices_partial)
            del og_matrices_partial
            D, sp_pairs = self.species_tree_distances(ogs, og_matrices)
            del og_matrices
            cmd_sp_tree, sp_tree_fn_ids = self.prepare_species_tree_command(D, sp_pairs, True)
            parallel_task_manager.run_command(cmd_sp_tree, True, False, raise_on_error=True)
        sp_tree_unrooted_fn = file_handler.handler.get_species_tree_unrooted_fn(True) 
        util.rename_tree_taxa(sp_tree_fn_ids, sp_tree_unrooted_fn, self.og_set.species_dict(), q_support=False, q_fix_negatives=True)  
        return sp_tree_fn_ids

# ==============================================================================================================================      
# DendroBlast   

def worker_og_matrices_read_blast_and_update_distances(cmd_queue, worker_status_queue, i_worker, og_matrices, n_genes, seqs_info,
                                                  blast_dir_list, ogs_per_species, q_double_blast):
    species_to_use = seqs_info.species_to_use
    og_genes = None   # per orthogroup and species: (gene iSeq array, matrix index array)
    with np.errstate(divide='ignore'):
        while True:
            try:
                ii_sp, sp1, n_seqs_sp1 = cmd_queue.get(True, 1)
                worker_status_queue.put(("start", i_worker, ii_sp))
                bs = [blast_file_processor.get_blast6_scores(seqs_info, blast_dir_list, sp1, sp2,
                                                         q_exclude_self_hits = False, q_double_blast=q_double_blast,
                                                         fmt="csr")
                      for sp2 in species_to_use]
                mins = np.ones((n_seqs_sp1, 1), dtype=np.float64)*9e99 
                maxes = np.zeros((n_seqs_sp1, 1), dtype=np.float64)
                for B in bs:
                    m0, m1 = csr_minmax(B)
                    mins = np.minimum(mins, m0)
                    maxes = np.maximum(maxes, m1)
                maxes_inv = 1./maxes
                # m[i][j] = 0.5*max(B[gi.iSeq, gj.iSeq], mins[gi.iSeq]) * maxes_inv[gi.iSeq]
                # for every gene pair of every orthogroup, computed with one
                # sparse lookup per species pair (all orthogroups together)
                # instead of one per gene pair.
                if og_genes is None:
                    og_genes = [
                        [(np.array([int(g.i_seq) for g, _ in og[s]], dtype=np.int64),
                          np.array([k for _, k in og[s]], dtype=np.intp))
                         for s in range(len(og))]
                        for og in ogs_per_species
                    ]
                for jj_sp, B in enumerate(bs):
                    blocks, rows_a, cols_b = [], [], []
                    for k, genes in enumerate(og_genes):
                        a, i_idx = genes[ii_sp]
                        b, j_idx = genes[jj_sp]
                        if a.size and b.size:
                            rows_a.append(np.repeat(a, b.size))
                            cols_b.append(np.tile(b, a.size))
                            blocks.append((k, i_idx, j_idx))
                    if not blocks:
                        continue
                    A = np.concatenate(rows_a)
                    scores = np.asarray(B[A, np.concatenate(cols_b)]).ravel()
                    values = 0.5 * np.maximum(scores, mins[A, 0]) * maxes_inv[A, 0]
                    pos = 0
                    for k, i_idx, j_idx in blocks:
                        block = values[pos:pos + i_idx.size * j_idx.size].reshape(i_idx.size, j_idx.size)
                        pos += i_idx.size * j_idx.size
                        m = og_matrices[k]
                        for i, row_values in zip(i_idx, block):
                            np.ctypeslib.as_array(m[i])[j_idx] = row_values
                del bs, B, mins, maxes, m0, m1, maxes_inv    # significantly reduces RAM usage
                worker_status_queue.put(("finish", i_worker, ii_sp))
            except queue.Empty:
                worker_status_queue.put(("empty", i_worker, None))
                return 

def get_ram_error_text():
    text = "ERROR: The computer ran out of RAM and killed OrthoFinder processes\n"
    text += "Try using a computer with more RAM. If you used the '-a' option\n"
    text += "it may be possible to complete the run by removing this option."
    return text

# ==============================================================================================================================

def csr_minmax(M):
    """Per-row minimum and maximum of the stored values (rows without values: 9e99 and 0)."""
    M = M.tocsr()
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    maxes = np.zeros((n, 1), dtype = np.float64)
    has = np.diff(M.indptr) > 0
    if has.any():
        starts = M.indptr[:-1][has]
        mins[has, 0] = np.minimum.reduceat(M.data, starts)
        maxes[has, 0] = np.maximum.reduceat(M.data, starts)
    return mins, maxes

def lil_minmax(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    maxes = np.zeros((n, 1), dtype = np.float64)
    for k_row in range(n):
        values=M.getrowview(k_row)
        if values.nnz == 0:
            continue
        mins[k_row] = min(values.data[0])
        maxes[k_row] = max(values.data[0])
    return mins, maxes

# ==============================================================================================================================    
# Species trees for two- & three-species analyses

def write_species_tree_ids_two_three(taxa, out_fn):
    """
    Get the unrooted species tree for two or three species
    Args:
        taxa - list of species names
    Returns:
    
    """
    t = tree.Tree()
    for s in taxa:
        t.add_child(tree.TreeNode(name=s))
    t.write(outfile=out_fn)
    
def get_species_tree_root_two_taxa(taxa):
    species_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
    t = tree.Tree("(%s,%s);" % (taxa[0], taxa[1]))  
    t.write(outfile=species_tree_fn_ids)
    return species_tree_fn_ids

 