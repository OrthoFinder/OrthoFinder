import os
import shutil
import traceback
import numpy as np
from scipy import sparse
import warnings
import numpy.core.numeric as numeric
from scipy.optimize import curve_fit
import multiprocessing as mp
from ..utils import util, file_handler, blast_file_processor, matrices, parallel_task_manager

try:
    import queue
except ImportError:
    import Queue as queue
try:
    from rich import print
except ImportError:
    ...


"""
scnorm
-------------------------------------------------------------------------------
"""


class ScNorm:
    @staticmethod
    def loglinear(x, a, b):
        return a * np.log10(x) + b

    @staticmethod
    def get_length_arrays_for_matrix(m, len_i, len_j):
        # Row-major order, as the previous lil-based version produced.
        m = m.tocsr()
        m.eliminate_zeros()
        m.sort_indices()
        I = np.repeat(np.arange(m.shape[0]), np.diff(m.indptr))
        J = m.indices
        scores = m.data
        li = np.array(len_i[I])
        lj = np.array(len_j[J])
        return li, lj, scores

    @staticmethod
    def get_top_percentile_of_scores(L, S, percentile_to_keep):
        # Get the top x% of hits at each length
        L = np.asarray(L)
        S = np.asarray(S, dtype=float)
        n_scores = len(S)
        # Sort by length, ties by original position (as sorted(zip(L, range(n))) did).
        indices = np.argsort(L, kind="stable")
        s_sorted = S[indices]
        l_sorted = L[indices]
        if n_scores < 100:
            # then we can't split them into bins, return all for fitting
            return l_sorted, s_sorted
        n_in_bins = 1000 if n_scores > 5000 else (200 if n_scores > 1000 else 20)
        n_bins, remainder = divmod(n_scores, n_in_bins)
        # Full bins only (the remainder is not used), one percentile per bin.
        n = n_bins * n_in_bins
        binned_scores = s_sorted[:n].reshape(n_bins, n_in_bins)
        cut_offs = np.percentile(binned_scores, percentile_to_keep, axis=1)
        keep = (binned_scores >= cut_offs[:, None]).ravel()
        return l_sorted[:n][keep], s_sorted[:n][keep]

    @staticmethod
    def calculate_fitting_parameters(lf, S):
        pars, covar = curve_fit(ScNorm.loglinear, lf, np.log10(S))
        return pars

    @staticmethod
    def scale_rows_cols(B, row_vals, col_vals, factor=1.0):
        """
        Return the CSR matrix with entries factor * row_vals[i] * B[i,j] * col_vals[j],
        equivalent to factor * diag(row_vals) * B * diag(col_vals) but computed
        directly on the non-zero values in O(nnz).
        """
        B = sparse.csr_matrix(B, copy=True)
        B.sum_duplicates()
        rows = np.repeat(np.arange(B.shape[0]), np.diff(B.indptr))
        B.data *= row_vals[rows]
        B.data *= col_vals[B.indices]
        if factor != 1.0:
            B.data *= factor
        B.eliminate_zeros()
        return B

    @staticmethod
    def normalise_scores_by_log_length_product(b, lq, lh, params):
        li_vals = lq ** (-params[0])
        lj_vals = lh ** (-params[0])
        return ScNorm.scale_rows_cols(b, li_vals, lj_vals, 10 ** (-params[1]))


"""
WaterfallMethod
-------------------------------------------------------------------------------
"""


class WaterfallMethod:
    @staticmethod
    def normalise_scores(B, lengths, i_species, j_species):
        li, lj, scores = ScNorm.get_length_arrays_for_matrix(
            B, lengths[i_species], lengths[j_species]
        )
        lf = li * lj
        top_lf, top_scores = ScNorm.get_top_percentile_of_scores(lf, scores, 95)
        if len(top_scores) > 1:
            fitting_parameters = ScNorm.calculate_fitting_parameters(top_lf, top_scores)
            return ScNorm.normalise_scores_by_log_length_product(
                B, lengths[i_species], lengths[j_species], fitting_parameters
            )
        elif i_species == j_species and B.nnz == 0:
            print(
                "WARNING: THIS IS UNCOMMON, there are zero hits when searching the genes in species %d against itself. Check the input proteome contains all the genes from that species and check the search program is working (default is diamond)."
                % i_species
            )
            return sparse.csr_matrix(B.get_shape())
        else:
            print(
                "WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored"
                % (i_species, j_species)
            )
            return sparse.csr_matrix(B.get_shape())

    @staticmethod
    def normalised_bit_score(B, lengths, i_species, j_species):
        """
        Args:
            B - sparse matrix
        Returns
            B' - CSR matrix
        """
        if B.nnz == 0:
            return sparse.csr_matrix(B)
        lq = lengths[i_species]
        lh = lengths[j_species]
        li_vals = lq ** (-0.5)
        lj_vals = lh ** (-0.5)
        return ScNorm.scale_rows_cols(B, li_vals, lj_vals)

    @staticmethod
    def process_blast_hits(
        seqs_info,
        blast_dir_list,
        lengths,
        i_species,
        d_pickle,
        q_double_blast,
        v2_scores,
        q_allow_empty=False,
    ):
        """
        iLimitNewSpecies: int - Only process fasta files with OrthoFidner ID >= iLimitNewSpecies
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            bi = []
            for j_species in range(seqs_info.n_species):
                bij = blast_file_processor.get_blast6_scores(
                    seqs_info,
                    blast_dir_list,
                    seqs_info.species_to_use[i_species],
                    seqs_info.species_to_use[j_species],
                    q_double_blast=q_double_blast,
                    q_allow_empty=q_allow_empty,
                    fmt="csr",
                )
                if v2_scores:
                    bij = WaterfallMethod.normalised_bit_score(
                        bij, lengths, i_species, j_species
                    )
                else:
                    bij = WaterfallMethod.normalise_scores(
                        bij, lengths, i_species, j_species
                    )
                bi.append(bij)
            matrices.dump_matrix_array("B", bi, i_species, d_pickle)
            BH = WaterfallMethod.get_bh_s(bi, seqs_info, i_species)
            matrices.dump_matrix_array("BH", BH, i_species, d_pickle)
            # util.print_time("Initial processing of species %d complete" % iSpecies)


    @staticmethod
    def worker_process_blast_hits_new(
        seqs_info,
        blast_dir_list,
        lengths,
        cmd_queue,
        d_pickle,
        q_double_blast,
        v2_scores,
        q_allow_empty,
        result_queue,
    ):  
        try:
            while True:
                try:
                    i_species = cmd_queue.get(True, 1)
                    if i_species is None:
                        break
                    WaterfallMethod.process_blast_hits(
                        seqs_info,
                        blast_dir_list,
                        lengths,
                        i_species,
                        d_pickle=d_pickle,
                        q_double_blast=q_double_blast,
                        v2_scores=v2_scores,
                        q_allow_empty=q_allow_empty,
                    )
                    result_queue.put((i_species, "success"))
                except queue.Empty:
                    continue
                except Exception:
                    result_queue.put(("error", traceback.format_exc()))
                    break
        finally:
            result_queue.put(None) 

    # @staticmethod
    # def worker_process_blast_hits_new(
    #     seqsInfo,
    #     blastDir_list,
    #     Lengths,
    #     iSpecies,
    #     d_pickle,
    #     qDoubleBlast,
    #     v2_scores,
    #     q_allow_empty,
    #     result_queue,
    # ):
    #     with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    #         try:
    #             WaterfallMethod.process_blast_hits(
    #                 seqsInfo,
    #                 blastDir_list,
    #                 Lengths, 
    #                 iSpecies,
    #                 d_pickle,
    #                 qDoubleBlast,
    #                 v2_scores,
    #                 q_allow_empty,
    #             )
    #             result_queue.put((iSpecies, "success"))
    #         except Exception as e:
    #             result_queue.put((iSpecies, e))

    @staticmethod
    def _row_max(W):
        """Row maxima of a CSR matrix, and which rows have any entry."""
        has = np.diff(W.indptr) > 0
        m = np.full(W.shape[0], -np.inf)
        if W.nnz:
            m[has] = np.maximum.reduceat(W.data, W.indptr[:-1][has])
        return m, has

    @staticmethod
    def _entries_above(W, threshold_per_row):
        """0/1 CSR matrix of the entries of W greater than their row's threshold."""
        rows = np.repeat(np.arange(W.shape[0]), np.diff(W.indptr))
        keep = W.data > threshold_per_row[rows]
        return sparse.csr_matrix(
            (np.ones(int(keep.sum())), (rows[keep], W.indices[keep])), shape=W.shape
        )

    @staticmethod
    def get_bh_s(pairwise_scores_matrices, seqs_info, i_species, tol=1e-3):
        """
        Best hits of each gene of iSpecies in every species (within tol of the
        row's best score); within its own species, the hits within tol of its
        best hit in any other species.

        Works on CSR arrays: a Python loop over every row of every species
        pair costs species^2 x genes calls.
        """
        n_seqs_i = seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_species]]
        best_hit_for_sequence = -1 * np.ones(n_seqs_i)
        H = [None for j in range(seqs_info.n_species)]
        for j in range(seqs_info.n_species):
            if i_species == j:
                # identify orthologs then come back to paralogs
                continue
            W = pairwise_scores_matrices[j].tocsr()
            m, has = WaterfallMethod._row_max(W)
            update = has & (m > best_hit_for_sequence)
            best_hit_for_sequence[update] = m[update]
            H[j] = WaterfallMethod._entries_above(W, m - tol)
        # now look for paralogs
        W = pairwise_scores_matrices[i_species].tocsr()
        H[i_species] = WaterfallMethod._entries_above(W, best_hit_for_sequence - tol)
        return H

    @staticmethod
    def connect_cognates(seqs_info, i_species, d_pickle, v2_scores=False):
        # calculate RBH for species i
        b_hix = matrices.load_matrix_array("BH", seqs_info, i_species, d_pickle)
        b_hxi = matrices.load_matrix_array("BH", seqs_info, i_species, d_pickle, row=False)
        rb_hi = matrices.matrices_and_tr_s(
            b_hix, b_hxi
        )  # twice as much work as before (only did upper triangular before)
        del b_hix, b_hxi
        B = matrices.load_matrix_array("B", seqs_info, i_species, d_pickle)
        connect = WaterfallMethod.connect_all_better_than_an_ortholog_s(
            rb_hi, B, seqs_info, i_species, v2_scores
        )
        matrices.dump_matrix_array("connect", connect, i_species, d_pickle)

    @staticmethod
    def worker_connect_cognates_new(
        cmd_queue, 
        result_queue, 
        d_pickle, 
        v2_scores=True
    ):
        try:
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    if args is None:
                        break
                    WaterfallMethod.connect_cognates(
                        *args, d_pickle=d_pickle, v2_scores=v2_scores
                    )

                    i_species = args[1]
                    result_queue.put((i_species, "success"))
                except queue.Empty:
                    continue
                except Exception:
                    result_queue.put(("error", traceback.format_exc()))
                    break
        finally:
            result_queue.put(None) 

    # @staticmethod
    # def worker_connect_cognates_new(seqsInfo, iSpecies, d_pickle, result_queue, v2_scores=True):
    #     with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    #         try:
    #             WaterfallMethod.connect_cognates(
    #                 seqsInfo, iSpecies, d_pickle=d_pickle, v2_scores=v2_scores
    #             )
    #             result_queue.put((iSpecies, "success"))
    #         except Exception as e:
    #             result_queue.put((iSpecies, e))
        

    @staticmethod
    def write_graph_parallel(func, seqs_info, n_process, i_unassigned=None):
        graph_fn = file_handler.handler.get_graph_filename(i_unassigned)
        # Should use PTM?
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(graph_fn, "w") as graph_file:
                graph_file.write(
                    "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"
                    % (seqs_info.n_seqs, seqs_info.n_seqs)
                )
                graph_file.write("\n(mclmatrix\nbegin\n\n")
            # On a failure (or a worker killed, e.g. out of memory) parallel_map
            # stops the other workers and raises the child's error.
            # Largest species first (load balancing); the parts are still
            # joined in species order below.
            n_seqs = [seqs_info.n_seqs_per_species[sp] for sp in seqs_info.species_to_use]
            args = [
                (seqs_info, graph_fn, i_spec, file_handler.handler.get_pickle_dir())
                for i_spec in sorted(range(seqs_info.n_species), key=lambda i: -n_seqs[i])
            ]
            parallel_task_manager.parallel_map(func, args, n_process)
            with open(graph_fn, "ab") as graph_file:
                for i_sp in range(seqs_info.n_species):
                    part_fn = graph_fn + "_%d" % i_sp
                    with open(part_fn, "rb") as part:
                        shutil.copyfileobj(part, graph_file, 16 * 1024 * 1024)
                    os.remove(part_fn)
            # Cleanup
            matrices.delete_matrices("B", file_handler.handler.get_pickle_dir())
            matrices.delete_matrices("connect", file_handler.handler.get_pickle_dir())
        return graph_fn

    @staticmethod
    def get_most_distant_s(RBH, B, seqs_info, i_spec):
        # most distant RBB - as cut-off for connecting to other genes
        most_distant = numeric.transpose(
            np.ones(seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]]) * 1e9
        )
        # Best hit in another species: species-specific paralogues will now be connected - closer than any gene in any other species
        best_hit = numeric.transpose(
            np.zeros(seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]])
        )
        for k_spec in range(seqs_info.n_species):
            B[k_spec] = B[k_spec].tocsr()
            if i_spec == k_spec:
                continue
            best_hit = np.maximum(best_hit, matrices.sparse_max_row(B[k_spec]))
            I, J = RBH[k_spec].nonzero()
            if len(I) > 0:
                most_distant[I] = np.minimum(B[k_spec][I, J], most_distant[I])
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = most_distant > 1e8
        most_distant[I] = (
            best_hit[I] + 1e-6
        )  # to connect to one in it's own species it must be closer than other species. We can deal with hits outside the species later
        return most_distant

    @staticmethod
    def get_most_distant_s_estimate_from_relative_rbhs(RBH, B, seqs_info, i_spec, p=90):
        """
        Get the cut-off score for each sequence
        Args:
            RBH - array of 0-1 RBH lil matrices
            B - array of corresponding lil score matrices (Similarity ~ 1/Distance)
            p - the percentile of times the estimate will include the most distant RBH
        Returns:
            mostDistant - the cut-off for each sequence in this species
        Implementation:
                        Stage 1: conversion factors for rbh from iSpec to furthest rbh
            - Get the relative similarity scores for RBHs from iSpec to species j
              versus species k where the two RBHs exist
            - Find the p^th percentile of what the rbh would be in species j given
              what the rbh is in species k. This is the conversion factor for
              similarity score of an rbh in j vs k
            - Take the column min: I.e. the expected score for the furthest RBH
            Stage 2: For each rbh observed, estimate what this would imply would be the cut-off
             for the furthest rbh
             - Take the furthest of these estimates (could potentially use a percentile, but I
               think the point is that we should take diverged genes as a clue we should look
               further and then allow the trees to sort it out.
        """
        # 1. Reduce to column vectors
        q_debug_gather = False
        RBH = [rbh.tocsr() for rbh in RBH]
        B = [b.tocsr() for b in B]
        RBH_B = [
            sparse.csr_matrix(rbh.multiply(b))
            for i, (rbh, b) in enumerate(zip(RBH, B))
            if i != i_spec
        ]
        # create a vector of these scores
        nseqi = seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]]

        def rbh_scores(rbh_b):
            # The score in rows with exactly one stored entry, else 0.
            rbh_b.sort_indices()
            counts = np.diff(rbh_b.indptr)
            z = np.zeros(nseqi)
            one = counts == 1
            z[one] = rbh_b.data[rbh_b.indptr[:-1][one]]
            return z

        # nsp-1 x nseqi
        Z = np.matrix([rbh_scores(rbh_b) for rbh_b in RBH_B])  # RBH if it exists else zero
        nsp_m1 = Z.shape[0]
        # Zr = 1.0 / Z
        Z = np.asarray(Z, dtype=float)
        # rs[isp, jsp]: the (100-p)th percentile of Z[isp]/Z[jsp] over the genes
        # with an RBH in both species (1.0 on the diagonal and where there are
        # none). All jsp are done together for each isp instead of one
        # np.percentile call per pair (species^2 calls for each species);
        # _row_percentiles gives exactly the values np.percentile gives.
        valid = Z > 0
        rs = np.ones((nsp_m1, nsp_m1))
        with np.errstate(divide="ignore", invalid="ignore"):
            for isp in range(nsp_m1):
                both = valid & valid[isp]
                both[isp] = False
                rows = np.flatnonzero(both.any(axis=1))
                for start in range(0, len(rows), 256):   # bounded memory
                    block = rows[start:start + 256]
                    ratio = Z[isp, :] / Z[block]
                    rs[isp, block] = WaterfallMethod._row_percentiles(ratio, both[block], 100 - p)
        C = np.matrix(rs)  # Conversion matrix:
        # To convert from a hit in species j to one in species i multiply by Cij
        # To convert from a hit in species j to the furthest hit, need to take
        # the column-wise minimum
        C = np.amin(
            C, axis=0
        )  # conversion from a RBH to the estmate of what the most distant RBH should be
        most_distant = numeric.transpose(
            np.ones(seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]]) * 1e9
        )

        best_hit = numeric.transpose(
            np.zeros(seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]])
        )
        j_conversion = 0
        for k_spec in range(seqs_info.n_species):
            B[k_spec] = B[k_spec].tocsr()
            if i_spec == k_spec:
                continue
            species_factor_to_most_distant = C[0, j_conversion]
            if q_debug_gather:
                print(species_factor_to_most_distant)
            j_conversion += 1
            best_hit = np.maximum(best_hit, matrices.sparse_max_row(B[k_spec]))
            I, J = RBH[k_spec].nonzero()
            if len(I) > 0:
                if q_debug_gather:
                    print(
                        "%d updated"
                        % np.count_nonzero(
                            species_factor_to_most_distant * B[k_spec][I, J]
                            < most_distant[I]
                        )
                    )
                most_distant[I] = np.minimum(
                    species_factor_to_most_distant * B[k_spec][I, J], most_distant[I]
                )
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = most_distant > 1e8
        most_distant[I] = (
            best_hit[I] + 1e-6
        )  # to connect to one in its own species it must be closer than other species. We can deal with hits outside the species later
        return most_distant

    @staticmethod
    def _row_percentiles(values, valid, q):
        """
        For each row r: np.percentile(values[r][valid[r]], q) (linear method),
        for all rows at once. Every row must have at least one valid value.
        Follows numpy's computation step by step, so the results are
        bit-identical to calling np.percentile on each row.
        """
        counts = valid.sum(axis=1)
        x = np.where(valid, values, np.nan)
        x.sort(axis=1)                       # NaNs sort to the end
        qf = np.true_divide(q, 100)
        v = (counts - 1) * qf                 # virtual index
        prev = np.floor(v)
        above = v >= counts - 1               # at/after the last value: take the last value
        i_prev = np.where(above, counts - 1, prev).astype(np.intp)
        i_next = np.where(above, counts - 1, prev + 1).astype(np.intp)
        r = np.arange(len(counts))
        a = x[r, i_prev]
        b = x[r, i_next]
        gamma = v - np.where(above, -1.0, prev)
        diff_b_a = b - a
        out = a + diff_b_a * gamma
        hi = gamma >= 0.5
        out[hi] = (b - diff_b_a * (1 - gamma))[hi]
        return out

    @staticmethod
    def connect_all_better_than_cutoff_s(B, most_distant, seqs_info, i_spec):
        connect = []
        n_seqs_i = seqs_info.n_seqs_per_species[seqs_info.species_to_use[i_spec]]
        for j_spec in range(seqs_info.n_species):
            M = B[j_spec].tocsr()
            rows = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))
            keep = M.data >= most_distant[rows]
            if i_spec == j_spec:
                keep &= rows != M.indices
            mat = sparse.csr_matrix(
                (np.ones(int(keep.sum())), (rows[keep], M.indices[keep])),
                shape=(n_seqs_i, seqs_info.n_seqs_per_species[seqs_info.species_to_use[j_spec]]),
            )
            connect.append(mat)
        return connect

    @staticmethod
    def connect_all_better_than_an_ortholog_s(RBH, B, seqs_info, i_spec, v2_scores):
        if v2_scores:
            most_distant = WaterfallMethod.get_most_distant_s_estimate_from_relative_rbhs(
                RBH, B, seqs_info, i_spec
            )
        else:
            most_distant = WaterfallMethod.get_most_distant_s(RBH, B, seqs_info, i_spec)
        connect = WaterfallMethod.connect_all_better_than_cutoff_s(
            B, most_distant, seqs_info, i_spec
        )
        return connect
