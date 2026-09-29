import os
import shutil
import traceback
import numpy as np
from scipy import sparse
import warnings
import numpy.core.numeric as numeric
from scipy.optimize import curve_fit
import multiprocessing as mp
from ..utils import util, files, blast_file_processor, matrices, parallel_task_manager

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


class scnorm:
    @staticmethod
    def loglinear(x, a, b):
        return a * np.log10(x) + b

    @staticmethod
    def GetLengthArraysForMatrix(m, len_i, len_j):
        # Row-major order, as the previous lil-based version produced.
        m = m.tocsr()
        m.eliminate_zeros()
        m.sort_indices()
        I = np.repeat(np.arange(m.shape[0]), np.diff(m.indptr))
        J = m.indices
        scores = m.data
        Li = np.array(len_i[I])
        Lj = np.array(len_j[J])
        return Li, Lj, scores

    @staticmethod
    def GetTopPercentileOfScores(L, S, percentileToKeep):
        # Get the top x% of hits at each length
        L = np.asarray(L)
        S = np.asarray(S, dtype=float)
        nScores = len(S)
        # Sort by length, ties by original position (as sorted(zip(L, range(n))) did).
        indices = np.argsort(L, kind="stable")
        s_sorted = S[indices]
        l_sorted = L[indices]
        if nScores < 100:
            # then we can't split them into bins, return all for fitting
            return l_sorted, s_sorted
        nInBins = 1000 if nScores > 5000 else (200 if nScores > 1000 else 20)
        nBins, remainder = divmod(nScores, nInBins)
        # Full bins only (the remainder is not used), one percentile per bin.
        n = nBins * nInBins
        binned_scores = s_sorted[:n].reshape(nBins, nInBins)
        cutOffs = np.percentile(binned_scores, percentileToKeep, axis=1)
        keep = (binned_scores >= cutOffs[:, None]).ravel()
        return l_sorted[:n][keep], s_sorted[:n][keep]

    @staticmethod
    def CalculateFittingParameters(Lf, S):
        pars, covar = curve_fit(scnorm.loglinear, Lf, np.log10(S))
        return pars

    @staticmethod
    def NormaliseScoresByLogLengthProduct(b, Lq, Lh, params):
        rangeq = list(range(len(Lq)))
        rangeh = list(range(len(Lh)))
        li_vals = Lq ** (-params[0])
        lj_vals = Lh ** (-params[0])
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.csr_matrix(10 ** (-params[1]) * li_matrix * b * lj_matrix)


"""
WaterfallMethod
-------------------------------------------------------------------------------
"""


class WaterfallMethod:
    @staticmethod
    def NormaliseScores(B, Lengths, iSpecies, jSpecies):
        Li, Lj, scores = scnorm.GetLengthArraysForMatrix(
            B, Lengths[iSpecies], Lengths[jSpecies]
        )
        Lf = Li * Lj
        topLf, topScores = scnorm.GetTopPercentileOfScores(Lf, scores, 95)
        if len(topScores) > 1:
            fittingParameters = scnorm.CalculateFittingParameters(topLf, topScores)
            return scnorm.NormaliseScoresByLogLengthProduct(
                B, Lengths[iSpecies], Lengths[jSpecies], fittingParameters
            )
        elif iSpecies == jSpecies and B.nnz == 0:
            print(
                "WARNING: THIS IS UNCOMMON, there are zero hits when searching the genes in species %d against itself. Check the input proteome contains all the genes from that species and check the search program is working (default is diamond)."
                % iSpecies
            )
            return sparse.csr_matrix(B.get_shape())
        else:
            print(
                "WARNING: Too few hits between species %d and species %d to normalise the scores, these hits will be ignored"
                % (iSpecies, jSpecies)
            )
            return sparse.csr_matrix(B.get_shape())

    @staticmethod
    def NormalisedBitScore(B, Lengths, iSpecies, jSpecies):
        """
        Args:
            B - sparse matrix
        Returns
            B' - CSR matrix
        """
        if B.nnz == 0:
            return sparse.csr_matrix(B)
        Lq = Lengths[iSpecies]
        Lh = Lengths[jSpecies]
        rangeq = list(range(len(Lq)))
        rangeh = list(range(len(Lh)))
        li_vals = Lq ** (-0.5)
        lj_vals = Lh ** (-0.5)
        li_matrix = sparse.csr_matrix((li_vals, (rangeq, rangeq)))
        lj_matrix = sparse.csr_matrix((lj_vals, (rangeh, rangeh)))
        return sparse.csr_matrix(li_matrix * B * lj_matrix)

    @staticmethod
    def ProcessBlastHits(
        seqsInfo,
        blastDir_list,
        Lengths,
        iSpecies,
        d_pickle,
        qDoubleBlast,
        v2_scores,
        q_allow_empty=False,
    ):
        """
        iLimitNewSpecies: int - Only process fasta files with OrthoFidner ID >= iLimitNewSpecies
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # process up to the best hits for each species
            Bi = []
            for jSpecies in range(seqsInfo.nSpecies):
                Bij = blast_file_processor.GetBLAST6Scores(
                    seqsInfo,
                    blastDir_list,
                    seqsInfo.speciesToUse[iSpecies],
                    seqsInfo.speciesToUse[jSpecies],
                    qDoubleBlast=qDoubleBlast,
                    q_allow_empty=q_allow_empty,
                    fmt="csr",
                )
                if v2_scores:
                    Bij = WaterfallMethod.NormalisedBitScore(
                        Bij, Lengths, iSpecies, jSpecies
                    )
                else:
                    Bij = WaterfallMethod.NormaliseScores(
                        Bij, Lengths, iSpecies, jSpecies
                    )
                Bi.append(Bij)
            matrices.DumpMatrixArray("B", Bi, iSpecies, d_pickle)
            BH = WaterfallMethod.GetBH_s(Bi, seqsInfo, iSpecies)
            matrices.DumpMatrixArray("BH", BH, iSpecies, d_pickle)
            # util.PrintTime("Initial processing of species %d complete" % iSpecies)


    @staticmethod
    def Worker_ProcessBlastHits(
        seqsInfo,
        blastDir_list,
        Lengths,
        cmd_queue,
        d_pickle,
        qDoubleBlast,
        v2_scores,
        q_allow_empty,
    ):
        while True:
            try:
                iSpecies = cmd_queue.get(True, 1)
                WaterfallMethod.ProcessBlastHits(
                    seqsInfo,
                    blastDir_list,
                    Lengths,
                    iSpecies,
                    d_pickle=d_pickle,
                    qDoubleBlast=qDoubleBlast,
                    v2_scores=v2_scores,
                    q_allow_empty=q_allow_empty,
                )
            except queue.Empty:
                return
            except Exception:
                i = seqsInfo.speciesToUse[iSpecies]
                print("ERROR: Error processing files Blast%d_*" % i)
                raise
    @staticmethod
    def Worker_ProcessBlastHits_New(
        seqsInfo,
        blastDir_list,
        Lengths,
        cmd_queue,
        d_pickle,
        qDoubleBlast,
        v2_scores,
        q_allow_empty,
        result_queue,
    ):  
        try:
            while True:
                try:
                    iSpecies = cmd_queue.get(True, 1)
                    if iSpecies is None:
                        break
                    WaterfallMethod.ProcessBlastHits(
                        seqsInfo,
                        blastDir_list,
                        Lengths,
                        iSpecies,
                        d_pickle=d_pickle,
                        qDoubleBlast=qDoubleBlast,
                        v2_scores=v2_scores,
                        q_allow_empty=q_allow_empty,
                    )
                    result_queue.put((iSpecies, "success"))
                except queue.Empty:
                    continue
                except Exception:
                    result_queue.put(("error", traceback.format_exc()))
                    break
        finally:
            result_queue.put(None) 

    # @staticmethod
    # def Worker_ProcessBlastHits_New(
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
    #             WaterfallMethod.ProcessBlastHits(
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
    def _RowMax(W):
        """Row maxima of a CSR matrix, and which rows have any entry."""
        has = np.diff(W.indptr) > 0
        m = np.full(W.shape[0], -np.inf)
        if W.nnz:
            m[has] = np.maximum.reduceat(W.data, W.indptr[:-1][has])
        return m, has

    @staticmethod
    def _EntriesAbove(W, threshold_per_row):
        """0/1 CSR matrix of the entries of W greater than their row's threshold."""
        rows = np.repeat(np.arange(W.shape[0]), np.diff(W.indptr))
        keep = W.data > threshold_per_row[rows]
        return sparse.csr_matrix(
            (np.ones(int(keep.sum())), (rows[keep], W.indices[keep])), shape=W.shape
        )

    @staticmethod
    def GetBH_s(pairwiseScoresMatrices, seqsInfo, iSpecies, tol=1e-3):
        """
        Best hits of each gene of iSpecies in every species (within tol of the
        row's best score); within its own species, the hits within tol of its
        best hit in any other species.

        Works on CSR arrays: a Python loop over every row of every species
        pair costs species^2 x genes calls.
        """
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpecies]]
        bestHitForSequence = -1 * np.ones(nSeqs_i)
        H = [None for j in range(seqsInfo.nSpecies)]
        for j in range(seqsInfo.nSpecies):
            if iSpecies == j:
                # identify orthologs then come back to paralogs
                continue
            W = pairwiseScoresMatrices[j].tocsr()
            m, has = WaterfallMethod._RowMax(W)
            update = has & (m > bestHitForSequence)
            bestHitForSequence[update] = m[update]
            H[j] = WaterfallMethod._EntriesAbove(W, m - tol)
        # now look for paralogs
        W = pairwiseScoresMatrices[iSpecies].tocsr()
        H[iSpecies] = WaterfallMethod._EntriesAbove(W, bestHitForSequence - tol)
        return H

    @staticmethod
    def ConnectCognates(seqsInfo, iSpecies, d_pickle, v2_scores=False):
        # calculate RBH for species i
        BHix = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle)
        BHxi = matrices.LoadMatrixArray("BH", seqsInfo, iSpecies, d_pickle, row=False)
        RBHi = matrices.MatricesAndTr_s(
            BHix, BHxi
        )  # twice as much work as before (only did upper triangular before)
        del BHix, BHxi
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpecies, d_pickle)
        connect = WaterfallMethod.ConnectAllBetterThanAnOrtholog_s(
            RBHi, B, seqsInfo, iSpecies, v2_scores
        )
        matrices.DumpMatrixArray("connect", connect, iSpecies, d_pickle)

    @staticmethod
    def Worker_ConnectCognates(cmd_queue, d_pickle, v2_scores=True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            while True:
                try:
                    args = cmd_queue.get(True, 1)
                    WaterfallMethod.ConnectCognates(
                        *args, d_pickle=d_pickle, v2_scores=v2_scores
                    )
                except queue.Empty:
                    return


    @staticmethod
    def Worker_ConnectCognates_New(
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
                    WaterfallMethod.ConnectCognates(
                        *args, d_pickle=d_pickle, v2_scores=v2_scores
                    )

                    iSpecies = args[1]
                    result_queue.put((iSpecies, "success"))
                except queue.Empty:
                    continue
                except Exception:
                    result_queue.put(("error", traceback.format_exc()))
                    break
        finally:
            result_queue.put(None) 

    # @staticmethod
    # def Worker_ConnectCognates_New(seqsInfo, iSpecies, d_pickle, result_queue, v2_scores=True):
    #     with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    #         try:
    #             WaterfallMethod.ConnectCognates(
    #                 seqsInfo, iSpecies, d_pickle=d_pickle, v2_scores=v2_scores
    #             )
    #             result_queue.put((iSpecies, "success"))
    #         except Exception as e:
    #             result_queue.put((iSpecies, e))
        

    @staticmethod
    def WriteGraphParallel(func, seqsInfo, nProcess, i_unassigned=None):
        graphFN = files.FileHandler.GetGraphFilename(i_unassigned)
        # Should use PTM?
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(graphFN, "w") as graphFile:
                graphFile.write(
                    "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"
                    % (seqsInfo.nSeqs, seqsInfo.nSeqs)
                )
                graphFile.write("\n(mclmatrix\nbegin\n\n")
            # On a failure (or a worker killed, e.g. out of memory) ParallelMap
            # stops the other workers and raises the child's error.
            args = [
                (seqsInfo, graphFN, iSpec, files.FileHandler.GetPickleDir())
                for iSpec in range(seqsInfo.nSpecies)
            ]
            parallel_task_manager.ParallelMap(func, args, nProcess)
            with open(graphFN, "ab") as graphFile:
                for iSp in range(seqsInfo.nSpecies):
                    part_fn = graphFN + "_%d" % iSp
                    with open(part_fn, "rb") as part:
                        shutil.copyfileobj(part, graphFile, 16 * 1024 * 1024)
                    os.remove(part_fn)
            # Cleanup
            matrices.DeleteMatrices("B", files.FileHandler.GetPickleDir())
            matrices.DeleteMatrices("connect", files.FileHandler.GetPickleDir())
        return graphFN

    @staticmethod
    def GetMostDistant_s(RBH, B, seqsInfo, iSpec):
        # most distant RBB - as cut-off for connecting to other genes
        mostDistant = numeric.transpose(
            np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]) * 1e9
        )
        # Best hit in another species: species-specific paralogues will now be connected - closer than any gene in any other species
        bestHit = numeric.transpose(
            np.zeros(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])
        )
        for kSpec in range(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            bestHit = np.maximum(bestHit, matrices.sparse_max_row(B[kSpec]))
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                mostDistant[I] = np.minimum(B[kSpec][I, J], mostDistant[I])
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = mostDistant > 1e8
        mostDistant[I] = (
            bestHit[I] + 1e-6
        )  # to connect to one in it's own species it must be closer than other species. We can deal with hits outside the species later
        return mostDistant

    @staticmethod
    def GetMostDistant_s_estimate_from_relative_rbhs(RBH, B, seqsInfo, iSpec, p=90):
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
            if i != iSpec
        ]
        # create a vector of these scores
        nseqi = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]

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
        rs = []
        for isp in range(nsp_m1):
            rs.append([])
            for jsp in range(nsp_m1):
                if isp == jsp:
                    rs[-1].append(1.0)
                    continue

                num = Z[isp, :]
                den = Z[jsp, :]
                mask = (num > 0) & (den > 0) 
                if not np.any(mask):
                    rs[-1].append(1.0)
                else:
                    ratio = num[mask] / den[mask]
                    r = np.percentile(ratio, 100 - p)
                    rs[-1].append(r)

                # ratios = np.multiply(Z[isp, :], Zr[jsp, :])
                # i_nonzeros = np.where(
                #     np.logical_and(np.isfinite(ratios), ratios > 0)
                # )  # only those which a score is availabe for both
                # if i_nonzeros[0].size == 0:
                #     rs[-1].append(1.0)  # no conversion between this pair
                # else:
                #     ratios = ratios[i_nonzeros]
                #     r = np.percentile(np.asarray(ratios), 100 - p)
                #     rs[-1].append(r)

                # calculate the n vectors and take the min
        C = np.matrix(rs)  # Conversion matrix:
        # To convert from a hit in species j to one in species i multiply by Cij
        # To convert from a hit in species j to the furthest hit, need to take
        # the column-wise minimum
        C = np.amin(
            C, axis=0
        )  # conversion from a RBH to the estmate of what the most distant RBH should be
        mostDistant = numeric.transpose(
            np.ones(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]) * 1e9
        )

        bestHit = numeric.transpose(
            np.zeros(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]])
        )
        j_conversion = 0
        for kSpec in range(seqsInfo.nSpecies):
            B[kSpec] = B[kSpec].tocsr()
            if iSpec == kSpec:
                continue
            species_factor_to_most_distant = C[0, j_conversion]
            if q_debug_gather:
                print(species_factor_to_most_distant)
            j_conversion += 1
            bestHit = np.maximum(bestHit, matrices.sparse_max_row(B[kSpec]))
            I, J = RBH[kSpec].nonzero()
            if len(I) > 0:
                if q_debug_gather:
                    print(
                        "%d updated"
                        % np.count_nonzero(
                            species_factor_to_most_distant * B[kSpec][I, J]
                            < mostDistant[I]
                        )
                    )
                mostDistant[I] = np.minimum(
                    species_factor_to_most_distant * B[kSpec][I, J], mostDistant[I]
                )
        # anything that doesn't have an RBB, set to distance to the closest gene in another species. I.e. it will hit just that gene and all genes closer to it in the its own species
        I = mostDistant > 1e8
        mostDistant[I] = (
            bestHit[I] + 1e-6
        )  # to connect to one in its own species it must be closer than other species. We can deal with hits outside the species later
        return mostDistant

    @staticmethod
    def ConnectAllBetterThanCutoff_s(B, mostDistant, seqsInfo, iSpec):
        connect = []
        nSeqs_i = seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]
        for jSpec in range(seqsInfo.nSpecies):
            M = B[jSpec].tocsr()
            rows = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))
            keep = M.data >= mostDistant[rows]
            if iSpec == jSpec:
                keep &= rows != M.indices
            mat = sparse.csr_matrix(
                (np.ones(int(keep.sum())), (rows[keep], M.indices[keep])),
                shape=(nSeqs_i, seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[jSpec]]),
            )
            connect.append(mat)
        return connect

    @staticmethod
    def ConnectAllBetterThanAnOrtholog_s(RBH, B, seqsInfo, iSpec, v2_scores):
        if v2_scores:
            mostDistant = WaterfallMethod.GetMostDistant_s_estimate_from_relative_rbhs(
                RBH, B, seqsInfo, iSpec
            )
        else:
            mostDistant = WaterfallMethod.GetMostDistant_s(RBH, B, seqsInfo, iSpec)
        connect = WaterfallMethod.ConnectAllBetterThanCutoff_s(
            B, mostDistant, seqsInfo, iSpec
        )
        return connect
