from __future__ import absolute_import

import os
import numpy as np
import subprocess
from scipy import sparse
import random
import time
import warnings
from collections import Counter, deque
import numpy.core.numeric as numeric
from scipy.optimize import curve_fit
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import queue
except ImportError:
    import Queue as queue

from ..tools import mcl, trees_msa, waterfall
from . import orthogroups_set
from ..utils import util, files, blast_file_processor, matrices, parallel_task_manager

try:
    from rich import print, progress
except ImportError:
    ...


def WriteGraph_perSpecies(args):
    seqsInfo, graphFN, iSpec, d_pickle = args
    # calculate the 2-way connections for one query species
    with open(graphFN + "_%d" % iSpec, "w") as graphFile:
        connect2 = []
        for jSpec in range(seqsInfo.nSpecies):
            m1 = matrices.LoadMatrix("connect", iSpec, jSpec, d_pickle)
            m2tr = numeric.transpose(
                matrices.LoadMatrix("connect", jSpec, iSpec, d_pickle)
            )
            connect2.append(m1 + m2tr)
            del m1, m2tr
        B = matrices.LoadMatrixArray("B", seqsInfo, iSpec, d_pickle)
        B_connect = matrices.MatricesAnd_s(connect2, B)
        del B, connect2

        W = [b.sorted_indices().tolil() for b in B_connect]
        del B_connect
        for query in range(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in range(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j, value in zip(row.rows[0], row.data[0]):
                    graphFile.write("%d:%.3f " % (j + jOffset, value))
            graphFile.write("$\n")
        if iSpec == (seqsInfo.nSpecies - 1):
            graphFile.write(")\n")
        util.PrintTime("Written final scores for species %d to graph file" % iSpec)


def WriteGraph_perSpecies_homology(args):
    seqsInfo, graphFN, iSpec, d_pickle = args
    # calculate the 2-way connections for one query species
    # W = [matrices.LoadMatrix("B", iSpec, jSpec, d_pickle).tolil() for jSpec in range(seqsInfo.nSpecies)]
    W = []
    for jSpec in range(seqsInfo.nSpecies):
        w1 = matrices.LoadMatrix("B", iSpec, jSpec, d_pickle)
        matrices.DumpMatrix("H", (w1 > 0).tolil(), iSpec, jSpec, d_pickle)
        w2tr = numeric.transpose(matrices.LoadMatrix("B", jSpec, iSpec, d_pickle))
        W.append((w1 + w2tr > 0).tolil())  # symmetrise
    # matrices.DumpMatrixArray("H", W, iSpec, d_pickle)
    with open(graphFN + "_%d" % iSpec, "w") as graphFile:
        for query in range(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSpec]]):
            offset = seqsInfo.seqStartingIndices[iSpec]
            graphFile.write("%d    " % (offset + query))
            for jSpec in range(seqsInfo.nSpecies):
                row = W[jSpec].getrowview(query)
                jOffset = seqsInfo.seqStartingIndices[jSpec]
                for j in row.rows[0]:
                    graphFile.write("%d:%.3f " % (j + jOffset, 1.0))
            graphFile.write("$\n")
        if iSpec == (seqsInfo.nSpecies - 1):
            graphFile.write(")\n")
        util.PrintTime("Written final scores for species %d to graph file" % iSpec)


def GetSequenceLengths(seqsInfo):
    sequenceLengths = []
    for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
        sequenceLengths.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
        fastaFilename = files.FileHandler.GetSpeciesFastaFN(iFasta)
        currentSequenceLength = 0
        iCurrentSequence = -1
        qFirstLine = True
        with open(fastaFilename) as infile:
            for row in infile:
                if len(row) > 1 and row[0] == ">":
                    if qFirstLine:
                        qFirstLine = False
                    else:
                        sequenceLengths[iSpecies][
                            iCurrentSequence
                        ] = currentSequenceLength
                        currentSequenceLength = 0
                    _, iCurrentSequence = util.GetIDPairFromString(row[1:])
                else:
                    currentSequenceLength += len(row.rstrip())
        sequenceLengths[iSpecies][iCurrentSequence] = currentSequenceLength
    return sequenceLengths


def DoOrthogroups(
    options,
    speciesInfoObj,
    seqsInfo,
    speciesNamesDict,
    speciesXML=None,
    i_unassigned=None,
):

    try:
        width = os.get_terminal_size().columns
    except OSError as e:
        width = 80



    # Run Algorithm, cluster and output cluster files with original accessions
    q_unassigned = i_unassigned is not None
    util.PrintUnderline(
        "Running OrthoFinder algorithm"
        + (" for clade-specific genes" if q_unassigned else "")
    )
    # it's important to free up the memory from python used for processing the genomes
    # before launching MCL because both use sizeable amounts of memory. The only
    # way I can find to do this is to launch the memory intensive python code
    # as separate process that exits before MCL is launched.

    Lengths = GetSequenceLengths(seqsInfo)  # Alternatively, self-self bit scores, but it amounts to the same thing

    # Process BLAST hits
    util.PrintTime("Initial processing of each species")
    cmd_queue = mp.Queue()
    blastDir_list = files.FileHandler.GetBlastResultsDir()
    if q_unassigned:
        blastDir_list = blastDir_list[:1]  # only use latet directory with unassigned gene searches
    for iSpeciesJob in range(seqsInfo.nSpecies):  # The i-th job, not the OrthoFinder species ID
        cmd_queue.put(iSpeciesJob)

    files.FileHandler.GetPickleDir()  # create the pickle directory before the parallel processing to prevent a race condition
    # Should use PTM?
    # args_list = [(seqsInfo, blastDir_list, Lengths, cmd_queue, files.FileHandler.GetPickleDir(), options.qDoubleBlast, options.v2_scores, q_unassigned)
    #              for i_ in range(options.nProcessAlg)]
    # parallel_task_manager.RunParallelMethods(WaterfallMethod.Worker_ProcessBlastHits, args_list, options.nProcessAlg)

    # runningProcesses = [
    #     mp.Process(
    #         target=waterfall.WaterfallMethod.Worker_ProcessBlastHits,
    #         args=(
    #             seqsInfo,
    #             blastDir_list,
    #             Lengths,
    #             cmd_queue,
    #             files.FileHandler.GetPickleDir(),
    #             options.qDoubleBlast,
    #             options.v2_scores,
    #             q_unassigned,
    #         ),
    #     )
    #     for i_ in range(options.nProcessAlg)
    # ]

    # for proc in runningProcesses:
    #     proc.start()
    # parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)

## ----------------------------------------------------------------------------

    gathering_progress = progress.Progress(
        progress.TextColumn("[progress.description]{task.description}"),
        progress.BarColumn(bar_width=width // 2),
        progress.SpinnerColumn(),
        progress.MofNCompleteColumn(),
        progress.TimeElapsedColumn(),
        transient=False,
        # progress.TextColumn("{task.completed}/{task.total}")
    )
    task = gathering_progress.add_task(
        "[yellow]Processing...", total=seqsInfo.nSpecies
    )
    update_cycle = 1

## ---------------------------------------------------------
    gathering_progress.start()
    with ProcessPoolExecutor(max_workers=options.nProcessAlg) as executor:
        # Submit the worker task n_processes times
        futures = {
            executor.submit(
            waterfall.WaterfallMethod.ProcessBlastHits,
            seqsInfo,
            blastDir_list,
            Lengths, 
            iSpeciesJob,
            files.FileHandler.GetPickleDir(),
            options.qDoubleBlast,
            options.v2_scores,
            q_unassigned,
            ): iSpeciesJob
            for iSpeciesJob in range(seqsInfo.nSpecies)
        }

        for i, future in enumerate(as_completed(futures)):
            iSpeciesJob = futures[future]
            try:
                future.result() 
            except Exception as e:
                for f in futures:
                    f.cancel()
                print("ERROR: Error processing files Blast%d_*" % seqsInfo.speciesToUse[iSpeciesJob])
                util.Fail()  
                break
            finally:
                if (i + 1) % update_cycle == 0:
                    gathering_progress.update(task, advance=update_cycle)
    gathering_progress.stop()
## --------------------------------------------------------------------

    if options.gathering_version < (3, 0):
        cmd_queue = mp.Queue()
        for iSpecies in range(seqsInfo.nSpecies):
            cmd_queue.put((seqsInfo, iSpecies))
        # args_list = [(cmd_queue, files.FileHandler.GetPickleDir(), options.v2_scores) for i_ in range(options.nProcessAlg)]
        # parallel_task_manager.RunParallelMethods(WaterfallMethod.Worker_ConnectCognates, args_list, options.nProcessAlg)
        runningProcesses = [
            mp.Process(
                target=waterfall.WaterfallMethod.Worker_ConnectCognates,
                args=(cmd_queue, files.FileHandler.GetPickleDir(), options.v2_scores),
            )
            for i_ in range(options.nProcessAlg)
        ]
        for proc in runningProcesses:
            proc.start()
        parallel_task_manager.ManageQueue(runningProcesses, cmd_queue)

        util.PrintTime("Connected putative homologues")
        graphFilename = waterfall.WaterfallMethod.WriteGraphParallel(
            WriteGraph_perSpecies, seqsInfo, options.nProcessAlg, i_unassigned
        )

        # 5b. MCL
        clustersFilename, clustersFilename_pairs = (
            files.FileHandler.CreateUnusedClustersFN(
                "_I%0.1f" % options.mclInflation, i_unassigned
            )
        )
        mcl.MCL.RunMCL(
            graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation
        )
        # If processing unassigned, then ignore all 'unclustered' genes - they will include any genes not included in this search
        mcl.ConvertSingleIDsToIDPair(
            seqsInfo, clustersFilename, clustersFilename_pairs, q_unassigned
        )
    elif options.gathering_version == (3, 2):
        graphFilename = waterfall.WaterfallMethod.WriteGraphParallel(
            WriteGraph_perSpecies_homology, seqsInfo, options.nProcessAlg, i_unassigned
        )
        clustersFilename, clustersFilename_pairs = (
            files.FileHandler.CreateUnusedClustersFN(
                "_I%0.1f" % options.mclInflation, i_unassigned
            )
        )
        mcl.MCL.RunMCL(
            graphFilename, clustersFilename, options.nProcessAlg, options.mclInflation
        )
        mcl.ConvertSingleIDsToIDPair(
            seqsInfo, clustersFilename, clustersFilename_pairs, q_unassigned
        )
    if not q_unassigned:
        post_clustering_orthogroups(
            clustersFilename_pairs,
            speciesInfoObj,
            seqsInfo,
            speciesNamesDict,
            options,
            speciesXML,
        )
    return clustersFilename_pairs


def post_clustering_orthogroups(
    clustersFilename_pairs,
    speciesInfoObj,
    seqsInfo,
    speciesNamesDict,
    options,
    speciesXML,
    q_incremental=False,
):
    """
    Write OGs & statistics to results files, write Fasta files.
    Args:
        q_incremental - These are not the final orthogroups, don't write results
    """
    ogs = mcl.GetPredictedOGs(clustersFilename_pairs)
    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    util.PrintUnderline("Writing orthogroups to file")
    idsDict = mcl.MCL.WriteOrthogroupFiles(
        ogs,
        [files.FileHandler.GetSequenceIDsFN()],
        resultsBaseFilename,
        clustersFilename_pairs,
    )
    if not q_incremental:
        mcl.MCL.CreateOrthogroupTable(
            ogs,
            idsDict,
            speciesNamesDict,
            speciesInfoObj.speciesToUse,
            resultsBaseFilename,
        )

    # Write Orthogroup FASTA files
    ogSet = orthogroups_set.OrthoGroupsSet(
        files.FileHandler.GetWorkingDirectory1_Read(),
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        idExtractor=util.FirstWordExtractor,
    )
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(
        files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse
    )
    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs):
        os.mkdir(d_seqs)
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idsDict, False)

    if not q_incremental:
        # stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
        if options.speciesXMLInfoFN:
            mcl.MCL.WriteOrthoXML(
                speciesXML,
                ogs,
                seqsInfo.nSeqsPerSpecies,
                idsDict,
                resultsBaseFilename + ".orthoxml",
                speciesInfoObj.speciesToUse,
            )
        print("")
        util.PrintTime("Done orthogroups")
        files.FileHandler.LogOGs()
