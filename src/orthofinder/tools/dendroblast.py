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

from . import stag, tree
from ..utils import util, files, parallel_task_manager, program_caller
from ..utils import blast_file_processor as BlastFileProcessor


class DendroBLASTTrees(object):
    def __init__(self, ogSet, nProcesses_alg, nProcess_std, qDoubleBlast):
        self.ogSet = ogSet
        self.nProcesses = nProcesses_alg
        self.nProcess_std = nProcess_std
        self.qDoubleBlast = qDoubleBlast
        # Check files exist
    
    def TreeFilename_IDs(self, iog):
        return files.FileHandler.GetOGsTreeFN(iog)
        
    def GetOGMatrices_FullParallel(self):
        """
        read the blast files as well, remove need for intermediate pickle and unpickle
        ogMatrices contains matrix M for each OG where:
            Mij = 0.5*max(Bij, Bmin_i)/Bmax_i
        """
        with warnings.catch_warnings():         
            warnings.simplefilter("ignore")
            ogs_all = self.ogSet.OGsAll()
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
            ogsPerSpecies = [[[(g, i) for i, g in enumerate(og) if g.iSp == iSp] for iSp in self.ogSet.seqsInfo.speciesToUse] for og in ogs]
            nGenes = [len(og) for og in ogs]
            nSeqs = self.ogSet.seqsInfo.nSeqsPerSpecies
            ogMatrices = [[mp.Array('d', n, lock=False) for _ in range(n)] for n in nGenes]
            blastDir_list = files.FileHandler.GetBlastResultsDir()
            cmd_queue = mp.Queue()
            for iiSp, sp1 in enumerate(self.ogSet.seqsInfo.speciesToUse):
                cmd_queue.put((iiSp, sp1, nSeqs[sp1]))
            worker_status_queue = mp.Queue()
            # Should use PTM?
            runningProcesses = [mp.Process(target=Worker_OGMatrices_ReadBLASTAndUpdateDistances,
                                           args=(cmd_queue, worker_status_queue, iWorker, ogMatrices, nGenes,
                                                 self.ogSet.seqsInfo, blastDir_list, ogsPerSpecies, self.qDoubleBlast))
                                for iWorker in range(self.nProcesses)]
            for proc in runningProcesses:
                proc.start()
            rota = [None for iWorker in range(self.nProcesses)]
            unfinished = []
            while True:
                # get process alive/dead
                time.sleep(1)
                alive = [proc.is_alive() for proc in runningProcesses]
                # read latest updates from queue, update rota
                try:
                    while True:
                        status, iWorker, iTask = worker_status_queue.get(True, 0.1)
                        if status == "start":
                            rota[iWorker] = iTask
                        elif status == "finish":
                            rota[iWorker] = None
                        elif status == "empty":
                            rota[iWorker] = "empty"
                except queue.Empty:
                    pass
                # if worker is dead but didn't finish task, issue warning
                for al, r in zip(alive, rota):
                    if (not al) and (r != "empty"):
                        text = GetRAMErrorText()
                        files.FileHandler.LogFailAndExit(text)
                        unfinished.append(r)
                if not any(alive):
                    break
                
            if len(unfinished) != 0:
                files.FileHandler.LogFailAndExit()
#                print("WARNING: Computer ran out of RAM and killed OrthoFinder processes")
#                print("OrthoFinder will attempt to run these processes once more. If it is")
#                print("unsuccessful again then it will have to exit. Consider using")
#                print("the option '-a 1' or running on a machine with more RAM")
            #ogMatrices = [np.matrix(m) for m in ogMatrices]
            return ogs, ogMatrices      
                   
    def CompleteOGMatrices(self, ogs, ogMatrices):
        newMatrices = []
        for iog, (og, m) in enumerate(zip(ogs, ogMatrices)):
            # dendroblast scores
            n = m.shape[0]
            m2 = np.zeros(m.shape)
            max_og = -9e99
            for i in range(n):
                for j in range(i):
                    m2[i, j] = -np.log(m[i,j] + m[j,i])  
                    m2[j, i] = m2[i, j]  
                    max_og = max(max_og, m2[i,j])
            newMatrices.append(m2)
        return newMatrices
        
    def CompleteAndWriteOGMatrices(self, ogs, ogMatrices):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        for iog, (og, m) in enumerate(zip(ogs, ogMatrices)):
            # dendroblast scores
            n = len(m)
            max_og = -9e99
            # Careful not to over-write a value and then attempt to try to use the old value
            for i in range(n):
                for j in range(i):
                    m[i][j] = -np.log(m[i][j] + m[j][i])  
                    m[j][i] = m[i][j]  
                    max_og = max(max_og, m[i][j])
            self.WritePhylipMatrix(m, [g.ToString() for g in og], files.FileHandler.GetOGsDistMatFN(iog), max_og)
        return ogMatrices
    
    @staticmethod
    def WritePhylipMatrix(m, names, outFN, max_og):
        """
        m - list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        max_og = 1.1*max_og
        sliver = 1e-6
        with open(outFN, 'w') as outfile:
            n = len(m)
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(names[i] + " ")
                # values could be -inf, these are the most distantly related so replace with max_og
                V = [0. + (0. if i==j else m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")
    
    def SpeciesTreeDistances(self, ogs, ogMatrices, method = 0):
        """
        ogMatrices - each matrix is a list of mp.Array  (so that each represents an nSeq x nSeq matrix
        """
        spPairs = list(itertools.combinations(self.ogSet.seqsInfo.speciesToUse, 2))
        D = [[] for _ in spPairs]
        if method == 0:
            """ closest distance for each species pair in each orthogroup"""
            for og, m in zip(ogs, ogMatrices):
                spDict = defaultdict(list)
                for i, g in enumerate(og):
                    spDict[g.iSp].append(i)
                for (sp1, sp2), d_list in zip(spPairs, D):
                    distances = [m[i][j] for i in spDict[sp1] for j in spDict[sp2]]
                    if len(distances) > 0: d_list.append(min(distances))
#                    d_list.append(min(distances) if len(distances) > 0 else None)
        return D, spPairs
    
    def PrepareSpeciesTreeCommand(self, D, spPairs, qPutInWorkingDir=False):
        n = len(self.ogSet.seqsInfo.speciesToUse)
        M = np.zeros((n, n))
        for (sp1, sp2), d in zip(spPairs, D):
            sp1 = self.ogSet.seqsInfo.speciesToUse.index(sp1)
            sp2 = self.ogSet.seqsInfo.speciesToUse.index(sp2)
            x = np.median(d)
            M[sp1, sp2] = x
            M[sp2, sp1] = x
        speciesMatrixFN = files.FileHandler.GetSpeciesTreeMatrixFN(qPutInWorkingDir)  
        sliver = 1e-6
        with open(speciesMatrixFN, 'w') as outfile:
            outfile.write("%d\n" % n)
            for i in range(n):
                outfile.write(str(self.ogSet.seqsInfo.speciesToUse[i]) + " ")
                V = [(0. + M[i,j]) for j in range(n)]  # hack to avoid printing out "-0"
                V = [sliver if 0 < v < sliver else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
                values = " ".join(["%.6f" % v for v in V])   
                outfile.write(values + "\n")       
        treeFN = files.FileHandler.GetSpeciesTreeUnrootedFN()
        cmd = " ".join(["fastme", "-i", speciesMatrixFN, "-o", treeFN, "-N", "-w", "O"] + (["-s"] if n < 1000 else []))
        return cmd, treeFN
                
    def PrepareGeneTreeCommand(self):
        cmds = []
        ogs = self.ogSet.OGsAll()
        for iog in range(len(ogs)):
            nTaxa = len(ogs[iog])
            if nTaxa < 4:
                continue
            cmds.append([" ".join(["fastme", "-i", files.FileHandler.GetOGsDistMatFN(iog), "-o",
                                   files.FileHandler.GetOGsTreeFN(iog), "-N", "-w", "O"] + (["-s"] if nTaxa < 1000 else []))])
        return cmds

    @staticmethod    
    def EnoughOGsForSTAG(ogs, speciesToUse):
        nSp = len(speciesToUse)
        nSp_perOG = [len(set([g.iSp for g in og])) for og in ogs]
        return (nSp_perOG.count(nSp) >= 100)
    
    def RunAnalysis(self, qSpeciesTree=True):
        """
        Args:
            qSpeciesTree - Bool: infer the species tree
        """
        util.PrintUnderline("Calculating gene distances")
        ogs, ogMatrices_partial = self.GetOGMatrices_FullParallel()
        ogMatrices = self.CompleteAndWriteOGMatrices(ogs, ogMatrices_partial)
        del ogMatrices_partial
        util.PrintTime("Done")
        cmds_trees = self.PrepareGeneTreeCommand()
        qLessThanFourSpecies = len(self.ogSet.seqsInfo.speciesToUse) < 4
        if not qSpeciesTree:
            qSTAG = False
        elif qLessThanFourSpecies:
            qSTAG = False
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            WriteSpeciesTreeIDs_TwoThree(self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        else:
            qSTAG = self.EnoughOGsForSTAG(ogs, self.ogSet.seqsInfo.speciesToUse)
            if not qSTAG:
                print("Using fallback species tree inference method")
                D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
                cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs)
                cmds_trees = [[cmd_spTree]] + cmds_trees
        del ogMatrices
        util.PrintUnderline("Inferring gene and species trees" if qSpeciesTree else "Inferring gene trees")
        program_caller.RunParallelCommands(self.nProcess_std, cmds_trees, qListOfList=True)
        if qSTAG:
            # Trees must have been completed
            print("")
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            stag.Run_ForOrthoFinder(files.FileHandler.GetOGsTreeDir(), files.FileHandler.GetWorkingDirectory_Write(),
                                    self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        if qSpeciesTree:
            util.RenameTreeTaxa(spTreeFN_ids, files.FileHandler.GetSpeciesTreeUnrootedFN(True), self.ogSet.SpeciesDict(),
                                qSupport=False, qFixNegatives=True)
            return spTreeFN_ids, qSTAG
        else:      
            return None, qSTAG
            
    def SpeciesTreeOnly(self):
        qLessThanFourSpecies = len(self.ogSet.seqsInfo.speciesToUse) < 4
        if qLessThanFourSpecies:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            WriteSpeciesTreeIDs_TwoThree(self.ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
        else:
            ogs, ogMatrices_partial = self.GetOGMatrices_FullParallel()
            ogMatrices = self.CompleteOGMatrices(ogs, ogMatrices_partial)
            del ogMatrices_partial
            D, spPairs = self.SpeciesTreeDistances(ogs, ogMatrices)
            del ogMatrices
            cmd_spTree, spTreeFN_ids = self.PrepareSpeciesTreeCommand(D, spPairs, True)
            parallel_task_manager.RunCommand(cmd_spTree, True, False)
        spTreeUnrootedFN = files.FileHandler.GetSpeciesTreeUnrootedFN(True) 
        util.RenameTreeTaxa(spTreeFN_ids, spTreeUnrootedFN, self.ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True)  
        return spTreeFN_ids

# ==============================================================================================================================      
# DendroBlast   

def Worker_OGMatrices_ReadBLASTAndUpdateDistances(cmd_queue, worker_status_queue, iWorker, ogMatrices, nGenes, seqsInfo,
                                                  blastDir_list, ogsPerSpecies, qDoubleBlast):
    speciesToUse = seqsInfo.speciesToUse
    with np.errstate(divide='ignore'):
        while True:
            try:
                iiSp, sp1, nSeqs_sp1 = cmd_queue.get(True, 1)
                worker_status_queue.put(("start", iWorker, iiSp))
                Bs = [BlastFileProcessor.GetBLAST6Scores(seqsInfo, blastDir_list, sp1, sp2,
                                                         qExcludeSelfHits = False, qDoubleBlast=qDoubleBlast)
                      for sp2 in speciesToUse]
                mins = np.ones((nSeqs_sp1, 1), dtype=np.float64)*9e99 
                maxes = np.zeros((nSeqs_sp1, 1), dtype=np.float64)
                for B in Bs:
                    m0, m1 = lil_minmax(B)
                    mins = np.minimum(mins, m0)
                    maxes = np.maximum(maxes, m1)
                maxes_inv = 1./maxes
                for jjSp, B  in enumerate(Bs):
                    for og, m in zip(ogsPerSpecies, ogMatrices):
                        for gi, i in og[iiSp]:
                            for gj, j in og[jjSp]:
                                    m[i][j] = 0.5*max(B[gi.iSeq, gj.iSeq], mins[gi.iSeq]) * maxes_inv[gi.iSeq]
                del Bs, B, mins, maxes, m0, m1, maxes_inv    # significantly reduces RAM usage
                worker_status_queue.put(("finish", iWorker, iiSp))
            except queue.Empty:
                worker_status_queue.put(("empty", iWorker, None))
                return 

def GetRAMErrorText():
    text = "ERROR: The computer ran out of RAM and killed OrthoFinder processes\n"
    text += "Try using a computer with more RAM. If you used the '-a' option\n"
    text += "it may be possible to complete the run by removing this option."
    return text

# ==============================================================================================================================

def lil_min(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        mins[kRow] = min(values.data[0])
    return mins 

def lil_max(M):
    n = M.shape[0]
    maxes = np.zeros((n, 1), dtype = np.float64)
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        maxes[kRow] = max(values.data[0])
    return maxes

def lil_minmax(M):
    n = M.shape[0]
    mins = np.ones((n, 1), dtype = np.float64) * 9e99
    maxes = np.zeros((n, 1), dtype = np.float64)
    for kRow in range(n):
        values=M.getrowview(kRow)
        if values.nnz == 0:
            continue
        mins[kRow] = min(values.data[0])
        maxes[kRow] = max(values.data[0])
    return mins, maxes

# ==============================================================================================================================    
# Species trees for two- & three-species analyses

def WriteSpeciesTreeIDs_TwoThree(taxa, outFN):
    """
    Get the unrooted species tree for two or three species
    Args:
        taxa - list of species names
    Returns:
    
    """
    t = tree.Tree()
    for s in taxa:
        t.add_child(tree.TreeNode(name=s))
    t.write(outfile=outFN)
    
def GetSpeciesTreeRoot_TwoTaxa(taxa):
    speciesTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
    t = tree.Tree("(%s,%s);" % (taxa[0], taxa[1]))  
    t.write(outfile=speciesTreeFN_ids)
    return speciesTreeFN_ids

 