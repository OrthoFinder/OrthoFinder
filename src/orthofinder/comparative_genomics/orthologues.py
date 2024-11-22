#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.comhor: david
import os
import sys
import csv
import time
import shutil
import numpy as np
import subprocess
from collections import Counter, defaultdict
import itertools
import multiprocessing as mp
import warnings
try: 
    import queue
except ImportError:
    import Queue as queue

from ..tools import tree, wrapper_phyldog, stag, stride, trees_msa, dendroblast
from ..tools import mcl as MCL
from ..gene_tree_inference import trees2ologs_dlcpar, trees2ologs_of, infer_trees

from ..utils import util, files, parallel_task_manager, program_caller
from ..orthogroups import accelerate, orthogroups_set
from . import stats
from .. import nThreadsDefault as nThreads
from .. import my_env

   

# ==============================================================================================================================      
# Main
            
def CheckUserSpeciesTree(speciesTreeFN, expSpecies):
    # File exists
    if not os.path.exists(speciesTreeFN):
        print(("Species tree file does not exist: %s" % speciesTreeFN))
        util.Fail()
    # Species in tree are unique
    try:
        t = tree.Tree(speciesTreeFN, format=1)
    except Exception as e:
        print("\nERROR: Incorrectly formated user-supplied species tree")
        print(str(e))
        util.Fail()
    actSpecies = (t.get_leaf_names())
    c = Counter(actSpecies)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names in species tree are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.Fail()
    # All required species are present
    actSpecies = set(actSpecies)
    ok = True
    for sp in expSpecies:
        if sp not in actSpecies:
            print(("ERROR: '%s' is missing from species tree" % sp))
            ok = False
    # expected species are unique
    c = Counter(expSpecies)
    if 1 != c.most_common()[0][1]:
        print("ERROR: Species names are not unique")
        for sp, n in c.most_common():
            if 1 != n:
                print(("Species '%s' appears %d times" % (sp, n)))
        util.Fail()
    expSpecies = set(expSpecies)
    for sp in actSpecies:
        if sp not in expSpecies:
            print(("ERROR: Additional species '%s' in species tree" % sp))
            ok = False
    if not ok: util.Fail()
    # Tree is rooted
    if len(t.get_children()) != 2:
        print("ERROR: Species tree is not rooted")
        util.Fail()

        
def WriteOrthologuesMatrix(fn, matrix, speciesToUse, speciesDict):
    with open(fn, util.csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([""] + [speciesDict[str(index)] for index in speciesToUse])
        for ii, iSp in enumerate(speciesToUse):
            overlap = [matrix[ii, jj] for jj, jSp in enumerate(speciesToUse)]
            writer.writerow([speciesDict[str(iSp)]] + overlap)   
    

def WriteOrthologuesStats(ogSet, nOrtho_sp):
    """
    nOrtho_sp is a util.nOrtho_sp object
    """
    speciesToUse = ogSet.speciesToUse
    speciesDict = ogSet.SpeciesDict()
    d = files.FileHandler.GetOGsStatsResultsDirectory()

    # print(nOrtho_sp.n)
    # print(nOrtho_sp.n_121)
    # print(nOrtho_sp.n_12m)
    # print( nOrtho_sp.n_m21)
    # print(nOrtho_sp.n_m2m)
    
    WriteOrthologuesMatrix(d + "OrthologuesStats_Totals.tsv", nOrtho_sp.n, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_one-to-one.tsv", nOrtho_sp.n_121, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_one-to-many.tsv", nOrtho_sp.n_12m, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_many-to-one.tsv", nOrtho_sp.n_m21, speciesToUse, speciesDict)
    WriteOrthologuesMatrix(d + "OrthologuesStats_many-to-many.tsv", nOrtho_sp.n_m2m, speciesToUse, speciesDict)
    
    # Duplications
    nodeCount = defaultdict(int)
    nodeCount_50 = defaultdict(int)
    ogCount = defaultdict(int)
    ogCount_50 = defaultdict(int)
    if not os.path.exists(files.FileHandler.GetDuplicationsFN()): return
    with open(files.FileHandler.GetDuplicationsFN(), util.csv_read_mode) as infile:
        reader = csv.reader(infile, delimiter="\t")
        next(reader)
        # for line in reader:
        #     try:
        #         og, node, _, support, _, _, _ = line
        #     except:
        #         print(line)
        #         raise
        for og, node, _, support, _, _, _ in reader:
            support = float(support)
            nodeCount[node] += 1
            ogCount[og] += 1
            if support >= 0.5:
                nodeCount_50[node] += 1
                ogCount_50[og] += 1
    with open(d + "Duplications_per_Species_Tree_Node.tsv", util.csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Species Tree Node", "Duplications (all)", "Duplications (50% support)"])
#        max_node = max([int(s[1:]) for s in nodeCount.keys()])    # Get largest node number
        for node in nodeCount:
            writer.writerow([node, nodeCount[node], nodeCount_50[node]])
    # Write on species tree
    in_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
    out_tree_fn = os.path.split(files.FileHandler.GetDuplicationsFN())[0] + "/SpeciesTree_Gene_Duplications_0.5_Support.txt"
    t = tree.Tree(in_tree_fn, format=1)
    for n in t.traverse():
        n.name = n.name + "_" + str(nodeCount_50[n.name])
    with open(out_tree_fn, 'w') as outfile:
        outfile.write(t.write(format=1)[:-1] + t.name + ";")
    with open(d + "Duplications_per_Orthogroup.tsv", util.csv_write_mode) as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Orthogroup", "Duplications (all)", "Duplications (50% support)"])
        if len(ogCount) > 0:
            max_og = max([int(s[2:]) for s in ogCount.keys()]) 
            pat = files.FileHandler.baseOgFormat 
            for i in range(max_og + 1):
                og = pat % i
                writer.writerow([og, ogCount[og], ogCount_50[og]])

def TwoAndThreeGeneHOGs(ogSet, st_rooted_labelled, hog_writer):
    ogs = ogSet.OGsAll()
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n < 2 or n > 3: continue
        og_name = "OG%07d" % iog
        sp_present = set([str(g.iSp) for g in og])
        stNode = trees2ologs_of.MRCA_node(st_rooted_labelled, sp_present)
        hogs_to_write = hog_writer.get_skipped_nodes(stNode, None)  
        if len(sp_present) > 1:
            # We don't create files for 'species specific HOGs'
            st_node = trees2ologs_of.MRCA_node(st_rooted_labelled, sp_present)
            hogs_to_write = hogs_to_write + [st_node.name]
        genes = [g.ToString() for g in og] # Inefficient as will convert back again, but trivial cost I think
        hog_writer.write_hog_genes(genes, hogs_to_write, og_name)

def TwoAndThreeGeneOrthogroups(ogSet, resultsDir, save_space, fewer_open_files):
    speciesDict = ogSet.SpeciesDict()
    sequenceDict = ogSet.SequenceDict()
    ogs = ogSet.OGsAll()
    nOrthologues_SpPair = util.nOrtho_sp(len(ogSet.speciesToUse))
    all_orthologues = []
    d_empty = defaultdict(list)
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n == 1: break
        elif n == 2:
            if og[0].iSp == og[1].iSp: continue
            # orthologues is a list of tuples of dictionaries
            # each dictionary is sp->list of genes in species
            d0 = defaultdict(list)
            d0[str(og[0].iSp)].append(str(og[0].iSeq))
            d1 = defaultdict(list)
            d1[str(og[1].iSp)].append(str(og[1].iSeq))
            orthologues = [(d0, d1, d_empty, d_empty)]  
        elif n == 3:
            sp = [g.iSp for g in og]
            c = Counter(sp) 
            nSp = len(c)
            if nSp == 3:
                g = [(str(g.iSp), str(g.iSeq)) for g in og]
                d0 = defaultdict(list)
                d0[g[0][0]].append(g[0][1])
                d1 = defaultdict(list)
                d1[g[1][0]].append(g[1][1])
                d1[g[2][0]].append(g[2][1])
                orthologues = [(d0, d1, d_empty, d_empty)]  
                d0 = defaultdict(list)
                d0[g[1][0]].append(g[1][1])
                d1 = defaultdict(list)
                d1[g[2][0]].append(g[2][1])
                orthologues.append((d0,d1, d_empty, d_empty))
            elif nSp == 2:             
                sp0, sp1 = list(c.keys())
                d0 = defaultdict(list)
                d0[str(sp0)] = [str(g.iSeq) for g in og if g.iSp == sp0]
                d1 = defaultdict(list)
                d1[str(sp1)] = [str(g.iSeq) for g in og if g.iSp == sp1]
                orthologues = [(d0, d1, d_empty, d_empty)]
            else: 
                continue # no orthologues
        elif n >= 4:
            continue
        all_orthologues.append((iog, orthologues))
    nspecies = len(ogSet.speciesToUse)
    sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}
    with trees2ologs_of.OrthologsFiles(resultsDir, speciesDict, ogSet.speciesToUse, nspecies,
                                       sp_to_index, save_space, fewer_open_files) as (olog_files_handles, suspect_genes_file_handles):
        olog_lines_tot = [["" for j in range(nspecies)] for i in range(nspecies)]
        olog_sus_lines_tot = ["" for i in range(nspecies)]
        nOrthologues_SpPair += trees2ologs_of.GetLinesForOlogFiles(all_orthologues, speciesDict, ogSet.speciesToUse, sequenceDict, 
                                                                   False, olog_lines_tot, olog_sus_lines_tot, fewer_open_files=fewer_open_files)
        # olog_sus_lines_tot will be empty
        lock_dummy = mp.Lock()
        for i in range(nspecies):
            for j in range(nspecies):
                trees2ologs_of.WriteOlogLinesToFile(olog_files_handles[i][j], olog_lines_tot[i][j], lock_dummy)
    return nOrthologues_SpPair
    
def ReconciliationAndOrthologues(
        recon_method,
        ogSet,
        nHighParallel,
        nLowParallel,
        iSpeciesTree=None,
        stride_dups=None,
        q_split_para_clades=False,
        fewer_open_files=False,
        save_space=False,
        old_version=False,
):
    """
    ogSet - info about the orthogroups, species etc.
    resultsDir - where the Orthologues top level results directory will go (should exist already)
    reconTreesRenamedDir - where to put the reconcilled trees that use the gene accessions
    iSpeciesTree - which of the potential roots of the species tree is this
    method - can be dlcpar, dlcpar_deep, of_recon
    """
    speciesTree_ids_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
    labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
    util.RenameTreeTaxa(speciesTree_ids_fn, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
    workingDir = files.FileHandler.GetWorkingDirectory_Write()    # workingDir - Orthologues working dir
    resultsDir_ologs = files.FileHandler.GetOrthologuesDirectory()
    reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)
    if "dlcpar" in recon_method:
        qDeepSearch = (recon_method == "dlcpar_convergedsearch")
        util.PrintTime("Starting DLCpar")
        dlcparResultsDir, dlcparLocusTreePat = trees2ologs_dlcpar.RunDlcpar(ogSet, speciesTree_ids_fn, workingDir, nHighParallel, qDeepSearch)
        util.PrintTime("Done DLCpar")
        spec_seq_dict = ogSet.Spec_SeqDict()
        ogs_all = ogSet.OGsAll()
        for iog, og in enumerate(ogs_all):
            if len(og) < 4:
                # For dlpar analysis can rely on ordered orthogroups
                break
            util.RenameTreeTaxa(dlcparResultsDir + dlcparLocusTreePat % iog, files.FileHandler.GetOGsReconTreeFN(iog),
                                spec_seq_dict, qSupport=False, qFixNegatives=False, inFormat=8, label='n')
    
        # Orthologue lists
        util.PrintUnderline("Inferring orthologues from gene trees" + (" (root %d)"%iSpeciesTree if iSpeciesTree != None else ""))
        pickleDir = files.FileHandler.GetPickleDir()
        nOrthologues_SpPair = trees2ologs_dlcpar.create_orthologue_lists(ogSet, resultsDir_ologs, dlcparResultsDir, pickleDir)  

    elif "phyldog" == recon_method:
        util.PrintTime("Starting Orthologues from Phyldog")
        nOrthologues_SpPair = trees2ologs_of.DoOrthologuesForOrthoFinder_Phyldog(ogSet, workingDir, trees2ologs_of.GeneToSpecies_dash,
                                                                                 resultsDir_ologs, reconTreesRenamedDir)
        util.PrintTime("Done Orthologues from Phyldog")
    else:
        start = time.time()
        util.PrintTime("Starting OF Orthologues")
        qNoRecon = ("only_overlap" == recon_method)
        # The next function should not create the HOG writer and label the species tree. This should be done here and passed as arguments
        species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
        # Label nodes of species tree
        species_tree_rooted_labelled.name = "N0"    
        iNode = 1
        node_names = [species_tree_rooted_labelled.name]
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = "N%d" % iNode
                node_names.append(n.name)
                iNode += 1
        # HOG Writer
        speciesDict = ogSet.SpeciesDict()
        SequenceDict = ogSet.SequenceDict()
        hog_writer = trees2ologs_of.HogWriter(species_tree_rooted_labelled, node_names, SequenceDict, speciesDict, ogSet.speciesToUse)
        nOrthologues_SpPair = trees2ologs_of.DoOrthologuesForOrthoFinder(ogSet, species_tree_rooted_labelled, trees2ologs_of.GeneToSpecies_dash, 
                                                                         stride_dups, qNoRecon, hog_writer, q_split_para_clades, nLowParallel,
                                                                         fewer_open_files, save_space, old_version=old_version)
        util.PrintTime("Done OF Orthologues")
        TwoAndThreeGeneHOGs(ogSet, species_tree_rooted_labelled, hog_writer)
        hog_writer.close_files()
    nOrthologues_SpPair += TwoAndThreeGeneOrthogroups(ogSet, resultsDir_ologs, save_space=save_space, fewer_open_files=fewer_open_files)
    if nLowParallel > 1 and "phyldog" != recon_method and "dlcpar" not in recon_method:
        trees2ologs_of.SortParallelFiles(nLowParallel, ogSet.speciesToUse, speciesDict, fewer_open_files)
    stop = time.time()
    # print("%fs for orthologs etc" % (stop-start))
    WriteOrthologuesStats(ogSet, nOrthologues_SpPair)
#    print("Identified %d orthologues" % nOrthologues)
        
                
def OrthologuesFromTrees(
        recon_method,
        nHighParallel,
        nLowParallel,
        userSpeciesTree_fn,
        qAddSpeciesToIDs,
        q_split_para_clades,
        fewer_open_files,
        old_version=False,
):
    """
    userSpeciesTree_fn - None if not supplied otherwise rooted tree using user species names (not orthofinder IDs)
    qUserSpTree - is the speciesTree_fn user-supplied
    
    Just infer orthologues from trees, don't do any of the preceeding steps.
    """
    speciesToUse, nSpAll, _ = util.GetSpeciesToUse(files.FileHandler.GetSpeciesIDsFN())    
    ogSet = orthogroups_set.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesToUse, nSpAll, qAddSpeciesToIDs,
                           idExtractor = util.FirstWordExtractor)
    if userSpeciesTree_fn != None:
        speciesDict = files.FileHandler.GetSpeciesDict()
        speciesToUseNames = [speciesDict[str(iSp)] for iSp in ogSet.speciesToUse]
        CheckUserSpeciesTree(userSpeciesTree_fn, speciesToUseNames)
        speciesTreeFN_ids = files.FileHandler.GetSpeciesTreeIDsRootedFN()
        infer_trees.ConvertUserSpeciesTree(userSpeciesTree_fn, speciesDict, speciesTreeFN_ids)
    util.PrintUnderline("Running Orthologue Prediction", True)
    util.PrintUnderline("Reconciling gene and species trees") 
    ReconciliationAndOrthologues(recon_method, ogSet, nHighParallel, nLowParallel, q_split_para_clades=q_split_para_clades,
                                 fewer_open_files=fewer_open_files, old_version=old_version)
    util.PrintUnderline("Writing results files")
    util.PrintTime("Writing results files")
    files.FileHandler.CleanWorkingDir2()


def OrthologuesWorkflow(speciesToUse, nSpAll, 
                       program_caller,
                       msa_method,
                       tree_method,
                       recon_method,
                       nHighParallel,
                       nLowParallel,
                       qDoubleBlast,
                       qAddSpeciesToIDs,
                       qTrim,
                       fewer_open_files,  # Open one ortholog file per species when analysing trees
                       cmd_order="ascending",
                       method_threads=None,
                       method_threads_large=None,
                       method_threads_small=None, 
                       threshold=None,
                       old_version=False,
                       userSpeciesTree = None,
                       qStopAfterSeqs = False,
                       qStopAfterAlign = False,
                       qStopAfterTrees = False, 
                       qMSA = False,
                       qPhyldog = False,
                       results_name = "",
                       q_split_para_clades=False,
                       save_space=False,
                       root_from_previous=False,
                       i_og_restart=0,
):
    ogSet = orthogroups_set.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesToUse, nSpAll,
                           qAddSpeciesToIDs, idExtractor = util.FirstWordExtractor)

    return_obj = infer_trees.InferGeneAndSpeciesTrees(ogSet,
                       program_caller, msa_method, tree_method,
                       nHighParallel, nLowParallel, qDoubleBlast, 
                       qAddSpeciesToIDs, qTrim, 
                       cmd_order,
                       method_threads,
                       method_threads_large,
                       method_threads_small, 
                       threshold,
                       old_version,
                       userSpeciesTree, qStopAfterSeqs, 
                       qStopAfterAlign, qMSA, qPhyldog,
                       results_name, root_from_previous, i_og_restart)
    if return_obj is None:
        return
    spTreeFN_ids, qSpeciesTreeSupports = return_obj

    return_obj = RootSpeciesTree(ogSet, spTreeFN_ids, qSpeciesTreeSupports,
                       nHighParallel, nLowParallel,
                       userSpeciesTree, qStopAfterSeqs, qStopAfterAlign, qStopAfterTrees, qMSA, qPhyldog,
                       results_name, q_split_para_clades, save_space, root_from_previous)
    if return_obj is None:
        return
    rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, stride_dups = return_obj

    InferOrthologs(ogSet, 
                   rooted_sp_tree, 
                   fn_rooted_sp_tree, 
                   q_multiple_roots, 
                   qSpeciesTreeSupports, 
                   stride_dups,
                   recon_method,
                    nHighParallel, 
                    nLowParallel, 
                    fewer_open_files,
                    userSpeciesTree, 
                    qPhyldog,
                    q_split_para_clades, 
                    save_space, 
                    root_from_previous,
                    old_version=old_version)

    fastaWriter = trees_msa.FastaWriter(files.FileHandler.GetSpeciesSeqsDir(), speciesToUse)
    ogs = accelerate.read_hogs(files.FileHandler.GetResultsDirectory1(), "N0")
    ogs = stats.add_unassigned_genes(ogs, ogSet.AllUsedSequenceIDs())
    species_dict = {int(k): v for k, v in ogSet.SpeciesDict().items()}
    ids_dict = ogSet.SequenceDict()
    stats.Stats(ogs, species_dict, speciesToUse, files.FileHandler.iResultsVersion, fastaWriter, ids_dict)


def RootSpeciesTree(ogSet, spTreeFN_ids, qSpeciesTreeSupports,
                             nHighParallel,
                             nLowParallel,
                             userSpeciesTree=None,
                             qStopAfterSeqs=False,
                             qStopAfterAlign=False,
                             qStopAfterTrees=False,
                             qMSA=False,
                             qPhyldog=False,
                             results_name="",
                             q_split_para_clades=False,
                             save_space=False,
                             root_from_previous=False,
                             ):
    """
    SpeciesTree
    spTreeFN_ids, or equivalently FileHandler.GetSpeciesTreeUnrootedFN() in all cases (user, inferred etc)
    Thus, we always have the species tree ids format
    
    With phyldog, we also have species_tree_ids_labelled_phyldog - with the node labels given by phyldog
    """    
         
    """ === 3 ===
    MSA:               RootSpeciesTree
    Phyldog:           RootSpeciesTree    
    Dendroblast:       RootSpeciesTree  
    MSA (ust):         ConvertSpeciesTreeIDs
    Phyldog (ust):     ConvertSpeciesTreeIDs
    Dendroblast (ust): ConvertSpeciesTreeIDs
    """
    """ B. Root species tree"""
    if qPhyldog:
        rootedSpeciesTreeFN = [spTreeFN_ids]
        roots = [None]
        qMultiple = False
        stride_dups = None
    elif userSpeciesTree:
        rootedSpeciesTreeFN = [spTreeFN_ids]
        roots = [None]
        qMultiple = False
        stride_dups = None
    elif len(ogSet.seqsInfo.speciesToUse) == 2:
        hardcodeSpeciesTree = dendroblast.GetSpeciesTreeRoot_TwoTaxa(ogSet.seqsInfo.speciesToUse)
        rootedSpeciesTreeFN = [hardcodeSpeciesTree]
        roots = [None]
        qMultiple = False
        stride_dups = None
    else:
        util.PrintUnderline("Best outgroup(s) for species tree") 
        util.PrintTime("Starting STRIDE")
        roots, clusters_counter, rootedSpeciesTreeFN, nSupport, _, _, stride_dups = \
            stride.GetRoot(spTreeFN_ids, files.FileHandler.GetOGsTreeDir(), stride.GeneToSpecies_dash, nHighParallel, qWriteRootedTree=True)
        util.PrintTime("Done STRIDE")
        nAll = sum(clusters_counter.values())
        nFP_mp = nAll - nSupport
        n_non_trivial = sum([v for k, v in clusters_counter.items() if len(k) > 1])
        if len(roots) > 1:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best roots and %d contradict them." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
            print("Best outgroups for species tree:")  
        else:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best root and %d contradict it." % (n_non_trivial, n_non_trivial-nFP_mp, nFP_mp)))
            print("Best outgroup for species tree:")  
        spDict = ogSet.SpeciesDict()
        for r in roots: print(("  " + (", ".join([spDict[s] for s in r]))  ))
        qMultiple = len(roots) > 1
    shutil.copy(rootedSpeciesTreeFN[0], files.FileHandler.GetSpeciesTreeIDsRootedFN())

    """
    SpeciesTree:
    We now have a list of rooted species trees: rootedSpeciesTreeFN (this should be recorded by the file handler)
    """
    if qStopAfterTrees:
        # root the gene trees using the species tree and write out their accessions - really I could remove the whole '-ot, -os, -oa' options, they are probably rarely used if ever.
        if userSpeciesTree:
            return
        # otherwise, root species tree
        resultsSpeciesTrees = []
        for i, (r, speciesTree_fn) in enumerate(zip(roots, rootedSpeciesTreeFN)):
            resultsSpeciesTrees.append(files.FileHandler.GetSpeciesTreeResultsFN(i, not qMultiple))
            util.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], ogSet.SpeciesDict(), qSupport=qSpeciesTreeSupports, qFixNegatives=True)
            labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
            util.RenameTreeTaxa(speciesTree_fn, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
        idDict = ogSet.Spec_SeqDict()
        qHaveSupport = None 
        for iog in ogSet.Get_iOGs4():
            infn = files.FileHandler.GetOGsTreeFN(iog)
            if os.path.exists(infn):
                if qHaveSupport is None: qHaveSupport = util.HaveSupportValues(infn)
                util.RenameTreeTaxa(infn, files.FileHandler.GetOGsTreeFN(iog, True), idDict, qSupport=qHaveSupport, qFixNegatives=True)       
        files.FileHandler.CleanWorkingDir2()
    if qMultiple:
        print("\nWARNING: Multiple potential species tree roots were identified, only one will be analyed.")
        for i, (r, speciesTree_fn) in enumerate(zip(roots, rootedSpeciesTreeFN)):
            unanalysedSpeciesTree = files.FileHandler.GetSpeciesTreeResultsFN(i, False)
            util.RenameTreeTaxa(speciesTree_fn, unanalysedSpeciesTree, ogSet.SpeciesDict(),
                                qSupport=qSpeciesTreeSupports, qFixNegatives=True, label='N')
    return roots[0], rootedSpeciesTreeFN[0], qMultiple, stride_dups


def InferOrthologs( ogSet, rooted_sp_tree, speciesTree_fn, qMultipleSpeciesTreeRoots, qSpeciesTreeSupports, stride_dups,
                    recon_method,
                    nHighParallel,
                    nLowParallel,
                    fewer_open_files,  # Open one ortholog file per species when analysing trees
                    userSpeciesTree=None,
                    qPhyldog=False,
                    q_split_para_clades=False,
                    save_space=False,
                    root_from_previous=False,
                    old_version=False,
                    ):
    """ C. Gene tree rooting & orthologs"""

    resultsSpeciesTrees = []
    i_rooted_sp_tree = 0
    util.PrintUnderline("Reconciling gene trees and species tree")         
    resultsSpeciesTrees.append(files.FileHandler.GetSpeciesTreeResultsFN(0, True))
    if (not userSpeciesTree) and (not qPhyldog) and len(ogSet.seqsInfo.speciesToUse) != 2:
        print(("Outgroup: " + (", ".join([ogSet.SpeciesDict()[s] for s in rooted_sp_tree]))))
    util.RenameTreeTaxa(speciesTree_fn, resultsSpeciesTrees[-1], ogSet.SpeciesDict(), qSupport=qSpeciesTreeSupports, qFixNegatives=True)
    util.PrintTime("Starting Recon and orthologues")
    ReconciliationAndOrthologues(recon_method, ogSet, nHighParallel, nLowParallel, i_rooted_sp_tree if qMultipleSpeciesTreeRoots else None,
                                 stride_dups=stride_dups, q_split_para_clades=q_split_para_clades,
                                 fewer_open_files=fewer_open_files, save_space=save_space, old_version=old_version)
    # util.PrintTime("Done Recon")

    files.FileHandler.CleanWorkingDir2()
    util.PrintUnderline("Writing results files", True)
