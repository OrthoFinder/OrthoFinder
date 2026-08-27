# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import sys
import csv
import glob
import resource
import operator
import itertools
import traceback
import multiprocessing as mp
from collections import Counter, defaultdict
import warnings
import queue
import time
import io
try:
    from rich import print
except ImportError:
    ...

from .tree_processor import TreeAnalyser, MRCA_node, GeneToSpecies_dash, GetLinesForOlogFiles
from .hog_processor import HogWriter
from .trees2ologs_processor import \
    (
        RunOrthologsParallel_Pipeline, 
        ParentOutputWriter, 
        SortNonHogOutputFiles, 
        ValidateHogWriterNoDuplicateIds,
        set_file_descriptor_limit,
        NonHogAppendWriter
    )
from ..tools import tree
from ..utils import util, files


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
        print_info=True,
        exist_msa=True,
        write_hog_tree=True,
        fix_files=True,
        write_to_rd=True,
        fd_limit=None,
        sort_non_hog_output=True,
        validate_hog_ids=True,
    ):

    speciesTree_ids_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
    labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()

    util.RenameTreeTaxa(
        speciesTree_ids_fn,
        labeled_tree_fn,
        ogSet.SpeciesDict(),
        qSupport=False,
        qFixNegatives=True,
        label="N"
    )

    workingDir = files.FileHandler.GetWorkingDirectory_Write()
    resultsDir_ologs = files.FileHandler.GetOrthologuesDirectory()
    reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)

    if print_info:
        util.PrintTime("Inferring orthologues from gene trees")

    qNoRecon = ("only_overlap" == recon_method)

    # Label species-tree internal nodes deterministically.
    species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
    species_tree_rooted_labelled.name = "N0"

    iNode = 1
    node_names = [species_tree_rooted_labelled.name]

    for n in species_tree_rooted_labelled.traverse():
        if (not n.is_leaf()) and (not n.is_root()):
            n.name = "N%d" % iNode
            node_names.append(n.name)
            iNode += 1

    speciesDict = ogSet.SpeciesDict()
    SequenceDict = ogSet.SequenceDict()

    if fd_limit is not None:
        set_file_descriptor_limit(fd_limit)

    hog_writer = HogWriter(
        species_tree_rooted_labelled,
        node_names,
        SequenceDict,
        speciesDict,
        ogSet.speciesToUse,
        write_to_rd=write_to_rd
    )

    try:
        nOrthologues_SpPair = DoOrthologuesForOrthoFinder(
            ogSet,
            species_tree_rooted_labelled,
            GeneToSpecies_dash,
            stride_dups,
            qNoRecon,
            hog_writer,
            q_split_para_clades,
            nLowParallel,
            fewer_open_files,
            save_space,
            old_version=old_version,
            print_info=print_info,
            exist_msa=exist_msa,
            write_hog_tree=write_hog_tree,
            fix_files=fix_files,
            fd_limit=fd_limit,
        )

        if print_info:
            util.PrintTime("Done of orthologues")

        TwoAndThreeGeneHOGs(
            ogSet,
            species_tree_rooted_labelled,
            hog_writer
        )

        hog_writer.close_files()

        if validate_hog_ids:
            ValidateHogWriterNoDuplicateIds(hog_writer)

        if not write_hog_tree or not fix_files:
            nOrthologues_SpPair += TwoAndThreeGeneOrthogroups(
                ogSet,
                resultsDir_ologs,
                save_space,
                fewer_open_files,
                write_hog_tree=write_hog_tree,
                fix_files=fix_files
            )

            if sort_non_hog_output:
                SortNonHogOutputFiles(
                    nLowParallel,
                    ogSet.speciesToUse,
                    speciesDict,
                    fewer_open_files,
                    save_space,
                    write_hog_tree,
                    fix_files,
                )

    finally:
        hog_writer.close_files()

    return nOrthologues_SpPair



def TwoAndThreeGeneHOGs(ogSet, st_rooted_labelled, hog_writer):
    ogs = ogSet.OGsAll()
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n < 2 or n > 3: continue
        og_name = "OG%07d" % iog
        sp_present = set([str(g.iSp) for g in og])
        stNode = MRCA_node(st_rooted_labelled, sp_present)
        hogs_to_write = hog_writer.get_skipped_nodes(stNode, None)  
        if len(sp_present) > 1:
            # We don't create files for 'species specific HOGs'
            st_node = MRCA_node(st_rooted_labelled, sp_present)
            hogs_to_write = hogs_to_write + [st_node.name]
        genes = [g.ToString() for g in og] # Inefficient as will convert back again, but trivial cost I think
        hog_writer.write_hog_genes(genes, hogs_to_write, og_name)


def TwoAndThreeGeneOrthogroups(
        ogSet,
        resultsDir,
        save_space,
        fewer_open_files,
        write_hog_tree=False,
        fix_files=False
    ):

    if write_hog_tree and fix_files:
        return util.nOrtho_sp(len(ogSet.speciesToUse))

    speciesDict = ogSet.SpeciesDict()
    sequenceDict = ogSet.SequenceDict()

    all_orthologues, nspecies, sp_to_index, _old_olog_lines, _old_sus_lines = AllOrthologues(ogSet)

    # Important:
    # In save_space mode, output is also one file per species,
    # so GetLinesForOlogFiles() must use dim2 = 1.
    use_fewer_open_files = fewer_open_files or save_space
    dim2 = 1 if use_fewer_open_files else nspecies

    olog_lines_tot = [
        ["" for _ in range(dim2)]
        for _ in range(nspecies)
    ]

    olog_sus_lines_tot = [
        ""
        for _ in range(nspecies)
    ]

    nOrthologues_SpPair = GetLinesForOlogFiles(
        all_orthologues,
        speciesDict,
        ogSet.speciesToUse,
        sequenceDict,
        qContainsSuspectOlogs=False,
        olog_lines=olog_lines_tot,
        olog_sus_lines=olog_sus_lines_tot,
        fewer_open_files=use_fewer_open_files,
    )

    writer = NonHogAppendWriter(
        dResultsOrthologues=resultsDir,
        speciesDict=speciesDict,
        speciesToUse=ogSet.speciesToUse,
        save_space=save_space,
        fewer_open_files=fewer_open_files,
        max_open=64,
    )

    try:
        writer.write_ortholog_lines(olog_lines_tot)
        writer.flush()

    finally:
        writer.close()

    return nOrthologues_SpPair


def AllOrthologues(ogSet):
    ogs = ogSet.OGsAll()
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
            
        else:
            continue
        if orthologues:
            all_orthologues.append((iog, orthologues))

    nspecies = len(ogSet.speciesToUse)
    sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}
    olog_lines_tot = [["" for j in range(nspecies)] for i in range(nspecies)]
    olog_sus_lines_tot = ["" for i in range(nspecies)]
    return all_orthologues, nspecies, sp_to_index, olog_lines_tot, olog_sus_lines_tot


def DoOrthologuesForOrthoFinder(
        ogSet,
        species_tree_rooted_labelled,
        GeneToSpecies,
        stride_dups,
        qNoRecon,
        hog_writer,
        q_split_paralogous_clades,
        n_parallel,
        fewer_open_files,
        save_space,
        old_version=False,
        print_info=True,
        exist_msa=True,
        write_hog_tree=True,
        fix_files=True,
        fd_limit=None,
    ):
    try:
        speciesDict = ogSet.SpeciesDict()
        SequenceDict = ogSet.SequenceDict()
        nspecies = len(ogSet.speciesToUse)

        dResultsOrthologues = files.FileHandler.GetOrthologuesDirectory()

        neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)
        iogs4 = list(ogSet.Get_iOGs4())
        reconTreesRenamedDir = files.FileHandler.GetOGsReconTreeDir(True)
        spec_seq_dict = ogSet.Spec_SeqDict()
        effective_fewer_open_files = fewer_open_files or save_space

        output_writer = ParentOutputWriter(
            dResultsOrthologues=dResultsOrthologues,
            speciesDict=speciesDict,
            speciesToUse=ogSet.speciesToUse,
            SequenceDict=SequenceDict,
            spec_seq_dict=spec_seq_dict,
            stride_dups=stride_dups,
            hog_writer=hog_writer,
            fewer_open_files=effective_fewer_open_files,
            save_space=save_space,
            write_hog_tree=write_hog_tree,
            fix_files=fix_files,
            max_open=64,
            flush_every_results=1000,
            flush_every_seconds=60,
        )

        ta = TreeAnalyser(
            len(iogs4),
            dResultsOrthologues,
            reconTreesRenamedDir,
            species_tree_rooted_labelled,
            ogSet.speciesToUse,
            GeneToSpecies,
            SequenceDict,
            speciesDict,
            spec_seq_dict,
            neighbours,
            qNoRecon,
            None,       
            stride_dups,
            None,  
            None,  
            hog_writer,
            q_split_paralogous_clades,
            fewer_open_files=effective_fewer_open_files,
            exist_msa=exist_msa,
            write_hog_tree=write_hog_tree,
            fix_files=fix_files
        )

        total_tasks = len(iogs4)

        if n_parallel == 1:
            progressbar, task = util.get_progressbar(total_tasks)
            progressbar.start()

            nOrthologues_SpPair = util.nOrtho_sp(nspecies)

            try:
                for completed_tasks, iog in enumerate(iogs4, start=1):
                    result = ta.AnalyseTree(iog)

                    if result is None:
                        progressbar.update(task, advance=1)
                        continue

                    nOrthologues_SpPair += result["n_orthologues"]
                    output_writer.write_result(result)

                    progressbar.update(task, advance=1)

            finally:
                progressbar.stop()
                output_writer.close()

            if print_info:
                util.PrintTime("Done writing orthologs")

        else:
            args_queue = mp.Queue()

            for iog in iogs4:
                args_queue.put(iog)

            nOrthologues_SpPair = RunOrthologsParallel_Pipeline(
                ta,
                len(ogSet.speciesToUse),
                args_queue,
                n_parallel,
                total_tasks,
                effective_fewer_open_files,
                output_writer,
                iogs_ordered=iogs4,
                n_ologs_cache=100,
                compatibility_mode=old_version,
                write_hog_tree=write_hog_tree,
                fix_files=fix_files,
                fd_limit=fd_limit,
                writer_queue_size=max(4 * n_parallel, 32),
            )


    except IOError as e:
        if str(e).startswith("[Errno 24] Too many open files"):
            util.number_open_files_exception_advice(len(ogSet.speciesToUse), True)
            util.Fail()
        else:
            raise

    return nOrthologues_SpPair


def GetSpeciesNeighbours(t):
    """
    Args: t = rooted species tree
    
    Returns:
    dict: species -> species_dict, such that species_dict: other_species -> toplogical_dist 
    """
    species = t.get_leaf_names()
    levels = {s:[] for s in species}
    for n in t.traverse('postorder'):
        if n.is_leaf(): continue
        children = n.get_children()
        leaf_sets = [set(ch.get_leaf_names()) for ch in children]
        not_i = [set.union(*[l for j, l in enumerate(leaf_sets) if j != i]) for i in range(len(children))]
        for l,n in zip(leaf_sets, not_i):
            for ll in l:
                levels[ll].append(n)
    neighbours = {sp:{other:n for n,others in enumerate(lev) for other in others} for sp, lev in levels.items()}
    return neighbours


