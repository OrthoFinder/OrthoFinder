# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologs
"""
from ..utils import logging as run_logging
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

from .tree_analyser import TreeAnalyser, mrca_node, gene_to_species_dash, get_lines_for_ortholog_files
from .hog_writer import HogWriter
from .parallel_pipeline import run_orthologs_parallel_pipeline, ParentOutputWriter, sort_non_hog_output_files, validate_hog_writer_no_duplicate_ids, set_file_descriptor_limit, NonHogAppendWriter
from ..phylo import tree
from ..utils import util, file_handler


@run_logging.RunLogger.stage("Infer orthologues")
def reconciliation_and_orthologs(
        recon_method,
        og_set,
        n_high_parallel,
        n_low_parallel,
        i_species_tree=None,
        stride_dups=None,
        q_split_para_clades=False,
        fewer_open_files=False,
        save_space=False,
        print_info=True,
        exist_msa=True,
        write_hog_tree=True,
        fix_files=True,
        write_to_rd=True,
        fd_limit=None,
        sort_non_hog_output=True,
        validate_hog_ids=True,
        write_named_hogs=True,
    ):

    species_tree_ids_fn = file_handler.handler.get_species_tree_ids_rooted_fn()
    labeled_tree_fn = file_handler.handler.get_species_tree_results_node_labels_fn()

    util.rename_tree_taxa(
        species_tree_ids_fn,
        labeled_tree_fn,
        og_set.species_dict(),
        q_support=False,
        q_fix_negatives=True,
        label="N"
    )

    working_dir = file_handler.handler.get_working_directory_write()
    results_dir_orthologs = file_handler.handler.get_orthologs_directory()
    recon_trees_renamed_dir = file_handler.handler.get_ogs_recon_tree_dir(True)

    if print_info:
        util.print_time("Inferring orthologues from gene trees")

    q_no_recon = ("only_overlap" == recon_method)

    # Label species-tree internal nodes deterministically.
    species_tree_rooted_labelled = tree.Tree(species_tree_ids_fn)
    species_tree_rooted_labelled.name = "N0"

    i_node = 1
    node_names = [species_tree_rooted_labelled.name]

    for n in species_tree_rooted_labelled.traverse():
        if (not n.is_leaf()) and (not n.is_root()):
            n.name = "N%d" % i_node
            node_names.append(n.name)
            i_node += 1

    species_dict = og_set.species_dict()
    sequence_dict = og_set.sequence_dict()

    if fd_limit is not None:
        set_file_descriptor_limit(fd_limit)

    hog_writer = HogWriter(
        species_tree_rooted_labelled,
        node_names,
        sequence_dict,
        species_dict,
        og_set.species_to_use,
        write_to_rd=write_to_rd,
        write_named=write_named_hogs,
    )

    try:
        n_orthologs_sp_pair = do_orthologs_for_ortho_finder(
            og_set,
            species_tree_rooted_labelled,
            gene_to_species_dash,
            stride_dups,
            q_no_recon,
            hog_writer,
            q_split_para_clades,
            n_low_parallel,
            fewer_open_files,
            save_space,
            print_info=print_info,
            exist_msa=exist_msa,
            write_hog_tree=write_hog_tree,
            fix_files=fix_files,
            fd_limit=fd_limit,
        )

        if print_info:
            util.print_time("Done of orthologues")

        two_and_three_gene_hogs(
            og_set,
            species_tree_rooted_labelled,
            hog_writer
        )

        hog_writer.close_files()

        if validate_hog_ids:
            validate_hog_writer_no_duplicate_ids(hog_writer)

        if not write_hog_tree or not fix_files:
            n_orthologs_sp_pair += two_and_three_gene_orthogroups(
                og_set,
                results_dir_orthologs,
                save_space,
                fewer_open_files,
                write_hog_tree=write_hog_tree,
                fix_files=fix_files
            )

            if sort_non_hog_output:
                sort_non_hog_output_files(
                    n_low_parallel,
                    og_set.species_to_use,
                    species_dict,
                    fewer_open_files,
                    save_space,
                    write_hog_tree,
                    fix_files,
                )

    finally:
        hog_writer.close_files()

    return n_orthologs_sp_pair



def two_and_three_gene_hogs(og_set, st_rooted_labelled, hog_writer):
    ogs = og_set.get_ogs_all()
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n < 2 or n > 3: continue
        og_name = "OG%07d" % iog
        sp_present = set([str(g.i_sp) for g in og])
        st_node = mrca_node(st_rooted_labelled, sp_present)
        hogs_to_write = hog_writer.get_skipped_nodes(st_node, None)  
        if len(sp_present) > 1:
            # We don't create files for 'species specific HOGs'
            st_node = mrca_node(st_rooted_labelled, sp_present)
            hogs_to_write = hogs_to_write + [st_node.name]
        # sorted so the row does not depend on the gene order of the orthogroup
        genes = sorted((g.to_string() for g in og), key=util.seq_id_key)
        hog_writer.write_hog_genes(genes, hogs_to_write, og_name)


def two_and_three_gene_orthogroups(
        og_set,
        results_dir,
        save_space,
        fewer_open_files,
        write_hog_tree=False,
        fix_files=False
    ):

    if write_hog_tree and fix_files:
        return util.NOrthoSp(len(og_set.species_to_use))

    species_dict = og_set.species_dict()
    sequence_dict = og_set.sequence_dict()

    all_orthologs, nspecies, sp_to_index, _old_ortholog_lines, _old_sus_lines = get_all_orthologs(og_set)

    # Important:
    # In save_space mode, output is also one file per species,
    # so get_lines_for_ortholog_files() must use dim2 = 1.
    use_fewer_open_files = fewer_open_files or save_space
    dim2 = 1 if use_fewer_open_files else nspecies

    ortholog_lines_tot = [
        ["" for _ in range(dim2)]
        for _ in range(nspecies)
    ]

    ortholog_sus_lines_tot = [
        ""
        for _ in range(nspecies)
    ]

    n_orthologs_sp_pair = get_lines_for_ortholog_files(
        all_orthologs,
        species_dict,
        og_set.species_to_use,
        sequence_dict,
        q_contains_suspect_orthologs=False,
        ortholog_lines=ortholog_lines_tot,
        ortholog_sus_lines=ortholog_sus_lines_tot,
        fewer_open_files=use_fewer_open_files,
    )

    writer = NonHogAppendWriter(
        d_results_orthologs=results_dir,
        species_dict=species_dict,
        species_to_use=og_set.species_to_use,
        save_space=save_space,
        fewer_open_files=fewer_open_files,
        max_open=64,
    )

    try:
        writer.write_ortholog_lines(ortholog_lines_tot)
        writer.flush()

    finally:
        writer.close()

    return n_orthologs_sp_pair


def get_all_orthologs(og_set):
    ogs = og_set.get_ogs_all()
    all_orthologs = []
    d_empty = defaultdict(list)
    for iog, og in enumerate(ogs):
        n = len(og) 
        if n == 1: break
        elif n == 2:
            if og[0].i_sp == og[1].i_sp: continue
            # orthologs is a list of tuples of dictionaries
            # each dictionary is sp->list of genes in species
            d0 = defaultdict(list)
            d0[str(og[0].i_sp)].append(str(og[0].i_seq))
            d1 = defaultdict(list)
            d1[str(og[1].i_sp)].append(str(og[1].i_seq))
            orthologs = [(d0, d1, d_empty, d_empty)]
        elif n == 3:
            sp = [g.i_sp for g in og]
            c = Counter(sp) 
            n_sp = len(c)
            if n_sp == 3:
                g = [(str(g.i_sp), str(g.i_seq)) for g in og]
                d0 = defaultdict(list)
                d0[g[0][0]].append(g[0][1])
                d1 = defaultdict(list)
                d1[g[1][0]].append(g[1][1])
                d1[g[2][0]].append(g[2][1])
                orthologs = [(d0, d1, d_empty, d_empty)]  
                d0 = defaultdict(list)
                d0[g[1][0]].append(g[1][1])
                d1 = defaultdict(list)
                d1[g[2][0]].append(g[2][1])
                orthologs.append((d0,d1, d_empty, d_empty))
            elif n_sp == 2:             
                sp0, sp1 = list(c.keys())
                d0 = defaultdict(list)
                d0[str(sp0)] = [str(g.i_seq) for g in og if g.i_sp == sp0]
                d1 = defaultdict(list)
                d1[str(sp1)] = [str(g.i_seq) for g in og if g.i_sp == sp1]
                orthologs = [(d0, d1, d_empty, d_empty)]
            else: 
                continue # no orthologs
            
        else:
            continue
        if orthologs:
            all_orthologs.append((iog, orthologs))

    nspecies = len(og_set.species_to_use)
    sp_to_index = {str(sp):i for i, sp in enumerate(og_set.species_to_use)}
    ortholog_lines_tot = [["" for j in range(nspecies)] for i in range(nspecies)]
    ortholog_sus_lines_tot = ["" for i in range(nspecies)]
    return all_orthologs, nspecies, sp_to_index, ortholog_lines_tot, ortholog_sus_lines_tot


def do_orthologs_for_ortho_finder(
        og_set,
        species_tree_rooted_labelled,
        GeneToSpecies,
        stride_dups,
        q_no_recon,
        hog_writer,
        q_split_paralogous_clades,
        n_parallel,
        fewer_open_files,
        save_space,
        print_info=True,
        exist_msa=True,
        write_hog_tree=True,
        fix_files=True,
        fd_limit=None,
    ):
    try:
        species_dict = og_set.species_dict()
        sequence_dict = og_set.sequence_dict()
        nspecies = len(og_set.species_to_use)

        d_results_orthologs = file_handler.handler.get_orthologs_directory()

        neighbours = get_species_neighbours(species_tree_rooted_labelled)
        iogs4 = list(og_set.get_i_ogs4())
        recon_trees_renamed_dir = file_handler.handler.get_ogs_recon_tree_dir(True)
        spec_seq_dict = og_set.spec_seq_dict()
        effective_fewer_open_files = fewer_open_files or save_space

        output_writer = ParentOutputWriter(
            d_results_orthologs=d_results_orthologs,
            species_dict=species_dict,
            species_to_use=og_set.species_to_use,
            sequence_dict=sequence_dict,
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
            d_results_orthologs,
            recon_trees_renamed_dir,
            species_tree_rooted_labelled,
            og_set.species_to_use,
            GeneToSpecies,
            sequence_dict,
            species_dict,
            spec_seq_dict,
            neighbours,
            q_no_recon,
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

            n_orthologs_sp_pair = util.NOrthoSp(nspecies)

            try:
                for completed_tasks, iog in enumerate(iogs4, start=1):
                    result = ta.analyse_tree(iog)

                    if result is None:
                        progressbar.update(task, advance=1)
                        continue

                    n_orthologs_sp_pair += result["n_orthologues"]
                    output_writer.write_result(result)

                    progressbar.update(task, advance=1)

            finally:
                progressbar.stop()
                output_writer.close()

            if print_info:
                util.print_time("Done writing orthologs")

        else:
            args_queue = mp.Queue()

            # Start the largest trees first so one big tree does not run alone
            # at the end. HOG rows are still committed in OG order.
            ogs_all = og_set.get_ogs_all()
            for iog in sorted(iogs4, key=lambda i: -len(ogs_all[i])):
                args_queue.put(iog)

            n_orthologs_sp_pair = run_orthologs_parallel_pipeline(
                ta,
                len(og_set.species_to_use),
                args_queue,
                n_parallel,
                total_tasks,
                effective_fewer_open_files,
                output_writer,
                iogs_ordered=iogs4,
                n_orthologs_cache=100,
                write_hog_tree=write_hog_tree,
                fix_files=fix_files,
                fd_limit=fd_limit,
                writer_queue_size=max(4 * n_parallel, 32),
                spill_dir=file_handler.handler.get_working_directory_write(),
            )


    except IOError as e:
        if str(e).startswith("[Errno 24] Too many open files"):
            util.number_open_files_exception_advice(len(og_set.species_to_use), True)
            util.fail()
        else:
            raise

    return n_orthologs_sp_pair


def get_species_neighbours(t):
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

