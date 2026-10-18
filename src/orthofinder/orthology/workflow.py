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
from ..utils import logging as run_logging
import os
import shutil
try:
    from rich import print
except ImportError:
    ...
from ..phylo import tree
from ..gene_tree import dendroblast
from ..species_tree import stride
from ..gene_tree import msa_trees
# from ..gene_tree_inference import trees2ologs_dlcpar
from . import reconciliation, hog_writer as hog_processor, parallel_pipeline
from ..gene_tree import infer_trees
from ..utils import util, file_handler
from ..mcl_orthogroups import orthogroups_set
from . import stats

from ..hog_orthogroups import finalise
from ..utils import printing
from ..utils import gene_ids
from ..utils import paths
from ..phylo import taxa
# from ..mcl_orthogroups import assign

# ==============================================================================================================================      
# Main
            

        
# def write_orthologs_matrix(fn, matrix, speciesToUse, speciesDict):
#     with open(fn, util.csv_write_mode) as outfile:
#         writer = csv.writer(outfile, delimiter="\t")
#         writer.writerow([""] + [speciesDict[str(index)] for index in speciesToUse])
#         for ii, iSp in enumerate(speciesToUse):
#             overlap = [matrix[ii, jj] for jj, jSp in enumerate(speciesToUse)]
#             writer.writerow([speciesDict[str(iSp)]] + overlap)   
    

# def write_orthologs_stats(ogSet, NOrthoSp):
#     """
#     NOrthoSp is a stats.NOrthoSp object
#     """
#     speciesToUse = ogSet.speciesToUse
#     speciesDict = ogSet.species_dict()
#     d = files.FileHandler.get_ogs_stats_results_directory()

#     # print(NOrthoSp.n)
#     # print(NOrthoSp.n_121)
#     # print(NOrthoSp.n_12m)
#     # print( NOrthoSp.n_m21)
#     # print(NOrthoSp.n_m2m)
    
#     write_orthologs_matrix(d + "OrthologsStats_Totals.tsv", NOrthoSp.n, speciesToUse, speciesDict)
#     write_orthologs_matrix(d + "OrthologsStats_one-to-one.tsv", NOrthoSp.n_121, speciesToUse, speciesDict)
#     write_orthologs_matrix(d + "OrthologsStats_one-to-many.tsv", NOrthoSp.n_12m, speciesToUse, speciesDict)
#     write_orthologs_matrix(d + "OrthologsStats_many-to-one.tsv", NOrthoSp.n_m21, speciesToUse, speciesDict)
#     write_orthologs_matrix(d + "OrthologsStats_many-to-many.tsv", NOrthoSp.n_m2m, speciesToUse, speciesDict)
    
#     # Duplications
#     nodeCount = defaultdict(int)
#     nodeCount_50 = defaultdict(int)
#     ogCount = defaultdict(int)
#     ogCount_50 = defaultdict(int)
#     if not os.path.exists(files.FileHandler.get_duplications_fn()): return
#     with open(files.FileHandler.get_duplications_fn(), util.csv_read_mode) as infile:
#         reader = csv.reader(infile, delimiter="\t")
#         next(reader)
#         # for line in reader:
#         #     try:
#         #         og, node, _, support, _, _, _ = line
#         #     except:
#         #         print(line)
#         #         raise
#         for og, node, _, support, _, _, _ in reader:
#             support = float(support)
#             nodeCount[node] += 1
#             ogCount[og] += 1
#             if support >= 0.5:
#                 nodeCount_50[node] += 1
#                 ogCount_50[og] += 1
#     with open(d + "Duplications_per_Species_Tree_Node.tsv", util.csv_write_mode) as outfile:
#         writer = csv.writer(outfile, delimiter="\t")
#         writer.writerow(["Species Tree Node", "Duplications (all)", "Duplications (50% support)"])
# #        max_node = max([int(s[1:]) for s in nodeCount.keys()])    # Get largest node number
#         for node in nodeCount:
#             writer.writerow([node, nodeCount[node], nodeCount_50[node]])
#     # Write on species tree
#     in_tree_fn = files.FileHandler.get_species_tree_results_node_labels_fn()
#     out_tree_fn = os.path.split(files.FileHandler.get_duplications_fn())[0] + "/SpeciesTree_Gene_Duplications_0.5_Support.txt"
#     t = tree.Tree(in_tree_fn, format=1)
#     for n in t.traverse():
#         n.name = n.name + "_" + str(nodeCount_50[n.name])
#     with open(out_tree_fn, 'w') as outfile:
#         outfile.write(t.write(format=1)[:-1] + t.name + ";")
#     with open(d + "Duplications_per_Orthogroup.tsv", util.csv_write_mode) as outfile:
#         writer = csv.writer(outfile, delimiter="\t")
#         writer.writerow(["Orthogroup", "Duplications (all)", "Duplications (50% support)"])
#         if len(ogCount) > 0:
#             max_og = max([int(s[2:]) for s in ogCount.keys()]) 
#             pat = files.FileHandler.baseOgFormat 
#             for i in range(max_og + 1):
#                 og = pat % i
#                 writer.writerow([og, ogCount[og], ogCount_50[og]])

def orthologs_workflow(
        seqs_info, 
        species_names_dict, 
        species_info_obj, 
        options,
        species_to_use, 
        n_sp_all, 
        program_caller,
        msa_method,
        tree_method,
        recon_method,
        n_high_parallel,
        n_low_parallel,
        q_double_blast,
        q_add_species_to_ids,
        q_trim,
        cmd_order="ascending",
        method_threads=None,
        user_species_tree = None,
        q_stop_after_seqs = False,
        q_stop_after_align = False,
        q_stop_after_trees = False, 
        q_stop_after_species_trees = False, 
        q_msa = False,
        q_phyldog = False,
        results_name = "",
        q_split_para_clades=False,
        save_space=False,
        root_from_previous=False,
        i_og_restart=0,
    ):
    og_set = orthogroups_set.OrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(), 
        species_to_use, 
        n_sp_all,
        q_add_species_to_ids, 
        tree_program=options.tree_program,
        id_extractor = gene_ids.FirstWordExtractor,
        mcl_inflation=options.mcl_inflation
    )

    return_obj = infer_trees.infer_gene_and_species_trees(
        og_set,
        program_caller, 
        msa_method, 
        tree_method,
        n_high_parallel, 
        n_low_parallel, 
        q_double_blast, 
        q_add_species_to_ids, 
        q_trim, 
        cmd_order,
        method_threads,
        user_species_tree, 
        q_stop_after_seqs, 
        q_stop_after_align, 
        q_msa, 
        q_phyldog,
        results_name, 
        root_from_previous, 
        i_og_restart,
        options.species_tree_program,
        n_skip=options.n_skip,
        msa_min_seq=options.msa_min_seq,
        species_tree_alignment=options.species_tree_alignment,
    )
    
    file_handler.handler.log_working_directory_trees()
    
    if q_stop_after_trees:
        print()
        print(f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] stopped after gene tree...")
        # (the empty results directories are removed as the run ends: finish.remove_empty_results_dirs)

        printing.print_run_end(file_handler.handler.get_results_directory1(), "-fgt", "its gene trees")
        util.success()
    
    if return_obj is None:
        return
    sp_tree_fn_ids, q_species_tree_supports = return_obj
    
    return_obj = root_species_tree(
        og_set, 
        sp_tree_fn_ids, 
        q_species_tree_supports,
        n_high_parallel, 
        n_low_parallel,
        user_species_tree, 
        q_stop_after_seqs, 
        q_stop_after_align, 
        q_stop_after_species_trees, 
        q_msa, 
        q_phyldog,
        results_name, 
        q_split_para_clades, 
        save_space, 
        root_from_previous
    )

    if q_stop_after_species_trees:
        print()
        print(f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] stopped after species tree...")
        # (the empty results directories are removed as the run ends: finish.remove_empty_results_dirs)
        printing.print_run_end(file_handler.handler.get_results_directory1(), "-fst", "its species tree")
        util.success()

    if return_obj is None:
        return
    rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, stride_dups = return_obj

    infer_orthologs_and_write_results(
        og_set, rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, q_species_tree_supports, stride_dups,
        seqs_info, species_names_dict, species_info_obj, options, species_to_use, recon_method,
        n_high_parallel, n_low_parallel, user_species_tree, q_phyldog,
        q_split_para_clades, root_from_previous, i_og_restart=i_og_restart,
    )


def infer_orthologs_and_write_results(
        og_set, rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, q_species_tree_supports, stride_dups,
        seqs_info, species_names_dict, species_info_obj, options, species_to_use, recon_method,
        n_high_parallel, n_low_parallel, user_species_tree, q_phyldog,
        q_split_para_clades, root_from_previous, i_og_restart=0, prev_wd=None,
    ):
    """
    From the rooted species tree to the results: orthologs and HOGs, then
    the results files updated to the HOG-based orthogroups and the statistics.
    The end of every analysis, whether run straight through (orthologs_workflow)
    or restarted from the species tree (-fst, orthologs_from_gene_species_trees),
    so that both give the same results.

    The first pass writes only to the working directory (the HOGs, from which
    the orthogroups are made); the second writes the results files, each once.
    """
    # --keep-legacy: the legacy files are written (WorkingDirectory/Legacy), else not at all
    file_handler.handler.keep_legacy = options.keep_legacy
    # First pass: the HOGs, in the working directory
    infer_orthologs(
        og_set, 
        rooted_sp_tree, 
        fn_rooted_sp_tree, 
        q_multiple_roots, 
        q_species_tree_supports, 
        stride_dups,
        recon_method,
        n_high_parallel, 
        n_low_parallel, 
        user_species_tree, 
        q_phyldog,
        q_split_para_clades, 
        root_from_previous,
        save_space=False,
        exist_msa=options.q_msa_trees,
        write_hog_tree=True,
        write_to_rd=False,
        fd_limit=options.fd_limit,
        # The update of the results files needs only N0_mcl.ids.tsv (n0_mcl_ids_fn);
        # the HOGs with gene names of this pass are legacy files (--keep-legacy: Legacy/HOGs)
        write_named_hogs=options.keep_legacy,
    )

    util.print_time("Converting MSA/Trees")
    og_set, new_ogs = finalise.update_output_files(
        og_set.species_dict(),
        og_set.sequence_dict(),
        og_set.species_to_use,
        og_set.all_used_sequence_ids(),
        species_info_obj,
        seqs_info,
        species_names_dict,
        options,
        options.n_process_alg,
        q_incremental=False,
        i_og_restart=i_og_restart,
        exist_msa=options.q_msa_trees,
        prev_wd=prev_wd,
    )
    util.print_time("Done updating MSA/Trees")
    paths.clear_dir(file_handler.handler.get_ogs_recon_tree_dir(q_results=True))
    # Second pass: the results files
    infer_orthologs(
        og_set,
        rooted_sp_tree,
        fn_rooted_sp_tree,
        q_multiple_roots,
        q_species_tree_supports,
        stride_dups,
        recon_method,
        n_high_parallel,
        n_low_parallel,
        user_species_tree,
        q_phyldog,
        q_split_para_clades,
        root_from_previous,
        save_space=options.save_space,
        print_info=False,
        exist_msa=options.q_msa_trees,
        write_hog_tree=False,
        write_to_rd=True,
        fd_limit=options.fd_limit
    )

    # (the sequences of the single-copy orthogroups only, read when they are known:
    # not every sequence of every species)
    fasta_writer = lambda genes: msa_trees.FastaWriter(file_handler.handler.get_species_seqs_dir(), species_to_use, genes)
    ogs = stats.add_unassigned_genes(new_ogs, og_set.all_used_sequence_ids())
    species_dict = {int(k): v for k, v in og_set.species_dict().items()}
    ids_dict = og_set.spec_seq_dict()   # (the single-copy FASTA files: <species>_<gene>)

    if not options.keep_legacy:
        os.remove(file_handler.handler.ogs_all_id_fn())
        os.remove(file_handler.handler.n0_mcl_ids_fn())
        shutil.rmtree(file_handler.handler.get_resolved_tree_id_dir())

    stats.stats(ogs, species_dict, species_to_use, file_handler.handler.i_results_version, fasta_writer, ids_dict)


@run_logging.RunLogger.stage("Infer species tree", require_result=True)
def root_species_tree(
        og_set, 
        sp_tree_fn_ids, 
        q_species_tree_supports,
        n_high_parallel,
        n_low_parallel,
        user_species_tree=None,
        q_stop_after_seqs=False,
        q_stop_after_align=False,
        q_stop_after_species_trees=False,
        q_msa=False,
        q_phyldog=False,
        results_name="",
        q_split_para_clades=False,
        save_space=False,
        root_from_previous=False,
    ):
    """
    SpeciesTree
    spTreeFN_ids, or equivalently FileHandler.get_species_tree_unrooted_fn() in all cases (user, inferred etc)
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
    if q_phyldog:
        rooted_species_tree_fn = [sp_tree_fn_ids]
        roots = [None]
        q_multiple = False
        stride_dups = None
    elif user_species_tree:
        rooted_species_tree_fn = [sp_tree_fn_ids]
        roots = [None]
        q_multiple = False
        stride_dups = None
    elif len(og_set.seqs_info.species_to_use) == 2:
        hardcode_species_tree = dendroblast.get_species_tree_root_two_taxa(og_set.seqs_info.species_to_use)
        rooted_species_tree_fn = [hardcode_species_tree]
        roots = [None]
        q_multiple = False
        stride_dups = None
    else:
        printing.print_underline("Best outgroup(s) for species tree") 
        util.print_time("Starting STRIDE")
        roots, clusters_counter, rooted_species_tree_fn, n_support, _, _, stride_dups = \
            stride.get_root(
                sp_tree_fn_ids, 
                file_handler.handler.get_ogs_tree_dir(), 
                stride.gene_to_species_dash, 
                n_high_parallel, 
                q_write_rooted_tree=True
            )
        util.print_time("Done STRIDE")
        n_all = sum(clusters_counter.values())
        n_fp_mp = n_all - n_support
        n_non_trivial = sum([v for k, v in clusters_counter.items() if len(k) > 1])
        if len(roots) > 1:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best roots and %d contradict them." % (n_non_trivial, n_non_trivial-n_fp_mp, n_fp_mp)))
            print("Best outgroups for species tree:")  
        else:
            print(("Observed %d well-supported, non-terminal duplications. %d support the best root and %d contradict it." % (n_non_trivial, n_non_trivial-n_fp_mp, n_fp_mp)))
            print("Best outgroup for species tree:")  
        sp_dict = og_set.species_dict()
        for r in roots: 
            print(("  " + (", ".join([sp_dict[s] for s in r]))  ))
        q_multiple = len(roots) > 1
    shutil.copy(rooted_species_tree_fn[0], file_handler.handler.get_species_tree_ids_rooted_fn())

    """
    SpeciesTree:
    We now have a list of rooted species trees: rootedSpeciesTreeFN (this should be recorded by the file handler)
    """
    if q_stop_after_species_trees:
        # root the gene trees using the species tree and write out their accessions - really I could remove the whole '-ot, -os, -oa' options, they are probably rarely used if ever.
        if user_species_tree:
            return
        # otherwise, root species tree
        results_species_trees = []
        for i, (r, species_tree_fn) in enumerate(zip(roots, rooted_species_tree_fn)):
            results_species_trees.append(file_handler.handler.get_species_tree_results_fn(i, not q_multiple))
            taxa.rename_tree_taxa(species_tree_fn, results_species_trees[-1], og_set.species_dict(), q_support=q_species_tree_supports, q_fix_negatives=True)
            labeled_tree_fn = file_handler.handler.get_species_tree_results_node_labels_fn()
            taxa.rename_tree_taxa(species_tree_fn, labeled_tree_fn, og_set.species_dict(), q_support=False, q_fix_negatives=True, label='N')
        file_handler.handler.copy_concat_orthogroups_to_results()
        id_dict = og_set.spec_seq_dict()
        q_have_support = None 
        for iog in og_set.get_i_ogs4():
            infn = file_handler.handler.get_ogs_tree_fn(iog)
            if os.path.exists(infn):
                if q_have_support is None: q_have_support = taxa.have_support_values(infn)
                taxa.rename_tree_taxa(infn, file_handler.handler.get_ogs_tree_fn(iog, True), id_dict, q_support=q_have_support, q_fix_negatives=True)       
        file_handler.handler.clean_working_dir2()
        
    if q_multiple:
        print("\nWARNING: Multiple potential species tree roots were identified, only one will be analyed.")
        for i, (r, species_tree_fn) in enumerate(zip(roots, rooted_species_tree_fn)):
            unanalysed_species_tree = file_handler.handler.get_species_tree_results_fn(i, False)
            taxa.rename_tree_taxa(species_tree_fn, unanalysed_species_tree, og_set.species_dict(),
                                q_support=q_species_tree_supports, q_fix_negatives=True, label='N')
    return roots[0], rooted_species_tree_fn[0], q_multiple, stride_dups


def infer_orthologs(
        og_set, 
        rooted_sp_tree, 
        species_tree_fn, 
        q_multiple_species_tree_roots, 
        q_species_tree_supports, 
        stride_dups,
        recon_method,
        n_high_parallel,
        n_low_parallel,
        user_species_tree=None,
        q_phyldog=False,
        q_split_para_clades=False,
        root_from_previous=False,
        save_space=False,
        print_info=True,
        exist_msa=True,
        write_hog_tree=True,
        working_dir="",
        write_to_rd=True,
        fd_limit=None,
        write_named_hogs=True,
    ):
    """ C. Gene tree rooting & orthologs
    write_to_rd: write the results files (the species tree's too), else only the working directory's
    write_named_hogs: write the HOG files with gene names (else only N0.ids.tsv)"""

    i_rooted_sp_tree = 0
    if print_info:
        printing.print_underline("Reconciling gene trees and species tree")
    if (not user_species_tree) and (not q_phyldog) and len(og_set.seqs_info.species_to_use) != 2:
        if print_info:
            print(("Outgroup: " + (", ".join([og_set.species_dict()[s] for s in rooted_sp_tree]))))
    if write_to_rd:
        file_handler.handler.copy_concat_orthogroups_to_results()
        taxa.rename_tree_taxa(
            species_tree_fn,
            file_handler.handler.get_species_tree_results_fn(0, True),
            og_set.species_dict(),
            q_support=q_species_tree_supports,
            q_fix_negatives=True
        )
    if print_info:
        util.print_time("Starting reconciliation and orthologs")

    n_orthologs_sp_pair = reconciliation.reconciliation_and_orthologs(
        recon_method, 
        og_set, 
        n_high_parallel, 
        n_low_parallel, 
        i_rooted_sp_tree if q_multiple_species_tree_roots else None,
        stride_dups=stride_dups, 
        q_split_para_clades=q_split_para_clades,
        save_space=save_space, 
        print_info=print_info,
        exist_msa=exist_msa,
        write_hog_tree=write_hog_tree,
        write_to_rd=write_to_rd,
        fd_limit=fd_limit,
        write_named_hogs=write_named_hogs,
    )
    # util.print_time("Done Recon")

    if not write_hog_tree:
        # print("%fs for orthologs etc" % (stop-start))
        stats.write_orthologs_stats(og_set, n_orthologs_sp_pair)
    #    print("Identified %d orthologs" % nOrthologs)

    file_handler.handler.clean_working_dir2(working_dir=working_dir)
    if print_info:
        printing.print_underline("Writing results files", True)
        # util.print_time("Start writing files")
        util.print_time("Start writing hierarchical orthogroups files")



def orthologs_from_gene_trees(
        seqs_info, 
        species_names_dict, 
        species_info_obj, 
        options,
        species_to_use, 
        n_sp_all, 
        recon_method,
        n_high_parallel,
        n_low_parallel,
        q_add_species_to_ids,
        user_species_tree = None,
        q_stop_after_seqs = False,
        q_stop_after_align = False,
        q_stop_after_trees = False, 
        q_msa = False,
        q_phyldog = False,
        results_name = "",
        q_split_para_clades=False,
        save_space=False,
        root_from_previous=False,
        i_og_restart=0,
):
    """
    Restart from the gene trees (-fgt): the gene trees of a run stopped after
    them (-ogt) are used. The species tree (inferred from them, or the one
    the run used, -s / --assign) is rooted as in a run straight through, then
    the analysis ends as one does (infer_orthologs_and_write_results).
    """
    og_set = orthogroups_set.OrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(), 
        species_to_use, 
        n_sp_all,
        q_add_species_to_ids, 
        tree_program=options.tree_program,
        id_extractor = gene_ids.FirstWordExtractor,
        mcl_inflation=options.mcl_inflation
    )

    sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
    prev_wd = os.path.dirname(sp_tree_fn_ids)
    file_handler.handler.log_working_directory_trees()
    if user_species_tree is not None:
        # as infer_gene_and_species_trees does: saved as 'unrooted', used as rooted
        infer_trees.convert_user_species_tree(user_species_tree, og_set.species_dict(), sp_tree_fn_ids)
    q_species_tree_supports = (
        False if (user_species_tree or q_msa or q_phyldog)
        else taxa.have_support_values(sp_tree_fn_ids)
    )

    return_obj = root_species_tree(
        og_set, 
        sp_tree_fn_ids, 
        q_species_tree_supports,
        n_high_parallel, 
        n_low_parallel,
        user_species_tree, 
        q_stop_after_seqs=q_stop_after_seqs,
        q_stop_after_align=q_stop_after_align,
        q_msa=q_msa,
        q_phyldog=q_phyldog,
        results_name=results_name,
        q_split_para_clades=q_split_para_clades,
        save_space=save_space,
        root_from_previous=root_from_previous,
    )
    if return_obj is None:
        return
    rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, stride_dups = return_obj
    infer_orthologs_and_write_results(
        og_set, rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, q_species_tree_supports, stride_dups,
        seqs_info, species_names_dict, species_info_obj, options, species_to_use, recon_method,
        n_high_parallel, n_low_parallel, user_species_tree, q_phyldog,
        q_split_para_clades, root_from_previous, i_og_restart=i_og_restart,
        prev_wd=prev_wd,
    )


def orthologs_from_gene_species_trees(
        seqs_info, 
        species_names_dict, 
        species_info_obj, 
        options,
        species_to_use, 
        n_sp_all, 
        recon_method,
        n_high_parallel,
        n_low_parallel,
        q_add_species_to_ids,
        user_species_tree_fn,
        q_split_para_clades=False,
        i_og_restart=0,
):
    """
    Restart from the species tree (-fst): the gene trees and the (unrooted)
    species tree of a run stopped after the species tree (-ost) are used.
    The species tree is rooted as in a run straight through (the rooting is
    deterministic, and STRIDE's duplications are needed for the results),
    then the analysis ends as one does (infer_orthologs_and_write_results).
    """
    og_set = orthogroups_set.OrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(), 
        species_to_use, 
        n_sp_all,
        q_add_species_to_ids, 
        tree_program=options.tree_program,
        id_extractor = gene_ids.FirstWordExtractor,
        mcl_inflation=options.mcl_inflation
    )
    
    sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
    prev_wd = os.path.dirname(sp_tree_fn_ids)

    if user_species_tree_fn is not None:
        # as infer_gene_and_species_trees does: saved as 'unrooted', used as rooted
        infer_trees.convert_user_species_tree(user_species_tree_fn, og_set.species_dict(), sp_tree_fn_ids)
    q_species_tree_supports = (
        False if (user_species_tree_fn or options.q_msa_trees or options.q_phyldog)
        else taxa.have_support_values(sp_tree_fn_ids)
    )
    return_obj = root_species_tree(
        og_set,
        sp_tree_fn_ids,
        q_species_tree_supports,
        n_high_parallel,
        n_low_parallel,
        user_species_tree_fn,
        q_msa=options.q_msa_trees,
        q_phyldog=options.q_phyldog,
        results_name=options.name,
        q_split_para_clades=q_split_para_clades,
        save_space=options.save_space,
    )
    if return_obj is None:
        return
    rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, stride_dups = return_obj
    infer_orthologs_and_write_results(
        og_set, rooted_sp_tree, fn_rooted_sp_tree, q_multiple_roots, q_species_tree_supports, stride_dups,
        seqs_info, species_names_dict, species_info_obj, options, species_to_use, recon_method,
        n_high_parallel, n_low_parallel, user_species_tree_fn, options.q_phyldog,
        q_split_para_clades, False, i_og_restart=i_og_restart, prev_wd=prev_wd,
    )
