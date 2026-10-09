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
# david_emms@hotmail.com

from __future__ import absolute_import

# first import parallel task manager to minimise RAM overhead for small processes
import multiprocessing as mp  # optional  (problems on OpenBSD)
import platform  # Y
import sys  # Y
from ..utils import printing
from ..utils import species_info
from ..utils import gene_ids
from ..utils import paths
from ..phylo import taxa
from ..orthology import ortholog_writers
from ..orthology import rooting

if __name__ == "__main__":
    if platform.system() == "Darwin":
        # https://github.com/davidemms/OrthoFinder/issues/570
        # https://github.com/davidemms/OrthoFinder/issues/663
        mp.set_start_method("fork")
    # else:
    #     # Should be more RAM efficient than fork and the time penalty
    #     # should be very small as we never try to create many processes
    #     mp.set_start_method('spawn')

import os  # Y
import glob

# os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised

import time
import copy  # Y
import time  # Y
import csv  # Y
import os.path  # Y
from ..utils import parallel_task_manager, file_handler, util, fasta_processor, file_io, run_state, orthoxml
from ..mcl_orthogroups import infer_orthogroups, orthogroups_set
from ..mcl_orthogroups import assign as assign
from orthofinder.utils import logging as run_logging


from ..mcl_orthogroups import mcl
from ..phylo import tree
from ..species_tree import astral
from ..orthology import reconciliation
from ..gene_tree import infer_trees
from . import process_args, check_dependencies, previous_run
from ..search import runner
from .. import orphan_genes_version, __version__, __location__
from ..orthology import workflow
from ..utils.printing import printer
from rich.markup import escape

try:
    from rich import print
except ImportError:
    ...

TEST_MODE = os.getenv("ORTHOFINDER_TEST_ISOLATE") == "1"

configfile_location = os.path.join(__location__, "run")
max_int = sys.maxsize
ok = False
while not ok:
    try:
        csv.field_size_limit(max_int)
        ok = True
    except OverflowError:
        max_int = int(max_int / 10)
sys.setrecursionlimit(10**6)

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)



"""
OrthoFinder
-------------------------------------------------------------------------------
"""


def orthoxml_species(options, species_info_obj):
    """
    With -x, the species information for Orthogroups.orthoxml and
    Orthogroups.html, one entry per species in the analysis; else None.
    """
    if not options.orthoxml:
        return None
    return orthoxml.species_for_analysis(
        file_handler.handler.get_species_ids_fn(), species_info_obj.species_to_use
    )


"""
Main
-------------------------------------------------------------------------------
"""


# 9
def get_orthologs(
    seqs_info,
    species_names_dict,
    species_info_obj,
    options,
    prog_caller,
    i_og_restart=0,
    species_xml=None,
):
    printing.print_underline("Analysing Orthogroups", True)
    workflow.orthologs_workflow(
        seqs_info,
        species_names_dict,
        species_info_obj,
        options,
        species_info_obj.species_to_use,
        species_info_obj.n_sp_all,
        prog_caller,
        options.msa_program,
        options.tree_program,
        options.recon_method,
        options.n_blast,
        options.n_process_alg,
        options.q_double_blast,
        options.q_add_species_to_ids,
        options.q_trim,
        options.fewer_open_files,
        options.cmd_order,
        options.method_threads,
        options.method_threads_large,
        options.method_threads_small,
        options.threshold,
        options.species_tree_fn,
        options.q_stop_after_seqs,
        options.q_stop_after_alignments,
        options.q_stop_after_trees,
        options.q_stop_after_species_trees,
        options.q_msa_trees,
        options.q_phyldog,
        options.name,
        options.q_split_para_clades,
        save_space=options.save_space,
        root_from_previous=False,
        i_og_restart=i_og_restart,
        species_xml=species_xml,
    )
    # util.print_time("Done writing files")


def warn_species_without_assigned_genes(ogs_new_species, species_info_obj, species_names_dict):
    """
    --assign: warn about new species none of whose genes were assigned to a
    core orthogroup (no hits in the orthogroup-profile search). Such a species
    is in no gene tree, so it cannot be placed in the species tree, and its
    genes stay unassigned: the run would otherwise end as if it had worked.
    """
    with_genes = {int(g.split("_", 1)[0]) for genes in ogs_new_species.values() for g in genes}
    missing = [i_sp for i_sp in range(species_info_obj.i_first_new_species, species_info_obj.n_sp_all)
               if i_sp not in with_genes]
    if missing:
        printer.print(
            "\nWARNING: No genes of these new species were assigned to the core orthogroups (no hits "
            "in the search against the orthogroup profiles): %s. They cannot be placed in the species "
            "tree, and their genes will be unassigned. Check that they are related to the core species "
            "and that the search program can compare these sequences (e.g. for nucleotide sequences, "
            "-d with -S blast or -S mmseqs).\n"
            % escape(", ".join(str(species_names_dict.get(i_sp, i_sp)) for i_sp in missing)),
            style="warning")


def between_core_orthogroups_workflow(
    continuation_dir,
    species_info_obj,
    seqs_info,
    options,
    prog_caller,
    species_names_dict,
    results_files,
    q_hogs,
    results_layout=None,
    only_missing=False,
):
    """
    Infer clade-specific orthogroups for the new species clades
    only_missing: continuing a run (-b): only the clade searches whose results are missing or incomplete
    results_layout: file_io.HitFormat - columns of the profile-search results
    n_unassigned: List[int] - number of unassigned genes per species
    """
    # Get current orthogroups - original orthogroups plus genes assigned to them
    if q_hogs:
        ogs = assign.read_hogs(continuation_dir, "N0")
    else:
        ogs = assign.get_original_orthogroups()
    i_og_restart = 0
    ogs_new_species, _ = assign.assign_genes(results_files, options.n_process_alg, layout=results_layout)
    warn_species_without_assigned_genes(ogs_new_species, species_info_obj, species_names_dict)
    clusters_filename_pairs = assign.write_all_orthogroups(
        ogs, ogs_new_species, []
    )  # this updates ogs

    og_set = orthogroups_set.OrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(),
        species_info_obj.species_to_use,
        species_info_obj.n_sp_all,
        options.q_add_species_to_ids,
        options.tree_program,
        id_extractor=gene_ids.FirstWordExtractor,
        mcl_inflation=options.mcl_inflation,
    )

    if options.q_stop_after_groups and options.species_tree_fn is None:
        # Can't infer clade-specific groups
        print("\nSpecies tree required for clade-speicfic orthogroups - skipping")
        return clusters_filename_pairs, i_og_restart

    n_unassigned = assign.write_unassigned_fasta(ogs, None, species_info_obj)

    # Get/Infer species tree
    if options.species_tree_fn is None:
        # Infer gene trees
        # We write orthogroup & stats results files in the following code, which we should avoid & only do once all OGs are done.
        infer_orthogroups.post_clustering_orthogroups(
            clusters_filename_pairs,
            species_info_obj,
            seqs_info,
            species_names_dict,
            options,
            species_xml=None,
            q_incremental=True,
        )

        infer_trees.infer_gene_and_species_trees(
            og_set,
            prog_caller,
            options.msa_program,
            options.tree_program,
            options.n_blast,
            options.n_process_alg,
            options.q_double_blast,
            options.q_add_species_to_ids,
            options.q_trim,
            cmd_order=options.cmd_order,
            method_threads=options.method_threads,
            method_threads_large=options.method_threads_large,
            method_threads_small=options.method_threads_small,
            threshold=options.threshold,
            user_species_tree=None,
            q_stop_after_seqs=False,
            q_stop_after_align=False,
            q_msa=options.q_msa_trees,
            q_phyldog=False,
            results_name=options.name,
            root_from_previous=True,
            n_skip=options.n_skip,
        )

        # Infer species tree
        astral_fn = file_handler.handler.get_astral_filename()
        astral.create_input_file(
            file_handler.handler.get_ogs_tree_dir(), astral_fn, n_skip=options.n_skip
        )
        species_tree_unrooted_fn = file_handler.handler.get_species_tree_unrooted_fn()
        parallel_task_manager.run_command(
            astral.get_astral_command(
                astral_fn, species_tree_unrooted_fn, options.n_blast
            ),
            # ASTRAL-Pro writes its progress log to stderr: only show it if it fails.
            q_print_stderr=False,
            raise_on_error=True,
        )

        # Root it
        core_rooted_species_tree = tree.Tree(
            file_handler.handler.get_core_species_tree_ids_rooted_fn(), format=1
        )
        species_to_speices_map = lambda x: x
        rooted_species_tree_ids, q_have_support = rooting.check_and_root_tree(
            species_tree_unrooted_fn, core_rooted_species_tree, species_to_speices_map
        )

        if rooted_species_tree_ids is None:
            print(
                "ERROR: Species tree inference failed. Please check for errors and check the species tree files: \n%s \n%s"
                % (species_tree_unrooted_fn, core_rooted_species_tree)
            )
            util.fail()

        rooted_species_tree_fn = file_handler.handler.get_species_tree_ids_rooted_fn()
        rooted_species_tree_ids.write(outfile=rooted_species_tree_fn)

        sp_tree_unrooted_fn = file_handler.handler.get_species_tree_results_fn(None, True)
        taxa.rename_tree_taxa(
            rooted_species_tree_ids,
            sp_tree_unrooted_fn,
            og_set.species_dict(),
            q_support=q_have_support,
            q_fix_negatives=True,
        )

        labeled_tree_fn = file_handler.handler.get_species_tree_results_node_labels_fn()
        taxa.rename_tree_taxa(
            rooted_species_tree_ids,
            labeled_tree_fn,
            og_set.species_dict(),
            q_support=False,
            q_fix_negatives=True,
            label="N",
        )
        i_og_restart = len(ogs)  # Need to process the clade-specific orthogroups only
    else:
        printing.print_underline("Using user-supplied species tree")
        sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
        infer_trees.convert_user_species_tree(
            options.species_tree_fn, og_set.species_dict(), sp_tree_fn_ids
        )
        rooted_species_tree_fn = sp_tree_fn_ids

    # Identify clades for clade-specific orthogroup inference
    i_species_core = set(species_info_obj.get_original_species())
    species_clades = assign.get_new_species_clades(rooted_species_tree_fn, i_species_core)
    printing.print_underline(
        "Identifying clade-specific orthogroups for the following clades:", q_heavy=True
    )
    species_dict = og_set.species_dict()
    for i, clade in enumerate(species_clades):
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
    print("")

    # The clade-specific step costs k^2 searches (and more for clustering) per
    # clade of k species: report its size and warn about very large clades.
    clade_rows = assign.clade_costs(
        [list(map(str, clade)) for clade in species_clades],
        map(str, i_species_core),
        n_genes={str(isp): n for isp, n in enumerate(n_unassigned)},
    )
    for warning in assign.report_clades(clade_rows, species_dict):
        run_logging.RunLogger.message(warning, level="WARNING")
    print("")

    if not species_clades:
        # Every new species falls in a part of the tree covered by the core:
        # there are no clade-specific orthogroups to infer.
        print("No new-species clades: skipping clade-specific orthogroup inference\n")
        clusters_filename_pairs, i_og_restart = assign.write_all_orthogroups(
            ogs, {}, [], restart_index=i_og_restart
        )
        return clusters_filename_pairs, i_og_restart

    # Clade-specific orthogroup inference
    runner.create_search_databases(
        species_info_obj, options, prog_caller, q_unassigned_genes=True
    )

    # provide list of clades, only run these searches (and only if the fasta files are non-empty)
    runner.run_search(
        options,
        species_info_obj,
        seqs_info,
        prog_caller,
        n_genes_per_species=n_unassigned,
        species_clades=species_clades,
        only_missing=only_missing,
    )
    # process the results files - only if they are present and non-empty
    options.v2_scores = True
    n_clades = len(species_clades)
    clusters_filename_pairs_unassigned_all = []
    for i_clade, clade in enumerate(species_clades):
        printing.print_underline(
            "OrthoFinder clutering on new species clade %d of %d"
            % (i_clade + 1, n_clades)
        )
        print(str(i_clade) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
        species_info_clade = copy.deepcopy(species_info_obj)
        species_info_clade.species_to_use = clade
        seqs_info_clade = species_info.seqs_info_recompute(seqs_info, clade)
        clusters_filename_pairs_unassigned = infer_orthogroups.do_orthogroups(
            options,
            species_info_clade,
            seqs_info_clade,
            species_names_dict,
            species_xml=None,
            i_unassigned=i_clade,
        )
        clusters_filename_pairs_unassigned_all.append(clusters_filename_pairs_unassigned)
    ogs_clade_specific_list = [
        mcl.get_predicted_ogs(filename)
        for filename in clusters_filename_pairs_unassigned_all
    ]

    # OGs have had assigned genes added to them already
    # Single-gene orthogroups removed here shift the positions of later ones,
    # so the restart index is adjusted along with them.
    clusters_filename_pairs, i_og_restart = assign.write_all_orthogroups(
        ogs, {}, ogs_clade_specific_list, restart_index=i_og_restart
    )
    return clusters_filename_pairs, i_og_restart


# def GetOrthologs_FromTrees(options):
#     orthologs.orthologs_from_trees(
#         options.min_seq,
#         options.recon_method,
#         options.nBlast,
#         options.nProcessAlg,
#         options.speciesTreeFN,
#         options.qAddSpeciesToIDs,
#         options.qSplitParaClades,
#         options.fewer_open_files,
#         exist_msa=options.qMSATrees,
#         fix_files=options.fix_files,
#         mclInflation=options.mclInflation
#     )


def main(args=None):
    file_handler.handler.reset()
    file_io.clear_caches()   # main() may run more than once in a process (e.g. tests)
    start = time.perf_counter()
    log = None
    current_step = "Initialisation"
    exit_code = 0
    try:
        if args is None:
            args = sys.argv[1:]
        input_args = args.copy()
        # Create PTM right at start
        ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        # prog_caller = get_program_caller()
        (
            prog_caller,
            options,
            fasta_dir,
            continuation_dir,
            results_dir_non_default,
            pickle_dir_non_default,
            user_specified_m,
        ) = process_args.process_args(args)

        printer.print(
            f"[bold dark_goldenrod]OrthoFinder[/bold dark_goldenrod] version [deep_sky_blue2]{__version__}[/deep_sky_blue2]",
            end="",
        )
        printer.print(
            " Copyright (C) 2014 [bold dark_goldenrod]David Emms[/bold dark_goldenrod]\n"
        )

        # A restart continues as the run would have gone on (run_state): with
        # its settings (programs, -d etc., unless given again). Restored first,
        # so that what follows (e.g. "-S blast" as blastn or blastp) and what is
        # reported and recorded agree with the settings the restart runs with.
        i_og_restart_state = 0
        state = {}
        changed = []
        q_restart = ((options.q_start_from_blast and not options.q_start_from_fasta) or options.q_start_from_groups
                     or options.q_start_from_trees or options.q_start_from_species_trees)
        if q_restart and continuation_dir is not None:
            _, state = run_state.find(continuation_dir)
            if not state:
                # no run_state.json (older versions, or the file removed): what
                # Log.txt records, read before this run rewrites it
                state = previous_run.state_from_log(continuation_dir)
            i_og_restart_state = state.get("assign", {}).get("i_og_restart", 0)
            changed = run_state.restore_settings(state, options, input_args)
            # the run's search program (restored above, or the same as the default)
            options.search_program_from_run = ("search_program" in state.get("settings", {})
                                               and not options.search_program_given)

        # "-S blast" is blastn or blastp (-d). Before the file handler is set
        # up, as it logs the search program.
        options = process_args.resolve_search_program(options, prog_caller)

        previous_run.initialise_file_handler(
            options,
            fasta_dir,
            continuation_dir,
            results_dir_non_default,
            pickle_dir_non_default,
        )


        # A restart keeps the record of the stages it continues from, not of
        # those it runs again (so running it again does not record them twice).
        orthogroup_stages = ("Infer orthogroups", "Infer clade-specific orthogroups")
        if options.q_start_from_blast and not options.q_start_from_fasta:
            file_io.trim_checkpoint(file_handler.handler.get_check_point_fn(), redo=orthogroup_stages)
        elif options.q_start_from_groups:
            file_io.trim_checkpoint(file_handler.handler.get_check_point_fn(), completed=[orthogroup_stages])
        elif options.q_start_from_trees or options.q_start_from_species_trees:
            file_io.trim_checkpoint(file_handler.handler.get_check_point_fn(),
                                    completed=[orthogroup_stages, ("Build alignments and gene trees",)])
        log = run_logging.Logger(
            file_handler.handler.get_check_point_fn(),
            fmt="%(asctime)s : %(message)s",
        )
        run_logging.RunLogger.set_logger(log)
        log.info(
            "Starting OrthoFinder v%s\n"
            "%d thread(s) for highly parallel tasks (BLAST searches etc.)\n"
            "%d thread(s) for OrthoFinder algorithm\n\n"
            "OrthoFinder version %s Copyright (C) 2014 David Emms\n\n"
            "Results directory:\n    %s\n",
            __version__,
            options.n_blast,
            options.n_process_alg,
            __version__,
            file_handler.handler.get_results_directory1(),
        )


        # print("Results directory: %s" % files.FileHandler.get_results_directory1())
        # printer.print("Results directory:", style="path")
        # printer.print(f"    [dark_cyan]{files.FileHandler.get_results_directory1()}")

        printer.print("[bold]Results directory:[/bold]")
        printer.print(f"    [dark_cyan]{file_handler.handler.get_results_directory1()}")

        # A restart also continues with the species tree the run used (unless
        # -s is given), and its point for restarting gene-tree inference
        # (--assign); its settings were restored above.
        results_dir = file_handler.handler.get_results_directory1()
        if q_restart:
            if changed:
                printer.print("Continuing with the settings of the run being continued: %s"
                              % escape(", ".join("%s=%s" % c for c in changed)))
                file_handler.handler.write_to_log("Settings of the run being continued: %s\n"
                                             % ", ".join("%s=%s" % c for c in changed))
            kept_tree = run_state.file(results_dir, state, "species_tree")
            if i_og_restart_state:
                # kept in this run's Log.txt too (a restart rewrites it)
                file_handler.handler.log_assign_first_orthogroup_to_infer(i_og_restart_state)
            if options.species_tree_fn is None and kept_tree is not None:
                options.species_tree_fn = kept_tree
                printer.print("Using the species tree of the run being continued:")
                printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(kept_tree))
        # -b on the results of an --assign analysis continues that analysis
        assign_resume = (options.q_start_from_blast and not options.q_start_from_fasta
                         and state.get("run", {}).get("analysis") == "assign")
        if assign_resume and not ("first_new_species" in state.get("assign", {})
                                  and "new_species" in state.get("assign", {})
                                  and "core_results" in state.get("dirs", {})):
            # stopped before the new species were added (e.g. while creating
            # the orthogroup profiles): there is nothing to continue from
            file_handler.handler.log_fail_and_exit(
                "ERROR: The '--assign' run in %s stopped before its new species were added, so it "
                "cannot be continued with -b. Run --assign again." % results_dir)
        run_state.record_run(results_dir, "assign" if options.q_fast_add else "core",
                             " ".join(["orthofinder"] + input_args))
        run_state.save_settings(results_dir, options)
        if options.species_tree_fn is not None:
            file_handler.handler.log_species_tree_used(run_state.keep_species_tree(results_dir, options.species_tree_fn))

        check_dependencies.check_dependencies(
            options,
            user_specified_m,
            prog_caller,
            file_handler.handler.get_working_directory1_read()[0],
            q_assign=assign_resume,
        )
        if continuation_dir is not None:
            printing.print_restart_info(options, continuation_dir)


        # if using previous Trees etc., check these are all present - Job for orthologs
        if options.q_start_from_blast and options.q_start_from_fasta:
            # 0. Check Files
            species_info_obj, species_to_use_names = previous_run.process_previous_files(
                file_handler.handler.get_working_directory1_read(), options.q_double_blast
            )
            # print(
            #     "\nAdding new species in %s to existing analysis in %s"
            #     % (fastaDir, continuationDir)
            # )
            printer.print(f"\nAdding new species in [dark_cyan]{fasta_dir}")
            printer.print(f"to existing analysis in [dark_cyan]{continuation_dir}")
            # 3.
            species_info_obj = fasta_processor.processes_new_fasta(
                fasta_dir, options.dna, species_info_obj, species_to_use_names
            )
            file_handler.handler.log_species()
            options = process_args.check_options(options, species_info_obj.species_to_use)
            # 4.
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )
            # 5.
            species_xml = orthoxml_species(options, species_info_obj)
            # 6.
            printing.print_underline("Dividing up work for BLAST for parallel processing")
            runner.create_search_databases(species_info_obj, options, prog_caller)
            # 7.
            current_step = "Sequence search"
            log.step(current_step, "Started")
            runner.run_search(options, species_info_obj, seqs_info, prog_caller)
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 8.
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            infer_orthogroups.do_orthogroups(
                options, species_info_obj, seqs_info, species_names_dict, species_xml
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9.
            if options.fix_files and not options.q_stop_after_mcl_groups:
                get_orthologs(
                    seqs_info,
                    species_names_dict,
                    species_info_obj,
                    options,
                    prog_caller,
                    species_xml=species_xml,
                )

        elif options.q_start_from_fasta:
            # 3.
            species_info_obj = None
            species_info_obj = fasta_processor.processes_new_fasta(fasta_dir, options.dna)
            file_handler.handler.log_species()
            options = process_args.check_options(options, species_info_obj.species_to_use)
            # 4
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )
            # 5.
            species_xml = orthoxml_species(options, species_info_obj)
            # 6.
            printing.print_underline("Dividing up work for BLAST for parallel processing")
            runner.create_search_databases(species_info_obj, options, prog_caller)
            # 7.
            current_step = "Sequence search"
            log.step(current_step, "Started")
            runner.run_search(options, species_info_obj, seqs_info, prog_caller)
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 8.
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            infer_orthogroups.do_orthogroups(
                options, species_info_obj, seqs_info, species_names_dict, species_xml
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9.4
            if options.fix_files and not options.q_stop_after_mcl_groups:
                get_orthologs(
                    seqs_info,
                    species_names_dict,
                    species_info_obj,
                    options,
                    prog_caller,
                    species_xml=species_xml,
                )

        elif options.q_start_from_blast and not assign_resume:
            working_dirs = file_handler.handler.get_working_directory1_read()

            species_info_obj, _ = previous_run.process_previous_files(
                working_dirs,
                options.q_double_blast,
                check_blast=False,
            )

            # A results directory of an --assign run has the searches of the new
            # species against the orthogroup profiles (Blast<i>_-1.txt), not
            # the all-versus-all searches that -b continues from.
            if glob.glob(os.path.join(glob.escape(working_dirs[0]), "Blast*_-1.txt*")):
                file_handler.handler.log_fail_and_exit(
                    "ERROR: %s is the working directory of an '--assign' run, which cannot be "
                    "restarted with -b (it has no all-versus-all search results). Run --assign "
                    "again to complete it."
                    % working_dirs[0])
            # Results that are missing, or incomplete because their search was
            # interrupted, are created by running their commands from
            # blast_commands.txt, writing to this working directory.
            to_create = previous_run.get_blast_results_to_create(
                species_info_obj.species_to_use, options.q_double_blast
            )
            if to_create:
                current_step = "Sequence search"
                log.step(current_step, "Started")
                runner.recreate_search_results(to_create, working_dirs[0], options, prog_caller)
                log.step(current_step, "Completed")
                current_step = "Workflow"
                species_info_obj, _ = previous_run.process_previous_files(
                    working_dirs, options.q_double_blast, check_blast=True
                )
                still = previous_run.get_blast_results_to_create(
                    species_info_obj.species_to_use, options.q_double_blast,
                    only=[fn for fn, _ in to_create],     # the others were checked already
                )
                if still:
                    file_handler.handler.log_fail_and_exit(
                        "ERROR: These search results are still missing or incomplete after running "
                        "their commands:\n%s" % "\n".join(path or fn for fn, path in still)
                    )
            file_handler.handler.log_species()

            if not to_create:
                printer.print("All required search results are present and complete.")
            printer.print(
                "\nUsing the sequence search results in [dark_cyan]%s[/dark_cyan]"
                % escape(file_handler.handler.get_working_directory1_read()[0])
            )
            options = process_args.check_options(options, species_info_obj.species_to_use)
            # 4.
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )
            # 5.
            species_xml = orthoxml_species(options, species_info_obj)
            # 8
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            infer_orthogroups.do_orthogroups(
                options, species_info_obj, seqs_info, species_names_dict, species_xml
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9
            if options.fix_files and not options.q_stop_after_mcl_groups:
                get_orthologs(
                    seqs_info,
                    species_names_dict,
                    species_info_obj,
                    options,
                    prog_caller,
                    species_xml=species_xml,
                )

        elif options.q_start_from_groups:
            # 0.
            check_blast = not options.q_msa_trees
            species_info_obj, _ = previous_run.process_previous_files(
                file_handler.handler.get_working_directory1_read(),
                options.q_double_blast,
                check_blast=check_blast,
            )
            file_handler.handler.log_species()
            options = process_args.check_options(options, species_info_obj.species_to_use)

            ### 9
            species_xml = orthoxml_species(options, species_info_obj)
            ### 8
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )
            
            # infer_orthogroups.do_orthogroups(
            #     options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML
            # )

            get_orthologs(
                seqs_info,
                species_names_dict,
                species_info_obj,
                options,
                prog_caller,
                i_og_restart=i_og_restart_state,
                species_xml=species_xml,
            )


        elif options.q_start_from_trees:
            species_info_obj, _ = previous_run.process_previous_files(
                file_handler.handler.get_working_directory1_read(),
                options.q_double_blast,
                check_blast=False,
            )
            file_handler.handler.log_species()
            options = process_args.check_options(options, species_info_obj.species_to_use)
            # GetOrthologs_FromTrees(options)

            # orthologs.orthologs_from_trees(
            #     options.min_seq,
            #     options.recon_method,
            #     options.nBlast,
            #     options.nProcessAlg,
            #     options.speciesTreeFN,
            #     options.qAddSpeciesToIDs,
            #     options.qSplitParaClades,
            #     options.fewer_open_files,
            #     exist_msa=options.qMSATrees,
            #     fix_files=options.fix_files,
            #     mclInflation=options.mclInflation
            # )

            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )

            workflow.orthologs_from_gene_trees(
                seqs_info,
                species_names_dict,
                species_info_obj,
                options,
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
                options.recon_method,
                options.n_blast,
                options.n_process_alg,
                options.q_add_species_to_ids,
                options.fewer_open_files,
                options.species_tree_fn,
                options.q_stop_after_seqs,
                options.q_stop_after_alignments,
                options.q_stop_after_trees,
                options.q_msa_trees,
                options.q_phyldog,
                options.name,
                options.q_split_para_clades,
                save_space=options.save_space,
                root_from_previous=False,
                i_og_restart=i_og_restart_state,
                species_xml=orthoxml_species(options, species_info_obj),
            )
        elif options.q_start_from_species_trees:
            species_info_obj, _ = previous_run.process_previous_files(
                file_handler.handler.get_working_directory1_read(),
                options.q_double_blast,
                check_blast=False,
            )
            file_handler.handler.log_species()
            options = process_args.check_options(options, species_info_obj.species_to_use)
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )

            workflow.orthologs_from_gene_species_trees(
                seqs_info,
                species_names_dict,
                species_info_obj,
                options,
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
                options.recon_method,
                options.n_blast,
                options.n_process_alg,
                options.q_add_species_to_ids,
                options.species_tree_fn,
                options.fewer_open_files,  # Open one ortholog file per species when analysing trees
                q_split_para_clades=options.q_split_para_clades,
                i_og_restart=i_og_restart_state,
                species_xml=orthoxml_species(options, species_info_obj),
            )

        elif options.q_fast_add or assign_resume:
            if assign_resume:
                # Continuing an --assign analysis stopped or interrupted
                # (-op, or during its searches): the core analysis and the new
                # species are those recorded (run_state). Its working directory
                # (the first base directory) already holds the new species.
                assign_state = run_state.assign_info(results_dir, state)
                state_dirs = run_state.dirs(results_dir, state)
                continuation_dir = state_dirs["core_results"]
                fasta_dir = assign_state["new_species"]
                first_new_species = assign_state["first_new_species"]
                printer.print("\nContinuing the --assign analysis of the new species in [dark_cyan]%s"
                              % escape(fasta_dir))
                printer.print("to the core analysis in [dark_cyan]%s" % escape(continuation_dir))
                # whether its profile search finished, before this run logs its own
                profile_search_done = not file_io.search_checkpoint(
                    file_handler.handler.get_working_directory1_read()[0], step="Search orthogroup profiles")[0]
                # The core species, as the run had them when it started (its
                # working directory had the core's SpeciesIDs.txt, the new
                # species were added after the profiles were made)
                species_info_obj = species_info.SpeciesInfo()
                species_info_obj.species_to_use = list(range(first_new_species))
                species_info_obj.n_sp_all = first_new_species
                wd_list = file_handler.handler.get_working_directory1_read()[1:]   # the core's
            else:
                # Prepare previous directory as database
                species_info_obj, species_to_use_names = previous_run.process_previous_files(
                    file_handler.handler.get_working_directory1_read(),
                    options.q_double_blast,
                    check_blast=False,
                )
                profile_search_done = False
                wd_list = file_handler.handler.get_working_directory1_read()
            # Check previous directory has been done with MSA trees
            if not assign.check_for_assign(continuation_dir, species_info_obj):
                util.fail()
            printing.print_underline("Creating orthogroup profiles")
            current_step = "Create orthogroup profiles"
            log.step(current_step, "Started")
            fn_diamond_db, q_hogs = assign.prepare_assign_database(
                options.min_seq,
                continuation_dir,
                wd_list,
                species_info_obj.n_sp_all,
                species_info_obj,
                options,
                prog_caller,
                tree_program=options.tree_program,
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # print(
            #     "\nAdding new species in %s to existing analysis in %s"
            #     % (fastaDir, continuationDir)
            # )
            if assign_resume:
                # the new species were added by the run being continued
                species_info_obj, _ = previous_run.process_previous_files(
                    file_handler.handler.get_working_directory1_read(),
                    options.q_double_blast,
                    check_blast=False,
                )
                species_info_obj.i_first_new_species = first_new_species
            else:
                printer.print(f"\nAdding new species in [dark_cyan]{fasta_dir}")
                printer.print(f"to existing analysis in [dark_cyan]{continuation_dir}")

                species_info_obj = fasta_processor.processes_new_fasta(
                    fasta_dir, options.dna, species_info_obj, species_to_use_names
                )
                # the core analysis and the new species, for a restart of this one
                run_state.record_assign(results_dir, core_results=continuation_dir, new_species=fasta_dir,
                                        first_new_species=species_info_obj.i_first_new_species)
                run_state.record_dirs(results_dir, core_working=file_handler.handler.get_working_directory1_read()[1])

            if options.search_program in ["mmseqs"]:
                print(f"Create {options.search_program} new species database")
                runner.create_search_databases(
                    species_info_obj, options, prog_caller, new_species=True
                )

            options = process_args.check_options(options, species_info_obj.species_to_use)
            seqs_info = species_info.get_seqs_info(
                file_handler.handler.get_working_directory1_read(),
                species_info_obj.species_to_use,
                species_info_obj.n_sp_all,
            )
            species_xml = orthoxml_species(options, species_info_obj)
            # Add genes to orthogroups
            current_step = "Search orthogroup profiles"
            log.step(current_step, "Started")
            results_files, results_layout = runner.run_search_assign(
                options, species_info_obj, fn_diamond_db, prog_caller,
                reuse_complete=profile_search_done,
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # Clade-specific genes
            species_names_dict = previous_run.species_name_dict(
                file_handler.handler.get_species_ids_fn()
            )
            # if orphan_genes_version == 1:
            #     # v1 - This is unsuitable, it does an all-v-all search of all unassigned genes. Although these should have
            #     # been depleted of all genes that are not clade-specific, the resulting search still takes too long.
            #     # clade_specific_orthogroups_v1 function is now inside the src/orthofinder/legacy/utils/clade_specific_orthogroups.py
            #     clustersFilename_pairs, i_og_restart = clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)
            #     raise Exception("If q_hjogs then should be reading the N0.tsv file, not the original clusters")
            #     infer_orthogroups.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML=None)
            if orphan_genes_version == 2:
                # v2 - Infer rooted species tree from new rooted gene trees, identify new species-clades & search within these
                current_step = "Infer clade-specific orthogroups"
                log.step(current_step, "Started")
                clusters_filename_pairs, i_og_restart = between_core_orthogroups_workflow(
                    continuation_dir,
                    species_info_obj,
                    seqs_info,
                    options,
                    prog_caller,
                    species_names_dict,
                    results_files,
                    q_hogs,
                    results_layout=results_layout,
                    only_missing=assign_resume,
                )
                log.step(current_step, "Completed")
                current_step = "Workflow"

                # Infer clade-specific orthogroup gene trees
                infer_orthogroups.post_clustering_orthogroups(
                    clusters_filename_pairs,
                    species_info_obj,
                    seqs_info,
                    species_names_dict,
                    options,
                    species_xml=species_xml,
                )
                if options.species_tree_fn is None:
                    # No user species tree, use the one we've just inferred
                    options.species_tree_fn = file_handler.handler.get_species_tree_results_fn(
                        None, True
                    )
                # kept for a restart (-fg, -fgt, -fst) after a stop (-og, -ogt, -ost)
                file_handler.handler.log_species_tree_used(run_state.keep_species_tree(results_dir, options.species_tree_fn))
                run_state.update(results_dir, "assign", i_og_restart=i_og_restart)
                file_handler.handler.log_assign_first_orthogroup_to_infer(i_og_restart)
            if options.fix_files and not options.q_stop_after_mcl_groups:
                get_orthologs(
                    seqs_info,
                    species_names_dict,
                    species_info_obj,
                    options,
                    prog_caller,
                    i_og_restart,
                    species_xml=species_xml,
                )
        else:
            raise NotImplementedError
            # ptm = parallel_task_manager.ParallelTaskManagerSingleton()
            ptm.stop()
        current_step = "Process output files"
        log.step(current_step, "Started")
        if not options.save_space:
            # split up the orthologs into one file per species-pair
            ortholog_writers.split_ortholog_files(file_handler.handler.get_orthologs_directory())

        ### ------------- Compress the Gene_Trees --------------
        gene_tree_dir = file_handler.handler.get_ogs_tree_dir(q_results=True)
        # usr_gene_tree_fn = files.FileHandler.GetUserTreeFN()
        # paths.compress_files(gene_tree_dir, usr_gene_tree_fn)

        ### ------------- Comprees the Resolved_Gene_Trees ------------
        resolved_gene_tree_dir = file_handler.handler.get_ogs_recon_tree_dir(q_results=True)
        usr_resolved_gene_tree_fn = file_handler.handler.get_user_resolved_tree_fn()
        paths.compress_files(resolved_gene_tree_dir, usr_resolved_gene_tree_fn)

        ### ---------- Clean up WorkingDirectory ---------------
        if options.rm_gene_trees:
            paths.cleanup_path(gene_tree_dir)

        if options.rm_resolved_gene_trees:
            paths.cleanup_path(resolved_gene_tree_dir)

        d_results = (
            os.path.normpath(file_handler.handler.get_results_directory1()) + os.path.sep
        )

        if options.fewer_open_files and options.save_space:
            for i in range(len(species_info_obj.species_to_use)):
                sp0 = species_info_obj.species_to_use[i]
                sp0_name = species_names_dict[sp0]
                sp_path = os.path.join(
                    file_handler.handler.get_orthologs_directory(), f"{sp0_name}.tsv"
                )
                # the uncompressed file left from writing the .tsv.gz; only a
                # regular file is removed, never e.g. a folder a user extracted to
                if os.path.isfile(sp_path):
                    os.remove(sp_path)

        if options.orthoxml:
            # -x: the HOGs and orthologs as OrthoXML and web pages, from the final
            # results files (the directories are not made if a run removed them)
            orthogroups_orthoxml = file_handler.handler.get_orthogroup_results_fn_base() + ".orthoxml"
            if os.path.exists(orthogroups_orthoxml):
                try:
                    orthoxml.write_hogs_and_orthologs(
                        orthogroups_orthoxml,
                        os.path.join(d_results, "Phylogenetic_Hierarchical_Orthogroups"),
                        os.path.join(d_results, "Orthologs"),
                    )
                except (OSError, KeyError, ValueError) as e:
                    print("WARNING: The HOGs and orthologs could not be written in OrthoXML format: "
                          + escape(repr(e)))

        # printer.print("\nResults:\n    %s" % d_results, style="path")
        if options.q_stop_after_mcl_groups:      # stopped after the orthogroups (-og)
            printing.print_run_end(d_results, "-fg", "its orthogroups")
        else:
            printing.print_run_end(d_results)
        file_handler.handler.write_to_log("OrthoFinder run completed\n", True)
        log.step(current_step, "Completed")
        current_step = "Workflow"
        log.info("OrthoFinder run completed in %.2f seconds", time.perf_counter() - start)

    except Exception as e:
        current_step = run_logging.RunLogger.failed_stage(current_step)
        exit_code = 1
        if log is not None:
            log.log("ERROR: Step failed",
                    step=current_step, level="ERROR", exc_info=True)
            for stream in ("stdout", "stderr"):
                output = getattr(e, stream, None)
                if output:
                    if isinstance(output, bytes):
                        output = output.decode("utf-8", errors="replace")
                    log.log("%s:\n%s", stream, output,
                            step=current_step, level="ERROR")
                    printer.console.print(f"{stream}:\n{output}", markup=False, highlight=False)
        import traceback
        traceback.print_exception(type(e), e, e.__traceback__)
        # ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()
        sys.exit(1)

    except KeyboardInterrupt:
        current_step = run_logging.RunLogger.failed_stage(current_step)
        exit_code = 1
        if log is not None:
            log.step(current_step, "Interrupted by user", level="WARNING")
        printer.print("\nProgram terminated by user.", style="error")
        sys.exit(1)

    except SystemExit as e:
        current_step = run_logging.RunLogger.failed_stage(current_step)
        exit_code = e.code
        if log is not None:
            if e.code is None or e.code == 0:
                log.step(current_step, "Workflow stopped normally before the full analysis finished")
            else:
                log.log("ERROR: Workflow exited with status %s", e.code,
                        step=current_step, level="ERROR", exc_info=True)
        raise

    finally:
        run_logging.RunLogger.set_logger(None)
        if log is not None:
            log.close()
        # ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()
        end = time.perf_counter()
        time_elapsed = end - start
        print()

        if len(input_args) == 0 or input_args[0] in [
            "--help",
            "-h",
            "-v",
            "--version",
            "-sm",
            "--scoring-matrix",
        ]:
            sys.exit()

        # printer.print(f"OrthoFinder finished in {time_elapsed:5f}s", end="\n" * 2, style="info")
        printer.print(
            f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] finished in ", end=""
        )
        printer.print(f"[green]{time_elapsed:5f}[/green]s", end="\n" * 2)
        file_handler.handler.reset()
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
