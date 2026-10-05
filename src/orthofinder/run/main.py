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
from ..utils import (
    parallel_task_manager,
    files,
    util,
    program_caller,
    fasta_processor,
    file_io,
    run_state,
)
from ..orthogroups import gathering, orthogroups_set
from ..orthogroups import accelerate as acc
from orthofinder.utils import logging as run_logging


from ..tools import astral, mcl, tree
from ..gene_tree_inference import trees2ologs_of, infer_trees, tree_processor
from . import process_args, check_dependencies, run_commands, species_info
from .. import orphan_genes_version, __version__, __location__
from ..comparative_genomics import orthologues
from ..utils.util import printer
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


def GetProgramCaller():
    config_file = os.path.join(configfile_location, "config.json")
    pc = program_caller.ProgramCaller(
        config_file if os.path.exists(config_file) else None
    )
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        pc.Add(pc_user)
    return pc


"""
Main
-------------------------------------------------------------------------------
"""


# 9
def GetOrthologues(
    seqsInfo,
    speciesNamesDict,
    speciesInfoObj,
    options,
    prog_caller,
    i_og_restart=0,
    speciesXML=None,
):
    util.PrintUnderline("Analysing Orthogroups", True)
    orthologues.OrthologuesWorkflow(
        seqsInfo,
        speciesNamesDict,
        speciesInfoObj,
        options,
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        prog_caller,
        options.msa_program,
        options.tree_program,
        options.recon_method,
        options.nBlast,
        options.nProcessAlg,
        options.qDoubleBlast,
        options.qAddSpeciesToIDs,
        options.qTrim,
        options.fewer_open_files,
        options.cmd_order,
        options.method_threads,
        options.method_threads_large,
        options.method_threads_small,
        options.threshold,
        options.speciesTreeFN,
        options.qStopAfterSeqs,
        options.qStopAfterAlignments,
        options.qStopAfterTrees,
        options.qStopAfterSpeciesTrees,
        options.qMSATrees,
        options.qPhyldog,
        options.name,
        options.qSplitParaClades,
        save_space=options.save_space,
        root_from_previous=False,
        i_og_restart=i_og_restart,
        speciesXML=speciesXML,
    )
    # util.PrintTime("Done writing files")


def WarnSpeciesWithoutAssignedGenes(ogs_new_species, speciesInfoObj, speciesNamesDict):
    """
    --assign: warn about new species none of whose genes were assigned to a
    core orthogroup (no hits in the orthogroup-profile search). Such a species
    is in no gene tree, so it cannot be placed in the species tree, and its
    genes stay unassigned: the run would otherwise end as if it had worked.
    """
    with_genes = {int(g.split("_", 1)[0]) for genes in ogs_new_species.values() for g in genes}
    missing = [iSp for iSp in range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll)
               if iSp not in with_genes]
    if missing:
        printer.print(
            "\nWARNING: No genes of these new species were assigned to the core orthogroups (no hits "
            "in the search against the orthogroup profiles): %s. They cannot be placed in the species "
            "tree, and their genes will be unassigned. Check that they are related to the core species "
            "and that the search program can compare these sequences (e.g. for nucleotide sequences, "
            "-d with -S blast or -S mmseqs).\n"
            % escape(", ".join(str(speciesNamesDict.get(iSp, iSp)) for iSp in missing)),
            style="warning")


def BetweenCoreOrthogroupsWorkflow(
    continuationDir,
    speciesInfoObj,
    seqsInfo,
    options,
    prog_caller,
    speciesNamesDict,
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
        ogs = acc.read_hogs(continuationDir, "N0")
    else:
        ogs = acc.get_original_orthogroups()
    i_og_restart = 0
    ogs_new_species, _ = acc.assign_genes(results_files, options.nProcessAlg, layout=results_layout)
    WarnSpeciesWithoutAssignedGenes(ogs_new_species, speciesInfoObj, speciesNamesDict)
    clustersFilename_pairs = acc.write_all_orthogroups(
        ogs, ogs_new_species, []
    )  # this updates ogs

    ogSet = orthogroups_set.OrthoGroupsSet(
        options.min_seq,
        files.FileHandler.GetWorkingDirectory1_Read(),
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        options.tree_program,
        idExtractor=util.FirstWordExtractor,
        mclInflation=options.mclInflation,
    )

    if options.qStopAfterGroups and options.speciesTreeFN is None:
        # Can't infer clade-specific groups
        print("\nSpecies tree required for clade-speicfic orthogroups - skipping")
        return clustersFilename_pairs, i_og_restart

    n_unassigned = acc.write_unassigned_fasta(ogs, None, speciesInfoObj)

    # Get/Infer species tree
    if options.speciesTreeFN is None:
        # Infer gene trees
        # We write orthogroup & stats results files in the following code, which we should avoid & only do once all OGs are done.
        gathering.post_clustering_orthogroups(
            clustersFilename_pairs,
            speciesInfoObj,
            seqsInfo,
            speciesNamesDict,
            options,
            speciesXML=None,
            q_incremental=True,
        )

        infer_trees.InferGeneAndSpeciesTrees(
            ogSet,
            prog_caller,
            options.msa_program,
            options.tree_program,
            options.nBlast,
            options.nProcessAlg,
            options.qDoubleBlast,
            options.qAddSpeciesToIDs,
            options.qTrim,
            cmd_order=options.cmd_order,
            method_threads=options.method_threads,
            method_threads_large=options.method_threads_large,
            method_threads_small=options.method_threads_small,
            threshold=options.threshold,
            userSpeciesTree=None,
            qStopAfterSeqs=False,
            qStopAfterAlign=False,
            qMSA=options.qMSATrees,
            qPhyldog=False,
            results_name=options.name,
            root_from_previous=True,
            n_skip=options.n_skip,
        )

        # Infer species tree
        astral_fn = files.FileHandler.GetAstralFilename()
        astral.create_input_file(
            files.FileHandler.GetOGsTreeDir(), astral_fn, n_skip=options.n_skip
        )
        species_tree_unrooted_fn = files.FileHandler.GetSpeciesTreeUnrootedFN()
        parallel_task_manager.RunCommand(
            astral.get_astral_command(
                astral_fn, species_tree_unrooted_fn, options.nBlast
            ),
            # ASTRAL-Pro writes its progress log to stderr: only show it if it fails.
            qPrintStderr=False,
            raise_on_error=True,
        )

        # Root it
        core_rooted_species_tree = tree.Tree(
            files.FileHandler.GetCoreSpeciesTreeIDsRootedFN(), format=1
        )
        species_to_speices_map = lambda x: x
        rooted_species_tree_ids, qHaveSupport = tree_processor.CheckAndRootTree(
            species_tree_unrooted_fn, core_rooted_species_tree, species_to_speices_map
        )

        if rooted_species_tree_ids is None:
            print(
                "ERROR: Species tree inference failed. Please check for errors and check the species tree files: \n%s \n%s"
                % (species_tree_unrooted_fn, core_rooted_species_tree)
            )
            util.Fail()

        rooted_species_tree_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
        rooted_species_tree_ids.write(outfile=rooted_species_tree_fn)

        spTreeUnrootedFN = files.FileHandler.GetSpeciesTreeResultsFN(None, True)
        util.RenameTreeTaxa(
            rooted_species_tree_ids,
            spTreeUnrootedFN,
            ogSet.SpeciesDict(),
            qSupport=qHaveSupport,
            qFixNegatives=True,
        )

        labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
        util.RenameTreeTaxa(
            rooted_species_tree_ids,
            labeled_tree_fn,
            ogSet.SpeciesDict(),
            qSupport=False,
            qFixNegatives=True,
            label="N",
        )
        i_og_restart = len(ogs)  # Need to process the clade-specific orthogroups only
    else:
        util.PrintUnderline("Using user-supplied species tree")
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        infer_trees.ConvertUserSpeciesTree(
            options.speciesTreeFN, ogSet.SpeciesDict(), spTreeFN_ids
        )
        rooted_species_tree_fn = spTreeFN_ids

    # Identify clades for clade-specific orthogroup inference
    iSpeciesCore = set(speciesInfoObj.get_original_species())
    species_clades = acc.get_new_species_clades(rooted_species_tree_fn, iSpeciesCore)
    util.PrintUnderline(
        "Identifying clade-specific orthogroups for the following clades:", qHeavy=True
    )
    species_dict = ogSet.SpeciesDict()
    for i, clade in enumerate(species_clades):
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
    print("")

    # The clade-specific step costs k^2 searches (and more for clustering) per
    # clade of k species: report its size and warn about very large clades.
    clade_rows = acc.clade_costs(
        [list(map(str, clade)) for clade in species_clades],
        map(str, iSpeciesCore),
        n_genes={str(isp): n for isp, n in enumerate(n_unassigned)},
    )
    for warning in acc.report_clades(clade_rows, species_dict):
        run_logging.RunLogger.message(warning, level="WARNING")
    print("")

    if not species_clades:
        # Every new species falls in a part of the tree covered by the core:
        # there are no clade-specific orthogroups to infer.
        print("No new-species clades: skipping clade-specific orthogroup inference\n")
        clustersFilename_pairs, i_og_restart = acc.write_all_orthogroups(
            ogs, {}, [], restart_index=i_og_restart
        )
        return clustersFilename_pairs, i_og_restart

    # Clade-specific orthogroup inference
    run_commands.CreateSearchDatabases(
        speciesInfoObj, options, prog_caller, q_unassigned_genes=True
    )

    # provide list of clades, only run these searches (and only if the fasta files are non-empty)
    run_commands.RunSearch(
        options,
        speciesInfoObj,
        seqsInfo,
        prog_caller,
        n_genes_per_species=n_unassigned,
        species_clades=species_clades,
        only_missing=only_missing,
    )
    # process the results files - only if they are present and non-empty
    options.v2_scores = True
    n_clades = len(species_clades)
    clustersFilename_pairs_unassigned_all = []
    for i_clade, clade in enumerate(species_clades):
        util.PrintUnderline(
            "OrthoFinder clutering on new species clade %d of %d"
            % (i_clade + 1, n_clades)
        )
        print(str(i_clade) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
        speciesInfo_clade = copy.deepcopy(speciesInfoObj)
        speciesInfo_clade.speciesToUse = clade
        seqsInfo_clade = util.SeqsInfoRecompute(seqsInfo, clade)
        clustersFilename_pairs_unassigned = gathering.DoOrthogroups(
            options,
            speciesInfo_clade,
            seqsInfo_clade,
            speciesNamesDict,
            speciesXML=None,
            i_unassigned=i_clade,
        )
        clustersFilename_pairs_unassigned_all.append(clustersFilename_pairs_unassigned)
    ogs_clade_specific_list = [
        mcl.GetPredictedOGs(filename)
        for filename in clustersFilename_pairs_unassigned_all
    ]

    # OGs have had assigned genes added to them already
    # Single-gene orthogroups removed here shift the positions of later ones,
    # so the restart index is adjusted along with them.
    clustersFilename_pairs, i_og_restart = acc.write_all_orthogroups(
        ogs, {}, ogs_clade_specific_list, restart_index=i_og_restart
    )
    return clustersFilename_pairs, i_og_restart


# def GetOrthologues_FromTrees(options):
#     orthologues.OrthologuesFromTrees(
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
    files.FileHandler.reset()
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
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        # prog_caller = GetProgramCaller()
        (
            prog_caller,
            options,
            fastaDir,
            continuationDir,
            resultsDir_nonDefault,
            pickleDir_nonDefault,
            user_specified_M,
        ) = process_args.ProcessArgs(args)

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
        q_restart = ((options.qStartFromBlast and not options.qStartFromFasta) or options.qStartFromGroups
                     or options.qStartFromTrees or options.qStartFromSpeciesTrees)
        if q_restart and continuationDir is not None:
            _, state = run_state.find(continuationDir)
            if not state:
                # no run_state.json (older versions, or the file removed): what
                # Log.txt records, read before this run rewrites it
                state = files.StateFromLog(continuationDir)
            i_og_restart_state = state.get("assign", {}).get("i_og_restart", 0)
            changed = run_state.restore_settings(state, options, input_args)
            # the run's search program (restored above, or the same as the default)
            options.search_program_from_run = ("search_program" in state.get("settings", {})
                                               and not options.search_program_given)

        # "-S blast" is blastn or blastp (-d). Before the file handler is set
        # up, as it logs the search program.
        options = process_args.ResolveSearchProgram(options, prog_caller)

        files.InitialiseFileHandler(
            options,
            fastaDir,
            continuationDir,
            resultsDir_nonDefault,
            pickleDir_nonDefault,
        )


        # A restart keeps the record of the stages it continues from, not of
        # those it runs again (so running it again does not record them twice).
        orthogroup_stages = ("Infer orthogroups", "Infer clade-specific orthogroups")
        if options.qStartFromBlast and not options.qStartFromFasta:
            file_io.trim_checkpoint(files.FileHandler.GetCheckPointFN(), redo=orthogroup_stages)
        elif options.qStartFromGroups:
            file_io.trim_checkpoint(files.FileHandler.GetCheckPointFN(), completed=[orthogroup_stages])
        elif options.qStartFromTrees or options.qStartFromSpeciesTrees:
            file_io.trim_checkpoint(files.FileHandler.GetCheckPointFN(),
                                    completed=[orthogroup_stages, ("Build alignments and gene trees",)])
        log = run_logging.Logger(
            files.FileHandler.GetCheckPointFN(),
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
            options.nBlast,
            options.nProcessAlg,
            __version__,
            files.FileHandler.GetResultsDirectory1(),
        )


        # print("Results directory: %s" % files.FileHandler.GetResultsDirectory1())
        # printer.print("Results directory:", style="path")
        # printer.print(f"    [dark_cyan]{files.FileHandler.GetResultsDirectory1()}")

        printer.print("[bold]Results directory:[/bold]")
        printer.print(f"    [dark_cyan]{files.FileHandler.GetResultsDirectory1()}")

        # A restart also continues with the species tree the run used (unless
        # -s is given), and its point for restarting gene-tree inference
        # (--assign); its settings were restored above.
        results_dir = files.FileHandler.GetResultsDirectory1()
        if q_restart:
            if changed:
                printer.print("Continuing with the settings of the run being continued: %s"
                              % escape(", ".join("%s=%s" % c for c in changed)))
                files.FileHandler.WriteToLog("Settings of the run being continued: %s\n"
                                             % ", ".join("%s=%s" % c for c in changed))
            kept_tree = run_state.file(results_dir, state, "species_tree")
            if i_og_restart_state:
                # kept in this run's Log.txt too (a restart rewrites it)
                files.FileHandler.LogAssignFirstOrthogroupToInfer(i_og_restart_state)
            if options.speciesTreeFN is None and kept_tree is not None:
                options.speciesTreeFN = kept_tree
                printer.print("Using the species tree of the run being continued:")
                printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(kept_tree))
        # -b on the results of an --assign analysis continues that analysis
        assign_resume = (options.qStartFromBlast and not options.qStartFromFasta
                         and state.get("run", {}).get("analysis") == "assign")
        if assign_resume and not ("first_new_species" in state.get("assign", {})
                                  and "new_species" in state.get("assign", {})
                                  and "core_results" in state.get("dirs", {})):
            # stopped before the new species were added (e.g. while creating
            # the orthogroup profiles): there is nothing to continue from
            files.FileHandler.LogFailAndExit(
                "ERROR: The '--assign' run in %s stopped before its new species were added, so it "
                "cannot be continued with -b. Run --assign again." % results_dir)
        run_state.record_run(results_dir, "assign" if options.qFastAdd else "core",
                             " ".join(["orthofinder"] + input_args))
        run_state.save_settings(results_dir, options)
        if options.speciesTreeFN is not None:
            files.FileHandler.LogSpeciesTreeUsed(run_state.keep_species_tree(results_dir, options.speciesTreeFN))

        check_dependencies.CheckDependencies(
            options,
            user_specified_M,
            prog_caller,
            files.FileHandler.GetWorkingDirectory1_Read()[0],
            q_assign=assign_resume,
        )
        if continuationDir is not None:
            util.PrintRestartInfo(options, continuationDir)


        # if using previous Trees etc., check these are all present - Job for orthologues
        if options.qStartFromBlast and options.qStartFromFasta:
            # 0. Check Files
            speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast
            )
            # print(
            #     "\nAdding new species in %s to existing analysis in %s"
            #     % (fastaDir, continuationDir)
            # )
            printer.print(f"\nAdding new species in [dark_cyan]{fastaDir}")
            printer.print(f"to existing analysis in [dark_cyan]{continuationDir}")
            # 3.
            speciesInfoObj = fasta_processor.ProcessesNewFasta(
                fastaDir, options.dna, speciesInfoObj, speciesToUse_names
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            # 5.
            # speciesXML = (
            #     species_info.GetXMLSpeciesInfo(speciesInfoObj, options)
            #     if options.speciesXMLInfoFN
            #     else None
            # )
            speciesXML = None
            # 6.
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7.
            current_step = "Sequence search"
            log.step(current_step, "Started")
            run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 8.
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            gathering.DoOrthogroups(
                options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9.
            if options.fix_files and not options.qStopAfterMCLGroups:
                GetOrthologues(
                    seqsInfo,
                    speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromFasta:
            # 3.
            speciesInfoObj = None
            speciesInfoObj = fasta_processor.ProcessesNewFasta(fastaDir, options.dna)
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            # 5.
            # speciesXML = (
            #     species_info.GetXMLSpeciesInfo(speciesInfoObj, options)
            #     if options.speciesXMLInfoFN
            #     else None
            # )
            speciesXML = None
            # 6.
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            # 7.
            current_step = "Sequence search"
            log.step(current_step, "Started")
            run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 8.
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            gathering.DoOrthogroups(
                options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9.4
            if options.fix_files and not options.qStopAfterMCLGroups:
                GetOrthologues(
                    seqsInfo,
                    speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromBlast and not assign_resume:
            working_dirs = files.FileHandler.GetWorkingDirectory1_Read()

            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                working_dirs,
                options.qDoubleBlast,
                check_blast=False,
            )

            # A results directory of an --assign run has the searches of the new
            # species against the orthogroup profiles (Blast<i>_-1.txt), not
            # the all-versus-all searches that -b continues from.
            if glob.glob(os.path.join(glob.escape(working_dirs[0]), "Blast*_-1.txt*")):
                files.FileHandler.LogFailAndExit(
                    "ERROR: %s is the working directory of an '--assign' run, which cannot be "
                    "restarted with -b (it has no all-versus-all search results). Run --assign "
                    "again to complete it."
                    % working_dirs[0])
            # Results that are missing, or incomplete because their search was
            # interrupted, are created by running their commands from
            # blast_commands.txt, writing to this working directory.
            to_create = species_info.GetBlastResultsToCreate(
                speciesInfoObj.speciesToUse, options.qDoubleBlast
            )
            if to_create:
                current_step = "Sequence search"
                log.step(current_step, "Started")
                run_commands.recreate_search_results(to_create, working_dirs[0], options, prog_caller)
                log.step(current_step, "Completed")
                current_step = "Workflow"
                speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                    working_dirs, options.qDoubleBlast, check_blast=True
                )
                still = species_info.GetBlastResultsToCreate(
                    speciesInfoObj.speciesToUse, options.qDoubleBlast,
                    only=[fn for fn, _ in to_create],     # the others were checked already
                )
                if still:
                    files.FileHandler.LogFailAndExit(
                        "ERROR: These search results are still missing or incomplete after running "
                        "their commands:\n%s" % "\n".join(path or fn for fn, path in still)
                    )
            files.FileHandler.LogSpecies()

            if not to_create:
                printer.print("All required search results are present and complete.")
            printer.print(
                "\nUsing the sequence search results in [dark_cyan]%s[/dark_cyan]"
                % escape(files.FileHandler.GetWorkingDirectory1_Read()[0])
            )
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            # 4.
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            # 5.
            # speciesXML = (
            #     species_info.GetXMLSpeciesInfo(speciesInfoObj, options)
            #     if options.speciesXMLInfoFN
            #     else None
            # )
            speciesXML = None
            # 8
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            current_step = "Infer orthogroups"
            log.step(current_step, "Started")
            gathering.DoOrthogroups(
                options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # 9
            if options.fix_files and not options.qStopAfterMCLGroups:
                GetOrthologues(
                    seqsInfo,
                    speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromGroups:
            # 0.
            check_blast = not options.qMSATrees
            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(),
                options.qDoubleBlast,
                check_blast=check_blast,
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)

            ### 9
            # speciesXML = (
            #     species_info.GetXMLSpeciesInfo(speciesInfoObj, options)
            #     if options.speciesXMLInfoFN
            #     else None
            # )
            speciesXML = None
            ### 8
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            
            # gathering.DoOrthogroups(
            #     options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML
            # )

            GetOrthologues(
                seqsInfo,
                speciesNamesDict,
                speciesInfoObj,
                options,
                prog_caller,
                i_og_restart=i_og_restart_state,
                speciesXML=speciesXML,
            )


        elif options.qStartFromTrees:
            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(),
                options.qDoubleBlast,
                check_blast=False,
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            # GetOrthologues_FromTrees(options)

            # orthologues.OrthologuesFromTrees(
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

            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )

            orthologues.OrthologuesFromGeneTrees(
                seqsInfo,
                speciesNamesDict,
                speciesInfoObj,
                options,
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
                options.recon_method,
                options.nBlast,
                options.nProcessAlg,
                options.qAddSpeciesToIDs,
                options.fewer_open_files,
                options.speciesTreeFN,
                options.qStopAfterSeqs,
                options.qStopAfterAlignments,
                options.qStopAfterTrees,
                options.qMSATrees,
                options.qPhyldog,
                options.name,
                options.qSplitParaClades,
                save_space=options.save_space,
                root_from_previous=False,
                i_og_restart=i_og_restart_state,
            )
        elif options.qStartFromSpeciesTrees:
            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(),
                options.qDoubleBlast,
                check_blast=False,
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )

            orthologues.OrthologuesFromGeneSpeciesTrees(
                seqsInfo,
                speciesNamesDict,
                speciesInfoObj,
                options,
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
                options.recon_method,
                options.nBlast,
                options.nProcessAlg,
                options.qAddSpeciesToIDs,
                options.speciesTreeFN,
                options.fewer_open_files,  # Open one ortholog file per species when analysing trees
                q_split_para_clades=options.qSplitParaClades,
                i_og_restart=i_og_restart_state,
                speciesXML=None,
            )

        elif options.qFastAdd or assign_resume:
            if assign_resume:
                # Continuing an --assign analysis stopped or interrupted
                # (-op, or during its searches): the core analysis and the new
                # species are those recorded (run_state). Its working directory
                # (the first base directory) already holds the new species.
                assign = run_state.assign_info(results_dir, state)
                state_dirs = run_state.dirs(results_dir, state)
                continuationDir = state_dirs["core_results"]
                fastaDir = assign["new_species"]
                first_new_species = assign["first_new_species"]
                printer.print("\nContinuing the --assign analysis of the new species in [dark_cyan]%s"
                              % escape(fastaDir))
                printer.print("to the core analysis in [dark_cyan]%s" % escape(continuationDir))
                # whether its profile search finished, before this run logs its own
                profile_search_done = not file_io.search_checkpoint(
                    files.FileHandler.GetWorkingDirectory1_Read()[0], step="Search orthogroup profiles")[0]
                # The core species, as the run had them when it started (its
                # working directory had the core's SpeciesIDs.txt, the new
                # species were added after the profiles were made)
                speciesInfoObj = util.SpeciesInfo()
                speciesInfoObj.speciesToUse = list(range(first_new_species))
                speciesInfoObj.nSpAll = first_new_species
                wd_list = files.FileHandler.GetWorkingDirectory1_Read()[1:]   # the core's
            else:
                # Prepare previous directory as database
                speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(
                    files.FileHandler.GetWorkingDirectory1_Read(),
                    options.qDoubleBlast,
                    check_blast=False,
                )
                profile_search_done = False
                wd_list = files.FileHandler.GetWorkingDirectory1_Read()
            # Check previous directory has been done with MSA trees
            if not acc.check_for_orthoxcelerate(continuationDir, speciesInfoObj):
                util.Fail()
            util.PrintUnderline("Creating orthogroup profiles")
            current_step = "Create orthogroup profiles"
            log.step(current_step, "Started")
            fn_diamond_db, q_hogs = acc.prepare_accelerate_database(
                options.min_seq,
                continuationDir,
                wd_list,
                speciesInfoObj.nSpAll,
                speciesInfoObj,
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
                speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                    files.FileHandler.GetWorkingDirectory1_Read(),
                    options.qDoubleBlast,
                    check_blast=False,
                )
                speciesInfoObj.iFirstNewSpecies = first_new_species
            else:
                printer.print(f"\nAdding new species in [dark_cyan]{fastaDir}")
                printer.print(f"to existing analysis in [dark_cyan]{continuationDir}")

                speciesInfoObj = fasta_processor.ProcessesNewFasta(
                    fastaDir, options.dna, speciesInfoObj, speciesToUse_names
                )
                # the core analysis and the new species, for a restart of this one
                run_state.record_assign(results_dir, core_results=continuationDir, new_species=fastaDir,
                                        first_new_species=speciesInfoObj.iFirstNewSpecies)
                run_state.record_dirs(results_dir, core_working=files.FileHandler.GetWorkingDirectory1_Read()[1])

            if options.search_program in ["mmseqs"]:
                print(f"Create {options.search_program} new species database")
                run_commands.CreateSearchDatabases(
                    speciesInfoObj, options, prog_caller, new_species=True
                )

            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            # Add genes to orthogroups
            current_step = "Search orthogroup profiles"
            log.step(current_step, "Started")
            results_files, results_layout = run_commands.RunSearch_accelerate(
                options, speciesInfoObj, fn_diamond_db, prog_caller,
                reuse_complete=profile_search_done,
            )
            log.step(current_step, "Completed")
            current_step = "Workflow"
            # Clade-specific genes
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            # if orphan_genes_version == 1:
            #     # v1 - This is unsuitable, it does an all-v-all search of all unassigned genes. Although these should have
            #     # been depleted of all genes that are not clade-specific, the resulting search still takes too long.
            #     # clade_specific_orthogroups_v1 function is now inside the src/orthofinder/legacy/utils/clade_specific_orthogroups.py
            #     clustersFilename_pairs, i_og_restart = clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)
            #     raise Exception("If q_hjogs then should be reading the N0.tsv file, not the original clusters")
            #     gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML=None)
            if orphan_genes_version == 2:
                # v2 - Infer rooted species tree from new rooted gene trees, identify new species-clades & search within these
                current_step = "Infer clade-specific orthogroups"
                log.step(current_step, "Started")
                clustersFilename_pairs, i_og_restart = BetweenCoreOrthogroupsWorkflow(
                    continuationDir,
                    speciesInfoObj,
                    seqsInfo,
                    options,
                    prog_caller,
                    speciesNamesDict,
                    results_files,
                    q_hogs,
                    results_layout=results_layout,
                    only_missing=assign_resume,
                )
                log.step(current_step, "Completed")
                current_step = "Workflow"

                # Infer clade-specific orthogroup gene trees
                gathering.post_clustering_orthogroups(
                    clustersFilename_pairs,
                    speciesInfoObj,
                    seqsInfo,
                    speciesNamesDict,
                    options,
                    speciesXML=None,
                )
                if options.speciesTreeFN is None:
                    # No user species tree, use the one we've just inferred
                    options.speciesTreeFN = files.FileHandler.GetSpeciesTreeResultsFN(
                        None, True
                    )
                # kept for a restart (-fg, -fgt, -fst) after a stop (-og, -ogt, -ost)
                files.FileHandler.LogSpeciesTreeUsed(run_state.keep_species_tree(results_dir, options.speciesTreeFN))
                run_state.update(results_dir, "assign", i_og_restart=i_og_restart)
                files.FileHandler.LogAssignFirstOrthogroupToInfer(i_og_restart)
            if options.fix_files and not options.qStopAfterMCLGroups:
                GetOrthologues(
                    seqsInfo,
                    speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    i_og_restart,
                    speciesXML=None,
                )
        else:
            raise NotImplementedError
            # ptm = parallel_task_manager.ParallelTaskManager_singleton()
            ptm.Stop()
        current_step = "Process output files"
        log.step(current_step, "Started")
        if not options.save_space:
            # split up the orthologs into one file per species-pair
            util.split_ortholog_files(files.FileHandler.GetOrthologuesDirectory())

        ### ------------- Compress the Gene_Trees --------------
        gene_tree_dir = files.FileHandler.GetOGsTreeDir(qResults=True)
        # usr_gene_tree_fn = files.FileHandler.GetUserTreeFN()
        # util.compress_files(gene_tree_dir, usr_gene_tree_fn)

        ### ------------- Comprees the Resolved_Gene_Trees ------------
        resolved_gene_tree_dir = files.FileHandler.GetOGsReconTreeDir(qResults=True)
        usr_resolved_gene_tree_fn = files.FileHandler.GetUserResolvedTreeFN()
        util.compress_files(resolved_gene_tree_dir, usr_resolved_gene_tree_fn)

        ### ---------- Clean up WorkingDirectory ---------------
        if options.rm_gene_trees:
            util.cleanup_path(gene_tree_dir)

        if options.rm_resolved_gene_trees:
            util.cleanup_path(resolved_gene_tree_dir)

        d_results = (
            os.path.normpath(files.FileHandler.GetResultsDirectory1()) + os.path.sep
        )

        if options.fewer_open_files and options.save_space:
            for i in range(len(speciesInfoObj.speciesToUse)):
                sp0 = speciesInfoObj.speciesToUse[i]
                sp0_name = speciesNamesDict[sp0]
                sp_path = os.path.join(
                    files.FileHandler.GetOrthologuesDirectory(), f"{sp0_name}.tsv"
                )
                # the uncompressed file left from writing the .tsv.gz; only a
                # regular file is removed, never e.g. a folder a user extracted to
                if os.path.isfile(sp_path):
                    os.remove(sp_path)

        # printer.print("\nResults:\n    %s" % d_results, style="path")
        if options.qStopAfterMCLGroups:      # stopped after the orthogroups (-og)
            util.PrintRunEnd(d_results, "-fg", "its orthogroups")
        else:
            util.PrintRunEnd(d_results)
        files.FileHandler.WriteToLog("OrthoFinder run completed\n", True)
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
        # ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
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
        # ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
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
        files.FileHandler.reset()
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
