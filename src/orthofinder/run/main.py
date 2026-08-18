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

import multiprocessing as mp
import platform
import sys

if __name__ == "__main__":
    if platform.system() == "Darwin":
        mp.set_start_method("fork")

import os
import time
import copy
import csv
import os.path

from ..utils import (
    parallel_task_manager,
    files,
    util,
    program_caller,
    fasta_processor,
)
from ..orthogroups import gathering, orthogroups_set
from ..orthogroups import accelerate as acc
from ..tools import astral, mcl, tree
from ..gene_tree_inference import trees2ologs_of, infer_trees
from . import process_args, check_dependencies, run_commands, species_info
from .. import orphan_genes_version, __version__, __location__
from ..comparative_genomics import orthologues
from ..utils.util import printer

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


def GetOrthologues(
        seqsInfo, speciesNamesDict,
        speciesInfoObj, options,
        prog_caller,
        i_og_restart=0,
        speciesXML=None,
    ):
    util.PrintUnderline("Analysing Orthogroups", True)
    orthologues.OrthologuesWorkflow(
        seqsInfo, speciesNamesDict,
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
        options.old_version,
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


def _get_singleton_orthogroups(ogs, ogs_unassigned_lists, species_ids):
    """Return genuine unassigned genes as one-gene orthogroups.

    MCL intentionally drops singleton rows during the unassigned-gene workflow
    because some can be implicit graph nodes. Recover singletons from the
    authoritative Unassigned.Species*.fa files instead.
    """
    clustered_unassigned = set()
    for ogs_unassigned in ogs_unassigned_lists:
        for og in ogs_unassigned:
            clustered_unassigned.update(og)

    existing_genes = set()
    for og in ogs:
        existing_genes.update(og)

    all_unassigned = set()
    for i_sp in species_ids:
        fw = fasta_processor.FastaWriter(
            files.FileHandler.GetSpeciesUnassignedFastaFN(i_sp)
        )
        all_unassigned.update(fw.SeqLists.keys())

    singleton_genes = all_unassigned - clustered_unassigned - existing_genes
    return [{gene} for gene in sorted(singleton_genes)]


def _write_single_copy_orthogroups(clusters_filename_pairs, speciesInfoObj):
    """Write the standard Orthogroups_SingleCopyOrthologues.txt file.

    An orthogroup is single-copy when it contains exactly one gene from every
    species in the current analysis. This is purely an orthogroup property and
    does not require MSAs, gene trees, a species tree, or orthologue inference.
    """
    ogs = mcl.GetPredictedOGs(clusters_filename_pairs)
    species_to_use = set(speciesInfoObj.speciesToUse)
    n_species = len(species_to_use)
    results_base = files.FileHandler.GetOrthogroupResultsFNBase()
    output_fn = results_base + "_SingleCopyOrthologues.txt"

    with open(output_fn, "w") as outfile:
        for i_og, og in enumerate(ogs):
            if len(og) != n_species:
                continue
            try:
                og_species = [int(gene.split("_", 1)[0]) for gene in og]
            except (ValueError, IndexError):
                continue
            if len(set(og_species)) == n_species and set(og_species) == species_to_use:
                outfile.write("OG%07d\n" % i_og)

    return output_fn


def AssignOnlyOrthogroupsWorkflow(
        continuationDir,
        speciesInfoObj,
        seqsInfo,
        options,
        prog_caller,
        speciesNamesDict,
        results_files,
        q_hogs,
    ):
    """Complete --assign orthogroup inference and stop before MSA/tree work.

    The new genes are first assigned to the existing/core orthogroups using the
    profile database. Remaining genes are clustered with MCL.

    If -s is supplied, the user tree is converted to internal species IDs and
    OrthoFinder's existing get_new_species_clades() logic is used, matching the
    clade-specific stage of the normal --assign workflow. No tree is inferred.

    If -s is not supplied, all newly added species are treated as one clade.
    Genuine genes left outside an MCL cluster are emitted as singleton OGs.
    """
    assert options.qFastAdd and options.qStopAfterGroups

    print("[OG-ONLY] -og recognised: no new MSAs, trees, or orthologues will be inferred")

    if q_hogs:
        ogs = acc.read_hogs(continuationDir, "N0")
    else:
        ogs = acc.get_original_orthogroups()

    # Add profile-assigned genes to the existing/core orthogroups.
    ogs_new_species, _ = acc.assign_genes(results_files)
    clustersFilename_pairs = acc.write_all_orthogroups(
        ogs, ogs_new_species, []
    )

    # write_all_orthogroups() updates `ogs` in place, so these FASTAs contain
    # genes not assigned to the existing/core OGs.
    n_unassigned = acc.write_unassigned_fasta(ogs, None, speciesInfoObj)

    core_species = set(speciesInfoObj.get_original_species())
    new_species = sorted(
        i_sp for i_sp in speciesInfoObj.speciesToUse
        if i_sp not in core_species
    )

    n_new_unassigned = sum(n_unassigned[i_sp] for i_sp in new_species)
    print("[OG-ONLY] %d unassigned gene(s) remain in the newly added species" % n_new_unassigned)

    if n_new_unassigned == 0:
        return clustersFilename_pairs

    util.PrintUnderline(
        "Inferring remaining orthogroups from unassigned genes (no MSA/tree inference)",
        qHeavy=True,
    )

    # Choose the search/clustering clades.
    if options.speciesTreeFN is not None:
        print("[OG-ONLY] Using user-supplied rooted species tree for clade-specific orthogroups")
        ogSet = orthogroups_set.OrthoGroupsSet(
            options.min_seq,
            files.FileHandler.GetWorkingDirectory1_Read(),
            speciesInfoObj.speciesToUse,
            speciesInfoObj.nSpAll,
            options.qAddSpeciesToIDs,
            options.tree_program,
            idExtractor=util.FirstWordExtractor,
        )
        species_tree_ids_fn = files.FileHandler.GetSpeciesTreeUnrootedFN()
        infer_trees.ConvertUserSpeciesTree(
            options.speciesTreeFN,
            ogSet.SpeciesDict(),
            species_tree_ids_fn,
        )

        # The supplied tree is expected to contain all current species.
        tree_species = set(map(int, tree.Tree(species_tree_ids_fn, format=1).get_leaf_names()))
        missing_species = set(speciesInfoObj.speciesToUse) - tree_species
        if missing_species:
            species_dict = ogSet.SpeciesDict()
            missing_names = [species_dict[str(i)] for i in sorted(missing_species)]
            print("ERROR: The species tree supplied with -s is missing species: %s" % ", ".join(missing_names))
            util.Fail()

        species_clades = acc.get_new_species_clades(
            species_tree_ids_fn, core_species
        )
        species_dict = ogSet.SpeciesDict()
        util.PrintUnderline(
            "Identifying clade-specific orthogroups for the following clades:",
            qHeavy=True,
        )
        for i_clade, clade in enumerate(species_clades):
            print(
                "%d: %s"
                % (i_clade, ", ".join(species_dict[str(i_sp)] for i_sp in clade))
            )
        print("")
    else:
        print("[OG-ONLY] No -s supplied: clustering unassigned genes from all new species together")
        species_clades = [new_species]

    # Build DBs from unassigned FASTAs and perform only the searches required
    # by the selected clades. This is the same clade-aware mechanism used by
    # the normal --assign workflow when a species tree is available.
    run_commands.CreateSearchDatabases(
        speciesInfoObj, options, prog_caller, q_unassigned_genes=True
    )
    run_commands.RunSearch(
        options,
        speciesInfoObj,
        seqsInfo,
        prog_caller,
        n_genes_per_species=n_unassigned,
        species_clades=species_clades,
    )

    options.v2_scores = True
    clusters_unassigned_all = []
    for i_clade, clade in enumerate(species_clades):
        if not any(n_unassigned[i_sp] for i_sp in clade):
            continue

        util.PrintUnderline(
            "OrthoFinder clustering on unassigned-gene clade %d of %d"
            % (i_clade + 1, len(species_clades))
        )
        speciesInfo_clade = copy.deepcopy(speciesInfoObj)
        speciesInfo_clade.speciesToUse = clade
        seqsInfo_clade = util.SeqsInfoRecompute(seqsInfo, clade)
        clusters_unassigned = gathering.DoOrthogroups(
            options,
            speciesInfo_clade,
            seqsInfo_clade,
            speciesNamesDict,
            speciesXML=None,
            i_unassigned=i_clade,
        )
        clusters_unassigned_all.append(clusters_unassigned)

    ogs_unassigned_lists = [
        mcl.GetPredictedOGs(filename)
        for filename in clusters_unassigned_all
    ]

    # MCL intentionally omits singleton rows in an unassigned-gene run.
    # Recover genuine singleton OGs from the authoritative unassigned FASTAs.
    singleton_ogs = _get_singleton_orthogroups(
        ogs, ogs_unassigned_lists, new_species
    )
    print("[OG-ONLY] Added %d singleton orthogroup(s)" % len(singleton_ogs))

    final_og_lists = list(ogs_unassigned_lists)
    if singleton_ogs:
        final_og_lists.append(singleton_ogs)

    clustersFilename_pairs = acc.write_all_orthogroups(
        ogs, {}, final_og_lists
    )
    return clustersFilename_pairs


def BetweenCoreOrthogroupsWorkflow(
        continuationDir,
        speciesInfoObj,
        seqsInfo,
        options,
        prog_caller,
        speciesNamesDict,
        results_files,
        q_hogs,
    ):
    """Infer clade-specific/new orthogroups for species added with --assign."""
    if q_hogs:
        ogs = acc.read_hogs(continuationDir, "N0")
    else:
        ogs = acc.get_original_orthogroups()

    i_og_restart = 0
    ogs_new_species, _ = acc.assign_genes(results_files)
    clustersFilename_pairs = acc.write_all_orthogroups(
        ogs, ogs_new_species, []
    )

    ogSet = orthogroups_set.OrthoGroupsSet(
        options.min_seq,
        files.FileHandler.GetWorkingDirectory1_Read(),
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        options.tree_program,
        idExtractor=util.FirstWordExtractor,
    )

    n_unassigned = acc.write_unassigned_fasta(ogs, None, speciesInfoObj)

    iSpeciesCore = set(speciesInfoObj.get_original_species())
    new_species = [
        i_sp for i_sp in speciesInfoObj.speciesToUse
        if i_sp not in iSpeciesCore
    ]

    # Custom --assign --og path. Infer remaining OGs directly from unassigned
    # genes and stop before any MSA, gene-tree, species-tree or orthologue work.
    if options.qStopAfterGroups and options.speciesTreeFN is None:
        if not any(n_unassigned[i_sp] for i_sp in new_species):
            print("No unassigned genes remain after assignment")
            return clustersFilename_pairs, i_og_restart

        util.PrintUnderline(
            "Inferring orthogroups for unassigned genes without tree inference",
            qHeavy=True,
        )

        run_commands.CreateSearchDatabases(
            speciesInfoObj, options, prog_caller, q_unassigned_genes=True
        )
        run_commands.RunSearch(
            options,
            speciesInfoObj,
            seqsInfo,
            prog_caller,
            n_genes_per_species=n_unassigned,
            species_clades=[new_species],
        )

        options.v2_scores = True
        speciesInfo_new = copy.deepcopy(speciesInfoObj)
        speciesInfo_new.speciesToUse = new_species
        seqsInfo_new = util.SeqsInfoRecompute(seqsInfo, new_species)
        clustersFilename_pairs_unassigned = gathering.DoOrthogroups(
            options,
            speciesInfo_new,
            seqsInfo_new,
            speciesNamesDict,
            speciesXML=None,
            i_unassigned=0,
        )
        ogs_unassigned = mcl.GetPredictedOGs(clustersFilename_pairs_unassigned)
        singleton_ogs = _get_singleton_orthogroups(
            ogs, [ogs_unassigned], new_species
        )
        print("Added %d singleton orthogroups" % len(singleton_ogs))

        clustersFilename_pairs = acc.write_all_orthogroups(
            ogs, {}, [ogs_unassigned, singleton_ogs]
        )
        return clustersFilename_pairs, i_og_restart

    # Normal --assign workflow below is unchanged, apart from singleton recovery
    # when the user supplies a species tree and also requests -og/--og.
    if options.speciesTreeFN is None:
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
            old_version=options.old_version,
            userSpeciesTree=None,
            qStopAfterSeqs=False,
            qStopAfterAlign=False,
            qMSA=options.qMSATrees,
            qPhyldog=False,
            results_name=options.name,
            root_from_previous=True,
            n_skip=options.n_skip
        )
        astral_fn = files.FileHandler.GetAstralFilename()
        astral.create_input_file(
            files.FileHandler.GetOGsTreeDir(), astral_fn, n_skip=options.n_skip
        )
        species_tree_unrooted_fn = files.FileHandler.GetSpeciesTreeUnrootedFN()
        parallel_task_manager.RunCommand(
            astral.get_astral_command(
                astral_fn, species_tree_unrooted_fn, options.nBlast
            )
        )
        core_rooted_species_tree = tree.Tree(
            files.FileHandler.GetCoreSpeciesTreeIDsRootedFN(), format=1
        )
        species_to_speices_map = lambda x: x
        rooted_species_tree_ids, qHaveSupport = trees2ologs_of.CheckAndRootTree(
            species_tree_unrooted_fn, core_rooted_species_tree, species_to_speices_map
        )
        if rooted_species_tree_ids is None:
            print(
                "ERROR: Species tree inference failed. Please check for errors "
                "and check the species tree files: \n%s \n%s"
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
        i_og_restart = len(ogs)
    else:
        util.PrintUnderline("Using user-supplied species tree")
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        infer_trees.ConvertUserSpeciesTree(
            options.speciesTreeFN, ogSet.SpeciesDict(), spTreeFN_ids
        )
        rooted_species_tree_fn = spTreeFN_ids

    species_clades = acc.get_new_species_clades(rooted_species_tree_fn, iSpeciesCore)
    util.PrintUnderline(
        "Identifying clade-specific orthogroups for the following clades:", qHeavy=True
    )
    species_dict = ogSet.SpeciesDict()
    for i, clade in enumerate(species_clades):
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
    print("")

    run_commands.CreateSearchDatabases(
        speciesInfoObj, options, prog_caller, q_unassigned_genes=True
    )
    run_commands.RunSearch(
        options,
        speciesInfoObj,
        seqsInfo,
        prog_caller,
        n_genes_per_species=n_unassigned,
        species_clades=species_clades,
    )

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

    if options.qStopAfterGroups:
        singleton_ogs = _get_singleton_orthogroups(
            ogs, ogs_clade_specific_list, new_species
        )
        print("Added %d singleton orthogroups" % len(singleton_ogs))
        ogs_clade_specific_list.append(singleton_ogs)

    clustersFilename_pairs = acc.write_all_orthogroups(
        ogs, {}, ogs_clade_specific_list
    )
    return clustersFilename_pairs, i_og_restart


def _check_dependencies(options, user_specified_M, prog_caller):
    """Run upstream dependency checks, with --assign --og special handling."""
    args = (
        options,
        user_specified_M,
        prog_caller,
        files.FileHandler.GetWorkingDirectory1_Read()[0],
    )

    if not (options.qFastAdd and options.qStopAfterGroups):
        check_dependencies.CheckDependencies(*args)
        return

    # The custom path performs MCL clustering, so explicitly require MCL.
    if not check_dependencies.CanRunMCL():
        util.Fail()

    # Upstream currently requires ASTRAL for every --assign run. The --og
    # branch returns before ASTRAL/tree inference, so suppress only that check.
    original_can_run_astral = check_dependencies.CanRunASTRAL
    check_dependencies.CanRunASTRAL = lambda: True
    try:
        check_dependencies.CheckDependencies(*args)
    finally:
        check_dependencies.CanRunASTRAL = original_can_run_astral


def main(args=None):
    files.FileHandler.reset()
    start = time.perf_counter()
    try:
        if args is None:
            args = sys.argv[1:]

        # Current OrthoFinder accepts -og and --only-groups. Add --og as a
        # convenient alias without requiring process_args.py to be replaced.
        args = ["-og" if arg == "--og" else arg for arg in args]
        input_args = args.copy()

        ptm = parallel_task_manager.ParallelTaskManager_singleton()
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
            f"[bold dark_goldenrod]OrthoFinder[/bold dark_goldenrod] version "
            f"[deep_sky_blue2]{__version__}[/deep_sky_blue2]",
            end="",
        )
        printer.print(
            " Copyright (C) 2014 [bold dark_goldenrod]David Emms[/bold dark_goldenrod]\n"
        )

        files.InitialiseFileHandler(
            options,
            fastaDir,
            continuationDir,
            resultsDir_nonDefault,
            pickleDir_nonDefault,
        )
        printer.print("Results directory:")
        printer.print(f"    [dark_cyan]{files.FileHandler.GetResultsDirectory1()}")

        _check_dependencies(options, user_specified_M, prog_caller)

        if options.qStartFromBlast and options.qStartFromFasta:
            speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast
            )
            printer.print(f"\nAdding new species in [dark_cyan]{fastaDir}")
            printer.print(f"to existing analysis in [dark_cyan]{continuationDir}")
            speciesInfoObj = fasta_processor.ProcessesNewFasta(
                fastaDir, options.dna, speciesInfoObj, speciesToUse_names
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            speciesXML = None
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            clustersFilename_pairs = gathering.DoOrthogroups(
                options,
                speciesInfoObj,
                seqsInfo,
                speciesNamesDict,
                speciesXML,
            )
            _write_single_copy_orthogroups(
                clustersFilename_pairs, speciesInfoObj
            )
            if options.fix_files or not options.qStopAfterGroups:
                GetOrthologues(
                    seqsInfo, speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromFasta:
            speciesInfoObj = fasta_processor.ProcessesNewFasta(fastaDir, options.dna)
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            speciesXML = None
            util.PrintUnderline("Dividing up work for BLAST for parallel processing")
            run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
            run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            clustersFilename_pairs = gathering.DoOrthogroups(
                options,
                speciesInfoObj,
                seqsInfo,
                speciesNamesDict,
                speciesXML,
            )
            _write_single_copy_orthogroups(
                clustersFilename_pairs, speciesInfoObj
            )
            if options.fix_files or not options.qStopAfterGroups:
                GetOrthologues(
                    seqsInfo,
                    speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromBlast:
            commands_fn = os.path.join(
                files.FileHandler.GetWorkingDirectory1_Read()[0], "blast_commands.txt"
            )
            if os.path.exists(commands_fn):
                commands = []
                with open(commands_fn) as reader:
                    for line in reader:
                        commands.append(line.strip())
                print("Using %d thread(s)" % options.nBlast)
                util.PrintTime("This may take some time...")
                program_caller.RunParallelCommands(
                    options.nBlast,
                    commands,
                    method_threads=options.method_threads,
                    method_threads_large=options.method_threads_large,
                    method_threads_small=options.method_threads_small,
                    threshold=options.threshold,
                    cmd_order=options.cmd_order,
                    tasksize=None,
                    qListOfList=False,
                    q_print_on_error=True,
                    q_always_print_stderr=False,
                    old_version=options.old_version,
                    dynamic_threads=options.dynamic_threads,
                )
            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast
            )
            files.FileHandler.LogSpecies()
            print(
                "Using previously calculated BLAST results in %s"
                % (files.FileHandler.GetWorkingDirectory1_Read()[0])
            )
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            speciesXML = None
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            clustersFilename_pairs = gathering.DoOrthogroups(
                options,
                speciesInfoObj,
                seqsInfo,
                speciesNamesDict,
                speciesXML,
            )
            _write_single_copy_orthogroups(
                clustersFilename_pairs, speciesInfoObj
            )
            if options.fix_files or not options.qStopAfterGroups:
                GetOrthologues(
                    seqsInfo, speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    speciesXML=speciesXML,
                )

        elif options.qStartFromGroups:
            check_blast = not options.qMSATrees
            speciesInfoObj, _ = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(),
                options.qDoubleBlast,
                check_blast=check_blast,
            )
            files.FileHandler.LogSpecies()
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            speciesXML = None
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            clustersFilename_pairs = gathering.DoOrthogroups(
                options,
                speciesInfoObj,
                seqsInfo,
                speciesNamesDict,
                speciesXML,
            )
            _write_single_copy_orthogroups(
                clustersFilename_pairs, speciesInfoObj
            )
            GetOrthologues(
                seqsInfo, speciesNamesDict,
                speciesInfoObj,
                options,
                prog_caller,
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
                options.old_version,
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
                options.fewer_open_files,
                old_version=options.old_version,
                q_split_para_clades=options.qSplitParaClades,
                i_og_restart=0,
                speciesXML=None,
            )

        elif options.qFastAdd:
            speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(
                files.FileHandler.GetWorkingDirectory1_Read(),
                options.qDoubleBlast,
                check_blast=False,
            )
            if not acc.check_for_orthoxcelerate(continuationDir, speciesInfoObj):
                util.Fail()
            util.PrintUnderline("Creating orthogroup profiles")
            wd_list = files.FileHandler.GetWorkingDirectory1_Read()
            fn_diamond_db, q_hogs = acc.prepare_accelerate_database(
                options.min_seq,
                continuationDir,
                wd_list,
                speciesInfoObj.nSpAll,
                tree_program=options.tree_program,
            )
            printer.print(f"\nAdding new species in [dark_cyan]{fastaDir}")
            printer.print(f"to existing analysis in [dark_cyan]{continuationDir}")
            speciesInfoObj = fasta_processor.ProcessesNewFasta(
                fastaDir, options.dna, speciesInfoObj, speciesToUse_names
            )
            options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
            seqsInfo = util.GetSeqsInfo(
                files.FileHandler.GetWorkingDirectory1_Read(),
                speciesInfoObj.speciesToUse,
                speciesInfoObj.nSpAll,
            )
            results_files = run_commands.RunSearch_accelerate(
                options, speciesInfoObj, fn_diamond_db, prog_caller
            )
            speciesNamesDict = species_info.SpeciesNameDict(
                files.FileHandler.GetSpeciesIDsFN()
            )
            if options.qStopAfterGroups:
                # HARD OG-ONLY PATH. Do not enter BetweenCoreOrthogroupsWorkflow,
                # because the normal helper contains MSA/gene-tree/species-tree
                # inference needed by full --assign runs.
                clustersFilename_pairs = AssignOnlyOrthogroupsWorkflow(
                    continuationDir,
                    speciesInfoObj,
                    seqsInfo,
                    options,
                    prog_caller,
                    speciesNamesDict,
                    results_files,
                    q_hogs,
                )
                # Internal unassigned-clade MCL runs create numbered temporary
                # cluster files. The final user-facing Orthogroups directory should
                # nevertheless match a fresh -f run and use the canonical filenames.
                files.FileHandler.iResultsVersion = 0
                gathering.post_clustering_orthogroups(
                    clustersFilename_pairs,
                    speciesInfoObj,
                    seqsInfo,
                    speciesNamesDict,
                    options,
                    speciesXML=None,
                )
                single_copy_fn = _write_single_copy_orthogroups(
                    clustersFilename_pairs, speciesInfoObj
                )
                print("[OG-ONLY] Wrote standard Orthogroups outputs including %s" % os.path.basename(single_copy_fn))
                print("[OG-ONLY] Orthogroups written. Stopping before MSA/tree/orthologue inference.")
                return

            if orphan_genes_version == 2:
                clustersFilename_pairs, i_og_restart = BetweenCoreOrthogroupsWorkflow(
                    continuationDir,
                    speciesInfoObj,
                    seqsInfo,
                    options,
                    prog_caller,
                    speciesNamesDict,
                    results_files,
                    q_hogs,
                )
                gathering.post_clustering_orthogroups(
                    clustersFilename_pairs,
                    speciesInfoObj,
                    seqsInfo,
                    speciesNamesDict,
                    options,
                    speciesXML=None,
                )
                if options.speciesTreeFN is None:
                    options.speciesTreeFN = files.FileHandler.GetSpeciesTreeResultsFN(
                        None, True
                    )
            if options.fix_files or not options.qStopAfterGroups:
                GetOrthologues(
                    seqsInfo, speciesNamesDict,
                    speciesInfoObj,
                    options,
                    prog_caller,
                    i_og_restart,
                    speciesXML=None,
                )

        else:
            raise NotImplementedError

        if not options.save_space and not options.qFastAdd:
            util.split_ortholog_files(
                files.FileHandler.GetOrthologuesDirectory()
            )

        # --assign --og deliberately creates no gene trees or reconciled trees.
        # Skip the normal tree collation/cleanup in this mode.
        if not (options.qFastAdd and options.qStopAfterGroups):
            gene_tree_dir = files.FileHandler.GetOGsTreeDir(qResults=True)
            resolved_gene_tree_dir = files.FileHandler.GetOGsReconTreeDir(qResults=True)
            usr_resolved_gene_tree_fn = files.FileHandler.GetUserResolvedTreeFN()
            util.compress_files(resolved_gene_tree_dir, usr_resolved_gene_tree_fn)

            if options.rm_gene_trees:
                util.cleanup_path(gene_tree_dir)
            if options.rm_resolved_gene_trees:
                util.cleanup_path(resolved_gene_tree_dir)

        d_results = (
            os.path.normpath(files.FileHandler.GetResultsDirectory1()) + os.path.sep
        )
        printer.print("\nResults directory:")
        printer.print(f"    [dark_cyan]{d_results}")
        util.PrintCitation(d_results)
        files.FileHandler.WriteToLog("OrthoFinder run completed\n", True)

    except Exception as e:
        print(str(e))
        util.print_traceback(e)
        ptm.Stop()
        sys.exit(1)

    except KeyboardInterrupt:
        printer.print("\nProgram terminated by user.", style="error")
        sys.exit(1)

    finally:
        ptm.Stop()
        end = time.perf_counter()
        time_elapsed = end - start
        print()
        if len(input_args) == 0 or input_args[0] in [
            "--help", "-h", "-v", "--version", "-sm", "--scoring-matrix"
        ]:
            sys.exit()

        printer.print(
            f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] finished in ", end=""
        )
        printer.print(f"[green]{time_elapsed:5f}[/green]s", end="\n" * 2)
        files.FileHandler.reset()
        sys.exit()


if __name__ == "__main__":
    main()
