"""--assign: adding species to the orthogroups of a core analysis."""

import copy  # Y

from rich.markup import escape

from orthofinder.utils import logging as run_logging

from .. import orphan_genes_version
from ..gene_tree import infer_trees
from ..mcl_orthogroups import assign as assign
from ..mcl_orthogroups import infer_orthogroups, mcl, orthogroups_set
from ..orthology import rooting
from ..phylo import taxa, tree
from ..search import runner
from ..species_tree import from_gene_trees
from ..utils import (
    fasta_processor,
    file_handler,
    file_io,
    gene_ids,
    printing,
    run_state,
    species_info,
    util,
)
from ..utils.printing import printer
from . import previous_run, process_args
from .workflows import get_orthologs

try:
    from rich import print
except ImportError:
    ...


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
    # The core orthogroups with the new species' genes, in memory (ogs is updated): the
    # clusters file is written once, with the clade-specific orthogroups (below)
    assign.write_all_orthogroups(ogs, ogs_new_species, [], write=False)

    og_set = orthogroups_set.OrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(),
        species_info_obj.species_to_use,
        species_info_obj.n_sp_all,
        options.q_add_species_to_ids,
        options.tree_program,
        id_extractor=gene_ids.FirstWordExtractor,
        mcl_inflation=options.mcl_inflation,
        ogs=list(ogs),   # (a copy of the list: the clade-specific orthogroups are added to ogs below)
    )

    n_unassigned = assign.write_unassigned_fasta(ogs, None, species_info_obj)

    # Get/Infer species tree
    if options.species_tree_fn is None:
        # Infer gene trees
        # We write orthogroup & stats results files in the following code, which we should avoid & only do once all OGs are done.
        infer_orthogroups.post_clustering_orthogroups(
            None,
            species_info_obj,
            seqs_info,
            species_names_dict,
            options,
            q_incremental=True,
            ogs=ogs,
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
        species_tree_unrooted_fn = file_handler.handler.get_species_tree_unrooted_fn()
        from_gene_trees.infer_species_tree(
            prog_caller, process_args.assign_species_tree_program(options), file_handler.handler.get_ogs_tree_dir(),
            file_handler.handler.get_astral_filename(), species_tree_unrooted_fn, options.n_blast,
            n_skip=options.n_skip)

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

        # The species tree the rest of the run uses (options.species_tree_fn), with the species
        # names, in the working directory: the results species tree files are written once,
        # with the orthologs (workflow.infer_orthologs)
        taxa.rename_tree_taxa(
            rooted_species_tree_ids,
            run_state.species_tree_copy_fn(file_handler.handler.rd1),
            og_set.species_dict(),
            q_support=q_have_support,
            q_fix_negatives=True,
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

def run_assign(ctx):
    """--assign: add the species in fasta_dir to the orthogroups of a core analysis
    (or continue such a run: -b on its results)"""
    options = ctx.options
    prog_caller = ctx.prog_caller
    fasta_dir = ctx.fasta_dir
    continuation_dir = ctx.continuation_dir
    results_dir = ctx.results_dir
    state = ctx.state
    assign_resume = ctx.assign_resume
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
    with ctx.step("Create orthogroup profiles"):
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
        printer.print(f"\nAdding new species in [dark_cyan]{escape(fasta_dir)}")
        printer.print(f"to existing analysis in [dark_cyan]{escape(continuation_dir)}")

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
    # Add genes to orthogroups
    with ctx.step("Search orthogroup profiles"):
        results_files, results_layout = runner.run_search_assign(
            options, species_info_obj, fn_diamond_db, prog_caller,
            reuse_complete=profile_search_done,
        )
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
        with ctx.step("Infer clade-specific orthogroups"):
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

        # Infer clade-specific orthogroup gene trees
        infer_orthogroups.post_clustering_orthogroups(
            clusters_filename_pairs,
            species_info_obj,
            seqs_info,
            species_names_dict,
            options,
        )
        q_inferred_species_tree = options.species_tree_fn is None
        if q_inferred_species_tree:
            # No user species tree, use the one we've just inferred
            options.species_tree_fn = run_state.species_tree_copy_fn(results_dir)
        # kept for a restart (-fg, -fgt, -fst) after a stop (-og, -ogt, -ost)
        file_handler.handler.log_species_tree_used(run_state.keep_species_tree(results_dir, options.species_tree_fn))
        run_state.update(results_dir, "assign", i_og_restart=i_og_restart)
        file_handler.handler.log_assign_first_orthogroup_to_infer(i_og_restart)
    if not options.q_stop_after_mcl_groups:
        get_orthologs(
            seqs_info,
            species_names_dict,
            species_info_obj,
            options,
            prog_caller,
            i_og_restart,
        )
    elif orphan_genes_version == 2 and q_inferred_species_tree:
        # Stopped after the orthogroups (-og): the species tree inferred to find the clades of
        # the clade-specific orthogroups, in the results, as the orthologs step writes it
        # (workflow.infer_orthologs)
        taxa.rename_tree_taxa(
            file_handler.handler.get_species_tree_ids_rooted_fn(),
            file_handler.handler.get_species_tree_results_fn(0, True),
            {str(i): name for i, name in species_names_dict.items()},
            q_support=False,
            q_fix_negatives=True,
        )
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict
