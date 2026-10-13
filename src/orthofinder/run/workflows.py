"""The run modes of OrthoFinder: each a function(ctx) -> (species_info_obj, species_names_dict),
for run.main to run and then finish (run.finish). --assign is in run.assign_run."""

import glob
import os  # Y
import os.path  # Y

from rich.markup import escape

from ..mcl_orthogroups import infer_orthogroups
from ..orthology import workflow
from ..search import runner
from ..utils import fasta_processor, file_handler, printing, species_info
from ..utils.printing import printer
from . import previous_run, process_args


def get_orthologs(
    seqs_info,
    species_names_dict,
    species_info_obj,
    options,
    prog_caller,
    i_og_restart=0,
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
    )

def add_species(ctx):
    """-f with -b: add the species in fasta_dir to the analysis being continued"""
    options = ctx.options
    species_info_obj, species_to_use_names = previous_run.process_previous_files(
        file_handler.handler.get_working_directory1_read(), options.q_double_blast
    )
    printer.print(f"\nAdding new species in [dark_cyan]{ctx.fasta_dir}")
    printer.print(f"to existing analysis in [dark_cyan]{ctx.continuation_dir}")
    species_info_obj = fasta_processor.processes_new_fasta(
        ctx.fasta_dir, options.dna, species_info_obj, species_to_use_names
    )
    return _from_new_fasta(ctx, species_info_obj)


def from_fasta(ctx):
    """-f: a new analysis of the species in fasta_dir"""
    species_info_obj = fasta_processor.processes_new_fasta(ctx.fasta_dir, ctx.options.dna)
    return _from_new_fasta(ctx, species_info_obj)


def _from_new_fasta(ctx, species_info_obj):
    """-f, -f -b: from the FASTA files of the species, now in the working directory"""
    file_handler.handler.log_species()
    options, seqs_info = _prepare(ctx.options, species_info_obj)
    _search(ctx, options, species_info_obj, seqs_info)
    species_names_dict = _orthogroups_and_orthologs(ctx, options, species_info_obj, seqs_info)
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict


def _prepare(options, species_info_obj):
    """The options for these species and their sequences."""
    options = process_args.check_options(options, species_info_obj.species_to_use)
    seqs_info = species_info.get_seqs_info(
        file_handler.handler.get_working_directory1_read(),
        species_info_obj.species_to_use,
        species_info_obj.n_sp_all,
    )
    return options, seqs_info


def _search(ctx, options, species_info_obj, seqs_info):
    """The sequence search databases and the searches."""
    printing.print_underline("Dividing up work for BLAST for parallel processing")
    runner.create_search_databases(species_info_obj, options, ctx.prog_caller)
    with ctx.step("Sequence search"):
        runner.run_search(options, species_info_obj, seqs_info, ctx.prog_caller)


def _orthogroups_and_orthologs(ctx, options, species_info_obj, seqs_info):
    """The orthogroups from the search results, then (unless -og, --no-fix-files) the orthologs."""
    species_names_dict = previous_run.species_name_dict(
        file_handler.handler.get_species_ids_fn()
    )
    with ctx.step("Infer orthogroups"):
        infer_orthogroups.do_orthogroups(
            options, species_info_obj, seqs_info, species_names_dict
        )
    if options.fix_files and not options.q_stop_after_mcl_groups:
        get_orthologs(
            seqs_info,
            species_names_dict,
            species_info_obj,
            options,
            ctx.prog_caller,
        )
    return species_names_dict


def from_blast(ctx):
    """-b: continue from the sequence search results of a previous run"""
    options = ctx.options
    prog_caller = ctx.prog_caller
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
        with ctx.step("Sequence search"):
            runner.recreate_search_results(to_create, working_dirs[0], options, prog_caller)
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
    options, seqs_info = _prepare(options, species_info_obj)
    species_names_dict = _orthogroups_and_orthologs(ctx, options, species_info_obj, seqs_info)
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict


def from_groups(ctx):
    """-fg: continue from the orthogroups of a previous run"""
    options = ctx.options
    prog_caller = ctx.prog_caller
    i_og_restart_state = ctx.i_og_restart_state
    check_blast = not options.q_msa_trees
    species_info_obj, _ = previous_run.process_previous_files(
        file_handler.handler.get_working_directory1_read(),
        options.q_double_blast,
        check_blast=check_blast,
    )
    file_handler.handler.log_species()
    options = process_args.check_options(options, species_info_obj.species_to_use)

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
    )
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict


def from_gene_trees(ctx):
    """-fgt: continue from the gene trees of a previous run"""
    options = ctx.options
    i_og_restart_state = ctx.i_og_restart_state
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
    )
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict


def from_species_trees(ctx):
    """-fst: continue from the species tree of a previous run"""
    options = ctx.options
    i_og_restart_state = ctx.i_og_restart_state
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
    )
    ctx.options = options  # check_options may make new options
    return species_info_obj, species_names_dict
