"""The start of a run: its options, results directory, log and restart state."""

from rich.markup import escape

from orthofinder.utils import logging as run_logging

from .. import __version__
from ..search import diamond_scoring
from ..utils import file_handler, file_io, printing, run_state
from ..utils.printing import printer
from . import check_dependencies, previous_run, process_args


def start_run(ctx, args, input_args):
    """
    Read the arguments; for a restart, restore the state of the run being
    continued; set up the results directory and the log; check the
    programs. Fills ctx.
    """
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
    state_dir = None  # the results directory of state, for its paths
    changed = []
    q_restart = ((options.q_start_from_blast and not options.q_start_from_fasta) or options.q_start_from_groups
                 or options.q_start_from_trees or options.q_start_from_species_trees)
    if q_restart and continuation_dir is not None:
        state_dir, state = run_state.find(continuation_dir)
        if not state:
            # no run_state.json (older versions, or the file removed): what
            # Log.txt records, read before this run rewrites it
            state = previous_run.state_from_log(continuation_dir)
        i_og_restart_state = state.get("assign", {}).get("i_og_restart", 0)
        changed = run_state.restore_settings(state, options, input_args, results_dir=state_dir)
        # the run's search program (restored above, or the same as the default)
        options.search_program_from_run = ("search_program" in state.get("settings", {})
                                           and not options.search_program_given)
    core_scoring = []
    if options.q_fast_add and continuation_dir is not None:
        # --assign: the new species' searches are scored as the core
        # analysis' were, unless --matrix etc. are given
        state_dir, state = run_state.find(continuation_dir)
        core_scoring = run_state.restore_settings(state, options, input_args, names=run_state.SCORING_SETTINGS,
                                                  results_dir=state_dir)
    # --matrix etc. that differ from the run being continued (or the core analysis)
    scoring_changes = run_state.scoring_differences(state, options, state_dir) if options.score_matrix_given else []

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
    if diamond_scoring.is_custom_matrix(options.score_matrix):
        # searched with a copy in the working directory, kept with the results
        options.score_matrix = diamond_scoring.copy_custom_matrix(
            options.score_matrix, file_handler.handler.get_working_directory_write())


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
    ctx.log = log
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
    printer.print(f"    [dark_cyan]{escape(file_handler.handler.get_results_directory1())}")

    # A restart also continues with the species tree the run used (unless
    # -s is given), and its point for restarting gene-tree inference
    # (--assign); its settings were restored above.
    results_dir = file_handler.handler.get_results_directory1()
    if core_scoring:
        printer.print("Scoring the searches as the core analysis did: %s"
                      % escape(", ".join("%s=%s" % c for c in core_scoring)))
        file_handler.handler.write_to_log("Settings of the core analysis: %s\n"
                                          % ", ".join("%s=%s" % c for c in core_scoring))
    if scoring_changes:
        text = ("NOTE: the searches of this run are scored differently from those of the %s: %s"
                % ("core analysis" if options.q_fast_add else "run being continued",
                   ", ".join("%s %s (was %s)" % (name, now, earlier) for name, earlier, now in scoring_changes)))
        printer.print(escape(text))
        file_handler.handler.write_to_log(text + "\n")
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
        recorded_programs=state.get("programs"),
    )
    if continuation_dir is not None:
        printing.print_restart_info(options, continuation_dir)
    ctx.options = options
    ctx.prog_caller = prog_caller
    ctx.fasta_dir = fasta_dir
    ctx.continuation_dir = continuation_dir
    ctx.results_dir = results_dir
    ctx.state = state
    ctx.i_og_restart_state = i_og_restart_state
    ctx.assign_resume = assign_resume
