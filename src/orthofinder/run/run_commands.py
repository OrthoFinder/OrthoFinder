from ..utils import util, program_caller, files, parallel_task_manager, file_io, fasta_processor
from ..utils.util import printer
from ..utils import logging as run_logging
import subprocess
import glob
import shutil
from . import run_info, species_info

# from .. import my_env
import os
import re
import time


def RunBlastDBCommand(command):
    capture = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=parallel_task_manager.my_env,
        shell=True,
    )
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stdout_lines = stdout.count("\n")
    n_stderr_lines = stderr.count("\n")
    nLines_success = 12
    if n_stdout_lines > nLines_success or n_stderr_lines > 0 or capture.returncode != 0:
        print("\nWARNING: Likely problem with input FASTA files")
        if capture.returncode != 0:
            print("makeblastdb returned an error code: %d" % capture.returncode)
        else:
            print("makeblastdb produced unexpected output")
        print("Command: %s" % " ".join(command))
        print("stdout:\n-------")
        print(stdout)
        if len(stderr) > 0:
            print("stderr:\n-------")
            print(stderr)


# 7
def RunSearch(
    options,
    speciessInfoObj,
    seqsInfo,
    prog_caller,
    q_new_species_unassigned_genes=False,
    n_genes_per_species=None,
    species_clades=None,
    only_missing=False,
):
    """
    n_genes_per_species: List[int] - optional, for use with unassigned genes. If a species has zero unassigned genes, don't search with/agaisnt it.
    only_missing: run only the searches whose results are missing or incomplete (continuing a run)
    """
    name_to_print = options.search_program
    if options.qStopAfterPrepare:
        util.PrintUnderline("%s commands that must be run" % name_to_print)
    elif species_clades is not None:
        util.PrintUnderline(
            "Running %s searches for new non-core species clades" % name_to_print
        )
    else:
        util.PrintUnderline("Running %s all-versus-all" % name_to_print, True)
    tasksizes = None
    if species_clades is None:
        # print("running GetOrderedSearchCommands")
        commands, tasksizes = run_info.GetOrderedSearchCommands(
            seqsInfo,
            speciessInfoObj,
            options,
            prog_caller,
            n_genes_per_species,
            q_new_species_unassigned_genes=q_new_species_unassigned_genes,
        )

    else:
        # print("running GetOrderedSearchCommands_clades")
        commands, tasksizes = run_info.GetOrderedSearchCommands_clades(
            seqsInfo,
            speciessInfoObj,
            options,
            prog_caller,
            n_genes_per_species,
            species_clades,
        )

    # if options.qStopAfterPrepare:
    # if options.save_blast_commands:
    # The columns of the results depend on the program and its output options:
    # check now that OrthoFinder can read them. They are read later from the
    # commands, which are kept next to the results in blast_commands.txt.
    fmt = check_search_output_format(commands[0], options.search_program, prog_caller) if commands else None
    commands_fn = files.FileHandler.GetBALSATCommandFN()
    if options.qStopAfterPrepare:
        # -op: the commands are run outside OrthoFinder, which would otherwise
        # fill in the thread count when running them
        saved = [program_caller.FillMethodThreads(c, options.method_threads)[1] for c in commands]
    else:
        saved = commands
    with open(commands_fn, "w") as writer:
        for line in file_io.commands_file_lines(saved, fmt):
            writer.write(line + "\n")
    if options.qStopAfterPrepare:
        # -op: stop here; the searches are run later (e.g. on a cluster), or
        # by restarting with -b, which runs the commands of missing results
        from rich.markup import escape
        printer.print("The %d search commands have been saved to:" % len(commands))
        printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(os.path.abspath(commands_fn)))
        printer.print("Run them yourself (e.g. on a cluster), or let OrthoFinder run them: a restart "
                      "with -b runs the commands of any search results that are missing or incomplete.")
        util.PrintRunEnd(files.FileHandler.GetResultsDirectory1(), "-b", "the sequence search")
        util.Success()
    if only_missing:
        # continuing a run (-b): only the searches whose results are missing
        # or incomplete (checkpoint.txt), e.g. not those run on a cluster
        to_create = dict(species_info.ResultsToCreate([results_written_by(c) for c in commands]))
        for fn, path in to_create.items():
            if path is not None:
                os.replace(path, path + ".incomplete")
                printer.print("Incomplete search results moved aside: %s.incomplete" % path)
        commands = [c for c in commands if results_written_by(c) in to_create]
        printer.print("%d of the searches are to be run (the others' results are complete)" % len(commands))
    # Logged only when the searches run (not with -op): a restart judges from
    # these whether the searches were interrupted (file_io.search_checkpoint).
    # The all-versus-all search is logged as a step by main.
    log_step = species_clades is not None and len(commands) > 0     # nothing to record if none are run
    if log_step:
        run_logging.RunLogger.step("Sequence search", "Started")
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
        tasksize=tasksizes,
        qListOfList=False,
        q_print_on_error=True,
        q_always_print_stderr=False,
        dynamic_threads=options.dynamic_threads,
        on_success=record_finished_search,
    )

    util.PrintTime("Done all-versus-all sequence search")
    if log_step:
        run_logging.RunLogger.step("Sequence search", "Completed")
    remove_search_only_files(
        [(options.search_program, files.FileHandler.GetSpeciesDatabaseN(i, options.search_program))
         for i in range(speciessInfoObj.nSpAll)],
        prog_caller,
    )
    remove_search_tmp_dirs(commands)     # MMseqs2's temporary folders


# 6
def CreateSearchDatabases(
    speciesInfoObj,
    options,
    prog_caller,
    q_unassigned_genes=False,
    core_infile=None,
    core_outfile=None,
    new_species=False,
):
    if new_species:
        iSpeciesToDo = list(
            range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll)
        )
        qForCreation = True
    else:
        iSpeciesToDo = range(max(speciesInfoObj.speciesToUse) + 1)
        qForCreation = False
    if core_infile is not None and core_outfile is not None:
        command = prog_caller.GetSearchMethodCommand_DB(
            options.search_program,
            core_infile,
            core_outfile,
            options.score_matrix,
            options.gapopen,
            options.gapextend,
            options.method_threads,
            nucleotide=options.dna,
        )
        ret_code = parallel_task_manager.RunCommand(
            command, qPrintOnError=True, qPrintStderr=False
        )
        if ret_code != 0:
            files.FileHandler.LogFailAndExit(
                f"ERROR: {options.search_program} makedb failed"
            )
    else:
        progressbar, task = util.get_progressbar(len(iSpeciesToDo))
        progressbar.update(task, description="[yellow]Creating %s databases" % options.search_program)
        progressbar.start()
        update_cycle = 1  # 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000

        for i, iSp in enumerate(iSpeciesToDo):
            fn_fasta = (
                files.FileHandler.GetSpeciesUnassignedFastaFN(iSp)
                if q_unassigned_genes
                else files.FileHandler.GetSpeciesFastaFN(iSp, qForCreation=qForCreation)
            )
            if os.stat(fn_fasta).st_size == 0:
                if (i + 1) % update_cycle == 0:
                    progressbar.update(task, advance=update_cycle)
                continue

            command = prog_caller.GetSearchMethodCommand_DB(
                options.search_program,
                fn_fasta,
                files.FileHandler.GetSpeciesDatabaseN(iSp, options.search_program),
                options.score_matrix,
                options.gapopen,
                options.gapextend,
                options.method_threads,
                nucleotide=options.dna,
            )
            # util.PrintTime("Creating %s database %d of %d" % (options.search_program, iSp + 1, len(iSpeciesToDo)))
            ret_code = parallel_task_manager.RunCommand(
                command, qPrintOnError=True, qPrintStderr=False
            )
            if ret_code != 0:
                files.FileHandler.LogFailAndExit("ERROR: %s database creation failed" % options.search_program)

            if (i + 1) % update_cycle == 0:
                progressbar.update(task, advance=update_cycle)

        progressbar.stop()


def CreateSearchDatabases_accelerate(speciesInfoObj, options, prog_caller):

    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))

    progressbar, task = util.get_progressbar(len(iSpeciesNew))
    progressbar.update(task, description="[yellow]Creating %s databases" % options.search_program)
    progressbar.start()
    update_cycle = (
        1  # 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000
    )

    for i, iSp in enumerate(iSpeciesNew):
        fn_fasta = files.FileHandler.GetSpeciesFastaFN(iSp, qForCreation=True)
        if os.stat(fn_fasta).st_size == 0:
            if (i + 1) % update_cycle == 0:
                progressbar.update(task, advance=update_cycle)
            continue

        command = prog_caller.GetSearchMethodCommand_DB(
            options.search_program,
            fn_fasta,
            files.FileHandler.GetSpeciesDatabaseN(iSp, options.search_program),
            options.score_matrix,
            options.gapopen,
            options.gapextend,
            options.method_threads,
            nucleotide=options.dna,
        )
        # util.PrintTime("Creating %s database %d of %d" % (options.search_program, iSp + 1, len(iSpeciesToDo)))
        ret_code = parallel_task_manager.RunCommand(
            command, qPrintOnError=True, qPrintStderr=False
        )
        if ret_code != 0:
            files.FileHandler.LogFailAndExit("ERROR: %s database creation failed" % options.search_program)

        if (i + 1) % update_cycle == 0:
            progressbar.update(task, advance=update_cycle)

    progressbar.stop()


# A results file written by a search command: ".../Blast<i>_<j>.txt" (not, e.g.,
# MMseqs2's temporary "/tmp/tmpBlast<i>_<j>.txt").
_RESULTS_IN_COMMAND = re.compile(r"""([^\s'"=<>|;]*?)\bBlast(\d+)_(-?\d+)\.txt""")


def _results_in_command(command):
    """The (directory, "Blast<i>_<j>.txt") results files named in a search command."""
    found = set()
    for m in _RESULTS_IN_COMMAND.finditer(command):
        directory = m.group(1)
        if directory == "" or directory.endswith(("/", os.sep)):
            found.add((directory, "Blast%s_%s.txt" % (m.group(2), m.group(3))))
    return found


def results_written_by(command):
    """The results file ".../Blast<i>_<j>.txt" a search command writes (None if not one)."""
    paths = {directory + name for directory, name in _results_in_command(command)}
    return paths.pop() if len(paths) == 1 else None


def record_finished_search(command):
    """Record in checkpoint.txt that a search finished, with its results file's size (see file_io)."""
    fn = results_written_by(command)
    if fn is None:
        return
    path = fn + ".gz" if os.path.exists(fn + ".gz") else fn
    if os.path.exists(path):
        run_logging.RunLogger.message(file_io.search_completed_message(path))


def search_commands_for(commands, needed_fns, working_dir):
    """
    The saved search commands that write the needed results files.

    Returns (commands, fns with no command, old results directories that were
    replaced by working_dir). The commands were saved with
    absolute paths; if the results directory has since been moved or copied,
    its old path in a command is replaced by working_dir, so the results are
    always written to the working directory being used, never to another run's.
    """
    wd = os.path.join(os.path.abspath(working_dir), "")
    writer_of = {}     # "Blast<i>_<j>.txt" -> (command, the directory in the command)
    for command in commands:
        outputs = _results_in_command(command)
        names = {name for _, name in outputs}
        if len(names) == 1:     # a command that writes one results file
            name = names.pop()
            writer_of.setdefault(name, (command, sorted(d for d, _ in outputs)[-1]))
    selected, no_command, moved_from = [], [], set()
    for fn in needed_fns:
        name = os.path.basename(fn)
        if name not in writer_of:
            no_command.append(fn)
            continue
        command, directory = writer_of[name]
        if directory and os.path.join(os.path.abspath(directory), "") != wd:
            command = command.replace(directory, wd)
            moved_from.add(directory)
        selected.append(command)
    return selected, no_command, sorted(moved_from)


# A species database in a search command: ".../<program>DBSpecies<j>", e.g.
# BlastDBSpecies3 (-S blast), diamondDBSpecies3 (.dmnd), mmseqsDBSpecies3 (.fa).
_DATABASE_IN_COMMAND = re.compile(r"""([^\s'"=<>|;]*?)\b([A-Za-z0-9_]+?)DBSpecies(\d+)\b""")


def recreate_missing_databases(commands, working_dir, prog_caller, options):
    """
    Create the species databases that search commands need and that no longer
    exist: BLAST+ databases are deleted once the searches are done
    (remove_search_only_files). A database is named <program>DBSpecies<j>
    and is made from Species<j>.fa, with makeblastdb for "Blast" (older
    versions' -S blast) or the program's db_cmd from the config. Returns
    the (program, database) pairs the commands use, to be cleaned up after
    the searches as in a normal run.
    """
    needed = {}
    used = set()
    for command in commands:
        for m in _DATABASE_IN_COMMAND.finditer(command):
            directory, program, iSp = m.group(1), m.group(2), int(m.group(3))
            if directory and not directory.endswith(("/", os.sep)):
                continue
            db = directory + "%sDBSpecies%d" % (program, iSp)
            used.add((program, db))
            # the database itself or its files (db.dmnd, db.pin, ...), not
            # those of another species (DBSpecies1 is not DBSpecies10)
            if not (os.path.exists(db) or glob.glob(glob.escape(db) + ".*")):
                needed[db] = (program, iSp)
    for db, (program, iSp) in sorted(needed.items()):
        # The species' FASTA may be in an earlier working directory (species
        # added to a previous analysis): look in all of them.
        try:
            fn_fasta = files.FileHandler.GetSpeciesFastaFN(iSp)
        except Exception:
            fn_fasta = os.path.join(working_dir, "Species%d.fa" % iSp)
        if not os.path.exists(fn_fasta):
            files.FileHandler.LogFailAndExit(
                "ERROR: Cannot create the search database %s: %s does not exist" % (db, fn_fasta))
        if program == "Blast":
            command = " ".join(["makeblastdb", "-dbtype", "prot", "-in", fn_fasta, "-out", db])
            printer.print("Creating the BLAST database for species %d (deleted after the searches)" % iSp)
            RunBlastDBCommand(command)
        else:
            # A restart (-b) is not given -d: the sequences tell the type
            nucleotide = options.dna or fasta_processor.sequence_type(fn_fasta) == "dna"
            command = prog_caller.GetSearchMethodCommand_DB(
                program, fn_fasta, db, options.score_matrix, options.gapopen,
                options.gapextend, options.method_threads, nucleotide=nucleotide)
            printer.print("Creating the %s database for species %d" % (program, iSp))
            if parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=False) != 0:
                files.FileHandler.LogFailAndExit("ERROR: could not create the %s database %s" % (program, db))
    return sorted(used)


# MMseqs2's temporary folder for one search, named in its command (config:
# "/tmp/tmpBASEOUTNAME"), e.g. /tmp/tmpBlast0_3.txt.
_TMP_IN_COMMAND = re.compile(r"""(/\S*/tmp\S*Blast\d+_-?\d+\.txt)(?=[\s;&|]|$)""")


def remove_search_tmp_dirs(commands):
    """Remove the temporary folders the search commands named (MMseqs2), once they have run."""
    for command in commands:
        for d in _TMP_IN_COMMAND.findall(command):
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)


def remove_search_only_files(databases, prog_caller):
    """
    Remove what is only needed while searching, once the searches are done:
      - BLAST+ databases (makeblastdb; "Blast" is older versions' -S blast);
      - MMseqs2 indexes (createindex, .idx files): a k-mer table of fixed
        size, ~0.9 GB per protein database however small. The database is
        kept; MMseqs2 searches it without an index, building one in memory.
    databases: (program, database path) pairs. Any search run later (a
    restart, or species added) creates what it needs again.
    """
    for program, db in databases:
        if program == "Blast" or prog_caller.UsesMakeBlastDB(program):
            remove = glob.glob(glob.escape(db) + ".*")      # not DBSpecies10.* for DBSpecies1
        elif prog_caller.UsesMMseqs(program):
            remove = glob.glob(glob.escape(db) + ".*idx*")
        else:
            continue
        for f in remove:
            if os.path.isfile(f):
                os.remove(f)


def recreate_search_results(to_create, working_dir, options, prog_caller):
    """
    Create the missing or incomplete search results of a restarted run (-b)
    by running their commands from blast_commands.txt in the working
    directory, writing to that directory. to_create: (fn, incomplete file or
    None) from species_info.GetBlastResultsToCreate. Exits with an
    explanation if a result cannot be created.
    """
    from rich.markup import escape
    n_incomplete = sum(1 for _, path in to_create if path is not None)
    util.PrintUnderline("Creating missing search results")
    printer.print("%d required search result(s) missing and %d incomplete (interrupted search) in:"
                  % (len(to_create) - n_incomplete, n_incomplete))
    printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(working_dir))
    commands_fn = os.path.join(working_dir, "blast_commands.txt")
    if not os.path.exists(commands_fn):
        files.FileHandler.LogFailAndExit(
            "ERROR: The search results below are missing or incomplete, and there is no %s "
            "with the commands to create them:\n%s"
            % (commands_fn, "\n".join(path or fn for fn, path in to_create)))
    with open(commands_fn) as reader:
        commands = [line.strip() for line in reader if line.strip() and not line.startswith("#")]
    selected, no_command, moved_from = search_commands_for(
        commands, [fn for fn, _ in to_create], working_dir)
    # The saved thread counts (-op, or edited for a cluster) are not used here:
    # OrthoFinder runs many searches at once, each with its own count (1 by default).
    selected = [program_caller.ResetSearchThreads(c) for c in selected]
    if no_command:
        files.FileHandler.LogFailAndExit(
            "ERROR: %s has no command to create these search results:\n%s"
            % (commands_fn, "\n".join(no_command)))
    if moved_from:
        printer.print("The results were created in [dark_cyan]%s[/dark_cyan], not in this directory: "
                      "the commands are run with this directory instead, so the results are written here."
                      % escape(", ".join(moved_from)), style="warning")
    databases = recreate_missing_databases(selected, working_dir, prog_caller, options)
    for fn, path in to_create:
        if path is not None:
            # Judged cut short (truncated, or no final line break): kept aside
            # rather than deleted, in case it was complete after all, and
            # replaced by the new search.
            os.replace(path, path + ".incomplete")
            printer.print("Incomplete search results moved aside: [dark_cyan]%s.incomplete[/dark_cyan]"
                          % escape(path))
    printer.print("Running [bold]%d[/bold] of the %d commands in [dark_cyan]%s[/dark_cyan]"
                  % (len(selected), len(commands), escape(commands_fn)))
    printer.print("\nUsing %d thread(s)" % options.nBlast)
    util.PrintTime("This may take some time...")
    program_caller.RunParallelCommands(
        options.nBlast,
        selected,
        method_threads=options.method_threads,
        method_threads_large=options.method_threads_large,
        method_threads_small=options.method_threads_small,
        threshold=options.threshold,
        cmd_order=options.cmd_order,
        tasksize=None,
        qListOfList=False,
        q_print_on_error=True,
        q_always_print_stderr=False,
        dynamic_threads=options.dynamic_threads,
        on_success=record_finished_search,
    )
    remove_search_only_files(databases, prog_caller)
    remove_search_tmp_dirs(selected)


def check_search_output_format(command, search_program, prog_caller,
                               required=file_io.REQUIRED_HIT_FIELDS):
    """The columns a search command writes (exits with an explanation if OrthoFinder cannot use them)."""
    try:
        fmt = file_io.search_hit_format(
            command, prog_caller.GetSearchOutputFields(search_program), required=required)
    except ValueError as e:
        print("\nERROR: OrthoFinder cannot read the output of the '%s' search: %s" % (search_program, e))
        print("Use a tab-separated output format with these columns, or declare the columns with "
              "\"output_fields\" in the config file entry for '%s'." % search_program)
        util.Fail()
    if fmt.assumed:
        print("Search program '%s': its output columns are assumed to be the standard BLAST "
              "tabular ones (%s)." % (search_program, " ".join(fmt.fields)))
    return fmt


def RunSearch_accelerate(
    options, speciessInfoObj, fn_diamond_db, prog_caller, q_one_query=False, reuse_complete=False
):
    """reuse_complete: the profile search of the run being continued finished (checkpoint.txt): its results are used if all present."""
    name_to_print = options.search_program
    util.PrintUnderline("Running %s profiles search" % name_to_print)
    commands, tasksizes, results_files = run_info.GetOrderedSearchCommands_accelerate(
        speciessInfoObj,
        fn_diamond_db,
        options,
        prog_caller,
        q_one_query=q_one_query,
        threads=options.nBlast,
    )
    if not commands:
        return [], file_io.STANDARD_HIT_FORMAT       # no new species to search
    # Columns of the results, passed on to their reader (they are read in this run).
    layout = check_search_output_format(
        commands[0], "diamond" if q_one_query else options.search_program, prog_caller,
        required=file_io.REQUIRED_PROFILE_HIT_FIELDS,
    )
    if reuse_complete and all(os.path.exists(fn) for fn in results_files):
        # continuing an --assign run (-b) whose profile search finished
        util.PrintTime("Using the profile search results of the run being continued")
        return results_files, layout
    if q_one_query:
        # Raises with the command's output (recorded in checkpoint.txt) on failure.
        parallel_task_manager.RunCommand(commands[0], qPrintOnError=True, raise_on_error=True)
        util.PrintTime("Done profiles search\n")
        return results_files, layout
    program_caller.RunParallelCommands(
        options.nBlast,
        commands,
        method_threads=options.method_threads,
        method_threads_large=options.method_threads_large,
        method_threads_small=options.method_threads_small,
        threshold=options.threshold,
        cmd_order=options.cmd_order,
        tasksize=tasksizes,
        qListOfList=False,
        q_print_on_error=True,
        q_always_print_stderr=False,
        dynamic_threads=options.dynamic_threads,
    )

    util.PrintTime("Done profiles search")
    # the profiles database is only searched here (BLAST+ files, MMseqs2 index)
    remove_search_only_files([(options.search_program, fn_diamond_db)], prog_caller)
    remove_search_tmp_dirs(commands)
    return results_files, layout
