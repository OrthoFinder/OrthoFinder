from ..utils import util, program_caller, file_handler, parallel_task_manager, run_state, programs
from ..gene_tree import msa_trees
from . import process_args
import contextlib
import io
import os
import re
import shutil
from .. import __location__
from ..utils.printing import printer
from ..utils import printing

_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_ERROR_WORDS = re.compile(r"error|not found|cannot|denied|illegal|segmentation|exec format|no such file", re.I)

# The programs whose chosen copy passed a test: not changed by a later test
_confirmed = set()


def _first_error(text):
    """The first line of a failed test's output that looks like the reason."""
    for line in _ANSI.sub("", text).splitlines():
        line = line.strip()
        if line and not line.startswith("Test can run") and _ERROR_WORDS.search(line):
            return line if len(line) <= 200 else line[:200] + " ..."
    return "its test failed"


def _same(a, b):
    return os.path.realpath(a) == os.path.realpath(b)


def _other_copies(names, chosen, path):
    """
    The choices to try after a test with these copies failed: all its
    programs from one directory (e.g. a whole installation that does not
    run on this computer), in the order of path, then each program's other
    copies on their own.
    """
    tries = []
    for d in path.split(os.pathsep):
        if not d:
            continue
        combo = {}
        for name in names:
            location = shutil.which(name, path=d)
            if location is None:
                break
            combo[name] = os.path.abspath(location)
        else:
            tries.append(combo)
    for name in names:
        for location in programs.candidates(name, path, chosen[name])[1:]:
            tries.append(dict(chosen, **{name: location}))
    unique = []
    seen = {tuple(os.path.realpath(chosen[n]) for n in names)}
    for combo in tries:
        key = tuple(os.path.realpath(combo[n]) for n in names)
        if key not in seen:
            seen.add(key)
            unique.append(combo)
    return unique


def run_test(test):
    """
    Run a dependency test, test() -> (success, cmd, result), cmd None if it
    ran nothing (a test turned off in the config), and choose the copy of each
    program cmd runs (programs.py): the one chosen already, else the first on
    PATH. If the test fails, try the other copies of the programs not chosen
    by an earlier test, quietly, and keep the first choice that passes.
    Returns test()'s (success, cmd, result) for the copies chosen.
    """
    success, cmd, result = test()
    if cmd is None:
        return success, cmd, result
    path = parallel_task_manager.my_env["PATH"]
    names = [p for p in programs.command_programs(cmd) if programs.can_choose(p)]
    for name in names:
        if name not in programs.chosen():
            found = programs.candidates(name, path)
            if found:
                programs.choose(name, found[0])
    names = [n for n in names if n in programs.chosen()]
    if success:
        _confirmed.update(names)
        return success, cmd, result
    free = [n for n in names if n not in _confirmed]
    chosen = {n: programs.chosen()[n] for n in free}
    tries = _other_copies(free, chosen, path) if free else []
    if not tries:
        return success, cmd, result
    printer.print("Trying the other copies of %s" % ", ".join(free))
    failed = []
    for combo in tries:
        for name, location in combo.items():
            programs.choose(name, location)
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            attempt = test()
        if attempt[0]:
            for name in free:
                if not _same(combo[name], chosen[name]):
                    printer.print("%s: the copy at %s does not work (see above), so OrthoFinder uses %s"
                                  % (name, chosen[name], combo[name]), style="warning", soft_wrap=True)
            _confirmed.update(names)
            printer.print(f'Test can run "[orange3]{" ".join(free)}[/orange3]" - [bold green]ok')
            return attempt
        failed.append((combo, _first_error(output.getvalue())))
    for name, location in chosen.items():
        programs.choose(name, location)
    printer.print("None of the other copies works:", style="error")
    for combo, reason in failed:
        copies = ", ".join(combo[n] for n in free if not _same(combo[n], chosen[n]))
        printer.print("    %s: %s" % (copies, reason), soft_wrap=True)
    return success, cmd, result


def choose_recorded(recorded):
    """
    Prefer the copies of the programs the run being continued used (if they
    are still there), so that it continues with the same versions.
    """
    _confirmed.clear()
    for name, location in (recorded or {}).items():
        if (programs.can_choose(name) and isinstance(location, str)
                and os.path.isfile(location) and os.access(location, os.X_OK)):
            programs.choose(name, os.path.abspath(location))


def warn_changed_programs(recorded, used):
    """A warning for each program a restart runs from a different copy than the run it continues."""
    for name, location in sorted(used.items()):
        before = (recorded or {}).get(name)
        if isinstance(before, str) and not _same(before, location):
            printer.print("WARNING: The run being continued used %s, this restart uses %s"
                          % (before, location), style="warning", soft_wrap=True)


def record_programs(commands, results_dir):
    """
    The full path of the copy of each program the commands ran (programs.py:
    the first that works, looked for in parallel_task_manager.program_dirs),
    in run_state.json: "programs": {name: path}.
    """
    chosen = programs.chosen()
    found = {}
    for command in commands:
        for program in programs.command_programs(command):
            if program in chosen:
                found[program] = chosen[program]
            else:
                location = shutil.which(program, path=parallel_task_manager.my_env["PATH"])
                if location is not None and programs.can_choose(program):
                    found[os.path.basename(program)] = os.path.abspath(location)
    if found and results_dir:
        run_state.update(results_dir, "programs", **found)
    return found


# The Python packages OrthoFinder needs (pyproject.toml): distribution name -> module name
PYTHON_PACKAGES = {
    "numpy": "numpy",
    "scipy": "scipy",
    "scikit-learn": "sklearn",
    "biopython": "Bio",
    "rich": "rich",
    "ete4": "ete4",
    "pandas": "pandas",
}


def record_python(results_dir):
    """
    The Python the run uses and its packages, in run_state.json: "python":
    {"executable", "version", "orthofinder": {"version", "path"}, "packages":
    {name: {"version", "path"} or null if it is not installed}}. The packages
    are found without importing them.
    """
    import importlib.metadata
    import importlib.util
    import platform
    import sys
    from .. import __version__

    packages = {}
    for dist, module in PYTHON_PACKAGES.items():
        try:
            spec = importlib.util.find_spec(module)
        except (ImportError, ValueError):
            spec = None
        if spec is None:
            packages[dist] = None
            continue
        try:
            version = importlib.metadata.version(dist)
        except importlib.metadata.PackageNotFoundError:
            version = None
        locations = list(spec.submodule_search_locations or [])
        path = locations[0] if locations else spec.origin
        packages[dist] = {"version": version, "path": os.path.realpath(path) if path else None}
    python = {
        "executable": os.path.abspath(sys.executable),   # (not resolved: e.g. a virtual environment's)
        "prefix": sys.prefix,                            # the environment (virtual or conda) or installation
        "version": platform.python_version(),
        "orthofinder": {"version": __version__, "path": os.path.realpath(__location__)},
        "packages": packages,
    }
    if results_dir:
        run_state.update(results_dir, "python", **python)
    return python


def record_system(results_dir):
    """The computer the run is on (utils.system_info), in run_state.json: "system"."""
    from ..utils import system_info
    system = system_info.system_info()
    if results_dir:
        run_state.update(results_dir, "system", **system)
    return system


def can_run_mcl(commands=None):
    command = "mcl -h"
    success, _, _ = run_test(lambda: (parallel_task_manager.can_run_command(command), command, None))
    if success:
        if commands is not None:
            commands.append(command)
        return True
    else:
        printer.print("ERROR: Cannot run MCL with the command \"%s\"" % command, style="error")
        program_caller.ProgramCaller.print_dependency_check_failure(command)
        printer.print("Please check MCL is installed and in the system path. See information above.\n", style="error")
        return False

def _method_test(result):
    """A ProgramCaller test's result as run_test takes it: True if the test is turned off in the config."""
    if result is True:
        return True, None, None
    success, stdout, stderr, cmd = result
    return success, cmd, None


def can_run_species_tree_program(prog_caller, program, commands=None):
    """Whether the species tree program (-ST) runs, on a few small gene trees."""
    if program not in prog_caller.list_species_tree_methods():
        # e.g. restored for a restart, from a config file no longer given
        printer.print('ERROR: The species tree program "%s" is not configured (program_type "species_tree" in '
                      'config.json or a file given with --config)\n' % program, style="error")
        return False
    d_deps_test = file_handler.handler.get_dependencies_check_dir()
    success, cmd, _ = run_test(lambda: _method_test(
        prog_caller.test_species_tree_method(program, d_deps_test, method_threads="2")))
    if cmd is None:   # not tested: "skip_check" in its config entry
        return True
    if success:
        if commands is not None:
            commands.append(cmd)
        return True
    printer.print("ERROR: Cannot run %s" % program, style="error")
    prog_caller.print_dependency_check_failure(cmd)
    printer.print("Please check %s is installed and that the executables are in the system path\n" % program, style="error")
    return False


def species_tree_program_to_check(options, q_assign):
    """The species tree program the analysis will run (-ST, or for --assign ASTRAL-Pro), or None."""
    if q_assign:
        return process_args.assign_species_tree_program(options)
    q_species_tree = not (options.q_stop_after_prepare or options.q_stop_after_mcl_groups
                          or options.q_start_from_trees or options.q_start_from_species_trees or options.species_tree_fn)
    if options.species_tree_program is not None and options.q_msa_trees and q_species_tree:
        return options.species_tree_program
    return None

def check_dependencies(options, user_specified_m, prog_caller, dir_for_temp_files, q_assign=False,
                       recorded_programs=None):
    """
    q_assign: an --assign analysis (also when continued with -b)
    recorded_programs: {name: path} of the programs the run being continued
        used (run_state.json), tried first
    """
    q_assign = q_assign or options.q_fast_add
    with run_state.stage(run_state.DEPENDENCY_CHECK):     # (its test commands noted as of it)
        _check_dependencies(options, user_specified_m, prog_caller, dir_for_temp_files, q_assign, recorded_programs)


def _check_dependencies(options, user_specified_m, prog_caller, dir_for_temp_files, q_assign, recorded_programs):
    choose_recorded(recorded_programs)
    commands = []     # those that ran: their programs are recorded in run_state.json
    printing.print_underline("Checking required programs are installed", True)
    if not user_specified_m:
        printer.print('Running with the recommended MSA tree inference by default. To revert to legacy method use "-M dendroblast".\n')
    if options.q_start_from_fasta or q_assign:
        d_deps_check = file_handler.handler.get_dependencies_check_dir()
        success, cmd, _ = run_test(lambda: _method_test(prog_caller.test_search_method(
            options.search_program, d_deps_check, scorematrix=options.score_matrix,
            gapopen=options.gapopen, gapextend=options.gapextend)))
        if not success:
            printer.print("\nERROR: Cannot run %s" % options.search_program, style="error")
            prog_caller.print_dependency_check_failure(cmd)
            printer.print("Please check %s is installed and that the executables are in the system path\n" % options.search_program,  style="error")
            util.fail()
        commands.append(cmd)
    if (options.q_start_from_fasta or options.q_start_from_blast) and not can_run_mcl(commands):
        util.fail()
    if not (options.q_stop_after_prepare or options.q_stop_after_mcl_groups or options.q_start_from_trees):
        if not can_run_ortholog_dependencies(dir_for_temp_files, 
                                            options.q_msa_trees, 
                                            options.q_phyldog, 
                                            options.q_stop_after_trees, 
                                            options.msa_program, 
                                            options.tree_program, 
                                            options.recon_method,
                                            prog_caller, 
                                            commands):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent ortholog inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.fail()
    species_tree_program = species_tree_program_to_check(options, q_assign)
    if species_tree_program is not None and not can_run_species_tree_program(prog_caller, species_tree_program, commands):
        util.fail()
    used = record_programs(commands, file_handler.handler.rd1)
    warn_changed_programs(recorded_programs, used)
    record_python(file_handler.handler.rd1)
    record_system(file_handler.handler.rd1)

def write_test_distances_file(test_fn):
    with open(test_fn, 'w', encoding="utf-8", newline="\n") as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
    return test_fn

def can_run_ortholog_dependencies(
        working_dir, 
        q_msa_gene_trees, 
        q_phyldog, 
        q_stop_after_trees, 
        msa_method, 
        tree_method,
        recon_method, 
        program_caller, 
        commands=None,
    ):
    """commands: if given, the test commands that ran are added to it"""
    if commands is None:
        commands = []
    d_deps_test = file_handler.handler.get_dependencies_check_dir()
    # FastME
    if not q_msa_gene_trees:
        test_fn = d_deps_test + "SimpleTest.phy"
        write_test_distances_file(test_fn)
        out_fn = d_deps_test + "SimpleTest.tre"
        cmd = "fastme -i %s -o %s" % (test_fn, out_fn)

        def fastme_test():
            if os.path.exists(out_fn):
                os.remove(out_fn)
            return parallel_task_manager.can_run_command(cmd, q_allow_stderr=False), cmd, None
        if not run_test(fastme_test)[0]:
            printer.print("ERROR: Cannot run fastme", style="error")
            program_caller.print_dependency_check_failure(cmd)
            printer.print("Please check FastME is installed and that the executables are in the system path.\n", style="error")
            return False
        commands.append(cmd)
    # DLCPar
    # if ("dlcpar" in recon_method) and not (qStopAfterTrees or qStopAfterAlignments):
    #     if not parallel_task_manager.can_run_command("dlcpar_search --version", qAllowStderr=False):
    #         printer.print("ERROR: Cannot run dlcpar_search", style="error")
    #         printer.print("Please check DLCpar is installed and that the executables are in the system path.\n", style="error")
    #         return False
        # if recon_method == "dlcpar_convergedsearch":
        #     capture = subprocess.Popen(
        #         "dlcpar_search --version", 
        #         shell=True, 
        #         stdout=subprocess.PIPE, 
        #         stderr=subprocess.PIPE, 
        #         env=parallel_task_manager.my_env
        #     )
        #     stdout = [x for x in capture.stdout]
        #     try:
        #         stdout = "".join([x.decode() for x in stdout])
        #     except (UnicodeDecodeError, AttributeError):
        #         stdout = "".join([x.encode() for x in stdout])
        #     version = stdout.split()[-1]
        #     tokens = list(map(int, version.split(".")))
        #     major, minor = tokens[:2]
        #     release = tokens[2] if len(tokens) > 2 else 0
        #     # require 1.0.1 or above            
        #     actual = (major, minor, release)
        #     required = [1,0,1]
        #     versionOK = True
        #     for r, a in zip(required, actual):
        #         if a > r:
        #             versionOK = True
        #             break
        #         elif a < r:
        #             versionOK = False
        #             break
        #         else:
        #             pass
        #             # need to check next level down
        #     if not versionOK:
        #         printer.print("ERROR: dlcpar_convergedsearch requires dlcpar_search version 1.0.1 or above", style="error")
        #         return False                   
    
    # FastTree & MAFFT
    if q_msa_gene_trees or q_phyldog:
        test_fn = msa_trees.write_test_file(d_deps_test)
        if msa_method is not None:
            success, cmd, _ = run_test(lambda: _method_test(program_caller.test_msa_method(msa_method, d_deps_test)))
            if not success:
                printer.print("ERROR: Cannot run MSA method '%s'" % msa_method, style="error")
                program_caller.print_dependency_check_failure(cmd)
                printer.print("Please check program is installed. If it is user-configured please check the configuration in the "
                      "orthofinder/config.json file\n", style="error")
                return False
            if cmd is not None:
                commands.append(cmd)
        if tree_method is not None:
            if q_msa_gene_trees:
                success, cmd, _ = run_test(lambda: _method_test(program_caller.test_tree_method(tree_method, d_deps_test)))
                if not success:
                   printer.print("ERROR: Cannot run tree method '%s'" % tree_method, style="error")
                   program_caller.print_dependency_check_failure(cmd)
                   printer.print("Please check program is installed. If it is user-configured please check the configuration in "
                         "the orthofinder/config.json file\n", style="error")
                   return False
                if cmd is not None:
                    commands.append(cmd)
            
    # if qPhyldog:
    #     if not parallel_task_manager.can_run_command("mpirun -np 1 phyldog", qAllowStderr=False):
    #         printer.print("ERROR: Cannot run mpirun -np 1 phyldog", style="error")
    #         printer.print("Please check phyldog is installed and that the executable is in the system path\n", style="error")
    #         return False
        
    return True    