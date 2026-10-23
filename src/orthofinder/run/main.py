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
import signal
import copy  # Y
import time  # Y
import os.path  # Y
from ..utils import parallel_task_manager, file_handler, util, fasta_processor, file_io, run_state, orthoxml
from ..mcl_orthogroups import infer_orthogroups, orthogroups_set
from ..mcl_orthogroups import assign as assign
from orthofinder.utils import logging as run_logging


from ..mcl_orthogroups import mcl
from ..phylo import tree, tree_pages
from ..orthology import reconciliation
from ..gene_tree import infer_trees
from . import process_args, check_dependencies, previous_run
from .context import RunContext
from . import assign_run, finish, start, workflows
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
sys.setrecursionlimit(10**6)

# Printing a name the console's encoding cannot hold (e.g. a Windows code page, or
# output redirected to a file) shows "?" for the character rather than failing
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except (AttributeError, ValueError):   # (not a text stream, e.g. replaced by a test)
        pass

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)



"""
OrthoFinder
-------------------------------------------------------------------------------
"""


class Terminated(KeyboardInterrupt):
    """SIGTERM (e.g. the time limit of a job): the run ends as after Ctrl-C."""


def _sigterm_handler(main_pid):
    def handler(signum, frame):
        if os.getpid() != main_pid:
            # a child process (it inherited the handler): end as by default
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            os.kill(os.getpid(), signal.SIGTERM)
            return
        raise Terminated()
    return handler


def main(args=None):
    try:
        previous_sigterm = signal.signal(signal.SIGTERM, _sigterm_handler(os.getpid()))
    except ValueError:   # not the main thread
        previous_sigterm = None
    file_handler.handler.reset()
    file_io.clear_caches()   # main() may run more than once in a process (e.g. tests)
    t_start = time.perf_counter()
    ctx = RunContext()  # ctx.current_step: "Initialisation"
    exit_code = 0
    try:
        if args is None:
            args = sys.argv[1:]
        input_args = args.copy()
        # Create PTM right at start
        ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        # prog_caller = get_program_caller()
        start.start_run(ctx, args, input_args)
        options = ctx.options

        if options.q_start_from_blast and options.q_start_from_fasta:
            run_workflow = workflows.add_species
        elif options.q_start_from_fasta:
            run_workflow = workflows.from_fasta
        elif options.q_start_from_blast and not ctx.assign_resume:
            run_workflow = workflows.from_blast
        elif options.q_start_from_groups:
            run_workflow = workflows.from_groups
        elif options.q_start_from_trees:
            run_workflow = workflows.from_gene_trees
        elif options.q_start_from_species_trees:
            run_workflow = workflows.from_species_trees
        elif options.q_fast_add or ctx.assign_resume:
            run_workflow = assign_run.run_assign
        else:
            raise NotImplementedError
        species_info_obj, species_names_dict = run_workflow(ctx)
        finish.process_output_files(ctx, species_info_obj, species_names_dict)
        ctx.log.info("OrthoFinder run completed in %.2f seconds", time.perf_counter() - t_start)

    except Exception as e:
        ctx.current_step = run_logging.RunLogger.failed_stage(ctx.current_step)
        exit_code = 1
        if ctx.log is not None:
            ctx.log.log("ERROR: Step failed",
                    step=ctx.current_step, level="ERROR", exc_info=True)
            for stream in ("stdout", "stderr"):
                output = getattr(e, stream, None)
                if output:
                    if isinstance(output, bytes):
                        output = output.decode("utf-8", errors="replace")
                    ctx.log.log("%s:\n%s", stream, output,
                            step=ctx.current_step, level="ERROR")
                    printer.console.print(f"{stream}:\n{output}", markup=False, highlight=False)
        import traceback
        traceback.print_exception(type(e), e, e.__traceback__)
        finish.clean_up_at_exit(ctx.options, q_failed=True)
        # ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()
        sys.exit(1)

    except KeyboardInterrupt as e:
        ctx.current_step = run_logging.RunLogger.failed_stage(ctx.current_step)
        exit_code = 1
        terminated = isinstance(e, Terminated)
        if ctx.log is not None:
            ctx.log.step(ctx.current_step, "Terminated (SIGTERM)" if terminated else "Interrupted by user",
                         level="WARNING")
        printer.print("\nProgram terminated (SIGTERM)." if terminated else "\nProgram terminated by user.",
                      style="error")
        finish.clean_up_at_exit(ctx.options, q_failed=True)
        sys.exit(1)

    except SystemExit as e:
        ctx.current_step = run_logging.RunLogger.failed_stage(ctx.current_step)
        exit_code = e.code
        q_failed = not (e.code is None or e.code == 0)
        if ctx.log is not None:
            if not q_failed:
                ctx.log.step(ctx.current_step, "Workflow stopped normally before the full analysis finished")
            else:
                ctx.log.log("ERROR: Workflow exited with status %s", e.code,
                        step=ctx.current_step, level="ERROR", exc_info=True)
        # a stop option (-ogt, -ost, ...) or a failure: the MCL orthogroups, no empty directories
        finish.clean_up_at_exit(ctx.options, q_failed=q_failed)
        raise

    finally:
        if previous_sigterm is not None:
            signal.signal(signal.SIGTERM, previous_sigterm)
        run_logging.RunLogger.set_logger(None)
        if ctx.log is not None:
            ctx.log.close()
        # ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()
        end = time.perf_counter()
        time_elapsed = end - t_start
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
        if exit_code not in (None, 0):   # (stopped by an error: not "finished")
            printer.print(f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] stopped after {time_elapsed:5f}s", end="\n" * 2)
        else:
            printer.print(
                f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] finished in ", end=""
            )
            printer.print(f"[green]{time_elapsed:5f}[/green]s", end="\n" * 2)
        file_handler.handler.reset()
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
