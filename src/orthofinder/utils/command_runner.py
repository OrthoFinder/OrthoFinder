"""Running external commands in parallel, sharing the CPU cores between them."""

import collections
import concurrent.futures
import os
import re
import shlex
import subprocess
import threading
import types

from . import parallel_task_manager, printing

try:
    from rich import print
except ImportError:
    ...


TOTAL_CORES = os.cpu_count() or 8


class CoreTokens(object):
    """
    Share TOTAL_CORES between concurrently running external commands.

    A command takes all the cores it needs at once, and requests are served
    in arrival order. Taking cores one at a time can deadlock: e.g. with 8
    cores and three commands needing 4 each, they can end up holding 3, 3 and
    2 and all wait forever.
    """

    def __init__(self, total):
        self.total = total
        self.available = total
        self.waiting = collections.deque()
        self.cond = threading.Condition()

    def acquire(self, n):
        n = max(1, min(int(n), self.total))
        ticket = object()
        with self.cond:
            self.waiting.append(ticket)
            try:
                while self.waiting[0] is not ticket or self.available < n:
                    self.cond.wait()
            except BaseException:
                self.waiting.remove(ticket)
                self.cond.notify_all()
                raise
            self.waiting.popleft()
            self.available -= n
            self.cond.notify_all()
        return n

    def release(self, n):
        with self.cond:
            self.available += n
            self.cond.notify_all()


TOKENS = CoreTokens(TOTAL_CORES)


_METHODTHREAD_RE = re.compile(r"(?<![A-Za-z0-9_])(METHODTHREADS?|METHODTHREAD)(?![A-Za-z0-9_])")


def fill_method_threads(command, method_threads):
    """
    (threads, command with its METHODTHREAD placeholder filled in): the
    threads one command is run with, -mt capped at the number of cores;
    DIAMOND always gets 1 (many of its commands run at once), as does a
    command without the placeholder.
    """
    try:
        prog = os.path.basename(shlex.split(command)[0]).lower()
    except Exception:
        prog = ""
    if _METHODTHREAD_RE.search(command) and not prog.startswith("diamond"):  # diamond makedb / blastp
        threads = max(1, min(int(method_threads), TOTAL_CORES))
    else:
        threads = 1
    return threads, _METHODTHREAD_RE.sub(str(threads), command)


# The thread-count options of the search programs, by program.
_SEARCH_THREAD_OPTIONS = (
    (re.compile(r"^diamond$"), ("-p", "--threads")),
    (re.compile(r"^mmseqs$"), ("--threads",)),
    (re.compile(r"^(t?blast[npx]|psiblast|deltablast)$"), ("-num_threads",)),
)


def reset_search_threads(command):
    """
    A saved search command with its thread counts set back to METHODTHREAD.

    Commands saved by -op get the -mt thread count (or more, when edited to
    run on a cluster); when OrthoFinder runs them itself (a -b restart) the
    count is filled in again when each command runs (fill_method_threads): 1
    by default, as many commands run at once (-t). Only the known thread
    options of DIAMOND, MMseqs2 and BLAST+ are changed.
    """
    parts = re.split(r"(;|&&|\|\||\|)", command)
    for i in range(0, len(parts), 2):     # the commands, not the separators
        words = parts[i].split()
        if not words:
            continue
        program = os.path.basename(words[0])
        for program_re, options in _SEARCH_THREAD_OPTIONS:
            if program_re.match(program):
                for option in options:
                    parts[i] = re.sub(r"(?<!\S)(%s)(\s+|=)\d+(?!\S)" % re.escape(option),
                                      r"\g<1>\g<2>METHODTHREAD", parts[i])
    return "".join(parts)


def run_command(command, method_threads, dynamic_threads=False, q_print_on_error=False, q_print_stderr=True, raise_on_error=False):
    """Run a single command with token gating."""
    
    # threads_needed = threads_from_cmd(command, method_threads)

    # if threads_needed > TOTAL_CORES:
    #     # if qPrintOnError:
    #     #     print(f"Capping threads from {threads_needed} to {TOTAL_CORES} for {prog} to avoid oversubscription.")
    #     threads_needed = TOTAL_CORES
    # if threads_needed < 1:
    #     threads_needed = 1

    threads_needed, command = fill_method_threads(command, method_threads)

    acquired = 0
    try:
        acquired = TOKENS.acquire(threads_needed)

        env = dict(parallel_task_manager.my_env) if hasattr(parallel_task_manager, "my_env") else os.environ.copy()
        if threads_needed > 1:
            env["OMP_NUM_THREADS"] = str(threads_needed)
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
        env.setdefault("MKL_NUM_THREADS", "1")
        env.setdefault("OMP_NESTED", "0")

        out = subprocess.PIPE if q_print_on_error or raise_on_error else subprocess.DEVNULL
        err = subprocess.PIPE if (q_print_on_error and q_print_stderr) or raise_on_error else subprocess.DEVNULL

        returncode, stdout, stderr = parallel_task_manager.run_monitored_command(
            command, env, stdout=out, stderr=err
        )

        if q_print_on_error or raise_on_error:
            if raise_on_error and returncode != 0:
                raise subprocess.CalledProcessError(returncode, command, output=stdout, stderr=stderr)
            if returncode != 0:
                print(f"\nERROR: external program returned code {returncode}")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout}")
                print(f"stderr:\n{stderr}")
            elif q_print_stderr and len(stderr) > 0 and not printing.stderr_exempt(stderr):
                print("\nWARNING: program produced output to stderr")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout}")
                print(f"stderr:\n{stderr}")
        return returncode
    finally:
        if acquired:
            TOKENS.release(acquired)


def run_parallel_commands(
    n_processes,
    commands,
    q_list_of_list,
    method_threads=None,
    method_threads_large=None,
    method_threads_small=None,
    threshold=None,
    cmd_order="descending",
    tasksize=None,
    q_trim=False,
    q_print_on_error=False,
    q_always_print_stderr=False,
    dynamic_threads=False,
    on_success=None,
):
    """
    on_success: called in this process with each command that succeeds (exit 0),
    e.g. to record finished searches; an error in it is reported, not raised.
    """
    if q_list_of_list:
        commands_and_no_filenames = [
            [(cmd, None) for cmd in cmd_list] for cmd_list in commands
        ]
    else:
        commands_and_no_filenames = [(cmd, None) for cmd in commands]
    run_parallel_commands_and_move_results_file(
        n_processes,
        commands_and_no_filenames,
        q_list_of_list,
        method_threads,
        method_threads_large,
        method_threads_small,
        threshold,
        cmd_order,
        tasksize,
        q_trim,
        q_print_on_error,
        q_always_print_stderr,
        dynamic_threads,
        on_success=on_success,
    )


def _run_parallel_commands_and_move_results_file(
    n_processes,
    commands_and_filenames,
    q_list_of_list,
    method_threads=None,
    method_threads_large=None,
    method_threads_small=None,
    threshold=None,
    cmd_order="descending",
    tasksize=None,
    q_trim=False,
    q_print_on_error=False,
    q_always_print_stderr=False,
    dynamic_threads=False,
    on_success=None,
):
    """
    Calls the commands in parallel and if required moves the results file to the required new filename
    Args:
        nProcess - the number of parallel process to use
        commands_and_filenames : tuple (cmd, actual_target_fn) where actual_target_fn = None if no move is required
                                 and actual_target_fn = (actual_fn, target_fn) is actual_fn is produced by cmd and this
                                 file should be moved to target_fn
        actual_target_fn - None if the cmd will save the results file to outfilename_proposed
                           otherwise (actual_fn, outfilename_proposed)
        qListOfList - if False then commands_and_filenames is a list of (cmd, actual_target_fn) tuples
                      if True then commands_and_filenames is a list of lists of (cmd, actual_target_fn) tuples where the elements
                      of the inner list need to be run in the order they appear.
        q_print_on_error - If error code returend print stdout & stederr
    """
    total_commands = len(commands_and_filenames)
    if method_threads is None:
        method_threads = "1"


    n_processes = max(1, min((n_processes, len(commands_and_filenames), TOTAL_CORES * 4)))

    progressbar, task = printing.get_progressbar(total_commands)
    progressbar.start()
    update_cycle = 1 #10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_processes) as executor:
        futures = {}
        for cmd_unit in commands_and_filenames:
            if cmd_unit is None:
                continue
            fut = executor.submit(
                worker_run_commands_and_move,
                cmd_unit,
                method_threads,
                q_list_of_list,
                q_print_on_error,
                q_always_print_stderr,
                dynamic_threads=dynamic_threads
            )
            futures[fut] = cmd_unit

        recorded = set()

        def record_success(future):
            recorded.add(future)
            try:
                on_success(futures[future][0] if not q_list_of_list else futures[future])
            except Exception as e:      # recording must not stop the run
                print("WARNING: could not record a finished command: %s" % e)

        try:
            for i, future in enumerate(concurrent.futures.as_completed(futures)):
                # The first failed command is raised straight away (with
                # its stdout/stderr, which main() writes to checkpoint.txt).
                result = future.result()
                if result != 0 and q_print_on_error:
                    print(f"ERROR occurred with command: {futures[future]}")
                if result == 0 and on_success is not None:
                    record_success(future)
                if (i + 1) % update_cycle == 0:
                    progressbar.update(task, advance=update_cycle)
        except BaseException:
            # A command failed or we were interrupted: don't start queued
            # commands and stop running ones, otherwise leaving the pool
            # would wait for all of them to finish.
            for fut in futures:
                fut.cancel()
            parallel_task_manager.kill_running_commands()
            progressbar.stop()
            if on_success is not None:
                # Commands that succeeded but were not yet recorded (e.g. they
                # finished while the batch was being stopped) are recorded too,
                # so a restart reuses their results rather than running them again.
                # Leaving the pool waits for the running ones anyway.
                concurrent.futures.wait(futures)
                for fut in futures:
                    if (fut not in recorded and not fut.cancelled()
                            and fut.exception() is None and fut.result() == 0):
                        record_success(fut)
            raise
    progressbar.stop()


def run_parallel_commands_and_move_results_file(*args, **kwargs):
    """
    See _run_parallel_commands_and_move_results_file. If a command fails, the other
    commands of the batch are stopped and the failure is raised. The abort
    flag that stops them is cleared once the whole batch has finished, so it
    never affects later commands.
    """
    parallel_task_manager.reset_command_abort()
    try:
        return _run_parallel_commands_and_move_results_file(*args, **kwargs)
    finally:
        parallel_task_manager.reset_command_abort()


q_print_first_traceback_0 = False


def worker_run_commands_and_move(
    command_fns_list, 
    method_threads, 
    q_list_of_lists, 
    q_print_on_error, 
    q_always_print_stderr, 
    dynamic_threads,
):
    """
    Continuously takes commands that need to be run from the cmd_and_filename_queue until the queue is empty. If required, moves
    the output filename produced by the cmd to a specified filename. The elements of the queue can be single cmd_filename tuples
    or an ordered list of tuples that must be run in the provided order.

    Args:
        cmd_and_filename_queue - queue containing (cmd, actual_target_fn) tuples (if qListOfLists is False) or a list of such
            tuples (if qListOfLists is True). Alternatively, 'cmd' can be a python fn and actual_target_fn the fn to call it on.
        nProcesses - the number of processes that are working on the queue.
        nToDo - The total number of elements in the original queue
        qListOfLists - Boolean, whether each element of the queue corresponds to a single command or a list of ordered commands
        qShell - Boolean, should a shell be used to run the command.

    Implementation:
        nProcesses and nToDo are used to print out the progress.
    """
    # while True:
    try:
        # i, command_fns_list = cmd_and_filename_queue.get(True, 1)

        # nDone = i - nProcesses + 1
        # if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
        #     util.print_time("Done %d of %d" % (nDone, nToDo))
        if not q_list_of_lists:
            command_fns_list = [command_fns_list]

        return_code = 1
        for command, fns in command_fns_list:
            if isinstance(command, types.FunctionType):
                # This will block the process, but it is ok for trimming, it takes minimal time
                fn = command
                fn(*fns)
                return_code = 0
            else:
                if not isinstance(command, str):
                    raise TypeError(f"Cannot run command: {command!r}")
                else:
                    return_code = run_command(
                        command,
                        method_threads,
                        q_print_on_error=q_print_on_error,
                        q_print_stderr=q_always_print_stderr,
                        raise_on_error=True,
                    )
                    if fns != None:
                        actual, target = fns
                        if os.path.exists(actual):
                            os.rename(actual, target)
        return return_code
    # except queue.Empty:
    #     return
    except BaseException:
        raise
