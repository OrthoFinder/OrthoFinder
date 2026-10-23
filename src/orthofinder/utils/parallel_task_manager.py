# -*- coding: utf-8 -*-
#
# Copyright 2017 David Emms
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
from . import logging as run_logging
import contextlib
import gc
import os
import sys
import platform
import time
import types
import signal
import threading
import datetime
import traceback
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, wait
from .. import __location__
from . import programs
from .printing import printer
try:
    from rich import print
except ImportError:
    ...
from . import printing
    
try:
    width = os.get_terminal_size().columns
except OSError as e:
    width = 80
try:
    import queue
except ImportError:
    import Queue as queue

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)


PROGRAM_DIRS_TEXT = (
    "the active conda environment, conda's base environment, ~/local/bin (where "
    "'make install_tools' installs them), PATH (system-wide installations included), "
    "then OrthoFinder's own bin directory (on Linux x86-64). The first copy of each "
    "program that works is used"
)


def bundled_binaries_run_here():
    """Whether the programs in OrthoFinder's own bin/ (built for Linux x86-64) can run on this computer."""
    return platform.system() == "Linux" and platform.machine().lower() in ("x86_64", "amd64")


def program_dirs(environ, location):
    """
    The directories searched for the programs OrthoFinder runs, in order of
    preference (those that exist, each once):
      1. the active conda environment (CONDA_PREFIX), if it is not base;
      2. conda's base environment (from CONDA_EXE, <base>/bin/conda, or
         MAMBA_ROOT_PREFIX): what the active environment does not have;
      3. ~/local/bin, where `make install_tools` installs them;
      4. PATH as the user has it: their other installations, system-wide
         ones included, in their order;
      5. OrthoFinder's own bin/ (and bin/mafft/bin), last: the programs there
         are prebuilt, so they may be missing or not run on this computer.
         Only on Linux x86-64, which they are built for.
    The dependency check uses the first copy of each program that works
    (programs.py), so one that does not run never hides a later one.
    """
    def conda_bin(prefix):
        return os.path.join(prefix, "Scripts" if os.name == "nt" else "bin")
    base = None
    if environ.get("CONDA_EXE"):
        base = os.path.dirname(os.path.dirname(environ["CONDA_EXE"]))
    elif environ.get("MAMBA_ROOT_PREFIX"):
        base = environ["MAMBA_ROOT_PREFIX"]
    dirs = []
    active = environ.get("CONDA_PREFIX")
    if active:
        dirs.append(conda_bin(active))
    if base:
        dirs.append(conda_bin(base))
    dirs.append(os.path.expanduser("~/local/bin"))
    dirs += environ.get("PATH", "").split(os.pathsep)
    if bundled_binaries_run_here():
        dirs += [os.path.join(location, "bin"), os.path.join(location, "bin", "mafft", "bin")]
    unique = []
    seen = set()
    for d in dirs:
        if d and os.path.isdir(d) and os.path.realpath(d) not in seen:
            seen.add(os.path.realpath(d))
            unique.append(d)
    return unique


def setup_environment():

    os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised
    os.environ["MKL_NUM_THREADS"] = "1"         # the same for numpy built with MKL (e.g. some conda channels)

    my_env = os.environ.copy()
    my_env["PATH"] = os.pathsep.join(program_dirs(my_env, __location__))
    return my_env

my_env = setup_environment()

# Fix LD_LIBRARY_PATH when using pyinstaller 
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''

def _reset_console_locks_in_child():
    """
    A process created by fork() inherits every lock in the state it was in.
    If another thread (e.g. a progress bar) was printing at that moment, the
    console lock is held by a thread that does not exist in the child, and the
    child's first print would block forever. Give the child fresh locks.
    """
    try:
        from rich import get_console
        consoles = [get_console(), printing.printer.console]
    except Exception:
        return
    for console in consoles:
        for name in ("_lock", "_record_buffer_lock"):
            if hasattr(console, name):
                setattr(console, name, threading.RLock())


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_reset_console_locks_in_child)


system = platform.system()
try:
    # (in OrthoFinder's own process: a process it started, e.g. with "spawn"
    # (ortholog_writers.split_ortholog_files), has its start method already)
    if system in ["Linux", "Darwin"] and mp.parent_process() is None:
        if mp.get_start_method(allow_none=True) != 'fork':
            mp.set_start_method('fork')
except RuntimeError as e:
    print(f"Multiprocessing context setting error on {system}: {e}")
    pass



_freeze_depth = 0
_freeze_lock = threading.Lock()


@contextlib.contextmanager
def frozen_for_fork():
    """
    While processes started by fork run (from before they start until they end): the
    objects there are left out of garbage collections (gc.freeze), the processes'
    and the parent's. A collection touches every object it looks at (e.g. the
    orthogroups and the gene names: GBs with 1000s of species), so a copy of each of
    their memory pages would be made: by each process that collects, and, if the
    parent collects, for the parent (the processes keep the old ones). The objects
    are still freed as soon as they are not used; only the cyclic collector leaves
    them until the processes have ended.
    """
    global _freeze_depth
    with _freeze_lock:
        if _freeze_depth == 0:
            gc.freeze()
        _freeze_depth += 1
    try:
        yield
    finally:
        with _freeze_lock:
            _freeze_depth -= 1
            if _freeze_depth == 0:
                gc.unfreeze()


def print_time(message):
    run_logging.RunLogger.message(message)
    printing.printer.print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message), style="default")
    sys.stdout.flush()


def print_no_new_line(text):
    printing.printer.print(text, end="")
    sys.stdout.flush()
    # sys.stdout.write(text)


class WorkerError(RuntimeError):
    """
    A child process or external command failed.

    The message carries the child's error (traceback, exit status or the
    command's output). Raising it, instead of calling Fail(), lets main()
    record the message in checkpoint.txt before OrthoFinder exits.
    """


def terminate_processes(processes, grace=5.0):
    """Stop all still-running child processes: terminate, then kill if needed."""
    processes = [p for p in processes if p is not None]
    for proc in processes:
        if proc.is_alive():
            proc.terminate()
    deadline = time.monotonic() + grace
    for proc in processes:
        proc.join(timeout=max(0.0, deadline - time.monotonic()))
    for proc in processes:
        if proc.is_alive():
            proc.kill()
            proc.join()


def wait_for_exit(processes, what="child processes", warn_after=10.0, warn_interval=600.0):
    """
    Wait for processes that have finished their work to exit.

    They are never stopped here: a process that is still running is only
    stopped when something has failed. A warning is printed if they take
    longer than warn_after seconds, and then every warn_interval seconds.
    Raises WorkerError if one exits with a non-zero status.
    """
    start = time.monotonic()
    next_warning = start + warn_after
    while True:
        alive = [p for p in processes if p.is_alive()]
        for proc in processes:
            if proc.exitcode not in (None, 0):
                raise WorkerError(
                    "%s (PID %s) exited with status %s" % (proc.name, proc.pid, proc.exitcode)
                )
        if not alive:
            return
        alive[0].join(timeout=0.5)
        now = time.monotonic()
        if now >= next_warning:
            print(
                "WARNING: waiting for %d %s to exit after finishing their work "
                "(%.0fs): %s" % (
                    len(alive), what, now - start,
                    ", ".join("PID %s" % p.pid for p in alive[:5]),
                )
            )
            next_warning = now + warn_interval


def parallel_map(function, args_list, n_processes, progress=None, initializer=None, initargs=()):
    """
    Run function(args) for each element of args_list in a process pool and
    return the results (in args_list order).

    If any call fails, or a worker is killed (e.g. out of memory), the other
    workers are stopped straight away and a WorkerError with the child's
    traceback is raised.

    progress, if given, is called as progress(n_completed): first with 0 once
    the workers have been started (start any progress bar then, not before,
    so no worker is forked while its refresh thread holds a lock), and after
    each completed task.

    initializer(*initargs), if given, is run once in each worker before its tasks:
    for data every task needs, which is then sent once per worker, not per task.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    args_list = list(args_list)
    if not args_list:
        return []
    pool = ProcessPoolExecutor(max_workers=max(1, min(n_processes, len(args_list))),
                               initializer=initializer, initargs=initargs)
    futures = {pool.submit(function, args): i for i, args in enumerate(args_list)}
    results = [None] * len(args_list)
    try:
        if progress is not None:
            progress(0)
        for n_done, future in enumerate(as_completed(futures), start=1):
            try:
                results[futures[future]] = future.result()
            except Exception as e:
                # The child's traceback is attached as the exception's cause.
                cause = e.__cause__
                detail = str(cause) if cause is not None else ""
                raise WorkerError(
                    "Worker failed: %s: %s\n%s" % (type(e).__name__, e, detail)
                ) from e
            if progress is not None:
                progress(n_done)
    except BaseException:
        for future in futures:
            future.cancel()
        # ProcessPoolExecutor has no public way to stop running workers.
        terminate_processes(list(getattr(pool, "_processes", {}).values()))
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    return results


# How often to check that a running external program is still using CPU (seconds).
# There is deliberately no time limit: a program that is still running is
# never stopped for being slow (a large alignment or tree can take hours), and
# running for a long time is not reported. Only a program that used no CPU over
# a whole interval is reported, as it may be stuck. Commands are only stopped
# when something has failed.
COMMAND_WARN_INTERVAL = 3600.0

_live_commands = set()
_live_commands_lock = threading.Lock()
# Set once a batch of commands has failed: no further command may start.
_commands_aborted = threading.Event()


def _format_duration(seconds):
    seconds = int(round(seconds))
    if seconds < 120:
        return "%ds" % seconds
    hours, minutes = divmod(seconds // 60, 60)
    return "%dh %02dmin" % (hours, minutes) if hours else "%dmin" % minutes


def _session_cpu_seconds(sid):
    """
    CPU seconds used so far by the processes in a session (Linux only).

    By session rather than process group, so a program that moves itself into
    its own process group (e.g. GNU timeout) is still counted. Each process's
    finished children (cutime, cstime) are included, so the total does not drop
    when a step of a multi-step program (e.g. MAFFT's script) exits.
    """
    if not os.path.isdir("/proc"):
        return None
    total_ticks = 0
    try:
        entries = os.listdir("/proc")
    except OSError:
        return None
    for entry in entries:
        if not entry.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % entry, encoding="utf-8") as f:
                stat = f.read()
        except OSError:
            continue
        # Fields after the "(comm)" field: state, ppid, pgrp, session, ...,
        # utime, stime, cutime, cstime
        fields = stat.rsplit(")", 1)[-1].split()
        try:
            if int(fields[3]) == sid:
                total_ticks += sum(int(f) for f in fields[11:15])
        except (IndexError, ValueError):
            continue
    return total_ticks / os.sysconf("SC_CLK_TCK")


def _wait_for_group_exit(popen, grace):
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        popen.poll()  # reap the shell so it does not linger as a zombie
        try:
            os.killpg(popen.pid, 0)
        except (ProcessLookupError, PermissionError):
            return True
        time.sleep(0.1)
    return False


def _kill_command(popen, grace=5.0):
    """Stop a command started by run_monitored_command, including its child processes."""
    if hasattr(os, "killpg"):
        # The command runs in its own process group (start_new_session), so this
        # reaches the real program, not just the shell that launched it.
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(popen.pid, sig)
            except (ProcessLookupError, PermissionError):
                break
            if _wait_for_group_exit(popen, grace):
                break
    elif popen.poll() is None:
        try:
            popen.kill()
        except OSError:
            pass
    try:
        popen.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass


def kill_running_commands():
    """
    Stop every external command that is still running, and refuse to start
    new ones until reset_command_abort() (used when a command fails or on abort).
    """
    with _live_commands_lock:
        _commands_aborted.set()
        running = list(_live_commands)
    for popen in running:
        _kill_command(popen)


def reset_command_abort():
    """Allow commands to start again (call before starting a new batch)."""
    _commands_aborted.clear()


def run_monitored_command(command, env, stdout=None, stderr=None,
                        warn_interval=None):
    """
    Run a shell command and return (returncode, stdout, stderr).

    Every warn_interval seconds the CPU time the command used in that interval
    is checked (Linux only). Taking long is normal and is not reported; a
    warning is printed only if it used no CPU in the interval, since it may be
    stuck. The command is never stopped for taking long; only if another
    command fails or OrthoFinder itself is stopped (kill_running_commands).
    """
    if warn_interval is None:
        warn_interval = COMMAND_WARN_INTERVAL
    use_group = hasattr(os, "killpg")

    if _commands_aborted.is_set():
        raise WorkerError("Not started because another command failed: %s" % command)
    from . import run_state
    run_state.note_command(command)      # (an example of each kind, in run_state.json)
    # Started outside the lock (under it, the threads' commands started one at a time:
    # ~3x slower for many short commands), then registered under it and the abort
    # checked again: kill_running_commands() sets the abort under the lock before it
    # lists the commands, so a command starting at the same moment is either in its
    # list or sees the abort here and is stopped.
    popen = subprocess.Popen(
        programs.with_paths(command), env=programs.with_env(env, command), shell=True,
        stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
        start_new_session=use_group,
    )
    with _live_commands_lock:
        _live_commands.add(popen)
        aborted = _commands_aborted.is_set()
    if aborted:
        _kill_command(popen)
        with _live_commands_lock:
            _live_commands.discard(popen)
        raise WorkerError("Stopped because another command failed: %s" % command)

    start = time.monotonic()
    cpu_prev = 0.0
    try:
        while True:
            try:
                # The timeout only sets how often we check it; it never stops the command.
                out, err = popen.communicate(timeout=warn_interval)
                return popen.returncode, out, err
            except subprocess.TimeoutExpired:
                pass

            # A program still using CPU is just slow, which is normal: say nothing.
            cpu_now = _session_cpu_seconds(popen.pid) if use_group else None
            if cpu_now is None:
                continue
            cpu_used = max(0.0, cpu_now - cpu_prev)
            cpu_prev = cpu_now
            if cpu_used >= 1.0:
                continue
            elapsed = time.monotonic() - start
            short_cmd = command if len(command) <= 300 else command[:300] + " ..."
            print(
                "WARNING: external command has used no CPU in the last %s (running for %s); "
                "it may be stuck, e.g. waiting on disk or network\n  %s"
                % (_format_duration(warn_interval), _format_duration(elapsed), short_cmd)
            )
    except BaseException:
        _kill_command(popen)
        raise
    finally:
        with _live_commands_lock:
            _live_commands.discard(popen)


def run_command(command, q_print_on_error=False, q_print_stderr=True, raise_on_error=False):
    """Run a single command"""
    capture = q_print_on_error or raise_on_error
    returncode, stdout, stderr = run_monitored_command(
        command,
        my_env,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.PIPE if capture else subprocess.DEVNULL,
    )
    if capture:
        if raise_on_error and returncode != 0:
            raise subprocess.CalledProcessError(
                returncode, command, output=stdout, stderr=stderr
            )
        if returncode != 0:
            print(
                (
                    "\nERROR: external program called by OrthoFinder returned an error code: %d"
                    % returncode
                )
            )
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
        elif q_print_stderr and len(stderr) > 0 and not printing.stderr_exempt(stderr):
            print("\nWARNING: program called by OrthoFinder produced output to stderr")
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
    return returncode


def can_run_command(
    command,
    q_allow_stderr=False,
    q_print=True,
    q_require_stdout=True,
    q_check_return_code=False,
):
    if q_print:
        print_no_new_line(f'Test can run "[orange3]{command.split()[0]}[/orange3]"')  # print without newline
    # communicate() rather than wait(): a program that prints more than the
    # pipe buffer would otherwise block forever.
    returncode, out, err = run_monitored_command(
        command, my_env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    stdout = out.splitlines(True)
    stderr = err.splitlines(True)
    if q_check_return_code:
        return_code_check = returncode == 0
    else:
        return_code_check = True
    if (
        (len(stdout) > 0 or not q_require_stdout)
        and (q_allow_stderr or len(stderr) == 0)
        and return_code_check
    ):
        if q_print:
            printing.printer.print(" - [bold green]ok")
        return True
    else:
        if q_print:
            printing.printer.print(" - [bold red]failed")
        if not return_code_check:
            printing.printer.print("Returned a non-zero code: %d" % returncode, style="error")
        print("\nstdout:")
        for l in stdout:
            print(l)
        print("\nstderr:")
        for l in stderr:
            print(l)
        return False


q_print_first_traceback_0 = False


q_print_first_traceback_1 = False


def worker_run_method(result_queue, function, args_queue):
    while True:
        try:
            args = args_queue.get(True, 0.1)
        except queue.Empty:
            return
        function(*args)
        result_queue.put((None, "success"))


def report_worker_failure(result_queue, function, *args):
    """Forward a worker's traceback to its parent without task-specific labels."""
    try:
        function(*args)
    except BaseException:
        result_queue.put(("error", traceback.format_exc()))
    finally:
        result_queue.put(None)

def run_method_parallel(
        function,
        args_queue,
        n_processes,
        total_tasks=0,
        show_progress=True,
    ):

    result_queue = mp.Queue()
    running_processes = [
        mp.Process(target=report_worker_failure,
                   args=(result_queue, worker_run_method, result_queue, function, args_queue))
        for i_ in range(n_processes)
    ]
    manage_queue_new(running_processes, total_tasks, n_processes, result_queue, show_progress=show_progress)

def frozen_while_running(function):
    """A function that starts processes by fork and waits for them, run frozen_for_fork."""
    import functools

    @functools.wraps(function)
    def frozen(*args, **kwargs):
        with frozen_for_fork():
            return function(*args, **kwargs)
    return frozen


@frozen_while_running
def manage_queue_new(
        running_processes,
        total_tasks,
        nprocess,
        result_queue,
        GRACE_PERIOD = 10.,
        STALL_TIMEOUT = 200.,
        show_progress = True,
    ):
    """
    Wait for workers, failing on reported errors or abnormal exits.

    STALL_TIMEOUT is the interval between warnings while no task completes,
    not a deadline: a single large task (e.g. sorting a big file) can
    legitimately take longer. Workers are only stopped when something fails.
    GRACE_PERIOD is how long to wait for finished workers to exit before
    printing a warning (they are still not stopped).
    """

    # Start the workers before the progress bar: its refresh thread holds the
    # console lock while drawing, and a process forked at that moment would
    # inherit the lock already held and could block on its first print.
    for proc in running_processes:
        proc.start()
    progressbar, task = printing.get_progressbar(total_tasks, visible=show_progress)
    update_cycle = 1
    if show_progress:
        progressbar.start()
    completed_tasks = 0
    active_workers = nprocess
    last_progress_time = time.time()
    last_warning_time = last_progress_time
    try:
        while completed_tasks < total_tasks or active_workers > 0:
            try:
                msg = result_queue.get(timeout=0.1)
            except queue.Empty:
                failed = [proc for proc in running_processes if proc.exitcode not in (None, 0)]
                if failed:
                    details = "; ".join(
                        f"{proc.name} (PID {proc.pid}) exited with status {proc.exitcode}"
                        for proc in failed
                    )
                    raise WorkerError(f"Worker exited without reporting a traceback: {details}")
                if all(proc.exitcode is not None for proc in running_processes):
                    raise WorkerError(
                        f"Workers exited before reporting completion "
                        f"({completed_tasks}/{total_tasks} tasks completed)."
                    )
                now = time.time()
                if now - last_warning_time > STALL_TIMEOUT:
                    print(
                        f"WARNING: No task has completed for {now - last_progress_time:.0f}s "
                        f"(completed {completed_tasks}/{total_tasks}); workers are still running."
                    )
                    last_warning_time = now
                continue

            if msg is None:
                active_workers -= 1
                continue

            if msg is False:
                raise WorkerError("Worker reported a fatal error without a traceback.")

            # Any worker can report its traceback without task-specific metadata.
            if isinstance(msg, tuple) and len(msg) == 2 and msg[0] == "error":
                raise WorkerError(f"Worker failed:\n{msg[1]}")

            if isinstance(msg, tuple) and len(msg) == 2 and msg[1] == "success":
                completed_tasks += 1
                if show_progress:
                    progressbar.update(task, advance=update_cycle)
                last_progress_time = time.time()
                last_warning_time = last_progress_time
                continue

            raise TypeError(f"Unexpected message from worker: {type(msg)} {msg!r}")

    except BaseException:
        # A failure: stop the other workers now instead of waiting for them.
        if show_progress:
            progressbar.stop()
        terminate_processes(running_processes)
        result_queue.cancel_join_thread()
        raise

    # All workers reported that they have finished: wait for them to exit.
    try:
        wait_for_exit(running_processes, "workers", GRACE_PERIOD)
    except BaseException:
        terminate_processes(running_processes)
        raise
    if show_progress:
        progressbar.stop()
    try:
        result_queue.close()
        result_queue.join_thread()
    except Exception:
        pass
    



def _i_spawn_processes(message_to_spawner, message_to_ptm):
    """
    Args:
        message_queue - for passing messages that a new queue of tasks should be started (PTM -> I_Space_Processes) or that the tasks are complete
        cmds_queue - queue containing tasks that should be done
    Use:
        A process should be started as early as possible (while RAM usage is low) with this method as its target.
        This is now a separate process with low RAM usage.
        Each time some parallel work is required then the queue for that is placed in the message_queue by the PTM.
        _i_spawn_processes - will spawn parallel processes when instructed by the message_queue in the message_queue and get them
        working on the queue. When the queue is empty it will wait for the next one. It can receive a special signal to exit - the None
        object
    """
    while True:
        try:
            # peak in qoq - it is the only method that tried to remove things from the queue
            message = message_to_spawner.get(timeout=0.1)
            if message is None:
                # Respond to request to terminate
                return
            # In which case, thread has been informed that there are tasks in the queue.
            func, args_list, n_parallel = message
            futures = []
            n_to_do = len(args_list)
            # Version 1: n worker threads for executing the method and a list of N arguments for calling the method
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
                # for i, _ in as_completed(futures):
                #     n_done = i+1
                #     if n_done >= 0 and divmod(n_done, 10 if n_done <= 200 else 100 if n_done <= 2000 else 1000)[1] == 0:
                #         print_time("Done %d of %d" % (n_done, n_to_do))
            # Version 2: launch n worker threads each executing a worker method that takes tasks from a queue
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
            wait(futures)
            message_to_ptm.put("Done")
            time.sleep(1)
        except queue.Empty:
            time.sleep(1)  # there wasn't anything this time, sleep then try again
    pass


class ParallelTaskManagerSingleton:
    """
    Creating new process requires forking parent process and can lea to very high RAM usage. One way to mitigate this is
    to create the pool of processes as early in execution as possible so that the memory footprint is low. The
    ParallelTaskManager takes care of that, and can be used by calling `RunParallelOrderedCommandLists` above.
    Apr 2023 Update:
    When running external programs there is no need to use multiprocessing, multithreading is sufficient since new process
    will be created anyway, so the SIL is no longer an issue.
    """

    class __Singleton(object):
        def __init__(self):
            """Implementation:
            Allocate a thread that will perform all the tasks
            Communicate with it using a queue.
            When provided with a list of commands it should fire up some workers and get them to run the commands and then exit.
            An alternative would be they should always stay alive - but then they could die for some reason? And I'd have to check how many there are.
            """
            self.message_to_spawner = mp.Queue()
            self.message_to_ptm = mp.Queue()
            # Orders/Messages:
            # None (PTM -> spawn_thread) - thread should return (i.e. exit)
            # 'Done' (spawn_thread -> PTM) - the cmds from the cmd queue have completed
            # Anything else = (nParallel, nTasks) (PTM -> spawn_thread) - cmds (nTasks of them) have been placed in the cmd queue,
            #   they should be executed using nParallel threads
            self.manager_process = mp.Process(
                target=_i_spawn_processes,
                args=(self.message_to_spawner, self.message_to_ptm),
            )
            self.manager_process.start()

    instance = None

    def __init__(self):
        if not ParallelTaskManagerSingleton.instance:
            ParallelTaskManagerSingleton.instance = (
                ParallelTaskManagerSingleton.__Singleton()
            )

    # def RunParallel(self, func, args_list, nParallel):
    #     """
    #     Args:
    #         cmd_list - list of commands or list of lists of commands (in which elements in inner list must be run in order)
    #         nParallel - number of parallel threads to use
    #         qShell - should the tasks be run in a shell
    #     """
    #     self.instance.message_to_spawner.put((func, args_list, nParallel))
    #     while True:
    #         try:
    #             signal = self.instance.message_to_PTM.get()
    #             if signal == "Done":
    #                 return
    #         except queue.Empty:
    #             pass
    #         time.sleep(1)

    def stop(self):
        """Warning, cannot be restarted"""
        self.instance.message_to_spawner.put(None)
        self.instance.manager_process.join()


# def RunParallelMethods(func, args_list, nProcesses):
#     """nProcesss - the number of processes to run in parallel
#     commands - list of lists of commands where the commands in the inner list are completed in order (the i_th won't run until
#     the i-1_th has finished).
#     """
#     ptm = ParallelTaskManagerSingleton()
#     ptm.RunParallel(func, args_list, nProcesses)


def success():
    ptm = ParallelTaskManagerSingleton()
    ptm.stop()
    sys.exit()


ISSUES_URL = "https://github.com/OrthoFinder/OrthoFinder/issues"


def _results_dir():
    """The results directory of the run, if it has one yet (else None)."""
    try:
        from . import file_handler      # (not at the top: file_handler imports this module)
        d = file_handler.handler.rd1
    except Exception:
        return None
    return d if d and os.path.isdir(d) else None


def fail():
    """Stop OrthoFinder after an error (already printed): what to do next."""
    sys.stderr.flush()
    ptm = ParallelTaskManagerSingleton()
    ptm.stop()
    printer.print("\nOrthoFinder stopped because of the error above.", style="error")
    results_dir = _results_dir()
    printer.print("If the message does not say how to fix it, search for it at %s, or open an issue "
                  "there with the command you ran and the error%s" % (ISSUES_URL, ", and attach this file:" if results_dir else "."))
    if results_dir:
        printer.print("    %s" % os.path.join(results_dir, "Log.txt"))
    sys.exit(1)
