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
from . import util
from ..utils.util import printer
try:
    from rich import print
except ImportError:
    ...
    
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


def setup_environment():

    os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised 

    my_env = os.environ.copy()
    # use orthofinder supplied executables by preference
    local_bin_dir = os.path.join(__location__, 'bin')
    bin_dirs = [
        "/opt/bin",
        "/usr/bin",
        "/usr/local/bin",
        os.path.expanduser("~/bin"),
        os.path.expanduser("~/.local/bin"),
        os.path.expanduser("~/local/bin"),
        local_bin_dir,
    ]
    for bin_dir in bin_dirs:
        my_env['PATH'] = bin_dir + os.pathsep + my_env['PATH']

    conda_prefix = my_env.get("CONDA_PREFIX")
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, "Scripts") if os.name == "nt" else os.path.join(conda_prefix, "bin")
        my_env["PATH"] = conda_bin + os.pathsep + my_env["PATH"]
    
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
        consoles = [get_console(), util.printer.console]
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
    if system in ["Linux", "Darwin"]:
        if mp.get_start_method(allow_none=True) != 'fork':
            mp.set_start_method('fork')
except RuntimeError as e:
    print(f"Multiprocessing context setting error on {system}: {e}")
    pass



def PrintTime(message):
    run_logging.RunLogger.message(message)
    util.printer.print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message), style="default")
    sys.stdout.flush()


def PrintNoNewLine(text):
    util.printer.print(text, end="")
    sys.stdout.flush()
    # sys.stdout.write(text)


class WorkerError(RuntimeError):
    """
    A child process or external command failed.

    The message carries the child's error (traceback, exit status or the
    command's output). Raising it, instead of calling Fail(), lets main()
    record the message in checkpoint.txt before OrthoFinder exits.
    """


def TerminateProcesses(processes, grace=5.0):
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


def WaitForExit(processes, what="child processes", warn_after=10.0, warn_interval=600.0):
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


def ParallelMap(function, args_list, nProcesses, progress=None):
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
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    args_list = list(args_list)
    if not args_list:
        return []
    pool = ProcessPoolExecutor(max_workers=max(1, min(nProcesses, len(args_list))))
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
        TerminateProcesses(list(getattr(pool, "_processes", {}).values()))
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    return results


def ManageQueue(runningProcesses, cmd_queue):
    """Manage a set of runningProcesses working through cmd_queue.
    If a process fails, stop the others straight away and raise WorkerError.
    Otherwise return when all work is complete.
    """
    try:
        while True:
            alive = False
            for proc in runningProcesses:
                if proc.is_alive():
                    alive = True
                elif proc.exitcode != 0:
                    raise WorkerError(
                        "%s (PID %s) exited with status %s"
                        % (proc.name, proc.pid, proc.exitcode)
                    )
            if not alive:
                return
            time.sleep(.1)
    except BaseException:
        TerminateProcesses(runningProcesses)
        cmd_queue.cancel_join_thread()
        raise

# not used
def RunCommand_Simple(command):
    subprocess.call(command, env=my_env, shell=True)


# How often to report an external program that is still running (seconds).
# There is deliberately no time limit: a program that is still running is
# never stopped for being slow (a large alignment or tree can take hours).
# Commands are only stopped when something has failed.
COMMAND_WARN_INTERVAL = 1800.0

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


def _process_group_cpu_seconds(pgid):
    """CPU seconds used so far by live processes in a process group (Linux only)."""
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
            with open("/proc/%s/stat" % entry) as f:
                stat = f.read()
        except OSError:
            continue
        # Fields after the "(comm)" field: state, ppid, pgrp, ..., utime, stime
        fields = stat.rsplit(")", 1)[-1].split()
        try:
            if int(fields[2]) == pgid:
                total_ticks += int(fields[11]) + int(fields[12])
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
    """Stop a command started by RunMonitoredCommand, including its child processes."""
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


def KillRunningCommands():
    """
    Stop every external command that is still running, and refuse to start
    new ones until ResetCommandAbort() (used when a command fails or on abort).
    """
    with _live_commands_lock:
        _commands_aborted.set()
        running = list(_live_commands)
    for popen in running:
        _kill_command(popen)


def ResetCommandAbort():
    """Allow commands to start again (call before starting a new batch)."""
    _commands_aborted.clear()


def RunMonitoredCommand(command, env, stdout=None, stderr=None,
                        warn_interval=None):
    """
    Run a shell command and return (returncode, stdout, stderr).

    While it runs, a warning is printed every warn_interval seconds with the
    elapsed time and the CPU time the command used in that interval, so a
    slow program (using CPU) can be told apart from a stuck one (not using
    CPU). The command is never stopped for taking long; only if another
    command fails or OrthoFinder itself is stopped (KillRunningCommands).
    """
    if warn_interval is None:
        warn_interval = COMMAND_WARN_INTERVAL
    use_group = hasattr(os, "killpg")

    # Start and register under the lock, so KillRunningCommands() cannot miss
    # a command that is starting at the same moment.
    with _live_commands_lock:
        if _commands_aborted.is_set():
            raise WorkerError("Not started because another command failed: %s" % command)
        popen = subprocess.Popen(
            command, env=env, shell=True,
            stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
            start_new_session=use_group,
        )
        _live_commands.add(popen)

    start = time.monotonic()
    cpu_prev = 0.0
    try:
        while True:
            try:
                # The timeout only sets how often we report; it never stops the command.
                out, err = popen.communicate(timeout=warn_interval)
                return popen.returncode, out, err
            except subprocess.TimeoutExpired:
                pass

            elapsed = time.monotonic() - start

            cpu_now = _process_group_cpu_seconds(popen.pid) if use_group else None
            if cpu_now is None:
                cpu_text = ""
            else:
                cpu_used = max(0.0, cpu_now - cpu_prev)
                cpu_prev = cpu_now
                cpu_text = (
                    "; CPU time used in the last %s: %s%s" % (
                        _format_duration(warn_interval), _format_duration(cpu_used),
                        " (not using CPU - it may be stuck, e.g. waiting on disk or network)"
                        if cpu_used < 1.0 else "",
                    )
                )
            short_cmd = command if len(command) <= 300 else command[:300] + " ..."
            print(
                "WARNING: external command still running after %s%s\n  %s"
                % (_format_duration(elapsed), cpu_text, short_cmd)
            )
    except BaseException:
        _kill_command(popen)
        raise
    finally:
        with _live_commands_lock:
            _live_commands.discard(popen)


def RunCommand(command, qPrintOnError=False, qPrintStderr=True, raise_on_error=False):
    """Run a single command"""
    capture = qPrintOnError or raise_on_error
    returncode, stdout, stderr = RunMonitoredCommand(
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
        elif qPrintStderr and len(stderr) > 0 and not util.stderr_exempt(stderr):
            print("\nWARNING: program called by OrthoFinder produced output to stderr")
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
    return returncode


def CanRunCommand(
    command,
    qAllowStderr=False,
    qPrint=True,
    qRequireStdout=True,
    qCheckReturnCode=False,
):
    if qPrint:
        PrintNoNewLine(f'Test can run "[orange3]{command.split()[0]}[/orange3]"')  # print without newline
    # communicate() rather than wait(): a program that prints more than the
    # pipe buffer would otherwise block forever.
    returncode, out, err = RunMonitoredCommand(
        command, my_env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    stdout = out.splitlines(True)
    stderr = err.splitlines(True)
    if qCheckReturnCode:
        return_code_check = returncode == 0
    else:
        return_code_check = True
    if (
        (len(stdout) > 0 or not qRequireStdout)
        and (qAllowStderr or len(stderr) == 0)
        and return_code_check
    ):
        if qPrint:
            util.printer.print(" - [bold green]ok")
        return True
    else:
        if qPrint:
            util.printer.print(" - [bold red]failed")
        if not return_code_check:
            util.printer.print("Returned a non-zero code: %d" % returncode, style="error")
        print("\nstdout:")
        for l in stdout:
            print(l)
        print("\nstderr:")
        for l in stderr:
            print(l)
        return False


q_print_first_traceback_0 = False


def Worker_RunCommands_And_Move(
    cmd_and_filename_queue,
    nProcesses,
    nToDo,
    qListOfLists,
    q_print_on_error,
    q_always_print_stderr,
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
    while True:
        try:
            i, command_fns_list = cmd_and_filename_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if (
                nDone >= 0
                and divmod(
                    nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000
                )[1]
                == 0
            ):
                PrintTime("Done %d of %d" % (nDone, nToDo))
            if not qListOfLists:
                command_fns_list = [command_fns_list]
            for command, fns in command_fns_list:
                if isinstance(command, types.FunctionType):
                    # This will block the process, but it is ok for trimming, it takes minimal time
                    fn = command
                    fn(*fns)
                else:
                    if not isinstance(command, str):
                        raise TypeError(f"Cannot run command: {command!r}")
                    else:
                        RunCommand(
                            command,
                            qPrintOnError=q_print_on_error,
                            qPrintStderr=q_always_print_stderr,
                            raise_on_error=True,
                        )
                        if fns != None:
                            actual, target = fns
                            if os.path.exists(actual):
                                os.rename(actual, target)
        except queue.Empty:
            return
        except BaseException:
            raise


q_print_first_traceback_1 = False


def Worker_RunMethod(result_queue, Function, args_queue):
    while True:
        try:
            args = args_queue.get(True, 0.1)
        except queue.Empty:
            return
        Function(*args)
        result_queue.put((None, "success"))


def ReportWorkerFailure(result_queue, function, *args):
    """Forward a worker's traceback to its parent without task-specific labels."""
    try:
        function(*args)
    except BaseException:
        result_queue.put(("error", traceback.format_exc()))
    finally:
        result_queue.put(None)

def RunMethodParallel(
        Function,
        args_queue,
        nProcesses,
        total_tasks=0,
        show_progress=True,
    ):

    result_queue = mp.Queue()
    runningProcesses = [
        mp.Process(target=ReportWorkerFailure,
                   args=(result_queue, Worker_RunMethod, result_queue, Function, args_queue))
        for i_ in range(nProcesses)
    ]
    ManageQueueNew(runningProcesses, total_tasks, nProcesses, result_queue, show_progress=show_progress)

def ManageQueueNew(
        runningProcesses,
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
    for proc in runningProcesses:
        proc.start()
    progressbar, task = util.get_progressbar(total_tasks, visible=show_progress)
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
                failed = [proc for proc in runningProcesses if proc.exitcode not in (None, 0)]
                if failed:
                    details = "; ".join(
                        f"{proc.name} (PID {proc.pid}) exited with status {proc.exitcode}"
                        for proc in failed
                    )
                    raise WorkerError(f"Worker exited without reporting a traceback: {details}")
                if all(proc.exitcode is not None for proc in runningProcesses):
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
        TerminateProcesses(runningProcesses)
        result_queue.cancel_join_thread()
        raise

    # All workers reported that they have finished: wait for them to exit.
    try:
        WaitForExit(runningProcesses, "workers", GRACE_PERIOD)
    except BaseException:
        TerminateProcesses(runningProcesses)
        raise
    if show_progress:
        progressbar.stop()
    try:
        result_queue.close()
        result_queue.join_thread()
    except Exception:
        pass
    



def _I_Spawn_Processes(message_to_spawner, message_to_PTM):
    """
    Args:
        message_queue - for passing messages that a new queue of tasks should be started (PTM -> I_Space_Processes) or that the tasks are complete
        cmds_queue - queue containing tasks that should be done
    Use:
        A process should be started as early as possible (while RAM usage is low) with this method as its target.
        This is now a separate process with low RAM usage.
        Each time some parallel work is required then the queue for that is placed in the message_queue by the PTM.
        _I_Spawn_Processes - will spawn parallel processes when instructed by the message_queue in the message_queue and get them
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
                #         PrintTime("Done %d of %d" % (n_done, n_to_do))
            # Version 2: launch n worker threads each executing a worker method that takes tasks from a queue
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
            wait(futures)
            message_to_PTM.put("Done")
            time.sleep(1)
        except queue.Empty:
            time.sleep(1)  # there wasn't anything this time, sleep then try again
    pass


class ParallelTaskManager_singleton:
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
            self.message_to_PTM = mp.Queue()
            # Orders/Messages:
            # None (PTM -> spawn_thread) - thread should return (i.e. exit)
            # 'Done' (spawn_thread -> PTM) - the cmds from the cmd queue have completed
            # Anything else = (nParallel, nTasks) (PTM -> spawn_thread) - cmds (nTasks of them) have been placed in the cmd queue,
            #   they should be executed using nParallel threads
            self.manager_process = mp.Process(
                target=_I_Spawn_Processes,
                args=(self.message_to_spawner, self.message_to_PTM),
            )
            self.manager_process.start()

    instance = None

    def __init__(self):
        if not ParallelTaskManager_singleton.instance:
            ParallelTaskManager_singleton.instance = (
                ParallelTaskManager_singleton.__Singleton()
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

    def Stop(self):
        """Warning, cannot be restarted"""
        self.instance.message_to_spawner.put(None)
        self.instance.manager_process.join()


# def RunParallelMethods(func, args_list, nProcesses):
#     """nProcesss - the number of processes to run in parallel
#     commands - list of lists of commands where the commands in the inner list are completed in order (the i_th won't run until
#     the i-1_th has finished).
#     """
#     ptm = ParallelTaskManager_singleton()
#     ptm.RunParallel(func, args_list, nProcesses)


def Success():
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    sys.exit()


def Fail():
    sys.stderr.flush()
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    printer.print(
        "ERROR: An error occurred, ***please review the error messages*** they may contain useful information about the problem.", style="error"
    )
    sys.exit(1)
