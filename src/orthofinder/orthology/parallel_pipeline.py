# -*- coding: utf-8 -*-

import os
import time
import sys
import resource
import traceback
import warnings

import queue
import multiprocessing as mp

try:
    from rich import print
except ImportError:
    ...

from ..utils import parallel_task_manager
from ..utils import printing
from . import stats
from .ortholog_writers import non_hog_writer_process, ordered_hog_writer_process, partition_non_hog_result


def get_n_writer_processes(
        nspecies,
        n_processes,
        fewer_open_files=False,
        max_writers=12,
    ):
    """Choose the number of non-HOG output writer processes."""
    nspecies = max(1, int(nspecies))
    n_processes = max(1, int(n_processes))

    species_per_writer = 96 if fewer_open_files else 48
    n_writers = (nspecies + species_per_writer - 1) // species_per_writer
    writer_cpu_cap = max(1, n_processes // 2)

    return max(
        1,
        min(
            int(n_writers),
            int(max_writers),
            writer_cpu_cap,
        )
    )


def worker_run_orthologs_method_pipeline(
        tree_analyser,
        nspecies,
        args_queue,
        hog_queue,
        non_hog_queues,
        progress_queue,
        fewer_open_files,
        need_ortholog_output,
        n_orthologs_cache=100,
        write_hog_tree=False,
        fix_files=False,
    ):
    """
    Analyse OGs and route HOG and non-HOG output independently.

    Ortholog counts are accumulated locally and sent once with
    "worker_done". Sending a dense nspecies x nspecies count object per OG
    is ~40 MB of pickled data per OG at 1,000 species.
    """
    n_writers = len(non_hog_queues)
    worker_pid = os.getpid()
    n_ortho_acc = stats.NOrthoSp(nspecies)

    while True:
        try:
            iog = args_queue.get(True, 0.1)

            if iog is None:
                break

            progress_queue.put(("start", worker_pid, iog))
            result = tree_analyser.analyse_tree(iog, n_ortho_acc=n_ortho_acc)

            if result is None:
                hog_queue.put(("skip", iog))
                progress_queue.put(("skip", iog))
                continue

            # Analysis is finished. Any further wait is on the output queues.
            progress_queue.put(("output", iog))

            hog_queue.put((
                "result",
                iog,
                result.get("cached_hogs", []),
            ))

            payloads = partition_non_hog_result(
                result,
                n_writers,
                tree_analyser,
                fewer_open_files,
                need_ortholog_output,
            )

            for owner, payload in payloads.items():
                non_hog_queues[owner].put(("result", payload))

            progress_queue.put(("result", iog))

        except queue.Empty:
            continue

        except Exception:
            tb = traceback.format_exc()
            # Report directly to the parent. A failed writer may have left a
            # full output queue, so even sending it an error could block here.
            progress_queue.put(("error", None, tb))
            return

    hog_queue.put(None)
    for writer_queue in non_hog_queues:
        writer_queue.put(None)
    progress_queue.put(("worker_done", worker_pid, n_ortho_acc))


def run_orthologs_parallel_pipeline(
        tree_analyser,
        nspecies,
        args_queue,
        n_processes,
        total_tasks,
        fewer_open_files,
        output_writer,
        iogs_ordered,
        n_orthologs_cache=100,
        write_hog_tree=False,
        fix_files=False,
        fd_limit=None,
        GRACE_PERIOD=10.0,
        STALL_TIMEOUT=120.0,
        writer_queue_size=None,
        n_writer_processes=None,
        spill_dir=None,
    ):
    """
    Parallel tree analysis with two output paths.

    HOG rows are committed in deterministic OG order. Other output files are
    written immediately by exclusive file owners and can be sorted afterward.

    STALL_TIMEOUT is the interval between long-wait warnings, not a deadline
    for a tree or an output batch. Silence alone cannot distinguish slow work
    from a hang. Child process failures are checked independently.
    """
    if fd_limit is not None:
        if sys.platform.startswith("linux") or sys.platform == "darwin":
            set_file_descriptor_limit(fd_limit)
        else:
            warnings.warn(
                "File descriptor limit adjustment is not supported on %s." %
                sys.platform
            )

    if n_writer_processes is None:
        n_writer_processes = get_n_writer_processes(
            nspecies,
            n_processes,
            fewer_open_files=fewer_open_files,
        )
    else:
        n_writer_processes = max(
            1,
            min(int(n_writer_processes), max(1, n_processes))
        )

    # print(
    #     "Output pipeline: %d analysis workers, 1 ordered HOG writer, "
    #     "%d non-HOG writer process%s (%d species, %s mode)." % (
    #         nProcesses,
    #         n_writer_processes,
    #         "" if n_writer_processes == 1 else "es",
    #         nspecies,
    #         "compact" if fewer_open_files else "pairwise",
    #     )
    # )

    if writer_queue_size is None:
        writer_queue_size = max(2 * n_processes, 16)

    hog_queue = mp.Queue(maxsize=max(4 * n_processes, 32))
    non_hog_queues = [
        mp.Queue(maxsize=writer_queue_size)
        for _ in range(n_writer_processes)
    ]
    progress_queue = mp.Queue(maxsize=max(4 * n_processes, 32))
    writer_status_queue = mp.Queue()

    for _ in range(n_processes):
        args_queue.put(None)

    hog_proc = mp.Process(
        target=ordered_hog_writer_process,
        args=(
            hog_queue,
            writer_status_queue,
            output_writer.hog_writer,
            iogs_ordered,
            n_processes,
            spill_dir,
        )
    )

    non_hog_procs = [
        mp.Process(
            target=non_hog_writer_process,
            args=(
                non_hog_queues[writer_id],
                writer_status_queue,
                output_writer,
                n_processes,
                writer_id,
            )
        )
        for writer_id in range(n_writer_processes)
    ]

    running_processes = [
        mp.Process(
            target=worker_run_orthologs_method_pipeline,
            args=(
                tree_analyser,
                nspecies,
                args_queue,
                hog_queue,
                non_hog_queues,
                progress_queue,
                fewer_open_files,
                output_writer.need_ortholog_output,
                n_orthologs_cache,
                write_hog_tree,
                fix_files,
            )
        )
        for _ in range(n_processes)
    ]

    hog_proc.start()
    for proc in non_hog_procs:
        proc.start()
    for proc in running_processes:
        proc.start()

    # Start the progress bar only after all children exist: its refresh thread
    # holds the console lock while drawing, which a forked child could inherit.
    progressbar, task = printing.get_progressbar(total_tasks)
    progressbar.start()

    n_orthologs_sp_pair = stats.NOrthoSp(nspecies)
    errors = []   # error text for the WorkerError raised at the end

    def report_error(text):
        print(text)
        errors.append(text)

    completed_tasks = 0
    skipped_tasks = 0
    active_workers = n_processes
    finished_worker_pids = set()
    in_flight = {}
    exited_without_status = {}
    fatal = False
    hog_done = False
    non_hog_done_ids = set()
    hog_counts = None
    last_progress_time = time.monotonic()
    last_writer_activity_time = None
    child_processes = running_processes + [hog_proc] + non_hog_procs

    try:
        while (
            completed_tasks < total_tasks
            or active_workers > 0
            or not hog_done
            or len(non_hog_done_ids) < n_writer_processes
        ):
            try:
                msg = progress_queue.get(timeout=0.1)
            except queue.Empty:
                msg = "__EMPTY__"

            if msg == "__EMPTY__":
                pass

            elif isinstance(msg, tuple) and msg[0] == "worker_done":
                _, worker_pid, worker_counts = msg
                finished_worker_pids.add(worker_pid)
                active_workers = n_processes - len(finished_worker_pids)
                n_orthologs_sp_pair += worker_counts

            elif isinstance(msg, tuple) and msg[0] == "start":
                _, worker_pid, iog = msg
                in_flight[iog] = (worker_pid, time.monotonic(), "analysing")

            elif isinstance(msg, tuple) and msg[0] == "output":
                iog = msg[1]
                if iog in in_flight:
                    worker_pid, started, _ = in_flight[iog]
                    in_flight[iog] = (worker_pid, started, "waiting for output queue")

            elif isinstance(msg, tuple) and msg[0] == "error":
                report_error("ERROR: worker error:\n%s" % msg[2])
                fatal = True
                break

            elif isinstance(msg, tuple) and msg[0] == "skip":
                in_flight.pop(msg[1], None)
                skipped_tasks += 1
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.monotonic()

            elif isinstance(msg, tuple) and msg[0] == "result":
                iog = msg[1]
                in_flight.pop(iog, None)
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.monotonic()

            else:
                fatal = True
                raise TypeError(
                    "Unexpected progress message: %s %r" %
                    (type(msg), msg)
                )

            while True:
                try:
                    wmsg = writer_status_queue.get_nowait()
                except queue.Empty:
                    break

                if last_writer_activity_time is not None:
                    last_writer_activity_time = time.monotonic()

                if isinstance(wmsg, tuple) and wmsg[0] == "hog_done":
                    _, hog_counts = wmsg
                    hog_done = True

                elif isinstance(wmsg, tuple) and wmsg[0] == "non_hog_done":
                    _, writer_id = wmsg
                    non_hog_done_ids.add(writer_id)

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "hog_progress", "non_hog_progress",
                }:
                    pass

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "hog_error",
                    "hog_close_error",
                    "non_hog_error",
                    "non_hog_close_error",
                }:
                    report_error("ERROR: writer error:\n%s" % wmsg[-1])
                    fatal = True
                    break

                else:
                    fatal = True
                    raise TypeError(
                        "Unexpected writer status message: %s %r" %
                        (type(wmsg), wmsg)
                    )

            if fatal:
                break

            now = time.monotonic()
            process_statuses = [
                ("analysis worker", proc, proc.pid in finished_worker_pids)
                for proc in running_processes
            ] + [("HOG writer", hog_proc, hog_done)] + [
                ("non-HOG writer %d" % writer_id, proc,
                 writer_id in non_hog_done_ids)
                for writer_id, proc in enumerate(non_hog_procs)
            ]

            for role, proc, reported_done in process_statuses:
                exitcode = proc.exitcode
                if exitcode is None:
                    continue
                if exitcode != 0:
                    report_error("ERROR: %s (pid=%d) exited with code %d." %
                          (role, proc.pid, exitcode))
                    fatal = True
                    break
                if not reported_done:
                    # Queue messages can still be waiting in the parent when
                    # a child exits normally. Drain them before starting the
                    # grace period for a missing completion message.
                    if role == "analysis worker" and msg != "__EMPTY__":
                        continue
                    first_seen = exited_without_status.setdefault(proc.pid, now)
                    if now - first_seen > GRACE_PERIOD:
                        report_error("ERROR: %s (pid=%d) exited without reporting completion." %
                              (role, proc.pid))
                        fatal = True
                        break

            if fatal:
                break

            if active_workers == 0 and completed_tasks < total_tasks:
                report_error("ERROR: all analysis workers finished, but only %d/%d "
                      "orthogroups were reported." % (completed_tasks, total_tasks))
                fatal = True
                break

            if (
                (completed_tasks < total_tasks or active_workers > 0)
                and now - last_progress_time > STALL_TIMEOUT
            ):
                pending = ", ".join(
                    "OG%07d (pid=%d, %.0fs, %s)" % (iog, pid, now - started, stage)
                    for iog, (pid, started, stage) in sorted(
                        in_flight.items(), key=lambda item: item[1][1]
                    )[:5]
                )
                if any(
                    stage != "analysing"
                    for _, _, stage in in_flight.values()
                ):
                    pending += (
                        " Workers waiting for the output queue means the "
                        "output writers (disk) are the bottleneck."
                    )
                print("WARNING: Still waiting for analysis workers "
                      "(completed %d/%d, active_workers=%d).%s" % (
                          completed_tasks, total_tasks, active_workers,
                          " Pending: " + pending if pending else "",
                      ))
                last_progress_time = now

            if (
                completed_tasks >= total_tasks
                and active_workers == 0
                and (
                    not hog_done
                    or len(non_hog_done_ids) < n_writer_processes
                )
            ):
                if last_writer_activity_time is None:
                    # Start the drain timer here, not at pipeline startup.
                    last_writer_activity_time = now
                elif now - last_writer_activity_time > STALL_TIMEOUT:
                    print("WARNING: Still waiting for output writers to finish; "
                          "analysis is complete.")
                    last_writer_activity_time = now

    finally:
        if fatal or sys.exc_info()[0] is not None:
            # A failed consumer can leave producers blocked on full queues.
            # Stop all children before joining any of them in the error path.
            parallel_task_manager.terminate_processes(child_processes)
            args_queue.cancel_join_thread()
        else:
            # Everything reported completion: wait for the children to exit.
            # They are not stopped for taking long, only warned about.
            try:
                parallel_task_manager.wait_for_exit(
                    child_processes, "ortholog pipeline processes", GRACE_PERIOD
                )
            except parallel_task_manager.WorkerError as e:
                report_error("ERROR: %s" % e)
                fatal = True
                parallel_task_manager.terminate_processes(child_processes)
        for proc in child_processes:
            proc.join()
            if proc.exitcode != 0:
                if not fatal:
                    report_error("ERROR: child process (pid=%d) exited with code %d." %
                          (proc.pid, proc.exitcode))
                fatal = True

        progressbar.stop()

        for q in [hog_queue] + non_hog_queues + [progress_queue, writer_status_queue]:
            try:
                q.close()
                q.join_thread()
            except Exception:
                pass

    if hog_counts is not None:
        output_writer.hog_writer.i_hog.clear()
        output_writer.hog_writer.i_hog.update(hog_counts)

    skip_rate = skipped_tasks / max(1, total_tasks)
    if skip_rate > 0.02:
        print(
            "WARNING: skipped %d/%d tasks (%.1f%%)." %
            (skipped_tasks, total_tasks, 100.0 * skip_rate)
        )

    if not fatal and not hog_done:
        report_error("ERROR: ordered HOG writer did not finish cleanly.")
        fatal = True

    if not fatal and len(non_hog_done_ids) != n_writer_processes:
        report_error(
            "ERROR: only %d/%d non-HOG writers finished cleanly." %
            (len(non_hog_done_ids), n_writer_processes)
        )
        fatal = True

    if fatal:
        # Raised (not Fail()) so main() records the child's error in checkpoint.txt.
        raise parallel_task_manager.WorkerError(
            "\n".join(errors) or "Ortholog inference failed."
        )

    return n_orthologs_sp_pair

def set_file_descriptor_limit(fd_limit) -> None:
    """
    Try to raise the soft open-file limit.

    This is only a convenience. The lazy writer should not depend on this.
    """
    try:
        if isinstance(fd_limit, (tuple, list)):
            requested_soft = int(fd_limit[0])
        else:
            requested_soft = int(fd_limit)

        if requested_soft <= 0:
            print("Ignoring non-positive file descriptor limit.")
            return

        soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)

        print(
            "Current file descriptor limits: soft=%s, hard=%s" %
            (soft_limit, hard_limit)
        )

        if soft_limit >= requested_soft:
            print(
                "File descriptor soft limit is already sufficient: %s" %
                soft_limit
            )
            return

        if hard_limit != resource.RLIM_INFINITY and hard_limit < requested_soft:
            print(
                "Cannot raise soft file descriptor limit to %s; "
                "hard limit is only %s. Increase the hard limit outside "
                "OrthoFinder if needed." %
                (requested_soft, hard_limit)
            )
            return

        resource.setrlimit(
            resource.RLIMIT_NOFILE,
            (requested_soft, hard_limit)
        )

        new_soft, new_hard = resource.getrlimit(resource.RLIMIT_NOFILE)

        print(
            "New file descriptor limits: soft=%s, hard=%s" %
            (new_soft, new_hard)
        )

    except AttributeError:
        print("File descriptor limit functions not available on this platform.")
    except Exception as e:
        print("Could not adjust file descriptor limit: %s" % e)


