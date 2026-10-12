import os
import subprocess

import pytest

from orthofinder.utils import command_runner


def test_commands_finishing_after_a_failure_are_recorded(tmp_path):
    # "late" ignores the SIGTERM sent when "exit 1" fails, so it succeeds only
    # after the batch has begun stopping; "killed" is stopped and must not be recorded.
    late = "trap '' TERM; sleep 1; echo done > %s" % (tmp_path / "late.txt")
    killed = "sleep 30; echo done > %s" % (tmp_path / "killed.txt")
    failing = "sleep 0.2; exit 1"
    recorded = []
    with pytest.raises(subprocess.CalledProcessError):
        command_runner.run_parallel_commands(
            3, [late, killed, failing], q_list_of_list=False, on_success=recorded.append)
    assert recorded == [late]
    assert (tmp_path / "late.txt").exists()
    assert not (tmp_path / "killed.txt").exists()


def test_follow_up_commands_run_before_queued_ones(tmp_path):
    # One at a time, so the order they run in is the order they were started in.
    log = tmp_path / "log.txt"
    units = [[("echo %s >> %s" % (name, log), None)] for name in ("a", "b", "c")]
    follow_up_unit = [("echo after_a >> %s" % log, None)]

    def follow_up(unit):
        return [follow_up_unit] if unit is units[0] else []

    command_runner.run_parallel_commands_and_move_results_file(1, units, True, follow_up=follow_up)
    assert log.read_text().split() == ["a", "after_a", "b", "c"]


def write_pid(fn):
    with open(fn, "w") as f:
        f.write(str(os.getpid()))


def fail_in_function(fn):
    raise ValueError("trimming failed")


def test_functions_run_in_a_separate_process(tmp_path):
    pid_fn = str(tmp_path / "pid.txt")
    done_fn = tmp_path / "done.txt"
    units = [[("echo done > %s" % done_fn, None), (write_pid, (pid_fn,))]]
    command_runner.run_parallel_commands_and_move_results_file(2, units, True)
    assert done_fn.exists()
    assert int(open(pid_fn).read()) != os.getpid()


def test_a_failing_function_stops_the_batch(tmp_path):
    later = tmp_path / "later.txt"
    units = [[(fail_in_function, (str(tmp_path / "x"),))],
             [("sleep 1; echo done > %s" % later, None)]]
    with pytest.raises(ValueError):
        command_runner.run_parallel_commands_and_move_results_file(1, units, True)
    assert not later.exists()
