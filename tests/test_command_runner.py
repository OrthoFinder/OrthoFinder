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
