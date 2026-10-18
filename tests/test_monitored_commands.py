"""run_monitored_command and kill_running_commands: no command is left running after an abort."""
import os
import subprocess
import threading
import time

import pytest

from orthofinder.utils import parallel_task_manager as ptm

ENV = dict(os.environ)


@pytest.fixture(autouse=True)
def reset():
    ptm.reset_command_abort()
    yield
    ptm.reset_command_abort()


def test_runs_command():
    rc, out, _ = ptm.run_monitored_command("echo ok", ENV, stdout=subprocess.PIPE)
    assert rc == 0 and out.strip() == b"ok"


def test_not_started_after_abort():
    ptm.kill_running_commands()
    with pytest.raises(ptm.WorkerError):
        ptm.run_monitored_command("true", ENV)


def test_abort_stops_running_command():
    result = {}

    def run():
        result["rc"] = ptm.run_monitored_command("sleep 30", ENV)[0]
    t = threading.Thread(target=run)
    start = time.monotonic()
    t.start()
    while not ptm._live_commands and time.monotonic() - start < 5:
        time.sleep(0.01)
    ptm.kill_running_commands()
    t.join(15)
    assert not t.is_alive() and result["rc"] != 0 and time.monotonic() - start < 15
    assert not ptm._live_commands


def test_abort_while_starting(monkeypatch):
    # The abort comes after the command has started but before it is registered:
    # it is stopped, not left running
    started = []
    real_popen = subprocess.Popen

    def popen(*args, **kwargs):
        p = real_popen(*args, **kwargs)
        started.append(p)
        ptm.kill_running_commands()     # (does not see p: not registered yet)
        return p
    monkeypatch.setattr(ptm.subprocess, "Popen", popen)
    with pytest.raises(ptm.WorkerError):
        ptm.run_monitored_command("sleep 30", ENV)
    assert started[0].poll() is not None and not ptm._live_commands
