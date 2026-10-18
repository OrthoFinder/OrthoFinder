"""Where OrthoFinder looks for the programs it runs, which copy it runs, and which it records."""
import os
import re
import subprocess

import pytest

from orthofinder import __location__
from orthofinder.run import check_dependencies
from orthofinder.utils import parallel_task_manager, programs


@pytest.fixture(autouse=True)
def no_choices():
    programs._chosen.clear()
    programs._env.clear()
    check_dependencies._confirmed.clear()
    yield
    programs._chosen.clear()
    programs._env.clear()
    check_dependencies._confirmed.clear()


def _env(tmp_path, monkeypatch, active=None, base=None, bundled=True):
    home = tmp_path / "home"
    (home / "local" / "bin").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    (tmp_path / "usr_bin").mkdir()
    env = {"PATH": str(tmp_path / "usr_bin")}
    if base is not None:
        (tmp_path / base / "bin").mkdir(parents=True)
        env["CONDA_EXE"] = str(tmp_path / base / "bin" / "conda")
    if active is not None:
        (tmp_path / active / "bin").mkdir(parents=True, exist_ok=True)
        env["CONDA_PREFIX"] = str(tmp_path / active)
    location = tmp_path / "orthofinder"
    (location / "bin" / "mafft" / "bin").mkdir(parents=True)
    monkeypatch.setattr(parallel_task_manager, "bundled_binaries_run_here", lambda: bundled)
    return env, str(location), home


def test_order_active_base_local_path_bundled(tmp_path, monkeypatch):
    env, location, home = _env(tmp_path, monkeypatch, active="envs/of3", base="miniconda3")
    assert parallel_task_manager.program_dirs(env, location) == [
        str(tmp_path / "envs/of3/bin"), str(tmp_path / "miniconda3/bin"), str(home / "local/bin"),
        str(tmp_path / "usr_bin"), os.path.join(location, "bin"), os.path.join(location, "bin", "mafft", "bin")]


def test_base_active_listed_once(tmp_path, monkeypatch):
    env, location, home = _env(tmp_path, monkeypatch, active="miniconda3", base="miniconda3")
    env["PATH"] = os.pathsep.join([str(tmp_path / "miniconda3/bin"), env["PATH"]])
    assert parallel_task_manager.program_dirs(env, location)[:3] == [
        str(tmp_path / "miniconda3/bin"), str(home / "local/bin"), str(tmp_path / "usr_bin")]


def test_bundled_only_where_they_run(tmp_path, monkeypatch):
    env, location, home = _env(tmp_path, monkeypatch, bundled=False)
    assert parallel_task_manager.program_dirs(env, location) == [str(home / "local/bin"), str(tmp_path / "usr_bin")]


def test_command_programs():
    assert programs.command_programs(
        "mmseqs createdb A B.fa ; mmseqs createindex B.fa B.tmp\nOMP_NUM_THREADS=1 FastTree x > y | gzip"
    ) == ["mmseqs", "FastTree", "gzip"]


def test_with_paths():
    programs.choose("diamond", "/opt/my tools/diamond")
    programs.choose("gzip", "/bin/gzip")
    programs.choose("cd", "/usr/bin/cd")      # (never chosen by the check: a shell word)
    assert programs.with_paths("diamond blastp -q diamond.fa | gzip > out.gz && cd x; A=1 B=2 gzip y") == (
        "'/opt/my tools/diamond' blastp -q diamond.fa | /bin/gzip > out.gz && /usr/bin/cd x; A=1 B=2 /bin/gzip y")
    assert programs.with_paths("mafft in > out") == "mafft in > out"
    assert not programs.can_choose("cd") and not programs.can_choose("/usr/bin/x") and programs.can_choose("mafft")


def test_bundled_mafft_needs_its_libexec():
    mafft = os.path.join(__location__, "bin", "mafft", "bin", "mafft")
    programs.choose("mafft", mafft)
    assert programs.with_env({"PATH": "x"}) == {
        "PATH": "x", "MAFFT_BINARIES": os.path.join(__location__, "bin", "mafft", "libexec")}
    programs.choose("mafft", "/usr/bin/mafft")
    env = {"PATH": "x"}
    assert programs.with_env(env) is env


def _out(capsys):
    return re.sub(r"\x1b\[[0-9;]*m", "", capsys.readouterr().out)


def _program(directory, name, works):
    directory.mkdir(parents=True, exist_ok=True)
    fn = directory / name
    fn.write_text("#!/bin/sh\n" + ("echo ok\n" if works else "echo 'error while loading shared libraries' >&2; exit 1\n"))
    fn.chmod(0o755)
    return str(fn)


def _test(command):
    def test():
        rc, _, _ = parallel_task_manager.run_monitored_command(
            command, parallel_task_manager.my_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return rc == 0, command, None
    return test


def _path(monkeypatch, *dirs):
    monkeypatch.setitem(parallel_task_manager.my_env, "PATH", os.pathsep.join(str(d) for d in dirs))


def test_candidates_recorded_first_each_copy_once(tmp_path):
    broken = _program(tmp_path / "a", "prog", False)
    works = _program(tmp_path / "b", "prog", True)
    (tmp_path / "c").mkdir()
    os.symlink(works, tmp_path / "c" / "prog")
    path = os.pathsep.join(str(tmp_path / d) for d in ("a", "b", "c", "missing"))
    assert programs.candidates("prog", path) == [broken, works]
    assert programs.candidates("prog", path, recorded=works) == [works, broken]


def test_broken_copy_is_passed_over(tmp_path, monkeypatch, capsys):
    broken = _program(tmp_path / "a", "prog", False)
    works = _program(tmp_path / "b", "prog", True)
    _path(monkeypatch, tmp_path / "a", tmp_path / "b")
    success, cmd, _ = check_dependencies.run_test(_test("prog -h"))
    assert success and programs.chosen() == {"prog": works}
    assert "the copy at %s does not work" % broken in _out(capsys)
    assert programs.with_paths("prog x") == works + " x"
    assert check_dependencies.record_programs([cmd], None) == {"prog": works}


def test_whole_installation_that_does_not_run(tmp_path, monkeypatch):
    # both programs of a command broken in the first directory: neither works on its own
    for name in ("makedb", "search"):
        _program(tmp_path / "a", name, False)
        _program(tmp_path / "b", name, True)
    _path(monkeypatch, tmp_path / "a", tmp_path / "b")
    assert check_dependencies.run_test(_test("makedb x && search y"))[0]
    assert programs.chosen() == {n: str(tmp_path / "b" / n) for n in ("makedb", "search")}


def test_none_works_keeps_first_and_says_why(tmp_path, monkeypatch, capsys):
    first = _program(tmp_path / "a", "prog", False)
    other = _program(tmp_path / "b", "prog", False)
    _path(monkeypatch, tmp_path / "a", tmp_path / "b")
    assert not check_dependencies.run_test(_test("prog -h"))[0]
    assert programs.chosen() == {"prog": first}
    out = _out(capsys)
    assert "None of the other copies works" in out and other in out


def test_confirmed_copy_not_changed_by_later_test(tmp_path, monkeypatch):
    first = _program(tmp_path / "a", "prog", True)
    _program(tmp_path / "b", "prog", True)
    _path(monkeypatch, tmp_path / "a", tmp_path / "b")
    assert check_dependencies.run_test(_test("prog -h"))[0]
    assert not check_dependencies.run_test(_test("prog -h && false"))[0]
    assert programs.chosen() == {"prog": first}


def test_restart_prefers_recorded_copy(tmp_path, monkeypatch, capsys):
    first = _program(tmp_path / "a", "prog", True)
    recorded = _program(tmp_path / "b", "prog", True)
    _path(monkeypatch, tmp_path / "a", tmp_path / "b")
    check_dependencies.choose_recorded({"prog": recorded, "gone": str(tmp_path / "gone")})
    assert check_dependencies.run_test(_test("prog -h"))[0]
    assert programs.chosen() == {"prog": recorded}
    check_dependencies.warn_changed_programs({"prog": first}, programs.chosen())
    assert "this restart uses" in _out(capsys)


def test_record_python(tmp_path):
    python = check_dependencies.record_python(None)     # (no results directory: only returned)
    assert os.path.isfile(python["executable"]) and python["version"].count(".") == 2
    assert os.path.isdir(python["orthofinder"]["path"])
    numpy = python["packages"]["numpy"]
    assert numpy["version"] and os.path.isfile(os.path.join(numpy["path"], "__init__.py"))
    assert set(python["packages"]) == set(check_dependencies.PYTHON_PACKAGES)


def test_record_system(tmp_path, monkeypatch):
    from orthofinder.utils import run_state, system_info
    system = check_dependencies.record_system(None)      # (no results directory: only returned)
    assert set(system) == {"os", "cpu", "memory_mb", "open_files", "container", "job"}
    assert system["os"]["system"] and system["cpu"]["cores_available"] >= 1
    # a value that cannot be read is None: nothing fails
    monkeypatch.setattr(system_info, "_cpu_model", lambda: 1 / 0)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    system = system_info.system_info()
    assert system["cpu"]["model"] is None and system["job"] == {"SLURM_JOB_ID": "123"}
    # written to run_state.json, as the other sections
    (tmp_path / run_state.WORKING_DIR_NAME).mkdir()
    check_dependencies.record_system(str(tmp_path))
    assert run_state.read(str(tmp_path))["system"]["os"]["system"] == system["os"]["system"]


def test_command_keys():
    from orthofinder.utils import run_state
    assert run_state.command_key("diamond blastp -d db -q q.fa") == "diamond blastp"
    assert run_state.command_key("/bin/diamond makedb --in a.fa -d db") == "diamond makedb"
    assert run_state.command_key("OMP_NUM_THREADS=1 FastTree -fastest a.fa > t") == "FastTree"
    assert run_state.command_key("mmseqs createdb a b ; mmseqs createindex b t") == "mmseqs createdb"
    assert run_state.command_key("FastTree /wd/Alignments_ids/SpeciesTreeAlignment.fa > t") == "FastTree (species tree)"
    assert run_state.command_key("mafft --localpair a.fa > b") == "mafft"
    assert run_state.command_key("  ") is None


def test_note_command_keeps_one_example_of_each(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from orthofinder.utils import file_handler, run_state
    (tmp_path / "WorkingDirectory").mkdir()
    monkeypatch.setattr(file_handler, "handler", SimpleNamespace(rd1=str(tmp_path) + "/"))
    monkeypatch.setattr(run_state, "_noted", set())
    monkeypatch.setattr(run_state, "_tested", set())
    commands = lambda: run_state.read(str(tmp_path))["commands"]
    run_state.record_run(str(tmp_path), "core", "orthofinder -f in")
    with run_state.stage(run_state.DEPENDENCY_CHECK):
        run_state.note_command("diamond blastp -d test")
    assert commands()["diamond blastp"]["stage"] == run_state.DEPENDENCY_CHECK
    # the run's first command replaces the test's; later ones, of any stage, do not
    with run_state.stage("Sequence search"):
        run_state.note_command("diamond blastp -d first")
        run_state.note_command("diamond blastp -d second")
    with run_state.stage("Search orthogroup profiles"):
        run_state.note_command("diamond blastp -d third")
    with run_state.stage(run_state.DEPENDENCY_CHECK):
        run_state.note_command("diamond blastp -d test2")
    assert list(commands()) == ["diamond blastp"]
    assert commands()["diamond blastp"]["command"].endswith("-d first")
    assert commands()["diamond blastp"]["stage"] == "Sequence search"
    # another process of the invocation (its own sets): the first stays
    monkeypatch.setattr(run_state, "_noted", set())
    with run_state.stage("Sequence search"):
        run_state.note_command("diamond blastp -d other_process")
    assert commands()["diamond blastp"]["command"].endswith("-d first")
    # a continuation replaces it
    monkeypatch.setattr(run_state, "_noted", set())
    run_state.record_run(str(tmp_path), "core", "orthofinder -b x")
    with run_state.stage("Sequence search"):
        run_state.note_command("diamond blastp -d continued")
    state = run_state.read(str(tmp_path))
    assert state["commands"]["diamond blastp"]["command"].endswith("-d continued")
    assert state["commands"]["diamond blastp"]["invocation"] == 2
    assert [i["number"] for i in state["run"]["invocations"]] == [1, 2]
    assert state["run"]["invocations"][1]["command"] == "orthofinder -b x"
