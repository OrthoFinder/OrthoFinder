"""
Which copy of each external program OrthoFinder runs.

The dependency check chooses, for each program its test commands run, the
first copy that works, looking in the directories of
parallel_task_manager.program_dirs in order (the PATH OrthoFinder runs
commands with). Every command OrthoFinder then runs starts the chosen copy:
with_paths puts its full path in place of the program's name, and with_env
adds the environment variables that copy needs (the bundled MAFFT's
MAFFT_BINARIES).
"""
import os
import re
import shlex
import shutil
import subprocess

from .. import __location__

# The separators of the commands in a command line (kept by split)
_SEPARATORS = re.compile(r"(\n|;|&&|\|\||\|)")
# The start of a command: spaces, VAR=value assignments, then the program
_COMMAND_START = re.compile(r"(\s*)((?:[A-Za-z_][A-Za-z0-9_]*=\S*\s+)*)(\S+)")

# Shell builtins and keywords: never replaced by a program found on PATH
SHELL_WORDS = {
    ".", ":", "[", "alias", "bg", "break", "case", "cd", "command", "continue", "do", "done",
    "echo", "elif", "else", "esac", "eval", "exec", "exit", "export", "false", "fg", "fi",
    "for", "function", "if", "jobs", "kill", "printf", "pwd", "read", "return", "set",
    "shift", "source", "test", "then", "time", "trap", "true", "type", "ulimit", "umask",
    "unset", "until", "wait", "while", "{", "}", "(", ")",
}

_chosen = {}        # program name -> full path of the copy to run
_env = {}           # program name -> {variable: value} that copy needs


def command_programs(command):
    """The programs a command line runs: the first word of each command in it (not VAR=value)."""
    programs = []
    for part in _SEPARATORS.split(command)[::2]:
        m = _COMMAND_START.match(part)
        if m and m.group(3) not in programs:
            programs.append(m.group(3))
    return programs


def can_choose(program):
    """Whether OrthoFinder chooses the copy of a program: a name, not a path or a shell word."""
    return program not in SHELL_WORDS and os.sep not in program and "/" not in program


def bundled_env(path):
    """The environment variables a copy needs: the bundled MAFFT needs its libexec directory."""
    bin_dir = os.path.join(__location__, "bin", "mafft", "bin")
    if os.path.basename(path) == "mafft" and os.path.dirname(os.path.abspath(path)) == bin_dir:
        return {"MAFFT_BINARIES": os.path.join(__location__, "bin", "mafft", "libexec")}
    return {}


def candidates(program, path, recorded=None):
    """
    The copies of a program, in the order they are tried: the one recorded
    for the run being continued (if it is still there), then one per
    directory of path (an os.pathsep-separated list), in its order.
    """
    found = []
    seen = set()
    first = [recorded] if recorded else []
    for location in first + [shutil.which(program, path=d) for d in path.split(os.pathsep) if d]:
        if location is None or not (os.path.isfile(location) and os.access(location, os.X_OK)):
            continue
        real = os.path.realpath(location)
        if real not in seen:
            seen.add(real)
            found.append(os.path.abspath(location))
    return found


def choose(program, location):
    """Run this copy of the program from now on."""
    _chosen[program] = location
    _env[program] = bundled_env(location)


def forget(program):
    """Run the program by its name again (the first copy on PATH)."""
    _chosen.pop(program, None)
    _env.pop(program, None)


def chosen():
    """{program: full path} of the copies chosen."""
    return dict(_chosen)


def _quote(path):
    return subprocess.list2cmdline([path]) if os.name == "nt" else shlex.quote(path)


def with_paths(command):
    """The command line with the full path of the chosen copy of each program it runs."""
    if not _chosen:
        return command
    parts = _SEPARATORS.split(command)
    for i in range(0, len(parts), 2):
        m = _COMMAND_START.match(parts[i])
        if m and m.group(3) in _chosen:
            parts[i] = parts[i][:m.start(3)] + _quote(_chosen[m.group(3)]) + parts[i][m.end(3):]
    return "".join(parts)


def with_env(env):
    """env with the variables the chosen copies need (env itself if there are none)."""
    extra = {k: v for values in _env.values() for k, v in values.items()}
    if not extra:
        return env
    env = dict(env)
    env.update(extra)
    return env
