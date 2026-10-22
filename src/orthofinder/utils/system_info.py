"""
The computer a run is on, for run_state.json ("system"): what explains most
problems with the programs and limits a run meets (a bundled program that does not
run with this glibc, too few open files, a cluster job given fewer cores or less
memory than the node has, WSL or Rosetta 2). Best effort: a value that cannot be
read is None, and nothing here fails a run. No host or user name: results
directories are shared.
"""
import os
import platform
import subprocess


def _safe(function, default=None):
    try:
        value = function()
    except Exception:
        return default
    return default if value in ("", None) else value


def _read(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def _os_release():
    """{"name", "version"} of the Linux distribution (/etc/os-release)."""
    fields = {}
    for line in _read("/etc/os-release").splitlines():
        key, sep, value = line.partition("=")
        if sep:
            fields[key.strip()] = value.strip().strip('"')
    return {"name": fields.get("PRETTY_NAME") or fields.get("NAME"),
            "version": fields.get("VERSION_ID") or fields.get("VERSION")}


def _sysctl(name):
    """A sysctl value on macOS (no shell; a short timeout)."""
    return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=2).stdout.strip()


def _cpu_model():
    system = platform.system()
    if system == "Linux":
        for line in _read("/proc/cpuinfo").splitlines():
            key, _, value = line.partition(":")
            # (x86: "model name"; ARM: "Model" or "CPU part" only)
            if key.strip() in ("model name", "Model", "Hardware"):
                return value.strip()
    if system == "Darwin":
        return _sysctl("machdep.cpu.brand_string")
    return platform.processor()


def _memory_mb():
    """{"total", "available"} of the computer in MB."""
    if platform.system() == "Linux":
        values = {}
        for line in _read("/proc/meminfo").splitlines():
            key, _, value = line.partition(":")
            if key in ("MemTotal", "MemAvailable"):
                values[key] = int(value.split()[0]) // 1024
        return {"total": values.get("MemTotal"), "available": values.get("MemAvailable")}
    total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // (1024 * 1024)
    return {"total": total, "available": None}


def _cgroup_memory_limit_mb():
    """The memory limit of this process's cgroup (a container or a cluster job) in MB, None if none."""
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        if os.path.exists(path):
            value = _read(path).strip()
            if value == "max" or not value.isdigit():
                return None
            limit = int(value)
            # (cgroup v1 writes a huge number for "no limit")
            return None if limit >= 2 ** 60 else limit // (1024 * 1024)
    return None


def _open_files_limit():
    """{"soft", "hard"} open-file limits now (after OrthoFinder raised the soft one); None: unlimited."""
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    unlimited = resource.RLIM_INFINITY
    return {"soft": None if soft == unlimited else soft, "hard": None if hard == unlimited else hard}


def _wsl():
    return "microsoft" in _read("/proc/version").lower()


def _rosetta():
    return _sysctl("sysctl.proc_translated") == "1"


def _container():
    """docker, podman, singularity / apptainer, or None."""
    if os.environ.get("APPTAINER_CONTAINER") or os.environ.get("APPTAINER_NAME"):
        return "apptainer"
    if os.environ.get("SINGULARITY_CONTAINER") or os.environ.get("SINGULARITY_NAME"):
        return "singularity"
    if os.path.exists("/.dockerenv"):
        return "docker"
    if os.path.exists("/run/.containerenv"):
        return "podman"
    return None


# The variables of a cluster job that explain its cores and memory
JOB_VARIABLES = (
    "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU",
    "PBS_JOBID", "PBS_NUM_PPN", "NCPUS", "LSB_JOBID", "LSB_DJOB_NUMPROC", "JOB_ID", "NSLOTS",
)


def system_info():
    """
    The computer the run is on: {"os": {"system", "release", "version", "machine",
    "distribution", "libc", "wsl", "rosetta2"}, "cpu": {"model", "cores", "cores_available"},
    "memory_mb": {"total", "available", "limit"}, "open_files": {"soft", "hard"},
    "container", "job"}.
    """
    from .. import n_cores
    system = platform.system()
    libc = _safe(platform.libc_ver, ("", ""))
    return {
        "os": {
            "system": _safe(platform.system),
            "release": _safe(platform.release),         # (the kernel's, on Linux)
            "version": _safe(platform.version),
            "machine": _safe(platform.machine),
            "distribution": _safe(_os_release) if system == "Linux" else (
                _safe(lambda: {"name": "macOS", "version": platform.mac_ver()[0]}) if system == "Darwin" else
                _safe(lambda: {"name": "Windows", "version": platform.win32_ver()[1]}) if system == "Windows" else None),
            "libc": {"name": libc[0], "version": libc[1]} if libc[0] else None,
            "wsl": _safe(_wsl, False) if system == "Linux" else False,
            "rosetta2": _safe(_rosetta, False) if system == "Darwin" else False,
        },
        "cpu": {
            "model": _safe(_cpu_model),
            "cores": _safe(os.cpu_count),
            "cores_available": n_cores,                 # (this process's CPU affinity: e.g. a job's)
        },
        "memory_mb": dict(_safe(_memory_mb, {"total": None, "available": None}),
                          limit=_safe(_cgroup_memory_limit_mb) if system == "Linux" else None),
        "open_files": _safe(_open_files_limit),
        "container": _safe(_container),
        "job": {name: os.environ[name] for name in JOB_VARIABLES if os.environ.get(name)} or None,
    }
