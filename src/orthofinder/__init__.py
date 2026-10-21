import multiprocessing as mp
import os
import sys

try:
    from ._version import __version__      # written by setuptools-scm when installed
except ImportError:
    __version__ = "unknown"                # run from the source without installing it


def _available_cores():
    """The cores this process may use (e.g. a job's allocation), else all of the machine's."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:   # not on every platform (e.g. macOS)
        return mp.cpu_count()


# The cores available, the most threads -t and -a can be
n_cores = _available_cores()
# The default number of search threads (-t): half the cores, so that the machine
# (and other programs on it) still has some
n_threads_default = max(1, n_cores // 2)

# MCL inflation parameter
g_mcl_inflation = 1.2

# Clade-specific genes
orphan_genes_version = 2

# Protocol vesion for pickling
# Updated to version 5, since the new OrthoFinder works on py38 and above
# version 5 add support for out-of-band data and speedup for in-band data.
pic_protocol = 5

# Get directory containing script/bundle
if getattr(sys, "frozen", False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )
