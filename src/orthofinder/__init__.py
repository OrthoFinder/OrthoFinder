import multiprocessing as mp
import os
import sys

try:
    from importlib.metadata import PackageNotFoundError, version
    __version__ = version(__name__)
except PackageNotFoundError:
    from ._version import __version__ as __version__

    
# Find the total number of threads on the host machine
nThreadsDefault = mp.cpu_count()

# MCL inflation parameter
g_mclInflation = 1.2

# Clade-specific genes
orphan_genes_version = 2

# Protocol vesion for pickling
# Updated to version 5, since the new OrthoFinder works on py38 and above
# version 5 add support for out-of-band data and speedup for in-band data.
picProtocol = 5

# Get directory containing script/bundle
if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
