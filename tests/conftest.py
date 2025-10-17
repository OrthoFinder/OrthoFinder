import os
import sys
import pytest
import helper
from typing import Optional, Union, List, Any
import multiprocessing as mp
from orthofinder.utils.util import CreateNewWorkingDirectory

## ------- Default ExampleData ---------
# HERE = Path(__file__).resolve()
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) # project root

EXAMPLEDATA = os.path.join(ROOT, "ExampleData")
EXAMPLE_RESULTS = os.path.join(ROOT, "ExampleData", "OrthoFinder")

ORTHOFINDER = os.path.join(ROOT, "src")
sys.path.insert(0, ORTHOFINDER)

ASSIGN = os.path.join(ROOT, "ExampleData", "AdditionalSpecies")
# EXPECTED_RESULTS = os.path.join(ROOT, "tests", "expected_results")

if EXAMPLE_RESULTS[-1] != os.sep:
    EXAMPLE_RESULTS += os.sep 

if ORTHOFINDER[-1] != os.sep:
    ORTHOFINDER += os.sep

# exampledata_baseDirectoryName = os.path.join(ROOT, "ExampleData", "OrthoFinder", "Results_")

# DEFAULT_PROJECTS_RESULTS = CreateNewWorkingDirectory(
#     exampledata_baseDirectoryName,
#     qDate=True,
#     search_program=None,
#     msa_program=None,
#     tree_program=None,
#     scorematrix=None,
#     gapopen=None,
#     gapextend=None,
#     extended_filename=False,
#     makedir=False
# )

DEFAULT_PROJECTS = EXAMPLEDATA
# DEFAULT_PROJECTS_RESULTS = _latest_output_dir()  # used if --projects not supplied
# working_dir = os.path.join(DEFAULT_PROJECTS_RESULTS, "WorkingDirectory")

## ------------ Default testdata --------------
TESTDATA = "proteome"
TESTDATA_DNA = "dna"

TESTDATA_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "data")) # tests root

TESTDATA_CORE = os.path.join(TESTDATA_ROOT, TESTDATA, "Core")

TESTDATA_ASSIGN = os.path.join(TESTDATA_ROOT, TESTDATA, "Assign")
TESTDATA_SPECIES_TREE_CORE = os.path.join(TESTDATA_ROOT, TESTDATA, "species_tree.txt")
TESTDATA_SPECIES_TREE_ASSIGN = os.path.join(TESTDATA_ROOT, TESTDATA, "species_tree_assign.txt")
TESTDATA_RESULTS = os.path.join(TESTDATA_CORE, "OrthoFinder")

TESTDATA_DNA_CORE = os.path.join(TESTDATA_ROOT, TESTDATA_DNA, "Genomes_Core")
# TESTDATA_DNA_ASSIGN = os.path.join(TESTDATA_ROOT, TESTDATA_DNA, "Genomes_Assign")
TESTDATA_DNA_RESULTS = os.path.join(TESTDATA_DNA_CORE, "OrthoFinder")

if TESTDATA_RESULTS[-1] != os.sep:
    TESTDATA_RESULTS += os.sep 

if TESTDATA_DNA_RESULTS[-1] != os.sep:
    TESTDATA_DNA_RESULTS += os.sep 


EXPECTED_RESULTS = os.path.join(TESTDATA_ROOT, "proteome_expected_results")


# exampledata_results_dir_list = sorted([
#     (entry.stat().st_mtime, entry.path) 
#     for entry in os.scandir(EXAMPLE_RESULTS)
# ])

# latest_output_dir = exampledata_results_dir_list[-1][1]
# working_dir = os.path.join(latest_output_dir, "WorkingDirectory")


# def _latest_output_dir() -> str:
#     try:
#         entries = sorted((e.stat().st_mtime, e.path) for e in os.scandir(EXAMPLE_RESULTS))
#         return entries[-1][1] if entries else ""
#     except FileNotFoundError:
#         return ""

DEFAULT_PROJECTS = {
    "exampledata":{
        "proteome":
        {
            "core": EXAMPLEDATA,
            "assign": ASSIGN,
        },
    },
    "testdata": {
        "proteome": {
            "core": TESTDATA_CORE,
            "assign": TESTDATA_ASSIGN,
            "species_tree_core": TESTDATA_SPECIES_TREE_CORE,
            "species_tree_assign": TESTDATA_SPECIES_TREE_ASSIGN
        },
        "dna":{
            "core": TESTDATA_DNA_CORE,
            # "assign": TESTDATA_DNA_ASSIGN,
        }
    }
}

def pytest_addoption(parser):
    parser.addoption(
        "--projects",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["proteome"]["core"],  
        help="Comma-separated list of project paths (defaults to latest ExampleData/OrthoFinder run).",
    )

    # parser.addoption(
    #     "--projects-results",
    #     action="store",
    #     default=DEFAULT_PROJECTS_RESULTS, 
    #     help="Comma-separated list of project result paths (defaults to latest ExampleData/OrthoFinder run).",
    # )

    parser.addoption(
        "--testdata",
        action="store_false",
        default=True,  
        help="Decide whether or not to use the test dataset.",
    )

    parser.addoption(
        "--dna",
        action="store_true",
        default=False,  
        help="Whether or not the input is DNA.",
    )

    parser.addoption(
        "--dna-projects",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["dna"]["core"],  
        help="Comma-separated list of DNA project paths (defaults to latest ExampleData/OrthoFinder run).",
    )

    parser.addoption(
        "--assign",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["proteome"]["assign"], 
        help="Comma-separated list of additional species paths.",
    )

    parser.addoption(
        "--expected-results",
        action="store",
        default=EXPECTED_RESULTS, 
        help="Comma-separated list of expected results paths.",
    )

    parser.addoption(
        "--cluster",
        action="store",
        default="mcl",  
        help="Clustering program",
    )

    parser.addoption(
        "--msa",
        action="store",
        default="famsa",  
        help="MSA program",
    )

    parser.addoption(
        "--gene-tree-method",
        action="store",
        default="fasttree",  
        help="Gene tree program",
    )

    parser.addoption("--dendroblast", action="store_true", default=False, help="Non-MSA program")

    parser.addoption(
        "--species-tree",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["proteome"]["species_tree_core"],  
        help="Species tree file",
    )

    parser.addoption(
        "--species-tree-assign",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["proteome"]["species_tree_assign"], 
        help="Species tree file for assign",
    )

    parser.addoption(
        "--species-tree-method",
        action="store",
        default="astral-pro",  
        help="Species tree program",
    )

    parser.addoption(
        "--recon-method",
        action="store",
        default="of_recon",  
        help="Rreconciliation method",
    )

    parser.addoption("--t", action="store", type=int, default=mp.cpu_count(),
                     help="Number of parallel sequence search threads")
    parser.addoption("--a", action="store", type=int,
                     default=min(16, max(1, int(mp.cpu_count() / 8))),
                     help="Number of parallel analysis threads")


def _to_list(value: Optional[Union[str, int, bool]]) -> List[Any]:
    """
    Normalize a CLI option value to a list:
      - str: split by commas (empty -> [])
      - int/bool: wrap as single-element list
      - None: []
      - list: return as-is (defensive)
    """
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        value = value.strip()
        return [x for x in value.split(",") if x] if value else []
    if isinstance(value, (int, bool)):
        return [value]
    # fallback
    return [value]

def _split_or_default(value: Optional[Union[str, int, bool]],
                      default_value: Optional[Union[str, int, bool]]) -> List[Any]:
    """
    Use value if present (after normalization); otherwise use default.
    Both paths return lists.
    """
    lst = _to_list(value)
    if lst:
        return lst
    return _to_list(default_value)

def _broadcast_to_len(lst: list, target_len: int, name: str):
    """
    Broadcast list `lst` to length `target_len`:
      - empty -> [None] * target_len
      - len==1 -> repeat
      - len==target_len -> as is
      - else -> error
    """
    if target_len <= 0:
        return lst
    if not lst:
        return [None] * target_len
    if len(lst) == 1:
        return lst * target_len
    if len(lst) == target_len:
        return lst
    raise pytest.usage_error(
        f"Length mismatch for '{name}': got {len(lst)} values, need {target_len}. "
        f"Provide 1 value to broadcast or {target_len} values to align."
    )

def pytest_generate_tests(metafunc):
    requested = [name for name in (
        "projects",
        "dna_projects",
        "testdata",
        "dna",
        # "projects_results",
        "assign",
        "expected_results",
        "cluster",
        "msa",
        "gene_tree_method",
        "dendroblast",
        "species_tree",
        "species_tree_assign",
        "species_tree_method",
        "recon_method",
        "t",
        "a",
    ) if name in metafunc.fixturenames]
    if not requested:
        return

    values_by_name = {
        "projects": _split_or_default(
            metafunc.config.getoption("--projects"), 
            default_value=DEFAULT_PROJECTS["testdata"]["proteome"]["core"]
        ),

        "dna_projects": _split_or_default(
            metafunc.config.getoption("--dna-projects"), 
            default_value=DEFAULT_PROJECTS["testdata"]["dna"]["core"]
        ),

        "testdata": _split_or_default(
            metafunc.config.getoption("--testdata"), default_value=True
        ),

        "dna": _split_or_default(
            metafunc.config.getoption("--dna"), default_value=False
        ),
        # "projects_results": _split_or_default(
        #     metafunc.config.getoption("--projects-results"), default_value=DEFAULT_PROJECTS_RESULTS
        # ),
        "assign": _split_or_default(
            metafunc.config.getoption("--assign"), default_value=DEFAULT_PROJECTS["testdata"]["proteome"]["assign"]
        ),
        "expected_results": _to_list(
            metafunc.config.getoption("--expected-results")
        ),
        "cluster": _split_or_default(
            metafunc.config.getoption("--cluster"), default_value="mcl"
        ),
        "msa": _split_or_default(
            metafunc.config.getoption("--msa"), default_value="famsa"
        ),
        "gene_tree_method": _split_or_default(
            metafunc.config.getoption("--gene-tree-method"), default_value="fasttree"
        ),
        "dendroblast": _split_or_default(
            metafunc.config.getoption("--dendroblast"), default_value=False  # ← keep as bool
        ),
        "species_tree_method": _split_or_default(
            metafunc.config.getoption("--species-tree-method"), default_value="astral-pro"
        ),

        "species_tree": _split_or_default(
            metafunc.config.getoption("--species-tree"), 
            default_value=DEFAULT_PROJECTS["testdata"]["proteome"]["species_tree_core"]
        ),
        "species_tree_assign": _split_or_default(
            metafunc.config.getoption("--species-tree-assign"), 
            default_value=DEFAULT_PROJECTS["testdata"]["proteome"]["species_tree_assign"]
        ),
        "recon_method": _split_or_default(
            metafunc.config.getoption("--recon-method"), default_value="of_recon"
        ),
        "sequence_search_threads": _split_or_default(
            metafunc.config.getoption("--t"), default_value=mp.cpu_count()
        ),
        "analysis_threads": _split_or_default(
            metafunc.config.getoption("--a"),
            default_value=min(16, max(1, int(mp.cpu_count() / 8)))
        ),
    }

    # decide base_len (first requested arg with >1, else first non-empty)
    base_len = 1
    for n in requested:
        if len(values_by_name[n]) > 1:
            base_len = len(values_by_name[n])
            break
    if base_len == 1:
        for n in requested:
            if values_by_name[n]:
                base_len = len(values_by_name[n])
                break

    aligned_lists = [_broadcast_to_len(values_by_name[n], base_len, n) for n in requested]
    params = list(zip(*aligned_lists))
    metafunc.parametrize(tuple(requested), params, scope="session")

### ------------------- CLI input --------------------------------

# @pytest.fixture(scope="session")
# def projects(request):
#     """Fixture giving the current project path for a test run"""
#     return request.config.getoption("--projects")

# @pytest.fixture(scope="session")
# def projects_results(projects):
#     """Fixture giving the projects results path for a test run"""

#     baseDirectoryName = os.path.join(projects, "OrthoFinder")
#     if baseDirectoryName[-1] == os.sep:
#         baseDirectoryName += "Results_"
#     else:
#         baseDirectoryName = baseDirectoryName + os.sep + "Results_"

#     results_dir = CreateNewWorkingDirectory(
#         baseDirectoryName,
#         qDate=True,
#         search_program=None,
#         msa_program=None,
#         tree_program=None,
#         scorematrix=None,
#         gapopen=None,
#         gapextend=None,
#         extended_filename=False,
#         makedir=False
#     )

#     return results_dir

# @pytest.fixture(scope="session")
# def projects_results(request, projects):
#     """Fixture giving the projects results path for a test run"""
#     return request.config.getoption("--projects-results")

@pytest.fixture(scope="session")
def projects(request):
    return request.param

@pytest.fixture(scope="session")
def dna_projects(request):
    return request.param

@pytest.fixture(scope="session")
def testdata(request):
    return request.param

@pytest.fixture(scope="session")
def dna(request):
    return request.param

@pytest.fixture(scope="session")
def assign(request):
    """Fixture giving the additional species path for a test run"""
    return request.param

@pytest.fixture(scope="session")
def expected_results(request):
    """Fixture giving the expected results path for a test run"""
    return request.param

@pytest.fixture(scope="session")
def cluster(request):
    """Fixture giving the clustering program used (e.g., mcl, etc.)"""
    return request.param

@pytest.fixture(scope="session")
def msa(request):
    """Fixture giving the MSA program (e.g., famsa, mafft, etc.)"""
    return request.config.getoption("--msa")

@pytest.fixture(scope="session")
def gene_tree_method(request):
    """Fixture giving the gene tree inference program (e.g., fasttree, iqtree, etc.)"""
    return request.param

@pytest.fixture(scope="session")
def dendroblast(request):
    """Fixture indicating whether dendroblast (non-MSA) mode is enabled"""
    return request.param

@pytest.fixture(scope="session")
def species_tree_method(request):
    """Fixture giving the species tree inference program (e.g., astral-pro, etc.)"""
    return request.param

@pytest.fixture(scope="session")
def species_tree(request):
    """Fixture giving the species tree file for core"""
    return request.param

@pytest.fixture(scope="session")
def species_tree_assign(request):
    """Fixture giving the species tree file for assign"""
    return request.param

@pytest.fixture(scope="session")
def recon_method(request):
    """Fixture giving the reconciliation program (e.g., of_recon, etc.)"""
    return request.param

@pytest.fixture(scope="session")
def sequence_search_threads(request):
    """Fixture giving the number of parallel sequence search threads"""
    return request.config.getoption("--t")

@pytest.fixture(scope="session")
def analysis_threads(request):
    """Fixture giving the number of analysis threads"""
    return request.config.getoption("--a")

# @pytest.fixture(scope="session")
# def dna_assign(testdata, assign):
#     if testdata:
#         return DEFAULT_PROJECTS["testdata"]["dna"]["assign"]
    
#     else:
#         return assign

### -----------------------------------------------------------------

# @pytest.fixture(scope="session")
# def project_overall_stats(projects_results):
#     return os.path.join(
#         helper.create_path(projects_results),
#         "Comparative_Genomics_Statistics",
#         "Statistics_Overall.tsv"
#     )

@pytest.fixture(scope="session")
def expected_overall_stats(expected_results):
    return os.path.join(
        helper.create_path(expected_results),
        "Statistics_Overall.tsv"
    )

# @pytest.fixture(scope="session")
# def orthogroups(projects_results):
#     return os.path.join(
#         helper.create_path(projects_results),
#         "Orthogroups",
#         "Orthogroups.txt"
#     )

