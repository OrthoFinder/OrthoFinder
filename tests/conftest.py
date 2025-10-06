import os
import sys
import pytest
import helper
from typing import Optional


# HERE = Path(__file__).resolve()
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
EXAMPLEDATA = os.path.join(ROOT, "ExampleData")
EXAMPLE_RESULTS = os.path.join(ROOT, "ExampleData", "OrthoFinder")
ORTHOFINDER = os.path.join(ROOT,"src")
sys.path.insert(0, ORTHOFINDER)

EXPECTED_RESULTS = os.path.join(ROOT, "tests", "expected_results")

if EXAMPLE_RESULTS[-1] != os.sep:
    EXAMPLE_RESULTS += os.sep 

if ORTHOFINDER[-1] != os.sep:
    ORTHOFINDER += os.sep

exampledata_results_dir_list = sorted([
    (entry.stat().st_mtime, entry.path) 
    for entry in os.scandir(EXAMPLE_RESULTS)
])

latest_output_dir = exampledata_results_dir_list[-1][1]
working_dir = os.path.join(latest_output_dir, "WorkingDirectory")


def _latest_output_dir() -> str:
    try:
        entries = sorted((e.stat().st_mtime, e.path) for e in os.scandir(EXAMPLE_RESULTS))
        return entries[-1][1] if entries else ""
    except FileNotFoundError:
        return ""

DEFAULT_PROJECTS = EXAMPLEDATA
DEFAULT_PROJECTS_RESULTS = _latest_output_dir()  # used if --projects not supplied


def pytest_addoption(parser):
    parser.addoption(
        "--projects",
        action="store",
        default=DEFAULT_PROJECTS,  
        help="Comma-separated list of project paths (defaults to latest ExampleData/OrthoFinder run).",
    )

    parser.addoption(
        "--projects-results",
        action="store",
        default=DEFAULT_PROJECTS_RESULTS, 
        help="Comma-separated list of project result paths (defaults to latest ExampleData/OrthoFinder run).",
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
        "--gene-tree",
        action="store",
        default="fasttree",  
        help="Gene tree program",
    )

    parser.addoption(
        "--dendroblast",
        action="store_true",
        default=False,  
        help="Non-MSA program",
    )

    parser.addoption(
        "--species-tree",
        action="store",
        default="astral-pro",  
        help="Species tree program",
    )

def _split_csv(s: str):
    return [x for x in s.split(",") if x] if s else []

def _split_or_default(s: str, default_value: Optional[str]):
    lst = _split_csv(s)
    if lst:
        return lst
    return [default_value] if default_value else []

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
    raise pytest.UsageError(
        f"Length mismatch for '{name}': got {len(lst)} values, need {target_len}. "
        f"Provide 1 value to broadcast or {target_len} values to align."
    )

def pytest_generate_tests(metafunc):
    """
    Dynamically parameterize tests that request any combination of:
    projects, projects_results, expected_results,
    cluster, msa, gene_tree, dendroblast, species_tree
    """
    # figure out what the test actually asks for
    requested = [name for name in (
        "projects",
        "projects_results",
        "expected_results",
        "cluster",
        "msa",
        "gene_tree",
        "dendroblast",
        "species_tree"
    ) if name in metafunc.fixturenames]
    if not requested:
        return

    # gather CLI options into lists
    values_by_name = {
        "projects": _split_or_default(
            metafunc.config.getoption("--projects"), default_value=DEFAULT_PROJECTS
        ),
        "projects_results": _split_or_default(
            metafunc.config.getoption("--projects-results"), default_value=DEFAULT_PROJECTS_RESULTS
        ),
        "expected_results": _split_csv(
            metafunc.config.getoption("--expected-results")
        ),
        "cluster": _split_or_default(
            metafunc.config.getoption("--cluster"), default_value="mcl"
        ),
        "msa": _split_or_default(
            metafunc.config.getoption("--msa"), default_value="famsa"
        ),
        "gene_tree": _split_or_default(
            metafunc.config.getoption("--gene-tree"), default_value="fasttree"
        ),
        "dendroblast": _split_or_default(
            # store_true returns bool, convert to list
            str(metafunc.config.getoption("--dendroblast")),
            default_value=str(False)
        ),
        "species_tree": _split_or_default(
            metafunc.config.getoption("--species-tree"), default_value="astral-pro"
        ),
    }

    # decide how many rows to generate:
    # pick first requested arg that has >1 elements or any non-empty
    base_len = 1
    for n in requested:
        if len(values_by_name[n]) > 1:
            base_len = len(values_by_name[n])
            break
    # if all requested are singletons or empty but projects has >1, pick that
    if base_len == 1:
        for n in requested:
            if values_by_name[n]:
                base_len = len(values_by_name[n])
                break

    # align each requested value list
    aligned_lists = [ _broadcast_to_len(values_by_name[n], base_len, n) for n in requested ]

    # parametrize with a single call
    params = list(zip(*aligned_lists))
    metafunc.parametrize(tuple(requested), params, scope="session")

### ------------------- CLI input --------------------------------

@pytest.fixture(scope="session")
def projects(request):
    """Fixture giving the current project path for a test run"""
    return request.config.getoption("--projects")

@pytest.fixture(scope="session")
def projects_results(request):
    """Fixture giving the projects results path for a test run"""
    return request.config.getoption("--projects-results")

@pytest.fixture(scope="session")
def expected_results(request):
    """Fixture giving the expected results path for a test run"""
    return request.config.getoption("--expected-results")

@pytest.fixture(scope="session")
def cluster(request):
    """Fixture giving the clustering program used (e.g., mcl, etc.)"""
    return request.config.getoption("--cluster")

@pytest.fixture(scope="session")
def msa(request):
    """Fixture giving the MSA program (e.g., famsa, mafft, etc.)"""
    return request.config.getoption("--msa")

@pytest.fixture(scope="session")
def gene_tree(request):
    """Fixture giving the gene tree inference program (e.g., fasttree, iqtree, etc.)"""
    return request.config.getoption("--gene-tree")

@pytest.fixture(scope="session")
def dendroblast(request):
    """Fixture indicating whether dendroblast (non-MSA) mode is enabled"""
    return request.config.getoption("--dendroblast")

@pytest.fixture(scope="session")
def species_tree(request):
    """Fixture giving the species tree inference program (e.g., astral-pro, etc.)"""
    return request.config.getoption("--species-tree")

### -----------------------------------------------------------------

@pytest.fixture(scope="session")
def project_overall_stats(projects_results):
    return os.path.join(
        helper.create_path(projects_results),
        "Comparative_Genomics_Statistics",
        "Statistics_Overall.tsv"
    )

@pytest.fixture(scope="session")
def expected_overall_stats(expected_results):
    return os.path.join(
        helper.create_path(expected_results),
        "Statistics_Overall.tsv"
    )

@pytest.fixture(scope="session")
def orthogroups(projects_results):
    return os.path.join(
        helper.create_path(projects_results),
        "Orthogroups",
        "Orthogroups.txt"
    )

