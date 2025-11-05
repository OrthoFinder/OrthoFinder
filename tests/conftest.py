import os
import sys
import pytest
import copy
from collections import deque 
import helper
from typing import Optional, Union, List, Any
import multiprocessing as mp
from orthofinder.utils.util import CreateNewWorkingDirectory
from orthofinder.utils import files

## ---------------- TEST OrthoFinder Commands ----------------
OF_BASELINE_OPTIONS = [
    "famsa_species_tree", 
    "famsa_species_tree_assign"
]

## ------- Default ExampleData ---------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) # project root
OF_USER_CONFIG = os.path.join(ROOT, "user_config.json")
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
COMMANDS_CONFIG = os.path.abspath(os.path.join(os.path.dirname(__file__), "of_run_commands.json")) # tests root


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
        "--run-all",
        action="store_true",
        default=False,  
        help="Decide wether or not to test all available OrthoFinder command options.",
    )
    
    parser.addoption(
        "--run-options",
        action="store",
        default=None,
        help="Comma-separated list of available OrthoFinder runs command options.",
    )
    
    parser.addoption(
        "--skip-runs",
        action="store",
        default=None,
        help="Comma-separated list of available OrthoFinder runs command options.",
    )
    
    parser.addoption(
        "--user-testconfig",
        action="store",
        default=None,
        help="OrthoFinder runs command options json file.",
    )

    parser.addoption(
        "--user-ofconfig",
        action="store",
        default=OF_USER_CONFIG,
        help="OrthoFinder external software command json file.",
    )
    

    parser.addoption(
        "--projects",
        action="store",
        default=DEFAULT_PROJECTS["testdata"]["proteome"]["core"],  
        help="Comma-separated list of project paths (defaults to latest ExampleData/OrthoFinder run).",
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

### ------------------- CLI input --------------------------------
@pytest.fixture(scope="session")
def projects(request):
    return request.config.getoption("--projects")

@pytest.fixture(scope="session")
def user_ofconfig(pytestconfig):
    return pytestconfig.getoption("--user-ofconfig")

@pytest.fixture(scope="session")
def dna_projects(request):
    return request.config.getoption("--dna-projects")

@pytest.fixture(scope="session")
def assign(request):
    """Fixture giving the additional species path for a test run"""
    return request.config.getoption("--assign")

@pytest.fixture(scope="session")
def expected_results(request):
    """Fixture giving the expected results path for a test run"""
    return request.config.getoption("--expected-results")

@pytest.fixture(scope="session")
def species_tree(request):
    """Fixture giving the species tree file for core"""
    return request.config.getoption("--species-tree")

@pytest.fixture(scope="session")
def species_tree_assign(request):
    """Fixture giving the species tree file for assign"""
    return request.config.getoption("--species-tree-assign")


@pytest.fixture(scope="session")
def baseline_options():
    return OF_BASELINE_OPTIONS 



def _load_of_command_dict_from_cli(config):
    user_config = config.getoption("--user-testconfig")
    run_all = bool(config.getoption("--run-all"))
    run_opts_s = config.getoption("--run-options") or ""
    run_opts = [x for x in run_opts_s.split(",") if x]
    skip_runs = config.getoption("--skip-runs") or ""
    skip_run_opts = [x for x in skip_runs.split(",") if x]
    

    of_commands = helper.read_config_file(COMMANDS_CONFIG, user_config)

    if run_all:
        return of_commands

    base = {cat: {n: cmd for n, cmd in cmds.items() if n in OF_BASELINE_OPTIONS}
            for cat, cmds in of_commands.items()}
    if run_opts:
        sel = {cat: {n: cmd for n, cmd in cmds.items() if n in run_opts}
               for cat, cmds in of_commands.items()}
        for cat, cmds in sel.items():
            base.setdefault(cat, {}).update(cmds)
    if skip_run_opts:
        sel = {cat: {n: cmd for n, cmd in cmds.items() if n not in skip_run_opts}
               for cat, cmds in of_commands.items()}
        for cat, cmds in sel.items():
            base.setdefault(cat, {}).update(cmds)
    return base

def pytest_generate_tests(metafunc):
    if "core_case" in metafunc.fixturenames or "assign_case" in metafunc.fixturenames:
        of_cmds = _load_of_command_dict_from_cli(metafunc.config)

        core = [(name, arg)  for cat, cmds in of_cmds.items()
                               if "assign" not in cat.lower()
                               for name, arg in cmds.items()]
        assign = [(name, arg) for cat, cmds in of_cmds.items()
                               if "core" not in cat.lower()
                               for name, arg in cmds.items()]

        if "core_case" in metafunc.fixturenames:
            metafunc.parametrize("core_case", core, ids=[f"core::{n}" for n, _ in core], scope="session")
        if "assign_case" in metafunc.fixturenames:
            metafunc.parametrize("assign_case", assign, ids=[f"assign::{n}" for n, _ in assign], scope="session")


@pytest.fixture(scope="session")
def baseline_arg_dict(baseline_options):
    of_commands = helper.read_config_file(COMMANDS_CONFIG)
    of_args_dict ={}
    for k1, v1 in of_commands.items():
        for k2, v2 in v1.items(): 
            if k2 in baseline_options:
                if len(v2) != 0 and isinstance(v2, str):
                    args = v2.strip().split()
                    if "-A" in args:
                        
                        of_args_dict["msa"] = args[args.index("-A") + 1]
                        
                    elif "-M" in args:
                        of_args_dict["msa"] = None 
                        
                    if "-T" in args:
                        of_args_dict["gene_tree_method"] = args[args.index("-T") + 1]

                    if "-ST" in args:
                        of_args_dict["species_tree_method"] = args[args.index("-ST") + 1]
    return of_args_dict
                    

@pytest.fixture(autouse=True)
def reset_orthofinder_state():
    """
    Reset OrthoFinder's FileHandler global state between tests to avoid
    'Changing WorkingDirectory1' errors when running multiple CLI calls
    in-process.
    """
    attrs_to_clear = ("wd_base", "wd1", "wd2", "base_dir")
    for attr in attrs_to_clear:
        if hasattr(files.FileHandler, attr):
            setattr(files.FileHandler, attr, "")
    yield 
    for attr in attrs_to_clear:
        if hasattr(files.FileHandler, attr):
            setattr(files.FileHandler, attr, "")
            
            
        
