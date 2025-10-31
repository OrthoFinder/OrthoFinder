import os
import pytest
from read_results import get_expected_results
from OF_funcs import OrthoFinderTestFuncs

@pytest.fixture(scope="class")
def expected_results_dict(expected_results):
    return get_expected_results(expected_results)


@pytest.fixture
def of_obj(projects, baseline_arg_dict, baseline_options, species_tree):
    return OrthoFinderTestFuncs(
        projects, 
        baseline_arg_dict,
        baseline_options, 
        species_tree=species_tree, 
    )

@pytest.fixture
def of_obj_assign(projects, baseline_arg_dict, baseline_options, assign, species_tree, species_tree_assign):
    return OrthoFinderTestFuncs(
        projects, 
        baseline_arg_dict,
        baseline_options, 
        assign,  
        species_tree, 
        species_tree_assign,
    )

@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_stats_overall(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["statistics_overall"]["stats_overall"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["statistics_overall"]["stats_overall"]

@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_ogs_nspecies(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["statistics_overall"]["og_nspecies"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["statistics_overall"]["og_nspecies"]


@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_recon_tree(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["resolved_gene_trees"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["resolved_gene_trees"]


@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_orthologues(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["orthologues"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["orthologues"]


@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_duplications(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["duplications"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["duplications"]


@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_orthogroups(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["orthogroups"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["orthogroups"]


@pytest.fixture(params=["core", "assign"], ids=["core", "assign"])
def of_and_expected_hogs(request, of_obj, of_obj_assign, expected_results_dict):
    if request.param == "core":
        return of_obj, expected_results_dict["core"]["hogs"]
    else:
        return of_obj_assign, expected_results_dict["assign"]["hogs"]
