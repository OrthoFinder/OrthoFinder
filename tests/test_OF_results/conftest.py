import os
import pytest
from read_results import get_overall_stats_info, get_og_nspecies_info, get_expected_results
from OF_funcs import OrthoFinderTestFuncs
import helper


@pytest.fixture(scope="class")
def expected_results_dict(expected_results):
    return get_expected_results(expected_results)


# @pytest.fixture(scope="class")
# def project_assign_results()

# @pytest.fixture(scope="class")
# def project_overall_stats_info(project_overall_stats):
#     return get_overall_stats_info(project_overall_stats)

# @pytest.fixture(scope="class")
# def expected_overall_stats_info(expected_overall_stats):
#     return get_overall_stats_info(expected_overall_stats)

# @pytest.fixture(scope="class")
# def project_og_nspecies_info(project_overall_stats):
#     return get_og_nspecies_info(project_overall_stats)

# @pytest.fixture(scope="class")
# def expected_og_nspecies_info(expected_overall_stats):
#     return get_og_nspecies_info(expected_overall_stats)

@pytest.fixture
def of_obj(projects, msa, gene_tree_method, recon_method, sequence_search_threads, analysis_threads, species_tree, species_tree_assign):
    return OrthoFinderTestFuncs(
        projects, 
        msa, 
        gene_tree_method, 
        recon_method, 
        sequence_search_threads, 
        analysis_threads,
        species_tree=species_tree, 
    )

@pytest.fixture
def of_obj_assign(projects, msa, gene_tree_method, recon_method, sequence_search_threads, analysis_threads, assign, species_tree, species_tree_assign):
    return OrthoFinderTestFuncs(
        projects, 
        msa, 
        gene_tree_method, 
        recon_method, 
        sequence_search_threads, 
        analysis_threads,
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




# @pytest.fixture
# def hogs(projects_results):
#     hog_dir = os.path.join(projects_results, "Phylogenetic_Hierarchical_Orthogroups")
#     hog_dict = {}
#     for file in os.listdir(hog_dir):
#         node = file.rsplit(".", 1)[0]
#         file_path = os.path.join(hog_dir, file)
#         hog_list = read_hog_file(file_path)
#         hog_dict[node] = hog_list

