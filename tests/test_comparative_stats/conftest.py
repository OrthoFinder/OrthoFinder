import os
import pytest
from read_stats_info import get_overall_stats_info, get_og_nspecies_info
from OF_funcs import OrthoFinderTestFuncs


@pytest.fixture(scope="class")
def project_overall_stats_info(project_overall_stats):
    return get_overall_stats_info(project_overall_stats)

@pytest.fixture(scope="class")
def expected_overall_stats_info(expected_overall_stats):
    return get_overall_stats_info(expected_overall_stats)

@pytest.fixture(scope="class")
def project_og_nspecies_info(project_overall_stats):
    return get_og_nspecies_info(project_overall_stats)

@pytest.fixture(scope="class")
def expected_og_nspecies_info(expected_overall_stats):
    return get_og_nspecies_info(expected_overall_stats)

@pytest.fixture
def of_obj(projects, msa, gene_tree_method, recon_method, sequence_search_threads, analysis_threads):
    return OrthoFinderTestFuncs(projects, msa, gene_tree_method, recon_method, sequence_search_threads, analysis_threads)

# @pytest.fixture
# def hogs(projects_results):
#     hog_dir = os.path.join(projects_results, "Phylogenetic_Hierarchical_Orthogroups")
#     hog_dict = {}
#     for file in os.listdir(hog_dir):
#         node = file.rsplit(".", 1)[0]
#         file_path = os.path.join(hog_dir, file)
#         hog_list = read_hog_file(file_path)
#         hog_dict[node] = hog_list