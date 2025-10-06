import os
import pytest
from read_stats_info import get_overall_stats_info, get_og_nspecies_info


from orthofinder.utils import fasta_processor
from orthofinder.utils import files, util
from orthofinder.run.process_args import Options
from orthofinder.file_updates.ogs import OrthoGroupsSet
from orthofinder.file_updates.file_updates import read_hog_file
from orthofinder.run.species_info import SpeciesNameDict

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
def sequence_info_obj(projects, projects_results):
    options = Options()
    working_dir = os.path.join(projects_results, "WorkingDirectory")
    species_id_fn = os.path.join(working_dir, "SpeciesIDs.txt")
    sequence_id_fn = os.path.join(working_dir, "SequenceIDs.txt")
    speciesInfoObj = None
    is_dna = False

    speciesInfoObj = fasta_processor.ProcessesNewFasta(
        projects, 
        is_dna, 
        species_id_fn=species_id_fn,
        sequence_id_fn=sequence_id_fn,
        working_dir=working_dir
        )
    sequenceInfoObj = util.GetSeqsInfo([working_dir], speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
    
    
    return sequenceInfoObj

@pytest.fixture
def ogset_obj(projects, projects_results, msa, gene_tree):

    options = Options()
    options.qStartFromFasta = True
    options.msa_program = msa
    options.tree_program = gene_tree
    working_dir = os.path.join(projects_results, "WorkingDirectory")
    species_id_fn = os.path.join(working_dir, "SpeciesIDs.txt")
    sequence_id_fn = os.path.join(working_dir, "SequenceIDs.txt")
    ogs_all_fn = os.path.join(working_dir, "OGsAll.tsv")

    speciesInfoObj = None
    is_dna = False
    # files.InitialiseFileHandler(
    #         options, 
    #         fastaDir=projects, 
    #         continuationDir=None, 
    #         resultsDir_nonDefault=None, 
    #         pickleDir_nonDefault=None,
    #         working_dir=working_dir,
    #     )
    
    speciesInfoObj = fasta_processor.ProcessesNewFasta(
        projects, 
        is_dna, 
        species_id_fn=species_id_fn,
        sequence_id_fn=sequence_id_fn,
        working_dir=working_dir
        )
    sequenceInfoObj = util.GetSeqsInfo([working_dir], speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
    
    ogSet = OrthoGroupsSet(
        options.min_seq,
        [working_dir], 
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        options.tree_program,
        idExtractor=util.FirstWordExtractor,
        species_id_fn=species_id_fn,
        sequence_id_fn=sequence_id_fn,
        ogs_all_fn=ogs_all_fn,
        results_dir=projects_results,
    )

    return ogSet

@pytest.fixture
def species_name_dict(projects_results):
    working_dir = os.path.join(projects_results, "WorkingDirectory")
    species_id_fn = os.path.join(working_dir, "SpeciesIDs.txt")
    return SpeciesNameDict(species_id_fn)


# @pytest.fixture
# def hogs(projects_results):
#     hog_dir = os.path.join(projects_results, "Phylogenetic_Hierarchical_Orthogroups")
#     hog_dict = {}
#     for file in os.listdir(hog_dir):
#         node = file.rsplit(".", 1)[0]
#         file_path = os.path.join(hog_dir, file)
#         hog_list = read_hog_file(file_path)
#         hog_dict[node] = hog_list