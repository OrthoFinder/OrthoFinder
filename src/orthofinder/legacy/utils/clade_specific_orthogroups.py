from ...utils import util
from ...orthogroups import accelerate as acc
from ...tools import mcl
from ...orthogroups import gathering

def clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files):
    ogs = acc.get_original_orthogroups()
    i_og_restart = len(ogs)
    ogs_new_species, _ = acc.assign_genes(results_files)
    n_unassigned = acc.write_unassigned_fasta(ogs, ogs_new_species, speciesInfoObj)
    CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=True)
    RunSearch(options, speciesInfoObj, seqsInfo, prog_caller, n_genes_per_species=n_unassigned, q_new_species_unassigned_genes=True)

    # Update info objects for clustering of only the new species
    speciesInfoObj_for_unassigned = copy.deepcopy(speciesInfoObj)
    speciesInfoObj_for_unassigned.speciesToUse = [iSp for iSp in speciesInfoObj.speciesToUse if iSp >= speciesInfoObj.iFirstNewSpecies]

    nSpecies_unassigned = len(speciesInfoObj_for_unassigned.speciesToUse)
    nSeqs_in_new_species = sum([seqsInfo.nSeqsPerSpecies[iSp] for iSp in speciesInfoObj_for_unassigned.speciesToUse])
    n_offset = seqsInfo.seqStartingIndices[-nSpecies_unassigned]
    starting_indicies = [n-n_offset for n in seqsInfo.seqStartingIndices[-nSpecies_unassigned:]]
    seqsInfo_for_unassigned = util.SequencesInfo(
        nSeqs = nSeqs_in_new_species,
        nSpecies = nSpecies_unassigned,
        speciesToUse = speciesInfoObj_for_unassigned.speciesToUse,
        seqStartingIndices = starting_indicies,
        nSeqsPerSpecies = seqsInfo.nSeqsPerSpecies
    )
    options.v2_scores = True
    clustersFilename_pairs_unassigned = gathering.DoOrthogroups(options, speciesInfoObj_for_unassigned, seqsInfo_for_unassigned,
                                                     speciesNamesDict, speciesXML=None, i_unassigned=0)
    ogs_clade_specific = mcl.GetPredictedOGs(clustersFilename_pairs_unassigned)
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, ogs_new_species, [ogs_clade_specific])
    return clustersFilename_pairs, i_og_restart