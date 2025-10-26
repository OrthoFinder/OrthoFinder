from collections import Counter
import pytest
import pytest_check as check
import numpy as np
from orthofinder.comparative_genomics.stats import add_unassigned_genes, OrthogroupsMatrix
from orthofinder.run.__main__ import main

# class TestComparativeStats:
#     def test_number_of_species(
#         self,
#         project_overall_stats_info,
#         expected_overall_stats_info,
#     ):
#         project_n = project_overall_stats_info["Number of species"]
#         expected_n = expected_overall_stats_info["Number of species"]
#         assert project_n == expected_n, (
#             f"Number of species from the project is {project_n}, "
#             f"but we expect {expected_n} species"
#         )


@pytest.mark.order(22)
class TestcComparativeStats:
    def test_overall_stats(self, of_and_expected_stats_overall):
        """
        Validate high-level statistics of an OrthoFinder run against the expected baseline.
        Uses pytest-check so all checks run even if one fails.
        """
        print("Validating high-level OrthoFinder statistics")
        of_obj, expected_overall_stats_info = of_and_expected_stats_overall

        # --- Number of species ---
        print("[1] Number of species")
        expected_n = int(expected_overall_stats_info["Number of species"])
        sequence_info_obj = of_obj.get_sequence_info_obj()
        ogset_obj = of_obj.get_og_obj()
        project_n = sequence_info_obj.nSpecies
        check.equal(
            project_n,
            expected_n,
            f"[1] Number of species should be {expected_n}, you got {project_n}"
        )

        # Prepare orthogroup data
        ogs = add_unassigned_genes(ogset_obj.AllOGs(), ogset_obj.AllUsedSequenceIDs())
        allOgs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
        properOGs = [og for og in allOgs if len(og) > 1]
        iogs_properOGs = [iog for iog, og in enumerate(allOgs) if len(og) > 1]
        allGenes = [g for og in allOgs for g in og]
        allGenesCounter = Counter([g[0] for g in allGenes])
        nGenes = sum(allGenesCounter.values())
        assignedGenesCounter = Counter([g[0] for og in properOGs for g in og])
        nAssigned = sum(assignedGenesCounter.values())

        speciesPresence = [set([g[0] for g in og]) for og in properOGs]
        nOgs = len(properOGs)
        speciesSpecificOGsCounter = Counter([next(iter(og_sp)) for og_sp in speciesPresence if len(og_sp) == 1])
        iSpecies = sequence_info_obj.speciesToUse

        # Species specific orthogroups - gene-based
        iSpeciesSpecificOGs = [i for i, og_sp in enumerate(speciesPresence) if len(og_sp) == 1]
        iSpSpecificOGsGeneCounts = [
            sum([len(properOGs[iog]) for iog in iSpeciesSpecificOGs if properOGs[iog][0][0] == iSp]) for iSp in
            iSpecies]

        # --- Number of orthogroups ---
        print("[2] Number of orthogroups")
        expected_ogs = int(expected_overall_stats_info["Number of orthogroups"])
        project_n = nOgs
        check.equal(
            project_n,
            expected_ogs,
            f"[2] Number of orthogroups should be {expected_ogs}, you got {project_n}"
        )

        # --- Total number of species-specific orthogroups ---
        print("[3] Number of species-specific orthogroups")
        expected_n = int(expected_overall_stats_info["Number of species-specific orthogroups"])
        project_n = sum(speciesSpecificOGsCounter.values())
        check.equal(
            project_n,
            expected_n,
            f"[3] Total number of species-specific orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Total number of genes ---
        print("[4] Total number of genes")
        expected_n = int(expected_overall_stats_info["Number of genes"])
        project_n = nGenes
        check.equal(
            project_n,
            expected_n,
            f"[4] Total number of genes should be {expected_n}, you got {project_n}"
        )

        # --- Total number of genes in species-specific orthogroups ---
        print("[5] Total Number of genes in species-specific orthogroups")
        expected_n = int(expected_overall_stats_info["Number of genes in species-specific orthogroups"])
        project_n = sum(iSpSpecificOGsGeneCounts)
        check.equal(
            project_n,
            expected_n,
            f"[5] Total number of genes in species-specific orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Number of genes in orthogroups ---
        print("[6] Total gene count in orthogroups")
        expected_n = int(expected_overall_stats_info["Number of genes in orthogroups"])
        project_n = nAssigned
        check.equal(
            project_n,
            expected_n,
            f"[6] Total number of genes in orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Number of unassigned genes ---
        print("[7] Number of unassigned genes")
        expected_n = int(expected_overall_stats_info["Number of unassigned genes"])
        project_n = nGenes - nAssigned
        check.equal(
            project_n,
            expected_n,
            f"[7] Total number of unassigned genes should be {expected_n}, you got {project_n}"
        )

        # --- Prepare percentage data ---
        percentFormat = "%0.1f"
        pAssigned = 100. * nAssigned / nGenes

        # --- Percentage of genes in orthogroups ---
        print("[8] Percentage of genes in orthogroups")
        expected_n = float(expected_overall_stats_info["Percentage of genes in orthogroups"])
        project_n = float(percentFormat % pAssigned)
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[8] Percentage of genes in orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Percentage of unassigned genes ---
        print("[9] Percentage of unassigned genes")
        expected_n = float(expected_overall_stats_info["Percentage of unassigned genes"])
        project_n = float(percentFormat % (100 * (nGenes - nAssigned) / nGenes))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[9] Percentage of unassigned genes should be {expected_n}, you got {project_n}"
        )


        # --- Percentage of genes in species-specific orthogroups ---
        print("[10] Percentage of genes in species-specific orthogroups")
        expected_n = float(expected_overall_stats_info["Percentage of genes in species-specific orthogroups"])
        project_n = float(percentFormat % (100. * sum(iSpSpecificOGsGeneCounts) / nGenes))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[10] Percentage of genes in species-specific orthogroups should be {expected_n}, you got {project_n}"
        )


        # 'averages'
        l = list(sorted(list(map(len, properOGs))))
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > nAssigned / 2)

        l2 = list(reversed(list(map(len, ogs))))
        L2 = np.cumsum(l2)
        j2, _ = next((i, x) for i, x in enumerate(L2) if x > nGenes / 2)
        G50 = l2[j2]
        O50 = len(l2) - j2

        # --- Mean orthogroup size ---
        print("[11] Mean orthogroup size")
        expected_n = float(expected_overall_stats_info["Mean orthogroup size"])
        project_n = float(np.round(np.mean(l), 1))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[11] Mean orthogroup size should be {expected_n}, you got {project_n}"
        )

        # --- Mean orthogroup size ---
        print("[12] Median orthogroup size")
        expected_n = float(expected_overall_stats_info["Median orthogroup size"])
        project_n = float(np.round(np.median(l), 1))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[12] Median orthogroup size should be {expected_n}, you got {project_n}"
        )

        # --- G50 (assigned genes) ---
        print("[13] G50 (assigned genes)")
        expected_n = int(expected_overall_stats_info["G50 (assigned genes)"])
        project_n = l[j]
        check.equal(
            project_n,
            expected_n,
            f"[13] G50 (assigned genes) should be {expected_n}, you got {project_n}"
        )

        # --- G50 (all genes) ---
        print("[14] G50 (assigned genes)")
        expected_n = int(expected_overall_stats_info["G50 (all genes)"])
        project_n = G50
        check.equal(
            project_n,
            expected_n,
            f"[14] G50 (all genes) should be {expected_n}, you got {project_n}"
        )

        # --- O50 (assigned genes) ---
        print("[15] O50 (assigned genes)")
        expected_n = int(expected_overall_stats_info["O50 (assigned genes)"])
        project_n = len(l) - j
        check.equal(
            project_n,
            expected_n,
            f"[15] O50 (assigned genes) should be {expected_n}, you got {project_n}"
        )

        # --- O50 (all genes) ---
        print("[15] O50 (all genes)")
        expected_n = int(expected_overall_stats_info["O50 (all genes)"])
        project_n = O50
        check.equal(
            project_n,
            expected_n,
            f"[15] O50 (all genes) should be {expected_n}, you got {project_n}"
        )

        # Single-copy orthogroups
        ogMatrix = OrthogroupsMatrix(iSpecies, properOGs)  # use iogs_properOGs for indexing
        nSpecies = len(iSpecies)
        nPresent = (ogMatrix > np.zeros((1, nSpecies))).sum(1)
        nCompleteOGs = list(nPresent).count(nSpecies)
        singleCopyOGs = (ogMatrix == np.ones((1, nSpecies))).all(1).nonzero()[0]
        nSingleCopy = len(singleCopyOGs)

        # --- O50 (all genes) ---
        print("[16] Number of orthogroups with all species present")
        expected_n = int(expected_overall_stats_info["Number of orthogroups with all species present"])
        project_n = nCompleteOGs
        check.equal(
            project_n,
            expected_n,
            f"[16] Number of orthogroups with all species present should be {expected_n}, you got {project_n}"
        )

        # --- O50 (all genes) ---
        print("[17] Number of single-copy orthogroups")
        expected_n = int(expected_overall_stats_info["Number of single-copy orthogroups"])
        project_n = nSingleCopy
        check.equal(
            project_n,
            expected_n,
            f"[17] Number of single-copy orthogroups should be {expected_n}, you got {project_n}"
        )

    def test_number_of_species_in_orthogroup(self, of_and_expected_ogs_nspecies):
        of_obj, expected_og_nspecies_info = of_and_expected_ogs_nspecies

        sequence_info_obj = of_obj.get_sequence_info_obj()
        ogset_obj = of_obj.get_og_obj()
        ogs = add_unassigned_genes(ogset_obj.AllOGs(), ogset_obj.AllUsedSequenceIDs())
        allOgs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
        properOGs = [og for og in allOgs if len(og) > 1]
        speciesPresence = [set([g[0] for g in og]) for og in properOGs]
        iSpecies = sequence_info_obj.speciesToUse
        nSp = len(iSpecies)
        # Species presence
        n = list(map(len, speciesPresence))
        print()
        for i in range(1, nSp + 1):
            print(f"Test {i} species in orthogroup")
            project_n = n.count(i)
            expected_n = expected_og_nspecies_info[i]
            check.equal(
                project_n,
                expected_n,
                f"Number of orthogroups with {i} species should be {expected_n}, you got {project_n}"
            )

