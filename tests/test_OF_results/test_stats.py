from collections import Counter
import pytest
import pytest_check as check
import numpy as np
from orthofinder.orthology.stats import add_unassigned_genes, orthogroups_matrix


@pytest.mark.skipif(
    "config.getoption('--skip-of-test')",
    reason="Skipping OrthoFinder results test"
)
@pytest.mark.order(12)
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
        project_n = sequence_info_obj.n_species
        check.equal(
            project_n,
            expected_n,
            f"[1] Number of species should be {expected_n}, you got {project_n}"
        )

        # Prepare orthogroup data
        ogs = add_unassigned_genes(ogset_obj.all_ogs(), ogset_obj.all_used_sequence_ids())
        all_ogs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
        proper_ogs = [og for og in all_ogs if len(og) > 1]
        iogs_proper_ogs = [iog for iog, og in enumerate(all_ogs) if len(og) > 1]
        all_genes = [g for og in all_ogs for g in og]
        all_genes_counter = Counter([g[0] for g in all_genes])
        n_genes = sum(all_genes_counter.values())
        assigned_genes_counter = Counter([g[0] for og in proper_ogs for g in og])
        n_assigned = sum(assigned_genes_counter.values())

        species_presence = [set([g[0] for g in og]) for og in proper_ogs]
        n_ogs = len(proper_ogs)
        species_specific_ogs_counter = Counter([next(iter(og_sp)) for og_sp in species_presence if len(og_sp) == 1])
        i_species = sequence_info_obj.species_to_use

        # Species specific orthogroups - gene-based
        i_species_specific_ogs = [i for i, og_sp in enumerate(species_presence) if len(og_sp) == 1]
        i_sp_specific_ogs_gene_counts = [
            sum([len(proper_ogs[iog]) for iog in i_species_specific_ogs if proper_ogs[iog][0][0] == i_sp]) for i_sp in
            i_species]

        # --- Number of orthogroups ---
        print("[2] Number of orthogroups")
        expected_ogs = int(expected_overall_stats_info["Number of orthogroups"])
        project_n = n_ogs
        check.equal(
            project_n,
            expected_ogs,
            f"[2] Number of orthogroups should be {expected_ogs}, you got {project_n}"
        )

        # --- Total number of species-specific orthogroups ---
        print("[3] Number of species-specific orthogroups")
        expected_n = int(expected_overall_stats_info["Number of species-specific orthogroups"])
        project_n = sum(species_specific_ogs_counter.values())
        check.equal(
            project_n,
            expected_n,
            f"[3] Total number of species-specific orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Total number of genes ---
        print("[4] Total number of genes")
        expected_n = int(expected_overall_stats_info["Number of genes"])
        project_n = n_genes
        check.equal(
            project_n,
            expected_n,
            f"[4] Total number of genes should be {expected_n}, you got {project_n}"
        )

        # --- Total number of genes in species-specific orthogroups ---
        print("[5] Total Number of genes in species-specific orthogroups")
        expected_n = int(expected_overall_stats_info["Number of genes in species-specific orthogroups"])
        project_n = sum(i_sp_specific_ogs_gene_counts)
        check.equal(
            project_n,
            expected_n,
            f"[5] Total number of genes in species-specific orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Number of genes in orthogroups ---
        print("[6] Total gene count in orthogroups")
        expected_n = int(expected_overall_stats_info["Number of genes in orthogroups"])
        project_n = n_assigned
        check.equal(
            project_n,
            expected_n,
            f"[6] Total number of genes in orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Number of unassigned genes ---
        print("[7] Number of unassigned genes")
        expected_n = int(expected_overall_stats_info["Number of unassigned genes"])
        project_n = n_genes - n_assigned
        check.equal(
            project_n,
            expected_n,
            f"[7] Total number of unassigned genes should be {expected_n}, you got {project_n}"
        )

        # --- Prepare percentage data ---
        percent_format = "%0.1f"
        p_assigned = 100. * n_assigned / n_genes

        # --- Percentage of genes in orthogroups ---
        print("[8] Percentage of genes in orthogroups")
        expected_n = float(expected_overall_stats_info["Percentage of genes in orthogroups"])
        project_n = float(percent_format % p_assigned)
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[8] Percentage of genes in orthogroups should be {expected_n}, you got {project_n}"
        )

        # --- Percentage of unassigned genes ---
        print("[9] Percentage of unassigned genes")
        expected_n = float(expected_overall_stats_info["Percentage of unassigned genes"])
        project_n = float(percent_format % (100 * (n_genes - n_assigned) / n_genes))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[9] Percentage of unassigned genes should be {expected_n}, you got {project_n}"
        )


        # --- Percentage of genes in species-specific orthogroups ---
        print("[10] Percentage of genes in species-specific orthogroups")
        expected_n = float(expected_overall_stats_info["Percentage of genes in species-specific orthogroups"])
        project_n = float(percent_format % (100. * sum(i_sp_specific_ogs_gene_counts) / n_genes))
        check.almost_equal(
            project_n,
            expected_n,
            rel=1e-6,
            msg=f"[10] Percentage of genes in species-specific orthogroups should be {expected_n}, you got {project_n}"
        )


        # 'averages'
        l = list(sorted(list(map(len, proper_ogs))))
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > n_assigned / 2)

        l2 = list(reversed(list(map(len, ogs))))
        L2 = np.cumsum(l2)
        j2, _ = next((i, x) for i, x in enumerate(L2) if x > n_genes / 2)
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
        og_matrix = orthogroups_matrix(i_species, proper_ogs)  # use iogs_properOGs for indexing
        n_species = len(i_species)
        n_present = (og_matrix > np.zeros((1, n_species))).sum(1)
        n_complete_ogs = list(n_present).count(n_species)
        single_copy_ogs = (og_matrix == np.ones((1, n_species))).all(1).nonzero()[0]
        n_single_copy = len(single_copy_ogs)

        # --- O50 (all genes) ---
        print("[16] Number of orthogroups with all species present")
        expected_n = int(expected_overall_stats_info["Number of orthogroups with all species present"])
        project_n = n_complete_ogs
        check.equal(
            project_n,
            expected_n,
            f"[16] Number of orthogroups with all species present should be {expected_n}, you got {project_n}"
        )

        # --- O50 (all genes) ---
        print("[17] Number of single-copy orthogroups")
        expected_n = int(expected_overall_stats_info["Number of single-copy orthogroups"])
        project_n = n_single_copy
        check.equal(
            project_n,
            expected_n,
            f"[17] Number of single-copy orthogroups should be {expected_n}, you got {project_n}"
        )

    def test_number_of_species_in_orthogroup(self, of_and_expected_ogs_nspecies):
        of_obj, expected_og_nspecies_info = of_and_expected_ogs_nspecies

        sequence_info_obj = of_obj.get_sequence_info_obj()
        ogset_obj = of_obj.get_og_obj()
        ogs = add_unassigned_genes(ogset_obj.all_ogs(), ogset_obj.all_used_sequence_ids())
        all_ogs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
        proper_ogs = [og for og in all_ogs if len(og) > 1]
        species_presence = [set([g[0] for g in og]) for og in proper_ogs]
        i_species = sequence_info_obj.species_to_use
        n_sp = len(i_species)
        # Species presence
        n = list(map(len, species_presence))
        print()
        for i in range(1, n_sp + 1):
            print(f"Test {i} species in orthogroup")
            project_n = n.count(i)
            expected_n = expected_og_nspecies_info[i]
            check.equal(
                project_n,
                expected_n,
                f"Number of orthogroups with {i} species should be {expected_n}, you got {project_n}"
            )

