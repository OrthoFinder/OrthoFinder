
import pytest
import helper

@pytest.mark.order(23)
def test_recon_tree(of_and_expected_recon_tree):
    of_obj, expected_recon_tree = of_and_expected_recon_tree
    of_gene_tree_dict = of_obj.get_recon_gene_trees()
    identifier = "_RECON_"
    leaf_names_list = []
    test_gene_list = []
    for _, gene_tree in of_gene_tree_dict.items():
        leaf_names = {n.name for n in gene_tree if n.is_leaf and n.name in identifier}
        if len(leaf_names) != 0:
            leaf_names_list.append(leaf_names)
            test_gene_list.append(gene_tree)

    assert len(leaf_names_list) == 0 or len(test_gene_list) == 0, "No reconciliation tree was found"
    for leaf_names, gene_tree in zip(leaf_names_list, test_gene_list):
        assert expected_recon_tree.issubset(leaf_names), "????"
        assert test_gene_tree.check_monophyly(expected_recon_tree)[0], "????"


@pytest.mark.order(24)
def test_orthologues(of_and_expected_orthologues):
    of_obj, expected = of_and_expected_orthologues
    of_orthologues = of_obj.get_orthologues()
    per_species_triplets, per_species_locations = helper._index_of_orthologues(of_orthologues)

    missing = []   
    relocations = []  

    for species, expected_ogs in expected.items():
        assert species in per_species_triplets, f"No orthologues found for species {species}."

        for og_expected, entries in expected_ogs.items():
            for e in entries:
                if e not in per_species_triplets[species]:
                    missing.append((species, og_expected, e))
                else:
                    actual_ogs = per_species_locations[species][e]
                    if og_expected not in actual_ogs:
                        relocations.append((species, og_expected, e, sorted(actual_ogs)))
    if relocations:
        print("\n[INFO] Expected triplets found under different OG IDs:")
        for species, og_exp, e, ogs in relocations:
            print(f"{species}: expected in {og_exp}, actually in {', '.join(ogs)} — {e}")

    if missing:
        msgs = []
        for species, og_exp, e in missing:
            partner, a, b = e
            msgs.append(
                f"{species} missing triplet from expected {og_exp}: "
                f"(partner={partner}, A={sorted(a)}, B={sorted(b)})"
            )
        raise AssertionError("Some expected orthologues were not found:\n" + "\n".join(msgs))

@pytest.mark.order(24)
def test_duplications(of_and_expected_duplications):
    of_obj, expected_dup_map = of_and_expected_duplications

    actual_dup_map = of_obj.get_duplications()
    actual_triplets, actual_locs = helper._index_duplications(actual_dup_map)
    expected_triplets, expected_locs = helper._index_duplications(expected_dup_map)
    missing = [t for t in expected_triplets if t not in actual_triplets]
    if missing:
        lines = []
        for label, A, B in missing:
            lines.append(f"Missing duplication: label={label}, Genes1={sorted(A)}, Genes2={sorted(B)}")
        raise AssertionError("Some expected duplications were not found:\n" + "\n".join(lines))

    reloc = []
    for t in expected_triplets:
        exp_ogs = expected_locs.get(t, set())
        act_ogs = actual_locs.get(t, set())
        if exp_ogs and act_ogs and exp_ogs != act_ogs:
            reloc.append((t, sorted(exp_ogs), sorted(act_ogs)))
    if reloc:
        print("\n[INFO] Duplications present under different OG IDs:")
        for (label, A, B), exp_ogs, act_ogs in reloc:
            print(f"label={label}, Genes1={sorted(A)}, Genes2={sorted(B)} | expected {exp_ogs}, found {act_ogs}")


@pytest.mark.order(25)
def test_orthogroups(of_and_expected_orthogroups):
    of_obj, expected_orthogrups = of_and_expected_orthogroups
    orthogroups_dict = of_obj.get_orthogroups()
    expected_orthogroups_count = len(expected_orthogrups)
    
    actual_ogs, actual_locs = helper._index_of_orthogroups(orthogroups_dict)
    expected_ogs, expected_locs = helper._index_of_orthogroups(expected_orthogrups)

    missing = [t for t in expected_ogs if t not in actual_ogs]
    if missing:
        lines = []
        for label, A, B in missing:
            lines.append(f"Missing orthogroups: label={label}, Genes1={sorted(A)}, Genes2={sorted(B)}")
        raise AssertionError("Some expected Orthogroups were not found:\n" + "\n".join(lines))

    reloc = []
    for k, t in expected_orthogrups.items():
        exp_ogs = expected_locs.get(t, set())
        act_ogs = actual_locs.get(t, set())
        if exp_ogs and act_ogs and exp_ogs != act_ogs:
            reloc.append((t, sorted(exp_ogs), sorted(act_ogs)))
    if reloc:
        print("\n[INFO] Orthogroups present under different OG IDs:")
        for (label, A, B), exp_ogs, act_ogs in reloc:
            print(f"label={label}, Genes1={sorted(A)}, Genes2={sorted(B)} | expected {exp_ogs}, found {act_ogs}")





