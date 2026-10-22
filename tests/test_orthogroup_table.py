from orthofinder.mcl_orthogroups.mcl import orthogroup_table_cells


def test_orthogroup_table_cells():
    # species 7, 3, 12 in that column order; the cells sorted, "" for a species without genes
    column = {7: 0, 3: 1, 12: 2}
    names = {"3_1": "b", "3_0": "c", "12_4": "a", "3_2": "a"}
    name_cells, id_cells, counts = orthogroup_table_cells(
        {"3_1", "3_0", "12_4", "3_2"}, names, column, 3
    )
    assert name_cells == ["", "a, b, c", "a"]
    assert id_cells == ["", "3_0, 3_1, 3_2", "12_4"]
    assert counts == [0, 3, 1]



def test_species_tree_orthogroups_none():
    # No orthogroup single-copy in enough species: none chosen (msa_trees then stops
    # before the alignments, advising -s or -ST); lowering the fraction only adds more
    import numpy as np
    from orthofinder.species_tree.concatenated_alignment import (
        determine_orthogroups_for_species_tree, get_orthogroup_occupancy_info)

    m = np.full((30, 40), 3)
    assert determine_orthogroups_for_species_tree(m, list(range(30)))[0] == []
    m[:5, :20] = 1
    _, n_ortho = get_orthogroup_occupancy_info(m)
    assert all(a <= b for a, b in zip(n_ortho, n_ortho[1:]))
