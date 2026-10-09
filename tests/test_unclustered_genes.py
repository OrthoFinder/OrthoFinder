"""The final MCL orthogroups (-og, --no-fix-files) account for every gene."""
from orthofinder.mcl_orthogroups.infer_orthogroups import add_unclustered_genes


def test_every_gene_clustered_nothing_added():
    ogs = [{"0_0", "1_0"}, {"0_1"}, {"1_1"}]
    assert add_unclustered_genes(ogs, [0, 1], {0: 2, 1: 2}) == 0
    assert ogs == [{"0_0", "1_0"}, {"0_1"}, {"1_1"}]


def test_unclustered_genes_added_at_the_end_in_order():
    # e.g. --assign: genes of new species 4 and 5 in no cluster
    ogs = [{"0_0", "4_0"}, {"0_1"}, {"5_1", "4_1"}]
    assert add_unclustered_genes(ogs, [0, 4, 5], {0: 2, 4: 3, 5: 2, 9: 5}) == 2
    # the orthogroups before keep their positions (their OG numbers); species 9 is not in the analysis
    assert ogs == [{"0_0", "4_0"}, {"0_1"}, {"5_1", "4_1"}, {"4_2"}, {"5_0"}]
