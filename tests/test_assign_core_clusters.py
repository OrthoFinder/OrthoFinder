from orthofinder.mcl_orthogroups import assign


def test_core_clusters_file_ignores_clade_clusters(tmp_path):
    # the working directory of an --assign analysis: its orthogroups and those of a new clade
    for name in ("clusters_OrthoFinder_unassigned_clade_0_I1.2.txt_id_pairs.txt",
                 "clusters_OrthoFinder_I1.2.txt_id_pairs.txt",
                 "clusters_OrthoFinder_unassigned_clade_1_I1.2.txt_id_pairs.txt"):
        (tmp_path / name).write_text("")
    wd = str(tmp_path) + "/"
    assert assign.core_clusters_file(wd) == wd + "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"
