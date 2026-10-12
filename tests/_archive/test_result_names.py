"""Continuing in a results directory of OrthoFinder before 3.2 renames its "Orthologue" names."""
from orthofinder.utils.file_handler import rename_earlier_result_names


def make_old_results(d):
    (d / "Orthologues" / "Orthologues_Sp_A").mkdir(parents=True)
    (d / "Orthologues" / "Orthologues_Sp_A" / "Sp_A__v__Sp_B.tsv").write_text("pairs")
    (d / "Orthologues" / "Sp_A.tsv.gz").write_text("orthologs")
    (d / "Single_Copy_Orthologue_Sequences").mkdir()
    (d / "Single_Copy_Orthologue_Sequences" / "OG0000001.fa").write_text(">a")
    (d / "Comparative_Genomics_Statistics").mkdir()
    for name in ("OrthologuesStats_Totals.tsv", "OrthologuesStats_one-to-one.tsv", "Statistics_Overall.tsv"):
        (d / "Comparative_Genomics_Statistics" / name).write_text(name)
    (d / "Orthogroups").mkdir()
    (d / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt").write_text("OG0000001")
    (d / "Orthogroups" / "Orthogroups.tsv").write_text("ogs")


def all_files(d):
    return {str(p.relative_to(d)): p.read_text() for p in d.rglob("*") if p.is_file()}


def test_old_names_renamed(tmp_path):
    make_old_results(tmp_path)
    renamed = rename_earlier_result_names(str(tmp_path))
    assert len(renamed) == 6
    assert all_files(tmp_path) == {
        "Orthologs/Orthologs_Sp_A/Sp_A__v__Sp_B.tsv": "pairs",
        "Orthologs/Sp_A.tsv.gz": "orthologs",
        "Single_Copy_Ortholog_Sequences/OG0000001.fa": ">a",
        "Comparative_Genomics_Statistics/OrthologsStats_Totals.tsv": "OrthologuesStats_Totals.tsv",
        "Comparative_Genomics_Statistics/OrthologsStats_one-to-one.tsv": "OrthologuesStats_one-to-one.tsv",
        "Comparative_Genomics_Statistics/Statistics_Overall.tsv": "Statistics_Overall.tsv",
        "Orthogroups/Orthogroups_SingleCopyOrthologs.txt": "OG0000001",
        "Orthogroups/Orthogroups.tsv": "ogs",
    }
    assert rename_earlier_result_names(str(tmp_path)) == []  # nothing left to rename


def test_new_names_not_overwritten(tmp_path):
    make_old_results(tmp_path)
    (tmp_path / "Orthologs").mkdir()
    (tmp_path / "Orthologs" / "Sp_A.tsv.gz").write_text("new")
    rename_earlier_result_names(str(tmp_path))
    assert (tmp_path / "Orthologs" / "Sp_A.tsv.gz").read_text() == "new"
    assert (tmp_path / "Orthologues" / "Sp_A.tsv.gz").read_text() == "orthologs"  # left as it was


def test_new_results_unchanged(tmp_path):
    (tmp_path / "Orthologs").mkdir()
    (tmp_path / "Orthogroups").mkdir()
    (tmp_path / "Orthogroups" / "Orthogroups.tsv").write_text("ogs")
    assert rename_earlier_result_names(str(tmp_path)) == []
    assert all_files(tmp_path) == {"Orthogroups/Orthogroups.tsv": "ogs"}
