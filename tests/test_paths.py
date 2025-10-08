import os 
import pytest
import helper
from orthofinder.run.__main__ import main

@pytest.mark.order(1)
def test_orthofinder_run(projects, msa, dendroblast, gene_tree):
    args = ["-f", os.path.basename(projects)]

    if msa != "famsa":
        args += ["-A", msa]

    if dendroblast:
        if "-A" in args:
            args.remove("-A")
        if msa in args:
            args.remove(msa)
        args += ["-M", "dendroblast"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises(SystemExit) as e:
        main(args)

    # Accept success as exit code 0 or None (both mean success)
    code = e.value.code
    assert code in (0, None)
    
@pytest.mark.order(2)
def test_project_dir(projects):
    assert os.path.exists(helper.create_path(projects))

@pytest.mark.order(3)
def test_project_results_dir(projects_results):
    assert os.path.exists(helper.create_path(projects_results))

@pytest.mark.order(4)
def test_expected_results_dirs(expected_results):
    assert os.path.exists(helper.create_path(expected_results))

@pytest.mark.order(5)
class TestProjectResultsFiles:
    def test_stats_file(self, project_overall_stats):
        assert os.path.exists(project_overall_stats)

    def test_orthogroups_file(self, orthogroups):
        assert os.path.exists(orthogroups)

@pytest.mark.order(6)
class TestExapectedResultsFiles:
    def test_stats_file(self, expected_overall_stats):
        assert os.path.exists(expected_overall_stats)
