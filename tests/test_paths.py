import os 
import pytest
import helper


def test_project_dir(projects):
    assert os.path.exists(helper.create_path(projects))

def test_project_results_dir(projects_results):
    assert os.path.exists(helper.create_path(projects_results))

def test_expected_results_dirs(expected_results):
    assert os.path.exists(helper.create_path(expected_results))

class TestProjectResultsFiles:
    def test_stats_file(self, project_overall_stats):
        assert os.path.exists(project_overall_stats)

    def test_orthogroups_file(self, orthogroups):
        assert os.path.exists(orthogroups)

class TestExapectedResultsFiles:
    def test_stats_file(self, expected_overall_stats):
        assert os.path.exists(expected_overall_stats)
