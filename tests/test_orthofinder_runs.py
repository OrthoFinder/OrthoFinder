import os
import re
from collections import deque
import pytest
import helper


@pytest.mark.order(1)
def test_orthofinder_direct_runs(direct_case, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd):
    name, argstr = direct_case
    helper.run_orthofinder_core_case(name, argstr, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd)

@pytest.mark.order(2)
def test_orthofinder_stopafter_runs(stopafter_case, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd):
    name, argstr = stopafter_case
    helper.run_orthofinder_core_case(name, argstr, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd)

@pytest.mark.order(3)
def test_orthofinder_restart_runs(restart_case, projects, example_data, assign, species_tree_assign, user_ofconfig, capfd):
    name, argstr = restart_case
    helper.run_orthofinder_assign_case(name, argstr, projects, example_data, assign, species_tree_assign, user_ofconfig, capfd)

@pytest.mark.order(4)
def test_orthofinder_core_runs(core_case, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd):
    name, argstr = core_case
    helper.run_orthofinder_core_case(name, argstr, projects, example_data, dna_projects, species_tree, user_ofconfig, capfd)

@pytest.mark.order(5)
def test_orthofinder_assign_runs(assign_case, projects, example_data, assign, species_tree_assign, user_ofconfig, capfd):
    name, argstr = assign_case
    helper.run_orthofinder_assign_case(name, argstr, projects, example_data, assign, species_tree_assign, user_ofconfig, capfd)
