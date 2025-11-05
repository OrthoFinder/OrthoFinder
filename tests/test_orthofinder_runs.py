import os
import re
from collections import deque
import pytest
import helper

@pytest.mark.order(1)
def test_orthofinder_core_runs(core_case, projects, dna_projects, species_tree, user_ofconfig, capfd):
    name, argstr = core_case
    helper.run_orthofinder_core_case(name, argstr, projects, dna_projects, species_tree, user_ofconfig, capfd)

@pytest.mark.order(2)
def test_orthofinder_assign_runs(assign_case, projects, assign, species_tree_assign, user_ofconfig, capfd):
    name, argstr = assign_case
    helper.run_orthofinder_assign_case(name, argstr, projects, assign, species_tree_assign, user_ofconfig, capfd)

