import pytest

@pytest.mark.regression
def test_unassigned_gens_left_diff(unassigned_gens_stats):
    left_diff_num = unassigned_gens_stats["left_diff"]
    left_label, right_label = unassigned_gens_stats["labels"]
    print()
    print(f"Test if {left_label} has distinct unassigned_gens from {right_label}.")
    assert left_diff_num == 0, f"The difference between {left_label} and {rihgt_label} is not empty!"

@pytest.mark.regression
def test_unassigned_gens_right_diff(unassigned_gens_stats):
    right_diff_num = unassigned_gens_stats["right_diff"]
    left_label, right_label = unassigned_gens_stats["labels"]
    print()
    print(f"Test if {right_label} has distinct unassigned_gens from {left_label}.")
    assert right_diff_num == 0, f"The difference between {right_label} and {left_label} is not empty!"

@pytest.mark.regression
def test_unassigned_gens_left_intersection(unassigned_gens_stats):
    intersection_num = unassigned_gens_stats["intersection"]
    left_diff_num = unassigned_gens_stats["left_diff"]
    left_num = intersection_num + left_diff_num
    left_label, right_label = unassigned_gens_stats["labels"]
    print()
    print(f"Test if the number of unassigned_gens in {left_label} is equal to the intersection with {right_label}.")
    assert intersection_num == left_num, f"The intersection is not equal to {left_label}!"

@pytest.mark.regression
def test_unassigned_gens_right_intersection(unassigned_gens_stats):
    intersection_num = unassigned_gens_stats["intersection"]
    right_diff_num = unassigned_gens_stats["right_diff"]
    right_num = intersection_num + right_diff_num
    left_label, right_label = unassigned_gens_stats["labels"]
    print()
    print(f"Test if the number of unassigned_gens in {right_label} is equal to the intersection with {left_label}.")
    assert intersection_num == right_num, f"The intersection is not equal to {right_label}!"