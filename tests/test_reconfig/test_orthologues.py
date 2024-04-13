import pytest

@pytest.mark.regression
def test_orthologues_left_diff(orthologues_stats):
    left_diff_num = orthologues_stats["left_diff"]
    left_label, right_label = orthologues_stats["labels"]
    print()
    print(f"Test if {left_label} has distinct orthologues from {right_label}.")
    assert left_diff_num == 0, f"The difference between {left_label} and {right_label} is not empty!"

@pytest.mark.regression
def test_orthologues_right_diff(orthologues_stats):
    right_diff_num = orthologues_stats["right_diff"]
    left_label, right_label = orthologues_stats["labels"]
    print()
    print(f"Test if {right_label} has distinct orthologues from {left_label}.")
    assert right_diff_num == 0, f"The difference between {right_label} and {left_label} is not empty!"

@pytest.mark.regression
def test_orthologues_left_intersection(orthologues_stats):
    intersection_num = orthologues_stats["intersection"]
    left_diff_num = orthologues_stats["left_diff"]
    left_num = intersection_num + left_diff_num
    left_label, right_label = orthologues_stats["labels"]
    print()
    print(f"Test if the number of orthologues in {left_label} is equal to the intersection with {right_label}.")
    assert intersection_num == left_num, f"The intersection is not equal to {left_label}!"

@pytest.mark.regression
def test_orthologues_right_intersection(orthologues_stats):
    intersection_num = orthologues_stats["intersection"]
    right_diff_num = orthologues_stats["right_diff"]
    right_num = intersection_num + right_diff_num
    left_label, right_label = orthologues_stats["labels"]
    print()
    print(f"Test if the number of orthologues in {right_label} is equal to the intersection with {left_label}.")
    assert intersection_num == right_num, f"The intersection is not equal to {right_label}!"