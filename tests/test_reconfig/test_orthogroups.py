import pytest

@pytest.mark.regression
def test_orthogroups_left_diff(orthogroup_stats):
    left_diff_num = orthogroup_stats["left_diff"]
    left_label, right_label = orthogroup_stats["labels"]
    print()
    print(f"Test if {left_label} has distinct orthogroups from {right_label}.")
    assert left_diff_num == 0, f"The difference between {left_label} and {rihgt_label} is not empty!"

@pytest.mark.regression
def test_orthogroups_right_diff(orthogroup_stats):
    right_diff_num = orthogroup_stats["right_diff"]
    left_label, right_label = orthogroup_stats["labels"]
    print()
    print(f"Test if {right_label} has distinct orthogroups from {left_label}.")
    assert right_diff_num == 0, f"The difference between {right_label} and {left_label} is not empty!"

@pytest.mark.regression
def test_orthogroups_left_intersection(orthogroup_stats):
    intersection_num = orthogroup_stats["intersection"]
    left_diff_num = orthogroup_stats["left_diff"]
    left_num = intersection_num + left_diff_num
    left_label, right_label = orthogroup_stats["labels"]
    print()
    print(f"Test if the number of orthogroups in {left_label} is equal to the intersection with {right_label}.")
    assert intersection_num == left_num, f"The intersection is not equal to {left_label}!"

@pytest.mark.regression
def test_orthogroups_right_intersection(orthogroup_stats):
    intersection_num = orthogroup_stats["intersection"]
    right_diff_num = orthogroup_stats["right_diff"]
    right_num = intersection_num + right_diff_num
    left_label, right_label = orthogroup_stats["labels"]
    print()
    print(f"Test if the number of orthogroups in {right_label} is equal to the intersection with {left_label}.")
    assert intersection_num == right_num, f"The intersection is not equal to {right_label}!"

@pytest.mark.regression
def test_orthogroups_similarity(orthogroup_stats):
    similarity = orthogroup_stats["similarity"]
    left_label, right_label = orthogroup_stats["labels"]
    print()
    print(f"Test if the similarity of the {left_label} and {right_label} falls into the range (99.0, 100.0)")
    print(f"Similarity {left_label} vs. {right_label}: {similarity}")
    print()
    assert (similarity > 99.0 and similarity <= 100.0), f"The similarity between {left_label} and {right_label} is outside of the range (99.0, 100.0)!"
