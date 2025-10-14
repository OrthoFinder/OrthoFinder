import os 
import pytest
import helper
from orthofinder.run.__main__ import main

@pytest.mark.order(1)
def test_orthofinder_dendroblast_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "-M", "dendroblast"]

    # if msa != "famsa":
    #     args += ["-A", msa]

    # if dendroblast:
    #     if "-A" in args:
    #         args.remove("-A")
    #     if msa in args:
    #         args.remove(msa)
    #     args += ["-M", "dendroblast"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)

@pytest.mark.order(2)
def test_orthofinder_mafft_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "-A", "mafft", "-t", "2", "-a", "4"]

    # if msa != "famsa":
    #     args += ["-A", msa]

    # if dendroblast:
    #     if "-A" in args:
    #         args.remove("-A")
    #     if msa in args:
    #         args.remove(msa)
    #     args += ["-M", "dendroblast"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)
    
 
@pytest.mark.order(3)
def test_orthofinder_assign_run(projects, assign, gene_tree):
    results_dir = os.path.join(projects, "OrthoFinder")
    projects_results = helper._latest_output_dir(results_dir)
    args = ["--core", projects_results, "--assign",  assign, "-nk", 2]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)


@pytest.mark.order(4)
def test_orthofinder_famsa_species_tree_run(projects, gene_tree, species_tree):
    args = ["-f", os.path.basename(projects)]
    # if msa != "famsa":
    #     args += ["-A", msa]

    # if dendroblast:
    #     if "-A" in args:
    #         args.remove("-A")
    #     if msa in args:
    #         args.remove(msa)
    #     args += ["-M", "dendroblast"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    if gene_tree != "fasttree":
        args += ["-s", species_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)

@pytest.mark.order(5)
def test_orthofinder_prepare_blast_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "-op"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)


# @pytest.mark.order(2)
# def test_orthofinder_start_from_blast_run(projects, gene_tree):

#     results_dir = os.path.join(projects, "OrthoFinder")
#     projects_results = helper._latest_output_dir(results_dir)
#     args = ["-f", os.path.basename(projects), "-b",  projects_results]

#     if gene_tree != "fasttree":
#         args += ["-T", gene_tree]

#     with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
#         main(args)
#     code = getattr(e.value, "code", 1)
#     assert code in (0, None)

@pytest.mark.order(6)
def test_orthofinder_old_version_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "--old-version"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)

@pytest.mark.order(7)
def test_orthofinder_no_file_fix_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "--no-fix-files"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)

@pytest.mark.order(8)
def test_orthofinder_blast_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects), "-S", "blast"]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)


@pytest.mark.order(9)
def test_orthofinder_famsa_run(projects, gene_tree):
    args = ["-f", os.path.basename(projects) + os.sep]

    if gene_tree != "fasttree":
        args += ["-T", gene_tree]

    with pytest.raises((SystemExit, RuntimeError, ValueError)) as e:
        main(args)
    code = getattr(e.value, "code", 1)
    assert code in (0, None)


    

# @pytest.mark.order(5)
# def test_project_dir(projects):
#     assert os.path.exists(helper.create_path(projects))

# @pytest.mark.order(4)
# def test_project_results_dir(projects_results):
#     assert os.path.exists(helper.create_path(projects_results))

# @pytest.mark.order(5)
# def test_expected_results_dirs(expected_results):
#     assert os.path.exists(helper.create_path(expected_results))

# @pytest.mark.order(6)
# class TestProjectResultsFiles:
#     def test_stats_file(self, project_overall_stats):
#         assert os.path.exists(project_overall_stats)

#     def test_orthogroups_file(self, orthogroups):
#         assert os.path.exists(orthogroups)

# @pytest.mark.order(6)
# class TestExapectedResultsFiles:
#     def test_stats_file(self, expected_overall_stats):
#         assert os.path.exists(expected_overall_stats)
