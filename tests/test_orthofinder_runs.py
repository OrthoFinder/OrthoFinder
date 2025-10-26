import os
import re
import pytest
import helper

# # if msa != "famsa":
# #     args += ["-A", msa]

# # if dendroblast:
# #     if "-A" in args:
# #         args.remove("-A")
# #     if msa in args:
# #         args.remove(msa)
# #     args += ["-M", "dendroblast"]

# # if gene_tree_method != "fasttree":
# #     args += ["-T", gene_tree_method]


class TestOrthoFinderRuns:

    # @pytest.mark.order(1)
    # def test_orthofinder_famsa_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(2)
    # def test_orthofinder_assign_run(self, projects, assign, gene_tree_method, capfd):
    #     results_dir = os.path.join(projects, "OrthoFinder")
    #     projects_results = helper._latest_output_dir(results_dir)
    #     args = ["--core", projects_results, "--assign",  assign, "-nk", 2]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text

    # @pytest.mark.order(3)
    # def test_orthofinder_dendroblast_run(self, projects, gene_tree_method, capfd):

    #     args = ["-f", projects, "-M", "dendroblast"]
        
    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
        
    #     assert not any(re.search(p, text) for p in fatal_patterns), text

    # @pytest.mark.order(4)
    # def test_orthofinder_mafft_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects, "-A", "mafft", "-t", "2", "-a", "4"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(5)
    # def test_orthofinder_prepare_blast_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects, "-op", "--save-blast-commands"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(6)
    # def test_orthofinder_start_from_blast_run(self, projects, gene_tree_method, capfd):

    #     results_dir = os.path.join(projects, "OrthoFinder")
    #     projects_results = helper._latest_output_dir(results_dir)
    #     args = ["-b",  projects_results, "--restart-of-blast"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(7)
    # def test_orthofinder_old_version_run(self, projects, capfd):
    #     args = ["-f", projects, "--old-version"]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(8)
    # def test_orthofinder_no_file_fix_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects, "--no-fix-files"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(9)
    # def test_orthofinder_blast_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects, "-S", "blast"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text



    # @pytest.mark.order(10)
    # def test_orthofinder_dna_run(self, dna_projects, capfd):
    #     args = ["-f", dna_projects, "-d", "--astral"]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text

    # @pytest.mark.skip(reason="Need to update the dataset")
    # @pytest.mark.order(11)
    # def test_orthofinder_dna_run(self, dna_projects, capfd):
    #     args = ["-f", dna_projects, "-d", "--astral"]
    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text



    # # # @pytest.mark.order(10)
    # # # def test_orthofinder_dna_assign_run(dna_projects, dna_assign, capfd):
    # # #     results_dir = os.path.join(dna_projects, "OrthoFinder")
    # # #     projects_results = helper._latest_output_dir(results_dir)
    # # #     args = ["--core", projects_results, "-d", "--assign", dna_assign, "-nk", 2]

    # # #     try:
    # # #         ret = main(args)
    # # #     except SystemExit as e:
    # # #         # success if code is 0 OR None
    # # #         code = e.code
    # # #         if isinstance(code, int):
    # # #             assert code == 0
    # # #         else:
    # # #             assert code is None  # sys.exit() or sys.exit(None)
    # # #     else:
    # # #         # success path where main returns normally
    # # #         assert ret is None
    # # #     out, err = capfd.readouterr()
    # # #     text = helper.clean_text(out + "\n" + err)  # strip ANSI, etc.
    # # #     err_msg = re.search(r"(?mi)^\s*ERROR:\s", text)
    # # #     assert "Traceback (most recent call last)" not in err, text
    # # #     assert not err_msg, text


    # @pytest.mark.order(12)
    # def test_orthofinder_oneway_search_run(self, projects, gene_tree_method, capfd):
    #     args = ["-f", projects, "-1", "-z", "-y"]

    #     if gene_tree_method != "fasttree":
    #         args += ["-T", gene_tree_method]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text

    # @pytest.mark.skip(reason="The test dataset is too small")
    # @pytest.mark.order(13)
    # def test_orthofinder_iqtree_run(self, projects, capfd):
    #     args = ["-f", projects, "-T", "iqtree3"]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.order(14)
    # def test_orthofinder_fasttree_fastest_run(self, projects, capfd):
    #     args = ["-f", projects, "-T", "fasttree_fastest"]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text


    # @pytest.mark.skip(reason="The test dataset is too small")
    # @pytest.mark.order(15)
    # def test_orthofinder_raxml_run(self, projects, capfd):
    #     args = ["-f", projects, "-T", "raxml"]

    #     code, out, err, text = helper._run_main(args, capfd)
    #     assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

    #     fatal_patterns = (
    #         r"(?mi)^\s*ERROR:\s",              
    #         r"(?m)^\s*Traceback \(most recent call last\):",    
    #         r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
    #         r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
    #     )
    #     assert not any(re.search(p, text) for p in fatal_patterns), text



    @pytest.mark.order(16)
    def test_orthofinder_famsa_species_tree_run(self, projects, gene_tree_method, species_tree, capfd):
        args = ["-f", projects, "-s", species_tree, "-rmgt", "-rmrgt", "-rmn0", "--save-space"]

        if gene_tree_method != "fasttree":
            args += ["-T", gene_tree_method]

        code, out, err, text = helper._run_main(args, capfd)
        assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

        fatal_patterns = (
            r"(?mi)^\s*ERROR:\s",              
            r"(?m)^\s*Traceback \(most recent call last\):",    
            r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
            r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
        )
        assert not any(re.search(p, text) for p in fatal_patterns), text





    @pytest.mark.order(17)
    def test_orthofinder_assign_with_species_tree_run(self, projects, assign, species_tree_assign, gene_tree_method, capfd):
        results_dir = os.path.join(projects, "OrthoFinder")
        projects_results = helper._latest_output_dir(results_dir, fileno=-1)
        args = ["--core", projects_results, "--assign",  assign, "-nk", 2, "-s", species_tree_assign, "-rmgt", "-rmrgt", "-rmn0", "--save-space"]

        if gene_tree_method != "fasttree":
            args += ["-T", gene_tree_method]

        code, out, err, text = helper._run_main(args, capfd)
        assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

        fatal_patterns = (
            r"(?mi)^\s*ERROR:\s",              
            r"(?m)^\s*Traceback \(most recent call last\):",    
            r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
            r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
        )
        assert not any(re.search(p, text) for p in fatal_patterns), text






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
