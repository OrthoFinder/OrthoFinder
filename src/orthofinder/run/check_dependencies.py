from ..utils import util, program_caller, file_handler, parallel_task_manager
from ..gene_tree import msa_trees
import os
import subprocess
from .. import __location__
from ..utils.printing import printer
from ..utils import printing

def can_run_mcl():
    command = "mcl -h"
    if parallel_task_manager.can_run_command(command):
        return True
    else:
        printer.print("ERROR: Cannot run MCL with the command \"%s\"" % command, style="error")
        program_caller.ProgramCaller.print_dependency_check_failure(command)
        printer.print("Please check MCL is installed and in the system path. See information above.\n", style="error")
        return False

def can_run_astral():
    d_deps_test = file_handler.handler.get_dependencies_check_dir()
    fn_test = d_deps_test + "astral.input.nwk"
    with open(fn_test, 'w') as outfile:
        for _ in range(20):
            outfile.write("(((A:1,B:1):2,(C:1,D:1):2,E:3):4);i\n")
    fn_output = d_deps_test + "astral.output.nwk"
    cmd = " ".join(["astral-pro", "-i", fn_test, "-o", fn_output, "-t", "2"])
    if parallel_task_manager.can_run_command(cmd, q_allow_stderr=True, q_require_stdout=False, q_check_return_code=True):
        return True
    else:
        printer.print("ERROR: Cannot run astral-pro", style="error")
        program_caller.ProgramCaller.print_dependency_check_failure(cmd)
        printer.print("Please check astral-pro is installed and that the executables are in the system path\n", style="error")
        return False

def check_dependencies(options, user_specified_m, prog_caller, dir_for_temp_files, q_assign=False):
    """q_assign: an --assign analysis (also when continued with -b)"""
    q_assign = q_assign or options.q_fast_add
    printing.print_underline("Checking required programs are installed", True)
    if not user_specified_m:
        printer.print('Running with the recommended MSA tree inference by default. To revert to legacy method use "-M dendroblast".\n')
    if options.q_start_from_fasta or q_assign:
        d_deps_check = file_handler.handler.get_dependencies_check_dir()
        success, stdout, stderr, cmd = prog_caller.test_search_method(options.search_program, 
                                                                    d_deps_check,
                                                                    scorematrix=options.score_matrix,
                                                                    gapopen=options.gapopen,
                                                                    gapextend=options.gapextend
                                                                    )
        if not success:
            printer.print("\nERROR: Cannot run %s" % options.search_program, style="error")
            prog_caller.print_dependency_check_failure(cmd)
            printer.print("Please check %s is installed and that the executables are in the system path\n" % options.search_program,  style="error")
            util.fail()
    if (options.q_start_from_fasta or options.q_start_from_blast) and not can_run_mcl():
        util.fail()
    if not (options.q_stop_after_prepare or options.q_stop_after_seqs or options.q_stop_after_groups or options.q_start_from_trees):
        if not can_run_ortholog_dependencies(dir_for_temp_files, 
                                            options.q_msa_trees, 
                                            options.q_phyldog, 
                                            options.q_stop_after_trees, 
                                            options.msa_program, 
                                            options.tree_program, 
                                            options.recon_method,
                                            prog_caller, 
                                            options.q_stop_after_alignments):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent orthologue inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.fail()
    if q_assign:
        if not can_run_astral():
            util.fail()

def write_test_distances_file(test_fn):
    with open(test_fn, 'w') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
    return test_fn

def can_run_ortholog_dependencies(
        working_dir, 
        q_msa_gene_trees, 
        q_phyldog, 
        q_stop_after_trees, 
        msa_method, 
        tree_method,
        recon_method, 
        program_caller, 
        q_stop_after_alignments
    ):
    d_deps_test = file_handler.handler.get_dependencies_check_dir()
    # FastME
    if not q_msa_gene_trees:
        test_fn = d_deps_test + "SimpleTest.phy"
        write_test_distances_file(test_fn)
        out_fn = d_deps_test + "SimpleTest.tre"
        if os.path.exists(out_fn): os.remove(out_fn)       
        cmd = "fastme -i %s -o %s" % (test_fn, out_fn)
        if not parallel_task_manager.can_run_command(cmd, q_allow_stderr=False):
            printer.print("ERROR: Cannot run fastme", style="error")
            program_caller.print_dependency_check_failure(cmd)
            printer.print("Please check FastME is installed and that the executables are in the system path.\n", style="error")
            return False
    # DLCPar
    # if ("dlcpar" in recon_method) and not (qStopAfterTrees or qStopAfterAlignments):
    #     if not parallel_task_manager.can_run_command("dlcpar_search --version", qAllowStderr=False):
    #         printer.print("ERROR: Cannot run dlcpar_search", style="error")
    #         printer.print("Please check DLCpar is installed and that the executables are in the system path.\n", style="error")
    #         return False
        # if recon_method == "dlcpar_convergedsearch":
        #     capture = subprocess.Popen(
        #         "dlcpar_search --version", 
        #         shell=True, 
        #         stdout=subprocess.PIPE, 
        #         stderr=subprocess.PIPE, 
        #         env=parallel_task_manager.my_env
        #     )
        #     stdout = [x for x in capture.stdout]
        #     try:
        #         stdout = "".join([x.decode() for x in stdout])
        #     except (UnicodeDecodeError, AttributeError):
        #         stdout = "".join([x.encode() for x in stdout])
        #     version = stdout.split()[-1]
        #     tokens = list(map(int, version.split(".")))
        #     major, minor = tokens[:2]
        #     release = tokens[2] if len(tokens) > 2 else 0
        #     # require 1.0.1 or above            
        #     actual = (major, minor, release)
        #     required = [1,0,1]
        #     versionOK = True
        #     for r, a in zip(required, actual):
        #         if a > r:
        #             versionOK = True
        #             break
        #         elif a < r:
        #             versionOK = False
        #             break
        #         else:
        #             pass
        #             # need to check next level down
        #     if not versionOK:
        #         printer.print("ERROR: dlcpar_convergedsearch requires dlcpar_search version 1.0.1 or above", style="error")
        #         return False                   
    
    # FastTree & MAFFT
    if q_msa_gene_trees or q_phyldog:
        test_fn = msa_trees.write_test_file(d_deps_test)
        if msa_method is not None:
            success, stdout, stderr, cmd = program_caller.test_msa_method(msa_method, d_deps_test)
            if not success and msa_method.startswith("mafft"):
                mafft_var = "MAFFT_BINARIES"
                print("Trying OrthoFinder packaged version of MAFFT\n")
                if mafft_var not in parallel_task_manager.my_env:
                    parallel_task_manager.my_env[mafft_var] = os.path.join(__location__, 'bin/mafft/libexec/')
                    parallel_task_manager.my_env["PATH"] = parallel_task_manager.my_env["PATH"] + ":" + \
                                                           os.path.join(__location__, 'bin/mafft/bin/')
                    success, stdout, stderr, cmd = program_caller.test_msa_method(msa_method, d_deps_test)
            if not success:
                printer.print("ERROR: Cannot run MSA method '%s'" % msa_method, style="error")
                program_caller.print_dependency_check_failure(cmd)
                printer.print("Please check program is installed. If it is user-configured please check the configuration in the "
                      "orthofinder/config.json file\n", style="error")
                return False
        if tree_method is not None:
            if q_msa_gene_trees and (not q_stop_after_alignments):
                success, stdout, stderr, cmd = program_caller.test_tree_method(tree_method, d_deps_test)
                if not success:
                   printer.print("ERROR: Cannot run tree method '%s'" % tree_method, style="error")
                   program_caller.print_dependency_check_failure(cmd)
                   printer.print("Please check program is installed. If it is user-configured please check the configuration in "
                         "the orthofinder/config.json file\n", style="error")
                   return False
            
    # if qPhyldog:
    #     if not parallel_task_manager.can_run_command("mpirun -np 1 phyldog", qAllowStderr=False):
    #         printer.print("ERROR: Cannot run mpirun -np 1 phyldog", style="error")
    #         printer.print("Please check phyldog is installed and that the executable is in the system path\n", style="error")
    #         return False
        
    return True    