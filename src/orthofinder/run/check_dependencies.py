from ..utils import util, program_caller, files, parallel_task_manager
from ..tools import trees_msa
import os
import subprocess
from .. import __location__, my_env

def CanRunBLAST():
    if parallel_task_manager.CanRunCommand("makeblastdb -help") and parallel_task_manager.CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        program_caller.ProgramCaller.PrintDependencyCheckFailure("makeblastdb -help\nblastp -help")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if parallel_task_manager.CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        program_caller.ProgramCaller.PrintDependencyCheckFailure(command)
        print("Please check MCL is installed and in the system path. See information above.\n")
        return False

def CanRunASTRAL():
    d_deps_test = files.FileHandler.GetDependenciesCheckDir()
    fn_test = d_deps_test + "astral.input.nwk"
    with open(fn_test, 'w') as outfile:
        for _ in range(20):
            outfile.write("(((A:1,B:1):2,(C:1,D:1):2,E:3):4);i\n")
    fn_output = d_deps_test + "astral.output.nwk"
    cmd = " ".join(["astral-pro", "-i", fn_test, "-o", fn_output, "-t", "2"])
    if parallel_task_manager.CanRunCommand(cmd, qAllowStderr=True, qRequireStdout=False, qCheckReturnCode=True):
        return True
    else:
        print("ERROR: Cannot run astral-pro")
        program_caller.ProgramCaller.PrintDependencyCheckFailure(cmd)
        print("Please check astral-pro is installed and that the executables are in the system path\n")
        return False

def CheckDependencies(options, user_specified_m, prog_caller, dirForTempFiles):
    util.PrintUnderline("Checking required programs are installed")
    if not user_specified_m:
        print("Running with the recommended MSA tree inference by default. To revert to legacy method use '-M dendroblast'.\n")
    if options.qStartFromFasta or options.qFastAdd:
        if options.search_program == "blast":
            if not CanRunBLAST(): util.Fail()
        else:
            d_deps_check = files.FileHandler.GetDependenciesCheckDir()
            success, stdout, stderr, cmd = prog_caller.TestSearchMethod(options.search_program, 
                                                                        d_deps_check,
                                                                        scorematrix=options.score_matrix,
                                                                        gapopen=options.gapopen,
                                                                        gapextend=options.gapextend
                                                                        )
            if not success:
                print("\nERROR: Cannot run %s" % options.search_program)
                prog_caller.PrintDependencyCheckFailure(cmd)
                print("Please check %s is installed and that the executables are in the system path\n" % options.search_program)
                util.Fail()
    if (options.qStartFromFasta or options.qStartFromBlast) and not CanRunMCL():
        util.Fail()
    if not (options.qStopAfterPrepare or options.qStopAfterSeqs or options.qStopAfterGroups or options.qStartFromTrees):
        if not CanRunOrthologueDependencies(dirForTempFiles, 
                                            options.qMSATrees, 
                                            options.qPhyldog, 
                                            options.qStopAfterTrees, 
                                            options.msa_program, 
                                            options.tree_program, 
                                            options.recon_method,
                                            prog_caller, 
                                            options.qStopAfterAlignments):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent orthologue inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.Fail()
    if options.qFastAdd:
        if not CanRunASTRAL():
            util.Fail()

def WriteTestDistancesFile(testFN):
    with open(testFN, 'w') as outfile:
        outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
    return testFN

def CanRunOrthologueDependencies(
        workingDir, 
        qMSAGeneTrees, 
        qPhyldog, 
        qStopAfterTrees, 
        msa_method, 
        tree_method,
        recon_method, 
        program_caller, 
        qStopAfterAlignments
    ):
    d_deps_test = files.FileHandler.GetDependenciesCheckDir()
    # FastME
    if not qMSAGeneTrees:
        testFN = d_deps_test + "SimpleTest.phy"
        WriteTestDistancesFile(testFN)
        outFN = d_deps_test + "SimpleTest.tre"
        if os.path.exists(outFN): os.remove(outFN)       
        cmd = "fastme -i %s -o %s" % (testFN, outFN)
        if not parallel_task_manager.CanRunCommand(cmd, qAllowStderr=False):
            print("ERROR: Cannot run fastme")
            program_caller.PrintDependencyCheckFailure(cmd)
            print("Please check FastME is installed and that the executables are in the system path.\n")
            return False
    # DLCPar
    if ("dlcpar" in recon_method) and not (qStopAfterTrees or qStopAfterAlignments):
        if not parallel_task_manager.CanRunCommand("dlcpar_search --version", qAllowStderr=False):
            print("ERROR: Cannot run dlcpar_search")
            print("Please check DLCpar is installed and that the executables are in the system path.\n")
            return False
        if recon_method == "dlcpar_convergedsearch":
            capture = subprocess.Popen("dlcpar_search --version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
            stdout = [x for x in capture.stdout]
            try:
                stdout = "".join([x.decode() for x in stdout])
            except (UnicodeDecodeError, AttributeError):
                stdout = "".join([x.encode() for x in stdout])
            version = stdout.split()[-1]
            tokens = list(map(int, version.split(".")))
            major, minor = tokens[:2]
            release = tokens[2] if len(tokens) > 2 else 0
            # require 1.0.1 or above            
            actual = (major, minor, release)
            required = [1,0,1]
            versionOK = True
            for r, a in zip(required, actual):
                if a > r:
                    versionOK = True
                    break
                elif a < r:
                    versionOK = False
                    break
                else:
                    pass
                    # need to check next level down
            if not versionOK:
                print("ERROR: dlcpar_convergedsearch requires dlcpar_search version 1.0.1 or above")
                return False                   
    
    # FastTree & MAFFT
    if qMSAGeneTrees or qPhyldog:
        testFN = trees_msa.WriteTestFile(d_deps_test)
        if msa_method is not None:
            success, stdout, stderr, cmd = program_caller.TestMSAMethod(msa_method, d_deps_test)
            if not success and msa_method.startswith("mafft"):
                mafft_var = "MAFFT_BINARIES"
                print("Trying OrthoFinder packaged version of MAFFT\n")
                if mafft_var not in parallel_task_manager.my_env:
                    parallel_task_manager.my_env[mafft_var] = os.path.join(__location__, 'bin/mafft/libexec/')
                    parallel_task_manager.my_env["PATH"] = parallel_task_manager.my_env["PATH"] + ":" + \
                                                           os.path.join(__location__, 'bin/mafft/bin/')
                    success, stdout, stderr, cmd = program_caller.TestMSAMethod(msa_method, d_deps_test)
            if not success:
                print("ERROR: Cannot run MSA method '%s'" % msa_method)
                program_caller.PrintDependencyCheckFailure(cmd)
                print("Please check program is installed. If it is user-configured please check the configuration in the "
                      "orthofinder/config.json file\n")
                return False
        if tree_method is not None:
            if qMSAGeneTrees and (not qStopAfterAlignments):
                success, stdout, stderr, cmd = program_caller.TestTreeMethod(tree_method, d_deps_test)
                if not success:
                   print("ERROR: Cannot run tree method '%s'" % tree_method)
                   program_caller.PrintDependencyCheckFailure(cmd)
                   print("Please check program is installed. If it is user-configured please check the configuration in "
                         "the orthofinder/config.json file\n")
                   return False
            
    if qPhyldog:
        if not parallel_task_manager.CanRunCommand("mpirun -np 1 phyldog", qAllowStderr=False):
            print("ERROR: Cannot run mpirun -np 1 phyldog")
            print("Please check phyldog is installed and that the executable is in the system path\n")
            return False
        
    return True    