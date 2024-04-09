from ..utils import util, program_caller, files, parallel_task_manager
import subprocess  
import glob     
import shutil
from . import run_info

def RunBlastDBCommand(command):
    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env, shell=True)
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stdout_lines = stdout.count("\n")
    n_stderr_lines = stderr.count("\n")
    nLines_success= 12
    if n_stdout_lines > nLines_success or n_stderr_lines > 0 or capture.returncode != 0:
        print("\nWARNING: Likely problem with input FASTA files")
        if capture.returncode != 0:
            print("makeblastdb returned an error code: %d" % capture.returncode)
        else:
            print("makeblastdb produced unexpected output")
        print("Command: %s" % " ".join(command))
        print("stdout:\n-------")
        print(stdout)
        if len(stderr) > 0:
            print("stderr:\n-------")
            print(stderr)


# 7
def RunSearch(options, speciessInfoObj, seqsInfo, prog_caller, 
              q_new_species_unassigned_genes=False, 
              n_genes_per_species=None, species_clades=None):
    """
    n_genes_per_species: List[int] - optional, for use with unassigned genes. If a species has zero unassigned genes, don't search with/agaisnt it.
    """
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    if options.qStopAfterPrepare:
        util.PrintUnderline("%s commands that must be run" % name_to_print)
    elif species_clades is not None:
        util.PrintUnderline("Running %s searches for new non-core species clades" % name_to_print)
    else:
        util.PrintUnderline("Running %s all-versus-all" % name_to_print)
    if species_clades is None:
        commands = run_info.GetOrderedSearchCommands(seqsInfo, 
                                                     speciessInfoObj,
                                                     options,
                                                     prog_caller, 
                                                     n_genes_per_species, 
                                                     q_new_species_unassigned_genes=q_new_species_unassigned_genes)
    else:
        commands = run_info.GetOrderedSearchCommands_clades(seqsInfo, 
                                                            speciessInfoObj,
                                                            options,
                                                            prog_caller, 
                                                            n_genes_per_species, 
                                                            species_clades)
    if options.qStopAfterPrepare:
        for command in commands:
            print(command)
        util.Success()
    print("Using %d thread(s)" % options.nBlast)
    util.PrintTime("This may take some time....")
    program_caller.RunParallelCommands(options.nBlast, commands, qListOfList=False,
                                       q_print_on_error=True, q_always_print_stderr=False)

    # remove BLAST databases
    util.PrintTime("Done all-versus-all sequence search")
    if options.search_program == "blast":
        for f in glob.glob(files.FileHandler.GetWorkingDirectory1_Read()[0] + "BlastDBSpecies*"):
            os.remove(f)
    if options.search_program == "mmseqs":
        for i in range(speciessInfoObj.nSpAll):
            for j in range(speciessInfoObj.nSpAll):
                tmp_dir = "/tmp/tmpBlast%d_%d.txt" % (i,j)
                if os.path.exists(tmp_dir):
                    try:
                        shutil.rmtree(tmp_dir)
                    except OSError:
                        time.sleep(1)
                        shutil.rmtree(tmp_dir, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted


# 6
def CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=False):
    iSpeciesToDo = range(max(speciesInfoObj.speciesToUse) + 1)
    for iSp in iSpeciesToDo:
        fn_fasta = files.FileHandler.GetSpeciesUnassignedFastaFN(iSp) if q_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iSp)
        if options.search_program == "blast":
            command = " ".join(["makeblastdb", "-dbtype", "prot", "-in", fn_fasta, "-out", files.FileHandler.GetSpeciesDatabaseN(iSp)])
            util.PrintTime("Creating Blast database %d of %d" % (iSp + 1, len(iSpeciesToDo)))
            run_commands.RunBlastDBCommand(command) 
        else:
            command = prog_caller.GetSearchMethodCommand_DB(options.search_program, 
                                                            fn_fasta, 
                                                            files.FileHandler.GetSpeciesDatabaseN(iSp, options.search_program),
                                                            options.score_matrix,
                                                            options.gapopen,
                                                            options.gapextend)
            util.PrintTime("Creating %s database %d of %d" % (options.search_program, iSp + 1, len(iSpeciesToDo)))
            ret_code = parallel_task_manager.RunCommand(command, qPrintOnError=True, qPrintStderr=False)
            if ret_code != 0:
                files.FileHandler.LogFailAndExit("ERROR: diamond makedb failed")


def RunSearch_accelerate(options, speciessInfoObj, fn_diamond_db, prog_caller, q_one_query=False):
    name_to_print = "BLAST" if options.search_program == "blast" else options.search_program
    util.PrintUnderline("Running %s profiles search" % name_to_print)
    commands, results_files = run_info.GetOrderedSearchCommands_accelerate(speciessInfoObj, 
                                                       fn_diamond_db, 
                                                       options,
                                                       prog_caller, 
                                                       q_one_query=q_one_query, 
                                                       threads=options.nBlast)
    if q_one_query:
        return_code = parallel_task_manager.RunCommand(commands[0], qPrintOnError=True)
        if return_code != 0:
            print("ERROR: DIAMOND search failed, see messages above")
            util.Fail()
        util.PrintTime("Done profiles search\n")
        return results_files
    program_caller.RunParallelCommands(options.nBlast, commands, qListOfList=False, q_print_on_error=True)
    util.PrintTime("Done profiles search")
    return results_files