import os
from ..utils import util, files
from orthofinder import __version__, g_mclInflation
from . import helpinfo
import shutil


# Control
class Options(object):#
    def __init__(self):
        self.nBlast = util.nThreadsDefault
        self.nProcessAlg = None
        self.qFastAdd = False  # Add species in near-linear time
        self.qStartFromBlast = False  # remove, just store BLAST to do
        self.qStartFromFasta = False  # local to argument checking
        self.qStartFromGroups = False
        self.qStartFromTrees = False
        self.qStopAfterPrepare = False
        self.qStopAfterGroups = False
        self.qStopAfterSeqs = False
        self.qStopAfterAlignments = False
        self.qStopAfterTrees = False
        self.qMSATrees = True  # Updated default
        self.qAddSpeciesToIDs = True
        self.qTrim = True
        self.gathering_version = (1, 0)    # < 3 is the original method
        self.search_program = "diamond"
        self.msa_program = "mafft"
        self.tree_program = "fasttree"
        self.recon_method = "of_recon"
        self.name = None   # name to identify this set of results
        self.qDoubleBlast = True
        self.qSplitParaClades = False
        self.qPhyldog = False
        self.speciesXMLInfoFN = None
        self.speciesTreeFN = None
        self.mclInflation = g_mclInflation
        self.dna = False
        self.fewer_open_files = True  # By default only open O(n) orthologs files at a time
        self.save_space = False  # On complete, have only one orthologs file per species
        self.v2_scores = False
        self.root_from_previous = False

    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)


def GetDirectoryArgument(arg, args):
    if len(args) == 0:
        print("Missing option for command line argument %s" % arg)
        util.Fail()
    directory = os.path.abspath(args.pop(0))
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
    return directory


def ProcessArgs(prog_caller, args):
    """ 
    Workflow
    | 1. Fasta Files | 2.  Prepare files    | 3.   Blast    | 4. Orthogroups    | 5.   Gene Trees     | 6.   Reconciliations/Orthologues   |

    Options
    Start from:
    -f: 1,2,..,6    (start from fasta files, --fasta)
    -b: 4,5,6       (start from blast results, --blast)
    -fg: 5,6         (start from orthogroups/do orthologue workflow, --from-groups)
    -ft: 6           (start from gene tree/do reconciliation, --from-trees)
    Stop at:
    -op: 2           (only prepare, --only-prepare)
    -og: 4           (orthogroups, --only-groups)
    """
    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.PrintHelp(prog_caller)
        util.Success() 

    options = Options()
    fastaDir = None
    continuationDir = None
    resultsDir_nonDefault = None
    pickleDir_nonDefault = None
    q_selected_msa_options = False
    q_selected_tree_options = False
    q_selected_search_option = False
    user_specified_M = False
    
    """
    -f: store fastaDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """    
    
    while len(args) > 0:
        arg = args.pop(0)

        if arg == "-f" or arg == "--fasta":
            if options.qStartFromFasta:
                print("Repeated argument: -f/--fasta\n")
                util.Fail()
            options.qStartFromFasta = True
            fastaDir = GetDirectoryArgument(arg, args)

        elif arg == "-b" or arg == "--blast":
            if options.qStartFromBlast:
                print("Repeated argument: -b/--blast\n")
                util.Fail()
            options.qStartFromBlast = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "--assign":
            options.qFastAdd = True
            fastaDir = GetDirectoryArgument(arg, args)

        elif arg == "--core":
            options.qFastAdd = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-fg" or arg == "--from-groups":
            if options.qStartFromGroups:
                print("Repeated argument: -fg/--from-groups\n")
                util.Fail()
            options.qStartFromGroups = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-ft" or arg == "--from-trees":
            if options.qStartFromTrees:
                print("Repeated argument: -ft/--from-trees\n")
                util.Fail()
            options.qStartFromTrees = True
            continuationDir = GetDirectoryArgument(arg, args)

        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nBlast = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()

        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
            except:
                print("Incorrect argument for number of BLAST threads: %s\n" % arg)
                util.Fail()  

        elif arg == "-1":
            options.qDoubleBlast = False

        elif arg == "-d" or arg == "--dna":
            options.dna = True
            if not q_selected_search_option:
                options.search_program = "blast_nucl"

        elif arg == "-X":
            options.qAddSpeciesToIDs = False

        elif arg == "-y":
            options.qSplitParaClades = True

        elif arg == "-z":
            options.qTrim = False

        elif arg == "-I" or arg == "--inflation":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.mclInflation = float(arg)
            except:
                print("Incorrect argument for MCL inflation parameter: %s\n" % arg)
                util.Fail()

        elif arg == "-c1":
            print("\nThe option 'c1' has been renamed '--c-homologs'")
            util.Fail()

        elif arg == "--c-homologs":
            options.gathering_version = (3, 2)

        elif arg == "--save-space":
            options.save_space = True

        elif arg == "-x" or arg == "--orthoxml":  
            if options.speciesXMLInfoFN:
                print("Repeated argument: -x/--orthoxml")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.speciesXMLInfoFN = args.pop(0)

        elif arg == "-n" or arg == "--name":  
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.name = args.pop(0)
            while options.name.endswith("/"): options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]): 
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()

        elif arg == "-o" or arg == "--output":  
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)

            while resultsDir_nonDefault.endswith("/"): 
                resultsDir_nonDefault = resultsDir_nonDefault[:-1]

            resultsDir_nonDefault += "/"
            if os.path.exists(resultsDir_nonDefault):
                print("ERROR: non-default output directory already exists: %s\n" % resultsDir_nonDefault)
                util.Fail()

            if " " in resultsDir_nonDefault:
                print("ERROR: non-default output directory cannot include spaces: %s\n" % resultsDir_nonDefault)
                util.Fail()
            checkDirName = resultsDir_nonDefault

            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]

            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print("ERROR: location '%s' for results directory '%s' does not exist.\n" % (path, newDir))
                util.Fail()

        elif arg == "-s" or arg == "--speciestree":  
            if options.speciesTreeFN:
                print("Repeated argument: -s/--speciestree")
                util.Fail()

            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.speciesTreeFN = args.pop(0)

        elif arg == "--scores-v2":
            options.v2_scores = True

        elif arg == "-S" or arg == "--search":
            choices = ['blast'] + prog_caller.ListSearchMethods()

            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0)
            if arg in choices:
                options.search_program = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-M" or arg == "--method":
            arg_M_or_msa = arg
            user_specified_M = True
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0)
            if arg == "msa": 
                options.qMSATrees = True
            elif arg == "phyldog": 
                options.qPhyldog = True
                options.recon_method = "phyldog"
                options.qMSATrees = False
            elif arg == "dendroblast":
                options.qMSATrees = False
            else:
                print("Invalid argument for option %s: %s" % (arg_M_or_msa, arg))
                print("Valid options are 'dendroblast' and 'msa'\n")
                util.Fail()

        elif arg == "-A" or arg == "--msa_program":
            choices = ['mafft'] + prog_caller.ListMSAMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.msa_program = arg
                q_selected_msa_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-T" or arg == "--tree_program":
            choices = ['fasttree'] + prog_caller.ListTreeMethods()
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.tree_program = arg
                q_selected_tree_options = True
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-R" or arg == "--recon_method":
            choices = ['of_recon', 'dlcpar', 'dlcpar_convergedsearch', 'only_overlap']
            switch_used = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            if arg in choices:
                options.recon_method = arg
            else:
                print("Invalid argument for option %s: %s" % (switch_used, arg))
                print("Valid options are: {%s}\n" % (", ".join(choices)))
                util.Fail()

        elif arg == "-p":
            pickleDir_nonDefault = GetDirectoryArgument(arg, args)

        elif arg == "-op" or arg == "--only-prepare":
            options.qStopAfterPrepare = True

        elif arg == "-og" or arg == "--only-groups":
            options.qStopAfterGroups = True

        elif arg == "-os" or arg == "--only-seqs":
            options.qStopAfterSeqs = True

        elif arg == "-oa" or arg == "--only-alignments":
            options.qStopAfterAlignments = True

        elif arg == "-ot" or arg == "--only-trees":
            options.qStopAfterTrees = True

        elif arg == "-h" or arg == "--help":
            helpinfo.PrintHelp(prog_caller)
            util.Success()

        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()

    # set a default for number of algorithm threads
    if options.nProcessAlg is None:
        options.nProcessAlg = min(16, max(1, int(options.nBlast/8)))

    if options.nBlast < 1:
        print("ERROR: Number of '-t' threads cannot be fewer than 1, got %d" % options.nBlast)
        util.Fail()    

    if options.nProcessAlg < 1:
        print("ERROR: Number of '-a' threads cannot be fewer than 1, got %d" % options.nProcessAlg)
        util.Fail()  

    # check argument combinations       
    if not (options.qStartFromFasta or options.qStartFromBlast or options.qStartFromGroups or options.qStartFromTrees or options.qFastAdd):
        print("ERROR: Please specify the input directory for OrthoFinder using one of the options: '-f', '-b', '-fg' or '-ft', '--assign'.")
        util.Fail()

    if options.qFastAdd:
        if (options.qStartFromFasta or options.qStartFromBlast or options.qStartFromGroups or options.qStartFromTrees):
            print("ERROR: Incompatible options used with --assign, cannot accept: '-f', '-b', '-fg' or '-ft'")
            util.Fail()
        if fastaDir is None:
            print("ERROR: '--assign' option also requires '--core' directory to be specified")
            util.Fail()
        if continuationDir is None:
            print("ERROR: '--core' option also requires '--assign' directory to be specified")
            util.Fail()
        if not options.qMSATrees:
            print("ERROR: --assign requires MSA trees, option '-M dendroblast' is invalid")
            util.Fail()

    if options.qStartFromFasta and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -f (start from fasta files) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()
        
    if options.qStartFromBlast and (options.qStartFromTrees or options.qStartFromGroups):
        print("ERROR: Incompatible arguments, -b (start from pre-calcualted BLAST results) and" + (" -fg (start from orthogroups)" if options.qStartFromGroups else " -ft (start from trees)"))
        util.Fail()      

    if options.qStartFromTrees and options.qStartFromGroups:
        print("ERROR: Incompatible arguments, -fg (start from orthogroups) and -ft (start from trees)")
        util.Fail()    

    if options.qStopAfterSeqs and (not options.qMSATrees):
        print("ERROR: Argument '-os' (stop after sequences) also requires option '-M msa'")
        util.Fail()   

    if options.qStopAfterAlignments and (not options.qMSATrees):
        print("ERROR: Argument '-oa' (stop after alignments) also requires option '-M msa'")
        util.Fail()     

    if (q_selected_msa_options or q_selected_tree_options) and (not options.qMSATrees and not options.qPhyldog):
        print("ERROR: Argument '-A' or '-T' (multiple sequence alignment/tree inference program) also requires option '-M msa'")
        util.Fail()       
        
    if options.qPhyldog and (not options.speciesTreeFN):
        print("ERROR: Phyldog currently needs a species tree to be provided")
        util.Fail()          

    if resultsDir_nonDefault != None and ((not options.qStartFromFasta) or options.qStartFromBlast):
        print("ERROR: Incompatible arguments, -o (non-default output directory) can only be used with a new OrthoFinder run using option '-f'")
        util.Fail()       
        
    if options.search_program not in (prog_caller.ListSearchMethods() + ['blast']):
        print("ERROR: Search program (%s) not configured in config.json file" % options.search_program)
        util.Fail()

    util.PrintTime("Starting OrthoFinder %s" % __version__)    
    print("%d thread(s) for highly parallel tasks (BLAST searches etc.)" % options.nBlast)
    print("%d thread(s) for OrthoFinder algorithm\n" % options.nProcessAlg)

    if options.qFastAdd and not q_selected_msa_options:
        print("INFO: For --assign defaulting to 'mafft --memsave' to reduce RAM usage\n")
        options.msa_program = "mafft_memsave"

    if options.qFastAdd and not q_selected_tree_options:
        print("INFO: For --assign defaulting to 'FastTree -fastest' to reduce RAM usage\n")
        options.tree_program = "fasttree_fastest"

    return options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault, user_specified_M


def DeleteDirectoryTree(d):
    if os.path.exists(d): 
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)   


def CheckOptions(options, speciesToUse):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = list(SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN()).values())
        orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any((options.qStopAfterPrepare, options.qStopAfterGroups, options.qStopAfterSeqs, options.qStopAfterAlignments, options.qStopAfterTrees))
    if q_do_orthologs:
        n_sp = len(speciesToUse)
        wd = files.FileHandler.GetWorkingDirectory_Write()
        wd_files_test = wd + "Files_test/"
        fh = []
        try:
            if not os.path.exists(wd_files_test):
                os.mkdir(wd_files_test)
            for i_sp in range(n_sp):
                di = wd_files_test + "Sp%d/" % i_sp
                if not os.path.exists(di):
                    os.mkdir(di)
                for j_sp in range(1):  # We only create a linear number of ortholog files now
                    fnij = di + "Sp%d.txt" % j_sp
                    fh.append(open(fnij, 'w'))
            # create a few extra files to be safe
            for i_extra in range(n_extra):
                fh.append(open(wd_files_test + "Extra%d.txt" % i_extra, 'w'))
            # close the files again and delete
            for fhh in fh:
                fhh.close()
            DeleteDirectoryTree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                util.number_open_files_exception_advice(len(speciesToUse), False)
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                util.Fail()
            else:
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                print("ERROR: Attempted to open required files for OrthoFinder run but an unexpected error occurred. \n\nStacktrace:")
                raise
    return options