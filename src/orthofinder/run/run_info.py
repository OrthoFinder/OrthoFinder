from ..utils import util, files
import itertools
import os 

def GetOrderedSearchCommands(seqsInfo, speciesInfoObj, 
                             options, prog_caller, 
                             n_genes_per_species=None, 
                             q_new_species_unassigned_genes=False):

    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    n_genes_per_species: List[int] - number of genes per species
    """
    iSpeciesPrevious = list(range(speciesInfoObj.iFirstNewSpecies))
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    if q_new_species_unassigned_genes:
        if n_genes_per_species is not None:
            exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
            iSpeciesNew = list(set(iSpeciesNew).difference(exclude))
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew)]
        taskSizes = [n_genes_per_species[i]*n_genes_per_species[j] for i,j in speciesPairs]
    else:
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew) if (options.qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious) if (options.qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew) if (options.qDoubleBlast or i <=j)]
        taskSizes = [seqsInfo.nSeqsPerSpecies[i]*seqsInfo.nSeqsPerSpecies[j] for i,j in speciesPairs]
    
    # Smaller files should be processed first 
    qLargestFirst = False
    if options.cmd_order != "ascending":
        qLargestFirst = True 
    
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, qLargestFirst=qLargestFirst)
    
    if options.old_version:
        method_threads = options.method_threads
    else:
        method_threads = None

    if options.search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(options.search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True),
                        scorematrix=options.score_matrix,
                        gapopen=options.gapopen,
                        gapextend=options.gapextend,
                        method_threads=method_threads
                        ) 
                        for iFasta, iDB in speciesPairs]
    return commands, taskSizes


def GetOrderedSearchCommands_clades(seqsInfo, 
                                    speciesInfoObj, 
                                    options, prog_caller,
                                    n_genes_per_species, 
                                    species_clades):
    """
    Search all species
    """
    if options.old_version:
        method_threads = options.method_threads
    else:
        method_threads = None

    exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
    speciesPairs = []
    for clade in species_clades:
        clade = list(set(clade).difference(exclude))
        speciesPairs.extend([(i, j) for i, j in itertools.product(clade, clade)])
    if options.search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) 
                    for iFasta, iDB in speciesPairs
                    ]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(options.search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True),
                        scorematrix=options.score_matrix,
                        gapopen=options.gapopen,
                        gapextend=options.gapextend,
                        method_threads=method_threads
                        ) 
                    for iFasta, iDB in speciesPairs
                    ]
    tasksize =  [
                os.stat(files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program)).st_size
                for iFasta, iDB in speciesPairs
                if os.path.isfile(files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program)) \
                    and os.path.exists(files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program))
                ]

    return  commands, tasksize


def GetOrderedSearchCommands_accelerate(speciesInfoObj, diamond_db, options, prog_caller, q_one_query, threads=1):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing
    the BLAST commands.
    """
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    tasksize = None

    if options.old_version:
        method_threads = options.method_threads
    else:
        method_threads = None
    
    if q_one_query:
        wd = files.FileHandler.GetSpeciesSeqsDir()[0]
        fn_single_fasta = "%sall_sequences.fa" % wd
        with open(fn_single_fasta, 'w') as outfile:
            for iFasta in iSpeciesNew:
                with open(files.FileHandler.GetSpeciesFastaFN(iFasta), 'r') as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write("\n")
        results = wd + "Blast_all_sequences.txt"
        # commands = [prog_caller.GetSearchMethodCommand_Search(search_program, fn_single_fasta, diamond_db, results,
                                                                # scorematrix=options.score_matrix, 
                                                                # gapopen=options.gapopen, 
                                                                # gapextend=options.gapextend,#
                                                                # method_threads=options.method_threads)]
        commands = ["diamond blastp --ignore-warnings -d %s -q %s -o %s --more-sensitive -p %d --quiet -e 0.001 --compress 1" % (diamond_db, fn_single_fasta, results, threads)]
        results_files = [results + ".gz"]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(
                        options.search_program,
                        files.FileHandler.GetSpeciesFastaFN(iFasta),
                        diamond_db,
                        files.FileHandler.GetBlastResultsFN(iFasta, -1, qForCreation=True),
                        scorematrix=options.score_matrix, 
                        gapopen=options.gapopen, 
                        gapextend=options.gapextend,#
                        method_threads=method_threads
            )
            for iFasta in iSpeciesNew
        ]
        results_files = [files.FileHandler.GetBlastResultsFN(iFasta, -1, qForCreation=True) + ".gz" for iFasta in iSpeciesNew]
        tasksize = [
            os.stat(files.FileHandler.GetSpeciesFastaFN(iFasta)).st_size
            for iFasta in iSpeciesNew
            if os.path.isfile(files.FileHandler.GetSpeciesFastaFN(iFasta)) \
                and os.path.exists(files.FileHandler.GetSpeciesFastaFN(iFasta))
        ]
    return commands, tasksize, results_files