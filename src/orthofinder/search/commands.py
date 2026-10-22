from ..utils import util, file_handler
import itertools
import os 

def get_ordered_search_commands(
        seqs_info, 
        species_info_obj, 
        options, 
        prog_caller, 
        n_genes_per_species=None, 
        q_new_species_unassigned_genes=False
    ):

    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    n_genes_per_species: List[int] - number of genes per species
    """
    i_species_previous = list(range(species_info_obj.i_first_new_species))
    i_species_new = list(range(species_info_obj.i_first_new_species, species_info_obj.n_sp_all))
    if q_new_species_unassigned_genes:
        if n_genes_per_species is not None:
            exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
            i_species_new = list(set(i_species_new).difference(exclude))
        species_pairs = [(i, j) for i, j in itertools.product(i_species_new, i_species_new)]
        task_sizes = [n_genes_per_species[i]*n_genes_per_species[j] for i,j in species_pairs]
    else:
        species_pairs = [(i, j) for i, j in itertools.product(i_species_new, i_species_new) if (options.q_double_blast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(i_species_new, i_species_previous) if (options.q_double_blast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(i_species_previous, i_species_new) if (options.q_double_blast or i <=j)]
        task_sizes = [seqs_info.n_seqs_per_species[i]*seqs_info.n_seqs_per_species[j] for i,j in species_pairs]
    
    # Smaller files should be processed first 
    q_largest_first = False
    if options.cmd_order != "ascending":
        q_largest_first = True 
    
    task_sizes, species_pairs = util.sort_array_pair_by_first(task_sizes, species_pairs, q_largest_first=q_largest_first)
    
    method_threads = None

    commands = [prog_caller.get_search_method_command_search(options.search_program,
                    file_handler.handler.get_species_unassigned_fasta_fn(i_fasta) if q_new_species_unassigned_genes else file_handler.handler.get_species_fasta_fn(i_fasta),
                    file_handler.handler.get_species_database_n(i_db, options.search_program),
                    file_handler.handler.get_blast_results_fn(i_fasta, i_db, q_for_creation=True),
                    scorematrix=options.score_matrix,
                    gapopen=options.gapopen,
                    gapextend=options.gapextend,
                    method_threads=method_threads,
                    nucleotide=options.dna,
                    ) 
                    for i_fasta, i_db in species_pairs]
    return commands, task_sizes


def get_ordered_search_commands_clades(
        seqs_info, 
        species_info_obj, 
        options, prog_caller,
        n_genes_per_species, 
        species_clades
    ):
    """
    Search all species
    """
    method_threads = None

    exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
    species_pairs = []
    for clade in species_clades:
        clade = list(set(clade).difference(exclude))
        species_pairs.extend([(i, j) for i, j in itertools.product(clade, clade)])
    commands = []
    for i_fasta, i_db in species_pairs:
        if options.search_program in ["mmseqs"]:
            unassigned = file_handler.handler.get_species_unassigned_fasta_fn(i_fasta)
            query = os.path.join(os.path.dirname(unassigned), os.path.basename(unassigned).split(".", 1)[-1])
        else:
            query = file_handler.handler.get_species_unassigned_fasta_fn(i_fasta)

        commands.append(
            prog_caller.get_search_method_command_search(
                options.search_program,
                query,
                file_handler.handler.get_species_database_n(i_db, options.search_program),
                file_handler.handler.get_blast_results_fn(i_fasta, i_db, q_for_creation=True),
                scorematrix=options.score_matrix,
                gapopen=options.gapopen,
                gapextend=options.gapextend,
                method_threads=method_threads,
                nucleotide=options.dna,
            )
        )
    # commands = [
    #     prog_caller.get_search_method_command_search(
    #         options.search_program,
    #         files.FileHandler.get_species_unassigned_fasta_fn(iFasta),
    #         files.FileHandler.get_species_database_n(iDB, options.search_program),
    #         files.FileHandler.get_blast_results_fn(iFasta, iDB, qForCreation=True),
    #         scorematrix=options.score_matrix,
    #         gapopen=options.gapopen,
    #         gapextend=options.gapextend,
    #         method_threads=method_threads
    #     ) 
    #     for iFasta, iDB in speciesPairs
    # ]
    tasksize =  [
        os.stat(file_handler.handler.get_species_database_n(i_db, options.search_program)).st_size
        for i_fasta, i_db in species_pairs
        if os.path.isfile(file_handler.handler.get_species_database_n(i_db, options.search_program)) \
            and os.path.exists(file_handler.handler.get_species_database_n(i_db, options.search_program))
    ]

    return  commands, tasksize


def get_ordered_search_commands_assign(species_info_obj, diamond_db, options, prog_caller, q_one_query, threads=1):
    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing
    the BLAST commands.
    """
    i_species_new = list(range(species_info_obj.i_first_new_species, species_info_obj.n_sp_all))
    tasksize = None

    method_threads = None
    
    if q_one_query:
        wd = file_handler.handler.get_species_seqs_dir()[0]
        fn_single_fasta = "%sall_sequences.fa" % wd
        with open(fn_single_fasta, 'w', encoding="utf-8", newline="\n") as outfile:
            for i_fasta in i_species_new:
                with open(file_handler.handler.get_species_fasta_fn(i_fasta), 'r', encoding="utf-8") as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write("\n")
        results = wd + "Blast_all_sequences.txt"
        # commands = [prog_caller.get_search_method_command_search(search_program, fn_single_fasta, diamond_db, results,
                                                                # scorematrix=options.score_matrix, 
                                                                # gapopen=options.gapopen, 
                                                                # gapextend=options.gapextend,#
                                                                # method_threads=options.method_threads)]
        commands = ["diamond blastp --ignore-warnings -d %s -q %s -o %s --more-sensitive -p %d --quiet -e 0.001 --compress 1" % (diamond_db, fn_single_fasta, results, threads)]
        results_files = [results + ".gz"]
    else:
        commands = [
            prog_caller.get_search_method_command_search(
                options.search_program,
                file_handler.handler.get_species_fasta_fn(i_fasta),
                diamond_db,
                file_handler.handler.get_blast_results_fn(i_fasta, -1, q_for_creation=True),
                scorematrix=options.score_matrix, 
                gapopen=options.gapopen, 
                gapextend=options.gapextend,#
                method_threads=method_threads,
                nucleotide=options.dna,
            )
            for i_fasta in i_species_new
        ]
        if options.search_program in ["mmseqs"]:
            results_files = [file_handler.handler.get_blast_results_fn(i_fasta, -1, q_for_creation=True) for i_fasta in i_species_new]
        else:
            results_files = [file_handler.handler.get_blast_results_fn(i_fasta, -1, q_for_creation=True) + ".gz" for i_fasta in i_species_new]
        tasksize = [
            os.stat(file_handler.handler.get_species_fasta_fn(i_fasta)).st_size
            for i_fasta in i_species_new
            if os.path.isfile(file_handler.handler.get_species_fasta_fn(i_fasta)) \
                and os.path.exists(file_handler.handler.get_species_fasta_fn(i_fasta))
        ]
    return commands, tasksize, results_files