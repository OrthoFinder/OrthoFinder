try:
    from rich import print
except ImportError:
    ...

from ..utils import logging as run_logging
from ..utils import util, file_handler
from ..phylo import tree
from . import dendroblast
from . import msa_trees
from ..utils import printing
from ..phylo import taxa

def convert_user_species_tree(species_tree_fn_in, species_dict, species_tree_fn_out):
    t = tree.Tree(species_tree_fn_in, format=1)  
    t.prune(t.get_leaf_names())
    rev_dict = {v:k for k,v in species_dict.items()}
    for sp in t:
        sp.name = rev_dict[sp.name]       
    t.write(outfile=species_tree_fn_out)

@run_logging.RunLogger.stage("Build alignments and gene trees", require_result=True)
def infer_gene_and_species_trees(
        og_set,
        program_caller,
        msa_method,
        tree_method,
        n_high_parallel,
        n_low_parallel,
        q_double_blast,
        q_add_species_to_ids,
        q_trim,
        cmd_order="descending",
        method_threads=None,
        method_threads_large=None,
        method_threads_small=None, 
        threshold=None,
        user_species_tree = None,
        q_stop_after_seqs = False,
        q_stop_after_align = False,
        q_msa = False,
        q_phyldog = False,
        results_name = "",
        root_from_previous = False,
        i_og_restart=0,
        fix_files=True,
        astral=True,
        dynamic_threads=False,
        n_skip=50,
        msa_min_seq=10
    ):
    """
    1. Setup:
        - ogSet, directories
        - DendroBLASTTress - object
    2. DendrobBLAST:
        - read scores
        - run_analysis: Get distance matrices, do trees
    3. Root species tree
    4. Reconciliation/Orthologs
    5. Clean up
    
    Variables:
    - ogSet - all the relevant information about the orthogroups, species etc.
    """
    tree_generation_method = "msa" if q_msa or q_phyldog else "dendroblast"
    stop_after = "seqs" if q_stop_after_seqs else "align" if q_stop_after_align else ""
    file_handler.handler.make_results_directory2(tree_generation_method, stop_after, results_name)    
    """ === 1 === ust = UserSpeciesTree
    MSA:               Sequences    Alignments                        GeneTrees    db    SpeciesTree
    Phyldog:           Sequences    Alignments                        GeneTrees    db    SpeciesTree  
    Dendroblast:                                  DistanceMatrices    GeneTrees    db    SpeciesTree
    MSA (ust):         Sequences    Alignments                        GeneTrees    db
    Phyldog (ust):     Sequences    Alignments                        GeneTrees    db      
    Dendroblast (ust):                            DistanceMatrices    GeneTrees    db        
    """
    q_db_species_tree = False
    if user_species_tree:
        if i_og_restart == 0: printing.print_underline("Using user-supplied species tree")
        sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()  # save it as 'unrooted' but is copied directly to 'rooted' filename
        convert_user_species_tree(user_species_tree, og_set.species_dict(), sp_tree_fn_ids)
    
    if q_msa or q_phyldog:
        """ A. MSA & Tree inference + unrooted species tree"""
        q_less_than_four_species = len(og_set.seqs_info.species_to_use) < 4
        tree_gen = msa_trees.TreesForOrthogroups(program_caller, msa_method, tree_method)
        if (not user_species_tree) and q_less_than_four_species:
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
            dendroblast.write_species_tree_ids_two_three(og_set.seqs_info.species_to_use, sp_tree_fn_ids)
            taxa.rename_tree_taxa(sp_tree_fn_ids, file_handler.handler.get_species_tree_unrooted_fn(True), og_set.species_dict(),
                                q_support=False, q_fix_negatives=True)
        q_do_msa_species_tree = (not q_less_than_four_species) and (not user_species_tree) and (not root_from_previous)
        util.print_time("Starting MSA/Trees")
        seqs_alignments_dirs = tree_gen.do_trees(
            og_set,
            og_set.spec_seq_dict(), 
            og_set.species_dict(), 
            og_set.species_to_use, 
            n_high_parallel, 
            q_stop_after_seqs, 
            q_stop_after_align or q_phyldog, 
            q_do_species_tree=q_do_msa_species_tree,
            q_trim=q_trim,
            i_og_restart=i_og_restart,
            cmd_order=cmd_order,
            method_threads=method_threads,
            method_threads_large=method_threads_large,
            method_threads_small=method_threads_small, 
            threshold=threshold,
            fix_files=fix_files,
            astral=astral,
            dynamic_threads=dynamic_threads,
            n_skip=n_skip,
            msa_min_seq=msa_min_seq
        )
        util.print_time("Done MSA/Trees")
        if q_do_msa_species_tree:
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
        if q_stop_after_seqs:
            print("")
            return
        elif q_stop_after_align:
            print("")
            return
        if q_db_species_tree and not user_species_tree and not q_less_than_four_species and not root_from_previous:
            db = dendroblast.DendroBLASTTrees(og_set, n_low_parallel, n_high_parallel, q_double_blast)
            printing.print_underline("Inferring species tree (calculating gene distances)")
            print("Loading BLAST scores")
            sp_tree_fn_ids = db.species_tree_only()

#         if qPhyldog:
# #            util.print_time("Do species tree for phyldog")
# #            spTreeFN_ids, spTreeUnrootedFN = db.species_tree_only()
#             if userSpeciesTree: 
#                 userSpeciesTree = convert_user_species_tree(userSpeciesTree, ogSet.species_dict(), files.FileHandler.get_species_tree_unrooted_fn())
#                 # not used for subsequent Phyldog steps
#             util.print_time("Starting phyldog")
#             species_tree_ids_labelled_phyldog = wrapper_phyldog.RunPhyldogAnalysis(files.FileHandler.GetPhyldogWorkingDirectory(),
#                                                                                    ogSet.get_i_ogs4(), ogSet.get_ogs_all(), speciesToUse, nHighParallel)
#             spTreeFN_ids = species_tree_ids_labelled_phyldog
    else:
        db = dendroblast.DendroBLASTTrees(og_set, n_low_parallel, n_high_parallel, q_double_blast)
        sp_tree_fn_ids, q_stag = db.run_analysis(user_species_tree == None)
        if user_species_tree != None:
            sp_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
    # if not fix_files:
    #     files.FileHandler.log_working_directory_trees()
    q_species_tree_supports = False if (user_species_tree or q_msa or q_phyldog) else q_stag

    return None if root_from_previous else sp_tree_fn_ids, q_species_tree_supports
