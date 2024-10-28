from ..utils import util, files
from ..tools import trees_msa, dendroblast, tree, wrapper_phyldog


def ConvertUserSpeciesTree(speciesTreeFN_in, speciesDict, speciesTreeFN_out):
    t = tree.Tree(speciesTreeFN_in, format=1)  
    t.prune(t.get_leaf_names())
    revDict = {v:k for k,v in speciesDict.items()}
    for sp in t:
        sp.name = revDict[sp.name]       
    t.write(outfile=speciesTreeFN_out)

def InferGeneAndSpeciesTrees(ogSet,
                       program_caller,
                       msa_method,
                       tree_method,
                       nHighParallel,
                       nLowParallel,
                       qDoubleBlast,
                       qAddSpeciesToIDs,
                       qTrim,
                       cmd_order="ascending",
                       method_threads=None,
                       method_threads_large=None,
                       method_threads_small=None, 
                       threshold=None,
                       old_version=False,
                       userSpeciesTree = None,
                       qStopAfterSeqs = False,
                       qStopAfterAlign = False,
                       qMSA = False,
                       qPhyldog = False,
                       results_name = "",
                       root_from_previous = False,
                       i_og_restart=0,
):
    """
    1. Setup:
        - ogSet, directories
        - DendroBLASTTress - object
    2. DendrobBLAST:
        - read scores
        - RunAnalysis: Get distance matrices, do trees
    3. Root species tree
    4. Reconciliation/Orthologues
    5. Clean up
    
    Variables:
    - ogSet - all the relevant information about the orthogroups, species etc.
    """
    tree_generation_method = "msa" if qMSA or qPhyldog else "dendroblast"
    stop_after = "seqs" if qStopAfterSeqs else "align" if qStopAfterAlign else ""
    files.FileHandler.MakeResultsDirectory2(tree_generation_method, stop_after, results_name)    
    """ === 1 === ust = UserSpeciesTree
    MSA:               Sequences    Alignments                        GeneTrees    db    SpeciesTree
    Phyldog:           Sequences    Alignments                        GeneTrees    db    SpeciesTree  
    Dendroblast:                                  DistanceMatrices    GeneTrees    db    SpeciesTree
    MSA (ust):         Sequences    Alignments                        GeneTrees    db
    Phyldog (ust):     Sequences    Alignments                        GeneTrees    db      
    Dendroblast (ust):                            DistanceMatrices    GeneTrees    db        
    """
    qDB_SpeciesTree = False
    if userSpeciesTree:
        if i_og_restart == 0: util.PrintUnderline("Using user-supplied species tree")
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()  # save it as 'unrooted' but is copied directly to 'rooted' filename
        ConvertUserSpeciesTree(userSpeciesTree, ogSet.SpeciesDict(), spTreeFN_ids)
    
    if qMSA or qPhyldog:
        """ A. MSA & Tree inference + unrooted species tree"""
        qLessThanFourSpecies = len(ogSet.seqsInfo.speciesToUse) < 4
        treeGen = trees_msa.TreesForOrthogroups(program_caller, msa_method, tree_method)
        if (not userSpeciesTree) and qLessThanFourSpecies:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
            dendroblast.WriteSpeciesTreeIDs_TwoThree(ogSet.seqsInfo.speciesToUse, spTreeFN_ids)
            util.RenameTreeTaxa(spTreeFN_ids, files.FileHandler.GetSpeciesTreeUnrootedFN(True), ogSet.SpeciesDict(),
                                qSupport=False, qFixNegatives=True)
        qDoMSASpeciesTree = (not qLessThanFourSpecies) and (not userSpeciesTree) and (not root_from_previous)
        util.PrintTime("Starting MSA/Trees")
        seqs_alignments_dirs = treeGen.DoTrees(ogSet,
                                               ogSet.Spec_SeqDict(), 
                                               ogSet.SpeciesDict(), 
                                               ogSet.speciesToUse, 
                                               nHighParallel, 
                                               qStopAfterSeqs, 
                                               qStopAfterAlign or qPhyldog, 
                                               qDoSpeciesTree=qDoMSASpeciesTree,
                                               qTrim=qTrim,
                                               i_og_restart=i_og_restart,
                                               cmd_order=cmd_order,
                                               method_threads=method_threads,
                                               method_threads_large=method_threads_large,
                                               method_threads_small=method_threads_small, 
                                               threshold=threshold,
                                               old_version=old_version
                                               )
        util.PrintTime("Done MSA/Trees")
        if qDoMSASpeciesTree:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        if qStopAfterSeqs:
            print("")
            return
        elif qStopAfterAlign:
            print("")
            return
        if qDB_SpeciesTree and not userSpeciesTree and not qLessThanFourSpecies and not root_from_previous:
            db = dendroblast.DendroBLASTTrees(ogSet, nLowParallel, nHighParallel, qDoubleBlast)
            util.PrintUnderline("Inferring species tree (calculating gene distances)")
            print("Loading BLAST scores")
            spTreeFN_ids = db.SpeciesTreeOnly()
        if qPhyldog:
#            util.PrintTime("Do species tree for phyldog")
#            spTreeFN_ids, spTreeUnrootedFN = db.SpeciesTreeOnly()
            if userSpeciesTree: 
                userSpeciesTree = ConvertUserSpeciesTree(userSpeciesTree, ogSet.SpeciesDict(), files.FileHandler.GetSpeciesTreeUnrootedFN())
                # not used for subsequent Phyldog steps
            util.PrintTime("Starting phyldog")
            species_tree_ids_labelled_phyldog = wrapper_phyldog.RunPhyldogAnalysis(files.FileHandler.GetPhyldogWorkingDirectory(),
                                                                                   ogSet.Get_iOGs4(), ogSet.OGsAll(), speciesToUse, nHighParallel)
            spTreeFN_ids = species_tree_ids_labelled_phyldog
    else:
        db = dendroblast.DendroBLASTTrees(ogSet, nLowParallel, nHighParallel, qDoubleBlast)
        spTreeFN_ids, qSTAG = db.RunAnalysis(userSpeciesTree == None)
        if userSpeciesTree != None:
            spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
    files.FileHandler.LogWorkingDirectoryTrees()
    qSpeciesTreeSupports = False if (userSpeciesTree or qMSA or qPhyldog) else qSTAG

    return None if root_from_previous else spTreeFN_ids, qSpeciesTreeSupports