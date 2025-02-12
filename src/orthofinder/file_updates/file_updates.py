from . import ogs, trees
import os
import tempfile
import ete3
import csv
import shutil
from ..utils import util, files
from ..tools import trees_msa


def update_output_files(
        working_dir,
        id_sequence_dict,
        single_ogs_list,
        speciesInfoObj,
        seqsInfo,
        speciesNamesDict,
        options,
        speciesXML,
        nprocess,
        q_incremental=False,
        i_og_restart=0,
        exist_msa=True,
    ):

    sequence_id_dict = {
        val: key
        for key, val in id_sequence_dict.items()
    }

    ## ------------------------ Fix OGs and OG Sequences -------------------------
    hog_n0_file = files.FileHandler.HierarchicalOrthogroupsFNN0()
    hogs_converter(hog_n0_file, sequence_id_dict)

    seq_id_dir = files.FileHandler.GetSeqsIDDir()
    seq_dir = files.FileHandler.GetResultsSeqsDir()

    ## Clean dirs 
    clear_dir(seq_id_dir)
    clear_dir(seq_dir)

    ogSet, treeGen, idDict, new_ogs, name_dictionary, species_names = post_hogs_processing(
        single_ogs_list,
        speciesInfoObj,
        seqsInfo,
        speciesNamesDict,
        options,
        speciesXML,
        q_incremental=q_incremental,
    )
    spec_seq_id_dict = {
        val: key
        for key, val in idDict.items()
    }

    # ## -------------------------- Fix Resolved Gene Trees -------------------------
    resolved_trees_working_dir = files.FileHandler.GetOGsReconTreeDir(qResults=True)
    trees_converter(resolved_trees_working_dir, spec_seq_id_dict)
    resolved_trees_working_dir2 = os.path.join(working_dir,  "Resolved_Gene_Trees2")
    os.makedirs(resolved_trees_working_dir2, exist_ok=True)
    shutil.copytree(resolved_trees_working_dir, resolved_trees_working_dir2, dirs_exist_ok=True)

    align_id_dir=None,
    align_id_dir2=None,
    if exist_msa:
        align_id_dir = files.FileHandler.GetAlignIDDir()
        align_id_dir2 = os.path.join(working_dir, "Alignments_ids2")
        os.makedirs(align_id_dir2, exist_ok=True)
        shutil.copytree(align_id_dir, align_id_dir2, dirs_exist_ok=True)

    clear_dir(resolved_trees_working_dir)
    clear_dir(align_id_dir)

    old_hog_n0 = read_hog_n0_file(hog_n0_file)
    hog_n0_over4genes = hog_file_over4genes(old_hog_n0)

    del old_hog_n0
    ## get list of unique OG
    unique_ogs = set(d['OG'] for d in hog_n0_over4genes)
    
    simplified_name_dict = {
        entry[1]: entry[0] 
        for hog_list in name_dictionary.values()
        for entry in hog_list
    }
    
    trees.post_ogs_processing(
        unique_ogs,
        resolved_trees_working_dir, 
        resolved_trees_working_dir2, 
        hog_n0_over4genes, 
        simplified_name_dict, 
        species_names, 
        nprocess,
        align_id_dir=align_id_dir,
        align_id_dir2=align_id_dir2,
    )
    
    shutil.rmtree(resolved_trees_working_dir2)
    shutil.rmtree(align_id_dir2)

    # ## -------------------------- Fix Gene Trees -------------------------
    tree_id_dir = files.FileHandler.GetOGsTreeDir(qResults=False)
    tree_dir = files.FileHandler.GetOGsTreeDir(qResults=True)

    ## Clean dirs 
    clear_dir(tree_id_dir)
    clear_dir(tree_dir)

    overwrite_gene_trees(
        resolved_trees_working_dir, 
        spec_seq_id_dict,
        tree_id_dir,
        tree_dir
    )
    clear_dir(resolved_trees_working_dir)

    ## ----------------------- Fix MSA Alignments --------------------------
    if exist_msa:
        align_dir = files.FileHandler.GetResultsAlignDir()
        clear_dir(align_dir)
        iogs_align = [i for i, og in enumerate(new_ogs) if len(og) >= 2 and i >= i_og_restart]
        # ids -> accessions
        alignmentFilesToUse = [treeGen.GetAlignmentFilename(i) for i in iogs_align]
        accessionAlignmentFNs = [treeGen.GetAlignmentFilename(i, True) for i in iogs_align]
        treeGen.RenameAlignmentTaxa(alignmentFilesToUse, accessionAlignmentFNs, idDict)

    return ogSet

def hogs_converter(hogs_n0_file, sequence_id_dict):

    with open(hogs_n0_file, newline='') as infile, \
        tempfile.NamedTemporaryFile(
            mode='w', delete=False, newline='', dir=os.path.dirname(hogs_n0_file)
        ) as temp_file:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = [fieldname.replace('.', '_') for fieldname in reader.fieldnames]
        writer = csv.DictWriter(temp_file, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for row in reader:
            writer.writerow({
                key: (
                    ", ".join([sequence_id_dict.get(gene, gene) for gene in str(val).split(", ")]) 
                    if val is not None and "," in str(val) else sequence_id_dict.get(str(val), str(val))
                ) if key not in {'OG', 'Gene Tree Parent Clade', 'HOG'} else val
                for key, val in row.items()
            })

    os.replace(temp_file.name, hogs_n0_file)


def trees_converter(
        resolved_trees_dir, 
        spec_seq_id_dict,
    ):
    for entry in os.scandir(resolved_trees_dir):
        if entry.is_file():
            try:
                with open(entry.path, "r") as infile:
                    newick_str = infile.read().strip()
                
                tree = ete3.Tree(newick_str, quoted_node_names=True, format=1)
                
                for leaf in tree.iter_leaves():
                    leaf.name = spec_seq_id_dict.get(leaf.name, leaf.name)
                modified_newick_str = tree.write(quoted_node_names=True, format=1)
    
                with open(entry.path, "w") as outfile:
                    outfile.write(modified_newick_str)
            
            except Exception as e:
                print(f"Error processing {entry.name}: {e}")

def overwrite_gene_trees(
        resolved_trees_dir, 
        spec_seq_id_dict,
        tree_id_dir,
        tree_dir,
    ):
    for entry in os.scandir(resolved_trees_dir):
        if entry.is_file():
            try:
                with open(entry.path, "r") as infile:
                    newick_str = infile.read().strip()
                
                tree = ete3.Tree(newick_str, quoted_node_names=True, format=1)
                original_newick_str = tree.write(format=5)
                
                for leaf in tree.iter_leaves():
                    leaf.name = spec_seq_id_dict.get(leaf.name, leaf.name)
                modified_newick_str = tree.write(format=5)

                tree_id_file = os.path.join(tree_id_dir, entry.name)
                tree_file = os.path.join(tree_dir, entry.name)
    
                with open(tree_id_file, "w") as outidfile, open(tree_file, "w") as outfile:
                    outidfile.write(modified_newick_str)
                    outfile.write(original_newick_str)
            
            except Exception as e:
                print(f"Error processing {entry.name}: {e}")

def post_hogs_processing(
        single_ogs_list,
        speciesInfoObj,
        seqsInfo,
        speciesNamesDict,
        options,
        speciesXML,
        q_incremental=False,
    ):
    """
    Write OGs & statistics to results files, write Fasta files.
    Args:
        q_incremental - These are not the final orthogroups, don't write results
    """
    new_ogs, name_dictionary, species_names = \
        ogs.update_ogs(files.FileHandler.HierarchicalOrthogroupsFNN0())
    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    # util.PrintUnderline("Writing orthogroups to file")
    new_ogs.extend(single_ogs_list)
    with open(files.FileHandler.OGsAllIDFN(), "w") as outfile:
        for og in new_ogs:
            outfile.write(", ".join(og) + "\n")
    
    idsDict = ogs.MCL.WriteOrthogroupFiles(
        new_ogs,
        [files.FileHandler.GetSequenceIDsFN()],
        resultsBaseFilename,
    )

    if not q_incremental:
        ogs.MCL.CreateOrthogroupTable(
            new_ogs,
            idsDict,
            speciesNamesDict,
            speciesInfoObj.speciesToUse,
            resultsBaseFilename,
        )

    # Write Orthogroup FASTA files
    ogSet = ogs.OrthoGroupsSet(
        files.FileHandler.GetWorkingDirectory1_Read(), 
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        idExtractor=util.FirstWordExtractor,
    )

    ## ------------------ Fix Orthogroup_Sequences and Sequences_ids --------------------
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(
        files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse
    )
    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs):
        os.mkdir(d_seqs)

    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idsDict, False)
    idDict = ogSet.Spec_SeqDict()
    idDict.update(ogSet.SpeciesDict()) # same code will then also convert concatenated alignment for species tree
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idDict, True)
    
    if not q_incremental:
        # stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
        if options.speciesXMLInfoFN:
            ogs.MCL.WriteOrthoXML(
                speciesXML,
                new_ogs,
                seqsInfo.nSeqsPerSpecies,
                idsDict,
                resultsBaseFilename + ".orthoxml",
                speciesInfoObj.speciesToUse,
            )
        # print("")
        # util.PrintTime("Done orthogroups")
        files.FileHandler.LogOGs()

    return ogSet, treeGen, idDict, new_ogs, name_dictionary, species_names

def read_hog_n0_file(hog_n0_file):
    hog_n0 = []
    with open(hog_n0_file, newline = '') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        reader.fieldnames = [
            fieldname.replace('.', '_') 
            for fieldname in reader.fieldnames
        ]
        hog_n0 = [row for row in reader]
    return hog_n0


## Function to get rid of hierarchical orthogroups with < 4 genes
def hog_file_over4genes(hog_n0):

    filtered_hog_n0 = []
    for row in hog_n0:
        if "n" in row["Gene Tree Parent Clade"]:
            genes = ', '.join([
                value 
                for key, value in row.items() 
                if key not in {'OG','Gene Tree Parent Clade', 'HOG'} and value
            ]).split(', ')

            if len(genes) >= 4:
                filtered_hog_n0.append(row)

    return filtered_hog_n0    

def clear_dir(of3_dir):
    with os.scandir(of3_dir) as entries:
        for entry in entries:
            try:
                if entry.is_file() or entry.is_symlink():
                    os.unlink(entry.path) 
                elif entry.is_dir():
                    shutil.rmtree(entry.path) 
            except Exception as e:
                print(f'Failed to delete {entry.path}. Reason: {e}')
