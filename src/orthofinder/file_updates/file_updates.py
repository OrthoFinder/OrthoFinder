from . import ogs, trees
import os
import tempfile
import csv
import shutil
from collections import defaultdict

from ..utils import files, util
from ..utils.util import printer

def update_output_files(
        working_dir,
        sp_ids,
        id_sequence_dict,
        species_to_use,
        all_seq_ids,
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
    
    sequence_id_dict = defaultdict(set)
    for key, value in id_sequence_dict.items():
        sequence_id_dict[value].add(key)

    iSps = list(map(str, sorted(species_to_use)))   # list of strings
    species_names = [sp_ids[i] for i in iSps]

    species_id_dict = {
        val: key 
        for key, val in sp_ids.items()
    }

    ## ------------------------ Fix OGs and OG Sequences -------------------------
    hog_n0_file = files.FileHandler.HierarchicalOrthogroupsFNN0()
    hogs_converter(hog_n0_file, sequence_id_dict, species_id_dict, species_names)

    # seq_id_dir = files.FileHandler.GetSeqsIDDir()
    seq_dir = files.FileHandler.GetResultsSeqsDir()
    
    # seq_id_dir2 = os.path.join(working_dir, "Sequences_ids2")
    # os.makedirs(seq_id_dir2, exist_ok=True)
    # shutil.copytree(seq_id_dir, seq_id_dir2, dirs_exist_ok=True)

    ## Clean dirs 
    # clear_dir(seq_id_dir)
    util.clear_dir(seq_dir)

    ogSet, treeGen, idDict, new_ogs, name_dictionary = ogs.post_hogs_processing(
        # single_ogs_list,
        all_seq_ids,
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
    # shutil.rmtree(seq_id_dir)
    # shutil.move(seq_id_dir2, seq_id_dir)

    # util.PrintTime("Updating MSA/Trees")

    # ## -------------------------- Fix Resolved Gene Trees and Gene Trees -------------------------
    resolved_trees_working_dir = files.FileHandler.GetOGsReconTreeDir(qResults=True)
    update_filenames(resolved_trees_working_dir, name_dictionary)

    resolved_trees_id_dir = files.FileHandler.GetResolvedTreeIDDir()
    update_filenames(resolved_trees_id_dir, name_dictionary)

    # hog_msa_dir = files.FileHandler.GetHOGMSADir()
    # align_id_dir = files.FileHandler.GetAlignIDDir()

    # shutil.move(align_id_dir, align_id_dir)
    # update_filenames(align_id_dir, name_dictionary)
   
    align_dir = files.FileHandler.GetResultsAlignDir()
    update_filenames(align_dir, name_dictionary)

    # align_id_dir = None
    # align_id_dir2 = None
    # if exist_msa:
    #     align_id_dir = files.FileHandler.GetAlignIDDir()
    #     align_id_dir2 = os.path.join(working_dir, "Alignments_ids2")
    #     os.makedirs(align_id_dir2, exist_ok=True)
    #     shutil.copytree(align_id_dir, align_id_dir2, dirs_exist_ok=True)

    # clear_dir(align_id_dir)

    # tree_id_dir = files.FileHandler.GetOGsTreeDir(qResults=False)
    # tree_dir = files.FileHandler.GetOGsTreeDir(qResults=True)

    # ## Clean dirs 
    # clear_dir(tree_id_dir)
    # clear_dir(tree_dir)

    # old_hog_n0 = read_hog_n0_file(hog_n0_file)
    # hog_n0_over4genes = hog_file_over4genes(old_hog_n0, options.min_seq)

    # del old_hog_n0
    # ## get list of unique OG
    # unique_ogs = set(d['OG'] for d in hog_n0_over4genes)
    # simplified_name_dict = {
    #     entry[1]: entry[0] 
    #     for hog_list in name_dictionary.values()
    #     for entry in hog_list
    # }
    
    # trees.post_ogs_processing(
    #     unique_ogs,
    #     resolved_trees_working_dir, 
    #     tree_id_dir,
    #     tree_dir,
    #     hog_n0_over4genes, 
    #     simplified_name_dict, 
    #     idDict,
    #     spec_seq_id_dict,
    #     species_names, 
    #     nprocess,
    #     align_id_dir=align_id_dir,
    #     align_id_dir2=align_id_dir2,
    # )
    
    # shutil.rmtree(align_id_dir)
    # shutil.move(align_id_dir2, align_id_dir)

    # clear_dir(resolved_trees_working_dir)

    ## ----------------------- Fix MSA Alignments --------------------------
    # if exist_msa:
    #     align_dir = files.FileHandler.GetResultsAlignDir()
    #     util.clear_dir(align_dir)
    #     iogs_align = [i for i, og in enumerate(new_ogs) if len(og) >= 2 and i >= i_og_restart]
    #     # ids -> accessions
    #     alignmentFilesToUse = [treeGen.GetAlignmentFilename(i) for i in iogs_align]
    #     accessionAlignmentFNs = [treeGen.GetAlignmentFilename(i, True) for i in iogs_align]
    #     treeGen.RenameAlignmentTaxa(alignmentFilesToUse, accessionAlignmentFNs, idDict)

    return ogSet

def hogs_converter(hogs_n0_file, sequence_id_dict, species_id_dict, species_names, rm_N0_ids=True):

    with open(hogs_n0_file, newline='') as infile, \
        tempfile.NamedTemporaryFile(
            mode='w', delete=False, newline='', dir=os.path.dirname(hogs_n0_file)
        ) as temp_file:
        reader = csv.DictReader(infile, delimiter='\t')

        fieldnames = ["HOG", "OG", "Gene Tree Parent Clade"] + species_names
        writer = csv.DictWriter(temp_file, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for row in reader:
            new_row = {
                key: (
                    ", ".join(
                        next(iter(
                            [s for s in sequence_id_dict.get(gene, set()) 
                             if s.split("_")[0] == species_id_dict[key]]
                        ), "")
                        for gene in str(val).split(", ")
                    )
                    if (val is not None and "," in str(val))
                    else next(
                        iter([s for s in sequence_id_dict.get(val, set()) 
                              if s.split("_")[0] == species_id_dict[key]]), ""
                    )
                ) if key not in {'OG', 'Gene Tree Parent Clade', 'HOG'} else val
                for key, val in row.items()
            }
            writer.writerow(new_row)
    

    os.replace(temp_file.name, hogs_n0_file)
    # shutil.copy(hogs_n0_file, os.path.join(os.path.dirname(hogs_n0_file), "N0_ids.tsv"))

def read_hog_n0_file(hog_n0_file):
    hog_n0 = []
    with open(hog_n0_file, newline = '') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        reader.fieldnames = [
            # fieldname.replace('.', '_') 
            fieldname
            for fieldname in reader.fieldnames
        ]
        hog_n0 = [row for row in reader]
    return hog_n0


## Function to get rid of hierarchical orthogroups with < 4 genes
def hog_file_over4genes(hog_n0, min_seq):

    filtered_hog_n0 = []
    for row in hog_n0:
        if "n" in row["Gene Tree Parent Clade"]:
            genes = ', '.join([
                value 
                for key, value in row.items() 
                if key not in {'OG','Gene Tree Parent Clade', 'HOG'} and value
            ]).split(', ')

            if len(genes) >= min_seq:
                filtered_hog_n0.append(row)

    return filtered_hog_n0    

def update_filenames(file_dir, name_dictionary):
 
    for entry in os.scandir(file_dir):
        filename, extension = entry.name.rsplit(".", 1)
        if "_" not in filename:
            continue
        old_name, node_name = filename.split("_")
        names = name_dictionary.get(old_name)
        if names is None:
            continue
        for i in names:
            if i[2] == node_name:
                new_filename = i[0] + "." + extension
                os.rename(entry.path, os.path.join(file_dir, new_filename))
                break




