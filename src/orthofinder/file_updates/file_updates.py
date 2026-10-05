from . import ogs, trees
import os
import itertools
import re
import shutil
import tempfile
from collections import defaultdict

from ..utils import files, util, file_io
from ..tools import newick

def update_output_files(
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
        prev_wd=None,
    ):
     
    iSps = list(map(str, sorted(species_to_use)))   # list of strings
    species_names = [sp_ids[i] for i in iSps]

    ## ------------------------ Fix OGs and OG Sequences -------------------------
    old_hog_n0_file = files.FileHandler.WDHierarchicalOrthogroupsFNN0()
    hog_n0_file = files.FileHandler.HierarchicalOrthogroupsFNN0()
    for file in os.listdir(files.FileHandler.GetResultHOGDir()):
        if file.startswith("N"):
            shutil.copy2(
                os.path.join(files.FileHandler.GetResultHOGDir(), file),  
                files.FileHandler.GetLegacyHOGDir()
            )

    hogs_converter(files.FileHandler.GetWorkingDirectory_Write() + "N0.ids.tsv", hog_n0_file, species_names)

    seq_dir = files.FileHandler.GetResultsSeqsDir()
    util.clear_dir(seq_dir)

    ogSet, idDict, name_dictionary, new_ogs = ogs.post_hogs_processing(
        all_seq_ids,
        speciesInfoObj,
        seqsInfo,
        speciesNamesDict,
        options,
        speciesXML,
        q_incremental=q_incremental,
    )


    util.PrintTime("Updating MSA/Trees")

    # ## -------------------------- Fix Resolved Gene Trees and Gene Trees -------------------------
    resolved_trees_working_dir = files.FileHandler.GetOGsReconTreeDir(qResults=True)
    resolved_trees_id_dir = files.FileHandler.GetResolvedTreeIDDir()
    if exist_msa:
        align_dir = files.FileHandler.GetResultsAlignDir()
        align_id_dir = files.FileHandler.GetAlignIDDir()
    else:
        align_dir = None
        align_id_dir = None

    if prev_wd is not None:
        align_id_dir = os.path.join(prev_wd, "Alignments_ids")
    old_hog_n0 = read_hog_file(hog_n0_file)
    # gene names of the trees -> IDs (a name of several genes: per orthogroup)
    spec_seq_id_dict = tree_name_to_id_dict(idDict, old_hog_n0, species_names)
    hog_n0_over4genes = hog_file_over4genes(old_hog_n0, 2)

    del old_hog_n0
    ## get list of unique OG
    unique_ogs = set(d['OG'] for d in hog_n0_over4genes)
    simplified_name_dict = {
        entry[1]: entry[0] 
        for hog_list in name_dictionary.values()
        for entry in hog_list
    }

    unique_ogs =  sorted(unique_ogs)
    tree_file_index = index_files(resolved_trees_working_dir, ".txt")
    # missing = sorted(set(unique_ogs) - set(tree_file_index))
    # if not options.qFastAdd:
    check_missing_ogs(unique_ogs, tree_file_index, resolved_trees_working_dir)

    fasta_file_index = index_files(align_id_dir, ".fa") if align_id_dir is not None else {}

    # hog_index = {
    #     unique_og: [row for row in hog_n0_over4genes if unique_og in row["OG"]]
    #     for unique_og in unique_ogs
    # }

    hog_index = build_hog_index(unique_ogs, hog_n0_over4genes)

    trees.post_ogs_processing(
        unique_ogs,
        resolved_trees_id_dir,
        hog_index, 
        simplified_name_dict, 
        idDict,
        spec_seq_id_dict,
        species_names, 
        nprocess,
        tree_file_index,
        fasta_file_index,  
        align_dir=align_dir,
        min_seq=options.min_seq,
        exist_msa=exist_msa
    )


    id_index = index_files(resolved_trees_id_dir, ".txt")

    expected_ogs = [
        simplified_name_dict[i["HOG"]]
        for i in hog_n0_over4genes
    ]
    # if not options.qFastAdd:
    check_missing_ogs(expected_ogs, id_index, resolved_trees_id_dir)

    ## ----------------------- Fix MSA Alignments --------------------------

    if exist_msa:
        CopyTinyAlignments(align_id_dir, align_dir, name_dictionary, idDict)

    return ogSet, new_ogs

_NEWICK_UNSAFE = re.compile("[" + newick._ILEGAL_NEWICK_CHARS + "]")


def tree_name_to_id_dict(idDict, hog_n0_rows=(), species_names=()):
    """
    The sequence IDs of the gene names on the gene trees' leaves
    (trees.TreeLeafNames), as names appear in the trees.

    The tree writer replaces characters that Newick cannot hold (; [ ] = ...)
    by "_", so "protein [E. coli]" is "protein _E. coli_" in a tree. Both
    forms map to the ID. A name (or tree form) shared by several genes is
    resolved per tree, from the IDs of each orthogroup in hog_n0_rows (the
    converted N0 HOGs: rows with "OG" and one column of IDs per species).
    """
    by_name = defaultdict(set)
    others = {}                      # entries that are not genes ("0" -> species name)
    for seq_id, name in idDict.items():
        if "_" not in seq_id:
            others[name] = seq_id
            continue
        by_name[name].add(seq_id)
        by_name[_NEWICK_UNSAFE.sub("_", name)].add(seq_id)
    name_to_id, ambiguous = {}, {}
    for name, seq_ids in by_name.items():
        if len(seq_ids) == 1:
            name_to_id[name] = next(iter(seq_ids))
        else:
            ambiguous[name] = sorted(seq_ids, key=util.seq_id_key)
    for name, seq_id in others.items():
        if name not in by_name:
            name_to_id[name] = seq_id
    # the orthogroups of the genes with an ambiguous name (only those are needed)
    ambiguous_ids = {seq_id for ids in ambiguous.values() for seq_id in ids}
    og_ids = defaultdict(set)
    if ambiguous_ids:
        for row in hog_n0_rows:
            for col in species_names:
                if row.get(col):
                    og_ids[row["OG"]].update(
                        g for g in (x.strip() for x in row[col].split(",")) if g in ambiguous_ids)
    return trees.TreeLeafNames(name_to_id, ambiguous, dict(og_ids))


_HOG_ID_COLUMNS = ("HOG", "OG", "Gene Tree Parent Clade")


def hogs_converter(hogs_n0_ids_file, hogs_n0_file, species_names):
    """
    Write the N0 HOG file with sequence IDs (hogs_n0_file), for the update of
    the results files, from N0.ids.tsv (hogs_n0_ids_file): the HogWriter
    writes it in the same run, with every HOG of N0 in sequence IDs and the
    same HOG IDs as N0.tsv. The IDs are exact (two genes of a species with the
    same name keep their own), so the HOG files with gene names need not be
    written for this (they are only kept as legacy files, -rmlg).

    The columns are put in the order of the results (species_names); rows are
    read as lists and only their non-empty cells are copied: with many species
    most cells are empty.
    """
    fieldnames = list(_HOG_ID_COLUMNS) + species_names
    temp_name = None
    try:
        with open(hogs_n0_ids_file, newline='') as infile, \
            tempfile.NamedTemporaryFile(
                mode='w', delete=False, newline='', dir=os.path.dirname(hogs_n0_file)
            ) as temp_file:
            temp_name = temp_file.name
            reader = file_io.unquoted_reader(infile)   # HOG files are written unquoted
            writer = file_io.writer(temp_file)
            header = next(reader, None)
            writer.writerow(fieldnames)
            if header is not None:
                unknown = set(header) - set(fieldnames)
                if unknown:
                    raise ValueError("Unexpected columns in %s: %s" % (hogs_n0_ids_file, sorted(unknown)[:5]))
                out_pos = [fieldnames.index(name) for name in header]
                for row in reader:
                    if not row:
                        continue
                    if len(row) > len(header):
                        raise ValueError("Row with more fields than the header in %s" % hogs_n0_ids_file)
                    out = [""] * len(fieldnames)
                    for i in itertools.compress(range(len(row)), row):   # non-empty cells only
                        out[out_pos[i]] = row[i]
                    writer.writerow(out)
    except BaseException:
        if temp_name is not None and os.path.exists(temp_name):
            os.remove(temp_name)
        raise
    os.replace(temp_name, hogs_n0_file)


def read_hog_file(hog_file):
    """
    The rows of a HOG file as dicts holding only the non-empty cells (the HOG,
    OG and Gene Tree Parent Clade columns are always included).

    With many species almost every species cell is empty; storing them all
    would need one dict entry per species per HOG (e.g. 1,200 x 100,000).
    Readers use row.get(species) and treat a missing cell like an empty one.
    """
    hog_n0 = []
    with open(hog_file, newline = '') as csvfile:
        reader = file_io.reader(csvfile)
        header = next(reader, None)
        if header is None:
            return hog_n0
        id_positions = [i for i, name in enumerate(header) if name in _HOG_ID_COLUMNS]
        for row in reader:
            if not row:
                continue                         # DictReader skipped blank lines too
            n = min(len(row), len(header))
            cells = dict.fromkeys(i for i in id_positions if i < n)
            cells.update(dict.fromkeys(itertools.compress(range(n), row)))
            hog_n0.append({header[i]: row[i] for i in sorted(cells)})
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

def CopyTinyAlignments(align_id_dir, align_dir, name_dictionary, idDict):
    for entry in os.scandir(align_id_dir):
        if entry.is_file():
            if '.' not in entry.name:
                continue
            
            old_name = entry.name.rsplit(".", 1)[0]
            entries = name_dictionary.get(old_name)
            if entries is None: 
                continue

            for row in entries:
                if len(row) >= 3 and row[2].strip() == '-':
                    genes_dict = trees.read_fasta(entry.path)
                    trees.write_fasta(align_dir, row[0], genes_dict, idDict)
                    break

def index_files(id_dir, extension=".fa"):
    file_index = {}
    if id_dir is None:
        return file_index
    if os.path.exists(id_dir):
        if os.path.isdir(id_dir):
            for entry in os.scandir(id_dir):
                if entry.is_file() and entry.name.endswith(extension) and entry.name.startswith("OG"):
                    key = entry.name[:-len(extension)]
                    file_index[key] = entry.path
        else:
            with open(id_dir) as reader:
                for line in reader:
                    key, val = line.strip().split(": ", 1)
                    file_index[key] = val
    
    return file_index


def build_hog_index(unique_ogs, hog_rows):
    unique_set = {str(og).strip() for og in unique_ogs}
    hog_index = {og: [] for og in unique_set}
    for row in hog_rows:
        og = str(row.get("OG", "")).strip()
        if og in hog_index:
            row["OG"] = og
            if "HOG" in row and row["HOG"] is not None:
                row["HOG"] = str(row["HOG"]).strip()
            if "Gene Tree Parent Clade" in row and row["Gene Tree Parent Clade"] is not None:
                row["Gene Tree Parent Clade"] = str(row["Gene Tree Parent Clade"]).strip()
            hog_index[og].append(row)
    return hog_index


def check_missing_ogs(ogs_list, index_dict, resolved_trees_dir):
    ogs_set = set(ogs_list)
    index_set = {
        file.rsplit("/", 1)[-1].partition(".")[0]
        for file in index_dict.values()
    }

    missing = ogs_set - index_set

    if missing:
        print(f"ERROR: {len(missing)} ID trees missing in {resolved_trees_dir}")
        for og in missing:
            print("Missing IDs:", og)
        util.Fail()