from . import convert, update_trees
import os
import itertools

from ..utils import file_handler, util, file_io
from ..utils import paths

def update_output_files(
        sp_ids,
        id_sequence_dict,
        species_to_use,
        all_seq_ids,
        species_info_obj,
        seqs_info,
        species_names_dict,
        options,
        nprocess,
        q_incremental=False,
        i_og_restart=0,
        exist_msa=True,
        prev_wd=None,
    ):
     
    i_sps = list(map(str, sorted(species_to_use)))   # list of strings
    species_names = [sp_ids[i] for i in i_sps]

    ## ------------------------ Fix OGs and OG Sequences -------------------------
    # The HOGs of N0 of the MCL orthogroups (the first pass), in sequence IDs: the
    # orthogroups are made from them. (The IDs are exact: two genes of a species with
    # the same name keep their own.)
    hog_n0_file = file_handler.handler.n0_mcl_ids_fn()

    seq_dir = file_handler.handler.get_results_seqs_dir()
    paths.clear_dir(seq_dir)

    og_set, id_dict, name_dictionary, new_ogs = convert.post_hogs_processing(
        all_seq_ids,
        species_info_obj,
        seqs_info,
        species_names_dict,
        options,
        q_incremental=q_incremental,
    )


    util.print_time("Updating MSA/Trees")

    # ## -------------------------- Fix Resolved Gene Trees and Gene Trees -------------------------
    resolved_trees_working_dir = file_handler.handler.get_ogs_recon_tree_dir(q_results=True)
    resolved_trees_id_dir = file_handler.handler.get_resolved_tree_id_dir()
    if exist_msa:
        align_dir = file_handler.handler.get_results_align_dir()
        align_id_dir = file_handler.handler.get_align_id_dir()
    else:
        align_dir = None
        align_id_dir = None

    if prev_wd is not None:
        align_id_dir = os.path.join(prev_wd, "Alignments_ids")
    old_hog_n0 = read_hog_file(hog_n0_file)
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

    update_trees.post_ogs_processing(
        unique_ogs,
        resolved_trees_id_dir,
        hog_index, 
        simplified_name_dict, 
        id_dict,
        None,       # (the trees have the gene IDs: tree_analyser.write_resolved_tree)
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
        copy_tiny_alignments(align_id_dir, align_dir, name_dictionary, id_dict)

    return og_set, new_ogs

_HOG_ID_COLUMNS = ("HOG", "OG", "Gene Tree Parent Clade")


def read_hog_file(hog_file):
    """
    The rows of a HOG file (written unquoted) as dicts holding only the non-empty
    cells (the HOG, OG and Gene Tree Parent Clade columns are always included).

    With many species almost every species cell is empty; storing them all
    would need one dict entry per species per HOG (e.g. 1,200 x 100,000).
    Readers use row.get(species) and treat a missing cell like an empty one.
    """
    hog_n0 = []
    with open(hog_file, newline = '', encoding="utf-8") as csvfile:
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

def copy_tiny_alignments(align_id_dir, align_dir, name_dictionary, id_dict):
    """
    The alignments of the orthogroups of 2 or 3 genes (no tree: gene tree parent
    clade "-") in the results, under their new names. Only the alignments,
    OG<i>.fa: the tree program may write its own files next to them (e.g.
    IQ-TREE's OG<i>.iqtree, .log, .mldist, .bionj), which are not alignments.
    """
    for old_name, path in sorted(index_files(align_id_dir, ".fa").items()):
        entries = name_dictionary.get(old_name)
        if entries is None:
            continue

        for row in entries:
            if len(row) >= 3 and row[2].strip() == '-':
                genes_dict = update_trees.trim_gap_columns(update_trees.read_fasta(path))
                update_trees.write_fasta(align_dir, row[0], genes_dict, id_dict)
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
            with open(id_dir, encoding="utf-8") as reader:
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
        util.fail()