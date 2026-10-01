import csv
import glob
import gzip
import os
import string
import random
from collections import defaultdict
from typing import Optional

try:
    from rich import print
except ImportError:
    ...
import numpy as np

from . import sample_genes
from ..tools import mcl, tree
from ..utils import util, files, parallel_task_manager, fasta_processor, blast_file_processor
from ..run import run_commands
from . import orthogroups_set


class XcelerateConfig(object):
    def __init__(self):
        self.n_for_profiles: Optional[int] = (
            10  # 10 is a good value if using this option
        )


xcelerate_config = XcelerateConfig()


def check_for_orthoxcelerate(input_dir, speciesInfoObj):
    # Add any specific checks required here
    if speciesInfoObj.speciesToUse != list(range(speciesInfoObj.nSpAll)):
        print(
            "ERROR: Removing species from 'core' results directory is not supported for an --assign analysis."
        )
        return False
    return True


def prepare_accelerate_database(
    min_seq,
    input_dir,
    wd_list,
    nSpAll,
    speciesInfoObj,
    options,
    prog_caller,
    tree_program="fasttree",
):
    if xcelerate_config.n_for_profiles is None:
        # create_shoot_db.create_full_database(input_dir, q_ids=True, subtrees_dir="")
        fn_diamond_db, q_hogs = create_profiles_database(
            min_seq,
            input_dir,
            wd_list,
            nSpAll,
            speciesInfoObj,
            options,
            prog_caller,
            selection="all",
            q_ids=True,
            subtrees_dir="",
            tree_program=tree_program,
        )
    else:
        fn_diamond_db, q_hogs = create_profiles_database(
            min_seq,
            input_dir,
            wd_list,
            nSpAll,
            speciesInfoObj,
            options,
            prog_caller,
            selection="kmeans",
            q_ids=True,
            n_for_profile=xcelerate_config.n_for_profiles,
            subtrees_dir="",
            tree_program=tree_program,
        )
    return fn_diamond_db, q_hogs


def _ogs_from_diamond_results_fast(fn_og_results_out, q_ignore_sub):
    """
    Same result as ogs_from_diamond_results (below), computed with arrays.
    Returns None if pandas is unavailable or the file has an unexpected layout.

    Lines are "<sp>_<seq>  <og>_<sp>_<seq>  ...  <evalue>  <bitscore>". For
    each query gene the best hit is the smallest (evalue, og, species), with og
    and species compared as text, exactly as the sorted() of tuples did.
    """
    try:
        n_fields = blast_file_processor.CountFieldsFirstLine(fn_og_results_out)
        if n_fields is None:
            return defaultdict(), defaultdict(lambda: defaultdict(int))
        if n_fields < 7:
            return None
        columns = [0, 1, 2, 3, n_fields - 2]
        # Orthogroup IDs are written as "%07d": read as numbers if they all are
        # (then numeric order is text order), otherwise as text.
        try:
            cols = blast_file_processor.ReadHitColumns(
                fn_og_results_out, n_fields, columns,
                [np.int64, np.int64, np.int64, np.int64, np.float64],
            )
            if cols is not None and len(cols[2]) and (cols[2].min() < 0 or cols[2].max() >= 10 ** 7):
                raise ValueError("orthogroup IDs not all 7 digits")
            og_numeric = True
        except ValueError:
            cols = blast_file_processor.ReadHitColumns(
                fn_og_results_out, n_fields, columns,
                [np.int64, np.int64, str, np.int64, np.float64],
            )
            og_numeric = False
    except (ValueError, OSError):
        return None
    if cols is None:
        return None
    q_sp, q_seq, og, h_sp, score = cols
    og_assignments = defaultdict()
    species_closest_hits = defaultdict(lambda: defaultdict(int))
    if len(score) == 0:
        return og_assignments, species_closest_hits

    # Rank orthogroup and species IDs by their text order (for tie-breaking).
    # Hash the IDs first (factorize) and sort only the distinct ones: sorting
    # one Python string per hit would cost more than the parsing saves.
    if og_numeric:
        og_values, og_index = np.unique(og, return_inverse=True)
        og_text = np.array(["%07d" % x for x in og_values], dtype=object)
    else:
        og_codes, og_distinct = blast_file_processor.pd.factorize(og, sort=False)
        og_distinct = np.asarray(og_distinct, dtype=object)
        if q_ignore_sub:
            og_distinct = np.array([x.split(".", 1)[0] for x in og_distinct], dtype=object)
        og_text, remap = np.unique(og_distinct, return_inverse=True)
        og_index = remap[og_codes]
    sp_values, sp_index = np.unique(h_sp, return_inverse=True)
    sp_text = np.array([str(x) for x in sp_values], dtype=object)
    sp_rank = np.empty(len(sp_text), dtype=np.int64)
    sp_rank[np.argsort(sp_text, kind="stable")] = np.arange(len(sp_text))

    gene_key = q_sp * (1 << 32) + q_seq
    order = np.lexsort((sp_rank[sp_index], og_index, score, gene_key))
    first_in_group = np.ones(len(order), dtype=bool)
    first_in_group[1:] = gene_key[order][1:] != gene_key[order][:-1]
    best = order[first_in_group]                     # one row per gene, by gene key
    # Report genes in the order they first appear in the file, as before.
    _, first_seen = np.unique(gene_key, return_index=True)
    for b in best[np.argsort(first_seen, kind="stable")]:
        gene = "%d_%d" % (q_sp[b], q_seq[b])
        og_assignments[gene] = og_text[og_index[b]]
        species_closest_hits[str(q_sp[b])][sp_text[sp_index[b]]] += 1
    return og_assignments, species_closest_hits


def ogs_from_diamond_results(fn_og_results_out, q_ignore_sub=False):
    """
    Get the OG based on the DIAMOND results
    Args:
        fn_og_results_out - fn of compressed DIAMOND results
        q_ignore_sub - ignore any subtrees and just look at overall OGs
    Returns:
        iog        : List[str] "int.int" or "int" of ordered, ambiguous assignments
        species_closest_hits : Dict[str, Dict[str, int]] search_species -> Dict[hit_species, n_closest_hits]
    Info:
        Hits to other groups are returned if np.log10 difference is less than 10
        and -np.log10 score > difference.
    """
    fast = _ogs_from_diamond_results_fast(fn_og_results_out, q_ignore_sub)
    if fast is not None:
        return fast

    ogs_sp_hits = defaultdict(list)
    scores_all_genes = defaultdict(list)

    if fn_og_results_out.endswith(".gz"):
        with gzip.open(fn_og_results_out, "rt") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for line in reader:
                gene = line[0]
                og, sp, _ = line[1].split("_")
                ogs_sp_hits[gene].append(
                    (og.split(".", 1)[0] if q_ignore_sub else og, sp)
                )
                scores_all_genes[gene].append(float(line[-2]))
    else:
        with open(fn_og_results_out, "rt") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for line in reader:
                gene = line[0]
                og, sp, _ = line[1].split("_")
                ogs_sp_hits[gene].append(
                    (og.split(".", 1)[0] if q_ignore_sub else og, sp)
                )
                scores_all_genes[gene].append(float(line[-2]))

    all_genes = list(ogs_sp_hits.keys())
    og_assignments = defaultdict()
    species_closest_hits = defaultdict(lambda: defaultdict(int))
    for gene in all_genes:
        ogs_sp = ogs_sp_hits[gene]
        scores = scores_all_genes[gene]
        sortedTuples = sorted(zip(scores, ogs_sp))
        # Only the best-scoring orthogroup is used below. (The previous code
        # also de-duplicated the orthogroup list with a quadratic scan and
        # computed log-scores that were never used.)
        if len(sortedTuples) > 0:
            og_assignments[gene] = sortedTuples[0][1][0]
            # print(sortedTuplesples)
            sp_hit = sortedTuples[0][1][1]
            # sys.exit()
            search_sp = gene.split("_")[0]
            species_closest_hits[search_sp][sp_hit] += 1
            # scores_all_genes[gene] = scores[0]
    return og_assignments, species_closest_hits


def get_original_orthogroups():
    wd_input_clusters = files.FileHandler.GetWorkingDirectory1_Read()[
        1
    ]  # first is the current working directory, second is the most recent of the input directories
    fn_clusters = glob.glob(wd_input_clusters + "clusters*_id_pairs.txt")
    if len(fn_clusters) == 0:
        print("ERROR: Couldn't find previous orthogroups in %s" % wd_input_clusters)
        util.Fail()
    elif len(fn_clusters) > 1:
        print("WARNING: Found multiple orthogroup files: %s" % fn_clusters)
    fn_clusters = fn_clusters[0]
    ogs = mcl.GetPredictedOGs(fn_clusters)
    return ogs


def _ogs_from_diamond_results_plain(fn):
    """ogs_from_diamond_results with plain dicts, so it can be returned from a worker."""
    ogs_all_genes, species_closest_hits = ogs_from_diamond_results(fn)
    return dict(ogs_all_genes), {k: dict(v) for k, v in species_closest_hits.items()}


def assign_genes(results_files, n_processes=1):
    """
    Returns OGs with the new species added + species group: Dict[query_species, closest_species]

    The results files (one per new species) are independent and are read in
    parallel; they are combined in the same order as before.
    """
    ogs = defaultdict(set)
    species_closest_hits_totals = defaultdict(lambda: defaultdict(int))
    if n_processes > 1 and len(results_files) > 1:
        per_file = parallel_task_manager.ParallelMap(
            _ogs_from_diamond_results_plain, results_files, n_processes
        )
    else:
        per_file = [_ogs_from_diamond_results_plain(fn) for fn in results_files]
    for ogs_all_genes, species_closest_hits in per_file:
        for query_species, hits in species_closest_hits.items():
            for hit_species, count in hits.items():
                species_closest_hits_totals[query_species][hit_species] += count
        for gene, og in ogs_all_genes.items():
            ogs[int(og)].add(gene)
    species_group = dict()
    for query_species, hits in species_closest_hits_totals.items():
        species_group[query_species] = max(hits, key=hits.get)
        # print((query_species, max(hits, key=hits.get), hits))
    return ogs, species_group


# def write_all_orthogroups(ogs: List[Set[str]], ogs_new_species: Dict[int, Set[str]], ogs_clade_specific: List[List[Set[str]]]):
def write_all_orthogroups(ogs, ogs_new_species, ogs_clade_specific_lists, restart_index=None):
    """
    Add the new species' genes and the clade-specific orthogroups to ogs (in
    place), drop single-gene orthogroups whose gene is now in a larger one, and
    write the clusters file.

    restart_index: the position in ogs (before this call) from which
    orthogroups still need trees. Removing single-gene orthogroups shifts the
    positions of everything after them, so if restart_index is given, the
    adjusted value is returned as well: (clusters filename, restart index).
    Without the adjustment the first clade-specific orthogroups would be
    skipped (no trees, no HOGs) and their genes reported as unassigned.
    """
    for iog, genes in ogs_new_species.items():
        ogs[iog].update(genes)

    for clade_ogs in ogs_clade_specific_lists:
        ogs.extend(clade_ogs)

    assigned_genes = {gene for og in ogs if len(og) > 1 for gene in og}
    keep = [len(og) != 1 or og.isdisjoint(assigned_genes) for og in ogs]
    if restart_index is not None:
        restart_index = sum(keep[:restart_index])
    ogs[:] = [og for og, k in zip(ogs, keep) if k]

    _, clustersFilename_pairs = files.FileHandler.CreateUnusedClustersFN()
    mcl.write_updated_clusters_file(ogs, clustersFilename_pairs)
    if restart_index is None:
        return clustersFilename_pairs
    return clustersFilename_pairs, restart_index


def _profile_genes_for_chunk(args):
    """
    Select the profile genes for a chunk of orthogroups (runs in a worker).
    Returns [(og_id_full, genes, q_subtrees)] in orthogroup order.
    """
    tasks, wd, subtrees_dir, pat_super, pat_sub_msa_glob, selection, n_for_profile_max = args
    out = []
    for iog, og in tasks:
        # Seeded per orthogroup: the result does not depend on the order or
        # the worker in which orthogroups are processed.
        rng = random.Random(iog)
        og_id = "%07d" % iog
        q_subtrees = bool(subtrees_dir) and os.path.exists(pat_super % iog)
        fn_msa = wd + "Alignments_ids/OG%07d.fa" % iog
        if q_subtrees:
            print("Subtrees: %d" % iog)
            fns_msa = list(glob.glob(pat_sub_msa_glob % iog))
        elif os.path.exists(fn_msa):
            fns_msa = [fn_msa]
        else:
            fns_msa = [wd + "Sequences_ids/OG%07d.fa" % iog]
        for fn in fns_msa:
            if not os.path.exists(fn):
                print("File does not exist, skipping: %s" % fn)
                continue
            i_part = os.path.basename(fn).rsplit(".", 2)[1] if q_subtrees else None
            fw_temp = fasta_processor.FastaWriter(fn)
            n_in_og = len(fw_temp.SeqLists)
            # Per orthogroup: must not lower the number used for later orthogroups.
            n_for_profile = min(n_in_og, n_for_profile_max)
            if selection == "kmeans" and len(fw_temp.SeqLists) > n_for_profile:
                fn_temp = None
                if q_subtrees:
                    # MSA needs to be modified
                    letters = string.ascii_lowercase
                    fn_temp = (
                        "/tmp/shoot_db_create"
                        + "".join(random.choice(letters) for i in range(6))
                        + os.path.basename(fn)
                    )
                    fw_temp.WriteSeqsToFasta(
                        [
                            g
                            for g in fw_temp.SeqLists
                            if not g.startswith("SHOOTOUTGROUP_")
                        ],
                        fn_temp,
                    )
                    fn = fn_temp
                # Don't trim as OrthoFinder has already trimmed by default
                s = sample_genes.select_from_aligned(fn, n_for_profile, q_trim=False, rng=rng)
                if fn_temp is not None:
                    os.remove(fn_temp)
            elif selection == "kmeans" or selection == "random":
                genes = og
                if q_subtrees:
                    genes = [
                        g
                        for g in fasta_processor.FastaWriter(fn).SeqLists
                        if not g.startswith("SHOOTOUTGROUP_")
                    ]
                s = sorted(genes)
                s = rng.sample(s, min(n_for_profile, len(s)))
            else:
                s = [
                    g
                    for g in fasta_processor.FastaWriter(fn).SeqLists.keys()
                    if not g.startswith("SHOOTOUTGROUP_")
                ]
            og_id_full = og_id + "." + i_part if q_subtrees else og_id
            out.append((og_id_full, s, q_subtrees))
    return out


def create_profiles_database(
    min_seq,
    din,
    wd_list,
    nSpAll,
    speciesInfoObj,
    options,
    prog_caller,
    selection="kmeans",
    n_for_profile=20,
    q_ids=True,
    subtrees_dir="",
    q_hogs=False,
    tree_program="fasttree",
):
    """
    Create a fasta file with profile genes from each orthogroup
    Args:
        din - Input OrthoFinder results directory
        selection - "kmeans|random|all"  - method to use for sampling genes from orthogroups
        n_for_profile - The number of genes to use from each orthogroup, when available
        q_ids - Convert subtrees (with user gene accessions) back to IDs for profiles database
        q_hogs - The use of HOGs is currently not recommended due to the divergence of the resulting
                 orthogroups from an analysis with the original method
    Notes:
    If the trees have been split into subtrees then profiles will be created for
    the subtrees instead.
    """
    if selection not in ("kmeans", "random", "all"):
        raise RuntimeError("selection method '%s' is not defined" % selection)
    wd = din + "WorkingDirectory/"
    if subtrees_dir:
        subtrees_label = "." + os.path.split(subtrees_dir)[1]
        pat_super = din + subtrees_dir + "/super/OG%07d.super.tre"
        pat_sub_msa_glob = din + subtrees_dir + "/msa_sub/OG%07d.*.fa"
    else:
        subtrees_label = ""
    fn_base = "profile_sequences.hogs" if q_hogs else "profile_sequences"
    if subtrees_label:
        if selection == "all":
            fn_fasta = wd + fn_base + ".%s.all.fa" % subtrees_label
        else:
            fn_fasta = (
                wd
                + fn_base
                + ".%s.%d_%s.fa" % (subtrees_label, n_for_profile, selection)
            )
    else:
        if selection == "all":
            fn_fasta = wd + fn_base + ".all.fa"
        else:
            fn_fasta = wd + fn_base + ".%d_%s.fa" % (n_for_profile, selection)

    fn_diamond_db = fn_fasta + ".dmnd"
    if os.path.exists(fn_diamond_db) and (
        options.search_program.split("_", 1)[0] in ["diamond", "blastp", "blastn"]
    ):
        # print("Profiles database already exists and will be reused: %s" % fn_diamond_db)
        print("Profiles database already exists and will be reused: ")
        print(f"[dark_cyan]{fn_diamond_db}[dark_cyan]")
        return fn_diamond_db, q_hogs

    og_set = orthogroups_set.OrthoGroupsSet(
        min_seq, wd_list, list(range(nSpAll)), nSpAll, True, tree_program=tree_program
    )
    ids = og_set.Spec_SeqDict()
    ids_rev = {v: k for k, v in ids.items()}
    if q_hogs:
        try:
            ids_simple = og_set.SequenceDict()
            ids_simple_rev = {v: k for k, v in ids_simple.items()}
            ogs = read_hogs(din, "N0", ids_simple_rev)
        except RuntimeError:
            print(
                "ERROR: Cannot read HOGs file, please report this error: https://github.com/davidemms/OrthoFinder/issues"
            )
            q_hogs = False
            print("WARNING: Using MCL-based orthogroups as a fall-back")
            # util.Fail()
    if not q_hogs:
        clusters_filename = list(glob.glob(wd + "clusters_OrthoFinder*id_pairs.txt"))
        if len(clusters_filename) == 0:
            print("ERROR: Can't find %s" % wd + "clusters_OrthoFinder*id_pairs.txt")
        ogs = mcl.GetPredictedOGs(clusters_filename[0])
        fn_fasta = fn_fasta[:-7] + ".fa"
    fw = fasta_processor.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    # print("WARNING: Check all gene names, can't start with '__'")
    # If there are subtrees then we need to convert their IDs in the profile file
    # back to internal IDs

    # Each orthogroup is independent, so they are processed in chunks in a
    # process pool; the results are joined in orthogroup order, so the profile
    # file is the same as with a single process.
    n_chunks_target = max(1, 8 * max(1, options.nProcessAlg))
    chunk_size = max(1, min(200, -(-len(ogs) // n_chunks_target)))
    chunks = [
        (
            [(iog, ogs[iog]) for iog in range(start, min(start + chunk_size, len(ogs)))],
            wd, subtrees_dir, pat_super if subtrees_dir else None,
            pat_sub_msa_glob if subtrees_dir else None, selection, n_for_profile,
        )
        for start in range(0, len(ogs), chunk_size)
    ]
    progressbar, task = util.get_progressbar(len(chunks))

    def progress(n_done):
        if n_done == 0:
            progressbar.start()   # only once the workers exist
        else:
            progressbar.update(task, completed=n_done)

    try:
        chunk_results = parallel_task_manager.ParallelMap(
            _profile_genes_for_chunk, chunks, options.nProcessAlg, progress=progress
        )
    finally:
        progressbar.stop()

    for result in chunk_results:
        for og_id_full, s, q_subtrees in result:
            if q_ids and q_subtrees:
                s = [ids_rev[ss] for ss in s]
            seq_write.extend(s)
            for ss in s:
                seq_convert[ss] = og_id_full + "_" + ss
    fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    # parallel_task_manager.RunCommand(" ".join(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db]), qPrintOnError=True, qPrintStderr=False)
    run_commands.CreateSearchDatabases(
        speciesInfoObj,
        options,
        prog_caller,
        core_infile=fn_fasta,
        core_outfile=fn_diamond_db,
    )

    return fn_diamond_db, q_hogs


class DummyIDs:
    def __getitem__(self, arg):
        return arg


def read_hogs(din, hog_name, ids_rev=None):
    fn_ids = din + "WorkingDirectory/%s.ids.tsv" % hog_name
    if os.path.exists(fn_ids):
        fn = fn_ids
        ids_rev = DummyIDs()
    elif ids_rev is None:
        return []
    else:
        fn = os.path.join(
            din, "Phylogenetic_Hierarchical_Orthogroups/%s.tsv" % hog_name
        )
    if not os.path.exists(fn):
        print("ERROR: %s does not exist" % fn)
        raise RuntimeError
    ogs = []
    with open(fn, util.csv_read_mode) as infile:
        reader = csv.reader(infile, delimiter="\t")
        next(reader)  # header
        for line in reader:
            ogs.append([])
            for species in line[3:]:
                genes = species.split(", ")
                genes = [ids_rev[g] for g in genes if g != ""]
                ogs[-1].extend(genes)
            ogs[-1] = set(ogs[-1])
    return ogs


def write_unassigned_fasta(ogs_orig_list, ogs_new_genes, speciesInfoObj):
    if ogs_new_genes is not None:
        assigned_genes = set.union(*ogs_new_genes.values())
        assigned_genes.update(set.union(*[og for og in ogs_orig_list if len(og) > 1]))
    else:
        assigned_genes = set.union(*[og for og in ogs_orig_list if len(og) > 1])
    # Write out files for all unassigned genes
    # iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    # for iSp in iSpeciesNew:
    n_unassigned = []
    for iSp in range(speciesInfoObj.nSpAll):
        fw = fasta_processor.FastaWriter(files.FileHandler.GetSpeciesFastaFN(iSp))
        unassigned = set(fw.SeqLists.keys()).difference(assigned_genes)
        n_unassigned.append(len(unassigned))
        fw.WriteSeqsToFasta(
            unassigned,
            files.FileHandler.GetSpeciesUnassignedFastaFN(iSp, qForCreation=True),
        )
    return n_unassigned


def sample_random(og, n_max):
    """
    Sample min(n_max, |og|) genes randomly from clade
    Args:
        og - set of strings
        n_max - max number of genes to sample
    Returns:
        genes - list of genes
    """
    return random.sample(sorted(og), min(n_max, len(og)))


# A new-species clade at least this large gets a warning at run time.
LARGE_CLADE_SPECIES = 50


def get_new_species_clades(rooted_species_tree_fn, core_species_ids, n_core_species=2):
    """
    New-species clades: the largest subtrees of the rooted species tree that
    contain at most n_core_species core species. Unassigned genes of the new
    species are clustered within each clade; every ordered pair of its species
    is searched (k^2 searches for k species), so the clade sizes, and with them
    the choice of core species, determine much of the runtime.

    Args:
        rooted_species_tree_fn - ids format
        core_species_ids: Set[str]
        n_core_species - maximum number of core species in a 'new clade of species'. 1 is the minimum, but 2 makes more
                         allowances for gene loss in one core species & still recovering associated orthogroup
    Returns:
        a list of sorted lists of species IDs. Groups containing only core
        species are not new-species clades and are left out.
    """
    # Compare IDs as strings throughout: the leaf names are strings, so the
    # core IDs must be too (an earlier version compared integers with strings,
    # so no group was ever left out).
    core_species_ids = set(map(str, core_species_ids))
    t = tree.Tree(rooted_species_tree_fn, format=1)

    def is_new_clade(node):
        return len(core_species_ids.intersection(node.get_leaf_names())) <= n_core_species

    species_clades = []
    for n in t.get_leaves(is_leaf_fn=is_new_clade):
        names = [str(name) for name in n.get_leaf_names()]
        if set(names).difference(core_species_ids):
            species_clades.append(list(sorted(map(int, names))))
    return species_clades


def clade_costs(clades, core_species, n_genes=None):
    """
    Per-clade size and number of species-pair searches, largest first.

    n_genes, if given, maps species -> number of unassigned genes; species
    without unassigned genes are not searched (as in the run). Without it the
    numbers are an upper bound.

    Returns a list of dicts: species, n_species, n_core, n_searches.
    """
    core_species = set(map(str, core_species))
    rows = []
    for clade in clades:
        searched = [sp for sp in clade if n_genes is None or n_genes.get(sp, 0) > 0]
        rows.append({
            "species": list(clade),
            "n_species": len(searched),
            "n_core": len(core_species.intersection(map(str, clade))),
            "n_searches": len(searched) ** 2,
        })
    rows.sort(key=lambda r: -r["n_species"])
    return rows


def report_clades(rows, species_names=None, large=None, max_listed=5, print_fn=print):
    """
    Print a summary of the clade-specific step and a warning for each clade of
    at least `large` species. Returns the list of warning messages.
    """
    if large is None:
        large = LARGE_CLADE_SPECIES

    def name(sp):
        return species_names.get(str(sp), str(sp)) if species_names else str(sp)

    total = sum(r["n_searches"] for r in rows)
    n_species = sum(r["n_species"] for r in rows)
    print_fn(
        "Clade-specific orthogroups: %d clade%s, %d species, %d species-pair searches"
        % (len(rows), "" if len(rows) == 1 else "s", n_species, total)
    )
    if rows:
        sizes = [r["n_species"] for r in rows]
        print_fn(
            "Clade sizes: largest %d, median %d; the largest clade accounts for %.0f%% of the searches"
            % (sizes[0], sorted(sizes)[len(sizes) // 2],
               100.0 * rows[0]["n_searches"] / max(1, total))
        )
    warnings = []
    for r in rows:
        if r["n_species"] < large:
            break
        members = [name(sp) for sp in r["species"]]
        listed = ", ".join(members[:max_listed])
        if len(members) > max_listed:
            listed += ", ... (%d more)" % (len(members) - max_listed)
        warnings.append(
            "WARNING: large new-species clade: %d species (%d core) -> %d searches "
            "and a clustering step that grows faster than that. Adding a core "
            "species within this part of the species tree would split it. "
            "Species: %s" % (r["n_species"], r["n_core"], r["n_searches"], listed)
        )
    for w in warnings:
        print_fn(w)
    return warnings
