"""Choosing orthogroups for the species tree and concatenating their alignments."""

from collections import Counter, defaultdict

import numpy as np

from ..utils import file_handler, printing, util

try:
    from rich import print
except ImportError:
    ...


def get_multicopy_cutoff(n_species, factor = 0.25, p_max = 0.25):
    """ Get the maximum number of species with multicopy genes for each number of non-single copy genes. A species is non-single
    copy if it is absent or if it is multicopy
    Args:
        nSpecies - number of species in the analysis
        factor - probability that a single copy gene is actually non-orthologous relative to the proportion of multi-copy species
        pMax - maximum allowed probability that one of the genes is non-orthologous
    Returns:
        maxMulticopy - array, i-th element is maximum number of allowed multicopy species if ispecies are non-single copy
    """        
    allowed_multicopy = []
    for i_non_single in range(n_species):
        n_single = n_species - i_non_single
        q_found_max = False
        for n_multiple in range(i_non_single+1):
            p_false = factor*n_multiple/float(n_species)
            p_any_false = 1.-(1-p_false)**n_single
            if p_any_false > p_max:
                allowed_multicopy.append(n_multiple - 1)
                q_found_max = True
                break
        if not q_found_max:
            allowed_multicopy.append(i_non_single)
    return allowed_multicopy  


def single_copy_with_probability_test(fraction, og_matrix):
    """ X fraction of species have a single copy of the gene, all other species can have whatever number of copies including 0"""
    n_species = og_matrix.shape[1]
    allowed_multicopy = get_multicopy_cutoff(n_species)
    multi = (og_matrix > 1 * np.ones((1, n_species))).sum(1) # multicopy
    excluded = n_species - (og_matrix == 1 * np.ones((1, n_species))).sum(1) # single-copy
    temp2 = (og_matrix == np.ones((1, n_species))).sum(1) / float(n_species)   # proportion of species with gene present
    c2 = temp2 > fraction 
    pass_frac_ogs = c2.nonzero()[0] 
    ok_ogs = [iog for iog in pass_frac_ogs if multi[iog] <= allowed_multicopy[excluded[iog]]]
    return ok_ogs 


def get_orthogroup_occupancy_info(m):
    """
    Args:
        m - orthogroup matrix
    """
    if getattr(m, "ndim", None) is None or m.ndim < 2:
        raise ValueError("m must be a 2D array/matrix")
    N = m.shape[1]
    if N == 0:
        raise ValueError("Matrix has zero columns; cannot compute occupancy.")
    eps = 1e-5
    fractions = [1.0 - n / N for n in range(N)]              # length N
    n_ortho = [float(len(single_copy_with_probability_test(F - eps, m)))
              for F in fractions]

    return fractions, n_ortho


def determine_orthogroups_for_species_tree(m, iogs4, n_ogs_min=100, n_sufficient=1000, increase_required=2.):
    """Orthogroups can be used if at least a fraction f of the species in the orthogroup are single copy, f is determined as described 
    below. Species that aren't single copy are allowed in the orthogroup as long as there aren't too many (criterion is quite strict). 
    The presence of species with multiple copies suggests that the single copy species might actually have arisen from duplication 
    and loss. Therefore, we use a probability model to determine how many of the excluded species can be multicopy before it is likely
    that the single copy genes are actually hidden paralogues.
    Args:
        m - the orthogroup matrix shape=(nOGs, nSpecies), each entry is number of genes from that species
        nOGsMin - the minimum number of orthogroups. Need to have at least this many before considering the next criteria
        p - proportionalIncreaseMin: Each decrease in the proportion of single copy species should bring about a relative increase 
        in the number of orthogroups that will be used of 'p', otherwise it is not worth lowering the bar
    """
    fractions, n_ortho = get_orthogroup_occupancy_info(m)
    if n_ortho[0] > n_sufficient:
        ogs_to_use = single_copy_with_probability_test(1.0-1e-5, m)
        return ogs_to_use, 1.0
    n_species = m.shape[1]
    for i in range(1, n_species):
        if n_ortho[i-1] < n_ogs_min: continue
        if n_ortho[i-1] > n_sufficient: break
        p = ((n_ortho[i] - n_ortho[i-1])/n_ortho[i-1]) / (-(fractions[i] - fractions[i-1])/fractions[i-1])
        if fractions[i] > 0.5:
            if p < increase_required: break
        else:
            # if fewer than half the species are in any orthogroup then set higher bar for reducing fraction further
            if p < 2*increase_required: break
    f = fractions[i-1]
    ogs_to_use = single_copy_with_probability_test(f-1e-5, m)
    ogs_to_use = [iogs4[i_] for i_ in ogs_to_use]
    return ogs_to_use, f   


class MSA(object):
    def __init__(self, msa_dict):
        self.seqs = msa_dict
        self.length = len(list(msa_dict.values())[0])


def read_alignment(fn):
    msa = dict()
    accession = None
    length = None
    seq = ""
    with open(fn, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if accession != None:
                    if length != None and len(seq) != length:
                        text  = "ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                        file_handler.handler.log_fail_and_exit(text)
                    msa[accession] = seq
                accession = line[1:]
                seq = ""
            else:
                seq += line
        if accession != None:
            if length != None and len(seq) != length:
                text = "Error: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                file_handler.handler.log_fail_and_exit(text)
            msa[accession] = seq
    return MSA(msa)


def orthogroup_matrix(ogs, min_seq):
    """
    (og_matrix, iogs4): the number of genes of each species (columns) in each
    orthogroup with at least min_seq genes (rows, the orthogroups iogs4).
    ogs: the orthogroups, lists of genes (with i_sp).
    """
    iogs4 = [i for i, og in enumerate(ogs) if len(og) >= min_seq]
    ogs = [ogs[i] for i in iogs4]
    i_species = sorted(set([gene.i_sp for og in ogs for gene in og]))
    species_index_dict = {i_sp: i_col for i_col, i_sp in enumerate(i_species)}
    # (i, j)-th entry of ogMatrix gives the number of genes from i in orthologous group j
    og_matrix = np.zeros((len(ogs), len(i_species)))
    for i_og, og in enumerate(ogs):
        for gene in og:
            og_matrix[i_og, species_index_dict[gene.i_sp]] += 1
    return og_matrix, iogs4


def orthogroup_pieces(fn, selected_seqs, all_species):
    """
    One orthogroup's piece of each species' row of the concatenated alignment.
    fn: its alignment; selected_seqs: the IDs of its single-copy genes;
    all_species: the IDs of all species.
    Returns ("ok", {species ID: aligned sequence}), or ("different lengths", None)
    or ("failed", None) if the alignment cannot be used.
    """
    try:
        alignment = read_alignment(fn)
        # This orthogroup's piece of each species' row, added only once the
        # whole alignment is known to be usable, so the rows stay in step.
        pieces = {}
        for name, al in alignment.seqs.items():
            if name.split()[0] in selected_seqs:
                i_sp = name.split("_")[0]
                pieces[i_sp] = al.replace('*', '-')  # this allows for the MSA method to have failed to put the sequence in the MSA
        if any(len(al) != alignment.length for al in pieces.values()):
            # A malformed MSA (rows of different lengths) would shift the
            # columns of every species after it.
            return "different lengths", None
        # now put blanks for the missing species, in a fixed order (set
        # order differs between runs, and sets the rows' order)
        for i_sp in sorted(all_species.difference(pieces), key=int):
            pieces[i_sp] = "-"*(alignment.length)
        return "ok", pieces
    except IndexError:
        # allow empty MSA (could fail for unknown reason)
        return "failed", None


def create_concatenated_alignment(ogs_to_use_ids, ogs, alignment_filename_function, output_filename, f_single_copy, f_max_gap=0.5,
                                  q_progress=False, species_names=None):
    """
    q_progress: show a progress bar
    species_names: {species ID: name} to name the rows (by default their species IDs)
    """
    all_species = {str(gene.i_sp) for og in ogs for gene in og}
    if q_progress:
        progressbar, task = printing.get_progressbar(len(ogs_to_use_ids))
        progressbar.start()
    results = []
    try:
        for i_og in ogs_to_use_ids:
            species_counts = Counter([gene.i_sp for gene in ogs[i_og]])
            selected_seqs = {gene.to_string() for gene in ogs[i_og] if species_counts[gene.i_sp] == 1}
            fn = alignment_filename_function(i_og)
            results.append((fn, orthogroup_pieces(fn, selected_seqs, all_species)))
            if q_progress:
                progressbar.update(task, advance=1)
    finally:
        if q_progress:
            progressbar.stop()
    # Each species' alignment is collected as a list of pieces and joined once:
    # appending to a growing string copies it every time.
    concatentaed_alignments = defaultdict(list)
    for fn, (status, pieces) in results:
        if status == "different lengths":
            print("WARNING: The sequences of an MSA have different lengths, it is not used for the species tree: %s" % fn)
        elif status == "failed":
            print("WARNING: An MSA failed for an unknown reason: %s" % fn)
            print("No tree or orthologs will be inferred for this orthogroup. To correct the issue, identify & correct the problematic gene sequence and rerun.")
        else:
            for i_sp, al in pieces.items():
                concatentaed_alignments[i_sp].append(al)
    concatentaed_alignments = {i_sp: "".join(parts) for i_sp, parts in concatentaed_alignments.items()}
    # Trim the completed alignment: to 50% of fraction of species present
    max_gap = (1.-f_max_gap*f_single_copy)*len(all_species)
    if len(concatentaed_alignments) == 0 or len(list(concatentaed_alignments.values())[0]) == 0:
        print("All MSAs for the concatenated multiple sequence alignment were empty.")
        print("Please correct the error and re-run.")
        util.fail()
    names = list(concatentaed_alignments.keys())
    row_name = (lambda i_sp: i_sp) if species_names is None else (lambda i_sp: species_names[i_sp])
    n_char = 80
    try:
        # one byte per position (a character array takes 4, and building it
        # from lists of characters 8 more): species x columns can be large
        length = len(concatentaed_alignments[names[0]])
        M = np.frombuffer("".join(concatentaed_alignments[name] for name in names).encode("latin-1"),
                          dtype=np.uint8).reshape(len(names), length)
    except (UnicodeEncodeError, ValueError):
        M = None      # unusual letters (rows all have the same length, see above): as before, below
    if M is not None:
        M = M[:, (M == ord("-")).sum(axis=0) <= max_gap]
        with open(output_filename, 'w') as outfile:
            for i_seq, name in enumerate(names):
                outfile.write(">%s\n" % row_name(name))
                seq = M[i_seq].tobytes().decode("latin-1")
                for i in range(0, len(seq), n_char):
                    outfile.write(seq[i:i+n_char] + "\n")
        return
    M = np.array([list(concatentaed_alignments[name]) for name in names])
    gap_counts = sum(M == "-")
    i_keep = np.where(gap_counts <= max_gap)
    M = M[:, i_keep]
    with open(output_filename, 'w') as outfile:
        for i_seq, name in enumerate(names):
            outfile.write(">%s\n" % row_name(name))
            seq = M[i_seq,:].tolist()[0]
            for i in range(0, len(seq), n_char):
                outfile.write("".join(seq[i:i+n_char]) + "\n")
