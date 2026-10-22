"""Choosing orthogroups for the species tree and concatenating their alignments."""

import functools
from collections import Counter, defaultdict

import numpy as np

from ..utils import printing, util

try:
    from rich import print
except ImportError:
    ...


def get_multicopy_cutoff(n_species, factor = 0.25, p_max = 0.25):
    """(a list: see _multicopy_cutoff)"""
    return list(_multicopy_cutoff(n_species, factor, p_max))


@functools.lru_cache(maxsize=None)
def _multicopy_cutoff(n_species, factor = 0.25, p_max = 0.25):
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
    return tuple(allowed_multicopy)   # (cached: species^2 steps, the same for every call)  


class Occupancy(object):
    """
    For each orthogroup (of the species tree candidates), the number of species
    with one gene (n_single) and with several (multi), of n_species: all that
    choosing them needs, without an orthogroups x species matrix (GBs with 1000s
    of species).
    """

    def __init__(self, n_single, multi, n_species):
        self.n_single = np.asarray(n_single)
        self.multi = np.asarray(multi)
        self.n_species = n_species

    @classmethod
    def of(cls, m):
        """An Occupancy, or one of an orthogroup matrix (orthogroups x species gene counts)."""
        if isinstance(m, cls):
            return m
        if getattr(m, "ndim", None) is None or m.ndim < 2:
            raise ValueError("m must be a 2D array/matrix")
        return cls((m == 1).sum(1), (m > 1).sum(1), m.shape[1])


def single_copy_with_probability_test(fraction, occupancy):
    """ X fraction of species have a single copy of the gene, all other species can have whatever number of copies including 0"""
    occupancy = Occupancy.of(occupancy)
    n_species = occupancy.n_species
    allowed_multicopy = get_multicopy_cutoff(n_species)
    multi = occupancy.multi # multicopy
    excluded = n_species - occupancy.n_single # single-copy
    temp2 = occupancy.n_single / float(n_species)   # proportion of species with gene present
    c2 = temp2 > fraction 
    pass_frac_ogs = c2.nonzero()[0] 
    ok_ogs = [iog for iog in pass_frac_ogs if multi[iog] <= allowed_multicopy[excluded[iog]]]
    return ok_ogs 


def get_orthogroup_occupancy_info(m):
    """
    Args:
        m - an Occupancy, or orthogroup matrix
    """
    occupancy = Occupancy.of(m)
    N = occupancy.n_species
    if N == 0:
        raise ValueError("Matrix has zero columns; cannot compute occupancy.")
    eps = 1e-5
    fractions = [1.0 - n / N for n in range(N)]              # length N
    # The counts of single_copy_with_probability_test(F - eps, m) for each F, from the
    # orthogroups' single- and multi-copy counts computed once (not a pass over the
    # whole matrix per fraction: species^2 x orthogroups), with the same comparisons
    n_single = occupancy.n_single
    # (allowed_multicopy[N]: no single-copy species, never above a fraction > 0)
    allowed_multicopy = np.asarray(get_multicopy_cutoff(N) + [0])
    ok = occupancy.multi <= allowed_multicopy[N - n_single]
    proportion_single = n_single / float(N)
    n_ortho = [float(np.count_nonzero(ok & (proportion_single > F - eps))) for F in fractions]

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
    m = Occupancy.of(m)
    n_species = m.n_species
    if n_species < 2:
        # the orthogroups large enough are of one species at most: none can be used
        return [], 1.0
    fractions, n_ortho = get_orthogroup_occupancy_info(m)
    if n_ortho[0] > n_sufficient:
        ogs_to_use = single_copy_with_probability_test(1.0-1e-5, m)
        return [iogs4[i_] for i_ in ogs_to_use], 1.0
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


class MSALengthError(ValueError):
    """The rows of an alignment have different lengths (read_alignment)."""


def read_alignment(fn):
    msa = dict()
    accession = None
    length = None
    seq = ""
    with open(fn, 'r', encoding="utf-8") as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if accession != None:
                    # (all the rows the length of the first)
                    if length is None:
                        length = len(seq)
                    elif len(seq) != length:
                        raise MSALengthError("Sequence length mismatch in MSA %s: %s & %d" % (fn, length, len(seq)))
                    msa[accession] = seq
                accession = line[1:]
                seq = ""
            else:
                seq += line
        if accession != None:
            if length is not None and len(seq) != length:
                raise MSALengthError("Sequence length mismatch in MSA %s: %s & %d" % (fn, length, len(seq)))
            msa[accession] = seq
    return MSA(msa)


def orthogroup_occupancy(ogs, min_seq):
    """
    (Occupancy, iogs4): of the orthogroups with at least min_seq genes (iogs4),
    the species counted those with genes in any of them (as the columns of an
    orthogroup matrix would be).
    ogs: the orthogroups, lists of genes (with i_sp).
    """
    iogs4 = [i for i, og in enumerate(ogs) if len(og) >= min_seq]
    n_single = np.zeros(len(iogs4), dtype=np.int64)
    multi = np.zeros(len(iogs4), dtype=np.int64)
    species = set()
    for row, i in enumerate(iogs4):
        counts = Counter(gene.i_sp for gene in ogs[i])
        species.update(counts)
        n_one = sum(1 for n in counts.values() if n == 1)
        n_single[row] = n_one
        multi[row] = len(counts) - n_one
    return Occupancy(n_single, multi, len(species)), iogs4


def orthogroup_pieces(fn, selected_seqs, all_species):
    """
    One orthogroup's piece of each species' row of the concatenated alignment.
    fn: its alignment; selected_seqs: the IDs of its single-copy genes;
    all_species: the IDs of all species.
    Returns ("ok", {species ID: aligned sequence}), or ("different lengths", None)
    or ("failed", None) if the alignment cannot be used.
    """
    try:
        try:
            alignment = read_alignment(fn)
        except MSALengthError:
            return "different lengths", None
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
    # Each species' alignment is collected as a list of pieces, each orthogroup's
    # added as it is read (not all kept first), then joined once: appending to a
    # growing string copies it every time.
    concatentaed_alignments = defaultdict(list)
    try:
        for i_og in ogs_to_use_ids:
            species_counts = Counter([gene.i_sp for gene in ogs[i_og]])
            selected_seqs = {gene.to_string() for gene in ogs[i_og] if species_counts[gene.i_sp] == 1}
            fn = alignment_filename_function(i_og)
            status, pieces = orthogroup_pieces(fn, selected_seqs, all_species)
            if status == "different lengths":
                print("WARNING: The sequences of an MSA have different lengths, it is not used for the species tree: %s" % fn)
            elif status == "failed":
                print("WARNING: An MSA failed for an unknown reason: %s" % fn)
                print("No tree or orthologs will be inferred for this orthogroup. To correct the issue, identify & correct the problematic gene sequence and rerun.")
            else:
                for i_sp, al in pieces.items():
                    concatentaed_alignments[i_sp].append(al)
            if q_progress:
                progressbar.update(task, advance=1)
    finally:
        if q_progress:
            progressbar.stop()
    for i_sp in list(concatentaed_alignments):
        # (one species at a time: its pieces are freed as its row is made)
        concatentaed_alignments[i_sp] = "".join(concatentaed_alignments[i_sp])
    concatentaed_alignments = dict(concatentaed_alignments)
    # Trim the completed alignment: to 50% of fraction of species present
    max_gap = (1.-f_max_gap*f_single_copy)*len(all_species)
    if len(concatentaed_alignments) == 0 or len(list(concatentaed_alignments.values())[0]) == 0:
        print("All MSAs for the concatenated multiple sequence alignment were empty.")
        print("Please correct the error and re-run.")
        util.fail()
    names = list(concatentaed_alignments.keys())
    row_name = (lambda i_sp: i_sp) if species_names is None else (lambda i_sp: species_names[i_sp])
    n_char = 80

    def row_bytes(name):
        # one byte per position (a character array takes 4): one species' row at a time,
        # not an array of them all (species x columns can be large)
        return np.frombuffer(concatentaed_alignments[name].encode("latin-1"), dtype=np.uint8)

    length = len(concatentaed_alignments[names[0]])
    try:
        if any(len(concatentaed_alignments[name]) != length for name in names):
            raise ValueError("rows of different lengths")
        gap_counts = np.zeros(length, dtype=np.int64)
        for name in names:
            gap_counts += row_bytes(name) == ord("-")
        keep = gap_counts <= max_gap
    except (UnicodeEncodeError, ValueError):
        keep = None      # unusual letters (rows all have the same length, see above): as before, below
    if keep is not None:
        with open(output_filename, 'w', encoding="utf-8", newline="\n") as outfile:
            for name in names:
                outfile.write(">%s\n" % row_name(name))
                seq = row_bytes(name)[keep].tobytes().decode("latin-1")
                for i in range(0, len(seq), n_char):
                    outfile.write(seq[i:i+n_char] + "\n")
        return
    M = np.array([list(concatentaed_alignments[name]) for name in names])
    gap_counts = sum(M == "-")
    i_keep = np.where(gap_counts <= max_gap)
    M = M[:, i_keep]
    with open(output_filename, 'w', encoding="utf-8", newline="\n") as outfile:
        for i_seq, name in enumerate(names):
            outfile.write(">%s\n" % row_name(name))
            seq = M[i_seq,:].tolist()[0]
            for i in range(0, len(seq), n_char):
                outfile.write("".join(seq[i:i+n_char]) + "\n")
