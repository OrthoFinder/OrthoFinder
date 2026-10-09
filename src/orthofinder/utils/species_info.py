"""The species and sequences of an analysis."""

import os
from collections import namedtuple


# Get Info from seqs IDs file?
def count_fasta_sequences(fasta_filename, block_size=1 << 24):
    """
    Number of lines starting with ">" in a FASTA file.

    Counts the byte pattern "<newline>>" in large blocks (C-level bytes.count)
    instead of looping over every line in Python; line ends "\n", "\r\n" and
    "\r" are all recognised, as with Python's text-mode line splitting.
    """
    n = 0
    previous = b"\n"            # the start of the file counts as a line start
    with open(fasta_filename, "rb") as infile:
        while True:
            block = infile.read(block_size)
            if not block:
                return n
            n += block.count(b"\n>") + block.count(b"\r>")
            if block[:1] == b">" and previous in (b"\n", b"\r"):
                n += 1             # a header right at the start of this block
            previous = block[-1:]


SequencesInfo = namedtuple(
    "SequencesInfo", "n_seqs n_species species_to_use seq_starting_indices n_seqs_per_species"
)


"""
n_seqs: int - number of sequences in species_to_use
n_species: int - number in species_to_use
species_to_use: List[int]
seq_starting_indices: List[int] of size |species_to_use|, reindexed from zero according to the speices present
n_seqs_per_species: Dict[int, int] - indexed by OrthoFinder species ID, info on all species, nto just those included in analysis
"""


def get_seqs_info(input_directory_list, species_to_use, n_sp_all):
    input_directory_list = list(input_directory_list)
    species_to_use = list(species_to_use)
    seq_starting_indices = [0]
    n_seqs = 0
    n_seqs_per_species = {}

    for i_fasta in range(n_sp_all):
        fasta_filename = None

        for d in input_directory_list:
            candidate = os.path.join(d, f"Species{i_fasta}.fa")
            if os.path.exists(candidate):
                fasta_filename = candidate
                break

        if fasta_filename is None:
            if i_fasta in species_to_use:
                raise FileNotFoundError(
                    f"Species{i_fasta}.fa not found in any of the input "
                    f"directories: {input_directory_list}"
                )
            continue
        n = count_fasta_sequences(fasta_filename)

        n_seqs_per_species[i_fasta] = n

        if i_fasta in species_to_use:
            n_seqs += n
            seq_starting_indices.append(n_seqs)

    if seq_starting_indices:
        seq_starting_indices = seq_starting_indices[:-1]

    n_species = len(species_to_use)

    return SequencesInfo(
        n_seqs=n_seqs,
        n_species=n_species,
        species_to_use=species_to_use,
        seq_starting_indices=seq_starting_indices,
        n_seqs_per_species=n_seqs_per_species,
    )


def seqs_info_recompute(seqs_info_orig, new_species_to_use):
    """
    Args:
        seqs_info_orig
        new_species_to_use: List[int]
    Recomputes SeqsInfo for a subset of species:
    - speciesToUse updated
    - seqStartingIndices recomputed
    - nSeqs updated
    - nSpecies updated
    """
    seq_starting_indices = [0]
    for isp in new_species_to_use:
        seq_starting_indices.append(
            seq_starting_indices[-1] + seqs_info_orig.n_seqs_per_species[isp]
        )
    n_seqs = seq_starting_indices[-1]
    seq_starting_indices = seq_starting_indices[:-1]
    return SequencesInfo(
        n_seqs=n_seqs,
        n_species=len(new_species_to_use),
        species_to_use=new_species_to_use,
        seq_starting_indices=seq_starting_indices,
        n_seqs_per_species=seqs_info_orig.n_seqs_per_species,
    )


class SpeciesInfo(object):
    def __init__(self):
        self.species_to_use = []  #       seqsInfo.iSpeciesToUse   - which to include for this analysis
        self.n_sp_all = None  #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.i_first_new_species = None  #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices

    def __str__(self):
        return str((self.species_to_use, self.n_sp_all, self.i_first_new_species))

    def get_original_species(self):
        if self.i_first_new_species is None:
            return self.species_to_use
        else:
            return [i_sp for i_sp in self.species_to_use if i_sp < self.i_first_new_species]


def get_species_to_use(species_ids_fn):
    """Returns species indices (int) to use and total number of species available"""
    species_to_use = []
    species_to_use_names = []
    n_skipped = 0
    with open(species_ids_fn, "r") as species_f:
        for line in species_f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                n_skipped += 1
            else:
                i_sp, sp_name = line.split(": ")
                if sp_name not in species_to_use_names:
                    species_to_use.append(int(i_sp))
                    species_to_use_names.append(sp_name)

    return species_to_use, len(species_to_use) + n_skipped, species_to_use_names


def seq_id_key(seq_id):
    """Sort key for "species_sequence" IDs in numeric order (0_61 before 0_317)."""
    i_sp, i_seq = seq_id.split("_")
    return int(i_sp), int(i_seq)
