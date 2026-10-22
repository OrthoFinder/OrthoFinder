# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com  

import os
import sys
import array
import numpy as np
from scipy import sparse
try:
    from rich import print
except ImportError:
    ...
from . import util, file_io, matrices

def get_blast6_scores(seqs_info, blast_dir_list, i_species, j_species, q_exclude_self_hits = True, sep = "_", q_double_blast=True, q_allow_empty=False, fmt="lil"):
    """
    Bit-score matrix (nSeqs_i x nSeqs_j) for the hits of species iSpecies
    against jSpecies: the highest score for each (query, hit) pair; hits with
    a score <= 0 are ignored. fmt is "lil" (default) or "csr".

    Hits are collected into arrays and the matrix is built once, rather than
    by setting one lil_matrix element per hit, which is several times slower.
    """
    q_same_species = i_species==j_species
    q_check_for_self_hits = q_exclude_self_hits and q_same_species
    if not q_double_blast:
        q_rev = (i_species > j_species)
    else:
        q_rev = False      
    if q_rev:
        i_species_open = j_species
        j_species_open = i_species
    else:        
        i_species_open = i_species
        j_species_open = j_species
    n_seqs_i = seqs_info.n_seqs_per_species[i_species]
    n_seqs_j = seqs_info.n_seqs_per_species[j_species]

    def result(I, J, S):
        # (in order of (I, J): _bit_score_matrix)
        B = matrices.csr_from_sorted(S, I, J, (n_seqs_i, n_seqs_j))
        return B.tolil() if fmt == "lil" else B

    empty = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
    row = ""
    # The results file, compressed or not, in the first directory that has it: looked
    # up once (one file per pair of species; on a network filesystem each lookup is slow)
    path = None
    for d in blast_dir_list:
        fn = os.path.join(d, "Blast%d_%d.txt" % (i_species_open, j_species_open))
        if os.path.exists(fn + ".gz"):
            path = fn + ".gz"
        elif os.path.exists(fn):
            path = fn
        if path is not None:
            break
    if path is None:
        if q_allow_empty:
            return result(*empty)
        path = fn       # (reading it reports the missing file)
    # The file's columns depend on the search program and its options (as
    # recorded when the search was run): read the IDs and bit score by name.
    layout = file_io.hit_format(path)
    # Fast path for large files.
    query, hit = ("subject_seq", "query_seq") if q_rev else ("query_seq", "subject_seq")
    cols = None
    if file_io.use_pandas(path):
        try:
            cols = file_io.read_hit_columns(path, (query, hit, "bitscore"), fmt=layout)
        except (ValueError, OSError):
            cols = None   # the line-by-line reader below reports any problem
    if cols is not None:
        return _bit_score_matrix(cols[query], cols[hit], cols["bitscore"], n_seqs_i, n_seqs_j,
                                 q_check_for_self_hits, i_species, j_species, result, empty)

    I = array.array("q")
    J = array.array("q")
    S = array.array("d")
    try:
        i_q, i_h = layout.index("qseqid"), layout.index("sseqid")
        if q_rev:
            i_q, i_h = i_h, i_q
        i_bitscore = layout.index("bitscore")
        with file_io.open_text(path) as blastfile:
            # Search output is never quoted, so a plain split gives the fields
            for line in file_io.hit_lines(blastfile, layout):
                row = line.rstrip("\r\n").split("\t")
                # Get hit and query IDs
                try:
                    sequence1_id = int(row[i_q].split(sep, 2)[1])
                    sequence2_id = int(row[i_h].split(sep, 2)[1])
                except (IndexError, ValueError):
                    sys.stderr.write("\nERROR: Query or hit sequence ID in BLAST results file was missing or incorrectly formatted.\n")
                    raise
                # Get bit score for pair
                try:
                    score = float(row[i_bitscore])
                except (IndexError, ValueError):
                    sys.stderr.write("\nERROR: field %d of the search results should be the bit score of the hit "
                                     "(columns: %s)\n" % (i_bitscore + 1, " ".join(layout.fields)))
                    raise
                I.append(sequence1_id)
                J.append(sequence2_id)
                S.append(score)
    except Exception:
        print("ERROR: Blast%d_%d.txt is corrupted" % (i_species, j_species))
        sys.stderr.write("Malformatted line in %sBlast%d_%d.txt\nOffending line was:\n" % (d, i_species, j_species))
        sys.stderr.write("\t".join(row) + "\n")
        raise

    return _bit_score_matrix(
        np.frombuffer(I, dtype=np.int64), np.frombuffer(J, dtype=np.int64),
        np.frombuffer(S, dtype=np.float64), n_seqs_i, n_seqs_j,
        q_check_for_self_hits, i_species, j_species, result, empty)


def _bit_score_matrix(I, J, S, n_seqs_i, n_seqs_j, q_check_for_self_hits, i_species, j_species, result, empty):
    """Max score per (query, hit) as a matrix; scores <= 0 and (optionally) self-hits dropped."""
    keep = S > 0     # the old element-wise update only stored scores above 0
    if q_check_for_self_hits:
        keep &= I != J
    I, J, S = I[keep], J[keep], S[keep]
    if I.size == 0:
        return result(*empty)

    bad = (I < 0) | (I >= n_seqs_i) | (J < 0) | (J >= n_seqs_j)
    if bad.any():
        k = int(np.argmax(bad))
        sequence1_id, sequence2_id = int(I[k]), int(J[k])
        def ord(n):
            return str(n)+("th" if 4<=n%100<=20 else {1:"st",2:"nd",3:"rd"}.get(n%10, "th"))
        sys.stderr.write("\nERROR: Inconsistent input files.\n")
        k_species, n_seqs_k, sequencek_id = (i_species,  n_seqs_i, sequence1_id) if not (0 <= sequence1_id < n_seqs_i) else (j_species,  n_seqs_j, sequence2_id)
        text = ("Blast%d_%d.txt is corrupted: Species%d.fa contains only %d sequences but found a query/hit "
                "for sequence %d_%d (i.e. %s sequence in species %d)."
                % (i_species, j_species, k_species, n_seqs_k, k_species, sequencek_id, ord(sequencek_id+1), k_species))
        print("ERROR: " + text)
        raise ValueError(text)

    # Keep the highest score for each (query, hit) pair: sort by pair, then
    # score, and take the last entry of each pair.
    key = I * n_seqs_j + J
    order = np.lexsort((S, key))
    key = key[order]
    last = np.ones(key.size, dtype=bool)
    last[:-1] = key[1:] != key[:-1]
    key = key[last]
    return result(key // n_seqs_j, key % n_seqs_j, S[order][last])
