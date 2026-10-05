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
import gzip
import array
import numpy as np
from scipy import sparse
try:
    from rich import print
except ImportError:
    ...
from . import util, file_io

PY2 = sys.version_info <= (3,)       
file_read_mode = 'rb' if PY2 else 'rt'

def GetBLAST6Scores(seqsInfo, blastDir_list, iSpecies, jSpecies, qExcludeSelfHits = True, sep = "_", qDoubleBlast=True, q_allow_empty=False, fmt="lil"):
    """
    Bit-score matrix (nSeqs_i x nSeqs_j) for the hits of species iSpecies
    against jSpecies: the highest score for each (query, hit) pair; hits with
    a score <= 0 are ignored. fmt is "lil" (default) or "csr".

    Hits are collected into arrays and the matrix is built once, rather than
    by setting one lil_matrix element per hit, which is several times slower.
    """
    qSameSpecies = iSpecies==jSpecies
    qCheckForSelfHits = qExcludeSelfHits and qSameSpecies
    if not qDoubleBlast:
        qRev = (iSpecies > jSpecies)
    else:
        qRev = False      
    if qRev:
        iSpeciesOpen = jSpecies
        jSpeciesOpen = iSpecies
    else:        
        iSpeciesOpen = iSpecies
        jSpeciesOpen = jSpecies
    nSeqs_i = seqsInfo.nSeqsPerSpecies[iSpecies]
    nSeqs_j = seqsInfo.nSeqsPerSpecies[jSpecies]

    def result(I, J, S):
        B = sparse.csr_matrix((S, (I, J)), shape=(nSeqs_i, nSeqs_j))
        return B.tolil() if fmt == "lil" else B

    empty = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
    row = ""
    for d in blastDir_list:
        if d[-1] != os.sep:
            d += os.sep
        fn = d + "Blast%d_%d.txt" % (iSpeciesOpen, jSpeciesOpen)
        if os.path.exists(fn) or os.path.exists(fn + ".gz"): break
    if q_allow_empty and not os.path.exists(fn) and not os.path.exists(fn + ".gz"):
        return result(*empty)
    path = fn + ".gz" if os.path.exists(fn + ".gz") else fn
    # The file's columns depend on the search program and its options (as
    # recorded when the search was run): read the IDs and bit score by name.
    layout = file_io.hit_format(path)
    # Fast path for large files.
    query, hit = ("subject_seq", "query_seq") if qRev else ("query_seq", "subject_seq")
    cols = None
    if file_io.use_pandas(path):
        try:
            cols = file_io.read_hit_columns(path, (query, hit, "bitscore"), fmt=layout)
        except (ValueError, OSError):
            cols = None   # the line-by-line reader below reports any problem
    if cols is not None:
        return _bit_score_matrix(cols[query], cols[hit], cols["bitscore"], nSeqs_i, nSeqs_j,
                                 qCheckForSelfHits, iSpecies, jSpecies, result, empty)

    I = array.array("q")
    J = array.array("q")
    S = array.array("d")
    try:
        iQ, iH = layout.index("qseqid"), layout.index("sseqid")
        if qRev:
            iQ, iH = iH, iQ
        iBitscore = layout.index("bitscore")
        with (gzip.open(fn + ".gz", file_read_mode) if os.path.exists(fn + ".gz") else open(fn, file_read_mode)) as blastfile:
            # Search output is never quoted, so a plain split gives the fields
            for line in file_io.hit_lines(blastfile, layout):
                row = line.rstrip("\r\n").split("\t")
                # Get hit and query IDs
                try:
                    sequence1ID = int(row[iQ].split(sep, 2)[1])
                    sequence2ID = int(row[iH].split(sep, 2)[1])
                except (IndexError, ValueError):
                    sys.stderr.write("\nERROR: Query or hit sequence ID in BLAST results file was missing or incorrectly formatted.\n")
                    raise
                # Get bit score for pair
                try:
                    score = float(row[iBitscore])
                except (IndexError, ValueError):
                    sys.stderr.write("\nERROR: field %d of the search results should be the bit score of the hit "
                                     "(columns: %s)\n" % (iBitscore + 1, " ".join(layout.fields)))
                    raise
                I.append(sequence1ID)
                J.append(sequence2ID)
                S.append(score)
    except Exception:
        print("ERROR: Blast%d_%d.txt is corrupted" % (iSpecies, jSpecies))
        sys.stderr.write("Malformatted line in %sBlast%d_%d.txt\nOffending line was:\n" % (d, iSpecies, jSpecies))
        sys.stderr.write("\t".join(row) + "\n")
        raise

    return _bit_score_matrix(
        np.frombuffer(I, dtype=np.int64), np.frombuffer(J, dtype=np.int64),
        np.frombuffer(S, dtype=np.float64), nSeqs_i, nSeqs_j,
        qCheckForSelfHits, iSpecies, jSpecies, result, empty)


def _bit_score_matrix(I, J, S, nSeqs_i, nSeqs_j, qCheckForSelfHits, iSpecies, jSpecies, result, empty):
    """Max score per (query, hit) as a matrix; scores <= 0 and (optionally) self-hits dropped."""
    keep = S > 0     # the old element-wise update only stored scores above 0
    if qCheckForSelfHits:
        keep &= I != J
    I, J, S = I[keep], J[keep], S[keep]
    if I.size == 0:
        return result(*empty)

    bad = (I < 0) | (I >= nSeqs_i) | (J < 0) | (J >= nSeqs_j)
    if bad.any():
        k = int(np.argmax(bad))
        sequence1ID, sequence2ID = int(I[k]), int(J[k])
        def ord(n):
            return str(n)+("th" if 4<=n%100<=20 else {1:"st",2:"nd",3:"rd"}.get(n%10, "th"))
        sys.stderr.write("\nERROR: Inconsistent input files.\n")
        kSpecies, nSeqs_k, sequencekID = (iSpecies,  nSeqs_i, sequence1ID) if not (0 <= sequence1ID < nSeqs_i) else (jSpecies,  nSeqs_j, sequence2ID)
        text = ("Blast%d_%d.txt is corrupted: Species%d.fa contains only %d sequences but found a query/hit "
                "for sequence %d_%d (i.e. %s sequence in species %d)."
                % (iSpecies, jSpecies, kSpecies, nSeqs_k, kSpecies, sequencekID, ord(sequencekID+1), kSpecies))
        print("ERROR: " + text)
        raise ValueError(text)

    # Keep the highest score for each (query, hit) pair: sort by pair, then
    # score, and take the last entry of each pair.
    key = I * nSeqs_j + J
    order = np.lexsort((S, key))
    key = key[order]
    last = np.ones(key.size, dtype=bool)
    last[:-1] = key[1:] != key[:-1]
    key = key[last]
    return result(key // nSeqs_j, key % nSeqs_j, S[order][last])
