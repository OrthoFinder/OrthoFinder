#!/usr/bin/env python3

import os
import shutil
# import argparse
import numpy as np
try:
    from rich import print
except ImportError:
    ...

class MSA(object):
    """
    An alignment as one array: a row of characters per sequence (np.frombuffer
    of the sequences' bytes), so the columns' counts and the trimming are
    single passes in numpy rather than loops over the characters.
    Gaps are "-" and "*" (as in the counts; "*" is written as "-").
    """
    def __init__(self, fn):
        self.names = []         # list of names
        seqs = []
        with open(fn, 'r', encoding="utf-8") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith(">"):
                    self.names.append(line[1:])
                    seqs.append([])
                else:
                    seqs[-1].append(line)
        seqs = ["".join(parts) for parts in seqs]
        self.n = len(self.names)    # number of sequences
        self.length = len(seqs[0]) if seqs else 0   # number of columns in MSA
        for i, seq in enumerate(seqs[1:], 2):
            if len(seq) != self.length:
                # (an error, not sys.exit(): that is a successful exit, and the run would
                # end as if stopped normally, without its results)
                raise ValueError("ERROR: Malformed alignment %s: sequence %d has length %d, previous was %d"
                                 % (fn, i, len(seq), self.length))
        # one byte a character; other characters (not Latin-1) as 4-byte code points
        try:
            self.encoding, dtype = "latin-1", np.uint8
            data = "".join(seqs).encode(self.encoding)
        except UnicodeEncodeError:
            self.encoding, dtype = "utf-32-le", np.uint32
            data = "".join(seqs).encode(self.encoding)
        self.M = np.frombuffer(data, dtype=dtype).reshape(self.n, self.length)
        self.residues = (self.M != ord("-")) & (self.M != ord("*"))   # non-gap characters

    def write_msa(self, i_cols, outfn, n_char = 80):
        kept = np.where(self.residues[:, i_cols], self.M[:, i_cols], ord("-")).astype(self.M.dtype)
        with open(outfn, 'w', encoding="utf-8", newline="\n") as outfile:
            for name, row in zip(self.names, kept):
                seq = row.tobytes().decode(self.encoding)
                outfile.write(">" + name + "\n")
                outfile.write("".join(seq[i:i+n_char] + "\n" for i in range(0, len(seq), n_char)))

def main(infn, outfn, msa_min_seq=10, f=0.1, n_min=500, c=0.75, exe=False):
    # if exe:
    #     run_exe(infn, outfn, f, n_min, c)
    # else:
    run_in_process(infn, outfn, msa_min_seq, f, n_min, c)

# def run_exe(infn, outfn, f, n_min, c):
#     """run trim in a separate process not forked from current process"""
#     subprocess.run(["python", __file__, "-c", str(c), infn, outfn, str(f), str(n_min)])


def run_in_process(infn, outfn, msa_min_seq, f, n_min, c):
    """
    Trim the alignment file. Auto reduce f parameter if required. OrthoFinder default
    parameters are used by default.
    Args:
        infn - input FASTA filename
        outfn - output FASTA filename
        f - Lower limit for fraction of meaningful characters
        n_min - The min trimmed length
        c - The min fraction of non-gap characters to conserve
    Notes:
        If f is too high then it will reduce f s.t. len >= n_min or the full 
        alignment is kept (if less than n_min).
        In some way fraction of gaps is unsatisfactory since as the number of sequences
        in the MSA gets larger a column could contain a high fraction of gaps and 
        yet still contain useful phylogenetic signal for (many) sequences that are
        closely related and don't have a gap in that column. For this reason the 
        c parameter is included which specifies the minimum total fraction of the 
        original non-gap characters that must be preserved.
        An alignment is not trimmed (it is copied) if a sequence would be left
        with no residues.
    """
    if not os.path.exists(infn):
        raise FileNotFoundError("ERROR: The alignment to trim does not exist: %s" % infn)
    try:
        # BEWARE: Input and output files may be the same
        with open(outfn, 'a', encoding="utf-8", newline="\n") as outfile:
            outfile.write(" ")
    except OSError as e:
        raise OSError("ERROR: Cannot write the trimmed alignment %s: %s" % (outfn, e)) from e

    msa = MSA(infn)   
    n = msa.n
    if msa.length <= n_min or n <= msa_min_seq:
        copy_input_to_output(infn, outfn)
        # print("Already below min length: %d" % msa.length)
        return

    length = msa.length
    max_gap = (1.-f)*n
    aa_counts = msa.residues.sum(axis=0)   # per column
    gap_counts = n - aa_counts
    aa_before = int(aa_counts.sum())
    i_keep = np.where(gap_counts <= max_gap)
    n_keep = i_keep[0].size     # it's an I, J tuple
    aa_after = aa_counts[i_keep].sum()
    # Check we have kept enough columns and characters
    if n_keep < n_min or aa_after < c*aa_before:
        # print("%0.3f: %d columns, %0.1f%% of characters retained" % (f, n_keep, 100.*aa_after/aa_before))
        f, i_keep = get_satifactory_f(gap_counts, aa_counts, n, f, n_min, c)

    to_keep = i_keep[0]
    # trimmed only if every sequence keeps a residue
    trim = bool(msa.residues[:, to_keep].any(axis=1).all())
    n_keep = to_keep.size
    if n_keep == length or not trim:
        copy_input_to_output(infn, outfn)
    if trim:
        msa.write_msa(to_keep, outfn)
    # print("%0.3f: %d->%d, %0.1f%% characters retained. Trimmed %s" % (f, length, to_keep.size, 100.*aa_after/aa_before, infn))

def get_satifactory_f(gap_counts, aa_counts, N, f_orig, n_min, c, tol = 0.001):
    """
    The f used was too large, get an f that gives n_min columns
    Args:
        gap_counts = np array of gaps per column
        aa_counts = np array of amino acids per column
        N - number of sequences
        f_orig - starting f for the search
        n_min - minimum alignment length
        c - minimum fraction of amino acids to keep
        tol - convergence tolerance
    Implementation:
        down - lower bound on non-gap percentage, giving larger alignment
        up - upper bound on non-gap character percentage
    """
    # binary search 
    aa_before = aa_counts.sum()
    up = f_orig
    down = 0.
    x_prev = f_orig
    x = 0.5*(down+up)
    while abs(x-x_prev) > tol:
        max_gap = (1.-x)*N
        i_keep = np.where(gap_counts <= max_gap)
        n_keep = i_keep[0].size
        aa_new = aa_counts[i_keep].sum()
        # print("%0.3f: %d columns, %0.1f%% of characters retained" % (x, n_keep, 100.*aa_new/aa_before))
        if n_keep < n_min or aa_new < c*aa_before:
            up = x
        elif (n_keep == n_min and aa_new >= c*aa_before) or (n_keep >= n_min and aa_new == c*aa_before):
            # one on the limit, the other satisfied - stop exactly here (x is used)
            down = x
            break
        else:
            down = x
        x_prev = x
        x = 0.5*(down+up)
    # Use the value which gave just above the min length
    x = down
    max_gap = (1.-x)*N
    i_keep = np.where(gap_counts <= max_gap)
    aa_new = aa_counts[i_keep].sum()
    # print("%0.3f: %d columns, %0.1f%% characters retained" % (x, n_keep, 100.*aa_new/aa_before))
    return x, i_keep


def copy_input_to_output(infn, outfn):      
    # if the outfn doesn't exist doesn't the check for being the same is redundant  
    if (not os.path.exists(outfn)) or (not os.path.samefile(infn, outfn)):
        shutil.copy(infn, outfn)

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("infn", help="Input alignment fasta file")
#     parser.add_argument("outfn", help="Output alignment fasta file")
#     parser.add_argument("gt", 
#                         type=float,
#                         help="Lower limit for fraction of non-gap characters")
#     parser.add_argument("min", type=int, help="Minimum number of columns to keep") 
#     parser.add_argument("-c", "--conserve", type=float, 
#                         default=0.0,
#                         help="Conserve at least this fraction of non-gap characters") 
#     args = parser.parse_args()
#     main(args.infn, args.outfn, args.gt, args.min, args.conserve)