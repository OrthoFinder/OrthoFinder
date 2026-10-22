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
import glob
import numpy as np
try:
    import cPickle as pic
except ImportError:
    import pickle as pic

from .. import pic_protocol

# Matrices are stored in blocks, one per species pair (iSpecies, jSpecies).
# dump_matrix_array writes all the blocks of one species into a single "pack"
# file: a header with the byte offset of every block, followed by the pickled
# blocks. That is n files per matrix instead of n^2 (e.g. ~1,270 instead of
# ~1.6 million at 1,270 species), while reading one block is still one open,
# one seek and one unpickle. dump_matrix still writes single-block files, and
# load_matrix reads either kind.

_PACK_MAGIC = b"OFPACK1\n"


def _pack_fn(name, i_species, d_pickle):
    return d_pickle + "%s%d.pack" % (name, i_species)


def _read_pack_index(f):
    if f.read(len(_PACK_MAGIC)) != _PACK_MAGIC:
        raise ValueError("Not a matrix pack file: %s" % f.name)
    n = int(np.frombuffer(f.read(8), dtype="<i8")[0])
    return np.frombuffer(f.read(8 * (n + 1)), dtype="<i8")


def dump_matrix(name, m, i_species, j_species, d_pickle):
    with open(d_pickle + "%s%d_%d.pic" % (name, i_species, j_species), 'wb') as pic_file:
        pic.dump(m, pic_file, protocol=pic_protocol)
    
def dump_matrix_array(name, matrix_array, i_species, d_pickle):
    fn = _pack_fn(name, i_species, d_pickle)
    n = len(matrix_array)
    offsets = np.zeros(n + 1, dtype="<i8")
    with open(fn + ".tmp", "wb") as f:
        f.write(b"\0" * (len(_PACK_MAGIC) + 8 * (n + 2)))   # header, filled in below
        for j_species, m in enumerate(matrix_array):
            offsets[j_species] = f.tell()
            pic.dump(m, f, protocol=pic_protocol)
        offsets[n] = f.tell()
        f.seek(0)
        f.write(_PACK_MAGIC)
        f.write(np.array([n], dtype="<i8").tobytes())
        f.write(offsets.tobytes())
    os.replace(fn + ".tmp", fn)   # readers never see a partly written file

def load_matrix(name, i_species, j_species, d_pickle): 
    fn = _pack_fn(name, i_species, d_pickle)
    if os.path.exists(fn):
        with open(fn, "rb") as f:
            # Read only this block's entry of the offset index (8 bytes), not
            # the whole index: column access reads one block from every pack.
            if f.read(len(_PACK_MAGIC)) != _PACK_MAGIC:
                raise ValueError("Not a matrix pack file: %s" % fn)
            n = int(np.frombuffer(f.read(8), dtype="<i8")[0])
            if not 0 <= j_species < n:
                raise IndexError("Block %d not in %s (%d blocks)" % (j_species, fn, n))
            f.seek(len(_PACK_MAGIC) + 8 + 8 * j_species)
            f.seek(int(np.frombuffer(f.read(8), dtype="<i8")[0]))
            return pic.load(f)
    with open(d_pickle + "%s%d_%d.pic" % (name, i_species, j_species), 'rb') as pic_file:  
        M = pic.load(pic_file)
    return M
        
def load_matrix_array(name, seqs_info, i_species, d_pickle, row=True):
    if row:
        fn = _pack_fn(name, i_species, d_pickle)
        if os.path.exists(fn):
            # The whole row is in one file: open it once, read the blocks in order.
            with open(fn, "rb") as f:
                offsets = _read_pack_index(f)
                if len(offsets) - 1 != seqs_info.n_species:
                    raise ValueError("%s holds %d blocks, expected %d" % (fn, len(offsets) - 1, seqs_info.n_species))
                f.seek(int(offsets[0]))
                return [pic.load(f) for _ in range(seqs_info.n_species)]
    matrix_array = []
    for j_species in range(seqs_info.n_species):
        if row == True:
            matrix_array.append(load_matrix(name, i_species, j_species, d_pickle))
        else:
            matrix_array.append(load_matrix(name, j_species, i_species, d_pickle))
    return matrix_array
              
def matrices_and_s(xarr, yarr):
    zarr = []
    for x, y in zip(xarr, yarr):
        zarr.append(x.multiply(y))
    return zarr
                
def matrices_and_tr_s(xarr, yarr):
    zarr = []
    for x, y in zip(xarr, yarr):
        zarr.append(x.multiply(y.transpose()))
    return zarr   
    
def delete_matrices(base_name, d_pickle):
    # Same name matching as before: e.g. "B" also matches the "BH" matrices.
    for pattern in (base_name + "*_*.pic", base_name + "*.pack", base_name + "*.pack.tmp"):
        for f in glob.glob(d_pickle + pattern):
            if os.path.exists(f): os.remove(f)

def csr_from_sorted(data, rows, cols, shape):
    """
    The CSR matrix of entries already in row order (and column order within a
    row, no duplicates): the same arrays as csr_matrix((data, (rows, cols))), a
    third of the time (it is made once or more for each pair of species). Entries
    not in that order are given to csr_matrix.
    """
    from scipy import sparse
    rows = np.asarray(rows)
    cols = np.asarray(cols)
    if rows.size > 1:
        key = rows.astype(np.int64) * shape[1] + cols
        if not (key[1:] > key[:-1]).all():
            return sparse.csr_matrix((data, (rows, cols)), shape=shape)
    indptr = np.zeros(shape[0] + 1, dtype=np.int64)
    np.cumsum(np.bincount(rows, minlength=shape[0]), out=indptr[1:])
    return sparse.csr_matrix((data, cols, indptr), shape=shape)


def sparse_max_row(csr_mat):
    ret = np.zeros(csr_mat.shape[0])
    ret[np.diff(csr_mat.indptr) != 0] = np.maximum.reduceat(csr_mat.data,csr_mat.indptr[:-1][np.diff(csr_mat.indptr)>0])
    return ret