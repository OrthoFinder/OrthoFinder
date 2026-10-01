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

from .. import picProtocol

# Matrices are stored in blocks, one per species pair (iSpecies, jSpecies).
# DumpMatrixArray writes all the blocks of one species into a single "pack"
# file: a header with the byte offset of every block, followed by the pickled
# blocks. That is n files per matrix instead of n^2 (e.g. ~1,270 instead of
# ~1.6 million at 1,270 species), while reading one block is still one open,
# one seek and one unpickle. DumpMatrix still writes single-block files, and
# LoadMatrix reads either kind.

_PACK_MAGIC = b"OFPACK1\n"


def _PackFN(name, iSpecies, d_pickle):
    return d_pickle + "%s%d.pack" % (name, iSpecies)


def _ReadPackIndex(f):
    if f.read(len(_PACK_MAGIC)) != _PACK_MAGIC:
        raise ValueError("Not a matrix pack file: %s" % f.name)
    n = int(np.frombuffer(f.read(8), dtype="<i8")[0])
    return np.frombuffer(f.read(8 * (n + 1)), dtype="<i8")


def DumpMatrix(name, m, iSpecies, jSpecies, d_pickle):
    with open(d_pickle + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'wb') as picFile:
        pic.dump(m, picFile, protocol=picProtocol)
    
def DumpMatrixArray(name, matrixArray, iSpecies, d_pickle):
    fn = _PackFN(name, iSpecies, d_pickle)
    n = len(matrixArray)
    offsets = np.zeros(n + 1, dtype="<i8")
    with open(fn + ".tmp", "wb") as f:
        f.write(b"\0" * (len(_PACK_MAGIC) + 8 * (n + 2)))   # header, filled in below
        for jSpecies, m in enumerate(matrixArray):
            offsets[jSpecies] = f.tell()
            pic.dump(m, f, protocol=picProtocol)
        offsets[n] = f.tell()
        f.seek(0)
        f.write(_PACK_MAGIC)
        f.write(np.array([n], dtype="<i8").tobytes())
        f.write(offsets.tobytes())
    os.replace(fn + ".tmp", fn)   # readers never see a partly written file

def LoadMatrix(name, iSpecies, jSpecies, d_pickle): 
    fn = _PackFN(name, iSpecies, d_pickle)
    if os.path.exists(fn):
        with open(fn, "rb") as f:
            # Read only this block's entry of the offset index (8 bytes), not
            # the whole index: column access reads one block from every pack.
            if f.read(len(_PACK_MAGIC)) != _PACK_MAGIC:
                raise ValueError("Not a matrix pack file: %s" % fn)
            n = int(np.frombuffer(f.read(8), dtype="<i8")[0])
            if not 0 <= jSpecies < n:
                raise IndexError("Block %d not in %s (%d blocks)" % (jSpecies, fn, n))
            f.seek(len(_PACK_MAGIC) + 8 + 8 * jSpecies)
            f.seek(int(np.frombuffer(f.read(8), dtype="<i8")[0]))
            return pic.load(f)
    with open(d_pickle + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'rb') as picFile:  
        M = pic.load(picFile)
    return M
        
def LoadMatrixArray(name, seqsInfo, iSpecies, d_pickle, row=True):
    if row:
        fn = _PackFN(name, iSpecies, d_pickle)
        if os.path.exists(fn):
            # The whole row is in one file: open it once, read the blocks in order.
            with open(fn, "rb") as f:
                offsets = _ReadPackIndex(f)
                if len(offsets) - 1 != seqsInfo.nSpecies:
                    raise ValueError("%s holds %d blocks, expected %d" % (fn, len(offsets) - 1, seqsInfo.nSpecies))
                f.seek(int(offsets[0]))
                return [pic.load(f) for _ in range(seqsInfo.nSpecies)]
    matrixArray = []
    for jSpecies in range(seqsInfo.nSpecies):
        if row == True:
            matrixArray.append(LoadMatrix(name, iSpecies, jSpecies, d_pickle))
        else:
            matrixArray.append(LoadMatrix(name, jSpecies, iSpecies, d_pickle))
    return matrixArray
              
def MatricesAnd_s(Xarr, Yarr):
    Zarr = []
    for x, y in zip(Xarr, Yarr):
        Zarr.append(x.multiply(y))
    return Zarr
                
def MatricesAndTr_s(Xarr, Yarr):
    Zarr = []
    for x, y in zip(Xarr, Yarr):
        Zarr.append(x.multiply(y.transpose()))
    return Zarr   
    
def DeleteMatrices(baseName, d_pickle):
    # Same name matching as before: e.g. "B" also matches the "BH" matrices.
    for pattern in (baseName + "*_*.pic", baseName + "*.pack", baseName + "*.pack.tmp"):
        for f in glob.glob(d_pickle + pattern):
            if os.path.exists(f): os.remove(f)

def sparse_max_row(csr_mat):
    ret = np.zeros(csr_mat.shape[0])
    ret[np.diff(csr_mat.indptr) != 0] = np.maximum.reduceat(csr_mat.data,csr_mat.indptr[:-1][np.diff(csr_mat.indptr)>0])
    return ret