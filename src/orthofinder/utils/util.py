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
import gzip
import os
import sys
import numpy as np
import datetime
from collections import namedtuple
from ..citation import citation, print_citation
from ..tools import tree
from . import parallel_task_manager
import shutil
import traceback

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'
csv_append_mode = 'ab' if PY2 else 'at'
csv_read_mode = 'rb' if PY2 else 'rt'

def print_traceback(e):
    if PY2:
        traceback.print_exc()
    else:
        traceback.print_tb(e.__traceback__)

def stderr_exempt(stderr):
    ok_line_starts = {"diamond v", "Licensed under the GNU GPL", "Check http://github.com/"}
    try:
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stderr = stderr.encode()
    lines = stderr.split("\n")
    for line in lines:
        if line.rstrip() == "": continue
        if any(line.startswith(x) for x in ok_line_starts): continue
        return False
    return True

"""
Utilities
-------------------------------------------------------------------------------
"""
SequencesInfo = namedtuple("SequencesInfo", "nSeqs nSpecies speciesToUse seqStartingIndices nSeqsPerSpecies")
"""
nSeqs: int - number of sequences in speciesToUse
nSpecies: int - number in speciesToUse
speciesToUse: List[int]
seqStartingIndices: List[int] of size |speciesToUse|, reindexed from zero according to the speices present
nSeqsPerSpecies: Dict[int, int] - indexed by OrthoFinder species ID, info on all species, nto just those included in analysis
"""

# Get Info from seqs IDs file?
def GetSeqsInfo(inputDirectory_list, speciesToUse, nSpAll):
    seqStartingIndices = [0]
    nSeqs = 0
    nSeqsPerSpecies = dict()
    for iFasta in range(nSpAll):
        for d in inputDirectory_list:
            fastaFilename = d + "Species%d.fa" % iFasta
            if os.path.exists(fastaFilename): break
        n = 0
        with open(fastaFilename) as infile:
            for line in infile:
                if len(line) > 1 and line[0] == ">":
                    n += 1
        nSeqsPerSpecies[iFasta] = n
        if iFasta in speciesToUse:
            nSeqs += n
            seqStartingIndices.append(nSeqs)
    seqStartingIndices = seqStartingIndices[:-1]
    nSpecies = len(speciesToUse)
    return SequencesInfo(nSeqs=nSeqs, nSpecies=nSpecies, speciesToUse=speciesToUse,
                         seqStartingIndices=seqStartingIndices, nSeqsPerSpecies=nSeqsPerSpecies)

def SeqsInfoRecompute(seqs_info_orig, new_species_to_use):
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
    seqStartingIndices = [0]
    for isp in new_species_to_use:
        seqStartingIndices.append(seqStartingIndices[-1] + seqs_info_orig.nSeqsPerSpecies[isp])
    nSeqs = seqStartingIndices[-1]
    seqStartingIndices = seqStartingIndices[:-1]
    return SequencesInfo(nSeqs=nSeqs, nSpecies=len(new_species_to_use), speciesToUse=new_species_to_use,
                         seqStartingIndices=seqStartingIndices, nSeqsPerSpecies=seqs_info_orig.nSeqsPerSpecies)


class SpeciesInfo(object):
    def __init__(self):
        self.speciesToUse = []           #       seqsInfo.iSpeciesToUse   - which to include for this analysis
        self.nSpAll = None               #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.iFirstNewSpecies = None     #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices
    def __str__(self):
        return str((self.speciesToUse, self.nSpAll, self.iFirstNewSpecies))
    def get_original_species(self):
        if self.iFirstNewSpecies is None:
            return self.speciesToUse
        else:
            return [iSp for iSp in  self.speciesToUse if iSp < self.iFirstNewSpecies]


# picProtocol = 1
# version = "3.0.1"
    
def PrintNoNewLine(text):
    parallel_task_manager.PrintNoNewLine(text)

def PrintTime(message):
    parallel_task_manager.PrintTime(message)   
      
"""
Directory and file management
-------------------------------------------------------------------------------
"""               
               
def GetDirectoryName(baseDirName, i, search_program, msa_program, tree_program,
                     scorematrix, gapopen, gapextend, extended_filename):
    
    if not extended_filename:
        if i == 0:
            return baseDirName + os.sep
        else:
            return baseDirName + ("_%d" % i) + os.sep
    else:
        if scorematrix and os.path.isfile(scorematrix):
            matrix_name = os.path.basename(scorematrix).split(".")[0]
        else:
            matrix_name = scorematrix if scorematrix else ""

        if gapopen and gapextend:
            gap_penalty = "-" + "-".join((gapopen, gapextend))
        else:
            gap_penalty = ""
        
        if len(matrix_name) != 0 or len(gap_penalty) != 0:
            extension = matrix_name + gap_penalty + "-" + "-".join((search_program, msa_program, tree_program))
        else:
            extension = "-".join((search_program, msa_program, tree_program))

        if i == 0:
            return baseDirName + "_" + extension + os.sep 
        else:
            return baseDirName + ("_%d" % i) + "_" + extension + os.sep

"""Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""
def CreateNewWorkingDirectory(baseDirectoryName, 
                              qDate=True,
                              search_program=None,
                              msa_program=None,
                              tree_program=None, 
                              scorematrix=None, 
                              gapopen=None, 
                              gapextend=None,
                              extended_filename=False):
    dateStr = datetime.date.today().strftime("%b%d") if qDate else ""
    iAppend = 0
    newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, 
                                        iAppend,
                                        search_program, 
                                        msa_program,
                                        tree_program, 
                                        scorematrix,
                                        gapopen, 
                                        gapextend,
                                        extended_filename)
    while os.path.exists(newDirectoryName):
        iAppend += 1
        newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, 
                                            iAppend, 
                                            search_program,
                                            msa_program,
                                            tree_program, 
                                            scorematrix, 
                                            gapopen, 
                                            gapextend,
                                            extended_filename)
    os.mkdir(newDirectoryName)
    return newDirectoryName

# Unused    
def CreateNewPairedDirectories(baseDirectoryName1, 
                               baseDirectoryName2,
                               search_program=None, 
                               msa_program=None, 
                               tree_program=None,                              
                               scorematrix=None, 
                               gapopen=None, 
                               gapextend=None,
                               extended_filename=False):
    dateStr = datetime.date.today().strftime("%b%d") 
    iAppend = 0
    newDirectoryName1 = GetDirectoryName(baseDirectoryName1 + dateStr,
                                        iAppend,
                                        search_program,
                                        msa_program,
                                        tree_program, 
                                        scorematrix, 
                                        gapopen, 
                                        gapextend,
                                        extended_filename)
    newDirectoryName2 = GetDirectoryName(baseDirectoryName2 + dateStr, 
                                        iAppend,
                                        search_program,
                                        msa_program,
                                        tree_program,  
                                        scorematrix, 
                                        gapopen, 
                                        gapextend,
                                        extended_filename)
    while os.path.exists(newDirectoryName1) or os.path.exists(newDirectoryName2):
        iAppend += 1
        newDirectoryName1 = GetDirectoryName(baseDirectoryName1 + dateStr, 
                                            iAppend,
                                            search_program,
                                            msa_program,
                                            tree_program,  
                                            scorematrix, 
                                            gapopen, 
                                            gapextend,
                                            extended_filename)
        newDirectoryName2 = GetDirectoryName(baseDirectoryName2 + dateStr, 
                                            iAppend, 
                                            search_program,
                                            msa_program,
                                            tree_program,  
                                            scorematrix, 
                                            gapopen, 
                                            gapextend,
                                            extended_filename)
    os.mkdir(newDirectoryName1)
    os.mkdir(newDirectoryName2)
    return newDirectoryName1, newDirectoryName2

def GetUnusedFilename(baseFilename, ext):
    iAppend = 0
    newFilename = baseFilename + ext
    while os.path.exists(newFilename):
        iAppend += 1
        newFilename = baseFilename + ("_%d" % iAppend) + ext
    return newFilename, iAppend
       
def SortArrayPairByFirst(useForSortAr, keepAlignedAr, qLargestFirst=False):
    sortedTuples = sorted(zip(useForSortAr, keepAlignedAr), reverse=qLargestFirst)
    useForSortAr = [i for i, j in sortedTuples]
    keepAlignedAr = [j for i, j in sortedTuples]
    return useForSortAr, keepAlignedAr

def GetSpeciesToUse(speciesIDsFN):
    """Returns species indices (int) to use and total number of species available """
    speciesToUse = []
    speciesToUse_names = []
    nSkipped = 0
    with open(speciesIDsFN, 'r') as speciesF:
        for line in speciesF:
            line = line.rstrip()
            if not line: continue
            if line.startswith("#"): nSkipped += 1
            else: 
                iSp, spName = line.split(": ")
                speciesToUse.append(int(iSp))
                speciesToUse_names.append(spName)
    return speciesToUse, len(speciesToUse) + nSkipped, speciesToUse_names
 
def Success():
    parallel_task_manager.Success()
   
def Fail():
    parallel_task_manager.Fail()
    
"""
IDExtractor
-------------------------------------------------------------------------------
"""

def GetIDPairFromString(line):
    return list(map(int, line.split("_")))

class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different 
    ways."""
    def GetIDToNameDict(self):
        raise NotImplementedError("Should not be implemented")

class FullAccession(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                line = line.rstrip()
                if not line: continue
#                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                id = id.strip()
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession 
                
    def GetIDToNameDict(self):
        return self.idToNameDict
                
class FirstWordExtractor(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        accs_in_species = []
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                iSp = int(id.split("_")[0])
                while len(accs_in_species) < iSp + 1:
                    accs_in_species.append(set())
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                # Only ensure the accessions are unique within the species, there's no need to ensure global uniqueness
                if accession in accs_in_species[iSp]:
                    raise RuntimeError("A duplicate accession was found using just first part: % s" % accession)
                accs_in_species[iSp].add(accession)
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession     
                
    def GetIDToNameDict(self):
        return self.idToNameDict

def HaveSupportValues(speciesTreeFN_ids):
    qHaveSupport = False
    try:
        tree.Tree(speciesTreeFN_ids, format=2)
        qHaveSupport = True
    except:
        pass
    return qHaveSupport

def RenameTreeTaxa(treeFN_or_tree, newTreeFilename, idsMap, qSupport, qFixNegatives=False, inFormat=None, label=None, qViaCopy=False):
    """
    qViaCopy - create a copy of the tree and edit this copy. I.e. don't make changes to the original 
    """
    if label != None: qSupport = False
    qHaveSupport = False
    try:
        if type(treeFN_or_tree) is tree.TreeNode:
            if qViaCopy:
                t = treeFN_or_tree.copy("newick")
            else:
                t = treeFN_or_tree
        else:
            qHaveSupport = False
            if inFormat == None:
                try:
                    t = tree.Tree(treeFN_or_tree, format=2)
                    qHaveSupport = True
                except:
                    t = tree.Tree(treeFN_or_tree)
            else:
                t = tree.Tree(treeFN_or_tree, format=inFormat)
        for node in t.get_leaves():
            node.name = idsMap[node.name]
        if qFixNegatives:
            tree_length = sum([n.dist for n in t.traverse() if n != t])
            sliver = tree_length * 1e-6
        iNode = 1
        for n in t.traverse():
            if qFixNegatives and n.dist < 0.0: n.dist = sliver
            if label != None:
                if (not n.is_leaf()) and (not n.is_root()):
                    n.name = label + ("%d" % iNode)
                    iNode += 1
        if label != None:
            with open(newTreeFilename, 'w') as outfile:
                outfile.write(t.write(format=3)[:-1] + label + "0;")  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        elif t.name == "N0" or t.name == "n0":
            with open(newTreeFilename, 'w') as outfile:
                outfile.write(t.write(format=3)[:-1] + t.name + ";")  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        else:
            if qSupport or qHaveSupport:
                t.write(outfile = newTreeFilename, format=2)  
            else:
                t.write(outfile = newTreeFilename, format=5)  
    except:
        pass
    
"""
Find results of previous run    
-------------------------------------------------------------------------------
"""

def GetSpeciesDirectory():
    # Confirms all required Sequence files and BLAST etc are present
    pass

def WriteCitation(d):
    with open(d + "Citation.txt", 'w') as outfile:
        outfile.write(citation)

def PrintCitation(d=None):
    if d is not None: WriteCitation(d)
    print()
    print(print_citation)  
    # print ("\nCITATION:")  
    # print (" When publishing work that uses OrthoFinder please cite:")
    # print (" Emms D.M. & Kelly S. (2019), Genome Biology 20:238\n")   

    # print (" If you use the species tree in your work then please also cite:")
    # print (" Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278")
    # print (" Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914")

def PrintUnderline(text, qHeavy=False):
    print(("\n" + text))
    n = len(text)
    if text.startswith("\n"): n -= 1
    print((("=" if qHeavy else "-") * n))

def FlowText(text, n=60):
    """Split text onto lines of no more that n characters long
    """
    lines = ""
    while len(text) > 0:
        if len(lines) > 0: lines += "\n"
        if len(text) > n:
            # split at no more than 60
            iEnd = n
            while iEnd > 0 and text[iEnd] != " ": iEnd-=1
            if iEnd == 0:
                # there was nowhere to split it at a blank, just have to split at 60
                lines += text[:n]
                text = text[n:]
            else:
                # split at blank
                lines += text[:iEnd]
                text = text[iEnd+1:]  # skip blank
        else:
            lines += text
            text = ""
    return lines

def number_open_files_exception_advice(n_species, q_at_trees):
    """
    Prints advice for user on "IOError: [Errno 24] Too many open files" exception
    Args:
        n_species - the number of species in the analysis
        q_at_trees - has this error occurred at the orthologs from trees stage
    """
    # parallel_task_manager.RunCommand("ulimit -Hn")    
    n_req = n_species + 100
    msg="\nERROR: The system limits on the number of files a process can open is too low. For %d species \
OrthoFinder needs to be able to open at least r=%d files. Please increase the limit and restart OrthoFinder\n\
1. Check the hard and soft limits on the number of open files for your system:\n\
    $ ulimit -Hn\n\
    $ ulimit -Sn\n\
2. If hard limit, h > r already, then you just need to increase the soft limit:\n\
    $ ulimit -n %d\n\
3. Alternatively, if h < r then you need to edit the file '/etc/security/limits.conf', this requires root privileges. \
To increase the limit to %d for user  called 'emms' add the lines:\n\
    emms hard nofile %d\n\
    emms soft nofile %d\n" % (n_species, n_req, n_req, n_req, n_req, n_req)   
    msg +="    (edit these lines to match your username)\n\
4. Check the limit has now been updated (if you changed the hard limit you'll need to open a new session and confirm it's updated):\n\
    $ ulimit -Sn" 

    if q_at_trees:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder 'from trees' using the '-ft' command"
    else:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder with the original command"
    msg_part_3 = "\nFor full details see: https://github.com/davidemms/OrthoFinder/issues/384"
    print(msg + "\n" + msg_part_2 + "\n" + msg_part_3 + "\n")
    print("Since this issue OrthoFinder has been updated to only require O(n) open files for n species rather than O(n^2).")
    print("Please follow the above advice to ensure that it is able to.")
"""
-------------------------------------------------------------------------------
"""

class nOrtho_sp(object):
    """ matrix of number of genes in species i that have orthologues/an orthologue in species j"""
    def __init__(self, nSp):
        self.n = np.zeros((nSp, nSp))
        self.n_121 = np.zeros((nSp, nSp))  # genes in i that have one orthologue in j
        self.n_12m = np.zeros((nSp, nSp))  # genes in i that have many orthologues in j
        self.n_m21 = np.zeros((nSp, nSp))  # genes in i that are in a many-to-one orthology relationship with genes in j
        self.n_m2m = np.zeros((nSp, nSp))  # genes in i that are in a many-to-many orthology relationship with genes in j
        
    def __iadd__(self, other):
        self.n += other.n
        self.n_121 += other.n_121
        self.n_12m += other.n_12m
        self.n_m21 += other.n_m21
        self.n_m2m += other.n_m2m
        return self


class nOrtho_cache(object):
    """ matrix of approx number of unwritten cached orthologs"""
    def __init__(self, nSp):
        self.n = np.zeros((nSp, nSp))
        
    def __iadd__(self, nOrtho_sp_obj):
        self.n += nOrtho_sp_obj.n
        return self

    def get_i_j_to_write(self, n_max_cache, fewer_open_files=False):
        if fewer_open_files:
            # all genes for species i are collated
            IJ = np.where(self.n.sum(axis=1) > n_max_cache)
            I = list(IJ[0])
            J = None
            for i in I:
                self.n[i, :] = 0
        else:
            IJ = np.where(self.n > n_max_cache)
            I = list(IJ[0])
            J = list(IJ[1])
            for i, j in zip(I,J):
                self.n[i, j] = 0
        return I, J
        
class Finalise(object):
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()

def writerow(fh, row):
    # CSV format specifies CRLF line endings: https://tools.ietf.org/html/rfc4180
    fh.write("\t".join(map(str, row)) + "\r\n")

def getrow(row):
    # CSV format specifies CRLF line endings: https://tools.ietf.org/html/rfc4180
    return "\t".join(map(str, row)) + "\r\n"

def version_parse_simple(sem_version):
    return list(map(int, sem_version.split(".")[:3]))


def file_open(filename_with_gz, mode, gz):
    if gz:
        return gzip.open(filename_with_gz + ".gz", mode)
    else:
        return open(filename_with_gz, mode)


def compress_files(src_dir, output_fn):
    try:
        with open(output_fn, "w") as writer:
            with os.scandir(src_dir) as it:
                for entry in it:
                    if entry.is_file():
                        fn = entry.name.partition("_")[0]
                        try:
                            with open(entry.path, "r") as infile:
                                og_tree = "".join(line.strip() for line in infile)
                            writer.write(f"{fn}: {og_tree}\n")
                        except IOError as e:
                            print(f"ERROR reading file {entry.name}: {e}")
    except IOError as e:
        print(f"ERROR writing to output file {output_fn}: {e}") 

def cleanup_path(path):
    try:
        if os.path.isfile(path) or os.path.islink(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))