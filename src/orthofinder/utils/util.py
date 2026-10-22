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
import sys


PY2 = sys.version_info <= (3,)
csv_write_mode = "wb" if PY2 else "wt"
csv_append_mode = "ab" if PY2 else "at"
csv_read_mode = "rb" if PY2 else "rt"


from . import parallel_task_manager


"""
Utilities
-------------------------------------------------------------------------------
"""
    # seqStartingIndices = [0]
    # nSeqs = 0
    # nSeqsPerSpecies = dict()
    # for iFasta in range(nSpAll):
    #     for d in inputDirectory_list:
    #         fastaFilename = os.path.join(d, "Species%d.fa" % iFasta)
    #         if os.path.exists(fastaFilename):
    #             break
    #     n = 0
    #     with open(fastaFilename) as infile:
    #         for line in infile:
    #             if len(line) > 1 and line[0] == ">":
    #                 n += 1
    #     nSeqsPerSpecies[iFasta] = n
    #     if iFasta in speciesToUse:
    #         nSeqs += n
    #         seqStartingIndices.append(nSeqs)

    # seqStartingIndices = seqStartingIndices[:-1]
    # nSpecies = len(speciesToUse)

    # return SequencesInfo(
    #     nSeqs=nSeqs,
    #     nSpecies=nSpecies,
    #     speciesToUse=speciesToUse,
    #     seqStartingIndices=seqStartingIndices,
    #     nSeqsPerSpecies=nSeqsPerSpecies,
    # )


# picProtocol = 1
# version = "3.0.1"


def print_no_new_line(text):
    parallel_task_manager.print_no_new_line(text)


def print_time(message):
    parallel_task_manager.print_time(message)


"""
Directory and file management
-------------------------------------------------------------------------------
"""


"""Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""


def sort_array_pair_by_first(use_for_sort_ar, keep_aligned_ar, q_largest_first=False):
    sorted_tuples = sorted(zip(use_for_sort_ar, keep_aligned_ar), reverse=q_largest_first)
    use_for_sort_ar = [i for i, j in sorted_tuples]
    keep_aligned_ar = [j for i, j in sorted_tuples]
    return use_for_sort_ar, keep_aligned_ar


def success():
    parallel_task_manager.success()


def fail():
    parallel_task_manager.fail()


"""
IDExtractor
-------------------------------------------------------------------------------
"""


"""
Find results of previous run    
-------------------------------------------------------------------------------
"""


# def FlowText(text, n=60):
#     """Split text onto lines of no more that n characters long"""
#     lines = ""
#     while len(text) > 0:
#         if len(lines) > 0:
#             lines += "\n"
#         if len(text) > n:
#             # split at no more than 60
#             iEnd = n
#             while iEnd > 0 and text[iEnd] != " ":
#                 iEnd -= 1
#             if iEnd == 0:
#                 # there was nowhere to split it at a blank, just have to split at 60
#                 lines += text[:n]
#                 text = text[n:]
#             else:
#                 # split at blank
#                 lines += text[:iEnd]
#                 text = text[iEnd + 1 :]  # skip blank
#         else:
#             lines += text
#             text = ""
#     return lines


"""
-------------------------------------------------------------------------------
"""


class Finalise(object):
    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()


def version_parse_simple(sem_version):
    return list(map(int, sem_version.split(".")[:3]))


