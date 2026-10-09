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
import glob
import numpy as np
import datetime
from collections import namedtuple
from . import logging as run_logging
from . import file_io
from ..citation import citation
from ..phylo import tree

import shutil
import traceback

try:
    from rich import print, progress
    from rich.console import Console
    from rich.theme import Theme

except ImportError:
    ...
try:
    width = os.get_terminal_size().columns
except OSError as e:
    width = 80

PY2 = sys.version_info <= (3,)
csv_write_mode = "wb" if PY2 else "wt"
csv_append_mode = "ab" if PY2 else "at"
csv_read_mode = "rb" if PY2 else "rt"


class OrthoPrinter:
    def __init__(self):
        self.console = Console(
            theme=Theme(
                {
                    "info": "deep_sky_blue2",
                    "path": "green4",
                    "orthofinder": "dark_goldenrod",
                    "other_methods": "dark_orange3",
                    "version": "dark_cyan",
                    "success": "bold green",
                    "error": "bold red",
                    "warning": "bold yellow",
                    "default": "default on default",  # bright_white
                }
            ),
            force_terminal=True,
        )

    def print(self, *messages, style: str = "default", sep: str = " ", end: str = "\n"):
        message = sep.join(map(str, messages))
        if style in ("error", "warning"):
            run_logging.RunLogger.message(message, level=style.upper())
        self.console.print(message, style=style, end=end)


printer = OrthoPrinter()

from . import parallel_task_manager


def print_traceback(e):
    if PY2:
        traceback.print_exc()
    else:
        traceback.print_tb(e.__traceback__)


def stderr_exempt(stderr):
    ok_line_starts = {
        "diamond v",
        "Licensed under the GNU GPL",
        "Check http://github.com/",
    }
    try:
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stderr = stderr.encode()
    lines = stderr.split("\n")
    for line in lines:
        if line.rstrip() == "":
            continue
        if any(line.startswith(x) for x in ok_line_starts):
            continue
        return False
    return True


"""
Utilities
-------------------------------------------------------------------------------
"""
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


def get_directory_name(
    base_dir_name,
    i,
    search_program,
    msa_program,
    tree_program,
    scorematrix,
    gapopen,
    gapextend,
    extended_filename,
):

    if not extended_filename:
        if i == 0:
            return base_dir_name + os.sep
        else:
            return base_dir_name + ("_%d" % i) + os.sep
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
            extension = (
                matrix_name
                + gap_penalty
                + "-"
                + "-".join((search_program, msa_program, tree_program))
            )
        else:
            extension = "-".join((search_program, msa_program, tree_program))

        if i == 0:
            return base_dir_name + "_" + extension + os.sep
        else:
            return base_dir_name + ("_%d" % i) + "_" + extension + os.sep


"""Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""


def create_new_working_directory(
    base_directory_name,
    q_date=True,
    search_program=None,
    msa_program=None,
    tree_program=None,
    scorematrix=None,
    gapopen=None,
    gapextend=None,
    extended_filename=False,
    makedir=True,
):
    date_str = datetime.date.today().strftime("%b%d") if q_date else ""
    i_append = 0
    new_directory_name = get_directory_name(
        base_directory_name + date_str,
        i_append,
        search_program,
        msa_program,
        tree_program,
        scorematrix,
        gapopen,
        gapextend,
        extended_filename,
    )
    if makedir:
        while os.path.exists(new_directory_name):
            i_append += 1
            new_directory_name = get_directory_name(
                base_directory_name + date_str,
                i_append,
                search_program,
                msa_program,
                tree_program,
                scorematrix,
                gapopen,
                gapextend,
                extended_filename,
            )
        os.mkdir(new_directory_name)
    return new_directory_name


# Unused
def create_new_paired_directories(
    base_directory_name1,
    base_directory_name2,
    search_program=None,
    msa_program=None,
    tree_program=None,
    scorematrix=None,
    gapopen=None,
    gapextend=None,
    extended_filename=False,
):
    date_str = datetime.date.today().strftime("%b%d")
    i_append = 0
    new_directory_name1 = get_directory_name(
        base_directory_name1 + date_str,
        i_append,
        search_program,
        msa_program,
        tree_program,
        scorematrix,
        gapopen,
        gapextend,
        extended_filename,
    )
    new_directory_name2 = get_directory_name(
        base_directory_name2 + date_str,
        i_append,
        search_program,
        msa_program,
        tree_program,
        scorematrix,
        gapopen,
        gapextend,
        extended_filename,
    )
    while os.path.exists(new_directory_name1) or os.path.exists(new_directory_name2):
        i_append += 1
        new_directory_name1 = get_directory_name(
            base_directory_name1 + date_str,
            i_append,
            search_program,
            msa_program,
            tree_program,
            scorematrix,
            gapopen,
            gapextend,
            extended_filename,
        )
        new_directory_name2 = get_directory_name(
            base_directory_name2 + date_str,
            i_append,
            search_program,
            msa_program,
            tree_program,
            scorematrix,
            gapopen,
            gapextend,
            extended_filename,
        )
    os.mkdir(new_directory_name1)
    os.mkdir(new_directory_name2)
    return new_directory_name1, new_directory_name2


def get_unused_filename(base_filename, ext, i_unassigned=None):
    i_append = 0
    new_filename = base_filename + ext
    if i_unassigned is not None:
        while os.path.exists(new_filename):
            i_append += 1
            new_filename = base_filename + ("_%d" % i_append) + ext
    return new_filename, i_append


def sort_array_pair_by_first(use_for_sort_ar, keep_aligned_ar, q_largest_first=False):
    sorted_tuples = sorted(zip(use_for_sort_ar, keep_aligned_ar), reverse=q_largest_first)
    use_for_sort_ar = [i for i, j in sorted_tuples]
    keep_aligned_ar = [j for i, j in sorted_tuples]
    return use_for_sort_ar, keep_aligned_ar


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


def print_continue_hint(option, directory, what):
    """After a run stopped early (-op, -og, -ogt, -ost): how to continue it."""
    from rich.markup import escape
    printer.print("\nTo continue this analysis from %s, run:" % what)
    printer.print("    [bold]orthofinder %s [dark_cyan]%s[/dark_cyan][/bold]" % (option, escape(directory)))


def print_run_end(d_results, continue_option=None, continue_from=None):
    """
    The end of a run, also of one stopped early (-op, -og, -ogh, -ogt, -ost):
    the results directory, how to continue a stopped run, and the citation.
    """
    from rich.markup import escape
    d_results = os.path.normpath(d_results) + os.path.sep
    printer.print("\n[bold]Results directory:[/bold]")
    printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(d_results))
    if continue_option is not None:
        print_continue_hint(continue_option, d_results, continue_from)
    print_citation(d_results)


def print_restart_info(options, directory):
    """At the start of a restarted run (-b, -fg, -fgt, -fst): what it continues from."""
    from rich.markup import escape
    if options.q_start_from_blast and not options.q_start_from_fasta:
        what = "the sequence search results"
    elif options.q_start_from_groups:
        what = "the orthogroups"
    elif options.q_start_from_trees:
        what = "the gene trees"
    elif options.q_start_from_species_trees:
        what = "the species tree"
    else:
        return
    printer.print("\n[bold]Restarting[/bold] from %s in:" % what)
    printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(directory))


def success():
    parallel_task_manager.success()


def fail():
    parallel_task_manager.fail()


"""
IDExtractor
-------------------------------------------------------------------------------
"""


def get_id_pair_from_string(line):
    return list(map(int, line.split("_")))


class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different
    ways."""

    def get_id_to_name_dict(self):
        raise NotImplementedError("Should not be implemented")


def clean_accession(accession):
    """
    A FASTA header as a gene name: characters that would break the output
    files are replaced by "_" (":" "," "(" ")" in trees; a tab, possible in a
    full header line, would split the name across columns in the unquoted
    HOG files).
    """
    return (
        accession.replace(":", "_")
        .replace(",", "_")
        .replace("(", "_")
        .replace(")", "_")
        .replace("\t", "_")
    )  # .replace(".", "_")


class GeneNames(object):
    """
    The gene names of one species, kept unique, also as they are written in
    the gene trees (where the characters newick cannot hold, ; [ ] = ..., are
    replaced by "_"). Two headers that differ only in such characters
    ("gene:1" and "gene,1", or "x[1]" and "x_1_") would give the same name,
    and the genes could not be told apart in the trees and orthogroup
    files: the later one gets a suffix ("gene_1_2"), the next free "_<n>" in
    the order of the input file.
    """

    _TREE_UNSAFE = None

    def __init__(self):
        self.used = set()            # the names as written in trees

    @classmethod
    def _tree_form(cls, name):
        if cls._TREE_UNSAFE is None:
            import re
            from ..phylo import newick
            cls._TREE_UNSAFE = re.compile("[" + newick._ILEGAL_NEWICK_CHARS + "]")
        return cls._TREE_UNSAFE.sub("_", name)

    def add(self, accession):
        """The unique name of the next gene with this header."""
        name = clean_accession(accession)
        tree_form = self._tree_form(name)
        if tree_form in self.used:
            n = 2
            while "%s_%d" % (tree_form, n) in self.used:
                n += 1
            name, tree_form = "%s_%d" % (name, n), "%s_%d" % (tree_form, n)
        self.used.add(tree_form)
        return name


class FullAccession(IDExtractor):
    def __init__(self, ids_filename):
        # only want the first part and nothing else (easy!)
        self.id_to_name_dict = dict()
        names_in_species = {}
        with open(ids_filename, "r") as ids_file:
            for line in ids_file:
                line = line.rstrip()
                if not line:
                    continue
                #                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                id = id.strip()
                if id in self.id_to_name_dict:
                    raise RuntimeError(
                        "ERROR: A duplicate id was found in the fasta files: % s" % id
                    )
                # Names are unique within a species (see GeneNames)
                i_sp = id.split("_")[0]
                if i_sp not in names_in_species:
                    names_in_species[i_sp] = GeneNames()
                self.id_to_name_dict[id] = names_in_species[i_sp].add(accession)

    def get_id_to_name_dict(self):
        return self.id_to_name_dict


class FirstWordExtractor(IDExtractor):
    def __init__(self, ids_filename):
        # only want the first part and nothing else (easy!)
        self.id_to_name_dict = dict()
        accs_in_species = []
        with open(ids_filename, "r") as ids_file:
            for line in ids_file:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                i_sp = int(id.split("_")[0])
                while len(accs_in_species) < i_sp + 1:
                    accs_in_species.append(set())
                # Replace problematic characters
                accession = (
                    accession.replace(":", "_")
                    .replace(",", "_")
                    .replace("(", "_")
                    .replace(")", "_")
                )  # .replace(".", "_")
                # Only ensure the accessions are unique within the species, there's no need to ensure global uniqueness
                if accession in accs_in_species[i_sp]:
                    raise RuntimeError(
                        "A duplicate accession was found using just first part: % s"
                        % accession
                    )
                accs_in_species[i_sp].add(accession)
                if id in self.id_to_name_dict:
                    raise RuntimeError(
                        "ERROR: A duplicate id was found in the fasta files: % s" % id
                    )
                self.id_to_name_dict[id] = accession

    def get_id_to_name_dict(self):
        return self.id_to_name_dict


def have_support_values(species_tree_fn_ids):
    q_have_support = False
    try:
        tree.Tree(species_tree_fn_ids, format=2)
        q_have_support = True
    except:
        pass
    return q_have_support


def rename_tree_taxa(
    tree_fn_or_tree,
    new_tree_filename,
    ids_map,
    q_support,
    q_fix_negatives=False,
    in_format=None,
    label=None,
    q_via_copy=False,
):
    """
    qViaCopy - create a copy of the tree and edit this copy. I.e. don't make changes to the original
    """
    if label != None:
        q_support = False
    q_have_support = False
    try:
        if type(tree_fn_or_tree) is tree.TreeNode:
            if q_via_copy:
                t = tree_fn_or_tree.copy("newick")
            else:
                t = tree_fn_or_tree
        else:
            q_have_support = False
            if in_format == None:
                try:
                    t = tree.Tree(tree_fn_or_tree, format=2)
                    q_have_support = True
                except:
                    t = tree.Tree(tree_fn_or_tree)
            else:
                t = tree.Tree(tree_fn_or_tree, format=in_format)
        for node in t.get_leaves():
            node.name = ids_map[node.name]
        if q_fix_negatives:
            tree_length = sum([n.dist for n in t.traverse() if n != t])
            sliver = tree_length * 1e-6
        i_node = 1
        for n in t.traverse():
            if q_fix_negatives and n.dist < 0.0:
                n.dist = sliver
            if label != None:
                if (not n.is_leaf()) and (not n.is_root()):
                    n.name = label + ("%d" % i_node)
                    i_node += 1
        if label != None:
            with open(new_tree_filename, "w") as outfile:
                outfile.write(
                    t.write(format=3)[:-1] + label + "0;"
                )  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        elif t.name == "N0" or t.name == "n0":
            with open(new_tree_filename, "w") as outfile:
                outfile.write(
                    t.write(format=3)[:-1] + t.name + ";"
                )  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        else:
            if q_support or q_have_support:
                t.write(outfile=new_tree_filename, format=2)
            else:
                t.write(outfile=new_tree_filename, format=5)
    except:
        pass


"""
Find results of previous run    
-------------------------------------------------------------------------------
"""


def get_species_directory():
    # Confirms all required Sequence files and BLAST etc are present
    pass


def write_citation(d):
    with open(d + "Citation.txt", "w") as outfile:
        outfile.write(citation)


def print_citation(d=None):
    if d is not None:
        write_citation(d)
    # printer.print(print_citation)
    printer.print("\n[bold]CITATION:[/bold]")
    printer.print(
        " When publishing work that uses [dark_goldenrod]OrthoFinder[/dark_goldenrod] please cite:"
    )
    printer.print(
        " Emms D.M., Liu Y., Belcher L., Holmes J. & Kelly S. (2026), "
        "Nature Methods 23:1327–1333, ",
        end="",
    )
    printer.print(
        "[dark_cyan]"
        "[link=https://doi.org/10.1038/s41592-026-03126-6]"
        "https://doi.org/10.1038/s41592-026-03126-6"
        "[/link]\n"
    )
    # printer.print(" Emms D.M. & Kelly S. (2019), Genome Biology 20:238\n")

    printer.print(" If you use the species tree in your work then please also cite:")
    printer.print(" Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278")
    printer.print(" Emms D.M. & Kelly S. (2018), bioRxiv ", end="")
    printer.print(
        "[dark_cyan][link=https://doi.org/10.1101/267914]https://doi.org/10.1101/267914[/link]"
    )


def print_underline(text, q_heavy=False):
    run_logging.RunLogger.message(text)
    n = len(text)
    if text.startswith("\n"):
        n -= 1
    print(("\n[bold]" + text + "[/bold]") if q_heavy else ("\n" + text))
    print((("=" if q_heavy else "-") * n))


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


def number_open_files_exception_advice(n_species, q_at_trees):
    """
    Prints advice for user on "IOError: [Errno 24] Too many open files" exception
    Args:
        n_species - the number of species in the analysis
        q_at_trees - has this error occurred at the orthologs from trees stage
    """
    # parallel_task_manager.run_command("ulimit -Hn")
    n_req = n_species + 100
    msg = (
        "\nERROR: The system limits on the number of files a process can open is too low. For %d species \
OrthoFinder needs to be able to open at least r=%d files. Please increase the limit and restart OrthoFinder\n\
1. Check the hard and soft limits on the number of open files for your system:\n\
    $ ulimit -Hn\n\
    $ ulimit -Sn\n\
2. If hard limit, h > r already, then you just need to increase the soft limit:\n\
    $ ulimit -n %d\n\
3. Alternatively, if h < r then you need to edit the file '/etc/security/limits.conf', this requires root privileges. \
To increase the limit to %d for user  called 'emms' add the lines:\n\
    emms hard nofile %d\n\
    emms soft nofile %d\n"
        % (n_species, n_req, n_req, n_req, n_req, n_req)
    )
    msg += "    (edit these lines to match your username)\n\
4. Check the limit has now been updated (if you changed the hard limit you'll need to open a new session and confirm it's updated):\n\
    $ ulimit -Sn"

    if q_at_trees:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder 'from trees' using the '-ft' command"
    else:
        msg_part_2 = (
            "5. Once the limit is updated restart OrthoFinder with the original command"
        )
    msg_part_3 = (
        "\nFor full details see: https://github.com/davidemms/OrthoFinder/issues/384"
    )
    print(msg + "\n" + msg_part_2 + "\n" + msg_part_3 + "\n")
    print(
        "Since this issue OrthoFinder has been updated to only require O(n) open files for n species rather than O(n^2)."
    )
    print("Please follow the above advice to ensure that it is able to.")


"""
-------------------------------------------------------------------------------
"""


class NOrthoSp(object):
    """matrix of number of genes in species i that have orthologues/an ortholog in species j"""

    def __init__(self, n_sp):
        self.n = np.zeros((n_sp, n_sp))
        self.n_121 = np.zeros((n_sp, n_sp))  # genes in i that have one ortholog in j
        self.n_12m = np.zeros((n_sp, n_sp))  # genes in i that have many orthologs in j
        self.n_m21 = np.zeros(
            (n_sp, n_sp)
        )  # genes in i that are in a many-to-one orthology relationship with genes in j
        self.n_m2m = np.zeros(
            (n_sp, n_sp)
        )  # genes in i that are in a many-to-many orthology relationship with genes in j

    def __iadd__(self, other):
        self.n += other.n
        self.n_121 += other.n_121
        self.n_12m += other.n_12m
        self.n_m21 += other.n_m21
        self.n_m2m += other.n_m2m
        return self

    def add_pair(self, i_l, i_r, n_l, n_r):
        """Record one orthology relationship between nL genes in iL and nR genes in iR."""
        self.n[i_l, i_r] += n_l
        self.n[i_r, i_l] += n_r

        if n_l == 1 and n_r == 1:
            self.n_121[i_l, i_r] += 1
            self.n_121[i_r, i_l] += 1

        elif n_l == 1:
            self.n_12m[i_l, i_r] += 1
            self.n_m21[i_r, i_l] += n_r

        elif n_r == 1:
            self.n_m21[i_l, i_r] += n_l
            self.n_12m[i_r, i_l] += 1

        else:
            self.n_m2m[i_l, i_r] += n_l
            self.n_m2m[i_r, i_l] += n_r


class NOrthoCache(object):
    """matrix of approx number of unwritten cached orthologs"""

    def __init__(self, n_sp):
        self.n = np.zeros((n_sp, n_sp))

    def __iadd__(self, n_ortho_sp_obj):
        self.n += n_ortho_sp_obj.n
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
            for i, j in zip(I, J):
                self.n[i, j] = 0
        return I, J


class Finalise(object):
    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManagerSingleton()
        ptm.stop()


def seq_id_key(seq_id):
    """Sort key for "species_sequence" IDs in numeric order (0_61 before 0_317)."""
    i_sp, i_seq = seq_id.split("_")
    return int(i_sp), int(i_seq)


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
            filenames = []
            with os.scandir(src_dir) as it:
                for entry in it:
                    if entry.is_file():
                        # fn = entry.name.partition("_")[0]
                        # fn = entry.name.partition(".")[0]
                        filenames.append((entry.name.partition(".")[0], entry))
            filenames = sorted(filenames)
            for fn, entry in filenames:
                try:
                    with open(entry.path, "r") as infile:
                        og_tree = "".join(line.strip() for line in infile)
                    writer.write(f"{fn}: {og_tree}\n")
                except IOError as e:
                    print(f"ERROR reading file {entry.name}: {e}")
    except IOError as e:
        printer.print(f"ERROR writing to output file {output_fn}: {e}", style="error")


def cleanup_path(path):
    try:
        if os.path.isfile(path) or os.path.islink(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))


def get_progressbar(len_task, visible=True):
    progressbar = progress.Progress(
        progress.TextColumn("[progress.description]{task.description}"),
        progress.BarColumn(bar_width=width // 2),
        progress.SpinnerColumn(),
        progress.MofNCompleteColumn(),
        progress.TimeElapsedColumn(),
        transient=False,
        # progress.TextColumn("{task.completed}/{task.total}")
    )

    task = progressbar.add_task(
            "[yellow]Processing...", total=len_task, visible=visible
    )
    return progressbar, task


def clear_dir(of3_dir):
    if os.path.exists(of3_dir):
        with os.scandir(of3_dir) as entries:
            for entry in entries:
                try:
                    if entry.is_file() or entry.is_symlink():
                        os.unlink(entry.path)
                    elif entry.is_dir():
                        shutil.rmtree(entry.path)

                except Exception as e:
                    printer.print(
                        f"Failed to delete {entry.path}. Reason: {e}", style="error"
                    )


def split_ortholog_files(d_orthologs, q_compress=False):
    if not d_orthologs.endswith("/"):
        d_orthologs += "/"
    filenames = list(glob.glob(d_orthologs + "*.tsv") + glob.glob(d_orthologs + "*.tsv.gz"))
    species = [
        os.path.splitext(os.path.splitext(os.path.basename(fn))[0])[0]
        if fn.endswith(".gz")
        else os.path.splitext(os.path.basename(fn))[0]
        for fn in filenames
    ]
    for fn, sp0 in zip(filenames, species):
        d_out = d_orthologs + "Orthologues_" + sp0 + "/"
        if not os.path.exists(d_out):
            os.mkdir(d_out)
        file_handles = []
        csv_writers = {}
        for sp1 in species:
            if sp0 == sp1:
                continue
            if q_compress:
                file_handles.append(
                    gzip.open(d_out + "%s__v__%s.tsv.gz" % (sp0, sp1), csv_write_mode)
                )
            else:
                file_handles.append(
                    open(d_out + "%s__v__%s.tsv" % (sp0, sp1), csv_write_mode)
                )
            csv_writers[sp1] = file_io.writer(file_handles[-1])
            csv_writers[sp1].writerow(("Orthogroup", sp0, sp1))
        with (
            gzip.open(fn, csv_read_mode)
            if fn.endswith(".gz")
            else open(fn, csv_read_mode) as infile
        ):
            reader = file_io.unquoted_reader(infile)   # rows are written unquoted
            next(reader)  # skip header
            for row in reader:
                if len(row) == 4:  # OG,species,genes1,genes2
                    csv_writers[row[1]].writerow(row[:1] + row[2:])
