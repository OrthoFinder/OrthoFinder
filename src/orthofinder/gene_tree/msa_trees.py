#!/usr/bin/env python3
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

"""
Created on Thu Sep 25 13:15:22 2014

@author: david
"""
import collections
import contextlib
import os
import glob
import pickle
import tempfile

import numpy as np
try:
    from rich import print
except ImportError:
    ...

from ..utils import util, parallel_task_manager
from ..utils import file_handler
from . import trim
from ..species_tree import from_gene_trees
from ..utils import printing
from ..phylo import taxa
from ..utils import command_runner
from ..species_tree.concatenated_alignment import create_concatenated_alignment, determine_orthogroups_for_species_tree


# The fewest sequences RAxML and RAxML-NG infer a tree of
RAXML_MIN_SEQS = 4


class FastaWriter(object):
    def __init__(self, fasta_file_dir_list, species_to_use, genes=None):
        """
        The list of previous directories could include species that are no longer used,
        but not ones that are conflicting
        genes: the sequence IDs to keep (None: all), e.g. those of a few orthogroups,
            so that every sequence of every species is not held in memory
        """
        self.seq_lists = dict()
        q_first = True
        accession = ""
        sequence = ""
        keep = True
        required_files = set(["Species%d.fa" % i for i in species_to_use])
        for d in fasta_file_dir_list:
            for fn in glob.glob(d + "Species*.fa"):
                if os.path.basename(fn) not in required_files: continue
                # each species from the first directory that has it, as
                # file_handler.get_species_fasta_fn finds it
                required_files.discard(os.path.basename(fn))
                with open(fn, 'r', encoding="utf-8") as fasta_file:
                    for line in fasta_file:
                        if line[0] == ">":
                            # deal with old sequence
                            if not q_first:
                                if keep:
                                    self.seq_lists[accession] = sequence
                                sequence = ""
                            q_first = False
                            # get id for new sequence
                            accession = line[1:].rstrip()
                            keep = genes is None or accession in genes
                        elif keep:
                            sequence += line
                    if keep:
                        self.seq_lists[accession] = sequence
    
    def write_seqs_to_fasta(self, seqs, out_filename):
        with open(out_filename, 'w', encoding="utf-8", newline="\n") as out_file:
            for seq in self.sort_seqs([s.to_string() for s in seqs]):
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % seq)
                    out_file.write(self.seq_lists[seq])
                else:
                    print(("ERROR: %s not found" % seq))
                                
    def write_seqs_to_fasta_with_new_accessions(self, seqs, out_filename, id_dict):
        try:
            seqs = [s.to_string() for s in seqs]
        except:
            pass
        with open(out_filename, 'w', encoding="utf-8", newline="\n") as out_file:
            for seq in self.sort_seqs(seqs):
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % id_dict[seq])
                    out_file.write(self.seq_lists[seq])
                    
    def sort_seqs(self, seqs):
        return sorted(seqs, key=lambda x: list(map(int, x.split("_"))))

# write_orthogroup_fasta_files: the sequences held in memory at a time (bytes of the
# species' FASTA files); more are sorted into ranges of orthogroups in temporary files
OG_FASTA_BUFFER_BYTES = 64 * 1024 * 1024


def _gene_numbers(gene):
    """(species, sequence) of a gene: a Seq or its "<species>_<sequence>" ID."""
    if isinstance(gene, str):
        a, b = gene.split("_")
        return int(a), int(b)
    return gene.i_sp, gene.i_seq


def _species_fasta_files(fasta_file_dir_list, species_to_use):
    """The species' FASTA files, each from the first directory that has it (as FastaWriter reads them)."""
    required_files = set(["Species%d.fa" % i for i in species_to_use])
    files = []
    for d in fasta_file_dir_list:
        for fn in glob.glob(d + "Species*.fa"):
            if os.path.basename(fn) not in required_files:
                continue
            required_files.discard(os.path.basename(fn))
            files.append(fn)
    return files


# _fasta_records reads a FASTA file in pieces of this many characters
FASTA_READ_CHARS = 16 * 1024 * 1024


def _fasta_records(fn):
    """
    (header, sequence text) of each record of a FASTA file, as FastaWriter keeps them
    (the header without ">" and its trailing spaces, the text its lines up to the next
    header, each with its newline): read in large pieces, each split at the starts of
    its records ("\\n>"), in the same text mode (newlines as FastaWriter reads them).
    Text before the first header is not in any record (FastaWriter would not keep it
    for a gene of this file either; OrthoFinder's species files start with one).
    """
    with open(fn, 'r', encoding="utf-8") as fasta_file:
        data = fasta_file.read(FASTA_READ_CHARS)
        while data and data[0] != ">":
            start = data.find("\n>")
            if start >= 0:
                data = data[start + 1:]
                break
            more = fasta_file.read(FASTA_READ_CHARS)
            data = data[-1:] + more if more else ""
        while data:
            more = fasta_file.read(FASTA_READ_CHARS)
            pieces = data.split("\n>")
            pieces[0] = pieces[0][1:]                  # (its ">")
            if more:
                last = pieces.pop()                    # (it may go on in the next piece)
                if last.endswith("\n") and more[0] == ">":
                    pieces.append(last[:-1])           # (complete: its "\\n>" across the pieces)
                    data = more
                else:
                    data = ">" + last + more
            else:
                data = ""
            for k, piece in enumerate(pieces):
                header, newline, text = piece.partition("\n")
                if newline and (more or k < len(pieces) - 1):
                    text += "\n"                     # (the newline before the next ">")
                yield header.rstrip(), text


def _gene_id(header):
    """
    (species, sequence) of a header that is a gene's ID "<species>_<sequence>" exactly,
    else None (e.g. not "01_2", which FastaWriter would not find for a gene either).
    """
    sp_text, _, seq_text = header.partition("_")
    if not (sp_text.isdigit() and seq_text.isdigit()) or \
            (len(sp_text) > 1 and sp_text[0] == "0") or (len(seq_text) > 1 and seq_text[0] == "0"):
        return None
    return int(sp_text), int(seq_text)


def _file_records(fn, og_of, extra):
    """(orthogroup k, (species, sequence), text) of the orthogroups' genes in a species' FASTA file, in its order."""
    # (the genes of SpeciesN.fa are "N_<sequence>": checked by that start)
    file_sp = int(os.path.basename(fn)[len("Species"):-len(".fa")])
    prefix = "%d_" % file_sp
    n_prefix = len(prefix)
    file_arr = og_of.get(file_sp)
    for header, text in _fasta_records(fn):
        seq_text = header[n_prefix:]
        if header.startswith(prefix) and seq_text.isdigit() and (seq_text[0] != "0" or len(seq_text) == 1):
            sp, seq, arr = file_sp, int(seq_text), file_arr
        else:
            numbers = _gene_id(header)
            if numbers is None:
                continue
            sp, seq = numbers
            arr = og_of.get(sp)
        if arr is None or seq >= len(arr):
            continue
        k = arr[seq]
        if k < 0:
            continue
        yield int(k), (sp, seq), text
        if extra:
            for k in extra.get((sp, seq), ()):
                yield k, (sp, seq), text


def _pack_fn(d):
    """
    The temporary file of this process's blocks of sequences (one a process, whatever
    the species' files it reads and the ranges: not a file a range and species' file,
    which with 1000s of species could be 100,000s of files, and file systems, above
    all of clusters, suffer many files).
    """
    return os.path.join(d, "pack_%d" % os.getpid())


# The data of the processes that sort the species' files into ranges of orthogroups
# (set once in each, by _set_route: too large to send with each file)
_ROUTE = None


def _set_route(route_fn):
    """
    The routing data of this process, from the file the main process wrote it to: not
    sent through the pipe that starts the process (MBs to GBs), where the main
    process would wait, without end, for a process that died as it started.
    """
    global _ROUTE
    with open(route_fn, "rb") as f:
        (og_of_fn, slices), extra, range_of, d, limit = pickle.load(f)
    # (the species' arrays, one file memory-mapped, read only: one copy, shared by the
    # processes through the file cache, not one each: 4 bytes a gene, ~80 MB with 20M
    # genes; each species' a view of its part)
    all_og_of = np.load(og_of_fn, mmap_mode="r")
    og_of = {sp: all_og_of[start:end] for sp, (start, end) in slices.items()}
    _ROUTE = (og_of, extra, range_of, d, limit)


def _route_file(fn, file_index, route=None):
    """
    Sort the orthogroups' sequences of one species' FASTA file by range of orthogroups,
    in blocks appended to this process's pack (_pack_fn), open while the file is read
    (one pack a process open, not opened for each block: 100s a file at scale).
    Returns the blocks: (range, file_index, block number, pack, offset, length), so
    that a range's are read back in the order of the files.
    route: (og_of, extra, range_of, d, limit), else that of the process (_set_route).
    """
    og_of, extra, range_of, d, limit = route if route is not None else _ROUTE
    pack = _pack_fn(d)
    blocks = []
    pending, held = {}, {}
    with open(pack, "ab") as pack_file:
        pack_file.seek(0, os.SEEK_END)

        def flush(rr):
            data = pickle.dumps(pending.pop(rr), pickle.HIGHEST_PROTOCOL)
            held.pop(rr)
            offset = pack_file.tell()
            pack_file.write(data)
            blocks.append((rr, file_index, len(blocks), pack, offset, len(data)))
        for k, gene, text in _file_records(fn, og_of, extra):
            rr = int(range_of[k])
            block = pending.get(rr)
            if block is None:
                block = pending[rr] = []
                held[rr] = 0
            block.append((k, gene, text))
            held[rr] += len(text) + 32
            if held[rr] >= limit:
                flush(rr)
        for rr in list(pending):
            flush(rr)
    return blocks


def _write_orthogroup_files(found_by_og, orthogroups, names, q_report_missing):
    """
    The FASTA file of each orthogroup (k, filename, its genes' (species, sequence)) from
    found_by_og {k: {(species, sequence): text}}: its genes sorted by species and
    sequence (as FastaWriter.sort_seqs), headed by their IDs or names[ID]. Returns the
    genes without a sequence, as (k, ID), if q_report_missing.
    """
    missing = []
    for k, filename, genes in orthogroups:
        found = found_by_og.get(k, {})
        parts = []
        for sp, seq in sorted(genes):
            text = found.get((sp, seq))
            gene = "%d_%d" % (sp, seq)
            if text is not None:
                parts.append(">%s\n" % (gene if names is None else names[gene]))
                parts.append(text)
            elif q_report_missing:
                missing.append((k, gene))
        with open(filename, 'w', encoding="utf-8", newline="\n") as out_file:
            out_file.write("".join(parts))
    return missing


def _write_range(blocks, orthogroups, names, q_report_missing):
    """
    The FASTA files of a range of orthogroups (_write_orthogroup_files), from its blocks
    (_route_file), taken in the order of the species' files, each read and decoded in
    turn (not all the range's serialized blocks held first): each pack opened once,
    at most one a process.
    """
    found_by_og = {}
    with contextlib.ExitStack() as stack:
        packs = {}
        for _, file_index, n, pack, offset, length in sorted(blocks, key=lambda b: (b[1], b[2])):
            f = packs.get(pack)
            if f is None:
                f = packs[pack] = stack.enter_context(open(pack, "rb"))
            f.seek(offset)
            for k, gene, text in pickle.loads(f.read(length)):
                found_by_og.setdefault(k, {})[gene] = text   # (a gene twice in the files: the last, as FastaWriter)
    return _write_orthogroup_files(found_by_og, orthogroups, names, q_report_missing)


def write_orthogroup_fasta_files(fasta_file_dir_list, species_to_use, ogs, filename_of, names=None,
                                 q_report_missing=True, buffer_bytes=None, tmp_dir=None, n_processes=1):
    """
    The FASTA file of each orthogroup (filename_of(i), unless it is there already,
    e.g. a restart): its genes' sequences, sorted by species and sequence, headed by
    their IDs or names[ID]; a gene without a sequence reported (q_report_missing) or
    left out. The same files as FastaWriter.write_seqs_to_fasta[_with_new_accessions]
    writes, without every sequence of every species in memory (GBs with 1000s of
    species): each gene's orthogroup in an array per species (4 bytes a gene), the
    species' FASTA files read once, the sequences in memory if at most buffer_bytes
    (OG_FASTA_BUFFER_BYTES), else sorted into ranges of orthogroups in temporary
    files (in tmp_dir), each then read on its own, and every orthogroup's file
    written once.
    n_processes: the species' files sorted, and the ranges written, that many at a
    time (new processes, sharing buffer_bytes); the files are the same.
    """
    if buffer_bytes is None:
        buffer_bytes = OG_FASTA_BUFFER_BYTES
    todo = [i for i in range(len(ogs)) if not os.path.exists(filename_of(i))]
    if not todo:
        return
    # each gene's orthogroup (its position in todo): an array per species
    def numbers_of(og):
        # (Seqs, as the orthogroup sets have them: their numbers without a call per gene)
        if og and not isinstance(og[0], str):
            return [(g.i_sp, g.i_seq) for g in og]
        return [_gene_numbers(g) for g in og]
    max_seq = {}
    for i in todo:
        for sp, seq in numbers_of(ogs[i]):
            if seq > max_seq.get(sp, -1):
                max_seq[sp] = seq
    # (int32: 4 bytes a gene, and each process sorting the files has a copy; int64
    # if there are too many orthogroups for it)
    og_dtype = np.int32 if len(todo) < 2 ** 31 else np.int64
    og_of = {sp: np.full(n + 1, -1, dtype=og_dtype) for sp, n in max_seq.items()}
    extra = {}          # (sp, seq) -> the other orthogroups of a gene in more than one
    for k, i in enumerate(todo):
        for sp, seq in numbers_of(ogs[i]):
            arr = og_of[sp]
            if arr[seq] < 0:
                arr[seq] = k
            elif arr[seq] != k:
                extra.setdefault((sp, seq), []).append(k)
    files = _species_fasta_files(fasta_file_dir_list, species_to_use)
    n_processes = max(1, int(n_processes or 1))
    # ranges of orthogroups (of about the same number of genes), each of whose sequences
    # fit in buffer_bytes (shared by the processes): as many as the species' files need
    total_bytes = sum(os.path.getsize(fn) for fn in files)
    range_bytes = max(1, buffer_bytes // n_processes)
    n_ranges = max(1, -(-total_bytes // range_bytes))
    n_genes = [len(ogs[i]) for i in todo]
    per_range = max(1, -(-sum(n_genes) // n_ranges))
    range_of = np.zeros(len(todo), dtype=og_dtype)
    r, in_range = 0, 0
    for k, n in enumerate(n_genes):
        if in_range and in_range + n > per_range:
            r, in_range = r + 1, 0
        range_of[k] = r
        in_range += n
    n_ranges = int(range_of[-1]) + 1
    ks_of_range = [[] for _ in range(n_ranges)]
    for k in range(len(todo)):
        ks_of_range[int(range_of[k])].append(k)

    def orthogroups_of(rr):
        return [(k, filename_of(todo[k]), numbers_of(ogs[todo[k]])) for k in ks_of_range[rr]]

    def names_of(orthogroups):
        if names is None or n_processes == 1:
            return names
        return {"%d_%d" % g: names["%d_%d" % g] for _, _, genes in orthogroups for g in genes if "%d_%d" % g in names}

    if n_ranges == 1:
        found_by_og = {}
        for fn in files:
            for k, gene, text in _file_records(fn, og_of, extra):
                found_by_og.setdefault(k, {})[gene] = text   # (a gene twice in the files: the last, as FastaWriter)
        missing = _write_orthogroup_files(found_by_og, orthogroups_of(0), names, q_report_missing)
    else:
        # the records of each range of orthogroups to blocks in temporary packs (one a
        # process), then each range read and its orthogroups' files written
        with tempfile.TemporaryDirectory(prefix="og_fasta_", dir=tmp_dir) as d:
            route = (og_of, extra, range_of, d, max(1, range_bytes // n_ranges))
            blocks_of_range = [[] for _ in range(n_ranges)]

            def add_blocks(blocks):
                for block in blocks:
                    blocks_of_range[block[0]].append(block)
            if n_processes == 1:
                for j, fn in enumerate(files):
                    add_blocks(_route_file(fn, j, route))
                missing = []
                for rr in range(n_ranges):
                    missing += _write_range(blocks_of_range[rr], orthogroups_of(rr), names, q_report_missing)
                    blocks_of_range[rr] = None
            else:
                import concurrent.futures
                import multiprocessing
                # (new processes, not forks: OrthoFinder may have threads running)
                context = multiprocessing.get_context("spawn")
                route_fn = os.path.join(d, "route.pickle")
                # (all species' arrays in one file, each one's part: one file, not one a species;
                # each copied into it, not first into one concatenated array)
                og_of_fn, slices, start = os.path.join(d, "og_of.npy"), {}, 0
                for sp, arr in og_of.items():
                    slices[sp] = (start, start + len(arr))
                    start += len(arr)
                all_og_of = np.lib.format.open_memmap(og_of_fn, mode="w+", dtype=og_dtype, shape=(start,))
                for sp, arr in og_of.items():
                    all_og_of[slices[sp][0]:slices[sp][1]] = arr
                all_og_of.flush()
                del all_og_of
                with open(route_fn, "wb") as f:
                    pickle.dump(((og_of_fn, slices),) + route[1:], f, pickle.HIGHEST_PROTOCOL)
                with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes, mp_context=context,
                                                            initializer=_set_route, initargs=(route_fn,)) as executor:
                    for future in [executor.submit(_route_file, fn, j) for j, fn in enumerate(files)]:
                        add_blocks(future.result())
                    # (at most 2 ranges a process waiting: the ranges' orthogroups are not all
                    # prepared, and held, at once)
                    missing, in_flight = [], collections.deque()
                    for rr in range(n_ranges):
                        if len(in_flight) >= 2 * n_processes:
                            missing += in_flight.popleft().result()
                        orthogroups = orthogroups_of(rr)
                        in_flight.append(executor.submit(_write_range, blocks_of_range[rr], orthogroups,
                                                         names_of(orthogroups), q_report_missing))
                        blocks_of_range[rr] = None
                    while in_flight:
                        missing += in_flight.popleft().result()
    # (as FastaWriter printed them, in the order of the orthogroups' files)
    for _, gene in sorted(missing, key=lambda m: m[0]) if n_processes > 1 else missing:
        print(("ERROR: %s not found" % gene))


def write_test_file(d):
    test_fn = d + "SimpleTest.fa"
    with open(test_fn, 'w', encoding="utf-8", newline="\n") as outfile:
        outfile.write(">a\nA\n>b\nA")
    return test_fn
 
""" 
-----------------------------------------------------------------------------
                    Identify orthogroups for species tree            
-----------------------------------------------------------------------------
"""

    # N = m.shape[1]
    # f = 1./N
    # fractions = []
    # nOrtho = []
    # for n in range(N):
    #     F = 1.-n*f
    #     fractions.append(F)
    #     nOrtho.append(len(single_copy_with_probability_test(F-1e-5, m)))
    # nOrtho = list(map(float, nOrtho))
    # return fractions, nOrtho
    
def trim_fn(fn, msa_min_seq):
    trim.main(fn, fn, msa_min_seq, 0.1, 500, 0.75, False)

""" 
-----------------------------------------------------------------------------
                             TreesForOrthogroups            
-----------------------------------------------------------------------------
"""    
     
class SpeciesTreeFirstProgress(object):
    """
    The progress of the alignments and trees, shown as two progress bars: one
    for the orthogroups for the species tree, then one for the rest. They run
    as one batch: the second bar starts with those of the rest already done.
    n_rest includes the tree of the concatenated alignment, added later.
    """

    def __init__(self, species_tree_units, n_rest):
        self.waiting_for = {id(unit) for unit in species_tree_units}
        self.species_tree_bar = command_runner.Progress(len(species_tree_units))
        self.rest_bar = None
        self.n_rest = n_rest
        self.n_rest_done = 0

    def start(self):
        self.species_tree_bar.start()

    def done(self, cmd_unit):
        if id(cmd_unit) in self.waiting_for:
            self.waiting_for.discard(id(cmd_unit))
            self.species_tree_bar.done(cmd_unit)
            if not self.waiting_for:
                self.species_tree_bar.stop()
        elif self.rest_bar is None:
            self.n_rest_done += 1
        else:
            self.rest_bar.done(cmd_unit)

    def species_tree_done(self):
        return not self.waiting_for

    def started_rest(self):
        return self.rest_bar is not None

    def start_rest(self):
        self.rest_bar = command_runner.Progress(self.n_rest)
        self.rest_bar.start()
        self.rest_bar.advance(self.n_rest_done)

    def added(self, n):
        pass    # the tree of the concatenated alignment is in n_rest already

    def stop(self):
        self.species_tree_bar.stop()
        if self.rest_bar is not None:
            self.rest_bar.stop()


class TreesForOrthogroups(object):
    def __init__(self, program_caller, msa_program, tree_program):
        self.program_caller = program_caller
        self.msa_program = msa_program
        self.tree_program = tree_program
    
    def get_fasta_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_seq_fn(i_og, q_results)
            
    def get_alignment_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_align_fn(i_og, q_results)
            
    def get_tree_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_tree_fn(i_og, q_results)
        
    # def write_fasta_files(self, fastaWriter, ogs, idDict, qBoth):
    #     # The results ones are now written by default after orthogroups, check they're not already there
    #     for iOg, og in enumerate(ogs):
    #         if not os.path.exists(self.get_fasta_filename(iOg, True)):
    #             fastaWriter.write_seqs_to_fasta_with_new_accessions(og, self.get_fasta_filename(iOg, True), idDict)
    #     if qBoth: 
    #         for iOg, og in enumerate(ogs):
    #             if not os.path.exists(self.get_fasta_filename(iOg)):
    #                 fastaWriter.write_seqs_to_fasta(og, self.get_fasta_filename(iOg))

    def write_fasta_files(self, fasta_file_dir_list, species_to_use, ogs, id_dict, q_id=True, q_results=False,
                          n_processes=1):
        """
        The FASTA file of each orthogroup (not those there already): in the results with
        the gene names (id_dict), and/or in the working directory with the IDs, from the
        species' FASTA files without holding all their sequences (write_orthogroup_fasta_files).
        """
        tmp_dir = file_handler.handler.get_working_directory_write()
        if q_results:
            write_orthogroup_fasta_files(fasta_file_dir_list, species_to_use, ogs,
                                         lambda i: self.get_fasta_filename(i, True), names=id_dict,
                                         q_report_missing=False, tmp_dir=tmp_dir, n_processes=n_processes)
        if q_id:
            write_orthogroup_fasta_files(fasta_file_dir_list, species_to_use, ogs, self.get_fasta_filename,
                                         tmp_dir=tmp_dir, n_processes=n_processes)

    def get_alignment_commands_and_new_filenames(self, ogs, i_og_restart=0, method_threads=None):
        iogs_align = [i for i, og in enumerate(ogs) if len(og) >= 2 and i >= i_og_restart]
        infn_list = [self.get_fasta_filename(i) for i in iogs_align]
        outfn_list = [self.get_alignment_filename(i) for i in iogs_align]
        id_list = ["OG%07d" % i for i in iogs_align]
        n_seqs = [len(ogs[i]) for i in iogs_align]
        return self.program_caller.get_msa_commands(self.msa_program, infn_list, outfn_list, id_list, n_seqs, method_threads=method_threads), iogs_align
        
    def get_tree_commands(self, alignments_for_tree, iogs_align, ogs, method_threads=None):
        # RAxML needs at least 4 sequences (the other programs 3)
        n_min = RAXML_MIN_SEQS if self.tree_program in ("raxml", "raxml-ng") else 3
        iogs_tree = [i for i in iogs_align if len(ogs[i]) >= n_min]
        alignments_for_tree = [alignments_for_tree[i] for i, iog in enumerate(iogs_align) if len(ogs[iog]) >= n_min]  # Remove those to align with tree filenames
        outfn_list = [self.get_tree_filename(i) for i in iogs_tree]
        id_list = ["OG%07d" % i for i in iogs_tree]
        n_seqs = [len(ogs[i]) for i in iogs_tree]
        return self.program_caller.get_tree_commands(self.tree_program, alignments_for_tree, outfn_list, id_list, n_seqs, method_threads=method_threads), iogs_tree
     
    def rename_alignment_taxa(self, ids_align_fns, acc_align_fns, ids_dict):
        for i, (align_fn, out_align_fn) in enumerate(zip(ids_align_fns, acc_align_fns)):
            if not os.path.exists(align_fn): continue
            with open(align_fn, 'r', encoding="utf-8") as infile, open(out_align_fn, 'w', encoding="utf-8", newline="\n") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(">" + ids_dict[line[1:].rstrip()] + "\n")
                    else:
                        outfile.write(line)

          
    def do_trees(
            self, 
            og_set, 
            id_dict, 
            species_id_dict, 
            species_to_use, 
            n_processes, 
            q_stop_after_seqs, 
            q_stop_after_alignments, 
            q_do_species_tree, 
            q_trim, 
            i_og_restart=0, 
            cmd_order="descending",
            method_threads=None,
            species_tree_program=None,
            n_skip=50,
            msa_min_seq=10,
            species_tree_alignment=False,
        ):
        """
        species_tree_program: the species tree from the gene trees with this
            program (-ST), instead of from the concatenated alignment
        species_tree_alignment: with species_tree_program, also write the
            concatenated alignment (it is always written without it)
        """
        q_from_gene_trees = species_tree_program is not None

        print_on_error = True
        id_dict.update(species_id_dict) # same code will then also convert concatenated alignment for species tree
        # 0       
        results_dirs_full_path = [
            file_handler.handler.get_results_seqs_dir(), 
            file_handler.handler.get_results_align_dir(), 
            file_handler.handler.get_results_trees_dir()
        ]
        
        # 1. (the orthogroups' FASTA files are written before: not all the sequences read again here)
        ogs = og_set.get_ogs_all()
        
        ## -------- Writing to WD/Sequence_ids is not is mcl_orthogroups/infer_orthogroups.py ---------------
        # self.write_fasta_files(fastaWriter, ogs, idDict, True)

        if q_stop_after_seqs: return results_dirs_full_path

        # 3
        # Get OGs to use for species tree
        if q_do_species_tree:
            occupancy, iogs4 = og_set.orthogroup_occupancy()
            i_ogs_for_species_tree, f_single_copy = determine_orthogroups_for_species_tree(occupancy, iogs4)
            concatenated_algn_fn = file_handler.handler.get_species_tree_concat_align_fn()
            if not i_ogs_for_species_tree and not q_from_gene_trees:
                # (before the alignments: else the run stops when they are done, maybe hours later)
                print("ERROR: No orthogroup is single-copy in enough species to infer the species tree from a "
                      "concatenated alignment of single-copy orthogroups (the default). This can happen with "
                      "many species, or many duplicated genes. Please either give the species tree (-s), or "
                      "infer it from the gene trees instead (-ST astral-pro).")
                util.fail()
        else:
            i_ogs_for_species_tree = []
        orig_align_commands_and_filenames, orig_iogs_align = \
            self.get_alignment_commands_and_new_filenames(
                ogs, 
                i_og_restart, 
                method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
            )

        # (references: neither list is changed, the originals are replaced by reversed copies)
        align_commands_and_filenames = orig_align_commands_and_filenames
        iogs_align = orig_iogs_align
        aln_tasksize = None
        if cmd_order == "ascending":
            # aln_tasksize = [len(og) for og in aln_filtered_ogs[::-1]]
            orig_align_commands_and_filenames = orig_align_commands_and_filenames[::-1]
            orig_iogs_align = orig_iogs_align[::-1]
            aln_tasksize = [
                os.stat(item).st_size
                for fn, _ in orig_align_commands_and_filenames
                for item in fn.split() if os.path.isfile(item) and os.path.exists(item)
            ]
        if q_stop_after_alignments:
            printing.print_underline("Inferring multiple sequence alignments")
            command_runner.run_parallel_commands_and_move_results_file(
                n_processes, 
                orig_align_commands_and_filenames, 
                False, 
                tasksize=aln_tasksize,
                method_threads=method_threads,
                cmd_order=cmd_order,
                q_trim=q_trim,
                q_print_on_error=print_on_error,
            )
            if q_do_species_tree:
                if not q_from_gene_trees:
                    create_concatenated_alignment(
                        i_ogs_for_species_tree, 
                        ogs, 
                        self.get_alignment_filename, 
                        concatenated_algn_fn, 
                        f_single_copy
                    )
                    # write OGs used to file
                    # in the working directory; copied to Species_Tree/ with the species tree
                    with open(file_handler.handler.get_concat_orthogroups_fn(), 'w', encoding="utf-8", newline="\n") as outfile:
                        for iog in i_ogs_for_species_tree: outfile.write("OG%07d\n" % iog)
            
            # ------------------ this section is not needed at this stage for the new procedure -------------
            # ids -> accessions
            alignment_files_to_use = [self.get_alignment_filename(i) for i, in orig_iogs_align]
            accession_alignment_fns = [self.get_alignment_filename(i, True) for i, in orig_iogs_align]
            if q_do_species_tree: 
                alignment_files_to_use.append(concatenated_algn_fn)
                accession_alignment_fns.append(file_handler.handler.get_species_tree_concat_align_fn(True))
            self.rename_alignment_taxa(alignment_files_to_use, accession_alignment_fns, id_dict)
            
            return results_dirs_full_path[:2]
        
        # Otherwise, alignments and trees, as one batch:
        # 1. Alignments (and trees) of the orthogroups for the species tree, first
        # 2. The other alignments and trees, largest orthogroups first
        # 3. When the orthogroups for the species tree are done, the concatenated
        #    alignment is created and its tree is started ahead of those still queued
        alignment_files_to_use = [self.get_alignment_filename(i) for i in iogs_align]

        orig_tree_commands_and_filenames, orig_iogs_tree = \
            self.get_tree_commands(
                alignment_files_to_use, 
                iogs_align, 
                ogs,
                method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
        )

        iog_to_align = dict(zip(iogs_align, align_commands_and_filenames))
        iog_to_tree = dict(zip(orig_iogs_tree, orig_tree_commands_and_filenames))
        iog_to_alignment_fn = dict(zip(iogs_align, alignment_files_to_use))

        def commands_for_og(iog):
            commands = [iog_to_align[iog]]
            if iog in iog_to_tree:
                if q_trim:
                    commands.append((trim_fn, (iog_to_alignment_fn[iog], msa_min_seq)))
                commands.append(iog_to_tree[iog])
            return commands

        i_ogs_for_species_tree_set = set(i_ogs_for_species_tree)
        species_tree_units = [commands_for_og(i) for i in i_ogs_for_species_tree if i in iog_to_align]
        commands_and_filenames = species_tree_units + \
            [commands_for_og(i) for i in iogs_align if i not in i_ogs_for_species_tree_set]
        if cmd_order == "ascending":
            commands_and_filenames = commands_and_filenames[::-1]

        follow_up = None
        progress = None
        if q_do_species_tree:
            if not q_from_gene_trees:
                print(("Species tree: Using %d orthogroups with minimum of %0.1f%% of species having single-copy genes in any orthogroup" % (len(i_ogs_for_species_tree), 100.*f_single_copy)))
            printing.print_underline("Inferring multiple sequence alignments for species tree") 
            species_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()

            def write_concatenated_alignment(q_progress=False):
                create_concatenated_alignment(i_ogs_for_species_tree, ogs, self.get_alignment_filename, concatenated_algn_fn, f_single_copy,
                                              q_progress=q_progress)
                # write OGs used to file
                # in the working directory; copied to Species_Tree/ with the species tree
                with open(file_handler.handler.get_concat_orthogroups_fn(), 'w', encoding="utf-8", newline="\n") as outfile:
                    for iog in i_ogs_for_species_tree: outfile.write("OG%07d\n" % iog)

            def species_tree_commands():
                """The commands for the tree of the concatenated alignment, which is created first (none with -ST)."""
                commands = []
                if not q_from_gene_trees:
                    write_concatenated_alignment()
                    util.print_time("Scheduling species tree inference from concatenated alignment")
                    commands = [
                        self.program_caller.get_tree_commands(
                                self.tree_program, 
                                [concatenated_algn_fn], 
                                [species_tree_fn_ids], 
                                ["SpeciesTree"],
                                method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
                            )
                        ]
                printing.print_underline("Inferring remaining multiple sequence alignments and gene trees") 
                return commands

            if species_tree_units:
                n_rest = len(commands_and_filenames) - len(species_tree_units) + (0 if q_from_gene_trees else 1)
                progress = SpeciesTreeFirstProgress(species_tree_units, n_rest)

                def follow_up(unit):
                    # Once, when the orthogroups for the species tree are all done
                    if not progress.species_tree_done() or progress.started_rest():
                        return []
                    commands = species_tree_commands()
                    progress.start_rest()
                    return commands
            else:
                commands_and_filenames = species_tree_commands() + commands_and_filenames
        else:
            printing.print_underline("Inferring multiple sequence alignments and gene trees")
        command_runner.run_parallel_commands_and_move_results_file(n_processes, 
                                                 commands_and_filenames, 
                                                 True, 
                                                 method_threads=method_threads,
                                                 cmd_order=cmd_order,
                                                 q_trim=q_trim,
                                                 q_print_on_error=print_on_error,
                                                 follow_up=follow_up,
                                                 progress=progress,
                                                 )
        
        # # Convert ids to accessions for MSA
        # accessionAlignmentFNs = [self.get_alignment_filename(i, True) for i in iogs_align]
        # # Add concatenated Alignment
        # if qDoSpeciesTree:
        #     alignmentFilesToUse.append(concatenated_algn_fn)
        #     accessionAlignmentFNs.append(files.FileHandler.get_species_tree_concat_align_fn(True))

        #     qHaveSupport = taxa.have_support_values(speciesTreeFN_ids)
        #     if os.path.exists(speciesTreeFN_ids):
        #         taxa.rename_tree_taxa(speciesTreeFN_ids, files.FileHandler.get_species_tree_unrooted_fn(True), idDict, qSupport=qHaveSupport, qFixNegatives=True)
        #     else:
        #         text = "ERROR: Species tree inference failed"
        #         files.FileHandler.log_fail_and_exit(text)


        # Add concatenated Alignment
        if q_do_species_tree:
            if q_from_gene_trees and species_tree_alignment:
                # Not needed for the species tree: after the alignments and trees
                printing.print_underline("Creating concatenated alignment for species tree")
                print(("Species tree: Using %d orthogroups with minimum of %0.1f%% of species having single-copy genes in any orthogroup" % (len(i_ogs_for_species_tree), 100.*f_single_copy)))
                # a separate step here (after the gene trees): with a progress bar
                write_concatenated_alignment(q_progress=True)
            if q_from_gene_trees:
                printing.print_underline("Inferring unrooted species tree using %s" % species_tree_program)
                species_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
                from_gene_trees.infer_species_tree(
                    self.program_caller, species_tree_program, file_handler.handler.get_ogs_tree_dir(),
                    file_handler.handler.get_core_astral_filename(), species_tree_fn_ids, n_processes, n_skip=n_skip)

            q_have_support = taxa.have_support_values(species_tree_fn_ids)
            if os.path.exists(species_tree_fn_ids):
                taxa.rename_tree_taxa(species_tree_fn_ids, file_handler.handler.get_species_tree_unrooted_fn(True), id_dict, q_support=q_have_support, q_fix_negatives=True)
                util.print_time("Done species tree inference")
            else:
                text = "ERROR: Species tree inference failed"
                file_handler.handler.log_fail_and_exit(text)
        
        return results_dirs_full_path[:2]
