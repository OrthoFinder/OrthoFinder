# -*- coding: utf-8 -*-

import os
import time
import gzip
import heapq
import pickle
import tempfile
import sys
import resource
import traceback
import warnings
from collections import OrderedDict

import queue
import multiprocessing as mp

try:
    from rich import print
except ImportError:
    ...

from ..utils import util, files, parallel_task_manager, file_io
from .tree_processor import IterOlogRow


class LazyFileCache(object):
    """
    Parent-process file cache.

    Keeps at most max_open file handles open. When another file is needed,
    the least-recently-used handle is closed.

    """

    def __init__(self, max_open=64):
        self.max_open = int(max_open)
        self.handles = OrderedDict()

    def get(self, path, mode, gz=False):
        key = (path, bool(gz))

        if key in self.handles:
            fh = self.handles.pop(key)
            self.handles[key] = fh
            return fh

        if len(self.handles) >= self.max_open:
            _, old_fh = self.handles.popitem(last=False)
            old_fh.close()

        fh = util.file_open(path, mode, gz=gz)
        self.handles[key] = fh
        return fh

    def write(self, path, text, mode=None, gz=False):
        if not text:
            return
        if mode is None:
            mode = util.csv_append_mode
        fh = self.get(path, mode, gz=gz)
        fh.write(text)

    def flush_all(self):
        for fh in self.handles.values():
            fh.flush()

    def close_all(self):
        for fh in self.handles.values():
            fh.close()
        self.handles.clear()

class ParentOutputWriter(object):
    """
    Owns all shared output files.

    Workers must not write orthologues, xenologues, duplications, suspect genes,
    or HOG files. They return data. This parent writer writes the data.
    """

    def __init__(
            self,
            dResultsOrthologues,
            speciesDict,
            speciesToUse,
            SequenceDict,
            spec_seq_dict,
            stride_dups,
            hog_writer,
            fewer_open_files,
            save_space,
            write_hog_tree,
            fix_files,
            max_open=64,
            flush_every_results=1000,
            flush_every_seconds=60,
        ):
        self.dResultsOrthologues = dResultsOrthologues
        self.speciesDict = speciesDict
        self.speciesToUse = list(speciesToUse)
        self.SequenceDict = SequenceDict
        self.spec_seq_dict = spec_seq_dict
        self.stride_dups = stride_dups
        self.hog_writer = hog_writer
        self.nspecies = len(self.speciesToUse)

        self.fewer_open_files = fewer_open_files or save_space
        self.save_space = save_space

        self.write_hog_tree = write_hog_tree
        self.fix_files = fix_files
        self.need_olog_output = (not write_hog_tree or not fix_files)

        self.cache = LazyFileCache(max_open=max_open)

        self.flush_every_results = int(flush_every_results)
        self.flush_every_seconds = float(flush_every_seconds)
        self.results_since_flush = 0
        self.last_flush_time = time.time()

        self.putative_xenolog_dir = files.FileHandler.GetPutativeXenelogsDir()
        self.suspect_genes_dir = files.FileHandler.GetSuspectGenesDir()
        self.duplications_path = files.FileHandler.GetDuplicationsFN()

        if self.need_olog_output:
            self.initialise_orthologue_outputs()

    def initialise_orthologue_outputs(self):
        for i in range(self.nspecies):
            sp0 = str(self.speciesToUse[i])
            sp0_name = self.speciesDict[sp0]

            if self.fewer_open_files:
                filename = os.path.join(
                    self.dResultsOrthologues,
                    "%s.tsv" % sp0_name
                )
                with util.file_open(filename, util.csv_write_mode, gz=self.save_space) as outfile:
                    writer = file_io.writer(outfile)
                    writer.writerow((
                        "Orthogroup",
                        "Species",
                        sp0_name,
                        "Orthologs"
                    ))
            else:
                d = os.path.join(
                    self.dResultsOrthologues,
                    "Orthologues_" + sp0_name
                )
                if not os.path.exists(d):
                    os.mkdir(d)

                for j in range(self.nspecies):
                    if j == i:
                        continue

                    sp1 = str(self.speciesToUse[j])
                    sp1_name = self.speciesDict[sp1]

                    fn = os.path.join(
                        d,
                        "%s__v__%s.tsv" % (sp0_name, sp1_name)
                    )

                    with open(fn, util.csv_write_mode) as outfile:
                        writer = file_io.writer(outfile)
                        writer.writerow(("Orthogroup", sp0_name, sp1_name))

        InitialiseSuspectGenesDirs(
            self.nspecies,
            self.speciesToUse,
            self.speciesDict
        )

        with open(self.duplications_path, util.csv_write_mode) as outfile:
            file_io.write_unquoted(
                outfile,
                [
                    "Orthogroup",
                    "Species Tree Node",
                    "Gene Tree Node",
                    "Support",
                    "Type",
                    "Genes 1",
                    "Genes 2"
                ]
            )

    def ortholog_path(self, i, j=None):
        sp0 = str(self.speciesToUse[i])
        sp0_name = self.speciesDict[sp0]

        if self.fewer_open_files:
            return os.path.join(
                self.dResultsOrthologues,
                "%s.tsv" % sp0_name
            )

        sp1 = str(self.speciesToUse[j])
        sp1_name = self.speciesDict[sp1]

        d = os.path.join(
            self.dResultsOrthologues,
            "Orthologues_" + sp0_name
        )

        return os.path.join(
            d,
            "%s__v__%s.tsv" % (sp0_name, sp1_name)
        )

    def xenolog_path(self, i):
        sp0 = str(self.speciesToUse[i])
        sp0_name = self.speciesDict[sp0]
        return os.path.join(self.putative_xenolog_dir, "%s.tsv" % sp0_name)

    def write_ortholog_lines(self, olog_lines):
        if self.fewer_open_files:
            for i in range(self.nspecies):
                text = olog_lines[i][0]
                if text:
                    self.cache.write(
                        self.ortholog_path(i),
                        text,
                        mode=util.csv_append_mode,
                        gz=self.save_space
                    )
        else:
            for i in range(self.nspecies):
                for j, text in IterOlogRow(olog_lines[i]):
                    if i == j:
                        continue
                    self.cache.write(
                        self.ortholog_path(i, j),
                        text,
                        mode=util.csv_append_mode,
                        gz=False
                    )

    def write_xenolog_lines(self, olog_sus_lines):
        for i in range(self.nspecies):
            text = olog_sus_lines[i]
            if text:
                self.cache.write(
                    self.xenolog_path(i),
                    text,
                    mode=util.csv_append_mode,
                    gz=False
                )

    def write_duplications(self, iog, duplications):
        if not duplications:
            return

        og_name = "OG%07d" % iog
        rows = DuplicationRows(
            og_name,
            duplications,
            self.speciesDict,
            self.spec_seq_dict,
            self.stride_dups
        )

        text = "".join(file_io.unquoted_line(row) for row in rows)

        self.cache.write(
            self.duplications_path,
            text,
            mode=util.csv_append_mode,
            gz=False
        )

    def write_suspect_genes(self, suspect_genes):
        if not suspect_genes:
            return

        species = list(map(str, self.speciesToUse))

        for i in range(self.nspecies):
            strsp0 = species[i]
            strsp0_ = strsp0 + "_"

            these_genes = [
                g for g in suspect_genes
                if g.startswith(strsp0_)
            ]

            if not these_genes:
                continue

            path = os.path.join(
                self.suspect_genes_dir,
                self.speciesDict[strsp0] + ".txt"
            )

            text = "\n".join(
                [self.SequenceDict[g] for g in these_genes]
            ) + "\n"

            self.cache.write(
                path,
                text,
                mode=util.csv_append_mode,
                gz=False
            )

    def write_hog_rows(self, cached_hogs):
        if cached_hogs:
            self.hog_writer.WriteCachedHOGs(cached_hogs, lock_hogs=None)


    def write_non_hog_result(self, result, do_flush=True):
        """
        Write non-HOG output for one analysed OG.

        These rows do not assign HOG IDs, so they can be streamed immediately.
        """
        if result is None:
            return

        if self.need_olog_output:
            self.write_duplications(
                result["iog"],
                result.get("duplications", [])
            )
            self.write_suspect_genes(
                result.get("suspect_genes", set())
            )
            self.write_ortholog_lines(
                result.get("olog_lines", [])
            )
            self.write_xenolog_lines(
                result.get("olog_sus_lines", [])
            )

        if do_flush:
            self.maybe_flush()


    def write_result(self, result):
        """
        Serial/safe ordered writer path.

        This keeps the serial parent-writer behaviour for n_parallel == 1 or for
        fallback parent-ordered mode.
        """
        if result is None:
            return

        self.write_hog_rows(result.get("cached_hogs", []))
        self.write_non_hog_result(result, do_flush=False)
        self.maybe_flush()

    def flush(self):
        self.cache.flush_all()

        if hasattr(self.hog_writer, "file_cache"):
            self.hog_writer.file_cache.flush_all()


    def maybe_flush(self):
        self.results_since_flush += 1

        now = time.time()

        q_by_count = (
            self.flush_every_results > 0
            and self.results_since_flush >= self.flush_every_results
        )

        q_by_time = (
            self.flush_every_seconds > 0
            and now - self.last_flush_time >= self.flush_every_seconds
        )

        if q_by_count or q_by_time:
            self.flush()
            self.results_since_flush = 0
            self.last_flush_time = now

    def close(self):
        """
        Final flush and close.

        close_all() would flush anyway, but explicit flush makes the intent clear.
        """
        try:
            self.flush()
        finally:
            self.cache.close_all()
            self.hog_writer.close_files()

# Out-of-order HOG batches held in memory (pickled) before spilling to disk.
HOG_PENDING_MEMORY_BYTES = 1024 * 1024 * 1024


class OrderedHogCommitter(object):
    """
    Commit HOG rows in deterministic OG order.

    Only cached_hogs are buffered. Full result payloads are not buffered here.
    Out-of-order batches are kept pickled in memory up to max_pending_bytes
    and spilled to a temporary file beyond that.
    """

    def __init__(
            self,
            hog_writer,
            iogs_ordered,
            max_pending_bytes=None,
            spill_dir=None,
        ):
        self.hog_writer = hog_writer
        self.iogs_ordered = list(sorted(iogs_ordered))
        # iog -> ("mem", pickled bytes) or ("disk", offset, length)
        self.pending_hogs = {}
        self.next_index = 0
        if max_pending_bytes is None:
            max_pending_bytes = HOG_PENDING_MEMORY_BYTES
        self.max_pending_bytes = max_pending_bytes
        self.pending_bytes = 0
        self.spill_dir = spill_dir
        self.spill_file = None

    def add_result(self, iog, cached_hogs):
        cached_hogs = cached_hogs or []
        if (
            self.next_index < len(self.iogs_ordered)
            and iog == self.iogs_ordered[self.next_index]
        ):
            # The common case: nothing to buffer.
            self._commit(cached_hogs)
            self.next_index += 1
        else:
            self._buffer(iog, cached_hogs)
        self._drain_ready()

    def add_skip(self, iog):
        self.add_result(iog, [])

    def _buffer(self, iog, cached_hogs):
        # Pickled bytes are far smaller than the live row objects. Past the
        # memory budget they go to a temporary file, so a slow early OG cannot
        # make the buffer grow without bound.
        data = pickle.dumps(cached_hogs, pickle.HIGHEST_PROTOCOL)
        if self.pending_bytes + len(data) <= self.max_pending_bytes:
            self.pending_hogs[iog] = ("mem", data)
            self.pending_bytes += len(data)
            return
        if self.spill_file is None:
            self.spill_file = tempfile.TemporaryFile(
                prefix="orthofinder_hog_spill_", dir=self.spill_dir
            )
        self.spill_file.seek(0, os.SEEK_END)
        offset = self.spill_file.tell()
        self.spill_file.write(data)
        self.pending_hogs[iog] = ("disk", offset, len(data))

    def _load(self, entry):
        if entry[0] == "mem":
            self.pending_bytes -= len(entry[1])
            return pickle.loads(entry[1])
        _, offset, length = entry
        self.spill_file.seek(offset)
        return pickle.loads(self.spill_file.read(length))

    def _commit(self, cached_hogs):
        if cached_hogs:
            self.hog_writer.WriteCachedHOGs(cached_hogs, lock_hogs=None)

    def _drain_ready(self):
        while self.next_index < len(self.iogs_ordered):
            next_iog = self.iogs_ordered[self.next_index]

            if next_iog not in self.pending_hogs:
                break

            self._commit(self._load(self.pending_hogs.pop(next_iog)))
            self.next_index += 1

    def close(self):
        if self.spill_file is not None:
            self.spill_file.close()
            self.spill_file = None

    def assert_finished(self):
        if self.pending_hogs:
            raise RuntimeError(
                "OrderedHogCommitter finished with %d pending HOG batches."
                % len(self.pending_hogs)
            )

        if self.next_index != len(self.iogs_ordered):
            raise RuntimeError(
                "OrderedHogCommitter stopped early: committed %d/%d OGs."
                % (self.next_index, len(self.iogs_ordered))
            )


class NonHogAppendWriter(object):
    """
    Append-only writer for non-HOG output.

    This does not initialise/truncate files.
    It is intended for small 2/3-gene orthogroups after the main
    tree-based orthologue-writing step has already created the files.
    """

    def __init__(
            self,
            dResultsOrthologues,
            speciesDict,
            speciesToUse,
            save_space=False,
            fewer_open_files=False,
            max_open=64,
        ):
        self.dResultsOrthologues = dResultsOrthologues
        self.speciesDict = speciesDict
        self.speciesToUse = list(speciesToUse)
        self.nspecies = len(self.speciesToUse)

        self.save_space = save_space
        self.fewer_open_files = fewer_open_files or save_space

        self.cache = LazyFileCache(max_open=max_open)

    def ortholog_path(self, i, j=None):
        sp0 = str(self.speciesToUse[i])
        sp0_name = self.speciesDict[sp0]

        if self.fewer_open_files:
            return os.path.join(
                self.dResultsOrthologues,
                "%s.tsv" % sp0_name
            )

        sp1 = str(self.speciesToUse[j])
        sp1_name = self.speciesDict[sp1]

        d = os.path.join(
            self.dResultsOrthologues,
            "Orthologues_" + sp0_name
        )

        return os.path.join(
            d,
            "%s__v__%s.tsv" % (sp0_name, sp1_name)
        )

    def write_ortholog_lines(self, olog_lines):
        if self.fewer_open_files:
            for i in range(self.nspecies):
                text = olog_lines[i][0]
                if text:
                    self.cache.write(
                        self.ortholog_path(i),
                        text,
                        mode=util.csv_append_mode,
                        gz=self.save_space
                    )
        else:
            for i in range(self.nspecies):
                for j, text in IterOlogRow(olog_lines[i]):
                    if i == j:
                        continue
                    self.cache.write(
                        self.ortholog_path(i, j),
                        text,
                        mode=util.csv_append_mode,
                        gz=False
                    )

    def flush(self):
        self.cache.flush_all()

    def close(self):
        self.cache.close_all()


def InitialiseSuspectGenesDirs(nspecies, speciesIDs, speciesDict):
    files.FileHandler.GetSuspectGenesDir()  # creates the directory
    dSuspectOrthologues = files.FileHandler.GetPutativeXenelogsDir()
    for index1 in range(nspecies):
        with open(dSuspectOrthologues + '%s.tsv' % speciesDict[str(speciesIDs[index1])], util.csv_write_mode) as outfile:
            writer1 = file_io.writer(outfile)
            writer1.writerow(("Orthogroup", speciesDict[str(speciesIDs[index1])], "Other"))

def WriteSuspectGenes(nspecies, speciesToUse, suspect_genes, speciesDict, SequenceDict):
    species = list(map(str, speciesToUse))
    dSuspectGenes = files.FileHandler.GetSuspectGenesDir()
    for index0 in range(nspecies):
        strsp0 = species[index0]
        strsp0_ = strsp0+"_"
        these_genes = [g for g in suspect_genes if g.startswith(strsp0_)]
        if len(these_genes) > 0:
            with open(dSuspectGenes + speciesDict[strsp0] + ".txt", util.csv_append_mode) as outfile:
                # not a CSV file so \n line endings are fine
                outfile.write("\n".join([SequenceDict[g] for g in these_genes]) + "\n")


def DuplicationRows(og_name, duplications, spIDs, seqIDs, stride_dups):
    """
    Convert duplication records into rows.

    This converts duplication records for the parent writer.
    """
    rows = []

    for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
        q_terminal = not sp_node_id.startswith("N")

        if stride_dups is None:
            isSTRIDE = "Terminal" if q_terminal else "Non-Terminal"
        else:
            if q_terminal:
                isSTRIDE = "Terminal"
            elif frozenset(genes0 + genes1) in stride_dups:
                isSTRIDE = "Non-Terminal: STRIDE"
            else:
                isSTRIDE = "Non-Terminal"

        gene_list0 = ", ".join([seqIDs[g] for g in genes0])
        gene_list1 = ", ".join([seqIDs[g] for g in genes1])

        rows.append([
            og_name,
            spIDs[sp_node_id] if q_terminal else sp_node_id,
            gene_node_name,
            str(frac),   # as text: all-str rows take the fast path
            isSTRIDE,
            gene_list0,
            gene_list1
        ])

    return rows


def WriteDuplications(dups_file_handle, og_name, duplications, spIDs, seqIDs, stride_dups):
    """
    Args:
        duplications - list of (sp_node_id, gene_node_name, fraction, genes0, genes1)
    """
    for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
        q_terminal = not sp_node_id.startswith("N")
        if stride_dups is None:
            isSTRIDE = "Terminal" if q_terminal else "Non-Terminal"
        else:
            isSTRIDE = "Terminal" if q_terminal else "Non-Terminal: STRIDE" if frozenset(genes0 + genes1) in stride_dups else "Non-Terminal"
        gene_list0 = ", ".join([seqIDs[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
        gene_list1 = ", ".join([seqIDs[g] for g in genes1])
        file_io.write_unquoted(dups_file_handle, [og_name, spIDs[sp_node_id] if q_terminal else sp_node_id, gene_node_name, str(frac), isSTRIDE, gene_list0, gene_list1]) 




def get_n_writer_processes(
        nspecies,
        n_processes,
        fewer_open_files=False,
        max_writers=12,
    ):
    """Choose the number of non-HOG output writer processes."""
    nspecies = max(1, int(nspecies))
    n_processes = max(1, int(n_processes))

    species_per_writer = 96 if fewer_open_files else 48
    n_writers = (nspecies + species_per_writer - 1) // species_per_writer
    writer_cpu_cap = max(1, n_processes // 2)

    return max(
        1,
        min(
            int(n_writers),
            int(max_writers),
            writer_cpu_cap,
        )
    )


def _empty_non_hog_payload(iog):
    return {
        "iog": iog,
        "olog_chunks": [],
        "xenolog_chunks": [],
        "suspect_chunks": [],
        "duplications": [],
    }


def PartitionNonHogResult(
        result,
        n_writers,
        tree_analyser,
        fewer_open_files,
        need_olog_output,
    ):
    """
    Split one analysed OG into file-owner payloads for non-HOG output.

    Orthologue, xenologue, and suspect-gene files are owned by source species.
    Duplications.tsv is owned by writer 0. Empty writer payloads are omitted.
    """
    if not need_olog_output:
        return {}

    iog = result["iog"]
    payloads = {}

    def payload_for(owner):
        payload = payloads.get(owner)
        if payload is None:
            payload = _empty_non_hog_payload(iog)
            payloads[owner] = payload
        return payload

    olog_lines = result.get("olog_lines", [])

    if fewer_open_files:
        for i, row in enumerate(olog_lines):
            if not row:
                continue
            text = row[0]
            if not text:
                continue
            owner = i % n_writers
            payload_for(owner)["olog_chunks"].append((i, None, text))
    else:
        for i, row in enumerate(olog_lines):
            owner = i % n_writers
            for j, text in IterOlogRow(row):
                if i == j:
                    continue
                payload_for(owner)["olog_chunks"].append((i, j, text))

    for i, text in enumerate(result.get("olog_sus_lines", [])):
        if not text:
            continue
        owner = i % n_writers
        payload_for(owner)["xenolog_chunks"].append((i, text))

    suspect_genes = result.get("suspect_genes", set())
    if suspect_genes:
        sp_id_to_index = {
            str(sp): i
            for i, sp in enumerate(tree_analyser.speciesToUse)
        }
        by_species = {}

        for g in suspect_genes:
            sp_id = g.split("_", 1)[0]
            i = sp_id_to_index.get(sp_id)
            if i is None:
                continue
            by_species.setdefault(i, []).append(g)

        for i, genes in by_species.items():
            genes.sort()
            text = "\n".join(
                tree_analyser.SequenceDict[g]
                for g in genes
            ) + "\n"
            owner = i % n_writers
            payload_for(owner)["suspect_chunks"].append((i, text))

    duplications = result.get("duplications", [])
    if duplications:
        payload_for(0)["duplications"] = duplications

    return payloads


def _write_non_hog_payload(output_writer, payload):
    for i, j, text in payload.get("olog_chunks", []):
        output_writer.cache.write(
            output_writer.ortholog_path(i, j),
            text,
            mode=util.csv_append_mode,
            gz=output_writer.save_space if output_writer.fewer_open_files else False,
        )

    for i, text in payload.get("xenolog_chunks", []):
        output_writer.cache.write(
            output_writer.xenolog_path(i),
            text,
            mode=util.csv_append_mode,
            gz=False,
        )

    for i, text in payload.get("suspect_chunks", []):
        sp_id = str(output_writer.speciesToUse[i])
        path = os.path.join(
            output_writer.suspect_genes_dir,
            output_writer.speciesDict[sp_id] + ".txt"
        )
        output_writer.cache.write(
            path,
            text,
            mode=util.csv_append_mode,
            gz=False,
        )

    duplications = payload.get("duplications", [])
    if duplications:
        output_writer.write_duplications(
            payload["iog"],
            duplications,
        )

    output_writer.maybe_flush()


def NonHogWriterProcess(
        writer_queue,
        writer_status_queue,
        output_writer,
        n_workers,
        writer_id,
    ):
    """Write only files exclusively owned by this non-HOG writer."""
    active_workers = n_workers
    error_text = None

    try:
        while active_workers > 0:
            msg = writer_queue.get()

            if msg is None:
                active_workers -= 1
                continue

            if not isinstance(msg, tuple) or len(msg) < 2:
                raise TypeError(
                    "Unexpected non-HOG writer message: %s %r" %
                    (type(msg), msg)
                )

            kind = msg[0]

            if kind == "result":
                _write_non_hog_payload(output_writer, msg[1])
                writer_status_queue.put(("non_hog_progress", writer_id))
            elif kind == "error":
                raise RuntimeError(msg[-1])
            else:
                raise TypeError(
                    "Unexpected non-HOG writer message kind: %r" % kind
                )

        output_writer.cache.flush_all()

    except Exception:
        error_text = traceback.format_exc()

    try:
        output_writer.cache.close_all()
    except Exception:
        close_error = traceback.format_exc()
        error_text = close_error if error_text is None else error_text + "\n" + close_error

    if error_text is None:
        writer_status_queue.put(("non_hog_done", writer_id))
    else:
        writer_status_queue.put((
            "non_hog_error",
            writer_id,
            error_text,
        ))


def OrderedHogWriterProcess(
        hog_queue,
        writer_status_queue,
        hog_writer,
        iogs_ordered,
        n_workers,
        spill_dir=None,
    ):
    """Commit HOG rows in deterministic OG order."""
    committer = OrderedHogCommitter(hog_writer, iogs_ordered, spill_dir=spill_dir)
    active_workers = n_workers
    error_text = None
    hog_counts = None

    try:
        while active_workers > 0:
            msg = hog_queue.get()

            if msg is None:
                active_workers -= 1
                continue

            if not isinstance(msg, tuple) or len(msg) < 2:
                raise TypeError(
                    "Unexpected HOG writer message: %s %r" %
                    (type(msg), msg)
                )

            kind = msg[0]

            if kind == "result":
                _, iog, cached_hogs = msg
                committer.add_result(iog, cached_hogs)
            elif kind == "skip":
                _, iog = msg
                committer.add_skip(iog)
            elif kind == "error":
                raise RuntimeError(msg[-1])
            else:
                raise TypeError(
                    "Unexpected HOG writer message kind: %r" % kind
                )

            writer_status_queue.put(("hog_progress",))

        committer.assert_finished()
        hog_writer.file_cache.flush_all()
        hog_counts = dict(hog_writer.iHOG)

    except Exception:
        error_text = traceback.format_exc()

    try:
        committer.close()
        hog_writer.close_files()
    except Exception:
        close_error = traceback.format_exc()
        error_text = close_error if error_text is None else error_text + "\n" + close_error

    if error_text is None:
        writer_status_queue.put(("hog_done", hog_counts))
    else:
        writer_status_queue.put((
            "hog_error",
            error_text,
        ))


def Worker_RunOrthologsMethod_Pipeline(
        tree_analyser,
        nspecies,
        args_queue,
        hog_queue,
        non_hog_queues,
        progress_queue,
        fewer_open_files,
        need_olog_output,
        n_ologs_cache=100,
        write_hog_tree=False,
        fix_files=False,
    ):
    """
    Analyse OGs and route HOG and non-HOG output independently.

    Orthologue counts are accumulated locally and sent once with
    "worker_done". Sending a dense nspecies x nspecies count object per OG
    is ~40 MB of pickled data per OG at 1,000 species.
    """
    n_writers = len(non_hog_queues)
    worker_pid = os.getpid()
    nOrtho_acc = util.nOrtho_sp(nspecies)

    while True:
        try:
            iog = args_queue.get(True, 0.1)

            if iog is None:
                break

            progress_queue.put(("start", worker_pid, iog))
            result = tree_analyser.AnalyseTree(iog, nOrtho_acc=nOrtho_acc)

            if result is None:
                hog_queue.put(("skip", iog))
                progress_queue.put(("skip", iog))
                continue

            # Analysis is finished. Any further wait is on the output queues.
            progress_queue.put(("output", iog))

            hog_queue.put((
                "result",
                iog,
                result.get("cached_hogs", []),
            ))

            payloads = PartitionNonHogResult(
                result,
                n_writers,
                tree_analyser,
                fewer_open_files,
                need_olog_output,
            )

            for owner, payload in payloads.items():
                non_hog_queues[owner].put(("result", payload))

            progress_queue.put(("result", iog))

        except queue.Empty:
            continue

        except Exception:
            tb = traceback.format_exc()
            # Report directly to the parent. A failed writer may have left a
            # full output queue, so even sending it an error could block here.
            progress_queue.put(("error", None, tb))
            return

    hog_queue.put(None)
    for writer_queue in non_hog_queues:
        writer_queue.put(None)
    progress_queue.put(("worker_done", worker_pid, nOrtho_acc))


def RunOrthologsParallel_Pipeline(
        tree_analyser,
        nspecies,
        args_queue,
        nProcesses,
        total_tasks,
        fewer_open_files,
        output_writer,
        iogs_ordered,
        n_ologs_cache=100,
        write_hog_tree=False,
        fix_files=False,
        fd_limit=None,
        GRACE_PERIOD=10.0,
        STALL_TIMEOUT=120.0,
        writer_queue_size=None,
        n_writer_processes=None,
        spill_dir=None,
    ):
    """
    Parallel tree analysis with two output paths.

    HOG rows are committed in deterministic OG order. Other output files are
    written immediately by exclusive file owners and can be sorted afterward.

    STALL_TIMEOUT is the interval between long-wait warnings, not a deadline
    for a tree or an output batch. Silence alone cannot distinguish slow work
    from a hang. Child process failures are checked independently.
    """
    if fd_limit is not None:
        if sys.platform.startswith("linux") or sys.platform == "darwin":
            set_file_descriptor_limit(fd_limit)
        else:
            warnings.warn(
                "File descriptor limit adjustment is not supported on %s." %
                sys.platform
            )

    if n_writer_processes is None:
        n_writer_processes = get_n_writer_processes(
            nspecies,
            nProcesses,
            fewer_open_files=fewer_open_files,
        )
    else:
        n_writer_processes = max(
            1,
            min(int(n_writer_processes), max(1, nProcesses))
        )

    # print(
    #     "Output pipeline: %d analysis workers, 1 ordered HOG writer, "
    #     "%d non-HOG writer process%s (%d species, %s mode)." % (
    #         nProcesses,
    #         n_writer_processes,
    #         "" if n_writer_processes == 1 else "es",
    #         nspecies,
    #         "compact" if fewer_open_files else "pairwise",
    #     )
    # )

    if writer_queue_size is None:
        writer_queue_size = max(2 * nProcesses, 16)

    hog_queue = mp.Queue(maxsize=max(4 * nProcesses, 32))
    non_hog_queues = [
        mp.Queue(maxsize=writer_queue_size)
        for _ in range(n_writer_processes)
    ]
    progress_queue = mp.Queue(maxsize=max(4 * nProcesses, 32))
    writer_status_queue = mp.Queue()

    for _ in range(nProcesses):
        args_queue.put(None)

    hog_proc = mp.Process(
        target=OrderedHogWriterProcess,
        args=(
            hog_queue,
            writer_status_queue,
            output_writer.hog_writer,
            iogs_ordered,
            nProcesses,
            spill_dir,
        )
    )

    non_hog_procs = [
        mp.Process(
            target=NonHogWriterProcess,
            args=(
                non_hog_queues[writer_id],
                writer_status_queue,
                output_writer,
                nProcesses,
                writer_id,
            )
        )
        for writer_id in range(n_writer_processes)
    ]

    runningProcesses = [
        mp.Process(
            target=Worker_RunOrthologsMethod_Pipeline,
            args=(
                tree_analyser,
                nspecies,
                args_queue,
                hog_queue,
                non_hog_queues,
                progress_queue,
                fewer_open_files,
                output_writer.need_olog_output,
                n_ologs_cache,
                write_hog_tree,
                fix_files,
            )
        )
        for _ in range(nProcesses)
    ]

    hog_proc.start()
    for proc in non_hog_procs:
        proc.start()
    for proc in runningProcesses:
        proc.start()

    # Start the progress bar only after all children exist: its refresh thread
    # holds the console lock while drawing, which a forked child could inherit.
    progressbar, task = util.get_progressbar(total_tasks)
    progressbar.start()

    nOrthologues_SpPair = util.nOrtho_sp(nspecies)
    errors = []   # error text for the WorkerError raised at the end

    def report_error(text):
        print(text)
        errors.append(text)

    completed_tasks = 0
    skipped_tasks = 0
    active_workers = nProcesses
    finished_worker_pids = set()
    in_flight = {}
    exited_without_status = {}
    fatal = False
    hog_done = False
    non_hog_done_ids = set()
    hog_counts = None
    last_progress_time = time.monotonic()
    last_writer_activity_time = None
    child_processes = runningProcesses + [hog_proc] + non_hog_procs

    try:
        while (
            completed_tasks < total_tasks
            or active_workers > 0
            or not hog_done
            or len(non_hog_done_ids) < n_writer_processes
        ):
            try:
                msg = progress_queue.get(timeout=0.1)
            except queue.Empty:
                msg = "__EMPTY__"

            if msg == "__EMPTY__":
                pass

            elif isinstance(msg, tuple) and msg[0] == "worker_done":
                _, worker_pid, worker_counts = msg
                finished_worker_pids.add(worker_pid)
                active_workers = nProcesses - len(finished_worker_pids)
                nOrthologues_SpPair += worker_counts

            elif isinstance(msg, tuple) and msg[0] == "start":
                _, worker_pid, iog = msg
                in_flight[iog] = (worker_pid, time.monotonic(), "analysing")

            elif isinstance(msg, tuple) and msg[0] == "output":
                iog = msg[1]
                if iog in in_flight:
                    worker_pid, started, _ = in_flight[iog]
                    in_flight[iog] = (worker_pid, started, "waiting for output queue")

            elif isinstance(msg, tuple) and msg[0] == "error":
                report_error("ERROR: worker error:\n%s" % msg[2])
                fatal = True
                break

            elif isinstance(msg, tuple) and msg[0] == "skip":
                in_flight.pop(msg[1], None)
                skipped_tasks += 1
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.monotonic()

            elif isinstance(msg, tuple) and msg[0] == "result":
                iog = msg[1]
                in_flight.pop(iog, None)
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.monotonic()

            else:
                fatal = True
                raise TypeError(
                    "Unexpected progress message: %s %r" %
                    (type(msg), msg)
                )

            while True:
                try:
                    wmsg = writer_status_queue.get_nowait()
                except queue.Empty:
                    break

                if last_writer_activity_time is not None:
                    last_writer_activity_time = time.monotonic()

                if isinstance(wmsg, tuple) and wmsg[0] == "hog_done":
                    _, hog_counts = wmsg
                    hog_done = True

                elif isinstance(wmsg, tuple) and wmsg[0] == "non_hog_done":
                    _, writer_id = wmsg
                    non_hog_done_ids.add(writer_id)

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "hog_progress", "non_hog_progress",
                }:
                    pass

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "hog_error",
                    "hog_close_error",
                    "non_hog_error",
                    "non_hog_close_error",
                }:
                    report_error("ERROR: writer error:\n%s" % wmsg[-1])
                    fatal = True
                    break

                else:
                    fatal = True
                    raise TypeError(
                        "Unexpected writer status message: %s %r" %
                        (type(wmsg), wmsg)
                    )

            if fatal:
                break

            now = time.monotonic()
            process_statuses = [
                ("analysis worker", proc, proc.pid in finished_worker_pids)
                for proc in runningProcesses
            ] + [("HOG writer", hog_proc, hog_done)] + [
                ("non-HOG writer %d" % writer_id, proc,
                 writer_id in non_hog_done_ids)
                for writer_id, proc in enumerate(non_hog_procs)
            ]

            for role, proc, reported_done in process_statuses:
                exitcode = proc.exitcode
                if exitcode is None:
                    continue
                if exitcode != 0:
                    report_error("ERROR: %s (pid=%d) exited with code %d." %
                          (role, proc.pid, exitcode))
                    fatal = True
                    break
                if not reported_done:
                    # Queue messages can still be waiting in the parent when
                    # a child exits normally. Drain them before starting the
                    # grace period for a missing completion message.
                    if role == "analysis worker" and msg != "__EMPTY__":
                        continue
                    first_seen = exited_without_status.setdefault(proc.pid, now)
                    if now - first_seen > GRACE_PERIOD:
                        report_error("ERROR: %s (pid=%d) exited without reporting completion." %
                              (role, proc.pid))
                        fatal = True
                        break

            if fatal:
                break

            if active_workers == 0 and completed_tasks < total_tasks:
                report_error("ERROR: all analysis workers finished, but only %d/%d "
                      "orthogroups were reported." % (completed_tasks, total_tasks))
                fatal = True
                break

            if (
                (completed_tasks < total_tasks or active_workers > 0)
                and now - last_progress_time > STALL_TIMEOUT
            ):
                pending = ", ".join(
                    "OG%07d (pid=%d, %.0fs, %s)" % (iog, pid, now - started, stage)
                    for iog, (pid, started, stage) in sorted(
                        in_flight.items(), key=lambda item: item[1][1]
                    )[:5]
                )
                if any(
                    stage != "analysing"
                    for _, _, stage in in_flight.values()
                ):
                    pending += (
                        " Workers waiting for the output queue means the "
                        "output writers (disk) are the bottleneck."
                    )
                print("WARNING: Still waiting for analysis workers "
                      "(completed %d/%d, active_workers=%d).%s" % (
                          completed_tasks, total_tasks, active_workers,
                          " Pending: " + pending if pending else "",
                      ))
                last_progress_time = now

            if (
                completed_tasks >= total_tasks
                and active_workers == 0
                and (
                    not hog_done
                    or len(non_hog_done_ids) < n_writer_processes
                )
            ):
                if last_writer_activity_time is None:
                    # Start the drain timer here, not at pipeline startup.
                    last_writer_activity_time = now
                elif now - last_writer_activity_time > STALL_TIMEOUT:
                    print("WARNING: Still waiting for output writers to finish; "
                          "analysis is complete.")
                    last_writer_activity_time = now

    finally:
        if fatal or sys.exc_info()[0] is not None:
            # A failed consumer can leave producers blocked on full queues.
            # Stop all children before joining any of them in the error path.
            parallel_task_manager.TerminateProcesses(child_processes)
            args_queue.cancel_join_thread()
        else:
            # Everything reported completion: wait for the children to exit.
            # They are not stopped for taking long, only warned about.
            try:
                parallel_task_manager.WaitForExit(
                    child_processes, "orthologue pipeline processes", GRACE_PERIOD
                )
            except parallel_task_manager.WorkerError as e:
                report_error("ERROR: %s" % e)
                fatal = True
                parallel_task_manager.TerminateProcesses(child_processes)
        for proc in child_processes:
            proc.join()
            if proc.exitcode != 0:
                if not fatal:
                    report_error("ERROR: child process (pid=%d) exited with code %d." %
                          (proc.pid, proc.exitcode))
                fatal = True

        progressbar.stop()

        for q in [hog_queue] + non_hog_queues + [progress_queue, writer_status_queue]:
            try:
                q.close()
                q.join_thread()
            except Exception:
                pass

    if hog_counts is not None:
        output_writer.hog_writer.iHOG.clear()
        output_writer.hog_writer.iHOG.update(hog_counts)

    skip_rate = skipped_tasks / max(1, total_tasks)
    if skip_rate > 0.02:
        print(
            "WARNING: skipped %d/%d tasks (%.1f%%)." %
            (skipped_tasks, total_tasks, 100.0 * skip_rate)
        )

    if not fatal and not hog_done:
        report_error("ERROR: ordered HOG writer did not finish cleanly.")
        fatal = True

    if not fatal and len(non_hog_done_ids) != n_writer_processes:
        report_error(
            "ERROR: only %d/%d non-HOG writers finished cleanly." %
            (len(non_hog_done_ids), n_writer_processes)
        )
        fatal = True

    if fatal:
        # Raised (not Fail()) so main() records the child's error in checkpoint.txt.
        raise parallel_task_manager.WorkerError(
            "\n".join(errors) or "Orthologue inference failed."
        )

    return nOrthologues_SpPair

def set_file_descriptor_limit(fd_limit) -> None:
    """
    Try to raise the soft open-file limit.

    This is only a convenience. The lazy writer should not depend on this.
    """
    try:
        if isinstance(fd_limit, (tuple, list)):
            requested_soft = int(fd_limit[0])
        else:
            requested_soft = int(fd_limit)

        if requested_soft <= 0:
            print("Ignoring non-positive file descriptor limit.")
            return

        soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)

        print(
            "Current file descriptor limits: soft=%s, hard=%s" %
            (soft_limit, hard_limit)
        )

        if soft_limit >= requested_soft:
            print(
                "File descriptor soft limit is already sufficient: %s" %
                soft_limit
            )
            return

        if hard_limit != resource.RLIM_INFINITY and hard_limit < requested_soft:
            print(
                "Cannot raise soft file descriptor limit to %s; "
                "hard limit is only %s. Increase the hard limit outside "
                "OrthoFinder if needed." %
                (requested_soft, hard_limit)
            )
            return

        resource.setrlimit(
            resource.RLIMIT_NOFILE,
            (requested_soft, hard_limit)
        )

        new_soft, new_hard = resource.getrlimit(resource.RLIMIT_NOFILE)

        print(
            "New file descriptor limits: soft=%s, hard=%s" %
            (new_soft, new_hard)
        )

    except AttributeError:
        print("File descriptor limit functions not available on this platform.")
    except Exception as e:
        print("Could not adjust file descriptor limit: %s" % e)


def SortNonHogOutputFiles(
        n_parallel,
        speciesToUse,
        speciesDict,
        fewer_open_files,
        save_space,
        write_hog_tree,
        fix_files,
    ):

    if write_hog_tree and fix_files:
        return

    species = [speciesDict[str(sp)] for sp in speciesToUse]
    dResultsOrthologues = files.FileHandler.GetOrthologuesDirectory()

    fns = []

    # Orthologue files.
    if fewer_open_files or save_space:
        for sp in species:
            fns.append((
                os.path.join(dResultsOrthologues, "%s.tsv" % sp),
                bool(save_space)
            ))
    else:
        for sp1 in species:
            d = os.path.join(dResultsOrthologues, "Orthologues_" + sp1)
            for sp2 in species:
                if sp1 == sp2:
                    continue
                fns.append((
                    os.path.join(d, "%s__v__%s.tsv" % (sp1, sp2)),
                    False
                ))

    # Xenolog files.
    dXenologs = files.FileHandler.GetPutativeXenelogsDir()
    for sp in species:
        fns.append((
            os.path.join(dXenologs, "%s.tsv" % sp),
            False
        ))

    # Duplications file.
    fns.append((
        files.FileHandler.GetDuplicationsFN(),
        False
    ))

    # Compressed outputs are written as fn + ".gz" (see util.file_open).
    # Sort the largest files first so one big file does not run alone at the end.
    existing = []
    for fn, gz in fns:
        path = fn + ".gz" if gz else fn
        if os.path.exists(path):
            existing.append((os.path.getsize(path), fn, gz))
    existing.sort(key=lambda x: -x[0])

    args_queue = mp.Queue()

    n_sort_tasks = 0
    for _, fn, gz in existing:
        args_queue.put((fn, gz))
        n_sort_tasks += 1

    parallel_task_manager.RunMethodParallel(
        SortFileByFirstColumnNoRepair,
        args_queue,
        n_parallel,
        total_tasks=n_sort_tasks,
        show_progress=False,
    )

    suspect_queue = mp.Queue()
    dSuspectGenes = files.FileHandler.GetSuspectGenesDir()
    suspect_fns = [
        os.path.join(dSuspectGenes, "%s.txt" % sp) for sp in species
    ]
    suspect_fns = sorted(
        (fn for fn in suspect_fns if os.path.exists(fn)),
        key=lambda fn: -os.path.getsize(fn),
    )
    n_suspect_tasks = 0
    for fn in suspect_fns:
        suspect_queue.put((fn,))
        n_suspect_tasks += 1

    parallel_task_manager.RunMethodParallel(
        SortPlainTextFile,
        suspect_queue,
        n_parallel,
        total_tasks=n_suspect_tasks,
        show_progress=False,
    )


# Above this many characters a file is sorted in chunks that are spilled to
# disk and merged, so memory use per sorting process stays bounded.
SORT_CHUNK_CHARS = 128 * 1024 * 1024


def _open_text(path, mode, gz):
    return gzip.open(path, mode) if gz else open(path, mode)


def _open_gz_text_as(path, name):
    """
    A gzipped text file written at path, but recording `name` as the file
    name in its gzip header: a temporary file that is renamed when complete
    keeps the final name when decompressed (gunzip -N, archive managers).
    """
    import io
    raw = open(path, "wb")
    try:
        gz = gzip.GzipFile(filename=name, mode="wb", fileobj=raw)
    except Exception:
        raw.close()
        raise
    gz.myfileobj = raw     # closed with the GzipFile, as gzip.open does
    return io.TextIOWrapper(gz)


def _write_sorted_run(lines, run_path, gz):
    # Runs from a compressed file are compressed too (fast level), so the
    # temporary space stays close to the compressed size, not the full text.
    lines.sort()
    if gz:
        outfile = gzip.open(run_path, util.csv_write_mode, compresslevel=1)
    else:
        outfile = open(run_path, util.csv_write_mode)
    with outfile:
        outfile.writelines(lines)


def SortLinesInFile(fn, gz=False, has_header=False, chunk_chars=SORT_CHUNK_CHARS):
    """
    Sort the lines of fn (fn + ".gz" if gz), keeping an optional header first.

    Small files are sorted in memory. Larger files are split into sorted runs
    next to the file and merged. The result is written to a temporary file
    and then renamed over the original, so an interrupted sort never leaves a
    truncated file.
    """
    path = fn + ".gz" if gz else fn
    tmp_path = path + ".sorting.tmp"
    run_paths = []
    run_files = []
    lines = []
    n_chars = 0
    header = None

    try:
        with _open_text(path, util.csv_read_mode, gz) as infile:
            if has_header:
                header = next(infile, None)
                if header is None:
                    return

            for line in infile:
                if not line.endswith("\n"):
                    line += "\n"
                lines.append(line)
                n_chars += len(line)
                if n_chars >= chunk_chars:
                    run_path = "%s.run%d.tmp" % (path, len(run_paths))
                    run_paths.append(run_path)
                    _write_sorted_run(lines, run_path, gz)
                    lines = []
                    n_chars = 0

        if not lines and not run_paths:
            return

        if run_paths:
            if lines:
                run_path = "%s.run%d.tmp" % (path, len(run_paths))
                run_paths.append(run_path)
                _write_sorted_run(lines, run_path, gz)
                lines = []
            for run_path in run_paths:
                run_files.append(_open_text(run_path, util.csv_read_mode, gz))
            sorted_lines = heapq.merge(*run_files)
        else:
            lines.sort()
            sorted_lines = lines

        # written as path's name (not the temporary name) in the gzip header
        outfile = (_open_gz_text_as(tmp_path, os.path.basename(path)) if gz
                   else open(tmp_path, util.csv_write_mode))
        with outfile:
            if header is not None:
                outfile.write(header)
            outfile.writelines(sorted_lines)

        os.replace(tmp_path, path)

    finally:
        for f in run_files:
            f.close()
        for p in run_paths + [tmp_path]:
            if os.path.exists(p):
                os.remove(p)


def SortFileByFirstColumnNoRepair(fn, gz=False):
    """
    Sort a TSV file by first column.

    This is for orthologues, xenologues, and duplications only.
    It must never be used for HOG files.

    Plain line order is the same as ordering by (first column, line): the
    first column is followed by a tab, which sorts before any character that
    can appear in an OG ID. Sorting whole lines avoids building a key per line.
    """
    SortLinesInFile(fn, gz=gz, has_header=True)


def SortPlainTextFile(fn):
    """Sort a plain-text output file deterministically."""
    SortLinesInFile(fn, gz=False, has_header=False)



def ValidateHogWriterNoDuplicateIds(hog_writer):

    if not getattr(hog_writer, "write_output", True):
        return

    paths = sorted(set(getattr(hog_writer, "hog_paths", {}).values()))

    bad_files = []

    for fn in paths:
        if not os.path.exists(fn):
            continue

        seen = set()
        duplicates = []
        malformed = []

        with open(fn, util.csv_read_mode) as infile:
            header = next(infile, None)

            for line_no, line in enumerate(infile, start=2):
                line = line.rstrip("\n")

                if not line:
                    continue

                parts = line.split("\t", 1)

                if len(parts) < 2:
                    malformed.append((line_no, line))
                    continue

                hog_id = parts[0]

                if ".HOG" not in hog_id:
                    malformed.append((line_no, line))
                    continue

                if hog_id in seen:
                    duplicates.append((line_no, hog_id))
                else:
                    seen.add(hog_id)

        if duplicates or malformed:
            bad_files.append((fn, duplicates[:20], malformed[:20]))

    if bad_files:
        msg = [
            "ERROR: HOG ID validation failed.",
            "No automatic renumbering was performed."
        ]

        for fn, duplicates, malformed in bad_files[:10]:
            msg.append("\nFile: %s" % fn)

            if duplicates:
                msg.append("Duplicate HOG IDs:")
                for line_no, hog_id in duplicates:
                    msg.append("  line %d: %s" % (line_no, hog_id))

            if malformed:
                msg.append("Malformed HOG rows:")
                for line_no, line in malformed:
                    msg.append("  line %d: %s" % (line_no, line[:200]))

        raise RuntimeError("\n".join(msg))
