# -*- coding: utf-8 -*-

import os
import time
import sys
import csv
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

from ..utils import util, files, parallel_task_manager


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
                    writer = csv.writer(outfile, delimiter="\t")
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
                        writer = csv.writer(outfile, delimiter="\t")
                        writer.writerow(("Orthogroup", sp0_name, sp1_name))

        InitialiseSuspectGenesDirs(
            self.nspecies,
            self.speciesToUse,
            self.speciesDict
        )

        with open(self.duplications_path, util.csv_write_mode) as outfile:
            util.writerow(
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
                for j in range(self.nspecies):
                    if i == j:
                        continue
                    text = olog_lines[i][j]
                    if text:
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

        text = "".join(util.getrow(row) for row in rows)

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

class OrderedHogCommitter(object):
    """
    Commit HOG rows in deterministic OG order.

    Only cached_hogs are buffered. Full result payloads are not buffered here.
    """

    def __init__(self, hog_writer, iogs_ordered):
        self.hog_writer = hog_writer
        self.iogs_ordered = list(sorted(iogs_ordered))
        self.pending_hogs = {}
        self.next_index = 0

    def add_result(self, iog, cached_hogs):
        self.pending_hogs[iog] = cached_hogs or []
        self._drain_ready()

    def add_skip(self, iog):
        self.pending_hogs[iog] = []
        self._drain_ready()

    def _drain_ready(self):
        while self.next_index < len(self.iogs_ordered):
            next_iog = self.iogs_ordered[self.next_index]

            if next_iog not in self.pending_hogs:
                break

            cached_hogs = self.pending_hogs.pop(next_iog)

            if cached_hogs:
                self.hog_writer.WriteCachedHOGs(cached_hogs, lock_hogs=None)

            self.next_index += 1

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
                for j in range(self.nspecies):
                    if i == j:
                        continue

                    text = olog_lines[i][j]
                    if text:
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
            writer1 = csv.writer(outfile, delimiter="\t")
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
            frac,
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
        util.writerow(dups_file_handle, [og_name, spIDs[sp_node_id] if q_terminal else sp_node_id, gene_node_name, frac, isSTRIDE, gene_list0, gene_list1]) 




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
            for j, text in enumerate(row):
                if i == j or not text:
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
    ):
    """Commit HOG rows in deterministic OG order."""
    committer = OrderedHogCommitter(hog_writer, iogs_ordered)
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

        committer.assert_finished()
        hog_writer.file_cache.flush_all()
        hog_counts = dict(hog_writer.iHOG)

    except Exception:
        error_text = traceback.format_exc()

    try:
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
    """Analyse OGs and route HOG and non-HOG output independently."""
    n_writers = len(non_hog_queues)

    while True:
        try:
            iog = args_queue.get(True, 0.1)

            if iog is None:
                break

            result = tree_analyser.AnalyseTree(iog)

            if result is None:
                hog_queue.put(("skip", iog))
                progress_queue.put(("skip", iog))
                continue

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

            progress_queue.put((
                "result",
                iog,
                result["n_orthologues"],
            ))

        except queue.Empty:
            continue

        except Exception:
            tb = traceback.format_exc()
            try:
                hog_queue.put(("error", tb))
            except Exception:
                pass
            for writer_queue in non_hog_queues:
                try:
                    writer_queue.put(("error", tb))
                except Exception:
                    pass
            progress_queue.put(("error", None, tb))
            break

    hog_queue.put(None)
    for writer_queue in non_hog_queues:
        writer_queue.put(None)
    progress_queue.put(None)


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
        compatibility_mode=False,
        write_hog_tree=False,
        fix_files=False,
        fd_limit=None,
        GRACE_PERIOD=10.0,
        STALL_TIMEOUT=120.0,
        writer_queue_size=None,
        n_writer_processes=None,
    ):
    """
    Parallel tree analysis with two output paths.

    HOG rows are committed in deterministic OG order. Other output files are
    written immediately by exclusive file owners and can be sorted afterward.
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

    progressbar, task = util.get_progressbar(total_tasks)
    progressbar.start()

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

    nOrthologues_SpPair = util.nOrtho_sp(nspecies)
    completed_tasks = 0
    skipped_tasks = 0
    active_workers = nProcesses
    fatal = False
    hog_done = False
    non_hog_done_ids = set()
    hog_counts = None
    last_progress_time = time.time()
    last_writer_activity_time = time.time()

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
                if (
                    completed_tasks < total_tasks
                    and time.time() - last_progress_time > STALL_TIMEOUT
                ):
                    print(
                        "ERROR: Stalled for %ss "
                        "(completed %d/%d, active_workers=%d)." % (
                            STALL_TIMEOUT,
                            completed_tasks,
                            total_tasks,
                            active_workers,
                        )
                    )
                    fatal = True
                    break

            elif msg is None:
                active_workers -= 1

            elif isinstance(msg, tuple) and msg[0] == "error":
                print("ERROR: worker error:")
                print(msg[2])
                fatal = True
                break

            elif isinstance(msg, tuple) and msg[0] == "skip":
                skipped_tasks += 1
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.time()

            elif isinstance(msg, tuple) and msg[0] == "result":
                _, iog, nOrtho = msg
                nOrthologues_SpPair += nOrtho
                completed_tasks += 1
                progressbar.update(task, advance=1)
                last_progress_time = time.time()

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

                last_writer_activity_time = time.time()

                if isinstance(wmsg, tuple) and wmsg[0] == "hog_done":
                    _, hog_counts = wmsg
                    hog_done = True

                elif isinstance(wmsg, tuple) and wmsg[0] == "non_hog_done":
                    _, writer_id = wmsg
                    non_hog_done_ids.add(writer_id)

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "hog_error",
                    "hog_close_error",
                    "non_hog_error",
                    "non_hog_close_error",
                }:
                    print("ERROR: writer error:")
                    print(wmsg[-1])
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

            if (
                completed_tasks >= total_tasks
                and active_workers == 0
                and (
                    not hog_done
                    or len(non_hog_done_ids) < n_writer_processes
                )
                and time.time() - last_writer_activity_time > STALL_TIMEOUT
            ):
                print("ERROR: output writers stalled after analysis completed.")
                fatal = True
                break

    finally:
        for proc in runningProcesses:
            proc.join(timeout=GRACE_PERIOD)
        for proc in runningProcesses:
            if proc.is_alive():
                proc.terminate()
        for proc in runningProcesses:
            proc.join()

        hog_proc.join(timeout=GRACE_PERIOD)
        if hog_proc.is_alive():
            hog_proc.terminate()
        hog_proc.join()

        for proc in non_hog_procs:
            proc.join(timeout=GRACE_PERIOD)
        for proc in non_hog_procs:
            if proc.is_alive():
                proc.terminate()
        for proc in non_hog_procs:
            proc.join()

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

    if fatal:
        util.Fail()

    if not hog_done:
        print("ERROR: ordered HOG writer did not finish cleanly.")
        util.Fail()

    if len(non_hog_done_ids) != n_writer_processes:
        print(
            "ERROR: only %d/%d non-HOG writers finished cleanly." %
            (len(non_hog_done_ids), n_writer_processes)
        )
        util.Fail()

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

    args_queue = mp.Queue()

    for fn, gz in fns:
        if os.path.exists(fn):
            args_queue.put((fn, gz))

    parallel_task_manager.RunMethodParallel(
        SortFileByFirstColumnNoRepair,
        args_queue,
        n_parallel
    )

    suspect_queue = mp.Queue()
    dSuspectGenes = files.FileHandler.GetSuspectGenesDir()
    for sp in species:
        fn = os.path.join(dSuspectGenes, "%s.txt" % sp)
        if os.path.exists(fn):
            suspect_queue.put((fn,))

    parallel_task_manager.RunMethodParallel(
        SortPlainTextFile,
        suspect_queue,
        n_parallel
    )


def SortFileByFirstColumnNoRepair(fn, gz=False):
    """
    Sort a TSV file by first column.

    This is for orthologues, xenologues, and duplications only.
    It must never be used for HOG files.
    """
    with util.file_open(fn, util.csv_read_mode, gz=gz) as infile:
        header = next(infile, None)
        if header is None:
            return

        lines = list(infile)

    if not lines:
        return

    lines.sort(key=lambda s: (s.split("\t", 1)[0], s))

    with util.file_open(fn, util.csv_write_mode, gz=gz) as outfile:
        outfile.write(header)
        outfile.write("".join(lines))


def SortPlainTextFile(fn):
    """Sort a plain-text output file deterministically."""
    with open(fn, util.csv_read_mode) as infile:
        lines = infile.readlines()

    if not lines:
        return

    lines.sort()

    with open(fn, util.csv_write_mode) as outfile:
        outfile.writelines(lines)



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
