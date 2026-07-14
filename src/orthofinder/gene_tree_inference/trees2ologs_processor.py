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

        This keeps the old parent-writer behaviour for n_parallel == 1 or for
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


class PipelineOutputWriter(object):

    def __init__(self, output_writer, iogs_ordered):
        self.output_writer = output_writer
        self.hog_committer = OrderedHogCommitter(
            output_writer.hog_writer,
            iogs_ordered
        )

    def write_result(self, result):
        if result is None:
            return

        iog = result["iog"]

        # Stream non-HOG output immediately.
        self.output_writer.write_non_hog_result(result, do_flush=False)

        # Commit HOGs only when OG order allows.
        self.hog_committer.add_result(
            iog,
            result.get("cached_hogs", [])
        )

        self.output_writer.maybe_flush()

    def write_skip(self, iog):
        self.hog_committer.add_skip(iog)

    def finish(self):
        self.hog_committer.assert_finished()
        self.output_writer.flush()

    def close(self):
        self.output_writer.close()


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

    This is the parent-writer version of WriteDuplications().
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



def OrthologPipelineWriterProcess(
        writer_queue,
        writer_status_queue,
        output_writer,
        iogs_ordered,
        n_workers,
    ):

    pipeline_writer = PipelineOutputWriter(
        output_writer,
        iogs_ordered
    )

    active_workers = n_workers
    fatal = False

    try:
        while active_workers > 0:
            msg = writer_queue.get()

            if msg is None:
                active_workers -= 1
                continue

            if not isinstance(msg, tuple):
                raise TypeError(
                    "Unexpected writer message: %s %r" %
                    (type(msg), msg)
                )

            kind = msg[0]

            if kind == "result":
                _, iog, result = msg

                pipeline_writer.write_result(result)

                writer_status_queue.put((
                    "written",
                    iog,
                    result["n_orthologues"]
                ))

            elif kind == "skip":
                _, iog = msg

                pipeline_writer.write_skip(iog)

                writer_status_queue.put((
                    "written_skip",
                    iog
                ))

            elif kind == "error":
                writer_status_queue.put(msg)
                fatal = True
                break

            else:
                raise TypeError(
                    "Unexpected writer message kind: %r" % kind
                )

        if not fatal:
            pipeline_writer.finish()

            writer_status_queue.put((
                "writer_done",
                dict(output_writer.hog_writer.iHOG)
            ))

    except Exception:
        writer_status_queue.put((
            "writer_error",
            traceback.format_exc()
        ))

    finally:
        try:
            pipeline_writer.close()
        except Exception:
            writer_status_queue.put((
                "writer_close_error",
                traceback.format_exc()
            ))


def Worker_RunOrthologsMethod_Pipeline(
        tree_analyser,
        nspecies,
        args_queue,
        writer_queue,
        progress_queue,
        fewer_open_files,
        n_ologs_cache=100,
        write_hog_tree=False,
        fix_files=False
    ):
    """
    Analysis worker.

    Heavy work happens here:
        result = tree_analyser.AnalyseTree(iog)

    It sends full results to writer_queue.
    It sends only worker errors/sentinels to progress_queue.
    """
    while True:
        try:
            iog = args_queue.get(True, 0.1)

            if iog is None:
                break

            result = tree_analyser.AnalyseTree(iog)

            if result is None:
                writer_queue.put(("skip", iog))
            else:
                writer_queue.put(("result", iog, result))

        except queue.Empty:
            continue

        except Exception:
            tb = traceback.format_exc()
            writer_queue.put(("error", None, tb))
            progress_queue.put(("error", None, tb))
            break

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
        old_version=False,
        write_hog_tree=False,
        fix_files=False,
        fd_limit=None,
        GRACE_PERIOD=10.0,
        STALL_TIMEOUT=120.0,
        writer_queue_size=None,
    ):
    """
    Pipeline architecture:

        analysis workers -> bounded writer_queue -> dedicated writer process

    HOG rows are ordered inside the writer process.
    Non-HOG rows are written immediately.
    """
    if old_version:
        print("WARNING: old_version parallel writer disabled; using pipeline writer mode.")

    if fd_limit is not None:
        if sys.platform.startswith("linux") or sys.platform == "darwin":
            set_file_descriptor_limit(fd_limit)
        else:
            warnings.warn(
                "File descriptor limit adjustment is not supported on %s."
                % sys.platform
            )

    if writer_queue_size is None:
        writer_queue_size = max(4 * nProcesses, 32)

    writer_queue = mp.Queue(maxsize=writer_queue_size)
    progress_queue = mp.Queue(maxsize=max(4 * nProcesses, 32))
    writer_status_queue = mp.Queue()

    progressbar, task = util.get_progressbar(total_tasks)
    progressbar.start()

    for _ in range(nProcesses):
        args_queue.put(None)

    writer_proc = mp.Process(
        target=OrthologPipelineWriterProcess,
        args=(
            writer_queue,
            writer_status_queue,
            output_writer,
            iogs_ordered,
            nProcesses,
        )
    )

    runningProcesses = [
        mp.Process(
            target=Worker_RunOrthologsMethod_Pipeline,
            args=(
                tree_analyser,
                nspecies,
                args_queue,
                writer_queue,
                progress_queue,
                fewer_open_files,
                n_ologs_cache,
                write_hog_tree,
                fix_files
            )
        )
        for _ in range(nProcesses)
    ]

    writer_proc.start()

    for proc in runningProcesses:
        proc.start()

    nOrthologues_SpPair = util.nOrtho_sp(nspecies)

    completed_tasks = 0
    skipped_tasks = 0
    active_workers = nProcesses
    fatal = False
    writer_done = False
    writer_hog_counts = None
    last_progress_time = time.time()

    try:
        while completed_tasks < total_tasks or active_workers > 0 or not writer_done:
            # Worker progress messages.
            try:
                msg = progress_queue.get(timeout=0.1)
            except queue.Empty:
                msg = "__EMPTY__"

            if msg == "__EMPTY__":
                if completed_tasks < total_tasks and time.time() - last_progress_time > STALL_TIMEOUT:
                    print(
                        "ERROR: Stalled for %ss "
                        "(completed %d/%d, active_workers=%d)." %
                        (
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
                _, iog = msg
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

            # Writer status messages.
            while True:
                try:
                    wmsg = writer_status_queue.get_nowait()
                except queue.Empty:
                    break

                if isinstance(wmsg, tuple) and wmsg[0] == "written":
                    _, iog, nOrtho = wmsg
                    nOrthologues_SpPair += nOrtho
                    completed_tasks += 1
                    progressbar.update(task, advance=1)
                    last_progress_time = time.time()

                elif isinstance(wmsg, tuple) and wmsg[0] == "written_skip":
                    _, iog = wmsg
                    skipped_tasks += 1
                    completed_tasks += 1
                    progressbar.update(task, advance=1)
                    last_progress_time = time.time()

                elif isinstance(wmsg, tuple) and wmsg[0] == "writer_done":
                    _, writer_hog_counts = wmsg
                    writer_done = True

                elif isinstance(wmsg, tuple) and wmsg[0] in {
                    "writer_error",
                    "writer_close_error",
                    "error",
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

            if active_workers == 0 and not writer_proc.is_alive() and not writer_done:
                print("ERROR: writer process exited without writer_done.")
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

        writer_proc.join(timeout=GRACE_PERIOD)

        if writer_proc.is_alive():
            writer_proc.terminate()

        writer_proc.join()

        progressbar.stop()

        for q in (writer_queue, progress_queue, writer_status_queue):
            try:
                q.close()
                q.join_thread()
            except Exception:
                pass

    if writer_hog_counts is not None:
        # Update parent copy before TwoAndThreeGeneHOGs() runs.
        output_writer.hog_writer.iHOG.clear()
        output_writer.hog_writer.iHOG.update(writer_hog_counts)

    skip_rate = skipped_tasks / max(1, total_tasks)

    if skip_rate > 0.02:
        print(
            "WARNING: skipped %d/%d tasks (%.1f%%)." %
            (skipped_tasks, total_tasks, 100.0 * skip_rate)
        )

    if fatal:
        util.Fail()

    if not writer_done:
        print("ERROR: writer process did not finish cleanly.")
        util.Fail()

    return nOrthologues_SpPair


def set_file_descriptor_limit(fd_limit) -> None:
    """
    Try to raise the soft open-file limit.

    This is only a convenience. The new lazy writer should not depend on this.
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
