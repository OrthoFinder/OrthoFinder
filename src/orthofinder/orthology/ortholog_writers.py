"""Writing the ortholog, HOG, duplication and suspect-gene output files."""

import glob
import gzip
import os
import pickle
import tempfile
import time
import traceback
from collections import OrderedDict

from ..utils import file_handler, file_io, util
from ..utils import paths as utils_paths
from ..utils.util import csv_read_mode, csv_write_mode
from .tree_analyser import iter_ortholog_row


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

        fh = utils_paths.file_open(path, mode, gz=gz)
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

    Workers must not write orthologs, xenologues, duplications, suspect genes,
    or HOG files. They return data. This parent writer writes the data.
    """

    def __init__(
            self,
            d_results_orthologs,
            species_dict,
            species_to_use,
            sequence_dict,
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
        self.d_results_orthologs = d_results_orthologs
        self.species_dict = species_dict
        self.species_to_use = list(species_to_use)
        self.sequence_dict = sequence_dict
        self.spec_seq_dict = spec_seq_dict
        self.stride_dups = stride_dups
        self.hog_writer = hog_writer
        self.nspecies = len(self.species_to_use)

        self.fewer_open_files = fewer_open_files or save_space
        self.save_space = save_space

        self.write_hog_tree = write_hog_tree
        self.fix_files = fix_files
        self.need_ortholog_output = (not write_hog_tree or not fix_files)

        self.cache = LazyFileCache(max_open=max_open)

        self.flush_every_results = int(flush_every_results)
        self.flush_every_seconds = float(flush_every_seconds)
        self.results_since_flush = 0
        self.last_flush_time = time.time()

        self.putative_xenolog_dir = file_handler.handler.get_putative_xenelogs_dir()
        self.suspect_genes_dir = file_handler.handler.get_suspect_genes_dir()
        self.duplications_path = file_handler.handler.get_duplications_fn()

        if self.need_ortholog_output:
            self.initialise_ortholog_outputs()

    def initialise_ortholog_outputs(self):
        for i in range(self.nspecies):
            sp0 = str(self.species_to_use[i])
            sp0_name = self.species_dict[sp0]

            if self.fewer_open_files:
                filename = os.path.join(
                    self.d_results_orthologs,
                    "%s.tsv" % sp0_name
                )
                with utils_paths.file_open(filename, util.csv_write_mode, gz=self.save_space) as outfile:
                    writer = file_io.writer(outfile)
                    writer.writerow((
                        "Orthogroup",
                        "Species",
                        sp0_name,
                        "Orthologs"
                    ))
            else:
                d = os.path.join(
                    self.d_results_orthologs,
                    "Orthologs_" + sp0_name
                )
                if not os.path.exists(d):
                    os.mkdir(d)

                for j in range(self.nspecies):
                    if j == i:
                        continue

                    sp1 = str(self.species_to_use[j])
                    sp1_name = self.species_dict[sp1]

                    fn = os.path.join(
                        d,
                        "%s__v__%s.tsv" % (sp0_name, sp1_name)
                    )

                    with open(fn, util.csv_write_mode) as outfile:
                        writer = file_io.writer(outfile)
                        writer.writerow(("Orthogroup", sp0_name, sp1_name))

        initialise_suspect_genes_dirs(
            self.nspecies,
            self.species_to_use,
            self.species_dict
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
        sp0 = str(self.species_to_use[i])
        sp0_name = self.species_dict[sp0]

        if self.fewer_open_files:
            return os.path.join(
                self.d_results_orthologs,
                "%s.tsv" % sp0_name
            )

        sp1 = str(self.species_to_use[j])
        sp1_name = self.species_dict[sp1]

        d = os.path.join(
            self.d_results_orthologs,
            "Orthologs_" + sp0_name
        )

        return os.path.join(
            d,
            "%s__v__%s.tsv" % (sp0_name, sp1_name)
        )

    def xenolog_path(self, i):
        sp0 = str(self.species_to_use[i])
        sp0_name = self.species_dict[sp0]
        return os.path.join(self.putative_xenolog_dir, "%s.tsv" % sp0_name)

    def write_ortholog_lines(self, ortholog_lines):
        if self.fewer_open_files:
            for i in range(self.nspecies):
                text = ortholog_lines[i][0]
                if text:
                    self.cache.write(
                        self.ortholog_path(i),
                        text,
                        mode=util.csv_append_mode,
                        gz=self.save_space
                    )
        else:
            for i in range(self.nspecies):
                for j, text in iter_ortholog_row(ortholog_lines[i]):
                    if i == j:
                        continue
                    self.cache.write(
                        self.ortholog_path(i, j),
                        text,
                        mode=util.csv_append_mode,
                        gz=False
                    )

    def write_xenolog_lines(self, ortholog_sus_lines):
        for i in range(self.nspecies):
            text = ortholog_sus_lines[i]
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
        rows = duplication_rows(
            og_name,
            duplications,
            self.species_dict,
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

        species = list(map(str, self.species_to_use))

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
                self.species_dict[strsp0] + ".txt"
            )

            text = "\n".join(
                [self.sequence_dict[g] for g in these_genes]
            ) + "\n"

            self.cache.write(
                path,
                text,
                mode=util.csv_append_mode,
                gz=False
            )

    def write_hog_rows(self, cached_hogs):
        if cached_hogs:
            self.hog_writer.write_cached_hogs(cached_hogs, lock_hogs=None)


    def write_non_hog_result(self, result, do_flush=True):
        """
        Write non-HOG output for one analysed OG.

        These rows do not assign HOG IDs, so they can be streamed immediately.
        """
        if result is None:
            return

        if self.need_ortholog_output:
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
            self.hog_writer.write_cached_hogs(cached_hogs, lock_hogs=None)

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
    tree-based ortholog-writing step has already created the files.
    """

    def __init__(
            self,
            d_results_orthologs,
            species_dict,
            species_to_use,
            save_space=False,
            fewer_open_files=False,
            max_open=64,
        ):
        self.d_results_orthologs = d_results_orthologs
        self.species_dict = species_dict
        self.species_to_use = list(species_to_use)
        self.nspecies = len(self.species_to_use)

        self.save_space = save_space
        self.fewer_open_files = fewer_open_files or save_space

        self.cache = LazyFileCache(max_open=max_open)

    def ortholog_path(self, i, j=None):
        sp0 = str(self.species_to_use[i])
        sp0_name = self.species_dict[sp0]

        if self.fewer_open_files:
            return os.path.join(
                self.d_results_orthologs,
                "%s.tsv" % sp0_name
            )

        sp1 = str(self.species_to_use[j])
        sp1_name = self.species_dict[sp1]

        d = os.path.join(
            self.d_results_orthologs,
            "Orthologs_" + sp0_name
        )

        return os.path.join(
            d,
            "%s__v__%s.tsv" % (sp0_name, sp1_name)
        )

    def write_ortholog_lines(self, ortholog_lines):
        if self.fewer_open_files:
            for i in range(self.nspecies):
                text = ortholog_lines[i][0]
                if text:
                    self.cache.write(
                        self.ortholog_path(i),
                        text,
                        mode=util.csv_append_mode,
                        gz=self.save_space
                    )
        else:
            for i in range(self.nspecies):
                for j, text in iter_ortholog_row(ortholog_lines[i]):
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


def initialise_suspect_genes_dirs(nspecies, species_ids, species_dict):
    file_handler.handler.get_suspect_genes_dir()  # creates the directory
    d_suspect_orthologs = file_handler.handler.get_putative_xenelogs_dir()
    for index1 in range(nspecies):
        with open(d_suspect_orthologs + '%s.tsv' % species_dict[str(species_ids[index1])], util.csv_write_mode) as outfile:
            writer1 = file_io.writer(outfile)
            writer1.writerow(("Orthogroup", species_dict[str(species_ids[index1])], "Other"))


def duplication_rows(og_name, duplications, sp_ids, seq_ids, stride_dups):
    """
    Convert duplication records into rows.

    This converts duplication records for the parent writer.
    """
    rows = []

    for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
        q_terminal = not sp_node_id.startswith("N")

        if stride_dups is None:
            is_stride = "Terminal" if q_terminal else "Non-Terminal"
        else:
            if q_terminal:
                is_stride = "Terminal"
            elif frozenset(genes0 + genes1) in stride_dups:
                is_stride = "Non-Terminal: STRIDE"
            else:
                is_stride = "Non-Terminal"

        gene_list0 = ", ".join([seq_ids[g] for g in genes0])
        gene_list1 = ", ".join([seq_ids[g] for g in genes1])

        rows.append([
            og_name,
            sp_ids[sp_node_id] if q_terminal else sp_node_id,
            gene_node_name,
            str(frac),   # as text: all-str rows take the fast path
            is_stride,
            gene_list0,
            gene_list1
        ])

    return rows


def _empty_non_hog_payload(iog):
    return {
        "iog": iog,
        "olog_chunks": [],
        "xenolog_chunks": [],
        "suspect_chunks": [],
        "duplications": [],
    }


def partition_non_hog_result(
        result,
        n_writers,
        tree_analyser,
        fewer_open_files,
        need_ortholog_output,
    ):
    """
    Split one analysed OG into file-owner payloads for non-HOG output.

    Ortholog, xenologue, and suspect-gene files are owned by source species.
    Duplications.tsv is owned by writer 0. Empty writer payloads are omitted.
    """
    if not need_ortholog_output:
        return {}

    iog = result["iog"]
    payloads = {}

    def payload_for(owner):
        payload = payloads.get(owner)
        if payload is None:
            payload = _empty_non_hog_payload(iog)
            payloads[owner] = payload
        return payload

    ortholog_lines = result.get("olog_lines", [])

    if fewer_open_files:
        for i, row in enumerate(ortholog_lines):
            if not row:
                continue
            text = row[0]
            if not text:
                continue
            owner = i % n_writers
            payload_for(owner)["olog_chunks"].append((i, None, text))
    else:
        for i, row in enumerate(ortholog_lines):
            owner = i % n_writers
            for j, text in iter_ortholog_row(row):
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
            for i, sp in enumerate(tree_analyser.species_to_use)
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
                tree_analyser.sequence_dict[g]
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
        sp_id = str(output_writer.species_to_use[i])
        path = os.path.join(
            output_writer.suspect_genes_dir,
            output_writer.species_dict[sp_id] + ".txt"
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


def non_hog_writer_process(
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


def ordered_hog_writer_process(
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
        hog_counts = dict(hog_writer.i_hog)

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
        d_out = d_orthologs + "Orthologs_" + sp0 + "/"
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
