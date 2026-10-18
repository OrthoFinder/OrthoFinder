"""Writing the ortholog, HOG, duplication and suspect-gene output files."""

import glob
import gzip
import io
import os
import pickle
import tempfile
import time
import traceback

from ..utils import file_handler, file_io, util
from ..utils.util import csv_read_mode, csv_write_mode


# The text of the gzip files (GzipBuffers) held in memory: all files together, shared by
# the writer processes; and per file, before it is written as one gzip member.
GZIP_BUFFER_BYTES = 256 * 1024 * 1024
GZIP_MEMBER_BYTES = 1024 * 1024
# The same for the other output files (TextBuffers), appended to in pieces of this size
TEXT_BUFFER_BYTES = 256 * 1024 * 1024
TEXT_MEMBER_BYTES = 1024 * 1024


class GzipBuffers(object):
    """
    The text of gzip files, compressed and appended to each file as gzip
    members: a file's text once it has GZIP_MEMBER_BYTES, and when the files
    together hold more than budget bytes, the largest first, down to half of it. Memory is bounded whatever
    the number of files (no compressor is kept per file), and the members are
    large enough to compress almost as well as one (~2% larger at 64 KB, <1%
    at 1 MB). A file of several members is a standard gzip file.

    headers: the text a file starts with (its first member then replaces the
    file); a file with a header is written even with no other text (close).
    """

    compress = True

    def __init__(self, budget=GZIP_BUFFER_BYTES, member_bytes=GZIP_MEMBER_BYTES):
        self.budget = budget
        self.member_bytes = member_bytes
        self.headers = {}      # path -> header text
        self.buffers = {}      # path -> [text]
        self.sizes = {}        # path -> characters held
        self.total = 0
        self.written = set()   # files to append to (their first member is written)

    def start(self, path):
        """The file is written (with its header) even if no text is added."""
        if path not in self.buffers:
            self.buffers[path] = []
            self.sizes[path] = 0

    def write(self, path, text):
        self.start(path)
        self.buffers[path].append(text)
        self.sizes[path] += len(text)
        self.total += len(text)
        if self.sizes[path] >= self.member_bytes:
            self._write_member(path)
        if self.total > self.budget:
            # the largest first, down to half the budget: one sort per half budget
            # of text (not a search of all the files on every write)
            for largest in sorted(self.sizes, key=self.sizes.get, reverse=True):
                if self.total <= self.budget // 2:
                    break
                self._write_member(largest)

    def _write_member(self, path):
        text = "".join(self.buffers.pop(path))
        self.total -= self.sizes.pop(path)
        first = path not in self.written
        if first and path in self.headers:
            text = self.headers[path] + text
        elif not text:
            return
        mode = "wb" if first and path in self.headers else "ab"
        data = text.encode(file_io.ENCODING)
        with open(path, mode) as outfile:
            outfile.write(gzip.compress(data, compresslevel=file_io.GZIP_LEVEL, mtime=0) if self.compress else data)
        self.written.add(path)

    def close(self):
        for path in list(self.buffers):
            self._write_member(path)


class TextBuffers(GzipBuffers):
    """
    The text of (uncompressed) files, appended to each file in pieces, as
    GzipBuffers: memory is bounded and no file is kept open, so with many files
    (e.g. a HOG file per species tree node) a file is opened about once per
    TEXT_MEMBER_BYTES of its text, not once per write (a file cache with fewer
    open files than are written to closes and reopens them constantly: slow,
    above all on network file systems). The files exist (with their header) before.
    """

    compress = False

    def __init__(self, budget=TEXT_BUFFER_BYTES, member_bytes=TEXT_MEMBER_BYTES):
        super(TextBuffers, self).__init__(budget, member_bytes)


def sorted_lines(text):
    """
    The lines of text (each ending in a newline) sorted. The output files are
    written in orthogroup order with the lines of each orthogroup sorted: the
    same as sorting the whole file, without writing it again.
    """
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    if len(lines) < 2:
        return text if text.endswith("\n") else text + "\n"
    lines.sort()
    return "\n".join(lines) + "\n"


class LazyFileCache(object):
    """
    The output files of a writer process: their text is buffered and appended in
    large pieces (TextBuffers; gz files: GzipBuffers, paths without .gz here,
    gz_headers: the header of a file, which its first member replaces), so no
    file is kept open and memory is bounded, whatever the number of files.
    The files are complete once close_all() is called.
    """

    def __init__(self):
        self.text = TextBuffers()
        self.gz = GzipBuffers()
        self.gz_headers = self.gz.headers   # path + ".gz" -> header text

    def set_budget(self, share):
        """This process's share (0-1) of the memory for the buffered text."""
        self.text.budget = int(TEXT_BUFFER_BYTES * share)
        self.gz.budget = int(GZIP_BUFFER_BYTES * share)

    def start_gz(self, path):
        """path + ".gz" is written (with its header) even if no text is added."""
        self.gz.start(path + ".gz")

    def write(self, path, text, gz=False):
        """Append text to path (the file exists, with its header, or gz_headers has it)."""
        if not text:
            return
        if gz:
            self.gz.write(path + ".gz", text)
        else:
            self.text.write(path, text)

    def close_all(self):
        self.text.close()
        self.gz.close()


class BufferedFile(object):
    """A file of a LazyFileCache, for code that writes to a file object (fh.write)."""

    def __init__(self, cache, path):
        self.cache = cache
        self.path = path

    def write(self, text):
        self.cache.write(self.path, text)


class ParentOutputWriter(object):
    """
    Owns all shared output files.

    Workers must not write orthologs, xenologues, duplications, suspect genes,
    or HOG files. They return data. This parent writer writes the data.

    The results are written in orthogroup order, the lines of each orthogroup
    sorted, so each file is written once, already sorted: the orthogroups
    analysed from their trees (write_result, or _write_non_hog_payload from
    the parallel writers, which commit them in order), then write_last: the
    orthologs of the orthogroups of 2 or 3 genes (last_ortholog_lines: the
    last orthogroups, those with the fewest genes) and the misplaced genes,
    sorted by name (held until then).
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
            save_space,
            write_hog_tree,
            last_ortholog_lines=None,
        ):
        self.d_results_orthologs = d_results_orthologs
        self.species_dict = species_dict
        self.species_to_use = list(species_to_use)
        self.sequence_dict = sequence_dict
        self.spec_seq_dict = spec_seq_dict
        self.stride_dups = stride_dups
        self.hog_writer = hog_writer
        self.nspecies = len(self.species_to_use)

        self.save_space = save_space

        self.write_hog_tree = write_hog_tree
        self.need_ortholog_output = not write_hog_tree

        self.cache = LazyFileCache()

        self.last_ortholog_lines = last_ortholog_lines
        self.suspect_genes = {}   # i_species -> misplaced gene names, written by write_last
        self._species_index = None   # species ID (str) -> i_species (write_suspect_genes)

        self.putative_xenolog_dir = None
        self.suspect_genes_dir = None
        self.duplications_path = None
        if self.need_ortholog_output:
            self.putative_xenolog_dir = file_handler.handler.get_putative_xenelogs_dir()
            self.suspect_genes_dir = file_handler.handler.get_suspect_genes_dir()
            self.duplications_path = file_handler.handler.get_duplications_fn()
            self.initialise_ortholog_outputs()

    def initialise_ortholog_outputs(self):
        for i in range(self.nspecies):
            sp0 = str(self.species_to_use[i])
            sp0_name = self.species_dict[sp0]

            filename = os.path.join(
                self.d_results_orthologs,
                "%s.tsv" % sp0_name
            )
            header = io.StringIO()
            file_io.writer(header).writerow((
                "Orthogroup",
                "Species",
                sp0_name,
                "Orthologs"
            ))
            if self.save_space:
                # written with the orthologs (GzipBuffers), by the writer that
                # owns it: even with no orthologs, by write_last
                self.cache.gz_headers[filename + ".gz"] = header.getvalue()
            else:
                with open(filename, util.csv_write_mode, encoding="utf-8", newline="\n") as outfile:
                    outfile.write(header.getvalue())

        initialise_suspect_genes_dirs(
            self.nspecies,
            self.species_to_use,
            self.species_dict
        )

        with open(self.duplications_path, util.csv_write_mode, encoding="utf-8", newline="\n") as outfile:
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

    def ortholog_path(self, i):
        """The ortholog file of species i, Orthologs/<species>.tsv (-pof: split by finish.py)."""
        sp0_name = self.species_dict[str(self.species_to_use[i])]
        return os.path.join(self.d_results_orthologs, "%s.tsv" % sp0_name)

    def xenolog_path(self, i):
        sp0 = str(self.species_to_use[i])
        sp0_name = self.species_dict[sp0]
        return os.path.join(self.putative_xenolog_dir, "%s.tsv" % sp0_name)

    def write_ortholog_lines(self, ortholog_lines):
        for i in range(self.nspecies):
            text = ortholog_lines[i]
            if text:
                self.cache.write(self.ortholog_path(i), sorted_lines(text), gz=self.save_space)

    def write_xenolog_lines(self, ortholog_sus_lines):
        for i in range(self.nspecies):
            text = ortholog_sus_lines[i]
            if text:
                self.cache.write(self.xenolog_path(i), sorted_lines(text))

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

        text = "".join(sorted(file_io.unquoted_line(row) for row in rows))

        self.cache.write(self.duplications_path, text)

    def write_suspect_genes(self, suspect_genes):
        """The misplaced genes of one orthogroup: held, and written sorted by write_last."""
        if not suspect_genes:
            return

        if self._species_index is None:
            self._species_index = {str(sp): i for i, sp in enumerate(self.species_to_use)}

        # one pass over the genes (not one per species: with many species most have none)
        for g in suspect_genes:
            i = self._species_index.get(g.split("_", 1)[0])
            if i is not None:
                self.suspect_genes.setdefault(i, []).append(self.sequence_dict[g])

    def write_last(self, owns=None):
        """
        After the orthogroups analysed from their trees: the orthologs of the
        orthogroups of 2 or 3 genes and the misplaced genes, of the species
        whose files this writer owns (owns(i_species), default all).
        """
        if not self.need_ortholog_output:
            return
        if owns is None:
            owns = lambda i: True
        if self.save_space:
            for i in range(self.nspecies):
                if owns(i):
                    self.cache.start_gz(self.ortholog_path(i))
        if self.last_ortholog_lines is not None:
            self.write_ortholog_lines([
                text if owns(i) else ""
                for i, text in enumerate(self.last_ortholog_lines)
            ])
        for i, names in sorted(self.suspect_genes.items()):
            if not owns(i):
                continue
            path = os.path.join(
                self.suspect_genes_dir,
                self.species_dict[str(self.species_to_use[i])] + ".txt"
            )
            with open(path, util.csv_write_mode, encoding="utf-8", newline="\n") as outfile:
                outfile.write(sorted_lines("\n".join(names) + "\n"))
        self.suspect_genes = {}

    def write_hog_rows(self, cached_hogs):
        if cached_hogs:
            self.hog_writer.write_cached_hogs(cached_hogs, lock_hogs=None)


    def write_result(self, result):
        """One analysed OG's output (n_parallel == 1: the results come in OG order)."""
        if result is None:
            return
        self.write_hog_rows(result.get("cached_hogs", []))
        if self.need_ortholog_output:
            self.write_duplications(result["iog"], result.get("duplications", []))
            self.write_suspect_genes(result.get("suspect_genes", set()))
            self.write_ortholog_lines(result.get("olog_lines", []))
            self.write_xenolog_lines(result.get("olog_sus_lines", []))

    def close(self):
        """Write the rest of the buffered text: the files are complete."""
        try:
            self.cache.close_all()
        finally:
            self.hog_writer.close_files()


# Out-of-order HOG batches held in memory (pickled) before spilling to disk.
HOG_PENDING_MEMORY_BYTES = 1024 * 1024 * 1024
# The same for the other output (orthologs etc.), shared by its writer processes.
NON_HOG_PENDING_MEMORY_BYTES = 1024 * 1024 * 1024


class OrderedCommitter(object):
    """
    Commit the output of each OG in deterministic OG order (iogs_ordered),
    whatever order the OGs are analysed in: commit(data) for each OG's data,
    in order (not for an OG without any, add_skip).

    Out-of-order data are kept pickled in memory up to max_pending_bytes
    and spilled to a temporary file beyond that.
    """

    def __init__(
            self,
            commit,
            iogs_ordered,
            max_pending_bytes=None,
            spill_dir=None,
        ):
        self.commit = commit
        self.iogs_ordered = list(sorted(iogs_ordered))
        # iog -> ("mem", pickled bytes) or ("disk", offset, length)
        self.pending = {}
        self.next_index = 0
        if max_pending_bytes is None:
            max_pending_bytes = HOG_PENDING_MEMORY_BYTES
        self.max_pending_bytes = max_pending_bytes
        self.pending_bytes = 0
        self.spill_dir = spill_dir
        self.spill_file = None

    def add_result(self, iog, data):
        if (
            self.next_index < len(self.iogs_ordered)
            and iog == self.iogs_ordered[self.next_index]
        ):
            # The common case: nothing to buffer.
            self._commit(data)
            self.next_index += 1
        else:
            self._buffer(iog, data)
        self._drain_ready()

    def add_skip(self, iog):
        self.add_result(iog, None)

    def _buffer(self, iog, data):
        # Pickled bytes are far smaller than the live row objects. Past the
        # memory budget they go to a temporary file, so a slow early OG cannot
        # make the buffer grow without bound.
        data = pickle.dumps(data, pickle.HIGHEST_PROTOCOL)
        if self.pending_bytes + len(data) <= self.max_pending_bytes:
            self.pending[iog] = ("mem", data)
            self.pending_bytes += len(data)
            return
        if self.spill_file is None:
            self.spill_file = tempfile.TemporaryFile(
                prefix="orthofinder_spill_", dir=self.spill_dir
            )
        self.spill_file.seek(0, os.SEEK_END)
        offset = self.spill_file.tell()
        self.spill_file.write(data)
        self.pending[iog] = ("disk", offset, len(data))

    def _load(self, entry):
        if entry[0] == "mem":
            self.pending_bytes -= len(entry[1])
            return pickle.loads(entry[1])
        _, offset, length = entry
        self.spill_file.seek(offset)
        return pickle.loads(self.spill_file.read(length))

    def _commit(self, data):
        if data:
            self.commit(data)

    def _drain_ready(self):
        while self.next_index < len(self.iogs_ordered):
            next_iog = self.iogs_ordered[self.next_index]

            if next_iog not in self.pending:
                break

            self._commit(self._load(self.pending.pop(next_iog)))
            self.next_index += 1

    def close(self):
        if self.spill_file is not None:
            self.spill_file.close()
            self.spill_file = None

    def assert_finished(self):
        name = type(self).__name__
        if self.pending:
            raise RuntimeError(
                "%s finished with %d pending OG batches." % (name, len(self.pending))
            )

        if self.next_index != len(self.iogs_ordered):
            raise RuntimeError(
                "%s stopped early: committed %d/%d OGs."
                % (name, self.next_index, len(self.iogs_ordered))
            )


class OrderedHogCommitter(OrderedCommitter):
    """Commit HOG rows (cached_hogs) in deterministic OG order."""

    def __init__(self, hog_writer, iogs_ordered, max_pending_bytes=None, spill_dir=None):
        self.hog_writer = hog_writer
        super(OrderedHogCommitter, self).__init__(
            lambda cached_hogs: hog_writer.write_cached_hogs(cached_hogs, lock_hogs=None),
            iogs_ordered,
            max_pending_bytes=max_pending_bytes,
            spill_dir=spill_dir,
        )


def initialise_suspect_genes_dirs(nspecies, species_ids, species_dict):
    file_handler.handler.get_suspect_genes_dir()  # creates the directory
    d_suspect_orthologs = file_handler.handler.get_putative_xenelogs_dir()
    for index1 in range(nspecies):
        with open(d_suspect_orthologs + '%s.tsv' % species_dict[str(species_ids[index1])], util.csv_write_mode, encoding="utf-8", newline="\n") as outfile:
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

    for i, text in enumerate(ortholog_lines):
        if text:
            payload_for(i % n_writers)["olog_chunks"].append((i, text))

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
    """One OG's output (committed in OG order: see ParentOutputWriter), the lines of each file sorted."""
    for i, text in payload.get("olog_chunks", []):
        output_writer.cache.write(output_writer.ortholog_path(i), sorted_lines(text), gz=output_writer.save_space)

    for i, text in payload.get("xenolog_chunks", []):
        output_writer.cache.write(output_writer.xenolog_path(i), sorted_lines(text))

    for i, text in payload.get("suspect_chunks", []):
        # held, and written sorted by write_last
        output_writer.suspect_genes.setdefault(i, []).extend(text.splitlines())

    duplications = payload.get("duplications", [])
    if duplications:
        output_writer.write_duplications(
            payload["iog"],
            duplications,
        )


# A writer's progress messages to the parent: at most one per this many seconds (the
# parent uses them only to tell a slow writer from a stuck one)
STATUS_INTERVAL = 1.0


def non_hog_writer_process(
        writer_queue,
        writer_status_queue,
        output_writer,
        n_workers,
        writer_id,
        n_writers,
        iogs_ordered,
        spill_dir=None,
    ):
    """
    Write only files exclusively owned by this non-HOG writer (those of the
    species i with i % n_writers == writer_id, and writer 0: Duplications.tsv),
    in OG order: each worker sends every writer a result or a skip for each OG.
    """
    committer = OrderedCommitter(
        lambda payload: _write_non_hog_payload(output_writer, payload),
        iogs_ordered,
        max_pending_bytes=NON_HOG_PENDING_MEMORY_BYTES // max(1, n_writers),
        spill_dir=spill_dir,
    )
    # the text held in memory: the writers share the budget
    output_writer.cache.set_budget(1.0 / max(1, n_writers))
    active_workers = n_workers
    error_text = None
    last_status = 0.0

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
                committer.add_result(msg[1]["iog"], msg[1])
            elif kind == "skips":
                for iog in msg[1]:
                    committer.add_skip(iog)
            elif kind == "error":
                raise RuntimeError(msg[-1])
            else:
                raise TypeError(
                    "Unexpected non-HOG writer message kind: %r" % kind
                )

            now = time.monotonic()
            if now - last_status >= STATUS_INTERVAL:
                writer_status_queue.put(("non_hog_progress", writer_id))
                last_status = now

        committer.assert_finished()
        output_writer.write_last(lambda i: i % n_writers == writer_id)

    except Exception:
        error_text = traceback.format_exc()

    try:
        committer.close()
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
    last_status = 0.0
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

            now = time.monotonic()
            if now - last_status >= STATUS_INTERVAL:
                writer_status_queue.put(("hog_progress",))
                last_status = now

        committer.assert_finished()
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


# split_ortholog_files: the text held for the pairwise files of a species, all together
SPLIT_BUFFER_BYTES = 64 * 1024 * 1024


def split_ortholog_files(d_orthologs, q_compress=False, buffer_bytes=None, n_processes=1):
    """
    The pairwise ortholog files (-pof), Orthologs_<sp0>/<sp0>__v__<sp1>.tsv, from the
    file of each species, <sp0>.tsv[.gz]: each species' file read once, its rows
    held for their pairwise files and appended to them all together once
    buffer_bytes (SPLIT_BUFFER_BYTES) are held (_split_species_file): memory
    bounded whatever the number of species or the size of the file, and no file
    kept open (one input, and an output while it is written to).
    n_processes: the species' files split at the same time (each into its own
    directory), sharing buffer_bytes.
    """
    if not d_orthologs.endswith("/"):
        d_orthologs += "/"
    filenames = list(glob.glob(d_orthologs + "*.tsv") + glob.glob(d_orthologs + "*.tsv.gz"))
    species = [
        os.path.basename(fn)[:-len(".tsv.gz")] if fn.endswith(".gz") else os.path.basename(fn)[:-len(".tsv")]
        for fn in filenames
    ]
    if buffer_bytes is None:
        buffer_bytes = SPLIT_BUFFER_BYTES
    n_processes = max(1, min(int(n_processes or 1), len(filenames)))
    tasks = []
    for fn, sp0 in zip(filenames, species):
        d_out = d_orthologs + "Orthologs_" + sp0 + "/"
        if not os.path.exists(d_out):
            os.mkdir(d_out)
        tasks.append((fn, d_out, sp0, [sp1 for sp1 in species if sp1 != sp0], q_compress,
                      max(1, buffer_bytes // n_processes)))
    if n_processes == 1:
        for task in tasks:
            _split_species_file(*task)
        return
    import concurrent.futures
    import multiprocessing
    # (new processes, not forks: OrthoFinder may have threads running, and a fork of
    # a process with threads can deadlock in the child)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes,
                                                mp_context=multiprocessing.get_context("spawn")) as executor:
        for future in [executor.submit(_split_species_file, *task) for task in tasks]:
            future.result()      # (an error in one is raised here)


# split_ortholog_files reads each species' file in pieces of this many characters
SPLIT_READ_CHARS = 4 * 1024 * 1024


def _split_species_file(fn, d_out, sp0, others, q_compress, buffer_bytes):
    """
    The pairwise files of sp0 with the species others, from its file fn, in one pass
    (split_ortholog_files). Per row only what it needs, with millions of rows: its
    text held in a list per file (appended to all together, once buffer_bytes are
    held), the species column cut out of it, and the file read in large pieces.
    """
    out_fns, headers, pending = {}, {}, {}
    for sp1 in others:
        out_fns[sp1] = d_out + "%s__v__%s.tsv" % (sp0, sp1) + (".gz" if q_compress else "")
        headers[sp1] = file_io.unquoted_line(("Orthogroup", sp0, sp1))
        pending[sp1] = []
    written = set()   # the files created (with their header)

    def write_out():
        # each file's rows in order: the first with its header (one open of each file
        # of a small species file), the others appended; gzip: a member each time
        for sp1, rows in pending.items():
            first = sp1 not in written
            if not rows and not first:
                continue
            text = (headers[sp1] if first else "") + "".join(rows)
            if q_compress:
                with open(out_fns[sp1], "wb" if first else "ab") as outfile:
                    outfile.write(gzip.compress(text.encode(file_io.ENCODING), compresslevel=file_io.GZIP_LEVEL, mtime=0))
            else:
                with open(out_fns[sp1], "w" if first else "a", encoding="utf-8", newline="\n") as outfile:
                    outfile.write(text)
            written.add(sp1)
            rows.clear()

    held = 0
    with file_io.open_text(fn, csv_read_mode) as infile:
        rest = infile.readline() and ""   # (the header)
        while True:
            block = infile.read(SPLIT_READ_CHARS)
            if not block:
                lines = [rest] if rest else []
            else:
                lines = (rest + block).split("\n")
                rest = lines.pop()
            for line in lines:
                # OG, species, genes1, genes2 -> OG, genes1, genes2
                i1 = line.find("\t")
                i2 = line.find("\t", i1 + 1)
                if i1 < 0 or i2 < 0:
                    continue
                rows = pending.get(line[i1 + 1:i2])
                if rows is None or line.count("\t") != 3:
                    continue
                if "\r" in line:
                    # (as file_io.unquoted_line writes the fields: those of a "\r\n" line)
                    row = line.rstrip("\r\n").split("\t")
                    text = file_io.unquoted_line((row[0], row[2], row[3]))
                else:
                    text = line[:i1] + line[i2:] + "\n"
                rows.append(text)
                held += len(text)
            if held > buffer_bytes:
                write_out()
                held = 0
            if not block:
                break
    write_out()
