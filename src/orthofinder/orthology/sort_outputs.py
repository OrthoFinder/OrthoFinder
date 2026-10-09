"""Sorting and checking the ortholog and HOG output files."""

import gzip
import heapq
import multiprocessing as mp
import os

from ..utils import file_handler, parallel_task_manager, util


def sort_non_hog_output_files(
        n_parallel,
        species_to_use,
        species_dict,
        fewer_open_files,
        save_space,
        write_hog_tree,
        fix_files,
    ):

    if write_hog_tree and fix_files:
        return

    species = [species_dict[str(sp)] for sp in species_to_use]
    d_results_orthologs = file_handler.handler.get_orthologs_directory()

    fns = []

    # Orthologue files.
    if fewer_open_files or save_space:
        for sp in species:
            fns.append((
                os.path.join(d_results_orthologs, "%s.tsv" % sp),
                bool(save_space)
            ))
    else:
        for sp1 in species:
            d = os.path.join(d_results_orthologs, "Orthologues_" + sp1)
            for sp2 in species:
                if sp1 == sp2:
                    continue
                fns.append((
                    os.path.join(d, "%s__v__%s.tsv" % (sp1, sp2)),
                    False
                ))

    # Xenolog files.
    d_xenologs = file_handler.handler.get_putative_xenelogs_dir()
    for sp in species:
        fns.append((
            os.path.join(d_xenologs, "%s.tsv" % sp),
            False
        ))

    # Duplications file.
    fns.append((
        file_handler.handler.get_duplications_fn(),
        False
    ))

    # Compressed outputs are written as fn + ".gz" (see paths.file_open).
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

    parallel_task_manager.run_method_parallel(
        sort_file_by_first_column_no_repair,
        args_queue,
        n_parallel,
        total_tasks=n_sort_tasks,
        show_progress=False,
    )

    suspect_queue = mp.Queue()
    d_suspect_genes = file_handler.handler.get_suspect_genes_dir()
    suspect_fns = [
        os.path.join(d_suspect_genes, "%s.txt" % sp) for sp in species
    ]
    suspect_fns = sorted(
        (fn for fn in suspect_fns if os.path.exists(fn)),
        key=lambda fn: -os.path.getsize(fn),
    )
    n_suspect_tasks = 0
    for fn in suspect_fns:
        suspect_queue.put((fn,))
        n_suspect_tasks += 1

    parallel_task_manager.run_method_parallel(
        sort_plain_text_file,
        suspect_queue,
        n_parallel,
        total_tasks=n_suspect_tasks,
        show_progress=False,
    )


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


# Above this many characters a file is sorted in chunks that are spilled to
# disk and merged, so memory use per sorting process stays bounded.
SORT_CHUNK_CHARS = 128 * 1024 * 1024


def sort_lines_in_file(fn, gz=False, has_header=False, chunk_chars=SORT_CHUNK_CHARS):
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


def sort_file_by_first_column_no_repair(fn, gz=False):
    """
    Sort a TSV file by first column.

    This is for orthologs, xenologues, and duplications only.
    It must never be used for HOG files.

    Plain line order is the same as ordering by (first column, line): the
    first column is followed by a tab, which sorts before any character that
    can appear in an OG ID. Sorting whole lines avoids building a key per line.
    """
    sort_lines_in_file(fn, gz=gz, has_header=True)


def sort_plain_text_file(fn):
    """Sort a plain-text output file deterministically."""
    sort_lines_in_file(fn, gz=False, has_header=False)


def validate_hog_writer_no_duplicate_ids(hog_writer):

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
