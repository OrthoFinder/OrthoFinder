#!/usr/bin/env python3

"""
Convert compact OrthoFinder ortholog results into pairwise species files.

OrthoFinder runs without ``-pof`` write one compressed file per species in the
``Orthologs`` directory. This utility recreates the pairwise layout produced
by ``-pof``::

    Orthologs/Orthologs_SpeciesA/SpeciesA__v__SpeciesB.tsv

INPUT may be either an OrthoFinder results directory or its ``Orthologs``
directory.

Examples
--------
split_ortholog_files Results_*/

split_ortholog_files Results_*/Orthologs/

split_ortholog_files Results_*/ -o PairwiseOrthologs/

split_ortholog_files Results_*/ --compress

split_ortholog_files Results_*/ --species Homo_sapiens Mus_musculus

split_ortholog_files Results_*/ --species Homo_sapiens --partners Mus_musculus Danio_rerio

With many species there is a pairwise file for each pair (n x (n-1) files): --species
and --partners make only those needed.
"""

import argparse
import gzip
import os
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack


def open_text(path, mode):
    """A text file (gzipped if .gz) as OrthoFinder writes them: UTF-8, "\n" line endings."""
    newline = None if "r" in mode else "\n"
    if path.endswith(".gz"):
        # (level 6: ~8x faster than gzip's default 9, files <1% larger)
        level = {} if "r" in mode else {"compresslevel": 6}
        return gzip.open(path, mode, encoding="utf-8", newline=newline, **level)
    return open(path, mode, encoding="utf-8", newline=newline)


class TabWriter(object):
    """The rows of a tab-separated file, as OrthoFinder writes them: tab-joined, unquoted."""

    def __init__(self, outfile):
        self.outfile = outfile

    def writerow(self, row):
        self.outfile.write("\t".join(row) + "\n")


def find_orthologs_directory(input_dir):
    input_dir = os.path.abspath(input_dir)

    if not os.path.isdir(input_dir):
        raise ValueError(f"input directory does not exist: {input_dir}")

    # "Orthologues" in the results of OrthoFinder versions before 3.2
    for name in ("Orthologs", "Orthologues"):
        nested = os.path.join(input_dir, name)
        if os.path.isdir(nested):
            return nested

    return input_dir


def compact_ortholog_files(orthologs_dir):
    files_by_species = {}

    for filename in sorted(os.listdir(orthologs_dir)):
        if filename.endswith(".tsv.gz"):
            species = filename[:-7]
        elif filename.endswith(".tsv"):
            species = filename[:-4]
        else:
            continue

        path = os.path.join(orthologs_dir, filename)
        if not os.path.isfile(path):
            continue

        existing = files_by_species.get(species)
        if existing:
            # Prefer the directly readable file when both representations are
            # present. This also makes reruns tolerant of a previously
            # decompressed compact input.
            if existing.endswith(".tsv") or path.endswith(".gz"):
                continue
        files_by_species[species] = path

    if len(files_by_species) < 2:
        raise ValueError(
            "expected at least two compact .tsv or .tsv.gz ortholog files in: "
            f"{orthologs_dir}"
        )

    return files_by_species


def pairwise_output_path(
    output_dir, source_species, target_species, compress=False
):
    suffix = ".tsv.gz" if compress else ".tsv"
    return os.path.join(
        output_dir,
        "Orthologs_" + source_species,
        f"{source_species}__v__{target_species}{suffix}",
    )


# A source file of at most this size (as text) is read into memory, its rows grouped by
# target species, and the pairwise files are then written one at a time: one pass, one
# file open. A larger one is read once per batch of pairwise files open at once (at most
# MAX_OPEN_PER_PROCESS: the open-file limit, often 1024, is shared by the processes).
IN_MEMORY_BYTES_PER_PROCESS = 128 * 1024 * 1024
MAX_OPEN_PER_PROCESS = 128


def text_size(path):
    """The size of a file as text (a .gz file: about 6x its size)."""
    size = os.path.getsize(path)
    return size * 6 if path.endswith(".gz") else size


def max_open_per_process(n_processes):
    """MAX_OPEN_PER_PROCESS, or fewer: the processes share half the open-file limit."""
    try:
        import resource
        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    except (ImportError, ValueError, OSError):
        return MAX_OPEN_PER_PROCESS
    if soft == resource.RLIM_INFINITY:
        return MAX_OPEN_PER_PROCESS
    return max(8, min(MAX_OPEN_PER_PROCESS, soft // (2 * max(1, n_processes))))


def split_species_file(task):
    source_species, source_file, species_names, targets, output_dir, compress, n_processes = task
    species_names = set(species_names)
    source_output_dir = os.path.join(output_dir, "Orthologs_" + source_species)
    if targets:
        os.makedirs(source_output_dir, exist_ok=True)

    output_count = 0
    row_count = 0
    expected_header = [
        "Orthogroup",
        "Species",
        source_species,
        "Orthologs",
    ]

    # (one pass over the source file, checking it, even with no targets)
    q_in_memory = text_size(source_file) <= IN_MEMORY_BYTES_PER_PROCESS
    max_open = len(targets) if q_in_memory else max_open_per_process(n_processes)
    batches = [
        targets[i:i + max_open]
        for i in range(0, len(targets), max(1, max_open))
    ] or [[]]
    for i_batch, batch in enumerate(batches):
        rows_of = {}    # q_in_memory: target species -> its rows
        with ExitStack() as stack:
            writers = {}
            for target_species in batch:
                if q_in_memory:
                    rows_of[target_species] = []
                    writers[target_species] = BufferedRows(rows_of[target_species])
                    continue
                output_file = pairwise_output_path(
                    output_dir, source_species, target_species, compress
                )
                outfile = stack.enter_context(open_text(output_file, "wt"))
                writer = TabWriter(outfile)
                writer.writerow(("Orthogroup", source_species, target_species))
                writers[target_species] = writer
                output_count += 1

            infile = stack.enter_context(open_text(source_file, "rt"))
            reader = (line.rstrip("\r\n").split("\t") for line in infile)
            header = next(reader, None)
            if header != expected_header:
                raise ValueError(
                    f"unexpected header in {source_file}: expected "
                    f"{expected_header!r}, found {header!r}"
                )

            for line_number, row in enumerate(reader, start=2):
                if i_batch == 0:
                    # (the file is checked in the first pass)
                    if len(row) != 4:
                        raise ValueError(
                            f"malformed row in {source_file} at line {line_number}: "
                            f"expected 4 columns, found {len(row)}"
                        )
                    if row[1] not in species_names or row[1] == source_species:
                        raise ValueError(
                            f"unknown target species {row[1]!r} in "
                            f"{source_file} at line {line_number}"
                        )

                writer = writers.get(row[1])
                if writer is not None:
                    writer.writerow((row[0], row[2], row[3]))
                    row_count += 1

        for target_species, rows in rows_of.items():
            output_file = pairwise_output_path(
                output_dir, source_species, target_species, compress
            )
            with open_text(output_file, "wt") as outfile:
                writer = TabWriter(outfile)
                writer.writerow(("Orthogroup", source_species, target_species))
                for row in rows:
                    writer.writerow(row)
            output_count += 1
        rows_of.clear()

    return output_count, row_count


class BufferedRows(object):
    """A writer that keeps the rows (in a list), to write them later."""

    def __init__(self, rows):
        self.writerow = rows.append


def selected_species(names, species_names, option):
    """The species of names (default all: species_names), in the order of species_names."""
    if not names:
        return list(species_names)
    unknown = sorted(set(names) - set(species_names))
    if unknown:
        raise ValueError(
            f"{option}: no ortholog file for {', '.join(unknown)}. The species are: "
            f"{', '.join(species_names[:20])}"
            + (f", ... ({len(species_names)} species)" if len(species_names) > 20 else "")
        )
    names = set(names)
    return [name for name in species_names if name in names]


def split_ortholog_files(
    input_dir,
    output_dir=None,
    compress=False,
    force=False,
    processes=None,
    species=None,
    partners=None,
):
    """
    species: the species whose pairwise files are made (default all); partners:
    the other species of those files (default all).
    """
    orthologs_dir = find_orthologs_directory(input_dir)
    files_by_species = compact_ortholog_files(orthologs_dir)
    species_names = sorted(files_by_species)
    sources = selected_species(species, species_names, "--species")
    targets_wanted = set(selected_species(partners, species_names, "--partners"))

    if output_dir is None:
        output_dir = orthologs_dir
    else:
        output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if processes is None:
        processes = min(4, len(sources), os.cpu_count() or 1)
    if processes < 1:
        raise ValueError("processes must be at least 1")

    tasks = []
    for source_species in sources:
        targets = []
        for target_species in species_names:
            if source_species == target_species or target_species not in targets_wanted:
                continue
            output_file = pairwise_output_path(
                output_dir, source_species, target_species, compress
            )
            if force or not os.path.exists(output_file):
                targets.append(target_species)
        tasks.append(
            (
                source_species,
                files_by_species[source_species],
                species_names,
                targets,
                output_dir,
                compress,
                processes,
            )
        )

    if processes == 1:
        results = map(split_species_file, tasks)
    else:
        with ProcessPoolExecutor(max_workers=processes) as executor:
            results = executor.map(split_species_file, tasks)
            results = list(results)

    output_count = 0
    row_count = 0
    for species_output_count, species_row_count in results:
        output_count += species_output_count
        row_count += species_row_count

    return output_count, row_count, output_dir


def main(args=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="OrthoFinder results directory or its Orthologs directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="output directory [default: the input Orthologs directory]",
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        help="write pairwise files as .tsv.gz",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="replace existing pairwise output files [default: skip them]",
    )
    parser.add_argument(
        "-t",
        "--processes",
        type=int,
        default=None,
        help="worker processes [default: up to 4, limited by CPU and species count]",
    )
    parser.add_argument(
        "-s",
        "--species",
        nargs="+",
        default=None,
        metavar="NAME",
        help="make the pairwise files of only these species (names as in Orthologs/<name>.tsv[.gz]) "
             "[default: all]",
    )
    parser.add_argument(
        "-p",
        "--partners",
        nargs="+",
        default=None,
        metavar="NAME",
        help="with only these other species [default: all]",
    )
    parsed = parser.parse_args(args)

    try:
        output_count, row_count, output_dir = split_ortholog_files(
            parsed.input,
            output_dir=parsed.output,
            compress=parsed.compress,
            force=parsed.force,
            processes=parsed.processes,
            species=parsed.species,
            partners=parsed.partners,
        )
    except (OSError, ValueError) as error:
        parser.exit(1, f"ERROR: {error}\n")

    print("Pairwise ortholog files written to:")
    print(output_dir)
    print(f"Created {output_count} files containing {row_count} ortholog rows.")


if __name__ == "__main__":
    main()
