#!/usr/bin/env python3

"""
Convert compact OrthoFinder orthologue results into pairwise species files.

OrthoFinder runs without ``-pof`` write one compressed file per species in the
``Orthologues`` directory. This utility recreates the pairwise layout produced
by ``-pof``::

    Orthologues/Orthologues_SpeciesA/SpeciesA__v__SpeciesB.tsv

INPUT may be either an OrthoFinder results directory or its ``Orthologues``
directory.

Examples
--------
split_ortholog_files Results_*/

split_ortholog_files Results_*/Orthologues/

split_ortholog_files Results_*/ -o PairwiseOrthologues/

split_ortholog_files Results_*/ --compress
"""

import argparse
import csv
import gzip
import os
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def find_orthologues_directory(input_dir):
    input_dir = os.path.abspath(input_dir)

    if not os.path.isdir(input_dir):
        raise ValueError(f"input directory does not exist: {input_dir}")

    nested = os.path.join(input_dir, "Orthologues")
    if os.path.isdir(nested):
        return nested

    return input_dir


def compact_orthologue_files(orthologues_dir):
    files_by_species = {}

    for filename in sorted(os.listdir(orthologues_dir)):
        if filename.endswith(".tsv.gz"):
            species = filename[:-7]
        elif filename.endswith(".tsv"):
            species = filename[:-4]
        else:
            continue

        path = os.path.join(orthologues_dir, filename)
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
            "expected at least two compact .tsv or .tsv.gz orthologue files in: "
            f"{orthologues_dir}"
        )

    return files_by_species


def pairwise_output_path(
    output_dir, source_species, target_species, compress=False
):
    suffix = ".tsv.gz" if compress else ".tsv"
    return os.path.join(
        output_dir,
        "Orthologues_" + source_species,
        f"{source_species}__v__{target_species}{suffix}",
    )


def split_species_file(task):
    source_species, source_file, species_names, targets, output_dir, compress = task
    source_output_dir = os.path.join(output_dir, "Orthologues_" + source_species)
    if targets:
        os.makedirs(source_output_dir, exist_ok=True)

    output_count = 0
    row_count = 0

    with ExitStack() as stack:
        writers = {}
        for target_species in targets:
            output_file = pairwise_output_path(
                output_dir, source_species, target_species, compress
            )
            outfile = stack.enter_context(open_text(output_file, "wt"))
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerow(("Orthogroup", source_species, target_species))
            writers[target_species] = writer
            output_count += 1

        with open_text(source_file, "rt") as infile:
            reader = csv.reader(infile, delimiter="\t")
            header = next(reader, None)
            expected_header = [
                "Orthogroup",
                "Species",
                source_species,
                "Orthologs",
            ]
            if header != expected_header:
                raise ValueError(
                    f"unexpected header in {source_file}: expected "
                    f"{expected_header!r}, found {header!r}"
                )

            for line_number, row in enumerate(reader, start=2):
                if len(row) != 4:
                    raise ValueError(
                        f"malformed row in {source_file} at line {line_number}: "
                        f"expected 4 columns, found {len(row)}"
                    )

                target_species = row[1]
                if target_species not in species_names or target_species == source_species:
                    raise ValueError(
                        f"unknown target species {target_species!r} in "
                        f"{source_file} at line {line_number}"
                    )

                writer = writers.get(target_species)
                if writer is not None:
                    writer.writerow((row[0], row[2], row[3]))
                    row_count += 1

    return output_count, row_count


def split_ortholog_files(
    input_dir,
    output_dir=None,
    compress=False,
    force=False,
    processes=None,
):
    orthologues_dir = find_orthologues_directory(input_dir)
    files_by_species = compact_orthologue_files(orthologues_dir)
    species_names = sorted(files_by_species)

    if output_dir is None:
        output_dir = orthologues_dir
    else:
        output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if processes is None:
        processes = min(4, len(species_names), os.cpu_count() or 1)
    if processes < 1:
        raise ValueError("processes must be at least 1")

    tasks = []
    for source_species in species_names:
        targets = []
        for target_species in species_names:
            if source_species == target_species:
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
        help="OrthoFinder results directory or its Orthologues directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="output directory [default: the input Orthologues directory]",
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
    parsed = parser.parse_args(args)

    try:
        output_count, row_count, output_dir = split_ortholog_files(
            parsed.input,
            output_dir=parsed.output,
            compress=parsed.compress,
            force=parsed.force,
            processes=parsed.processes,
        )
    except (OSError, ValueError) as error:
        parser.exit(1, f"ERROR: {error}\n")

    print("Pairwise orthologue files written to:")
    print(output_dir)
    print(f"Created {output_count} files containing {row_count} orthologue rows.")


if __name__ == "__main__":
    main()
