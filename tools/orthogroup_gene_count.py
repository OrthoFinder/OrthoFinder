#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert OrthoFinder orthogroup files into gene-count tables.

INPUT
  Input may be:
    - Orthogroups.tsv
    - Orthogroups.txt
    - N*.tsv hierarchical orthogroup tables
    - a directory containing any of the above

OUTPUT
  For tabular files:
      INPUT.GeneCount.tsv

  For Orthogroups.txt:
      Orthogroups.GeneCount.tsv

Examples
--------
orthogroup_gene_count Orthogroups.tsv

orthogroup_gene_count Orthogroups.txt

orthogroup_gene_count N0.tsv

orthogroup_gene_count Results_*/

orthogroup_gene_count Results_*/ -o GeneCountTables/
"""

import os
import csv
import argparse


def is_orthogroup_table(path):
    if not os.path.isfile(path):
        return False

    basename = os.path.basename(path)

    if basename in {"Orthogroups.tsv", "Orthogroups.txt"}:
        return True

    if basename.startswith("N") and basename.endswith(".tsv"):
        return True

    return False


def is_colon_orthogroups_file(path):
    return os.path.basename(path) == "Orthogroups.txt"


def default_output_file(input_file):
    root, ext = os.path.splitext(input_file)

    if os.path.basename(input_file) == "Orthogroups.txt":
        return root + ".txt.GeneCount.tsv"

    return root + ".GeneCount" + ext


def convert_tabular_gene_counts(input_file, output_file=None):
    if output_file is None:
        output_file = default_output_file(input_file)

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        header = next(reader)

        # HOG tables contain:
        # HOG, Level, Description, Species...
        #
        # Orthogroups.tsv contains:
        # Orthogroup, Species...
        n_col_skip = 3 if header[0] == "HOG" else 1

        writer.writerow(header)

        for line in reader:
            counts = [
                0 if cell == "" else len(cell.split(", "))
                for cell in line[n_col_skip:]
            ]

            writer.writerow(line[:n_col_skip] + counts)

    print("Orthogroup gene count table written to:")
    print(output_file)


def convert_colon_gene_counts(input_file, output_file=None):
    if output_file is None:
        output_file = default_output_file(input_file)

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Orthogroup", "GeneCount"])

        for line in infile:
            line = line.strip()

            if not line:
                continue

            if ":" not in line:
                continue

            orthogroup, genes = line.split(":", 1)
            genes = genes.strip().split()

            writer.writerow([orthogroup, len(genes)])

    print("Orthogroup gene count table written to:")
    print(output_file)


def convert_gene_counts(input_file, output_file=None):
    if is_colon_orthogroups_file(input_file):
        convert_colon_gene_counts(input_file, output_file)
    else:
        convert_tabular_gene_counts(input_file, output_file)


def get_input_files(input_path):
    if os.path.isfile(input_path):
        if is_orthogroup_table(input_path):
            return [input_path]

        print("ERROR: input file is not a supported orthogroup file:")
        print(input_path)
        raise SystemExit(1)

    if os.path.isdir(input_path):
        files = []

        for root, _, filenames in os.walk(input_path):
            for filename in sorted(filenames):
                path = os.path.join(root, filename)

                if is_orthogroup_table(path):
                    files.append(path)

        return sorted(files)

    print("ERROR: input path does not exist:")
    print(input_path)
    raise SystemExit(1)


def convert_directory(input_dir, output_dir=None):
    input_dir = os.path.abspath(input_dir)

    files_to_process = get_input_files(input_dir)

    if not files_to_process:
        print("ERROR: no orthogroup files found in directory:")
        print(input_dir)
        raise SystemExit(1)

    if output_dir is None:
        output_dir = os.path.join(input_dir, "GeneCountTables")
    else:
        output_dir = os.path.abspath(output_dir)

    os.makedirs(output_dir, exist_ok=True)

    for input_file in files_to_process:
        rel_path = os.path.relpath(input_file, input_dir)
        output_file = os.path.join(output_dir, rel_path)
        output_file = default_output_file(output_file)

        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        print("\nProcessing:")
        print(input_file)

        convert_gene_counts(input_file, output_file)


def main(args=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="Orthogroups.tsv, Orthogroups.txt, N*.tsv file, or directory",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Output filename for file input, or output directory for directory input "
            "[default: INPUT.GeneCount.tsv or INPUT/GeneCountTables]"
        ),
    )

    parsed = parser.parse_args(args)

    input_path = os.path.abspath(parsed.input)

    if os.path.isdir(input_path):
        convert_directory(input_path, parsed.output)
    else:
        convert_gene_counts(input_path, parsed.output)


if __name__ == "__main__":
    main()