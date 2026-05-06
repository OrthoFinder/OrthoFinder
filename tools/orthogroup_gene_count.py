#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert an OrthoFinder orthogroup table into a gene-count table.

INPUT
  Input may be:
    - Orthogroups.tsv
    - N*.tsv hierarchical orthogroup tables

The script replaces comma-separated gene lists with the number
of genes per species.

Examples
--------
orthogroup_gene_count Orthogroups.tsv

orthogroup_gene_count N0.tsv

orthogroup_gene_count Orthogroups.tsv -o Orthogroups.GeneCount.tsv
"""

import os
import csv
import argparse


def convert_gene_counts(input_file, output_file=None):

    if output_file is None:
        root, ext = os.path.splitext(input_file)
        output_file = root + ".GeneCount" + ext

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


def main(args=None):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="Orthogroups.tsv or hierarchical orthogroup .tsv file",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output filename [default: INPUT.GeneCount.tsv]",
    )

    parsed = parser.parse_args(args)

    convert_gene_counts(parsed.input, parsed.output)


if __name__ == "__main__":
    main()