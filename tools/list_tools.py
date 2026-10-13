# tools/list_tools.py

#!/usr/bin/env python3

def main():
    print(
"""
OrthoFinder utility scripts
===========================

primary_transcript
    Extract longest transcript isoform per gene.

ncbi_primary_transcript
    Extract longest transcripts from NCBI genome ZIP downloads.

make_ultrametric
    Convert rooted trees to ultrametric trees.

convert_orthofinder_tree_ids
    Convert between OrthoFinder IDs and gene accessions.

create_hog_fastas
    Create FASTA files for hierarchical orthogroups (HOGs).

orthogroup_gene_count
    Convert orthogroup tables into gene-count tables.

split_ortholog_files
    Convert compact ortholog results into pairwise species files.

view_results
    Write the pages to view the results of a finished run in a web browser
    (as a run with --view: gene trees, statistics, orthogroups and HOGs).

species_tree_alignment
    Write the concatenated alignment of the orthogroups chosen for a species
    tree, from the results of a finished run (e.g. one with -ST astral-pro).

mcl_orthogroup_sequences
    Write the sequences of the MCL orthogroups (before the HOGs, e.g. of a
    run with -og), a FASTA file for each, from the working directory.

Run any command with --help for detailed usage. From the OrthoFinder
directory, without installing OrthoFinder, each can also be run as a
script: python tools/<command>.py (e.g. python tools/make_ultrametric.py --help),
and view_results as python tools/view_results/view_results.py. The
tools/view_results directory also works on its own, e.g. copied or
downloaded next to the results: it needs Python and ete4, not OrthoFinder.

Examples
--------
primary_transcript --help
make_ultrametric --help
split_ortholog_files --help
view_results --help
species_tree_alignment --help
mcl_orthogroup_sequences --help
"""
    )

if __name__ == "__main__":
    main()
