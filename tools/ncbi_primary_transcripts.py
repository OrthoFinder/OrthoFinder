#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:30:34 2024

@author: user
"""

import zipfile
import os
import sys
from Bio import SeqIO
from collections import defaultdict

def unzip_file(zip_path, extract_to):
    """Unzip the given ZIP file."""
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

def find_files(root_dir):
    """Find the protein FASTA and GFF files in the extracted directory."""
    protein_fasta = None
    gff_file = None
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith("protein.faa"):
                protein_fasta = os.path.join(root, file)
            elif file.endswith("genomic.gff"):
                gff_file = os.path.join(root, file)
    return protein_fasta, gff_file

def parse_fasta(fasta_file):
    """Parse the FASTA file and store the sequences by protein ID."""
    protein_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_dict[record.id] = str(record.seq)
    return protein_dict

def parse_gff(gff_file):
    """Parse the GFF file to map protein IDs to gene names."""
    gene_map = {}
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attributes = fields[8]
            if "protein_id=" in attributes and "gene=" in attributes:
                protein_id = attributes.split("protein_id=")[-1].split(";")[0]
                gene = attributes.split("gene=")[-1].split(";")[0]
                gene_map[protein_id] = gene
    return gene_map

def get_longest_transcripts(protein_dict, gene_map):
    """Find the longest protein sequence for each gene."""
    gene_transcripts = defaultdict(list)
    
    # Group protein sequences by gene
    for protein_id, sequence in protein_dict.items():
        if protein_id in gene_map:
            gene = gene_map[protein_id]
            gene_transcripts[gene].append((protein_id, sequence))

    # For each gene, select the longest transcript
    longest_transcripts = {}
    for gene, transcripts in gene_transcripts.items():
        longest_protein = max(transcripts, key=lambda x: len(x[1]))
        longest_transcripts[gene] = longest_protein

    return longest_transcripts

def write_longest_fasta(output_file, longest_transcripts):
    """Write the longest transcript for each gene to a new FASTA file."""
    with open(output_file, 'w') as out_fasta:
        for gene, (protein_id, sequence) in longest_transcripts.items():
            out_fasta.write(f">{protein_id} gene={gene}\n")
            out_fasta.write(f"{sequence}\n")

def process_zip_file(zip_file_path, output_dir):
    """Process a single ZIP file."""
    extract_to = zip_file_path.replace('.zip', '')
    unzip_file(zip_file_path, extract_to)

    # Step 2: Find the protein FASTA and GFF files
    protein_fasta, gff_file = find_files(extract_to)

    if not protein_fasta or not gff_file:
        print(f"Could not find required files in {zip_file_path}")
        return

    # Step 3: Parse the protein FASTA file
    protein_dict = parse_fasta(protein_fasta)

    # Step 4: Parse the GFF file to map protein IDs to genes
    gene_map = parse_gff(gff_file)

    # Step 5: Get the longest transcript for each gene
    longest_transcripts = get_longest_transcripts(protein_dict, gene_map)

    # Step 6: Write the longest transcripts to a new FASTA file
    output_fasta_path = os.path.join(output_dir, f"{os.path.basename(zip_file_path).replace('.zip', '_lt.fasta')}")
    write_longest_fasta(output_fasta_path, longest_transcripts)
    print(f"Longest transcripts saved to {output_fasta_path}")


if __name__ == "__main__":
    # Check if there are zip files in the current directory
    current_directory = os.getcwd()

    # Define the output directory: a subfolder 'primary_transcripts' in the current directory
    output_directory = os.path.join(current_directory, 'primary_transcripts')

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Process each ZIP file in the current directory
    for filename in os.listdir(current_directory):
        if filename.endswith(".zip"):
            zip_file_path = os.path.join(current_directory, filename)
            process_zip_file(zip_file_path, output_directory)