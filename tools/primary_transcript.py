
"""
Create primary transcript FASTA files by retaining the longest isoform per gene.

Supports:
- Single FASTA files
- Directories (searched recursively)
- .zip archives
- .tar.gz / .tgz / .tar archives

Gene identification
-------------------
By default the script automatically detects:

1. Ensembl-style FASTA headers
   using:
       gene:
       gene=
       locus:
       locus=

2. NCBI protein FASTA headers
   using isoform annotations.

Optional transcript parsing modes
---------------------------------
These modes define how transcript IDs are grouped into genes.

last_dot
    Remove the final dot suffix.

    Example:
        AT1G01010.1 -> AT1G01010
        AT1G01010.2 -> AT1G01010

space
    Use only the first whitespace-separated token.

    Example:
        transcript1 extra_info -> transcript1

last_dot_before_first_space
    Apply both rules:
    - take first token before whitespace
    - remove final dot suffix

    Example:
        AT1G01010.1 description text
        -> AT1G01010

Examples
--------
Process a FASTA file:
    primary_transcript proteins.fa

Process a directory recursively:
    primary_transcript proteomes/

Process an archive:
    primary_transcript proteomes.tar.gz

Use custom transcript grouping:
    primary_transcript proteins.fa last_dot
"""

FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fas", ".fsa", ".pep"}


import os
import re
import sys
import argparse
import tarfile
import zipfile
import tempfile
from collections import Counter, defaultdict

# Use the 'all' version rather than ab initio

def is_archive(fn):
    return zipfile.is_zipfile(fn) or tarfile.is_tarfile(fn)


def is_fasta_file(fn):
    return os.path.splitext(fn)[1].lower() in FASTA_EXTENSIONS

def process_fasta_file(fn, dout, gene_name_function_name=None):
    if not CheckFile(fn):
        return

    if gene_name_function_name is not None:
        q_use_original_accession_line = True
        ScanTags_with_fn(fn, gene_name_function_name)
    else:
        if IsNCBI(fn):
            print("Identified as NCBI file")
            gene_name_function_name = GetGeneName_NCBI
            q_use_original_accession_line = True
        else:
            gene_name_function_name = GetGeneName_Ensembl
            q_use_original_accession_line = False
            print('Looking for "gene=" or "gene:" to identify isoforms of same gene')

    CreatePrimaryTranscriptsFile(
        fn,
        dout,
        gene_name_function_name,
        q_use_original_accession_line,
    )

def process_archive(fn, dout, gene_name_function=None):

    print("\nExtracting archive: %s" % fn)

    with tempfile.TemporaryDirectory() as tmpdir:

        if zipfile.is_zipfile(fn):

            with zipfile.ZipFile(fn, "r") as zf:
                zf.extractall(tmpdir)

        else:

            with tarfile.open(fn, "r:*") as tf:
                tf.extractall(tmpdir)

        process_input(tmpdir, dout, gene_name_function)

def CheckFile(fn):
    """
    Checks for:
    - Duplicated accession lines
    """
    accs = set()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                a = l.rstrip()[1:]
                if a in accs:
                    print("\nERROR: duplicated sequence accession:\n%s" % a)
                    print("\nPlease correct this and then rerun the script.\n")
                    return False
                accs.add(a)
    return True

def ScanTags(fn):
    """
    For ensembl genomes, look for tag:id and count repeated ids
    :param fn:
    :return:
    """
    tags = set()
    tokens = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            tokens.append([t.split(":", 1) for t in line.rstrip().split() if ":" in t])
            tags.update([t[0] for t in tokens[-1]])
    for this_tag in tags:
        print(this_tag)
        # print(tokens[-1])
        c = Counter([idd for acc in tokens for t, idd in acc if t == this_tag])
        print(c.most_common(5))
        print("")

def ScanTags_NCBI(fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(line[1:].split(".", 1)[0])
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))

def ScanTags_with_fn(fn, gene_name_fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(gene_name_fn(line))
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))
    # print(genes[0])
    # print(sorted(genes)[:10])

def GetGeneName_Ensembl(acc_line):
    tokens = [(t.split("=") if "=" in t else t.split(":"))[1] for t in acc_line.rstrip().split() if ("gene:" in t or "gene=" in t or "locus:" in t or "locus=" in t)]
    if len(tokens) != 1: return None
    return tokens[0]

def IsNCBI(fn):
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                l = l.rstrip()
                if l.startswith(">NP_") and l.endswith("]"): return True
                elif l.startswith(">XP_") and l.endswith("]"): return True
                elif l.startswith(">YP_") and l.endswith("]"): return True
                elif l.startswith(">WP_") and l.endswith("]"): return True
                return False
    return False

def GetGeneName_NCBI(acc_line):
    acc_line = acc_line[1:]
    original = acc_line
    # look for "isoform X[:d]+" or "isoform [:d]+"
    acc_line = re.sub("isoform [0-9, A-Z]+ ", "", acc_line)
    acc_line = re.sub("isoform X[0-9, A-Z]+ ", "", acc_line)
    # This last step is nasty. These are the same gene:
    # >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # as the name is the same and they same 'isoform ...'
    # Whereas these are not the same gene, even though the names are identical
    # because they don't say isoform:
    # >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
    # >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
    #
    # To deal with that, we remove the ID (e.g. XP_024356342.1) if it says 'isoform'
    # so that the lines are identical, but not when it doesn't say 'isoform'
    # so that the lines are different. If I were writting the script from scratch
    # for NCBI files I'd do it a different way, but this is a way to handle it so 
    # that it works with the existing logic in the file.
    if original != acc_line:
        acc_line = acc_line.split(None, 1)[-1]
    return acc_line

def CreatePrimaryTranscriptsFile(fn, dout, gene_name_fn, q_use_original_accession_line):
    # Get genes and lengths
    max_gene_lens = defaultdict(int)
    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1
    nAcc = 0
    nGeneUnidentified = 0
    acc_to_use = defaultdict(str)
    iLine = -1
    while iLine < N:
        iLine += 1
        line = lines[iLine]
        if not line.startswith(">"): continue
        nAcc += 1
        iLineAcc = iLine
        gene = gene_name_fn(line)
        if gene is None:
            nGeneUnidentified += 1
            continue
        # get length
        l = 0
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if line.startswith(">"):
                iLine -= 1
                break
            l += len(line.rstrip())
        if l > max_gene_lens[gene]:
            max_gene_lens[gene] = l
            acc_to_use[gene] = iLineAcc
    print("Found %d accessions, %d genes, %d unidentified transcripts" % (nAcc, len(max_gene_lens), nGeneUnidentified))
    # print(gene)
    # print(sorted(max_gene_lens.keys())[:10])
    # print(len(set(max_gene_lens.keys())))

    # Get longest version for each gene
    # Parse file second time and only write out sequences that are longest variant
    nGenesWritten = 0
    outfn = os.path.join(dout, os.path.basename(fn))
    with open(outfn, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"): continue
            gene = gene_name_fn(line)
            # transcripts not identifying the gene should be written
            if gene != None and iLine != acc_to_use[gene]: continue
            if q_use_original_accession_line or gene == None:
                acc_line_out = line + "\n"
            else:
                 acc_line_out = ">%s\n" % gene
            nGenesWritten += 1
            outfile.write(acc_line_out)
            while iLine < N:
                iLine += 1
                line = lines[iLine]
                if line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(line + "\n")
    print("Wrote %d genes" % nGenesWritten)
    if nGenesWritten != len(max_gene_lens) + nGeneUnidentified:
        print("ERROR")
        raise Exception
    print(outfn)


def last_dot(text):
    return text[1:].rstrip().rsplit(".", 1)[0]

def space(text):
    return text[1:].rstrip().split(None, 1)[0]

def last_dot_before_first_space(text):
    return text[1:].rstrip().split(None, 1)[0].rstrip().rsplit(".", 1)[0]

def process_input(path, dout, gene_name_function=None):

    if os.path.isdir(path):

        for root, _, files in os.walk(path):
            for filename in sorted(files):
                process_input(
                    os.path.join(root, filename),
                    dout,
                    gene_name_function,
                )

    elif is_archive(path):

        process_archive(path, dout, gene_name_function)

    elif is_fasta_file(path):

        print("\nProcessing FASTA: %s" % path)
        process_fasta_file(path, dout, gene_name_function)

def main(args=None):

    function_dict = {"last_dot":last_dot, "space":space, "last_dot_before_first_space":last_dot_before_first_space}

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="Input FASTA file, directory, or archive",
    )

    parser.add_argument(
        "transcript_id_mode",
        nargs="?",
        metavar="mode",
        choices=sorted(function_dict),
        help="Transcript grouping mode",
    )

    parsed = parser.parse_args(args)

    fn = parsed.input

    gene_name_function = None
    if parsed.transcript_id_mode:
        gene_name_function = function_dict[parsed.transcript_id_mode]

    if os.path.isdir(fn):
        dout = os.path.join(os.path.abspath(fn), "primary_transcripts")
    else:
        dout = os.path.join(
            os.path.dirname(os.path.abspath(fn)),
            "primary_transcripts",
        )

    os.makedirs(dout, exist_ok=True)

    process_input(fn, dout, gene_name_function)

if __name__ == "__main__":
    main()