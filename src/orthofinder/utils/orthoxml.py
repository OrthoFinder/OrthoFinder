"""
With --orthoxml: the orthogroups and HOGs in OrthoXML format
(http://orthoxml.org, version 0.3), for other tools:

    Orthogroups/Orthogroups.orthoxml
    Phylogenetic_Hierarchical_Orthogroups/N<i>.orthoxml (one per N<i>.tsv)

(To look at the orthogroups and HOGs in a web browser: --groups-view,
orthology/groups_page.py.)

A gene has the same id in all the OrthoXML files. The species information
comes from the input FASTA files. OrthoXML requires an NCBI taxonomy ID and a
database for each species, which OrthoFinder does not know: the taxonomy ID
is 0 (unknown) and the database is the FASTA file.
"""

import csv
import glob
import os
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape, quoteattr

from .. import __version__

GENE_SEP = ", "

NAMESPACE = "http://orthoXML.org/2011/"
SCHEMA_LOCATION = "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd"
UNKNOWN_TAX_ID = "0"
UNKNOWN_VERSION = "NA"


def species_for_analysis(species_ids_fn, species_to_use):
    """
    (species name, FASTA filename) for each species in species_to_use, from
    SpeciesIDs.txt ("<ID>: <FASTA filename>" lines). The species name is the
    FASTA filename without its extension, as in the other results files.
    """
    fasta_files = {}
    with open(species_ids_fn) as infile:
        for line in infile:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            i_species, fasta_filename = line.split(": ", 1)
            fasta_files[int(i_species)] = fasta_filename
    return [(fasta_files[i].rsplit(".", 1)[0], fasta_files[i]) for i in species_to_use]


# --- Orthogroups ---------------------------------------------------------------

def write_results(species, ogs, n_seqs_per_species, id_dict, species_to_use, results_base_filename):
    """
    Write <results_base_filename>.orthoxml.

    species: from species_for_analysis, one entry per species in species_to_use
    ogs: the orthogroups, sets of sequence IDs ("<species>_<sequence>"); orthogroup
         i is OG<i> as in Orthogroups.tsv. Single-gene orthogroups are not groups.
    n_seqs_per_species: {species ID: number of sequences}
    id_dict: {sequence ID: gene name}
    """
    # every gene of every species, numbered across species in this order
    gene_number, species_genes = {}, []
    for i_species, (name, fasta_filename) in zip(species_to_use, species):
        genes = []
        for i_seq in range(n_seqs_per_species[i_species]):
            seq_id = f"{i_species}_{i_seq}"
            gene_number[seq_id] = len(gene_number)
            genes.append((gene_number[seq_id], id_dict[seq_id]))
        species_genes.append((name, fasta_filename, genes))
    multi = [(i_og, og) for i_og, og in enumerate(ogs) if len(og) >= 2]
    groups = ((f"OG{i_og:07d}", sorted(gene_number[s] for s in og), ()) for i_og, og in multi)
    write_orthoxml(results_base_filename + ".orthoxml", species_genes, groups)


# --- HOGs, from the results files ---------------------------------------------------

def write_hog_files(orthogroups_orthoxml, hog_dir):
    """
    From the final results files: an OrthoXML file for each
    Phylogenetic_Hierarchical_Orthogroups/N<i>.tsv. The species and the gene
    ids are those of Orthogroups.orthoxml. Returns the files written.
    """
    if not hog_dir or not os.path.isdir(hog_dir):
        return []
    species, gene_ids = read_genes(orthogroups_orthoxml)
    return write_hogs(species, gene_ids, hog_dir)


def read_genes(orthoxml_fn):
    """[(species name, FASTA file)], {species name: {gene name: gene id}} of an OrthoXML file."""
    species, gene_ids = [], {}
    for _, el in ET.iterparse(orthoxml_fn):
        if el.tag == f"{{{NAMESPACE}}}species":
            db = el.find(f"{{{NAMESPACE}}}database")
            name = el.get("name")
            species.append((name, db.get("name")))
            gene_ids[name] = {g.get("geneId"): int(g.get("id")) for g in db.iter(f"{{{NAMESPACE}}}gene")}
            el.clear()
    return species, gene_ids


def split_genes(cell):
    return cell.split(GENE_SEP) if cell else []


def node_number(fn):
    n = os.path.basename(fn)[1:-4]
    return int(n) if n.isdigit() else -1


def write_hogs(species, gene_ids, hog_dir):
    written = []
    for fn in sorted(glob.glob(os.path.join(hog_dir, "N*.tsv")), key=node_number):
        node = os.path.basename(fn)[:-4]
        with open(fn, newline="") as infile:
            table = list(csv.reader(infile, delimiter="\t"))
        header, rows = table[0], table[1:]
        if not rows:  # OrthoXML needs at least one group
            continue
        sp_cols = header[3:]
        # the species of this node's clade: those with genes in its HOGs
        clade = [j for j in range(len(sp_cols)) if any(row[3 + j] for row in rows)]
        clade_names = [sp_cols[j] for j in clade]
        groups = []
        for row in rows:
            genes = sorted(gene_ids[sp_cols[j]][g] for j in clade for g in split_genes(row[3 + j]))
            if len(genes) >= 2:  # an OrthoXML group has at least two genes
                groups.append((row[0], genes, (("orthogroup", row[1]), ("gene tree parent clade", row[2]))))
        if groups:
            out = fn[:-4] + ".orthoxml"
            write_orthoxml(
                out, genes_of(species, gene_ids, groups, clade_names), groups,
                notes=f"Hierarchical orthogroups (HOGs) at species tree node {node} "
                      f"({', '.join(clade_names)}), as in {os.path.basename(fn)}, "
                      f"except the HOGs of a single gene.")
            written.append(out)
    return written


def genes_of(species, gene_ids, groups, species_names):
    """(name, FASTA file, [(gene id, gene name)]) of the given species, for the genes in the groups."""
    used = {g for _, genes, _ in groups for g in genes}
    fasta = dict(species)
    out = []
    for name in species_names:
        genes = sorted((i, g) for g, i in gene_ids[name].items() if i in used)
        if genes:
            out.append((name, fasta[name], genes))
    return out


# --- writing -----------------------------------------------------------------

def write_orthoxml(filename, species_genes, groups, notes=None):
    """
    species_genes: [(species name, FASTA file, [(gene id, gene name)])]
    groups: [(group id, [gene ids], [(property name, value)])]
    """
    with open(filename, "w") as out:
        out.write('<?xml version="1.0" encoding="utf-8"?>\n')
        out.write(f'<orthoXML xmlns={quoteattr(NAMESPACE)} xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation={quoteattr(SCHEMA_LOCATION)} origin="OrthoFinder" '
                  f'originVersion={quoteattr(__version__)} version="0.3">\n')
        if notes:
            out.write(f"  <notes>{escape(notes)}</notes>\n")
        for name, fasta_filename, genes in species_genes:
            out.write(f'  <species name={quoteattr(name)} NCBITaxId="{UNKNOWN_TAX_ID}">\n')
            out.write(f'    <database name={quoteattr(fasta_filename)} version="{UNKNOWN_VERSION}">\n')
            out.write("      <genes>\n")
            out.writelines(f'        <gene id="{i}" geneId={quoteattr(g)}/>\n' for i, g in genes)
            out.write("      </genes>\n")
            out.write("    </database>\n")
            out.write("  </species>\n")
        out.write("  <groups>\n")
        for group_id, genes, properties in groups:
            out.write(f"    <orthologGroup id={quoteattr(group_id)}>\n")
            out.writelines(f"      <property name={quoteattr(p)} value={quoteattr(v)}/>\n" for p, v in properties)
            out.writelines(f'      <geneRef id="{i}"/>\n' for i in genes)
            out.write("    </orthologGroup>\n")
        out.write("  </groups>\n")
        out.write("</orthoXML>\n")
