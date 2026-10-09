"""
The orthogroups in OrthoXML format (http://orthoxml.org, version 0.3), written
when the user gives -x <file>.

The file gives, for each species, one line of five tab-separated fields:
    FASTA filename, species name, NCBI taxonomy ID, source database, database version
e.g.
    Mycoplasma_genitalium.faa	Mycoplasma genitalium	2097	Ensembl	Ensembl 46
"""

import os
from xml.sax.saxutils import quoteattr

from .. import __version__

N_FIELDS = 5
FIELDS = "FASTA filename, species name, NCBI taxonomy ID, source database, database version"
NAMESPACE = "http://orthoXML.org/2011/"
SCHEMA_LOCATION = "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd"


class SpeciesFileError(Exception):
    pass


def species_key(fasta_filename):
    """The species name OrthoFinder uses for a FASTA file (as in SpeciesIDs.txt, without the extension)."""
    return os.path.basename(fasta_filename).rsplit(".", 1)[0]


def read_species_file(species_file):
    """
    The lines of a -x file, as {species key: (species name, NCBI taxonomy ID,
    database, database version)}. Raises SpeciesFileError if the file cannot
    be read or a line is not in the format.
    """
    try:
        with open(species_file) as infile:
            lines = infile.read().splitlines()
    except OSError as e:
        raise SpeciesFileError(f"Cannot read the OrthoXML species file {species_file}: {e.strerror}")
    species = {}
    for i_line, line in enumerate(lines, 1):
        if not line.strip():
            continue
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) != N_FIELDS or not all(fields):
            raise SpeciesFileError(
                f"Line {i_line} of the OrthoXML species file {species_file} is not in the format.\n"
                f"Each line should have {N_FIELDS} tab-separated fields: {FIELDS}")
        fasta_filename, name, tax_id, database, version = fields
        if not tax_id.isdigit():
            raise SpeciesFileError(
                f"Line {i_line} of the OrthoXML species file {species_file}: "
                f"the NCBI taxonomy ID should be a number, not '{tax_id}'")
        key = species_key(fasta_filename)
        if key in species:
            raise SpeciesFileError(
                f"The OrthoXML species file {species_file} has more than one line for {fasta_filename}")
        species[key] = (name, tax_id, database, version)
    if not species:
        raise SpeciesFileError(f"The OrthoXML species file {species_file} is empty")
    return species


def species_for_analysis(species_file, species_names_dict, species_to_use):
    """
    The -x information for each species in the analysis, in the order of
    species_to_use. species_names_dict: {species ID: name} (from SpeciesIDs.txt).
    Raises SpeciesFileError if a species is missing from the file.
    """
    species = read_species_file(species_file)
    missing = [species_names_dict[i] for i in species_to_use if species_names_dict[i] not in species]
    if missing:
        missing_list = "\n  ".join(missing)
        raise SpeciesFileError(f"The OrthoXML species file {species_file} has no line for:\n  {missing_list}")
    return [species[species_names_dict[i]] for i in species_to_use]


def write_orthoxml(species_xml, ogs, n_seqs_per_species, id_dict, species_to_use, filename):
    """
    Write the orthogroups to filename in OrthoXML.

    species_xml: from species_for_analysis, one entry per species in species_to_use
    ogs: the orthogroups, sets of sequence IDs ("<species>_<sequence>"); orthogroup
         i is OG<i> as in Orthogroups.tsv. Single-gene orthogroups are not groups.
    n_seqs_per_species: {species ID: number of sequences}
    id_dict: {sequence ID: gene name}
    """
    gene_number = {}  # sequence ID -> the gene's id in the file (unique across species)
    with open(filename, "w") as out:
        out.write('<?xml version="1.0" encoding="utf-8"?>\n')
        out.write(f'<orthoXML xmlns={quoteattr(NAMESPACE)} xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation={quoteattr(SCHEMA_LOCATION)} origin="OrthoFinder" '
                  f'originVersion={quoteattr(__version__)} version="0.3">\n')
        for i_species, (name, tax_id, database, version) in zip(species_to_use, species_xml):
            out.write(f"  <species name={quoteattr(name)} NCBITaxId={quoteattr(tax_id)}>\n")
            out.write(f"    <database name={quoteattr(database)} version={quoteattr(version)}>\n")
            out.write("      <genes>\n")
            for i_seq in range(n_seqs_per_species[i_species]):
                seq_id = f"{i_species}_{i_seq}"
                gene_number[seq_id] = len(gene_number)
                out.write(f'        <gene id="{gene_number[seq_id]}" geneId={quoteattr(id_dict[seq_id])}/>\n')
            out.write("      </genes>\n")
            out.write("    </database>\n")
            out.write("  </species>\n")
        out.write("  <groups>\n")
        for i_og, og in enumerate(ogs):
            if len(og) < 2:
                continue
            out.write(f'    <orthologGroup id="OG{i_og:07d}">\n')
            out.writelines(f'      <geneRef id="{i_gene}"/>\n' for i_gene in sorted(gene_number[seq_id] for seq_id in og))
            out.write("    </orthologGroup>\n")
        out.write("  </groups>\n")
        out.write("</orthoXML>\n")
