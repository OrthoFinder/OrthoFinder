"""
With --orthoxml: the orthogroups in OrthoXML format (http://orthoxml.org), for
other tools: Orthogroups/Orthogroups.orthoxml, written at the end of the run
from its final results files (write_results). With the HOGs (when the run
infers them), each orthogroup has its hierarchical orthogroups nested in it,
the species tree as the taxonomy, as OMA and orthoxml-tools write and read
HOGs; without, e.g. with -og, it is a group for each orthogroup.

(To look at the orthogroups and HOGs in a web browser: --groups-view,
orthology/groups_page.py.)

The species information comes from the input FASTA files. OrthoXML requires an
NCBI taxonomy ID and a database for each species, which OrthoFinder does not
know: the taxonomy ID is 0 (unknown) and the database is the FASTA file.
"""

import os
from collections import defaultdict
from xml.sax.saxutils import escape, quoteattr

from .. import __version__
from ..phylo import tree as tree_lib
from . import file_io
from .gene_ids import clean_species_name

GENE_SEP = ", "

NAMESPACE = "http://orthoXML.org/2011/"
UNKNOWN_TAX_ID = "0"
UNKNOWN_VERSION = "NA"


def species_for_analysis(species_ids_fn, species_to_use):
    """
    (species name, FASTA filename) for each species in species_to_use, from
    SpeciesIDs.txt ("<ID>: <FASTA filename>" lines). The species name is that of
    the other results files (clean_species_name).
    """
    fasta_files = {}
    with file_io.open_text_any(species_ids_fn) as infile:
        for line in infile:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            i_species, fasta_filename = line.split(": ", 1)
            fasta_files[int(i_species)] = fasta_filename
    return [(clean_species_name(fasta_files[i]), fasta_files[i]) for i in species_to_use]


# --- Orthogroups ---------------------------------------------------------------

def write_results(results_dir, species):
    """
    Write Orthogroups/Orthogroups.orthoxml, once, at the end of the run, from its
    final results files. Returns the file, or None if the results have no
    Orthogroups/Orthogroups.tsv (e.g. -ogt, -ost).

    With the HOGs (the run inferred them): the orthogroups with their HOGs
    nested in them (hog_groups). Without (e.g. -og): a group for each orthogroup
    of Orthogroups.tsv (the MCL orthogroups with -og). Single-gene orthogroups
    are not groups.

    species: [(species name, FASTA file)], from species_for_analysis
    """
    og_dir = os.path.join(results_dir, "Orthogroups")
    og_fn = os.path.join(og_dir, "Orthogroups.tsv")
    if not os.path.isfile(og_fn):
        return None
    out_fn = os.path.join(og_dir, "Orthogroups.orthoxml")
    gene_ids = read_genes([og_fn, os.path.join(og_dir, "Orthogroups_UnassignedGenes.tsv")])
    hogs = hog_groups(results_dir, gene_ids)
    if hogs is None:
        return write_flat(og_fn, gene_ids, out_fn, species)
    top_level, species_tree, taxon_ids = hogs
    write_orthoxml(out_fn, species_genes_of(gene_ids, species), top_level,
                   notes="The orthogroups of Orthogroups.tsv, each with its hierarchical "
                         "orthogroups (HOGs) nested in it.",
                   version=HOGS_VERSION, taxonomy=(species_tree, taxon_ids), score_defs=[COMPLETENESS])
    return out_fn


def write_flat(og_fn, gene_ids, out_fn, species):
    """
    out_fn: a group for each orthogroup of og_fn (an Orthogroups.tsv-like table,
    e.g. MCL_Orthogroups.tsv). Single-gene orthogroups are not groups.
    gene_ids: read_genes of the table and its unassigned genes.
    """
    groups = []
    for og_id, _, sp_genes in read_table(og_fn, 1, gene_ids):
        genes = sorted(g for ids in sp_genes.values() for g in ids)
        if len(genes) >= 2:
            groups.append((og_id, genes, ()))
    write_orthoxml(out_fn, species_genes_of(gene_ids, species), groups)
    return out_fn


def species_genes_of(gene_ids, species):
    """[(species name, FASTA file, [(gene id, gene name)])] for write_orthoxml."""
    fasta = dict(species)
    return [(name, fasta.get(name, name), sorted((i, g) for g, i in genes.items()))
            for name, genes in gene_ids.items()]


def read_genes(tables):
    """
    {species name: {gene name: gene id}}, the species in the column order of the
    first table: every gene of the tables (Orthogroups.tsv-like: an id column,
    then one per species), e.g. Orthogroups.tsv and Orthogroups_UnassignedGenes.tsv,
    numbered across species.
    """
    genes = {}
    for path in tables:
        if not os.path.isfile(path):
            continue
        with file_io.open_text(path) as infile:
            reader = file_io.reader(infile)
            sp_cols = next(reader)[1:]
            for sp in sp_cols:
                genes.setdefault(sp, {})
            for row in reader:
                for sp, cell in zip(sp_cols, row[1:]):
                    for g in split_genes(cell):
                        genes[sp].setdefault(g, None)
    n = 0
    for sp_genes in genes.values():
        for g in sp_genes:
            sp_genes[g] = n
            n += 1
    return genes


def read_table(fn, n_id_cols, gene_ids):
    """
    [(group, gene tree parent clade, {species: [gene id]})] of Orthogroups.tsv
    (n_id_cols 1) or an N<i>.tsv (n_id_cols 3: HOG, OG, gene tree parent clade).
    """
    with file_io.open_text(fn) as infile:
        reader = file_io.reader(infile)
        sp_cols = next(reader)[n_id_cols:]
        return [(row[0], row[2] if n_id_cols == 3 else None,
                 {sp: [gene_ids[sp][g] for g in split_genes(cell)]
                  for sp, cell in zip(sp_cols, row[n_id_cols:]) if cell})
                for row in reader]


def split_genes(cell):
    return cell.split(GENE_SEP) if cell else []


# --- With the HOGs --------------------------------------------------------------------

HOGS_VERSION = "0.5"  # the first version with a taxonomy
COMPLETENESS = ("CompletenessScore",
                ("The fraction of the species of the group's taxon (its node of the species tree) "
                 "with a gene in the group"))


def hog_groups(results_dir, gene_ids):
    """
    The orthogroups with their HOGs nested in them, from Orthogroups.tsv,
    Phylogenetic_Hierarchical_Orthogroups/N<i>.tsv and
    Species_Tree/SpeciesTree_rooted_node_labels.txt: (top-level groups, species
    tree, {node: taxon id}), or None if the run did not infer the HOGs.

    A top-level group for each orthogroup (the HOGs of N0), its id the
    orthogroup. In the group of a HOG at a node, for each child of the node:
      - a HOG of the child node within it: its group; several: a paralogGroup
        of their groups (a duplication on the branch to the child)
      - a species' gene: the gene; several: a paralogGroup of them
    A gene of a HOG in no HOG of the child node of its species (e.g. an
    orthogroup without a gene tree, or a copy whose other genes were lost) is a
    copy at that child, as its HOGs are: several such genes of a species are a
    paralogGroup (never orthologs of each other), and those of several species an
    orthologGroup. These are OrthoFinder's orthologs (Orthologs/), checked with
    orthoxml-tools export-pairs, except where the gene tree and the HOGs differ
    (a gene placed against the species tree). A group of a single member is
    replaced by that member (an OrthoXML group has two or more), e.g. a HOG of
    the genes of one child: at the top level, that member then has the id of
    the orthogroup, and its own as the property "HOG". A paralogGroup is not a
    top-level group: a root HOG is always an orthologGroup, as OMA and
    orthoxml-tools (e.g. to-nhx) expect, although the OrthoXML schema requires
    two members: orthoxml-tools validate reports these groups, but its converters
    (to-csv, to-nhx) skip a top-level paralogGroup.

    Each orthologGroup has, as in OMA: the taxonId of its node in the taxonomy
    (the species tree), the node's name as the property TaxRange, and the score
    CompletenessScore.
    """
    species_tree_fn = os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    hog_dir = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    if not os.path.isfile(species_tree_fn) or not os.path.isdir(hog_dir):
        return None
    species_tree = tree_lib.Tree(species_tree_fn, format=1)
    preorder = list(species_tree.traverse("preorder"))
    taxon_ids = {node: i for i, node in enumerate(preorder, 1)}
    n_species = {node: len(node.get_leaf_names()) for node in preorder if not node.is_leaf()}

    def hog_table(node):
        """The HOGs at node (read_table)."""
        if node.is_root():
            return read_table(os.path.join(results_dir, "Orthogroups", "Orthogroups.tsv"), 1, gene_ids)
        fn = os.path.join(hog_dir, node.name + ".tsv")
        return read_table(fn, 3, gene_ids) if os.path.isfile(fn) else []

    # 1. What is in each HOG (node, i): the HOGs of the child nodes in it, and its
    # other genes (with the child species they are from, or None), and its species
    hogs = {}                    # (node, i) -> (HOG id, gene tree parent clade, {species})
    n_hogs = {}                  # node -> number of HOGs
    kids = defaultdict(list)     # (node, i) -> [(child node, j)]
    genes = defaultdict(list)    # (node, i) -> [(child node of its species, species, gene id)]: in no HOG of the child
    tables = {species_tree: hog_table(species_tree)}
    for node in preorder:
        if node.is_leaf():
            continue
        table = tables.pop(node)
        n_hogs[node] = len(table)
        hog_of_gene = {}
        for i, (hog_id, parent_clade, sp_genes) in enumerate(table):
            hogs[(node, i)] = (hog_id, parent_clade, set(sp_genes))
            hog_of_gene.update((g, i) for ids in sp_genes.values() for g in ids)
        for child in node.children:
            if child.is_leaf():
                continue
            tables[child] = hog_table(child)
            for j, (_, _, sp_genes) in enumerate(tables[child]):
                child_genes = [g for ids in sp_genes.values() for g in ids]
                kids[(node, hog_of_gene[child_genes[0]])].append((child, j))
                for g in child_genes:
                    del hog_of_gene[g]
        # (the child node each species is in: the species itself, or the clade it is in)
        species_child = {leaf: child for child in node.children for leaf in child.get_leaf_names()}
        for i, (_, _, sp_genes) in enumerate(table):
            for sp, ids in sp_genes.items():
                genes[(node, i)].extend((species_child.get(sp), sp, g) for g in ids if g in hog_of_gene)

    # 2. The group of each HOG, children first: ("gene", id) or
    # (tag, id, taxonId, [(score, value)], [(property, value)], [members])
    built = {}
    for node in reversed(preorder):
        if node.is_leaf():
            continue
        for key in ((node, i) for i in range(n_hogs[node])):
            members = []
            key_kids, key_genes = kids.pop(key, []), genes.pop(key, [])
            for child in node.children:
                parts = [built.pop(k) for k in key_kids if k[0] is child]
                # The genes of the child's species in no HOG of the child: each species'
                # genes one copy (several: a paralogGroup, never orthologs of each other),
                # a copy at the child as a HOG of it is (e.g. one whose other genes were
                # lost, or an orthogroup without a gene tree); the species' copies
                # orthologs of each other, as OrthoFinder has them without a gene tree
                by_species = defaultdict(list)
                for c, sp, g in key_genes:
                    if c is child:
                        by_species[sp].append(("gene", g))
                copies = [ids[0] if len(ids) == 1 else ("paralogGroup", None, None, [], [], ids)
                          for ids in by_species.values()]
                if len(copies) == 1:
                    parts.append(copies[0])
                elif copies:
                    parts.append(("orthologGroup", None, taxon_ids[child], [], [("TaxRange", child.name)], copies))
                if len(parts) == 1:
                    members.append(parts[0])
                elif parts:
                    members.append(("paralogGroup", None, None, [], [], parts))
            hog_id, parent_clade, hog_species = hogs.pop(key)
            scores = [(COMPLETENESS[0], round(len(hog_species) / n_species[node], 4))]
            properties = [("TaxRange", node.name)]
            if not node.is_root():
                properties.append(("gene tree parent clade", parent_clade))
            if len(members) > 1 or (node.is_root() and members[0][0] == "paralogGroup"):
                built[key] = ("orthologGroup", hog_id, taxon_ids[node], scores, properties, members)
            elif node.is_root() and members[0][0] == "orthologGroup":
                _, child_hog_id, taxon_id, child_scores, child_properties, sub = members[0]
                built[key] = ("orthologGroup", hog_id, taxon_id, child_scores,
                              child_properties + [("HOG", child_hog_id)], sub)
            else:
                built[key] = members[0]

    top_level = [built[(species_tree, i)] for i in range(n_hogs[species_tree])
                 if built[(species_tree, i)][0] != "gene"]
    return top_level, species_tree, taxon_ids


# --- writing -----------------------------------------------------------------

def write_orthoxml(filename, species_genes, groups, notes=None, version="0.3", taxonomy=None, score_defs=()):
    """
    species_genes: [(species name, FASTA file, [(gene id, gene name)])]
    groups: [(group id, [gene ids], [(property name, value)])], or (nested, with
            version 0.5) [(tag, group id, taxonId, [(score id, value)], [(property name, value)],
            [member])], a member a group or ("gene", gene id)
    taxonomy: (species tree, {node: taxon id}), with version 0.5
    score_defs: [(score id, description)]
    """
    with open(filename, "w", encoding="utf-8", newline="\n") as out:
        out.write('<?xml version="1.0" encoding="utf-8"?>\n')
        # (no xsi:schemaLocation, as in FastOMA: orthoxml-tools filter etc. write it wrongly)
        out.write(f'<orthoXML xmlns={quoteattr(NAMESPACE)} origin="OrthoFinder" '
                  f'originVersion={quoteattr(__version__)} version="{version}">\n')
        if notes:
            out.write(f"  <notes>{escape(notes)}</notes>\n")
        leaf_taxa = {node.name: i for node, i in taxonomy[1].items() if node.is_leaf()} if taxonomy else {}
        for name, fasta_filename, genes in species_genes:
            taxon = f' taxonId="{leaf_taxa[name]}"' if name in leaf_taxa else ""
            out.write(f'  <species name={quoteattr(name)} NCBITaxId="{UNKNOWN_TAX_ID}"{taxon}>\n')
            out.write(f'    <database name={quoteattr(fasta_filename)} version="{UNKNOWN_VERSION}">\n')
            out.write("      <genes>\n")
            out.writelines(f'        <gene id="{i}" geneId={quoteattr(g)}/>\n' for i, g in genes)
            out.write("      </genes>\n")
            out.write("    </database>\n")
            out.write("  </species>\n")
        if taxonomy:
            write_taxonomy(out, *taxonomy)
        if score_defs:
            out.write("  <scores>\n")
            out.writelines(f"    <scoreDef id={quoteattr(i)} desc={quoteattr(d)}/>\n" for i, d in score_defs)
            out.write("  </scores>\n")
        out.write("  <groups>\n")
        for group in groups:
            if len(group) == 3:  # flat
                group_id, genes, properties = group
                group = ("orthologGroup", group_id, None, [], properties, [("gene", i) for i in genes])
            write_group(out, group, 2)
        out.write("  </groups>\n")
        out.write("</orthoXML>\n")


def write_taxonomy(out, species_tree, taxon_ids):
    """The species tree as the taxonomy (a taxon for each node)."""
    out.write("  <taxonomy>\n")
    stack = [(species_tree, 2)]
    while stack:
        node, depth = stack.pop()
        if node is None:            # closes a taxon
            out.write("  " * depth + "</taxon>\n")
            continue
        pad = "  " * depth
        if node.is_leaf():
            out.write(f'{pad}<taxon id="{taxon_ids[node]}" name={quoteattr(node.name)}/>\n')
            continue
        out.write(f'{pad}<taxon id="{taxon_ids[node]}" name={quoteattr(node.name)}>\n')
        stack.append((None, depth))
        stack.extend((child, depth + 1) for child in reversed(node.children))
    out.write("  </taxonomy>\n")


def write_group(out, group, depth):
    """A group and its members, nested (without recursion: species trees can be deep)."""
    stack = [(group, depth)]
    while stack:
        item, depth = stack.pop()
        pad = "  " * depth
        if isinstance(item, str):   # closes a group
            out.write(f"{pad}</{item}>\n")
        elif item[0] == "gene":
            out.write(f'{pad}<geneRef id="{item[1]}"/>\n')
        else:
            tag, group_id, taxon_id, scores, properties, members = item
            attrs = (f" id={quoteattr(group_id)}" if group_id is not None else "") + \
                    (f' taxonId="{taxon_id}"' if taxon_id is not None else "")
            out.write(f"{pad}<{tag}{attrs}>\n")
            out.writelines(f'{pad}  <score id={quoteattr(i)} value="{v}"/>\n' for i, v in scores)
            out.writelines(f"{pad}  <property name={quoteattr(p)} value={quoteattr(v)}/>\n" for p, v in properties)
            stack.append((tag, depth))
            stack.extend((m, depth + 1) for m in reversed(members))
