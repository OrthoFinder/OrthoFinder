# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologs
"""
import os
import sys
import itertools
from collections import defaultdict 

import multiprocessing as mp

try:
    from rich import print
except ImportError:
    ...
from ..utils import file_handler, file_io
from ..phylo import taxa as phylo_taxa
from . import stats
from .rooting import check_and_root_tree
from .resolve import resolve_tree


debug = False   # HOGs


def gene_to_species_dash(g):
  return g.split("_", 1)[0]
  
ortho_finder_ids = gene_to_species_dash


def species_and_gene_dash(g):
  return g.split("_", 1)

  
species_and_gene_lookup = {
    gene_to_species_dash:species_and_gene_dash, 
}


class TreeAnalyser(object):
    def __init__(
            self, 
            n_ogs, 
            d_results_orthologs, 
            recon_trees_renamed_dir, 
            species_tree_rooted_labelled, 
            species_to_use, 
            GeneToSpecies, 
            sequence_dict, 
            species_dict, 
            spec_seq_dict, 
            neighbours, 
            q_no_recon, 
            dups_file_handle, 
            stride_dups, 
            orthologs_files_handles, 
            putative_xenolog_file_handles, 
            hog_writer, 
            q_split_paralogous_clades, 
            fewer_open_files=False,
            exist_msa=True,
            write_hog_tree=True,
            fix_files=True,
    ):
        self.n_ogs = n_ogs
        self.d_results_orthologs = d_results_orthologs
        self.recon_trees_renamed_dir = recon_trees_renamed_dir
        self.species_tree_rooted_labelled = species_tree_rooted_labelled
        self.species_to_use = species_to_use
        self.nspecies = len(self.species_to_use)
        self.GeneToSpecies = GeneToSpecies
        self.sequence_dict = sequence_dict
        self.species_dict = species_dict
        self.spec_seq_dict = spec_seq_dict
        self.spec_seq_dict[">%d genes" % (5*self.nspecies)] = ">%d genes" % (5*self.nspecies)   # used in Duplications
        self.neighbours = neighbours
        self.q_no_recon = q_no_recon
        self.dups_file_handle = dups_file_handle
        self.stride_dups = stride_dups
        self.orthologs_files_handles = orthologs_files_handles
        self.putative_xenolog_file_handles = putative_xenolog_file_handles
        self.hog_writer = hog_writer
        self.q_split_paralogous_clades = q_split_paralogous_clades
        self.fewer_open_files = fewer_open_files
        self.lock_orthologs = [mp.Lock() for i in range(self.nspecies)]   # lock the larger of the two species index
        self.lock_dups = mp.Lock()
        self.lock_suspect = mp.Lock()
        self.lock_hogs = mp.Lock()
        self.exist_msa = exist_msa
        self.write_hog_tree = write_hog_tree
        self.fix_files = fix_files


    def analyse_tree(self, iog, n_ortho_acc=None):
        """
        Analyse one gene tree.

        It returns a payload to the parent process.

        If nOrtho_acc (a stats.NOrthoSp) is given, this tree's ortholog
        counts are added to it once the whole tree has been analysed, and the
        payload's "n_orthologues" is None. Otherwise it holds a dense
        NOrthoSp for this tree.
        """
        og_name = "OG%07d" % iog
        n_species = len(self.species_to_use)

        try:
            if self.write_hog_tree or not self.fix_files:
                tree_file = file_handler.handler.get_ogs_tree_fn(iog)
            else:
                tree_file = file_handler.handler.get_resolved_tree_id_dir() + "OG%07d.txt" % iog

            if not os.path.exists(tree_file):
                return None

            rooted_tree_ids, q_have_support = check_and_root_tree(
                tree_file,
                self.species_tree_rooted_labelled,
                self.GeneToSpecies
            )

            if rooted_tree_ids is None:
                return None

            orthologs, recon_tree, suspect_genes, dups = get_orthologs_from_tree(
                iog,
                rooted_tree_ids,
                self.species_tree_rooted_labelled,
                self.GeneToSpecies,
                self.neighbours,
                q_get_dups=True,
                q_no_recon=self.q_no_recon
            )

            ortholog_lines = empty_ortholog_lines(self.nspecies, self.fewer_open_files)
            ortholog_sus_lines = ["" for _ in range(self.nspecies)]

            pair_counts = [] if n_ortho_acc is not None else None
            n_orthologs_sp_pair = get_lines_for_ortholog_files(
                [(iog, orthologs)],
                self.species_dict,
                self.species_to_use,
                self.sequence_dict,
                len(suspect_genes) > 0,
                ortholog_lines,
                ortholog_sus_lines,
                fewer_open_files=self.fewer_open_files,
                pair_counts=pair_counts,
            )

            cached_hogs = get_hogs_from_tree(
                iog,
                recon_tree,
                self.hog_writer,
                self.q_split_paralogous_clades,
            )

            phylo_taxa.rename_tree_taxa(
                recon_tree,
                self.recon_trees_renamed_dir + "OG%07d.txt" % iog,
                self.spec_seq_dict,
                q_support=False,
                q_fix_negatives=True
            )

            if n_ortho_acc is not None:
                # Only count trees that were analysed completely.
                for counts in pair_counts:
                    n_ortho_acc.add_pair(*counts)

            return {
                "iog": iog,
                "n_orthologues": n_orthologs_sp_pair,
                "olog_lines": ortholog_lines,
                "olog_sus_lines": ortholog_sus_lines,
                "duplications": dups,
                "suspect_genes": suspect_genes,
                "cached_hogs": cached_hogs,
            }

        except Exception as e:
            print(str(e))
            print("WARNING: Unknown error analysing tree %s" % og_name)

            return {
                "iog": iog,
                "n_orthologues": (
                    stats.NOrthoSp(n_species) if n_ortho_acc is None else None
                ),
                "olog_lines": empty_ortholog_lines(self.nspecies, self.fewer_open_files),
                "olog_sus_lines": ["" for _ in range(self.nspecies)],
                "duplications": [],
                "suspect_genes": set(),
                "cached_hogs": [],
            }

def empty_ortholog_lines(nspecies, fewer_open_files):
    """
    Per-species containers for one tree's ortholog rows.

    In pairwise mode each row is a dict {j: text} holding only the partner
    species that have orthologs, instead of an nspecies x nspecies grid
    (1,000,000 cells per tree at 1,000 species). Use iter_ortholog_row to read it.
    """
    if fewer_open_files:
        return [[""] for _ in range(nspecies)]
    return [{} for _ in range(nspecies)]


def iter_ortholog_row(row):
    """(j, text) for the non-empty cells of an olog_lines row (list or dict)."""
    items = row.items() if isinstance(row, dict) else enumerate(row)
    return ((j, text) for j, text in items if text)


class SpeciesSets(object):
    """
    Species below every node of a gene tree, computed once in one postorder pass.

    Each set is stored as an int bitmask (one bit per species), so memory is
    about nspecies/8 bytes per node rather than a Python set per node, and a
    union is a single OR. Valid while the tree topology does not change.
    """

    def __init__(self, tree, GeneToSpecies):
        self.species_names = []
        self.species_bit = {}
        self.masks = {}
        self.leaves = {}
        self._set_cache = {}
        for n in tree.traverse("postorder"):
            if n.is_leaf():
                sp = GeneToSpecies(n.name)
                bit = self.species_bit.get(sp)
                if bit is None:
                    bit = 1 << len(self.species_names)
                    self.species_bit[sp] = bit
                    self.species_names.append(sp)
                self.masks[id(n)] = bit
                self.leaves.setdefault(n.name, n)
            else:
                mask = 0
                for c in n.get_children():
                    mask |= self.masks[id(c)]
                self.masks[id(n)] = mask

    def mask(self, node):
        return self.masks[id(node)]

    def species(self, node):
        """A new set of the species below node (callers may modify it)."""
        mask = self.masks[id(node)]
        cached = self._set_cache.get(mask)
        if cached is None:
            names = self.species_names
            out = []
            m = mask
            while m:
                low = m & -m
                out.append(names[low.bit_length() - 1])
                m ^= low
            cached = frozenset(out)
            self._set_cache[mask] = cached
        return set(cached)

    def leaf(self, name):
        return self.leaves[name]


def get_orthologs_from_tree(
        iog,
        tree,
        species_tree_rooted,
        GeneToSpecies,
        neighbours,
        q_get_dups=False,
        q_no_recon=False,
    ):
    """Infer orthologs and duplications from a reconciled gene tree."""
    og_name = "OG%07d" % iog
    q_prune = False

    species_and_gene = species_and_gene_lookup[GeneToSpecies]

    orthologs = []
    duplications = []

    if not q_no_recon:
        tree = resolve_tree(tree, GeneToSpecies)

    if q_prune:
        tree.prune(tree.get_leaf_names())

    if len(tree) == 1:
        return set(orthologs), tree, set(), duplications

    i_node = 1
    tree.name = "n0"

    suspect_genes = set()
    empty_set = set()

    # The topology is fixed from here on, so species sets and leaf lookups can
    # be computed once instead of re-walking subtrees at every node.
    sp_sets = SpeciesSets(tree, GeneToSpecies)
    mrca_cache = {}

    def species_tree_mrca(mask, species):
        node = mrca_cache.get(mask)
        if node is None:
            node = mrca_node(species_tree_rooted, species)
            mrca_cache[mask] = node
        return node

    for n in tree.traverse("preorder"):
        if n.is_leaf():
            continue

        if not n.is_root():
            n.name = "n%d" % i_node
            i_node += 1

        ch = n.get_children()

        if len(ch) == 2:
            o_size, overlap, sp0, sp1 = overlap_size(
                n,
                GeneToSpecies,
                suspect_genes,
                sp_sets=sp_sets,
            )

            sp_present = sp0.union(sp1)
            if len(sp_present) == len(sp_sets.species(n)):
                # No species removed as suspect: the set is exactly the species below n.
                st_node = species_tree_mrca(sp_sets.mask(n), sp_present)
            else:
                st_node = mrca_node(species_tree_rooted, sp_present)

            n.add_feature("sp_node", st_node.name)

            if o_size != 0:
                q_resolved, misplaced_genes = resolve_overlap(
                    overlap,
                    sp0,
                    sp1,
                    ch,
                    tree,
                    neighbours,
                    GeneToSpecies,
                    sp_sets=sp_sets,
                )

                for g in misplaced_genes:
                    nn = sp_sets.leaf(g)
                    nn.add_feature("X", True)

            else:
                misplaced_genes = empty_set

            dup = o_size != 0 and not q_resolved

            n.add_feature("dup", dup)

            if dup:
                n.add_feature("dup_type", "species_overlap_duplication")

                if q_get_dups:
                    genes0 = ch[0].get_leaf_names()
                    genes1 = ch[1].get_leaf_names()

                    duplications.append((
                        st_node.name,
                        n.name,
                        float(o_size) / len(st_node),
                        genes0,
                        genes1
                    ))

            else:
                n.add_feature("dup_type", "speciation")

                orthologs.append(
                    orthologs_and_suspect(
                        ch,
                        suspect_genes,
                        misplaced_genes,
                        species_and_gene
                    )
                )

                suspect_genes.update(misplaced_genes)

        elif len(ch) > 2:
            species = [sp_sets.species(child) for child in ch]
            all_species = set.union(*species)
            st_node = species_tree_mrca(sp_sets.mask(n), all_species)
            n.add_feature("sp_node", st_node.name)

            if len(all_species) == 1:
                genes = n.get_leaf_names()
                duplications.append((st_node.name, n.name, 1.0, genes, []))
                n.add_feature("dup", True)
            else:
                dups = []
                for (n0, s0), (n1, s1) in itertools.combinations(
                        zip(ch, species), 2):
                    if len(s0.intersection(s1)) == 0:
                        orthologs.append(
                            orthologs_and_suspect(
                                (n0, n1),
                                suspect_genes,
                                empty_set,
                                species_and_gene
                            )
                        )
                        dups.append(False)
                    else:
                        dups.append(True)

                if all(dups):
                    genes = n.get_leaf_names()
                    duplications.append((st_node.name, n.name, 1.0, genes, []))

                n.add_feature("dup", all(dups))

    return orthologs, tree, suspect_genes, duplications


def mrca_node(t_rooted, taxa):
    return (t_rooted & next(taxon for taxon in taxa)) if len(taxa) == 1 else t_rooted.get_common_ancestor(taxa)


def overlap_size(node, GeneToSpecies, suspect_genes, sp_sets=None):
    if sp_sets is None:
        descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()}.difference(suspect_genes) for n in node.get_children()]
    else:
        descendents = [sp_sets.species(n).difference(suspect_genes) for n in node.get_children()]
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]

def resolve_overlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies, rel_overlap_cutoff=4, sp_sets=None):
    """
    Is an overlap suspicious and if so can it be resolved by identifying genes that are out of place?
    Args:
        overlap - the species with genes in both clades
        sp0 - the species below ch[0]
        sp1 - the species below ch[1]
        ch - the two child nodes
        tree - the gene tree
        neighbours - dictionary species->neighbours, where neighbours is a list of the sets of species observed at successive topological distances from the species
    Returns:
        qSuccess - has the overlap been resolved
        genes_removed - the out-of-place genes that have been removed so as to resolve the overlap
    
    Implementation:
        - The number of species in the overlap must be a 5th or less of the number of species in each clade - What if it's a single gene that's out of place? Won't make a difference then to the orthologs!
        - for each species with genes in both clades: the genes in one clade must all be more out of place (according to the 
          species tree) than all the gene from that species in the other tree
    """
    o_size = len(overlap)
    lsp0 = len(sp0)
    lsp1 = len(sp1)
    if (o_size == lsp0 or o_size == lsp1) or (rel_overlap_cutoff*o_size >= lsp0 and rel_overlap_cutoff*o_size >= lsp1): 
        return False, []
    # The overlap looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    n_a_removed = 0
    n_b_removed = 0
    q_resolved = True
    for sp in overlap:
        A = [g for g in ch[0].get_leaf_names() if GeneToSpecies(g) == sp]
        B = [g for g in ch[1].get_leaf_names() if GeneToSpecies(g) == sp]
        a_levels = []
        b_levels = []
        for X, level in zip((A,B),(a_levels, b_levels)):
            for g in X:
                if sp_sets is None:
                    gene_node = tree & g
                    species_below = lambda node: set([GeneToSpecies(gg) for gg in node.get_leaf_names()])
                else:
                    gene_node = sp_sets.leaf(g)
                    species_below = sp_sets.species
                r = gene_node.up
                next_species = species_below(r)
                # having a gene from the same species isn't enough?? No, but we add to the count I think.
                while len(next_species) == 1:
                    r = r.up
                    next_species = species_below(r)
                next_species.remove(sp)
                # get the level
                # the sum of the closest and furthest expected distance topological distance for the closest genes in the gene tree (based on species tree topology)
                neigh = neighbours[sp]
                observed = [neigh[n_sp] for n_sp in next_species]
                level.append(min(observed) + max(observed))
        q_remove_a = max(b_levels) + 2 < min(a_levels)   # if the clade is one step up the tree further way (min=max) then this gives +2. There's no way this is a problem                        
        q_remove_b = max(a_levels) + 2 < min(b_levels)                           
        if q_remove_a and rel_overlap_cutoff*o_size < len(sp0):
            n_a_removed += len(a_levels)
            genes_removed.extend(A)
        elif q_remove_b and rel_overlap_cutoff*o_size < len(sp1):
            n_b_removed += len(b_levels)
            genes_removed.extend(B)
        else:
            q_resolved = False
            break
    if q_resolved:
        return True, set(genes_removed)
    else:
        return False, set()
          
def orthologs_and_suspect(ch, suspect_genes, misplaced_genes, species_and_gene):
    """
    ch - the two child nodes that are orthologous
    suspect_genes - genes already identified as misplaced at lower levels
    misplaced_genes - genes identified as misplaced at this level

    Returns the tuple (o_0, o_1, os_0, os_1) where each element is a dictionary from species to genes from that species,
    the o are orthologs, the os are 'suspect' orthologs because the gene was previously identified as suspect
    """
    d = [defaultdict(list) for _ in range(2)]
    d_sus = [defaultdict(list) for _ in range(2)] 
    for node, di, d_susi in zip(ch, d, d_sus):
        for g in [g for g in node.get_leaf_names() if g not in misplaced_genes]:
            sp, seq = species_and_gene(g)
            if g in suspect_genes:
                d_susi[sp].append(seq)
            else:
                di[sp].append(seq)
    return d[0], d[1], d_sus[0], d_sus[1]


def get_lines_for_ortholog_files(
        orthologs_alltrees,
        species_dict,
        i_species_to_use,
        sequence_dict,
        q_contains_suspect_orthologs,
        ortholog_lines,
        ortholog_sus_lines,
        fewer_open_files,
        pair_counts=None,
    ):
    """
    Add ortholog rows to olog_lines/olog_sus_lines.

    Returns a dense NOrthoSp for these trees. If pair_counts is a list, the
    (iL, iR, nL, nR) records are appended to it instead and None is returned,
    which avoids allocating five nspecies x nspecies matrices per tree.
    """
    nsp = len(i_species_to_use)
    n_ortho = stats.NOrthoSp(nsp) if pair_counts is None else None

    # Rows are collected in lists and joined once at the end. Appending with
    # olog_lines[i][j] += row copies the whole accumulated string every time
    # (Python can't append in place to a string held in a list), which is
    # quadratic for large trees.
    ortholog_parts = defaultdict(list)   # (i, j) -> rows
    sus_parts = defaultdict(list)    # i -> rows

    sp_to_i = {str(sp): i for i, sp in enumerate(i_species_to_use)}
    sp_label = [species_dict[str(sp)] for sp in i_species_to_use]
    getrow = file_io.unquoted_line

    seq_cache = {}

    def seq_name(sp, g):
        key = str(sp) + "_" + str(g)

        v = seq_cache.get(key)

        if v is None:
            v = sequence_dict[key]
            seq_cache[key] = v

        return v

    def join_gene_names(sp, genes):
        return ", ".join(
            sorted(seq_name(sp, g) for g in genes)
        )

    def add_stats(i_l, i_r, n_l, n_r):
        if pair_counts is None:
            n_ortho.add_pair(i_l, i_r, n_l, n_r)
        else:
            pair_counts.append((i_l, i_r, n_l, n_r))

    def add_ortholog_row(og, sp_l, genes_l, sp_r, genes_r):
        """
        Add one species-pair ortholog row in both directions.
        """
        sp_l = str(sp_l)
        sp_r = str(sp_r)

        if sp_l == sp_r:
            return

        if not genes_l or not genes_r:
            return

        i_l = sp_to_i[sp_l]
        i_r = sp_to_i[sp_r]

        n_l = len(genes_l)
        n_r = len(genes_r)

        text_l = join_gene_names(sp_l, genes_l)
        text_r = join_gene_names(sp_r, genes_r)

        if fewer_open_files:
            ortholog_parts[(i_l, 0)].append(getrow((
                og,
                sp_label[i_r],
                text_l,
                text_r
            )))

            ortholog_parts[(i_r, 0)].append(getrow((
                og,
                sp_label[i_l],
                text_r,
                text_l
            )))

        else:
            ortholog_parts[(i_l, i_r)].append(getrow((
                og,
                text_l,
                text_r
            )))

            ortholog_parts[(i_r, i_l)].append(getrow((
                og,
                text_r,
                text_l
            )))

        add_stats(i_l, i_r, n_l, n_r)

    def add_suspect_row(og, sp_l, genes_l, sp_r, genes_r):
        sp_l = str(sp_l)
        sp_r = str(sp_r)

        if sp_l == sp_r:
            return

        if not genes_l or not genes_r:
            return

        i_l = sp_to_i[sp_l]
        i_r = sp_to_i[sp_r]

        text_l = join_gene_names(sp_l, genes_l)
        text_r = join_gene_names(sp_r, genes_r)

        sus_parts[i_l].append(getrow((og, text_l, text_r)))
        sus_parts[i_r].append(getrow((og, text_r, text_l)))

    for iog, orthologs_onetree in orthologs_alltrees:
        og = "OG%07d" % iog

        for leaves_l, leaves_r, leaves_l_sus, leaves_r_sus in orthologs_onetree:
            for sp_l in sorted(leaves_l.keys(), key=lambda sp: sp_to_i[str(sp)]):
                genes_l = leaves_l[sp_l]

                if not genes_l:
                    continue

                for sp_r in sorted(leaves_r.keys(), key=lambda sp: sp_to_i[str(sp)]):
                    genes_r = leaves_r[sp_r]

                    if not genes_r:
                        continue

                    add_ortholog_row(
                        og,
                        sp_l,
                        genes_l,
                        sp_r,
                        genes_r
                    )

            if not q_contains_suspect_orthologs:
                continue

            leaves_sus = (leaves_l_sus, leaves_r_sus)
            leaves_norm = (leaves_l, leaves_r)

            for i_pair in range(2):
                leaves_sus0 = leaves_sus[i_pair]
                leaves_sus1 = leaves_sus[1 - i_pair]
                leaves_norm1 = leaves_norm[1 - i_pair]

                partner_species = dict.fromkeys(
                    list(leaves_norm1.keys()) + list(leaves_sus1.keys())
                )

                for sp0 in sorted(leaves_sus0.keys(), key=lambda sp: sp_to_i[str(sp)]):
                    genes0 = leaves_sus0[sp0]

                    if not genes0:
                        continue

                    for sp1 in sorted(partner_species.keys(), key=lambda sp: sp_to_i[str(sp)]):
                        if str(sp0) == str(sp1):
                            continue

                        genes1 = (
                            list(leaves_norm1.get(sp1, [])) +
                            list(leaves_sus1.get(sp1, []))
                        )

                        if not genes1:
                            continue

                        add_suspect_row(
                            og,
                            sp0,
                            genes0,
                            sp1,
                            genes1
                        )

    for (i, j), parts in ortholog_parts.items():
        row = ortholog_lines[i]
        if isinstance(row, dict):
            # Sparse row: only species pairs with orthologs are stored.
            row[j] = row.get(j, "") + "".join(parts)
        else:
            row[j] += "".join(parts)
    for i, parts in sus_parts.items():
        ortholog_sus_lines[i] += "".join(parts)

    return n_ortho


def get_hogs_from_tree(
        iog,
        tree,
        hog_writer,
        q_split_paralogous_clades,
    ):
    """
    Generate HOG rows from one gene tree.

    Important:
        This function no longer writes to files.
        It returns cached_hogs to the parent process.
    """
    og_name = "OG%07d" % iog
    if debug:
        print("\n===== %s =====" % og_name)

    try:
        tree = hog_writer.mark_dups_below(tree)
        cached_hogs = []

        for n in tree.traverse("preorder"):
            cached_hogs.extend(
                hog_writer.write_clade_v2(
                    n,
                    og_name,
                    q_split_paralogous_clades
                )
            )

        return cached_hogs

    except Exception:
        print("WARNING: HOG analysis for %s failed" % og_name)
        print(
            "Please report to https://github.com/davidemms/OrthoFinder/issues "
            "including SpeciesTree_rooted_ids.txt and Trees_ids/%s_tree_id.txt "
            "from WorkingDirectory/" % og_name
        )
        raise
