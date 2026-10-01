# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import sys
import operator
import itertools
from collections import defaultdict 

import multiprocessing as mp

try:
    from rich import print
except ImportError:
    ...
from ..tools import tree as tree_lib
from . import resolve
from ..utils import util, files


debug = False   # HOGs


def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
  
OrthoFinderIDs = GeneToSpecies_dash


def SpeciesAndGene_dash(g):
  return g.split("_", 1)

  
SpeciesAndGene_lookup = {
    GeneToSpecies_dash:SpeciesAndGene_dash, 
}


class RootMap(object):
    def __init__(self, setA, setB, GeneToSpecies):
        self.setA = setA
        self.setB = setB
        self.GeneToSpecies = GeneToSpecies
        
    def GeneMap(self, gene_name):
        sp = self.GeneToSpecies(gene_name)
        if sp in self.setA: return True 
        elif sp in self.setB: return False
        else: return None

class TreeAnalyser(object):
    def __init__(
            self, 
            nOgs, 
            dResultsOrthologues, 
            reconTreesRenamedDir, 
            species_tree_rooted_labelled, 
            speciesToUse, 
            GeneToSpecies, 
            SequenceDict, 
            speciesDict, 
            spec_seq_dict, 
            neighbours, 
            qNoRecon, 
            dups_file_handle, 
            stride_dups, 
            ologs_files_handles, 
            putative_xenolog_file_handles, 
            hog_writer, 
            q_split_paralogous_clades, 
            fewer_open_files=False,
            exist_msa=True,
            write_hog_tree=True,
            fix_files=True,
    ):
        self.nOgs = nOgs
        self.dResultsOrthologues = dResultsOrthologues
        self.reconTreesRenamedDir = reconTreesRenamedDir
        self.species_tree_rooted_labelled = species_tree_rooted_labelled
        self.speciesToUse = speciesToUse
        self.nspecies = len(self.speciesToUse)
        self.GeneToSpecies = GeneToSpecies
        self.SequenceDict = SequenceDict
        self.speciesDict = speciesDict
        self.spec_seq_dict = spec_seq_dict
        self.spec_seq_dict[">%d genes" % (5*self.nspecies)] = ">%d genes" % (5*self.nspecies)   # used in Duplications
        self.neighbours = neighbours
        self.qNoRecon = qNoRecon
        self.dups_file_handle = dups_file_handle
        self.stride_dups = stride_dups
        self.ologs_files_handles = ologs_files_handles
        self.putative_xenolog_file_handles = putative_xenolog_file_handles
        self.hog_writer = hog_writer
        self.q_split_paralogous_clades = q_split_paralogous_clades
        self.fewer_open_files = fewer_open_files
        self.lock_ologs = [mp.Lock() for i in range(self.nspecies)]   # lock the larger of the two species index
        self.lock_dups = mp.Lock()
        self.lock_suspect = mp.Lock()
        self.lock_hogs = mp.Lock()
        self.exist_msa = exist_msa
        self.write_hog_tree = write_hog_tree
        self.fix_files = fix_files


    def AnalyseTree(self, iog, nOrtho_acc=None):
        """
        Analyse one gene tree.

        It returns a payload to the parent process.

        If nOrtho_acc (a util.nOrtho_sp) is given, this tree's orthologue
        counts are added to it once the whole tree has been analysed, and the
        payload's "n_orthologues" is None. Otherwise it holds a dense
        nOrtho_sp for this tree.
        """
        og_name = "OG%07d" % iog
        n_species = len(self.speciesToUse)

        try:
            if self.write_hog_tree or not self.fix_files:
                tree_file = files.FileHandler.GetOGsTreeFN(iog)
            else:
                tree_file = files.FileHandler.GetResolvedTreeIDDir() + "OG%07d.txt" % iog

            if not os.path.exists(tree_file):
                return None

            rooted_tree_ids, qHaveSupport = CheckAndRootTree(
                tree_file,
                self.species_tree_rooted_labelled,
                self.GeneToSpecies
            )

            if rooted_tree_ids is None:
                return None

            ologs, recon_tree, suspect_genes, dups = GetOrthologues_from_tree(
                iog,
                rooted_tree_ids,
                self.species_tree_rooted_labelled,
                self.GeneToSpecies,
                self.neighbours,
                q_get_dups=True,
                qNoRecon=self.qNoRecon
            )

            olog_lines = EmptyOlogLines(self.nspecies, self.fewer_open_files)
            olog_sus_lines = ["" for _ in range(self.nspecies)]

            pair_counts = [] if nOrtho_acc is not None else None
            nOrthologues_SpPair = GetLinesForOlogFiles(
                [(iog, ologs)],
                self.speciesDict,
                self.speciesToUse,
                self.SequenceDict,
                len(suspect_genes) > 0,
                olog_lines,
                olog_sus_lines,
                fewer_open_files=self.fewer_open_files,
                pair_counts=pair_counts,
            )

            cached_hogs = GetHOGs_from_tree(
                iog,
                recon_tree,
                self.hog_writer,
                self.q_split_paralogous_clades,
            )

            util.RenameTreeTaxa(
                recon_tree,
                self.reconTreesRenamedDir + "OG%07d.txt" % iog,
                self.spec_seq_dict,
                qSupport=False,
                qFixNegatives=True
            )

            if nOrtho_acc is not None:
                # Only count trees that were analysed completely.
                for counts in pair_counts:
                    nOrtho_acc.add_pair(*counts)

            return {
                "iog": iog,
                "n_orthologues": nOrthologues_SpPair,
                "olog_lines": olog_lines,
                "olog_sus_lines": olog_sus_lines,
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
                    util.nOrtho_sp(n_species) if nOrtho_acc is None else None
                ),
                "olog_lines": EmptyOlogLines(self.nspecies, self.fewer_open_files),
                "olog_sus_lines": ["" for _ in range(self.nspecies)],
                "duplications": [],
                "suspect_genes": set(),
                "cached_hogs": [],
            }

def EmptyOlogLines(nspecies, fewer_open_files):
    """
    Per-species containers for one tree's orthologue rows.

    In pairwise mode each row is a dict {j: text} holding only the partner
    species that have orthologues, instead of an nspecies x nspecies grid
    (1,000,000 cells per tree at 1,000 species). Use IterOlogRow to read it.
    """
    if fewer_open_files:
        return [[""] for _ in range(nspecies)]
    return [{} for _ in range(nspecies)]


def IterOlogRow(row):
    """(j, text) for the non-empty cells of an olog_lines row (list or dict)."""
    items = row.items() if isinstance(row, dict) else enumerate(row)
    return ((j, text) for j, text in items if text)


def CheckAndRootTree(treeFN, species_tree_rooted, GeneToSpecies):
    """
    Check that the tree can be analysed and rooted
    Root tree
    Returns None if this fails, i.e. checks: exists, has more than one gene, can be rooted
    """
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0:
        return None, False
    qHaveSupport = False
    try:
        tree = tree_lib.Tree(treeFN, format=2)
        qHaveSupport = True
    except:
        try:
            tree = tree_lib.Tree(treeFN)
        except:
            tree = tree_lib.Tree(treeFN, format=3)
    if len(tree) == 1: 
        return None, False
    root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    if root == None: 
        return None, False
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)
    return tree, qHaveSupport

def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = GetRoots(tree, species_tree_rooted, GeneToSpecies)
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
        else:
            return None # single species tree
        
def GetRoots(tree, species_tree_rooted, GeneToSpecies):
    """
    Allow non-binary gene or species trees.
    (A,B,C) => consider splits A|BC, B|AC, C|AB - this applies to gene and species tree
    If a clean ingroup/outgroup split cannot be found then score root by geometric mean of fraction of expected species actually 
    observed for the two splits
    """
    speciesObserved = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
    if len(speciesObserved) == 1:
        return [next(n for n in tree)] # arbitrary root if all genes are from the same species
    
    # use species tree to find correct outgroup according to what species are present in the gene tree
    n = species_tree_rooted
    children = n.get_children()
    leaves = [set(ch.get_leaf_names()) for ch in children]
    have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
    while sum(have) < 2:
        n = children[have.index(True)]
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]

    # Get splits to look for
    roots_list = []
    scores_list = []   # the fraction completeness of the two clades
#    roots_set = set()
    # Species below/above each node: the same for every split, so computed once
    # (the tree is not modified in this loop).
    StoreSpeciesSets(tree, GeneToSpecies)   # sets of species
    for i in range(len(leaves)):
        t1 = leaves[i]
        t2 = set.union(*[l for j,l in enumerate(leaves) if j!=i])
        # G - set of species in gene tree
        # First relevant split in species tree is (A,B), such that A \cap G \neq \emptyset and A \cap G \neq \emptyset
        # label all nodes in gene tree according the whether subsets of A, B or both lie below node
        root_mapper = RootMap(t1, t2, GeneToSpecies)    
        sett1 = set(t1)
        sett2 = set(t2)
        nt1 = float(len(t1))
        nt2 = float(len(t2))
        N_recip = 1./(nt1*nt1*nt2*nt2)
        GeneMap = root_mapper.GeneMap
        StoreSpeciesSets(tree, GeneMap, "inout_") # ingroup/outgroup identification
        # find all possible locations in the gene tree at which the root should be

        T = {True,}
        F = {False,}
        TF = set([True, False])
        for m in tree.traverse('postorder'):
            if m.is_leaf(): 
                if len(m.inout_up) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
            else:
                if len(m.inout_up) == 1 and len(m.inout_down) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
                nodes = m.get_children() if m.is_root() else [m] + m.get_children()
                clades = [ch.inout_down for ch in nodes] if m.is_root() else ([m.inout_up] + [ch.inout_down for ch in m.get_children()])
                # do we have the situation A | B or (A,B),S?
                if len(nodes) == 3:
                    if all([len(c) == 1 for c in clades]) and T in clades and F in clades:
                        # unique root
                        if clades.count(T) == 1:
                            return [nodes[clades.index(T)]]
                        else:
                            return [nodes[clades.index(F)]]
                    elif T in clades and F in clades:
                        #AB-(A,B) or B-(AB,A)
                        ab = [c == TF for c in clades]
                        i = ab.index(True)
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                    elif clades.count(TF) >= 2:  
                        # (A,A,A)-excluded, (A,A,AB)-ignore as want A to be bigest without including B, (A,AB,AB), (AB,AB,AB) 
                        i = 0
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                elif T in clades and F in clades:
                    roots_list.append(m)
                    scores_list.append(0)  # last choice
    # If we haven't found a unique root then use the scores for completeness of ingroup/outgroup to root
    if len(roots_list) == 0: 
        return [] # This shouldn't occur
    return [sorted(zip(scores_list, roots_list), key=lambda x: x[0], reverse=True)[0][1]]


def StoreSpeciesSets(t, GeneMap, tag="sp_"):
    if t is None:
        raise ValueError("StoreSpeciesSets got t=None (tree is None)")
    tag_up = tag + "up"
    tag_down = tag + "down"  
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node_type = GeneMap(node.name)
            node.add_feature(tag_down, set() if node_type is None else {node_type})
        elif node.is_root():
            continue
        else:
            node.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in node.get_children()]))
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature(tag_up, set())
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                node.add_feature(tag_up, set.union(*[other.__getattribute__(tag_down) for other in others]))
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                sp_downs = set.union(*[other.__getattribute__(tag_down) for other in others])
                node.add_feature(tag_up, parent.__getattribute__(tag_up).union(sp_downs))
    t.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in t.get_children()]))

def OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, n1, n2):
    f_dup = len(sp_up.intersection(sett1)) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * len(sp_down.intersection(sett2)) * N_recip
    f_a = len(sp_up.intersection(sett1)) * (n2-len(sp_up.intersection(sett2))) * (n1-len(sp_down.intersection(sett1))) * len(sp_down.intersection(sett2)) * N_recip
    f_b = (n1-len(sp_up.intersection(sett1))) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * (n2-len(sp_down.intersection(sett2))) * N_recip
    choice = (f_dup, f_a, f_b)
    return max(choice)



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


def GetOrthologues_from_tree(
        iog,
        tree,
        species_tree_rooted,
        GeneToSpecies,
        neighbours,
        q_get_dups=False,
        qNoRecon=False,
    ):
    """Infer orthologues and duplications from a reconciled gene tree."""
    og_name = "OG%07d" % iog
    qPrune = False

    SpeciesAndGene = SpeciesAndGene_lookup[GeneToSpecies]

    orthologues = []
    duplications = []

    if not qNoRecon:
        tree = Resolve(tree, GeneToSpecies)

    if qPrune:
        tree.prune(tree.get_leaf_names())

    if len(tree) == 1:
        return set(orthologues), tree, set(), duplications

    iNode = 1
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
            node = MRCA_node(species_tree_rooted, species)
            mrca_cache[mask] = node
        return node

    for n in tree.traverse("preorder"):
        if n.is_leaf():
            continue

        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1

        ch = n.get_children()

        if len(ch) == 2:
            oSize, overlap, sp0, sp1 = OverlapSize(
                n,
                GeneToSpecies,
                suspect_genes,
                sp_sets=sp_sets,
            )

            sp_present = sp0.union(sp1)
            if len(sp_present) == len(sp_sets.species(n)):
                # No species removed as suspect: the set is exactly the species below n.
                stNode = species_tree_mrca(sp_sets.mask(n), sp_present)
            else:
                stNode = MRCA_node(species_tree_rooted, sp_present)

            n.add_feature("sp_node", stNode.name)

            if oSize != 0:
                qResolved, misplaced_genes = ResolveOverlap(
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

            dup = oSize != 0 and not qResolved

            n.add_feature("dup", dup)

            if dup:
                n.add_feature("dup_type", "species_overlap_duplication")

                if q_get_dups:
                    genes0 = ch[0].get_leaf_names()
                    genes1 = ch[1].get_leaf_names()

                    duplications.append((
                        stNode.name,
                        n.name,
                        float(oSize) / len(stNode),
                        genes0,
                        genes1
                    ))

            else:
                n.add_feature("dup_type", "speciation")

                orthologues.append(
                    Orthologs_and_Suspect(
                        ch,
                        suspect_genes,
                        misplaced_genes,
                        SpeciesAndGene
                    )
                )

                suspect_genes.update(misplaced_genes)

        elif len(ch) > 2:
            species = [sp_sets.species(child) for child in ch]
            all_species = set.union(*species)
            stNode = species_tree_mrca(sp_sets.mask(n), all_species)
            n.add_feature("sp_node", stNode.name)

            if len(all_species) == 1:
                genes = n.get_leaf_names()
                duplications.append((stNode.name, n.name, 1.0, genes, []))
                n.add_feature("dup", True)
            else:
                dups = []
                for (n0, s0), (n1, s1) in itertools.combinations(
                        zip(ch, species), 2):
                    if len(s0.intersection(s1)) == 0:
                        orthologues.append(
                            Orthologs_and_Suspect(
                                (n0, n1),
                                suspect_genes,
                                empty_set,
                                SpeciesAndGene
                            )
                        )
                        dups.append(False)
                    else:
                        dups.append(True)

                if all(dups):
                    genes = n.get_leaf_names()
                    duplications.append((stNode.name, n.name, 1.0, genes, []))

                n.add_feature("dup", all(dups))

    return orthologues, tree, suspect_genes, duplications


def Resolve(tree, GeneToSpecies):
    if tree is None:
        raise ValueError("Resolve got tree=None")
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        new_tree = resolve.resolve(n, GeneToSpecies)
        if new_tree is not None:
            tree = new_tree
    return tree


def MRCA_node(t_rooted, taxa):
    return (t_rooted & next(taxon for taxon in taxa)) if len(taxa) == 1 else t_rooted.get_common_ancestor(taxa)


def OverlapSize(node, GeneToSpecies, suspect_genes, sp_sets=None):
    if sp_sets is None:
        descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()}.difference(suspect_genes) for n in node.get_children()]
    else:
        descendents = [sp_sets.species(n).difference(suspect_genes) for n in node.get_children()]
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]

def ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies, relOverlapCutoff=4, sp_sets=None):
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
    oSize = len(overlap)
    lsp0 = len(sp0)
    lsp1 = len(sp1)
    if (oSize == lsp0 or oSize == lsp1) or (relOverlapCutoff*oSize >= lsp0 and relOverlapCutoff*oSize >= lsp1): 
        return False, []
    # The overlap looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    nA_removed = 0
    nB_removed = 0
    qResolved = True
    for sp in overlap:
        A = [g for g in ch[0].get_leaf_names() if GeneToSpecies(g) == sp]
        B = [g for g in ch[1].get_leaf_names() if GeneToSpecies(g) == sp]
        A_levels = []
        B_levels = []
        for X, level in zip((A,B),(A_levels, B_levels)):
            for g in X:
                if sp_sets is None:
                    gene_node = tree & g
                    species_below = lambda node: set([GeneToSpecies(gg) for gg in node.get_leaf_names()])
                else:
                    gene_node = sp_sets.leaf(g)
                    species_below = sp_sets.species
                r = gene_node.up
                nextSpecies = species_below(r)
                # having a gene from the same species isn't enough?? No, but we add to the count I think.
                while len(nextSpecies) == 1:
                    r = r.up
                    nextSpecies = species_below(r)
                nextSpecies.remove(sp)
                # get the level
                # the sum of the closest and furthest expected distance topological distance for the closest genes in the gene tree (based on species tree topology)
                neigh = neighbours[sp]
                observed = [neigh[nSp] for nSp in nextSpecies]
                level.append(min(observed) + max(observed))
        qRemoveA = max(B_levels) + 2 < min(A_levels)   # if the clade is one step up the tree further way (min=max) then this gives +2. There's no way this is a problem                        
        qRemoveB = max(A_levels) + 2 < min(B_levels)                           
        if qRemoveA and relOverlapCutoff*oSize < len(sp0):
            nA_removed += len(A_levels)
            genes_removed.extend(A)
        elif qRemoveB and relOverlapCutoff*oSize < len(sp1):
            nB_removed += len(B_levels)
            genes_removed.extend(B)
        else:
            qResolved = False
            break
    if qResolved:
        return True, set(genes_removed)
    else:
        return False, set()
          
def Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene):
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
            sp, seq = SpeciesAndGene(g)
            if g in suspect_genes:
                d_susi[sp].append(seq)
            else:
                di[sp].append(seq)
    return d[0], d[1], d_sus[0], d_sus[1]


def GetLinesForOlogFiles(
        orthologues_alltrees,
        speciesDict,
        iSpeciesToUse,
        sequenceDict,
        qContainsSuspectOlogs,
        olog_lines,
        olog_sus_lines,
        fewer_open_files,
        pair_counts=None,
    ):
    """
    Add orthologue rows to olog_lines/olog_sus_lines.

    Returns a dense nOrtho_sp for these trees. If pair_counts is a list, the
    (iL, iR, nL, nR) records are appended to it instead and None is returned,
    which avoids allocating five nspecies x nspecies matrices per tree.
    """
    nsp = len(iSpeciesToUse)
    nOrtho = util.nOrtho_sp(nsp) if pair_counts is None else None

    # Rows are collected in lists and joined once at the end. Appending with
    # olog_lines[i][j] += row copies the whole accumulated string every time
    # (Python can't append in place to a string held in a list), which is
    # quadratic for large trees.
    olog_parts = defaultdict(list)   # (i, j) -> rows
    sus_parts = defaultdict(list)    # i -> rows

    sp_to_i = {str(sp): i for i, sp in enumerate(iSpeciesToUse)}
    sp_label = [speciesDict[str(sp)] for sp in iSpeciesToUse]
    getrow = util.getrow

    seq_cache = {}

    def seq_name(sp, g):
        key = str(sp) + "_" + str(g)

        v = seq_cache.get(key)

        if v is None:
            v = sequenceDict[key]
            seq_cache[key] = v

        return v

    def join_gene_names(sp, genes):
        return ", ".join(
            sorted(seq_name(sp, g) for g in genes)
        )

    def add_stats(iL, iR, nL, nR):
        if pair_counts is None:
            nOrtho.add_pair(iL, iR, nL, nR)
        else:
            pair_counts.append((iL, iR, nL, nR))

    def add_olog_row(og, spL, genesL, spR, genesR):
        """
        Add one species-pair orthologue row in both directions.
        """
        spL = str(spL)
        spR = str(spR)

        if spL == spR:
            return

        if not genesL or not genesR:
            return

        iL = sp_to_i[spL]
        iR = sp_to_i[spR]

        nL = len(genesL)
        nR = len(genesR)

        textL = join_gene_names(spL, genesL)
        textR = join_gene_names(spR, genesR)

        if fewer_open_files:
            olog_parts[(iL, 0)].append(getrow((
                og,
                sp_label[iR],
                textL,
                textR
            )))

            olog_parts[(iR, 0)].append(getrow((
                og,
                sp_label[iL],
                textR,
                textL
            )))

        else:
            olog_parts[(iL, iR)].append(getrow((
                og,
                textL,
                textR
            )))

            olog_parts[(iR, iL)].append(getrow((
                og,
                textR,
                textL
            )))

        add_stats(iL, iR, nL, nR)

    def add_suspect_row(og, spL, genesL, spR, genesR):
        spL = str(spL)
        spR = str(spR)

        if spL == spR:
            return

        if not genesL or not genesR:
            return

        iL = sp_to_i[spL]
        iR = sp_to_i[spR]

        textL = join_gene_names(spL, genesL)
        textR = join_gene_names(spR, genesR)

        sus_parts[iL].append(getrow((og, textL, textR)))
        sus_parts[iR].append(getrow((og, textR, textL)))

    for iog, orthologues_onetree in orthologues_alltrees:
        og = "OG%07d" % iog

        for leavesL, leavesR, leavesL_sus, leavesR_sus in orthologues_onetree:
            for spL in sorted(leavesL.keys(), key=lambda sp: sp_to_i[str(sp)]):
                genesL = leavesL[spL]

                if not genesL:
                    continue

                for spR in sorted(leavesR.keys(), key=lambda sp: sp_to_i[str(sp)]):
                    genesR = leavesR[spR]

                    if not genesR:
                        continue

                    add_olog_row(
                        og,
                        spL,
                        genesL,
                        spR,
                        genesR
                    )

            if not qContainsSuspectOlogs:
                continue

            leaves_sus = (leavesL_sus, leavesR_sus)
            leaves_norm = (leavesL, leavesR)

            for iPair in range(2):
                leaves_sus0 = leaves_sus[iPair]
                leaves_sus1 = leaves_sus[1 - iPair]
                leaves_norm1 = leaves_norm[1 - iPair]

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

    for (i, j), parts in olog_parts.items():
        row = olog_lines[i]
        if isinstance(row, dict):
            # Sparse row: only species pairs with orthologues are stored.
            row[j] = row.get(j, "") + "".join(parts)
        else:
            row[j] += "".join(parts)
    for i, parts in sus_parts.items():
        olog_sus_lines[i] += "".join(parts)

    return nOrtho


def GetHOGs_from_tree(
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
