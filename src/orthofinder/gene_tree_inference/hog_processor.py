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
import itertools
from collections import defaultdict

import multiprocessing as mp

try:
    from rich import print
except ImportError:
    ...

from .trees2ologs_processor import LazyFileCache
from ..utils import util, files


PY2 = sys.version_info <= (3,)

debug = False   # HOGs

class HogWriter(object):
    def __init__(
            self, 
            species_tree, 
            species_tree_node_names, 
            seq_ids, 
            sp_ids, 
            species_to_use,
            write_to_rd,
            write_output=True
        ):
        """
        Prepare files, get ready to write.
        species_tree_node_names - list of species tree nodes
        seq_ids - dict of sequence ids
        sp_ids - dict of species ids
        species_to_use - list of ints
        """
        self.sp_ids = sp_ids
        self.seq_ids = seq_ids
        if write_to_rd:
            q_results = True
        else:
            q_results = False

        self.fhs = dict()
        self.iSps = list(map(str, sorted(species_to_use)))   # list of strings
        self.i_sp_to_index = {int(isp):i_col for i_col, isp in enumerate(self.iSps)}
        self.iHOG = defaultdict(int)
        self.lock_iHOG = mp.Lock()
        self.species_tree = species_tree
        species_names = [sp_ids[i] for i in self.iSps]

        self.write_output = write_output
        self.q_results = q_results
        self.node_names = list(species_tree_node_names)
        self.species_names = species_names
        self.hog_paths = {}
        self.file_cache = LazyFileCache(max_open=64)

        if write_output:
            d = os.path.dirname(
                files.FileHandler.GetHierarchicalOrthogroupsFN(
                    "N0",
                    q_results=q_results
                )
            )
            if not os.path.exists(d):
                os.mkdir(d)

            for name in species_tree_node_names + ["N0.ids"]:
                if not name.endswith(".ids"):
                    fn = files.FileHandler.GetHierarchicalOrthogroupsFN(
                        name,
                        q_results=q_results
                    )
                else:
                    fn = files.FileHandler.GetHierarchicalOrthogroupsFN(
                        name,
                        q_results=q_results,
                        extension=".ids"
                    )

                self.hog_paths[name] = fn
                
                with open(fn, util.csv_write_mode) as fh:
                    util.writerow(
                        fh,
                        ["HOG", "OG", "Gene Tree Parent Clade"] + species_names
                    )

        # Map from HOGs to genes that must be contained in them
        self.hog_contents = dict()  # sp_node_name = hog_name-> list of contents fo hog (internal nodes and leaves)
        for n in species_tree.traverse():
            desc = n.get_descendants()
            self.hog_contents[n.name] = set([int(nn.name) if nn.is_leaf() else nn.name for nn in desc])
        self.comp_nodes = self.get_comparable_nodes(self.species_tree)


    def _get_hog_handle(self, hog_name):
        if not self.write_output:
            return None
        fn = self.hog_paths[hog_name]
        return self.file_cache.get(fn, util.csv_append_mode, gz=False)


    def close_files(self):
        if hasattr(self, "file_cache"):
            self.file_cache.close_all()

    def get_hog_index(self, hog_name):
        with self.lock_iHOG:
            i = self.iHOG[hog_name]
            self.iHOG[hog_name] += 1
            return i
        
    def write_hog_genes(self, genes, sp_node_name_list, og_name):
        """
        Write two & three gene HOGs.

        This is now lazy: it opens only the target HOG file currently being written.
        """
        if not self.write_output:
            return

        if len(sp_node_name_list) == 0:
            return

        genes_per_species = defaultdict(list)
        genes_per_species_ids = defaultdict(list)

        for g in genes:
            isp, _ = g.split("_")
            genes_per_species[isp].append(self.seq_ids[g])
            genes_per_species_ids[isp].append(g)

        i_hogs = [
            self.get_hog_index(sp_node_name)
            for sp_node_name in sp_node_name_list
        ]

        row_genes = [
            ", ".join(genes_per_species[isp])
            for isp in self.iSps
        ]

        for i_hog, sp_node_name in zip(i_hogs, sp_node_name_list):
            hog_id = "%s.HOG%07d" % (sp_node_name, i_hog)
            fh = self._get_hog_handle(sp_node_name)
            util.writerow(
                fh,
                [hog_id, og_name, "-"] + row_genes
            )

            if sp_node_name == "N0":
                row_genes_ids = [
                    ", ".join(genes_per_species_ids[isp])
                    for isp in self.iSps
                ]
                fh_ids = self._get_hog_handle("N0.ids")
                util.writerow(
                    fh_ids,
                    [hog_id, og_name, "-"] + row_genes_ids
                )


    def write_clade_v2(self, n, og_name, split_paralogous_clades_from_same_hog=False):
        """
        Original HOG behaviour.

        Do not use polytomy_boundary here.
        HOG writing should depend on the original dup/sp_node/dups_below logic.
        """
        if n.is_leaf():
            return []

        if debug:
            print("\nTree node: %s" % n.name)
            print(n.sp_node)

        if getattr(n, "dup", False) and n.sp_node == "N0":
            n.add_feature("done", set())

        ch = n.get_children()

        if (
                split_paralogous_clades_from_same_hog
                and getattr(n, "dup", False)
                and len(ch) == 2
                and ch[0].sp_node == ch[1].sp_node
            ):
            hogs_to_write = (
                set()
                if n.sp_node.startswith("N")
                else self.comp_nodes[n.sp_node][0].copy()
            )
        else:
            hogs_to_write = self.comp_nodes[n.sp_node][0].copy()

        genes_ids_per_species_id = self.get_descendant_genes(n)

        if debug:
            print("Dups below: " + str(n.dups_below))

        stop_at_dups = lambda nn: nn.name in n.dups_below

        sp_node = self.species_tree & n.sp_node

        hogs_to_write.update({
            nn.name
            for nn in sp_node.traverse("preorder", is_leaf_fn=stop_at_dups)
            if (not nn.is_leaf()) and (nn.name not in n.dups_below)
        })

        if not n.is_root():
            hogs_to_write.difference_update(n.up.done)

        n.add_feature(
            "done",
            hogs_to_write if n.is_root() else n.up.done.union(hogs_to_write)
        )

        if len(hogs_to_write) == 0:
            return []

        if debug:
            print(hogs_to_write)

        return self.get_hog_file_entries(
            hogs_to_write,
            genes_ids_per_species_id,
            og_name,
            n.name
        )


    
    # def write_clade_v2(self, n, og_name, split_paralogous_clades_from_same_hog=False):
    #     """
    #     Look at parent node to know when to start, look at dups below to know when 
    #     to stop.
    #     - Current MRCA could be excluded either because it's already been done or
    #       because of a duplication below

    #     - Species-specific clades could still be HOGs if they are all that remain
    #       of that clade 
    #       Args:
    #         n - gene tree node
    #         og_name - name to use in file output
    #         split_paralogous_clades_from_same_hog - should clades which are within 
    #             the same HOG but are paralogous be split up
    #     """
    #     if n.is_leaf():
    #         return []
    #     if debug: print("\nTree node: %s" % n.name)
    #     if debug: print(n.sp_node)
    #     if (n.dup and n.sp_node == "N0"): 
    #         n.add_feature("done", set())
    #     # self.comp_nodes[n.sp_node] is the set of HOGs relevant to this node

    #     # Only skip doing HOGs for above if it is a dup and want to split paralogous clades from same HOG
    #     ch = n.get_children()
    #     if (split_paralogous_clades_from_same_hog and n.dup and (ch[0].sp_node == ch[1].sp_node)):
    #         # continue to record single-species orthogroups
    #         hogs_to_write = set() if n.sp_node.startswith("N") else self.comp_nodes[n.sp_node][0].copy()
    #     else:
    #         hogs_to_write = self.comp_nodes[n.sp_node][0].copy()
        
    #     # get scl & remove HOGs that can't be written yet due to duplications
    #     # 0. Get the scl units below this node in the gene tree
    #     # I.e. get genes (referenced by species ID) below each scl (the relevant gene tree nodes)
    #     genes_ids_per_species_id = self.get_descendant_genes(n)

    #     # scl_mrca = {nn.sp_node for nn in scl if not nn.is_leaf()}
    #     if debug: print("Dups below: " + str(n.dups_below))
    #     stop_at_dups = lambda nn : nn.name in n.dups_below
    #     sp_node = self.species_tree & (n.sp_node)
    #     # don't need skip for dups, that's recorded in dups_below
    #     # traverse the species tree from the current node and record all nodes before hitting a duplication node from the gene tree
    #     hogs_to_write.update({nn.name for nn in sp_node.traverse('preorder', is_leaf_fn = stop_at_dups) if (not nn.is_leaf()) and (not nn.name in n.dups_below)})
        
    #     if not n.is_root():
    #         hogs_to_write.difference_update(n.up.done)

    #     # 2. Write HOGs
    #     n.add_feature("done", hogs_to_write if n.is_root() else n.up.done.union(hogs_to_write))
    #     if len(hogs_to_write) == 0:
    #         return []

    #     if debug: print(hogs_to_write)
    #     return self.get_hog_file_entries(hogs_to_write, genes_ids_per_species_id, og_name, n.name)




    def get_descendant_genes(self, n):
        """
        Attempt at a simplified replacement to get_scl_units as shouldn't need to 
        care about which scl a gene belongs to.
        Args:
            n - node under consideration
        Returns:
            dict:sp_id (int)->string of genes ids, comma separated
        """
        genes_per_species = defaultdict(list) # iCol (before 'name' columns) -> text string of genes
        genes = n.get_leaves()
        q_have_legitimate_gene = False # may be misplaced genes
        for g in genes:
            if "X" in g.features: continue
            q_have_legitimate_gene = True
            isp = int(g.name.split("_")[0])
            genes_per_species[isp].append(g.name)
        for k, v in genes_per_species.items():
            genes_per_species[k] = ", ".join(v)
        return genes_per_species

    def get_hog_file_entries(self, hogs_to_write, genes_ids_per_species_id, og_name, gt_node_name):
        """
        Write the HOGs that can be determined from this gene tree node.
        Args:
            hogs_to_write - list of HOG names
            genes_ids_per_species_id - dict:sp_id (int)->str, comma separated list of gene ids
            og_name - OG name
            gt_node_name - gene tree node name
        Implementation:
            - We have the HOGs that need writing plus knowledge of what scl units 
              each HOG should contain. For each hog take the intersection of what 
              we have with what the hog should contain.
        """
        ret = []
        for h in hogs_to_write:
            # print("HOG: " + h)
            q_empty = True
            # 2. We know the scl, these are the 'taxonomic units' available (clades or individual species in species tree for this node of the gene tree)
            # Note there can be at most one of each. Only a subset of these will fall under this HOG.
            units = self.hog_contents[h].intersection(genes_ids_per_species_id.keys())
            # print("Units: " + str(units))
            genes_row_ids = ["" for _ in self.iSps]
            genes_row = ["" for _ in self.iSps]
            # put the units into the row
            for isp in units:
                # translate the species ID to the species column it should be in
                # after accounting for removed species
                genes_row_ids[self.i_sp_to_index[isp]] = genes_ids_per_species_id[isp]
                genes_row[self.i_sp_to_index[isp]] = ", ".join(sorted([self.seq_ids[g] for g in genes_ids_per_species_id[isp].split(", ")]))
                q_empty = False
            if not q_empty: 
                # print((h, genes_row))
                ret.append((h, [og_name, gt_node_name] + genes_row))
                if h == "N0":
                    ret.append(("N0.ids", [og_name, gt_node_name] + genes_row_ids))
                # self.writers[h].writerow(["%s.HOG%07d" % (h, self.get_hog_index(h)),  og_name, gt_node_name] + genes_row)
        return ret

    # def close_files(self):
    #     for fh in self.fhs.values():
    #         fh.close()

    @staticmethod
    def get_skipped_nodes(n_sp_this, n_above_name, n_stop = None, n_gene=None):
        """
        Get the HOGs for the series of skipped species tree nodes above the current 
        node. 
        Args:
            n_sp_this - ete3 node from species tree 
            n_above_name - MRCA species tree node name for the gene tree node above
            n_stop - a HOG name above the MRCA that should be the last one added
            n_gene - ete3 node from the gene tree
        Implementation/Questions:
            - This is only used by the OGs with fewer than 4 taxa (and therefore 
              no tree) now
        """
        n = n_sp_this
        missed_sp_node_names = []
        if n.name == n_above_name:
            return missed_sp_node_names
        n = n.up
        while n is not None and n.name != n_above_name:
            missed_sp_node_names.append(n.name)
            if n.name == n_stop:
                break
            n = n.up
        # if above node is a duplication then we won't have written out a HOG for that, pass the next node 
        if n_stop is None:
            if n_gene is not None and n_gene.up is not None and n_gene.up.dup and n is not None:
                # then we also need to write the HOG above
                missed_sp_node_names.append(n.name)
        return missed_sp_node_names



    def mark_dups_below(self, tree):
        """
        Original HOG duplication logic, with one safety fix:
        every internal node always gets dups_below.
        """
        for n in tree.traverse("postorder"):
            if n.is_leaf():
                n.sp_node = n.name.split("_")[0]
                continue

            q_dup_evidenced = False

            if getattr(n, "dup", False):
                mrcas = [
                    ch.name.split("_")[0] if ch.is_leaf() else ch.sp_node
                    for ch in n.get_children()
                ]

                if len(set(mrcas)) == 1 and len(mrcas) > 1:
                    n.add_feature("dup_level", mrcas[0])
                    q_dup_evidenced = True
                else:
                    l = self.get_evidenced_dup_level(mrcas)

                    if l is None:
                        # Do not crash and do not continue.
                        # This node is not a valid evidenced duplication.
                        n.dup = False
                    else:
                        n.add_feature("dup_level", l)
                        q_dup_evidenced = True

            dups_below = set()

            for ch in n.get_children():
                if ch.is_leaf():
                    continue
                dups_below.update(getattr(ch, "dups_below", set()))

            if (
                    getattr(n, "dup", False)
                    and q_dup_evidenced
                    and getattr(n, "dup_level", "").startswith("N")
                ):
                dups_below.add(n.dup_level)

            n.add_feature("dups_below", dups_below)

        return tree


    def get_evidenced_dup_level(self, mrcas):
        """
        Implementation
        - V3.1: Currently, we need a representative from X and Y in both ch1 & ch2.
        - What if we asked for evidence from each descendant clade of evidence of 
          a duplication?
            - Version that would be too stringent for duplication identification: 
              genetree =(ch1, ch2)n, sptree = (X,Y) and ask for a
              single representative of X that is in both ch1 & ch2 and similarly 
              for Y in ch1 & ch2.
            - Better version: V3.1 criterion (X & Y seen in ch1 & 2) plus two copies 
              of a gene from X and two copies of a gene from Y in clade n (don't 
              have to be correct topology such that they fall correctly in ch1 & 
              ch2, just evidence of duplicated genes).
                - Note, these two copies could arise from a separate & well-evidenced
                  lower duplication. The idea was that this would strike the right
                  balance, but can we do better? We'd have to identify the genes
                  below that go into each component of the interpretation of the 
                  tree. This is for another time, too much for now.

        - a. get MRCA for ch1 & ch2 and then the lower of the two
        - b. get all multi-copy species, get MRCA
        - Get lower of a & b
        """
        # try each in turn and see if it is supported
        attested = set()
        for l1, l2 in itertools.combinations(mrcas, 2):
            if l1 == l2:
                attested.add(l1)
            elif l2 in self.comp_nodes[l1][1]:
                attested.add(l2)
            elif l1 in self.comp_nodes[l2][1]:
                attested.add(l1)
        if len(attested) == 1:
            return attested.pop()
        elif len(attested) == 0:
            print("WARNING: Unexpected gene tree topology")
            print(mrcas)
            # raise Exception()
            return None
        else:
            # get the highest in the tree
            attested = list(attested)
            ancestor_lists = [(self.species_tree & a).get_ancestors() for a in attested]
            x = len(attested)
            for i in range(x):
                if all(attested[i] in ancestor_lists[j] for j in range(x) if j!=i):
                    return attested[i]
        print("WARNING: Unexpected gene tree topology 2")
        print(mrcas)
        # raise Exception()
        return None

    def WriteCachedHOGs(
            self,
            cached_hogs,
            lock_hogs=None,
        ):

        if not self.write_output:
            return

        if not cached_hogs:
            return

        d = defaultdict(list)
        for h, row in cached_hogs:
            d[h].append(row)

        if lock_hogs is not None:
            lock_hogs.acquire()

        try:
            for h, hog_rows in d.items():
                fh = self._get_hog_handle(h)
                for r in hog_rows:
                    hog_id = "%s.HOG%07d" % (h, self.get_hog_index(h))
                    util.writerow(fh, [hog_id] + r)
        finally:
            if lock_hogs is not None:
                lock_hogs.release()

    @staticmethod
    def is_hog_boundary(n):
        return bool(
            getattr(n, "dup", False) or
            getattr(n, "polytomy_boundary", False)
        )

    @staticmethod
    def scl_fn(n):
        return n.is_leaf() or getattr(n, "dup", False)

    def get_comparable_nodes(self, sp_tree):
        """
        Return a dictionary of comaprable nodes
        Node NX < NY if NX is on the path between NY and the root.
        If a node is not <, =, > another then they are incomparable
        Args:
            sp_tree - sp_tree with labelled nodes
        Returns:
            comp_nodes - dict:NX -> ( {n|n<NX}, {n|n>NX} ) i.e. (higher_nodes, lower_nodes)
        """
        comp_nodes = dict()
        for n in sp_tree.traverse('postorder'):
            nodes_below = set()
            if not n.is_leaf():
                for ch in n.get_children():
                    if not ch.is_leaf():
                        nodes_below.update(ch.nodes_below)
                    nodes_below.add(ch.name)
            above = set([nn.name for nn in n.get_ancestors()])
            n.add_feature('nodes_below', nodes_below)
            comp_nodes[n.name] = (above, nodes_below, above.union(nodes_below.union(set(n.name))))
        return comp_nodes


