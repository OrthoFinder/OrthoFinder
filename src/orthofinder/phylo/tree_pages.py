"""
--tree-view: a page to look at the species tree and the
resolved gene trees in a web browser, with each gene tree's alignment:

    Resolved_Gene_Trees/Resolved_Gene_Trees.html
    Resolved_Gene_Trees/Tree_Data/<orthogroup>.js

The page draws the trees itself (in plain JavaScript, nothing to install, no
internet needed): a click on a node collapses or expands its clade, and the
tree can be saved as an SVG picture as it is shown. Below each gene tree is its
alignment (MultipleSequenceAlignments/<orthogroup>.fa), its rows in the order
of the tree's genes. The data of each gene tree is written as a script the page
loads when the tree is shown: a page opened from disk may load scripts, but
not read other files.

Genes are labelled "<species>: <gene>" and coloured by species, the same
colour everywhere; a branch is in a species' colour when all the genes below it
are of that species. Gene duplications (Gene_Duplication_Events) are marked on
their nodes. The internal nodes are labelled as in the results: N0, N1, ... in
the species tree (the HOG files), n0, n1, ... in the gene trees (their "Gene
tree parent clade").

The species tree explorer (a second tab) opens the species tree node by node,
in three ways: the HOGs at each node (Phylogenetic_Hierarchical_Orthogroups),
each with its clade of the gene tree and the HOGs it contains at the nodes
below; an orthogroup in the species tree (its gene tree mapped onto it: the
gene copies on each branch, the duplications and the losses); and each
species' genes, by orthogroup. Their data is in Tree_Data/HOGs_<node>.js and
Tree_Data/Species_<k>.js, loaded when opened.

All genes (a third tab): the species tree, with each orthogroup's gene tree
attached at the node where the orthogroup was gained, opened group by group
and orthogroup by orthogroup.
"""

import csv
import json
import os

from ete4 import Tree
from ete4.parser.newick import NewickError

from .. import __version__
from ..utils import view_pages

# colour-blind friendly, distinguishable on white (Okabe-Ito, then Tableau 10)
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255", "#117733",
           "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#B07AA1", "#9C755F", "#BAB0AC"]
DATA_DIR = "Tree_Data"


def species_and_genes(orthogroups_dir):
    """The species (the Orthogroups.tsv columns) and {gene: species} of all genes."""
    species, gene_species = [], {}
    for fn in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        path = os.path.join(orthogroups_dir, fn)
        if not os.path.exists(path):
            continue
        with open(path, newline="", encoding="utf-8", errors="replace") as infile:
            rows = list(csv.reader(infile, delimiter="\t"))
        species = rows[0][1:]
        for row in rows[1:]:
            for sp, cell in zip(species, row[1:]):
                for gene in cell.split(", ") if cell else []:
                    gene_species[gene] = sp
    return species, gene_species


def genes_by_species(orthogroups_dir, n_species):
    """
    For each species, its genes in each orthogroup ([[orthogroup, [genes]], ...])
    and its unassigned genes, from Orthogroups.tsv and Orthogroups_UnassignedGenes.tsv.
    """
    result = [{"orthogroups": [], "unassigned": []} for _ in range(n_species)]
    for fn, unassigned in (("Orthogroups.tsv", False), ("Orthogroups_UnassignedGenes.tsv", True)):
        path = os.path.join(orthogroups_dir, fn)
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as infile:
            next(infile, None)  # the header: the same species, in the same order
            for line in infile:
                row = line.rstrip("\r\n").split("\t")
                for k, cell in enumerate(row[1:n_species + 1]):
                    if cell:
                        genes = cell.split(", ")
                        if unassigned:
                            result[k]["unassigned"].extend(genes)
                        else:
                            result[k]["orthogroups"].append([row[0], genes])
    return result


def read_hogs(filename, species_number):
    """
    The HOGs of a species tree node: [[HOG, orthogroup, gene tree parent clade,
    [[species number, [genes]], ...]], ...], from Phylogenetic_Hierarchical_Orthogroups/<node>.tsv.
    """
    hogs = []
    with open(filename, encoding="utf-8", errors="replace") as infile:
        header = infile.readline().rstrip("\r\n").split("\t")
        columns = [(k, species_number[sp]) for k, sp in enumerate(header) if sp in species_number]
        for line in infile:
            row = line.rstrip("\r\n").split("\t")
            if len(row) < 3:
                continue
            genes = [[sp, row[k].split(", ")] for k, sp in columns if k < len(row) and row[k].strip()]
            hogs.append([row[0], row[1], row[2], genes])
    return hogs


def count_orthogroups(orthogroups_dir):
    """The number of orthogroups (Orthogroups.tsv), or None if it is not there."""
    path = os.path.join(orthogroups_dir, "Orthogroups.tsv")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8", errors="replace") as infile:
        return max(0, sum(1 for line in infile if line.strip()) - 1)


def write_explorer_data(results_dir, data_dir, species, species_number):
    """
    The data of the species tree explorer: the HOGs of each species tree node
    (HOGs_<node>.js) and the genes of each species (Species_<k>.js), each
    loaded by the page when it is opened. Returns [[node, HOGs], ...] and the
    genes of each species.
    """
    hog_nodes = []
    hog_dir = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    for fn in sorted(os.listdir(hog_dir)) if os.path.isdir(hog_dir) else []:
        if not fn.endswith(".tsv"):
            continue
        node = fn[:-4]
        hogs = read_hogs(os.path.join(hog_dir, fn), species_number)
        hog_nodes.append([node, len(hogs)])
        with open(os.path.join(data_dir, f"HOGs_{node}.js"), "w", encoding="utf-8") as out:
            out.write(f"hogsLoaded({json.dumps({'node': node, 'hogs': hogs}, separators=(',', ':'))});\n")
    n_genes = []
    for k, genes in enumerate(genes_by_species(os.path.join(results_dir, "Orthogroups"), len(species))):
        n_genes.append(sum(len(g) for _, g in genes["orthogroups"]) + len(genes["unassigned"]))
        with open(os.path.join(data_dir, f"Species_{k}.js"), "w", encoding="utf-8") as out:
            out.write(f"speciesLoaded({json.dumps({'k': k, **genes}, separators=(',', ':'))});\n")
    return hog_nodes, n_genes


def leaf_species(leaf_name, species_by_length, gene_species):
    """(species, gene name) of a gene tree leaf: "<species>_<gene>", or "<gene>" (-X)."""
    for sp in species_by_length:
        if leaf_name.startswith(sp + "_"):
            return sp, leaf_name[len(sp) + 1:]
    return gene_species.get(leaf_name), leaf_name


def read_alignment(filename):
    """{sequence name: aligned sequence} of a FASTA alignment."""
    seqs, name = {}, None
    with open(filename, encoding="utf-8", errors="replace") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            elif line and name is not None:
                seqs[name].append(line)
    return {name: "".join(parts) for name, parts in seqs.items()}


def read_duplications(filename):
    """
    {orthogroup: [[gene tree node, support, species tree node], ...]} from
    Duplications.tsv ({} if there is none).
    """
    dups = {}
    if not os.path.exists(filename):
        return dups
    with open(filename, encoding="utf-8", errors="replace") as infile:
        next(infile, None)  # the header
        for line in infile:
            # only the first fields: the gene lists can be very long
            og, species_node, node, support = line.split("\t", 4)[:4]
            dups.setdefault(og, []).append([node, float(support), species_node])
    return dups


def flat_tree(tree):
    """
    A tree as a list of its nodes in preorder (a parent before its children):
    [parent's index (-1 for the root), branch length, the leaf's number for a
    leaf or the name of an internal node], and the names of the leaves. A list,
    not nested lists, so that deep trees need no recursion, here or in the page.
    """
    index, nodes, leaf_names = {}, [], []
    for node in tree.traverse("preorder"):
        index[id(node)] = len(nodes)
        parent = index[id(node.up)] if node.up is not None else -1
        if node.is_leaf:
            nodes.append([parent, node.dist or 0, len(leaf_names)])
            leaf_names.append(node.name)
        else:
            nodes.append([parent, node.dist or 0, node.name or ""])
    return nodes, leaf_names


def alignment_stats(aln):
    """
    Of an alignment: its columns, the percentage of gaps, the shortest, median
    and longest sequence (without gaps) and the percentage of identical columns
    (the same residue in every sequence, no gaps).
    """
    seqs = [seq.upper() for seq in aln.values()]
    n_cols = max((len(seq) for seq in seqs), default=0)
    if not n_cols:
        return None
    lengths = sorted(len(seq) - seq.count("-") - seq.count(".") for seq in seqs)
    n_gaps = sum(n_cols - n for n in lengths)
    identical = sum(1 for col in zip(*seqs) if col[0] not in "-." and col.count(col[0]) == len(col))
    return {"columns": n_cols,
            "gaps": round(100 * n_gaps / (n_cols * len(seqs)), 1),
            "lengths": [lengths[0], lengths[len(lengths) // 2], lengths[-1]],
            "identical": round(100 * identical / n_cols, 1)}


def tree_stats(nodes, genes, duplications):
    """Of a gene tree: its genes per species ([species number, genes]), length and duplications."""
    per_species = {}
    for _, sp in genes:
        per_species[sp] = per_species.get(sp, 0) + 1
    return {"species": sorted(per_species.items()),
            "length": round(sum(node[1] for node in nodes), 4),
            "duplications": [len(duplications), sum(1 for d in duplications if d[1] >= 0.5)]}


def write_tree_data(filename, data):
    """The data of one gene tree, as a script that gives it to the page."""
    with open(filename, "w", encoding="utf-8") as out:
        out.write(f"treeLoaded({json.dumps(data, separators=(',', ':'))});\n")


class SpeciesTreeMap:
    """A species tree (flat, as flat_tree) for mapping gene trees onto it."""

    def __init__(self, nodes, leaf_species):
        """leaf_species: the species number of each of its leaves (-1 if not known)."""
        self.n = len(nodes)
        self.parent = [node[0] for node in nodes]
        self.kids = [[] for _ in nodes]
        self.level = [0] * self.n
        for k in range(1, self.n):
            self.kids[self.parent[k]].append(k)
            self.level[k] = self.level[self.parent[k]] + 1
        self.leaf_of = {leaf_species[node[2]]: k for k, node in enumerate(nodes)
                        if not self.kids[k] and leaf_species[node[2]] >= 0}

    def lca(self, a, b):
        while self.level[a] > self.level[b]:
            a = self.parent[a]
        while self.level[b] > self.level[a]:
            b = self.parent[b]
        while a != b:
            a, b = self.parent[a], self.parent[b]
        return a


def reconcile(gene_nodes, gene_species, smap):
    """
    A gene tree (flat) mapped onto the species tree, each of its nodes to the
    last common ancestor of its genes' species (gene_species: the species
    number of each gene). Returns, for each species tree node, the gene copies
    on the branch above it and the gene losses there, and the node its root is
    mapped to (None if no gene's species is in the species tree). As the
    page's reconcile() (the species tree explorer), for one orthogroup.
    """
    n = len(gene_nodes)
    kids = [[] for _ in range(n)]
    for k in range(1, n):
        kids[gene_nodes[k][0]].append(k)
    m = [-1] * n
    for k in range(n - 1, -1, -1):
        if not kids[k]:
            m[k] = smap.leaf_of.get(gene_species[gene_nodes[k][2]], -1)
            continue
        for c in kids[k]:
            if m[c] >= 0:
                m[k] = m[c] if m[k] < 0 else smap.lca(m[k], m[c])
    copies, losses = [0] * smap.n, [0] * smap.n
    for k in range(n):
        if m[k] < 0:
            continue
        p = gene_nodes[k][0]
        if p < 0:  # the lineage above the gene tree's root
            y = m[k]
            while y >= 0:
                copies[y] += 1
                y = smap.parent[y]
            continue
        if m[p] < 0 or m[p] == m[k]:
            continue
        y = m[k]
        while True:
            copies[y] += 1
            x = smap.parent[y]
            if x == m[p]:
                if any(m[c] == m[p] for c in kids[p]):  # a duplication at p: its other copy is lost below x
                    for z in smap.kids[x]:
                        if z != y:
                            losses[z] += 1
                break
            for z in smap.kids[x]:  # the lineage passes x without a copy in z
                if z != y:
                    losses[z] += 1
            y = x
    return copies, losses, (m[0] if m[0] >= 0 else None)


class AllOrthogroups:
    """The gene trees of all orthogroups mapped onto the species tree, added up node by node."""

    def __init__(self, smap):
        self.smap = smap
        self.copies, self.losses, self.gains = [0] * smap.n, [0] * smap.n, [0] * smap.n
        self.dups = [[0, 0] for _ in range(smap.n)]  # all, support >= 0.5
        # the orthogroups behind the numbers: gained, [orthogroup, losses], [orthogroup, duplications]
        self.gained = [[] for _ in range(smap.n)]
        self.lost = [[] for _ in range(smap.n)]
        self.duplicated = [[] for _ in range(smap.n)]
        self.n_trees = 0

    def add(self, og, gene_nodes, gene_species, og_dups, node_index):
        copies, losses, root = reconcile(gene_nodes, gene_species, self.smap)
        if root is None:
            return
        self.n_trees += 1
        self.gains[root] += 1
        self.gained[root].append(og)
        for k in range(self.smap.n):
            self.copies[k] += copies[k]
            if losses[k]:
                self.losses[k] += losses[k]
                self.lost[k].append([og, losses[k]])
        here = {}
        for _, support, species_node in og_dups:
            k = node_index.get(species_node)
            if k is not None:
                self.dups[k][0] += 1
                self.dups[k][1] += support >= 0.5
                here[k] = here.get(k, 0) + 1
        for k, count in here.items():
            self.duplicated[k].append([og, count])

    def summary(self, n_orthogroups):
        return {"n_trees": self.n_trees, "n_orthogroups": n_orthogroups, "copies": self.copies,
                "losses": self.losses, "gains": self.gains, "dups": self.dups}

    def lists(self):
        by_count = lambda rows: sorted(rows, key=lambda r: (-r[1], r[0]))  # the most first
        return {"gained": self.gained, "lost": [by_count(r) for r in self.lost],
                "duplicated": [by_count(r) for r in self.duplicated]}


def write_tree_pages(results_dir, written=None):
    """
    Write the page of the trees of a results directory and the data of each
    gene tree. Returns the number of gene trees (0 if the directory has no
    resolved gene trees). written: the pages a run is writing, to link to
    (view_pages.present; None: the pages that are there).
    """
    gene_trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    species_tree_fn = os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    if not os.path.exists(gene_trees_fn):
        return 0
    with open(gene_trees_fn, encoding="utf-8", errors="replace") as infile:
        entries = [line.strip().split(": ", 1) for line in infile if line.strip()]
    if not entries:  # e.g. -og: the file is there but has no trees
        return 0
    species, gene_species = species_and_genes(os.path.join(results_dir, "Orthogroups"))
    colour = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(species)}
    species_by_length = sorted(species, key=len, reverse=True)
    species_number = {sp: k for k, sp in enumerate(species)}
    msa_dir = os.path.join(results_dir, "MultipleSequenceAlignments")
    duplications = read_duplications(os.path.join(results_dir, "Gene_Duplication_Events", "Duplications.tsv"))
    data_dir = os.path.join(results_dir, "Resolved_Gene_Trees", DATA_DIR)

    def label_of(name):
        """["<species>: <gene>" or "<gene>", species number or -1] of a gene."""
        sp, gene = leaf_species(name, species_by_length, gene_species)
        return [f"{sp}: {gene}" if sp else gene, species_number.get(sp, -1)]

    species_tree, all_ogs = None, None
    if os.path.exists(species_tree_fn):
        with open(species_tree_fn, encoding="utf-8", errors="replace") as infile:
            tree = Tree(infile.read().strip(), parser=1)
        nodes, names = flat_tree(tree)
        species_tree = {"genes": [[name, species_number.get(name, -1)] for name in names], "nodes": nodes}
        all_ogs = AllOrthogroups(SpeciesTreeMap(nodes, [species_number.get(name, -1) for name in names]))
        node_index = {(names[node[2]] if isinstance(node[2], int) else node[2]): k for k, node in enumerate(nodes)}

    trees = []  # [orthogroup, [gene names], stats]
    n_failed = 0
    for og, newick in entries:
        try:
            tree = Tree(newick, parser=1)
        except (NewickError, ValueError):  # a tree that cannot be read
            n_failed += 1
            continue
        nodes, leaf_names = flat_tree(tree)
        genes = [label_of(name) for name in leaf_names]
        og_dups = duplications.get(og, [])
        data = {"og": og, "genes": genes, "nodes": nodes, "dups": og_dups}
        stats = tree_stats(nodes, genes, og_dups)
        msa_fn = os.path.join(msa_dir, og + ".fa")
        if os.path.exists(msa_fn):
            try:
                aln = read_alignment(msa_fn)
            except (ValueError, OSError):  # the tree is still shown, without its alignment
                aln = None
            if aln:
                # a sequence for each gene of the tree (None if it has none), then
                # the alignment's other sequences, which are not in the tree
                in_tree = set(leaf_names)
                extra = [name for name in aln if name not in in_tree]
                data["genes"] = genes + [label_of(name) for name in extra]
                data["seqs"] = [aln.get(name) for name in leaf_names] + [aln[name] for name in extra]
                stats["alignment"] = alignment_stats(aln)
        os.makedirs(data_dir, exist_ok=True)
        write_tree_data(os.path.join(data_dir, og + ".js"), data)
        trees.append([og, sorted({leaf_species(n, species_by_length, gene_species)[1] for n in leaf_names}), stats])
        if all_ogs:
            all_ogs.add(og, nodes, [sp for _, sp in genes], og_dups, node_index)
    if n_failed:
        print(f"WARNING: {n_failed} gene trees could not be read")
    if not trees:
        return 0

    hog_nodes, n_genes = write_explorer_data(results_dir, data_dir, species, species_number)
    explorer = {"hog_nodes": hog_nodes, "species_genes": n_genes}
    if all_ogs:
        explorer["all_orthogroups"] = all_ogs.summary(count_orthogroups(os.path.join(results_dir, "Orthogroups")))
        with open(os.path.join(data_dir, "All_Orthogroups.js"), "w", encoding="utf-8") as out:
            out.write(f"allOrthogroupsLoaded({json.dumps(all_ogs.lists(), separators=(',', ':'))});\n")
    explorer.update(nav=view_pages.nav(results_dir, "trees", written),
                    stats_page=view_pages.link(results_dir, "stats", written),
                    groups_page=view_pages.link(results_dir, "groups", written))
    write_page(os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.html"),
               species, colour, trees, species_tree, explorer)
    return len(trees)


def write_page(filename, species, colour, trees, species_tree, extra=None):
    """extra: the explorer's data and the links to the other pages (view_pages)."""
    data = {"version": __version__, "data": DATA_DIR, "species_tree": species_tree,
            "species": [[sp, colour[sp]] for sp in species], "trees": trees, **(extra or {})}
    with open(filename, "w", encoding="utf-8") as out:
        out.write(HTML.replace("__NAV_CSS__", view_pages.NAV_CSS).replace("__NAV_HTML__", view_pages.NAV_HTML)
                  .replace("__NAV_JS__", view_pages.NAV_JS).replace("__THEME_JS__", view_pages.THEME_JS)
                  .replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gene trees</title>
__THEME_JS__
<style>
:root { color-scheme:light; --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e8; --head:#f5f5f7; --accent:#2f6fde;
  --ink-fff:#fff; --ink-555:#555; --ink-333:#333; --ink-777:#777; --ink-2b3a67:#2b3a67; --ink-fff1a8:#fff1a8; --ink-c8c8c8:#c8c8c8; --ink-666:#666; --ink-999:#999; --ink-222:#222; --ink-444:#444; --ink-888:#888; --ink-aaa:#aaa; --ink-eeeef2:#eeeef2; --ink-ccc:#ccc; --ink-bbb:#bbb; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff;
  --ink-fff:#1b1b1e; --ink-555:#a0a0a8; --ink-333:#d0d0d6; --ink-777:#8f8f98; --ink-2b3a67:#9fb2ea; --ink-fff1a8:#5c4d00; --ink-c8c8c8:#45454b; --ink-666:#a6a6ad; --ink-999:#7d7d85; --ink-222:#e6e6ea; --ink-444:#c4c4ca; --ink-888:#8a8a92; --ink-aaa:#6f6f77; --ink-eeeef2:#2a2a2f; --ink-ccc:#3a3a40; --ink-bbb:#55555c; } }
:root[data-theme="dark"] { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff;
  --ink-fff:#1b1b1e; --ink-555:#a0a0a8; --ink-333:#d0d0d6; --ink-777:#8f8f98; --ink-2b3a67:#9fb2ea; --ink-fff1a8:#5c4d00; --ink-c8c8c8:#45454b; --ink-666:#a6a6ad; --ink-999:#7d7d85; --ink-222:#e6e6ea; --ink-444:#c4c4ca; --ink-888:#8a8a92; --ink-aaa:#6f6f77; --ink-eeeef2:#2a2a2f; --ink-ccc:#3a3a40; --ink-bbb:#55555c; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width:1400px; margin:0 auto; padding:24px clamp(16px, 2vw, 32px) 48px; }
h1 { font-size:32px; margin:0 0 6px; } h2 { font-size:19px; margin:26px 0 10px; }
#species > details > summary { margin:6px 0 8px; } #species > details > summary b { font-size:19px; }
.sub, .muted { color:var(--muted); font-size:12px; }
.legend { display:flex; flex-wrap:wrap; gap:6px 14px; margin:12px 0; font-size:13px; }
.legend span::before { content:""; display:inline-block; width:10px; height:10px; border-radius:2px; margin-right:6px; background:var(--c); }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:8px; }
input, select { font:inherit; padding:7px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
input { width:min(380px,100%); } select { max-width:100%; }
button { font:inherit; padding:6px 12px; border:1px solid var(--line); border-radius:6px; background:var(--head); color:var(--fg); cursor:pointer; }
button:disabled { opacity:.4; cursor:default; }
.tree { border:1px solid var(--line); border-radius:8px; padding:8px; background:var(--ink-fff); overflow:auto; }
.tree svg { display:block; }
.tree.capped { max-height:60vh; }
.tree [data-k] { cursor:pointer; }
.seg { display:inline-flex; } .seg button { border-radius:0; } .seg button:first-child { border-radius:6px 0 0 6px; }
.seg button:last-child { border-radius:0 6px 6px 0; border-left:0; } .seg button[aria-pressed="true"] { background:var(--accent); color:#fff; border-color:var(--accent); }
.key { font-size:12px; color:var(--muted); margin:0 0 6px; }
.key svg { vertical-align:-1px; }
.stats { display:grid; grid-template-columns:max-content 1fr; gap:4px 14px; margin:4px 0 12px; font-size:13px; }
.stats dt { color:var(--muted); } .stats dd { margin:0; }
.chip { display:inline-block; margin:0 10px 2px 0; white-space:nowrap; }
.chip::before { content:""; display:inline-block; width:9px; height:9px; border-radius:2px; margin-right:5px; background:var(--c); }
.chip.none { opacity:.45; }
.badge { display:inline-block; padding:0 7px; border-radius:9px; background:var(--head); border:1px solid var(--line); font-size:12px; }
.aln { border:1px solid var(--line); border-radius:8px; background:var(--ink-fff); overflow:auto; max-height:70vh; }
.aln canvas { position:sticky; top:0; left:0; display:block; }
#alnhover { font-family:ui-monospace, Menlo, Consolas, monospace; }
#maintabs { display:flex; gap:4px; margin:14px 0 4px; border-bottom:1px solid var(--line); }
#maintabs button { border:0; border-radius:0; background:none; padding:8px 12px; color:var(--muted); border-bottom:2px solid transparent; margin-bottom:-1px; }
#maintabs button[aria-selected="true"] { color:var(--fg); border-bottom-color:var(--accent); font-weight:600; }
section[hidden] { display:none; }
#ex-modes label { cursor:pointer; margin-right:10px; white-space:nowrap; }
#ex-modes input { width:auto; padding:0; margin:0 4px 0 0; vertical-align:-1px; }
.outline details { margin-left:16px; border-left:1px solid var(--line); padding-left:10px; }
.outline > details { margin-left:0; border-left:0; padding-left:0; }
.outline summary { cursor:pointer; padding:3px 0; }
.outline .leaf { margin-left:16px; padding:3px 0 3px 10px; border-left:1px solid var(--line); }
.outline .body { padding-bottom:4px; }
.outline .items > details, .outline .items > .gene { margin-left:0; }
.outline .gene { font-size:13px; padding:1px 0; }
.outline .sub { margin:4px 0 4px 0; }
.outline .more { margin:4px 0 8px; }
.outline .tools { margin:6px 0; }
.hogtree, .below { margin:6px 0; font-size:13px; }
.clade { overflow:auto; background:var(--ink-fff); border:1px solid var(--line); border-radius:6px; padding:4px; margin:4px 0; max-height:60vh; }
.clade svg { display:block; } .clade [data-k] { cursor:pointer; }
__NAV_CSS__
</style></head><body><main>
__NAV_HTML__
<h1>Gene trees</h1>
<div class="sub" id="sub"></div>
<div class="legend" id="legend"></div>
<nav id="maintabs" role="tablist"><button role="tab" data-tab="trees" aria-selected="true">Gene trees</button><button role="tab" data-tab="explore" aria-selected="false">Species tree explorer</button><button role="tab" data-tab="allgenes" aria-selected="false">All genes</button></nav>
<section id="tab-trees">
<div id="species"></div>
<h2>Resolved gene trees</h2>
<div class="tools">
  <input id="q" placeholder="Search an orthogroup ID or gene name" autocomplete="off">
  <select id="og"></select>
  <button id="prev">Previous</button><button id="next">Next</button>
  <span class="muted" id="count"></span>
  <a id="og-stats" hidden title="This orthogroup in the statistics page">in the statistics</a>
  <a id="og-groups" hidden title="The genes of this orthogroup and its HOGs">its genes</a>
</div>
<dl class="stats" id="stats"></dl>
<div class="tools">
  <span class="seg" role="group" aria-label="Layout"><button id="lay-rect" aria-pressed="true">Rectangular</button><button id="lay-circ" aria-pressed="false">Circular</button></span>
  <span class="seg" role="group" aria-label="Branch lengths" title="Branch lengths: to scale, or every branch the same length (to see the branching order when lengths crowd the tree)"><button id="len-scaled" aria-pressed="true">Scaled</button><button id="len-equal" aria-pressed="false">Equal</button></span>
  <button id="tout" title="Smaller">&minus;</button><button id="tin" title="Larger">+</button>
  <button id="expand" disabled>Expand all</button>
  <button id="save">Save as SVG</button>
  <button id="whole" hidden></button>
</div>
<div class="key" id="key"></div>
<div class="tree capped" id="treebox"></div>
<div id="alnbox" hidden>
<h2>Alignment <span class="muted" id="alninfo"></span></h2>
<div class="tools">
  <button id="zout" title="Zoom out">&minus;</button><button id="zin" title="Zoom in">+</button>
  <span class="muted" id="alnhover">Point at a residue to see its position</span>
</div>
<div class="aln" id="alnview"><div id="alnspace"><canvas id="cv"></canvas></div></div>
</div>
</section>
<section id="tab-allgenes" hidden>
<p class="muted" id="ag-help" style="margin-top:12px"></p>
<form class="tools" id="ag-find"><input id="ag-q" placeholder="Find a gene or orthogroup" autocomplete="off"><button type="submit">Find</button>
  <button type="button" id="ag-open">Open all groups</button><button type="button" id="ag-close">Close all</button><span class="muted" id="ag-msg"></span></form>
<div class="key">Click a group to open its orthogroups, an orthogroup to open its gene tree, a node to close or open it. Genes and branches in their species' colour;
  <svg width="10" height="10"><rect x="1" y="1" width="8" height="8" fill="#d62728"/></svg> gene duplication (support &ge; 0.5). Each node on its own row, its children below it.</div>
<div class="tree capped" id="ag-tree"></div>
</section>
<section id="tab-explore" hidden>
<div class="tools" id="ex-modes" style="margin-top:12px">
  <label><input type="radio" name="exmode" value="hogs" checked> HOGs at each node</label>
  <label><input type="radio" name="exmode" value="og"> Orthogroups in the species tree</label>
  <label><input type="radio" name="exmode" value="species"> Species and their genes</label>
</div>
<p class="muted" id="ex-help"></p>
<div id="ex-view"></div>
</section>
</main>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const num = (x, d = 0) => Number(x).toLocaleString(undefined, {maximumFractionDigits: d});
const MIXED = "var(--ink-555)", DUP = "#d62728", NODE = "var(--ink-2b3a67)";
const spColour = s => s >= 0 ? `var(--sp${s})` : MIXED;
// The drawings' colours are variables, so that they follow the light or dark colours: --ink-<colour on white>
// (the style) and --sp<species> (here, made lighter on dark until they read: contrast 4.5 with the dark panels).
// Save as SVG puts back the colours on white (lightSvg).
function lighten(hex) {
  if (!/^#[0-9a-f]{6}$/i.test(hex)) return hex;
  const rgb = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
  const lum = c => c.map(v => { v /= 255; return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4; })
    .reduce((a, v, i) => a + v * [0.2126, 0.7152, 0.0722][i], 0);
  const bg = lum([27, 27, 30]);  // --ink-fff, dark
  for (let t = 0; t <= 1; t += 0.05) {
    const c = rgb.map(v => Math.round(v + (255 - v) * t));
    if ((lum(c) + 0.05) / (bg + 0.05) >= 4.5) return "#" + c.map(v => v.toString(16).padStart(2, "0")).join("");
  }
  return "#ffffff";
}
{
  const light = D.species.map(([, c], k) => `--sp${k}:${c};`).join(""), dark = D.species.map(([, c], k) => `--sp${k}:${lighten(c)};`).join("");
  const st = document.createElement("style");
  st.textContent = `:root { ${light} } @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { ${dark} } } :root[data-theme="dark"] { ${dark} }`;
  document.head.appendChild(st);
}
const rootStyle = getComputedStyle(document.documentElement);
const ink = name => rootStyle.getPropertyValue("--ink-" + name).trim();  // for the alignment (a canvas)
const cssColour = c => c.startsWith("var(") ? rootStyle.getPropertyValue(c.slice(4, -1)).trim() : c;
const lightSvg = svg => svg.replace(/var\(--ink-([0-9a-f]+)\)/g, "#$1").replace(/var\(--sp(\d+)\)/g, (_, k) => D.species[+k][1]);
$("sub").textContent = "OrthoFinder " + D.version + ". The trees of Resolved_Gene_Trees.txt, genes coloured by species; "
  + "n0, n1, ... are the gene tree nodes (the HOG files' gene tree parent clade). "
  + "Below each tree, its alignment (MultipleSequenceAlignments).";
__NAV_JS__
$("legend").innerHTML = D.species.map(([s, c]) => `<span style="--c:${c}">${esc(s)}</span>`).join("");
$("key").innerHTML = "Click a node to collapse or expand its clade. A branch is in a species' colour when all its genes are of that species. "
  + `<svg width="10" height="10"><rect x="1" y="1" width="8" height="8" fill="${DUP}"/></svg> gene duplication (support &ge; 0.5), `
  + `<svg width="10" height="10"><rect x="1.5" y="1.5" width="7" height="7" fill="var(--ink-fff)" stroke="${DUP}" stroke-width="1.5"/></svg> support &lt; 0.5.`;

// ---- trees: drawn here, as SVG ----
// A tree is {genes: [[label, species number]], nodes: [[parent, branch length, gene number or node name]]},
// its nodes in preorder (a parent before its children).
const ROWH = 18, FONT = "12px system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif", TREE_ZOOMS = [0.35, 0.5, 0.7, 1, 1.4, 2, 3, 4.5, 7, 10, 15, 20];
const measurer = document.createElement("canvas").getContext("2d");
function textWidth(s, font = FONT) { measurer.font = font; return measurer.measureText(s).width; }

function prepareTree(T) {
  const n = T.nodes.length, kids = Array.from({length: n}, () => []);
  const depth = new Float64Array(n), sp = new Int32Array(n).fill(-2), nGenes = new Int32Array(n), deepest = new Float64Array(n);
  const udepth = new Float64Array(n), udeepest = new Float64Array(n);  // every branch of length 1
  const unit = T.nodes.every(([p, d]) => p < 0 || !d);  // no branch lengths: draw the topology
  for (let k = 0; k < n; k++) {
    const p = T.nodes[k][0];
    if (p >= 0) { kids[p].push(k); depth[k] = depth[p] + (unit ? 1 : T.nodes[k][1]); udepth[k] = udepth[p] + 1; }
  }
  for (let k = n - 1; k >= 0; k--) {
    if (!kids[k].length) { nGenes[k] = 1; sp[k] = T.genes[T.nodes[k][2]][1]; deepest[k] = depth[k]; udeepest[k] = udepth[k]; }
    const p = T.nodes[k][0];
    if (p < 0) continue;
    nGenes[p] += nGenes[k]; deepest[p] = Math.max(deepest[p], deepest[k]); udeepest[p] = Math.max(udeepest[p], udeepest[k]);
    sp[p] = sp[p] === -2 ? sp[k] : (sp[p] === sp[k] ? sp[p] : -1);
  }
  const dups = new Map((T.dups || []).map(([name, support]) => [name, support]));
  return Object.assign(T, {kids, depth, sp, nGenes, deepest, udepth, udeepest, unitLengths: unit, dupAt: dups, collapsed: new Set()});
}

function layoutTree(G) {
  // the rows (genes and collapsed clades) from top to bottom, and each shown node's row
  const rows = [], order = [], y = new Float64Array(G.nodes.length), stack = [0];
  while (stack.length) {
    const k = stack.pop();
    order.push(k);
    if (G.kids[k].length && !G.collapsed.has(k)) for (let j = G.kids[k].length - 1; j >= 0; j--) stack.push(G.kids[k][j]);
    else { y[k] = rows.length; rows.push(k); }
  }
  for (let j = order.length - 1; j >= 0; j--) {
    const k = order[j], c = G.kids[k];
    if (c.length && !G.collapsed.has(k)) y[k] = (y[c[0]] + y[c[c.length - 1]]) / 2;
  }
  return {rows, order, y};
}

function cladeLabel(G, k) {
  const name = G.nodes[k][2] || "clade", n = G.nGenes[k];
  if (G.sp[k] >= 0) return `${name}: ${num(n)} ${D.species[G.sp[k]][0]} gene${n === 1 ? "" : "s"}`;
  const species = new Set(), stack = [k];
  while (stack.length) { const j = stack.pop(); if (G.kids[j].length) stack.push(...G.kids[j]); else species.add(G.sp[j]); }
  return `${name}: ${num(n)} genes, ${species.size} species`;
}
function rowLabel(G, k) { return G.kids[k].length ? "▸ " + cladeLabel(G, k) : G.genes[G.nodes[k][2]][0]; }

function niceLength(x) {
  const p = Math.pow(10, Math.floor(Math.log10(x))), m = x / p;
  return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p;
}

function treeSvg(G, opts) {
  const dep = opts.equal ? G.udepth : G.depth, deep = opts.equal ? G.udeepest : G.deepest;
  const L = layoutTree(G), {rows, order, y} = L;
  const labels = rows.map(k => rowLabel(G, k));
  const labelW = Math.max(40, ...labels.map(s => textWidth(s)));
  const maxDepth = deep[0] || 1;
  const treeW = Math.max(120, Math.min(1400, opts.width - labelW - 60)) * TREE_ZOOMS[opts.zoom ?? 3];
  const sx = treeW / maxDepth, left = 14, top = 6;
  const X = k => left + dep[k] * sx, Y = r => top + r * ROWH + ROWH / 2;
  const q = (opts.highlight || "").toLowerCase();
  const out = [];
  // branches: horizontal to each node, vertical across each expanded node's children
  for (const k of order) {
    const p = G.nodes[k][0], col = spColour(G.sp[k]), yk = Y(y[k]);
    out.push(`<line x1="${p >= 0 ? X(p) : left - 10}" y1="${yk}" x2="${X(k)}" y2="${yk}" stroke="${col}" stroke-width="1.6"/>`);
    const c = G.kids[k];
    if (c.length && !G.collapsed.has(k))
      out.push(`<line x1="${X(k)}" y1="${Y(y[c[0]]) - 0.8}" x2="${X(k)}" y2="${Y(y[c[c.length - 1]]) + 0.8}" stroke="${col}" stroke-width="1.6"/>`);
  }
  // genes and collapsed clades
  rows.forEach((k, r) => {
    const col = spColour(G.sp[k]), yk = Y(r), label = labels[r];
    let x = X(k);
    if (G.kids[k].length) {
      const x2 = left + deep[k] * sx;
      out.push(`<g data-k="${k}"><title>${esc(cladeLabel(G, k))} (click to expand)</title>`
        + `<polygon points="${x},${yk} ${x2},${yk - 6} ${x2},${yk + 6}" fill="${col}" fill-opacity="0.35" stroke="${col}"/></g>`);
      x = x2;
    } else {
      out.push(`<circle cx="${x}" cy="${yk}" r="2.5" fill="${col}"/>`);
    }
    if (q && label.toLowerCase().includes(q))
      out.push(`<rect x="${x + 4}" y="${yk - ROWH / 2 + 1}" width="${textWidth(label) + 4}" height="${ROWH - 2}" fill="var(--ink-fff1a8)"/>`);
    out.push(`<text x="${x + 6}" y="${yk}" dy="0.35em" fill="${col}" font-size="12" font-family="system-ui, sans-serif"`
      + (G.kids[k].length ? ` data-k="${k}"` : "") + `>${esc(label)}</text>`);
  });
  // internal nodes: their names, and a marker (a square for a gene duplication) to click
  for (const k of order) {
    if (!G.kids[k].length || G.collapsed.has(k)) continue;
    const x = X(k), yk = Y(y[k]), name = G.nodes[k][2], support = G.dupAt.get(name);
    if (name) out.push(`<text x="${x - 3}" y="${yk - 4}" text-anchor="end" fill="var(--ink-777)" font-size="9" font-family="system-ui, sans-serif">${esc(name)}</text>`);
    const tip = `${cladeLabel(G, k)}` + (support !== undefined ? `; gene duplication, support ${support}` : "") + (opts.clickable ? " (click to collapse)" : "");
    const mark = support === undefined ? `<circle cx="${x}" cy="${yk}" r="3" fill="${NODE}"/>`
      : `<rect x="${x - 4}" y="${yk - 4}" width="8" height="8" fill="${support >= 0.5 ? DUP : "var(--ink-fff)"}" stroke="${DUP}" stroke-width="1.5"/>`;
    out.push(opts.clickable ? `<g data-k="${k}"><title>${esc(tip)}</title><circle cx="${x}" cy="${yk}" r="8" fill="transparent"/>${mark}</g>`
                            : `<g><title>${esc(tip)}</title>${mark}</g>`);
  }
  // scale bar
  const height = top + rows.length * ROWH + 34, bar = niceLength(maxDepth / 5);
  if (!G.unitLengths && !opts.equal) {
    const by = height - 18;
    out.push(`<line x1="${left}" y1="${by}" x2="${left + bar * sx}" y2="${by}" stroke="var(--ink-333)"/>`
      + `<line x1="${left}" y1="${by - 4}" x2="${left}" y2="${by + 4}" stroke="var(--ink-333)"/><line x1="${left + bar * sx}" y1="${by - 4}" x2="${left + bar * sx}" y2="${by + 4}" stroke="var(--ink-333)"/>`
      + `<text x="${left}" y="${by + 14}" font-size="11" font-family="system-ui, sans-serif" fill="var(--ink-333)">${num(bar, 6)}</text>`);
  }
  const width = Math.ceil(left + treeW + 12 + labelW + 16);
  return {svg: `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`
    + `<rect width="100%" height="100%" fill="var(--ink-fff)"/>${out.join("")}</svg>`, rows};
}

// The same tree in a circle: the rows (genes and collapsed clades) clockwise from the top,
// a node's distance from the centre its distance from the root.
function circularTreeSvg(G, opts) {
  const dep = opts.equal ? G.udepth : G.depth, deep = opts.equal ? G.udeepest : G.deepest;
  const L = layoutTree(G), {rows, order, y} = L, n = Math.max(1, rows.length);
  const labels = rows.map(k => rowLabel(G, k));
  const labelW = Math.max(40, ...labels.map(s => textWidth(s)));
  const maxDepth = deep[0] || 1;
  const spaced = n * ROWH / (2 * Math.PI);  // room for a label in each row at the outside
  const R = Math.max(90, spaced, Math.min(1400, opts.width) / 2 - labelW - 24) * TREE_ZOOMS[opts.zoom ?? 3];
  const c = R + labelW + 24, sr = R / maxDepth;  // the centre (the picture is cropped to what is drawn)
  const rad = k => dep[k] * sr, ang = row => -Math.PI / 2 + 2 * Math.PI * row / n;
  const P = (r, t) => `${(c + r * Math.cos(t)).toFixed(2)},${(c + r * Math.sin(t)).toFixed(2)}`;
  const arc = (r, t1, t2) => `M${P(r, t1)}A${r.toFixed(2)},${r.toFixed(2)} 0 ${t2 - t1 > Math.PI ? 1 : 0} 1 ${P(r, t2)}`;
  const q = (opts.highlight || "").toLowerCase(), out = [];
  let x0 = c - 20, x1 = c + 20, y0 = c - 20, y1 = c + 20;  // how far the drawing reaches
  const reach = (r, t) => { const x = c + r * Math.cos(t), yy = c + r * Math.sin(t);
    x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, yy); y1 = Math.max(y1, yy); };
  for (const k of order) {  // branches: out from the centre to each node, round each expanded node's children
    const p = G.nodes[k][0], col = spColour(G.sp[k]), t = ang(y[k]);
    if (p >= 0) out.push(`<path d="M${P(rad(p), t)}L${P(rad(k), t)}" stroke="${col}" stroke-width="1.6" fill="none"/>`);
    const ch = G.kids[k];
    if (ch.length && !G.collapsed.has(k) && rad(k) > 0)
      out.push(`<path d="${arc(rad(k), ang(y[ch[0]]), ang(y[ch[ch.length - 1]]))}" stroke="${col}" stroke-width="1.6" fill="none"/>`);
  }
  const dt = Math.min(0.4 * 2 * Math.PI / n, 0.3);
  rows.forEach((k, r) => {  // genes and collapsed clades, their labels turned to their angle
    const col = spColour(G.sp[k]), t = ang(r), label = labels[r];
    let rEnd = rad(k);
    if (G.kids[k].length) {
      rEnd = deep[k] * sr;
      out.push(`<g data-k="${k}"><title>${esc(cladeLabel(G, k))} (click to expand)</title>`
        + `<path d="M${P(rad(k), t)}L${P(rEnd, t - dt)}A${rEnd},${rEnd} 0 0 1 ${P(rEnd, t + dt)}Z" fill="${col}" fill-opacity="0.35" stroke="${col}"/></g>`);
    } else out.push(`<circle cx="${(c + rEnd * Math.cos(t)).toFixed(2)}" cy="${(c + rEnd * Math.sin(t)).toFixed(2)}" r="2.5" fill="${col}"/>`);
    const left = Math.cos(t) < 0, deg = t * 180 / Math.PI + (left ? 180 : 0), rl = rEnd + 6;
    reach(rl + textWidth(label) + 4, t); reach(rl + textWidth(label) + 4, t + 0.02); reach(rl + textWidth(label) + 4, t - 0.02);
    const hit = q && label.toLowerCase().includes(q);
    out.push(`<text transform="translate(${P(rl, t)}) rotate(${deg.toFixed(2)})" dy="0.35em" text-anchor="${left ? "end" : "start"}" fill="${col}"`
      + ` font-size="12" font-family="system-ui, sans-serif"${hit ? ` stroke="var(--ink-fff1a8)" stroke-width="7" paint-order="stroke" font-weight="600"` : ""}`
      + (G.kids[k].length ? ` data-k="${k}"` : "") + `>${esc(label)}</text>`);
  });
  for (const k of order) {  // internal nodes: a marker (a square for a gene duplication) to click
    if (!G.kids[k].length || G.collapsed.has(k)) continue;
    const t = ang(y[k]), x = c + rad(k) * Math.cos(t), yy = c + rad(k) * Math.sin(t), name = G.nodes[k][2], support = G.dupAt.get(name);
    const tip = `${name ? name + ": " : ""}${cladeLabel(G, k)}` + (support !== undefined ? `; gene duplication, support ${support}` : "") + (opts.clickable ? " (click to collapse)" : "");
    const mark = support === undefined ? `<circle cx="${x}" cy="${yy}" r="3" fill="${NODE}"/>`
      : `<rect x="${x - 4}" y="${yy - 4}" width="8" height="8" fill="${support >= 0.5 ? DUP : "var(--ink-fff)"}" stroke="${DUP}" stroke-width="1.5"/>`;
    out.push(opts.clickable ? `<g data-k="${k}"><title>${esc(tip)}</title><circle cx="${x}" cy="${yy}" r="8" fill="transparent"/>${mark}</g>`
                            : `<g><title>${esc(tip)}</title>${mark}</g>`);
    if (name && n <= 200)  // its name (N0, N1, ... in the species tree; n0, n1, ... in a gene tree), unless there are too many
      out.push(`<text x="${(x + 5).toFixed(1)}" y="${(yy - 5).toFixed(1)}" fill="var(--ink-777)" font-size="9" font-family="system-ui, sans-serif">${esc(name)}</text>`);
  }
  // the picture: what the drawing reaches, and below it the scale bar
  x0 = Math.floor(x0 - 10); y0 = Math.floor(y0 - 10); x1 = Math.ceil(x1 + 10); y1 = Math.ceil(y1 + 10);
  if (!G.unitLengths && !opts.equal) {
    const bar = niceLength(maxDepth / 5), bx = x0 + 14, by = y1 + 6;
    out.push(`<line x1="${bx}" y1="${by}" x2="${bx + bar * sr}" y2="${by}" stroke="var(--ink-333)"/>`
      + `<line x1="${bx}" y1="${by - 4}" x2="${bx}" y2="${by + 4}" stroke="var(--ink-333)"/><line x1="${bx + bar * sr}" y1="${by - 4}" x2="${bx + bar * sr}" y2="${by + 4}" stroke="var(--ink-333)"/>`
      + `<text x="${bx}" y="${by + 14}" font-size="11" font-family="system-ui, sans-serif" fill="var(--ink-333)">${num(bar, 6)}</text>`);
    y1 += 30;
  }
  const w = x1 - x0, h = y1 - y0;
  return {svg: `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="${x0} ${y0} ${w} ${h}" style="margin:0 auto">`
    + `<rect x="${x0}" y="${y0}" width="${w}" height="${h}" fill="var(--ink-fff)"/>${out.join("")}</svg>`, rows};
}

if (D.species_tree) {
  // its own layout switch (remembered apart from the gene trees')
  const S = prepareTree(D.species_tree);
  let layout = "rect";
  try { if (localStorage.getItem("orthofinder-species-tree-layout") === "circ") layout = "circ"; } catch (e) { /* not available */ }
  $("species").innerHTML = `<details open><summary><b>Species tree</b> <span class="muted">(N0, N1, ... as in the HOG files)</span></summary>`
    + `<div class="tools" style="margin-top:6px"><span class="seg" role="group" aria-label="Layout of the species tree">`
    + `<button id="sp-rect">Rectangular</button><button id="sp-circ">Circular</button></span></div>`
    + `<div class="tree" id="species-tree"></div></details>`;
  const drawSpeciesTree = () => {
    const draw = layout === "circ" ? circularTreeSvg : treeSvg;
    $("species-tree").innerHTML = draw(S, {width: Math.min(700, $("species").clientWidth), clickable: false}).svg;
    $("sp-rect").setAttribute("aria-pressed", layout === "rect"); $("sp-circ").setAttribute("aria-pressed", layout === "circ");
  };
  const setSpeciesLayout = l => {
    layout = l;
    try { localStorage.setItem("orthofinder-species-tree-layout", l); } catch (e) { /* not available */ }
    drawSpeciesTree();
  };
  $("sp-rect").onclick = () => setSpeciesLayout("rect");
  $("sp-circ").onclick = () => setSpeciesLayout("circ");
  drawSpeciesTree();
}

// ---- the gene tree shown ----
const loaded = {};
let shown = D.trees, i = 0, G = null, wanted = null, treeZoom = 3, wholeTree = false, treeLayout = "rect";
let treeEqual = false;  // every branch the same length
try { if (localStorage.getItem("orthofinder-tree-layout") === "circ") treeLayout = "circ"; } catch (e) { /* not available */ }
try { treeEqual = localStorage.getItem("orthofinder-tree-lengths") === "equal"; } catch (e) { /* not available */ }
window.treeLoaded = T => { loaded[T.og] = T; gotData("tree:" + T.og, T); if (T.og === wanted) setTree(T); };
function show() {
  $("og").innerHTML = shown.map((t, k) => `<option value="${k}">${esc(t[0])} (${t[1].length} genes)</option>`).join("");
  $("og").value = String(i);
  const t = shown[i];
  $("count").textContent = shown.length ? `${i + 1} of ${shown.length.toLocaleString()}` : "no matches";
  $("prev").disabled = i <= 0;
  $("next").disabled = i >= shown.length - 1;
  showStats(t);
  $("og-stats").hidden = !(t && D.stats_page);
  if (t && D.stats_page) $("og-stats").href = `${D.stats_page}#orthogroups/${encodeURIComponent(t[0])}`;
  $("og-groups").hidden = !(t && D.groups_page);
  if (t && D.groups_page) $("og-groups").href = `${D.groups_page}#orthogroups/${encodeURIComponent(t[0])}`;
  if (!t) { wanted = null; G = null; $("treebox").innerHTML = ""; $("alnbox").hidden = true; return; }
  wanted = t[0];
  if (loaded[t[0]]) return setTree(loaded[t[0]]);
  $("treebox").innerHTML = `<p class="muted">Loading ${esc(t[0])}...</p>`;
  const s = document.createElement("script");
  s.src = `${D.data}/${encodeURIComponent(t[0])}.js`;
  s.onerror = () => { if (wanted === t[0]) $("treebox").innerHTML = `<p class="muted">The data of ${esc(t[0])} could not be loaded (${esc(s.src)}). The page needs its folder ${esc(D.data)} next to it: copy (or download) them together, e.g. the whole results directory.</p>`; };
  document.head.appendChild(s);
}
function setTree(T) {
  G = T.kids ? T : prepareTree(T);
  $("treebox").scrollTo(0, 0);
  drawTree(true);
}
function drawTree(newTree) {
  if (!G) return;
  const draw = treeLayout === "circ" ? circularTreeSvg : treeSvg;
  const box = $("treebox"), {svg, rows} = draw(G, {width: box.clientWidth, zoom: treeZoom, clickable: true, highlight: $("q").value.trim(), equal: treeEqual});
  $("len-scaled").setAttribute("aria-pressed", !treeEqual); $("len-equal").setAttribute("aria-pressed", treeEqual);
  $("lay-rect").setAttribute("aria-pressed", treeLayout === "rect"); $("lay-circ").setAttribute("aria-pressed", treeLayout === "circ");
  box.innerHTML = svg;
  if (treeLayout === "circ" && newTree) centreTree();
  $("expand").disabled = !G.collapsed.size;
  $("tout").disabled = treeZoom === 0; $("tin").disabled = treeZoom === TREE_ZOOMS.length - 1;
  treeCap();
  if (G.seqs) setAlignment(alignmentRows(G, rows), newTree); else { aln = null; $("alnbox").hidden = true; }
}
$("treebox").addEventListener("click", e => {
  const el = e.target.closest("[data-k]");
  if (!el || !G) return;
  const k = +el.dataset.k;
  if (G.collapsed.has(k)) G.collapsed.delete(k); else G.collapsed.add(k);
  drawTree(false);
});
function setLayout(layout) {
  treeLayout = layout;
  try { localStorage.setItem("orthofinder-tree-layout", layout); } catch (e) { /* not available */ }
  $("treebox").scrollTo(0, 0);
  drawTree(false);
  if (layout === "circ") centreTree();
}
function centreTree() {  // a circular tree larger than its box: show its middle
  const box = $("treebox");
  box.scrollTo((box.scrollWidth - box.clientWidth) / 2, (box.scrollHeight - box.clientHeight) / 2);
}
$("lay-rect").onclick = () => setLayout("rect");
function setLengths(equal) {
  treeEqual = equal;
  try { localStorage.setItem("orthofinder-tree-lengths", equal ? "equal" : "scaled"); } catch (e) { /* not available */ }
  drawTree(false);
  if (treeLayout === "circ") centreTree();
}
$("len-scaled").onclick = () => setLengths(false);
$("len-equal").onclick = () => setLengths(true);
$("lay-circ").onclick = () => setLayout("circ");
$("expand").onclick = () => { if (G) { G.collapsed.clear(); drawTree(false); } };
$("tin").onclick = () => { treeZoom = Math.min(TREE_ZOOMS.length - 1, treeZoom + 1); drawTree(false); };
$("tout").onclick = () => { treeZoom = Math.max(0, treeZoom - 1); drawTree(false); };
$("save").onclick = () => {
  const svg = $("treebox").querySelector("svg");
  if (!svg || !G) return;
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + lightSvg(svg.outerHTML)], {type: "image/svg+xml"}));
  a.download = `${G.og}.svg`;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
};

// A tall tree is shown in a box of 60% of the window's height, so that its
// alignment is close; "Show the whole tree" lifts the limit (for all trees).
function treeCap() {
  const svg = $("treebox").querySelector("svg");
  $("treebox").classList.toggle("capped", !wholeTree);
  $("whole").hidden = !svg || svg.getBoundingClientRect().height + 18 <= window.innerHeight * 0.6;
  $("whole").textContent = wholeTree ? "Limit the tree's height" : "Show the whole tree";
}
window.addEventListener("resize", treeCap);
$("whole").onclick = () => { wholeTree = !wholeTree; treeCap(); };

function showStats(t) {
  if (!t) { $("stats").innerHTML = ""; return; }
  const st = t[2], counts = Object.fromEntries(st.species), n = t[1].length;
  const present = D.species.filter((_, k) => counts[k]).length;
  const singleCopy = present === D.species.length && D.species.every((_, k) => counts[k] === 1);
  const chips = D.species.map(([sp, c], k) =>
    `<span class="chip${counts[k] ? "" : " none"}" style="--c:${c}">${esc(sp)} ${counts[k] || 0}</span>`).join("");
  const [dups, supported] = st.duplications;
  let html = `<dt>Orthogroup</dt><dd>${num(n)} genes, ${present} of ${D.species.length} species`
    + (singleCopy ? ` <span class="badge">single-copy in all species</span>` : "") + `</dd>`
    + `<dt>Genes per species</dt><dd>${chips}</dd>`
    + `<dt>Gene tree</dt><dd>total branch length ${num(st.length, 3)}; `
    + `${num(dups)} gene duplication${dups === 1 ? "" : "s"} (${num(supported)} with support &ge; 0.5)</dd>`;
  const a = st.alignment;
  if (a) html += `<dt>Alignment</dt><dd>${num(a.columns)} columns, ${num(a.gaps, 1)}% gaps, `
    + `${num(a.identical, 1)}% identical columns; sequence lengths ${num(a.lengths[0])} to ${num(a.lengths[2])} (median ${num(a.lengths[1])})</dd>`;
  $("stats").innerHTML = html;
}

// ---- the alignment of the tree shown: its rows follow the tree's rows ----
function alignmentRows(G, treeRows) {
  const rows = treeRows.map(k => G.kids[k].length
    ? [rowLabel(G, k), G.sp[k], null]
    : [...G.genes[G.nodes[k][2]], G.seqs[G.nodes[k][2]]]);
  const nLeaves = G.nodes.filter(nd => typeof nd[2] === "number").length;
  for (let g = nLeaves; g < G.genes.length; g++) rows.push([...G.genes[g], G.seqs[g]]);  // not in the tree
  return rows;
}
const ROW = 16, RULER = 20, ZOOMS = [1, 2, 3, 5, 8, 11, 14, 18, 24];
const START_SPACE = 14, END_SPACE = 14;  // the space before the first column and after the last (the same)
const AA = {}, NT = {};
for (const [cs, col] of [["AILMFWV", "#b7d0ff"], ["KR", "#ffb8b8"], ["DE", "#f3c0ff"], ["NQST", "#bdf0bd"],
                         ["C", "#ffd9e8"], ["G", "#ffd9a6"], ["P", "#fff0a0"], ["HY", "#aef0ee"]])
  for (const c of cs) AA[c] = col;
for (const [c, col] of [["A", "#a8e6a8"], ["C", "#a9cdff"], ["G", "#ffd598"], ["T", "#ffacac"], ["U", "#ffacac"]]) NT[c] = col;
let aln = null, zoom = 4, nameW = 0, numW = 0, scheme = AA, drawPending = false;  // numW: the row numbers' column
function setAlignment(rows, reset) {
  const seqs = rows.map(r => r[2]).filter(s => s);
  aln = {rows, ncol: Math.max(0, ...seqs.map(s => s.length))};
  const residues = seqs.slice(0, 20).join("").replace(/[-.]/g, "").toUpperCase();
  scheme = residues && /^[ACGTUN]*$/.test(residues) ? NT : AA;
  $("alnbox").hidden = false;
  $("alninfo").textContent = `(${seqs.length} sequences, ${aln.ncol.toLocaleString()} columns; rows in the order of the tree${treeLayout === "circ" ? ", clockwise from the top" : ""})`;
  numW = textWidth(String(rows.length), "11px system-ui, sans-serif") + 14;
  nameW = numW + Math.min(Math.max(...rows.map(r => textWidth(r[0])), 40) + 16, Math.max(160, $("alnview").clientWidth * 0.45));
  if (reset) $("alnview").scrollTo(0, 0);
  layout();
}
function layout() {
  if (!aln) return;
  const cw = ZOOMS[zoom];
  $("alnspace").style.width = (nameW + START_SPACE + aln.ncol * cw + END_SPACE) + "px";
  $("alnspace").style.height = (RULER + aln.rows.length * ROW) + "px";
  $("zout").disabled = zoom === 0; $("zin").disabled = zoom === ZOOMS.length - 1;
  draw();
}
function fit(ctx, text, width) {
  if (ctx.measureText(text).width <= width) return text;
  let lo = 0, hi = text.length;
  while (lo < hi) { const m = (lo + hi + 1) >> 1; if (ctx.measureText(text.slice(0, m) + "…").width <= width) lo = m; else hi = m - 1; }
  return text.slice(0, lo) + "…";
}
function draw() {
  drawPending = false;
  if (!aln) return;
  const view = $("alnview"), cv = $("cv"), space = $("alnspace"), dpr = window.devicePixelRatio || 1;
  const w = Math.min(view.clientWidth, space.offsetWidth), h = Math.min(view.clientHeight, space.offsetHeight);
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    cv.style.width = w + "px"; cv.style.height = h + "px";
  }
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = ink("fff"); ctx.fillRect(0, 0, w, h);
  const x0 = view.scrollLeft, y0 = view.scrollTop, cw = ZOOMS[zoom], rows = aln.rows, ncol = aln.ncol;
  const left = nameW + START_SPACE - x0, colX = c => left + c * cw;  // where column c (from 0) starts
  const c0 = Math.max(0, Math.floor((nameW - left) / cw)), c1 = Math.min(ncol, Math.ceil((w - left) / cw) + 1);
  const r0 = Math.floor(y0 / ROW), r1 = Math.min(rows.length, Math.ceil((y0 + h - RULER) / ROW) + 1);
  ctx.textAlign = "center"; ctx.textBaseline = "middle";
  ctx.font = `${Math.min(12, cw + 1)}px ui-monospace, Menlo, Consolas, monospace`;
  for (let r = r0; r < r1; r++) {
    const seq = rows[r][2], y = RULER + r * ROW - y0;
    if (seq === null) {  // a collapsed clade, or a gene of the tree without a sequence
      ctx.fillStyle = ink("eeeef2"); ctx.fillRect(colX(c0), y + 1, (c1 - c0) * cw, ROW - 2);
      continue;
    }
    for (let c = c0; c < c1; c++) {
      const ch = seq[c];
      if (ch === "-" || ch === "." || ch === undefined) continue;
      const x = colX(c);
      ctx.fillStyle = scheme[ch.toUpperCase()] || "#e4e4e4";
      ctx.fillRect(x, y + 1, cw, ROW - 2);
      if (cw >= 8) { ctx.fillStyle = "#222"; ctx.fillText(ch, x + cw / 2, y + ROW / 2); }
    }
  }
  // the ruler (column numbers) and the gene names stay in place when scrolling
  ctx.fillStyle = ink("fff"); ctx.fillRect(0, 0, w, RULER);
  ctx.strokeStyle = ink("ccc"); ctx.beginPath(); ctx.moveTo(0, RULER - 0.5); ctx.lineTo(w, RULER - 0.5); ctx.stroke();
  // the column numbers: as many as fit (every column, 2, 5, 10, ... as the zoom allows), and the first and last
  ctx.font = "10px system-ui, sans-serif"; ctx.textAlign = "center";
  const labelW = ctx.measureText(String(ncol)).width + 8;
  const step = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].find(s => s * cw >= labelW) || 10000;
  const mid = c => colX(c - 1) + cw / 2;  // the middle of column c (from 1)
  ctx.fillStyle = ink("999");
  if (cw >= 5) for (let c = c0 + 1; c <= c1; c++) if (colX(c - 1) >= nameW) ctx.fillRect(Math.round(mid(c)), RULER - 3, 1, 2);
  ctx.fillStyle = ink("666");
  const label = c => { const x = mid(c); if (x >= nameW) { ctx.fillRect(Math.round(x), RULER - 6, 1, 5); ctx.fillText(String(c), x, RULER / 2 - 1); } };
  for (let c = Math.max(step, Math.ceil((c0 + 1) / step) * step); c <= c1; c += step)
    if ((c - 1) * cw >= labelW && (ncol - c) * cw >= labelW) label(c);  // not over the first or last column's number
  label(1);
  if (ncol > 1) label(ncol);
  ctx.fillStyle = ink("bbb");  // the start and end of the alignment
  if (colX(0) > nameW) ctx.fillRect(Math.round(colX(0)) - 2, RULER, 1, h - RULER);
  if (colX(ncol) < w) ctx.fillRect(Math.round(colX(ncol)) + 1, RULER, 1, h - RULER);
  ctx.textAlign = "left";
  ctx.fillStyle = ink("fff"); ctx.fillRect(0, 0, nameW, h);
  ctx.font = FONT;
  for (let r = r0; r < r1; r++) {  // the row numbers, then the names
    const [label, sp] = rows[r], y = RULER + r * ROW - y0;
    ctx.font = "11px system-ui, sans-serif"; ctx.fillStyle = ink("999"); ctx.textAlign = "right";
    ctx.fillText(String(r + 1), numW - 8, y + ROW / 2);
    ctx.font = FONT; ctx.fillStyle = cssColour(spColour(sp)); ctx.textAlign = "left";
    ctx.fillText(fit(ctx, label, nameW - numW - 14), numW, y + ROW / 2);
  }
  ctx.strokeStyle = ink("ccc"); ctx.beginPath(); ctx.moveTo(nameW - 4.5, 0); ctx.lineTo(nameW - 4.5, h); ctx.stroke();
}
function redraw() { if (!drawPending) { drawPending = true; requestAnimationFrame(draw); } }
$("alnview").addEventListener("scroll", redraw);
window.addEventListener("resize", redraw);
window.addEventListener("orthofinder-theme", redraw);  // Colours: Light, Dark or Auto
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
$("zin").onclick = () => { zoomTo(zoom + 1); };
$("zout").onclick = () => { zoomTo(zoom - 1); };
function zoomTo(z) {
  if (!aln || z < 0 || z >= ZOOMS.length) return;
  // keep the column in the middle of the view in the middle
  const view = $("alnview"), centre = (nameW + view.clientWidth) / 2;
  const col = (view.scrollLeft + centre - nameW - START_SPACE) / ZOOMS[zoom];
  zoom = z; layout();
  view.scrollLeft = nameW + START_SPACE + col * ZOOMS[zoom] - centre;
}
$("cv").addEventListener("mousemove", e => {
  if (!aln) return;
  const view = $("alnview"), cw = ZOOMS[zoom];
  const c = Math.floor((e.offsetX - nameW - START_SPACE + view.scrollLeft) / cw), r = Math.floor((e.offsetY - RULER + view.scrollTop) / ROW);
  const row = aln.rows[r];
  if (e.offsetX < nameW || e.offsetY < RULER || !row || !row[2] || c < 0 || c >= row[2].length) {
    $("alnhover").textContent = row && e.offsetY >= RULER ? `row ${r + 1} \u00B7 ${row[0]}` : "Point at a residue to see its position";
    return;
  }
  const ch = row[2][c], n = row[2].slice(0, c + 1).replace(/[-.]/g, "").length;
  $("alnhover").textContent = `row ${r + 1} · ${row[0]} · column ${c + 1} · ` + (ch === "-" || ch === "." ? "gap" : `residue ${n}: ${ch}`);
});

// ---- choosing the tree ----
$("q").addEventListener("input", e => {
  const q = e.target.value.trim().toLowerCase();
  const cur = shown[i] && shown[i][0];
  shown = q ? D.trees.filter(([og, genes]) => og.toLowerCase().includes(q) || genes.some(g => g.toLowerCase().includes(q))) : D.trees;
  i = Math.max(0, shown.findIndex(t => t[0] === cur));
  if (shown[i] && shown[i][0] === cur) { drawTree(false); $("count").textContent = `${i + 1} of ${shown.length.toLocaleString()}`; }
  else show();
});
$("og").addEventListener("change", e => { i = +e.target.value; show(); });
$("prev").onclick = () => { i--; show(); };
$("next").onclick = () => { i++; show(); };
// ---- the two tabs: the gene trees, and the species tree explorer ----
function openMainTab(id) {
  for (const t of ["trees", "explore", "allgenes"]) {
    $("tab-" + t).hidden = t !== id;
    document.querySelector(`#maintabs [data-tab="${t}"]`).setAttribute("aria-selected", t === id);
  }
  if (id === "explore") renderExplorer();
  else if (id === "allgenes") renderAllGenes();
  else if (G) drawTree(false);  // drawn for the width it now has
}
$("maintabs").onclick = e => { const b = e.target.closest("[data-tab]"); if (b) openMainTab(b.dataset.tab); };
function showGeneTree(og) {  // the Gene trees tab, at this orthogroup
  $("q").value = ""; shown = D.trees; i = Math.max(0, D.trees.findIndex(t => t[0] === og));
  openMainTab("trees"); show(); window.scrollTo(0, $("tab-trees").offsetTop);
}

// ---- data files, each loaded once, when first needed ----
const cache = {}, waiting = {};
function fetchData(key, src) {
  if (key in cache) return Promise.resolve(cache[key]);
  return new Promise((resolve, reject) => {
    if (!waiting[key]) {
      waiting[key] = [];
      const s = document.createElement("script");
      s.src = src;
      s.onerror = () => { (waiting[key] || []).forEach(w => w.reject(new Error(`could not be loaded (${src})`))); delete waiting[key]; };
      document.head.appendChild(s);
    }
    waiting[key].push({resolve, reject});
  });
}
function gotData(key, value) { cache[key] = value; (waiting[key] || []).forEach(w => w.resolve(value)); delete waiting[key]; }
window.hogsLoaded = d => gotData("hogs:" + d.node, d.hogs);
window.speciesLoaded = d => gotData("species:" + d.k, d);
const TREE_OGS = new Set(D.trees.map(t => t[0]));
const loadTree = og => loaded[og] ? Promise.resolve(loaded[og]) : fetchData("tree:" + og, `${D.data}/${encodeURIComponent(og)}.js`);
const loadHogs = node => fetchData("hogs:" + node, `${D.data}/HOGs_${encodeURIComponent(node)}.js`);
const loadSpecies = k => fetchData("species:" + k, `${D.data}/Species_${k}.js`);

// ---- the species tree ----
const ST = (() => {
  if (!D.species_tree) return null;
  const T = D.species_tree, n = T.nodes.length, kids = Array.from({length: n}, () => []), level = new Int32Array(n);
  for (let k = 1; k < n; k++) { kids[T.nodes[k][0]].push(k); level[k] = level[T.nodes[k][0]] + 1; }
  const isLeaf = k => !kids[k].length, sp = k => isLeaf(k) ? T.genes[T.nodes[k][2]][1] : -1;
  const name = k => isLeaf(k) ? T.genes[T.nodes[k][2]][0] : T.nodes[k][2];
  const leafOf = {}, byName = {};
  for (let k = 0; k < n; k++) { if (isLeaf(k) && sp(k) >= 0) leafOf[sp(k)] = k; byName[name(k)] = k; }
  const parent = k => T.nodes[k][0];
  function lca(a, b) { while (level[a] > level[b]) a = parent(a); while (level[b] > level[a]) b = parent(b); while (a !== b) { a = parent(a); b = parent(b); } return a; }
  function speciesBelow(k) { const out = new Set(), st = [k]; while (st.length) { const j = st.pop(); if (isLeaf(j)) out.add(sp(j)); else st.push(...kids[j]); } return out; }
  return {T, n, kids, level, isLeaf, sp, name, leafOf, byName, parent, lca, speciesBelow};
})();
const HOG_COUNT = Object.fromEntries(D.hog_nodes || []);
const geneName = g => g[1] >= 0 && g[0].startsWith(D.species[g[1]][0] + ": ") ? g[0].slice(D.species[g[1]][0].length + 2) : g[0];
const geneKey = (sp, name) => sp + "\t" + name;
const spChip = (k, text) => `<span class="chip" style="--c:${D.species[k][1]}">${esc(text)}</span>`;

// a gene tree cut down to some of its genes (keep(gene)), below one of its nodes (or all of it)
function pruneTree(T, keep, cladeName) {
  const n = T.nodes.length, kids = Array.from({length: n}, () => []);
  for (let k = 1; k < n; k++) kids[T.nodes[k][0]].push(k);
  let root = 0;
  if (cladeName) { const r = T.nodes.findIndex(nd => typeof nd[2] === "string" && nd[2] === cladeName); if (r >= 0) root = r; }
  const has = new Uint8Array(n);
  for (let k = n - 1; k >= 0; k--) has[k] = kids[k].length ? +kids[k].some(c => has[c]) : +!!keep(T.genes[T.nodes[k][2]]);
  const nodes = [], genes = [], stack = [[root, -1]];
  while (stack.length) {
    let [k, np] = stack.pop(), d = np < 0 ? 0 : T.nodes[k][1];
    if (!has[k]) continue;
    let live = kids[k].filter(c => has[c]);
    while (live.length === 1) { k = live[0]; d += T.nodes[k][1]; live = kids[k].filter(c => has[c]); }  // no nodes with one child
    if (np < 0) d = 0;
    const idx = nodes.length;
    if (!kids[k].length) { nodes.push([np, d, genes.length]); genes.push(T.genes[T.nodes[k][2]]); }
    else { nodes.push([np, d, T.nodes[k][2]]); for (let j = live.length - 1; j >= 0; j--) stack.push([live[j], idx]); }
  }
  return genes.length ? {og: T.og, genes, nodes, dups: T.dups} : null;
}
// a small gene tree drawn in the explorer; a click on a node collapses or expands it
const smallTrees = [];
function smallTreeHtml(T) {
  const id = smallTrees.length;
  smallTrees.push(prepareTree(T));
  return `<div class="clade" data-gid="${id}"></div>`;
}
function drawSmallTrees(root) {
  for (const el of root.querySelectorAll(".clade[data-gid]:not([data-drawn])")) {
    el.dataset.drawn = 1;
    el.innerHTML = treeSvg(smallTrees[+el.dataset.gid], {width: Math.max(360, el.clientWidth), zoom: 2, clickable: true}).svg;
  }
}
$("tab-explore").addEventListener("click", e => {
  const box = e.target.closest(".clade[data-gid]"), k = e.target.closest("[data-k]");
  if (box && k) {
    const T = smallTrees[+box.dataset.gid], j = +k.dataset.k;
    if (T.collapsed.has(j)) T.collapsed.delete(j); else T.collapsed.add(j);
    box.innerHTML = treeSvg(T, {width: Math.max(360, box.clientWidth), zoom: 2, clickable: true}).svg;
    return;
  }
  const og = e.target.closest("[data-open-og]");
  if (og) { e.preventDefault(); showGeneTree(og.dataset.openOg); }
});
const openOgLink = og => TREE_OGS.has(og) ? `<a href="#" data-open-og="${esc(og)}">gene tree and alignment</a>` : "";

// <details> whose content is made when first opened
function lazy(summaryHtml, fill, {open = false, cls = ""} = {}) {
  const d = document.createElement("details");
  if (cls) d.className = cls;
  d.innerHTML = `<summary>${summaryHtml}</summary><div class="body"></div>`;
  let filled = false;
  const doFill = () => { if (d.open && !filled) { filled = true; fill(d.querySelector(".body")); } };
  d.addEventListener("toggle", doFill);
  if (open) { d.open = true; setTimeout(doFill); }
  return d;
}
// a long list shown a part at a time, with a search
function longList(el, items, render, {placeholder, match, step = 30}) {
  el.innerHTML = `<div class="tools"><input placeholder="${esc(placeholder)}" autocomplete="off"><span class="muted"></span></div><div class="items"></div>`
    + `<button class="more" hidden>Show more</button>`;
  const input = el.querySelector("input"), list = el.querySelector(".items"), more = el.querySelector(".more"), count = el.querySelector(".muted");
  let shownItems = items, n = 0;
  const add = () => {
    const end = Math.min(shownItems.length, n + step);
    for (; n < end; n++) list.appendChild(render(shownItems[n]));
    more.hidden = n >= shownItems.length;
    count.textContent = `${n.toLocaleString()} of ${shownItems.length.toLocaleString()}`;
  };
  input.oninput = () => { const q = input.value.trim().toLowerCase(); shownItems = q ? items.filter(x => match(x, q)) : items; list.innerHTML = ""; n = 0; add(); };
  more.onclick = add;
  add();
}

// ---- the explorer ----
const MODES = {
  hogs: "Open a node of the species tree to see its hierarchical orthogroups (HOGs); open a HOG to see its genes "
    + "(its clade of the gene tree, coloured by species) and the HOGs it contains at the nodes below.",
  og: "All orthogroups, or one, in the species tree: the width of each branch shows how many gene copies it carries, "
    + "■ the gene duplications at each node (Duplications.tsv) and ✕ the gene losses (from mapping the gene tree onto the species tree). "
    + "Click a node or a species to see its genes.",
  species: "Open a species to see its genes, in each orthogroup and unassigned.",
};
let exMode = "hogs", exDrawn = null;
function renderExplorer() {
  const mode = document.querySelector('input[name="exmode"]:checked').value;
  if (exDrawn === mode) return;
  exDrawn = exMode = mode;
  $("ex-help").textContent = MODES[mode];
  const view = $("ex-view");
  view.innerHTML = "";
  if (!ST) { view.innerHTML = `<p class="muted">There is no species tree in these results (Species_Tree/SpeciesTree_rooted_node_labels.txt).</p>`; return; }
  if (mode === "og") return renderOgMode(view);
  const outline = document.createElement("div");
  outline.className = "outline";
  outline.appendChild(speciesTreeNode(0, mode, true));
  view.appendChild(outline);
}
$("ex-modes").addEventListener("change", renderExplorer);

function speciesTreeNode(k, mode, open = false) {
  if (ST.isLeaf(k)) {
    const sp = ST.sp(k), n = sp >= 0 ? (D.species_genes || [])[sp] : undefined;
    const label = (sp >= 0 ? spChip(sp, ST.name(k)) : esc(ST.name(k))) + (n !== undefined ? ` <span class="muted">${num(n)} genes</span>` : "");
    if (mode !== "species" || sp < 0) { const div = document.createElement("div"); div.className = "leaf"; div.innerHTML = label; return div; }
    return lazy(label, body => fillSpecies(body, sp));
  }
  const name = ST.name(k), nSp = ST.speciesBelow(k).size;
  const label = `<b>${esc(name)}</b> <span class="muted">${nSp} species` + (mode === "hogs" ? `; ${HOG_COUNT[name] !== undefined ? num(HOG_COUNT[name]) + " HOGs" : "no HOG file"}` : "") + `</span>`;
  return lazy(label, body => {
    if (mode === "hogs" && HOG_COUNT[name] !== undefined)
      body.appendChild(lazy(`HOGs at ${esc(name)} <span class="muted">(${num(HOG_COUNT[name])})</span>`, b => fillHogs(b, name), {cls: "hoglist"}));
    for (const c of ST.kids[k]) body.appendChild(speciesTreeNode(c, mode));
  }, {open});
}

function hogSummary(h) {
  const n = h[3].reduce((s, [, g]) => s + g.length, 0);
  return `<b>${esc(h[0])}</b> <span class="muted">${esc(h[1])} · ${num(n)} gene${n === 1 ? "" : "s"}, ${h[3].length} species</span> `
    + h[3].map(([k, g]) => spChip(k, String(g.length))).join("");
}
function fillHogs(body, node) {
  body.innerHTML = `<p class="muted">Loading the HOGs of ${esc(node)}...</p>`;
  loadHogs(node).then(hogs => longList(body, hogs, h => hogDetails(h, node), {
    placeholder: "Search a HOG, orthogroup or gene", match: (h, q) => h[0].toLowerCase().includes(q) || h[1].toLowerCase().includes(q)
      || h[3].some(([, g]) => g.some(x => x.toLowerCase().includes(q)))}),
    err => { body.innerHTML = `<p class="muted">The HOGs of ${esc(node)} ${esc(err.message)}.</p>`; });
}
function hogDetails(h, node) {
  return lazy(hogSummary(h), body => {
    const genes = new Set(h[3].flatMap(([k, g]) => g.map(x => geneKey(k, x))));
    body.innerHTML = `<div class="muted">${esc(h[0])}: the genes of ${esc(h[1])} below the node ${esc(h[2])} of its gene tree. ${openOgLink(h[1])}</div><div class="hogtree"></div><div class="below"></div>`;
    const treeEl = body.querySelector(".hogtree");
    const listGenes = () => { treeEl.innerHTML = h[3].map(([k, g]) => `<div>${spChip(k, D.species[k][0])}${g.map(esc).join(", ")}</div>`).join(""); };
    if (TREE_OGS.has(h[1])) {
      treeEl.innerHTML = `<p class="muted">Loading the gene tree of ${esc(h[1])}...</p>`;
      loadTree(h[1]).then(T => {
        const P = pruneTree(T, g => genes.has(geneKey(g[1], geneName(g))), h[2]);
        if (!P) return listGenes();
        treeEl.innerHTML = smallTreeHtml(P);
        drawSmallTrees(treeEl);
      }, listGenes);
    } else listGenes();
    // the HOGs it contains at the nodes below, and its genes in the species below
    const below = body.querySelector(".below"), k = ST.byName[node];
    if (k === undefined) return;
    for (const c of ST.kids[k]) {
      if (ST.isLeaf(c)) {
        const sp = ST.sp(c), g = (h[3].find(([s]) => s === sp) || [, []])[1];
        below.insertAdjacentHTML("beforeend", `<div class="leaf">${sp >= 0 ? spChip(sp, ST.name(c)) : esc(ST.name(c)) + ": "}`
          + (g.length ? g.map(esc).join(", ") : `<span class="muted">no genes</span>`) + `</div>`);
        continue;
      }
      const cname = ST.name(c), holder = document.createElement("div");
      holder.className = "sub";
      below.appendChild(holder);
      if (HOG_COUNT[cname] === undefined) { holder.innerHTML = `<span class="muted">${esc(cname)}: no HOG file</span>`; continue; }
      holder.innerHTML = `<span class="muted">Loading the HOGs at ${esc(cname)}...</span>`;
      loadHogs(cname).then(hogs => {
        const inside = hogs.filter(x => x[1] === h[1] && x[3].every(([s, g]) => g.every(y => genes.has(geneKey(s, y)))));
        holder.innerHTML = `<div class="muted">At ${esc(cname)}: ${inside.length ? `${inside.length} HOG${inside.length === 1 ? "" : "s"}` : "no HOG (no genes of its species, or lost)"}</div>`;
        for (const x of inside) holder.appendChild(hogDetails(x, cname));
      }, err => { holder.innerHTML = `<span class="muted">The HOGs at ${esc(cname)} ${esc(err.message)}.</span>`; });
    }
  }, {cls: "hog"});
}

function fillSpecies(body, sp) {
  body.innerHTML = `<p class="muted">Loading the genes of ${esc(D.species[sp][0])}...</p>`;
  loadSpecies(sp).then(d => {
    body.innerHTML = "";
    if (d.unassigned.length)
      body.appendChild(lazy(`Unassigned genes <span class="muted">(${num(d.unassigned.length)}, in no orthogroup)</span>`,
        b => longList(b, d.unassigned, g => { const x = document.createElement("div"); x.className = "gene"; x.textContent = g; return x; },
                      {placeholder: "Search a gene", match: (g, q) => g.toLowerCase().includes(q), step: 100})));
    const list = document.createElement("div");
    body.appendChild(list);
    const ogs = [...d.orthogroups].sort((a, b) => b[1].length - a[1].length || (a[0] < b[0] ? -1 : 1));
    longList(list, ogs, ([og, genes]) => lazy(`<b>${esc(og)}</b> <span class="muted">${num(genes.length)} gene${genes.length === 1 ? "" : "s"}</span>`,
      b => { b.innerHTML = `<div>${genes.map(esc).join(", ")}</div><div>${openOgLink(og)}</div>`; }),
      {placeholder: "Search an orthogroup or gene", match: ([og, genes], q) => og.toLowerCase().includes(q) || genes.some(g => g.toLowerCase().includes(q))});
  }, err => { body.innerHTML = `<p class="muted">The genes of ${esc(D.species[sp][0])} ${esc(err.message)}.</p>`; });
}

// An orthogroup in the species tree: its gene tree mapped onto the species tree
// (each gene tree node to the last common ancestor of its species).
function reconcile(T) {
  const n = T.nodes.length, kids = Array.from({length: n}, () => []), M = new Int32Array(n).fill(-1);
  for (let k = 1; k < n; k++) kids[T.nodes[k][0]].push(k);
  for (let k = n - 1; k >= 0; k--) {
    if (!kids[k].length) { const s = ST.leafOf[T.genes[T.nodes[k][2]][1]]; M[k] = s === undefined ? -1 : s; continue; }
    for (const c of kids[k]) if (M[c] >= 0) M[k] = M[k] < 0 ? M[c] : ST.lca(M[k], M[c]);
  }
  const copies = new Int32Array(ST.n), losses = new Int32Array(ST.n), at = Array.from({length: ST.n}, () => []), genesOf = {};
  for (const nd of T.nodes) if (typeof nd[2] === "number") { const sp = T.genes[nd[2]][1]; genesOf[sp] = (genesOf[sp] || 0) + 1; }
  const isDup = k => kids[k].some(c => M[c] === M[k]);
  for (let k = 0; k < n; k++) {
    if (M[k] < 0) continue;
    at[M[k]].push(k);
    const p = T.nodes[k][0];
    if (p < 0) { for (let y = M[k]; y >= 0; y = ST.parent(y)) copies[y]++; continue; }  // the lineage above the gene tree's root
    if (M[p] < 0 || M[p] === M[k]) continue;
    let y = M[k];
    while (true) {
      copies[y]++;
      const x = ST.parent(y);
      if (x === M[p]) { if (isDup(p)) for (const z of ST.kids[x]) if (z !== y) losses[z]++; break; }
      for (const z of ST.kids[x]) if (z !== y) losses[z]++;  // the lineage passes x without a copy in z
      y = x;
    }
  }
  const dupsAt = Array.from({length: ST.n}, () => []);
  for (const [node, support, spNode] of T.dups || []) if (ST.byName[spNode] !== undefined) dupsAt[ST.byName[spNode]].push([node, support]);
  return {kids, M, copies, losses, at, dupsAt, genesOf};
}
function renderOgMode(view) {
  view.innerHTML = `<div class="tools"><input id="ex-og-q" placeholder="Search an orthogroup ID or gene name" autocomplete="off"><select id="ex-og"></select></div>`
    + `<div id="ex-og-sum" class="muted"></div><div class="tree" id="ex-og-tree"></div><div id="ex-og-node"></div>`;
  const fillSelect = (list, all) => { $("ex-og").innerHTML = (all && D.all_orthogroups ? `<option value="${ALL}">All orthogroups</option>` : "")
    + list.slice(0, 2000).map(t => `<option>${esc(t[0])}</option>`).join(""); };
  fillSelect(D.trees, true);
  $("ex-og-q").oninput = e => {
    const q = e.target.value.trim().toLowerCase();
    fillSelect(q ? D.trees.filter(([og, genes]) => og.toLowerCase().includes(q) || genes.some(g => g.toLowerCase().includes(q))) : D.trees, !q);
    drawOg();
  };
  $("ex-og").onchange = drawOg;
  drawOg();
}
let ogShown = null;
const ALL = "*all*";  // "All orthogroups" in the choice (no orthogroup ID has a "*")
function drawOg() {
  const og = $("ex-og").value;
  if (og === ALL) return drawAll();
  if (!og) { $("ex-og-tree").innerHTML = ""; $("ex-og-sum").textContent = "No orthogroup matches."; return; }
  ogShown = og;
  $("ex-og-tree").innerHTML = `<p class="muted">Loading ${esc(og)}...</p>`;
  $("ex-og-node").innerHTML = "";
  loadTree(og).then(T => {
    if (ogShown !== og) return;
    const R = reconcile(T), nGenes = T.nodes.filter(nd => typeof nd[2] === "number").length;
    const nSp = new Set(T.genes.slice(0, nGenes).map(g => g[1])).size;
    $("ex-og-sum").innerHTML = `${esc(og)}: ${num(nGenes)} genes in ${nSp} of ${D.species.length} species; `
      + `${num(R.dupsAt.reduce((s, a) => s + a.length, 0))} gene duplications, ${num(R.losses.reduce((s, x) => s + x, 0))} gene losses. ${openOgLink(og)}`;
    $("ex-og-tree").innerHTML = reconciledSvg(T, R);
    $("ex-og-tree").onclick = e => { const g = e.target.closest("[data-s]"); if (g) showSpeciesNode(T, R, +g.dataset.s); };
  }, err => { $("ex-og-tree").innerHTML = `<p class="muted">${esc(og)} ${esc(err.message)}.</p>`; });
}
function speciesTreeLayout(W, labelExtra, left = 30) {  // left: room for the root's branch and label  // the species tree drawn left to right, a row for each species
  const n = ST.n, kids = ST.kids, depth = new Float64Array(n), unit = ST.T.nodes.every(([p, d]) => p < 0 || !d);
  for (let k = 1; k < n; k++) depth[k] = depth[ST.parent(k)] + (unit ? 1 : ST.T.nodes[k][1]);
  const y = new Float64Array(n), order = [], st = [0];
  let row = 0;
  while (st.length) { const k = st.pop(); order.push(k); if (kids[k].length) for (let j = kids[k].length - 1; j >= 0; j--) st.push(kids[k][j]); else y[k] = row++; }
  for (let j = order.length - 1; j >= 0; j--) { const k = order[j]; if (kids[k].length) y[k] = (y[kids[k][0]] + y[kids[k][kids[k].length - 1]]) / 2; }
  const rowH = 44, top = 26, labelW = Math.min(300, Math.max(...ST.T.genes.map(g => textWidth(g[0])))) + labelExtra;
  const maxD = Math.max(...depth) || 1, sx = (W - labelW - left - 20) / maxD;
  return {order, kids, height: top + row * rowH, X: k => left + depth[k] * sx, Y: k => top + y[k] * rowH};
}
function reconciledSvg(T, R) {
  const W = Math.max(420, $("ex-og-tree").clientWidth - 20), {order, kids, height, X, Y} = speciesTreeLayout(W, 110);
  const width = c => c ? 1.5 + 1.4 * Math.min(c, 8) : 1, out = [];
  for (const k of order) {
    const p = ST.parent(k), c = R.copies[k], col = ST.isLeaf(k) && ST.sp(k) >= 0 ? spColour(ST.sp(k)) : "var(--ink-555)";
    const x0 = p >= 0 ? X(p) : 10;
    out.push(`<line x1="${x0}" x2="${X(k)}" y1="${Y(k)}" y2="${Y(k)}" stroke="${c ? col : "var(--ink-c8c8c8)"}" stroke-width="${width(c)}" stroke-linecap="round"/>`);
    if (kids[k].length) out.push(`<line x1="${X(k)}" x2="${X(k)}" y1="${Y(kids[k][0])}" y2="${Y(kids[k][kids[k].length - 1])}" stroke="var(--ink-555)" stroke-width="1.2"/>`);
    const mid = (x0 + X(k)) / 2;
    if (c) out.push(`<text x="${mid}" y="${Y(k) - 6}" text-anchor="middle" font-size="11" fill="var(--ink-333)">×${c}</text>`);
    if (R.losses[k]) out.push(`<text x="${mid}" y="${Y(k) + 15}" text-anchor="middle" font-size="11" fill="${DUP}">✕ ${R.losses[k]} lost</text>`);
  }
  for (const k of order) {
    const x = X(k), yk = Y(k), d = R.dupsAt[k], leaf = ST.isLeaf(k), nAt = R.at[k].length;
    const tip = `${ST.name(k)}: ${R.copies[k]} gene cop${R.copies[k] === 1 ? "y" : "ies"} on the branch above it` + (d.length ? `, ${d.length} duplication${d.length === 1 ? "" : "s"} here` : "")
      + (leaf ? `, ${R.genesOf[ST.sp(k)] || 0} genes` : "")
      + (R.losses[k] ? `, ${R.losses[k]} lost` : "") + ". Click for its genes.";
    out.push(`<g data-s="${k}" style="cursor:pointer"><title>${esc(tip)}</title><circle cx="${x}" cy="${yk}" r="12" fill="transparent"/>`
      + (d.length ? `<rect x="${x - 5}" y="${yk - 5}" width="10" height="10" fill="${d.some(z => z[1] >= 0.5) ? DUP : "var(--ink-fff)"}" stroke="${DUP}" stroke-width="1.5"/>`
                  : `<circle cx="${x}" cy="${yk}" r="4" fill="${nAt || R.copies[k] ? NODE : "var(--ink-c8c8c8)"}"/>`) + `</g>`);
    if (d.length) out.push(`<text x="${x - 8}" y="${yk - 7}" text-anchor="end" font-size="11" fill="${DUP}">${d.length}</text>`);
    if (leaf) {
      const sp = ST.sp(k), g = R.genesOf[sp] || 0;
      out.push(`<g data-s="${k}" style="cursor:pointer"><text x="${x + 10}" y="${yk}" dy="0.35em" font-size="12" fill="${sp >= 0 ? spColour(sp) : "var(--ink-333)"}">${esc(ST.name(k))}`
        + ` <tspan fill="var(--ink-666)">(${g} gene${g === 1 ? "" : "s"})</tspan></text></g>`);
    } else out.push(`<text x="${x - 6}" y="${yk + 15}" text-anchor="end" font-size="10" fill="var(--ink-777)">${esc(ST.name(k))}</text>`);
  }
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${height}"><rect width="100%" height="100%" fill="var(--ink-fff)"/>${out.join("")}</svg>`;
}
// All orthogroups: their gene trees mapped onto the species tree (when the page was written), added up
function drawAll() {
  const A = D.all_orthogroups;
  ogShown = ALL;
  $("ex-og-node").innerHTML = "";
  const nDups = A.dups.reduce((s, d) => s + d[0], 0), nLost = A.losses.reduce((s, x) => s + x, 0);
  $("ex-og-sum").innerHTML = `All orthogroups: the ${num(A.n_trees)} with a resolved gene tree`
    + (A.n_orthogroups !== null ? ` (of ${num(A.n_orthogroups)}; gene trees are made for the orthogroups with at least 4 genes)` : "")
    + `, mapped onto the species tree and added up: ${num(nDups)} gene duplications, ${num(nLost)} gene losses. `
    + `The width of each branch shows the gene copies on it (in these orthogroups); + the orthogroups gained there (all their genes are of the species below it). Click a node or a species for its orthogroups.`;
  const branchLabel = k => `${num(A.copies[k])} genes` + (A.gains[k] ? `<tspan fill="${NODE}" font-weight="600">  +${num(A.gains[k])} gained</tspan>` : "");
  const rootW = textWidth(`${num(A.copies[0])} genes  +${num(A.gains[0])} gained`) + 34;
  const W = Math.max(420, $("ex-og-tree").clientWidth - 20), {order, kids, height, X, Y} = speciesTreeLayout(W, 130, rootW);
  const maxC = Math.max(1, ...A.copies), width = c => c ? 1.5 + 9 * Math.sqrt(c / maxC) : 1, out = [];
  for (const k of order) {
    const p = ST.parent(k), c = A.copies[k], col = ST.isLeaf(k) && ST.sp(k) >= 0 ? spColour(ST.sp(k)) : "var(--ink-555)";
    const x0 = p >= 0 ? X(p) : X(k) - 16, mid = (x0 + X(k)) / 2;
    out.push(`<line x1="${x0}" x2="${X(k)}" y1="${Y(k)}" y2="${Y(k)}" stroke="${c ? col : "var(--ink-c8c8c8)"}" stroke-width="${width(c)}"/>`);
    if (kids[k].length) out.push(`<line x1="${X(k)}" x2="${X(k)}" y1="${Y(kids[k][0])}" y2="${Y(kids[k][kids[k].length - 1])}" stroke="var(--ink-555)" stroke-width="1.2"/>`);
    out.push(p >= 0 ? `<text x="${mid}" y="${Y(k) - 8}" text-anchor="middle" font-size="11" fill="var(--ink-333)">${branchLabel(k)}</text>`
                    : `<text x="${x0 - 4}" y="${Y(k)}" dy="0.35em" text-anchor="end" font-size="11" fill="var(--ink-333)">${branchLabel(k)}</text>`);
    if (A.losses[k]) out.push(`<text x="${mid}" y="${Y(k) + 17}" text-anchor="middle" font-size="11" fill="${DUP}">\u2715 ${num(A.losses[k])} lost</text>`);
  }
  for (const k of order) {
    const x = X(k), yk = Y(k), d = A.dups[k], leaf = ST.isLeaf(k);
    const tip = `${ST.name(k)}: ${num(A.copies[k])} gene copies on the branch above it, ${num(A.gains[k])} orthogroups gained there, `
      + `${num(A.losses[k])} genes lost; ${num(d[0])} gene duplications here (${num(d[1])} with support \u2265 0.5). Click for its orthogroups.`;
    out.push(`<g data-s="${k}" style="cursor:pointer"><title>${esc(tip)}</title><circle cx="${x}" cy="${yk}" r="12" fill="transparent"/>`
      + (d[0] ? `<rect x="${x - 5}" y="${yk - 5}" width="10" height="10" fill="${d[1] ? DUP : "var(--ink-fff)"}" stroke="${DUP}" stroke-width="1.5"/>` : `<circle cx="${x}" cy="${yk}" r="4" fill="${NODE}"/>`) + `</g>`);
    if (d[0]) out.push(`<text x="${x - 8}" y="${yk - 7}" text-anchor="end" font-size="11" fill="${DUP}">${num(d[0])}</text>`);
    if (leaf) {
      const sp = ST.sp(k);
      out.push(`<g data-s="${k}" style="cursor:pointer"><text x="${x + 10}" y="${yk}" dy="0.35em" font-size="12" fill="${sp >= 0 ? spColour(sp) : "var(--ink-333)"}">${esc(ST.name(k))}</text></g>`);
    } else out.push(`<text x="${x - 6}" y="${yk + 15}" text-anchor="end" font-size="10" fill="var(--ink-777)">${esc(ST.name(k))}</text>`);
  }
  $("ex-og-tree").innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${height}"><rect width="100%" height="100%" fill="var(--ink-fff)"/>${out.join("")}</svg>`;
  $("ex-og-tree").onclick = e => { const g = e.target.closest("[data-s]"); if (g) showAllNode(+g.dataset.s); };
}
window.allOrthogroupsLoaded = d => gotData("all", d);
function showAllNode(k) {
  const A = D.all_orthogroups, el = $("ex-og-node"), d = A.dups[k];
  el.innerHTML = `<h3 style="margin:12px 0 4px">${esc(ST.name(k))}${ST.isLeaf(k) ? "" : ` <span class="muted">(${ST.speciesBelow(k).size} species)</span>`}</h3>`
    + `<div class="muted node-note">${num(A.copies[k])} gene copies on the branch above it; ${num(d[0])} gene duplications here (${num(d[1])} with support \u2265 0.5). Loading its orthogroups...</div>`;
  fetchData("all", `${D.data}/All_Orthogroups.js`).then(L => {
    if (!el.isConnected) return;
    const note = el.querySelector(".node-note");
    note.textContent = note.textContent.replace(" Loading its orthogroups...", "");
    const outline = document.createElement("div");
    outline.className = "outline";
    const section = (title, rows, render) => {
      if (!rows.length) return;
      outline.appendChild(lazy(`${title} <span class="muted">(${num(rows.length)})</span>`, b => longList(b, rows, render,
        {placeholder: "Search an orthogroup", match: (r, q) => (Array.isArray(r) ? r[0] : r).toLowerCase().includes(q)})));
    };
    const ogRow = (og, note) => { const x = document.createElement("div"); x.className = "gene";
      x.innerHTML = `<b>${esc(og)}</b> ${note ? `<span class="muted">${esc(note)}</span> ` : ""}${openOgLink(og)}`; return x; };
    section("Orthogroups gained here", L.gained[k], og => ogRow(og));
    section("Orthogroups with gene duplications here", L.duplicated[k], ([og, c]) => ogRow(og, `${c} duplication${c === 1 ? "" : "s"}`));
    section("Orthogroups with gene losses here", L.lost[k], ([og, c]) => ogRow(og, `${c} lost`));
    if (!outline.children.length) outline.innerHTML = `<p class="muted">No orthogroups gained, duplicated or lost here.</p>`;
    el.appendChild(outline);
  }, err => { el.insertAdjacentHTML("beforeend", `<p class="muted">The orthogroups ${esc(err.message)}.</p>`); });
}
function showSpeciesNode(T, R, k) {
  const below = ST.speciesBelow(k), el = $("ex-og-node"), d = R.dupsAt[k];
  const P = pruneTree(T, g => below.has(g[1]));
  el.innerHTML = `<h3 style="margin:12px 0 4px">${esc(ST.name(k))}${ST.isLeaf(k) ? "" : ` <span class="muted">(${below.size} species)</span>`}</h3>`
    + `<div class="muted">${R.copies[k]} gene cop${R.copies[k] === 1 ? "y" : "ies"} on the branch above it`
    + (d.length ? `; gene duplications here: ${d.map(([node, s]) => `${esc(node)} (support ${s})`).join(", ")}` : "")
    + (R.losses[k] ? `; ${R.losses[k]} lost` : "") + `. ${P ? "Its genes in the gene tree:" : "No genes of its species in this orthogroup."}</div>`
    + (P ? smallTreeHtml(P) : "");
  drawSmallTrees(el);
}

// ---- All genes: the species tree, each orthogroup's gene tree attached where it was gained ----
// A tree built again at each change from what is open: species tree nodes, a group of the
// orthogroups gained at each, the orthogroups (closed, or their gene tree) and the genes.
const agOpen = new Set(), agClosed = new Set();  // groups and orthogroups opened; species and gene tree nodes closed
let agGainedAt = null, agHighlight = "", agWanted = null;
const AG_STEP = 20, AG_ROW = 20;
const ogInfo = Object.fromEntries(D.trees.map(t => [t[0], t]));
function renderAllGenes() {
  if (!ST || !D.all_orthogroups) {
    $("ag-view").innerHTML = `<p class="muted">This needs the species tree and the orthogroups mapped onto it (Species_Tree/SpeciesTree_rooted_node_labels.txt).</p>`;
    return;
  }
  const A = D.all_orthogroups;
  $("ag-help").textContent = `The species tree, with each orthogroup's gene tree attached at the node where the orthogroup was gained `
    + `(all its genes are of the species below it): the ${num(A.n_trees)} orthogroups with a resolved gene tree`
    + (A.n_orthogroups !== null ? ` (of ${num(A.n_orthogroups)}; gene trees are made for the orthogroups with at least 4 genes)` : "")
    + `. Click a group to open its orthogroups, an orthogroup to open its gene tree, a node to close or open it. `
    + `Drawn without branch lengths: those of the species tree and the gene trees are not on the same scale.`;
  if (agGainedAt) return drawAllGenes();
  $("ag-tree").innerHTML = `<p class="muted">Loading the orthogroups...</p>`;
  fetchData("all", `${D.data}/All_Orthogroups.js`).then(L => {
    agGainedAt = {};
    L.gained.forEach((ogs, k) => ogs.forEach(og => { agGainedAt[og] = k; }));
    agGained = L.gained;
    drawAllGenes();
  }, err => { $("ag-tree").innerHTML = `<p class="muted">The orthogroups ${esc(err.message)}.</p>`; });
}
let agGained = null;
function agBuild() {
  const nodes = [];
  const add = (node, parent) => { node.parent = parent; node.kids = []; node.idx = nodes.length; nodes.push(node); if (parent >= 0) nodes[parent].kids.push(node.idx); return node.idx; };
  const addOg = (og, parent) => {
    const t = ogInfo[og], key = "o:" + og, open = agOpen.has(key), T = loaded[og] || cache["tree:" + og];
    const sps = t[2].species.map(([k]) => k);
    const me = add({key, kind: "og", og, open, sp: sps.length === 1 ? sps[0] : -1, label: og, note: `${num(t[1].length)} genes, ${sps.length} species`, species: t[2].species}, parent);
    if (!open) return;
    if (!T) { nodes[me].note += " (loading)"; agLoad(og); return; }
    const GT = T.kids ? T : prepareTree(T), map = [me], dups = GT.dupAt;
    for (let j = 1; j < GT.nodes.length; j++) {
      const leaf = !GT.kids[j].length, key2 = `g:${og}:${j}`, name = GT.nodes[j][2];
      if (map[GT.nodes[j][0]] === undefined) continue;  // below a closed node
      const g = leaf ? GT.genes[name] : null;
      const idx = add({key: key2, kind: leaf ? "gene" : "gnode", label: leaf ? g[0] : name, sp: GT.sp[j],
                       dup: leaf ? undefined : dups.get(name), n: GT.nGenes[j]}, map[GT.nodes[j][0]]);
      if (leaf || !agClosed.has(key2)) map[j] = idx;
    }
    const rootDup = dups.get(GT.nodes[0][2]);
    if (rootDup !== undefined) nodes[me].dup = rootDup;
  };
  const addSpecies = (k, parent) => {
    const key = "s:" + k, leaf = ST.isLeaf(k);
    const me = add({key, kind: "snode", label: ST.name(k), sp: leaf ? ST.sp(k) : -1, leaf}, parent);
    if (agClosed.has(key)) { nodes[me].closed = true; return; }
    const gained = agGained[k] || [];
    if (gained.length) {
      const gkey = "b:" + k, b = add({key: gkey, kind: "group", open: agOpen.has(gkey),
        label: `${num(gained.length)} orthogroup${gained.length === 1 ? "" : "s"} gained ${leaf ? "in" : "at"} ${ST.name(k)}`}, me);
      if (agOpen.has(gkey)) for (const og of gained) addOg(og, b);
    }
    if (!leaf) for (const c of ST.kids[k]) addSpecies(c, me);
  };
  addSpecies(0, -1);
  return nodes;
}
function agLoad(og) {
  if (agWanted && agWanted.has(og)) return;
  (agWanted = agWanted || new Set()).add(og);
  loadTree(og).then(() => { agWanted.delete(og); drawAllGenes(); }, () => { agWanted.delete(og); });
}
function drawAllGenes() {
  // each node on its own row, its children below it, indented: labels never meet
  const nodes = agBuild(), level = new Int32Array(nodes.length), row = new Int32Array(nodes.length);
  const order = [], st = [0];
  while (st.length) { const k = st.pop(); order.push(k); for (let j = nodes[k].kids.length - 1; j >= 0; j--) st.push(nodes[k].kids[j]); }
  order.forEach((k, r) => { row[k] = r; if (k) level[k] = level[nodes[k].parent] + 1; });
  const X = k => 14 + level[k] * AG_STEP, Y = k => 8 + row[k] * AG_ROW + AG_ROW / 2, q = agHighlight.toLowerCase();
  const out = [], dupMark = (x, y, support) =>
    `<rect x="${x - 4}" y="${y - 4}" width="8" height="8" fill="${support >= 0.5 ? DUP : "var(--ink-fff)"}" stroke="${DUP}" stroke-width="1.5"/>`;
  let width = 0, focus = null;
  for (const k of order) {  // the branches: down from the parent, then across to the node
    const nd = nodes[k], p = nd.parent;
    if (p < 0) continue;
    const col = nd.kind === "snode" ? (nd.leaf && nd.sp >= 0 ? spColour(nd.sp) : "var(--ink-555)") : nd.kind === "group" ? "var(--ink-999)" : spColour(nd.sp);
    out.push(`<path d="M${X(p)},${Y(p)}V${Y(k)}H${X(k)}" fill="none" stroke="${col}" stroke-width="${nd.kind === "snode" ? 2.6 : 1.4}"/>`);
  }
  for (const k of order) {
    const nd = nodes[k], x = X(k), yk = Y(k), hit = `<rect x="${x - 9}" y="${yk - AG_ROW / 2}" width="18" height="${AG_ROW}" fill="transparent"/>`;
    let mark = "", parts = [], tip = "", clickable = true;
    if (nd.kind === "snode") {
      mark = `<circle cx="${x}" cy="${yk}" r="4.5" fill="${nd.closed ? "var(--ink-fff)" : (nd.leaf && nd.sp >= 0 ? spColour(nd.sp) : "var(--ink-333)")}" stroke="var(--ink-333)" stroke-width="1.2"/>`;
      parts = [[(nd.closed ? "▸ " : "") + nd.label, nd.leaf && nd.sp >= 0 ? spColour(nd.sp) : "var(--ink-222)", 600], [nd.leaf ? "" : "  species tree node", "var(--ink-999)", 400]];
      tip = `${nd.label}: click to ${nd.closed ? "open" : "close"}`;
    } else if (nd.kind === "group") {
      mark = `<rect x="${x - 4}" y="${yk - 4}" width="8" height="8" rx="2" fill="var(--ink-999)"/>`;
      parts = [[(nd.open ? "▾ " : "▸ ") + nd.label, "var(--ink-444)", 600]];
      tip = `click to ${nd.open ? "close" : "open"}`;
    } else if (nd.kind === "og") {
      mark = nd.dup !== undefined ? dupMark(x, yk, nd.dup) : `<circle cx="${x}" cy="${yk}" r="3.5" fill="${NODE}"/>`;
      parts = [[(nd.open ? "▾ " : "▸ ") + nd.label, spColour(nd.sp), 600], ["  " + nd.note, "var(--ink-777)", 400]];
      tip = `${nd.label}: ${nd.note}; click to ${nd.open ? "close" : "open its gene tree"}`;
      if (q && nd.label.toLowerCase() === q) focus = yk;
    } else if (nd.kind === "gnode") {
      mark = nd.dup !== undefined ? dupMark(x, yk, nd.dup) : `<circle cx="${x}" cy="${yk}" r="3" fill="${NODE}"/>`;
      parts = [[(nd.kids.length ? "" : "▸ ") + (nd.label || "clade"), "var(--ink-888)", 400], [`  ${nd.n} genes` + (nd.dup !== undefined ? `, gene duplication (support ${nd.dup})` : ""), "var(--ink-aaa)", 400]];
      tip = `${nd.label}: ${nd.n} genes; click to ${nd.kids.length ? "close" : "open"}`;
    } else {  // a gene
      mark = `<circle cx="${x}" cy="${yk}" r="2.5" fill="${spColour(nd.sp)}"/>`;
      parts = [[nd.label, spColour(nd.sp), 400]];
      clickable = false;
      if (q && nd.label.toLowerCase().includes(q)) {
        out.push(`<rect x="${x + 5}" y="${yk - AG_ROW / 2 + 1}" width="${textWidth(nd.label) + 4}" height="${AG_ROW - 2}" fill="var(--ink-fff1a8)"/>`);
        focus = focus ?? yk;
      }
    }
    const label = parts.filter(([t]) => t).map(([t, c, w]) => `<tspan fill="${c}" font-weight="${w}">${esc(t)}</tspan>`).join("");
    const text = `<text x="${x + 8}" y="${yk}" dy="0.35em" font-size="12">${label}</text>`;
    out.push(clickable ? `<g data-key="${esc(nd.key)}" style="cursor:pointer"><title>${esc(tip)}</title>${hit}${mark}${text}</g>` : mark + text);
    let w = x + 8 + textWidth(parts.map(([t]) => t).join("")) + 10;
    if (nd.kind === "og" && !nd.open) {
      for (const [s2, c] of nd.species) { out.push(`<rect x="${w}" y="${yk - 4.5}" width="9" height="9" rx="2" fill="${D.species[s2][1]}"><title>${esc(D.species[s2][0])}: ${c}</title></rect>`); w += 12; }
    }
    width = Math.max(width, w + 10);
  }
  const height = 16 + order.length * AG_ROW;
  $("ag-tree").innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" width="${Math.ceil(width)}" height="${height}"><rect width="100%" height="100%" fill="var(--ink-fff)"/>${out.join("")}</svg>`;
  if (focus !== null) { const box = $("ag-tree"); box.scrollTop = Math.max(0, focus - box.clientHeight / 2); }
}
$("ag-tree").addEventListener("click", e => {
  const el = e.target.closest("[data-key]");
  if (!el) return;
  const key = el.dataset.key;
  if (key.startsWith("b:") || key.startsWith("o:")) { if (agOpen.has(key)) agOpen.delete(key); else agOpen.add(key); }
  else if (agClosed.has(key)) agClosed.delete(key); else agClosed.add(key);
  drawAllGenes();
});
$("ag-open").onclick = () => { if (agGained) { agGained.forEach((ogs, k) => { if (ogs.length) agOpen.add("b:" + k); }); drawAllGenes(); } };
$("ag-close").onclick = () => { agOpen.clear(); agClosed.clear(); agHighlight = ""; drawAllGenes(); };
$("ag-find").onsubmit = e => {
  e.preventDefault();
  const q = $("ag-q").value.trim(), ql = q.toLowerCase();
  if (!q || !agGainedAt) return;
  let og = TREE_OGS.has(q) ? q : null;
  if (!og) { const t = D.trees.find(t => t[1].some(g => g.toLowerCase().includes(ql))); og = t && t[0]; }
  if (!og || agGainedAt[og] === undefined) { $("ag-msg").textContent = `Not found in the gene trees: ${q}`; return; }
  $("ag-msg").textContent = TREE_OGS.has(q) ? "" : `In ${og}`;
  for (let k = agGainedAt[og]; k >= 0; k = ST.parent(k)) agClosed.delete("s:" + k);
  agOpen.add("b:" + agGainedAt[og]);
  agHighlight = q;
  if (og !== q) agOpen.add("o:" + og);
  drawAllGenes();
};

const want = decodeURIComponent(location.hash.slice(1));  // Resolved_Gene_Trees.html#OG0000010, #explore or #allgenes
if (want && want !== "explore" && want !== "allgenes") { const k = D.trees.findIndex(t => t[0] === want); if (k >= 0) i = k; }
show();
if (want === "explore" || want === "allgenes") openMainTab(want);
window.addEventListener("hashchange", () => {  // e.g. sent here from another page in OrthoFinder_Results.html
  const w = decodeURIComponent(location.hash.slice(1));
  if (w === "explore" || w === "allgenes") openMainTab(w);
  else if (D.trees.some(t => t[0] === w)) showGeneTree(w);
});
</script></body></html>
"""


def main(args=None):
    """python -m orthofinder.phylo.tree_pages <results directory>: write the page of its trees."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.phylo.tree_pages",
        description="Write Resolved_Gene_Trees/Resolved_Gene_Trees.html for OrthoFinder results: "
                    "the species tree and the resolved gene trees, with each gene tree's alignment, "
                    "to look at in a web browser.")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    results_dir = parser.parse_args(args).results_dir
    if not os.path.isdir(results_dir):
        parser.error(f"not a directory: {results_dir}")
    n_trees = write_tree_pages(results_dir)
    if not n_trees:
        print(f"No resolved gene trees in {results_dir}")
        return 1
    print(f"Gene tree page ({n_trees} trees): "
          f"{os.path.join(results_dir, 'Resolved_Gene_Trees', 'Resolved_Gene_Trees.html')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
