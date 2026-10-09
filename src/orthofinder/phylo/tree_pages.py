"""
--tree-view (a hidden option): a page to look at the species tree and the
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
"""

import csv
import json
import os

from ete4 import Tree
from ete4.parser.newick import NewickError

from .. import __version__

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
        with open(path, newline="") as infile:
            rows = list(csv.reader(infile, delimiter="\t"))
        species = rows[0][1:]
        for row in rows[1:]:
            for sp, cell in zip(species, row[1:]):
                for gene in cell.split(", ") if cell else []:
                    gene_species[gene] = sp
    return species, gene_species


def leaf_species(leaf_name, species_by_length, gene_species):
    """(species, gene name) of a gene tree leaf: "<species>_<gene>", or "<gene>" (-X)."""
    for sp in species_by_length:
        if leaf_name.startswith(sp + "_"):
            return sp, leaf_name[len(sp) + 1:]
    return gene_species.get(leaf_name), leaf_name


def read_alignment(filename):
    """{sequence name: aligned sequence} of a FASTA alignment."""
    seqs, name = {}, None
    with open(filename) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            elif line and name is not None:
                seqs[name].append(line)
    return {name: "".join(parts) for name, parts in seqs.items()}


def read_duplications(filename):
    """{orthogroup: [[gene tree node, support], ...]} from Duplications.tsv ({} if there is none)."""
    dups = {}
    if not os.path.exists(filename):
        return dups
    with open(filename) as infile:
        next(infile, None)  # the header
        for line in infile:
            # only the first fields: the gene lists can be very long
            og, _, node, support = line.split("\t", 4)[:4]
            dups.setdefault(og, []).append([node, float(support)])
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
            "duplications": [len(duplications), sum(1 for _, support in duplications if support >= 0.5)]}


def write_tree_data(filename, data):
    """The data of one gene tree, as a script that gives it to the page."""
    with open(filename, "w") as out:
        out.write(f"treeLoaded({json.dumps(data, separators=(',', ':'))});\n")


def write_tree_pages(results_dir):
    """
    Write the page of the trees of a results directory and the data of each
    gene tree. Returns the number of gene trees (0 if the directory has no
    resolved gene trees).
    """
    gene_trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    species_tree_fn = os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    if not os.path.exists(gene_trees_fn):
        return 0
    with open(gene_trees_fn) as infile:
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
    if n_failed:
        print(f"WARNING: {n_failed} gene trees could not be read")
    if not trees:
        return 0

    species_tree = None
    if os.path.exists(species_tree_fn):
        with open(species_tree_fn) as infile:
            tree = Tree(infile.read().strip(), parser=1)
        nodes, names = flat_tree(tree)
        species_tree = {"genes": [[name, species_number.get(name, -1)] for name in names], "nodes": nodes}
    write_page(os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.html"),
               species, colour, trees, species_tree)
    return len(trees)


def write_page(filename, species, colour, trees, species_tree):
    data = {"version": __version__, "data": DATA_DIR, "species_tree": species_tree,
            "species": [[sp, colour[sp]] for sp in species], "trees": trees}
    with open(filename, "w") as out:
        out.write(HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gene trees</title>
<style>
:root { --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e8; --head:#f5f5f7; --accent:#2f6fde; }
@media (prefers-color-scheme: dark) { :root { --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff; } }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width:1200px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:20px; margin:0 0 4px; } h2 { font-size:15px; margin:24px 0 8px; }
.sub, .muted { color:var(--muted); font-size:12px; }
.legend { display:flex; flex-wrap:wrap; gap:6px 14px; margin:12px 0; font-size:13px; }
.legend span::before { content:""; display:inline-block; width:10px; height:10px; border-radius:2px; margin-right:6px; background:var(--c); }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:8px; }
input, select { font:inherit; padding:7px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
input { width:min(380px,100%); } select { max-width:100%; }
button { font:inherit; padding:6px 12px; border:1px solid var(--line); border-radius:6px; background:var(--head); color:var(--fg); cursor:pointer; }
button:disabled { opacity:.4; cursor:default; }
.tree { border:1px solid var(--line); border-radius:8px; padding:8px; background:#fff; overflow:auto; }
.tree svg { display:block; }
.tree.capped { max-height:60vh; }
.tree [data-k] { cursor:pointer; }
.key { font-size:12px; color:var(--muted); margin:0 0 6px; }
.key svg { vertical-align:-1px; }
.stats { display:grid; grid-template-columns:max-content 1fr; gap:4px 14px; margin:4px 0 12px; font-size:13px; }
.stats dt { color:var(--muted); } .stats dd { margin:0; }
.chip { display:inline-block; margin:0 10px 2px 0; white-space:nowrap; }
.chip::before { content:""; display:inline-block; width:9px; height:9px; border-radius:2px; margin-right:5px; background:var(--c); }
.chip.none { opacity:.45; }
.badge { display:inline-block; padding:0 7px; border-radius:9px; background:var(--head); border:1px solid var(--line); font-size:12px; }
.aln { border:1px solid var(--line); border-radius:8px; background:#fff; overflow:auto; max-height:70vh; }
.aln canvas { position:sticky; top:0; left:0; display:block; }
#alnhover { font-family:ui-monospace, Menlo, Consolas, monospace; }
</style></head><body><main>
<h1>Gene trees</h1>
<div class="sub" id="sub"></div>
<div class="legend" id="legend"></div>
<div id="species"></div>
<h2>Resolved gene trees</h2>
<div class="tools">
  <input id="q" placeholder="Search an orthogroup ID or gene name" autocomplete="off">
  <select id="og"></select>
  <button id="prev">Previous</button><button id="next">Next</button>
  <span class="muted" id="count"></span>
</div>
<dl class="stats" id="stats"></dl>
<div class="tools">
  <button id="tout" title="Narrower">&minus;</button><button id="tin" title="Wider">+</button>
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
</main>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const num = (x, d = 0) => Number(x).toLocaleString(undefined, {maximumFractionDigits: d});
const MIXED = "#555", DUP = "#d62728", NODE = "#2b3a67";
const spColour = s => s >= 0 ? D.species[s][1] : MIXED;
$("sub").textContent = "OrthoFinder " + D.version + ". The trees of Resolved_Gene_Trees.txt, genes coloured by species; "
  + "n0, n1, ... are the gene tree nodes (the HOG files' gene tree parent clade). "
  + "Below each tree, its alignment (MultipleSequenceAlignments).";
$("legend").innerHTML = D.species.map(([s, c]) => `<span style="--c:${c}">${esc(s)}</span>`).join("");
$("key").innerHTML = "Click a node to collapse or expand its clade. A branch is in a species' colour when all its genes are of that species. "
  + `<svg width="10" height="10"><rect x="1" y="1" width="8" height="8" fill="${DUP}"/></svg> gene duplication (support &ge; 0.5), `
  + `<svg width="10" height="10"><rect x="1.5" y="1.5" width="7" height="7" fill="#fff" stroke="${DUP}" stroke-width="1.5"/></svg> support &lt; 0.5.`;

// ---- trees: drawn here, as SVG ----
// A tree is {genes: [[label, species number]], nodes: [[parent, branch length, gene number or node name]]},
// its nodes in preorder (a parent before its children).
const ROWH = 18, FONT = "12px system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif", TREE_ZOOMS = [0.35, 0.5, 0.7, 1, 1.4, 2, 3, 4.5];
const measurer = document.createElement("canvas").getContext("2d");
function textWidth(s, font = FONT) { measurer.font = font; return measurer.measureText(s).width; }

function prepareTree(T) {
  const n = T.nodes.length, kids = Array.from({length: n}, () => []);
  const depth = new Float64Array(n), sp = new Int32Array(n).fill(-2), nGenes = new Int32Array(n), deepest = new Float64Array(n);
  const unit = T.nodes.every(([p, d]) => p < 0 || !d);  // no branch lengths: draw the topology
  for (let k = 0; k < n; k++) {
    const p = T.nodes[k][0];
    if (p >= 0) { kids[p].push(k); depth[k] = depth[p] + (unit ? 1 : T.nodes[k][1]); }
  }
  for (let k = n - 1; k >= 0; k--) {
    if (!kids[k].length) { nGenes[k] = 1; sp[k] = T.genes[T.nodes[k][2]][1]; deepest[k] = depth[k]; }
    const p = T.nodes[k][0];
    if (p < 0) continue;
    nGenes[p] += nGenes[k]; deepest[p] = Math.max(deepest[p], deepest[k]);
    sp[p] = sp[p] === -2 ? sp[k] : (sp[p] === sp[k] ? sp[p] : -1);
  }
  const dups = new Map((T.dups || []).map(([name, support]) => [name, support]));
  return Object.assign(T, {kids, depth, sp, nGenes, deepest, unitLengths: unit, dupAt: dups, collapsed: new Set()});
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
  const L = layoutTree(G), {rows, order, y} = L;
  const labels = rows.map(k => rowLabel(G, k));
  const labelW = Math.max(40, ...labels.map(s => textWidth(s)));
  const maxDepth = G.deepest[0] || 1;
  const treeW = Math.max(120, Math.min(900, opts.width - labelW - 60)) * TREE_ZOOMS[opts.zoom ?? 3];
  const sx = treeW / maxDepth, left = 14, top = 6;
  const X = k => left + G.depth[k] * sx, Y = r => top + r * ROWH + ROWH / 2;
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
      const x2 = left + G.deepest[k] * sx;
      out.push(`<g data-k="${k}"><title>${esc(cladeLabel(G, k))} (click to expand)</title>`
        + `<polygon points="${x},${yk} ${x2},${yk - 6} ${x2},${yk + 6}" fill="${col}" fill-opacity="0.35" stroke="${col}"/></g>`);
      x = x2;
    } else {
      out.push(`<circle cx="${x}" cy="${yk}" r="2.5" fill="${col}"/>`);
    }
    if (q && label.toLowerCase().includes(q))
      out.push(`<rect x="${x + 4}" y="${yk - ROWH / 2 + 1}" width="${textWidth(label) + 4}" height="${ROWH - 2}" fill="#fff1a8"/>`);
    out.push(`<text x="${x + 6}" y="${yk}" dy="0.35em" fill="${col}" font-size="12" font-family="system-ui, sans-serif"`
      + (G.kids[k].length ? ` data-k="${k}"` : "") + `>${esc(label)}</text>`);
  });
  // internal nodes: their names, and a marker (a square for a gene duplication) to click
  for (const k of order) {
    if (!G.kids[k].length || G.collapsed.has(k)) continue;
    const x = X(k), yk = Y(y[k]), name = G.nodes[k][2], support = G.dupAt.get(name);
    if (name) out.push(`<text x="${x - 3}" y="${yk - 4}" text-anchor="end" fill="#777" font-size="9" font-family="system-ui, sans-serif">${esc(name)}</text>`);
    const tip = `${cladeLabel(G, k)}` + (support !== undefined ? `; gene duplication, support ${support}` : "") + (opts.clickable ? " (click to collapse)" : "");
    const mark = support === undefined ? `<circle cx="${x}" cy="${yk}" r="3" fill="${NODE}"/>`
      : `<rect x="${x - 4}" y="${yk - 4}" width="8" height="8" fill="${support >= 0.5 ? DUP : "#fff"}" stroke="${DUP}" stroke-width="1.5"/>`;
    out.push(opts.clickable ? `<g data-k="${k}"><title>${esc(tip)}</title><circle cx="${x}" cy="${yk}" r="8" fill="transparent"/>${mark}</g>`
                            : `<g><title>${esc(tip)}</title>${mark}</g>`);
  }
  // scale bar
  const height = top + rows.length * ROWH + 34, bar = niceLength(maxDepth / 5);
  if (!G.unitLengths) {
    const by = height - 18;
    out.push(`<line x1="${left}" y1="${by}" x2="${left + bar * sx}" y2="${by}" stroke="#333"/>`
      + `<line x1="${left}" y1="${by - 4}" x2="${left}" y2="${by + 4}" stroke="#333"/><line x1="${left + bar * sx}" y1="${by - 4}" x2="${left + bar * sx}" y2="${by + 4}" stroke="#333"/>`
      + `<text x="${left}" y="${by + 14}" font-size="11" font-family="system-ui, sans-serif" fill="#333">${num(bar, 6)}</text>`);
  }
  const width = Math.ceil(left + treeW + 12 + labelW + 16);
  return {svg: `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`
    + `<rect width="100%" height="100%" fill="#fff"/>${out.join("")}</svg>`, rows};
}

if (D.species_tree) {
  const S = prepareTree(D.species_tree);
  $("species").innerHTML = `<details open><summary><b>Species tree</b> <span class="muted">(N0, N1, ... as in the HOG files)</span></summary>`
    + `<div class="tree">${treeSvg(S, {width: Math.min(700, $("species").clientWidth), clickable: false}).svg}</div></details>`;
}

// ---- the gene tree shown ----
const loaded = {};
let shown = D.trees, i = 0, G = null, wanted = null, treeZoom = 3, wholeTree = false;
window.treeLoaded = T => { loaded[T.og] = T; if (T.og === wanted) setTree(T); };
function show() {
  $("og").innerHTML = shown.map((t, k) => `<option value="${k}">${esc(t[0])} (${t[1].length} genes)</option>`).join("");
  $("og").value = String(i);
  const t = shown[i];
  $("count").textContent = shown.length ? `${i + 1} of ${shown.length.toLocaleString()}` : "no matches";
  $("prev").disabled = i <= 0;
  $("next").disabled = i >= shown.length - 1;
  showStats(t);
  if (!t) { wanted = null; G = null; $("treebox").innerHTML = ""; $("alnbox").hidden = true; return; }
  wanted = t[0];
  if (loaded[t[0]]) return setTree(loaded[t[0]]);
  $("treebox").innerHTML = `<p class="muted">Loading ${esc(t[0])}...</p>`;
  const s = document.createElement("script");
  s.src = `${D.data}/${encodeURIComponent(t[0])}.js`;
  s.onerror = () => { if (wanted === t[0]) $("treebox").innerHTML = `<p class="muted">The data of ${esc(t[0])} could not be loaded (${esc(s.src)}).</p>`; };
  document.head.appendChild(s);
}
function setTree(T) {
  G = T.kids ? T : prepareTree(T);
  $("treebox").scrollTo(0, 0);
  drawTree(true);
}
function drawTree(newTree) {
  if (!G) return;
  const box = $("treebox"), {svg, rows} = treeSvg(G, {width: box.clientWidth, zoom: treeZoom, clickable: true, highlight: $("q").value.trim()});
  box.innerHTML = svg;
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
$("expand").onclick = () => { if (G) { G.collapsed.clear(); drawTree(false); } };
$("tin").onclick = () => { treeZoom = Math.min(TREE_ZOOMS.length - 1, treeZoom + 1); drawTree(false); };
$("tout").onclick = () => { treeZoom = Math.max(0, treeZoom - 1); drawTree(false); };
$("save").onclick = () => {
  const svg = $("treebox").querySelector("svg");
  if (!svg || !G) return;
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + svg.outerHTML], {type: "image/svg+xml"}));
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
const ROW = 16, RULER = 20, ZOOMS = [1, 2, 3, 5, 8, 11, 14];
const AA = {}, NT = {};
for (const [cs, col] of [["AILMFWV", "#b7d0ff"], ["KR", "#ffb8b8"], ["DE", "#f3c0ff"], ["NQST", "#bdf0bd"],
                         ["C", "#ffd9e8"], ["G", "#ffd9a6"], ["P", "#fff0a0"], ["HY", "#aef0ee"]])
  for (const c of cs) AA[c] = col;
for (const [c, col] of [["A", "#a8e6a8"], ["C", "#a9cdff"], ["G", "#ffd598"], ["T", "#ffacac"], ["U", "#ffacac"]]) NT[c] = col;
let aln = null, zoom = 4, nameW = 0, scheme = AA, drawPending = false;
function setAlignment(rows, reset) {
  const seqs = rows.map(r => r[2]).filter(s => s);
  aln = {rows, ncol: Math.max(0, ...seqs.map(s => s.length))};
  const residues = seqs.slice(0, 20).join("").replace(/[-.]/g, "").toUpperCase();
  scheme = residues && /^[ACGTUN]*$/.test(residues) ? NT : AA;
  $("alnbox").hidden = false;
  $("alninfo").textContent = `(${seqs.length} sequences, ${aln.ncol.toLocaleString()} columns; rows in the order of the tree)`;
  nameW = Math.min(Math.max(...rows.map(r => textWidth(r[0])), 40) + 16, Math.max(160, $("alnview").clientWidth * 0.45));
  if (reset) $("alnview").scrollTo(0, 0);
  layout();
}
function layout() {
  if (!aln) return;
  const cw = ZOOMS[zoom];
  $("alnspace").style.width = (nameW + aln.ncol * cw) + "px";
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
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h);
  const x0 = view.scrollLeft, y0 = view.scrollTop, cw = ZOOMS[zoom], rows = aln.rows, ncol = aln.ncol;
  const c0 = Math.floor(x0 / cw), c1 = Math.min(ncol, Math.ceil((x0 + w - nameW) / cw) + 1);
  const r0 = Math.floor(y0 / ROW), r1 = Math.min(rows.length, Math.ceil((y0 + h - RULER) / ROW) + 1);
  ctx.textAlign = "center"; ctx.textBaseline = "middle";
  ctx.font = `${Math.min(12, cw + 1)}px ui-monospace, Menlo, Consolas, monospace`;
  for (let r = r0; r < r1; r++) {
    const seq = rows[r][2], y = RULER + r * ROW - y0;
    if (seq === null) {  // a collapsed clade, or a gene of the tree without a sequence
      ctx.fillStyle = "#eeeef2"; ctx.fillRect(nameW + c0 * cw - x0, y + 1, (c1 - c0) * cw, ROW - 2);
      continue;
    }
    for (let c = c0; c < c1; c++) {
      const ch = seq[c];
      if (ch === "-" || ch === "." || ch === undefined) continue;
      const x = nameW + c * cw - x0;
      ctx.fillStyle = scheme[ch.toUpperCase()] || "#e4e4e4";
      ctx.fillRect(x, y + 1, cw, ROW - 2);
      if (cw >= 8) { ctx.fillStyle = "#222"; ctx.fillText(ch, x + cw / 2, y + ROW / 2); }
    }
  }
  // the ruler (column numbers) and the gene names stay in place when scrolling
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, RULER);
  ctx.strokeStyle = "#ccc"; ctx.beginPath(); ctx.moveTo(0, RULER - 0.5); ctx.lineTo(w, RULER - 0.5); ctx.stroke();
  const step = [10, 20, 50, 100, 200, 500, 1000].find(s => s * cw >= 60) || 1000;
  ctx.fillStyle = "#666"; ctx.font = "10px system-ui, sans-serif"; ctx.textAlign = "left";
  for (let c = Math.max(step, Math.ceil((c0 + 1) / step) * step); c <= c1; c += step) {
    const x = nameW + (c - 1) * cw - x0 + cw / 2;
    ctx.fillRect(x, RULER - 5, 1, 4);
    ctx.fillText(String(c), x + 2, RULER / 2 - 1);
  }
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, nameW, h);
  ctx.font = FONT;
  for (let r = r0; r < r1; r++) {
    const [label, sp] = rows[r], y = RULER + r * ROW - y0;
    ctx.fillStyle = spColour(sp);
    ctx.fillText(fit(ctx, label, nameW - 14), 6, y + ROW / 2);
  }
  ctx.strokeStyle = "#ccc"; ctx.beginPath(); ctx.moveTo(nameW - 4.5, 0); ctx.lineTo(nameW - 4.5, h); ctx.stroke();
}
function redraw() { if (!drawPending) { drawPending = true; requestAnimationFrame(draw); } }
$("alnview").addEventListener("scroll", redraw);
window.addEventListener("resize", redraw);
$("zin").onclick = () => { zoomTo(zoom + 1); };
$("zout").onclick = () => { zoomTo(zoom - 1); };
function zoomTo(z) {
  if (!aln || z < 0 || z >= ZOOMS.length) return;
  const view = $("alnview"), mid = (view.scrollLeft + (view.clientWidth - nameW) / 2) / ZOOMS[zoom];
  zoom = z; layout();
  view.scrollLeft = mid * ZOOMS[zoom] - (view.clientWidth - nameW) / 2;
}
$("cv").addEventListener("mousemove", e => {
  if (!aln) return;
  const view = $("alnview"), cw = ZOOMS[zoom];
  const c = Math.floor((e.offsetX - nameW + view.scrollLeft) / cw), r = Math.floor((e.offsetY - RULER + view.scrollTop) / ROW);
  const row = aln.rows[r];
  if (e.offsetX < nameW || e.offsetY < RULER || !row || !row[2] || c < 0 || c >= row[2].length) {
    $("alnhover").textContent = row && e.offsetY >= RULER ? row[0] : "Point at a residue to see its position";
    return;
  }
  const ch = row[2][c], n = row[2].slice(0, c + 1).replace(/[-.]/g, "").length;
  $("alnhover").textContent = `${row[0]} · column ${c + 1} · ` + (ch === "-" || ch === "." ? "gap" : `residue ${n}: ${ch}`);
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
const want = decodeURIComponent(location.hash.slice(1));  // Resolved_Gene_Trees.html#OG0000010
if (want) { const k = D.trees.findIndex(t => t[0] === want); if (k >= 0) i = k; }
show();
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
