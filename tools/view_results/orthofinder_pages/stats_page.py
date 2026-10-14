# Copied from src/orthofinder/orthology/stats_page.py by tools/view_results/make_copies.py: edit that file, not this one.
"""
--stats-view: a page to explore the statistics of a run in a
web browser, with interactive plots:

    Comparative_Genomics_Statistics/Statistics.html
    Comparative_Genomics_Statistics/Statistics_Data/*.js

Overview (the headline numbers, orthogroup sizes, species per orthogroup),
Species (genes in orthogroups, orthogroups and orthologs shared by each pair of
species), Orthogroups (every orthogroup: genes per species, duplications;
searchable, sortable, filterable), Duplications (on the species tree, per
orthogroup) and HOGs (per species tree node, each HOG).

It is made from the results files only (Comparative_Genomics_Statistics,
Orthogroups.GeneCount.tsv, Phylogenetic_Hierarchical_Orthogroups, the species
tree), so it can be written again for existing results:

    python -m orthofinder.orthology.stats_page <results directory>

The page draws its plots itself (plain JavaScript, nothing to install, no
internet needed). It holds a summary of each orthogroup; the genes of each
species in them, the orthogroups of each species, the species x species
matrices and the HOGs of each species tree node are in data files the page
loads when they are needed (view_pages.write_data), so that it stays small and
quick for thousands of species: a page opened from disk may load scripts, but
not read other files.
"""

import json
import os
import re

from ete4 import Tree

from . import __version__
from .tree_pages import PALETTE, flat_tree
from . import view_pages

STATS_DIR = "Comparative_Genomics_Statistics"
DATA_DIR = "Statistics_Data"
ORTHOLOG_KINDS = ["Totals", "one-to-one", "one-to-many", "many-to-one", "many-to-many"]


def read_rows(filename):
    """The rows of a tab-separated file, split only on tabs (cells may be very long)."""
    with open(filename, encoding="utf-8", errors="replace") as infile:
        return [line.rstrip("\r\n").split("\t") for line in infile]


def number(text):
    """A number from a statistics file ("81.0", "2,733", "'1"), or the text if it is not one."""
    t = text.strip().lstrip("'").replace(",", "")
    try:
        return int(t)
    except ValueError:
        try:
            return float(t)
        except ValueError:
            return text.strip()


def first_block(rows):
    """The rows before the first empty row."""
    block = []
    for row in rows:
        if not any(cell.strip() for cell in row):
            break
        block.append(row)
    return block


def read_overall(filename):
    """[[statistic, value], ...] from Statistics_Overall.tsv (its first block)."""
    return [[row[0], number(row[1])] for row in first_block(read_rows(filename)) if len(row) >= 2]


def read_per_species(filename, species):
    """[[statistic, [value for each species]], ...] from Statistics_PerSpecies.tsv (its first block)."""
    block = first_block(read_rows(filename))
    column = {sp: k for k, sp in enumerate(block[0][1:], 1)}
    return [[row[0], [number(row[column[sp]]) if sp in column and column[sp] < len(row) else None
                      for sp in species]]
            for row in block[1:]]


COUNT = re.compile(r"(?:^|\t)([1-9][0-9]*)(?=[\t\r\n]|$)")  # a cell of a number other than 0


def read_gene_counts(filename):
    """
    The species and [[orthogroup, [species, genes, species, genes, ...]], ...]
    (the species with genes) from Orthogroups.GeneCount.tsv: only its cells
    other than 0 are read, as most are 0 for many species.
    """
    with open(filename, encoding="utf-8", errors="replace") as infile:
        header = infile.readline().rstrip("\r\n").split("\t")
        n = len(header) - 2 if header[-1] == "Total" else len(header) - 1
        rows = []
        for line in infile:
            cut = line.find("\t")
            if cut <= 0:
                continue
            pairs, col, last = [], 0, 0
            for m in COUNT.finditer(line, cut):
                col += line.count("\t", last, m.start(1))
                last = m.start(1)
                if col <= n:
                    pairs += [col - 1, int(m.group(1))]
            rows.append([line[:cut], pairs])
    return header[1:1 + n], rows


def read_matrix(filename, species):
    """A species x species matrix (OrthologsStats_*.tsv), in the order of species (12.0 as 12: smaller for the page)."""
    rows = read_rows(filename)
    column = {sp: k for k, sp in enumerate(rows[0][1:], 1)}
    by_row = {row[0]: row for row in rows[1:] if row and row[0]}
    whole = lambda v: int(v) if isinstance(v, float) and v.is_integer() else v
    return [[whole(number(by_row[a][column[b]])) if a in by_row and b in column else None for b in species]
            for a in species]


def shared_orthogroups(counts, n_species):
    """
    The number of orthogroups with genes of both species, for each pair of
    species: from a bitmask of the orthogroups of each species, so the time
    goes with species x species, not with the species of each orthogroup squared.
    """
    bits = [[] for _ in range(n_species)]  # the orthogroups (row numbers) of each species
    for i, (_, pairs) in enumerate(counts):
        for k in pairs[::2]:
            bits[k].append(i)
    masks = [int.from_bytes(bitmap(b), "little") for b in bits]
    shared = [[0] * n_species for _ in range(n_species)]
    for a in range(n_species):
        for b in range(a, n_species):
            shared[a][b] = shared[b][a] = (masks[a] & masks[b]).bit_count()
    return shared


def bitmap(numbers):
    """The bytes of a bitmask with the bits of numbers set (little-endian)."""
    out = bytearray(max(numbers, default=-1) // 8 + 1)
    for i in numbers:
        out[i >> 3] |= 1 << (i & 7)
    return bytes(out)


def node_number(name):
    m = re.match(r"N(\d+)$", name)
    return int(m.group(1)) if m else 10 ** 9


def read_hogs(filename, species_index):
    """
    The HOGs of a species tree node ([HOG, orthogroup, gene tree parent clade,
    [species, genes, species, genes, ...]]), from Phylogenetic_Hierarchical_
    Orthogroups/<node>.tsv: only its non-empty cells are read.
    """
    hogs = []
    with open(filename, encoding="utf-8", errors="replace") as infile:
        header = infile.readline().rstrip("\r\n").split("\t")
        column = {k: species_index[sp] for k, sp in enumerate(header) if sp in species_index}
        for line in infile:
            fields = line.split("\t", 3)
            if len(fields) < 3:
                continue
            pairs = []
            for k, text in view_pages.sparse_cells(line):
                if k in column and k >= 3:
                    pairs += [column[k], text.count(", ") + 1]
            hogs.append([fields[0], fields[1], fields[2].rstrip("\r\n"), pairs])
    return hogs


def hog_summary(node, hogs):
    sizes = [sum(h[3][1::2]) for h in hogs]
    n_species = [len(h[3]) // 2 for h in hogs]
    return {"node": node, "hogs": len(hogs), "genes": sum(sizes),
            "multi_species": sum(1 for k in n_species if k > 1),
            "sizes": sizes_histogram(sizes)}


SIZE_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, None)]


def sizes_histogram(sizes):
    """[[bin label, count], ...] of sizes in SIZE_BINS."""
    counts = [0] * len(SIZE_BINS)
    for s in sizes:
        for k, (lo, hi) in enumerate(SIZE_BINS):
            if s >= lo and (hi is None or s <= hi):
                counts[k] += 1
                break
    labels = [str(lo) if lo == hi else (f"{lo}-{hi}" if hi else f">{lo - 1}") for lo, hi in SIZE_BINS]
    return [[label, c] for label, c in zip(labels, counts)]


def write_stats_page(results_dir, written=None):
    """
    Write the statistics page of a results directory and its data files.
    written: the pages a run is writing, to link to (view_pages.present;
    None: the pages that are there). Returns the page's filename, or None if
    the results have no orthogroups yet.
    """
    stats_dir = os.path.join(results_dir, STATS_DIR)
    gene_count_fn = view_pages.orthogroups_file(results_dir, ".GeneCount.tsv")   # or the MCL orthogroups
    if gene_count_fn is None:
        return None
    species, counts = read_gene_counts(gene_count_fn)
    species_index = {sp: k for k, sp in enumerate(species)}
    data_dir = os.path.join(stats_dir, DATA_DIR)
    os.makedirs(data_dir, exist_ok=True)
    data = {"version": __version__, "species": [[sp, PALETTE[k % len(PALETTE)]] for k, sp in enumerate(species)],
            "data": DATA_DIR}

    def stats_file(name):
        path = os.path.join(stats_dir, name)
        return path if os.path.exists(path) else None

    if fn := stats_file("Statistics_Overall.tsv"):
        data["overall"] = read_overall(fn)
    if fn := stats_file("Statistics_PerSpecies.tsv"):
        data["per_species"] = read_per_species(fn, species)

    dups = {}
    if fn := stats_file("Duplications_per_Orthogroup.tsv"):
        for row in read_rows(fn)[1:]:
            if len(row) >= 3:
                dups[row[0]] = [number(row[1]), number(row[2])]
    with_trees = set()
    trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    if view_pages.link(results_dir, "trees", written) and os.path.exists(trees_fn):
        with open(trees_fn, encoding="utf-8", errors="replace") as infile:
            with_trees = {line[:line.find(":")] for line in infile if ":" in line}
        data["tree_page"] = view_pages.link(results_dir, "trees", written)  # --tree-view
    data["groups_page"] = view_pages.link(results_dir, "groups", written)  # --groups-view
    data["nav"] = view_pages.nav(results_dir, "stats", written)

    # each orthogroup in columns (its genes, species, duplications, those with support >= 0.5, 1 if
    # it has a tree page), and the genes of each species in it in chunks: Counts_<c>
    ids = [og for og, _ in counts]
    og_dups = [dups.get(og, [0, 0]) for og in ids]
    data["ogs"] = {"ids": view_pages.compact_ids(ids), "genes": [sum(p[1::2]) for _, p in counts],
                   "species": [len(p) // 2 for _, p in counts],
                   "dups": [d[0] for d in og_dups], "dups05": [d[1] for d in og_dups],
                   "tree": [int(og in with_trees) for og in ids]}
    data["max_genes"] = max((max(p[1::2]) for _, p in counts if p), default=1)  # of a species in an orthogroup
    data["counts"] = view_pages.write_chunks(data_dir, "Counts", [p for _, p in counts], len)
    data["og_sizes"] = sizes_histogram(data["ogs"]["genes"])
    of_species = [[] for _ in species]  # the orthogroups of each species: Species_<k>
    for i, (_, pairs) in enumerate(counts):
        for k in pairs[::2]:
            of_species[k].append(i)
    for k, rows in enumerate(of_species):
        view_pages.write_data(data_dir, f"Species_{k}", rows)

    # the species x species matrices: Matrix_<i>
    matrices = [("Orthogroups shared", shared_orthogroups(counts, len(species)))]
    for kind in ORTHOLOG_KINDS:
        if fn := stats_file(f"OrthologsStats_{kind}.tsv"):
            matrices.append(("Orthologs: " + ("all" if kind == "Totals" else kind), read_matrix(fn, species)))
    for i, (_, m) in enumerate(matrices):
        view_pages.write_data(data_dir, f"Matrix_{i}", m)
    data["matrices"] = [name for name, _ in matrices]

    if fn := stats_file("Duplications_per_Species_Tree_Node.tsv"):
        data["node_duplications"] = {row[0]: [number(row[1]), number(row[2])]
                                     for row in read_rows(fn)[1:] if len(row) >= 3}
    species_tree_fn = os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    if os.path.exists(species_tree_fn):
        with open(species_tree_fn, encoding="utf-8", errors="replace") as infile:
            nodes, names = flat_tree(Tree(infile.read().strip(), parser=1))
        data["species_tree"] = {"genes": [[name, species_index.get(name, -1)] for name in names], "nodes": nodes}

    # the HOGs of each species tree node: HOGs_<node> (in columns) and HOGCounts_<node>_<c>
    hog_dir = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    hog_files = sorted((fn for fn in os.listdir(hog_dir) if fn.endswith(".tsv")) if os.path.isdir(hog_dir) else [],
                       key=lambda fn: (node_number(fn[:-4]), fn))
    og_row = {og: i for i, og in enumerate(ids)}
    data["hogs"] = []
    for fn in hog_files:
        node = fn[:-4]
        hogs = read_hogs(os.path.join(hog_dir, fn), species_index)
        data["hogs"].append(hog_summary(node, hogs))
        view_pages.write_data(data_dir, f"HOGs_{node}", {
            "ids": view_pages.compact_ids([h[0] for h in hogs]),
            "ogs": [og_row.get(h[1], h[1]) for h in hogs],  # the orthogroup's row (its ID if not one)
            "clades": [h[2] for h in hogs], "genes": [sum(h[3][1::2]) for h in hogs],
            "species": [len(h[3]) // 2 for h in hogs],
            "max_genes": max((max(h[3][1::2]) for h in hogs if h[3]), default=1),
            "chunks": view_pages.write_chunks(data_dir, f"HOGCounts_{node}", [h[3] for h in hogs], len)})

    page = os.path.join(stats_dir, "Statistics.html")
    with open(page, "w", encoding="utf-8") as out:
        out.write(HTML.replace("__NAV_CSS__", view_pages.NAV_CSS).replace("__NAV_HTML__", view_pages.NAV_HTML)
                  .replace("__NAV_JS__", view_pages.NAV_JS).replace("__THEME_JS__", view_pages.THEME_JS)
                  .replace("__DATA_JS__", view_pages.DATA_JS)
                  .replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
    return page


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>OrthoFinder statistics</title>
__THEME_JS__
<style>
:root { color-scheme:light; --page:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9;
  --axis:#c3c2b7; --ring:rgba(11,11,11,0.10); --s1:#2a78d6; --s2:#eb6834; --heat-lo:#cde2fb; --heat-hi:#104281; --zero:#f0efec;
  --accent:#2a78d6; --head:#f1f0ec; --pagenav-line:var(--grid); }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme="light"])) { color-scheme:dark; --page:#0d0d0d; --surface:#1a1a19; --ink:#ffffff;
  --ink2:#c3c2b7; --muted:#898781; --grid:#2c2c2a; --axis:#383835; --ring:rgba(255,255,255,0.10); --s1:#3987e5; --s2:#d95926;
  --heat-lo:#184f95; --heat-hi:#b7d3f6; --zero:#2a2a28; --accent:#6da7ec; --head:#242422; } }
:root[data-theme="dark"] { color-scheme:dark; --page:#0d0d0d; --surface:#1a1a19; --ink:#ffffff; --ink2:#c3c2b7; --muted:#898781;
  --grid:#2c2c2a; --axis:#383835; --ring:rgba(255,255,255,0.10); --s1:#3987e5; --s2:#d95926; --heat-lo:#184f95; --heat-hi:#b7d3f6;
  --zero:#2a2a28; --accent:#6da7ec; --head:#242422; }
* { box-sizing:border-box; }
body { margin:0; background:var(--page); color:var(--ink); font:14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width:1400px; margin:0 auto; padding:24px clamp(16px, 2vw, 32px) 48px; }
h1 { font-size:32px; margin:0 0 6px; } h2 { font-size:18px; margin:0 0 2px; } h3 { font-size:15px; margin:0; font-weight:600; }
.sub, .muted { color:var(--ink2); font-size:12px; }
nav#tabs { display:flex; gap:4px; flex-wrap:wrap; margin:16px 0; border-bottom:1px solid var(--grid); }
nav#tabs button { border:0; background:none; color:var(--ink2); font:inherit; padding:8px 12px; cursor:pointer; border-bottom:2px solid transparent; margin-bottom:-1px; }
nav#tabs button[aria-selected="true"] { color:var(--ink); border-bottom-color:var(--accent); font-weight:600; }
section[hidden] { display:none; }
.tiles { display:grid; grid-template-columns:repeat(auto-fill, minmax(170px, 1fr)); gap:10px; margin-bottom:14px; }
.tile { background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:12px 14px; }
.tile .label { color:var(--ink2); font-size:12px; } .tile .value { font-size:24px; font-weight:600; margin-top:2px; }
.tile .note { color:var(--muted); font-size:12px; }
.grid2 { display:grid; grid-template-columns:repeat(auto-fit, minmax(min(100%, 460px), 1fr)); gap:14px; margin-bottom:14px; }
.card { background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:14px; margin-bottom:14px; min-width:0; }
.card header { display:flex; justify-content:space-between; align-items:flex-start; gap:10px; margin-bottom:8px; }
.card .chart { width:100%; overflow-x:auto; }
.chart svg { display:block; }
.chart text { font-family:system-ui, -apple-system, "Segoe UI", sans-serif; }
.linkbtn { border:0; background:none; color:var(--accent); font:inherit; font-size:12px; cursor:pointer; padding:0; white-space:nowrap; }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:10px; }
input, select { font:inherit; padding:6px 10px; border:1px solid var(--axis); border-radius:6px; background:var(--surface); color:var(--ink); }
button.btn { font:inherit; padding:5px 12px; border:1px solid var(--axis); border-radius:6px; background:var(--head); color:var(--ink); cursor:pointer; }
button.btn:disabled { opacity:.4; cursor:default; }
table { border-collapse:collapse; width:100%; font-size:13px; }
th, td { text-align:left; padding:5px 8px; border-bottom:1px solid var(--grid); vertical-align:middle; }
th { color:var(--ink2); font-weight:600; white-space:nowrap; position:sticky; top:0; background:var(--surface); }
th[data-sort] { cursor:pointer; } th[data-sort]:hover { color:var(--ink); }
td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; }
.tablewrap { overflow:auto; max-height:70vh; }
.sw { display:inline-block; width:9px; height:9px; border-radius:2px; margin-right:6px; vertical-align:0; }
.tag { display:inline-block; padding:0 6px; border-radius:8px; border:1px solid var(--axis); color:var(--ink2); font-size:11px; margin-right:4px; white-space:nowrap; }
a { color:var(--accent); }
.legend { display:flex; gap:14px; flex-wrap:wrap; font-size:12px; color:var(--ink2); margin-bottom:6px; align-items:center; }
.genes-link { font-size:11px; margin-left:4px; }
.species-key { font-size:13px; color:var(--ink); margin:12px 0 0; }
.species-key details { width:100%; } .species-key summary { color:var(--ink2); }
.species-key details .legend { max-height:30vh; overflow:auto; margin-top:6px; }
.scroll { max-height:70vh; overflow:auto; }
canvas.strip { display:block; }
.ramp { display:inline-block; width:90px; height:9px; border-radius:2px; vertical-align:-1px; margin:0 6px; }
#tip { position:fixed; z-index:10; pointer-events:none; background:var(--surface); color:var(--ink); border:1px solid var(--ring);
  border-radius:8px; padding:7px 10px; font-size:12px; box-shadow:0 4px 14px rgba(0,0,0,.15); max-width:320px; }
#tip .v { font-weight:600; font-size:13px; } #tip .l { color:var(--ink2); white-space:pre-line; }
.hov { opacity:.75; }
.heat-cv.hov { opacity:1; }
details summary { cursor:pointer; color:var(--ink2); font-size:13px; margin:4px 0; }
.empty { color:var(--muted); font-size:13px; padding:12px 0; }
__NAV_CSS__
</style></head><body><main>
__NAV_HTML__
<h1>OrthoFinder statistics</h1>
<div class="sub" id="sub"></div>
<div class="legend species-key" id="species-key"></div>
<nav id="tabs" role="tablist"></nav>
<section id="tab-overview"></section>
<section id="tab-species" hidden></section>
<section id="tab-orthogroups" hidden></section>
<section id="tab-duplications" hidden></section>
<section id="tab-hogs" hidden></section>
</main>
<div id="tip" hidden><div class="v"></div><div class="l"></div></div>
<script>
const D = __DATA__;
__DATA_JS__
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const fmt = (x, d = 1) => typeof x === "number" ? x.toLocaleString(undefined, {maximumFractionDigits: d}) : (x ?? "");
const compact = x => x >= 1e6 ? fmt(x / 1e6) + "M" : x >= 1e5 ? fmt(x / 1e3, 0) + "K" : fmt(x, 1);
const tok = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const S = D.species.map(s => s[0]), NS = S.length;
const swatch = k => k >= 0 ? `<span class="sw" style="background:${D.species[k][1]}"></span>` : "";
const sum = a => a.reduce((x, y) => x + y, 0);
// Math.max(...a) throws for more than about 125k values (e.g. the orthogroups): a loop
const maxOf = (a, start = -Infinity) => a.reduce((x, y) => y > x ? y : x, start);
const minOf = (a, start = Infinity) => a.reduce((x, y) => y < x ? y : x, start);
$("sub").textContent = `OrthoFinder ${D.version}. Made from the results files (Comparative_Genomics_Statistics, Orthogroups, `
  + `Phylogenetic_Hierarchical_Orthogroups, Species_Tree). Point at a plot for its numbers; "table" shows them as a table.`;

// the species' colours, as on the tree page (and in the genes per species cells); for many, a list to open
{
  const key = D.species.map(([sp, c]) => `<span><span class="sw" style="background:${c}"></span>${esc(sp)}</span>`).join("");
  $("species-key").innerHTML = NS > 40 ? `<details><summary>The species and their colours (${NS.toLocaleString()})</summary><div class="legend">${key}</div></details>` : key;
}
// the species in the order of the species tree's leaves (any not in it after them): for plots of many species
const TREE_ORDER = (() => {
  const T = D.species_tree, seen = new Set(), order = [];
  if (T) for (const nd of T.nodes) if (typeof nd[2] === "number") { const k = T.genes[nd[2]][1]; if (k >= 0 && !seen.has(k)) { seen.add(k); order.push(k); } }
  for (let k = 0; k < NS; k++) if (!seen.has(k)) order.push(k);
  return order;
})();
__NAV_JS__

// ---- tooltip: any element with data-tv (the value) and data-tl (what it is) ----
const tip = $("tip");
let hovered = null;
function showTip(t, x, y) {
  if (hovered !== t) { hovered?.classList.remove("hov"); hovered = t; t.classList.add("hov"); }
  tip.querySelector(".v").textContent = t.dataset.tv;
  tip.querySelector(".l").textContent = t.dataset.tl || "";
  tip.hidden = false;
  const r = tip.getBoundingClientRect();
  tip.style.left = Math.min(window.innerWidth - r.width - 8, x + 14) + "px";
  tip.style.top = (y + r.height + 20 > window.innerHeight ? y - r.height - 10 : y + 14) + "px";
}
function hideTip() { tip.hidden = true; hovered?.classList.remove("hov"); hovered = null; }
document.addEventListener("pointermove", e => { const t = e.target.closest && e.target.closest("[data-tv]"); t ? showTip(t, e.clientX, e.clientY) : hideTip(); });
document.addEventListener("focusin", e => { const t = e.target.closest && e.target.closest("[data-tv]");
  if (t) { const r = t.getBoundingClientRect(); showTip(t, r.right, r.top); } else hideTip(); });
const tipAttrs = (value, label) => `data-tv="${esc(value)}" data-tl="${esc(label)}" tabindex="0"`;

// ---- shared pieces ----
function niceStep(x) { const p = Math.pow(10, Math.floor(Math.log10(x || 1))), m = x / p; return (m > 5 ? 10 : m > 2 ? 5 : m > 1 ? 2 : 1) * p; }
function ticks(max, n = 4) { const step = niceStep(max / n), top = Math.max(step, Math.ceil(max / step) * step), t = []; for (let v = 0; v <= top + step / 2; v += step) t.push(v); return t; }
const measure = document.createElement("canvas").getContext("2d");
function textW(s, size = 12) { measure.font = `${size}px system-ui, -apple-system, "Segoe UI", sans-serif`; return measure.measureText(s).width; }
function hex(c) { c = c.replace("#", ""); return [0, 2, 4].map(i => parseInt(c.slice(i, i + 2), 16)); }
function mix(a, b, t) { const A = hex(a), B = hex(b); return "#" + A.map((v, i) => Math.round(v + (B[i] - v) * t).toString(16).padStart(2, "0")).join(""); }
function heat(v, min, max) { return v ? mix(tok("--heat-lo"), tok("--heat-hi"), max > min ? (v - min) / (max - min) : 1) : tok("--zero"); }
function inkOn(bg) { const [r, g, b] = hex(bg).map(v => { v /= 255; return v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.35 ? "#0b0b0b" : "#ffffff"; }
function roundedBar(x, y, w, h, r, dir) {  // 4px rounded data end, square at the baseline
  r = Math.min(r, dir === "up" ? h : w, w / 2, h / 2);
  if (dir === "up") return `M${x},${y + h}V${y + r}Q${x},${y} ${x + r},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h}Z`;
  return `M${x},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h - r}Q${x + w},${y + h} ${x + w - r},${y + h}H${x}Z`;
}
function card(parent, title, note, id) {
  const c = document.createElement("div"); c.className = "card";
  c.innerHTML = `<header><div><h2>${esc(title)}</h2>${note ? `<div class="muted">${esc(note)}</div>` : ""}</div>`
    + `<button class="linkbtn" data-table>table</button></header><div class="chart"${id ? ` id="${id}"` : ""}></div><div class="tableview" hidden></div>`;
  parent.appendChild(c);
  return c;
}
function setTable(c, head, rows) {  // the table view of a chart: head [[label, numeric]], rows of cells
  // built when first shown (that of a species x species matrix can be large)
  const tv = c.querySelector(".tableview");
  tv.build = () => { tv.build = null; tv.innerHTML = `<div class="tablewrap"><table><thead><tr>${head.map(([h, n]) => `<th${n ? ' class="num"' : ""}>${esc(h)}</th>`).join("")}</tr></thead><tbody>`
    + rows.map(r => `<tr>${r.map((v, k) => `<td${head[k][1] ? ' class="num"' : ""}>${typeof v === "number" ? fmt(v) : v}</td>`).join("")}</tr>`).join("") + `</tbody></table></div>`; };
  if (!tv.hidden) tv.build();
}
document.addEventListener("click", e => {
  const b = e.target.closest("[data-table]");
  if (!b) return;
  const tv = b.closest(".card").querySelector(".tableview"), ch = b.closest(".card").querySelector(".chart");
  tv.hidden = !tv.hidden; if (!tv.hidden && tv.build) tv.build(); ch.hidden = !tv.hidden; b.textContent = tv.hidden ? "table" : "chart";
});

// a column chart of counts in ordered bins: bins [[label, value]]
function columns(el, bins, {unit, xTitle}) {
  const W = Math.max(300, el.clientWidth), H = 210, left = 46, right = 8, top = 18, bottom = 40;
  const pw = W - left - right, ph = H - top - bottom, max = Math.max(1, ...bins.map(b => b[1]));
  const t = ticks(max), ymax = t[t.length - 1], band = pw / bins.length, bw = Math.min(24, band * 0.7);
  const Y = v => top + ph - v / ymax * ph, out = [];
  for (const v of t) out.push(`<line x1="${left}" x2="${W - right}" y1="${Y(v)}" y2="${Y(v)}" stroke="${v ? "var(--grid)" : "var(--axis)"}"/>`
    + `<text x="${left - 6}" y="${Y(v)}" dy="0.32em" text-anchor="end" font-size="11" fill="var(--muted)">${compact(v)}</text>`);
  const showValues = bins.length <= 14, labelEvery = Math.ceil(bins.length / Math.max(1, Math.floor(pw / 34)));
  bins.forEach(([label, v], k) => {
    const cx = left + band * (k + 0.5), h = ph * v / ymax;
    out.push(`<g ${tipAttrs(`${fmt(v)} ${unit}`, `${xTitle}: ${label}`)}><rect x="${cx - band / 2}" y="${top}" width="${band}" height="${ph}" fill="transparent"/>`
      + (v ? `<path d="${roundedBar(cx - bw / 2, Y(v), bw, h, 4, "up")}" fill="var(--s1)"/>` : "") + `</g>`);
    if (showValues && v) out.push(`<text x="${cx}" y="${Y(v) - 5}" text-anchor="middle" font-size="11" fill="var(--ink2)">${compact(v)}</text>`);
    if (k % labelEvery === 0) out.push(`<text x="${cx}" y="${H - bottom + 16}" text-anchor="middle" font-size="11" fill="var(--muted)">${esc(label)}</text>`);
  });
  out.push(`<text x="${left + pw / 2}" y="${H - 6}" text-anchor="middle" font-size="11" fill="var(--ink2)">${esc(xTitle)}</text>`);
  el.innerHTML = `<svg width="${W}" height="${H}" role="img">${out.join("")}</svg>`;
}

// horizontal bars: rows [{label, k (species number or -1), value, tl}]
function hbars(el, rows, {unit}) {
  const W = Math.max(300, el.clientWidth), rowH = 26, top = 4;
  const labelW = Math.min(240, Math.max(60, ...rows.map(r => textW(r.label) + (r.k >= 0 ? 16 : 0)))) + 10;
  const max = Math.max(1, ...rows.map(r => r.value)), valW = textW(compact(max)) + 12, pw = W - labelW - valW - 8;
  const out = [], H = top + rows.length * rowH + 4;
  rows.forEach((r, i) => {
    const y = top + i * rowH, w = pw * r.value / max;
    let x = 0;
    if (r.k >= 0) { out.push(`<rect x="0" y="${y + 9}" width="9" height="9" rx="2" fill="${D.species[r.k][1]}"/>`); x = 15; }
    out.push(`<text x="${x}" y="${y + rowH / 2}" dy="0.32em" font-size="12" fill="var(--ink)">${esc(fit(r.label, labelW - x - 8))}</text>`
      + `<g ${tipAttrs(`${fmt(r.value)} ${unit}`, r.tl || r.label)}><rect x="${labelW}" y="${y}" width="${pw + valW}" height="${rowH}" fill="transparent"/>`
      + (r.value ? `<path d="${roundedBar(labelW, y + 6, Math.max(1, w), rowH - 12, 4, "right")}" fill="var(--s1)"/>` : "") + `</g>`
      + `<text x="${labelW + w + 6}" y="${y + rowH / 2}" dy="0.32em" font-size="11" fill="var(--ink2)">${compact(r.value)}</text>`);
  });
  out.push(`<line x1="${labelW}" x2="${labelW}" y1="${top}" y2="${H - 4}" stroke="var(--axis)"/>`);
  el.innerHTML = `<svg width="${W}" height="${H}" role="img">${out.join("")}</svg>`;
}
function fit(s, w, size = 12) {
  if (textW(s, size) <= w) return s;
  let lo = 0, hi = s.length;
  while (lo < hi) { const m = (lo + hi + 1) >> 1; if (textW(s.slice(0, m) + "…", size) <= w) lo = m; else hi = m - 1; }
  return s.slice(0, lo) + "…";
}

// stacked horizontal bars of two parts, as percentages: rows [{label, k, parts: [a, b]}], names of the parts
function stacked(el, rows, names) {
  const W = Math.max(300, el.clientWidth), rowH = 28, top = 2, cols = ["var(--s1)", "var(--s2)"];
  const labelW = Math.min(240, Math.max(60, ...rows.map(r => textW(r.label) + 16))) + 10, pw = W - labelW - 20;
  const out = [], H = top + rows.length * rowH + 22;
  rows.forEach((r, i) => {
    const y = top + i * rowH, total = sum(r.parts) || 1;
    out.push(`<rect x="0" y="${y + 9}" width="9" height="9" rx="2" fill="${D.species[r.k][1]}"/>`
      + `<text x="15" y="${y + rowH / 2}" dy="0.32em" font-size="12" fill="var(--ink)">${esc(fit(r.label, labelW - 23))}</text>`);
    let x = labelW;
    r.parts.forEach((v, j) => {
      const w = pw * v / total, pct = 100 * v / total, last = j === r.parts.length - 1;
      if (w > 0) {
        const ww = Math.max(0.5, w - (last ? 0 : 2));  // a 2px surface gap between the parts
        out.push(`<g ${tipAttrs(`${fmt(v)} genes (${fmt(pct)}%)`, `${r.label}: ${names[j]}`)}>`
          + `<path d="${last ? roundedBar(x, y + 5, ww, rowH - 10, 4, "right") : `M${x},${y + 5}h${ww}v${rowH - 10}h${-ww}Z`}" fill="${cols[j]}"/></g>`);
        const label = `${fmt(pct, 0)}%`;
        if (textW(label, 11) + 10 < ww) out.push(`<text x="${x + 6}" y="${y + rowH / 2}" dy="0.32em" font-size="11" fill="${inkOn(tok(j ? "--s2" : "--s1"))}" pointer-events="none">${label}</text>`);
      }
      x += w;
    });
  });
  for (const p of [0, 25, 50, 75, 100]) {
    const x = labelW + pw * p / 100;
    out.push(`<text x="${x}" y="${H - 6}" text-anchor="${p === 0 ? "start" : p === 100 ? "end" : "middle"}" font-size="11" fill="var(--muted)">${p}%</text>`);
  }
  el.innerHTML = `<div class="legend">${names.map((n, j) => `<span><span class="sw" style="background:${cols[j]}"></span>${esc(n)}</span>`).join("")}</div>`
    + `<svg width="${W}" height="${H}" role="img">${out.join("")}</svg>`;
}

// a species x species heatmap: in SVG, or above HEAT_SVG_MAX species (an SVG of NS x NS cells
// is too slow) on a canvas, a cell for each pair of species, named by pointing at it
const HEAT_SVG_MAX = 80;
function heatmap(el, M, {unit, diagonal}) {
  const vals = M.flatMap((r, a) => r.filter((v, b) => (a !== b || diagonal) && typeof v === "number" && v > 0));
  const min = vals.length ? minOf(vals) : 0, max = vals.length ? maxOf(vals) : 1;
  const legend = `<div class="legend">fewer<span class="ramp" style="background:linear-gradient(90deg, ${tok("--heat-lo")}, ${tok("--heat-hi")})"></span>more `
    + `(${fmt(min)} to ${fmt(max)} ${esc(unit)}; 0 in grey; rows → columns${NS > HEAT_SVG_MAX ? `; ${D.species_tree ? "in the order of the species tree; " : ""}point at a cell for its species` : ""})</div>`;
  if (NS > HEAT_SVG_MAX) return heatCanvas(el, M, {unit, diagonal, min, max, legend});
  const W = Math.max(300, el.clientWidth), labelW = Math.min(220, Math.max(...S.map(s => textW(s) + 16))) + 8;
  const cell = Math.max(12, Math.min(52, Math.floor((W - labelW - 20) / NS))), colH = Math.min(150, Math.max(...S.map(s => textW(s, 11))) * 0.72 + 14);
  const showValues = cell >= 34, out = [], top = colH;
  S.forEach((s, b) => {
    const x = labelW + b * cell + cell / 2;
    out.push(`<text transform="translate(${x + 3},${top - 6}) rotate(-45)" font-size="11" fill="var(--ink2)">${esc(fit(s, 190, 11))}</text>`);
  });
  M.forEach((row, a) => {
    const y = top + a * cell;
    out.push(`<rect x="0" y="${y + cell / 2 - 4.5}" width="9" height="9" rx="2" fill="${D.species[a][1]}"/>`
      + `<text x="15" y="${y + cell / 2}" dy="0.32em" font-size="12" fill="var(--ink)">${esc(fit(S[a], labelW - 23))}</text>`);
    row.forEach((v, b) => {
      const x = labelW + b * cell, onDiag = a === b && !diagonal;
      const bg = onDiag || typeof v !== "number" ? tok("--zero") : heat(v, min, max);
      out.push(`<g ${tipAttrs(onDiag ? "—" : `${fmt(v)} ${unit}`, `${S[a]} → ${S[b]}`)}>`
        + `<rect x="${x + 1}" y="${y + 1}" width="${cell - 2}" height="${cell - 2}" rx="2" fill="${bg}"/></g>`);
      if (showValues && !onDiag && typeof v === "number")
        out.push(`<text x="${x + cell / 2}" y="${y + cell / 2}" dy="0.32em" text-anchor="middle" font-size="11" fill="${inkOn(bg)}" pointer-events="none">${compact(v)}</text>`);
    });
  });
  const H = top + NS * cell + 6;
  el.innerHTML = legend + `<svg width="${Math.max(W, labelW + NS * cell + 10)}" height="${H}" role="img">${out.join("")}</svg>`;
}
function heatCanvas(el, M, {unit, diagonal, min, max, legend}) {
  // a pixel for each cell, scaled up unsmoothed; the species' colours along the top and the left;
  // the species in the order of the species tree (P: position -> species), so that clades are blocks
  const P = TREE_ORDER, W = Math.max(300, el.clientWidth), bar = 6, at = bar + 2;
  const cell = Math.max(1, Math.min(12, Math.floor((W - at - 4) / NS))), side = at + NS * cell;
  const lo = hex(tok("--heat-lo")), hi = hex(tok("--heat-hi")), zero = hex(tok("--zero"));
  const img = new ImageData(NS, NS), px = img.data;
  for (let a = 0; a < NS; a++) for (let b = 0; b < NS; b++) {
    const v = M[P[a]][P[b]], i = 4 * (a * NS + b), t = max > min ? (v - min) / (max - min) : 1;
    const c = (a === b && !diagonal) || typeof v !== "number" || v <= 0 ? null : t;
    for (let j = 0; j < 3; j++) px[i + j] = c === null ? zero[j] : Math.round(lo[j] + (hi[j] - lo[j]) * c);
    px[i + 3] = 255;
  }
  const small = document.createElement("canvas");
  small.width = small.height = NS; small.getContext("2d").putImageData(img, 0, 0);
  const dpr = window.devicePixelRatio || 1, cv = document.createElement("canvas");
  cv.width = cv.height = Math.round(side * dpr); cv.style.width = cv.style.height = side + "px";
  cv.className = "heat-cv"; cv.tabIndex = 0; cv.setAttribute("role", "img");
  const ctx = cv.getContext("2d");
  ctx.scale(dpr, dpr); ctx.imageSmoothingEnabled = false;
  P.forEach((sp, k) => { ctx.fillStyle = D.species[sp][1]; ctx.fillRect(at + k * cell, 0, cell, bar); ctx.fillRect(0, at + k * cell, bar, cell); });
  ctx.drawImage(small, at, at, NS * cell, NS * cell);
  // the tooltip (pointermove on the document, after this): the cell's (or the colour's) species
  cv.addEventListener("pointermove", e => {
    const r = cv.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const b = Math.floor((x - at) / cell), a = Math.floor((y - at) / cell), ok = k => k >= 0 && k < NS;
    let tv = null, tl = "";
    if (ok(a) && ok(b)) { const v = M[P[a]][P[b]]; tv = (a === b && !diagonal) || typeof v !== "number" ? "—" : `${fmt(v)} ${unit}`; tl = `${S[P[a]]} → ${S[P[b]]}`; }
    else if (y < bar && ok(b)) [tv, tl] = [S[P[b]], "column"];
    else if (x < bar && ok(a)) [tv, tl] = [S[P[a]], "row"];
    if (tv === null) { delete cv.dataset.tv; delete cv.dataset.tl; } else Object.assign(cv.dataset, {tv, tl});
  });
  el.innerHTML = legend + `<div class="tablewrap"></div>`;
  el.lastChild.appendChild(cv);
}

// the genes of each species of an orthogroup or HOG (pairs: [species, genes, species, genes, ...], those with
// genes) as a strip, on a canvas: a cell for each species in its colour, stronger for more genes (grey: none),
// with the number when it fits; for many species, each pixel some species (in the order of the species tree),
// in the colour of the one with the most genes. Pointing at it names the species there.
const STRIP_NUMBERS = NS <= 20, STRIP_CW = STRIP_NUMBERS ? 22 : NS <= 60 ? 6 : NS <= 160 ? 3 : 0;
const STRIP_W = STRIP_CW ? NS * STRIP_CW : 480, STRIP_H = STRIP_NUMBERS ? 16 : 12;
const STRIP_ORDER = STRIP_CW ? S.map((_, k) => k) : TREE_ORDER, STRIP_POS = new Int32Array(NS);
STRIP_ORDER.forEach((k, p) => { STRIP_POS[k] = p; });
const stripX = k => STRIP_CW ? k * STRIP_CW : Math.floor(STRIP_POS[k] * STRIP_W / NS);
const stripCanvas = () => `<canvas class="strip" width="${STRIP_W}" height="${STRIP_H}" style="width:${STRIP_W}px;height:${STRIP_H}px" tabindex="0"></canvas>`;
function cell(ctx, x, w) { ctx.beginPath(); if (ctx.roundRect && w > 3) ctx.roundRect(x, 0, w, STRIP_H, 2); else ctx.rect(x, 0, w, STRIP_H); ctx.fill(); }
function drawStrip(cv, pairs, max) {  // pairs null: the heading, the species' colours
  const dpr = window.devicePixelRatio || 1, ctx = cv.getContext("2d"), surface = tok("--surface"), zero = tok("--zero");
  cv.width = Math.round(STRIP_W * dpr); cv.height = Math.round(STRIP_H * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  Object.assign(cv, {pairs, max});
  const colour = (k, c) => pairs ? mix(surface, D.species[k][1], 0.3 + 0.7 * Math.log1p(c) / Math.log1p(max)) : D.species[k][1];
  const w = STRIP_CW ? STRIP_CW - 2 : 1;
  if (STRIP_CW) { ctx.fillStyle = zero; for (let k = 0; k < NS; k++) cell(ctx, stripX(k), w); }
  else { ctx.fillStyle = pairs ? zero : surface; ctx.fillRect(0, 0, STRIP_W, STRIP_H); }
  const at = new Map();  // pixel -> [genes, species] of the species there with the most genes
  if (pairs) for (let j = 0; j < pairs.length; j += 2) { const x = stripX(pairs[j]), b = at.get(x); if (!b || pairs[j + 1] > b[0]) at.set(x, [pairs[j + 1], pairs[j]]); }
  else for (let k = NS - 1; k >= 0; k--) at.set(stripX(k), [1, k]);
  ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.font = "10px system-ui, sans-serif";
  for (const [x, [c, k]] of at) {
    const bg = colour(k, c);
    ctx.fillStyle = bg; cell(ctx, x, w);
    if (STRIP_NUMBERS && pairs) { ctx.fillStyle = inkOn(bg); ctx.fillText(c > 999 ? compact(c) : String(c), x + w / 2, STRIP_H / 2 + 0.5); }
  }
}
// pointing at a strip: the species there (before the tooltip reads it: the capture phase)
document.addEventListener("pointermove", e => {
  const cv = e.target;
  if (!(cv instanceof HTMLCanvasElement) || !cv.classList.contains("strip")) return;
  const x = e.clientX - cv.getBoundingClientRect().left;
  let ks = [];
  if (STRIP_CW) { const k = Math.floor(x / STRIP_CW); if (k >= 0 && k < NS) ks = [k]; }
  else { const p0 = Math.max(0, Math.floor(x * NS / STRIP_W)), p1 = Math.min(NS, Math.max(p0 + 1, Math.floor((x + 1) * NS / STRIP_W))); ks = STRIP_ORDER.slice(p0, p1); }
  if (!cv.pairs) { cv.dataset.tv = ks.map(k => S[k]).join(", "); cv.dataset.tl = "the species' colour"; return; }
  const genes = new Map();
  for (let j = 0; j < cv.pairs.length; j += 2) genes.set(cv.pairs[j], cv.pairs[j + 1]);
  cv.dataset.tv = `${fmt(sum(cv.pairs.filter((_, j) => j % 2)))} genes, ${fmt(cv.pairs.length / 2)} species`;
  cv.dataset.tl = ks.map(k => `${S[k]}: ${genes.get(k) || 0}`).join("\n");
}, true);
// the column's heading: the species of the cells, in their colours (point at one for its name)
const stripHead = () => `Genes per species${STRIP_CW ? "" : " (species tree order)"}<br>` + stripCanvas().replace('class="strip"', 'class="strip" data-head');
function drawStrips(el, pairsOf, max) {
  for (const cv of el.querySelectorAll("canvas.strip")) drawStrip(cv, cv.hasAttribute("data-head") ? null : pairsOf(+cv.dataset.r), max);
}

// a table of many rows, sortable and in pages; prepare(rows of the page) loads what they need, after(el) draws
function dataTable(el, {cols, rows, pageSize = 50, sortCol = null, sortDesc = true, prepare, after}) {
  let page = 0, order = rows, seq = 0;
  const sortBy = (c, desc) => { sortCol = c; sortDesc = desc; const key = cols[c].key;
    order = [...rows].sort((a, b) => { const x = key(a), y = key(b); return (x < y ? -1 : x > y ? 1 : 0) * (desc ? -1 : 1); }); page = 0; };
  if (sortCol !== null) sortBy(sortCol, sortDesc);
  async function render() {
    const my = ++seq, n = order.length, pages = Math.max(1, Math.ceil(n / pageSize)), p0 = page * pageSize, view = order.slice(p0, p0 + pageSize);
    if (prepare) {
      try { await prepare(view); } catch (e) { if (my === seq) el.innerHTML = `<div class="empty">${esc(e.message)}</div>`; return; }
      if (my !== seq) return;  // another page asked for meanwhile
    }
    el.innerHTML = `<div class="tablewrap"><table><thead><tr>` + cols.map((c, k) =>
      `<th${c.num ? ' class="num"' : ""}${c.key ? ` data-sort="${k}"` : ""}>${c.head || esc(c.label)}${sortCol === k ? (sortDesc ? " ▾" : " ▴") : ""}</th>`).join("")
      + `</tr></thead><tbody>` + (view.length ? view.map(r => `<tr>${cols.map(c => `<td${c.num ? ' class="num"' : ""}>${c.html(r)}</td>`).join("")}</tr>`).join("")
        : `<tr><td colspan="${cols.length}" class="empty">Nothing matches.</td></tr>`)
      + `</tbody></table></div><div class="tools" style="margin-top:8px"><button class="btn" data-p="-1"${page ? "" : " disabled"}>Previous</button>`
      + `<button class="btn" data-p="1"${page < pages - 1 ? "" : " disabled"}>Next</button>`
      + `<span class="muted">${n ? `${fmt(p0 + 1)}–${fmt(Math.min(n, p0 + pageSize))} of ${fmt(n)}` : "0"}</span></div>`;
    if (after) after(el);
  }
  el.onclick = e => {
    const th = e.target.closest("th[data-sort]"), b = e.target.closest("[data-p]");
    if (th) { const c = +th.dataset.sort; sortBy(c, sortCol === c ? !sortDesc : true); render(); }
    if (b) { page += +b.dataset.p; render(); }
  };
  render();
  return {setRows(r) { rows = r; if (sortCol !== null) sortBy(sortCol, sortDesc); else { order = r; page = 0; } render(); }, redraw: render};
}
const treeLink = og => D.tree_page ? `<a href="${esc(D.tree_page)}#${encodeURIComponent(og)}" title="Its gene tree and alignment">${esc(og)}</a>` : esc(og);
// an orthogroup's ID: a link to its tree (if it has one), and to its genes (--groups-view)
const ogLinks = (og, hasTree) => (hasTree ? treeLink(og) : esc(og))
  + (D.groups_page ? ` <a class="genes-link" href="${esc(D.groups_page)}#orthogroups/${encodeURIComponent(og)}" title="Its genes">genes</a>` : "");
function saveText(text, filename) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], {type: "text/tab-separated-values"}));
  a.download = filename; document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// ---- the tabs ----
const TABS = [["overview", "Overview"], ["species", "Species"], ["orthogroups", "Orthogroups"], ["duplications", "Duplications"], ["hogs", "HOGs"]];
const drawn = {}, renderers = {}, onResize = {};  // onResize: {plot: draw it again}
$("tabs").innerHTML = TABS.map(([id, name]) => `<button role="tab" data-tab="${id}">${name}</button>`).join("");
function openTab(id) {
  if (!renderers[id]) id = "overview";
  hideTip();
  for (const [t] of TABS) { $("tab-" + t).hidden = t !== id; document.querySelector(`[data-tab="${t}"]`).setAttribute("aria-selected", t === id); }
  if (!drawn[id]) { drawn[id] = true; renderers[id](); }
  history.replaceState(null, "", "#" + id);
}
$("tabs").onclick = e => { const b = e.target.closest("[data-tab]"); if (b) openTab(b.dataset.tab); };

// the orthogroups, in columns (D.ogs): a row number for each; the genes of each species in it are in chunks (Counts_<c>)
const O = D.ogs, OGIDS = expandIds(O.ids), NOG = OGIDS.length, OG_ROWS = Array.from({length: NOG}, (_, r) => r);
const isAll = r => O.species[r] === NS, isSingle = r => isAll(r) && O.genes[r] === NS, isOne = r => O.species[r] === 1;

// Overview
renderers.overview = () => {
  const sec = $("tab-overview"), OV = Object.fromEntries(D.overall || []);
  const nAll = OG_ROWS.filter(isAll).length, nSingle = OG_ROWS.filter(isSingle).length;
  const tiles = [
    ["Species", NS], ["Genes", OV["Number of genes"]],
    ["Genes in orthogroups", OV["Percentage of genes in orthogroups"], "%", OV["Number of genes in orthogroups"] !== undefined ? `${fmt(OV["Number of genes in orthogroups"])} genes` : ""],
    ["Orthogroups", NOG], ["With all species", nAll, "", `${fmt(100 * nAll / (NOG || 1))}% of orthogroups`],
    ["Single-copy orthogroups", nSingle, "", "one gene in each species"],
    ["Species-specific orthogroups", OV["Number of species-specific orthogroups"]],
    ["Unassigned genes", OV["Number of unassigned genes"]],
  ].filter(t => t[1] !== undefined);
  sec.innerHTML = `<div class="tiles">${tiles.map(([l, v, u, n]) => `<div class="tile"><div class="label">${esc(l)}</div>`
    + `<div class="value">${typeof v === "number" ? (u === "%" ? fmt(v) + "%" : compact(v)) : esc(v)}</div>${n ? `<div class="note">${esc(n)}</div>` : ""}</div>`).join("")}</div><div class="grid2" id="ov-grid"></div>`;
  const g = $("ov-grid");
  const c1 = card(g, "Orthogroup sizes", "Orthogroups by their number of genes");
  const c2 = card(g, "Species per orthogroup", "Orthogroups by the number of species with genes in them");
  const bySp = Array.from({length: NS}, (_, k) => [String(k + 1), 0]);
  for (const n of O.species) if (n) bySp[n - 1][1]++;
  const draw = () => { columns(c1.querySelector(".chart"), D.og_sizes, {unit: "orthogroups", xTitle: "Genes"});
                       columns(c2.querySelector(".chart"), bySp, {unit: "orthogroups", xTitle: "Species"}); };
  draw(); onResize["overview"] = draw;
  setTable(c1, [["Genes", 0], ["Orthogroups", 1]], D.og_sizes.map(([l, v]) => [esc(l), v]));
  setTable(c2, [["Species", 0], ["Orthogroups", 1]], bySp.map(([l, v]) => [l, v]));
  if (D.overall) {
    const c = card(sec, "All overall statistics", "Statistics_Overall.tsv");
    c.querySelector("[data-table]").remove();
    c.querySelector(".chart").innerHTML = `<div class="tablewrap"><table><tbody>${D.overall.map(([k, v]) => `<tr><td>${esc(k)}</td><td class="num">${typeof v === "number" ? fmt(v) : esc(v)}</td></tr>`).join("")}</tbody></table></div>`;
  }
};

// Species
const BIG_TABLE = 100;  // above this many species, a species x species matrix is saved as a file, not shown as a table
renderers.species = () => {
  const sec = $("tab-species"), P = Object.fromEntries((D.per_species || []).map(([k, v]) => [k, v]));
  if (P["Number of genes in orthogroups"] && P["Number of unassigned genes"]) {
    const c = card(sec, "Genes in orthogroups", "Of each species' genes, those in orthogroups and those not assigned to any");
    if (NS > 30) c.querySelector(".chart").classList.add("scroll");
    const rows = S.map((s, k) => ({label: s, k, parts: [P["Number of genes in orthogroups"][k] || 0, P["Number of unassigned genes"][k] || 0]}));
    const draw = () => stacked(c.querySelector(".chart"), rows, ["In orthogroups", "Unassigned"]);
    draw(); onResize["species-stacked"] = draw;
    setTable(c, [["Species", 0], ["In orthogroups", 1], ["Unassigned", 1], ["% in orthogroups", 1]],
      rows.map(r => [swatch(r.k) + esc(r.label), r.parts[0], r.parts[1], 100 * r.parts[0] / (sum(r.parts) || 1)]));
  }
  const c = card(sec, "Shared between species", "");
  c.querySelector("header > div").insertAdjacentHTML("beforeend", `<div style="margin-top:6px"><select id="mat">${D.matrices.map((n, i) => `<option value="${i}">${esc(n)}</option>`).join("")}</select></div>`);
  let M = null, name = "";
  if (NS > BIG_TABLE) {
    const b = c.querySelector("[data-table]");
    b.removeAttribute("data-table"); b.textContent = "save as TSV"; b.title = "The matrix as a tab-separated file";
    b.onclick = () => { if (M) saveText(["", ...S].join("\t") + "\n" + M.map((r, a) => [S[a], ...r.map(v => v ?? "")].join("\t")).join("\n") + "\n", name.replace(/[^\w-]+/g, "_") + ".tsv"); };
  }
  const drawM = () => {
    const i = +$("mat").value, og = i === 0;
    name = D.matrices[i];
    c.querySelector(".chart").innerHTML = `<div class="empty">Loading ${esc(name)}...</div>`;
    loadData(`Matrix_${i}`).then(m => {
      if (+$("mat").value !== i) return;
      M = m;
      heatmap(c.querySelector(".chart"), M, {unit: og ? "orthogroups" : "orthologs", diagonal: og});
      if (NS <= BIG_TABLE) setTable(c, [["", 0], ...S.map(s => [s, 1])], M.map((r, a) => [swatch(a) + esc(S[a]), ...r.map((v, b) => a === b && !og ? "—" : v)]));
    }, err => { c.querySelector(".chart").innerHTML = `<div class="empty">${esc(name)} ${esc(err.message)}.</div>`; });
  };
  $("mat").onchange = drawM; drawM(); onResize["species-matrix"] = () => { if (M) heatmap(c.querySelector(".chart"), M, {unit: +$("mat").value === 0 ? "orthogroups" : "orthologs", diagonal: +$("mat").value === 0}); };
  sec.insertAdjacentHTML("beforeend", `<div class="muted" style="margin:-6px 0 14px">Orthogroups shared: orthogroups with genes of both species (on the diagonal, those with genes of the species). `
    + `Orthologs: as in Comparative_Genomics_Statistics/OrthologsStats_*.tsv, from the row's species to the column's.</div>`);
  if (D.per_species) {
    const t = card(sec, "All per-species statistics", "Statistics_PerSpecies.tsv");
    t.querySelector("[data-table]").remove();
    t.querySelector(".chart").innerHTML = `<div class="tablewrap"><table><thead><tr><th></th>${S.map((s, k) => `<th class="num">${swatch(k)}${esc(s)}</th>`).join("")}</tr></thead><tbody>`
      + D.per_species.map(([k, v]) => `<tr><td>${esc(k)}</td>${v.map(x => `<td class="num">${typeof x === "number" ? fmt(x) : esc(x ?? "")}</td>`).join("")}</tr>`).join("") + `</tbody></table></div>`;
  }
};

// Orthogroups
renderers.orthogroups = () => {
  const sec = $("tab-orthogroups");
  sec.innerHTML = `<div class="tools"><input id="og-q" placeholder="Search an orthogroup ID" autocomplete="off">`
    + `<select id="og-kind"><option value="">All orthogroups</option><option value="all">With all species</option><option value="single">Single-copy (one gene in each species)</option>`
    + `<option value="one">Species-specific (one species)</option><option value="dup">With gene duplications</option><option value="tree">With a gene tree</option></select>`
    + `<select id="og-sp"><option value="-1">Any species</option>${S.map((s, k) => `<option value="${k}">With ${esc(s)}</option>`).join("")}</select>`
    + `<button class="btn" id="og-dl">Download as TSV</button><span class="muted" id="og-msg"></span></div><div class="card"><div id="og-table"></div></div>`
    + `<div class="muted">Genes per species: ${STRIP_CW ? `a cell for each species (${S.length}), in its colour (as in the column's heading)`
      : `the ${fmt(NS)} species in the order of the species tree, a few in each pixel, in the colour of the one with the most genes`}, `
    + `stronger for more genes, grey for none${STRIP_NUMBERS ? "; the number of genes in it" : ""}. Point at it for the species there and their genes.</div>`;
  let counts = new Map();
  const cols = [
    {label: "Orthogroup", key: r => OGIDS[r], html: r => ogLinks(OGIDS[r], O.tree[r])},
    {label: "Genes", key: r => O.genes[r], html: r => fmt(O.genes[r]), num: true},
    {label: "Species", key: r => O.species[r], html: r => `${fmt(O.species[r])} of ${fmt(NS)}`, num: true},
    {label: "Genes per species", head: stripHead(), html: r => stripCanvas().replace('class="strip"', `class="strip" data-r="${r}"`)},
    {label: "Duplications (support ≥ 0.5)", key: r => O.dups05[r], html: r => fmt(O.dups05[r]), num: true},
    {label: "", html: r => (isSingle(r) ? `<span class="tag">single-copy</span>` : isAll(r) ? `<span class="tag">all species</span>` : "")
      + (isOne(r) ? `<span class="tag">species-specific</span>` : "")},
  ];
  const table = dataTable($("og-table"), {cols, rows: OG_ROWS,
    prepare: view => chunkRows("Counts", D.counts, view).then(m => { counts = m; }),
    after: el => drawStrips(el, r => counts.get(r), D.max_genes)});
  onResize["orthogroups"] = () => table.redraw();
  let filtered = OG_ROWS, applying = 0;
  const apply = async () => {
    const my = ++applying, q = $("og-q").value.trim().toLowerCase(), kind = $("og-kind").value, sp = +$("og-sp").value;
    let withSp = null;
    if (sp >= 0) {
      $("og-msg").textContent = `Loading the orthogroups of ${S[sp]}...`;
      try { withSp = new Set(await loadData(`Species_${sp}`)); } catch (e) { $("og-msg").textContent = `The orthogroups of ${S[sp]} ${e.message}.`; return; }
      if (my !== applying) return;
      $("og-msg").textContent = "";
    }
    filtered = OG_ROWS.filter(r => (!q || OGIDS[r].toLowerCase().includes(q)) && (!withSp || withSp.has(r))
      && (!kind || (kind === "all" ? isAll(r) : kind === "single" ? isSingle(r) : kind === "one" ? isOne(r) : kind === "dup" ? O.dups05[r] > 0 : O.tree[r])));
    table.setRows(filtered);
  };
  for (const id of ["og-q", "og-kind", "og-sp"]) $(id).addEventListener("input", apply);
  $("og-dl").onclick = async () => {
    // the genes of each species: a column for each, or for many species, the species with genes in one column
    $("og-msg").textContent = `Preparing ${fmt(filtered.length)} orthogroups...`;
    let m;
    try { m = await chunkRows("Counts", D.counts, filtered); } catch (e) { $("og-msg").textContent = `The genes per species ${e.message}.`; return; }
    $("og-msg").textContent = "";
    const wide = NS <= 200, lines = [["Orthogroup", ...(wide ? S : ["Genes per species"]), "Total", "Duplications", "Duplications (support >= 0.5)"].join("\t")];
    for (const r of filtered) {
      const p = m.get(r), cells = [];
      if (wide) { const g = new Array(NS).fill(0); for (let j = 0; j < p.length; j += 2) g[p[j]] = p[j + 1]; cells.push(...g); }
      else { const parts = []; for (let j = 0; j < p.length; j += 2) parts.push(`${S[p[j]]}: ${p[j + 1]}`); cells.push(parts.join(", ")); }
      lines.push([OGIDS[r], ...cells, O.genes[r], O.dups[r], O.dups05[r]].join("\t"));
    }
    saveText(lines.join("\n") + "\n", "Orthogroups_selected.tsv");
  };
};

// Duplications
function speciesTree(el, T, values) {  // values: {node name: [all, support >= 0.5]}
  const n = T.nodes.length, kids = Array.from({length: n}, () => []), depth = new Float64Array(n);
  const unit = T.nodes.every(([p, d]) => p < 0 || !d);
  for (let k = 1; k < n; k++) { const p = T.nodes[k][0]; kids[p].push(k); depth[k] = depth[p] + (unit ? 1 : T.nodes[k][1]); }
  const name = k => kids[k].length ? T.nodes[k][2] : T.genes[T.nodes[k][2]][0];
  const y = new Float64Array(n); let row = 0;
  const order = [], st = [0];
  while (st.length) { const k = st.pop(); order.push(k); if (kids[k].length) for (let j = kids[k].length - 1; j >= 0; j--) st.push(kids[k][j]); else y[k] = row++; }
  for (let j = order.length - 1; j >= 0; j--) { const k = order[j]; if (kids[k].length) y[k] = (y[kids[k][0]] + y[kids[k][kids[k].length - 1]]) / 2; }
  // for many species, closer rows and smaller circles, without the numbers (pointing at a node gives them)
  const rowH = row > 150 ? 16 : row > 60 ? 24 : 42, roomy = rowH >= 42, rMax = roomy ? 12 : rowH / 2 - 1;
  const W = Math.max(320, el.clientWidth), top = roomy ? 36 : 12, labelW = Math.min(260, maxOf(T.genes.map(g => textW(g[0]) + 16), 0)) + 70;
  const maxD = maxOf(depth, 0) || 1, sx = (W - labelW - 40) / maxD, X = k => 24 + depth[k] * sx, Y = k => top + y[k] * rowH;
  const vmax = maxOf(Object.values(values).map(v => v[1]), 1), out = [];
  for (const k of order) {
    const p = T.nodes[k][0];
    out.push(`<line x1="${p >= 0 ? X(p) : 8}" x2="${X(k)}" y1="${Y(k)}" y2="${Y(k)}" stroke="var(--axis)" stroke-width="1.5"/>`);
    if (kids[k].length) out.push(`<line x1="${X(k)}" x2="${X(k)}" y1="${Y(kids[k][0])}" y2="${Y(kids[k][kids[k].length - 1])}" stroke="var(--axis)" stroke-width="1.5"/>`);
  }
  for (const k of order) {
    const v = values[name(k)] || [0, 0], r = Math.min(rMax, 2 + (rMax - 2) * Math.sqrt(v[1] / vmax)), x = X(k), yk = Y(k), leaf = !kids[k].length;
    out.push(`<g ${tipAttrs(`${fmt(v[1])} duplications with support ≥ 0.5`, `${name(k)}${leaf ? "" : " (the species below it)"}: ${fmt(v[0])} in all`)}>`
      + `<circle cx="${x}" cy="${yk}" r="${Math.max(r, rowH / 2)}" fill="transparent"/><circle cx="${x}" cy="${yk}" r="${r}" fill="var(--s1)" fill-opacity="0.22" stroke="var(--s1)" stroke-width="1.5"/></g>`
      + (roomy ? `<text x="${x}" y="${yk - r - 4}" text-anchor="middle" font-size="11" fill="var(--ink2)" pointer-events="none">${compact(v[1])}</text>` : ""));
    if (leaf) {
      const sp = T.genes[T.nodes[k][2]][1], fs = Math.min(12, rowH - 4);
      out.push((sp >= 0 ? `<rect x="${x + rMax + 8}" y="${yk - 4.5}" width="9" height="9" rx="2" fill="${D.species[sp][1]}"/>` : "")
        + `<text x="${x + rMax + 23}" y="${yk}" dy="0.32em" font-size="${fs}" fill="var(--ink)">${esc(name(k))}</text>`);
    } else if (roomy) out.push(`<text x="${x - 6}" y="${yk + 14}" text-anchor="end" font-size="10" fill="var(--muted)">${esc(name(k))}</text>`);
  }
  el.innerHTML = `<svg width="${W}" height="${top + row * rowH + (roomy ? 0 : 8)}" role="img">${out.join("")}</svg>`;
}
renderers.duplications = () => {
  const sec = $("tab-duplications");
  if (D.species_tree && D.node_duplications) {
    const c = card(sec, "Gene duplications on the species tree", "Duplications with support ≥ 0.5 at each node (the circle's area); at a species, those within it");
    if (NS > 30) c.querySelector(".chart").classList.add("scroll");
    const draw = () => speciesTree(c.querySelector(".chart"), D.species_tree, D.node_duplications);
    draw(); onResize["duplications-tree"] = draw;
    setTable(c, [["Species tree node", 0], ["Duplications", 1], ["With support ≥ 0.5", 1]],
      Object.entries(D.node_duplications).map(([k, v]) => [esc(k), v[0], v[1]]));
  } else sec.insertAdjacentHTML("beforeend", `<div class="empty">No duplications per species tree node in these results (Duplications_per_Species_Tree_Node.tsv).</div>`);
  const bins = [[0, 0], [1, 1], [2, 2], [3, 3], [4, 5], [6, 10], [11, 20], [21, 50], [51, Infinity]].map(([lo, hi]) =>
    [lo === hi ? String(lo) : hi === Infinity ? `>${lo - 1}` : `${lo}-${hi}`, O.dups05.filter(d => d >= lo && d <= hi).length]);
  const c2 = card(sec, "Duplications per orthogroup", "Orthogroups by their number of gene duplications with support ≥ 0.5");
  const draw2 = () => columns(c2.querySelector(".chart"), bins, {unit: "orthogroups", xTitle: "Duplications"});
  draw2(); onResize["duplications-hist"] = draw2;
  setTable(c2, [["Duplications", 0], ["Orthogroups", 1]], bins.map(([l, v]) => [esc(l), v]));
  const top = [...OG_ROWS].sort((a, b) => O.dups05[b] - O.dups05[a]).slice(0, 15).filter(r => O.dups05[r]);
  if (top.length) {
    const c3 = card(sec, "Orthogroups with the most duplications", "Duplications with support ≥ 0.5");
    const draw3 = () => hbars(c3.querySelector(".chart"), top.map(r => ({label: OGIDS[r], k: -1, value: O.dups05[r], tl: `${OGIDS[r]}: ${fmt(O.genes[r])} genes, ${fmt(O.dups[r])} duplications in all`})), {unit: "duplications"});
    draw3(); onResize["duplications-top"] = draw3;
    setTable(c3, [["Orthogroup", 0], ["Genes", 1], ["Duplications", 1], ["With support ≥ 0.5", 1]], top.map(r => [O.tree[r] ? treeLink(OGIDS[r]) : esc(OGIDS[r]), O.genes[r], O.dups[r], O.dups05[r]]));
  }
};

// HOGs: those of each species tree node in a data file (HOGs_<node>, in columns), the genes of each species
// in them in chunks (HOGCounts_<node>_<c>)
let hogWanted = null;
renderers.hogs = () => {
  const sec = $("tab-hogs");
  if (!D.hogs.length) { sec.innerHTML = `<div class="empty">No hierarchical orthogroups (HOGs) in these results (Phylogenetic_Hierarchical_Orthogroups).</div>`; return; }
  const c = card(sec, "HOGs at each species tree node", "Hierarchical orthogroups (Phylogenetic_Hierarchical_Orthogroups/<node>.tsv)");
  if (D.hogs.length > 30) c.querySelector(".chart").classList.add("scroll");
  const draw = () => hbars(c.querySelector(".chart"), D.hogs.map(h => ({label: h.node, k: -1, value: h.hogs, tl: `${h.node}: ${fmt(h.genes)} genes, ${fmt(h.multi_species)} HOGs with more than one species`})), {unit: "HOGs"});
  draw(); onResize["hogs-nodes"] = draw;
  setTable(c, [["Node", 0], ["HOGs", 1], ["Genes", 1], ["With more than one species", 1]], D.hogs.map(h => [esc(h.node), h.hogs, h.genes, h.multi_species]));
  sec.insertAdjacentHTML("beforeend", `<div class="tools"><label class="muted" for="hog-node">Species tree node</label>`
    + `<select id="hog-node">${D.hogs.map(h => `<option>${esc(h.node)}</option>`).join("")}</select></div><div id="hog-node-view"></div>`);
  $("hog-node").onchange = () => loadHogs($("hog-node").value);
  loadHogs(D.hogs[0].node);
};
function loadHogs(node) {
  hogWanted = node;
  const h = D.hogs.find(x => x.node === node), v = $("hog-node-view");
  v.innerHTML = `<div class="tiles"><div class="tile"><div class="label">HOGs at ${esc(node)}</div><div class="value">${compact(h.hogs)}</div></div>`
    + `<div class="tile"><div class="label">Genes in them</div><div class="value">${compact(h.genes)}</div></div>`
    + `<div class="tile"><div class="label">With more than one species</div><div class="value">${compact(h.multi_species)}</div><div class="note">${fmt(100 * h.multi_species / (h.hogs || 1))}% of the HOGs</div></div></div>`
    + `<div class="grid2" id="hog-grid"></div><div class="tools"><input id="hog-q" placeholder="Search a HOG or orthogroup ID" autocomplete="off"></div>`
    + `<div class="card"><div id="hog-table"><div class="empty">Loading the HOGs of ${esc(node)}...</div></div></div>`;
  const c = card($("hog-grid"), `HOG sizes at ${node}`, "HOGs by their number of genes");
  card($("hog-grid"), `Species per HOG at ${node}`, "HOGs by the number of species with genes in them");
  const draw = () => columns(c.querySelector(".chart"), h.sizes, {unit: "HOGs", xTitle: "Genes"});
  draw(); onResize["hogs-sizes"] = draw;
  setTable(c, [["Genes", 0], ["HOGs", 1]], h.sizes.map(([l, n]) => [esc(l), n]));
  loadData(`HOGs_${node}`).then(H => { if (hogWanted === node) showHogs(node, H); },
    err => { if (hogWanted === node) $("hog-table").innerHTML = `<div class="empty">The HOGs of ${esc(node)} ${esc(err.message)}.</div>`; });
}
function showHogs(node, H) {
  const ids = expandIds(H.ids), rows = ids.map((_, r) => r);
  const ogOf = r => typeof H.ogs[r] === "number" ? OGIDS[H.ogs[r]] : H.ogs[r], treeOf = r => typeof H.ogs[r] === "number" && O.tree[H.ogs[r]];
  const bySp = Array.from({length: NS}, (_, k) => [String(k + 1), 0]);
  for (const n of H.species) if (n) bySp[n - 1][1]++;
  const c2 = $("hog-grid").lastElementChild;
  const draw = () => columns(c2.querySelector(".chart"), bySp, {unit: "HOGs", xTitle: "Species"});
  draw(); onResize["hogs-species"] = draw;
  setTable(c2, [["Species", 0], ["HOGs", 1]], bySp);
  let counts = new Map();
  const cols = [
    {label: "HOG", key: r => ids[r], html: r => esc(ids[r])},
    {label: "Orthogroup", key: r => ogOf(r), html: r => ogLinks(ogOf(r), treeOf(r))},
    {label: "Gene tree clade", key: r => H.clades[r], html: r => esc(H.clades[r])},
    {label: "Genes", key: r => H.genes[r], html: r => fmt(H.genes[r]), num: true},
    {label: "Species", key: r => H.species[r], html: r => fmt(H.species[r]), num: true},
    {label: "Genes per species", head: stripHead(), html: r => stripCanvas().replace('class="strip"', `class="strip" data-r="${r}"`)},
  ];
  const table = dataTable($("hog-table"), {cols, rows,
    prepare: view => chunkRows(`HOGCounts_${node}`, H.chunks, view).then(m => { counts = m; }),
    after: el => drawStrips(el, r => counts.get(r), H.max_genes)});
  onResize["hogs-table"] = () => table.redraw();
  $("hog-q").oninput = () => { const q = $("hog-q").value.trim().toLowerCase(); table.setRows(q ? rows.filter(r => ids[r].toLowerCase().includes(q) || ogOf(r).toLowerCase().includes(q)) : rows); };
}

// redraw the plots for the window's width and the colour scheme
let resizeTimer = null;
const redrawAll = () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(() => Object.values(onResize).forEach(f => f()), 150); };
window.addEventListener("resize", redrawAll);
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll);
window.addEventListener("orthofinder-theme", redrawAll);  // Colours: Light, Dark or Auto
function goToHash() {  // Statistics.html#species, #orthogroups/OG0000003
  const [tab, og] = decodeURIComponent(location.hash.slice(1)).split("/");
  openTab(tab || "overview");
  if (tab === "orthogroups" && og) { $("og-q").value = og; $("og-q").dispatchEvent(new Event("input")); }
}
goToHash();
window.addEventListener("hashchange", goToHash);  // e.g. sent here from another page in OrthoFinder_Results.html
</script></body></html>
"""


def main(args=None):
    """python -m orthofinder.orthology.stats_page <results directory>: write the statistics page."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.orthology.stats_page",
        description="Write Comparative_Genomics_Statistics/Statistics.html for OrthoFinder results: "
                    "the statistics of the run, with plots, to look at in a web browser.")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    results_dir = parser.parse_args(args).results_dir
    if not os.path.isdir(results_dir):
        parser.error(f"not a directory: {results_dir}")
    page = write_stats_page(results_dir)
    if page is None:
        print(f"No orthogroups in {results_dir}")
        return 1
    print(f"Statistics page: {page}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
