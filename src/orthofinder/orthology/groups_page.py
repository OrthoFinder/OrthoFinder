"""
--groups-view: a page to look at the genes of the orthogroups and HOGs in a web
browser, a tab for each:

    Orthogroups/Groups.html
    Orthogroups/Groups_Data/*.js

Made from the results files (Orthogroups.tsv, Phylogenetic_Hierarchical_
Orthogroups/N<i>.tsv), so it can be written again for existing results:

    python -m orthofinder.orthology.groups_page <results directory>

The page holds only a summary. Each table (the orthogroups, the HOGs of a
species tree node) is a data file of its rows (their IDs and sizes), and their
genes are in chunks, loaded for the rows shown; genes are found with an index
of their names (view_pages.write_gene_index). So the page stays small and quick
for thousands of species: a page opened from disk may load scripts, but not
read other files.
"""

import json
import os
import re

from .. import __version__
from ..phylo.tree_pages import PALETTE
from ..utils import view_pages

DATA_DIR = "Groups_Data"


def node_number(name):
    m = re.match(r"N(\d+)$", name)
    return int(m.group(1)) if m else 10 ** 9


def read_groups(filename, species_index, first=1):
    """
    The header and rows of Orthogroups.tsv or a HOG file: [[the first cells
    (up to first), [[species, [genes]], ...]], ...] with the species of the
    header's columns from first on (species_index: {name: number}), read
    without splitting each row at every tab (most cells are empty).
    """
    rows = []
    with open(filename, encoding="utf-8", errors="replace") as infile:
        header = infile.readline().rstrip("\r\n").split("\t")
        column = {k: species_index[sp] for k, sp in enumerate(header) if k >= first and sp in species_index}
        for line in infile:
            fields = line.rstrip("\r\n").split("\t", first)
            if not fields[0] or len(fields) < first:
                continue
            genes = sorted((column[k], text.split(", ")) for k, text in view_pages.sparse_cells(line) if k in column)
            rows.append([fields[:first], genes])
    return header, rows


def table_file(data_dir, key, rows, extra):
    """
    A table for the page: <key> (its rows' IDs, genes and species, extra, and
    the first row of each chunk) and <key>_<c> (the genes of each species in its
    rows, [[species, [genes]], ...]). Returns its number of rows.
    """
    view_pages.write_data(data_dir, key, {
        **extra, "genes": [sum(len(g) for _, g in genes) for _, genes in rows],
        "species": [len(genes) for _, genes in rows],
        "chunks": view_pages.write_chunks(data_dir, key, [genes for _, genes in rows],
                                          lambda genes: sum(len(g) for _, g in genes))})
    return len(rows)


def write_groups_page(results_dir, written=None):
    """
    Write the groups page of a results directory and its data files. written:
    the pages a run is writing (view_pages.present).
    Returns the page's filename, or None if the results have no orthogroups.
    """
    og_dir = os.path.join(results_dir, "Orthogroups")
    og_fn = view_pages.orthogroups_file(results_dir)   # or the MCL orthogroups (e.g. -og)
    if og_fn is None:
        return None
    q_mcl = os.path.basename(og_fn).startswith("MCL_")
    with open(og_fn, encoding="utf-8", errors="replace") as infile:
        species = infile.readline().rstrip("\r\n").split("\t")[1:]
    species_index = {sp: k for k, sp in enumerate(species)}
    _, og_rows = read_groups(og_fn, species_index)
    data_dir = os.path.join(og_dir, DATA_DIR)
    os.makedirs(data_dir, exist_ok=True)
    tabs = []

    # the orthogroups (Orthogroups, which the page loads first: the IDs of the orthogroups, for the HOGs)
    ids = [first[0] for first, _ in og_rows]
    og_row = {og: i for i, og in enumerate(ids)}
    tree_page = view_pages.link(results_dir, "trees", written)
    with_trees = set()
    trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    if tree_page and os.path.exists(trees_fn):
        with open(trees_fn, encoding="utf-8", errors="replace") as infile:
            with_trees = {line[:line.find(":")] for line in infile if ":" in line}
    n = table_file(data_dir, "Orthogroups", og_rows, {"ids": view_pages.compact_ids(ids), "cols": list(range(len(species))),
                                                      "tree": [int(og in with_trees) for og in ids]})
    tabs.append({"id": "orthogroups", "name": "MCL orthogroups" if q_mcl else "Orthogroups",
                 "note": f"The orthogroups of {os.path.basename(og_fn)}.",
                 "views": [{"name": "All orthogroups", "key": "Orthogroups", "rows": n}]})
    gene_index = view_pages.write_gene_index(data_dir, ((gene, i) for i, (_, genes) in enumerate(og_rows)
                                                        for _, gs in genes for gene in gs))
    del og_rows

    # the HOGs, a table for each species tree node (the species of its clade)
    hog_dir = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    views = []
    hog_files = sorted((f for f in os.listdir(hog_dir) if f.endswith(".tsv")) if os.path.isdir(hog_dir) else [],
                       key=lambda f: (node_number(f[:-4]), f))
    for fn in hog_files:
        node = fn[:-4]
        _, hog_rows = read_groups(os.path.join(hog_dir, fn), species_index, first=3)
        clade = sorted({k for _, genes in hog_rows for k, _ in genes})
        n = table_file(data_dir, f"HOGs_{node}", hog_rows, {
            "ids": view_pages.compact_ids([first[0] for first, _ in hog_rows]),
            "ogs": [og_row.get(first[1], first[1]) for first, _ in hog_rows],  # the orthogroup's row (its ID if not one)
            "clades": [first[2] for first, _ in hog_rows], "cols": clade})
        names = [species[k] for k in clade]
        views.append({"name": f"{node}: " + (", ".join(names) if len(names) <= 4 else
                                             f"{len(names)} species ({', '.join(names[:3])}, ...)"),
                      "key": f"HOGs_{node}", "rows": n})
    if views:
        tabs.append({"id": "hogs", "name": "HOGs", "views": views,
                     "note": "The hierarchical orthogroups of Phylogenetic_Hierarchical_Orthogroups/N<i>.tsv, "
                             "one species tree node at a time (the species of its clade)."})

    data = {"version": __version__, "data": DATA_DIR, "tabs": tabs, "genes": gene_index,
            "species": [[sp, PALETTE[k % len(PALETTE)]] for k, sp in enumerate(species)],
            "nav": view_pages.nav(results_dir, "groups", written),
            "tree_page": tree_page if with_trees else None,
            "stats_page": view_pages.link(results_dir, "stats", written)}
    page = os.path.join(og_dir, "Groups.html")
    with open(page, "w", encoding="utf-8") as out:
        out.write(HTML.replace("__NAV_CSS__", view_pages.NAV_CSS).replace("__NAV_HTML__", view_pages.NAV_HTML)
                  .replace("__NAV_JS__", view_pages.NAV_JS).replace("__THEME_JS__", view_pages.THEME_JS)
                  .replace("__DATA_JS__", view_pages.DATA_JS)
                  .replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
    return page


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Orthogroups and HOGs</title>
__THEME_JS__
<style>
:root { color-scheme:light; --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e8; --head:#f5f5f7; --accent:#2f6fde; --chip:#eef3fd; --mark:#fff1a8; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff; --chip:#1e2a40; --mark:#5c4d00; } }
:root[data-theme="dark"] { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff; --chip:#1e2a40; --mark:#5c4d00; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width:1400px; margin:0 auto; padding:24px clamp(16px, 2vw, 32px) 48px; }
h1 { font-size:32px; margin:0 0 6px; } h2 { font-size:19px; margin:20px 0 8px; }
.sub, .muted { color:var(--muted); font-size:12px; }
.legend { display:flex; flex-wrap:wrap; gap:6px 14px; margin:12px 0; font-size:13px; }
.sw { display:inline-block; width:10px; height:10px; border-radius:2px; margin-right:6px; vertical-align:0; }
.tabs { display:flex; gap:4px; flex-wrap:wrap; margin:16px 0 4px; border-bottom:1px solid var(--line); }
.tabs button { border:0; border-radius:0; background:none; padding:8px 12px; color:var(--muted); border-bottom:2px solid transparent; margin-bottom:-1px; }
.tabs button[aria-selected="true"] { color:var(--fg); border-bottom-color:var(--accent); font-weight:600; }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin:10px 0 8px; }
input, select { font:inherit; padding:7px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
input { width:min(380px,100%); } select { max-width:100%; }
button { font:inherit; padding:6px 12px; border:1px solid var(--line); border-radius:6px; background:var(--head); color:var(--fg); cursor:pointer; }
button:disabled { opacity:.4; cursor:default; }
.wrap { overflow:auto; border:1px solid var(--line); border-radius:8px; max-height:75vh; }
table { border-collapse:collapse; width:100%; font-size:13px; }
th, td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--line); vertical-align:top; }
th { background:var(--head); font-weight:600; white-space:nowrap; position:sticky; top:0; z-index:1; }
td.num, th.num { font-variant-numeric:tabular-nums; text-align:right; }
.id { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; white-space:nowrap; }
.links { display:block; font-size:11px; } .links a { margin-right:8px; }
.gene { display:inline-block; background:var(--chip); border-radius:4px; padding:1px 6px; margin:1px 3px 1px 0;
        font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:12px; }
details { display:inline; } summary { display:inline; cursor:pointer; color:var(--accent); font-size:12px; white-space:nowrap; }
a { color:var(--accent); }
mark { background:var(--mark); color:inherit; }
.spg { margin:1px 0 3px; } .spn { font-size:12px; color:var(--muted); margin-right:4px; }
.legend > details { display:block; width:100%; } .legend > details > summary { display:list-item; color:var(--muted); font-size:13px; }
.legend details .legend { max-height:30vh; overflow:auto; margin-top:6px; }
__NAV_CSS__
</style></head><body><main>
__NAV_HTML__
<h1>Orthogroups and HOGs</h1>
<div class="sub" id="sub"></div>
<div class="legend" id="legend"></div>
<nav class="tabs" id="tabs" role="tablist"></nav>
<p class="muted" id="note"></p>
<div class="tools">
  <select id="view"></select>
  <input id="q" placeholder="Search an ID or gene name" autocomplete="off">
  <button id="prev">Previous</button><button id="next">Next</button><span class="muted" id="count"></span>
</div>
<div class="wrap"><table id="table"></table></div>
</main>
<script>
const D = __DATA__;
__DATA_JS__
__NAV_JS__
const PAGE = 100, SHOW = 6, SHOW_SPECIES = 6, COLUMNS_MAX = 20;  // a column for each species up to COLUMNS_MAX of them
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const NS = D.species.length;
$("sub").textContent = `OrthoFinder ${D.version}. The genes of each orthogroup and HOG, from the results files. `
  + `Each table is loaded when it is opened, and the genes of the rows shown.`;
{
  const key = D.species.map(([s, c]) => `<span><span class="sw" style="background:${c}"></span>${esc(s)}</span>`).join("");
  $("legend").innerHTML = NS > 40 ? `<details><summary>The species and their colours (${NS.toLocaleString()})</summary><div class="legend">${key}</div></details>` : key;
}
$("tabs").innerHTML = D.tabs.map(t => `<button role="tab" data-tab="${t.id}">${esc(t.name)}</button>`).join("");
const message = text => { $("table").innerHTML = `<tr><td class="muted">${esc(text)}</td></tr>`; };

const q = () => $("q").value.trim().toLowerCase();
const mark = s => { const t = q(), i = t ? String(s).toLowerCase().indexOf(t) : -1;
  return i < 0 ? esc(s) : esc(String(s).slice(0, i)) + `<mark>${esc(String(s).slice(i, i + t.length))}</mark>` + esc(String(s).slice(i + t.length)); };
const chips = gs => gs.map(g => `<span class="gene">${mark(g)}</span>`).join("");
const hasHit = gs => { const t = q(); return t && gs.some(g => g.toLowerCase().includes(t)); };
// the hidden part of a long list: made when opened (a large orthogroup can have thousands of genes), or at
// once, open, if the search matches in it
let more = [];
function hidden(summary, html, hit) {
  if (hit) return `<details open><summary>${summary}</summary>${html()}</details>`;
  more.push(html);
  return `<details data-more="${more.length - 1}"><summary>${summary}</summary></details>`;
}
document.addEventListener("toggle", e => {
  const d = e.target;
  if (d.open && d.dataset && d.dataset.more !== undefined && more[+d.dataset.more]) {
    d.insertAdjacentHTML("beforeend", more[+d.dataset.more]()); delete d.dataset.more;
  }
}, true);
function genes(gs) {
  if (gs.length <= SHOW) return chips(gs);
  const rest = gs.slice(SHOW);
  return chips(gs.slice(0, SHOW)) + hidden(`+${rest.length} more`, () => chips(rest), hasHit(rest));
}
// the genes of many species in one cell: each species' genes after its name; past SHOW_SPECIES species, the rest to open
function bySpecies(groups) {
  if (q()) groups = [...groups.filter(([, gs]) => hasHit(gs)), ...groups.filter(([, gs]) => !hasHit(gs))];  // those matching the search first
  const one = ([k, gs]) => `<div class="spg"><span class="sw" style="background:${D.species[k][1]}"></span><span class="spn">${esc(D.species[k][0])}</span> ${genes(gs)}</div>`;
  const rest = groups.slice(SHOW_SPECIES);
  return groups.slice(0, SHOW_SPECIES).map(one).join("")
    + (rest.length ? hidden(`+${rest.length} more species`, () => rest.map(one).join(""), rest.some(([, gs]) => hasHit(gs))) : "");
}
function ogCell(og, hasTree) {
  const links = [];
  if (D.tree_page && hasTree) links.push(`<a href="${esc(D.tree_page)}#${encodeURIComponent(og)}" title="Its gene tree and alignment">tree</a>`);
  if (D.stats_page) links.push(`<a href="${esc(D.stats_page)}#orthogroups/${encodeURIComponent(og)}" title="In the statistics">statistics</a>`);
  return `<td class="id">${mark(og)}${links.length ? `<span class="links">${links.join("")}</span>` : ""}</td>`;
}
const num = x => Number(x).toLocaleString();

// a table (a data file): its rows' IDs (and for HOGs their orthogroups and gene tree clades), genes and species,
// the species of its columns (cols), and its genes in chunks, loaded for the rows shown
let OG = null, OGIDS = null;  // the orthogroups, loaded first
function prepare(t) {
  if (t.ready) return t;
  t.ids = expandIds(t.ids);
  t.rowIds = t.ogs ? t.ogs.map(x => typeof x === "number" ? OGIDS[x] : x) : t.ids;  // its orthogroup's ID
  t.columns = t.cols.length <= COLUMNS_MAX;  // a column for each species, else one for all
  t.ready = true;
  return t;
}
function headHtml(t) {
  const sp = k => `<th><span class="sw" style="background:${D.species[k][1]}"></span>${esc(D.species[k][0])}</th>`;
  return "<thead><tr>" + (t.ogs ? "<th>HOG</th><th>Orthogroup</th><th>Gene tree clade</th>" : "<th>Orthogroup</th>")
    + `<th class="num">Genes</th>` + (t.columns ? t.cols.map(sp).join("") : `<th>Genes by species <span class="muted">(${num(t.cols.length)} species)</span></th>`) + "</tr></thead>";
}
function rowHtml(t, r, groups) {
  const og = t.ogs && typeof t.ogs[r] === "number" ? t.ogs[r] : (t.ogs ? -1 : r), hasTree = og >= 0 && OG.tree[og];
  let cells;
  if (t.columns) { const of = new Map(groups); cells = t.cols.map(k => `<td>${genes(of.get(k) || [])}</td>`).join(""); }
  else cells = `<td>${bySpecies(groups)}</td>`;
  return "<tr>" + (t.ogs ? `<td class="id">${mark(t.ids[r])}</td>${ogCell(t.rowIds[r], hasTree)}<td>${mark(t.clades[r])}</td>` : ogCell(t.ids[r], hasTree))
    + `<td class="num">${num(t.genes[r])}<div class="muted" style="white-space:nowrap">${num(t.species[r])} species</div></td>${cells}</tr>`;
}

let tab = D.tabs[0], table = null, rows = [], page = 0, seq = 0;
function openTab(id, view = 0) {
  tab = D.tabs.find(t => t.id === id) || D.tabs[0];
  for (const b of document.querySelectorAll("#tabs button")) b.setAttribute("aria-selected", b.dataset.tab === tab.id);
  $("note").textContent = tab.note;
  $("view").innerHTML = tab.views.map((v, i) => `<option value="${i}">${esc(v.name)} (${v.rows.toLocaleString()})</option>`).join("");
  $("view").style.display = tab.views.length > 1 ? "" : "none";
  $("view").value = String(view);
  openView(view);
}
function openView(i) {
  const key = tab.views[i].key;
  table = null;
  message(`Loading ${tab.views[i].name}...`);
  loadData(key).then(t => { if (tab.views[+$("view").value].key === key) { table = prepare(t); table.key = key; filter(); } },
    err => message(`The table ${key} ${err.message}.`));
}
// the rows that match the search: by their IDs, and by their genes (the start of a gene's name, with the index)
async function filter() {
  if (!table) return;
  const my = ++seq, t = q(), T = table;
  if (!t) { rows = T.ids.map((_, r) => r); page = 0; return draw(); }
  const hits = new Set();
  T.ids.forEach((id, r) => { if (id.toLowerCase().includes(t) || (T.ogs && T.rowIds[r].toLowerCase().includes(t))) hits.add(r); });
  $("count").textContent = "searching...";
  try {
    const found = await findGenes(D.genes, t, 5000);  // [[gene, orthogroup's row]]
    if (my !== seq) return;
    const ogs = new Set(found.map(g => g[1]));
    if (!T.ogs) ogs.forEach(r => hits.add(r));
    else {  // the HOGs of those orthogroups with the gene (an orthogroup can have several HOGs at a node)
      const cand = T.ogs.map((x, r) => [x, r]).filter(([x]) => typeof x === "number" && ogs.has(x)).map(([, r]) => r);
      const of = await chunkRows(T.key, T.chunks, cand);
      if (my !== seq) return;
      for (const r of cand) if (of.get(r).some(([, gs]) => gs.some(g => g.toLowerCase().startsWith(t)))) hits.add(r);
    }
  } catch (e) { if (my === seq) message(`The index of the genes ${e.message}.`); return; }
  rows = [...hits].sort((a, b) => a - b);
  page = 0; draw();
}
async function draw() {
  const my = ++seq, T = table, start = page * PAGE, shown = rows.slice(start, start + PAGE);
  let of;
  try { of = await chunkRows(T.key, T.chunks, shown); } catch (e) { if (my === seq) message(`The genes ${e.message}.`); return; }
  if (my !== seq) return;
  more = [];
  $("table").innerHTML = headHtml(T) + "<tbody>" + shown.map(r => rowHtml(T, r, of.get(r))).join("") + "</tbody>";
  $("count").textContent = rows.length
    ? `${(start + 1).toLocaleString()}–${Math.min(start + PAGE, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}` : "no matches";
  $("prev").disabled = page === 0;
  $("next").disabled = start + PAGE >= rows.length;
}
$("tabs").onclick = e => { const b = e.target.closest("[data-tab]"); if (b) { openTab(b.dataset.tab); history.replaceState(null, "", "#" + b.dataset.tab); } };
$("view").onchange = e => openView(+e.target.value);
let typing = null;
$("q").addEventListener("input", () => { clearTimeout(typing); typing = setTimeout(filter, 200); });
$("prev").onclick = () => { page--; draw(); };
$("next").onclick = () => { page++; draw(); };
// Groups.html#orthogroups/OG0000003: that orthogroup; #hogs: the HOGs
function goToHash() {
  const [tab, search] = decodeURIComponent(location.hash.slice(1)).split("/");
  $("q").value = search || "";
  openTab(tab || D.tabs[0].id);
}
message("Loading the orthogroups...");
loadData("Orthogroups").then(t => {
  OG = t; OGIDS = expandIds(t.ids);
  goToHash();
  window.addEventListener("hashchange", goToHash);  // e.g. sent here from another page in OrthoFinder_Results.html
}, err => message(`The orthogroups ${err.message}.`));
</script></body></html>
"""


def main(args=None):
    """python -m orthofinder.orthology.groups_page <results directory>: write the groups page."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.orthology.groups_page",
        description="Write Orthogroups/Groups.html for OrthoFinder results: the genes of each "
                    "orthogroup and HOG, to look at in a web browser.")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    a = parser.parse_args(args)
    if not os.path.isdir(a.results_dir):
        parser.error(f"not a directory: {a.results_dir}")
    page = write_groups_page(a.results_dir)
    if page is None:
        print(f"No orthogroups in {a.results_dir}")
        return 1
    print(f"Orthogroups and HOGs page: {page}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
