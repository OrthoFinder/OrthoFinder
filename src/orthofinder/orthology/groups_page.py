"""
--groups-view: a page to look at the genes of the orthogroups and HOGs in a web
browser, a tab for each:

    Orthogroups/Groups.html
    Orthogroups/Groups_Data/<table>.js

Made from the results files (Orthogroups.tsv, Phylogenetic_Hierarchical_
Orthogroups/N<i>.tsv), so it can be written again for existing results:

    python -m orthofinder.orthology.groups_page <results directory>

The page holds only a summary; each table (the orthogroups, the HOGs of a
species tree node) is a script the page
loads when it is opened, so that the page stays small for large runs (a page
opened from disk may load scripts, but not read other files).
"""

import json
import os
import re

from .. import __version__
from ..phylo.tree_pages import PALETTE
from ..utils import view_pages

DATA_DIR = "Groups_Data"


def split_genes(cell):
    return cell.split(", ") if cell else []


def read_table(filename):
    """The header and rows of a tab-separated file, split only on tabs (cells may be very long)."""
    with open(filename, encoding="utf-8", errors="replace") as infile:
        rows = [line.rstrip("\r\n").split("\t") for line in infile]
    return (rows[0], [r for r in rows[1:] if r and r[0]]) if rows else ([], [])


def node_number(name):
    m = re.match(r"N(\d+)$", name)
    return int(m.group(1)) if m else 10 ** 9


def write_table(data_dir, key, columns, rows):
    """A table for the page, as a script that gives it to the page when loaded."""
    with open(os.path.join(data_dir, key + ".js"), "w", encoding="utf-8") as out:
        table = {"key": key, "columns": columns, "rows": rows}
        out.write(f"tableLoaded({json.dumps(table, separators=(',', ':'))});\n")


def write_groups_page(results_dir, written=None):
    """
    Write the groups page of a results directory and its tables. written: the
    pages a run is writing (view_pages.present).
    Returns the page's filename, or None if the results have no orthogroups.
    """
    og_dir = os.path.join(results_dir, "Orthogroups")
    og_fn = os.path.join(og_dir, "Orthogroups.tsv")
    if not os.path.exists(og_fn):
        return None
    header, og_rows = read_table(og_fn)
    species = header[1:]
    n = len(species)
    data_dir = os.path.join(og_dir, DATA_DIR)
    os.makedirs(data_dir, exist_ok=True)
    tabs = []

    def genes_column(k):
        return {"name": species[k], "kind": "genes", "species": k}

    # the orthogroups
    rows = [[r[0], sum(len(split_genes(c)) for c in r[1:n + 1])] + [split_genes(c) for c in (r[1:n + 1] + [""] * n)[:n]]
            for r in og_rows]
    write_table(data_dir, "Orthogroups", [{"name": "Orthogroup", "kind": "og"}, {"name": "Genes", "kind": "num"}]
                + [genes_column(k) for k in range(n)], rows)
    tabs.append({"id": "orthogroups", "name": "Orthogroups", "note": "The orthogroups of Orthogroups.tsv.",
                 "views": [{"name": "All orthogroups", "key": "Orthogroups", "rows": len(rows)}]})

    # the HOGs, a table for each species tree node (the species of its clade)
    hog_dir = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    views = []
    hog_files = sorted((f for f in os.listdir(hog_dir) if f.endswith(".tsv")) if os.path.isdir(hog_dir) else [],
                       key=lambda f: (node_number(f[:-4]), f))
    for fn in hog_files:
        node = fn[:-4]
        head, hog_rows = read_table(os.path.join(hog_dir, fn))
        cols = [(i, species.index(sp)) for i, sp in enumerate(head) if i >= 3 and sp in species]
        clade = [(i, k) for i, k in cols if any(i < len(r) and r[i] for r in hog_rows)]
        rows = [[r[0], r[1], r[2], sum(len(split_genes(r[i])) for i, _ in clade if i < len(r))]
                + [split_genes(r[i]) if i < len(r) else [] for i, _ in clade] for r in hog_rows]
        write_table(data_dir, f"HOGs_{node}", [{"name": "HOG", "kind": "id"}, {"name": "Orthogroup", "kind": "og"},
                                               {"name": "Gene tree clade", "kind": "text"}, {"name": "Genes", "kind": "num"}]
                    + [genes_column(k) for _, k in clade], rows)
        views.append({"name": f"{node}: {', '.join(species[k] for _, k in clade)}", "key": f"HOGs_{node}", "rows": len(rows)})
    if views:
        tabs.append({"id": "hogs", "name": "HOGs", "views": views,
                     "note": "The hierarchical orthogroups of Phylogenetic_Hierarchical_Orthogroups/N<i>.tsv, "
                             "one species tree node at a time (the species of its clade)."})

    tree_page = view_pages.link(results_dir, "trees", written)
    tree_ogs = []
    trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    if tree_page and os.path.exists(trees_fn):
        with open(trees_fn, encoding="utf-8", errors="replace") as infile:
            tree_ogs = [line.split(":", 1)[0] for line in infile if ":" in line]
    data = {"version": __version__, "data": DATA_DIR, "tabs": tabs,
            "species": [[sp, PALETTE[k % len(PALETTE)]] for k, sp in enumerate(species)],
            "nav": view_pages.nav(results_dir, "groups", written),
            "tree_page": tree_page if tree_ogs else None, "tree_ogs": tree_ogs,
            "stats_page": view_pages.link(results_dir, "stats", written)}
    page = os.path.join(og_dir, "Groups.html")
    with open(page, "w", encoding="utf-8") as out:
        out.write(HTML.replace("__NAV_CSS__", view_pages.NAV_CSS).replace("__NAV_HTML__", view_pages.NAV_HTML)
                  .replace("__NAV_JS__", view_pages.NAV_JS).replace("__THEME_JS__", view_pages.THEME_JS)
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
__NAV_JS__
const PAGE = 100, SHOW = 6;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const TREE_OGS = new Set(D.tree_ogs);
$("sub").textContent = `OrthoFinder ${D.version}. The genes of each orthogroup and HOG, from the results files. `
  + `Each table is loaded when it is opened.`;
$("legend").innerHTML = D.species.map(([s, c]) => `<span><span class="sw" style="background:${c}"></span>${esc(s)}</span>`).join("");
$("tabs").innerHTML = D.tabs.map(t => `<button role="tab" data-tab="${t.id}">${esc(t.name)}</button>`).join("");

// tables, each loaded once, when first opened
const tables = {}, waiting = {};
window.tableLoaded = t => { tables[t.key] = t; (waiting[t.key] || []).forEach(f => f(t)); delete waiting[t.key]; };
function loadTable(key, then) {
  if (tables[key]) return then(tables[key]);
  if (!waiting[key]) {
    waiting[key] = [];
    const s = document.createElement("script");
    s.src = `${D.data}/${encodeURIComponent(key)}.js`;
    s.onerror = () => { delete waiting[key]; $("table").innerHTML = `<tr><td class="muted">The table ${esc(key)} could not be loaded (${esc(s.src)}). The page needs its folder ${esc(D.data)} next to it: copy (or download) them together, e.g. the whole results directory.</td></tr>`; };
    document.head.appendChild(s);
  }
  waiting[key].push(then);
}

const q = () => $("q").value.trim().toLowerCase();
const mark = s => { const t = q(), i = t ? String(s).toLowerCase().indexOf(t) : -1;
  return i < 0 ? esc(s) : esc(String(s).slice(0, i)) + `<mark>${esc(String(s).slice(i, i + t.length))}</mark>` + esc(String(s).slice(i + t.length)); };
const chips = gs => gs.map(g => `<span class="gene">${mark(g)}</span>`).join("");
function genes(gs) {
  const t = q(), hit = t && gs.slice(SHOW).some(g => g.toLowerCase().includes(t));
  return gs.length <= SHOW ? chips(gs) : chips(gs.slice(0, SHOW))
    + `<details${hit ? " open" : ""}><summary>+${gs.length - SHOW} more</summary>${chips(gs.slice(SHOW))}</details>`;
}
function ogCell(og) {
  const links = [];
  if (D.tree_page && TREE_OGS.has(og)) links.push(`<a href="${esc(D.tree_page)}#${encodeURIComponent(og)}" title="Its gene tree and alignment">tree</a>`);
  if (D.stats_page) links.push(`<a href="${esc(D.stats_page)}#orthogroups/${encodeURIComponent(og)}" title="In the statistics">statistics</a>`);
  return `<td class="id">${mark(og)}${links.length ? `<span class="links">${links.join("")}</span>` : ""}</td>`;
}
const th = c => c.kind === "genes" ? `<th><span class="sw" style="background:${D.species[c.species][1]}"></span>${esc(c.name)}</th>`
  : `<th${c.kind === "num" ? ' class="num"' : ""}>${esc(c.name)}</th>`;
const td = (c, v) => c.kind === "genes" ? `<td>${genes(v)}</td>` : c.kind === "og" ? ogCell(v)
  : c.kind === "num" ? `<td class="num">${Number(v).toLocaleString()}</td>`
  : `<td class="${c.kind === "id" ? "id" : ""}">${mark(v)}</td>`;
const text = r => r.map(v => Array.isArray(v) ? v.join(" ") : String(v)).join(" ").toLowerCase();

let tab = D.tabs[0], table = null, rows = [], page = 0;
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
  $("table").innerHTML = `<tr><td class="muted">Loading ${esc(tab.views[i].name)}...</td></tr>`;
  loadTable(key, t => { if (tab.views[+$("view").value].key === key) { table = t; filter(); } });
}
function filter() {
  if (!table) return;
  const t = q();
  rows = t ? table.rows.filter(r => text(r).includes(t)) : table.rows;
  page = 0; draw();
}
function draw() {
  const start = page * PAGE;
  $("table").innerHTML = "<thead><tr>" + table.columns.map(th).join("") + "</tr></thead><tbody>"
    + rows.slice(start, start + PAGE).map(r => "<tr>" + r.map((v, i) => td(table.columns[i], v)).join("") + "</tr>").join("") + "</tbody>";
  $("count").textContent = rows.length
    ? `${(start + 1).toLocaleString()}–${Math.min(start + PAGE, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}` : "no matches";
  $("prev").disabled = page === 0;
  $("next").disabled = start + PAGE >= rows.length;
}
$("tabs").onclick = e => { const b = e.target.closest("[data-tab]"); if (b) { openTab(b.dataset.tab); history.replaceState(null, "", "#" + b.dataset.tab); } };
$("view").onchange = e => openView(+e.target.value);
$("q").addEventListener("input", filter);
$("prev").onclick = () => { page--; draw(); };
$("next").onclick = () => { page++; draw(); };
// Groups.html#orthogroups/OG0000003: that orthogroup; #hogs: the HOGs
function goToHash() {
  const [tab, search] = decodeURIComponent(location.hash.slice(1)).split("/");
  $("q").value = search || "";
  openTab(tab || D.tabs[0].id);
}
goToHash();
window.addEventListener("hashchange", goToHash);  // e.g. sent here from another page in OrthoFinder_Results.html
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
