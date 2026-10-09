"""
With -x: OrthoFinder's results in OrthoXML format (http://orthoxml.org,
version 0.3), for other tools, and as pages to look at in a web browser:

    Orthogroups/Orthogroups.orthoxml, Orthogroups.html
    Phylogenetic_Hierarchical_Orthogroups/N<i>.orthoxml (one per N<i>.tsv), HOGs.html
    Orthologs/Orthologs.orthoxml, Orthologs.html (only with -pof: with many species they are very large)

A gene has the same id in all the OrthoXML files. The species information
comes from the input FASTA files. OrthoXML requires an NCBI taxonomy ID and a
database for each species, which OrthoFinder does not know: the taxonomy ID
is 0 (unknown) and the database is the FASTA file.
"""

import csv
import glob
import gzip
import json
import os
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape, quoteattr

from .. import __version__

NAMESPACE = "http://orthoXML.org/2011/"
SCHEMA_LOCATION = "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd"
UNKNOWN_TAX_ID = "0"
UNKNOWN_VERSION = "NA"
GENE_SEP = ", "


def species_for_analysis(species_ids_fn, species_to_use):
    """
    (species name, FASTA filename) for each species in species_to_use, from
    SpeciesIDs.txt ("<ID>: <FASTA filename>" lines). The species name is the
    FASTA filename without its extension, as in the other results files.
    """
    fasta_files = {}
    with open(species_ids_fn) as infile:
        for line in infile:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            i_species, fasta_filename = line.split(": ", 1)
            fasta_files[int(i_species)] = fasta_filename
    return [(fasta_files[i].rsplit(".", 1)[0], fasta_files[i]) for i in species_to_use]


# --- Orthogroups ---------------------------------------------------------------

def write_results(species, ogs, n_seqs_per_species, id_dict, species_to_use, results_base_filename):
    """
    Write <results_base_filename>.orthoxml and <results_base_filename>.html.

    species: from species_for_analysis, one entry per species in species_to_use
    ogs: the orthogroups, sets of sequence IDs ("<species>_<sequence>"); orthogroup
         i is OG<i> as in Orthogroups.tsv. Single-gene orthogroups are not groups.
    n_seqs_per_species: {species ID: number of sequences}
    id_dict: {sequence ID: gene name}
    """
    # every gene of every species, numbered across species in this order
    gene_number, species_genes = {}, []
    for i_species, (name, fasta_filename) in zip(species_to_use, species):
        genes = []
        for i_seq in range(n_seqs_per_species[i_species]):
            seq_id = f"{i_species}_{i_seq}"
            gene_number[seq_id] = len(gene_number)
            genes.append((gene_number[seq_id], id_dict[seq_id]))
        species_genes.append((name, fasta_filename, genes))
    multi = [(i_og, og) for i_og, og in enumerate(ogs) if len(og) >= 2]
    groups = ((f"OG{i_og:07d}", sorted(gene_number[s] for s in og), ()) for i_og, og in multi)
    write_orthoxml(results_base_filename + ".orthoxml", species_genes, groups)

    column = {i_species: i for i, i_species in enumerate(species_to_use)}
    rows = []
    n_assigned = [0] * len(species)
    for i_og, og in multi:
        per_species = [[] for _ in species]
        for seq_id in sorted(og, key=lambda s: tuple(map(int, s.split("_")))):
            i_col = column[int(seq_id.split("_")[0])]
            per_species[i_col].append(id_dict[seq_id])
            n_assigned[i_col] += 1
        rows.append([f"OG{i_og:07d}", len(og)] + per_species)
    n_genes = [n_seqs_per_species[i] for i in species_to_use]
    write_page(
        results_base_filename + ".html",
        title="Orthogroups",
        note="The orthogroups of Orthogroups.tsv (also in Orthogroups.orthoxml).",
        stats=[(len(species), "species"), (len(rows), "orthogroups"), (sum(n_genes), "genes"),
               (percent(sum(n_assigned), sum(n_genes)), "genes in orthogroups")],
        info={"title": "Species", "columns": ["Species", "FASTA file", "#Genes", "#In orthogroups"],
              "rows": [[name, fasta, n, percent(a, n)] for (name, fasta), n, a in zip(species, n_genes, n_assigned)]},
        views=[{"name": "Orthogroups",
                "columns": ["@Orthogroup", "#Genes"] + ["*" + name for name, _ in species],
                "rows": rows}],
    )


# --- HOGs and orthologs, from the results files ----------------------------------

def write_hogs_and_orthologs(orthogroups_orthoxml, hog_dir, orthologs_dir):
    """
    From the final results files: an OrthoXML file for each
    Phylogenetic_Hierarchical_Orthogroups/N<i>.tsv and one for the orthologs
    (Orthologs/<species>.tsv[.gz]), and a page for each. The species and
    the gene ids are those of Orthogroups.orthoxml. Returns the files written.
    """
    species, gene_ids = read_genes(orthogroups_orthoxml)
    written = []
    if hog_dir and os.path.isdir(hog_dir):
        written += write_hogs(species, gene_ids, hog_dir)
    if orthologs_dir and os.path.isdir(orthologs_dir):
        written += write_orthologs(species, gene_ids, orthologs_dir)
    return written


def read_genes(orthoxml_fn):
    """[(species name, FASTA file)], {species name: {gene name: gene id}} of an OrthoXML file."""
    species, gene_ids = [], {}
    for _, el in ET.iterparse(orthoxml_fn):
        if el.tag == f"{{{NAMESPACE}}}species":
            db = el.find(f"{{{NAMESPACE}}}database")
            name = el.get("name")
            species.append((name, db.get("name")))
            gene_ids[name] = {g.get("geneId"): int(g.get("id")) for g in db.iter(f"{{{NAMESPACE}}}gene")}
            el.clear()
    return species, gene_ids


def split_genes(cell):
    return cell.split(GENE_SEP) if cell else []


def node_number(fn):
    n = os.path.basename(fn)[1:-4]
    return int(n) if n.isdigit() else -1


def write_hogs(species, gene_ids, hog_dir):
    written, views = [], []
    for fn in sorted(glob.glob(os.path.join(hog_dir, "N*.tsv")), key=node_number):
        node = os.path.basename(fn)[:-4]
        with open(fn, newline="") as infile:
            table = list(csv.reader(infile, delimiter="\t"))
        header, rows = table[0], table[1:]
        if not rows:  # OrthoXML needs at least one group
            continue
        sp_cols = header[3:]
        # the species of this node's clade: those with genes in its HOGs
        clade = [j for j in range(len(sp_cols)) if any(row[3 + j] for row in rows)]
        clade_names = [sp_cols[j] for j in clade]
        groups, page_rows = [], []
        for row in rows:
            per_species = [split_genes(row[3 + j]) for j in clade]
            genes = sorted(gene_ids[sp_cols[j]][g] for j, gs in zip(clade, per_species) for g in gs)
            if len(genes) >= 2:  # an OrthoXML group has at least two genes
                groups.append((row[0], genes, (("orthogroup", row[1]), ("gene tree parent clade", row[2]))))
            page_rows.append([row[0], row[1], row[2], len(genes)] + per_species)
        if groups:
            out = fn[:-4] + ".orthoxml"
            write_orthoxml(
                out, genes_of(species, gene_ids, groups, clade_names), groups,
                notes=f"Hierarchical orthogroups (HOGs) at species tree node {node} "
                      f"({', '.join(clade_names)}), as in {os.path.basename(fn)}, "
                      f"except the HOGs of a single gene.")
            written.append(out)
        views.append({"name": f"{node}: {', '.join(clade_names)}",
                      "columns": ["@HOG", "Orthogroup", "Gene tree parent clade", "#Genes"]
                                 + ["*" + name for name in clade_names],
                      "rows": page_rows})
    if views:
        out = os.path.join(hog_dir, "HOGs.html")
        write_page(
            out, title="Hierarchical orthogroups (HOGs)",
            note="The HOGs of the N<i>.tsv files, one species tree node at a time (also in N<i>.orthoxml, "
                 "without the HOGs of a single gene). The orthogroups are the HOGs at the root, N0.",
            stats=[(len(views), "species tree nodes"), (sum(len(v["rows"]) for v in views), "HOGs")],
            views=views)
        written.append(out)
    return written


def write_orthologs(species, gene_ids, orthologs_dir):
    names = [name for name, _ in species]
    index = {name: i for i, name in enumerate(names)}
    pairs = {}  # (i, j), i < j -> rows [orthogroup, genes of i, genes of j]
    for i, name in enumerate(names):
        fn = os.path.join(orthologs_dir, name + ".tsv.gz")  # the results file, or with -pof
        if not os.path.exists(fn):
            fn = fn[:-3]
        if not os.path.exists(fn):
            continue
        with (gzip.open(fn, "rt", newline="") if fn.endswith(".gz") else open(fn, newline="")) as infile:
            reader = csv.reader(infile, delimiter="\t")
            next(reader, None)  # Orthogroup, Species, <this species>, Orthologs
            for og, other, genes, orthologs in reader:
                j = index[other]
                if i < j:  # each relationship is in both species' files: keep it once
                    pairs.setdefault((i, j), []).append([og, split_genes(genes), split_genes(orthologs)])
    if not pairs:
        return []
    groups, views = [], []
    for (i, j), rows in sorted(pairs.items()):
        for og, genes_i, genes_j in rows:
            genes = sorted([gene_ids[names[i]][g] for g in genes_i] + [gene_ids[names[j]][g] for g in genes_j])
            groups.append((str(len(groups) + 1), genes, (("orthogroup", og),)))
        views.append({"name": f"{names[i]} - {names[j]}",
                      "columns": ["@Orthogroup", "Type", "*" + names[i], "*" + names[j]],
                      "rows": [[og, relationship(a, b), a, b] for og, a, b in rows]})
    out_xml = os.path.join(orthologs_dir, "Orthologs.orthoxml")
    used = sorted({names[k] for pair in pairs for k in pair}, key=index.get)
    write_orthoxml(out_xml, genes_of(species, gene_ids, groups, used), groups,
                   notes="Orthologs: each group is one ortholog relationship of the Orthologs files, "
                         "between the genes of two species (genes of the same species are in-paralogs).")
    out_html = os.path.join(orthologs_dir, "Orthologs.html")
    all_rows = [r for v in views for r in v["rows"]]
    write_page(
        out_html, title="Orthologs",
        note="The orthologs of the Orthologs files, one species pair at a time (also in Orthologs.orthoxml).",
        stats=[(len(views), "species pairs"), (len(all_rows), "ortholog relationships"),
               (sum(len(r[2]) * len(r[3]) for r in all_rows), "ortholog gene pairs"),
               (sum(r[1] == "one-to-one" for r in all_rows), "one-to-one")],
        views=views)
    return [out_xml, out_html]


def relationship(genes_a, genes_b):
    return f"{'one' if len(genes_a) == 1 else 'many'}-to-{'one' if len(genes_b) == 1 else 'many'}"


def genes_of(species, gene_ids, groups, species_names):
    """(name, FASTA file, [(gene id, gene name)]) of the given species, for the genes in the groups."""
    used = {g for _, genes, _ in groups for g in genes}
    fasta = dict(species)
    out = []
    for name in species_names:
        genes = sorted((i, g) for g, i in gene_ids[name].items() if i in used)
        if genes:
            out.append((name, fasta[name], genes))
    return out


# --- writing -----------------------------------------------------------------

def percent(a, b):
    return f"{100 * a / b:.1f}%" if b else "-"


def write_orthoxml(filename, species_genes, groups, notes=None):
    """
    species_genes: [(species name, FASTA file, [(gene id, gene name)])]
    groups: [(group id, [gene ids], [(property name, value)])]
    """
    with open(filename, "w") as out:
        out.write('<?xml version="1.0" encoding="utf-8"?>\n')
        out.write(f'<orthoXML xmlns={quoteattr(NAMESPACE)} xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation={quoteattr(SCHEMA_LOCATION)} origin="OrthoFinder" '
                  f'originVersion={quoteattr(__version__)} version="0.3">\n')
        if notes:
            out.write(f"  <notes>{escape(notes)}</notes>\n")
        for name, fasta_filename, genes in species_genes:
            out.write(f'  <species name={quoteattr(name)} NCBITaxId="{UNKNOWN_TAX_ID}">\n')
            out.write(f'    <database name={quoteattr(fasta_filename)} version="{UNKNOWN_VERSION}">\n')
            out.write("      <genes>\n")
            out.writelines(f'        <gene id="{i}" geneId={quoteattr(g)}/>\n' for i, g in genes)
            out.write("      </genes>\n")
            out.write("    </database>\n")
            out.write("  </species>\n")
        out.write("  <groups>\n")
        for group_id, genes, properties in groups:
            out.write(f"    <orthologGroup id={quoteattr(group_id)}>\n")
            out.writelines(f"      <property name={quoteattr(p)} value={quoteattr(v)}/>\n" for p, v in properties)
            out.writelines(f'      <geneRef id="{i}"/>\n' for i in genes)
            out.write("    </orthologGroup>\n")
        out.write("  </groups>\n")
        out.write("</orthoXML>\n")


def write_page(filename, title, note, stats, views, info=None):
    """
    A standalone page (the data is in the page) for a web browser.

    stats: [(value, label)] shown at the top
    info: an optional table {"title", "columns", "rows"} shown before the views
    views: [{"name", "columns", "rows"}], one shown at a time (a selector if
           there are several). A column label's first character gives its kind:
           "@" an id, "#" a number, "*" a list of genes; otherwise text.
    """
    data = {"title": title, "note": note, "version": __version__,
            "stats": [[v if isinstance(v, str) else f"{v:,}", label] for v, label in stats],
            "info": info, "views": views}
    with open(filename, "w") as out:
        out.write(HTML.replace("__TITLE__", escape(title))
                  .replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))


HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root { --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e8; --head:#f5f5f7; --accent:#2f6fde; --chip:#eef3fd; }
@media (prefers-color-scheme: dark) { :root { --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --head:#1f1f23; --accent:#7aa7ff; --chip:#1e2a40; } }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width:1200px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:20px; margin:0 0 4px; } h2 { font-size:15px; margin:28px 0 8px; }
.sub { color:var(--muted); font-size:12px; }
.stats { display:flex; gap:12px; flex-wrap:wrap; margin:16px 0; }
.stat { border:1px solid var(--line); border-radius:8px; padding:10px 14px; min-width:120px; }
.stat b { display:block; font-size:20px; font-variant-numeric:tabular-nums; } .stat span { color:var(--muted); font-size:12px; }
.wrap { overflow-x:auto; border:1px solid var(--line); border-radius:8px; }
table { border-collapse:collapse; width:100%; }
th, td { text-align:left; padding:7px 10px; border-bottom:1px solid var(--line); vertical-align:top; }
th { background:var(--head); font-weight:600; white-space:nowrap; }
td.num, th.num { font-variant-numeric:tabular-nums; text-align:right; }
tr:last-child td { border-bottom:none; }
.id { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; color:var(--accent); white-space:nowrap; }
td.text { white-space:nowrap; }
.gene { display:inline-block; background:var(--chip); border-radius:4px; padding:1px 6px; margin:1px 3px 1px 0;
        font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:12px; }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:8px; }
input, select { font:inherit; padding:7px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
input { width:min(380px,100%); } select { max-width:100%; }
button { font:inherit; padding:6px 12px; border:1px solid var(--line); border-radius:6px; background:var(--head); color:var(--fg); cursor:pointer; }
button:disabled { opacity:.4; cursor:default; } .muted { color:var(--muted); font-size:12px; }
details { display:inline; } summary { display:inline; cursor:pointer; color:var(--accent); font-size:12px; white-space:nowrap; }
</style></head><body><main>
<h1 id="title"></h1>
<div class="sub" id="sub"></div>
<div class="stats" id="stats"></div>
<div id="info"></div>
<h2 id="views_title"></h2>
<div class="tools">
  <select id="view"></select>
  <input id="q" placeholder="Search an ID or gene name" autocomplete="off">
  <button id="prev">Previous</button><button id="next">Next</button><span class="muted" id="count"></span>
</div>
<div class="wrap"><table id="table"></table></div>
</main>
<script>
const D = __DATA__;
const PAGE = 100, SHOW = 6;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const kind = c => ({"@": "id", "#": "num", "*": "genes"})[c[0]] || "text";
const label = c => kind(c) === "text" ? c : c.slice(1);
const chips = gs => gs.map(g => `<span class="gene">${esc(g)}</span>`).join("");
const genes = gs => gs.length <= SHOW ? chips(gs)
  : chips(gs.slice(0, SHOW)) + `<details><summary>+${gs.length - SHOW} more</summary>${chips(gs.slice(SHOW))}</details>`;
const th = c => `<th${kind(c) === "num" ? " class='num'" : ""}>${esc(label(c))}</th>`;
const td = (c, v) => kind(c) === "genes" ? `<td>${genes(v)}</td>`
  : kind(c) === "num" ? `<td class='num'>${typeof v === "number" ? v.toLocaleString() : esc(v)}</td>`
  : kind(c) === "id" ? `<td class='id'>${esc(v)}</td>` : `<td class='text'>${esc(v)}</td>`;
$("title").textContent = D.title;
$("sub").textContent = "OrthoFinder " + D.version + ". " + D.note;
$("stats").innerHTML = D.stats.map(([v, l]) => `<div class="stat"><b>${esc(v)}</b><span>${esc(l)}</span></div>`).join("");
if (D.info) $("info").innerHTML = `<h2>${esc(D.info.title)}</h2><div class="wrap"><table>` +
  "<tr>" + D.info.columns.map(th).join("") + "</tr>" +
  D.info.rows.map(r => "<tr>" + r.map((v, i) => td(D.info.columns[i], v)).join("") + "</tr>").join("") + "</table></div>";
$("views_title").textContent = D.views.length > 1 ? "Show" : D.views[0].name;
$("view").innerHTML = D.views.map((v, i) => `<option value="${i}">${esc(v.name)} (${v.rows.length.toLocaleString()})</option>`).join("");
$("view").style.display = D.views.length > 1 ? "" : "none";
const text = r => r.map(v => Array.isArray(v) ? v.join(" ") : String(v)).join(" ").toLowerCase();
let view = D.views[0], rows = view.rows, page = 0;
function draw() {
  const start = page * PAGE;
  $("table").innerHTML = "<tr>" + view.columns.map(th).join("") + "</tr>" +
    rows.slice(start, start + PAGE).map(r => "<tr>" + r.map((v, i) => td(view.columns[i], v)).join("") + "</tr>").join("");
  $("count").textContent = rows.length
    ? `${(start + 1).toLocaleString()}-${Math.min(start + PAGE, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}` : "no matches";
  $("prev").disabled = page === 0;
  $("next").disabled = start + PAGE >= rows.length;
}
function filter() {
  const q = $("q").value.trim().toLowerCase();
  rows = q ? view.rows.filter(r => text(r).includes(q)) : view.rows;
  page = 0; draw();
}
$("view").addEventListener("change", e => { view = D.views[+e.target.value]; filter(); });
$("q").addEventListener("input", filter);
$("prev").onclick = () => { page--; draw(); };
$("next").onclick = () => { page++; draw(); };
draw();
</script></body></html>
"""
