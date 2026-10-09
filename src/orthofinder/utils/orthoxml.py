"""
With -x: the orthogroups in OrthoXML format (http://orthoxml.org, version 0.3),
Orthogroups.orthoxml, and as a page to look at in a web browser,
Orthogroups.html, both next to Orthogroups.tsv.

The species information comes from the input FASTA files. OrthoXML requires
an NCBI taxonomy ID and a database for each species, which OrthoFinder does
not know: the taxonomy ID is 0 (unknown) and the database is the FASTA file.
"""

import json
from xml.sax.saxutils import quoteattr

from .. import __version__

NAMESPACE = "http://orthoXML.org/2011/"
SCHEMA_LOCATION = "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd"
UNKNOWN_TAX_ID = "0"
UNKNOWN_VERSION = "NA"


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


def write_results(species, ogs, n_seqs_per_species, id_dict, species_to_use, results_base_filename):
    """
    Write <results_base_filename>.orthoxml and <results_base_filename>.html.

    species: from species_for_analysis, one entry per species in species_to_use
    ogs: the orthogroups, sets of sequence IDs ("<species>_<sequence>"); orthogroup
         i is OG<i> as in Orthogroups.tsv. Single-gene orthogroups are not groups.
    n_seqs_per_species: {species ID: number of sequences}
    id_dict: {sequence ID: gene name}
    """
    write_orthoxml(species, ogs, n_seqs_per_species, id_dict, species_to_use, results_base_filename + ".orthoxml")
    write_html(species, ogs, n_seqs_per_species, id_dict, species_to_use, results_base_filename + ".html")


def write_orthoxml(species, ogs, n_seqs_per_species, id_dict, species_to_use, filename):
    gene_number = {}  # sequence ID -> the gene's id in the file (unique across species)
    with open(filename, "w") as out:
        out.write('<?xml version="1.0" encoding="utf-8"?>\n')
        out.write(f'<orthoXML xmlns={quoteattr(NAMESPACE)} xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation={quoteattr(SCHEMA_LOCATION)} origin="OrthoFinder" '
                  f'originVersion={quoteattr(__version__)} version="0.3">\n')
        for i_species, (name, fasta_filename) in zip(species_to_use, species):
            out.write(f'  <species name={quoteattr(name)} NCBITaxId="{UNKNOWN_TAX_ID}">\n')
            out.write(f'    <database name={quoteattr(fasta_filename)} version="{UNKNOWN_VERSION}">\n')
            out.write("      <genes>\n")
            for i_seq in range(n_seqs_per_species[i_species]):
                seq_id = f"{i_species}_{i_seq}"
                gene_number[seq_id] = len(gene_number)
                out.write(f'        <gene id="{gene_number[seq_id]}" geneId={quoteattr(id_dict[seq_id])}/>\n')
            out.write("      </genes>\n")
            out.write("    </database>\n")
            out.write("  </species>\n")
        out.write("  <groups>\n")
        for i_og, og in enumerate(ogs):
            if len(og) < 2:
                continue
            out.write(f'    <orthologGroup id="OG{i_og:07d}">\n')
            out.writelines(f'      <geneRef id="{i_gene}"/>\n' for i_gene in sorted(gene_number[seq_id] for seq_id in og))
            out.write("    </orthologGroup>\n")
        out.write("  </groups>\n")
        out.write("</orthoXML>\n")


def write_html(species, ogs, n_seqs_per_species, id_dict, species_to_use, filename):
    """A standalone page (the data is in the page) showing the orthogroups."""
    column = {i_species: i for i, i_species in enumerate(species_to_use)}
    groups = []
    n_assigned = [0] * len(species_to_use)
    for i_og, og in enumerate(ogs):
        if len(og) < 2:
            continue
        per_species = [[] for _ in species_to_use]
        for seq_id in sorted(og, key=lambda s: tuple(map(int, s.split("_")))):
            i_col = column[int(seq_id.split("_")[0])]
            per_species[i_col].append(id_dict[seq_id])
            n_assigned[i_col] += 1
        groups.append([f"OG{i_og:07d}", per_species])
    data = {
        "version": __version__,
        "species": [{"name": name, "fasta": fasta, "n_genes": n_seqs_per_species[i], "n_assigned": n}
                    for (name, fasta), i, n in zip(species, species_to_use, n_assigned)],
        "groups": groups,
    }
    with open(filename, "w") as out:
        out.write(HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))


HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>OrthoFinder Orthogroups</title>
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
.og { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; color:var(--accent); white-space:nowrap; }
.gene { display:inline-block; background:var(--chip); border-radius:4px; padding:1px 6px; margin:1px 3px 1px 0;
        font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:12px; }
.tools { display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:8px; }
input { font:inherit; padding:7px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); width:min(420px,100%); }
button { font:inherit; padding:6px 12px; border:1px solid var(--line); border-radius:6px; background:var(--head); color:var(--fg); cursor:pointer; }
button:disabled { opacity:.4; cursor:default; } .muted { color:var(--muted); font-size:12px; }
details { display:inline; } summary { display:inline; cursor:pointer; color:var(--accent); font-size:12px; white-space:nowrap; }
</style></head><body><main>
<h1>Orthogroups</h1>
<div class="sub" id="sub"></div>
<div class="stats" id="stats"></div>
<h2>Species</h2>
<div class="wrap"><table id="species"></table></div>
<h2>Orthogroups</h2>
<div class="tools">
  <input id="q" placeholder="Search an orthogroup ID or gene name" autocomplete="off">
  <button id="prev">Previous</button><button id="next">Next</button><span class="muted" id="info"></span>
</div>
<div class="wrap"><table id="groups"></table></div>
</main>
<script>
const D = __DATA__;
const PAGE = 100, SHOW = 6;
const chips = gs => gs.map(g => `<span class="gene">${esc(g)}</span>`).join("");
const cell = gs => gs.length <= SHOW ? chips(gs)
  : chips(gs.slice(0, SHOW)) + `<details><summary>+${gs.length - SHOW} more</summary>${chips(gs.slice(SHOW))}</details>`;
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const pct = (a, b) => b ? (100 * a / b).toFixed(1) + "%" : "-";
document.getElementById("sub").textContent = "OrthoFinder " + D.version + ". The same orthogroups as Orthogroups.tsv and Orthogroups.orthoxml.";
const nGenes = D.species.reduce((a, s) => a + s.n_genes, 0);
const nAssigned = D.species.reduce((a, s) => a + s.n_assigned, 0);
document.getElementById("stats").innerHTML = [
  [D.species.length, "species"], [D.groups.length.toLocaleString(), "orthogroups"],
  [nGenes.toLocaleString(), "genes"], [pct(nAssigned, nGenes), "genes in orthogroups"]
].map(([v, l]) => `<div class="stat"><b>${v}</b><span>${l}</span></div>`).join("");
document.getElementById("species").innerHTML =
  "<tr><th>Species</th><th>FASTA file</th><th class='num'>Genes</th><th class='num'>In orthogroups</th></tr>" +
  D.species.map(s => `<tr><td>${esc(s.name)}</td><td>${esc(s.fasta)}</td>` +
    `<td class="num">${s.n_genes.toLocaleString()}</td><td class="num">${pct(s.n_assigned, s.n_genes)}</td></tr>`).join("");
const head = "<tr><th>Orthogroup</th><th class='num'>Genes</th>" + D.species.map(s => `<th>${esc(s.name)}</th>`).join("") + "</tr>";
let rows = D.groups, page = 0;
function draw() {
  const start = page * PAGE, shown = rows.slice(start, start + PAGE);
  document.getElementById("groups").innerHTML = head + shown.map(([id, per]) =>
    `<tr><td class="og">${esc(id)}</td><td class="num">${per.reduce((a, g) => a + g.length, 0)}</td>` +
    per.map(gs => `<td>${cell(gs)}</td>`).join("") + "</tr>").join("");
  document.getElementById("info").textContent = rows.length
    ? `${(start + 1).toLocaleString()}-${Math.min(start + PAGE, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}` : "no matches";
  document.getElementById("prev").disabled = page === 0;
  document.getElementById("next").disabled = start + PAGE >= rows.length;
}
document.getElementById("q").addEventListener("input", e => {
  const q = e.target.value.trim().toLowerCase();
  rows = q ? D.groups.filter(([id, per]) => id.toLowerCase().includes(q) || per.some(gs => gs.some(g => g.toLowerCase().includes(q)))) : D.groups;
  page = 0; draw();
});
document.getElementById("prev").onclick = () => { page--; draw(); };
document.getElementById("next").onclick = () => { page++; draw(); };
draw();
</script></body></html>
"""
