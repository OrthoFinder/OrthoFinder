"""
The pages to view a run's results in a web browser, and the links between
them: --tree-view, --stats-view and --groups-view (or --view for all three).
Each page is one directory down from the results directory, so a link from
one page to another is "../<its path>".

With two or more pages, a run also writes OrthoFinder_Results.html in the
results directory: a page with a tab for each, shown in it (an iframe), so
that there is one page to open. A page shown in it hides its own bar of links
to the others, and its links to them open the other page's tab (messages to
the parent page: a page opened from disk cannot reach into another).
"""

import heapq
import json
import os
import pickle
import re
import sys
import tempfile

from .. import __version__

PAGES = {  # key: (name, path in the results directory)
    "trees": ("Gene trees", "Resolved_Gene_Trees/Resolved_Gene_Trees.html"),
    "stats": ("Statistics", "Comparative_Genomics_Statistics/Statistics.html"),
    "groups": ("Orthogroups & HOGs", "Orthogroups/Groups.html"),
}


def orthogroups_file(results_dir, suffix=".tsv"):
    """
    Orthogroups/Orthogroups<suffix>, or without it (e.g. after -og)
    MCL_Orthogroups<suffix>, the MCL orthogroups; None if neither is there.
    """
    for name in ("Orthogroups", "MCL_Orthogroups"):
        path = os.path.join(results_dir, "Orthogroups", name + suffix)
        if os.path.exists(path):
            return path
    return None


def present(results_dir, written=None):
    """
    The pages of a results directory: written, the keys of the pages a run is
    writing (they may not be there yet), or if None, the pages that are there.
    """
    if written is not None:
        return set(written)
    return {key for key, (_, path) in PAGES.items() if os.path.exists(os.path.join(results_dir, path))}


def link(results_dir, key, written=None):
    """The link to a page from another page, or None if it is not there."""
    return "../" + PAGES[key][1] if key in present(results_dir, written) else None


def nav(results_dir, current, written=None):
    """[[name, link, 1 if it is this page, key]] of the pages, for the bar at the top of each."""
    keys = present(results_dir, written) | {current}
    return [[name, "../" + path, int(key == current), key] for key, (name, path) in PAGES.items() if key in keys]


# ---- the pages' data ----
# A page holds a summary; the rest is in data files, each a script the page loads when it is needed (a page
# opened from disk may load scripts, but not read other files): <data dir>/<name>.js, which calls ofData(name,
# value) (DATA_JS). Large tables are split into chunks, so that a page loads the rows it shows, and genes are
# found with an index (write_gene_index), so that thousands of species need not be loaded to find one gene.

def write_data(data_dir, name, value):
    """<data_dir>/<name>.js: a data file, which gives value to the page (loadData(name), DATA_JS)."""
    with open(os.path.join(data_dir, name + ".js"), "w", encoding="utf-8", newline="\n") as out:
        out.write(f"ofData({json.dumps(name)},{json.dumps(value, separators=(',', ':'))});\n")


NON_EMPTY = re.compile(r"[^\t\r\n]+")


def sparse_cells(line):
    """
    (column, text) of the non-empty cells of a line of a tab-separated file
    (column 0 the first), found without splitting the line at every tab: a
    row of thousands of species (Orthogroups.tsv, the HOG files) is mostly
    empty cells.
    """
    col, last = 0, 0
    for m in NON_EMPTY.finditer(line):
        col += line.count("\t", last, m.start())
        last = m.start()
        yield col, m.group()


def compact_ids(ids):
    """
    IDs for a page: {"prefix", "width", "n"} if they are <prefix><0, 1, 2, ...
    zero-padded to width> in order (OG0000000, OG0000001, ...; N3.HOG0000000,
    ...), which the page makes again (expandIds, DATA_JS), else the list.
    """
    if ids:
        m = re.match(r"(.*?)(\d+)$", ids[0])
        if m and int(m.group(2)) == 0:
            prefix, width = m.group(1), len(m.group(2))
            if all(x == f"{prefix}{i:0{width}d}" for i, x in enumerate(ids)):
                return {"prefix": prefix, "width": width, "n": len(ids)}
    return list(ids)


def write_chunks(data_dir, prefix, rows, weight, limit=50000):
    """
    rows in data files <prefix>_<c>, each of rows of about limit in weight
    (weight(row): e.g. its genes), so that the page loads the rows it shows
    (chunkRows, DATA_JS). Returns the first row of each chunk. rows can be
    an iterator: each chunk is written once it is complete, so only one is
    held at a time.
    """
    starts, total, chunk = [0], 0, []
    i = -1
    for i, row in enumerate(rows):
        w = weight(row) + 1
        if total and total + w > limit:
            write_data(data_dir, f"{prefix}_{len(starts) - 1}", chunk)
            starts.append(i)
            total, chunk = 0, []
        total += w
        chunk.append(row)
    write_data(data_dir, f"{prefix}_{len(starts) - 1}", chunk)
    return starts


# The genes of the gene index sorted in memory at a time (write_gene_index): more are
# sorted in runs of this many, kept in temporary files, and merged
GENE_INDEX_RUN = 1000000


def _gene_key(entry):
    return entry[0].lower()


# Temporary files of entries (genes): pickled in blocks of this many, read back block by block
_BLOCK = 10000


class EntryFile(object):
    """
    A temporary file of entries, written then read back in order (in blocks of
    _BLOCK, pickled: much quicker than a line of text each). Use as a context
    manager; add(entry) then entries().
    """

    def __init__(self, directory):
        self.file = tempfile.TemporaryFile(dir=directory)
        self.block = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.file.close()

    def add(self, entry):
        self.block.append(entry)
        if len(self.block) >= _BLOCK:
            self._flush()

    def _flush(self):
        if self.block:
            pickle.dump(self.block, self.file, pickle.HIGHEST_PROTOCOL)
            self.block = []

    def entries(self):
        self._flush()
        self.file.seek(0)
        while True:
            try:
                block = pickle.load(self.file)
            except EOFError:
                return
            yield from block


def _sorted_genes(genes, tmp_dir, run=GENE_INDEX_RUN):
    """
    The (gene, value) entries sorted by gene ignoring case, as sorted() would
    sort them (equal names in the order given), in bounded memory: in memory if
    there are at most run of them, else sorted runs written to tmp_dir and
    merged (both sorts are stable, and heapq.merge takes equal names from the
    earlier run first).
    """
    genes = iter(genes)
    first = []
    for entry in genes:
        first.append(entry)
        if len(first) >= run:
            break
    else:
        yield from sorted(first, key=_gene_key)
        return
    runs = []
    batch = first
    try:
        while batch:
            batch.sort(key=_gene_key)
            run_file = EntryFile(tmp_dir)
            for entry in batch:
                run_file.add(entry)
            runs.append(run_file)
            batch = []
            for entry in genes:
                batch.append(entry)
                if len(batch) >= run:
                    break
        yield from heapq.merge(*(r.entries() for r in runs), key=_gene_key)
    finally:
        for r in runs:
            r.file.close()


def write_gene_index(data_dir, genes, shard=20000, run=GENE_INDEX_RUN):
    """
    An index of gene names, to find genes by the start of their name
    (findGenes, DATA_JS): genes [(gene, value)] (or an iterator of them), e.g.
    the orthogroup each is in, sorted (ignoring case) into data files Genes_<i>
    of shard genes each, in bounded memory (_sorted_genes: tens of millions of
    genes are not all held at once).
    Returns the first gene (in lower case) of each, for the page.
    """
    index = []
    with tempfile.TemporaryDirectory(prefix="gene_index_", dir=data_dir) as tmp_dir:
        batch = []
        for entry in _sorted_genes(genes, tmp_dir, run):
            batch.append(list(entry))
            if len(batch) == shard:
                write_data(data_dir, f"Genes_{len(index)}", batch)
                index.append(batch[0][0].lower())
                batch = []
        if batch:
            write_data(data_dir, f"Genes_{len(index)}", batch)
            index.append(batch[0][0].lower())
    return index


# The script for the data files: in a page, after D (its data, with D.data its data directory) is set.
DATA_JS = r"""// ---- data files (view_pages.write_data): scripts, each loaded once, when first needed ----
const dataCache = {}, dataWaiting = {};
window.ofData = (name, value) => { dataCache[name] = value; (dataWaiting[name] || []).forEach(w => w.resolve(value)); delete dataWaiting[name]; };
function loadData(name) {
  if (name in dataCache) return Promise.resolve(dataCache[name]);
  return new Promise((resolve, reject) => {
    if (!dataWaiting[name]) {
      dataWaiting[name] = [];
      const s = document.createElement("script");
      s.src = `${D.data}/${encodeURIComponent(name)}.js`;
      s.onerror = () => {
        const e = new Error(`could not be loaded (${s.src}). The page needs its folder ${D.data} next to it: copy (or download) them together, e.g. the whole results directory`);
        (dataWaiting[name] || []).forEach(w => w.reject(e)); delete dataWaiting[name];
      };
      document.head.appendChild(s);
    }
    dataWaiting[name].push({resolve, reject});
  });
}
// IDs (compact_ids): the list, or {prefix, width, n} for <prefix>0000, <prefix>0001, ...
const expandIds = ids => Array.isArray(ids) ? ids : Array.from({length: ids.n}, (_, i) => ids.prefix + String(i).padStart(ids.width, "0"));
// rows in chunks (write_chunks): the chunk of a row, and the data of some rows (a Map: row -> data)
function chunkOf(starts, row) { let lo = 0, hi = starts.length - 1; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (starts[m] <= row) lo = m; else hi = m - 1; } return lo; }
async function chunkRows(prefix, starts, rows) {
  const chunks = [...new Set(rows.map(r => chunkOf(starts, r)))];
  await Promise.all(chunks.map(c => loadData(`${prefix}_${c}`)));
  return new Map(rows.map(r => { const c = chunkOf(starts, r); return [r, dataCache[`${prefix}_${c}`][r - starts[c]]]; }));
}
// genes (write_gene_index; index: the first gene of each file): [[gene, value], ...] of the genes whose
// name starts with q (in any case), at most limit
async function findGenes(index, q, limit = 1000) {
  q = q.toLowerCase();
  if (!index || !index.length || !q) return [];
  let s = 0;
  while (s + 1 < index.length && index[s + 1] < q) s++;  // (a gene's name can be in two files: in two species)
  const out = [];
  for (; s < index.length; s++) {
    for (const e of await loadData(`Genes_${s}`)) {
      const g = e[0].toLowerCase();
      if (g.startsWith(q)) { out.push(e); if (out.length >= limit) return out; }
      else if (g > q) return out;
    }
  }
  return out;
}"""


HUB = "OrthoFinder_Results.html"


def write_hub(results_dir, keys):
    """OrthoFinder_Results.html: a tab for each page (keys, of PAGES). Returns its filename."""
    pages = [[key, name, path] for key, (name, path) in PAGES.items() if key in keys]
    data = {"version": __version__, "pages": pages, "results": os.path.basename(os.path.normpath(results_dir))}
    fn = os.path.join(results_dir, HUB)
    with open(fn, "w", encoding="utf-8", newline="\n") as out:
        out.write(HUB_HTML.replace("__THEME_JS__", THEME_JS).replace("__THEME_HTML__", THEME_HTML)
                  .replace("__THEME_SELECT_JS__", THEME_SELECT_JS)
                  .replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
    return fn


# Light or dark colours: the viewer's choice, else the system's. Each page's style has its dark colours under
# @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) } and :root[data-theme="dark"].
# THEME_JS goes in the <head> (so that a page opens in its colours): it sets <html data-theme> from the choice
# kept in the browser, or from the parent page (OrthoFinder_Results.html), and then fires "orthofinder-theme".
THEME_JS = """<script>
const THEME_KEY = "orthofinder-theme";
function applyTheme(t) {
  if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
  else delete document.documentElement.dataset.theme;
  window.dispatchEvent(new Event("orthofinder-theme"));
}
try { applyTheme(localStorage.getItem(THEME_KEY)); } catch (e) { /* not available */ }
window.addEventListener("message", e => {
  if (e.source === window.parent && e.source !== window && e.data && e.data.orthofinder === "theme") applyTheme(e.data.theme);
});
</script>"""
THEME_HTML = ('<label class="theme" title="Light or dark colours (Auto: as the system\'s)">Colours '
              '<select id="theme-select"><option value="">Auto</option><option value="light">Light</option>'
              '<option value="dark">Dark</option></select></label>')
THEME_SELECT_JS = """const themeSelect = document.getElementById("theme-select");
try { themeSelect.value = localStorage.getItem(THEME_KEY) || ""; } catch (e) { /* not available */ }
themeSelect.onchange = () => {
  const t = themeSelect.value;
  try { if (t) localStorage.setItem(THEME_KEY, t); else localStorage.removeItem(THEME_KEY); } catch (e) { /* not available */ }
  applyTheme(t);
};"""

# The bar at the top of each page: its style, element and script (the page's data has "nav"), with the colours'
# choice at its right (in OrthoFinder_Results.html, the parent page has it instead).
NAV_CSS = (".pagebar { position:relative; } "
           ".pagenav { display:flex; gap:6px 18px; flex-wrap:wrap; align-items:center; font-size:13px; "
           "margin:0 0 14px; padding:0 140px 8px 0; border-bottom:1px solid var(--pagenav-line, var(--line, #e3e3e8)); } "
           ".pagenav[hidden] { display:none; } "
           ".pagenav b { font-weight:600; } .pagenav .what { color:var(--pagenav-muted, var(--muted, #6b6b70)); } "
           ".pagebar .theme { position:absolute; top:0; right:0; font-size:12px; color:var(--pagenav-muted, var(--muted, #6b6b70)); } "
           ".pagebar .theme[hidden] { display:none; } "
           ".theme select { font:inherit; font-size:12px; padding:2px 6px; width:auto; }")
NAV_HTML = '<div class="pagebar"><nav class="pagenav" id="pagenav" hidden></nav>' + THEME_HTML + '</div>'
NAV_JS = """const IN_HUB = window.parent !== window;  // shown in OrthoFinder_Results.html
if (IN_HUB) document.querySelector(".pagebar .theme").hidden = true;
else {
""" + "\n".join("  " + line for line in THEME_SELECT_JS.split("\n")) + """
}
if (D.nav && D.nav.length > 1 && !IN_HUB) {
  const el = document.getElementById("pagenav"), e = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
  el.hidden = false;
  el.innerHTML = `<span class="what">OrthoFinder results:</span>` + D.nav.map(([name, href, here]) =>
    here ? `<b>${e(name)}</b>` : `<a href="${e(href)}">${e(name)}</a>`).join("");
}
if (IN_HUB && D.nav) {
  // a link to another page opens its tab in the parent page, at the place linked to (#...)
  document.addEventListener("click", ev => {
    const a = ev.target.closest && ev.target.closest("a[href]");
    if (!a || ev.ctrlKey || ev.metaKey || ev.shiftKey) return;
    const [path, hash] = a.getAttribute("href").split("#");
    const page = D.nav.find(p => p[1] === path && !p[2]);
    if (!page) return;
    ev.preventDefault();
    window.parent.postMessage({orthofinder: "open", key: page[3], hash: hash || ""}, "*");
  });
}"""


# The parent page: a tab for each page, each shown when first opened and kept, and a link to open it on its own.
HUB_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>OrthoFinder results</title>
__THEME_JS__
<style>
:root { color-scheme:light; --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e8; --accent:#2f6fde; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --accent:#7aa7ff; } }
:root[data-theme="dark"] { color-scheme:dark; --bg:#161618; --fg:#ececf0; --muted:#9a9aa2; --line:#2e2e33; --accent:#7aa7ff; }
* { box-sizing:border-box; }
html, body { margin:0; height:100%; }
body { background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; display:flex; flex-direction:column; }
header { padding:10px clamp(16px, 2vw, 32px) 0; border-bottom:1px solid var(--line); }
.title { display:flex; gap:12px; align-items:baseline; flex-wrap:wrap; }
h1 { font-size:20px; margin:0; } .muted { color:var(--muted); font-size:12px; }
nav { display:flex; gap:4px; flex-wrap:wrap; margin-top:6px; }
nav button { font:inherit; border:0; background:none; color:var(--muted); padding:8px 14px; cursor:pointer; border-bottom:2px solid transparent; margin-bottom:-1px; }
nav button[aria-selected="true"] { color:var(--fg); border-bottom-color:var(--accent); font-weight:600; }
nav a { margin-left:auto; align-self:center; font-size:12px; color:var(--accent); }
.theme { margin-left:auto; color:var(--muted); font-size:12px; }
.theme select { font:inherit; padding:2px 6px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
#frames { flex:1; position:relative; }
iframe { position:absolute; inset:0; width:100%; height:100%; border:0; background:var(--bg); }
iframe[hidden] { display:none; }
</style></head><body>
<header>
  <div class="title"><h1>OrthoFinder results</h1><span class="muted" id="sub"></span>__THEME_HTML__</div>
  <nav id="tabs" role="tablist"></nav>
</header>
<div id="frames"></div>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
$("sub").textContent = `${D.results} \u00B7 OrthoFinder ${D.version}`;
$("tabs").innerHTML = D.pages.map(([key, name]) => `<button role="tab" data-key="${key}">${esc(name)}</button>`).join("")
  + `<a id="alone" target="_blank" title="This page on its own, in a new tab">open on its own</a>`;
const frames = {};
function open(key, hash) {
  const page = D.pages.find(p => p[0] === key) || D.pages[0];
  key = page[0];
  for (const b of document.querySelectorAll("#tabs button")) b.setAttribute("aria-selected", b.dataset.key === key);
  for (const [k, f] of Object.entries(frames)) f.hidden = k !== key;
  const url = page[2] + (hash ? "#" + hash : "");
  if (!frames[key]) {  // made when first opened, then kept (it keeps its place)
    const f = document.createElement("iframe");
    f.title = page[1]; f.src = url;
    f.addEventListener("load", () => sendTheme(f));
    $("frames").appendChild(f);
    frames[key] = f;
  } else if (hash) frames[key].src = url;  // the page goes to it (#...)
  $("alone").href = url;
  history.replaceState(null, "", "#" + key + (hash ? "/" + hash : ""));
}
$("tabs").onclick = e => { const b = e.target.closest("[data-key]"); if (b) open(b.dataset.key); };
window.addEventListener("message", e => {
  const m = e.data;
  if (m && m.orthofinder === "open") open(m.key, m.hash);
});
// the colours chosen here, for the pages shown in it
__THEME_SELECT_JS__
const sendTheme = f => f.contentWindow && f.contentWindow.postMessage({orthofinder: "theme", theme: document.documentElement.dataset.theme || ""}, "*");
window.addEventListener("orthofinder-theme", () => Object.values(frames).forEach(sendTheme));
// OrthoFinder_Results.html#stats, #trees/OG0000003, #groups/orthogroups/OG0000003: that page, at that place
const start = decodeURIComponent(location.hash.slice(1)), cut = start.indexOf("/");
open(cut < 0 ? start : start.slice(0, cut), cut < 0 ? "" : start.slice(cut + 1));
</script></body></html>
"""


# ---- serving the pages ----
# For results on another computer (e.g. a cluster): the pages served over HTTP, to open in a browser on this
# one, e.g. through an ssh tunnel. The pages are the same files; the server only sends them, compressed (the
# data files are text that compresses ~5-10x), and only the files of the results directory.

COMPRESS = {".html", ".js", ".json", ".css", ".svg", ".tsv", ".txt"}


def start_page(results_dir):
    """The page to open first: OrthoFinder_Results.html, else the one page there; None if there are none."""
    if os.path.exists(os.path.join(results_dir, HUB)):
        return HUB
    return next((path for _, path in PAGES.values() if os.path.exists(os.path.join(results_dir, path))), None)


def make_server(results_dir, host="127.0.0.1", port=8000, tries=20):
    """
    An HTTP server of the files of results_dir (not a list of them), at the first port free from port on
    (tries of them; port 0: any). Call its serve_forever().
    """
    import gzip
    import http.server
    from collections import OrderedDict

    root = os.path.abspath(results_dir)
    first = start_page(root)
    compressed = OrderedDict()  # (path, modified, size): gzip, the most recent last

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=root, **kwargs)

        def send_head(self):
            if self.path.split("?")[0] in ("", "/"):
                self.send_response(302)
                self.send_header("Location", "/" + (first or ""))
                self.end_headers()
                return None
            if self.path == "/favicon.ico":  # asked for by browsers: none
                self.send_response(204)
                self.end_headers()
                return None
            path = self.translate_path(self.path)
            if os.path.isdir(path):
                self.send_error(404, "Not a page")  # no lists of the files
                return None
            ext = os.path.splitext(path)[1].lower()
            if ext not in COMPRESS or "gzip" not in self.headers.get("Accept-Encoding", "") or not os.path.isfile(path):
                return super().send_head()
            st = os.stat(path)
            key = (path, st.st_mtime_ns, st.st_size)
            body = compressed.pop(key, None)
            if body is None:
                with open(path, "rb") as f:
                    body = gzip.compress(f.read(), compresslevel=5)
            compressed[key] = body
            while len(compressed) > 200:
                compressed.popitem(last=False)
            self.send_response(200)
            self.send_header("Content-Type", self.guess_type(path))
            self.send_header("Content-Encoding", "gzip")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-cache")  # pages written again are seen at once
            self.end_headers()
            if self.command == "GET":
                self.wfile.write(body)
            return None

        def log_request(self, code="-", size="-"):  # quiet: only the requests that failed
            if isinstance(code, int) and code >= 400:
                super().log_request(code, size)

    for p in [0] if port == 0 else range(port, port + tries):
        try:
            return http.server.ThreadingHTTPServer((host, p), Handler)
        except OSError:
            continue
    raise OSError(f"no free port from {port} to {port + tries - 1}")


def serve(results_dir, host="127.0.0.1", port=8000):
    """Serve the pages of results_dir until stopped (Ctrl-C), with how to open them."""
    import socket
    server = make_server(results_dir, host, port)
    port = server.server_address[1]
    page = start_page(results_dir) or ""
    print(f"Serving the pages of {results_dir} (Ctrl-C to stop):")
    if host in ("127.0.0.1", "localhost"):
        print(f"  on this computer:     http://localhost:{port}/{page}")
        print(f"  from another one:     ssh -N -L {port}:localhost:{port} {os.environ.get('USER', 'you')}@{socket.gethostname()}"
              f"\n                        then open http://localhost:{port}/{page}")
    else:
        print(f"  http://{socket.gethostname()}:{port}/{page}")
        print("  WARNING: anyone who can reach this computer on this port can read these results")
    sys.stdout.flush()  # before waiting (e.g. the output going to a file)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopped")
    finally:
        server.server_close()


def add_serve_options(parser):
    """--serve, --port and --host, for the commands that write the pages."""
    parser.add_argument("--serve", action="store_true",
                        help="then serve the pages over HTTP, to open them in a browser (e.g. results on a "
                             "cluster, through an ssh tunnel) until stopped with Ctrl-C")
    parser.add_argument("--port", type=int, default=8000, help="with --serve: the port (default 8000, or the next free one)")
    parser.add_argument("--host", default="127.0.0.1",
                        help="with --serve: the address to serve on (default 127.0.0.1, only this computer; "
                             "0.0.0.0: every computer that can reach it, with no password)")


def main(args=None):
    """python -m orthofinder.utils.view_pages <results directory>: OrthoFinder_Results.html for its pages."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.utils.view_pages",
        description="Write OrthoFinder_Results.html for OrthoFinder results: a tab for each of the pages "
                    "that are there (written with --tree-view, --stats-view, --groups-view or their commands).")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    add_serve_options(parser)
    a = parser.parse_args(args)
    keys = present(a.results_dir)
    if not keys:
        print(f"No pages in {a.results_dir} (see --tree-view, --stats-view, --groups-view)")
        return 1
    print(f"OrthoFinder results page ({', '.join(PAGES[k][0] for k in PAGES if k in keys)}): "
          f"{write_hub(a.results_dir, keys)}")
    if a.serve:
        serve(a.results_dir, a.host, a.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
