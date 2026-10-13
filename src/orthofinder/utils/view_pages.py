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

import json
import os

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


HUB = "OrthoFinder_Results.html"


def write_hub(results_dir, keys):
    """OrthoFinder_Results.html: a tab for each page (keys, of PAGES). Returns its filename."""
    pages = [[key, name, path] for key, (name, path) in PAGES.items() if key in keys]
    data = {"version": __version__, "pages": pages, "results": os.path.basename(os.path.normpath(results_dir))}
    fn = os.path.join(results_dir, HUB)
    with open(fn, "w", encoding="utf-8") as out:
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


def main(args=None):
    """python -m orthofinder.utils.view_pages <results directory>: OrthoFinder_Results.html for its pages."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.utils.view_pages",
        description="Write OrthoFinder_Results.html for OrthoFinder results: a tab for each of the pages "
                    "that are there (written with --tree-view, --stats-view, --groups-view or their commands).")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    results_dir = parser.parse_args(args).results_dir
    keys = present(results_dir)
    if not keys:
        print(f"No pages in {results_dir} (see --tree-view, --stats-view, --groups-view)")
        return 1
    print(f"OrthoFinder results page ({', '.join(PAGES[k][0] for k in PAGES if k in keys)}): "
          f"{write_hub(results_dir, keys)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
