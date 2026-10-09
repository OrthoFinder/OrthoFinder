"""
--tree-view (a hidden option): pictures of the species tree and of the
resolved gene trees, drawn with ete4 (whose drawing needs PyQt6), and a page
to look at them in a web browser:

    Species_Tree/SpeciesTree_rooted_node_labels.svg
    Resolved_Gene_Trees/Tree_Images/<orthogroup>.svg
    Resolved_Gene_Trees/Resolved_Gene_Trees.html

Genes are coloured by species, the same colour in all the pictures. The
internal nodes are labelled as in the results: N0, N1, ... in the species tree
(the HOG files), n0, n1, ... in the gene trees (their "Gene tree parent clade").
"""

import contextlib
import csv
import io
import json
import os

from .. import __version__

# colour-blind friendly, distinguishable on white (Okabe-Ito, then Tableau 10)
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255", "#117733",
           "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#B07AA1", "#9C755F", "#BAB0AC"]
NODE_LABEL_COLOUR = "#777777"
IMAGES_DIR = "Tree_Images"


class TreeViewUnavailable(Exception):
    pass


def missing_pyqt6_message(results_dir):
    """Why --tree-view did not draw the trees, and how to draw them for these results."""
    return (
        "--tree-view: the trees were not drawn. OrthoFinder draws them with ete4, whose drawing "
        "uses the Qt graphics library (PyQt6). PyQt6 is not installed with OrthoFinder, as it is "
        "large and only this option needs it. To draw the trees, install it:\n"
        "    pip install PyQt6      (or: pip install \"orthofinder[trees]\")\n"
        "and then draw the trees of these results, without running the analysis again:\n"
        f"    python -m orthofinder.phylo.tree_pages {results_dir}"
    )


def load_treeview():
    """ete4's drawing classes; raises TreeViewUnavailable without PyQt6."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")  # draw to files, without a display
    try:
        from ete4 import Tree
        from ete4.parser.newick import NewickError
        from ete4.treeview import TextFace, TreeStyle, faces
    except ImportError as e:
        raise TreeViewUnavailable(str(e))
    return Tree, NewickError, TreeStyle, TextFace, faces


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


def write_tree_pages(results_dir):
    """
    Draw the trees of a results directory and write the page. Returns the
    number of gene trees drawn (0 if the directory has no resolved gene trees).
    Raises TreeViewUnavailable without PyQt6.
    """
    gene_trees_fn = os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    species_tree_fn = os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    if not os.path.exists(gene_trees_fn):
        return 0
    with open(gene_trees_fn) as infile:
        entries = [line.strip().split(": ", 1) for line in infile if line.strip()]
    if not entries:  # e.g. -og: the file is there but has no trees
        return 0
    Tree, NewickError, TreeStyle, TextFace, faces = load_treeview()
    species, gene_species = species_and_genes(os.path.join(results_dir, "Orthogroups"))
    colour = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(species)}
    species_by_length = sorted(species, key=len, reverse=True)

    def draw(tree, out_fn, leaf_label):
        def layout(node):
            if node.is_leaf:
                text, col = leaf_label(node.name)
                faces.add_face_to_node(TextFace(" " + text, fgcolor=col, fsize=9), node, column=0)
            elif node.name:
                faces.add_face_to_node(TextFace(node.name, fgcolor=NODE_LABEL_COLOUR, fsize=7), node,
                                       column=0, position="branch-top")
        style = TreeStyle()
        style.show_leaf_name = False
        style.show_scale = True
        style.layout_fn = layout
        with contextlib.redirect_stdout(io.StringIO()):  # ete4 prints "Wrote file: ..."
            tree.render(out_fn, tree_style=style)

    images_dir = os.path.join(results_dir, "Resolved_Gene_Trees", IMAGES_DIR)
    os.makedirs(images_dir, exist_ok=True)
    trees = []  # [orthogroup, [gene names]]
    n_failed = 0
    for og, newick in entries:
        try:
            tree = Tree(newick, parser=1)
            genes = []

            def leaf_label(name, genes=genes):
                sp, gene = leaf_species(name, species_by_length, gene_species)
                genes.append(gene)
                return gene, colour.get(sp, "#000000")

            draw(tree, os.path.join(images_dir, og + ".svg"), leaf_label)
            trees.append([og, sorted(set(genes))])
        except (NewickError, ValueError, OSError):  # a tree that cannot be read or drawn
            n_failed += 1
    if n_failed:
        print(f"WARNING: {n_failed} gene trees could not be drawn")
    if not trees:
        os.rmdir(images_dir)  # empty: nothing was drawn
        return 0

    species_tree_img = None
    if os.path.exists(species_tree_fn):
        with open(species_tree_fn) as infile:
            species_tree = Tree(infile.read().strip(), parser=1)
        species_tree_img = os.path.join(os.path.dirname(species_tree_fn), "SpeciesTree_rooted_node_labels.svg")
        draw(species_tree, species_tree_img, lambda name: (name, colour.get(name, "#000000")))
    write_page(os.path.join(results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.html"),
               species, colour, trees,
               os.path.relpath(species_tree_img, os.path.join(results_dir, "Resolved_Gene_Trees"))
               if species_tree_img else None)
    return len(trees)


def write_page(filename, species, colour, trees, species_tree_img):
    data = {"version": __version__, "images": IMAGES_DIR, "species_tree": species_tree_img,
            "species": [[sp, colour[sp]] for sp in species], "trees": trees}
    with open(filename, "w") as out:
        out.write(HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))


HTML = """<!doctype html>
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
a { color:var(--accent); }
.tree { border:1px solid var(--line); border-radius:8px; padding:12px; background:#fff; overflow:auto; }
.tree img { display:block; max-width:none; }
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
  <span class="muted" id="count"></span> <a id="open" target="_blank">open the picture</a>
</div>
<div class="tree"><img id="img" alt=""></div>
</main>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
$("sub").textContent = "OrthoFinder " + D.version + ". The trees of Resolved_Gene_Trees.txt, genes coloured by species; "
  + "n0, n1, ... are the gene tree nodes (the HOG files' gene tree parent clade).";
$("legend").innerHTML = D.species.map(([s, c]) => `<span style="--c:${c}">${esc(s)}</span>`).join("");
if (D.species_tree) $("species").innerHTML = `<details open><summary><b>Species tree</b> <span class="muted">(N0, N1, ... as in the HOG files)</span></summary>`
  + `<div class="tree"><img src="${esc(D.species_tree)}" alt="Species tree"></div></details>`;
let shown = D.trees, i = 0;
function show() {
  $("og").innerHTML = shown.map((t, k) => `<option value="${k}">${esc(t[0])} (${t[1].length} genes)</option>`).join("");
  $("og").value = String(i);
  const t = shown[i];
  $("img").src = t ? `${D.images}/${t[0]}.svg` : "";
  $("img").alt = t ? `Gene tree of ${t[0]}` : "no gene tree";
  $("open").href = t ? `${D.images}/${t[0]}.svg` : "#";
  $("count").textContent = shown.length ? `${i + 1} of ${shown.length.toLocaleString()}` : "no matches";
  $("prev").disabled = i <= 0;
  $("next").disabled = i >= shown.length - 1;
}
$("q").addEventListener("input", e => {
  const q = e.target.value.trim().toLowerCase();
  shown = q ? D.trees.filter(([og, genes]) => og.toLowerCase().includes(q) || genes.some(g => g.toLowerCase().includes(q))) : D.trees;
  i = 0; show();
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
    """python -m orthofinder.phylo.tree_pages <results directory>: draw its trees and write the page."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m orthofinder.phylo.tree_pages",
        description="Draw the species tree and the resolved gene trees of OrthoFinder results "
                    "(ete4, which needs PyQt6) and write Resolved_Gene_Trees/Resolved_Gene_Trees.html.")
    parser.add_argument("results_dir", help="an OrthoFinder results directory (Results_<name>)")
    results_dir = parser.parse_args(args).results_dir
    if not os.path.isdir(results_dir):
        parser.error(f"not a directory: {results_dir}")
    try:
        n_trees = write_tree_pages(results_dir)
    except TreeViewUnavailable:
        print(missing_pyqt6_message(results_dir))
        return 1
    if not n_trees:
        print(f"No resolved gene trees in {results_dir}")
        return 1
    print(f"Gene tree pictures ({n_trees}): "
          f"{os.path.join(results_dir, 'Resolved_Gene_Trees', 'Resolved_Gene_Trees.html')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
