"""
The pages to view a run's results in a web browser (--tree-view, --stats-view,
--groups-view, or --view for all three), written from the results files: at
the end of a run, or later for existing results (tools/view_results).
"""

import os

from ..orthology import groups_page, stats_page
from ..phylo import tree_pages
from ..utils import view_pages


def write_views(results_dir, trees=False, stats=False, groups=False):
    """
    Write the pages asked for, linked to each other and to those already there
    (e.g. one written before with tools/view_results, which is not written
    again): each is told which of the others there are, also those written
    after it.

    Returns the pages written ([(what, path)]), OrthoFinder_Results.html (a tab
    for each page there) if there are two or more pages, else None, and all
    the pages there ([(what, path)]).
    """
    there = view_pages.present(results_dir)
    stats = stats and view_pages.orthogroups_file(results_dir, ".GeneCount.tsv") is not None
    groups = groups and view_pages.orthogroups_file(results_dir) is not None
    linked = there | ({"stats"} if stats else set()) | ({"groups"} if groups else set())
    written = {}  # key: (what, path)
    if trees:
        n_trees = tree_pages.write_tree_pages(results_dir, linked | {"trees"})
        if n_trees:
            linked.add("trees")
            written["trees"] = (f"Gene trees ({n_trees})", os.path.join(results_dir, view_pages.PAGES["trees"][1]))
        else:
            linked.discard("trees")
    if stats:
        page = stats_page.write_stats_page(results_dir, linked)
        if page:
            written["stats"] = ("Statistics", page)
    if groups:
        page = groups_page.write_groups_page(results_dir, written=linked)
        if page:
            written["groups"] = ("Orthogroups & HOGs", page)
    every = [written.get(key) or (name, os.path.join(results_dir, path))
             for key, (name, path) in view_pages.PAGES.items() if key in written or key in there]
    hub = view_pages.write_hub(results_dir, [k for k in view_pages.PAGES if k in written or k in there]) \
        if written and len(every) >= 2 else None
    return list(written.values()), hub, every
