#!/usr/bin/env python3

"""
Write the pages to view OrthoFinder results in a web browser, from the results
files of a run that has finished (the same pages as a run with --view,
--tree-view, --stats-view or --groups-view):

    Resolved_Gene_Trees/Resolved_Gene_Trees.html       (--tree-view)
    Comparative_Genomics_Statistics/Statistics.html    (--stats-view)
    Orthogroups/Groups.html                            (--groups-view)
    OrthoFinder_Results.html                           (with two or more: a tab for each)

With no option, all three. INPUT is an OrthoFinder results directory.

Examples
--------
view_results Results_*/

view_results Results_*/ --tree-view

view_results Results_*/ --stats-view --groups-view
"""

import argparse
import os
import sys

from orthofinder.run.views import write_views
from orthofinder.utils.printing import print_pages


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="view_results",
        description="Write the pages to view OrthoFinder results in a web browser, from the results "
                    "files of a finished run (as a run with --view). With no option, all three pages.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  view_results Results_*/\n  view_results Results_*/ --tree-view\n"
               "  view_results Results_*/ --stats-view --groups-view")
    parser.add_argument("results_dir", metavar="INPUT", help="an OrthoFinder results directory (Results_<name>)")
    parser.add_argument("--tree-view", action="store_true", help="the trees: Resolved_Gene_Trees/Resolved_Gene_Trees.html")
    parser.add_argument("--stats-view", action="store_true",
                        help="the statistics, with plots: Comparative_Genomics_Statistics/Statistics.html")
    parser.add_argument("--groups-view", action="store_true",
                        help="the genes of the orthogroups and HOGs: Orthogroups/Groups.html")
    parser.add_argument("--view", action="store_true", help="all three pages (the default)")
    a = parser.parse_args(args)
    results_dir = os.path.normpath(a.results_dir)
    if not os.path.isdir(os.path.join(results_dir, "Orthogroups")):
        parser.error(f"not an OrthoFinder results directory (no Orthogroups directory): {a.results_dir}")
    every = a.view or not (a.tree_view or a.stats_view or a.groups_view)
    pages, hub, there = write_views(results_dir, trees=every or a.tree_view, stats=every or a.stats_view,
                                    groups=every or a.groups_view)
    if not pages:
        print(f"No pages written: {results_dir} has no resolved gene trees or orthogroups yet", file=sys.stderr)
        return 1
    print_pages(there if hub else pages, hub)  # with OrthoFinder_Results.html: all its pages
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
