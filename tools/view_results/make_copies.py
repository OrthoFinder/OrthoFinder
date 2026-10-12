#!/usr/bin/env python3

"""
Copy OrthoFinder's page code (src/orthofinder) into orthofinder_pages/, so that
this directory works on its own: view_results.py and the results are all
someone needs (with Python and ete4), not OrthoFinder.

Run it after changing the pages in src/orthofinder (a test checks that the
copies are current):

    python tools/view_results/make_copies.py

and to also write the directory as a zip file, to download:

    python tools/view_results/make_copies.py --zip dist/view_results.zip
"""

import argparse
import ast
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "src", "orthofinder")
PACKAGE = "orthofinder_pages"

# copy: the module in src/orthofinder it is copied from
MODULES = {
    "views.py": "run/views.py",
    "tree_pages.py": "phylo/tree_pages.py",
    "stats_page.py": "orthology/stats_page.py",
    "groups_page.py": "orthology/groups_page.py",
    "view_pages.py": "utils/view_pages.py",
}
# their imports of the rest of OrthoFinder, and those of the copies (each other)
IMPORTS = [
    ("from ..orthology import groups_page, stats_page\n", "from . import groups_page, stats_page\n"),
    ("from ..phylo import tree_pages\n", "from . import tree_pages\n"),
    ("from ..utils import view_pages\n", "from . import view_pages\n"),
    ("from ..phylo.tree_pages import ", "from .tree_pages import "),
    ("from .. import __version__\n", "from . import __version__\n"),
]
# from utils/printing.py: where the pages are and how to open them
PRINTING_FUNCTIONS = ("open_commands", "print_pages")

HEADER = "# Copied from src/orthofinder/%s by tools/view_results/make_copies.py: edit that file, not this one.\n"

PRINTING = '''# Copied from src/orthofinder/utils/printing.py by tools/view_results/make_copies.py: edit that file, not this one.
"""Where the pages are and how to open them, as OrthoFinder prints it (in colour with rich, if installed)."""

import os

try:
    from rich.console import Console
    from rich.markup import escape
    printer = Console()
except ImportError:
    import re

    def escape(text):
        return text

    class _Printer:
        def print(self, text="", soft_wrap=False):
            print(re.sub(r"\\[/?[a-z_0-9 ]+\\]", "", text))

    printer = _Printer()
'''

INIT = '''"""
Copies of OrthoFinder's page code, made by tools/view_results/make_copies.py
from src/orthofinder: edit those files, not these, and run make_copies.py.
"""

__version__ = "%s"   # of the OrthoFinder they were copied from
'''


def orthofinder_version():
    """The version of the OrthoFinder copied from (src/orthofinder/_version.py is made when it is built or installed)."""
    sys.path.insert(0, os.path.join(SRC, ".."))
    try:
        from orthofinder import __version__
        return __version__
    except ImportError:
        return "unknown"
    finally:
        sys.path.pop(0)


def copies(version=None):
    """{file name in orthofinder_pages: its text}."""
    out = {}
    for name, source in MODULES.items():
        with open(os.path.join(SRC, source), encoding="utf-8") as infile:
            text = infile.read()
        for old, new in IMPORTS:
            text = text.replace(old, new)
        left = re.findall(r"^\s*from \.\.\S* import .*$|^\s*import orthofinder.*$|^\s*from orthofinder\S* import .*$",
                          text, re.M)
        if left:
            raise ValueError("%s imports other parts of OrthoFinder; add them to make_copies.py: %s" % (source, left))
        out[name] = HEADER % source + text
    with open(os.path.join(SRC, "utils", "printing.py"), encoding="utf-8") as infile:
        printing = infile.read()
    functions = {node.name: ast.get_source_segment(printing, node) for node in ast.parse(printing).body
                 if isinstance(node, ast.FunctionDef) and node.name in PRINTING_FUNCTIONS}
    text = PRINTING
    for name in PRINTING_FUNCTIONS:
        body = re.sub(r"^\s*from rich\.markup import escape\n", "", functions[name] + "\n", flags=re.M)
        text += "\n\n" + body
    out["printing.py"] = text
    out["__init__.py"] = INIT % (version or orthofinder_version())
    return out


def write(directory=os.path.join(HERE, PACKAGE)):
    os.makedirs(directory, exist_ok=True)
    for name, text in copies().items():
        with open(os.path.join(directory, name), "w", encoding="utf-8") as outfile:
            outfile.write(text)


def write_zip(fn):
    """This directory as fn, in a directory view_results (without Python's caches)."""
    os.makedirs(os.path.dirname(os.path.abspath(fn)), exist_ok=True)
    with zipfile.ZipFile(fn, "w", zipfile.ZIP_DEFLATED) as z:
        for root, dirs, files in os.walk(HERE):
            dirs[:] = sorted(d for d in dirs if d != "__pycache__")
            for f in sorted(files):
                if not f.endswith(".pyc"):
                    path = os.path.join(root, f)
                    z.write(path, os.path.join("view_results", os.path.relpath(path, HERE)))


def main(args=None):
    parser = argparse.ArgumentParser(description="Copy OrthoFinder's page code into %s/" % PACKAGE)
    parser.add_argument("--zip", metavar="FILE", help="also write this directory as a zip file, e.g. dist/view_results.zip")
    a = parser.parse_args(args)
    if not os.path.isdir(SRC):
        sys.exit("make_copies.py copies from the OrthoFinder source code (src/orthofinder), which is not here: "
                 "run it in the OrthoFinder directory")
    write()
    print("Copied the page code to %s" % os.path.join(HERE, PACKAGE))
    if a.zip:
        write_zip(a.zip)
        print("Written %s" % a.zip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
