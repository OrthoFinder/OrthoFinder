# Copied from src/orthofinder/utils/printing.py by tools/view_results/make_copies.py: edit that file, not this one.
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
            print(re.sub(r"\[/?[a-z_0-9 ]+\]", "", text))

    printer = _Printer()


def open_commands(path, system=None):
    """
    Commands to open a file in a web browser on this system (or system, as
    platform.system() gives it), with the path quoted for its shell. On Linux:
    the default browser (xdg-open), or Chrome or Firefox, their messages hidden
    (a browser that is already open prints errors, e.g. "Broken pipe", that
    are not about the page).
    """
    import platform
    import shlex
    system = system or platform.system()
    if system == "Windows":
        return ['start "" "%s"' % path]
    quoted = shlex.quote(path)
    if system == "Darwin":
        return ["open %s" % quoted]
    return ["xdg-open %s" % quoted, "google-chrome %s 2>/dev/null" % quoted, "firefox %s 2>/dev/null" % quoted]


def print_pages(pages, hub=None):
    """
    The web pages a run wrote (--tree-view, --stats-view, --groups-view),
    [(what, path)], and how to open them in a web browser (the commands printed
    whole, to copy). hub: OrthoFinder_Results.html, a tab for each page (with two
    or more): the one page to open.
    """
    from pathlib import Path
    if not pages:
        return
    if hub:
        names = escape(", ".join(what for what, _ in pages))
        printer.print(f"\n[bold]The results in a web browser[/bold] ({names}, a tab for each):")
        path = hub
    else:
        printer.print(f"\n[bold]The results in a web browser[/bold] ({escape(pages[0][0])}):")
        path = pages[0][1]
    path = os.path.abspath(path)
    printer.print(f"    [dark_cyan]{escape(path)}[/dark_cyan]", soft_wrap=True)
    printer.print("Open it in a web browser, or paste this address into one:")
    printer.print("    " + escape(Path(path).as_uri()), soft_wrap=True)
    printer.print("or open it with, e.g.:")
    for command in open_commands(path):
        printer.print("    " + escape(command), soft_wrap=True)
    printer.print("To look at them on another computer (e.g. Windows or macOS), copy the whole results directory: "
                  "each page reads its data from the folders next to it.")
