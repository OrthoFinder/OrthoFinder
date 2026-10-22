"""Printing to the console and the run log: the printer, run summaries and progress bars."""

import os

from ..citation import citation
from . import logging as run_logging

try:
    from rich import print, progress
    from rich.console import Console
    from rich.theme import Theme

except ImportError:
    ...
try:
    width = os.get_terminal_size().columns
except OSError as e:
    width = 80


class OrthoPrinter:
    def __init__(self):
        self.console = Console(
            theme=Theme(
                {
                    "info": "deep_sky_blue2",
                    "path": "green4",
                    "orthofinder": "dark_goldenrod",
                    "other_methods": "dark_orange3",
                    "version": "dark_cyan",
                    "success": "bold green",
                    "error": "bold red",
                    "warning": "bold yellow",
                    "default": "default on default",  # bright_white
                }
            ),
            force_terminal=True,
        )

    def print(self, *messages, style: str = "default", sep: str = " ", end: str = "\n", soft_wrap: bool = False):
        """soft_wrap: do not break the lines (e.g. so that a command can be copied whole)"""
        message = sep.join(map(str, messages))
        if style in ("error", "warning"):
            run_logging.RunLogger.message(message, level=style.upper())
        self.console.print(message, style=style, end=end, soft_wrap=soft_wrap)


printer = OrthoPrinter()


def stderr_exempt(stderr):
    ok_line_starts = {
        "diamond v",
        "Licensed under the GNU GPL",
        "Check http://github.com/",
    }
    if isinstance(stderr, bytes):
        # (a program's output need not be UTF-8, e.g. in a Windows console's code page)
        stderr = stderr.decode(errors="replace")
    lines = stderr.split("\n")
    for line in lines:
        if line.rstrip() == "":
            continue
        if any(line.startswith(x) for x in ok_line_starts):
            continue
        return False
    return True


def print_continue_hint(option, directory, what):
    """After a run stopped early (-op, -og, -ogt, -ost): how to continue it."""
    from rich.markup import escape
    printer.print("\nTo continue this analysis from %s, run:" % what)
    printer.print("    [bold]orthofinder %s [dark_cyan]%s[/dark_cyan][/bold]" % (option, escape(directory)))


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
    from rich.markup import escape
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


def print_run_end(d_results, continue_option=None, continue_from=None):
    """
    The end of a run, also of one stopped early (-op, -og, -ogt, -ost):
    the results directory, how to continue a stopped run, and the citation.
    """
    from rich.markup import escape
    d_results = os.path.normpath(d_results) + os.path.sep
    printer.print("\n[bold]Results directory:[/bold]")
    printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(d_results))
    if continue_option is not None:
        print_continue_hint(continue_option, d_results, continue_from)
    print_citation(d_results)


def print_restart_info(options, directory):
    """At the start of a restarted run (-b, -fg, -fgt, -fst): what it continues from."""
    from rich.markup import escape
    if options.q_start_from_blast and not options.q_start_from_fasta:
        what = "the sequence search results"
    elif options.q_start_from_groups:
        what = "the orthogroups"
    elif options.q_start_from_trees:
        what = "the gene trees"
    elif options.q_start_from_species_trees:
        what = "the species tree"
    else:
        return
    printer.print("\n[bold]Restarting[/bold] from %s in:" % what)
    printer.print("    [dark_cyan]%s[/dark_cyan]" % escape(directory))


def write_citation(d):
    with open(d + "Citation.txt", "w", encoding="utf-8", newline="\n") as outfile:
        outfile.write(citation)


def print_citation(d=None):
    if d is not None:
        write_citation(d)
    # printer.print(print_citation)
    printer.print("\n[bold]CITATION:[/bold]")
    printer.print(
        " When publishing work that uses [dark_goldenrod]OrthoFinder[/dark_goldenrod] please cite:"
    )
    printer.print(
        " Emms D.M., Liu Y., Belcher L., Holmes J. & Kelly S. (2026), "
        "Nature Methods 23:1327–1333, ",
        end="",
    )
    printer.print(
        "[dark_cyan]"
        "[link=https://doi.org/10.1038/s41592-026-03126-6]"
        "https://doi.org/10.1038/s41592-026-03126-6"
        "[/link]\n"
    )
    # printer.print(" Emms D.M. & Kelly S. (2019), Genome Biology 20:238\n")

    printer.print(" If you use the species tree in your work then please also cite:")
    printer.print(" Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278")
    printer.print(" Emms D.M. & Kelly S. (2018), bioRxiv ", end="")
    printer.print(
        "[dark_cyan][link=https://doi.org/10.1101/267914]https://doi.org/10.1101/267914[/link]"
    )


def print_underline(text, q_heavy=False):
    run_logging.RunLogger.message(text)
    n = len(text)
    if text.startswith("\n"):
        n -= 1
    print(("\n[bold]" + text + "[/bold]") if q_heavy else ("\n" + text))
    print((("=" if q_heavy else "-") * n))


# The open files one OrthoFinder process needs, whatever the number of species (the output
# writers keep a bounded number open), and per thread (the pipes of its worker processes
# and commands): measured, 64 files for 8 threads and 128 for 64, with this headroom
OPEN_FILES_BASE = 128
OPEN_FILES_PER_THREAD = 2


def open_files_needed(n_threads):
    return OPEN_FILES_BASE + OPEN_FILES_PER_THREAD * max(1, int(n_threads))


def number_open_files_exception_advice(n_species, q_at_trees, n_req=None):
    """
    Prints advice for user on "IOError: [Errno 24] Too many open files" exception
    Args:
        n_species - the number of species in the analysis (the files needed do not depend on it)
        q_at_trees - has this error occurred at the orthologs from trees stage
        n_req - the files needed (open_files_needed), if known
    """
    if n_req is None:
        n_req = open_files_needed(1)
    msg = (
        "\nERROR: The system limit on the number of files a process can open is too low. OrthoFinder "
        "needs to be able to open at least r=%d files (with these threads, for any number of species: "
        "%d here). Please increase the limit and restart OrthoFinder\n"
        "1. Check the hard and soft limits on the number of open files for your system:\n"
        "    $ ulimit -Hn\n"
        "    $ ulimit -Sn\n"
        "2. If the hard limit h >= r, increase the soft limit:\n"
        "    $ ulimit -n %d\n"
        "3. Otherwise ask your system administrator to raise the hard limit (on Linux, in "
        "'/etc/security/limits.conf', which needs root privileges), e.g. for a user called 'user':\n"
        "    user hard nofile %d\n"
        "    user soft nofile %d\n"
        "4. Check the limit has been updated (after changing the hard limit, in a new session):\n"
        "    $ ulimit -Sn"
        % (n_req, n_species, n_req, n_req, n_req)
    )
    if q_at_trees:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder 'from trees' using the '-ft' command"
    else:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder with the original command"
    print(msg + "\n" + msg_part_2 + "\n")


def get_progressbar(len_task, visible=True):
    progressbar = progress.Progress(
        progress.TextColumn("[progress.description]{task.description}"),
        progress.BarColumn(bar_width=width // 2),
        progress.SpinnerColumn(),
        progress.MofNCompleteColumn(),
        progress.TimeElapsedColumn(),
        transient=False,
        # progress.TextColumn("{task.completed}/{task.total}")
    )

    task = progressbar.add_task(
            "[yellow]Processing...", total=len_task, visible=visible
    )
    return progressbar, task
