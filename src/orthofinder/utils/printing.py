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

    def print(self, *messages, style: str = "default", sep: str = " ", end: str = "\n"):
        message = sep.join(map(str, messages))
        if style in ("error", "warning"):
            run_logging.RunLogger.message(message, level=style.upper())
        self.console.print(message, style=style, end=end)


printer = OrthoPrinter()


def stderr_exempt(stderr):
    ok_line_starts = {
        "diamond v",
        "Licensed under the GNU GPL",
        "Check http://github.com/",
    }
    try:
        stderr = stderr.decode()
    except (UnicodeDecodeError, AttributeError):
        stderr = stderr.encode()
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


def print_run_end(d_results, continue_option=None, continue_from=None):
    """
    The end of a run, also of one stopped early (-op, -og, -ogh, -ogt, -ost):
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
    with open(d + "Citation.txt", "w") as outfile:
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


def number_open_files_exception_advice(n_species, q_at_trees):
    """
    Prints advice for user on "IOError: [Errno 24] Too many open files" exception
    Args:
        n_species - the number of species in the analysis
        q_at_trees - has this error occurred at the orthologs from trees stage
    """
    # parallel_task_manager.run_command("ulimit -Hn")
    n_req = n_species + 100
    msg = (
        "\nERROR: The system limits on the number of files a process can open is too low. For %d species \
OrthoFinder needs to be able to open at least r=%d files. Please increase the limit and restart OrthoFinder\n\
1. Check the hard and soft limits on the number of open files for your system:\n\
    $ ulimit -Hn\n\
    $ ulimit -Sn\n\
2. If hard limit, h > r already, then you just need to increase the soft limit:\n\
    $ ulimit -n %d\n\
3. Alternatively, if h < r then you need to edit the file '/etc/security/limits.conf', this requires root privileges. \
To increase the limit to %d for user  called 'emms' add the lines:\n\
    emms hard nofile %d\n\
    emms soft nofile %d\n"
        % (n_species, n_req, n_req, n_req, n_req, n_req)
    )
    msg += "    (edit these lines to match your username)\n\
4. Check the limit has now been updated (if you changed the hard limit you'll need to open a new session and confirm it's updated):\n\
    $ ulimit -Sn"

    if q_at_trees:
        msg_part_2 = "5. Once the limit is updated restart OrthoFinder 'from trees' using the '-ft' command"
    else:
        msg_part_2 = (
            "5. Once the limit is updated restart OrthoFinder with the original command"
        )
    msg_part_3 = (
        "\nFor full details see: https://github.com/davidemms/OrthoFinder/issues/384"
    )
    print(msg + "\n" + msg_part_2 + "\n" + msg_part_3 + "\n")
    print(
        "Since this issue OrthoFinder has been updated to only require O(n) open files for n species rather than O(n^2)."
    )
    print("Please follow the above advice to ensure that it is able to.")


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
