"""The help (-h): a section for each group of options (cli_options)."""
try:
    from rich import print
except ImportError:
    ...
from rich.table import Table
from rich.console import Console

from ..utils import printing
from . import cli_options

width = 22  # of the column of the options


def option_rows(option, prog_caller):
    """The rows of an option in the help: (option, text); the values it takes listed under it."""
    lines = option["help"].split("\n")
    rows = [(cli_options.help_label(option), lines[0])] + [("", line) for line in lines[1:]]
    if "choices" in option:
        rows.append(("", "Options: " + ", ".join(option["choices"](prog_caller))))
    return rows


def print_help(prog_caller):
    console = Console()

    print("")
    print("[bold]SIMPLE USAGE:[/bold]")
    print(
        " Run full [dark_goldenrod]OrthoFinder[/dark_goldenrod] analysis on [red1]FASTA[/red1] format proteomes in <dir>"
    )
    print("   orthofinder \\[options] -f <dir>")
    print("")
    print(
        " To assign species from <dir1> to existing [dark_goldenrod]OrthoFinder[/dark_goldenrod] orthogroups in <dir2>"
    )
    print("   orthofinder \\[options] --assign <dir1> --core <dir2>")
    print("")

    for group in cli_options.OPTION_GROUPS:
        if group.get("hidden"):   # e.g. the development options
            continue
        shown = [o for o in group["options"] if o["help"] is not None and cli_options.is_available(o)]
        if not shown:   # (e.g. only preview options)
            continue
        table = Table(show_header=False, box=None, expand=False)
        table.add_column("Option", justify="left", width=width, no_wrap=False, overflow="fold")
        table.add_column("Description", justify="left", overflow="fold")
        for option in shown:
            for row in option_rows(option, prog_caller):
                table.add_row(*row)
        console.print("[bold]%s:[/bold]" % group["title"])
        console.print(table)
        console.print()

    console.print("[bold]LICENSE:[/bold]")
    print(
        " Distributed under the [dodger_blue1]GNU General Public License (GPLv3)[/dodger_blue1]. See LICENSE"
    )
    printing.print_citation()
