import os
from ..utils import util
from ..utils.util import printer
from .. import g_mclInflation, nThreadsDefault

try:
    from rich import print
except ImportError:
    ...

from rich.table import Table
from rich.console import Console
from rich.text import Text

try:
    width = os.get_terminal_size().columns
except OSError as e:
    width = 80

def PrintHelp(prog_caller):

    console = Console()
    table_options = Table(show_header=False, box=None, expand=False)
    table_options.add_column("Option", justify="left", no_wrap=True)
    table_options.add_column("Description", justify="left", overflow="fold")

    table_stoping = Table(show_header=False, box=None)
    table_stoping.add_column("Option", justify="left", no_wrap=True)
    table_stoping.add_column("Description", justify="left", overflow="fold")

    table_restarting = Table(show_header=False, box=None)
    table_restarting.add_column("Option", justify="left", no_wrap=True)
    table_restarting.add_column("Description", justify="left", overflow="fold")

    table_version = Table(show_header=False, box=None)
    table_version.add_column("Option", justify="left", no_wrap=True)
    table_version.add_column("Description", justify="left", overflow="fold")

    # Add two columns: one for the option, one for the description
    # "no_wrap=True" on the first column so it doesn't break mid-option.
    # "overflow='fold'" on the second column to wrap text if it's long.

    # table.add_column("", justify="left", no_wrap=True)
    # table.add_column("", justify="left", overflow="fold")

    msa_ops = prog_caller.ListMSAMethods()
    tree_ops = prog_caller.ListTreeMethods()
    search_ops = prog_caller.ListSearchMethods()

    print("")
    print("SIMPLE USAGE:")
    print("Run full [dark_goldenrod]OrthoFinder[/dark_goldenrod] analysis on [red1]FASTA[/red1] format proteomes in <dir>")
    print("  orthofinder \\[options] -f <dir>")
    # print("")
    # print("Add new species in <dir1> to previous run in <dir2> and run new analysis")
    # print("  orthofinder [options] -f <dir1> -b <dir2>")
    print("")
    print("To assign species from <dir1> to existing [dark_goldenrod]OrthoFinder[/dark_goldenrod] orthogroups in <dir2>")
    print("  orthofinder \\[options] --assign <dir1> --core <dir2>")
    print("")

    # print("OPTIONS:")
    # print(
    #     " -t <int>                Number of parallel sequence search threads [Default = %d]"
    #     % nThreadsDefault
    # )
    table_options.add_row(
        "-t <[bright_magenta]int[/bright_magenta]>",
        f"Number of parallel sequence search threads [Default = [deep_sky_blue2]{nThreadsDefault}[/deep_sky_blue2]]"
    )
    table_options.add_row(
        "-a <[bright_magenta]int[/bright_magenta]>",
        "Number of parallel analysis threads"
    )
    # # print(" -a <int>                Number of parallel analysis threads")

    # print(" -d                      Input is DNA sequences")                ### is this still an option??
    table_options.add_row(
        "-d",
        "Input is DNA sequences"
    )
    # print(
    #     ' -M <txt>                Method for gene tree inference. Options "dendroblast" & "msa"'
    # )
    table_options.add_row(
        "-M <[bright_magenta]txt[/bright_magenta]>",
        'Method for gene tree inference. Options [dark_cyan]"dendroblast"[/dark_cyan] & [dark_cyan]"msa"[/dark_cyan] [Default = [dark_cyan]msa[/dark_cyan]]'
    )
    
    # print("                         [Default = msa]")

    # print(" -S <txt>                Sequence search program [Default = diamond]")
    # print("                         Options: " + ", ".join(["blast"] + search_ops))
    table_options.add_row(
        "-S <[bright_magenta]txt[/bright_magenta]>",
        "Sequence search program [Default = [dark_cyan]diamond[/dark_cyan]]",
    )
    table_options.add_row(
        "",
        "Options: " + ", ".join(["blast"] + search_ops)
    )

    # print(' -A <txt>                MSA program, requires "-M msa" [Default = FAMSA]')      ##edited famsa
    # print("                         Options: " + ", ".join(msa_ops))
    table_options.add_row(
        "-A <[bright_magenta]txt[/bright_magenta]>",
        'MSA program, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]FAMSA[/dark_cyan]]',
    )
    table_options.add_row(
        "",
        "Options: " + ", ".join(msa_ops)
    )

    # print(
    #     ' -T <txt>                Tree inference method, requires "-M msa" [Default = FastTree]'    ## corrected to FastTree
    # )
    # print("                         Options: " + ", ".join(tree_ops))

    table_options.add_row(
        "-T <[bright_magenta]txt[/bright_magenta]>",
        'Tree inference method, requires [dark_cyan]"-M msa"[/dark_cyan] [Default = [dark_cyan]FastTree[/dark_cyan]]',
    )
    table_options.add_row(
        "",
        "Options: " + ", ".join(tree_ops)
    )

    # #    print(" -R <txt>                Tree reconciliation method [Default = of_recon]")
    # #    print("                      Options: of_recon, dlcpar, dlcpar_convergedsearch")
    
    # print(" -s <file>               User-specified rooted species tree")

    table_options.add_row(
        "-s <[bright_magenta]file[/bright_magenta]>",
        'Tree reconciliation method [Default = [dark_cyan]of_recon[/dark_cyan]]',
    )
    # print(
    #     " -I <int>                MCL inflation parameter [Default = %0.1f]"
    #     % g_mclInflation
    # )
    table_options.add_row(
        "-I <[bright_magenta]int[/bright_magenta]>",
        f'MCL inflation parameter [Default = [deep_sky_blue2]{g_mclInflation:0.1f}[/deep_sky_blue2]]',
    )

    # print(
    #     " --save-space            Only create one compressed orthologs file per species"           ## is this an option? i thought it did this already?
    # )
    table_options.add_row(
        "--save-space",
        "Only create one compressed orthologs file per species",
    )

    # print(" -x <file>               Info for outputting results in OrthoXML format")
    table_options.add_row(
        "-x <[bright_magenta]file[/bright_magenta]>",
        "Info for outputting results in OrthoXML format",
    )
    # print(" -p <dir>                Write the temporary pickle files to <dir>")
    table_options.add_row(
        "-p <[bright_magenta]dir[/bright_magenta]>",
        "Write the temporary pickle files to <[bright_magenta]dir[/bright_magenta]>",
    )

    # print(" -1                      Only perform one-way sequence search")

    table_options.add_row(
        "[deep_sky_blue2]-1[/deep_sky_blue2]",
        "Only perform one-way sequence search",
    )

    # print(" -X                      Don't add species names to sequence IDs")
    # print(
    #     " -y                      Split paralogous clades below root of a HOG into separate HOGs"
    # )
    # print(
    #     " -z                      Don't trim MSAs (columns>=90% gap, min. alignment length 500)"
    # )
    # print(" -n <txt>                Name to append to the results directory")
    # print(" -o <txt>                Non-default results directory")
    # print(" -h                      Print this help text")
    # print(
    #     " -efn                    Extend the output directory name with the name of the scoring matrix, gap penalties, search program, MSA program and tree program"
    # )
    table_options.add_row(
        "[bold]-efn[/bold]",
        (
            "Extend the output directory name with the name of the scoring matrix, "
            "gap penalties, search program, MSA program and tree program"
        ),
    )
    # print(" --matrix <txt>          Scoring matrix allowed by DIAMOND")
    # print(" --custom-matrix <txt>   Custom scoring matrix")

    console.print("[bold]OPTIONS:[/bold]")
    console.print(table_options)
    console.print()


    print("")
    print("WORKFLOW STOPPING OPTIONS:")
    print(" -op                     Stop after preparing input files for BLAST")
    print(" -og                     Stop after inferring orthogroups")                              ### I think we only go as far as doing blast??
    #print(" -os                     Stop after writing sequence files for orthogroups")
    #print("                         (requires '-M msa')")
    #print(" -oa                     Stop after inferring alignments for orthogroups")
    #print("                         (requires '-M msa')")
    #print(" -ot                     Stop after inferring gene trees for orthogroups ")

    print("")
    print("WORKFLOW RESTART COMMANDS:")
    print(
        " -b  <dir>               Start [dark_goldenrod]OrthoFinder[/dark_goldenrod] from pre-computed BLAST results in <dir>"         ### Im pretty sure this is the only intermediate we want to use
    )
    #print(
        #" -fg <dir>               Start OrthoFinder from pre-computed orthogroups in <dir>"
    #)
    #print(
        #" -ft <dir>               Start OrthoFinder from pre-computed gene trees in <dir>"
    #)

    print("")
    print("VERSION:")
    print(" -v                      Show the current version number")

    print("")
    print("LICENSE:")
    print(" Distributed under the [dodger_blue1]GNU General Public License (GPLv3)[/dodger_blue1]. See License.md")
    util.PrintCitation()
