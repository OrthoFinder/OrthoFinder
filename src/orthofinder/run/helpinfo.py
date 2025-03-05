from ..utils import util
from ..utils.util import printer
from .. import g_mclInflation, nThreadsDefault

try:
    from rich import print
except ImportError:
    ...


def PrintHelp(prog_caller):
    msa_ops = prog_caller.ListMSAMethods()
    tree_ops = prog_caller.ListTreeMethods()
    search_ops = prog_caller.ListSearchMethods()

    print("")
    print("SIMPLE USAGE:")
    print("Run full [dark_goldenrod]OrthoFinder[/dark_goldenrod] analysis on [orange3]FASTA[/orange3] format proteomes in <dir>")
    print("  orthofinder [options] -f <dir>")
    # print("")
    # print("Add new species in <dir1> to previous run in <dir2> and run new analysis")
    # print("  orthofinder [options] -f <dir1> -b <dir2>")
    print("")
    print("To assign species from <dir1> to existing [dark_goldenrod]OrthoFinder[/dark_goldenrod] orthogroups in <dir2>")
    print("  orthofinder [options] --assign <dir1> --core <dir2>")
    print("")

    print("OPTIONS:")
    print(
        " -t <int>                Number of parallel sequence search threads [Default = %d]"
        % nThreadsDefault
    )
    print(" -a <int>                Number of parallel analysis threads")
    print(" -d                      Input is DNA sequences")                ### is this still an option??
    print(
        ' -M <txt>                Method for gene tree inference. Options "dendroblast" & "msa"'
    )
    print("                         [Default = msa]")
    print(" -S <txt>                Sequence search program [Default = diamond]")
    print("                         Options: " + ", ".join(["blast"] + search_ops))
    print(' -A <txt>                MSA program, requires "-M msa" [Default = FAMSA]')      ##edited famsa
    print("                         Options: " + ", ".join(msa_ops))
    print(
        ' -T <txt>                Tree inference method, requires "-M msa" [Default = FastTree]'    ## corrected to FastTree
    )
    print("                         Options: " + ", ".join(tree_ops))
    #    print(" -R <txt>                Tree reconciliation method [Default = of_recon]")
    #    print("                      Options: of_recon, dlcpar, dlcpar_convergedsearch")
    print(" -s <file>               User-specified rooted species tree")
    print(
        " -I <int>                MCL inflation parameter [Default = %0.1f]"
        % g_mclInflation
    )
    print(
        " --save-space            Only create one compressed orthologs file per species"           ## is this an option? i thought it did this already?
    )
    print(" -x <file>               Info for outputting results in OrthoXML format")
    print(" -p <dir>                Write the temporary pickle files to <dir>")
    print(" -1                      Only perform one-way sequence search")
    print(" -X                      Don't add species names to sequence IDs")
    print(
        " -y                      Split paralogous clades below root of a HOG into separate HOGs"
    )
    print(
        " -z                      Don't trim MSAs (columns>=90% gap, min. alignment length 500)"
    )
    print(" -n <txt>                Name to append to the results directory")
    print(" -o <txt>                Non-default results directory")
    print(" -h                      Print this help text")
    print(
        " -efn                    Extend the output directory name with the name of the scoring matrix, gap penalties, search program, MSA program and tree program"
    )
    print(" --matrix <txt>          Scoring matrix allowed by DIAMOND")
    print(" --custom-matrix <txt>   Custom scoring matrix")

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
