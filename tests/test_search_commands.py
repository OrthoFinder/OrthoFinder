"""The results files a search command writes, found in the command (search.runner)."""
import pytest

from orthofinder.search import runner


@pytest.mark.parametrize("command, found", [
    ("diamond blastp -d /w/diamondDBSpecies1 -q /w/Species0.fa -o /w/Blast0_1.txt", {("/w/", "Blast0_1.txt")}),
    ("blastp -query /w/Species0.fa -db /w/BlastDBSpecies1 -out Blast0_1.txt", {("", "Blast0_1.txt")}),
    ("diamond blastp --out=/w/Blast5_-1.txt --compress 1", {("/w/", "Blast5_-1.txt")}),
    ("blastp -outfmt 6 -query /w/Species3.fa | gzip > '/a dir/w/Blast3_4.txt.gz'", {("dir/w/", "Blast3_4.txt")}),
    # MMseqs2: the temporary folder tmpBlast<i>_<j>.txt is not a results file
    ("mmseqs search /w/mmseqsDBSpecies0 /w/mmseqsDBSpecies3.fa /w/Blast0_3.txt.db /w/tmpBlast0_3.txt ; "
     "mmseqs convertalis /w/mmseqsDBSpecies0 /w/mmseqsDBSpecies3.fa /w/Blast0_3.txt.db /w/Blast0_3.txt",
     {("/w/", "Blast0_3.txt")}),
    ("x -o /w/NotBlast0_1.txt", set()),
])
def test_results_in_command(command, found):
    assert runner._results_in_command(command) == found


def test_search_commands_for_moved_results(tmp_path):
    # Commands saved with the results directory's old path write to the one being used
    wd = str(tmp_path) + "/"
    commands = ["diamond blastp -d /old/WorkingDirectory/diamondDBSpecies1 -o /old/WorkingDirectory/Blast0_1.txt",
                "diamond blastp -d /old/WorkingDirectory/diamondDBSpecies0 -o /old/WorkingDirectory/Blast1_0.txt"]
    selected, missing, moved = runner.search_commands_for(commands, [wd + "Blast1_0.txt", wd + "Blast2_0.txt"], wd)
    assert selected == ["diamond blastp -d %sdiamondDBSpecies0 -o %sBlast1_0.txt" % (wd, wd)]
    assert missing == [wd + "Blast2_0.txt"] and list(moved) == ["/old/WorkingDirectory/"]
