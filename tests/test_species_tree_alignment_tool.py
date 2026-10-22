"""tools/species_tree_alignment.py writes the same files as the run that inferred the species tree."""

import filecmp
import glob
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")


@pytest.mark.slow
def test_species_tree_alignment_tool_reproduces_the_run(tmp_path):
    fasta_dir = tmp_path / "ExampleData"
    fasta_dir.mkdir()
    for fn in glob.glob(os.path.join(ROOT, "ExampleData", "*.faa")):
        shutil.copy(fn, fasta_dir)
    subprocess.run([sys.executable, os.path.join(ROOT, "orthofinder.py"), "-f", str(fasta_dir), "-t", "2", "-n", "t"],
                   check=True, capture_output=True)
    results = fasta_dir / "OrthoFinder" / "Results_t"
    out = tmp_path / "out"
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "species_tree_alignment.py"), str(results), "-o", str(out)],
                   check=True, capture_output=True)
    assert filecmp.cmp(results / "WorkingDirectory" / "Alignments_ids" / "SpeciesTreeAlignment.fa",
                       out / "SpeciesTreeAlignment.fa", shallow=False)
    with open(results / "Species_Tree" / "Orthogroups_for_concatenated_alignment.txt") as a, \
            open(out / "Orthogroups_for_concatenated_alignment.txt") as b:
        assert sorted(a) == sorted(b)
