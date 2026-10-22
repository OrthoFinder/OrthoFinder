"""Removing a species from an analysis: commented out in SpeciesIDs.txt, then -b on its results."""

import filecmp
import glob
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
REMOVED = "Mycoplasma_genitalium"


def run_orthofinder(*args):
    subprocess.run([sys.executable, os.path.join(ROOT, "orthofinder.py"), *args, "-t", "2"],
                   check=True, capture_output=True, text=True)


@pytest.mark.slow
def test_removing_a_species_gives_the_results_of_the_other_species(tmp_path):
    """
    The -b restart infers the orthogroups again, in the results it is given: the files of
    the earlier orthogroups (sequences, alignments, gene trees) must not be used for the
    new ones, which with the same numbers are other genes.
    """
    all_species, others = tmp_path / "all", tmp_path / "others"
    all_species.mkdir()
    others.mkdir()
    for fn in glob.glob(os.path.join(ROOT, "ExampleData", "*.faa")):
        shutil.copy(fn, all_species)
        if REMOVED not in fn:
            shutil.copy(fn, others)
    run_orthofinder("-f", str(all_species), "-n", "all")
    run_orthofinder("-f", str(others), "-n", "others")
    results = all_species / "OrthoFinder" / "Results_all"

    species_ids = results / "WorkingDirectory" / "SpeciesIDs.txt"
    lines = species_ids.read_text().splitlines(True)
    species_ids.write_text("".join("#" + l if REMOVED in l else l for l in lines))
    run_orthofinder("-b", str(results))

    for name in ("Orthogroups/Orthogroups.tsv", "Orthogroups/Orthogroups_UnassignedGenes.tsv",
                 "Gene_Duplication_Events/Duplications.tsv"):
        assert filecmp.cmp(results / name, others / "OrthoFinder" / "Results_others" / name, shallow=False), name
    assert REMOVED not in (results / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").read_text()
