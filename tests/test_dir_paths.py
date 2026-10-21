"""Directory paths are kept ending in the separator, wherever they come from."""

import os
import types

from orthofinder.run import previous_run
from orthofinder.utils import file_handler, paths


def test_as_dir():
    assert paths.as_dir("/a/b") == "/a/b" + os.sep
    assert paths.as_dir("/a/b" + os.sep) == "/a/b" + os.sep


def test_file_handler_directories_end_in_separator():
    h = file_handler._FileHandler()
    h.wd_current = "/r/Results_x/WorkingDirectory"
    h.wd_trees = "/r/Results_x/WorkingDirectory"
    h.rd1 = "/r/Results_x"
    h.wd_base = ["/r/Results_x/WorkingDirectory", "/r/Results_y/WorkingDirectory/"]
    assert h.wd_current == "/r/Results_x/WorkingDirectory/"
    assert h.wd_trees == "/r/Results_x/WorkingDirectory/"
    assert h.rd1 == "/r/Results_x/"
    assert h.wd_base == ["/r/Results_x/WorkingDirectory/", "/r/Results_y/WorkingDirectory/"]
    assert h.get_species_fasta_fn(0, q_for_creation=True) == "/r/Results_x/WorkingDirectory/Species0.fa"


def test_log_working_directory_without_separator(tmp_path):
    """
    Results of an older version (no run_state.json), moved, whose Log.txt
    records WorkingDirectory_Base without the final separator: -b on the
    results directory gave .../WorkingDirectorySpeciesIDs.txt
    """
    results = tmp_path / "OrthoFinder" / "Results_Aug27"
    wd = results / "WorkingDirectory"
    wd.mkdir(parents=True)
    (wd / "SpeciesIDs.txt").write_text("0: a.fa\n")
    clusters = "clusters_OrthoFinder_I1.5.txt_id_pairs.txt"
    (wd / clusters).write_text("")
    (results / "Log.txt").write_text(
        "WorkingDirectory_Base: /old/place/Results_Aug27/WorkingDirectory\n"
        "FN_Orthogroups: /old/place/Results_Aug27/WorkingDirectory/%s\n" % clusters
    )
    options = types.SimpleNamespace(q_start_from_fasta=True, q_start_from_blast=True)
    pfl = previous_run.PreviousFilesLocatorNew(options, str(results))
    assert pfl.wd_base_prev == [str(wd) + os.sep]
    assert pfl.wd_trees == str(wd) + os.sep
