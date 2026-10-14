"""tools/view_results: the pages to view the results of a finished run (as a run with --view)."""
import json
import os
import re
import subprocess
import sys

from test_tree_pages import make_results

from tools.view_results import make_copies, view_results


def page_data(path):
    return json.loads(re.search(r"const D = (.*);\n", path.read_text()).group(1))


def results(d):
    make_results(d)  # orthogroups, gene trees, HOGs, alignments, duplications
    (d / "Orthogroups" / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\tSp_A\tSp_B\tTotal\nOG0000000\t2\t1\t3\nOG0000001\t1\t2\t3\n")


def test_all_pages(tmp_path, capsys):
    results(tmp_path)
    assert view_results.main([str(tmp_path)]) == 0  # no option: all three
    for page in ("Resolved_Gene_Trees/Resolved_Gene_Trees.html", "Comparative_Genomics_Statistics/Statistics.html",
                 "Orthogroups/Groups.html", "OrthoFinder_Results.html"):
        assert (tmp_path / page).exists(), page
    assert "OrthoFinder_Results.html" in capsys.readouterr().out
    hub = page_data(tmp_path / "OrthoFinder_Results.html")
    assert [k for k, _, _ in hub["pages"]] == ["trees", "stats", "groups"]


def test_one_page_linked_to_those_there(tmp_path, capsys):
    results(tmp_path)
    view_results.main([str(tmp_path), "--tree-view"])
    assert not (tmp_path / "OrthoFinder_Results.html").exists()  # one page: nothing to put in tabs
    tree_page = tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html"
    written = tree_page.stat().st_mtime_ns
    view_results.main([str(tmp_path), "--stats-view"])
    assert tree_page.stat().st_mtime_ns == written  # the page there is not written again
    stats = page_data(tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html")
    assert stats["tree_page"] == "../Resolved_Gene_Trees/Resolved_Gene_Trees.html"  # but linked to
    hub = page_data(tmp_path / "OrthoFinder_Results.html")  # two pages now: a tab for each
    assert [k for k, _, _ in hub["pages"]] == ["trees", "stats"]
    assert "Gene trees, Statistics" in capsys.readouterr().out.replace("\n", " ")


def test_not_results(tmp_path, capsys):
    import pytest
    with pytest.raises(SystemExit):
        view_results.main([str(tmp_path)])
    assert "not an OrthoFinder results directory" in capsys.readouterr().err


def test_pages_with_the_default_encoding_of_windows(tmp_path, monkeypatch):
    """
    On Windows, a file opened as text without an encoding is in the code page
    (cp1252), not UTF-8: the pages are written as UTF-8 (as they say they are),
    from the results (UTF-8) and their text that cp1252 does not have (e.g. ▸).
    """
    import builtins
    results(tmp_path)
    for fn in tmp_path.rglob("*"):   # a species name with an accent
        if fn.is_file():
            fn.write_text(fn.read_text(encoding="utf-8").replace("Sp_A", "Sp_Á"), encoding="utf-8")
    real_open = builtins.open

    def windows_open(file, mode="r", buffering=-1, encoding=None, *args, **kwargs):
        if "b" not in mode and encoding is None:
            encoding = "cp1252"
        return real_open(file, mode, buffering, encoding, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", windows_open)
    assert view_results.main([str(tmp_path)]) == 0
    monkeypatch.undo()
    for page in ("Resolved_Gene_Trees/Resolved_Gene_Trees.html", "Comparative_Genomics_Statistics/Statistics.html",
                 "Orthogroups/Groups.html", "OrthoFinder_Results.html"):
        text = (tmp_path / page).read_text(encoding="utf-8")   # UTF-8, as its <meta charset> says
        assert '<meta charset="utf-8">' in text
    assert "▸" in (tmp_path / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.html").read_text(encoding="utf-8")
    assert "Sp_Á" in json.dumps(page_data(tmp_path / "Orthogroups" / "Groups.html"), ensure_ascii=False)


def test_commands_to_open_the_pages():
    from orthofinder.utils.printing import open_commands
    # quoted for each system's shell: paths with spaces are usual on macOS and Windows
    assert open_commands(r"C:\Users\Jane Doe\Results\OrthoFinder_Results.html", "Windows") == \
        [r'start "" "C:\Users\Jane Doe\Results\OrthoFinder_Results.html"']
    assert open_commands("/Users/jane doe/Results/OrthoFinder_Results.html", "Darwin") == \
        ["open '/Users/jane doe/Results/OrthoFinder_Results.html'"]
    assert open_commands("/home/jane/Results/OrthoFinder_Results.html", "Linux") == \
        ["xdg-open /home/jane/Results/OrthoFinder_Results.html",
         "google-chrome /home/jane/Results/OrthoFinder_Results.html 2>/dev/null",
         "firefox /home/jane/Results/OrthoFinder_Results.html 2>/dev/null"]


def test_copies_of_the_page_code_are_current():
    """tools/view_results/orthofinder_pages: if not, run python tools/view_results/make_copies.py"""
    directory = os.path.join(make_copies.HERE, make_copies.PACKAGE)
    with open(os.path.join(directory, "__init__.py")) as infile:
        version = re.search(r'__version__ = "(.*)"', infile.read()).group(1)
    for name, text in make_copies.copies(version).items():
        with open(os.path.join(directory, name), encoding="utf-8") as infile:
            assert infile.read() == text, "%s is not current: run python tools/view_results/make_copies.py" % name


def test_works_without_orthofinder(tmp_path):
    """The directory on its own (e.g. downloaded as a zip): it does not import OrthoFinder."""
    results(tmp_path / "Results")
    zip_fn = tmp_path / "view_results.zip"
    make_copies.write_zip(str(zip_fn))
    import zipfile
    zipfile.ZipFile(zip_fn).extractall(tmp_path / "download")
    script = tmp_path / "download" / "view_results" / "view_results.py"
    # Python without the installed OrthoFinder (an editable install is a finder on sys.meta_path)
    code = ("import runpy, sys\n"
            "sys.meta_path = [f for f in sys.meta_path if 'editable' not in repr(f)]\n"
            "sys.path = [p for p in sys.path if not os.path.isdir(os.path.join(p, 'orthofinder'))]\n"
            "sys.path.insert(0, %r)\n"
            "sys.argv = [%r, %r]\n"
            "runpy.run_path(%r, run_name='__main__')\n"
            "assert 'orthofinder' not in sys.modules\n") % (str(script.parent), str(script), str(tmp_path / "Results"),
                                                               str(script))
    run = subprocess.run([sys.executable, "-c", "import os\n" + code], capture_output=True, text=True,
                         cwd=str(tmp_path))
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "Results" / "OrthoFinder_Results.html").exists()
