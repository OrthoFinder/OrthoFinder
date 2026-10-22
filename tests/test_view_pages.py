"""The pages to view the results (--view) and OrthoFinder_Results.html, a tab for each."""
import json
import re
import shutil
import subprocess

import pytest

from orthofinder.utils import view_pages


def test_present_and_links(tmp_path):
    assert view_pages.present(str(tmp_path)) == set()
    (tmp_path / "Orthogroups").mkdir()
    (tmp_path / "Orthogroups" / "Groups.html").write_text("")
    assert view_pages.present(str(tmp_path)) == {"groups"}  # the pages that are there
    assert view_pages.present(str(tmp_path), {"trees", "stats"}) == {"trees", "stats"}  # those a run is writing
    assert view_pages.link(str(tmp_path), "groups") == "../Orthogroups/Groups.html"
    assert view_pages.link(str(tmp_path), "stats") is None


def test_hub(tmp_path):
    results = tmp_path / "Results_Oct10"
    results.mkdir()
    fn = view_pages.write_hub(str(results), {"groups", "trees"})
    assert fn == str(results / "OrthoFinder_Results.html")
    page = (results / "OrthoFinder_Results.html").read_text()
    data = json.loads(re.search(r"const D = (.*);\n", page).group(1))
    # a tab for each page, in the same order always, with its path from the results directory
    assert data["pages"] == [["trees", "Gene trees", "Resolved_Gene_Trees/Resolved_Gene_Trees.html"],
                             ["groups", "Orthogroups & HOGs", "Orthogroups/Groups.html"]]
    assert data["results"] == "Results_Oct10"


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to check the scripts")
def test_scripts_are_valid_javascript(tmp_path):
    view_pages.write_hub(str(tmp_path), {"trees", "stats"})
    js = tmp_path / "hub.js"
    js.write_text((tmp_path / "OrthoFinder_Results.html").read_text().split("<script>")[1].split("</script>")[0])
    nav = tmp_path / "nav.js"
    nav.write_text("const D = {nav: []};\n" + view_pages.NAV_JS)
    for f in (js, nav):
        result = subprocess.run(["node", "--check", str(f)], capture_output=True, text=True, check=False)
        assert result.returncode == 0, result.stderr


def test_command(tmp_path, capsys):
    assert view_pages.main([str(tmp_path)]) == 1  # no pages
    (tmp_path / "Comparative_Genomics_Statistics").mkdir()
    (tmp_path / "Comparative_Genomics_Statistics" / "Statistics.html").write_text("")
    assert view_pages.main([str(tmp_path)]) == 0
    assert "OrthoFinder results page (Statistics)" in capsys.readouterr().out
    assert (tmp_path / "OrthoFinder_Results.html").exists()


def test_serve(tmp_path):
    # the pages over HTTP (e.g. results on a cluster): the files of the results directory, compressed
    import gzip
    import threading
    import urllib.error
    import urllib.request
    (tmp_path / "Orthogroups").mkdir()
    (tmp_path / "Orthogroups" / "Groups.html").write_text("<p>groups</p>" * 100)
    (tmp_path / "secret.txt").write_text("")
    server = view_pages.make_server(str(tmp_path / "Orthogroups" / ".."), port=0)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        with urllib.request.urlopen(url + "/") as r:  # to the page there (one: no OrthoFinder_Results.html)
            assert r.url == url + "/Orthogroups/Groups.html"
        with urllib.request.urlopen(urllib.request.Request(url + "/Orthogroups/Groups.html",
                                                           headers={"Accept-Encoding": "gzip"})) as r:
            assert r.headers["Content-Encoding"] == "gzip"
            assert gzip.decompress(r.read()).decode() == "<p>groups</p>" * 100
        for path in ("/Orthogroups/", "/../" + tmp_path.name + "/secret.txt", "/missing.js"):
            with pytest.raises(urllib.error.HTTPError) as e:  # no lists of files, nothing outside it
                urllib.request.urlopen(url + path)
            assert e.value.code == 404, path
        with urllib.request.urlopen(url + "/favicon.ico") as r:
            assert r.status == 204
        # the port taken: the next one
        other = view_pages.make_server(str(tmp_path), port=server.server_address[1])
        assert other.server_address[1] == server.server_address[1] + 1
        other.server_close()
    finally:
        server.shutdown()
        server.server_close()
