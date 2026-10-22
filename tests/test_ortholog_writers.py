import gzip
import random
import subprocess

from orthofinder.orthology.hog_writer import HogWriter
from orthofinder.orthology.ortholog_writers import (
    GzipBuffers, LazyFileCache, OrderedCommitter, OrderedHogCommitter, sorted_lines,
)
from orthofinder.phylo import tree


def test_gzip_buffers_members_and_budget(tmp_path):
    # small members and budget, as for a large run: several members per file
    gz = GzipBuffers(budget=3000, member_bytes=1000)
    paths = [str(tmp_path / ("sp%d.tsv.gz" % i)) for i in range(5)]
    for p in paths:
        open(p, "wb").write(b"replaced")       # a file the first member replaces
        gz.headers[p] = "Orthogroup\tSp\n"
        gz.start(p)                             # as write_last does for the files a writer owns
    expected = {p: "Orthogroup\tSp\n" for p in paths}
    rng = random.Random(0)
    for i in range(2000):
        p = rng.choice(paths[:4])               # paths[4] gets no text
        line = "OG%07d\tx\t%s\n" % (i, "g" * (i % 50))
        gz.write(p, line)
        expected[p] += line
        assert gz.total <= gz.budget
    gz.close()
    for p in paths:
        assert gzip.decompress(open(p, "rb").read()).decode() == expected[p]
        assert subprocess.run(["gzip", "-t", p]).returncode == 0
    assert open(paths[0], "rb").read().count(b"\x1f\x8b\x08") > 1   # several members


def test_gzip_buffers_append_without_header(tmp_path):
    path = str(tmp_path / "b.tsv.gz")
    with gzip.open(path, "wt") as f:
        f.write("header\n")
    gz = GzipBuffers()
    gz.write(path, "row\n")
    gz.close()
    assert gzip.open(path, "rt").read() == "header\nrow\n"


def test_cache_gz_uses_header(tmp_path):
    fn = str(tmp_path / "c.tsv")
    cache = LazyFileCache()
    cache.gz_headers[fn + ".gz"] = "h\n"
    cache.write(fn, "b\n", gz=True)
    cache.write(str(tmp_path / "plain.txt"), "p\n")
    cache.write(fn, "c\n", gz=True)
    cache.close_all()
    assert gzip.open(fn + ".gz", "rt").read() == "h\nb\nc\n"
    # a file with a header and no text is still written
    cache = LazyFileCache()
    cache.gz_headers[fn + ".gz"] = "only\n"
    cache.start_gz(fn)
    cache.close_all()
    assert gzip.open(fn + ".gz", "rt").read() == "only\n"


def test_sorted_lines():
    assert sorted_lines("b\na\n") == "a\nb\n"
    assert sorted_lines("a\n") == "a\n"
    assert sorted_lines("a") == "a\n"


def test_ordered_committer_any_order(tmp_path):
    committed = []
    iogs = list(range(50))
    committer = OrderedCommitter(committed.append, iogs, max_pending_bytes=200, spill_dir=str(tmp_path))
    order = iogs[:]
    random.Random(1).shuffle(order)
    for iog in order:
        if iog % 7 == 0:
            committer.add_skip(iog)
        else:
            committer.add_result(iog, {"iog": iog, "text": "x" * iog})
    committer.assert_finished()
    committer.close()
    assert [d["iog"] for d in committed] == [i for i in iogs if i % 7 != 0]


def test_hog_ids_unique_and_in_og_order(tmp_path):
    """
    HOG IDs are given by one counter per species-tree node in the one process
    that writes the HOGs, in OG order (OrderedHogCommitter), whatever order
    the OGs are analysed in: unique and consecutive, and N0.ids has the IDs
    of N0.tsv row by row. (No check of the files is needed during a run.)
    """
    species_tree = tree.Tree("((0,1)N1,2)N0;", format=1)
    hog_writer = HogWriter(
        species_tree, ["N0", "N1"], {}, {"0": "A", "1": "B", "2": "C"}, [0, 1, 2],
        write_to_rd=True, write_output=False,
    )
    # the files of a run, without a results directory
    hog_writer.write_output = True
    for name in ("N0", "N1", "N0.ids"):
        hog_writer.hog_paths[name] = str(tmp_path / (name + ".tsv"))

    iogs = list(range(40))
    committer = OrderedHogCommitter(hog_writer, iogs, max_pending_bytes=300, spill_dir=str(tmp_path))
    order = iogs[:]
    random.Random(2).shuffle(order)
    for iog in order:
        og = "OG%07d" % iog
        if iog % 5 == 0:
            committer.add_skip(iog)
            continue
        n_n1 = iog % 3     # 0 to 2 HOGs at N1 in this OG
        cached_hogs = [("N0", [og, "n0", "A_%d" % iog, "", "C_%d" % iog]),
                       ("N0.ids", [og, "n0", "0_%d" % iog, "", "2_%d" % iog])]
        cached_hogs += [("N1", [og, "n%d" % k, "A_%d_%d" % (iog, k), "B_%d_%d" % (iog, k), ""])
                        for k in range(n_n1)]
        committer.add_result(iog, cached_hogs)
    committer.assert_finished()
    committer.close()
    hog_writer.close_files()

    def rows(name):
        with open(tmp_path / (name + ".tsv")) as f:
            return [line.rstrip("\n").split("\t") for line in f]

    n0, n0_ids, n1 = rows("N0"), rows("N0.ids"), rows("N1")
    with_hogs = [i for i in iogs if i % 5 != 0]
    assert [r[0] for r in n0] == ["N0.HOG%07d" % i for i in range(len(with_hogs))]
    assert [r[1] for r in n0] == ["OG%07d" % i for i in with_hogs]
    assert [r[0] for r in n0_ids] == [r[0] for r in n0]
    assert [r[1] for r in n0_ids] == [r[1] for r in n0]
    n1_ogs = [i for i in with_hogs for _ in range(i % 3)]
    assert [r[0] for r in n1] == ["N1.HOG%07d" % i for i in range(len(n1_ogs))]
    assert [r[1] for r in n1] == ["OG%07d" % i for i in n1_ogs]


def test_text_buffers_many_files_bounded(tmp_path):
    # with many files (e.g. a HOG file per species tree node) the text is appended in
    # pieces, in the order written, with at most budget characters held
    from orthofinder.orthology.ortholog_writers import TextBuffers

    paths = [str(tmp_path / ("N%d.tsv" % i)) for i in range(300)]
    for p in paths:
        with open(p, "w") as f:
            f.write("HOG\n")
    tb = TextBuffers(budget=2000, member_bytes=500)
    expected = {p: "HOG\n" for p in paths}
    rng = random.Random(2)
    for k in range(5000):
        p = paths[rng.randrange(len(paths))]
        line = "row%d\n" % k
        tb.write(p, line)
        expected[p] += line
        assert tb.total <= 2000
    tb.close()
    for p in paths:
        assert open(p).read() == expected[p]


def test_buffers_over_budget_write_largest_down_to_half(tmp_path):
    from orthofinder.orthology.ortholog_writers import TextBuffers
    buffers = TextBuffers(budget=100, member_bytes=1000)
    small, large = str(tmp_path / "small.txt"), str(tmp_path / "large.txt")
    for fn in (small, large):
        open(fn, "w").close()
    buffers.write(small, "s" * 10)
    buffers.write(large, "l" * 95)     # 105 > 100: the largest is written, leaving 10 <= 50
    assert buffers.total == 10 and open(large).read() == "l" * 95 and open(small).read() == ""
    buffers.close()
    assert open(small).read() == "s" * 10


def test_resolved_tree_not_written_stops_the_run(tmp_path, monkeypatch):
    # (the tree's orthologs and HOGs must not be dropped with only a warning)
    import pytest
    from orthofinder.orthology import tree_analyser
    from orthofinder.utils.parallel_task_manager import WorkerError
    ta = object.__new__(tree_analyser.TreeAnalyser)
    ta.recon_trees_renamed_dir = str(tmp_path) + "/"
    ta.spec_seq_dict = {}
    ta.write_hog_tree = False

    def fails(*args, **kwargs):
        raise OSError(28, "No space left on device")
    monkeypatch.setattr(tree_analyser.phylo_taxa, "rename_tree_taxa", fails)
    with pytest.raises(WorkerError, match="OG0000005"):
        ta.write_resolved_tree(None, 5)
