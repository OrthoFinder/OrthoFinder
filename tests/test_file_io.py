import csv
import gzip
import io
import random

import numpy as np
import pytest

from orthofinder.utils import file_io

ALPHABET = ["a", "Z", "_", ".", ",", " ", ", ", "\t", '"', "\r", "\n", "\r\n", "é", "\x00", "#", "\\"]


def random_value(rng):
    r = rng.random()
    if r < 0.6:
        return "".join(rng.choice(ALPHABET[:7] if rng.random() < 0.8 else ALPHABET)
                       for _ in range(rng.randrange(0, 8)))
    return rng.choice([None, 0, -3, 1.5, 0.1, 1e-300, True, np.int64(7), np.float64(2.25),
                       np.float32(0.5), float("nan"), ""])


def random_rows(seed, n=3000):
    rng = random.Random(seed)
    return [[random_value(rng) for _ in range(rng.randrange(0, 6))] for _ in range(n)]


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("terminator", ["\r\n", "\n"])
def test_writer_matches_csv(seed, terminator):
    rows = random_rows(seed)
    expected, actual = io.StringIO(), io.StringIO()
    w = csv.writer(expected, delimiter="\t", lineterminator=terminator)
    for row in rows:
        w.writerow(row)
    file_io.writer(actual, lineterminator=terminator).writerows(rows)
    assert actual.getvalue() == expected.getvalue()


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("terminator", ["\r\n", "\n"])
@pytest.mark.parametrize("newline", [None, ""])
@pytest.mark.parametrize("quotes", [True, False])
def test_reader_matches_csv(tmp_path, seed, terminator, newline, quotes):
    rows = random_rows(seed)
    if not quotes:
        rows = [[v for v in row if not (isinstance(v, str) and set(v) & set('\t"\r\n'))] for row in rows]
    fn = tmp_path / "t.tsv"
    with open(fn, "w", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator=terminator).writerows(rows)
    with open(fn, newline=newline) as f:
        expected = list(csv.reader(f, delimiter="\t"))
    with open(fn, newline=newline) as f:
        actual = list(file_io.reader(f))
    assert actual == expected


def _old_getrow(row):
    # util.getrow, which file_io.unquoted_line replaced
    cleaned = [str(x).rstrip("\r\n") for x in row]
    return "\t".join(cleaned) + "\n"


@pytest.mark.parametrize("seed", range(5))
def test_unquoted_line_matches_old_getrow(seed):
    rows = random_rows(seed)
    assert [file_io.unquoted_line(r) for r in rows] == [_old_getrow(r) for r in rows]
    assert [file_io.unquoted_line(iter(r)) for r in rows] == [_old_getrow(r) for r in rows]
    out = io.StringIO()
    for r in rows:
        file_io.write_unquoted(out, r)
    assert out.getvalue() == "".join(_old_getrow(r) for r in rows)


def test_unquoted_reader_keeps_quotes(tmp_path):
    # Rows as written to HOG files and Duplications.tsv: tab-joined, never quoted.
    rows = [["N0.HOG0000000", "OG0000000", "n0", '"geneA, geneB', "geneC"],
            ["N0.HOG0000001", "OG0000001", "-", 'g"x', ""],
            []]
    fn = tmp_path / "N0.tsv"
    fn.write_text("".join("\t".join(r) + "\n" for r in rows))
    with open(fn) as f:
        assert list(file_io.unquoted_reader(f)) == rows
    with open(fn) as f:   # what a quote-aware reader makes of the same line
        assert next(csv.reader(f, delimiter="\t")) != rows[0]


def test_writer_default_line_ending_is_lf():
    out = io.StringIO()
    w = file_io.writer(out)
    w.writerow(["Orthogroup", "Sp0"])
    w.writerow(["OG0000000", 'gene"x'])     # quoted row, written by csv
    assert out.getvalue() == 'Orthogroup\tSp0\nOG0000000\t"gene""x"\n'


def test_reader_gzip_and_next(tmp_path):
    fn = tmp_path / "t.tsv.gz"
    with gzip.open(fn, "wt") as f:
        w = file_io.writer(f)
        w.writerow(["Orthogroup", "Sp0", "Sp1"])
        w.writerow(["OG0000000", "g1, g2", ""])
        w.writerow(["OG0000001", "", 'gene"x'])
    with gzip.open(fn, "rt") as f:
        r = file_io.reader(f)
        assert next(r) == ["Orthogroup", "Sp0", "Sp1"]
        assert list(r) == [["OG0000000", "g1, g2", ""], ["OG0000001", "", 'gene"x']]


def test_file_info_reports_line_endings(tmp_path):
    fn = tmp_path / "t.tsv"
    fn.write_bytes(b"a\tb\r\nc\td\n")
    info = file_io.file_info(str(fn))
    assert info["lines"] == 2 and info["fields per line"] == "2"
    assert info["line endings"].startswith("mixed")


@pytest.mark.parametrize("ext", ["", ".gz"])
def test_convert_round_trip(tmp_path, ext):
    rows = [["Orthogroup", "Sp0", "Sp1"], ["OG0000000", '"gene', 'g"x'], ["OG0000001", "", "y"]]
    unquoted = tmp_path / ("in.tsv" + ext)
    with file_io.open_text(str(unquoted), "wt") as f:
        for r in rows:
            file_io.write_unquoted(f, r)
    quoted = tmp_path / ("q.tsv" + ext)
    back = tmp_path / ("back.tsv" + ext)
    assert file_io.convert(str(unquoted), str(quoted), "unquoted", "quoted") == 3
    with file_io.open_text(str(quoted)) as f:
        assert list(file_io.reader(f)) == rows
    file_io.convert(str(quoted), str(back), "quoted", "unquoted")
    with file_io.open_text(str(back)) as f:
        assert list(file_io.unquoted_reader(f)) == rows


def test_convert_to_unquoted_rejects_tab_in_field(tmp_path):
    fn = tmp_path / "t.tsv"
    fn.write_text('h\tx\n1\t"a\tb"\n')
    with pytest.raises(ValueError):
        file_io.convert(str(fn), str(tmp_path / "out.tsv"), "quoted", "unquoted")


def test_standalone_copy_imports_without_orthofinder(tmp_path):
    import shutil, subprocess, sys
    shutil.copy(file_io.__file__, tmp_path / "file_io.py")
    code = "import file_io, sys; assert not any(m.startswith('orthofinder') for m in sys.modules)"
    subprocess.run([sys.executable, "-c", code], cwd=tmp_path, check=True)
    out = subprocess.run([sys.executable, "file_io.py", "--help"], cwd=tmp_path, check=True,
                         capture_output=True, text=True).stdout
    assert "info" in out and "convert" in out


_HITS = ("3_17\t5_2\t60.8\t120\t40\t1\t1\t120\t3\t122\t5.99e-111\t315\n"
         "3_18\t5_9\t99.0\t50\t0\t0\t2\t51\t1\t50\t1e-05\t88.5\n")
_STD = file_io.STANDARD_HIT_FORMAT


def test_read_hit_columns_by_name(tmp_path):
    pytest.importorskip("pandas")
    fn = tmp_path / "Blast3_5.txt"
    fn.write_text(_HITS)
    names = ("query_species", "query_seq", "subject_seq", "pident", "qstart", "evalue", "bitscore")
    cols = file_io.read_hit_columns(str(fn), names, fmt=_STD)
    assert list(cols) == list(names)
    assert cols["query_species"].tolist() == [3, 3] and cols["query_seq"].tolist() == [17, 18]
    assert cols["subject_seq"].tolist() == [2, 9] and cols["qstart"].tolist() == [1, 2]
    assert cols["pident"].tolist() == [60.8, 99.0] and cols["evalue"].tolist() == [5.99e-111, 1e-05]
    assert cols["bitscore"].tolist() == [315.0, 88.5]


def test_read_hit_columns_custom_id_parts(tmp_path):
    pytest.importorskip("pandas")
    fn = tmp_path / "Blast4_-1.txt"
    fn.write_text("4_7\t0000012_1_33\t50\t10\t5\t0\t1\t10\t1\t10\t2e-30\t120\n")
    parts = {"qseqid": ("query_species", "query_seq"), "sseqid": ("og", "subject_species", "subject_seq")}
    cols = file_io.read_hit_columns(str(fn), ("og", "subject_species", "evalue"), fmt=_STD,
                                    id_parts=parts, dtypes={"og": str})
    assert cols["og"].tolist() == ["0000012"] and cols["subject_species"].tolist() == [1]
    assert cols["evalue"].tolist() == [2e-30]


def test_read_hit_columns_errors(tmp_path):
    pytest.importorskip("pandas")
    fn = tmp_path / "Blast3_5.txt"
    fn.write_text(_HITS)
    with pytest.raises(KeyError):           # unknown column: a programming error
        file_io.read_hit_columns(str(fn), ("query_seq", "no_such_column"), fmt=_STD)
    fn.write_text("3_17\t5_2\t60.8\n")      # not the expected layout: callers fall back
    with pytest.raises(ValueError):
        file_io.read_hit_columns(str(fn), ("query_seq", "bitscore"), fmt=_STD)
    text = file_io.HitFormat(("qseqid", "sseqid", "stitle", "bitscore"))   # text column may hold "_"
    fn.write_text("3_17\t5_2\tsome_title\t315\n")
    with pytest.raises(ValueError):
        file_io.read_hit_columns(str(fn), ("query_seq", "bitscore"), fmt=text)


@pytest.mark.parametrize("cmd, fields, header", [
    ("diamond blastp -d DB -q IN -o OUT --more-sensitive --compress 1", file_io.BLAST6_FIELDS, False),
    ("blastp -outfmt 6 -query IN -db DB | gzip > OUT.gz", file_io.BLAST6_FIELDS, False),
    ("mmseqs search A B C tmp ; mmseqs convertalis A B C OUT", file_io.MMSEQS_DEFAULT_FIELDS, False),
    ('blastp -outfmt "6 qseqid sseqid evalue bitscore" -query IN', ("qseqid", "sseqid", "evalue", "bitscore"), False),
    ("blastp -outfmt '7 std qlen' -query IN", file_io.BLAST6_FIELDS + ("qlen",), False),
    ("diamond blastp -d DB -q IN -o OUT --outfmt 6 sseqid qseqid bitscore -p 1", ("sseqid", "qseqid", "bitscore"), False),
    ("diamond blastp -d DB -q IN -o OUT -f 6 qseqid sseqid bitscore", ("qseqid", "sseqid", "bitscore"), False),
    ('mmseqs convertalis A B C OUT --format-output "query,target,evalue,bits" --format-mode 4',
     ("qseqid", "sseqid", "evalue", "bitscore"), True),
])
def test_hit_format_from_command(cmd, fields, header):
    fmt = file_io.search_hit_format(cmd)
    assert fmt.fields == tuple(fields) and fmt.header == header and not fmt.assumed


@pytest.mark.parametrize("cmd", [
    "blastp -query IN -db DB -out OUT",                              # -outfmt 0: pairwise text
    "blastp -outfmt 10 -query IN",                                   # CSV
    "diamond blastp -d DB -q IN -o OUT --outfmt 5",                  # XML
    "diamond blastp -d DB -q IN -o OUT --outfmt 6 qseqid sseqid evalue",   # no bit score
    "mmseqs convertalis A B C OUT --format-mode 1",                  # SAM
])
def test_hit_format_unreadable_or_incomplete(cmd):
    with pytest.raises(ValueError):
        file_io.search_hit_format(cmd)


def test_hit_format_declared_and_assumed():
    fmt = file_io.search_hit_format("mytool -i IN -o OUT", declared_fields="sseqid qseqid bitscore")
    assert fmt.fields == ("sseqid", "qseqid", "bitscore")
    fmt = file_io.search_hit_format("mytool -i IN -o OUT")
    assert fmt.fields == file_io.BLAST6_FIELDS and fmt.assumed


def test_hit_format_from_blast_commands(tmp_path):
    """The layout is read from blast_commands.txt next to the results, which every search writes."""
    def commands_dir(name, lines):
        d = tmp_path / name
        d.mkdir()
        (d / "blast_commands.txt").write_text("".join(l + "\n" for l in lines))
        return d
    d = commands_dir("blast", ["blastp -outfmt '6 sseqid qseqid bitscore' -query a.fa -db x | gzip > Blast0_0.txt.gz"])
    assert file_io.hit_format(str(d / "Blast0_0.txt.gz")).fields == ("sseqid", "qseqid", "bitscore")
    d = commands_dir("mmseqs", ['mmseqs convertalis A B C Blast0_0.txt --format-output "query,target,bits" --format-mode 4'])
    fmt = file_io.hit_format(str(d / "Blast0_0.txt"))
    assert fmt.fields == ("qseqid", "sseqid", "bitscore") and fmt.header
    # columns declared in the config ("output_fields") come first, as a comment
    declared = file_io.search_hit_format("mytool -i IN -o OUT", declared_fields="sseqid qseqid bitscore")
    lines = file_io.commands_file_lines(["mytool -i a -o Blast0_0.txt"], declared)
    assert lines[0] == "# output_fields: sseqid qseqid bitscore"
    d = commands_dir("custom", lines)
    assert file_io.hit_format(str(d / "Blast0_0.txt")).fields == ("sseqid", "qseqid", "bitscore")
    assert file_io.commands_file_lines(["blastp -outfmt 6"], file_io.STANDARD_HIT_FORMAT) == ["blastp -outfmt 6"]
    # results made outside OrthoFinder: the standard columns
    assert file_io.hit_format(str(tmp_path / "none" / "Blast0_0.txt")).assumed


@pytest.mark.parametrize("use_pandas", [True, False])
def test_blast_scores_with_custom_layout(tmp_path, monkeypatch, use_pandas):
    """Same matrix from a standard file and from one with other columns, order, comments and a header."""
    pytest.importorskip("pandas")
    import types
    from orthofinder.utils import blast_file_processor as bfp
    std_dir, custom_dir = tmp_path / "std", tmp_path / "custom"
    std_dir.mkdir(), custom_dir.mkdir()
    (std_dir / "Blast3_5.txt").write_text(_HITS)
    (std_dir / "blast_commands.txt").write_text("diamond blastp -d DB -q a.fa -o Blast3_5.txt\n")
    # same hits as "evalue bitscore sseqid qseqid", with a comment and a header line
    lines = ["# a comment\n", "evalue\tbits\ttarget\tquery\n"]
    for l in _HITS.splitlines():
        r = l.split("\t")
        lines.append("\t".join([r[10], r[11], r[1], r[0]]) + "\n")
    (custom_dir / "Blast3_5.txt").write_text("".join(lines))
    (custom_dir / "blast_commands.txt").write_text(
        'mmseqs convertalis A B C Blast3_5.txt --format-output "evalue,bits,target,query" --format-mode 4\n')
    monkeypatch.setattr(file_io, "use_pandas", lambda p: use_pandas)
    info = types.SimpleNamespace(nSeqsPerSpecies={3: 20, 5: 20})
    a = bfp.GetBLAST6Scores(info, [str(std_dir)], 3, 5, fmt="csr")
    b = bfp.GetBLAST6Scores(info, [str(custom_dir)], 3, 5, fmt="csr")
    assert a.nnz == 2 and (a != b).nnz == 0
    assert a[17, 2] == 315.0 and a[18, 9] == 88.5


def test_search_checkpoint(tmp_path):
    hits = tmp_path / "Blast0_1.txt.gz"
    hits.write_bytes(b"x" * 123)
    line = file_io.search_completed_message(str(hits))
    assert line == "Search completed: Blast0_1.txt.gz (123 bytes)"

    def status(*lines):
        (tmp_path / "checkpoint.txt").write_text("".join("2026-10-04 10:00:00 : %s\n" % l for l in lines))
        return file_io.search_checkpoint(str(tmp_path))

    started, completed = "[Sequence search] Started", "[Sequence search] Completed"
    # records, newest wins
    assert status("Search completed: Blast0_1.txt.gz (100 bytes)", line,
                  "Search completed: Blast1_0.txt (55 bytes)")[1] == {"Blast0_1.txt.gz": 123, "Blast1_0.txt": 55}
    assert status(started, line, completed)[0] is False
    assert status(started, line)[0] is True                                  # killed: no more lines
    assert status(started, "[Sequence search] Interrupted by user")[0] is True
    assert status(started, "[ERROR] [Sequence search] ERROR: Workflow exited with status 1")[0] is True
    assert status(started, "[Sequence search] Workflow stopped normally before the full analysis finished")[0] is False  # -op
    assert status(started, completed, started)[0] is True                   # a later restart interrupted
    assert status(started, started, completed)[0] is False                  # ... then finished
    assert status("Starting OrthoFinder")[0] is False                        # no search logged
    assert file_io.search_checkpoint(str(tmp_path / "missing")) == (False, {})
