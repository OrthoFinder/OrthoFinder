# -*- coding: utf-8 -*-
"""
OrthoFinder's file readers and writers, kept in one place.

Add new readers/writers here rather than in the modules that use them, so
each file format has one implementation and its reader and writer match.

Text files are UTF-8, with "\n" line endings on every OS (open_text), so that
results are the same bytes on Linux, macOS and Windows and are read the same
everywhere. Readers accept "\n", "\r\n" and "\r".

Tab-separated files  reader, writer, unquoted_line, write_unquoted
    Orthogroups.tsv and the other orthogroup tables, the HOG files, the
    statistics files, the ortholog and xenolog files, Duplications.tsv. Fields
    are joined by tabs as they are, never quoted: the names in them cannot hold
    a tab, line break or quote (gene_ids.clean_name), so the files read the
    same with csv readers, spreadsheets, pandas, R or str.split("\t"). A line
    break at the end of a field is dropped.

Search hits  hit_format, read_hit_columns, hit_lines, search_hit_format
    Search results (Blast*_*.txt[.gz]), written by DIAMOND, BLAST+, MMseqs2
    or a user-configured program. Their columns depend on the program and its
    output options, so the layout is worked out from the search command
    (kept in blast_commands.txt next to the results) and columns are read
    by name. Large files are read with pandas (optional), small ones
    line by line by the caller (blast_file_processor, orthogroups.assign).

Standalone use
    The module needs only the Python standard library (pandas and numpy are
    optional, for read_hit_columns), so it can be copied into other projects
    and imported as `import file_io`. It also runs as a command:

        python file_io.py info FILE...
        python -m orthofinder.utils.file_io info Orthogroups.tsv

    Files ending in .gz are read gzipped.
"""
import argparse
import gzip
import io
import os
import re


# ------------------------------------------------------------ text files ----

ENCODING = "utf-8"
# The gzip compression level of the files written: 6 (zlib's default) compresses
# ~8x faster than 9 (gzip's default) and the files are <1% larger (measured on
# ortholog files: 55 vs 6.4 MB/s, 16.8% of the text either way)
GZIP_LEVEL = 6


def open_text(path, mode="rt", errors=None):
    """
    A text file (gzipped if path ends in .gz): UTF-8, its lines written ending in
    "\n" on every OS, and read with any line ending (as "\n").
    """
    if "b" in mode:
        raise ValueError("open_text opens text files: mode %r" % mode)
    if "t" not in mode:
        mode += "t"
    newline = None if "r" in mode else "\n"
    if path.endswith(".gz"):
        level = {} if "r" in mode else {"compresslevel": GZIP_LEVEL}
        return gzip.open(path, mode, encoding=ENCODING, errors=errors, newline=newline, **level)
    return open(path, mode, encoding=ENCODING, errors=errors, newline=newline)


class open_text_any(object):
    """
    The lines of a text file that an older version or a user wrote, on any OS: older
    versions wrote SpeciesIDs.txt, SequenceIDs.txt, Log.txt and the species trees in
    the system's encoding (e.g. cp1252 on Windows), this one in UTF-8, and a user's
    file (e.g. a species tree, -s) may be either, with a byte order mark. Each line
    is read as UTF-8, or else as cp1252, so that e.g. a restart from older results
    (-b) reads them. A file object for reading (for line in f, next(f), read(),
    with ... as f).
    """

    def __init__(self, path):
        self._f = open(path, encoding="utf-8-sig", errors="surrogateescape")

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self._f)
        if not line.isascii():      # (nearly every line is: nothing else done)
            try:
                line.encode(ENCODING)
            except UnicodeEncodeError:
                line = line.encode(ENCODING, "surrogateescape").decode("cp1252", errors="replace")
        return line

    def read(self):
        return "".join(self)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._f.close()

    def close(self):
        self._f.close()


# ---------------------------------------------------------- tab-separated ----

def reader(infile):
    """The rows of a tab-separated file (lists of str; [] for an empty line)."""
    for line in infile:
        line = line.rstrip("\r\n")
        yield line.split("\t") if line else []


class writer(object):
    """Writes the rows of a tab-separated file (unquoted_line), like csv.writer."""

    def __init__(self, outfile):
        self._write = outfile.write

    def writerow(self, row):
        return self._write(unquoted_line(row))

    def writerows(self, rows):
        for row in rows:
            self._write(unquoted_line(row))


def unquoted_line(row):
    """
    One line of a tab-separated file: the fields as str (None as "None"), each
    without trailing line breaks, joined by tabs, ending in "\n".
    """
    if type(row) is not list and type(row) is not tuple:
        row = list(row)
    try:
        line = "\t".join(row)                     # rows are usually all str
    except TypeError:
        row = [x if type(x) is str else str(x) for x in row]
        line = "\t".join(row)
    if "\n" in line or "\r" in line:
        return "\t".join([x.rstrip("\r\n") for x in row]) + "\n"
    return line + "\n"


def write_unquoted(fh, row):
    """Write one line of a tab-separated file (see unquoted_line)."""
    fh.write(unquoted_line(row))


# ----------------------------------------------------- search-program hits ----
#
# Search results (Blast<i>_<j>.txt[.gz]) are tab-separated, one hit per line,
# but which columns they have depends on the search program and on the output
# options in its command, which users can change in the config file. The
# layout is worked out from the command (hit_format_from_command) before the
# search runs, and again from blast_commands.txt when the results are read
# (hit_format), so a column is always read by its name.

# The standard columns of BLAST tabular output ("std"): DIAMOND's default,
# blastp/blastn -outfmt 6.
BLAST6_FIELDS = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# MMseqs2's default (convertalis / easy-search --format-output), in the names
# used here: column 3 is "fident", a fraction (0.304) rather than a percentage.
MMSEQS_DEFAULT_FIELDS = ("qseqid", "sseqid", "fident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# MMseqs2 column names -> the BLAST names used here (others keep their name).
_MMSEQS_NAMES = {"query": "qseqid", "target": "sseqid", "bits": "bitscore",
                 "alnlen": "length", "tstart": "sstart", "tend": "send", "tlen": "slen"}

# Columns OrthoFinder reads: all searches, and the orthogroup-profile searches.
REQUIRED_HIT_FIELDS = ("qseqid", "sseqid", "bitscore")
REQUIRED_PROFILE_HIT_FIELDS = ("qseqid", "sseqid", "evalue")

# OrthoFinder's sequence IDs ("3_17": species 3, sequence 17) are split at "_"
# into these named parts by read_hit_columns.
SEQ_ID_PARTS = {"qseqid": ("query_species", "query_seq"),
                "sseqid": ("subject_species", "subject_seq")}

# Numeric columns, with the type they are read as. Other columns are text,
# which may contain "_" (stitle, ...), so read_hit_columns does not read them.
_HIT_DTYPES = {
    "query_species": "int64", "query_seq": "int64",
    "subject_species": "int64", "subject_seq": "int64",
    "length": "int64", "mismatch": "int64", "gapopen": "int64", "gaps": "int64",
    "nident": "int64", "positive": "int64", "score": "int64", "raw": "int64",
    "qstart": "int64", "qend": "int64", "sstart": "int64", "send": "int64",
    "qlen": "int64", "slen": "int64",
    "pident": "float64", "fident": "float64", "ppos": "float64",
    "evalue": "float64", "bitscore": "float64",
    "qcovs": "float64", "qcovhsp": "float64", "qcovus": "float64",
    "qcov": "float64", "tcov": "float64",
}

_BLAST_PROGRAMS = {"blastp", "blastn", "blastx", "tblastn", "tblastx", "psiblast",
                   "deltablast", "rpsblast", "rpstblastn"}


class HitFormat(object):
    """
    The layout of a search-results file: its columns (names as in
    BLAST6_FIELDS), whether its first line is a header (MMseqs2
    --format-mode 4), and whether it was declared in the config
    ("output_fields") or assumed (a program whose options OrthoFinder cannot
    read) rather than read from the search command.
    Lines starting with "#" (BLAST -outfmt 7) are always skipped.
    """

    def __init__(self, fields, header=False, assumed=False, declared=False):
        self.fields = tuple(fields)
        self.header = header
        self.assumed = assumed
        self.declared = declared

    def index(self, field):
        """Position of a column (ValueError if the file does not have it)."""
        return self.fields.index(field)

    def missing(self, required):
        return [f for f in required if f not in self.fields]

    def __eq__(self, other):
        return (self.fields, self.header) == (other.fields, other.header)

    def __repr__(self):
        return "HitFormat(%r, header=%r, assumed=%r)" % (self.fields, self.header, self.assumed)


STANDARD_HIT_FORMAT = HitFormat(BLAST6_FIELDS)


def _command_segments(cmd):
    """The simple commands of a shell command line (split at | ; && ||), as token lists."""
    import shlex
    lexer = shlex.shlex(cmd, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    segment = []
    for token in lexer:
        if token in ("|", ";", "&&", "||", "&"):
            if segment:
                yield segment
            segment = []
        else:
            segment.append(token)
    if segment:
        yield segment


def _option_value(tokens, names):
    """The value of the last of the options `names` in a token list (also --opt=value), or None."""
    value = None
    for i, t in enumerate(tokens):
        for name in names:
            if t == name and i + 1 < len(tokens):
                value = (i + 1, tokens[i + 1])
            elif t.startswith(name + "="):
                value = (i, t[len(name) + 1:])
    return value


_FIELD_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")   # e.g. qseqid, full_sseq


def _blast_fields(words):
    """Expand BLAST/DIAMOND field words ("std" -> the 12 standard columns)."""
    fields = []
    for w in words:
        fields.extend(BLAST6_FIELDS if w == "std" else (w,))
    return tuple(fields) if fields else BLAST6_FIELDS


def hit_format_from_command(cmd):
    """
    The layout of the results a search command writes, from its output
    options: BLAST+ (-outfmt 6/7 [fields]), DIAMOND (--outfmt/-f 6 [fields])
    and MMseqs2 (convertalis/easy-search --format-output, --format-mode 0/4).
    Returns None if the command runs none of these programs.
    Raises ValueError if the output is not tab-separated text (e.g. BLAST
    -outfmt 0 or 10, DIAMOND --outfmt 5), which OrthoFinder cannot read.
    """
    fmt = None
    for tokens in _command_segments(cmd):
        program = os.path.basename(tokens[0])
        if program in _BLAST_PROGRAMS:
            found = _option_value(tokens, ("-outfmt",))
            words = found[1].split() if found else ["0"]
            if words[0] not in ("6", "7"):
                raise ValueError(
                    "%s -outfmt %s does not write tab-separated hits: use -outfmt 6 (or 7), "
                    "e.g. -outfmt 6 or -outfmt \"6 qseqid sseqid evalue bitscore\"" % (program, words[0]))
            fmt = HitFormat(_blast_fields(words[1:]))
        elif program == "diamond" and len(tokens) > 1 and tokens[1] in ("blastp", "blastx"):
            found = _option_value(tokens, ("--outfmt", "-f"))
            if found is None:
                fmt = STANDARD_HIT_FORMAT
                continue
            i, first = found
            words = first.split()
            if words[0] != "6":
                raise ValueError(
                    "diamond --outfmt %s does not write tab-separated hits: use --outfmt 6, "
                    "e.g. --outfmt 6 qseqid sseqid evalue bitscore" % words[0])
            words = words[1:]
            if len(first.split()) == 1:     # fields as separate arguments
                for t in tokens[i + 1:]:
                    # up to the next option, or anything that is not a field
                    # name (e.g. a redirection: "2> err.log")
                    if not _FIELD_NAME_RE.match(t):
                        break
                    words.append(t)
            fmt = HitFormat(_blast_fields(words))
        elif program == "mmseqs" and len(tokens) > 1 and tokens[1] in (
                "convertalis", "easy-search", "easy-linsearch"):
            mode = _option_value(tokens, ("--format-mode",))
            mode = mode[1] if mode else "0"
            if mode not in ("0", "4"):
                raise ValueError(
                    "mmseqs --format-mode %s does not write tab-separated hits: use "
                    "--format-mode 0 (the default) or 4" % mode)
            out = _option_value(tokens, ("--format-output",))
            fields = (tuple(_MMSEQS_NAMES.get(f, f) for f in out[1].split(","))
                      if out else MMSEQS_DEFAULT_FIELDS)
            fmt = HitFormat(fields, header=(mode == "4"))
    return fmt


def search_hit_format(cmd, declared_fields=None, required=REQUIRED_HIT_FIELDS):
    """
    The layout of a search command's results: from "output_fields" in the
    config if declared, else from the command's options, else the standard
    12 columns (assumed). Raises ValueError, with an explanation, if the
    output cannot be read or lacks a column OrthoFinder needs.
    """
    if declared_fields:
        fmt = HitFormat(declared_fields.replace(",", " ").split(), declared=True)
    else:
        fmt = hit_format_from_command(cmd) or HitFormat(BLAST6_FIELDS, assumed=True)
    missing = fmt.missing(required)
    if missing:
        raise ValueError(
            "the search results would not have the column(s) %s, which OrthoFinder needs "
            "(columns: %s). Command: %s" % (", ".join(missing), " ".join(fmt.fields), cmd))
    return fmt


# blast_commands.txt (written next to the results of every search, with the
# commands) starts with this comment when the columns were declared in the
# config ("output_fields"), as they cannot be read from the command then.
OUTPUT_FIELDS_COMMENT = "# output_fields:"


def commands_file_lines(commands, fmt):
    """The lines of blast_commands.txt: the commands, after the declared columns if any."""
    lines = []
    if fmt is not None and fmt.declared:
        lines.append("%s %s" % (OUTPUT_FIELDS_COMMENT, " ".join(fmt.fields)))
    return lines + list(commands)


_hit_format_cache = {}


def clear_caches():
    """Forget the layouts read so far (a new run may reuse a directory)."""
    _hit_format_cache.clear()


def hit_format(path):
    """
    The layout of a search-results file (Blast<i>_<j>.txt), from the search
    command in blast_commands.txt in the same directory, which every search
    writes next to its results (also for a run restarted from them); the
    standard 12 columns if there is no such file (e.g. results made outside
    OrthoFinder). The layout of the orthogroup-profile searches is passed to
    their reader instead (orthogroups.assign), as they are read in the same run.
    """
    directory = os.path.dirname(os.path.abspath(path))
    if directory not in _hit_format_cache:
        fmt = None
        commands_fn = os.path.join(directory, "blast_commands.txt")
        if os.path.exists(commands_fn):
            with open(commands_fn, encoding="utf-8") as infile:
                for line in infile:
                    if line.startswith(OUTPUT_FIELDS_COMMENT):
                        fmt = HitFormat(line[len(OUTPUT_FIELDS_COMMENT):].split(), declared=True)
                        break
                    if line.strip() and not line.startswith("#"):
                        try:
                            fmt = hit_format_from_command(line.strip())
                        except ValueError:
                            fmt = None
                        break
        _hit_format_cache[directory] = fmt or HitFormat(BLAST6_FIELDS, assumed=True)
    return _hit_format_cache[directory]


def hit_lines(infile, fmt):
    """The hit lines of an open results file: without "#" comment lines and the header line."""
    skip_header = fmt.header
    for line in infile:
        if not line.strip() or line.startswith("#"):
            continue
        if skip_header:
            skip_header = False
            continue
        yield line


# Hit files are parsed in blocks of this many bytes (bounded memory).
HIT_READ_BLOCK_BYTES = 64 * 1024 * 1024

# Smallest hit files (on disk) parsed with pandas. Each pandas call costs
# ~0.3 ms however small the file, so below ~5,000 lines the line-by-line
# reader is faster; above it pandas is, up to ~1.5x for 1M lines.
PANDAS_MIN_BYTES = 256 * 1024
PANDAS_MIN_BYTES_GZ = 64 * 1024

_pd = False   # pandas module, None if not installed; imported on first use


def _pandas():
    global _pd
    if _pd is False:
        try:
            import pandas
            _pd = pandas
        except ImportError:   # optional: without it the line-by-line readers are used
            _pd = None
    return _pd


def use_pandas(path):
    """Whether pandas is installed and the hit file is large enough to gain from it."""
    if _pandas() is None:
        return False
    limit = PANDAS_MIN_BYTES_GZ if path.endswith(".gz") else PANDAS_MIN_BYTES
    return os.path.getsize(path) >= limit


def _hit_file_blocks(path, fmt):
    """Complete hit lines of a (gzipped) text file, in blocks of bytes, without comment/header lines."""
    import re
    comment = re.compile(rb"(?m)^#[^\n]*\n")
    skip_header = fmt.header
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as infile:
        rest = b""
        while True:
            block = infile.read(HIT_READ_BLOCK_BYTES)
            if not block:
                data = rest if rest.endswith(b"\n") or not rest.strip() else rest + b"\n"
            else:
                data = rest + block
                cut = data.rfind(b"\n") + 1
                if cut == 0:
                    rest = data
                    continue
                data, rest = data[:cut], data[cut:]
            if b"#" in data:
                data = comment.sub(b"", data)
            if skip_header and data.strip():
                data = data.lstrip(b"\n")
                data = data[data.find(b"\n") + 1:]
                skip_header = False
            if data.strip():
                yield data
            if not block:
                return


def hit_column_names(fields=BLAST6_FIELDS, id_parts=SEQ_ID_PARTS):
    """The columns of a hit file once the ID fields are split into their parts."""
    names = []
    for field in fields:
        names.extend(id_parts.get(field, (field,)))
    return names


def read_hit_columns(path, columns, fmt=None, id_parts=SEQ_ID_PARTS, dtypes=None):
    """
    Read named columns of a search-results file with pandas' C parser.

    columns: names from hit_column_names(fmt.fields, id_parts), e.g.
        ("query_seq", "subject_seq", "bitscore"); any numeric column of the
        file can be asked for, e.g. "evalue", "pident" or "qstart".
    fmt: the file's layout (default: hit_format(path), as recorded when the
        search was run).
    id_parts: the fields holding "_"-joined IDs, and the names of their parts
        (default: SEQ_ID_PARTS, e.g. "qseqid" "3_17" -> query_species 3,
        query_seq 17).
    dtypes: {column: dtype} to override _HIT_DTYPES; str reads a column as
        text (returned as an object array).

    "_" is treated as a field separator while parsing, which is how the ID
    parts become separate integer columns, so every other column of the file
    must be numeric (a text column such as stitle may contain "_"). Much faster
    than a Python loop over the lines: no Python object is created per field,
    and numbers are parsed into arrays directly. Floats are parsed with
    float_precision="round_trip", which is bit-identical to float().

    Returns {column: numpy array}, or None if pandas is not installed. Raises
    KeyError for a column the layout does not have (a programming error), and
    ValueError if the file has a text column, a line does not fit the layout
    or a value is missing; callers then fall back to their line-by-line reader,
    which reports any problem.
    """
    pd = _pandas()
    if pd is None:
        return None
    import numpy as np
    if fmt is None:
        fmt = hit_format(path)
    names = hit_column_names(fmt.fields, id_parts)
    for c in columns:
        if c not in names:
            raise KeyError("no column %r in %s" % (c, names))
    dtypes = dtypes or {}
    text = [n for n in names if n not in _HIT_DTYPES and n not in dtypes]
    if text:
        raise ValueError("text column(s) %s: read line by line" % ", ".join(text))
    n_fields = len(names)
    wanted = [dtypes.get(c, _HIT_DTYPES.get(c)) for c in columns]
    read_dtypes = [str if d in (str, object) else np.dtype(d) for d in wanted]
    out_dtypes = [object if d in (str, object) else np.dtype(d) for d in wanted]
    columns = [names.index(c) for c in columns]
    parts = [[] for _ in columns]
    checked = False
    for data in _hit_file_blocks(path, fmt):
        data = data.replace(b"_", b"\t")
        if not checked:
            first = data.lstrip(b"\n")
            first = first[:first.find(b"\n")]
            if first.count(b"\t") + 1 != n_fields:
                raise ValueError("unexpected number of fields")
            checked = True
        try:
            df = pd.read_csv(
                io.BytesIO(data), sep="\t", header=None, engine="c",
                names=list(range(n_fields)), usecols=columns,
                dtype=dict(zip(columns, read_dtypes)),
                float_precision="round_trip", skip_blank_lines=True,
                na_filter=True, keep_default_na=False, na_values=[""],
            )
        except Exception as e:
            raise ValueError("unexpected file layout: %s" % e)
        for i, (col, dtype) in enumerate(zip(columns, out_dtypes)):
            values = df[col]
            if values.isna().any():
                raise ValueError("missing values")
            parts[i].append(values.to_numpy(dtype=dtype))
    return {
        names[col]: np.concatenate(p) if p else np.zeros(0, dtype=dtype)
        for col, p, dtype in zip(columns, parts, out_dtypes)
    }


# A restart (-b) decides from checkpoint.txt (the run log in the working
# directory, appended to by every run there) which search results it can use:
#   - each search that finishes is recorded with the size of its results file
#     ("Search completed: Blast0_1.txt.gz (1234 bytes)");
#   - the search phase is logged as "[Sequence search] Started", then
#     "Completed" when every search finished. If the last one started has no
#     "Completed" (the run was stopped, failed or was killed during the
#     searches), it was interrupted.
# Not interrupted (or no search phase logged at all: -op, searches run outside
# OrthoFinder, older versions): every results file that exists is complete.
# Interrupted: a file is complete only if recorded with its current size.
SEARCH_COMPLETED = "Search completed:"
_SEARCH_COMPLETED_RE = re.compile(r"Search completed: (\S+) \((\d+) bytes\)")
SEARCH_STEP = "[Sequence search] "


def search_completed_message(path):
    """The checkpoint.txt line for a finished search that wrote path."""
    return "%s %s (%d bytes)" % (SEARCH_COMPLETED, os.path.basename(path), os.path.getsize(path))


def search_checkpoint(directory, step="Sequence search"):
    """
    (interrupted, {results file name: size when its search finished}) from
    checkpoint.txt in directory: whether the last search phase logged there
    was interrupted, and the searches recorded as finished (see above).
    step: the name of the search phase ("Search orthogroup profiles" for the
    --assign profile search).
    """
    SEARCH_STEP = "[%s] " % step
    interrupted = False
    done = {}
    fn = os.path.join(directory, "checkpoint.txt")
    if os.path.exists(fn):
        with open(fn, errors="replace", encoding="utf-8") as infile:
            for line in infile:
                if SEARCH_COMPLETED in line:
                    m = _SEARCH_COMPLETED_RE.search(line)
                    if m:
                        done[m.group(1)] = int(m.group(2))
                elif SEARCH_STEP in line:
                    status = line.split(SEARCH_STEP, 1)[1].strip()
                    if status == "Started":
                        interrupted = True
                    elif status == "Completed" or status.startswith("Workflow stopped normally"):
                        # finished, or -op stopping before any search ran
                        interrupted = False
    return interrupted, done


RUN_HEADER = "Starting OrthoFinder v"


def _is_search_record(line):
    return SEARCH_COMPLETED in line or SEARCH_STEP in line


def trim_checkpoint(path, completed=(), redo=()):
    """
    Trim checkpoint.txt (the record of every run of a results directory) back
    to the point a restart continues from, so that running the same restart
    again does not record the same stages twice:
      completed: the stages that must have completed before that point, in
        order, each as a tuple of alternative names (e.g. ("Infer
        orthogroups", "Infer clade-specific orthogroups"));
      redo: the stages the restart runs again; if not given, any stage
        started after those completed (or a later run's start).
    The lines from the first redone stage on are removed, except the search
    records ("[Sequence search] ...", "Search completed: ..."), which a -b
    restart uses to decide which search results are complete. Nothing is
    changed if the stages are not found (e.g. a checkpoint of an older version).
    """
    if not os.path.exists(path):
        return
    with open(path, errors="replace", encoding="utf-8") as infile:
        lines = infile.readlines()
    start = 0
    for alternatives in completed:
        markers = ["[%s] Completed" % name for name in alternatives]
        i = next((k for k in range(start, len(lines)) if any(m in lines[k] for m in markers)), None)
        if i is None:
            return
        start = i + 1
    if redo:
        starts = ["[%s] Started" % name for name in redo]
        cut = next((k for k in range(start, len(lines)) if any(m in lines[k] for m in starts)), None)
    else:
        cut = next((k for k in range(start, len(lines))
                    if "] Started" in lines[k] or RUN_HEADER in lines[k]), None)
    if cut is None:
        return
    kept = lines[:cut] + [line for line in lines[cut:] if _is_search_record(line)]
    # Of later runs (restarts) before the cut, only their search records are
    # kept (e.g. of searches a -b ran), so that repeating a restart records
    # nothing twice; the first run, which made the results directory, is kept.
    headers = [k for k in range(len(kept)) if RUN_HEADER in kept[k]]
    drop = set()
    for i, h in enumerate(headers):
        end = headers[i + 1] if i + 1 < len(headers) else len(kept)
        if h > 0 and h >= start:
            drop.update(k for k in range(h, end) if not _is_search_record(kept[k]))
    kept = [line for k, line in enumerate(kept) if k not in drop]
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as outfile:
        outfile.writelines(kept)
    os.replace(tmp, path)


def hit_file_is_empty(path, fmt=None):
    """Whether a search-results file has no hit lines (fmt: its layout, default hit_format(path))."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as infile:
        return next(hit_lines(infile, fmt or hit_format(path)), None) is None


# ------------------------------------------------------------ command line ----

def file_info(path):
    """
    Encoding, line endings, row and field counts of a tab-separated file, and
    how many of its lines have quotes (files of older versions quoted some fields).
    """
    opener = gzip.open if path.endswith(".gz") else open
    n_lf = n_crlf = n_lines = n_quote_lines = 0
    n_fields = set()
    utf8 = True
    with opener(path, "rb") as infile:
        for line in infile:
            n_lines += 1
            if line.endswith(b"\r\n"):
                n_crlf += 1
            elif line.endswith(b"\n"):
                n_lf += 1
            if b'"' in line:
                n_quote_lines += 1
            if utf8:
                try:
                    line.decode(ENCODING)
                except UnicodeDecodeError:
                    utf8 = False
            line = line.rstrip(b"\r\n")
            if line:
                n_fields.add(line.count(b"\t") + 1)
    if n_crlf and n_lf:
        endings = "mixed (%d LF, %d CRLF)" % (n_lf, n_crlf)
    else:
        endings = "CRLF (Windows)" if n_crlf else "LF" if n_lf else "none"
    return {
        "path": path,
        "encoding": "UTF-8" if utf8 else "not UTF-8",
        "lines": n_lines,
        "line endings": endings,
        "fields per line": ("%d" % min(n_fields) if len(n_fields) == 1 else
                            "%d-%d" % (min(n_fields), max(n_fields)) if n_fields else "-"),
        "lines with quotes": n_quote_lines,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Inspect tab-separated files (OrthoFinder's file_io).")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("info", help="encoding, line endings, rows, fields and quotes of each file")
    p.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    for path in args.files:
        for key, value in file_info(path).items():
            print("%-18s %s" % (key + ":", value))
        print()


if __name__ == "__main__":
    main()
