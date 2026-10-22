"""
The scoring matrices of DIAMOND searches (--matrix, --custom-matrix) and their
gap penalties (--gapopen, --gapextend).

DIAMOND itself checks them too (the dependency check runs a search with them),
but checking them here says which values are allowed before anything runs.
"""
import os

# The matrices DIAMOND supports: (default gap open, default gap extend), and
# for each gap extend penalty, the gap open penalties allowed with it. As
# DIAMOND 2.0.13 accepts them (tests/test_diamond_scoring.py checks them
# against the DIAMOND OrthoFinder runs); the defaults are NCBI's.
DIAMOND_MATRICES = {
    "BLOSUM45": ((14, 2), {1: range(16, 20), 2: range(12, 17), 3: range(10, 14)}),
    "BLOSUM50": ((13, 2), {1: range(15, 20), 2: range(12, 17), 3: range(9, 14)}),
    "BLOSUM62": ((11, 1), {1: range(9, 14), 2: range(6, 12)}),
    "BLOSUM80": ((10, 1), {1: range(9, 12), 2: (6, 7, 8, 9, 13, 25)}),
    "BLOSUM90": ((10, 1), {1: range(9, 12), 2: range(6, 10)}),
    "PAM250": ((14, 2), {1: range(17, 22), 2: range(13, 18), 3: range(11, 16)}),
    "PAM70": ((10, 1), {1: range(9, 12), 2: range(6, 9)}),
    "PAM30": ((9, 1), {1: range(8, 11), 2: range(5, 8)}),
}
DEFAULT_MATRIX = "BLOSUM62"

# The search methods of config.json for a custom matrix (--custom-matrix SCOREMATRIX)
CUSTOM_MATRIX_METHODS = {"diamond": "diamond_custom", "diamond_ultra_sens": "diamond_ultra_sens_custom"}


# The amino acids a custom matrix must have scores for: DIAMOND reads a matrix
# without one of them, but then finds no hits.
STANDARD_AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

CUSTOM_MATRIX_FORMAT = (
    "A custom matrix is a text file in the format of NCBI's matrices "
    "(https://ftp.ncbi.nlm.nih.gov/blast/matrices/): a line of one-letter amino acid codes, "
    "then for each code, in the same order, a line of the code and its integer score with each "
    "code of the first line. It needs the 20 standard amino acids (others, e.g. B, Z, X and *, "
    "are optional). Lines starting with # are comments.")


def read_custom_matrix(fn):
    """
    The lines of a custom matrix file as DIAMOND needs them (DIAMOND rejects a
    line with spaces at its end, as NCBI's matrix files have), and warnings.
    Raises ValueError saying what is wrong with the file.
    """
    if not os.path.isfile(fn):
        raise ValueError("The custom scoring matrix file does not exist: %s" % fn)

    def wrong(i_line, text):
        return ValueError("Custom scoring matrix %s, line %d: %s\n%s" % (fn, i_line, text, CUSTOM_MATRIX_FORMAT))

    with open(fn, errors="replace", encoding="utf-8") as infile:
        lines = [line.rstrip() for line in infile]
    rows = [(i + 1, line.split()) for i, line in enumerate(lines) if line.strip() and not line.startswith("#")]
    if not rows:
        raise ValueError("Custom scoring matrix %s has no scores\n%s" % (fn, CUSTOM_MATRIX_FORMAT))
    i_header, codes = rows[0]
    codes = [c.upper() for c in codes]
    if any(len(c) != 1 for c in codes):
        raise wrong(i_header, "the first line should be the one-letter amino acid codes")
    if len(set(codes)) != len(codes):
        raise wrong(i_header, "an amino acid code is repeated")
    missing = [a for a in STANDARD_AMINO_ACIDS if a not in codes]
    if missing:
        raise wrong(i_header, "no scores for %s (DIAMOND would find no hits)" % ", ".join(missing))
    if len(rows) - 1 != len(codes):
        raise ValueError("Custom scoring matrix %s has %d lines of scores, for %d amino acids\n%s"
                         % (fn, len(rows) - 1, len(codes), CUSTOM_MATRIX_FORMAT))
    scores = {}
    for (i_line, row), code in zip(rows[1:], codes):
        if row[0].upper() != code:
            raise wrong(i_line, "the line for %s was expected (the lines in the order of the codes of the first "
                                "line), got %s" % (code, row[0]))
        if len(row) - 1 != len(codes):
            raise wrong(i_line, "%d scores, for %d amino acids" % (len(row) - 1, len(codes)))
        try:
            scores[code] = [int(v) for v in row[1:]]
        except ValueError:
            raise wrong(i_line, "the scores should be integers")
    warnings = []
    asymmetric = [(a, b) for i, a in enumerate(codes) for j, b in enumerate(codes) if i < j and scores[a][j] != scores[b][i]]
    if asymmetric:
        a, b = asymmetric[0]
        warnings.append("The custom scoring matrix is not symmetric (e.g. %s,%s and %s,%s), which DIAMOND accepts: %s"
                        % (a, b, b, a, fn))
    return [line + "\n" for line in lines], warnings


CUSTOM_MATRIX_FN = "custom_scoring_matrix.txt"


def copy_custom_matrix(fn, working_dir):
    """
    The copy of a custom matrix in the working directory that the searches use
    (as DIAMOND reads it), so it is kept with the results for restarts.
    """
    copy_fn = os.path.join(working_dir, CUSTOM_MATRIX_FN)
    if os.path.abspath(fn) != os.path.abspath(copy_fn):
        lines, _ = read_custom_matrix(fn)
        with open(copy_fn, "w", encoding="utf-8", newline="\n") as outfile:
            outfile.writelines(lines)
    return copy_fn


def is_custom_matrix(score_matrix):
    """Whether a scoring matrix is a file (--custom-matrix) rather than one of DIAMOND's."""
    return score_matrix.upper() not in DIAMOND_MATRICES


def _values(values):
    """The allowed values, as text: ranges as first-last."""
    values = list(values)
    if values == list(range(values[0], values[-1] + 1)) and len(values) > 2:
        return "%d-%d" % (values[0], values[-1])
    return ", ".join(map(str, values))


def scoring(matrix=None, custom_matrix=None, gapopen=None, gapextend=None):
    """
    (score_matrix, gapopen, gapextend), as text for the search commands, from
    those given (None: not given): a missing gap penalty is the matrix's
    default, or with only --gapextend, the largest gap open allowed with it.
    Raises ValueError, saying what is allowed, for values DIAMOND does not accept.
    """
    if matrix is not None and custom_matrix is not None:
        raise ValueError("--matrix and --custom-matrix cannot be used together")
    for name, value in (("--gapopen", gapopen), ("--gapextend", gapextend)):
        if value is not None and value < 0:
            raise ValueError("%s must be 0 or more, got %d" % (name, value))
    if custom_matrix is not None:
        if gapopen is None or gapextend is None:
            raise ValueError("--custom-matrix also needs --gapopen and --gapextend (DIAMOND has no default "
                             "gap penalties for a custom matrix)")
        return os.path.abspath(custom_matrix), str(gapopen), str(gapextend)

    name = (matrix or DEFAULT_MATRIX).upper()
    if name not in DIAMOND_MATRICES:
        raise ValueError("Unknown scoring matrix for DIAMOND: %s\nValid matrices are: %s (or a matrix file with "
                         "--custom-matrix)" % (matrix, ", ".join(DIAMOND_MATRICES)))
    (default_open, default_extend), allowed = DIAMOND_MATRICES[name]
    extend = default_extend if gapextend is None else gapextend
    if extend not in allowed:
        raise ValueError("Gap extend penalty %d is not supported by DIAMOND for %s. Valid gap extend penalties: %s"
                         % (extend, name, ", ".join(map(str, allowed))))
    if gapopen is None:
        gapopen = default_open if extend == default_extend else max(allowed[extend])
    if gapopen not in allowed[extend]:
        raise ValueError("Gap open penalty %d is not supported by DIAMOND for %s with gap extend %d. Valid gap "
                         "open penalties: %s" % (gapopen, name, extend, _values(allowed[extend])))
    return name, str(gapopen), str(extend)
