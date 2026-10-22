"""The default number of analysis threads (-a): a quarter of the search threads (-t), at least 1."""

import pytest

from orthofinder.run import process_args


@pytest.mark.parametrize("n_blast, expected", [(1, 1), (3, 1), (4, 1), (8, 2), (12, 3), (64, 16), (128, 32)])
def test_default_alg_threads(n_blast, expected):
    assert process_args.default_alg_threads(n_blast) == expected


def test_default_from_t_and_a_given(tmp_path):
    prog_caller = process_args.program_caller_with_configs([])
    options = process_args._read_options(["-f", str(tmp_path), "-t", "8"], prog_caller)[0]
    assert options.n_process_alg == 2
    options = process_args._read_options(["-f", str(tmp_path), "-t", "8", "-a", "5"], prog_caller)[0]
    assert options.n_process_alg == 5          # (given: not the default)
