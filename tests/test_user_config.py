"""user_config.json, the example of a user's own config (~/config_orthofinder_user.json)."""

import json
import os

from orthofinder.utils import program_caller

ROOT = os.path.join(os.path.dirname(__file__), "..")
USER_CONFIG = os.path.join(ROOT, "user_config.json")
CONFIG = os.path.join(ROOT, "src", "orthofinder", "run", "config.json")


def entries(fn):
    """The program entries of a config file (its "__comment" keys repeat, so read as pairs)."""
    with open(fn, encoding="utf-8") as infile:
        pairs = json.load(infile, object_pairs_hook=lambda p: p)
    return {name: dict(value) for name, value in pairs if isinstance(value, list)}


def test_user_config_uses_only_known_options():
    """A misspelt option (e.g. "ouput_filename") is ignored, so the program's results are not found."""
    known = {key for entry in entries(CONFIG).values() for key in entry}
    for name, entry in entries(USER_CONFIG).items():
        assert set(entry) <= known, "%s: unknown options %s" % (name, sorted(set(entry) - known))


def test_user_config_raxml_is_the_default_one():
    """The RAxML entry names the program as the bioconda package and the Makefile install it."""
    assert entries(USER_CONFIG)["raxml"] == entries(CONFIG)["raxml"]


def test_user_config_loads():
    prog_caller = program_caller.ProgramCaller(CONFIG)
    prog_caller.add(program_caller.ProgramCaller(USER_CONFIG))
    assert "raxml" in prog_caller.list_tree_methods()
