"""The options of features not yet released (cli_options.PREVIEW_OPTIONS): in the help and
usable on the development branch, neither in the released OrthoFinder (main)."""

import pytest

from orthofinder.run import cli_options, helpinfo, process_args

PREVIEW = ("--tree-view", "--stats-view", "--groups-view", "--view",
           "--matrix", "--custom-matrix", "-go", "--gapopen", "-ge", "--gapextend",
           "--orthoxml", "-x")


def test_the_preview_options():
    marked = {name for name, option in cli_options.OPTIONS_BY_NAME.items() if option.get("preview")}
    assert marked == set(PREVIEW)


def help_text(capsys):
    helpinfo.print_help(process_args.program_caller_with_configs([]))
    return capsys.readouterr().out


@pytest.mark.parametrize("on", [True, False])
def test_help(monkeypatch, capsys, on):
    monkeypatch.setattr(cli_options, "PREVIEW_OPTIONS", on)
    text = help_text(capsys)
    for name in ("--view", "--tree-view", "--matrix", "--custom-matrix", "--gapopen", "--gapextend", "--orthoxml"):
        assert (name in text) == on, name     # (also in the help of other options: -mog)
    assert ("ADDITIONAL OUTPUT OPTIONS" in text) == on     # (only preview options)
    assert "-y" in text and "-pof" in text                 # the others


@pytest.mark.parametrize("args", [["--view"], ["--stats-view"], ["--matrix", "BLOSUM45"], ["-go", "11"],
                                  ["--orthoxml"], ["-x"]])
def test_released_rejects_them(monkeypatch, tmp_path, args):
    monkeypatch.setattr(cli_options, "PREVIEW_OPTIONS", False)
    with pytest.raises(process_args.ArgumentError, match="not available in this version"):
        process_args._read_options(["-f", str(tmp_path)] + args, process_args.program_caller_with_configs([]))
