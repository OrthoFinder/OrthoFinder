import pytest 
from orthofinder.run import helpinfo

@pytest.mark.skipif(
    "config.getoption('--skip-of-test')",
    reason="Skipping OrthoFinder results test"
)
@pytest.mark.order(6)
@pytest.mark.unit
def test_print_help_info(of_obj, capsys):
    helpinfo.PrintHelp(of_obj.prog_caller)
    out, err = capsys.readouterr()
    assert "usage" in out.lower() or "help" in out.lower()