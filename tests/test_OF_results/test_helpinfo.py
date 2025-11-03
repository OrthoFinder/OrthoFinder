import pytest 
from orthofinder.run import helpinfo

@pytest.mark.order(3)
@pytest.mark.unit
def test_print_help_info(of_obj, capsys):
    helpinfo.PrintHelp(of_obj.prog_caller)
    out, err = capsys.readouterr()
    assert "usage" in out.lower() or "help" in out.lower()