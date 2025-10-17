import os 
import sys
import re
import traceback
import pytest
from orthofinder.run.__main__ import main

FATAL_MARKERS = (
    "Traceback (most recent call last)",
    "CRITICAL",
    "FATAL",
    "Segmentation fault",
    "Specified directory doesn't exist",
    "NCBI C++ Exception",  # BLAST hard failure
)

NON_FATAL_PATTERNS = (
    r"\b0 errors?\b",
    r"\bno errors?\b",
    r"\berror rate\b",
    r"\bERRORS?:\s*0\b",
)

def clean_text(s: str) -> str:
    # strip ANSI color codes
    return re.sub(r"\x1b\[[0-9;]*m", "", s)


def _run_main(args, capfd):
    """
    Call main(args) in-process and return (exit_code, out, err, text).
    - Success: exit_code == 0
    - If main() returns None, we treat it as 0
    - If main() calls sys.exit(None), treat as 0
    """
    try:
        ret = main(args)
        code = 0 if ret is None else int(ret)
    except SystemExit as e:
        # sys.exit(0) or sys.exit(None) are both "success"
        code = e.code if isinstance(e.code, int) else 0
    except Exception:
        # Unexpected crash: fail with traceback + whatever was printed so far
        out, err = capfd.readouterr()
        tb = traceback.format_exc()
        pytest.fail(
            "Unexpected exception from main()\n\n"
            f"{tb}\n"
            f"--- Captured stdout ---\n{out}\n"
            f"--- Captured stderr ---\n{err}\n"
        )

    out, err = capfd.readouterr()
    text = out + err
    return code, out, err, text



def create_path(arg):
    filepath = os.path.abspath(arg)
    if not os.path.isfile(filepath) and filepath[-1] != os.sep:
        filepath += os.sep
    return filepath

def _latest_output_dir(results_dir, fileno=-1) -> str:
    if isinstance(results_dir, list):
        results_dir = results_dir[0]

    results_dir = os.path.abspath(results_dir)

    try:
        entries = [
            (os.stat(os.path.join(results_dir, name)).st_mtime,
             os.path.join(results_dir, name))
            for name in os.listdir(results_dir)
        ]
        return sorted(entries)[fileno][1] if entries else ""
    except FileNotFoundError:
        return ""
    except NotADirectoryError:
        return results_dir






# def get_dir_path(arg):
#     directory = os.path.abspath(arg)
#     if not os.path.isfile(directory) and directory[-1] != os.sep:
#         directory += os.sep
#     if not os.path.exists(directory):
#         print("Specified directory doesn't exist: %s" % directory)
#         sys.exit(1)
#     return directory

# def get_file_path(arg):
#     file_path = os.path.abspath(arg)
#     directory = os.path.dirname(file_path)
#     if not os.path.exists(directory):
#         print("Directory points to the file doesn't exist: %s" % directory)
#         sys.exit(1)
#     if not os.path.isfile(file_path):
#         print("Specified file doesn't exist: %s" % file_path)
#         sys.exit(1)
#     return file_path
