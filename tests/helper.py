import os 
import sys

def create_path(arg):
    filepath = os.path.abspath(arg)
    if not os.path.isfile(filepath) and filepath[-1] != os.sep:
        filepath += os.sep
    return filepath

def _latest_output_dir(results_dir) -> str:
    if isinstance(results_dir, list):
        results_dir = results_dir[0]

    results_dir = os.path.abspath(results_dir)

    try:
        entries = [
            (os.stat(os.path.join(results_dir, name)).st_mtime,
             os.path.join(results_dir, name))
            for name in os.listdir(results_dir)
        ]
        return sorted(entries)[-1][1] if entries else ""
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
