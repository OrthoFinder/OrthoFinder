"""Directory and file helpers: working directories, cleaning up, compressing."""

import datetime
import os
import shutil

try:
    from rich import print
except ImportError:
    ...

from .printing import printer
from . import file_io


def get_directory_name(
    base_dir_name,
    i,
    search_program,
    msa_program,
    tree_program,
    scorematrix,
    gapopen,
    gapextend,
    extended_filename,
):

    if not extended_filename:
        if i == 0:
            return base_dir_name + os.sep
        else:
            return base_dir_name + ("_%d" % i) + os.sep
    else:
        if scorematrix and os.path.isfile(scorematrix):
            matrix_name = os.path.basename(scorematrix).split(".")[0]
        else:
            matrix_name = scorematrix if scorematrix else ""

        if gapopen and gapextend:
            gap_penalty = "-" + "-".join((gapopen, gapextend))
        else:
            gap_penalty = ""

        if len(matrix_name) != 0 or len(gap_penalty) != 0:
            extension = (
                matrix_name
                + gap_penalty
                + "-"
                + "-".join((search_program, msa_program, tree_program))
            )
        else:
            extension = "-".join((search_program, msa_program, tree_program))

        if i == 0:
            return base_dir_name + "_" + extension + os.sep
        else:
            return base_dir_name + ("_%d" % i) + "_" + extension + os.sep


def create_new_working_directory(
    base_directory_name,
    q_date=True,
    search_program=None,
    msa_program=None,
    tree_program=None,
    scorematrix=None,
    gapopen=None,
    gapextend=None,
    extended_filename=False,
    makedir=True,
):
    date_str = datetime.date.today().strftime("%b%d") if q_date else ""
    i_append = 0
    new_directory_name = get_directory_name(
        base_directory_name + date_str,
        i_append,
        search_program,
        msa_program,
        tree_program,
        scorematrix,
        gapopen,
        gapextend,
        extended_filename,
    )
    if makedir:
        while os.path.exists(new_directory_name):
            i_append += 1
            new_directory_name = get_directory_name(
                base_directory_name + date_str,
                i_append,
                search_program,
                msa_program,
                tree_program,
                scorematrix,
                gapopen,
                gapextend,
                extended_filename,
            )
        os.mkdir(new_directory_name)
    return new_directory_name


# Unused
def get_unused_filename(base_filename, ext, i_unassigned=None):
    i_append = 0
    new_filename = base_filename + ext
    if i_unassigned is not None:
        while os.path.exists(new_filename):
            i_append += 1
            new_filename = base_filename + ("_%d" % i_append) + ext
    return new_filename, i_append


def file_open(filename_with_gz, mode, gz):
    """A text file (file_io.open_text): UTF-8, "\\n" line endings; filename_with_gz + ".gz" if gz."""
    return file_io.open_text(filename_with_gz + ".gz" if gz else filename_with_gz, mode)


def compress_files(src_dir, output_fn):
    try:
        with open(output_fn, "w", encoding="utf-8", newline="\n") as writer:
            filenames = []
            with os.scandir(src_dir) as it:
                for entry in it:
                    if entry.is_file():
                        # fn = entry.name.partition("_")[0]
                        # fn = entry.name.partition(".")[0]
                        filenames.append((entry.name.partition(".")[0], entry))
            filenames = sorted(filenames)
            for fn, entry in filenames:
                try:
                    with open(entry.path, "r", encoding="utf-8") as infile:
                        og_tree = "".join(line.strip() for line in infile)
                    writer.write(f"{fn}: {og_tree}\n")
                except IOError as e:
                    print(f"ERROR reading file {entry.name}: {e}")
    except IOError as e:
        printer.print(f"ERROR writing to output file {output_fn}: {e}", style="error")


def cleanup_path(path):
    try:
        if os.path.isfile(path) or os.path.islink(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))


def clear_dir(of3_dir):
    if os.path.exists(of3_dir):
        with os.scandir(of3_dir) as entries:
            for entry in entries:
                try:
                    if entry.is_file() or entry.is_symlink():
                        os.unlink(entry.path)
                    elif entry.is_dir():
                        shutil.rmtree(entry.path)

                except Exception as e:
                    printer.print(
                        f"Failed to delete {entry.path}. Reason: {e}", style="error"
                    )
