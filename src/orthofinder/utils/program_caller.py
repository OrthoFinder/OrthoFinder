#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2016 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.comhor: david
import os
import json
import numpy as np
import subprocess
# from .. import my_env
import traceback
import re 
try:
    import queue
except ImportError:
    import Queue as queue

from . import util, parallel_task_manager
from .printing import printer
try:
    from rich import print
except ImportError:
    ...

import multiprocessing as mp

# def _to_text_or_none(x):
#     if isinstance(x, str):
#         return x
#     if isinstance(x, bytes):
#         try:
#             return x.decode("utf-8", "ignore")
#         except Exception:
#             return None
#     return None 

# def threads_from_cmd(cmd_maybe, method_threads: int) -> int:
#     s = _to_text_or_none(cmd_maybe)
#     if not s:
#         return 1
#     return int(method_threads) if _METHODTHREAD_RE.search(s) else 1

# def _detect_cores():
#     try:
#         return len(os.sched_getaffinity(0))
#     except Exception:
#         return mp.cpu_count() or 1


_SPLIT_CMDS_RE = re.compile(r"\s*(?:;|&&|\|)\s*")
_INPUT_FLAGS = {"-s", "--msa", "-q", "--query", "--in", "--input", "-align"}
_FASTA_EXTS = (".fa", ".fasta", ".faa", ".aln", ".fa.gz", ".fasta.gz", ".faa.gz", ".aln.gz")

# def _tokens(cmd):
#     try:
#         return shlex.split(cmd, posix=True)
#     except Exception:
#         return cmd.split()

# def _first_existing_path(tokens):
#     for t in tokens:
#         if t in {">", "1>", "2>"}:
#             break
#         if ("/" in t or t.lower().endswith(_FASTA_EXTS)) and os.path.exists(t):
#             return t
#     return None

# def extract_input_path(cmd: str):
#     """Return the most likely INPUT path for this (possibly compound) command."""
#     for sub in _SPLIT_CMDS_RE.split(cmd):
#         toks = _tokens(sub)
#         if not toks:
#             continue
#         for i, t in enumerate(toks[:-1]):
#             if t in _INPUT_FLAGS and os.path.exists(toks[i+1]):
#                 return toks[i+1]
#         prog = os.path.basename(toks[0]).lower()
#         if prog in {"famsa", "fasttree", "fasttreemp", "mafft", "mafft-memsave", "muscle"}:
#             inp = _first_existing_path(toks)
#             if inp: return inp
#     return None


# def _threads_by_size(path: str, base_threads: int) -> int:
#     try:
#         sz = os.path.getsize(path)
#     except Exception:
#         return 1
#     if sz <= 10000:        # < ~200 KB
#         return 1
#     return max(4, base_threads) 

# try:
#     longer_file = impresources.files(test_sequences) / "longer.txt"
#     shorter_file = impresources.files(test_sequences) / "shorter.txt"
#     with longer_file.open("rt") as f1, shorter_file.open(
#         "rt"
#     ) as f2:  # "rt" as text file with universal newlines
#         longer_sequence = f1.read()
#         shorter_sequence = f2.read()
# except AttributeError:
#     # Python < PY3.9, fall back to method deprecated in PY3.11.
#     longer_sequence = impresources.read_text(test_sequences, "longer.txt")
#     shorter_sequence = impresources.read_text(test_sequences, "shorter.txt")
# except FileNotFoundError as e:
#     print(f"File not found: {e.filename}")


longer_sequence = \
""">0_0
MNINSPNDKEIALKSYTETFLDILRQELGDQMLYKNFFANFEIKDVSKIGHITIGTTNVTPNSQYVIRAY
ESSIQKSLDETFERKCTFSFVLLDSAVKKKVKRERKEAAIENIELSNREVDKTKTFENYVEGNFNKEAIR
IAKLIVEGEEDYNPIFIYGKSGIGKTHLLNAICNELLKKEVSVKYINANSFTRDISYFLQENDQRKLKQI
RNHFDNADIVMFDDFQSYGIGNKKATIELIFNILDSRINQKRTTIICSDRPIYSLQNSFDARLISRLSMG
LQLSIDEPQKADLLKILDYMIDINKMTPELWEDDAKNFIVKNYANSIRSLIGAVNRLRFYNSEIVKTNSR
YTLAIVNSILKDIQQVKEKVTPDVIIEYVAKYYKLSRSEILGKSRRKDVVLARHIAIWIVKKQLDLSLEQ
IGRFFGNRDHSTIINAVRKIEKETEQSDITFKRTISEISNEIFKKN
>1_2
MKTKLKRFLEEISVHFNEANSELLDAFVHSIDFVFEENDNIYIYFESPYFFNEFKNKLNHLINVENAVVF
NDYLSLEWKKIIKENKRVNLLNKKEADTLKEKLATLKKQEKYKINPLSKGIKEKYNFGNYLVFEFNKEAV
YLAKQIANKTTHSNWNPIIIEGKPGYGKSHLLQAIANERQKLFPEEKICVLSSDDFGSEFLKSVIAPDPT
HIESFKSKYKDYDLLMIDDVQIISNRPKTNETFFTIFNSLVDQKKTIVITLDCKIEEIQDKLTARMISRF
QKGINVRINQPNKNEIIQIFKQKFKENNLEKYMDDHVIEEISDFDEGDIRKIEGSVSTLVFMNQMYGSTK
TKDQILKSFIEKVTNRKNLILSKDPKYVFDKIKYHFNVSEDVLKSSKRKKEIVQARHICMYVLKNVYNKN
LSQIGKLLRKDHTTVRHGIDKVEEELENDPNLKSFLDLFKN"""

shorter_sequence = \
""">A
MSKVIELKGIYAKYNKKSDYILEDLNLNVESGEFIAIIGPSGVGKSTLFKVIVNALEISKGSVRLFGQNI
>B
MLKLLSKFPLKVKLMALFAVILSTLHPFLSILIPTVTRQLITYLANSNINSEVSVYIFKSSWIIGSFSYA
>C
MQITVKDLVHTFLAKTPYELNAIDNINVTIKQGEFVGVIGQTGSGKTTFIEHLNALLLPSAGSVEWVFEN
>D
MIKVTDLMFKYPSAQANAIEKLNLEIESGKYVAILGHNGSGKSTFSKLLVALYKPADGKIELDGTTISKE"""


class InvalidEntryException(Exception):
    pass


# MMseqs2 steps that must be told nucleotide sequences are to be compared as
# nucleotides (--search-type 3, like blastn); on DNA they otherwise stop
# and ask for --search-type 2 (translated) or 3 (nucleotide).
MMSEQS_SEARCH_TYPE_RE = re.compile(r"(^|[\s;&|/])mmseqs\s+(createindex|search|easy-search|easy-linsearch|linsearch)\b")
_MMSEQS_CREATEINDEX_RE = re.compile(r"(^|[\s;&|/])mmseqs\s+createindex\b")
_KMER_OPTION_RE = re.compile(r"(^|\s)-k(\s|=|$)")
_SEARCH_TYPE_RE = re.compile(r"--search-type[\s=]+(\S+)")

# k-mer length for MMseqs2 nucleotide searches. At MMseqs2's default (15) every
# search allocates a k-mer table of ~8.4 GB, however small the input; 13 takes
# ~0.5 GB and found as many hits in tests (see the comment in config.json).
MMSEQS_NUCLEOTIDE_KMER = 13


def add_nucleotide_search_type(command):
    """
    command for nucleotide sequences: "--search-type 3" and
    "-k MMSEQS_NUCLEOTIDE_KMER" are added to each MMseqs2 search step that
    does not set them already (the k-mer length only for nucleotide search,
    type 3), and MMseqs2 createindex steps are left out. A nucleotide index holds a k-mer table
    of fixed size (8 GB on disk at the default k-mer length, however small
    the input), so each search builds its index in memory instead.
    Other programs' commands are returned unchanged.
    """
    parts = re.split(r"(;|&&|\|\||\|)", command)
    kept = []                  # [separator before it, command] pairs
    for i in range(0, len(parts), 2):     # the commands, not the separators
        part, sep = parts[i], parts[i - 1] if i else ""
        if _MMSEQS_CREATEINDEX_RE.search(part):
            continue
        if MMSEQS_SEARCH_TYPE_RE.search(part):
            stripped = part.rstrip()
            options = ""
            search_type = _SEARCH_TYPE_RE.search(part)
            if search_type is None:
                options += " --search-type 3"
            # The k-mer length is for nucleotide k-mers: not for a search type
            # the user chose (e.g. 2, translated, has amino-acid k-mers).
            nucleotide = search_type is None or search_type.group(1) == "3"
            if nucleotide and not _KMER_OPTION_RE.search(part):
                options += " -k %d" % MMSEQS_NUCLEOTIDE_KMER
            part = stripped + options + part[len(stripped):]
        kept.append((sep, part))
    if not kept:
        return command
    return "".join([kept[0][1].lstrip()] + [sep + part for sep, part in kept[1:]]).rstrip()


# One codon per amino acid, to make nucleotide test sequences from the protein ones.
_TEST_CODONS = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA", "E": "GAA",
    "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT",
    "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG", "Y": "TAT", "V": "GTT",
}


def _nucleotide_test_sequences():
    """longer_sequence back-translated (header lines kept)."""
    return "\n".join(
        line if line.startswith(">") else "".join(_TEST_CODONS.get(c, "NNN") for c in line.strip())
        for line in longer_sequence.split("\n")
    )


def _has_hit_lines(path):
    """Whether a (gzipped) search-results file has a hit line (not empty or only "#" comments)."""
    from . import file_io
    try:
        return not file_io.hit_file_is_empty(path, file_io.STANDARD_HIT_FORMAT)
    except (OSError, EOFError, UnicodeDecodeError):
        return False


class Method(object):
    def __init__(self, name, config_dict):
        self.skip_check = False
        if "cmd_line" in config_dict:
            self.cmd = config_dict["cmd_line"]
        else:
            print(
                ("WARNING: Incorrectly formatted configuration file entry: %s" % name)
            )
            print("'cmd_line' entry is missing")
            raise InvalidEntryException
        if "output_filename" in config_dict:
            self.non_default_outfn = config_dict["output_filename"]
        else:
            self.non_default_outfn = None
        # Non-advertised methods, switch to a faster method if number of sequences is greater than X
        if "cmd_line_fast" in config_dict and "n_seqs_use_fast" in config_dict:
            self.cmd_fast = config_dict["cmd_line_fast"]
            self.n_seqs_use_fast = int(config_dict["n_seqs_use_fast"])
        else:
            self.cmd_fast = None
            self.n_seqs_use_fast = None
        if "skip_check" in config_dict:
            if config_dict["skip_check"] == "true":
                self.skip_check = True


class ProgramCaller(object):
    def __init__(self, configure_file):
        self.msa = dict()
        self.tree = dict()
        self.search_db = dict()
        self.search_search = dict()
        # Optional, for a search program OrthoFinder cannot read the output
        # options of: the columns its results have, e.g. "qseqid sseqid bitscore".
        self.search_output_fields = dict()
        # Add default methods
        # self.msa["mafft"] = Method(
        #     "mafft",
        #     {
        #         "cmd_line": "mafft --localpair --maxiterate 1000 --anysymbol INPUT > OUTPUT",
                # "cmd_line_fast": "mafft --anysymbol INPUT > OUTPUT",
                # "n_seqs_use_fast": "500",
        #     },
        # )

        # self.tree["fasttree"] = Method(
        #     "fasttree", {"cmd_line": "FastTree INPUT > OUTPUT"}
        # )
        
        if configure_file == None:
            return
        if not os.path.exists(configure_file):
            print(("WARNING: Configuration file, '%s', does not exist. No user-confgurable multiple sequence alignment or tree inference methods have been added.\n" % configure_file))
            return
        with open(configure_file, "r") as infile:
            try:
                d = json.load(infile)
            except ValueError:
                print(("WARNING: Incorrectly formatted configuration file %s" % configure_file))
                print("File is not in .json format. No user-confgurable multiple sequence alignment or tree inference methods have been added.\n")
                return
            for name, v in d.items():
                if name == "__comment":
                    continue
                if " " in name:
                    print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                    print(("No space is allowed in name: '%s'" % name))
                    continue

                if "program_type" not in v:
                    print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                    print("'program_type' entry is missing")
                try:
                    if v["program_type"] == "msa":
                        if name in self.msa:
                            print(("Multiple sequence alignment method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.msa[name] = Method(name, v)
                    elif v["program_type"] == "tree":
                        if name in self.tree:
                            print(("Tree inference method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.tree[name] = Method(name, v)
                    elif v["program_type"] == "search":
                        if ("db_cmd" not in v) or ("search_cmd" not in v):
                            print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                            print("'cmd_line' entry is missing")
                            raise InvalidEntryException
                        if name in self.search_db:
                            print(("Sequence search method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.search_db[name] = Method(name, {"cmd_line": v["db_cmd"]})
                            self.search_search[name] = Method(name, {"cmd_line": v["search_cmd"]})
                            if "output_fields" in v:
                                self.search_output_fields[name] = v["output_fields"]
                            if "output_filename" in v:
                                print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                                print(
                                    "'output_filename' option is not supported for 'program_type' 'search'"
                                )
                    else:
                        print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                        print(("'program_type' should be 'msa' or 'tree', got '%s'" % v["program_type"]))
                except InvalidEntryException:
                    pass

    def add(self, other):
        self.msa.update(other.msa)
        self.tree.update(other.tree)
        self.search_db.update(other.search_db)
        self.search_search.update(
            other.search_search
        )  # search_db & search_search are only added together
        self.search_output_fields.update(other.search_output_fields)

    def list_msa_methods(self):
        return [key for key in self.msa]

    def list_tree_methods(self):
        return [key for key in self.tree]

    def get_search_output_fields(self, method_name):
        """The columns declared with "output_fields" for a search method, or None."""
        return self.search_output_fields.get(method_name)

    def list_search_methods(self):
        return [key for key in self.search_db]

    def get_msa_method_command(
        self,
        method_name,
        infilename,
        outfilename_proposed,
        identifier,
        n_seqs=None,
        method_threads=None,
    ):
        return self._get_command(
            "msa",
            method_name,
            infilename,
            outfilename_proposed,
            identifier,
            n_seqs=n_seqs,
            method_threads=method_threads,
        )

    def get_tree_method_command(
        self,
        method_name,
        infilename,
        outfilename_proposed,
        identifier,
        n_seqs=None,
        method_threads=None,
    ):
        return self._get_command(
            "tree",
            method_name,
            infilename,
            outfilename_proposed,
            identifier,
            n_seqs=n_seqs,
            method_threads=method_threads,
        )

    def get_search_method_command_db(
        self,
        method_name,
        infilename,
        outfilename,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
        nucleotide=False,
    ):
        """nucleotide: the sequences are DNA (see add_nucleotide_search_type)."""
        command = self._get_command(
            "search_db",
            method_name,
            infilename,
            outfilename,
            scorematrix=scorematrix,
            gapopen=gapopen,
            gapextend=gapextend,
            method_threads=method_threads,
        )[
            0
        ]  # output filename isn't returned
        return add_nucleotide_search_type(command) if nucleotide else command

    def get_search_method_command_search(
        self,
        method_name,
        queryfilename,
        dbfilename,
        outfilename,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
        nucleotide=False,
    ):
        """nucleotide: the sequences are DNA (see add_nucleotide_search_type)."""
        command = self._get_command(
            "search_search",
            method_name,
            queryfilename,
            outfilename,
            None,
            dbfilename,
            scorematrix=scorematrix,
            gapopen=gapopen,
            gapextend=gapextend,
            method_threads=method_threads,
        )[
            0
        ]  # output filename isn't returned
        return add_nucleotide_search_type(command) if nucleotide else command

    def uses_diamond(self, method_name):
        """Whether a configured search method runs DIAMOND (which compares protein sequences only)."""
        if method_name not in self.search_search:
            return False
        return re.search(r"(^|[\s;&|/])diamond\s", self.search_search[method_name].cmd) is not None

    def uses_make_blast_db(self, method_name):
        """Whether a configured search method builds BLAST+ databases (makeblastdb)."""
        if method_name not in self.search_db:
            return False
        return re.search(r"(^|[\s;&|/])makeblastdb\s", self.search_db[method_name].cmd) is not None

    def uses_mmseqs(self, method_name):
        """Whether a configured search method runs MMseqs2."""
        if method_name not in self.search_db:
            return False
        return any(MMSEQS_SEARCH_TYPE_RE.search(m.cmd)
                   for m in (self.search_db[method_name], self.search_search[method_name]))

    def get_msa_commands(
        self,
        method_name,
        infn_list,
        outfn_list,
        id_list,
        n_seqs=None,
        method_threads=None,
    ):
        if n_seqs == None:
            return [
                self.get_msa_method_command(
                    method_name, infn, outfn, ident, method_threads=method_threads
                )
                for infn, outfn, ident in zip(infn_list, outfn_list, id_list)
            ]
        else:
            return [
                self.get_msa_method_command(
                    method_name, infn, outfn, ident, n, method_threads=method_threads
                )
                for infn, outfn, ident, n in zip(infn_list, outfn_list, id_list, n_seqs)
            ]

    def get_tree_commands(
        self,
        method_name,
        infn_list,
        outfn_list,
        id_list,
        n_seqs=None,
        method_threads=None,
    ):
        if n_seqs == None:
            return [
                self.get_tree_method_command(
                    method_name, infn, outfn, ident, method_threads=method_threads
                )
                for infn, outfn, ident in zip(infn_list, outfn_list, id_list)
            ]
        else:
            return [
                self.get_tree_method_command(
                    method_name, infn, outfn, ident, n, method_threads=method_threads
                )
                for infn, outfn, ident, n in zip(infn_list, outfn_list, id_list, n_seqs)
            ]

    # # not used
    # def GetSearchCommands_DB(
    #     self,
    #     method_name,
    #     infn_list,
    #     outfn_list,
    #     scorematrix=None,
    #     gapopen=None,
    #     gapextend=None,
    #     method_threads=None,
    # ):
    #     return [
    #         self.get_search_method_command_db(
    #             method_name,
    #             infn,
    #             outfn,
    #             scorematrix=scorematrix,
    #             gapopen=gapopen,
    #             gapextend=gapextend,
    #             method_threads=method_threads,
    #         )
    #         for infn, outfn in zip(infn_list, outfn_list)
    #     ]

    # # not used
    # def GetSearchCommands_Search(
    #     self,
    #     method_name,
    #     querryfn_list,
    #     dblist,
    #     outfn_list,
    #     scorematrix=None,
    #     gapopen=None,
    #     gapextend=None,
    #     method_threads=None,
    # ):
    #     return [
    #         self.get_search_method_command_search(
    #             method_name,
    #             querryfn,
    #             dbname,
    #             outfn,
    #             scorematrix=scorematrix,
    #             gapopen=gapopen,
    #             gapextend=gapextend,
    #             method_threads=method_threads,
    #         )
    #         for querryfn, dbname, outfn in zip(querryfn_list, dblist, outfn_list)
    #     ]

    # # not used
    # def CallMSAMethod(
    #     self,
    #     method_name,
    #     infilename,
    #     outfilename,
    #     identifier,
    #     nSeqs=None,
    #     method_threads=None,
    # ):
    #     return self._call_method(
    #         "msa",
    #         method_name,
    #         infilename,
    #         outfilename,
    #         identifier,
    #         nSeqs=nSeqs,
    #         method_threads=method_threads,
    #     )

    # # not used
    # def CallTreeMethod(
    #     self,
    #     method_name,
    #     infilename,
    #     outfilename,
    #     identifier,
    #     nSeqs=None,
    #     method_threads=None,
    # ):
    #     return self._call_method(
    #         "tree",
    #         method_name,
    #         infilename,
    #         outfilename,
    #         identifier,
    #         nSeqs=nSeqs,
    #         method_threads=method_threads,
    #     )

    def call_search_method_db(
        self,
        method_name,
        infilename,
        outfilename,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
    ):
        return self._call_method(
            "search_db",
            method_name,
            infilename,
            outfilename,
            scorematrix=scorematrix,
            gapopen=gapopen,
            gapextend=gapextend,
            method_threads=method_threads,
        )

    def call_search_method_search(
        self,
        method_name,
        queryfilename,
        dbfilename,
        outfilename,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
    ):
        return self._call_method(
            "search_search",
            method_name,
            queryfilename,
            outfilename,
            dbname=dbfilename,
            scorematrix=scorematrix,
            gapopen=gapopen,
            gapextend=gapextend,
            method_threads=method_threads,
        )

    def test_msa_method(self, method_name, d_test, method_threads="1"):
        return self._test_method(
            "msa", method_name, d_test, method_threads=method_threads
        )

    def test_tree_method(self, method_name, d_test, method_threads="1"):
        return self._test_method(
            "tree", method_name, d_test, method_threads=method_threads
        )

    def test_search_method(
        self,
        method_name,
        d_deps_check,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads="1",
    ):
        """
        Run the search program's database and search commands on test
        sequences searched against themselves: protein sequences, then the
        same as nucleotides (for a nucleotide program such as blastn). It
        works if either search finds hits.
        """
        stdout_all, stderr_all = [], []
        for fasta_text in (longer_sequence, _nucleotide_test_sequences()):
            fasta = self._write_test_sequence_longer(d_deps_check, fasta_text)
            dbname = d_deps_check + method_name + "DBSpecies0"
            stdout_db, stderr_db, cmd_db = self.call_search_method_db(
                method_name,
                fasta,
                dbname,
                scorematrix=scorematrix,
                gapopen=gapopen,
                gapextend=gapextend,
                method_threads=method_threads,
            )
            # it doesn't matter what file(s) it writes out the database to, only that we can use the database
            resultsfn = d_deps_check + "test_search_results.txt"
            for fn in (resultsfn, resultsfn + ".gz"):
                if os.path.exists(fn):
                    os.remove(fn)        # from the previous test sequences
            stdout_s, stderr_s, cmd_s = self.call_search_method_search(
                method_name,
                fasta,
                dbname,
                resultsfn,
                scorematrix=scorematrix,
                gapopen=gapopen,
                gapextend=gapextend,
                method_threads=method_threads,
            )
            stdout_all += stdout_db + stdout_s
            stderr_all += stderr_db + stderr_s
            # A results file alone is not enough: with "... | gzip > OUTPUT.gz"
            # the shell creates it even if the program is not installed.
            success = any(
                os.path.exists(fn) and _has_hit_lines(fn) for fn in (resultsfn, resultsfn + ".gz")
            )
            if success:
                break
        if not success:
            print("%s produced the following output:" % method_name)
            print("\n".join(stdout_all))
            print("\n".join(stderr_all))
        # The test database and results are not needed again (an MMseqs2
        # index alone is ~0.9 GB, whatever the size of the test sequences).
        import glob
        for fn in ([dbname] + glob.glob(glob.escape(dbname) + ".*")
                   + [resultsfn, resultsfn + ".gz"]):
            if os.path.isfile(fn):
                os.remove(fn)
        cmd = cmd_db + "\n" + cmd_s
        return success, stdout_all, stderr_all, cmd

    def _call_method(
        self,
        method_type,
        method_name,
        infilename,
        outfilename,
        identifier=None,
        dbname=None,
        n_seqs=None,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
    ):

        cmd, actual_target_fns = self._get_command(
            method_type,
            method_name,
            infilename,
            outfilename,
            identifier,
            dbname,
            n_seqs,
            scorematrix=scorematrix,
            gapopen=gapopen,
            gapextend=gapextend,
            method_threads=method_threads,
        )
        # Read stdout and stderr together: reading one to the end first can
        # deadlock if the program fills the other pipe's buffer.
        _, out, err = parallel_task_manager.run_monitored_command(
            cmd,
            parallel_task_manager.my_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout = out.splitlines(True)
        stderr = err.splitlines(True)
        try:
            stdout = [x.decode() for x in stdout]
            stderr = [x.decode() for x in stderr]
        except (UnicodeDecodeError, AttributeError):
            stdout = [x.encode() for x in stdout]
            stderr = [x.encode() for x in stderr]
        if actual_target_fns != None:
            actual, target = actual_target_fns
            if os.path.exists(actual):
                os.rename(actual, target)
        return stdout, stderr, cmd

    def _should_skip_test(self, method_type, method_name):
        if method_type == "msa":
            dictionary = self.msa
        elif method_type == "tree":
            dictionary = self.tree
        elif method_type == "search_db":
            dictionary = self.search_db
        elif method_type == "search_search":
            dictionary = self.search_search
        else:
            raise NotImplementedError
        if method_name not in dictionary:
            raise Exception(
                "No %s method called '%s'"
                % (self._get_method_type_name(method_type), method_name)
            )
        method_parameters = dictionary[method_name]
        return method_parameters.skip_check

    def _test_method(self, method_type, method_name, d_test, method_threads=None):
        util.print_no_new_line(f'Test can run "[orange3]{method_name.split()[0]}[/orange3]"')
        # util.print_no_new_line('Test can run "%s"' % method_name)
        if self._should_skip_test(method_type, method_name):
            printer.print(" - test has been manually over-ridden", style="warning")
            return True
        infn = self._write_test_sequence(d_test)
        propossed_outfn = infn + f".{method_type}.{method_name}.output.txt"
        stdout, stderr, cmd = self._call_method(
            method_type,
            method_name,
            infn,
            propossed_outfn,
            f"test",
            method_threads=method_threads,
        )
        success = (
            os.path.exists(propossed_outfn) and os.stat(propossed_outfn).st_size > 0
        )
        if success:
            printer.print(" - [bold green]ok")
        else:
            printer.print(" - [bold red]failed")
            print("".join(stdout))
            print("".join(stderr))
        return success, stdout, stderr, cmd

    def _replace_variables(
        self,
        instring,
        infilename,
        outfilename,
        identifier=None,
        dbname=None,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
    ):

        path, basename = os.path.split(infilename)
        path_out, basename_out = os.path.split(outfilename)
        infilename = os.path.abspath(infilename)
        outfilename = os.path.abspath(outfilename)
        
        path = os.path.abspath(path)# + os.path.sep

        outstring = (
            instring.replace("INPUT", infilename)
            .replace("OUTPUT", outfilename)
            .replace("BASENAME", basename)
            .replace("PATH", path)
            .replace("BASEOUTNAME", basename_out)
        )

        if "diamond" in instring:
            if scorematrix and gapopen and gapextend:
                outstring = (
                    outstring.replace("SCOREMATRIX", scorematrix)
                    .replace("GAPOPEN", gapopen)
                    .replace("GAPEXTEND", gapextend)
                )

        if identifier is not None:
            outstring = outstring.replace("IDENTIFIER", identifier)
        if dbname is not None:
            outstring = outstring.replace("DATABASE", dbname)
        if method_threads is not None and "METHODTHREAD" in instring:
            outstring = outstring.replace("METHODTHREAD", method_threads)

        return outstring

    def _get_method_type_name(self, method_type):
        if method_type == "msa":
            return "multiple sequence alignment"
        elif method_type == "tree":
            return "tree"
        elif method_type == "search_db" or method_type == "search_search":
            return "alignment search"
        else:
            raise NotImplementedError

    def _get_command(
        self,
        method_type,
        method_name,
        infilename,
        outfilename_proposed,
        identifier=None,
        dbname=None,
        n_seqs=None,
        scorematrix=None,
        gapopen=None,
        gapextend=None,
        method_threads=None,
    ):
        """
        Returns:
            cmd, actual_target_fn
            Where:
                cmd - The command line that should be called
                actual_target_fn - None if the cmd will save the results file to outfilename_proposed
                                   otherwise (actual_fn, outfilename_proposed)
        """

        if method_type == "msa":
            dictionary = self.msa
        elif method_type == "tree":
            dictionary = self.tree
        elif method_type == "search_db":
            dictionary = self.search_db
        elif method_type == "search_search":
            dictionary = self.search_search
        else:
            raise NotImplementedError
        if method_name not in dictionary:
            raise Exception(
                "No %s method called '%s'"
                % (self._get_method_type_name(method_type), method_name)
            )
        method_parameters = dictionary[method_name]

        if (
            n_seqs != None
            and method_parameters.cmd_fast != None
            and n_seqs >= method_parameters.n_seqs_use_fast
        ):
            cmd = self._replace_variables(
                method_parameters.cmd_fast,
                infilename,
                outfilename_proposed,
                identifier,
                dbname,
                scorematrix=scorematrix,
                gapopen=gapopen,
                gapextend=gapextend,
                method_threads=method_threads,
            )
        else:
            cmd = self._replace_variables(
                method_parameters.cmd,
                infilename,
                outfilename_proposed,
                identifier,
                dbname,
                scorematrix=scorematrix,
                gapopen=gapopen,
                gapextend=gapextend,
                method_threads=method_threads,
            )
        actual_target_fn = None
        if method_parameters.non_default_outfn:
            actual_fn = self._replace_variables(
                method_parameters.non_default_outfn,
                infilename,
                outfilename_proposed,
                identifier,
                scorematrix=scorematrix,
                gapopen=gapopen,
                gapextend=gapextend,
                method_threads=method_threads,
            )
            target_fn = outfilename_proposed
            actual_target_fn = (actual_fn, target_fn)
        # print((cmd, actual_target_fn))
        return cmd, actual_target_fn

    def _write_test_sequence(self, working_dir):
        fn = working_dir + "Species0.fa"
        with open(fn, "w") as outfile:
            outfile.write(shorter_sequence)
        return fn

    def _write_test_sequence_longer(self, working_dir, fasta_text=None):
        fn = working_dir + "Species0.fa"
        with open(fn, "w") as outfile:
            outfile.write(longer_sequence if fasta_text is None else fasta_text)
        return fn

    @staticmethod
    def print_dependency_check_failure(cmd):
        """Error message including environment & PATH if a dependency check fails"""
        print("\nEnvironment:")
        print(
            {
                key: parallel_task_manager.my_env[key]
                for key in sorted(parallel_task_manager.my_env)
            }
        )
        print("\nCommand:")
        print("export PATH=%s:" % parallel_task_manager.my_env["PATH"])
        print(cmd)
        print(
            "\nResolve any issues so that you can successfully run the above commands for OrthoFinder's dependencies on your computer and then re-run OrthoFinder."
        )


# ========================================================================================================================


# def run_parallel_commands_and_move_results_file(nProcesses, commands_and_filenames, qListOfList, q_print_on_error=False,
#                                           q_always_print_stderr=False):
#     """
#     Calls the commands in parallel and if required moves the results file to the required new filename
#     Args:
#         nProcess - the number of parallel process to use
#         commands_and_filenames : tuple (cmd, actual_target_fn) where actual_target_fn = None if no move is required
#                                  and actual_target_fn = (actual_fn, target_fn) is actual_fn is produced by cmd and this
#                                  file should be moved to target_fn
#         actual_target_fn - None if the cmd will save the results file to outfilename_proposed
#                            otherwise (actual_fn, outfilename_proposed)
#         qListOfList - if False then commands_and_filenames is a list of (cmd, actual_target_fn) tuples
#                       if True then commands_and_filenames is a list of lists of (cmd, actual_target_fn) tuples where the elements
#                       of the inner list need to be run in the order they appear.
#         q_print_on_error - If error code returend print stdout & stederr
#     """
#     cmd_queue = Queue()
#     i = -1
#     for i, cmd in enumerate(commands_and_filenames):
#         cmd_queue.put((i, cmd))

#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         futures = [executor.submit(parallel_task_manager.worker_run_commands_and_move,
#                                    cmd_queue,
#                                    nProcesses,
#                                    i+1,
#                                    qListOfList,
#                                    q_print_on_error,
#                                    q_always_print_stderr=q_always_print_stderr)
#                    for _ in range(nProcesses)]
#     concurrent.futures.wait(futures)


