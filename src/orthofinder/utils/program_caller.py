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
import concurrent.futures
import os
import sys
import json
import numpy as np
import subprocess
from .. import my_env 
import types 
import traceback
try:
    import queue
except ImportError:
    import Queue as queue     

from . import util, parallel_task_manager

try:
    from importlib import resources as impresources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as impresources

from .. import test_sequences


import multiprocessing as mp



try:
    longer_file = (impresources.files(test_sequences) / 'longer.txt')
    shorter_file = (impresources.files(test_sequences) / 'shorter.txt')
    with longer_file.open("rt") as f1, shorter_file.open("rt") as f2:  # "rt" as text file with universal newlines
        longer_sequence = f1.read()
        shorter_sequence = f2.read()
except AttributeError:
    # Python < PY3.9, fall back to method deprecated in PY3.11.
    longer_sequence = impresources.read_text(test_sequences, 'longer.txt')
    shorter_sequence = impresources.read_text(test_sequences, 'shorter.txt')
except FileNotFoundError as e:
        print(f"File not found: {e.filename}")


class InvalidEntryException(Exception):
    pass


class Method(object):
    def __init__(self, name, config_dict):
        self.skip_check = False
        if 'cmd_line' in config_dict:
            self.cmd = config_dict['cmd_line']
        else:
            print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
            print("'cmd_line' entry is missing")
            raise InvalidEntryException
        if 'ouput_filename' in config_dict:
            self.non_default_outfn = config_dict['ouput_filename'] 
        else:
            self.non_default_outfn = None
        # Non-advertised methods, switch to a faster method if number of sequences is greater than X
        if "cmd_line_fast" in config_dict and "n_seqs_use_fast" in config_dict:
            self.cmd_fast = config_dict['cmd_line_fast']
            self.n_seqs_use_fast = int(config_dict['n_seqs_use_fast'])
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
        # Add default methods
        self.msa['mafft'] = Method('mafft', {"cmd_line": "mafft --localpair --maxiterate 1000 --anysymbol INPUT > OUTPUT", "cmd_line_fast": "mafft --anysymbol INPUT > OUTPUT", "n_seqs_use_fast" : "500"})
        self.tree['fasttree'] = Method('fasttree', {"cmd_line": "FastTree INPUT > OUTPUT"})
        if configure_file == None:
            return
        if not os.path.exists(configure_file):
            print(("WARNING: Configuration file, '%s', does not exist. No user-confgurable multiple sequence alignment or tree inference methods have been added.\n" % configure_file))
            return
        with open(configure_file, 'r') as infile:
            try:
                d = json.load(infile)
            except ValueError:
                print(("WARNING: Incorrectly formatted configuration file %s" % configure_file))
                print("File is not in .json format. No user-confgurable multiple sequence alignment or tree inference methods have been added.\n")
                return
            for name, v in d.items():
                if name == "__comment": continue
                if " " in name:
                    print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                    print(("No space is allowed in name: '%s'" % name))
                    continue
                    
                if 'program_type' not in v:
                    print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                    print("'program_type' entry is missing")
                try:
                    if v['program_type'] == 'msa':
                        if name in self.msa:
                            print(("Multiple sequence alignment method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.msa[name] = Method(name, v)
                    elif v['program_type'] == 'tree':
                        if name in self.tree:
                            print(("Tree inference method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.tree[name] = Method(name, v)
                    elif v['program_type'] == 'search':
                        if ('db_cmd' not in v) or ('search_cmd' not in v):
                            print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                            print("'cmd_line' entry is missing")
                            raise InvalidEntryException
                        if name in self.search_db:
                            print(("Sequence search method '%s' has already been defined, skipping config file entry." % name))
                        else:
                            self.search_db[name] = Method(name, {'cmd_line':v['db_cmd']})
                            self.search_search[name] = Method(name, {'cmd_line':v['search_cmd']})
                            if 'ouput_filename' in v:
                                print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                                print("'ouput_filename' option is not supported for 'program_type' 'search'")
                    else:
                        print(("WARNING: Incorrectly formatted configuration file entry: %s" % name))
                        print(("'program_type' should be 'msa' or 'tree', got '%s'" % v['program_type']))
                except InvalidEntryException:
                    pass
    
    def Add(self, other):
        self.msa.update(other.msa)
        self.tree.update(other.tree)
        self.search_db.update(other.search_db)
        self.search_search.update(other.search_search)  # search_db & search_search are only added together
    
    def ListMSAMethods(self):
        return [key for key in self.msa]
    def ListTreeMethods(self):
        return [key for key in self.tree]
    def ListSearchMethods(self):
        return [key for key in self.search_db]

    def GetMSAMethodCommand(self, method_name, infilename, outfilename_proposed, identifier, nSeqs=None, method_threads=None):
        return self._GetCommand('msa', method_name, infilename, outfilename_proposed, 
                                identifier, nSeqs=nSeqs, method_threads=method_threads)
    
    def GetTreeMethodCommand(self, method_name, infilename, outfilename_proposed, identifier, nSeqs=None, method_threads=None):
        return self._GetCommand('tree', method_name, infilename, outfilename_proposed, 
                                identifier, nSeqs=nSeqs, method_threads=method_threads)
    
    def GetSearchMethodCommand_DB(self, method_name, infilename, outfilename, 
                                  scorematrix=None, gapopen=None, gapextend=None, method_threads=None):
        return self._GetCommand('search_db', method_name, infilename, outfilename, 
                                 scorematrix=scorematrix, gapopen=gapopen, 
                                 gapextend=gapextend, method_threads=method_threads)[0]  # output filename isn't returned
    
    def GetSearchMethodCommand_Search(self, method_name, queryfilename, dbfilename, outfilename, \
                                      scorematrix=None, gapopen=None, gapextend=None, method_threads=None):
        return self._GetCommand('search_search', method_name, queryfilename, 
                                 outfilename, None, dbfilename, 
                                 scorematrix=scorematrix, gapopen=gapopen, 
                                 gapextend=gapextend, method_threads=method_threads)[0]  # output filename isn't returned
    
    def GetMSACommands(self, method_name, infn_list, outfn_list, id_list, nSeqs=None, method_threads=None):
        if nSeqs == None:        
            return [self.GetMSAMethodCommand(method_name, infn, outfn, ident, method_threads=method_threads) for infn, outfn, ident in zip(infn_list, outfn_list, id_list)]
        else:        
            return [self.GetMSAMethodCommand(method_name, infn, outfn, ident, n, method_threads=method_threads) for infn, outfn, ident, n in zip(infn_list, outfn_list, id_list, nSeqs)]
       
    def GetTreeCommands(self, method_name, infn_list, outfn_list, id_list, nSeqs=None, method_threads=None):        
        if nSeqs == None:        
            return [self.GetTreeMethodCommand(method_name, infn, outfn, ident, method_threads=method_threads) for infn, outfn, ident in zip(infn_list, outfn_list, id_list)]
        else:
            return [self.GetTreeMethodCommand(method_name, infn, outfn, ident, n, method_threads=method_threads) for infn, outfn, ident, n in zip(infn_list, outfn_list, id_list, nSeqs)]
   
    # not used      
    def GetSearchCommands_DB(self, method_name,  infn_list, outfn_list, 
                             scorematrix=None, gapopen=None, gapextend=None, method_threads=None):        
        return [self.GetSearchMethodCommand_DB(method_name, infn, outfn, 
                                               scorematrix=scorematrix, 
                                               gapopen=gapopen, 
                                               gapextend=gapextend, 
                                               method_threads=method_threads) for infn, outfn in zip(infn_list, outfn_list)]
    # not used 
    def GetSearchCommands_Search(self, method_name, querryfn_list, dblist, outfn_list,
                                 scorematrix=None, gapopen=None, gapextend=None, method_threads=None):        
        return [self.GetSearchMethodCommand_Search(method_name, querryfn, 
                                                   dbname, outfn,
                                                   scorematrix=scorematrix, 
                                                   gapopen=gapopen, 
                                                   gapextend=gapextend,
                                                   method_threads=method_threads) for querryfn, dbname, outfn in zip(querryfn_list, dblist, outfn_list)]
    
    # not used
    def CallMSAMethod(self, method_name, infilename, outfilename, identifier, nSeqs=None, method_threads=None):
        return self._CallMethod('msa', method_name, infilename, outfilename, identifier, nSeqs=nSeqs, method_threads=method_threads)
    # not used    
    def CallTreeMethod(self, method_name, infilename, outfilename, identifier, nSeqs=None, method_threads=None):
        return self._CallMethod('tree', method_name, infilename, outfilename, identifier, nSeqs=nSeqs, method_threads=method_threads)   
        
    def CallSearchMethod_DB(self, method_name, infilename, outfilename, 
                            scorematrix=None, gapopen=None, gapextend=None, method_threads=None):
        return self._CallMethod('search_db', method_name, infilename, outfilename, 
                                scorematrix=scorematrix, gapopen=gapopen, gapextend=gapextend,
                                method_threads=method_threads)     
        
    def CallSearchMethod_Search(self, method_name, queryfilename, 
                                dbfilename, outfilename,
                                scorematrix=None, gapopen=None, gapextend=None, method_threads=None):
        return self._CallMethod('search_search', method_name, queryfilename, 
                                outfilename, dbname=dbfilename, 
                                scorematrix=scorematrix, gapopen=gapopen, gapextend=gapextend,
                                method_threads=method_threads)      

    def TestMSAMethod(self, method_name, d_test, method_threads="1"):
        return self._TestMethod('msa', method_name, d_test, method_threads=method_threads)
        
    def TestTreeMethod(self, method_name, d_test, method_threads="1"):
        return self._TestMethod('tree', method_name, d_test, method_threads=method_threads)
        
    def TestSearchMethod(self, method_name, d_deps_check, 
                         scorematrix=None, gapopen=None, gapextend=None, method_threads="1"):
        success = False
        fasta = self._WriteTestSequence_Longer(d_deps_check)
        dbname = d_deps_check + method_name + "DBSpecies0"
        stdout_db, stderr_db, cmd_db = self.CallSearchMethod_DB(method_name, fasta, dbname,
                                                                 scorematrix=scorematrix, 
                                                                 gapopen=gapopen,
                                                                 gapextend=gapextend,
                                                                 method_threads=method_threads
                                                                )
        # it doesn't matter what file(s) it writes out the database to, only that we can use the database
        resultsfn = d_deps_check + "test_search_results.txt"
        stdout_s, stderr_s, cmd_s = self.CallSearchMethod_Search(method_name, fasta, 
                                                                 dbname, resultsfn,
                                                                 scorematrix=scorematrix,
                                                                 gapopen=gapopen,
                                                                 gapextend=gapextend,
                                                                 method_threads=method_threads)

        success = os.path.exists(resultsfn) or os.path.exists(resultsfn + ".gz")
        if not success:
            print("%s produced the following output:" % method_name)
            print("\n".join(stdout_db))
            print("\n".join(stderr_db))
            print("\n".join(stdout_s))
            print("\n".join(stderr_s))
        cmd = cmd_db + "\n" + cmd_s
        return success, stdout_db  + stdout_s, stderr_db + stderr_s, cmd
    
    def _CallMethod(self, method_type, method_name, infilename,
                    outfilename, identifier=None, dbname=None, nSeqs=None, 
                    scorematrix=None, gapopen=None, gapextend=None, method_threads=None):

        cmd, actual_target_fns = self._GetCommand(method_type, method_name, infilename, 
                                                  outfilename, identifier, dbname, nSeqs, 
                                                  scorematrix=scorematrix,
                                                  gapopen=gapopen, 
                                                  gapextend=gapextend,
                                                  method_threads=method_threads
                                                  )
        capture = subprocess.Popen(cmd, shell=True, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   env=parallel_task_manager.my_env)
        stdout = [x for x in capture.stdout]
        stderr = [x for x in capture.stderr]
        try:
            stdout = [x.decode() for x in stdout]
            stderr = [x.decode() for x in stderr]
        except (UnicodeDecodeError, AttributeError):
            stdout = [x.encode() for x in stdout]
            stderr = [x.encode() for x in stderr]
        capture.communicate()
        if actual_target_fns != None:
            actual, target = actual_target_fns
            if os.path.exists(actual): 
                os.rename(actual, target)
        return stdout, stderr, cmd
    
    def _ShouldSkipTest(self, method_type, method_name):
        if method_type == 'msa':
            dictionary = self.msa
        elif method_type == 'tree':
            dictionary = self.tree
        elif method_type == 'search_db':
            dictionary = self.search_db
        elif method_type == 'search_search':
            dictionary = self.search_search
        else:
            raise NotImplementedError
        if method_name not in dictionary:
            raise Exception("No %s method called '%s'" % (self._GetMethodTypeName(method_type), method_name))
        method_parameters = dictionary[method_name]
        return method_parameters.skip_check

    def _TestMethod(self, method_type, method_name, d_test, method_threads=None):
        util.PrintNoNewLine("Test can run \"%s\"" % method_name) 
        if self._ShouldSkipTest(method_type, method_name):
            print(" - test has been manually over-ridden")
            return True
        infn = self._WriteTestSequence(d_test)
        propossed_outfn = infn + ".output.txt"
        stdout, stderr, cmd = self._CallMethod(method_type, method_name, infn, propossed_outfn, "test", method_threads=method_threads)
        success = os.path.exists(propossed_outfn) and os.stat(propossed_outfn).st_size > 0 
        if success:
            print(" - ok")
        else:
            print(" - failed")
            print("".join(stdout))
            print("".join(stderr))
        return success, stdout, stderr, cmd

    def _ReplaceVariables(self, instring, infilename, outfilename, 
                          identifier=None, dbname=None, scorematrix=None, 
                          gapopen=None, gapextend=None, method_threads=None):

        path, basename = os.path.split(infilename)
        path_out, basename_out = os.path.split(outfilename)
        outstring = instring.replace("INPUT", infilename).\
                    replace("OUTPUT", outfilename).\
                    replace("BASENAME", basename).\
                    replace("PATH", path).\
                    replace("BASEOUTNAME", basename_out)

        if "diamond" in instring:
            if scorematrix and gapopen and gapextend:
                outstring = outstring.replace("SCOREMATRIX", scorematrix).\
                            replace("GAPOPEN", gapopen).\
                            replace("GAPEXTEND", gapextend)
        
        if identifier is not None:
            outstring = outstring.replace("IDENTIFIER", identifier)
        if dbname is not  None:
            outstring = outstring.replace("DATABASE", dbname)
        if method_threads is not None and "METHODTHREAD" in instring:
            outstring = outstring.replace("METHODTHREAD", method_threads)

        return outstring

    def _GetMethodTypeName(self, method_type):
        if method_type == 'msa':
            return "multiple sequence alignment"
        elif method_type == 'tree':
            return "tree"
        elif method_type == 'search_db' or method_type == 'search_search':
            return "alignment search"
        else:
            raise NotImplementedError
        
    def _GetCommand(self, method_type, method_name, 
                    infilename, outfilename_proposed, 
                    identifier=None, dbname=None,
                    nSeqs=None, scorematrix=None,
                    gapopen=None, gapextend=None, method_threads=None):
        """
        Returns:
            cmd, actual_target_fn
            Where:
                cmd - The command line that should be called
                actual_target_fn - None if the cmd will save the results file to outfilename_proposed 
                                   otherwise (actual_fn, outfilename_proposed)
        """
        
        if method_type == 'msa':
            dictionary = self.msa
        elif method_type == 'tree':
            dictionary = self.tree
        elif method_type == 'search_db':
            dictionary = self.search_db
        elif method_type == 'search_search':
            dictionary = self.search_search
        else:
            raise NotImplementedError
        if method_name not in dictionary:
            raise Exception("No %s method called '%s'" % (self._GetMethodTypeName(method_type), method_name))
        method_parameters = dictionary[method_name]

        if nSeqs != None and method_parameters.cmd_fast != None and nSeqs >= method_parameters.n_seqs_use_fast:
            cmd = self._ReplaceVariables(method_parameters.cmd_fast, infilename, 
                                         outfilename_proposed, identifier, dbname,
                                         scorematrix=scorematrix, 
                                         gapopen=gapopen,
                                         gapextend=gapextend,
                                         method_threads=method_threads
                                         )
        else:
            cmd = self._ReplaceVariables(method_parameters.cmd, infilename, 
                                         outfilename_proposed, identifier, dbname,
                                         scorematrix=scorematrix, 
                                         gapopen=gapopen,
                                         gapextend=gapextend,
                                         method_threads=method_threads
                                         )
        actual_target_fn = None
        if method_parameters.non_default_outfn:
            actual_fn = self._ReplaceVariables(method_parameters.non_default_outfn, 
                                               infilename, outfilename_proposed, identifier,
                                               scorematrix=scorematrix, 
                                               gapopen=gapopen,
                                               gapextend=gapextend,
                                               method_threads=method_threads
                                               )

            target_fn = outfilename_proposed
            actual_target_fn = (actual_fn, target_fn)
        # print((cmd, actual_target_fn))
        return cmd, actual_target_fn
        
    def _WriteTestSequence(self, working_dir):
        fn = working_dir + "Species0.fa"
        with open(fn, 'w') as outfile:
            outfile.write(shorter_sequence)
        return fn
        
    def _WriteTestSequence_Longer(self, working_dir):
        fn = working_dir + "Species0.fa"
        with open(fn, 'w') as outfile:
            outfile.write(longer_sequence)
        return fn
   
    @staticmethod
    def PrintDependencyCheckFailure(cmd):
        """Error message including environment & PATH if a dependency check fails
        """
        print("\nEnvironment:")
        print({key: parallel_task_manager.my_env[key] for key in sorted(parallel_task_manager.my_env)})
        print("\nCommand:")
        print("export PATH=%s:" % parallel_task_manager.my_env['PATH'])
        print(cmd)
        print("\nResolve any issues so that you can successfully run the above commands for OrthoFinder's dependencies on your computer and then re-run OrthoFinder.")


# ========================================================================================================================

def RunParallelCommands(nProcesses, commands, qListOfList,
                        method_threads=None,
                        method_threads_large=None,
                        method_threads_small=None, 
                        threshold=None, 
                        cmd_order="ascending",
                        tasksize=None,
                        qTrim=False,
                        q_print_on_error=False, 
                        q_always_print_stderr=False,
                        old_version=False):
    
    if qListOfList:
        commands_and_no_filenames = [[(cmd, None) for cmd in cmd_list] for cmd_list in commands]
    else:
        commands_and_no_filenames = [(cmd, None) for cmd in commands]
    RunParallelCommandsAndMoveResultsFile(nProcesses, 
                                          commands_and_no_filenames, 
                                          qListOfList,
                                          method_threads,
                                          method_threads_large,
                                          method_threads_small, 
                                          threshold,
                                          cmd_order,
                                          tasksize,
                                          qTrim,
                                          q_print_on_error,
                                          q_always_print_stderr,
                                          old_version
                                          )


# def RunParallelCommandsAndMoveResultsFile(nProcesses, commands_and_filenames, qListOfList, q_print_on_error=False,
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
#         futures = [executor.submit(parallel_task_manager.Worker_RunCommands_And_Move,
#                                    cmd_queue,
#                                    nProcesses,
#                                    i+1,
#                                    qListOfList,
#                                    q_print_on_error,
#                                    q_always_print_stderr=q_always_print_stderr)
#                    for _ in range(nProcesses)]
#     concurrent.futures.wait(futures)


def RunParallelCommandsAndMoveResultsFile(nProcesses, commands_and_filenames, qListOfList,
                                          method_threads=None,
                                          method_threads_large=None,
                                          method_threads_small=None, 
                                          threshold=None,
                                          cmd_order="ascending",
                                          tasksize=None,
                                          qTrim=False,
                                          q_print_on_error=False,
                                          q_always_print_stderr=False,
                                          old_version=False
                                          ):
    """
    Calls the commands in parallel and if required moves the results file to the required new filename
    Args:
        nProcess - the number of parallel process to use
        commands_and_filenames : tuple (cmd, actual_target_fn) where actual_target_fn = None if no move is required 
                                 and actual_target_fn = (actual_fn, target_fn) is actual_fn is produced by cmd and this 
                                 file should be moved to target_fn
        actual_target_fn - None if the cmd will save the results file to outfilename_proposed 
                           otherwise (actual_fn, outfilename_proposed)
        qListOfList - if False then commands_and_filenames is a list of (cmd, actual_target_fn) tuples
                      if True then commands_and_filenames is a list of lists of (cmd, actual_target_fn) tuples where the elements 
                      of the inner list need to be run in the order they appear.
        q_print_on_error - If error code returend print stdout & stederr
    """
    if old_version:
        print("\n*** You are running David's version of Multiprocessing ***\n")
        cmd_queue = queue.Queue()
        i = -1
        for i, cmd in enumerate(commands_and_filenames):
            # print(cmd)
            cmd_queue.put((i, cmd))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(parallel_task_manager.Worker_RunCommands_And_Move,
                                    cmd_queue,
                                    nProcesses,
                                    i+1,
                                    qListOfList,
                                    q_print_on_error,
                                    q_always_print_stderr=q_always_print_stderr)
                    for _ in range(nProcesses)]
        concurrent.futures.wait(futures)

    else:
        completed_count = 0
        total_commands = len(commands_and_filenames)
        if threshold is not None and tasksize is not None and len(tasksize) != 0:
            cutoff = np.percentile(tasksize, threshold)#.astype(int)

            if not qListOfList:
                if "METHODTHREAD" in commands_and_filenames[0][0]:
                    if method_threads_large is not None and method_threads_small is not None:
                        commands_and_filenames = [
                            (cmd[0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[1]) 
                            for cmd, size in zip(commands_and_filenames, tasksize)
                        ]
                    elif method_threads is not None:
                        commands_and_filenames = [
                            (cmd[0].replace("METHODTHREAD", method_threads), cmd[1]) 
                            for cmd in commands_and_filenames
                        ]
                    else:
                        commands_and_filenames = [
                            (cmd[0].replace("METHODTHREAD", "1"), cmd[1]) 
                            for cmd in commands_and_filenames
                        ]

            else:
                if "METHODTHREAD" in commands_and_filenames[0][0][0]:
                    if qTrim:
                        if method_threads_large is not None and method_threads_small is not None:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[0][1]), 
                                (cmd[1][0], cmd[1][1].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small)), 
                                (cmd[2][0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[2][1])] 
                                if len(cmd) > 1 
                                else [(cmd[0][0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[0][1])]
                                for cmd, size in zip(commands_and_filenames, tasksize)
                            ]
                        elif method_threads is not None:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1]), 
                                (cmd[1][0], cmd[1][1].replace("METHODTHREAD", method_threads)), 
                                (cmd[2][0].replace("METHODTHREAD", method_threads), cmd[2][1])] 
                                if len(cmd) > 1 
                                else [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1])] 
                                for cmd in commands_and_filenames
                            ]
                        else:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", "1"), cmd[0][1]), 
                                (cmd[1][0], cmd[1][1].replace("METHODTHREAD", "1")), 
                                (cmd[2][0].replace("METHODTHREAD", "1"), cmd[2][1])]
                                if len(cmd) > 1 
                                else [(cmd[0][0].replace("METHODTHREAD", "1"), cmd[0][1])]
                                for cmd in commands_and_filenames
                            ]
                    else:
                        if method_threads_large is not None and method_threads_small is not None:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[0][1]), 
                                (cmd[1][0].replace("METHODTHREAD", method_threads_large if size >= cutoff else method_threads_small), cmd[1][1])] 
                                for cmd, size in zip(commands_and_filenames, tasksize)
                            ]
                        elif method_threads is not None:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1]),  
                                (cmd[1][0].replace("METHODTHREAD", method_threads), cmd[1][1])] 
                                for cmd in commands_and_filenames
                            ]

                        else:
                            commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", "1"), cmd[0][1]), 
                                (cmd[1][0].replace("METHODTHREAD", "1"), cmd[1][1])]  
                                for cmd in commands_and_filenames
                            ]

            large_file_commands = [cmd for cmd, size in zip(commands_and_filenames, tasksize) if size >= cutoff]
            small_file_commands = [cmd for cmd, size in zip(commands_and_filenames, tasksize) if size < cutoff]
            
            nProcesses_small_files = nProcesses
            nProcesses_large_files = nProcesses

            if method_threads_small is not None and method_threads_large is not None:
                if nProcesses_small_files * int(method_threads_small) > mp.cpu_count():
                    nProcesses_small_files = mp.cpu_count() // int(method_threads_small)
                
                if nProcesses_large_files * int(method_threads_large) > mp.cpu_count():
                    nProcesses_large_files = mp.cpu_count() // int(method_threads_large)

            if method_threads is not None:
                if nProcesses * int(method_threads) > mp.cpu_count():
                    nProcesses_small_files = mp.cpu_count() // int(method_threads)
                    nProcesses_large_files = mp.cpu_count() // int(method_threads)

            nProcesses_files = np.amin((nProcesses_small_files, nProcesses_large_files))
            if cmd_order == "ascending":
                split_files_MP(nProcesses_files, 
                    small_file_commands + large_file_commands, 
                    completed_count, 
                    total_commands,
                    qListOfList,
                    q_print_on_error,
                    q_always_print_stderr)
                
            else:
                split_files_MP(nProcesses_files, 
                            large_file_commands + small_file_commands, 
                            completed_count, 
                            total_commands,
                            qListOfList,
                            q_print_on_error,
                            q_always_print_stderr)
                    

                
            # if small_file_commands:
                # with concurrent.futures.ThreadPoolExecutor(max_workers=nProcesses_small_files) as executor:
                #     futures = {executor.submit(Worker_RunCommands_And_Move, 
                #                             cmd, 
                #                             qListOfList,
                #                             q_print_on_error,
                #                             q_always_print_stderr
                #                             ): cmd for cmd in small_file_commands if cmd is not None}
                #     for future in concurrent.futures.as_completed(futures):
                #         try:
                #             result = future.result()

                #             completed_count += 1
                #             if divmod(completed_count, 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000)[1] == 0:
                #                 util.PrintTime("Done %d of %d" % (completed_count, total_commands))

                #             if result != 0 and q_print_on_error:
                #                 print(f"ERROR occurred with small-file command: {futures[future]}")
                #         except Exception as e:
                #             print(f"Exception with small-file command {futures[future]}: {e}")

            # if large_file_commands:
                # with concurrent.futures.ThreadPoolExecutor(max_workers=nProcesses_large_files) as executor:
                #     futures = {executor.submit(Worker_RunCommands_And_Move, 
                #                             cmd, 
                #                             qListOfList,
                #                             q_print_on_error,
                #                             q_always_print_stderr
                #                             ): cmd for cmd in large_file_commands if cmd is not None}
                #     for future in concurrent.futures.as_completed(futures):
                #         try:
                #             result = future.result()

                #             completed_count += 1
                #             if divmod(completed_count, 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000)[1] == 0:
                #                 util.PrintTime("Done %d of %d" % (completed_count, total_commands))

                #             if result != 0 and q_print_on_error:
                #                 print(f"ERROR occurred with large-file command: {futures[future]}")
                #         except Exception as e:
                #             print(f"Exception with large-file command {futures[future]}: {e}")

        else:
            if method_threads is None:
                method_threads = "1"

            if method_threads is not None:
                if nProcesses * int(method_threads) > mp.cpu_count():
                    nProcesses = mp.cpu_count() // int(method_threads)
                
            if not qListOfList:
                if "METHODTHREAD" in commands_and_filenames[0][0]:
                    commands_and_filenames = [(cmd[0].replace("METHODTHREAD", method_threads), cmd[1]) 
                                            for cmd in commands_and_filenames]

            else:
                if "METHODTHREAD" in commands_and_filenames[0][0][0]:
                    if qTrim:
                        commands_and_filenames = [
                            [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1]), 
                            (cmd[1][0], cmd[1][1].replace("METHODTHREAD", method_threads)), 
                            (cmd[2][0].replace("METHODTHREAD", method_threads), cmd[2][1])]
                            if len(cmd) > 1 
                            else [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1])] 
                            for cmd in commands_and_filenames
                        ]

                    else:
                        commands_and_filenames = [
                                [(cmd[0][0].replace("METHODTHREAD", method_threads), cmd[0][1]), 
                                (cmd[1][0].replace("METHODTHREAD", method_threads), cmd[1][1])]
                                for cmd in commands_and_filenames
                            ]
   
            with concurrent.futures.ThreadPoolExecutor(max_workers=nProcesses) as executor:
                futures = {executor.submit(Worker_RunCommands_And_Move, 
                                        cmd, 
                                        qListOfList, 
                                        q_print_on_error, 
                                        q_always_print_stderr): cmd for cmd in commands_and_filenames if cmd is not None}

                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        completed_count += 1 
                        if divmod(completed_count, 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000)[1] == 0:
                            util.PrintTime("Done %d of %d" % (completed_count, total_commands))

                        if result != 0 and q_print_on_error:
                            print(f"ERROR occurred with command: {futures[future]}")
                    except Exception as e:
                        print(f"Exception with command {futures[future]}: {e}")


def split_files_MP(nthreads, 
                   file_commands, 
                   completed_count, 
                   total_commands,
                    qListOfList,
                    q_print_on_error,
                    q_always_print_stderr):

    with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
        futures = {executor.submit(Worker_RunCommands_And_Move, 
                                cmd, 
                                qListOfList,
                                q_print_on_error,
                                q_always_print_stderr
                                ): cmd for cmd in file_commands if cmd is not None}
        
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()

                completed_count += 1
                if divmod(completed_count, 10 if total_commands <= 200 else 100 if total_commands <= 2000 else 1000)[1] == 0:
                    util.PrintTime("Done %d of %d" % (completed_count, total_commands))

                if result != 0 and q_print_on_error:
                    print(f"ERROR occurred with large-file command: {futures[future]}")
            except Exception as e:
                print(f"Exception with large-file command {futures[future]}: {e}")



q_print_first_traceback_0 = False
def Worker_RunCommands_And_Move(command_fns_list, qListOfLists, q_print_on_error, q_always_print_stderr):
    """
    Continuously takes commands that need to be run from the cmd_and_filename_queue until the queue is empty. If required, moves 
    the output filename produced by the cmd to a specified filename. The elements of the queue can be single cmd_filename tuples
    or an ordered list of tuples that must be run in the provided order.
  
    Args:
        cmd_and_filename_queue - queue containing (cmd, actual_target_fn) tuples (if qListOfLists is False) or a list of such 
            tuples (if qListOfLists is True). Alternatively, 'cmd' can be a python fn and actual_target_fn the fn to call it on.
        nProcesses - the number of processes that are working on the queue.
        nToDo - The total number of elements in the original queue
        qListOfLists - Boolean, whether each element of the queue corresponds to a single command or a list of ordered commands
        qShell - Boolean, should a shell be used to run the command.
        
    Implementation:
        nProcesses and nToDo are used to print out the progress.
    """
    # while True:
    try:
        # i, command_fns_list = cmd_and_filename_queue.get(True, 1)
        
        # nDone = i - nProcesses + 1
        # if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
        #     util.PrintTime("Done %d of %d" % (nDone, nToDo))
        if not qListOfLists:
            command_fns_list = [command_fns_list]
    
        return_code = 1
        for command, fns in command_fns_list:
            if isinstance(command, types.FunctionType):
                # This will block the process, but it is ok for trimming, it takes minimal time
                fn = command
                fn(fns)
                return_code = 0
            else:
                if not isinstance(command, str):
                    print("ERROR: Cannot run command: " + str(command))
                    print("Please report this issue.")
                else:
                    return_code = RunCommand(command, qPrintOnError=q_print_on_error, qPrintStderr=q_always_print_stderr)
                    if fns != None:
                        actual, target = fns
                        if os.path.exists(actual):
                            os.rename(actual, target)
        return return_code            
    # except queue.Empty:
    #     return     
    except Exception as e:
        print("WARNING: ")
        print(str(e))
        global q_print_first_traceback_0
        if not q_print_first_traceback_0:
            util.print_traceback(e)
            q_print_first_traceback_0 = True
    except:
        print("WARNING: Unknown caught unknown exception")


def RunCommand(command, qPrintOnError=False, qPrintStderr=True):
    """ Run a single command """
    popen = subprocess.Popen(command, env=my_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if qPrintOnError:
        stdout, stderr = popen.communicate()
        if popen.returncode != 0:
            print(("\nERROR: external program called by OrthoFinder returned an error code: %d" % popen.returncode))
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
        elif qPrintStderr and len(stderr) > 0 and not util.stderr_exempt(stderr):
            print("\nWARNING: program called by OrthoFinder produced output to stderr")
            print(("\nCommand: %s" % command))
            print(("\nstdout:\n%s" % stdout))
            print(("stderr:\n%s" % stderr))
        return popen.returncode
    else:
        popen.communicate()
        return popen.returncode