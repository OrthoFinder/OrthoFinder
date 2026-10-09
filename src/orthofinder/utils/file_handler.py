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


"""
Handles location of all input and output files
Users:
1. Call initialise_file_handler
2. Interact with FileHandler object (which is initialised by the above method)
    - it is an instance of _FileHandler

The code also contains the help class: PreviousFilesLocator (and child classes of it)
"""

from . import logging as run_logging
import os
import sys
import glob
import time
import shutil
import datetime
try:
    from rich import print
except ImportError:
    ...
from . import util
from . import run_state
from .. import __version__

# Log.txt line of the species tree an analysis uses (FileHandler.log_species_tree_used).
SPECIES_TREE_LOG_KEY = "FN_SpeciesTree: "
# --assign: the first orthogroup whose gene tree was not inferred to place the
# new species (the trees of those before it are reused by a restart).
ASSIGN_FIRST_OG_LOG_KEY = "Assign_FirstOrthogroupToInfer: "

# Log.txt lines of the programs a run used (start_log), as the run_state settings.
PROGRAM_LOG_KEYS = (("Search program: ", "search_program"), ("MSA program: ", "msa_program"),
                    ("Tree program: ", "tree_program"))


def state_from_log(continuation_dir):
    """
    What a restart continues with, from the Log.txt of the results directory
    it is given (or of the one holding that WorkingDirectory), in the form of
    a run_state state, for results directories without run_state.json (older
    versions, or the file removed): the programs used, the species tree
    (FN_SpeciesTree) and, for --assign, the first orthogroup whose gene tree
    was not inferred to place the new species. The last record of each counts.
    Read before the restart writes its own Log.txt. {} if there is no Log.txt.
    """
    d = os.path.abspath(continuation_dir)
    log_fn = next((os.path.join(r, "Log.txt") for r in (d, os.path.dirname(d))
                   if os.path.exists(os.path.join(r, "Log.txt"))), None)
    if log_fn is None:
        return {}
    state = {}
    with open(log_fn, errors="replace") as infile:
        for line in infile:
            line = line.rstrip("\r\n")
            for key, name in PROGRAM_LOG_KEYS:
                if line.startswith(key):
                    state.setdefault("settings", {})[name] = line[len(key):].strip()
            if line.startswith(SPECIES_TREE_LOG_KEY):
                state.setdefault("files", {})["species_tree"] = PreviousFilesLocatorNew._in_this_results_dir(
                    line[len(SPECIES_TREE_LOG_KEY):], log_fn, 2)
            if line.startswith(ASSIGN_FIRST_OG_LOG_KEY):
                state.setdefault("assign", {})["i_og_restart"] = int(line[len(ASSIGN_FIRST_OG_LOG_KEY):])
    return state


# Log.txt lines of the earlier working directories an analysis reads from,
# after its own (WorkingDirectory_Base): those of the core analysis (--assign),
# or of the analysis it added species to (-b with -f).
EARLIER_WD_KEYS = ("WorkingDirectory_Core", "WorkingDirectory_Previous")

class _FileHandler(object):    
    def __init__(self):
        self.reset() 

    def reset(self):
        self.base_og_format = "OG%07d"
        self.wd_base = []               # Base: blast, species & sequence IDs, species fasta files - should not request this and then write here
        self.wd_current = None          # Location to write out any new files
        self.wd_trees = None            # Location of working dir containing tree files 
        self.rd1 = None
        self.file_identifier_string = "OrthoFinder"
        self.clusters_filename = None
        self.i_results_version = 0
        self.nondefault_pickle_dir = None
        self.species_tree_rooted_ids_fn = None
        self.species_un_tree_rooted_ids_fn = None
        self.multiple_rooted_species_trees_dir = None
        self.species_ids_corrected = None
        self.q_assign = False
     
    """ ========================================================================================== """
    # RefactorDS - FileHandler
    def create_output_dir_from_start_new(self, fasta_dir, base, 
                                     user_name = None, 
                                     old_wd_base_list=None,
                                     search_program=None,
                                     msa_program=None, 
                                     tree_program=None,
                                     scorematrix=None,
                                     gapopen=None,
                                     gapextend=None,
                                     extended_filename=False):
        """
        The initial difference will be that results will go in OrthoFinder/Results_DATE or USER_SPECIFIED/RESULTS_DATE
        whereas before they went in Results_DATE or USER_SPECIFIED.
        
        If this is a composite analysis (-f + -b) then old_wd_base_list != None 
        
        old_wd_base_list - first item is the WD from a previous analysis to be extended. If this extended other
          ones itself then there will be other items in the list.
        """
        if user_name == None:
            self.rd1 = util.create_new_working_directory(base + "Results_", 
                                                      search_program=search_program,
                                                      msa_program=msa_program,
                                                      tree_program=tree_program,
                                                      scorematrix=scorematrix,
                                                      gapopen=gapopen,
                                                      gapextend=gapextend,
                                                      extended_filename=extended_filename)
        else:
            self.rd1 = util.create_new_working_directory(base + "Results_" + user_name, 
                                                      q_date=False,
                                                      search_program=search_program,
                                                      msa_program=msa_program,
                                                      tree_program=tree_program,
                                                      scorematrix=scorematrix,
                                                      gapopen=gapopen,
                                                      gapextend=gapextend,
                                                      extended_filename=extended_filename)
        self.wd_current = os.path.join(self.rd1, "WorkingDirectory") + os.sep
        os.mkdir(self.wd_current)
        self.wd_base = [self.wd_current]        
        if old_wd_base_list != None:
            shutil.copy(old_wd_base_list[0] + "SpeciesIDs.txt", self.wd_current + "SpeciesIDs.txt")
            shutil.copy(old_wd_base_list[0] + "SequenceIDs.txt", self.wd_current + "SequenceIDs.txt")
            # The earlier working directories are recorded in run_state.json and
            # Log.txt (start_log); older versions wrote previous_wd.txt instead.
            self.wd_base.extend(old_wd_base_list)
        self.wd_trees = self.wd_current
        self.start_log(search_program=search_program, msa_program=msa_program, tree_program=tree_program,
                      scorematrix=scorematrix, gapopen=gapopen, gapextend=gapextend)
        
    # RefactorDS - PreviousFilesLocator
    def start_from_orthogroups_or_sequence_search(
            self, 
            wd_base_list, 
            base, 
            clusters_filename_pairs=None, 
            user_name = None, 
            user_species_tree=None,
            search_program=None,
            msa_program=None, 
            tree_program=None,
            scorematrix=None,
            gapopen=None,
            gapextend=None,
            extended_filename=False,
            use_base = False,
        ):
        """
        NEed to initialise:
        wd_base
        wd_trees
        wd_current
        """
        if len(self.wd_base) != 0: raise Exception("Changing WorkingDirectory1")
        self.wd_base = wd_base_list
        if use_base:
            self.wd_current = wd_base_list[0]
            if clusters_filename_pairs != None: self.clusters_filename = clusters_filename_pairs[:-len("_id_pairs.txt")]
            self.rd1 = os.path.dirname(self.wd_current[:-1]) + os.sep
            os.makedirs(self.wd_current, exist_ok=True)
        else:    
            if clusters_filename_pairs != None: self.clusters_filename = clusters_filename_pairs[:-len("_id_pairs.txt")]
            if user_name == None:
                self.rd1 = util.create_new_working_directory(base + "Results_",
                                                        search_program=search_program,
                                                        msa_program=msa_program,
                                                        tree_program=tree_program, 
                                                        scorematrix=scorematrix,
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            else:
                self.rd1 = util.create_new_working_directory(base + "Results_" + user_name,
                                                        search_program=search_program, 
                                                        q_date=False,
                                                        scorematrix=scorematrix,
                                                        msa_program=msa_program,
                                                        tree_program=tree_program, 
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            self.wd_current = os.path.join(self.rd1, "WorkingDirectory")  + os.sep
            os.mkdir(self.wd_current)
        self.wd_trees = self.wd_current
        self.start_log(search_program=search_program, msa_program=msa_program, tree_program=tree_program,
                      scorematrix=scorematrix, gapopen=gapopen, gapextend=gapextend)
    
    
    def start_from_trees(self, 
                       wd1_list, 
                       wd2,
                       base, 
                       clusters_filename_pairs,
                       species_tree_fn, 
                       q_is_user_species_tree,
                       q_start_from_species_trees=False,
                       user_name=None,
                       search_program=None,
                       msa_program=None,
                       tree_program=None,
                       scorematrix=None,
                       gapopen=None,
                       gapextend=None,
                       extended_filename=False,
                       use_base=False):
        """
        Convert user species tree here if necessary
        For OF species tree copy it to location given by FileHandler
        For user species tree, this must be done immediately by OF code
        """
        self.wd_base = wd1_list
        self.wd_trees = wd2

        if use_base:
            self.wd_current = wd1_list[0]
            if clusters_filename_pairs != None: self.clusters_filename = clusters_filename_pairs[:-len("_id_pairs.txt")]
            self.rd1 = os.path.dirname(self.wd_current[:-1]) + os.sep
            os.makedirs(self.wd_current, exist_ok=True)
        else:    
            if user_name == None:
                self.rd1 = util.create_new_working_directory(base + "Results_",
                                                        search_program=search_program, 
                                                        msa_program=msa_program,
                                                        tree_program=tree_program,
                                                        scorematrix=scorematrix,
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            else:
                self.rd1 = util.create_new_working_directory(base + "Results_" + user_name, 
                                                        q_date=False,
                                                        search_program=search_program,
                                                        msa_program=msa_program,
                                                        tree_program=tree_program,
                                                        scorematrix=scorematrix,
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            self.wd_current = os.path.join(self.rd1, "WorkingDirectory") + os.sep
            os.mkdir(self.wd_current)
        self.clusters_filename = clusters_filename_pairs[:-len("_id_pairs.txt")]
        self.start_log(search_program=search_program, msa_program=msa_program, tree_program=tree_program,
                      scorematrix=scorematrix, gapopen=gapopen, gapextend=gapextend)
        if q_start_from_species_trees:
            if not q_is_user_species_tree:
                shutil.copy(species_tree_fn, self.get_species_tree_ids_rooted_fn())
            self.write_to_log("Species Tree: %s\n" % species_tree_fn)
        else:
            self.write_to_log("Unrooted Species Tree: %s\n" % species_tree_fn)
        self.log_working_directory_trees()
                                         
    def create_output_directories(self, options, previous_files_locator, base_dir, fasta_dir=None):
        # --assign; a restart reads it from run_state.json (start_log), or Log.txt
        self.q_assign = bool(options.q_fast_add) or getattr(previous_files_locator, "q_assign", False)
        if (options.q_start_from_fasta and options.q_start_from_blast) or options.q_fast_add:
            wd1 = previous_files_locator.get_start_from_blast()
            self.create_output_dir_from_start_new(fasta_dir, base_dir, 
                                              user_name=options.name, 
                                              old_wd_base_list = wd1,
                                              search_program=options.search_program,
                                              msa_program=options.msa_program,
                                              tree_program=options.tree_program,
                                              scorematrix=options.score_matrix,
                                              gapopen=options.gapopen,
                                              gapextend=options.gapextend,
                                              extended_filename=options.extended_filename)

        elif options.q_start_from_fasta:
            self.create_output_dir_from_start_new(fasta_dir, base_dir, 
                                              user_name=options.name,
                                              search_program=options.search_program,
                                              scorematrix=options.score_matrix,
                                              msa_program=options.msa_program,
                                              tree_program=options.tree_program,
                                              gapopen=options.gapopen,
                                              gapextend=options.gapextend,
                                              extended_filename=options.extended_filename)
            
        elif options.q_start_from_blast: 
            wd1 = previous_files_locator.get_start_from_blast()
            self.start_from_orthogroups_or_sequence_search(wd1, 
                                                      base_dir,
                                                      user_name=options.name,
                                                      search_program=options.search_program,
                                                      msa_program=options.msa_program,
                                                      tree_program=options.tree_program,
                                                      scorematrix=options.score_matrix,
                                                      gapopen=options.gapopen,
                                                      gapextend=options.gapextend,
                                                      extended_filename=options.extended_filename,
                                                      use_base=options.use_base)  
            
        elif options.q_start_from_groups:
            wd1, clusters_filename_pairs = previous_files_locator.get_start_from_ogs()
            self.start_from_orthogroups_or_sequence_search(wd1, 
                                                      base_dir,
                                                      clusters_filename_pairs, 
                                                      user_name=options.name,
                                                      search_program=options.search_program,
                                                      msa_program=options.msa_program,
                                                      tree_program=options.tree_program,
                                                      scorematrix=options.score_matrix,
                                                      gapopen=options.gapopen,
                                                      gapextend=options.gapextend,
                                                      extended_filename=options.extended_filename,
                                                      use_base=options.use_base)
        elif options.q_start_from_trees or options.q_start_from_species_trees:
            wd1, clusters_filename_pairs, wd_trees, species_tree_fn, unrootedspecies_tree_fn = previous_files_locator.get_start_from_trees()
            if options.species_tree_fn != None:
                q_is_user_species_tree = True
                species_tree_fn = options.species_tree_fn
            elif species_tree_fn != None:
                q_is_user_species_tree = False
                if options.q_start_from_trees:
                    species_tree_fn = unrootedspecies_tree_fn
            else:
                if options.q_start_from_trees:
                    print("ERROR: Could not find gene tree")
                else:
                    print("ERROR: Could not find species tree")
                util.fail()
            self.start_from_trees(wd1, 
                                wd_trees,
                                base_dir, 
                                clusters_filename_pairs,
                                species_tree_fn, 
                                q_is_user_species_tree,
                                q_start_from_species_trees=options.q_start_from_species_trees,
                                user_name=options.name,
                                search_program=options.search_program,
                                msa_program=options.msa_program,
                                tree_program=options.tree_program,
                                scorematrix=options.score_matrix,
                                gapopen=options.gapopen,
                                gapextend=options.gapextend,
                                extended_filename=options.extended_filename,
                                use_base=options.use_base)
            
        # elif :
            # wd1, clustersFilename_pairs, wd_trees, speciesTreeFN = previous_files_locator.get_start_from_trees()
            # if options.speciesTreeFN != None:
            #     qIsUserSpeciesTree = True
            #     speciesTreeFN = options.speciesTreeFN
            # elif speciesTreeFN != None:
            #     qIsUserSpeciesTree = False
            # else:
            #     print("ERROR: Could not find species tree")
            #     util.Fail()
            # self.start_from_trees(wd1, 
            #                     wd_trees,
            #                     base_dir, 
            #                     clustersFilename_pairs,
            #                     speciesTreeFN, 
            #                     qIsUserSpeciesTree,
            #                     user_name=options.name,
            #                     search_program=options.search_program,
            #                     msa_program=options.msa_program,
            #                     tree_program=options.tree_program,
            #                     scorematrix=options.score_matrix,
            #                     gapopen=options.gapopen,
            #                     gapextend=options.gapextend,
            #                     extended_filename=options.extended_filename)
            

        if (options.q_start_from_groups or options.q_start_from_trees or options.q_start_from_species_trees) and previous_files_locator.species_ids_lines != None:
            # In only these cases, it's possible that the SpeciesIDs.txt file is out of sync and the version in the previous log should be used instead
            self.create_corrected_species_ids_file(previous_files_locator.species_ids_lines)

    def create_corrected_species_ids_file(self, species_ids_lines):
        self.species_ids_corrected = self.wd_current + "SpeciesIDs.txt"
        with open(self.species_ids_corrected, 'w') as outfile:
            outfile.write(species_ids_lines)
                                                      
    """ ========================================================================================== """
       
       
    # # RefactorDS - FileHandler
    # def SetNondefaultPickleDir(self, d):
    #     self.pickleDir = d
   
    def get_pickle_dir(self):
        if self.nondefault_pickle_dir != None: 
            d = self.pickleDir
        else:
            d = os.path.join(self.wd_current, "pickle") + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d
            
    """ Standard Methods
        ========================================================================================== """                       
    def log_species(self):
        text = "\nSpecies used: \n"
        fn = self.get_species_ids_fn()
        with open(fn, 'r') as infile:
            species_ids_lines = "".join(infile.readlines())
        text += species_ids_lines
        self.write_to_log(text + "\n")
        # the species of this analysis (SpeciesIDs.txt may later be edited, e.g.
        # to leave out species with '#', for another run): -fg etc. use these
        run_state.update(self.rd1, "species", ids_lines=species_ids_lines)
        
    """ Standard Directories
        ========================================================================================== """
    
    def get_working_directory1_read(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        return self.wd_base 
        
    def get_working_directory_write(self):
        if self.wd_current == None: raise Exception("No wd_current")
        return self.wd_current 
        
    def get_results_directory1(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1
        
    def get_orthologs_directory(self):
        """"Where the directories of species orthologs are"""
        if self.rd1 == None: raise Exception("No rd1")
        d = os.path.join(self.rd1, "Orthologues") + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    """ Orthogroups files 
        ========================================================================================== """
    
    def get_species_ids_fn(self):
        if self.species_ids_corrected != None:
            return self.species_ids_corrected
        return os.path.join(self.wd_base[0], "SpeciesIDs.txt")
        
    def get_sequence_ids_fn(self):
        # It is always in the first of the 'extension' directories (as this is the relevant one)
        return os.path.join(self.wd_base[0], "SequenceIDs.txt")
    
#    def get_species_ids_fn(self):
#        if self.species_ids_corrected != None:
#            return self.species_ids_corrected
#        if len(self.wd_base) == 0: raise Exception("No wd1")
#        for d in self.wd_base:
#            fn = d + "SpeciesIDs.txt"
#            if os.path.exists(fn): return fn
#        raise Exception(fn + " not found")
#        
#    def get_sequence_ids_fn(self):
#        if len(self.wd_base) == 0: raise Exception("No wd1")
#        for d in self.wd_base:
#            fn = d + "SequenceIDs.txt"
#            if os.path.exists(fn): return fn
#        raise Exception(fn + " not found")
        
    def get_species_seqs_dir(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        return self.wd_base 
        
    def get_species_fasta_fn(self, i_species, q_for_creation=False):
        """
        qForCreation: A path is required at which the file should be created (don't search for it)
        """
        if len(self.wd_base) == 0: raise Exception("No wd1")
        if q_for_creation:
            if self.wd_base[0][-1] != os.sep:
                self.wd_base[0] += os.sep
            return "%sSpecies%d.fa" % (self.wd_base[0], i_species)
        for d in self.wd_base:
            if d[-1] != os.sep:
                d += os.sep
            fn = "%sSpecies%d.fa" % (d, i_species)
            if os.path.exists(fn): return fn
        raise Exception(fn + " not found")

    def get_species_unassigned_fasta_fn(self, i_species, q_for_creation=False):
        """
        qForCreation: A path is required at which the file should be created (don't search for it)
        """
        if len(self.wd_base) == 0: raise Exception("No wd1")
        if q_for_creation:
            return "%sUnassigned.Species%d.fa" % (self.wd_base[0], i_species)
        for d in self.wd_base:
            fn = "%sUnassigned.Species%d.fa" % (d, i_species)
            if os.path.exists(fn): return fn
        raise Exception(fn + " not found")
        
    def get_sorted_species_fasta_files(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        fasta_filenames = []
        for d in self.wd_base:
            fasta_filenames.extend(glob.glob(os.path.join(d, "Species*.fa")))
        species_indices = []
        for f in fasta_filenames:
            start = f.rfind("Species")
            species_indices.append(int(f[start+7:-3]))
        indices, sorted_fasta = util.sort_array_pair_by_first(species_indices, fasta_filenames)
        return sorted_fasta  
        
    def get_species_database_n(self, i_species, program="Blast"):
        return "%s%sDBSpecies%d" % (self.wd_current, program, i_species)
    
    def get_blast_results_dir(self):
        return self.wd_base
        
    def get_blast_results_fn(self, i_species_search, j_species_db, q_for_creation=False, raise_exception=True):

        if len(self.wd_base) == 0: raise Exception("No wd1")
        if self.wd_base[0][-1] != os.sep:
            self.wd_base[0] += os.sep
        if q_for_creation: 
            return "%sBlast%d_%d.txt" % (self.wd_base[0], i_species_search, j_species_db)     
        for d in self.wd_base:
            if d[-1] != os.sep:
                d += os.sep
            fn = "%sBlast%d_%d.txt" % (d, i_species_search, j_species_db)
            if os.path.exists(fn) or os.path.exists(fn + ".gz"): return fn
        if raise_exception:
            raise Exception(fn + " not found")
        else:
            return "%sBlast%d_%d.txt" % (self.wd_base[0], i_species_search, j_species_db)
        
    def get_graph_filename(self, i_unassigned=None):
        if self.wd_current == None: raise Exception("No wd_current")
        identifier = self.file_identifier_string + ("" if i_unassigned is None else "_unassigned_clade_%d" % i_unassigned)
        return self.wd_current + "%s_graph.txt" % identifier
        
    def create_unused_clusters_fn(self, mcl_inflation=1.2, i_unassigned=None):
        mcl_inflation_str = "_I%0.1f" % mcl_inflation
        if self.wd_current == None: raise Exception("No wd_current")
        identifier = self.file_identifier_string + ("" if i_unassigned is None else "_unassigned_clade_%d" % i_unassigned)
        filename_suggestion = self.wd_current + "clusters_%s%s" % (identifier, mcl_inflation_str)
        self.clusters_filename, self.i_results_version = util.get_unused_filename(filename_suggestion, ".txt")
        return self.clusters_filename, self.clusters_filename + "_id_pairs.txt"
        
    # def SetClustersFN(self, pairsFN):
    #     self.clustersFilename = pairsFN[:-len("_id_pairs.txt")]
    #     log = "Orthogroups used: %s\n\n" % self.clustersFilename
    #     self.write_to_log(log)
        
    def get_clusters_fn(self, mcl_inflation=1.2, i_unassigned=None):
        mcl_inflation_str = "_I%0.1f" % mcl_inflation
        if self.clusters_filename is not None:
            return self.clusters_filename + "_id_pairs.txt"
        else:
            if self.wd_current == None: raise Exception("No wd_current")
            identifier = self.file_identifier_string + ("" if i_unassigned is None else "_unassigned_clade_%d" % i_unassigned)
            filename_suggestion = self.wd_current + "clusters_%s%s" % (identifier, mcl_inflation_str)
            self.clusters_filename, self.i_results_version = util.get_unused_filename(filename_suggestion, ".txt", i_unassigned)

            return self.clusters_filename + "_id_pairs.txt"
        
    """ Orthologues files
        ========================================================================================== """
    
    def get_results_seqs_dir_single_copy(self):
        d = self.rd1 + "Single_Copy_Orthologue_Sequences/"
        if not os.path.exists(d): os.mkdir(d)
        return d        
        
    def get_results_seqs_dir(self):
        d = self.rd1 + "Orthogroup_Sequences/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def get_seqs_id_dir(self):
        return self.wd_current + "Sequences_ids/"
    
    def get_legacy_dir(self):
        d = os.path.join(self.wd_current, "Legacy")
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def get_legacy_hog_dir(self):
        legacy_dir = os.path.join(self.get_legacy_dir(), "HOGs")
        if not os.path.exists(legacy_dir): os.mkdir(legacy_dir)
        return legacy_dir
        
    def get_results_align_dir(self):
        msa_dir = os.path.join(self.rd1, "MultipleSequenceAlignments/")
        if not os.path.exists(msa_dir): os.mkdir(msa_dir)
        return msa_dir
    
    def get_align_id_dir(self):
        return os.path.join(self.wd_current,  "Alignments_ids/")
        
    def get_results_trees_dir(self):
        # return self.rd1 + "Gene_Trees/"
        return os.path.join(self.wd_trees, "Gene_Trees/")

    def get_user_resolved_tree_dir(self):
        d = self.rd1 + "Resolved_Gene_Trees/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def get_resolved_tree_id_dir(self):
        d = self.wd_current + "Resolved_Gene_Trees_ids/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def get_results_species_tree_dir(self):
        d = self.rd1 + "Species_Tree/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    # def GetUserTreeFN(self):
    #     return self.get_user_tree_dir() + "Gene_Trees.txt"
    
    def get_user_resolved_tree_fn(self):
        return self.get_user_resolved_tree_dir() + "Resolved_Gene_Trees.txt"
    
    def get_ogs_seq_fn(self, i_og, q_results=False):
        if q_results:
            return self.rd1 + "Orthogroup_Sequences/" + (self.base_og_format % i_og) + ".fa"
        else:
            return self.wd_current + "Sequences_ids/" + (self.base_og_format % i_og) + ".fa"
            
    def get_ogs_align_fn(self, i_og, q_results=False):
        if q_results:
            return self.rd1 + "MultipleSequenceAlignments/" + (self.base_og_format % i_og) + ".fa"
        else:
            return self.wd_current + "Alignments_ids/" + (self.base_og_format % i_og) + ".fa"
            
    def get_ogs_tree_fn(self, i_og, q_results=False):
        if q_results:
            # return self.rd1 + "Gene_Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
            # return self.wd_trees + "Gene_Trees/" + (self.baseOgFormat % iOG) + "_tree.txt"
            return self.wd_trees + "Gene_Trees/" + (self.base_og_format % i_og) + ".txt"
        else:
            # return self.wd_trees + "Trees_ids/" + (self.baseOgFormat % iOG) + "_tree_id.txt"   
            return self.wd_trees + "Trees_ids/" + (self.base_og_format % i_og) + ".txt"   
        
    def get_species_tree_concat_align_fn(self, q_results=False):
        if q_results:
            return self.rd1 + "MultipleSequenceAlignments/" + "SpeciesTreeAlignment.fa"
        else:
            return self.wd_current + "Alignments_ids/SpeciesTreeAlignment.fa"  
        
    def get_species_tree_matrix_fn(self, q_put_in_working_dir = False):
        if q_put_in_working_dir:
            return self.wd_current + "SpeciesMatrix.phy"
        else:
            return self.wd_current + "Distances/SpeciesMatrix.phy"
            
    def get_species_tree_unrooted_fn(self, q_accessions=False):
        if q_accessions:
            return os.path.join(self.wd_trees, "SpeciesTree_unrooted.txt")
        else: 
            return os.path.join(self.wd_trees, "SpeciesTree_unrooted_ids.txt")
                        
    def get_species_tree_ids_rooted_fn(self):
        return self.wd_current + "SpeciesTree_rooted_ids.txt"

    def get_core_species_tree_ids_rooted_fn(self):
        return self.wd_base[1] + "SpeciesTree_rooted_ids.txt"
            
    def get_species_tree_results_fn(self, i, q_unique):
        """
        The results species tree (rooted, accessions, support values)
        i: index for species tree, starting at 0
        qUnique: bool, has a unique root been identified (as it may not be known exatly which branch the root belongs on)
        E.g. if there were just one species tree, the correct call would be get_species_tree_results_fn(0,True)
        """
        d = self.rd1 + "Species_Tree/"
        if not os.path.exists(d): os.mkdir(d)
        if q_unique:
            return d + "SpeciesTree_rooted.txt"
        else:
            if not self.multiple_rooted_species_trees_dir:
                self.multiple_rooted_species_trees_dir = d + "Potential_Rooted_Species_Trees/"
                if not os.path.exists(self.multiple_rooted_species_trees_dir): os.mkdir(self.multiple_rooted_species_trees_dir)
            return self.multiple_rooted_species_trees_dir + "SpeciesTree_rooted_at_outgroup_%d.txt" % i    
            
    def get_species_tree_results_node_labels_fn(self):
        return self.get_species_tree_results_fn(0, True)[:-4] + "_node_labels.txt"
        
    def get_ogs_dist_mat_fn(self, i_og):
        return self.wd_current + "Distances/OG%07d.phy" % i_og
        
    def get_species_dict(self):
        d = util.FullAccession(self.get_species_ids_fn()).get_id_to_name_dict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}

    def get_hierarchical_orthogroups_fn(self, sp_node_name, q_results=True, extension=""):
        if q_results and not extension:
            hog_dir = os.path.join(self.rd1, "Phylogenetic_Hierarchical_Orthogroups/")
            if not os.path.exists(hog_dir):
                os.makedirs(hog_dir, exist_ok=True)
            return hog_dir + "%s.tsv" % sp_node_name
        elif extension == ".ids":
            return self.get_working_directory_write() + "%s.tsv" % sp_node_name
        else:
            wd_hog_dir = os.path.join(self.rd1, "Phylogenetic_Hierarchical_Orthogroups/")
            if not os.path.exists(wd_hog_dir):
                os.makedirs(wd_hog_dir, exist_ok=True)
            return wd_hog_dir + "%s.tsv" % sp_node_name
     
    def get_result_hog_dir(self):
        hog_dir = os.path.join(self.rd1, "Phylogenetic_Hierarchical_Orthogroups/")
        if not os.path.exists(hog_dir):
            os.makedirs(hog_dir, exist_ok=True)
        return hog_dir
    
    def wd_hierarchical_orthogroups_fn_n0(self):
        return os.path.join(self.get_legacy_dir(), "N0.tsv")
    
    def hierarchical_orthogroups_fn_n0(self):
        hog_dir = os.path.join(self.rd1, "Phylogenetic_Hierarchical_Orthogroups/")
        if not os.path.exists(hog_dir):
            os.makedirs(hog_dir, exist_ok=True)

        return os.path.join(hog_dir, "N0.tsv")
    
    def ogs_all_id_fn(self):
        return self.wd_current + "OGsAll.tsv"
    
    def ogs_id_fn(self, fixed=False, extension=".tsv"):
        if fixed:
            return os.path.join(self.get_legacy_dir(), "Fixed_Orthogroups_id" + extension)
        else:
            return os.path.join(self.get_legacy_dir(), "Orthogroups_id" + extension)

        

    """ ========================================================================================== """
    
    # def GetHOGsTreeDir(self):
    #     hog_tree_dir = self.wd_current + "HOG_Trees"
    #     if not os.path.exists(hog_tree_dir): 
    #         os.mkdir(hog_tree_dir)
    #     return hog_tree_dir

    # def GetHOGMSADir(self):
    #     hog_msa_dir = self.wd_current + "HOG_Alignments"
    #     if not os.path.exists(hog_msa_dir): 
    #         os.mkdir(hog_msa_dir)
    #     return hog_msa_dir
    
    def get_check_point_fn(self):
        check_point_fn = os.path.join(self.wd_current, "checkpoint.txt")
        return check_point_fn


    def get_balsat_command_fn(self):
        commands_fn = os.path.join(self.wd_current, "blast_commands.txt")
        return commands_fn

    def get_ogs_tree_dir(self, q_results=False):
        if q_results:
            # return self.rd1 + "Gene_Trees/" 
            return self.wd_trees + "Gene_Trees/" 
        else:
            return self.wd_trees + "Trees_ids/"

    def get_astral_filename(self):
        return self.wd_trees + "astral_input.nk"

    def get_core_astral_filename(self):
        return self.wd_trees + "astral_input_core.nk"
            
    def get_ogs_recon_tree_dir(self, q_results=False):
        if q_results:
            d = self.wd_trees + "Resolved_Gene_Trees/" 
            # d = self.rd1 + "Resolved_Gene_Trees/"
            if not os.path.exists(d): os.mkdir(d)
            return d
        else:
            raise NotImplementedError()
    
    # def GetWDOGsReconTreeDir(self):
    #     d = self.wd_trees + "Resolved_Gene_Trees_1/" 
    #     if not os.path.exists(d): os.mkdir(d)
    #     return d

    # def GetPhyldogWorkingDirectory(self):
    #     d = self.wd_current + "phyldog/"
    #     if not os.path.exists(d): os.mkdir(d)
    #     return d
            
    # def GetPhyldogOGResultsTreeFN(self, i):
    #     return self.wd_current + "phyldog/OG%07d.ReconciledTree.txt" % i
    
    """ ========================================================================================== """
         
    def clean_working_dir2(self, working_dir=""):
        if not working_dir:
            working_dir = self.wd_current
        dirs = ['Distances/']
        for d in dirs:
            d_full = working_dir + d
            if os.path.exists(d_full): 
                try:
                    shutil.rmtree(d_full)
                except OSError:
                    time.sleep(1)
                    shutil.rmtree(d_full, True)  # shutil / NFS bug - ignore errors, it's less crucial that the files are deleted
                    
    """ ************************************************************************************************************************* """
            
# RefactorDS - FileHandler 
    """ Standard Methods ========================================================================================== """  
    def log_fail_and_exit(self, text=""):
        if text:
            run_logging.RunLogger.message(text, level="ERROR")
        if text != "": print(text)
        self.write_to_log("\nERROR: An error occurred\n" + text)
        util.fail()
             
    def write_to_log(self, text, q_with_time=False):
        prepend = ""
        if q_with_time:
            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + " : "
        with open(os.path.join(self.rd1, "Log.txt"), 'a') as outfile:
            outfile.write(prepend + text)
    
    def start_log(self, search_program=None, msa_program=None, tree_program=None, 
                 scorematrix=None, gapopen=None, gapextend=None):
        """
        Start Log.txt afresh, for a restart too: it records this run only. A
        restart reads what it needs from the previous Log.txt before this
        (state_from_log, PreviousFilesLocatorNew, when there is no
        run_state.json), and records the directories, species tree etc. again.
        """
        with open(os.path.join(self.rd1, "Log.txt"), "w"):
            pass
        self.write_to_log("Started OrthoFinder version " + __version__ + "\n", True)
        text = "Command Line: " + " ".join(sys.argv) + "\n\n"

        if search_program:
            text += f"Search program: {search_program}\n"
        if msa_program:
            text += f"MSA program: {msa_program}\n"
        if tree_program:
            text += f"Tree program: {tree_program}\n"

        if scorematrix:
            text += f"Scoring matrix: {scorematrix}\n"
        if gapextend:
            text += f"Gap extend: {gapextend}\n"
        if gapopen:
            text += f"Gap open: {gapopen}\n"

        text += "\nWorkingDirectory_Base: %s\n" % self.wd_base[0]
        # The earlier working directories the analysis also reads from (a
        # restart without run_state.json reads them from here): those of the
        # core analysis for --assign, otherwise those of the analysis that
        # species were added to (-b with -f).
        if len(self.wd_base) > 1:
            q_assign = self.q_assign or run_state.read(self.rd1).get("run", {}).get("analysis") == "assign"
            if q_assign:
                text += "ResultsDirectory_Core: %s\n" % (os.path.dirname(os.path.normpath(self.wd_base[1])) + os.sep)
            key = "WorkingDirectory_Core" if q_assign else "WorkingDirectory_Previous"
            text += "".join("%s: %s\n" % (key, d) for d in self.wd_base[1:])
        self.write_to_log(text)
        # the directories, as the restart options read them (run_state)
        run_state.record_dirs(self.rd1, working=self.wd_current, base=self.wd_base, trees=self.wd_trees)
        if self.clusters_filename != None:self.log_ogs()
    
    def log_ogs(self):
        self.write_to_log("FN_Orthogroups: %s\n" % (self.clusters_filename + "_id_pairs.txt"))
        run_state.record_file(self.rd1, "orthogroups", self.clusters_filename + "_id_pairs.txt")
    
    def log_assign_first_orthogroup_to_infer(self, i_og_restart):
        """--assign: see ASSIGN_FIRST_OG_LOG_KEY (a restart without run_state.json reads it from here, state_from_log)."""
        self.write_to_log("%s%d\n" % (ASSIGN_FIRST_OG_LOG_KEY, i_og_restart))

    def log_species_tree_used(self, species_tree_fn):
        """
        The species tree the analysis uses (-s, or for --assign the one inferred
        to place the new species): a restart without run_state.json reads it
        from here (state_from_log).
        """
        self.write_to_log("%s%s\n" % (SPECIES_TREE_LOG_KEY, species_tree_fn))

    def log_working_directory_trees(self):
        self.write_to_log("WorkingDirectory_Trees: %s\n" % self.wd_trees)
        run_state.record_dirs(self.rd1, trees=self.wd_trees)
        
    def make_results_directory2(self, tree_generation_method, stop_after="", append_name=""):
        """
        Args
        tree_method: msa, dendroblast, phyldog (determines the directory structure that will be created)
        stop_after: seqs, align
        """
        # RefactorDS - need to change where it puts things
        if self.rd1 == None: raise Exception("No rd1") 
        self.wd_trees = self.wd_current
        d = self.rd1 + "Orthologues/"
        if not os.path.exists(d): os.mkdir(d)
        if tree_generation_method == "msa":
            for i, d in enumerate([self.get_results_seqs_dir(), self.wd_current + "Sequences_ids/", self.get_results_align_dir(), self.wd_current + "Alignments_ids/", self.get_results_trees_dir(), self.wd_current + "Trees_ids/"]):
                if stop_after == "seqs" and i == 2: break 
                if stop_after == "align" and i == 4: break 
                if not os.path.exists(d): os.mkdir(d)
        elif tree_generation_method == "dendroblast":
            for i, d in enumerate([self.wd_current + "Distances/", self.get_results_trees_dir(), self.wd_current + "Trees_ids/"]):
                if not os.path.exists(d): os.mkdir(d)
    
    def get_orthogroup_results_fn_base(self):
        if self.rd1 == None: 
            raise Exception("No rd1")
        if self.i_results_version == None:
            raise Exception("Base results identifier has not been created")
        d = self.rd1 + "Orthogroups/"
        if not os.path.exists(d): os.mkdir(d)
        return d + "Orthogroups" + ("" if self.i_results_version == 0 else "_%d" % self.i_results_version)
    
    def get_orthogroup_results_dir(self):
        d = self.rd1 + "Orthogroups/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def get_ogs_stats_results_directory(self):
        d = self.rd1 + "Comparative_Genomics_Statistics/"
        if not os.path.exists(d): os.mkdir(d)
        return d
    
    def get_results_duplications_dir(self):
        d = self.rd1 + "Gene_Duplication_Events/"
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    def get_duplications_fn(self):
        # d = self.rd1 + "Gene_Duplication_Events/"
        # if not os.path.exists(d): os.mkdir(d)
        return os.path.join(self.get_results_duplications_dir(), "Duplications.tsv")
        # return d + "Duplications.tsv"
    
    def get_suspect_genes_dir(self):
        d = self.rd1 + "Phylogenetically_Misplaced_Genes/"
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    def get_putative_xenelogs_dir(self):
        d = self.rd1 + "Putative_Xenologs/"
        if not os.path.exists(d): os.mkdir(d)
        return d    

    def get_dependencies_check_dir(self):
        d = self.get_working_directory_write() + "dependencies/"
        if not os.path.exists(d):
            os.mkdir(d)
        return d

handler = _FileHandler()
                    
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """

class Unprocessable(Exception):
    pass

class PreviousFilesLocator(object):
    def __init__(self):
        self.wd_base_prev = []
        self.clusters_filename_pairs = None
        self.wd_trees = None
        self.home_for_results = None
        self.species_tree_rooted_ids_fn = None
        self.species_tree_unrooted_ids_fn = None
        self.species_ids_lines = None
        self.q_assign = False       # an --assign analysis (Log.txt lists WorkingDirectory_Core)
                
    def get_home_for_results(self):
        return self.home_for_results

    def get_start_from_blast(self):
        return self.wd_base_prev

    def get_start_from_ogs(self):
        return self.wd_base_prev, self.clusters_filename_pairs

    def get_start_from_trees(self):
        return self.wd_base_prev, self.clusters_filename_pairs, self.wd_trees, self.species_tree_unrooted_ids_fn, self.species_tree_rooted_ids_fn
        
""" ************************************************************************************************************************* """

class PreviousFilesLocatorNew(PreviousFilesLocator):
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        if not continuation_dir.endswith("/"): continuation_dir += "/"
        # self.home_for_results = continuationDir + "../"
        self.home_for_results = os.path.abspath(continuation_dir + "..") + os.sep
        if (options.q_start_from_fasta and not options.q_start_from_blast):
            # there are no files to find
            return
        if not self._is_new_dir_structure(continuation_dir): raise Unprocessable("Input directory structure is not processable as new structure")
        self._process_log(os.path.join(continuation_dir, "Log.txt"))
    
    def _is_new_dir_structure(self, input_dir):
        return os.path.exists(os.path.join(input_dir, "Log.txt"))
    
    def _process_log(self, log_fn):
        """
        Get all relevant data from log file. 
        Checks the paths saved do exist still
        Should work with relevant paths to allow directory to move
        Other methods can then check that the data required for a particular run is available
        """
        wd_base_anchor = None
        earlier_wds = []     # WorkingDirectory_Core/_Previous lines after the last WorkingDirectory_Base
        with open(log_fn, 'r') as infile:
            for line in infile:
                for key in EARLIER_WD_KEYS:
                    if line.startswith(key + ": "):
                        earlier_wds.append(line.rstrip()[len(key) + 2:])
                        self.q_assign = key == "WorkingDirectory_Core"

                if line.startswith("Species used:"):
                    self.species_ids_lines = ""
                    line = next(infile)
                    while line.rstrip() != "":
                        self.species_ids_lines += line
                        line = next(infile)
                wd_base_str = "WorkingDirectory_Base: "
                wd_trees_str = "WorkingDirectory_Trees: "
                clusters_str = "FN_Orthogroups: "
                if line.startswith(wd_base_str): 
                    wd_base_anchor = self._in_this_results_dir(line.rstrip()[len(wd_base_str):], log_fn, 1)
                    if not os.path.exists(wd_base_anchor):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(wd_base_anchor[:-1])
                        path, d_res = os.path.split(path)
                        # wd_base_anchor = os.path.split(logFN)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        wd_base_anchor = os.path.abspath(os.path.split(log_fn)[0] + ("/../%s/%s/" % (d_res, d_wd)))
                        if not os.path.exists(wd_base_anchor):
                            print("ERROR: Missing directory: %s" % wd_base_anchor)
                            util.fail()
                    earlier_wds = []
                    self.wd_base_prev = [wd_base_anchor]
                    self.wd_trees = wd_base_anchor
                if line.startswith(clusters_str): 
                    clusters_fn_full_path = line.rstrip()[len(clusters_str):]
                    self.clusters_filename_pairs = self._in_this_results_dir(clusters_fn_full_path, log_fn, 2)
                    if not os.path.exists(self.clusters_filename_pairs):
                        # try to see if it's a relative directory to current one
                        path, clusters_fn = os.path.split(self.clusters_filename_pairs)
                        path, d_wd = os.path.split(path)
                        path, d_res = os.path.split(path)
                        # self.clustersFilename_pairs = os.path.split(logFN)[0] + ("/../%s/%s/%s" %(d_res, d_wd, clusters_fn))
                        self.clusters_filename_pairs = os.path.abspath(os.path.split(log_fn)[0] + ("/../%s/%s/%s" %(d_res, d_wd, clusters_fn)))
                        if not os.path.exists(self.clusters_filename_pairs):
                            print("ERROR: Missing orthogroups file: %s or %s" % (self.clusters_filename_pairs, clusters_fn_full_path))
                            util.fail()
#                    self._get_ogs_file(wd_ogs_path)
                if line.startswith(wd_trees_str): 
                    self.wd_trees = self._in_this_results_dir(line.rstrip()[len(wd_trees_str):], log_fn, 1)
                    if not os.path.exists(self.wd_trees):
                        # try to see if it's a relative directory to current one
                        path, d_wd = os.path.split(self.wd_trees[:-1])
                        path, d_res = os.path.split(path)
                        self.wd_trees = os.path.split(log_fn)[0] + ("/../%s/%s/" % (d_res, d_wd))
                        self.wd_trees = os.path.abspath(self.wd_trees)
                        if not os.path.exists(self.wd_trees):
                            print("ERROR: Missing directory: %s" % self.wd_trees)
                            util.fail()
                    self.species_tree_rooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_rooted_ids.txt") 
                    self.species_tree_unrooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_unrooted_ids.txt")

        if wd_base_anchor is not None:
            # The earlier working directories the analysis reads from: listed in
            # Log.txt (this version), or followed from previous_wd.txt (older
            # versions, which wrote one in each working directory).
            if earlier_wds:
                for d in earlier_wds:
                    if not os.path.isdir(d):
                        print("ERROR: Missing directory: %s" % d)
                        util.fail()
                self.wd_base_prev = [wd_base_anchor] + earlier_wds
            else:
                self.wd_base_prev = self.get_wd_base_chain(wd_base_anchor)

    @staticmethod
    def _in_this_results_dir(recorded, log_fn, n_parts):
        """
        A path recorded in Log.txt, in the results directory being restarted
        if it has it there: the last n_parts of the recorded path (e.g.
        "WorkingDirectory/") under the directory of logFN. Paths are recorded
        absolute, so for a results directory that was copied or moved they
        name the original location, which must not be read or written then
        (whatever the copy is called).
        A working directory is used here only if it is a complete one (it has
        the run's SpeciesIDs.txt): older versions restarted into a new results
        directory whose WorkingDirectory held only the new files, and whose
        Log.txt names the results directory it continued from.
        Otherwise the recorded path (checked by the caller).
        """
        trailing = recorded.endswith(os.sep)
        parts = recorded.rstrip(os.sep).split(os.sep)[-n_parts:]
        local = os.path.join(os.path.dirname(os.path.abspath(log_fn)), *parts)
        if n_parts == 1 and not os.path.exists(os.path.join(local, "SpeciesIDs.txt")):
            return recorded
        if os.path.exists(local):
            return local + (os.sep if trailing else "")
        return recorded

    @staticmethod           
    def get_wd_base_chain(wd_base_anchor):
        chain = [wd_base_anchor]
        while os.path.exists(chain[-1] + "previous_wd.txt"):
            with open(chain[-1] + "previous_wd.txt", 'r') as infile:
                wd = infile.readline().rstrip()
                if not os.path.exists(wd):
                    # try to see if it's a relative directory to current one
                    path, d_wd = os.path.split(wd[:-1])
                    path, d_res = os.path.split(path)
                    wd = wd_base_anchor + ("/../../%s/%s/" % (d_res, d_wd))
                    wd = os.path.abspath(wd)
                    if wd[-1] != os.sep:
                        wd += os.sep
                chain.append(wd)
        return chain
                
            
""" ************************************************************************************************************************* """

class PreviousFilesLocatorState(PreviousFilesLocator):
    """
    The files of a results directory from its run_state.json (see run_state):
    the working directories it reads from (base), its orthogroups file, its
    trees directory. Raises Unprocessable if the directory has no state
    (e.g. made by an older version): Log.txt is then read instead.
    """
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        results_dir, state = run_state.find(continuation_dir)
        if not state.get("dirs", {}).get("base"):
            raise Unprocessable("No run_state.json in %s" % continuation_dir)
        self.results_dir = results_dir
        self.home_for_results = os.path.abspath(os.path.join(results_dir, "..")) + os.sep
        if options.q_start_from_fasta and not options.q_start_from_blast:
            return
        dirs = run_state.dirs(results_dir, state)
        for d in dirs["base"]:
            if not os.path.isdir(d):
                print("ERROR: A working directory of this analysis is missing: %s" % d)
                print("(recorded in %s)" % run_state.state_path(results_dir))
                util.fail()
        self.wd_base_prev = dirs["base"]
        self.wd_trees = dirs.get("trees") or self.wd_base_prev[0]
        self.clusters_filename_pairs = run_state.file(results_dir, state, "orthogroups")
        self.species_ids_lines = state.get("species", {}).get("ids_lines")
        self.species_tree_rooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_rooted_ids.txt")
        self.species_tree_unrooted_ids_fn = os.path.join(self.wd_trees, "SpeciesTree_unrooted_ids.txt")


class PreviousFilesLocatorOld(PreviousFilesLocator):
    def __init__(self, options, continuation_dir):
        PreviousFilesLocator.__init__(self)
        if not continuation_dir.endswith("/"): continuation_dir += "/"
        self.home_for_results = continuation_dir + "OrthoFinder/"
        if options.q_start_from_groups or options.q_start_from_trees:
            # User can specify it using clusters_id_pairs file, process this first to get the workingDirectory
            ogs_dir = os.path.abspath(continuation_dir + "../" if options.q_start_from_trees else continuation_dir)

            self.wd_base_prev, self.orthofinder_results_dir, self.clusters_filename_pairs = self._get_ogs_file(ogs_dir)
            if options.q_start_from_trees:
                self._find_from_trees(continuation_dir, options.species_tree_fn, unrooted=True)

            elif options.q_start_from_species_trees:
                self._find_from_trees(continuation_dir, options.species_tree_fn)
        elif options.q_start_from_blast:
            if self._is_working_directory(continuation_dir): 
                self.wd_base_prev = continuation_dir
            elif self._is_working_directory(os.path.join(continuation_dir, "WorkingDirectory") + os.sep):
                self.wd_base_prev = os.path.join(continuation_dir, "WorkingDirectory") + os.sep
            else:
                self.wd_base_prev = continuation_dir   # nothing much to do, set this as the one to try and fail later
        if self.wd_base_prev:
            self.wd_base_prev = [self.wd_base_prev]
        else:
            self.wd_base_prev = []
                
    def _get_ogs_file(self, user_arg):
        """returns the WorkingDirectory, ResultsDirectory and clusters_id_pairs filename"""
        q_specified_results_file = False
        if user_arg == None:
            print("ERROR: orthofinder_results_directory has not been specified")
            util.fail()
        if os.path.isfile(user_arg):
            fn = os.path.split(user_arg)[1]
            if ("clusters_OrthoFinder_" not in fn) or ("txt_id_pairs.txt" not in fn):
                print("ERROR:\n    %s\nis neither a directory or a clusters_OrthoFinder_*.txt_id_pairs.txt file." % user_arg)
                util.fail()
            q_specified_results_file = True
            # user has specified specific results file
        elif user_arg[-1] != os.path.sep: 
            user_arg += os.path.sep
        
        # find required files
        if q_specified_results_file:
            orthofinder_working_dir = os.path.split(user_arg)[0] + os.sep
            if not self._is_working_directory(orthofinder_working_dir):
                print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s" % orthofinder_working_dir)
                util.fail()
        else:
            orthofinder_working_dir = os.path.split(user_arg)[0] if q_specified_results_file else user_arg
            if not self._is_working_directory(orthofinder_working_dir):
                orthofinder_working_dir = user_arg + "WorkingDirectory" + os.sep   
                if not self._is_working_directory(orthofinder_working_dir):
                    print("ERROR: cannot find files from OrthoFinder run in directory:\n   %s\nor\n   %s\n" % (user_arg, orthofinder_working_dir))
                    util.fail()
                
        if q_specified_results_file:
            print("\nUsing orthogroups in file:\n    %s" % user_arg)
            return orthofinder_working_dir, orthofinder_working_dir, user_arg
        else:     
            # identify orthogroups file
            clusters_files = glob.glob(orthofinder_working_dir + "clusters_OrthoFinder_*.txt_id_pairs.txt")
            orthogroup_files = glob.glob(orthofinder_working_dir + "OrthologousGroups*.txt") + glob.glob(orthofinder_working_dir + "Orthogroups*.txt")
            if orthofinder_working_dir != user_arg:
                orthogroup_files += glob.glob(user_arg + "OrthologousGroups*.txt")
                orthogroup_files += glob.glob(user_arg + "Orthogroups*.txt")
            # User may have specified a WorkingDirectory and results could be in directory above
            if len(orthogroup_files) < len(clusters_files):
                # orthogroupFiles += glob.glob(userArg + ".." + os.sep + "OrthologousGroups*.txt")
                # orthogroupFiles += glob.glob(userArg + ".." + os.sep + "Orthogroups*.txt")
                og_path = os.path.abspath(user_arg + "..") 
                orthogroup_files += glob.glob(os.path.join(og_path, "OrthologousGroups*.txt"))
                orthogroup_files += glob.glob(os.path.join(og_path,"Orthogroups*.txt"))
            clusters_files = sorted(clusters_files)
            orthogroup_files = sorted(orthogroup_files)
            if len(clusters_files) > 1 or len(orthogroup_files) > 1:
                print("ERROR: Results from multiple OrthoFinder runs found\n")
                print("Tab-delimiter Orthogroups*.txt/OrthologousGroups*.txt files:")
                for fn in orthogroup_files:
                    print("    " + fn)
                print("With corresponding cluster files:")
                for fn in clusters_files:
                    print("    " + fn)
                print("\nPlease run with only one set of results in directories or specifiy the specific clusters_OrthoFinder_*.txt_id_pairs.txt file on the command line")
                util.fail()        
                
            if len(clusters_files) != 1 or len(orthogroup_files) != 1:
                print("ERROR: Results not found in <orthofinder_results_directory> or <orthofinder_results_directory>/WorkingDirectory")
                print("\nCould not find:\n    Orthogroups*.txt/OrthologousGroups*.txt\nor\n    clusters_OrthoFinder_*.txt_id_pairs.txt")
                util.fail()
                
            print("\nUsing orthogroups in file:\n    %s" % orthogroup_files[0])
            print("and corresponding clusters file:\n    %s" % clusters_files[0])
            return orthofinder_working_dir, user_arg, clusters_files[0]
            
    def _is_working_directory(self, orthofinder_working_dir):
        ok = True
        ok = ok and len(glob.glob(orthofinder_working_dir + "clusters_OrthoFinder_*.txt_id_pairs.txt")) > 0
        ok = ok and len(glob.glob(orthofinder_working_dir + "Species*.fa")) > 0
        return ok
        
    def _find_from_trees(self, orthologs_dir, user_species_tree, unrooted=False):
        """
        if userSpeciesTree == None: Use existing tree
        """
        print("\nFind from trees:")
        print((orthologs_dir, user_species_tree))
        self.wd_trees = os.path.join(orthologs_dir, "WorkingDirectory") + os.sep
        # Find species tree
        if user_species_tree == None:
            possibilities = ["SpeciesTree_ids_0_rooted.txt", "SpeciesTree_ids_1_rooted.txt", "SpeciesTree_user_ids.txt", "SpeciesTree_unrooted_0_rooted.txt", "STAG_SpeciesTree_ids_0_rooted.txt"] # etc (only need to determine if unique)
            n_trees = 0
            for p in possibilities:
                for d in [self.wd_trees, self.wd_trees + "Trees_ids/"]:
                    fn = d + p
                    if os.path.exists(fn): 
                        n_trees += 1
                        species_tree_fn = fn
            if n_trees == 0:
                print("\nERROR: There is a problem with the specified directory. The rooted species tree %s or %s is not present." % (possibilities[0], possibilities[2]))
                print("Please rectify the problem or alternatively use the -s option to specify the species tree to use.\n")
                util.fail()
            if n_trees > 1:
                print("\nERROR: There is more than one rooted species tree in the specified directory structure. Please use the -s option to specify which species tree should be used\n")
                util.fail()
            if unrooted and "unrooted" in fn:
                self.species_tree_unrooted_ids_fn = species_tree_fn
            else:
                self.species_tree_rooted_ids_fn = species_tree_fn
        else:
            if not os.path.exists(user_species_tree):
                print("\nERROR: %s does not exist\n" % user_species_tree)
                util.fail()
            if unrooted:
                self.species_tree_unrooted_ids_fn = species_tree_fn
            else:
                self.species_tree_rooted_ids_fn = species_tree_fn
            
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """

def initialise_file_handler(
        options, 
        fasta_dir=None, 
        continuation_dir=None, 
        results_dir_non_default=None, 
        pickle_dir_non_default=None,
        working_dir="",
    ):
    """
    Creates a file handler object which will determine the location of all the files:
    Results will be under the user specified directory of the default results location. Defaults:
        - New, from start: 
                FastaDir/OrthoFinder/Results_Date
            or
                resultsDir_nonDefault/Results_Date
        - New, continuation: Existing_OrthoFinder_Dir/Results_Date
        - Old, continuation: 
                ContinuationDir/OrthoFinder/Results_Date
            or 
                resultsDir_nonDefault/Results_Date
        
    
    Implementation
    1. Working out if an old directory structure is being used
    2. Construct and appropriate PreviousFilesLocator if necessary - this locates all required files
    3. Pass this to FileHandler - this creates the directory structure required for this run
    4. if error: print and exit
    5. Return FileHandler
    
    Tasks:
    - Switch this round, I can tell if it's and old or new directory right from the start - read log and check info present,
    perhaps just psss it to the new file handler and let it decide if everything is there
    """
    # 1 & 2
    # If starting from scratch, no need for a PreviousFileLocator
    if fasta_dir is not None and fasta_dir[-1] != os.sep:
        fasta_dir += os.sep
    if options.q_start_from_fasta and not options.q_start_from_blast:
        pfl = None
        base_dir = results_dir_non_default if results_dir_non_default != None else fasta_dir + "OrthoFinder/"
    else:
        try:
            try:
                # The locations recorded in the results directory (run_state.json)
                pfl = PreviousFilesLocatorState(options, continuation_dir)
            except Unprocessable:
                # Otherwise from its Log.txt (results of older versions)
                pfl = PreviousFilesLocatorNew(options, continuation_dir)
            # don't create any new directory, it already exists
            base_dir = pfl.get_home_for_results()

        except Unprocessable:
            pfl = PreviousFilesLocatorOld(options, continuation_dir)
            if not pfl.wd_base_prev:
                print("ERROR: Couldn't process as an OrthoFinder results directory: %s" % continuation_dir)
                util.fail()
            base_dir = results_dir_non_default if results_dir_non_default is not None else pfl.get_home_for_results()
    if not os.path.exists(base_dir): os.mkdir(base_dir)
    if base_dir[-1] != os.sep:
        base_dir += os.sep
    # 3 
    # RefactorDS - this might be suitable as a constructor now
    # base_dir - should now exist
    
    """The previous file locator should decide where the output directory should be rooted? 
    Rules:
    - If starting from Fasta then Fasta/OrthoFinder/Results_Date
          - Or, SpecifiedDirectory/Results_Date if user specified
    - If starting from a previous new-structure directory (TopLevel/Results_X) then TopLevel/Results_Date
    - If starting from a previous old-structure directory then, as high up as we can go and still be in the directory structure:
        - Fasta/Results_OldDate/OrthoFinder/Results_Date
    """
    if not working_dir:
        handler.create_output_directories(options, pfl, base_dir, fasta_dir)
    else:
        if working_dir[-1] != os.sep:
            working_dir += os.sep
        handler.wd_base = [working_dir]
        handler.wd_current = working_dir
        handler.rd1 = os.path.dirname(working_dir)
        handler.wd_trees = working_dir
        
