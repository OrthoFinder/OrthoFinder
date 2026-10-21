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
from . import gene_ids, file_io
from . import paths

# Log.txt line of the species tree an analysis uses (FileHandler.log_species_tree_used).
SPECIES_TREE_LOG_KEY = "FN_SpeciesTree: "
# --assign: the first orthogroup whose gene tree was not inferred to place the
# new species (the trees of those before it are reused by a restart).
ASSIGN_FIRST_OG_LOG_KEY = "Assign_FirstOrthogroupToInfer: "

# Log.txt lines of the programs a run used (start_log), as the run_state settings.
PROGRAM_LOG_KEYS = (("Search program: ", "search_program"), ("MSA program: ", "msa_program"),
                    ("Tree program: ", "tree_program"))


# Log.txt lines of the earlier working directories an analysis reads from,
# after its own (WorkingDirectory_Base): those of the core analysis (--assign),
# or of the analysis it added species to (-b with -f).
EARLIER_WD_KEYS = ("WorkingDirectory_Core", "WorkingDirectory_Previous")


def rename_earlier_result_names(results_dir):
    """
    The results of OrthoFinder before 3.2 used "Orthologue" in some names
    (Orthologues/, OrthologuesStats_*.tsv, ...). A run that continues in such
    a results directory (-b, -fg, -fgt, -fst) renames them first, so it
    updates them in place, as before, and the directory has one set of names.
    Nothing is renamed if the new name is already there. Returns
    [(old path, new path)] of what was renamed.
    """
    renamed = []

    def rename(old, new):
        if os.path.exists(old) and not os.path.exists(new):
            os.rename(old, new)
            renamed.append((old, new))

    orthologs_dir = os.path.join(results_dir, "Orthologs")
    rename(os.path.join(results_dir, "Orthologues"), orthologs_dir)
    if os.path.isdir(orthologs_dir):
        for name in os.listdir(orthologs_dir):  # -pof: Orthologues_<species>/
            if name.startswith("Orthologues_"):
                rename(os.path.join(orthologs_dir, name), os.path.join(orthologs_dir, "Orthologs_" + name[len("Orthologues_"):]))
    rename(os.path.join(results_dir, "Single_Copy_Orthologue_Sequences"),
           os.path.join(results_dir, "Single_Copy_Ortholog_Sequences"))
    for subdir, old_part, new_part in (("Comparative_Genomics_Statistics", "OrthologuesStats_", "OrthologsStats_"),
                                       ("Orthogroups", "_SingleCopyOrthologues", "_SingleCopyOrthologs")):
        d = os.path.join(results_dir, subdir)
        if os.path.isdir(d):
            for name in os.listdir(d):
                if old_part in name:
                    rename(os.path.join(d, name), os.path.join(d, name.replace(old_part, new_part)))
    return renamed


def _dir_attribute(name):
    """An attribute holding a directory path, always kept ending in the separator (paths.as_dir)."""
    key = "_" + name
    def get(self):
        return getattr(self, key)
    def set(self, value):
        setattr(self, key, None if value is None else paths.as_dir(value))
    return property(get, set)


class _FileHandler(object):
    # Directories: file names are appended to them, so they always end in the separator
    wd_current = _dir_attribute("wd_current")
    wd_trees = _dir_attribute("wd_trees")
    rd1 = _dir_attribute("rd1")

    @property
    def wd_base(self):
        return self._wd_base

    @wd_base.setter
    def wd_base(self, dirs):
        self._wd_base = [paths.as_dir(d) for d in dirs]

    def __init__(self):
        self.reset() 

    def reset(self):
        self.base_og_format = "OG%07d"
        self.wd_base = []               # Base: blast, species & sequence IDs, species fasta files - should not request this and then write here
        self.wd_current = None          # Location to write out any new files
        self.wd_trees = None            # Location of working dir containing tree files
        self.keep_legacy = False        # --keep-legacy: the legacy files are written (WorkingDirectory/Legacy)
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
            self.rd1 = paths.create_new_working_directory(base + "Results_", 
                                                      search_program=search_program,
                                                      msa_program=msa_program,
                                                      tree_program=tree_program,
                                                      scorematrix=scorematrix,
                                                      gapopen=gapopen,
                                                      gapextend=gapextend,
                                                      extended_filename=extended_filename)
        else:
            self.rd1 = paths.create_new_working_directory(base + "Results_" + user_name, 
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
            shutil.copy(os.path.join(old_wd_base_list[0], "SpeciesIDs.txt"), self.wd_current + "SpeciesIDs.txt")
            shutil.copy(os.path.join(old_wd_base_list[0], "SequenceIDs.txt"), self.wd_current + "SequenceIDs.txt")
            # The earlier working directories are recorded in run_state.json and
            # Log.txt (start_log); older versions wrote previous_wd.txt instead.
            self.wd_base = self.wd_base + list(old_wd_base_list)
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
            self.rd1 = os.path.dirname(os.path.normpath(self.wd_current))
            os.makedirs(self.wd_current, exist_ok=True)
            rename_earlier_result_names(self.rd1)  # continuing in an existing results directory
        else:    
            if clusters_filename_pairs != None: self.clusters_filename = clusters_filename_pairs[:-len("_id_pairs.txt")]
            if user_name == None:
                self.rd1 = paths.create_new_working_directory(base + "Results_",
                                                        search_program=search_program,
                                                        msa_program=msa_program,
                                                        tree_program=tree_program, 
                                                        scorematrix=scorematrix,
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            else:
                self.rd1 = paths.create_new_working_directory(base + "Results_" + user_name,
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
            self.rd1 = os.path.dirname(os.path.normpath(self.wd_current))
            os.makedirs(self.wd_current, exist_ok=True)
            rename_earlier_result_names(self.rd1)  # continuing in an existing results directory
        else:    
            if user_name == None:
                self.rd1 = paths.create_new_working_directory(base + "Results_",
                                                        search_program=search_program, 
                                                        msa_program=msa_program,
                                                        tree_program=tree_program,
                                                        scorematrix=scorematrix,
                                                        gapopen=gapopen,
                                                        gapextend=gapextend,
                                                        extended_filename=extended_filename)
            else:
                self.rd1 = paths.create_new_working_directory(base + "Results_" + user_name, 
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
        with open(self.species_ids_corrected, 'w', encoding="utf-8", newline="\n") as outfile:
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
        with file_io.open_text_any(fn) as infile:
            species_ids_lines = "".join(infile)
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
        d = os.path.join(self.rd1, "Orthologs") + os.sep
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
            return "%sSpecies%d.fa" % (self.wd_base[0], i_species)
        for d in self.wd_base:
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
        if q_for_creation: 
            return "%sBlast%d_%d.txt" % (self.wd_base[0], i_species_search, j_species_db)     
        for d in self.wd_base:
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
        self.clusters_filename, self.i_results_version = paths.get_unused_filename(filename_suggestion, ".txt")
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
            self.clusters_filename, self.i_results_version = paths.get_unused_filename(filename_suggestion, ".txt", i_unassigned)

            return self.clusters_filename + "_id_pairs.txt"
        
    """ Orthologs files
        ========================================================================================== """
    
    def get_results_seqs_dir_single_copy(self):
        d = self.rd1 + "Single_Copy_Ortholog_Sequences/"
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
        
    CONCAT_ORTHOGROUPS_FN = "Orthogroups_for_concatenated_alignment.txt"

    def get_concat_orthogroups_fn(self):
        """The orthogroups of the concatenated alignment for the species tree (working directory)."""
        return os.path.join(self.wd_current, "Alignments_ids", self.CONCAT_ORTHOGROUPS_FN)

    def copy_concat_orthogroups_to_results(self):
        """
        With the results species tree: Species_Tree/Orthogroups_for_concatenated_alignment.txt,
        from the working directory of the run that made the alignment (also for a restart).
        Not from the base working directories: those of the core run for --assign, whose
        species tree is not from a concatenated alignment.
        """
        for d in [self.wd_current, self.wd_trees]:
            if d is None:
                continue
            fn = os.path.join(d, "Alignments_ids", self.CONCAT_ORTHOGROUPS_FN)
            if os.path.isfile(fn):
                out_dir = os.path.join(self.rd1, "Species_Tree")
                os.makedirs(out_dir, exist_ok=True)
                shutil.copy2(fn, os.path.join(out_dir, self.CONCAT_ORTHOGROUPS_FN))
                return

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
        return gene_ids.species_names(self.get_species_ids_fn())

    def get_hierarchical_orthogroups_fn(self, sp_node_name, q_results=True, extension=""):
        """
        A HOG file of the ortholog inference, each written once. The HOGs of the
        two passes are named by the orthogroups they are of:
        q_results (the second pass, of the final orthogroups): N<i>.tsv in the
            results (not N0: the orthogroups are Orthogroups.tsv) and N0.ids.tsv in
            the working directory;
        else (the first pass, of the MCL orthogroups): N0_mcl.ids.tsv, from which the
            orthogroups are updated (n0_mcl_ids_fn), and the HOGs with gene names,
            N<i>_mcl.tsv, which are only written as legacy files (--keep-legacy:
            Legacy/HOGs).
        """
        if extension == ".ids":   # sp_node_name "N0.ids"
            if q_results:
                return self.get_working_directory_write() + "%s.tsv" % sp_node_name
            return self.n0_mcl_ids_fn()
        if q_results:
            return self.get_result_hog_dir() + "%s.tsv" % sp_node_name
        return os.path.join(self.get_legacy_hog_dir(), "%s_mcl.tsv" % sp_node_name)

    def n0_mcl_ids_fn(self):
        """
        N0_mcl.ids.tsv: the HOGs of N0 of the MCL orthogroups (the first pass), from
        which the orthogroups are updated. A legacy file (--keep-legacy: Legacy/HOGs),
        else in the working directory, removed at the end of the run.
        """
        if self.keep_legacy:
            return os.path.join(self.get_legacy_hog_dir(), "N0_mcl.ids.tsv")
        return self.get_working_directory_write() + "N0_mcl.ids.tsv"

    def get_result_hog_dir(self):
        hog_dir = os.path.join(self.rd1, "Phylogenetic_Hierarchical_Orthogroups/")
        if not os.path.exists(hog_dir):
            os.makedirs(hog_dir, exist_ok=True)
        return hog_dir

    def ogs_all_id_fn(self):
        return self.wd_current + "OGsAll.tsv"

    def ogs_id_fn(self, fixed=False, extension=".tsv"):
        """The orthogroups as sequence IDs, before (MCL) and after (fixed) the update: legacy files (--keep-legacy) only."""
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
        with open(os.path.join(self.rd1, "Log.txt"), 'a', encoding="utf-8", newline="\n") as outfile:
            outfile.write(prepend + text)
    
    def start_log(self, search_program=None, msa_program=None, tree_program=None, 
                 scorematrix=None, gapopen=None, gapextend=None):
        """
        Start Log.txt afresh, for a restart too: it records this run only. A
        restart reads what it needs from the previous Log.txt before this
        (state_from_log, PreviousFilesLocatorNew, when there is no
        run_state.json), and records the directories, species tree etc. again.
        """
        with open(os.path.join(self.rd1, "Log.txt"), "w", encoding="utf-8", newline="\n"):
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
        d = self.rd1 + "Orthologs/"
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
    
    def get_mcl_orthogroup_results_fn_base(self):
        """Orthogroups/MCL_Orthogroups: the MCL orthogroups (run.finish.write_mcl_orthogroups)."""
        d, name = os.path.split(self.get_orthogroup_results_fn_base())
        return os.path.join(d, "MCL_" + name)

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

""" ************************************************************************************************************************* """

""" ************************************************************************************************************************* """

""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """
""" ************************************************************************************************************************* """

