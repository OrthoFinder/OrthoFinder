#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
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
# david_emms@hotmail.com

"""
Created on Thu Sep 25 13:15:22 2014

@author: david
"""
import os
import glob
import copy
try:
    from rich import print
except ImportError:
    ...

from ..utils import util, parallel_task_manager
from ..utils import file_handler
from . import trim
from ..species_tree import astral as astral_pro
from ..utils import printing
from ..phylo import taxa
from ..utils import command_runner
from ..species_tree.concatenated_alignment import create_concatenated_alignment, determine_orthogroups_for_species_tree


class FastaWriter(object):
    def __init__(self, fasta_file_dir_list, species_to_use):
        """
        The list of previous directories could include species that are no longer used,
        but not ones that are conflicting
        """
        self.seq_lists = dict()
        q_first = True
        accession = ""
        sequence = ""
        required_files = set(["Species%d.fa" % i for i in species_to_use])
        for d in fasta_file_dir_list:
            for fn in glob.glob(d + "Species*.fa"):
                if os.path.basename(fn) not in required_files: continue
                with open(fn, 'r') as fasta_file:
                    for line in fasta_file:
                        if line[0] == ">":
                            # deal with old sequence
                            if not q_first:
                                self.seq_lists[accession] = sequence
                                sequence = ""
                            q_first = False
                            # get id for new sequence
                            accession = line[1:].rstrip()
                        else:
                            sequence += line
                    self.seq_lists[accession] = sequence
    
    def write_seqs_to_fasta(self, seqs, out_filename):
        with open(out_filename, 'w') as out_file:
            for seq in self.sort_seqs([s.to_string() for s in seqs]):
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % seq)
                    out_file.write(self.seq_lists[seq])
                else:
                    print(("ERROR: %s not found" % seq))
                                
    def write_seqs_to_fasta_with_new_accessions(self, seqs, out_filename, id_dict):
        try:
            seqs = [s.to_string() for s in seqs]
        except:
            pass
        with open(out_filename, 'w') as out_file:
            for seq in self.sort_seqs(seqs):
                if seq in self.seq_lists:
                    out_file.write(">%s\n" % id_dict[seq])
                    out_file.write(self.seq_lists[seq])
                    
    def sort_seqs(self, seqs):
        return sorted(seqs, key=lambda x: list(map(int, x.split("_"))))

def write_test_file(d):
    test_fn = d + "SimpleTest.fa"
    with open(test_fn, 'w') as outfile:
        outfile.write(">a\nA\n>b\nA")
    return test_fn
 
""" 
-----------------------------------------------------------------------------
                    Identify orthogroups for species tree            
-----------------------------------------------------------------------------
"""

    # N = m.shape[1]
    # f = 1./N
    # fractions = []
    # nOrtho = []
    # for n in range(N):
    #     F = 1.-n*f
    #     fractions.append(F)
    #     nOrtho.append(len(single_copy_with_probability_test(F-1e-5, m)))
    # nOrtho = list(map(float, nOrtho))
    # return fractions, nOrtho
    
def trim_fn(fn, msa_min_seq):
    trim.main(fn, fn, msa_min_seq, 0.1, 500, 0.75, False)

""" 
-----------------------------------------------------------------------------
                             TreesForOrthogroups            
-----------------------------------------------------------------------------
"""    
     
class TreesForOrthogroups(object):
    def __init__(self, program_caller, msa_program, tree_program):
        self.program_caller = program_caller
        self.msa_program = msa_program
        self.tree_program = tree_program
    
    def get_fasta_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_seq_fn(i_og, q_results)
            
    def get_alignment_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_align_fn(i_og, q_results)
            
    def get_tree_filename(self, i_og, q_results=False):
        return file_handler.handler.get_ogs_tree_fn(i_og, q_results)
        
    # def write_fasta_files(self, fastaWriter, ogs, idDict, qBoth):
    #     # The results ones are now written by default after orthogroups, check they're not already there
    #     for iOg, og in enumerate(ogs):
    #         if not os.path.exists(self.get_fasta_filename(iOg, True)):
    #             fastaWriter.write_seqs_to_fasta_with_new_accessions(og, self.get_fasta_filename(iOg, True), idDict)
    #     if qBoth: 
    #         for iOg, og in enumerate(ogs):
    #             if not os.path.exists(self.get_fasta_filename(iOg)):
    #                 fastaWriter.write_seqs_to_fasta(og, self.get_fasta_filename(iOg))

    def write_fasta_files(self, fasta_writer, ogs, id_dict, q_id=True, q_results=False):
        # The results ones are now written by default after orthogroups, check they're not already there
        if q_results:
            for i_og, og in enumerate(ogs):
                if not os.path.exists(self.get_fasta_filename(i_og, True)):
                    fasta_writer.write_seqs_to_fasta_with_new_accessions(og, self.get_fasta_filename(i_og, True), id_dict)
        if q_id: 
            for i_og, og in enumerate(ogs):
                if not os.path.exists(self.get_fasta_filename(i_og)):
                    fasta_writer.write_seqs_to_fasta(og, self.get_fasta_filename(i_og))

    def get_alignment_commands_and_new_filenames(self, ogs, i_og_restart=0, method_threads=None):
        iogs_align = [i for i, og in enumerate(ogs) if len(og) >= 2 and i >= i_og_restart]
        infn_list = [self.get_fasta_filename(i) for i in iogs_align]
        outfn_list = [self.get_alignment_filename(i) for i in iogs_align]
        id_list = ["OG%07d" % i for i in iogs_align]
        n_seqs = [len(ogs[i]) for i in iogs_align]
        return self.program_caller.get_msa_commands(self.msa_program, infn_list, outfn_list, id_list, n_seqs, method_threads=method_threads), iogs_align
        
    def get_tree_commands(self, alignments_for_tree, iogs_align, ogs, method_threads=None):
        iogs_tree = [i for i in iogs_align if len(ogs[i]) >= 3]
        alignments_for_tree = [alignments_for_tree[i] for i, iog in enumerate(iogs_align) if len(ogs[iog]) >= 3]  # Remove those to align with tree filenames
        outfn_list = [self.get_tree_filename(i) for i in iogs_tree]
        id_list = ["OG%07d" % i for i in iogs_tree]
        n_seqs = [len(ogs[i]) for i in iogs_tree]
        return self.program_caller.get_tree_commands(self.tree_program, alignments_for_tree, outfn_list, id_list, n_seqs, method_threads=method_threads), iogs_tree
     
    def rename_alignment_taxa(self, ids_align_fns, acc_align_fns, ids_dict):
        for i, (align_fn, out_align_fn) in enumerate(zip(ids_align_fns, acc_align_fns)):
            if not os.path.exists(align_fn): continue
            with open(align_fn, 'r') as infile, open(out_align_fn, 'w') as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(">" + ids_dict[line[1:].rstrip()] + "\n")
                    else:
                        outfile.write(line)

          
    def do_trees(
            self, 
            og_set, 
            id_dict, 
            species_id_dict, 
            species_to_use, 
            n_processes, 
            q_stop_after_seqs, 
            q_stop_after_alignments, 
            q_do_species_tree, 
            q_trim, 
            i_og_restart=0, 
            cmd_order="descending",
            method_threads=None,
            method_threads_large=None,
            method_threads_small=None, 
            threshold=None,
            fix_files=True,
            astral=False,
            dynamic_threads=False,
            n_skip=50,
            msa_min_seq=10,
        ):

        print_on_error = True
        id_dict.update(species_id_dict) # same code will then also convert concatenated alignment for species tree
        # 0       
        results_dirs_full_path = [
            file_handler.handler.get_results_seqs_dir(), 
            file_handler.handler.get_results_align_dir(), 
            file_handler.handler.get_results_trees_dir()
        ]
        
        # 1.
        fasta_writer = FastaWriter(file_handler.handler.get_species_seqs_dir(), species_to_use)
        ogs = og_set.get_ogs_all()
        
        ## -------- Writing to WD/Sequence_ids is not is mcl_orthogroups/infer_orthogroups.py ---------------
        # self.write_fasta_files(fastaWriter, ogs, idDict, True)

        if q_stop_after_seqs: return results_dirs_full_path

        # 3
        # Get OGs to use for species tree
        if q_do_species_tree:
            og_matrix, iogs4 = og_set.orthogroup_matrix()
            i_ogs_for_species_tree, f_single_copy = determine_orthogroups_for_species_tree(og_matrix, iogs4)
            concatenated_algn_fn = file_handler.handler.get_species_tree_concat_align_fn()
        else:
            i_ogs_for_species_tree = []
        orig_align_commands_and_filenames, orig_iogs_align = \
            self.get_alignment_commands_and_new_filenames(
                ogs, 
                i_og_restart, 
                method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
            )

        align_commands_and_filenames = copy.deepcopy(orig_align_commands_and_filenames)
        iogs_align = copy.deepcopy(orig_iogs_align)
        aln_tasksize = None
        if cmd_order == "ascending":
            # aln_tasksize = [len(og) for og in aln_filtered_ogs[::-1]]
            orig_align_commands_and_filenames = orig_align_commands_and_filenames[::-1]
            orig_iogs_align = orig_iogs_align[::-1]
            aln_tasksize = [
                os.stat(item).st_size
                for fn, _ in orig_align_commands_and_filenames
                for item in fn.split() if os.path.isfile(item) and os.path.exists(item)
            ]
        if q_stop_after_alignments:
            printing.print_underline("Inferring multiple sequence alignments")
            command_runner.run_parallel_commands_and_move_results_file(
                n_processes, 
                orig_align_commands_and_filenames, 
                False, 
                tasksize=aln_tasksize,
                method_threads=method_threads,
                method_threads_large=method_threads_large,
                method_threads_small=method_threads_small, 
                threshold=threshold,
                cmd_order=cmd_order,
                q_trim=q_trim,
                q_print_on_error=print_on_error,
                dynamic_threads=dynamic_threads
            )
            if q_do_species_tree:
                if not astral:
                    create_concatenated_alignment(
                        i_ogs_for_species_tree, 
                        ogs, 
                        self.get_alignment_filename, 
                        concatenated_algn_fn, 
                        f_single_copy
                    )
                    # write OGs used to file
                    d_species_tree = os.path.split(file_handler.handler.get_species_tree_results_fn(0, True))[0] + "/"
                    with open(d_species_tree + "Orthogroups_for_concatenated_alignment.txt", 'w') as outfile:
                        for iog in i_ogs_for_species_tree: outfile.write("OG%07d\n" % iog)
            
            # ------------------ this section is not needed at this stage for the new procedure -------------
            # ids -> accessions
            alignment_files_to_use = [self.get_alignment_filename(i) for i, in orig_iogs_align]
            accession_alignment_fns = [self.get_alignment_filename(i, True) for i, in orig_iogs_align]
            if q_do_species_tree: 
                alignment_files_to_use.append(concatenated_algn_fn)
                accession_alignment_fns.append(file_handler.handler.get_species_tree_concat_align_fn(True))
            self.rename_alignment_taxa(alignment_files_to_use, accession_alignment_fns, id_dict)
            
            return results_dirs_full_path[:2]
        
        # Otherwise, alignments and trees
        # Strategy is
        # 1. Do alignments (and trees) require for species tree
        # 2. Create concatenated alignment
        # 3. Create second list of commands [speciestree] + [remaining alignments and trees]
        alignment_files_to_use = [self.get_alignment_filename(i) for i in iogs_align]

        orig_tree_commands_and_filenames, orig_iogs_tree = \
            self.get_tree_commands(
                alignment_files_to_use, 
                iogs_align, 
                ogs,
                method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
        )

        tree_commands_and_filenames = copy.deepcopy(orig_tree_commands_and_filenames)
        iogs_tree = copy.deepcopy(orig_iogs_tree)
        commands_and_filenames = []
        tree_tasksize = None
        if q_do_species_tree:
            if not astral:
                print(("Species tree: Using %d orthogroups with minimum of %0.1f%% of species having single-copy genes in any orthogroup" % (len(i_ogs_for_species_tree), 100.*f_single_copy)))
            printing.print_underline("Inferring multiple sequence alignments for species tree") 
            # util.print_time("This may take some time...")
            # Do required alignments and trees
            species_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
            orig_commands_and_filenames = []
            if q_trim:
                for i in i_ogs_for_species_tree:
                    orig_commands_and_filenames.append([align_commands_and_filenames[i], 
                                                (trim_fn, (alignment_files_to_use[i], msa_min_seq)), 
                                                tree_commands_and_filenames[i]])
            else:
                for i in i_ogs_for_species_tree:
                    orig_commands_and_filenames.append([align_commands_and_filenames[i],
                                                   tree_commands_and_filenames[i]])
            commands_and_filenames = copy.deepcopy(orig_commands_and_filenames)        
            if cmd_order == "ascending":
                # tree_tasksize = [len(og) for og in tree_filtered_ogs[::-1]]
                orig_commands_and_filenames = orig_commands_and_filenames[::-1]
                tree_tasksize = [
                        os.stat(item).st_size
                        for cmd in orig_commands_and_filenames
                        for item in cmd[0][0].split() if os.path.isfile(item) and os.path.exists(item)
                    ]
                
            command_runner.run_parallel_commands_and_move_results_file(
                n_processes, 
                orig_commands_and_filenames, 
                True, 
                tasksize=tree_tasksize,
                method_threads=method_threads,
                method_threads_large=method_threads_large,
                method_threads_small=method_threads_small, 
                threshold=threshold,
                cmd_order=cmd_order,
                q_trim=q_trim,
                q_print_on_error=print_on_error,
                dynamic_threads=dynamic_threads
            )
            create_concatenated_alignment(i_ogs_for_species_tree, ogs, self.get_alignment_filename, concatenated_algn_fn, f_single_copy)
            # write OGs used to file
            d_species_tree = os.path.split(file_handler.handler.get_species_tree_results_fn(0, True))[0] + "/"
            with open(d_species_tree + "Orthogroups_for_concatenated_alignment.txt", 'w') as outfile:
                for iog in i_ogs_for_species_tree: outfile.write("OG%07d\n" % iog)
            # Add species tree to list of commands to run
            if not astral:
                util.print_time("Scheduling species tree inference from concatenated alignment")
                commands_and_filenames = [
                    self.program_caller.get_tree_commands(
                            self.tree_program, 
                            [concatenated_algn_fn], 
                            [species_tree_fn_ids], 
                            ["SpeciesTree"],
                            method_threads=None  # filled in when run, which reserves that many cores (command_runner.fill_method_threads)
                        )
                    ]
            printing.print_underline("Inferring remaining multiple sequence alignments and gene trees") 
        else:
            printing.print_underline("Inferring multiple sequence alignments and gene trees")
        # util.print_time("This may take some time...")
        # Now continue as before
        i_ogs_for_species_tree = set(i_ogs_for_species_tree)
        iog_to_align_index = {iog: index for index, iog in enumerate(iogs_align)}
        for i, iog in enumerate(iogs_tree):
            if iog in i_ogs_for_species_tree: 
                continue
            if q_trim:
                commands_and_filenames.append([align_commands_and_filenames[iog_to_align_index[iog]],
                                              (trim_fn, (alignment_files_to_use[iog_to_align_index[iog]], msa_min_seq)),
                                              tree_commands_and_filenames[i]])
            else:
                commands_and_filenames.append([align_commands_and_filenames[iog_to_align_index[iog]],
                                              tree_commands_and_filenames[i]])
        for iog in set(iogs_align).difference(iogs_tree):
            if iog in i_ogs_for_species_tree: 
                continue
            commands_and_filenames.append([align_commands_and_filenames[iog_to_align_index[iog]]])
        
        if cmd_order == "ascending":
            commands_and_filenames = commands_and_filenames[::-1]
            tree_tasksize = [
                        os.stat(item).st_size
                        for cmd in orig_commands_and_filenames
                        for item in cmd[0][0].split() if os.path.isfile(item) and os.path.exists(item)
                    ]
        command_runner.run_parallel_commands_and_move_results_file(n_processes, 
                                                 commands_and_filenames, 
                                                 True, 
                                                 tasksize=tree_tasksize,
                                                 method_threads=method_threads,
                                                 method_threads_large=method_threads_large,
                                                 method_threads_small=method_threads_small, 
                                                 threshold=threshold,
                                                 cmd_order=cmd_order,
                                                 q_trim=q_trim,
                                                 q_print_on_error=print_on_error,
                                                 dynamic_threads=dynamic_threads
                                                 )
        
        # # Convert ids to accessions for MSA
        # accessionAlignmentFNs = [self.get_alignment_filename(i, True) for i in iogs_align]
        # # Add concatenated Alignment
        # if qDoSpeciesTree:
        #     alignmentFilesToUse.append(concatenated_algn_fn)
        #     accessionAlignmentFNs.append(files.FileHandler.get_species_tree_concat_align_fn(True))

        #     qHaveSupport = taxa.have_support_values(speciesTreeFN_ids)
        #     if os.path.exists(speciesTreeFN_ids):
        #         taxa.rename_tree_taxa(speciesTreeFN_ids, files.FileHandler.get_species_tree_unrooted_fn(True), idDict, qSupport=qHaveSupport, qFixNegatives=True)
        #     else:
        #         text = "ERROR: Species tree inference failed"
        #         files.FileHandler.log_fail_and_exit(text)


        # Add concatenated Alignment
        if q_do_species_tree:
            if astral:
                printing.print_underline("Inferring unrooted species tree using Astral-Pro.") 
                astral_fn = file_handler.handler.get_core_astral_filename()
                astral_pro.create_input_file(file_handler.handler.get_ogs_tree_dir(), astral_fn, n_skip=n_skip)
                species_tree_fn_ids = file_handler.handler.get_species_tree_unrooted_fn()
                parallel_task_manager.run_command(
                    astral_pro.get_astral_command(
                        astral_fn, species_tree_fn_ids, n_processes
                    ),
                    # ASTRAL-Pro writes its progress log to stderr: only show it if it fails.
                    q_print_stderr=False,
                    raise_on_error=True,
                )

            q_have_support = taxa.have_support_values(species_tree_fn_ids)
            if os.path.exists(species_tree_fn_ids):
                taxa.rename_tree_taxa(species_tree_fn_ids, file_handler.handler.get_species_tree_unrooted_fn(True), id_dict, q_support=q_have_support, q_fix_negatives=True)
                util.print_time("Done species tree inference")
            else:
                text = "ERROR: Species tree inference failed"
                file_handler.handler.log_fail_and_exit(text)


         # ------------------ this section is not needed at this stage for the new procedure -------------
        if not fix_files:
            # Convert ids to accessions for MSA
            accession_alignment_fns = [self.get_alignment_filename(i, True) for i in iogs_align]
            # Add concatenated Alignment
            if q_do_species_tree:
                alignment_files_to_use.append(concatenated_algn_fn)
                accession_alignment_fns.append(file_handler.handler.get_species_tree_concat_align_fn(True))

            self.rename_alignment_taxa(alignment_files_to_use, accession_alignment_fns, id_dict)
        
        return results_dirs_full_path[:2]
