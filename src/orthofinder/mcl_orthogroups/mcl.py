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
from __future__ import absolute_import
import os
import sys
from typing import List, Set
try:
    from rich import print
except ImportError:
    ...

from ..utils import parallel_task_manager, util, file_handler, file_io

from collections import defaultdict, Counter
from xml.etree.ElementTree import SubElement    # Y
from .. import __version__
from ..utils import gene_ids

def get_predicted_ogs(clusters_filename):
    """
    Returns:
        ogs: List[Set[str]] - List of sets of genes ids (as strings)
    """
    predicted_ogs = []
    n_ogs_string = ""
    q_contains_profiles = False
    with open(clusters_filename, 'r') as cluster_file:
        header = True
        og = set()
        for line in cluster_file:
            if header:
                if line.count("begin"):
                    header = False
            else:
                if line.find(")") != -1:
                    break
                if line[-2] == "$":
                    line = line[:-3]
                if line[0] == " ":
                    # continuation of group
                    x = line.split()
                    y = [x_ for x_ in x if not x_.startswith('Prof')]
                    og = og.union(y)
                else:
                    # new OG
                    if len(og) != 0:
                        predicted_ogs.append(og)
                    n_ogs_string, line = line.split(" ", 1)
                    x = line.split()
                    y = [x_ for x_ in x if not x_.startswith('Prof')]
                    if len(x) != len(y):
                        q_contains_profiles = True
                    og = set(y)
        if len(og) > 0:
            predicted_ogs.append(og)
    return predicted_ogs
    
def get_id_pair(species_starting_indices, single_id, species_to_use):   
    for i, starting_index in enumerate(species_starting_indices):
        if starting_index > single_id:
            return "%d_%d" % (species_to_use[i-1], single_id - species_starting_indices[i-1])
    return "%d_%d" % (species_to_use[-1], single_id - species_starting_indices[len(species_starting_indices)-1]) 

def convert_single_ids_to_id_pair(seqs_info, clusters_filename, new_filename, q_unassigned=False):
    with open(clusters_filename, 'r') as cluster_file, open(new_filename, "w") as output:
        header = True
        for line in cluster_file:
            append_dollar = False
            initial_text = ""
            ids_string = ""
            ids = []
            if header:
                output.write(line)
                if line.count("begin"):
                    header = False
            else:
                if line.find(")") != -1:
                    output.write(line)
                    break
                if q_unassigned and not line.startswith(" ") and len(line.split()) == 3:
                    # This is a single gene, not necessarily from this MCL clustering but here implicitly
                    continue
                if line[-2] == "$":
                    line = line[:-3]
                    append_dollar = True
                if line[0] != " ":
                    initial_text, line = line.split(None, 1)
                # continuation of group
                ids = line.split()
                for id in ids:
                    ids_string += get_id_pair(seqs_info.seq_starting_indices, int(id), seqs_info.species_to_use) + " "
                output.write(initial_text + "      " + ids_string)
                if append_dollar:
                    output.write("$\n")
                else:
                    output.write("\n")


def write_updated_clusters_file(ogs: List[Set[str]], clusters_filename_pairs):
    with open(clusters_filename_pairs, 'w') as outfile:
        outfile.write("(mclmatrix")
        outfile.write("begin\n")
        for iog, og in enumerate(ogs):
            outfile.write("%d      %s $\n" % (iog, " ".join(sorted(og))))
        outfile.write(")\n")


def id_full_dict(ids_filenames, func=gene_ids.FirstWordExtractor):

    full_dict = dict()
    for ids_filename in ids_filenames:
        id_extract = func(ids_filename)
        id_dict = id_extract.get_id_to_name_dict()
        full_dict.update(id_dict)

    return full_dict

class MCL:
    @staticmethod
    def create_ogs(predicted_ogs, output_filename, orthogroups_id_fn, id_dict):
        """Orthogroups.txt-like files: output_filename with the gene names, orthogroups_id_fn
        with the sequence IDs (either may be None: not written)."""
        if output_filename is not None:
            with open(output_filename, 'w') as output_file:
                for i_og, og in enumerate(predicted_ogs):
                    output_file.write("OG%07d: " % i_og)
                    output_file.write(" ".join(sorted([id_dict[seq] for seq in og])))
                    output_file.write("\n")
        if orthogroups_id_fn is not None:
            with open(orthogroups_id_fn, "w") as outputid_file:
                for i_og, og in enumerate(predicted_ogs):
                    outputid_file.write("OG%07d: " % i_og)
                    outputid_file.write(" ".join(sorted([seq for seq in og])))
                    outputid_file.write("\n")

    # @staticmethod

    @staticmethod
    def run_mcl(graph_filename, clusters_filename, n_processes, inflation):
        command = " ".join(["mcl", graph_filename, "-I", str(inflation), "-o", clusters_filename, "-te", str(n_processes), "-V", "all"])
        parallel_task_manager.run_command(command, q_print_on_error=True, raise_on_error=True)
        util.print_time("Ran MCL")

    @staticmethod
    def write_orthogroup_files(
            ogs, 
            ids_filenames, 
            results_base_filename, 
            clusters_filename_pairs,
            write_ids=True,
        ):
        """
        <results_base_filename>.txt (None: not written) and, with write_ids, the
        orthogroups as sequence IDs in the working directory. Returns {sequence ID: gene name}.
        """
        output_fn = None if results_base_filename is None else results_base_filename + ".txt"
        orthogroups_id_fn = file_handler.handler.ogs_id_fn(extension=".txt") if write_ids else None
        try:
            # fullDict = dict()
            # for idsFilename in idsFilenames:
            #     idExtract = gene_ids.FirstWordExtractor(idsFilename)
            #     idDict = idExtract.get_id_to_name_dict()
            #     fullDict.update(idDict)
            full_dict = id_full_dict(ids_filenames, func=gene_ids.FirstWordExtractor)
            MCL.create_ogs(ogs, output_fn, orthogroups_id_fn, full_dict)
        except KeyError as e:
            sys.stderr.write("ERROR: Sequence ID not found in %s\n" % ids_filenames)
            sys.stderr.write(str(e) + "\n")
            file_handler.handler.log_fail_and_exit(("ERROR: Sequence ID not found in %s\n" % ids_filenames) + str(e) + "\n")
        except RuntimeError as error:
            print(str(error))
            if str(error).startswith("ERROR"):
                err_text = "ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. If %s was prepared manually then please check the IDs are correct. " % (ids_filenames, clusters_filename_pairs, ids_filenames)
                file_handler.handler.log_fail_and_exit(err_text)
            else:
                print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")
                try:
                    # fullDict = dict()
                    # for idsFilename in idsFilenames:
                    #     idExtract = gene_ids.FullAccession(idsFilename)
                    #     idDict = idExtract.get_id_to_name_dict()
                    #     fullDict.update(idDict)
                    full_dict = id_full_dict(ids_filenames, func=gene_ids.FullAccession)
                    MCL.create_ogs(ogs, output_fn, orthogroups_id_fn, full_dict)
                except:
                    err_text = "ERROR: %s contains a duplicate ID. The IDs for the orthogroups in %s will not be replaced with the sequence accessions. This is probably because the same accession was used more than once in your input FASTA files. However, if %s was prepared manually then you may need to check that file instead." % (ids_filenames, clusters_filename_pairs, ids_filenames)
                    file_handler.handler.log_fail_and_exit(err_text)
        return full_dict


    @staticmethod
    def create_orthogroup_table(
        ogs,
        id_to_name_dict,
        species_names_dict,
        species_to_use,
        results_base_filename,
        write_ids=True,
    ):
        """<results_base_filename>.tsv, .GeneCount.tsv, _UnassignedGenes.tsv and, with
        write_ids, the orthogroups as sequence IDs in the working directory."""

        n_species = len(species_names_dict)
        
        ogs_ids = [[seq for seq in og] for og in ogs]
        ogs_names = [[id_to_name_dict[seq] for seq in og] for og in ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in ogs]

        # write out
        orthogroups_id_fn = file_handler.handler.ogs_id_fn() if write_ids else os.devnull
        output_filename = results_base_filename + ".tsv"
        output_filename_counts = results_base_filename + ".GeneCount.tsv"
        single_gene_filename = results_base_filename + "_UnassignedGenes.tsv"
        with open(output_filename, util.csv_write_mode) as output_file, \
            open(orthogroups_id_fn, util.csv_write_mode) as ogidfile, \
            open(single_gene_filename, util.csv_write_mode) as single_gene_file, \
            open(output_filename_counts, util.csv_write_mode) as out_file_counts:
            
            ogid_filewriter = file_io.writer(ogidfile)
            file_writer = file_io.writer(output_file)
            
            file_writer_counts = file_io.writer(out_file_counts)
            single_gene_writer = file_io.writer(single_gene_file)
            for writer in [ogid_filewriter, file_writer, single_gene_writer]:
                row = ["Orthogroup"] + [species_names_dict[index] for index in species_to_use]
                writer.writerow(row)
            file_writer_counts.writerow(row + ['Total'])

            for i_og, (og, og_names, og_ids) in enumerate(zip(ogs_ints, ogs_names, ogs_ids)):
                og_dict = defaultdict(list)
                og_id_dict = defaultdict(list)
                row = ["OG%07d" % i_og]
                id_row = ["OG%07d" % i_og]
                this_output_writer = file_writer
                # separate it into sequences from each species
                if len(og) == 1:
                    row.extend(['' for x in range(n_species)])
                    row[species_to_use.index(og[0][0]) + 1] = og_names[0]
                    this_output_writer = single_gene_writer
                else:
                    for (i_species, i_sequence), name, ids in zip(og, og_names, og_ids):
                        og_dict[species_to_use.index(i_species)].append(name)
                        og_id_dict[species_to_use.index(i_species)].append(ids)
                    for i_species in range(n_species):
                        row.append(", ".join(sorted(og_dict[i_species])))
                        id_row.append(", ".join(sorted(og_id_dict[i_species])))
                    ogid_filewriter.writerow(id_row)
                    counts = Counter([i_species for i_species, _ in og])
                    counts_row = [counts[i_species] for i_species in species_to_use]
                    file_writer_counts.writerow(row[:1] + counts_row + [sum(counts_row)])
                this_output_writer.writerow(row)
