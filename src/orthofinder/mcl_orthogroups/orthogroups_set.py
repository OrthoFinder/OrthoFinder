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
# david_emms@hotmail.comhor: david

import numpy as np

try:
    from rich import print
except ImportError:
    ...
from . import mcl as MCL
from ..utils import util, file_handler


class Seq(object):
    def __init__(self, seq_input):
        """Constructor takes sequence in any format and returns generators the
        Seq object accordingly. If performance is really important then can write
        individual an @classmethod to do that without the checks"""
        if type(seq_input) is str:
            a, b = seq_input.split("_")
            self.i_sp = int(a)
            self.i_seq = int(b)
        elif len(seq_input) == 2:
            if seq_input[0] is str:
                self.i_sp, self.i_seq = list(map(int, seq_input))
            else:
                self.i_sp = seq_input[0]
                self.i_seq = seq_input[1]
        else:
            raise NotImplementedError

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.to_string()

    def to_string(self):
        return "%d_%d" % (self.i_sp, self.i_seq)


# ==============================================================================================================================


class OrthogroupSet(object):
    def __init__(
        self,
        min_seq,
        orthofinder_working_dir_list,
        species_to_use,
        n_sp_all,
        q_add_species_to_ids,
        tree_program="fasttree",
        id_extractor=util.FirstWordExtractor,
        mcl_inflation=1.2,
    ):

        self.species_ids_ex = util.FullAccession(file_handler.handler.get_species_ids_fn())
        self._spec_seq_ids = None
        self._extractor = id_extractor
        self.seq_ids_ex = None
        self.ogs_all = None
        self.ogs_single = None
        self.i_ogs4 = None
        self.species_to_use = species_to_use  # list of ints
        self.seqs_info = util.get_seqs_info(
            orthofinder_working_dir_list, self.species_to_use, n_sp_all
        )
        self.id_to_og = None
        self.q_add_species_to_ids = q_add_species_to_ids
        self.cached_seq_ids_dict = None
        self.min_seq = min_seq
        self.tree_program = tree_program
        self.mcl_inflation = mcl_inflation

    def sequence_dict(self):
        """returns Dict[str, str]"""
        if self.cached_seq_ids_dict is not None:
            return self.cached_seq_ids_dict
        if self.seq_ids_ex == None:
            try:
                self.seq_ids_ex = self._extractor(file_handler.handler.get_sequence_ids_fn())
            except RuntimeError as error:
                print(str(error))
                if str(error).startswith("ERROR"):
                    file_handler.handler.log_fail_and_exit()
                else:
                    print(
                        "Tried to use only the first part of the accession in order to list the sequences in each orthogroup"
                    )
                    print(
                        "more concisely but these were not unique. The full accession line will be used instead.\n"
                    )
                    self.seq_ids_ex = util.FullAccession(
                        file_handler.handler.get_sequence_ids_fn()
                    )
        self.cached_seq_ids_dict = self.seq_ids_ex.get_id_to_name_dict()
        return self.cached_seq_ids_dict

    def species_dict(self):
        """returns Dict[str, str]"""
        d = self.species_ids_ex.get_id_to_name_dict()
        return {k: v.rsplit(".", 1)[0] for k, v in d.items()}

    def spec_seq_dict(self):
        """returns Dict[str, str]"""
        if self._spec_seq_ids != None:
            return self._spec_seq_ids
        seqs = self.sequence_dict()
        seqs = {
            k: v for k, v in seqs.items() if int(k.split("_")[0]) in self.species_to_use
        }
        if not self.q_add_species_to_ids:
            self._spec_seq_ids = seqs
            return seqs
        specs = self.species_dict()
        specs_ed = {k: v.replace(".", "_").replace(" ", "_") for k, v in specs.items()}
        self._spec_seq_ids = {
            seq_id: specs_ed[seq_id.split("_")[0]] + "_" + name
            for seq_id, name in seqs.items()
        }
        return self._spec_seq_ids

    def get_i_ogs4(self):
        if self.i_ogs4 is None:
            ogs = self.get_ogs_all()
            self.i_ogs4 = [i for i, og in enumerate(ogs) if len(og) >= self.min_seq]
        return self.i_ogs4

    def get_ogs_all(self, i_unassigned=None):
        if self.ogs_all is None:
            ogs = MCL.get_predicted_ogs(file_handler.handler.get_clusters_fn(self.mcl_inflation, i_unassigned))
            if self.tree_program in ["raxml", "raxml-ng"]:
                self.ogs_all = [
                    [Seq(g) for g in og] for og in ogs if len(og) >= self.min_seq
                ]
            else:
                self.ogs_all = [[Seq(g) for g in og] for og in ogs]
            # self.ogs_all = sorted(self.ogs_all, key=len, reverse=True)
        return self.ogs_all

    # def get_ogs_single(self):
    #     if self.ogs_single is None:
    #         ogs = MCL.get_predicted_ogs(files.FileHandler.get_clusters_fn())
    #         # self.ogs_single = [og for og in ogs if len(og) == 1]
    #         all_assigned = {g for og in ogs if len(og) != 1 for g in og}
    #         unassigned = set(self.all_used_sequence_ids()).difference(all_assigned)
    #         self.ogs_single = [{g,} for g in unassigned]

    #     return self.ogs_single

    def get_ogs_single(self):
        if self.ogs_single is None:
            ogs = MCL.get_predicted_ogs(file_handler.handler.get_clusters_fn())
            self.ogs_single = [og for og in ogs if len(og) == 1]
        return self.ogs_single

    # def OGs4AssumeOrdered(self):
    #     ogs_all = self.get_ogs_all()
    #     iogs4 = self.get_i_ogs4()
    #     return [ogs_all[i] for i in iogs4]

    def orthogroup_matrix(self):
        """qReduce give a matrix with only as many columns as species for cases when
        clustering has been performed on a subset of species"""
        ogs = self.get_ogs_all()
        iogs4 = self.get_i_ogs4()
        ogs = [ogs[i] for i in iogs4]
        i_species = sorted(set([gene.i_sp for og in ogs for gene in og]))
        species_index_dict = {i_sp: i_col for i_col, i_sp in enumerate(i_species)}
        n_species = len(i_species)
        n_groups = len(ogs)
        # (i, j)-th entry of ogMatrix gives the number of genes from i in orthologous group j
        og_matrix = np.zeros((n_groups, n_species))
        for i_og, og in enumerate(ogs):
            for gene in og:
                og_matrix[i_og, species_index_dict[gene.i_sp]] += 1
        return og_matrix, iogs4

    def id_to_og_dict(self):
        if self.id_to_og != None:
            return self.id_to_og
        # Maybe shouldn't include unclustered genes:
        self.id_to_og = {
            g.to_string(): iog for iog, og in enumerate(self.get_ogs_all()) for g in og
        }
        return self.id_to_og

    def all_used_sequence_ids(self):
        ids_dict = self.sequence_dict()
        species_to_use_strings = list(map(str, self.species_to_use))
        all_ids = [
            s for s in ids_dict.keys() if s.split("_")[0] in species_to_use_strings
        ]
        return all_ids
