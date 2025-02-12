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
from ..tools import mcl as MCL
from ..utils import util, files


class Seq(object):
    def __init__(self, seqInput):
        """ Constructor takes sequence in any format and returns generators the 
        Seq object accordingly. If performance is really important then can write 
        individual an @classmethod to do that without the checks"""
        if type(seqInput) is str:
            a,b = seqInput.split("_")
            self.iSp = int(a)
            self.iSeq = int(b)
        elif len(seqInput) == 2:
            if seqInput[0] is str:
                self.iSp, self.iSeq = list(map(int, seqInput))
            else:
                self.iSp= seqInput[0]
                self.iSeq = seqInput[1]
        else:
            raise NotImplementedError
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not self.__eq__(other)         
        
    def __repr__(self):
        return self.ToString()
    
    def ToString(self):
        return "%d_%d" % (self.iSp, self.iSeq)

# ==============================================================================================================================
        
class OrthoGroupsSet(object):
    def __init__(
            self, 
            orthofinderWorkingDir_list, 
            speciesToUse, 
            nSpAll, 
            qAddSpeciesToIDs, 
            idExtractor = util.FirstWordExtractor
        ):
        
        self.speciesIDsEx = util.FullAccession(files.FileHandler.GetSpeciesIDsFN())
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.seqIDsEx = None
        self.ogs_all = None
        self.ogs_single = None
        self.iOgs4 = None
        self.speciesToUse = speciesToUse     # list of ints
        self.seqsInfo = util.GetSeqsInfo(orthofinderWorkingDir_list, self.speciesToUse, nSpAll)
        self.id_to_og = None
        self.qAddSpeciesToIDs = qAddSpeciesToIDs
        self.cached_seq_ids_dict = None

    def SequenceDict(self):
        """returns Dict[str, str]"""
        if self.cached_seq_ids_dict is not None:
            return self.cached_seq_ids_dict
        if self.seqIDsEx == None:
            try:
                self.seqIDsEx = self._extractor(files.FileHandler.GetSequenceIDsFN())
            except RuntimeError as error:
                print(str(error))
                if str(error).startswith("ERROR"): 
                    files.FileHandler.LogFailAndExit()
                else:
                    print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup")
                    print("more concisely but these were not unique. The full accession line will be used instead.\n")
                    self.seqIDsEx = util.FullAccession(files.FileHandler.GetSequenceIDsFN())
        self.cached_seq_ids_dict = self.seqIDsEx.GetIDToNameDict()
        return self.cached_seq_ids_dict
        
    def SpeciesDict(self):
        """returns Dict[str, str]"""
        d = self.speciesIDsEx.GetIDToNameDict()
        return {k: v.rsplit(".", 1)[0] for k, v in d.items()}
        
    def Spec_SeqDict(self):
        """returns Dict[str, str]"""
        if self._Spec_SeqIDs != None:
            return self._Spec_SeqIDs
        seqs = self.SequenceDict()
        seqs = {k:v for k,v in seqs.items() if int(k.split("_")[0]) in self.speciesToUse}
        if not self.qAddSpeciesToIDs:
            self._Spec_SeqIDs = seqs
            return seqs
        specs = self.SpeciesDict()
        specs_ed = {k:v.replace(".", "_").replace(" ", "_") for k,v in specs.items()}
        self._Spec_SeqIDs = {seqID:specs_ed[seqID.split("_")[0]] + "_" + name for seqID, name in seqs.items()}
        return self._Spec_SeqIDs

    def Get_iOGs4(self):
        if self.iOgs4 is None:
            ogs = self.OGsAll()
            self.iOgs4 = [i for i, og in enumerate(ogs) if len(og) >= 4]
        return self.iOgs4

    def OGsAll(self):
        if self.ogs_all is None:
            ogs = MCL.GetPredictedOGs(files.FileHandler.GetClustersFN())
            self.ogs_all = [[Seq(g) for g in og] for og in ogs]
        return self.ogs_all
    
    def OGsSingle(self):
        if self.ogs_single is None:
            ogs = MCL.GetPredictedOGs(files.FileHandler.GetClustersFN())
            self.ogs_single = [og for og in ogs if len(og) == 1]
        return self.ogs_single

    # def OGs4AssumeOrdered(self):
    #     ogs_all = self.OGsAll()
    #     iogs4 = self.Get_iOGs4()
    #     return [ogs_all[i] for i in iogs4]

    def OrthogroupMatrix(self):
        """ qReduce give a matrix with only as many columns as species for cases when
        clustering has been performed on a subset of species"""
        ogs = self.OGsAll()
        iogs4 = self.Get_iOGs4()
        ogs = [ogs[i] for i in iogs4]
        iSpecies = sorted(set([gene.iSp for og in ogs for gene in og]))
        speciesIndexDict = {iSp:iCol for iCol, iSp in enumerate(iSpecies)}
        nSpecies = len(iSpecies)
        nGroups = len(ogs)
        # (i, j)-th entry of ogMatrix gives the number of genes from i in orthologous group j
        ogMatrix = np.zeros((nGroups, nSpecies)) 
        for i_og, og in enumerate(ogs):
            for gene in og:
                ogMatrix[i_og, speciesIndexDict[gene.iSp]] += 1
        return ogMatrix, iogs4
        
    def ID_to_OG_Dict(self):
        if self.id_to_og != None:
            return self.id_to_og
        # Maybe shouldn't include unclustered genes:
        self.id_to_og = {g.ToString():iog for iog, og in enumerate(self.OGsAll()) for g in og}
        return self.id_to_og

    def AllUsedSequenceIDs(self):
        ids_dict = self.SequenceDict()
        species_to_use_strings = list(map(str, self.speciesToUse))
        all_ids = [s for s in ids_dict.keys() if s.split("_")[0] in species_to_use_strings]
        return all_ids
