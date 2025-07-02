import sys
import os 
import csv 
from operator import itemgetter
from collections import defaultdict, Counter

from ..tools import trees_msa
from ..tools import mcl as MCL
from ..utils import util, files

import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom
from .. import __version__


def post_hogs_processing(
        all_seq_ids,
        speciesInfoObj,
        seqsInfo,
        speciesNamesDict,
        options,
        speciesXML,
        q_incremental=False,
    ):
    """
    Write OGs & statistics to results files, write Fasta files.
    Args:
        q_incremental - These are not the final orthogroups, don't write results
    """
    new_ogs, name_dictionary = \
        update_ogs(files.FileHandler.HierarchicalOrthogroupsFNN0())
    resultsBaseFilename = files.FileHandler.GetOrthogroupResultsFNBase()
    # util.PrintUnderline("Writing orthogroups to file")
    all_assigned = set([g for og in new_ogs for g in og])
    unassigned = set(all_seq_ids).difference(all_assigned)
    single_ogs_list = [{g,} for g in unassigned]
    new_ogs.extend(single_ogs_list)

    with open(files.FileHandler.OGsAllIDFN(), "w") as outfile:
        for og in new_ogs:
            outfile.write(", ".join(og) + "\n")

    idsDict = MCL.WriteOrthogroupFiles(
        new_ogs,
        [files.FileHandler.GetSequenceIDsFN()],
        resultsBaseFilename,
    )

    if not q_incremental:
        MCL.CreateOrthogroupTable(
            new_ogs,
            idsDict,
            speciesNamesDict,
            speciesInfoObj.speciesToUse,
            resultsBaseFilename,
        )

    # Write Orthogroup FASTA files
    ogSet = OrthoGroupsSet(
        options.min_seq,
        files.FileHandler.GetWorkingDirectory1_Read(), 
        speciesInfoObj.speciesToUse,
        speciesInfoObj.nSpAll,
        options.qAddSpeciesToIDs,
        options.tree_program,
        idExtractor=util.FirstWordExtractor,
    )

    ## ------------------ Fix Orthogroup_Sequences and Sequences_ids --------------------
    treeGen = trees_msa.TreesForOrthogroups(None, None, None)
    fastaWriter = trees_msa.FastaWriter(
        files.FileHandler.GetSpeciesSeqsDir(), speciesInfoObj.speciesToUse
    )

    d_seqs = files.FileHandler.GetResultsSeqsDir()
    if not os.path.exists(d_seqs):
        os.mkdir(d_seqs)

    # Update Orthogroup_Sequneces
    treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idsDict, qID=False, qResults=True)

    idDict = ogSet.Spec_SeqDict()
    idDict.update(ogSet.SpeciesDict()) # same code will then also convert concatenated alignment for species tree

    # Update Orthogroup_Sequneces and Sequences_ids
    # treeGen.WriteFastaFiles(fastaWriter, ogSet.OGssAll(), idDict, True)
    # treeGen.WriteFastaFiles(fastaWriter, ogSet.OGsAll(), idDict, False) # Set to False, only update the Orthogroup_Sequneces 
    
    if not q_incremental:
        # stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
        if options.speciesXMLInfoFN:
            MCL.WriteOrthoXML(
                speciesXML,
                new_ogs,
                seqsInfo.nSeqsPerSpecies,
                idsDict,
                resultsBaseFilename + ".orthoxml",
                speciesInfoObj.speciesToUse,
            )
        # print("")
        # util.PrintTime("Done orthogroups")
        files.FileHandler.LogOGs()

    return ogSet, idDict, name_dictionary

def update_ogs(input_path):
    sorted_matrix = read_hogs_to_matrix(input_path)   
    name_dictionary = {}     
    new_og_list = []
    # For each line in sorted HOG order replace HOG name with index OG name (based on length of enumerate so HOG.N0 + 0*x + number)
    for pos, line in enumerate(sorted_matrix):
        new_og_name = "OG%07d" % pos
        key = line[2]
        if key in name_dictionary:
            additional = [new_og_name] + [line[1]] + [line[3]]           
            new_value = name_dictionary[key]
            new_value.append(additional)
            name_dictionary.update({key: new_value})
        else:
            name_dictionary[key] = [[new_og_name] + [line[1]] + [line[3]]]

        new_og_set = set(", ".join(line[4:]).replace("\n", "").split(", "))
        new_og_list.append({gene for gene in new_og_set if len(gene) != 0})
    return  new_og_list, name_dictionary

def read_hogs_to_matrix(input_path):
    #holds lines to write to new output file
    matrix = []
    with open(input_path) as input_file:
        for i, line in enumerate(input_file):
            line_split = line.strip().split("\t")
            if i == 0:
                # species_names = line_split[3:]
                continue
            # Count number of genes in each line as a count.
            count = len(list(filter(None, (", ".join(line_split[3:]).split(", ")))))
            # Add line with gene count to matrix variable.
            matrix.append([count] + line_split)

    sorted_matrix = sorted(matrix, key=itemgetter(0), reverse=True)
    return sorted_matrix  

def GetSingleID(speciesStartingIndices, seq, speciesToUse): 
    a, b = seq.split("_")
    iSpecies = int(a)
    iSeq = int(b)
    offset = speciesStartingIndices[speciesToUse.index(iSpecies)]
    return iSeq + offset  


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
            min_seq,
            orthofinderWorkingDir_list, 
            speciesToUse, 
            nSpAll, 
            qAddSpeciesToIDs, 
            tree_prgram = "fasttree",
            idExtractor = util.FirstWordExtractor
        ):
        
        self.speciesIDsEx = util.FullAccession(files.FileHandler.GetSpeciesIDsFN())
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.seqIDsEx = None
        self.ogs_all = None
        self.iOgs4 = None
        self.speciesToUse = speciesToUse     # list of ints
        self.seqsInfo = util.GetSeqsInfo(orthofinderWorkingDir_list, self.speciesToUse, nSpAll)
        self.id_to_og = None
        self.qAddSpeciesToIDs = qAddSpeciesToIDs
        self.cached_seq_ids_dict = None
        self.min_seq = min_seq
        self.tree_program = tree_prgram

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
            self.iOgs4 = [i for i, og in enumerate(ogs) if len(og) >= self.min_seq]
        return self.iOgs4

    def OGsAll(self):
        if self.ogs_all is None:
            with open(files.FileHandler.OGsAllIDFN()) as infile:
                ogs = [og.strip().split(", ") for og in infile]
            if self.tree_program == "raxml":
                self.ogs_all = [[Seq(g) for g in og]  for og in ogs if len(og) >= self.min_seq]
            else:
                self.ogs_all = [[Seq(g) for g in og] for og in ogs]

            # self.ogs_all = sorted(self.ogs_all, key=len, reverse=True)
        return self.ogs_all
    
    def AllOGs(self):
        with open(files.FileHandler.OGsAllIDFN()) as infile:
                ogs = [og.strip().split(", ") for og in infile]
        if self.tree_program == "raxml":
            all_ogs = [[g for g in og]  for og in ogs if len(og) >= self.min_seq]
        else:
            all_ogs = [[g for g in og] for og in ogs]
        return all_ogs
        
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


class MCL:
    @staticmethod
    def CreateOGs(predictedOGs, outputFilename, idDict):
        with open(outputFilename, 'w') as outputFile:
            for iOg, og in enumerate(predictedOGs):
                outputFile.write("OG%07d: " % iOg)
                accessions = sorted([idDict[seq] for seq in og if idDict.get(seq) is not None])
                outputFile.write(" ".join(accessions))
                outputFile.write("\n")

    @staticmethod
    def prettify(elem):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    @staticmethod
    def WriteOrthoXML(speciesInfo, predictedOGs, nSequencesDict, idDict, orthoxmlFilename, speciesToUse):
        """ speciesInfo: ordered array for which each element has
            fastaFilename, speciesName, NCBITaxID, sourceDatabaseName, databaseVersionFastaFile
        """
        # Write OrthoXML file
        root = ET.Element("orthoXML")
        root.set('xsi:schemaLocation', "http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd")
        root.set('originVersion', __version__)
        root.set('origin', 'OrthoFinder')
        root.set('version', "0.3")
        root.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")
        #notes = SubElement(root, 'notes')

        # Species: details of source of genomes and sequences they contain
        speciesStartingIndices = []
        iGene_all = 0
        for iPos, thisSpeciesInfo in enumerate(speciesInfo):
            iSpecies = speciesToUse[iPos]
            nSeqs = nSequencesDict[iSpecies]
            speciesNode = SubElement(root, 'species')
            speciesNode.set('NCBITaxId', thisSpeciesInfo[2])           # required
            speciesNode.set('name', thisSpeciesInfo[1])                # required
            speciesDatabaseNode = SubElement(speciesNode, "database")
            speciesDatabaseNode.set('name', thisSpeciesInfo[3])            # required
            speciesDatabaseNode.set('version', thisSpeciesInfo[4])         # required
    #            speciesDatabaseNode.set('geneLink', "")        # skip
    #            speciesDatabaseNode.set('protLink', "")        # skip
    #            speciesDatabaseNode.set('transcriptLink', "")  # skip
            allGenesNode = SubElement(speciesDatabaseNode, "genes")
            speciesStartingIndices.append(iGene_all)
            for iGene_species in range(nSeqs):
                geneNode = SubElement(allGenesNode, 'gene')
                geneNode.set("geneId", idDict["%d_%d" % (iSpecies , iGene_species)])
                geneNode.set('id', str(iGene_all))       # required
    #                geneNode.set("protID", "")  # skip
                iGene_all += 1

        # Scores tag - unused
    #            scoresNode = SubElement(root, 'scores')        # skip

        # Orthogroups
        allGroupsNode = SubElement(root, 'groups')
        for iOg, og in enumerate(predictedOGs):
            groupNode = SubElement(allGroupsNode, 'orthologGroup')
            groupNode.set('id', str(iOg))
    #                groupScoreNode = SubElement(groupNode, 'score')    # skip
    #                groupScoreNode.set('id', "")                       # skip
    #                groupScoreNode.set('value', "")                    # skip
    #                SubElement(groupNode, 'property')                  # skip
            for seq in og:
                geneNode = SubElement(groupNode, 'geneRef')
                geneNode.set('id', str(GetSingleID(speciesStartingIndices, seq, speciesToUse)))
    #                    SubElement(geneNode, 'score')                  # skip
        with open(orthoxmlFilename, 'w') as orthoxmlFile:
    #            ET.ElementTree(root).write(orthoxmlFile)
            orthoxmlFile.write(MCL.prettify(root))
        print("Orthogroups have been written to orthoxml file:\n   %s" % orthoxmlFilename)

    @staticmethod
    def WriteOrthogroupFiles(
            ogs, 
            idsFilenames, 
            resultsBaseFilename, 
        ):
        outputFN = resultsBaseFilename + ".txt"
        try:
            fullDict = dict()
            for idsFilename in idsFilenames:
                idExtract = util.FirstWordExtractor(idsFilename)
                idDict = idExtract.GetIDToNameDict()
                fullDict.update(idDict)
            MCL.CreateOGs(ogs, outputFN, fullDict)
        except KeyError as e:
            sys.stderr.write("ERROR: Sequence ID not found in %s\n" % idsFilename)
            sys.stderr.write(str(e) + "\n")
            files.FileHandler.LogFailAndExit(("ERROR: Sequence ID not found in %s\n" % idsFilename) + str(e) + "\n")
        except RuntimeError as error:
            print(str(error))
            if str(error).startswith("ERROR"):
                err_text = "ERROR: %s contains a duplicate ID. " % (idsFilename)
                files.FileHandler.LogFailAndExit(err_text)
            else:
                print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")
                try:
                    fullDict = dict()
                    for idsFilename in idsFilenames:
                        idExtract = util.FullAccession(idsFilename)
                        idDict = idExtract.GetIDToNameDict()
                        fullDict.update(idDict)
                    MCL.CreateOGs(ogs, outputFN, fullDict)
                except:
                    err_text = "ERROR: %s contains a duplicate ID. " % (idsFilename)
                    files.FileHandler.LogFailAndExit(err_text)
        return fullDict


    @staticmethod
    def CreateOrthogroupTable(
        ogs,
        idToNameDict,
        speciesNamesDict,
        speciesToUse,
        resultsBaseFilename
    ):

        nSpecies = len(speciesNamesDict)

        ogs_names = [[idToNameDict[seq] for seq in og] for og in ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in ogs]

        # write out
        outputFilename = resultsBaseFilename + ".tsv"
        outputFilename_counts = resultsBaseFilename + ".GeneCount.tsv"
        singleGeneFilename = resultsBaseFilename + "_UnassignedGenes.tsv"
        with open(outputFilename, util.csv_write_mode) as outputFile, \
            open(singleGeneFilename, util.csv_write_mode) as singleGeneFile, \
                open(outputFilename_counts, util.csv_write_mode) as outFile_counts:
            fileWriter = csv.writer(outputFile, delimiter="\t")
            fileWriter_counts = csv.writer(outFile_counts, delimiter="\t")
            singleGeneWriter = csv.writer(singleGeneFile, delimiter="\t")
            for writer in [fileWriter, singleGeneWriter]:
                row = ["Orthogroup"] + [speciesNamesDict[index] for index in speciesToUse]
                writer.writerow(row)
            fileWriter_counts.writerow(row + ['Total'])
            for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
                ogDict = defaultdict(list)
                row = ["OG%07d" % iOg]
                thisOutputWriter = fileWriter
                # separate it into sequences from each species
                if len(og) == 1:
                    row.extend(['' for x in range(nSpecies)])
                    row[speciesToUse.index(og[0][0]) + 1] = og_names[0]
                    thisOutputWriter = singleGeneWriter
                else:
                    for (iSpecies, iSequence), name in zip(og, og_names):
                        ogDict[speciesToUse.index(iSpecies)].append(name)
                    for iSpecies in range(nSpecies):
                        row.append(", ".join(sorted(ogDict[iSpecies])))
                    counts = Counter([iSpecies for iSpecies, _ in og])
                    counts_row = [counts[iSpecies] for iSpecies in speciesToUse]
                    fileWriter_counts.writerow(row[:1] + counts_row + [sum(counts_row)])
                thisOutputWriter.writerow(row)


    @staticmethod
    def SingleGeneWriter(
        ogs,
        idToNameDict,
        speciesNamesDict,
        speciesToUse,
        resultsBaseFilename
    ):

        nSpecies = len(speciesNamesDict)

        ogs_names = [[idToNameDict[seq] for seq in og] for og in ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in ogs]

        # write out
        outputFilename = resultsBaseFilename + ".tsv"
        outputFilename_counts = resultsBaseFilename + ".GeneCount.tsv"
        singleGeneFilename = resultsBaseFilename + "_UnassignedGenes.tsv"
        with open(outputFilename, util.csv_write_mode) as outputFile, \
            open(singleGeneFilename, util.csv_write_mode) as singleGeneFile, \
                open(outputFilename_counts, util.csv_write_mode) as outFile_counts:
            fileWriter = csv.writer(outputFile, delimiter="\t")
            fileWriter_counts = csv.writer(outFile_counts, delimiter="\t")
            singleGeneWriter = csv.writer(singleGeneFile, delimiter="\t")
            for writer in [fileWriter, singleGeneWriter]:
                row = ["Orthogroup"] + [speciesNamesDict[index] for index in speciesToUse]
                writer.writerow(row)
            fileWriter_counts.writerow(row + ['Total'])
            for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
                ogDict = defaultdict(list)
                row = ["OG%07d" % iOg]
                thisOutputWriter = fileWriter
                # separate it into sequences from each species
                if len(og) == 1:
                    row.extend(['' for x in range(nSpecies)])
                    row[speciesToUse.index(og[0][0]) + 1] = og_names[0]
                    thisOutputWriter = singleGeneWriter
                else:
                    for (iSpecies, iSequence), name in zip(og, og_names):
                        ogDict[speciesToUse.index(iSpecies)].append(name)
                    for iSpecies in range(nSpecies):
                        row.append(", ".join(sorted(ogDict[iSpecies])))
                    counts = Counter([iSpecies for iSpecies, _ in og])
                    counts_row = [counts[iSpecies] for iSpecies in speciesToUse]
                    fileWriter_counts.writerow(row[:1] + counts_row + [sum(counts_row)])
                thisOutputWriter.writerow(row)