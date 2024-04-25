from .utils import util, files 

def BitScore(sequence):
    pass

def GetMaxBitscores(seqsInfo):
    bit_scores = []
    for iSpecies, iFasta in enumerate(seqsInfo.speciesToUse):
        bit_scores.append(np.zeros(seqsInfo.nSeqsPerSpecies[iFasta]))
        fastaFilename = files.FileHandler.GetSpeciesFastaFN(iFasta)
        current_sequence = ""
        qFirstLine = True
        with open(fastaFilename) as infile:
            for row in infile:
                if len(row) > 1 and row[0] == ">":
                    if qFirstLine:
                        qFirstLine = False
                    else:
                        bit_scores[iSpecies][iCurrentSequence] = BitScore(current_sequence)
                        current_sequence = ""
                    _, iCurrentSequence = util.GetIDPairFromString(row[1:])
                else:
                    current_sequence += row.rstrip()
        bit_scores[iSpecies][iCurrentSequence] = BitScore(current_sequence)
    return bit_scores