import DataFile
import ArgLib
import Constants

# set both leucine and isoleucine codons to translate to isoleucine
codonDict = {

'TTT': 'F',
'TTC': 'F',
'TTA': 'L',
'TTG': 'L',
'CTT': 'L',
'CTC': 'L',
'CTA': 'L',
'CTG': 'L',
'ATT': 'I',
'ATC': 'I',
'ATA': 'I',
'ATG': 'M',
'GTT': 'V',
'GTC': 'V',
'GTA': 'V',
'GTG': 'V',
'TCT': 'S',
'TCC': 'S',
'TCA': 'S',
'TCG': 'S',
'CCT': 'P',
'CCC': 'P',
'CCA': 'P',
'CCG': 'P',
'ACT': 'T',
'ACC': 'T',
'ACA': 'T',
'ACG': 'T',
'GCT': 'A',
'GCC': 'A',
'GCA': 'A',
'GCG': 'A',
'TAT': 'Y',
'TAC': 'Y',
'TAA': '*',
'TAG': '*',
'CAT': 'H',
'CAC': 'H',
'CAA': 'Q',
'CAG': 'Q',
'AAT': 'N',
'AAC': 'N',
'AAA': 'K',
'AAG': 'K',
'GAT': 'D',
'GAC': 'D',
'GAA': 'E',
'GAG': 'E',
'TGT': 'C',
'TGC': 'C',
'TGA': '*',
'TGG': 'W',
'CGT': 'R',
'CGC': 'R',
'CGA': 'R',
'CGG': 'R',
'AGT': 'S',
'AGC': 'S',
'AGA': 'R',
'AGG': 'R',
'GGT': 'G',
'GGC': 'G',
'GGA': 'G',
'GGG': 'G',

}



# mapping onto antisense strand, includes N:N to map undetermined nucleotides onto undetermined nucleotides
antisenseMap = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', 'N': 'N'}

def calculateRandomMatchPeptideProbabilities(fastaFileName, pepLengths=range(5,20), ILEqual=False):

    probDict = {}
    numAA = 20
    if ILEqual:
        numAA -= 1

    for pepLength in pepLengths:
        seqGen = sequenceGenerator(fastaFileName)
        uniquePeptSet = set()
        for seqName, sequence in seqGen:
            startInd = 0
            while startInd <= len(sequence) - pepLength:
                subSeq = sequence[startInd:startInd+pepLength]
                if all(aa in Constants.aminoacids for aa in subSeq):
                    uniquePeptSet.add(subSeq)
                startInd += 1

        probDict[pepLength] = '%.5g' % ((float(len(uniquePeptSet))/(numAA**pepLength)),)

    return probDict

# frames are 1,2,3,-1,-2,-3
def getTranslation(sequence, frame):
    if frame < 0:
        # get antisense strand and reverse orientation
        sequence = ''.join([antisenseMap[nucleotide] for nucleotide in sequence])[::-1]
        frame = frame * -1

    translatedSequence = []
    for i in range(frame-1, len(sequence), 3):
        try:
            translatedSequence += [codonDict[sequence[i:i+3]]]
        except KeyError:
            translatedSequence += ['X']

    return ''.join(translatedSequence)
    
    
def sequenceGenerator(fastaFileName):
    fastaFile = open(fastaFileName)

    seqName = None
    for line in fastaFile:
        if line[0] == '>':
            if seqName != None:
                yield seqName, ''.join(sequence)

            seqName = line.strip()
            sequence = []
        else:
            sequence += [line.strip()]

    yield seqName, ''.join(sequence)




# seqStart and seqEnd are defined as distance away from start of chromosome with first position indexed as one, seqStart is inclusive and seqEnd is exclusive
def chunkAndWriteSequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, lineLength=60, chunkSize=None, dontReadThrough=['*', 'X'], minPeptLength=10, overhang=100):
    i = 0
    while i < len(fullSequence):
        if fullSequence[i] in dontReadThrough:
            i += 1
        else:
            lowestInd = min(fullSequence.find(forbidAA,i) for forbidAA in dontReadThrough)
            if lowestInd == -1:
                lowestInd = len(fullSequence)
                for forbidAA in dontReadThrough:
                    if fullSequence.find(forbidAA,i) != -1:
                        lowestInd = min(lowestInd, fullSequence.find(forbidAA,i))
            
            if lowestInd - i > minPeptLength and (not chunkSize or lowestInd - i <= chunkSize):
                prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, i, lowestInd, lineLength=lineLength)
            elif chunkSize and lowestInd - i > chunkSize:
                for startInd in range(i, lowestInd, chunkSize-overhang):
                    prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, startInd, min(startInd+chunkSize, lowestInd), lineLength=60)
            
            i = lowestInd
                    
                
def prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, startInd, endInd, lineLength=60):
    seqStart = 1 + abs(frame) + 3*startInd
    seqEnd = seqStart + 3*(endInd - startInd)
    if frame < 0:
        seqEnd = ntSeqLength - seqStart + 2
        seqStart = seqEnd - 3*(endInd - startInd)

    seqName = getTransSeqNameForGFY(baseSeqName, frame, '_start_Pos%i_end_Pos%i' % (seqStart, seqEnd))
    writeSequence(outFile, seqName, fullSequence[startInd:endInd], lineLength=lineLength)

def writeSequence(outFile, seqName, sequence, lineLength=60):
    outFile.write(seqName + '\n')
    for i in range(0, len(sequence), lineLength):
        outFile.write(sequence[i:i+lineLength] + '\n')

def getTransSeqNameForGFY(seqName, frame, addStringToBase=''):
    seqNameList = seqName.split(' ')
    seqNameList[0] = seqNameList[0] + ('_+' if frame > 0 else '_') + str(frame) + addStringToBase
    return ' '.join(seqNameList) + ' frame=%i' %(frame,)
    

if __name__ == '__main__':
    print 'This program will take a FASTA file from the --lads argument and output a six-frame translation of the file to output. Number refers to maximum size of sequence in resulting FASTA file. If a chromosomal region exceeds this length with no stop codons, the sequence will be chunked with a 100 aa overhang at each edge. Minimum Length of peptide in FASTA file is 5.'
    options = ArgLib.parse(['lads', 'output', 'number'])

    outFile = open(options.output, 'w')
    #chunkSize = int(options.number)
    
    for seqName, sequence in sequenceGenerator(options.lads):
        for frame in [1, 2, 3, -1, -2, -3]:
        #for frame in [1,2,3]:
            #print seqName, frame
            transSeq = getTranslation(sequence.upper(), frame)
            chunkAndWriteSequence(outFile, seqName, transSeq, frame, len(sequence), lineLength=60, chunkSize=2000, dontReadThrough=['X'], minPeptLength=10, overhang=50)
            #transSeqName = seqName + ('_+' if frame > 0 else '_') + str(frame)
            #transSeqName = getTransSeqNameForGFY(seqName, frame)
            #writeSequence(outFile, transSeqName, transSeq)
            
    outFile.close()
            
