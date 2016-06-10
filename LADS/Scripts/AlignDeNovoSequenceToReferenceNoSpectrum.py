import sys
import os
import networkx as nx


sys.path.insert(1, os.path.abspath('../'))

import numpy as np
import time
import subprocess
import shlex
import pickle
import itertools
from collections import defaultdict, deque

import ArgLib
import DataFile
import Analytics as An
import SixFrameDNATranslator as SFDT
import Constants

TAG_NODE_SCORE = 1
ISO_CHANGE_PENALTY = -0.5
MOD_PENALTY = -1
UNDEF_MOD_PENALTY = -2

###################
# KNOWN BUGS
# Adding arginine mass mod (+156) instead of matching an isobaric sub to arginine (due to scoring matching AA +1), penalize mods more?
# Allowing insertion and deletion of AAs instead of matching isobaric sub (due to scoring matching AA + 1), penalize mods more?
####################


# modHash is dictionary mapping mod name to its frequency
# Filtering Criteria: terminal specificity, peptide length, alignment score, mod types, longest dissimilarity
def isValidMod(modInfo, modHash = None, modThreshold = 1000, termAA = ['K']):

    mods = eval(modInfo['Modifications'])
    alignScore = eval(modInfo['Alignment Score'])
    refSeqs = eval(modInfo['DB Peptide'])

    def termCompatible(refSeqs):
        for peptide in list(itertools.chain.from_iterable(refSeqs.values())):
            if peptide[0] in termAA or peptide[-3] in termAA:
                return True
    
    def allModsIsobaric(mods):
        if any([mod[0] != 'Isobaric Substitution' for mod in mods]):
            return False
        else:
            return True

    def modsPrevalent(mods):
        if modHash == None:
            return False
        else:
            for mod in mods:
                if mod[0] in ['Undefined Mass Shift', 'Deletion', 'Insertion'] or '->' in mod[0] or modHash[mod[0]] < modThreshold:
                    return False

            return True
    
    if alignScore == None and modInfo['Accuracy'] != 'None':
        return True
    elif allModsIsobaric(mods) and (alignScore > 6 or (alignScore > 3 and termCompatible(refSeqs))):
        return True
    elif alignScore > 6 and (modsPrevalent(mods) or termCompatible(refSeqs)): 
        return True
    else:
        return False
    


def grabGreedyTags(diags):
    acceptedTagRegions = []
    for diag in diags:
        diagX, diagY = zip(*diag)
        
        acceptTag = True
        for tag in acceptedTagRegions:
            tagX, tagY = zip(*tag)
            print (not set(tagX) & set(diagX)), (not set(tagY) & set(diagY))
            if set(tagX) & set(diagX) or set(tagY) & set(diagY):
                acceptTag = False

        if acceptTag:
            acceptedTagRegions += [diag]

    return acceptedTagRegions
    
                                                 

def assembleSequenceTags(indicesDict):
    overlapRegions = []
    tagInds = sorted(indicesDict.keys())


    # Get overlapping sequence tag match regions
    diagHash = defaultdict(deque)
    for tagInd in tagInds:
        for match in indicesDict[tagInd]:
            diagHash[tagInd - match] += [(tagInd, match)]

    #print 'Diag Hash', diagHash
    #Break Apart non-contiguous diagonals
    diags = []
    for diag in diagHash.values():
        tempDiag = [diag.popleft()]
        while diag:
            pair = diag.popleft()
            if tempDiag[-1][1] + 1 == pair[1]:
                tempDiag += [pair]
            else:
                diags += [tempDiag]
                tempDiag = [pair]

        diags += [tempDiag]
            

    return diags
            
                
def generateSequenceTags(deNovoPept, dbPept, tagLength = 3, ambigAA = 'X'):
    indicesDict = {}

    for i in range(len(deNovoPept) - tagLength + 1):
        tag = deNovoPept[i:i+tagLength]
        if ambigAA not in tag:
            matches = list(findAll(tag, dbPept))
            if matches:
                indicesDict[i] = matches

    #print 'Indices Dict', indicesDict
    # Assign tags to peptide positions in a greedy fashion
    acceptedTagRegions = assembleSequenceTags(indicesDict)

    #print 'Accepted Tag Regions', acceptedTagRegions
    acceptedTags = []
    for tagRegion in sorted(acceptedTagRegions, key = lambda t: t[0][0]):
        acceptedTags += [((tagRegion[0][0], tagRegion[-1][0] + tagLength), (tagRegion[0][1], tagRegion[-1][1] + tagLength))]

    return acceptedTags

def greedyTagsToAlignment(deNovoPept, dbPept, acceptedTags):
    alignment = [deNovoPept[:acceptedTags[0][0][0]], dbPept[:acceptedTags[0][0][1]]]

    for i, tag in enumerate(acceptedTags):
        alignment[0] += '*%s*' % (deNovoPept[tag[0][0]:tag[1][0]],)
        alignment[1] += '*%s*' % (dbPept[tag[0][1]:tag[1][1]],)

        if i < len(acceptedTags) - 1:
            alignment[0] += deNovoPept[tag[1][0]:acceptedTags[i+1][0][0]]
            alignment[1] += dbPept[tag[1][1]:acceptedTags[i+1][0][1]]

    alignment[0] += deNovoPept[tag[1][0]:]
    alignment[1] += dbPept[tag[1][1]:]

    return alignment


def generateStartAndEndTags(deNovoPept, dbPept):
    startTags = []
    endTags = []

    for i in range(len(dbPept) + 1):
        startTags += (((0,0), (i,i)),)
        endTags += (((len(deNovoPept), len(deNovoPept)), (i,i)),)

    return startTags, endTags

def getSequenceTagGraph(startTags, endTags, sequenceTags):
    tagGraph = nx.DiGraph()

    #print 'Sequence Tags', sequenceTags
    for tag in startTags:
        tagGraph.add_node(tag, position="start")
    for tag in endTags:
        tagGraph.add_node(tag, position="end")
    for tag in sequenceTags:
        tagGraph.add_node(tag, position="internal")

    for startTag in startTags:
        for internalTag in sequenceTags:
            if startTag[0][1] <= internalTag[0][0] and startTag[1][1] <= internalTag[1][0]:
                tagGraph.add_edge(startTag, internalTag)

    for endTag in endTags:
        for internalTag in sequenceTags:
            if internalTag[0][1] <= endTag[0][0] and internalTag[1][1] <= endTag[1][0]:
                tagGraph.add_edge(internalTag, endTag)

    for i in range(len(sequenceTags)):
        for j in range(i):
            tag1 = sequenceTags[i]
            tag2 = sequenceTags[j]

            if tag1[0][1] <= tag2[0][0] and tag1[1][1] <= tag2[1][0]:
                tagGraph.add_edge(tag1, tag2)
            elif tag2[0][1] <= tag1[0][0] and tag2[1][1] <= tag1[1][0]:
                tagGraph.add_edge(tag2, tag1)

    return tagGraph

def alignDeNovoToDBSequence(deNovoPeptWithMods, deNovoPept, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = None, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-2, undefModPenalty=-3, defaultScore=0):
    deNovoPRMLadder = An.getPRMLadder(deNovoPeptWithMods, ambigEdges = deNovoAmbigEdges, addEnds=True)
    #print deNovoPRMLadder

    print 'De Novo', deNovoPept
    print 'DB', dbPept
    
    dbPRMLadder = An.getPRMLadder(dbPept, addEnds=True)

    startTags, endTags = generateStartAndEndTags(deNovoPept, dbPept)
    sequenceTags = generateSequenceTags(deNovoPept, dbPept, tagLength=tagLength)

    tagGraph = getSequenceTagGraph(startTags, endTags, sequenceTags)

    maxScore = None
    maxScoringTag = None
    
    #print sorted(tagGraph.nodes(data=True))
    #print sorted(tagGraph.edges(data=True))
    for tag in nx.topological_sort(tagGraph):
        nodeScore = tag[0][1] - tag[0][0]
        #print 'Tag', tag
        for prevTag in tagGraph.predecessors(tag):
            nModSymbol = None
            # Define terminus of peptide for modification annotation
            if tagGraph.node[prevTag]['position'] == 'start':
                term = 'N-term'
            elif tagGraph.node[tag]['position'] == 'end':
                term = 'C-term'
            else:
                term = None

            
            refMass = dbPRMLadder[tag[1][0]] - dbPRMLadder[prevTag[1][1]]
            deNovoMass = deNovoPRMLadder[tag[0][0]] - deNovoPRMLadder[prevTag[0][1]]
            refSubSeq = dbPept[prevTag[1][1]:tag[1][0]]
            deNovoSubSeq = deNovoPept[prevTag[0][1]:tag[0][0]]

            mods = resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hashedUnimodDict, unimodDict, paramsDict, term=term, nModSymbol=nModSymbol)
            modPenalty = defModPenalty
            for mod in mods:
                if 'Isobaric Substitution' == mod[0]:
                    modPenalty = isobaricPenalty
                elif 'Insertion' == mod[0] or 'Deletion' == mod[0]:
                    modPenalty = inDelPenalty
                elif 'Undefined Mass Shift' == mod[0]:
                    modPenalty = undefModPenalty
            if not mods:
                modPenalty = 0

            tagGraph.edge[prevTag][tag]['edgeScore'] = nodeScore + modPenalty
            tagGraph.edge[prevTag][tag]['mods'] = mods

            print prevTag, tag, deNovoSubSeq, refSubSeq, mods
            
            if 'score' not in tagGraph.node[prevTag]:
                tagGraph.node[prevTag]['score'] = defaultScore

            try:
                tagGraph.node[tag]['score'] = max(tagGraph.node[tag]['score'], tagGraph.node[prevTag]['score'] + nodeScore + modPenalty)
            except KeyError:
                tagGraph.node[tag]['score'] = tagGraph.node[prevTag]['score'] + nodeScore + modPenalty

            if tagGraph.node[tag]['position'] == 'end' and tagGraph.node[tag]['score'] > maxScore:
                maxScore = tagGraph.node[tag]['score']
                maxScoringTag = tag

    if maxScoringTag != None:
        return getBestAlignment(tagGraph, dbPept, maxScore, maxScoringTag)
    else:
        return None, None, None
    

def calculateMassDeltaFrequencies(precMassArr, epsilon=0.01):
    precMassMatrix = np.reshape(np.tile(precMassArr, precMassArr.size), (precMassArr.size, precMassArr.size))
    deltaMatrix = precMassMatrix - precMassMatrix.transpose()

    freqHash = defaultdict(int)
    for entry in deltaMatrix.flat:
        hMass = hashMass(entry, epsilon)
        for m in range(hMass-1, hMass+2):
            freqHash[m] += 1

    return freqHash

def getScoringModel(freqHash):
    return 0

    
def getBestAlignment(tagGraph, dbPept, maxScore, maxScoringTag):
    modListReversed = tuple()

    endInd = maxScoringTag[1][0]
    currentNode = maxScoringTag
    currentScore = maxScore
    while True:
        for node in tagGraph.predecessors(currentNode):
            if currentScore - tagGraph.edge[node][currentNode]['edgeScore'] == tagGraph.node[node]['score']:
                modListReversed += tagGraph.edge[node][currentNode]['mods']
                if tagGraph.node[node]['position'] == 'start':
                    return dbPept[node[1][0]:endInd], modListReversed[::-1], maxScore

                currentNode = node
                currentScore = tagGraph.node[node]['score']
                break
            
    

def getAlignment(deNovoPept, dbPept, AAMap, scoreMatrix, gapOpenPen=-5, gapExtendPen=0):
    alignment = An.alignSequences(deNovoPept, dbPept, AAMap, scoreMatrix, gapOpenPen, gapExtendPen)[1][0]
    return alignment

def findAll(sub, string):
    index = -1
    try:
        while True:
            index = string.index(sub, index + 1)
            yield index
    except ValueError:
        pass

def updateSeqCount(seqCountArr, seqModArr, proteinSeq, peptide, modList = None, countUp = 1, maxModLength=5):
    #print proteinSeq, peptide
    startPos = proteinSeq.index(peptide)
    #print startPos, len(proteinSeq)
    for i in range(startPos, startPos + len(peptide)):
        seqCountArr[i] += countUp

    
    if modList != None:
        for mod in modList:
            if 'Isobaric' not in mod[0][0] and len(mod[2]) < maxModLength:
                seqModArr[startPos + peptide.index(mod[2])] += 1
    

def getConnectedDisagreementRegions(disagreeArr, trueVal = 0):
    intervals = []

    startInd = -1
    endInd = -1
    for i, elem in enumerate(disagreeArr):
        if elem == trueVal:
            if startInd == -1:
                startInd = i
            endInd = i+1
        else:
            if startInd != -1:
                intervals += [(startInd, endInd)]
            startInd = endInd = -1

    if startInd != -1:
        intervals += [(startInd, endInd)]
        
    return intervals
                

def hashUnimodDict(unimodDict, epStep=0.0025, maxEp=0.1):
    hUnimodDict = {}

    maxIntEp = int(np.ceil(maxEp/epStep))
    minIntEp = int(np.floor(-maxEp/epStep))

    for mod in unimodDict:
        massKey = np.round(unimodDict[mod]['mass']/epStep)

        for intEp in range(minIntEp, maxIntEp+1):
            try:
                hUnimodDict[intEp + massKey][mod] = epStep*intEp
            except KeyError:
                hUnimodDict[intEp + massKey] = {mod: epStep*intEp}

    return hUnimodDict

def hashMass(modMass, epStep=0.0025):
    return int(np.round(modMass/epStep))

def getAlignedIndsMap(alignment):
    alignedIndsMap = {'De Novo': {}, 'Ref': {}}
    numAADeNovo = 0
    numAARef = 0
    
    for i in range(len(alignment[0])):
        alignedIndsMap['De Novo'][i] = numAADeNovo
        alignedIndsMap['Ref'][i] = numAARef
        if alignment[0][i] != '-':
            numAADeNovo += 1
        if alignment[1][i] != '-':
            numAARef += 1

    alignedIndsMap['De Novo'][len(alignment[0])] = numAADeNovo
    alignedIndsMap['Ref'][len(alignment[0])] = numAARef
    return alignedIndsMap

def isValidModLocation(unmodSeq, term, location):
    if location[0] in unmodSeq:
        if location[1] == 'Anywhere':
            return True
        elif location[1] == term:
            if term == 'N-term' and unmodSeq[0] == location[0]:
                return True
            elif term == 'C-term' and unmodSeq[-1] == location[0]:
                return True
            else:
                return False
        else:
            return False
    elif location[0] == term:
        return True
    else:
        return False

def resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSeq, unmodSeq, paramsDict, term = None, nModSymbol=None, epsilon=0.02):

    possUnmodMasses = [refMass]
    for modEntry in paramsDict['Static Mods']:
        if modEntry[1] in unmodSeq:
            possUnmodMasses += [refMass - float(paramsDict['Static Mods'][modEntry])]

    for modEntry in paramsDict['Diff Mods'].values():
        if modEntry[1] in unmodSeq or (term != None and modEntry[1].lower() == term.lower()):
            possUnmodMasses += [refMass + float(modEntry[2])]
            
    matchedMod = (None, None, 10*epsilon)
    for unmodMass in possUnmodMasses:

        # see if iso substitution was masked by diff mod
        if abs(unmodMass - deNovoMass) < epsilon:
            return ('Isobaric Substitution', 0, deNovoMass - unmodMass, deNovoMass, deNovoSeq, unmodSeq)
        
        try:
            possMods = hUnimodDict[hashMass(deNovoMass - unmodMass)]
        except KeyError:
            continue

        for mod in possMods:
            if abs(possMods[mod]) < epsilon:
                locations = unimodDict[mod]['locations']
                for location in locations:
                    if isValidModLocation(unmodSeq, term, location):
                        if abs(possMods[mod]) < abs(matchedMod[2]):
                            matchedMod = (mod, unimodDict[mod]['mass'], possMods[mod], deNovoMass, deNovoSeq, unmodSeq)
                        break

    return matchedMod

# Returns a tuple of modifications found in the given interval
def resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hUnimodDict, unimodDict, paramsDict, term=None, nModSymbol=None, epsilon=0.02):
    if not refSubSeq and not deNovoSubSeq:
        return tuple()
    elif not refSubSeq:
        return (('Insertion', deNovoMass, 0, deNovoMass, deNovoSubSeq, refSubSeq),)
    elif not deNovoSubSeq:
        return (('Deletion', -refMass, 0, deNovoMass, deNovoSubSeq, refSubSeq),)
    elif abs(deNovoMass - refMass) < epsilon:
        return (('Isobaric Substitution', 0, deNovoMass - refMass, deNovoMass, deNovoSubSeq, refSubSeq),)
    else:
        matchedMod = resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSubSeq, refSubSeq, paramsDict, term = term, nModSymbol=nModSymbol, epsilon=epsilon)
        if matchedMod[0] != None:
            return (matchedMod,)
        else:
            return (('Undefined Mass Shift', None, deNovoMass - refMass, deNovoMass, deNovoSubSeq, refSubSeq),)

def getAccAndPrecForModRefPeptide(modList, newRefEndInds, deNovoSeq, deNovoUnmodSeq, refSeq, alignedIndsMap, deNovoAmbigEdges=[]):
    prevIntervalStart = newRefEndInds['start']
    #print 'End Inds', newRefEndInds
    tempSeq = ''
    tempAmbigEdges = []
    for interval in sorted(modList):
        if 'Isobaric' not in modList[interval][0][0] and not ('Insertion' in modList[interval][0][0] and (alignedIndsMap['De Novo'][interval[0]] < 2 or (len(deNovoUnmodSeq) - alignedIndsMap['De Novo'][interval[1]]) < 2)):
            tempSeq += refSeq[prevIntervalStart:alignedIndsMap['Ref'][interval[0]]] + 'X'
            #print 'Mod list interval', modList[interval]
            tempAmbigEdges += [(0, modList[interval][0][3])]
            #print 'temp ambig edges', tempAmbigEdges
            prevIntervalStart = alignedIndsMap['Ref'][interval[1]]
    #print 'TempSeq', tempSeq, tempAmbigEdges
    tempSeq += refSeq[prevIntervalStart:(len(refSeq) + newRefEndInds['end'])]
    #print deNovoSeq, refSeq, tempSeq, deNovoAmbigEdges, tempAmbigEdges
    comp = An.comparePeptideResults(deNovoSeq, tempSeq, ambigEdges1=deNovoAmbigEdges, ambigEdges2=tempAmbigEdges, ppm=10)
    return comp[0], comp[1]

def compareSequences(deNovoPep, deNovoUnmodPep, refPep, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = [], epsilon = 0.02):
    if 'X' in refPep:
        refPep = refPep.translate(None, 'X')

    # KLUDGE: REMOVE WHEN REWRITE
    #deNovoPep = An.stripModifications(deNovoPep, noRemove=['#', '*'])
    
    alignment = getAlignment(deNovoUnmodPep, refPep, AAMap, scoreMatrix)
    alignedIndsMap = getAlignedIndsMap(alignment)
    

    disagreeArr = [1 if alignment[0][i] == alignment[1][i] else 0 for i in range(len(alignment[0]))]
    intervals = getConnectedDisagreementRegions(disagreeArr)

    try:
        refPRMLadder = An.getPRMLadder(refPep)
    except KeyError:
        return None
    
    deNovoPRMLadder = An.getPRMLadder(deNovoPep, ambigEdges=deNovoAmbigEdges)

    allResolved = True
    modList = {}
    newRefEndInds = {'start': 0, 'end': 0}

    # rough check of whether or not intervals can be easily explained
    for interval in intervals:
        deNovoSubSeq = deNovoUnmodPep[alignedIndsMap['De Novo'][interval[0]]:alignedIndsMap['De Novo'][interval[1]]]
        refSubSeq = refPep[alignedIndsMap['Ref'][interval[0]]:alignedIndsMap['Ref'][interval[1]]]

        if alignedIndsMap['De Novo'][interval[0]] == 0:
                term = 'N-term'
        elif alignedIndsMap['De Novo'][interval[1]] == len(deNovoUnmodPep):
            term = 'C-term'
        else:
            term = None

        if deNovoSubSeq != '' and refSubSeq != '':
            deNovoMass = deNovoPRMLadder[alignedIndsMap['De Novo'][interval[1]]] - deNovoPRMLadder[alignedIndsMap['De Novo'][interval[0]]]
            if term == None:
                refMass = refPRMLadder[alignedIndsMap['Ref'][interval[1]]] - refPRMLadder[alignedIndsMap['Ref'][interval[0]]]
                modList[interval] = resolveInterval(refMass, deNovoMass, refSubSeq, hashedUnimodDict, unimodDict, paramsDict, term=term, epsilon=epsilon), deNovoSubSeq, refSubSeq
            else:
                minSizedMod = ((None, None, 10000000,),)
                for i in range(len(refSubSeq)):
                    if term == 'N-term':
                        refMass = refPRMLadder[alignedIndsMap['Ref'][interval[1]]] - refPRMLadder[alignedIndsMap['Ref'][interval[0]] + i]
                        subRefSubSeq = refSubSeq[i:]
                    else:
                        refMass = refPRMLadder[alignedIndsMap['Ref'][interval[1]] - i] - refPRMLadder[alignedIndsMap['Ref'][interval[0]]]
                        subRefSubSeq = refSubSeq[:-i]
                    mod = resolveInterval(refMass, deNovoMass, subRefSubSeq, hashedUnimodDict, unimodDict, paramsDict, term=term, epsilon=epsilon)
                    if 'TX' in deNovoUnmodPep:
                        print deNovoSubSeq, refSubSeq, subRefSubSeq, mod
                    if (abs(minSizedMod[0][2]) > abs(mod[2]) and (minSizedMod[0][0] == None or 'Isobaric' not in minSizedMod[0][0])) or 'Isobaric' in mod[0]:
                        if mod[1] != None or (mod[1] == None and minSizedMod[0][1] == None) or ('Isobaric' in mod[0] and 'Isobaric' not in minSizedMod[0][0]):
                            minSizedMod = mod, deNovoSubSeq, subRefSubSeq
                            if term == 'N-term':
                                newRefEndInds['start'] = i
                            else:
                                newRefEndInds['end'] = -i
                modList[interval] = minSizedMod
                    
        else:
            # Make sure that lack of sequence is due to overhang of reference peptide
            if alignedIndsMap['De Novo'][interval[1]] == 0:
                newRefEndInds['start'] = len(refSubSeq)
            elif alignedIndsMap['De Novo'][interval[0]] == len(deNovoUnmodPep):
                newRefEndInds['end'] = -len(refSubSeq)
#            elif term != None:
#                raise ValueError('Not enough reference sequence provided for resoluton of terminal discrepancies. De Novo: %s, Reference %s' % (deNovoPep, refPep))
            elif term == None:
                if deNovoSubSeq == '':
                    refMass = refPRMLadder[alignedIndsMap['Ref'][interval[1]]] - refPRMLadder[alignedIndsMap['Ref'][interval[0]]]
                    modList[interval] = ('Deletion', refMass, 0, -refMass), deNovoSubSeq, refSubSeq
                else:
                    deNovoMass = deNovoPRMLadder[alignedIndsMap['De Novo'][interval[1]]] - deNovoPRMLadder[alignedIndsMap['De Novo'][interval[0]]]
                    modList[interval] = ('Insertion', deNovoMass, 0, deNovoMass), deNovoSubSeq, refSubSeq

    #print 'Mod List: ', modList
    acc, prec =  getAccAndPrecForModRefPeptide(modList, newRefEndInds, deNovoPep, deNovoUnmodPep, refPep, alignedIndsMap, deNovoAmbigEdges)
    
    return modList, newRefEndInds, alignment, acc, prec
            
def loadFASTAs(fastaFile, fromBLAST=True):
    seqDict = {}
    seqGen = SFDT.sequenceGenerator(fastaFile)

    for seqName, seq in seqGen:
        seqNameTruncated = seqName.split(' ')[0][1:]
        if fromBLAST and len(seqNameTruncated) >= 63:
            seqNameTruncated = seqNameTruncated[:63] + '...'
        seqDict[seqNameTruncated] = seq

    return seqDict


def getExtendedSequence(seqDict, reference, nTermMass=0, cTermMass=0):
    protSeq = seqDict[reference[0]]
    baseSeq = protSeq[reference[1]-1:reference[2]]
    print baseSeq, nTermMass, cTermMass, reference
    
    if nTermMass > 0:
        nTermSubSeq, nTermSubMass = '', 0
        for aa in removeNoncanonicalAminoAcids(protSeq[:reference[1]-1])[::-1]:
            nTermSubMass += Constants.aminoacids[aa][2]
            if nTermSubMass >= nTermMass:
                break
            else:
                nTermSubSeq = aa + nTermSubSeq
                print nTermSubSeq, nTermSubMass
                
        baseSeq = nTermSubSeq + baseSeq

    if cTermMass > 0:
        cTermSubSeq, cTermSubMass = '', 0
        for aa in removeNoncanonicalAminoAcids(protSeq[reference[2]:]):
            cTermSubMass += Constants.aminoacids[aa][2]
            if cTermSubMass >= cTermMass:
                break
            else:
                cTermSubSeq = cTermSubSeq + aa
                    
        baseSeq = baseSeq + cTermSubSeq

    return baseSeq
                                                                                                                                                                                        


def getReferenceSequence(proteinSeq, peptide, start=None, end=None):
    if start == None:
        start = list(findAll(peptide,proteinSeq))
        end = [startInd + len(peptide) for startInd in start]
    proteinSeq = '-' + proteinSeq + '-'
    return [proteinSeq[start[i]] + '.' + proteinSeq[start[i]+1:end[i]+1] + '.' + proteinSeq[end[i]+1] for i in range(len(start))]

def getStartAndEndInds(proteinSeq, peptide):
    startInd = proteinSeq.index(peptide)
    return startInd, startInd+len(peptide)

def removeNoncanonicalAminoAcids(subjectSeq):
    return subjectSeq.replace('Z', 'Q').replace('B', 'N').replace('X', '').replace('*', '').replace('U', 'C').replace('L', 'I')

if __name__ == '__main__':
    print 'This program will take a tdv file of  BLAST results indexed by scanF and attempt to explain the discrepancies using the unimod modification dictionary. Mainprogname is a dict mapping the name of the score and peptide field to their corresponding URIs. Fields are ScanF, Peptide, References, Ambig Edges (optional), Ref Peptide, and Num Identical'
    options = ArgLib.parse(['init', 'output', 'ppmstd', 'comp', 'unimoddict', 'mainprogname'], [{'opts': ('-f', '--fasta'), 'attrs': {'type': 'string', 'dest': 'fasta', 'help': 'Location of reference fasta containing reference proteins (same file used to generate BLAST DB).'}}])

    paramsDict = ArgLib.parseInitFile(options.init, options)
    seqDict = loadFASTAs(options.fasta)

    infoDict = eval(options.mainprogname)

    with open(options.unimoddict) as fin:
        unimodDict = pickle.load(fin)
    hashedUnimodDict = hashUnimodDict(unimodDict)

    outFile = open(options.output, 'w')
    cols = ['ScanF', 'Score', 'Peptide', 'Unmod Peptide', 'References', 'Modifications', 'DB Peptide', 'Alignment Score']
    if 'Ambig Edges' in infoDict:
        cols.insert(2, 'Ambig Edges')
        
    outFile.write('\t'.join([col for col in cols]) + '\n')

    for entry in DataFile.getScanInfo(options.comp, delimiter='\t'):
        scanData = {}
        scanData['ScanF'] = entry[infoDict['ScanF']]
        scanData['Peptide'] = entry[infoDict['Peptide']]
        scanData['Unmod Peptide'] = An.stripModifications(scanData['Peptide'], noRemove=[])
        scanData['Score'] = entry[infoDict['Score']]
        scanData['Alignment Score'] = None
        
        if 'Ambig Edges' in infoDict:
            ambigEdges = eval(entry[infoDict['Ambig Edges']])
            scanData['Ambig Edges'] = ambigEdges
        else:
            ambigEdges = []
        deNovoPRMLadder = An.getPRMLadder(scanData['Peptide'], ambigEdges=ambigEdges)
                
        refList = eval(entry[infoDict['References']])
        subjSequence = seqDict[refList[0][0]][refList[0][1]-1:refList[0][2]]

        if scanData['Unmod Peptide'] == subjSequence:
            scanData['Modifications'] = []
            
            refSeqDict = {}
            for reference in refList:
                refSeqDict[reference] = getReferenceSequence(seqDict[reference[0]], subjSequence, start= [reference[1]-1], end= [reference[2]])
                
            scanData['DB Peptide'] = refSeqDict
            scanData['References'] = [ref[0] for ref in refList]

        else:
            modLists = {}
            refDict = {}            

            for reference in refList:

                nTermMass, cTermMass = 0, 0
                if int(entry['Query Start']) > 1:
                    nTermMass = deNovoPRMLadder[int(entry['Query Start'])-1] + 250
                if int(entry['Query End']) < len(scanData['Unmod Peptide']):
                    cTermMass = deNovoPRMLadder[-1] - deNovoPRMLadder[int(entry['Query End'])] + 250

                print scanData['ScanF'], subjSequence, scanData['Unmod Peptide'], nTermMass, cTermMass
                extendedSequence = getExtendedSequence(seqDict, reference, nTermMass, cTermMass)
                if extendedSequence not in modLists:
                    dbPept, modList, maxScore = alignDeNovoToDBSequence(scanData['Peptide'], scanData['Unmod Peptide'], extendedSequence, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = ambigEdges, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-3, undefModPenalty=-3, defaultScore=0)

                    modLists[extendedSequence] = [dbPept, modList, maxScore]
                    refDict[extendedSequence] = [reference]
                else:
                    refDict[extendedSequence] += [reference]


            if len(modLists) == 0:
                continue

            extendedSequence = max(modLists.iterkeys(), key=lambda k: modLists[k][2])

            if modLists[extendedSequence][0] == None:
                continue
            
            scanData['References'] = [ref[0] for ref in refDict[extendedSequence]]

            refSeqDict = {}
            for reference in refDict[extendedSequence]:
                refSeqDict[reference] = getReferenceSequence(seqDict[reference[0]], modLists[extendedSequence][0])
                
                                   
            scanData['DB Peptide'] = refSeqDict
            scanData['Modifications'] = modLists[extendedSequence][1]
            scanData['Alignment Score'] = modLists[extendedSequence][2]

        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

    outFile.close()
            
