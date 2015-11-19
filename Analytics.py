'''
Created on Jul 5, 2011

@author: Arun
'''

from collections import deque
import pickle
import copy
import numpy as np
import time
import svmutil

import Constants
import SpectraAlignment as SA
import DataFile
import ProbNetwork as PN

searchprogs = ['LADS', 'SEQUEST', 'MASCOT', 'SORCERER', 'PepNovo', 'NovoHMM', 'MaxQuant', 'LuteFisk', 'X!Tandem', 'SPIDER', 'PEAKS', 'OMSSA', 'InSpect', 'pNovo', 'Combined']
avgAAMass = 110

def getAlignmentRatios(scanInfoFName, dtaDir, delta, epsilon=0.02):
    scanInfo = DataFile.getScanInfo(scanInfoFName)
    dtaNames = DataFile.getDTAFNamesInDir(dtaDir)
    
    scansToUse = scanInfo
    """
    for i in range(len(scanInfo) - 1):
        if (int(scanInfo[i][0]) + 1 == int(scanInfo[i+1][0])):
            if (scanInfo[i][1] == scanInfo[i+1][1]):
                scansToUse += [scanInfo[i]]
        else:
            scansToUse += [scanInfo[i]]
    """
    ratios = []
    goodRatios = []
    for i in range(len(scansToUse)):
        for j in range(i + 1, len(scansToUse)):
            if j == i + 1:
                print '%s percent done' % str(float(i) / len(scansToUse))
            if np.abs(np.abs(float(scansToUse[i][1]) - float(scansToUse[j][1])) - delta) < epsilon:
                dta1 = '244.%(scanF)04i.%(scanF)04i.1.dta' % {'scanF': int(scansToUse[i][0])}
                dta2 = '244.%(scanF)04i.%(scanF)04i.1.dta' % {'scanF': int(scansToUse[j][0])}
                spec1 = DataFile.getMassIntPairs(dtaDir + dta1)
                spec2 = DataFile.getMassIntPairs(dtaDir + dta2)
                ratio = SA.getSharedPeaksRatio(float(scansToUse[i][1]), spec1, float(scansToUse[j][1]), spec2, epsilon)
                print ratio, scansToUse[i], scansToUse[j]
                ratios.extend([(ratio, scansToUse[i], scansToUse[j])])

    with open('heavylightpairs.txt', 'w') as fout:
        pickle.dump(ratios, fout)
    return ratios

def findDeltaPairs(dtaList, delta, ppm=5, intEp=20):
    precMassArr = np.zeros((len(dtaList), 2))
    pairs = []
    
    for i in range(len(dtaList)):
        precMassArr[i] = [DataFile.getPrecMassAndCharge(dtaList[i])[0], DataFile.getScanNum(dtaList[i])]

    maxPrecMass = np.max(precMassArr, 0)[0]
    epsilon = ppm * 10**-6 * maxPrecMass
    resolution = epsilon/intEp
    
    hPrecMassArr = copy.copy(precMassArr)
    hPrecMassArr[:,0] = np.round(hPrecMassArr[:,0]/resolution)
    hashedDict = {}
    for elem in hPrecMassArr:
        hInd = int(elem[0])
        for hMass in range(hInd-intEp, hInd+intEp+1):
            try:
                hashedDict[hMass] += [(hMass-hInd, elem[1])]
            except KeyError:
                hashedDict[hMass] = [(hMass-hInd, elem[1])]
    
    shiftHashDict = copy.copy(precMassArr)
    shiftHashDict[:,0] = np.round((shiftHashDict[:,0] - delta)/resolution)
    for i, elem in enumerate(shiftHashDict):
        hInd = int(elem[0])
        if hInd in hashedDict:
            for possiblePair in hashedDict[hInd]:
                if abs(possiblePair[0]) * resolution * 10**6/precMassArr[i][0] < ppm:
                    pairs += [(int(possiblePair[1]), int(elem[1]))]

    return pairs

def findDeltaPairsClusters(clusters, scanFDict, delta, ppm=5, intEp=20):
    pairs = []
    precMassArr = np.zeros(len(clusters))
    for i, cluster in enumerate(clusters):
        precMassArr[i] = np.average(np.array([scanFDict[scanF]['precMass'] for scanF in cluster]))
    
    maxPrecMass = np.max(precMassArr)
    epsilon = ppm * 10**-6 * maxPrecMass
    resolution = epsilon/intEp
    
    hPrecMassArr = copy.copy(precMassArr)
    hPrecMassArr = np.round(hPrecMassArr/resolution)
    hashedDict = {}
    for i, elem in enumerate(hPrecMassArr):
        hInd = int(elem)
        for hMass in range(hInd-intEp, hInd+intEp+1):
            try:
                hashedDict[hMass] += [(hMass-hInd, i)]
            except KeyError:
                hashedDict[hMass] = [(hMass-hInd, i)]
    
    shiftHashDict = copy.copy(precMassArr)
    shiftHashDict = np.round((shiftHashDict - delta)/resolution)
    for i, elem in enumerate(shiftHashDict):
        hInd = int(elem)
        if hInd in hashedDict:
            for possiblePair in hashedDict[hInd]:
                if abs(possiblePair[0]) * resolution * 10**6/precMassArr[i] < ppm:
                    pairs += [(possiblePair[1], i)]
    
    return pairs
                 
    
def findSamePrecMassClusters(dtaList, ppm=5):
    precMassArr = np.zeros((len(dtaList), 2))
    for i in range(len(dtaList)):
        precMassArr[i] = [DataFile.getPrecMassAndCharge(dtaList[i])[0], DataFile.getScanNum(dtaList[i])]
        
    precMassArr = precMassArr[np.argsort(precMassArr[:,0])]
    
    clusters = [[i] for i in range(precMassArr.shape[0])]
    
    i = 0
    while i < len(clusters):
        mergeClusters = False
        epsilon = ppm * 10**-6 * precMassArr[clusters[i][0]][0]
        for precMassInd1 in clusters[i]:
            for precMassInd2 in clusters[i - 1]:
                if (np.abs(precMassArr[precMassInd1][0] - precMassArr[precMassInd2][0]) < epsilon):
                    mergeClusters = True
                    break
            
        if mergeClusters:
            clusters[i - 1].extend(clusters[i])
            del clusters[i]
        else:
            i = i + 1
    
    scanFClusters = []
    for cluster in clusters:
        scanFClusters += [[precMassArr[i][1] for i in cluster]]
        
    return scanFClusters

def getSamePeptideClusters(precMassClusters, scanFDict, svmModel, svmRange, ppmSTD=5, cutOff=0):
    trueClusters = []
    
    for cluster in precMassClusters:
        if len(cluster) == 1:
            trueClusters += [cluster]
        else:
            pairIndex = []
            xVals = []
            specs = []
            for i in range(len(cluster)):
                specs +=  [DataFile.getMassIntPairs(scanFDict[cluster[i]]['dta'])]
                
            dMatrix = np.ones((len(cluster), len(cluster))) * -2
            for i in range(len(cluster)):
                for j in range(i+1, len(cluster)):
                    epSTD = ppmSTD * 10 ** -6 * scanFDict[cluster[i]]['precMass']
            
                    SVMClassificationInfo = SA.getSpectraPairInfoForSVMClassification(specs[i], specs[j], scanFDict[cluster[i]]['precMass'], NMod=0, CMod=0, epsilon=2*epSTD)
                    xVals += [SVMClassificationInfo]
                    pairIndex += [(i, j)]
            
            xValsNorm = svmutil.normalize_instances(xVals, svmRange)
            pLabs = svmutil.svm_predict([0]*len(xValsNorm), xValsNorm, svmModel)[0]
#            print pLabs
            for i, pLab in enumerate(pLabs):
            # Scale distances by 4: totalTICRatio, 1: TotalSharedPeaksRatio
                dMatrix[pairIndex[i][0]][pairIndex[i][1]] =  dMatrix[pairIndex[i][1]][pairIndex[i][0]] = xVals[i][1] if pLab==1 else -1

            trueClusters += heirarchicalClusteringAverageLinkage([[scanF] for scanF in cluster], dMatrix, cutOff=cutOff)
    
    return trueClusters
            
def heirarchicalClusteringAverageLinkage(clusters, dMatrix, cutOff=0):
    maxInd = np.where(dMatrix == dMatrix.max())
#    print 'dMatrix', dMatrix
#    print maxInd, dMatrix[maxInd][0], cutOff, dMatrix[maxInd][0] < cutOff
    if dMatrix[maxInd][0] < cutOff:
        return clusters
    else:
        c1, c2 = maxInd[0][0], maxInd[1][0]
#        print 'c1 c2', c1, c2
        dMatrix[c1, :] = (dMatrix[c1, :]*len(clusters[c1]) + dMatrix[c2, :]*len(clusters[c2]))/(len(clusters[c1]) + len(clusters[c2]))
        dMatrix[:, c1] = dMatrix[c1, :]
        
        dMatrix = np.delete(np.delete(dMatrix, c2, 0), c2, 1)
        clusters[c1] = clusters[c1] + clusters[c2] 
        dMatrix[c1, c1] = -2
        del clusters[c2]
        
        return heirarchicalClusteringAverageLinkage(clusters, dMatrix, cutOff=cutOff)
        
        
def getPairDiffDistribution(pairs):
    scanDiffVec = np.array([pair[0]-pair[1] for pair in pairs])
    hist, bin_edges = np.histogram(scanDiffVec)
    print 'Histogram: ', hist
    print 'Bins: ', bin_edges

def getImprVec(shiftLightSpec, heavySpec):
    heavySpecLen = heavySpec.shape[0]
    closestInds = np.searchsorted(heavySpec[:,0], shiftLightSpec[:,0])
    imprVec = np.zeros(shiftLightSpec.shape)
    for i, ind in enumerate(closestInds):
        if ind == 0:
            imprVec[i][1] = np.abs(heavySpec[0][1] - shiftLightSpec[i][1])
            imprVec[i][0] = heavySpec[0][1]
        elif ind == heavySpecLen:
            imprVec[i][1] = np.abs(heavySpec[-1][1] - shiftLightSpec[i][1])
            imprVec[i][0] = heavySpec[-1][1]
        else:
            botOff = np.abs(heavySpec[ind-1][1] - shiftLightSpec[i][1])
            topOff = np.abs(heavySpec[ind][1] - shiftLightSpec[i][1])
            if botOff > topOff:
                imprVec[i] = [heavySpec[ind-1][1], botOff]
            else:
                imprVec[i] = [heavySpec[ind][1], topOff]

    return imprVec

def intensityDotProduct(lightSpec, heavySpec, Nmod, Cmod, epsilon):
    lightSpec[:,1] = lightSpec[:,1]/np.sum(lightSpec[:,1])
    heavySpec[:,1] = heavySpec[:,1]/np.sum(heavySpec[:,1])
    NShiftSpec = copy.copy(lightSpec)
    NShiftSpec[:,0] += Nmod
    CShiftSpec = copy.copy(lightSpec)
    CShiftSpec[:,0] += Cmod
    
    NImprVec = getImprVec(NShiftSpec, heavySpec)
    CImprVec = getImprVec(CShiftSpec, heavySpec)
    diffImp = NImprVec[:,1] - CImprVec[:,1]
    medDiffImp = np.median(diffImp)
    
    heavyPrimeVec = np.zeros(lightSpec.shape[0])
    for i in range(lightSpec.shape[0]):
        if diffImp[i] < medDiffImp:
            heavyPrimeVec[i] = NImprVec[i][0]
        else:
            heavyPrimeVec[i] = CImprVec[i][0]
    
    lightVec = lightSpec[:,1]
    return np.dot(lightVec, heavyPrimeVec)/(np.linalg.norm(lightVec) * np.linalg.norm(heavyPrimeVec))

#Assumption-->all amino acids denoted by single letters
def stripModifications(seq, ambigAA='X', noRemove=[]):
    stripSeq = []
    for aa in seq:
        if aa in Constants.aminoacids or aa == ambigAA or aa in noRemove:
            stripSeq += [aa]
    
    return ''.join(stripSeq)

def writeFASTAFile(compInfo, LADSprogName, infoMap, FASTAout, ambigEdgeCutoff=1):
    outFile = open(FASTAout, 'w')
    for scan in compInfo:
        if scan[LADSprogName + ' ' + infoMap['LADS']['Peptide']] != 'None' and int(scan[LADSprogName + ' Num Ambig Edges']) <= ambigEdgeCutoff:
            outFile.write('>' + scan['ScanF'] + '\n')
            outFile.write(stripModifications(scan[LADSprogName + ' ' + infoMap['LADS']['Peptide']]) + '\n')
    outFile.close() 

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = DataFile.getScanNum(dta)
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'sequenced': False}
    
    return scanFDict

#returns a number from 0 to 1 indicating how confident LADS is in the edge
def getAAConfidence(G, ambigAAMap={'X': '-'}, prevNode=None, nextNode=None):
    fromNodeConf = False
    toNodeConf = False
    seqMap = {}
    for aa in Constants.aminoacids.keys():
        seqMap[aa] = aa
    seqMap.update(ambigAAMap)
        
    if prevNode != None and G.node[prevNode['prm']]['memo'][(seqMap[prevNode['formAA']], seqMap[prevNode['lattAA']])] > 0:
        fromNodeConf = True
    elif prevNode == None:
        fromNodeConf = True
    
    if nextNode != None and G.node[nextNode['prm']]['memo'][(seqMap[nextNode['formAA']], seqMap[nextNode['lattAA']])] > 0:
        toNodeConf = True
    elif nextNode == None:
        toNodeConf = True
    
    return 0.5*fromNodeConf + 0.5*toNodeConf
        

def getSharedPeaksRatio(lightPath, heavyPath, pairConfig, epsilon):
    lightPairs = DataFile.getMassIntPairs(lightPath)
    heavyPairs = DataFile.getMassIntPairs(heavyPath)
    N, C = SA.getNandCIons(lightPairs, heavyPairs, pairConfig['NMod'], pairConfig['CMod'], epsilon=epsilon)
    return SA.getSharedPeaksRatio(lightPairs, heavyPairs, N, C)

def getPairedAndUnpairedSpectra(dtaDir, dtaList, Nmod, Cmod, ppm=5, cutOff=0.1, verbose=False):
    specPairs = []
    unpairedSpecs = []
    delta = Nmod + Cmod
    for i in range(len(dtaList)):
        paired = False
        precMass1 = DataFile.getPrecMassAndCharge(dtaList[i])[0]
        spec1 = DataFile.getMassIntPairs(dtaList[i])
        for j in range(i + 1, len(dtaList)):
            precMass2 = DataFile.getPrecMassAndCharge(dtaList[j])[0]
            epsilon = ppm * 10 ** -6 * max(precMass1, precMass2)
            if np.abs(np.abs(precMass1 - precMass2) - delta) < epsilon:
                spec2 = DataFile.getMassIntPairs(dtaList[j])
                if precMass1 < precMass2:
                    N, C = SA.getNandCIons(spec1, spec2, Nmod, Cmod, epsilon=epsilon)
                    ratio = SA.getSharedPeaksRatio(spec1, spec2, N, C)
                else:
                    N, C = SA.getNandCIons(spec2, spec1, Nmod, Cmod, epsilon=epsilon)
                    ratio = SA.getSharedPeaksRatio(spec2, spec1, N, C)
                if ratio > cutOff:
                    if verbose:
                        print 'Pair found', dtaList[i], dtaList[j]
                    paired = True
                    specs = (dtaList[i], dtaList[j])
                    lightInd = int(precMass2 < precMass1)
                    specPairs.extend([(ratio, specs[lightInd], specs[1 - lightInd])])
        
        if not paired:
            unpairedSpecs.extend([dtaList[i]])
            if verbose:
                print 'No pairs for', dtaList[i]
    return specPairs, unpairedSpecs

def getAllAAs(seq, ambigAA='X', ambigEdges=None):
    AAs = []
    nodeGen = Constants.nodeInfoGen(seq, considerTerminalMods=True, ambigEdges=ambigEdges)
    for node in nodeGen:
        AAs.extend([node['formAA']])
   
    AAs.extend([node['lattAA']])
    return AAs

def getFilteredPRMLadder(seq, specs, filterRule, ambigAA='X', addEnds=True, ambigEdges=None):
    prmLadder = []
    nodeGen = Constants.nodeInfoGen(seq, considerTerminalMods=True, addTerminalNodes=addEnds, ambigEdges=ambigEdges, ambigAA=ambigAA)
    for node in nodeGen:
        addToLadder = False
        for spec in specs:
            if node['prm'] == 0 or node['prm'] == spec.pm or filterRule(node, spec):
                addToLadder = True
            
        if addToLadder:
            prmLadder += [node['prm']]
    
    return prmLadder

def getPRMLadder(seq, ambigAA='X', addEnds=True, ambigEdges=None):
    prmLadder = []
    nodeGen = Constants.nodeInfoGen(seq, considerTerminalMods=True, addTerminalNodes=addEnds, ambigEdges=ambigEdges, ambigAA=ambigAA)
    for node in nodeGen:
        prmLadder.extend([node['prm']])
    
    return prmLadder

def getAllInds(seq, char):
    return [i for i, x in enumerate(seq) if x == char]

def preprocessSequence(seq, seqMap, replaceExistingTerminalMods=False, ambigAA='X', ambigEdges=None):
    s = list(seq)
    replNTerm = True
    replCTerm = True
    replaceDict = {'Mods':{}, 'AAs':{}}
    for mod in seqMap['Mods']:
        repInds = getAllInds(s, mod)
        if repInds:
            if seqMap['Mods'][mod] in Constants.NTermMods and not replaceExistingTerminalMods:
                replNTerm = False
            if seqMap['Mods'][mod] in Constants.CTermMods and not replaceExistingTerminalMods:
                replCTerm = False
            replaceDict['Mods'][mod] = repInds
    for aa in seqMap['AAs']:
        repInds = getAllInds(s, aa)
        if repInds:
            replaceDict['AAs'][aa] = repInds
    
    for charType in replaceDict:
        for repChar in replaceDict[charType]:
            for ind in replaceDict[charType][repChar]:
                s[ind] = seqMap[charType][repChar]

    s = list(''.join(s))
    AAs = getAllAAs(''.join(s), ambigAA=ambigAA, ambigEdges=ambigEdges)
    if 'N-Term' in seqMap['Mods'] and replNTerm:
        if s[len(AAs[0])] in Constants.NTermMods:
            del s[len(AAs[0])]
        s.insert(len(AAs[0]), seqMap['Mods']['N-Term'])
    
    if 'C-Term' in seqMap['Mods'] and replCTerm:
        if s[-1] in Constants.CTermMods:
            del s[-1]
        s.extend(seqMap['Mods']['C-Term'])

    return ''.join(s)

def getAffinePathGraph(seq1, seq2, AAMap, scoreMatrix, gapOpenPenalty= -5, gapExtendPenalty= -2):
    size = (7, len(seq1) + 2, len(seq2) + 2)
    negInf = -(2 ** 30)
    
    #not really bitArrays, but the closest you can get with numpy (i.e., byte arrays)
    bitArrays = np.zeros(size, dtype=np.dtype('int8'))
    bitArrays[2][-1][-1] = 1
    
    P = np.zeros((size[1] - 1, size[2] - 1), dtype=np.dtype('int'))
    P[0, :] = negInf

    Q = np.zeros((size[1] - 1, size[2] - 1), dtype=np.dtype('int'))
    Q[:, 0] = negInf
    
    R = np.zeros((size[1] - 1, size[2] - 1), dtype=np.dtype('int'))
    R[1:, 0] = [gapOpenPenalty + gapExtendPenalty * i for i in range(1, size[1] - 1)]
    R[0, 1:] = [gapOpenPenalty + gapExtendPenalty * i for i in range(1, size[2] - 1)]
    
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):
            
            if i:
                P[i][j] = gapExtendPenalty + max(P[i - 1][j], R[i - 1][j] + gapOpenPenalty)
                if P[i][j] == P[i - 1][j] + gapExtendPenalty:
                    bitArrays[3][i - 1][j] = 1
                if P[i][j] == R[i - 1][j] + gapExtendPenalty + gapOpenPenalty:
                    bitArrays[4][i - 1][j] = 1
            
            if j:   
                Q[i][j] = gapExtendPenalty + max(Q[i][j - 1], R[i][j - 1] + gapOpenPenalty)
                if Q[i][j] == Q[i][j - 1] + gapExtendPenalty:
                    bitArrays[5][i][j - 1] = 1
                if Q[i][j] == R[i][j - 1] + gapExtendPenalty + gapOpenPenalty:
                    bitArrays[6][i][j - 1] = 1

            if i and j:    
                diagScore = R[i - 1][j - 1] + scoreMatrix[AAMap[seq1[i - 1]]][AAMap[seq2[j - 1]]]
            else:
                diagScore = negInf
            
            if i or j:
                R[i][j] = max(P[i][j], Q[i][j], diagScore)
                if R[i][j] == P[i][j]:
                    bitArrays[0][i][j] = 1
                if R[i][j] == Q[i][j]:
                    bitArrays[1][i][j] = 1
                if R[i][j] == diagScore:
                    bitArrays[2][i][j] = 1

    for i in range(len(seq1) + 1)[::-1]:
        for j in range(len(seq2) + 1)[::-1]:
            if not ((bitArrays[0][i + 1][j] and bitArrays[4][i][j]) or (bitArrays[1][i][j + 1] and bitArrays[6][i][j]) or bitArrays[2][i + 1][j + 1]):
                for k in range(3):
                    bitArrays[k][i][j] = 0

            if bitArrays[0][i + 1][j] or bitArrays[1][i][j + 1] or bitArrays[2][i + 1][j + 1]:
                if bitArrays[0][i + 1][j] and bitArrays[3][i][j]:
                    bitArrays[3][i + 1][j] = 1 - bitArrays[4][i][j]
                    bitArrays[4][i][j] = 1 - bitArrays[0][i][j]
                    bitArrays[0][i][j] = 1
                else:
                    bitArrays[3][i + 1][j] = bitArrays[4][i][j] = 0
            
                if bitArrays[1][i][j + 1] and bitArrays[5][i][j]:
                    bitArrays[5][i][j + 1] = 1 - bitArrays[6][i][j]
                    bitArrays[6][i][j] = 1 - bitArrays[1][i][j]
                    bitArrays[1][i][j] = 1
                else:
                    bitArrays[5][i][j + 1] = bitArrays[6][i][j] = 0

    return R[-1][-1], bitArrays

                
                
def DFSAffinePathGraph(bitArrays):
    queue = deque()
    
    if bitArrays[0][-2][-2]:
        queue.extend([[(bitArrays.shape[1] - 3, bitArrays.shape[2] - 2), (bitArrays.shape[1] - 2, bitArrays.shape[2] - 2)]])
    if bitArrays[1][-2][-2]:
        queue.extend([[(bitArrays.shape[1] - 2, bitArrays.shape[2] - 3), (bitArrays.shape[1] - 2, bitArrays.shape[2] - 2)]])
    if bitArrays[2][-2][-2]:
        queue.extend([[(bitArrays.shape[1] - 3, bitArrays.shape[2] - 3), (bitArrays.shape[1] - 2, bitArrays.shape[2] - 2)]])      
    
    while queue:
        indices = queue.pop()
        yield indices
        
        currInd = indices[0]
        prevInd = indices[1]
        prevEdgeType = tuple(np.array(prevInd) - np.array(currInd))
        
        if bitArrays[0][currInd[0]][currInd[1]]:
            if not (prevEdgeType != (1, 0) and bitArrays[4][currInd[0]][currInd[1]]) and not (prevEdgeType == (0, 1) and bitArrays[5][prevInd[0]][prevInd[1]]):
                queue.extend([[(currInd[0] - 1, currInd[1])] + indices])
        
        if bitArrays[1][currInd[0]][currInd[1]]:
            if not (prevEdgeType != (0, 1) and bitArrays[6][currInd[0]][currInd[1]]) and not (prevEdgeType == (1, 0) and bitArrays[3][prevInd[0]][prevInd[1]]):
                queue.extend([[(currInd[0], currInd[1] - 1)] + indices])
        
        if bitArrays[2][currInd[0]][currInd[1]]:
            if not (prevEdgeType == (0, 1) and bitArrays[5][prevInd[0]][prevInd[1]]) and not (prevEdgeType == (1, 0) and bitArrays[3][prevInd[0]][prevInd[1]]):
                queue.extend([[(currInd[0] - 1, currInd[1] - 1)] + indices])
            
    
    
def alignSequences(seq1, seq2, AAMap, scoreMatrix, gapOpenPenalty= -4, gapExtendPenalty= -2):

    score, bitArrays = getAffinePathGraph(seq1, seq2, AAMap, scoreMatrix, gapOpenPenalty, gapExtendPenalty)
    
    optimalSeqs = []
    for indexSeq in DFSAffinePathGraph(bitArrays):
        if indexSeq[0] == (0, 0):
            optimalSeqs.extend([getSeqsFromIndex(seq1, seq2, indexSeq)])
    
    return score, optimalSeqs
    

def getSeqsFromIndex(seq1, seq2, indices):
    aSeq1 = ''
    aSeq2 = ''
    for i in range(len(indices) - 1):
        if (indices[i + 1][0] - indices[i][0]):
            aSeq1 += seq1[indices[i][0]]
        else:
            aSeq1 += '-'
        
        if (indices[i + 1][1] - indices[i][1]):
            aSeq2 += seq2[indices[i][1]]
        else:
            aSeq2 += '-'
    
    return (aSeq1, aSeq2)
    
def preprocessLADSScanInfo(scanInfo, seqMap, pairConfigurations, optFieldMap={}):
    
    def insertLADSScanDict(scanNum, otherScanNum, seq, score, scan, optFieldMap={}):
        info = {}
        info['Peptide'] = seq
        info['PScore'] = score
        info['Paired Spectrum'] = otherScanNum
        for field in optFieldMap.keys():
            try:
                info[optFieldMap[field]] = scan[field]
            except KeyError:
                info[optFieldMap[field]] = 'None'

        addToInfo = True
        try:
            AAs = getAllAAs(seq, ambigEdges=scan['AMBIGUOUS EDGES'])
        except KeyError:
            addToInfo = False
            
        if addToInfo:    
            try:   
                if (not LADSInfo[scanNum]['Paired Spectrum'] or otherScanNum) and score > LADSInfo[scanNum]['PScore']:
                    LADSInfo[scanNum] = info
            except KeyError:
                LADSInfo[scanNum] = info

    heavySeqMaps = {}
    for confName in pairConfigurations:
        heavySeqMaps[confName] = copy.deepcopy(seqMap)
        heavySeqMaps[confName]['Mods']['N-Term'] = pairConfigurations[confName]['NModSymbol']
        heavySeqMaps[confName]['Mods']['C-Term'] = pairConfigurations[confName]['CModSymbol']
        
    LADSInfo = {}
    for scan in scanInfo:
        try:
            scan['AMBIGUOUS EDGES'] = eval(scan['AMBIGUOUS EDGES'])
        except NameError:
            continue
        lightScan = int(float(scan['LIGHT SCAN']))
        try:
            heavyScan = int(float(scan['HEAVY SCAN']))
        except ValueError:
            heavyScan = False
        
        score = float(scan['SCORE'])
        seq = scan['SEQ']
        try:
            insertLADSScanDict(lightScan, heavyScan, preprocessSequence(seq, seqMap, ambigEdges=scan['AMBIGUOUS EDGES']), score, scan, optFieldMap)
            if heavyScan:
                insertLADSScanDict(heavyScan, lightScan, preprocessSequence(seq, heavySeqMaps[scan['PAIR CONFIGURATION'].lower()], replaceExistingTerminalMods=True, ambigEdges=scan['AMBIGUOUS EDGES']), score, scan, optFieldMap)
        except KeyError:
            print 'KEY ERROR', scan
        
    return LADSInfo
        
def preprocessDatabaseScanInfo(scanInfo, seqMap, fieldMap, scanNumKey='ScanF', srchID=None, seqDelimLen=2):
    DBScanInfo = {}
    for scan in scanInfo:
        #print scan
        scanNum = int(scan[scanNumKey])
        info = {}
        if not (srchID == None or int(scan['SrchID']) == srchID):
            continue

        if 'Ambiguous Edges' in scan:
            scan['Ambiguous Edges'] = eval(scan['Ambiguous Edges'])
            ambigEdges = scan['Ambiguous Edges']
        else:
            ambigEdges = None
            
        try:
            if seqDelimLen > 0:
                procSeq = preprocessSequence(scan['Peptide'][seqDelimLen:-seqDelimLen], seqMap, ambigEdges = ambigEdges)
            else:
                procSeq = preprocessSequence(scan['Peptide'], seqMap, ambigEdges = ambigEdges)
        #except KeyError:
        #    print 'ERROR: Could not process sequence %s with current sequence map' % (scan['Peptide'][seqDelimLen:-seqDelimLen],)
        #    continue
        except IndexError:
            if scan['Peptide'] == 'NULL':
                continue
            else:
                raise IndexError
            
        info['Peptide'] = procSeq
        for field in fieldMap.keys():
            info[fieldMap[field]] = scan[field]

        addToInfo = True
        try:
            AAs = getAllAAs(procSeq, ambigEdges=ambigEdges)
        except KeyError:
            raise KeyError
        DBScanInfo[scanNum] = info
    
    return DBScanInfo

def preprocessPepnovoScanInfo(scanInfo, seqMap, fieldMap, scanNumKey='ScanNum'):
    PNovoScanInfo = {}
    for scan in scanInfo:
        scanNum = int(scan[scanNumKey])
        info = {}
        peptide = scan['RnkScrSequence']
        ambigEdges = []
        PM = float(scan['[M+H]']) - Constants.mods['H2O'] - Constants.mods['H+']
        if float(scan['N-Gap']):
            peptide = 'X' + peptide
            ambigEdges += [(0, float(scan['N-Gap']))]
        if float(scan['C-Gap']):
            peptide = peptide + 'X'
            ambigEdges += [(PM - float(scan['C-Gap']), PM)]


        print peptide, scanNum, 'seq:', scan['RnkScrSequence'], scan['N-Gap'], scan['C-Gap'], scan['RnkScr']
        if len(peptide) == 1:
            continue
        procSeq = preprocessSequence(peptide, seqMap, ambigEdges=ambigEdges, ambigAA='X', replaceExistingTerminalMods=True)
        AAs = getAllAAs(procSeq, ambigEdges=ambigEdges)
        
        for field in fieldMap.keys():
            info[fieldMap[field]] = scan[field]

        info['Peptide'] = procSeq
        info['Ambiguous Edges'] = ambigEdges
        PNovoScanInfo[scanNum] = info

    return PNovoScanInfo
        
                
def comparePeptideResults(seq1, seq2, ambigEdges1=None, ambigEdges2=None, ambigAA='X', ppm=50):
    L1 = np.array(getPRMLadder(seq1, ambigEdges=ambigEdges1, addEnds=False))
    L2 = np.array(getPRMLadder(seq2, ambigEdges=ambigEdges2, addEnds=False))
    epsilon = 1 * 10 ** -6 * L1[-1] * ppm 

    sharedPRMs = getSharedPRMs(L1, L2, epsilon)
    if sharedPRMs:
        accuracy = float(len(sharedPRMs)) / (len(L2))
        precision = float(len(sharedPRMs)) / (len(L1))
        consensus = alignWithPRMs(seq1, seq2, ambigEdges1, ambigEdges2, sharedPRMs, epsilon=epsilon)
        return accuracy, precision, consensus
    else:
        return 0,0, (seq1, seq2)

def getSharedPRMs(prmLadder1, prmLadder2, epsilon=0.5):
    hashTable = {}
    for i in range(prmLadder1.size):
        key = np.round(prmLadder1[i] / epsilon)
        hashTable[key] = [(i, prmLadder1[i])]
    
    temp = np.zeros((prmLadder2.size, 2))
    temp[:, 0] = prmLadder2
    pairedIonData = SA.getPairedIons(hashTable, temp, delta=0.0, epsilon=epsilon)
    sharedPRMs = []
    for key in sorted(pairedIonData.keys()):
        sharedPRMs += [zip(*pairedIonData[key])[1]]

    if sharedPRMs:
        return zip(*sharedPRMs)[0]
    else:
        return []

def getScanComparisonInfo(scansInfo, infoMap, progDict, scanFields=['Score'], compFields=['Precision', 'Accuracy'], specificColDict={}, minNumProgs=2):
    compInfo = {}
    keys = []
    progNames = progDict.keys()
    for field in scanFields:
        keys.extend([progName + ' ' + infoMap[progDict[progName]][field] for progName in progNames])
    
    for cField in compFields:
        for i in range(len(progNames)):
            for j in range(i + 1, len(progNames)):
                try:
                    scansInfo[0][progNames[i] + ' ' + progNames[j] + ' ' + cField]
                    keys.extend([progNames[i] + ' ' + progNames[j] + ' ' + cField])
                except KeyError:
                    keys.extend([progNames[j] + ' ' + progNames[i] + ' ' + cField])
    
    for prog, cols in specificColDict.items():
        for progName, progNprog in progDict.items():
            if progNprog == prog:
                keys.extend([progName + ' ' + col for col in cols])
    
    for scan in scansInfo:
        numProgs = 0
        for progName in progNames:
            if scan[progName + ' ' + infoMap[progDict[progName]]['Peptide']] != 'None':
                numProgs += 1
        
        if numProgs >= minNumProgs:
            scanInfo = {}
            for key in keys:
                try:
                    scanInfo[key] = scan[key]
                except KeyError:
                    scanInfo[key] = None
            compInfo[int(scan['ScanF'])] = scanInfo
    
    return compInfo

def alignWithPRMs(seq1, seq2, ambigEdges1, ambigEdges2, alignedPRMs, epsilon=0.02):
    consensus = [list(seq1), list(seq2)]
    i=0
    j=0
    subMass1=0
    subMass2=0
    for k, prm in enumerate(alignedPRMs):
        while i < len(seq1) and np.abs(subMass1 - prm) > epsilon:
            i += 1
            try:
                subMass1 = getPRMLadder(seq1[:i], ambigEdges=ambigEdges1)[-1]
            except KeyError:
                pass
        
        while j < len(seq2) and np.abs(subMass2 - prm) > epsilon:
            j += 1
            try:
                subMass2 = getPRMLadder(seq2[:j], ambigEdges=ambigEdges2)[-1]
            except KeyError:
                pass
        
        consensus[0].insert(i+k, '|')
        consensus[1].insert(j+k, '|')
    
    return (''.join(consensus[0]), ''.join(consensus[1])) 

def getLADSPScore(seq, dtaPath, PNet, ppm=5, ambigEdges=None, ambigAA='X', ambigPenalty=20):
    pairs = DataFile.getMassIntPairs(dtaPath)
    precMass = DataFile.getPrecMassAndCharge(dtaPath)[0]
    epsilon = ppm * precMass * 10 ** -6
    spec = PN.Spectrum(PNet, precMass, Nmod=0, Cmod=0, epsilon=epsilon, spectrum=pairs)
    spec.initializeNoiseModel()
    nodeGen = Constants.nodeInfoGen(seq, considerTerminalMods=True, ambigEdges=ambigEdges)
    pScore = 0
    node = nodeGen.next()
    pScore += spec.getNodeScore(**node)
    pScore += spec.getPriorScore(prm=0, formAA=None, lattAA=node['formAA'])
    if node['formAA'] == ambigAA:
        pScore -= ambigPenalty
        
    for node in nodeGen:
        pScore += spec.getNodeScore(**node)
        if node['formAA'] == ambigAA:
            pScore -= ambigPenalty
            
    pScore += spec.getPriorScore(prm=precMass- Constants.mods['H+'] - Constants.mods['H2O'], formAA=node['lattAA'], lattAA=None)
    if node['lattAA'] == ambigAA:
        pScore -= ambigPenalty
    
    return pScore  

def findPairsInTrueClusters(trueClusters, pairConfig, results, infoMap, progDict, progName=None, isComp=True, ppm=5):
    pairs = {}
    delta = pairConfig['NMod']+pairConfig['CMod']
    if isComp:
        preceder = progName + ' '
    else:
        preceder = ''
    
    for i in range(len(trueClusters)):
        try:
            PM1 = float(results[list(trueClusters[i])[0]][preceder + infoMap[progDict[progName]]['Obs M+H']])
            seq1 = results[list(trueClusters[i])[0]][preceder + infoMap[progDict[progName]]['Peptide']]
        except ValueError:
            continue
        for j in range(i + 1, len(trueClusters)):
            try:
                PM2 = float(results[list(trueClusters[j])[0]][preceder + infoMap[progDict[progName]]['Obs M+H']])
                seq2 = results[list(trueClusters[j])[0]][preceder + infoMap[progDict[progName]]['Peptide']]
            except ValueError:
                continue
            
            epsilon = ppm * 10 ** -6 * max(PM1, PM2)
            if np.abs(np.abs(PM1 - PM2) - delta) < epsilon:
                pair = [(i, PM1, seq1), (j, PM2, seq2)]
                lightInd = PM2 < PM1
                lAAs = getAllAAs(pair[lightInd][2])
                hAAs = getAllAAs(pair[1 - lightInd][2])
                if lAAs[1:-1] == hAAs[1:-1] and lAAs[0]+pairConfig['NModSymbol'] == pair[1-lightInd][2][:(len(lAAs[0])+len(pairConfig['NModSymbol']))]:
                    pairs[(pair[lightInd][0], pair[1 - lightInd][0])] = 0
    
    return pairs
    
def findPairsInSearchResults(results, infoMap, progDict, pairConfig, progName=None, isComp=True, ppm=5):    
    scans = sorted(results.keys())
    pairs = {}
    delta = pairConfig['NMod']+pairConfig['CMod']
    if isComp:
        preceder = progName + ' '
    else:
        preceder = ''
    for i in range(len(scans)):
        try:
            PM1 = float(results[scans[i]][preceder + infoMap[progDict[progName]]['Obs M+H']])
            seq1 = results[scans[i]][preceder + infoMap[progDict[progName]]['Peptide']]
        except ValueError:
            continue
        for j in range(i + 1, len(scans)):
            try:
                PM2 = float(results[scans[j]][preceder + infoMap[progDict[progName]]['Obs M+H']])
                seq2 = results[scans[j]][preceder + infoMap[progDict[progName]]['Peptide']]
            except ValueError:
                continue
            epsilon = ppm * 10 ** -6 * max(PM1, PM2)
            if np.abs(np.abs(PM1 - PM2) - delta) < epsilon:
                pair = [(scans[i], PM1, seq1), (scans[j], PM2, seq2)]
                lightInd = PM2 < PM1
                lAAs = getAllAAs(pair[lightInd][2])
                hAAs = getAllAAs(pair[1 - lightInd][2])
                if lAAs[1:-1] == hAAs[1:-1] and lAAs[0]+pairConfig['NModSymbol'] == pair[1-lightInd][2][:(len(lAAs[0])+len(pairConfig['NModSymbol']))]:
                    pairs[(pair[lightInd][0], pair[1 - lightInd][0])] = 0
    
    return pairs

def getTrueClusters(results, infoMap, progDict, progName=None, isComp=True, ppm=5):
    pairConfig = {'NMod': 0, 'CMod': 0, 'NStatic': 0, 'CStatic': 0, 'NModSymbol': '', 'CModSymbol': ''}
    pairs = findPairsInSearchResults(results, infoMap, progDict, pairConfig, progName=progName, isComp=isComp, ppm=ppm)

    trueClusters = [set(pair) for pair in pairs]
    prevNumClusters = len(trueClusters)
    curNumClusters = 0
    while prevNumClusters != curNumClusters:
        newTrueClusters = []
        mergedInds = set()
        for i in range(len(trueClusters)):
            if i not in mergedInds:
                mergedCluster = trueClusters[i]
                for j in range(i+1, len(trueClusters)):
                    if len(mergedCluster & trueClusters[j]) > 0:
                        mergedCluster = mergedCluster | trueClusters[j]
                        mergedInds.add(j)
                        
                newTrueClusters += [mergedCluster]
                        
        prevNumClusters = curNumClusters
        curNumClusters = len(newTrueClusters)
        trueClusters = newTrueClusters

    """
    trueClusters = []
    for pair in pairs:
        addedToCluster = False
        print pair
        for cluster in trueClusters:
            if pair[0] in cluster or pair[1] in cluster:
                cluster.add(pair[0])
                cluster.add(pair[1])
                addedToCluster = True
                break
       
        prevNumClusters = curNumClusters
        curNumClusters = len(newTrueClusters)
        trueClusters = newTrueClusters
            
    """     
#    print trueClusters
    allClusteredScanFs = set()
    for cluster in trueClusters:
        allClusteredScanFs = allClusteredScanFs | cluster
    notClusteredScanFs = set(results.keys()) - allClusteredScanFs

    
    for scanF in notClusteredScanFs:
        trueClusters += [set([scanF])]
    
    return trueClusters

def checkPairsDict(pairsDict, scanF1, scanF2, pairConfigName):
    inPairsDict = False
    pairConfigName = pairConfigName.lower()
    if pairConfigName == 'n/a':
        return None
    for progName in pairsDict[pairConfigName].keys():
        if (scanF1, scanF2) in pairsDict[pairConfigName][progName]:
            pairsDict[pairConfigName][progName][(scanF1, scanF2)] = 1
            inPairsDict = True
        elif (scanF2, scanF1) in pairsDict[pairConfigName][progName]:
            pairsDict[pairConfigName][progName][(scanF2, scanF1)] = 1
            inPairsDict = True
    
    return inPairsDict

def determinePairType(pairsDict, scanInfo, progDict, infoMap, LADSprogName):
    scanF1 = int(scanInfo['ScanF'])
    try:
        scanInfo[LADSprogName + ' Paired Spectrum']
    except KeyError:
        return None
    
    hasInfo = False
    for progName in progDict.keys():
        if progName == LADSprogName:
            continue
        if scanInfo[progName + ' ' + infoMap[progDict[progName]]['Score']] != 'None':
            hasInfo = True
    
    if not hasInfo:
        return None


    try:
        return checkPairsDict(pairsDict, scanF1, int(scanInfo[LADSprogName + ' Paired Spectrum']), scanInfo[LADSprogName + ' Pair Configuration'])
    except ValueError:
        return None
        
def getAccuracyAndPrecisionNames(prog1, prog2, searchRow):
    try:
        searchRow[prog1 + ' ' + prog2 + ' Accuracy']
        accName = prog1 + ' ' + prog2 + ' Accuracy'
        precName = prog1 + ' ' + prog2 + ' Precision'
    except KeyError:
        try:
            searchRow[prog2 + ' ' + prog1 + ' Accuracy']
            accName = prog2 + ' ' + prog1 + ' Precision'
            precName = prog2 + ' ' + prog1 + ' Accuracy'
        except KeyError:
            raise KeyError('ERROR: Accuracy comparison between %s and %s not present in CompSearch output' % (prog1, prog2))
    
    return accName, precName
    

def getCompStats(compSearchPath, mainProgName, progDict, infoMap, paramsDict, mainProgFields=['PScore', 'Num Ambig Edges'], getPairStats=True):
    compSearchInfo = DataFile.getScanInfo(compSearchPath, delimiter='\t')
    unpaired = {}
    other = {}
    stats = {}

    for progName, prog in progDict.items():
        if progName == mainProgName:
            continue
        
        unpaired[progName] = {'accuracyVec': np.array([]), 'precisionVec': np.array([]), 'numScans': 0}
        accName, precName = getAccuracyAndPrecisionNames(progName, mainProgName, compSearchInfo[0])
        stats[progName] = {}
        stats[progName]['accName'] = accName
        stats[progName]['precName'] = precName
        for progfield in mainProgFields:
            unpaired[progName][progfield] = np.array([])
        other[progName] = copy.deepcopy(unpaired[progName])

    pairsDict = {}
    if getPairStats:
        truePairs = {}
        falsePairs = {}
        compInfo = getScanComparisonInfo(compSearchInfo, infoMap, progDict, scanFields=['Score', 'Peptide', 'Obs M+H'], specificColDict={'LADS': ['Num Ambig Edges', 'Paired Spectrum', 'Pair Configuration']})
        for pairConfigName in paramsDict['Pair Configurations']:
            truePairs[pairConfigName] = {}
            falsePairs[pairConfigName] = {}
            pairsDict[pairConfigName] = {}
            for progName in progDict.keys():
                if progName == mainProgName:
                    continue
                pairsDict[pairConfigName][progName] = findPairsInSearchResults(compInfo, infoMap, progDict, paramsDict['Pair Configurations'][pairConfigName], progName=progName)
                truePairs[pairConfigName][progName] = copy.deepcopy(unpaired[progName])
                falsePairs[pairConfigName][progName] = copy.deepcopy(unpaired[progName])

    print 'Compiling stats'
    for scan in compSearchInfo:
        scanF1 = int(scan['ScanF'])
        pairType = determinePairType(pairsDict, scan, progDict, infoMap, mainProgName)
        if pairType == None:
            temp = unpaired
        elif pairType:
            temp = truePairs[scan[mainProgName + ' Pair Configuration'].lower()]
        else:
            temp = falsePairs[scan[mainProgName + ' Pair Configuration'].lower()]
        for progName in stats.keys():
            try:
                if scan[progName + ' ' + infoMap[progDict[progName]]['Score']] != 'None':
                    temp[progName]['numScans'] += 1
                temp[progName]['accuracyVec'] = np.append(temp[progName]['accuracyVec'], float(scan[stats[progName]['accName']]))
                temp[progName]['precisionVec'] = np.append(temp[progName]['precisionVec'], float(scan[stats[progName]['precName']]))
                for progfield in mainProgFields:
                    temp[progName][progfield] = np.append(temp[progName][progfield], float(scan[mainProgName + ' ' + progfield]))
            except ValueError:
                other[progName]['numScans'] += 1
                for progfield in mainProgFields:
                    try:
                        other[progName][progfield] = np.append(other[progName][progfield], float(scan[mainProgName + ' ' + progfield]))
                    except ValueError:
                        print 'ERROR in getting main %s data for scan %s, peptide %s, %s %s' % (mainProgName, scan['ScanF'], scan[mainProgName + ' ' + infoMap[progDict[mainProgName]]['Peptide']], progfield, scan[mainProgName + ' ' + progfield])
                        pass
    
    for progName in stats.keys():
        if getPairStats:
            stats[progName]['truepairs'] = {}
            stats[progName]['falsepairs'] = {}
            stats[progName]['pairsDict'] = {}
            stats[progName]['unpaired'] = unpaired[progName]
            stats[progName]['other'] = other[progName]
            stats[progName]['composite'] = {}
            for pairConfigName in truePairs:
                stats[progName]['truepairs'][pairConfigName] = truePairs[pairConfigName][progName]
                stats[progName]['falsepairs'][pairConfigName] = falsePairs[pairConfigName][progName]
                stats[progName]['pairsDict'][pairConfigName] = pairsDict[pairConfigName][progName]
                
            for field in stats[progName]['unpaired']:
                try:
                    truePairsComp = np.concatenate([stats[progName]['truepairs'][pairConfigName][field] for pairConfigName in stats[progName]['truepairs']])
                    falsePairsComp = np.concatenate([stats[progName]['falsepairs'][pairConfigName][field] for pairConfigName in stats[progName]['falsepairs']])
                    stats[progName]['composite'][field] = np.concatenate((truePairsComp, falsePairsComp, stats[progName]['unpaired'][field]))
                except ValueError:
                    pass
            
            numTruePairs = np.sum([stats[progName]['truepairs'][pairConfigName]['numScans'] for pairConfigName in stats[progName]['truepairs']])
            numFalsePairs = np.sum([stats[progName]['falsepairs'][pairConfigName]['numScans'] for pairConfigName in stats[progName]['falsepairs']])
            stats[progName]['composite']['numScans'] = numTruePairs + numFalsePairs + stats[progName]['unpaired']['numScans']
        else:
            stats[progName]['other'] = other[progName]
            stats[progName]['composite'] = unpaired[progName]

    return stats
        
if __name__ == '__main__':
    #scanInfoF = 'adevabhaktuni_1310166306.csv'
    #dirPath = '//lab//user_data//adevabhaktuni//LF2_short_HCD+CID_ath001862_244//'
    #print getAlignmentRatios(scanInfoF, dirPath, delta=Constants.mods['*'])
    """
    seq2 = 'HLELNHLELNPDAAAA'
    seq1 = 'HLELNLLLLLLLLLLLLLLLLLLLLLLLLLPDEHATK'
    path = 'BLOSUM90.txt'
    
    seqMap = {}
    for aa in Constants.aminoacids.keys():
        seqMap[aa] = aa
    
    seqMap['I'] = 'I'
    seqMap['L'] = 'I'
    seqMap['J'] = 'I'
    seqMap['X'] = 'J'
    seqMap['K*'] = 'J'
    seqMap['M#'] = 'O'
    seqMap['O'] = 'O'
    seqMap['-'] = 'X'
    dir = 'C:\\Users\\Arun\\Proteomics Files\\'
    AAMap, scoreMatrix = DataFile.getScoringMatrix(path)
    #print alignSequences(seq1, seq2, AAMap, scoreMatrix, gapOpenPenalty=-5, gapExtendPenalty=-2)
    """
    """
    LADSfields = ['LIGHT SCAN', 'HEAVY SCAN', 'SEQ', 'SCORE', 'AMBIGUOUS EDGES', 'M+H', 'EPSILON', 'SHARED PEAKS RATIO', 'SEQUENCING TIME', 'NUM PATHS', 'AA EDGES']
    LADSScanInfo = DataFile.getScanInfo(dir + 'ath001862UPen10KPen15LRRestrictTest.tdv', LADSfields, delimiter='\t')
    optFieldMap = {'AMBIGUOUS EDGES': 'Ambiguous Edges', 'M+H': 'Obs M+H', 'EPSILON': 'Epsilon', 'SHARED PEAKS RATIO': 'Shared Peaks Ratio', 'SEQUENCING TIME': 'Sequencing Time'}
    LADSInfo = preprocessLADSScanInfo(LADSScanInfo, seqMap, optFieldMap)
    for key in sorted(LADSInfo.keys()):
        print key, LADSInfo[key]
    """
    #MASCOTFields = ['ScanF', 'Obs M+H', 'PPM', 'Ion Score', 'h-factor', 'Reference', 'Peptide']
    #MASCOTData = DataFile.getScanInfo(dir + 'LF2_short_CID+HCD_ath001862_MASCOT.csv', MASCOTFields, delimiter=',')
    #MASCOTFieldMap = {'Obs M+H': 'Obs M+H', 'Ion Score': 'Ion Score', 'h-factor': 'h-factor', 'PPM': 'PPM', 'Reference': 'Reference'}
    
    #MASCOTInfo = preprocessDatabaseScanInfo(MASCOTData, seqMap, MASCOTFieldMap)
    #for key in sorted(MASCOTInfo.keys()):
    #    print key, MASCOTInfo[key]
    """
    seq = '-K*M#--'
    ambigEdges=[(0,1000),(0,2000),(0,4000)]
    paramsDict = DataFile.parseParams('./Misc/LADS_SILAC_Trypsin.ini')
    seqMap = DataFile.generateSeqMap(['SEQUEST', 'MASCOT', 'LADS'], paramsDict)
    Constants.NTermMods['['] = 10000
    Constants.CTermMods[']'] = 20000
    Constants.NTermMods['^'] = 40000
    seqMap['LADS']['Mods']['N-Term'] = '['
    seqMap['LADS']['Mods']['$'] = '^'
    seqMap['LADS']['Mods']['C-Term'] = ']'
    #nodeGen = Constants.nodeInfoGen(seq, considerTerminalMods=True, addTerminalNodes=True, ambigEdges=ambigEdges)
    #for node in nodeGen:
        #print node
    newSeq = preprocessSequence(seq, seqMap['LADS'], replaceExistingTerminalMods=True, ambigEdges=copy.copy(ambigEdges))
    print newSeq
    print getPRMLadder(newSeq, ambigAA='X', addEnds=True, ambigEdges=copy.copy(ambigEdges))
    """
    """
    #print comparePeptideResults('AAKKIKK', 'KAAIKKK')
    dirPath = '/home/arun/Proteomics_Data/LF2_short_HCD+CID_ath001862_244'
    heavyPath = dirPath + '244.3383.3383.1.dta'
    lightPath = dirPath + '3760.0160.0160.3.dta'
    """
    dbInfo = DataFile.getDBInfo('./Misc/searchformatparams.pag')
    compInfo = DataFile.getScanInfo('/home/arun/Downloads/ath009552_ppm5_alpha0.90_min300_max500_PC0.05_SGCut300_SEQUEST_CompareSearches.py_UnitTest.tdv', delimiter='\t')
    writeFASTAFile(compInfo, 'LADS Unit Test', dbInfo['infoMap'], 'test.txt', ambigEdgeCutoff=1)

    
