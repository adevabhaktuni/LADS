'''
Created on Sep 26, 2011

@author: Arun
'''
import os
import sys

sys.path.insert(1, os.path.abspath(os.pardir))


import glob
import numpy as np
import pickle
import sys
from collections import defaultdict

import Analytics as An
import ArgLib
import DataFile
import CompareSearches as CS
import SpectraAlignment as SA

def reverseDict(dictToReverse):
    return dict((dictToReverse[i], i) for i in dictToReverse)

#both keys and values must be tuples
def addDicts(dict1, dict2):
    combDict = copy.deepcopy(dict1)
    for key in dict2:
        if key in dict1:
            combDict[key] = dict1[key] + dict2[key]
        else:
            combDict[key] = dict2[key]
            
    return combDict

def prepIonTableForAddition(ionTable, pairType):
    Inds = []
    for key in ionTable:
        Inds += [zip(*ionTable[key])[0]]
    ionDict = dict((((IndPair[0], pairType[0]),),((IndPair[1], pairType[1]),)) for IndPair in Inds)
    return ionDict
                                                                        
def getFunction(funcname):
    try:
        modname, funcname = funcname.split('.')
        mod = __import__(modname)
        return getattr(mod, funcname)
    except ValueError:
        return getattr(An, funcname)

def countNumberUniquePeptides(results, searchProg, infoMap):
    seqDict = {}
    for scanF in results:
        seq = results[scanF][infoMap[searchProg]['Peptide']]
        if seq in seqDict:
            seqDict[seq] += 1
        else:
            seqDict[seq] = 1

    print 'Number unique peptides', len(seqDict)
    print 'Number total peptides', sum(seqDict[seq] for seq in seqDict)

def getScanFDiffVec(pairs):
    vec = []
    for pair in pairs:
        vec += [pair[0] - pair[1]]

    return np.array(vec)

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = DataFile.getScanNum(dta)
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'sequenced': False}
        
    return scanFDict
                            
def getSharedPeaksRatio(lightPairs, heavyPairs, pairConfig, epsilon):
    N, C = SA.getNandCIons(lightPairs, heavyPairs, pairConfig['NMod'], pairConfig['CMod'], epsilon=epsilon)
    return SA.getSharedPeaksRatio(lightPairs, heavyPairs, N, C)

def getClusterRecallAndPrecision(results, testClusters, trueClusters, scanFDict, hashEp=0.02):
    hashedTrueClusters = {}
    for cluster in trueClusters:
        hMass = np.round(scanFDict[list(cluster)[0]]['precMass']/hashEp)
        for ind in hMass-1, hMass, hMass+1:
            if ind in hashedTrueClusters:
                hashedTrueClusters[ind] += [cluster]
            else:
                hashedTrueClusters[ind] = [cluster]

    recSum, precSum, nonSingRecSum, nonSingPrecSum = 0, 0, 0, 0
    numRelevantClusters = 0
    numNonSingularClusters = 0
#    restClusterScanFs = set()
    for cluster in testClusters:
        restCluster = set()
        for scanF in cluster:
            if scanF in results:
                restCluster.add(scanF)
        try:
            trueClusters = hashedTrueClusters[np.round(scanFDict[cluster[0]]['precMass']/hashEp)]
        except KeyError:
            continue
        
        if len(restCluster) > 0:
#            print 'testCluster', cluster
#            print 'restCluster', restCluster
#            restClusterScanFs = restClusterScanFs | restCluster
            numRelevantClusters += 1
            recSum += max([len(restCluster & trueCluster)/float(len(trueCluster)) for trueCluster in trueClusters])
            precSum += max([len(restCluster & trueCluster)/float(len(restCluster)) for trueCluster in trueClusters])

        # Now look at clusters which have more than one element (gets rid of score boost for correctly guessing clusters with one element in them)
        if len(restCluster) > 1:
            numNonSingularClusters += 1
            nonSingRecSum += max([len(restCluster & trueCluster)/float(len(trueCluster)) for trueCluster in trueClusters])
            nonSingPrecSum += max([len(restCluster & trueCluster)/float(len(restCluster)) for trueCluster in trueClusters])
#    print 'number relevant clusters', numRelevantClusters
#    print 'number relevant scanFs in test clusters', len(restClusterScanFs)
    return recSum/numRelevantClusters, precSum/numRelevantClusters, nonSingRecSum/numNonSingularClusters, nonSingPrecSum/numNonSingularClusters

def getClusterPairwiseStats(results, precMassClusters, testClusters, trueClusters):
    trueClusterDict = {}
    for cluster in trueClusters:
        lCluster = list(cluster)
        for i in range(len(lCluster)):
            for j in range(i+1, len(lCluster)):
                trueClusterDict[(lCluster[i], lCluster[j])] = 0

    numPrecMassClusterPairs = 0
    for cluster in precMassClusters:
        numRelScanFs = sum([int(scanF in results) for scanF in cluster])
        numPrecMassClusterPairs += numRelScanFs * (numRelScanFs - 1) / 2

    print len(trueClusterDict), numPrecMassClusterPairs
    TP, FP, TN, FN = 0, 0, 0, 0
    for cluster in testClusters:
        restCluster = []
        for scanF in cluster:
            if scanF in results:
                restCluster += [scanF]

        for i in range(len(restCluster)):
            for j in range(i+1, len(restCluster)):
                if (restCluster[i], restCluster[j]) in trueClusterDict or (restCluster[j], restCluster[i]) in trueClusterDict:
                    TP += 1
                else:
                    print 'FP', restCluster[i], restCluster[j]
                    print results[restCluster[i]], results[restCluster[j]]
                    FP += 1

    FN = len(trueClusterDict) - TP
    TN = numPrecMassClusterPairs - FN - TP - FP

    return TP, FP, TN, FN
        

def writePairsDistributionInfo(outFile, truePairs, testPairs, name, numBins=20):
    truePairsVec = np.array(getScanFDiffVec(truePairs), dtype=np.float64)
    testPairsVec = np.array(getScanFDiffVec(testPairs), dtype=np.float64)
    print 'truePairs: ', truePairsVec
    print 'testPairs: ', testPairsVec
    minDiff = testPairsVec.min()
    maxDiff = testPairsVec.max()
    print 'testpairs: ', minDiff, maxDiff, 'truepairs: ', truePairsVec.min(), truePairsVec.max()
    bins = np.floor((np.arange(numBins) * (maxDiff-minDiff)/numBins + minDiff)).astype(int)
    binsDict = {}
    for i in range(numBins):
        binsDict[i] = [0,0,0]

    truePairsVec = np.floor((truePairsVec-minDiff)/(maxDiff-minDiff) * numBins).astype(int)
    testPairsVec = np.floor((testPairsVec-minDiff)/(maxDiff-minDiff) * numBins).astype(int)
    for bin in testPairsVec:
        if bin >= numBins:
            binsDict[numBins-1][0] += 1
        else:
            binsDict[bin][0] += 1

    for bin in truePairsVec:
        if bin >= numBins:
            binsDict[numBins-1][1] += 1
        else:
            binsDict[bin][1] += 1

    for bin in binsDict:
        binsDict[bin][2] = binsDict[bin][0] - binsDict[bin][1]

    outFile.write('\n%s Scan Number Difference Distribution. Max Diff: %i' % (name, maxDiff) + '\n')
    outFile.write('\t'.join(['Diff Bin', 'Test Pairs', 'True Pairs', 'False Pairs']) + '\n')
    for i in range(numBins):
        outFile.write('\t'.join([str(elem) for elem in [bins[i], binsDict[i][0], binsDict[i][1], binsDict[i][2]]]) + '\n')
    
    
if __name__ == '__main__':
    print 'Model refers to svmmodel used'
    options = ArgLib.parse(['dtadir', 'combined', 'sequest', 'mascot', 'database', 'output', 'ppmstd', 'init', 'symbolmap'])
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    progDict = ArgLib.getProgDict(An.searchprogs, options)

    dbDict = DataFile.getDBInfo(options.database)
    infoMap = dbDict['infoMap']
    
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    
    seqMap = DataFile.generateSeqMap(progDict, symbolMap, paramsDict)

    processedInfo = {} 
    
    if options.mascot:
        MASCOTdict = eval(options.mascot)
        processedInfo.update(CS.parseScans(MASCOTdict, 'MASCOT', seqMap, dbDict))     
    if options.sequest:
        SEQUESTdict = eval(options.sequest)
        processedInfo.update(CS.parseScans(SEQUESTdict, 'SEQUEST', seqMap, dbDict))
    if options.combined:
        combinedDict = eval(options.combined)
        processedInfo.update(CS.parseScans(combinedDict, 'Combined', seqMap, dbDict, delimiter='\t', seqDelimLen=0))
    
    if len(processedInfo.keys()) > 1:
        print 'ERROR: Can only compare results to one database search output. Exiting...'
        sys.exit(-1)
    
    progName = processedInfo.keys()[0]

    dtaList = glob.glob(options.dtadir + '/*.dta')
    scanFDict = getScanFDict(dtaList)
    
    precMassClusters = An.findSamePrecMassClusters(dtaList, ppm=options.ppmstd)

    clusterOut = open(options.output + '_cluster.txt', 'w')

    
    for cluster in precMassClusters:
        if len(cluster) == 1:
            continue

        specs = []
        for scanF in cluster:
            specs += [DataFile.getMassIntPairs(scanFDict[scanF]['dta'])]

        for i in range(len(cluster)):
            for j in range(i+1, len(cluster)):
                if cluster[i] in processedInfo[progName] and cluster[j] in processedInfo[progName]:
                    epSTD = options.ppmstd * 10 ** -6 * scanFDict[cluster[i]]['precMass']
                
                    SVMClassificationInfo = SA.getSpectraPairInfoForSVMClassification(specs[i], specs[j], scanFDict[cluster[i]]['precMass'], NMod=0, CMod=0, epsilon=2*epSTD)
                    seq1 = processedInfo[progName][cluster[i]][infoMap[progDict[progName]]['Peptide']]
                    seq2 = processedInfo[progName][cluster[j]][infoMap[progDict[progName]]['Peptide']]

                    xVal = 1 if seq1 == seq2 else -1
                    clusterOut.write(' '.join([str(xVal)] + ['%i:%f' % (key, SVMClassificationInfo[key]) for key in sorted(SVMClassificationInfo)]) + ' # Scans %s, %i - %s, %i\n' % (processedInfo[progName][cluster[i]][infoMap[progDict[progName]]['Peptide']], cluster[i], processedInfo[progName][cluster[j]][infoMap[progDict[progName]]['Peptide']], cluster[j]))

    clusterOut.close()

    
    for pairConfigName in paramsDict['Pair Configurations']:
        pairConfig = paramsDict['Pair Configurations'][pairConfigName]

        delta = pairConfig['NMod'] + pairConfig['CMod']
        deltaPairs = An.findDeltaPairsClusters(precMassClusters, scanFDict, delta, ppm=options.ppmstd)

        pairsOut = open(options.output + '_' + pairConfigName + '.txt', 'w')

        for pair in deltaPairs:

            epSTD = options.ppmstd * 10 ** -6 * scanFDict[precMassClusters[pair[0]][0]]['precMass']
            
            # Get all possible true pairings from database search results
            uniquePeptideDict = defaultdict(lambda: {'light': [], 'heavy': []})
            for scanF in precMassClusters[pair[0]]:
                if scanF in processedInfo[progName]:
                    uniquePeptideDict[An.stripModifications(processedInfo[progName][scanF][infoMap[progDict[progName]]['Peptide']], noRemove=['#'])]['light'] += [scanF]

            for scanF in precMassClusters[pair[1]]:
                if scanF in processedInfo[progName]:
                    uniquePeptideDict[An.stripModifications(processedInfo[progName][scanF][infoMap[progDict[progName]]['Peptide']], noRemove=['#'])]['heavy'] += [scanF]


            mergedSpecDictLight = {}
            mergedSpecDictHeavy = {}
            for peptide in uniquePeptideDict:
                if len(uniquePeptideDict[peptide]['light']) > 0:
                    mergedSpecDictLight[peptide] = SA.mergeSpectra([DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in uniquePeptideDict[peptide]['light']], epsilon = 2*epSTD)
                if len(uniquePeptideDict[peptide]['heavy']) > 0:
                    mergedSpecDictHeavy[peptide] = SA.mergeSpectra([DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in uniquePeptideDict[peptide]['heavy']], epsilon = 2*epSTD)                

            for lightPept in mergedSpecDictLight:
                for heavyPept in mergedSpecDictHeavy:

                    xVal = 1 if lightPept == heavyPept else -1
                    SVMClassificationInfo = SA.getSpectraPairInfoForSVMClassification(mergedSpecDictLight[lightPept], mergedSpecDictHeavy[heavyPept], scanFDict[precMassClusters[pair[0]][0]]['precMass'], pairConfig['NMod'], pairConfig['CMod'], epsilon=2*epSTD)

                    pairsOut.write(' '.join([str(xVal)] + ['%i:%f' % (key, SVMClassificationInfo[key]) for key in sorted(SVMClassificationInfo)]) + ' # light peptide: %s, %s - heavy scans: %s, %s\n' % (lightPept, str(uniquePeptideDict[lightPept]['light']), heavyPept, str(uniquePeptideDict[heavyPept]['heavy'])))


        pairsOut.close()
                                                                                                                                 

