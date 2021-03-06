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
import time
import svmutil
import os
import copy

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
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    
    seqMap = DataFile.generateSeqMap(progDict, symbolMap, paramsDict)
    outFile = open(options.output, 'w')

    print options.dtadir
    dtaList = glob.glob(options.dtadir + '/*.dta')
    scanFDict = getScanFDict(dtaList)

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

    outFile.write('Scan information fetched. Total number of scans: %i. Number of scans considered for validation: %i' % (len(scanFDict), len(processedInfo[progName])))

    progPairs = {}
    for pairConfigName in paramsDict['Pair Configurations']:
        progPairs[pairConfigName] = An.findPairsInSearchResults(processedInfo[progName], dbDict['infoMap'], progDict, paramsDict['Pair Configurations'][pairConfigName], progName=progName, isComp=False, ppm=options.ppmstd)

    pairs = {}
    times = {}
    
    t1 = time.time()
    print 'Getting Clusters'
    clusterSVMModel = svmutil.svm_load_model(paramsDict['Cluster Configuration']['model'])
    clusterSVMRanges = svmutil.load_ranges(os.path.splitext((paramsDict['Cluster Configuration']['model']))[0] + '.range')
    
    precMassClusters = An.findSamePrecMassClusters(dtaList, ppm=options.ppmstd)
#    print 'precMassClusters', precMassClusters
    samePeptideClusters = An.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=float(paramsDict['Cluster Configuration']['cutoff']))
#    samePeptideClusters = An.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=4)
    outFile.write('\n%i Clusters Found. Time taken: %f\n' % (len(samePeptideClusters), time.time()-t1))
    trueClusters = An.getTrueClusters(processedInfo[progName], dbDict['infoMap'], progDict, progName=progName, isComp=False, ppm=options.ppmstd)
    """
    clusteredScanFs = set()
    nonUniqueClusteredScanFs = []
    for cluster in trueClusters:
        if len(cluster & clusteredScanFs) > 0:
            print 'ERROR: Duplicate ScanFs in true clusters at cluster', cluster
        clusteredScanFs = clusteredScanFs | cluster
        nonUniqueClusteredScanFs += list(cluster)
    print 'number of unique scanFs in true clusters: ', len(clusteredScanFs)
    print 'number of total scanFs in true clusters: ', len(nonUniqueClusteredScanFs)
    """
    clustRecall, clustPrecision, nonSingClustRecall, nonSingClustPrecision = getClusterRecallAndPrecision(processedInfo[progName], samePeptideClusters, trueClusters, scanFDict)
    """
    infoMap = dbDict['infoMap']
    for cluster in trueClusters:
        lClust = list(cluster)
        for i in range(len(cluster)):
            for j in range(i+1, len(cluster)):
                seq1 = processedInfo[progName][lClust[i]][infoMap[progDict[progName]]['Peptide']]
                seq2 = processedInfo[progName][lClust[j]][infoMap[progDict[progName]]['Peptide']]
                if seq1 != seq2:
                    print 'ERROR', seq1, seq2, 'do not match for', lClust[i], lClust[j], 'in cluster', cluster
    print 'Finished testing true clusters'
    countNumberUniquePeptides(processedInfo[progName], progDict[progName], infoMap)
    """            
    outFile.write('\nCluster Statistics\n')
    outFile.write('\t'.join(['Type', 'Number', 'Max', 'Min', 'Precision', 'Recall', 'NonSing_Prec', 'NonSing_Rec']) + '\n')
    outFile.write('\t'.join([str(attr) for attr in ['True', len(trueClusters), max([len(cluster) for cluster in trueClusters]), min([len(cluster) for cluster in trueClusters]), 'N/A', 'N/A', 'N/A', 'N/A']]) + '\n')
    outFile.write('\t'.join([str(attr) for attr in ['Test', len(samePeptideClusters), max([len(cluster) for cluster in samePeptideClusters]), min([len(cluster) for cluster in samePeptideClusters]), clustPrecision, clustRecall, nonSingClustPrecision, nonSingClustRecall]]) + '\n')

    # Measure pairwise stats for cluster quality (similar to assessing pairs)
    TP, FP, TN, FN = getClusterPairwiseStats(processedInfo[progName], precMassClusters, samePeptideClusters, trueClusters)
    try:
        sensitivity = int(float(TP)/(TP + FN) * 10000)/float(100)
    except ZeroDivisionError:
        sensitivity = 100
        
    try:
        precision = int(float(TP)/(TP + FP) * 10000)/float(100)
    except ZeroDivisionError:
        precision = 100
    outFile.write('\n' + '\t'.join(['TP', 'FP', 'TN', 'FN', 'Sensitivity', 'Precision']) + '\n')
    outFile.write('\t'.join([str(val) for val in [TP, FP, TN, FN, sensitivity, precision]]) + '\n')
    
    for pairConfigName in paramsDict['Pair Configurations']:
        pairs[pairConfigName] = {}
        pairConfig = paramsDict['Pair Configurations'][pairConfigName]
        startTime = time.time()

        svmModel = svmutil.svm_load_model(pairConfig['Model'])
        svmRange = svmutil.load_ranges(os.path.splitext(pairConfig['Model'])[0] + '.range')

        delta = pairConfig['NMod'] + pairConfig['CMod']
        deltaPairs = An.findDeltaPairsClusters(samePeptideClusters, scanFDict, delta, ppm=options.ppmstd)
        outFile.write('\nTotal number of cluster pairs considered for pair Config %s (including pairs not reported by database search results): %i\n' % (pairConfigName, len(deltaPairs)))
        x, y = [], []
 #       testedDeltaPairedScanFs = set()
 #       testedPairs = set()
 #       possPairsList = []
        testedDeltaPairs = []
        for pair in deltaPairs:
            if not (any([scanF in processedInfo[progName] for scanF in samePeptideClusters[pair[0]]]) and any([scanF in processedInfo[progName] for scanF in samePeptideClusters[pair[1]]])):
                continue
            print 'getting data for pair', pair
            #testedPairs += [pair]
            testedDeltaPairs += [pair]
#            for scanF in list(samePeptideClusters[pair[0]]) + list(samePeptideClusters[pair[1]]):
#                testedDeltaPairedScanFs.add(scanF)
#            for possPair in [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[pair[0]] for heavyScanF in samePeptideClusters[pair[1]]]:
#                testedPairs.add(possPair)

            lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
            heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
            lightPrecMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))
 #           heavyPrecMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[1]]]))
 #           print lightPrecMass, heavyPrecMass
 #           print samePeptideClusters[pair[0]], samePeptideClusters[pair[1]]
            epSTD = options.ppmstd * 10 ** -6 * lightPrecMass

            lightMergedSpec = SA.mergeSpectra(lightSpecs, epsilon=2*epSTD)
            heavyMergedSpec = SA.mergeSpectra(heavySpecs, epsilon=2*epSTD)

            """
            NTermTable, CTermTable = SA.getNandCIons(lightMergedSpec, heavyMergedSpec, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD)
            NCrossTable, CCrossTable = SA.getCrossPairedIons(lightMergedSpec, heavyMergedSpec, lightPrecMass, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD)

            NTermIonDict = prepIonTableForAddition(NTermTable, ['b', 'b'])
            CTermIonDict = prepIonTableForAddition(CTermTable, ['y', 'y'])
            NCrossIonDict = prepIonTableForAddition(NCrossTable, ['y', 'b'])
            CCrossIonDict = prepIonTableForAddition(CCrossTable, ['b', 'y'])
            
            allPairedIonsDict = addDicts(reverseDict(addDicts(NTermIonDict, CCrossIonDict)), reverseDict(addDicts(NCrossIonDict, CTermIonDict)))
            symLightInds = set()
            symHeavyInds = set()
            totalLightInds = set()
            totalHeavyInds = set()
            doubleSymLightInds = set()
            doubleSymHeavyInds = set()
            
            for heavyIons in allPairedIonsDict:
                if len(heavyIons) == 2 and len(allPairedIonsDict[heavyIons]) == 2:
                    for ion in heavyIons:
                        doubleSymHeavyInds.add(ion[0])
                    for ion in allPairedIonsDict[heavyIons]:
                        doubleSymLightInds.add(ion[0])
                if len(heavyIons) == 2 or len(allPairedIonsDict[heavyIons]) ==2:
                    for ion in heavyIons:
                        symHeavyInds.add(ion[0])
                    for ion in allPairedIonsDict[heavyIons]:
                        symLightInds.add(ion[0])
                        
                for ion in heavyIons:
                    totalHeavyInds.add(ion[0])
                for ion in allPairedIonsDict[heavyIons]:
                    totalLightInds.add(ion[0])

            totalNumPeaks = float(lightMergedSpec.shape[0] + heavyMergedSpec.shape[0])
            totalSharedPeaksRatio = (len(totalLightInds) + len(totalHeavyInds))/totalNumPeaks
            singleSymSharedPeaksRatio = (len(symLightInds) + len(symHeavyInds))/totalNumPeaks
            doubleSymSharedPeaksRatio = (len(doubleSymLightInds) + len(doubleSymHeavyInds))/totalNumPeaks
            """
            SVMClassificationInfo = SA.getSpectraPairInfoForSVMClassification(lightMergedSpec, heavyMergedSpec, lightPrecMass, pairConfig['NMod'], pairConfig['CMod'], epsilon=2*epSTD)
            
            
            possPairs = [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[pair[0]] for heavyScanF in samePeptideClusters[pair[1]]]
#            possPairsList += [set(possPairs)]
            y += [1 if any([pair in progPairs[pairConfigName] for pair in possPairs]) else -1]
            #x += [{1: totalSharedPeaksRatio, 2: singleSymSharedPeaksRatio, 3: scanFDict[pair[0]]['precMass']}]
            x += [SVMClassificationInfo]

#            pairs[pairConfigName][pair] = getSharedPeaksRatio(lightSpec, heavySpec, pairConfig, epsilon)
#            print pair, pairs[pairConfigName][pair]


        x = svmutil.normalize_instances(x, svmRange)
        pLab = svmutil.svm_predict(y, x, svmModel)[0]
        pairs[pairConfigName] = {'test labels': pLab, 'true labels': y, 'pairs': testedDeltaPairs}
        times[pairConfigName] = time.time() - startTime

#    for i, pair in enumerate(testedPairs):
#        print pairs['same']['test labels'][i], pairs['same']['true labels'][i], pairs['same']['tested pairs'][i]
#        print processedInfo[progName][pairs['same']['tested pairs'][i][0]]['Peptide'], processedInfo[progName][pairs['same']['tested pairs'][i][1]]['Peptide']
        
    for pairConfigName in paramsDict['Pair Configurations']:
        truePairedScanFs = set()
        for pair in progPairs[pairConfigName]:
            truePairedScanFs.add(pair[0])
            truePairedScanFs.add(pair[1])

#        print 'number of true paired scanFs', len(truePairedScanFs)
#        print 'number of tested delta paired scanFs', len(testedDeltaPairedScanFs)
#        print 'intersection', len(testedDeltaPairedScanFs & truePairedScanFs)
#        print 'number of true pairs', len(progPairs[pairConfigName])
#        print 'number of tested pairs', len(testedPairs)
#        print 'intersection', len(set(progPairs[pairConfigName]) & testedPairs)
        
        outFile.write('\n%s. Time taken %f\n' % (pairConfigName, times[pairConfigName]))
#        writePairsDistributionInfo(outFile, progPairs[pairConfigName], pairs[pairConfigName]['pairs'], pairConfigName)
        outFile.write('\n' + '\t'.join(['Cutoff', 'TP', 'FP', 'TN', 'FN', 'Sensitivity', 'Precision']) + '\n')
        TP, FN, FP, TN = 0, 0, 0, 0
        pairsConsidered = set()
#        numPairsConsidered = 0
#        numPairsPassed = 0
        print len(pairs[pairConfigName]['test labels'])
        print len(pairs[pairConfigName]['pairs'])
        for i in range(len(pairs[pairConfigName]['test labels'])):
            clustPair = pairs[pairConfigName]['pairs'][i]
#            print clustPair
            possPairs = [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[clustPair[0]] for heavyScanF in samePeptideClusters[clustPair[1]]]
#            if len(set(possPairs) | pairs[pairConfigName]['possPairsList'][i]) != len(possPairs):
#                print 'ERROR: discrepancy in possible pairs for cluster pair', i, clustPair
#                print 'previous possible pairs', pairs[pairConfigName]['possPairsList'][i]
#                print 'current possible pairs', set(possPairs)
#            print  samePeptideClusters[clustPair[0]], samePeptideClusters[clustPair[1]], possPairs
            for pair in possPairs:
#                numPairsConsidered += 1
                if not (pair[0] in processedInfo[progName] and pair[1] in processedInfo[progName]):
                    continue
#                numPairsPassed += 1
                pairsConsidered.add(pair)
                if pairs[pairConfigName]['test labels'][i] == -1:
                    if pair not in progPairs[pairConfigName]:
                        TN += 1
                    else:
                        FN += 1
                else:
                    if pair in progPairs[pairConfigName]:
                        TP += 1
                    else:
                        FP += 1
                    
        cutoff = 'Cl_SVM'

        print 'length of progPairs', len(progPairs[pairConfigName])
#        print 'number of possible pairs', numPairsConsidered
#        print 'number of pairs tallied', numPairsPassed
        for pair in progPairs[pairConfigName]:
            if pair not in pairsConsidered:
                FN += 1
                
        try:
            sensitivity = int(float(TP)/(TP + FN) * 10000)/float(100)
        except ZeroDivisionError:
            sensitivity = 100
            
        try:
            precision = int(float(TP)/(TP + FP) * 10000)/float(100)
        except ZeroDivisionError:
            precision = 100       
        outFile.write('\t'.join([str(val) for val in [cutoff, TP, FP, TN, FN, sensitivity, precision]]) + '\n')
                    
    """
    for pairConfigName in paramsDict['Pair Configurations']:
        outFile.write('\n%s. Time taken %f\n' % (pairConfigName, times[pairConfigName]))
        print pairConfigName
        writePairsDistributionInfo(outFile, progPairs[pairConfigName], pairs[pairConfigName], pairConfigName)
        outFile.write('\n' + '\t'.join(['Cutoff', 'TP', 'FP', 'TN', 'FN', 'Sensitivity', 'Precision']) + '\n')
        for cutoff in np.arange(40, dtype=np.float64)/200:
            TP, FP, FN, TN = 0, 0, 0, 0
            for pair in pairs[pairConfigName]:
                if pair in progPairs[pairConfigName]:
                    if pairs[pairConfigName][pair] >= cutoff:
                        TP += 1
                    else:
                        FN += 1
                else:
                    if pairs[pairConfigName][pair] >= cutoff:
                        FP += 1
                    else:
                        TN += 1
            for pair in progPairs[pairConfigName]:
                if pair not in pairs[pairConfigName]:
                FN += 1
            try:
                sensitivity = int(float(TP)/(TP + FN) * 10000)/float(100)
            except ZeroDivisionError:
                sensitivity = 100

            try:
                precision = int(float(TP)/(TP + FP) * 10000)/float(100)
            except ZeroDivisionError:
                precision = 100
                
            outFile.write('\t'.join([str(val) for val in [cutoff, TP, FP, TN, FN, sensitivity, precision]]) + '\n')
    """
