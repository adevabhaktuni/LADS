from collections import defaultdict
import itertools

import numpy as np
import Analytics as An
import copy


def getSensitivityAndPrecisionArrs(scoreArr, critArr, critCutOff=0.9, numPoints=100):
    truePosInds = np.where(critArr >= critCutOff)
    trueScores = scoreArr[truePosInds]
    falseScores = np.delete(scoreArr, truePosInds)

    print scoreArr
    maxScore = np.amax(scoreArr)
    minScore = np.amin(scoreArr)
    scoreCutoffs = np.linspace(minScore, maxScore, num=numPoints)

    sens, prec = [], []
    for cutoff in scoreCutoffs:
        TP = len(np.where(trueScores >= cutoff)[0])
        FP = len(np.where(falseScores >= cutoff)[0])
        sens += [float(TP)/len(trueScores)]
        try:
            prec += [float(TP)/(TP + FP)]
        except ZeroDivisionError:
            prec += [1]

    return sens, prec, len(trueScores), len(falseScores)

def writeSensitivityAndPrecisionToFile(scoreArr, sensArr, precArr, outFileName, allScoresArr = None):
    outFile = open(outFileName, 'w')

    cols = ['Score Cutoff', 'Sensitivity', 'Precision', 'Num Results']
    if allScoresArr != None:
        cols += ['Num Total Results']
    outFile.write('\t'.join([col for col in cols]) + '\n')
    cutOffs = np.linspace(np.amin(scoreArr), np.amax(scoreArr), num=len(sensArr))
    
    for i in range(cutOffs.size):
        outFile.write('\t'.join([str(cutOffs[i]), str(sensArr[i]), str(precArr[i]), str( np.where(scoreArr >= cutOffs[i])[0].size ), str( np.where(allScoresArr >= cutOffs[i])[0].size ) if allScoresArr != None else 'None']) + '\n')
    
    outFile.close()

def getTrueAndFalseScores(scoreArr, critArr, critCutOff=0.9):
    truePosInds = np.where(critArr >= critCutOff)
    trueScores = scoreArr[truePosInds]
    falseScores = np.delete(scoreArr, truePosInds)

    return trueScores, falseScores

def calculateFDRArray(scoreArr, critArr, critCutOff=0.9, numPoints=100):
    trueScores, falseScores = getTrueAndFalseScores(scoreArr, critArr, critCutOff)

    FDRs = []
    scoreCutoffs = np.linspace(np.amin(scoreArr), np.amax(scoreArr), num=numPoints)
    for cutoff in scoreCutoffs:
        TP = len(np.where(trueScores >= cutoff)[0])
        FP = len(np.where(falseScores >= cutoff)[0])
        try:
            FDR = float(FP)/(TP + FP)
        except ZeroDivisionError:
            if FP > 0:
                FDR = 1
            else:
                FDR = 0

        FDRs += [FDR]
    
    return scoreCutoffs, FDRs

def writeFDRArr(scoreCutoffs, FDRArr, outFileName):
    outFile = open(outFileName, 'w')
    
    cols = ['Score Cutoff', 'FDR']
    outFile.write('\t'.join([col for col in cols]) + '\n')

    for i in range(scoreCutoffs.size):
        outFile.write('\t'.join([str(scoreCutoffs[i]), str(FDRArr[i])]) + '\n')

    outFile.close()


def getUniquePeptDict(scanDict, scoreKey, peptideKey, scanKey = 'ScanF', nullVal = 'None', noStrip=['#'], datasets=None):
    scanFDict = defaultdict(lambda: dict([(dataset, []) for dataset in datasets]))
    uniquePeptDict = {}
    if datasets == None:
        datasets = scanDict.keys()

    for dataset in datasets:
        for item in scanDict[dataset]:
            if item[peptideKey] == nullVal:
                continue
            
            strippedPept = An.stripModifications(item[peptideKey], noRemove=noStrip)
            
            if strippedPept in uniquePeptDict and float(item[scoreKey]) > float(uniquePeptDict[strippedPept][scoreKey]):
                uniquePeptDict[strippedPept] = item
            elif strippedPept not in uniquePeptDict:
                uniquePeptDict[strippedPept] = item


            scanFDict[strippedPept][dataset] += [item[scanKey]]

    return uniquePeptDict, scanFDict


def indexByFoundIn(scanFDict, uniquePeptDict):
    numDatasets = len(scanFDict.values()[0])

    indexedByFoundIn = dict([(i, []) for i in range(1, numDatasets + 1)])

    for peptide in scanFDict:
        numFound = 0

        for dataset in scanFDict[peptide]:
            if len(scanFDict[peptide][dataset]) > 0:
                numFound += 1

        indexedByFoundIn[numFound] += [uniquePeptDict[peptide]]

    return indexedByFoundIn


# ROCsByFoundIn is list of ROCs, each of whcih corresponds to a different 'found in' bin
def optimizationOverFoundIn(ROCsByFoundIn, precBin=200):

    # Data structure initialization
    numCategories = len(ROCsByFoundIn)

    optimalIndsAtDesiredPrecs = defaultdict(lambda: {'Inds': [], 'Precision': 0, 'Sensitivity': 0, 'Num Results': 0, 'Num Total Results': 0})

    totalNumCorrect = 0
    for ROCs in ROCsByFoundIn:
        if len(ROCs) > 0:
            totalNumCorrect += float(ROCs[0]['Precision']) * int(ROCs[0]['Num Results']) / float(ROCs[0]['Sensitivity'])

    indList = [range(len(item)) for item in ROCsByFoundIn]

    ROCslength = [len(item) for item in ROCsByFoundIn]
    counts = np.zeros((numCategories, max(ROCslength), 2))
    
    for i, ROCs in enumerate(ROCsByFoundIn):
        for j, item in enumerate(ROCs):
            counts[i][j] = [float(item['Precision']) * float(item['Num Results']), float(item['Num Results'])]

    num = 0
    for inds in itertools.product(*indList):
        prec, sens, numResults = calculateStats(counts, inds, totalNumCorrect)

        binnedPrec = round(prec * precBin)
        if optimalIndsAtDesiredPrecs[binnedPrec]['Sensitivity'] < sens:
            optimalIndsAtDesiredPrecs[binnedPrec] = {'Inds': inds, 'Precision': prec, 'Sensitivity': sens, 'Num Results': numResults}

        if num % 1000000 == 0:
            print num
            print optimalIndsAtDesiredPrecs
        num += 1
        
    return optimalIndsAtDesiredPrecs

def calculateStats(counts, inds, totalCorrect):
    numResults = 0
    numCorrect = 0
    for i, ind in enumerate(inds):
        numResults += counts[i][ind][1]
        numCorrect += counts[i][ind][0]

    return numCorrect/numResults, numCorrect/totalCorrect, numResults
    

def getUniqueClusters(postScores, epsilon=0.1, delta=6.03766):
    currentItem = postScores[0]
    uniquesByCluster = []
    scanFs = []
    
    for item in postScores:
        obsMH = float(item['M+H'])
        if abs(obsMH - float(currentItem['M+H'])) < epsilon or abs(abs(obsMH - float(currentItem['M+H'])) - delta) < epsilon:
            scanFs += [item['ScanF']]
            if float(item['LADS Post Score']) > float(currentItem['LADS Post Score']):
                currentItem = item
        else:
            currentItem['ScanFs'] = scanFs
            uniquesByCluster += [currentItem]
            
            currentItem = item
            scanFs = [item['ScanF']]
            
    currentItem['ScanFs'] = scanFs
    uniquesByCluster += [currentItem]
                
    return uniquesByCluster
            
        
            
