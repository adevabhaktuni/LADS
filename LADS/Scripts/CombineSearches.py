'''
Created on Nov 8, 2011

@author: Arun
'''

import ArgLib
import Analytics as An
import DataFile

import numpy as np
import pickle

def getAllScanF(processedInfo):
    scanFs = np.array([], dtype=np.dtype('int16'))
    for progName in processedInfo.keys():
        scanFs = np.append(scanFs, np.array(processedInfo[progName].keys(), dtype=np.dtype('int16')))

    return np.unique(scanFs)

def getPeptideData(progName, progDict, scanF):
    try:
        seq = processedInfo[progName][scanF]['Peptide']
        if progDict[progName] == 'LADS':
            return (seq, processedInfo[progName][scanF]['Ambiguous Edges'])
        else:
            return (seq, None)
    except KeyError:
        return False

def parseDBScans(fDict, prog, seqMap, dbDict):
    processedInfo = {}
    for csvfile in fDict.keys():
        MASCOTData = DataFile.getScanInfo(csvfile, dbDict[prog]['fields'], delimiter=',')
        processedInfo[fDict[csvfile]] = An.preprocessDatabaseScanInfo(MASCOTData, seqMap[fDict[csvfile]], dbDict[prog]['fieldmap'])
    
    return processedInfo

#Number argument refers to minimum number of search prog results which have the same peptide for it to be included in the final output
if __name__== '__main__':
    options = ArgLib.parse(['init', 'sequest', 'lads', 'mascot', 'output', 'database', 'symbolmap', 'number'])
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    dbDict = DataFile.getDBInfo(options.database)
    progDict = ArgLib.getProgDict(An.searchprogs, options)
    
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap(progDict, symbolMap, paramsDict)
    
    if hasattr(options, 'number'):
        minNumScans = int(options.number)
    else:
        minNumScans = 1
        
    processedInfo = {}  
    if options.lads:
        LADSdict = eval(options.lads)
        for tdvfile in LADSdict.keys():
            LADSScanInfo = DataFile.getScanInfo(tdvfile, dbDict['LADS']['fields'], delimiter='\t')
            processedInfo[LADSdict[tdvfile]] = An.preprocessLADSScanInfo(LADSScanInfo, seqMap[LADSdict[tdvfile]], paramsDict['LADS Parameters']['pair configurations'], dbDict['LADS']['fieldmap'])
        
    if options.mascot:
        MASCOTdict = eval(options.mascot)
        processedInfo.update(parseDBScans(MASCOTdict, 'MASCOT', seqMap, dbDict))
        
    if options.sequest:
        SEQUESTdict = eval(options.sequest)
        processedInfo.update(parseDBScans(SEQUESTdict, 'SEQUEST', seqMap, dbDict))
        
    cols = ['ScanF']
    progNames = processedInfo.keys()
    cols.extend([val for val in dbDict[progDict[progNames[0]]]['cols']])
    
    outFile = open(options.output, 'w')
    outFile.write(','.join([col for col in cols]) + '\n')
    
    for i in getAllScanF(processedInfo):
        scanData = {}
        scanData['ScanF'] = i
        peptides = {}
        for progName in progNames:
            peptide = getPeptideData(progName, progDict, i)
            if peptide:
                try:
                    peptides[peptide] += [progName]
                except KeyError:
                    peptides[peptide] = [progName]
        
        maxScore = None
        chosenProg = None
        for peptide in peptides:
            if len(peptides[peptide]) >= minNumScans:
                for progName in peptides[peptide]:
                    if not (maxScore == None) or processedInfo[progName][i][dbDict['infoMap'][progDict[progName]]['Score']] > maxScore:
                        maxScore = processedInfo[progName][i][dbDict['infoMap'][progDict[progName]]['Score']]
                        chosenProg = progName
        
        if chosenProg:
            for col in dbDict[progDict[chosenProg]]['cols']:
                try:
                    scanData[col] = processedInfo[progName][i][col]
                except KeyError:
                    scanData[col] = None
            scanData['Peptide'] = 'X.' + scanData['Peptide'] + '.X'
            outFile.write(','.join([str(scanData[col]) for col in cols]) + '\n')
