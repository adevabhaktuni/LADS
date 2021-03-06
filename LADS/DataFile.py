'''
Created on Jun 1, 2011

@author: Arun
'''

import numpy as np
import os
import csv
import pickle
import anydbm as dbm
import ConfigParser
import glob

import Constants

#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.2499.2499.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3367.3367.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3368.3368.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3383.3383.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3384.3384.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3876.3876.1.dta"
#"C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3877.3877.1.dta"


'''
This method will return a tuple containing the precursor mass and charge read
from the first line of a .dta file
@param absPath: absolute path to .dta file
@return: tuple, first element is precursor mass, last element is charge
@author: Arun
'''
def getPrecMassAndCharge(absPath):
    fin = open(absPath, 'r')
    
    firstline = fin.readline()
    length = len(firstline)
    for i in range(length):
        if firstline[i] == '\t':
            precMass = float(firstline[0:i])
            charge = float(firstline[i + 1:length])
    
    fin.close()  
    return precMass, charge

def getDTAFNamesInDir(dirPath):
    dtas = []
    dirList = os.listdir(dirPath)
    for file in dirList:
        ext = os.path.splitext(file)[1]
        if ext == '.dta':
            dtas.extend([dirPath + file])
    
    return dtas

def getDTAByScanNum(dtaDir, scanNum):
    return glob.glob(dtaDir + '/*.%(scanNum)04i.*dta' % {'scanNum': scanNum})

def indexDataByKey(data, key='ScanF', overrideKey=None, dtyper=int):
    dataDict = {}
    for datum in data:
        if dtyper(datum[key]) in dataDict and overrideKey != None:
            if float(datum[overrideKey]) > float(dataDict[dtyper(datum[key])][overrideKey]):
                dataDict[dtyper(datum[key])] = datum
        else:
            dataDict[dtyper(datum[key])] = datum
            
    return dataDict  


def combineDatafiles(data1, data2, key='ScanF', overrideKey1=None, overrideKey2=None, dtyper=int):
    indexedData1 = indexDataByKey(data1, key=key, overrideKey=overrideKey1, dtyper=dtyper)
    indexedData2 = indexDataByKey(data2, key=key, overrideKey=overrideKey2, dtyper=dtyper)

    cols1 = indexedData1.values()[0].keys()
    cols2 = indexedData2.values()[0].keys()
    try:
        cols1.remove(key)
        cols2.remove(key)
    except ValueError:
        pass

    combinedDataArr = []
    allKeyVals = set(indexedData1.keys()) | set(indexedData2.keys())
    for keyVal in allKeyVals:
        scanData = {}
        scanData[key] = keyVal

        if keyVal in indexedData1:
            for col in cols1:
                scanData[col] = indexedData1[keyVal][col]
        else:
            for col in cols1:
                scanData[col] = None

        if keyVal in indexedData2:
            for col in cols2:
                scanData[col] = indexedData2[keyVal][col]
        else:
            for col in cols2:
                scanData[col] = None

        combinedDataArr += [scanData]

    return combinedDataArr

def getScanInfo(fname, fields=None, delimiter=',', skipLines=0):
    runInfo = []
    rec = csv.reader(open(fname, 'r'), delimiter=delimiter)
    for i in range(skipLines):
        rec.next()
    
    cols = rec.next()
    
    print 'cols', cols
    if not fields:
        fields = [col for col in cols]
                
    fieldInds = []
    for i in range(len(cols)):
        if cols[i] in fields:
            fieldInds += [i]
    
    for row in rec:
        rowInfo = {}
        for ind in fieldInds:
            rowInfo[cols[ind]] = row[ind]
            
        runInfo += [rowInfo]
    
    return runInfo
    
    
def getScanNum(dtaFName):
    fname = os.path.basename(dtaFName)
    base = os.path.splitext(fname)[0]
    base = os.path.splitext(base)[0]
    base = os.path.splitext(base)[0]
    scanNum = os.path.splitext(base)[1]
    return int(scanNum[1:])
    

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getMassIntPairs(absPath):
    fin = open(absPath, 'r')

    massIntPairs = np.zeros((file_len(absPath) - 1 , 2))
    fin.readline()
    
    for i, line in enumerate(fin): 
        length = len(line)
        for j in range(length):
            if line[j] == '\t':
                firstnum = float(line[0:j])
                secnum = float(line[j + 1:length])
                massIntPairs[i, 0] = firstnum
                massIntPairs[i, 1] = secnum
                break

    fin.close()
    return massIntPairs

def getScoringMatrix(absPath):
    fin = open(absPath, 'r')
    line = fin.readline()
    while line[0] == '#':
        line = fin.readline()
    
    i = 0
    AAMap = {}
    for char in line:
        if char != ' ' and char != '\n':
            AAMap[char] = i
            i += 1
    
    dim = max(AAMap.values()) + 1
    matrix = np.zeros((dim, dim), dtype=np.dtype('int'))
    for i in range(dim):
        line = fin.readline()
        char = line[0]
        if AAMap[char] != i:
            print 'ERROR: scoring matrix not properly formatted for AA %s' % char
        else:
            data = line[2:]
            k = j = 0
            while j < dim:
                if data[k] == '-':
                    matrix[i][j] = int(data[k:k + 2])
                    k += 1
                    j += 1
                elif data[k] != ' ':
                    matrix[i][j] = int(data[k])
                    j += 1
                
                k += 1
    
    return AAMap, matrix

def getDBInfo(db, key=None):
    db = dbm.open(db, 'r')
    if key:
        return pickle.loads(db[key])
    else:
        dbDict = {}
        for key in db.keys():
            dbDict[key] = pickle.loads(db[key])
        return dbDict

#    modDict = {'SEQUEST': {'silac(k)': '#', 'oxidation': '*'},
#                'MASCOT': {'silac(k)': '*', 'oxidation': '#'},
#               'LADS': {'silac(k)': '*', 'oxidation': '#'}}
#    AADict = {'SEQUEST': {'X': 'I'}, 'MASCOT': {}, 'LADS': {'-': 'X'}}
def generateSeqMap(progDict, symbolMap, paramsDict):
    seqMap = {}
    for prog in progDict:
        seqMap[prog] = {'AAs':{}, 'Mods': {}}
        if progDict[prog] != 'LADS':
            seqMap[prog]['AAs']['L'] = 'I'
        elif progDict[prog] == 'LADS':
            seqMap[prog]['AAs']['-'] = 'X'
        for aa in Constants.origAAs.keys():
            seqMap[prog]['AAs'][aa] = aa
        try:
            seqMap[prog]['AAs'].update(symbolMap[prog]['AAs'])
        except KeyError:
            print 'AAs dictionary for program %s not present in symbol map' % (prog,)
            
    for mod in paramsDict['Diff Mods'].keys():
        modType = paramsDict['Diff Mods'][mod][0]
        for prog in progDict:
            try:
                progMod = symbolMap[prog]['Mods'][modType]
                if symbolMap[prog]['Mods'][modType] != 'static':
                    seqMap[prog]['Mods'][progMod] = mod
                else:
                    # modify all amino acids with static mod
                    seqMap[prog]['AAs'][paramsDict['Diff Mods'][mod][1]] = seqMap[prog]['AAs'][paramsDict['Diff Mods'][mod][1]] + mod
            except KeyError:
                print 'Modification of type %s unaccounted for in modDict for program %s' % (modType, prog)
            
    for modData in paramsDict['Static Mods'].keys():
        modType = modData[0]
        for prog in progDict:
            try:
                progMod = symbolMap[prog]['Mods'][modType]
                seqMap[prog]['Mods'][progMod] = ''
            except KeyError:
                print 'Modification of type %s unaccounted for in modDict for program %s' % (modType, prog)
        
    return seqMap

paramHandler = {

'Static Mods': lambda datum, args: Constants.addStaticMod(datum, args),
'Diff Mods': lambda datum, args: Constants.addDiffMod(datum, args),
'Models': lambda datum, args: parseModelInfo(datum, args),
'Enzyme': lambda datum, args: parseEnzymeInfo(datum, args),
'Amino Acids': lambda datum, args: Constants.addAA(datum, args),
'LADS Parameters': lambda datum, args: parseLADSParametersInfo(datum, args),
'Pair Configurations': lambda datum, args: parsePairConfiguration(datum, args),
'Cluster Configuration': lambda datum, args: parseClusterConfiguration(datum, args),
'Columns': lambda datum, args: parseColumns(datum, args),
'Alpha': lambda datum, args: parseColumns(datum, args),
'Minedge': lambda datum, args: parseColumns(datum, args),
'Maxedge': lambda datum, args: parseColumns(datum, args),
'Ambiguity Penalty': lambda datum, args: parseColumns(datum, args),
'PPM Penalty': lambda datum, args: parseColumns(datum, args),
'ppmstd': lambda datum, args: parseColumns(datum, args),
'Symbol Map': lambda datum, args: parseColumns(datum, args),
'Config File': lambda datum, args: parseColumns(datum, args),
'Model': lambda datum, args: parseColumns(datum, args)
}

def parseParams(fname):
    paramsDict = {}
    params = ConfigParser.ConfigParser()
    params.read(fname)
    for section in params.sections():
        #print "params section " + section
        paramsDict[section] = {}
        for datum in params.options(section):
            try:
                paramsDict[section].update(paramHandler[section](datum, params.get(section, datum).split(' ')))
            except KeyError:
                raise KeyError('%s not a valid category. Valid categories are: %s' % (section, str(paramHandler.keys())))
    return paramsDict

def parseColumns(datum, args):
    return {datum: args[0]}

def parseModelInfo(datum, args):
    return {datum: {'config': args[0], 'model': args[1]}}

def parseEnzymeInfo(datum, args):
    if datum == 'specificity':
        return {datum: parseEnzymeSpecificity(args)}
    else:
        return {datum: args[0]}

def parseLADSParametersInfo(datum, args):
    return {datum: args[0]}
    
def parseEnzymeSpecificity(motifs):
    flankingPairs = []
    for motif in motifs:
        cleavageInd = motif.index(';')
        before = motif[:cleavageInd]
        after = motif[cleavageInd + 1:]
        if before != '.*':
            before = before.split('|')
        else:
            before = (None,)
        
        if after != '.*':
            after = after.split('|')
        else:
            after = (None,)
        
        flankingPairs.extend([(bef, aft) for bef in before for aft in after])
    
    return flankingPairs


def parsePairConfiguration(datum, pairConfig):
    pairConfiguration = {'NStatic': np.array([float(pairConfig[0])]), 'CStatic': np.array([float(pairConfig[1])]), 'NMod': float(pairConfig[2]), 'CMod': float(pairConfig[3]), 'NModSymbol': pairConfig[4], 'CModSymbol': pairConfig[5], 'Model': pairConfig[6]}
    if pairConfiguration['NModSymbol'] in ['None', '0', 'False']:
        pairConfiguration['NModSymbol'] = ''
    if pairConfiguration['CModSymbol'] in ['None', '0', 'False']:
        pairConfiguration['CModSymbol'] = ''
        
    return {datum: pairConfiguration}

def parseClusterConfiguration(datum, args):
    return {datum: args[0]}

"""
def plotSpectra(massIntPairs, precMass = 'None', charge = 'Undefined'):
    masses = massIntPairs[:,0]
    intensities = massIntPairs[:,1] 
    minMass = masses.min()
    maxMass = masses.max()
    maxInt = intensities.max()

    for i in range(len(masses)):
        pl.axvline(x = masses[i], ymin = 0, ymax = intensities[i]/(maxInt * 1.1))

    pl.ylim([0, maxInt * 1.1])
    pl.xlabel( 'M/Z, NumSpectra: ' + str(np.size(masses)) )
    pl.ylabel('Intensity')
    pl.title('Precursor Mass: ' + str(precMass) + ', Charge: ' + str(charge))
    pl.xlim([minMass * 0.8,maxMass * 1.2])
    pl.show()
"""
if __name__ == '__main__':
    #print getScanInfo('adevabhaktuni_1310166306.csv')
    """
    dirPath = 'C:\\Users\\Arun\\Pythonprojects\\DeNovoSequencing\\src\\SpectraCorrelation\\LF2_short_HCD+CID_ath001862_244\\'
    names = getDTAFNamesInDir(dirPath)
    fname = 'C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3383.3383.1.dta'
    massIntPairs = getMassIntPairs(fname)
    precMass, charge = getPrecMassAndCharge(fname)
    plotSpectra(massIntPairs, precMass, charge)
    """
    #fields = ['LIGHT SCAN', 'HEAVY SCAN', 'SEQ', 'SCORE', 'AMBIGUOUS EDGES', 'M+H', 'EPSILON', 'SHARED PEAKS RATIO', 'SEQUENCING TIME']
    #print getScanInfo('C:\\Users\\Arun\\Proteomics Files\\ath001862UPen10KPen15LRRestrictTest.tdv', fields, delimiter='\t')
    paramsDict = parseParams('./Misc/LADS_SILAC_Trypsin.ini')
    print paramsDict
    seqMap = generateSeqMap(['SEQUEST', 'MASCOT'], paramsDict)
    print seqMap
    
