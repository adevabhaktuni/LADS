'''
Created on Sep 29, 2011

@author: Arun
'''
import os
import sys
sys.path.insert(1, os.path.abspath(os.pardir))

import shlex, subprocess
import time
import numpy as np

import ArgLib
import DataFile
import Analytics
import CompareSearches
import Constants
#In its current incarnation, the false positives consist of true FPs, as well as cases in which the first member of each pair was found by one of the compared datasets (not LADS), and the other was not present in the same dataset
#the true positives undercount because each pair only lists its best possible match, excluding other pairs to that scan
#false negatives overcount for this aforementioned reason as well
def getArgs(options, argList):
    args = []
    for arg in argList:
        try:
            opt = getattr(options, arg)
            if opt != None:
                args.extend(['--' + arg, '\"' + str(opt) + '\"'])
        except AttributeError:
            pass
    
    return args

def executeProcess(interpreter, prog, args, outputBase, TDVoutput=True, logStdout=True):
    if TDVoutput:
        args.extend(['--output', getOutputName(outputBase, prog, '.tdv')])
    if logStdout:
        logFile = open(getOutputName(outputBase, prog, '.log'), 'w')

    cmd = ' '.join([interpreter] + [prog] + args)
    print cmd
    try:
        proc = subprocess.Popen(shlex.split(cmd), stdout=logFile, stderr=logFile)
    except NameError:
        proc = subprocess.Popen(shlex.split(cmd))
    
    print proc.communicate()
    if logStdout:
        logFile.close()

def getOutputName(outputBase, prog, ext):
    return outputBase + '_' + prog + '_UnitTest' + ext

def writeHistogram(outFile, vector, name, numBins=15):
#    hist, bins = np.histogram(vector, bins=numBins)
    hist, bins = np.histogram(vector, bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    outFile.write('\n%s Distribution. Max %s: %f \n' % (name, name, bins[-1]))
    outFile.write('\t'.join([name, 'Count']) + '\n')
    for i in range(hist.size):
        outFile.write('\t'.join([str(bins[i]), str(hist[i])]) + '\n')

def getNumFalseNegatives(pairs):
    numFN = 0
    for key in pairs.keys():
        if not pairs[key]:
            numFN += 1
    
    return numFN

def getVecScore(vector, numScans):
    return np.sum(vector) / numScans

def writeCategoryInfo(stats, outFile, categoryList, name=''):
    cols = ['Prog', 'Num Scans', 'Accuracy Score', 'Precision Score']
    outFile.write('\t'.join(cols) + '\n')
    for progName in stats.keys():
        dataDict = stats[progName]
        for category in categoryList:
            dataDict = dataDict[category]
            
        try:
            numScans = dataDict['numScans']
            outFile.write('\t'.join([progName, str(numScans), str(getVecScore(dataDict['accuracyVec'], numScans)), str(getVecScore(dataDict['precisionVec'], numScans))]) + '\n')
        except KeyError:
            outFile.write('\t'.join([progName, str(dataDict['numScans']), 'N/A', 'N/A']) + '\n')

    for progName in stats.keys():
        dataDict = stats[progName]
        for category in categoryList:
            dataDict = dataDict[category]
        try:
            outFile.write('\n%s %s Distributions\n' % (progName, name))
            writeHistogram(outFile, dataDict['precisionVec'], '%s %s Precision' % (progName, name), numBins=10)
            writeHistogram(outFile, dataDict['accuracyVec'], '%s %s Accuracy' % (progName, name), numBins=10)
        except KeyError:
            pass
        outFile.write('\n%s %s Ambiguous Edge Distribution\n' % (progName, name))
        outFile.write('\t'.join(['Num Ambig Edges', 'Count']) + '\n')
        try:
            for i in range(int(dataDict['Num Ambig Edges'].max())):
                outFile.write('\t'.join([str(i), str(len(np.where(dataDict['Num Ambig Edges'] == i)[0]))]) + '\n')
        except ValueError:
            outFile.write('\t'.join(['N/A', 'N/A']) + '\n')
        except KeyError:
            print 'ERROR: No data for Ambiguous Edges available'

        try:
            writeHistogram(outFile, dataDict['PScore'], '%s %s PScore' % (progName, name))
        except KeyError:
            print 'ERROR: No Data for PScore available'
    
unitTestName = 'LADS Unit Test'
if __name__ == '__main__':
    options = ArgLib.parse(['init', 'denovoscript', 'dtadir', 'config', 'model', 'output', 'columns', 'verbose', 'paircutoff', 'ppmsyserror', 'ppmstd', 'ppmpenalty', 'ambigpenalty', 'minedge', 'maxedge', 'alpha', 'lads', 'sequest', 'mascot', 'pepnovo', 'database', 'mainprogname', 'progdict', 'comp', 'subgraphcut', 'symbolmap', 'pnovo', 'peaks', 'combined', 'srchid'])
    
    outBase = os.path.splitext(options.output)[0]
    paramsDict = ArgLib.parseInitFile(options.init, options)

    interpreter = 'python2.6'
    if options.denovoscript:
        DNSprog = options.denovoscript
        progDict = {unitTestName: 'LADS'}
    else:
        progDict = {}

    if options.progdict:
        print options.progdict
        progDict = eval(options.progdict)
    else:
        ArgLib.getProgDict(Analytics.searchprogs, options, progDict=progDict)

    try:
        args = getArgs(options, ['init', 'dtadir', 'config', 'model', 'columns', 'paircutoff', 'ppmsyserror', 'ppmstd', 'ppmpenalty', 'ambigpenalty', 'minedge', 'maxedge', 'alpha', 'subgraphcut', 'symbolmap'])
        start = time.time()
        executeProcess(interpreter, DNSprog, args, outBase)
        DNSTime = time.time() - start
    except NameError:
        pass

    if not options.comp:
        args = getArgs(options, ['init', 'sequest', 'mascot', 'pepnovo', 'database', 'symbolmap', 'pnovo', 'peaks', 'combined', 'srchid'])
        if options.lads:
            ladsDict = eval(options.lads)
        else:
            ladsDict = {}

        if options.denovoscript:
            ladsDict.update({getOutputName(outBase, DNSprog, '.tdv'): unitTestName})
            args.extend(['--lads', '\"' + str(ladsDict) + '\"'])
        else:
            args.extend(['--lads', '\"' + str(ladsDict) + '\"'])
            
        executeProcess(interpreter, 'CompareSearches.py', args, outBase)

    infoMap = DataFile.getDBInfo(options.database, key='infoMap')
    if progDict[options.mainprogname] == 'LADS':
        getPairStats = True
        mainProgFields = ['PScore', 'Num Ambig Edges']
    else:
        getPairStats = False
        mainProgFields = [infoMap[progDict[options.mainprogname]]['Score']]

    if options.denovoscript:
        stats = Analytics.getCompStats(getOutputName(outBase, 'CompareSearches.py', '.tdv'), unitTestName, progDict, infoMap, paramsDict)
    elif options.mainprogname:
        unitTestName = options.mainprogname
        if not options.comp:
            stats = Analytics.getCompStats(getOutputName(outBase, 'CompareSearches.py', '.tdv'), options.mainprogname, progDict, infoMap, paramsDict, getPairStats=getPairStats, mainProgFields=mainProgFields)
        else:
            stats = Analytics.getCompStats(options.comp, options.mainprogname, progDict, infoMap, paramsDict, getPairStats=getPairStats, mainProgFields=mainProgFields)
        
    outFile = open(options.output, 'w')

    outFile.write('\nOverall Comparison Statistics\n')
    writeCategoryInfo(stats, outFile, ['composite'], name='Composite')
    outFile.write('\nTrue Pairs\n')
    if getPairStats:
        for pairConfigName in paramsDict['Pair Configurations']:
            writeCategoryInfo(stats, outFile, ['truepairs', pairConfigName], name='%s True Pairs' % (pairConfigName,))
        outFile.write('\nFalse Pairs\n')
        for pairConfig in paramsDict['Pair Configurations']:
            writeCategoryInfo(stats, outFile, ['falsepairs', pairConfigName], name='%s False Pairs' % (pairConfigName,))
        outFile.write('\nUnpaired\n')
        writeCategoryInfo(stats, outFile, ['unpaired'], name='Unpaired')

    outFile.write('\nOther scans (no comparison available)\n')
    writeCategoryInfo(stats, outFile, ['other'], name='Other')
    

    if getPairStats:
        outFile.write('\n' + '\t'.join(['Program'] + [progName for progName in stats.keys()]) + '\n')
        for pairConfigName in paramsDict['Pair Configurations']:
            outFile.write('\t'.join(['TP %s (underestimated)' % (pairConfigName,)] + [str(stats[progName]['truepairs'][pairConfigName]['numScans']) for progName in stats.keys()]) + '\n')
            outFile.write('\t'.join(['FP %s (overestimated)' % (pairConfigName,)] + [str(stats[progName]['falsepairs'][pairConfigName]['numScans']) for progName in stats.keys()]) + '\n')
            outFile.write('\t'.join(['FN %s (overestimated)' % (pairConfigName,)] + [str(getNumFalseNegatives(stats[progName]['pairsDict'][pairConfigName])) for progName in stats.keys()]) + '\n')
    
    try:
        outFile.write('\n\nSequencing time: %f' % (DNSTime,) + '\n')
    except NameError:
        pass
    
    outFile.close()
        
        
        
    
    
    
    
    
    
    
