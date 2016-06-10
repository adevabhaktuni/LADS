'''
Created on Sep 7, 2011

@author: Arun
'''

import sys
import os

sys.path.insert(1, os.path.abspath(os.pardir))
libpath = os.path.abspath(os.pardir) + "/libsvm-3.12/python"
sys.path.append(libpath)

import glob
import time
import pickle
import numpy as np
import copy
import svmutil

import ArgLib
import Constants
import DeNovoSequencer as DNS
import SpectraAlignment as SA
import Analytics
import DataFile
import ProbNetwork as PN
import Queue
import multiprocessing
from multiprocessing import Process, current_process, Pool, Pipe
import threading
from threading import Lock, BoundedSemaphore
from time import sleep
import commands
import fcntl

print_lock = Lock()
spectrum_lock = Lock()
out_lock = Lock()

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = DataFile.getScanNum(dta)
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'sequenced': False}
    
    return scanFDict

def getSharedPeaksRatio(lightPath, heavyPath, epsilon):
    lightPairs = DataFile.getMassIntPairs(lightPath)
    heavyPairs = DataFile.getMassIntPairs(heavyPath)
    N, C = SA.getNandCIons(lightPairs, heavyPairs, pairConfig['NMod'], pairConfig['CMod'], epsilon=epsilon)
    return SA.getSharedPeaksRatio(lightPairs, heavyPairs, N, C)

def validateHeavySequence(seq, heavySeqMap, ambigEdges):
    try:
        if seq != '-':
            heavySeq = Analytics.preprocessSequence(seq, heavySeqMap, replaceExistingTerminalMods=True, ambigEdges=ambigEdges)
            AAs = Analytics.getAllAAs(heavySeq, ambigEdges=ambigEdges)
            return True
        else:
            return False
    except KeyError:
        return False

def getPairs(pairs, xVals):
    for pair in pairs:
        lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
        heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
        lightPrecMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))

        epSTD = (float(paramsDict['ppmstd']['value'])) * 10 ** -6 * lightPrecMass

        lightMergedSpec = SA.mergeSpectra(lightSpecs, epsilon=2*epSTD)
        heavyMergedSpec = SA.mergeSpectra(heavySpecs, epsilon=2*epSTD)

        svmClassificationData = SA.getSpectraPairInfoForSVMClassification(lightMergedSpec, heavyMergedSpec, lightPrecMass, NMod=pairConfig['NMod'], CMod=pairConfig['CMod'], epsilon=2*epSTD)
        xVals.put([svmClassificationData])
    
    return xVals

def getPairsThread(pairs):
    xVals = []
    res = Queue()
    processors = multiprocessing.cpu_count()
    width = len(pairs)/(processors -1)
    processes = []
    print "Initial pairs len " + str(len(pairs)) + " " + str(width)
    ##processes = [Process(target=getPairs, args=(pairs, xVals)) for i in xrange(processors)]

    available = processors -1
    for i in xrange(available):
        start = width * i
        end = ((i + 1) * width)
        p = Process(target=getPairs, args=(pairs[start:end], res))
        processes.append(p)

    if (width * processors < len(pairs)):

        p = Process(target=getPairs, args=(pairs[(width*processors):(len(pairs))], res))
        processes.append(p)
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    xVals = res.get()
    for r in range(res.qsize()):
        xVals += res.get()

    return xVals

def getSequencing(pair, sharedPeaks, paramsDict, outFile, res):
    global print_lock, spectrum_lock

    result = []

    scanData = {}
    lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
    heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
    precMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))
    
    epMean = options.ppmsyserror * precMass * 10**-6
    epSTD = options.ppmstd * precMass * 10**-6
                
    scanData['shared peaks ratio'] = sharedPeaks

    s1 = time.time()
    sharedInfo, starts, ends, deltas, G = DNS.prepPairedSpectrumGraph(lightSpecs, heavySpecs, precMass, addEnds, ppmSTD=options.ppmstd, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], verbose=options.verbose)
    scanData['M+H'] = precMass
    
    specs = []
    for massIntPairs in lightSpecs:
        specs += [PN.Spectrum(PNet, precMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=massIntPairs)]
    for massIntPairs in heavySpecs:
        specs += [PN.Spectrum(PNet, precMass + pairConfig['NMod'] + pairConfig['CMod'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD, spectrum=massIntPairs)]
    for spec in specs:
        spec.initializeNoiseModel()

    # with spectrum_lock:
    temp = DNS.getSpectrumGraphDataThread(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, termModHash=termModHash, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose)
    temp_scan = temp[0]
    peps = temp[1]
    scanData.update(temp_scan)
    
    scanData['pair configuration'] = pairConfigName

    with print_lock:
        print 'Now sequencing light scan(s) %s, heavy scan(s) %s with shared peaks ratio %f \n' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), scanData['shared peaks ratio'])
        # out.append('Now sequencing light scan(s) ' + str(samePeptideClusters[pair[0]]) + ', heavy scan(s) ' + str(samePeptideClusters[pair[1]]) + ' with shared peaks ratio ' + str(scanData['shared peaks ratio']) + ' \n' )
        Ord = np.argsort(-1 * np.array(scanData['over_scores']))
        if scanData['blind'] == 0:
            for i in range(min(Ord.size, 10)):
                try:
                    print 'Score: ', peps[0][Ord[i]], 'Seq: ', ''.join(peps[1][Ord[i]])
                    # out.append('Score: ' + str(peps[0][Ord[i]]) + ' Seq: ' + ''.join(peps[1][Ord[i]]))
                except TypeError:
                    print 'Score: ', peps[0][Ord[i]], 'Seq: ', peps[1][Ord[i]]
                    # out.append('Score: ' + str(peps[0][Ord[i]]) + ' Seq: ' + str(peps[1][Ord[i]]))
        elif scanData['blind'] == 1:
            for i in range(min(Ord.size, maxNum)):
                try:
                    print 'Score: ', peps[0][Ord[i]], 'Seq: ', ''.join(peps[1][Ord[i]][0]), 'Mod Names: ', peps[2][Ord[i]][1]
                    # out.append('Score: ' + str(peps[0][Ord[i]]) + ' Seq: ' + ''.join(peps[1][Ord[i]][0]) + ' Mod Names: ' + peps[2][Ord[i]][1])
                except TypeError:
                    print 'Score: ', peps[0][Ord[i]], 'Seq: ', peps[1][Ord[i]][0], 'Mod Names: ', peps[2][1]
                    # out.append('Score: ' + str(peps[0][Ord[i]]) + ' Seq: ' + peps[1][Ord[i]][0] +  ' Mod Names: ' + peps[2][1])
        
        scanData['sequencing time'] = time.time() - s1
        print '\nTime Taken:', time.time() - s1, '\n'    
    # out.append('\nTime Taken: ' + str(time.time() - s1) + '\n')

    if validateHeavySequence(scanData['seq'], heavySeqMap, scanData['ambiguous edges']):
        for scanF in samePeptideClusters[pair[0]] + samePeptideClusters[pair[1]]:
            scanFDict[scanF]['sequenced'] = True
        if options.output:
            for pair in [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[pair[0]] for heavyScanF in samePeptideClusters[pair[1]]]:
                scanData['light scan'] = int(pair[0])
                scanData['heavy scan'] = int(pair[1])                  
                # outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
                # print str(scanData[col])
                res.append([str(scanData[col]) for col in cols])
        else:
            print 'WARNING: Invalid sequence! Unsuccessful sequencing of %s and %s with pair configuration %s' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), pairConfigName)

    exit(0)

def getSequencingThread(pairs, xVals, paramsDict, outFile, cols, pLab):    
    processes = []
    results = []
    res = multiprocessing.Queue()
    man = multiprocessing.Manager()
    L = man.list()

    for i, pair in enumerate(pairs) :
        if pLab[i] == -1:
            continue
            
        sharedPeaks = xVals[i][1]
        p = Process(target=getSequencing, args=((pair), sharedPeaks, paramsDict, outFile, L))
        processes.append(p)
    
    curr = 0
    sent = 0
    for p in processes:
        count = commands.getoutput('ps -a | grep python | wc -l')
        while (int(count) > int(multiprocessing.cpu_count() + 4) and curr < (len(processes)-4)):
            count = commands.getoutput('ps -a | grep python | wc -l')
            prev = processes[curr]
            if(prev.exitcode == None and curr < (len(processes) - multiprocessing.cpu_count())):
                sleep(0.1)
            else:
                prev.terminate() 
                curr += 1
        
        p.start()

    for p in processes:
        p.join()

    for l in L:
        for j in l:
            outFile.write(str(j) + '\t')
        outFile.write('\n')

if __name__ == '__main__' :
    options = ArgLib.parse(['init', 'dtadir', 'config', 'model', 'output', 'columns', 'verbose', 'paircutoff', 'ppmsyserror', 'ppmstd', 'ppmpenalty', 'ambigpenalty', 'minedge', 'maxedge', 'alpha', 'subgraphcut', 'symbolmap'])
    epStep = 0.00025
    maxEp = 0.1
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap({'LADS Unit Test': 'LADS'}, symbolMap, paramsDict)
    
    if options.columns:
        with open(options.columns) as fin:
            cols = pickle.load(fin)
    else:
        print 'Using default cols'
        cols = ['light scan', 'heavy scan', 'pair configuration', 'M+H', 'score', 'seq', 'epsilon', 'ambiguous edges', 'num ambig edges']
    
    if options.output:
        outFile = open(options.output, 'w')
        outFile.write('\t'.join([col.upper() for col in cols]) + '\n')
    
    PNet = PN.ProbNetwork(options.config, options.model)

    dtaList = glob.glob(options.dtadir + '/*.dta')
    scanFDict = getScanFDict(dtaList)
    
    aas = Constants.addPepsToAADict(300)
    hashedAAs = Constants.hashAAsEpsilonRange(aas, epStep, maxEp)
    
    ambigOpenPenalty = 0
    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(options.minedge, ambigOpenPenalty, options.ambigpenalty)
    ppmPenaltyFun = DNS.getPPMPenaltyFun(options.ppmstd, hashedAAs, options.minedge, options.ppmpenalty, options.ppmsyserror, epStep)
    
    print 'Getting Clusters'
    parent = os.path.abspath(os.pardir)
    clusterSVMModel = svmutil.svm_load_model(parent + paramsDict['Cluster Configuration']['model'])
    clusterSVMRanges = svmutil.load_ranges(parent + os.path.splitext((paramsDict['Cluster Configuration']['model']))[0] + '.range')

    precMassClusters = Analytics.findSamePrecMassClusters(dtaList, ppm=options.ppmstd)
#    print 'precMassClusters', precMassClusters                                                                                                                                                                      
    samePeptideClusters = Analytics.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=float(paramsDict['Cluster Configuration']['cutoff']))
#    samePeptideClusters = Analytics.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=4)
#    samePeptideClusters = An.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=4)

    # To test without any clustering
    #samePeptideClusters = [[scanF] for scanF in scanFDict]
    
    for pairConfigName in paramsDict['Pair Configurations']:
        
        print 'Getting heavy-light pairs for %s' % (pairConfigName,)
        t1 = time.time()

        pairConfig = paramsDict['Pair Configurations'][pairConfigName]
        pairs = Analytics.findDeltaPairsClusters(samePeptideClusters, scanFDict, pairConfig['NMod']+pairConfig['CMod'], ppm=options.ppmstd)
        addEnds = DNS.getSpectrumGraphEndpointInitFunction(pairConfig['NStatic'], pairConfig['CStatic'], paramsDict['Enzyme']['specificity'])
        termModHash = Constants.getTermModHashForPairConfig(pairConfig)
        
        svmModel = svmutil.svm_load_model(parent + pairConfig['Model'])
        svmRange = svmutil.load_ranges(parent + os.path.splitext(pairConfig['Model'])[0] + '.range')
        
        xVals = []
        # xVals = getPairsThread(pairs)
        for pair in pairs:
            lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
            heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
            lightPrecMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))

            epSTD = options.ppmstd * 10 ** -6 * lightPrecMass

            lightMergedSpec = SA.mergeSpectra(lightSpecs, epsilon=2*epSTD)
            heavyMergedSpec = SA.mergeSpectra(heavySpecs, epsilon=2*epSTD)


            svmClassificationData = SA.getSpectraPairInfoForSVMClassification(lightMergedSpec, heavyMergedSpec, lightPrecMass, NMod=pairConfig['NMod'], CMod=pairConfig['CMod'], epsilon=2*epSTD)
            xVals += [svmClassificationData]
        
        
        xValsNorm = svmutil.normalize_instances(xVals, svmRange)
        pLab = svmutil.svm_predict([0]*len(xValsNorm), xValsNorm, svmModel)[0]
        
        print 'Pairs found. Time taken:', time.time() - t1, '\n'
        heavySeqMap = copy.deepcopy(seqMap['LADS Unit Test'])
        heavySeqMap['Mods']['N-Term'] = paramsDict['Pair Configurations'][pairConfigName]['NModSymbol']
        heavySeqMap['Mods']['C-Term'] = paramsDict['Pair Configurations'][pairConfigName]['CModSymbol']
        
#        hyperParameters = PNet.getHyperParameters(pairConfigName)
#        ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(hyperParameters['minedge'], hyperParameters['ambigopen'], hyperParameters['ambigextend'])
#        ppmPenaltyFun = DNS.getPPMPenaltyFun(hyperParameters['ppmstd'], hashedAAs, hyperParameters['minedge'], hyperParameters['ppmpen'], 0, epStep)

        getSequencingThread(pairs, xVals, paramsDict, outFile, cols, pLab)
        # for i, pair in enumerate(pairs):
        #     if pLab[i] == -1:
        #         continue
            
        #     scanData = {}
        #     lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
        #     heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
        #     precMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))
            
        #     epMean = options.ppmsyserror * precMass * 10**-6
        #     epSTD = options.ppmstd * precMass * 10**-6
                        
        #     scanData['shared peaks ratio'] = xVals[i][1]
        
        #     print 'Now sequencing light scan(s) %s, heavy scan(s) %s with shared peaks ratio %f \n' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), xVals[i][1])
        #     s1 = time.time()
        #     sharedInfo, starts, ends, deltas, G = DNS.prepPairedSpectrumGraph(lightSpecs, heavySpecs, precMass, addEnds, ppmSTD=options.ppmstd, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], verbose=options.verbose)
        #     scanData['M+H'] = precMass
            
        #     specs = []
        #     for massIntPairs in lightSpecs:
        #         specs += [PN.Spectrum(PNet, precMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=massIntPairs)]
        #     for massIntPairs in heavySpecs:
        #         specs += [PN.Spectrum(PNet, precMass + pairConfig['NMod'] + pairConfig['CMod'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD, spectrum=massIntPairs)]
        #     for spec in specs:
        #         spec.initializeNoiseModel()

        #     scanData.update(DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, termModHash=termModHash, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose))
        #     scanData['sequencing time'] = time.time() - s1
        #     scanData['pair configuration'] = pairConfigName
        #     if validateHeavySequence(scanData['seq'], heavySeqMap, scanData['ambiguous edges']):
        #         for scanF in samePeptideClusters[pair[0]] + samePeptideClusters[pair[1]]:
        #             scanFDict[scanF]['sequenced'] = True
        #         if options.output:
        #             for pair in [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[pair[0]] for heavyScanF in samePeptideClusters[pair[1]]]:
        #                 scanData['light scan'] = int(pair[0])
        #                 scanData['heavy scan'] = int(pair[1])
        #                 outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
        #     else:
        #         print 'WARNING: Invalid sequence! Unsuccessful sequencing of %s and %s with pair configuration %s' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), pairConfigName)
        #     print '\nTime Taken:', time.time() - s1, '\n'     

    if len(paramsDict['Pair Configurations'].keys()) == 0:
        unpaired = glob.glob(options.dtadir + '/*.dta')
    
    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])  
    termModHash = Constants.createTermModHashAAs(N=Constants.NTermMods, C=Constants.CTermMods) 

#    hyperParameters = PNet.getHyperParameters(None)
#    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(hyperParameters['minedge'], hyperParameters['ambigopen'], hyperParameters['ambigextend'])
#    ppmPenaltyFun = DNS.getPPMPenaltyFun(hyperParameters['ppmstd'], hashedAAs, hyperParameters['minedge'], hyperParameters['ppmpen'], 0, epStep)

    for cluster in samePeptideClusters:
        if scanFDict[cluster[0]]['sequenced'] or True:
            continue
        
        scanData = {}
        specsInfo = [DataFile.getMassIntPairs(scanFDict[scanF]['dta']) for scanF in cluster]
        precMass = np.average(np.array([scanFDict[scanF]['precMass'] for scanF in cluster]))
        
        scanData['light scan'] = None
        scanData['heavy scan'] = 'N/A'
        scanData['shared peaks ratio'] = 'N/A'
        
        print 'Now sequencing unpaired scan(s) %s \n' % (str(cluster),)
        s1 = time.time()
        
        sharedInfo, starts, ends, deltas, G = DNS.prepSingleSpectrumGraph(specsInfo, precMass, addEnds, options.ppmstd, 0, 0, verbose=options.verbose)
        scanData['M+H'] = precMass
        
        epMean = options.ppmsyserror * precMass * 10**-6
        epSTD = options.ppmstd * precMass * 10**-6
        specs = [PN.Spectrum(PNet, precMass, Nmod=0, Cmod=0, epsilon=2*epSTD, spectrum=massIntPairs) for massIntPairs in specsInfo]
        for spec in specs:
            spec.initializeNoiseModel()
        
        scanData.update(DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, termModHash=termModHash, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose))
        scanData['sequencing time'] = time.time() - s1
        scanData['pair configuration'] = 'N/A'
        if options.output:
            for scanF in cluster:
                scanData['light scan'] = int(scanF)
                outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

        print '\nTime Taken:', time.time() - s1, '\n'
    
    if options.output:    
        outFile.close()
    print 'Finished Sequencing. Total Time Taken:', time.time() - t1
