'''
Created on Sep 7, 2011

@author: Arun
'''

import sys
import os

sys.path.insert(1, os.path.abspath(os.pardir))

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

if __name__ == '__main__' :
    options = ArgLib.parse(['init', 'dtadir', 'config', 'model', 'output', 'columns', 'verbose', 'paircutoff', 'ppmsyserror', 'ppmstd', 'ppmpenalty', 'ambigpenalty', 'minedge', 'maxedge', 'alpha', 'subgraphcut', 'symbolmap', 'modalpha', 'unimoddict'])
    epStep = 0.00025
    maxEp = 0.1

    # Any peptides generated using this dictionary will not have any static modifications
    with open(options.unimoddict, 'r') as fin:
        unimodDict = pickle.load(fin)
    modAAsDict = Constants.parseModifications(unimodDict)
    unimodPeptDict = Constants.getUnimodPeptDict(200, modAAsDict)
    hashedMods = Constants.hashUnimodModAAsEpsilonRange(unimodPeptDict, epStep, maxEp)
    modAAsDict['C-term'].update(modAAsDict['Anywhere'])
    modAAsDict['N-term'].update(modAAsDict['Anywhere'])

    paramsDict = ArgLib.parseInitFile(options.init, options)
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap({'LADS Unit Test': 'LADS'}, symbolMap, paramsDict)

    
    if options.columns:
        with open(options.columns) as fin:
            cols = pickle.load(fin)
    else:
        print 'Using default cols'
        cols = ['light scan', 'heavy scan', 'pair configuration', 'M+H', 'score', 'seq', 'epsilon', 'ambiguous edges', 'num ambig edges', 'mod names']
    
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
    
    clusterSVMModel = svmutil.svm_load_model(paramsDict['Cluster Configuration']['model'])
    clusterSVMRanges = svmutil.load_ranges(os.path.splitext((paramsDict['Cluster Configuration']['model']))[0] + '.range')

    precMassClusters = Analytics.findSamePrecMassClusters(dtaList, ppm=options.ppmstd)
#    print 'precMassClusters', precMassClusters                                                                                                                                                                      
    samePeptideClusters = Analytics.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=float(paramsDict['Cluster Configuration']['cutoff']))
#    samePeptideClusters = Analytics.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=4)
#    samePeptideClusters = An.getSamePeptideClusters(precMassClusters, scanFDict, clusterSVMModel, clusterSVMRanges, ppmSTD=options.ppmstd, cutOff=4)                                                   
    
    for pairConfigName in paramsDict['Pair Configurations']:
        
        print 'Getting heavy-light pairs for %s' % (pairConfigName,)
        t1 = time.time()

        pairConfig = paramsDict['Pair Configurations'][pairConfigName]
        pairs = Analytics.findDeltaPairsClusters(samePeptideClusters, scanFDict, pairConfig['NMod']+pairConfig['CMod'], ppm=options.ppmstd)
        addEnds = DNS.getSpectrumGraphEndpointInitFunction(pairConfig['NStatic'], pairConfig['CStatic'], paramsDict['Enzyme']['specificity'])
        termModHash = Constants.getTermModHashForPairConfig(pairConfig)
        
        svmModel = svmutil.svm_load_model(pairConfig['Model'])
        svmRange = svmutil.load_ranges(os.path.splitext(pairConfig['Model'])[0] + '.range')
        
        xVals = []
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

        for i, pair in enumerate(pairs):
            if pLab[i] == -1:
                continue
            
            scanData = {}
            lightSpecs = [DataFile.getMassIntPairs(scanFDict[lightScanF]['dta']) for lightScanF in samePeptideClusters[pair[0]]]
            heavySpecs = [DataFile.getMassIntPairs(scanFDict[heavyScanF]['dta']) for heavyScanF in samePeptideClusters[pair[1]]]
            precMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in samePeptideClusters[pair[0]]]))
            
            epMean = options.ppmsyserror * precMass * 10**-6
            epSTD = options.ppmstd * precMass * 10**-6
                        
            scanData['shared peaks ratio'] = xVals[i][1]
        
            print 'Now sequencing light scan(s) %s, heavy scan(s) %s with shared peaks ratio %f \n' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), xVals[i][1])
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

            scanData.update(DNS.getSpectrumGraphDataBlindMods(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, hashedMods, modAAsDict, termModHash=termModHash, modAlpha=options.modalpha, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose))
            scanData['sequencing time'] = time.time() - s1
            scanData['pair configuration'] = pairConfigName
            if validateHeavySequence(scanData['seq'], heavySeqMap, scanData['ambiguous edges']):
                for scanF in samePeptideClusters[pair[0]] + samePeptideClusters[pair[1]]:
                    scanFDict[scanF]['sequenced'] = True
                if options.output:
                    for pair in [(lightScanF, heavyScanF) for lightScanF in samePeptideClusters[pair[0]] for heavyScanF in samePeptideClusters[pair[1]]]:
                        scanData['light scan'] = int(pair[0])
                        scanData['heavy scan'] = int(pair[1])
                        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
            else:
                print 'WARNING: Invalid sequence! Unsuccessful sequencing of %s and %s with pair configuration %s' % (str(samePeptideClusters[pair[0]]), str(samePeptideClusters[pair[1]]), pairConfigName)
            print '\nTime Taken:', time.time() - s1, '\n'   

    if len(paramsDict['Pair Configurations'].keys()) == 0:
        unpaired = glob.glob(options.dtadir + '/*.dta')
    
    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])  
    termModHash = Constants.createTermModHashAAs(N=Constants.NTermMods, C=Constants.CTermMods) 

    if options.output:    
        outFile.close()
    print 'Finished Sequencing. Total Time Taken:', time.time() - t1

#    hyperParameters = PNet.getHyperParameters(None)
#    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(hyperParameters['minedge'], hyperParameters['ambigopen'], hyperParameters['ambigextend'])
#    ppmPenaltyFun = DNS.getPPMPenaltyFun(hyperParameters['ppmstd'], hashedAAs, hyperParameters['minedge'], hyperParameters['ppmpen'], 0, epStep)
"""
    for cluster in samePeptideClusters:
        if scanFDict[cluster[0]]['sequenced']:
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
""" 

