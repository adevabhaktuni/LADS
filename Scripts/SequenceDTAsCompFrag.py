'''
Created on Sep 7, 2011

@author: Arun
'''

import glob
import time
import pickle
import numpy as np
import copy

import ArgLib
import Constants
import DeNovoSequencer as DNS
import SpectraAlignment as SA
import Analytics
import DataFile
import ProbNetwork as PN

def getCIDETDPairs(scanFDict, ppm=5):
    pairs = []
    prevScanF = 0
    for scanF in sorted(scanFDict.keys()):
        #use a tolerance of 10 * ppm for precursor mass to cast a wide net
        if scanF - prevScanF == 1 and np.abs(scanFDict[scanF]['precMass'] - scanFDict[prevScanF]['precMass']) <  10**-5*ppm*scanFDict[scanF]['precMass'] and (not pairs or pairs[-1][1] != prevScanF):
            pairs += [(prevScanF, scanF)]
        prevScanF = scanF

    return pairs

def getSharedPeaksRatio(lightPath, heavyPath, epsilon):
    lightPairs = DataFile.getMassIntPairs(lightPath)
    heavyPairs = DataFile.getMassIntPairs(heavyPath)
    N, C = SA.getNandCIons(lightPairs, heavyPairs, 17.0265, -16.0187, epsilon=epsilon)
    return SA.getSharedPeaksRatio(lightPairs, heavyPairs, N, C)

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = DataFile.getScanNum(dta)
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'sequenced': False}
    
    return scanFDict


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
    
    t1 = time.time()
    print 'Configuring LADS for sequencing...'
    ETDPNet = PN.ProbNetwork(paramsDict['Models']['etd']['config'], paramsDict['Models']['etd']['model'])
    HCDPNet = PN.ProbNetwork(paramsDict['Models']['hcd']['config'], paramsDict['Models']['hcd']['model'])
    dtaList = glob.glob(options.dtadir + '/*.dta')
    scanFDict = getScanFDict(dtaList)
    aas = Constants.addPepsToAADict(options.minedge)
    hashedAAs = Constants.hashAAsEpsilonRange(aas, epStep, maxEp)

    ambigOpenPenalty = 0
    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(options.minedge, ambigOpenPenalty, options.ambigpenalty)
    ppmPenaltyFun = DNS.getPPMPenaltyFun(options.ppmstd, hashedAAs, options.minedge, options.ppmpenalty, options.ppmsyserror, epStep)

    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
    termModHash = Constants.createTermModHashAAs(N=copy.deepcopy(Constants.NTermMods), C=copy.deepcopy(Constants.CTermMods)) 
    print 'Getting Pairs...'
    pairs = getCIDETDPairs(scanFDict)

    for pair in pairs:    
        scanData = {}
        precMass = scanFDict[pair[0]]['precMass']
            
        epMean = options.ppmsyserror * precMass * 10**-6
        epSTD = options.ppmstd * precMass * 10**-6
            
        #(lightPath, heavyPath) = (scanFDict[pair[0]]['dta'], scanFDict[pair[1]]['dta'])
        #scanData['shared peaks ratio'] = getSharedPeaksRatio(lightPath, heavyPath, epSTD)
        scanData['shared peaks ratio'] = 'N/A'
        
        #skip if no common peaks are found
        #if scanData['shared peaks ratio'] == 0:
        #    continue
            
        scanData['light scan'] = pair[0]
        scanData['heavy scan'] = pair[1]

        lightSpecs = [DataFile.getMassIntPairs(scanFDict[pair[0]]['dta'])]
        heavySpecs = [DataFile.getMassIntPairs(scanFDict[pair[1]]['dta'])]
        
        print 'Now sequencing %i %i\n' % (pair[0], pair[1])
        s1 = time.time()
        sharedInfo, starts, ends, deltas, G = DNS.prepPairedSpectrumGraph(lightSpecs, heavySpecs, precMass, addEnds, ppmSTD=options.ppmstd, Nmod=17.0265, Cmod=-16.0187, verbose=options.verbose)
                    
        scanData['M+H'] = precMass
        lightSpec = PN.Spectrum(HCDPNet, precMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=lightSpecs[0])
        lightSpec.initializeNoiseModel()
        heavySpec = PN.Spectrum(ETDPNet, scanFDict[pair[1]]['precMass'], Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=heavySpecs[0])
        heavySpec.initializeNoiseModel()
        specs = [lightSpec, heavySpec]
        
        scanData.update(DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, termModHash=termModHash, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose))
        scanData['sequencing time'] = time.time() - s1
        scanData['pair configuration'] = 'N/A'
        scanFDict[pair[0]]['sequenced'] = True
        scanFDict[pair[1]]['sequenced'] = True
        if options.output:
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

        print '\nTime Taken:', time.time() - s1, '\n'   

    """
    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
    termModHash = Constants.createTermModHashAAs(N=copy.deepcopy(Constants.NTermMods), C=copy.deepcopy(Constants.CTermMods)) 
    for scanF in scanFDict:
        if scanFDict[scanF]['sequenced']:
            continue
        
        scanData = {}
        precMass = scanFDict[scanF]['precMass']
        
        scanData['light scan'] = scanF
        scanData['heavy scan'] = 'N/A'
        scanData['shared peaks ratio'] = 'N/A'
        
        print 'Now sequencing unpaired spectrum %s \n' % scanFDict[scanF]['dta']
        s1 = time.time()
        
        sharedInfo, starts, ends, deltas, G = DNS.prepSingleSpectrumGraph(scanFDict[scanF]['dta'], addEnds, options.ppmstd, 0, 0, options.verbose)
        precMass = sharedInfo['precMass']
        scanData['M+H'] = precMass
        
        epMean = options.ppmsyserror * precMass * 10**-6
        epSTD = options.ppmstd * precMass * 10**-6
        specs = [PN.Spectrum(PNet, precMass, Nmod=0, Cmod=0, epsilon=2*epSTD, spectrum=sharedInfo['pairs'])]
        specs[0].initializeNoiseModel()
        
        scanData.update(DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, termModHash=termModHash, maxEdge=options.maxedge, minEdge=options.minedge, subGraphCut=options.subgraphcut, subAlpha=0.3, alpha=options.alpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=options.verbose))
        scanData['sequencing time'] = time.time() - s1
        scanData['pair configuration'] = 'N/A'
        if options.output:
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
    
        print '\nTime Taken:', time.time() - s1, '\n'

    """
    
    if options.output:    
        outFile.close()
    print 'Finished Sequencing. Total Time Taken:', time.time() - t1
