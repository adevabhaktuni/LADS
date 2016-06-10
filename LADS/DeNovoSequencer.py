'''
Created on Jun 1, 2011

@author: Arun
'''

import numpy as np
import Constants

import DataFile
import networkx as nx
from collections import deque
import copy
import SpectraAlignment as SA
import time
import ProbNetwork as PN

SGMass = Constants.aminoacids['S'][2] + Constants.aminoacids['G'][2]

def getSpectrumGraphEndpointInitFunction(startOffsets, endOffsets, specificity):
    minusFromCTerm = np.array([])
    addToNTerm = np.array([])
    for flankingPair in specificity:
        if flankingPair[0]:
            PRMladder = PN.getPRMLadder(flankingPair[0], startMass=0, addEnds=True, considerTerminalMods=False)
            for end in endOffsets:
                minusFromCTerm = np.append(minusFromCTerm, [end + PRM for PRM in PRMladder])
        else:
            minusFromCTerm = endOffsets
        
        if flankingPair[1]:
            PRMladder = PN.getPRMLadder(flankingPair[1], startMass=0, addEnds=True, considerTerminalMods=False)
            for start in startOffsets:
                addToNTerm = np.append(addToNTerm, [start + PRM for PRM in PRMladder])

        else:
            addToNTerm = startOffsets
    
    minusFromCTerm = np.unique(minusFromCTerm)
    addToNTerm = np.unique(addToNTerm)
    def addEndsToSpectrumGraph(G, PM):
        G.add_nodes_from(PM - minusFromCTerm)
        G.add_nodes_from(addToNTerm)
        
        return (startOffsets, PM - endOffsets)
    
    return addEndsToSpectrumGraph
        
def getAmbigEdgePenaltyFunction(minEdge, ambigEdgeOpenPenalty, ambigEdgeExtendPenalty):
    const = float(ambigEdgeExtendPenalty)/minEdge
    def penaltyFun(fromNode, toNode):
        return ambigEdgeOpenPenalty + const * (toNode-fromNode)
        #return ambigEdgePenalty
    return penaltyFun

def getPPMPenaltyFun(ppmSTD, hashedAAsWithEpRange, minEdge, maxPPMPenalty=15, ppmSysError=0, epStep=0.0005):
    const1 = maxPPMPenalty*(2*ppmSTD-ppmSTD)
    const2 = 2*ppmSTD-maxPPMPenalty*(ppmSTD)
    def PPMPenalty(epDiff, precMass):
        ppmDiff = np.abs(epDiff/precMass * 10**6 - ppmSysError)
        if ppmDiff < ppmSTD:
            return 0
        else:
            return maxPPMPenalty - const1/((maxPPMPenalty-1)*ppmDiff+const2)
    
    def penaltyFun(edgeMass, edgeSeq, precMass):
        hashedEdge = np.round(edgeMass/epStep)
        if edgeSeq == '-':
            try:
                edgeEpDiff = hashedAAsWithEpRange[hashedEdge]['min']
                return PPMPenalty(edgeEpDiff, precMass)
            except KeyError:
            # Penalize edges which aren't supposed to correspond to unmodified peptide sequences
                return maxPPMPenalty

        else:
            try:
                edgeEpDiff = hashedAAsWithEpRange[hashedEdge]['seqs'][edgeSeq]
                return PPMPenalty(edgeEpDiff, precMass)
            except KeyError:
                return maxPPMPenalty
                
    return penaltyFun
    

def getSubGraphScore(deltas, startMass, endMass, specs, alpha=0.9, epsilon=0.02):
    sDeltas = sliceDeltas(deltas, startMass, endMass)
    subG = getSpectrumGraph(startMass, endMass, sDeltas, maxNodeNum=100)
    labelEdgesWithSequence(subG, maxMassCutoff=500, minMassCutoff=minEdge)
    probNetScoreFunc(subG, specs, unknownPenalty=unknownPenalty)
    maxScore = insertLScores(subG)
    zeroScore = insertRScores(subG)
    print maxScore, zeroScore
    if np.abs(maxScore - zeroScore) < 1:
        subG = getSubOptimalSubgraph(subG, maxScore, alpha)
        clearLAndRScores(subG)
        print 'Resolving ambiguous edges'
        subG = resolveAmbiguousEdges(subG, sDeltas, specs, keepAmbigEdges=True)
    
    print 'Ambiguous Edges resolved'
    probNetScoreFunc(subG, specs, useExistingMemo=True, unknownPenalty=unknownPenalty)
    destScore = insertLScores(subG)
    oriScore = insertRScores(subG)
    print destScore, oriScore
    paths = getSubOptimalPaths(subG, alpha, endMass, oriScore)
    
    seqs = []
    for path in paths:
        seqs.extend([getSequenceFromNodes(subG, path[1])])
    
    scores = list(zip(*paths)[0])
    printSequences(seqs, scores, maxNum=50)

def sequenceSingleSpectrum(pairs, pepMass, PNet, unknownPenalty=5, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, Nmod=0.0, Cmod=0.0, alpha=0.95, ppm=5, numSeqs=10, aas=Constants.aminoacids, verbose=False):
    PM = pepMass + Constants.mods['H+'] + Constants.mods['H2O']
    epsilon = ppm * PM * 10 ** -6
    specs = [PN.Spectrum(PNet, PM, Nmod=Nmod, Cmod=Cmod, epsilon=epsilon, spectrum=pairs)]
    specs[0].initializeNoiseModel()
    
    deltas = getDeltaIons(pairs, pepMass, Nmod, Cmod)
    
    G = getSpectrumGraph(0, pepMass, deltas, maxNodeNum=50)
    G.add_node(pepMass - Constants.aminoacids['K'][2], intensity=0)
    G.add_node(pepMass - Constants.aminoacids['X'][2], intensity=0)
    G = mergeNodes(G, epsilon=2 * epsilon)
    
    scores, seqs = solveSpectrumGraph(G, deltas, pepMass, specs, unknownPenalty=unknownPenalty, maxEdge=maxEdge, minEdge=minEdge, subGraphCut=subGraphCut, subAlpha=subAlpha, alpha=alpha, epsilon=epsilon, aas=aas, verbose=verbose)
    if seqs:
        printSequences(seqs, scores, maxNum=10)
    else:
        print 'No Sequences Found'

def prepSingleSpectrumGraph(specsInfo, precMass, addEndsFn, ppmSTD=5, Nmod=0, Cmod=0, verbose=False, maxNodeNum=80):
    
    epSTD = ppmSTD * precMass * 10 ** -6
    mergedPairs = SA.mergeSpectra(specsInfo, 2*epSTD)
    
    NTermTable = SA.getNandCIons(mergedPairs, mergedPairs, Nmod=0, Cmod=0, epsilon=2*epSTD)[0]
    CCrossTable = SA.getCrossPairedIons(mergedPairs, mergedPairs, precMass, Nmod=0, Cmod=0, epsilon=2*epSTD)[1]
    NTermIonDict = SA.prepIonTableForAddition(NTermTable, ['b', 'b'])
    CCrossIonDict = SA.prepIonTableForAddition(CCrossTable, ['b', 'y'])

    symIonsDict = SA.reverseDict(SA.addDicts(NTermIonDict, CCrossIonDict))
    
    sharedInfo = {'precMass': precMass, 'pairs': mergedPairs}
    PM = precMass - Constants.mods['H+'] - Constants.mods['H2O']
    
    
    deltas = getUnpairedRankedDeltas(mergedPairs, PM, symIonsDict)
    G = getSpectrumGraph(0, PM, deltas, maxNodeNum=maxNodeNum, addEndMasses=False)
    starts, ends = addEndsFn(G, PM)
    return sharedInfo, starts, ends, deltas, mergeNodes(G, epsilon=2 * epSTD )
       
        
def sequencePairedSpectra(NInd, CInd, lightPairs, heavyPairs, pepMass, PNet, Nmod, Cmod, unknownPenalty=5, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, numSeqs=10, alpha=0.95, ppm=5, aas=Constants.aminoacids, verbose=False):
    G = nx.DiGraph()
    maxInt = max(lightPairs[:, 1].max(), heavyPairs[:, 1].max())
    G.add_node(0, intensity=maxInt)
    G.add_node(pepMass, intensity=maxInt)
    
    PM = pepMass + Constants.mods['H+'] + Constants.mods['H2O']
    epsilon = ppm * PM * 10 ** -6
    lightSpec = PN.Spectrum(PNet, PM, Nmod=0.0, Cmod=0.0, epsilon=epsilon, spectrum=lightPairs)
    lightSpec.initializeNoiseModel()
    heavySpec = PN.Spectrum(PNet, PM + Nmod + Cmod, Nmod=Nmod, Cmod=Cmod, epsilon=epsilon, spectrum=heavyPairs)
    heavySpec.initializeNoiseModel()
    specs = [lightSpec, heavySpec]
    
    CTermIons = lightPairs[CInd[0]]
    CTermIons[:, 0] = PM - CTermIons[:, 0]
    addMassesToSpectrumGraph(G, CTermIons, pepMass)

    NTermIons = lightPairs[NInd[0]]
    NTermIons[:, 0] = NTermIons[:, 0] - Constants.mods['H+']
    addMassesToSpectrumGraph(G, NTermIons, pepMass)
    
    G.add_node(pepMass - Constants.aminoacids['K'][2], intensity=0)
    
    deltas = prepareDeltas(NInd, CInd, lightPairs, heavyPairs, pepMass)
    
    G = mergeNodes(G, epsilon=2 * epsilon)
    
    scores, seqs = solveSpectrumGraph(G, deltas, pepMass, specs, unknownPenalty=unknownPenalty, maxEdge=maxEdge, minEdge=minEdge, subGraphCut=subGraphCut, subAlpha=subAlpha, alpha=alpha, epsilon=epsilon, aas=aas, verbose=verbose)
    if seqs:
        printSequences(seqs, scores, maxNum=10)
    else:
        print 'No Sequences Found'

def getInitialNumPairedNodes(allPairedIonsDict, defaultNum=50):
    if len(allPairedIonsDict) == 0:
        return defaultNum

    # only put nodes with b and y ions present in both spectra into initial graph
    numNodes = 0
    for ion in allPairedIonsDict:
        if len(ion) == 2 and len(allPairedIonsDict[ion]) == 2:
            numNodes += 1
            
    return numNodes if numNodes > 0 else len(allPairedIonsDict)

def prepPairedSpectrumGraphAlt(lightSpecs, heavySpecs, precMass, addEndsFn, Nmod, Cmod, ppmSTD=5, maxNodeNum=80, verbose=False):
    epSTD = ppmSTD * precMass * 10 ** -6^M
                                          
    
    #Minimum number of nodes in graph in precursor mass divided by average molecular weight of amino acid
    sharedInfo = getPairedSpectraInfoForSequencing(lightSpecs, heavySpecs, precMass, Nmod=Nmod, Cmod=Cmod, epsilon=2*epSTD, verbose=verbose)^M
    
    PM = precMass - Constants.mods['H+'] - Constants.mods['H2O']
    
    deltas = getPairedRankedDeltas(sharedInfo['lightPairs'], sharedInfo['heavyPairs'], sharedInfo['pairedIonsDict'], PM, Nmod, Cmod)
    #    print deltas
    #    print maxNodeNum, len(sharedInfo['pairedIonsDict'])
    G = getSpectrumGraph(0, PM, deltas, maxNodeNum=min(getInitialNumPairedNodes(sharedInfo['pairedIonsDict'], defaultNum=maxNodeNum), maxNodeNum), addEndMasses=False)
    starts, ends = addEndsFn(G, PM)
                                   
    return sharedInfo, starts, ends, deltas, mergeNodes(G, epsilon=2 * epSTD)
    


def prepPairedSpectrumGraph(lightSpecs, heavySpecs, precMass, addEndsFn, Nmod, Cmod, ppmSTD=5, maxNodeNum=50, verbose=False): 
    epSTD = ppmSTD * precMass * 10 ** -6
    

    #Minimum number of nodes in graph in precursor mass divided by average molecular weight of amino acid
    sharedInfo = getPairedSpectraInfoForSequencing(lightSpecs, heavySpecs, precMass, Nmod=Nmod, Cmod=Cmod, epsilon=2*epSTD, verbose=verbose)
    
    PM = precMass - Constants.mods['H+'] - Constants.mods['H2O'] 

    deltas = getPairedRankedDeltas(sharedInfo['lightPairs'], sharedInfo['heavyPairs'], sharedInfo['pairedIonsDict'], PM, Nmod, Cmod)
#    print deltas
#    print maxNodeNum, len(sharedInfo['pairedIonsDict'])
    G = getSpectrumGraph(0, PM, deltas, maxNodeNum=min(max(np.floor(precMass/110), len(sharedInfo['pairedIonsDict'])), maxNodeNum), addEndMasses=False)
    starts, ends = addEndsFn(G, PM)
    
    return sharedInfo, starts, ends, deltas, mergeNodes(G, epsilon=2 * epSTD)

                                                                                
def appendToDeltas(newDeltas, deltas=np.transpose(np.array([[],[]]))):
    newDeltas = np.array(newDeltas)
    try:
        newDeltas = newDeltas[np.argsort(-1*newDeltas[:,1])]
        newDeltas = newDeltas[np.argsort(-1*newDeltas[:,2], kind='mergesort')]
#        print 'deltas', deltas
#        print 'newdeltas', newDeltas
        return np.append(deltas, np.array(zip(newDeltas[:,0], range(deltas.shape[0], deltas.shape[0]+newDeltas.shape[0]))), axis=0)
    except IndexError:
        return deltas
    
def getUnpairedBAndYDeltas(ions, PM, NMod=0, CMod=0):
    try:
        bDeltas = copy.copy(ions)
        bDeltas[:,0] = bDeltas[:,0] - Constants.mods['H+']
        yDeltas = copy.copy(ions)
        yDeltas[:,0] = PM + Constants.mods['H2O'] + Constants.mods['H+'] - yDeltas[:,0]
        return np.append(bDeltas, yDeltas, axis=0)
    except IndexError:
        return np.transpose(np.array([],[],[]))
    
def getPairedRankedDeltas(lightMergedSpec, heavyMergedSpec, allPairedIonsDict, lightPM, NMod, CMod):
    pairedLightInds = set()
    pairedHeavyInds = set()
    doubleSymDeltas = []
    singleSymDeltas = []
    pairedIonDeltas = []
        
    for heavyIons in allPairedIonsDict:
        deltaMass = 0
        deltaIntensity = 0
        deltaCount = 0
#        print heavyIons, allPairedIonsDict[heavyIons]
        for ion in allPairedIonsDict[heavyIons]:
#            print 'light', ion, lightMergedSpec[ion[0]][0]
            pairedLightInds.add(ion[0])
            deltaMass += PN.ProbNetwork.deltaRules[ion[1]](lightPM, lightMergedSpec[ion[0]][0], 0, 0)
            deltaIntensity += lightMergedSpec[ion[0]][1]
            deltaCount += lightMergedSpec[ion[0]][2]
        for ion in heavyIons:
#            print 'heavy', ion, heavyMergedSpec[ion[0]][0] 
            pairedHeavyInds.add(ion[0])
            deltaMass += PN.ProbNetwork.deltaRules[ion[1]](lightPM, heavyMergedSpec[ion[0]][0], NMod, CMod)
            deltaIntensity += heavyMergedSpec[ion[0]][1]
            deltaCount += heavyMergedSpec[ion[0]][2]
            
        if len(heavyIons) == 2 and len(allPairedIonsDict[heavyIons]) == 2:
#            print 'doubleSym'
            doubleSymDeltas += [[deltaMass/4, deltaIntensity, deltaCount]]
        elif len(heavyIons) == 2 or len(allPairedIonsDict[heavyIons]) == 2:
#            print 'singleSym'
            singleSymDeltas += [[deltaMass/3, deltaIntensity, deltaCount]]
        else:
#            print 'paired'
            pairedIonDeltas += [[deltaMass/2, deltaIntensity, deltaCount]]
    
    deltas = appendToDeltas(doubleSymDeltas)
#    print deltas
    deltas = appendToDeltas(singleSymDeltas, deltas)
#    print deltas
    deltas = appendToDeltas(pairedIonDeltas, deltas)
#    print deltas
    
    unpairedLightIons = lightMergedSpec[list(set(range(lightMergedSpec.shape[0])) - pairedLightInds)]
    unpairedLightDeltas = getUnpairedBAndYDeltas(unpairedLightIons, lightPM, NMod=0, CMod=0)
    
    unpairedHeavyIons = heavyMergedSpec[list(set(range(heavyMergedSpec.shape[0])) - pairedHeavyInds)]
    unpairedHeavyDeltas = getUnpairedBAndYDeltas(unpairedHeavyIons, lightPM, NMod=NMod, CMod=CMod)
    
    deltas = appendToDeltas(np.append(unpairedLightDeltas, unpairedHeavyDeltas, axis=0), deltas)
    return deltas[np.argsort(deltas[:,0])]
    
def getUnpairedRankedDeltas(mergedSpec, PM, symIonsDict):
    symInds = set()
    
    symDeltas = []
    for ions in symIonsDict:
        if len(ions) == 2:
#            print ions, symIonsDict[ions], mergedSpec[ions[0][0]], mergedSpec[ions[1][0]]
            deltaMass = 0
            deltaIntensity = 0
            deltaCount = 0
            for ion in ions:
                symInds.add(ion[0])
#                print PN.ProbNetwork.deltaRules[ion[1]](PM, mergedSpec[ion[0]][0], 0, 0)
                deltaMass += PN.ProbNetwork.deltaRules[ion[1]](PM, mergedSpec[ion[0]][0], 0, 0)
                deltaIntensity += mergedSpec[ion[0]][1]
                deltaCount += mergedSpec[ion[0]][2]
            
#            print [deltaMass/2, deltaIntensity, deltaCount]
            symDeltas += [[deltaMass/2, deltaIntensity, deltaCount]]
    
    deltas = appendToDeltas(symDeltas)
    unpairedIons = mergedSpec[list(set(range(mergedSpec.shape[0])) - symInds)]
    unpairedDeltas = getUnpairedBAndYDeltas(unpairedIons, PM)
    
    deltas = appendToDeltas(unpairedDeltas, deltas)
    return deltas[np.argsort(deltas[:,0])]
    
def adjustStartsAndEnds(G, starts, ends):
    for i, start in enumerate(starts):
        if start not in G.node:
            adjNodeInd = np.argmin(np.abs(np.array(G.nodes())-start))
            starts[i] = G.nodes()[adjNodeInd]
    
    for i, end in enumerate(ends):
        if end not in G.node:
            adjNodeInd = np.argmin(np.abs(np.array(G.nodes())-end))
            ends[i] = G.nodes()[adjNodeInd]

def getSpectrumGraphData(G, deltas, specs, starts, ends, PM, ambigPenaltyFun, ppmPenaltyFun, hashedAAs, epMean=0, epSTD=0.02, termModHash={'NTerm': {}, 'CTerm': {}}, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epStep=0.0005, maxNumPaths=10, verbose=False):
    graphData = {}
    adjustStartsAndEnds(G, starts, ends)
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    seqEdges = 0
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['seq'] != '-':
            seqEdges += 1
    graphData['AA edges'] = seqEdges
    graphData['start nodes'] = len(G.nodes())
    graphData['start edges'] = len(G.edges())
    graphData['epsilon'] = epSTD
      
    k = 0
    while not pathExists(G, starts, ends):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    probNetScoreFunc(G, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
#    print '\nInitial spectrum graph edges: ', sorted(G.edges())
    #pairs = getSymPairs(G, PM, epsilon)
    #print [pair for pair in sorted(pairs.keys())]
    matG = getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=epSTD, symPenalty=50, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    #print 'Done getting matG'
    maxScore = insertLScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    zeroScore = insertRScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    #print matGmaxScore, matGzeroScore
    subG = getSubOptimalSubgraphMatG(matG, G, maxScore, matG.graph['ends'], starts, ends, alpha=alpha)
#    print 'Suboptimal subgraph nodes from matG: ', str(sorted(subG.nodes()))
#    print 'Suboptimal subgraph edges from matG: ', str(sorted(subG.edges()))
    
    #intensityScoreFunc(G, lightSpec, heavySpec)
    #maxScore = insertLScores(G, starts, ends)
    #zeroScore = insertRScores(G, starts, ends)
    #listEdges(G)
    #listNodes(G)
    #print maxScore, zeroScore
    #subG = getSubOptimalSubgraph(G, maxScore, alpha)
    if verbose:
        print '\nInitial spectrum graph nodes:', sorted(G.nodes()) 
        print '\nSuboptimal Subgraph nodes at alpha factor %f using lScore %f and rScore %f: ' % (alpha, maxScore, zeroScore) + str(sorted(subG.nodes()))
        print '\nSuboptimal Subgraph edges: ', str(sorted(subG.edges()))
        print '\nResolving Ambiguous Edges'
    
    subG = resolveAmbiguousEdges(subG, deltas, specs, PM, hashedAAs, subGraphCut=subGraphCut, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
    #print 'After resolving and merging:', sorted(subG.nodes())
    #print 'Edges:', sorted(subG.edges())
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    clearLAndRScores(subG)
#    print sorted(subG.edges())
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)
    
    graphData['end nodes'] = len(subG.nodes())
    graphData['end edges'] = len(subG.edges())
    graphData['max score from ori'] = oriScore
    graphData['max score from dest'] = destScore

    prevNumPaths = 0
    pathAlphas = [0.90, 0.95, 0.99]
    paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    while (prevNumPaths != len(paths) and len(paths) < maxNumPaths and pathAlphas):
        prevNumPaths = len(paths)
        paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    
    graphData['num paths'] = len(paths)
    graphData['alpha'] = alpha
    if paths:
        seqs = []
        for path in paths:
            seqs.extend([getSequenceFromNodes(subG, path[1], PM, termModHash)])
    
        scores = list(zip(*paths)[0])

        printSequences(seqs, scores, maxNum=maxNumPaths)
        
        Ord = np.argsort(-1 * np.array(scores))
        #for i in range(10):
        #    print paths[Ord[i]]
        graphData['score'] = scores[Ord[0]]
        
        ambigEdges = []
        numAmbig = 0
        for i in range(len(seqs[Ord[0]])):
            if type(seqs[Ord[0]][i]) == tuple:
                ambigEdges.extend([seqs[Ord[0]][i]])
                numAmbig += 1
                seqs[Ord[0]][i] = '-'

        graphData['seq'] = ''.join(seqs[Ord[0]])
        graphData['ambiguous edges'] = ambigEdges
        graphData['num ambig edges'] = numAmbig
    
    else:
        graphData['seq'] = '-'
        graphData['score'] = '-'
        graphData['ambiguous edges'] = 'N/A'
        graphData['num ambig edges'] = 'N/A'
        print 'No Sequences Found'
        
    return graphData


def getSpectrumGraphDataBlindMods(G, deltas, specs, starts, ends, PM, ambigPenaltyFun, ppmPenaltyFun, hashedAAs, hashedMods, modAAsDict, epMean=0, epSTD=0.02, termModHash={'NTerm': {}, 'CTerm': {}}, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epStep=0.0005, maxNumPaths=10, verbose=False, modAlpha=0.99):
    graphData = {}
    adjustStartsAndEnds(G, starts, ends)
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    seqEdges = 0
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['seq'] != '-':
            seqEdges += 1
    graphData['AA edges'] = seqEdges
    graphData['start nodes'] = len(G.nodes())
    graphData['start edges'] = len(G.edges())
    graphData['epsilon'] = epSTD
      
    k = 0
    while not pathExists(G, starts, ends):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    probNetScoreFunc(G, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
#    print '\nInitial spectrum graph edges: ', sorted(G.edges())
    #pairs = getSymPairs(G, PM, epsilon)
    #print [pair for pair in sorted(pairs.keys())]
    matG = getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=epSTD, symPenalty=30, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    #print 'Done getting matG'
    maxScore = insertLScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    zeroScore = insertRScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    #print matGmaxScore, matGzeroScore
    subG = getSubOptimalSubgraphMatG(matG, G, maxScore, matG.graph['ends'], starts, ends, alpha=alpha)
#    print 'Suboptimal subgraph nodes from matG: ', str(sorted(subG.nodes()))
#    print 'Suboptimal subgraph edges from matG: ', str(sorted(subG.edges()))
    
    #intensityScoreFunc(G, lightSpec, heavySpec)
    #maxScore = insertLScores(G, starts, ends)
    #zeroScore = insertRScores(G, starts, ends)
    #listEdges(G)
    #listNodes(G)
    #print maxScore, zeroScore
    #subG = getSubOptimalSubgraph(G, maxScore, alpha)
    if verbose:
        print '\nInitial spectrum graph nodes:', sorted(G.nodes()) 
        print '\nSuboptimal Subgraph nodes at alpha factor %f using lScore %f and rScore %f: ' % (alpha, maxScore, zeroScore) + str(sorted(subG.nodes()))
        print '\nSuboptimal Subgraph edges: ', str(sorted(subG.edges()))
        print '\nResolving Ambiguous Edges'
    
    subG = resolveAmbiguousEdges(subG, deltas, specs, PM, hashedAAs, subGraphCut=subGraphCut, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
#    print 'After resolving and merging:', sorted(subG.nodes())
    #print 'Edges:', sorted(subG.edges())
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    clearLAndRScores(subG)
#    print sorted(subG.edges())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)
    
    subG = getSubOptimalSubgraph(subG, oriScore, modAlpha)
#    print 'After modalpha cutoff', sorted(subG.nodes())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    
    subG = resolveAmbiguousEdgesBlindMods(subG, deltas, specs, PM, hashedMods, subGraphCut=200, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
#    print 'After resolve ambiguous edges blind mods', sorted(subG.nodes())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    
    clearLAndRScores(subG)
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)

    graphData['end nodes'] = len(subG.nodes())
    graphData['end edges'] = len(subG.edges())
    graphData['max score from ori'] = oriScore
    graphData['max score from dest'] = destScore

    prevNumPaths = 0
    pathAlphas = [0.90, 0.95, 0.99]
    paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    while (prevNumPaths != len(paths) and len(paths) < maxNumPaths and pathAlphas):
        prevNumPaths = len(paths)
        paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    
    graphData['num paths'] = len(paths)
    graphData['alpha'] = alpha
    if paths:
        seqs = []
        for path in paths:
            seqs.extend([getSequenceFromNodesBlindMods(subG, path[1], PM, hashedMods, modAAsDict, termModHash, epStep=epStep, epsilon=0.02)])
    
        scores = list(zip(*paths)[0])

        printSequencesBlindMods(seqs, scores, maxNum=maxNumPaths)
        
        Ord = np.argsort(-1 * np.array(scores))
        #for i in range(10):
        #    print paths[Ord[i]]
        graphData['score'] = scores[Ord[0]]
        
        ambigEdges = []
        numAmbig = 0
        for i in range(len(seqs[Ord[0]][0])):
            if type(seqs[Ord[0]][0][i]) == tuple:
                ambigEdges.extend([seqs[Ord[0]][0][i]])
                numAmbig += 1
                seqs[Ord[0]][0][i] = '-'
            elif seqs[Ord[0]][0][i][0] == 'X':
                ambigEdges.extend([(paths[Ord[0]][1][i], paths[Ord[0]][1][i+1])])

 
        graphData['seq'] = ''.join(seqs[Ord[0]][0])
        graphData['ambiguous edges'] = ambigEdges
        graphData['num ambig edges'] = numAmbig
        graphData['mod names'] = seqs[Ord[0]][1]
    
    else:
        graphData['seq'] = '-'
        graphData['score'] = '-'
        graphData['ambiguous edges'] = 'N/A'
        graphData['num ambig edges'] = 'N/A'
        graphData['mod names'] = 'N/A'
        print 'No Sequences Found'
        
    return graphData

def getSpectrumGraphDataThread(G, deltas, specs, starts, ends, PM, ambigPenaltyFun, ppmPenaltyFun, hashedAAs, epMean=0, epSTD=0.02, termModHash={'NTerm': {}, 'CTerm': {}}, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epStep=0.0005, maxNumPaths=10, verbose=False):
    graphData = {}
    adjustStartsAndEnds(G, starts, ends)
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    seqEdges = 0
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['seq'] != '-':
            seqEdges += 1
    graphData['AA edges'] = seqEdges
    graphData['start nodes'] = len(G.nodes())
    graphData['start edges'] = len(G.edges())
    graphData['epsilon'] = epSTD
      
    k = 0
    while not pathExists(G, starts, ends):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    probNetScoreFunc(G, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
#    print '\nInitial spectrum graph edges: ', sorted(G.edges())
    #pairs = getSymPairs(G, PM, epsilon)
    #print [pair for pair in sorted(pairs.keys())]
    matG = getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=epSTD, symPenalty=50, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    #print 'Done getting matG'
    maxScore = insertLScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    zeroScore = insertRScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    #print matGmaxScore, matGzeroScore
    subG = getSubOptimalSubgraphMatG(matG, G, maxScore, matG.graph['ends'], starts, ends, alpha=alpha)
#    print 'Suboptimal subgraph nodes from matG: ', str(sorted(subG.nodes()))
#    print 'Suboptimal subgraph edges from matG: ', str(sorted(subG.edges()))
    
    #intensityScoreFunc(G, lightSpec, heavySpec)
    #maxScore = insertLScores(G, starts, ends)
    #zeroScore = insertRScores(G, starts, ends)
    #listEdges(G)
    #listNodes(G)
    #print maxScore, zeroScore
    #subG = getSubOptimalSubgraph(G, maxScore, alpha)
    if verbose:
        print '\nInitial spectrum graph nodes:', sorted(G.nodes()) 
        print '\nSuboptimal Subgraph nodes at alpha factor %f using lScore %f and rScore %f: ' % (alpha, maxScore, zeroScore) + str(sorted(subG.nodes()))
        print '\nSuboptimal Subgraph edges: ', str(sorted(subG.edges()))
        print '\nResolving Ambiguous Edges'
    
    subG = resolveAmbiguousEdges(subG, deltas, specs, PM, hashedAAs, subGraphCut=subGraphCut, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
    #print 'After resolving and merging:', sorted(subG.nodes())
    #print 'Edges:', sorted(subG.edges())
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    clearLAndRScores(subG)
#    print sorted(subG.edges())
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)
    
    graphData['end nodes'] = len(subG.nodes())
    graphData['end edges'] = len(subG.edges())
    graphData['max score from ori'] = oriScore
    graphData['max score from dest'] = destScore

    prevNumPaths = 0
    pathAlphas = [0.90, 0.95, 0.99]
    paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    while (prevNumPaths != len(paths) and len(paths) < maxNumPaths and pathAlphas):
        prevNumPaths = len(paths)
        paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    
    graphData['num paths'] = len(paths)
    graphData['alpha'] = alpha
    if paths:
        seqs = []
        for path in paths:
            seqs.extend([getSequenceFromNodes(subG, path[1], PM, termModHash)])
    
        scores = list(zip(*paths)[0])

        ret_print = []
        ret_print = printSequencesThread(seqs, scores, maxNum=maxNumPaths)
        graphData['over_scores'] = ret_print[0]
        graphData['blind'] = 0

        Ord = np.argsort(-1 * np.array(scores))
        #for i in range(10):
        #    print paths[Ord[i]]
        graphData['score'] = scores[Ord[0]]
        
        ambigEdges = []
        numAmbig = 0
        for i in range(len(seqs[Ord[0]])):
            if type(seqs[Ord[0]][i]) == tuple:
                ambigEdges.extend([seqs[Ord[0]][i]])
                numAmbig += 1
                seqs[Ord[0]][i] = '-'

        
        graphData['seq'] = ''.join(seqs[Ord[0]])
        graphData['ambiguous edges'] = ambigEdges
        graphData['num ambig edges'] = numAmbig
    
    else:
        graphData['seq'] = '-'
        graphData['score'] = '-'
        graphData['ambiguous edges'] = 'N/A'
        graphData['num ambig edges'] = 'N/A'
        print 'No Sequences Found'
        
    return graphData, ret_print


def getSpectrumGraphDataBlindModsThread(G, deltas, specs, starts, ends, PM, ambigPenaltyFun, ppmPenaltyFun, hashedAAs, hashedMods, modAAsDict, epMean=0, epSTD=0.02, termModHash={'NTerm': {}, 'CTerm': {}}, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epStep=0.0005, maxNumPaths=10, verbose=False, modAlpha=0.99):
    graphData = {}
    adjustStartsAndEnds(G, starts, ends)
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    seqEdges = 0
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['seq'] != '-':
            seqEdges += 1
    graphData['AA edges'] = seqEdges
    graphData['start nodes'] = len(G.nodes())
    graphData['start edges'] = len(G.edges())
    graphData['epsilon'] = epSTD
      
    k = 0
    while not pathExists(G, starts, ends):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epMean=epMean, epSTD=epSTD, epStep=epStep, hashedAAs=hashedAAs)
    
    probNetScoreFunc(G, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
#    print '\nInitial spectrum graph edges: ', sorted(G.edges())
    #pairs = getSymPairs(G, PM, epsilon)
    #print [pair for pair in sorted(pairs.keys())]
    matG = getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=epSTD, symPenalty=30, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    #print 'Done getting matG'
    maxScore = insertLScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    zeroScore = insertRScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    #print matGmaxScore, matGzeroScore
    subG = getSubOptimalSubgraphMatG(matG, G, maxScore, matG.graph['ends'], starts, ends, alpha=alpha)
#    print 'Suboptimal subgraph nodes from matG: ', str(sorted(subG.nodes()))
#    print 'Suboptimal subgraph edges from matG: ', str(sorted(subG.edges()))
    
    #intensityScoreFunc(G, lightSpec, heavySpec)
    #maxScore = insertLScores(G, starts, ends)
    #zeroScore = insertRScores(G, starts, ends)
    #listEdges(G)
    #listNodes(G)
    #print maxScore, zeroScore
    #subG = getSubOptimalSubgraph(G, maxScore, alpha)
    if verbose:
        print '\nInitial spectrum graph nodes:', sorted(G.nodes()) 
        print '\nSuboptimal Subgraph nodes at alpha factor %f using lScore %f and rScore %f: ' % (alpha, maxScore, zeroScore) + str(sorted(subG.nodes()))
        print '\nSuboptimal Subgraph edges: ', str(sorted(subG.edges()))
        print '\nResolving Ambiguous Edges'
    
    subG = resolveAmbiguousEdges(subG, deltas, specs, PM, hashedAAs, subGraphCut=subGraphCut, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
#    print 'After resolving and merging:', sorted(subG.nodes())
    #print 'Edges:', sorted(subG.edges())
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    clearLAndRScores(subG)
#    print sorted(subG.edges())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)
    
    subG = getSubOptimalSubgraph(subG, oriScore, modAlpha)
#    print 'After modalpha cutoff', sorted(subG.nodes())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    
    subG = resolveAmbiguousEdgesBlindMods(subG, deltas, specs, PM, hashedMods, subGraphCut=200, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epMean=epMean, epSTD=epSTD, epStep=epStep, verbose=verbose)
#    print 'After resolve ambiguous edges blind mods', sorted(subG.nodes())
#    for edge in sorted(subG.edges()):
#        print edge[0], edge[1], subG.edge[edge[0]][edge[1]]['seq']
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun)
    
    clearLAndRScores(subG)
    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)

    graphData['end nodes'] = len(subG.nodes())
    graphData['end edges'] = len(subG.edges())
    graphData['max score from ori'] = oriScore
    graphData['max score from dest'] = destScore

    prevNumPaths = 0
    pathAlphas = [0.90, 0.95, 0.99]
    paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    while (prevNumPaths != len(paths) and len(paths) < maxNumPaths and pathAlphas):
        prevNumPaths = len(paths)
        paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    
    graphData['num paths'] = len(paths)
    graphData['alpha'] = alpha
    if paths:
        seqs = []
        for path in paths:
            seqs.extend([getSequenceFromNodesBlindMods(subG, path[1], PM, hashedMods, modAAsDict, termModHash, epStep=epStep, epsilon=0.02)])
    
        scores = list(zip(*paths)[0])

        ret_print = []
        ret_print = printSequencesBlindModsThread(seqs, scores, maxNum=maxNumPaths)
        graphData['over_scores'] = ret_print[0]
        graphData['blind'] = 1

        Ord = np.argsort(-1 * np.array(scores))
        #for i in range(10):
        #    print paths[Ord[i]]
        graphData['score'] = scores[Ord[0]]
        
        ambigEdges = []
        numAmbig = 0
        for i in range(len(seqs[Ord[0]][0])):
            if type(seqs[Ord[0]][0][i]) == tuple:
                ambigEdges.extend([seqs[Ord[0]][0][i]])
                numAmbig += 1
                seqs[Ord[0]][0][i] = '-'
            elif seqs[Ord[0]][0][i][0] == 'X':
                ambigEdges.extend([(paths[Ord[0]][1][i], paths[Ord[0]][1][i+1])])
        
        
        graphData['seq'] = ''.join(seqs[Ord[0]][0])
        graphData['ambiguous edges'] = ambigEdges
        graphData['num ambig edges'] = numAmbig
        graphData['mod names'] = seqs[Ord[0]][1]
        
    
    else:
        graphData['seq'] = '-'
        graphData['score'] = '-'
        graphData['ambiguous edges'] = 'N/A'
        graphData['num ambig edges'] = 'N/A'
        graphData['mod names'] = 'N/A'
        print 'No Sequences Found'
        
    return graphData, ret_print

#Returns the index location of a mass in the matrix spectrum graph
#x vals (first half of mass list) are denoted by 0
#y vals (second half of mass list) are denoted by 1
def getSymPairIndexer(sortedSymList):
    x, y = zip(*sortedSymList)
    indHolder = {}
    for i, mass in enumerate(x):
        indHolder[mass] = [i, 0]
    for i, mass in enumerate(y):
        indHolder[mass] = [i, 1]
    
    return indHolder
    

#assumes that all starts occur in first half of peptide mass, and all ends occur in last half
def getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=0.02, symPenalty=30, ambigPenaltyFun=lambda fromNode, toNode: 20, ppmPenaltyFun=lambda edgeMass, edgeSeq, precMass: 0):
    matG = nx.DiGraph()
    precMass = PM + Constants.mods['H2O'] + Constants.mods['H+']
    symPairs = getSymPairs(G, PM, epsilon)
    sortedSymList = sorted(symPairs.keys())
    symListSize = len(sortedSymList)
    indHolder = getSymPairIndexer(sortedSymList)
    
    for i in range(len(sortedSymList)):
        for j in range(len(sortedSymList)):
            matG.add_node((i, j), penalty=0)
    
    matG.graph['symList'] = {}
    matG.graph['symList'][0], matG.graph['symList'][1] = zip(*sortedSymList)
    termNodes = []
    for edge in G.edges():
        indFrom = indHolder[edge[0]]
        indTo = indHolder[edge[1]]
        if indFrom[1] == 1 and indTo[1] == 1:
            for xInd in range(indFrom[0]+1):
                if symPairs[sortedSymList[xInd]][0]:
                    matG.add_edge((xInd, indTo[0]), (xInd, indFrom[0]), seq=G.edge[edge[0]][edge[1]]['seq'])
                    matG.edge[(xInd, indTo[0])][(xInd, indFrom[0])]['penalty'] = (ambigPenaltyFun(edge[0], edge[1]) if G.edge[edge[0]][edge[1]]['seq'] == '-' else 0) + ppmPenaltyFun(edge[1]-edge[0], G.edge[edge[0]][edge[1]]['seq'], precMass)
        elif indFrom[1] == 0 and indTo[1] == 0:
            for yInd in range(indTo[0]+1):
                if symPairs[sortedSymList[yInd]][1]:
                    matG.add_edge((indFrom[0], yInd), (indTo[0], yInd), seq=G.edge[edge[0]][edge[1]]['seq'])
                    matG.edge[(indFrom[0], yInd)][(indTo[0], yInd)]['penalty'] = (ambigPenaltyFun(edge[0], edge[1]) if G.edge[edge[0]][edge[1]]['seq'] == '-' else 0) + ppmPenaltyFun(edge[1]-edge[0], G.edge[edge[0]][edge[1]]['seq'], precMass)
        else:
            termNodes += [(indFrom[0], indTo[0])]
            matG.node[(indFrom[0], indTo[0])]['seq'] = G.edge[edge[0]][edge[1]]['seq']
            matG.node[(indFrom[0], indTo[0])]['edgePenalty'] = (ambigPenaltyFun(edge[0], edge[1]) if G.edge[edge[0]][edge[1]]['seq'] == '-' else 0) + ppmPenaltyFun(edge[1]-edge[0], G.edge[edge[0]][edge[1]]['seq'], precMass)
    

   
    minNs = [(indHolder[start][0], indHolder[end][0]) for start in starts for end in ends]
    matG.graph['scores'] = {0: {}, 1: {}}
    matG.graph['starts'] = minNs
    matG.graph['ends'] = termNodes
    #print 'minNs: ', minNs
    #print 'maxNs: ', termNodes
    for node in G.nodes():
        indNode = indHolder[node]
        if G.node[node]['memo']:
            if indNode[1] == 1:
                matG.graph['scores'][1][indNode[0]] = max(G.node[node]['memo'].values())
            elif indNode[1] == 0:
                matG.graph['scores'][0][indNode[0]] = max(G.node[node]['memo'].values())
        else:
            matG.graph['scores'][indNode[1]][indNode[0]] = 0
    
    for i in range(symListSize):
        matG.node[(i, i)]['penalty'] += symPenalty
                        
    return matG

def getMatGDimensionRelationship(node1, node2):
    if node1[0] == node2[0]:
        return 1
    elif node1[1] == node2[1]:
        return 0
    else:
        raise ValueError('Error, nodes %s %s differ along more than one dimension' % (str(node1), str(node2)))

def insertLScoresMatG(matG, minNs, maxNs):
    top = nx.topological_sort(matG)
    for minN in minNs:
        matG.node[minN]['lScore'] = 0
      
    maxLScores = {}
    for node in maxNs:
        maxLScores[node] = None
    for pred in top: 
        for node in matG.successors_iter(pred):
            try:
                prevLScore = matG.node[pred]['lScore']
            except KeyError:
                continue
            
            dimDiff = getMatGDimensionRelationship(pred, node)
            newScore = prevLScore+matG.graph['scores'][dimDiff][node[dimDiff]]-matG.node[node]['penalty']-matG.edge[pred][node]['penalty']
            if 'lScore' in matG.node[node]:
                matG.node[node]['lScore'] = max(matG.node[node]['lScore'], newScore)
            else:
                matG.node[node]['lScore'] = newScore
            
            if node in maxLScores:
                maxLScores[node] = matG.node[node]['lScore']-matG.node[node]['edgePenalty']
    
    return max(maxLScores.values())

def insertRScoresMatG(matG, minNs, maxNs):
    revTop = nx.topological_sort(matG)[::-1]
    for maxN in maxNs:
        matG.node[maxN]['rScore'] = matG.graph['scores'][0][maxN[0]]+matG.graph['scores'][1][maxN[1]]-matG.node[maxN]['penalty']-matG.node[maxN]['edgePenalty']
    
    maxRScores = {}
    for node in minNs:
        maxRScores[node] = None
    for succ in revTop:
        for node in matG.predecessors_iter(succ):
            try:
                prevRScore = matG.node[succ]['rScore']
            except KeyError:
                continue
            
            dimDiff = getMatGDimensionRelationship(succ, node)
            newScore = prevRScore+matG.graph['scores'][dimDiff][node[dimDiff]]-matG.node[node]['penalty']-matG.edge[node][succ]['penalty']
            if 'rScore' in matG.node[node]:
                matG.node[node]['rScore'] = max(matG.node[node]['rScore'], newScore)
            else:
                matG.node[node]['rScore'] = newScore
        
            if node in maxRScores:
                maxRScores[node] = matG.node[node]['rScore']
            
    return max(maxRScores.values())

def getSubOptimalSubgraphMatG(matG, G, maxScore, maxNs, starts, ends, alpha=0.9):
    subNodes = {}
    includeEdges = set()
    removeEdges = set()
    if maxScore == None:
        return nx.DiGraph()
    
    cutOff = 1.01 * maxScore if maxScore < 0 else alpha * maxScore
#    cutOff = maxScore - 20
    for edge in matG.edges():
        dimDiff = getMatGDimensionRelationship(edge[0], edge[1])
        try:
#            if (matG.graph['symList'][dimDiff][edge[1][dimDiff]] < 475 and matG.graph['symList'][dimDiff][edge[1][dimDiff]] > 474) and (matG.graph['symList'][dimDiff][edge[0][dimDiff]] < 234 and matG.graph['symList'][dimDiff][edge[0][dimDiff]] > 230):
#                print 'subotpimal path score', (matG.graph['symList'][dimDiff][edge[0][dimDiff]], matG.graph['symList'][dimDiff][edge[1][dimDiff]]), matG.node[edge[0]]['lScore']+matG.node[edge[1]]['rScore']-matG.edge[edge[0]][edge[1]]['penalty']-matG.graph['scores'][1-dimDiff][edge[0][1-dimDiff]], maxScore, cutOff
            if matG.node[edge[0]]['lScore']+matG.node[edge[1]]['rScore']-matG.edge[edge[0]][edge[1]]['penalty']-matG.graph['scores'][1-dimDiff][edge[0][1-dimDiff]] >= cutOff:
                subNodes[matG.graph['symList'][dimDiff][edge[0][dimDiff]]] = 0
                subNodes[matG.graph['symList'][dimDiff][edge[1][dimDiff]]] = 0
#                print matG.node[edge[0]]['lScore']+matG.node[edge[1]]['rScore']-matG.edge[edge[0]][edge[1]]['penalty']-matG.graph['scores'][1-dimDiff][edge[0][1-dimDiff]], maxScore, cutOff
                if matG.edge[edge[0]][edge[1]]['seq'] == '-':
                    includeEdges.add((matG.graph['symList'][dimDiff][edge[dimDiff][dimDiff]], matG.graph['symList'][dimDiff][edge[1-dimDiff][dimDiff]]))
            elif matG.edge[edge[0]][edge[1]]['seq'] == '-':  
                removeEdges.add((matG.graph['symList'][dimDiff][edge[dimDiff][dimDiff]], matG.graph['symList'][dimDiff][edge[1-dimDiff][dimDiff]]))
        except KeyError:
            pass
    
#    print sorted(G.edges())
    for removeEdge in (removeEdges - includeEdges):
        G.remove_edge(*removeEdge)
    
#    print sorted(G.edges())
    for maxN in maxNs:
        try:
            if matG.node[maxN]['lScore']-matG.node[maxN]['edgePenalty'] < cutOff:
#                print 'removing', maxN, matG.node[maxN]['lScore']-matG.node[maxN]['edgePenalty'], cutOff, maxScore
                G.remove_edge(matG.graph['symList'][0][maxN[0]], matG.graph['symList'][1][maxN[1]])
            elif not (matG.graph['symList'][0][maxN[0]] in starts and matG.graph['symList'][1][maxN[1]] in ends):
                subNodes[matG.graph['symList'][0][maxN[0]]] = 0
                subNodes[matG.graph['symList'][1][maxN[1]]] = 0
        except KeyError:
            pass
    
    return nx.subgraph(G, subNodes.keys())

def solveSpectrumGraph(G, deltas, pepMass, specs, unknownPenalty=5, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epsilon=0.02, aas=Constants.aminoacids, verbose=False):
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epsilon=epsilon, aas=aas)
    G = mergeNodes(G, epsilon=2 * epsilon, keepEdges=True, keepAttributes=True)
    
    k = 0
    while (not pathExists(G, min(G.nodes()), max(G.nodes()))):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epsilon=epsilon, aas=aas)
        G = mergeNodes(G, epsilon=2 * epsilon, keepEdges=True, keepAttributes=True)
    
    probNetScoreFunc(G, specs, unknownPenalty=unknownPenalty / 2)
    #intensityScoreFunc(G, lightSpec, heavySpec)
    maxScore = insertLScores(G, 0, pepMass)
    zeroScore = insertRScores(G, 0, pepMass)
    #listEdges(G)
    #listNodes(G)
    print maxScore, zeroScore
    
    subG = getSubOptimalSubgraph(G, maxScore, alpha)
    #print sorted(subG.nodes())
    #print sorted(subG.edges(), key = lambda a: a[0])
    #print 'SubG nodes:', sorted(subG.nodes())
    if verbose:
        print 'Resolving Ambiguous Edges'

    subG = resolveAmbiguousEdges(subG, deltas, specs, subGraphCut=subGraphCut, unknownPenalty=unknownPenalty, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epsilon=epsilon, aas=aas, verbose=verbose)
    
    probNetScoreFunc(subG, specs, useExistingMemo=True, unknownPenalty=unknownPenalty)
    clearLAndRScores(subG)
    destScore = insertLScores(subG, 0, pepMass)
    oriScore = insertRScores(subG, 0, pepMass)
    print oriScore, destScore
    paths = getSubOptimalPaths(subG, alpha, 0, pepMass, oriScore)
    if paths:
        seqs = []
        for path in paths:
            seqs.extend([getSequenceFromNodes(subG, path[1])])
    
    
        scores = list(zip(*paths)[0])
        return scores, seqs
    else:
        return None, None

def getSymPairs(G, PM, epsilon, intEp=3):
    masses = np.sort(np.array(G.nodes()))
    antiMasses = PM + Constants.mods['H2O'] - masses
    halfIndex = np.searchsorted(masses, (PM + Constants.mods['H2O'])/2)
    
    res = epsilon/intEp

    hAntiMasses = np.round(antiMasses/res)
    antiDict = {}
    for i, elem in enumerate(hAntiMasses):
        hInd = int(elem)
        for hMass in range(hInd-intEp, hInd+intEp+1):
            antiDict[hMass] = masses[i]
    
    hMasses = np.round(masses/res)
    pairs = {}
    for i in range(halfIndex):
        if hMasses[i] in antiDict:
            pairs[(masses[i], antiDict[hMasses[i]])] = [True, True]
        else:
            pairs[(masses[i], antiMasses[i])] = [True, False]
    
    for i in range(halfIndex, len(masses)):
        if hMasses[i] in antiDict:
            pairs[(antiDict[hMasses[i]], masses[i])] = [True, True]
        else:
            pairs[(antiMasses[i], masses[i])] = [False, True]
    
    return pairs


    
def clearLAndRScores(G):
    for node in G.nodes():
        try:
            del G.node[node]['lScores']
        except KeyError:
            pass
        
        try:
            del G.node[node]['rScores']
        except KeyError:
            pass

def addSeqToGraph(G, seq, startMass=0, endMass=None):
    nodeGen = PN.nodeInfoGen(seq, startMass=startMass)
    prevNode = startMass
    for node in nodeGen:
        G.add_edge(prevNode, node['prm'], seq=node['formAA'])
        prevNode = node['prm']
    
    if endMass:
        try:
            G.add_edge(prevNode, endMass, seq=node['lattAA'])
        except UnboundLocalError:
            if len(seq) > 0:
                G.add_edge(prevNode, endMass, seq=seq)
    else:
        G.add_edge(prevNode, prevNode + Constants.aminoacids[node['lattAA']][2], seq=node['lattAA'])

def addSeqToGraphBlindMods(G, seq, startMass, endMass):
    nodeGen = Constants.nodeInfoGen(seq[0], startMass=startMass, ambigEdges=seq[1])
    prevNode = startMass

#    print 'Adding seq', seq
    for node in nodeGen:
#        print 'Adding Edge', prevNode, node['prm'], node['formAA']
        G.add_edge(prevNode, node['prm'], seq=node['formAA'])
        prevNode = node['prm']

    try:
        G.add_edge(prevNode, endMass, seq=node['lattAA'])
    except UnboundLocalError:
        if len(seq) > 0:
            G.add_edge(prevNode, endMass, seq=seq)
    

def resolveAmbiguousEdges(G, deltas, specs, PM, hashedAAs, maxEdge=500, minEdge=SGMass, ambigPenaltyFun=lambda fromNode, toNode: 20, ppmPenaltyFun=lambda edgeMass, edgeSeq, precMass: 0, keepAmbigEdges=True, subGraphCut=300, subAlpha=0.3, epStep=0.0005, epMean=0, epSTD=0.02, verbose=False):
    tolerance = 1.5 * epSTD
    if keepAmbigEdges:
        ambigEdges = []
    for edge in G.edges():
        pred = edge[0]
        succ = edge[1]
        if G.edge[pred][succ]['seq'] == '-':
            G.remove_edge(*edge)
            if verbose:
                print 'Ambiguous Edge: ', edge
            if succ - pred < subGraphCut:
                try:
                    subSeqs = hashedAAs[np.round((succ - pred)/epStep)]['seqs']
                except KeyError:
                    if succ - pred > minEdge:
                        subSeqs = Constants.getCandidatePeptides(succ - pred, tolerance)
                    else:
                        subSeqs = None
                    
                if subSeqs:
                    for seq in subSeqs:
                        try:
#                            print seq, epMean, 2*epSTD
                            if np.abs(subSeqs[seq]-epMean) < 2*epSTD:
#                                print 'adding seq', seq
                                addSeqToGraph(G, seq, startMass=pred, endMass=succ)
                        except TypeError:
                            addSeqToGraph(G, seq, startMass=pred, endMass=succ)
                    if verbose:
                        print 'Edge Resolved.', len(subSeqs), 'subsequences.', subSeqs
                elif keepAmbigEdges:
                    if verbose:
                        print 'Edge ambiguous'
                    ambigEdges += [(pred, succ)]
            else:
                #print 'creating subgraph'
                sDeltas = sliceDeltas(deltas, pred, succ)
                #print sDeltas
                subG = getSpectrumGraph(pred, succ, sDeltas, maxNodeNum=40)
                subG = mergeNodes(subG, epsilon=epSTD)
                labelEdgesWithSequence(subG, maxMassCutoff=min(succ - pred - 1, maxEdge), epMean=epMean, epStep=epStep, epSTD=epSTD, minMassCutoff=minEdge, hashedAAs=hashedAAs)
                #print 'Initial subgraph nodes:', sorted(subG.nodes())
                probNetScoreFunc(subG, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, calculateEndPrior=False)
                subMaxScore = insertLScores(subG, [pred], [succ])
                subOriScore = insertRScores(subG, [pred], [succ])
                subG = getSubOptimalSubgraph(subG, subMaxScore, alpha=subAlpha)
                #print 'Suboptimal nodes at L,R scores of %f %f and alpha of %f: ' % (subMaxScore, subOriScore, subAlpha) , str(sorted(subG.nodes()))
                #print 'Suboptimal edges: ', str(sorted(subG.edges()))
                clearLAndRScores(subG)
                if pathExists(subG, [pred], [succ]):
                    subG = resolveAmbiguousEdges(subG, sDeltas, specs, PM, hashedAAs, maxEdge=maxEdge, minEdge=minEdge, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, subGraphCut=subGraphCut, subAlpha=subAlpha, epStep=epStep, epMean=epMean, epSTD=epSTD, keepAmbigEdges=True, verbose=verbose)
                    if pathExists(subG, [pred], [succ]):
                        if verbose:
                            print 'Edge resolved. Subgraph optimization succeeded'
                        G = nx.compose(G, subG)
                    elif keepAmbigEdges:
                        if verbose:
                            print 'Edge ambiguous'
                        ambigEdges += [(pred, succ)]
                elif keepAmbigEdges:
                    if verbose:
                        print 'Edge ambiguous'
                    ambigEdges += [(pred, succ)]

    if keepAmbigEdges:
        for edge in ambigEdges:
            G.add_edge(edge[0], edge[1], seq='-')

    return mergeNodes(G, clusterMethod='single', epsilon=epSTD, keepEdges=True, keepAttributes=True)

def getCorrectHashedModDict(hashedMods, pred, succ, PM):
    
    if pred == 0.0:
#        print pred, succ, PM, 'N-term'
        return hashedMods['N-term']
    elif succ == PM:
#        print pred, succ, PM, 'C-term'
        return hashedMods['C-term']
    else:
#        print pred, succ, PM, 'Anywhere'
        return hashedMods['Anywhere']

def resolveAmbiguousEdgesBlindMods(G, deltas, specs, PM, hashedMods, maxEdge=500, minEdge=SGMass, ambigPenaltyFun=lambda fromNode, toNode: 20, ppmPenaltyFun=lambda edgeMass, edgeSeq, precMass: 0, keepAmbigEdges=True, subGraphCut=300, subAlpha=0.3, epStep=0.0005, epMean=0, epSTD=0.02, verbose=False):
    tolerance = 1.5 * epSTD
    maxNode = max(G.nodes())
#    print 'Entering resolveAmbigousEdgesBlindMods'
    if keepAmbigEdges:
        ambigEdges = []
    for edge in G.edges():
        pred = edge[0]
        succ = edge[1]
        if G.edge[pred][succ]['seq'] == '-':
            print 'Ambiguous edge', edge
            G.remove_edge(*edge)
            hashedAAs = getCorrectHashedModDict(hashedMods, pred, succ, maxNode)                
            try:
                print succ-pred, hashedAAs[np.round((succ-pred)/epStep)]['seqs']
                subSeqs = hashedAAs[np.round((succ - pred)/epStep)]['seqs']
            except KeyError:
                subSeqs = None

            if verbose:
                print 'Ambiguous Edge: ', edge
            if succ - pred < subGraphCut:    
                if subSeqs:
                    for seq in subSeqs:
                        try:
#                            print seq, epMean, 2*epSTD
                            if np.abs(subSeqs[seq]-epMean) < 2*epSTD:
#                                print 'adding seq', seq
                                addSeqToGraphBlindMods(G, seq, startMass=pred, endMass=succ)
                        except TypeError:
                            addSeqToGraphBlindMods(G, seq, startMass=pred, endMass=succ)
                    if verbose:
                        print 'Edge Resolved.', len(subSeqs), 'subsequences.', subSeqs
                elif keepAmbigEdges:
                    if verbose:
                        print 'Edge ambiguous'
                    ambigEdges += [(pred, succ)]
            else:
                addedUnimodEdge = False
                if subSeqs:
                    for seq in subSeqs:
                        if len(seq[0]) == 1 and np.abs(subSeqs[seq]-epMean) < 2*epSTD:
                            G.add_edge(pred, succ, seq='X')
                            addedUnimodEdge = True
                            break
                            
                print 'creating subgraph'
                sDeltas = sliceDeltas(deltas, pred, succ)
                print sDeltas
                subG = getSpectrumGraph(pred, succ, sDeltas, maxNodeNum=40)
                subG = mergeNodes(subG, epsilon=epSTD)
                labelEdgesWithSequenceBlindMods(subG, maxNode, hashedMods, maxMassCutoff=maxEdge, epMean=epMean, epStep=epStep, epSTD=epSTD, minMassCutoff=minEdge)
                print 'Initial subgraph nodes:', sorted(subG.nodes())
                probNetScoreFunc(subG, specs, PM, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, calculateEndPrior=False)
                subMaxScore = insertLScores(subG, [pred], [succ])
                subOriScore = insertRScores(subG, [pred], [succ])
                subG = getSubOptimalSubgraph(subG, subMaxScore, alpha=subAlpha)
                print 'Suboptimal nodes at L,R scores of %f %f and alpha of %f: ' % (subMaxScore, subOriScore, subAlpha) , str(sorted(subG.nodes()))
                print 'Suboptimal edges: ', str(sorted(subG.edges()))
                clearLAndRScores(subG)
                if pathExists(subG, [pred], [succ]):
                    subG = resolveAmbiguousEdgesBlindMods(subG, sDeltas, specs, PM, hashedMods, maxEdge=min(succ - pred - 1, maxEdge), minEdge=minEdge, ambigPenaltyFun=ambigPenaltyFun, ppmPenaltyFun=ppmPenaltyFun, subGraphCut=subGraphCut, subAlpha=subAlpha, epStep=epStep, epMean=epMean, epSTD=epSTD, keepAmbigEdges=True, verbose=verbose)
                    print 'subG nodes after resolveAmbiguousEdgesBlindMods', sorted(subG.nodes())
                    print 'subG edges after resolveAmbiguousEdgesBlindMods', sorted(subG.edges())
                    if pathExists(subG, [pred], [succ]):
                        if verbose:
                            print 'Edge resolved. Subgraph optimization succeeded'
                        G = nx.compose(G, subG)
                    elif keepAmbigEdges and not addedUnimodEdge:
                        if verbose:
                            print 'Edge ambiguous'
                        ambigEdges += [(pred, succ)]
                elif keepAmbigEdges and not addedUnimodEdge:
                    if verbose:
                        print 'Edge ambiguous'
                    ambigEdges += [(pred, succ)]

    if keepAmbigEdges:
        for edge in ambigEdges:
            G.add_edge(edge[0], edge[1], seq='-')

#    print 'Exiting resolveAmbiguousEdgesBlindMods'
    return mergeNodes(G, clusterMethod='single', epsilon=epSTD, keepEdges=True, keepAttributes=True)

def pathExists(G, starts, ends):
    localG = copy.deepcopy(G)
    top = nx.topological_sort(localG)
    for minN in starts:
        try:
            for node in localG.successors(minN):
                if node in ends:
                    return True
                localG.node[node]['pathExists'] = True
        except nx.NetworkXError:
            pass
            
    for pred in top:
        if 'pathExists' not in localG.node[pred]:
            continue
        for node in G.successors_iter(pred):
            if node in ends:
                return True
            localG.node[node]['pathExists'] = True
            
    return False

def printSequences(seqs, scores, maxNum=10):
    Ord = np.argsort(-1 * np.array(scores))
    for i in range(min(Ord.size, maxNum)):
        try:
            print 'Score: ', scores[Ord[i]], 'Seq: ', ''.join(seqs[Ord[i]])
        except TypeError:
            print 'Score: ', scores[Ord[i]], 'Seq: ', seqs[Ord[i]]
        

def printSequencesBlindMods(seqs, scores, maxNum=10):
    Ord = np.argsort(-1 * np.array(scores))
    for i in range(min(Ord.size, maxNum)):
        try:
            print 'Score: ', scores[Ord[i]], 'Seq: ', ''.join(seqs[Ord[i]][0]), 'Mod Names: ', seqs[Ord[i]][1]
        except TypeError:
            print 'Score: ', scores[Ord[i]], 'Seq: ', seqs[Ord[i]][0], 'Mod Names: ', seqs[Ord[i]][1]

def printSequencesThread(seqs, scores, maxNum=10):
    scores_ret = []
    sequence_ret = []

    Ord = np.argsort(-1 * np.array(scores))
    for i in range(min(Ord.size, maxNum)):
        try:
            score = copy.deepcopy(scores[Ord[i]])
            pep = copy.deepcopy(seqs[Ord[i]])
            scores_ret.append(score)
            sequence_ret.append(pep)
        except TypeError:
            score = copy.deepcopy(scores[Ord[i]])
            pep = copy.deepcopy(seqs[Ord[i]])
            scores_ret.append(copy.deepcopy(scores[Ord[i]]))
            sequence_ret.append(copy.deepcopy(seqs[Ord[i]]))

    return scores_ret, sequence_ret

def printSequencesBlindModsThread(seqs, scores, maxNum=10):
    scores_ret = []
    sequence_ret = []
    mods = []

    Ord = np.argsort(-1 * np.array(scores))   
    for i in range(min(Ord.size, maxNum)):
        try:
            score = copy.deepcopy(scores[Ord[i]])
            pep = copy.deepcopy(seqs[Ord[i]][0])
            mod = copy.deepcopy(seqs[Ord[i]][1])
            scores_ret.append(score)
            sequence_ret.append(pep)
            mods.append(mod)
        except TypeError:
            score = copy.deepcopy(scores[Ord[i]])
            pep = copy.deepcopy(seqs[Ord[i]][0])
            mod = copy.deepcopy(seqs[Ord[i]][1])
            scores_ret.append(score)
            sequence_ret.append(pep)
            mods.append(mod)

    return scores_ret, sequence_ret, mods

def updateScoresAndSeqs(G, scores, seqs, subSeqs, specs):
    i = 0
    while i < len(scores):   
        numMods = 0
        for j in range(len(seqs[i])):
            try:
                pairs = subSeqs[seqs[i][j]]
                for subScore, subSeq in pairs:
                    newScore, newSeq = connectSubsequenceToBackBone(G, scores[i], seqs[i], j, subScore, subSeq, specs)
                    scores.append(newScore)
                    seqs.append(newSeq)
                
                if not pairs:
                    newSeq = copy.deepcopy(seqs[i])
                    newSeq[j] = ('None', seqs[i][j])
                    newScore = scores[i] - 5 ** (numMods + 1)
                    scores.append(newScore)
                    seqs.append(newSeq)
                
                del scores[i]
                del seqs[i]
                i -= 1
                break
            except KeyError:
                if type(seqs[i][j]) != str:
                    numMods += 1

        i += 1     

def connectSubsequenceToBackBone(G, score, seq, pos, subScore, subSeq, specs):
    interval = seq[pos]
    #print subSeq, 'subScore: ', subScore
    newScore = score + subScore
    newSeq = copy.deepcopy(seq)
    
    if pos > 0:
        if type(seq[pos - 1]) == str:
            formAA = seq[pos - 1][-1]
        else:
            formAA = '-'    
        newScore -= G.node[interval[0]]['memo'][(formAA, '-')]
        #print 'Subtracting: ', G.node[interval[0]]['memo'][(formAA, '-')], 'AAs: ', (formAA, '-')
        try:
            newScore += G.node[interval[0]]['memo'][(formAA, subSeq[0])]
        except KeyError:
            G.node[interval[0]]['memo'][(formAA, subSeq[0])] = max(spec.getNodeScore(prm=interval[0], formAA=formAA, lattAA=subSeq[0]) for spec in specs)
            newScore += G.node[interval[0]]['memo'][(formAA, subSeq[0])]
        
        #print 'Adding: ', G.node[interval[0]]['memo'][(formAA, subSeq[0])], 'AAs: ', (formAA, subSeq[0])
            
    
    if pos < len(seq) - 1:
        if type(seq[pos + 1]) == str:
            lattAA = seq[pos + 1][0]
        else:
            lattAA = '-'
        newScore -= G.node[interval[1]]['memo'][('-', lattAA)]
        #print 'Subtracting: ', G.node[interval[1]]['memo'][('-', lattAA)], 'AAs: ', ('-', lattAA)
        try:
            newScore += G.node[interval[1]]['memo'][(subSeq[-1], lattAA)]
        except KeyError:
            G.node[interval[1]]['memo'][(subSeq[-1], lattAA)] = max(spec.getNodeScore(prm=interval[1], formAA=subSeq[-1], lattAA=lattAA) for spec in specs)
            newScore += G.node[interval[1]]['memo'][(subSeq[-1], lattAA)]
        
        #print 'Adding: ', G.node[interval[1]]['memo'][(subSeq[-1], lattAA)], 'AAs: ', (subSeq[-1], lattAA)
    
    newSeq[pos] = subSeq
    return newScore, newSeq
    
    
            
def scoreSubsequence(seq, startMass, specs):
    nodeGen = PN.nodeInfoGen(seq, startMass)
    score = 0
    for node in nodeGen:
        nodeScore = max(spec.getNodeScore(**node) for spec in specs)
        score += nodeScore

    return score


def getSubOptimalSubgraph(G, maxScore, alpha=0.9):
    subNodes = np.array([])
    for edge in G.edges():
        pred = edge[0]
        succ = edge[1]
        if isSubOptimalMember(G, edge, alpha, maxScore):
            subNodes = np.union1d(subNodes, np.array([pred, succ]))
        else:
            if G.edge[pred][succ]['seq'] == '-':
                G.remove_edge(pred, succ)
    
    return nx.subgraph(G, subNodes)

def getLRScore(G, edge):
    try:
        lScore = G.node[edge[0]]['lScores'][edge[1]]
    except KeyError:
        lScore = 0
            
    try:
        rScore = G.node[edge[1]]['rScores'][edge[0]]
    except KeyError:
        rScore = 0
    
    #if edge[0] > 300 and edge[1] < 650:
    #    print 'LR SCORE', edge, lScore, rScore, lScore + rScore
    return lScore + rScore      

def isSubOptimalMember(G, edge, alpha, maxScore):
    score = getLRScore(G, edge)
    cutOff = 1.01 * maxScore if maxScore < 0 else alpha * maxScore

    if score and score >= cutOff:
        return True
    else:
        return False
            
    
def findNode(G, mass, epsilon=0.01):
    nodes = np.array(G.nodes())
    nodes.sort()
    ind = np.searchsorted(nodes, mass)
    resid = np.array([nodes[ind] - mass, mass - nodes[ind - 1]])
    minInd = np.argmin(resid)
    if (resid[minInd] < epsilon):
        return nodes[ind - minInd]
    else:
        return None

def listEdges(G):
    edges = G.edges()
    for edge in edges:
        #if G.edge[edge[0]][edge[1]]['seq'] != '-':
        print 'Edge: ', edge, ', AA: ', G.edge[edge[0]][edge[1]]['seq']
        
def listNodes(G):
    nodes = sorted(G.nodes())
    for node in nodes:
        print 'Node: ', node, ', intensity: ', G.node[node]['intensity']
        
#Place-holder scoring function, will replace with probability network
def intensityScoreFunc(G, spec1, spec2):
    nodes = G.nodes()
    for n in nodes:
        print n, G.successors(n)
        scores = {}
        G.node[n]['scores'] = scores
        succ = G.successors(n)
        pred = G.predecessors(n)
        if succ and pred:
            for p in pred:
                for s in succ:
                    scores[(p, s)] = G.node[n]['intensity']

def probNetScoreFunc(G, specs, PM, useExistingMemo=False, ppmPenaltyFun=lambda edgeMass, edgeSeq, precMass: 0, ambigPenaltyFun=lambda fromNode, toNode: 0, calculateEndPrior=True):
    nodes = G.nodes()
    precMass = PM + Constants.mods['H2O'] + Constants.mods['H+']
    for n in nodes:
        if not useExistingMemo:
            G.node[n]['memo'] = {}
        try:
            memo = G.node[n]['memo']
        except KeyError:
            memo = {}
            G.node[n]['memo'] = memo
        
        scores = {}
        G.node[n]['scores'] = scores
        succ = G.successors(n)
        pred = G.predecessors(n)
        if succ and pred:
            for p in pred:
                for s in succ:
                    fAA = G.edge[p][n]['seq']
                    lAA = G.edge[n][s]['seq']
                    try:
                        scores[(p, s)] = memo[(fAA, lAA)]
                    except KeyError:
                        #print 'AAs: prev: %s, latt: %s' % (G.edge[p][n]['seq'], G.edge[n][s]['seq'])
                        #scores[(p, s)] = max(spec.getNodeScore(prm=n, formAA=fAA, lattAA=lAA) for spec in specs)
                        scores[(p, s)] = sum(spec.getNodeScore(prm=n, formAA=fAA, lattAA=lAA) for spec in specs)
                        memo[(fAA, lAA)] = scores[(p, s)]
                    
                    scores[(p,s)] -= ppmPenaltyFun(n-p, fAA, precMass)
#                    print 'ppm penalty', p, n, s, fAA, ppmPenaltyFun(n-p, fAA, precMass)
                    if fAA == '-':
                        scores[(p, s)] -= ambigPenaltyFun(p, n)
                    if not G.successors(s):
                        scores[(p,s)] -= ppmPenaltyFun(s-n, lAA, precMass)
                        if lAA == '-':
                            scores[(p, s)] -= ambigPenaltyFun(n, s)
                    
                    if calculateEndPrior:
                        if not G.predecessors(p):
                            scores[(p, s)] += len(specs)*specs[0].getPriorScore(prm=0, formAA=None, lattAA=fAA)
                        elif not G.successors(s):
                            scores[(p, s)] += len(specs)*specs[0].getPriorScore(prm=PM, formAA=lAA, lattAA=None)
                    
def DFSSubOptimalPaths(G, origins, alpha, maxScore):
    cutOff = alpha * maxScore if maxScore > 0 else 1.01 * maxScore
    queue = deque()
    for origin in origins:
        if origin in G and G.node[origin]['rScores'] > cutOff:
            for n in G.successors(origin):
                for s in G.successors(n):
                    try:
                        if G.node[n]['scores'][(origin, s)] + G.node[s]['rScores'][n] > cutOff:
                            queue.extend([(G.node[n]['scores'][(origin, s)], [origin, n, s])])
                    except KeyError:
                        if G.node[n]['scores'][(origin, s)] > cutOff:
                            queue.extend([(G.node[n]['scores'][(origin, s)], [origin, n, s])])
        
    while queue:
        path = queue.pop()
        yield path
        for s in G.successors_iter(path[1][-1]):     
            n = path[1][-1]
            newScore = path[0] + G.node[n]['scores'][(path[1][-2], s)]
            try:
                if newScore + G.node[s]['rScores'][n] > cutOff:
                    queue.extend([(newScore, path[1] + [s])])
            except KeyError:
                if newScore > cutOff:
                    queue.extend([(newScore, path[1] + [s])])
                
def getSubOptimalPaths(G, alpha, sources, destinations, maxScore):
    paths = []
    gen = DFSSubOptimalPaths(G, sources, alpha, maxScore)
    for item in gen:
        if item[1][-1] in destinations:         
            paths.extend([item])
    
    return paths

def getSources(G):
    sources = []
    for node in G.nodes():
        if not G.predecessors(node) and G.successors(node):
            sources.extend([node])
    
    return sources

def getSinks(G):
    sinks = []
    for node in G.nodes():
        if not G.successors(node) and G.predecessors(node):
            sinks.extend([node])
    
    return sinks 

def insertLScores(G, minNs, maxNs):
    top = nx.topological_sort(G)
    for minN in minNs:
        try:
            for node in G.successors(minN):
                if node in maxNs:
                    continue
                if 'lScores' not in G.node[node]:
                    G.node[node]['lScores'] = {}
                for succ in G.successors(node):
                    try:
                        G.node[node]['lScores'][succ] = max(G.node[node]['scores'][(minN, succ)], G.node[node]['lScores'][succ])
                    except KeyError:
                        G.node[node]['lScores'][succ] = G.node[node]['scores'][(minN, succ)]
                        
        except nx.NetworkXError:
            pass
    
    for pred in top: 
        for node in G.successors_iter(pred):
            try:
                prevLScore = G.node[pred]['lScores'][node]
            except KeyError:
                continue
            
            succ = G.successors(node)
            if succ:
                try: 
                    lScores = G.node[node]['lScores']
                except KeyError:
                    lScores = {}
                    G.node[node]['lScores'] = lScores
                    
                for out in succ:
                    newScore = prevLScore + G.node[node]['scores'][(pred, out)]
                    try:
                        lScores[out] = max(lScores[out], newScore)
                    except KeyError:
                        lScores[out] = newScore
            else:
                try:
                    G.node[node]['lScores'] = max(G.node[node]['lScores'], prevLScore)
                except KeyError:
                    G.node[node]['lScores'] = prevLScore
    
    maxLScores = []
    for node in maxNs:
        try:
            maxLScores.extend([G.node[node]['lScores']])
        except KeyError:
            pass
    
    try:
        return max(maxLScores)
    except ValueError:
        return 0
    
def insertRScores(G, minNs, maxNs):
    revTop = nx.topological_sort(G)[::-1]
    for maxN in maxNs:
        try:
            for node in G.predecessors(maxN):
                if node in minNs:
                    continue
                if 'rScores' not in G.node[node]:
                    G.node[node]['rScores'] = {}
                for pred in G.predecessors(node):
                    try:
                        G.node[node]['rScores'][pred] = max(G.node[node]['rScores'][pred], G.node[node]['scores'][(pred, maxN)])
                    except KeyError:
                        G.node[node]['rScores'][pred] = G.node[node]['scores'][(pred, maxN)]
        except nx.NetworkXError:
            pass

    for succ in revTop:
        for node in G.predecessors_iter(succ):
            try:
                prevRScore = G.node[succ]['rScores'][node]
            except KeyError:
                continue
            
            pred = G.predecessors(node)
            if pred:
                try: 
                    rScores = G.node[node]['rScores']
                except KeyError:
                    rScores = {}
                    G.node[node]['rScores'] = rScores
                    
                for inNode in pred:
                    newScore = prevRScore + G.node[node]['scores'][(inNode, succ)]
                    try:
                        rScores[inNode] = max(rScores[inNode], newScore)
                    except KeyError:
                        rScores[inNode] = newScore
            else:
                try:
                    G.node[node]['rScores'] = max(G.node[node]['rScores'], prevRScore)
                except KeyError:
                    G.node[node]['rScores'] = prevRScore
    
    maxRScores = []
    for node in minNs:
        try:
            maxRScores.extend([G.node[node]['rScores']])
        except KeyError:
            pass
    
    try:
        return max(maxRScores)
    except ValueError:
        return 0

def optimizeSequence(G, seq, path, deltas, epsilon):
    posSeq = []
    if path[0] != 0:
        print 'connecting zero mass with maximally connection graph'
        subG = getSpectrumGraph(0, path[0], deltas)
        mergeNodes(subG)
        labelEdgesWithSequence(subG, epsilon=epsilon, massCutoff=300)
        subPaths = findAllPaths(subG, min(subG.nodes()), max(subG.nodes()))
        posSeq.extend([[getSequenceFromEdges(subG, p) for p in subPaths]])
        
    for i in range(len(seq)):
        if seq[i] == '-':
            minMass = path[i]
            maxMass = path[i + 1]
            subG = getSpectrumGraph(minMass, maxMass, deltas)
            mergeNodes(subG)
            labelEdgesWithSequence(subG, epsilon=epsilon, massCutoff=300)
            subPaths = findAllPaths(subG, min(subG.nodes()), max(subG.nodes()))
            if subPaths:
                posSeq.extend([[getSequenceFromEdges(subG, p) for p in subPaths]])
            else:
                posSeq.extend([seq[i]])
        else:
            posSeq.extend([seq[i]])
                        
    if (path[-1] < max(G.nodes())):
        print 'connecting peptide mass with maximally connected graph'
        subG = getSpectrumGraph(path[-1], max(G.nodes()), deltas)
        mergeNodes(subG)
        labelEdgesWithSequence(subG, epsilon=epsilon, massCutoff=300)
        subPaths = findAllPaths(subG, min(subG.nodes()), max(subG.nodes()))
        posSeq.extend([[getSequenceFromEdges(subG, p) for p in subPaths]])
        
    return posSeq

def sliceDeltas(deltas, minMass, maxMass, offset=Constants.aminoacids['G'][2] - 1):
    minInd = np.searchsorted(deltas[:, 0], minMass + offset)
    maxInd = np.searchsorted(deltas[:, 0], maxMass - offset)
    return deltas[minInd:maxInd]

def getSpectrumGraph(minMass, maxMass, deltas, maxNodeNum=10, addEndMasses=True):
    subG = nx.DiGraph()
    
    if maxNodeNum != None and deltas.shape[0] > maxNodeNum:
        intOrd = np.argsort(deltas[:, 1])
        addNodes = deltas[intOrd[:maxNodeNum]]
    else:
        addNodes = deltas
    
    if addEndMasses:
        subG.add_node(minMass)
        subG.add_node(maxMass)
    
    addMassesToSpectrumGraph(subG, addNodes, maxMass, minMass=minMass)
    return subG

def isCompleteSequence(seq, mass, epsilon=0.02):
    partMass = 0
    for aa in seq:
        if len(aa) > 1:
            return False
        else:
            try:
                partMass += Constants.aminoacids[aa[0]][2]
            except KeyError:
                return False
    
    return (np.abs(partMass - mass) < epsilon)

def prepareDeltas(NInd, CInd, lightPairs, heavyPairs, pepMass, Nmod, Cmod):
    nonNCInd = np.setdiff1d(np.arange(lightPairs.shape[0]), np.append(CInd[0], NInd[0]), assume_unique=True)
    lUnknownMasses = lightPairs[nonNCInd]
    lDeltas = getDeltaIons(lUnknownMasses, pepMass)
    
    nonNCInd = np.setdiff1d(np.arange(heavyPairs.shape[0]), np.append(CInd[1], NInd[1]), assume_unique=True)
    hUnknownMasses = heavyPairs[nonNCInd]
    hDeltas = getDeltaIons(hUnknownMasses, pepMass, nTermMod=Nmod, cTermMod=Cmod)
    
    deltas = np.append(lDeltas, hDeltas, axis=0)
    Ord = np.argsort(deltas[:, 0])
    deltas = deltas[Ord]
    return deltas
    
def getDeltaIons(massIntList, pepMassUnMod, nTermMod=0.0, cTermMod=0.0):
    H = Constants.mods['H+']
    H2O = Constants.mods['H2O']
    NH3 = Constants.mods['NH3']
    numDelt = 2
    addNodesList = np.zeros((massIntList.shape[0] * numDelt, 2))
    for i in range(massIntList.shape[0]):
        mH = massIntList[i, 0] - H
        #Have not yet added in a-ions, or double neutral loss b- and y-ions
        #NMasses = np.array([mH, mH+H2O, mH+NH3]) - nTermMod
        #CMasses = pepMassUnMod - (np.array([mH-H2O, mH, mH-H2O+NH3]) - cTermMod)
        NMasses = np.array([mH]) - nTermMod
        CMasses = pepMassUnMod - (np.array([mH - H2O]) - cTermMod)
        addMasses = np.array(zip(np.append(NMasses, CMasses), [massIntList[i, 1]] * numDelt))
        addNodesList[numDelt * i:numDelt * i + numDelt] = addMasses

    return addNodesList
        
def addMassesToSpectrumGraph(G, massIntList, maxMass, minMass=0, spec1=None, spec2=None, cutOff=5):
    for i in range(massIntList.shape[0]):
        if massIntList[i, 0] > (Constants.aminoacids['G'][2] - 1 + minMass) and massIntList[i, 0] < maxMass + 1 - Constants.aminoacids['G'][2]:
            try:
                score1 = spec1.getNodeScore(prm=massIntList[i, 0], formAA='-', lattAA='-')
                score2 = spec2.getNodeScore(prm=massIntList[i, 0], formAA='-', lattAA='-')
                print massIntList[i, 0], score1, score2
                if max(score1, score2) > cutOff:
                    G.add_node(massIntList[i, 0], intensity=massIntList[i, 1])
            except AttributeError:
                G.add_node(massIntList[i, 0], intensity=massIntList[i, 1]) 

def DFSAllowRepeats(G, origin, AAEdgesOnly=False):
    queue = deque()
    for succ in G.successors_iter(origin):
        if (not AAEdgesOnly) or G.edge[origin][succ]['seq'] != '-':
            queue.extend([(origin, succ)])
    while queue:
        node = queue.pop()
        yield node
        for successor in G.successors_iter(node[1]): 
            if (not AAEdgesOnly) or G.edge[origin][succ]['seq'] != '-':
                queue.extend([(node[1], successor)])

def findAllPaths(G, origin, destination):
    print origin, destination
    dfs = DFSAllowRepeats(G, origin)
    visited = deque()
    allPaths = []
    
    for edge in dfs:
        while visited and (visited[-1][1] != edge[0]):
            visited.pop()
            
        visited.extend([edge])
                
        if visited[-1][1] == destination:
            allPaths.extend([copy.copy(visited)])
                
    return allPaths       
    
def longestPaths(G):
    top = np.array(nx.topological_sort(G))
    lengthTo = np.zeros(top.size)
    pathTo = [[]] * top.size
    
    for vInd in range(top.size):
        v = top[vInd]
        if not pathTo[vInd]:
            pathTo[vInd] = [top[vInd]]
        for w in G.successors_iter(v):
            wInd = (np.where(top == w))[0]
            if lengthTo[wInd] <= lengthTo[vInd] + G.node[w]['intensity']:
                lengthTo[wInd] = lengthTo[vInd] + G.node[w]['intensity']
                pathTo[wInd] = pathTo[vInd] + [w]
    
    maxInds = (np.where(lengthTo == lengthTo.max()))[0]
    print "max lengths: ", [lengthTo[Ind] for Ind in maxInds]
    return [pathTo[Ind] for Ind in maxInds]
    
def getSequenceFromEdges(G, path):
    seq = []
    for edge in path: 
        AA = G.edge[edge[0]][edge[1]]['seq']
        if AA == '-':
            newSeq = Constants.getCandidatePeptides(edge[1] - edge[0], tolerance=0.02)
            if newSeq:
                seq = seq + [newSeq]
            else:
                seq = []
                break
        else:
            seq = seq + [AA]
    
    return seq
                
    
def getSequenceFromNodes(G, path, PM, termModHash, epsilon=0.02):
    seq = []
    for i in range(1, len(path)):
        AA = G.edge[path[i - 1]][path[i]]['seq']
        if AA == '-':
            seq += [(path[i - 1], path[i])]
        else:
            seq = seq + [AA]
        
    try:
        symb = termModHash['NTerm'][np.round(path[0] / epsilon)]
        try:
            seq[0] += symb[0]
        except TypeError:
            seq.insert(1, symb[0])
    except KeyError:
        pass

    try:
        symb = termModHash['CTerm'][np.round((PM - path[-1]) / epsilon)]
        try:
            seq[-1] += symb[0]
        except TypeError:
            seq.insert(len(seq), symb[0])
    except KeyError:
        pass
    
    return seq

def getModDescription(hashedAAs, modAAsDict, edgeMass, epStep=0.0005):
    modAA = ''
    minEp = 10
    hMass = np.round(edgeMass/epStep)

    for seq in hashedAAs[hMass]['seqs']:
        if hashedAAs[hMass]['seqs'][seq] < minEp and len(seq[0]) == 1:
            print 'mod description', seq
            
            minEp = hashedAAs[hMass]['seqs'][seq]
            print np.round(seq[1][0][1], decimals=8), modAAsDict[np.round(seq[1][0][1], decimals=8)]
            modAA = modAAsDict[np.round(seq[1][0][1], decimals=8)]

    return modAA

def getSequenceFromNodesBlindMods(G, path, PM, hashedMods, modAAsDict, termModHash, epStep=0.0005, epsilon=0.02):
    seq = []
    modNames = []
    for i in range(1, len(path)):
        AA = G.edge[path[i - 1]][path[i]]['seq']
        if AA == '-':
            seq += [(path[i - 1], path[i])]
        elif AA == 'X':
            hMass = np.round((path[i] - path[i-1])/epStep)
            hashedAAs = getCorrectHashedModDict(hashedMods, path[i-1], path[i], PM)
            seq = seq + [AA]
            modNames += [getModDescription(hashedAAs, getCorrectHashedModDict(modAAsDict, path[i-1], path[i], PM), path[i] - path[i-1], epStep)]
            
        else:
            seq = seq + [AA]
        
    try:
        symb = termModHash['NTerm'][np.round(path[0] / epsilon)]
        try:
            seq[0] += symb[0]
        except TypeError:
            seq.insert(1, symb[0])
    except KeyError:
        pass

    try:
        symb = termModHash['CTerm'][np.round((PM - path[-1]) / epsilon)]
        try:
            seq[-1] += symb[0]
        except TypeError:
            seq.insert(len(seq), symb[0])
    except KeyError:
        pass
    
    return seq, modNames

#complete is default, and more conservative clustering method
def mergeNodes(G, epsilon=0.02, clusterMethod='complete', keepEdges=False, keepAttributes=True):
    #print 'Merging Nodes'
    massNodes = np.array(G.nodes())
    massNodes.sort()
    if clusterMethod == 'complete':
        clusters = getClustersCompleteLinkage(massNodes, epsilon)
    else:
        clusters = getClustersSingleLinkage(massNodes, epsilon)
        
    massMap = {}
    #print 'Clusters: ', clusters
    for cluster in clusters:
        newMass = sum(cluster) / len(cluster)
        for mass in cluster:
            massMap[mass] = newMass
    
    #print 'Massmap: ',[(key, massMap[key]) for key in sorted(massMap.keys())]
    newG = nx.DiGraph()
    for mass in massMap.keys():
        newG.add_node(massMap[mass])      
        if keepAttributes:
            attributes = G.node[mass]
            for attr in attributes.keys():
                try:
                    #print attr, newG.node[massMap[mass]][attr], G.node[mass][attr], mass, massMap[mass]
                    newG.node[massMap[mass]][attr] += G.node[mass][attr]
                except KeyError:
                    newG.node[massMap[mass]][attr] = G.node[mass][attr]
                except TypeError:
                    try:
                        newG.node[massMap[mass]][attr].update(G.node[mass][attr])
                    except AttributeError:
                        continue
                    except TypeError:
                        continue
                    
        
        if keepEdges:
            for pred in G.predecessors(mass):
                newG.add_edge(massMap[pred], massMap[mass], seq=G.edge[pred][mass]['seq'])
            for succ in G.successors(mass):
                newG.add_edge(massMap[mass], massMap[succ], seq=G.edge[mass][succ]['seq']) 
    
    #print 'Before:', len(G.nodes()), 'nodes'
    #print 'After:', len(newG.nodes()), 'nodes'  
    return newG

def getClustersSingleLinkage(nodes, epsilon=0.02):
    clusters = []
    for mass in nodes:
        clusters.extend([[mass]])
    
    i = 0
    while i < len(clusters):
        mergeClusters = False
        for node1 in clusters[i]:
            for node2 in clusters[i - 1]:
                if (np.abs(node1 - node2) < epsilon):
                    mergeClusters = True
                    break
                if mergeClusters:
                    break
            
        if mergeClusters:
            clusters[i - 1].extend(clusters[i])
            del clusters[i]
        else:
            i = i + 1
           
    return clusters

#greedy complete clustering, does not attempt to minimize spread of each completely linked heirarchical cluster
def getClustersCompleteLinkage(nodes, epsilon=0.02):
    clusters = []
    for mass in nodes:
        clusters.extend([[mass]])
    
    i = 0
    while i < len(clusters):
        mergeClusters = True
        for node1 in clusters[i]:
            for node2 in clusters[i - 1]:
                if (np.abs(node1 - node2) > epsilon):
                    mergeClusters = False
                    break
                if not mergeClusters:
                    break
            
        if mergeClusters:
            clusters[i - 1].extend(clusters[i])
            del clusters[i]
        else:
            i = i + 1
           
    return clusters

SGMass = Constants.aminoacids['G'][2] + Constants.aminoacids['S'][2]

def labelEdgesWithSequence(G, hashedAAs, epMean=0, epSTD=0.02, epStep=0.0005, maxMassCutoff=200, minMassCutoff=SGMass):
    massNodes = np.array(G.nodes())
    massNodes.sort()
    for i in range(massNodes.size):
        for j in range(i):
            massDiff = massNodes[i] - massNodes[j]
            if (massDiff < maxMassCutoff):
                testKey = np.round(massDiff/epStep)
                try: 
                    testEp = hashedAAs[testKey]['min']
                    possSeqs = hashedAAs[testKey]['seqs'].keys()
                except KeyError:
                    if massDiff > minMassCutoff:
                        testEp, possSeqs = epMean, []
                    else:
                        testEp, possSeqs = None, None
                
#               if (massNodes[j] < 716 and massNodes[j]) > 714 or (massNodes[j] < 552 and massNodes[j] > 548):
#                    print 'maybe label', massNodes[j], massNodes[i], testEp, possSeqs, epMean, epSTD
                
                if testEp != None and np.abs(testEp-epMean) < 2*epSTD:
#                    if (massNodes[j] < 716 and massNodes[j]) > 714 or (massNodes[j] < 552 and massNodes[j] > 548):
#                        print 'will label', massNodes[j], massNodes[i], testEp, possSeqs
                    if len(possSeqs) == 1 and possSeqs[0] in Constants.aminoacids:
                        G.add_edge(massNodes[j], massNodes[i], seq=possSeqs[0])
                    else:
                        G.add_edge(massNodes[j], massNodes[i], seq='-')


def labelEdgesWithSequenceBlindMods(G, PM, hashedPeps, epMean=0, epSTD=0.02, epStep=0.0005, maxMassCutoff=200, minMassCutoff=SGMass):
    massNodes = np.array(G.nodes())
    massNodes.sort()
    
    for i in range(massNodes.size):
        for j in range(i):
            massDiff = massNodes[i] - massNodes[j]
            if (massDiff < maxMassCutoff):
                hashedAAs = getCorrectHashedModDict(hashedPeps, massNodes[j], massNodes[i], PM)
                testKey = np.round(massDiff/epStep)
                try: 
                    testEp = hashedAAs[testKey]['min']
                    possSeqs = hashedAAs[testKey]['seqs'].keys()
                except KeyError:
                    if massDiff > minMassCutoff:
                        testEp, possSeqs = epMean, []
                    else:
                        testEp, possSeqs = None, None
                
                if testEp != None and np.abs(testEp-epMean) < 2*epSTD:
                    if len(possSeqs) == 1 and len(possSeqs[0][0]) == 1:
                        G.add_edge(massNodes[j], massNodes[i], seq=possSeqs[0][0])
                    else:
                        G.add_edge(massNodes[j], massNodes[i], seq='-')


def getCompareMatrix(lightPairs, heavyPairs, delta, epsilon=0.01):
    compMatrix = np.zeros((lightPairs.shape[0], heavyPairs.shape[0]))
    print (lightPairs.shape[0], heavyPairs.shape[0])
    for i in xrange(lightPairs.shape[0]):
        for j in xrange(heavyPairs.shape[0]):
            compMatrix[i][j] = (np.abs(heavyPairs[j][0] - lightPairs[i][0] - delta) < epsilon)
    return compMatrix
"""
def getCTermIons(lightPairs, heavyPairs):
    compMatrix = getCompareMatrix(lightPairs, heavyPairs, delta=Constants.mods['*'])
    Ind = np.where(compMatrix == True)
    return Ind

def getNTermIons(lightPairs, heavyPairs):
    compMatrix = getCompareMatrix(lightPairs, heavyPairs, delta=0.0)
    Ind = np.where(compMatrix == True)
    return Ind

"""
def uniqueMasses(masslist, decimel=1):
    a, indices = np.unique(masslist.round(decimel), return_index=True)
    return indices
    

def pairYToBIons(YIons, lightTotMass, heavyTotMass, lightMasses, heavyMasses):
    lightYIons = YIons[:, 0]
    heavyYIons = YIons[:, 1]
    
    def getBYIonPairs(YIons, totalMass, masses, epsilon=0.1):
        YIonsMat = np.tile(YIons, (np.size(masses), 1))
        YIonsMat = np.abs(YIonsMat + masses.reshape((-1, 1)) - totalMass - Constants.mods['H+'])
        IonPairsInd = np.where(YIonsMat < epsilon)
        return zip(masses[IonPairsInd[0]], YIons[IonPairsInd[1]])

    lightBYPairs = getBYIonPairs(lightYIons, lightTotMass, lightMasses)
    heavyBYPairs = getBYIonPairs(heavyYIons, heavyTotMass, heavyMasses)
    
    return (lightBYPairs, heavyBYPairs)
     
def getAASequence(massseq, tolerance=0.02):
    aa = (Constants.getCandidatePeptides(massseq[0] - Constants.mods['Hyd'], tolerance))[0]
    aaseq = [aa]
    for i in xrange(massseq.size - 1):
        aamass = massseq[i + 1] - massseq[i]
        aa = (Constants.getCandidatePeptides(aamass, tolerance))[0]
        aaseq = aaseq + [aa]
    
    return ''.join(aaseq)        

def getPairedSpectraInfoForSequencing(lightSpecs, heavySpecs, lightPrecMass, Nmod, Cmod, epsilon=0.02, verbose=False):
    mergedHeavySpec = SA.mergeSpectra(heavySpecs, epsilon=epsilon)
    mergedLightSpec = SA.mergeSpectra(lightSpecs, epsilon=epsilon)
    
    NTermTable, CTermTable = SA.getNandCIons(mergedLightSpec, mergedHeavySpec, Nmod=Nmod, Cmod=Cmod, epsilon=epsilon)
    NCrossTable, CCrossTable = SA.getCrossPairedIons(mergedLightSpec, mergedHeavySpec, lightPrecMass, Nmod=Nmod, Cmod=Cmod, epsilon=epsilon)

    NTermIonDict = SA.prepIonTableForAddition(NTermTable, ['b', 'b'])
    CTermIonDict = SA.prepIonTableForAddition(CTermTable, ['y', 'y'])
    NCrossIonDict = SA.prepIonTableForAddition(NCrossTable, ['y', 'b'])
    CCrossIonDict = SA.prepIonTableForAddition(CCrossTable, ['b', 'y'])

    allPairedIonsDict = SA.addDicts(SA.reverseDict(SA.addDicts(NTermIonDict, CCrossIonDict)), SA.reverseDict(SA.addDicts(NCrossIonDict, CTermIonDict)))
    
    return {'lightPairs': mergedLightSpec, 'heavyPairs': mergedHeavySpec, 'pairedIonsDict': allPairedIonsDict}
    
def initializeSpectrumGraph(PNet, paramsDict, lightPath, heavyPath=None, ppm=5, usePaired=False, pairConfigName=None, verbose=False):
    if usePaired == True:
        pairConfig = paramsDict['Pair Configurations'][pairConfigName]
        addEnds = getSpectrumGraphEndpointInitFunction(pairConfig['NStatic'], pairConfig['CStatic'], paramsDict['Enzyme']['specificity'])
        termModHash = Constants.getTermModHashForPairConfig(pairConfig)
        sharedInfo, starts, ends, deltas, G = prepPairedSpectrumGraph(lightPath, heavyPath, addEnds, ppm=ppm, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], verbose=verbose)
        precMass = sharedInfo['lightPrecMass']
        
        epsilon = ppm * precMass * 10 ** -6
        lightSpec = PN.Spectrum(PNet, precMass, Nmod=0.0, Cmod=0.0, epsilon=epsilon, spectrum=sharedInfo['lightPairs'])
        lightSpec.initializeNoiseModel()
        heavySpec = PN.Spectrum(PNet, precMass + pairConfig['NMod'] + pairConfig['CMod'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=epsilon, spectrum=sharedInfo['heavyPairs'])
        heavySpec.initializeNoiseModel()
        specs = [lightSpec, heavySpec]
    else:
        addEnds = getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
        termModHash = Constants.createTermModHashAAs()   
        sharedInfo, starts, ends, deltas, G = prepSingleSpectrumGraph(lightPath, addEnds, ppm, 0, 0, verbose=verbose)
        precMass = sharedInfo['precMass']
        
        epsilon = ppm * precMass * 10 ** -6
        specs = [PN.Spectrum(PNet, precMass, Nmod=0, Cmod=0, epsilon=epsilon, spectrum=sharedInfo['pairs'])]
        specs[0].initializeNoiseModel()
    
    return sharedInfo, starts, ends, deltas, termModHash, specs, G

def getSpectrumGraphPaths(G, deltas, specs, starts, ends, PM, termModHash={'NTerm': {}, 'CTerm': {}}, unknownPenalty=5, maxEdge=500, minEdge=SGMass, subGraphCut=300, subAlpha=0.3, alpha=0.95, epsilon=0.02, aas=Constants.aminoacids, maxNumPaths=10, verbose=False):
    hashedAAs = Constants.hashAAs(aas, epsilon)
    labelEdgesWithSequence(G, maxMassCutoff=maxEdge, minMassCutoff=minEdge, epsilon=epsilon, hashedAAs=hashedAAs)
    
    seqEdges = 0
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['seq'] != '-':
            seqEdges += 1
      
    k = 0
    while (not any([pathExists(G, start, end) for start in starts for end in ends])):
        k += 1
        nodes = sorted(G.nodes())
        newDiff = max(nodes[i + k] - nodes[i] for i in range(len(nodes) - k))
        if verbose:
            print 'No Path Found, expanding maxMassCutoff with cutoff', newDiff
        labelEdgesWithSequence(G, maxMassCutoff=np.ceil(newDiff), minMassCutoff=minEdge, epsilon=epsilon, hashedAAs=hashedAAs)
    
    ambigPenaltyFun = getAmbigEdgePenaltyFunction(minEdge, unknownPenalty)
    probNetScoreFunc(G, specs, PM, ambigPenaltyFun=ambigPenaltyFun)
    
    pairs = getSymPairs(G, PM, epsilon)
    #print [pair for pair in sorted(pairs.keys())]
    matG = getMatrixSpectrumGraph(G, PM, starts, ends, epsilon=epsilon, symPenalty=30, ambigPenaltyFun=ambigPenaltyFun)
    #print 'Done getting matG'
    maxScore = insertLScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    zeroScore = insertRScoresMatG(matG, matG.graph['starts'], matG.graph['ends'])
    #print matGmaxScore, matGzeroScore
    subG = getSubOptimalSubgraphMatG(matG, G, maxScore, matG.graph['ends'], alpha=0.90)
    #print 'Suboptimal subgraph nodes from matG: ', str(sorted(subG.nodes()))
    #print 'Suboptimal subgraph edges from matG: ', str(sorted(subG.edges()))
    
    #intensityScoreFunc(G, lightSpec, heavySpec)
    #maxScore = insertLScores(G, starts, ends)
    #zeroScore = insertRScores(G, starts, ends)
    #listEdges(G)
    #listNodes(G)
    #print maxScore, zeroScore
    #subG = getSubOptimalSubgraph(G, maxScore, alpha)
    if verbose:
        print '\nInitial spectrum graph nodes:', sorted(G.nodes()) 
        print '\nSuboptimal Subgraph nodes at alpha factor %f using lScore %f and rScore %f: ' % (alpha, maxScore, zeroScore) + str(sorted(subG.nodes()))
        print '\nSuboptimal Subgraph edges: ', str(sorted(subG.edges()))
        print '\nResolving Ambiguous Edges'
    
    subG = resolveAmbiguousEdges(subG, deltas, specs, PM, hashedAAs, subGraphCut=subGraphCut, ambigPenaltyFun=ambigPenaltyFun, maxEdge=maxEdge, minEdge=minEdge, subAlpha=subAlpha, epsilon=epsilon, verbose=verbose)
    #print 'After resolving and merging:', sorted(subG.nodes())
    #print 'Edges:', sorted(subG.edges())
    probNetScoreFunc(subG, specs, PM, useExistingMemo=True, ambigPenaltyFun=ambigPenaltyFun)
    clearLAndRScores(subG)

    oriScore = insertLScores(subG, starts, ends)
    destScore = insertRScores(subG, starts, ends)

    prevNumPaths = 0
    pathAlphas = [0.90, 0.95, 0.99]
    paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    while (prevNumPaths != len(paths) and len(paths) < maxNumPaths and pathAlphas):
        prevNumPaths = len(paths)
        paths = getSubOptimalPaths(subG, pathAlphas.pop(), starts, ends, oriScore)
    
    return paths, subG

if __name__ == '__main__':
    minEdge = 300
    unknownPenalty = 4
    maxPPMPenalty=12
    seqCutoff = 10000
    alpha = 0.9
    ppmSTD = 5
    ppmSysError = 0
    epStep = 0.00025
    maxEp = 0.1
    verbose = True
    
    paramsDict = DataFile.parseParams('./Misc/LADS_LysC_guanStatic_lDimethDiff_hDimethDiff_cluster_light_heavy.ini')
    pairConfigName = 'lightdimethyl_heavydimethyl'
    aas = {}
    #aas = Constants.getPeptsOfMaxLength(2)
    aas.update(Constants.addPepsToAADict(300))
    hashedAAs = Constants.hashAAsEpsilonRange(aas, epStep, maxEp)
    ambigPenaltyFun = getAmbigEdgePenaltyFunction(minEdge, 0, unknownPenalty)
    ppmPenaltyFun = getPPMPenaltyFun(ppmSTD, hashedAAs, minEdge, maxPPMPenalty, ppmSysError, epStep)
    
    PNet = PN.ProbNetwork('./Scoring_Functions/LysC_HCD_likelihood_b_y_a_pos_AAClass_singleNeutralLoss_IntBin5_config.txt', './Scoring_Functions/ath013833_ath013837_LysC_guanDiMeth_LogTIC_intBin5_bya_singNeutLoss_AAClass_pos_likelihood.model')
    
    t1 = time.time()
    dirPath = '/lab/core_data/adevabhaktuni/ath019105_EllenPoopLADSFracs_2921/_deiso/'
    lightPath = dirPath + '2921.9456.8463.5.dta'
    heavyPath = dirPath + '2921.9415.8425.4.dta'
    usePaired = True

    #Constants.massLadder('PFSNSHNTQK')
    
    if usePaired == True:
        pairConfig = paramsDict['Pair Configurations'][pairConfigName]
        addEnds = getSpectrumGraphEndpointInitFunction(pairConfig['NStatic'], pairConfig['CStatic'], paramsDict['Enzyme']['specificity'])
        termModHash = Constants.getTermModHashForPairConfig(pairConfig)

        lightSpecs = [DataFile.getMassIntPairs(lightPath)]
        heavySpecs = [DataFile.getMassIntPairs(heavyPath)]

        pMass = DataFile.getPrecMassAndCharge(lightPath)[0]
        
        epMean = ppmSysError * pMass * 10**-6
        epSTD = ppmSTD * pMass * 10**-6

        sharedInfo, starts, ends, deltas, G = prepPairedSpectrumGraph(lightSpecs, heavySpecs, pMass, addEnds, ppmSTD=ppmSTD, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], verbose=verbose)

        print 'precMass:', pMass, 'starts and ends:', starts, ends

        lightSpec = PN.Spectrum(PNet, pMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=lightSpecs[0])
        lightSpec.initializeNoiseModel()
        heavySpec = PN.Spectrum(PNet, pMass + pairConfig['NMod'] + pairConfig['CMod'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD, spectrum=heavySpecs[0])
        heavySpec.initializeNoiseModel()
        specs = [lightSpec, heavySpec]
    else:
        addEnds = getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
        termModHash = Constants.createTermModHashAAs(N=Constants.NTermMods, C=Constants.CTermMods)   
        s1 = time.time()
        sharedInfo, starts, ends, deltas, G = prepSingleSpectrumGraph(lightPath, addEnds, ppmSTD, 0, 0, verbose=verbose)
        pMass = sharedInfo['precMass']
        
        print 'precMass:', pMass, 'starts and ends:', starts, ends
        epMean = ppmSysError * pMass * 10**-6
        epSTD = ppmSTD * pMass * 10**-6
        specs = [PN.Spectrum(PNet, pMass, Nmod=0, Cmod=0, epsilon=2*epSTD, spectrum=sharedInfo['pairs'])]
        specs[0].initializeNoiseModel()
    
    print getSpectrumGraphData(G, deltas, specs, starts, ends, pMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, epMean=epMean, epSTD=epSTD, termModHash=termModHash, maxEdge=seqCutoff, minEdge=minEdge, subGraphCut=300, subAlpha=0.3, alpha=alpha, epStep=epStep, verbose=verbose)

    t2 = time.time()
    print 'Time taken: ', t2 - t1
    """
    Constants.massLadder('VFJENVJRDAVTYTEHAK')
    start = 359.220880
    end = 1356.740120
    pairs = DataFile.getMassIntPairs(lightPath)
    precMass, charge = DataFile.getPrecMassAndCharge(lightPath)
    spec = PN.Spectrum(PNet, precMass, Nmod=0.0, Cmod=0.0, spectrum=pairs)
    spec.initializeNoiseModel()
    deltas = getDeltaIons(pairs, precMass - Constants.mods['H+'] - Constants.mods['H2O'], 0.0, 0.0)
    getSubGraphScore(deltas, start, end, [spec])
    
    print Constants.getCandidatePeptides(end-start, tolerance=0.03)
    """
    
