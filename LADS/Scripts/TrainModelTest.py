'''
Created on Jul 8, 2011

@author: Arun
'''

import os
import sys
#make sure that all imports refer to LADS-rev120 version of LADS
sys.path.insert(1, os.path.abspath(os.pardir))

import numpy as np
import glob
import pickle
import random
import copy
#import scipy.optimize as sciopt

import DataFile
import ProbNetwork as PN
import Constants
import Analytics as An
import ArgLib
import DeNovoSequencer as DNS

def parseScans(fDict, prog, seqMap, dbDict, delimiter=',', seqDelimLen=2):
    processedInfo = {}
    for csvfile in fDict.keys():
        data = DataFile.getScanInfo(csvfile, dbDict[prog]['fields'], delimiter=delimiter)
        processedInfo[csvfile] = An.preprocessDatabaseScanInfo(data, seqMap[fDict[csvfile]], dbDict[prog]['fieldmap'], seqDelimLen=seqDelimLen)
    
    return processedInfo

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = DataFile.getScanNum(dta)
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'sequenced': False}
        
    return scanFDict

def getHyperParameterTrainingSet(scanFDicts, processedInfo, pairConfig=None, numSpectra=100, massRange=None, clusterSize=None):
    trainingSet = []
    for run in processedInfo:
        trueClusters = An.getTrueClusters(processedInfo[run], dbDict['infoMap'], {run: 'SEQUEST'}, run, isComp=False, ppm=options.ppmstd)
        if pairConfig != None:
            truePairs = An.findPairsInTrueClusters(trueClusters, pairConfig, processedInfo[run], dbDict['infoMap'], {run: 'SEQUEST'}, run, isComp=False, ppm=options.ppmstd)
            trainingSet += [(run, scanFDicts[run][list(trueClusters[pair[0]])[0]]['precMass'], (trueClusters[pair[0]], trueClusters[pair[1]])) for pair in truePairs]
        else:
            trainingSet += [(run, scanFDicts[run][list(cluster)[0]]['precMass'], (cluster,)) for cluster in trueClusters]   
    
    if massRange != None:
        trainingSet = [item for item in trainingSet if (item[1] >= massRange[1] and item[1] < massRange[0])]
    if clusterSize != None:
        trainingSet = [item for item in trainingSet if (max([len(cluster) for cluster in item[2]]) >= clusterSize)]
    trainingSet = random.sample(trainingSet, min(numSpectra, len(trainingSet)))
    if len(trainingSet) < numSpectra:
        print 'ERROR: Number of spectra/clusters/pairs available %i is less than number requested for optimization %i' % (len(trainingSet), numSpectra)

    dtaList, trueSeqs =  [], []
    for item in trainingSet:
        run, precMass, clusters = item

        dtaList += [tuple([tuple([scanFDicts[run][scanF]['dta'] for scanF in cluster]) for cluster in clusters])]
        trueSeqs += [processedInfo[run][list(clusters[0])[0]]['Peptide']]

    return dtaList, trueSeqs
    

def getMaxInterval(attr, testVals, initHyperParameters, dtaList, trueSeqs):
    attrArr = np.zeros((len(testVals),2))
    attrArr[:,0] = np.array(testVals)
    hyperParameters = copy.deepcopy(initHyperParameters)
    for i in range(attrArr.shape[0]):
        print 'testing %s: %f' % (attr, float(attrArr[i][0]))
        hyperParameters[attr] = attrArr[i][0]
        acc, prec = getAccAndPrec(dtaList, trueSeqs, hyperParameters)
        print acc, prec
        attrArr[i][1] = acc

    Ord = np.argsort(-1*attrArr[:,1])
    return [attrArr[Ord[0]][0], attrArr[Ord[1]][0]]
    
def accuracyCostFunctionGenerator(scanFDicts, processedInfo, hyperParameters, epsilons, paramsToOptimize=['ambigopen', 'ambigextend', 'ppmpen', 'alpha', 'ppmstd', 'maxnodenum'], numSpectra=100, pairConfig=None, massRange=None, clusterSize=None):
    dtaList, trueSeqs = getHyperParameterTrainingSet(scanFDicts, processedInfo, numSpectra=numSpectra, pairConfig=pairConfig, massRange=massRange, clusterSize=clusterSize)
    defaultParams = copy.copy(hyperParameters)
    
    #minimize 1 - accuracy
    def costFunc(x):
        for i, param in enumerate(paramsToOptimize):
            hyperParameters[param] = defaultParams[param] + x[i] * epsilons[param]
        hyperParameters['maxnodenum'] = int(np.round(hyperParameters['maxnodenum']))
        print hyperParameters
        if pairConfig == None:
            acc, prec = getAccAndPrec(dtaList, trueSeqs, hyperParameters)
        else:
            acc, prec = getAccAndPrecPaired(dtaList, trueSeqs, hyperParameters, pairConfig)
        logFile.write('\t'.join([str(hyperParameters[param]) for param in hyperParameters] +  [str(acc), str(prec)]) + '\n')
        print acc, prec
        return 1 - acc

    return costFunc

def getAccAndPrecPaired(dtaList, trueSeqs, hyperParameters, pairConfig):
    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(hyperParameters['minedge'], hyperParameters['ambigopen'], hyperParameters['ambigextend'])
    ppmPenaltyFun = DNS.getPPMPenaltyFun(hyperParameters['ppmstd'], hashedAAs, hyperParameters['minedge'], hyperParameters['ppmpen'], 0, epStep)

    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
    termModHash = Constants.createTermModHashAAs(N=copy.deepcopy(Constants.NTermMods), C=copy.deepcopy(Constants.CTermMods))
    hyperParameterTrainSeqs = []
    
    for pair in dtaList:
        lightSpecs = [DataFile.getMassIntPairs(dta) for dta in pair[0]]
        heavySpecs = [DataFile.getMassIntPairs(dta) for dta in pair[1]]
        lightPrecMass = np.average([DataFile.getPrecMassAndCharge(dta)[0] for dta in pair[0]])

        epSTD = hyperParameters['ppmstd'] * lightPrecMass * 10**-6
        sharedInfo, starts, ends, deltas, G = DNS.prepPairedSpectrumGraph(lightSpecs, heavySpecs, lightPrecMass, addEnds, ppmSTD=hyperParameters['ppmstd'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'])
        
        specs = []
        for massIntPairs in lightSpecs:
            specs += [PN.Spectrum(trainNet, lightPrecMass, Nmod=0, Cmod=0, epsilon=2*epSTD, spectrum=massIntPairs) for massIntPairs in lightSpecs]
        for massIntPairs in heavySpecs:
            specs += [PN.Spectrum(trainNet, lightPrecMass+pairConfig['NMod']+pairConfig['CMod'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD, spectrum=massIntPairs) for massIntPairs in heavySpecs]
        [spec.initializeNoiseModel() for spec in specs]
        graphData =  DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, lightPrecMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, epMean=0, epSTD=epSTD, termModHash=termModHash, maxEdge=500, minEdge=hyperParameters['minedge'], subGraphCut=hyperParameters['subgraphcut'], subAlpha=hyperParameters['subalpha'], alpha=hyperParameters['alpha'], epStep=epStep)
        
        seq = An.preprocessSequence(graphData['seq'], seqMap['LADS Unit Test'], ambigEdges=graphData['ambiguous edges'])
        hyperParameterTrainSeqs += [(seq, graphData['ambiguous edges'])]

    totPrec, totAcc = 0,0
    for i, trueSeq in enumerate(trueSeqs):
        prec, acc = An.comparePeptideResults(hyperParameterTrainSeqs[i][0], trueSeq, ambigEdges1=hyperParameterTrainSeqs[i][1], ambigEdges2=[])[:2]
        totPrec += prec
        totAcc += acc

    return (totAcc/len(dtaList), totPrec/len(dtaList))

    
def getAccAndPrec(dtaList, trueSeqs, hyperParameters):
    ambigPenaltyFun = DNS.getAmbigEdgePenaltyFunction(hyperParameters['minedge'], hyperParameters['ambigopen'], hyperParameters['ambigextend'])
    ppmPenaltyFun = DNS.getPPMPenaltyFun(hyperParameters['ppmstd'], hashedAAs, hyperParameters['minedge'], hyperParameters['ppmpen'], 0, epStep)
    
    addEnds = DNS.getSpectrumGraphEndpointInitFunction(np.array(Constants.NTermMods.values()), np.array(Constants.CTermMods.values()), paramsDict['Enzyme']['specificity'])
    termModHash = Constants.createTermModHashAAs(N=copy.deepcopy(Constants.NTermMods), C=copy.deepcopy(Constants.CTermMods))
    hyperParameterTrainSeqs = []
    for cluster in dtaList:
        
        specsInfo = [DataFile.getMassIntPairs(dta) for dta in cluster[0]]
        pMass = np.average([DataFile.getPrecMassAndCharge(dta)[0] for scanF in cluster[0]])
#        pMass = DataFile.getPrecMassAndCharge(dta)[0]
#        pairs = DataFile.getMassIntPairs(dta)
        sharedInfo, starts, ends, deltas, G = DNS.prepSingleSpectrumGraph(specsInfo, pMass, addEnds, hyperParameters['ppmstd'], 0, 0, maxNodeNum=hyperParameters['maxnodenum'])
        
        epSTD = hyperParameters['ppmstd'] * pMass * 10**-6
        specs = [PN.Spectrum(trainNet, pMass, Nmod=0, Cmod=0, epsilon=2*epSTD, spectrum=massIntPairs) for massIntPairs in specsInfo]
        [spec.initializeNoiseModel() for spec in specs]
        
        graphData = DNS.getSpectrumGraphData(G, deltas, specs, starts, ends, pMass - Constants.mods['H+'] - Constants.mods['H2O'], ambigPenaltyFun, ppmPenaltyFun, hashedAAs, epMean=0, epSTD=epSTD, termModHash=termModHash, maxEdge=500, minEdge=hyperParameters['minedge'], subGraphCut=hyperParameters['subgraphcut'], subAlpha=hyperParameters['subalpha'], alpha=hyperParameters['alpha'], epStep=epStep)
        
        seq = An.preprocessSequence(graphData['seq'], seqMap['LADS Unit Test'], ambigEdges=graphData['ambiguous edges'])
        hyperParameterTrainSeqs += [(seq, graphData['ambiguous edges'])]

    totPrec, totAcc = 0,0
    for i, trueSeq in enumerate(trueSeqs):
        prec, acc = An.comparePeptideResults(hyperParameterTrainSeqs[i][0], trueSeq, ambigEdges1=hyperParameterTrainSeqs[i][1], ambigEdges2=[])[:2]
        totPrec += prec
        totAcc += acc

    return (totAcc/len(dtaList), totPrec/len(dtaList))
                        

if __name__ == '__main__':
    print 'Note: search program arguments are of the form \"{\'path/to/results\':  \'path/to/dtadir\'}\"\nmainprogname argument refers to the type of SEQUEST/MASCOT (must edit code to specify) search done to get the search results'
    options = ArgLib.parse(['init', 'config', 'output', 'mascot', 'sequest', 'ppmstd', 'database', 'mainprogname', 'symbolmap', 'minedge', 'alpha'], [{'opts': ('-P', '--pairconfig'), 'attrs': {'type': 'string', 'dest': 'pairconfig', 'help': 'Optional argument indicating pair configuration of results for optimization'}}, {'opts': ('-M', '--massrange'), 'attrs': {'type': 'string', 'dest': 'massrange', 'help': 'Optional argument specifying precursor mass range of spectra to optimize over. Format: [low_mass, high_mass]'}}, {'opts': ('-C', '--clustersize'), 'attrs': {'type': 'string', 'dest': 'clustersize', 'help': 'Minimum size of clusters to optimize over'}}, {'opts': ('-a', '--parameters'), 'attrs': {'type': 'string', 'dest': 'parameters', 'help': 'Parameters to optimize over using bfgs'}}])
    
    trainNet = PN.ProbNetwork(options.config)
    dbDict = DataFile.getDBInfo(options.database)
    paramsDict = DataFile.parseParams(options.init)

    clusterSize = eval(options.clustersize) if options.clustersize != None else None
    parametersToOptimize = eval(options.parameters) if options.parameters else ['ambigopen', 'ambigextend', 'ppmpen', 'maxnodenum', 'alpha']
    massRange = eval(options.massrange) if options.massrange != None else None
    pairConfig = paramsDict['Pair Configurations'][options.pairconfig] if options.pairconfig != None else None
    
    progDict = {options.mainprogname: 'MASCOT', 'LADS Unit Test': 'LADS'}

    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap(progDict, symbolMap, paramsDict)

    processedInfo = {}
    scanFDicts = {}
    if options.sequest:
        dtaDict = eval(options.sequest)
        SEQUESTdict = {}
        for progPath in dtaDict:
            SEQUESTdict[progPath] = options.mainprogname
            
            dtaList = glob.glob(dtaDict[progPath] + '/*dta')
            scanFDicts[progPath] = getScanFDict(dtaList)
        print SEQUESTdict
        processedInfo.update(parseScans(SEQUESTdict, 'SEQUEST', seqMap, dbDict))

    if options.mascot:
        dtaDict = eval(options.mascot)
        MASCOTdict = {}
        for progPath in dtaDict:
            MASCOTdict[progPath] = options.mainprogname
            
            dtaList = glob.glob(dtaDict[progPath] + '/*dta')
            scanFDicts[progPath] = getScanFDict(dtaList)
        print MASCOTdict
        processedInfo.update(parseScans(MASCOTdict, 'MASCOT', seqMap, dbDict))

#    print Constants.NTermMods
    """
    processedInfo = {}    
    progDict = {}
    if options.mascot:
        MASCOTdict = eval(options.mascot)
        processedInfo.update(CS.parseDBScans(MASCOTdict, 'MASCOT', seqMap, dbDict))
        CS.updateProgDict(MASCOTdict, progDict, 'MASCOT')
    elif options.sequest:
        SEQUESTdict = eval(options.sequest)
        processedInfo.update(CS.parseDBScans(SEQUESTdict, 'SEQUEST', seqMap, dbDict))
        CS.updateProgDict(SEQUESTdict, progDict, 'SEQUEST')
    """
    numScans = 0
#    hyperParameterTrainSet = []
#    trueSeqsB = []
    for run in processedInfo:
        print run
        for scanF in processedInfo[run]:
            dtaFile = scanFDicts[run][scanF]['dta']
            precMass = DataFile.getPrecMassAndCharge(dtaFile)[0]
            print 'Current Scan Num:', scanF, dtaFile
            massIntPairs = DataFile.getMassIntPairs(dtaFile)
            numScans += 1
#            if random.randint(1, 25) == 1:
#                hyperParameterTrainSet += [dtaFile]
#                trueSeqs += [processedInfo[run][scanF]['Peptide']]
            epSTD = options.ppmstd * precMass * 10 ** -6
            S = PN.Spectrum(trainNet, precMass, 0.0, 0.0, massIntPairs, epsilon=2*epSTD, sequence=processedInfo[run][scanF]['Peptide'])
            S.useToTrain()
    
    model = trainNet.getModelFromCounts(hyperParameters={'ambig_open': 20, 'ambig_extend': 20, 'ppm_pen': 20}, comment='Trained from configuration file %s on %i spectra using datasets %s' % (options.config, numScans, str(dtaDict.values())))
    epStep = 0.00025
    print 'prepping for de novo sequencing'
    aas = Constants.addPepsToAADict(300)
    hashedAAs = Constants.hashAAsEpsilonRange(aas, epStep, 0.1)
    # Initialize hyperParameters to reasonable default values

    hyperParameters = {'minedge': 300, 'ambigopen': 0, 'ambigextend': 4, 'ppmpen': 12, 'ppmstd': 5, 'maxnodenum': 50, 'subgraphcut': 300, 'subalpha': 0.3, 'alpha': 0.9}
    """
    #define epsilons as desired step sizes when  epsilon given to fmin_bfgs is 0.1 (i.e., epsilon['attr'] = 0.1*desired step size)
    epsilons = {'minedge': 200, 'ambigopen': 4, 'ambigextend': 4, 'ppmpen': 4, 'ppmstd': 1, 'maxnodenum': 50, 'subgraphcut': 200, 'subalpha': 0.1, 'alpha': 0.2}
    
    
    logFile = open(options.output + '.log', 'w')
    logFile.write('\t'.join([param for param in hyperParameters] + ['Accuracy', 'Precision']) + '\n')
    
    # gradient calculation is done as percentage deviation from initial value
    x0 = [0] * len(parametersToOptimize)
    
    costFunc = accuracyCostFunctionGenerator(scanFDicts, processedInfo, hyperParameters, epsilons, paramsToOptimize=parametersToOptimize, numSpectra=100, pairConfig=pairConfig, massRange=massRange, clusterSize=clusterSize)
    optParams = sciopt.fmin_bfgs(costFunc, np.array(x0), gtol=0.1, epsilon=0.1, maxiter=4)
    
    logFile.close()
    print optParams
    """
    
    """
    dtaList, trueSeqs = getHyperParameterTrainingSet(scanFDicts, processedInfo, numSpectra=300)
    coarseInterval = getMaxInterval('ambigopen', [0, 10, 20, 30, 40, 50], hyperParameters, dtaList, trueSeqs)
    fineInterval = getMaxInterval('ambigopen', [coarseInterval[0] + step*(coarseInterval[1] - coarseInterval[0])/5 for step in range(6)], hyperParameters, dtaList, trueSeqs)
    hyperParameters['ambigopen'] = fineInterval[0]
    print hyperParameters

    coarseInterval = getMaxInterval('ambigextend', [0, 10, 20, 30, 40, 50], hyperParameters, dtaList, trueSeqs)
    fineInterval = getMaxInterval('ambigextend', [coarseInterval[0] + step*(coarseInterval[1] - coarseInterval[0])/5 for step in range(6)], hyperParameters, dtaList, trueSeqs)
    hyperParameters['ambigextend'] = fineInterval[0]
    print hyperParameters
    
    coarseInterval = getMaxInterval('ppmpen', [0, 10, 20, 30, 40, 50], hyperParameters, dtaList, trueSeqs)
    fineInterval = getMaxInterval('ppmpen', [coarseInterval[0] + step*(coarseInterval[1] - coarseInterval[0])/5 for step in range(6)], hyperParameters, dtaList, trueSeqs)
    hyperParameters['ppmpen'] = fineInterval[0]
    print hyperParameters

    coarseInterval = getMaxInterval('alpha', [0.86, 0.88, 0.90, 0.92, 0.94], hyperParameters, dtaList, trueSeqs)
    hyperParameters['alpha'] = coarseInterval[0]
    
    model['hyper_parameters'] = hyperParameters
    """

    with open(options.output, 'w') as fout:
        pickle.dump(model, fout)
    
            

        
