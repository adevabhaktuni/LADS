import os
import sys

sys.path.insert(1, os.path.abspath('../'))
libpath = os.path.abspath(os.pardir) + "/libsvm-3.12/python"
sys.path.append(libpath)
                
import svmutil

import ArgLib
import Analytics as An
import GenerateSVMFormDiscriminantTrainingData as GLFD
import DataFile
import Discriminator
import SequenceDTAsTDV as SDTDV
import ProbNetwork as PN
import Constants

import numpy as np
import glob
import pickle
import copy

def getScanScoreDictClusterNormScore(LADSSeqInfo, seqEntry):
    scanScoreDict = {}
    clusterLen = len(seqEntry[0]) + len(seqEntry[1])
    PSM = LADSSeqInfo[seqEntry][0]

    PSMData = {'Seq': (PSM[1], PSM[2]), 'Raw Score': PSM[0], 'Post Score': PSM[0]/clusterLen}
    for scan in seqEntry[0]:
        scanScoreDict[scan] = PSMData
    for scan in seqEntry[1]:
        scanScoreDict[scan] = PSMData

    return scanScoreDict

    

def getSpectrumAndPSMFeatureDict(LADSSeqInfo, seqEntry, scanFDict, pairConfig, PNet):

    featureList = []
    lightScans = seqEntry[0]
    heavyScans = seqEntry[1]
    
    lightSpecs = [DataFile.getMassIntPairs(scanFDict[int(lightScanF)]['dta']) for lightScanF in lightScans]
    heavySpecs = [DataFile.getMassIntPairs(scanFDict[int(heavyScanF)]['dta']) for heavyScanF in heavyScans]
    avgLightPrecMass = np.average(np.array([scanFDict[lightScanF]['precMass'] for lightScanF in lightScans]))
    
    epSTD = options.ppmstd * 10**-6 * avgLightPrecMass
    
    specs = []
    for i, massIntPairs in enumerate(lightSpecs):
        specs += [PN.Spectrum(PNet, scanFDict[lightScans[i]]['precMass'], Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=massIntPairs)]
    for i, massIntPairs in enumerate(heavySpecs):
        specs += [PN.Spectrum(PNet, scanFDict[heavyScans[i]]['precMass'], Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD, spectrum=massIntPairs)]
    for spec in specs:
        spec.initializeNoiseModel()
                                                                                                                                                    
    clusterPairingStats = Discriminator.getClusterPairingStats(lightSpecs, heavySpecs, avgLightPrecMass, pairConfig, epSTD=epSTD)
    GLFD.addClusterPairingStatsToFeatureList(clusterPairingStats, featureList)

    scoreStats = {}
    truePMs = {}
    prmLadders = {}
    for PSM in LADSSeqInfo[seqEntry]:
        lightSeq = An.preprocessSequence(PSM[1], seqMap, ambigEdges=PSM[2])
        scoreStats[PSM[:2]] = Discriminator.getScoreStats(specs, lightSeq, ambigEdges=PSM[2])

        prmLadderWithEnds = An.getPRMLadder(lightSeq, ambigEdges=PSM[2], addEnds=True)
        truePMs[PSM[:2]] = prmLadderWithEnds[-1]
        prmLadders[PSM[:2]] = prmLadderWithEnds[1:-1]
        
    PSMList = scoreStats.keys()
    spectrumOrderedScoreStats, clusterScoreStats = GLFD.compileScoreStats(scoreStats, specs, PSMList)

    spectrumAndPSMSpecificFeatureDict = {}
        
    PSMIndexDict = dict([(PSM, i) for i, PSM in enumerate(PSMList)])
    for i, PSM in enumerate(LADSSeqInfo[seqEntry]):
        PSMSpecificFeatureList = copy.copy(featureList)

        peptLength = len(prmLadders[PSM[:2]]) + 1

        # Add LADS PScore (and normalized variants)  and delta rank, delta score (LADS PScore) to feature list
        PSMSpecificFeatureList += [PSM[0], PSM[0]/peptLength, PSM[0]/len(specs), -i, PSM[0]-LADSSeqInfo[seqEntry][0][0]]
        # Add Total Path Score (and normalized variants) and delta rank, delta score (total path score)  and total minimum node score to feature list
        totalPathScore = scoreStats[PSM[:2]]['Total Path Score']
        PSMSpecificFeatureList += [totalPathScore, totalPathScore/peptLength, totalPathScore/len(specs), -clusterScoreStats['PSM Rankings'][PSMIndexDict[PSM[:2]]], totalPathScore-clusterScoreStats['Max Cluster Path Score'], scoreStats[PSM[:2]]['Total Minimum Node Score']]
        
        # Add minimum path score, maximum path score, (and normalized variants) and minimum score/maximum score for cluster to feature list
        PSMSpecificFeatureList += [scoreStats[PSM[:2]]['Minimum Path Score'], scoreStats[PSM[:2]]['Minimum Path Score']/peptLength, scoreStats[PSM[:2]]['Maximum Path Score'], scoreStats[PSM[:2]]['Maximum Path Score']/peptLength, scoreStats[PSM[:2]]['Minimum Path Score']/scoreStats[PSM[:2]]['Maximum Path Score']]
        
        # Add difference between minimum and maximum ranking for PSM across cluster to feature list
        rankingsForPSM = [spectrumOrderedScoreStats[i]['PSM Rankings'][PSMIndexDict[PSM[:2]]] for i in spectrumOrderedScoreStats]
        PSMSpecificFeatureList += [min(rankingsForPSM) - max(rankingsForPSM)]
        
        #Add Number forbidden node pairs (and normalized variants) to feature list
        numForbiddenPairs = Discriminator.getNumForbiddenPairs(prmLadders[PSM[:2]], avgLightPrecMass)
        PSMSpecificFeatureList += [numForbiddenPairs, 2.0*numForbiddenPairs/(peptLength-1)]

        # Add number of ambiguous edges to feature list
        PSMSpecificFeatureList += [len(PSM[2])]
        
        # Add stats for PRM Evidence over cluster (and normalized variants) to feature list
        PSMSpecificFeatureList += [scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['All Evidence'], scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['All Evidence']/float(peptLength-1), scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['Majority Evidence'], scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['Majority Evidence']/float(peptLength-1), scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['None Evidence'], scoreStats[PSM[:2]]['Aggregate PRM Score Statistics']['None Evidence']/float(peptLength-1)]

        # Add stats for paired PRMs and their corresponding ion types to feature list
        pairedPRMStats = Discriminator.getPairedPRMStats(prmLadders[PSM[:2]], clusterPairingStats['Light Merged Spec'], clusterPairingStats['Heavy Merged Spec'], lightSpecs, heavySpecs, clusterPairingStats['Cluster Paired PRM Information'], epSTD=epSTD)
        GLFD.addPairedPRMStatsToFeatureList(pairedPRMStats, PSMSpecificFeatureList, len(prmLadders[PSM[:2]]))

        pairedPRMLadder = pairedPRMStats['Paired PRM Ladder']        
    
        for i, scan in enumerate(lightScans):
            spectrumSpecificFeatureList = copy.copy(PSMSpecificFeatureList)
            # Add path score (and normalized variants), delta rank, delta score, number of negative PRMs, and minimum node score for spectrum to feature list
            pathScore = spectrumOrderedScoreStats[i]['Path Scores'][PSMIndexDict[PSM[:2]]]
            numNegativePRMs = spectrumOrderedScoreStats[i]['Num Negative PRMs'][PSMIndexDict[PSM[:2]]]
            spectrumSpecificFeatureList += [pathScore, pathScore/peptLength, pathScore/scoreStats[PSM[:2]]['Maximum Path Score'], -spectrumOrderedScoreStats[i]['PSM Rankings'][PSMIndexDict[PSM[:2]]], spectrumOrderedScoreStats[i]['Delta Scores'][PSMIndexDict[PSM[:2]]], numNegativePRMs, numNegativePRMs/float(peptLength-1), spectrumOrderedScoreStats[i]['Min Node Scores'][PSMIndexDict[PSM[:2]]]]
            
            # Add mass deviation from true peptide mass to feature list
            precMass = scanFDict[scan]['precMass']
            spectrumSpecificFeatureList += [abs(truePMs[PSM[:2]] + Constants.mods['H2O'] + Constants.mods['H+'] - precMass)]
        
            peakAnnotationMassOffsetStats = Discriminator.getPeakAnnotationAndMassOffsetStats(DataFile.getMassIntPairs(scanFDict[scan]['dta']), specs[i], prmLadders[PSM[:2]], pairedPRMLadder, PNet)
            GLFD.addPeakAnnotationStatsToFeatureList(PNet, peakAnnotationMassOffsetStats, spectrumSpecificFeatureList, peptLength)
            GLFD.addMassOffsetStatsToFeatureList(peakAnnotationMassOffsetStats, spectrumSpecificFeatureList)
        
            spectrumSpecificFeatureList += [precMass, GLFD.getChargeStateFromDTAFName(scanFDict[scan]['dta']), peptLength]
            spectrumAndPSMSpecificFeatureDict[(scan, PSM[:2])] = spectrumSpecificFeatureList

        for j, scan in enumerate(heavyScans):
            i = j + len(lightScans)
            
            spectrumSpecificFeatureList = copy.copy(PSMSpecificFeatureList)
            # Add path score (and normalized variants), delta rank, delta score, number of negative PRMs, and minimum node score for spectrum to feature list
            pathScore = spectrumOrderedScoreStats[i]['Path Scores'][PSMIndexDict[PSM[:2]]]
            numNegativePRMs = spectrumOrderedScoreStats[i]['Num Negative PRMs'][PSMIndexDict[PSM[:2]]]
            spectrumSpecificFeatureList += [pathScore, pathScore/peptLength, pathScore/scoreStats[PSM[:2]]['Maximum Path Score'], -spectrumOrderedScoreStats[i]['PSM Rankings'][PSMIndexDict[PSM[:2]]], spectrumOrderedScoreStats[i]['Delta Scores'][PSMIndexDict[PSM[:2]]], numNegativePRMs, numNegativePRMs/float(peptLength-1), spectrumOrderedScoreStats[i]['Min Node Scores'][PSMIndexDict[PSM[:2]]]]
            
            # Add mass deviation from true peptide mass to feature list
            precMass = scanFDict[scan]['precMass']
            spectrumSpecificFeatureList += [abs(truePMs[PSM[:2]] + pairConfig['NMod'] + pairConfig['CMod'] + Constants.mods['H2O'] + Constants.mods['H+'] - precMass)]
            
            peakAnnotationMassOffsetStats = Discriminator.getPeakAnnotationAndMassOffsetStats(DataFile.getMassIntPairs(scanFDict[scan]['dta']), specs[i], prmLadders[PSM[:2]], pairedPRMLadder, PNet)
            GLFD.addPeakAnnotationStatsToFeatureList(PNet, peakAnnotationMassOffsetStats, spectrumSpecificFeatureList, peptLength)
            GLFD.addMassOffsetStatsToFeatureList(peakAnnotationMassOffsetStats, spectrumSpecificFeatureList)
            
            spectrumSpecificFeatureList += [precMass, GLFD.getChargeStateFromDTAFName(scanFDict[scan]['dta']), peptLength]
            spectrumAndPSMSpecificFeatureDict[(scan, PSM[:2])] = spectrumSpecificFeatureList

    return spectrumAndPSMSpecificFeatureDict
    
def getScanScoreDictLinearDiscriminator(LADSSeqInfo, seqEntry, scanFDict, linearDiscriminatorVec, pairConfig, PNet):
    scanScoreDict = {}
    spectrumAndPSMSpecificFeatureDict = getSpectrumAndPSMFeatureDict(LADSSeqInfo, seqEntry, scanFDict, pairConfig, PNet)
    # Now get PSM with highest rank score for each scan
    fullPSMList = LADSSeqInfo[seqEntry]

    for scan in lightScans + heavyScans:
        discScoreList = []
        for PSM in fullPSMList:
            featureList = spectrumAndPSMSpecificFeatureDict[(scan, PSM[:2])]
            discScore = calculateDiscScoreFromLinearDiscriminator(linearDiscriminatorVec, np.array(featureList))
            discScoreList += [discScore]

        highestDiscScoreInd = np.argmax(discScoreList)
        scanScoreDict[scan] = {'Seq': (fullPSMList[highestDiscScoreInd][1], fullPSMList[highestDiscScoreInd][2]), 'Raw Score': fullPSMList[highestDiscScoreInd][0], 'Post Score': discScoreList[highestDiscScoreInd]}

    return scanScoreDict

def getScanScoreDictSVM(LADSSeqInfo, seqEntry, scanFDict, svmModel, svmRange, pairConfig, PNet, desired_feats=None):
    scanScoreDict = {}

    spectrumAndPSMSpecificFeatureDict = getSpectrumAndPSMFeatureDict(LADSSeqInfo, seqEntry, scanFDict, pairConfig, PNet)
    # Now get PSM with highest rank score for each scan
    fullPSMList = LADSSeqInfo[seqEntry]
    for scan in lightScans + heavyScans:
        xVals = []
        for PSM in fullPSMList:
            featureList = spectrumAndPSMSpecificFeatureDict[(scan, PSM[:2])]
            if desired_feats != None:
                xVals += [dict((i+1, featureList[desired_feats[i] - 1]) for i in range(len(desired_feats)))]
            else:
                xVals += [dict((i+1, featureList[i]) for i in range(len(featureList)))]

        xValsNorm = svmutil.normalize_instances(xVals, svmRange)
        
        probs = zip(*svmutil.svm_predict([0] * len(xValsNorm), xValsNorm, svmModel, '-b 1')[2])[0]
        #probs = zip(*svmutil.svm_predict([0] * len(xValsNorm), xValsNorm, svmModel, '-b 1')[2])[1]
            
        highestProbInd = np.argmax(probs)
        scanScoreDict[scan] = {'Seq': (fullPSMList[highestProbInd][1], fullPSMList[highestProbInd][2]), 'Raw Score': fullPSMList[highestProbInd][0], 'Post Score': probs[highestProbInd]}

    return scanScoreDict

def getScanScoreDictRankBoost(LADSSeqInfo, seqEntry, scanFDict, rankingModel, pairConfig, PNet):
    scanScoreDict = {}

    spectrumAndPSMSpecificFeatureDict = getSpectrumAndPSMFeatureDict(LADSSeqInfo, seqEntry, scanFDict, pairConfig, PNet)
    # Now get PSM with highest rank score for each scan
    fullPSMList = LADSSeqInfo[seqEntry]
    for scan in lightScans + heavyScans:
        rankScoreList = []
        for PSM in fullPSMList:
            featureList = spectrumAndPSMSpecificFeatureDict[(scan, PSM[:2])]
            rankScore = calculateRankScoreFromRankBoostModel(rankingModel, featureList)
            rankScoreList += [rankScore]
            
        highestRankScoreInd = np.argmax(rankScoreList)
        scanScoreDict[scan] = {'Seq': (fullPSMList[highestRankScoreInd][1], fullPSMList[highestRankScoreInd][2]), 'Raw Score': fullPSMList[highestRankScoreInd][0], 'Post Score': rankScoreList[highestRankScoreInd]}

    return scanScoreDict

def calculateDiscScoreFromLinearDiscriminator(wVec, featureVec):
    return np.dot(-1*wVec, featureVec)

def calculateRankScoreFromRankBoostModel(rankingModel, featureList):
    rankScore = 0
    for learner in rankingModel:
    # each entry is a tuple of the form (featureLabel, cutoff, weight)
        rankScore += learner[2] if featureList[learner[0] - 1] > learner[1] else 0

    return rankScore
    
def loadRankingModel(rankingModelLoc):
    rankingModel = []
    modelFile = open(rankingModelLoc)
    for line in modelFile:
        if line[0] == '#':
            continue

        weakLearners = line.split(' ')
        for learner in weakLearners:
            learnerTuple = learner.split(':')
            rankingModel += [(int(learnerTuple[0]), float(learnerTuple[1]), float(learnerTuple[2]))]

    return rankingModel

            

if __name__ == '__main__':
    print 'This program generates a results file containing Raw lads output postscored with the algorithm of choice. The discmodel is a supplied model, if necessary for the postscoring algorithm'
    options = ArgLib.parse(['init', 'ppmstd', 'dtadir', 'lads', 'sequest', 'config', 'model', 'output', 'symbolmap'], optArgs=[{'opts': ('-D', '--discmodel'), 'attrs': {'type': 'string', 'dest': 'discmodel', 'help': 'Model used to calculate discriminant score'}}, {'opts': ('-P', '--pairconfig'), 'attrs': {'type': 'string', 'dest': 'pairconfig', 'help': 'Name of LADS Pair Configuration'}}, {'opts': ('-F', '--featurelist'), 'attrs': {'type': 'string', 'dest': 'featurelist', 'help': 'File containing pickled list of desired features (optional)'}}])
    parent = os.path.abspath(os.pardir)
                           
    PNet = PN.ProbNetwork(options.config, options.model)
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    pairConfigurations = paramsDict['Pair Configurations']

    LADSSeqInfo = GLFD.parseSequenceDTAsLogfile(options.lads)

    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap({'LADS Unit Test': 'LADS'}, symbolMap, paramsDict)
    seqMap = seqMap['LADS Unit Test']

    if options.featurelist:
        with open(options.featurelist) as fin:
            desired_feats = pickle.load(fin)
    else:
        desired_feats = None

    heavySeqMaps = {}
    for confName in pairConfigurations:
        heavySeqMaps[confName] = copy.deepcopy(seqMap)
        heavySeqMaps[confName]['Mods']['N-Term'] = pairConfigurations[confName]['NModSymbol']
        heavySeqMaps[confName]['Mods']['C-Term'] = pairConfigurations[confName]['CModSymbol']

    if options.pairconfig:
        pairConfigName = options.pairconfig
    else:
        pairConfigName = pairConfigurations.keys()[0]
        
    dtaList = glob.glob(options.dtadir + '/*.dta')
    scanFDict = SDTDV.getScanFDict(dtaList)

    if options.sequest:
        SEQUESTMASCOTResults = DataFile.indexDataByKey(DataFile.getScanInfo(options.sequest, delimiter='\t'))
    else:
        SEQUESTMASCOTResults = {}

    outFile = open(options.output, 'w')
                           
#    rankModel = loadRankingModel(options.discmodel)
#    linearDiscVec = np.loadtxt(options.discmodel)
    svmModel = svmutil.svm_load_model(options.discmodel)
    svmRange = svmutil.load_ranges(os.path.splitext(options.discmodel)[0] + '.range')
    """
    cols = ['ScanF', 'M+H', 'PEAKS ALC', 'PEAKS TLC', 'PEAKS Sequence', 'SEQUEST XCorr', 'MASCOT Ion Score', 'SEQUEST MASCOT Sequence', 'Accuracy', 'Precision']
    outFile.write('\t'.join([col for col in cols]) + '\n')
    PEAKSInfo = DataFile.indexDataByKey(DataFile.getScanInfo(options.lads, delimiter='\t'), key='ScanF')

    for scan in PEAKSInfo:
        if scan not in SEQUESTMASCOTResults:
            continue

        scanData = {'ScanF': scan}
        scanData['M+H'] = scanFDict[scan]['precMass']
        scanData['PEAKS ALC'] = PEAKSInfo[scan]['ALC (%)']
        scanData['PEAKS TLC'] = PEAKSInfo[scan]['TLC']
        scanData['PEAKS Sequence'] = PEAKSInfo[scan]['Peptide']
        
        comp = An.comparePeptideResults(PEAKSInfo[scan]['Peptide'], SEQUESTMASCOTResults[scan]['Peptide'], ambigEdges1=[], ambigEdges2=[], ppm=20)
        scanData['SEQUEST XCorr'] = SEQUESTMASCOTResults[scan]['SEQUEST XCorr']
        scanData['MASCOT Ion Score'] = SEQUESTMASCOTResults[scan]['MASCOT Ion Score']
        scanData['SEQUEST MASCOT Sequence'] = SEQUESTMASCOTResults[scan]['Peptide']
        scanData['Accuracy'] = comp[0]
        scanData['Precision'] = comp[1]

        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
    
    outFile.close()
    """

    cols = ['ScanF', 'M+H', 'LADS Raw Score', 'LADS Sequence', 'LADS Ambig Edges', 'LADS Post Score', 'SEQUEST XCorr', 'MASCOT Ion Score', 'SEQUEST MASCOT Sequence', 'Accuracy', 'Precision']
    outFile.write('\t'.join([col for col in cols]) + '\n')

    for seqEntry in LADSSeqInfo:
        lightScans = seqEntry[0]
        heavyScans = seqEntry[1]

        scanScoreDict = getScanScoreDictSVM(LADSSeqInfo, seqEntry, scanFDict, svmModel, svmRange, pairConfigurations[pairConfigName], PNet, desired_feats = desired_feats)
        
#        scanScoreDict = getScanScoreDictRankBoost(LADSSeqInfo, seqEntry, scanFDict, rankModel, pairConfigurations['lightdimethyl_heavydimethyl'], PNet)
#        scanScoreDict = getScanScoreDictClusterNormScore(LADSSeqInfo, seqEntry)

        for i, scan in enumerate(lightScans):

            scanData = {'ScanF': scan}
                        
            lightSeq = An.preprocessSequence(scanScoreDict[scan]['Seq'][0], seqMap, ambigEdges=scanScoreDict[scan]['Seq'][1])
            scanData['LADS Sequence'] = lightSeq
            scanData['LADS Ambig Edges'] = scanScoreDict[scan]['Seq'][1]
            scanData['LADS Raw Score'] = scanScoreDict[scan]['Raw Score']
            scanData['LADS Post Score'] = scanScoreDict[scan]['Post Score']
            scanData['M+H'] = scanFDict[scan]['precMass']

            try:
                comp = An.comparePeptideResults(lightSeq, SEQUESTMASCOTResults[scan]['Peptide'], ambigEdges1=scanScoreDict[scan]['Seq'][1], ambigEdges2=[], ppm=20)            
                scanData['SEQUEST XCorr'] = SEQUESTMASCOTResults[scan]['SEQUEST XCorr']
                scanData['MASCOT Ion Score'] = SEQUESTMASCOTResults[scan]['MASCOT Ion Score']
                scanData['SEQUEST MASCOT Sequence'] = SEQUESTMASCOTResults[scan]['Peptide']
                scanData['Accuracy'] = comp[0]
                scanData['Precision'] = comp[1]
            except KeyError:
                scanData['SEQUEST XCorr'] = None
                scanData['MASCOT Ion Score'] = None
                scanData['SEQUEST MASCOT Sequence'] = None
                scanData['Accuracy'] = None
                scanData['Precision'] = None
                
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
            
        for i, scan in enumerate(heavyScans):

            scanData = {'ScanF': scan}

            heavySeq = An.preprocessSequence(scanScoreDict[scan]['Seq'][0], heavySeqMaps[pairConfigName], replaceExistingTerminalMods=True, ambigEdges=scanScoreDict[scan]['Seq'][1])
            scanData['LADS Sequence'] = heavySeq
            scanData['LADS Ambig Edges'] = scanScoreDict[scan]['Seq'][1]
            scanData['LADS Raw Score'] = scanScoreDict[scan]['Raw Score']
            scanData['LADS Post Score'] = scanScoreDict[scan]['Post Score']
            scanData['M+H'] = scanFDict[scan]['precMass']
            
            try:
                comp = An.comparePeptideResults(heavySeq, SEQUESTMASCOTResults[scan]['Peptide'], ambigEdges1=scanScoreDict[scan]['Seq'][1], ambigEdges2=[], ppm=20)
                scanData['SEQUEST XCorr'] = SEQUESTMASCOTResults[scan]['SEQUEST XCorr']
                scanData['MASCOT Ion Score'] = SEQUESTMASCOTResults[scan]['MASCOT Ion Score']
                scanData['SEQUEST MASCOT Sequence'] = SEQUESTMASCOTResults[scan]['Peptide']
                scanData['Accuracy'] = comp[0]
                scanData['Precision'] = comp[1]
            except KeyError:
                scanData['SEQUEST XCorr'] = None
                scanData['MASCOT Ion Score'] = None
                scanData['SEQUEST MASCOT Sequence'] = None
                scanData['Accuracy'] = None
                scanData['Precision'] = None

            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

