
import os
import sys

sys.path.insert(1, os.path.abspath('../'))


import ArgLib
import DeNovoSequencer as DNS
import Analytics as An
import DataFile
import CompareSearches as CS
import SequenceDTAsTDV as SDTDV
import Constants
import ProbNetwork as PN
import Discriminator

import numpy
import glob
import numpy as np
import os
import pickle
import re
import copy
"""
Other possible features:
    - average path score (normalized by peptide length) for [original scoring, consensus rescore, considered spectrum]
    - normalized path score (normalized by cluster size) for [original scoring, consensus rescore]
    - minimum path score (of all members of cluster)
    - minimal cleavage scores (score of lowest scoring PRM)
    - lowest score / highest score (of all members in cluster)
    - norm score - score/highest score (of all members in cluster)
    - number of PRMS with negative score
    - % PRMs with negative score (normalized by peptide length)
    - delta mass - abs(predicted - observed peptide mass)
    - delta rank (change in rank from highest scoring peptide) for [original scoring, consensus rescore, considered spectrum]
    - delta score (change in score from highest scoring peptide) for [original scoring, consensus rescore, considered spectrum]    
"""
def compileScoreStats(scoreStats, specs, PSMList):
    spectrumStats = {}
    totalPathStats = {'PSM Rankings': [], 'Path Scores': [], 'Max Cluster Path Score': 0}

    for i in range(len(specs)):
        spectrumStats[i] = {'PSM Rankings': [], 'Num Negative PRMs': [], 'Path Scores': [], 'Delta Scores': [], 'Min Node Scores': []}

    for PSM in PSMList:
        totalPathStats['Path Scores'] += [scoreStats[PSM]['Total Path Score']]
        
        for i in range(len(specs)):
            spectrumStats[i]['Path Scores'] += [scoreStats[PSM]['Spectrum Specific Score Statistics'][i]['Score']]
            spectrumStats[i]['Min Node Scores'] += [scoreStats[PSM]['Spectrum Specific Score Statistics'][i]['Min Node Score']]
            spectrumStats[i]['Num Negative PRMs'] += [scoreStats[PSM]['Spectrum Specific Score Statistics'][i]['Negative']]

    totalPathStats['PSM Rankings'] = np.argsort(np.argsort(-1*np.array(totalPathStats['Path Scores'])))
    totalPathStats['Max Cluster Path Score'] = max(totalPathStats['Path Scores'])
    
    for i in spectrumStats:
        spectrumStats[i]['PSM Rankings'] = np.argsort(np.argsort(-1*np.array(spectrumStats[i]['Path Scores'])))
        spectrumStats[i]['Delta Scores'] = np.array(spectrumStats[i]['Path Scores']) - max(spectrumStats[i]['Path Scores'])

    return spectrumStats, totalPathStats

def getChargeStateFromDTAFName(dta):
    dtaBase, ext = os.path.splitext(dta)
    return int(os.path.splitext(dtaBase)[1][1:])

def parseSequenceDTAsLogfile(seqLogfile):
    seqLogfilehandler = open(seqLogfile)

    line = None
    seqData = {}
    while line != '':
        line = seqLogfilehandler.readline()
        if len(re.findall('shared peaks ratio', line)) > 0:
            tempScanFs = re.findall('\[.*?\]', line)
            scanFs = tuple([tuple(eval(scans)) for scans in tempScanFs])

            seqData[scanFs] = []
            seqLogfilehandler.readline()
            line = seqLogfilehandler.readline()
            while len(re.findall('Score', line)) > 0:
                print line
                PSMArr = line.split('  ')
                try:
                    score = float(PSMArr[1][:-4])

                    if PSMArr[-1][0] != '[':
                        seq = PSMArr[-1][:-1]
                    else:
                        seq = eval(PSMArr[-1][:-1])

                    ambigEdges = []
                    if type(seq) == list:
                        for i in range(len(seq)):
                            if type(seq[i]) == tuple:
                                ambigEdges.extend([seq[i]])
                                seq[i] = '-'
                        seq = ''.join(seq)

                    seqData[scanFs] += [(score, seq, ambigEdges)]
                except ValueError:
                    pass

                line = seqLogfilehandler.readline()

            # Make sure sequences were found
            if len(seqData[scanFs]) == 0:
                del seqData[scanFs]
            
    return seqData
                    
def addClusterPairingStatsToFeatureList(clusterPairingStats, featureList):                                                        
    featureList += [clusterPairingStats['Shared Peaks Ratio']]

    numPairedIons = sum([clusterPairingStats['Pair Type Stats'][pairedIonType] for pairedIonType in Discriminator.pairTypes])
    featureList += list(np.array([clusterPairingStats['Pair Type Stats'][pairedIonType] for pairedIonType in Discriminator.pairTypes], dtype=np.float)/numPairedIons)
    
def addPairedPRMStatsToFeatureList(pairedPRMStats, featureList, prmLadderLen):
    prmLadderLen = float(prmLadderLen)
    featureList += [len(pairedPRMStats['Paired PRM Ladder']), len(pairedPRMStats['Paired PRM Ladder'])/prmLadderLen]
    featureList += [pairedPRMStats['Num Paired PRMs With Majority Evidence Light'], pairedPRMStats['Num Paired PRMs With Majority Evidence Light']/prmLadderLen, pairedPRMStats['Num Paired PRMs With Majority Evidence Heavy'], pairedPRMStats['Num Paired PRMs With Majority Evidence Heavy']/prmLadderLen]

    featureList += list(np.array([pairedPRMStats['Pair Type Stats for Paired PRMs'][pairedIonType] for pairedIonType in Discriminator.pairTypes], dtype=np.float)/prmLadderLen)

"""
-# annotated peaks in top 25, 50 peaks
-# annotated peaks in top 25, 50 peaks explained by paired PRMs
-% explained intensity
- # of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc.
- % of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc. (normalized by peptide length)
- # peak annotations with 'only neutral losses' (i.e., b-H2O found with no b)
"""
def addPeakAnnotationStatsToFeatureList(PNet, peakAnnotationStats, featureList, peptLength):
    featureList += [peakAnnotationStats['Percent Explained Intensity'], peakAnnotationStats['Annotated Peaks 25'], peakAnnotationStats['Annotated Peaks 50'], peakAnnotationStats['Annotated Peaks Paired 25'], peakAnnotationStats['Annotated Peaks Paired 50']]

    for ion in list(PNet._ions) + ['imm', 'internal', 'only neut']:
        featureList += [peakAnnotationStats['Fragment Ion Counts'][ion], peakAnnotationStats['Fragment Ion Counts'][ion]/float(peptLength-1)]

"""
- average mass offset for fragment x = b, y
- max mass offset for fragment x = b, y
- max relative offset for fragment x = b, y (i.e., b_n+1 - b_n - mass(intervening polypeptide))
"""
def addMassOffsetStatsToFeatureList(massOffsetStats, featureList):
   for ion in ['b', 'y']:
       featureList += [massOffsetStats['Average Mass Offsets'][ion], massOffsetStats['Max Mass Offsets'][ion], massOffsetStats['Max Rel Mass Offsets'][ion]]

def generateFeatureNames(PNet):
    featureNames = []

    featureNames += ['Shared Peaks Ratio']
    featureNames += ['Proportion %s Ion Pairs in allPairedIonsDict (Potential PRMs normalized)' % (pairedIonType,) for pairedIonType in Discriminator.pairTypes]
    featureNames += ['LADS PScore', 'LADS PScore (Length Normalized)', 'LADS PScore (Cluster size normalized)', 'LADS PScore Delta Rank', 'LADS PScore Delta Score']
    featureNames += ['Cluster Total Path Score', 'Cluster Total Path Score (Length Normalized)', 'Cluster Total Path Score (Cluster size normalized)', 'Cluster Total Path Score Delta Rank', 'Cluster Total Path Score Delta Score', 'Cluster Total Path Minimum Node Score']
    featureNames += ['Minimum Path Score', 'Minimum Path Score (Length Normalized)', 'Maximum Path Score', 'Maximum Path Score (Length Normalized)', 'Minimum Path Score/Maximum Path Score']
    featureNames += ['Highest Rank - Lowest Rank for PSM Over Cluster']
    featureNames += ['Num Forbidden Node Pairs', 'Proportion Forbidden Node PRMs (PRM Ladder Length Normalized)']
    featureNames += ['Number of Ambiguous Edges']
    featureNames += ['# PRMs with evidence from all members of cluster', '% PRMs with evidence from all members of cluster (PRM Ladder Length Normalized)', '# PRMs with evidence from majority of cluster', '% PRMs with evidence from majority of cluster (PRM Ladder Length Normalized)', '# PRMs with no evidence from any members of cluster', '% PRMs with no evidence from any members of cluster (PRM Ladder Length Normalized)']
    featureNames += ['Num Paired PRMs', 'Num Paired PRMs (PRM Ladder Length Normalized)']
    featureNames += ['Num Paired PRMs With Majority Peaks Light', 'Num Paired PRMs With Majority Peaks Light (PRM Ladder Length Normalized)', 'Number Paired PRMs With Majority Peaks Heavy', 'Number Paired PRMs With Majority Peaks Heavy (PRM Ladder Length Normalized)']
    featureNames += ['Proportion %s-based PRMs in Paired PRM Ladder (PRM Ladder Length Normalized)' % (pairedIonType,) for pairedIonType in Discriminator.pairTypes]
    featureNames += ['Spectrum Path Score', 'Spectrum Path Score (Peptide Length Normalized)', 'Spectrum Path Score/Maximum Path Score Over Cluster', 'Spectrum Path Score Delta Rank', 'Spectrum Path Score Delta Score', 'Spectrum Path Score Number Negative Scoring PRMs', 'Spectrum Path Score Proportion Negative Scoring PRMs (PRM Ladder Length Normalized)', 'Spectrum Path Minimum Node Score']
    featureNames += ['Absolute Deviation Observed - Predicted Precursor Ion Mass']
    featureNames += ['% Explained TIC', 'Num Annotated Peaks of Top 25', 'Num Annotated Peaks of Top 50', 'Num Annotated Peaks of Top 25 By Paired PRMs', 'Num Annotated Peaks of Top 50 by Paired PRMs']
    for ion in list(PNet._ions) + ['imm', 'internal', 'only neut']:
        featureNames += ['# Peak annotations for fragment ion %s' % (ion,), 'Proportion Peak Annotations for fragment ion %s (Normalized by PRM Ladder Length)' %(ion,)]
    for ion in ['b', 'y']:
        featureNames += ['Average Mass Deviation for %s-ions' % (ion,), 'Maximum Amplitude Mass Deviation For %s-ions' % (ion,), 'Maximum Amplitude Relative Mass Deviation for %s-ions' % (ion,)]
    featureNames += ['Precursor Ion Mass', 'Charge State', 'Peptide Length']

    return featureNames

def writeFeatures(featureList, rank, qid, outFile, comment=""):
    outFile.write('%i ' % (rank,) + ' '.join(['%i:%f' % (i+1, feature) for i, feature in enumerate(featureList)]) + '# %s\n' % (comment,))

def printFeatureNames(featureNames):
    for i, feature in enumerate(featureNames):
        print '%i. %s' % (i+1, feature)

def printFeatures(featureNames, featureList):
    for i, feature in enumerate(featureNames):
        print '%i. %s: %f' % (i+1, feature, featureList[i])

if __name__ == '__main__':
    print 'This program generates LETOR format training data for the training of a discriminator. dtadir is of the formate {/loc of dtadir: (loc of LADS SequenceDTAsTDV.py LOG file, loc of combined SEQUEST-MASCOT database results'
    options = ArgLib.parse(['init', 'dtadir', 'ppmstd', 'symbolmap', 'output', 'model', 'config'])

    paramsDict = ArgLib.parseInitFile(options.init, options)
    pairConfigurations = paramsDict['Pair Configurations']
    ppm = float(options.ppmstd)

    dtadirInfo = eval(options.dtadir)

    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap({'LADS Unit Test': 'LADS'}, symbolMap, paramsDict)
    seqMap = seqMap['LADS Unit Test']

    PNet = PN.ProbNetwork(options.config, options.model)
    outFile = open(options.output, 'w')

    featureNames = generateFeatureNames(PNet)
    #printFeatureNames(featureNames)

    heavySeqMaps = {}
    for confName in pairConfigurations:
        heavySeqMaps[confName] = copy.deepcopy(seqMap)
        heavySeqMaps[confName]['Mods']['N-Term'] = pairConfigurations[confName]['NModSymbol']
        heavySeqMaps[confName]['Mods']['C-Term'] = pairConfigurations[confName]['CModSymbol']


    pairConfig = paramsDict['Pair Configurations']['silac_light_heavy']
    for dtadir in dtadirInfo:
        dtaList = glob.glob(dtadir + '/*.dta')
        print 'Num DTAs', len(dtaList)
        scanFDict = SDTDV.getScanFDict(dtaList)
        print 'Getting LADS Seq Info for %s' % (dtadir,)
        LADSSeqInfo = parseSequenceDTAsLogfile(dtadirInfo[dtadir][0])
        print 'Got LADS Seq Info for %s %i' % (dtadir, len(LADSSeqInfo))
        SEQUESTMASCOTResults = DataFile.indexDataByKey(DataFile.getScanInfo(dtadirInfo[dtadir][1], delimiter='\t'))
        print 'Num DB Results', len(SEQUESTMASCOTResults)
        for seqEntry in LADSSeqInfo:
            print seqEntry
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
            addClusterPairingStatsToFeatureList(clusterPairingStats, featureList)
            
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
            spectrumOrderedScoreStats, clusterScoreStats = compileScoreStats(scoreStats, specs, PSMList)

            PSMIndexDict = dict([(PSM, i) for i, PSM in enumerate(PSMList)])
            for i, PSM in enumerate(LADSSeqInfo[seqEntry]):
                PSMSpecificFeatureList = copy.copy(featureList)
                lightSeq = An.preprocessSequence(PSM[1], seqMap, ambigEdges=PSM[2])
                heavySeq = An.preprocessSequence(PSM[1], heavySeqMaps['silac_light_heavy'], replaceExistingTerminalMods=True, ambigEdges=PSM[2])
                
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
                addPairedPRMStatsToFeatureList(pairedPRMStats, PSMSpecificFeatureList, len(prmLadders[PSM[:2]]))

                pairedPRMLadder = pairedPRMStats['Paired PRM Ladder']                
                for i, scan in enumerate(lightScans):
                    if int(scan) not in SEQUESTMASCOTResults:
                        continue
                    
                    spectrumSpecificFeatureList = copy.copy(PSMSpecificFeatureList)

                    # Add path score (and normalized variants), delta rank, delta score, number of negative PRMs, and minimum node score for spectrum to feature list
                    pathScore = spectrumOrderedScoreStats[i]['Path Scores'][PSMIndexDict[PSM[:2]]]
                    numNegativePRMs = spectrumOrderedScoreStats[i]['Num Negative PRMs'][PSMIndexDict[PSM[:2]]]
                    spectrumSpecificFeatureList += [pathScore, pathScore/peptLength, pathScore/scoreStats[PSM[:2]]['Maximum Path Score'], -spectrumOrderedScoreStats[i]['PSM Rankings'][PSMIndexDict[PSM[:2]]], spectrumOrderedScoreStats[i]['Delta Scores'][PSMIndexDict[PSM[:2]]], numNegativePRMs, numNegativePRMs/float(peptLength-1), spectrumOrderedScoreStats[i]['Min Node Scores'][PSMIndexDict[PSM[:2]]]]
                    
                    # Add mass deviation from true peptide mass to feature list
                    precMass = scanFDict[scan]['precMass']
                    spectrumSpecificFeatureList += [abs(truePMs[PSM[:2]] + Constants.mods['H2O'] + Constants.mods['H+'] - precMass)]
                    
                    peakAnnotationMassOffsetStats = Discriminator.getPeakAnnotationAndMassOffsetStats(DataFile.getMassIntPairs(scanFDict[scan]['dta']), specs[i], prmLadders[PSM[:2]], pairedPRMLadder, PNet)
                    addPeakAnnotationStatsToFeatureList(PNet, peakAnnotationMassOffsetStats, spectrumSpecificFeatureList, peptLength)
                    addMassOffsetStatsToFeatureList(peakAnnotationMassOffsetStats, spectrumSpecificFeatureList)

                    spectrumSpecificFeatureList += [precMass, getChargeStateFromDTAFName(scanFDict[scan]['dta']), peptLength]
                    
                    comp = An.comparePeptideResults(lightSeq, SEQUESTMASCOTResults[scan]['Peptide'], ambigEdges1=PSM[2], ambigEdges2=[], ppm=5)
                    acc, prec = comp[0], comp[1]

#                    print lightSeq, SEQUESTMASCOTResults[scan]['Peptide']
#                    print "acc, prec", acc, prec

                    if prec < 1:
                        rank = -1
                    else:
                        rank = 1

                    writeFeatures(spectrumSpecificFeatureList, rank, 1, outFile, comment="Scan %i From DTA Directory %s" % (int(scan), dtadir))
#                    printFeatures(featureNames, spectrumSpecificFeatureList)

                for j, scan in enumerate(heavyScans):
                    if int(scan) not in SEQUESTMASCOTResults:
                        continue
                    #adjust offset in specs array
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
                    addPeakAnnotationStatsToFeatureList(PNet, peakAnnotationMassOffsetStats, spectrumSpecificFeatureList, peptLength)
                    addMassOffsetStatsToFeatureList(peakAnnotationMassOffsetStats, spectrumSpecificFeatureList)

                    spectrumSpecificFeatureList += [precMass, getChargeStateFromDTAFName(scanFDict[scan]['dta']), peptLength]
                    
                    try:
                        comp = An.comparePeptideResults(heavySeq, SEQUESTMASCOTResults[scan]['Peptide'], ambigEdges1=PSM[2], ambigEdges2=[], ppm=5)
                    except KeyError:
                        continue
                    
                    acc, prec = comp[0], comp[1]

                    if prec < 1:
                        rank = -1
                    else:
                        rank = 1

                    writeFeatures(spectrumSpecificFeatureList, rank, 1, outFile, comment="Scan %i From DTA Directory %s" % (int(scan), dtadir))
                    #printFeatures(featureNames, spectrumSpecificFeatureList)

    outFile.close()
                    
