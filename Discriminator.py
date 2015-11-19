

"""
This module calculates various features of putative peptide-spectral matches and uses a discriminator to rescore these PSMs

Possible pair/cluster-based features:
-shared peaks ratio - # shared peaks / total number of peaks in cluster - Done
-num paired PRMs - Done
-% paired PRMs - Done
-# paired PRMs with evidence from majority of light cluster - Done
-# paired PRMs with evidence from majority of heavy cluster - Done
-% paired PRMs with evidence from  majority of light, heavy cluster - Done
- num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc. num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc in light, heavy pair (normalized by number of paired ions, also include raw counts?) - Done
-num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc. num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc. in paired PRM ladder (normalized by peptide length, also include raw counts?) - Done

-# PRMs with evidence from all members of cluster - Done
- % PRMs with evidence from all members of cluster - Done
- # PRMs with evidence from majority of cluster - Done
- % PRMs with evidence from majority of cluster - Done
- # PRMs with no evidence from any member of cluster - Done
- % PRMs with no evidence from any member of cluster - Done

Other possible features:
- path score for [original scoring - Done, consensus rescore - Done, considered spectrum - Done]
- average path score (normalized by peptide length) for [original scoring - Done, consensus rescore - Done, considered spectrum - Done]
- normalized path score (normalized by cluster size) for [original scoring - Done, consensus rescore - Done]
- minimum path score (of all members of cluster) (normalized by peptide length) - Done, Done
- maximum path score (of all members of cluster) (normalized by peptide length) - Done, Done
- minimal cleavage scores (score of lowest scoring PRM) for [consensus - Done, considered spectrum - Done]
- lowest score / highest score (of all members in cluster) - Done
- norm score - score/highest score (of all members in cluster) - Done
- number of PRMS with negative score for considered spectrum - Done
- % PRMs with negative score for considered spectrum - Done
- delta mass - abs(predicted - observed peptide mass) - Done
- delta rank (change in rank from highest scoring peptide) for [original scoring - Done, consensus rescore - Done, considered spectrum - Done]
- delta score (change in score from highest scoring peptide) for [original scoring - Done, consensus rescore - Done, considered spectrum - Done]
- highest rank - lowest rank (of peptide for all spectra in cluster after rescoring) - Done
- number forbidden node pairs (normalized by peptide length) - Done, Done

-# annotated peaks in top 25, 50 peaks - Done, Done
-# annotated peaks in top 25, 50 peaks explained by paired PRMs - Done, Done
-% explained intensity - Done
- # of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc. - Done
- % of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc. (normalized by peptide length) - Done
- # peak annotations with 'only neutral losses' (i.e., b-H2O found with no b) - Done
- average mass offset for fragment x = b, y - Done
- max mass offset for fragment x = b, y - Done
- max relative offset for fragment x = b, y (i.e., b_n+1 - b_n - mass(intervening polypeptide)) - Done


-average/minimal triplet category
-category of amino acid triplets at N/C-terminus (adapt to specificity?)

-number ambiguous edges - Done
-observed peptide mass (Frank trains different RankBoost model for different peptide masses)
-observed peptide charge state
-peptide length

Note: for all attributes which assume only one spectrum, take all answers from 'best' spectrum (How to determine?)
"""


import numpy as np
import networkx as nx

import DeNovoSequencer as DNS
import DataFile
import ProbNetwork as PN
import Analytics as An
import SpectraAlignment as SA
import Constants

pairTypes = ['by_by', 'by_b', 'by_y', 'b_by', 'y_by', 'b_y', 'b_b', 'y_b', 'y_y']

def closestInSortedArray(val, array):
    #    print val, array
    ind = np.searchsorted(array, val)
    try:
        leftTol = abs(val - array[ind-1])
    except IndexError:
        leftTol = 100000
        
    try:
        rightTol = abs(array[ind] - val)
    except IndexError:
        rightTol = 100000
        
    return array[ind - int(leftTol < rightTol)]
                                                        

def getAllPairedIonsDict(lightMergedSpec, heavyMergedSpec, lightPrecMass, pairConfig, epSTD=0.01):
 
    NTermTable, CTermTable = SA.getNandCIons(lightMergedSpec, heavyMergedSpec, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD)
    NCrossTable, CCrossTable = SA.getCrossPairedIons(lightMergedSpec, heavyMergedSpec, lightPrecMass, Nmod=pairConfig['NMod'], Cmod=pairConfig['CMod'], epsilon=2*epSTD)
    
    NTermIonDict = SA.prepIonTableForAddition(NTermTable, ['b', 'b'])
    CTermIonDict = SA.prepIonTableForAddition(CTermTable, ['y', 'y'])
    NCrossIonDict = SA.prepIonTableForAddition(NCrossTable, ['y', 'b'])
    CCrossIonDict = SA.prepIonTableForAddition(CCrossTable, ['b', 'y'])

    return SA.addDicts(SA.reverseDict(SA.addDicts(NTermIonDict, CCrossIonDict)), SA.reverseDict(SA.addDicts(NCrossIonDict, CTermIonDict)))

"""
This method calculates the following features:
-shared peaks ratio - # shared peaks / total number of peaks in cluster
-num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc. for spectra
-Also calculates information required by getPairedPRMStats to calculate pairedPRMs for a prmLadder and their associated pair type and cluster stats
"""
def getClusterPairingStats(lightSpecs, heavySpecs, lightPrecMass, pairConfig, epSTD = 0.01):
    lightMergedSpec = SA.mergeSpectra(lightSpecs, epsilon=2*epSTD)
    heavyMergedSpec = SA.mergeSpectra(heavySpecs, epsilon=2*epSTD)

    allPairedIonsDict = getAllPairedIonsDict(lightMergedSpec, heavyMergedSpec, lightPrecMass, pairConfig, epSTD)

    specPairedPRMs = {}
    pairTypeCount = {}
    for pairType in pairTypes:
        specPairedPRMs[pairType] = []
        pairTypeCount[pairType] = 0
        
    numLightInds = 0
    numHeavyInds = 0
    for heavyIons in allPairedIonsDict:
        deltaMasses = []
        pairType = {'light': [], 'heavy': []}
        for ion in allPairedIonsDict[heavyIons]:
            pairType['light'] += [ion[1]]
            deltaMasses += [PN.ProbNetwork.deltaRules[ion[1]](lightPrecMass-Constants.mods['H+']-Constants.mods['H2O'], lightMergedSpec[ion[0]][0], 0, 0)]
        for ion in heavyIons:
            pairType['heavy'] += [ion[1]]
            deltaMasses += [PN.ProbNetwork.deltaRules[ion[1]](lightPrecMass-Constants.mods['H+']-Constants.mods['H2O'], heavyMergedSpec[ion[0]][0], pairConfig['NMod'], pairConfig['CMod'])]
            
        pairTypeString = ''.join(pairType['light']) + '_' + ''.join(pairType['heavy'])
        specPairedPRMs[pairTypeString] += [(sum(deltaMasses)/len(deltaMasses), (allPairedIonsDict[heavyIons], heavyIons))]
        pairTypeCount[pairTypeString] += 1
        
        numLightInds += len(allPairedIonsDict[heavyIons])
        numHeavyInds += len(heavyIons)

    sharedPeaksRatio = float(numLightInds + numHeavyInds)/(lightMergedSpec.shape[0] + heavyMergedSpec.shape[0])

    return {'Cluster Paired PRM Information': specPairedPRMs, 'Shared Peaks Ratio': sharedPeaksRatio, 'Pair Type Stats': pairTypeCount, 'Light Merged Spec': lightMergedSpec, 'Heavy Merged Spec': heavyMergedSpec, 'Num Paired Ions': numLightInds + numHeavyInds}

"""
This method calculates the following features:
-num paired PRMs
-# paired PRMs with evidence from majority of light cluster
-# paired PRMs with evidence from majority of heavy cluster
-% paired PRMs with evidence from  majority of light, heavy cluster
- num by_by pairs, by_y pairs, b_by pairs, b_b, y_y, etc. for paired PRMs
"""
def getPairedPRMStats(prmLadder, lightMergedSpec, heavyMergedSpec, lightSpecs, heavySpecs, specPairedPRMs, epSTD=0.01):
    pairTypeResults = {}
    pairedPRMs = []
    pairedPRMsEvidenceLight = 0
    pairedPRMsEvidenceHeavy = 0
    for pairType in specPairedPRMs:
        pairTypeResults[pairType] = 0
        for prmInfo in specPairedPRMs[pairType]:
            arrayVal = closestInSortedArray(prmInfo[0], prmLadder)
            if abs(arrayVal - prmInfo[0]) < 2*epSTD:
                pairedPRMs += [arrayVal]
                pairTypeResults[pairType] += 1
                for ion in prmInfo[1][0]:
                    if 2*lightMergedSpec[ion[0]][2] >= len(lightSpecs):
                        pairedPRMsEvidenceLight += 1
                        break
                for ion in prmInfo[1][1]:
                    if 2*heavyMergedSpec[ion[0]][2] >= len(heavySpecs):
                        pairedPRMsEvidenceHeavy += 1
                        break
                        
    return {'Paired PRM Ladder': pairedPRMs, 'Pair Type Stats for Paired PRMs': pairTypeResults, 'Num Paired PRMs With Majority Evidence Light': pairedPRMsEvidenceLight, 'Num Paired PRMs With Majority Evidence Heavy': pairedPRMsEvidenceHeavy}

"""
This method calculates the following features:
-# annotated peaks in top 25, 50 peaks
-# annotated peaks in top 25, 50 peaks explained by paired PRMs
-% explained intensity
- # of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc.
- % of peak annotations for fragment x = b; y; a; y+2; y-H2O, internal, immonium, etc. (normalized by peptide length)
- # peak annotations with 'only neutral losses' (i.e., b-H2O found with no b)
- average mass offset for fragment x = b, y
- max mass offset for fragment x = b, y
- max relative offset for fragment x = b, y (i.e., b_n+1 - b_n - mass(intervening polypeptide))
"""
def getPeakAnnotationAndMassOffsetStats(spectrumData, spec, PRMLadder, pairedPRMLadder, PNet, ppmstd=5):
    explainedPeakInds = set()
    explainedPeakIndsByPairedPRMs = set()

    indsByInt = np.argsort(-1*spectrumData[:,1])
    top25Inds = set(indsByInt[:25])
    top50Inds = set(indsByInt[:50])

    epSTD = ppmstd*spec._mH*10**-6

    ionCounts = {}
    for ion in PNet._ions:
        ionCounts[ion] = 0
    ionCounts['internal'] = 0
    ionCounts['imm'] = 0
    ionCounts['only neut'] = 0
    
    massOffsets = {'b': [], 'y': []}

    pairedPRMInds = set([i for i in range(len(PRMLadder)) if PRMLadder[i] in pairedPRMLadder])

    #calculate statistics for fragment ions in model
    ionIntsByType = {}
    for i, prm in enumerate(PRMLadder):
        ionInts = PNet.getIonIntensitiesForPRM(spec, prm)
        for j, ion in enumerate(PNet._ions):
            ionIntsByType[ion] = ionInts[j]
            if ionInts[j] > 0:
                ionMass = PNet.ionRules[PNet._ions[j]](spec.pm, prm, spec.Nmod, spec.Cmod)
                closestInd = np.argmin(np.abs(spectrumData[:,0]-ionMass))
                if abs(spectrumData[closestInd][0] - ionMass) < 5*epSTD:
                    ionCounts[PNet._ions[j]] += 1

                    explainedPeakInds.add(closestInd)

                    if i in pairedPRMInds:
                        explainedPeakIndsByPairedPRMs.add(closestInd)
                        
                    if PNet._ions[j] == 'b' or PNet._ions[j] == 'y':
                        massOffsets[PNet._ions[j]] += [spectrumData[closestInd][0] - ionMass]

        for ion in ['b','y']:
            if ionIntsByType[ion] == 0 and ionIntsByType[ion+'-H2O'] > 0:
                ionCounts['only neut'] += 1
                break

    #Now calculate statistics for immonium ions and internal fragments
    for i in range(len(PRMLadder)):
        for j in range(i+1, len(PRMLadder)):
            intFragMass = PRMLadder[j] - PRMLadder[i] + Constants.mods['H+']
            intScore = spec.getIntScore(intFragMass)
            closestInd = np.argmin(np.abs(spectrumData[:,0]-intFragMass))
            if abs(spectrumData[closestInd][0] - intFragMass) < 5*epSTD:
                ionCounts['internal'] += 1
                explainedPeakInds.add(closestInd)
                
                if i in pairedPRMInds and j in pairedPRMInds:
                    explainedPeakIndsByPairedPRMs.add(closestInd)

    PRMLadderWithEnds = [28.0313] + list(PRMLadder) + [spec.pm]
    for i in range(len(PRMLadderWithEnds)-1):
        imIonMass = PRMLadderWithEnds[i+1] - PRMLadderWithEnds[i] - Constants.mods['CO'] + Constants.mods['H+']
        intScore = spec.getIntScore(imIonMass)
        closestInd = np.argmin(np.abs(spectrumData[:,0]-imIonMass))
        if abs(spectrumData[closestInd][0] - imIonMass) < 5*epSTD:
            ionCounts['imm'] += 1
            explainedPeakInds.add(closestInd)
                    
            if i in pairedPRMInds and j in pairedPRMInds:
                explainedPeakIndsByPairedPRMs.add(closestInd)

    
    for ion in massOffsets:
        if len(massOffsets[ion]) == 0:
            massOffsets[ion] = [2*epSTD]
    
    relMassOffsets = {'b': np.array(massOffsets['b'])[1:] - np.array(massOffsets['b'])[:-1], 'y': np.array(massOffsets['y'])[1:] - np.array(massOffsets['y'])[:-1]}
    if relMassOffsets['b'].size == 0:
        relMassOffsets['b'] = [2*epSTD]
    if relMassOffsets['y'].size == 0:
        relMassOffsets['y'] = [2*epSTD]
    
    explainedInt = np.sum(spectrumData[list(explainedPeakInds), 1])/np.sum(spectrumData[:,1])
    return {'Percent Explained Intensity': explainedInt, 'Average Mass Offsets': {'b': sum(massOffsets['b'])/len(massOffsets['b']), 'y': sum(massOffsets['y'])/len(massOffsets['y'])}, 'Max Mass Offsets': {'b': max(np.abs(massOffsets['b'])), 'y': max(np.abs(massOffsets['y']))}, 'Max Rel Mass Offsets': {'b': max(np.abs(relMassOffsets['b'])), 'y': max(np.abs(relMassOffsets['y']))}, 'Fragment Ion Counts': ionCounts, 'Annotated Peaks 25': len(top25Inds & explainedPeakInds), 'Annotated Peaks 50': len(top50Inds & explainedPeakInds), 'Annotated Peaks Paired 25': len(top25Inds & explainedPeakIndsByPairedPRMs), 'Annotated Peaks Paired 50': len(top50Inds & explainedPeakIndsByPairedPRMs)}


"""
This function computes the following features for all spectra in the cluster:
- path score (of spectrum relative to peptide)
- minimum path score (of all members of cluster)
- maximum path score (of all members of cluster)
- minimal cleavage scores (score of lowest scoring PRM)
- number of PRMS with negative score (per spectrum)
- delta - abs(predicted - observed peptide mass)
- total path score (dependent on cluster size)
- # PRMs with evidence from all members of cluster
- % PRMs with evidence from all members of cluster
- # PRMs with evidence from majority of cluster
- % PRMs with evidence from majority of cluster
- # PRMs with no evidence from any member of cluster
- % PRMs with no evidence from any member of cluster

NOTE: Score is calculated using just prior and posterior from Pnet, no ambiguous edge or ppm penalties levied!
"""
def getScoreStats(specs, peptide, ambigEdges=None, ppmstd=5):

    prmStats = {'All Evidence': 0, 'Majority Evidence': 0, 'None Evidence': 0}
    nodeGen = Constants.nodeInfoGen(peptide, considerTerminalMods=True, ambigEdges=ambigEdges)
    prevNodeMass = 0
    totalMinNodeScore = 100000

    specScoreStats = {}
    for i in range(len(specs)):
        specScoreStats[i] = {'Score': 0, 'Negative': 0, 'Min Node Score': 100000}
        
    for node in nodeGen:
        numPos = 0
        totalNodeScore = 0
        
        for i, spec in enumerate(specs):
            if prevNodeMass == 0:
                specScoreStats[i]['Score'] += spec.getPriorScore(prm=0, formAA=None, lattAA=node['formAA'])

            score = spec.getNodeScore(prm=node['prm'], formAA=node['formAA'], lattAA=node['lattAA'])
            specScoreStats[i]['Score'] += score
            totalNodeScore += score
            
            if score > 0:
                numPos += 1
            if score < 0:
                specScoreStats[i]['Negative'] += 1
            if score < specScoreStats[i]['Min Node Score']:
                specScoreStats[i]['Min Node Score'] = score

        if numPos == len(specs):
            prmStats['All Evidence'] += 1
        if 2*numPos >= len(specs):
            prmStats['Majority Evidence'] += 1
        if numPos == 0:
            prmStats['None Evidence'] += 1
        if totalNodeScore < totalMinNodeScore:
            totalMinNodeScore = totalNodeScore

        prevNodeMass = node['prm']

    for i, spec in enumerate(specs):
        specScoreStats[i]['Score'] += spec.getPriorScore(prm=spec.pm, formAA=node['lattAA'], lattAA=None)

    scores = [specScoreStats[i]['Score'] for i in specScoreStats]
    totalPathScore = sum(scores)
    maxPathScore = max(scores)
    minPathScore = min(scores)
        
    return {'Total Path Score': totalPathScore, 'Total Minimum Node Score': totalMinNodeScore, 'Minimum Path Score': minPathScore, 'Maximum Path Score': maxPathScore, 'Aggregate PRM Score Statistics': prmStats, 'Spectrum Specific Score Statistics': specScoreStats}




"""
This function computes the following features:
- number forbidden node pairs
"""
def getNumForbiddenPairs(prmLadder, precMass, ppmstd=5):
    antiLadder = precMass - Constants.mods['H+'] - prmLadder

#    print 'prmLadder', prmLadder
#    print 'antiLadder', antiLadder
    epsilon = 2*ppmstd*precMass*10**-6
    res = epsilon/5

    hPRMLadder = np.round(prmLadder/res)
    hDict = {}
    for elem in hPRMLadder:
        for hMass in range(elem-5, elem+6):
            hDict[hMass] = 0
            
    numAntisymNodes = 0
    hAntiLadder = np.round(antiLadder/res)
    for hMass in hAntiLadder:
        if hMass in hDict:
            numAntisymNodes += 1

#    print 'Num Anitsym: ', numAntisymNodes/2

    return numAntisymNodes/2
    
    
