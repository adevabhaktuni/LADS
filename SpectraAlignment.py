'''
Created on Jun 27, 2011

@author: Arun
'''

import DataFile
import Constants
import DeNovoSequencer as DNS
import numpy as np
import copy
import time
from pprint import pprint

def buildSpectrumTable(precMass1, precMass2, spectra1, spectra2, resolution=0.5):
    binnedSpec1 = binSpectrum(precMass1, spectra1, resolution)
    binnedSpec2 = binSpectrum(precMass2, spectra2, resolution)
    
    specTable = np.outer(binnedSpec2, binnedSpec1)
    return specTable
    
def binSpectrum(precMass, spectra, resolution=0.5):
    numBins = np.round(precMass / resolution) + 1
    binnedSpec = np.zeros(numBins)
    for pair in spectra:
        mass = pair[0]
        intMass = np.round(mass / resolution)
        binnedSpec[intMass] = 1
    
    return binnedSpec

def parallelogramMethod(specTable, resolution=0.5):
    intH2O = np.round(Constants.mods['H2O'] / resolution)
    intH = np.round(Constants.mods['H+'] / resolution)
    
    delta = specTable.shape[0] - specTable.shape[1]
    mass = min(specTable.shape[0] - 1 - intH2O - intH, specTable.shape[1] - 1 - intH2O - intH)
    
    count = 0
    offset = mass + delta + intH2O + 2 * intH
    density = []
    for i in xrange(specTable.shape[1]):
        try:            
            if specTable[i][i] > 0:
                count += 1
            if specTable[delta + offset - i][offset - i]:
                count += 1
        except IndexError:
            continue       
        density.extend([float(count)])
        
    return density

def hashSpectrum(spectrum, specID=1, prevTable=None, epsilon=0.04):
    if prevTable:
        hashTable = copy.deepcopy(prevTable)
    else:
        hashTable = {}
        
    for mass in spectrum[:, 0]:
        intMass = np.round(mass / epsilon)
        try:
            hashTable[intMass] = hashTable[intMass] + [(specID, mass)]
        except KeyError:
            hashTable[intMass] = [(specID, mass)]
        
    return hashTable   

def parallelogramMethodCont(precMass1, precMass2, spectra1, spectra2, epsilon=0.04):
    delta = precMass1 - precMass2
    baseMass = min(precMass1, precMass2) - Constants.mods['H2O'] - Constants.mods['H+']
    
    offset = baseMass + delta + Constants.mods['H2O'] + 2 * Constants.mods['H+']    
    intOffset = np.round(offset / epsilon)
    
    hashTable1 = hashSpectrum(spectra1, specID=1)
    
    hashDiag1 = hashSpectrum(spectra2, specID=2, prevTable=copy.deepcopy(hashTable1))
    for key in hashDiag1.keys():
        if len(hashDiag1[key]) <= 1:
            del hashDiag1[key]
    
    spectra2[:, 0] = spectra2[:, 0] + delta
    hashDiag2 = hashSpectrum(spectra2, specID=2, prevTable=copy.deepcopy(hashTable1))
    for key in hashDiag2.keys():
        if len(hashDiag2[key]) <= 1:
            del hashDiag2[key]
            
    zeroDiagMasses = np.sort(np.array(hashDiag1.keys()))
    deltaDiagMasses = np.sort(np.array(hashDiag2.keys()))
    
    density = []
    for i in range(zeroDiagMasses.size):
        count = i + 1
        recipMass = intOffset - zeroDiagMasses[i]
        count += (np.where(deltaDiagMasses >= recipMass))[0].size
        density.extend([count / zeroDiagMasses[i]])
    
    return density
    
def mergeSpectra(specs, epsilon=0.02):
    mergedDict = {}
    numSpecs = len(specs)
    for i in range(numSpecs):
        hMasses = np.round(specs[i][:,0]/epsilon)
        for j, hMass in enumerate(hMasses):
            if hMass in mergedDict:
                mergedDict[hMass]['pairs'] += [specs[i][j]]
                mergedDict[hMass]['present'][i] = 1
            else:
                mergedDict[hMass] = {'pairs': [specs[i][j]], 'present': np.zeros(numSpecs)}
                mergedDict[hMass]['present'][i] = 1
    
    mergedSpec = np.zeros((len(mergedDict), 3))
    for i, hMass in enumerate(sorted(mergedDict)):
        mergedSpec[i][:2] = np.average(np.array(mergedDict[hMass]['pairs']), 0)
        mergedSpec[i][2] = np.sum(mergedDict[hMass]['present'])
    
    #scale intensities by number of spectra in which they are found (basically the intensity of a merged peak is the sum of the intensities of all peaks merged into it)
    mergedSpec[:,1] = mergedSpec[:,1] * mergedSpec[:,2]
    return mergedSpec

def getSharedPeaksRatios(lightMergedSpec, heavyMergedSpec, lightPrecMass, NMod=0, CMod=0, epsilon=0.02):
    NTermTable, CTermTable = getNandCIons(lightMergedSpec, heavyMergedSpec, Nmod=NMod, Cmod=CMod, epsilon=epsilon)
    NCrossTable, CCrossTable = getCrossPairedIons(lightMergedSpec, heavyMergedSpec, lightPrecMass, Nmod=NMod, Cmod=CMod, epsilon=epsilon)

    NTermIonDict = prepIonTableForAddition(NTermTable, ['b', 'b'])
    CTermIonDict = prepIonTableForAddition(CTermTable, ['y', 'y'])
    NCrossIonDict = prepIonTableForAddition(NCrossTable, ['y', 'b'])
    CCrossIonDict = prepIonTableForAddition(CCrossTable, ['b', 'y'])
            
    allPairedIonsDict = addDicts(reverseDict(addDicts(NTermIonDict, CCrossIonDict)), reverseDict(addDicts(NCrossIonDict, CTermIonDict)))
    symLightInds = set()
    symHeavyInds = set()
    totalLightInds = set()
    totalHeavyInds = set()
    doubleSymLightInds = set()
    doubleSymHeavyInds = set()
            
    for heavyIons in allPairedIonsDict:
        if len(heavyIons) == 2 and len(allPairedIonsDict[heavyIons]) == 2:
            for ion in heavyIons:
                doubleSymHeavyInds.add(ion[0])
            for ion in allPairedIonsDict[heavyIons]:
                doubleSymLightInds.add(ion[0])
        if len(heavyIons) == 2 or len(allPairedIonsDict[heavyIons]) ==2:
            for ion in heavyIons:
                symHeavyInds.add(ion[0])
            for ion in allPairedIonsDict[heavyIons]:
                symLightInds.add(ion[0])
                        
        for ion in heavyIons:
            totalHeavyInds.add(ion[0])
        for ion in allPairedIonsDict[heavyIons]:
            totalLightInds.add(ion[0])

    totalNumPeaks = float(lightMergedSpec.shape[0] + heavyMergedSpec.shape[0])
    totalSharedPeaksRatio = (len(totalLightInds) + len(totalHeavyInds))/totalNumPeaks
    singleSymSharedPeaksRatio = (len(symLightInds) + len(symHeavyInds))/totalNumPeaks
    doubleSymSharedPeaksRatio = (len(doubleSymLightInds) + len(doubleSymHeavyInds))/totalNumPeaks
    
    return totalSharedPeaksRatio, singleSymSharedPeaksRatio, doubleSymSharedPeaksRatio


def getSpectraPairInfoForSVMClassification(lightMergedSpec, heavyMergedSpec, lightPrecMass, NMod=0, CMod=0, epsilon=0.02):
    NTermTable, CTermTable = getNandCIons(lightMergedSpec, heavyMergedSpec, Nmod=NMod, Cmod=CMod, epsilon=epsilon)
    NCrossTable, CCrossTable = getCrossPairedIons(lightMergedSpec, heavyMergedSpec, lightPrecMass, Nmod=NMod, Cmod=CMod, epsilon=epsilon)

    NTermIonDict = prepIonTableForAddition(NTermTable, ['b', 'b'])
    CTermIonDict = prepIonTableForAddition(CTermTable, ['y', 'y'])
    NCrossIonDict = prepIonTableForAddition(NCrossTable, ['y', 'b'])
    CCrossIonDict = prepIonTableForAddition(CCrossTable, ['b', 'y'])
            
    allPairedIonsDict = addDicts(reverseDict(addDicts(NTermIonDict, CCrossIonDict)), reverseDict(addDicts(NCrossIonDict, CTermIonDict)))
    symLightInds = set()
    symHeavyInds = set()
    totalLightInds = set()
    totalHeavyInds = set()
    doubleSymLightInds = set()
    doubleSymHeavyInds = set()
            
    for heavyIons in allPairedIonsDict:
        if len(heavyIons) == 2 and len(allPairedIonsDict[heavyIons]) == 2:
            for ion in heavyIons:
                doubleSymHeavyInds.add(ion[0])
            for ion in allPairedIonsDict[heavyIons]:
                doubleSymLightInds.add(ion[0])
        if len(heavyIons) == 2 or len(allPairedIonsDict[heavyIons]) ==2:
            for ion in heavyIons:
                symHeavyInds.add(ion[0])
            for ion in allPairedIonsDict[heavyIons]:
                symLightInds.add(ion[0])
                        
        for ion in heavyIons:
            totalHeavyInds.add(ion[0])
        for ion in allPairedIonsDict[heavyIons]:
            totalLightInds.add(ion[0])

    totalNumPeaks = float(lightMergedSpec.shape[0] + heavyMergedSpec.shape[0])
    totalSharedPeaksRatio = (len(totalLightInds) + len(totalHeavyInds))/totalNumPeaks
    singleSymSharedPeaksRatio = (len(symLightInds) + len(symHeavyInds))/totalNumPeaks
    doubleSymSharedPeaksRatio = (len(doubleSymLightInds) + len(doubleSymHeavyInds))/totalNumPeaks

    lightTIC = np.sum(lightMergedSpec[:,1])
    heavyTIC = np.sum(heavyMergedSpec[:,1])
    sharedLightTIC = np.sum(lightMergedSpec[list(totalLightInds), 1])
    sharedHeavyTIC = np.sum(heavyMergedSpec[list(totalHeavyInds), 1])
    singleSymLightTIC = np.sum(lightMergedSpec[list(symLightInds), 1])
    singleSymHeavyTIC = np.sum(heavyMergedSpec[list(symHeavyInds), 1])
    doubleSymLightTIC = np.sum(lightMergedSpec[list(doubleSymLightInds), 1])
    doubleSymHeavyTIC = np.sum(heavyMergedSpec[list(doubleSymHeavyInds), 1])

    totalTICRatio = (sharedLightTIC + sharedHeavyTIC)/(lightTIC + heavyTIC)
    singleSymTICRatio = (singleSymLightTIC + singleSymHeavyTIC)/(lightTIC + heavyTIC)
    doubleSymTICRatio = (doubleSymLightTIC + doubleSymHeavyTIC)/(lightTIC + heavyTIC)

#    return {1: totalSharedPeaksRatio, 2: singleSymSharedPeaksRatio, 3: doubleSymSharedPeaksRatio, 4: lightPrecMass}
    return {1: totalSharedPeaksRatio, 2: singleSymSharedPeaksRatio, 3: doubleSymSharedPeaksRatio, 4: totalTICRatio, 5: singleSymTICRatio, 6: doubleSymTICRatio, 7: lightPrecMass}
    
def getPairedIons(hashTable, spectra, delta, epsilon=0.02):
    table = copy.deepcopy(hashTable)
    for i in range(spectra.shape[0]):
        m = spectra[i, 0] - delta
        hMass = np.round(m / epsilon)
        try:
            table[hMass] += [(i, spectra[i, 0])]
        except KeyError:
            try:
                table[hMass - 1] += [(i, spectra[i, 0])]
            except KeyError:
                try:
                    table[hMass + 1] += [(i, spectra[i, 0])]
                except KeyError:
                    pass
            
    for mass in table.keys():
        l = len(table[mass])
        if l < 2:
            del table[mass]
        if l > 2:
            table[mass] = table[mass][:2]
    
    return table

def reverseDict(dictToReverse):
    return dict((dictToReverse[i], i) for i in dictToReverse)

#both keys and values must be tuples                                                                                                                                                                                 
def addDicts(dict1, dict2):
    combDict = copy.deepcopy(dict1)
    for key in dict2:
        if key in dict1:
            combDict[key] = dict1[key] + dict2[key]
        else:
            combDict[key] = dict2[key]

    return combDict

def prepIonTableForAddition(ionTable, pairType):
    Inds = []
    for key in ionTable:
        Inds += [zip(*ionTable[key])[0]]
    ionDict = dict((((IndPair[0], pairType[0]),),((IndPair[1], pairType[1]),)) for IndPair in Inds)
    return ionDict

  

def getSymmetricIons(spec, PM, epsilon=0.02):
    hashTable = {}
    symSpecMasses = PM + Constants.mods['H2O'] + 2*Constants.mods['H+'] - spec[:,0]
    for i in range(symSpecMasses.size):
        key = np.round(symSpecMasses[i]/epsilon)
        hashTable[key-1] = True
        hashTable[key] = True
        hashTable[key+1] = True
        
    symMasses = {}
    for i in range(spec.shape[0]):
        if np.round(spec[i][0]/epsilon) in hashTable:
            symMasses[i] = spec[i][0]
            
    return symMasses

def getCrossPairedIons(lightSpec, heavySpec, lightPrecMass, Nmod, Cmod, epsilon=0.02):
    hashTable = {}
    for i in range(lightSpec.shape[0]):
        mass = lightSpec[i, 0]
        key = np.round(mass / epsilon)
        hashTable[key] = [(i, mass)]
        
    revHeavySpec = copy.copy(heavySpec)
    revHeavySpec[:,0] = lightPrecMass + Constants.mods['H+'] - heavySpec[:,0]
    NCrossTable = getPairedIons(hashTable, revHeavySpec, delta=-Nmod, epsilon=epsilon)
    CCrossTable = getPairedIons(hashTable, revHeavySpec, delta=-Cmod, epsilon=epsilon)
    return NCrossTable, CCrossTable                                            

def getNandCIons(spectra1, spectra2, Nmod, Cmod, epsilon=0.02):
    hashTable = {}
    for i in range(spectra1.shape[0]):
        mass = spectra1[i, 0]
        key = np.round(mass / epsilon)
        hashTable[key] = [(i, mass)]
    
    # Get NTerm Ions
    NTermTable = getPairedIons(hashTable, spectra2, delta=Nmod, epsilon=epsilon)
    # Get NTerm Ions
    CTermTable = getPairedIons(hashTable, spectra2, delta=Cmod, epsilon=epsilon)
    return NTermTable, CTermTable
    
def getSharedPeaksRatio(spectra1, spectra2, N, C):
    return float(len(N.keys()) + len(C.keys())) / (spectra1.shape[0] + spectra2.shape[0])
           
    
def getHashedMassesForAlignment(precMass, spec, epsilon=0.08):
    hashSpec = hashSpectrum(spec, epsilon=epsilon)
    hshMasses = np.sort(np.array(hashSpec.keys()))
    intPrecMass = np.round(precMass / epsilon)
    hshMasses = np.insert(hshMasses, 0, 0)
    hshMasses = np.append(hshMasses, intPrecMass)
    
    return hashSpec, hshMasses

def diagonalHash(masses1, masses2):
    diagHash = {}
    for i in range(masses1.size + masses2.size - 1):
        for j in range(i + 1):
            if j < masses1.size and i - j < masses2.size:
                m = masses1[j]
                n = masses2[i - j]
                try:
                    diagHash[m - n] = diagHash[m - n] + [(j, i - j)]
                except KeyError:
                    diagHash[m - n] = [(j, i - j)]
    
    return diagHash
            
def getMassesForDiagHash(precMass, spec):
    masses = spec[:, 0]
    masses = np.insert(masses, 0, 0)
    masses = np.append(masses, precMass)
    return masses

def diagonalHash2(masses1, masses2, epsilon=0.02):
    diagHash = {}
    for i in range(masses1.size + masses2.size - 1):
        for j in range(i + 1):
            if j < masses1.size and i - j < masses2.size:
                m = masses1[j]
                n = masses2[i - j]
                try: 
                    ind = np.floor((m - n) / epsilon)
                    diagHash[ind] = diagHash[ind] + [(j, i - j)]
                except KeyError:
                    try:
                        ind = np.ceil((m - n) / epsilon)
                        diagHash[ind] = diagHash[ind] + [(j, i - j)]
                    except KeyError:
                        ind = np.round((m - n) / epsilon)
                        diagHash[ind] = [(j, i - j)]
    
    return diagHash

def getD0Matrix(masses1, masses2, diagHash):
    D0 = np.ones((masses1.size, masses2.size))
    D0[0][0] = 0
    for key in diagHash.keys():
        count = 0
        inds = diagHash[key]
        for pair in inds:
            D0[pair[0], pair[1]] = count
            count += 1
            
    return D0

def getMMatrix(D):
    M = np.zeros(D.shape)
    for i in range(D.shape[0] + D.shape[1] - 1):
        for j in range(1, i):
            if j < D.shape[0] and i - j < D.shape[1]:
                Dval = D[j, i - j]
                Mi = M[j - 1, i - j]
                Mj = M[j, i - j - 1]
                M[j][i - j] = max(Dval, Mi, Mj)
    
    return M

def getDkMatrix(masses1, masses2, Mprev, diagHash, epsilon=0.02):
    Dk = np.zeros(Mprev.shape)
    for i in range(Dk.shape[0] + Dk.shape[1] - 1):
        for j in range(1, i):
            if j < Dk.shape[0] and i - j < Dk.shape[1]:
                m = masses1[j]
                n = masses2[i - j]
                coDiag = hashMasses(m, n, (j, i - j), diagHash, epsilon)
                if len(coDiag) > 1:
                    temp = coDiag[0]
                    for elem in coDiag:
                        if (elem[0] == j and elem[1] == i - j):
                            break
                        else:
                            temp = elem   
                                               
                    Dk[j, i - j] = max(Dk[temp[0], temp[1]], Mprev[j - 1, i - j - 1]) + 1
                else:
                    Dk[j, i - j] = Mprev[j - 1, i - j - 1] + 1
    
    return Dk


def MaxPath(k, Dstack, masses1, masses2, d, sub=None, N=None, epsilon=0.02):
    path = []
    diagHash = copy.deepcopy(d)
    
    if not sub:
        sub = (Dstack.shape[0], Dstack.shape[1])
    k = min(k, Dstack.shape[2] - 1)
    
    Dcurrent = Dstack[0:sub[0], 0:sub[1], k] 
    
    if not N:
        N = np.amax(Dcurrent)
        
    inds = np.where(Dcurrent == N)
    curInd = (0, 0)
    for i in range(inds[0].size):
        if (inds[0][i] + inds[1][i]) > (curInd[0] + curInd[1]):
            curInd = (inds[0][i], inds[1][i])
    
    m = masses1[curInd[0]]
    n = masses2[curInd[1]]
    curDiag = hashMasses(m, n, (curInd[0], curInd[1]), diagHash, epsilon)
    
    diagInd = curDiag[-1]
    while diagInd != curInd:
        curDiag.pop()
        diagInd = curDiag[-1]
      
    while N > 0:
        path.extend([(masses1[curInd[0]], masses2[curInd[1]])])
        curDiag.pop()
        if curDiag:
            newInd = curDiag[-1]
            if Dcurrent[newInd[0], newInd[1]] == N - 1:
                N -= 1
                curInd = newInd
            else:
                path.extend(MaxPath(k - 1, Dstack, masses1, masses2, diagHash, sub=(curInd[0], curInd[1]), N=N - 1, epsilon=epsilon))
                break 
        else:
            path.extend(MaxPath(k - 1, Dstack, masses1, masses2, diagHash, sub=(curInd[0], curInd[1]), N=N - 1, epsilon=epsilon))
            break 
    
    return path

def hashMasses(mass1, mass2, inds, d, epsilon=0.02):
    try:
        coDiag = d[np.floor((mass1 - mass2) / epsilon)]
        if coDiag[0] <= inds:
            return coDiag
    except KeyError:
        pass
    
    return d[np.ceil((mass1 - mass2) / epsilon)]
                    
def getKMaxPath(k, masses1, masses2, epsilon=0.02):    
    d = diagonalHash2(masses1, masses2, epsilon)

    Dstack = np.zeros((masses1.size, masses2.size, k + 1))
    Dstack[:, :, 0] = getD0Matrix(masses1, masses2, d)
    
    for i in range(1, k + 1):
        M = getMMatrix(Dstack[:, :, i - 1])
        Dstack[:, :, i] = getDkMatrix(masses1, masses2, M, d, epsilon=epsilon)
    
    path = MaxPath(k, Dstack, masses1, masses2, d, epsilon=epsilon)
    return path[::-1]

def getDiagMasses(diag, masses1, masses2):
    diagMasses = []
    for pair in diag:
        diagMasses += [(masses1[pair[0]], masses2[pair[1]])]
    
    return diagMasses

def getSymmetrizedSpectrum(masses, precMH, epsilon=0.04):
    revMasses = precMH - (masses - Constants.mods['H+'])
    return combineSpectraCompleteLinkage(masses, revMasses, epsilon)
    
def combineSpectraCompleteLinkage(masses1, masses2, epsilon=0.04):
    combMasses = np.append(masses1, masses2)
    combMasses = np.sort(combMasses)
    clusters = DNS.getClustersCompleteLinkage(combMasses, epsilon)
    combMasses = []
    for cluster in clusters:
        if len(cluster) > 1:
            combMasses += [sum(cluster) / len(cluster)]
        else:
            combMasses += cluster
    
    return np.sort(np.array(combMasses))
    
if __name__ == '__main__':
    dirPath = 'C:\\Users\\Arun\\Pythonprojects\\DeNovoSequencing\\LF2_short_HCD+CID_ath001862_244\\'
    
    ppm = 5
    heavyPath = dirPath + '244.3611.3611.1.dta'
    lightPath = dirPath + '244.3619.3619.1.dta'
    heavyPairs = DataFile.getMassIntPairs(heavyPath)
    lightPairs = DataFile.getMassIntPairs(lightPath)
    heavyPrecMass, heavyCharge = DataFile.getPrecMassAndCharge(heavyPath) 
    lightPrecMass, lightCharge = DataFile.getPrecMassAndCharge(lightPath)   
    
    print ppm * 10 ** -6 * heavyPrecMass
    print getSharedPeaksRatio(lightPairs, heavyPairs, Nmod=0, Cmod=Constants.mods['*'], epsilon=ppm * heavyPrecMass * 10 ** -6)
    
     
    """
    tPath = dirPath + '244.0855.0855.1.dta'
    tMass = DataFile.getPrecMassAndCharge(tPath)[0] 
    tPairs = DataFile.getMassIntPairs(tPath)
    tIons = tPairs[:,0]
    tIons = np.insert(tIons, 0, 0)
    tIons = np.append(tIons, tMass)
    tIons = getSymmetrizedSpectrum(tIons, tMass)
    print tIons

    NTermTable, CTermTable = getNandCIons(lightPrecMass, lightPairs, heavyPrecMass, heavyPairs)
    
    NIons = zip(*NTermTable.values())[0]
    NIons = np.array(zip(*NIons)[1])
    NIons = np.append(NIons, [0, lightPrecMass])
    #NIons.sort()
    

    CIons = zip(*CTermTable.values())[0]
    CIons = np.array(zip(*CIons)[1])
    CIons = np.append(CIons, [0, lightPrecMass])
    #CIons.sort()
    
    refIons = combineSpectraCompleteLinkage(NIons, lightPrecMass - (CIons - Constants.mods['H+']))
    print getKMaxPath(1, refIons, tIons)
    """
    """
    print unModSeq
    bIons = Constants.getBIons(unModSeq)
    bIons = np.array([0] + list(bIons) + [Constants.getMH(unModSeq)])
    yIons = Constants.getYIons(unModSeq)
    yIons = np.array([0] + list(yIons))
    
    print bIons
    print yIons
    print 'Max Path from bIons: '
    print getKMaxPath(1, bIons, tIons)
    """
    
    #print getSharedPeaksRatio(lightPrecMass, lightPairs, heavyPrecMass, heavyPairs)
    #path = getKMaxPath(1, lightPrecMass, lightPairs, heavyPrecMass, heavyPairs)
    #print path
    #print len(path)
    """
    t=[]
    for key in hashDiags[0].keys():
        t+= [hashDiags[0][key][0][1]]
    print 'length zero diag: ', len(t)
    print sorted(t)
    
    t=[]
    h=[]
    for pair in zeroDiag:
        try:
            t += [masses1[pair[0]]]
            h += [masses2[pair[1]]]
        except KeyError:
            continue
    print 'length zero diag: ', len(t)
    print sorted(t)
    print sorted(h)
    
    t=[]
    for key in hashDiags[1].keys():
        t+= [hashDiags[1][key][0][1]]
    print 'length heavy diag: ', len(t)
    print sorted(t)
    
    t=[]
    h=[]
    for pair in heavyDiag:
        try:
            t += [masses1[pair[0]]]
            h += [masses2[pair[1]]]
        except KeyError:
            continue
    print 'length heavy diag: ', len(t)
    print sorted(t)
    print sorted(h)
    #D0 = getD0Matrix(masses1[1], masses2[1], d)
    #print np.amax(D0)
    #print getAlignmentRatios("C:\\Users\\Arun\\Dropbox\\SpectraCorrelation\\SpectraCorrelation\\spectrafnames.txt")
    #plt.hist(counts, bins=100)
    #plt.show()
    """
