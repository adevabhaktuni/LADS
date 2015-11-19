'''
Created on Jun 6, 2011

@author: Arun
'''
import numpy as np
from collections import deque
import copy

#Dict of amino acid masses and formulas, first mass is monoisotopic mass, second mass is average mass
#obtained from http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html

aminoacids = {
              
'A': ('Ala', 'C3H5ON', 71.03711, 71.0788),
'R': ('Arg', 'C6H12ON4', 156.10111, 156.1875),
'N': ('Asn', 'C4H6O2N2', 114.04293, 114.1038),
'D': ('Asp', 'C4H5O3N', 115.02694, 115.0886),
'C': ('Cys', 'C3H5ONS', 103.00919, 103.1388),
'E': ('Glu', 'C5H7O3N', 129.04259, 129.1155),
'Q': ('Gln', 'C5H8O2N2', 128.05858, 128.1307),
'G': ('Gly', 'C2H3ON', 57.02146, 57.0519),
'H': ('His', 'C6H7ON3', 137.05891, 137.1411),
'I': ('Ile', 'C6H11ON', 113.08406, 113.1594),
#'L': ('Leu', 'C6H11ON', 113.08406, 113.1594),
'K': ('Lys', 'C6H12ON2', 128.09496, 128.1741),
'M': ('Met', 'C5H9ONS', 131.04049, 131.1926),
'F': ('Phe', 'C9H9ON', 147.06841, 147.1766),
'P': ('Pro', 'C5H7ON', 97.05276, 97.1167),
'S': ('Ser', 'C3H5O2N', 87.03203, 87.0782),
'T': ('Thr', 'C4H7O2N', 101.04768, 101.1051),
'W': ('Trp', 'C11H10ON2', 186.07931, 186.2132),
'Y': ('Tyr', 'C9H9O2N', 163.06333, 163.1760),
'V': ('Val', 'C5H9ON', 99.06841, 99.1326)                         
}

origAAs = copy.copy(aminoacids)

mods = {  
'H+': 1.00727638,
'Hyd': 1.0078250321,
'H2O': 18.0105647,
'OH': 17.0032883,
'NH3': 17.0265491,
'NH2': 16.0187240679,
'CO': 27.9949146
}


NTermMods = {'static': 0}
CTermMods = {'static': 0}

diffModSymbols = ['#', '*', '^', '$', '?', '~', '<', '>', '&', '!', '@', '[', ']', '(', ')']

def hashAAsWithPPM(aas, precMass, ppmSTD, ppmSysError=0, ppmStep=0.5, maxNoDevs=3):
    ppmStepEp = ppmStep * precMass * 10**-6
    maxIntPPM = np.ceil((ppmSTD*maxNoDevs+ppmSysError)/ppmStep)
    minIntPPM = np.floor((-ppmSTD*maxNoDevs+ppmSysError)/ppmStep)

    hashed = {}
    for key in aas:
        newKey = np.round(aas[key][2]/ppmStepEp)
        for intPPM in range(minIntPPM, maxIntPPM+1):
            try:
                hashed[intPPM+newKey]['seqs'][key] = intPPM
                if intPPM < hashed[intPPM+newKey]['min']:
                    hashed[intPPM+newKey]['min'] = intPPM
            except KeyError:
                hashed[intPPM+newKey] = {'min': intPPM, 'seqs': {key: intPPM}}

    return hashed
                
def hashAAsEpsilonRange(aas, epStep=0.0005, maxEp=0.08):
    maxIntEp = int(np.ceil(maxEp/epStep))
    minIntEp = int(np.floor(-maxEp/epStep))
    
    hashed = {}
    for key in aas:
        newKey = np.round(aas[key][2]/epStep)
        for intEp in range(minIntEp, maxIntEp+1):
            try:
                hashed[intEp+newKey]['seqs'][key] = epStep*intEp
                if np.abs(epStep*intEp) < np.abs(hashed[intEp+newKey]['min']):
                    hashed[intEp+newKey]['min'] = epStep*intEp
            except KeyError:
                hashed[intEp+newKey] = {'min': epStep*intEp, 'seqs': {key: epStep*intEp}}

    return hashed

def hashUnimodModAAsEpsilonRange(unimodPeptDict, epStep = 0.0005, maxEp=0.08):
    maxIntEp = int(np.ceil(maxEp/epStep))
    minIntEp = int(np.floor(-maxEp/epStep))
    
    hashed = {}
    modList = {}
    for terminus in unimodPeptDict:
        hashed[terminus] = {}
#        modList[terminus] = zip(*unimodPeptDict[terminus])[0]
        for i, peptide in enumerate(unimodPeptDict[terminus]):
            newKey = np.round(peptide[1]/epStep)
#            print peptide, epStep, newKey
            for intEp in range(minIntEp, maxIntEp+1):
                try:
                    hashed[terminus][intEp+newKey]['seqs'][peptide[0]] = epStep*intEp
                    if np.abs(epStep*intEp) < np.abs(hashed[intEp+newKey]['min']):
                        hashed[terminus][intEp+newKey]['min'] = epStep*intEp
                except KeyError:
                    hashed[terminus][intEp+newKey] = {'min': epStep*intEp, 'seqs': {peptide[0]: epStep*intEp}}

    return hashed

def hashAAs(aas, epsilon):
    hashed = {}
    for key in aas.keys():
        try:
            newKey = np.round(aas[key][2] / epsilon)
        except TypeError:
            newKey = np.round(aas[key] / epsilon)
            
        for hKey in [newKey - 1, newKey, newKey + 1]:
            try:
                hashed[hKey] += [key]
            except KeyError:
                hashed[hKey] = [key]
                      
    return hashed

def createTermModHashAAs(epsilon=0.02, N=copy.deepcopy(NTermMods), C=copy.deepcopy(CTermMods)):
    try:
        del N['static']
    except KeyError:
        pass
    
    try:
        del C['static']
    except KeyError:
        pass
    
    termModHash = {}
    termModHash['NTerm'] = hashAAs(N, epsilon)
    termModHash['CTerm'] = hashAAs(C, epsilon)
    return termModHash

def addStaticMod(datum, args):
    try:
        AAData = list(aminoacids[args[0]])
        AAData[2] += float(args[1])
        aminoacids[args[0]] = tuple(AAData)
    except KeyError:
        if args[0] == 'N-Term':
            NTermMods['static'] = float(args[1])
        elif args[0] == 'C-Term':
            CTermMods['static'] = float(args[1])
        else:
            raise KeyError('%s not a valid amino acid. Must add AAs to dictionary before modifying them.' % (args[0],))
    
#    mods[datum] = float(args[1])
    return {(datum, args[0]): args[1]}
    
def addDiffMod(datum, args):
    print datum, args
    try:
        symb = args[3]
    except IndexError:
        symb = None
        for modSymbol in diffModSymbols:
            if modSymbol not in mods.keys():
                symb = modSymbol
                break
        if not symb:    
            symb = '+%i' % (int(args[1]),)
    
    mods[symb] = float(args[1])
    try:
        if bool(int(args[2])):
            AAdata = list(origAAs[args[0]])
        else:
            AAdata = list(aminoacids[args[0]])
            
        AAdata[2] += float(args[1])
        print AAdata
        aminoacids[args[0] + symb] = tuple(AAdata)
    except KeyError:
        if args[0] == 'N-Term':
            NTermMods[symb] = float(args[1]) + (0 if bool(int(args[2])) else NTermMods['static'])
        elif args[0] == 'C-Term':
            CTermMods[symb] = float(args[1]) + (0 if bool(int(args[2])) else NTermMods['static'])
        else:
            raise KeyError('%s not a valid amino acid. Must add AAs to dictionary before modifying them.' % (args[0],))

    return {symb: (datum, args[0], args[1])}
            
def getAA(mass, tolerance=.1):
    residual = tolerance
    aa = None
    for aaKey in aminoacids.keys():
        newResid = np.abs(aminoacids[aaKey][2] - mass)
        if newResid < residual:
            residual = newResid
            aa = aaKey
    
    return aa

def addAA(datum, args):
    try:
        aminoacids[args[0]]
        raise ValueError("%s already exists in amino acid dictionary. Choose a character different than the following: %s" % (args[0], aminoacids.keys()))
    except KeyError:
        aminoacids[args[0]] = tuple(args[1], args[2], float(args[3]), float(args[4]))
        return {args[0]:  tuple(args[1:])}

# each entry is of the form ((modAA, ambigEdge), modMass)
def getModList(unimodDict):
    modSet = set()
    for modMass in unimodDict:
        modSet.add(('X' if unimodDict[modMass][1] != 'Unmod' else unimodDict[modMass][0], ((0.0, modMass),) if unimodDict[modMass][1] != 'Unmod' else ()))

    return modSet

# Entries in the queue are of the form [(seq, ambigEdges), mass]
def blindModPeptDFS(maxPolyPepMass, addModList, maxPolyPepLength = 3, seedList=None):
    queue = deque()
    for item in seedList:
        if item[0] in aminoacids:
            mass = aminoacids[item[0]][2]
        else:
            mass = item[1][0][1]
        queue.extend([(item, mass)])

    while queue:
        pept, peptMass = queue.pop()
        yield (pept, peptMass)
        for modEntry in addModList:
            if modEntry[0] in aminoacids:
                mass = aminoacids[modEntry[0]][2]
            else:
                mass = modEntry[1][0][1]
            if peptMass + mass <= maxPolyPepMass and len(pept[0]) + 1 <= maxPolyPepLength:
                queue.extend([((pept[0] + modEntry[0], pept[1] + modEntry[1]), peptMass + mass)])

def peptDFS(mass, aminoacids=aminoacids):
    queue = deque([(aaKey, aminoacids[aaKey][2]) for aaKey in aminoacids.keys()])
    while queue:
        pept, peptMass = queue.pop()
        yield (pept, peptMass)
        for aaKey in aminoacids.keys():
            aaMass = aminoacids[aaKey][2]
            if (peptMass + aaMass <= mass):
                queue.extend([(pept + aaKey, peptMass + aaMass)])
    
#Depth-first search for candidate peptides of mass M
def getCandidatePeptides(mass, tolerance=0.02):
    peptGen = peptDFS(mass + tolerance)
    candidateList = []
    for pair in peptGen:
        peptide = pair[0]
        peptMass = pair[1]
        if np.abs(peptMass - mass) < tolerance:
            candidateList.extend([peptide])
    
    return candidateList

def getPeptidesLessThan(mass):
    peptGen = peptDFS(mass)
    candidateList = []
    for pair in peptGen:
        peptide = pair[0]
        peptMass = pair[1]
        if peptMass < mass:
            candidateList.extend([peptide])
    
    return candidateList

def parseModifications(modDict, aminoacids=aminoacids, epsilon=0.02, minModAAMass = 43, precision=8):
    intEp = 0.0025
    hashedAAs = hashAAsEpsilonRange(aminoacids, intEp, epsilon)

    modFormattedAADict = {}
    for aa in aminoacids:
        modFormattedAADict[aminoacids[aa][2]] = (aa, 'Unmod')
        
    modAAsDict = {'N-term': {}, 'Anywhere': {}, 'C-term': {}}
    modAAsDict['Anywhere'] = copy.deepcopy(modFormattedAADict)

    
    for mod in modDict:
        for location in modDict[mod]['locations']:
            aa = location[0] if not location[0] == 'L' else 'I'
            if location[1] == 'Anywhere' and location[0] != 'N-term' and location[0] != 'C-term':
                modMass = np.round(aminoacids[aa][2] + modDict[mod]['mass'], decimals=precision)
                if np.round(modMass/intEp) not in hashedAAs and modMass > minModAAMass:
                    modAAsDict['Anywhere'][modMass] = (aa, mod)
                    
    # Now do terminal mods
    acetylTermModMasses = []
    for mod in modDict:
        for location in modDict[mod]['locations']:
            if not (location[1] == 'Anywhere' and location[0] != 'N-term' and location[0] != 'C-term'):
                aa = location[0] if not location[0] == 'L' else 'I'
                if mod == 'Acetyl':
                    print aa
                for modMass in modAAsDict['Anywhere']:
                    # Don't consider double modified terminal modifications
                    if 'Unmod' == modAAsDict['Anywhere'][modMass][1] and (modAAsDict['Anywhere'][modMass][0] == aa or aa == 'N-term' or aa == 'C-term'):
                        term = location[1]
                        if aa == 'N-term' or aa == 'C-term':
                            term = aa
                        termModMass = np.round(modMass + modDict[mod]['mass'], decimals=precision)
                        if np.round(termModMass/intEp) not in hashedAAs and termModMass > minModAAMass:
                            if not(termModMass in modAAsDict[term] and 'Unmod' in modAAsDict[term][termModMass][1]):
                                modAAsDict[term][termModMass] = (modAAsDict['Anywhere'][modMass][0], modAAsDict['Anywhere'][modMass][1] + ' ' + term + ' ' + mod)
                            if mod == 'Acetyl' and 'Unmod' in modAAsDict['Anywhere'][modMass][1]:
                                acetylTermModMasses += [termModMass]


    return modAAsDict

def getUnimodPeptDict(maxPolyPepMass, unimodDict, maxPolyPepLength=3):
    unimodPeptDict = {}
    
    anywhereSeedSet = getModList(unimodDict['Anywhere']) 
    ntermSeedSet = getModList(unimodDict['N-term']) | anywhereSeedSet
    ctermSeedSet = getModList(unimodDict['C-term']) | anywhereSeedSet

    print 'Anywhere'
    unimodPeptDict['Anywhere'] = []
    for peptide in blindModPeptDFS(maxPolyPepMass, list(anywhereSeedSet), maxPolyPepLength = 3, seedList=list(anywhereSeedSet)):
        unimodPeptDict['Anywhere'] += [peptide]

    print 'C-term'
    unimodPeptDict['C-term'] = []
    for peptide in blindModPeptDFS(maxPolyPepMass, list(anywhereSeedSet), maxPolyPepLength = 3, seedList=list(ctermSeedSet)):
        unimodPeptDict['C-term'] += [((peptide[0][0][::-1], peptide[0][1][::-1]), peptide[1])]

    print 'N-term'
    unimodPeptDict['N-term'] = []
    for peptide in blindModPeptDFS(maxPolyPepMass, list(anywhereSeedSet), maxPolyPepLength = 3, seedList=list(ntermSeedSet)):
        unimodPeptDict['N-term'] += [peptide]

    return unimodPeptDict

def addPepsToAADict(mass, aminoacids=aminoacids):
    aas = copy.copy(aminoacids)
    addList = getPeptidesLessThan(mass)
    for pept in addList:
        try:
            aas[pept]
        except KeyError:
            peptData = ['', '', 0, 0]
            for aa in AAGen(pept):
                for i in range(len(peptData)):
                    peptData[i] += aas[aa][i]
                
                aas[pept] = tuple(peptData)
    
    return aas

def getPeptsOfMaxLength(length=3):
    aas = copy.copy(aminoacids)
    currLengthPepts = copy.copy(aas)
    peptLengths = 1
    while peptLengths < length:
        newLengthPepts = {}
        for pept in currLengthPepts:
            for aa in aminoacids:
                newLengthPepts[pept + aa] = [currLengthPepts[pept][i]+aminoacids[aa][i] for i in range(4)]
        
        aas.update(newLengthPepts)
        currLengthPepts = newLengthPepts
        peptLengths += 1
        
    return aas        

def getMass(peptide):
    mass = 0
    for aa in list(peptide):
        mass += aminoacids[aa][2]
    
    return mass

def AAGen(seq, aaDict=aminoacids):
    nodeGen = nodeInfoGen(seq, aaDict=aaDict, considerTerminalMods=True)
    for node in nodeGen:
        yield node['formAA']

    yield node['lattAA']

def nodeInfoGen(seq, startMass=0, aaDict=aminoacids, addTerminalNodes=False, considerTerminalMods=False, ambigAA='X', ambigEdges=None):

    i = 0
    AAs = []
    aa = ''
    if ambigEdges:
        ambigEdges = copy.copy(list(ambigEdges))
    prm = startMass + (NTermMods['static'] if considerTerminalMods else 0)
    
    while i < len(seq):
        try:
            while i < len(seq):
                aamass = aaDict[aa + seq[i]][2]
                aa += seq[i]
                i += 1
        except KeyError:
            if aa == '' and seq[i] == ambigAA:
                aa = ambigAA
                i += 1
                aamass = ambigEdges[0][1] - ambigEdges[0][0]
                del ambigEdges[0]
            if len(AAs) > 0:
                if aa != '':
                    yield {'prm': prm, 'lattAA': aa, 'formAA': AAs[-1]}
                elif considerTerminalMods and seq[i] in NTermMods and len(AAs) == 1:
                    prm += NTermMods[seq[i]] - NTermMods['static']
                    i += 1
                elif considerTerminalMods and seq[i] in CTermMods and seq[i] == seq[-1]:
                    break
                else:
                    raise KeyError('%s is not a valid amino acid in seq %s' % (aa + seq[i], seq))
            elif addTerminalNodes and aa != '':
                yield {'prm': startMass, 'lattAA': aa, 'formAA': None}
            elif aa == '':
                raise KeyError('%s is not a valid amino acid in seq %s' % (aa + seq[i], seq))
            
            if aa:
                AAs += [aa]
                prm += aamass
            aa = ''
            aamass = 0
    
    
    if i < len(seq):
        aamass += CTermMods[seq[i]]
        aa = AAs[-1]
    elif aa == '':
        aa = AAs[-1] 
    elif len(AAs) > 0:   
        yield {'prm': prm, 'lattAA': aa, 'formAA': AAs[-1]}
        if considerTerminalMods:
            aamass += CTermMods['static']
    
    if addTerminalNodes and aa:
        if len(AAs) == 0:
            yield {'prm': startMass, 'lattAA': aa, 'formAA': None}
            if considerTerminalMods:
                aamass += CTermMods['static']
            
        yield {'prm': prm + aamass, 'lattAA': None, 'formAA': aa}

def getTermModHashForPairConfig(pairConfig):
    N = {}
    C = {}
    for mod in mods.keys():
        for NMod in pairConfig['NStatic']:
            if np.abs(mods[mod] - NMod) < 0.001:
                N[mod] = mods[mod]
        for CMod in pairConfig['CStatic']:
            if np.abs(mods[mod] - CMod) < 0.001:
                C[mod] = mods[mod]
            
    return createTermModHashAAs(N=N, C=C)

def massLadder(seq, startMass=0):
    nodes = nodeInfoGen(seq, startMass)
    for node in nodes:
        print 'Prm: %f, AA: %s' % (node['prm'], node['formAA'])
        
    print 'Prm: %f, AA: %s' % (node['prm'] + aminoacids[node['lattAA']][2], node['lattAA'])
    
if __name__ == '__main__':
    aas = getPeptsOfMaxLength(3)
    for aa in aas:
        print aa, aas[aa]
    print 'length', len(aas)
    aas = addPepsToAADict(300)
    print len(aas)

#    nodeGen = nodeInfoGen('AAAXGHYX', addTerminalNodes=True, considerTerminalMods=True, ambigEdges=[(10,20), (20,40)])
#    for node in nodeGen:
#        print node
