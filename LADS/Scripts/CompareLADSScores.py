'''
Created on Sep 23, 2011

@author: Arun
'''

import Analytics as An
import Constants
import DataFile
import ProbNetwork as PN

import ArgLib
import glob
import anydbm as dbm
import pickle
import os

pScoreErr = 0.5

if __name__ == '__main__':
    options = ArgLib.parse(['dtadir', 'output', 'config', 'model', 'comp', 'database', 'ppm', 'ambigpenalty', 'progdict'])
    
    Constants.aminoacids['I'] = Constants.aminoacids['J']
    Constants.aminoacids['J'] = (Constants.aminoacids['K'][0], Constants.aminoacids['K'][1], Constants.aminoacids['K'][2] + Constants.mods['*'], Constants.aminoacids['K'][3])
    Constants.aminoacids['O'] = (Constants.aminoacids['M'][0], Constants.aminoacids['M'][1] + 'O', Constants.aminoacids['M'][2] + Constants.mods['#'], Constants.aminoacids['M'][3])
    Constants.aminoacids['C'] = (Constants.aminoacids['C'][0], Constants.aminoacids['C'][1], Constants.aminoacids['C'][2] + Constants.mods['Carbamidomethyl'], Constants.aminoacids['C'][3])
    
    revMap = {}
    for aa in Constants.aminoacids.keys():
        revMap[aa] = aa
    
    revMap['I'] = 'J'
    revMap['J'] = 'K'
    revMap['X'] = 'X'
    revMap['O'] = 'M'
    
    PNet = PN.ProbNetwork(options.config, options.model)
    
    db = dbm.open(options.database, 'r')
    infoMap = pickle.loads(db['infoMap'])
    db.close()
    
    scansInfo = DataFile.getScanInfo(options.comp, delimiter='\t')
    
    progs = []
    for key in scansInfo[0].keys():
        prog = key.split(' ', 1)[0]
        if prog in An.searchprogs and prog not in progs:
            progs.extend([prog])
            
    compInfo = An.getScanComparisonInfo(scansInfo, infoMap, progs=progs, scanFields=['Score', 'Peptide'], compFields=[], specificCols=['LADS Ambiguous Edges'])
    del scansInfo
    
    cols = compInfo[compInfo.keys()[0]].keys()
    
    stats = {}
    for prog in progs:
        cols.extend([prog + ' Result PScore'])
        stats[prog + ' better'] = 0
        stats[prog + ' worse'] = 0
        stats[prog + ' equal'] = 0
        stats[prog + ' none'] = 0
        
    cols.sort()
    cols = ['ScanF'] + cols
    
    if options.output:
        outFile = open(options.output, 'w')
        outFile.write('\t'.join([col for col in cols]) + '\n')
    
    
    for scanF in sorted(compInfo.keys()):
        scanData = {}
        dtaPath = glob.glob(options.dtadir + '*' + '%04d' % scanF + '*.dta')[0]
        scanData['ScanF'] = scanF
        for col in cols:
            try:
                scanData[col] = compInfo[scanF][col]
            except KeyError:
                pass
        
        print '\n Scan Number %i \n' % scanF
        print '\t'.join(['Program', 'Peptide', 'Original Score', 'LADS Score'])
        for prog in progs:
            if prog == 'LADS':
                ambigEdges = eval(scanData['LADS Ambiguous Edges'])
            else:
                ambigEdges = None

            if scanData[prog + ' ' + infoMap[prog]['Peptide']] != 'None':    
                pScore = An.getLADSPScore(seq=scanData[prog + ' ' + infoMap[prog]['Peptide']], dtaPath=dtaPath, PNet=PNet, ppm=options.ppm, ambigEdges=ambigEdges, revMap=revMap, ambigAAPen=options.ambigpenalty)
            else:
                pScore = None

            scanData[prog + ' Peptide PScore'] = pScore
            if not pScore:
                stats[prog + ' none'] += 1
            elif pScore > float(scanData['LADS ' + infoMap['LADS']['Score']]) + pScoreErr:
                stats[prog + ' better'] += 1
            elif pScore < float(scanData['LADS ' + infoMap['LADS']['Score']]) - pScoreErr:
                stats[prog + ' worse'] += 1
            else:
                stats[prog + ' equal'] += 1
            
            peptide = str(scanData[prog + ' ' + infoMap[prog]['Peptide']])
            oldScore = infoMap[prog]['Score'] + ': ' + str(scanData[prog + ' ' + infoMap[prog]['Score']])
            print '\t'.join([prog, peptide, oldScore, str(pScore)])
        
        if options.output:    
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
    
    print '\n Stats \n'
    print '\t'.join(['Program', 'Better', 'Equal', 'Worse', 'No Result'])
    for prog in progs:
        print '\t'.join([prog, str(stats[prog + ' better']), str(stats[prog + ' equal']), str(stats[prog + ' worse']), str(stats[prog + ' none'])])
        
        
    
    
