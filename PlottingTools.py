'''
Created on Sep 19, 2011

@author: Arun
'''

try:
    import pylab as plt
except ImportError:
    pass

import DataFile
import anydbm as dbm
import pickle
import numpy as np

def getScanComparisonInfo(scansInfo, infoMap, progs=['LADS', 'MASCOT', 'SEQUEST'], scanFields=['Score'], compFields=['Precision', 'Accuracy']):
    compInfo = {}
    
    keys = []
    for field in scanFields:
        keys.extend([prog + ' ' + infoMap[prog][field] for prog in progs])
    
    for cField in compFields:
        for i in range(len(progs)):
            for j in range(i + 1, len(progs)):
                try:
                    scansInfo[0][progs[i] + ' ' + progs[j] + ' ' + cField]
                    keys.extend([progs[i] + ' ' + progs[j] + ' ' + cField])
                except KeyError:
                    keys.extend([progs[j] + ' ' + progs[i] + ' ' + cField])
    
    for scan in scansInfo:
        numProgs = 0
        for prog in progs:
            if scan[prog + ' ' + infoMap[prog]['Peptide']] != 'None':
                numProgs += 1
        
        if numProgs > 1:
            scanInfo = {}
            for key in keys:
                scanInfo[key] = scan[key]
            compInfo[int(scan['ScanF'])] = scanInfo
    
    return compInfo

def scatterPlot(compInfo, axis1, axis2):
    axis1Vals = []
    axis2Vals = []
    for scanF in compInfo.keys():
        if compInfo[scanF][axis1] != 'None' and compInfo[scanF][axis2] != 'None':
            axis1Vals.extend([float(compInfo[scanF][axis1])])
            axis2Vals.extend([float(compInfo[scanF][axis2])])
    
    plt.scatter(axis1Vals, axis2Vals)
    plt.xlabel(axis1)
    plt.ylabel(axis2)
    plt.show()
    
if __name__ == '__main__':
    scansfName = 'compareSearches_MASCOT_LADSUPen10KPen15All_SEQUEST_ath001862.tdv'
    scansInfo = DataFile.getScanInfo(scansfName, delimiter='\t')
    
    infoMap = {'MASCOT': {'Score': 'Ion Score', 'Peptide': 'Peptide', 'Reference': 'Reference'},
               'SEQUEST': {'Score': 'XCorr', 'Peptide': 'Peptide', 'Reference': 'Reference'},
               'LADS': {'Score': 'PScore', 'Peptide': 'Peptide', 'Reference': 'Reference'}}
    
    compInfo = getScanComparisonInfo(scansInfo, infoMap)
    scatterPlot(compInfo, 'SEQUEST XCorr', 'LADS PScore')
    
    
    
