import sys
import os
sys.path.insert(1, os.path.abspath('../LADS-rev120/'))

import numpy as np
import os
import shlex, subprocess
import time

import ArgLib
import DataFile
import Analytics as An


def parseBLASTOutput(blastOutput):
    cols = ['qid', 'sid', 'E', 'N', 'Sprime', 'S', 'alignlen', 'nident', 'npos', 'nmism', 'pcident', 'pcpos', 'qgaps', 'qgaplen', 'sgaps', 'sgaplen', 'qframe', 'qstart', 'qend', 'sframe', 'sstart', 'send']

    blastOut = open(blastOutput)

    blastDict = {}

    for line in blastOut:
        scanDict = dict( zip(cols, line.strip().split('\t')) )

        qid = int(scanDict['qid'])
        if qid not in blastDict:
            blastDict[qid] = { 'references': [( scanDict['sid'], int(scanDict['sstart']), int(scanDict['send']) )], 'min e-value': float(scanDict['E']), 'align score': int(scanDict['S']), 'align length': int(scanDict['alignlen']), 'q gap length': int(scanDict['qgaplen']), 'qstart': int(scanDict['qstart']), 'qend': int(scanDict['qend']), 'nident': int(scanDict['nident']), 'npos': int(scanDict['npos']), 's gap length': int(scanDict['sgaplen']) }
        else:
            alignScore = int(scanDict['S'])

            if alignScore == blastDict[qid]['align score']:
                blastDict[qid]['references'] += [( scanDict['sid'], int(scanDict['sstart']), int(scanDict['send']) )]
                blastDict[qid]['min e-value'] = min(( float(scanDict['E']), blastDict[qid]['min e-value']))
                
            elif alignScore > blastDict[qid]['align score']:
                blastDict[qid] = { 'references': [( scanDict['sid'], int(scanDict['sstart']), int(scanDict['send']) )], 'min e-value': float(scanDict['E']), 'align score': int(scanDict['S']), 'align length': int(scanDict['alignlen']), 'q gap length': int(scanDict['qgaplen']), 'qstart': int(scanDict['qstart']), 'qend': int(scanDict['qend']), 'nident': int(scanDict['nident']), 'npos': int(scanDict['npos']), 's gap length': int(scanDict['sgaplen']) }                        

    return blastDict
            
def writeBLASTResults(parsedBLASTOutput, peptideDict, outFileName):
    outFile = open(outFileName, 'w')

    cols = ['QID', 'Peptide', 'Unmod Peptide', 'Ambig Edges', 'References', 'Alignment Score', 'Min E-Value', 'Query Start', 'Query End', 'Num Identical', 'Num Positive', 'Query Gap Length', 'Subject Gap Length', 'Alignment Length']

    outFile.write('\t'.join(cols) + '\n')
    for qid in parsedBLASTOutput:
        print qid
        item = parsedBLASTOutput[qid]
        print item
        scanData = {'QID': qid, 'References': item['references'], 'Peptide': peptideDict[qid][0], 'Unmod Peptide': peptideDict[qid][1], 'Ambig Edges': peptideDict[qid][2], 'Alignment Score': item['align score'], 'Min E-Value': item['min e-value'], 'Query Start': item['qstart'], 'Query End': item['qend'], 'Num Identical': item['nident'], 'Num Positive': item['npos'], 'Alignment Length': item['align length'], 'Query Gap Length': item['q gap length'], 'Subject Gap Length': item['s gap length'] }
        
        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

    outFile.close()

def compareWithDBSearchResults(blastResults, indexedDBResults, dbPeptideKey, dbReferenceKey, outFileName):
    cols = ['ScanF', 'De Novo Peptide', 'DB Peptide', 'Longest Shared Substring Length', 'DB Reference', 'BLAST References', 'BLAST Num Identical', 'BLAST DB Agree?']
    outFile = open(outFileName, 'w')
    outFile.write('\t'.join(cols) + '\n')

    numAgree, numDisagree = 0, 0
    for scanF in blastResults:
        if scanF in indexedDBResults:
            dbItem = indexedDBResults[scanF]
            if dbItem[dbPeptideKey] != 'None':
                scanData = {'ScanF': scanF, 'De Novo Peptide': blastResults[scanF]['Peptide'], 'DB Peptide': dbItem[dbPeptideKey], 'DB Reference': dbItem[dbReferenceKey], 'BLAST References': blastResults[scanF]['References'], 'BLAST Num Identical': blastResults[scanF]['Num Identical']}

                if dbItem[dbReferenceKey].replace('|', '_') in [ref[0] for ref in eval(blastResults[scanF]['References'])]:
                    scanData['BLAST DB Agree?'] = True
                    numAgree += 1
                else:
                    scanData['BLAST DB Agree?'] = False
                    numDisagree += 1

                scanData['Longest Shared Substring Length'] = len(long_substr([blastResults[scanF]['Unmod Peptide'], An.stripModifications(dbItem[dbPeptideKey])]))

                outFile.write('\t'.join([str(scanData[col]) for col in cols])  + '\n')

    outFile.close()



# Taken from stack overflow url http://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings-python
def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]

    return substr
                
def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False

    return True
                
def compileProteinStats(BLASTresults, procCompInfo, infoMap, mainProgName):
    prots = {}
    for scan in BLASTresults:
        scanF = int(scan['ScanF'])
        protList = eval(scan['Desc'])
        for prot in protList:
            try:
                prots[prot]['Max P-Score'] = max(prots[prot]['Max P-Score'], float(procCompInfo[scanF][mainProgName + ' ' + infoMap['LADS']['Score']]))
                prots[prot]['Min E-val'] = min(prots[prot]['Min E-val'], float(scan['E-Value']))
                prots[prot]['Recip E-val Sum'] += 1/float(scan['E-Value'])
                prots[prot]['Peptides'] += 1
            except KeyError:
                prots[prot] = {}
                prots[prot]['Max P-Score'] = float(procCompInfo[scanF][mainProgName + ' ' + infoMap['LADS']['Score']])
                prots[prot]['Min E-val'] = float(scan['E-Value'])
                prots[prot]['Recip E-val Sum'] = 1/float(scan['E-Value'])
                prots[prot]['Peptides'] = 1

    return prots

def getDBProteinStats(procCompInfo, infoMap, progDict, mainProgName):
    dbProts = {}
    for progName in progDict:
        if progName != mainProgName:
            dbProts[progName] = {}
            
    for scanF in procCompInfo:
        for progName in dbProts:
            if procCompInfo[scanF][progName + ' ' + infoMap[progDict[progName]]['Score']] != 'None':
                prot = procCompInfo[scanF][progName + ' ' + infoMap[progDict[progName]]['Reference']]
                if prot in dbProts[progName]:
                    dbProts[progName][prot]['Max ' + infoMap[progDict[progName]]['Score']] = max(dbProts[progName][prot]['Max ' + infoMap[progDict[progName]]['Score']], float(procCompInfo[scanF][progName + ' ' + infoMap[progDict[progName]]['Score']]))
                    dbProts[progName][prot]['Peptides'] += 1
                else:
                    dbProts[progName][prot] = {}
                    dbProts[progName][prot]['Max ' + infoMap[progDict[progName]]['Score']] = float(procCompInfo[scanF][progName + ' ' + infoMap[progDict[progName]]['Score']])
                    dbProts[progName][prot]['Peptides'] = 1

    return dbProts

def getNonRedundantReferences(prots, dbProts):
    nrProts = {}
    for prot in prots:
        nrProts[prot] = 0

    for progName in dbProts:
        for prot in dbProts[progName]:
            nrProts[prot] = 0

    return nrProts.keys()

def writeFASTAFile(runInfo, peptideKey, FASTAout, ambigEdgesKey=None, ambigEdgeCutoff=1):
    outFile = open(FASTAout, 'w')
    fastaDict = {}

    for scanF in runInfo:
        scan = runInfo[scanF]
        
        if ambigEdgesKey != None:
            ambigEdgesLength = eval(scan[ambigEdgesKey])
            if type(ambigEdgesLength) != int:
                ambigEdgesLength = len(ambigEdgesLength)

        unmodPept = An.stripModifications(scan[peptideKey])
        if scan[peptideKey] != 'None' and (ambigEdgesKey == None or ambigEdgesLength <= ambigEdgeCutoff):
            outFile.write('>' + str(scanF) + '\n')
            outFile.write(unmodPept + '\n')

            fastaDict[scanF] = [scan[peptideKey], unmodPept, [] if ambigEdgesKey == None else scan[ambigEdgesKey]]
        
    outFile.close()
    
    return fastaDict


if __name__ == '__main__':
    print 'This program will take peptides from a tabular formatted file and BLAST them against the desired proteome using AB-BLAST. LADS option is input dataset, init option is name of the peptide column in dataset. sequest option is name of ambig edge key in dataset (leave blank if none)'
    options = ArgLib.parse(['output', 'lads', 'init', 'sequest'], [
        {'opts': ('-b', '--blastout'), 'attrs': {'type': 'string', 'dest': 'blastout', 'help': 'AB BLAST output file. If it exists, can compile protein stats without re-computing BLAST alignments.'}},
        {'opts': ('-d', '--blastdb'), 'attrs': {'type': 'string', 'dest': 'blastdb', 'help': 'Database to BLAST against.'}},
        {'opts': ('-m', '--matrix'), 'attrs': {'type': 'string', 'dest': 'matrix', 'help': 'Scoring matrix to use for BLAST'}},
        {'opts': ('-Q', '--gapopen'), 'attrs': {'type': 'int', 'dest': 'gapopen', 'help': 'gap open penalty for BLAST'}},
        {'opts': ('-R', '--gapextend'), 'attrs': {'type': 'int', 'dest': 'gapextend', 'help': 'gap extend penalty for BLAST.'}},
        {'opts': ('-T', '--threshold'), 'attrs': {'type': 'int', 'dest': 'threshold', 'help': 'word score threshold for BLAST'}}])
    
                         
    
#    progDict = eval(options.progdict)
    compInfo = DataFile.getScanInfo(options.lads, delimiter='\t')
    compInfo = DataFile.indexDataByKey(compInfo, 'ScanF', dtyper=lambda num: int(float(num)))
    #compInfo = DataFile.indexDataByKey(compInfo, 'ScanF', overrideKey = 'LADS Post Score', dtyper=float)

#    procCompInfo = An.getScanComparisonInfo(compInfo, dbInfo['infoMap'], progDict, scanFields=['Reference', 'Score'], compFields = [], specificColDict={}, minNumProgs=1)
    curDir = os.getcwd()
    if options.output != None:
        outFile = open(options.output, 'w')
        FASTAout = os.path.splitext(options.output)[0] + '.fasta'
#        writeFASTAFile(compInfo, 'LADS Sequence', FASTAout, ambigEdgesKey='LADS Ambig Edges', ambigEdgeCutoff=5)
        peptideDict = writeFASTAFile(compInfo, options.init, FASTAout, ambigEdgesKey=options.sequest, ambigEdgeCutoff=2)
#        writeFASTAFile(compInfo, 'LADS Unit Test Peptide', FASTAout, ambigEdgesKey='LADS Unit Test Ambiguous Edges', ambigEdgeCutoff=2)
        
        t1 = time.time()

        os.chdir('/lab/seq_db/user/adevabhaktuni/ABblastpDBs/')
        blastOut = os.path.splitext(options.output)[0] + '_ABBLAST.txt'
        #cmd = '../ABBLAST/blastp %s %s -B 100 -V 1 -E 100000 -spoutmax 1 -sort_by_highscore -mformat 2 -matrix %s -T %i -Q %i -R %i -o %s' % (options.blastdb, FASTAout, options.matrix, options.threshold, options.gapopen, options.gapextend, blastOut)
        # For bigger databases
        cmd = '../ABBLAST/blastp %s %s -B 1000 -V 1 -E 100000 -W 5 -spoutmax 1 -sort_by_highscore -mformat 2 -matrix %s -T %i -Q %i -R %i -o %s -nogaps -wstrict' % (options.blastdb, FASTAout, options.matrix, options.threshold, options.gapopen, options.gapextend, blastOut)
        print cmd

        proc = subprocess.Popen(shlex.split(cmd))
        print proc.communicate()
        print 'Time taken: ', time.time()-t1
        
    elif options.blastout != None:
        blastOut = options.blastout

    
    BLASTresults = parseBLASTOutput(blastOut)
    writeBLASTResults(BLASTresults, peptideDict, options.output)
