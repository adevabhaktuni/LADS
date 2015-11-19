'''
Created on Jan 24, 2012

@author: arun
'''

from Tkinter import *
import tkFileDialog

import glob
import numpy as np
import time
import pickle
import copy

import ArgLib
import DataFile
import Analytics as An
import ProbNetwork as PN
import DeNovoSequencer as DNS
import Constants

from tkColorChooser import askcolor


class appmanager_tk():
    
    def getHexString(self, rgbProps):
        rgbVals = np.floor(255*rgbProps)
        rVal = hex(int(rgbVals[0]))[2:]
        gVal = hex(int(rgbVals[1]))[2:]
        bVal = hex(int(rgbVals[2]))[2:]
        return '#'+rVal.zfill(2)+gVal.zfill(2)+bVal.zfill(2)             
    
    def curryFileBrowser(self, stringVar, **kwargs):
        def browseFiles():
            stringVar.set(tkFileDialog.askopenfilename(parent=self._root, **kwargs))
        return browseFiles
    
    def curryDirBrowser(self, stringVar, **kwargs):
        def browseDirs():
            stringVar.set(tkFileDialog.askdirectory(parent=self._root, **kwargs))
        return browseDirs
        
    def loadInit(self):
        self._paramsDict = DataFile.parseParams(self._selectedInitFile.get())
        with open('../Misc/symbolmap.txt', 'r') as fin:
            symbolMap = pickle.load(fin)
        self._seqMap = DataFile.generateSeqMap({'LADS Unit Test': 'LADS'}, symbolMap, self._paramsDict)['LADS Unit Test']
        self._aas = Constants.addPepsToAADict(self._minedge)
        
        
    def listDTAs(self):
        self._dtaList = glob.glob(self._selectedDir.get() + '/*.dta')
        if not self._dtaList:
            self._selectedDir.set('No DTAs in selected directory!')
        else:
            self._scanFDict = An.getScanFDict(self._dtaList)
            self._indexedScanFList = np.zeros(len(self._scanFDict))
            for i, scanF in enumerate(sorted(self._scanFDict.keys())):
                self._scanFListbox.insert(END, str(scanF))
                self._indexedScanFList[i] = scanF
    
    def __init__(self, master):
        self._ppm = ArgLib.A['ppm']['attrs']['default']
        self._paircutoff = ArgLib.A['paircutoff']['attrs']['default']
        self._pnet = PN.ProbNetwork('../Scoring_Functions/LysC_likelihood_prior_config.txt', '../Scoring_Functions/LysC_likelihood_prior_model.txt')
        self._ambigpenalty = ArgLib.A['ambigpenalty']['attrs']['default']
        self._minedge = ArgLib.A['minedge']['attrs']['default']
        self._maxedge = ArgLib.A['maxedge']['attrs']['default']
        self._subgraphcut = ArgLib.A['subgraphcut']['attrs']['default']
        self._alpha = ArgLib.A['alpha']['attrs']['default']
        
        self._root = master
        self._numseq = 10
        self.initialize()
        self._root.mainloop()
        
    def getBrowseLoadFrame(self, storeStringVar, browseText, browseCallback, loadText, loadCallback, master=None):
        if master == None:
            selectFrame = Frame(self._root)
        else:
            selectFrame = Frame(master)
        nameEntry = Entry(selectFrame, textvariable=storeStringVar)
        browseButton = Button(selectFrame, text=browseText, command=browseCallback)
        loadButton = Button(selectFrame, text=loadText, command=loadCallback)
        nameEntry.pack(side=LEFT)
        browseButton.pack(side=LEFT)
        loadButton.pack(side=LEFT)
        return selectFrame
    
    def getListBoxWScrollbar(self, listFrame):
        scroller = Scrollbar(listFrame, orient=VERTICAL)
        lbox = Listbox(listFrame, yscrollcommand=scroller.set)
        scroller.config(command=lbox.yview)
        scroller.pack(side=RIGHT, fill=Y)
        lbox.pack(side=LEFT, fill=BOTH, expand=1)
        return lbox
        
    def getPairs(self):
        pairs = {}
        for scanF in self._scanFDict:
            self._scanFDict[scanF]['paired scans'] = {}
        for pairConfigName in self._paramsDict['Pair Configurations']:
            pairConfig = self._paramsDict['Pair Configurations'][pairConfigName]
            pairs[pairConfigName] = An.findDeltaPairs(self._dtaList, pairConfig['NMod']+pairConfig['CMod'], ppm=self._ppm)
            for pair in pairs[pairConfigName]:
                pairData = {'light': pair[0], 'heavy': pair[1], 'pair configuration': pairConfigName, 'pair score': None, 'light precmass': self._scanFDict[pair[0]]['precMass'], 'heavy precmass': self._scanFDict[pair[1]]['precMass']}
                self._scanFDict[pair[0]]['paired scans'][pair[1]] = pairData
                self._scanFDict[pair[1]]['paired scans'][pair[0]] = pairData
    
    def listPairs(self, event):
        curScanF = self._indexedScanFList[int(self._scanFListbox.curselection()[0])]
        self._indexedPairData = {}
        self._pairedScanListbox.delete(0, END)
        self._pairedScanListbox.insert(END, '<Unpaired>')
        self._indexedPairData[0] = {'light': curScanF, 'light precmass': self._scanFDict[curScanF]['precMass'], 'heavy precmass': 'N/A', 'heavy': 'N/A', 'pair configuration': 'N/A', 'pair score': 'N/A',}
        if 'paired scans' in self._scanFDict[curScanF]:
            for i, pairedScan in enumerate(sorted(self._scanFDict[curScanF]['paired scans'].keys())):
                self._indexedPairData[i+1] = self._scanFDict[curScanF]['paired scans'][pairedScan]
                self._pairedScanListbox.insert(END, str(pairedScan))
        
    def updatePairInfo(self, event):
        curPairedScanData = self._indexedPairData[int(self._pairedScanListbox.curselection()[0])]
        if curPairedScanData['pair score'] == None:
            curPairedScanData['pair score'] = An.getSharedPeaksRatio(self._scanFDict[curPairedScanData['light']]['dta'], self._scanFDict[curPairedScanData['heavy']]['dta'], self._paramsDict['Pair Configurations'][curPairedScanData['pair configuration']], epsilon=self._ppm * 10**-6 * curPairedScanData['light precmass'])
        for labelDatum in self._pairInfoLabelVars:
            self._pairInfoLabelVars[labelDatum].set(str(curPairedScanData[labelDatum]))
        
        if curPairedScanData['pair score'] != 'N/A':
            if curPairedScanData['pair score'] > self._paircutoff:
                self._pairScoreLabel.config(fg='dark green')
            else:
                self._pairScoreLabel.config(fg='red')
        else:
            self._pairScoreLabel.config(fg='black')
        
    def initializePairInfoDisplay(self):
        pairedScansFrame = LabelFrame(self._root, text='Paired Scans', padx=5, pady=5)
        
        pairedScanListFrame = Frame(pairedScansFrame)
        self._pairedScanListbox = self.getListBoxWScrollbar(pairedScanListFrame)
        pairedScanListFrame.grid(rowspan=2)
        self._pairedScanListbox.bind('<Double-Button-1>', self.updatePairInfo)
        
        pairInfoFrame = LabelFrame(pairedScansFrame, text='Pair Information', fg='blue')
        self._pairInfoLabelVars = {'light': 0, 'heavy': 0, 'pair configuration': 0, 'pair score': 0, 'light precmass': 0, 'heavy precmass': 0}
        for labelDatum in self._pairInfoLabelVars:
            self._pairInfoLabelVars[labelDatum] = StringVar()
            self._pairInfoLabelVars[labelDatum].set('N/A')
        
        for i, labelText in enumerate(['Light ScanF: ', 'Light Precursor Mass: ', 'Heavy ScanF: ', 'Heavy Precursor Mass: ', 'Pair Configuration: ', 'Shared Peaks Ratio: ']):
            Label(pairInfoFrame, text=labelText).grid(row=i, sticky=W)
        
        for i, labelDatum in enumerate(['light', 'light precmass', 'heavy', 'heavy precmass', 'pair configuration']):
            Label(pairInfoFrame, textvariable=self._pairInfoLabelVars[labelDatum]).grid(row=i, column=1, sticky=E)
        
        self._pairScoreLabel = Label(pairInfoFrame, textvariable=self._pairInfoLabelVars['pair score'])
        self._pairScoreLabel.grid(row=len(self._pairInfoLabelVars)-1, column=1, sticky=E)
        
        seqButton = Button(pairedScansFrame, text="Sequence!", command=self.sequenceDTAs)
        
        pairInfoFrame.grid(column=1, row=0)
        seqButton.grid(column=1, row=1)
        pairedScansFrame.pack()
    
    def sequenceDTAs(self):
        curPairedScanData = self._indexedPairData[int(self._pairedScanListbox.curselection()[0])]
        t1 = time.time()
        if curPairedScanData['heavy'] != 'N/A':
            heavySeqMap = copy.deepcopy(self._seqMap)
            heavySeqMap['Mods']['N-Term'] = self._paramsDict['Pair Configurations'][curPairedScanData['pair configuration']]['NModSymbol']
            heavySeqMap['Mods']['C-Term'] = self._paramsDict['Pair Configurations'][curPairedScanData['pair configuration']]['CModSymbol']
            sharedInfo, starts, ends, deltas, termModHash, specs, G = DNS.initializeSpectrumGraph(self._pnet, self._paramsDict, self._scanFDict[curPairedScanData['light']]['dta'], heavyPath=self._scanFDict[curPairedScanData['heavy']]['dta'], ppm=self._ppm, usePaired=True, pairConfigName=curPairedScanData['pair configuration'], verbose=False)
            precMass = sharedInfo['lightPrecMass']
        else:
            sharedInfo, starts, ends, deltas, termModHash, specs, G = DNS.initializeSpectrumGraph(self._pnet, self._paramsDict, self._scanFDict[curPairedScanData['light']]['dta'], ppm=self._ppm, verbose=False)
            precMass = sharedInfo['precMass']
        
        epsilon = self._ppm * precMass * 10 ** -6
        paths, subG = DNS.getSpectrumGraphPaths(G, deltas, specs, starts, ends, precMass - Constants.mods['H+'] - Constants.mods['H2O'], termModHash=termModHash, unknownPenalty=self._ambigpenalty, maxEdge=self._maxedge, minEdge=self._minedge, subGraphCut=self._subgraphcut, subAlpha=0.3, alpha=self._alpha, epsilon=epsilon, aas=self._aas, verbose=False)
        seqTime = time.time() - t1
        if paths:
            seqs = []
            for path in paths:
                seqs.extend([DNS.getSequenceFromNodes(subG, path[1], precMass - Constants.mods['H+'] - Constants.mods['H2O'], termModHash)])
    
            scores = list(zip(*paths)[0])
            Ord = np.argsort(-1 * np.array(scores))
            
            ambigEdges = []
            numAmbig = 0
            for j in range(self._numseq):
                try:
                    for i in range(len(seqs[Ord[j]])):
                        if type(seqs[Ord[j]][i]) == tuple:
                            ambigEdges.extend([seqs[Ord[j]][i]])
                            numAmbig += 1
                            seqs[Ord[j]][i] = '-'
                
                    curSeq = ''.join(seqs[Ord[j]])
                    curSeq = An.preprocessSequence(curSeq, self._seqMap, ambigEdges=ambigEdges)
                    if j == 0 and curPairedScanData['heavy'] != 'N/A':
                        try:
                            curHeavySeq = An.preprocessSequence(curSeq, heavySeqMap, replaceExistingTerminalMods=True, ambigEdges=ambigEdges)
                            AAs = An.getAllAAs(curHeavySeq, ambigEdges=ambigEdges)
                            self._seqStatus.set('Paired Sequencing Successful! Heavy Sequence: %s. Time taken: %f seconds' % (curHeavySeq, seqTime))
                        except KeyError:
                            self._seqStatus.set('ERROR: Heavy Sequence %s is not a valid sequence! Time wasted: %f seconds' % (curHeavySeq, seqTime))
                    elif j == 0:
                        self._seqStatus.set('Unpaired Sequencing Successful! Time taken: %f seconds' % (seqTime))
                    
                    for labelInst in self._seqScoreData[j]['seq'].children.values():
                        labelInst.destroy()
                    self.displayConfColoredSequence(subG, self._seqScoreData[j]['seq'], paths[Ord[j]][1], curSeq, ambigEdges=ambigEdges)
                    self._seqScoreData[j]['score'].set(str(scores[Ord[j]]))
                except IndexError:
                    for labelInst in self._seqScoreData[j]['seq'].children.values():
                        labelInst.destroy()
                    self._seqScoreData[j]['score'].set('')
        else:
            self._seqStatus.set('ERROR: No Sequences Found! Time wasted: %f seconds' % seqTime)
                        
    def displayConfColoredSequence(self, G, masterFrame, path, seq, ambigEdges=None):
        nodeGen = Constants.nodeInfoGen(seq, addTerminalNodes=False, considerTerminalMods=True, ambigEdges=ambigEdges)
        prevNode = None
        for i, node in enumerate(nodeGen):
            print node, path[i+1]
            node['prm'] = path[i+1]
            confScore = An.getAAConfidence(G, prevNode=prevNode, nextNode=node)
            prevNode = node
            hexColor = self.getHexString(np.array([1-confScore, confScore, 0]))
            if prevNode == None and seq[len(node['formAA'])] in Constants.NTermMods:
                Label(masterFrame, text=node['formAA']+seq[len(node['formAA'])], fg='white', bg=hexColor).pack(side=LEFT)
            else:
                Label(masterFrame, text=node['formAA'], fg='white', bg=hexColor).pack(side=LEFT)
        
        confScore = An.getAAConfidence(G, prevNode=prevNode, nextNode=None)
        hexColor = self.getHexString(np.array([1-confScore, confScore, 0]))
        if seq[-1] in Constants.CTermMods:
            Label(masterFrame, text=node['lattAA']+seq[-1], fg='white', bg=hexColor).pack(side=LEFT)
        else:
            Label(masterFrame, text=node['lattAA'], fg='white', bg=hexColor).pack(side=LEFT)
        
    def initializeSeqInfoDisplay(self):
        seqInfoFrame = LabelFrame(self._root, text='Sequence Information', padx=5, pady=5)
        
        seqResultsFrame = LabelFrame(seqInfoFrame, text='Sequence Results')
        Label(seqResultsFrame, text='Score').grid(row=0, sticky=W)
        Label(seqResultsFrame, text='Sequence').grid(row=0, column=1, sticky=E)
        self._seqScoreData = {}
        for i in range(self._numseq):
            self._seqScoreData[i] = {'seq': Frame(seqResultsFrame), 'score': StringVar()}
            Label(seqResultsFrame, textvariable=self._seqScoreData[i]['score'], fg='blue').grid(row=i+1, column=0, sticky=W)
            self._seqScoreData[i]['seq'].grid(row=i+1, column=1, sticky=E)
        
        self._seqStatus = StringVar()
        seqStatusFrame = LabelFrame(seqInfoFrame, text='Status', fg='Red')
        Label(seqStatusFrame, textvariable=self._seqStatus, anchor=W, justify=LEFT).pack()
        seqStatusFrame.pack(side=BOTTOM)
        seqResultsFrame.pack()
        seqInfoFrame.pack(side=BOTTOM)
        
        
            
    
    def initialize(self):
        self._selectedInitFile = StringVar()
        initSelectFrame = self.getBrowseLoadFrame(self._selectedInitFile, 'Browse', self.curryFileBrowser(self._selectedInitFile, defaultextension='ini'), 'Load Init', self.loadInit)
        initSelectFrame.pack(side=TOP)
        
        self._selectedDir = StringVar()
        dtaDirSelectFrame = self.getBrowseLoadFrame(self._selectedDir, 'Browse', self.curryDirBrowser(self._selectedDir, mustexist=True), 'List DTAs', self.listDTAs)
        dtaDirSelectFrame.pack(side=TOP)
        
        labelDTAListFrame = LabelFrame(self._root, text='DTA List', padx=5, pady=5)
        scanFListFrame = Frame(labelDTAListFrame, padx=5, pady=5)
        self._scanFListbox = self.getListBoxWScrollbar(scanFListFrame)
        scanFListFrame.pack(side=TOP)
        b = Button(labelDTAListFrame, text="Get Pairs!", command=self.getPairs)
        b.pack(side=BOTTOM)
        self._scanFListbox.bind('<Double-Button-1>', self.listPairs)
        labelDTAListFrame.pack()
        
        self.initializePairInfoDisplay()
        self.initializeSeqInfoDisplay()
        
    
    
if __name__ == '__main__':
    root = Tk()
    appmanager_tk(root)