'''
Created on Sep 2, 2011

@author: Arun
'''
#Take dta.tgz file, find paired spectra, and spit out sequences

import glob
from optparse import OptionParser
import os
import time

import Constants
import DeNovoSequencer as DNS
import Analytics
import ProbNetwork as PN
import DataFile

parser = OptionParser(usage="%prog [options]")
parser.add_option('-d', '--dir', type='string', dest='dtaDir', help='Directory of .dta files or dta.tgz file')
parser.add_option('-c', '--config', type='string', dest='config', help='Path to model configuration file')
parser.add_option('-m', '--model', type='string', dest='model', help='Path to probabilistic model file')

parser.add_option('-v', '--verbose', action='store_true', default=False, dest='verbose', help='Print status updates to std out')
parser.add_option('-P', '--pairCutoff', type='float', dest='pairCutoff', default=0.1, help='Cutoff shared peaks ratio when determining paired spectra')
parser.add_option('-N', '--Nmod', type='float', dest='Nmod', default=0, help='N-terminal mass difference between paired spectra')
parser.add_option('-C', '--Cmod', type='float', dest='Cmod', default=Constants.mods['*'], help='C-terminal mass difference between paired spectra')
parser.add_option('-p', '--ppm', type='float', dest='ppm', default=5, help='PPM tolerance in mass accuracy')
parser.add_option('-a', '--ambigEdgePenalty', type='float', dest='ambigEdgePenalty', default=5, help='Score penalty for ambiguous edges in spectrum graph')
parser.add_option('-e', '--minEdge', type='float', dest='minEdge', default=160, help='Minimum length of ambiguous edges in spectrum graph')
parser.add_option('-E', '--maxEdge', type='float', dest='maxEdge', default=500, help='Maximum length of ambiguous edges in spectrum graph')
parser.add_option('-l', '--alpha', type='float', dest='alpha', default=0.95, help='Minimum score ratio between suboptimal and optimal paths')

if __name__ == '__main__' :
    (options, args) = parser.parse_args()
    
    if not options.dtaDir or not options.model or not options.config:
        print 'ERROR: missing model, config, or dtaDir'
        exit(-1)
    
    Constants.aminoacids['C'] = (Constants.aminoacids['C'][0], Constants.aminoacids['C'][1], Constants.aminoacids['C'][2] + Constants.mods['Carbamidomethyl'], Constants.aminoacids['C'][3])
    Constants.aminoacids['O'] = (Constants.aminoacids['M'][0], Constants.aminoacids['M'][1] + 'O', Constants.aminoacids['M'][2] + Constants.mods['#'], Constants.aminoacids['M'][3])
    if options.Cmod == Constants.mods['*']:
        Constants.aminoacids['X'] = (Constants.aminoacids['K'][0], Constants.aminoacids['K'][1], Constants.aminoacids['K'][2] + Constants.mods['*'], Constants.aminoacids['K'][3])
    
    PNet = PN.ProbNetwork(options.config, options.model)
    if options.verbose:
        t1 = time.time()
        print 'Getting heavy-light pairs'

    dtaList = glob.glob(options.dtaDir + '/*.dta')
    (paired, unpaired) = Analytics.getPairedAndUnpairedSpectra(options.dtaDir, dtaList, delta=(options.Nmod + options.Cmod), ppm=options.ppm, cutOff=options.pairCutoff)
    if options.verbose:
        t2 = time.time()
        print 'Finished getting paired spectra. Time taken: ', t2 - t1
        print 'Starting Sequencing'
    
    aas = Constants.addPepsToAADict(options.minEdge)
    for pair in paired:
        (lightSpec, heavySpec) = pair[1:]
        if options.verbose:
            print 'Now sequencing %s %s with shared peaks ratio %f' % (lightSpec, heavySpec, pair[0])
            s1 = time.time()
            
        heavyPath = heavySpec
        lightPath = lightSpec
        sharedInfo = DNS.getPairedSpectraInfoForSequencing(lightPath, heavyPath, options.verbose)
        DNS.sequencePairedSpectra(sharedInfo['NInd'], sharedInfo['CInd'], sharedInfo['lightPairs'], sharedInfo['heavyPairs'], sharedInfo['lightPrecMass'] - Constants.mods['H+'] - Constants.mods['H2O'], PNet, alpha=options.alpha, unknownPenalty=options.ambigEdgePenalty, maxEdge=options.maxEdge, minEdge=options.minEdge, Nmod=options.Nmod, Cmod=options.Cmod, aas=aas, verbose=options.verbose)
        
        if options.verbose:
            print 'Time taken:', time.time() - s1 
    
    for spec in unpaired:
        if options.verbose:
            print 'Now sequencing unpaired spectrum %s' % spec
            s1 = time.time()
            
        precMass = DataFile.getPrecMassAndCharge(spec)[0]
        pairs = DataFile.getMassIntPairs(spec)
        DNS.sequenceSingleSpectrum(pairs, precMass - Constants.mods['H+'] - Constants.mods['H2O'], PNet, alpha=options.alpha, unknownPenalty=options.ambigEdgePenalty, maxEdge=options.maxEdge, minEdge=options.minEdge, aas=aas, verbose=options.verbose)
        
        if options.verbose:
            print 'Time taken:', time.time() - s1
    
    if options.verbose:
        print 'Finished sequencing. Time taken: ', time.time() - t2
        print 'Total time taken for program: ', time.time() - t1
            
            
