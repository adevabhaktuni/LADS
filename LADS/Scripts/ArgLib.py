'''
Created on Sep 23, 2011

@author: Arun
'''
import sys
from optparse import OptionParser

import Constants
import DataFile

#{'opts': (), 'attrs': {}}
A = {
     
'dtadir': {'opts': ('-d', '--dtadir'), 'attrs': {'type': 'string', 'dest': 'dtadir', 'help': 'Directory of .dta files or dta.tgz file'}},

'config': {'opts': ('-c', '--config'), 'attrs': {'type': 'string', 'dest': 'config', 'help': 'Path to model configuration file. If none is provided, will check the params file.'}},

'model': {'opts': ('-m', '--model'), 'attrs': {'type': 'string', 'dest': 'model', 'help': 'Path to probabilistic model file. If none is provided, will check the params file.'}},

'output': {'opts': ('-o', '--output'), 'attrs': {'type': 'string', 'dest': 'output', 'help': 'Name of output file (unless absolute path, defaults to local directory)'}},

'database': {'opts': ('-D', '--database'), 'attrs': {'type': 'string', 'dest': 'database', 'help': 'Database containing formatting parameters of search programs'}},

'comp': {'opts': ('-T', '--comp'), 'attrs': {'type': 'string', 'dest': 'comp', 'help': 'Location of compareSearches.py output (.tdv)'}},

'columns': {'opts': ('-k', '--columns'), 'attrs': {'type': 'string', 'dest': 'columns', 'help': 'Path to pickled python tuple listing columns to include in LADS .tdv output, default is light scan number, heavy scan number (where applicable), [M+H], score, sequence'}},

'verbose': {'opts': ('-v', '--verbose'), 'attrs': {'action': 'store_true', 'default': False, 'dest': 'verbose', 'help': 'Print status updates to std out'}},

'paircutoff': {'opts': ('-P', '--paircutoff'), 'attrs': {'type': 'float', 'dest': 'paircutoff', 'default': 0.05, 'help': 'Cutoff shared peaks ratio when determining paired spectra'}},

#'Nmod': {'opts': ('-N', '--Nmod'), 'attrs': {'type': 'float', 'dest': 'Nmod', 'default': 0, 'help': 'N-terminal mass difference between paired spectra'}},

#'Cmod': {'opts': ('-C', '--Cmod'), 'attrs': {'type': 'float', 'dest': 'Cmod', 'default': Constants.mods['*'], 'help': 'C-terminal mass difference between paired spectra'}},

'ppmstd': {'opts': ('-p', '--ppmstd'), 'attrs': {'type': 'float', 'dest': 'ppmstd', 'default': 5, 'help': 'Expected standard deviation in ppm error distributions of spectrum graph edges'}},

'ppmsyserror': {'opts': ('-M', '--ppmsyserror'), 'attrs': {'type': 'float', 'dest': 'ppmsyserror', 'default': 0, 'help': 'Expected systematic error in the ppm error of the edges from true sequence masses'}},

'ambigpenalty': {'opts': ('-a', '--ambigpenalty'), 'attrs': {'type': 'float', 'dest': 'ambigpenalty', 'default': 20, 'help': 'Score penalty for ambiguous edges in spectrum graph'}},

'ppmpenalty': {'opts': ('-A', '--ppmpenalty'), 'attrs': {'type': 'float', 'dest': 'ppmpenalty', 'default': 20, 'help': 'Maximum score penalty for ppm deviation from true sequence mass in edges of spectrm graph'}},

'minedge': {'opts': ('-e', '--minedge'), 'attrs': {'type': 'float', 'dest': 'minedge', 'default': 300, 'help': 'Minimum length of ambiguous edges in spectrum graph'}},

'maxedge': {'opts': ('-E', '--maxedge'), 'attrs': {'type': 'float', 'dest': 'maxedge', 'default': 500, 'help': 'Maximum length of ambiguous edges in spectrum graph'}},

'alpha': {'opts': ('-l', '--alpha'), 'attrs': {'type': 'float', 'dest': 'alpha', 'default': 0.90, 'help': 'Minimum score ratio between suboptimal and optimal paths'}},

#'matrix': {'opts': ('-M', '--matrix'), 'attrs': {'type': 'string', 'dest': 'matrix', 'help': 'Location of scoring matrix for sequence alignment'}},

#'gapopen': {'opts': ('-g', '--gapopen'), 'attrs': {'type': 'int', 'dest': 'gapopen', 'default':-4, 'help': 'gap initiation penalty for sequence alignment'}},

#'gapextend': {'opts': ('-G', '--gapextend'), 'attrs': {'type': 'int', 'dest': 'gapextend', 'default':-2, 'help': 'gap extension penalty for sequence alignment'}},

'lads': {'opts': ('-L', '--lads'), 'attrs': {'type': 'string', 'dest': 'lads', 'help': 'Dictionary mapping location of LADS search results (.tdv) to its name in the output'}},

'sequest':{'opts': ('-s', '--sequest'), 'attrs': {'type': 'string', 'dest': 'sequest', 'help': 'Dictionary mapping location of SEQUEST search results (.csv) to its name in the output'}},

'mascot': {'opts': ('-t', '--mascot'), 'attrs': {'type': 'string', 'dest': 'mascot', 'help': 'Dictionary mapping location of MASCOT search results (.csv) to its name in the output'}},

'pepnovo': {'opts': ('-O', '--pepnovo'), 'attrs': {'type': 'string', 'dest': 'pepnovo', 'help': 'Location of pepNOVO search results (.tdv)'}},

'pnovo': {'opts': ('-N', '--pnovo'), 'attrs': {'type': 'string', 'dest': 'pnovo', 'help': 'Location of pNovo search results (.tdv)'}},

'peaks': {'opts': ('-K', '--peaks'), 'attrs': {'type': 'string', 'dest': 'peaks', 'help': 'Location of PEAKS search results (.tdv)'}},

'pairfinder': {'opts': ('-f', '--pairfinder'), 'attrs': {'type': 'string', 'dest': 'pairfinder', 'help': 'Function used to find heavy light pairs. String input can be either module.function or function (module defaults to Analytics). Functon must take arguments of form func(lightspec, heavyspec, nmod, cmod, epsilon)'}},

'testvals':  {'opts': ('-V', '--testvals'), 'attrs': {'type': 'string', 'dest': 'testvals', 'help': 'Input of the form a,b,c, etc. which evals to a tuple. Used to test a program over multiple values of a parameter'}},

'progdict': {'opts': ('-r', '--progdict'), 'attrs': {'type': 'string', 'dest': 'progdict', 'help': 'Dictionary which maps names of searches in CompareSearches.py output to the program which generated them'}},

'denovoscript': {'opts': ('-n', '--denovoscript'), 'attrs': {'type': 'string', 'dest': 'denovoscript', 'help': 'name of de novo program script to unit test'}},

'mainprogname': {'opts': ('-R', '--mainprogname'), 'attrs': {'type': 'string', 'dest': 'mainprogname', 'help': 'name of output to compare against others when analyzing compareSearches.py output'}},

'subgraphcut': {'opts': ('-u', '--subgraphcut'), 'attrs': {'type': 'float', 'dest': 'subgraphcut', 'default': 300, 'help': 'maximum edge length to be taken for complete sequencing (as opposed to creating a sub-spectrum graph) in resolving ambiguous edges'}},

'number': {'opts': ('-b', '--number'), 'attrs': {'type': 'float', 'dest': 'number', 'help': 'Just a number, use varies with program'}},

'init': {'opts': ('-i', '--init'), 'attrs': {'type': 'string', 'dest': 'init', 'help': 'Initiation file used to configure LADS'}},

'symbolmap': {'opts': ('-S', '--symbolmap'), 'attrs': {'type': 'string', 'dest': 'symbolmap', 'help': 'Key mapping program symbols from search program output the their respective modification types or amino acids(path to pickled python dict).'}},

'modalpha': {'opts': ('-q', '--modalpha'), 'attrs': {'type': 'float', 'dest': 'modalpha', 'help': 'suboptimal pathing cutoff alpha for modification amino acids'}},

'unimoddict': {'opts': ('-Q', '--unimoddict'), 'attrs': {'type': 'string', 'dest': 'unimoddict', 'help': 'location of pickled unimod dictionary'}},

'combined': {'opts': ('-C', '--combined'), 'attrs': {'type': 'string', 'dest': 'combined', 'help': 'Dictionary which maps location of combined sequest/mascot results to their name in the output'}},

'srchid': {'opts': ('-H', '--srchid'), 'attrs': {'type': 'int', 'dest': 'srchid', 'default': None, 'help': 'Dictionary which maps name of database search in output to the searchID of interest in aggregated exported database search file (.csv), optional'}}
}

def parse(arglist, optArgs=[]):
    parser = OptionParser(usage="%prog [options]")
        
    for arg in arglist:
        parser.add_option(*A[arg]['opts'], **A[arg]['attrs'])

    for optArg in optArgs:
        parser.add_option(*optArg['opts'], **optArg['attrs'])
        
    (options, args) = parser.parse_args()

    return options

def parseInitFile(init, options):
    paramsDict = DataFile.parseParams(init)
    for param in paramsDict['LADS Parameters'].keys():
        try:
            paramType = A[param]['attrs']['type']
            val = init['LADS Parameters'][param]
            if paramType != 'string':
                val = getattr('__builtin__', paramType)(val)
            
            setattr(options, param, val)
        except KeyError:
            pass
    
    return paramsDict

def getProgDict(progs, options, progDict={}):
    for prog in progs:
        if hasattr(options, prog.lower()) and getattr(options, prog.lower()) != None:
            fDict = eval(getattr(options, prog.lower()))
            for name in fDict.values():
                progDict[name] = prog
    
    return progDict
        
