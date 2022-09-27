#! /usr/bin/env python
import os, pickle, sys, ROOT
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser
from myutils import BetterConfigParser, copytree, copytreePSI, ParseInfo
import utils
import glob
import os
argv = sys.argv

#------------------------------------------------------------------------------
# prepare sample .txt files from ntuple folders
#------------------------------------------------------------------------------

#get files info from config
parser = OptionParser()
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="directory config")
(opts, args) = parser.parse_args(argv)
config = BetterConfigParser()
config.read(opts.config)

inputPath = config.get('Directories', 'samplepath')
outputPath = config.get('Directories', 'samplefiles')
print 'INPUT: ', inputPath
print 'OUTPUT:', outputPath
print '-'*80
if '/pnfs/' in inputPath:
    inputPath = '/pnfs/' + inputPath.split('/pnfs/')[-1]

sampleNameReplacementRules = {
        'RunIISummer16MiniAODv2_': '',
        }
sampleFileReplacementRules = {
        '/pnfs/psi.ch/cms/trivcat': ''
        }

# find all folders which contain ntuples
inputTextFileNames = glob.glob(inputPath + '/*/')
for inputTextFileName in inputTextFileNames:
    sampleFolderName = inputTextFileName.strip('/').split('/')[-1]
    print '-> '+sampleFolderName
    sampleName = sampleFolderName
    for k,v in sampleNameReplacementRules.iteritems():
        sampleName = sampleName.replace(k, v)
    print '  -> name: ' + sampleName

    # find all .root trees in this folder
    sampleFiles = glob.glob(inputPath + '/{sampleFolderName}/*/*/*.root'.format(sampleFolderName=sampleFolderName))
    print '  -> files:'
    treeNames = []
    filesToWrite = []
    for sampleFile in sampleFiles:
        sampleFileFinal = sampleFile
        for k,v in sampleFileReplacementRules.iteritems():
            sampleFileFinal = sampleFileFinal.replace(k, v)
        print '      ',sampleFileFinal
        treeName = sampleFileFinal.split('/')[-2] + '/' + sampleFileFinal.split('/')[-1]
        if treeName in treeNames:
            print '\x1b[31mERROR: duplicate tree!!\x1b[0m'
        else:
            treeNames.append(treeName)
            filesToWrite.append(sampleFileFinal)

    # write output .txt file
    outputFileName = outputPath + '/' + sampleName + '.txt'
    with open(outputFileName, 'w') as outputFile:
        for fileToWrite in filesToWrite:
            outputFile.write(fileToWrite.strip() + '\n')
    print ' => wrote ' + outputFileName
