#!/usr/bin/env python
import ROOT
import sys
import numpy as np
import datetime
import os

# RUN:
# ./SignalEwkQcdCorrection.py input.root output.root applyCorrection(bool)
#

# **********************************************************************************************************************
#  Class to apply corrections text file
# **********************************************************************************************************************
class TreeCopierWithCorrectionFromFile:
    def __init__(self):
        self.correctionFiles = {}
        self.correctionFileIndexColumns = {}
        self.delimiter = '  '
        self.commentsMarker = '#'
        self.outputBranches = []
        self.corrections = {}
        self.inputBranchName = None

    def loadCorrectionFiles(self, corrections):
        for correction in corrections:
            with open(correction['file'], 'rb') as correctionFile:
                self.correctionFiles[correction['name']] = np.genfromtxt(correctionFile,
                                 delimiter=correction['delimiter'] if 'delimiter' in correction else self.delimiter,
                                 comments=correction['comments'] if 'comments' in correction else self.commentsMarker)
                self.correctionFileIndexColumns[correction['name']] = self.correctionFiles[correction['name']][:, 0]

    def setInputBranch(self, inputBranchName):
        self.inputBranchName = inputBranchName

    def addOutputBranch(self, branchName, branchFormula):
        self.outputBranches.append({'name': branchName, 'formula': branchFormula})

    def getCorrection(self, v):
        corr = {}

        # select corresponding line from the table with correct Vpt
        for correctionName, correctionFile in self.correctionFiles.iteritems():
            index = self.correctionFileIndexColumns[correctionName].searchsorted(v)
            index -= 1
            if index < 0:
                index = 0
            if index > len(correctionFile) - 1:
                index = len(correctionFile) - 1
            corr[correctionName] = correctionFile[index]

        return corr

    # ------------------------------------------------------------------------------------------------------------------
    # copy the tree and add new branch with weight from text file
    # ------------------------------------------------------------------------------------------------------------------
    def copy(self, inputFileName, outputFileName, applyNNLO):

        ifile = ROOT.TFile.Open(inputFileName, "READ")
        ofile = ROOT.TFile.Open(outputFileName, "RECREATE")

        # copy other file objects
        ifile.cd()
        obj = ROOT.TObject
        for key in ROOT.gDirectory.GetListOfKeys():
            ifile.cd()
            obj = key.ReadObj()
            if obj.GetName() == 'tree':
                continue
            ofile.cd()
            obj.Write(key.GetName())

        ifile.cd()
        tree = ifile.Get('tree')
        ofile.cd()

        # clone TTree
        otree = tree.CloneTree(0)

        # --------------------------------------------------------------------------------------------------------------
        # set branches
        # --------------------------------------------------------------------------------------------------------------

        # set branch addresses for input variable (1D at the moment!)
        inputBranchValue = np.array([0], dtype='f')
        tree.SetBranchAddress(self.inputBranchName, inputBranchValue)

        # define additional output branches
        outputBranchValues = {}
        for outputBranch in self.outputBranches:
            outputBranchValues[outputBranch['name']] = np.array([0], dtype=float)
            obranch = otree.Branch(outputBranch['name'], outputBranchValues[outputBranch['name']],
                                   outputBranch['name'] + '/D')

        # --------------------------------------------------------------------------------------------------------------
        # apply corrections
        # --------------------------------------------------------------------------------------------------------------
        nEntries = tree.GetEntries()
        print "nEntries = ", nEntries
        for entry in range(nEntries):
            if entry % 10000 == 0:
                print "processing entry: %d"%entry

            tree.GetEntry(entry)

            # add all output branches
            for outputBranch in self.outputBranches:
                corrData = self.getCorrection(inputBranchValue[0])
                outputBranchValues[outputBranch['name']][0] = outputBranch['formula'](corrData) if applyNNLO else 1.0

            otree.Fill()

        # --------------------------------------------------------------------------------------------------------------
        # write and close files
        # --------------------------------------------------------------------------------------------------------------
        print "writing to file..."
        ofile.cd()
        otree.Write()

        ofile.Close()
        ifile.Close()

        print "done."

    # prepare the file with multiplicative weights
    def prepare(self, silent=False, inputFileName = None, vptHistogramFileName = None):
        outputTable = []
        histograms = {}
        histogramsCorr = {}
        ifileVpt = None

        xbbPath = '/'.join(os.path.realpath(__file__).split('/')[:-3]) + '/'
        tempHistogramName = xbbPath + 'python/weights/Zll_temp_histograms.root'
        print "creating temporary .root file:", tempHistogramName
        rootFile = ROOT.TFile(tempHistogramName, 'RECREATE')


        for outputBranch in self.outputBranches:
            histograms[outputBranch['name']] = ROOT.TH1D(outputBranch['name'], outputBranch['name'], 100, 0.0, 500.0)
            histogramsCorr[outputBranch['name']] = ROOT.TH1D("Vpt_corr_" + outputBranch['name'], "Vpt_corr_" + outputBranch['name'], 100, 0.0, 500.0)

        # create a 1D function with the corrected cross section (saved in TH1D)
        #  and some debug output
        for Vpt in range(0, 500, 5):
            outputString = "%d "%Vpt
            for outputBranch in self.outputBranches:
                if not silent:
                    print "Vpt:",Vpt
                corrData = self.getCorrection(Vpt)
                if not silent:
                    print "corrdata:", corrData
                corrValue = outputBranch['formula'](corrData)
                histograms[outputBranch['name']].Fill(Vpt, corrValue)
                outputString += "%0.7f "%corrValue
            if not silent:
                outputTable.append(outputString)

        if inputFileName:

            VptHistogram = None

            # read Vpt distribution as a histogram from .root file
            if vptHistogramFileName:
                ifileVpt = ROOT.TFile.Open(vptHistogramFileName, "READ")

                # read Vpt distribution as a histogram
                VptHistogram = ifileVpt.Get('VptIncl_AjetIncl__V_pT')
                if VptHistogram:
                    VptHistogram.SetDirectory(0)
                    print "Vpt distribution for normalization is taken from HISTO:", vptHistogramFileName

            # signal sample Vpt histogram
            ifile = ROOT.TFile.Open(inputFileName, "READ")
            VptHistogramSignal = ROOT.TH1D("Vpt", "Vpt", 100, 0.0, 500.0)
            VptHistogramSignal.SetDirectory(0)

            tree = ifile.Get('tree')
            inputBranchValue = np.array([0], dtype='f')
            tree.SetBranchAddress("V_pt", inputBranchValue)
            nEntries = tree.GetEntries()
            print "nEntries = ", nEntries
            for entry in range(nEntries):
                if entry % 10000 == 0:
                    print "processing entry: %d" % entry
                tree.GetEntry(entry)
                VptHistogramSignal.Fill(inputBranchValue[0])
            ifile.Close()

            # or fill it from tree if it is not found
            if not VptHistogram:
                print "Vpt distribution for normalization is taken from TREE:", inputFileName
                VptHistogram = VptHistogramSignal.Clone("Vpt_SIGNAL")

            # compute relative corrections by dividing "nnloQCD * (1 + delta_EWK)" over "powheg"
            ratioHistograms = {}
            for outputBranch in self.outputBranches:
                print "add ratio:", outputBranch['name']
                ratioHistograms[outputBranch['name']] = histograms[outputBranch['name']].Clone("ratio_" + outputBranch['name'])
                ratioHistograms[outputBranch['name']].Divide(VptHistogram)
                ratioHistograms[outputBranch['name']].SetDirectory(rootFile)

            # apply relative corrections (normalization is still wrong)
            ifile = ROOT.TFile.Open(inputFileName, "READ")
            tree = ifile.Get('tree')
            inputBranchValue = np.array([0], dtype='f')
            tree.SetBranchAddress("V_pt", inputBranchValue)
            nEntries = tree.GetEntries()
            print "nEntries = ", nEntries
            for entry in range(nEntries):
                if entry % 10000 == 0:
                    print "processing entry: %d" % entry

                tree.GetEntry(entry)
                for outputBranch in self.outputBranches:
                    histogramsCorr[outputBranch['name']].Fill(inputBranchValue[0], ratioHistograms[outputBranch['name']].GetBinContent(ratioHistograms[outputBranch['name']].GetXaxis().FindBin(inputBranchValue[0])))

            # compute normalization and apply it to histograms
            #  take Vpt histogram from SIGNAL here
            ifile.Close()
            for outputBranch in self.outputBranches:
                scaleFactor = 1.0 * VptHistogramSignal.Integral() / histogramsCorr[outputBranch['name']].Integral()
                print "normalization factor:", scaleFactor
                histogramsCorr[outputBranch['name']].Scale(scaleFactor)
                ratioHistograms[outputBranch['name']].Scale(scaleFactor)

        # write text file with multiplicative and normalized weights
        for outputBranch in self.outputBranches:
            outputFileName = xbbPath + 'python/weights/Zll_weight_' + outputBranch['name'] + '.dat'
            with open(outputFileName, 'w') as outfile:
                outfile.write("# this file has been automatically created on " + str(datetime.datetime.now()) + "\n")
                outfile.write("# command: %s\n"%(' '.join(sys.argv)))
                outfile.write("# V_pt  correction\n")
                for i in range(1, ratioHistograms[outputBranch['name']].GetXaxis().GetNbins()+1):
                    outfile.write("%f %f\n"%(ratioHistograms[outputBranch['name']].GetBinLowEdge(i),ratioHistograms[outputBranch['name']].GetBinContent(i)))
            print "weights written to:", outputFileName

        rootFile.Write()
        rootFile.Close()
        if ifileVpt:
            ifileVpt.Close()

# **********************************************************************************************************************
#  main
# **********************************************************************************************************************

try:
    xbbPath = '/'.join(os.path.realpath(__file__).split('/')[:-3]) + '/'
except:
    print "ERROR: CAN'T GET XBB PATH!!!"
    sys.exit(9)

theTreeCopier = TreeCopierWithCorrectionFromFile()
theTreeCopier.setInputBranch('V_pt')

# --------------------------------------------------------------------------------------------------------------
# prepare: reads the (N)NLO cross sections from the files and computes the normalization factors.
# It writes a file of multiplicative weights in the end.
# --------------------------------------------------------------------------------------------------------------
if len(sys.argv) >= 3 and sys.argv[1] == 'prepare':
    print "running prepare:"

    # define text files which include the corrections
    correctionFiles = [
        {
            'name': 'deltaEWK',
            'file': xbbPath + 'python/weights/Zll/distributions/dat.ptv',
            'delimiter': '  '  # two spaces here!
        },
        {
            'name': 'nnloQCD',
            'file': xbbPath + 'python/weights/Zll/nnlo_qcd/ze.ptv',
            'delimiter': ' '  # just one space here!
        }
    ]
    theTreeCopier.loadCorrectionFiles(correctionFiles)

    # define new branches
    sigmaTotal = 0.10155859E+02  # doesn't really matter, will be normalized again later
    COLUMN_EWK_NLO = 5
    COLUMN_EWK_BORN = 1
    COLUMN_QCD_NNLO = 2

    # see https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGVHYR4
    theTreeCopier.addOutputBranch(
        branchName='SignalCS_nloEWK_nnloQCD',
        branchFormula=lambda x: sigmaTotal * x['nnloQCD'][COLUMN_QCD_NNLO] * (
            1.0 + x['deltaEWK'][COLUMN_EWK_NLO] / x['deltaEWK'][COLUMN_EWK_BORN])
    )

    theTreeCopier.prepare(silent=False, inputFileName=sys.argv[2], vptHistogramFileName=sys.argv[3] if len(sys.argv) > 3 else None)

# --------------------------------------------------------------------------------------------------------------
# add a branch containing the multiplicative weights calculated above
# --------------------------------------------------------------------------------------------------------------
elif len(sys.argv) == 4:
    applyNNLO = sys.argv[3].lower().strip() == 'true'

    # define text files which include the corrections
    correctionFiles = [
        {
            'name': 'signalWeight',
            'file': xbbPath + 'python/weights/Zll_weight_SignalCS_nloEWK_nnloQCD.dat',
            'delimiter': ' '
        }
    ]
    theTreeCopier.loadCorrectionFiles(correctionFiles)

    theTreeCopier.addOutputBranch(
        branchName='SignalWeight_nloEWK_nnloQCD',
        branchFormula=lambda x: x['signalWeight'][1]
    )

    # now run the tree copier again and fill the correctly normalized weight
    theTreeCopier.copy(inputFileName=sys.argv[1], outputFileName=sys.argv[2], applyNNLO=applyNNLO)

else:
    print "syntax(1): {file} prepare full_phase_space_input.root".format(file=sys.argv[0] if len(sys.argv) > 0 else './file.py')
    print "syntax(2): {file} input.root output.root useNNLO(bool)".format(file=sys.argv[0] if len(sys.argv) > 0 else './file.py')

