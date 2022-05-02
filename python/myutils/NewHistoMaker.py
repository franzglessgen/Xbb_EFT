from __future__ import print_function
import ROOT
import TdrStyles
import os
import sys
import array

from NewTreeCache import TreeCache as TreeCache
from sampleTree import SampleTree

class NewHistoMaker:

    instanceCounter = 0

    def __init__(self, config, sample, sampleTree, histogramOptions):
        NewHistoMaker.instanceCounter += 1
        self.debug = 'XBBDEBUG' in os.environ
        self.config = config
        self.sample = sample
        self.sampleTree = sampleTree
        self.histogramOptions = histogramOptions
        self.histogram = None
        self.EvalCut = config.get('Cuts', 'EvalCut')
        self.originalBins = []
        TdrStyles.tdrStyle()

    def create1Dhistogram(self, histogramName, histogramOptions):
        if 'binList' in histogramOptions:
            print("DEBUG: creating histogram with variable bins:", histogramOptions['binList'])
            binArray = array.array('d', histogramOptions['binList']) 
            return ROOT.TH1F(histogramName, histogramName, len(histogramOptions['binList']) - 1, binArray)
        else:
            print("DEBUG: creating histogram with equidistant bins:", histogramOptions['nBinsX'], histogramOptions['minX'], histogramOptions['maxX'])
            return ROOT.TH1F(histogramName, histogramName, histogramOptions['nBinsX'], histogramOptions['minX'], histogramOptions['maxX'])

    def create2Dhistogram(self, histogramName, histogramOptions):
        return ROOT.TH2F(histogramName, histogramName, histogramOptions['nBinsX'], histogramOptions['minX'], histogramOptions['maxX'], histogramOptions['nBinsY'], histogramOptions['minY'], histogramOptions['maxY'])

    def createProfilehistogram(self, histogramName, histogramOptions):
        return ROOT.TProfile(histogramName, histogramName, histogramOptions['nBins'], histogramOptions['minX'], histogramOptions['maxX'])
    

    def getHistogramName(self):
        histogramName = (self.histogramOptions['name']+'_') if 'name' in self.histogramOptions else ''
        histogramName += self.sample.name + ('_' + self.histogramOptions['var'] if 'var' in self.histogramOptions else '')

        # add unique instance counter to avoid same name for histogram and ROOT complaining
        if 'uniqueid' in self.histogramOptions and self.histogramOptions['uniqueid']:
            histogramName += '_instance%d'%NewHistoMaker.instanceCounter
        return histogramName

    def initializeHistogram(self):
        self.histogramName = self.getHistogramName()

        is2D = ':' in self.histogramOptions['treeVar'].replace('::', '')
        isProfile = 'customProfile' in self.histogramOptions and self.histogramOptions['customProfile']

	#self.histogram = self.create2Dhistogram(self.histogramName, self.histogramOptions) if is2D else self.create1Dhistogram(self.histogramName, self.histogramOptions)
       
	if is2D:
		self.histogram = self.create2Dhistogram(self.histogramName, self.histogramOptions) 
	elif isProfile: 
		self.histogram = self.createProfilehistogram(self.histogramName, self.histogramOptions) 
	else:
		self.histogram = self.create1Dhistogram(self.histogramName, self.histogramOptions) 
		

        self.histogram.Sumw2()
        self.histogram.SetTitle(self.sample.name)

        return self.histogram

    def scaleHistogram(self):
        if self.sample.type != 'DATA':
            ScaleFactor = self.sampleTree.getScale(self.sample)
            if ScaleFactor != 0:
                self.histogram.Scale(ScaleFactor)
            else:
                print ("\x1b[31mWARNING: histogram scaling factor is 0!\x1b[0m")

    def getHistogram(self, cut='1'):
        if self.initializeHistogram():
            # apply weights only to MC and not to DATA
            if 'group' in self.histogramOptions and self.histogramOptions['group'] == 'DATA':
                weightF = '1'
            else:
                weightF = "({weight})".format(weight=self.histogramOptions['weight'] if ('weight' in self.histogramOptions and self.histogramOptions['weight']) else '1') 

            # drop training events and rescale MC by 2 for BDT plots
            if ('BDT' in self.histogramOptions['treeVar'] or 'DNN' in self.histogramOptions['treeVar']) and self.sample.type != 'DATA':
                cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=self.EvalCut)
                weightF = '(({weight1})*({weight2}))'.format(weight1=weightF, weight2='2.0')
                print("INFO: training events removed for \x1b[32m", self.histogramOptions['treeVar'], "\x1b[0m plot with additional cut \x1b[35m", self.EvalCut, "\x1b[0m, MC rescaled by \x1b[36m2.0\x1b[0m")
            # add blind cut on data if specified for this variable
            if 'blindCut' in self.histogramOptions and self.histogramOptions['group'] == 'DATA':
                cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=self.histogramOptions['blindCut'])
                #print('cut is', cut)
                #sys.exit()

            # per sample special weight
            if self.config.has_option('Weights', 'useSpecialWeight') and eval(self.config.get('Weights', 'useSpecialWeight')):
                specialweight = self.sample.specialweight
                weightF = "(({weight})*({specialweight}))".format(weight=weightF, specialweight=specialweight)
                print ("INFO: use specialweight: {specialweight}".format(specialweight=specialweight))

            # (deprecated) additional blinding cut applied to DATA only in BDT plots
            if 'BDT' in self.histogramOptions['treeVar'] and self.sample.type == 'DATA':
                if self.config.has_section('Blinding') and self.config.has_option('Blinding', 'BlindBDTinSR'):
                    self.blindingCut = self.config.get('Blinding', 'BlindBDTinSR')
                    cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=self.blindingCut)
                    print ("\x1b[31mUSE BLINDING CUT FOR DATA in BDT!!", self.blindingCut ,"\x1b[0m")

            # (deprecated) additional mass blinding
            if 'mass' in self.histogramOptions['treeVar'] and self.sample.type == 'DATA':
                if self.config.has_section('Blinding') and self.config.has_option('Blinding', 'BlindDataInMassPlots') and eval(self.config.get('Blinding', 'BlindDataInMassPlots')):
                    self.blindingCut = '0'
                    cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=self.blindingCut)
                    print ("\x1b[31mUSE BLINDING CUT FOR DATA in BDT!!", self.blindingCut ,"\x1b[0m")

            # region/var specific blinding cut
            if 'blindCut' in self.histogramOptions and self.sample.type == 'DATA':
                cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=self.histogramOptions['blindCut'])
                print("INFO: found region/var specific blinding cut in histogramOptions:\n--->{cut}".format(cut=self.histogramOptions['blindCut']))

            # global cut affecting everything!! (but not caching, so can be used to further tighten selection in plots without re-caching)
            if self.config.has_option('Cuts', 'additionalPlottingCut'):
                print("\x1b[31mINFO: there is a global additional cut used (defined in Cuts->additionalPlottingCut)!\x1b[0m")
                globalCut = self.config.get('Cuts', 'additionalPlottingCut')
                print(globalCut)
                cut = '(({cut1})&&({cut2}))'.format(cut1=cut, cut2=globalCut)
            if self.config.has_option('Cuts', 'additionalPlottingWeight'):
                print("\x1b[31mINFO: there is a global additional weight used, ALSO FOR DATA! (defined in Cuts->additionalPlottingWeight)\x1b[0m")
                globalCut = self.config.get('Cuts', 'additionalPlottingWeight')
                weightF = '({cut2})*({cut1})'.format(cut1=weightF, cut2=globalCut)
                print("INFO: ->", weightF)

            # add tree cut 
            if self.config.has_option('Plot_general', 'drawWeight'):
                weightF = '({w})*{s}'.format(w=weightF, s=self.sampleTree.getScale(self.sample))
                weightF_original = weightF
                if self.config.get('Plot_general', 'drawWeight').startswith('log'):
                    weightF = 'TMath::Log10(TMath::Abs({w}))'.format(w=weightF)
                elif self.config.get('Plot_general', 'drawWeight').startswith('ln'):
                    weightF = 'TMath::Log(TMath::Abs({w}))'.format(w=weightF)
                if self.config.get('Plot_general', 'drawWeight').endswith('weighted'):
                    selection = "({weight})*({cut})".format(weight=weightF_original, cut=cut) 
                    nEvents = self.sampleTree.tree.Draw('{var}>>{histogramName}'.format(var=weightF, histogramName=self.histogramName), selection)
                else:
                    nEvents = self.sampleTree.tree.Draw('{var}>>{histogramName}'.format(var=weightF, histogramName=self.histogramName), cut)
            else:
                selection = "({weight})*({cut})".format(weight=weightF, cut=cut) 
 
        	if 'customProfile' in self.histogramOptions and self.histogramOptions['customProfile']:
			print(">>>>>>>>>>>>>>> CUSTOMPROFILE", weightF)

			nEvents = self.sampleTree.tree.Draw('{var}>>{histogramName}'.format(var=self.histogramOptions['treeVar'] + "[]:" + self.histogramOptions['treeVar']+"_index[]", histogramName=self.histogramName), selection, "prof")
			
		

		else:
			nEvents = self.sampleTree.tree.Draw('{var}>>{histogramName}'.format(var=self.histogramOptions['treeVar'], histogramName=self.histogramName), selection)
                	self.scaleHistogram()
            if nEvents < 0:
                print ("\x1b[31mERROR: error in TTree:Draw! returned {nEvents}\x1b[0m".format(nEvents=nEvents))
            if self.debug:
                print(selection, " => # events:", nEvents, "weighted:", self.histogram.Integral())
        else:
            print ("ERROR: initialization of histogram failed!")
            raise Exception("HistoMakerInitializationError")

        if 'addOverFlow' in self.histogramOptions and self.histogramOptions['addOverFlow']:
            uFlow = self.histogram.GetBinContent(0)+self.histogram.GetBinContent(1)
            oFlow = self.histogram.GetBinContent(self.histogram.GetNbinsX()+1)+self.histogram.GetBinContent(self.histogram.GetNbinsX())
            uFlowErr = ROOT.TMath.Sqrt(ROOT.TMath.Power(self.histogram.GetBinError(0),2)+ROOT.TMath.Power(self.histogram.GetBinError(1),2))
            oFlowErr = ROOT.TMath.Sqrt(ROOT.TMath.Power(self.histogram.GetBinError(self.histogram.GetNbinsX()),2)+ROOT.TMath.Power(self.histogram.GetBinError(self.histogram.GetNbinsX()+1),2))
            self.histogram.SetBinContent(1,uFlow)
            self.histogram.SetBinContent(self.histogram.GetNbinsX(),oFlow)
            self.histogram.SetBinError(1,uFlowErr)
            self.histogram.SetBinError(self.histogram.GetNbinsX(),oFlowErr)


        if 'plotEqualSize' in self.histogramOptions and 'binList' in self.histogramOptions:
            self.histogramOptionsEqualBins = self.histogramOptions.copy()
            self.histogramOptionsEqualBins['nBinsX'] = len(self.histogramOptionsEqualBins['binList']) - 1
            del self.histogramOptionsEqualBins['binList']

            equalHistogram = self.create1Dhistogram(self.histogramName+'_eqBins', self.histogramOptionsEqualBins)
            nBinsSource = self.histogram.GetXaxis().GetNbins()
            nBinsDest = equalHistogram.GetXaxis().GetNbins()
            if nBinsSource != nBinsDest:
                print("ERROR: bin list not compatiable!!", nBinsSource, nBinsDest)
                raise Exception("HistogramDefinitionError")
            for i in range(1, nBinsSource+1):
                self.originalBins.append([self.histogram.GetXaxis().GetBinLowEdge(i), self.histogram.GetXaxis().GetBinUpEdge(i)])
                equalHistogram.SetBinContent(i, self.histogram.GetBinContent(i))
                equalHistogram.SetBinError(i, self.histogram.GetBinError(i))
            self.histogram = equalHistogram

        return self.histogram

