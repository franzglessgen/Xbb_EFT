#!/usr/bin/env python
import ROOT
from WeightInfo                   import WeightInfo
from HyperPoly                    import HyperPoly

#from ROOT import TVector3, TLorentzVector
#from math import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
import numpy as np
#import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class EFT_params(AddCollectionsModule):


    def __init__(self, branchName='EFT_params'):
        super(EFT_params, self).__init__()
        self.branchName = branchName


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
       
        self.addVectorBranch(self.branchName + '_weight', length=255)
        self.addVectorBranch(self.branchName + '_coeff', length=255)
        self.addIntegerBranch(self.branchName + '_np')
        self.addBranch(self.branchName + '_chi2_ndof')
        self.addIntegerBranch(self.branchName + '_nWC')

        reweight_pkl = "/work/fglessge/CMSSW_10_1_0/src/Xbb_EFT/python/myutils/SMEFTsim_VH_reweight_card.pkl"        # this has the information on WCs written in the reweight card (particularly, which WCs are used and in what combination)
        weightInfo = WeightInfo( reweight_pkl )
        ref_point_coordinates = [ weightInfo.ref_point_coordinates[var] for var in weightInfo.variables ]

        interpolationOrder = int(2)                              # because we're using a second-order polynomial
        self.hyperPoly       = HyperPoly( interpolationOrder )

        weightInfo_data = list(weightInfo.data.iteritems())
        weightInfo_data.sort( key = lambda w: w[1] )
        basepoint_coordinates = map( lambda d: [d[v] for v in weightInfo.variables] , map( lambda w: self.interpret_weight(w[0]), weightInfo_data) )

        self.hyperPoly.initialize( basepoint_coordinates, ref_point_coordinates )

        self.sampleTree = initVars['sampleTree']


    def processEvent(self, tree):


        
	# if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b(self.branchName + '_weight')[0] = -99997.0
            self._b(self.branchName + '_coeff')[0] = -99997.0
            self._b(self.branchName + '_np')[0] = -99997
            self._b(self.branchName + '_chi2_ndof')[0] = -99997.0
            self._b(self.branchName + '_nWC')[0] = -99997

            #fl = ROOT.TFile.Open("NanoAOD_ZH.root")
            #tree = fl.Get("Events")
            #tree.GetEntry(0)

            weights  = [ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]
            


	    coeff     = self.hyperPoly.get_parametrization( weights )
            np        = self.hyperPoly.ndof                                 # no. of coefficients in polynomial
            chi2_ndof = self.hyperPoly.chi2_ndof( coeff, weights )          # chi^2/ndof of the polynomial fit
            nWC       = self.hyperPoly.nvar                                 # no. of Wilson coefficients

            
	    for i in range(np): self._b(self.branchName + '_coeff')[i] = coeff[i]
            self._b(self.branchName + '_np')[0] = np
            self._b(self.branchName + '_chi2_ndof')[0] = chi2_ndof
            self._b(self.branchName + '_nWC')[0] = nWC
            
            WCs = []
            for ic in range(nWC):
                WCs.append(0)
            for n in range(len(WCs)):
                WCs[n] = 1
                self._b(self.branchName + '_weight')[n] = self.hyperPoly.eval(coeff, *WCs)
		WCs[n] = 0


    def interpret_weight(self, weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res


if __name__ == '__main__':


	reweight_pkl = "/work/fglessge/CMSSW_10_1_0/src/Xbb_EFT/python/myutils/SMEFTsim_VH_reweight_card.pkl"        # this has the information on WCs written in the reweight card (particularly, which WCs are used and in what combination)
        weightInfo = WeightInfo( reweight_pkl )
        ref_point_coordinates = [ weightInfo.ref_point_coordinates[var] for var in weightInfo.variables ]


	print(ref_point_coordinates)
        weightInfo_data = list(weightInfo.data.iteritems())
        
	print(weightInfo_data[0][1])
	print(weightInfo_data[0])

	
	weightInfo_data.sort( key = lambda w: w[1] )

	for i in range(len(weightInfo_data)):
		print(weightInfo_data[i][1])
        
	
    	def Interpret_weight(weight_id):
		str_s = weight_id.split('_')
		res={}
		for i in range(len(str_s)/2):
		    res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
		return res
	

	#print(weightInfo_data[0])
	#print(Interpret_weight(weightInfo_data[0][0]))
        
	basepoint_coordinates = map( lambda d: [d[v] for v in weightInfo.variables] , map( lambda w: Interpret_weight(w[0]), weightInfo_data) )

	print(basepoint_coordinates)
	print(ref_point_coordinates)

#Examples to get weight for few combinations of WCs are below

#for ieft in range(nWC):
#    WCs = []
#    for ic in range(hyperPoly.nvar):
#        WCs.append(0)
#    if(ieft>0):
#        WCs[ieft-1] = 1.
#    EFT_wt =  hyperPoly.eval(coeff,*WCs)
#    print("Weight id:",ieft+1,"weight:",EFT_wt)

