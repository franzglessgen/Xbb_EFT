#!/usr/bin/env python
import ROOT
from sklearn.preprocessing import PolynomialFeatures

#from ROOT import TVector3, TLorentzVector
#from math import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
import numpy as np
#import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class EFT_interpolation(AddCollectionsModule):


    def __init__(self, branchName='EFT_interpolation'):
        super(EFT_interpolation, self).__init__()
        self.branchName = branchName


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
       
        self.addVectorBranch(self.branchName + '_weight', length=255)
        #self.addVectorBranch(self.branchName + '_coeff', length=255)
        #self.addIntegerBranch(self.branchName + '_np')
        #self.addBranch(self.branchName + '_chi2_ndof')
        self.addIntegerBranch(self.branchName + '_nWC')


	self.nWC = int(self.config.get('WCGeneral', 'nbofWC'))
	self.ndof = 2*self.nWC + self.nWC*(self.nWC -1)/2 + 1        

	#This hold the inverse Vandermonde matrix for computation of the interpolation
	self.VandermondePath = self.config.get('WCGeneral', 'SimWC')
	self.invVandermonde = np.load(self.VandermondePath)
	self.polyfeatures = PolynomialFeatures(2)		
	self.WCtoeval = self.polyfeatures.fit_transform(np.diag(np.ones(self.nWC)))

	#print(">>>>>>>>>>>>>>>>>>>>>>>> PARAMS : ")

        self.sampleTree = initVars['sampleTree']

    def processEvent(self, tree):



        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b(self.branchName + '_weight')[0] = -99997.0
            #self._b(self.branchName + '_coeff')[0] = -99997.0
            #self._b(self.branchName + '_np')[0] = -99997
            #self._b(self.branchName + '_chi2_ndof')[0] = -99997.0
            self._b(self.branchName + '_nWC')[0] = -99997

            #fl = ROOT.TFile.Open("NanoAOD_ZH.root")
            #tree = fl.Get("Events")
            #tree.GetEntry(0)

            weights  = np.array([ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]).reshape(tree.nLHEReweightingWeight,1)
            


	    coeff     = np.matmul(self.invVandermonde, weights).flatten()
	
           
	    """ 
	    for i in range(self.ndof): self._b(self.branchName + '_coeff')[i] = coeff[i]
            self._b(self.branchName + '_np')[0] = self.ndof
	    """           
            self._b(self.branchName + '_nWC')[0] = self.nWC
 
	    EvaluatedWC = np.matmul(self.WCtoeval, coeff).flatten()

            for n in range(self.nWC):
                self._b(self.branchName + '_weight')[n] = EvaluatedWC[n]


    def interpret_weight(self, weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res


if __name__ == '__main__':

	print("Interpolation module")

