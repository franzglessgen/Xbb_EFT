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

class EFT_correlation(AddCollectionsModule):


    def __init__(self, branchName='EFT_correlation'):
        super(EFT_correlation, self).__init__()
        self.branchName = branchName


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
      

	self.nWC = int(self.config.get('WCGeneral', 'nbofWC'))
	self.ndof = 2*self.nWC + self.nWC*(self.nWC -1)/2 + 1        

	#This hold the inverse Vandermonde matrix for computation of the interpolation
	self.VandermondePath = self.config.get('WCGeneral', 'SimWC')
	self.invVandermonde = np.load(self.VandermondePath)
	self.polyfeatures = PolynomialFeatures(2)		
			
	self.WCindices = self.config.get('Correlation', 'WCindices').split(",")
	self.nWCtocorrelate = len(self.WCindices)	

	#self.RangeStep = int(self.config.get('Scaling', 'NbofPoints'))
        #self.ScalingRange = np.linspace(int(self.ScalingParam[0]), int(self.ScalingParam[1]), self.RangeStep) 
	self.NbEvents = 0
        self.sampleTree = initVars['sampleTree']
	self.WCtoeval = 0

	self.branchsuffix = ""
	for i in range(self.nWCtocorrelate):
		self.branchsuffix+="c" + self.WCindices[i] + "_"

	print("BRANCHNAME", self.branchsuffix)

        #self.addBranch(self.branchName + '_weight_' + self.branchsuffix)
        self.addVectorBranch(self.branchName + '_weight_' + self.branchsuffix, length=1)


    def processEvent(self, tree):


	self.NbEvents+=1
	print("Event ", self.NbEvents, " out of ", tree.GetEntries())

        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
	    
	    #Compute the polynomial weights
	    weights  = np.array([ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]).reshape(tree.nLHEReweightingWeight,1)

	    coeff     = np.matmul(self.invVandermonde, weights).flatten()
            

	    self.WCtoeval = np.zeros((self.nWC, 1))
	    
	    for i in range(self.nWCtocorrelate):
		#Set the needed WC to 1 to get the total mixed contribution
	    	self.WCtoeval[int(self.WCindices[i]), 0]	= 1
	  

	    self.WCtoeval = self.WCtoeval.T
	    self.WCtoeval = self.polyfeatures.fit_transform(self.WCtoeval)

	    EvaluatedWC = np.matmul(self.WCtoeval, coeff).flatten()

	    self._b(self.branchName + '_weight_' + self.branchsuffix)[0] = EvaluatedWC[0]



if __name__ == '__main__':

	print("WC correlation module, simulate two WC at a time")
	print("Produces the weights for the case of n WC = 1")
