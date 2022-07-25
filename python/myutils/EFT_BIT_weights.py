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

class EFT_BIT_weights(AddCollectionsModule):


    def __init__(self, branchName='EFT_BIT_weights'):
        super(EFT_BIT_weights, self).__init__()
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
                

        #self.RangeStep = int(self.config.get('Scaling', 'NbofPoints'))
            #self.ScalingRange = np.linspace(int(self.ScalingParam[0]), int(self.ScalingParam[1]), self.RangeStep) 
        self.NbEvents = 0
        self.sampleTree = initVars['sampleTree']

        

        self.addVectorBranch(self.branchName , length=self.ndof)


    def processEvent(self, tree):


	self.NbEvents+=1
	print("Event ", self.NbEvents, " out of ", tree.GetEntries())

        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
	    
	    #Compute the polynomial weights
	    weights  = np.array([ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]).reshape(tree.nLHEReweightingWeight,1)

	    coeff     = np.matmul(self.invVandermonde, weights).flatten()
            
        for s in range(len(coeff)):

	        self._b(self.branchName)[s] = coeff[s]



if __name__ == '__main__':

	print("WC correlation module, simulate two WC at a time")
	print("Produces the weights for the case of n WC = 1")
