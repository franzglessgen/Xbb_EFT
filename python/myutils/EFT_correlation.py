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
        self.sampleTree = initVars['sampleTree']
        
        self.nWC = int(self.config.get('WCGeneral', 'nbofWC'))
        self.ndof = 2*self.nWC + self.nWC*(self.nWC -1)/2 + 1        
        
        #This hold the inverse Vandermonde matrix for computation of the interpolation
        self.VandermondePath = self.config.get('WCGeneral', 'SimWC')
        self.invVandermonde = np.load(self.VandermondePath)
        self.polyfeatures = PolynomialFeatures(2)		
        self.NbEvents = 0
        self.WCtoeval = np.zeros((self.nWC, self.ndof - 1))
        for i in range(self.nWC):
            #SM + LIN + QUAD needs to evaluate the poly at c = 1
            self.WCtoeval[i, i]	= 1
	        #QUAD needs to evaluate the poly at c = -1 also
            self.WCtoeval[i, i + self.nWC]	= -1
            for j in range(i + 1, self.nWC):
		        #Set the needed WCs to 1 to get the total mixed contribution
                self.WCtoeval[i, j - i -1 + 2*self.nWC + (self.nWC - 1)*i - i*(i -1)/2]	= 1
                self.WCtoeval[j, j - i -1 + 2*self.nWC + (self.nWC - 1)*i - i*(i -1)/2]	= 1
	             
    
        self.WCtoeval = self.WCtoeval.T
        #print("SM+LIN+QUAD", self.WCtoeval[:,:self.nWC])
        #print("QUAD", self.WCtoeval[self.nWC:2*self.nWC, :])
        #print("MIXED", self.WCtoeval[2*self.nWC:2*self.nWC + 18, :])
        self.WCtoeval = self.polyfeatures.fit_transform(self.WCtoeval)
        self.addVectorBranch(self.branchName + '_weight', length=self.ndof)


    def processEvent(self, tree):

        self.NbEvents+=1
        print("Event ", self.NbEvents, " out of ", tree.GetEntries())

        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
	    
	    #Compute the polynomial weights
	    weights  = np.array([ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]).reshape(tree.nLHEReweightingWeight,1)

	    coeff     = np.matmul(self.invVandermonde, weights).flatten()
            
	    EvaluatedWC = np.matmul(self.WCtoeval, coeff).flatten()
        
        self._b(self.branchName + '_weight')[0] = 1
        #SM + LIN + QUAD
        for i in range(1, self.nWC+1):
            #print("LIN", i-1)
            self._b(self.branchName + '_weight')[i] = EvaluatedWC[i-1]

        #QUAD
        for i in range(self.nWC+1, 2*self.nWC + 1):
            #print("QUAD", i-1-self.nWC, i-1)
            quad_weight = 0.5*(EvaluatedWC[i - 1 - self.nWC] + EvaluatedWC[i - 1] - 2)
            self._b(self.branchName + '_weight')[i] = quad_weight
        #MIXED
        for i in range(2*self.nWC + 1, self.ndof):
            #print("MIXED", i-1)
            self._b(self.branchName + '_weight')[i] = EvaluatedWC[i-1]



if __name__ == '__main__':

	print("WC correlation module, simulate two WC at a time")
	print("Produces the weights for the case of n WC = 1")
