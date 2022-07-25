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

class EFT_scaling(AddCollectionsModule):


    def __init__(self, branchName='EFT_scaling'):
        super(EFT_scaling, self).__init__()
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
                
        self.ScalingParam = self.config.get('Scaling', 'Range').split(",")
        

        self.RangeStep = int(self.config.get('Scaling', 'NbofPoints'))
        self.ScalingRange = np.linspace(int(self.ScalingParam[0]), int(self.ScalingParam[1]), self.RangeStep) 
        self.NbEvents = 0
        self.sampleTree = initVars['sampleTree']
        self.WCtoeval = 0

        for i in range(self.nWC):
                self.addVectorBranch(self.branchName + '_weight_c' + str(i), length=self.RangeStep)
                self.addVectorBranch(self.branchName + '_weight_c' + str(i) +"_index", length=self.RangeStep)


    def processEvent(self, tree):


	self.NbEvents+=1
	print("Event ", self.NbEvents, " out of ", tree.GetEntries())

        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
	    
	    #Compute the polynomial weights
	    weights  = np.array([ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]).reshape(tree.nLHEReweightingWeight,1)

	    coeff     = np.matmul(self.invVandermonde, weights).flatten()
            
	    for i in range(self.nWC):
            	    #self._b(self.branchName + '_weight_c' + str(i))[0] = -99997.0

		    self.WCtoeval = np.zeros((self.nWC, self.RangeStep))
	            self.WCtoeval[i, :]	= self.ScalingRange
	            self.WCtoeval	= self.WCtoeval.T
		    self.WCtoeval = self.polyfeatures.fit_transform(self.WCtoeval)

		    EvaluatedWC = np.matmul(self.WCtoeval, coeff).flatten()

		    for n in range(self.RangeStep):
            	    	self._b(self.branchName + '_weight_c' + str(i))[n] = EvaluatedWC[n]
            	    	self._b(self.branchName + '_weight_c' + str(i) +"_index")[n] = n+1


    def interpret_weight(self, weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res


if __name__ == '__main__':

	print("Scaling module")

