#!/usr/bin/env python
import ROOT
from WeightInfo                   import WeightInfo
from HyperPoly                    import HyperPoly

#from ROOT import TVector3, TLorentzVector
#from math import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
import numpy
#import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class EFT_reweight(AddCollectionsModule):


    def __init__(self, branchName='EFT_reweight', wc_choice = 0,  wc_start = -1.0, wc_end = 2.0, wc_step = 0.1):
        super(EFT_reweight, self).__init__()
        self.branchName = branchName
        self.wc_choice = wc_choice
        self.wc_start = wc_start
        self.wc_end = wc_end
        self.wc_step = wc_step


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']
       
     
        self.addVectorBranch(self.branchName + '_weight', length=255)
        #self.addBranch(self.branchName + '_np')
        #self.addBranch(self.branchName + '_chi2_ndof')
        #self.addBranch(self.branchName + '_nWC')

        reweight_pkl = "/work/vperovic/CMSSW_10_1_0/src/Xbb/python/myutils/SMEFTsim_VH_reweight_card.pkl"        # this has the information on WCs written in the reweight card (particularly, which WCs are used and in what combination)
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
            #self._b(self.branchName + '_np')[0] = -99997.0
            #self._b(self.branchName + '_chi2_ndof')[0] = -99997.0
            #self._b(self.branchName + '_nWC')[0] = -99997.0

            #fl = ROOT.TFile.Open("NanoAOD_ZH.root")
            #tree = fl.Get("Events")
            #tree.GetEntry(0)

            #weights  = [ tree.LHEReweightingWeight[i] for i in range(tree.nLHEReweightingWeight) ]
            #coeff     = self.hyperPoly.get_parametrization( weights )
            #np        = self.hyperPoly.ndof                                 # no. of coefficients in polynomial
            #chi2_ndof = self.hyperPoly.chi2_ndof( coeff, weights )          # chi^2/ndof of the polynomial fit
            #nWC       = self.hyperPoly.nvar                                 # no. of Wilson coefficients

            #coeff = [tree.EFT_params_coeff[i] for i in range(int(tree.EFT_params_np))]
            coeff = [tree.EFT_params_coeff[i] for i in range(tree.EFT_params_np)]
            np = tree.EFT_params_np
            chi2_ndof = tree.EFT_params_chi2_ndof
            nWC = tree.EFT_params_nWC


            #for i in range(np): self._b(self.branchName + '_coeff')[i] = coeff[i]
            #self._b(self.branchName + '_np')[0] = np
            #self._b(self.branchName + '_chi2_ndof')[0] = chi2_ndof
            #self._b(self.branchName + '_nWC')[0] = nWC

            
            #WC_start = -1.0
            #WC_end = 2.0
            #WC_step = 0.1
            #target_WC=1 #according to the reweight card 0 is cHj1, 1 is cHj3 etc
            WCs = []
            #for ic in range(int(nWC)):
            for ic in range(nWC):
                WCs.append(0)
            WC_vals = numpy.arange(self.wc_start, self.wc_end + self.wc_step, self.wc_step)
            for cnt, WC_val in enumerate(WC_vals):
                WCs[self.wc_choice] = WC_val
                self._b(self.branchName + '_weight')[cnt] = self.hyperPoly.eval(coeff, *WCs)
                
            

    def interpret_weight(self, weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res




#Examples to get weight for few combinations of WCs are below

#for ieft in range(nWC):
#    WCs = []
#    for ic in range(hyperPoly.nvar):
#        WCs.append(0)
#    if(ieft>0):
#        WCs[ieft-1] = 1.
#    EFT_wt =  hyperPoly.eval(coeff,*WCs)
#    print("Weight id:",ieft+1,"weight:",EFT_wt)

