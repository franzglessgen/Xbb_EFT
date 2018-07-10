#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import os

class METtriggerFromWS(AddCollectionsModule):

    def __init__(self, workspace="data/met/vhbb_metsf.root", fName="data_nominal_110OR170"):
        self.debug = 'XBBDEBUG' in os.environ
        self.workspaceFileName = workspace
        self.fName = fName
        super(METtriggerFromWS, self).__init__()

    def loadWorkspace(self, wsName):
        rootFile = ROOT.TFile.Open(wsName, "READ")
        wsp = rootFile.Get("w")
        wsp.Print()
        rootFile.Close()
        self.met_trigger_sf = wsp.function(self.fName).functor(ROOT.RooArgList(wsp.argSet("met_mht")))

    def customInit(self, initVars):
        self.isData = initVars['sample'].isData()
        if not self.isData:
            self.addBranch('weight_mettrigSF', 1.0)
            self.loadWorkspace(self.workspaceFileName)
    
    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree) and not self.isData:
            self.markProcessed(tree)

            met_mht = array.array('d', [min(min(tree.MET_pt, tree.MHT_pt), 500.0)])
            self.setBranch('weight_mettrigSF', self.met_trigger_sf.eval(met_mht[0]) if met_mht[0] >= 120 else 0.0)

        return True