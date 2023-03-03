#!/usr/bin/env python
from __future__ import print_function
import ROOT
import numpy as np
from Jet import Jet
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import sys
import os
import json
from myutils.XbbTools import XbbTools
from myutils.XbbTools import XbbMvaInputsList

# tfZllDNN repository has to be cloned inside the python folder
sys.path.append("..")
sys.path.append("../LGBMvhbb/")
from LGBMvhbb.LGBMregressor import LGBMRegressor

# reads tensorflow checkpoints and applies DNN output to ntuples,
# including variations from systematic uncertainties
# needs: tensorflow >=1.4
class LGBMEvaluator(AddCollectionsModule):

    def __init__(self, mvaName, nano=False, condition=None):
        self.mvaName = mvaName
        self.nano = nano
        self.debug = False
        self.condition = condition
        self.fixInputs = []
        super(LGBMEvaluator, self).__init__()

    def customInit(self, initVars):
        self.config = initVars['config']
        self.sampleTree = initVars['sampleTree']
        self.sample = initVars['sample']

        if self.condition:
            self.sampleTree.addFormula(self.condition)

        self.hJidx = self.config.get('General', 'hJidx') if self.config.has_option('General', 'hJidx') else 'hJidx'


        #PATH and NAMES
        self.checkpointpath = self.config.get(self.mvaName, 'checkpointpath') if self.config.has_option(self.mvaName, 'checkpointpath') else None

        print("PATH", self.checkpointpath)
        
        self.checkpointname = self.config.get(self.mvaName, 'checkpointprefix') if self.config.has_option(self.mvaName, 'checkpointprefix') else None

        print("Name", self.checkpointname)

        #if self.config.has_option(self.mvaName, 'branchName'):
        #    self.branchName = self.config.get(self.mvaName, 'branchName')
        #else:
        #    print("\x1b[31mERROR: 'branchName' option missing for MVA config section [%s]! .../model.ckpt has to be specified to be able to restore classifier.\x1b[0m"%self.mvaName)
        #    raise Exception("BranchNameError")

        if self.checkpointpath is None or self.checkpointname is None:
            print("\x1b[31mERROR: 'checkpointname' or 'checkpointpath' option missing for MVA config section [%s]! .../model.ckpt has to be specified to be able to restore classifier.\x1b[0m"%self.mvaName)
            raise Exception("CheckpointError")

        #Components of WC to add
        #String of WC indices to consider, separated by commas in training.ini

        self.WC = self.config.get(self.mvaName, 'WC') if self.config.has_option(self.mvaName, 'WC') else None
        if self.WC is None:
            print("\x1b[31mERROR: 'WC' option missing for MVA config section [%s]! .../model.ckpt has to be specified to be able to restore classifier.\x1b[0m"%self.mvaName)
            raise Exception("WCError")

        #Create the list of components:
        self.quadOnly = self.config.get(self.mvaName, 'quadOnly') if self.config.has_option(self.mvaName, 'quadOnly') else None

        self.regressors = []
        self.checkpoints = []
        self.WCpostfix = self.WC.split(",")
        
        #The quadOnly parameter is for testing purposes. If it is true, only the Quad parts of the WC will be considered (no mixing, no Lin)
        #If it is false, only the Lin part wil be considered (no mixing, no Quad)
        #If it is not present in the config, all the components will be evaluated (all Lin and Quad and mixed terms for the WC list). 

        for i in range(len(self.WCpostfix)):
            #Quad
            if self.quadOnly:
                fullname, _= LGBMRegressor.define_names(WCs = self.WCpostfix[i], fullname = self.checkpointname, quad = True)
                regressor = LGBMRegressor.load(fullname = fullname, path = self.checkpointpath)
                self.regressors.append(regressor) 
                self.checkpoints.append(fullname) 
            #Lin
            elif self.quadOnly is not None:
                fullname, _= LGBMRegressor.define_names(WCs = self.WCpostfix[i], fullname = self.checkpointname, quad = False)
                regressor = LGBMRegressor.load(fullname = fullname, path = self.checkpointpath)
                self.regressors.append(regressor) 
                self.checkpoints.append(fullname) 
            #Lin, Mixed and Quad
            else: 
                fullname, _= LGBMRegressor.define_names(WCs = self.WCpostfix[i], fullname = self.checkpointname, quad = False)
                regressor = LGBMRegressor.load(fullname = fullname, path = self.checkpointpath)
                self.regressors.append(regressor) 
                self.checkpoints.append(fullname) 
                
                fullname, _= LGBMRegressor.define_names(WCs = self.WCpostfix[i], fullname = self.checkpointname, quad = True)
                regressor = LGBMRegressor.load(fullname = fullname, path = self.checkpointpath)
                self.regressors.append(regressor) 
                self.checkpoints.append(fullname) 
                
                for j in range(i+1, len(self.WCpostfix)):
                    fullname, _= LGBMRegressor.define_names(WCs = self.WCpostfix[i] + "," + self.WCpostfix[j], fullname = self.checkpointname, quad = True)
                    regressor = LGBMRegressor.load(fullname = fullname, path = self.checkpointpath)
                    self.regressors.append(regressor) 
                    self.checkpoints.append(fullname) 
                           

        print("POSTFIX", self.checkpoints)



        # FEATURES
        self.featuresConfig = None
        self.featuresCheckpoint = None
        self.featuresConfigElectrons = None
        self.featuresConfigMuons = None
        

        try:
            self.featuresConfig = self.config.get(self.config.get(self.mvaName, "treeVarSet"), "Nominal").strip().split(" ")
            self.featuresConfigElectrons = self.config.get(self.config.get(self.mvaName, "treeVarSet"), "OnlyElectrons").strip().split(" ")
            self.featuresConfigMuons = self.config.get(self.config.get(self.mvaName, "treeVarSet"), "OnlyMuons").strip().split(" ")
        except Exception as e:
            print("WARNING: could not get treeVarSet from config:", e)
       
        print("FEATURES", self.featuresConfig)
        print("Only Electrons", self.featuresConfigElectrons)
        print("Only Muons", self.featuresConfigMuons)


        


 

        #Add features from checkpoint

        self.features = self.featuresConfig 
        self.featuresElectrons = self.featuresConfigElectrons 
        self.featuresMuons = self.featuresConfigMuons

        self.featureList = XbbMvaInputsList(self.features, config=self.config)
        self.featureListElectrons = XbbMvaInputsList(self.featuresElectrons, config=self.config)
        self.featureListMuons = XbbMvaInputsList(self.featuresMuons, config=self.config)
        self.nFeatures   = self.featureList.length()
        self.nFeaturesElectrons   = self.featureListElectrons.length()
        self.nFeaturesMuons   = self.featureListMuons.length()


        assert self.nFeaturesElectrons == self.nFeaturesMuons, "Leptons features number do not match"


        # SIGNAL definition
        self.signalIndex = 0
        if self.config.has_option(self.mvaName, 'signalIndex'):
            self.signalIndex = eval(self.config.get(self.mvaName, 'signalIndex'))
        self.signalClassIds = [self.signalIndex]

        # systematics
        self.systematics = self.config.get(self.mvaName, 'systematics').split(' ') if self.config.has_option(self.mvaName, 'systematics') else self.config.get('systematics', 'systematics').split(' ')

        # create output branches
        self.bitCollections = []
        for name in self.checkpoints:
            self.bitCollection = Collection(name, self.systematics, leaves=True)
            self.addCollection(self.bitCollection)
            self.bitCollections.append(self.bitCollection)



        # create formulas for input variables
        self.inputVariables = {}
        self.inputVariablesElectrons = {}
        self.inputVariablesMuons = {}
        for syst in self.systematics:
            systBase, UD              = XbbTools.splitSystVariation(syst, sample=self.sample)
            self.inputVariables[syst] = [XbbTools.sanitizeExpression(self.featureList.get(i, syst=systBase, UD=UD), self.config, debug=self.debug) for i in range(self.nFeatures)]
            self.inputVariablesElectrons[syst] = [XbbTools.sanitizeExpression(self.featureListElectrons.get(i, syst=systBase, UD=UD), self.config, debug=self.debug) for i in range(self.nFeaturesElectrons)]
            self.inputVariablesMuons[syst] = [XbbTools.sanitizeExpression(self.featureListMuons.get(i, syst=systBase, UD=UD), self.config, debug=self.debug) for i in range(self.nFeaturesMuons)]
            for var in self.inputVariables[syst]:
                self.sampleTree.addFormula(var)
            for var in self.inputVariablesElectrons[syst]:
                self.sampleTree.addFormula(var)
            for var in self.inputVariablesMuons[syst]:
                self.sampleTree.addFormula(var)


        self.Vtype_index = self.featuresConfig.index("Vtype")

        print("Vtype index is", self.Vtype_index)
        print("Number of features : ", self.nFeatures + self.nFeaturesElectrons) 


        #self.inputVariables = []       
        #self.inputVariablesElectrons = []       
        #self.inputVariablesMuons = []       

 
        #for name in self.checkpoints:
        #    self.addBranch(name)

        # create formulas for input variables
        #print([self.featureList.get(i) for i in range(self.nFeatures)])
        
        #for i in range(self.nFeatures):
        #    print("feature", i, self.featureList.get(i))
        #    self.sampleTree.addFormula(self.featureList.get(i))
        #    self.inputVariables.append(self.featureList.get(i))

        ##Split in electrons and muons 

        #for i in range(self.nFeaturesElectrons):
        #    print("feature electron", i, self.featureListElectrons.get(i))
        #    self.sampleTree.addFormula(self.featureListElectrons.get(i))
        #    self.inputVariablesElectrons.append(self.featureListElectrons.get(i))


        #for i in range(self.nFeaturesMuons):
        #    print("feature muon", i, self.featureListMuons.get(i))
        #    self.sampleTree.addFormula(self.featureListMuons.get(i))
        #    self.inputVariablesMuons.append(self.featureListMuons.get(i))

        # create tensorflow graph
        #self.ev = TensorflowDNNEvaluator(checkpoint=self.checkpoint, scaler=self.scalerDump)
        
        self.NbEvents = 0
      



    # called from main loop for every event
    def processEvent(self, tree):
        
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            self.NbEvents+=1
            print("Event ", self.NbEvents, " out of ", tree.GetEntries())
            
            if self.condition is None or self.sampleTree.evaluate(self.condition):
                # fill input variables
                inputs = np.full((len(self.systematics), self.nFeatures + self.nFeaturesElectrons), 0.0, dtype=np.float32)
                for j, syst in enumerate(self.systematics):
                    for i, var in enumerate(self.inputVariables[syst]):
                        inputs[j,i] = self.sampleTree.evaluate(var)
                    
                    #Choose electron or muon feature
                    for i, var in enumerate(self.inputVariablesElectrons[syst]):
                        Vtype = inputs[0,self.Vtype_index]

                        #Muons
                        if not Vtype%2:
                            thisvar = self.inputVariablesMuons[syst][i]                    
                            inputs[j,i + self.nFeatures] = self.sampleTree.evaluate(thisvar)
                        #Electrons
                        else:
                            thisvar = self.inputVariablesElectrons[syst][i]                    
                            inputs[j,i + self.nFeatures] = self.sampleTree.evaluate(thisvar)

                for i,bitCollection in enumerate(self.bitCollections):
                    for j, syst in enumerate(self.systematics):
                        bitCollection[bitCollection.name][j] = self.regressors[i].gbm.predict(inputs[j:j+1,:], num_threads = 10) 



                 
            else:
                print("The condition is not met, check your config")
                raise Exception("CheckpointError")

 
        return True


