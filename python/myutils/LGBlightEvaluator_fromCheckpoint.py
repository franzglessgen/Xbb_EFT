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

        if self.config.has_option(self.mvaName, 'branchName'):
            self.branchName = self.config.get(self.mvaName, 'branchName')
        else:
            print("\x1b[31mERROR: 'branchName' option missing for MVA config section [%s]! .../model.ckpt has to be specified to be able to restore classifier.\x1b[0m"%self.mvaName)
            raise Exception("BranchNameError")

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
        try:
            self.featuresConfig = self.config.get(self.config.get(self.mvaName, "treeVarSet"), "Nominal").strip().split(" ")
        except Exception as e:
            print("WARNING: could not get treeVarSet from config:", e)
       
        print("FEATURES", self.featuresConfig)

 
        #TOADD: save the list of variable names in checkpoint
        #if 'variables' in self.info:
        #    self.featuresCheckpoint = self.info['variables']

        #Add features from checkpoint

        self.features = self.featuresConfig 

        #self.features = [self.features[i].replace("VHbb::", "ROOT.VHbb.") for i in range(len(self.features))]  
        self.featureList = XbbMvaInputsList(self.features, config=self.config)
        self.nFeatures   = self.featureList.length()


        # SIGNAL definition
        self.signalIndex = 0
        if self.config.has_option(self.mvaName, 'signalIndex'):
            self.signalIndex = eval(self.config.get(self.mvaName, 'signalIndex'))
        self.signalClassIds = [self.signalIndex]

        # systematics
        #self.systematics = self.config.get(self.mvaName, 'systematics').split(' ') if self.config.has_option(self.mvaName, 'systematics') else self.config.get('systematics', 'systematics').split(' ')

        # create output branches
        self.inputVariables = []       

 
        for name in self.checkpoints:
            self.addBranch(name)

        # create formulas for input variables
        #print([self.featureList.get(i) for i in range(self.nFeatures)])
        
        for i in range(self.nFeatures):
            print("feature", i, self.featureList.get(i))
            self.sampleTree.addFormula(self.featureList.get(i))
            self.inputVariables.append(self.featureList.get(i))

        #Reorder inputVariables to match training
         




        # create tensorflow graph
        #self.ev = TensorflowDNNEvaluator(checkpoint=self.checkpoint, scaler=self.scalerDump)
        
        self.NbEvents = 0
      



    # called from main loop for every event
    def processEvent(self, tree):

        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
           
            self.NbEvents+=1
            print("Event ", self.NbEvents, " out of ", tree.GetEntries())
            #fake_features = np.ones((1,44))
            #pred = Regressor.gbm.predict(fake_features)
            #print(pred)
    

            if self.condition is None or self.sampleTree.evaluate(self.condition):
                # fill input variables
                inputs = np.full((1, self.nFeatures), 0.0, dtype=np.float32)
                for i, var in enumerate(self.inputVariables):
                        if var in self.fixInputs:
                            inputs[0,i] = self.fixInputs[var]
                        else:
                            inputs[0,i] = self.sampleTree.evaluate(var)

                # use TensorflowDNNEvaluator
                
                #Softmax probas, replace with BDT evaluation
                
                for i in range(len(self.checkpoints)):
                     
                    self._b(self.checkpoints[i])[0] = self.regressors[i].gbm.predict(inputs) 

                 
            else:
                print("The condition is not met, check your config")
                raise Exception("CheckpointError")

 
        return True


