#!/usr/bin/env python
from __future__ import print_function
import ROOT
ROOT.gROOT.SetBatch(True)
import NewTreeCache as TreeCache
from sampleTree import SampleTree as SampleTree
from Datacard import Datacard
from sample_parser import ParseInfo
import resource
import sys
import pickle
import glob
import shutil
import numpy as np
import pandas as pd
import os

class SampleTreesToDataFrameConverter(object):

    #TODO: take regions as argument
    def __init__(self, config,config_name="DC"):
        self.config = config
        self.config_name = config_name
        if config.has_option('LimitGeneral', 'List_for_rebinner'):
            self.regions = [x.strip() for x in config.get('LimitGeneral', 'List_for_rebinner').split(',') if len(x.strip()) > 0]
        else:
            self.regions = [x.strip() for x in config.get('LimitGeneral', 'List').split(',') if len(x.strip()) > 0 and config.get('dc:%s'%(x.strip()),'type').lower()!='cr' ]
        
        print("Convert regions:")
        print(", ".join(self.regions))
        self.dcMakers = {region: Datacard(config=self.config, region=region) for region in self.regions}
        self.dfs = {region: {} for region in self.regions} 
        

    def loadSamples(self,safe_hdf=False,path="dumps/",force=False):
        for region in self.regions:
            print("\n\n==========================================\nINFO: Load samples for %s"%region)
            dcMaker = self.dcMakers[region]
            for sampleType in ['SIG','BKG']:
                treevar = dcMaker.treevar
                hdf_filename = "%s_%s_%s.hdf"%(self.config_name,region,treevar.split(".")[0])
                if not force:
                    try:
                        df = pd.read_hdf("%s/%s"%(path,hdf_filename),sampleType)
                        print("INFO: HDF file for %s %s %s %s existing. Skip loading"%(self.config_name,region,treevar,sampleType))
                        self.dfs[region][sampleType] = df
                        continue
                    except:
                        print("INFO: No HDF file found matching %s %s %s %s. Load samples from cache."%(self.config_name,region,treevar,sampleType))
                        pass
                df = pd.DataFrame({treevar: [], 'weight': []})         
                samples = dcMaker.samples[sampleType]
                weightF = dcMaker.weightF
                EvalCut = dcMaker.EvalCut
                try:
                    for i, sample in enumerate(samples): 
                        print("INFO: Add ",sampleType," sample ", i, " of ", len(samples))
                        # get sample tree from cache
                        tc = TreeCache.TreeCache(
                                sample=sample,
                                cutList=dcMaker.getCacheCut(sample),
                                inputFolder=dcMaker.path,
                                config=self.config,
                                debug=False
                            )
                        if not tc.isCached():
                            print("\x1b[31m:ERROR not cached! run cachedc step again\x1b[0m")
                            raise "Tree not cached"
                        sampleTree = tc.getTree()
               
                        if sampleTree:
                            treeScale = sampleTree.getScale(sample) * 2.0
                            print ('scale:', treeScale)
                            
                            # initialize numpy array
                            nSamples = sampleTree.GetEntries()
                            input_treevar = []
                            input_weight = []
                            sampleTree.addFormula(treevar)
                            sampleTree.addFormula(weightF)
                            sampleTree.addFormula(EvalCut)
                            
                            # fill numpy array from ROOT tree
                            for i, event in enumerate(sampleTree):
                                if sampleTree.evaluate(EvalCut):
                                    input_treevar.append(sampleTree.evaluate(treevar))
                                    input_weight.append(sampleTree.evaluate(weightF)*treeScale)
                                    # total weight comes from weightF (btag, lepton sf, ...) and treeScale to scale MC to x-section

                            sampleDf = pd.DataFrame({treevar: input_treevar, 'weight': input_weight})
                            df = df.append(sampleDf,ignore_index=True)
                            
                        else:
                            print ("\x1b[31mERROR: TREE NOT FOUND:", sample.name, " -> not cached??\x1b[0m")
                            raise "Tree not cached"
                        
                    if safe_hdf:
                        print("INFO: safe %s in %s/%s"%(sampleType,path,hdf_filename))
                        if not os.path.exists(path):
                                os.makedirs(path)
                        df.to_hdf("%s/%s"%(path,hdf_filename),sampleType)
                    
                    self.dfs[region][sampleType] = df
                except:
                    print("\x1b[31mERROR: reading %s samples for %s failed"%(sampleType,region))

    def getDataFrames(self):
        return self.dfs

    def makeHDF(self,path="./"):
        for region, dataframes in self.dfs.iteritems():
            for dataType, df in dataframes.iteritems():
                df.to_hdf(path+"/"+region+".hdf",dataType)



