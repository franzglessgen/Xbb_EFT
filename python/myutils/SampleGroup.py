#!/usr/bin/env python
from __future__ import print_function
from myutils import ParseInfo
from BranchTools import AddCollectionsModule

class SampleGroup(AddCollectionsModule):

    def __init__(self, groupDict=None, prefix="is_", eventCounts=None):
        super(SampleGroup, self).__init__()
        self.groupDict = groupDict
        self.prefix = prefix
        if eventCounts:
            with open(eventCounts, 'r') as inFile:
                self.eventCountsDict = eval(inFile.read())
        else:
            self.eventCountsDict = None

        self.branches = []
        self.eventNumberOffset = 0L

    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.sampleTree = initVars['sampleTree']
        self.config = initVars['config']
        self.samplesInfo = ParseInfo(samples_path=self.config.get('Directories', 'dcSamples'), config=self.config)
	self.subsamples = [x for x in self.samplesInfo if x.identifier == self.sample.identifier and x.subsample]
        print("INFO: subsamples/cut")
        for s in self.subsamples:
            print(" >", s.name, s.subcut)
            self.sampleTree.addFormula(s.subcut)


        if not self.groupDict:
            self.groupDict = eval(self.config.get('LimitGeneral','Group'))

        self.groupNames = list(set(self.groupDict.values()))
        self.groups = {k: [x for x,y in self.groupDict.items() if y==k] for k in self.groupNames}

        for groupName, sampleNames in self.groups.items():
            self.branches.append({'name': self.prefix + groupName, 'formula': self.isInGroup, 'arguments': groupName}) 


        self.branches.append({'name': 'sampleIndex', 'formula': self.getSampleIndex, 'type': 'i'})

        if self.eventCountsDict:
            self.branches.append({'name': 'event_unique', 'formula': self.getEventNumber, 'type': 'l'})

            if len(self.sampleTree.sampleFileNames) != 1:
                print("ERROR: adding unique event numbers for chains is not implemented!")
                raise Exception("SampleGroup__customInit__not_implemented")
            self.eventNumberOffset = self.eventCountsDict[self.sample.identifier][self.sampleTree.sampleFileNames[0]]

    def getBranches(self):
        return self.branches

    def isInGroup(self, event, arguments):
        groupName = arguments
        if len(self.subsamples) > 0:
            for s in self.subsamples:
                if self.sampleTree.evaluate(s.subcut):
                    return s.name in self.groups[groupName]
        else:
            return self.sample.name in self.groups[groupName]

    def getSampleIndex(self, event):
        if len(self.subsamples) > 0:
            for s in self.subsamples:
                if self.sampleTree.evaluate(s.subcut):
                    return s.index 
            #print("\x1b[31mERROR: none of the subcuts applies for event in ", self.sample.identifier, "! => no sample index can be defined!\x1b[0m")
            #raise Exception("SampleSubcutDefinitionError")
            return self.sample.index
        else:
            return self.sample.index

    def getEventNumber(self, event):
        return self.eventNumberOffset + event.GetReadEntry()
