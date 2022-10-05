#!/usr/bin/env python
import ROOT
import sys
ROOT.gROOT.SetBatch(True)
from myutils import BetterConfigParser, ParseInfo
import glob
from multiprocessing import Process


# 1) ./getExtWeights.py configname
# 2) ./getExtWeight.py configname verify

list_broken_file = []

def countEvents(rootFileName, cut = "1"):

    #print 'counting events in', rootFileName
    #print 'cuts is', cut
    f = ROOT.TFile.Open(rootFileName)
    if not f or f.IsZombie():
        print "\x1b[31mnot found:",rootFileName,"\x1b[0m"
        list_broken_file.append(rootFileName)
        nevents = 0
    else: 
        tree = f.Get("Events")
        nevents = 1.* tree.Draw("",cut)
        f.Close()
    #nevents = 1.* tree.GetEntries(cut)
    #print 'number of events is', nevents
    return nevents

def getEventCount(Sample, region="1"):
    sysOut = config.get('Directories','SYSout').strip()
    #sysOut = config.get('Directories','PREPout').strip()
    t3proto = 'root://t3dcachedb.psi.ch:1094'
    sysOutMountedPath = sysOut.replace(t3proto,'').replace('root://t3dcachedb03.psi.ch:1094','')
    #SampleCounts = {}
    count = 0
    #print 'Sample is', Sample
    for sample in Sample:
        cut = ''
        if sample in SampleCuts:
            #print'taking additional', SampleCuts[sample], ' cut for', sample
            cut = region +'&('+SampleCuts[sample]+')'
        else:
            cut = region
        #print 'sample is', sample
        #print 'cut is', cut
        fileMask = "{path}/{sample}/{tree}.root".format(path=sysOutMountedPath, sample=sample, tree='*')
        #print fileMask
        SampleFiles = glob.glob(fileMask)
        #print sample,"=>",len(SampleFiles),"files"
        nevents = 0
        for SampleFile in SampleFiles:
            nevents += countEvents(t3proto + '/' + SampleFile, cut)
            #print 'nevents is', nevents
        #Jprint '-------------'
        #Jprint 'count for file', fileMask, 'is'
        #Jprint nevents
        count += nevents
    #print 'Total count is', count
    return count
    # get total count
    #totalCount = sum([n for sampleName,n in SampleCounts.iteritems()])

    # relative weights
    #for sampleName,n in SampleCounts.iteritems():
    #    print sampleName,":",(1.0*n/totalCount if totalCount > 0 else '-')

def getStichWeight(string, Sample1, Sample2, region):
    nevent_1 = getEventCount(Sample1, region)
    nevent_2 = getEventCount(Sample2, region)
    ##return nevent_1/(nevent_1+nevent_2)
    #print 'sample1 is', Sample1
    #print 'sample2 is', Sample2
    #print 'nevent1 is', nevent_1
    #print 'nevent2 is', nevent_2
    print 'weight for ', string, 'is', nevent_1/(nevent_1+nevent_2)
    print 'list of broken files'
    print list_broken_file

def getExtWeights(config, extParts):

    sysOut = config.get('Directories','SYSout').strip()
    #sysOut = config.get('Directories','PREPout').strip()
    t3proto = 'root://t3dcachedb.psi.ch:1094'
    sysOutMountedPath = sysOut.replace(t3proto,'').replace('root://t3dcachedb03.psi.ch:1094','')
    extPartCounts = {}
    for extPart in extParts:
        fileMask = "{path}/{sample}/{tree}.root".format(path=sysOutMountedPath, sample=extPart, tree='*')
        #print fileMask
        extPartFiles = glob.glob(fileMask)
        #print extPart,"=>",len(extPartFiles),"files"
        extPartCounts[extPart] = 0
        for extPartFile in extPartFiles:
            extPartCounts[extPart] += countEvents(t3proto + '/' + extPartFile)
            #root://t3dcachedb03.psi.ch:1094
    # get total count
    totalCount = sum([n for sampleName,n in extPartCounts.iteritems()])

    # relative weights
    for sampleName,n in extPartCounts.iteritems():
        print sampleName,":",(1.0*n/totalCount if totalCount > 0 else '-')

config = BetterConfigParser()
configPath = sys.argv[1] + '/samples_nosplit.ini'
config.read(configPath)
configPath = sys.argv[1] + '/paths.ini'
config.read(configPath)

sampleDict = {}
sampleWeights = {}
verify = len(sys.argv) > 2 and sys.argv[2]=='verify'

##############
# Cuts
##############

preselectionCut =  "(( (Vtype == 2 && <!General|muTrig!>) || (Vtype == 3 && <!General|eTrig!> )) && <!General|BasicCuts!>)"

SampleCuts = {"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":"(LHE_HT<100)","DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":"(LHE_HT<100)"}

####
#Phase-space cuts
####

BJets           = "(LHE_Nb>0)"
BGenFilerJets   = "(LHE_Nb==0) && Sum$(GenPart_status == 2 && (abs((GenPart_pdgId)/100)%10 == 5 ) || (abs(GenPart_pdgId/1000)%10 == 5)) > 0"

VPT0              = "(LHE_Vpt<100)"
VPT100            = "(LHE_Vpt>100 && LHE_Vpt<200)"
VPT200            = "(LHE_Vpt>200)"

##############
# Samples
##############




WlvjetsHTbinned  = ['WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8']


#Benriched
WlvBjets         = ["WBJetsToLNu_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","WBJetsToLNu_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]

#BGenFilter
WlvBGenFilterjets         = ["WJetsToLNu_BGenFilter_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","WJetsToLNu_BGenFilter_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]

########
#DY+Jets
########

DYjetsHTbinned  = ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8']

DYBjets     = ['DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8']

DYBjetsIncl = ['DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8']
DYBjetsVpt = ['DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8']

DYBGenFilterjets         = ["DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]

def runInParallel(func, arglist):
    proc = []
    for arg_ in arglist:
        #p = Process(target=getStichWeight, args=arg_)
        p = Process(target=func, args=arg_)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

########
## W+jets
########

#
## for b-enriched
#print 'B-enriched'
#arglist = [
#    ('VPT100', WlvjetsHTbinned, WlvBjets,VPT100 +"&&"+BJets,),
#    ('VPT200', WlvjetsHTbinned, WlvBjets,VPT200 +"&&"+BJets,),
#        ]
#runInParallel(getStichWeight,arglist)
#
#print 'BGenFilter'
## for b-enriched
#arglist = [
#    ('VPT100',WlvjetsHTbinned, WlvBGenFilterjets,VPT100 +"&&"+BGenFilerJets,),
#    ('VPT200',WlvjetsHTbinned, WlvBGenFilterjets,VPT200 +"&&"+BGenFilerJets,),
#        ]
#runInParallel(getStichWeight,arglist)

#######
# Z+jets
#######

# for b-enriched
print 'B-enriched'
arglist = [
    ('VPT0', DYjetsHTbinned, DYBjetsIncl, VPT0 +"&&"+BJets,),
    ('VPT100', DYjetsHTbinned, DYBjetsVpt, VPT100 +"&&"+BJets,),
    ('VPT200', DYjetsHTbinned, DYBjetsVpt, VPT200 +"&&"+BJets,),
        ]
runInParallel(getStichWeight,arglist)

## inclusive vs vpt
#print 'B-enriched: inclusive vs Vpt'
#arglist = [
#    ('VPT0', DYBjetsIncl, DYBjetsVpt, VPT0 +"&&"+BJets,),
#    ('VPT100', DYBjetsIncl, DYBjetsVpt, VPT100 +"&&"+BJets,),
#    ('VPT200', DYBjetsIncl, DYBjetsVpt, VPT200 +"&&"+BJets,),
#        ]
#runInParallel(getStichWeight,arglist)

#print 'BGenFilter'
## for b-enriched
#arglist = [
#    ('1', DYjetsHTbinned, DYBGenFilterjets, BGenFilerJets,),
#        ]
#runInParallel(getStichWeight,arglist)

print 'all broken files', list_broken_file
