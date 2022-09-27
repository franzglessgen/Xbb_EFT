#!/usr/bin/env python
import ROOT
import sys
ROOT.gROOT.SetBatch(True)
from myutils import BetterConfigParser, ParseInfo
import glob
from multiprocessing import Process


# 1) ./getExtWeights.py configname
# 2) ./getExtWeight.py configname verify

def countEvents(rootFileName, cut = "1"):

    f = ROOT.TFile.Open(rootFileName)
    if not f or f.IsZombie():
        print "\x1b[31mnot found:",rootFileName,"\x1b[0m"
    tree = f.Get("tree")

    nevents = 1.* tree.Draw("",cut)
    f.Close()
    return nevents

def getEventCount(Sample, region="1"):
    sysOut = config.get('Directories','SYSout').strip()
    t3proto = 'root://t3dcachedb.psi.ch:1094'
    sysOutMountedPath = sysOut.replace(t3proto,'').replace('root://t3dcachedb03.psi.ch:1094','')
    #SampleCounts = {}
    count = 0
    for sample in Sample:
        cut = ''
        if sample in SampleCuts:
            #print'taking additional', SampleCuts[sample], ' cut for', sample
            cut = region +'&('+SampleCuts[sample]+')'
        else:
            cut = region
        fileMask = "{path}/{sample}/{tree}.root".format(path=sysOutMountedPath, sample=sample, tree='*')
        #print fileMask
        SampleFiles = glob.glob(fileMask)
        print sample,"=>",len(SampleFiles),"files"
        nevents = 0
        for SampleFile in SampleFiles:
            nevents += countEvents(t3proto + '/' + SampleFile, cut)
        #print 'count for file', fileMask, 'is', nevents
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
    #return nevent_1/(nevent_1+nevent_2)
    print 'weight for ', string, 'is', nevent_1/(nevent_1+nevent_2)

#def getExtWeights2(Sample):
#    Totalcount = 0
#    countDic = {}
#    for sample_ in Sample:
#        sample = []
#        sample.append(sample_)
#        print 'sample is', sample
#        count = getEventCount(sample)
#        print 'count is', count
#        Totalcount += count
#        countDic[sample] = count
#    for sample in Sample:
#        print 'weight for ', sample, 'is', countDic[sample]/Totalcount

def getExtWeights(config, extParts):

    sysOut = config.get('Directories','SYSout').strip()
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
#Specialweight
##############

ZLLIncl = ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1","DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext2"]

ZLLjetsHT0       = ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
#ZLLjetsHT70      = ["DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
ZLLjetsHT100     = ["DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
ZLLjetsHT200     = ["DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
ZLLjetsHT400     = ["DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
ZLLjetsHT600     = ["DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
ZLLjetsHT800     = ["DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
ZLLjetsHT1200    = ["DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
ZLLjetsHT2500    = ["DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]

ZlljetsHTbinned = [
        "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1",
#        "DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1",
        "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1",
        "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1",
        "DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]

##V25 only
#ZLLBjets         = ["DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
##V25 + v25b
ZLLBjets         = ["DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1","DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
#
SampleCuts = {"DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":"lheV_pt<100","DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1":"lheHT<100"}

ZLLBjetsVPT0 = ["DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]
ZLLBjetsVPT100= ["DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]
ZLLBjetsVPT200 = ["DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"]

ggZH = ["ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8_ext1","ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8_ext2"]
ZH = ["ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8","ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8_ext1"]

####
#Phase-space cuts
####

DYBJets           = "(lheNb>0)"

HT0               = "(lheHT<100)"
#HT70              = "(lheHT>70&& lheHT<100)"
HT100             = "(lheHT>100 && lheHT<200)"
HT200             = "(lheHT>200 && lheHT<400)"
HT400             = "(lheHT>400 && lheHT<600)"
HT600             = "(lheHT>600 && lheHT<800)"
HT800             = "(lheHT>800 && lheHT<1200)"
HT1200            = "(lheHT>1200 && lheHT<2500)"
HT2500            = "(lheHT>2500)"

VPT0              = "(lheV_pt<100)"
VPT100            = "(lheV_pt>100 && lheV_pt<200)"
VPT200            = "(lheV_pt>200)"


def runInParallel(func, arglist):
    proc = []
    for arg_ in arglist:
        #p = Process(target=getStichWeight, args=arg_)
        p = Process(target=func, args=arg_)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

#arglist = [
#    ('HT0'   ,ZLLjetsHT0   , ZLLBjets, HT0   +"&&"+DYBJets,),
##    ('HT70'  ,ZLLjetsHT70  , ZLLBjets, HT70  +"&&"+DYBJets,),
#    ('HT100' ,ZLLjetsHT100 , ZLLBjets, HT100 +"&&"+DYBJets,),
#    ('HT200' ,ZLLjetsHT200 , ZLLBjets, HT200 +"&&"+DYBJets,),
#    ('HT400' ,ZLLjetsHT400 , ZLLBjets, HT400 +"&&"+DYBJets,),
#    ('HT600' ,ZLLjetsHT600 , ZLLBjets, HT600 +"&&"+DYBJets,),
#    ('HT800' ,ZLLjetsHT800 , ZLLBjets, HT800 +"&&"+DYBJets,),
#    ('HT1200',ZLLjetsHT1200, ZLLBjets, HT1200+"&&"+DYBJets,),
#    ('HT2500',ZLLjetsHT2500, ZLLBjets, HT2500+"&&"+DYBJets,)
#        ]
##runInParallel(arglist)
#runInParallel(getStichWeight,arglist)

#arglist = [
#    ('VPT0'  ,ZlljetsHTbinned, ZLLBjets,VPT0   +"&&"+DYBJets,),
#    ('VPT100',ZlljetsHTbinned, ZLLBjets,VPT200 +"&&"+DYBJets,),
#    ('VPT200',ZlljetsHTbinned, ZLLBjets,VPT100 +"&&"+DYBJets,),
#        ]
#runInParallel(getStichWeight,arglist)

#arglist = [
##    (config, ZLLjetsHT0,),
##    (config, ZLLjetsHT70,),
#    (config, ZLLjetsHT100,),
#    (config, ZLLjetsHT200,),
#    (config, ZLLjetsHT400,)
##    (config, ZLLjetsHT600,),
##    (config, ZLLjetsHT800,),
##    (config, ZLLjetsHT1200,),
##    (config, ZLLjetsHT2500,)
#        ]
#runInParallel(getExtWeights, arglist)

#arglist = [
#    (config, ZLLBjetsVPT100,),
#    (config, ZLLBjetsVPT200,)
#        ]
#runInParallel(getExtWeights, arglist)

#arglist = [
#    (config, ggZH,),
#    (config, ZH,)
#        ]
#runInParallel(getExtWeights, arglist)

arglist = [(config, ZLLIncl,)]
runInParallel(getExtWeights, arglist)

#print 'Weights for HT 0 are',    getStichWeight(ZLLjetsHT0   , ZLLBjets, HT0   +"&&"+DYBJets)
#print 'Weights for HT 70 are',   getStichWeight(ZLLjetsHT70  , ZLLBjets, HT70  +"&&"+DYBJets)
#print 'Weights for HT 100 are',  getStichWeight(ZLLjetsHT100 , ZLLBjets, HT100 +"&&"+DYBJets)
#print 'Weights for HT 200 are',  getStichWeight(ZLLjetsHT200 , ZLLBjets, HT200 +"&&"+DYBJets)
#print 'Weights for HT 400 are',  getStichWeight(ZLLjetsHT400 , ZLLBjets, HT400 +"&&"+DYBJets)
#print 'Weights for HT 600 are',  getStichWeight(ZLLjetsHT600 , ZLLBjets, HT600 +"&&"+DYBJets)
#print 'Weights for HT 800 are',  getStichWeight(ZLLjetsHT800 , ZLLBjets, HT800 +"&&"+DYBJets)
#print 'Weights for HT 1200 are', getStichWeight(ZLLjetsHT1200, ZLLBjets, HT1200+"&&"+DYBJets)
#print 'Weights for HT 2500 are', getStichWeight(ZLLjetsHT2500, ZLLBjets, HT2500+"&&"+DYBJets)


#for section in config.sections():
#    #try:
#    #    sampleName = config.get(section, 'sampleName')
#    #except:
#    sampleName = section
#    sampleNameShort = sampleName
#    if '_ext' in sampleName:
#        sampleNameShort = sampleName.split('_ext')[0]
#    elif '_backup' in sampleName:
#        sampleNameShort = sampleName.split('_backup')[0]
#    if sampleNameShort in sampleDict:
#        sampleDict[sampleNameShort].append(sampleName)
#    else:
#        sampleDict[sampleNameShort] = [sampleName]
#    if verify:
#        extweight = 0
#        try:
#            extweight = float(config.get(section, 'extweight'))
#        except:
#            pass
#        sampleWeights[sampleNameShort] = float(sampleWeights[sampleNameShort]) + extweight if sampleNameShort in sampleWeights else extweight
#
#if verify:
#    for sampleNameShort,totalWeight in sampleWeights.iteritems():
#        if sampleNameShort in sampleDict and len(sampleDict[sampleNameShort])>1:
#            print sampleNameShort,":",totalWeight
#else:
#    for sample,extParts in sampleDict.iteritems():
#        if len(extParts) > 1:
#            print '-'*80
#            print sample,":"
#            getExtWeights(config,extParts)
#
