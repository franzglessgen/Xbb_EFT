#! /usr/bin/env python
from __future__ import print_function
import ROOT
import os

fileNameList = '''
logs_Zvv2017//plots_V11_kfactors_NbGt0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zvv2017//plots_V11_kfactors_NbGt0_DYB/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zvv2017//plots_V11_kfactors_Nb0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zvv2017//plots_V11_kfactors_Nb0_BGenFilter/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Wlv2017//plots_V11_kfactors_NbGt0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Wlv2017//plots_V11_kfactors_NbGt0_DYB/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Wlv2017//plots_V11_kfactors_Nb0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Wlv2017//plots_V11_kfactors_Nb0_BGenFilter/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_V11_kfactors_NbGt0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_V11_kfactors_NbGt0_DYB/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_V11_kfactors_Nb0/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_V11_kfactors_Nb0_BGenFilter/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017/plots_V11_comp2018V3_HT/Plots/Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017/plots_V11_comp2018V3_DYBPlusBGenFilter/Plots/Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017/plots_V11_comp2018V3_HT/Plots/Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017/plots_V11_comp2018V3_inclOnly/Plots//Inclusive__LHE_Vpt_.shapes.root
'''

fileNames  = [x.strip() for x in fileNameList.split('\n') if len(x.strip()) > 0]
rootFiles  = [ROOT.TFile.Open(x) for x in fileNames]
histograms = [f.Get("summedMcHistograms") for f in rootFiles]

print("-"*80)
print("k-factors (to be multiplied with b-enriched samples to get back to yields of HTbinned samples)")
print(" ZBJets:            %1.3f"%(histograms[0].Integral()/histograms[1].Integral()))
print(" ZJets_BGenFilter:  %1.3f"%(histograms[2].Integral()/histograms[3].Integral()))
print(" WBJets:            %1.3f"%(histograms[4].Integral()/histograms[5].Integral()))
print(" WJets_BGenFilter:  %1.3f"%(histograms[6].Integral()/histograms[7].Integral()))
print(" DYBJets:           %1.3f"%(histograms[8].Integral()/histograms[9].Integral()))
print(" DYJets_BGenFilter: %1.3f"%(histograms[10].Integral()/histograms[11].Integral()))
print("-"*80)

def histogramRatio(h1,h2):
    hr = h1.Clone()
    hr.Divide(h2)
    return hr

ratioHistograms     = [histogramRatio(histograms[i],histograms[i+1]) for i in [0,2,4,6,8,10,12,14]]
#ratioHistograms     = [histogramRatio(histograms[i+1],histograms[i]) for i in [0,2,4,6,8,10]]
ratioHistogramNames = ['ZBJets','ZJets_BGenFilter','WBJets','WJets_BGenFilter','DYBJets','DYJets_BGenFilter','DYBJets_plus_DYJets_BGenFilter','DYjets_inclusive']

try:
    os.makedirs('results/benriched/')
except:
    pass

rangeMin = 100
rangeMax = 800

for i in range(len(ratioHistograms)):
    c1 = ROOT.TCanvas("c1","c1",800,500)
    fcn = ROOT.TF1("f1", "pol2", rangeMin, rangeMax)
    r = ratioHistograms[i].Fit(fcn, "RS", "", rangeMin, rangeMax) 
    print(ratioHistogramNames[i])
    print(" chi2/#df =", fcn.GetChisquare()/fcn.GetNDF())
    print((" (%1.3e + %1.3e*min(LHE_Vpt,800) + %1.3e*min(LHE_Vpt,800)**2)"%(r.Parameter(0),r.Parameter(1), r.Parameter(2))).replace('+ -','- '))
    #outputFileName = 'results/benriched/' + ratioHistogramNames[i] + '.png'
    outputFileName = 'results/benriched/' + ratioHistogramNames[i] + '.root'
    c1.SaveAs(outputFileName)



