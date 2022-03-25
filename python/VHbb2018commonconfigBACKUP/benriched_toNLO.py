#! /usr/bin/env python
from __future__ import print_function
import ROOT
import os

ratioHistogramNames = ['DY2b','DY1b','DY0b']
fileNameList = '''
logs_Zll2017//plots_NLO2_NLO/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_NLO2/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_NLO1_NLO/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_NLO1/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_NLO0_NLO/Plots//Inclusive__LHE_Vpt_.shapes.root
logs_Zll2017//plots_NLO0/Plots//Inclusive__LHE_Vpt_.shapes.root
'''
rangeMin = 50
rangeMax = 500


fileNames  = [x.strip() for x in fileNameList.split('\n') if len(x.strip()) > 0]
rootFiles  = [ROOT.TFile.Open(x) for x in fileNames]
histograms = [f.Get("summedMcHistograms") for f in rootFiles]

def histogramRatio(h1,h2,rebin=4):
    hr = h1.Clone()
    h2c = h2.Clone()
    hr.Rebin(rebin)
    hr.Scale(1.0/rebin)
    h2c.Rebin(rebin)
    h2c.Scale(1.0/rebin)
    print(hr.GetXaxis().GetNbins(), h2c.GetXaxis().GetNbins())
    hr.Divide(h2c)
    return hr

ratioHistograms     = [histogramRatio(histograms[2*i],histograms[2*i+1]) for i in range(len(ratioHistogramNames))]

print("-"*80)
print("NLO/LO normalization:")
for i in range(len(ratioHistogramNames)):
    print(" " + ratioHistogramNames[i] + " :            %1.3f"%(histograms[2*i].Integral()/histograms[2*i+1].Integral()))
print("-"*80)

try:
    os.makedirs('results/benriched/')
except:
    pass

for i in range(len(ratioHistograms)):
    c1 = ROOT.TCanvas("c1","c1",800,500)
    
    if True:
        fitModel = "[0]-[1]*(x-%1.3e)+[2]*(x-%1.3e)**2"%(rangeMin, rangeMin)
        fcn = ROOT.TF1("f1", fitModel, rangeMin, rangeMax)
        r = ratioHistograms[i].Fit(fcn, "RS", "", rangeMin, rangeMax) 

        print(ratioHistogramNames[i])
        print(" chi2/#df =", fcn.GetChisquare()/fcn.GetNDF())
        print(fitModel.replace("[0]","%1.4e"%r.Parameter(0)).replace("[1]","%1.4e"%r.Parameter(1)).replace("[2]","%1.4e"%r.Parameter(2)).replace('+ -','- '))

    else:
        fitModel = "[1]-(x-%1.3e)*[0]"%rangeMin
        fcn = ROOT.TF1("f0", fitModel, rangeMin, rangeMax)
        r = ratioHistograms[i].Fit(fcn, "RS", "", rangeMin, rangeMax) 

        fitModel = fitModel.replace("[1]","%1.5e"%r.Parameter(1))
        p1 = r.Parameter(1)
        fcn = ROOT.TF1("f1", fitModel, rangeMin, rangeMax)
        r = ratioHistograms[i].Fit(fcn, "RS", "", rangeMin, rangeMax) 

        print(ratioHistogramNames[i])
        print(" chi2/#df =", fcn.GetChisquare()/fcn.GetNDF())
        print(("%1.3e-(%1.3e*(max(min(LHE_Vpt,%1.3e),%1.3e)-%1.3e))"%(p1, r.Parameter(0),rangeMax,rangeMin,rangeMin)).replace('+ -','- '))

    outputFileName = 'results/benriched/NLO_over_LO_' + ratioHistogramNames[i] + '.root'
    c1.SaveAs(outputFileName)
    ROOT.gPad.Update()
    raw_input()



