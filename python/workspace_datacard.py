#!/usr/bin/env python
import os
import sys
import ROOT 
from ROOT import TFile
from array import array
from math import sqrt
from copy import copy
#suppres the EvalInstace conversion warning bug
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from BetterConfigParser import BetterConfigParser
from samplesclass import sample
from mvainfos import mvainfo
import pickle
from progbar import progbar
from printcolor import printc
from gethistofromtree import getHistoFromTree, orderandadd


#CONFIGURE
#load config
config = BetterConfigParser()
config.read('./config')
#get locations:
Wdir=config.get('Directories','Wdir')
#systematics
systematics=config.get('systematics','systematics')
systematics=systematics.split(' ')
#TreeVar Array
#Niklas: Not needed?
#MVA_Vars={}
#for systematic in systematics:
#    MVA_Vars[systematic]=config.get('treeVars',systematic)
#    MVA_Vars[systematic]=MVA_Vars[systematic].split(' ')
weightF=config.get('Weights','weightF')
path=sys.argv[1]
var=sys.argv[2]
plot=config.get('Limit',var)
infofile = open(path+'/samples.info','r')
info = pickle.load(infofile)
infofile.close()
options = plot.split(',')
if len(options) < 12:
    print "You have to choose option[11]: either Mjj or BDT"
    sys.exit("You have to choose option[11]: either Mjj or BDT")
name=options[1]
title = options[2]
nBins=int(options[3])
xMin=float(options[4])
xMax=float(options[5])
mass=options[9]
data=options[10]
anType=options[11]
RCut=options[7]
setup=config.get('Limit','setup')
setup=setup.split(',')
ROOToutname = options[6]
outpath=config.get('Directories','limits')
outfile = ROOT.TFile(outpath+'vhbb_TH_'+ROOToutname+'.root', 'RECREATE')


##############################
# MAYBE EDIT THIS:
discr_names = ['ZjLF','ZjHF', 'TT','VV', 's_Top', 'VH', 'WjLF', 'WjHF', 'QCD'] #corresponding to setup
data_name = ['data_obs']
#systematicsnaming={'JER':'CMS_res_j','JES':'CMS_scale_j','beff':'CMS_eff_b','bmis':'CMS_fake_b'}
systematicsnaming={'JER':'CMS_res_j','JES':'CMS_scale_j','beff':'CMS_eff_b','bmis':'CMS_fake_b_7TeV'}
if '8TeV' in options[10]:
        systematicsnaming={'JER':'CMS_res_j','JES':'CMS_scale_j','beff':'CMS_eff_b','bmis':'CMS_fake_b_8TeV'}
#### rescaling by factor 4
scaling=True
if 'RTight' in RCut:
    Datacradbin=options[10]
elif 'RMed' in RCut:
    Datacradbin=options[10]
else:
    Datacradbin=options[10]

    #EDIT!
MC_rescale_factor=1.0
if anaTag=='7TeV': MC_rescale_factor=2.0

#############################

WS = ROOT.RooWorkspace('%s'%Datacradbin,'%s'%Datacradbin) #Zee
print 'WS initialized'
disc= ROOT.RooRealVar(name,name,xMin,xMax)
obs = ROOT.RooArgList(disc)
ROOT.gROOT.SetStyle("Plain")
datas = []
datatyps =[]
histos = []
typs = []
statUps=[]
statDowns=[]
blind=eval(config.get('Limit','blind'))
if blind: 
    print 'I AM BLINDED!'
counter=0
for job in info:
    if job.type == 'BKG':
        #print 'MC'
        hTemp, typ = getHistoFromTree(job,options,MC_rescale_factor)
        histos.append(hTemp)
        typs.append(typ)
        if counter == 0:
            hDummy = copy(hTemp)
        else:
            hDummy.Add(hTemp)
        counter += 1
    elif job.type == 'SIG' and job.name == mass:
        hTemp, typ = getHistoFromTree(job,options,MC_rescale_factor)
        histos.append(hTemp)
        typs.append(typ)    
    elif job.name in data:
        #print 'DATA'
        hTemp, typ = getHistoFromTree(job,options)
        datas.append(hTemp)
        datatyps.append(typ)
    if job.type == 'SIG' and job.name == 'ZH125':
        hSigInjec, typen = getHistoFromTree(job,options,MC_rescale_factor)
        if counter == 0:
            hDummy = copy(hSigInjec)
        else:
            hDummy.Add(hSigInjec)
        counter += 1


MC_integral=0
MC_entries=0
for histo in histos:
    MC_integral+=histo.Integral()
printc('green','', 'MC integral = %s'%MC_integral)  
#order and add together
histos, typs = orderandadd(histos,typs,setup)
rescaleSqrtN = False

for i in range(0,len(histos)):
    histos[i].SetName(discr_names[i])
    #histos[i].SetDirectory(outfile)
    outfile.cd()
    histos[i].Write()
    statUps.append(histos[i].Clone())
    statDowns.append(histos[i].Clone())
    statUps[i].SetName('%sCMS_vhbb_stats_%s_%sUp'%(discr_names[i],discr_names[i],options[10]))
    statDowns[i].SetName('%sCMS_vhbb_stats_%s_%sDown'%(discr_names[i],discr_names[i],options[10]))
    #statUps[i].Sumw2()
    #statDowns[i].Sumw2()
    errorsum=0
    total=0
    for j in range(histos[i].GetNbinsX()+1):
        errorsum=errorsum+(histos[i].GetBinError(j))**2
        
    errorsum=sqrt(errorsum)
    total=histos[i].Integral()
    
    #shift up and down with statistical error
    for j in range(histos[i].GetNbinsX()+1):
	if rescaleSqrtN:
            statUps[i].SetBinContent(j,statUps[i].GetBinContent(j)+statUps[i].GetBinError(j)/total*errorsum)
        else:
            statUps[i].SetBinContent(j,statUps[i].GetBinContent(j)+statUps[i].GetBinError(j))
        if rescaleSqrtN:
            statDowns[i].SetBinContent(j,statDowns[i].GetBinContent(j)-statDowns[i].GetBinError(j)/total*errorsum)
        else:
            statDowns[i].SetBinContent(j,statDowns[i].GetBinContent(j)-statDowns[i].GetBinError(j))
        #statUps[i].SetBinContent(j,statUps[i].GetBinContent(j)+statUps[i].GetBinError(j))
        #statDowns[i].SetBinContent(j,statDowns[i].GetBinContent(j)-statDowns[i].GetBinError(j))

    '''
    ######################
    #trying some crazy shifting:
    anzahlBins=histos[i].GetNbinsX()
    contentarray=[]
    errorarray=[]
    indexarray=[]
    for j in range(0,anzahlBins):
        Ncontent=histos[i].GetBinContent(j)
        Nerror=histos[i].GetBinError(j)
        if Ncontent>0:
            contentarray.append(Ncontent)
            errorarray.append(Nerror)
            indexarray.append(j)
    nonzeroBins=len(contentarray)
    ungerade=nonzeroBins%2
    half=(nonzeroBins-ungerade)/2
    newarray=[0]*nonzeroBins
    if ungerade:
        #factor=-1
        for m in range(0,half):
            newarray[m]=(half-m)*(-1)*errorarray[m]/half
            newarray[m+half+1]=(m)*(+1)*errorarray[m+half+1]/half
        newarray[half+1]=0
    else:
        #factor=-1
        for m in range(0,half):
            newarray[m]=(half-m)*(-1)*errorarray[m]/half
            newarray[m+half]=(m)*(+1)*errorarray[m+half]/half
    for j in range(0,anzahlBins):
        if j in indexarray:
            whereisit=indexarray.index(j)
            statUps[i].SetBinContent(j,contentarray[whereisit]+newarray[whereisit])
            statDowns[i].SetBinContent(j,contentarray[whereisit]-newarray[whereisit])
        else:
            statUps[i].SetBinContent(j,0)
            statDowns[i].SetBinContent(j,0)
    ###################
    '''

    statUps[i].Write()
    statDowns[i].Write()
    histPdf = ROOT.RooDataHist(discr_names[i],discr_names[i],obs,histos[i])
    #UP stats of MCs
    RooStatsUp = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sUp'%(discr_names[i],discr_names[i],options[10]),'%sCMS_vhbb_stats_%s_%sUp'%(discr_names[i],discr_names[i],options[10]),obs, statUps[i])
    #DOWN stats of MCs
    RooStatsDown = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sDown'%(discr_names[i],discr_names[i],options[10]),'%sCMS_vhbb_stats_%s_%sDown'%(discr_names[i],discr_names[i],options[10]),obs, statDowns[i])
    getattr(WS,'import')(histPdf)
    getattr(WS,'import')(RooStatsUp)
    getattr(WS,'import')(RooStatsDown)

#dunnmies - only to fill in empty histos for QCD and Wj
#Wlight,Wbb,QCD
for i in range(6,9):
    dummy = ROOT.TH1F(discr_names[i], 'discriminator', nBins, xMin, xMax)
    outfile.cd()
    dummy.Write()
    #nominal
    histPdf = ROOT.RooDataHist(discr_names[i],discr_names[i],obs,dummy)
    #UP stats of MCs
    RooStatsUp = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sUp'%(discr_names[i],discr_names[i],options[10]),'%sCMS_vhbb_stats_%s_%sUp'%(discr_names[i],discr_names[i],options[10]),obs, dummy)
    #DOWN stats of MCs
    RooStatsDown = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sDown'%(discr_names[i],discr_names[i],options[10]),'%sCMS_vhbb_stats_%s_%sDown'%(discr_names[i],discr_names[i],options[10]),obs, dummy)
    getattr(WS,'import')(histPdf)
    getattr(WS,'import')(RooStatsUp)
    getattr(WS,'import')(RooStatsDown)

#HISTOGRAMM of DATA    
d1 = ROOT.TH1F('d1','d1',nBins,xMin,xMax)
for i in range(0,len(datas)):
    d1.Add(datas[i],1)
printc('green','','\nDATA integral = %s\n'%d1.Integral())
flow = d1.GetEntries()-d1.Integral()
if flow > 0:
    printc('red','','U/O flow: %s'%flow)
d1.SetName(data_name[0])
outfile.cd()
d1.Write()
if blind:
    hDummy.SetName(data_name[0])
    histPdf = ROOT.RooDataHist('data_obs','data_obs',obs,hDummy)
    #rooDummy = ROOT.RooDataHist('data_obs','data_obs',obs,hDummy)
    #toy = ROOT.RooHistPdf('data_obs','data_obs',ROOT.RooArgSet(obs),rooDummy)
    #rooDataSet = toy.generate(ROOT.RooArgSet(obs),int(d1.Integral()))
    #histPdf = ROOT.RooDataHist('data_obs','data_obs',ROOT.RooArgSet(obs),rooDataSet.reduce(ROOT.RooArgSet(obs)))
else:
    histPdf = ROOT.RooDataHist('data_obs','data_obs',obs,d1)
#ROOT.RooAbsData.plotOn(histPdf,frame)
#frame.Draw()
#IMPORT
getattr(WS,'import')(histPdf)

#SYSTEMATICS:
UD = ['Up','Down']
systhistosarray=[]
Coco=0 #iterates over (all systematics) * (up,down)

bdt = False
mjj = False

#print str(anType)
#print len(options)
if str(anType) == 'BDT':
    bdt = True
    systematics = ['JER','JES','beff','bmis']
elif str(anType) == 'Mjj':
    mjj = True
    systematics = ['JER','JES']

nominalShape = options[0]
    
for sys in systematics:
    for Q in UD:
        ff=options[0].split('.')
        if bdt == True:
            ff[1]='%s_%s'%(sys,Q.lower())
            options[0]=nominalShape.replace('.nominal','.%s_%s'%(sys,Q.lower()))
        elif mjj == True:
            ff[0]='H_%s'%(sys)
            ff[1]='mass_%s'%(Q.lower())
            options[0]='.'.join(ff)
        print options[0]


        print '\n'
        printc('blue','','\t--> doing systematic %s %s'%(sys,Q.lower())) 

        systhistosarray.append([])
        typsX = []

        for job in info:
            #print job.name
            if job.type == 'BKG':
                #print 'MC'
                hTemp, typ = getHistoFromTree(job,options,MC_rescale_factor)
                systhistosarray[Coco].append(hTemp)
                typsX.append(typ)
            elif job.type == 'SIG' and job.name == mass:
                #print 'MC'
                hTemp, typ = getHistoFromTree(job,options,MC_rescale_factor)
                systhistosarray[Coco].append(hTemp)
                typsX.append(typ)

        MC_integral=0
        for histoX in systhistosarray[Coco]:
            MC_integral+=histoX.Integral()
        printc('green','', 'MC integral = %s'%MC_integral)  
        systhistosarray[Coco], typsX = orderandadd(systhistosarray[Coco],typsX,setup)
        '''
        # do the linear fit blabla
        for i in range(0,len(systhistosarray[Coco])):
            #systhistosarray[Coco][i]
            #histos[i]
            for bin in range(0,histos[i].GetSize()):
                A=systhistosarray[Coco][i].GetBinContent(bin)
                B=histos[i].GetBinContent(bin)
                systhistosarray[Coco][i].SetBinContent(bin,A-B)
            #Fit:
            FitFunction=ROOT.TF1('FitFunction','pol1')
            systhistosarray[Coco][i].Fit('FitFunction')
        '''
        if scaling:
            #or now i try some rescaling by 4:
            for i in range(0,len(systhistosarray[Coco])):
                #systhistosarray[Coco][i]
                #histos[i]
                for bin in range(0,histos[i].GetSize()):
                    A=systhistosarray[Coco][i].GetBinContent(bin)
                    B=histos[i].GetBinContent(bin)
                    systhistosarray[Coco][i].SetBinContent(bin,B+((A-B)/4.))
        # finaly lpop over histos
        for i in range(0,len(systhistosarray[Coco])):
            systhistosarray[Coco][i].SetName('%s%s%s'%(discr_names[i],systematicsnaming[sys],Q))
            outfile.cd()
            systhistosarray[Coco][i].Write()            
            histPdf = ROOT.RooDataHist('%s%s%s'%(discr_names[i],systematicsnaming[sys],Q),'%s%s%s'%(discr_names[i],systematicsnaming[sys],Q),obs,systhistosarray[Coco][i])
            getattr(WS,'import')(histPdf)
        Coco+=1
        #print Coco
WS.writeToFile(outpath+'vhbb_WS_'+ROOToutname+'.root')
   #WS.writeToFile("testWS.root")
   
   
   

#write DATAcard:
if '8TeV' in options[10]:
    pier = open(Wdir+'/pier8TeV.txt','r')
else:
    pier = open(Wdir+'/pier.txt','r')
scalefactors=pier.readlines()
pier.close()
f = open(outpath+'vhbb_DC_'+ROOToutname+'.txt','w')
f.write('imax\t1\tnumber of channels\n')
f.write('jmax\t8\tnumber of backgrounds (\'*\' = automatic)\n')
f.write('kmax\t*\tnumber of nuisance parameters (sources of systematical uncertainties)\n\n')
if bdt==True:
    f.write('shapes * * vhbb_WS_%s.root $CHANNEL:$PROCESS $CHANNEL:$PROCESS$SYSTEMATIC\n\n'%ROOToutname)
else:
    f.write('shapes * * vhbb_TH_%s.root $PROCESS $PROCESS$SYSTEMATIC\n\n'%ROOToutname)
f.write('bin\t%s\n\n'%Datacradbin)
if blind:
    f.write('observation\t%s\n\n'%(hDummy.Integral()))
else:
    f.write('observation\t%s\n\n'%(int(d1.Integral())))
f.write('bin\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(Datacradbin,Datacradbin,Datacradbin,Datacradbin,Datacradbin,Datacradbin,Datacradbin,Datacradbin,Datacradbin))
f.write('process\tVH\tWjLF\tWjHF\tZjLF\tZjHF\tTT\ts_Top\tVV\tQCD\n')

f.write('process\t0\t1\t2\t3\t4\t5\t6\t7\t8\n')
f.write('rate\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(histos[5].Integral(),0,0,histos[0].Integral(),histos[1].Integral(),histos[2].Integral(),histos[4].Integral(),histos[3].Integral(),0)) #\t1.918\t0.000 0.000\t135.831  117.86  18.718 1.508\t7.015\t0.000
if '7TeV' in options[10]:
    f.write('lumi_7TeV\tlnN\t1.022\t-\t-\t-\t-\t-\t1.022\t1.022\t1.022\n')
if '8TeV' in options[10]:
    f.write('lumi_8TeV\tlnN\t1.05\t-\t-\t-\t-\t-\t1.05\t1.05\t1.05\n')
f.write('pdf_qqbar\tlnN\t1.01\t-\t-\t-\t-\t-\t-\t1.01\t-\n')
f.write('pdf_gg\tlnN\t-\t-\t-\t-\t-\t-\t1.01\t-\t1.01\n')
f.write('QCDscale_VH\tlnN\t1.04\t-\t-\t-\t-\t-\t-\t-\t-\n')
f.write('QCDscale_ttbar\tlnN\t-\t-\t-\t-\t-\t-\t1.06\t-\t-\n')
f.write('QCDscale_VV\tlnN\t-\t-\t-\t-\t-\t-\t-\t1.04\t-\n')
f.write('QCDscale_QCD\tlnN\t-\t-\t-\t-\t-\t-\t-\t-\t1.30\n')
f.write('CMS_vhbb_boost_EWK\tlnN\t1.05\t-\t-\t-\t-\t-\t-\t-\t-\n')
f.write('CMS_vhbb_boost_QCD\tlnN\t1.10\t-\t-\t-\t-\t-\t-\t-\t-\n')
f.write('CMS_vhbb_ST\tlnN\t-\t-\t-\t-\t-\t-\t1.29\t-\t-\n')
f.write('CMS_vhbb_VV\tlnN\t-\t-\t-\t-\t-\t-\t-\t1.30\t-\n')
if '7TeV' in options[10]:
	f.write('CMS_vhbb_ZjLF_ex\tlnN\t-\t-\t-\t1.05\t-\t-\t-\t-\t-\n')
	f.write('CMS_vhbb_ZjHF_ex\tlnN\t-\t-\t-\t-\t1.05\t-\t-\t-\t-\n')
	f.write('CMS_vhbb_TT_ex\tlnN\t-\t-\t-\t-\t-\t1.05\t-\t-\t-\n')
if '8TeV' in options[10]:
	f.write('CMS_vhbb_ZjLF_ex_8TeV\tlnN\t-\t-\t-\t1.05\t-\t-\t-\t-\t-\n')
	f.write('CMS_vhbb_ZjHF_ex_8TeV\tlnN\t-\t-\t-\t-\t1.05\t-\t-\t-\t-\n')
	f.write('CMS_vhbb_TT_ex_8TeV\tlnN\t-\t-\t-\t-\t-\t1.05\t-\t-\t-\n')
for line in scalefactors:
    f.write(line)
if 'Zee' in options[10]:
    f.write('CMS_eff_m lnN\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
    f.write('CMS_eff_e lnN\t1.04\t-\t-\t-\t-\t-\t1.04\t1.04\t1.04\n')
    #f.write('CMS_trigger_m\tlnN\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
    #f.write('CMS_trigger_e\tlnN\t1.02\t-\t-\t-\t-\t-\t1.02\t1.02\t-\n')
if 'Zmm' in options[10]:
    f.write('CMS_eff_e lnN\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
    f.write('CMS_eff_m lnN\t1.04\t-\t-\t-\t-\t-\t1.04\t1.04\t1.04\n')
    #f.write('CMS_trigger_e\tlnN\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
    #f.write('CMS_trigger_m\tlnN\t1.01\t-\t-\t-\t-\t-\t1.01\t1.01\t-\n')

f.write('CMS_vhbb_trigger_MET\tlnN\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
f.write('CMS_vhbb_stats_%s_%s\tshape\t1.0\t-\t-\t-\t-\t-\t-\t-\t-\n'%(discr_names[5], options[10]))
f.write('CMS_vhbb_stats_%s_%s\tshape\t-\t-\t-\t1.0\t-\t-\t-\t-\t-\n'%(discr_names[0], options[10]))
f.write('CMS_vhbb_stats_%s_%s\tshape\t-\t-\t-\t-\t1.0\t-\t-\t-\t-\n'%(discr_names[1], options[10]))
f.write('CMS_vhbb_stats_%s_%s\tshape\t-\t-\t-\t-\t-\t1.0\t-\t-\t-\n'%(discr_names[2], options[10]))
f.write('CMS_vhbb_stats_%s_%s\tshape\t-\t-\t-\t-\t-\t-\t1.0\t-\t-\n'%(discr_names[4], options[10]))
f.write('CMS_vhbb_stats_%s_%s\tshape\t-\t-\t-\t-\t-\t-\t-\t1.0\t-\n'%(discr_names[3], options[10]))
#SYST
if bdt==True:
    if scaling:
        f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['JER'])
        f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['JES'])
        f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['beff'])
        f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['bmis'])
    else:
        #SYST4
        f.write('%s\tshape\t0.25\t-\t-\t0.25\t0.25\t0.25\t0.25\t0.25\t-\n'%systematicsnaming['JER'])
        f.write('%s\tshape\t0.25\t-\t-\t0.25\t0.25\t0.25\t0.25\t0.25\t-\n'%systematicsnaming['JES'])
        f.write('%s\tshape\t0.25\t-\t-\t0.25\t0.25\t0.25\t0.25\t0.25\t-\n'%systematicsnaming['beff'])
        f.write('%s\tshape\t0.25\t-\t-\t0.25\t0.25\t0.25\t0.25\t0.25\t-\n'%systematicsnaming['bmis'])
else:
    f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['JER'])
    f.write('%s\tshape\t1.0\t-\t-\t1.0\t1.0\t1.0\t1.0\t1.0\t-\n'%systematicsnaming['JES'])
f.close()

outfile.Write()
outfile.Close()
