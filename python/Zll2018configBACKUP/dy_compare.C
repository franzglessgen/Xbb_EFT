#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <vector>
#include "RooPolynomial.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFormulaVar.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include "TColor.h"
#include "TAxis.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TLegend.h>
#include "cmath" 


using namespace RooFit;
using namespace RooStats;
using namespace std;

void plot_ratio(TString title, TH1D* hist1,TH1D* hist2, TString s1, TString s2);


TString s1 = "/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src/Xbb/python/Wlv2018config/PU_hist/MyDataPileupHistogram_mbxsec_";
TString s2 = ".root";

void dy_compare()
{   
    TFile* krunal_F = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-21jun-withb/Plots/IR_F__LHEVpt_FB_.shapes.root");
    TFile* krunal_B = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-21jun-withb/Plots/IR_B__LHEVpt_FB_.shapes.root");

    TFile* krunal_noF = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-21jun-nob/Plots/IR_F__LHEVpt_FB_.shapes.root");
    TFile* krunal_noB = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-21jun-nob/Plots/IR_B__LHEVpt_FB_.shapes.root");

    TH1D *kF = (TH1D*)krunal_F->Get("summedMcHistograms");
    TH1D *kB = (TH1D*)krunal_B->Get("summedMcHistograms");

    TH1D *knoF = (TH1D*)krunal_noF->Get("summedMcHistograms");
    TH1D *knoB = (TH1D*)krunal_noB->Get("summedMcHistograms");


    TFile* pirmin_F = new TFile("/mnt/t3nfs01/data01/shome/berger_p2/VHbb/CMSSW_10_1_0/src//Xbb/python/logs_Zll2017//plots_V11_comp2018_Nb0_BGenFilter/Plots//Inclusive__LHE_Vpt_.shapes.root");
    TFile* pirmin_B = new TFile("/mnt/t3nfs01/data01/shome/berger_p2/VHbb/CMSSW_10_1_0/src//Xbb/python/logs_Zll2017//plots_V11_comp2018_NbGt0_DYB/Plots//Inclusive__LHE_Vpt_.shapes.root");

    TFile* pirmin_noF = new TFile("/mnt/t3nfs01/data01/shome/berger_p2/VHbb/CMSSW_10_1_0/src//Xbb/python/logs_Zll2017//plots_V11_comp2018_Nb0/Plots//Inclusive__LHE_Vpt_.shapes.root");
    TFile* pirmin_noB = new TFile("/mnt/t3nfs01/data01/shome/berger_p2/VHbb/CMSSW_10_1_0/src//Xbb/python/logs_Zll2017//plots_V11_comp2018_NbGt0/Plots//Inclusive__LHE_Vpt_.shapes.root");

    TH1D *pF = (TH1D*)pirmin_F->Get("summedMcHistograms");
    TH1D *pB = (TH1D*)pirmin_B->Get("summedMcHistograms");
   
    TH1D *pnoF = (TH1D*)pirmin_noF->Get("summedMcHistograms");
    TH1D *pnoB = (TH1D*)pirmin_noB->Get("summedMcHistograms");
/*
    plot_ratio("bGenFilter",kF,pF,"2017","2018");   
    plot_ratio("DYB",kB,pB,"2017","2018");

    plot_ratio("no_bGenFilter",knoF,pnoF,"2017","2018");
    plot_ratio("no_DYB",knoB,pnoB,"2017","2018");
*/
    //plot_ratio("kb",kB,kB,"2017","2018");
    
    TFile* krunal_F1 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-withb/Plots/INC_F__LHEVpt_FB_.shapes.root");
    TFile* krunal_B1 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-withb/Plots/INC_B__LHEVpt_FB_.shapes.root");

    TFile* krunal_noF1 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-nob/Plots/INC_F__LHEVpt_FB_.shapes.root");
    TFile* krunal_noB1 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-nob/Plots/INC_B__LHEVpt_FB_.shapes.root");

    TH1D *kF1 = (TH1D*)krunal_F1->Get("summedMcHistograms");
    TH1D *kB1 = (TH1D*)krunal_B1->Get("summedMcHistograms");

    TH1D *knoF1 = (TH1D*)krunal_noF1->Get("summedMcHistograms");
    TH1D *knoB1 = (TH1D*)krunal_noB1->Get("summedMcHistograms");

    plot_ratio(kF1,knoF1,"HT","bGenFilter");

}

void plot_ratio(TString title, TH1D* hist1,TH1D* hist2, TString s1, TString s2)
{
    TCanvas *g2 = new TCanvas(title,title,800,600) ;
    g2->Range(0,0,1,1);
    g2->SetFillColor(0);
    g2->SetBorderMode(0);
    g2->SetBorderSize(2);
    g2->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    g2->Divide(1,2);
    g2->cd(1);
    gPad->Range(-69.03766,-0.0002509835,628.8703,0.06681207);
    gPad->SetGridx();
    gPad->SetTopMargin(0.006574267);
    gPad->SetBottomMargin(0.001912046);
    hist1->SetLineColor(kBlack);
    hist2->SetLineColor(kRed);
    gPad->SetPad(0,0.33,1,0.9378806);
       
    THStack* hs = new THStack(title,title);
    hist1->SetFillColor(kWhite);
    hist2->SetFillColor(kWhite);
    hist1->Scale(1/hist1->Integral());
    hist2->Scale(1/hist2->Integral());
    hs->Add(hist1,"HIST");
    hs->Add(hist2,"HIST");
    hs->Draw("nostack");
    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.03);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* la = new TLegend(0.6516547,0.7688758,0.857482,1.292674,NULL,"brNDC");
    TLegend* l1a;
    l1a = la;
    l1a->AddEntry(hist1,"2017","l");
    l1a->AddEntry(hist2,"2018","l");
    l1a->SetTextFont(12);
    l1a->SetTextSize(0.06237499);
    l1a->SetBorderSize(0);
    l1a->SetFillColor(0);
    l1a->Draw("same");



    g2->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3536693);
    TH1D *ratio_1 = new TH1D("ratio","ratio",80,0,800);
    ratio_1->Sumw2();
    ratio_1->Divide(hist1,hist2);
    ratio_1->SetStats(0);
    ratio_1->SetLineColor(kRed);
    ratio_1->SetTitle("");
    ratio_1->GetXaxis()->SetLabelSize(0.12);
    ratio_1->GetXaxis()->SetTitle("LHE_Vpt");
    ratio_1->GetXaxis()->SetTitleSize(0.15);
    ratio_1->GetXaxis()->SetTitleOffset(0.99);
    ratio_1->GetXaxis()->SetTickLength(0.03);
    ratio_1->GetYaxis()->SetTitle("2017/2018");
    ratio_1->GetYaxis()->CenterTitle();
    ratio_1->GetYaxis()->SetTitleSize(0.115);
    ratio_1->GetYaxis()->SetLabelSize(0.105);
    ratio_1->GetYaxis()->SetTitleOffset(0.35);
    ratio_1->GetYaxis()->SetNdivisions(406);
    ratio_1->GetYaxis()->SetRangeUser(0.5,1.5);
    //ratio->Sumw2();
    ratio_1->Draw("E");
   // ratio->Sumw2();

    Double_t low = ratio_1->GetXaxis()->GetXmin();
    Double_t high = ratio_1->GetXaxis()->GetXmax();
    TLine *line = new TLine(low, 1, high, 1);
    line->SetLineColor(1);
    line->SetLineStyle(7);
    line->Draw();
    gPad->Modified();
    g2->Update();
    g2->Draw();

}

