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

//TH1D* getMC_hist(TString inputfile, TH1D *h1);
TH1D* getdata_hist(TString mbxsec);
//void plot_ratio(TH1D* data_pu,TH1D* mc_pu_1, TH1D* mc_pu_2, TH1D* mc_pu_3, TH1D* mc_pu_4, TH1D* data_pu_5, TH1D* data_pu_6, TString s1, TString s2, TString s3, TString s4, TString s5, TString s6);
void plot_660(TH1D* mc_pu,TH1D* mc_pu_rwt_660, TH1D* data_pu_660);


TString s1 = "/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src/Xbb/python/Wlv2018config/PU_hist/MyDataPileupHistogram_mbxsec_";
TString s2 = ".root";
TString SR_hpt = "Zll_BDT_highpt_";
TString IR_F = "INC_F_";
TString IR_B = "INC_B_";
TString IR_FB = "IR_FB_";
TString VHF_hpt = "Zll_CRZb_highpt_";

//TString region = VHF_hpt;
//TString region = SR_hpt;
//TString region = IR_F;
TString region = IR_B;
//TString region = IR_FB;

void dy_study()
{   

    //TH1D *data_pu_300 = getdata_hist("30000"); 
    //21jun-nob
    TFile* file1 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-nob/Plots/"+ region +"_LHEVpt_FB_.shapes.root");
    //file2->ls();

    TFile* file2 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-18july-withb/Plots/"+region+"_LHEVpt_FB_.shapes.root");
    //TFile* file2 = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src//Xbb/python/logs_Zll2018//runplot-21jun-withbcorr/Plots/"+ region +"_LHEVpt_FB_.shapes.root");
    //file3->ls();


                //*(&h1) = (TH1D*)file2->Get("pileup");
    TH1D *h1 = (TH1D*)file1->Get("summedMcHistograms");
    TH1D *h2 = (TH1D*)file2->Get("summedMcHistograms");
    TH1D *ratio_1 = new TH1D("ratio","HT/DYB",80,0,800);
    //TH1D *ratio_1 = new TH1D("ratio","ratio of corr_with by w/o b-enriched DY samples",80,0,800);
    //TH1D *ratio_2 = new TH1D("ratio","ratio",80,0,800);
    ratio_1->Sumw2();
    ratio_1->Divide(h1,h2);
    ratio_1->GetYaxis()->SetRangeUser(0.0,3.0);
    TF1 *func = new TF1("func","pol2",100,800);
    ratio_1->Fit(func);
    double bin_width = h1->GetXaxis()->GetBinWidth(0);
    std::ostringstream strs;
    strs << bin_width;
    std::string str = strs.str();
//    str.erase (str.find_last_not_of('0') + 1, std::string::npos);
    TString yaxis = "events/" + str;
    ratio_1->GetYaxis()->SetTitle(yaxis);
    ratio_1->GetXaxis()->SetTitle("LHE_Vpt");
    //ratio_1->GetYaxis()->SetTitle("Entries/"+binwidth);
    //ratio_1->Divide(h2,h1);    
    ratio_1->Draw();

    TAxis *xaxis = h1->GetXaxis();
    TAxis *xaxis2 = h2->GetXaxis();
//TAxis *yaxis = h->GetYaxis();
    //Int_t binx = xaxis->FindBin(100);
    //Double_t sum_100_200;
    //Double_t sum_200_400;
    Double_t sum_100_200 = h1->Integral(xaxis->FindBin(100),xaxis->FindBin(200));
    Double_t sum_200_400 = h1->Integral(xaxis->FindBin(200),xaxis->FindBin(400));
    Double_t sum2_100_200 = h2->Integral(xaxis2->FindBin(100),xaxis2->FindBin(200));
    Double_t sum2_200_400 = h2->Integral(xaxis2->FindBin(200),xaxis2->FindBin(400));

    cout<<sum_100_200/sum2_100_200<<endl;
    cout<<sum_200_400/sum2_200_400<<endl;
//
/*
    TCanvas *closure_660 = new TCanvas("closure_660","closure_660",800,600);
    mc_pu_rwt_660->SetLineColor(kBlue);
    data_pu_660->SetLineColor(kRed);
    mc_pu->SetLineColor(kBlack);
    mc_pu_rwt_660->Draw("HIST");
    data_pu_660->Draw("E SAME");
    mc_pu->Draw("HIST SAME");
    closure_660->Draw();
*/

//    plot_ratio(mc_pu, data_pu_600, data_pu_620, data_pu_640, data_pu_660, data_pu_680, data_pu_670, "_60_0", "_68_0", "_69_0", "_70_0", "_71_0", "_72_0");

   /* TH1D* ndata_pu_660 = data_pu_660;
    TH1D* nmc_pu_rwt_660 = mc_pu_rwt_660;
    TH1D* nmc_pu = mc_pu;
  
    plot_660(nmc_pu,nmc_pu_rwt_660,ndata_pu_660);
*/
    //plot_ratio(mc_pu, data_pu_600, data_pu_620, data_pu_640, data_pu_660, data_pu_680, data_pu_692, "_60_0", "_62_0", "_64_0", "_66_0", "_68_0", "_69_2");
    //plot_ratio(mc_pu, data_pu_700, data_pu_720, data_pu_740, data_pu_760, data_pu_780, data_pu_800, "_70_0", "_72_0", "_74_0", "_76_0", "_78_0", "_80_0");
}
/*
TH1D* getdata_hist(TString mbxsec)
{
    TFile* file2 = new TFile(s1+mbxsec+s2);
    file2->ls();

    //(&h1) = (TH1D*)file2->Get("pileup");
    TH1D *h1 = (TH1D*)file2->Get("pileup");
    h1->Scale(1/h1->Integral());

    TFile* file_data_pu = new TFile("/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src/Xbb/python/Wlv2018config/pu_hist/data_pu_" + mbxsec + ".root","RECREATE");
    h1->Write();
    file_data_pu->Close();
    file_data_pu->Delete();

    return h1;
}
    g->Clear();
    g_up->Clear();
    g_down->Clear();




    TCanvas *g1 = new TCanvas("g1","g1",800,600);
    TH1D *data_pu = data->ProjectionX("pileup","e");
    data_pu->Draw("E");
    g1->Draw();

    TH1D *data_pu_up = data_up->ProjectionX("pileup","e");
    TH1D *data_pu_down = data_down->ProjectionX("pileup","e");    


    TFile* file_data_pu = new TFile("data_pu.root","RECREATE");
    data_pu->Write();
    file_data_pu->Close();
    file_data_pu->Delete();

    TFile* file_data_pu_up = new TFile("data_pu_up.root","RECREATE");
    data_pu_up->Write();
    file_data_pu_up->Close();
    file_data_pu_up->Delete();

    TFile* file_data_pu_down = new TFile("data_pu_down.root","RECREATE");
    data_pu_down->Write();
    file_data_pu_down->Close();
    file_data_pu_down->Delete();




//    TCanvas *g2 = new TCanvas("g2","g2",800,600) ;
//    TH1D *mc_pu = mc->ProjectionX();
*/

/*    g2->Divide(1,2);
    g2->cd(1);
    gPad->SetPad(0,0.33,1,0.9378806);
    data_pu->Draw();
    mc_pu->Draw("same");
    g2->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    ratio->Draw("hist");
    g2->Draw();   
*/

/*
    double r[100];
    TH1D *ratio = new TH1D("ratio","ratio",100,0,100);
    ratio->GetYaxis()->SetRangeUser(0,2);    
    for(int i=0; i<100; i++)
    {
	r[i] = data_rho[i]/y[i];
        ratio->Fill(i+0.5,r[i]);
        cout<<i<<" for value "<<r[i]<<endl;
    }
*/	
/*    
    TH1D *mc_pu = mc->ProjectionX("pileup","e");

    TFile* file_mc_pu = new TFile("mc_pu.root","RECREATE");
    mc_pu->Write();
    file_mc_pu->Close();
    file_mc_pu->Delete();



    TCanvas *g3 = new TCanvas("g3","g3",800,600) ;
    TH1D *ratio = new TH1D("ratio","ratio",100,0,100);
    ratio->Sumw2();
    ratio->Divide(data_pu,mc_pu);
    ratio->Draw("E");
    cout<<ratio->GetBinError(20);
     
   //data_pu->Divide(mc_pu); 
    //data_pu->Draw();
    //g3->Draw();

*/


void plot_660(TH1D* nmc_pu,TH1D* nmc_pu_rwt_660, TH1D* ndata_pu_660) 
{
    TCanvas *g1 = new TCanvas("g1","g1",800,600) ;
    g1->Range(0,0,1,1);
    g1->SetFillColor(0);
    g1->SetBorderMode(0);
    g1->SetBorderSize(2);
    g1->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    g1->Divide(1,2);
    g1->cd(1);
    gPad->Range(-69.03766,-0.0002509835,628.8703,0.06681207);
    gPad->SetGridx();
    gPad->SetTopMargin(0.006574267);
    gPad->SetBottomMargin(0.001912046);
    ndata_pu_660->SetLineColor(kBlack);
    ndata_pu_660->SetLineStyle(7);
    ndata_pu_660->SetLineWidth(3);    
    nmc_pu_rwt_660->SetLineColor(kRed);
    nmc_pu->SetLineColor(kGreen+2);
    gPad->SetPad(0,0.33,1,0.9378806);
       
    //mc_pu->Sumw2();
    //data_pu_1->Sumw2();
    THStack* hs = new THStack("hs","Closure test for 66.0 mb min bias xsec of data");
    ndata_pu_660->SetLineWidth(3);
    //data_pu_660->Sumw2();
    hs->Add(nmc_pu,"HIST");
    hs->Add(ndata_pu_660,"HIST");
    hs->Add(nmc_pu_rwt_660,"HIST");
    hs->Draw("nostack");
    TText T; T.SetTextFont(42); T.SetTextAlign(21);
    //T.DrawTextNDC(.5,.95,"Default drawing option");
    //hs->GetHistogram()->SetTitle("Closure test for data_nPU = 66.0 mb");
    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.03);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* la = new TLegend(0.6516547,0.7688758,0.857482,1.12674,NULL,"brNDC");
    TLegend* l1a;
    l1a = la;
    l1a->AddEntry(nmc_pu,"MC_nPU_66.0","l");
    l1a->AddEntry(nmc_pu_rwt_660,"MC_nPU_rwt_66.0","l");
    l1a->AddEntry(ndata_pu_660,"data_nPU_66.0","l");
    l1a->SetTextFont(12);
    l1a->SetTextSize(0.06237499);
    l1a->SetBorderSize(0);
    l1a->SetFillColor(0);
    l1a->Draw("same");



    g1->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3536693);
    TH1D *ratio_1 = new TH1D("ratio","ratio",100,0,100);
    TH1D *ratio_2 = new TH1D("ratio","ratio",100,0,100);    
    ratio_1->Divide(ndata_pu_660,nmc_pu);
    ratio_2->Divide(ndata_pu_660,nmc_pu_rwt_660);
    ratio_1->SetStats(0);
    ratio_1->SetLineColor(kGreen+2);
    ratio_2->SetLineColor(kRed);
    ratio_1->SetTitle("");
    ratio_1->GetXaxis()->SetLabelSize(0.12);
    ratio_1->GetXaxis()->SetTitle("nPU");
    ratio_1->GetXaxis()->SetTitleSize(0.15);
    ratio_1->GetXaxis()->SetTitleOffset(0.99);
    ratio_1->GetXaxis()->SetTickLength(0.03);
    ratio_1->GetYaxis()->SetTitle("data/MC");
    ratio_1->GetYaxis()->CenterTitle();
    ratio_1->GetYaxis()->SetTitleSize(0.115);
    ratio_1->GetYaxis()->SetLabelSize(0.105);
    ratio_1->GetYaxis()->SetTitleOffset(0.35);
    ratio_1->GetYaxis()->SetNdivisions(406);
    ratio_1->GetYaxis()->SetRangeUser(0.5,1.5);
    //ratio->Sumw2();
    ratio_1->Draw("HIST");
    ratio_2->Draw("HIST SAME");
   // ratio->Sumw2();

    Double_t low = ratio_1->GetXaxis()->GetXmin();
    Double_t high = ratio_1->GetXaxis()->GetXmax();
    TLine *line = new TLine(low, 1, high, 1);
    line->SetLineColor(1);
    line->SetLineStyle(7);
    line->SetLineWidth(3);
    line->Draw();
    gPad->Modified();
    g1->Update();
    g1->Draw();

}
