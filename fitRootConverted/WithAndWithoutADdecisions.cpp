#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"



void WithAndWithoutADdecisions(){

  TFile* fileListWith    = new TFile("AnalysisResultsLHC18qr_CMUP6_15062020.root");   // with    AD decisions
  TFile* fileListWithout = new TFile("AnalysisResultsLHC18qr_CMUP6_noADflags.root");  // without AD decisions
  TDirectory* dirWith    = fileListWith   ->GetDirectory("MyTask");
  TDirectory* dirWithout = fileListWithout->GetDirectory("MyTask");
  TList* listingsWith;
  TList* listingsWithout;
  dirWith   ->GetObject("ADcheck_", listingsWith);
  dirWithout->GetObject("ADcheck_", listingsWithout);
  TH1F* Pt[4][2];
  Pt[0][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAzeroH");
  Pt[1][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAzeroH");
  Pt[2][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAanyH");
  Pt[3][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAanyH");
  Pt[0][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAzeroH");
  Pt[1][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAzeroH");
  Pt[2][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAanyH");
  Pt[3][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAanyH");
  TH1F* PtTwo[4][2];
  // PtTwo[0][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAzeroH");
  // PtTwo[1][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAzeroH");
  // PtTwo[2][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAanyH");
  // PtTwo[3][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAanyH");
  // PtTwo[0][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAzeroH");
  // PtTwo[1][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAzeroH");
  // PtTwo[2][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAanyH");
  // PtTwo[3][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAanyH");
  PtTwo[0][0] = (TH1F*)Pt[0][0]->Clone("PtTwo_00");
  PtTwo[1][0] = (TH1F*)Pt[1][0]->Clone("PtTwo_10");
  PtTwo[2][0] = (TH1F*)Pt[2][0]->Clone("PtTwo_20");
  PtTwo[3][0] = (TH1F*)Pt[3][0]->Clone("PtTwo_30");
  PtTwo[0][1] = (TH1F*)Pt[0][1]->Clone("PtTwo_01");
  PtTwo[1][1] = (TH1F*)Pt[1][1]->Clone("PtTwo_11");
  PtTwo[2][1] = (TH1F*)Pt[2][1]->Clone("PtTwo_21");
  PtTwo[3][1] = (TH1F*)Pt[3][1]->Clone("PtTwo_31");
  TH1F* PtV3[4][2];
  PtV3[0][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3");
  PtV3[1][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3");
  PtV3[2][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3");
  PtV3[3][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3");
  PtV3[0][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3");
  PtV3[1][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3");
  PtV3[2][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3");
  PtV3[3][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3");
  TH1F* PtTwoV3[4][2];
  // PtTwoV3[0][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3");
  // PtTwoV3[1][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3");
  // PtTwoV3[2][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3");
  // PtTwoV3[3][0] = (TH1F*)listingsWith   ->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3");
  // PtTwoV3[0][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3");
  // PtTwoV3[1][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3");
  // PtTwoV3[2][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3");
  // PtTwoV3[3][1] = (TH1F*)listingsWithout->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3");
  PtTwoV3[0][0] = (TH1F*)PtV3[0][0]->Clone("PtTwoV3_00");
  PtTwoV3[1][0] = (TH1F*)PtV3[1][0]->Clone("PtTwoV3_10");
  PtTwoV3[2][0] = (TH1F*)PtV3[2][0]->Clone("PtTwoV3_20");
  PtTwoV3[3][0] = (TH1F*)PtV3[3][0]->Clone("PtTwoV3_30");
  PtTwoV3[0][1] = (TH1F*)PtV3[0][1]->Clone("PtTwoV3_01");
  PtTwoV3[1][1] = (TH1F*)PtV3[1][1]->Clone("PtTwoV3_11");
  PtTwoV3[2][1] = (TH1F*)PtV3[2][1]->Clone("PtTwoV3_21");
  PtTwoV3[3][1] = (TH1F*)PtV3[3][1]->Clone("PtTwoV3_31");


  for (size_t i = 0; i < 4; i++) {
    Pt[i][0]     ->SetLineWidth(3);
    PtTwo[i][0]  ->SetLineWidth(3);
    PtV3[i][0]   ->SetLineWidth(3);
    PtTwoV3[i][0]->SetLineWidth(3);
    Pt[i][1]     ->SetLineWidth(3);
    PtTwo[i][1]  ->SetLineWidth(3);
    PtV3[i][1]   ->SetLineWidth(3);
    PtTwoV3[i][1]->SetLineWidth(3);
    Pt[i][0]     ->SetLineColor(kRed);
    PtTwo[i][0]  ->SetLineColor(kRed);
    PtV3[i][0]   ->SetLineColor(kRed);
    PtTwoV3[i][0]->SetLineColor(kRed);
    Pt[i][0]     ->Sumw2();
    PtTwo[i][0]  ->Sumw2();
    PtV3[i][0]   ->Sumw2();
    PtTwoV3[i][0]->Sumw2();
    Pt[i][1]     ->Sumw2();
    PtTwo[i][1]  ->Sumw2();
    PtV3[i][1]   ->Sumw2();
    PtTwoV3[i][1]->Sumw2();
    Pt[i][0]     ->GetXaxis()->SetTitleOffset(1.25);
    Pt[i][0]     ->GetYaxis()->SetTitleOffset(1.45);
    Pt[i][0]     ->GetXaxis()->SetTitleSize(0.045);
    Pt[i][0]     ->GetYaxis()->SetTitleSize(0.045);
    Pt[i][0]     ->GetXaxis()->SetLabelSize(0.045);
    Pt[i][0]     ->GetYaxis()->SetLabelSize(0.045);
    Pt[i][0]     ->GetXaxis()->SetTitleFont(42);
    Pt[i][0]     ->GetYaxis()->SetTitleFont(42);
    Pt[i][0]     ->GetXaxis()->SetLabelFont(42);
    Pt[i][0]     ->GetYaxis()->SetLabelFont(42);
    Pt[i][0]     ->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    Pt[i][0]     ->GetYaxis()->SetTitle("Counts");
    PtTwo[i][0]  ->GetXaxis()->SetTitleOffset(1.25);
    PtTwo[i][0]  ->GetYaxis()->SetTitleOffset(1.45);
    PtTwo[i][0]  ->GetXaxis()->SetTitleSize(0.045);
    PtTwo[i][0]  ->GetYaxis()->SetTitleSize(0.045);
    PtTwo[i][0]  ->GetXaxis()->SetLabelSize(0.045);
    PtTwo[i][0]  ->GetYaxis()->SetLabelSize(0.045);
    PtTwo[i][0]  ->GetXaxis()->SetTitleFont(42);
    PtTwo[i][0]  ->GetYaxis()->SetTitleFont(42);
    PtTwo[i][0]  ->GetXaxis()->SetLabelFont(42);
    PtTwo[i][0]  ->GetYaxis()->SetLabelFont(42);
    PtTwo[i][0]  ->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    PtTwo[i][0]  ->GetYaxis()->SetTitle("Ratio");
    PtTwo[i][1]  ->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    PtTwo[i][1]  ->GetYaxis()->SetTitle("Ratio");
    PtV3[i][0]   ->GetXaxis()->SetTitleOffset(1.25);
    PtV3[i][0]   ->GetYaxis()->SetTitleOffset(1.45);
    PtV3[i][0]   ->GetXaxis()->SetTitleSize(0.045);
    PtV3[i][0]   ->GetYaxis()->SetTitleSize(0.045);
    PtV3[i][0]   ->GetXaxis()->SetLabelSize(0.045);
    PtV3[i][0]   ->GetYaxis()->SetLabelSize(0.045);
    PtV3[i][0]   ->GetXaxis()->SetTitleFont(42);
    PtV3[i][0]   ->GetYaxis()->SetTitleFont(42);
    PtV3[i][0]   ->GetXaxis()->SetLabelFont(42);
    PtV3[i][0]   ->GetYaxis()->SetLabelFont(42);
    PtV3[i][0]   ->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    PtV3[i][0]   ->GetYaxis()->SetTitle("Counts");
    PtTwoV3[i][0]->GetXaxis()->SetTitleOffset(1.25);
    PtTwoV3[i][0]->GetXaxis()->SetTitleOffset(1.45);
    PtTwoV3[i][0]->GetXaxis()->SetTitleSize(0.045);
    PtTwoV3[i][0]->GetXaxis()->SetTitleSize(0.045);
    PtTwoV3[i][0]->GetXaxis()->SetLabelSize(0.045);
    PtTwoV3[i][0]->GetXaxis()->SetLabelSize(0.045);
    PtTwoV3[i][0]->GetXaxis()->SetTitleFont(42);
    PtTwoV3[i][0]->GetXaxis()->SetTitleFont(42);
    PtTwoV3[i][0]->GetXaxis()->SetLabelFont(42);
    PtTwoV3[i][0]->GetXaxis()->SetLabelFont(42);
    PtTwoV3[i][0]->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    PtTwoV3[i][0]->GetYaxis()->SetTitle("Ratio");
    PtTwoV3[i][1]->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    PtTwoV3[i][1]->GetYaxis()->SetTitle("Ratio");
    Pt[i][0]     ->Rebin(5);
    PtTwo[i][0]  ->Rebin(5);
    Pt[i][0]     ->GetXaxis()->SetRangeUser(0, 4);
    PtTwo[i][0]  ->GetXaxis()->SetRangeUser(0, 4);
    Pt[i][1]     ->Rebin(5);
    PtTwo[i][1]  ->Rebin(5);
    Pt[i][1]     ->GetXaxis()->SetRangeUser(0, 4);
    PtTwo[i][1]  ->GetXaxis()->SetRangeUser(0, 4);
    PtTwoV3[i][0]->GetYaxis()->SetRangeUser(0.75, 2);
    PtTwoV3[i][1]->GetYaxis()->SetRangeUser(0.75, 2);
    PtTwoV3[i][0]->GetXaxis()->SetRangeUser(0., 0.6);
    PtTwoV3[i][1]->GetXaxis()->SetRangeUser(0., 0.6);
    PtV3[i][0]   ->GetXaxis()->SetRangeUser(0., 0.6);
    PtV3[i][1]   ->GetXaxis()->SetRangeUser(0., 0.6);

  }


  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);



  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(4,2);
  c1->cd(1);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  Pt[0][0]->Draw();
  Pt[0][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c1->cd(2);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  Pt[1][0]->Draw();
  Pt[1][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c1->cd(3);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  Pt[2][0]->Draw();
  Pt[2][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c1->cd(4);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  Pt[3][0]->Draw();
  Pt[3][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c1->cd(5);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwo[0][0]->Divide(PtTwo[0][1]);
  PtTwo[0][0]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c1->cd(6);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwo[1][0]->Divide(PtTwo[1][1]);
  PtTwo[1][0]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c1->cd(7);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwo[2][0]->Divide(PtTwo[2][1]);
  PtTwo[2][0]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c1->cd(8);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwo[3][0]->Divide(PtTwo[3][1]);
  PtTwo[3][0]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");



  TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->Divide(4,2);
  c2->cd(1);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  PtV3[0][0]->Draw();
  PtV3[0][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c2->cd(2);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  PtV3[1][0]->Draw();
  PtV3[1][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c2->cd(3);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  PtV3[2][0]->Draw();
  PtV3[2][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c2->cd(4);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  PtV3[3][0]->Draw();
  PtV3[3][1]->Draw("same");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Counts");
  c2->cd(5);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwoV3[0][1]->Divide(PtTwoV3[0][0]);
  PtTwoV3[0][1]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c2->cd(6);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwoV3[1][1]->Divide(PtTwoV3[1][0]);
  PtTwoV3[1][1]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c2->cd(7);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwoV3[2][1]->Divide(PtTwoV3[2][0]);
  PtTwoV3[2][1]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  c2->cd(8);
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");
  PtTwoV3[3][1]->Divide(PtTwoV3[3][0]);
  PtTwoV3[3][1]->Draw("ep");
  gPad->SetTitle(  ";p_{T} [GeV/#it{c}];Ratio");


}
