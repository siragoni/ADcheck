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


//_____________________________________________________________________________
/* -
 * - Ratios of AD multiplicities.
 * -
 */
void RatiosADmultiplicities(const char* AnalysisName){
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("ADcheck_", listings);



  TH1F *fADAmultiplicityTotalH          = (TH1F*)listings->FindObject("fADAmultiplicityTotalH");
  TH1F *fADAmultiplicity0N0NclassTotalH = (TH1F*)listings->FindObject("fADAmultiplicity0N0NclassTotalH");
  TH1F *fADAmultiplicity0NXNclassTotalH = (TH1F*)listings->FindObject("fADAmultiplicity0NXNclassTotalH");
  TH1F *fADAmultiplicityXN0NclassTotalH = (TH1F*)listings->FindObject("fADAmultiplicityXN0NclassTotalH");
  TH1F *fADAmultiplicityXNXNclassTotalH = (TH1F*)listings->FindObject("fADAmultiplicityXNXNclassTotalH");

  TH1F *fADCmultiplicityTotalH          = (TH1F*)listings->FindObject("fADCmultiplicityTotalH");
  TH1F *fADCmultiplicity0N0NclassTotalH = (TH1F*)listings->FindObject("fADCmultiplicity0N0NclassTotalH");
  TH1F *fADCmultiplicity0NXNclassTotalH = (TH1F*)listings->FindObject("fADCmultiplicity0NXNclassTotalH");
  TH1F *fADCmultiplicityXN0NclassTotalH = (TH1F*)listings->FindObject("fADCmultiplicityXN0NclassTotalH");
  TH1F *fADCmultiplicityXNXNclassTotalH = (TH1F*)listings->FindObject("fADCmultiplicityXNXNclassTotalH");


  // fADAmultiplicityTotalH         ->Rebin(10);
  // fADAmultiplicity0N0NclassTotalH->Rebin(10);
  // fADAmultiplicity0NXNclassTotalH->Rebin(10);
  // fADAmultiplicityXN0NclassTotalH->Rebin(10);
  // fADAmultiplicityXNXNclassTotalH->Rebin(10);
  //
  // fADCmultiplicityTotalH         ->Rebin(10);
  // fADCmultiplicity0N0NclassTotalH->Rebin(10);
  // fADCmultiplicity0NXNclassTotalH->Rebin(10);
  // fADCmultiplicityXN0NclassTotalH->Rebin(10);
  // fADCmultiplicityXNXNclassTotalH->Rebin(10);



  // fADAmultiplicityTotalH         ->Rebin(20);
  // fADAmultiplicity0N0NclassTotalH->Rebin(20);
  // fADAmultiplicity0NXNclassTotalH->Rebin(20);
  // fADAmultiplicityXN0NclassTotalH->Rebin(20);
  // fADAmultiplicityXNXNclassTotalH->Rebin(20);
  //
  // fADCmultiplicityTotalH         ->Rebin(20);
  // fADCmultiplicity0N0NclassTotalH->Rebin(20);
  // fADCmultiplicity0NXNclassTotalH->Rebin(20);
  // fADCmultiplicityXN0NclassTotalH->Rebin(20);
  // fADCmultiplicityXNXNclassTotalH->Rebin(20);


  fADAmultiplicityTotalH         ->Rebin(40);
  fADAmultiplicity0N0NclassTotalH->Rebin(40);
  fADAmultiplicity0NXNclassTotalH->Rebin(40);
  fADAmultiplicityXN0NclassTotalH->Rebin(40);
  fADAmultiplicityXNXNclassTotalH->Rebin(40);

  fADCmultiplicityTotalH         ->Rebin(40);
  fADCmultiplicity0N0NclassTotalH->Rebin(40);
  fADCmultiplicity0NXNclassTotalH->Rebin(40);
  fADCmultiplicityXN0NclassTotalH->Rebin(40);
  fADCmultiplicityXNXNclassTotalH->Rebin(40);


  fADAmultiplicityTotalH         ->Sumw2();
  fADAmultiplicity0N0NclassTotalH->Sumw2();
  fADAmultiplicity0NXNclassTotalH->Sumw2();
  fADAmultiplicityXN0NclassTotalH->Sumw2();
  fADAmultiplicityXNXNclassTotalH->Sumw2();

  fADCmultiplicityTotalH         ->Sumw2();
  fADCmultiplicity0N0NclassTotalH->Sumw2();
  fADCmultiplicity0NXNclassTotalH->Sumw2();
  fADCmultiplicityXN0NclassTotalH->Sumw2();
  fADCmultiplicityXNXNclassTotalH->Sumw2();


  TH1F* fADAmultiplicityTotalClone          = (TH1F*) fADAmultiplicityTotalH         ->Clone("fADAmultiplicityTotalClone");
  TH1F* fADAmultiplicity0N0NclassTotalClone = (TH1F*) fADAmultiplicity0N0NclassTotalH->Clone("fADAmultiplicity0N0NclassTotalClone");
  TH1F* fADAmultiplicity0NXNclassTotalClone = (TH1F*) fADAmultiplicity0NXNclassTotalH->Clone("fADAmultiplicity0NXNclassTotalClone");
  TH1F* fADAmultiplicityXN0NclassTotalClone = (TH1F*) fADAmultiplicityXN0NclassTotalH->Clone("fADAmultiplicityXN0NclassTotalClone");
  TH1F* fADAmultiplicityXNXNclassTotalClone = (TH1F*) fADAmultiplicityXNXNclassTotalH->Clone("fADAmultiplicityXNXNclassTotalClone");

  TH1F* fADCmultiplicityTotalClone          = (TH1F*) fADCmultiplicityTotalH         ->Clone("fADCmultiplicityTotalClone");
  TH1F* fADCmultiplicity0N0NclassTotalClone = (TH1F*) fADCmultiplicity0N0NclassTotalH->Clone("fADCmultiplicity0N0NclassTotalClone");
  TH1F* fADCmultiplicity0NXNclassTotalClone = (TH1F*) fADCmultiplicity0NXNclassTotalH->Clone("fADCmultiplicity0NXNclassTotalClone");
  TH1F* fADCmultiplicityXN0NclassTotalClone = (TH1F*) fADCmultiplicityXN0NclassTotalH->Clone("fADCmultiplicityXN0NclassTotalClone");
  TH1F* fADCmultiplicityXNXNclassTotalClone = (TH1F*) fADCmultiplicityXNXNclassTotalH->Clone("fADCmultiplicityXNXNclassTotalClone");


  TH1F* fADA_0NXNvsXN0N = (TH1F*) fADAmultiplicity0NXNclassTotalH->Clone("fADA_0NXNvsXN0N");
  fADA_0NXNvsXN0N->Divide(fADAmultiplicityXN0NclassTotalClone);

  TH1F* fADC_0NXNvsXN0N = (TH1F*) fADCmultiplicity0NXNclassTotalH->Clone("fADC_0NXNvsXN0N");
  fADC_0NXNvsXN0N->Divide(fADCmultiplicityXN0NclassTotalClone);

  TH1F* fADA_0NXNvs_ADC_0NXN = (TH1F*) fADAmultiplicity0NXNclassTotalH->Clone("fADA_0NXNvs_ADC_0NXN");
  fADA_0NXNvs_ADC_0NXN->Divide(fADCmultiplicity0NXNclassTotalClone);

  TH1F* fADA_XN0Nvs_ADC_XN0N = (TH1F*) fADAmultiplicityXN0NclassTotalH->Clone("fADA_XN0Nvs_ADC_XN0N");
  fADA_XN0Nvs_ADC_XN0N->Divide(fADCmultiplicityXN0NclassTotalClone);

  TH1F* fADA_0N0Nvs_ADC_XNXN = (TH1F*) fADAmultiplicity0N0NclassTotalH->Clone("fADA_0N0Nvs_ADC_XNXN");
  fADA_0N0Nvs_ADC_XNXN->Divide(fADCmultiplicityXNXNclassTotalClone);

  TH1F* fADA_XNXNvs_ADC_0N0N = (TH1F*) fADAmultiplicityXNXNclassTotalH->Clone("fADA_XNXNvs_ADC_0N0N");
  fADA_XNXNvs_ADC_0N0N->Divide(fADCmultiplicity0N0NclassTotalClone);




  TCanvas* PtDistrCanvas = new TCanvas( "RatiosADmulti", "RatiosADmulti", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  /* - Beautifying is starting now.
     -
   */
  fADA_0NXNvsXN0N->GetXaxis()->SetTitleOffset(1.25);
  // fADA_0NXNvsXN0N->GetYaxis()->SetTitleOffset(1.25);
  fADA_0NXNvsXN0N->GetYaxis()->SetTitleOffset(1.45);
  fADA_0NXNvsXN0N->GetXaxis()->SetTitleSize(0.045);
  fADA_0NXNvsXN0N->GetYaxis()->SetTitleSize(0.045);
  fADA_0NXNvsXN0N->GetXaxis()->SetLabelSize(0.045);
  fADA_0NXNvsXN0N->GetYaxis()->SetLabelSize(0.045);
  fADA_0NXNvsXN0N->GetXaxis()->SetTitleFont(42);
  fADA_0NXNvsXN0N->GetYaxis()->SetTitleFont(42);
  fADA_0NXNvsXN0N->GetXaxis()->SetLabelFont(42);
  fADA_0NXNvsXN0N->GetYaxis()->SetLabelFont(42);
  // fADA_0NXNvsXN0N->GetXaxis()->SetNdivisions(408);
  // fADA_0NXNvsXN0N->GetYaxis()->SetRangeUser(10, fDimuonPtDistributionDataH->GetMaximum()*10.);
  fADA_0NXNvsXN0N->GetYaxis()->SetRangeUser(fADA_0NXNvsXN0N->GetMaximum()*0.0001, fADA_0NXNvsXN0N->GetMaximum()*100.);
  fADA_0NXNvsXN0N->GetXaxis()->SetRangeUser(0, 200);
  // fADA_0NXNvsXN0N->GetXaxis()->SetRangeUser(0, 3);
  // gPad ->SetLogy();
  fADA_0NXNvsXN0N     ->Draw("PEsame");
  fADC_0NXNvsXN0N     -> SetLineColor(kRed);
  fADA_0NXNvs_ADC_0NXN-> SetLineColor(kMagenta);
  fADA_XN0Nvs_ADC_XN0N-> SetLineColor(kYellow+1);
  fADA_0N0Nvs_ADC_XNXN-> SetLineColor(kCyan);
  fADA_XNXNvs_ADC_0N0N-> SetLineColor(kYellow);
  // fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  // fTwoGammaToMuMedium -> SetLineColor(kGreen);
  // fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  // fHighPtTail         -> SetLineColor(kGreen+1);
  fADA_0NXNvsXN0N     -> SetLineWidth(3);
  fADC_0NXNvsXN0N     -> SetLineWidth(3);
  fADA_0NXNvs_ADC_0NXN-> SetLineWidth(3);
  fADA_XN0Nvs_ADC_XN0N-> SetLineWidth(3);
  fADA_0N0Nvs_ADC_XNXN-> SetLineWidth(3);
  fADA_XNXNvs_ADC_0N0N-> SetLineWidth(3);
  // fADC_0NXNvsXN0N     -> Draw("same");
  // fADA_0NXNvs_ADC_0NXN-> Draw("same");
  // fADA_XN0Nvs_ADC_XN0N-> Draw("same");
  // fADA_0N0Nvs_ADC_XNXN-> Draw("same");
  fADA_XNXNvs_ADC_0N0N-> Draw("same");


  // TLatex* latex = new TLatex();
  // latex->SetTextSize(0.05);
  // latex->SetTextFont(42);
  // latex->SetTextAlign(11);
  // latex->SetNDC();
  // // latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  // latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  // latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));

  TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  l->AddEntry(  fADA_0NXNvsXN0N,      "fADA_0NXNvsXN0N");
  l->AddEntry(  fADC_0NXNvsXN0N,      "fADC_0NXNvsXN0N");
  l->AddEntry(  fADA_0NXNvs_ADC_0NXN, "fADA_0NXNvs_ADC_0NXN");
  l->AddEntry(  fADA_XN0Nvs_ADC_XN0N, "fADA_XN0Nvs_ADC_XN0N");
  l->AddEntry(  fADA_0N0Nvs_ADC_XNXN, "fADA_0N0Nvs_ADC_XNXN");
  l->AddEntry(  fADA_XNXNvs_ADC_0N0N, "fADA_XNXNvs_ADC_0N0N");
  l->Draw();

  gPad->SaveAs("pngResults/RatiosAD.png", "RECREATE");




}
