// c++ headers
#include <iostream>
#include <fstream>
#include <iomanip>

// root headers
#include <Rtypes.h>
#include <TMath.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;









void DrawSuppressionFactor(){
  

  TGraphErrors *SuppressionFactor2;
  TH1F* SuppressionFactor = new TH1F("SuppressionFactor", "SuppressionFactor", 10000, -0.5, 999.5);
  const Double_t ALICE_photo_W[9]    = { 19.130700000, 24.5643000, 31.541200000, 493.3870000, 633.5210000,   813.4570000,            97.1536,    124.748,      160.179      };
  // const Double_t ALICE_photo_sig[9]  = {  0.00882072,  0.0138693,  0.0168523,     0.0473069,    0.0523713,     0.0628193,             0.0218977,   0.0230963,    0.0246617  };
  // const Double_t ALICE_photo_sta[9]  = {  0.000260183, 0.000307948,0.000593475,   0.00775179,   0.00799204,    0.024678,              0.00101229,  0.00136407,   0.00506907 };
  const Double_t ALICE_photo_sig[9]  = {  0.00882072,  0.0138693,  0.0168523,     0.0473069,   0.0523713,      0.0628193,       0.0218977,  0.0246131,  0.0246617};
  const Double_t ALICE_photo_sta[9]  = {  0.000260183, 0.000307948,0.000593475,   0.00775179,   0.00799204,    0.024678,        0.00461081, 0.000587666,  0.00708194 };
  const Double_t ALICE_photo_sys[9]  = {  0.00069,    0.00056,    0.00135,       0.00695,   0.00899,           0.01334,         0.00175, 0.00141, 0.00485};
  const Double_t IA[9]               = {  1.3259E-02, 1.7995E-02, 2.2492E-02,     1.4200E-01,   1.6707E-01,    1.9655E-01,            4.9208E-02,  5.7978E-02,   6.8268E-02 };
  Double_t SuppressionFactors[9]     = { 0,0,0, 0,0,0, 0,0,0 };
  Double_t SuppressionFactors_sta[9] = { 0,0,0, 0,0,0, 0,0,0 };
  Double_t SuppressionFactors_sys[9] = { 0,0,0, 0,0,0, 0,0,0 };
  for (Int_t i = 0; i < 9; i++) {
    SuppressionFactors[i]     = TMath::Sqrt(ALICE_photo_sig[i]/IA[i]);
    SuppressionFactors_sta[i] = ALICE_photo_sta[i]/(2.*TMath::Sqrt(ALICE_photo_sig[i]*IA[i]));
    SuppressionFactors_sys[i] = ALICE_photo_sys[i]/(2.*TMath::Sqrt(ALICE_photo_sig[i]*IA[i]));
    cout << "Supp +/- Stat. +/- Sys. = " << SuppressionFactors[i] << " +/- " << SuppressionFactors_sta[i] << " +/- " << SuppressionFactors_sys[i] << endl;
  }
  SuppressionFactor2  = new TGraphErrors(9, ALICE_photo_W, SuppressionFactors, 0, SuppressionFactors_sta);
  TBox *Box[9];
  Double_t xBoxL[9]     = { 0,0,0, 0,0,0, 0,0,0 };
  Double_t xBoxR[9]     = { 0,0,0, 0,0,0, 0,0,0 };
  Double_t yBoxL[9]     = { 0,0,0, 0,0,0, 0,0,0 };
  Double_t yBoxR[9]     = { 0,0,0, 0,0,0, 0,0,0 };
  for (size_t i = 0; i < 9; i++) {
    xBoxL[i] = ALICE_photo_W[i]*0.95;
    xBoxR[i] = ALICE_photo_W[i]*1.05;
    yBoxL[i] = SuppressionFactors[i]-SuppressionFactors_sys[i];
    yBoxR[i] = SuppressionFactors[i]+SuppressionFactors_sys[i];
  }
  for (size_t i = 0; i < 9; i++) {
    Box[i] = new TBox(xBoxL[i],yBoxL[i],xBoxR[i],yBoxR[i]);
  }
  for (Int_t i = 0; i < 9; i++) {
    Int_t bin = SuppressionFactor->GetXaxis()->FindBin(ALICE_photo_W[i]);
    cout << "bin = " << bin << endl;
    SuppressionFactor->SetBinContent( SuppressionFactor->GetXaxis()->FindBin(ALICE_photo_W[i]),  SuppressionFactors[i]);
    SuppressionFactor->SetBinError(   SuppressionFactor->GetXaxis()->FindBin(ALICE_photo_W[i]),  SuppressionFactors_sta[i] );
  }
  new TCanvas;
  SuppressionFactor->SetMarkerStyle(20);
  SuppressionFactor->SetMarkerColor(2);
  SuppressionFactor2->SetLineColor(2);
  SuppressionFactor2->SetMaximum(1.5);
  SuppressionFactor2->SetMinimum(0.000000);
  SuppressionFactor2->GetXaxis()->SetLimits(10.,1000.);

    






  gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTopMargin(0.14);
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  // gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
  SuppressionFactor->SetTitle("");
  SuppressionFactor->GetXaxis()->SetTitleOffset(1.25);
  SuppressionFactor->GetYaxis()->SetTitleOffset(1.25);
  SuppressionFactor->GetXaxis()->SetTitleSize(0.045);
  SuppressionFactor->GetYaxis()->SetTitleSize(0.045);
  SuppressionFactor->GetXaxis()->SetLabelSize(0.045);
  SuppressionFactor->GetYaxis()->SetLabelSize(0.045);
  SuppressionFactor->GetXaxis()->SetTitleFont(42);
  SuppressionFactor->GetYaxis()->SetTitleFont(42);
  SuppressionFactor->GetXaxis()->SetLabelFont(42);
  SuppressionFactor->GetYaxis()->SetLabelFont(42);

  SuppressionFactor->GetXaxis()->SetTitle("W_{#gammaPb} [GeV]");
  SuppressionFactor->GetYaxis()->SetTitle("S_{Pb} [a.u.]");
  SuppressionFactor->GetYaxis()->SetRangeUser(0.,1.5);
  SuppressionFactor->GetXaxis()->SetRangeUser(10.,1000.);
  SuppressionFactor->SetLineWidth(2);
  SuppressionFactor->SetLineColor(2);
  SuppressionFactor->Draw("same");
  for (size_t i = 0; i < 9; i++) {
    Box[i]->SetFillStyle(0);
    Box[i]->SetLineWidth(2);
    Box[i]->SetLineColor(2);
    Box[i]->Draw("same");
  }

  Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  Double_t bxmin = TMath::Power((mjpsi/1000.),2.);
  Double_t bxmax = TMath::Power((mjpsi/10.),2.);
  TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  fbx->SetParameter(0, mjpsi);

  TGaxis *axis2 = new TGaxis(1000., 1.5, 10., 1.5, "fbx", 510, "+G");
  axis2->SetTextFont(42);
  axis2->SetLabelFont(42);
  Double_t siz = 0.045;
  axis2->SetTitleSize(siz); axis2->SetLabelSize(siz);
  axis2->SetLabelOffset(-0.035);
  axis2->SetTitleOffset();
  axis2->SetTitle("Bjorken-#it{x}");
  axis2->Draw("same");

  TLatex *bxtit2 = new TLatex();
  bxtit2->SetTextFont(42);
  bxtit2->SetTextSize(siz);
  bxtit2->SetTextAlign(31);
  bxtit2->DrawLatex(1000., 1.5+450, "Bjorken-#it{x}");
  TLatex* latex6 = new TLatex();
  latex6->SetTextSize(0.045);
  latex6->SetTextFont(42);
  latex6->SetTextAlign(11);
  latex6->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.7,"This thesis");
  // latex5->DrawLatex(0.31,0.8,"Coherent J/#psi");
  latex6->DrawLatex(0.4,0.78,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
//   latex6->DrawLatex(0.4,0.7,"This thesis");
  latex6->DrawLatex(0.4,0.62,"Coherent J/#psi");

  gPad->SaveAs("SuppressionFactor-refined.pdf", "recreate");

}