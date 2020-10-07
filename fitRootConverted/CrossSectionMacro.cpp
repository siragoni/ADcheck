#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TString.h"
using namespace std;
#include <math.h>
#include <vector>


//_____________________________________________________________________________
/* -
   -
 */
void CrossSection(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *Coherent0N0N;
  TGraphErrors *Coherent0NXN;
  TGraphErrors *CoherentXN0N;
  TGraphErrors *CoherentXNXN;
  TGraphErrors *Coherent0N0N_;
  TGraphErrors *Coherent0NXN_;
  TGraphErrors *CoherentXN0N_;
  TGraphErrors *CoherentXNXN_;


  // Double_t DSigmaDy0N0N[3]  = { 0,0,0};
  // Double_t DSigmaDy0NXN[3]  = { 0,0,0};
  // Double_t DSigmaDyXN0N[3]  = { 0,0,0};
  // Double_t DSigmaDyXNXN[3]  = { 0,0,0};
  Double_t DSigmaDy0N0N[3]  = { 1.644916975, 2.310923447, 2.702658714 };
  Double_t DSigmaDy0NXN[3]  = { 0.057679491, 0.116341016, 0.160207017 };
  Double_t DSigmaDyXN0N[3]  = { 0.075702188, 0.133863829, 0.210089027 };
  Double_t DSigmaDyXNXN[3]  = { 0.048983822, 0.095930944, 0.169544318 };

  Double_t ey0N0N[3] = { 0.039134952, 0.033247312, 0.060306836 };
  Double_t ey0NXN[3] = { 0.008066679, 0.009370019, 0.018885340 };
  Double_t eyXN0N[3] = { 0.011973722, 0.053616776, 0.016704030 };
  Double_t eyXNXN[3] = { 0.006641558, 0.006633638, 0.013797879 };

  Double_t ex[3]  = { 0., 0., 0. };
  Double_t ex_[3] = { 0.25, 0.25, 0.25 };


  Double_t ey0N0N_[3] = { DSigmaDy0N0N[0]*0.1, DSigmaDy0N0N[1]*0.1, DSigmaDy0N0N[2]*0.1 };
  Double_t ey0NXN_[3] = { DSigmaDy0NXN[0]*0.1, DSigmaDy0NXN[1]*0.1, DSigmaDy0NXN[2]*0.1 };
  Double_t eyXN0N_[3] = { DSigmaDyXN0N[0]*0.1, DSigmaDyXN0N[1]*0.1, DSigmaDyXN0N[2]*0.1 };
  Double_t eyXNXN_[3] = { DSigmaDyXNXN[0]*0.1, DSigmaDyXNXN[1]*0.1, DSigmaDyXNXN[2]*0.1 };


  Double_t x1[3]          = { -4+1*(4-2.5)/6, -4+3*(4-2.5)/6, -4+5*(4-2.5)/6    };


  Coherent0N0N  = new TGraphErrors(3, x1, DSigmaDy0N0N, ex, ey0N0N);
  Coherent0NXN  = new TGraphErrors(3, x1, DSigmaDy0NXN, ex, ey0NXN);
  CoherentXN0N  = new TGraphErrors(3, x1, DSigmaDyXN0N, ex, eyXN0N);
  CoherentXNXN  = new TGraphErrors(3, x1, DSigmaDyXNXN, ex, eyXNXN);
  Coherent0N0N_ = new TGraphErrors(3, x1, DSigmaDy0N0N, ex_, ey0N0N_);
  Coherent0NXN_ = new TGraphErrors(3, x1, DSigmaDy0NXN, ex_, ey0NXN_);
  CoherentXN0N_ = new TGraphErrors(3, x1, DSigmaDyXN0N, ex_, eyXN0N_);
  CoherentXNXN_ = new TGraphErrors(3, x1, DSigmaDyXNXN, ex_, eyXNXN_);






  Coherent0N0N->SetMarkerStyle(20);
  Coherent0N0N->SetMarkerSize(2);
  Coherent0N0N->SetLineWidth(3);
  Coherent0N0N->SetMarkerColor(2);
  Coherent0N0N->SetLineColor(2);
  Coherent0NXN->SetMarkerStyle(20);
  Coherent0NXN->SetMarkerSize(2);
  Coherent0NXN->SetLineWidth(3);
  Coherent0NXN->SetMarkerColor(3);
  Coherent0NXN->SetLineColor(3);
  CoherentXN0N->SetMarkerStyle(20);
  CoherentXN0N->SetMarkerSize(2);
  CoherentXN0N->SetLineWidth(3);
  CoherentXN0N->SetMarkerColor(4);
  CoherentXN0N->SetLineColor(4);
  CoherentXNXN->SetMarkerStyle(20);
  CoherentXNXN->SetMarkerColor(6);
  CoherentXNXN->SetLineColor(6);
  CoherentXNXN->SetMarkerSize(2);
  CoherentXNXN->SetLineWidth(3);




  Coherent0N0N_->SetMarkerStyle(20);
  Coherent0N0N_->SetMarkerSize(2);
  Coherent0N0N_->SetLineWidth(3);
  Coherent0N0N_->SetMarkerColor(2);
  Coherent0N0N_->SetLineColor(2);
  Coherent0NXN_->SetMarkerStyle(20);
  Coherent0NXN_->SetMarkerSize(2);
  Coherent0NXN_->SetLineWidth(3);
  Coherent0NXN_->SetMarkerColor(3);
  Coherent0NXN_->SetLineColor(3);
  CoherentXN0N_->SetMarkerStyle(20);
  CoherentXN0N_->SetMarkerSize(2);
  CoherentXN0N_->SetLineWidth(3);
  CoherentXN0N_->SetMarkerColor(4);
  CoherentXN0N_->SetLineColor(4);
  CoherentXNXN_->SetMarkerStyle(20);
  CoherentXNXN_->SetMarkerColor(6);
  CoherentXNXN_->SetLineColor(6);
  CoherentXNXN_->SetMarkerSize(2);
  CoherentXNXN_->SetLineWidth(3);

  Coherent0N0N_->SetFillColor(0);
  Coherent0N0N_->SetFillStyle(0);
  Coherent0NXN_->SetFillColor(0);
  Coherent0NXN_->SetFillStyle(0);
  CoherentXN0N_->SetFillColor(0);
  CoherentXN0N_->SetFillStyle(0);
  CoherentXNXN_->SetFillColor(0);
  CoherentXNXN_->SetFillStyle(0);



  TCanvas *c = new TCanvas("c","c",800,800);


  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gPad->SetTitle(  ";y;d#sigma/dy [mb]" );

  Coherent0N0N->GetYaxis()->SetRangeUser(0., 6.);
  Coherent0N0N->GetXaxis()->SetRangeUser(-5.0, -1.5);
  Coherent0NXN->GetYaxis()->SetRangeUser(0., 6.);
  Coherent0NXN->GetXaxis()->SetRangeUser(-5.0, -1.5);
  CoherentXN0N->GetYaxis()->SetRangeUser(0., 6.);
  CoherentXN0N->GetXaxis()->SetRangeUser(-5.0, -1.5);
  CoherentXNXN->GetYaxis()->SetRangeUser(0., 6.);
  CoherentXNXN->GetXaxis()->SetRangeUser(-5.0, -1.5);


  TMultiGraph *mg = new TMultiGraph();
  mg->Add(Coherent0N0N, "AP");
  mg->Add(Coherent0NXN, "AP");
  mg->Add(CoherentXN0N, "AP");
  mg->Add(CoherentXNXN, "AP");
  mg->Add(Coherent0N0N_, "2");
  mg->Add(Coherent0NXN_, "2");
  mg->Add(CoherentXN0N_, "2");
  mg->Add(CoherentXNXN_, "2");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
  mg->GetYaxis()->SetRangeUser(0., 9.);
  mg->GetXaxis()->SetRangeUser(-5.0, -1.5);
  mg->Draw("APL");


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset");




  TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
  leg_pt->SetFillStyle(0);
  leg_pt->SetBorderSize(0);
  leg_pt->SetTextSize(0.04);
  leg_pt->AddEntry(Coherent0N0N,"0N0N", "PI");
  leg_pt->AddEntry(Coherent0NXN,"0NXN", "PI");
  leg_pt->AddEntry(CoherentXN0N,"XN0N", "PI");
  leg_pt->AddEntry(CoherentXNXN,"XNXN", "PI");
  leg_pt->Draw();


  // Coherent0N0N->Draw("AP");
  // // Coherent0NXN->Draw("APsame");
  // // CoherentXN0N->Draw("APsame");
  // // CoherentXNXN->Draw("APsame");
  // Coherent0NXN->Draw("AP");
  // CoherentXN0N->Draw("AP");
  // CoherentXNXN->Draw("AP");
  // //
  // Coherent0N0N_->Draw("2");
  // // Coherent0NXN_->Draw("2same");
  // // CoherentXN0N_->Draw("2same");
  // // CoherentXNXN_->Draw("2same");
  // // Coherent0NXN_->Draw("2");
  // // CoherentXN0N_->Draw("2");
  // // CoherentXNXN_->Draw("2");
  //
  // // Coherent0N0N->Draw("AP");
  // // Coherent0NXN->Draw("APsame");
  // // CoherentXN0N->Draw("APsame");
  // // CoherentXNXN->Draw("APsame");




  // Coherent0N0N->Draw("AP");
  // Coherent0N0N_->Draw("2");
  //
  //
  // // Coherent0NXN->Draw("APsame");
  // // CoherentXN0N->Draw("APsame");
  // // CoherentXNXN->Draw("APsame");
  // Coherent0NXN->Draw("APsame");
  // Coherent0NXN_->Draw("2same");

  // CoherentXN0N->Draw("AP");
  // CoherentXN0N_->Draw("2");
  //
  // CoherentXNXN->Draw("AP");
  // CoherentXNXN_->Draw("2");

  // //
  // Coherent0N0N_->Draw("2");
  // Coherent0NXN_->Draw("2same");
  // CoherentXN0N_->Draw("2same");
  // CoherentXNXN_->Draw("2same");
  // Coherent0NXN_->Draw("2");
  // CoherentXN0N_->Draw("2");
  // CoherentXNXN_->Draw("2");





}
