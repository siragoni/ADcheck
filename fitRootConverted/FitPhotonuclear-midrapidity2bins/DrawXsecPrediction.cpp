#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* -
 * - Original macro:
 * - https://root.cern/doc/v610/graphShade_8C.html
 */
void DrawXsecPrediction(Int_t selectionFlag = 0)
{
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);





  TGraphErrors *Coherent0N0N;
  TGraphErrors *Coherent0NXN;
  TGraphErrors *CoherentXN0N;
  TGraphErrors *CoherentXNXN;
  TGraphErrors *Coherent0N0N2;
  TGraphErrors *Coherent0NXN2;
  TGraphErrors *CoherentXN0N2;
  TGraphErrors *CoherentXNXN2;
  Double_t DSigmaDy0N0N[3]  = { 1.5860 * 0.5,2.3150 * 0.5,2.6570 * 0.5};
  Double_t DSigmaDy0NXN[3]  = { 0.0430,0.1092,0.2302};
  Double_t DSigmaDyXN0N[3]  = { 0.0993 * 2.0,0.1681 * 2.0,0.2359 * 2.0};
  Double_t DSigmaDyXNXN[3]  = { 0.0830,0.1693,0.2681};

  Double_t x1Error[3]       = {  0.25, 0.25, 0.25 };
  // Double_t x1[3]            = { -3.75+(4-2.5)/6, -3.25+(4-2.5)/6, -2.75+(4-2.5)/6};
  Double_t x1[3]            = { -3.75, -3.25, -2.75};

  Double_t y1Error0N0N[3]   = { 0.0489 * 0.5,0.0333 * 0.5,0.0760 * 0.5};
  Double_t y1Error0NXN[3]   = { 0.0035,0.0053,0.0126};
  Double_t y1ErrorXN0N[3]   = { 0.0078 * 2.0,0.0098 * 2.0,0.0208 * 2.0};
  Double_t y1ErrorXNXN[3]   = { 0.0119,0.0120,0.0208};

  Double_t y2Error0N0N[3]  = { 0.0844,0.0844,0.0844};
  Double_t y2Error0NXN[3]  = { 0.3329,0.1910,0.1167};
  Double_t y2ErrorXN0N[3]  = { 0.0857,0.0885,0.0917};
  Double_t y2ErrorXNXN[3]  = { 0.1412,0.1427,0.1401};


  TBox *Box0N0N2[3];
  TBox *Box0NXN2[3];
  TBox *BoxXN0N2[3];
  TBox *BoxXNXN2[3];
  Double_t xBoxL[3]     = { x1[0]-x1Error[0],x1[1]-x1Error[1],x1[2]-x1Error[2] };
  Double_t xBoxR[3]     = { x1[0]+x1Error[0],x1[1]+x1Error[1],x1[2]+x1Error[2] };
  Double_t yBox0N0NL[3] = { DSigmaDy0N0N[0]*(1.-y2Error0N0N[0]), DSigmaDy0N0N[1]*(1.-y2Error0N0N[1]), DSigmaDy0N0N[2]*(1.-y2Error0N0N[2]) };
  Double_t yBox0N0NR[3] = { DSigmaDy0N0N[0]*(1.+y2Error0N0N[0]), DSigmaDy0N0N[1]*(1.+y2Error0N0N[1]), DSigmaDy0N0N[2]*(1.+y2Error0N0N[2]) };
  Double_t yBox0NXNL[3] = { DSigmaDy0NXN[0]*(1.-y2Error0NXN[0]), DSigmaDy0NXN[1]*(1.-y2Error0NXN[1]), DSigmaDy0NXN[2]*(1.-y2Error0NXN[2]) };
  Double_t yBox0NXNR[3] = { DSigmaDy0NXN[0]*(1.+y2Error0NXN[0]), DSigmaDy0NXN[1]*(1.+y2Error0NXN[1]), DSigmaDy0NXN[2]*(1.+y2Error0NXN[2]) };
  Double_t yBoxXN0NL[3] = { DSigmaDyXN0N[0]*(1.-y2ErrorXN0N[0]), DSigmaDyXN0N[1]*(1.-y2ErrorXN0N[1]), DSigmaDyXN0N[2]*(1.-y2ErrorXN0N[2]) };
  Double_t yBoxXN0NR[3] = { DSigmaDyXN0N[0]*(1.+y2ErrorXN0N[0]), DSigmaDyXN0N[1]*(1.+y2ErrorXN0N[1]), DSigmaDyXN0N[2]*(1.+y2ErrorXN0N[2]) };
  Double_t yBoxXNXNL[3] = { DSigmaDyXNXN[0]*(1.-y2ErrorXNXN[0]), DSigmaDyXNXN[1]*(1.-y2ErrorXNXN[1]), DSigmaDyXNXN[2]*(1.-y2ErrorXNXN[2]) };
  Double_t yBoxXNXNR[3] = { DSigmaDyXNXN[0]*(1.+y2ErrorXNXN[0]), DSigmaDyXNXN[1]*(1.+y2ErrorXNXN[1]), DSigmaDyXNXN[2]*(1.+y2ErrorXNXN[2]) };

  for (size_t i = 0; i < 3; i++) {
    // Box0N0N2[i] = new TBox(xBoxL[i]+0.1,yBox0N0NL[i],xBoxR[i]-0.1,yBox0N0NR[i]);
    // Box0NXN2[i] = new TBox(xBoxL[i]+0.1,yBox0NXNL[i],xBoxR[i]-0.1,yBox0NXNR[i]);
    // BoxXN0N2[i] = new TBox(xBoxL[i]+0.1,yBoxXN0NL[i],xBoxR[i]-0.1,yBoxXN0NR[i]);
    // BoxXNXN2[i] = new TBox(xBoxL[i]+0.1,yBoxXNXNL[i],xBoxR[i]-0.1,yBoxXNXNR[i]);
    Box0N0N2[i] = new TBox(xBoxL[i]+0.0,yBox0N0NL[i],xBoxR[i]-0.0,yBox0N0NR[i]);
    Box0NXN2[i] = new TBox(xBoxL[i]+0.0,yBox0NXNL[i],xBoxR[i]-0.0,yBox0NXNR[i]);
    BoxXN0N2[i] = new TBox(xBoxL[i]+0.0,yBoxXN0NL[i],xBoxR[i]-0.0,yBoxXN0NR[i]);
    BoxXNXN2[i] = new TBox(xBoxL[i]+0.0,yBoxXNXNL[i],xBoxR[i]-0.0,yBoxXNXNR[i]);
  }

  // Coherent0N0N  = new TGraphErrors(3, x1, DSigmaDy0N0N, x1Error, y1Error0N0N);
  // Coherent0NXN  = new TGraphErrors(3, x1, DSigmaDy0NXN, x1Error, y1Error0NXN);
  // CoherentXN0N  = new TGraphErrors(3, x1, DSigmaDyXN0N, x1Error, y1ErrorXN0N);
  // CoherentXNXN  = new TGraphErrors(3, x1, DSigmaDyXNXN, x1Error, y1ErrorXNXN);
  Coherent0N0N  = new TGraphErrors(3, x1, DSigmaDy0N0N, 0, y1Error0N0N);
  Coherent0NXN  = new TGraphErrors(3, x1, DSigmaDy0NXN, 0, y1Error0NXN);
  CoherentXN0N  = new TGraphErrors(3, x1, DSigmaDyXN0N, 0, y1ErrorXN0N);
  CoherentXNXN  = new TGraphErrors(3, x1, DSigmaDyXNXN, 0, y1ErrorXNXN);

  Coherent0N0N2 = new TGraphErrors(3, x1, DSigmaDy0N0N, x1Error, y2Error0N0N);
  Coherent0NXN2 = new TGraphErrors(3, x1, DSigmaDy0NXN, x1Error, y2Error0NXN);
  CoherentXN0N2 = new TGraphErrors(3, x1, DSigmaDyXN0N, x1Error, y2ErrorXN0N);
  CoherentXNXN2 = new TGraphErrors(3, x1, DSigmaDyXNXN, x1Error, y2ErrorXNXN);

  TAxis *axis = Coherent0N0N->GetXaxis();
  axis->SetLimits(-5.,1.);                 // along X
  Coherent0N0N->GetHistogram()->SetMaximum(2.);   // along
  Coherent0N0N->GetHistogram()->SetMinimum(0.);  //   Y
  TAxis *axis2 = Coherent0NXN->GetXaxis();
  axis2->SetLimits(-5.,1.);                 // along X
  Coherent0NXN->GetHistogram()->SetMaximum(2.);   // along
  Coherent0NXN->GetHistogram()->SetMinimum(0.);  //   Y
  TAxis *axis3 = CoherentXN0N->GetXaxis();
  axis3->SetLimits(-5.,1.);                 // along X
  CoherentXN0N->GetHistogram()->SetMaximum(2.);   // along
  CoherentXN0N->GetHistogram()->SetMinimum(0.);  //   Y
  TAxis *axis4 = CoherentXNXN->GetXaxis();
  axis4->SetLimits(-5.,1.);                 // along X
  CoherentXNXN->GetHistogram()->SetMaximum(2.);   // along
  CoherentXNXN->GetHistogram()->SetMinimum(0.);  //   Y

  TMultiGraph *mg = new TMultiGraph();
  Coherent0N0N->SetMarkerStyle(20);
  Coherent0N0N->SetMarkerColor(2);
  Coherent0N0N->SetLineColor(2);
  Coherent0N0N2->SetMarkerStyle(25);
  Coherent0N0N2->SetMarkerColor(2);
  Coherent0N0N2->SetLineColor(2);
  if ( selectionFlag == 0 || selectionFlag == 1 ){
    mg->Add(Coherent0N0N);
  }
  // mg->Add(Coherent0N0N2);
  Coherent0NXN->SetMarkerStyle(20);
  Coherent0NXN->SetMarkerColor(3);
  Coherent0NXN->SetLineColor(3);
  Coherent0NXN2->SetMarkerStyle(25);
  Coherent0NXN2->SetMarkerColor(3);
  Coherent0NXN2->SetLineColor(3);
  if ( selectionFlag == 2 ){
    mg->Add(Coherent0NXN);
  }
  // mg->Add(Coherent0NXN2);
  CoherentXN0N->SetMarkerStyle(20);
  CoherentXN0N->SetMarkerColor(4);
  CoherentXN0N->SetLineColor(4);
  CoherentXN0N2->SetMarkerStyle(25);
  CoherentXN0N2->SetMarkerColor(4);
  CoherentXN0N2->SetLineColor(4);
  if ( selectionFlag == 0 || selectionFlag == 3 ){
    mg->Add(CoherentXN0N);
  }
  // mg->Add(CoherentXN0N2);
  CoherentXNXN->SetMarkerStyle(20);
  CoherentXNXN->SetMarkerColor(6);
  CoherentXNXN->SetLineColor(6);
  CoherentXNXN2->SetMarkerStyle(25);
  CoherentXNXN2->SetMarkerColor(6);
  CoherentXNXN2->SetLineColor(6);
  if ( selectionFlag == 0 || selectionFlag == 4 ){
    mg->Add(CoherentXNXN);
  }
  // mg->Add(CoherentXNXN2);

  Coherent0N0N->SetTitle("J/#Psi 0N0N");
  Coherent0NXN->SetTitle("J/#Psi 0NXN");
  CoherentXN0N->SetTitle("J/#Psi XN0N");
  CoherentXNXN->SetTitle("J/#Psi XNXN");

  mg->SetMaximum(4.0);
  mg->SetMinimum(0.0000001);
  mg->GetXaxis()->SetLimits(-5,1.2);
  // mg->GetYaxis()->SetLimits(0.,1.);
  mg->GetHistogram()->SetMaximum(3.);   // along
  mg->GetHistogram()->SetMinimum(0.);  //   Y
  mg->Draw("APsame");
  // mg->Draw("same");
  // TAxis *axis = Coherent0N0N->GetXaxis();
  //
  // axis->SetLimits(-5.,1.);                 // along X
  // Coherent0N0N->GetHistogram()->SetMaximum(1.);   // along
  // Coherent0N0N->GetHistogram()->SetMinimum(0.);  //   Y

  // Coherent0N0N->Draw("APsame");
  // Coherent0NXN->Draw("APsame");
  // CoherentXN0N->Draw("APsame");
  // CoherentXNXN->Draw("APsame");


  // for (size_t i = 0; i < 3; i++) {
  //   Box0N0N2[i]->SetFillStyle(0);
  //   Box0NXN2[i]->SetFillStyle(0);
  //   BoxXN0N2[i]->SetFillStyle(0);
  //   BoxXNXN2[i]->SetFillStyle(0);
  //
  //   Box0N0N2[i]->SetLineWidth(2);
  //   Box0NXN2[i]->SetLineWidth(2);
  //   BoxXN0N2[i]->SetLineWidth(2);
  //   BoxXNXN2[i]->SetLineWidth(2);
  //
  //   Box0N0N2[i]->SetLineColor(2);
  //   Box0NXN2[i]->SetLineColor(3);
  //   BoxXN0N2[i]->SetLineColor(4);
  //   BoxXNXN2[i]->SetLineColor(6);
  //   Box0N0N2[i]->Draw("same");
  //   Box0NXN2[i]->Draw("same");
  //   BoxXN0N2[i]->Draw("same");
  //   BoxXNXN2[i]->Draw("same");
  //
  // }

  // mg->Draw("APL");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
  // Change the axis limits
  // gPad->BuildLegend();












  // TFile* file = new TFile("Michal-Broz-xsec/xSection_Cent.root");
  TFile* file = new TFile("../../Michal-Broz-xsec/xSection_Cent2.root");
  TGraphErrors* Cent0N0N = (TGraphErrors*) file->Get("gXsection_Cent_Stat_0n0n");
  TGraphAsymmErrors* Cent0N0Nerr = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0n0n");
  for (int i=0;i<Cent0N0N->GetN();i++) Cent0N0N->GetY()[i] *= 0.5;
  for (int i=0;i<Cent0N0N->GetN();i++) Cent0N0N->SetPointError(i,0.,Cent0N0N->GetErrorY(i) * 0.5);
  for (int i=0;i<Cent0N0Nerr->GetN();i++) Cent0N0Nerr->SetPointError(i,Cent0N0Nerr->GetErrorX(i) * 0.5,Cent0N0Nerr->GetErrorX(i) * 0.5, Cent0N0Nerr->GetErrorY(i) * 0.25,Cent0N0Nerr->GetErrorY(i) * 0.25);
  for (int i=0;i<Cent0N0Nerr->GetN();i++) Cent0N0Nerr->GetY()[i] *= 0.5;
  TGraphErrors* Cent0NXN = (TGraphErrors*) file->Get("gXsection_Cent_Stat_0nXn");
  TGraphAsymmErrors* Cent0NXNerr = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0nXn");
  for (int i=0;i<Cent0NXN->GetN();i++) Cent0NXN->GetY()[i] *= 1.;
  for (int i=0;i<Cent0NXN->GetN();i++) Cent0NXN->SetPointError(i,0.,Cent0NXN->GetErrorY(i) * 1.);
  for (int i=0;i<Cent0NXNerr->GetN();i++) Cent0NXNerr->SetPointError(i,Cent0NXNerr->GetErrorX(i) * 1.,Cent0NXNerr->GetErrorX(i) * 1., Cent0NXNerr->GetErrorY(i) * 1., Cent0NXNerr->GetErrorY(i) * 1.);
  for (int i=0;i<Cent0NXNerr->GetN();i++) Cent0NXNerr->GetY()[i] *= 1.;
  // for (int i=0;i<Cent0NXN->GetN();i++) Cent0NXN->GetY()[i] *= 0.5;
  // for (int i=0;i<Cent0NXN->GetN();i++) Cent0NXN->SetPointError(i,0.,Cent0NXN->GetErrorY(i) * 0.5);
  // for (int i=0;i<Cent0NXNerr->GetN();i++) Cent0NXNerr->SetPointError(i,Cent0NXNerr->GetErrorX(i) * 0.5,Cent0NXNerr->GetErrorX(i) * 0.5, Cent0NXNerr->GetErrorY(i) * 0.5, Cent0NXNerr->GetErrorY(i) * 0.5);
  // for (int i=0;i<Cent0NXNerr->GetN();i++) Cent0NXNerr->GetY()[i] *= 0.5;

  // TGraphAsymmErrors* CentXN0N = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_Xn0n");
  TGraphErrors* CentXNXN = (TGraphErrors*) file->Get("gXsection_Cent_Stat_XnXn");
  TGraphAsymmErrors* CentXNXNerr = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_XnXn");
  Cent0N0N->SetMarkerStyle(20);
  Cent0N0N->SetMarkerColor(2);
  Cent0N0N->SetLineColor(2);
  if ( selectionFlag == 0 || selectionFlag == 1 ){
    mg->Add(Cent0N0N);
  }
  // mg->Add(Cent0N0N);
  Cent0NXN->SetMarkerStyle(20);
  Cent0NXN->SetMarkerColor(4);
  Cent0NXN->SetLineColor(4);
  if ( selectionFlag == 0 || selectionFlag == 2 || selectionFlag == 3 ){
    mg->Add(Cent0NXN);
  }
  // CentXN0N->SetMarkerStyle(20);
  // CentXN0N->SetMarkerColor(4);
  // CentXN0N->SetLineColor(4);
  // mg->Add(CentXN0N);
  CentXNXN->SetMarkerStyle(20);
  CentXNXN->SetMarkerColor(6);
  CentXNXN->SetLineColor(6);
  if ( selectionFlag == 0 || selectionFlag == 4 ){
    mg->Add(CentXNXN);
  }
  // mg->Add(CentXNXN);
  mg->Draw("APsame");










  TFile* file2 = new TFile("../../Michal-Broz-xsec/predictions_with_neutrons.root");

  // LTA onon
  TGraph* lta_onon_h = (TGraph*) file2->Get("glta00_h");
  TGraph* lta_onon_l = (TGraph*) file2->Get("glta00_l");
  for (int i=0;i<lta_onon_h->GetN();i++) lta_onon_h->GetY()[i] *= 0.5;
  for (int i=0;i<lta_onon_l->GetN();i++) lta_onon_l->GetY()[i] *= 0.5;
  TGraph* lta_onon   = new TGraph(lta_onon_l->GetN() +lta_onon_h->GetN());

  for(Int_t i = 0;                  i < lta_onon_h->GetN();                       i++) {
    lta_onon->SetPoint(i, lta_onon_h->GetX()[i], lta_onon_h->GetY()[i]);
    cout << "[" << i << "] = " << lta_onon_h->GetX()[i] << endl;

  }
  cout << "stop" << endl;
  for(Int_t i = lta_onon_h->GetN(); i < (lta_onon_l->GetN() +lta_onon_h->GetN()); i++) {
    // lta_onon->SetPoint(i, lta_onon_l->GetX()[i], lta_onon_l->GetY()[i]);
    lta_onon->SetPoint(i, lta_onon_l->GetX()[lta_onon_l->GetN() +lta_onon_h->GetN()-i], lta_onon_l->GetY()[lta_onon_l->GetN() +lta_onon_h->GetN()-i]);
    cout << "[" << i << "] = " << lta_onon_l->GetX()[i-lta_onon_h->GetN()] << endl;

  }

  lta_onon->SetFillColorAlpha(kRed-2, 0.6);
  //lta_onon->SetFillColor(kOrange-2);
  lta_onon->SetLineColor(kRed-2);
  lta_onon->SetLineWidth(1);

  lta_onon->SetFillStyle(3008);

  lta_onon_h->SetLineColor(kRed-2);
  lta_onon_h->SetLineWidth(1);
  lta_onon_l->SetLineColor(kRed-2);
  lta_onon_l->SetLineWidth(1);
  static TGraph *lta_onon_coll[3] = {lta_onon, lta_onon_l, lta_onon_h};

  if ( selectionFlag == 0 || selectionFlag == 1 ){
    lta_onon_coll[0]->Draw("fsame");
    lta_onon_coll[1]->Draw("lsame");
    lta_onon_coll[2]->Draw("lsame");
  }

  // lta_onon_coll[0]->Draw("fsame");
  // lta_onon_coll[1]->Draw("lsame");
  // lta_onon_coll[2]->Draw("lsame");





  // LTA onxn
  TGraph* lta_onxn_h = (TGraph*) file2->Get("glta0X_h");
  TGraph* lta_onxn_l = (TGraph*) file2->Get("glta0X_l");
  for (int i=0;i<lta_onxn_h->GetN();i++) lta_onxn_h->GetY()[i] *= 1.;
  for (int i=0;i<lta_onxn_l->GetN();i++) lta_onxn_l->GetY()[i] *= 1.;
  // for (int i=0;i<lta_onxn_h->GetN();i++) lta_onxn_h->GetY()[i] *= 0.5;
  // for (int i=0;i<lta_onxn_l->GetN();i++) lta_onxn_l->GetY()[i] *= 0.5;
  TGraph* lta_onxn   = new TGraph(lta_onxn_l->GetN() +lta_onxn_h->GetN());

  for(Int_t i = 0;                  i < lta_onxn_h->GetN();                       i++) {
    lta_onxn->SetPoint(i, lta_onxn_h->GetX()[i], lta_onxn_h->GetY()[i]);
    // cout << "[" << i << "] = " << lta_onxn_h->GetX()[i] << endl;

  }
  // cout << "stop" << endl;
  for(Int_t i = lta_onxn_h->GetN(); i < (lta_onxn_l->GetN() +lta_onxn_h->GetN()); i++) {
    // lta_onon->SetPoint(i, lta_onon_l->GetX()[i], lta_onon_l->GetY()[i]);
    lta_onxn->SetPoint(i, lta_onxn_l->GetX()[lta_onxn_l->GetN() +lta_onxn_h->GetN()-i], lta_onxn_l->GetY()[lta_onxn_l->GetN() +lta_onxn_h->GetN()-i]);
    // cout << "[" << i << "] = " << lta_onon_l->GetX()[i-lta_onon_h->GetN()] << endl;

  }

  lta_onxn->SetFillColorAlpha(kBlue-2, 0.6);
  //lta_onon->SetFillColor(kOrange-2);
  lta_onxn->SetLineColor(kBlue-2);
  lta_onxn->SetLineWidth(1);

  lta_onxn->SetFillStyle(3008);

  lta_onxn_h->SetLineColor(kBlue-2);
  lta_onxn_h->SetLineWidth(1);
  lta_onxn_l->SetLineColor(kBlue-2);
  lta_onxn_l->SetLineWidth(1);
  static TGraph *lta_onxn_coll[3] = {lta_onxn, lta_onxn_l, lta_onxn_h};

  if ( selectionFlag == 0 || selectionFlag == 2 || selectionFlag == 3 ){
    lta_onxn_coll[0]->Draw("fsame");
    lta_onxn_coll[1]->Draw("lsame");
    lta_onxn_coll[2]->Draw("lsame");
  }

  // lta_onxn_coll[0]->Draw("fsame");
  // lta_onxn_coll[1]->Draw("lsame");
  // lta_onxn_coll[2]->Draw("lsame");









  // LTA xnxn
  TGraph* lta_xnxn_h = (TGraph*) file2->Get("gltaXX_h");
  TGraph* lta_xnxn_l = (TGraph*) file2->Get("gltaXX_l");
  for (int i=0;i<lta_xnxn_h->GetN();i++) lta_xnxn_h->GetY()[i] *= 1;
  for (int i=0;i<lta_xnxn_l->GetN();i++) lta_xnxn_l->GetY()[i] *= 1;
  TGraph* lta_xnxn   = new TGraph(lta_xnxn_l->GetN() +lta_xnxn_h->GetN());

  for(Int_t i = 0;                  i < lta_xnxn_h->GetN();                       i++) {
    lta_xnxn->SetPoint(i, lta_xnxn_h->GetX()[i], lta_xnxn_h->GetY()[i]);
    // cout << "[" << i << "] = " << lta_onxn_h->GetX()[i] << endl;

  }
  // cout << "stop" << endl;
  for(Int_t i = lta_xnxn_h->GetN(); i < (lta_xnxn_l->GetN() +lta_xnxn_h->GetN()); i++) {
    // lta_onon->SetPoint(i, lta_onon_l->GetX()[i], lta_onon_l->GetY()[i]);
    lta_xnxn->SetPoint(i, lta_xnxn_l->GetX()[lta_xnxn_l->GetN() +lta_xnxn_h->GetN()-i], lta_xnxn_l->GetY()[lta_xnxn_l->GetN() +lta_xnxn_h->GetN()-i]);
    // cout << "[" << i << "] = " << lta_onon_l->GetX()[i-lta_onon_h->GetN()] << endl;

  }

  lta_xnxn->SetFillColorAlpha(kMagenta-2, 0.6);
  //lta_onon->SetFillColor(kOrange-2);
  lta_xnxn->SetLineColor(kMagenta-2);
  lta_xnxn->SetLineWidth(1);

  lta_xnxn->SetFillStyle(3008);

  lta_xnxn_h->SetLineColor(kMagenta-2);
  lta_xnxn_h->SetLineWidth(1);
  lta_xnxn_l->SetLineColor(kMagenta-2);
  lta_xnxn_l->SetLineWidth(1);
  static TGraph *lta_xnxn_coll[3] = {lta_xnxn, lta_xnxn_l, lta_xnxn_h};

  if ( selectionFlag == 0 || selectionFlag == 4 ){
    lta_xnxn_coll[0]->Draw("fsame");
    lta_xnxn_coll[1]->Draw("lsame");
    lta_xnxn_coll[2]->Draw("lsame");
  }

  // lta_xnxn_coll[0]->Draw("fsame");
  // lta_xnxn_coll[1]->Draw("lsame");
  // lta_xnxn_coll[2]->Draw("lsame");


  gPad->Modified();
  gPad->Update();












  TBox *Box0N0Nmid[3];
  TBox *Box0NXNmid[3];
  TBox *BoxXNXNmid[3];
  Double_t xBoxLmid[3]     = { -0.8,-0.2,0.2 };
  Double_t xBoxRmid[3]     = { -0.2,0.2,0.8 };
  Double_t yBox0N0NLmid[3] = {0,0,0};
  Double_t yBox0N0NRmid[3] = {0,0,0};
  Double_t yBox0NXNLmid[3] = {0,0,0};
  Double_t yBox0NXNRmid[3] = {0,0,0};
  Double_t yBoxXNXNLmid[3] = {0,0,0};
  Double_t yBoxXNXNRmid[3] = {0,0,0};
  for (Int_t i = 0; i < 3; i++)
  {
    yBox0N0NLmid[i] = (Cent0N0N->GetY()[i]-Cent0N0Nerr->GetErrorY(i));
    yBox0N0NRmid[i] = (Cent0N0N->GetY()[i]+Cent0N0Nerr->GetErrorY(i));
    yBox0NXNLmid[i] = (Cent0NXN->GetY()[i]-Cent0NXNerr->GetErrorY(i));
    yBox0NXNRmid[i] = (Cent0NXN->GetY()[i]+Cent0NXNerr->GetErrorY(i));
    yBoxXNXNLmid[i] = (CentXNXN->GetY()[i]-CentXNXNerr->GetErrorY(i));
    yBoxXNXNRmid[i] = (CentXNXN->GetY()[i]+CentXNXNerr->GetErrorY(i));
  }



  for (size_t i = 0; i < 3; i++) {
    Box0N0Nmid[i] = new TBox(xBoxLmid[i]+0.0,yBox0N0NLmid[i],xBoxRmid[i]-0.0,yBox0N0NRmid[i]);
    Box0NXNmid[i] = new TBox(xBoxLmid[i]+0.0,yBox0NXNLmid[i],xBoxRmid[i]-0.0,yBox0NXNRmid[i]);
    BoxXNXNmid[i] = new TBox(xBoxLmid[i]+0.0,yBoxXNXNLmid[i],xBoxRmid[i]-0.0,yBoxXNXNRmid[i]);
  }

















  for (size_t i = 0; i < 3; i++) {
    Box0N0N2[i]->SetFillStyle(0);
    Box0NXN2[i]->SetFillStyle(0);
    BoxXN0N2[i]->SetFillStyle(0);
    BoxXNXN2[i]->SetFillStyle(0);
    Box0N0Nmid[i]->SetFillStyle(0);
    Box0NXNmid[i]->SetFillStyle(0);
    BoxXNXNmid[i]->SetFillStyle(0);

    Box0N0N2[i]->SetLineWidth(2);
    Box0NXN2[i]->SetLineWidth(2);
    BoxXN0N2[i]->SetLineWidth(2);
    BoxXNXN2[i]->SetLineWidth(2);
    Box0N0Nmid[i]->SetLineWidth(2);
    Box0NXNmid[i]->SetLineWidth(2);
    BoxXNXNmid[i]->SetLineWidth(2);

    Box0N0N2[i]->SetLineColor(2);
    Box0NXN2[i]->SetLineColor(3);
    BoxXN0N2[i]->SetLineColor(4);
    BoxXNXN2[i]->SetLineColor(6);
    Box0N0Nmid[i]->SetLineColor(2);
    Box0NXNmid[i]->SetLineColor(4);
    BoxXNXNmid[i]->SetLineColor(6);
    // Box0N0N2[i]->Draw("same");
    // Box0NXN2[i]->Draw("same");
    // BoxXN0N2[i]->Draw("same");
    // BoxXNXN2[i]->Draw("same");
    if ( selectionFlag == 0 || selectionFlag == 1 ){
      Box0N0N2[i]->Draw("same");
      Box0N0Nmid[i]->Draw("same");
    }
    if (  selectionFlag == 2 ){
      Box0NXN2[i]->Draw("same");
      Box0NXNmid[i]->Draw("same");
    }
    if ( selectionFlag == 0 || selectionFlag == 3 ){
      BoxXN0N2[i]->Draw("same");
      Box0NXNmid[i]->Draw("same");
    }
    if ( selectionFlag == 0 || selectionFlag == 4 ){
      BoxXNXN2[i]->Draw("same");
      BoxXNXNmid[i]->Draw("same");
    }


  }

  // TLegend *leg_pt = new TLegend(0.5,0.45,0.85,0.79);
  TLegend *leg_pt = new TLegend(0.2,0.55,0.5,0.89);
  leg_pt->SetFillStyle(0);
  leg_pt->SetBorderSize(0);
  leg_pt->SetTextSize(0.03);
  leg_pt->AddEntry(Coherent0N0N,"Coh J/#psi 0N0N x0.5", "EP");
  // leg_pt->AddEntry(Coherent0NXN,"Coh J/#psi 0NXN", "EP");
  leg_pt->AddEntry(CoherentXN0N,"Coh J/#psi XN0N", "EP");
  leg_pt->AddEntry(CoherentXNXN,"Coh J/#psi XNXN", "EP");
  leg_pt->Draw();

  TLegend *leg_theo = new TLegend(0.5,0.55,0.85,0.89);
  leg_theo->SetFillStyle(0);
  leg_theo->SetBorderSize(0);
  leg_theo->SetTextSize(0.03);
  leg_theo->AddEntry(lta_onon_coll[0],"LTA 0N0N x0.5", "F");
  leg_theo->AddEntry(lta_onxn_coll[0],"LTA 0NXN", "F");
  leg_theo->AddEntry(lta_xnxn_coll[0],"LTA XNXN", "F");
  leg_theo->Draw();



  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{s_{NN}} = 5.02 TeV");


  gPad->Modified();
  gPad->Update();
  // gPad->SaveAs("pngResults/xsec.pdf", "recreate");
  gPad->SaveAs(Form("xsec_%d.pdf", selectionFlag), "recreate");


}
