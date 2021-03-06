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
/* - Drawing function for the invariant mass distributions' retrieved valued,
   - after the fit in terms of CosTheta bins.
   -
 */
void RatioSigmaWithWithoutAD(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *Coherent0N0N;
  TGraphErrors *Coherent0NXN;
  TGraphErrors *CoherentXN0N;
  TGraphErrors *CoherentXNXN;
  Double_t DSigmaDy0N0N[3]  = { 0,0,0};
  Double_t DSigmaDy0NXN[3]  = { 0,0,0};
  Double_t DSigmaDyXN0N[3]  = { 0,0,0};
  Double_t DSigmaDyXNXN[3]  = { 0,0,0};
  Double_t fI[3]            = { 0.0595, 0.056,  0.041 };
  Double_t fI0N0N[3]        = { 0.0430, 0.027,  0.043 };
  Double_t fI0NXN[3]        = { 0.0840, 0.084,  0.084 };
  Double_t fIXN0N[3]        = { 5.6430, 2.335,  1.511 };
  Double_t fIXNXN[3]        = { 0.5000, 0.485,  0.500 };
  Double_t fD               = 0.055;
  Double_t eJPsi[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI             = 532.926;
  Double_t BR               = 0.05961;
  Double_t NumOfJPsi0N0N[3] = { 2355, 6951, 2869 };
  Double_t NumOfJPsi0NXN[3] = {  180,  590,  303 };
  Double_t NumOfJPsiXN0N[3] = {  298,  839,  417 };
  Double_t NumOfJPsiXNXN[3] = {   91,  357,  204 };

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ){
      // DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI0N0N[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI0NXN[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fIXN0N[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fIXNXN[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
  }

  Double_t x1[3]          = { -4+1*(4-2.5)/6, -4+3*(4-2.5)/6, -4+5*(4-2.5)/6    };
  Double_t y1Error0N0N[3] = { // 83 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 119/( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 98 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              72 /( (1+fI0N0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              119/( (1+fI0N0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              98 /( (1+fI0N0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1Error0NXN[3] = { // 18 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 34 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 32 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              18 /( (1+fI0NXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              34 /( (1+fI0NXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              32 /( (1+fI0NXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXN0N[3] = { // 29 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 52 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 29 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              29 /( (1+fIXN0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              52 /( (1+fIXN0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              29 /( (1+fIXN0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXNXN[3] = { // 16 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 25 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 24 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              16 /( (1+fIXNXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              25 /( (1+fIXNXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              24 /( (1+fIXNXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t x1Error[3]     = {  (4-2.5)/6, (4-2.5)/6, (4-2.5)/6 };




  Double_t DSigmaDy0N0N_withoutAD[3]  = { 0,0,0};
  Double_t DSigmaDy0NXN_withoutAD[3]  = { 0,0,0};
  Double_t DSigmaDyXN0N_withoutAD[3]  = { 0,0,0};
  Double_t DSigmaDyXNXN_withoutAD[3]  = { 0,0,0};
  Double_t fI_withoutAD[3]            = { 0.0595, 0.056,  0.041 };
  Double_t fI0N0N_withoutAD[3]        = { 0.0430, 0.027,  0.043 };
  Double_t fI0NXN_withoutAD[3]        = { 0.0840, 0.084,  0.084 };
  Double_t fIXN0N_withoutAD[3]        = { 5.6430, 2.335,  1.511 };
  Double_t fIXNXN_withoutAD[3]        = { 0.5000, 0.485,  0.500 };
  Double_t fD_withoutAD               = 0.055;
  Double_t eJPsi_withoutAD[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI_withoutAD             = 532.926;
  Double_t BR_withoutAD               = 0.05961;
  Double_t NumOfJPsi0N0N_withoutAD[3] = { 2361, 6988, 2872 };
  Double_t NumOfJPsi0NXN_withoutAD[3] = {  185,  630,  315 };
  Double_t NumOfJPsiXN0N_withoutAD[3] = {  327,  937,  468 };
  Double_t NumOfJPsiXNXN_withoutAD[3] = {  109,  444,  242 };

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ){
      // DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDy0N0N_withoutAD[iLoop] = NumOfJPsi0N0N_withoutAD[iLoop]/( (1+fI0N0N_withoutAD[iLoop]+fD_withoutAD)*eJPsi_withoutAD[iLoop]*BR_withoutAD*LUMI_withoutAD*0.5*0.95*1000 );
      DSigmaDy0NXN_withoutAD[iLoop] = NumOfJPsi0NXN_withoutAD[iLoop]/( (1+fI0NXN_withoutAD[iLoop]+fD_withoutAD)*eJPsi_withoutAD[iLoop]*BR_withoutAD*LUMI_withoutAD*0.5*0.95*1000 );
      DSigmaDyXN0N_withoutAD[iLoop] = NumOfJPsiXN0N_withoutAD[iLoop]/( (1+fIXN0N_withoutAD[iLoop]+fD_withoutAD)*eJPsi_withoutAD[iLoop]*BR_withoutAD*LUMI_withoutAD*0.5*0.95*1000 );
      DSigmaDyXNXN_withoutAD[iLoop] = NumOfJPsiXNXN_withoutAD[iLoop]/( (1+fIXNXN_withoutAD[iLoop]+fD_withoutAD)*eJPsi_withoutAD[iLoop]*BR_withoutAD*LUMI_withoutAD*0.5*0.95*1000 );
  }

  Double_t x1_withoutAD[3]          = { -4+1*(4-2.5)/6, -4+3*(4-2.5)/6, -4+5*(4-2.5)/6    };
  Double_t y1Error0N0N_withoutAD[3] = { // 83 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 119/( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 98 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              72 /( (1+fI0N0N_withoutAD[0]+fD_withoutAD)*eJPsi_withoutAD[0]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              119/( (1+fI0N0N_withoutAD[1]+fD_withoutAD)*eJPsi_withoutAD[1]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              98 /( (1+fI0N0N_withoutAD[2]+fD_withoutAD)*eJPsi_withoutAD[2]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 )
                              };
  Double_t y1Error0NXN_withoutAD[3] = { // 18 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 34 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 32 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              18 /( (1+fI0NXN_withoutAD[0]+fD_withoutAD)*eJPsi_withoutAD[0]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              34 /( (1+fI0NXN_withoutAD[1]+fD_withoutAD)*eJPsi_withoutAD[1]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              32 /( (1+fI0NXN_withoutAD[2]+fD_withoutAD)*eJPsi_withoutAD[2]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXN0N_withoutAD[3] = { // 29 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 52 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 29 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              29 /( (1+fIXN0N_withoutAD[0]+fD_withoutAD)*eJPsi_withoutAD[0]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              52 /( (1+fIXN0N_withoutAD[1]+fD_withoutAD)*eJPsi_withoutAD[1]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              29 /( (1+fIXN0N_withoutAD[2]+fD_withoutAD)*eJPsi_withoutAD[2]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXNXN_withoutAD[3] = { // 16 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 25 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 24 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              16 /( (1+fIXNXN_withoutAD[0]+fD_withoutAD)*eJPsi_withoutAD[0]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              25 /( (1+fIXNXN_withoutAD[1]+fD_withoutAD)*eJPsi_withoutAD[1]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 ),
                              24 /( (1+fIXNXN_withoutAD[2]+fD_withoutAD)*eJPsi_withoutAD[2]*BR_withoutAD*LUMI_withoutAD*1.5*0.95*1000 )
                              };
  Double_t x1Error_withoutAD[3]     = {  (4-2.5)/6, (4-2.5)/6, (4-2.5)/6 };




  Coherent0N0N = new TGraphErrors(3, x1, DSigmaDy0N0N, x1Error, y1Error0N0N);
  Coherent0NXN = new TGraphErrors(3, x1, DSigmaDy0NXN, x1Error, y1Error0NXN);
  CoherentXN0N = new TGraphErrors(3, x1, DSigmaDyXN0N, x1Error, y1ErrorXN0N);
  CoherentXNXN = new TGraphErrors(3, x1, DSigmaDyXNXN, x1Error, y1ErrorXNXN);

  TMultiGraph *mg = new TMultiGraph();
  Coherent0N0N->SetMarkerStyle(20);
  Coherent0N0N->SetMarkerColor(2);
  Coherent0N0N->SetLineColor(2);
  mg->Add(Coherent0N0N);
  Coherent0NXN->SetMarkerStyle(20);
  Coherent0NXN->SetMarkerColor(3);
  Coherent0NXN->SetLineColor(3);
  mg->Add(Coherent0NXN);
  CoherentXN0N->SetMarkerStyle(20);
  CoherentXN0N->SetMarkerColor(4);
  CoherentXN0N->SetLineColor(4);
  mg->Add(CoherentXN0N);
  CoherentXNXN->SetMarkerStyle(20);
  CoherentXNXN->SetMarkerColor(6);
  CoherentXNXN->SetLineColor(6);
  mg->Add(CoherentXNXN);

  Coherent0N0N->SetTitle("J/#Psi 0N0N");
  Coherent0NXN->SetTitle("J/#Psi 0NXN");
  CoherentXN0N->SetTitle("J/#Psi XN0N");
  CoherentXNXN->SetTitle("J/#Psi XNXN");

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
  // Change the axis limits
  gPad->BuildLegend();
  gPad->Modified();
  // mg->GetXaxis()->SetLimits(-1., 1.);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(6000.);
  // c2->Print("pngResults/MultiGraph1Dview.png");
  // mg->Draw("a fb l3d");
  // c2->Print("pngResults/MultiGraph2Dview.png");
  // return c2;

  cout << "0N0N : " << endl << DSigmaDy0N0N[0] << endl << DSigmaDy0N0N[1] << endl << DSigmaDy0N0N[2] << endl;
  cout << "0NXN : " << endl << DSigmaDy0NXN[0] << endl << DSigmaDy0NXN[1] << endl << DSigmaDy0NXN[2] << endl;
  cout << "XN0N : " << endl << DSigmaDyXN0N[0] << endl << DSigmaDyXN0N[1] << endl << DSigmaDyXN0N[2] << endl;
  cout << "XNXN : " << endl << DSigmaDyXNXN[0] << endl << DSigmaDyXNXN[1] << endl << DSigmaDyXNXN[2] << endl;

  cout << "0N0N E: " << endl << y1Error0N0N[0] << endl << y1Error0N0N[1] << endl << y1Error0N0N[2] << endl;
  cout << "0NXN E: " << endl << y1Error0NXN[0] << endl << y1Error0NXN[1] << endl << y1Error0NXN[2] << endl;
  cout << "XN0N E: " << endl << y1ErrorXN0N[0] << endl << y1ErrorXN0N[1] << endl << y1ErrorXN0N[2] << endl;
  cout << "XNXN E: " << endl << y1ErrorXNXN[0] << endl << y1ErrorXNXN[1] << endl << y1ErrorXNXN[2] << endl;



  TCanvas *Ratios = new TCanvas("Ratios","Ratios",600,400);
  Double_t Ratio0N0N[3]  = { 0,0,0 };
  Double_t Error0N0N[3]  = { 0,0,0 };
  Double_t Ratio0NXN[3]  = { 0,0,0 };
  Double_t Error0NXN[3]  = { 0,0,0 };
  Double_t RatioXN0N[3]  = { 0,0,0 };
  Double_t ErrorXN0N[3]  = { 0,0,0 };
  Double_t RatioXNXN[3]  = { 0,0,0 };
  Double_t ErrorXNXN[3]  = { 0,0,0 };
  for (size_t i = 0; i < 3; i++) {
    Ratio0N0N[i] = DSigmaDy0N0N[i] / DSigmaDy0N0N_withoutAD[i];
    Error0N0N[i] = y1Error0N0N[i]  + y1Error0N0N_withoutAD[i];
    Ratio0NXN[i] = DSigmaDy0NXN[i] / DSigmaDy0NXN_withoutAD[i];
    Error0NXN[i] = y1Error0NXN[i]  + y1Error0NXN_withoutAD[i];
    RatioXN0N[i] = DSigmaDyXN0N[i] / DSigmaDyXN0N_withoutAD[i];
    ErrorXN0N[i] = y1ErrorXN0N[i]  + y1ErrorXN0N_withoutAD[i];
    RatioXNXN[i] = DSigmaDyXNXN[i] / DSigmaDyXNXN_withoutAD[i];
    ErrorXNXN[i] = y1ErrorXNXN[i]  + y1ErrorXNXN_withoutAD[i];
  }
  // TGraphErrors* RatioPlot = new TGraphErrors(3, x1, Ratio, x1Error, Error);
  // RatioPlot->Draw("APL");
  // RatioPlot->GetXaxis()->SetTitle("y");
  // RatioPlot->GetYaxis()->SetTitle("Ratio [au]");
  // // Change the axis limits
  // gPad->BuildLegend();
  // gPad->Modified();

  TGraphErrors *Ratios0N0N;
  TGraphErrors *Ratios0NXN;
  TGraphErrors *RatiosXN0N;
  TGraphErrors *RatiosXNXN;
  Ratios0N0N = new TGraphErrors(3, x1, Ratio0N0N, x1Error, Error0N0N);
  Ratios0NXN = new TGraphErrors(3, x1, Ratio0NXN, x1Error, Error0NXN);
  RatiosXN0N = new TGraphErrors(3, x1, RatioXN0N, x1Error, ErrorXN0N);
  RatiosXNXN = new TGraphErrors(3, x1, RatioXNXN, x1Error, ErrorXNXN);

  TMultiGraph *mg2 = new TMultiGraph();
  Ratios0N0N->SetMarkerStyle(20);
  Ratios0N0N->SetMarkerColor(2);
  Ratios0N0N->SetLineColor(2);
  mg2->Add(Ratios0N0N);
  Ratios0NXN->SetMarkerStyle(20);
  Ratios0NXN->SetMarkerColor(3);
  Ratios0NXN->SetLineColor(3);
  mg2->Add(Ratios0NXN);
  RatiosXN0N->SetMarkerStyle(20);
  RatiosXN0N->SetMarkerColor(4);
  RatiosXN0N->SetLineColor(4);
  mg2->Add(RatiosXN0N);
  RatiosXNXN->SetMarkerStyle(20);
  RatiosXNXN->SetMarkerColor(6);
  RatiosXNXN->SetLineColor(6);
  mg2->Add(RatiosXNXN);

  Ratios0N0N->SetTitle("J/#Psi 0N0N");
  Ratios0NXN->SetTitle("J/#Psi 0NXN");
  RatiosXN0N->SetTitle("J/#Psi XN0N");
  RatiosXNXN->SetTitle("J/#Psi XNXN");

  mg2->Draw("APL");
  mg2->GetXaxis()->SetTitle("y");
  mg2->GetYaxis()->SetTitle("Ratios [a.u.]");
  // Change the axis limits
  gPad->BuildLegend();
  gPad->Modified();


}
