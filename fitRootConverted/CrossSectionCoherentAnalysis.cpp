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
void CrossSectionCoherentAnalysis(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *CoherentMine;
  TGraphErrors *CoherentTomas;
  TGraphErrors *CoherentApril;
  TGraphErrors *CoherentMicha;
  TGraphErrors *CoherentLHCb;
  Double_t DSigmaDyCMUP6[6]   = { 0,0,0,0,0,0 };
  Double_t DSigmaDyTomas[6]   = { 0,0,0,0,0,0 };
  Double_t DSigmaDyLHCb[5]    = { 1.10,1.73,2.28,2.60,3.0 };
  Double_t DSigmaDyCentralB[5]= { 3.5, 3.8, 3.75, 3.8, 3.5 };
  Double_t DSigmaDyCMUP11[6]  = { 1.615,1.938,2.377,2.831,3.018,3.531 };
  // Double_t fICMUP6[6]         = { 0.154, 0.121, 0.094, 0.074, 0.073, 0.078 }; // Different valeus for high pt model
  Double_t fICMUP6[6]         = { 0.083, 0.070, 0.058, 0.065, 0.051, 0.044 }; // Different valeus for high pt model
  // Double_t fICMUP6[6]         = { 0.064, 0.058, 0.060, 0.052, 0.049, 0.049 }; // from coherent paper

  Double_t errfICMUP6[6]      = { 0.018,0.008,0.005,0.012,0.011,0.010 };

  Double_t fD                 = 0.055;
  Double_t eJPsiEvgeny[6]     = { 0.051, 0.140, 0.204, 0.191, 0.119, 0.029 };
  // Double_t eJPsiCMUPC[6]      = { 0.155, 0.145, 0.209, 0.193, 0.120, 0.029 };
  // Double_t fICMUP11[6]        = { 0.155, 0.116, 0.096, 0.075, 0.075, 0.084 };
  Double_t eJPsiCMUP6[6]      = { 0.053, 0.145, 0.209, 0.193, 0.120, 0.029 };
  Double_t LUMI               = 532.926;
  Double_t BR                 = 0.05961;
  Double_t NumOfJPsi[6]       = { 684.0, 2302.0, 4324.0, 4670.0, 3095.0, 817.0 };

  for( Int_t iLoop = 0; iLoop < 6; iLoop++ ){
      DSigmaDyCMUP6[iLoop] = NumOfJPsi[iLoop]/( (1+fICMUP6[iLoop]+fD)*eJPsiCMUP6[iLoop]*BR*LUMI*0.25*0.95*1000 );
  }

  Double_t x1[6]           = { -4.0+0*(4.0-2.5)/6.0+0.125, -4.0+1*(4.0-2.5)/6.0+0.125, -4.0+2*(4.0-2.5)/6.0+0.125, -4.0+3*(4.0-2.5)/6.0+0.125, -4.0+4*(4.0-2.5)/6.0+0.125,-4.0+5*(4.0-2.5)/6.0+0.125    };
  Double_t y1ErrorCMUP6[6] = { TMath::Sqrt( ( 38.0*38.0/( 684.0* 684.0)) + (errfICMUP6[0]*errfICMUP6[0]/((1+fICMUP6[0]+fD)*(1+fICMUP6[0]+fD))) )*DSigmaDyCMUP6[0],
                               TMath::Sqrt( ( 71.0*71.0/(2302.0*2302.0)) + (errfICMUP6[1]*errfICMUP6[1]/((1+fICMUP6[1]+fD)*(1+fICMUP6[1]+fD))) )*DSigmaDyCMUP6[1],
                               TMath::Sqrt( ( 98.0*98.0/(4324.0*4324.0)) + (errfICMUP6[2]*errfICMUP6[2]/((1+fICMUP6[2]+fD)*(1+fICMUP6[2]+fD))) )*DSigmaDyCMUP6[2],
                               TMath::Sqrt( ( 83.0*83.0/(4670.0*4670.0)) + (errfICMUP6[3]*errfICMUP6[3]/((1+fICMUP6[3]+fD)*(1+fICMUP6[3]+fD))) )*DSigmaDyCMUP6[3],
                               TMath::Sqrt( ( 84.0*84.0/(3095.0*3095.0)) + (errfICMUP6[4]*errfICMUP6[4]/((1+fICMUP6[4]+fD)*(1+fICMUP6[4]+fD))) )*DSigmaDyCMUP6[4],
                               TMath::Sqrt( ( 44.0*44.0/( 817.0* 817.0)) + (errfICMUP6[5]*errfICMUP6[5]/((1+fICMUP6[5]+fD)*(1+fICMUP6[5]+fD))) )*DSigmaDyCMUP6[5]
                               };
  Double_t y1Error[6]       = { 0.060, 0.042, 0.040, 0.047, 0.061, 0.137 };
  Double_t y1ErrorCentral[5]= { 0.1, 0.1, 0.1, 0.1, 0.1 };
  Double_t y1ErrorLHCb[5]   = { 0.22, 0.15, 0.15, 0.19, 0.4 };
  // Double_t x1Error[6]     = {  (4-2.5)/6, (4-2.5)/6, (4-2.5)/6 };
  Double_t x1Error[6]       = {  0,0,0,0,0,0 };
  Double_t x1ErrorCentral[5]= {  0,0,0,0,0 };
  Double_t x1Central[5]     = { -0.55, -0.25, 0.0, 0.25,  0.55    };
  Double_t x1LHCb[5]        = { -4.25, -3.75, -3.25, -2.75, -2.25   };



  CoherentMine  = new TGraphErrors(6, x1,        DSigmaDyCMUP6,    x1Error,          y1ErrorCMUP6       );
  CoherentApril = new TGraphErrors(6, x1,        DSigmaDyCMUP11,   x1Error,          y1Error            );
  CoherentMicha = new TGraphErrors(5, x1Central, DSigmaDyCentralB, x1ErrorCentral,   y1ErrorCentral     );
  CoherentLHCb  = new TGraphErrors(5, x1LHCb,    DSigmaDyLHCb,     x1ErrorCentral,   y1ErrorLHCb       );

  TMultiGraph *mg = new TMultiGraph();
  CoherentMine->SetMarkerStyle(20);
  CoherentMine->SetMarkerColor(2);
  CoherentMine->SetLineColor(2);
  mg->Add(CoherentMine);
  CoherentApril->SetMarkerStyle(20);
  CoherentApril->SetMarkerColor(3);
  CoherentApril->SetLineColor(3);
  mg->Add(CoherentApril);
  CoherentMicha->SetMarkerStyle(20);
  CoherentMicha->SetMarkerColor(4);
  CoherentMicha->SetLineColor(4);
  mg->Add(CoherentMicha);
  CoherentLHCb->SetMarkerStyle(20);
  CoherentLHCb->SetMarkerColor(6);
  CoherentLHCb->SetLineColor(6);
  mg->Add(CoherentLHCb);

  CoherentMine ->SetTitle("J/#Psi CMUP6");
  CoherentApril->SetTitle("J/#Psi CMUP11");
  CoherentMicha->SetTitle("J/#Psi Central Barrel");
  CoherentLHCb ->SetTitle("J/#Psi LHCb");

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

  cout << "CMUP6_0 = " << endl << DSigmaDyCMUP6[0] << " +/- " << y1ErrorCMUP6[0] << endl;
  cout << "CMUP6_1 = " << endl << DSigmaDyCMUP6[1] << " +/- " << y1ErrorCMUP6[1] << endl;
  cout << "CMUP6_2 = " << endl << DSigmaDyCMUP6[2] << " +/- " << y1ErrorCMUP6[2] << endl;
  cout << "CMUP6_3 = " << endl << DSigmaDyCMUP6[3] << " +/- " << y1ErrorCMUP6[3] << endl;
  cout << "CMUP6_4 = " << endl << DSigmaDyCMUP6[4] << " +/- " << y1ErrorCMUP6[4] << endl;
  cout << "CMUP6_5 = " << endl << DSigmaDyCMUP6[5] << " +/- " << y1ErrorCMUP6[5] << endl;

}
