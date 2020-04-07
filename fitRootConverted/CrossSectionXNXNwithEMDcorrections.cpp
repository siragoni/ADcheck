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
void CrossSectionXNXNwithEMDcorrections(){
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
  // Double_t fI0N0N[3]        = { 0.0080, 0.0100, 0.0100 };
  // Double_t fI0NXN[3]        = { 0.0340, 0.0400, 0.0640 };
  // Double_t fIXN0N[3]        = { 1.7440, 0.7910, 0.4790 };
  // Double_t fIXNXN[3]        = { 0.2080, 0.1320, 0.0800 };
  // Double_t errfI0N0N[3]     = { 0.0046, 0.0025, 0.0038 };
  // Double_t errfI0NXN[3]     = { 0.0025, 0.0109, 0.0216 };
  // Double_t errfIXN0N[3]     = { 0.4339, 0.1007, 0.0730 };
  // Double_t errfIXNXN[3]     = { 0.0632, 0.0349, 0.0359 };

  // Perfect
  // Double_t fI0N0N[3]        = { 0.0060, 0.0080, 0.0090 };
  // Double_t fI0NXN[3]        = { 0.0360, 0.0460, 0.0690 };
  // Double_t fIXN0N[3]        = { 2.5840, 1.0650, 0.5280 };
  // Double_t fIXNXN[3]        = { 0.1960, 0.1630, 0.1050 };
  // Double_t errfI0N0N[3]     = { 0.0019, 0.0017, 0.0025 };
  // Double_t errfI0NXN[3]     = { 0.0216, 0.0104, 0.0194 };
  // Double_t errfIXN0N[3]     = { 1.2846, 0.1399, 0.0923 };
  // Double_t errfIXNXN[3]     = { 0.0881, 0.0337, 0.0336 };


  // Sidebands
  Double_t fI0N0N[3]        = { 0.0060, 0.0060, 0.0060 };
  Double_t fI0NXN[3]        = { 0.0210, 0.0320, 0.0520 };
  Double_t fIXN0N[3]        = { 1.6020, 0.6840, 0.3570 };
  Double_t fIXNXN[3]        = { 0.1490, 0.1120, 0.0773 };
  Double_t errfI0N0N[3]     = { 0.0019, 0.0005, 0.0009 };
  Double_t errfI0NXN[3]     = { 0.0191, 0.0249, 0.0174 };
  Double_t errfIXN0N[3]     = { 0.5834, 0.0738, 0.1116 };
  Double_t errfIXNXN[3]     = { 0.0793, 0.0271, 0.0288 };


  Double_t fD               = 0.055;
  Double_t eJPsi[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI             = 532.926;
  Double_t BR               = 0.05961;
  Double_t BeforeEmd0N0N[3] = { 2355, 6966, 2871 };
  Double_t BeforeEmd0NXN[3] = {  185,  630,  315 };
  Double_t BeforeEmdXN0N[3] = {  327,  937,  468 };
  Double_t BeforeEmdXNXN[3] = {  109,  444,  242 };


  Double_t BeforeEmd[3][4];
  BeforeEmd[0][0] = 2355;
  BeforeEmd[0][1] =  185;
  BeforeEmd[0][2] =  327;
  BeforeEmd[0][3] =  109;
  BeforeEmd[1][0] = 6966;
  BeforeEmd[1][1] =  630;
  BeforeEmd[1][2] =  937;
  BeforeEmd[1][3] =  444;
  BeforeEmd[2][0] = 2871;
  BeforeEmd[2][1] =  315;
  BeforeEmd[2][2] =  468;
  BeforeEmd[2][3] =  242;


  Double_t EMDcorrection[4][4];
  EMDcorrection[0][0] =  1.08;
  EMDcorrection[0][1] = -0.07809;
  EMDcorrection[0][2] = -0.07792;
  EMDcorrection[0][3] =  0.005637;
  EMDcorrection[1][0] = -0.03926;
  EMDcorrection[1][1] =  1.115;
  EMDcorrection[1][2] =  0.002934;
  EMDcorrection[1][3] = -0.08072;
  EMDcorrection[2][0] = -0.04157;
  EMDcorrection[2][1] =  0.03119;
  EMDcorrection[2][2] =  1.113;
  EMDcorrection[2][3] = -0.08074;
  EMDcorrection[3][0] =  0.001512;
  EMDcorrection[3][1] = -0.04306;
  EMDcorrection[3][2] = -0.04059;
  EMDcorrection[3][3] =  1.156;


  Double_t NumOfJPsi0N0N[3] = { 0,0,0 };
  Double_t NumOfJPsi0NXN[3] = { 0,0,0 };
  Double_t NumOfJPsiXN0N[3] = { 0,0,0 };
  Double_t NumOfJPsiXNXN[3] = { 0,0,0 };
  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ){
    for ( Int_t j = 0; j < 4; j++ ) {
      NumOfJPsi0N0N[iLoop] += EMDcorrection[0][j] * BeforeEmd[iLoop][j];
      NumOfJPsi0NXN[iLoop] += EMDcorrection[1][j] * BeforeEmd[iLoop][j];
      NumOfJPsiXN0N[iLoop] += EMDcorrection[2][j] * BeforeEmd[iLoop][j];
      NumOfJPsiXNXN[iLoop] += EMDcorrection[3][j] * BeforeEmd[iLoop][j];
    }
  }


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
  // Double_t y1Error0N0N[3] = { // 83 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 119/( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 98 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             72 /( (1+fI0N0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             100/( (1+fI0N0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             82 /( (1+fI0N0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             };
  // Double_t y1Error0NXN[3] = { // 18 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 34 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 32 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             18 /( (1+fI0NXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             34 /( (1+fI0NXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             32 /( (1+fI0NXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             };
  // Double_t y1ErrorXN0N[3] = { // 29 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 52 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 29 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             29 /( (1+fIXN0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             52 /( (1+fIXN0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             29 /( (1+fIXN0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             };
  // Double_t y1ErrorXNXN[3] = { // 16 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 25 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             // 24 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             16 /( (1+fIXNXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
  //                             25 /( (1+fIXNXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
  //                             24 /( (1+fIXNXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
  //                             };
  Double_t y1Error0N0N[3] = { TMath::Sqrt( ( 72.0* 72.0/(2355.0*2355.0)) + (errfI0N0N[0]*errfI0N0N[0]/((1+fI0N0N[0]+fD)*(1+fI0N0N[0]+fD))) )*DSigmaDy0N0N[0],
                              TMath::Sqrt( (100.0*100.0/(6966.0*6966.0)) + (errfI0N0N[1]*errfI0N0N[1]/((1+fI0N0N[1]+fD)*(1+fI0N0N[1]+fD))) )*DSigmaDy0N0N[1],
                              TMath::Sqrt( ( 82.0* 82.0/(2871.0*2871.0)) + (errfI0N0N[2]*errfI0N0N[2]/((1+fI0N0N[2]+fD)*(1+fI0N0N[2]+fD))) )*DSigmaDy0N0N[2]
                              };
  Double_t y1Error0NXN[3] = { TMath::Sqrt( ( 18.0* 18.0/( 185.0* 185.0)) + (errfI0NXN[0]*errfI0NXN[0]/((1+fI0NXN[0]+fD)*(1+fI0NXN[0]+fD))) )*DSigmaDy0NXN[0],
                              TMath::Sqrt( ( 34.0* 34.0/( 630.0* 630.0)) + (errfI0NXN[1]*errfI0NXN[1]/((1+fI0NXN[1]+fD)*(1+fI0NXN[1]+fD))) )*DSigmaDy0NXN[1],
                              TMath::Sqrt( ( 32.0* 32.0/( 315.0* 315.0)) + (errfI0NXN[2]*errfI0NXN[2]/((1+fI0NXN[2]+fD)*(1+fI0NXN[2]+fD))) )*DSigmaDy0NXN[2]
                              };
  Double_t y1ErrorXN0N[3] = { TMath::Sqrt( ( 29.0* 29.0/( 327.0* 327.0)) + (errfIXN0N[0]*errfIXN0N[0]/((1+fIXN0N[0]+fD)*(1+fIXN0N[0]+fD))) )*DSigmaDyXN0N[0],
                              TMath::Sqrt( ( 52.0* 52.0/( 937.0* 937.0)) + (errfIXN0N[1]*errfIXN0N[1]/((1+fIXN0N[1]+fD)*(1+fIXN0N[1]+fD))) )*DSigmaDyXN0N[1],
                              TMath::Sqrt( ( 29.0* 29.0/( 468.0* 468.0)) + (errfIXN0N[2]*errfIXN0N[2]/((1+fIXN0N[2]+fD)*(1+fIXN0N[2]+fD))) )*DSigmaDyXN0N[2]
                              };
  Double_t y1ErrorXNXN[3] = { TMath::Sqrt( ( 16.0* 16.0/( 109.0* 109.0)) + (errfIXNXN[0]*errfIXNXN[0]/((1+fIXNXN[0]+fD)*(1+fIXNXN[0]+fD))) )*DSigmaDyXNXN[0],
                              TMath::Sqrt( ( 25.0* 25.0/( 444.0* 444.0)) + (errfIXNXN[1]*errfIXNXN[1]/((1+fIXNXN[1]+fD)*(1+fIXNXN[1]+fD))) )*DSigmaDyXNXN[1],
                              TMath::Sqrt( ( 24.0* 24.0/( 242.0* 242.0)) + (errfIXNXN[2]*errfIXNXN[2]/((1+fIXNXN[2]+fD)*(1+fIXNXN[2]+fD))) )*DSigmaDyXNXN[2]
                              };
  Double_t x1Error[3]     = {  (4.0-2.5)/6.0, (4.0-2.5)/6.0, (4.0-2.5)/6.0 };

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



  new TCanvas;
  Double_t Ratio[3]  = { 0,0,0 };
  Double_t Error[3]  = { 0,0,0 };
  for (size_t i = 0; i < 3; i++) {
    Ratio[i] = DSigmaDy0NXN[i] / DSigmaDyXN0N[i];
    Error[i] = y1Error0NXN[i]  + y1ErrorXN0N[i];
  }
  TGraphErrors* RatioPlot = new TGraphErrors(3, x1, Ratio, x1Error, Error);
  RatioPlot->Draw("APL");
  RatioPlot->GetXaxis()->SetTitle("y");
  RatioPlot->GetYaxis()->SetTitle("Ratio [au]");
  // Change the axis limits
  gPad->BuildLegend();
  gPad->Modified();

}
