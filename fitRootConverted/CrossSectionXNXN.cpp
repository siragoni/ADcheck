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
void CrossSectionXNXN(){
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
  // Double_t fI0N0N[3]        = { 0.0430, 0.027,  0.043 };
  // Double_t fI0NXN[3]        = { 0.0840, 0.084,  0.084 };
  // Double_t fIXN0N[3]        = { 5.6430, 2.335,  1.511 };
  // Double_t fIXNXN[3]        = { 0.5000, 0.485,  0.500 };
  // Double_t fI0N0N[3]        = { 0.008, 0.009, 0.010 };
  // Double_t fI0NXN[3]        = { 0.033, 0.039, 0.061 };
  // Double_t fIXN0N[3]        = { 1.648, 0.817, 0.424 };
  // Double_t fIXNXN[3]        = { 0.198, 0.129, 0.087 };
  // Double_t errfI0N0N[3]     = { 0.0000016, 0.0000022, 0.0000024 };
  // Double_t errfI0NXN[3]     = { 0.0000087, 0.0000102, 0.0000174 };
  // Double_t errfIXN0N[3]     = { 0.0002170, 0.0001043, 0.0000425 };
  // Double_t errfIXNXN[3]     = { 0.0000198, 0.0000185, 0.0000133 };
  Double_t fI0N0N[3]        = { 0.0060, 0.0080, 0.0090 };
  Double_t fI0NXN[3]        = { 0.0360, 0.0460, 0.0690 };
  Double_t fIXN0N[3]        = { 2.5840, 1.0650, 0.5280 };
  Double_t fIXNXN[3]        = { 0.1960, 0.1630, 0.1050 };
  Double_t errfI0N0N[3]     = { 0.0019, 0.0017, 0.0025 };
  Double_t errfI0NXN[3]     = { 0.0216, 0.0104, 0.0194 };
  Double_t errfIXN0N[3]     = { 1.2846, 0.1399, 0.0923 };
  Double_t errfIXNXN[3]     = { 0.0881, 0.0337, 0.0336 };

  Double_t fD               = 0.055;
  Double_t eJPsi[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI             = 532.926;
  Double_t BR               = 0.05961;
  Double_t NumOfJPsi0N0N[3] = { 2355, 6966, 2871 };
  Double_t NumOfJPsi0NXN[3] = {  185,  630,  315 };
  Double_t NumOfJPsiXN0N[3] = {  327,  937,  468 };
  Double_t NumOfJPsiXNXN[3] = {  109,  444,  242 };

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
  Double_t y1Error0N0N[3] = { TMath::Sqrt( ( 72* 72/(2355*2355)) + (errfI0N0N[0]*errfI0N0N[0]/((1+fI0N0N[0]+fD)*(1+fI0N0N[0]+fD))) )*DSigmaDy0N0N[0],
                              TMath::Sqrt( (100*100/(6966*6966)) + (errfI0N0N[1]*errfI0N0N[1]/((1+fI0N0N[1]+fD)*(1+fI0N0N[1]+fD))) )*DSigmaDy0N0N[1],
                              TMath::Sqrt( ( 82* 82/(2871*2871)) + (errfI0N0N[2]*errfI0N0N[2]/((1+fI0N0N[2]+fD)*(1+fI0N0N[2]+fD))) )*DSigmaDy0N0N[2]
                              };
  Double_t y1Error0NXN[3] = { TMath::Sqrt( ( 18* 18/( 185* 185)) + (errfI0NXN[0]*errfI0NXN[0]/((1+fI0NXN[0]+fD)*(1+fI0NXN[0]+fD))) )*DSigmaDy0NXN[0],
                              TMath::Sqrt( ( 34* 34/( 630* 630)) + (errfI0NXN[1]*errfI0NXN[1]/((1+fI0NXN[1]+fD)*(1+fI0NXN[1]+fD))) )*DSigmaDy0NXN[1],
                              TMath::Sqrt( ( 32* 32/( 315* 315)) + (errfI0NXN[2]*errfI0NXN[2]/((1+fI0NXN[2]+fD)*(1+fI0NXN[2]+fD))) )*DSigmaDy0NXN[2]
                              };
  Double_t y1ErrorXN0N[3] = { TMath::Sqrt( ( 29* 29/( 327* 327)) + (errfIXN0N[0]*errfIXN0N[0]/((1+fIXN0N[0]+fD)*(1+fIXN0N[0]+fD))) )*DSigmaDyXN0N[0],
                              TMath::Sqrt( ( 52* 52/( 937* 937)) + (errfIXN0N[1]*errfIXN0N[1]/((1+fIXN0N[1]+fD)*(1+fIXN0N[1]+fD))) )*DSigmaDyXN0N[1],
                              TMath::Sqrt( ( 29* 29/( 468* 468)) + (errfIXN0N[2]*errfIXN0N[2]/((1+fIXN0N[2]+fD)*(1+fIXN0N[2]+fD))) )*DSigmaDyXN0N[2]
                              };
  Double_t y1ErrorXNXN[3] = { TMath::Sqrt( ( 16* 16/( 109* 109)) + (errfIXNXN[0]*errfIXNXN[0]/((1+fIXNXN[0]+fD)*(1+fIXNXN[0]+fD))) )*DSigmaDyXNXN[0],
                              TMath::Sqrt( ( 25* 25/( 444* 444)) + (errfIXNXN[1]*errfIXNXN[1]/((1+fIXNXN[1]+fD)*(1+fIXNXN[1]+fD))) )*DSigmaDyXNXN[1],
                              TMath::Sqrt( ( 24* 24/( 242* 242)) + (errfIXNXN[2]*errfIXNXN[2]/((1+fIXNXN[2]+fD)*(1+fIXNXN[2]+fD))) )*DSigmaDyXNXN[2]
                              };
  Double_t x1Error[3]     = {  (4-2.5)/6, (4-2.5)/6, (4-2.5)/6 };

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
