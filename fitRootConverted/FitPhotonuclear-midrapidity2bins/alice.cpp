
/*

  plot of sigma(gamma-p) for semiforward and mid-rapidity paper, Alice data,
  other experiments and models

  compile:

g++ alice.cpp `root-config --libs --cflags` -lEG -lMinuit -o alice

  run:

./alice

*/

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

// local headers
// #include "alice.h"
#include "alice2.h"






// read model from PHYSICAL REVIEW C 97, 024901 (2018), GG-hs (red solid line)
void readHS(TGraph *gr)
{
  // square of jpsi mass
  const double jpsiMass2 = 3.0969*3.0969;
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;

  // reading input file
  ifstream ifs("predictions/fig2_sigma_GGhs_xBj.txt");
  const int nlines = 41;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> x;
    ifs >> sig;
    // convert x into W
    double W = std::sqrt(jpsiMass2/x);
    // fill TGraph
    gr->SetPoint(gr->GetN(),W,sig/1000.);
  }
  ifs.close();
}
//_____________________________________________________________________________
// read model from Phys.Lett.B 817 (2021) 136306, b-BK-A and b-BK--GG
void readBK(TGraph *grA, TGraph *grGG)
{
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;
  double W = 0;
  
  // reading first file
  const int nlines = 68;
  ifstream ifs("predictions/Pb_GG_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> x;
    ifs >> W;
    ifs >> sig;
    // fill TGraph
    grGG->SetPoint(grGG->GetN(),W,sig/1000.);
  }
  ifs.close();

  // reading first file
  ifstream ifsA("predictions/Pb_A_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifsA >> x;
    ifsA >> W;
    ifsA >> sig;
    // fill TGraph
    grA->SetPoint(grA->GetN(),W,sig/1000.);
  }
  ifsA.close();

}
//_____________________________________________________________________________
// read model from STARlight Impulse Approximation
void readImpulseApproximation(TGraph *gr)
{
  // square of jpsi mass
  const double jpsiMass2 = 3.0969*3.0969;
  
  // tmp storage of input data
  double y = 0;
  double W = 0;
  double flux = 0;
  double sig = 0;
  double dsig = 0;
  double W2 = 0;
  double flux2 = 0;
  double sig2 = 0;
  double dsig2 = 0;
  double boh = 0;

  // reading input file
  ifstream ifs("predictions/ImpulseApproximation.txt");
  const int nlines = 810;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> y;
    ifs >> W;
    ifs >> flux;
    ifs >> sig;
    ifs >> dsig;
    ifs >> W2;
    ifs >> flux2;
    ifs >> sig2;
    ifs >> dsig2;
    ifs >> boh;
    // fill TGraph
    if ( i < 3 ){ 
    cout << "y     = " << y     << endl;
    cout << "W     = " << W     << endl;
    cout << "flux  = " << flux  << endl;
    cout << "sig   = " << sig   << endl;
    cout << "dsig  = " << dsig  << endl;
    cout << "W2    = " << W2    << endl;
    cout << "flux2 = " << flux2 << endl;
    cout << "sig2  = " << sig2  << endl;
    cout << "dsig2 = " << dsig2 << endl;
    cout << "boh   = " << dsig2 << endl;
    cout << "=================" << endl;
    }
    gr->SetPoint(gr->GetN(),W,sig);
  }
  ifs.close();
}
//_____________________________________________________________________________
// read model STARlight 
void readSTARlight(TGraph *gr)
{
  // square of jpsi mass
  const double jpsiMass2 = 3.0969*3.0969;
  
  // tmp storage of input data
  double y = 0;
  double W = 0;
  double flux = 0;
  double sig = 0;
  double dsig = 0;
  double W2 = 0;
  double flux2 = 0;
  double sig2 = 0;
  double dsig2 = 0;
  double boh = 0;

  // reading input file
  ifstream ifs("predictions/STARlight-XNXN.txt");
  const int nlines = 810;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> y;
    ifs >> W;
    ifs >> flux;
    ifs >> sig;
    ifs >> dsig;
    ifs >> W2;
    ifs >> flux2;
    ifs >> sig2;
    ifs >> dsig2;
    ifs >> boh;
    // fill TGraph
    if ( i < 3 ){ 
    cout << "y     = " << y     << endl;
    cout << "W     = " << W     << endl;
    cout << "flux  = " << flux  << endl;
    cout << "sig   = " << sig   << endl;
    cout << "dsig  = " << dsig  << endl;
    cout << "W2    = " << W2    << endl;
    cout << "flux2 = " << flux2 << endl;
    cout << "sig2  = " << sig2  << endl;
    cout << "dsig2 = " << dsig2 << endl;
    cout << "boh   = " << dsig2 << endl;
    cout << "=================" << endl;
    }
    gr->SetPoint(gr->GetN(),W,sig/1000);
  }
  ifs.close();
}
//_____________________________________________________________________________
int main(void) {

  //configure layout of the plot
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02,"Y");

  TCanvas* c1 = new TCanvas("c1","c1",1000,700);

  TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  frame1->GetXaxis()->SetMoreLogLabels();
  frame1->GetXaxis()->SetTitleOffset(1.25);
  frame1->GetYaxis()->SetTitleOffset(1.25);
  frame1->GetXaxis()->SetTitleSize(0.045);
  frame1->GetYaxis()->SetTitleSize(0.045);
  frame1->GetXaxis()->SetLabelSize(0.045);
  frame1->GetYaxis()->SetLabelSize(0.045);
  frame1->GetXaxis()->SetTitleFont(42);
  frame1->GetYaxis()->SetTitleFont(42);
  frame1->GetXaxis()->SetLabelFont(42);
  frame1->GetYaxis()->SetLabelFont(42);

  // gPad->SetLeftMargin(0.09);
  // gPad->SetRightMargin(0.01);
  // gPad->SetTopMargin(0.11);//0.025
  // gPad->SetBottomMargin(0.13);
  // gPad->SetLogx();
  // gPad->SetLogy();
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTopMargin(0.14);
  // gPad->SetGridx();
  // gPad->SetGridy();
  gStyle->SetOptStat(0);


  gStyle->SetOptStat("");
  gStyle->SetPalette(1);
  gStyle->SetLineWidth(2);      //axis line
  gStyle->SetFrameLineWidth(2); //frame line
  TGaxis::SetMaxDigits(3);

  frame1->SetTitle(";W_{#gammaPb} (GeV);#sigma(#gammaPb) (mb)");
  Float_t siz = 0.045;
  frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
  frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
  // frame1->GetXaxis()->SetTitleOffset(1.4);

  //upper horizontal axis with Bjorken-x
  Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  Double_t bxmin = TMath::Power((mjpsi/gxmax),2.);
  Double_t bxmax = TMath::Power((mjpsi/gxmin),2.);
  TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  fbx->SetParameter(0, mjpsi);

  TGaxis *axis = new TGaxis(gxmax, gymax, gxmin, gymax, "fbx", 510, "+G");
  axis->SetTextFont(42);
  axis->SetLabelFont(42);
  axis->SetTitleSize(siz); axis->SetLabelSize(siz);
  axis->SetLabelOffset(-0.035);
  axis->SetTitleOffset();
  axis->SetTitle("Bjorken-#it{x}");
  axis->Draw("same");

  // TLatex *bxtit = new TLatex();
  // bxtit->SetTextFont(42);
  // bxtit->SetTextSize(siz);
  // bxtit->SetTextAlign(31);
  // bxtit->DrawLatex(gxmax, gymax+450, "Bjorken-#it{x}");
  //end of layout configuration

  //alice fit
  alice *alic = alice::instance();

  TGraph *gShade = alic->fit();
  // gShade->Draw("fsame");


  // TGraphErrors *STARLIGHTpred      = alic->getStarlightWithNuclear();
  // STARLIGHTpred->Draw("psame");


  //alice data
  TGraphErrors *aliceData      = alic->getData();
  aliceData->Draw("psame");


  // TGraphErrors *STARLIGHTpred2      = alic->getStarlight();
  // STARLIGHTpred2->Draw("psame");





  TGraph *predHS = new TGraph(); // prediction from the hot-spot model
  readHS(predHS);
  TGraph *predBKA = new TGraph(); // prediction from the b-BK equation
  TGraph *predBKGG = new TGraph(); // prediction from the b-BK equation  
  readBK(predBKA,predBKGG);
  // plot prediction from the hot-spot model
  predHS->SetLineColor(kRed);
  predHS->SetLineWidth(2);
  predHS->Draw("c,same");
  // plot prediction from the b-BK-A model
  predBKA->SetLineColor(kBlue);
  predBKA->SetLineWidth(2);
  predBKA->Draw("c,same");
  // plot prediction from the b-BK-GG model
  predBKGG->SetLineColor(kCyan);
  predBKGG->SetLineWidth(2);
  predBKGG->Draw("c,same");




  TGraph *predIA = new TGraph(); // prediction from the hot-spot model
  readImpulseApproximation(predIA);
  // plot prediction from STARlight
  predIA->SetLineColor(kMagenta);
  predIA->SetLineWidth(2);
  predIA->Draw("c,same");
  TGraph *predSL = new TGraph(); // prediction from the hot-spot model
  readSTARlight(predSL);
  // plot prediction from STARlight
  predSL->SetLineColor(kBlack);
  predSL->SetLineWidth(2);
  predSL->Draw("c,same");



  //legend for experimental data
  Double_t xl = 0.15, dxl = 0.2;
  Double_t yl = 0.52, dyl = 0.3;
  TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(siz-0.005);
  leg1->AddEntry(aliceData,"ALICE data","p");
  // leg1->AddEntry(STARLIGHTpred,"STARlight","p");
  leg1->AddEntry(predSL,"STARlight","l");
  // leg1->AddEntry(STARLIGHTpred2,"Impulse Approximation","p");
  leg1->AddEntry(predIA,"Impulse Approximation","l");
  leg1->AddEntry(predHS,"Hot-spot Model (PRC97 (2018) 024901)","l");
  leg1->AddEntry(predBKA,"b-BK-A (PLB817(2021)136306)","l");
  leg1->AddEntry(predBKGG,"b-BK-GG (PLB817(2021)136306)","l");  
  leg1->Draw("same");


  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.43,0.80,"LHC18qr, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.43,0.33,"N = 0.0218 #pm 0.0006");
  // latex5->DrawLatex(0.43,0.27,"#delta = 0.3349 #pm 0.0334");
  // latex5->DrawLatex(0.43,0.21,"#chi^{2}/NDF = 15.92 / (9 - 2) = 2.27");



  //put vertical scale atop the models
  gPad->RedrawAxis("Y");

  c1->SaveAs("fit.pdf");

  return 0;

}//main
