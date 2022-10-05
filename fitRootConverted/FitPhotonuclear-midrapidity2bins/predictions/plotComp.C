//
// this program compares the cross sections from Simone with our predictions
//

// -----------------------------------------------------------------
// all headers are defined here

// c++ headers
#include <iostream>
#include <fstream>
#include <cmath>

// root headers
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TPad.h>
#include <TStyle.h>

// -----------------------------------------------------------------
// data and models

// alice data
const int ALICE_photo_n = 9;
Double_t ALICE_photo_W[ALICE_photo_n]   = { 19.130700000, 24.5643000, 31.541200000, 493.3870000, 633.5210000, 813.4570000, 97.1536, 124.748, 160.179 };
Double_t ALICE_photo_sig[ALICE_photo_n] = {  0.00882072,  0.0138693,  0.0168523,   0.0473069,   0.0523713,   0.0628193, 0.0218977, 0.0230963,  0.0246617};
Double_t ALICE_photo_sta[ALICE_photo_n]  = {  0.000260183,0.000307948,0.000593475,   0.00775179,   0.00799204,   0.024678, 0.00101229, 0.00136407, 0.00506907};


// read model from PHYSICAL REVIEW C 97, 024901 (2018), GG-hs (red solid line)
void readHS(TGraph *gr)
{
  // square of jpsi mass
  const double jpsiMass2 = 3.0969*3.0969;
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;

  // reading input file
  ifstream ifs("fig2_sigma_GGhs_xBj.txt");
  const int nlines = 41;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> x;
    ifs >> sig;
    // convert x into W
    double W = std::sqrt(jpsiMass2/x);
    // fill TGraph
    gr->SetPoint(gr->GetN(),W,sig);
  }
  ifs.close();
}

// read model from Phys.Lett.B 817 (2021) 136306, b-BK-A and b-BK--GG
void readBK(TGraph *grA, TGraph *grGG)
{
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;
  double W = 0;
  
  // reading first file
  const int nlines = 68;
  ifstream ifs("Pb_GG_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> x;
    ifs >> W;
    ifs >> sig;
    // fill TGraph
    grGG->SetPoint(grGG->GetN(),W,sig);
  }
  ifs.close();

  // reading first file
  ifstream ifsA("Pb_A_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifsA >> x;
    ifsA >> W;
    ifsA >> sig;
    // fill TGraph
    grA->SetPoint(grA->GetN(),W,sig);
  }
  ifsA.close();

}

// -----------------------------------------------------------------
// plotting program

void plotComp()
{
  // convert alice data to mub
  for(int i=0;i<ALICE_photo_n;i++) {
    ALICE_photo_sig[i] *= 1000;
    ALICE_photo_sta[i] *= 1000;
  };
  
  // set graphs
  TGraphErrors *dataALICE = new TGraphErrors(ALICE_photo_n,ALICE_photo_W,ALICE_photo_sig,NULL,ALICE_photo_sta);
  TGraph *predHS = new TGraph(); // prediction from the hot-spot model
  readHS(predHS);
  TGraph *predBKA = new TGraph(); // prediction from the b-BK equation
  TGraph *predBKGG = new TGraph(); // prediction from the b-BK equation  
  readBK(predBKA,predBKGG);
  
  // define canvas and pad position/attributes
  TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  TH1F* frame1 = gPad->DrawFrame(10,1,1000,1000);
  gPad->SetLeftMargin(0.09);  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.11);  gPad->SetBottomMargin(0.13);
  gPad->SetLogx();  gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  
  frame1->SetTitle(";W_{#gammaPb} (GeV);#sigma(#gamma+p #rightarrow J/#psi+p) (#mub)");
  Float_t siz = 0.045;
  frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
  frame1->SetTitleSize(siz, "X");  frame1->SetLabelSize(siz, "X");
  frame1->GetXaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetTitleOffset(1.0);  

  // plot data
  dataALICE->SetMarkerStyle(kFullCircle);
  dataALICE->SetMarkerColor(kBlack);  
  dataALICE->Draw("p,e1,same");

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

  // add legends
  TLegend* leg1 = new TLegend(0.2, 0.6, 0.4, 0.8);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(siz);
  leg1->AddEntry(dataALICE,"ALICE","pe");
  leg1->AddEntry(predHS,"Hot-spot Model (PRC97 (2018) 024901)","l");
  leg1->AddEntry(predBKA,"b-BK-A (PLB817(2021)136306)","l");
  leg1->AddEntry(predBKGG,"b-BK-GG (PLB817(2021)136306)","l");  
  leg1->Draw("same");

}
