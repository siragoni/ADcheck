// #include "Riostream.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVectorD.h"

Double_t FittingPhotonuclear(Int_t element, Int_t mode = 0, Int_t rap = 0);
void     Computing          ();
void     BeautifyPad          ();
void     BeautifyHisto      (TH1* histogram);

TFile* file = new TFile("../../Michal-Broz-xsec/xSection_Cent.root");
TGraphAsymmErrors* Cent0N0N = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0n0n");
TGraphAsymmErrors* Cent0NXN = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0nXn");
TGraphAsymmErrors* CentXNXN = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_XnXn");
TFile* file2 = new TFile("xSection_Cent_2bins.root");
TGraphAsymmErrors* Cent0N0N_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0n0n");
TGraphAsymmErrors* Cent0NXN_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0nXn");
TGraphAsymmErrors* CentXNXN_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_XnXn");

//______________________________________________
void BeautifyPad(){
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
}
//______________________________________________
void BeautifyHisto(TH1* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");
}

// https://root-forum.cern.ch/t/tdecompsvd-not-working/22313
//_____________________________________
Double_t FittingPhotonuclear(Int_t element, Int_t mode = 0, Int_t rap = 0)
{

  //==========================
  // SYSTEMATIC UNCERTAINTIES
  //==========================
  // Correlated ones
  //--------------------------
  Double_t Lumi = 2.7;
  Double_t BR   = 0.6;
  Double_t Trk  = 3.0;
  Double_t Trg  = 6.2;
  Double_t Mtch = 1.0;
  Double_t fD   = 0.7;
  Double_t PFlux= 2.0;
  Double_t Correlated = TMath::Sqrt(Lumi*Lumi+PFlux*PFlux+BR*BR+Trg*Trg+Trk*Trk+Mtch*Mtch+fD*fD);
  //--------------------------
  // Yields
  Double_t SignalExtr4to25[3] = {0.36, 0.35, 0.74};
  Double_t SignalExtr4to35[3] = {0.54, 0.53, 0.90};
  Double_t SignalExtr35to3[3] = {0.10, 0.63, 0.67};
  Double_t SignalExtr3to25[3] = {0.20, 1.26, 0.81};
  //--------------------------
  // Incoherent
  Double_t fI4to25[3] = {0.41, 0.85, 3.07};
  Double_t fI4to35[3] = {0.37, 0.53, 2.23};
  Double_t fI35to3[3] = {0.41, 0.90, 3.32};
  Double_t fI3to25[3] = {0.38, 0.61, 1.59};
  //--------------------------
  // Pile-up
  Double_t PileUp4to25[3] = {0.12,-0.68,-0.11};
  Double_t PileUp4to35[3] = {0.12,-0.82,-0.11};
  Double_t PileUp35to3[3] = {0.12,-0.70,-0.11};
  Double_t PileUp3to25[3] = {0.12,-0.56,-0.12};
  //--------------------------
  // ZDC eff
  Double_t eZdc4to25[3] = {0.71,-2.56,-11.10};
  Double_t eZdc4to35[3] = {0.52, 1.33,-11.13};
  Double_t eZdc35to3[3] = {0.68, 2.49,-11.11};
  Double_t eZdc3to25[3] = {0.96, 3.37,-11.09};


  //==========================
  // Final computation
  //--------------------------
  Double_t PartialSys4to25[3] = {Correlated, Correlated, Correlated};
  Double_t PartialSys4to35[3] = {Correlated, Correlated, Correlated};
  Double_t PartialSys35to3[3] = {Correlated, Correlated, Correlated};
  Double_t PartialSys3to25[3] = {Correlated, Correlated, Correlated};
  for (Int_t i = 0; i < 3; i++) {
    PartialSys4to25[i] = TMath::Sqrt(Correlated*Correlated + SignalExtr4to25[i]*SignalExtr4to25[i] + fI4to25[i]*fI4to25[i]);
    PartialSys4to35[i] = TMath::Sqrt(Correlated*Correlated + SignalExtr4to35[i]*SignalExtr4to35[i] + fI4to35[i]*fI4to35[i]);
    PartialSys35to3[i] = TMath::Sqrt(Correlated*Correlated + SignalExtr35to3[i]*SignalExtr35to3[i] + fI35to3[i]*fI35to3[i]);
    PartialSys3to25[i] = TMath::Sqrt(Correlated*Correlated + SignalExtr3to25[i]*SignalExtr3to25[i] + fI3to25[i]*fI3to25[i]);
  }




  //==========================
  // Renormalising them
  //--------------------------
  for (Int_t i = 0; i < 3; i++) {

    PartialSys4to25[i] = PartialSys4to25[i] * 0.01;
    PartialSys4to35[i] = PartialSys4to35[i] * 0.01;
    PartialSys35to3[i] = PartialSys35to3[i] * 0.01;
    PartialSys3to25[i] = PartialSys3to25[i] * 0.01;

    PileUp4to25[i] = PileUp4to25[i] * 0.01;
    PileUp4to35[i] = PileUp4to35[i] * 0.01;
    PileUp35to3[i] = PileUp35to3[i] * 0.01;
    PileUp3to25[i] = PileUp3to25[i] * 0.01;

    eZdc4to25[i] = eZdc4to25[i] * 0.01;
    eZdc4to35[i] = eZdc4to35[i] * 0.01;
    eZdc35to3[i] = eZdc35to3[i] * 0.01;
    eZdc3to25[i] = eZdc3to25[i] * 0.01;

  }





  //==========================
  // Solving linear equations
  // with SVD technique
  //==========================
  //
  // We have an equation like
  // b = Ax
  // where b are the measured
  // cross sections,
  // A is the flux matrix
  // x is our unknown vector
  // of photonuclear sigmas
  // for positive and
  // negative rapidity


  // 4 neutron emission classes
  // 2 unknown photonuclear xsecs
  // TMatrixD Amatrix(4,2);
  // TMatrixD Amatrix0(4,2);
  // TMatrixD Amatrix1(4,2);
  // TMatrixD Amatrix2(4,2);
  TMatrixD Amatrix(3,2);
  TMatrixD Amatrix0(3,2);
  TMatrixD Amatrix1(3,2);
  TMatrixD Amatrix2(3,2);
  Amatrix.Zero();
  Amatrix0.Zero();
  Amatrix1.Zero();
  Amatrix2.Zero();
  //=============================
  // MIDRAPIDITY CASE
  //-----------------------------
  // 3 neutron emission classes
  // 2 unknown photonuclear xsecs
  TMatrixD Amatrix_mid(3,2);
  TMatrixD Amatrix0_mid(3,2);
  TMatrixD Amatrix1_mid(3,2);
  TMatrixD Amatrix2_mid(3,2);
  Amatrix_mid.Zero();
  Amatrix0_mid.Zero();
  Amatrix1_mid.Zero();
  Amatrix2_mid.Zero();


  TMatrixD Amatrix0_mid_twobins(3,2);
  TMatrixD Amatrix1_mid_twobins(3,2);
  Amatrix0_mid_twobins.Zero();
  Amatrix1_mid_twobins.Zero();


  // Fluxes from Michal
  /*
  y = 0 0n0n:  Low = 61.537 High = 61.537
  y = 0 0nXn:  Low = 16.1841 High = 16.1841
  y = 0 XnXn:  Low = 5.03313 High = 5.03313
  y = 0.25 0n0n:  Low = 68.7589 High = 54.2806
  y = 0.25 0nXn:  Low = 16.3828 High = 15.8818
  y = 0.25 XnXn:  Low = 5.05204 High = 4.99478
  y = 0.575 0n0n:  Low = 78.3588 High = 45.2674
  y = 0.575 0nXn:  Low = 16.5604 High = 15.3431
  y = 0.575 XnXn:  Low = 5.06777 High = 4.92338
  y = 2.75 0n0n:  Low = 142.948 High = 3.66191
  y = 2.75 0nXn:  Low = 16.7132 High = 4.28895
  y = 2.75 XnXn:  Low = 5.04786 High = 2.04812
  y = 3.25 0n0n:  Low = 157.789 High = 1.07192
  y = 3.25 0nXn:  Low = 16.6849 High = 1.70239
  y = 3.25 XnXn:  Low = 5.03916 High = 0.949336
  y = 3.75 0n0n:  Low = 172.551 High = 0.200768
  y = 3.75 0nXn:  Low = 16.6506 High = 0.421108
  y = 3.75 XnXn:  Low = 5.02952 High = 0.275626
  */


  // -4 < y < -2.5
  Amatrix(0,0) = 162.75;
  Amatrix(0,1) = 1.0973;
  Amatrix(1,0) = 18.313*0.5;
  Amatrix(1,1) = 1.9936*0.5;
  Amatrix(2,0) = 6.5209;
  Amatrix(2,1) = 1.3546;
  // -4 < y < -3.5
  Amatrix0(0,0) = 178.28;
  Amatrix0(0,1) = 0.2149;
  Amatrix0(1,0) = 18.312*0.5;
  Amatrix0(1,1) = 0.52656*0.5;
  Amatrix0(2,0) = 6.5199;
  Amatrix0(2,1) = 0.425;
  // -3.5 < y < -3.0
  Amatrix1(0,0) = 162.75;
  Amatrix1(0,1) = 1.0973;
  Amatrix1(1,0) = 18.313*0.5;
  Amatrix1(1,1) = 1.9936*0.5;
  Amatrix1(2,0) = 6.5209;
  Amatrix1(2,1) = 1.3546;
  // -3.0 < y < -2.5
  Amatrix2(0,0) = 147.23;
  Amatrix2(0,1) = 3.6998;
  Amatrix2(1,0) = 18.311*0.5;
  Amatrix2(1,1) = 4.8253*0.5;
  Amatrix2(2,0) = 6.5223;
  Amatrix2(2,1) = 2.7674;


  TDecompSVD svd(Amatrix,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd0(Amatrix0,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd1(Amatrix1,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd2(Amatrix2,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  Bool_t ok;
  Bool_t ok0;
  Bool_t ok1;
  Bool_t ok2;








  //=============================
  // MIDRAPIDITY MEASUREMENT
  //-----------------------------
  // y = 0 0n0n:  Low = 61.537 High = 61.537
  // y = 0 0nXn:  Low = 16.1841 High = 16.1841
  // y = 0 XnXn:  Low = 5.03313 High = 5.03313
  // y = 0.25 0n0n:  Low = 68.7589 High = 54.2806
  // y = 0.25 0nXn:  Low = 16.3828 High = 15.8818
  // y = 0.25 XnXn:  Low = 5.05204 High = 4.99478
  // y = 0.575 0n0n:  Low = 78.3588 High = 45.2674
  // y = 0.575 0nXn:  Low = 16.5604 High = 15.3431
  // y = 0.575 XnXn:  Low = 5.06777 High = 4.92338
  // y = 0.575
  Amatrix0_mid(0,0) = 78.3588;
  Amatrix0_mid(0,1) = 45.2674;
  Amatrix0_mid(1,0) = 16.5604;
  Amatrix0_mid(1,1) = 15.3431;
  Amatrix0_mid(2,0) = 5.06777;
  Amatrix0_mid(2,1) = 4.92338;
  // y = 0.25
  Amatrix1_mid(0,0) = 68.7589;
  Amatrix1_mid(0,1) = 54.2806;
  Amatrix1_mid(1,0) = 16.3828;
  Amatrix1_mid(1,1) = 15.8818;
  Amatrix1_mid(2,0) = 5.05204;
  Amatrix1_mid(2,1) = 4.99478;
  // y = 0.
  Amatrix2_mid(0,0) = 61.537;
  Amatrix2_mid(0,1) = 61.537;
  Amatrix2_mid(1,0) = 16.1841;
  Amatrix2_mid(1,1) = 16.1841;
  Amatrix2_mid(2,0) = 5.03313;
  Amatrix2_mid(2,1) = 5.03313;


  // TDecompSVD svd_mid(Amatrix_mid,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd0_mid(Amatrix0_mid,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd1_mid(Amatrix1_mid,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd2_mid(Amatrix2_mid,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  // Bool_t ok_mid;
  Bool_t ok0_mid;
  Bool_t ok1_mid;
  Bool_t ok2_mid;



  // Fluxes from Michal 2
  // y = 0 0n0n:  Low = 61.537 High = 61.537
  // y = 0 0nXn:  Low = 16.1841 High = 16.1841
  // y = 0 XnXn:  Low = 5.03313 High = 5.03313
  // y = 0.5 0n0n:  Low = 76.1357 High = 47.3637
  // y = 0.5 0nXn:  Low = 16.5274 High = 15.494
  // y = 0.5 XnXn:  Low = 5.06535 High = 4.94607
  // y = 0.5
  Amatrix0_mid_twobins(0,0) = 77.709;
  Amatrix0_mid_twobins(0,1) = 47.873;
  Amatrix0_mid_twobins(1,0) = 17.962;
  Amatrix0_mid_twobins(1,1) = 16.789;
  Amatrix0_mid_twobins(2,0) = 6.4902;
  Amatrix0_mid_twobins(2,1) = 6.3229;

  // y = 0.
  Amatrix1_mid_twobins(0,1) = 62.374;
  Amatrix1_mid_twobins(0,0) = 62.374;
  Amatrix1_mid_twobins(1,1) = 17.565;
  Amatrix1_mid_twobins(1,0) = 17.565;
  Amatrix1_mid_twobins(2,1) = 6.4383;
  Amatrix1_mid_twobins(2,0) = 6.4383;
  TDecompSVD svd0_mid_2(Amatrix0_mid_twobins,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd1_mid_2(Amatrix1_mid_twobins,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  Bool_t ok0_mid_2;
  Bool_t ok1_mid_2;





  //====================
  // Corrected sigmas
  // -4 < y < -2.5
  //====================
  // CorrSigma  stat.
  // 2.2390    0.0307
  // 0.1184    0.0041
  // 0.1667    0.0074
  // 0.1695    0.0088
  TVectorD xsec4to25orig(3);
  xsec4to25orig.Zero();
  xsec4to25orig(0) = 2.2390;
  xsec4to25orig(1) = 0.1667;
  xsec4to25orig(2) = 0.1695;
  Double_t Total4to25 = 2.772690;
  Double_t StatTotal4to25 = 0.022*TMath::Sqrt(1.5);
  TRandom R;
  TVectorD xsec4to25(3);
  xsec4to25.Zero();
  if(         mode == 0 ){
    xsec4to25(0) = xsec4to25orig(0);
    xsec4to25(1) = xsec4to25orig(1);
    xsec4to25(2) = xsec4to25orig(2);
  } else if ( mode == 1) {
    xsec4to25(0) = xsec4to25orig(0)*(1. + PartialSys4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. + PartialSys4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. + PartialSys4to25[2]);
  } else if ( mode == 2) {
    xsec4to25(0) = xsec4to25orig(0)*(1. - PartialSys4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. - PartialSys4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. - PartialSys4to25[2]);
  } else if ( mode == 3) {
    xsec4to25(0) = xsec4to25orig(0)*(1. + PileUp4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. + PileUp4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. + PileUp4to25[2]);
  } else if ( mode == 4) {
    xsec4to25(0) = xsec4to25orig(0)*(1. - PileUp4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. - PileUp4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. - PileUp4to25[2]);
  } else if ( mode == 5) {
    xsec4to25(0) = xsec4to25orig(0)*(1. + eZdc4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. + eZdc4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. + eZdc4to25[2]);
  } else if ( mode == 6) {
    xsec4to25(0) = xsec4to25orig(0)*(1. - eZdc4to25[0]);
    xsec4to25(1) = xsec4to25orig(1)*(1. - eZdc4to25[1]);
    xsec4to25(2) = xsec4to25orig(2)*(1. - eZdc4to25[2]);
  } else {
    xsec4to25(0) = xsec4to25orig(0);
    xsec4to25(1) = xsec4to25orig(1);
    xsec4to25(2) = xsec4to25orig(2);
  }








  TVectorD xsec4to35orig(3);
  xsec4to35orig.Zero();
  xsec4to35orig(0) = 1.5860;
  // xsec4to35orig(1) = 0.0430;
  xsec4to35orig(1) = 0.0993;
  xsec4to35orig(2) = 0.0830;
  Double_t Total4to35 = 1.882012;
  Double_t StatTotal4to35 = 0.0307*TMath::Sqrt(1.5);
  TVectorD xsec4to35(3);
  xsec4to35.Zero();
  if(         mode == 0 ){
    xsec4to35(1) = xsec4to35orig(1);
    xsec4to35(2) = xsec4to35orig(2);
    xsec4to35(0) = xsec4to35orig(0);
  } else if ( mode == 1) {
    xsec4to35(0) = xsec4to35orig(0)*(1. + PartialSys4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. + PartialSys4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. + PartialSys4to35[2]);
  } else if ( mode == 2) {
    xsec4to35(0) = xsec4to35orig(0)*(1. - PartialSys4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. - PartialSys4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. - PartialSys4to35[2]);
  } else if ( mode == 3) {
    xsec4to35(0) = xsec4to35orig(0)*(1. + PileUp4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. + PileUp4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. + PileUp4to35[2]);
  } else if ( mode == 4) {
    xsec4to35(0) = xsec4to35orig(0)*(1. - PileUp4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. - PileUp4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. - PileUp4to35[2]);
  } else if ( mode == 5) {
    xsec4to35(0) = xsec4to35orig(0)*(1. + eZdc4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. + eZdc4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. + eZdc4to35[2]);
  } else if ( mode == 6) {
    xsec4to35(0) = xsec4to35orig(0)*(1. - eZdc4to35[0]);
    xsec4to35(1) = xsec4to35orig(1)*(1. - eZdc4to35[1]);
    xsec4to35(2) = xsec4to35orig(2)*(1. - eZdc4to35[2]);
  } else {
    xsec4to35(1) = xsec4to35orig(1);
    xsec4to35(2) = xsec4to35orig(2);
    xsec4to35(0) = xsec4to35orig(0);
  }







  TVectorD xsec35to3orig(3);
  xsec35to3orig.Zero();
  xsec35to3orig(0) = 2.3150;
  // xsec35to3orig(1) = 0.1092;
  xsec35to3orig(1) = 0.1681;
  xsec35to3orig(2) = 0.1693;
  Double_t Total35to3 = 2.856231;
  Double_t StatTotal35to3 = 0.0307*TMath::Sqrt(1.5);
  TVectorD xsec35to3(3);
  xsec35to3.Zero();
  if( mode == 1 ){
    xsec35to3(1) = xsec35to3orig(1);
    xsec35to3(2) = xsec35to3orig(2);
    xsec35to3(0) = xsec35to3orig(0);
  } else if ( mode == 1) {
    xsec35to3(0) = xsec35to3orig(0)*(1. + PartialSys35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. + PartialSys35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. + PartialSys35to3[2]);
  } else if ( mode == 2) {
    xsec35to3(0) = xsec35to3orig(0)*(1. - PartialSys35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. - PartialSys35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. - PartialSys35to3[2]);
  } else if ( mode == 3) {
    xsec35to3(0) = xsec35to3orig(0)*(1. + PileUp35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. + PileUp35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. + PileUp35to3[2]);
  } else if ( mode == 4) {
    xsec35to3(0) = xsec35to3orig(0)*(1. - PileUp35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. - PileUp35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. - PileUp35to3[2]);
  } else if ( mode == 5) {
    xsec35to3(0) = xsec35to3orig(0)*(1. + eZdc35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. + eZdc35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. + eZdc35to3[2]);
  } else if ( mode == 6) {
    xsec35to3(0) = xsec35to3orig(0)*(1. - eZdc35to3[0]);
    xsec35to3(1) = xsec35to3orig(1)*(1. - eZdc35to3[1]);
    xsec35to3(2) = xsec35to3orig(2)*(1. - eZdc35to3[2]);
  } else {
    xsec35to3(1) = xsec35to3orig(1);
    xsec35to3(2) = xsec35to3orig(2);
    xsec35to3(0) = xsec35to3orig(0);
  }











  TVectorD xsec3to25orig(3);
  xsec3to25orig.Zero();
  xsec3to25orig(0) = 2.6570;
  // xsec3to25orig(1) = 0.2302;
  xsec3to25orig(1) = 0.2359;
  xsec3to25orig(2) = 0.2681;
  Double_t Total3to25 = 3.479975;
  Double_t StatTotal3to25 = 0.0307*TMath::Sqrt(1.5);
  // TRandom R;
  TVectorD xsec3to25(3);
  xsec3to25.Zero();
  if(         mode == 0 ){
    xsec3to25(1) = xsec3to25orig(1);
    xsec3to25(2) = xsec3to25orig(2);
    xsec3to25(0) = xsec3to25orig(0);
  } else if ( mode == 1) {
    xsec3to25(0) = xsec3to25orig(0)*(1. + PartialSys3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. + PartialSys3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. + PartialSys3to25[2]);
  } else if ( mode == 2) {
    xsec3to25(0) = xsec3to25orig(0)*(1. - PartialSys3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. - PartialSys3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. - PartialSys3to25[2]);
  } else if ( mode == 3) {
    xsec3to25(0) = xsec3to25orig(0)*(1. + PileUp3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. + PileUp3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. + PileUp3to25[2]);
  } else if ( mode == 4) {
    xsec3to25(0) = xsec3to25orig(0)*(1. - PileUp3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. - PileUp3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. - PileUp3to25[2]);
  } else if ( mode == 5) {
    xsec3to25(0) = xsec3to25orig(0)*(1. + eZdc3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. + eZdc3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. + eZdc3to25[2]);
  } else if ( mode == 6) {
    xsec3to25(0) = xsec3to25orig(0)*(1. - eZdc3to25[0]);
    xsec3to25(1) = xsec3to25orig(1)*(1. - eZdc3to25[1]);
    xsec3to25(2) = xsec3to25orig(2)*(1. - eZdc3to25[2]);
  } else {
    xsec3to25(1) = xsec3to25orig(1);
    xsec3to25(2) = xsec3to25orig(2);
    xsec3to25(0) = xsec3to25orig(0);
  }






  const TVectorD photonuclear4to25 = svd.Solve(xsec4to25,ok);
  const TVectorD photonuclear4to35 = svd0.Solve(xsec4to35,ok0);
  const TVectorD photonuclear35to3 = svd1.Solve(xsec35to3,ok1);
  const TVectorD photonuclear3to25 = svd2.Solve(xsec3to25,ok2);
  cout << "Was it fine? " << ok << ok0 << ok1 << ok2 << endl;
  photonuclear4to25.Print();
  photonuclear4to35.Print();
  photonuclear35to3.Print();
  photonuclear3to25.Print();















  //=======================
  // MIDRAPIDITY
  //-----------------------
  // TGraphAsymmErrors* Cent0N0N = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0n0n");
  // TGraphAsymmErrors* Cent0NXN = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_0nXn");
  // TGraphAsymmErrors* CentXNXN = (TGraphAsymmErrors*) file->Get("gXsection_Cent_Syst_XnXn");
  TVectorD xsec0575(3);
  xsec0575.Zero();
  xsec0575(0) = Cent0N0N->GetY()[0];
  xsec0575(1) = Cent0NXN->GetY()[0];
  xsec0575(2) = CentXNXN->GetY()[0];
  TVectorD xsec025(3);
  xsec025.Zero();
  xsec025(0) = Cent0N0N->GetY()[1];
  xsec025(1) = Cent0NXN->GetY()[1];
  xsec025(2) = CentXNXN->GetY()[1];
  TVectorD xsec000(3);
  xsec000.Zero();
  xsec000(0) = Cent0N0N->GetY()[2];
  xsec000(1) = Cent0NXN->GetY()[2];
  xsec000(2) = CentXNXN->GetY()[2];
  if( mode == 1 ){
    xsec0575(0) = xsec0575(0) + Cent0N0N->GetErrorY(0);
    xsec0575(1) = xsec0575(1) + Cent0NXN->GetErrorY(0);
    xsec0575(2) = xsec0575(2) + CentXNXN->GetErrorY(0);
    xsec025(0)  = xsec025(0)  + Cent0N0N->GetErrorY(1);
    xsec025(1)  = xsec025(1)  + Cent0NXN->GetErrorY(1);
    xsec025(2)  = xsec025(2)  + CentXNXN->GetErrorY(1);
    xsec000(0)  = xsec000(0)  + Cent0N0N->GetErrorY(2);
    xsec000(1)  = xsec000(1)  + Cent0NXN->GetErrorY(2);
    xsec000(2)  = xsec000(2)  + CentXNXN->GetErrorY(2);
  } else if (mode == 2){
    xsec0575(0) = xsec0575(0) - Cent0N0N->GetErrorY(0);
    xsec0575(1) = xsec0575(1) - Cent0NXN->GetErrorY(0);
    xsec0575(2) = xsec0575(2) - CentXNXN->GetErrorY(0);
    xsec025(0)  = xsec025(0)  - Cent0N0N->GetErrorY(1);
    xsec025(1)  = xsec025(1)  - Cent0NXN->GetErrorY(1);
    xsec025(2)  = xsec025(2)  - CentXNXN->GetErrorY(1);
    xsec000(0)  = xsec000(0)  - Cent0N0N->GetErrorY(2);
    xsec000(1)  = xsec000(1)  - Cent0NXN->GetErrorY(2);
    xsec000(2)  = xsec000(2)  - CentXNXN->GetErrorY(2);
  }
  // file->Close();
  const TVectorD photonuclear0575 = svd0_mid.Solve(xsec0575,ok0_mid);
  const TVectorD photonuclear025  = svd1_mid.Solve(xsec025, ok1_mid);
  const TVectorD photonuclear000  = svd2_mid.Solve(xsec000, ok2_mid);
  cout << "Was it fine? " <<  ok0_mid << ok1_mid << ok2_mid << endl;
  photonuclear0575.Print();
  photonuclear025.Print();
  photonuclear000.Print();











  // TFile* file2 = new TFile("xSection_Cent_2bins.root");
  // TGraphAsymmErrors* Cent0N0N_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0n0n");
  // TGraphAsymmErrors* Cent0NXN_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0nXn");
  // TGraphAsymmErrors* CentXNXN_2 = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_XnXn");

  TVectorD xsec05(3);
  xsec05.Zero();
  xsec05(0) = Cent0N0N_2->GetY()[0];
  xsec05(1) = Cent0NXN_2->GetY()[0];
  xsec05(2) = CentXNXN_2->GetY()[0];
  TVectorD xsec0(3);
  xsec0.Zero();
  xsec0(0) = Cent0N0N_2->GetY()[2];
  xsec0(1) = Cent0NXN_2->GetY()[2];
  xsec0(2) = CentXNXN_2->GetY()[2];

  if( mode == 1 ){
    xsec05(0) = xsec05(0) + Cent0N0N_2->GetErrorY(0);
    xsec05(1) = xsec05(1) + Cent0NXN_2->GetErrorY(0);
    xsec05(2) = xsec05(2) + CentXNXN_2->GetErrorY(0);
    xsec0(0)  = xsec0(0)  + Cent0N0N_2->GetErrorY(2);
    xsec0(1)  = xsec0(1)  + Cent0NXN_2->GetErrorY(2);
    xsec0(2)  = xsec0(2)  + CentXNXN_2->GetErrorY(2);
  } else if (mode == 2){
    xsec05(0) = xsec05(0) - Cent0N0N_2->GetErrorY(0);
    xsec05(1) = xsec05(1) - Cent0NXN_2->GetErrorY(0);
    xsec05(2) = xsec05(2) - CentXNXN_2->GetErrorY(0);
    xsec0(0)  = xsec0(0)  - Cent0N0N_2->GetErrorY(2);
    xsec0(1)  = xsec0(1)  - Cent0NXN_2->GetErrorY(2);
    xsec0(2)  = xsec0(2)  - CentXNXN_2->GetErrorY(2);
  }
  const TVectorD photonuclear05 = svd0_mid_2.Solve(xsec05,ok0_mid_2);
  const TVectorD photonuclear0  = svd1_mid_2.Solve(xsec0, ok1_mid_2);
  cout << "Was it fine? " <<  ok0_mid_2 << ok1_mid_2 << endl;
  photonuclear05.Print();
  photonuclear0.Print();














  Double_t Result = -999.;
  if( rap == 0){
    if (element == 0) {
      Result = photonuclear4to25(0);
    } else {
      Result = photonuclear4to25(1);
    }
  } else if (rap == 1) {
    if (element == 0) {
      Result = photonuclear4to35(0);
    } else {
      Result = photonuclear4to35(1);
    }
  } else if (rap == 2) {
    if (element == 0) {
      Result = photonuclear35to3(0);
    } else {
      Result = photonuclear35to3(1);
    }
  } else if (rap == 3) {
    if (element == 0) {
      Result = photonuclear3to25(0);
    } else {
      Result = photonuclear3to25(1);
    }
  } else if (rap == 4) {
    if (element == 0) {
      Result = photonuclear0575(0);
    } else {
      Result = photonuclear0575(1);
    }
  } else if (rap == 5) {
    if (element == 0) {
      Result = photonuclear025(0);
    } else {
      Result = photonuclear025(1);
    }
  } else if (rap == 6) {
    if (element == 0) {
      Result = photonuclear000(0);
    } else {
      Result = photonuclear000(1);
    }
  } else if (rap == 7) {
    if (element == 0) {
      Result = photonuclear05(0);
    } else {
      Result = photonuclear05(1);
    }
  } else if (rap == 8) {
    if (element == 0) {
      Result = photonuclear0(0);
    } else {
      Result = photonuclear0(1);
    }
  } else {
    Result = -999.;
  }
  return Result;
}
//_____________________________________________________________________________
void Computing(){
  Double_t HighSolutionx[6]  = {FittingPhotonuclear(0,1,0),FittingPhotonuclear(0,2,0),FittingPhotonuclear(0,3,0),FittingPhotonuclear(0,4,0),FittingPhotonuclear(0,5,0),FittingPhotonuclear(0,6,0)};
  Double_t LowSolutionx[6]   = {FittingPhotonuclear(1,1,0),FittingPhotonuclear(1,2,0),FittingPhotonuclear(1,3,0),FittingPhotonuclear(1,4,0),FittingPhotonuclear(1,5,0),FittingPhotonuclear(1,6,0)};

  Double_t HighSolution_corr   = 0.5*(HighSolutionx[0]-HighSolutionx[1]);
  Double_t HighSolution_pileup = 0.5*(HighSolutionx[2]-HighSolutionx[3]);
  Double_t HighSolution_ezdc   = 0.5*(HighSolutionx[4]-HighSolutionx[5]);
  Double_t LowSolution_corr    = 0.5*(LowSolutionx[0]-LowSolutionx[1]);
  Double_t LowSolution_pileup  = 0.5*(LowSolutionx[2]-LowSolutionx[3]);
  Double_t LowSolution_ezdc    = 0.5*(LowSolutionx[4]-LowSolutionx[5]);




  Double_t HighSolution0x[6]  = {FittingPhotonuclear(0,1,1),FittingPhotonuclear(0,2,1),FittingPhotonuclear(0,3,1),FittingPhotonuclear(0,4,1),FittingPhotonuclear(0,5,1),FittingPhotonuclear(0,6,1)};
  Double_t LowSolution0x[6]   = {FittingPhotonuclear(1,1,1),FittingPhotonuclear(1,2,1),FittingPhotonuclear(1,3,1),FittingPhotonuclear(1,4,1),FittingPhotonuclear(1,5,1),FittingPhotonuclear(1,6,1)};

  Double_t HighSolution0_corr   = 0.5*(HighSolution0x[0]-HighSolution0x[1]);
  Double_t HighSolution0_pileup = 0.5*(HighSolution0x[2]-HighSolution0x[3]);
  Double_t HighSolution0_ezdc   = 0.5*(HighSolution0x[4]-HighSolution0x[5]);
  Double_t LowSolution0_corr    = 0.5*(LowSolution0x[0]-LowSolution0x[1]);
  Double_t LowSolution0_pileup  = 0.5*(LowSolution0x[2]-LowSolution0x[3]);
  Double_t LowSolution0_ezdc    = 0.5*(LowSolution0x[4]-LowSolution0x[5]);



  Double_t HighSolution1x[6]  = {FittingPhotonuclear(0,1,2),FittingPhotonuclear(0,2,2),FittingPhotonuclear(0,3,2),FittingPhotonuclear(0,4,2),FittingPhotonuclear(0,5,2),FittingPhotonuclear(0,6,2)};
  Double_t LowSolution1x[6]   = {FittingPhotonuclear(1,1,2),FittingPhotonuclear(1,2,2),FittingPhotonuclear(1,3,2),FittingPhotonuclear(1,4,2),FittingPhotonuclear(1,5,2),FittingPhotonuclear(1,6,2)};

  Double_t HighSolution1_corr   = 0.5*(HighSolution1x[0]-HighSolution1x[1]);
  Double_t HighSolution1_pileup = 0.5*(HighSolution1x[2]-HighSolution1x[3]);
  Double_t HighSolution1_ezdc   = 0.5*(HighSolution1x[4]-HighSolution1x[5]);
  Double_t LowSolution1_corr    = 0.5*(LowSolution1x[0]-LowSolution1x[1]);
  Double_t LowSolution1_pileup  = 0.5*(LowSolution1x[2]-LowSolution1x[3]);
  Double_t LowSolution1_ezdc    = 0.5*(LowSolution1x[4]-LowSolution1x[5]);





  Double_t HighSolution2x[6]  = {FittingPhotonuclear(0,1,3),FittingPhotonuclear(0,2,3),FittingPhotonuclear(0,3,3),FittingPhotonuclear(0,4,3),FittingPhotonuclear(0,5,3),FittingPhotonuclear(0,6,3)};
  Double_t LowSolution2x[6]   = {FittingPhotonuclear(1,1,3),FittingPhotonuclear(1,2,3),FittingPhotonuclear(1,3,3),FittingPhotonuclear(1,4,3),FittingPhotonuclear(1,5,3),FittingPhotonuclear(1,6,3)};

  Double_t HighSolution2_corr   = 0.5*(HighSolution2x[0]-HighSolution2x[1]);
  Double_t HighSolution2_pileup = 0.5*(HighSolution2x[2]-HighSolution2x[3]);
  Double_t HighSolution2_ezdc   = 0.5*(HighSolution2x[4]-HighSolution2x[5]);
  Double_t LowSolution2_corr    = 0.5*(LowSolution2x[0]-LowSolution2x[1]);
  Double_t LowSolution2_pileup  = 0.5*(LowSolution2x[2]-LowSolution2x[3]);
  Double_t LowSolution2_ezdc    = 0.5*(LowSolution2x[4]-LowSolution2x[5]);











  Double_t HighSolution  = FittingPhotonuclear(0,0,0);
  Double_t LowSolution   = FittingPhotonuclear(1,0,0);



  Double_t HighSolution0  = FittingPhotonuclear(0,0,1);
  Double_t LowSolution0   = FittingPhotonuclear(1,0,1);




  Double_t HighSolution1  = FittingPhotonuclear(0,0,2);
  Double_t LowSolution1   = FittingPhotonuclear(1,0,2);




  Double_t HighSolution2  = FittingPhotonuclear(0,0,3);
  Double_t LowSolution2   = FittingPhotonuclear(1,0,3);

  cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << ", " << HighSolution_corr << ", " << HighSolution_pileup << ", " << HighSolution_ezdc << endl;
  cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << ", " << LowSolution_corr  << ", " << LowSolution_pileup  << ", " << LowSolution_ezdc  << endl;

  cout << "Element 0 (-4 < y < -3.5)= " << HighSolution0 << ", " << HighSolution0_corr << ", " << HighSolution0_pileup << ", " << HighSolution0_ezdc << endl;
  cout << "Element 1 (-4 < y < -3.5)= " << LowSolution0  << ", " << LowSolution0_corr  << ", " << LowSolution0_pileup  << ", " << LowSolution0_ezdc  << endl;

  cout << "Element 0 (-3.5 < y < -3)= " << HighSolution1 << ", " << HighSolution1_corr << ", " << HighSolution1_pileup << ", " << HighSolution1_ezdc << endl;
  cout << "Element 1 (-3.5 < y < -3)= " << LowSolution1  << ", " << LowSolution1_corr  << ", " << LowSolution1_pileup  << ", " << LowSolution1_ezdc  << endl;

  cout << "Element 0 (-3 < y < -2.5)= " << HighSolution2 << ", " << HighSolution2_corr << ", " << HighSolution2_pileup << ", " << HighSolution2_ezdc << endl;
  cout << "Element 1 (-3 < y < -2.5)= " << LowSolution2  << ", " << LowSolution2_corr  << ", " << LowSolution2_pileup  << ", " << LowSolution2_ezdc  << endl;


  // //==========================
  // // MIDRAPIDITY
  // //--------------------------
  Double_t HighSolution0575  = FittingPhotonuclear(0,0,4);
  Double_t HighSolution0575U = FittingPhotonuclear(0,1,4)-FittingPhotonuclear(0,0,4);
  Double_t HighSolution0575D = FittingPhotonuclear(0,2,4)-FittingPhotonuclear(0,0,4);
  Double_t LowSolution0575   = FittingPhotonuclear(1,0,4);
  Double_t LowSolution0575U  = FittingPhotonuclear(1,1,4)-FittingPhotonuclear(1,0,4);
  Double_t LowSolution0575D  = FittingPhotonuclear(1,2,4)-FittingPhotonuclear(1,0,4);

  Double_t HighSolution025  = FittingPhotonuclear(0,0,5);
  Double_t HighSolution025U = FittingPhotonuclear(0,1,5)-FittingPhotonuclear(0,0,5);
  Double_t HighSolution025D = FittingPhotonuclear(0,2,5)-FittingPhotonuclear(0,0,5);
  Double_t LowSolution025   = FittingPhotonuclear(1,0,5);
  Double_t LowSolution025U  = FittingPhotonuclear(1,1,5)-FittingPhotonuclear(1,0,5);
  Double_t LowSolution025D  = FittingPhotonuclear(1,2,5)-FittingPhotonuclear(1,0,5);

  Double_t HighSolution000  = FittingPhotonuclear(0,0,6);
  Double_t HighSolution000U = FittingPhotonuclear(0,1,6)-FittingPhotonuclear(0,0,6);
  Double_t HighSolution000D = FittingPhotonuclear(0,2,6)-FittingPhotonuclear(0,0,6);
  Double_t LowSolution000   = FittingPhotonuclear(1,0,6);
  Double_t LowSolution000U  = FittingPhotonuclear(1,1,6)-FittingPhotonuclear(1,0,6);
  Double_t LowSolution000D  = FittingPhotonuclear(1,2,6)-FittingPhotonuclear(1,0,6);

  cout << "Element 0 (y = 0.575)= " << HighSolution0575 << " + " << HighSolution0575U << " - " << HighSolution0575D << endl;
  cout << "Element 1 (y = 0.575)= " << LowSolution0575  << " + " << LowSolution0575U  << " - " << LowSolution0575D  << endl;

  cout << "Element 0 (y = 0.25)= "  << HighSolution025 << " + " << HighSolution025U << " - " << HighSolution025D << endl;
  cout << "Element 1 (y = 0.25)= "  << LowSolution025  << " + " << LowSolution025U  << " - " << LowSolution025D  << endl;

  cout << "Element 0 (y = 0.)= "    << HighSolution000 << " + " << HighSolution000U << " - " << HighSolution000D << endl;
  cout << "Element 1 (y = 0.)= "    << LowSolution000  << " + " << LowSolution000U  << " - " << LowSolution000D  << endl;

  Double_t HighSolution05  = FittingPhotonuclear(0,0,7);
  Double_t HighSolution05U = FittingPhotonuclear(0,1,7)-FittingPhotonuclear(0,0,7);
  Double_t HighSolution05D = FittingPhotonuclear(0,2,7)-FittingPhotonuclear(0,0,7);
  Double_t LowSolution05   = FittingPhotonuclear(1,0,7);
  Double_t LowSolution05U  = FittingPhotonuclear(1,1,7)-FittingPhotonuclear(1,0,7);
  Double_t LowSolution05D  = FittingPhotonuclear(1,2,7)-FittingPhotonuclear(1,0,7);

  Double_t HighSolution0X  = FittingPhotonuclear(0,0,8);
  Double_t HighSolution0XU = FittingPhotonuclear(0,1,8)-FittingPhotonuclear(0,0,8);
  Double_t HighSolution0XD = FittingPhotonuclear(0,2,8)-FittingPhotonuclear(0,0,8);
  Double_t LowSolution0X   = FittingPhotonuclear(1,0,8);
  Double_t LowSolution0XU  = FittingPhotonuclear(1,1,8)-FittingPhotonuclear(1,0,8);
  Double_t LowSolution0XD  = FittingPhotonuclear(1,2,8)-FittingPhotonuclear(1,0,8);

  cout << "Element 0 (y = 0.5)= " << HighSolution05 << " + " << HighSolution05U << " - " << HighSolution05D << endl;
  cout << "Element 1 (y = 0.5)= " << LowSolution05  << " + " << LowSolution05U  << " - " << LowSolution05D  << endl;

  cout << "Element 0 (y = 0.)= "  << HighSolution0X << " + " << HighSolution0XU << " - " << HighSolution0XD << endl;
  cout << "Element 1 (y = 0.)= "  << LowSolution0X  << " + " << LowSolution0XU  << " - " << LowSolution0XD  << endl;





  Double_t rapidity[12] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75, -0.5, 0., 0., 0.5 };
  Double_t Wgp[12];
  for (Int_t i = 0; i < 12; i++) {
    Wgp[i] = TMath::Sqrt(2.*2510*3.1*TMath::Exp(rapidity[i]));
  }
  Double_t y[12]    = {HighSolution,  HighSolution0,  HighSolution1,  HighSolution2,  LowSolution,  LowSolution0,  LowSolution1,  LowSolution2,  HighSolution05,  HighSolution0X,  LowSolution0X,  LowSolution05  };
  Double_t ey[12]   = {HighSolution_corr, HighSolution0_corr, HighSolution1_corr, HighSolution2_corr, LowSolution_corr, LowSolution0_corr, LowSolution1_corr, LowSolution2_corr, HighSolution05U, HighSolution0XU, LowSolution0XU, LowSolution05U };
  Double_t ey2[12]  = {HighSolution_pileup, HighSolution0_pileup, HighSolution1_pileup, HighSolution2_pileup, LowSolution_pileup, LowSolution0_pileup, LowSolution1_pileup, LowSolution2_pileup, HighSolution05U, HighSolution0XU, LowSolution0XU, LowSolution05U };
  Double_t ey3[12]  = {HighSolution_ezdc, HighSolution0_ezdc, HighSolution1_ezdc, HighSolution2_ezdc, LowSolution_ezdc, LowSolution0_ezdc, LowSolution1_ezdc, LowSolution2_ezdc, HighSolution05U, HighSolution0XU, LowSolution0XU, LowSolution05U };





  Double_t TotalEy[12];
  for (size_t i = 0; i < 12; i++) {
    TotalEy[i] = TMath::Sqrt( ey[i]*ey[i] + ey2[i]*ey2[i] + ey3[i]*ey3[i] );
  }












  for (Int_t i = 0; i < 8; i++) {
    cout << "============" << endl;
    cout << "Wgp["   << i << "] = " << Wgp[i] << endl;
    cout << "sigma[" << i << "] = " << y[i]   << endl;
    cout << "corr.["<< i << "]  = " << ey[i]  << endl;
    cout << "pileup["<< i << "] = " << ey2[i]  << endl;
    cout << "eZdc["<< i << "]   = " << ey3[i]  << endl;
    cout << "TOTAL["<< i << "]   = "<< TotalEy[i]  << endl;
  }

}
