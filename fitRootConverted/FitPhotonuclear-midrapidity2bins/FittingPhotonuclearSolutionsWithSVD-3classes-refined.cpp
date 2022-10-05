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
  Amatrix(0,0) = 157.789;
  Amatrix(0,1) = 1.07192;
  Amatrix(1,0) = 16.6849*0.5;
  Amatrix(1,1) = 1.70239*0.5;
  // Amatrix(2,0) = 16.6849*0.5;
  // Amatrix(2,1) = 1.70239*0.5;
  Amatrix(2,0) = 5.03916;
  Amatrix(2,1) = 0.949336;
  // -4 < y < -3.5
  Amatrix0(0,0) = 172.551;
  Amatrix0(0,1) = 0.200768;
  Amatrix0(1,0) = 16.6506*0.5;
  Amatrix0(1,1) = 0.421108*0.5;
  // Amatrix0(2,0) = 16.6506*0.5;
  // Amatrix0(2,1) = 0.421108*0.5;
  Amatrix0(2,0) = 5.02952;
  Amatrix0(2,1) = 0.275626;
  // -3.5 < y < -3.0
  Amatrix1(0,0) = 157.789;
  Amatrix1(0,1) = 1.07192;
  Amatrix1(1,0) = 16.6849*0.5;
  Amatrix1(1,1) = 1.70239*0.5;
  // Amatrix1(2,0) = 16.6849*0.5;
  // Amatrix1(2,1) = 1.70239*0.5;
  Amatrix1(2,0) = 5.03916;
  Amatrix1(2,1) = 0.949336;
  // -3.0 < y < -2.5
  Amatrix2(0,0) = 142.948;
  Amatrix2(0,1) = 3.66191;
  Amatrix2(1,0) = 16.7132*0.5;
  Amatrix2(1,1) = 4.28895*0.5;
  // Amatrix2(2,0) = 16.7132*0.5;
  // Amatrix2(2,1) = 4.28895*0.5;
  Amatrix2(2,0) = 5.04786;
  Amatrix2(2,1) = 2.04812;


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
  Amatrix0_mid_twobins(0,0) = 76.1357;
  Amatrix0_mid_twobins(0,1) = 47.3637;
  Amatrix0_mid_twobins(1,0) = 16.5274;
  Amatrix0_mid_twobins(1,1) = 15.494;
  Amatrix0_mid_twobins(2,0) = 5.06535;
  Amatrix0_mid_twobins(2,1) = 4.94607;
  // y = 0.
  Amatrix1_mid_twobins(0,0) = 61.537;
  Amatrix1_mid_twobins(0,1) = 61.537;
  Amatrix1_mid_twobins(1,0) = 16.1841;
  Amatrix1_mid_twobins(1,1) = 16.1841;
  Amatrix1_mid_twobins(2,0) = 5.03313;
  Amatrix1_mid_twobins(2,1) = 5.03313;
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
  // if( mode == 1 ){
  //   xsec4to25(0) = 2.2390 +TMath::Sqrt( 0.0307*0.0307 +  0.0855*0.0855*2.2390*2.2390 );
  //   xsec4to25(1) = 0.1667 +TMath::Sqrt( 0.0074*0.0074 +  0.1085*0.1085*0.1667*0.1667 );
  //   xsec4to25(2) = 0.1695 +TMath::Sqrt( 0.0088*0.0088 +  0.0924*0.0924*0.1695*0.1695 );
  // } else if (mode == 2){
  //   xsec4to25(0) = 2.2390 -TMath::Sqrt( 0.0307*0.0307 +  0.0855*0.0855*2.2390*2.2390 );
  //   xsec4to25(1) = 0.1667 -TMath::Sqrt( 0.0074*0.0074 +  0.1085*0.1085*0.1667*0.1667 );
  //   xsec4to25(2) = 0.1695 -TMath::Sqrt( 0.0088*0.0088 +  0.0924*0.0924*0.1695*0.1695 );
  //
  // }
  TRandom R;
  TVectorD xsec4to25(3);
  xsec4to25.Zero();
  xsec4to25(1) = gRandom->Gaus(xsec4to25orig(1), 0.0074);
  xsec4to25(2) = gRandom->Gaus(xsec4to25orig(2), 0.0088);
  xsec4to25(0) = gRandom->Gaus(Total4to25, StatTotal4to25) - 2.*xsec4to25(1) - xsec4to25(2);









  TVectorD xsec4to35orig(3);
  xsec4to35orig.Zero();
  xsec4to35orig(0) = 1.5860;
  // xsec4to35orig(1) = 0.0430;
  xsec4to35orig(1) = 0.0993;
  xsec4to35orig(2) = 0.0830;
  // if( mode == 1 ){
  //   xsec4to35(0) = 1.5860 +TMath::Sqrt( 0.0307*0.0307 +  0.0690*0.0690*1.5860*1.5860 );
  //   // xsec4to35(1) = 0.0430 +TMath::Sqrt( 0.0489*0.0489 +  0.3872*0.3872*0.0430*0.0430 );
  //   xsec4to35(1) = 0.0993 +TMath::Sqrt( 0.0084*0.0084 +  0.1092*0.1092*0.0993*0.0993 );
  //   xsec4to35(2) = 0.0830 +TMath::Sqrt( 0.0119*0.0119 +  0.0761*0.0761*0.0830*0.0830 );
  // } else if (mode == 2){
  //   xsec4to35(0) = 1.5860 -TMath::Sqrt( 0.0307*0.0307 +  0.0690*0.0690*1.5860*1.5860 );
  //   // xsec4to35(1) = 0.0430 -TMath::Sqrt( 0.0489*0.0489 +  0.3872*0.3872*0.0430*0.0430 );
  //   xsec4to35(1) = 0.0993 -TMath::Sqrt( 0.0084*0.0084 +  0.1092*0.1092*0.0993*0.0993 );
  //   xsec4to35(2) = 0.0830 -TMath::Sqrt( 0.0119*0.0119 +  0.0761*0.0761*0.0830*0.0830 );
  //
  // }
  Double_t Total4to35 = 1.882012;
  Double_t StatTotal4to35 = 0.0307*TMath::Sqrt(1.5);
  // TRandom R;
  TVectorD xsec4to35(3);
  xsec4to35.Zero();
  // xsec4to35(0) = gRandom->Gaus(Total4to35, StatTotal4to35) - gRandom->Gaus(xsec4to35orig(1), 0.0489) - gRandom->Gaus(xsec4to35orig(2), 0.0084) - gRandom->Gaus(xsec4to35orig(3), 0.0119);
  // xsec4to35(1) = gRandom->Gaus(xsec4to35orig(1), 0.0035);
  xsec4to35(1) = gRandom->Gaus(xsec4to35orig(1), 0.0078);
  xsec4to35(2) = gRandom->Gaus(xsec4to35orig(2), 0.0119);
  xsec4to35(0) = gRandom->Gaus(Total4to35, StatTotal4to35) - 2.*  xsec4to35(1) - xsec4to35(2);







  TVectorD xsec35to3orig(3);
  xsec35to3orig.Zero();
  xsec35to3orig(0) = 2.3150;
  // xsec35to3orig(1) = 0.1092;
  xsec35to3orig(1) = 0.1681;
  xsec35to3orig(2) = 0.1693;
  // if( mode == 1 ){
  //   xsec35to3(0) = 2.3150 +TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.3150*2.3150 );
  //   // xsec35to3(1) = 0.1092 +TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.1092*0.1092 );
  //   xsec35to3(1) = 0.1681 +TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.1681*0.1681 );
  //   xsec35to3(2) = 0.1693 +TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.1693*0.1693 );
  // } else if (mode == 2){
  //   xsec35to3(0) = 2.3150 -TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.3150*2.3150 );
  //   // xsec35to3(1) = 0.1092 -TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.1092*0.1092 );
  //   xsec35to3(1) = 0.1681 -TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.1681*0.1681 );
  //   xsec35to3(2) = 0.1693 -TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.1693*0.1693 );
  //
  // }
  Double_t Total35to3 = 2.856231;
  Double_t StatTotal35to3 = 0.0307*TMath::Sqrt(1.5);
  TVectorD xsec35to3(3);
  xsec35to3.Zero();
  // xsec35to3(1) = gRandom->Gaus(xsec35to3orig(1), 0.0053);
  xsec35to3(1) = gRandom->Gaus(xsec35to3orig(1), 0.0098);
  xsec35to3(2) = gRandom->Gaus(xsec35to3orig(2), 0.0120);
  xsec35to3(0) = gRandom->Gaus(Total35to3, StatTotal35to3) -2.* xsec35to3(1) - xsec35to3(2);










  TVectorD xsec3to25orig(3);
  xsec3to25orig.Zero();
  xsec3to25orig(0) = 2.6570;
  // xsec3to25orig(1) = 0.2302;
  xsec3to25orig(1) = 0.2359;
  xsec3to25orig(2) = 0.2681;
  // if( mode == 1 ){
  //   xsec3to25(0) = 2.6570 +TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.6570*2.6570 );
  //   // xsec3to25(1) = 0.2302 +TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.2302*0.2302 );
  //   xsec3to25(1) = 0.2359 +TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.2359*0.2359 );
  //   xsec3to25(2) = 0.2681 +TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.2681*0.2681 );
  // } else if (mode == 2){
  //   xsec3to25(0) = 2.6570 -TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.6570*2.6570 );
  //   // xsec3to25(1) = 0.2302 -TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.2302*0.2302 );
  //   xsec3to25(1) = 0.2359 -TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.2359*0.2359 );
  //   xsec3to25(2) = 0.2681 -TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.2681*0.2681 );
  //
  // }
  Double_t Total3to25 = 3.479975;
  Double_t StatTotal3to25 = 0.0307*TMath::Sqrt(1.5);
  // TRandom R;
  TVectorD xsec3to25(3);
  xsec3to25.Zero();
  // xsec3to25(1) = gRandom->Gaus(xsec3to25orig(1), 0.0126);
  xsec3to25(1) = gRandom->Gaus(xsec3to25orig(1), 0.0208);
  xsec3to25(2) = gRandom->Gaus(xsec3to25orig(2), 0.0249);
  xsec3to25(0) = gRandom->Gaus(Total3to25, StatTotal3to25) - 2.*xsec3to25(1) - xsec3to25(2);






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
  const TVectorD photonuclear05 = svd0_mid.Solve(xsec05,ok0_mid_2);
  const TVectorD photonuclear0  = svd2_mid.Solve(xsec0, ok1_mid_2);
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
  Double_t HighSolutionx  = -1.;
  Double_t LowSolutionx   = -1.;
  TH1F* histoHigh4to25 = new TH1F("histoHigh4to25", "histoHigh4to25", 100000, 0, 10);
  TH1F* histoLow4to25  = new TH1F("histoLow4to25",  "histoLow4to25",  100000, 0, 10);
  for (Int_t i = 0; i < 10000; i++) {
    HighSolutionx  = FittingPhotonuclear(0,0,0);
    LowSolutionx   = FittingPhotonuclear(1,0,0);
    histoHigh4to25->Fill(HighSolutionx);
    histoLow4to25->Fill(LowSolutionx);
  }
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoHigh4to25);
  histoHigh4to25->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoHigh4to25->GetYaxis()->SetTitle("Counts [a.u.]");
  histoHigh4to25->GetXaxis()->SetRangeUser(0.,10);
  histoHigh4to25->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"High solution");
  latex10->DrawLatex(0.2,0.70,"-4 < y < -2.5");
  TF1* gHigh4to25 = new TF1("gHigh4to25", "gaus", 0, 2);
  gHigh4to25->SetParameter(0,1);
  gHigh4to25->SetParameter(1,0.3);
  gHigh4to25->SetParameter(2,0.3);
  histoHigh4to25->Fit(gHigh4to25);
  Double_t HighSolutionU = gHigh4to25->GetParameter(2);
  Double_t HighSolutionD = gHigh4to25->GetParameter(1);
  // gPad->SaveAs(Form("SignalExtractionCoarse/Fitting/residual-refolded-data-%d.pdf", Iterations), "recreate");
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLow4to25);
  histoLow4to25->Rebin(10);
  histoLow4to25->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoLow4to25->GetYaxis()->SetTitle("Counts [a.u.]");
  histoLow4to25->GetXaxis()->SetRangeUser(0.,10);
  histoLow4to25->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Low solution");
  latex10->DrawLatex(0.2,0.70,"-4 < y < -2.5");
  TF1* gLow4to25 = new TF1("gLow4to25", "gaus", 0, 2);
  gLow4to25->SetParameter(0,1);
  gLow4to25->SetParameter(1,0.3);
  gLow4to25->SetParameter(2,0.3);
  histoLow4to25->Fit(gLow4to25);
  Double_t LowSolutionU = gLow4to25->GetParameter(2);
  Double_t LowSolutionD = gLow4to25->GetParameter(1);










  Double_t HighSolution0x  = -1.;
  Double_t LowSolution0x   = -1.;
  TH1F* histoHigh4to35 = new TH1F("histoHigh4to35", "histoHigh4to35", 100000, 0, 10);
  TH1F* histoLow4to35  = new TH1F("histoLow4to35",  "histoLow4to35",  100000, 0, 10);
  for (Int_t i = 0; i < 10000; i++) {
    HighSolution0x  = FittingPhotonuclear(0,0,1);
    LowSolution0x   = FittingPhotonuclear(1,0,1);
    histoHigh4to35->Fill(HighSolution0x);
    histoLow4to35->Fill(LowSolution0x);
  }
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoHigh4to35);
  histoHigh4to35->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoHigh4to35->GetYaxis()->SetTitle("Counts [a.u.]");
  histoHigh4to35->GetXaxis()->SetRangeUser(0.,10);
  histoHigh4to35->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"High solution");
  latex10->DrawLatex(0.2,0.70,"-4 < y < -3.5");
  TF1* gHigh4to35 = new TF1("gHigh4to35", "gaus", 0, 2);
  gHigh4to35->SetParameter(0,1);
  gHigh4to35->SetParameter(1,0.3);
  gHigh4to35->SetParameter(2,0.3);
  histoHigh4to35->Fit(gHigh4to35);
  Double_t HighSolutionU0 = gHigh4to35->GetParameter(2);
  Double_t HighSolutionD0 = gHigh4to35->GetParameter(1);
  // gPad->SaveAs(Form("SignalExtractionCoarse/Fitting/residual-refolded-data-%d.pdf", Iterations), "recreate");
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLow4to35);
  histoLow4to35->Rebin(10);
  histoLow4to35->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoLow4to35->GetYaxis()->SetTitle("Counts [a.u.]");
  histoLow4to35->GetXaxis()->SetRangeUser(0.,10);
  histoLow4to35->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Low solution");
  latex10->DrawLatex(0.2,0.70,"-4 < y < -3.5");
  TF1* gLow4to35 = new TF1("gLow4to35", "gaus", 0, 2);
  gLow4to35->SetParameter(0,1);
  gLow4to35->SetParameter(1,0.3);
  gLow4to35->SetParameter(2,0.3);
  histoLow4to35->Fit(gLow4to35);
  Double_t LowSolutionU0 = gLow4to35->GetParameter(2);
  Double_t LowSolutionD0 = gLow4to35->GetParameter(1);









  Double_t HighSolution1x  = -1.;
  Double_t LowSolution1x   = -1.;
  TH1F* histoHigh35to3 = new TH1F("histoHigh35to3", "histoHigh35to3", 100000, 0, 10);
  TH1F* histoLow35to3  = new TH1F("histoLow35to3",  "histoLow35to3",  100000, 0, 10);
  for (Int_t i = 0; i < 10000; i++) {
    HighSolution1x  = FittingPhotonuclear(0,0,2);
    LowSolution1x   = FittingPhotonuclear(1,0,2);
    histoHigh35to3->Fill(HighSolution1x);
    histoLow35to3->Fill(LowSolution1x);
  }
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoHigh35to3);
  histoHigh35to3->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoHigh35to3->GetYaxis()->SetTitle("Counts [a.u.]");
  histoHigh35to3->GetXaxis()->SetRangeUser(0.,10);
  histoHigh35to3->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"High solution");
  latex10->DrawLatex(0.2,0.70,"-3.5 < y < -3");
  TF1* gHigh35to3 = new TF1("gHigh35to3", "gaus", 0, 2);
  gHigh35to3->SetParameter(0,1);
  gHigh35to3->SetParameter(1,0.3);
  gHigh35to3->SetParameter(2,0.3);
  histoHigh35to3->Fit(gHigh35to3);
  Double_t HighSolutionU1 = gHigh35to3->GetParameter(2);
  Double_t HighSolutionD1 = gHigh35to3->GetParameter(1);
  // gPad->SaveAs(Form("SignalExtractionCoarse/Fitting/residual-refolded-data-%d.pdf", Iterations), "recreate");
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLow35to3);
  histoLow35to3->Rebin(10);
  histoLow35to3->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoLow35to3->GetYaxis()->SetTitle("Counts [a.u.]");
  histoLow35to3->GetXaxis()->SetRangeUser(0.,10);
  histoLow35to3->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Low solution");
  latex10->DrawLatex(0.2,0.70,"-3.5 < y < -3");
  TF1* gLow35to3 = new TF1("gLow35to3", "gaus", 0, 2);
  gLow35to3->SetParameter(0,1);
  gLow35to3->SetParameter(1,0.3);
  gLow35to3->SetParameter(2,0.3);
  histoLow35to3->Fit(gLow35to3);
  Double_t LowSolutionU1 = gLow35to3->GetParameter(2);
  Double_t LowSolutionD1 = gLow35to3->GetParameter(1);








  Double_t HighSolution2x  = -1.;
  Double_t LowSolution2x   = -1.;
  TH1F* histoHigh3to25 = new TH1F("histoHigh3to25", "histoHigh3to25", 100000, 0, 10);
  TH1F* histoLow3to25  = new TH1F("histoLow3to25",  "histoLow3to25",  100000, 0, 10);
  for (Int_t i = 0; i < 10000; i++) {
    HighSolution2x  = FittingPhotonuclear(0,0,3);
    LowSolution2x   = FittingPhotonuclear(1,0,3);
    histoHigh3to25->Fill(HighSolution2x);
    histoLow3to25->Fill(LowSolution2x);
  }
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoHigh3to25);
  histoHigh3to25->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoHigh3to25->GetYaxis()->SetTitle("Counts [a.u.]");
  histoHigh3to25->GetXaxis()->SetRangeUser(0.,10);
  histoHigh3to25->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"High solution");
  latex10->DrawLatex(0.2,0.70,"-3 < y < -2.5");
  TF1* gHigh3to25 = new TF1("gHigh3to25", "gaus", 0, 2);
  gHigh3to25->SetParameter(0,1);
  gHigh3to25->SetParameter(1,0.3);
  gHigh3to25->SetParameter(2,0.3);
  histoHigh3to25->Fit(gHigh3to25);
  Double_t HighSolutionU2 = gHigh3to25->GetParameter(2);
  Double_t HighSolutionD2 = gHigh3to25->GetParameter(1);
  // gPad->SaveAs(Form("SignalExtractionCoarse/Fitting/residual-refolded-data-%d.pdf", Iterations), "recreate");
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLow3to25);
  histoLow3to25->Rebin(10);
  histoLow3to25->GetXaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histoLow3to25->GetYaxis()->SetTitle("Counts [a.u.]");
  histoLow3to25->GetXaxis()->SetRangeUser(0.,10);
  histoLow3to25->Draw("");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Low solution");
  latex10->DrawLatex(0.2,0.70,"-3 < y < -2.5");
  TF1* gLow3to25 = new TF1("gLow3to25", "gaus", 0, 2);
  gLow3to25->SetParameter(0,1);
  gLow3to25->SetParameter(1,0.3);
  gLow3to25->SetParameter(2,0.3);
  histoLow3to25->Fit(gLow3to25);
  Double_t LowSolutionU2 = gLow3to25->GetParameter(2);
  Double_t LowSolutionD2 = gLow3to25->GetParameter(1);














  Double_t HighSolution  = FittingPhotonuclear(0,0,0);
  Double_t LowSolution   = FittingPhotonuclear(1,0,0);



  Double_t HighSolution0  = FittingPhotonuclear(0,0,1);
  Double_t LowSolution0   = FittingPhotonuclear(1,0,1);




  Double_t HighSolution1  = FittingPhotonuclear(0,0,2);
  Double_t LowSolution1   = FittingPhotonuclear(1,0,2);




  Double_t HighSolution2  = FittingPhotonuclear(0,0,3);
  Double_t LowSolution2   = FittingPhotonuclear(1,0,3);

  cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << " + " << HighSolutionU << " - " << HighSolutionD << endl;
  cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << " + " << LowSolutionU  << " - " << LowSolutionD  << endl;

  cout << "Element 0 (-4 < y < -3.5)= " << HighSolution0 << " + " << HighSolutionU0 << " - " << HighSolutionD0 << endl;
  cout << "Element 1 (-4 < y < -3.5)= " << LowSolution0  << " + " << LowSolutionU0  << " - " << LowSolutionD0  << endl;

  cout << "Element 0 (-3.5 < y < -3)= " << HighSolution1 << " + " << HighSolutionU1 << " - " << HighSolutionD1 << endl;
  cout << "Element 1 (-3.5 < y < -3)= " << LowSolution1  << " + " << LowSolutionU1  << " - " << LowSolutionD1  << endl;

  cout << "Element 0 (-3 < y < -2.5)= " << HighSolution2 << " + " << HighSolutionU2 << " - " << HighSolutionD2 << endl;
  cout << "Element 1 (-3 < y < -2.5)= " << LowSolution2  << " + " << LowSolutionU2  << " - " << LowSolutionD2  << endl;


  //==========================
  // MIDRAPIDITY
  //--------------------------
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





  TH1F* histo = new TH1F("histo", "histo", 1000, -0.5, 999.5);
  // Double_t rapidity[14] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75, -0.575, -0.25, 0., 0., 0.25, 0.575};
  // Double_t Wgp[14];
  // for (Int_t i = 0; i < 14; i++) {
  //   Wgp[i] = TMath::Sqrt(2.*2510*3.1*TMath::Exp(rapidity[i]));
  // }
  // Double_t y[14]   = {HighSolution,  HighSolution0,  HighSolution1,  HighSolution2,  LowSolution,  LowSolution0,  LowSolution1,  LowSolution2,  HighSolution0575,  HighSolution025,  HighSolution000,  LowSolution000,  LowSolution025,  LowSolution0575  };
  // Double_t ey[14]  = {HighSolutionU, HighSolutionU0, HighSolutionU1, HighSolutionU2, LowSolutionU, LowSolutionU0, LowSolutionU1, LowSolutionU2, HighSolution0575U, HighSolution025U, HighSolution000U, LowSolution000U, LowSolution025U, LowSolution0575U };
  Double_t rapidity[12] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75, -0.5, 0., 0., 0.5 };
  Double_t Wgp[12];
  for (Int_t i = 0; i < 12; i++) {
    Wgp[i] = TMath::Sqrt(2.*2510*3.1*TMath::Exp(rapidity[i]));
  }
  Double_t y[12]   = {HighSolution,  HighSolution0,  HighSolution1,  HighSolution2,  LowSolution,  LowSolution0,  LowSolution1,  LowSolution2,  HighSolution05,  HighSolution0X,  LowSolution0X,  LowSolution05  };
  Double_t ey[12]  = {HighSolutionU, HighSolutionU0, HighSolutionU1, HighSolutionU2, LowSolutionU, LowSolutionU0, LowSolutionU1, LowSolutionU2, HighSolution05U, HighSolution0XU, LowSolution0XU, LowSolution05U };
  for (Int_t i = 0; i < 12; i++) {
    Int_t bin = histo->GetXaxis()->FindBin(Wgp[i]);
    cout << "bin = " << bin << endl;
    histo->SetBinContent( histo->GetXaxis()->FindBin(Wgp[i]),  y[i]);
    histo->SetBinError(   histo->GetXaxis()->FindBin(Wgp[i]), ey[i]);
  }
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTopMargin(0.14);
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  // gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
  histo->SetTitle("");
  histo->GetXaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->GetXaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.045);
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetYaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);

  histo->GetXaxis()->SetTitle("W_{#gammaPb} [GeV]");
  histo->GetYaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histo->GetYaxis()->SetRangeUser(0.002,0.2);
  // histo->GetYaxis()->SetRangeUser(-0.1,0.2);
  histo->GetXaxis()->SetRangeUser(10.,1000.);
  histo->SetLineWidth(5);
  histo->SetLineColor(2);
  histo->Draw("same");
  Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  Double_t bxmin = TMath::Power((mjpsi/1000.),2.);
  Double_t bxmax = TMath::Power((mjpsi/10.),2.);
  TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  fbx->SetParameter(0, mjpsi);

  TGaxis *axis = new TGaxis(1000., 0.2, 10., 0.2, "fbx", 510, "+G");
  axis->SetTextFont(42);
  axis->SetLabelFont(42);
  Double_t siz = 0.045;
  axis->SetTitleSize(siz); axis->SetLabelSize(siz);
  axis->SetLabelOffset(-0.035);
  axis->SetTitleOffset();
  axis->SetTitle("Bjorken-#it{x}");
  axis->Draw("same");

  TLatex *bxtit = new TLatex();
  bxtit->SetTextFont(42);
  bxtit->SetTextSize(siz);
  bxtit->SetTextAlign(31);
  bxtit->DrawLatex(1000., 0.2+450, "Bjorken-#it{x}");
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.7,"This thesis");
  // latex5->DrawLatex(0.31,0.8,"Coherent J/#psi");
  latex5->DrawLatex(0.31,0.78,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.7,"This thesis");
  latex5->DrawLatex(0.31,0.62,"Coherent J/#psi");

  gPad->SaveAs("solutions_mid.pdf", "recreate");














  // ==================
  // Suppression factor
  // ==================
  TH1F* SuppressionFactor = new TH1F("SuppressionFactor", "SuppressionFactor", 1000, -0.5, 999.5);
  // Double_t rapidity[8] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75};
  // // Double_t rapidity[8] = {3.25, 3.75, 3.25, 2.75, -3.25, -3.75, -3.25, -2.75};
  // Double_t Wgp[8];
  // for (Int_t i = 0; i < 8; i++) {
  //   Wgp[i] = TMath::Sqrt(2.*2760*3.1*TMath::Exp(rapidity[i]));
  // }
  // Double_t y[8]   = {HighSolution, HighSolution0, HighSolution1, HighSolution2, LowSolution, LowSolution0, LowSolution1, LowSolution2};
  // Double_t ey[8]  = {HighSolutionU, HighSolutionU0, HighSolutionU1, HighSolutionU2, LowSolutionU, LowSolutionU0, LowSolutionU1, LowSolutionU2};
  // ========================================
  // Wgp(1) (GeV) k*(dn_1/dk)  sigma_2(gam+A)
  // ----------------------------------------
  // IA
  // 2.1082E+01   1.9869E+02   1.5119E-02
  // 2.7060E+01   1.8317E+02   1.9706E-02
  // 3.4050E+01   1.6889E+02   2.3864E-02
  // 5.1928E+02   1.0599E+01   1.4670E-01
  // 6.6677E+02   4.0781E+00   1.7259E-01
  // 8.5615E+02   1.0018E+00   2.0305E-01

  // 1.9176E+01   2.0459E+02   1.3259E-02
  // 2.4610E+01   1.8907E+02   1.7995E-02
  // 3.1591E+01   1.7355E+02   2.2492E-02
  // 4.9396E+02   1.2312E+01   1.4200E-01
  // 6.3425E+02   5.0905E+00   1.6707E-01
  // 8.1439E+02   1.4145E+00   1.9655E-01


  // 97.1536
  // 124.748
  // 124.748
  // 160.179

  // STARlight
  // 9.9235E+01   1.0240E+02   4.9860E-02
  // 1.1645E+02   9.2484E+01   5.5380E-02
  // 1.3196E+02   8.4752E+01   6.0102E-02
  // 1.3196E+02   8.4752E+01   6.0102E-02
  // 1.4878E+02   7.7358E+01   6.5006E-02
  // 1.7547E+02   6.7252E+01   7.2396E-02

  // 9.3456E+01   1.0613E+02   4.7930E-02
  // 1.1022E+02   9.5890E+01   5.3419E-02
  // 1.2489E+02   8.8152E+01   5.7978E-02
  // 1.2489E+02   8.8152E+01   5.7978E-02
  // 1.4152E+02   8.0435E+01   6.2917E-02
  // 1.6608E+02   7.0611E+01   6.9845E-02



  // IA
  // 9.7270E+01   1.0364E+02   4.9208E-02
  // 1.2489E+02   8.8152E+01   5.7978E-02
  // 1.2489E+02   8.8152E+01   5.7978E-02
  // 1.6037E+02   7.2753E+01   6.8268E-02

  // STARlight
  // 9.7270E+01   1.0364E+02   3.2432E-02
  // 1.2489E+02   8.8152E+01   3.6961E-02
  // 1.2489E+02   8.8152E+01   3.6961E-02
  // 1.6037E+02   7.2753E+01   4.1996E-02

  Double_t IA[12]  = { 1.7995E-02, 1.3259E-02, 1.7995E-02, 2.2492E-02, 1.6707E-01, 1.9655E-01, 1.6707E-01, 1.4200E-01,             4.9208E-02, 5.7978E-02, 5.7978E-02, 6.8268E-02 };
  Double_t SuppressionFactorValue[12];
  Double_t SuppressionFactorError[12];
  for (Int_t i = 0; i < 14; i++) {
    Int_t bin = histo->GetXaxis()->FindBin(Wgp[i]);
    cout << "bin = " << bin << endl;
    SuppressionFactor->SetBinContent( SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  TMath::Sqrt(y[i]/IA[i]));
    // SuppressionFactor->SetBinError(   SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  TMath::Sqrt(ey[i]/IA[i]));
    SuppressionFactor->SetBinError(   SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  ey[i]/(2.*TMath::Sqrt(y[i]*IA[i])) );
    SuppressionFactorValue[i] = TMath::Sqrt(y[i]/IA[i]);
    SuppressionFactorError[i] = ey[i]/(2.*TMath::Sqrt(y[i]*IA[i]));
  }
  new TCanvas;
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
  SuppressionFactor->GetYaxis()->SetRangeUser(-0.3,1.5);
  SuppressionFactor->GetXaxis()->SetRangeUser(10.,1000.);
  SuppressionFactor->SetLineWidth(5);
  SuppressionFactor->SetLineColor(2);
  SuppressionFactor->Draw("same");
  // Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  // Double_t bxmin = TMath::Power((mjpsi/1000.),2.);
  // Double_t bxmax = TMath::Power((mjpsi/10.),2.);
  // TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  // fbx->SetParameter(0, mjpsi);

  TGaxis *axis2 = new TGaxis(1000., 1.5, 10., 1.5, "fbx", 510, "+G");
  axis2->SetTextFont(42);
  axis2->SetLabelFont(42);
  // Double_t siz = 0.045;
  axis2->SetTitleSize(siz); axis->SetLabelSize(siz);
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
  latex6->DrawLatex(0.4,0.7,"This thesis");
  latex6->DrawLatex(0.4,0.62,"Coherent J/#psi");

  gPad->SaveAs("SuppressionFactor_mid.pdf", "recreate");



  for (Int_t i = 0; i < 12; i++) {
    cout << "============" << endl;
    cout << "Wgp["   << i << "] = " << Wgp[i] << endl;
    cout << "sigma[" << i << "] = " << y[i]   << endl;
    cout << "Uncert["<< i << "] = " << ey[i]  << endl;
    cout << "S["     << i << "] = " << SuppressionFactorValue[i]   << endl;
    cout << "eS["    << i << "] = " << SuppressionFactorError[i]  << endl;
  }

}
