#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <vector>



//_____________________________________________________________________________
/* - Histograms to be used for the fit.
 * - What happens is that we will interpolate the many points together...
 * -
 */
TH1F* fTwoGammaToMuMedium;
TH1F* fTwoGammaToMuMediumOrig;

TH1F* fTwoGammaToMuMediumLowerSide;
TH1F* fTwoGammaToMuMediumOrigLowerSide;

TH1F* fTwoGammaToMuMediumHigherSide;
TH1F* fTwoGammaToMuMediumOrigHigherSide;

TH1F* Data0N0Nlower;
TH1F* Data0N0Nhigher;

TH1F* Data0NXNlower;
TH1F* Data0NXNhigher;

TH1F* DataXN0Nlower;
TH1F* DataXN0Nhigher;

TH1F* DataXNXNlower;
TH1F* DataXNXNhigher;





TH1F* fTwoGammaToMuMediumC;
TH1F* fTwoGammaToMuMediumOrigC;

TH1F* fTwoGammaToMuMediumLowerSideC;
TH1F* fTwoGammaToMuMediumOrigLowerSideC;

TH1F* fTwoGammaToMuMediumHigherSideC;
TH1F* fTwoGammaToMuMediumOrigHigherSideC;

TH1F* Data0N0NlowerC;
TH1F* Data0N0NhigherC;

TH1F* Data0NXNlowerC;
TH1F* Data0NXNhigherC;

TH1F* DataXN0NlowerC;
TH1F* DataXN0NhigherC;

TH1F* DataXNXNlowerC;
TH1F* DataXNXNhigherC;


//_____________________________________________________________________________
/* - Sidebands templates.
 * -
 */
void Sidebands(){
  TFile* fileList = new TFile("AnalysisResultsLHC18qr_19032020.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  // TFile* fileMC   = new TFile("MCtrainResults/2019-09-17/kTwoGammaToMuMedium/AnalysisResults.root");
  TFile* fileMC   = new TFile("AnalysisResultsSidebandsMC.root");
  TDirectory* dirMC = fileMC->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir->GetObject("ADcheck_", listings);
  TList* listingsMC;
  dirMC->GetObject("MyOutputContainer", listingsMC);
  /* - We now do the same as before to ascertain if the TList was there and
   * - to try to retrieve the plots. Result:
   *   listings->ls()
   *     OBJ: TList	  MyOutputContainer	          Doubly linked list          : 0
   *     OBJ: TH1F	  fNumberMuonsH	              fNumberMuonsH               : 0 at: 0x5a145f0
   *     OBJ: TH1F	  fCounterH	                  fCounterH                   : 0 at: 0x5a3b570
   *     OBJ: TH1F	  fEtaMuonH	                  fEtaMuonH                   : 0 at: 0x5a3ba80
   *     OBJ: TH1F	  fRAbsMuonH	                fRAbsMuonH                  : 0 at: 0x5a3c0c0
   *     OBJ: TH1F	  fInvariantMassDistributionH	fInvariantMassDistributionH : 0 at: 0x5a3c720
   */
  fTwoGammaToMuMediumOrig = (TH1F*)listingsMC->FindObject("fTemplatePtDistributionH");
  fTwoGammaToMuMediumOrig -> Rebin(5);

  fTwoGammaToMuMediumOrigLowerSide = (TH1F*)listingsMC->FindObject("fTemplatePtDistributionHLowerSide");
  fTwoGammaToMuMediumOrigLowerSide -> Rebin(5);

  fTwoGammaToMuMediumOrigHigherSide = (TH1F*)listingsMC->FindObject("fTemplatePtDistributionHHigherSide");
  fTwoGammaToMuMediumOrigHigherSide -> Rebin(5);



  Float_t PtBins[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
                          0.200, 0.225, 0.250, 0.275, 0.350, 0.425, 0.500, 0.575,
                          0.650, 0.725,
                          0.800, 0.875, 0.950, 1.100, 1.250, 1.400, 1.600, 1.800,
                          2.000, 2.500, 3.000, 3.500, 4.000, 5.000
                        };
  Int_t   PtBinNumber = sizeof(PtBins)/sizeof(Float_t) - 1; // or just = 9
  fTwoGammaToMuMedium           = new TH1F("fTwoGammaToMuMedium",           "fTwoGammaToMuMedium",           PtBinNumber, PtBins );
  fTwoGammaToMuMediumLowerSide  = new TH1F("fTwoGammaToMuMediumLowerSide",  "fTwoGammaToMuMediumLowerSide",  PtBinNumber, PtBins );
  fTwoGammaToMuMediumHigherSide = new TH1F("fTwoGammaToMuMediumHigherSide", "fTwoGammaToMuMediumHigherSide", PtBinNumber, PtBins );

  Double_t BinCenter = 0;
  for ( Int_t ibin = 1; ibin < fTwoGammaToMuMediumOrig->GetNbinsX(); ibin++ ) {
    BinCenter = ((TAxis*)fTwoGammaToMuMediumOrig->GetXaxis())->GetBinCenter(ibin);
    if ( BinCenter > PtBins[PtBinNumber-1] ) continue;
    cout << "BinCenter = " << BinCenter << endl;
    for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
      if ( BinCenter < PtBins[ibinVariable+1] ){
        fTwoGammaToMuMedium->SetBinContent(ibinVariable+1, fTwoGammaToMuMedium->GetBinContent(ibinVariable+1) + (fTwoGammaToMuMediumOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
        break;
      }
    }
    for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
      if ( BinCenter < PtBins[ibinVariable+1] ){
        fTwoGammaToMuMediumLowerSide->SetBinContent(ibinVariable+1, fTwoGammaToMuMediumLowerSide->GetBinContent(ibinVariable+1) + (fTwoGammaToMuMediumOrigLowerSide->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
        break;
      }
    }
    for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
      if ( BinCenter < PtBins[ibinVariable+1] ){
        fTwoGammaToMuMediumHigherSide->SetBinContent(ibinVariable+1, fTwoGammaToMuMediumHigherSide->GetBinContent(ibinVariable+1) + (fTwoGammaToMuMediumOrigHigherSide->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
        break;
      }
    }
  }

  /* - Firstly we normalize the histograms.
     - Remember to always Sumw2()!!
     -
   */
  fTwoGammaToMuMedium -> Sumw2();
  Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );

  fTwoGammaToMuMediumLowerSide -> Sumw2();
  Double_t Integral_fTwoGammaToMuMediumLowerSide = fTwoGammaToMuMediumLowerSide -> Integral();
  fTwoGammaToMuMediumLowerSide -> Scale( 1/Integral_fTwoGammaToMuMediumLowerSide );

  fTwoGammaToMuMediumHigherSide -> Sumw2();
  Double_t Integral_fTwoGammaToMuMediumHigherSide = fTwoGammaToMuMediumHigherSide -> Integral();
  fTwoGammaToMuMediumHigherSide -> Scale( 1/Integral_fTwoGammaToMuMediumHigherSide );




  Data0N0Nlower  = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide");
  Data0N0Nhigher = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide");
  Data0N0Nlower  -> Sumw2();
  Data0N0Nhigher -> Sumw2();
  Double_t Integral_0N0Nlower  = Data0N0Nlower  -> Integral();
  Double_t Integral_0N0Nhigher = Data0N0Nhigher -> Integral();
  Data0N0Nlower  -> Scale( 1/Integral_0N0Nlower  );
  Data0N0Nhigher -> Scale( 1/Integral_0N0Nhigher );

  // Guillermo's way.
  Data0NXNlower  = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide");
  Data0NXNhigher = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide");
  Data0NXNlower  -> Sumw2();
  Data0NXNhigher -> Sumw2();
  Double_t Integral_0NXNlower  = Data0NXNlower  -> Integral();
  Double_t Integral_0NXNhigher = Data0NXNhigher -> Integral();
  Data0NXNlower  -> Scale( 1/Integral_0NXNlower  );
  Data0NXNhigher -> Scale( 1/Integral_0NXNhigher );

  DataXN0Nlower  = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide");
  DataXN0Nhigher = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide");
  DataXN0Nlower  -> Sumw2();
  DataXN0Nhigher -> Sumw2();
  Double_t Integral_XN0Nlower  = DataXN0Nlower  -> Integral();
  Double_t Integral_XN0Nhigher = DataXN0Nhigher -> Integral();
  DataXN0Nlower  -> Scale( 1/Integral_XN0Nlower  );
  DataXN0Nhigher -> Scale( 1/Integral_XN0Nhigher );

  DataXNXNlower  = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3LowerSide");
  DataXNXNhigher = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv3HigherSide");
  DataXNXNlower  -> Sumw2();
  DataXNXNhigher -> Sumw2();
  Double_t Integral_XNXNlower  = DataXNXNlower  -> Integral();
  Double_t Integral_XNXNhigher = DataXNXNhigher -> Integral();
  DataXNXNlower  -> Scale( 1/Integral_XNXNlower  );
  DataXNXNhigher -> Scale( 1/Integral_XNXNhigher );




  fTwoGammaToMuMediumC           = (TH1F*) fTwoGammaToMuMedium          -> Clone("fTwoGammaToMuMediumC");
  fTwoGammaToMuMediumLowerSideC  = (TH1F*) fTwoGammaToMuMediumLowerSide -> Clone("fTwoGammaToMuMediumLowerSideC");
  fTwoGammaToMuMediumHigherSideC = (TH1F*) fTwoGammaToMuMediumHigherSide-> Clone("fTwoGammaToMuMediumHigherSideC");
  Data0N0NlowerC                 = (TH1F*) Data0N0Nlower                -> Clone("Data0N0NlowerC");
  Data0N0NhigherC                = (TH1F*) Data0N0Nhigher               -> Clone("Data0N0NhigherC");
  Data0NXNlowerC                 = (TH1F*) Data0NXNlower                -> Clone("Data0NXNlowerC");
  Data0NXNhigherC                = (TH1F*) Data0NXNhigher               -> Clone("Data0NXNhigherC");
  DataXN0NlowerC                 = (TH1F*) DataXN0Nlower                -> Clone("DataXN0NlowerC");
  DataXN0NhigherC                = (TH1F*) DataXN0Nhigher               -> Clone("DataXN0NhigherC");
  DataXNXNlowerC                 = (TH1F*) DataXNXNlower                -> Clone("DataXNXNlowerC");
  DataXNXNhigherC                = (TH1F*) DataXNXNhigher               -> Clone("DataXNXNhigherC");

  // TCanvas* c1 = new TCanvas("0N0N low ratio high", "0N0N low ratio high", 900, 800);
  // Data0N0NlowerC->Divide(Data0N0NhigherC);
  // Data0N0NlowerC->Draw();
  // TCanvas* c2 = new TCanvas("0NXN low ratio high", "0NXN low ratio high", 900, 800);
  // Data0NXNlowerC->Divide(Data0NXNhigherC);
  // Data0NXNlowerC->Draw();
  // TCanvas* c3 = new TCanvas("XN0N low ratio high", "XN0N low ratio high", 900, 800);
  // DataXN0NlowerC->Divide(DataXN0NhigherC);
  // DataXN0NlowerC->Draw();
  // TCanvas* c4 = new TCanvas("XNXN low ratio high", "XNXN low ratio high", 900, 800);
  // DataXNXNlowerC->Divide(DataXNXNhigherC);
  // DataXNXNlowerC->Draw();
  // TCanvas* c5 = new TCanvas("MC low ratio high", "MC low ratio high", 900, 800);
  // fTwoGammaToMuMediumLowerSideC->Divide(fTwoGammaToMuMediumHigherSideC);
  // fTwoGammaToMuMediumLowerSideC->Draw();
  // TCanvas* c6 = new TCanvas("XN0N low ratio MC lower", "XN0N low ratio MC lower", 900, 800);
  // DataXN0Nlower->Divide(fTwoGammaToMuMedium);
  // DataXN0Nlower->Draw();

  TCanvas* c1 = new TCanvas("0N0N low ratio high", "0N0N low ratio high", 900, 800);
  gPad->SetLogy();
  Data0N0NlowerC ->SetLineColor(kBlue);
  Data0N0NlowerC ->SetLineStyle(kSolid);
  Data0N0NlowerC ->SetLineWidth(3);
  Data0N0NlowerC ->Draw();
  Data0N0NhigherC->SetLineColor(kRed);
  Data0N0NhigherC->SetLineStyle(kSolid);
  Data0N0NhigherC->SetLineWidth(3);
  Data0N0NhigherC->Draw("same");
  TCanvas* c2 = new TCanvas("0NXN low ratio high", "0NXN low ratio high", 900, 800);
  gPad->SetLogy();
  Data0NXNlowerC ->SetLineColor(kBlue);
  Data0NXNlowerC ->SetLineStyle(kSolid);
  Data0NXNlowerC ->SetLineWidth(3);
  Data0NXNlowerC ->Draw();
  Data0NXNhigherC->SetLineColor(kRed);
  Data0NXNhigherC->SetLineStyle(kSolid);
  Data0NXNhigherC->SetLineWidth(3);
  Data0NXNhigherC->Draw("same");
  TCanvas* c3 = new TCanvas("XN0N low ratio high", "XN0N low ratio high", 900, 800);
  gPad->SetLogy();
  DataXN0NlowerC ->SetLineColor(kBlue);
  DataXN0NlowerC ->SetLineStyle(kSolid);
  DataXN0NlowerC ->SetLineWidth(3);
  DataXN0NlowerC ->Draw();
  DataXN0NhigherC->SetLineColor(kRed);
  DataXN0NhigherC->SetLineStyle(kSolid);
  DataXN0NhigherC->SetLineWidth(3);
  DataXN0NhigherC->Draw("same");
  TCanvas* c4 = new TCanvas("XNXN low ratio high", "XNXN low ratio high", 900, 800);
  gPad->SetLogy();
  DataXNXNlowerC ->SetLineColor(kBlue);
  DataXNXNlowerC ->SetLineStyle(kSolid);
  DataXNXNlowerC ->SetLineWidth(3);
  DataXNXNlowerC ->Draw();
  DataXNXNhigherC->SetLineColor(kRed);
  DataXNXNhigherC->SetLineStyle(kSolid);
  DataXNXNhigherC->SetLineWidth(3);
  DataXNXNhigherC->Draw("same");
  TCanvas* c5 = new TCanvas("MC low ratio high", "MC low ratio high", 900, 800);
  gPad->SetLogy();
  DataXN0NlowerC ->SetLineColor(kBlue);
  DataXN0NlowerC ->SetLineStyle(kSolid);
  DataXN0NlowerC ->SetLineWidth(3);
  DataXN0NlowerC ->Draw();
  DataXN0NhigherC->SetLineColor(kRed);
  DataXN0NhigherC->SetLineStyle(kSolid);
  DataXN0NhigherC->SetLineWidth(3);
  DataXN0NhigherC->Draw("same");
  TCanvas* c6 = new TCanvas("MC higher ratio MC lower", "XN0N low ratio MC lower", 900, 800);
  gPad->SetLogy();
  fTwoGammaToMuMediumLowerSideC ->SetLineColor(kBlue);
  fTwoGammaToMuMediumLowerSideC ->SetLineStyle(kSolid);
  fTwoGammaToMuMediumLowerSideC ->SetLineWidth(3);
  fTwoGammaToMuMediumLowerSideC ->Draw();
  fTwoGammaToMuMediumHigherSideC->SetLineColor(kRed);
  fTwoGammaToMuMediumHigherSideC->SetLineStyle(kSolid);
  fTwoGammaToMuMediumHigherSideC->SetLineWidth(3);
  fTwoGammaToMuMediumHigherSideC->Draw("same");




  // fDimuonPtDistributionDataH->SetLineColor(kBlue);
  // fDimuonPtDistributionDataH->SetLineStyle(kSolid);
  // fDimuonPtDistributionDataH->SetLineWidth(3);
  // fDimuonPtDistributionDataH->SetMarkerStyle(kFullCircle);
  // fDimuonPtDistributionDataH->SetMarkerColor(kBlue);
  // fDimuonPtDistributionDataH->SetMarkerSize(1);
  // fDimuonPtDistributionDataH->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
  //                                                         fDimuonPtDistributionDataH->GetXaxis()->GetBinWidth(1)
  //                                                       )
  //                                                   );
  // fDimuonPtDistributionDataH->SetTitle("");
  // TCanvas* PtDistrCanvas = new TCanvas( "PtDistrCanvas", "PtDistrCanvas", 900, 800 );
  // gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  // gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  // /* - Beautifying is starting now.
  //    -
  //  */
  // fDimuonPtDistributionDataH->GetXaxis()->SetTitleOffset(1.25);
  // // fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.25);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.45);
  // fDimuonPtDistributionDataH->GetXaxis()->SetTitleSize(0.045);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleSize(0.045);
  // fDimuonPtDistributionDataH->GetXaxis()->SetLabelSize(0.045);
  // fDimuonPtDistributionDataH->GetYaxis()->SetLabelSize(0.045);
  // fDimuonPtDistributionDataH->GetXaxis()->SetTitleFont(42);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleFont(42);
  // fDimuonPtDistributionDataH->GetXaxis()->SetLabelFont(42);
  // fDimuonPtDistributionDataH->GetYaxis()->SetLabelFont(42);
  // fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(fDimuonPtDistributionDataH->GetMaximum()*0.0001, fDimuonPtDistributionDataH->GetMaximum()*1000.);
  // // fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 5.5);
  // fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 3);
  // fDimuonPtDistributionDataH->Draw("PEsame");
  // fTwoGammaToMuMedium -> SetLineColor(kGreen);
  // fTwoGammaToMuMedium -> SetLineWidth(3);
  // TH1F* fTwoGammaToMuMediumC = (TH1F*) fTwoGammaToMuMedium -> Clone("fTwoGammaToMuMediumC");
  // fTwoGammaToMuMediumC -> Draw("HISTsame");
  //
  // TLatex* latex = new TLatex();
  // latex->SetTextSize(0.05);
  // latex->SetTextFont(42);
  // latex->SetTextAlign(11);
  // latex->SetNDC();
  // latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  // latex->DrawLatex(0.17,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f #pm %.7f", int1i, intun, int1c, f_I, errfI));
  // latex->SetTextSize(0.045);
  //
  // TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  // l->SetMargin(0.1);
  // l->SetBorderSize(0);
  // // l->AddEntry(  fDimuonPtDistributionDataH, "ALICE data 2018");
  // if      ( selectionFlag == 0 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2");
  // else if ( selectionFlag == 1 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, 0N0N");
  // else if ( selectionFlag == 2 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, 0NXN");
  // else if ( selectionFlag == 3 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, XN0N");
  // else if ( selectionFlag == 4 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, XNXN");
  // else                           l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2");
  // l->AddEntry(  FitPtDistr, Form( "Fit: #chi^{2}/NDF = %.2f / %.2d = %.2f  ",
  //                                  FitPtDistr->GetChisquare(),
  //                                  FitPtDistr->GetNDF(),
  //                                  FitPtDistr->GetChisquare()/FitPtDistr->GetNDF()
  //                                  )
  //                                 );
  // l->AddEntry(  fTwoGammaToMuMediumC, "Continuum  #gamma#gamma #rightarrow #mu#mu");
  // l->Draw();
  //
  // gPad->SaveAs("pngResults/fitPtDistrALL_v3.png",    "RECREATE");


}
