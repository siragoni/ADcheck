#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;
TH1F* hPt;
TH1F* h1c;
TH1F* h1i;
TH1F* hfc;
TH1F* hfi;
TH1F* hgl;
TH1F* hun;
TF1* fgu;

TH1F* fCohJpsiToMuOrig;
TH1F* fCohPsi2sToMuOrig;
TH1F* fCohPsi2sToMuPiOrig;
TH1F* fIncohJpsiToMuOrig;
TH1F* fIncohPsi2sToMuOrig;
TH1F* fIncohPsi2sToMuPiOrig;
TH1F* fTwoGammaToMuMediumOrig;
TH1F* fTwoGammaToMuHighOrig;



void fitPtROOFit( const char* AnalysisName, const int selectionFlag, const int selectionFlag2, Double_t rSigma2s1s = 0.150 )
{
  // ---------------------------------------------------
  // I m p o r t i n g   R O O T   h i s t o g r a m s
  // ===================================================
  // I m p o r t   T H 1   i n t o   a   R o o D a t a H i s t
  // ---------------------------------------------------------
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TFile* fileMC[8];
  // fileMC[0] = new TFile("MCtrainResults/2020-06-19/kCohJpsiToMu/AnalysisResults.root");
  // fileMC[1] = new TFile("MCtrainResults/2020-06-19/kCohPsi2sToMu/AnalysisResults.root");
  // fileMC[2] = new TFile("MCtrainResults/2020-06-19/kCohPsi2sToMuPi/AnalysisResults.root");
  // fileMC[3] = new TFile("MCtrainResults/2020-06-19/kIncohJpsiToMu/AnalysisResults.root");
  // fileMC[4] = new TFile("MCtrainResults/2020-06-19/kIncohPsi2sToMu/AnalysisResults.root");
  // fileMC[5] = new TFile("MCtrainResults/2020-06-19/kIncohPsi2sToMuPi/AnalysisResults.root");
  // fileMC[6] = new TFile("MCtrainResults/2020-06-19/kTwoGammaToMuMedium/AnalysisResults.root");
  // fileMC[7] = new TFile("MCtrainResults/2020-06-19/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[0] = new TFile("MCtrainResults/2020-06-26/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2020-06-26/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2020-06-26/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2020-06-26/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2020-06-26/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2020-06-26/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2020-06-26/kTwoGammaToMuMedium/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2020-06-26/kTwoGammaToMuHigh/AnalysisResults.root");
  TDirectory* dirMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir->GetObject("ADcheck_", listings);
  TList* listingsMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory]->GetObject("MyOutputContainer", listingsMC[iDirectory]);
  }
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
  // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fDimuonPtDistributionH");
  // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fDimuonPtDistributionH");
  // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fDimuonPtDistributionH");
  // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fDimuonPtDistributionH");
  // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fDimuonPtDistributionH");
  // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fDimuonPtDistributionH");
  // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fDimuonPtDistributionH");
  // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fDimuonPtDistributionH");
  if        ( selectionFlag2 == 1 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    fCohJpsiToMuOrig       ->Add((TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fCohPsi2sToMuOrig      ->Add((TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fCohPsi2sToMuPiOrig    ->Add((TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fIncohJpsiToMuOrig     ->Add((TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fIncohPsi2sToMuOrig    ->Add((TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fIncohPsi2sToMuPiOrig  ->Add((TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    // fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    // fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_1"));
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_0");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 2 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_2");
    fCohJpsiToMuOrig       ->Add((TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fCohPsi2sToMuOrig      ->Add((TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fCohPsi2sToMuPiOrig    ->Add((TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fIncohJpsiToMuOrig     ->Add((TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fIncohPsi2sToMuOrig    ->Add((TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fIncohPsi2sToMuPiOrig  ->Add((TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    // fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    // fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_3"));
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_1");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_1");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_1");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 3 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fCohJpsiToMuOrig       ->Add((TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fCohPsi2sToMuOrig      ->Add((TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fCohPsi2sToMuPiOrig    ->Add((TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fIncohJpsiToMuOrig     ->Add((TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fIncohPsi2sToMuOrig    ->Add((TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fIncohPsi2sToMuPiOrig  ->Add((TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    // fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    // fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fTwoGammaToMuMediumOrig->Add((TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));
    fTwoGammaToMuHighOrig  ->Add((TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_5"));

    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 0 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_3");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 4 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_4");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 5 ) {
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityMoreH_5");
    // fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    // fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    // fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    // fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    // fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_2");
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  }
  /* - Rebin
     -
   */
  // fCohJpsiToMu        -> Rebin(4);
  // fCohPsi2sToMu       -> Rebin(4);
  // fCohPsi2sToMuPi     -> Rebin(4);
  // fIncohJpsiToMu      -> Rebin(4);
  // fIncohPsi2sToMu     -> Rebin(4);
  // fIncohPsi2sToMuPi   -> Rebin(4);
  // fTwoGammaToMuMedium -> Rebin(4);
  // fTwoGammaToMuHigh   -> Rebin(4);
  // fCohJpsiToMuOrig        -> Rebin(5);  // -> 25 MeV bins
  // fCohPsi2sToMuOrig       -> Rebin(5);
  // fCohPsi2sToMuPiOrig     -> Rebin(5);
  // fIncohJpsiToMuOrig      -> Rebin(5);
  // fIncohPsi2sToMuOrig     -> Rebin(5);
  // fIncohPsi2sToMuPiOrig   -> Rebin(5);
  // fTwoGammaToMuMediumOrig -> Rebin(5);
  // fTwoGammaToMuHighOrig   -> Rebin(5);
  fCohJpsiToMuOrig        -> Sumw2();
  fCohPsi2sToMuOrig       -> Sumw2();
  fCohPsi2sToMuPiOrig     -> Sumw2();
  fIncohJpsiToMuOrig      -> Sumw2();
  fIncohPsi2sToMuOrig     -> Sumw2();
  fIncohPsi2sToMuPiOrig   -> Sumw2();
  fTwoGammaToMuMediumOrig -> Sumw2();
  fTwoGammaToMuHighOrig   -> Sumw2();



  h1c = (TH1F*) fCohJpsiToMuOrig        -> Clone("fCohJpsiToMu");
  // hfc = (TH1F*) fCohPsi2sToMuOrig       -> Clone("fCohPsi2sToMu");
  hfc     = (TH1F*) fCohPsi2sToMuPiOrig     -> Clone("fCohPsi2sToMuPi");
  h1i = (TH1F*) fIncohJpsiToMuOrig      -> Clone("fIncohJpsiToMu");
  // hfi = (TH1F*) fIncohPsi2sToMuOrig     -> Clone("fIncohPsi2sToMu");
  hfi   = (TH1F*) fIncohPsi2sToMuPiOrig   -> Clone("fIncohPsi2sToMuPi");
  // h1i   = (TH1F*) fIncohPsi2sToMuPiOrig   -> Clone("fIncohPsi2sToMuPi");
  hgl = (TH1F*) fTwoGammaToMuMediumOrig -> Clone("fTwoGammaToMuMedium");
  // fTwoGammaToMuHigh   = (TH1F*) fTwoGammaToMuHighOrig   -> Clone("fTwoGammaToMuHigh");

  if      ( selectionFlag == 0 ) {
      if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
      else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_1");
      else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_2");
      else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_3");
      else if ( selectionFlag2 == 4 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_4");
      else if ( selectionFlag2 == 5 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_5");
      else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
  // } else if   ( selectionFlag  == 1 ) {
  //     if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
  //     else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_1");
  //     else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_2");
  //     // else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_3");
  //     // else if ( selectionFlag2 == 4 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_4");
  //     // else if ( selectionFlag2 == 5 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_5");
  //     else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
  }
  else if ( selectionFlag == 1 ) {
       if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv2");
       else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_2");
       else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv2");
  }
  else if ( selectionFlag == 2 ) {
       if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv2");
       else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_2");
       else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv2");
  }
  else if ( selectionFlag == 3 ) {
       if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv2");
       else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_2");
       else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv2");
  }
  else if ( selectionFlag == 4 ) {
       if      ( selectionFlag2 == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv2");
       else if ( selectionFlag2 == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_2");
       else                            hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv2");
  }

  else                                hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  // hPt->Rebin(5);
  // hPt->Draw("PE");


  hPt->Rebin(5);
  h1c->Rebin(5);
  h1i->Rebin(5);
  hfc->Rebin(5);
  hfi->Rebin(5);
  hgl->Rebin(5);

  // hPt->Rebin(2);
  // h1c->Rebin(2);
  // h1i->Rebin(2);
  // hfc->Rebin(2);
  // hfi->Rebin(2);
  // hgl->Rebin(2);
  //
  // hPt->Rebin(2);
  // h1c->Rebin(2);
  // h1i->Rebin(2);
  // // hfc->Rebin(2);
  // hfi->Rebin(2);
  // hgl->Rebin(2);








  Int_t nBinsX = hPt->GetNbinsX();

  TF1* fun = new TF1("fun","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0,20);
  fun->SetParameter(0,1);
  //  fun->SetParameter(1,debug==4 ? 1.25 : 1.);
  //  fun->SetParameter(2,debug==4 ? 6.1 : 1.);
  // fun->SetParameter(1,debug==4 ? 1.6 : 1.79);
  // fun->SetParameter(2,debug==4 ? 3.58 : 3.58);
  fun->SetParameter(1, 1.79);
  fun->SetParameter(2, 3.58);
  fun->SetNpx(nBinsX);
  // hun = (TH1D*) fun->GetHistogram()->Clone("hun");
  hun = (TH1F*) fun->GetHistogram()->Clone("hun");
  for (Int_t ibin=1;ibin<=hun->GetNbinsX();ibin++) hun->SetBinError(ibin,0);




  // Declare observable x
  RooRealVar   pT("pT",  "pT",   0, 20);
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist rPt("rPt", "rPt", pT, Import(*hPt));
  RooDataHist r1c("r1c", "r1c", pT, Import(*h1c));
  RooDataHist r1i("r1i", "r1i", pT, Import(*h1i));
  RooDataHist rfc("rfc", "rfc", pT, Import(*hfc));
  RooDataHist rfi("rfi", "rfi", pT, Import(*hfi));
  RooDataHist rgl("rgl", "rgl", pT, Import(*hgl));
  RooDataHist run("run", "run", pT, Import(*hun));
  // P l o t   a n d   f i t   a   R o o D a t a H i s t
  // ---------------------------------------------------
  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  // RooPlot *frame = pT.frame(Title("Dimuon pT distribution"));
  RooPlot *frame = pT.frame(0.0, 3.0);
  rPt.plotOn(frame);
  // // Fit a Gaussian p.d.f to the data
  // RooRealVar mean("mean", "mean", 0, -10, 10);
  // RooRealVar sigma("sigma", "sigma", 3, 0.1, 10);
  // RooGaussian gauss("gauss", "gauss", x, mean, sigma);
  // gauss.fitTo(dh);
  // gauss.plotOn(frame);






  RooHistPdf pdf1c("pdf1c", "pdf1c", pT, r1c, 0);
  RooHistPdf pdf1i("pdf1i", "pdf1i", pT, r1i, 0);
  RooHistPdf pdffc("pdffc", "pdffc", pT, rfc, 0);
  RooHistPdf pdffi("pdffi", "pdffi", pT, rfi, 0);
  RooHistPdf pdfgl("pdfgl", "pdfgl", pT, rgl, 0);
  RooHistPdf pdfun("pdfun", "pdfun", pT, run, 0);








  Double_t ngl = 0.0;
  if        ( selectionFlag  == 0 ) { // COHERENT analysis
    if      ( selectionFlag2 == 0 ) ngl = 138 ;
    else if ( selectionFlag2 == 1 ) ngl = 637 ;
    else if ( selectionFlag2 == 2 ) ngl = 1296 ;
    else if ( selectionFlag2 == 3 ) ngl = 1551 ;
    else if ( selectionFlag2 == 4 ) ngl = 1092 ;
    else if ( selectionFlag2 == 5 ) ngl = 303 ;
  } else if ( selectionFlag  == 1 ) { // 0N0N
    if      ( selectionFlag2 == 0 ) ngl = 911 ;
    else if ( selectionFlag2 == 1 ) ngl = 911 ;
    else if ( selectionFlag2 == 2 ) ngl = 3378 ;
    else if ( selectionFlag2 == 3 ) ngl = 1638 ;
    else                            ngl = 911 ;
  } else if ( selectionFlag  == 2 ) { // 0NXN
    if      ( selectionFlag2 == 0 ) ngl =  66 ;
    else if ( selectionFlag2 == 1 ) ngl =  66 ;
    else if ( selectionFlag2 == 2 ) ngl =  234 ;
    else if ( selectionFlag2 == 3 ) ngl =  151 ;
    else                            ngl =  66 ;
  } else if ( selectionFlag  == 3 ) { // XN0N
    if      ( selectionFlag2 == 0 ) ngl =   82 ;
    else if ( selectionFlag2 == 1 ) ngl =   82 ;
    else if ( selectionFlag2 == 2 ) ngl =  307 ;
    else if ( selectionFlag2 == 3 ) ngl =  146 ;
    else                            ngl =   82 ;
  } else if ( selectionFlag  == 4 ) { // XNXN
    if      ( selectionFlag2 == 0 ) ngl =   19 ;
    else if ( selectionFlag2 == 1 ) ngl =   19 ;
    else if ( selectionFlag2 == 2 ) ngl =  112 ;
    else if ( selectionFlag2 == 3 ) ngl =   51 ;
    else                            ngl =   19 ;
  } else {
    ngl =  592 ;
  }




  // RooRealVar ggN     ( "ggN",    "number of gg",        ngl); // fixed by fit to mass peak
  RooRealVar ggN     ( "ggN",    "number of gg",        ngl, ngl*0.8, ngl*1.2 ); // fixed by fit to mass peak
  // RooRealVar ggN     ( "ggN",    "number of gg",        1000, 0, 10000); // fixed by fit to mass peak
  // RooRealVar cohN    ( "cohN",   "number of coh",       10,  0, 20 ); // free
  RooRealVar cohN    ( "cohN",   "number of coh",       1000, 0, 100000 ); // free
  RooRealVar icN     ( "icN",    "number of ic",        1000, 0, 10000); //free
  // RooFormulaVar p2scohN("p2scohN","cohN*0.05", RooArgList(cohN));
  RooFormulaVar p2scohN("p2scohN","cohN*0.045", RooArgList(cohN));
  // RooRealVar p2scohN ( "p2scohN","number of coh p2s",   nPsiCoh); // fixed to STARLIGHT ratio
  // RooRealVar p2sicN  ( "p2sicN", "number of incoh p2s", 24,  20, 28);
  RooFormulaVar p2sicN("p2sicN","icN*0.045", RooArgList(icN));
  RooRealVar dissoc  ( "dissoc", "number of dissoc",    1000, 0, 10000 );
  RooAddPdf  sum     ( "sum",    "extended sum of all",
                       RooArgList(pdf1i, pdf1c, pdffc,   pdffi,    pdfgl, pdfun),
                       RooArgList(icN,   cohN,  p2scohN, p2sicN,   ggN,   dissoc)
                       );




  // RooFitResult* r = sum.fitTo(rPt,Extended(kTRUE),Save());
  RooFitResult* r = sum.fitTo(rPt,Range(0., 3.),Extended(kTRUE),Save());




  sum.plotOn (frame,LineColor(kBlack), Range(0.,3.0) ) ;
  // sum.paramOn(frame);
  // sum.paramOn(frame,rPt);
  sum.plotOn (frame,Components(pdf1c), LineColor(kGreen),   Range(0.,3.0)  );
  sum.plotOn (frame,Components(pdf1i), LineColor(kRed),     Range(0.,3.0)  );
  sum.plotOn (frame,Components(pdffc), LineColor(kYellow),  Range(0.,3.0)  );
  sum.plotOn (frame,Components(pdffi), LineColor(kBlue),    Range(0.,3.0)  );
  sum.plotOn (frame,Components(pdfgl), LineColor(kCyan),    Range(0.,3.0)  );
  sum.plotOn (frame,Components(pdfun), LineColor(kMagenta), Range(0.,3.0)  );
  // sum.paramOn(frame,rPt);
  frame->Draw();



  // new TCanvas("rf201_composite", "rf201_composite", 600, 600);
  // gPad->SetLeftMargin(0.15);
  // xframe->GetYaxis()->SetTitleOffset(1.4);
  // xframe->Draw();




  // get fractions of PDF in jRecPt range [0,Pt cut]
  pT.setRange("signal",0.0,0.25);
  RooAbsReal* cohI    = pdf1c.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* icI     = pdf1i.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* p2scohI = pdffc.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* p2sicI  = pdffi.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* ggI     = pdfgl.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* funI    = pdfun.createIntegral(pT,NormSet(pT),Range("signal")) ;
  RooAbsReal* sI      = sum.createIntegral  (pT,NormSet(pT),Range("signal")) ;





  //_______________________________
  // TOMAS SNIPPET
  RooRealVar ncoh(    "ncoh",    "ncoh",    cohI->getVal(),    "" );
  RooRealVar ncoh2(   "ncoh2",   "ncoh2",   p2scohI->getVal(), "" );
  RooRealVar nincoh(  "nincoh",  "nincoh",  icI->getVal(),     "" );
  RooRealVar nincoh2( "nincoh2", "nincoh2", p2sicI->getVal(),  "" );
  RooRealVar ndisoc(  "ndisoc",  "ndisoc",  funI->getVal(),    "" );

  RooFormulaVar var_ncoh(     "var_ncoh",    "cohN*ncoh",      RooArgList(cohN,ncoh)      );
  RooFormulaVar var_ncoh2(    "var_ncoh2",   "p2scohN*ncoh2",  RooArgList(p2scohN,ncoh2)  );
  RooFormulaVar var_nincoh(   "var_nincoh",  "icN*nincoh",     RooArgList(icN,nincoh)     );
  RooFormulaVar var_nincoh2(  "var_nincoh2", "p2sicN*nincoh2", RooArgList(p2sicN,nincoh2) );
  RooFormulaVar var_ndisoc(   "var_ndisoc",  "dissoc*ndisoc",  RooArgList(dissoc,ndisoc)  );
  RooFormulaVar f_I_extended( "f_I_extended",
                              "(var_nincoh+var_nincoh2+var_ndisoc)/(var_ncoh+var_ncoh2)",
                              RooArgList(var_nincoh, var_nincoh2, var_ndisoc, var_ncoh, var_ncoh2)
                              );

  //_______________________________





  //Calculate chi^2
  Double_t total = sum.expectedEvents(pT);
  cout << "total = " << total << endl;
  Double_t chi2 = 0;
  for( Int_t i = 0; i < 120; i++){
  // for( Int_t i = 0; i < nBinsX; i++){
    // pT.setRange( "bin", i*(1.0/((Double_t)nBinsX)), (i+1)*(1.0/((Double_t)nBinsX)) );
    // RooAbsReal* binI   = sum.createIntegral( pT, NormSet(pT), Range("bin") );
    pT.setRange( Form("bin%i", i), i*0.025, (i+1)*0.025 );
    RooAbsReal* binI   = sum.createIntegral( pT, NormSet(pT), Range(Form("bin%i", i)) );
   	Double_t fBin      = binI->getVal();
   	Double_t sumPoint  = fBin * total;
   	cout << "sumPoint" << sumPoint << " ";
   	Double_t dataPoint = hPt->GetBinContent(i+1);
   	cout <<  dataPoint << endl;
   	Double_t sqDiff    = (sumPoint - dataPoint)*(sumPoint-dataPoint);
   	cout <<  sqDiff    << endl;
   	Double_t sqError   = sumPoint;
   	Double_t chiBin    = sqDiff/sqError;
   	// cout <<  chiBin    << endl;
   	chi2+=chiBin;
  }

  cout << "chi^2 =     " << chi2 << endl;
  // cout << "chi^2/dof = " << chi2/((Double_t)nBinsX) << endl;
  cout << "chi^2/dof = " << chi2/((Double_t)120. ) << endl;

  Double_t fCoh     = cohI   ->getVal();
  Double_t fICoh    = icI    ->getVal();
  // Double_t fPsiCoh  = p2scohI->getVal();
  Double_t fPsiICoh = p2sicI ->getVal();
  Double_t fGG      = ggI    ->getVal();
  Double_t fdisso   = funI   ->getVal();


 cout << "For jRectPt<Pt cut GeV/c" << endl;
 cout << "Number of gg integrated under Pt cut: "              << ggN.getValV()    * fGG      << " +/- " << ggN.getError()    * fGG      << endl;
 // cout << "Number of coh feed-down integrated under Pt cut: " << p2scohN.getValV() * fPsiCoh << " +/- " << p2scohN.getError() * fPsiCoh << endl;
 // cout << "Number of incoh feed-down integrated under Pt cut: " << p2sicN.getValV() * fPsiICoh << " +/- " << p2sicN.getError() * fPsiICoh << endl;
 cout << "Number of incoherent integrated under Pt cut: "      << icN.getValV()    * fICoh    << " +/- " << icN.getError()    * fICoh    << endl;
 cout << "Number of coherent integrated under Pt cut: "        << cohN.getValV()   * fCoh     << " +/- " << cohN.getError()   * fCoh     << endl;
 // Double_t N_coh2      = cohN.getVal()    * fCoh; //Number of coherent J/psi according to fit
 Double_t N_coh2      = cohN.getVal()    * fCoh; //Number of coherent J/psi according to fit
 // Double_t N_coh2      = cohN.getValV()    * fCoh; //Number of coherent J/psi according to fit
 Double_t N_cohError  = cohN.getError()   * fCoh; //Error on coh J/psi according to fit
 // Double_t N_I         =      fICoh; //Number of incoherent J/psi below Pt cut
 Double_t N_I         = icN.getVal()     * fICoh; //Number of incoherent J/psi below Pt cut
 // Double_t N_I         = icN.getValV()     * fICoh; //Number of incoherent J/psi below Pt cut
 Double_t N_IError    = icN.getError()    * fICoh; //Error on inc J/psi below Pt cut
 // Double_t N_cohFD = p2scohN.getValV() * fPsiCoh; //Number of coh FD below Pt cut
 // Double_t N_cohFDError = p2scohN.getError() * fPsiCoh; //Error on coh FD < Pt cut
 // Double_t N_icFD      = p2sicN.getValV()  * fPsiICoh; //Number of ic FD below Pt cut
 // Double_t N_icFDError = p2sicN.getError() * fPsiICoh; //Error on ic FD below Pt cut
 // Double_t N_diss      = fdisso; //Number of ic FD below Pt cut
 Double_t N_diss      = dissoc.getVal()  * fdisso; //Number of ic FD below Pt cut
 // Double_t N_diss      = dissoc.getValV()  * fdisso; //Number of ic FD below Pt cut
 Double_t N_dissError = dissoc.getError() * fdisso; //Error on ic FD below Pt cut


 cout << "coh = " << cohN.getVal() << endl;
 cout << "ich = " << icN.getVal() << endl;
 cout << "dissoc = " << dissoc.getVal() << endl;
 //When cross sections and efficiencies are given for pt<Pt cut:
 //cout << "fFDcoh = " << p2scohN.getValV() * fPsiCoh/(cohN.getValV() * fCoh) << endl;
 //cout << "fFDincoh = " << p2sicN.getValV() * fPsiICoh/(cohN.getValV() * fCoh) << endl;
 //When cross sections and effeciences are given for all pt:
 // fracFDcoh = p2scohN.getValV()/cohN.getValV();
 //  fracFDic = p2sicN.getValV()/icN.getValV();
 //
 //  cout << "fFDcoh = " << fracFDcoh << ";  Target: 0.105" << endl;
 //  cout << "fFDincoh = " << fracFDic << ";  Target: 0.0993" << endl;

 cout << "****************************** " << endl;
 cout << "Fractions of the templates under Pt=Pt cut GeV/c" << endl;
 cout << "Coherent fraction " << fCoh << endl;
 cout << "Incoherent fraction " << fICoh << endl;
 // cout << "Coherent Feed-down fraction " << fPsiCoh << endl;
 cout << "Incoherent Feed-down fraction " << fPsiICoh << endl;
 cout << "Gamma+Gamma fraction " << fGG << endl;

 Double_t f_I   =  (N_I     + N_diss     ) / (N_coh2);
 Double_t ErrfI = sqrt((N_IError + N_dissError)*(N_IError + N_dissError) / ((N_I + N_diss)*(N_I + N_diss)) + (N_cohError*N_cohError)/(N_coh2*N_coh2) )*f_I;
 // fFD = ((p2scohN.getValV() * fPsiCoh) + (p2sicN.getValV() * fPsiICoh)) / (cohN.getValV() * fCoh);








 Double_t lumi = 11.92;

 TLatex* latex = new TLatex();
 latex->SetNDC();
 latex->SetTextSize(0.045);
 latex->SetTextAlign(22);
 latex->DrawLatex(0.56,0.95,Form("ALICE, Pb#minusp #sqrt{#it{s}_{NN}} = 8.16 TeV"));
 // latex->DrawLatex(0.45,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f #pm %.3f", N_I, N_diss, N_coh2, f_I, ErrfI));
 latex->DrawLatex(0.45,0.86,Form("f_{I} =  %.3f #pm %.3f", f_I_extended.getVal(), f_I_extended.getPropagatedError(*r)));
 if ( selectionFlag == 7 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -2.500));
 if ( selectionFlag == 8 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.250));
 if ( selectionFlag == 9 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.250,  -2.500));
 if ( selectionFlag == 10) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.500));
 if ( selectionFlag == 11) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.500,  -3.000));
 if ( selectionFlag == 12) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.000,  -2.500));
 if ( selectionFlag == 13) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.625));
 if ( selectionFlag == 14) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.625,  -3.250));
 if ( selectionFlag == 15) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.250,  -2.875));
 if ( selectionFlag == 16) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-2.875,  -2.500));
 if ( selectionFlag == 17) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.700));
 if ( selectionFlag == 18) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.700,  -3.400));
 if ( selectionFlag == 19) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.400,  -3.100));
 if ( selectionFlag == 20) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.100,  -2.800));
 if ( selectionFlag == 21) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-2.800,  -2.500));
 latex->SetTextSize(0.033);
 latex->SetTextAlign(12);
 // latex->DrawLatex(0.48,0.89,Form("UPC, L_{#lower[-0.3]{int}} = %.0f #pm %.0f #mub^{-1}",lumi,lumi*0.05));
 // latex->DrawLatex(0.48,0.85,Form("#minus%.2f < #it{y} < #minus%.2f",-gYMin[iy],-gYMax[iy]));
 // latex->DrawLatex(0.48,0.81,Form("%.2f < #it{m}_{#mu#mu} < %.2f GeV/#it{c}^{2}",gMMin[im],gMMax[im]));

 if (0) {
//    TLegend* l = new TLegend(0.45,0.47,0.98,0.77);
// //    TLegend* l = new TLegend(0.55,0.50,0.98,0.84);
//    l->SetMargin(0.1);
//    l->AddEntry(hPt,"ALICE data","p");
//    l->AddEntry(fsum,Form("Fit: #chi^{2}/NDF=%.1f\n",chi2ndf),"l");
//    l->AddEntry(h1c,Form("Coherent J/#psi: %.0f #pm %.0f",n1c,n1c_err));
//    l->AddEntry(h1i,Form("Incoherent J/#psi: %.0f #pm %.0f",n1i,n1i_err));
//    l->AddEntry(hun,Form("Incoherent dissocitive J/#psi: %.0f",nun),"l");
//    // l->AddEntry(hfc,"Coherent #psi' feeddown");
//    l->AddEntry(hfi,"Incoherent #psi' feeddown");
//    l->AddEntry(hgl,Form("Continuum #gamma#gamma #rightarrow #mu#mu: %.0f",ngl));
//    l->Draw();
 }  else {
   // TLegend* l = new TLegend(0.47,0.48,0.985,0.78);
   // l->SetMargin(0.09);
   // l->AddEntry(hPt,"ALICE data","p");
   // l->AddEntry(pdf1c,"Coherent J/#psi");
   // l->AddEntry(pdf1i,"Incoherent J/#psi");
   // l->AddEntry(pdfun,"Incoherent J/#psi with nucleon dissociation");
   // // l->AddEntry(hfc,"Coherent J/#psi from #psi' decay");
   // l->AddEntry(pdffi,"Incoherent J/#psi from #psi' decay");
   // l->AddEntry(pdfgl,"Continuum #gamma#gamma #rightarrow #mu#mu");
   // l->AddEntry(sum,"Sum");
   // l->Draw();
 }





 for (int i=0; i<frame->numItems(); i++) {
   TString obj_name=frame->nameOf(i); if (obj_name=="") continue;
   cout << Form("%d. '%s'\n",i,obj_name.Data());
 }



 TString names[] = {
   "h_rPt",
   "sum_Norm[pT]_Range[0.000000_1.000000]",
   "sum_Norm[pT]_Comp[pdf1c]_Range[0.000000_1.600000]",
   "sum_Norm[pT]_Comp[pdf1i]_Range[0.000000_1.600000]",
   "sum_Norm[pT]_Comp[pdffi]_Range[0.000000_1.600000]",
   "sum_Norm[pT]_Comp[pdfgl]_Range[0.000000_1.600000]",
   "sum_Norm[pT]_Comp[pdfun]_Range[0.000000_1.600000]",
   ""
 };

 TString signs[] = {
   "ALICE data",
   "Total fit",
   "#gamma+Pb",
   "Exclusive",
   "#gamma#gamma",
   "Non-exclusive"
 };


 TLegend* l = new TLegend(0.47,0.48,0.985,0.78);
 l->SetMargin(0.09);
 Int_t i=-1;
 while ( names[++i] != "" ) {
   TObject *obj = frame->findObject(names[i].Data());
   if (!obj) {
     Warning("fitBi4",Form("Can't find item = %s in the frame2!\n",names[i].Data()));
     continue;
   }
   l->AddEntry(obj,signs[i],"l");
 }
 l->Draw();




 // ------------------------------------------------------------------------------------------------------------------------------------
 // Pt plot
 //---Create pt canvas
 TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
 // cPt->Divide(2);

 // //-----------------------------------------------------------------------------------
 // // Draw Correlation Matrix
 // cPt->cd(2);
 // TPad *cPt2 = new TPad("cPt2", "cPt2",0.001,0.001,0.999,0.999);
 // // ---Pad option
 // cPt2->SetRightMargin(0.15);
 // cPt2->SetLeftMargin(0.15);
 // cPt2->Draw();
 // cPt2->cd();
 // // ---TH2D correlation matrix
 // TH2* hCorr_pt = fit_pt->correlationHist();
 // hCorr_pt->SetMarkerSize(1.6);
 // hCorr_pt->GetXaxis()->SetLabelSize(0.045);
 // hCorr_pt->GetYaxis()->SetLabelSize(0.045);
 // hCorr_pt->Draw("zcol,text");
 //
 //-----------------------------------------------------------------------------------
 // Draw pt Histogram
 // cPt->cd(1);
 TPad *cPt1 = new TPad("cPt1", "cPt1",0.001,0.001,0.999,0.999);
 // ---Pad option
 cPt1->SetLogy();
 cPt1->SetLeftMargin(0.14);
 cPt1->SetRightMargin(0.01);
 cPt1->SetBottomMargin(0.12);
 cPt1->Draw();
 cPt1->cd();
 // ---TH1 pt spectrum
 // RooPlot* frame_pt = pT.frame(Title("Pt fit")) ;
 RooPlot* frame_pt = pT.frame(Title(" ")) ;
 rPt.plotOn(frame_pt,Name("rPt")/*,Binning(bin_pt),*/,MarkerStyle(20),MarkerSize(0.9));
 sum.plotOn(frame_pt,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
 // ---Psi(2S)
 // Pt_fit_func.plotOn(frame_pt,Name("PDF_CohPsi2sToMu"), Components(PDF_CohPsi2sToMu), LineColor(kBlue+1), LineWidth(2));
 // Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohPsi2sToMu"),Components(PDF_IncohPsi2sToMu), LineColor(kRed+1), LineWidth(2));
 // ---J/Psi
 sum.plotOn(frame_pt,Name("PDF_CohJpsiToMu"),       Components(pdf1c), LineColor(kBlue+1), LineWidth(2));
 sum.plotOn(frame_pt,Name("PDF_IncohJpsiToMu"),     Components(pdf1i), LineColor(kRed+1), LineWidth(2));
 sum.plotOn(frame_pt,Name("PDF_CohPsi2sToMuPi"),    Components(pdffc), LineColor(kCyan+1), LineWidth(2));
 sum.plotOn(frame_pt,Name("PDF_IncohPsi2sToMuPi"),  Components(pdffi), LineColor(kOrange), LineWidth(2));
 sum.plotOn(frame_pt,Name("PDF_IncohJpsiToX"),      Components(pdfun), LineColor(kMagenta+1), LineWidth(2));
 sum.plotOn(frame_pt,Name("PDF_TwoGammaToMuMedium"),Components(pdfgl), LineColor(kGreen+2), LineWidth(2));
 // ---Frame pt option
 frame_pt->SetAxisRange(0,3,"X");
 frame_pt->SetAxisRange(0.025,11000000,"Y");
 frame_pt->GetXaxis()->SetTitle("Dimuon #it{p}_{T} (GeV/#it{c})");
 frame_pt->GetYaxis()->SetTitle("Counts per  0.25 GeV/#it{c} ");
 // frame_pt->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (pt_range_max-pt_range_min)/n_bins_pt*1000));
 frame_pt->GetYaxis()->SetTitleOffset(1.3);
 frame_pt->Draw();
 //---Calculate chi2 pt
 // Double_t chi2_pt = frame_pt->chiSquare("sum","rPt",r->floatParsFinal().getSize());
 // Double_t chi2_pt = frame_pt->chiSquare("sum","rPt",sum.getParameters(rPt)->selectByAttrib("Constant",kFALSE)->getSize());
 Double_t chi2_pt = frame->chiSquare("sum","rPt",sum.getParameters(rPt)->selectByAttrib("Constant",kFALSE)->getSize());
 // ---pt plot description
 TLatex * text_pt1 = new TLatex (0.1,0.38*frame_pt->GetMaximum(),"Preliminary");
 text_pt1->SetTextSize(0.035);
 text_pt1->Draw();
 TLatex * text_pt2 = new TLatex (0.85,0.38*frame_pt->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
 text_pt2->SetTextSize(0.035);
 text_pt2->Draw();
 TLatex * text_pt3 = 0x0;
 if      ( selectionFlag == 1 ) {
    if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -2.5");
    else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -3.5");
    else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -3.5 < y < -3.0");
    else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -3.0 < y < -2.5");
    else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -2.5");
 }
 else if ( selectionFlag == 2 ) {
    if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -2.5");
    else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -3.5");
    else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -3.5 < y < -3.0");
    else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -3.0 < y < -2.5");
    else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -2.5");
 }
 else if ( selectionFlag == 3 ) {
   if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -2.5");
   else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -3.5");
   else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -3.5 < y < -3.0");
   else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -3.0 < y < -2.5");
   else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -2.5");
 }
 else if ( selectionFlag == 4 ) {
   if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
   else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -3.5");
   else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -3.5 < y < -3.0");
   else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -3.0 < y < -2.5");
   else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
 }
 else                              text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
 text_pt3->SetTextSize(0.035);
 text_pt3->Draw();

 // TLatex * text_pt3 = new TLatex (0.85,0.11*frame_pt->GetMaximum(),Form("#bf{%.1f < y < %.1f, %.2f < #it{m_{#mu^{+}#mu^{-}}} < %.2f}",y_min,y_max, m_min_jpsi, m_max_jpsi));
 // text_pt3->Draw();
 // TLatex * text_pt4 = new TLatex (0.1,0.11*frame_pt->GetMaximum(),Form("#bf{f_{I} = %.1f%% #pm %.1f%%}",100*f_I.getVal(), 100*f_I.getPropagatedError(*fit_pt) ));
 // text_pt4->Draw();
 // ---Legend pt
 TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
 leg_pt->SetFillStyle(0);
 leg_pt->SetBorderSize(0);
 leg_pt->SetTextSize(0.04);
 leg_pt->AddEntry("rPt","ALICE Data", "P");
 leg_pt->AddEntry("sum",Form("Fit: #chi^{2}/dof= %.3f",chi2/((Double_t)120. )),"L");
 // ------Psi(2S)
 // leg_pt->AddEntry("PDF_CohPsi2sToMu","Coherent #psi(2S)", "L");
 // leg_pt->AddEntry("PDF_IncohPsi2sToMu","Incoherent #psi(2S)", "L");
 // ------J/Psi
 leg_pt->AddEntry("PDF_CohJpsiToMu","Coherent J/#psi", "L");
 leg_pt->AddEntry("PDF_IncohJpsiToMu","Incoherent J/#psi", "L");
 leg_pt->AddEntry("PDF_CohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
 leg_pt->AddEntry("PDF_IncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
 leg_pt->AddEntry("PDF_IncohJpsiToX","Nucleon disoc.", "L");
 leg_pt->AddEntry("PDF_TwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
 // leg_pt->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2_pt),"");
 leg_pt->Draw();





 gPad->SaveAs(Form("pngResults/fitPtROOFit_%d.png", selectionFlag),  "RECREATE");

}
