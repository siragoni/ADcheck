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
TH1F* fCohJpsiToMu;
TH1F* fCohPsi2sToMu;
TH1F* fCohPsi2sToMuPi;
TH1F* fIncohJpsiToMu;
TH1F* fIncohPsi2sToMu;
TH1F* fIncohPsi2sToMuPi;
TH1F* fTwoGammaToMuMedium;
TH1F* fTwoGammaToMuHigh;
TH1F* fHighPtTail;
TH1F* fHighPtTailPreliminary;
//_____________________________________________________________________________
/* - Fit function for the Pt plots.
 * - I am using simple ROOT to make gaussian fits to the plot.
 */
Double_t fPtDistr(Double_t* x,Double_t* par)
{
  /* - Par 0, 1, 2:   coherent.
     - Par 3, 4, 5:   incoherent.
     - Par 6      :   gamma+gamma.
     -
   */
  Double_t val = 0;
  val += par[0]* ( fCohJpsiToMu       ->Interpolate(x[0]) );   //needed
  // val += par[1]* ( fCohPsi2sToMu      ->Interpolate(x[0]) );
  val += par[0]* ( fCohPsi2sToMuPi    ->Interpolate(x[0]) ) * par[4];   //needed
  val += par[1]* ( fIncohJpsiToMu     ->Interpolate(x[0]) );   //needed
  // val += par[4]* ( fIncohPsi2sToMu    ->Interpolate(x[0]) );
  val += par[1]* ( fIncohPsi2sToMuPi  ->Interpolate(x[0]) ) * par[4];   //needed
  val += par[2]* ( fTwoGammaToMuMedium->Interpolate(x[0]) );   //needed
  // val += par[2]* ( fTwoGammaToMuHigh  ->Interpolate(x[0]) );   //needed
  val += par[3]* ( fHighPtTail        ->Interpolate(x[0]) );   //needed

  return val;
}
//_____________________________________________________________________________
/* - Fit function for the Pt plots.
 * - I am using simple ROOT to make gaussian fits to the plot.
 */
Double_t fPtDistrPreliminary(Double_t* x,Double_t* par)
{
  Double_t val = 0;
  val += par[0]* ( fHighPtTailPreliminary->Interpolate(x[0]) );   //needed
  return val;
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitPtDistr(const char* AnalysisName, const int selectionFlag, const int selectionFlag2){
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TFile* fileMC[8];
  fileMC[0] = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-09-17/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-09-17/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-09-17/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-09-17/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-09-17/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-09-17/kTwoGammaToMuMedium/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-09-17/kTwoGammaToMuHigh/AnalysisResults.root");
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
  if        ( selectionFlag2 == 0 ) {
    fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 1 ) {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_0");
    fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 2 ) {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_1");
    fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 3 ) {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_2");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_2");
    fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
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
  fCohJpsiToMu        -> Rebin(5);  // -> 25 MeV bins
  fCohPsi2sToMu       -> Rebin(5);
  fCohPsi2sToMuPi     -> Rebin(5);
  fIncohJpsiToMu      -> Rebin(5);
  fIncohPsi2sToMu     -> Rebin(5);
  fIncohPsi2sToMuPi   -> Rebin(5);
  fTwoGammaToMuMedium -> Rebin(5);
  fTwoGammaToMuHigh   -> Rebin(5);

  /* - Firstly we normalize the histograms.
     - Remember to always Sumw2()!!
     -
   */
  fCohJpsiToMu        -> Sumw2();
  fCohPsi2sToMu       -> Sumw2();
  fCohPsi2sToMuPi     -> Sumw2();
  fIncohJpsiToMu      -> Sumw2();
  fIncohPsi2sToMu     -> Sumw2();
  fIncohPsi2sToMuPi   -> Sumw2();
  fTwoGammaToMuMedium -> Sumw2();
  fTwoGammaToMuHigh   -> Sumw2();
  // Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral(0, 20);
  // Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral(0, 20);
  // Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral(0, 20);
  // Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral(0, 20);
  // Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral(0, 20);
  // Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral(0, 20);
  // Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral(0, 20);
  // Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral(0, 20);
  Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral();
  Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral();
  Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral();
  Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral();
  Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral();
  Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral();
  fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  fCohPsi2sToMu       -> Scale( 1/Integral_fCohPsi2sToMu       );
  fCohPsi2sToMuPi     -> Scale( 1/Integral_fCohPsi2sToMuPi     );
  fIncohJpsiToMu      -> Scale( 1/Integral_fIncohJpsiToMu      );
  fIncohPsi2sToMu     -> Scale( 1/Integral_fIncohPsi2sToMu     );
  fIncohPsi2sToMuPi   -> Scale( 1/Integral_fIncohPsi2sToMuPi   );
  fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );
  fTwoGammaToMuHigh   -> Scale( 1/Integral_fTwoGammaToMuHigh   );

  /* - High Pt-tail, with HERA's data.
     -
   */
  TF1* fModelForHighPtTail = new TF1("fModelForHighPtTail","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0,4);
  fModelForHighPtTail->SetParameter(0,1);
//  fModelForHighPtTail->SetParameter(1,debug==4 ? 1.25 : 1.);
//  fModelForHighPtTail->SetParameter(2,debug==4 ? 6.1 : 1.);
  fModelForHighPtTail->SetParameter(1, 1.6/*1.79*/);
  fModelForHighPtTail->SetParameter(2, 3.58);
  fModelForHighPtTail->SetNpx( fCohJpsiToMu->GetNbinsX() );
  fHighPtTail = (TH1F*) fModelForHighPtTail->GetHistogram()->Clone("fHighPtTail");
  for (Int_t ibin=1; ibin<=fHighPtTail->GetNbinsX(); ibin++) {
    fHighPtTail->SetBinError(ibin,0);
  }



  TH1F *fDimuonPtDistributionDataH = 0x0;
  if      ( selectionFlag == 0 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  }
  else if ( selectionFlag == 1 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv2");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_2");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroHv2");
  }
  else if ( selectionFlag == 2 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv2");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_2");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyHv2");
  }
  else if ( selectionFlag == 3 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv2");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_2");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroHv2");
  }
  else if ( selectionFlag == 4 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv2");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_0");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_1");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_2");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyHv2");
  }
  else                                 fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  fDimuonPtDistributionDataH->Rebin(5);
  fDimuonPtDistributionDataH->Draw("PE");


  /* - Preliminary fit to the unknown contribution.
   * -
   */
  TF1* fModelForHighPtTailPreliminary = new TF1("fModelForHighPtTailPreliminary","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0.5,5.);
  fModelForHighPtTailPreliminary->SetParameter(0,1);
  // fModelForHighPtTailPreliminary->SetParameter(1,debug==4 ? 1.25 : 1.);
  // fModelForHighPtTailPreliminary->SetParameter(2,debug==4 ? 6.1 : 1.);
  fModelForHighPtTailPreliminary->SetParameter(1, 1.6/*1.79*/);
  fModelForHighPtTailPreliminary->SetParameter(2, 3.58);
  fModelForHighPtTailPreliminary->SetNpx( fCohJpsiToMu->GetNbinsX() );
  fHighPtTailPreliminary = (TH1F*) fModelForHighPtTailPreliminary->GetHistogram()->Clone("fHighPtTailPreliminary");
  for (Int_t ibin=1; ibin<=fHighPtTailPreliminary->GetNbinsX(); ibin++) {
    fHighPtTailPreliminary->SetBinError(ibin,0);
  }
  TCanvas* UtilityCanvas = new TCanvas("UtilityCanvas", "UtilityCanvas", 800, 900);
  TH1F* fDimuonPtDistributionDataCloneH = (TH1F*) fDimuonPtDistributionDataH->Clone("fDimuonPtDistributionDataCloneH");
  TF1* FitPtDistrPreliminary = new TF1(  "FitPtDistrPreliminary",
                                         fPtDistrPreliminary,
                                         1.5, 3,
                                         /*7*/1
                                         );
  FitPtDistrPreliminary->SetNpx(1000);
  fDimuonPtDistributionDataCloneH->Fit(FitPtDistrPreliminary, "","", 1.5, 3.);
  // fDimuonPtDistributionDataCloneH->Fit(fModelForHighPtTailPreliminary, "","", 1.5, 3.);
  gStyle->SetOptFit(100);
  // UtilityCanvas->SaveAs("UtilityCanvas.png");


  new TCanvas;
  TF1* FitPtDistr = new TF1(  "FitPtDistr",
                              fPtDistr,
                              0, 3,
                              /*7*/5
                              );
  FitPtDistr->SetNpx(1000);




  Double_t kFeedDownCoherent = 0.05;     // neutral element
  Double_t kError            = 0.01;
  FitPtDistr->SetParameter(4, kFeedDownCoherent);
  FitPtDistr->SetParLimits(4, FitPtDistr->GetParameter(4)*(1-kError), FitPtDistr->GetParameter(4)*(1+kError));


  FitPtDistr->SetParLimits(0, 0.00000000000000001, 99999999999);
  FitPtDistr->SetParLimits(1, 0.00000000000000001, 99999999999);

  // FEED-DOWN COHERENT
  // Double_t kFeedDownIncoherent = 1705;     // neutral element
  // Double_t kErrorInc           = 0.175;
  // FitPtDistr->SetParameter(5, kFeedDownIncoherent);
  // FitPtDistr->SetParLimits(5, FitPtDistr->GetParameter(5)*(1-kErrorInc), FitPtDistr->GetParameter(5)*(1+kErrorInc));

  // Gamma+Gamma Medium
  // FitPtDistr->SetParameter(2, 4800);
  // FitPtDistr->SetParameter(2, 6531);
  if        ( selectionFlag == 0 ) {
    // FitPtDistr->SetParameter(2, 6531);  //2018
    // FitPtDistr->SetParameter(2, 9035);  //2018+2015 with SPD
    FitPtDistr->SetParameter(2, 11000);  //2018+2015 no SPD
  } else if ( selectionFlag == 1 ) { // 0N0N
    // FitPtDistr->SetParameter(2, 5213);  //2018
    // FitPtDistr->SetParameter(2, 7754);  //2018+2015 with SPD
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 6460);          //2018 no SPD, no AD
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 6460 * 0.165);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 6460 * 0.572);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 6460 * 0.263);  //2018 no SPD, no AD
    else                            FitPtDistr->SetParameter(2, 6460);          //2018 no SPD, no AD
  } else if ( selectionFlag == 2 ) { // 0NXN
    // FitPtDistr->SetParameter(2,  370);  //2018
    // FitPtDistr->SetParameter(2,  439);  //2018+2015 with SPD
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 512);          //2018 no SPD, no AD
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 512 * 0.151);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 512 * 0.528);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 512 * 0.321);  //2018 no SPD, no AD
    else                            FitPtDistr->SetParameter(2, 512);          //2018 no SPD, no AD
  } else if ( selectionFlag == 3 ) {
    // FitPtDistr->SetParameter(2,  470);  //2018
    // FitPtDistr->SetParameter(2,  543);  //2018+2015 with SPD
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 647);          //2018 no SPD, no AD
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 647 * 0.165);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 647 * 0.573);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 647 * 0.262);  //2018 no SPD, no AD
    else                            FitPtDistr->SetParameter(2, 647);          //2018 no SPD, no AD
  } else if ( selectionFlag == 4 ) {
    // FitPtDistr->SetParameter(2,  140);  //2018
    // FitPtDistr->SetParameter(2,  145);  //2018+2015 with SPD
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 258);          //2018 no SPD, no AD
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 258 * 0.107);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 258 * 0.613);  //2018 no SPD, no AD
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 258 * 0.280);  //2018 no SPD, no AD
    else                            FitPtDistr->SetParameter(2, 258);          //2018 no SPD, no AD
  }

  // FitPtDistr->FixParameter(2, 0);
  FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*0.95, FitPtDistr->GetParameter(2)*1.05);

  // Unknown component
  // if        ( selectionFlag == 0 ) {
  //   // FitPtDistr->SetParameter(3, 743);
  //   // FitPtDistr->SetParLimits(3, 720, 745);
  //   FitPtDistr->SetParameter(3, 580);            // 2018+2015
  //   FitPtDistr->SetParLimits(3, 550, 600);       // 2018+2015
  // } else if ( selectionFlag == 1 ) {
  //   // FitPtDistr->SetParameter(3, 125);
  //   // FitPtDistr->SetParLimits(3, 100, 150);
  //   FitPtDistr->SetParameter(3, 125);            // 2018+2015
  //   FitPtDistr->SetParLimits(3, 120, 130);       // 2018+2015
  // } else if ( selectionFlag == 2 ) {
  //   FitPtDistr->SetParameter(3, 20);
  //   FitPtDistr->SetParLimits(3, 0, 50);
  // } else if ( selectionFlag == 3 ) {
  //   // FitPtDistr->SetParameter(3, 540);
  //   // FitPtDistr->SetParLimits(3, 500, 600);
  //   FitPtDistr->SetParameter(3, 400);            // 2018+2015
  //   FitPtDistr->SetParLimits(3, 380, 420);       // 2018+2015
  // } else if ( selectionFlag == 4 ) {
  //   // FitPtDistr->SetParameter(3, 40);
  //   // FitPtDistr->SetParLimits(3, 0, 100);
  //   FitPtDistr->SetParameter(3, 25);             // 2018+2015
  //   FitPtDistr->SetParLimits(3, 0, 100);         // 2018+2015
  // }

  FitPtDistr->SetParameter(3, FitPtDistrPreliminary->GetParameter(0));
  // FitPtDistr->SetParLimits(3, FitPtDistrPreliminary->GetParameter(0)*0.9, FitPtDistrPreliminary->GetParameter(0)*1.1);
  // FitPtDistr->SetParameter(3, fModelForHighPtTailPreliminary->GetParameter(0));
  // FitPtDistr->SetParLimits(3, fModelForHighPtTailPreliminary->GetParameter(0)*0.9, fModelForHighPtTailPreliminary->GetParameter(0)*1.1);
  // FitPtDistr->SetParLimits(7, 740, 745);

  fDimuonPtDistributionDataH->SetLineColor(kBlue);
  fDimuonPtDistributionDataH->SetLineStyle(kSolid);
  fDimuonPtDistributionDataH->SetLineWidth(3);
  fDimuonPtDistributionDataH->SetMarkerStyle(kFullCircle);
  fDimuonPtDistributionDataH->SetMarkerColor(kBlue);
  fDimuonPtDistributionDataH->SetMarkerSize(1);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  fDimuonPtDistributionDataH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
                                                          fDimuonPtDistributionDataH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fDimuonPtDistributionDataH->SetTitle("");
  fDimuonPtDistributionDataH->Fit( FitPtDistr,"","", 0./*1.2*/, 3. );
  TCanvas* PtDistrCanvas = new TCanvas( "PtDistrCanvas", "PtDistrCanvas", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  /* - Beautifying is starting now.
     -
   */
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleOffset(1.25);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.25);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.45);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelFont(42);
  // fDimuonPtDistributionDataH->GetXaxis()->SetNdivisions(408);
  // fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(10, fDimuonPtDistributionDataH->GetMaximum()*10.);
  fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(fDimuonPtDistributionDataH->GetMaximum()*0.0001, fDimuonPtDistributionDataH->GetMaximum()*10.);
  // fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 5.5);
  fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 3);
  gPad ->SetLogy();
  fDimuonPtDistributionDataH->Draw("PEsame");
  fCohJpsiToMu        -> SetLineColor(kRed);
  fCohPsi2sToMu       -> SetLineColor(kMagenta);
  fCohPsi2sToMuPi     -> SetLineColor(kYellow+1);
  fIncohJpsiToMu      -> SetLineColor(kCyan);
  fIncohPsi2sToMu     -> SetLineColor(kYellow);
  fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  fTwoGammaToMuMedium -> SetLineColor(kGreen);
  fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  fHighPtTail         -> SetLineColor(kGreen+1);
  fCohJpsiToMu        -> SetLineWidth(3);
  fCohPsi2sToMu       -> SetLineWidth(3);
  fCohPsi2sToMuPi     -> SetLineWidth(3);
  fIncohJpsiToMu      -> SetLineWidth(3);
  fIncohPsi2sToMu     -> SetLineWidth(3);
  fIncohPsi2sToMuPi   -> SetLineWidth(3);
  fTwoGammaToMuMedium -> SetLineWidth(3);
  fTwoGammaToMuHigh   -> SetLineWidth(3);
  fHighPtTail         -> SetLineWidth(3);
  TH1F* fCohJpsiToMuC        = (TH1F*) fCohJpsiToMu        -> Clone("fCohJpsiToMuC");
  TH1F* fCohPsi2sToMuC       = (TH1F*) fCohPsi2sToMu       -> Clone("fCohPsi2sToMuC");
  TH1F* fCohPsi2sToMuPiC     = (TH1F*) fCohPsi2sToMuPi     -> Clone("fCohPsi2sToMuPiC");
  TH1F* fIncohJpsiToMuC      = (TH1F*) fIncohJpsiToMu      -> Clone("fIncohJpsiToMuC");
  TH1F* fIncohPsi2sToMuC     = (TH1F*) fIncohPsi2sToMu     -> Clone("fIncohPsi2sToMuC");
  TH1F* fIncohPsi2sToMuPiC   = (TH1F*) fIncohPsi2sToMuPi   -> Clone("fIncohPsi2sToMuPiC");
  TH1F* fTwoGammaToMuMediumC = (TH1F*) fTwoGammaToMuMedium -> Clone("fTwoGammaToMuMediumC");
  TH1F* fTwoGammaToMuHighC   = (TH1F*) fTwoGammaToMuHigh   -> Clone("fTwoGammaToMuHighC");
  TH1F* fHighPtTailC         = (TH1F*) fHighPtTail         -> Clone("fHighPtTailC");
  fCohJpsiToMuC        -> Scale( FitPtDistr->GetParameter(0) );
  fCohPsi2sToMuC       -> Scale( FitPtDistr->GetParameter(0) * FitPtDistr->GetParameter(4) );
  fCohPsi2sToMuPiC     -> Scale( FitPtDistr->GetParameter(0) * FitPtDistr->GetParameter(4) );
  fIncohJpsiToMuC      -> Scale( FitPtDistr->GetParameter(1) );
  fIncohPsi2sToMuC     -> Scale( FitPtDistr->GetParameter(1) * FitPtDistr->GetParameter(4) );
  fIncohPsi2sToMuPiC   -> Scale( FitPtDistr->GetParameter(1) * FitPtDistr->GetParameter(4) );
  fTwoGammaToMuMediumC -> Scale( FitPtDistr->GetParameter(2) );
  fTwoGammaToMuHighC   -> Scale( FitPtDistr->GetParameter(2) );
  fHighPtTailC         -> Scale( FitPtDistr->GetParameter(3) );
  fCohJpsiToMuC        -> Draw("HISTsame");
  // fCohPsi2sToMuC       -> Draw("same");
  fCohPsi2sToMuPiC     -> Draw("HISTsame");
  fIncohJpsiToMuC      -> Draw("HISTsame");
  // fIncohPsi2sToMuC     -> Draw("same");
  fIncohPsi2sToMuPiC   -> Draw("HISTsame");
  fTwoGammaToMuMediumC -> Draw("HISTsame");
  // fTwoGammaToMuHighC   -> Draw("same");
  fHighPtTailC         -> Draw("Esame");
  // fCohJpsiToMu        -> Scale( FitPtDistr->GetParameter(0) );
  // fCohPsi2sToMu       -> Scale( FitPtDistr->GetParameter(1) );
  // fCohPsi2sToMuPi     -> Scale( FitPtDistr->GetParameter(2) );
  // fIncohJpsiToMu      -> Scale( FitPtDistr->GetParameter(3) );
  // fIncohPsi2sToMu     -> Scale( FitPtDistr->GetParameter(4) );
  // fIncohPsi2sToMuPi   -> Scale( FitPtDistr->GetParameter(5) );
  // fTwoGammaToMuMedium -> Scale( FitPtDistr->GetParameter(6) );
  // fTwoGammaToMuHigh   -> Scale( FitPtDistr->GetParameter(0) );
  // fCohJpsiToMu        -> Draw("same");
  // fCohPsi2sToMu       -> Draw("same");
  // fCohPsi2sToMuPi     -> Draw("same");
  // fIncohJpsiToMu      -> Draw("same");
  // fIncohPsi2sToMu     -> Draw("same");
  // fIncohPsi2sToMuPi   -> Draw("same");
  // fTwoGammaToMuMedium -> Draw("same");
  // fTwoGammaToMuHigh   -> Draw("same");


  /* - COMPUTING $f_I$
   * -
   * - This snippet is for computing the
   * - $f_I$ fraction for correcting the
   * - cross sections later on.
   * - Basically the idea is that
   * -
   * - f_I = ( Incoh + Incoh_Dissociative )/ Coherent
   * -
   */
  Double_t int1c = fCohJpsiToMuC  ->Integral( 1/*fCohJpsiToMuC  ->GetXaxis()->FindBin(0.000001)*/, fCohJpsiToMuC  ->GetXaxis()->FindBin(0.2499999) );
  Double_t int1i = fIncohJpsiToMuC->Integral( 1/*fIncohJpsiToMuC->GetXaxis()->FindBin(0.000001)*/, fIncohJpsiToMuC->GetXaxis()->FindBin(0.2499999) );
  Double_t intun = fHighPtTailC   ->Integral( 1/*fHighPtTailC   ->GetXaxis()->FindBin(0.000001)*/, fHighPtTailC   ->GetXaxis()->FindBin(0.2499999) );
  Double_t f_I   = (int1i + intun ) / int1c;

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.17,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f ", int1i, intun, int1c, f_I));
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  // latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  /* - Chi square computation.
     -
   */
  // latex->DrawLatex(0.55,0.68,Form( "  #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
  //                                  FitPtDistr->GetChisquare(),
  //                                  FitPtDistr->GetNDF(),
  //                                  FitPtDistr->GetChisquare()/FitPtDistr->GetNDF()
  //                                 )
  //                                );

  TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  // l->AddEntry(  fDimuonPtDistributionDataH, "ALICE data 2018");
  if      ( selectionFlag == 0 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2");
  else if ( selectionFlag == 1 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, 0N0N");
  else if ( selectionFlag == 2 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, 0NXN");
  else if ( selectionFlag == 3 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, XN0N");
  else if ( selectionFlag == 4 ) l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2, XNXN");
  else                           l->AddEntry(  fDimuonPtDistributionDataH, "ALICE Run 2");
  l->AddEntry(  FitPtDistr, Form( "Fit: #chi^{2}/NDF = %.2f / %.2d = %.2f  ",
                                   FitPtDistr->GetChisquare(),
                                   FitPtDistr->GetNDF(),
                                   FitPtDistr->GetChisquare()/FitPtDistr->GetNDF()
                                   )
                                  );
  l->AddEntry(  fCohJpsiToMuC,        "Coherent   J/#psi");
  l->AddEntry(  fIncohJpsiToMuC,      "Incoherent J/#psi");
  l->AddEntry(  fHighPtTailC,         "Incoherent dissociative J/#psi");
  l->AddEntry(  fCohPsi2sToMuPiC,     "Coherent   #psi(2S) feeddown");
  l->AddEntry(  fIncohPsi2sToMuPiC,   "Incoherent #psi(2S) feeddown");
  l->AddEntry(  fTwoGammaToMuMediumC, "Continuum  #gamma#gamma #rightarrow #mu#mu");
  l->Draw();

  // gPad->SaveAs("pngResults/fitPtDistr.png", "RECREATE");
  if        ( selectionFlag == 0 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");
  } else if ( selectionFlag == 1 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistr0N0N.png",   "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistr0N0N_0.png", "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistr0N0N_1.png", "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistr0N0N_2.png", "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistr0N0N.png",   "RECREATE");
  } else if ( selectionFlag == 2 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistr0NXN.png",   "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistr0NXN_0.png", "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistr0NXN_1.png", "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistr0NXN_2.png", "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistr0NXN.png",   "RECREATE");
  } else if ( selectionFlag == 3 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistrXN0N.png",   "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistrXN0N_0.png", "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistrXN0N_1.png", "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistrXN0N_2.png", "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistrXN0N.png",   "RECREATE");
  } else if ( selectionFlag == 4 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistrXNXN.png",   "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistrXNXN_0.png", "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistrXNXN_1.png", "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistrXNXN_2.png", "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistrXN0N.png",   "RECREATE");
  } else                               gPad->SaveAs("pngResults/fitPtDistrALL.png",    "RECREATE");


}
