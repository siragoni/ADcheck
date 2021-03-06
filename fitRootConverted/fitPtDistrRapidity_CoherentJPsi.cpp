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

TH1F* fCohJpsiToMuOrig;
TH1F* fCohPsi2sToMuOrig;
TH1F* fCohPsi2sToMuPiOrig;
TH1F* fIncohJpsiToMuOrig;
TH1F* fIncohPsi2sToMuOrig;
TH1F* fIncohPsi2sToMuPiOrig;
TH1F* fTwoGammaToMuMediumOrig;
TH1F* fTwoGammaToMuHighOrig;

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
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 1 ) {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_0");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_0");
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 2 ) {
    // fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionRapidityH_1");
    // fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionRapidityH_1");
    fCohJpsiToMuOrig        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuOrig       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
    fCohPsi2sToMuPiOrig     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
    fIncohJpsiToMuOrig      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuOrig     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
    fIncohPsi2sToMuPiOrig   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuMediumOrig = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
    fTwoGammaToMuHighOrig   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  } else if ( selectionFlag2 == 3 ) {
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
  fCohJpsiToMuOrig        -> Rebin(5);  // -> 25 MeV bins
  fCohPsi2sToMuOrig       -> Rebin(5);
  fCohPsi2sToMuPiOrig     -> Rebin(5);
  fIncohJpsiToMuOrig      -> Rebin(5);
  fIncohPsi2sToMuOrig     -> Rebin(5);
  fIncohPsi2sToMuPiOrig   -> Rebin(5);
  fTwoGammaToMuMediumOrig -> Rebin(5);
  fTwoGammaToMuHighOrig   -> Rebin(5);


  fCohJpsiToMu        = (TH1F*) fCohJpsiToMuOrig        -> Clone("fCohJpsiToMu");
  fCohPsi2sToMu       = (TH1F*) fCohPsi2sToMuOrig       -> Clone("fCohPsi2sToMu");
  fCohPsi2sToMuPi     = (TH1F*) fCohPsi2sToMuPiOrig     -> Clone("fCohPsi2sToMuPi");
  fIncohJpsiToMu      = (TH1F*) fIncohJpsiToMuOrig      -> Clone("fIncohJpsiToMu");
  fIncohPsi2sToMu     = (TH1F*) fIncohPsi2sToMuOrig     -> Clone("fIncohPsi2sToMu");
  fIncohPsi2sToMuPi   = (TH1F*) fIncohPsi2sToMuPiOrig   -> Clone("fIncohPsi2sToMuPi");
  fTwoGammaToMuMedium = (TH1F*) fTwoGammaToMuMediumOrig -> Clone("fTwoGammaToMuMedium");
  fTwoGammaToMuHigh   = (TH1F*) fTwoGammaToMuHighOrig   -> Clone("fTwoGammaToMuHigh");



  // Float_t PtBins[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
  //                         0.200, 0.225, 0.250, 0.275, 0.350, 0.425, 0.500, 0.575,
  //                         0.650, 0.725,
  //                         0.800, 0.875, 0.950, 1.100, 1.250, 1.400, 1.600, 1.800,
  //                         2.000, 2.500, 3.000, 3.500, 4.000, 5.000
  //                       };
  // Int_t   PtBinNumber = sizeof(PtBins)/sizeof(Float_t) - 1; // or just = 9
  // fCohJpsiToMu        = new TH1F("fCohJpsiToMu",        "fCohJpsiToMu",        PtBinNumber, PtBins );
  // fCohPsi2sToMu       = new TH1F("fCohPsi2sToMu",       "fCohPsi2sToMu",       PtBinNumber, PtBins );
  // fCohPsi2sToMuPi     = new TH1F("fCohPsi2sToMuPi",     "fCohPsi2sToMuPi",     PtBinNumber, PtBins );
  // fIncohJpsiToMu      = new TH1F("fIncohJpsiToMu",      "fIncohJpsiToMu",      PtBinNumber, PtBins );
  // fIncohPsi2sToMu     = new TH1F("fIncohPsi2sToMu",     "fIncohPsi2sToMu",     PtBinNumber, PtBins );
  // fIncohPsi2sToMuPi   = new TH1F("fIncohPsi2sToMuPi",   "fIncohPsi2sToMuPi",   PtBinNumber, PtBins );
  // fTwoGammaToMuMedium = new TH1F("fTwoGammaToMuMedium", "fTwoGammaToMuMedium", PtBinNumber, PtBins );
  // fTwoGammaToMuHigh   = new TH1F("fTwoGammaToMuHigh",   "fTwoGammaToMuHigh",   PtBinNumber, PtBins );
  //
  // Double_t BinCenter = 0;
  // for ( Int_t ibin = 1; ibin < fCohJpsiToMuOrig->GetNbinsX(); ibin++ ) {
  //   BinCenter = ((TAxis*)fCohJpsiToMuOrig->GetXaxis())->GetBinCenter(ibin);
  //   if ( BinCenter > PtBins[PtBinNumber-1] ) continue;
  //   cout << "BinCenter = " << BinCenter << endl;
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fCohJpsiToMu->Fill( fCohJpsiToMuOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fCohJpsiToMu->SetBinContent(ibinVariable+1, fCohJpsiToMu->GetBinContent(ibinVariable+1) + (fCohJpsiToMuOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fCohPsi2sToMu->Fill( fCohPsi2sToMuOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fCohPsi2sToMu->SetBinContent(ibinVariable+1, fCohPsi2sToMu->GetBinContent(ibinVariable+1) + (fCohPsi2sToMuOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fCohPsi2sToMuPi->Fill( fCohPsi2sToMuPiOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fCohPsi2sToMuPi->SetBinContent(ibinVariable+1, fCohPsi2sToMuPi->GetBinContent(ibinVariable+1) + (fCohPsi2sToMuPiOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fIncohJpsiToMu->Fill( fIncohJpsiToMuOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fIncohJpsiToMu->SetBinContent(ibinVariable+1, fIncohJpsiToMu->GetBinContent(ibinVariable+1) + (fIncohJpsiToMuOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fIncohPsi2sToMu->Fill( fIncohPsi2sToMuOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fIncohPsi2sToMu->SetBinContent(ibinVariable+1, fIncohPsi2sToMu->GetBinContent(ibinVariable+1) + (fIncohPsi2sToMuOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fIncohPsi2sToMuPi->Fill( fIncohPsi2sToMuPiOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fIncohPsi2sToMuPi->SetBinContent(ibinVariable+1, fIncohPsi2sToMuPi->GetBinContent(ibinVariable+1) + (fIncohPsi2sToMuPiOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fTwoGammaToMuMedium->Fill( fTwoGammaToMuMediumOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fTwoGammaToMuMedium->SetBinContent(ibinVariable+1, fTwoGammaToMuMedium->GetBinContent(ibinVariable+1) + (fTwoGammaToMuMediumOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  //   for( Int_t ibinVariable = 0; ibinVariable < PtBinNumber-1; ibinVariable++ ) {
  //     if ( BinCenter < PtBins[ibinVariable+1] ){
  //       // fTwoGammaToMuHigh->Fill( fTwoGammaToMuHighOrig->GetBinContent(ibin), 1./(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       fTwoGammaToMuHigh->SetBinContent(ibinVariable+1, fTwoGammaToMuHigh->GetBinContent(ibinVariable+1) + (fTwoGammaToMuHighOrig->GetBinContent(ibin))/(PtBins[ibinVariable+1]-PtBins[ibinVariable]) );
  //       break;
  //     }
  //   }
  // }
  // new TCanvas;
  // fCohJpsiToMu        -> Draw();
  // new TCanvas;
  // fCohPsi2sToMu        -> Draw();
  // new TCanvas;
  // fCohPsi2sToMuPi        -> Draw();
  // new TCanvas;
  // fIncohJpsiToMu        -> Draw();
  // new TCanvas;
  // fIncohPsi2sToMu        -> Draw();
  // new TCanvas;
  // fIncohPsi2sToMuPi        -> Draw();
  // new TCanvas;
  // fTwoGammaToMuMedium        -> Draw();
  // new TCanvas;
  // fTwoGammaToMuHigh        -> Draw();
  // new TCanvas;

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
  TF1* fModelForHighPtTail = new TF1("fModelForHighPtTail","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0.,4);
  fModelForHighPtTail->SetParameter(0,1);
//  fModelForHighPtTail->SetParameter(1,debug==4 ? 1.25 : 1.);
//  fModelForHighPtTail->SetParameter(2,debug==4 ? 6.1 : 1.);
  fModelForHighPtTail->SetParameter(1, 1.6/*1.79*/);
  fModelForHighPtTail->SetParameter(2, 3.58);
  // fModelForHighPtTail->SetParameter(2, 2.);
  // fModelForHighPtTail->SetNpx( fCohJpsiToMu->GetNbinsX() );
  fModelForHighPtTail->SetNpx( 4.0/0.025 );
  fHighPtTail = (TH1F*) fModelForHighPtTail->GetHistogram()->Clone("fHighPtTail");
  for (Int_t ibin=1; ibin<=fHighPtTail->GetNbinsX(); ibin++) {
    fHighPtTail->SetBinError(ibin,0);
  }
  Double_t Integral_fHighPtTail = fHighPtTail->Integral();
  fHighPtTail->Scale( 1/Integral_fHighPtTail );
  // Float_t PtBins2[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
  //                         0.200, 0.225, 0.250, 0.275, 0.350, 0.425, 0.500, 0.575,
  //                         0.650, 0.725,
  //                         0.800, 0.875, 0.950, 1.100, 1.250, 1.400, 1.600, 1.800,
  //                         2.000, 2.500, 3.000, 3.500, 4.000, 5.000
  //                       };
  // Int_t   PtBinNumber2 = sizeof(PtBins2)/sizeof(Float_t) - 1; // or just = 9
  // fHighPtTail = new TH1F( "fHighPtTail", "fHighPtTail", PtBinNumber2, PtBins2 );
  // for (Int_t ibin=1; ibin<=fHighPtTail->GetNbinsX(); ibin++) {
  //   fHighPtTail->SetBinError(ibin,0);
  //   fHighPtTail->SetBinContent(ibin,  fModelForHighPtTail->Integral(PtBins2[ibin-1], PtBins2[ibin])/(PtBins2[ibin]-PtBins2[ibin-1]));
  // }
  // Double_t Integral_fHighPtTail = fHighPtTail->Integral();
  // fHighPtTail->Scale( 1/Integral_fHighPtTail );



  TH1F *fDimuonPtDistributionDataH = 0x0;
  if      ( selectionFlag == 0 ) {
       if      ( selectionFlag2 == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
       else if ( selectionFlag2 == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_1");
       else if ( selectionFlag2 == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_2");
       else if ( selectionFlag2 == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_3");
       else if ( selectionFlag2 == 4 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_4");
       else if ( selectionFlag2 == 5 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_5");
       else                            fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionRapidityH_0");
  }
  else                                 fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  fDimuonPtDistributionDataH->Rebin(5);
  fDimuonPtDistributionDataH->Draw("PE");


  /* - Preliminary fit to the unknown contribution.
   * -
   */
  TF1* fModelForHighPtTailPreliminary = new TF1("fModelForHighPtTailPreliminary","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0.0000000001,5.);
  fModelForHighPtTailPreliminary->SetParameter(0,1);
  // fModelForHighPtTailPreliminary->SetParameter(1,debug==4 ? 1.25 : 1.);
  // fModelForHighPtTailPreliminary->SetParameter(2,debug==4 ? 6.1 : 1.);
  // fModelForHighPtTailPreliminary->SetParameter(1, 1.6);
  fModelForHighPtTailPreliminary->SetParameter(1, 1.79);
  fModelForHighPtTailPreliminary->SetParameter(2, 3.58);
  // fModelForHighPtTailPreliminary->SetParameter(2, 6.1);
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
  // UtilityCanvas->SaveAs("UtilityCanvas_v3.png");


  new TCanvas;
  TF1* FitPtDistr = new TF1(  "FitPtDistr",
                              fPtDistr,
                              0, 3,
                              /*7*/5
                              );
  FitPtDistr->SetNpx(1000);




  Double_t kFeedDownCoherent = 0.05;     // neutral element
  Double_t kError            = 0.01;
  // if ( selectionFlag  == 1 ) {
  //   FitPtDistr->FixParameter(4, kFeedDownCoherent * 2);
  // } else {
    FitPtDistr->FixParameter(4, kFeedDownCoherent);
  // }
  // FitPtDistr->SetParameter(4, kFeedDownCoherent);
  // FitPtDistr->SetParLimits(4, FitPtDistr->GetParameter(4)*(1-kError), FitPtDistr->GetParameter(4)*(1+kError));


  // // FitPtDistr->SetParLimits(0, 0.00000000000000001, 99999999999);
  // // FitPtDistr->SetParLimits(1, 0.00000000000000001, 99999999999);
  // if ( selectionFlag == 1 || selectionFlag == 3 ) {
  //   FitPtDistr->SetParLimits(0, 10, 10000);
  // } else {
  //   FitPtDistr->SetParLimits(0, 1, 1000);
  // }
  // // FitPtDistr->SetParLimits(1, 0.1, 10000);



  if        ( selectionFlag  == 0 ) { // COHERENT analysis
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(0, 11716*0.043 );
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(0, 11716*0.144 );
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(0, 11716*0.271 );
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(0, 11716*0.293 );
    else if ( selectionFlag2 == 4 ) FitPtDistr->SetParameter(0, 11716*0.194 );
    else if ( selectionFlag2 == 5 ) FitPtDistr->SetParameter(0, 11716*0.051 );
  } else if ( selectionFlag  == 1 ) { // 0N0N
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(0, 9450 );
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(0, 1856 );
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(0, 5336 );
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(0, 2230 );
    else                            FitPtDistr->SetParameter(0, 9450 );
  } else if ( selectionFlag  == 2 ) { // 0NXN
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(0,  864 );
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(0,  135 );
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(0,  484 );
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(0,  232 );
    else                            FitPtDistr->SetParameter(0,  864 );
  } else if ( selectionFlag  == 3 ) {
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(0,  810 );
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(0,   97 );
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(0,  430 );
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(0,  265 );
    else                            FitPtDistr->SetParameter(0,  810 );
  } else if ( selectionFlag  == 4 ) {
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(0,  592 );
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(0,   65 );
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(0,  316 );
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(0,  187 );
    else                            FitPtDistr->SetParameter(0,  592 );
  }

  FitPtDistr->SetParLimits(0, FitPtDistr->GetParameter(0)*0.8, FitPtDistr->GetParameter(0)*1.4);

  if        ( selectionFlag  == 1 ) { // 0N0N
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1, 110, 300);
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,  10,  90);
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,  50, 150);
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,  25,  80);
    else                            FitPtDistr->SetParLimits(1, 0.1, 100);
  } else if ( selectionFlag  == 2 ) { // 0NXN
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1,  70, 100);
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,   1,  30);
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,  30,  70);
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,  10,  50);
    else                            FitPtDistr->SetParLimits(1, 0.1, 100);
  } else if ( selectionFlag  == 3 ) {
    // if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1,  90, 150);
    // else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,  40,  70);
    // else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,  70, 1000);
    // else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,  25,  45);
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1,  90, 800);
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,  40, 300);
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,  70,1000);
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,  25, 150);
    else                            FitPtDistr->SetParLimits(1, 0.1, 100);
  } else if ( selectionFlag  == 4 ) {
    // if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1,   1,  40);
    // else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,   1,  40);
    // else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,   1,  40);
    // else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,   1,  40);
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParLimits(1,   1, 120);
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParLimits(1,   1,  60);
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParLimits(1,   1,  80);
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParLimits(1,   1,  60);
    else                            FitPtDistr->SetParLimits(1, 0.1, 100);
  }




  // FEED-DOWN COHERENT
  // Double_t kFeedDownIncoherent = 1705;     // neutral element
  // Double_t kErrorInc           = 0.175;
  // FitPtDistr->SetParameter(5, kFeedDownIncoherent);
  // FitPtDistr->SetParLimits(5, FitPtDistr->GetParameter(5)*(1-kErrorInc), FitPtDistr->GetParameter(5)*(1+kErrorInc));

  // Gamma+Gamma Medium
  // FitPtDistr->SetParameter(2, 4800);
  // FitPtDistr->SetParameter(2, 6531);
  // if        ( selectionFlag == 0 ) {
  //   // FitPtDistr->SetParameter(2, 6531);  //2018
  //   // FitPtDistr->SetParameter(2, 9035);  //2018+2015 with SPD
  //   FitPtDistr->SetParameter(2, 11000);  //2018+2015 no SPD
  // } else if ( selectionFlag == 1 ) { // 0N0N
  //   // FitPtDistr->SetParameter(2, 5213);  //2018
  //   // FitPtDistr->SetParameter(2, 7754);  //2018+2015 with SPD
  //   if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 6460);          //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 6460 * 0.165);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 6460 * 0.572);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 6460 * 0.263);  //2018 no SPD, no AD
  //   else                            FitPtDistr->SetParameter(2, 6460);          //2018 no SPD, no AD
  // } else if ( selectionFlag == 2 ) { // 0NXN
  //   // FitPtDistr->SetParameter(2,  370);  //2018
  //   // FitPtDistr->SetParameter(2,  439);  //2018+2015 with SPD
  //   if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 512);          //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 512 * 0.151);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 512 * 0.528);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 512 * 0.321);  //2018 no SPD, no AD
  //   else                            FitPtDistr->SetParameter(2, 512);          //2018 no SPD, no AD
  // } else if ( selectionFlag == 3 ) {
  //   // FitPtDistr->SetParameter(2,  470);  //2018
  //   // FitPtDistr->SetParameter(2,  543);  //2018+2015 with SPD
  //   if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 647);          //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 647 * 0.165);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 647 * 0.573);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 647 * 0.262);  //2018 no SPD, no AD
  //   else                            FitPtDistr->SetParameter(2, 647);          //2018 no SPD, no AD
  // } else if ( selectionFlag == 4 ) {
  //   // FitPtDistr->SetParameter(2,  140);  //2018
  //   // FitPtDistr->SetParameter(2,  145);  //2018+2015 with SPD
  //   if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 258);          //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 258 * 0.107);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 258 * 0.613);  //2018 no SPD, no AD
  //   else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 258 * 0.280);  //2018 no SPD, no AD
  //   else                            FitPtDistr->SetParameter(2, 258);          //2018 no SPD, no AD
  // }

  // // FitPtDistr->FixParameter(2, 0);
  // if ( selectionFlag == 1 ) {
  //   FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*0.6, FitPtDistr->GetParameter(2)*1.4);
  // } else {
  //   FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*0.6, FitPtDistr->GetParameter(2)*1.4);
  // }

  if        ( selectionFlag == 0 ) {
    // FitPtDistr->SetParameter(2, 6531);  //2018
    // FitPtDistr->SetParameter(2, 9035);  //2018+2015 with SPD
    // if      ( selectionFlag2 == 0 ) FitPtDistr->FixParameter(2, 7877 * 0.030);
    // else if ( selectionFlag2 == 1 ) FitPtDistr->FixParameter(2, 7877 * 0.132);
    // else if ( selectionFlag2 == 2 ) FitPtDistr->FixParameter(2, 7877 * 0.263);
    // else if ( selectionFlag2 == 3 ) FitPtDistr->FixParameter(2, 7877 * 0.308);
    // else if ( selectionFlag2 == 4 ) FitPtDistr->FixParameter(2, 7877 * 0.211);
    // else if ( selectionFlag2 == 5 ) FitPtDistr->FixParameter(2, 7877 * 0.059);
    // else                            FitPtDistr->FixParameter(2, 7877);
    if      ( selectionFlag2 == 0 ) FitPtDistr->SetParameter(2, 7877 * 0.030);
    else if ( selectionFlag2 == 1 ) FitPtDistr->SetParameter(2, 7877 * 0.132);
    else if ( selectionFlag2 == 2 ) FitPtDistr->SetParameter(2, 7877 * 0.263);
    else if ( selectionFlag2 == 3 ) FitPtDistr->SetParameter(2, 7877 * 0.308);
    else if ( selectionFlag2 == 4 ) FitPtDistr->SetParameter(2, 7877 * 0.211);
    else if ( selectionFlag2 == 5 ) FitPtDistr->SetParameter(2, 7877 * 0.059);
    else                            FitPtDistr->SetParameter(2, 7877);
  }
  FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*0.4, FitPtDistr->GetParameter(2)*1.4);



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
  fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(fDimuonPtDistributionDataH->GetMaximum()*0.0001, fDimuonPtDistributionDataH->GetMaximum()*1000.);
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
  // Double_t int1c = fCohJpsiToMuC  ->Integral( 1/*fCohJpsiToMuC  ->GetXaxis()->FindBin(0.000001)*/, fCohJpsiToMuC  ->GetXaxis()->FindBin(0.2499999) );
  // Double_t int1i = fIncohJpsiToMuC->Integral( 1/*fIncohJpsiToMuC->GetXaxis()->FindBin(0.000001)*/, fCohJpsiToMuC->GetXaxis()->FindBin(0.2499999) );
  // Double_t intun = fHighPtTailC   ->Integral( 1/*fHighPtTailC   ->GetXaxis()->FindBin(0.000001)*/, fCohJpsiToMuC   ->GetXaxis()->FindBin(0.2499999) );
  // Double_t f_I   = (int1i + intun ) / int1c;

  Double_t err1c = 0;
  Double_t err1i = 0;
  Double_t errun = 0;
  Double_t errfI = 0;
  Double_t int1c = fCohJpsiToMuC  ->IntegralAndError( 1, fCohJpsiToMuC->GetXaxis()->FindBin(0.2499999), err1c, "" );
  Double_t int1i = fIncohJpsiToMuC->IntegralAndError( 1, fCohJpsiToMuC->GetXaxis()->FindBin(0.2499999), err1i, "" );
  Double_t intun = fHighPtTailC   ->IntegralAndError( 1, fCohJpsiToMuC->GetXaxis()->FindBin(0.2499999), errun, "" );
  Double_t f_I   = (int1i + intun ) / int1c;

  // err1c *= (int1c * (FitPtDistr->GetParError(0)));
  // err1i *= (int1i * (FitPtDistr->GetParError(1)));
  // errun *= (intun * (FitPtDistr->GetParError(3)));
  err1c += int1c * (FitPtDistr->GetParError(0) / FitPtDistr->GetParameter(0));
  err1i += int1i * (FitPtDistr->GetParError(1) / FitPtDistr->GetParameter(1));
  errun += intun * (FitPtDistr->GetParError(3) / FitPtDistr->GetParameter(3));

  // errfI          = err1c + err1i + errun;
  Double_t RatioErrfIOverfI        = 0;
  Double_t RatioErr1cOver1c        = 0;
  Double_t RatioErrIncohOverIncoh  = 0;
  RatioErr1cOver1c                 = err1c / int1c;
  RatioErrIncohOverIncoh           = ( err1i + errun ) / ( int1i + intun );
  RatioErrfIOverfI                 = RatioErrIncohOverIncoh + RatioErr1cOver1c;
  errfI                            = RatioErrfIOverfI * f_I;
  // errfI          = (RatioErrfIOverfI + (FitPtDistr->GetParError(0) / FitPtDistr->GetParameter(0)) + (FitPtDistr->GetParError(1) / FitPtDistr->GetParameter(1)) + (FitPtDistr->GetParError(3) / FitPtDistr->GetParameter(3))) * f_I;
  cout << "err1c = " << err1c;
  cout << "err1i = " << err1i;
  cout << "errun = " << errun;
  cout << "errfI = " << errfI;

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  // latex->DrawLatex(0.17,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f ", int1i, intun, int1c, f_I));
  latex->DrawLatex(0.17,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f #pm %.7f", int1i, intun, int1c, f_I, errfI));
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

  // gPad->SaveAs("pngResults/fitPtDistr_v3.png", "RECREATE");
  if        ( selectionFlag == 0 )  {
       if      ( selectionFlag2 == 0 ) gPad->SaveAs("pngResults/fitPtDistrALL_0.png",    "RECREATE");
       else if ( selectionFlag2 == 1 ) gPad->SaveAs("pngResults/fitPtDistrALL_1.png",    "RECREATE");
       else if ( selectionFlag2 == 2 ) gPad->SaveAs("pngResults/fitPtDistrALL_2.png",    "RECREATE");
       else if ( selectionFlag2 == 3 ) gPad->SaveAs("pngResults/fitPtDistrALL_3.png",    "RECREATE");
       else if ( selectionFlag2 == 4 ) gPad->SaveAs("pngResults/fitPtDistrALL_4.png",    "RECREATE");
       else if ( selectionFlag2 == 5 ) gPad->SaveAs("pngResults/fitPtDistrALL_5.png",    "RECREATE");
       else                            gPad->SaveAs("pngResults/fitPtDistrALL.png",      "RECREATE");
  } else                               gPad->SaveAs("pngResults/fitPtDistrALL.png",      "RECREATE");






  /* -
   * - Here we compute the Chi Square of the fit.
   * - This is needed to understand which bins
   * - may lead to problems...
   * -
   * - This is especially the case for the XN0N
   * - plots and the 0N0N neutron emission class.
   * -
   * -
   */
  // Double_t ChiSquareSingleContribution[PtBinNumber2];
  // TCanvas* ChiSquareCanvas = new TCanvas( "ChiSquareCanvas", "ChiSquareCanvas", 900, 800 );
  // TH1F* ChiSquareH  = new TH1F( "ChiSquareH" ,"ChiSquareH" , PtBinNumber2, PtBins2);
  // TH1F* ChiSquareH2 = new TH1F( "ChiSquareH2","ChiSquareH2", PtBinNumber2, PtBins2);
  // TH1F* ChiSquareH3 = new TH1F( "ChiSquareH3","ChiSquareH3", PtBinNumber2, PtBins2);
  // // TH1F* fTwoGammaFromModelH = (TH1F*) GammaGammaFit->GetHistogram()->Clone("fTwoGammaFromModelH");
  // // for (Int_t ibin=1; ibin<=fTwoGammaFromModelH->GetNbinsX(); ibin++) {
  // //   fTwoGammaFromModelH->SetBinError(ibin,0);
  // // }
  // Int_t whichBin  = fDimuonPtDistributionDataH->GetXaxis()->FindBin(0.001);
  // Int_t whichBin2 = ChiSquareH                ->GetXaxis()->FindBin(0.001);
  // for( Int_t iLoop = 0; iLoop < PtBinNumber2 - 3; iLoop++ ) {
  //   // Double_t binContent = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t binContent = fDimuonPtDistributionDataH->GetBinContent(           whichBin + iLoop);
  //   Double_t binError   = fDimuonPtDistributionDataH->GetBinError(             whichBin + iLoop);
  //   Double_t binCenter  = fDimuonPtDistributionDataH->GetXaxis()->GetBinCenter(whichBin + iLoop);
  //   // Double_t fitValue   = fFitInvMass->Eval(2.225+iLoop*0.5);
  //   Double_t fitValue   = FitPtDistr->Eval(binCenter);
  //   cout << "binContent = " << binContent << endl;
  //   cout << "fitValue   = " << fitValue   << endl;
  //   // Double_t pull = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5)) - fFitInvMass->Eval(2.225+iLoop*0.5);
  //   Double_t pull = binContent - fitValue;
  //   cout << pull << endl;
  //   // Double_t help = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   // Double_t help = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t pullHelp;
  //   if( pull != 0 && binContent != 0 ) pullHelp = pull /binContent;
  //   else pullHelp = 0;
  //   cout << "pullHelp = " << pullHelp << endl;
  //   Double_t pullHelp2;
  //   if( pull != 0 && binError != 0 ) pullHelp2 = pull /binError;
  //   else pullHelp2 = 0;
  //   cout << "pullHelp2 = " << pullHelp2 << endl;
  //   cout << "binCenter = " << ChiSquareH                ->GetBinCenter( whichBin2 + iLoop ) << endl;
  //   cout << "binCenter = " << fDimuonPtDistributionDataH->GetBinCenter( whichBin  + iLoop ) << endl;
  //
  //   ChiSquareH->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pullHelp
  //                              // pull
  //                              );
  //   ChiSquareH2->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pullHelp2
  //                              // pull
  //                              );
  //   ChiSquareH3->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pull*pull/(binError*PtBinNumber2)
  //                              // pull
  //                              );
  //   ChiSquareSingleContribution[iLoop] = pull*pull/(binError);
  //
  //
  // }
  // Double_t AllChiSquare = 0;
  // for( Int_t i = 0; i < PtBinNumber2 - 3; i++){
  //   cout << "ChiSquareSingleContribution = " << ChiSquareSingleContribution[i] << endl;
  //   AllChiSquare += ChiSquareSingleContribution[i];
  // }
  // cout << "AllChiSquare     = " << AllChiSquare << endl;
  // cout << "AllChiSquare/DOF = " << (AllChiSquare/PtBinNumber2) << endl;
  // ChiSquareH->Draw();
  // TCanvas* ChiSquareCanvas2 = new TCanvas( "ChiSquareCanvas2", "ChiSquareCanvas2", 900, 800 );
  // ChiSquareH2->Draw();
  // TCanvas* ChiSquareCanvas3 = new TCanvas( "ChiSquareCanvas3", "ChiSquareCanvas3", 900, 800 );
  // ChiSquareH3->Draw();


}
