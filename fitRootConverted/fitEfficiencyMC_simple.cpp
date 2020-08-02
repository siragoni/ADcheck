#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"



//_____________________________________________________________________________
/* - Computes the efficiency of the MC as a
 * - simple division.
 * -
 */
void fitEfficiencyMC(){

  // TFile* fileList = new TFile("MCtrainResults/2020-06-26/kCohJpsiToMu/AnalysisResults.root");  // same settings as CheckAD
  // TFile* fileList = new TFile("MCtrainResults/2020-06-26/kCohJpsiToMu/AnalysisResultsCoherentCMUP6.root");  // same settings as CheckAD
  TFile* fileList = new TFile("AnalysisResultsCohLHC16b2.root");  // same settings as CheckAD
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
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
  TH1F* fEfficiencyPerRunH[6]   = { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 };
  TH1F* fMCEfficiencyPerRunH[6] = { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 };
  for (Int_t i = 0; i < 6; i++) {
    // fEfficiencyPerRunH[i]   = (TH1F*)listings->FindObject( Form("fEfficiencyPerRunRapidityH_%d",   i) );
    fEfficiencyPerRunH[i]   = (TH1F*)listings->FindObject( Form("fInvariantMassDistributionCoherentRapidityBinsH_%d",   i) );
    fMCEfficiencyPerRunH[i] = (TH1F*)listings->FindObject( Form("fMCEfficiencyPerRunRapidityH_%d", i) );
    fEfficiencyPerRunH[i]   ->Sumw2();
    fMCEfficiencyPerRunH[i] ->Sumw2();
  }

  Double_t RECnumber[6]   = { 0,0,0,0,0,0 };
  Double_t GENnumber[6]   = { 0,0,0,0,0,0 };
  Double_t efficiency[6]  = { 0,0,0,0,0,0 };
  Double_t RECint = 0;
  Double_t GENint = 0;
  Double_t EFFint = 0;
  for (Int_t i = 0; i < 6; i++) {
    RECnumber[i]  = fEfficiencyPerRunH[i]  ->GetEntries();
    GENnumber[i]  = fMCEfficiencyPerRunH[i]->GetEntries();
    efficiency[i] = RECnumber[i] / GENnumber[i];
    RECint += RECnumber[i];
    GENint += GENnumber[i];
  }

  EFFint = RECint / GENint;


  for (Int_t i = 0; i < 6; i++) {
    cout << "RECnumber[ " << i << "] = " << RECnumber[i]  << endl;
    cout << "GENnumber[ " << i << "] = " << GENnumber[i]  << endl;
    cout << "efficiency[" << i << "] = " << efficiency[i] << endl;
  }


  cout << "RECint = " << RECint  << endl;
  cout << "GENint = " << GENint  << endl;
  cout << "EFFint = " << EFFint  << endl;




}
