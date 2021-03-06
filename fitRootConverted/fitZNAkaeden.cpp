#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
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
/* - Fit function for the ZNC plots.
 * - I am using simple ROOT to make gaussian fits to the plot.
 */
Double_t fGaussiansTogether(Double_t* x,Double_t* par)
{
  /* - Par 0, 3, 6, 9 :   normalizations.
     - Par 1, 4, 7, 10:   means.
     - Par 2, 5, 8, 11:   widths.
     -
   */
  Double_t ZerothPeak = par[0]*TMath::Exp( -(x[0]-par[1])*(x[0]-par[1])   / (2*par[2]*par[2])   );
  Double_t FirstPeak  = par[3]*TMath::Exp( -(x[0]-par[4])*(x[0]-par[4])   / (2*par[5]*par[5])   );
  Double_t SecondPeak = par[6]*TMath::Exp( -(x[0]-par[7])*(x[0]-par[7])   / (2*par[8]*par[8])   );
  Double_t ThirdPeak  = par[9]*TMath::Exp( -(x[0]-par[10])*(x[0]-par[10]) / (2*par[11]*par[11]) );

  return ( ZerothPeak + FirstPeak + SecondPeak + ThirdPeak );
}
//_____________________________________________________________________________
/* - Single Components of the 4 Gaussians fits.
   - To change from one to the other a "set parameter" on every parameter
   - is highly needed...
   -
 */
Double_t fSingleGaussian(Double_t* x,Double_t* par)
{
  Double_t SinglePeak = par[0]*TMath::Exp( -(x[0]-par[1])*(x[0]-par[1])   / (2*par[2]*par[2])   );
  return SinglePeak;
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitZNA(const char* AnalysisName){
  /* - There are three cases for the selectionFlag:
     - 1) = 0 ; this implies the traditional pt-integrated plot;
     - 2) = 1 ; this is instead the coherent component;
     - 3) = 2 ; this is the incoherent component;
     - 4) = 3 ; ******************* ;
     -
   */
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir->GetObject("ADcheck_", listings);
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
  TH1F *fZNAEnergyAgainstEntriesH = (TH1F*)listings->FindObject("fZNAEnergyAgainstEntriesH");
  fZNAEnergyAgainstEntriesH->Rebin(100);
  // fZNAEnergyAgainstEntriesH->Draw();


  TF1* FourGaussians = new TF1( "FourGaussians",
                                fGaussiansTogether,
                                -10000, 10000,
                                12
                               );
  FourGaussians->SetNpx(1000);
  /* - Setting lhe limits for the many parameters of the fit.
     - This values are based on what had already been used for the
     - RooFit case.
     -
   */
  // ZEROTH PEAK
  FourGaussians->SetParameter(0, 1);
  FourGaussians->SetParameter(1, 0);
  FourGaussians->SetParLimits(1, -30, +30);
  FourGaussians->SetParameter(2, 400);
  FourGaussians->SetParLimits(2, 200, 800);
  // FIRST PEAK
  FourGaussians->SetParameter(3, 1);
  FourGaussians->SetParameter(4, 2300);
  FourGaussians->SetParLimits(4, 2000, 2600);
  FourGaussians->SetParameter(5, 500);
  FourGaussians->SetParLimits(5, 200, 800);
  // SECOND PEAK
  FourGaussians->SetParameter(6,1);
  FourGaussians->SetParameter(7, 2*(FourGaussians->GetParameter(4)) );
  FourGaussians->SetParLimits(7, 4500, 5500);
  FourGaussians->SetParameter(8, TMath::Sqrt( 2*( FourGaussians->GetParameter(5)*FourGaussians->GetParameter(5) -
                                                  FourGaussians->GetParameter(2)*FourGaussians->GetParameter(2) )
                                              +   FourGaussians->GetParameter(2)*FourGaussians->GetParameter(2)
                                             )
                                    );
  FourGaussians->SetParLimits(8, 200, 1000);
  // THIRD PEAK
  FourGaussians->SetParameter(9,  1);
  FourGaussians->SetParameter(10, 3*(FourGaussians->GetParameter(4)) );
  FourGaussians->SetParLimits(10, 7000, 8000);
  FourGaussians->SetParameter(11, TMath::Sqrt( 3*( FourGaussians->GetParameter(5)*FourGaussians->GetParameter(5) -
                                                   FourGaussians->GetParameter(2)*FourGaussians->GetParameter(2) )
                                               +   FourGaussians->GetParameter(2)*FourGaussians->GetParameter(2)
                                              )
                                    );
  FourGaussians->SetParLimits(11, 200, 1000);
  fZNAEnergyAgainstEntriesH->SetLineColor(kBlue);
  fZNAEnergyAgainstEntriesH->SetLineStyle(kSolid);
  fZNAEnergyAgainstEntriesH->SetLineWidth(3);
  fZNAEnergyAgainstEntriesH->SetMarkerStyle(kFullCircle);
  fZNAEnergyAgainstEntriesH->SetMarkerSize(1);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetTitle("ZNA Energy [a.u.]");
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitle( Form( "Counts / (%.0f a.u.)",
                                                         fZNAEnergyAgainstEntriesH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fZNAEnergyAgainstEntriesH->SetTitle("");
  fZNAEnergyAgainstEntriesH->Fit( FourGaussians,"","", -2000., 9000 );
  TCanvas* ZNAEnergy = new TCanvas( "ZNAEnergy", "ZNAEnergy", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);




  /* - Beautifying is starting now.
     -
   */
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetTitleOffset(1.25);
  // fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitleOffset(1.25);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitleOffset(1.45);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetTitleSize(0.045);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitleSize(0.045);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetLabelSize(0.045);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetLabelSize(0.045);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetTitleFont(42);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitleFont(42);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetLabelFont(42);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetLabelFont(42);
  fZNAEnergyAgainstEntriesH->GetXaxis()->SetNdivisions(408);
  fZNAEnergyAgainstEntriesH->GetYaxis()->SetRangeUser(5, fZNAEnergyAgainstEntriesH->GetMaximum()*10.);
  gPad ->SetLogy();
  fZNAEnergyAgainstEntriesH->Draw("SAME");
  /* - Plotting the single Components of the fits in the same TCanvas.
     -
   */
  TF1* SingleGaussian = new TF1( "SingleGaussian",
                                 fSingleGaussian,
                                 -10000, 10000,
                                 3
                                );
  SingleGaussian->SetNpx(1000);
  // DRAWING the ZEROTH PEAK
  SingleGaussian->SetParameter(0, FourGaussians->GetParameter(0));
  SingleGaussian->SetParameter(1, FourGaussians->GetParameter(1));
  SingleGaussian->SetParameter(2, FourGaussians->GetParameter(2));
  SingleGaussian->SetLineColor(kMagenta);
  SingleGaussian->SetLineStyle(kSolid);
  SingleGaussian->SetLineWidth(3);
  SingleGaussian->Draw("ehistsame");
  Double_t ZerothPeakNorm = SingleGaussian->Integral(-2000,9000);
  cout << "ZerothPeakNorm " << ZerothPeakNorm << endl;
  // SingleGaussian->Print();
  // DRAWING the FIRST PEAK
  SingleGaussian->SetParameter(0, FourGaussians->GetParameter(3));
  SingleGaussian->SetParameter(1, FourGaussians->GetParameter(4));
  SingleGaussian->SetParameter(2, FourGaussians->GetParameter(5));
  SingleGaussian->SetLineColor(kRed);
  SingleGaussian->SetLineStyle(kSolid);
  SingleGaussian->SetLineWidth(3);
  SingleGaussian->Draw("SAME");
  Double_t FirstPeakNorm  = SingleGaussian->Integral(-2000,9000);
  cout << "FirstPeakNorm " << FirstPeakNorm << endl;
  // DRAWING the ZEROTH PEAK
  SingleGaussian->SetParameter(0, FourGaussians->GetParameter(6));
  SingleGaussian->SetParameter(1, FourGaussians->GetParameter(7));
  SingleGaussian->SetParameter(2, FourGaussians->GetParameter(8));
  SingleGaussian->SetLineColor(kGreen+1);
  SingleGaussian->SetLineStyle(kDashed);
  SingleGaussian->SetLineWidth(3);
  SingleGaussian->Draw("SAME");
  Double_t SecondPeakNorm = SingleGaussian->Integral(-2000,9000);
  cout << "SecondPeakNorm " << SecondPeakNorm << endl;
  // DRAWING the ZEROTH PEAK
  SingleGaussian->SetParameter(0, FourGaussians->GetParameter(9));
  SingleGaussian->SetParameter(1, FourGaussians->GetParameter(10));
  SingleGaussian->SetParameter(2, FourGaussians->GetParameter(11));
  SingleGaussian->SetLineColor(kGreen+2);
  SingleGaussian->SetLineStyle(kDashed);
  SingleGaussian->SetLineWidth(3);
  SingleGaussian->Draw("SAME");
  Double_t ThirdPeakNorm  = SingleGaussian->Integral(-2000,9000);
  cout << "ThirdPeakNorm " << ThirdPeakNorm << endl;

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex->DrawLatex(0.55,0.72,Form("N_{0 neutrons} = %.0f #pm %.0f", FourGaussians->GetParameter(0), FourGaussians->GetParError(0)));
  latex->DrawLatex(0.55,0.66,Form("N_{1 neutrons} = %.0f #pm %.0f", FourGaussians->GetParameter(3), FourGaussians->GetParError(3)));
  latex->DrawLatex(0.55,0.60,Form("N_{2 neutrons} = %.0f #pm %.0f", FourGaussians->GetParameter(6), FourGaussians->GetParError(6)));
  latex->DrawLatex(0.55,0.54,Form("N_{3 neutrons} = %.0f #pm %.0f", FourGaussians->GetParameter(9), FourGaussians->GetParError(9)));
  // latex->DrawLatex(0.55,0.58,Form("m_{0 neutrons} = %.0f #pm %.0f", mean0Neutrons.getVal(), mean0Neutrons.getError()));
  // latex->DrawLatex(0.55,0.54,Form("m_{1 neutron } = %.0f #pm %.0f",  mean1Neutron.getVal(),  mean1Neutron.getError()));
  // latex->DrawLatex(0.55,0.50,Form("m_{2 neutrons} = %.0f #pm %.0f", mean2Neutrons.getVal(), mean2Neutrons.getError()));
  // latex->DrawLatex(0.55,0.46,Form("m_{3 neutrons} = %.0f #pm %.0f", mean3Neutrons.getVal(), mean3Neutrons.getError()));
  // latex->DrawLatex(0.55,0.42,Form("s_{0 neutrons} = %.0f #pm %.0f", sigma0Neutrons.getVal(), sigma0Neutrons.getError()));
  // latex->DrawLatex(0.55,0.38,Form("s_{1 neutron } = %.0f #pm %.0f",  sigma1Neutron.getVal(),  sigma1Neutron.getError()));
  // latex->DrawLatex(0.55,0.34,Form("s_{2 neutrons} = %.0f #pm %.0f", sigma2Neutrons.getVal(), sigma2Neutrons.getError()));
  // latex->DrawLatex(0.55,0.30,Form("s_{3 neutrons} = %.0f #pm %.0f", sigma3Neutrons.getVal(), sigma3Neutrons.getError()));
  /* - Chi square computation.
     -
   */
  latex->DrawLatex(0.55,0.18,Form( "      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                   FourGaussians->GetChisquare(),
                                   FourGaussians->GetNDF(),
                                   FourGaussians->GetChisquare()/FourGaussians->GetNDF()
                                  )
                                 );
  gPad->SaveAs("pngResults/fZNAEnergyKaeden.png", "RECREATE");





}
