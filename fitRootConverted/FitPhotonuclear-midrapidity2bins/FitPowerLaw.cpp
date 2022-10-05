// #include "Riostream.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVectorD.h"

// ==========================
// Global data
// --------------------------
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsY;
std::vector< Double_t > values;
std::vector< Double_t > errors;



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
//_____________________________________________________________________________
inline void alice::fcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
// chi2 function to fit data according to prescription in
// Eur.Phys.J. C63 (2009) 625-678, section 9, eq (31)

  // ch2 = sum_i [m_i-mu_i-Sij]^2/D + Sbj
  // Sij = sum_j g_ij*m_i*b_j
  // Sbj = sum_j (b_j)^2
  // D = d_i_stat^2*mu_i*(m_i-Sij)+(d_i,unc*m_i)^2
  // mu_i is the measurement
  // g_ij relative normalization uncertainty at point i from source j
  // d_i = relative uncertainty (either stat or uncorr)

  //Double_t bj[nCorA] = {par[2],par[3],par[4],par[5]};
  Double_t bj[nCorA];
  for(Int_t i=0; i<nCorA; i++) {
    bj[i] = par[i+2];
  }

  Double_t Sbj = 0;
  for(Int_t j=0;j<nCorA;j++) Sbj += (bj[j]*bj[j]);
  Double_t chi2 = 0;
  for (Int_t i=0; i<ALICE_2013_n;i++) {
    Double_t mu_i = ALICE_2013_sig[i];
    Double_t m_i = par[0]*TMath::Power(ALICE_2013_W[i]/90,par[1]);
    Double_t Sij = 0;
    for(Int_t j=0;j<nCorA;j++) Sij += m_i*bj[j]*ALICE_2013_corr_A[j][i];
    Double_t d_stat_i =  ALICE_2013_sta[i]/ALICE_2013_sig[i];
    Double_t d_unc_i = ALICE_2013_unc[i];
    Double_t D = (d_stat_i*d_stat_i*mu_i*(m_i-Sij))+((d_unc_i*m_i)*(d_unc_i*m_i));
    chi2 += ((m_i-mu_i-Sij)*(m_i-mu_i-Sij)/D);
  }
  chi2 += Sbj;
  f = chi2;
}
//______________________________________________
void alice::makeFit()
{

  cout << endl<< endl<< endl;
  if (opt == 0)   cout << " - o - o - o - o - o - o - o - o -   Model A   - o - o - o - o - o - o - o - o -" << endl;
  cout << endl<< endl<< endl;

  // initialize minuit with a maximum of parameters
  Int_t nPar = 0;
  if (opt == 0) nPar = 2+nCorA;
  TMinuit *myMinuit = new TMinuit(nPar);

  // set the function with the minimization process
  if (opt == 0) myMinuit->SetFCN(ifcnChi2ModelA);

  // define parameters
  myMinuit->DefineParameter(0,"N",70,1,0.0,1000.);
  myMinuit->DefineParameter(1,"#delta",0.7,0.1,0.1,10.);
  Char_t ParName[120];
  Int_t nCor = 0;
  if (opt == 0) nCor = nCorA;
  for(Int_t j=0;j<nCor;j++) {
    sprintf(ParName,"b_{%d}",j);
    myMinuit->DefineParameter(2+j,ParName,0.001,0.0001,-1,1.);
  }
  // migrad
  myMinuit->SetMaxIterations(500);
  myMinuit->Migrad();

  ///////////////////////////////////////////////////////////////////////////////
  // get results
  ///////////////////////////////////////////////////////////////////////////////
  Double_t Cov[nPar*nPar];
  myMinuit->mnemat(Cov,nPar);
  Double_t N = 0;  // normalization
  Double_t Nerr = 0;  // normalization error
  myMinuit->GetParameter(0,N,Nerr);
  Double_t D = 0;  // exponent
  Double_t Derr = 0;  // exponent error
  myMinuit->GetParameter(1,D,Derr);
  Double_t covNN = Cov[0];
  Double_t covDD = Cov[nPar+1];
  Double_t covND = Cov[1];
  Double_t rho = covND/TMath::Sqrt(covDD*covNN);
  //status of the fit
  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  myMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

  cout << "#####  Fit results  #####" << endl;
  cout << " N     = " << fixed << setprecision(2) << N << " +/- " << Nerr << endl;
  cout << " delta = " << D << " +/- " << Derr << endl;
  cout << " rho   = " << fixed << setprecision(2) << rho << endl;
  cout << " chi2  = " << fixed << setprecision(2) << fmin << endl;
  cout << " stat: " << istat << endl;

  delete myMinuit;

  Nval = N;
  Dval = D;

  fitGraph = GetShade(N,D,covNN,covDD,covND);
}
//______________________________________________
void ifcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //interface to chi2 function

  alice::instance()->fcnChi2ModelA(npar, gin, f, par, iflag);

}//interface to chi2 function
//______________________________________________
void FitPowerLaw(){
  TFile* SavingFile = new TFile("photonuclear.root");
  TH1F* photonuclear = (TH1F*) SavingFile->Get("histo");
  new TCanvas;
  BeautifyPad();
  gPad->SetLogx();
  gPad->SetLogy();
  photonuclear->GetYaxis()->SetRangeUser(0.002,0.2);
  photonuclear->GetXaxis()->SetRangeUser(10.,1000.);
  photonuclear->Draw();
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
  latex5->DrawLatex(0.31,0.78,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.7,"This thesis");
  latex5->DrawLatex(0.31,0.62,"Coherent J/#psi");




  for (Int_t i = 1; i < photonuclear->GetNbinsX(); i++) {
    if(photonuclear->GetBinContent(i) > 0.001){
      coordsX.push_back(photonuclear->GetBinCenter(i));
      coordsY.push_back(photonuclear->GetBinContent(i));
      errors.push_back(photonuclear->GetBinError(i));
    }
  }


}
