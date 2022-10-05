
//global constants with range of the plot
// const Double_t gxmin = 10., gymin = 0.001, gxmax = 1000., gymax = 1*1000.;
const Double_t gxmin = 10., gymin = 0.001, gxmax = 1000., gymax = 2.;

//_____________________________________________________________________________
void SetGraphStyle(TGraph* g, Color_t mcolor, Style_t mstyle, Size_t msize, Color_t lcolor, Style_t lstyle, Width_t lwidth){
  g->SetMarkerColor(mcolor);
  g->SetMarkerSize(msize);
  g->SetMarkerStyle(mstyle);
  g->SetLineColor(lcolor);
  g->SetLineStyle(lstyle);
  g->SetLineWidth(lwidth);
}//SetGraphStyle

void ifcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);//declaration of interface to chi2 function

//_____________________________________________________________________________
class alice {
  Int_t opt;//option for the fit
  // static const Int_t ALICE_photo_n = 4;
  // static const Int_t ALICE_photo_n = 5;
  static const Int_t ALICE_photo_n = 9;
  Double_t ALICE_chi2[ALICE_photo_n]   = { 0,0,0,0,0,0    ,0,0,0, };
  const Double_t ALICE_photo_W[ALICE_photo_n]   = { 19.130700000, 24.5643000, 31.541200000, 493.3870000, 633.5210000, 813.4570000,             97.1536, 124.748, 160.179 };
  // const Double_t ALICE_photo_sig[ALICE_photo_n] = {  0.00882072,  0.0138693,  0.0168523,     0.0473069,   0.0523713,     0.0628193,            0.0218977, 0.0230963,  0.0246617};
  // const Double_t ALICE_photo_sta[ALICE_photo_n]  = {  0.000260183,0.000307948,0.000593475,   0.00775179,   0.00799204,   0.024678,           0.00101229, 0.00136407, 0.00506907};
  // const Double_t ALICE_photo_sys[ALICE_photo_n]  = {  0.00069,    0.00056,    0.00135,       0.00695,   0.00899,         0.01676,           0.00101229, 0.00136407, 0.00506907};
  const Double_t ALICE_photo_sig[ALICE_photo_n]  = {  0.00882072,  0.0138693,  0.0168523,     0.0473069,   0.0523713,      0.0628193,       0.0218977,  0.0246131,  0.0246617};
  const Double_t ALICE_photo_sta[ALICE_photo_n]  = {  0.000260183, 0.000307948,0.000593475,   0.00775179,   0.00799204,    0.024678,        0.00461081, 0.000587666,  0.00708194 };
  const Double_t ALICE_photo_sys[ALICE_photo_n]  = {  0.00069,    0.00056,    0.00135,       0.00695,   0.00899,           0.01334,         0.00175, 0.00141, 0.00485};
  
  
  
  // model A, considers MUON and trigger efficiencyes uncorrelated
  // static const Int_t nUncA = 10; //  sources of uncorrelated error
  // const Double_t ALICE_2013_unc_A[nUncA][ALICE_2013_n] = // relative uncertainty!
  // {
  //   {0, 0, 0, 0.057, 0.012, 0.009, 0.008, 0.033, 0}, // TPC track selection
  //   {0, 0, 0, 0, 0, 0.013, 0.006, 0, 0}, // PID efficiency
  //   {0.04, 0.04, 0.04, 0.02, 0.02, 0, 0, 0.03, 0.06}, // muon track eff
  //   {0.01, 0.01, 0.01, 0.005, 0.005, 0, 0, 0.005, 0.01}, // muon matching eff
  //   {0.028, 0.028, 0.028, 0.02, 0.02, 0, 0 , 0.016, 0.032}, // muon trig eff
  //   {0, 0, 0, 0, 0, 0.101, 0.101, 0, 0}, // central trigger eff
  //   {0, 0, 0, 0, 0, 0, 0, 0.034, 0.035}, // V0C efficiency
  //   {0.066, 0.037, 0.08, 0.036, 0.022, 0.021, 0.019, 0.03, 0.06}, // signal extraction, last number should be -0,+0.06
  //   {0.031, 0.02, 0.02, 0.013, 0.01, 0.01, 0.01, 0.014, 0.031}, // feed-down
  //   //lumi at mid-rapidity by quadratic average:
  //   //sqrt((3.3%*2.1)**2 + (3%*4.8)**2) / (2.1 + 4.8) = 2.3%
  //   {0.033, 0.033, 0.033, 0.033, 0.033, 0.023, 0.023, 0.03, 0.03} // lumi uncorrelated
  // };
  static const Int_t nUncA = 1; //  sources of uncorrelated error
  const Double_t ALICE_photo_unc_A[nUncA][ALICE_photo_n] = // relative uncertainty!
  {
    // {0, 0, 0, 0}
    // {0.00038/ALICE_photo_sig[0],0.00029/ALICE_photo_sig[1],0.00071/ALICE_photo_sig[2],0.00232/ALICE_photo_sig[3],0.00165/ALICE_photo_sig[4],0.00357/ALICE_photo_sig[5], 0, 0,  0}
    {0.00069/ALICE_photo_sig[0],0.00054/ALICE_photo_sig[1],0.00131/ALICE_photo_sig[2],   0.00385/ALICE_photo_sig[3],0.00236/ALICE_photo_sig[4],0.00628/ALICE_photo_sig[5], 0, 0,  0}
    // {0, 0, 0, 0, 0}
  };
  Double_t ALICE_photo_unc[ALICE_photo_n]; // to store total uncorrelated uncertainty
  // static const Int_t nCorA = 3; // sources of correlated error
  // const Double_t ALICE_2013_corr_A[nCorA][ALICE_2013_n] = // relative uncertainty!
  // {
  //   {0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016}, // lumi
  //   {0.006, 0.006, 0.006, 0.006, 0.006, 0.004, 0.004, 0.006, 0.006}, // br
  //   {0.02, 0.02, 0.02, 0.027, 0.035, 0.021, 0.021, 0.012, 0.02}//, // v0a veto
  //   //{0, 0, 0, 0, 0, 0.09, 0.09, 0, 0} // TOF trigger, Eur. Phys. J. C (2013) 73:2617, now in central trigger eff
  // };
  static const Int_t nCorA = 2; // sources of correlated error
  const Double_t ALICE_photo_corr_A[nCorA][ALICE_photo_n] = // relative uncertainty!
  {
    // {0,0,0,0}
    {0.00001/ALICE_photo_sig[0],0.00002/ALICE_photo_sig[1],0.00003/ALICE_photo_sig[2],-0.00039/ALICE_photo_sig[3],-0.00063/ALICE_photo_sig[4],-0.00125/ALICE_photo_sig[5], 0, 0,  0},
    {0.00007/ALICE_photo_sig[0],0.00016/ALICE_photo_sig[1],0.00032/ALICE_photo_sig[2],-0.00577/ALICE_photo_sig[3],-0.00865/ALICE_photo_sig[4],-0.01170/ALICE_photo_sig[5], 0, 0,  0}
    // {0,0,0,0,0}
  };
  Double_t ALICE_photo_tot[ALICE_photo_n]; // to store total error
  Double_t ALICE_photo_statsys[ALICE_photo_n]; // to store stat + total sys error for the plot
  Double_t Nval, Dval; // power-law parameters for further use
  TGraph *fitGraph; // to store fit result
  void TotalUncA();
  void TotalErrorA();
  void StatSysError();
  TGraph* GetShade(Double_t N, Double_t d, Double_t cov00, Double_t cov11, Double_t cov01);
  void makeFit();
  static alice alic;
  alice(Int_t oo=0);
  alice(const alice&);
  Double_t Nexp   = 0.0000001;
  Double_t rhoExp = 0.0000001;
public:
  static alice* instance() {return &alic;}
  void fcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  TGraph *fit() const {return fitGraph;}
  Double_t getPowerLaw(Double_t wval) const {return Nval*pow(wval/90.,Dval);}
  void getPowerLaw(Double_t &nn, Double_t &dd) const {nn = Nval; dd = Dval;}
  TGraphErrors *getData();
  TGraphErrors *getStarlight();
  TGraphErrors *getStarlightWithNuclear();
  void printData();
};//alice
alice alice::alic;

alice::alice(Int_t oo) : opt(oo), Nval(-1), Dval(-1), fitGraph(0) {

  for(Int_t i=0; i<ALICE_photo_n; i++) {ALICE_photo_unc[i]=0.; ALICE_photo_tot[i]=0.; ALICE_photo_statsys[i]=0.;}

  if (opt == 0) {
    TotalUncA();
    // TotalErrorA();
    StatSysError();
    makeFit();
  } else {
    cout << " Option not recognized. Bye!" << endl;
    exit(1);
  }

}//alice

void alice::TotalUncA()
{
// compute total uncorrelated uncertainty
  for(Int_t j=0;j<ALICE_photo_n;j++) {
    Double_t unc = 0;
    for(Int_t i=0;i<nUncA;i++) unc+=TMath::Power(ALICE_photo_unc_A[i][j],2);
    ALICE_photo_unc[j]=TMath::Sqrt(unc);
    //  cout << " Total uncorrelated uncertainty for measurement " << j << " is " << ALICE_photo_unc[j] << endl;
  }
}
void alice::TotalErrorA()
{
  // to print also sys err on differential cross section
  // Double_t dSigmaDy[] = {6.9, 8.7, 10.6, 10.0, 7.1};
  // compute total  uncertainty
  for(Int_t j=0;j<ALICE_photo_n;j++) {
    Double_t tot = (ALICE_photo_unc[j]*ALICE_photo_unc[j]);
    for(Int_t i=0;i<nCorA;i++) tot+=TMath::Power(ALICE_photo_corr_A[i][j],2);
    ALICE_photo_tot[j]=TMath::Sqrt(tot)*ALICE_photo_sig[j];
    cout << " Total uncertainty for measurement " << j << " is ";
    cout << fixed << setprecision(1) << ALICE_photo_tot[j] << "  (" << 100*TMath::Sqrt(tot) << "%)";
    // if(j > 2 && j < 8) {
    //   cout << ",   Delta dSigma/dy = " << TMath::Sqrt(tot)*dSigmaDy[j-3];
    // }
    cout << endl;
  }
}
void alice::StatSysError()
{
// quadratic sum of statistical uncertainty and total systematic uncertainty, to be used for the plot
  for(Int_t i=0; i<ALICE_photo_n; i++) {
    Double_t err = 0.;
    // err = TMath::Power(ALICE_photo_sta[i], 2) + TMath::Power(ALICE_photo_tot[i], 2);
    err = TMath::Power(ALICE_photo_sta[i], 2) + TMath::Power(ALICE_photo_sys[i], 2);
    ALICE_photo_statsys[i] = TMath::Sqrt(err);
  }
}
TGraph* alice::GetShade(Double_t N, Double_t d, Double_t cov00, Double_t cov11, Double_t cov01)
{
// Define graph with error band from results of the fit
  // Float_t wmin=gxmin;
  // Float_t wmax=gxmax;
  Float_t wmin=18.;
  Float_t wmax=1018;
  Int_t nw = 200;
  Color_t color = 29;
  Double_t dw = (wmax-wmin)/nw;
  TGraph *gshade = new TGraph(2*nw+2);
  for (Int_t i=0;i<=nw;i++){
    Double_t w = wmin+i*dw;
    Double_t mJpsiParticle = TDatabasePDG::Instance()->GetParticle(443)->Mass();
    Double_t mProton       = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    Double_t mLead         = 207.2 * 0.9315;
    // Double_t y = N*pow(w/90,d);
    // Double_t y = N*pow(w/90,d)*(1. - (mLead+mJpsiParticle)*(mLead+mJpsiParticle)/(w*w)  );
    // Double_t y = N*pow(w/90,d)*(1. - (10.*mProton+mJpsiParticle)*(10.*mProton+mJpsiParticle)/(w*w)  );
    Double_t thresholdfactor = (1. - (1.*mProton+mJpsiParticle)*(1.*mProton+mJpsiParticle)/(w*w)  );
    // Double_t y = N*pow(w/90,d)*(1. - (1.*mProton+mJpsiParticle)*(1.*mProton+mJpsiParticle)/(w*w)  );
    Double_t y = N*pow(w/90,d)*thresholdfactor*thresholdfactor;
    Double_t dydN = y/N;
    Double_t dydd = y*log(w/90);
    Double_t dy = sqrt(dydN*dydN*cov00 + dydd*dydd*cov11 + 2*dydN*dydd*cov01);
    gshade->SetPoint(i       ,w,y+dy);
    gshade->SetPoint(2*nw-i+1,w,y-dy);
  }
  gshade->SetFillColor(color);
  gshade->SetLineColor(color);
  return gshade;
}
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

  Double_t mJpsiParticle = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  Double_t mProton       = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mLead         = 207.2 * 0.9315;


  //Double_t bj[nCorA] = {par[2],par[3],par[4],par[5]};
  Double_t bj[nCorA];
  for(Int_t i=0; i<nCorA; i++) {
    bj[i] = par[i+2];
  }

  Double_t Sbj = 0;
  for(Int_t j=0;j<nCorA;j++) Sbj += (bj[j]*bj[j]);
  Double_t chi2 = 0;
  for (Int_t i=0; i<ALICE_photo_n;i++) {
    Double_t mu_i = ALICE_photo_sig[i];
    // Double_t m_i = par[0]*TMath::Power(ALICE_photo_W[i]/90,par[1]);
    // Double_t m_i = par[0]*(1. - (10.*mProton+mJpsiParticle)*(10.*mProton+mJpsiParticle)/(ALICE_photo_W[i]*ALICE_photo_W[i])  )*TMath::Power(ALICE_photo_W[i]/90,par[1]);
    Double_t thresholdfactor = (1. - (1.*mProton+mJpsiParticle)*(1.*mProton+mJpsiParticle)/(ALICE_photo_W[i]*ALICE_photo_W[i])  );
    // Double_t m_i = par[0]*(1. - (1.*mProton+mJpsiParticle)*(1.*mProton+mJpsiParticle)/(ALICE_photo_W[i]*ALICE_photo_W[i])  )*TMath::Power(ALICE_photo_W[i]/90,par[1]);
    Double_t m_i = par[0]*thresholdfactor*thresholdfactor*TMath::Power(ALICE_photo_W[i]/90,par[1]);
    // Double_t m_i = par[0]*(1. - (mLead+mJpsiParticle)*(mLead+mJpsiParticle)/(ALICE_photo_W[i]*ALICE_photo_W[i])  )*TMath::Power(ALICE_photo_W[i]/90,par[1]);
    Double_t Sij = 0;
    for(Int_t j=0;j<nCorA;j++) Sij += m_i*bj[j]*ALICE_photo_corr_A[j][i];
    Double_t d_stat_i =  ALICE_photo_sta[i]/ALICE_photo_sig[i];
    Double_t d_unc_i = ALICE_photo_unc[i];
    Double_t D = (d_stat_i*d_stat_i*mu_i*(m_i-Sij))+((d_unc_i*m_i)*(d_unc_i*m_i));
    chi2 += ((m_i-mu_i-Sij)*(m_i-mu_i-Sij)/D);
    ALICE_chi2[i] = ((m_i-mu_i-Sij)*(m_i-mu_i-Sij)/D);
  }
  chi2 += Sbj;
  f = chi2;
}
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
  cout << " N     = " << fixed << setprecision(4) << N << " +/- " << Nerr << endl;
  cout << " delta = " << D << " +/- " << Derr << endl;
  cout << " rho   = " << fixed << setprecision(2) << rho << endl;
  cout << " chi2  = " << fixed << setprecision(2) << fmin << endl;
  cout << " stat: " << istat << endl;

  delete myMinuit;

  Nval = N;
  Dval = D;

  fitGraph = GetShade(N,D,covNN,covDD,covND);
  Nexp = N;
  rhoExp = D;
  cout << "=========================" << endl;
  cout << "Single Chi2 contributions" << endl;
  cout << "=========================" << endl;
  cout << "First  = " << ALICE_chi2[0] << endl;
  cout << "Second = " << ALICE_chi2[1] << endl;
  cout << "Third  = " << ALICE_chi2[2] << endl;
  cout << "Fourth = " << ALICE_chi2[3] << endl;
  cout << "Fifth  = " << ALICE_chi2[4] << endl;
  cout << "Sixth  = " << ALICE_chi2[5] << endl;


}

void ifcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //interface to chi2 function

  alice::instance()->fcnChi2ModelA(npar, gin, f, par, iflag);

}//interface to chi2 function

TGraphErrors *alice::getData() {
  //forward data
  const Int_t   np       = 11;
  // const Double_t w[np]   = {ALICE_photo_W[0], ALICE_photo_W[1], ALICE_photo_W[2], ALICE_photo_W[8]};
  // const Double_t sig[np] = {ALICE_photo_sig[0], ALICE_photo_sig[1], ALICE_photo_sig[2], ALICE_photo_sig[8]};
  // const Double_t err[np] = {ALICE_photo_statsys[0], ALICE_photo_statsys[1], ALICE_photo_statsys[2], ALICE_photo_statsys[8]};

  TGraphErrors *gr = new TGraphErrors(np, ALICE_photo_W, ALICE_photo_sig, NULL, ALICE_photo_statsys);
  SetGraphStyle(gr ,kRed  ,kFullDiamond     ,2.,kRed   , 1,2);

  return gr;

}//getForward



TGraphErrors *alice::getStarlight() {
  //forward data
  const Int_t   np       = 6;
  const Double_t w[np]   = {ALICE_photo_W[0], ALICE_photo_W[1], ALICE_photo_W[2], ALICE_photo_W[3], ALICE_photo_W[4], ALICE_photo_W[5]};
  const Double_t sig[np] = {1.3259E-02, 1.7995E-02, 2.2492E-02, 1.4200E-01, 1.6707E-01, 1.9655E-01};
  // const Double_t err[np] = {ALICE_photo_statsys[0], ALICE_photo_statsys[1], ALICE_photo_statsys[2], ALICE_photo_statsys[8]};
  // 1.9176E+01   2.0459E+02   1.3259E-02
  // 2.4610E+01   1.8907E+02   1.7995E-02
  // 3.1591E+01   1.7355E+02   2.2492E-02
  // 4.9396E+02   1.2312E+01   1.4200E-01
  // 6.3425E+02   5.0905E+00   1.6707E-01
  // 8.1439E+02   1.4145E+00   1.9655E-01



  TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, NULL);
  SetGraphStyle(gr ,kMagenta  ,kFullCircle     ,2.,kMagenta   , 1,2);

  return gr;

}//getForward

TGraphErrors *alice::getStarlightWithNuclear() {
  //forward data
  const Int_t   np       = 9;
  // 9.7270E+01   1.0364E+02   3.2432E-02
  // 1.2489E+02   8.8152E+01   3.6961E-02
  // 1.6037E+02   7.2753E+01   4.1996E-02
  const Double_t w[np]   = {ALICE_photo_W[0], ALICE_photo_W[1], ALICE_photo_W[2], ALICE_photo_W[3], ALICE_photo_W[4], ALICE_photo_W[5],       ALICE_photo_W[6], ALICE_photo_W[7], ALICE_photo_W[8] };
  const Double_t sig[np] = {1.0410E-02, 1.3791E-02, 1.6819E-02, 7.1943E-02, 8.0387E-02, 8.9516E-02,         3.2432E-02, 3.6961E-02, 4.1996E-02};
  const Double_t err[np] = {ALICE_photo_statsys[0], ALICE_photo_statsys[1], ALICE_photo_statsys[2], ALICE_photo_statsys[8]};
  // // IA
  // 1.9176E+01   2.0459E+02   1.3259E-02
  // 2.4610E+01   1.8907E+02   1.7995E-02
  // 3.1591E+01   1.7355E+02   2.2492E-02
  // 4.9396E+02   1.2312E+01   1.4200E-01
  // 6.3425E+02   5.0905E+00   1.6707E-01
  // 8.1439E+02   1.4145E+00   1.9655E-01

  // // With Nuclear
  // 1.9176E+01   2.0459E+02   1.0410E-02
  // 2.4610E+01   1.8907E+02   1.3791E-02
  // 3.1591E+01   1.7355E+02   1.6819E-02
  // 4.9396E+02   1.2312E+01   7.1943E-02
  // 6.3425E+02   5.0905E+00   8.0387E-02
  // 8.1439E+02   1.4145E+00   8.9516E-02

  // 9.3456E+01   1.0613E+02   3.1752E-02
  // 1.1022E+02   9.5890E+01   3.4636E-02
  // 1.2489E+02   8.8152E+01   3.6961E-02
  // 1.2489E+02   8.8152E+01   3.6961E-02
  // 1.4152E+02   8.0435E+01   3.9413E-02
  // 1.6608E+02   7.0611E+01   4.2744E-02


  TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, NULL);
  SetGraphStyle(gr ,kBlue  ,kFullCircle     ,2.,kBlue   , 1,2);

  return gr;

}//getForward









//_____________________________________________________________________________
TGraph *getRatio(TGraph *g) {

  //ratio of model from TGraph argument over the data

  Double_t wval, sigma;
  Int_t i = 0;

  //output graph
  TGraph *gout = (TGraph*) g->Clone();

  //alice data
  alice *alic = alice::instance();

  //input graph loop
  while( g->GetPoint(i, wval, sigma) >=0 ) {

    gout->SetPoint(i, wval, sigma/alic->getPowerLaw(wval) );

    i++;
  }//input graph loop

  return gout;

}//getRatio
