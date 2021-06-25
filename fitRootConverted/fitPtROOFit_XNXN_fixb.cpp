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
TH1F* hun;
TF1* fgu;



void fitPtROOFit( const char* AnalysisName, Int_t selectionFlag = 0 )
{
  // ---------------------------------------------------
  // I m p o r t i n g   R O O T   h i s t o g r a m s
  // ===================================================
  // I m p o r t   T H 1   i n t o   a   R o o D a t a H i s t
  // ---------------------------------------------------------
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("ADcheck_", listings);

  if      ( selectionFlag == 0 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeH");
  else if ( selectionFlag == 1 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeZNCzeroZNAzeroH");
  else if ( selectionFlag == 2 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeZNCzeroZNAanyH");
  else if ( selectionFlag == 3 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeZNCanyZNAzeroH");
  else if ( selectionFlag == 4 ) hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeZNCanyZNAanyH");
  else                           hPt = (TH1F*)listings->FindObject("fDimuonPtDistributionDissociativeH");


  hPt->Rebin(10);








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




  // ---Create Incoherent Disocitation PDF
  // Values taken from https://arxiv.org/abs/1304.5162
  RooRealVar b("b","b",1.79, 0., 2.5); //1.67, 1.91);
  RooRealVar negativeb("negativeb","negativeb",-1.79, -2.5, -0.5); //1.67, 1.91);
  RooRealVar n("n","n",3.58, 2.5,4.5); //3.43,3.73);

  // b.setConstant(kTRUE);
  n.setConstant(kTRUE);

  RooGenericPdf  PDF_IncohJpsiToX("PDF_IncohJpsiToX","pT*pow((1+pow(pT,2)*b/n),-n)",RooArgSet(pT, b, n));







  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist rPt("rPt", "rPt", pT, Import(*hPt));
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






  RooHistPdf pdfun("pdfun", "pdfun", pT, run, 0);












  RooRealVar dissoc  ( "dissoc", "number of dissoc",    1000, 0, 10000 );
  RooAddPdf  sum     ( "sum",    "extended sum of all",
                       // RooArgList(pdfun),
                       RooArgList(PDF_IncohJpsiToX),
                       RooArgList(dissoc)
                       );




  // RooFitResult* r = sum.fitTo(rPt,Extended(kTRUE),Save());
  RooFitResult* r = sum.fitTo(rPt,Range(0.2, 3.),Extended(kTRUE),Save());




  // sum.plotOn (frame,LineColor(kBlack), Range(0.,3.0) ) ;
  // sum.paramOn(frame);
  // sum.paramOn(frame,rPt);
  // sum.plotOn (frame,Components(pdfun), LineColor(kMagenta), Range(0.,3.0)  );
  sum.plotOn (frame,Components(PDF_IncohJpsiToX), LineColor(kMagenta), Range(0.,3.0)  );
  // sum.paramOn(frame,rPt);
  frame->Draw();















//
//   // get fractions of PDF in jRecPt range [0,Pt cut]
//   pT.setRange("signal",0.0,0.25);
//   RooAbsReal* cohI    = pdf1c.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* icI     = pdf1i.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* p2scohI = pdffc.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* p2sicI  = pdffi.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* ggI     = pdfgl.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* funI    = pdfun.createIntegral(pT,NormSet(pT),Range("signal")) ;
//   RooAbsReal* sI      = sum.createIntegral  (pT,NormSet(pT),Range("signal")) ;
//
//
//
//
//
//   //_______________________________
//   // TOMAS SNIPPET
//   RooRealVar ncoh(    "ncoh",    "ncoh",    cohI->getVal(),    "" );
//   RooRealVar ncoh2(   "ncoh2",   "ncoh2",   p2scohI->getVal(), "" );
//   RooRealVar nincoh(  "nincoh",  "nincoh",  icI->getVal(),     "" );
//   RooRealVar nincoh2( "nincoh2", "nincoh2", p2sicI->getVal(),  "" );
//   RooRealVar ndisoc(  "ndisoc",  "ndisoc",  funI->getVal(),    "" );
//
//   RooFormulaVar var_ncoh(     "var_ncoh",    "cohN*ncoh",      RooArgList(cohN,ncoh)      );
//   RooFormulaVar var_ncoh2(    "var_ncoh2",   "p2scohN*ncoh2",  RooArgList(p2scohN,ncoh2)  );
//   RooFormulaVar var_nincoh(   "var_nincoh",  "icN*nincoh",     RooArgList(icN,nincoh)     );
//   RooFormulaVar var_nincoh2(  "var_nincoh2", "p2sicN*nincoh2", RooArgList(p2sicN,nincoh2) );
//   RooFormulaVar var_ndisoc(   "var_ndisoc",  "dissoc*ndisoc",  RooArgList(dissoc,ndisoc)  );
//   RooFormulaVar f_I_extended( "f_I_extended",
//                               "(var_nincoh+var_nincoh2+var_ndisoc)/(var_ncoh+var_ncoh2)",
//                               RooArgList(var_nincoh, var_nincoh2, var_ndisoc, var_ncoh, var_ncoh2)
//                               );
//
//   //_______________________________
//
//
//
//
//
//   //Calculate chi^2
//   Double_t total = sum.expectedEvents(pT);
//   cout << "total = " << total << endl;
//   Double_t chi2 = 0;
//   for( Int_t i = 0; i < 120; i++){
//   // for( Int_t i = 0; i < nBinsX; i++){
//     // pT.setRange( "bin", i*(1.0/((Double_t)nBinsX)), (i+1)*(1.0/((Double_t)nBinsX)) );
//     // RooAbsReal* binI   = sum.createIntegral( pT, NormSet(pT), Range("bin") );
//     pT.setRange( Form("bin%i", i), i*0.025, (i+1)*0.025 );
//     RooAbsReal* binI   = sum.createIntegral( pT, NormSet(pT), Range(Form("bin%i", i)) );
//    	Double_t fBin      = binI->getVal();
//    	Double_t sumPoint  = fBin * total;
//    	cout << "sumPoint" << sumPoint << " ";
//    	Double_t dataPoint = hPt->GetBinContent(i+1);
//    	cout <<  dataPoint << endl;
//    	Double_t sqDiff    = (sumPoint - dataPoint)*(sumPoint-dataPoint);
//    	cout <<  sqDiff    << endl;
//    	Double_t sqError   = sumPoint;
//    	Double_t chiBin    = sqDiff/sqError;
//    	// cout <<  chiBin    << endl;
//    	chi2+=chiBin;
//   }
//
//   cout << "chi^2 =     " << chi2 << endl;
//   // cout << "chi^2/dof = " << chi2/((Double_t)nBinsX) << endl;
//   cout << "chi^2/dof = " << chi2/((Double_t)120. ) << endl;
//
//   Double_t fCoh     = cohI   ->getVal();
//   Double_t fICoh    = icI    ->getVal();
//   // Double_t fPsiCoh  = p2scohI->getVal();
//   Double_t fPsiICoh = p2sicI ->getVal();
//   Double_t fGG      = ggI    ->getVal();
//   Double_t fdisso   = funI   ->getVal();
//
//
//  cout << "For jRectPt<Pt cut GeV/c" << endl;
//  cout << "Number of gg integrated under Pt cut: "              << ggN.getValV()    * fGG      << " +/- " << ggN.getError()    * fGG      << endl;
//  // cout << "Number of coh feed-down integrated under Pt cut: " << p2scohN.getValV() * fPsiCoh << " +/- " << p2scohN.getError() * fPsiCoh << endl;
//  // cout << "Number of incoh feed-down integrated under Pt cut: " << p2sicN.getValV() * fPsiICoh << " +/- " << p2sicN.getError() * fPsiICoh << endl;
//  cout << "Number of incoherent integrated under Pt cut: "      << icN.getValV()    * fICoh    << " +/- " << icN.getError()    * fICoh    << endl;
//  cout << "Number of coherent integrated under Pt cut: "        << cohN.getValV()   * fCoh     << " +/- " << cohN.getError()   * fCoh     << endl;
//  // Double_t N_coh2      = cohN.getVal()    * fCoh; //Number of coherent J/psi according to fit
//  Double_t N_coh2      = cohN.getVal()    * fCoh; //Number of coherent J/psi according to fit
//  // Double_t N_coh2      = cohN.getValV()    * fCoh; //Number of coherent J/psi according to fit
//  Double_t N_cohError  = cohN.getError()   * fCoh; //Error on coh J/psi according to fit
//  // Double_t N_I         =      fICoh; //Number of incoherent J/psi below Pt cut
//  Double_t N_I         = icN.getVal()     * fICoh; //Number of incoherent J/psi below Pt cut
//  // Double_t N_I         = icN.getValV()     * fICoh; //Number of incoherent J/psi below Pt cut
//  Double_t N_IError    = icN.getError()    * fICoh; //Error on inc J/psi below Pt cut
//  // Double_t N_cohFD = p2scohN.getValV() * fPsiCoh; //Number of coh FD below Pt cut
//  // Double_t N_cohFDError = p2scohN.getError() * fPsiCoh; //Error on coh FD < Pt cut
//  // Double_t N_icFD      = p2sicN.getValV()  * fPsiICoh; //Number of ic FD below Pt cut
//  // Double_t N_icFDError = p2sicN.getError() * fPsiICoh; //Error on ic FD below Pt cut
//  // Double_t N_diss      = fdisso; //Number of ic FD below Pt cut
//  Double_t N_diss      = dissoc.getVal()  * fdisso; //Number of ic FD below Pt cut
//  // Double_t N_diss      = dissoc.getValV()  * fdisso; //Number of ic FD below Pt cut
//  Double_t N_dissError = dissoc.getError() * fdisso; //Error on ic FD below Pt cut
//
//
//  cout << "coh = " << cohN.getVal() << endl;
//  cout << "ich = " << icN.getVal() << endl;
//  cout << "dissoc = " << dissoc.getVal() << endl;
//  //When cross sections and efficiencies are given for pt<Pt cut:
//  //cout << "fFDcoh = " << p2scohN.getValV() * fPsiCoh/(cohN.getValV() * fCoh) << endl;
//  //cout << "fFDincoh = " << p2sicN.getValV() * fPsiICoh/(cohN.getValV() * fCoh) << endl;
//  //When cross sections and effeciences are given for all pt:
//  // fracFDcoh = p2scohN.getValV()/cohN.getValV();
//  //  fracFDic = p2sicN.getValV()/icN.getValV();
//  //
//  //  cout << "fFDcoh = " << fracFDcoh << ";  Target: 0.105" << endl;
//  //  cout << "fFDincoh = " << fracFDic << ";  Target: 0.0993" << endl;
//
//  cout << "****************************** " << endl;
//  cout << "Fractions of the templates under Pt=Pt cut GeV/c" << endl;
//  cout << "Coherent fraction " << fCoh << endl;
//  cout << "Incoherent fraction " << fICoh << endl;
//  // cout << "Coherent Feed-down fraction " << fPsiCoh << endl;
//  cout << "Incoherent Feed-down fraction " << fPsiICoh << endl;
//  cout << "Gamma+Gamma fraction " << fGG << endl;
//
//  Double_t f_I   =  (N_I     + N_diss     ) / (N_coh2);
//  Double_t ErrfI = sqrt((N_IError + N_dissError)*(N_IError + N_dissError) / ((N_I + N_diss)*(N_I + N_diss)) + (N_cohError*N_cohError)/(N_coh2*N_coh2) )*f_I;
//  // fFD = ((p2scohN.getValV() * fPsiCoh) + (p2sicN.getValV() * fPsiICoh)) / (cohN.getValV() * fCoh);
//
//
//
//
//
//
//
//
//  Double_t lumi = 11.92;
//
//  TLatex* latex = new TLatex();
//  latex->SetNDC();
//  latex->SetTextSize(0.045);
//  latex->SetTextAlign(22);
//  latex->DrawLatex(0.56,0.95,Form("ALICE, Pb#minusp #sqrt{#it{s}_{NN}} = 8.16 TeV"));
//  // latex->DrawLatex(0.45,0.86,Form("f_{I} = #frac{%.3f + %.3f}{%.3f} = %.3f #pm %.3f", N_I, N_diss, N_coh2, f_I, ErrfI));
//  latex->DrawLatex(0.45,0.86,Form("f_{I} =  %.3f #pm %.3f", f_I_extended.getVal(), f_I_extended.getPropagatedError(*r)));
//  if ( selectionFlag == 7 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -2.500));
//  if ( selectionFlag == 8 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.250));
//  if ( selectionFlag == 9 ) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.250,  -2.500));
//  if ( selectionFlag == 10) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.500));
//  if ( selectionFlag == 11) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.500,  -3.000));
//  if ( selectionFlag == 12) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.000,  -2.500));
//  if ( selectionFlag == 13) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.625));
//  if ( selectionFlag == 14) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.625,  -3.250));
//  if ( selectionFlag == 15) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.250,  -2.875));
//  if ( selectionFlag == 16) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-2.875,  -2.500));
//  if ( selectionFlag == 17) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-4.000,  -3.700));
//  if ( selectionFlag == 18) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.700,  -3.400));
//  if ( selectionFlag == 19) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.400,  -3.100));
//  if ( selectionFlag == 20) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-3.100,  -2.800));
//  if ( selectionFlag == 21) latex->DrawLatex(0.80,0.20,Form("%.3f < y < %.3f",-2.800,  -2.500));
//  latex->SetTextSize(0.033);
//  latex->SetTextAlign(12);
//  // latex->DrawLatex(0.48,0.89,Form("UPC, L_{#lower[-0.3]{int}} = %.0f #pm %.0f #mub^{-1}",lumi,lumi*0.05));
//  // latex->DrawLatex(0.48,0.85,Form("#minus%.2f < #it{y} < #minus%.2f",-gYMin[iy],-gYMax[iy]));
//  // latex->DrawLatex(0.48,0.81,Form("%.2f < #it{m}_{#mu#mu} < %.2f GeV/#it{c}^{2}",gMMin[im],gMMax[im]));
//
//  if (0) {
// //    TLegend* l = new TLegend(0.45,0.47,0.98,0.77);
// // //    TLegend* l = new TLegend(0.55,0.50,0.98,0.84);
// //    l->SetMargin(0.1);
// //    l->AddEntry(hPt,"ALICE data","p");
// //    l->AddEntry(fsum,Form("Fit: #chi^{2}/NDF=%.1f\n",chi2ndf),"l");
// //    l->AddEntry(h1c,Form("Coherent J/#psi: %.0f #pm %.0f",n1c,n1c_err));
// //    l->AddEntry(h1i,Form("Incoherent J/#psi: %.0f #pm %.0f",n1i,n1i_err));
// //    l->AddEntry(hun,Form("Incoherent dissocitive J/#psi: %.0f",nun),"l");
// //    // l->AddEntry(hfc,"Coherent #psi' feeddown");
// //    l->AddEntry(hfi,"Incoherent #psi' feeddown");
// //    l->AddEntry(hgl,Form("Continuum #gamma#gamma #rightarrow #mu#mu: %.0f",ngl));
// //    l->Draw();
//  }  else {
//    // TLegend* l = new TLegend(0.47,0.48,0.985,0.78);
//    // l->SetMargin(0.09);
//    // l->AddEntry(hPt,"ALICE data","p");
//    // l->AddEntry(pdf1c,"Coherent J/#psi");
//    // l->AddEntry(pdf1i,"Incoherent J/#psi");
//    // l->AddEntry(pdfun,"Incoherent J/#psi with nucleon dissociation");
//    // // l->AddEntry(hfc,"Coherent J/#psi from #psi' decay");
//    // l->AddEntry(pdffi,"Incoherent J/#psi from #psi' decay");
//    // l->AddEntry(pdfgl,"Continuum #gamma#gamma #rightarrow #mu#mu");
//    // l->AddEntry(sum,"Sum");
//    // l->Draw();
//  }
//
//
//
//
//
//  for (int i=0; i<frame->numItems(); i++) {
//    TString obj_name=frame->nameOf(i); if (obj_name=="") continue;
//    cout << Form("%d. '%s'\n",i,obj_name.Data());
//  }
//
//
//
//  TString names[] = {
//    "h_rPt",
//    "sum_Norm[pT]_Range[0.000000_1.000000]",
//    "sum_Norm[pT]_Comp[pdf1c]_Range[0.000000_1.600000]",
//    "sum_Norm[pT]_Comp[pdf1i]_Range[0.000000_1.600000]",
//    "sum_Norm[pT]_Comp[pdffi]_Range[0.000000_1.600000]",
//    "sum_Norm[pT]_Comp[pdfgl]_Range[0.000000_1.600000]",
//    "sum_Norm[pT]_Comp[pdfun]_Range[0.000000_1.600000]",
//    ""
//  };
//
//  TString signs[] = {
//    "ALICE data",
//    "Total fit",
//    "#gamma+Pb",
//    "Exclusive",
//    "#gamma#gamma",
//    "Non-exclusive"
//  };
//
//
//  TLegend* l = new TLegend(0.47,0.48,0.985,0.78);
//  l->SetMargin(0.09);
//  Int_t i=-1;
//  while ( names[++i] != "" ) {
//    TObject *obj = frame->findObject(names[i].Data());
//    if (!obj) {
//      Warning("fitBi4",Form("Can't find item = %s in the frame2!\n",names[i].Data()));
//      continue;
//    }
//    l->AddEntry(obj,signs[i],"l");
//  }
//  l->Draw();
//
//
//
//
//  // ------------------------------------------------------------------------------------------------------------------------------------
//  // Pt plot
//  //---Create pt canvas
//  TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
//  // cPt->Divide(2);
//
//  // //-----------------------------------------------------------------------------------
//  // // Draw Correlation Matrix
//  // cPt->cd(2);
//  // TPad *cPt2 = new TPad("cPt2", "cPt2",0.001,0.001,0.999,0.999);
//  // // ---Pad option
//  // cPt2->SetRightMargin(0.15);
//  // cPt2->SetLeftMargin(0.15);
//  // cPt2->Draw();
//  // cPt2->cd();
//  // // ---TH2D correlation matrix
//  // TH2* hCorr_pt = fit_pt->correlationHist();
//  // hCorr_pt->SetMarkerSize(1.6);
//  // hCorr_pt->GetXaxis()->SetLabelSize(0.045);
//  // hCorr_pt->GetYaxis()->SetLabelSize(0.045);
//  // hCorr_pt->Draw("zcol,text");
//  //
//  //-----------------------------------------------------------------------------------
//  // Draw pt Histogram
//  // cPt->cd(1);
//  TPad *cPt1 = new TPad("cPt1", "cPt1",0.001,0.001,0.999,0.999);
//  // ---Pad option
//  cPt1->SetLogy();
//  cPt1->SetLeftMargin(0.14);
//  cPt1->SetRightMargin(0.01);
//  cPt1->SetBottomMargin(0.12);
//  cPt1->Draw();
//  cPt1->cd();
//  // ---TH1 pt spectrum
//  // RooPlot* frame_pt = pT.frame(Title("Pt fit")) ;
//  RooPlot* frame_pt = pT.frame(Title(" ")) ;
//  rPt.plotOn(frame_pt,Name("rPt")/*,Binning(bin_pt),*/,MarkerStyle(20),MarkerSize(0.9));
//  sum.plotOn(frame_pt,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
//  // ---Psi(2S)
//  // Pt_fit_func.plotOn(frame_pt,Name("PDF_CohPsi2sToMu"), Components(PDF_CohPsi2sToMu), LineColor(kBlue+1), LineWidth(2));
//  // Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohPsi2sToMu"),Components(PDF_IncohPsi2sToMu), LineColor(kRed+1), LineWidth(2));
//  // ---J/Psi
//  sum.plotOn(frame_pt,Name("PDF_CohJpsiToMu"),       Components(pdf1c), LineColor(kBlue+1), LineWidth(2));
//  sum.plotOn(frame_pt,Name("PDF_IncohJpsiToMu"),     Components(pdf1i), LineColor(kRed+1), LineWidth(2));
//  sum.plotOn(frame_pt,Name("PDF_CohPsi2sToMuPi"),    Components(pdffc), LineColor(kCyan+1), LineWidth(2));
//  sum.plotOn(frame_pt,Name("PDF_IncohPsi2sToMuPi"),  Components(pdffi), LineColor(kOrange), LineWidth(2));
//  sum.plotOn(frame_pt,Name("PDF_IncohJpsiToX"),      Components(pdfun), LineColor(kMagenta+1), LineWidth(2));
//  sum.plotOn(frame_pt,Name("PDF_TwoGammaToMuMedium"),Components(pdfgl), LineColor(kGreen+2), LineWidth(2));
//  // ---Frame pt option
//  frame_pt->SetAxisRange(0,3,"X");
//  frame_pt->SetAxisRange(0.025,11000000,"Y");
//  frame_pt->GetXaxis()->SetTitle("Dimuon #it{p}_{T} (GeV/#it{c})");
//  frame_pt->GetYaxis()->SetTitle("Counts per  0.25 GeV/#it{c} ");
//  // frame_pt->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (pt_range_max-pt_range_min)/n_bins_pt*1000));
//  frame_pt->GetYaxis()->SetTitleOffset(1.3);
//  frame_pt->Draw();
//  //---Calculate chi2 pt
//  // Double_t chi2_pt = frame_pt->chiSquare("sum","rPt",r->floatParsFinal().getSize());
//  // Double_t chi2_pt = frame_pt->chiSquare("sum","rPt",sum.getParameters(rPt)->selectByAttrib("Constant",kFALSE)->getSize());
//  Double_t chi2_pt = frame->chiSquare("sum","rPt",sum.getParameters(rPt)->selectByAttrib("Constant",kFALSE)->getSize());
//  // ---pt plot description
//  TLatex * text_pt1 = new TLatex (0.1,0.38*frame_pt->GetMaximum(),"Preliminary");
//  text_pt1->SetTextSize(0.035);
//  text_pt1->Draw();
//  TLatex * text_pt2 = new TLatex (0.85,0.38*frame_pt->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
//  text_pt2->SetTextSize(0.035);
//  text_pt2->Draw();
//  TLatex * text_pt3 = 0x0;
//  if      ( selectionFlag == 1 ) {
//     if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -2.5");
//     else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -3.5");
//     else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -3.5 < y < -3.0");
//     else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -3.0 < y < -2.5");
//     else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0N0N, -4.0 < y < -2.5");
//  }
//  else if ( selectionFlag == 2 ) {
//     if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -2.5");
//     else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -3.5");
//     else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -3.5 < y < -3.0");
//     else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -3.0 < y < -2.5");
//     else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"0NXN, -4.0 < y < -2.5");
//  }
//  else if ( selectionFlag == 3 ) {
//    if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -2.5");
//    else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -3.5");
//    else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -3.5 < y < -3.0");
//    else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -3.0 < y < -2.5");
//    else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XN0N, -4.0 < y < -2.5");
//  }
//  else if ( selectionFlag == 4 ) {
//    if      ( selectionFlag2 == 0 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
//    else if ( selectionFlag2 == 1 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -3.5");
//    else if ( selectionFlag2 == 2 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -3.5 < y < -3.0");
//    else if ( selectionFlag2 == 3 ) text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -3.0 < y < -2.5");
//    else                            text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
//  }
//  else                              text_pt3 = new TLatex (0.1,0.18*frame_pt->GetMaximum(),"XNXN, -4.0 < y < -2.5");
//  text_pt3->SetTextSize(0.035);
//  text_pt3->Draw();
//
//  // TLatex * text_pt3 = new TLatex (0.85,0.11*frame_pt->GetMaximum(),Form("#bf{%.1f < y < %.1f, %.2f < #it{m_{#mu^{+}#mu^{-}}} < %.2f}",y_min,y_max, m_min_jpsi, m_max_jpsi));
//  // text_pt3->Draw();
//  // TLatex * text_pt4 = new TLatex (0.1,0.11*frame_pt->GetMaximum(),Form("#bf{f_{I} = %.1f%% #pm %.1f%%}",100*f_I.getVal(), 100*f_I.getPropagatedError(*fit_pt) ));
//  // text_pt4->Draw();
//  // ---Legend pt
//  TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
//  leg_pt->SetFillStyle(0);
//  leg_pt->SetBorderSize(0);
//  leg_pt->SetTextSize(0.04);
//  leg_pt->AddEntry("rPt","ALICE Data", "P");
//  leg_pt->AddEntry("sum",Form("Fit: #chi^{2}/dof= %.3f",chi2/((Double_t)120. )),"L");
//  // ------Psi(2S)
//  // leg_pt->AddEntry("PDF_CohPsi2sToMu","Coherent #psi(2S)", "L");
//  // leg_pt->AddEntry("PDF_IncohPsi2sToMu","Incoherent #psi(2S)", "L");
//  // ------J/Psi
//  leg_pt->AddEntry("PDF_CohJpsiToMu","Coherent J/#psi", "L");
//  leg_pt->AddEntry("PDF_IncohJpsiToMu","Incoherent J/#psi", "L");
//  leg_pt->AddEntry("PDF_CohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
//  leg_pt->AddEntry("PDF_IncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
//  leg_pt->AddEntry("PDF_IncohJpsiToX","Nucleon disoc.", "L");
//  leg_pt->AddEntry("PDF_TwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
//  // leg_pt->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2_pt),"");
//  leg_pt->Draw();
//
//
//
//
//
//  gPad->SaveAs(Form("pngResults/fitPtROOFit_%d.png", selectionFlag),  "RECREATE");

}
