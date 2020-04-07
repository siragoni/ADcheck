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



Double_t fLumiPerRun = 0;
Int_t    fRunNum     = 0;
Double_t GlobalAxE   = 0;
// Double_t fLumi       = 10853.457;
// Double_t fLumi       = 11919.684;
Double_t fLumi       = 532.926230;

//_____________________________________________________________________________
/* - There were plenty of ways to do this...
 * - However, recently the STL libraries were
 * - creating confusion on the LEGO framework
 * - (they didn't fire at all).
 * - This problem was not found on local, where
 * - things were working properly...
 * - So I am using the most barbaric C-style
 * - arrays/for...
 */
void SetLuminosityCap()
{
  fLumiPerRun = 0;
  /* - Here I am rounding up the number for 10k,
   * - so that I don't have to do tedious conversions...
   * - I am adding 1 entry to the number obtained by 40k,
   * - so that I am not missing anything...
   * -
   */
  if      ( fRunNum == 295585 ) { fLumiPerRun = 0.0793352; }
  else if ( fRunNum == 295586 ) { fLumiPerRun = 0.238599; }
  else if ( fRunNum == 295587 ) { fLumiPerRun = 0.109518; }
  else if ( fRunNum == 295588 ) { fLumiPerRun = 0.135709; }
  else if ( fRunNum == 295589 ) { fLumiPerRun = 0.281897; }
  else if ( fRunNum == 295612 ) { fLumiPerRun = 0.448985; }
  else if ( fRunNum == 295615 ) { fLumiPerRun = 0.0565828; }
  else if ( fRunNum == 295665 ) { fLumiPerRun = 0.334733; }
  else if ( fRunNum == 295666 ) { fLumiPerRun = 0.323941; }
  else if ( fRunNum == 295667 ) { fLumiPerRun = 0.0970128; }
  else if ( fRunNum == 295668 ) { fLumiPerRun = 0.130088; }
  else if ( fRunNum == 295671 ) { fLumiPerRun = 0.325985; }
  else if ( fRunNum == 295673 ) { fLumiPerRun = 0.312556; }
  else if ( fRunNum == 295675 ) { fLumiPerRun = 0.13204; }
  else if ( fRunNum == 295676 ) { fLumiPerRun = 0.321284; }
  else if ( fRunNum == 295677 ) { fLumiPerRun = 0.265538; }
  else if ( fRunNum == 295714 ) { fLumiPerRun = 0.345554; }
  else if ( fRunNum == 295716 ) { fLumiPerRun = 0.33778; }
  else if ( fRunNum == 295717 ) { fLumiPerRun = 0.287941; }
  else if ( fRunNum == 295718 ) { fLumiPerRun = 0.257395; }
  else if ( fRunNum == 295719 ) { fLumiPerRun = 0.294713; }
  else if ( fRunNum == 295723 ) { fLumiPerRun = 0.506379; }
  else if ( fRunNum == 295725 ) { fLumiPerRun = 0.8453; }
  else if ( fRunNum == 295753 ) { fLumiPerRun = 0.384579; }
  else if ( fRunNum == 295754 ) { fLumiPerRun = 1.12405; }
  else if ( fRunNum == 295755 ) { fLumiPerRun = 1.01645; }
  else if ( fRunNum == 295758 ) { fLumiPerRun = 1.62041; }
  else if ( fRunNum == 295759 ) { fLumiPerRun = 0.532728; }
  else if ( fRunNum == 295762 ) { fLumiPerRun = 0.274725; }
  else if ( fRunNum == 295763 ) { fLumiPerRun = 1.0236; }
  else if ( fRunNum == 295786 ) { fLumiPerRun = 0.746701; }
  else if ( fRunNum == 295788 ) { fLumiPerRun = 3.06545; }
  else if ( fRunNum == 295791 ) { fLumiPerRun = 0.85713; }
  else if ( fRunNum == 295816 ) { fLumiPerRun = 1.20498; }
  else if ( fRunNum == 295818 ) { fLumiPerRun = 0.145167; }
  else if ( fRunNum == 295819 ) { fLumiPerRun = 2.74121; }
  else if ( fRunNum == 295822 ) { fLumiPerRun = 2.2132; }
  else if ( fRunNum == 295825 ) { fLumiPerRun = 0.255836; }
  else if ( fRunNum == 295826 ) { fLumiPerRun = 1.48001; }
  else if ( fRunNum == 295829 ) { fLumiPerRun = 0.934546; }
  else if ( fRunNum == 295831 ) { fLumiPerRun = 0.779334; }
  else if ( fRunNum == 295854 ) { fLumiPerRun = 1.73457; }
  else if ( fRunNum == 295855 ) { fLumiPerRun = 1.74464; }
  else if ( fRunNum == 295856 ) { fLumiPerRun = 1.73755; }
  else if ( fRunNum == 295859 ) { fLumiPerRun = 1.05037; }
  else if ( fRunNum == 295860 ) { fLumiPerRun = 0.834142; }
  else if ( fRunNum == 295861 ) { fLumiPerRun = 1.06703; }
  else if ( fRunNum == 295863 ) { fLumiPerRun = 0.728227; }
  else if ( fRunNum == 295881 ) { fLumiPerRun = 0.711494; }
  else if ( fRunNum == 295908 ) { fLumiPerRun = 2.92589; }
  else if ( fRunNum == 295909 ) { fLumiPerRun = 0.787373; }
  else if ( fRunNum == 295910 ) { fLumiPerRun = 3.19503; }
  else if ( fRunNum == 295913 ) { fLumiPerRun = 3.12937; }
  else if ( fRunNum == 295936 ) { fLumiPerRun = 1.47357; }
  else if ( fRunNum == 295937 ) { fLumiPerRun = 0.405657; }
  else if ( fRunNum == 295941 ) { fLumiPerRun = 1.67616; }
  else if ( fRunNum == 295942 ) { fLumiPerRun = 1.92337; }
  else if ( fRunNum == 295943 ) { fLumiPerRun = 1.67424; }
  else if ( fRunNum == 295945 ) { fLumiPerRun = 2.03661; }
  else if ( fRunNum == 295947 ) { fLumiPerRun = 2.69035; }
  else if ( fRunNum == 296061 ) { fLumiPerRun = 1.29676; }
  else if ( fRunNum == 296062 ) { fLumiPerRun = 1.80833; }
  else if ( fRunNum == 296063 ) { fLumiPerRun = 2.69974; }
  else if ( fRunNum == 296065 ) { fLumiPerRun = 2.45097; }
  else if ( fRunNum == 296066 ) { fLumiPerRun = 0.733886; }
  else if ( fRunNum == 296068 ) { fLumiPerRun = 1.98122; }
  else if ( fRunNum == 296123 ) { fLumiPerRun = 0.506486; }
  else if ( fRunNum == 296128 ) { fLumiPerRun = 0.445452; }
  else if ( fRunNum == 296132 ) { fLumiPerRun = 2.10822; }
  else if ( fRunNum == 296133 ) { fLumiPerRun = 1.73213; }
  else if ( fRunNum == 296134 ) { fLumiPerRun = 3.91041; }
  else if ( fRunNum == 296135 ) { fLumiPerRun = 2.3412; }
  else if ( fRunNum == 296142 ) { fLumiPerRun = 1.78935; }
  else if ( fRunNum == 296143 ) { fLumiPerRun = 0.534028; }
  else if ( fRunNum == 296191 ) { fLumiPerRun = 5.051; }
  else if ( fRunNum == 296192 ) { fLumiPerRun = 0.497364; }
  else if ( fRunNum == 296194 ) { fLumiPerRun = 2.84426; }
  else if ( fRunNum == 296195 ) { fLumiPerRun = 0.737647; }
  else if ( fRunNum == 296196 ) { fLumiPerRun = 2.42215; }
  else if ( fRunNum == 296197 ) { fLumiPerRun = 2.07279; }
  else if ( fRunNum == 296198 ) { fLumiPerRun = 0.813921; }
  else if ( fRunNum == 296241 ) { fLumiPerRun = 0.845868; }
  else if ( fRunNum == 296242 ) { fLumiPerRun = 0.95166; }
  else if ( fRunNum == 296243 ) { fLumiPerRun = 1.56742; }
  else if ( fRunNum == 296244 ) { fLumiPerRun = 8.37179; }
  else if ( fRunNum == 296246 ) { fLumiPerRun = 1.82994; }
  else if ( fRunNum == 296247 ) { fLumiPerRun = 1.1763; }
  else if ( fRunNum == 296269 ) { fLumiPerRun = 3.8392; }
  else if ( fRunNum == 296270 ) { fLumiPerRun = 1.51158; }
  else if ( fRunNum == 296273 ) { fLumiPerRun = 7.23985; }
  else if ( fRunNum == 296279 ) { fLumiPerRun = 0.405692; }
  else if ( fRunNum == 296280 ) { fLumiPerRun = 1.61122; }
  else if ( fRunNum == 296303 ) { fLumiPerRun = 2.01567; }
  else if ( fRunNum == 296304 ) { fLumiPerRun = 6.11939; }
  else if ( fRunNum == 296307 ) { fLumiPerRun = 2.93104; }
  else if ( fRunNum == 296309 ) { fLumiPerRun = 2.10255; }
  else if ( fRunNum == 296312 ) { fLumiPerRun = 2.12718; }
  else if ( fRunNum == 296377 ) { fLumiPerRun = 6.08141; }
  else if ( fRunNum == 296378 ) { fLumiPerRun = 5.34678; }
  else if ( fRunNum == 296379 ) { fLumiPerRun = 2.09698; }
  else if ( fRunNum == 296380 ) { fLumiPerRun = 2.88934; }
  else if ( fRunNum == 296381 ) { fLumiPerRun = 1.44081; }
  else if ( fRunNum == 296383 ) { fLumiPerRun = 1.51276; }
  else if ( fRunNum == 296414 ) { fLumiPerRun = 4.8759; }
  else if ( fRunNum == 296419 ) { fLumiPerRun = 2.74845; }
  else if ( fRunNum == 296420 ) { fLumiPerRun = 1.40973; }
  else if ( fRunNum == 296423 ) { fLumiPerRun = 1.57873; }
  else if ( fRunNum == 296424 ) { fLumiPerRun = 0.386229; }
  else if ( fRunNum == 296433 ) { fLumiPerRun = 4.46674; }
  else if ( fRunNum == 296472 ) { fLumiPerRun = 0.862237; }
  else if ( fRunNum == 296509 ) { fLumiPerRun = 2.99727; }
  else if ( fRunNum == 296510 ) { fLumiPerRun = 9.07207; }
  else if ( fRunNum == 296511 ) { fLumiPerRun = 2.5694; }
  else if ( fRunNum == 296514 ) { fLumiPerRun = 0.48925; }
  else if ( fRunNum == 296516 ) { fLumiPerRun = 0.613001; }
  else if ( fRunNum == 296547 ) { fLumiPerRun = 1.08314; }
  else if ( fRunNum == 296548 ) { fLumiPerRun = 1.37636; }
  else if ( fRunNum == 296549 ) { fLumiPerRun = 4.86114; }
  else if ( fRunNum == 296550 ) { fLumiPerRun = 3.96203; }
  else if ( fRunNum == 296551 ) { fLumiPerRun = 2.02168; }
  else if ( fRunNum == 296552 ) { fLumiPerRun = 0.483545; }
  else if ( fRunNum == 296553 ) { fLumiPerRun = 0.707355; }
  else if ( fRunNum == 296615 ) { fLumiPerRun = 1.56721; }
  else if ( fRunNum == 296616 ) { fLumiPerRun = 0.53916; }
  else if ( fRunNum == 296618 ) { fLumiPerRun = 2.28102; }
  else if ( fRunNum == 296619 ) { fLumiPerRun = 1.5574; }
  else if ( fRunNum == 296622 ) { fLumiPerRun = 0.703628; }
  else if ( fRunNum == 296623 ) { fLumiPerRun = 2.16132; }
  else if ( fRunNum == 296690 ) { fLumiPerRun = 6.55512; }
  else if ( fRunNum == 296691 ) { fLumiPerRun = 0.651063; }
  else if ( fRunNum == 296694 ) { fLumiPerRun = 5.02114; }
  else if ( fRunNum == 296749 ) { fLumiPerRun = 9.27577; }
  else if ( fRunNum == 296750 ) { fLumiPerRun = 8.21552; }
  else if ( fRunNum == 296781 ) { fLumiPerRun = 0.817883; }
  else if ( fRunNum == 296784 ) { fLumiPerRun = 3.0019; }
  else if ( fRunNum == 296785 ) { fLumiPerRun = 1.9085; }
  else if ( fRunNum == 296786 ) { fLumiPerRun = 0.753734; }
  else if ( fRunNum == 296787 ) { fLumiPerRun = 3.24915; }
  else if ( fRunNum == 296791 ) { fLumiPerRun = 0.757278; }
  else if ( fRunNum == 296793 ) { fLumiPerRun = 1.39522; }
  else if ( fRunNum == 296794 ) { fLumiPerRun = 3.13783; }
  else if ( fRunNum == 296799 ) { fLumiPerRun = 2.77608; }
  else if ( fRunNum == 296836 ) { fLumiPerRun = 1.51136; }
  else if ( fRunNum == 296838 ) { fLumiPerRun = 0.542885; }
  else if ( fRunNum == 296839 ) { fLumiPerRun = 2.95419; }
  else if ( fRunNum == 296848 ) { fLumiPerRun = 2.19145; }
  else if ( fRunNum == 296849 ) { fLumiPerRun = 11.8282; }
  else if ( fRunNum == 296850 ) { fLumiPerRun = 2.8369; }
  else if ( fRunNum == 296851 ) { fLumiPerRun = 0.890021; }
  else if ( fRunNum == 296852 ) { fLumiPerRun = 0.948014; }
  else if ( fRunNum == 296890 ) { fLumiPerRun = 6.96187; }
  else if ( fRunNum == 296894 ) { fLumiPerRun = 4.47163; }
  else if ( fRunNum == 296899 ) { fLumiPerRun = 2.11682; }
  else if ( fRunNum == 296900 ) { fLumiPerRun = 2.78312; }
  else if ( fRunNum == 296903 ) { fLumiPerRun = 1.03903; }
  else if ( fRunNum == 296930 ) { fLumiPerRun = 1.42317; }
  else if ( fRunNum == 296931 ) { fLumiPerRun = 0.528865; }
  else if ( fRunNum == 296932 ) { fLumiPerRun = 1.19931; }
  else if ( fRunNum == 296934 ) { fLumiPerRun = 2.62967; }
  else if ( fRunNum == 296935 ) { fLumiPerRun = 4.49271; }
  else if ( fRunNum == 296938 ) { fLumiPerRun = 1.65053; }
  else if ( fRunNum == 296941 ) { fLumiPerRun = 2.93132; }
  else if ( fRunNum == 296966 ) { fLumiPerRun = 3.37409; }
  else if ( fRunNum == 296967 ) { fLumiPerRun = 0.804596; }
  else if ( fRunNum == 296968 ) { fLumiPerRun = 3.18308; }
  else if ( fRunNum == 296969 ) { fLumiPerRun = 1.88784; }
  else if ( fRunNum == 296971 ) { fLumiPerRun = 0.6895; }
  else if ( fRunNum == 296975 ) { fLumiPerRun = 7.47159; }
  else if ( fRunNum == 296976 ) { fLumiPerRun = 1.11722; }
  else if ( fRunNum == 296979 ) { fLumiPerRun = 1.09943; }
  else if ( fRunNum == 297029 ) { fLumiPerRun = 7.20294; }
  else if ( fRunNum == 297031 ) { fLumiPerRun = 5.94515; }
  else if ( fRunNum == 297035 ) { fLumiPerRun = 1.30617; }
  else if ( fRunNum == 297085 ) { fLumiPerRun = 0.977753; }
  else if ( fRunNum == 297117 ) { fLumiPerRun = 2.34892; }
  else if ( fRunNum == 297118 ) { fLumiPerRun = 2.44282; }
  else if ( fRunNum == 297119 ) { fLumiPerRun = 2.68704; }
  else if ( fRunNum == 297123 ) { fLumiPerRun = 3.3714; }
  else if ( fRunNum == 297124 ) { fLumiPerRun = 0.639463; }
  else if ( fRunNum == 297128 ) { fLumiPerRun = 2.41074; }
  else if ( fRunNum == 297129 ) { fLumiPerRun = 2.82995; }
  else if ( fRunNum == 297132 ) { fLumiPerRun = 2.81789; }
  else if ( fRunNum == 297133 ) { fLumiPerRun = 1.16976; }
  else if ( fRunNum == 297193 ) { fLumiPerRun = 7.64123; }
  else if ( fRunNum == 297194 ) { fLumiPerRun = 9.86729; }
  else if ( fRunNum == 297196 ) { fLumiPerRun = 2.1255; }
  else if ( fRunNum == 297218 ) { fLumiPerRun = 6.39259; }
  else if ( fRunNum == 297219 ) { fLumiPerRun = 9.29989; }
  else if ( fRunNum == 297221 ) { fLumiPerRun = 2.83193; }
  else if ( fRunNum == 297222 ) { fLumiPerRun = 1.69325; }
  else if ( fRunNum == 297278 ) { fLumiPerRun = 0.601609; }
  else if ( fRunNum == 297310 ) { fLumiPerRun = 0.670071; }
  else if ( fRunNum == 297312 ) { fLumiPerRun = 2.40205; }
  else if ( fRunNum == 297315 ) { fLumiPerRun = 7.93229; }
  else if ( fRunNum == 297317 ) { fLumiPerRun = 4.31559; }
  else if ( fRunNum == 297363 ) { fLumiPerRun = 1.89669; }
  else if ( fRunNum == 297366 ) { fLumiPerRun = 2.05394; }
  else if ( fRunNum == 297367 ) { fLumiPerRun = 3.11285; }
  else if ( fRunNum == 297372 ) { fLumiPerRun = 3.22421; }
  else if ( fRunNum == 297379 ) { fLumiPerRun = 6.92989; }
  else if ( fRunNum == 297380 ) { fLumiPerRun = 1.46125; }
  else if ( fRunNum == 297405 ) { fLumiPerRun = 0.51899; }
  else if ( fRunNum == 297408 ) { fLumiPerRun = 4.05969; }
  else if ( fRunNum == 297413 ) { fLumiPerRun = 2.99189; }
  else if ( fRunNum == 297414 ) { fLumiPerRun = 2.21433; }
  else if ( fRunNum == 297415 ) { fLumiPerRun = 6.76725; }
  else if ( fRunNum == 297441 ) { fLumiPerRun = 4.92805; }
  else if ( fRunNum == 297442 ) { fLumiPerRun = 1.98748; }
  else if ( fRunNum == 297446 ) { fLumiPerRun = 8.1332; }
  else if ( fRunNum == 297450 ) { fLumiPerRun = 1.95205; }
  else if ( fRunNum == 297451 ) { fLumiPerRun = 1.33275; }
  else if ( fRunNum == 297452 ) { fLumiPerRun = 1.15124; }
  else if ( fRunNum == 297479 ) { fLumiPerRun = 7.71966; }
  else if ( fRunNum == 297481 ) { fLumiPerRun = 9.55402; }
  else if ( fRunNum == 297483 ) { fLumiPerRun = 1.99555; }
  else if ( fRunNum == 297512 ) { fLumiPerRun = 1.5838; }
  else if ( fRunNum == 297537 ) { fLumiPerRun = 1.80636; }
  else if ( fRunNum == 297540 ) { fLumiPerRun = 0.628094; }
  else if ( fRunNum == 297541 ) { fLumiPerRun = 3.96702; }
  else if ( fRunNum == 297542 ) { fLumiPerRun = 1.5516; }
  else if ( fRunNum == 297544 ) { fLumiPerRun = 7.29158; }
  else if ( fRunNum == 297558 ) { fLumiPerRun = 0.478978; }
  else if ( fRunNum == 297588 ) { fLumiPerRun = 5.26723; }
  else if ( fRunNum == 297590 ) { fLumiPerRun = 3.14808; }
  else if ( fRunNum == 297595 ) { fLumiPerRun = 0.0; }
  else                          { fLumiPerRun = 0.0; }


}
//_____________________________________________________________________________
/* - Computes the efficiency of the MC on a
 * - run-by-run basis.
 * -
 */
void fitEfficiencyMC(){

  // TFile* fileList = new TFile("AnalysisResultsMCcohJPsiLHC18l7.root");        // same settings as polarisation analysis
  TFile* fileList = new TFile("AnalysisResultsCoherentMC_ADrestrictions.root");  // same settings as CheckAD
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
  // TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_0");
  // TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_0");
  TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_1");
  TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_1");
  // TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_2");
  // TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_2");
  // TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_3");
  // TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_3");
  // TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_4");
  // TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_4");
  // TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunRapidityH_5");
  // TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunRapidityH_5");

  fEfficiencyPerRunH  ->Sumw2();
  fMCEfficiencyPerRunH->Sumw2();

  TH1F* RealEfficiency = (TH1F*) fEfficiencyPerRunH->Clone("RealEfficiency");
  TH1F* MCEfficiency = (TH1F*) fMCEfficiencyPerRunH->Clone("MCEfficiency");
  // RealEfficiency->Divide(fMCEfficiencyPerRunH);
  // RealEfficiency->Draw("ep");

  Double_t ComputedEfficiency[300];
  Double_t ValueEfficiency[300];
  Double_t ValueEfficiency2[300];
  Double_t Error1[300];
  Double_t Error2[300];
  TString  runLabel[300];
  for ( Int_t i = 0; i < 300; i++ ) ComputedEfficiency[i] = 0;
  for ( Int_t i = 0; i < 300; i++ ) ValueEfficiency[i] = 0;
  for ( Int_t i = 0; i < 300; i++ ) ValueEfficiency2[i] = 0;

  Int_t counter = 0;
  cout << "fEfficiencyPerRunH : " << endl;
  for ( Int_t iLoop = 1; iLoop <= RealEfficiency->GetNbinsX(); iLoop++ ) {
    if ( iLoop > 229 ) break;
    TString label          = RealEfficiency->GetXaxis()->GetBinLabel(iLoop);
    runLabel[iLoop]        = RealEfficiency->GetXaxis()->GetBinLabel(iLoop);
    // if ( atoi( label.Data() ) == 266615) continue;
    ValueEfficiency[iLoop] = RealEfficiency->GetBinContent(iLoop);
    if ( ValueEfficiency[iLoop] != 0 ) {
      Error1[iLoop]          = RealEfficiency->GetBinError(iLoop) / ValueEfficiency[iLoop];
    } else {
      Error1[iLoop]          = 0;
    }
    // ValueEfficiency[iLoop] = fEfficiencyPerRunH->GetXaxis()->GetBinContent(iLoop);
    counter++;
    cout <<   atoi( label.Data() ) << endl;
    fRunNum = atoi( label.Data() );
    SetLuminosityCap();
    Double_t weight = fLumiPerRun / fLumi;
    for ( Int_t iLoop2 = 1; iLoop2 <= MCEfficiency->GetNbinsX(); iLoop2++ ) {
      if ( iLoop2 > 280 ) break;
      TString label2 = MCEfficiency->GetXaxis()->GetBinLabel(iLoop2);
      // if ( atoi(label2.Data()) == atoi(label.Data()) ) {
      if ( atoi( MCEfficiency->GetXaxis()->GetBinLabel(iLoop2) ) == atoi( RealEfficiency->GetXaxis()->GetBinLabel(iLoop) ) ) {
        cout << atoi( MCEfficiency->GetXaxis()->GetBinLabel(iLoop2) ) << endl;
        cout << atoi( RealEfficiency->GetXaxis()->GetBinLabel(iLoop) ) << endl;
        // cout << "OK : " << counter << endl;
        // if ( ValueEfficiency != 0 && fMCEfficiencyPerRunH->GetXaxis()->GetBinContent(iLoop2) != 0 ) {
        if ( ValueEfficiency[iLoop] != 0 && MCEfficiency->GetBinContent(iLoop2) != 0 ) {
          // ComputedEfficiency[iLoop] = (Double_t)ValueEfficiency[iLoop]; // / (Double_t)MCEfficiency->GetBinContent(iLoop2);
          ValueEfficiency2[iLoop] = MCEfficiency->GetBinContent(iLoop2);
          Error2[iLoop]           = MCEfficiency->GetBinError(iLoop2) / ValueEfficiency2[iLoop];
          ComputedEfficiency[iLoop] = (Double_t)ValueEfficiency[iLoop] / (Double_t)ValueEfficiency2[iLoop];
          // ComputedEfficiency[iLoop] = ValueEfficiency[iLoop] / fMCEfficiencyPerRunH->GetXaxis()->GetBinContent(iLoop2);
        } else {
          ComputedEfficiency[iLoop] = 0;
        }
      } else {
        continue;
      }

      // cout << label2.Data() << endl;
      // counter++;
    }


    GlobalAxE += (weight * ComputedEfficiency[iLoop]);
    cout << "Global AxE   = " << GlobalAxE << endl;
    cout << "Contribution ="  << (weight * ComputedEfficiency[iLoop]) << endl;

  }
  cout << "counter = " << counter << endl;


  // for ( Int_t i = 0; i < 200 ; i++ ) {
  //   cout << "Run " << i << " = " << ComputedEfficiency[i] << endl;
  // }
  //
  // for ( Int_t i = 0; i < 200 ; i++ ) {
  //   cout << "Run tr1 " << i << " = " << ValueEfficiency[i] << endl;
  //   cout << "Run tr2 " << i << " = " << ValueEfficiency2[i] << endl;
  // }


  // counter = 0;
  // cout << "fMCEfficiencyPerRunH : " << endl;
  // for ( Int_t iLoop = 1; iLoop <= fMCEfficiencyPerRunH->GetNbinsX(); iLoop++ ) {
  //   TString label = fMCEfficiencyPerRunH->GetXaxis()->GetBinLabel(iLoop);
  //   cout << label.Data() << endl;
  //   // counter++;
  // }
  // cout << "counter = " << counter << endl;
  TCanvas* EffCanvas = new TCanvas("EffCanvas","EffCanvas",900,800);
  TH1F* Eff = new TH1F( "eff" , "eff", 200, -0.5, 199.5 );
  // Eff->SetStats(0);
  Eff->SetFillColor(38);
  Eff->LabelsDeflate();
  for ( Int_t i = 0; i < 229; i++) {
    // Eff->Fill( i, ComputedEfficiency[i] );
    Eff->Fill( runLabel[i].Data() , ComputedEfficiency[i] );
    Eff->SetBinError( i+1 , (Error1[i] + Error2[i]) * ComputedEfficiency[i]);
  }
  Eff->Draw("ep");

  cout << "Global AxE = " << GlobalAxE << endl;

  // TFile f("pngResults/efficiency.root", "recreate");
  // // RealEfficiency->Write();
  // f.Close();
}
