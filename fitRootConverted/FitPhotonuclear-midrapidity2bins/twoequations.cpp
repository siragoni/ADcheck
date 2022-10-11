// c++ headers
#include <iostream>
#include <fstream>
#include <iomanip>

// root headers
#include <Rtypes.h>
#include <TMath.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;

void twoequations(){

    double xnxn[3] = {0.0830, 0.1693, 0.2681};
    double xnon[3] = {0.0993, 0.1681, 0.2359};
    double onxn[3] = {0.0430, 0.1092, 0.2302};
    double onon[3] = {1.5860, 2.3150, 2.6570};

    double xnxn_stat[3] = {0.0119, 0.0120, 0.0249};
    double xnon_stat[3] = {0.0078, 0.0098, 0.0208};
    double onxn_stat[3] = {0.0035, 0.0053, 0.0126};
    double onon_stat[3] = {0.0489, 0.0333, 0.0760};

    double xnxn_sysperc[3] = {0.1348, 0.1364, 0.1336};
    double xnon_sysperc[3] = {0.0746, 0.0778, 0.0814};
    double onxn_sysperc[3] = {0.3302, 0.1863, 0.1088};
    double onon_sysperc[3] = {0.0731, 0.0731, 0.0734};

    double xnxn_sys[3] = {0., 0., 0.};
    double xnon_sys[3] = {0., 0., 0.};
    double onxn_sys[3] = {0., 0., 0.};
    double onon_sys[3] = {0., 0., 0.};

    for (size_t i = 0; i < 3; i++)
    {
        xnxn_sys[i] = xnxn_sysperc[i]*xnxn[i];
        xnon_sys[i] = xnon_sysperc[i]*xnon[i];
        onxn_sys[i] = onxn_sysperc[i]*onxn[i];
        onon_sys[i] = onon_sysperc[i]*onon[i];
    }






    double photon_one_xnxn[3] = {6.5199, 6.5209, 6.5225};
    double photon_two_xnxn[3] = {0.4250, 1.3546, 2.7674};

    double photon_one_xnon[3] = {18.31200*0.5, 18.3120*0.5, 18.3120*0.5};
    double photon_two_xnon[3] = { 0.52656*0.5,  1.9936*0.5,  4.8253*0.5};
    
    double photon_one_onxn[3] = {18.31200*0.5, 18.3120*0.5, 18.3120*0.5};
    double photon_two_onxn[3] = { 0.52656*0.5,  1.9936*0.5,  4.8253*0.5};

    double photon_one_onon[3] = {178.2800, 162.7500, 147.2300};
    double photon_two_onon[3] = {  0.2149,   1.0973,   3.6998};



    double sigma_onon_xnon_plus[3]  = {0,0,0};
    double sigma_onon_xnon_minus[3] = {0,0,0};

    double sigma_onon_xnxn_plus[3]  = {0,0,0};
    double sigma_onon_xnxn_minus[3] = {0,0,0};

    double sigma_onon_onxn_plus[3]  = {0,0,0};
    double sigma_onon_onxn_minus[3] = {0,0,0};




    double sigma_onon_xnon_plus_stat[3]  = {0,0,0};
    double sigma_onon_xnon_minus_stat[3] = {0,0,0};

    double sigma_onon_xnxn_plus_stat[3]  = {0,0,0};
    double sigma_onon_xnxn_minus_stat[3] = {0,0,0};

    double sigma_onon_onxn_plus_stat[3]  = {0,0,0};
    double sigma_onon_onxn_minus_stat[3] = {0,0,0};


    double sigma_onon_xnon_plus_sys[3]  = {0,0,0};
    double sigma_onon_xnon_minus_sys[3] = {0,0,0};

    double sigma_onon_xnxn_plus_sys[3]  = {0,0,0};
    double sigma_onon_xnxn_minus_sys[3] = {0,0,0};

    double sigma_onon_onxn_plus_sys[3]  = {0,0,0};
    double sigma_onon_onxn_minus_sys[3] = {0,0,0};



    for (size_t i = 0; i < 3; i++)
    {
        sigma_onon_xnon_plus[i] = ( onon[i]/photon_two_onon[i] - xnon[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus[i] = ( onon[i]/photon_two_onon[i] - xnxn[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus[i] = ( onon[i]/photon_two_onon[i] - onxn[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );
    
        sigma_onon_xnon_minus[i] = ( onon[i]/photon_one_onon[i] - xnon[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus[i] = ( onon[i]/photon_one_onon[i] - xnxn[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus[i] = ( onon[i]/photon_one_onon[i] - onxn[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );

    
    


        sigma_onon_xnon_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + xnon_stat[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + xnxn_stat[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + onxn_stat[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );

        sigma_onon_xnon_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + xnon_stat[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + xnxn_stat[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + onxn_stat[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );




        sigma_onon_xnon_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - xnon_sys[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - xnxn_sys[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - onxn_sys[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );
    
        sigma_onon_xnon_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - xnon_sys[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - xnxn_sys[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - onxn_sys[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );

    
    }
    

    
    for (size_t i = 0; i < 3; i++)
    {
        cout << "sigma_onon_xnon_plus[" << i << "]  = " <<  sigma_onon_xnon_plus[i]  << " +/- " << sigma_onon_xnon_plus_stat[i]  << " (stat.) +/- " << sigma_onon_xnon_plus_sys[i]  << " (sys.)" << endl;
        cout << "sigma_onon_xnon_minus[" << i << "] = " <<  sigma_onon_xnon_minus[i] << " +/- " << sigma_onon_xnon_minus_stat[i] << " (stat.) +/- " << sigma_onon_xnon_minus_sys[i] << " (sys.)" << endl;
    }
    for (size_t i = 0; i < 3; i++)
    {
        cout << "sigma_onon_xnxn_plus[" << i << "]  = " <<  sigma_onon_xnxn_plus[i]  << " +/- " << sigma_onon_xnxn_plus_stat[i]  << " (stat.) +/- " << sigma_onon_xnxn_plus_sys[i]  << " (sys.)" << endl;
        cout << "sigma_onon_xnxn_minus[" << i << "] = " <<  sigma_onon_xnxn_minus[i] << " +/- " << sigma_onon_xnxn_minus_stat[i] << " (stat.) +/- " << sigma_onon_xnxn_minus_sys[i] << " (sys.)" << endl;
    }
    for (size_t i = 0; i < 3; i++)
    {
        cout << "sigma_onon_onxn_plus[" << i << "]  = " <<  sigma_onon_onxn_plus[i]  << " +/- " << sigma_onon_onxn_plus_stat[i] << " (stat.) +/- " << sigma_onon_onxn_plus_sys[i]  << " (sys.)" << endl;
        cout << "sigma_onon_onxn_minus[" << i << "] = " <<  sigma_onon_onxn_minus[i] << " +/- " << sigma_onon_onxn_plus_sys[i]  << " (stat.) +/- " << sigma_onon_onxn_minus_sys[i] << " (sys.)" << endl;
    }
  



}








//____________
void twoequations-midrapidity(){

    TFile* file2 = new TFile("../../Michal-Broz-xsec/xSection_Cent2.root");
    TGraphAsymmErrors* Cent0N0N = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0n0n");
    TGraphAsymmErrors* Cent0NXN = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_0nXn");
    TGraphAsymmErrors* CentXNXN = (TGraphAsymmErrors*) file2->Get("gXsection_Cent_Syst_XnXn");
    TGraphErrors* Cent0N0N_2 = (TGraphErrors*) file2->Get("gXsection_Cent_Stat_0n0n");
    TGraphErrors* Cent0NXN_2 = (TGraphErrors*) file2->Get("gXsection_Cent_Stat_0nXn");
    TGraphErrors* CentXNXN_2 = (TGraphErrors*) file2->Get("gXsection_Cent_Stat_XnXn");




    double xnxn[2] = {0,0};
    double xnon[2] = {0,0};
    double onon[2] = {0,0};

    double xnxn_stat[2] = {0,0};
    double xnon_stat[2] = {0,0};
    double onon_stat[2] = {0,0};

    // double xnxn_sysperc[2] = {0,0};
    // double xnon_sysperc[2] = {0,0};
    // double onon_sysperc[2] = {0,0};

    double xnxn_sys[2] = {0,0};
    double xnon_sys[2] = {0,0};
    double onon_sys[2] = {0,0};

    for (size_t i = 0; i < 2; i++)
    {
        // xnxn_sys[i] = xnxn_sysperc[i]*xnxn[i];
        // xnon_sys[i] = xnon_sysperc[i]*xnon[i];
        // onon_sys[i] = onon_sysperc[i]*onon[i];
        xnxn[i] = CentXNXN_2->GetY()[i];
        xnon[i] = Cent0NXN_2->GetY()[i];
        onon[i] = Cent0N0N_2->GetY()[i];

        xnxn_stat[i] = CentXNXN_2->GetErrorY(i);
        xnon_stat[i] = Cent0NXN_2->GetErrorY(i);
        onon_stat[i] = Cent0N0N_2->GetErrorY(i);

    }






    double photon_one_xnxn[2] = {6.5199, 6.5209, 6.5225};
    double photon_two_xnxn[2] = {0.4250, 1.3546, 2.7674};

    double photon_one_xnon[2] = {18.31200*0.5, 18.3120*0.5, 18.3120*0.5};
    double photon_two_xnon[2] = { 0.52656*0.5,  1.9936*0.5,  4.8253*0.5};
    
    double photon_one_onon[2] = {178.2800, 162.7500, 147.2300};
    double photon_two_onon[2] = {  0.2149,   1.0973,   3.6998};



    double sigma_onon_xnon_plus[2]  = {0,0};
    double sigma_onon_xnon_minus[2] = {0,0};

    double sigma_onon_xnxn_plus[2]  = {0,0};
    double sigma_onon_xnxn_minus[2] = {0,0};





    double sigma_onon_xnon_plus_stat[2]  = {0,0};
    double sigma_onon_xnon_minus_stat[2] = {0,0};

    double sigma_onon_xnxn_plus_stat[2]  = {0,0};
    double sigma_onon_xnxn_minus_stat[2] = {0,0};

    double sigma_onon_onxn_plus_stat[2]  = {0,0};
    double sigma_onon_onxn_minus_stat[2] = {0,0};


    double sigma_onon_xnon_plus_sys[2]  = {0,0};
    double sigma_onon_xnon_minus_sys[2] = {0,0};

    double sigma_onon_xnxn_plus_sys[2]  = {0,0};
    double sigma_onon_xnxn_minus_sys[2] = {0,0};

    double sigma_onon_onxn_plus_sys[2]  = {0,0};
    double sigma_onon_onxn_minus_sys[2] = {0,0};



    for (size_t i = 0; i < 2; i++)
    {
        sigma_onon_xnon_plus[i] = ( onon[i]/photon_two_onon[i] - xnon[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus[i] = ( onon[i]/photon_two_onon[i] - xnxn[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus[i] = ( onon[i]/photon_two_onon[i] - onxn[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );
    
        sigma_onon_xnon_minus[i] = ( onon[i]/photon_one_onon[i] - xnon[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus[i] = ( onon[i]/photon_one_onon[i] - xnxn[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus[i] = ( onon[i]/photon_one_onon[i] - onxn[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );

    
    


        sigma_onon_xnon_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + xnon_stat[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + xnxn_stat[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus_stat[i] = ( onon_stat[i]/photon_two_onon[i] + onxn_stat[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );

        sigma_onon_xnon_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + xnon_stat[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + xnxn_stat[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus_stat[i] = ( onon_stat[i]/photon_one_onon[i] + onxn_stat[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );




        sigma_onon_xnon_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - xnon_sys[i]/photon_two_xnon[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnon[i]/photon_two_xnon[i] );
        sigma_onon_xnxn_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - xnxn_sys[i]/photon_two_xnxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_xnxn[i]/photon_two_xnxn[i] );
        sigma_onon_onxn_plus_sys[i] = ( onon_sys[i]/photon_two_onon[i] - onxn_sys[i]/photon_two_onxn[i] ) / ( photon_one_onon[i]/photon_two_onon[i] - photon_one_onxn[i]/photon_two_onxn[i] );
    
        sigma_onon_xnon_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - xnon_sys[i]/photon_one_xnon[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnon[i]/photon_one_xnon[i] );
        sigma_onon_xnxn_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - xnxn_sys[i]/photon_one_xnxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_xnxn[i]/photon_one_xnxn[i] );
        sigma_onon_onxn_minus_sys[i] = ( onon_sys[i]/photon_one_onon[i] - onxn_sys[i]/photon_one_onxn[i] ) / ( photon_two_onon[i]/photon_one_onon[i] - photon_two_onxn[i]/photon_one_onxn[i] );

    
    }
    

    
    for (size_t i = 0; i < 2; i++)
    {
        cout << "sigma_onon_xnon_plus[" << i << "]  = " <<  sigma_onon_xnon_plus[i]  << " +/- " << sigma_onon_xnon_plus_stat[i]  << " (stat.) +/- " << sigma_onon_xnon_plus_sys[i]  << " (sys.)" << endl;
        cout << "sigma_onon_xnon_minus[" << i << "] = " <<  sigma_onon_xnon_minus[i] << " +/- " << sigma_onon_xnon_minus_stat[i] << " (stat.) +/- " << sigma_onon_xnon_minus_sys[i] << " (sys.)" << endl;
    }
    for (size_t i = 0; i < 2; i++)
    {
        cout << "sigma_onon_xnxn_plus[" << i << "]  = " <<  sigma_onon_xnxn_plus[i]  << " +/- " << sigma_onon_xnxn_plus_stat[i]  << " (stat.) +/- " << sigma_onon_xnxn_plus_sys[i]  << " (sys.)" << endl;
        cout << "sigma_onon_xnxn_minus[" << i << "] = " <<  sigma_onon_xnxn_minus[i] << " +/- " << sigma_onon_xnxn_minus_stat[i] << " (stat.) +/- " << sigma_onon_xnxn_minus_sys[i] << " (sys.)" << endl;
    }
  



}




