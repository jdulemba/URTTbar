#include "Analyses/URTTbar/interface/Alpha_Corrections.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "URAnalysis/AnalysisFW/interface/Wards.h"


Alpha_Corrections::Alpha_Corrections(){

    URParser &parser = URParser::instance();
    parser.addCfgParameter<string>("alpha_correction", "fit", "");
    parser.parseArguments();

        //get filename that has fit parameters
    string fname = parser.getCfgPar<string>("alpha_correction", "fit");
    TFile alpha_file(DataFile(fname).path().c_str());

    HistoOwnershipWard hward;
    DirectoryWard dward;

    THad_E_Mtt_1d_ = get_from<TH2D>(alpha_file, "THad_E_Mtt_1d", "THad_E_Mtt_1d_clone");
    THad_E_Mtt_2d_ = get_from<TH2D>(alpha_file, "THad_E_Mtt_2d", "THad_E_Mtt_2d_clone");
    THad_P_Mtt_1d_ = get_from<TH2D>(alpha_file, "THad_P_Mtt_1d", "THad_P_Mtt_1d_clone");
    THad_P_Mtt_2d_ = get_from<TH2D>(alpha_file, "THad_P_Mtt_2d", "THad_P_Mtt_2d_clone");
    THad_M_Mtt_1d_ = get_from<TH2D>(alpha_file, "THad_M_Mtt_1d", "THad_M_Mtt_1d_clone");
    THad_M_Mtt_2d_ = get_from<TH2D>(alpha_file, "THad_M_Mtt_2d", "THad_M_Mtt_2d_clone");
    THad_E_All_1d_ = get_from<TH1D>(alpha_file, "THad_E_All_1d", "THad_E_All_1d_clone");
    THad_E_All_2d_ = get_from<TH1D>(alpha_file, "THad_E_All_2d", "THad_E_All_2d_clone");
    THad_P_All_1d_ = get_from<TH1D>(alpha_file, "THad_P_All_1d", "THad_P_All_1d_clone");
    THad_P_All_2d_ = get_from<TH1D>(alpha_file, "THad_P_All_2d", "THad_P_All_2d_clone");
    THad_M_All_1d_ = get_from<TH1D>(alpha_file, "THad_M_All_1d", "THad_M_All_1d_clone");
    THad_M_All_2d_ = get_from<TH1D>(alpha_file, "THad_M_All_2d", "THad_M_All_2d_clone");

}

Alpha_Corrections::~Alpha_Corrections()
{}

double Alpha_Corrections::alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range)
{

    if( fit_var != "E" && fit_var != "P" && fit_var!= "M" && fit_var != "None" ){
        Logger::log().error() << "Only energy (E), momentum (P), and mass (M) variables are supported!" << endl;
        throw 42;
    }

    double alpha = 1.; // reset alpha

    if( fit_var != "None" ){
        if( fit_degree != "1d" && fit_degree != "2d" ){
            Logger::log().error() << "Only linear (1D) and quadratic (2D) fits are supported!" << endl;
            throw 42;
        }
        if( fit_range != "All" && fit_range != "Mtt" ){
            Logger::log().error() << "Mtt fit range not supported!" << endl;
            throw 42;
        }

        string hname = "THad_"+fit_var+"_"+fit_range+"_"+fit_degree;

        //cout << "hist name: " << hname << endl;
       
            // find alpha from TH1D hists 
        if( hname == "THad_E_All_1d" ) alpha = get_1D_corr( THad_E_All_1d_, 173.1/perm.THad().M() );
        else if( hname == "THad_E_All_2d" ) alpha = get_1D_corr( THad_E_All_2d_, 173.1/perm.THad().M() );
        else if( hname == "THad_P_All_1d" ) alpha = get_1D_corr( THad_P_All_1d_, 173.1/perm.THad().M() );
        else if( hname == "THad_P_All_2d" ) alpha = get_1D_corr( THad_P_All_2d_, 173.1/perm.THad().M() );
        else if( hname == "THad_M_All_1d" ) alpha = get_1D_corr( THad_M_All_1d_, 173.1/perm.THad().M() );
        else if( hname == "THad_M_All_2d" ) alpha = get_1D_corr( THad_M_All_2d_, 173.1/perm.THad().M() );

            // find alpha from TH2D hists 
        else if( hname == "THad_E_Mtt_1d" ) alpha = get_2D_corr( THad_E_Mtt_1d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else if( hname == "THad_E_Mtt_2d" ) alpha = get_2D_corr( THad_E_Mtt_2d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else if( hname == "THad_P_Mtt_1d" ) alpha = get_2D_corr( THad_P_Mtt_1d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else if( hname == "THad_P_Mtt_2d" ) alpha = get_2D_corr( THad_P_Mtt_2d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else if( hname == "THad_M_Mtt_1d" ) alpha = get_2D_corr( THad_M_Mtt_1d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else if( hname == "THad_M_Mtt_2d" ) alpha = get_2D_corr( THad_M_Mtt_2d_, 173.1/perm.THad().M(), perm.LVect().M() );
        else{
            Logger::log().error() << "Combination of input variables not supported! Hist name that crashed: " << hname  << endl;
            throw 42;
        }
    }

    return alpha;

}

double Alpha_Corrections::get_1D_corr( std::shared_ptr<TH1> h, double mthad ) const
{
    int binx = h->GetXaxis()->FindFixBin( mthad );
    binx = TMath::Max(binx, 1);
    binx = TMath::Min(binx, h->GetNbinsX());

    double binval = h->GetBinContent( binx );

    return binval;
}

double Alpha_Corrections::get_2D_corr( std::shared_ptr<TH2> h, double mthad, double mtt ) const
{
    int binx = h->GetXaxis()->FindFixBin( mthad );
    binx = TMath::Max(binx, 1);
    binx = TMath::Min(binx, h->GetNbinsX());

    int biny = h->GetYaxis()->FindFixBin( mtt );
    biny = TMath::Max(biny, 1);
    biny = TMath::Min(biny, h->GetNbinsY());

    double binval = h->GetBinContent( binx, biny );

    return binval;
}


TLorentzVector Alpha_Corrections::Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range) // scaled reco thad 4-vec
{
    double alpha_corr = alpha(perm, fit_degree, fit_var, fit_range);
    TLorentzVector alpha_thad(alpha_corr*perm.THad().Px(), alpha_corr*perm.THad().Py(), alpha_corr*perm.THad().Pz(), alpha_corr*perm.THad().E());
    return alpha_thad;
}


double Alpha_Corrections::alpha_thad_cthstar(Permutation &perm, string fit_degree, string fit_var, string fit_range)
{
    TLorentzVector alpha_thad = Alpha_THad(perm, fit_degree, fit_var, fit_range);

    TLorentzVector reco_corr_ttang = alpha_thad+perm.TLep();
    TLorentzVector reco_ttcm_corr_thad = alpha_thad;
    reco_ttcm_corr_thad.Boost(-1*reco_corr_ttang.BoostVector());
    double reco_corr_thad_cth = reco_corr_ttang.Vect().Unit().Dot(reco_ttcm_corr_thad.Vect().Unit());
    //double reco_corr_labframe_thad_cth = reco_corr_ttang.Vect().Unit().Dot(alpha_thad.Vect().Unit());

    return reco_corr_thad_cth;
}

//double Alpha_Corrections::median_alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range)
//{
////    if( fit_var == "None" ){
////        //Logger::log().debug() << "Alpha correction set to 1." << endl;
//        median_alpha_ = 1.;
////    }
////    else{
////        if( fit_var != "E" && fit_var != "P" && fit_var!= "M" && fit_var != "None" ){
////            Logger::log().error() << "Only energy (E), momentum (P), and mass (M) variables are supported!" << endl;
////            throw 42;
////        }
////    
////        string mtt_range;
////        double MTHad;
////    
////        if( fit_range == "Mtt" ){
////            if( 200. <= perm.LVect().M() && perm.LVect().M() < 350. ){
////                mtt_range = "Mttbar200to350";
////                    //values taken from http://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2016Data/JetpTcut20/LeadPt0/Only_Alpha_Correction/Full/ttJets/3J_Event_Plots/Lost_BP/Alpha_Correction/THad_E/y_proj/Alpha_THad_E_Mttbar200to350_yprojections.png
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.05;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.19;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.35;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.51;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.69;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.89;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 2.07;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 2.29;
////            }
////            else if( 350. <= perm.LVect().M() && perm.LVect().M() < 400. ){
////                mtt_range = "Mttbar350to400";
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.03;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.15;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.27;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.39;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.49;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.59;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.67;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.73;
////            }
////            else if( 400. <= perm.LVect().M() && perm.LVect().M() < 500. ){
////                mtt_range = "Mttbar400to500";
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.01;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.11;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.21;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.31;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.39;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.43;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.47;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.51;
////            }
////            else if( 500. <= perm.LVect().M() && perm.LVect().M() < 700. ){
////                mtt_range = "Mttbar500to700";
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.07;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.17;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.23;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.29;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.33;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.35;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.37;
////            }
////            else if( 700. <= perm.LVect().M() && perm.LVect().M() < 1000. ){
////                mtt_range = "Mttbar700to1000";
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.05;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.13;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.17;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.23;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.23;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.27;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.29;
////            }
////            else if( perm.LVect().M() >= 1000. ){
////                mtt_range = "Mttbar1000toInf";
////                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
////                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.03;
////                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.09;
////                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.15;
////                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.19;
////                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.23;
////                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.21;
////                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.23;
////            }
////        }
////        else if( fit_range == "All" ){
////            mtt_range = fit_range;
////                // values taken from http://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2016Data/JetpTcut20/LeadPt0/Only_Alpha_Correction/Full/ttJets/3J_Event_Plots/Lost_BP/Alpha_Correction/THad_E/y_proj/Alpha_THad_E_yprojections.png
////            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.01;
////            if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.11;
////            if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.23;
////            if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.35;
////            if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.47;
////            if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.57;
////            if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.65;
////            if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.71;
////        }
////        else{
////            Logger::log().error() << "Fit range not supported!" << endl;
////            throw 42;
////        }
////    
////        vector<double> pars = jsonFile_["THad_"+fit_var][mtt_range][fit_degree];
////   
////        if( fit_degree == "1d" ) median_alpha_ = pars[0]*( MTHad ) + pars[1];
////        else if ( fit_degree == "2d" ) median_alpha_ = pars[0]*( pow(MTHad,2) ) + pars[1]*( MTHad ) + pars[2];
////        else{
////            Logger::log().error() << "Only linear (1D) and quadratic (2D) fits are supported!" << endl;
////            throw 42;
////        }
////    }
////
//    return median_alpha_;
//
//}
//
//TLorentzVector Alpha_Corrections::median_Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range) // scaled reco thad 4-vec
//{
//    median_alpha_ = median_alpha(perm, fit_degree, fit_var, fit_range);
//    Logger::log().debug() << "median alpha: " << median_alpha_ << endl;
//    TLorentzVector median_alpha_thad(median_alpha_*perm.THad().Px(), median_alpha_*perm.THad().Py(), median_alpha_*perm.THad().Pz(), median_alpha_*perm.THad().E());
//    return median_alpha_thad;
//}

