#include "Analyses/URTTbar/interface/Alpha_Corrections.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"

Alpha_Corrections::Alpha_Corrections(){

    URParser &parser = URParser::instance();
    //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
    parser.addCfgParameter<string>("alpha_correction", "fit", "");
    parser.parseArguments();

        //get filename that has fit parameters
    string fname = parser.getCfgPar<string>("alpha_correction", "fit");

        // open file and create dict from json
    std::ifstream ifs(DataFile(fname).path().c_str());
    ifs >> jsonFile_;
}

Alpha_Corrections::~Alpha_Corrections()
{}

double Alpha_Corrections::alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range)
{
    if( fit_var == "None" ){
        //Logger::log().debug() << "Alpha correction set to 1." << endl;
        alpha_ = 1.;
    }
    else{
        if( fit_var != "E" && fit_var != "P" && fit_var != "None" ){
            Logger::log().error() << "Only energy (E) and momentum (P) variables are supported!" << endl;
            throw 42;
        }
    
        string mtt_range;
    
        if( fit_range == "Mtt" ){
            if( 200. <= perm.LVect().M() && perm.LVect().M() < 350. ) mtt_range = "Mttbar200to350";
            else if( 350. <= perm.LVect().M() && perm.LVect().M() < 400. ) mtt_range = "Mttbar350to400";
            else if( 400. <= perm.LVect().M() && perm.LVect().M() < 500. ) mtt_range = "Mttbar400to500";
            else if( 500. <= perm.LVect().M() && perm.LVect().M() < 700. ) mtt_range = "Mttbar500to700";
            else if( 700. <= perm.LVect().M() && perm.LVect().M() < 1000. ) mtt_range = "Mttbar700to1000";
            else if( perm.LVect().M() >= 1000. ) mtt_range = "Mttbar1000toInf";
        }
        else if( fit_range == "All" ) mtt_range = fit_range;
        else{
            Logger::log().error() << "Fit range not supported!" << endl;
            throw 42;
        }
    
        vector<double> pars = jsonFile_["THad_"+fit_var][mtt_range][fit_degree+"_g0.9"]["Pars"];
    
        if( fit_degree == "1D" ) alpha_ = pars[0]*( 173.1/perm.THad().M() ) + pars[1];
        else if ( fit_degree == "2D" ) alpha_ = pars[0]*( pow(173.1/perm.THad().M(),2) ) + pars[1]*( 173.1/perm.THad().M() ) + pars[2];
        else{
            Logger::log().error() << "Only linear (1D) and quadratic (2D) fits are supported!" << endl;
            throw 42;
        }
    }

    return alpha_;

}

TLorentzVector Alpha_Corrections::Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range) // scaled reco thad 4-vec
{
    alpha_ = alpha(perm, fit_degree, fit_var, fit_range);
    TLorentzVector alpha_thad(alpha_*perm.THad().Px(), alpha_*perm.THad().Py(), alpha_*perm.THad().Pz(), alpha_*perm.THad().E());
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

double Alpha_Corrections::median_alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range)
{
    if( fit_var == "None" ){
        //Logger::log().debug() << "Alpha correction set to 1." << endl;
        alpha_ = 1.;
    }
    else{
        if( fit_var != "E" && fit_var != "P" && fit_var != "None" ){
            Logger::log().error() << "Only energy (E) and momentum (P) variables are supported!" << endl;
            throw 42;
        }
    
        string mtt_range;
        double MTHad;
    
        if( fit_range == "Mtt" ){
            if( 200. <= perm.LVect().M() && perm.LVect().M() < 350. ){
                mtt_range = "Mttbar200to350";
                    //values taken from http://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2016Data/JetpTcut20/LeadPt0/Only_Alpha_Correction/Full/ttJets/3J_Event_Plots/Lost_BP/Alpha_Correction/THad_E/y_proj/Alpha_THad_E_Mttbar200to350_yprojections.png
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.05;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.19;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.35;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.51;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.69;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.89;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 2.07;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 2.29;
            }
            else if( 350. <= perm.LVect().M() && perm.LVect().M() < 400. ){
                mtt_range = "Mttbar350to400";
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.03;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.15;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.27;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.39;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.49;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.59;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.67;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.73;
            }
            else if( 400. <= perm.LVect().M() && perm.LVect().M() < 500. ){
                mtt_range = "Mttbar400to500";
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.01;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.11;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.21;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.31;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.39;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.43;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.47;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.51;
            }
            else if( 500. <= perm.LVect().M() && perm.LVect().M() < 700. ){
                mtt_range = "Mttbar500to700";
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.07;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.17;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.23;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.29;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.33;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.35;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.37;
            }
            else if( 700. <= perm.LVect().M() && perm.LVect().M() < 1000. ){
                mtt_range = "Mttbar700to1000";
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.05;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.13;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.17;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.23;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.23;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.27;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.29;
            }
            else if( perm.LVect().M() >= 1000. ){
                mtt_range = "Mttbar1000toInf";
                if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 0.99;
                if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.03;
                if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.09;
                if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.15;
                if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.19;
                if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.23;
                if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.21;
                if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.23;
            }
        }
        else if( fit_range == "All" ){
            mtt_range = fit_range;
                // values taken from http://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2016Data/JetpTcut20/LeadPt0/Only_Alpha_Correction/Full/ttJets/3J_Event_Plots/Lost_BP/Alpha_Correction/THad_E/y_proj/Alpha_THad_E_yprojections.png
            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) MTHad = 1.01;
            if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) MTHad = 1.11;
            if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) MTHad = 1.23;
            if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) MTHad = 1.35;
            if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) MTHad = 1.47;
            if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) MTHad = 1.57;
            if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) MTHad = 1.65;
            if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) MTHad = 1.71;
        }
        else{
            Logger::log().error() << "Fit range not supported!" << endl;
            throw 42;
        }
    
        vector<double> pars = jsonFile_["THad_"+fit_var][mtt_range][fit_degree+"_g0.9"]["Pars"];
   
        if( fit_degree == "1D" ) median_alpha_ = pars[0]*( MTHad ) + pars[1];
        else if ( fit_degree == "2D" ) median_alpha_ = pars[0]*( pow(MTHad,2) ) + pars[1]*( MTHad ) + pars[2];
        else{
            Logger::log().error() << "Only linear (1D) and quadratic (2D) fits are supported!" << endl;
            throw 42;
        }
    }

    return median_alpha_;

}

TLorentzVector Alpha_Corrections::median_Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range) // scaled reco thad 4-vec
{
    median_alpha_ = median_alpha(perm, fit_degree, fit_var, fit_range);
    Logger::log().debug() << "median alpha: " << median_alpha_ << endl;
    TLorentzVector median_alpha_thad(median_alpha_*perm.THad().Px(), median_alpha_*perm.THad().Py(), median_alpha_*perm.THad().Pz(), median_alpha_*perm.THad().E());
    return median_alpha_thad;
}

