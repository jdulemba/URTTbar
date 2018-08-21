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
    if( fit_var != "E" && fit_var != "P" ){
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
    else if ( fit_degree == "2D" ) alpha_ = pars[0] + pars[1]*( 173.1/perm.THad().M() ) + pars[2]*( pow(173.1/perm.THad().M(),2) );
    else{
        Logger::log().error() << "Only linear (1D) and quadratic (2D) fits are supported!" << endl;
        throw 42;
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

