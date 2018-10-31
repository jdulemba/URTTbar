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
            Logger::log().error() << "Only linear (1d) and quadratic (2d) fits are supported!" << endl;
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

