#include "Analyses/URTTbar/interface/QCD_WeightProducer.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "URAnalysis/AnalysisFW/interface/Wards.h"


QCD_WeightProducer::QCD_WeightProducer(){

    URParser &parser = URParser::instance();
    parser.addCfgParameter<string>("qcd_correction", "fit", "");
    parser.parseArguments();

        //get filename that has fit parameters
    string fname = parser.getCfgPar<string>("qcd_correction", "fit");
    TFile qcd_file(DataFile(fname).path().c_str());

    HistoOwnershipWard hward;
    DirectoryWard dward;

    frac_b_CSV_p0_ = get_from<TH1D>(qcd_file, "Fractions_b_CSV_Pt/p0", "frac_b_CSV_pt_p0_clone");
    frac_b_nonCSV_p0_ = get_from<TH1D>(qcd_file, "Fractions_b_nonCSV_Pt/p0", "frac_b_nonCSV_pt_p0_clone");
    eff_b_p0_ = get_from<TH1D>(qcd_file, "Efficiencies_b_Pt_fit/p0", "eff_b_pt_p0_clone");
    eff_b_p1_ = get_from<TH1D>(qcd_file, "Efficiencies_b_Pt_fit/p1", "eff_b_pt_p1_clone");
    eff_b_p2_ = get_from<TH1D>(qcd_file, "Efficiencies_b_Pt_fit/p2", "eff_b_pt_p2_clone");
    eff_p_p0_ = get_from<TH1D>(qcd_file, "Efficiencies_p_Pt_fit/p0", "eff_p_pt_p0_clone");
    eff_p_p1_ = get_from<TH1D>(qcd_file, "Efficiencies_p_Pt_fit/p1", "eff_p_pt_p1_clone");
    eff_p_p2_ = get_from<TH1D>(qcd_file, "Efficiencies_p_Pt_fit/p2", "eff_p_pt_p2_clone");
    r_data_nonCSV_p0_ = get_from<TH1D>(qcd_file, "r_data_nonCSV_Pt_fit/p0", "r_data_nonCSV_pt_p0_clone");
    r_data_nonCSV_p1_ = get_from<TH1D>(qcd_file, "r_data_nonCSV_Pt_fit/p1", "r_data_nonCSV_pt_p1_clone");
    r_data_nonCSV_p2_ = get_from<TH1D>(qcd_file, "r_data_nonCSV_Pt_fit/p2", "r_data_nonCSV_pt_p2_clone");

}

QCD_WeightProducer::~QCD_WeightProducer()
{}

double QCD_WeightProducer::qcd_weight( double lep_kvar, string fit_variable, string fit_variation)
//double QCD_WeightProducer::qcd_weight(Permutation &perm, string fit_variable, string fit_variation)
{

    //if( fit_variable != "Pt" && fit_variable != "Eta" ){
    //    Logger::log().error() << "Only lep pT and eta are supported for qcd weight estimation!" << endl;
    //    throw 42;
    //}
    if( fit_variable != "Pt" ){
        Logger::log().error() << "Only lep pT is supported for qcd weight estimation!" << endl;
        throw 42;
    }
    if( fit_variation != "nominal" && fit_variation != "up" && fit_variation != "down" ){
        Logger::log().error() << "Nominal, up, or down variation must be chosen!" << endl;
        throw 42;
    }

    //double lep_kvar = -1.;
    //if( fit_variable == "Pt" ) lep_kvar = perm.L()->Pt();
    //else lep_kvar = perm.L()->Eta();

    double qcd_weight = 1.; // reset weight

    double r_data =  SIGMOID( get_par( r_data_nonCSV_p0_, fit_variation ), get_par( r_data_nonCSV_p1_, fit_variation ),
                             get_par( r_data_nonCSV_p2_, fit_variation ), lep_kvar );
    double eff_b = SIGMOID( get_par( eff_b_p0_, fit_variation ), get_par( eff_b_p1_, fit_variation ),
                            get_par( eff_b_p2_, fit_variation ), lep_kvar ); 

    double eff_p = SIGMOID( get_par( eff_p_p0_, fit_variation ), get_par( eff_p_p1_, fit_variation ),
                            get_par( eff_p_p2_, fit_variation ), lep_kvar ); 

    double frac_b_CSV = get_par( frac_b_CSV_p0_, fit_variation );
    double frac_b_nonCSV = get_par( frac_b_nonCSV_p0_, fit_variation );

    double r_QCD = frac_b_nonCSV*eff_b + (1. - frac_b_nonCSV)*eff_p;
    double sf = r_data/r_QCD;

    qcd_weight = ( frac_b_CSV*sf*eff_b )/( 1. - sf*eff_b ) + ( ( 1. - frac_b_CSV )*sf*eff_p )/( 1. - sf*eff_p );

    return qcd_weight;

}

double QCD_WeightProducer::get_par( std::shared_ptr<TH1> h, string variation ) const
{
    double binval = h->GetBinContent( 1 );
    double binerr = h->GetBinError( 1 );

    if( variation == "nominal" ) return binval;
    else if( variation == "up" ) return binval+binerr;
    else if( variation == "down" ) return binval-binerr;
    else {
        Logger::log().error() << "Nominal, up, or down variation must be chosen!" << endl;
        throw 42;
    }
}


double QCD_WeightProducer::SIGMOID( double p0, double p1, double p2, double X ) const
{

    double sigval =  p0/(1. +TMath::Exp( -p1*( X - p2 ) ) );

    return sigval;
}

