#include "TROOT.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "Analyses/URTTbar/interface/DR_TTGenMatcher.h"
#include "Analyses/URTTbar/interface/Alpha_Corrections.h"
#include "Analyses/URTTbar/interface/QCD_WeightProducer.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/GenObject.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include <algorithm>
#include <list>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include <boost/algorithm/string/predicate.hpp>
#include "TMath.h"
#include "TRandom3.h"
#include <set>
#include "Analyses/URTTbar/interface/IDMet.h"
#include "TUUID.h"
#include "Analyses/URTTbar/interface/systematics.h"
//#include "Analyses/URTTbar/interface/TTBarPlots.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"

#include "TTree.h"
//#include <string>

using namespace TMath;
using namespace std;

class ttbar_post_alpha_reco : public AnalyzerBase
{

    public:
        //enum TTNaming_Merged {RIGHT, MERGE_SWAP, MERGE, WRONG};

    private:
        //counters
        unsigned long evt_idx_ = 0; //event index

        // bins for hists
        int mass_bins_ = 180;
        int pt_bins_ = 100;
        int eta_bins_ = 100;
        int costh_bins_ = 200;

        double pt_min_ = 0.;
        double pt_max_ = 1000.;
        double eta_min_ = -2.5;
        double eta_max_ = 2.5;
        double costh_min_ = -1.;
        double costh_max_ = 1.;

        // mass disc
        int massdisc_bins_ = 120;
        double massdisc_min_ = -10.;
        double massdisc_max_ = 20.;

        // ns disc
        int nsdisc_bins_ = 80;
        double nsdisc_min_ = -5.;
        double nsdisc_max_ = 15.;

        // total disc
        int combdisc_bins_ = 160;
        double combdisc_min_ = -10.;
        double combdisc_max_ = 30.;

        //histograms
        map< string, map < string, map< string, map< string, map< string, RObject > > > > > histos_merged_;
        map< string, map < string, map< string, map< string, map< string, map< string, RObject > > > > > > histos_disc_;
        //unordered_map<string, map< string, RObject> > histos_;
        //unordered_map<string, map< string, RObject> > histos_;
        map<string, map< string, RObject> > histos_;

        CutFlowTracker tracker_; //tracks how many events pass cuts

        //switches
        bool isData_, isTTbar_, isSignal_;

        //selectors and helpers
        TTGenParticleSelector genp_selector_; //selects generator level objects
        TTGenMatcher matcher_; //matches particles on generator level
        DR_TTGenMatcher dr_matcher_;
        Alpha_Corrections alpha_corr_;
        QCD_WeightProducer qcd_reweight_;
        TTPermutator permutator_;
        TTObjectSelector object_selector_; //selects ttbar objects
        float evt_weight_;
        TRandom3 randomizer_;// = TRandom3(98765);
        MCWeightProducer mc_weights_;
        BTagSFProducer btag_sf_;
        TTBarSolver solver_; //solves ttbar events

        // Lepton scale factors
        LeptonSF electron_sf_, muon_sf_;

        vector<systematics::SysShifts> systematics_;

        // btag cuts
        IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;
        //IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
        //IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
        //IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;

        float cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_;

    public:
        ttbar_post_alpha_reco(const std::string output_filename):
            AnalyzerBase("ttbar_post_alpha_reco", output_filename),
            tracker_(),
            //genp_selector_(TTGenParticleSelector::SelMode::LHE),
            //genp_selector_(),
            matcher_(),
            permutator_(),
            object_selector_(),
            mc_weights_(),
            evt_weight_(1.),
            electron_sf_("electron_sf", false),
            muon_sf_("muon_sf"),
            randomizer_(),
            btag_sf_("permutations.tightb", "permutations.looseb"),
            solver_(true),
            //pdf_uncs_(254),
            systematics_()
            //sync_(false),
            //sync_tree_(0),
            //sync_info_(){

            //genp_selector_(TTGenParticleSelector::SelMode::LHE), //for ttJets files
    {
        Logger::log().debug() << "ttbar_alpha_reco ctor" << endl;
        cut_tight_b_ = btag_sf_.tight_cut();
        cut_loose_b_ = btag_sf_.loose_cut();

        // get parameters from cfg file
        URParser& parser = URParser::instance();

        parser.addCfgParameter<float>("jets", "ptmin", "minimum pt");
        parser.addCfgParameter<float>("jets", "etamax", "maximum eta");
        parser.addCfgParameter<float>("jets", "lead_ptmin", "minimum leading jet pt");
        parser.parseArguments();

        cut_jet_ptmin_ = parser.getCfgPar<float>("jets", "ptmin" );
        cut_jet_etamax_ = parser.getCfgPar<float>("jets", "etamax");
        cut_leadjet_ptmin_ = parser.getCfgPar<float>("jets", "lead_ptmin" );


        //find out which sample we're running on
        opts::variables_map &values = URParser::instance().values();
        string output_file = values["output"].as<std::string>();
        string sample = systematics::get_sample(output_file);
        isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
        isTTbar_ = boost::starts_with(sample, "ttJets");
        isData_  = boost::starts_with(sample, "data");

        if( !(isTTbar_ || isSignal_) ) {
            Logger::log().error() << "This analyzer is only supposed to run on ttbar and signal samples!" << endl;
            throw 49;
        }

        if( isTTbar_ ) genp_selector_ = TTGenParticleSelector(TTGenParticleSelector::SelMode::LHE);
        else genp_selector_ = TTGenParticleSelector();

        //set tracker
        tracker_.use_weight(&evt_weight_);
        object_selector_.set_tracker(&tracker_);

        Logger::log().debug() << "NOSYS: " << values.count("nosys") << endl;
        //choose systematics to run based on sample
        systematics_ = {systematics::SysShifts::NOSYS};
        if( !(isTTbar_ || isSignal_) ) {
            Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
            throw 49;
        }

        //Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << ", isSignal_: " << isSignal_ << endl;
        if(isData_) {
            if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
            else object_selector_.lepton_type(1);
        }

        if(!isData_) mc_weights_.init(sample);

    };

        TDirectory* getDir(string path){
            TDirectory* out = (TDirectory*) outFile_.Get(path.c_str());
            if(out) return out;
            outFile_.mkdir(path.c_str());
            return (TDirectory*) outFile_.Get(path.c_str());
        }

        template <class H, typename ... Args>
            void book(string folder, string name, Args ... args)
            {
                getDir(folder)->cd();
                histos_[folder][name] = RObject::book<H>(name.c_str(), args ...);
            }

        //This method is called once per job at the beginning of the analysis
        //book here your histograms/tree and run every initialization needed

        /// book and fill plots after hadronic top mass corrections applied to lost-jet events
        void book_reconstruction_plots( string folder, string cat ){

            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};

                //book uncorrected plots
            book<TH1F>(folder+"/Reconstruction/Mass"+cat+"/Uncorrected", "THad", "", 500, 0., 500.);
            book<TH1F>(folder+"/Reconstruction/Mass"+cat+"/Uncorrected", "TTbar", "", mass_bins_, 200., 2000.);
            book<TH2D>(folder+"/Reconstruction/Mass"+cat+"/Uncorrected", "Reco_vs_Gen_TTbar", ";Gen M(ttbar); Reco M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);

            //book<TH1F>(folder+"/Reconstruction/Costh"+cat, "THad_Labframe", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh"+cat+"/Uncorrected", "THad", "", costh_bins_, -1., 1.);
            book<TH2D>(folder+"/Reconstruction/Costh"+cat+"/Uncorrected", "Reco_vs_Gen_THad", ";Gen Cos(theta)(t_{h}); Reco Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);

                //book corrected plots
            for( auto fit_degree : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = "/Corrected_"+fit_degree+"_"+fit_var+"_"+fit_range;

                        book<TH1F>(folder+"/Reconstruction/Mass"+cat+hname, "THad", "", 500, 0., 500.);
                        book<TH1F>(folder+"/Reconstruction/Mass"+cat+hname, "TTbar", "", mass_bins_, 200., 2000.);
                        book<TH2D>(folder+"/Reconstruction/Mass"+cat+hname, "Reco_vs_Gen_TTbar", ";Gen M(ttbar); Reco M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);

                        //book<TH1F>(folder+"/Reconstruction/Costh"+cat, "THad_Labframe", "", costh_bins_, -1., 1.);
                        book<TH1F>(folder+"/Reconstruction/Costh"+cat+hname, "THad", "", costh_bins_, -1., 1.);
                        book<TH2D>(folder+"/Reconstruction/Costh"+cat+hname, "Reco_vs_Gen_THad", ";Gen Cos(theta)(t_{h}); Reco Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);
                    }
                }
            }
        }

        void book_resolution_plots( string folder, string cat ){

            vector<string> shoulder_parts = {"L30", "G30"};
            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};

                //book uncorreted plots
            book<TH1F>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "TLep", "", 500, -1000., 500.);
            book<TH1F>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "THad", "", 500, -1000., 500.);
            book<TH1F>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "TTbar", "", 500, -1000., 1000.);
            //book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Alpha_vs_Gen_MTTbar", "", 80, 0., 2000., 100, 0., 5.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_MTHad_vs_Gen_MTTbar", "", 80, 0., 2000., 500, -1000., 500.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_MTHad_vs_Gen_THadPt", "", 100, 0., 1000., 500, -1000., 500.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_MTHad_vs_Uncorrected_MTHad", ";173.1/Uncorrected M(t_{h}); M(t_{h}) Resolution", 8, 0.9, 2.5, 500, -1000., 500.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_THadPt_vs_Gen_THadPt", "", 100, 0., 1000., 200., -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_THadEta_vs_Gen_THadEta", "", 100, -3., 3., 200., -6., 6.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_EtaTTbar_vs_Gen_EtaTTbar", "", 100, -3., 3., 200., -6., 6.);
            //book<TH1F>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Frac_THad", "", mass_bins_, -2., 2.);
            //book<TH1F>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Frac_TTbar", "", mass_bins_, -2., 2.);
            //book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Frac_TTbar_vs_Gen_THadPt", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_TLepPt_vs_Gen_TLepPt", "", 100, 0., 1000., 200., -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass"+cat+"/Uncorrected", "Reso_TLepEta_vs_Gen_TLepEta", "", 100, -3., 3., 200., -6., 6.);

            book<TH1F>(folder+"/Resolution/Costh"+cat+"/Uncorrected", "THad", "", costh_bins_, -2., 2.);

                //book corrected plots
            for( auto fit_degree : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = "/Corrected_"+fit_degree+"_"+fit_var+"_"+fit_range;

                        book<TH1F>(folder+"/Resolution/Mass"+cat+hname, "THad", "", 500, -1000., 500.);
                        book<TH1F>(folder+"/Resolution/Mass"+cat+hname, "TTbar", "", 500, -1000., 1000.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Alpha_vs_Gen_MTTbar", "", 80, 0., 2000., 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Median_Alpha_vs_Gen_MTTbar", "", 80, 0., 2000., 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Alpha_Diff_vs_Gen_MTTbar", ";;#alpha computed-median", 80, 0., 2000., 200, -1., 1.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Computed_Alpha_vs_Uncorrected_MTHad", "", 8, 0.9, 2.5, 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Median_Alpha_vs_Uncorrected_MTHad", "", 8, 0.9, 2.5, 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Alpha_Diff_vs_Uncorrected_MTHad", ";;#alpha computed-median", 8, 0.9, 2.5, 200, -1., 1.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_MTHad_vs_Gen_MTTbar", "", 80, 0., 2000., 500, -1000., 500.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_MTHad_vs_Gen_THadPt", "", 100, 0., 1000., 500, -1000., 500.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_MTHad_vs_Uncorrected_MTHad", ";173.1/Uncorrected M(t_{h}); M(t_{h}) Resolution", 8, 0.9, 2.5, 500, -1000., 500.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_THadPt_vs_Gen_THadPt", "", 100, 0., 1000., 200., -1000., 1000.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_THadEta_vs_Gen_THadEta", "", 100, -3., 3., 200., -6., 6.);
                        book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Reso_EtaTTbar_vs_Gen_EtaTTbar", "", 100, -3., 3., 200., -6., 6.);
                        //book<TH1F>(folder+"/Resolution/Mass"+cat+hname, "Frac_THad"+hname, "", mass_bins_, -2., 2.);
                        //book<TH1F>(folder+"/Resolution/Mass"+cat+hname, "Frac_TTbar"+hname, "", mass_bins_, -2., 2.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname, "Frac_TTbar_vs_Gen_THadPt"+hname, "", 100, 0., 1000., mass_bins_, -2., 2.);

                        book<TH1F>(folder+"/Resolution/Costh"+cat+hname, "THad", "", costh_bins_, -2., 2.);

                        for( auto shoulder_part : shoulder_parts ){
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_Mass_Uncorrected", "", 500, 0., 500.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_Mass_Corrected", "", 500, 0., 500.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_Pt_Uncorrected", "", 1000, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_Pt_Corrected", "", 1000, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_P_Uncorrected", "", 500, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_P_Corrected", "", 500, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_E_Uncorrected", "", 500, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "THad_E_Corrected", "", 500, 0., 1000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "TTbar_Mass_Uncorrected", "", mass_bins_, 0., 2000.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "TTbar_Mass_Corrected", "", mass_bins_, 0., 2000.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Alpha_vs_Gen_MTTbar", "", 80, 0., 2000., 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Median_Alpha_vs_Gen_MTTbar", "", 80, 0., 2000., 300, 0., 3.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Alpha_Diff_vs_Gen_MTTbar", ";;#alpha computed-median", 80, 0., 2000., 200, -1., 1.);
                            book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Computed_Alpha_vs_Uncorrected_MTHad", "", 8, 0.9, 2.5, 300, 0., 3.);
                            //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Median_Alpha_vs_Uncorrected_MTHad", "", 8, 0.9, 2.5, 300, 0., 3.);
                            //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Alpha_Diff_vs_Uncorrected_MTHad", ";;#alpha computed-median", 8, 0.9, 2.5, 200, -1., 1.);
                            book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "AlphaE_vs_Uncorrected_MTHad", "; 173.1/Uncorrected M(t_{h});#alpha_{E}", 8, 0.9, 2.5, 100, 0., 5.);
                            book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "AlphaP_vs_Uncorrected_MTHad", "; 173.1/Uncorrected M(t_{h});#alpha_{P}", 8, 0.9, 2.5, 100, 0., 5.);

                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_over_Reco_BHadE", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_over_Reco_BHadPt", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Match_over_Reco_BHadE", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Match_over_Reco_BHadPt", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_WJa_over_Reco_WJetE", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_WJa_over_Reco_WJetPt", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_WJb_over_Reco_WJetE", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Gen_WJb_over_Reco_WJetPt", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Match_over_Reco_WJetE", "", 50, 0., 5.);
                            book<TH1F>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Match_over_Reco_WJetPt", "", 50, 0., 5.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_MTHad_vs_Gen_MTTbar", "", 80, 0., 2000., 500, -1000., 500.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_MTHad_vs_Gen_THadPt", "", 100, 0., 1000., 500, -1000., 500.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_MTHad_vs_Uncorrected_MTHad", ";173.1/Uncorrected M(t_{h}); M(t_{h}) Resolution", 8, 0.9, 2.5, 500, -1000., 500.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_THadPt_vs_Gen_THadPt", "", 100, 0., 1000., 200., -1000., 1000.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_THadEta_vs_Gen_THadEta", "", 100, -3., 3., 200., -6., 6.);
                        //book<TH2D>(folder+"/Resolution/Mass"+cat+hname+"/"+shoulder_part, "Reso_EtaTTbar_vs_Gen_EtaTTbar", "", 100, -3., 3., 200., -6., 6.);

                        }
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_LCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_GCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_LCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_GCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_LCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_GCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_LCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                        //book<TH2D>(folder+"/Resolution"+cat+hname, "Reso_MTTbar_GCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                    }
                }
            }

        }

        void book_parton_accept_plots( string folder, string cat ){

            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};
            vector<string> npartons = { "All", "2Partons", "3Partons", "4Partons" };

                //book uncorrected plots
            for( auto nparts : npartons ){
                book<TH1F>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTTbar_"+nparts, "", mass_bins_, -1000., 1000.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTTbar_vs_Gen_MTTbar_"+nparts, "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTHad_vs_Gen_MTTbar_"+nparts, "", 80, 0., 2000., 500, -1000., 500.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTHad_vs_Uncorrected_MTHad_"+nparts, ";173.1/Uncorrected M(t_{h}); M(t_{h}) Resolution", 8, 0.9, 2.5, 500, -1000., 500.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_THadPt_vs_Gen_THadPt_"+nparts, "", 100, 0., 1000., 200., -1000., 1000.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_THadEta_vs_Gen_THadEta_"+nparts, "", 100, -3., 3., 200., -6., 6.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_EtaTTbar_vs_Gen_EtaTTbar_"+nparts, "", 100, -3., 3., 200., -6., 6.);
                book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Frac_MTTbar_vs_Gen_THadPt_"+nparts, "", 100, 0., 1000., mass_bins_, -2., 2.);
            }

                //book corrected plots
            for( auto fit_degree : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = "_Corrected_"+fit_degree+"_"+fit_var+"_"+fit_range;

                        for( auto nparts : npartons ){
                            book<TH1F>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTTbar_"+nparts+hname, "", mass_bins_, -1000., 1000.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTTbar_vs_Gen_MTTbar_"+nparts+hname, "", 80, 0., 2000., mass_bins_, -1000., 1000.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTHad_vs_Gen_MTTbar_"+nparts+hname, "", 80, 0., 2000., 500, -1000., 500.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_MTHad_vs_Uncorrected_MTHad_"+nparts+hname, ";173.1/Uncorrected M(t_{h}); M(t_{h}) Resolution", 8, 0.9, 2.5, 500, -1000., 500.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_THadPt_vs_Gen_THadPt_"+nparts+hname, "", 100, 0., 1000., 200., -1000., 1000.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_THadEta_vs_Gen_THadEta_"+nparts+hname, "", 100, -3., 3., 200., -6., 6.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Reso_EtaTTbar_vs_Gen_EtaTTbar_"+nparts+hname, "", 100, -3., 3., 200., -6., 6.);
                            book<TH2D>(folder+"/Resolution/Parton_Acceptance"+cat, "Frac_MTTbar_vs_Gen_THadPt_"+nparts+hname, "", 100, 0., 1000., mass_bins_, -2., 2.);
                        }
                    }
                }
            }

        }

        void book_plots( string folder, string cat ){
            book_reconstruction_plots(folder, cat);
            book_resolution_plots(folder, cat);
            //book_parton_accept_plots(folder, cat);
        }


        void fill_reconstruction_plots( string folder, string cat, GenTTBar &ttbar, Permutation &perm ){
            auto uncorr_reco_mass_dir = histos_.find(folder+"/Reconstruction/Mass"+cat+"/Uncorrected");
            auto uncorr_reco_costh_dir = histos_.find(folder+"/Reconstruction/Costh"+cat+"/Uncorrected");


            // costheta variables
            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >

            // costhetastar for original perm thad
            double reco_thad_cth = reco_thad_cthstar(perm);

            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};
            for( auto fit_deg : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = "/Corrected_"+fit_deg+"_"+fit_var+"_"+fit_range;
                        auto corr_reco_mass_dir = histos_.find(folder+"/Reconstruction/Mass"+cat+hname);
                        auto corr_reco_costh_dir = histos_.find(folder+"/Reconstruction/Costh"+cat+hname);

                            // get corrected thad 4-vec and costhetastar vals
                        TLorentzVector Alpha_THad = alpha_corr_.Alpha_THad(perm, fit_deg, fit_var, fit_range);
                        double alpha_thad_cth = alpha_corr_.alpha_thad_cthstar(perm, fit_deg, fit_var, fit_range);

                        //double nominal_qcd_weight = qcd_reweight_.qcd_weight( perm, "Pt", "nominal" );
                        //double up_qcd_weight = qcd_reweight_.qcd_weight( perm, "Pt", "up" );
                        //double down_qcd_weight = qcd_reweight_.qcd_weight( perm, "Pt", "down" );
                        //cout << "Qcd nominal: " << nominal_qcd_weight << ", up: " << up_qcd_weight << ", down: " << down_qcd_weight << endl << endl;

                            // fill plots for corrected thad
                        corr_reco_mass_dir->second["THad"].fill( Alpha_THad.M(), evt_weight_ );
                        corr_reco_mass_dir->second["TTbar"].fill( (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        corr_reco_mass_dir->second["Reco_vs_Gen_TTbar"].fill( ttbar.M(), (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        corr_reco_costh_dir->second["THad"].fill( alpha_thad_cth, evt_weight_ );
                        corr_reco_costh_dir->second["Reco_vs_Gen_THad"].fill( gen_cosths.first, alpha_thad_cth, evt_weight_ );
                    }
                }
            }

              // fill plots for original thad
            uncorr_reco_mass_dir->second["THad"].fill( perm.THad().M(), evt_weight_ );
            uncorr_reco_mass_dir->second["TTbar"].fill( perm.LVect().M(), evt_weight_ );
            uncorr_reco_mass_dir->second["Reco_vs_Gen_TTbar"].fill( ttbar.M(), perm.LVect().M(), evt_weight_ );
            uncorr_reco_costh_dir->second["THad"].fill( reco_thad_cth, evt_weight_ );
            uncorr_reco_costh_dir->second["Reco_vs_Gen_THad"].fill(gen_cosths.first, reco_thad_cth, evt_weight_);
        }


        void fill_shoulder_resolution_plots( string folder, GenTTBar &ttbar, Permutation &perm, TLorentzVector &corr_thad, string fit_deg, string fit_var, string fit_range){

            auto dir = histos_.find(folder);

            dir->second["THad_Mass_Uncorrected"].fill( perm.THad().M(), evt_weight_ );
            dir->second["THad_Mass_Corrected"].fill( corr_thad.M(), evt_weight_ );
            dir->second["THad_Pt_Uncorrected"].fill( perm.THad().Pt(), evt_weight_ );
            dir->second["THad_Pt_Corrected"].fill( corr_thad.Pt(), evt_weight_ );
            dir->second["THad_P_Uncorrected"].fill( perm.THad().P(), evt_weight_ );
            dir->second["THad_P_Corrected"].fill( corr_thad.P(), evt_weight_ );
            dir->second["THad_E_Uncorrected"].fill( perm.THad().E(), evt_weight_ );
            dir->second["THad_E_Corrected"].fill( corr_thad.E(), evt_weight_ );

            dir->second["TTbar_Mass_Uncorrected"].fill( perm.LVect().M(), evt_weight_ );
            dir->second["TTbar_Mass_Corrected"].fill( (corr_thad + perm.TLep()).M(), evt_weight_ );

            dir->second["Computed_Alpha_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
            //dir->second["Median_Alpha_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
            //dir->second["Alpha_Diff_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range)-alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
            dir->second["AlphaE_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), ttbar.had_top()->E()/perm.THad().E(), evt_weight_ );
            dir->second["AlphaP_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), ttbar.had_top()->P()/perm.THad().P(), evt_weight_ );

            dir->second["Gen_over_Reco_BHadE"].fill( ttbar.had_b()->E()/perm.BHad()->E(), evt_weight_ );
            dir->second["Gen_over_Reco_BHadPt"].fill( ttbar.had_b()->Pt()/perm.BHad()->Pt(), evt_weight_ );

            if( perm.BHad()->match() ){
                dir->second["Match_over_Reco_BHadE"].fill( perm.BHad()->match()->E()/perm.BHad()->E(), evt_weight_ );
                dir->second["Match_over_Reco_BHadPt"].fill( perm.BHad()->match()->Pt()/perm.BHad()->Pt(), evt_weight_ );

                //cout << "Perm bhad E: " << perm.BHad()->E() << ", bhad match E: " << perm.BHad()->match()->E() << endl;
                //cout << "gen bhad E: " << ttbar.had_b()->E() << ", gen blep E: " << ttbar.lep_b()->E() << endl;
                //cout << "Perm bhad Pt: " << perm.BHad()->Pt() << ", bhad match Pt: " << perm.BHad()->match()->Pt() << endl;
                //cout << "gen bhad Pt: " << ttbar.had_b()->Pt() << ", gen blep Pt: " << ttbar.lep_b()->Pt() << endl << endl;
            }

            dir->second["Gen_WJa_over_Reco_WJetE"].fill( ttbar.had_W()->first->E()/perm.WJa()->E(), evt_weight_ );
            dir->second["Gen_WJa_over_Reco_WJetPt"].fill( ttbar.had_W()->first->Pt()/perm.WJa()->Pt(), evt_weight_ );
            dir->second["Gen_WJb_over_Reco_WJetE"].fill( ttbar.had_W()->second->E()/perm.WJa()->E(), evt_weight_ );
            dir->second["Gen_WJb_over_Reco_WJetPt"].fill( ttbar.had_W()->second->Pt()/perm.WJa()->Pt(), evt_weight_ );

            if( perm.WJa()->match() ){
                dir->second["Match_over_Reco_WJetE"].fill( perm.WJa()->match()->E()/perm.WJa()->E(), evt_weight_ );
                dir->second["Match_over_Reco_WJetPt"].fill( perm.WJa()->match()->Pt()/perm.WJa()->Pt(), evt_weight_ );
                //cout << "Perm wjet E: " << perm.WJa()->E() << ", wjet match: " << perm.WJa()->match()->E() << endl;
                //cout << "gen wja E: " << ttbar.had_W()->first->E() << ", gen wjb E: " << ttbar.had_W()->second->E() << endl;
                //cout << "Perm wjet Pt: " << perm.WJa()->Pt() << ", wjet match: " << perm.WJa()->match()->Pt() << endl;
                //cout << "gen wja Pt: " << ttbar.had_W()->first->Pt() << ", gen wjb Pt: " << ttbar.had_W()->second->Pt() << endl << endl;
            }
            //else cout << folder << endl;
            //cout << "gen wja E: " << ttbar.had_W()->first->E() << ", gen wjb E: " << ttbar.had_W()->second->E() << endl;

        }



        void fill_resolution_plots( string folder, string cat, GenTTBar &ttbar, Permutation &perm ){
            //auto reso_dir = histos_.find(folder+"/Resolution"+cat);
            auto uncorr_reso_costh_dir = histos_.find(folder+"/Resolution/Costh"+cat+"/Uncorrected");
            auto uncorr_reso_mass_dir = histos_.find(folder+"/Resolution/Mass"+cat+"/Uncorrected");

//cout << endl;
//cout << "Event: " << evt_idx_ << endl;
//cout << "   Original THad M(thad) Res: " << ttbar.had_top()->M() - perm.THad().M() << ", M(tt) Res: " << ttbar.M() - perm.LVect().M() << ", 173.1/Reco M(thad): " << 173.1/perm.THad().M() << ", Reco M(tt): " << perm.LVect().M() << ", pT: " << perm.THad().Pt() << ", E: " << perm.THad().E() << ", P: " << perm.THad().P() << endl;

            // costheta variables
            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            // costhetastar for original perm thad
            double reco_thad_cth = reco_thad_cthstar(perm);

            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};
            for( auto fit_deg : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = "/Corrected_"+fit_deg+"_"+fit_var+"_"+fit_range;
                        auto corr_reso_costh_dir = histos_.find(folder+"/Resolution/Costh"+cat+hname);
                        auto corr_reso_mass_dir = histos_.find(folder+"/Resolution/Mass"+cat+hname);

                            // get corrected thad 4-vec and costhetastar vals
                        TLorentzVector Alpha_THad = alpha_corr_.Alpha_THad(perm, fit_deg, fit_var, fit_range);
                        double alpha_thad_cth = alpha_corr_.alpha_thad_cthstar(perm, fit_deg, fit_var, fit_range);

                        //TLorentzVector median_Alpha_THad = alpha_corr_.median_Alpha_THad(perm, fit_deg, fit_var, fit_range);

                        //cout << "   " << hname << ", alpha = " << alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range) << ", median alpha: " << alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range) << endl;
                        //cout << "       Corrected M(thad) Res: " << ttbar.had_top()->M()-Alpha_THad.M() << ", M(tt) Res: " << ttbar.M()-(Alpha_THad + perm.TLep()).M() << ", Reco M(thad): " << Alpha_THad.M() << ", Reco M(tt): " << (Alpha_THad + perm.TLep()).M() << ", pT: " << Alpha_THad.Pt() << ", Mass: " << Alpha_THad.M() << ", E: " << Alpha_THad.E() << ", P: " << Alpha_THad.P() << endl;
                        //cout << "       Median Corrected M(thad) Res: " << ttbar.had_top()->M()-median_Alpha_THad.M() << ", M(tt) Res: " << ttbar.M()-(median_Alpha_THad + perm.TLep()).M() << ", Reco M(thad): " << median_Alpha_THad.M() << ", Reco M(tt): " << (median_Alpha_THad + perm.TLep()).M() << ", pT: " << median_Alpha_THad.Pt() << ", Mass: " << median_Alpha_THad.M() << ", E: " << median_Alpha_THad.E() << ", P: " << median_Alpha_THad.P() << endl;


                            // investigate shoulder in Mthad resolution dists
                        if( ttbar.had_top()->M() - Alpha_THad.M() < 30 ){
                            fill_shoulder_resolution_plots( folder+"/Resolution/Mass"+cat+hname+"/L30", ttbar, perm, Alpha_THad, fit_deg, fit_var, fit_range );
                        }
                        else{
                            fill_shoulder_resolution_plots( folder+"/Resolution/Mass"+cat+hname+"/G30", ttbar, perm, Alpha_THad, fit_deg, fit_var, fit_range );
                        }

                            // fill plots for corrected thad
                        corr_reso_mass_dir->second["THad"].fill( ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        corr_reso_mass_dir->second["TTbar"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        corr_reso_mass_dir->second["Alpha_vs_Gen_MTTbar"].fill( ttbar.M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        //corr_reso_mass_dir->second["Median_Alpha_vs_Gen_MTTbar"].fill( ttbar.M(), alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        //corr_reso_mass_dir->second["Alpha_Diff_vs_Gen_MTTbar"].fill( ttbar.M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range)-alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        corr_reso_mass_dir->second["Computed_Alpha_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        //corr_reso_mass_dir->second["Median_Alpha_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        //corr_reso_mass_dir->second["Alpha_Diff_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), alpha_corr_.alpha(perm, fit_deg, fit_var, fit_range)-alpha_corr_.median_alpha(perm, fit_deg, fit_var, fit_range), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_MTTbar_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_MTHad_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_MTHad_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_MTHad_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_THadPt_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - Alpha_THad.Pt(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_THadEta_vs_Gen_THadEta"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - Alpha_THad.Eta(), evt_weight_ );
                        corr_reso_mass_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar"].fill( ttbar.Eta(), ttbar.Eta() - (Alpha_THad +perm.TLep()).Eta(), evt_weight_ );
                        //reso_mass_dir->second["Frac_THad_Corrected_"+hname].fill( (ttbar.had_top()->M() - Alpha_THad.M())/ttbar.had_top()->M(), evt_weight_ );
                        //reso_mass_dir->second["Frac_TTbar_Corrected_"+hname].fill( (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        //reso_mass_dir->second["Frac_TTbar_vs_Gen_THadPt_Corrected_"+hname].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        corr_reso_costh_dir->second["THad"].fill( gen_cosths.first - alpha_thad_cth, evt_weight_);

                        //if( ttbar.M() - (Alpha_THad + perm.TLep()).M() > ttbar.M()/3. ){
                        //    reso_dir->second["Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz(), evt_weight_ );
                        //}
                        //else{
                        //    reso_dir->second["Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt(), evt_weight_ );
                        //    reso_dir->second["Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz(), evt_weight_ );
                        //}
                    }
                }
            }

                // fill plots for original thad
            uncorr_reso_mass_dir->second["THad"].fill( ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["TTbar"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_MTTbar_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_MTHad_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_MTHad_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_MTHad_vs_Uncorrected_MTHad"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_THadPt_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_THadEta_vs_Gen_THadEta"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar"].fill( ttbar.Eta(), ttbar.Eta() - perm.LVect().Eta(), evt_weight_ );
            //uncorr_reso_mass_dir->second["Frac_THad"].fill( (ttbar.had_top()->M() - perm.THad().M())/ttbar.had_top()->M(), evt_weight_ );
            //uncorr_reso_mass_dir->second["Frac_TTbar"].fill( (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            //uncorr_reso_mass_dir->second["Frac_TTbar_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

            uncorr_reso_mass_dir->second["TLep"].fill( ttbar.lep_top()->M() - perm.TLep().M(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_TLepPt_vs_Gen_TLepPt"].fill( ttbar.lep_top()->Pt(), ttbar.lep_top()->Pt() - perm.TLep().Pt(), evt_weight_ );
            uncorr_reso_mass_dir->second["Reso_TLepEta_vs_Gen_TLepEta"].fill( ttbar.lep_top()->Eta(), ttbar.lep_top()->Eta() - perm.TLep().Eta(), evt_weight_ );

            uncorr_reso_costh_dir->second["THad"].fill( gen_cosths.first - reco_thad_cth, evt_weight_ );

        }

        void fill_parton_accept_plots( string folder, string cat, GenTTBar &ttbar, Permutation &perm ){
            auto reso_part_dir = histos_.find(folder+"/Resolution/Parton_Acceptance"+cat);


            vector<string> fit_degrees = {"1d", "2d"};
            vector<string> fit_vars    = {"E", "P", "M"};
            vector<string> fit_ranges    = {"All", "Mtt"};
            for( auto fit_deg : fit_degrees ){
                for( auto fit_var : fit_vars ){
                    for( auto fit_range : fit_ranges ){
                        string hname = fit_deg+"_"+fit_var+"_"+fit_range;

                            // get corrected thad 4-vec and costhetastar vals
                        TLorentzVector Alpha_THad = alpha_corr_.Alpha_THad(perm, fit_deg, fit_var, fit_range);
                        //double alpha_thad_cth = alpha_corr_.alpha_thad_cthstar(perm, fit_deg, fit_var);

                            // fill plots for corrected thad
                        reso_part_dir->second["Reso_MTTbar_All_Corrected_"+hname].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_All_Corrected_"+hname].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                        reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_All_Corrected_"+hname].fill( ttbar.M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_All_Corrected_"+hname].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                        reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_All_Corrected_"+hname].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - Alpha_THad.Pt(), evt_weight_ );
                        reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_All_Corrected_"+hname].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - Alpha_THad.Eta(), evt_weight_ );
                        reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_All_Corrected_"+hname].fill( ttbar.Eta(), ttbar.Eta() - (Alpha_THad + perm.TLep()).Eta(), evt_weight_ );
                        reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_All_Corrected_"+hname].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        if( ttbar.two_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                            reso_part_dir->second["Reso_MTTbar_2Partons_Corrected_"+hname].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_2Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_2Partons_Corrected_"+hname].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_2Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - Alpha_THad.Pt(), evt_weight_ );
                            reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_2Partons_Corrected_"+hname].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - Alpha_THad.Eta(), evt_weight_ );
                            reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_2Partons_Corrected_"+hname].fill( ttbar.Eta(), ttbar.Eta() - (Alpha_THad + perm.TLep()).Eta(), evt_weight_ );
                            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_2Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        }
                        if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                            reso_part_dir->second["Reso_MTTbar_3Partons_Corrected_"+hname].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_3Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_3Partons_Corrected_"+hname].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_3Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - Alpha_THad.Pt(), evt_weight_ );
                            reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_3Partons_Corrected_"+hname].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - Alpha_THad.Eta(), evt_weight_ );
                            reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_3Partons_Corrected_"+hname].fill( ttbar.Eta(), ttbar.Eta() - (Alpha_THad + perm.TLep()).Eta(), evt_weight_ );
                            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_3Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        }
                        if( ttbar.all_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                            reso_part_dir->second["Reso_MTTbar_4Partons_Corrected_"+hname].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_4Partons_Corrected_"+hname].fill( ttbar.M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_4Partons_Corrected_"+hname].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );
                            reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_4Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - Alpha_THad.Pt(), evt_weight_ );
                            reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_4Partons_Corrected_"+hname].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - Alpha_THad.Eta(), evt_weight_ );
                            reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_4Partons_Corrected_"+hname].fill( ttbar.Eta(), ttbar.Eta() - (Alpha_THad + perm.TLep()).Eta(), evt_weight_ );
                            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_4Partons_Corrected_"+hname].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
                        }
                    }
                }
            }


                    // reso breaking up by number of partons in acceptance
                        // combined categores
            reso_part_dir->second["Reso_MTTbar_All"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_All"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
            reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_All"].fill( ttbar.M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_All"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_All"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_ );
            reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_All"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_ );
            reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_All"].fill( ttbar.Eta(), ttbar.Eta() - perm.LVect().Eta(), evt_weight_ );
            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_All"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

            if( ttbar.two_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_2Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_2Partons"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_2Partons"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_ );
                reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_2Partons"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_ );
                reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_2Partons"].fill( ttbar.Eta(), ttbar.Eta() - perm.LVect().Eta(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_2Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            }
            if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_3Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_3Partons"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_3Partons"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_ );
                reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_3Partons"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_ );
                reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_3Partons"].fill( ttbar.Eta(), ttbar.Eta() - perm.LVect().Eta(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_3Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            }
            if( ttbar.all_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_4Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTHad_vs_Uncorrected_MTHad_4Partons"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
                reso_part_dir->second["Reso_THadPt_vs_Gen_THadPt_4Partons"].fill( ttbar.had_top()->Pt(), ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_ );
                reso_part_dir->second["Reso_THadEta_vs_Gen_THadEta_4Partons"].fill( ttbar.had_top()->Eta(), ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_ );
                reso_part_dir->second["Reso_EtaTTbar_vs_Gen_EtaTTbar_4Partons"].fill( ttbar.Eta(), ttbar.Eta() - perm.LVect().Eta(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_4Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            }

        }
        /// book and fill plots for hadronic top mass corrections for lost-jet events

        void fill_plots( string folder, string cat, GenTTBar &ttbar, Permutation &perm ){
            fill_reconstruction_plots( folder, cat, ttbar, perm );
            fill_resolution_plots( folder, cat, ttbar, perm );
            //fill_parton_accept_plots( folder, cat, ttbar, perm );
        }

        /// book and fill plots for matched permutations
        string matched_perm_categories( Permutation &mp ){
        //void fill_matched_perm_plots( Permutation &mp, Permutation &bp, GenTTBar &ttbar ){

            stringstream mp_evt_type;

            if( mp.unique_matches() == 0 ){
                mp_evt_type << "0";
            }
            else if( mp.unique_matches() == 1 ){
                mp_evt_type << "1";
            }
            else if( mp.unique_matches() == 2 ){
                mp_evt_type << "2";
                //cout << "mp: " << mp << endl;
            }
            else if( mp.unique_matches() == 3 ){
                mp_evt_type << "3/";

                if( mp.Merged_Event() ){
                    mp_evt_type << "Merged/";
                    if( mp.Merged_BHadBLep() ) mp_evt_type << "BHadBLep";
                    else if( mp.Merged_BHadWJa() ) mp_evt_type << "BHadWJa";
                    else if( mp.Merged_BHadWJb() ) mp_evt_type << "BHadWJb";
                    else if( mp.Merged_BLepWJa() ) mp_evt_type << "BLepWJa";
                    else if( mp.Merged_BLepWJb() ) mp_evt_type << "BLepWJb";
                    else{
                        mp_evt_type << "WJets";
                    }
                }
                else if( mp.Lost_Event() ){
                    mp_evt_type << "Lost/";
                    if( mp.Lost_BHad() ) mp_evt_type << "BHad";
                    else if( mp.Lost_BLep() ) mp_evt_type << "BLep";
                    else if( mp.Lost_WJa() ) mp_evt_type << "WJa";
                    else{
                        mp_evt_type << "WJb";
                    }
                }
                else{
                    Logger::log().fatal() << "This should be a lost or merged perm!" << endl;
                }
            }
            else{
                Logger::log().fatal() << "There are " << mp.unique_matches() << " for the matched permutation. This shouldn't be possible!" << endl;
            }

            //cout << mp_evt_type.str() << endl;
            return mp_evt_type.str();
            //return "3J_Event_Plots/Matched_Perm/Post_Alpha_Correction/"+mp_evt_type.str();
        }


        void book_gen_plots( string folder ){

            vector<string> gen_obs = { "BHad", "BLep", "WJa", "WJb" };

            for( auto gen_ob : gen_obs ){
                book<TH2D>(folder, "Gen_"+gen_ob+"_Pt_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen p_{T}", 20, 0., 2000., 100, 0., 1000.);
                book<TH2D>(folder, "Gen_"+gen_ob+"_Eta_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen |#eta|", 20, 0., 2000., 100, 0., 5.);
                book<TH2D>(folder, "Gen_"+gen_ob+"_Eta_vs_Pt", ";Gen p_{T};Gen |#eta|", 100, 0., 1000., 100, 0., 5.);
                book<TH2D>(folder+"/In_Acceptance", "Gen_"+gen_ob+"_Pt_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen p_{T}", 20, 0., 2000., 100, 0., 1000.);
                book<TH2D>(folder+"/In_Acceptance", "Gen_"+gen_ob+"_Eta_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen |#eta|", 20, 0., 2000., 100, 0., 5.);
                book<TH2D>(folder+"/In_Acceptance", "Gen_"+gen_ob+"_Eta_vs_Pt", ";Gen p_{T};Gen |#eta|", 100, 0., 1000., 100, 0., 5.);
                book<TH2D>(folder+"/Out_Acceptance", "Gen_"+gen_ob+"_Pt_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen p_{T}", 20, 0., 2000., 100, 0., 1000.);
                book<TH2D>(folder+"/Out_Acceptance", "Gen_"+gen_ob+"_Eta_vs_Gen_MTTbar", ";Gen M(t#bar{t});Gen |#eta|", 20, 0., 2000., 100, 0., 5.);
                book<TH2D>(folder+"/Out_Acceptance", "Gen_"+gen_ob+"_Eta_vs_Pt", ";Gen p_{T};Gen |#eta|", 100, 0., 1000., 100, 0., 5.);
            }

        }

        void fill_gen_plots( GenTTBar &ttbar ){

            auto gen_dir = histos_.find( "3J/Gen_Plots" );
            auto gen_in_dir = histos_.find( "3J/Gen_Plots/In_Acceptance" );
            auto gen_out_dir = histos_.find( "3J/Gen_Plots/Out_Acceptance" );

                // All Events
            gen_dir->second["Gen_BHad_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_b()->Pt(), evt_weight_ );
            gen_dir->second["Gen_BHad_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_b()->Eta()), evt_weight_ );
            gen_dir->second["Gen_BHad_Eta_vs_Pt"].fill( ttbar.had_b()->Pt(), Abs(ttbar.had_b()->Eta()), evt_weight_ );

            gen_dir->second["Gen_BLep_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.lep_b()->Pt(), evt_weight_ );
            gen_dir->second["Gen_BLep_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );
            gen_dir->second["Gen_BLep_Eta_vs_Pt"].fill( ttbar.lep_b()->Pt(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );

            gen_dir->second["Gen_WJa_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->first->Pt(), evt_weight_ );
            gen_dir->second["Gen_WJa_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );
            gen_dir->second["Gen_WJa_Eta_vs_Pt"].fill( ttbar.had_W()->first->Pt(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );

            gen_dir->second["Gen_WJb_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->second->Pt(), evt_weight_ );
            gen_dir->second["Gen_WJb_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );
            gen_dir->second["Gen_WJb_Eta_vs_Pt"].fill( ttbar.had_W()->second->Pt(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );

                // bhad in/out of acceptance
            if( ttbar.is_bhad_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ) {
                gen_in_dir->second["Gen_BHad_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_b()->Pt(), evt_weight_ );
                gen_in_dir->second["Gen_BHad_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_b()->Eta()), evt_weight_ );
                gen_in_dir->second["Gen_BHad_Eta_vs_Pt"].fill( ttbar.had_b()->Pt(), Abs(ttbar.had_b()->Eta()), evt_weight_ );
            }
            else{
                gen_out_dir->second["Gen_BHad_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_b()->Pt(), evt_weight_ );
                gen_out_dir->second["Gen_BHad_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_b()->Eta()), evt_weight_ );
                gen_out_dir->second["Gen_BHad_Eta_vs_Pt"].fill( ttbar.had_b()->Pt(), Abs(ttbar.had_b()->Eta()), evt_weight_ );
            }

                // blep in/out of acceptance
            if( ttbar.is_blep_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ) {
                gen_in_dir->second["Gen_BLep_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.lep_b()->Pt(), evt_weight_ );
                gen_in_dir->second["Gen_BLep_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );
                gen_in_dir->second["Gen_BLep_Eta_vs_Pt"].fill( ttbar.lep_b()->Pt(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );
            }
            else{
                gen_out_dir->second["Gen_BLep_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.lep_b()->Pt(), evt_weight_ );
                gen_out_dir->second["Gen_BLep_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );
                gen_out_dir->second["Gen_BLep_Eta_vs_Pt"].fill( ttbar.lep_b()->Pt(), Abs(ttbar.lep_b()->Eta()), evt_weight_ );
            }

                // wja in/out of acceptance
            if( ttbar.is_wja_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ) {
                gen_in_dir->second["Gen_WJa_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->first->Pt(), evt_weight_ );
                gen_in_dir->second["Gen_WJa_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );
                gen_in_dir->second["Gen_WJa_Eta_vs_Pt"].fill( ttbar.had_W()->first->Pt(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );
            }
            else{
                gen_out_dir->second["Gen_WJa_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->first->Pt(), evt_weight_ );
                gen_out_dir->second["Gen_WJa_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );
                gen_out_dir->second["Gen_WJa_Eta_vs_Pt"].fill( ttbar.had_W()->first->Pt(), Abs(ttbar.had_W()->first->Eta()), evt_weight_ );
            }

                // wjb in/out of acceptance
            if( ttbar.is_wjb_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ) {
                gen_in_dir->second["Gen_WJb_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->second->Pt(), evt_weight_ );
                gen_in_dir->second["Gen_WJb_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );
                gen_in_dir->second["Gen_WJb_Eta_vs_Pt"].fill( ttbar.had_W()->second->Pt(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );
            }
            else{
                gen_out_dir->second["Gen_WJb_Pt_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.had_W()->second->Pt(), evt_weight_ );
                gen_out_dir->second["Gen_WJb_Eta_vs_Gen_MTTbar"].fill( ttbar.M(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );
                gen_out_dir->second["Gen_WJb_Eta_vs_Pt"].fill( ttbar.had_W()->second->Pt(), Abs(ttbar.had_W()->second->Eta()), evt_weight_ );
            }

        }



        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;

            // LOST 3-jet events
            vector<string> lost_evt_type_categories = {"NO_MP", // categories for best perm
                "CORRECT_WJET_CORRECT_Bs", "CORRECT_WJET_SWAPPED_Bs", "CORRECT_WJET_CORRECT_BHAD", "CORRECT_WJET_CORRECT_BLEP", "CORRECT_WJET_WRONG_Bs",
                "WRONG_WJET_CORRECT_Bs", "WRONG_WJET_SWAPPED_Bs", "WRONG_WJET_CORRECT_BHAD", "WRONG_WJET_CORRECT_BLEP", "WRONG_WJET_WRONG_Bs",

                //"Lost_Jets", "Merged_Jets"
            };

            // matched perm
            vector<string> matched_perm_nunique_solutions = {
                "0", "1", "2", "3"
            };
            vector<string> matched_perm_3matches_evt_tyes = {
                "Lost/BHad", "Lost/BLep", "Lost/WJa", "Lost/WJb",
                "Merged/BHadBLep", "Merged/BHadWJa", "Merged/BHadWJb",
                "Merged/BLepWJa", "Merged/BLepWJb", "Merged/WJets"
            };

            for( auto& sys : systematics_ ){
                string sys_name = systematics::shift_to_name.at(sys);
                string lost_bp_dname = "3J/"+sys_name+"/Lost_BP";

                // plots for lost_bp from 3-jet events
                book_plots(lost_bp_dname+"/Post_Alpha_Correction", "" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                    book_plots(lost_bp_dname+"/Post_Alpha_Correction/"+lost_evt_type_cat, "" );
                }

                // plots for matched perm
                string mp_dname = "3J/"+sys_name+"/Matched_Perm";
                for( auto mp_solution : matched_perm_nunique_solutions ){
                    // plots for lost_bp from 3-jet events
                    if( mp_solution == "3" ){
                        for( auto evt_type : matched_perm_3matches_evt_tyes ){
                            book_plots(mp_dname+"/Post_Alpha_Correction", "/"+mp_solution+"/"+evt_type );
                        }
                    }
                    else{
                        book_plots(mp_dname+"/Post_Alpha_Correction", "/"+mp_solution );
                    }

                }
            }

            string gen_dir = "3J/Gen_Plots";
            book_gen_plots( gen_dir );

            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }


        // find costheta values for thad and tlep
        //reco
        std::pair< double, double > reco_costh_tops( Permutation &perm ){

            hyp::TTbar reco_ttang(perm); // doesn't work because not all wjets defined
            auto reco_ttcm = reco_ttang.to_CM();
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            double reco_tlep_cth = reco_ttang.unit3D().Dot(reco_ttcm.tlep().unit3D());

            return std::make_pair(reco_thad_cth, reco_tlep_cth);
        }

        //gen
        std::pair< double, double > gen_costh_tops( GenTTBar &ttbar ){

            hyp::TTbar gen_ttang(ttbar);
            auto gen_ttcm = gen_ttang.to_CM();
            double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());
            double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());

            return std::make_pair(gen_thad_cth, gen_tlep_cth);
        }

        // compute costhetastar for permutation thad
        double reco_thad_cthstar( Permutation &perm ){
            hyp::TTbar reco_ttang(perm); 
            auto reco_ttcm = reco_ttang.to_CM();
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            //double reco_thad_labframe_cth = reco_ttang.unit3D().Dot(perm.THad().Vect().Unit());

            return reco_thad_cth;
        }

        // top pT reweighting formula based on section 5.5 of AN-16-272
        double top_pt_reweighting( string sys_var ){
            double top_pt1, top_pt2;

            if( sys_var == "nominal" ){top_pt1 = 0.; top_pt2 = 0.;}
            else if( sys_var == "up_nom" ){top_pt1 = 1.; top_pt2 = 0.;}
            else if( sys_var == "down_nom" ){top_pt1 = -1.; top_pt2 = 0.;}
            else if( sys_var == "nom_up" ){top_pt1 = 0.; top_pt2 = 1.;}
            else if( sys_var == "nom_down" ){top_pt1 = 0.; top_pt2 = -1.;}
            else{
                Logger::log().error() << sys_var << " isn't a valid option for the top pt reweighting!" << endl;
                throw 40;
            }

            double p0 =  6.15025*pow(10., -2.) + 3.243*pow(10., -2.)*top_pt1 - 4.353*pow(10., -7.)*top_pt2;
            double p1 = -5.17833*pow(10., -4.) - 1.404*pow(10., -4.)*top_pt1 - 1.005*pow(10., -4.)*top_pt2;

            double weight = exp(p0+p1*( (genp_selector_.top()->Pt()+genp_selector_.tbar()->Pt())/2 ));
            return weight;
        }


        // find best perm assuming it's merged
        Permutation merged_best_perm( Permutation &perm, string disc_type ){

            Permutation merged_bp;
            double merged_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Merged(test_perm);

                double tp_discval = 3e10;

                if( disc_type == "Total" ){ tp_discval = test_perm.Prob(); }
                else if( disc_type == "NS" ){ tp_discval = test_perm.NuDiscr(); }
                else if( disc_type == "Mass" ){ tp_discval = test_perm.MassDiscr(); }
                else{
                    Logger::log().debug() << "Not a valid discriminant type" << endl;
                    throw 42;
                }

                if( tp_discval < merged_lowest_Totaldisc_3J ){
                    merged_lowest_Totaldisc_3J = tp_discval;
                    merged_bp = test_perm;
                }
            }
            return merged_bp;
        }

        // find best perm assuming it's lost 
        Permutation lost_best_perm( Permutation &perm, string disc_type ){

            Permutation lost_bp;
            double lost_lowest_discval = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Lost(test_perm);

                double tp_discval = 3e10;

                if( disc_type == "Total" ){ tp_discval = test_perm.Prob(); }
                else if( disc_type == "NS" ){ tp_discval = test_perm.NuDiscr(); }
                else if( disc_type == "Mass" ){ tp_discval = test_perm.MassDiscr(); }
                else{
                    Logger::log().debug() << "Not a valid discriminant type" << endl;
                    throw 42;
                }

                if( tp_discval < lost_lowest_discval ){
                    lost_lowest_discval = tp_discval;
                    lost_bp = test_perm;
                }
            }

            return lost_bp;
        }


        void lost_bp_cats( GenTTBar &ttbar, Permutation &lost_bp, string dname ){

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                tracker_.track("Not semilep events");
                return;
            }
            tracker_.track("semilep");

            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());

                // matched perm plots
            string mp_cat = matched_perm_categories( mp );
            fill_plots( dname+"/Matched_Perm/Post_Alpha_Correction", "/"+mp_cat, ttbar, lost_bp);

            string lost_perm_status;
            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) ) lost_perm_status = "NO_MP"; // check if mp exists, has bhad and blep and at least one wjet
            if( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ){ // lost_bp wjet matched to at least one of the mp wjets
                // check if mp and lost_bp misidentified b's
                if( mp.AreBsSame(lost_bp)  ) lost_perm_status = "CORRECT_WJET_CORRECT_Bs";
                else if( mp.AreBsFlipped(lost_bp) ) lost_perm_status = "CORRECT_WJET_SWAPPED_Bs";
                else if( mp.IsBHadCorrect(lost_bp) && !mp.IsBLepCorrect(lost_bp) ) lost_perm_status = "CORRECT_WJET_CORRECT_BHAD";
                else if( !mp.IsBHadCorrect(lost_bp) && mp.IsBLepCorrect(lost_bp) ) lost_perm_status = "CORRECT_WJET_CORRECT_BLEP";
                else lost_perm_status = "CORRECT_WJET_WRONG_Bs";
            }
            else{ // lost_bp wjet not matched to at least one of the mp wjets
                // check if mp and lost_bp misidentified b's
                if( mp.AreBsSame(lost_bp)  ) lost_perm_status = "WRONG_WJET_CORRECT_Bs";
                else if( mp.AreBsFlipped(lost_bp) ) lost_perm_status = "WRONG_WJET_SWAPPED_Bs";
                else if( mp.IsBHadCorrect(lost_bp) && !mp.IsBLepCorrect(lost_bp) ) lost_perm_status = "WRONG_WJET_CORRECT_BHAD";
                else if( !mp.IsBHadCorrect(lost_bp) && mp.IsBLepCorrect(lost_bp) ) lost_perm_status = "WRONG_WJET_CORRECT_BLEP";
                else lost_perm_status = "WRONG_WJET_WRONG_Bs";
            }

            //if( lost_bp.THad().M() > 180.0 ) return; // throw out events in which M(thad) from original thad is unphysical

            // lost_bp plots
            fill_plots( dname+"/Lost_BP/Post_Alpha_Correction", "", ttbar, lost_bp);

            // break dists up by classification comparing best_perm and matched_perm
            fill_plots( dname+"/Lost_BP/Post_Alpha_Correction/"+lost_perm_status, "", ttbar, lost_bp);



            fill_gen_plots( ttbar );
            //// break up dists based on matched perm being lost or merged
            //if( mp.Lost_Event() ){
            //    fill_plots( dname+"/Lost_BP/Post_Alpha_Correction/Lost_Jets", "", ttbar, lost_bp);
            //}
            //if( mp.Merged_Event() ){
            //    fill_plots( dname+"/Lost_BP/Post_Alpha_Correction/Merged_Jets", "", ttbar, lost_bp);
            //}

        } // end of lost_bp_cats


        void process_evt( systematics::SysShifts shift, URStreamer &event ){
            permutator_.reset_3J();

            //Fill truth of resp matrix
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            //select reco objects
            if( !object_selector_.select(event, shift) ) return;
            tracker_.track("obj selection");
            if( object_selector_.clean_jets().size() != 3 ) return;
            tracker_.track("3 jets selection");

            bool lep_is_tight = (object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTMU ||
                    object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTEL);

            //MC Weight for lepton selection
            float lep_weight=1;
            if(!isData_) {
                if(object_selector_.tight_muons().size() == 1) {
                    lep_weight = muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta());
                    // if(-1.3 < object_selector_.muon()->Eta() && object_selector_.muon()->Eta() < -1)
                    //  cout << "Mu " << *object_selector_.muon() << " weight: " << lep_weight << " prev: " << evt_weight_ << endl;
                }
                if(object_selector_.tight_electrons().size() == 1) lep_weight = electron_sf_.get_sf(object_selector_.electron()->Pt(), object_selector_.electron()->etaSC());
            }
            evt_weight_ *= lep_weight;
            tracker_.track("MC weights");

            //cut on btag
            auto &clean_jets = object_selector_.clean_jets();
            sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });
            if(!clean_jets[0]->BTagId(cut_tight_b_)) return;
            if(!clean_jets[1]->BTagId(cut_loose_b_)) return;

            // preselection
            bool preselection_pass = permutator_.preselection(
                    object_selector_.clean_jets(), object_selector_.lepton(),
                    object_selector_.met(), object_selector_.lepton_charge(),
                    event.rho().value(), lep_is_tight
                    );
            tracker_.track("permutation pre-selection done (not applied)");

            if( !preselection_pass ) return;
            tracker_.track("perm preselection");

            //find mc weight for btag
            double bweight = 1;
            if(!isData_){
                bweight *= btag_sf_.scale_factor(clean_jets, shift, false);
                //tracker_.track("not data and not sync");
            }
            evt_weight_ *= bweight;

            //find mc weight for top pt for ttbar events
            if( isTTbar_ ){
                double top_pt_weight = top_pt_reweighting("nominal");
                evt_weight_ *= top_pt_weight;
            }

            IDJet* wj2 = 0;
            Permutation event_perm(clean_jets[2], wj2, clean_jets[0], clean_jets[1], object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wja, wjb, bhad, blep
            Permutation best_perm = lost_best_perm( event_perm, "Total" );

            bool reco_3J_success = !best_perm.IsEmpty() && best_perm.Prob() < 1e9;
            if( !reco_3J_success ) return;

//            fill_gen_plots( ttbar );

            if( best_perm.THad().M() > 180.0 ) return; // throw out events in which M(thad) from original thad is unphysical

            string dname = "3J/"+systematics::shift_to_name.at(shift);
            lost_bp_cats( ttbar, best_perm, dname );

        }// end of process_evt




        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            opts::variables_map &values = URParser::instance().values();
            int report = values["report"].as<int>();
            Logger::log().debug() << "-- DONE -- reporting every -- " << report << " out of " << tree_->GetEntries() << endl;

            //int limit = 1000;
            int limit = -100;

            URStreamer event(tree_);

            while( event.next() )
            {

                if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << " run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;

                //if(evt_idx_ < 20) Logger::log().debug() << "Beginning event " << evt_idx_ << " eventnumber: " << event.evt << endl;

                if(limit > 0 && evt_idx_ > limit) {
                    return;
                }
                evt_idx_++;

                tracker_.track("All Events");

                //long and time consuming
                if(isTTbar_){
                    bool selection =    genp_selector_.select(event);
                    tracker_.track("gen selection");
                    if(!selection) {
                        Logger::log().error() << "Error: TTGenParticleSelector was not able to find all the generated top decay products in event " << evt_idx_ << endl <<
                            "run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
                        continue;
                    }
                }

                tracker_.deactivate();
                for(auto shift : systematics_){
                    evt_weight_ = mc_weights_.evt_weight(event, shift);
                    if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
                    //Logger::log().debug() << "processing: " << shift << endl;

                    process_evt(shift, event);
                    tracker_.group("");
                    if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
                }

            } // end of event loop

            Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
        } // end of analyze()

        //this method is called at the end of the job, by default saves
        //every histogram/tree produced, override it if you need something more
        virtual void end()
        {
            outFile_.Write();
            tracker_.writeTo(outFile_);
            Logger::log().debug() << "End of end() " << evt_idx_ << endl;
        }

        //do you need command-line or cfg options? If so implement this
        //method to book the options you need. CLI parsing is provided
        //by AnalysisFW/interface/URParser.h and uses boost::program_options
        //look here for a quickstart tutorial:
        //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html

        static void setOptions()
        {
            URParser &parser = URParser::instance();
            opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
            opts.add_options()
                ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                //("report", opts::value<int>()->default_value(1), "report every in debug mode");
                ("report", opts::value<int>()->default_value(10000), "report every in debug mode");
        }
};

//make it executable
int main(int argc, char *argv[])
{
    URParser &parser = URParser::instance(argc, argv);
    URDriver<ttbar_post_alpha_reco> test;
    int thing = test.run();

    opts::variables_map &values = parser.values();
    string output_file = values["output"].as<std::string>();
    Logger::log().debug() << "  Output File: " << output_file << endl;

    return thing;
}
