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
        float alpha_correction_slope_, alpha_correction_yint_;

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
        parser.addCfgParameter<float>("alpha_correction", "slope", "slope to be used in alpha correction");
        parser.addCfgParameter<float>("alpha_correction", "yint", "y-int to be used in alpha correction");
        parser.parseArguments();

        cut_jet_ptmin_ = parser.getCfgPar<float>("jets", "ptmin" );
        cut_jet_etamax_ = parser.getCfgPar<float>("jets", "etamax");
        cut_leadjet_ptmin_ = parser.getCfgPar<float>("jets", "lead_ptmin" );

        alpha_correction_slope_ = parser.getCfgPar<float>("alpha_correction", "slope" );
        alpha_correction_yint_ = parser.getCfgPar<float>("alpha_correction", "yint" );


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
        void book_post_alpha_correction_plots( string folder ){
            // reco corrected and uncorrected plots
            book<TH1F>(folder+"/Reconstruction/Mass", "THad", "", 500, 0., 500.);
            book<TH1F>(folder+"/Reconstruction/Mass", "THad_Corrected", "", 500, 0., 500.);
            book<TH1F>(folder+"/Reconstruction/Mass", "TTbar", "", mass_bins_, 200., 2000.);
            book<TH1F>(folder+"/Reconstruction/Mass", "TTbar_Corrected", "", mass_bins_, 200., 2000.);
            book<TH2D>(folder+"/Reconstruction/Mass", "Reco_vs_Gen_TTbar", ";Gen M(ttbar); Reco M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);
            book<TH2D>(folder+"/Reconstruction/Mass", "Reco_vs_Gen_TTbar_Corrected", ";Gen M(ttbar); Reco Corrected M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);

            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Labframe", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Labframe_Corrected", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Corrected", "", costh_bins_, -1., 1.);
            book<TH2D>(folder+"/Reconstruction/Costh", "Reco_vs_Gen_THad", ";Gen Cos(theta)(t_{h}); Reco Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);
            book<TH2D>(folder+"/Reconstruction/Costh", "Reco_vs_Gen_THad_Corrected", ";Gen Cos(theta)(t_{h}); Reco Corrected Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);

            // reso corrected and uncorrected plots
            book<TH1F>(folder+"/Resolution/Mass", "THad", "", 500, -1000., 500.);
            book<TH1F>(folder+"/Resolution/Mass", "THad_Corrected", "", 500, -1000., 500.);
            book<TH1F>(folder+"/Resolution/Mass", "TTbar", "", 500, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "TTbar_Corrected", "", 500, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass", "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass", "Reso_MTTbar_vs_Gen_MTTbar_Corrected", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_THad", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_THad_Corrected", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_TTbar", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_TTbar_Corrected", "", mass_bins_, -2., 2.);
            book<TH2D>(folder+"/Resolution/Mass", "Frac_TTbar_vs_Gen_THadPt", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH2D>(folder+"/Resolution/Mass", "Frac_TTbar_vs_Gen_THadPt_Corrected", "", 100, 0., 1000., mass_bins_, -2., 2.);

            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_All", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_All", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_All", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_2Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_2Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_2Partons", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_3Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_3Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_3Partons", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_4Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_4Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_4Partons", "", 100, 0., 1000., mass_bins_, -2., 2.);

            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_All_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_All_Corrected", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_All_Corrected", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_2Partons_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_2Partons_Corrected", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_2Partons_Corrected", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_3Partons_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_3Partons_Corrected", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_3Partons_Corrected", "", 100, 0., 1000., mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_4Partons_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_4Partons_Corrected", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Frac_MTTbar_vs_Gen_THadPt_4Partons_Corrected", "", 100, 0., 1000., mass_bins_, -2., 2.);

            book<TH1F>(folder+"/Resolution/Costh", "THad", "", costh_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Costh", "THad_Corrected", "", costh_bins_, -2., 2.);


            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution", "Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);

        }

        void fill_post_alpha_correction_plots( string folder, GenTTBar &ttbar, Permutation &perm ){
            auto reco_mass_dir = histos_.find(folder+"/Reconstruction/Mass");
            auto reco_costh_dir = histos_.find(folder+"/Reconstruction/Costh");

            auto reso_dir = histos_.find(folder+"/Resolution");
            auto reso_costh_dir = histos_.find(folder+"/Resolution/Costh");
            auto reso_mass_dir = histos_.find(folder+"/Resolution/Mass");
            auto reso_part_dir = histos_.find(folder+"/Resolution/Parton_Acceptance");

                //alphas found by fitting Alpha_THad_P/E hists
                    //values taken from 1degree vals ttp://home.fnal.gov/~jdulemba/Plots/ttbar_post_alpha_reco/2018/Compare_Lost_Merged_Jets/Full/ttJetsM0/3J_Event_Plots/Final_Reco/Clear_and_MassCut_Classes/Class_Lost/Alpha_Correction/fit_parameters.json
            // y = mx +b -> m is first element in json, b is second

            double alpha_E = alpha_correction_slope_*( 173.1/perm.THad().M() ) + alpha_correction_yint_;
            //double alpha_E = 0.4019*( 173.1/perm.THad().M() ) + 0.5834; // only alpha_E used because it's more consistent over mtt spectrum 
            //double alpha_P = 0.1544*( 173.1/perm.THad().M() ) + 0.8599;

            //TLorentzVector Alpha_THad(alpha_P*perm.THad().Px(), alpha_P*perm.THad().Py(), alpha_P*perm.THad().Pz(), alpha_E*perm.THad().E());
            TLorentzVector Alpha_THad(alpha_E*perm.THad().Px(), alpha_E*perm.THad().Py(), alpha_E*perm.THad().Pz(), alpha_E*perm.THad().E());


            // costheta variables
            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
                // perm thad
            hyp::TTbar reco_ttang(perm); // doesn't work because not all wjets defined
            auto reco_ttcm = reco_ttang.to_CM();
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            double reco_thad_labframe_cth = reco_ttang.unit3D().Dot(perm.THad().Vect().Unit());
            //cout << "reco LF cth: " << reco_thad_labframe_cth << endl;

                // corrected thad
            TLorentzVector reco_corr_ttang = Alpha_THad+perm.TLep();
            TLorentzVector reco_ttcm_corr_thad = Alpha_THad;
            reco_ttcm_corr_thad.Boost(-1*reco_corr_ttang.BoostVector());
            double reco_corr_thad_cth = reco_corr_ttang.Vect().Unit().Dot(reco_ttcm_corr_thad.Vect().Unit());
            double reco_corr_labframe_thad_cth = reco_corr_ttang.Vect().Unit().Dot(Alpha_THad.Vect().Unit());
            //cout << "LF corr cth: " << reco_corr_labframe_cth << endl;

              // reco plots
                  // mass plots
            reco_mass_dir->second["THad"].fill( perm.THad().M(), evt_weight_ );
            reco_mass_dir->second["THad_Corrected"].fill( Alpha_THad.M(), evt_weight_ );
            reco_mass_dir->second["TTbar"].fill( perm.LVect().M(), evt_weight_ );
            reco_mass_dir->second["TTbar_Corrected"].fill( (Alpha_THad + perm.TLep()).M(), evt_weight_ );
            reco_mass_dir->second["Reco_vs_Gen_TTbar"].fill( ttbar.M(), perm.LVect().M(), evt_weight_ );
            reco_mass_dir->second["Reco_vs_Gen_TTbar_Corrected"].fill( ttbar.M(), (Alpha_THad + perm.TLep()).M(), evt_weight_ );

                  // costh plots
            reco_costh_dir->second["THad_Labframe"].fill( reco_thad_labframe_cth, evt_weight_ );
            reco_costh_dir->second["THad_Labframe_Corrected"].fill( reco_corr_labframe_thad_cth, evt_weight_ );
            reco_costh_dir->second["THad"].fill( reco_thad_cth, evt_weight_ );
            reco_costh_dir->second["THad_Corrected"].fill( reco_corr_thad_cth, evt_weight_ );
            reco_costh_dir->second["Reco_vs_Gen_THad"].fill(gen_cosths.first, reco_thad_cth, evt_weight_);
            reco_costh_dir->second["Reco_vs_Gen_THad_Corrected"].fill( gen_cosths.first, reco_corr_thad_cth, evt_weight_ );

                // reso plots
                // reso plots
            reso_mass_dir->second["THad"].fill( ttbar.had_top()->M() - perm.THad().M(), evt_weight_ );
            reso_mass_dir->second["THad_Corrected"].fill( ttbar.had_top()->M() - Alpha_THad.M(), evt_weight_ );

            reso_mass_dir->second["TTbar"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
            reso_mass_dir->second["Reso_MTTbar_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );

            reso_mass_dir->second["TTbar_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
            reso_mass_dir->second["Reso_MTTbar_vs_Gen_MTTbar_Corrected"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );

            reso_mass_dir->second["Frac_THad"].fill( (ttbar.had_top()->M() - perm.THad().M())/ttbar.had_top()->M(), evt_weight_ );
            reso_mass_dir->second["Frac_THad_Corrected"].fill( (ttbar.had_top()->M() - Alpha_THad.M())/ttbar.had_top()->M(), evt_weight_ );
            reso_mass_dir->second["Frac_TTbar"].fill( (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            reso_mass_dir->second["Frac_TTbar_Corrected"].fill( (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
            reso_mass_dir->second["Frac_TTbar_vs_Gen_THadPt"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );
            reso_mass_dir->second["Frac_TTbar_vs_Gen_THadPt_Corrected"].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );

                    // reso breaking up by number of partons in acceptance
                        // combined categores
            reso_part_dir->second["Reso_MTTbar_All"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_All"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_All"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

            reso_part_dir->second["Reso_MTTbar_All_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
            reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_All_Corrected"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
            reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_All_Corrected"].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );

            if( ttbar.two_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_2Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_2Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

                reso_part_dir->second["Reso_MTTbar_2Partons_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons_Corrected"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_2Partons_Corrected"].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
            }
            if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_3Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_3Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

                reso_part_dir->second["Reso_MTTbar_3Partons_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons_Corrected"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_3Partons_Corrected"].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
            }
            if( ttbar.all_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_4Partons"].fill( ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_4Partons"].fill( ttbar.had_top()->Pt(), (ttbar.M() - perm.LVect().M())/ttbar.M(), evt_weight_ );

                reso_part_dir->second["Reso_MTTbar_4Partons_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons_Corrected"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M(), evt_weight_ );
                reso_part_dir->second["Frac_MTTbar_vs_Gen_THadPt_4Partons_Corrected"].fill( ttbar.had_top()->Pt(), (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M(), evt_weight_ );
            }
            reso_costh_dir->second["THad"].fill( gen_cosths.first - reco_thad_cth, evt_weight_ );
            reso_costh_dir->second["THad_Corrected"].fill( gen_cosths.first - reco_corr_thad_cth, evt_weight_);


            if( ttbar.M() - (Alpha_THad + perm.TLep()).M() > ttbar.M()/3. ){
                reso_dir->second["Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz(), evt_weight_ );
            }
            else{
                reso_dir->second["Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt(), evt_weight_ );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz(), evt_weight_ );
            }

        }
        /// book and fill plots for hadronic top mass corrections for lost-jet events



        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;

            // LOST 3-jet events
            vector<string> lost_evt_type_categories = {"CORRECT_B", "WRONG_B", "OTHER", // categories for best perm and matched perm
                "CORRECT_BHAD", "CORRECT_BLEP", "CORRECT_Bs", "SWAPPED_Bs", "OTHER_MATCH" // categories for perm object matches and gen objects
            };
            vector<string> lost_best_perm_solutions = {
                "Lost_BP",
            };

            for( auto lost_bp_solution : lost_best_perm_solutions ){

                // plots for lost_bp from 3-jet events
                book_post_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Post_Alpha_Correction" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                    book_post_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Post_Alpha_Correction/"+lost_evt_type_cat );
                }
            }

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


        void lost_bp_cats( GenTTBar &ttbar, Permutation &lost_bp, string bp_solution ){

            string lost_perm_status;

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

            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) ) lost_perm_status = "OTHER"; // check if mp exists, has bhad and blep and at least one wjet
            // check if mp and lost_bp misidentified b's
            else if( mp.AreBsSame(lost_bp)  ) lost_perm_status = "CORRECT_B";
            else if( mp.AreBsFlipped(lost_bp) ) lost_perm_status = "WRONG_B";
            else lost_perm_status = "OTHER";


            string gen_match_status;
            if( lost_bp.BHad()->match() == ttbar.had_b() && lost_bp.BLep()->match() != ttbar.lep_b() ) gen_match_status = "CORRECT_BHAD"; // only perm bhad match is same as gen bhad
            else if( lost_bp.BHad()->match() != ttbar.had_b() && lost_bp.BLep()->match() == ttbar.lep_b() ) gen_match_status = "CORRECT_BLEP"; // only perm blep match is same as gen blep
            else if( lost_bp.BHad()->match() == ttbar.had_b() && lost_bp.BLep()->match() == ttbar.lep_b() ) gen_match_status = "CORRECT_Bs"; // both perm b matches same as gen bs
            else if( lost_bp.BHad()->match() == ttbar.lep_b() && lost_bp.BLep()->match() == ttbar.had_b() ) gen_match_status = "SWAPPED_Bs"; // perm bhad match same as gen blep, perm blep match same as gen bhad
            else gen_match_status = "OTHER_MATCH";


            // lost_bp plots
            fill_post_alpha_correction_plots("3J_Event_Plots/"+bp_solution+"/Post_Alpha_Correction", ttbar, lost_bp);

            // break dists up by classification comparing best_perm and matched_perm
            fill_post_alpha_correction_plots("3J_Event_Plots/"+bp_solution+"/Post_Alpha_Correction/"+lost_perm_status, ttbar, lost_bp);

            // break dists up by classification comparing best_perm matches and gen objects
            fill_post_alpha_correction_plots("3J_Event_Plots/"+bp_solution+"/Post_Alpha_Correction/"+gen_match_status, ttbar, lost_bp);

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
                    event.rho().value()
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
            lost_bp_cats( ttbar, best_perm, "Lost_BP" );

        }// end of process_evt




        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            opts::variables_map &values = URParser::instance().values();
            int report = values["report"].as<int>();
            Logger::log().debug() << "-- DONE -- reporting every -- " << report << " out of " << tree_->GetEntries() << endl;

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
