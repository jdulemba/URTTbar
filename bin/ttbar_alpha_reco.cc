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
using namespace systematics;
typedef SysShifts Sys;

class ttbar_alpha_reco : public AnalyzerBase
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
        bool isData_, isTTbar_, isSignal_, runsys_;

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

        float cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_;

    public:
        ttbar_alpha_reco(const std::string output_filename):
            AnalyzerBase("ttbar_alpha_reco", output_filename),
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

        //choose systematics to run based on sample
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

        Logger::log().debug() << "NOSYS: " << values.count("nosys") << endl;

        runsys_ = !( values.count("nosys") || isData_ );
        if( !runsys_ )
            systematics_ = {systematics::SysShifts::NOSYS};
        else{
            systematics_ = {
                Sys::NOSYS,
                Sys::JES_UP,  Sys::JES_DW,
                Sys::JER_UP,  Sys::JER_DW,
                Sys::MET_UP,  Sys::MET_DW,
                Sys::PU_UP,   Sys::PU_DW,
                Sys::BEFF_UP, Sys::BEFF_DW,
                Sys::BFAKE_UP, Sys::BFAKE_DW,
                Sys::LEPEFF_UP, Sys::LEPEFF_DW,
            };
            if(isTTbar_) {
                systematics_.push_back(Sys::HDAMP_UP);
                systematics_.push_back(Sys::HDAMP_DW);
                systematics_.push_back(Sys::RENORM_UP);
                systematics_.push_back(Sys::RENORM_DW);
                systematics_.push_back(Sys::FACTOR_UP);
                systematics_.push_back(Sys::FACTOR_DW);
                systematics_.push_back(Sys::RENFACTOR_UP);
                systematics_.push_back(Sys::RENFACTOR_DW);
            }
        }
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

        /// book and fill plots at gen-level
        void book_gen_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", 100, 200., 2000.);// nbins, mass_min, mass_max
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, 150., 200.);
            book<TH1F>(folder+"/Mass", "TLep", "", mass_bins_, 150., 200.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, costh_min_, costh_max_);
            book<TH2D>(folder+"/Costh", "Mmerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, mass_bins_, 0., 300.);
            book<TH2D>(folder+"/Costh", "Ptmerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(folder+"/Costh", "Etamerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -1.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, eta_min_, eta_max_);          
            //book<TH1F>(folder+"/Eta", "TTbar", "", 100, -2.5, 2.5);

        }
        void fill_gen_plots( string folder, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(ttbar.M(), evt_weight_);
            pt_dir->second["TTbar"].fill(ttbar.Pt(), evt_weight_);
            //eta_dir->second["TTbar"].fill(ttbar.Eta(), evt_weight_);

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            if( ttbar.had_top() ){
                double merged_parton_mass = -10.;
                double merged_parton_pt = -10.;
                double merged_parton_eta = -10.;
                if( ttbar.only_merged_bhadwja(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->first).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->first).Eta();
                }
                if( ttbar.only_merged_bhadwjb(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->second).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->second).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->second).Eta();
                }
                if( ttbar.only_merged_wjawjb(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->first).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->first).Eta();
                }

                mass_dir->second["THad"].fill(ttbar.had_top()->M(), evt_weight_);
                pt_dir->second["THad"].fill(ttbar.had_top()->Pt(), evt_weight_);
                pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P(), evt_weight_);
                pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P(), evt_weight_);
                eta_dir->second["THad"].fill(ttbar.had_top()->Eta(), evt_weight_);
                costh_dir->second["THad"].fill(gen_cosths.first, evt_weight_);
                costh_dir->second["Mmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_mass, evt_weight_);
                costh_dir->second["Ptmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_pt, evt_weight_);
                costh_dir->second["Etamerged_vs_Costh"].fill(gen_cosths.first, merged_parton_eta, evt_weight_);

            }
            if( ttbar.lep_top() ){
                //cout << "tlep m: " << ttbar.lep_top()->M() << endl;
                mass_dir->second["TLep"].fill(ttbar.lep_top()->M(), evt_weight_);
                pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt(), evt_weight_);
                eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta(), evt_weight_);
                costh_dir->second["TLep"].fill(gen_cosths.second, evt_weight_);

            }
        }
        /// book and fill plots at gen-level

        /// book and fill plots at reco-level
        void book_reco_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", mass_bins_, 200., 2000.);
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, 50., 250.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -1.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TTbar", "",eta_bins_, eta_min_, eta_max_);
        }
        void fill_reco_plots( string folder, Permutation &perm ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(perm.LVect().M(), evt_weight_);
            pt_dir->second["TTbar"].fill(perm.LVect().Pt(), evt_weight_);
            eta_dir->second["TTbar"].fill(perm.LVect().Eta(), evt_weight_);

            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(perm.THad().M(), evt_weight_);
            pt_dir->second["THad"].fill(perm.THad().Pt(), evt_weight_);
            pt_dir->second["THad_PT_P"].fill(perm.THad().Pt()/perm.THad().P(), evt_weight_);
            pt_dir->second["THad_PZ_P"].fill(perm.THad().Pz()/perm.THad().P(), evt_weight_);
            eta_dir->second["THad"].fill(perm.THad().Eta(), evt_weight_);
            costh_dir->second["THad"].fill(reco_cosths.first, evt_weight_);
            pt_dir->second["TLep"].fill(perm.TLep().Pt(), evt_weight_);
            eta_dir->second["TLep"].fill(perm.TLep().Eta(), evt_weight_);
            costh_dir->second["TLep"].fill(reco_cosths.second, evt_weight_);
        }
        /// book and fill plots at reco-level

        /// book and fill resolution plots
        void book_reso_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, -1000., 500.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, -2*costh_max_, 2*costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, -2*costh_max_, 2*costh_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.2, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -0.5, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, -2*eta_max_, 2*eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, -2*eta_max_, 2*eta_max_);
            //book<TH1F>(folder+"/Eta", "TTbar", "",eta_bins_, -2*eta_max_, 2*eta_max_);
        }
        void fill_reso_plots( string folder, Permutation &perm, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(ttbar.M() - perm.LVect().M(), evt_weight_);
            pt_dir->second["TTbar"].fill(ttbar.Pt() - perm.LVect().Pt(), evt_weight_);
            //eta_dir->second["TTbar"].fill(ttbar.Eta() - perm.LVect().Eta(), evt_weight_);

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(ttbar.had_top()->M() - perm.THad().M(), evt_weight_);
            pt_dir->second["THad"].fill(ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_);
            pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P() - perm.THad().Pt()/perm.THad().P(), evt_weight_);
            pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P() - perm.THad().Pz()/perm.THad().P(), evt_weight_);
            eta_dir->second["THad"].fill(ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_);
            costh_dir->second["THad"].fill(gen_cosths.first - reco_cosths.first, evt_weight_);
            pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt() - perm.TLep().Pt(), evt_weight_);
            eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta() - perm.TLep().Eta(), evt_weight_);
            costh_dir->second["TLep"].fill(gen_cosths.second - reco_cosths.second, evt_weight_);
        }
        /// book and fill resolution plots


        /// book and fill plots for permutation discriminants
        void book_disc_plots( string folder ){
            book<TH1F>(folder, "Massdisc", "", massdisc_bins_, massdisc_min_, massdisc_max_);
            book<TH1F>(folder, "NSdisc", "", nsdisc_bins_, nsdisc_min_, nsdisc_max_);
            book<TH1F>(folder, "Totaldisc", "", combdisc_bins_, combdisc_min_, combdisc_max_);
        }
        void fill_disc_plots( string folder, Permutation &perm ){
            auto dir = histos_.find(folder);

            dir->second["Massdisc"].fill(perm.MassDiscr(), evt_weight_);
            dir->second["NSdisc"].fill(perm.NuDiscr(), evt_weight_);
            dir->second["Totaldisc"].fill(perm.Prob(), evt_weight_);
        }
        /// book and fill plots for permutation discriminants


        /// book and fill plots for hadronic top mass corrections for lost-jet events
        void book_alpha_correction_plots( string folder ){

                // 3D plots of ( x=173.1/Mth, y=Mtt, z=alpha )
            book<TH3D>(folder+"/THad_P", "Alpha_THad_P_Mtt_vs_Mthad_vs_Alpha", ";173.1/M(t_{h}); M(t#bar{t}); #alpha_{P}= Gen P(t_{h})/Reco P(t_{h})", 32, 0.9, 2.5, 36, 200., 2000., 500, 0., 10.);
            book<TH3D>(folder+"/THad_E", "Alpha_THad_E_Mtt_vs_Mthad_vs_Alpha", ";173.1/M(t_{h}); M(t#bar{t}); #alpha_{E}= Gen E(t_{h})/Reco E(t_{h})", 32, 0.9, 2.5, 36, 200., 2000., 500, 0., 10.);
            book<TH3D>(folder+"/THad_M", "Alpha_THad_M_Mtt_vs_Mthad_vs_Alpha", ";173.1/M(t_{h}); M(t#bar{t}); #alpha_{M}= Gen M(t_{h})/Reco M(t_{h})", 32, 0.9, 2.5, 36, 200., 2000., 500, 0., 10.);


                // entire mass spectrum
            book<TH2D>(folder+"/THad_P", "Alpha_THad_P", ";173.1/M(t_{h}); #alpha_{P}= Gen P(t_{h})/Reco P(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);
            book<TH2D>(folder+"/THad_E", "Alpha_THad_E", ";173.1/M(t_{h}); #alpha_{E}= Gen E(t_{h})/Reco E(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);
            book<TH2D>(folder+"/THad_M", "Alpha_THad_M", ";173.1/M(t_{h}); #alpha_{M}= Gen M(t_{h})/Reco M(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);

                // parts of mass spectrum
            vector<string> mttbar_ranges = { "200to350", "350to400", "400to500", "500to700", "700to1000", "1000toInf"};
            for( auto m_range : mttbar_ranges ){
                book<TH2D>(folder+"/THad_P", "Alpha_THad_P_Mttbar"+m_range, ";173.1/M(t_{h}); #alpha_{P}= Gen P(t_{h})/Reco P(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);
                book<TH2D>(folder+"/THad_E", "Alpha_THad_E_Mttbar"+m_range, ";173.1/M(t_{h}); #alpha_{E}= Gen E(t_{h})/Reco E(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);
                book<TH2D>(folder+"/THad_M", "Alpha_THad_M_Mttbar"+m_range, ";173.1/M(t_{h}); #alpha_{M}= Gen M(t_{h})/Reco M(t_{h})", 32, 0.9, 2.5, 500, 0., 10.);
            }

            if( !boost::contains(folder, "nosys") ) return;

            vector<string> mthad_ranges = { "0.9to1.1", "1.1to1.3", "1.3to1.5", "1.5to1.7", "1.7to1.9", "1.9to2.1", "2.1to2.3", "2.3to2.5"};

            // gen vs reco plots for bins of 173.1/reco M(thad)
            for( auto m_range : mthad_ranges ){
                book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_"+m_range, ";Reco E(t_{h}) [GeV]; Gen E(t_{h}) [GeV]", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
                book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_"+m_range, ";Reco P(t_{h}) [GeV]; Gen P(t_{h}) [GeV]", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
                book<TH2D>(folder+"/THad_M", "Gen_vs_Reco_THadM_"+m_range, ";Reco M(t_{h}) [GeV]; Gen M(t_{h}) [GeV]", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
            }

        }

        void fill_alpha_correction_plots( string folder, GenTTBar &ttbar, Permutation &perm ){
            auto thad_E_dir = histos_.find(folder+"/THad_E");
            auto thad_P_dir = histos_.find(folder+"/THad_P");
            auto thad_M_dir = histos_.find(folder+"/THad_M");

            if( perm.THad().M() > 180.0 ) return;
            //dir->second["Reco_vs_Gen_MTTbar"].fill( ttbar.M(), perm.LVect().M(), evt_weight_ );

            thad_P_dir->second["Alpha_THad_P"].fill( 173.1/perm.THad().M(), ttbar.had_top()->P()/perm.THad().P(), evt_weight_ );//mthad taken from pdg
            thad_E_dir->second["Alpha_THad_E"].fill( 173.1/perm.THad().M(), ttbar.had_top()->E()/perm.THad().E(), evt_weight_ );
            thad_M_dir->second["Alpha_THad_M"].fill( 173.1/perm.THad().M(), ttbar.had_top()->M()/perm.THad().M(), evt_weight_ );

            thad_P_dir->second["Alpha_THad_P_Mtt_vs_Mthad_vs_Alpha"].fill( 173.1/perm.THad().M(), perm.LVect().M(), ttbar.had_top()->P()/perm.THad().P(), evt_weight_ );
            thad_E_dir->second["Alpha_THad_E_Mtt_vs_Mthad_vs_Alpha"].fill( 173.1/perm.THad().M(), perm.LVect().M(), ttbar.had_top()->E()/perm.THad().E(), evt_weight_ ); 
            thad_M_dir->second["Alpha_THad_M_Mtt_vs_Mthad_vs_Alpha"].fill( 173.1/perm.THad().M(), perm.LVect().M(), ttbar.had_top()->M()/perm.THad().M(), evt_weight_ ); 

            string mtt_range;
            if( 200. <= perm.LVect().M() && perm.LVect().M() < 350. ) mtt_range = "200to350";
            else if( 350. <= perm.LVect().M() && perm.LVect().M() < 400. ) mtt_range = "350to400";
            else if( 400. <= perm.LVect().M() && perm.LVect().M() < 500. ) mtt_range = "400to500";
            else if( 500. <= perm.LVect().M() && perm.LVect().M() < 700. ) mtt_range = "500to700";
            else if( 700. <= perm.LVect().M() && perm.LVect().M() < 1000. ) mtt_range = "700to1000";
            else if( perm.LVect().M() >= 1000. ) mtt_range = "1000toInf";

            if( perm.LVect().M() >= 200. ){
                thad_P_dir->second["Alpha_THad_P_Mttbar"+mtt_range].fill( 173.1/perm.THad().M(), ttbar.had_top()->P()/perm.THad().P(), evt_weight_ );
                thad_E_dir->second["Alpha_THad_E_Mttbar"+mtt_range].fill( 173.1/perm.THad().M(), ttbar.had_top()->E()/perm.THad().E(), evt_weight_ );
                thad_M_dir->second["Alpha_THad_M_Mttbar"+mtt_range].fill( 173.1/perm.THad().M(), ttbar.had_top()->M()/perm.THad().M(), evt_weight_ );
            }

            if( !boost::contains(folder, "nosys") ) return;

            string mthad_range;
            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ) mthad_range = "0.9to1.1";
            if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ) mthad_range = "1.1to1.3";
            if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ) mthad_range = "1.3to1.5"; 
            if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ) mthad_range = "1.5to1.7"; 
            if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ) mthad_range = "1.7to1.9"; 
            if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ) mthad_range = "1.9to2.1"; 
            if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ) mthad_range = "2.1to2.3"; 
            if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ) mthad_range = "2.3to2.5"; 

            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_"+mthad_range].fill( perm.THad().E(), ttbar.had_top()->E(), evt_weight_ );
                thad_P_dir->second["Gen_vs_Reco_THadP_"+mthad_range].fill( perm.THad().P(), ttbar.had_top()->P(), evt_weight_ );
                thad_M_dir->second["Gen_vs_Reco_THadM_"+mthad_range].fill( perm.THad().M(), ttbar.had_top()->M(), evt_weight_ );
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
                "WRONG_WJET_CORRECT_Bs", "WRONG_WJET_SWAPPED_Bs", "WRONG_WJET_CORRECT_BHAD", "WRONG_WJET_CORRECT_BLEP", "WRONG_WJET_WRONG_Bs"
            };

            for( auto& sys : systematics_ ){
                string sys_name = systematics::shift_to_name.at(sys);
                string dname = "3J/"+sys_name;
                book_alpha_correction_plots(dname+"/Alpha_Correction" );

                if( sys != Sys::NOSYS ) continue;

                book_disc_plots(dname+"/Discr" );
                book_gen_plots( dname+"/Gen" );
                book_reco_plots(dname+"/Reconstruction" );
                book_reso_plots(dname+"/Resolution" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                    // plots for lost_bp from 3-jet events
                    book_disc_plots(dname+"/Discr/"+lost_evt_type_cat );
                    book_gen_plots( dname+"/Gen/"+lost_evt_type_cat );
                    book_reco_plots(dname+"/Reconstruction/"+lost_evt_type_cat );
                    book_reso_plots(dname+"/Resolution/"+lost_evt_type_cat );
                    book_alpha_correction_plots(dname+"/Alpha_Correction/"+lost_evt_type_cat );
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


            // lost_bp plots
            fill_alpha_correction_plots(dname+"/Alpha_Correction", ttbar, lost_bp); // only fill alpha correction plots for systematics

            if( !boost::contains(dname, "nosys") ) return;

            fill_disc_plots(dname+"/Discr", lost_bp );
            fill_gen_plots( dname+"/Gen", ttbar );
            fill_reco_plots(dname+"/Reconstruction", lost_bp );
            fill_reso_plots(dname+"/Resolution", lost_bp, ttbar );

            // break dists up by classification comparing best_perm and matched_perm
            fill_disc_plots(dname+"/Discr/"+lost_perm_status, lost_bp );
            fill_gen_plots( dname+"/Gen/"+lost_perm_status, ttbar );
            fill_reco_plots(dname+"/Reconstruction/"+lost_perm_status, lost_bp );
            fill_reso_plots(dname+"/Resolution/"+lost_perm_status, lost_bp, ttbar );
            fill_alpha_correction_plots(dname+"/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp);

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
                ("nosys", opts::value<int>()->default_value(1), "do not run systematics")
                ("report", opts::value<int>()->default_value(10000), "report every in debug mode");
        }
};

//make it executable
int main(int argc, char *argv[])
{
    URParser &parser = URParser::instance(argc, argv);
    URDriver<ttbar_alpha_reco> test;
    int thing = test.run();

    opts::variables_map &values = parser.values();
    string output_file = values["output"].as<std::string>();
    Logger::log().debug() << "  Output File: " << output_file << endl;

    return thing;
}
