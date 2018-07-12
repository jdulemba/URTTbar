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

class ttbar_final_reco_3J : public AnalyzerBase
{

    public:
        //enum TTNaming_Merged {RIGHT, MERGE_SWAP, MERGE, WRONG};

    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        //	    double njets_ = 0;

        double disc_cut_ = 2.;


        // bins for hists
        int mass_bins_ = 100;
        int pt_bins_ = 50;
        int eta_bins_ = 50;
        int costh_bins_ = 50;

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
        bool isData_, isTTbar_;

        //selectors and helpers
        TTObjectSelector object_selector_; //selects ttbar objects
        TTGenParticleSelector genp_selector_; //selects generator level objects
        TTBarSolver solver_; //solves ttbar events
        TTGenMatcher matcher_; //matches particles on generator level
        TTPermutator permutator_;

        DR_TTGenMatcher dr_matcher_;

        float evt_weight_;
        MCWeightProducer mc_weights_;
        IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
        IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
        IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;

        float cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_;


    public:
        ttbar_final_reco_3J(const std::string output_filename):
            AnalyzerBase("ttbar_final_reco_3J", output_filename),
            tracker_(),
            object_selector_(),
            genp_selector_(TTGenParticleSelector::SelMode::LHE), //for ttJets files
            //genp_selector_(), // for other files
            solver_(),
            evt_weight_(1.),
            mc_weights_()
    {
        //set tracker
        tracker_.use_weight(&evt_weight_);
        //find out which sample we're running on

        opts::variables_map &values = URParser::instance().values();
        string output_file = values["output"].as<std::string>();
        string sample = systematics::get_sample(output_file);
        bool isSignal = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
        isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal;

        isData_  = boost::starts_with(sample, "data");
        Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << ", isSignal: " << isSignal << endl;
        if(isData_) {
            if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
            else object_selector_.lepton_type(1);
        }
        if(!isData_) mc_weights_.init(sample);

        // get parameters from cfg file
        URParser& parser = URParser::instance();

        parser.addCfgParameter<float>("gen_jets", "ptmin", "minimum pt");
        parser.addCfgParameter<float>("gen_jets", "etamax", "maximum eta");
        parser.addCfgParameter<float>("gen_jets", "lead_ptmin", "minimum leading jet pt");
        parser.parseArguments();

        cut_jet_ptmin_ = parser.getCfgPar<float>("gen_jets", "ptmin" );
        cut_jet_etamax_ = parser.getCfgPar<float>("gen_jets", "etamax");
        cut_leadjet_ptmin_ = parser.getCfgPar<float>("gen_jets", "lead_ptmin" );

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

            mass_dir->second["TTbar"].fill(ttbar.M());
            pt_dir->second["TTbar"].fill(ttbar.Pt());
            //eta_dir->second["TTbar"].fill(ttbar.Eta());

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

                mass_dir->second["THad"].fill(ttbar.had_top()->M());
                pt_dir->second["THad"].fill(ttbar.had_top()->Pt());
                pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P());
                pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P());
                eta_dir->second["THad"].fill(ttbar.had_top()->Eta());
                costh_dir->second["THad"].fill(gen_cosths.first);
                costh_dir->second["Mmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_mass);
                costh_dir->second["Ptmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_pt);
                costh_dir->second["Etamerged_vs_Costh"].fill(gen_cosths.first, merged_parton_eta);

            }
            if( ttbar.lep_top() ){
                //cout << "tlep m: " << ttbar.lep_top()->M() << endl;
                mass_dir->second["TLep"].fill(ttbar.lep_top()->M());
                pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt());
                eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta());
                costh_dir->second["TLep"].fill(gen_cosths.second);

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

            mass_dir->second["TTbar"].fill(perm.LVect().M());
            pt_dir->second["TTbar"].fill(perm.LVect().Pt());
            eta_dir->second["TTbar"].fill(perm.LVect().Eta());

            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(perm.THad().M());
            pt_dir->second["THad"].fill(perm.THad().Pt());
            pt_dir->second["THad_PT_P"].fill(perm.THad().Pt()/perm.THad().P());
            pt_dir->second["THad_PZ_P"].fill(perm.THad().Pz()/perm.THad().P());
            eta_dir->second["THad"].fill(perm.THad().Eta());
            costh_dir->second["THad"].fill(reco_cosths.first);
            pt_dir->second["TLep"].fill(perm.TLep().Pt());
            eta_dir->second["TLep"].fill(perm.TLep().Eta());
            costh_dir->second["TLep"].fill(reco_cosths.second);
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

            mass_dir->second["TTbar"].fill(ttbar.M() - perm.LVect().M());
            pt_dir->second["TTbar"].fill(ttbar.Pt() - perm.LVect().Pt());
            //eta_dir->second["TTbar"].fill(ttbar.Eta() - perm.LVect().Eta());

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(ttbar.had_top()->M() - perm.THad().M());
            pt_dir->second["THad"].fill(ttbar.had_top()->Pt() - perm.THad().Pt());
            pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P() - perm.THad().Pt()/perm.THad().P());
            pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P() - perm.THad().Pz()/perm.THad().P());
            eta_dir->second["THad"].fill(ttbar.had_top()->Eta() - perm.THad().Eta());
            costh_dir->second["THad"].fill(gen_cosths.first - reco_cosths.first);
            pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt() - perm.TLep().Pt());
            eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta() - perm.TLep().Eta());
            costh_dir->second["TLep"].fill(gen_cosths.second - reco_cosths.second);
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

            dir->second["Massdisc"].fill(perm.MassDiscr());
            dir->second["NSdisc"].fill(perm.NuDiscr());
            dir->second["Totaldisc"].fill(perm.Prob());
        }
        /// book and fill plots for permutation discriminants


        /// book and fill plots for hadronic top mass corrections for lost-jet events
        void book_alpha_correction_plots( string folder ){
            book<TH2D>(folder+"/THad_P", "Alpha_THad_P", "", 32, 0.9, 2.5, 500, 0., 50.);
            book<TH2D>(folder+"/THad_E", "Alpha_THad_E", "", 32, 0.9, 2.5, 500, 0., 50.);

            // gen vs reco plots for bins of 173.1/reco M(thad)
                // Energy
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_0.9to1.1", "", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_1.1to1.3", "", 75, 0., 1500., 100, 0., 2000.); // 1.1 < 173.1/reco M(thad) < 1.3
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_1.3to1.5", "", 75, 0., 1500., 100, 0., 2000.); // 1.3 < 173.1/reco M(thad) < 1.5
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_1.5to1.7", "", 75, 0., 1500., 100, 0., 2000.); // 1.5 < 173.1/reco M(thad) < 1.7
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_1.7to1.9", "", 75, 0., 1500., 100, 0., 2000.); // 1.7 < 173.1/reco M(thad) < 1.9
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_1.9to2.1", "", 75, 0., 1500., 100, 0., 2000.); // 1.9 < 173.1/reco M(thad) < 2.1
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_2.1to2.3", "", 75, 0., 1500., 100, 0., 2000.); // 2.1 < 173.1/reco M(thad) < 2.3
            book<TH2D>(folder+"/THad_E", "Gen_vs_Reco_THadE_2.3to2.5", "", 75, 0., 1500., 100, 0., 2000.); // 2.3 < 173.1/reco M(thad) < 2.5
               // Momentum
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_0.9to1.1", "", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_1.1to1.3", "", 75, 0., 1500., 100, 0., 2000.); // 1.1 < 173.1/reco M(thad) < 1.3
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_1.3to1.5", "", 75, 0., 1500., 100, 0., 2000.); // 1.3 < 173.1/reco M(thad) < 1.5
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_1.5to1.7", "", 75, 0., 1500., 100, 0., 2000.); // 1.5 < 173.1/reco M(thad) < 1.7
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_1.7to1.9", "", 75, 0., 1500., 100, 0., 2000.); // 1.7 < 173.1/reco M(thad) < 1.9
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_1.9to2.1", "", 75, 0., 1500., 100, 0., 2000.); // 1.9 < 173.1/reco M(thad) < 2.1
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_2.1to2.3", "", 75, 0., 1500., 100, 0., 2000.); // 2.1 < 173.1/reco M(thad) < 2.3
            book<TH2D>(folder+"/THad_P", "Gen_vs_Reco_THadP_2.3to2.5", "", 75, 0., 1500., 100, 0., 2000.); // 2.3 < 173.1/reco M(thad) < 2.5

        }

        void fill_alpha_correction_plots( string folder, GenTTBar &ttbar, Permutation &perm ){
            auto thad_E_dir = histos_.find(folder+"/THad_E");
            auto thad_P_dir = histos_.find(folder+"/THad_P");

            if( perm.THad().M() > 180.0 ) return;
            //dir->second["Reco_vs_Gen_MTTbar"].fill( ttbar.M(), perm.LVect().M() );

            thad_P_dir->second["Alpha_THad_P"].fill( 173.1/perm.THad().M(), ttbar.had_top()->P()/perm.THad().P() );//mthad taken from pdg
            thad_E_dir->second["Alpha_THad_E"].fill( 173.1/perm.THad().M(), ttbar.had_top()->E()/perm.THad().E() );

            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_0.9to1.1"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_0.9to1.1"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_1.1to1.3"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_1.1to1.3"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_1.3to1.5"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_1.3to1.5"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_1.5to1.7"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_1.5to1.7"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_1.7to1.9"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_1.7to1.9"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_1.9to2.1"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_1.9to2.1"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_2.1to2.3"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_2.1to2.3"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ){
                thad_E_dir->second["Gen_vs_Reco_THadE_2.3to2.5"].fill( perm.THad().E(), ttbar.had_top()->E() );
                thad_P_dir->second["Gen_vs_Reco_THadP_2.3to2.5"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }

        }
        /// book and fill plots for hadronic top mass corrections for lost-jet events


        /// book and fill plots after hadronic top mass corrections applied to lost-jet events
        void book_post_alpha_correction_plots( string folder ){
            // reco corrected and uncorrected plots
            book<TH1F>(folder+"/Reconstruction/Mass", "THad", "", 50, 50., 250.);
            book<TH1F>(folder+"/Reconstruction/Mass", "THad_Corrected", "", 100, 50., 450.);
            book<TH1F>(folder+"/Reconstruction/Mass", "TTbar", "", mass_bins_, 200., 2000.);
            book<TH1F>(folder+"/Reconstruction/Mass", "TTbar_Corrected", "", mass_bins_, 200., 2000.);
            book<TH2D>(folder+"/Reconstruction/Mass", "Reco_vs_Gen_TTbar", ";Gen M(ttbar); Reco M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);
            book<TH2D>(folder+"/Reconstruction/Mass", "Reco_vs_Gen_Corrected_TTbar", ";Gen M(ttbar); Reco Corrected M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);

            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Labframe", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Labframe_Corrected", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad", "", costh_bins_, -1., 1.);
            book<TH1F>(folder+"/Reconstruction/Costh", "THad_Corrected", "", costh_bins_, -1., 1.);
            book<TH2D>(folder+"/Reconstruction/Costh", "Reco_vs_Gen_Corrected_THad", ";Gen Cos(theta)(t_{h}); Reco Corrected Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);
            book<TH2D>(folder+"/Reconstruction/Costh", "Reco_vs_Gen_THad", ";Gen Cos(theta)(t_{h}); Reco Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);

            // reso corrected and uncorrected plots
            book<TH1F>(folder+"/Resolution/Mass", "THad", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "THad_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "TTbar", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "TTbar_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass", "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Mass", "Reso_MTTbar_Corrected_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_THad", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_THad_Corrected", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_TTbar", "", mass_bins_, -2., 2.);
            book<TH1F>(folder+"/Resolution/Mass", "Frac_TTbar_Corrected", "", mass_bins_, -2., 2.);

            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_2Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_2Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_3Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_3Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_4Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_vs_Gen_MTTbar_4Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_2Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_vs_Gen_MTTbar_2Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_3Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_vs_Gen_MTTbar_3Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_4Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder+"/Resolution/Parton_Acceptance", "Reso_MTTbar_Corrected_vs_Gen_MTTbar_4Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);

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
                    //values taken from 1degree vals ttp://home.fnal.gov/~jdulemba/Plots/ttbar_final_reco_3J/2018/Compare_Lost_Merged_Jets/Full/ttJetsM0/3J_Event_Plots/Final_Reco/Clear_and_MassCut_Classes/Class_Lost/Alpha_Correction/fit_parameters.json
            double alpha_E = 0.6356*( 173.1/perm.THad().M() ) + 0.2862; // only alpha_E used because it's more consistent over mtt spectrum 
            //double alpha_P = 0.3853*( 173.1/perm.THad().M() ) + 0.5502;

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
            reco_mass_dir->second["THad"].fill( perm.THad().M() );
            reco_mass_dir->second["THad_Corrected"].fill( Alpha_THad.M() );
            reco_mass_dir->second["TTbar"].fill( perm.LVect().M() );
            reco_mass_dir->second["TTbar_Corrected"].fill( (Alpha_THad + perm.TLep()).M() );
            reco_mass_dir->second["Reco_vs_Gen_TTbar"].fill( ttbar.M(), perm.LVect().M() );
            reco_mass_dir->second["Reco_vs_Gen_Corrected_TTbar"].fill( ttbar.M(), (Alpha_THad + perm.TLep()).M() );

                  // costh plots
            reco_costh_dir->second["THad_Labframe"].fill( reco_thad_labframe_cth );
            reco_costh_dir->second["THad_Labframe_Corrected"].fill( reco_corr_labframe_thad_cth );
            reco_costh_dir->second["THad"].fill( reco_thad_cth );
            reco_costh_dir->second["THad_Corrected"].fill( reco_corr_thad_cth );
            reco_costh_dir->second["Reco_vs_Gen_Corrected_THad"].fill( gen_cosths.first, reco_corr_thad_cth );
            reco_costh_dir->second["Reco_vs_Gen_THad"].fill(gen_cosths.first, reco_thad_cth);

                // reso plots
            reso_mass_dir->second["THad"].fill( ttbar.had_top()->M() - perm.THad().M() );
            reso_mass_dir->second["THad_Corrected"].fill( ttbar.had_top()->M() - Alpha_THad.M() );

            reso_mass_dir->second["TTbar"].fill( ttbar.M() - perm.LVect().M() );
            reso_mass_dir->second["Reso_MTTbar_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

            reso_mass_dir->second["TTbar_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            reso_mass_dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );

            reso_mass_dir->second["Frac_THad"].fill( (ttbar.had_top()->M() - perm.THad().M())/ttbar.had_top()->M() );
            reso_mass_dir->second["Frac_THad_Corrected"].fill( (ttbar.had_top()->M() - Alpha_THad.M())/ttbar.had_top()->M() );
            reso_mass_dir->second["Frac_TTbar"].fill( (ttbar.M() - perm.LVect().M())/ttbar.M() );
            reso_mass_dir->second["Frac_TTbar_Corrected"].fill( (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M() );


            if( ttbar.two_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_2Partons"].fill( ttbar.M() - perm.LVect().M() );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                reso_part_dir->second["Reso_MTTbar_Corrected_2Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                reso_part_dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }
            if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_3Partons"].fill( ttbar.M() - perm.LVect().M() );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                reso_part_dir->second["Reso_MTTbar_Corrected_3Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                reso_part_dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }
            if( ttbar.all_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){
                reso_part_dir->second["Reso_MTTbar_4Partons"].fill( ttbar.M() - perm.LVect().M() );
                reso_part_dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                reso_part_dir->second["Reso_MTTbar_Corrected_4Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                reso_part_dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }

            reso_costh_dir->second["THad"].fill( gen_cosths.first - reco_thad_cth );
            reso_costh_dir->second["THad_Corrected"].fill( gen_cosths.first - reco_corr_thad_cth );


            if( ttbar.M() - (Alpha_THad + perm.TLep()).M() > ttbar.M()/3. ){
                reso_dir->second["Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt() );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M() );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt() );
                reso_dir->second["Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz() );
            }
            else{
                reso_dir->second["Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt() );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M() );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt() );
                reso_dir->second["Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz() );
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


            vector<string> alpha_status = {"Pre_Alpha", "Post_Alpha"};

        //// plots for 3-jet events using event jets as perm
                // MERGED 3-jet events
            vector<string> merged_evt_type_categories = {"RIGHT", "MERGED_SWAP", "MERGED", "WRONG"};
            vector<string> merged_best_perm_solutions = {
                                                         "Clear_Classification/Class_Merged",
                                                         "Final_Classification/Class_Merged",
                                                        };


            for( auto merged_bp_solution : merged_best_perm_solutions ){

                    // plots for merged_bp from 3-jet events
                book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr" );
                book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen" );
                book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction" );
                book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution" );

                for( auto merged_evt_type_cat : merged_evt_type_categories ){
                        // plots for merged_bp from 3-jet events
                    book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr/"+merged_evt_type_cat );
                    book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/"+merged_evt_type_cat );
                    book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/"+merged_evt_type_cat );
                    book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/"+merged_evt_type_cat );

                }
            }


                // LOST 3-jet events
            vector<string> lost_evt_type_categories = {"RIGHT", "LOST_SWAP", "LOST", "WRONG"};
            vector<string> lost_best_perm_solutions = {
                                                        "Clear_Classification/Class_Lost",
                                                        "Final_Classification/Class_Lost",
                                                      };

            for( auto lost_bp_solution : lost_best_perm_solutions ){

                    // plots for lost_bp from 3-jet events
                book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr" );
                book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen" );
                book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction" );
                book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution" );
                book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction" );
                book_post_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Post_Alpha_Correction" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                        // plots for lost_bp from 3-jet events
                    book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/"+lost_evt_type_cat );
                    book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/"+lost_evt_type_cat );
                    book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/"+lost_evt_type_cat );
                    book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/"+lost_evt_type_cat );
                    book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction/"+lost_evt_type_cat );
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

        // find best perm assuming it's merged
        Permutation merged_best_perm( Permutation &perm ){

            Permutation merged_bp;
            double merged_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Merged(test_perm);
                if( test_perm.Prob() < merged_lowest_Totaldisc_3J ){
                    merged_lowest_Totaldisc_3J = test_perm.Prob();
                    merged_bp = test_perm;
                }
            }
            return merged_bp;
        }

        // find best perm assuming it's lost 
        Permutation lost_best_perm( Permutation &perm ){

            Permutation lost_bp;
            double lost_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Lost(test_perm);
                if( test_perm.Prob() < lost_lowest_Totaldisc_3J ){
                    lost_lowest_Totaldisc_3J = test_perm.Prob();
                    lost_bp = test_perm;
                }
            }
            return lost_bp;
        }



        // events with 3 jets
        std::pair< Permutation, Permutation >  process_all3J_evt(URStreamer &event)
            //std::pair< Permutation, Permutation >  process_all3J_evt(URStreamer &event, Permutation &matched_perm )
        {
            //initialize permutation objects
            //jets
            IDJet* wj1 = 0;
            IDJet* wj2 = 0;
            IDJet* bj1 = 0;
            IDJet* bj2 = 0;
            Permutation empty_perm; // perm.WJa(), WJb(), BHad(), BLep()

            vector<IDJet*> jets_vector;
            for( auto jet : object_selector_.clean_jets() ){
                jets_vector.push_back(jet);
            }
            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });

            tracker_.track("All Num BTag");

            if( !jets_vector[1]->BTagId(IDJet::BTag::CSVMEDIUM) ) return std::make_pair(empty_perm, empty_perm); // require at least 2 jets that pass btag

            tracker_.track("Num BTag >= 2");

            bj1 = jets_vector[0];
            bj2 = jets_vector[1];
            wj1 = jets_vector[2];

            // create perm from the event and find best perms from merged/lost TTBarSolver
            Permutation event_perm(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wjat, wjb, bhad, blep
            Permutation merged_bp = merged_best_perm( event_perm );
            Permutation lost_bp = lost_best_perm( event_perm );

            return std::make_pair(merged_bp, lost_bp);
        } // end of process_all3J_evt()



        void best_perm_categories( URStreamer &event ){
            permutator_.reset_3J();

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            tracker_.track("Before best perm"); 

            string best_perm_solution;

            // get best perm from event determined by likelihood disc
            std::pair< Permutation, Permutation > best_perms = process_all3J_evt(event); // < merged_bp, lost_bp >
            Permutation merged_bp = best_perms.first;
            Permutation lost_bp = best_perms.second;
            if( merged_bp.IsEmpty() && lost_bp.IsEmpty() ){ // neither best_perm has solution
                tracker_.track("best perm doesn't exist");
                best_perm_solution = "No_BP_Solutions";
                return; // skip to next event if perm is empty
            }
            else if( !merged_bp.IsEmpty() && lost_bp.IsEmpty() ){ // only merged best_perm has solution
                tracker_.track("only merged_bp exists");
                best_perm_solution = "Only_Merged";
                merged_bp_cats( ttbar, merged_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }
            else if( merged_bp.IsEmpty() && !lost_bp.IsEmpty() ){ // only lost best_perm has solution
                tracker_.track("only lost_bp exists");
                best_perm_solution = "Only_Lost";
                lost_bp_cats( ttbar, lost_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }
            else{ // both best_perms have solution
                tracker_.track("both bp exists");
                best_perm_solution = "Both_BP";
                both_bp_cats( ttbar, merged_bp, lost_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }

        } // end of best_perm_categories



        void both_bp_cats( GenTTBar &ttbar, Permutation &merged_bp, Permutation &lost_bp, string bp_solution ){

            //if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
            //    tracker_.track("Not semilep events");
            //    return;
            //}
            //tracker_.track("semilep");

            string Merged_BP_Class;
            string Lost_BP_Class;

            // values for classifications taken from slide 3 of https://indico.cern.ch/event/714185/contributions/2936743/attachments/1617897/2572311/2018_March15_Update.pdf
            // TreatMerged LDA equation: y = 0.013230*x + -2.481506
            // TreatLost LDA equation: y = -0.014565*x + 15.402863

            double TreatMerged_equation = merged_bp.Prob()-0.013230*merged_bp.THad().Pt();
            double TreatLost_equation = lost_bp.Prob()+0.014565*lost_bp.THad().Pt();

            if( TreatMerged_equation <= -2.481506 ) Merged_BP_Class = "Merged_BP/Class_Merged";
            else Merged_BP_Class = "Merged_BP/Class_Lost";
            if( TreatLost_equation <= 15.402863 ) Lost_BP_Class = "Lost_BP/Class_Lost";
            else Lost_BP_Class = "Lost_BP/Class_Merged";

            // event is classified as merged for both discriminants
            if( Merged_BP_Class == "Merged_BP/Class_Merged" && Lost_BP_Class == "Lost_BP/Class_Merged" ) merged_bp_cats( ttbar, merged_bp, "Both_BP/Class_Merged" );

            // event is classified as lost for both discriminants
            else if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Lost" ) lost_bp_cats( ttbar, lost_bp, "Both_BP/Class_Lost" );
            //if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Lost" ) lost_bp_cats( ttbar, lost_bp, "Both_BP/Class_Lost" );

            // event is classified as lost for one discriminant and merged for the other
            else{

                string merged_bp_thad_mass_cut;
                string lost_bp_thad_mass_cut;
                double thad_mass_cut;
                //if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Merged" ) thad_mass_cut = 170.;
                //else thad_mass_cut = 140.;
                if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Merged" ) thad_mass_cut = 170.;
                else thad_mass_cut = 160.;

                if( merged_bp.THad().M() > thad_mass_cut ) merged_bp_thad_mass_cut = "/THad_Mass/GCut"; // expected to have better resolution
                else merged_bp_thad_mass_cut = "/THad_Mass/LCut"; // expected to have worse resolution
                if( lost_bp.THad().M() > thad_mass_cut ) lost_bp_thad_mass_cut = "/THad_Mass/GCut"; // expected to have worse resolution
                else lost_bp_thad_mass_cut = "/THad_Mass/LCut"; // expected to have better resolution

                merged_bp_cats( ttbar, merged_bp, "Both_BP/Opposite_Class/"+Merged_BP_Class+merged_bp_thad_mass_cut );
                lost_bp_cats( ttbar, lost_bp, "Both_BP/Opposite_Class/"+Lost_BP_Class+lost_bp_thad_mass_cut );

            }

        }


        void merged_bp_cats( GenTTBar &ttbar, Permutation &merged_bp, string bp_solution ){

            string merged_perm_status;

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                tracker_.track("Not semilep events");
                //merged_perm_status = "NOTSEMILEP";
                //cout << merged_perm_status << ": " << evt_idx_ << endl;
                return;
            }
            tracker_.track("semilep");

            // get matched perm from SEMILEP events
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());


            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Merged_Event() ) merged_perm_status = "WRONG";
            if( mp.Merged_Event() ){
                if( mp.AreBsSame(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "RIGHT";
                else if( mp.AreBsFlipped(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "MERGED_SWAP";
                else merged_perm_status = "MERGED";
            }
            else merged_perm_status = "WRONG";

            // merged_bp plots

                // filling reco and reso plots in cases where the classification is obvious
            if( bp_solution == "Both_BP/Class_Merged"  || bp_solution == "Only_Merged" ){
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Merged/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Resolution", merged_bp, ttbar );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Merged/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Resolution/"+merged_perm_status, merged_bp, ttbar );
            }

                // filling reco and reso plots in for final classifications
            if( bp_solution == "Only_Merged" || bp_solution == "Both_BP/Class_Merged" ||
                        bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged/THad_Mass/GCut" ||
                        bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost/THad_Mass/GCut" ){

                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Merged/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Merged/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Merged/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Merged/Resolution", merged_bp, ttbar );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Merged/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Merged/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Merged/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Merged/Resolution/"+merged_perm_status, merged_bp, ttbar );
            }

        } // end of merged_bp_cats


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


            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Lost_Event() ) lost_perm_status = "WRONG";
            if( mp.Lost_Event() ){
                if( mp.AreBsSame(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "RIGHT";
                else if( mp.AreBsFlipped(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "LOST_SWAP";
                else lost_perm_status = "LOST";
            }
            else lost_perm_status = "WRONG";

            // lost_bp plots

                // filling reco and reso plots in cases where the classification is obvious
            if( bp_solution == "Both_BP/Class_Lost" || bp_solution == "Only_Lost" ){
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Lost/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Resolution", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Alpha_Correction", ttbar, lost_bp);
                fill_post_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Post_Alpha_Correction", ttbar, lost_bp);

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Lost/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Resolution/"+lost_perm_status, lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp);
                fill_post_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Post_Alpha_Correction/"+lost_perm_status, ttbar, lost_bp);
            }

                // filling reco and reso plots in for final classifications
            if( bp_solution == "Only_Lost" || bp_solution == "Both_BP/Class_Lost" ||
                        bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost/THad_Mass/LCut" ||
                        bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged/THad_Mass/LCut" ){

                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Lost/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Lost/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Lost/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Lost/Resolution", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Alpha_Correction", ttbar, lost_bp );
                fill_post_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Post_Alpha_Correction", ttbar, lost_bp );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Lost/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Lost/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Lost/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Lost/Resolution/"+lost_perm_status, lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp );
                fill_post_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Post_Alpha_Correction/"+lost_perm_status, ttbar, lost_bp );
            }

        } // end of lost_bp_cats


        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            //while( event.next() && evt_idx_ < 30000 )
            while( event.next() )
            {
                evt_idx_++;
                if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                tracker_.track("All Events");


                //generator selection
                bool selection = genp_selector_.select(event);
                if( !selection ){
                    Logger::log().debug() << "event has no gen selection " << endl;
                    continue;
                }
                tracker_.track("gen selection");
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                //if( ttbar.M() > 700 ) continue;

                int njets = 0;
                if( object_selector_.select(event) ) njets = object_selector_.clean_jets().size();
                if( njets < 3 ) continue;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets << endl;
                tracker_.track("njet cuts");

                /// 3 jet events
                if( njets == 3 ){ 
                    tracker_.track("njets = 3");
                    best_perm_categories(event);
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
                ("report,s", opts::value<int>()->default_value(1), "report every");
        }
};

//make it executable
int main(int argc, char *argv[])
{
    URParser &parser = URParser::instance(argc, argv);
    URDriver<ttbar_final_reco_3J> test;
    int thing = test.run();
    return thing;
}
