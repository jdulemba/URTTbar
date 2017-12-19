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

class ttbar_reco_3J : public AnalyzerBase
{

    public:
        //enum TTNaming_Merged {RIGHT, MERGE_SWAP, MERGE, WRONG};

    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        //	    double njets_ = 0;

        double disc_cut_ = 2.;
        double M_TM_THad_Pt = -1e10;

        //histograms
        map< string, map < string, map< string, map< string, map< string, RObject > > > > > histos_merged_;
        map< string, map < string, map< string, map< string, RObject > > > > histos_disc_;
        //unordered_map<string, map< string, RObject> > histos_;
        unordered_map<string, map< string, RObject> > histos_;

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


        // TTrees for multidim analysis
        TTree *Merged_TreatMerged = 0;
//        TTree *Merged_TreatLost = 0;
//        TTree *Lost_TreatMerged = 0;
//        TTree *Lost_TreatLost = 0;

    public:
        ttbar_reco_3J(const std::string output_filename):
            AnalyzerBase("ttbar_reco_3J", output_filename),
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
        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;

            //int nbins;
            double mass_max;//, mass_min;

            int mass_bins = 450;
            int pt_bins = 200;
            int eta_bins = 200;
            int costh_bins = 100;

            double pt_min = 0.;
            double pt_max = 1000.;
            double eta_min = -2.4;
            double eta_max = 2.4;
            double costh_min = -1.;
            double costh_max = 1.;

            // mass disc
            int massdisc_bins = 120;
            double massdisc_min = -10.;
            double massdisc_max = 20.;

            // ns disc
            int nsdisc_bins = 80;
            double nsdisc_min = -5.;
            double nsdisc_max = 15.;

            // total disc
            int combdisc_bins = 160;
            double combdisc_min = -10.;
            double combdisc_max = 30.;




            if( sample == "ttJetsM0" ){
                //nbins = 180;
                ////			m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                //mass_min = 200.;
                mass_max = 2000.;
            }
            else if( sample == "ttJetsM700" ){
                //nbins = 30;
                ////			m_bins[nbins+1] = {700, 800, 900, 1000};
                //mass_min = 700.;
                mass_max = 1000.;
            }
            else if( sample == "ttJetsM1000" ){
                //nbins = 100;
                ////			m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                //mass_min = 1000.;
                mass_max = 2000.;
            }
            else{
                //nbins = 80;
                //mass_min = 200.;
                mass_max = 1000.;
            }

            string merged_expected = "Expected_Plots/MERGED";
            book<TH1F>(merged_expected, "Expected_Event_Categories_3J", "", 5, 0.5, 5.5);

            string sys_plots = "System_Plots";
            book<TH1F>(sys_plots, "TTbar_Mass", "", mass_bins, 200., 2000.);

            //attempting to make mass, pt, eta, and costheta hists for merged events
            vector<string> plot_cats = {"Merged_Plots"};
            vector<string> plot_types = {"Reconstruction", "Resolution", "Gen"};
            vector<string> evt_types = {"RIGHT", "MERGED_SWAP", "MERGED", "WRONG"};
            //vector<TTNaming_Merged> evt_types = {RIGHT, MERGE_SWAP, MERGE, WRONG};
            vector<string> disc_cats = {"LCut", "GCut"};

            for( auto plot_cat : plot_cats ){
                TDirectory* dir_plot_cat = outFile_.mkdir(plot_cat.c_str());
                dir_plot_cat->cd();
                for( auto plot_type : plot_types ){
                    TDirectory* dir_plot_type = dir_plot_cat->mkdir(plot_type.c_str());
                    dir_plot_type->cd();
                    for( auto evt_type : evt_types ){
                        TDirectory* dir_evt_type = dir_plot_type->mkdir(evt_type.c_str());
                        dir_plot_type->cd();
                        for( auto disc_cat : disc_cats ){
                            TDirectory* dir_disc = dir_evt_type->mkdir(disc_cat.c_str());
                            dir_disc->cd();

                            if( plot_type == "Reconstruction" || plot_type == "Gen" ){
                                    /// Mass
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Mass"] = RObject::book<TH1F>("THad_Mass", "", mass_bins, 100., 250.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Mass"] = RObject::book<TH1F>("TLep_Mass", "", mass_bins, 100., 250.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Mass"] = RObject::book<TH1F>("TTbar_Mass", "", mass_bins, 200., 2000.);

                                    /// Pt
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Pt"] = RObject::book<TH1F>("THad_Pt", "", pt_bins, pt_min, pt_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Pt"] = RObject::book<TH1F>("TLep_Pt", "", pt_bins, pt_min, pt_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Pt"] = RObject::book<TH1F>("TTbar_Pt", "", pt_bins, pt_min, pt_max);

                                    /// Eta
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Eta"] = RObject::book<TH1F>("THad_Eta", "", eta_bins, eta_min, eta_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Eta"] = RObject::book<TH1F>("TLep_Eta", "", eta_bins, eta_min, eta_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Eta"] = RObject::book<TH1F>("TTbar_Eta", "", eta_bins, eta_min, eta_max);

                                    /// Costh 
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Costh"] = RObject::book<TH1F>("THad_Costh", "", costh_bins, costh_min, costh_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Costh"] = RObject::book<TH1F>("TLep_Costh", "", costh_bins, costh_min, costh_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Costh"] = RObject::book<TH1F>("TTbar_Costh", "", costh_bins, costh_min, costh_max);
                            }
                            if( plot_type == "Resolution" ){
                                    /// Mass
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Mass"] = RObject::book<TH1F>("THad_Mass", "", mass_bins, -1000., 500.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Mass"] = RObject::book<TH1F>("TLep_Mass", "", mass_bins, -500., 500.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Mass"] = RObject::book<TH1F>("TTbar_Mass", "", mass_bins, -1000., mass_max);

                                    /// Pt
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Pt"] = RObject::book<TH1F>("THad_Pt", "", pt_bins, -500., 500.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Pt"] = RObject::book<TH1F>("TLep_Pt", "", pt_bins, -1000., 1000.);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Pt"] = RObject::book<TH1F>("TTbar_Pt", "", pt_bins, -1000., 500.);

                                    /// Eta
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Eta"] = RObject::book<TH1F>("THad_Eta", "", eta_bins, -2*eta_max, 2*eta_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Eta"] = RObject::book<TH1F>("TLep_Eta", "", eta_bins, -2*eta_max, 2*eta_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Eta"] = RObject::book<TH1F>("TTbar_Eta", "", eta_bins, -2*eta_max, 2*eta_max);

                                    /// Costh 
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["THad_Costh"] = RObject::book<TH1F>("THad_Costh", "", costh_bins, -2*costh_max, 2*costh_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TLep_Costh"] = RObject::book<TH1F>("TLep_Costh", "", costh_bins, -2*costh_max, 2*costh_max);
                                histos_merged_[plot_cat][plot_type][evt_type][disc_cat]["TTbar_Costh"] = RObject::book<TH1F>("TTbar_Costh", "", costh_bins, -2*costh_max, 2*costh_max);
                            }
                        }
                    }
                }
            }



            // Discriminant values for best perms
            vector<string> disc_plot_types = {"Best_Perm_Disc_Plots"};
            vector<string> event_types = {"MERGED", "LOST"};
            for( auto plot_type : disc_plot_types ){
                TDirectory* dir_plot_type = outFile_.mkdir(plot_type.c_str());
                dir_plot_type->cd();
                for( auto evt_type : event_types ){
                    TDirectory* dir_evt_type = dir_plot_type->mkdir(evt_type.c_str());
                    dir_evt_type->cd();

                    histos_disc_[plot_type][evt_type][""]["3J_bp_Massdisc"] = RObject::book<TH1F>("3J_bp_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
                    histos_disc_[plot_type][evt_type][""]["3J_bp_NSdisc"] = RObject::book<TH1F>("3J_bp_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
                    histos_disc_[plot_type][evt_type][""]["3J_bp_Totaldisc"] = RObject::book<TH1F>("3J_bp_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);

                    histos_disc_[plot_type][evt_type][""]["3J_bp_Totaldisc_vs_LeadTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_LeadTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);
                    histos_disc_[plot_type][evt_type][""]["3J_bp_Totaldisc_vs_SubleadTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_SubleadTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);
                    histos_disc_[plot_type][evt_type][""]["3J_bp_Totaldisc_vs_AverageTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_AverageTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);

                    vector<string> evt_cats = {"RIGHT", evt_type+"_SWAP", evt_type, "WRONG"};
                    for( auto evt_cat : evt_cats ){
                        TDirectory* dir_evt_cat = dir_evt_type->mkdir(evt_cat.c_str());
                        dir_evt_cat->cd();
                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_Massdisc"] = RObject::book<TH1F>("3J_bp_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_NSdisc"] = RObject::book<TH1F>("3J_bp_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_Totaldisc"] = RObject::book<TH1F>("3J_bp_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);

                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_Totaldisc_vs_LeadTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_LeadTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);
                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_Totaldisc_vs_SubleadTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_SubleadTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);
                        histos_disc_[plot_type][evt_type][evt_cat]["3J_bp_Totaldisc_vs_AverageTopPt"] = RObject::book<TH2D>("3J_bp_Totaldisc_vs_AverageTopPt", "", pt_bins, pt_min, pt_max, combdisc_bins, combdisc_min, combdisc_max);
                    }
                }
            }

            // Plots for matched perms
            vector<string> mp_plot_dirs = {"Matched_Perm_Plots"};
            vector<string> event_types_lower = {"Merged", "Lost"};
            vector<string> event_treatments = {"TreatMerged", "TreatLost"};
            for( auto mp_plot_dir : mp_plot_dirs ){ // Matched_Perm_Plots
                TDirectory* dir_mp_plot_dir = outFile_.mkdir(mp_plot_dir.c_str());
                dir_mp_plot_dir->cd();
                for( auto evt_type : event_types_lower ){ // Merged/Lost
                    TDirectory* dir_evt_type = dir_mp_plot_dir->mkdir(evt_type.c_str());
                    dir_evt_type->cd();
                    for( auto evt_treat : event_treatments ){ // Treat event as Merged or Lost
                        TDirectory* dir_evt_treat = dir_evt_type->mkdir(evt_treat.c_str());
                        dir_evt_treat->cd();

                        histos_disc_[mp_plot_dir][evt_type][evt_treat]["3J_THad_Pt"] = RObject::book<TH1F>("3J_THad_Pt", "", pt_bins, pt_min, pt_max);
                        histos_disc_[mp_plot_dir][evt_type][evt_treat]["3J_Massdisc"] = RObject::book<TH1F>("3J_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
                        histos_disc_[mp_plot_dir][evt_type][evt_treat]["3J_NSdisc"] = RObject::book<TH1F>("3J_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
                        histos_disc_[mp_plot_dir][evt_type][evt_treat]["3J_Totaldisc"] = RObject::book<TH1F>("3J_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);
                    }
                }
            }


        // test perm jets
            // merged events
            string merged_disc_tp = "Test_Perm_Disc_Plots/MERGED";

                // discriminant values
            book<TH1F>(merged_disc_tp, "3J_tp_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(merged_disc_tp, "3J_tp_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
            book<TH1F>(merged_disc_tp, "3J_tp_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);

            book<TH1F>(merged_disc_tp, "Massdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(merged_disc_tp, "NSdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(merged_disc_tp, "Totaldisc_cat", "", 3, 0.5, 3.5);

            // lost events
            string lost_disc_tp = "Test_Perm_Disc_Plots/LOST";

                // discriminant values
            book<TH1F>(lost_disc_tp, "3J_tp_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(lost_disc_tp, "3J_tp_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
            book<TH1F>(lost_disc_tp, "3J_tp_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);

            book<TH1F>(lost_disc_tp, "Massdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(lost_disc_tp, "NSdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(lost_disc_tp, "Totaldisc_cat", "", 3, 0.5, 3.5);


            // TTrees for Lost and Merged events
            Merged_TreatMerged = new TTree("Merged_TreatMerged", "Tree with merged events using merged discr.");
            Merged_TreatMerged->Branch("M_TM_THad_Pt", &M_TM_THad_Pt, "M_TM_THad_Pt/D");


            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        // events with 3 jets
        std::pair< Permutation, Permutation >  process_all3J_evt(URStreamer &event, Permutation &matched_perm )
        {
            auto merged_disc_tp_dir = histos_.find("Test_Perm_Disc_Plots/MERGED");
            auto lost_disc_tp_dir = histos_.find("Test_Perm_Disc_Plots/LOST");
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


            // initialize best_perms for merged and lost jet solvers
                // merged perm
            Permutation merged_bp;
            double merged_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
                solver_.Solve_3J_Merged(test_perm);

                merged_disc_tp_dir->second["3J_tp_Massdisc"].fill(test_perm.Merged3JDiscr());
                merged_disc_tp_dir->second["3J_tp_NSdisc"].fill(test_perm.NuDiscr());
                merged_disc_tp_dir->second["3J_tp_Totaldisc"].fill(test_perm.Prob());

                if( test_perm.Merged3JDiscr() < 10. ) merged_disc_tp_dir->second["Massdisc_cat"].fill(1); // 1 if mass disc is in range
                else merged_disc_tp_dir->second["Massdisc_cat"].fill(3); // 3 if mass disc isn't in range
                if( test_perm.NuDiscr() < 10. && test_perm.NuDiscr() > 0. ) merged_disc_tp_dir->second["NSdisc_cat"].fill(1); // 1 if NS disc is in range
                else merged_disc_tp_dir->second["NSdisc_cat"].fill(3); // 3 if NS disc isn't in range
                if( test_perm.Prob() < 15. ) merged_disc_tp_dir->second["Totaldisc_cat"].fill(1); // 1 if total disc is in range
                else merged_disc_tp_dir->second["Totaldisc_cat"].fill(3); // 3 if total disc isn't in range

                if( test_perm.Prob() < merged_lowest_Totaldisc_3J ){
                    merged_lowest_Totaldisc_3J = test_perm.Prob();
                    merged_bp = test_perm;
                }
            }

                // lost perm
            Permutation lost_bp;
            double lost_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
                solver_.Solve_3J_Lost(test_perm);

                lost_disc_tp_dir->second["3J_tp_Massdisc"].fill(test_perm.Lost3JDiscr());
                lost_disc_tp_dir->second["3J_tp_NSdisc"].fill(test_perm.NuDiscr());
                lost_disc_tp_dir->second["3J_tp_Totaldisc"].fill(test_perm.Prob());

                //if( test_perm.Lost3JDiscr() < 10. ) lost_disc_tp_dir->second["Massdisc_cat"].fill(1); // 1 if mass disc is in range
                //else lost_disc_tp_dir->second["Massdisc_cat"].fill(3); // 3 if mass disc isn't in range
                //if( test_perm.NuDiscr() < 10. && test_perm.NuDiscr() > 0. ) lost_disc_tp_dir->second["NSdisc_cat"].fill(1); // 1 if NS disc is in range
                //else lost_disc_tp_dir->second["NSdisc_cat"].fill(3); // 3 if NS disc isn't in range
                //if( test_perm.Prob() < 15. ) lost_disc_tp_dir->second["Totaldisc_cat"].fill(1); // 1 if total disc is in range
                //else lost_disc_tp_dir->second["Totaldisc_cat"].fill(3); // 3 if total disc isn't in range

                if( test_perm.Prob() < lost_lowest_Totaldisc_3J ){
                    lost_lowest_Totaldisc_3J = test_perm.Prob();
                    lost_bp = test_perm;
                }
            }

            return std::make_pair(merged_bp, lost_bp);
        } // end of process_3J_evt()


        void matched_perm_info( URStreamer &event ){
            permutator_.reset_3J();

            auto mp_dir = histos_disc_.find("Matched_Perm_Plots")->second;

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) return; // skip to next event if it's not semilep

            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());

            if( mp.IsEmpty() || !(mp.BHad() && mp.BLep()) || !(mp.WJa() || mp.WJb()) ) return; //make sure matched_perm has BHad, BLep, and at least one WJet
            if( mp.unique_matches() < 3 ) return; // require perm to have at least 3 distinct objects
            if( !(mp.Merged_Event() || mp.Lost_Event()) ) return;// require event to be merged or lost
            if( mp.Merged_BHadBLep() ) return; // event must have separate bhad and blep

            string event_perm_type;
            //string event_perm_treatement;

            if( mp.Merged_Event() ) event_perm_type = "Merged";
            if( mp.Lost_Event() ) event_perm_type = "Lost";

            M_TM_THad_Pt = -1e10;

                // merged perm
            Permutation merged_bp;
            double merged_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(mp.WJa(), mp.WJb(), mp.BHad(), mp.BLep(), mp.L(), mp.MET(), mp.LepCharge()) ){
                solver_.Solve_3J_Merged(test_perm);

                if( test_perm.Prob() < merged_lowest_Totaldisc_3J ){
                    merged_lowest_Totaldisc_3J = test_perm.Prob();
                    merged_bp = test_perm;
                }
            }

                // lost perm
            Permutation lost_bp;
            double lost_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(mp.WJa(), mp.WJb(), mp.BHad(), mp.BLep(), mp.L(), mp.MET(), mp.LepCharge()) ){
                solver_.Solve_3J_Lost(test_perm);

                if( test_perm.Prob() < lost_lowest_Totaldisc_3J ){
                    lost_lowest_Totaldisc_3J = test_perm.Prob();
                    lost_bp = test_perm;
                }
            }

                // initialize THad objects for merged and lost best perms
            TLorentzVector merged_bp_THad;
            TLorentzVector lost_bp_THad;
            if( merged_bp.Merged_BHadWJa() || merged_bp.Merged_BLepWJa() || merged_bp.Merged_WJets() || merged_bp.Lost_WJa() ){
                merged_bp_THad = *merged_bp.BHad()+*merged_bp.WJb();
            }
            if( merged_bp.Merged_BHadWJb() || merged_bp.Merged_BLepWJb() || merged_bp.Lost_WJb() ){
                merged_bp_THad = *merged_bp.BHad()+*merged_bp.WJa();
            }

            if( lost_bp.Merged_BHadWJa() || lost_bp.Merged_BLepWJa() || lost_bp.Merged_WJets() || lost_bp.Lost_WJa() ){
                lost_bp_THad = *lost_bp.BHad()+*lost_bp.WJb();
            }
            if( lost_bp.Merged_BHadWJb() || lost_bp.Merged_BLepWJb() || lost_bp.Lost_WJb() ){
                lost_bp_THad = *lost_bp.BHad()+*lost_bp.WJa();
            }

            // fill plots
            mp_dir[event_perm_type]["TreatMerged"]["3J_THad_Pt"].fill(merged_bp_THad.Pt());
            mp_dir[event_perm_type]["TreatMerged"]["3J_Massdisc"].fill(merged_bp.Merged3JDiscr());
            mp_dir[event_perm_type]["TreatMerged"]["3J_NSdisc"].fill(merged_bp.NuDiscr());
            mp_dir[event_perm_type]["TreatMerged"]["3J_Totaldisc"].fill(merged_bp.Prob());

            mp_dir[event_perm_type]["TreatLost"]["3J_THad_Pt"].fill(lost_bp_THad.Pt());
            mp_dir[event_perm_type]["TreatLost"]["3J_Massdisc"].fill(lost_bp.Lost3JDiscr());
            mp_dir[event_perm_type]["TreatLost"]["3J_NSdisc"].fill(lost_bp.NuDiscr());
            mp_dir[event_perm_type]["TreatLost"]["3J_Totaldisc"].fill(lost_bp.Prob());

            if( mp.Merged_Event() ){
                M_TM_THad_Pt = merged_bp_THad.Pt();
                //cout << "Merged treat merged thad pt: " << M_TM_THad_Pt << endl;
                if( M_TM_THad_Pt > 0 ) Merged_TreatMerged->Fill();
            }

        }// end of matched_perm_info()

        //void merged_bp_cats( URStreamer &event ){
        //    permutator_.reset_3J();

        //    // subdirectories
        //    auto merged_exp_dir = histos_.find("Expected_Plots/MERGED");
        //    auto merged_plots_dir = histos_merged_.find("Merged_Plots")->second;
        //    auto merged_bp_disc_dir = histos_disc_.find("Best_Perm_Disc_Plots")->second;

        //    //generator selection
        //    bool selection = genp_selector_.select(event);
        //    if( !selection ){
        //        Logger::log().debug() << "event has no gen selection " << endl;
        //        return;
        //    }
        //    GenTTBar &ttbar = genp_selector_.ttbar_system();

        //    tracker_.track("Before best perm"); 

        //    // get best perm from event determined by likelihood disc
        //    std::pair< Permutation, Permutation > best_perms = process_all3J_evt(event); // < merged_bp, lost_bp >
        //    Permutation merged_bp = best_perms.first;
        //    if( merged_bp.IsEmpty() ){
        //        tracker_.track("best perm doesn't exist");
        //        return; // skip to next event if perm is empty
        //    }
        //    tracker_.track("best perm exists");

        //    TLorentzVector leadingtop = ( merged_bp.THad().Pt() > merged_bp.TLep().Pt() ) ? merged_bp.THad() : merged_bp.TLep();
        //    TLorentzVector subleadingtop = ( merged_bp.THad().Pt() > merged_bp.TLep().Pt() ) ? merged_bp.TLep() : merged_bp.THad();

        //    string merged_perm_status;
        //    string merged_perm_disc_cut_status;
        //    if( merged_bp.Prob() < disc_cut_ ){
        //        merged_perm_disc_cut_status = "LCut";
        //    }
        //    else{
        //        merged_perm_disc_cut_status = "GCut";
        //    }

        //        // fill discriminant plots
        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_Massdisc"].fill(merged_bp.Merged3JDiscr());
        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_NSdisc"].fill(merged_bp.NuDiscr());
        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_Totaldisc"].fill(merged_bp.Prob());

        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), merged_bp.Prob());
        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), merged_bp.Prob());
        //    merged_bp_disc_dir["MERGED"][""]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), merged_bp.Prob());

        //        // variables for costh*
        //    hyp::TTbar gen_ttang(ttbar);
        //    auto gen_ttcm = gen_ttang.to_CM();
        //    double gen_ttbar_cth = gen_ttang.unit3D().Dot(gen_ttcm.unit3D());

        //    hyp::TTbar reco_ttang(merged_bp); // doesn't work because not all wjets defined
        //    auto reco_ttcm = reco_ttang.to_CM();
        //    double reco_ttbar_cth = reco_ttang.unit3D().Dot(reco_ttcm.unit3D());
        //    double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
        //    double reco_tlep_cth = reco_ttang.unit3D().Dot(reco_ttcm.tlep().unit3D());

        //    merged_exp_dir->second["Expected_Event_Categories_3J"].fill(1);// tot expected events == 1

        //    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty

        //        merged_exp_dir->second["Expected_Event_Categories_3J"].fill(5);// expected other events == 5
        //        //Logger::log().debug() << "not semilep event " << evt_idx_ << endl;

        //        // Reco Plots
        //        merged_plots_dir["Reconstruction"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Mass"].fill(merged_bp.LVect().M());
        //        merged_plots_dir["Reconstruction"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Pt"].fill(merged_bp.LVect().Pt());
        //        merged_plots_dir["Reconstruction"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Eta"].fill(merged_bp.LVect().Eta());
        //        merged_plots_dir["Reconstruction"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Costh"].fill(reco_ttbar_cth);
        //        // Resolution Plots
        //        merged_plots_dir["Resolution"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Mass"].fill(ttbar.M() - merged_bp.LVect().M());
        //        merged_plots_dir["Resolution"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Pt"].fill(ttbar.Pt() - merged_bp.LVect().Pt());
        //        merged_plots_dir["Resolution"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Eta"].fill(ttbar.Eta() - merged_bp.LVect().Eta());
        //        merged_plots_dir["Resolution"]["WRONG"][merged_perm_disc_cut_status]["TTbar_Costh"].fill(gen_ttbar_cth - reco_ttbar_cth);

        //        // Discriminant Plots 
        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_Massdisc"].fill(merged_bp.Merged3JDiscr());
        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_NSdisc"].fill(merged_bp.NuDiscr());
        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_Totaldisc"].fill(merged_bp.Prob());

        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), merged_bp.Prob());
        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), merged_bp.Prob());
        //        merged_bp_disc_dir["MERGED"]["WRONG"]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), merged_bp.Prob());

        //        tracker_.track("Not semilep events");

        //        return;
        //    }
        //    tracker_.track("semilep");

        //    if( ttbar.partial_hadronic_merged(0.4) ){ // gen partons merged
        //        merged_exp_dir->second["Expected_Event_Categories_3J"].fill(3);// expected merged events == 3
        //    }
        //    else{ // gen partons not merged
        //        merged_exp_dir->second["Expected_Event_Categories_3J"].fill(5);// expected other events == 5
        //    }

        //    double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());
        //    double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());

        //    // get matched perm from event
        //    Permutation mp = dr_matcher_.dr_match(
        //            genp_selector_.ttbar_final_system(),
        //            object_selector_.clean_jets(),
        //            object_selector_.lepton(),
        //            object_selector_.met(),
        //            object_selector_.lepton_charge());


        //    if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Merged_Event() ) merged_perm_status = "WRONG";
        //    if( mp.Merged_Event() ){
        //        if( mp.AreBsSame(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "RIGHT";
        //        else if( mp.AreBsFlipped(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "MERGED_SWAP";
        //        else merged_perm_status = "MERGED";
        //    }
        //    else merged_perm_status = "WRONG";

        //    // Gen Plots
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Mass"].fill(ttbar.M());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Pt"].fill(ttbar.Pt());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Eta"].fill(ttbar.Eta());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Costh"].fill(gen_ttbar_cth);

        //        // Mass
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["THad_Mass"].fill(ttbar.had_top()->M());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Mass"].fill(ttbar.lep_top()->M());
        //        // Pt
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["THad_Pt"].fill(ttbar.had_top()->Pt());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Pt"].fill(ttbar.lep_top()->Pt());
        //        // Eta
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["THad_Eta"].fill(ttbar.had_top()->Eta());
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Eta"].fill(ttbar.lep_top()->Eta());
        //        // Costh
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["THad_Costh"].fill(gen_thad_cth);
        //    merged_plots_dir["Gen"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Costh"].fill(gen_tlep_cth);
        //    
        //    // Reco Plots
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Mass"].fill(merged_bp.LVect().M());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Pt"].fill(merged_bp.LVect().Pt());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Eta"].fill(merged_bp.LVect().Eta());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Costh"].fill(reco_ttbar_cth);

        //        // Mass
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["THad_Mass"].fill(merged_bp.THad().M());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Mass"].fill(merged_bp.TLep().M());
        //        // Pt
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["THad_Pt"].fill(merged_bp.THad().Pt());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Pt"].fill(merged_bp.TLep().Pt());
        //        // Eta
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["THad_Eta"].fill(merged_bp.THad().Eta());
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Eta"].fill(merged_bp.TLep().Eta());
        //        // Costh
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["THad_Costh"].fill(reco_thad_cth);
        //    merged_plots_dir["Reconstruction"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Costh"].fill(reco_tlep_cth);
        //    
        //    // Resolution Plots
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Mass"].fill(ttbar.M() - merged_bp.LVect().M());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Pt"].fill(ttbar.Pt() - merged_bp.LVect().Pt());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Eta"].fill(ttbar.Eta() - merged_bp.LVect().Eta());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TTbar_Costh"].fill(gen_ttbar_cth - reco_ttbar_cth);

        //        // Mass
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["THad_Mass"].fill(ttbar.had_top()->M() - merged_bp.THad().M());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Mass"].fill(ttbar.lep_top()->M() - merged_bp.TLep().M());
        //        // Pt
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["THad_Pt"].fill(ttbar.had_top()->Pt() - merged_bp.THad().Pt());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Pt"].fill(ttbar.lep_top()->Pt() - merged_bp.TLep().Pt());
        //        // Eta
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["THad_Eta"].fill(ttbar.had_top()->Eta() - merged_bp.THad().Eta());
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Eta"].fill(ttbar.lep_top()->Eta() - merged_bp.TLep().Eta());
        //        // Costh
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["THad_Costh"].fill(gen_thad_cth - reco_thad_cth);
        //    merged_plots_dir["Resolution"][merged_perm_status][merged_perm_disc_cut_status]["TLep_Costh"].fill(gen_tlep_cth - reco_tlep_cth);
        //   
        //    // Discriminant Plots 
        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_Massdisc"].fill(merged_bp.Merged3JDiscr());
        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_NSdisc"].fill(merged_bp.NuDiscr());
        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_Totaldisc"].fill(merged_bp.Prob());

        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), merged_bp.Prob());
        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), merged_bp.Prob());
        //    merged_bp_disc_dir["MERGED"][merged_perm_status]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), merged_bp.Prob());

        //} // end of merged_bp_cats


        //void lost_bp_cats( URStreamer &event ){
        //    permutator_.reset_3J();

        //    // subdirectories
        //    auto lost_bp_disc_dir = histos_disc_.find("Best_Perm_Disc_Plots")->second;

        //    //generator selection
        //    bool selection = genp_selector_.select(event);
        //    if( !selection ){
        //        Logger::log().debug() << "event has no gen selection " << endl;
        //        return;
        //    }
        //    GenTTBar &ttbar = genp_selector_.ttbar_system();

        //    tracker_.track("Before lost best perm"); 

        //    // get best perm from event determined by likelihood disc
        //    std::pair< Permutation, Permutation > best_perms = process_all3J_evt(event); // < merged_bp, lost_bp >
        //    Permutation lost_bp = best_perms.second;
        //    if( lost_bp.IsEmpty() ){
        //        tracker_.track("lost best perm doesn't exist");
        //        return; // skip to next event if perm is empty
        //    }

        //    tracker_.track("lost best perm exists");

        //    TLorentzVector leadingtop = ( lost_bp.THad().Pt() > lost_bp.TLep().Pt() ) ? lost_bp.THad() : lost_bp.TLep();
        //    TLorentzVector subleadingtop = ( lost_bp.THad().Pt() > lost_bp.TLep().Pt() ) ? lost_bp.TLep() : lost_bp.THad();

        //    // Discriminant Plots 
        //    lost_bp_disc_dir["LOST"][""]["3J_bp_Massdisc"].fill(lost_bp.Lost3JDiscr());
        //    lost_bp_disc_dir["LOST"][""]["3J_bp_NSdisc"].fill(lost_bp.NuDiscr());
        //    lost_bp_disc_dir["LOST"][""]["3J_bp_Totaldisc"].fill(lost_bp.Prob());

        //    lost_bp_disc_dir["LOST"][""]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), lost_bp.Prob());
        //    lost_bp_disc_dir["LOST"][""]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), lost_bp.Prob());
        //    lost_bp_disc_dir["LOST"][""]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), lost_bp.Prob());


        //    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
        //        // Discriminant Plots 
        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_Massdisc"].fill(lost_bp.Lost3JDiscr());
        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_NSdisc"].fill(lost_bp.NuDiscr());
        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_Totaldisc"].fill(lost_bp.Prob());

        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), lost_bp.Prob());
        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), lost_bp.Prob());
        //        lost_bp_disc_dir["LOST"]["WRONG"]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), lost_bp.Prob());

        //        tracker_.track("Not semilep events");
        //        return;
        //    }

        //    // get matched perm from event
        //    Permutation mp = dr_matcher_.dr_match(
        //            genp_selector_.ttbar_final_system(),
        //            object_selector_.clean_jets(),
        //            object_selector_.lepton(),
        //            object_selector_.met(),
        //            object_selector_.lepton_charge());

        //    string lost_perm_status;

        //    if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Lost_Event() ) lost_perm_status = "WRONG";
        //    if( mp.Lost_Event() ){
        //        if( mp.AreBsSame(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "RIGHT";
        //        else if( mp.AreBsFlipped(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "LOST_SWAP";
        //        else lost_perm_status = "LOST";
        //    }
        //    else lost_perm_status = "WRONG";

        //    // Discriminant Plots 
        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_Massdisc"].fill(lost_bp.Lost3JDiscr());
        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_NSdisc"].fill(lost_bp.NuDiscr());
        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_Totaldisc"].fill(lost_bp.Prob());

        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_Totaldisc_vs_LeadTopPt"].fill(leadingtop.Pt(), lost_bp.Prob());
        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_Totaldisc_vs_SubleadTopPt"].fill(subleadingtop.Pt(), lost_bp.Prob());
        //    lost_bp_disc_dir["LOST"][lost_perm_status]["3J_bp_Totaldisc_vs_AverageTopPt"].fill(0.5*(leadingtop.Pt()+subleadingtop.Pt()), lost_bp.Prob());

        //}// end of lost_bp_cats


        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            while(event.next() /*&& evt_idx_ < 30000*/)
            {
                evt_idx_++;
                if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                tracker_.track("All Events");

                auto sys_dir = histos_.find("System_Plots")->second;

                //generator selection
                bool selection = genp_selector_.select(event);
                if( !selection ){
                    Logger::log().debug() << "event has no gen selection " << endl;
                    continue;
                }
                tracker_.track("gen selection");
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                sys_dir["TTbar_Mass"].fill( ttbar.M() );
                //if( ttbar.M() > 700 ) continue;

                int njets = 0;
                if( object_selector_.select(event) ) njets = object_selector_.clean_jets().size();
                if( njets < 3 ) continue;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets << endl;
                tracker_.track("njet cuts");

                /// 3 jet events
                if( njets == 3 ){ 
                    tracker_.track("njets = 3");
                    //merged_bp_cats(event);
                    //lost_bp_cats(event);
                    matched_perm_info(event);
                }

                //Merged_TreatMerged->Fill();

            } // end of event loop

            Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
        } // end of analyze()

        //this method is called at the end of the job, by default saves
        //every histogram/tree produced, override it if you need something more
        virtual void end()
        {
            //Merged_TreatMerged->Write();
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
    URDriver<ttbar_reco_3J> test;
    int thing = test.run();
    return thing;
}
