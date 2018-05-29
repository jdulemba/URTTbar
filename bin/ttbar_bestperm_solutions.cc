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

class ttbar_bestperm_solutions : public AnalyzerBase
{

    public:

    private:
        //counters
        unsigned long evt_idx_ = 0; //event index

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


    public:
        ttbar_bestperm_solutions(const std::string output_filename):
            AnalyzerBase("ttbar_bestperm_solutions", output_filename),
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


        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;


        //// plots for 3-jet events using event jets as perm
            vector<string> merged_evt_type_categories = {"RIGHT", "MERGED_SWAP", "MERGED", "WRONG"};
            vector<string> merged_best_perm_solutions = {
                                                         "Only_Merged", "Merged_BP", "Both_BP/Merged_BP/Class_Merged",
                                                         "Both_BP/Merged_BP/Class_Lost", "Both_BP/Class_Merged",
                                                         "Both_BP/Opposite_Class/Merged_BP/Class_Merged",
                                                         "Both_BP/Opposite_Class/Merged_BP/Class_Lost"
                                                        };


                // MERGED 3-jet events
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
                                                        "Only_Lost", "Lost_BP", "Both_BP/Lost_BP/Class_Merged",
                                                        "Both_BP/Lost_BP/Class_Lost", "Both_BP/Class_Lost",
                                                        "Both_BP/Opposite_Class/Lost_BP/Class_Merged",
                                                        "Both_BP/Opposite_Class/Lost_BP/Class_Lost"
                                                      };

            for( auto lost_bp_solution : lost_best_perm_solutions ){

                    // plots for lost_bp from 3-jet events
                book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr" );
                book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen" );
                book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction" );
                book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                        // plots for lost_bp from 3-jet events
                    book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/"+lost_evt_type_cat );
                    book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/"+lost_evt_type_cat );
                    book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/"+lost_evt_type_cat );
                    book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/"+lost_evt_type_cat );

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

            if( !merged_bp.IsEmpty() ) merged_bp_cats( ttbar, merged_bp, "Merged_BP" );
            if( !lost_bp.IsEmpty() ) lost_bp_cats( ttbar, lost_bp, "Lost_BP" );

            //cout << best_perm_solution << endl;
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

            merged_bp_cats( ttbar, merged_bp, "Both_BP/"+Merged_BP_Class );
            lost_bp_cats( ttbar, lost_bp, "Both_BP/"+Lost_BP_Class );

            // event is classified as merged for both discriminants
            if( Merged_BP_Class == "Merged_BP/Class_Merged" && Lost_BP_Class == "Lost_BP/Class_Merged" ) merged_bp_cats( ttbar, merged_bp, "Both_BP/Class_Merged" );

            // event is classified as lost for both discriminants
            else if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Lost" ) lost_bp_cats( ttbar, lost_bp, "Both_BP/Class_Lost" );

            // event is classified as lost for one discriminant and merged for the other
            else{
                merged_bp_cats( ttbar, merged_bp, "Both_BP/Opposite_Class/"+Merged_BP_Class );
                lost_bp_cats( ttbar, lost_bp, "Both_BP/Opposite_Class/"+Lost_BP_Class );

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

                    // not breaking up by perm comparison
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution", merged_bp, ttbar );

                    // not breaking up by disc val
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+merged_perm_status, merged_bp, ttbar );


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

                    // not breaking up by comparison
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution", lost_bp, ttbar );


                    // breaking up by perm comp
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+lost_perm_status, lost_bp, ttbar );

        } // end of lost_bp_cats


        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            //while(event.next() && evt_idx_ < 30000)
            while(event.next() )
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
    URDriver<ttbar_bestperm_solutions> test;
    int thing = test.run();
    return thing;
}
