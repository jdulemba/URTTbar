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

#include "TROOT.h"
#include "TTree.h"
//#include <string>

using namespace TMath;
//using namespace std;

class jet_perm_disc : public AnalyzerBase
{
    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        double njets_ = 0;

//        const char *DRnames_[2] = {"DRP4", "DRP8"};
//        double DR_[2] = {0.4, 0.8};
        //            const char *DRnames_[4] = {"DRP4", "DRP5", "DRP6", "DRP8"};
        //	    double DR_[4] = {0.4, 0.5, 0.6, 0.8};
        //            const char *DRnames_[1] = {"DRP4"};
        //	    double DR_[1] = {0.4};

//        double disc_cut_ = 2.;


        //histograms
        unordered_map<string, map< string, RObject> > histos_;
        //	    map<TTNaming, string> > naming_;

        CutFlowTracker tracker_; //tracks how many events pass cuts

        // initialize bin ranges
        int nbins;
        double mass_max_, mass_min_;

        int mass_bins_ = 500;
        int pt_bins_ = 200;
        int DR_bins_ = 200;
        int eta_bins_ = 200;
        int costh_bins_ = 100;

        double DR_min_ = 0.;
        double DR_max_ = 8.;
        double pt_min_ = 0.;
        double pt_max_ = 2000.;
        double eta_min_ = -2.4;
        double eta_max_ = 2.4;
        double costh_min_ = -1.;
        double costh_max_ = 1.;

        // initialize variables for m_3j vs m_2j dist

            // m_3j vars
        double Mass_BHadWJaWJb_; // correct combo
        double Mass_BLepWJaWJb_; // wrong combo

            // m_2j vars
        double Mass_WJaWJb_; // correct combo
        double Mass_BHadWJa_, Mass_BHadWJb_, Mass_BLepWJa_, Mass_BLepWJb_; // wrong combo
        


        //values for ttree

        // 3 Jets
        double Unmerged_BLep_mass_3J = -1.;
        double Unmerged_BLep_comb_mass_3J = -1.;
        double Max_Unmerged_BLep_mass_3J = -1.;

        double Merged_BLep_mass_3J = -1.;
        double Merged_BLep_comb_mass_3J = -1.;
        double Max_Merged_BLep_mass_3J = -1.;

        double Unmerged_BHad_mass_3J = -1.;
        double Unmerged_BHad_comb_mass_3J = -1.;
        double Max_Unmerged_BHad_mass_3J = -1.;

        double Merged_BHad_mass_3J = -1.;
        double Merged_BHad_comb_mass_3J = -1.;
        double Max_Merged_BHad_mass_3J = -1.;

        double Correct_merged_wjet_mass_3J = -1.;
        double Correct_merged_wjet_comb_mass_3J = -1.;
        double Max_Correct_merged_wjet_mass_3J = -1.;

        double Wrong_merged_wjet_mass_3J = -1.;
        double Wrong_merged_wjet_comb_mass_3J = -1.;
        double Max_Wrong_merged_wjet_mass_3J = -1.;




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
        jet_perm_disc(const std::string output_filename):
            AnalyzerBase("jet_perm_disc", output_filename),
            tracker_(),
            object_selector_(),
            genp_selector_(TTGenParticleSelector::SelMode::LHE),
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
        Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
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

            if( sample == "ttJetsM0" ){
                nbins = 22;
                //			m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min_ = 250.;
                mass_max_ = 2000.;
            }
            if( sample == "ttJetsM700" ){
                nbins = 10;
                //			m_bins[nbins+1] = {700, 800, 900, 1000};
                mass_min_ = 700.;
                mass_max_ = 1000.;
            }
            if( sample == "ttJetsM1000" ){
                nbins = 10;
                //			m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min_ = 1000.;
                mass_max_ = 2000.;
            }

//            // perm disc hists
//            int massdisc_bins = 20;
//            double massdisc_min = -10.;
//            double massdisc_max = 10.;
//
//            // ns disc hists
//            int ns_disc_bins = 20;
//            double ns_disc_min = 0.;
//            double ns_disc_max = 10.;
//
//            // ns chi hists
//            int nschi_bins = 100;
//            double nschi_min = 0.;
//            double nschi_max = 1000.;
//
//            // combined disc hists
//            int comb_disc_bins = 20;
//            double comb_disc_min = -5.;
//            double comb_disc_max = 15.;

            // Gen Mttbar Plots
            string gen_plots = "Gen_Plots";
            book<TH1F>(gen_plots, "Mttbar", "", 18, 200., 2000.);


            // Resolved 4+ jets plots
            string disc_correct_4PJ = "Correct_4PJ_Disc_Plots";
            string disc_wrong_4PJ = "Wrong_4PJ_Disc_Plots";

            book<TH2D>(disc_correct_4PJ, "4PJ_M3j_vs_M2j_correct", "", 100, 0., 1000., 100, 0, 1000.); 
            book<TH2D>(disc_wrong_4PJ, "4PJ_M3j_vs_M2j_wrong", "", 100, 0., 1000., 100, 0, 1000.); 

            // Merged jets discriminant plots
            string disc_right = "Correct_Disc_Plots";
            string disc_wrong = "Wrong_Disc_Plots";

            // 3 jets
            book<TH2D>(disc_right, "3J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
            book<TH2D>(disc_wrong, "3J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}


            // Neutrino Solver discriminant plots
            book<TH1F>(disc_right, "3J_nusolver_chi2_correct", "", 1000, 0., 1000.); // ns chi2 val hist for correct perms
            book<TH1F>(disc_wrong, "3J_nusolver_chi2_wrong", "", 1000, 0., 1000.); // ns chi2 val hist for wrong perms



        Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }


        void matched_4PJ_perms(URStreamer &event, Permutation &mp){// perms have all jets and are unmerged

            auto correct_4PJ_disc_dir = histos_.find("Correct_4PJ_Disc_Plots");
            auto wrong_4PJ_disc_dir = histos_.find("Wrong_4PJ_Disc_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();

//            cout << "M3j: " << (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() << endl;
//            cout << "M2j: " << (*mp.WJa()+*mp.WJb()).M() << endl;

            // correct jet combinations
            correct_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_correct"].fill( (*mp.WJa()+*mp.WJb()).M(), (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() );

            // wrong jet combos
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BHad()+*mp.WJa()).M(), (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() );
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BHad()+*mp.WJb()).M(), (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() );

            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BLep()+*mp.WJa()).M(), (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() );
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BLep()+*mp.WJb()).M(), (*mp.BHad()+*mp.WJa()+*mp.WJb()).M() );
            
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BHad()+*mp.WJa()).M(), (*mp.BLep()+*mp.WJa()+*mp.WJb()).M() );
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BHad()+*mp.WJb()).M(), (*mp.BLep()+*mp.WJa()+*mp.WJb()).M() );

            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BLep()+*mp.WJa()).M(), (*mp.BLep()+*mp.WJa()+*mp.WJb()).M() );
            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.BLep()+*mp.WJb()).M(), (*mp.BLep()+*mp.WJa()+*mp.WJb()).M() );

            wrong_4PJ_disc_dir->second["4PJ_M3j_vs_M2j_wrong"].fill( (*mp.WJa()+*mp.WJb()).M(), (*mp.BLep()+*mp.WJa()+*mp.WJb()).M() );

        }


        // right/wrong matched perm plots
        void matched_perm_3J_plots(URStreamer &event, Permutation &matched_perm){

            auto disc_wrong_dir = histos_.find("Wrong_Disc_Plots");
            auto disc_correct_dir = histos_.find("Correct_Disc_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();

//            GenTop* thad = ttbar.had_top();


            solver_.Solve_3J_Merged(matched_perm);


            if( !matched_perm.Merged_Event() ){ // no jets merged
                if( matched_perm.WJa() ){ // if wja exists
                    // blep paired with wja
                    Unmerged_BLep_mass_3J = matched_perm.BLep()->M();
                    Unmerged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                    Max_Unmerged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                    //                      t1->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // bhad paired with wja
                    Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                    Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                    Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
                    //                t3->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                }

                if( matched_perm.WJb() ){ // if wjb exists
                    // blep paired with wjb
                    Unmerged_BLep_mass_3J = matched_perm.BLep()->M();
                    Unmerged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJb()).M();
                    Max_Unmerged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJb()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJb()->M();
                    //              t1->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // bhad paired with wjb
                    Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                    Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                    Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
                    //                t3->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                }
            } // end of unmerged events

            // 2 jets merged together
            //wrong perm
            if( matched_perm.Merged_BLepWJa() && matched_perm.WJb() ){ // only blep and wja merged and wjb exists
                Merged_BLep_mass_3J = matched_perm.BLep()->M(); 
                Merged_BLep_comb_mass_3J =  (*matched_perm.BLep()+*matched_perm.WJb()).M();
                Max_Merged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJb()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJb()->M();
                //            t2->Fill();
                disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);
                disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

            }

            //wrong perm
            if( matched_perm.Merged_BLepWJb() && matched_perm.WJa() ){ // only blep and wjb merged and wja exists
                Merged_BLep_mass_3J = matched_perm.BLep()->M(); 
                Merged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                Max_Merged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                //            t2->Fill();
                disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);
                disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

            }

            // correct perm
            if( matched_perm.Merged_BHadWJa() && matched_perm.WJb() ){ // only bhad and wja merged and wjb exists
                Merged_BHad_mass_3J = matched_perm.BHad()->M();
                Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
                //            t4->Fill();
                disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

            }

            // correct perm
            if( matched_perm.Merged_BHadWJb() && matched_perm.WJa() ){ // only bhad and wjb merged and wja exists
                Merged_BHad_mass_3J =  matched_perm.BHad()->M();
                Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
                //            t4->Fill();
                disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

            }

            // merged wjets
            if( matched_perm.Merged_WJets() ){ // only wjets merged
                // correct perm
                Correct_merged_wjet_mass_3J = matched_perm.BHad()->M();
                Correct_merged_wjet_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                Max_Correct_merged_wjet_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
                //			t5->Fill();
                disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Correct_merged_wjet_mass_3J, Correct_merged_wjet_comb_mass_3J);
                disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());


                //wrong perm
                Wrong_merged_wjet_mass_3J = matched_perm.BLep()->M();
                Wrong_merged_wjet_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                Max_Wrong_merged_wjet_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                //			t6->Fill();
                disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Wrong_merged_wjet_mass_3J, Wrong_merged_wjet_comb_mass_3J);
                //            disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());


            } // end of merged W jets loop

            permutator_.reset_3J();

        } // end of matched_perm_3J_plots()




        // events with 4+ jets
        Permutation process_4PJ_evt(URStreamer &event){
            //auto like_dir = histos_.find("Likelihood_Plots");
            //auto test_dir = histos_.find("Test_Perm_Plots");

            //initialize permutation objects
            //jets
            IDJet* wj1 = 0;
            IDJet* wj2 = 0;
            IDJet* bj1 = 0;
            IDJet* bj2 = 0;
            Permutation empty_perm; // perm.WJa(), WJb(), BHad(), BLep()

            vector<IDJet*> jets_vector;
            int num_btag = 0;

            for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
                jets_vector.push_back(*jets);

                if( (*jets)->BTagId(cut_medium_b_) ) ++num_btag;
            }

            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->btagCSVV2() > B->btagCSVV2() ); }); // CSVv2 btagger

            if( num_btag < 2 ) return empty_perm; // require at least 2 b-tagged jets

            bj1 = jets_vector[0];
            bj2 = jets_vector[1];
            wj1 = jets_vector[2];
            wj2 = jets_vector[3];

//            // test perm plots
//
//            // mass
//            test_dir->second["3J_bj1wj1_Mass"].fill( (*bj1+*wj1).M() );
//            test_dir->second["3J_bj2wj1_Mass"].fill( (*bj2+*wj1).M() );
//            test_dir->second["3J_bj1_Mass"].fill( bj1->M() );
//            test_dir->second["3J_bj2_Mass"].fill( bj2->M() );
//            test_dir->second["3J_wj1_Mass"].fill( wj1->M() );
//
//            // MVA
//            test_dir->second["3J_bj1_MVA"].fill( bj1->btagCMVA() );
//            test_dir->second["3J_bj2_MVA"].fill( bj2->btagCMVA() );
//            test_dir->second["3J_wj1_MVA"].fill( wj1->btagCMVA() );
//
//            // CvsL
//            test_dir->second["3J_bj1_CvsL"].fill( bj1->CvsLtag() );
//            test_dir->second["3J_bj2_CvsL"].fill( bj2->CvsLtag() );
//            test_dir->second["3J_wj1_CvsL"].fill( wj1->CvsLtag() );
//
//            //Mult 
//            test_dir->second["3J_bj1_Mult"].fill( bj1->nConstituents() );
//            test_dir->second["3J_bj2_Mult"].fill( bj2->nConstituents() );
//            test_dir->second["3J_wj1_Mult"].fill( wj1->nConstituents() );
//
//
//
//
            Permutation best_perm(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());
//            double lowest_massdisc_3J = 1e10;
//            double lowest_NSchi_3J = 1e10;
//            double lowest_NSdisc_3J = 1e10;
//            double lowest_Totaldisc_3J = 1e10;
//
//            for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
//                solver_.Solve_3J_Merged(test_perm);
//
//                if( test_perm.Merged3JDiscr() < lowest_massdisc_3J ){
//                    lowest_massdisc_3J = test_perm.Merged3JDiscr();
//                    //                best_perm = test_perm;
//                }
//                like_dir->second["3J_Event_All_Massdisc"].fill(test_perm.Merged3JDiscr()); // all perm disc values
//
//                if( test_perm.NuChisq() < lowest_NSchi_3J && test_perm.NuChisq() > 0. ){
//                    lowest_NSchi_3J = test_perm.NuChisq();
//                    //                cout << "NS chi: " << lowest_NSchi_3J << endl;
//                }
//                like_dir->second["3J_Event_All_NSchi"].fill(test_perm.NuChisq()); // all neutrino solver chi2 values
//
//                if( test_perm.NuDiscr() < lowest_NSdisc_3J ){
//                    lowest_NSdisc_3J = test_perm.NuDiscr();
//                    //                cout << "NS disc: " << lowest_NSdisc_3J << endl;
//                }
//                like_dir->second["3J_Event_All_NSdisc"].fill(test_perm.NuDiscr()); // all NS disc values
//
//                if( test_perm.Prob() < lowest_Totaldisc_3J ){
//                    lowest_Totaldisc_3J = test_perm.Prob();
//                    best_perm = test_perm;
//                    //                cout << "NS disc: " << lowest_NSdisc_3J << endl;
//                }
//                like_dir->second["3J_Event_All_Totaldisc"].fill(test_perm.Prob()); // all total (perm + NS) disc values
//
//            }
//
//            // lowest values from event
//            like_dir->second["3J_Event_Lowest_Massdisc"].fill(lowest_massdisc_3J); // lowest perm disc value
//            like_dir->second["3J_Event_Lowest_NSchi"].fill(lowest_NSchi_3J); // lowest neutrino solver chi2 value
//            like_dir->second["3J_Event_Lowest_NSdisc"].fill(lowest_NSdisc_3J); // lowest neutrino solver disc value
//            like_dir->second["3J_Event_Lowest_Totaldisc"].fill(lowest_Totaldisc_3J); // lowest total (perm+NS) disc value
//
//
//            if( !best_perm.IsEmpty() ){
//                // values from best perm
//                like_dir->second["3J_Event_Best_Perm_Massdisc"].fill(best_perm.Merged3JDiscr()); // best perm perm disc value
//                like_dir->second["3J_Event_Best_Perm_NSchi"].fill(best_perm.NuChisq()); // best perm neutrino solver chi2 value
//                like_dir->second["3J_Event_Best_Perm_NSdisc"].fill(best_perm.NuDiscr()); // best perm neutrino solver disc value
//                like_dir->second["3J_Event_Best_Perm_Totaldisc"].fill(best_perm.Prob()); // best perm total (perm+NS) disc value
//            }
//
//            //        cout << "" << endl; 
            return best_perm;
        }



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


                // subdirectories
                auto gen_dir = histos_.find("Gen_Plots");

                //generator selection
                bool selection = genp_selector_.select(event);
                if( !selection ){
                    Logger::log().debug() << "event has no gen selection " << endl;
                    continue;
                }
                tracker_.track("gen selection");
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                //if( ttbar.M() > 700 ) continue;
                gen_dir->second["Mttbar"].fill(ttbar.M());

                njets_ = 0;
                if( object_selector_.select(event) ) njets_ = object_selector_.clean_jets().size();
                if( njets_ < 3 ) continue;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets_ << endl;
                tracker_.track("njet cuts");

                if( njets_  >= 4 ){
                    permutator_.reset();

                    // get matched perm from event
                    Permutation matched_perm_4PJ = matcher_.match(
                            genp_selector_.ttbar_final_system(),
                            object_selector_.clean_jets(),
                            object_selector_.veto_electrons(),
                            object_selector_.veto_muons());

                    // get rid of events that are merged or don't have b's and wjets
                    if( matched_perm_4PJ.Merged_Event() || !(matched_perm_4PJ.IsTHadComplete() && matched_perm_4PJ.BLep()) ) continue;

                    matched_4PJ_perms(event, matched_perm_4PJ);

//                    if( !matched_perm_4PJ.Merged_Event() ){ 
//                        cout << "Unmerged 4+ event" << endl;
//                    }

                    Permutation tp = process_4PJ_evt(event);
                    if( tp.IsEmpty() ) continue;
//                    cout << "bj1 pT: " << tp.BHad()->Pt() << endl;
//                    cout << "bj2 pT: " << tp.BLep()->Pt() << endl;
//                    cout << "wj1 pT: " << tp.WJa()->Pt() << endl;
//                    cout << "wj2 pT: " << tp.WJb()->Pt() << endl;
                }


                /// 3 jet events

                if( njets_ == 3 ){ // Likelihood value computation

                    cout << "njets = " << njets_ << endl;

                    //likelihood_3J_calc(event);

                    //if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue; // skip non-semilep events

                    //// get matched perm from event
                    //Permutation matched_perm = dr_matcher_.dr_match(
                    //        genp_selector_.ttbar_final_system(),
                    //        object_selector_.clean_jets(),
                    //        object_selector_.lepton(),
                    //        object_selector_.met(),
                    //        object_selector_.lepton_charge());

                    //if( matched_perm.IsEmpty() ){ // skip to next event if perm is empty
                    //    permutator_.reset_3J();
                    //    continue;
                    //}

                    //matched_perm_3J_plots(event, matched_perm);
                    //if( matched_perm.Merged_Event() ) merged_3J_evt(event, matched_perm);
                    //else unmerged_3J_evt(event, matched_perm);

                    //permutator_.reset_3J();
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
    URDriver<jet_perm_disc> test;
    int thing = test.run();
    return thing;
}
