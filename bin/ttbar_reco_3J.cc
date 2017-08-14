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

class ttbar_reco_3J : public AnalyzerBase
{
    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        //	    double njets_ = 0;

        double disc_cut_ = 2.;


        //histograms
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

            int nbins;
            double mass_max, mass_min;

            int mass_bins = 500;
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
            int massdisc_bins = 20;
            double massdisc_min = -10.;
            double massdisc_max = 10.;

            // ns disc
            int nsdisc_bins = 20;
            double nsdisc_min = 0.;
            double nsdisc_max = 10.;

            // total disc
            int combdisc_bins = 20;
            double combdisc_min = -5.;
            double combdisc_max = 15.;




            if( sample == "ttJetsM0" ){
                nbins = 22;
                //			m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min = 250.;
                mass_max = 2000.;
            }
            else if( sample == "ttJetsM700" ){
                nbins = 10;
                //			m_bins[nbins+1] = {700, 800, 900, 1000};
                mass_min = 700.;
                mass_max = 1000.;
            }
            else if( sample == "ttJetsM1000" ){
                nbins = 10;
                //			m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min = 1000.;
                mass_max = 2000.;
            }
            else{
                nbins = 10;
                mass_min = 200.;
                mass_max = 2000.;
            }

            string exp = "Expected_Plots";

            book<TH1F>(exp, "Expected_Event_Categories_3J", "", 5, 0.5, 5.5);
//            book<TH1F>(exp, "Expected_Event_Categories_4J", "", 5, 0.5, 5.5);
//            book<TH1F>(exp, "Expected_Event_Categories_5PJ", "", 5, 0.5, 5.5);

            //        const char *Objects[3] = {"thad", "tlep", "ttbar"};
            //        const char *kinvars[4] = {"Mass", "Pt", "Eta", "Costh"};
            const char *categories[4] = {"RIGHT", "MERGE_SWAP", "MERGE", "WRONG"};
            const char *cut_cats[2] = {"LCut", "GCut"};

            // Gen Plots
            string gen = "Gen_Plots";

            /// Mass
            book<TH1F>(gen, "thad_Mass", "", mass_bins, 125., 225.);
            book<TH1F>(gen, "tlep_Mass", "", mass_bins, 125., 225.);
            book<TH1F>(gen, "ttbar_Mass", "", nbins, mass_min, mass_max);

            ///Pt 
            book<TH1F>(gen, "thad_Pt", "", pt_bins, pt_min, pt_max);
            book<TH1F>(gen, "tlep_Pt", "", pt_bins, pt_min, pt_max);
            book<TH1F>(gen, "ttbar_Pt", "", pt_bins, pt_min, pt_max);

            ///Eta 
            book<TH1F>(gen, "thad_Eta", "", eta_bins, eta_min, eta_max);
            book<TH1F>(gen, "tlep_Eta", "", eta_bins, eta_min, eta_max);
            book<TH1F>(gen, "ttbar_Eta", "", eta_bins, eta_min, eta_max);

            ///Costh 
            book<TH1F>(gen, "thad_Costh", "", costh_bins, costh_min, costh_max);
            book<TH1F>(gen, "tlep_Costh", "", costh_bins, costh_min, costh_max);
            book<TH1F>(gen, "ttbar_Costh", "", costh_bins, costh_min, costh_max);



            // Reco jet plots
            string reco = "Reco_Plots";

            for( int i = 0; i < 2; ++i ){
                for( int j = 0; j < 4; ++j ){
                    /// Mass
                    book<TH1F>(reco, std::string("THad_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, 100., 250.);
                    book<TH1F>(reco, std::string("TLep_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, 100., 250.);
                    book<TH1F>(reco, std::string("TTbar_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, 200., 2000.);

                    ///Pt 
                    book<TH1F>(reco, std::string("THad_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, pt_min, pt_max);
                    book<TH1F>(reco, std::string("TLep_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, pt_min, pt_max);
                    book<TH1F>(reco, std::string("TTbar_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, pt_min, pt_max);

                    ///Eta 
                    book<TH1F>(reco, std::string("THad_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, eta_min, eta_max);
                    book<TH1F>(reco, std::string("TLep_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, eta_min, eta_max);
                    book<TH1F>(reco, std::string("TTbar_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, eta_min, eta_max);

                    ///Costh 
                    book<TH1F>(reco, std::string("THad_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, costh_min, costh_max);
                    book<TH1F>(reco, std::string("TLep_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, costh_min, costh_max);
                    book<TH1F>(reco, std::string("TTbar_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, costh_min, costh_max);
                }
            }


            // Jet Resolution plots (gen-reco)
            string resolution = "Resolution_Plots";

            for( int i = 0; i < 2; ++i ){
                for( int j = 0; j < 4; ++j ){
                    /// Mass
                    book<TH1F>(resolution, std::string("thad_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, -1000., 500.);
                    book<TH1F>(resolution, std::string("tlep_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, -500., 500.);
                    book<TH1F>(resolution, std::string("ttbar_Mass_")+cut_cats[i]+"_"+categories[j], "", mass_bins, -1000., mass_max);

                    ///Pt 
                    book<TH1F>(resolution, std::string("thad_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, -500., 500.);
                    book<TH1F>(resolution, std::string("tlep_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, -1000., 1000.);
                    book<TH1F>(resolution, std::string("ttbar_Pt_")+cut_cats[i]+"_"+categories[j], "", pt_bins, -1000, 500.);

                    ///Eta 
                    book<TH1F>(resolution, std::string("thad_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, -2*eta_max, 2*eta_max);
                    book<TH1F>(resolution, std::string("tlep_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, -2*eta_max, 2*eta_max);
                    book<TH1F>(resolution, std::string("ttbar_Eta_")+cut_cats[i]+"_"+categories[j], "", eta_bins, -2*eta_max, 2*eta_max);

                    ///Costh solution
                    book<TH1F>(resolution, std::string("thad_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, -2*costh_max, 2*costh_max);
                    book<TH1F>(resolution, std::string("tlep_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, -2*costh_max, 2*costh_max);
                    book<TH1F>(resolution, std::string("ttbar_Costh_")+cut_cats[i]+"_"+categories[j], "", costh_bins, -2*costh_max, 2*costh_max);
                }
            }

            // top tagger input vars

            // test perm jets
            string test_per = "Test_Perm_Plots";

            // Mass
            book<TH1F>(test_per, "3J_bj1wj1_Mass", "", mass_bins, 0., 400.);
            book<TH1F>(test_per, "3J_bj2wj1_Mass", "", mass_bins, 0., 400.);
            book<TH1F>(test_per, "3J_bj1_Mass", "", 50, 0., 200.);
            book<TH1F>(test_per, "3J_bj2_Mass", "", 50, 0., 200.);
            book<TH1F>(test_per, "3J_wj1_Mass", "", 50, 0., 200.);

            // CSV
            book<TH1F>(test_per, "3J_bj1_CSV", "", 50, 0., 1.);
            book<TH1F>(test_per, "3J_bj2_CSV", "", 50, 0., 1.);
            book<TH1F>(test_per, "3J_wj1_CSV", "", 50, 0., 1.);

            // CvsL
            book<TH1F>(test_per, "3J_bj1_CvsL", "", 20, -0.5, 1.);
            book<TH1F>(test_per, "3J_bj2_CvsL", "", 20, -0.5, 1.);
            book<TH1F>(test_per, "3J_wj1_CvsL", "", 20, -0.5, 1.);

            // Multiplicity
            book<TH1F>(test_per, "3J_bj1_Mult", "", 101, 0., 100.);
            book<TH1F>(test_per, "3J_bj2_Mult", "", 101, 0., 100.);
            book<TH1F>(test_per, "3J_wj1_Mult", "", 101, 0., 100.);

            // num btagged
            book<TH1F>(test_per, "Num_BTag", "", 4, -0.5, 3.5);

            // discriminant values
            book<TH1F>(test_per, "3J_tp_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(test_per, "3J_tp_NSdisc", "", nsdisc_bins, nsdisc_min, nsdisc_max);
            book<TH1F>(test_per, "3J_tp_Totaldisc", "", combdisc_bins, combdisc_min, combdisc_max);

            book<TH1F>(test_per, "Massdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(test_per, "NSdisc_cat", "", 3, 0.5, 3.5);
            book<TH1F>(test_per, "Totaldisc_cat", "", 3, 0.5, 3.5);

            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        // events with 3 jets
        Permutation process_3J_evt(URStreamer &event){
            auto test_dir = histos_.find("Test_Perm_Plots");
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

            test_dir->second["Num_BTag"].fill(num_btag);

            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });

            tracker_.track("All Num BTag");

            if( num_btag < 2 ) return empty_perm; // require at least 2 b-tagged jets

            tracker_.track("Num BTag >= 2");

            bj1 = jets_vector[0];
            bj2 = jets_vector[1];
            wj1 = jets_vector[2];


            // test perm plots
            // mass
            test_dir->second["3J_bj1wj1_Mass"].fill( (*bj1+*wj1).M() );
            test_dir->second["3J_bj2wj1_Mass"].fill( (*bj2+*wj1).M() );
            test_dir->second["3J_bj1_Mass"].fill( bj1->M() );
            test_dir->second["3J_bj2_Mass"].fill( bj2->M() );
            test_dir->second["3J_wj1_Mass"].fill( wj1->M() );

            // CSV
            test_dir->second["3J_bj1_CSV"].fill( bj1->csvIncl() );
            test_dir->second["3J_bj2_CSV"].fill( bj2->csvIncl() );
            test_dir->second["3J_wj1_CSV"].fill( wj1->csvIncl() );

            // CvsL
            test_dir->second["3J_bj1_CvsL"].fill( bj1->CvsLtag() );
            test_dir->second["3J_bj2_CvsL"].fill( bj2->CvsLtag() );
            test_dir->second["3J_wj1_CvsL"].fill( wj1->CvsLtag() );

            //Mult
            test_dir->second["3J_bj1_Mult"].fill( bj1->numberOfDaughters() );
            test_dir->second["3J_bj2_Mult"].fill( bj2->numberOfDaughters() );
            test_dir->second["3J_wj1_Mult"].fill( wj1->numberOfDaughters() );



            Permutation best_perm;
            //        double lowest_massdisc_3J = 1e10;
            //        double lowest_NSdisc_3J = 1e10;
            double lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
                solver_.Solve_3J(test_perm);

                test_dir->second["3J_tp_Massdisc"].fill(test_perm.PermDiscr());
                test_dir->second["3J_tp_NSdisc"].fill(test_perm.NuDiscr());
                test_dir->second["3J_tp_Totaldisc"].fill(test_perm.Prob());

                if( test_perm.PermDiscr() < 10. ) test_dir->second["Massdisc_cat"].fill(1); // 1 if mass disc is in range
                else test_dir->second["Massdisc_cat"].fill(3); // 3 if mass disc isn't in range
                if( test_perm.NuDiscr() < 10. && test_perm.NuDiscr() > 0. ) test_dir->second["NSdisc_cat"].fill(1); // 1 if NS disc is in range
                else test_dir->second["NSdisc_cat"].fill(3); // 3 if NS disc isn't in range
                if( test_perm.Prob() < 15. ) test_dir->second["Totaldisc_cat"].fill(1); // 1 if total disc is in range
                else test_dir->second["Totaldisc_cat"].fill(3); // 3 if total disc isn't in range

                if( test_perm.Prob() < lowest_Totaldisc_3J ){
                    lowest_Totaldisc_3J = test_perm.Prob();
                    best_perm = test_perm;
                }

            }
            //        cout << "" << endl; 
            return best_perm;
        } // end of process_3J_evt()


        void best_perm_cats( URStreamer &event ){
            permutator_.reset_3J();
            // subdirectories
            auto gen_dir = histos_.find("Gen_Plots");
            auto reco_dir = histos_.find("Reco_Plots");
            auto res_dir = histos_.find("Resolution_Plots");
            auto exp_dir = histos_.find("Expected_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            tracker_.track("Before best perm"); 

            // get best perm from event determined by likelihood disc
            Permutation best_perm = process_3J_evt(event);
            if( best_perm.IsEmpty() ){
                tracker_.track("best perm doesn't exist");
                return; // skip to next event if perm is empty
            }

            tracker_.track("best perm exists");

            hyp::TTbar gen_ttang(ttbar);

            auto gen_ttcm = gen_ttang.to_CM();

            double gen_ttbar_cth = gen_ttang.unit3D().Dot(gen_ttcm.unit3D());

            hyp::TTbar reco_ttang(best_perm); // doesn't work because not all wjets defined
            auto reco_ttcm = reco_ttang.to_CM();

            double reco_ttbar_cth = reco_ttang.unit3D().Dot(reco_ttcm.unit3D());
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            double reco_tlep_cth = reco_ttang.unit3D().Dot(reco_ttcm.tlep().unit3D());

            //Gen Plots 
            gen_dir->second["ttbar_Mass"].fill(ttbar.M());
            gen_dir->second["ttbar_Pt"].fill(ttbar.Pt());
            gen_dir->second["ttbar_Eta"].fill(ttbar.Eta());
            gen_dir->second["ttbar_Costh"].fill(gen_ttbar_cth);

            exp_dir->second["Expected_Event_Categories_3J"].fill(1);// tot expected events == 1

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty

                exp_dir->second["Expected_Event_Categories_3J"].fill(5);// expected other events == 5
                //Logger::log().debug() << "not semilep event " << evt_idx_ << endl;

                if( best_perm.Prob() < disc_cut_ ){
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_LCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_LCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_LCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_LCut_WRONG"].fill(reco_ttbar_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_LCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_LCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_LCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_LCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);
                }
                else{
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_GCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_GCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_GCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_GCut_WRONG"].fill(reco_ttbar_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_GCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_GCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_GCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_GCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);
                }

                tracker_.track("Not semilep events");

                return;
            }
            tracker_.track("semilep");

            if( ttbar.merged_bhadw_partons(0.4) || ttbar.merged_w_partons(0.4) ){ // gen partons merged
                exp_dir->second["Expected_Event_Categories_3J"].fill(3);// expected merged events == 3
            }
            else{ // gen partons not merged
                exp_dir->second["Expected_Event_Categories_3J"].fill(5);// expected other events == 5
            }

            double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());
            double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());

            //Gen Plots
            // Mass
            gen_dir->second["thad_Mass"].fill(ttbar.had_top()->M());
            gen_dir->second["tlep_Mass"].fill(ttbar.lep_top()->M());

            //Pt 
            gen_dir->second["thad_Pt"].fill(ttbar.had_top()->Pt());
            gen_dir->second["tlep_Pt"].fill(ttbar.lep_top()->Pt());

            //Eta 
            gen_dir->second["thad_Eta"].fill(ttbar.had_top()->Eta());
            gen_dir->second["tlep_Eta"].fill(ttbar.lep_top()->Eta());

            //Costh 
            gen_dir->second["thad_Costh"].fill(gen_thad_cth);
            gen_dir->second["tlep_Costh"].fill(gen_tlep_cth);



            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());

            if( mp.IsEmpty() ){ // skip to next event if perm is empty
                if( best_perm.Prob() < disc_cut_ ){
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_LCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_LCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_LCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_LCut_WRONG"].fill(reco_ttbar_cth);

                    reco_dir->second["THad_Mass_LCut_WRONG"].fill(best_perm.THad().M()); // Mass
                    reco_dir->second["TLep_Mass_LCut_WRONG"].fill(best_perm.TLep().M());
                    reco_dir->second["THad_Pt_LCut_WRONG"].fill(best_perm.THad().Pt()); // Pt
                    reco_dir->second["TLep_Pt_LCut_WRONG"].fill(best_perm.TLep().Pt());
                    reco_dir->second["THad_Eta_LCut_WRONG"].fill(best_perm.THad().Eta()); // Eta
                    reco_dir->second["TLep_Eta_LCut_WRONG"].fill(best_perm.TLep().Eta());
                    reco_dir->second["THad_Costh_LCut_WRONG"].fill(reco_thad_cth); // Costh
                    reco_dir->second["TLep_Costh_LCut_WRONG"].fill(reco_tlep_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_LCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_LCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_LCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_LCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);

                    res_dir->second["thad_Mass_LCut_WRONG"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                    res_dir->second["tlep_Mass_LCut_WRONG"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                    res_dir->second["thad_Pt_LCut_WRONG"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                    res_dir->second["tlep_Pt_LCut_WRONG"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                    res_dir->second["thad_Eta_LCut_WRONG"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                    res_dir->second["tlep_Eta_LCut_WRONG"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                    res_dir->second["thad_Costh_LCut_WRONG"].fill(gen_thad_cth - reco_thad_cth); // Costh
                    res_dir->second["tlep_Costh_LCut_WRONG"].fill(gen_tlep_cth - reco_tlep_cth);
                }
                else{
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_GCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_GCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_GCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_GCut_WRONG"].fill(reco_ttbar_cth);

                    reco_dir->second["THad_Mass_GCut_WRONG"].fill(best_perm.THad().M()); // Mass
                    reco_dir->second["TLep_Mass_GCut_WRONG"].fill(best_perm.TLep().M());
                    reco_dir->second["THad_Pt_GCut_WRONG"].fill(best_perm.THad().Pt()); // Pt
                    reco_dir->second["TLep_Pt_GCut_WRONG"].fill(best_perm.TLep().Pt());
                    reco_dir->second["THad_Eta_GCut_WRONG"].fill(best_perm.THad().Eta()); // Eta
                    reco_dir->second["TLep_Eta_GCut_WRONG"].fill(best_perm.TLep().Eta());
                    reco_dir->second["THad_Costh_GCut_WRONG"].fill(reco_thad_cth); // Costh
                    reco_dir->second["TLep_Costh_GCut_WRONG"].fill(reco_tlep_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_GCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_GCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_GCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_GCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);

                    res_dir->second["thad_Mass_GCut_WRONG"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                    res_dir->second["tlep_Mass_GCut_WRONG"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                    res_dir->second["thad_Pt_GCut_WRONG"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                    res_dir->second["tlep_Pt_GCut_WRONG"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                    res_dir->second["thad_Eta_GCut_WRONG"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                    res_dir->second["tlep_Eta_GCut_WRONG"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                    res_dir->second["thad_Costh_GCut_WRONG"].fill(gen_thad_cth - reco_thad_cth); // Costh
                    res_dir->second["tlep_Costh_GCut_WRONG"].fill(gen_tlep_cth - reco_tlep_cth);
                }

                tracker_.track("matched perm dne");

                return;
            }

            tracker_.track("matched perm exists");
            if( mp.Merged_Event() ){ // matched perm has merged objects

                // objects in best_perm and matched_perm are the same and in same position (RIGHT)
                if( mp.AreBsSame(best_perm) && ( best_perm.WJa() == mp.WJa() || best_perm.WJa() == mp.WJb() ) ){
                    if( best_perm.Prob() < disc_cut_ ){
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_LCut_RIGHT"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_LCut_RIGHT"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_LCut_RIGHT"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_LCut_RIGHT"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_LCut_RIGHT"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_LCut_RIGHT"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_LCut_RIGHT"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_LCut_RIGHT"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_LCut_RIGHT"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_LCut_RIGHT"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_LCut_RIGHT"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_LCut_RIGHT"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_LCut_RIGHT"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_LCut_RIGHT"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_LCut_RIGHT"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_LCut_RIGHT"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_LCut_RIGHT"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_LCut_RIGHT"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_LCut_RIGHT"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_LCut_RIGHT"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_LCut_RIGHT"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_LCut_RIGHT"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_LCut_RIGHT"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_LCut_RIGHT"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                    else{
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_GCut_RIGHT"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_GCut_RIGHT"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_GCut_RIGHT"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_GCut_RIGHT"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_GCut_RIGHT"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_GCut_RIGHT"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_GCut_RIGHT"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_GCut_RIGHT"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_GCut_RIGHT"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_GCut_RIGHT"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_GCut_RIGHT"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_GCut_RIGHT"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_GCut_RIGHT"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_GCut_RIGHT"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_GCut_RIGHT"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_GCut_RIGHT"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_GCut_RIGHT"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_GCut_RIGHT"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_GCut_RIGHT"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_GCut_RIGHT"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_GCut_RIGHT"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_GCut_RIGHT"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_GCut_RIGHT"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_GCut_RIGHT"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                }

                // objects in best_perm and matched_perm are the same but in differnt orders (MERGE_SWAP)
                else if( mp.AreBsFlipped(best_perm) && ( best_perm.WJa() == mp.WJa() || best_perm.WJa() == mp.WJb() ) ){
                    if( best_perm.Prob() < disc_cut_ ){
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_LCut_MERGE_SWAP"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_LCut_MERGE_SWAP"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_LCut_MERGE_SWAP"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_LCut_MERGE_SWAP"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_LCut_MERGE_SWAP"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_LCut_MERGE_SWAP"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_LCut_MERGE_SWAP"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_LCut_MERGE_SWAP"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_LCut_MERGE_SWAP"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_LCut_MERGE_SWAP"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_LCut_MERGE_SWAP"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_LCut_MERGE_SWAP"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_LCut_MERGE_SWAP"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_LCut_MERGE_SWAP"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_LCut_MERGE_SWAP"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_LCut_MERGE_SWAP"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_LCut_MERGE_SWAP"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_LCut_MERGE_SWAP"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_LCut_MERGE_SWAP"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_LCut_MERGE_SWAP"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_LCut_MERGE_SWAP"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_LCut_MERGE_SWAP"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_LCut_MERGE_SWAP"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_LCut_MERGE_SWAP"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                    else{
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_GCut_MERGE_SWAP"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_GCut_MERGE_SWAP"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_GCut_MERGE_SWAP"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_GCut_MERGE_SWAP"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_GCut_MERGE_SWAP"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_GCut_MERGE_SWAP"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_GCut_MERGE_SWAP"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_GCut_MERGE_SWAP"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_GCut_MERGE_SWAP"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_GCut_MERGE_SWAP"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_GCut_MERGE_SWAP"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_GCut_MERGE_SWAP"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_GCut_MERGE_SWAP"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_GCut_MERGE_SWAP"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_GCut_MERGE_SWAP"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_GCut_MERGE_SWAP"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_GCut_MERGE_SWAP"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_GCut_MERGE_SWAP"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_GCut_MERGE_SWAP"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_GCut_MERGE_SWAP"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_GCut_MERGE_SWAP"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_GCut_MERGE_SWAP"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_GCut_MERGE_SWAP"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_GCut_MERGE_SWAP"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                }

                // objects in best_perm and matched_perm are not all the same (MERGE)
                else{
                    if( best_perm.Prob() < disc_cut_ ){
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_LCut_MERGE"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_LCut_MERGE"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_LCut_MERGE"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_LCut_MERGE"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_LCut_MERGE"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_LCut_MERGE"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_LCut_MERGE"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_LCut_MERGE"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_LCut_MERGE"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_LCut_MERGE"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_LCut_MERGE"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_LCut_MERGE"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_LCut_MERGE"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_LCut_MERGE"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_LCut_MERGE"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_LCut_MERGE"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_LCut_MERGE"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_LCut_MERGE"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_LCut_MERGE"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_LCut_MERGE"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_LCut_MERGE"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_LCut_MERGE"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_LCut_MERGE"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_LCut_MERGE"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                    else{
                        //Reco Plots
                        reco_dir->second["TTbar_Mass_GCut_MERGE"].fill(best_perm.LVect().M());
                        reco_dir->second["TTbar_Pt_GCut_MERGE"].fill(best_perm.LVect().Pt());
                        reco_dir->second["TTbar_Eta_GCut_MERGE"].fill(best_perm.LVect().Eta());
                        reco_dir->second["TTbar_Costh_GCut_MERGE"].fill(reco_ttbar_cth);

                        reco_dir->second["THad_Mass_GCut_MERGE"].fill(best_perm.THad().M()); // Mass
                        reco_dir->second["TLep_Mass_GCut_MERGE"].fill(best_perm.TLep().M());
                        reco_dir->second["THad_Pt_GCut_MERGE"].fill(best_perm.THad().Pt()); // Pt
                        reco_dir->second["TLep_Pt_GCut_MERGE"].fill(best_perm.TLep().Pt());
                        reco_dir->second["THad_Eta_GCut_MERGE"].fill(best_perm.THad().Eta()); // Eta
                        reco_dir->second["TLep_Eta_GCut_MERGE"].fill(best_perm.TLep().Eta());
                        reco_dir->second["THad_Costh_GCut_MERGE"].fill(reco_thad_cth); // Costh
                        reco_dir->second["TLep_Costh_GCut_MERGE"].fill(reco_tlep_cth);

                        //Res Plots 
                        res_dir->second["ttbar_Mass_GCut_MERGE"].fill(ttbar.M() - best_perm.LVect().M());
                        res_dir->second["ttbar_Pt_GCut_MERGE"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                        res_dir->second["ttbar_Eta_GCut_MERGE"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                        res_dir->second["ttbar_Costh_GCut_MERGE"].fill(gen_ttbar_cth - reco_ttbar_cth);

                        res_dir->second["thad_Mass_GCut_MERGE"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                        res_dir->second["tlep_Mass_GCut_MERGE"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                        res_dir->second["thad_Pt_GCut_MERGE"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                        res_dir->second["tlep_Pt_GCut_MERGE"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                        res_dir->second["thad_Eta_GCut_MERGE"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                        res_dir->second["tlep_Eta_GCut_MERGE"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                        res_dir->second["thad_Costh_GCut_MERGE"].fill(gen_thad_cth - reco_thad_cth); // Costh
                        res_dir->second["tlep_Costh_GCut_MERGE"].fill(gen_tlep_cth - reco_tlep_cth);
                    }
                }

                tracker_.track("matched perm merged");

            }// end of merged events loop

            else{ // events not merged

                if( best_perm.Prob() < disc_cut_ ){
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_LCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_LCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_LCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_LCut_WRONG"].fill(reco_ttbar_cth);

                    reco_dir->second["THad_Mass_LCut_WRONG"].fill(best_perm.THad().M()); // Mass
                    reco_dir->second["TLep_Mass_LCut_WRONG"].fill(best_perm.TLep().M());
                    reco_dir->second["THad_Pt_LCut_WRONG"].fill(best_perm.THad().Pt()); // Pt
                    reco_dir->second["TLep_Pt_LCut_WRONG"].fill(best_perm.TLep().Pt());
                    reco_dir->second["THad_Eta_LCut_WRONG"].fill(best_perm.THad().Eta()); // Eta
                    reco_dir->second["TLep_Eta_LCut_WRONG"].fill(best_perm.TLep().Eta());
                    reco_dir->second["THad_Costh_LCut_WRONG"].fill(reco_thad_cth); // Costh
                    reco_dir->second["TLep_Costh_LCut_WRONG"].fill(reco_tlep_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_LCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_LCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_LCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_LCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);

                    res_dir->second["thad_Mass_LCut_WRONG"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                    res_dir->second["tlep_Mass_LCut_WRONG"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                    res_dir->second["thad_Pt_LCut_WRONG"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                    res_dir->second["tlep_Pt_LCut_WRONG"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                    res_dir->second["thad_Eta_LCut_WRONG"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                    res_dir->second["tlep_Eta_LCut_WRONG"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                    res_dir->second["thad_Costh_LCut_WRONG"].fill(gen_thad_cth - reco_thad_cth); // Costh
                    res_dir->second["tlep_Costh_LCut_WRONG"].fill(gen_tlep_cth - reco_tlep_cth);
                }
                else{
                    //Reco Plots
                    reco_dir->second["TTbar_Mass_GCut_WRONG"].fill(best_perm.LVect().M());
                    reco_dir->second["TTbar_Pt_GCut_WRONG"].fill(best_perm.LVect().Pt());
                    reco_dir->second["TTbar_Eta_GCut_WRONG"].fill(best_perm.LVect().Eta());
                    reco_dir->second["TTbar_Costh_GCut_WRONG"].fill(reco_ttbar_cth);

                    reco_dir->second["THad_Mass_GCut_WRONG"].fill(best_perm.THad().M()); // Mass
                    reco_dir->second["TLep_Mass_GCut_WRONG"].fill(best_perm.TLep().M());
                    reco_dir->second["THad_Pt_GCut_WRONG"].fill(best_perm.THad().Pt()); // Pt
                    reco_dir->second["TLep_Pt_GCut_WRONG"].fill(best_perm.TLep().Pt());
                    reco_dir->second["THad_Eta_GCut_WRONG"].fill(best_perm.THad().Eta()); // Eta
                    reco_dir->second["TLep_Eta_GCut_WRONG"].fill(best_perm.TLep().Eta());
                    reco_dir->second["THad_Costh_GCut_WRONG"].fill(reco_thad_cth); // Costh
                    reco_dir->second["TLep_Costh_GCut_WRONG"].fill(reco_tlep_cth);

                    //Res Plots 
                    res_dir->second["ttbar_Mass_GCut_WRONG"].fill(ttbar.M() - best_perm.LVect().M());
                    res_dir->second["ttbar_Pt_GCut_WRONG"].fill(ttbar.Pt() - best_perm.LVect().Pt());
                    res_dir->second["ttbar_Eta_GCut_WRONG"].fill(ttbar.Eta() - best_perm.LVect().Eta());
                    res_dir->second["ttbar_Costh_GCut_WRONG"].fill(gen_ttbar_cth - reco_ttbar_cth);

                    res_dir->second["thad_Mass_GCut_WRONG"].fill(ttbar.had_top()->M() - best_perm.THad().M()); // Mass
                    res_dir->second["tlep_Mass_GCut_WRONG"].fill(ttbar.lep_top()->M() - best_perm.TLep().M());
                    res_dir->second["thad_Pt_GCut_WRONG"].fill(ttbar.had_top()->Pt() - best_perm.THad().Pt()); // Pt
                    res_dir->second["tlep_Pt_GCut_WRONG"].fill(ttbar.lep_top()->Pt() - best_perm.TLep().Pt());
                    res_dir->second["thad_Eta_GCut_WRONG"].fill(ttbar.had_top()->Eta() - best_perm.THad().Eta()); // Eta
                    res_dir->second["tlep_Eta_GCut_WRONG"].fill(ttbar.lep_top()->Eta() - best_perm.TLep().Eta());
                    res_dir->second["thad_Costh_GCut_WRONG"].fill(gen_thad_cth - reco_thad_cth); // Costh
                    res_dir->second["tlep_Costh_GCut_WRONG"].fill(gen_tlep_cth - reco_tlep_cth);
                }

            } // end of unmerged events

            tracker_.track("matched perm unmerged");

        } // end of best_perm_cats


        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            while(event.next() /*&& evt_idx_ < 30000*/)
            {
//                auto exp_dir = histos_.find("Expected_Plots");

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

                //            if( ttbar.M() > 700 ) continue;

                int njets = 0;
                if( object_selector_.select(event) ) njets = object_selector_.clean_jets().size();
                if( njets < 3 ) continue;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets << endl;
                tracker_.track("njet cuts");

                /// 3 jet events
                if( njets == 3 ){ 
                    tracker_.track("njets = 3");
                    best_perm_cats(event);
                }

//                /// 4 jet events
//                else if( njets == 4 ){
//                    tracker_.track("njets = 4");
//                    exp_dir->second["Expected_Event_Categories_4J"].fill(1);// tot expeced events == 1
//                    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
//                        exp_dir->second["Expected_Event_Categories_4J"].fill(5);// expected other events == 5
//                        continue;
//                    }
//
//                    if( ttbar.merged_bhadw_partons(0.4) || ttbar.merged_w_partons(0.4) ){ // gen partons merged
//                        exp_dir->second["Expected_Event_Categories_4J"].fill(3);// expected merged events == 3
//                    }
//                    else{ // gen partons not merged
//                        exp_dir->second["Expected_Event_Categories_4J"].fill(5);// expected other events == 5
//                    }
//                }
//
//                /// 5+ jet events
//                else if( njets > 4 ){
//                    tracker_.track("njets = 5+");
//                    exp_dir->second["Expected_Event_Categories_5PJ"].fill(1);// tot expeced events == 1
//                    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
//                        exp_dir->second["Expected_Event_Categories_5PJ"].fill(5);// expected other events == 5
//                        continue;
//                    }
//
//                    if( ttbar.merged_bhadw_partons(0.4) || ttbar.merged_w_partons(0.4) ){ // gen partons merged
//                        exp_dir->second["Expected_Event_Categories_5PJ"].fill(3);// expected merged events == 3
//                    }
//                    else{ // gen partons not merged
//                        exp_dir->second["Expected_Event_Categories_5PJ"].fill(5);// expected other events == 5
//                    }
//                }

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
    URDriver<ttbar_reco_3J> test;
    int thing = test.run();
    return thing;
}
