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

        const char *DRnames_[2] = {"DRP4", "DRP8"};
        double DR_[2] = {0.4, 0.8};
        //            const char *DRnames_[4] = {"DRP4", "DRP5", "DRP6", "DRP8"};
        //	    double DR_[4] = {0.4, 0.5, 0.6, 0.8};
        //            const char *DRnames_[1] = {"DRP4"};
        //	    double DR_[1] = {0.4};

        double disc_cut_ = 2.;


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


        // Likelihood Discriminant
        double PermDiscr_3J = -1.;


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
        IDJet::BTag cut_tight_b_ = IDJet::BTag::MVATIGHT;
        IDJet::BTag cut_medium_b_ = IDJet::BTag::MVAMEDIUM;
        IDJet::BTag cut_loose_b_ = IDJet::BTag::MVALOOSE;


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

            // perm disc hists
            int massdisc_bins = 20;
            double massdisc_min = -10.;
            double massdisc_max = 10.;

            // ns disc hists
            int ns_disc_bins = 20;
            double ns_disc_min = 0.;
            double ns_disc_max = 10.;

            // ns chi hists
            int nschi_bins = 100;
            double nschi_min = 0.;
            double nschi_max = 1000.;

            // combined disc hists
            int comb_disc_bins = 20;
            double comb_disc_min = -5.;
            double comb_disc_max = 15.;

            // Gen Mttbar Plots
            string gen_plots = "Gen_Plots";
            book<TH1F>(gen_plots, "Mttbar", "", 18, 200., 2000.);


            // Likelihood Plots
            string likelihood = "Likelihood_Plots";

            /// Mass disc

            // 3 jets

            book<TH1F>(likelihood, "3J_Event_All_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(likelihood, "3J_Event_Lowest_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(likelihood, "3J_Event_Best_Perm_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);

            book<TH1F>(likelihood, "3J_Best_Perm_RIGHT_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_SWAP_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_WRONG_Massdisc", "", massdisc_bins, massdisc_min, massdisc_max);

            /// Neutrino solver 

            // 3 jets
            book<TH1F>(likelihood, "3J_Event_All_NSchi", "", nschi_bins, nschi_min, nschi_max);
            book<TH1F>(likelihood, "3J_Event_Lowest_NSchi", "", nschi_bins, nschi_min, nschi_max);
            book<TH1F>(likelihood, "3J_Event_Best_Perm_NSchi", "", nschi_bins, nschi_min, nschi_max);

            book<TH1F>(likelihood, "3J_Best_Perm_RIGHT_NSchi", "", nschi_bins, nschi_min, nschi_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_SWAP_NSchi", "", nschi_bins, nschi_min, nschi_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_NSchi", "", nschi_bins, nschi_min, nschi_max);
            book<TH1F>(likelihood, "3J_Best_Perm_WRONG_NSchi", "", nschi_bins, nschi_min, nschi_max);

            book<TH1F>(likelihood, "3J_Event_All_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);
            book<TH1F>(likelihood, "3J_Event_Lowest_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);
            book<TH1F>(likelihood, "3J_Event_Best_Perm_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);

            book<TH1F>(likelihood, "3J_Best_Perm_RIGHT_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_SWAP_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_WRONG_NSdisc", "", ns_disc_bins, ns_disc_min, ns_disc_max);

            /// Total disc 

            // 3 jets
            book<TH1F>(likelihood, "3J_Event_All_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);
            book<TH1F>(likelihood, "3J_Event_Lowest_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);
            book<TH1F>(likelihood, "3J_Event_Best_Perm_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);

            book<TH1F>(likelihood, "3J_Best_Perm_RIGHT_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_SWAP_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_MERGE_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);
            book<TH1F>(likelihood, "3J_Best_Perm_WRONG_Totaldisc", "", comb_disc_bins, comb_disc_min, comb_disc_max);


            // Best Perm
            string best_per = "Best_Perm_Plots";

            //perms w total disc values < 2
            // Mass
            book<TH1F>(best_per, "3J_BHadWJet_Mass_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHad_Mass_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_LCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHadWJet_Mass_RIGHT_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_MERGE_SWAP_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_MERGE_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_WRONG_LCut", "", mass_bins_, 0., 400.);

            book<TH1F>(best_per, "3J_BLepWJet_Mass_RIGHT_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_MERGE_SWAP_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_MERGE_LCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_WRONG_LCut", "", mass_bins_, 0., 400.);

            book<TH1F>(best_per, "3J_BHad_Mass_RIGHT_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_MERGE_SWAP_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_MERGE_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_WRONG_LCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BLep_Mass_RIGHT_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_MERGE_SWAP_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_MERGE_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_WRONG_LCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_WJet_Mass_RIGHT_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_MERGE_SWAP_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_MERGE_LCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_WRONG_LCut", "", 50, 0., 200.);

            // Pt 
            book<TH1F>(best_per, "3J_BHad_Pt_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_LCut", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(best_per, "3J_BHad_Pt_RIGHT_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_MERGE_SWAP_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_MERGE_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_WRONG_LCut", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(best_per, "3J_BLep_Pt_RIGHT_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_MERGE_SWAP_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_MERGE_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_WRONG_LCut", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(best_per, "3J_WJet_Pt_RIGHT_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_MERGE_SWAP_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_MERGE_LCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_WRONG_LCut", "", pt_bins_, pt_min_, pt_max_);

            // Eta 
            book<TH1F>(best_per, "3J_BHad_Eta_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_LCut", "", eta_bins_, eta_min_, eta_max_);

            book<TH1F>(best_per, "3J_BHad_Eta_RIGHT_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_MERGE_SWAP_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_MERGE_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_WRONG_LCut", "", eta_bins_, eta_min_, eta_max_);

            book<TH1F>(best_per, "3J_BLep_Eta_RIGHT_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_MERGE_SWAP_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_MERGE_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_WRONG_LCut", "", eta_bins_, eta_min_, eta_max_);

            book<TH1F>(best_per, "3J_WJet_Eta_RIGHT_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_MERGE_SWAP_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_MERGE_LCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_WRONG_LCut", "", eta_bins_, eta_min_, eta_max_);

            // Costh 
            book<TH1F>(best_per, "3J_BHad_Costh_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_LCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_BHad_Costh_RIGHT_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_MERGE_SWAP_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_MERGE_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_WRONG_LCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_BLep_Costh_RIGHT_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_MERGE_SWAP_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_MERGE_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_WRONG_LCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_WJet_Costh_RIGHT_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_MERGE_SWAP_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_MERGE_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_WRONG_LCut", "", pt_bins_, -1., 1.);

            // Mult
            book<TH1F>(best_per, "3J_BHad_Mult_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Mult_LCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Mult_LCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_BHad_Mult_RIGHT_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_MERGE_SWAP_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_MERGE_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_WRONG_LCut", "", 101, 0., 100.);

            book<TH1F>(best_per, "3J_BLep_Mult_RIGHT_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_MERGE_SWAP_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_MERGE_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_WRONG_LCut", "", 101, 0., 100.);

            book<TH1F>(best_per, "3J_WJet_Mult_RIGHT_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_MERGE_SWAP_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_MERGE_LCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_WRONG_LCut", "", 101, 0., 100.);


            //perms w total disc values > 2
            // Mass
            book<TH1F>(best_per, "3J_BHadWJet_Mass_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHad_Mass_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHadWJet_Mass_RIGHT_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_MERGE_SWAP_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_MERGE_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BHadWJet_Mass_WRONG_GCut", "", mass_bins_, 0., 400.);

            book<TH1F>(best_per, "3J_BLepWJet_Mass_RIGHT_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_MERGE_SWAP_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_MERGE_GCut", "", mass_bins_, 0., 400.);
            book<TH1F>(best_per, "3J_BLepWJet_Mass_WRONG_GCut", "", mass_bins_, 0., 400.);

            book<TH1F>(best_per, "3J_BHad_Mass_RIGHT_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_MERGE_SWAP_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_MERGE_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BHad_Mass_WRONG_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BLep_Mass_RIGHT_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_MERGE_SWAP_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_MERGE_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mass_WRONG_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_WJet_Mass_RIGHT_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_MERGE_SWAP_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_MERGE_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mass_WRONG_GCut", "", 50, 0., 200.);

            // Pt 
            book<TH1F>(best_per, "3J_BHad_Pt_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Pt_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Pt_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHad_Pt_RIGHT_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_MERGE_SWAP_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_MERGE_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BHad_Pt_WRONG_GCut", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(best_per, "3J_BLep_Pt_RIGHT_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_MERGE_SWAP_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_MERGE_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_BLep_Pt_WRONG_GCut", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(best_per, "3J_WJet_Pt_RIGHT_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_MERGE_SWAP_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_MERGE_GCut", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(best_per, "3J_WJet_Pt_WRONG_GCut", "", pt_bins_, pt_min_, pt_max_);

            // Eta 
            book<TH1F>(best_per, "3J_BHad_Eta_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Eta_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Eta_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHad_Eta_RIGHT_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_MERGE_SWAP_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_MERGE_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BHad_Eta_WRONG_GCut", "", eta_bins_, eta_min_, eta_max_);

            book<TH1F>(best_per, "3J_BLep_Eta_RIGHT_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_MERGE_SWAP_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_MERGE_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_BLep_Eta_WRONG_GCut", "", eta_bins_, eta_min_, eta_max_);

            book<TH1F>(best_per, "3J_WJet_Eta_RIGHT_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_MERGE_SWAP_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_MERGE_GCut", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(best_per, "3J_WJet_Eta_WRONG_GCut", "", eta_bins_, eta_min_, eta_max_);

            // Costh 
            book<TH1F>(best_per, "3J_BHad_Costh_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Costh_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Costh_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHad_Costh_RIGHT_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_MERGE_SWAP_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_MERGE_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BHad_Costh_WRONG_GCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_BLep_Costh_RIGHT_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_MERGE_SWAP_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_MERGE_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_BLep_Costh_WRONG_GCut", "", pt_bins_, -1., 1.);

            book<TH1F>(best_per, "3J_WJet_Costh_RIGHT_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_MERGE_SWAP_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_MERGE_GCut", "", pt_bins_, -1., 1.);
            book<TH1F>(best_per, "3J_WJet_Costh_WRONG_GCut", "", pt_bins_, -1., 1.);

            // Mult
            book<TH1F>(best_per, "3J_BHad_Mult_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_BLep_Mult_GCut", "", 50, 0., 200.);
            book<TH1F>(best_per, "3J_WJet_Mult_GCut", "", 50, 0., 200.);

            book<TH1F>(best_per, "3J_BHad_Mult_RIGHT_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_MERGE_SWAP_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_MERGE_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BHad_Mult_WRONG_GCut", "", 101, 0., 100.);

            book<TH1F>(best_per, "3J_BLep_Mult_RIGHT_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_MERGE_SWAP_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_MERGE_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_BLep_Mult_WRONG_GCut", "", 101, 0., 100.);

            book<TH1F>(best_per, "3J_WJet_Mult_RIGHT_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_MERGE_SWAP_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_MERGE_GCut", "", 101, 0., 100.);
            book<TH1F>(best_per, "3J_WJet_Mult_WRONG_GCut", "", 101, 0., 100.);




            // Merged jets discriminant plots
            string disc_right = "Correct_Disc_Plots";
            string disc_wrong = "Wrong_Disc_Plots";

            // 3 jets
            book<TH2D>(disc_right, "3J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
        book<TH2D>(disc_wrong, "3J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}


        // Neutrino Solver discriminant plots
        book<TH1F>(disc_right, "3J_nusolver_chi2_correct", "", 1000, 0., 1000.); // ns chi2 val hist for correct perms
        book<TH1F>(disc_wrong, "3J_nusolver_chi2_wrong", "", 1000, 0., 1000.); // ns chi2 val hist for wrong perms


        // Delta (gen-reco) Plots
        string delt = "Delta_Plots";


        // 3 jets
        // objects that pass preselection
        book<TH1F>(delt, "3J_Preselection_Correct_Perms_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
        book<TH1F>(delt, "3J_Preselection_Correct_Perms_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);
        book<TH1F>(delt, "3J_Preselection_Wrong_Perms_Delta_mass", "" , mass_bins_, -2*mass_max_, 2*mass_max_);
        book<TH1F>(delt, "3J_Preselection_Wrong_Perms_Delta_costh", "" , costh_bins_, -2*costh_max_, 2*costh_max_);

        // objects that pass the square cut
        book<TH1F>(delt, "3J_Square_Cut_Correct_Perms_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
        book<TH1F>(delt, "3J_Square_Cut_Correct_Perms_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);
        book<TH1F>(delt, "3J_Square_Cut_Wrong_Perms_Delta_mass", "" , mass_bins_, -2*mass_max_, 2*mass_max_);
        book<TH1F>(delt, "3J_Square_Cut_Wrong_Perms_Delta_costh", "" , costh_bins_, -2*costh_max_, 2*costh_max_);


        // top-tagger input vars

        // test perm jets
        string test_per = "Test_Perm_Plots";

        // Mass
        book<TH1F>(test_per, "3J_bj1wj1_Mass", "", mass_bins_, 0., 400.);
        book<TH1F>(test_per, "3J_bj2wj1_Mass", "", mass_bins_, 0., 400.);
        book<TH1F>(test_per, "3J_bj1_Mass", "", 50, 0., 200.);
        book<TH1F>(test_per, "3J_bj2_Mass", "", 50, 0., 200.);
        book<TH1F>(test_per, "3J_wj1_Mass", "", 50, 0., 200.);

        // MVA
        book<TH1F>(test_per, "3J_bj1_MVA", "", 50, -1., 1.);
        book<TH1F>(test_per, "3J_bj2_MVA", "", 50, -1., 1.);
        book<TH1F>(test_per, "3J_wj1_MVA", "", 50, -1., 1.);

        // CvsL
        book<TH1F>(test_per, "3J_bj1_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(test_per, "3J_bj2_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(test_per, "3J_wj1_CvsL", "", 20, -0.5, 1.);

        // Multiplicity
        book<TH1F>(test_per, "3J_bj1_Mult", "", 101, 0., 100.);
        book<TH1F>(test_per, "3J_bj2_Mult", "", 101, 0., 100.);
        book<TH1F>(test_per, "3J_wj1_Mult", "", 101, 0., 100.);


        // partially merged jets
        string partial_merge = "Partial_Merge_Matched_Perm_Plots";

        // Mass
        book<TH1F>(partial_merge, "3J_BHadWJet_Mass", "", mass_bins_, 0., 500.);
        book<TH1F>(partial_merge, "3J_BLepWJet_Mass", "", mass_bins_, 0., 500.);
        book<TH1F>(partial_merge, "3J_BHad_Mass", "", 50, 0., 200.);
        book<TH1F>(partial_merge, "3J_BLep_Mass", "", 50, 0., 200.);
        book<TH1F>(partial_merge, "3J_WJet_Mass", "", 50, 0., 200.);

        // MVA
        book<TH1F>(partial_merge, "3J_BHad_MVA", "", 50, -1., 1.);
        book<TH1F>(partial_merge, "3J_BLep_MVA", "", 50, -1., 1.);
        book<TH1F>(partial_merge, "3J_WJet_MVA", "", 50, -1., 1.);

        // CvsL
        book<TH1F>(partial_merge, "3J_BHad_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(partial_merge, "3J_BLep_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(partial_merge, "3J_WJet_CvsL", "", 20, -0.5, 1.);

        // Multiplicity
        book<TH1F>(partial_merge, "3J_BHad_Mult", "", 101, 0., 100.);
        book<TH1F>(partial_merge, "3J_BLep_Mult", "", 101, 0., 100.);
        book<TH1F>(partial_merge, "3J_WJet_Mult", "", 101, 0., 100.);


        // unmerged jets
        string unmerged = "Unmerged_Matched_Perm_Plots";

        // Mass
        book<TH1F>(unmerged, "3J_BHadWJet_Mass", "", mass_bins_, 0., 500.);
        book<TH1F>(unmerged, "3J_BLepWJet_Mass", "", mass_bins_, 0., 500.);
        book<TH1F>(unmerged, "3J_BHad_Mass", "", 50, 0., 200.);
        book<TH1F>(unmerged, "3J_BLep_Mass", "", 50, 0., 200.);
        book<TH1F>(unmerged, "3J_WJet_Mass", "", 50, 0., 200.);

        // MVA
        book<TH1F>(unmerged, "3J_BHad_MVA", "", 50, -1., 1.);
        book<TH1F>(unmerged, "3J_BLep_MVA", "", 50, -1., 1.);
        book<TH1F>(unmerged, "3J_WJet_MVA", "", 50, -1., 1.);

        // CvsL
        book<TH1F>(unmerged, "3J_BHad_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(unmerged, "3J_BLep_CvsL", "", 20, -0.5, 1.);
        book<TH1F>(unmerged, "3J_WJet_CvsL", "", 20, -0.5, 1.);

        // Multiplicity
        book<TH1F>(unmerged, "3J_BHad_Mult", "", 101, 0., 100.);
        book<TH1F>(unmerged, "3J_BLep_Mult", "", 101, 0., 100.);
        book<TH1F>(unmerged, "3J_WJet_Mult", "", 101, 0., 100.);

        Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        // events with 3 jets
        Permutation process_3J_evt(URStreamer &event){
            auto like_dir = histos_.find("Likelihood_Plots");
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

            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->CombinedMVA() > B->CombinedMVA() ); });

            if( num_btag < 2 ) return empty_perm; // require at least 2 b-tagged jets

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

            // MVA
            test_dir->second["3J_bj1_MVA"].fill( bj1->CombinedMVA() );
            test_dir->second["3J_bj2_MVA"].fill( bj2->CombinedMVA() );
            test_dir->second["3J_wj1_MVA"].fill( wj1->CombinedMVA() );

            // CvsL
            test_dir->second["3J_bj1_CvsL"].fill( bj1->CvsLtag() );
            test_dir->second["3J_bj2_CvsL"].fill( bj2->CvsLtag() );
            test_dir->second["3J_wj1_CvsL"].fill( wj1->CvsLtag() );

            //Mult 
            test_dir->second["3J_bj1_Mult"].fill( bj1->numberOfDaughters() );
            test_dir->second["3J_bj2_Mult"].fill( bj2->numberOfDaughters() );
            test_dir->second["3J_wj1_Mult"].fill( wj1->numberOfDaughters() );




            Permutation best_perm;
            double lowest_massdisc_3J = 1e10;
            double lowest_NSchi_3J = 1e10;
            double lowest_NSdisc_3J = 1e10;
            double lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
                solver_.Solve_3J(test_perm);

                if( test_perm.PermDiscr() < lowest_massdisc_3J ){
                    lowest_massdisc_3J = test_perm.PermDiscr();
                    //                best_perm = test_perm;
                }
                like_dir->second["3J_Event_All_Massdisc"].fill(test_perm.PermDiscr()); // all perm disc values

                if( test_perm.NuChisq() < lowest_NSchi_3J && test_perm.NuChisq() > 0. ){
                    lowest_NSchi_3J = test_perm.NuChisq();
                    //                cout << "NS chi: " << lowest_NSchi_3J << endl;
                }
                like_dir->second["3J_Event_All_NSchi"].fill(test_perm.NuChisq()); // all neutrino solver chi2 values

                if( test_perm.NuDiscr() < lowest_NSdisc_3J ){
                    lowest_NSdisc_3J = test_perm.NuDiscr();
                    //                cout << "NS disc: " << lowest_NSdisc_3J << endl;
                }
                like_dir->second["3J_Event_All_NSdisc"].fill(test_perm.NuDiscr()); // all NS disc values

                if( test_perm.Prob() < lowest_Totaldisc_3J ){
                    lowest_Totaldisc_3J = test_perm.Prob();
                    best_perm = test_perm;
                    //                cout << "NS disc: " << lowest_NSdisc_3J << endl;
                }
                like_dir->second["3J_Event_All_Totaldisc"].fill(test_perm.Prob()); // all total (perm + NS) disc values

            }

            // lowest values from event
            like_dir->second["3J_Event_Lowest_Massdisc"].fill(lowest_massdisc_3J); // lowest perm disc value
            like_dir->second["3J_Event_Lowest_NSchi"].fill(lowest_NSchi_3J); // lowest neutrino solver chi2 value
            like_dir->second["3J_Event_Lowest_NSdisc"].fill(lowest_NSdisc_3J); // lowest neutrino solver disc value
            like_dir->second["3J_Event_Lowest_Totaldisc"].fill(lowest_Totaldisc_3J); // lowest total (perm+NS) disc value


            if( !best_perm.IsEmpty() ){
                // values from best perm
                like_dir->second["3J_Event_Best_Perm_Massdisc"].fill(best_perm.PermDiscr()); // best perm perm disc value
                like_dir->second["3J_Event_Best_Perm_NSchi"].fill(best_perm.NuChisq()); // best perm neutrino solver chi2 value
                like_dir->second["3J_Event_Best_Perm_NSdisc"].fill(best_perm.NuDiscr()); // best perm neutrino solver disc value
                like_dir->second["3J_Event_Best_Perm_Totaldisc"].fill(best_perm.Prob()); // best perm total (perm+NS) disc value
            }

            //        cout << "" << endl; 
            return best_perm;
        }


        // likelihood value calculation/categorization
        void likelihood_3J_calc(URStreamer &event){
            auto like_dir = histos_.find("Likelihood_Plots");
            auto best_dir = histos_.find("Best_Perm_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            // get best perm from event determined by likelihood disc
            Permutation best_perm = process_3J_evt(event);
            if( best_perm.IsEmpty() ){
                permutator_.reset_3J();
                return; // skip to next event if perm is empty
            }

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                like_dir->second["3J_Best_Perm_WRONG_Massdisc"].fill(best_perm.PermDiscr());
                like_dir->second["3J_Best_Perm_WRONG_NSchi"].fill(best_perm.NuChisq());
                like_dir->second["3J_Best_Perm_WRONG_NSdisc"].fill(best_perm.NuDiscr());
                like_dir->second["3J_Best_Perm_WRONG_Totaldisc"].fill(best_perm.Prob());

                if( best_perm.Prob() < disc_cut_ ){
                    best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                    best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                    best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                    best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                    best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                    best_dir->second["3J_BHadWJet_Mass_WRONG_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

                else{
                    best_dir->second["3J_BHadWJet_Mass_WRONG_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_GCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_GCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_GCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_GCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_GCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_GCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_GCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_GCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_GCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

                permutator_.reset_3J();
                return;
            }
            tracker_.track("3J semilep events");

            // get matched perm from event
            Permutation matched_perm = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());


            if( matched_perm.IsEmpty() ){ // skip to next event if perm is empty
                like_dir->second["3J_Best_Perm_WRONG_Massdisc"].fill(best_perm.PermDiscr());
                like_dir->second["3J_Best_Perm_WRONG_NSchi"].fill(best_perm.NuChisq());
                like_dir->second["3J_Best_Perm_WRONG_NSdisc"].fill(best_perm.NuDiscr());
                like_dir->second["3J_Best_Perm_WRONG_Totaldisc"].fill(best_perm.Prob());

                if( best_perm.Prob() < disc_cut_ ){
                    best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                    best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                    best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                    best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                    best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                    best_dir->second["3J_BHadWJet_Mass_WRONG_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

                else{
                    best_dir->second["3J_BHadWJet_Mass_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_GCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_GCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_GCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_GCut"].fill( best_perm.BHad()->Pt() ); // pt
                    best_dir->second["3J_BLep_Pt_GCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_GCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_GCut"].fill( best_perm.BHad()->Eta() ); // eta
                    best_dir->second["3J_BLep_Eta_GCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_GCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_GCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                    best_dir->second["3J_BLep_Costh_GCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_GCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                    best_dir->second["3J_BLep_Mult_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_GCut"].fill( best_perm.WJa()->numberOfDaughters() );

                    best_dir->second["3J_BHadWJet_Mass_WRONG_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_GCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_GCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_GCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_GCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_GCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_GCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_GCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_GCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_GCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

                permutator_.reset_3J();
                return;
            }
            tracker_.track("matched events");


            if( matched_perm.Merged_Event() ){ // matched perm has merged objets

                // objects in best_perm and matched_perm are the same and in same position
                if( matched_perm.AreBsSame(best_perm) && ( best_perm.WJa() == matched_perm.WJa() || best_perm.WJa() == matched_perm.WJb() ) ){ 
                    like_dir->second["3J_Best_Perm_RIGHT_Massdisc"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_Best_Perm_RIGHT_NSchi"].fill(best_perm.NuChisq());
                    like_dir->second["3J_Best_Perm_RIGHT_NSdisc"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Best_Perm_RIGHT_Totaldisc"].fill(best_perm.Prob());

                    if( best_perm.Prob() < disc_cut_ ){
                        best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_RIGHT_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_RIGHT_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_RIGHT_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_RIGHT_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_RIGHT_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_RIGHT_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_RIGHT_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_RIGHT_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_RIGHT_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_RIGHT_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_RIGHT_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_RIGHT_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_RIGHT_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_RIGHT_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_RIGHT_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_RIGHT_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_RIGHT_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                    else{
                        best_dir->second["3J_BHadWJet_Mass_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_GCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_GCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_GCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_GCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_RIGHT_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_RIGHT_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_RIGHT_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_RIGHT_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_RIGHT_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_RIGHT_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_RIGHT_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_RIGHT_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_RIGHT_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_RIGHT_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_RIGHT_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_RIGHT_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_RIGHT_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_RIGHT_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_RIGHT_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_RIGHT_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_RIGHT_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                }

                // objects in best_perm and matched_perm are the same but in differnt orders
                else if( matched_perm.AreBsFlipped(best_perm) && ( best_perm.WJa() == matched_perm.WJa() || best_perm.WJa() == matched_perm.WJb() ) ){
                    like_dir->second["3J_Best_Perm_MERGE_SWAP_Massdisc"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_Best_Perm_MERGE_SWAP_NSchi"].fill(best_perm.NuChisq());
                    like_dir->second["3J_Best_Perm_MERGE_SWAP_NSdisc"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Best_Perm_MERGE_SWAP_Totaldisc"].fill(best_perm.Prob());

                    if( best_perm.Prob() < disc_cut_ ){
                        best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_MERGE_SWAP_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_MERGE_SWAP_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_MERGE_SWAP_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_MERGE_SWAP_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_MERGE_SWAP_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_MERGE_SWAP_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_MERGE_SWAP_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_MERGE_SWAP_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_MERGE_SWAP_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_MERGE_SWAP_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_MERGE_SWAP_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_MERGE_SWAP_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_MERGE_SWAP_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_MERGE_SWAP_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_MERGE_SWAP_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_MERGE_SWAP_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_MERGE_SWAP_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                    else{
                        best_dir->second["3J_BHadWJet_Mass_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_GCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_GCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_GCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_GCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_MERGE_SWAP_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_MERGE_SWAP_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_MERGE_SWAP_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_MERGE_SWAP_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_MERGE_SWAP_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_MERGE_SWAP_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_MERGE_SWAP_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_MERGE_SWAP_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_MERGE_SWAP_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_MERGE_SWAP_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_MERGE_SWAP_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_MERGE_SWAP_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_MERGE_SWAP_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_MERGE_SWAP_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_MERGE_SWAP_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_MERGE_SWAP_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_MERGE_SWAP_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                }

                // objects in best_perm and matched_perm are not all the same
                else{
                    like_dir->second["3J_Best_Perm_MERGE_Massdisc"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_Best_Perm_MERGE_NSchi"].fill(best_perm.NuChisq());
                    like_dir->second["3J_Best_Perm_MERGE_NSdisc"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Best_Perm_MERGE_Totaldisc"].fill(best_perm.Prob());

                    if( best_perm.Prob() < disc_cut_ ){
                        best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_MERGE_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_MERGE_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_MERGE_LCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_MERGE_LCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_MERGE_LCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_MERGE_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_MERGE_LCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_MERGE_LCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_MERGE_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_MERGE_LCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_MERGE_LCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_MERGE_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_MERGE_LCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_MERGE_LCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_MERGE_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_MERGE_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_MERGE_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                    else{
                        best_dir->second["3J_BHadWJet_Mass_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_GCut"].fill( best_perm.BHad()->Pt() ); // pt
                        best_dir->second["3J_BLep_Pt_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_GCut"].fill( best_perm.BHad()->Eta() ); // eta
                        best_dir->second["3J_BLep_Eta_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_GCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                        best_dir->second["3J_BLep_Costh_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                        best_dir->second["3J_BLep_Mult_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_GCut"].fill( best_perm.WJa()->numberOfDaughters() );

                        best_dir->second["3J_BHadWJet_Mass_MERGE_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                        best_dir->second["3J_BLepWJet_Mass_MERGE_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                        best_dir->second["3J_BHad_Mass_MERGE_GCut"].fill( best_perm.BHad()->M() );
                        best_dir->second["3J_BLep_Mass_MERGE_GCut"].fill( best_perm.BLep()->M() );
                        best_dir->second["3J_WJet_Mass_MERGE_GCut"].fill( best_perm.WJa()->M() );

                        best_dir->second["3J_BHad_Pt_MERGE_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                        best_dir->second["3J_BLep_Pt_MERGE_GCut"].fill( best_perm.BLep()->Pt() );
                        best_dir->second["3J_WJet_Pt_MERGE_GCut"].fill( best_perm.WJa()->Pt() );

                        best_dir->second["3J_BHad_Eta_MERGE_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                        best_dir->second["3J_BLep_Eta_MERGE_GCut"].fill( best_perm.BLep()->Eta() );
                        best_dir->second["3J_WJet_Eta_MERGE_GCut"].fill( best_perm.WJa()->Eta() );

                        best_dir->second["3J_BHad_Costh_MERGE_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                        best_dir->second["3J_BLep_Costh_MERGE_GCut"].fill( best_perm.BLep()->CosTheta() );
                        best_dir->second["3J_WJet_Costh_MERGE_GCut"].fill( best_perm.WJa()->CosTheta() );

                        best_dir->second["3J_BHad_Mult_MERGE_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                        best_dir->second["3J_BLep_Mult_MERGE_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                        best_dir->second["3J_WJet_Mult_MERGE_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                    }

                }
            }

            else{ // events not merged
                like_dir->second["3J_Best_Perm_WRONG_Massdisc"].fill(best_perm.PermDiscr());
                like_dir->second["3J_Best_Perm_WRONG_NSchi"].fill(best_perm.NuChisq());
                like_dir->second["3J_Best_Perm_WRONG_NSdisc"].fill(best_perm.NuDiscr());
                like_dir->second["3J_Best_Perm_WRONG_Totaldisc"].fill(best_perm.Prob());

                if( best_perm.Prob() < disc_cut_ ){
                    best_dir->second["3J_BHadWJet_Mass_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_LCut"].fill( best_perm.BHad()->Pt() ); // pt
                    best_dir->second["3J_BLep_Pt_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_LCut"].fill( best_perm.BHad()->Eta() ); // eta
                    best_dir->second["3J_BLep_Eta_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_LCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                    best_dir->second["3J_BLep_Costh_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                    best_dir->second["3J_BLep_Mult_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_LCut"].fill( best_perm.WJa()->numberOfDaughters() );

                    best_dir->second["3J_BHadWJet_Mass_WRONG_LCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_LCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_LCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_LCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_LCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_LCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_LCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_LCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_LCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_LCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_LCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_LCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_LCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_LCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_LCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_LCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_LCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

                else{
                    best_dir->second["3J_BHadWJet_Mass_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_GCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_GCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_GCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_GCut"].fill( best_perm.BHad()->Pt() ); // pt
                    best_dir->second["3J_BLep_Pt_GCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_GCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_GCut"].fill( best_perm.BHad()->Eta() ); // eta
                    best_dir->second["3J_BLep_Eta_GCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_GCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_GCut"].fill( best_perm.BHad()->CosTheta() ); // costh
                    best_dir->second["3J_BLep_Costh_GCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_GCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); // Mult
                    best_dir->second["3J_BLep_Mult_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_GCut"].fill( best_perm.WJa()->numberOfDaughters() );

                    best_dir->second["3J_BHadWJet_Mass_WRONG_GCut"].fill( (*best_perm.BHad()+*best_perm.WJa()).M() ); //mass
                    best_dir->second["3J_BLepWJet_Mass_WRONG_GCut"].fill( (*best_perm.BLep()+*best_perm.WJa()).M() );
                    best_dir->second["3J_BHad_Mass_WRONG_GCut"].fill( best_perm.BHad()->M() );
                    best_dir->second["3J_BLep_Mass_WRONG_GCut"].fill( best_perm.BLep()->M() );
                    best_dir->second["3J_WJet_Mass_WRONG_GCut"].fill( best_perm.WJa()->M() );

                    best_dir->second["3J_BHad_Pt_WRONG_GCut"].fill( best_perm.BHad()->Pt() ); //pt
                    best_dir->second["3J_BLep_Pt_WRONG_GCut"].fill( best_perm.BLep()->Pt() );
                    best_dir->second["3J_WJet_Pt_WRONG_GCut"].fill( best_perm.WJa()->Pt() );

                    best_dir->second["3J_BHad_Eta_WRONG_GCut"].fill( best_perm.BHad()->Eta() ); //eta
                    best_dir->second["3J_BLep_Eta_WRONG_GCut"].fill( best_perm.BLep()->Eta() );
                    best_dir->second["3J_WJet_Eta_WRONG_GCut"].fill( best_perm.WJa()->Eta() );

                    best_dir->second["3J_BHad_Costh_WRONG_GCut"].fill( best_perm.BHad()->CosTheta() ); //costh
                    best_dir->second["3J_BLep_Costh_WRONG_GCut"].fill( best_perm.BLep()->CosTheta() );
                    best_dir->second["3J_WJet_Costh_WRONG_GCut"].fill( best_perm.WJa()->CosTheta() );

                    best_dir->second["3J_BHad_Mult_WRONG_GCut"].fill( best_perm.BHad()->numberOfDaughters() ); //mult
                    best_dir->second["3J_BLep_Mult_WRONG_GCut"].fill( best_perm.BLep()->numberOfDaughters() );
                    best_dir->second["3J_WJet_Mult_WRONG_GCut"].fill( best_perm.WJa()->numberOfDaughters() );
                }

            } 

            permutator_.reset_3J();

        } // end of likelihood_3J_calc()


        // right/wrong matched perm plots
        void matched_perm_3J_plots(URStreamer &event, Permutation &matched_perm){

            auto disc_wrong_dir = histos_.find("Wrong_Disc_Plots");
            auto disc_correct_dir = histos_.find("Correct_Disc_Plots");
            auto delt_dir = histos_.find("Delta_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            GenTop* thad = ttbar.had_top();


            solver_.Solve_3J(matched_perm);


            if( !matched_perm.Merged_Event() ){ // no jets merged
                if( matched_perm.WJa() ){ // if wja exists
                    // blep paired with wja

                    Unmerged_BLep_mass_3J = matched_perm.BLep()->M();
                    Unmerged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                    Max_Unmerged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                    //                      t1->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());

                    // bhad paired with wja

                    Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                    Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                    Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
                    //                t3->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());

                }

                if( matched_perm.WJb() ){ // if wjb exists
                    // blep paired with wjb

                    Unmerged_BLep_mass_3J = matched_perm.BLep()->M();
                    Unmerged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJb()).M();
                    Max_Unmerged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJb()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJb()->M();
                    //              t1->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());

                    // bhad paired with wjb

                    Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                    Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                    Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
                    //                t3->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());

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

                // preselection
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());

            }

            //wrong perm
            if( matched_perm.Merged_BLepWJb() && matched_perm.WJa() ){ // only blep and wjb merged and wja exists
                Merged_BLep_mass_3J = matched_perm.BLep()->M(); 
                Merged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                Max_Merged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                //            t2->Fill();
                disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);
                disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                // preselection
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());

            }

            // correct perm
            if( matched_perm.Merged_BHadWJa() && matched_perm.WJb() ){ // only bhad and wja merged and wjb exists
                Merged_BHad_mass_3J = matched_perm.BHad()->M();
                Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
                //            t4->Fill();
                disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

                // preselection
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());

            }

            // correct perm
            if( matched_perm.Merged_BHadWJb() && matched_perm.WJa() ){ // only bhad and wjb merged and wja exists
                Merged_BHad_mass_3J =  matched_perm.BHad()->M();
                Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
                //            t4->Fill();
                disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

                // preselection
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());

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

                // preselection
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());


                //wrong perm
                Wrong_merged_wjet_mass_3J = matched_perm.BLep()->M();
                Wrong_merged_wjet_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                Max_Wrong_merged_wjet_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
                //			t6->Fill();
                disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Wrong_merged_wjet_mass_3J, Wrong_merged_wjet_comb_mass_3J);
                //            disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                // preselection
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());


            } // end of merged W jets loop

            permutator_.reset_3J();

        } // end of matched_perm_3J_plots()



        void merged_3J_evt(URStreamer &event, Permutation &hyp){

            auto merged_dir = histos_.find("Partial_Merge_Matched_Perm_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            if( !(hyp.WJa() && hyp.WJb()) ) return; // require all perm objects to exist (2 of them the same)
            tracker_.track("Merged 3J all events");

            if( hyp.Merged_BHadWJb() && hyp.WJa() ){ // bhad and wjb merged
                // invariant mass
                merged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJa()).M() );
                merged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJa()).M() );

                merged_dir->second["3J_WJet_Mass"].fill( hyp.WJa()->M() ); // jet mass
                merged_dir->second["3J_WJet_MVA"].fill( hyp.WJa()->CombinedMVA() ); // jet csv
                merged_dir->second["3J_WJet_CvsL"].fill( hyp.WJa()->CvsLtag() ); // jet CvsL
                merged_dir->second["3J_WJet_Mult"].fill( hyp.WJa()->numberOfDaughters() ); // jet mult

                tracker_.track("Merged 3J category");
            }

            if( hyp.Merged_BHadWJa() && hyp.WJb() ){ // bhad and wja merged
                // invariant mass
                merged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJb()).M() );
                merged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJb()).M() );

                merged_dir->second["3J_WJet_Mass"].fill( hyp.WJb()->M() ); // jet mass
                merged_dir->second["3J_WJet_MVA"].fill( hyp.WJb()->CombinedMVA() ); // jet csv
                merged_dir->second["3J_WJet_CvsL"].fill( hyp.WJb()->CvsLtag() ); // jet CvsL
                merged_dir->second["3J_WJet_Mult"].fill( hyp.WJb()->numberOfDaughters() ); // jet mult

                tracker_.track("Merged 3J category");
            }

            if( hyp.Merged_WJets() ){ // wjets merged
                // invariant mass
                merged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJa()).M() );
                merged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJa()).M() );

                merged_dir->second["3J_WJet_Mass"].fill( hyp.WJa()->M() ); // jet mass
                merged_dir->second["3J_WJet_MVA"].fill( hyp.WJa()->CombinedMVA() ); // jet csv
                merged_dir->second["3J_WJet_CvsL"].fill( hyp.WJa()->CvsLtag() ); // jet CvsL
                merged_dir->second["3J_WJet_Mult"].fill( hyp.WJa()->numberOfDaughters() ); // jet mult

                tracker_.track("Merged 3J category");
            }

            if( hyp.Merged_BLepWJb() && hyp.WJa() ){ // blep and wjb merged
                // invariant mass
                merged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJa()).M() );
                merged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJa()).M() );

                merged_dir->second["3J_WJet_Mass"].fill( hyp.WJa()->M() ); // jet mass
                merged_dir->second["3J_WJet_MVA"].fill( hyp.WJa()->CombinedMVA() ); // jet csv
                merged_dir->second["3J_WJet_CvsL"].fill( hyp.WJa()->CvsLtag() ); // jet CvsL
                merged_dir->second["3J_WJet_Mult"].fill( hyp.WJa()->numberOfDaughters() ); // jet mult

                tracker_.track("Merged 3J category");
            }

            if( hyp.Merged_BLepWJa() && hyp.WJb() ){ // blep and wja merged
                // invariant mass
                merged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJb()).M() );
                merged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJb()).M() );

                merged_dir->second["3J_WJet_Mass"].fill( hyp.WJb()->M() ); // jet mass
                merged_dir->second["3J_WJet_MVA"].fill( hyp.WJb()->CombinedMVA() ); // jet csv
                merged_dir->second["3J_WJet_CvsL"].fill( hyp.WJb()->CvsLtag() ); // jet CvsL
                merged_dir->second["3J_WJet_Mult"].fill( hyp.WJb()->numberOfDaughters() ); // jet mult

                tracker_.track("Merged 3J category");
            }


            // BHad plots        
            merged_dir->second["3J_BHad_Mass"].fill( hyp.BHad()->M() );
            merged_dir->second["3J_BHad_MVA"].fill( hyp.BHad()->CombinedMVA() );
            merged_dir->second["3J_BHad_CvsL"].fill( hyp.BHad()->CvsLtag() );
            merged_dir->second["3J_BHad_Mult"].fill( hyp.BHad()->numberOfDaughters() );

            // BLep plots        
            merged_dir->second["3J_BLep_Mass"].fill( hyp.BLep()->M() );
            merged_dir->second["3J_BLep_MVA"].fill( hyp.BLep()->CombinedMVA() );
            merged_dir->second["3J_BLep_CvsL"].fill( hyp.BLep()->CvsLtag() );
            merged_dir->second["3J_BLep_Mult"].fill( hyp.BLep()->numberOfDaughters() );


        } // end of merged_3J_evt


        void unmerged_3J_evt(URStreamer &event, Permutation &hyp){

            auto unmerged_dir = histos_.find("Unmerged_Matched_Perm_Plots");

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            tracker_.track("gen selection");
            GenTTBar &ttbar = genp_selector_.ttbar_system();


            if( hyp.WJa() ){
                // invariant mass
                unmerged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJa()).M() );
                unmerged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJa()).M() );

                unmerged_dir->second["3J_WJet_Mass"].fill( hyp.WJa()->M() ); // jet mass
                unmerged_dir->second["3J_WJet_MVA"].fill( hyp.WJa()->CombinedMVA() ); // jet csv
                unmerged_dir->second["3J_WJet_CvsL"].fill( hyp.WJa()->CvsLtag() ); // jet CvsL
                unmerged_dir->second["3J_WJet_Mult"].fill( hyp.WJa()->numberOfDaughters() ); // jet mult
            }

            if( hyp.WJb() ){
                // invariant mass
                unmerged_dir->second["3J_BHadWJet_Mass"].fill( (*hyp.BHad()+*hyp.WJb()).M() );
                unmerged_dir->second["3J_BLepWJet_Mass"].fill( (*hyp.BLep()+*hyp.WJb()).M() );

                unmerged_dir->second["3J_WJet_Mass"].fill( hyp.WJb()->M() ); // jet mass
                unmerged_dir->second["3J_WJet_MVA"].fill( hyp.WJb()->CombinedMVA() ); // jet csv
                unmerged_dir->second["3J_WJet_CvsL"].fill( hyp.WJb()->CvsLtag() ); // jet CvsL
                unmerged_dir->second["3J_WJet_Mult"].fill( hyp.WJb()->numberOfDaughters() ); // jet mult
            }

            // BHad plots        
            unmerged_dir->second["3J_BHad_Mass"].fill( hyp.BHad()->M() );
            unmerged_dir->second["3J_BHad_MVA"].fill( hyp.BHad()->CombinedMVA() );
            unmerged_dir->second["3J_BHad_CvsL"].fill( hyp.BHad()->CvsLtag() );
            unmerged_dir->second["3J_BHad_Mult"].fill( hyp.BHad()->numberOfDaughters() );

            // BLep plots        
            unmerged_dir->second["3J_BLep_Mass"].fill( hyp.BLep()->M() );
            unmerged_dir->second["3J_BLep_MVA"].fill( hyp.BLep()->CombinedMVA() );
            unmerged_dir->second["3J_BLep_CvsL"].fill( hyp.BLep()->CvsLtag() );
            unmerged_dir->second["3J_BLep_Mult"].fill( hyp.BLep()->numberOfDaughters() );

        } // end of unmerged_3J_evts




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

                /// 3 jet events

                if( njets_ == 3 ){ // Likelihood value computation

                    //                cout << "gen ttbar Pt: " << ttbar.Pt() << endl;

                    likelihood_3J_calc(event);

                    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue; // skip non-semilep events

                    // get matched perm from event
                    Permutation matched_perm = dr_matcher_.dr_match(
                            genp_selector_.ttbar_final_system(),
                            object_selector_.clean_jets(),
                            object_selector_.lepton(),
                            object_selector_.met(),
                            object_selector_.lepton_charge());

                    if( matched_perm.IsEmpty() ){ // skip to next event if perm is empty
                        permutator_.reset_3J();
                        continue;
                    }

                    matched_perm_3J_plots(event, matched_perm);
                    if( matched_perm.Merged_Event() ) merged_3J_evt(event, matched_perm);
                    else unmerged_3J_evt(event, matched_perm);

                    permutator_.reset_3J();
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
