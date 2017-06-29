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

class new_jet_perm_disc : public AnalyzerBase
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

            double maxmjet_cut_3J_ = 40.;
            double maxmjet_cut_4J_ = 50.;
//            double b_cut_ = 50.;
            double b_and_jet_upper_cut_ = 200.;
            double b_and_jet_lower_cut_ = 100.;


        //histograms
	    unordered_map<string, map< string, RObject> > histos_;
//	    map<TTNaming, string> > naming_;

	    CutFlowTracker tracker_; //tracks how many events pass cuts

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

			// 4 Jets
	    double Merged_BLep_mass_4J = -1.;
	    double Merged_BLep_comb_mass_4J = -1.;
	    double Max_Merged_BLep_mass_4J = -1.;

	    double Wrong_merged_wjet_mass_4J = -1.;
	    double Wrong_merged_wjet_comb_mass_4J = -1.;
	    double Max_Wrong_merged_wjet_mass_4J = -1.;

	    double Unmerged_BLep_mass_4J = -1.;
	    double Unmerged_BLep_comb_mass_4J = -1.;
	    double Max_Unmerged_BLep_mass_4J = -1.;

	    double Unmerged_BHad_mass_4J = -1.;
	    double Unmerged_BHad_comb_mass_4J = -1.;
	    double Max_Unmerged_BHad_mass_4J = -1.;

	    double Merged_BHad_mass_4J = -1.;
	    double Merged_BHad_comb_mass_4J = -1.;
	    double Max_Merged_BHad_mass_4J = -1.;

	    double Correct_merged_wjet_mass_4J = -1.;
	    double Correct_merged_wjet_comb_mass_4J = -1.;
	    double Max_Correct_merged_wjet_mass_4J = -1.;

            // Likelihood Discriminant
        double PermDiscr_3J = -1.;
        double PermDiscr_4J = -1.;


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

//		// 3 jets
//	    TTree *t1 = 0; // Unmerged_BLep_3J
//	    TTree *t2 = 0; // Merged_BLep_3J
//	    TTree *t3 = 0; // Unmerged_BHad_3J
//	    TTree *t4 = 0; // Merged_BHad_3J
//	    TTree *t5 = 0; // Correct_merged_wjet_3J
//	    TTree *t6 = 0; // Wrong_merged_wjet_3J
//
//		// 4 jets
//	    TTree *t7 = 0; // Merged_BLep_4J
//	    TTree *t8 = 0; // Wrong_merged_wjet_4J
//	    TTree *t9 = 0; // Unmerged_BLep_4J
//	    TTree *t10 = 0; // Unmerged_BHad_4J
//	    TTree *t11 = 0; // Merged_BHad_4J
//	    TTree *t12 = 0; // Correct_merged_wjet_4J
//
//        // Likelihood Discriminant
//        TTree *t13 = 0; // PermDiscr_3J
//        TTree *t14 = 0; // PermDiscr_4J

	public:
//        enum TTNaming {RIGHT, MERGE_SWAP, MERGE, WRONG};
	    new_jet_perm_disc(const std::string output_filename):
		AnalyzerBase("new_jet_perm_disc", output_filename),
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
//			// 3 jets
//		t1 = new TTree("Unmerged_BLep_3J", "Tree with Unmerged_BLep_3J.");
//		t1->Branch("Unmerged_BLep_mass_3J", &Unmerged_BLep_mass_3J, "Unmerged_BLep_mass_3J/D");
//		t1->Branch("Unmerged_BLep_comb_mass_3J", &Unmerged_BLep_comb_mass_3J, "Unmerged_BLep_comb_mass_3J/D");
//		t1->Branch("Max_Unmerged_BLep_mass_3J", &Max_Unmerged_BLep_mass_3J, "Max_Unmerged_BLep_mass_3J/D");
//
//		t2 = new TTree("Merged_BLep_3J", "Tree with Merged_BLep_3J.");
//		t2->Branch("Merged_BLep_mass_3J", &Merged_BLep_mass_3J, "Merged_BLep_mass_3J/D");
//		t2->Branch("Merged_BLep_comb_mass_3J", &Merged_BLep_comb_mass_3J, "Merged_BLep_comb_mass_3J/D");
//		t2->Branch("Max_Merged_BLep_mass_3J", &Max_Merged_BLep_mass_3J, "Max_Merged_BLep_mass_3J/D");
//
//		t3 = new TTree("Unmerged_BHad_3J", "Tree with Unmerged_BHad_3J.");
//		t3->Branch("Unmerged_BHad_mass_3J", &Unmerged_BHad_mass_3J, "Unmerged_BHad_mass_3J/D");
//		t3->Branch("Unmerged_BHad_comb_mass_3J", &Unmerged_BHad_comb_mass_3J, "Unmerged_BHad_comb_mass_3J/D");
//		t3->Branch("Max_Unmerged_BHad_mass_3J", &Max_Unmerged_BHad_mass_3J, "Max_Unmerged_BHad_mass_3J/D");
//
//		t4 = new TTree("Merged_BHad_3J", "Tree with Merged_BHad_3J.");
//		t4->Branch("Merged_BHad_mass_3J", &Merged_BHad_mass_3J, "Merged_BHad_mass_3J/D");
//		t4->Branch("Merged_BHad_comb_mass_3J", &Merged_BHad_comb_mass_3J, "Merged_BHad_comb_mass_3J/D");
//		t4->Branch("Max_Merged_BHad_mass_3J", &Max_Merged_BHad_mass_3J, "Max_Merged_BHad_mass_3J/D");
//
//		t5 = new TTree("Correct_merged_wjet_3J", "Tree with Correct_merged_wjet_3J.");
//		t5->Branch("Correct_merged_wjet_mass_3J", &Correct_merged_wjet_mass_3J, "Correct_merged_wjet_mass_3J/D");
//		t5->Branch("Correct_merged_wjet_comb_mass_3J", &Correct_merged_wjet_comb_mass_3J, "Correct_merged_wjet_comb_mass_3J/D");
//		t5->Branch("Max_Correct_merged_wjet_mass_3J", &Max_Correct_merged_wjet_mass_3J, "Max_Correct_merged_wjet_mass_3J/D");
//
//		t6 = new TTree("Wrong_merged_wjet_3J", "Tree with Wrong_merged_wjet_3J.");
//		t6->Branch("Wrong_merged_wjet_mass_3J", &Wrong_merged_wjet_mass_3J, "Wrong_merged_wjet_mass_3J/D");
//		t6->Branch("Wrong_merged_wjet_comb_mass_3J", &Wrong_merged_wjet_comb_mass_3J, "Wrong_merged_wjet_comb_mass_3J/D");
//		t6->Branch("Max_Wrong_merged_wjet_mass_3J", &Max_Wrong_merged_wjet_mass_3J, "Max_Wrong_merged_wjet_mass_3J/D");
//
//
//			// 4 jets
//		t7 = new TTree("Merged_BLep_4J", "Tree with Merged_BLep_4J.");
//		t7->Branch("Merged_BLep_mass_4J", &Merged_BLep_mass_4J, "Merged_BLep_mass_4J/D");
//		t7->Branch("Merged_BLep_comb_mass_4J", &Merged_BLep_comb_mass_4J, "Merged_BLep_comb_mass_4J/D");
//		t7->Branch("Max_Merged_BLep_mass_4J", &Max_Merged_BLep_mass_4J, "Max_Merged_BLep_mass_4J/D");
//
//		t8 = new TTree("Wrong_merged_wjet_4J", "Tree with Wrong_merged_wjet_4J.");
//		t8->Branch("Wrong_merged_wjet_mass_4J", &Wrong_merged_wjet_mass_4J, "Wrong_merged_wjet_mass_4J/D");
//		t8->Branch("Wrong_merged_wjet_comb_mass_4J", &Wrong_merged_wjet_comb_mass_4J, "Wrong_merged_wjet_comb_mass_4J/D");
//		t8->Branch("Max_Wrong_merged_wjet_mass_4J", &Max_Wrong_merged_wjet_mass_4J, "Max_Wrong_merged_wjet_mass_4J/D");
//
//		t9 = new TTree("Unmerged_BLep_4J", "Tree with Unmerged_BLep_4J.");
//		t9->Branch("Unmerged_BLep_mass_4J", &Unmerged_BLep_mass_4J, "Unmerged_BLep_mass_4J/D");
//		t9->Branch("Unmerged_BLep_comb_mass_4J", &Unmerged_BLep_comb_mass_4J, "Unmerged_BLep_comb_mass_4J/D");
//		t9->Branch("Max_Unmerged_BLep_mass_4J", &Max_Unmerged_BLep_mass_4J, "Max_Unmerged_BLep_mass_4J/D");
//
//		t10 = new TTree("Unmerged_BHad_4J", "Tree with Unmerged_BHad_4J.");
//		t10->Branch("Unmerged_BHad_mass_4J", &Unmerged_BHad_mass_4J, "Unmerged_BHad_mass_4J/D");
//		t10->Branch("Unmerged_BHad_comb_mass_4J", &Unmerged_BHad_comb_mass_4J, "Unmerged_BHad_comb_mass_4J/D");
//		t10->Branch("Max_Unmerged_BHad_mass_4J", &Max_Unmerged_BHad_mass_4J, "Max_Unmerged_BHad_mass_4J/D");
//
//		t11 = new TTree("Merged_BHad_4J", "Tree with Merged_BHad_4J.");
//		t11->Branch("Merged_BHad_mass_4J", &Merged_BHad_mass_4J, "Merged_BHad_mass_4J/D");
//		t11->Branch("Merged_BHad_comb_mass_4J", &Merged_BHad_comb_mass_4J, "Merged_BHad_comb_mass_4J/D");
//		t11->Branch("Max_Merged_BHad_mass_4J", &Max_Merged_BHad_mass_4J, "Max_Merged_BHad_mass_4J/D");
//
//		t12 = new TTree("Correct_merged_wjet_4J", "Tree with Correct_merged_wjet_4J.");
//		t12->Branch("Correct_merged_wjet_mass_4J", &Correct_merged_wjet_mass_4J, "Correct_merged_wjet_mass_4J/D");
//		t12->Branch("Correct_merged_wjet_comb_mass_4J", &Correct_merged_wjet_comb_mass_4J, "Correct_merged_wjet_comb_mass_4J/D");
//		t12->Branch("Max_Correct_merged_wjet_mass_4J", &Max_Correct_merged_wjet_mass_4J, "Max_Correct_merged_wjet_mass_4J/D");
//
//
//
//            // Likelihood Discriminant
//        t13 = new TTree("PermDiscr_3J", "Tree with likelihood discriminant 3J.");
//        t13->Branch("PermDiscr_3J", &PermDiscr_3J, "PermDiscr_3J/D");
//
//        t14 = new TTree("PermDiscr_4J", "Tree with likelihood discriminant 4J.");
//        t14->Branch("PermDiscr_4J", &PermDiscr_4J, "PermDiscr_4J/D");



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
        int permdisc_bins = 20;
        double permdisc_min = -10.;
        double permdisc_max = 10.;

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

        /// Permutation disc

            // 3 jets
//        book<TH1F>(likelihood, "3J_Likelihood_Ratio", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Preselection_Likelihood_Lowest_Merged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Preselection_Likelihood_Lowest_Unmerged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Preselection_Likelihood_Lowest_Correct", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Preselection_Likelihood_Lowest_Wrong", "", 100, -10., 10.);
//
//        book<TH1F>(likelihood, "3J_Square_Cut_Likelihood_Lowest_Merged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Square_Cut_Likelihood_Lowest_Unmerged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Square_Cut_Likelihood_Lowest_Correct", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "3J_Square_Cut_Likelihood_Lowest_Wrong", "", 100, -10., 10.);

        book<TH1F>(likelihood, "3J_Permdisc_All", "", permdisc_bins, permdisc_min, permdisc_max);
        book<TH1F>(likelihood, "3J_Permdisc_Lowest", "", permdisc_bins, permdisc_min, permdisc_max);
        book<TH1F>(likelihood, "3J_Permdisc_Best_Perm", "", permdisc_bins, permdisc_min, permdisc_max);

        book<TH1F>(likelihood, "3J_Permdisc_Best_Perm_RIGHT", "", permdisc_bins, permdisc_min, permdisc_max);
        book<TH1F>(likelihood, "3J_Permdisc_Best_Perm_MERGE_SWAP", "", permdisc_bins, permdisc_min, permdisc_max);
        book<TH1F>(likelihood, "3J_Permdisc_Best_Perm_MERGE", "", permdisc_bins, permdisc_min, permdisc_max);
        book<TH1F>(likelihood, "3J_Permdisc_Best_Perm_WRONG", "", permdisc_bins, permdisc_min, permdisc_max);

        /// Neutrino solver 

            // 3 jets
        book<TH1F>(likelihood, "3J_NSchi_All", "", nschi_bins, nschi_min, nschi_max);
        book<TH1F>(likelihood, "3J_NSchi_Lowest", "", nschi_bins, nschi_min, nschi_max);
        book<TH1F>(likelihood, "3J_NSchi_Best_Perm", "", nschi_bins, nschi_min, nschi_max);

        book<TH1F>(likelihood, "3J_NSchi_Best_Perm_RIGHT", "", nschi_bins, nschi_min, nschi_max);
        book<TH1F>(likelihood, "3J_NSchi_Best_Perm_MERGE_SWAP", "", nschi_bins, nschi_min, nschi_max);
        book<TH1F>(likelihood, "3J_NSchi_Best_Perm_MERGE", "", nschi_bins, nschi_min, nschi_max);
        book<TH1F>(likelihood, "3J_NSchi_Best_Perm_WRONG", "", nschi_bins, nschi_min, nschi_max);

        book<TH1F>(likelihood, "3J_NSdisc_All", "", ns_disc_bins, ns_disc_min, ns_disc_max);
        book<TH1F>(likelihood, "3J_NSdisc_Lowest", "", ns_disc_bins, ns_disc_min, ns_disc_max);
        book<TH1F>(likelihood, "3J_NSdisc_Best_Perm", "", ns_disc_bins, ns_disc_min, ns_disc_max);

        book<TH1F>(likelihood, "3J_NSdisc_Best_Perm_RIGHT", "", ns_disc_bins, ns_disc_min, ns_disc_max);
        book<TH1F>(likelihood, "3J_NSdisc_Best_Perm_MERGE_SWAP", "", ns_disc_bins, ns_disc_min, ns_disc_max);
        book<TH1F>(likelihood, "3J_NSdisc_Best_Perm_MERGE", "", ns_disc_bins, ns_disc_min, ns_disc_max);
        book<TH1F>(likelihood, "3J_NSdisc_Best_Perm_WRONG", "", ns_disc_bins, ns_disc_min, ns_disc_max);

        /// Total disc 

            // 3 jets
        book<TH1F>(likelihood, "3J_Totaldisc_All", "", comb_disc_bins, comb_disc_min, comb_disc_max);
        book<TH1F>(likelihood, "3J_Totaldisc_Lowest", "", comb_disc_bins, comb_disc_min, comb_disc_max);
        book<TH1F>(likelihood, "3J_Totaldisc_Best_Perm", "", comb_disc_bins, comb_disc_min, comb_disc_max);

        book<TH1F>(likelihood, "3J_Totaldisc_Best_Perm_RIGHT", "", comb_disc_bins, comb_disc_min, comb_disc_max);
        book<TH1F>(likelihood, "3J_Totaldisc_Best_Perm_MERGE_SWAP", "", comb_disc_bins, comb_disc_min, comb_disc_max);
        book<TH1F>(likelihood, "3J_Totaldisc_Best_Perm_MERGE", "", comb_disc_bins, comb_disc_min, comb_disc_max);
        book<TH1F>(likelihood, "3J_Totaldisc_Best_Perm_WRONG", "", comb_disc_bins, comb_disc_min, comb_disc_max);

            // 4 jets
//        book<TH1F>(likelihood, "4J_Likelihood_Ratio", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Preselection_Likelihood_Lowest_Merged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Preselection_Likelihood_Lowest_Unmerged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Preselection_Likelihood_Lowest_Correct", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Preselection_Likelihood_Lowest_Wrong", "", 100, -10., 10.);
//
//        book<TH1F>(likelihood, "4J_Square_Cut_Likelihood_Lowest_Merged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Square_Cut_Likelihood_Lowest_Unmerged", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Square_Cut_Likelihood_Lowest_Correct", "", 100, -10., 10.);
//        book<TH1F>(likelihood, "4J_Square_Cut_Likelihood_Lowest_Wrong", "", 100, -10., 10.);

//        book<TH1F>(likelihood, "4J_Likelihood_All", "", 20, -10., 10.);
//        book<TH1F>(likelihood, "4J_Likelihood_Lowest", "", 20, -10., 10.);

//        book<TH2D>(likelihood, "3J_vs_4J_Preselection_Likelihood_Lowest_Merged", "", 100, -10., 10., 100, -10., 10.);
//        book<TH2D>(likelihood, "3J_vs_4J_Preselection_Likelihood_Lowest_Unmerged", "", 100, -10., 10., 100, -10., 10.);
        


    // Merged jets discriminant plots
        string disc_right = "Correct_Disc_Plots";
        string disc_wrong = "Wrong_Disc_Plots";

            // 3 jets
        book<TH2D>(disc_right, "3J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
        book<TH2D>(disc_wrong, "3J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}

            // 4 jets
        book<TH2D>(disc_right, "4J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
        book<TH2D>(disc_wrong, "4J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}

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

//            // 4 jets
//                        // objects that pass preselection
//        book<TH1F>(delt, "4J_Preselection_Correct_Perms_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
//        book<TH1F>(delt, "4J_Preselection_Correct_Perms_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);
//        book<TH1F>(delt, "4J_Preselection_Wrong_Perms_Delta_mass", "" , mass_bins_, -2*mass_max_, 2*mass_max_);
//        book<TH1F>(delt, "4J_Preselection_Wrong_Perms_Delta_costh", "" , costh_bins_, -2*costh_max_, 2*costh_max_);
//
//                        // objects that pass the square cut
//        book<TH1F>(delt, "4J_Square_Cut_Correct_Perms_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
//        book<TH1F>(delt, "4J_Square_Cut_Correct_Perms_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);
//        book<TH1F>(delt, "4J_Square_Cut_Wrong_Perms_Delta_mass", "" , mass_bins_, -2*mass_max_, 2*mass_max_);
//        book<TH1F>(delt, "4J_Square_Cut_Wrong_Perms_Delta_costh", "" , costh_bins_, -2*costh_max_, 2*costh_max_);
//

		Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
    }

// events with 3 jets
    Permutation process_3J_evt(URStreamer &event){
        auto like_dir = histos_.find("Likelihood_Plots");
        tracker_.track("njets = 3");

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
        
        sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });

        if( num_btag < 2 ) return empty_perm; // require at least 2 b-tagged jets
        else{ 
            if( jets_vector[0]->BTagId(cut_medium_b_) == 0 ) swap(jets_vector[0],jets_vector[2]); // if (wj, bj, bj)
            else if( jets_vector[1]->BTagId(cut_medium_b_) == 0 ) swap(jets_vector[0],jets_vector[2]); // if (bj, wj, bj) 

            bj1 = jets_vector[0];
            bj2 = jets_vector[1];
            wj1 = jets_vector[2];
        }

        Permutation best_perm;
        double lowest_permdisc_3J = 1e10;
        double lowest_NSchi_3J = 1e10;
        double lowest_NSdisc_3J = 1e10;
        double lowest_Totaldisc_3J = 1e10;

        for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
            solver_.Solve_3J(test_perm);

            if( test_perm.PermDiscr() < lowest_permdisc_3J ){
                lowest_permdisc_3J = test_perm.PermDiscr();
//                best_perm = test_perm;
            }
            like_dir->second["3J_Permdisc_All"].fill(test_perm.PermDiscr()); // all perm disc values

            if( test_perm.NuChisq() < lowest_NSchi_3J && test_perm.NuChisq() > 0. ){
                lowest_NSchi_3J = test_perm.NuChisq();
//                cout << "NS chi: " << lowest_NSchi_3J << endl;
            }
            like_dir->second["3J_NSchi_All"].fill(test_perm.NuChisq()); // all neutrino solver chi2 values

            if( test_perm.NuDiscr() < lowest_NSdisc_3J ){
                lowest_NSdisc_3J = test_perm.NuDiscr();
//                cout << "NS disc: " << lowest_NSdisc_3J << endl;
            }
            like_dir->second["3J_NSdisc_All"].fill(test_perm.NuDiscr()); // all NS disc values

            if( test_perm.Prob() < lowest_Totaldisc_3J ){
                lowest_Totaldisc_3J = test_perm.Prob();
                best_perm = test_perm;
//                cout << "NS disc: " << lowest_NSdisc_3J << endl;
            }
            like_dir->second["3J_Totaldisc_All"].fill(test_perm.Prob()); // all total (perm + NS) disc values

        }

        // lowest values from event
        like_dir->second["3J_Permdisc_Lowest"].fill(lowest_permdisc_3J); // lowest perm disc value
        like_dir->second["3J_NSchi_Lowest"].fill(lowest_NSchi_3J); // lowest neutrino solver chi2 value
        like_dir->second["3J_NSdisc_Lowest"].fill(lowest_NSdisc_3J); // lowest neutrino solver disc value
        like_dir->second["3J_Totaldisc_Lowest"].fill(lowest_Totaldisc_3J); // lowest total (perm+NS) disc value


        if( !best_perm.IsEmpty() ){
            // values from best perm
            like_dir->second["3J_Permdisc_Best_Perm"].fill(best_perm.PermDiscr()); // best perm perm disc value
            like_dir->second["3J_NSchi_Best_Perm"].fill(best_perm.NuChisq()); // best perm neutrino solver chi2 value
            like_dir->second["3J_NSdisc_Best_Perm"].fill(best_perm.NuDiscr()); // best perm neutrino solver disc value
            like_dir->second["3J_Totaldisc_Best_Perm"].fill(best_perm.Prob()); // best perm total (perm+NS) disc value

//            cout << "   in function: " << lowest_NSchi_3J << endl;
        }

//        cout << "" << endl; 
        return best_perm;
    }

//// events with 4 jets
//    Permutation process_4J_evt(URStreamer &event){
////        auto like_dir = histos_.find("Likelihood_Plots");
//        tracker_.track("njets = 4");
//
//		//initialize permutation objects
//			//jets
//		IDJet* wj1 = 0;
//		IDJet* wj2 = 0;
//		IDJet* bj1 = 0;
//		IDJet* bj2 = 0;
//        Permutation empty_perm; // perm.WJa(), WJb(), BHad(), BLep()
////
//		vector<IDJet*> jets_vector;
//        int num_btag = 0;
////        int num_wjet = 0;
//    	for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
//    		jets_vector.push_back(*jets);
//
//            if( (*jets)->BTagId(cut_medium_b_) ) ++num_btag;
//        }
//        
//        sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });
//
//        if( num_btag < 2 ) return empty_perm; // require at least 2 b-tagged jets
//        else{ 
//            if( jets_vector[0]->BTagId(cut_medium_b_) == 0 && jets_vector[1]->BTagId(cut_medium_b_) == 0 ){ // if (wj, wj, bj, bj)
//                swap(jets_vector[0], jets_vector[2]);
//                swap(jets_vector[1], jets_vector[3]);
//            }
//            else if( jets_vector[0]->BTagId(cut_medium_b_) == 0 && jets_vector[2]->BTagId(cut_medium_b_) == 0 ){ // if (wj, bj, wj, bj)
//                swap(jets_vector[0], jets_vector[3]);
//            }
//            else if( jets_vector[0]->BTagId(cut_medium_b_) == 0 && jets_vector[3]->BTagId(cut_medium_b_) == 0 ){ // if (wj, bj, bj, wj)
//                swap(jets_vector[0], jets_vector[2]);
//            }
//            else if( jets_vector[1]->BTagId(cut_medium_b_) == 0 && jets_vector[2]->BTagId(cut_medium_b_) == 0 ){ // if (bj, wj, wj, bj)
//                swap(jets_vector[1], jets_vector[3]);
//            }
//            else if( jets_vector[1]->BTagId(cut_medium_b_) == 0 && jets_vector[3]->BTagId(cut_medium_b_) == 0 ){ // if (bj, wj, bj, wj)
//                swap(jets_vector[1], jets_vector[2]);
//            }
//
//            bj1 = jets_vector[0];
//            bj2 = jets_vector[1];
//            wj1 = ( jets_vector[2]->Pt() > jets_vector[3]->Pt() ) ? jets_vector[2] : jets_vector[3];
//            wj2 = ( jets_vector[2]->Pt() > jets_vector[3]->Pt() ) ? jets_vector[3] : jets_vector[2];
//        }
//
//        cout << "wj1 pt: " << wj1->Pt() << endl;
//        cout << "  wj2 pt: " << wj2->Pt() << endl;
//        if( !(bj1 && bj2) ) cout << "hi" << endl;
//
//        Permutation best_perm;
////        double lowest_permdisc_3J = 1e10;
////        for( auto test_perm : permutator_.permutations_3J(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge()) ){
////            solver_.Solve_3J(test_perm);
////
////            if( test_perm.PermDiscr() < lowest_permdisc_3J ){
////                lowest_permdisc_3J = test_perm.PermDiscr();
////                best_perm = test_perm;
////            }
////            like_dir->second["3J_Permdisc_All"].fill(test_perm.PermDiscr()); // all likelihood values
////
////
////        }
//////cout << "bp isEmpty: " << best_perm.IsEmpty() << endl;
////        if( !best_perm.IsEmpty() ){
////            like_dir->second["3J_Permdisc_Lowest"].fill(best_perm.PermDiscr()); // lowest event value
//////            cout << "bp disc: " << best_perm.PermDiscr() << endl;
////        }
//
//        return best_perm;
//    }



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
            auto like_dir = histos_.find("Likelihood_Plots");
//			auto delt_dir = histos_.find("Delta_Plots");
//            auto disc_correct_dir = histos_.find("Correct_Disc_Plots");
//            auto disc_wrong_dir = histos_.find("Wrong_Disc_Plots");
            auto gen_dir = histos_.find("Gen_Plots");

			//generator selection
			bool selection = genp_selector_.select(event);
            if( !selection ){
                    Logger::log().debug() << "event has no gen selection " << endl;
                    continue;
            }
            tracker_.track("gen selection");
			GenTTBar &ttbar = genp_selector_.ttbar_system();

//            if( ttbar.M() > 700 ) continue;
            gen_dir->second["Mttbar"].fill(ttbar.M());


            njets_ = 0;
            if( object_selector_.select(event) ) njets_ = object_selector_.clean_jets().size();
			if( njets_ < 3 ) continue;
//            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets_ << endl;
            tracker_.track("njet cuts");

        /// 3 jet events

            if( njets_ == 3 ){ // Likelihood value computation

//                cout << "Event number: " << evt_idx_ << endl;

            // get best perm from event
                Permutation best_perm = process_3J_evt(event);
                if( best_perm.IsEmpty() ){
                    permutator_.reset_3J();
                    continue; // skip to next event if perm is empty
                }

			    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                    like_dir->second["3J_Permdisc_Best_Perm_WRONG"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_NSchi_Best_Perm_WRONG"].fill(best_perm.NuChisq());
                    like_dir->second["3J_NSdisc_Best_Perm_WRONG"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Totaldisc_Best_Perm_WRONG"].fill(best_perm.Prob());
                    permutator_.reset_3J();
                    continue;
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
                    like_dir->second["3J_Permdisc_Best_Perm_WRONG"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_NSchi_Best_Perm_WRONG"].fill(best_perm.NuChisq());
                    like_dir->second["3J_NSdisc_Best_Perm_WRONG"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Totaldisc_Best_Perm_WRONG"].fill(best_perm.Prob());
                    permutator_.reset_3J();
                    continue;
                }
                tracker_.track("matched events");


                if( matched_perm.Merged_Event() ){ // matched perm has merged objets

            // objects in best_perm and matched_perm are the same and in same position
                    if( matched_perm.AreBsSame(best_perm) && ( best_perm.WJa() == matched_perm.WJa() || best_perm.WJa() == matched_perm.WJb() ) ){ 
                        like_dir->second["3J_Permdisc_Best_Perm_RIGHT"].fill(best_perm.PermDiscr());
                        like_dir->second["3J_NSchi_Best_Perm_RIGHT"].fill(best_perm.NuChisq());
                        like_dir->second["3J_NSdisc_Best_Perm_RIGHT"].fill(best_perm.NuDiscr());
                        like_dir->second["3J_Totaldisc_Best_Perm_RIGHT"].fill(best_perm.Prob());
                    }

            // objects in best_perm and matched_perm are the same but in differnt orders
                    else if( matched_perm.AreBsFlipped(best_perm) && ( best_perm.WJa() == matched_perm.WJa() || best_perm.WJa() == matched_perm.WJb() ) ){
                        like_dir->second["3J_Permdisc_Best_Perm_MERGE_SWAP"].fill(best_perm.PermDiscr());
                        like_dir->second["3J_NSchi_Best_Perm_MERGE_SWAP"].fill(best_perm.NuChisq());
                        like_dir->second["3J_NSdisc_Best_Perm_MERGE_SWAP"].fill(best_perm.NuDiscr());
                        like_dir->second["3J_Totaldisc_Best_Perm_MERGE_SWAP"].fill(best_perm.Prob());
                    }

           // objects in best_perm and matched_perm are not all the same
                    else{
                        like_dir->second["3J_Permdisc_Best_Perm_MERGE"].fill(best_perm.PermDiscr());
                        like_dir->second["3J_NSchi_Best_Perm_MERGE"].fill(best_perm.NuChisq());
                        like_dir->second["3J_NSdisc_Best_Perm_MERGE"].fill(best_perm.NuDiscr());
                        like_dir->second["3J_Totaldisc_Best_Perm_MERGE"].fill(best_perm.Prob());
                    }
                }

                else{ // events not merged
//                else if( !matched_perm.Merged_Event() ){
                    like_dir->second["3J_Permdisc_Best_Perm_WRONG"].fill(best_perm.PermDiscr());
                    like_dir->second["3J_NSchi_Best_Perm_WRONG"].fill(best_perm.NuChisq());
                    like_dir->second["3J_NSdisc_Best_Perm_WRONG"].fill(best_perm.NuDiscr());
                    like_dir->second["3J_Totaldisc_Best_Perm_WRONG"].fill(best_perm.Prob());
                } 

                permutator_.reset_3J();

            }


            if( njets_ == 3 ){ // matched perm plots

                auto disc_wrong_dir = histos_.find("Wrong_Disc_Plots");
                auto disc_correct_dir = histos_.find("Correct_Disc_Plots");
    			auto delt_dir = histos_.find("Delta_Plots");

                if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue; // skip non-semilep events

                GenTop* thad = ttbar.had_top();

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

                        // square cut
                        if( Max_Unmerged_BLep_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());
                        }

                            // bhad paired with wja

                        Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                        Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                        Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
//                        t3->Fill();
                        disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                        disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                        // preselection
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());

                        // square cut
                        if( Max_Unmerged_BHad_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());
                        }
                    }

                    if( matched_perm.WJb() ){ // if wjb exists
                            // blep paired with wjb

                        Unmerged_BLep_mass_3J = matched_perm.BLep()->M();
                        Unmerged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJb()).M();
                        Max_Unmerged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJb()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJb()->M();
  //                      t1->Fill();
                        disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);
                        disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                        // preselection
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());

                        // square cut
                        if( Max_Unmerged_BLep_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BLep()+*matched_perm.WJb()).M() < b_and_jet_upper_cut_ && (*matched_perm.BLep()+*matched_perm.WJb()).M() > b_and_jet_lower_cut_ ){
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());
                        }                             

                            // bhad paired with wjb

                        Unmerged_BHad_mass_3J = matched_perm.BHad()->M();
                        Unmerged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                        Max_Unmerged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
//                        t3->Fill();
                        disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);
                        disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                        // preselection
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());

                        // square cut
                        if( Max_Unmerged_BHad_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BHad()+*matched_perm.WJb()).M() < b_and_jet_upper_cut_ && (*matched_perm.BHad()+*matched_perm.WJb()).M() > b_and_jet_lower_cut_ ){
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                            delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());
                        }
                    }
                } // end of unmerged events

        // 2 jets merged together
                //wrong perm
                if( matched_perm.Merged_BLepWJa() && matched_perm.WJb() ){ // only blep and wja merged and wjb exists
                    Merged_BLep_mass_3J = matched_perm.BLep()->M(); 
                    Merged_BLep_comb_mass_3J =  (*matched_perm.BLep()+*matched_perm.WJb()).M();
                    Max_Merged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJb()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJb()->M();
//                    t2->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());

                    // square cut
                    if( Max_Merged_BLep_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BLep()+*matched_perm.WJb()).M() < b_and_jet_upper_cut_ && (*matched_perm.BLep()+*matched_perm.WJb()).M() > b_and_jet_lower_cut_ ){
                        delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJb()).M());
                        delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJb()).CosTheta());
                    }
                }

                    //wrong perm
                if( matched_perm.Merged_BLepWJb() && matched_perm.WJa() ){ // only blep and wjb merged and wja exists
                    Merged_BLep_mass_3J = matched_perm.BLep()->M(); 
                    Merged_BLep_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
                    Max_Merged_BLep_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
//                    t2->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);
                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                    delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());

                    // square cut
                    if( Max_Merged_BLep_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
                        delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
                        delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());
                    }
                }

                    // correct perm
                if( matched_perm.Merged_BHadWJa() && matched_perm.WJb() ){ // only bhad and wja merged and wjb exists
                    Merged_BHad_mass_3J = matched_perm.BHad()->M();
                    Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJb()).M();
                    Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJb()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJb()->M();
//                    t4->Fill();
                    disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                    disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                    delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());

                    // square cut
                    if( Max_Merged_BHad_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BHad()+*matched_perm.WJb()).M() < b_and_jet_upper_cut_ && (*matched_perm.BHad()+*matched_perm.WJb()).M() > b_and_jet_lower_cut_ ){
                        delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJb()).M());
                        delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJb()).CosTheta());
                    }
                }

                // correct perm
                if( matched_perm.Merged_BHadWJb() && matched_perm.WJa() ){ // only bhad and wjb merged and wja exists
                    Merged_BHad_mass_3J =  matched_perm.BHad()->M();
                    Merged_BHad_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
                    Max_Merged_BHad_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
//                    t4->Fill();
                    disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);
                    disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

                    // preselection
                    delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                    delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());

                    // square cut
                    if( Max_Merged_BHad_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
                        delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                        delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());
                    }
                }

				// merged wjets
				if( matched_perm.Merged_WJets() ){ // only wjets merged
				// correct perm
					Correct_merged_wjet_mass_3J = matched_perm.BHad()->M();
					Correct_merged_wjet_comb_mass_3J = (*matched_perm.BHad()+*matched_perm.WJa()).M();
					Max_Correct_merged_wjet_mass_3J = ( matched_perm.BHad()->M() > matched_perm.WJa()->M() ) ? matched_perm.BHad()->M() : matched_perm.WJa()->M();
//					t5->Fill();
                    disc_correct_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Correct_merged_wjet_mass_3J, Correct_merged_wjet_comb_mass_3J);
                    disc_correct_dir->second["3J_nusolver_chi2_correct"].fill(matched_perm.NuChisq());

                   // preselection
                   delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                   delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());

                   // square cut
                   if( Max_Correct_merged_wjet_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BHad()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
                       delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BHad()+*matched_perm.WJa()).M());
                       delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BHad()+*matched_perm.WJa()).CosTheta());
                   }

				//wrong perm
					Wrong_merged_wjet_mass_3J = matched_perm.BLep()->M();
					Wrong_merged_wjet_comb_mass_3J = (*matched_perm.BLep()+*matched_perm.WJa()).M();
					Max_Wrong_merged_wjet_mass_3J = ( matched_perm.BLep()->M() > matched_perm.WJa()->M() ) ? matched_perm.BLep()->M() : matched_perm.WJa()->M();
//					t6->Fill();
                    disc_wrong_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Wrong_merged_wjet_mass_3J, Wrong_merged_wjet_comb_mass_3J);
//                    disc_wrong_dir->second["3J_nusolver_chi2_wrong"].fill(matched_perm.NuChisq());

					// preselection
					delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
					delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());
	
					// square cut
					if( Max_Wrong_merged_wjet_mass_3J > maxmjet_cut_3J_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() < b_and_jet_upper_cut_ && (*matched_perm.BLep()+*matched_perm.WJa()).M() > b_and_jet_lower_cut_ ){
					    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*matched_perm.BLep()+*matched_perm.WJa()).M());
					    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*matched_perm.BLep()+*matched_perm.WJa()).CosTheta());
					}

				} // end of merged W jets loop

            } // end of 3 jets


//        /// 4 jet events
//            if( njets_ == 4 ) process_4J_evt(event);

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
	URDriver<new_jet_perm_disc> test;
	int thing = test.run();
	return thing;
}
