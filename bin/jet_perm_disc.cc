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
#include "Analyses/URTTbar/interface/TTPermutator.h"
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
#include "Analyses/URTTbar/interface/TTBarPlots.h"
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

            double b_cut_ = 50.;
            double b_and_jet_upper_cut_ = 200.;
            double b_and_jet_lower_cut_ = 150.;


        //histograms
	    unordered_map<string, map< string, RObject> > histos_;
	    //map<TTNaming, string> naming_;
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



	//switches
	    bool isData_, isTTbar_;

	//selectors and helpers
            TTObjectSelector object_selector_; //selects ttbar objects
            TTGenParticleSelector genp_selector_; //selects generator level objects
	    TTBarSolver solver_; //solves ttbar events
	    TTGenMatcher matcher_; //matches particles on generator level
            TTPermutator permutator_;

	    float evt_weight_;
	    MCWeightProducer mc_weights_;
	    IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
	    IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
	    IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;

		// 3 jets
	    TTree *t1 = 0; // Unmerged_BLep_3J
	    TTree *t2 = 0; // Merged_BLep_3J
	    TTree *t3 = 0; // Unmerged_BHad_3J
	    TTree *t4 = 0; // Merged_BHad_3J
	    TTree *t5 = 0; // Correct_merged_wjet_3J
	    TTree *t6 = 0; // Wrong_merged_wjet_3J

		// 4 jets
	    TTree *t7 = 0; // Merged_BLep_4J
	    TTree *t8 = 0; // Wrong_merged_wjet_4J
	    TTree *t9 = 0; // Unmerged_BLep_4J
	    TTree *t10 = 0; // Unmerged_BHad_4J
	    TTree *t11 = 0; // Merged_BHad_4J
	    TTree *t12 = 0; // Correct_merged_wjet_4J


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

			//Init Solver
			string filename = "prob_ttJets.permProbComputer.test.root";
                        Logger::log().debug() << "solver file: " << filename << endl;
                        TFile probfile(DataFile(filename).path().c_str());
                        TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
			solver_.Init(td, false, true, true, true, true); //probfile, btag, nusolv,massdis,anghad,anglep

//			cut_tight_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.tightb"));
//			cut_loose_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.looseb"));
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
			// 3 jets
		t1 = new TTree("Unmerged_BLep_3J", "Tree with Unmerged_BLep_3J.");
		t1->Branch("Unmerged_BLep_mass_3J", &Unmerged_BLep_mass_3J, "Unmerged_BLep_mass_3J/D");
		t1->Branch("Unmerged_BLep_comb_mass_3J", &Unmerged_BLep_comb_mass_3J, "Unmerged_BLep_comb_mass_3J/D");
		t1->Branch("Max_Unmerged_BLep_mass_3J", &Max_Unmerged_BLep_mass_3J, "Max_Unmerged_BLep_mass_3J/D");

		t2 = new TTree("Merged_BLep_3J", "Tree with Merged_BLep_3J.");
		t2->Branch("Merged_BLep_mass_3J", &Merged_BLep_mass_3J, "Merged_BLep_mass_3J/D");
		t2->Branch("Merged_BLep_comb_mass_3J", &Merged_BLep_comb_mass_3J, "Merged_BLep_comb_mass_3J/D");
		t2->Branch("Max_Merged_BLep_mass_3J", &Max_Merged_BLep_mass_3J, "Max_Merged_BLep_mass_3J/D");

		t3 = new TTree("Unmerged_BHad_3J", "Tree with Unmerged_BHad_3J.");
		t3->Branch("Unmerged_BHad_mass_3J", &Unmerged_BHad_mass_3J, "Unmerged_BHad_mass_3J/D");
		t3->Branch("Unmerged_BHad_comb_mass_3J", &Unmerged_BHad_comb_mass_3J, "Unmerged_BHad_comb_mass_3J/D");
		t3->Branch("Max_Unmerged_BHad_mass_3J", &Max_Unmerged_BHad_mass_3J, "Max_Unmerged_BHad_mass_3J/D");

		t4 = new TTree("Merged_BHad_3J", "Tree with Merged_BHad_3J.");
		t4->Branch("Merged_BHad_mass_3J", &Merged_BHad_mass_3J, "Merged_BHad_mass_3J/D");
		t4->Branch("Merged_BHad_comb_mass_3J", &Merged_BHad_comb_mass_3J, "Merged_BHad_comb_mass_3J/D");
		t4->Branch("Max_Merged_BHad_mass_3J", &Max_Merged_BHad_mass_3J, "Max_Merged_BHad_mass_3J/D");

		t5 = new TTree("Correct_merged_wjet_3J", "Tree with Correct_merged_wjet_3J.");
		t5->Branch("Correct_merged_wjet_mass_3J", &Correct_merged_wjet_mass_3J, "Correct_merged_wjet_mass_3J/D");
		t5->Branch("Correct_merged_wjet_comb_mass_3J", &Correct_merged_wjet_comb_mass_3J, "Correct_merged_wjet_comb_mass_3J/D");
		t5->Branch("Max_Correct_merged_wjet_mass_3J", &Max_Correct_merged_wjet_mass_3J, "Max_Correct_merged_wjet_mass_3J/D");

		t6 = new TTree("Wrong_merged_wjet_3J", "Tree with Wrong_merged_wjet_3J.");
		t6->Branch("Wrong_merged_wjet_mass_3J", &Wrong_merged_wjet_mass_3J, "Wrong_merged_wjet_mass_3J/D");
		t6->Branch("Wrong_merged_wjet_comb_mass_3J", &Wrong_merged_wjet_comb_mass_3J, "Wrong_merged_wjet_comb_mass_3J/D");
		t6->Branch("Max_Wrong_merged_wjet_mass_3J", &Max_Wrong_merged_wjet_mass_3J, "Max_Wrong_merged_wjet_mass_3J/D");


			// 4 jets
		t7 = new TTree("Merged_BLep_4J", "Tree with Merged_BLep_4J.");
		t7->Branch("Merged_BLep_mass_4J", &Merged_BLep_mass_4J, "Merged_BLep_mass_4J/D");
		t7->Branch("Merged_BLep_comb_mass_4J", &Merged_BLep_comb_mass_4J, "Merged_BLep_comb_mass_4J/D");
		t7->Branch("Max_Merged_BLep_mass_4J", &Max_Merged_BLep_mass_4J, "Max_Merged_BLep_mass_4J/D");

		t8 = new TTree("Wrong_merged_wjet_4J", "Tree with Wrong_merged_wjet_4J.");
		t8->Branch("Wrong_merged_wjet_mass_4J", &Wrong_merged_wjet_mass_4J, "Wrong_merged_wjet_mass_4J/D");
		t8->Branch("Wrong_merged_wjet_comb_mass_4J", &Wrong_merged_wjet_comb_mass_4J, "Wrong_merged_wjet_comb_mass_4J/D");
		t8->Branch("Max_Wrong_merged_wjet_mass_4J", &Max_Wrong_merged_wjet_mass_4J, "Max_Wrong_merged_wjet_mass_4J/D");

		t9 = new TTree("Unmerged_BLep_4J", "Tree with Unmerged_BLep_4J.");
		t9->Branch("Unmerged_BLep_mass_4J", &Unmerged_BLep_mass_4J, "Unmerged_BLep_mass_4J/D");
		t9->Branch("Unmerged_BLep_comb_mass_4J", &Unmerged_BLep_comb_mass_4J, "Unmerged_BLep_comb_mass_4J/D");
		t9->Branch("Max_Unmerged_BLep_mass_4J", &Max_Unmerged_BLep_mass_4J, "Max_Unmerged_BLep_mass_4J/D");

		t10 = new TTree("Unmerged_BHad_4J", "Tree with Unmerged_BHad_4J.");
		t10->Branch("Unmerged_BHad_mass_4J", &Unmerged_BHad_mass_4J, "Unmerged_BHad_mass_4J/D");
		t10->Branch("Unmerged_BHad_comb_mass_4J", &Unmerged_BHad_comb_mass_4J, "Unmerged_BHad_comb_mass_4J/D");
		t10->Branch("Max_Unmerged_BHad_mass_4J", &Max_Unmerged_BHad_mass_4J, "Max_Unmerged_BHad_mass_4J/D");

		t11 = new TTree("Merged_BHad_4J", "Tree with Merged_BHad_4J.");
		t11->Branch("Merged_BHad_mass_4J", &Merged_BHad_mass_4J, "Merged_BHad_mass_4J/D");
		t11->Branch("Merged_BHad_comb_mass_4J", &Merged_BHad_comb_mass_4J, "Merged_BHad_comb_mass_4J/D");
		t11->Branch("Max_Merged_BHad_mass_4J", &Max_Merged_BHad_mass_4J, "Max_Merged_BHad_mass_4J/D");

		t12 = new TTree("Correct_merged_wjet_4J", "Tree with Correct_merged_wjet_4J.");
		t12->Branch("Correct_merged_wjet_mass_4J", &Correct_merged_wjet_mass_4J, "Correct_merged_wjet_mass_4J/D");
		t12->Branch("Correct_merged_wjet_comb_mass_4J", &Correct_merged_wjet_comb_mass_4J, "Correct_merged_wjet_comb_mass_4J/D");
		t12->Branch("Max_Correct_merged_wjet_mass_4J", &Max_Correct_merged_wjet_mass_4J, "Max_Correct_merged_wjet_mass_4J/D");



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


    // Merged jets discriminant plots
        string disc = "Disc_Plots";

            // 3 jets
        book<TH2D>(disc, "3J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
        book<TH2D>(disc, "3J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}

            // 4 jets
        book<TH2D>(disc, "4J_mbpjet_vs_maxmjet_correct", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}
        book<TH2D>(disc, "4J_mbpjet_vs_maxmjet_wrong", "", 100, 0., 200., 200, 0., 2000.); // hist for m_{b+jet} vs max m_{jet}


	// Delta (gen-reco) Plots
        string delt = "Delta_Plots";

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


	// Hadronic side event categories
        for( int i = 0; i < 2; i++ ){
                        // 3 jets
                book<TH1F>(DRnames_[i], "3J_Had_Resolved_vs_mttbar", "", nbins, mass_min_, mass_max_); // 3 hadronic jets resolved
                book<TH1F>(DRnames_[i], "3J_Merged_BHadWJet_vs_mttbar", "", nbins, mass_min_, mass_max_); // BHad and WJa or WJb merged
                book<TH1F>(DRnames_[i], "3J_Merged_WJets_vs_mttbar", "", nbins, mass_min_, mass_max_); // jets from W merged
                book<TH1F>(DRnames_[i], "3J_Merged_THad_vs_mttbar", "", nbins, mass_min_, mass_max_); // all jets on hadronic side merged
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_vs_mttbar", "", nbins, mass_min_, mass_max_); // 1 or more partons not matched

                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_BHad_mass", "", mass_bins_, 0, 200);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_BHad_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_BHad_eta", "", eta_bins_, -3, 3);

                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJa_mass", "", mass_bins_, 0, 200);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJa_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJa_eta", "", eta_bins_, -3, 3);

                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJb_mass", "", mass_bins_, 0, 200);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJb_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "3J_Non_Reconstructable_WJb_eta", "", eta_bins_, -3, 3);

        				// 4 jets
                book<TH1F>(DRnames_[i], "4J_Had_Resolved_vs_mttbar", "", nbins, mass_min_, mass_max_); // 3 hadronic jets resolved
                book<TH1F>(DRnames_[i], "4J_Merged_BHadWJet_vs_mttbar", "", nbins, mass_min_, mass_max_); // BHad and WJa or WJb merged
                book<TH1F>(DRnames_[i], "4J_Merged_WJets_vs_mttbar", "", nbins, mass_min_, mass_max_); // jets from W merged
                book<TH1F>(DRnames_[i], "4J_Merged_THad_vs_mttbar", "", nbins, mass_min_, mass_max_); // all jets on hadronic side merged
                book<TH1F>(DRnames_[i], "4J_Non_Reconstructable_vs_mttbar", "", nbins, mass_min_, mass_max_); // 1 or more partons not matched

                        // 5 or more jets
                book<TH1F>(DRnames_[i], "5PJ_Had_Resolved_vs_mttbar", "", nbins, mass_min_, mass_max_); // 3 hadronic jets resolved
                book<TH1F>(DRnames_[i], "5PJ_Merged_BHadWJet_vs_mttbar", "", nbins, mass_min_, mass_max_); // BHad and WJa or WJb merged
                book<TH1F>(DRnames_[i], "5PJ_Merged_WJets_vs_mttbar", "", nbins, mass_min_, mass_max_); // jets from W merged
                book<TH1F>(DRnames_[i], "5PJ_Merged_THad_vs_mttbar", "", nbins, mass_min_, mass_max_); // all jets on hadronic side merged
                book<TH1F>(DRnames_[i], "5PJ_Non_Reconstructable_vs_mttbar", "", nbins, mass_min_, mass_max_); // 1 or more partons not matched

        }

		Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
    }

        //This method is called once every file, contains the event loop
        ///run your proper analysis here
	virtual void analyze()
	{
		Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

        URStreamer event(tree_);

		while(event.next() /*&& evt_idx_ < 30*/)
        {
            evt_idx_++;
            if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

			auto delt_dir = histos_.find("Delta_Plots");
            auto disc_dir = histos_.find("Disc_Plots");

			//generator selection
			bool selection = genp_selector_.select(event);
            tracker_.track("gen selection");
            if( !selection ){
                    Logger::log().debug() << "event has no selection " << endl;
                    continue;
            }
			GenTTBar &ttbar = genp_selector_.ttbar_system();

            Permutation matched_perm;
            if( object_selector_.select(event) ){
                    if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
                            matched_perm = matcher_.match(
                                    genp_selector_.ttbar_final_system(),
                                    object_selector_.clean_jets(),
                                    object_selector_.loose_electrons(),
                                    object_selector_.loose_muons());
                    }
                    matched_perm.SetMET(object_selector_.met());
            }

			if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue;

			njets_ = object_selector_.clean_jets().size();

			if( njets_ < 3 ) continue;


			//initialize gen jets
//			const GenObject* lepton = 0;
			const GenObject* BLep = 0;
			const GenObject* BHad = 0;
			const GenObject* WJa = 0;
			const GenObject* WJb = 0;
//			const GenTop* tlep = 0;
			const GenTop* thad = 0;

			//Generator level definitions
//			if( ttbar.lep_top() ) tlep = ttbar.lep_top();
			if( ttbar.had_top() ) thad = ttbar.had_top();
//			if( ttbar.lepton() ) lepton = ttbar.lepton();
			if( ttbar.lep_b() ) BLep = ttbar.lep_b();
			if( ttbar.had_b() ) BHad = ttbar.had_b();
			if( ttbar.had_W() ){
				if( ttbar.had_W()->first && !ttbar.had_W()->second ) WJa = ttbar.had_W()->first;
				if( ttbar.had_W()->second && !ttbar.had_W()->first ) WJb = ttbar.had_W()->second;
				if( ttbar.had_W()->first && ttbar.had_W()->second ){
					WJa = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
					WJb = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
				}
			}

			if( !thad ) continue;


		//////////////// Perm Matching
			for( int i = 0; i < 1; i++ ){

				auto reco_dir = histos_.find(DRnames_[i]);

				//initialize permutation objects
					//jets
				const IDJet* best_BHad = 0;
				const IDJet* best_BLep = 0;
				const IDJet* best_WJa = 0;
				const IDJet* best_WJb = 0;


				list<IDJet*> bhad_list;
				float bhad_DR = 1e10;// itialize dr to high number
				for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
					if( !(BHad == 0) ){
						if( (*jets)->DeltaR(*BHad) > DR_[i] ) continue;
						if( (*jets)->DeltaR(*BHad) < bhad_DR ){
							bhad_DR = (*jets)->DeltaR(*BHad);
							bhad_list.push_back(*jets);
							best_BHad = (bhad_list.back());
						}
					}
				}
				list<IDJet*> blep_list;
				float blep_DR = 1e10;// itialize dr to high number
				for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
					if( !(BLep == 0) ){
						if( (*jets)->DeltaR(*BLep) > DR_[i] ) continue;
						if( (*jets)->DeltaR(*BLep) < blep_DR ){
							blep_DR = (*jets)->DeltaR(*BLep);
							blep_list.push_back(*jets);
							best_BLep = (blep_list.back());
						}
					}
				}
				list<IDJet*> wja_list;
				float wja_DR = 1e10;// itialize dr to high number
				for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
					if( !(WJa == 0) ){
						if( (*jets)->DeltaR(*WJa) > DR_[i] ) continue;
						if( (*jets)->DeltaR(*WJa) < wja_DR ){
							wja_DR = (*jets)->DeltaR(*WJa);
							wja_list.push_back(*jets);
							best_WJa = (wja_list.back());
						}
					}
				}
				list<IDJet*> wjb_list;
				float wjb_DR = 1e10;// itialize dr to high number
				for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
					if( !(WJb == 0) ){
						if( (*jets)->DeltaR(*WJb) > DR_[i] ) continue;
						if( (*jets)->DeltaR(*WJb) < wjb_DR ){
							wjb_DR = (*jets)->DeltaR(*WJb);
							wjb_list.push_back(*jets);
							best_WJb = (wjb_list.back());
						}
					}
				}

		///////// fill trees for discriminant and preselection/signal cuts
        
   /// 3 jets
				if( njets_ == 3 && i == 0 ){ // 3 jets in events and DR = 0.4
					if( !(best_BHad && best_BLep) ) continue; // require both b's
					if( !(best_WJa || best_WJb) ) continue; // require one wjet
					if( best_BHad == best_BLep ) continue; // both b's can't be merged
					if( (best_BHad == best_WJa && best_BHad == best_WJb) || (best_BLep == best_WJa && best_BLep == best_WJb) ) continue; // three jets can't be merged

            // 0 jets merged together
                    //all wrong perms
                    if( !(best_BLep == best_WJa) && !(best_BLep == best_WJb) && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_WJa == best_WJb) ){
                        if( best_WJa ){ // if wja exists
                                // blep paired with wja
                            Unmerged_BLep_mass_3J = best_BLep->M();
                            Unmerged_BLep_comb_mass_3J = (*best_BLep+*best_WJa).M();
                            Max_Unmerged_BLep_mass_3J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
                            t1->Fill();
                            disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);

                            // preselection
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());

                            // square cut
                            if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJa).M() > b_and_jet_lower_cut_ ){
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
                            }

                                // bhad paired with wja
                            Unmerged_BHad_mass_3J = best_BHad->M();
                            Unmerged_BHad_comb_mass_3J = (*best_BHad+*best_WJa).M();
                            Max_Unmerged_BHad_mass_3J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
                            t3->Fill();
                            disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);

                            // preselection
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());

                            // square cut
                            if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJa).M() > b_and_jet_lower_cut_ ){
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
                            }
                        }

                        if( best_WJb ){ // if wjb exists
                                // blep paired with wjb
                            Unmerged_BLep_mass_3J = best_BLep->M();
                            Unmerged_BLep_comb_mass_3J = (*best_BLep+*best_WJb).M();
                            Max_Unmerged_BLep_mass_3J = ( best_BLep->M() > best_WJb->M() ) ? best_BLep->M() : best_WJb->M();
                            t1->Fill();
                            disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_3J, Unmerged_BLep_comb_mass_3J);

                            // preselection
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJb).M());
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJb).CosTheta());

                            // square cut
                            if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJb).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJb).M() > b_and_jet_lower_cut_ ){
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJb).M());
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJb).CosTheta());
                            }

                                // bhad paired with wjb
                            Unmerged_BHad_mass_3J = best_BHad->M();
                            Unmerged_BHad_comb_mass_3J = (*best_BHad+*best_WJb).M();
                            Max_Unmerged_BHad_mass_3J = ( best_BHad->M() > best_WJb->M() ) ? best_BHad->M() : best_WJb->M();
                            t3->Fill();
                            disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_3J, Unmerged_BHad_comb_mass_3J);

                            // preselection
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJb).M());
                            delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJb).CosTheta());

                            // square cut
                            if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJb).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJb).M() > b_and_jet_lower_cut_ ){
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJb).M());
                                    delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJb).CosTheta());
                            }
                        }
                    }

            // 2 jets merged together
                    //wrong perm
                    if( (best_BLep == best_WJa) && !(best_BLep == best_WJb) && !(best_BHad == best_WJb) && best_WJb ){ // only blep and wja merged and wjb exists
                        Merged_BLep_mass_3J = best_BLep->M(); 
                        Merged_BLep_comb_mass_3J =  (*best_BLep+*best_WJb).M();
                        Max_Merged_BLep_mass_3J = ( best_BLep->M() > best_WJb->M() ) ? best_BLep->M() : best_WJb->M();
                        t2->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);

                        // preselection
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJb).M());
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJb).CosTheta());

                        // square cut
                        if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJb).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJb).M() > b_and_jet_lower_cut_ ){
                                delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJb).M());
                                delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJb).CosTheta());
                        }
                    }

                    //wrong perm
                    if( (best_BLep == best_WJb) && !(best_BLep == best_WJa) && !(best_BHad == best_WJa) && best_WJa ){ // only blep and wjb merged and wja exists
                        Merged_BLep_mass_3J = best_BLep->M(); 
                        Merged_BLep_comb_mass_3J = (*best_BLep+*best_WJa).M();
                        Max_Merged_BLep_mass_3J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
                        t2->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_3J, Merged_BLep_comb_mass_3J);

                        // preselection
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
                        delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());

                        // square cut
                        if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJa).M() > b_and_jet_lower_cut_ ){
                                delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
                                delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
                        }
                    }

                    // correct perm
                    if( (best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BLep == best_WJb) && best_WJb ){ // only bhad and wja merged and wjb exists
                        Merged_BHad_mass_3J = best_BHad->M();
                        Merged_BHad_comb_mass_3J = (*best_BHad+*best_WJb).M();
                        Max_Merged_BHad_mass_3J = ( best_BHad->M() > best_WJb->M() ) ? best_BHad->M() : best_WJb->M();
                        t4->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);

                        // preselection
                        delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJb).M());
                        delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJb).CosTheta());

                        // square cut
                        if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJb).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJb).M() > b_and_jet_lower_cut_ ){
                                delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJb).M());
                                delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJb).CosTheta());
                        }
                    }

                    // correct perm
                    if( (best_BHad == best_WJb) && !(best_BHad == best_WJa) && !(best_BLep == best_WJa) && best_WJa ){ // only bhad and wjb merged and wja exists
                        Merged_BHad_mass_3J =  best_BHad->M();
                        Merged_BHad_comb_mass_3J = (*best_BHad+*best_WJa).M();
                        Max_Merged_BHad_mass_3J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
                        t4->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_3J, Merged_BHad_comb_mass_3J);

                        // preselection
                        delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                        delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());

                        // square cut
                        if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJa).M() > b_and_jet_lower_cut_ ){
                                delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                                delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
                        }
                    }

					// merged wjets
					if( (best_WJa == best_WJb) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){ // only wjets merged
					// correct perm
						Correct_merged_wjet_mass_3J = best_BHad->M();
						Correct_merged_wjet_comb_mass_3J = (*best_BHad+*best_WJa).M();
						Max_Correct_merged_wjet_mass_3J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
						t5->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_correct"].fill(Max_Correct_merged_wjet_mass_3J, Correct_merged_wjet_comb_mass_3J);

                       // preselection
                       delt_dir->second["3J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                       delt_dir->second["3J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());

                       // square cut
                       if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJa).M() > b_and_jet_lower_cut_ ){
                               delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                               delt_dir->second["3J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
                       }

					//wrong perm
						Wrong_merged_wjet_mass_3J = best_BLep->M();
						Wrong_merged_wjet_comb_mass_3J = (*best_BLep+*best_WJa).M();
						Max_Wrong_merged_wjet_mass_3J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
						t6->Fill();
                        disc_dir->second["3J_mbpjet_vs_maxmjet_wrong"].fill(Max_Wrong_merged_wjet_mass_3J, Wrong_merged_wjet_comb_mass_3J);

						// preselection
						delt_dir->second["3J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
						delt_dir->second["3J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
	
						// square cut
						if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJa).M() > b_and_jet_lower_cut_ ){
							delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
							delt_dir->second["3J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
						}

					} // end of merged W jets loop

				} // end of njets = 3 loop for disc and presel/signal



    /// 4 jets
				if( njets_ == 4 && i == 0 ){ // 4 jets in events and DR = 0.4
					if( !(best_BHad && best_BLep) ) continue; // require both b's
					if( !(best_WJa || best_WJb) ) continue; // require one wjet
					if( best_BHad == best_BLep ) continue; // both b's can't be merged
					if( (best_BHad == best_WJa && best_BHad == best_WJb) || (best_BLep == best_WJa && best_BLep == best_WJb) ) continue; // three jets can't be merged

            // 0 jets merged together
                    //all wrong perms
                    if( !(best_BLep == best_WJa) && !(best_BLep == best_WJb) && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_WJa == best_WJb) ){
                        if( best_WJa ){ // if wja exists
                                // blep paired with wja
                            Unmerged_BLep_mass_4J = best_BLep->M();
                            Unmerged_BLep_comb_mass_4J = (*best_BLep+*best_WJa).M();
                            Max_Unmerged_BLep_mass_4J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
                            t9->Fill();
                            disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_4J, Unmerged_BLep_comb_mass_4J);

                                // bhad paired with wja
                            Unmerged_BHad_mass_4J = best_BHad->M();
                            Unmerged_BHad_comb_mass_4J = (*best_BHad+*best_WJa).M();
                            Max_Unmerged_BHad_mass_4J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
                            t10->Fill();
                            disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_4J, Unmerged_BHad_comb_mass_4J);
                        }

                        if( best_WJb ){ // if wjb exists
                                // blep paired with wjb
                            Unmerged_BLep_mass_4J = best_BLep->M();
                            Unmerged_BLep_comb_mass_4J = (*best_BLep+*best_WJb).M();
                            Max_Unmerged_BLep_mass_4J = ( best_BLep->M() > best_WJb->M() ) ? best_BLep->M() : best_WJb->M();
                            t9->Fill();
                            disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BLep_mass_4J, Unmerged_BLep_comb_mass_4J);

                                // bhad paired with wjb
                            Unmerged_BHad_mass_4J = best_BHad->M();
                            Unmerged_BHad_comb_mass_4J = (*best_BHad+*best_WJb).M();
                            Max_Unmerged_BHad_mass_4J = ( best_BHad->M() > best_WJb->M() ) ? best_BHad->M() : best_WJb->M();
                            t10->Fill();
                            disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Unmerged_BHad_mass_4J, Unmerged_BHad_comb_mass_4J);
                        }
                    }

            // 2 jets merged together
                    //wrong perm
                    if( (best_BLep == best_WJa) && !(best_BLep == best_WJb) && !(best_BHad == best_WJb) && best_WJb ){ // only blep and wja merged and wjb exists
                        Merged_BLep_mass_4J = best_BLep->M(); 
                        Merged_BLep_comb_mass_4J =  (*best_BLep+*best_WJb).M();
                        Max_Merged_BLep_mass_4J = ( best_BLep->M() > best_WJb->M() ) ? best_BLep->M() : best_WJb->M();
                        t7->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_4J, Merged_BLep_comb_mass_4J);
                    }

                    //wrong perm
                    if( (best_BLep == best_WJb) && !(best_BLep == best_WJa) && !(best_BHad == best_WJa) && best_WJa ){ // only blep and wjb merged and wja exists
                        Merged_BLep_mass_4J = best_BLep->M(); 
                        Merged_BLep_comb_mass_4J = (*best_BLep+*best_WJa).M();
                        Max_Merged_BLep_mass_4J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
                        t7->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Merged_BLep_mass_4J, Merged_BLep_comb_mass_4J);
                    }

                    // correct perm
                    if( (best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BLep == best_WJb) && best_WJb ){ // only bhad and wja merged and wjb exists
                        Merged_BHad_mass_4J = best_BHad->M();
                        Merged_BHad_comb_mass_4J = (*best_BHad+*best_WJb).M();
                        Max_Merged_BHad_mass_4J = ( best_BHad->M() > best_WJb->M() ) ? best_BHad->M() : best_WJb->M();
                        t11->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_4J, Merged_BHad_comb_mass_4J);
                    }

                    // correct perm
                    if( (best_BHad == best_WJb) && !(best_BHad == best_WJa) && !(best_BLep == best_WJa) && best_WJa ){ // only bhad and wjb merged and wja exists
                        Merged_BHad_mass_4J =  best_BHad->M();
                        Merged_BHad_comb_mass_4J = (*best_BHad+*best_WJa).M();
                        Max_Merged_BHad_mass_4J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
                        t11->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_correct"].fill(Max_Merged_BHad_mass_4J, Merged_BHad_comb_mass_4J);
                    }

					// merged wjets
					if( (best_WJa == best_WJb) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){ // only wjets merged
					// correct perm
						Correct_merged_wjet_mass_4J = best_BHad->M();
						Correct_merged_wjet_comb_mass_4J = (*best_BHad+*best_WJa).M();
						Max_Correct_merged_wjet_mass_4J = ( best_BHad->M() > best_WJa->M() ) ? best_BHad->M() : best_WJa->M();
						t12->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_correct"].fill(Max_Correct_merged_wjet_mass_4J, Correct_merged_wjet_comb_mass_4J);

//                       // preselection
//                       delt_dir->second["4J_Preselection_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
//                       delt_dir->second["4J_Preselection_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
//
//                       // square cut
//                       if( best_BHad->M() > b_cut_ && (*best_BHad+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BHad+*best_WJa).M() > b_and_jet_lower_cut_ ){
//                               delt_dir->second["4J_Square_Cut_Correct_Perms_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
//                               delt_dir->second["4J_Square_Cut_Correct_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
//                       }

					// wrong perm
						Wrong_merged_wjet_mass_4J = best_BLep->M();
						Wrong_merged_wjet_comb_mass_4J = (*best_BLep+*best_WJa).M();
						Max_Wrong_merged_wjet_mass_4J = ( best_BLep->M() > best_WJa->M() ) ? best_BLep->M() : best_WJa->M();
						t8->Fill();
                        disc_dir->second["4J_mbpjet_vs_maxmjet_wrong"].fill(Max_Wrong_merged_wjet_mass_4J, Wrong_merged_wjet_comb_mass_4J);

//						// preselection
//						delt_dir->second["4J_Preselection_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
//						delt_dir->second["4J_Preselection_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
//	
//						// square cut
//						if( best_BLep->M() > b_cut_ && (*best_BLep+*best_WJa).M() < b_and_jet_upper_cut_ && (*best_BLep+*best_WJa).M() > b_and_jet_lower_cut_ ){
//							delt_dir->second["4J_Square_Cut_Wrong_Perms_Delta_mass"].fill(thad->M()-(*best_BLep+*best_WJa).M());
//							delt_dir->second["4J_Square_Cut_Wrong_Perms_Delta_costh"].fill(thad->CosTheta()-(*best_BLep+*best_WJa).CosTheta());
//						}

					} // end of merged W jets loop

				} // end of njets = 4 loop for disc and presel/signal

		//////////////// Hadronic side event categories
				if( njets_ == 3 ){
                                /// Hadronic side resolved
                                        if( best_BHad && best_WJa && best_WJb && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["3J_Had_Resolved_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only BHad and one WJet merged
                                        // BHad and WJa merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["3J_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                        // BHad and WJb merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJb) && !(best_BHad == best_WJa) && !(best_BHad == best_BLep) && !(best_WJa == best_BLep) ){
                                                reco_dir->second["3J_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only jets from W are merged
                                        if( best_BHad && best_WJa && best_WJb && (best_WJa == best_WJb) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["3J_Merged_WJets_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only three hadronic jets merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && (best_BHad == best_WJb) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["3J_Merged_THad_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// non-reconstructable events (one or more hadronic jets missing)
                                        if( !best_BHad || !best_WJa || !best_WJb ){
                                                reco_dir->second["3J_Non_Reconstructable_vs_mttbar"].fill(ttbar.M());

                                                if( best_BHad ){
                                                        reco_dir->second["3J_Non_Reconstructable_BHad_mass"].fill(best_BHad->M());
                                                        reco_dir->second["3J_Non_Reconstructable_BHad_pt"].fill(best_BHad->Pt());
                                                        reco_dir->second["3J_Non_Reconstructable_BHad_eta"].fill(best_BHad->Eta());
                                                }
                                                if( best_WJa ){
                                                        reco_dir->second["3J_Non_Reconstructable_WJa_mass"].fill(best_WJa->M());
                                                        reco_dir->second["3J_Non_Reconstructable_WJa_pt"].fill(best_WJa->Pt());
                                                        reco_dir->second["3J_Non_Reconstructable_WJa_eta"].fill(best_WJa->Eta());
                                                }
                                                if( best_WJb ){
                                                        reco_dir->second["3J_Non_Reconstructable_WJb_mass"].fill(best_WJb->M());
                                                        reco_dir->second["3J_Non_Reconstructable_WJb_pt"].fill(best_WJb->Pt());
                                                        reco_dir->second["3J_Non_Reconstructable_WJb_eta"].fill(best_WJb->Eta());
                                                }
                                        }
				} // end of njets = 3 loop for hadronic events

                                if( njets_ == 4 ){
                                /// Hadronic side resolved
                                        if( best_BHad && best_WJa && best_WJb && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["4J_Had_Resolved_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only BHad and one WJet merged
                                        // BHad and WJa merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["4J_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                        // BHad and WJb merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJb) && !(best_BHad == best_WJa) && !(best_BHad == best_BLep) && !(best_WJa == best_BLep) ){
                                                reco_dir->second["4J_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only jets from W are merged
                                        if( best_BHad && best_WJa && best_WJb && (best_WJa == best_WJb) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["4J_Merged_WJets_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only three hadronic jets merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && (best_BHad == best_WJb) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["4J_Merged_THad_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// non-reconstructable events (one or more hadronic jets missing)
                                        if( !best_BHad || !best_WJa || !best_WJb ){
                                                reco_dir->second["4J_Non_Reconstructable_vs_mttbar"].fill(ttbar.M());
                                        }

                                } // end of njets = 4 loop for hadronic events

                                if( njets_ > 4 ){

                                /// Hadronic side resolved
                                        if( best_BHad && best_WJa && best_WJb && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["5PJ_Had_Resolved_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only BHad and one WJet merged
                                        // BHad and WJa merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) && !(best_WJb == best_BLep) ){
                                                reco_dir->second["5PJ_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                        // BHad and WJb merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJb) && !(best_BHad == best_WJa) && !(best_BHad == best_BLep) && !(best_WJa == best_BLep) ){
                                                reco_dir->second["5PJ_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only jets from W are merged
                                        if( best_BHad && best_WJa && best_WJb && (best_WJa == best_WJb) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["5PJ_Merged_WJets_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// only three hadronic jets merged
                                        if( best_BHad && best_WJa && best_WJb && (best_BHad == best_WJa) && (best_BHad == best_WJb) && !(best_BHad == best_BLep) ){
                                                reco_dir->second["5PJ_Merged_THad_vs_mttbar"].fill(ttbar.M());
                                        }

                                /// non-reconstructable events (one or more hadronic jets missing)
                                        if( !best_BHad || !best_WJa || !best_WJb ){
                                                reco_dir->second["5PJ_Non_Reconstructable_vs_mttbar"].fill(ttbar.M());
                                        }

                                } // end of njets > 4 loop for hadronic events

			} // end of perm matching loop for different DR values

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
