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

class jet_effs : public AnalyzerBase
{
	private:
        //counters
	    unsigned long evt_idx_ = 0; //event index
	    double njets_ = 0;
	    double jets3_ = 0;
	    double matched_BHadWJa_ = 0;
	    double matched_BHadWJb_ = 0;
	    double matched_WJaWJb_ = 0;
	    double jet_pt_min_ = 25;
	    double jet_eta_cut_ = 2.4;

            const char *DRnames_[4] = {"DRP4", "DRP5", "DRP6", "DRP8"};
	    double DR_[4] = {0.4, 0.5, 0.6, 0.8};
            //const char *DRnames_[1] = {"DRP4"};
	    //double DR_[1] = {0.4};

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

	    double DR_min_ = 0.;
	    double DR_max_ = 8.;
	    double pt_min_ = 0.;
	    double pt_max_ = 2000.;
	    double eta_min_ = -8.;
	    double eta_max_ = 8.;


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

	public:
	    jet_effs(const std::string output_filename):
		AnalyzerBase("jet_effs", output_filename),
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

		opts::variables_map &values = URParser::instance().values();
		string output_file = values["output"].as<std::string>();
                string sample = systematics::get_sample(output_file);
		Logger::log().debug() << "		" << sample << endl;

		if( sample == "ttJetsM0" ){
			nbins = 22;
//                	m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
			mass_min_ = 250.;
			mass_max_ = 2000.;
		}
		if( sample == "ttJetsM700" ){
			nbins = 10;
//                	m_bins[nbins+1] = {700, 800, 900, 1000};
			mass_min_ = 700.;
			mass_max_ = 1000.;
		}
		if( sample == "ttJetsM1000" ){
			nbins = 10;
//                	m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
			mass_min_ = 1000.;
			mass_max_ = 2000.;
		}

		string folder = "Gen_Plots";

	//Gen Objects
			//DR
		book<TH2D>(folder, "DR_LepBHad_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DR_LepBLep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DR_LepWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DR_LepWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DR_BHadBLep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                book<TH2D>(folder, "DR_BHadWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                book<TH2D>(folder, "DR_BHadWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                book<TH2D>(folder, "DR_BLepWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DR_BLepWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                book<TH2D>(folder, "DR_WJaWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);

		book<TH1F>(folder, "DR_LepBHad", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_LepBLep", "", DR_bins_, DR_min_, DR_max_) ;
		book<TH1F>(folder, "DR_LepWJa", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_LepWJb", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_BHadBLep", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_BHadWJa", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_BHadWJb", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_BLepWJa", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_BLepWJb", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DR_WJaWJb", "", DR_bins_, DR_min_, DR_max_);

		book<TH2D>(folder, "DRmin_thad_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DRmin_tlep_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DRmin_thad_vs_ptthad", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DRmin_tlep_vs_pttlep", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
		
		book<TH1F>(folder, "DRmin_thad", "", DR_bins_, DR_min_, DR_max_);
		book<TH1F>(folder, "DRmin_tlep", "", DR_bins_, DR_min_, DR_max_);

		book<TH2D>(folder, "DRmax_thad_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
		book<TH2D>(folder, "DRmax_thad_vs_ptthad", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);

                book<TH1F>(folder, "DRmax_thad", "", DR_bins_, DR_min_, DR_max_);

		for( int i = 0; i < 4; i++ ){
			book<TH1F>(folder, std::string("DRmin_thad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmin_thad_l")+DRnames_[i]+"_vs_ptthad", "", nbins, pt_min_, pt_max_);
	                book<TH1F>(folder, std::string("DRmin_tlep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmin_tlep_l")+DRnames_[i]+"_vs_pttlep", "", nbins, pt_min_, pt_max_);
	
	                book<TH1F>(folder, std::string("DRmin_thad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmin_thad_g")+DRnames_[i]+"_vs_ptthad", "", nbins, pt_min_, pt_max_);
	                book<TH1F>(folder, std::string("DRmin_tlep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmin_tlep_g")+DRnames_[i]+"_vs_pttlep", "", nbins, pt_min_, pt_max_);

	                book<TH1F>(folder, std::string("DRmax_thad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmax_thad_l")+DRnames_[i]+"_vs_ptthad", "", nbins, pt_min_, pt_max_);
	
	                book<TH1F>(folder, std::string("DRmax_thad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DRmax_thad_g")+DRnames_[i]+"_vs_ptthad", "", nbins, pt_min_, pt_max_);
	
	                book<TH1F>(folder, std::string("DR_LepBHad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepBLep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadBLep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BLepWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BLepWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_WJaWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	
	                book<TH1F>(folder, std::string("DR_LepBHad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepBLep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_LepWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadBLep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BHadWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BLepWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_BLepWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
	                book<TH1F>(folder, std::string("DR_WJaWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
		}

			//Pt
                book<TH2D>(folder, "Pt_Lep_vs_Mtt", "m_{t#bar t} [GeV]; l p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH2D>(folder, "Pt_BLep_vs_Mtt", "m_{t#bar t} [GeV]; b_{l} p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH2D>(folder, "Pt_BHad_vs_Mtt", "m_{t#bar t} [GeV]; b_{h} p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH2D>(folder, "Pt_WJa_vs_Mtt", "m_{t#bar t} [GeV]; WJa p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH2D>(folder, "Pt_WJb_vs_Mtt", "m_{t#bar t} [GeV]; WJb p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

                book<TH1F>(folder, "Pt_Lep", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_BLep", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_BHad", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_WJa", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_WJb", "", pt_bins_, pt_min_, pt_max_);

                book<TH1F>(folder, "Pt_ttbar", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_thad", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(folder, "Pt_tlep", "", pt_bins_, pt_min_, pt_max_);

			//Eta
                book<TH2D>(folder, "Eta_Lep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH2D>(folder, "Eta_BLep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH2D>(folder, "Eta_BHad_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH2D>(folder, "Eta_WJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH2D>(folder, "Eta_WJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);

                book<TH1F>(folder, "Eta_Lep", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(folder, "Eta_BLep", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(folder, "Eta_BHad", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(folder, "Eta_WJa", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(folder, "Eta_WJb", "", eta_bins_, eta_min_, eta_max_);

			//System Plots
                book<TH1F>(folder, "Mass_ttbar", "", mass_bins_, mass_min_, mass_max_);
                book<TH1F>(folder, "Mass_thad", "", mass_bins_, 125., 250.);
                book<TH1F>(folder, "Mass_tlep", "", mass_bins_, 125., 250.);
                book<TH1F>(folder, "nJets", "", 13, 2.5, 15.5);
                book<TH1F>(folder, "nMatched_objects_3J", "", 4, 0.5, 4.5);
                book<TH1F>(folder, "nMatched_objects_4J", "", 4, 0.5, 4.5);
                book<TH1F>(folder, "nMatched_objects_5PJ", "", 4, 0.5, 4.5);

	//Reco Jets
		for( int i = 0; i < 4; i++ ){
				//Pt
        	        book<TH1F>(DRnames_[i], "Matched_perm_BHad_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_BLep_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_pt", "", pt_bins_, pt_min_, pt_max_);

      				//Eta
			book<TH1F>(DRnames_[i], "Matched_perm_BHad_eta", "", eta_bins_, eta_min_, eta_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_BLep_eta", "", eta_bins_, eta_min_, eta_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_eta", "", eta_bins_, eta_min_, eta_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_eta", "", eta_bins_, eta_min_, eta_max_);

        	               //break ttJetsM0 into 700-1000 and > 1000
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM700_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM700_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM1000_pt", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM1000_pt", "", pt_bins_, pt_min_, pt_max_);

        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM700_frac_p", "", 10, 0., 5.);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM700_frac_p", "", 10, 0., 5.);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM1000_frac_p", "", 10, 0., 5.);
        	        book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM1000_frac_p", "", 10, 0., 5.);

				//Matched objects
        	        book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);

        	        book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);

        	        book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
        	        book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);

				//BTagging
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BHad_tight_fail", "", nbins, mass_min_, mass_max_);

        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_less_BLep_tight_fail", "", nbins, mass_min_, mass_max_);

        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BHad_tight_fail", "", nbins, mass_min_, mass_max_);

        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
        	        book<TH1F>(DRnames_[i], "BTag_great_BLep_tight_fail", "", nbins, mass_min_, mass_max_);

				//Merged jets
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_and_WJb_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_WJb_mass", "", 100, 0., 200.);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJa_massDivpt", "", 50, 0., 0.5);

        	        book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_and_WJa_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_WJa_mass", "", 100, 0., 200.);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Merged_BHadWJb_massDivpt", "", 50, 0., 0.5);

        	        book<TH2D>(DRnames_[i], "Merged_BHadWJet_mass_vs_BHad_mass","",100, 0.,300.,100,0., 300.);

        	        book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_and_WJb_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJa_massDivpt", "", 50, 0., 0.5);

        	        book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_and_WJa_mass", "", 100, 0., 300.);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Merged_BLepWJb_massDivpt", "", 50, 0., 0.5);

        	        book<TH2D>(DRnames_[i], "Merged_BLepWJet_mass_vs_BLep_mass","",100, 0.,300.,100,0., 300.);

				//Unmerged jets
        	        book<TH1F>(DRnames_[i], "Unmerged_BHad_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJa_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJb_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Unmerged_BHad_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJa_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJb_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Unmerged_BHad_massDivpt", "", 50, 0., 0.5);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJa_massDivpt", "", 50, 0., 0.5);
        	        book<TH1F>(DRnames_[i], "Unmerged_WJb_massDivpt", "", 50, 0., 0.5);

        	        book<TH2D>(DRnames_[i], "Unmerged_BHadWJet_mass_vs_BHad_mass","",100, 0., 300., 100, 0., 300.);
        	        book<TH2D>(DRnames_[i], "Unmerged_BHadWJet_highest_mass_vs_BHad_mass", "", 100, 0., 300., 100, 0., 300.);

        	        book<TH1F>(DRnames_[i], "Unmerged_BLep_mass", "", 100, 0., 100.);
        	        book<TH1F>(DRnames_[i], "Unmerged_BLep_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
        	        book<TH1F>(DRnames_[i], "Unmerged_BLep_massDivpt", "", 50, 0., 0.5);

        	        book<TH2D>(DRnames_[i], "Unmerged_BLepWJet_mass_vs_BLep_mass","",100, 0., 300., 100, 0., 300.);
        	        book<TH2D>(DRnames_[i], "Unmerged_BLepWJet_highest_mass_vs_BLep_mass", "", 100, 0., 300., 100, 0., 300.);
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

			auto gen_dir = histos_.find("Gen_Plots");

			//generator selection
			bool selection = genp_selector_.select(event);
                        tracker_.track("gen selection");
                        if( !selection ){
                                Logger::log().debug() << "event has no selection " << endl;
                                continue;
                        }
			GenTTBar &ttbar = genp_selector_.ttbar_system();
			gen_dir->second["Pt_ttbar"].fill(ttbar.Pt());
                        gen_dir->second["Mass_ttbar"].fill(ttbar.M());

			//Gen level matching
			Permutation matched_perm;
	    		//double nmatched_objects = 0;

			if( object_selector_.select(event) ){
				gen_dir->second["nJets"].fill(object_selector_.clean_jets().size());
				if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
					matched_perm = matcher_.match(
						genp_selector_.ttbar_final_system(),
						object_selector_.clean_jets(),
						object_selector_.loose_electrons(),
						object_selector_.loose_muons());
				}
				matched_perm.SetMET(object_selector_.met());
//				if( object_selector_.clean_jets().size() == 3 ){
//					for(vector<IDJet*>::iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){ 
////						cout << "jet E: " << (*jets)->E() << endl;
////						cout << "jet |p|: " << (*jets)->P() << endl;
////						cout << "jet M: " << (*jets)->M() << endl;
//					}
//				}
			}
//			if( matched_perm.L() || matched_perm.BHad() || matched_perm.BLep() || matched_perm.WJa() || matched_perm.WJb() ) cout << matched_perm << endl;
//			cout << matched_perm << endl;

			//initialize kinematic vars
			double DR_LepBHad = -1;
			double DR_LepBLep = -1;
			double DR_LepWJa = -1;
			double DR_LepWJb = -1;
			double DR_BHadBLep = -1;
			double DR_BHadWJa = -1;
			double DR_BHadWJb = -1;
			double DR_BLepWJa = -1;
			double DR_BLepWJb = -1;
			double DR_WJaWJb = -1;
			
			double Pt_Lep = -1;
			double Pt_BHad = -1;
			double Pt_BLep = -1;
			double Pt_WJa = -1;
			double Pt_WJb = -1;
			
			double Eta_Lep = 1e3;
			double Eta_BHad = 1e3;
			double Eta_BLep = 1e3;
			double Eta_WJa = 1e3;
			double Eta_WJb = 1e3;

			//initialize gen jets
      			const GenObject* lepton = 0;
			const GenObject* BLep = 0;
      			const GenObject* BHad = 0;
			const GenObject* WJa = 0;
			const GenObject* WJb = 0;
			const GenTop* tlep = 0;
			const GenTop* thad = 0;

			//Generator level definitions
			if( ttbar.lep_top() ){ //check if lep top exists
				tlep = ttbar.lep_top();
				gen_dir->second["Pt_tlep"].fill(tlep->Pt());
				gen_dir->second["Mass_tlep"].fill(tlep->M());
			}
			if( ttbar.had_top() ){ //check if had top exists
				thad = ttbar.had_top();
				gen_dir->second["Pt_thad"].fill(thad->Pt());
				gen_dir->second["Mass_thad"].fill(thad->M());
			}
			if( ttbar.lepton() ){
				lepton = ttbar.lepton();
				Pt_Lep = lepton->Pt();
                        	gen_dir->second["Pt_Lep_vs_Mtt"].fill(ttbar.M(), Pt_Lep);
                        	gen_dir->second["Pt_Lep"].fill(Pt_Lep);

                        	Eta_Lep = lepton->Eta();
				gen_dir->second["Eta_Lep_vs_Mtt"].fill(ttbar.M(), Eta_Lep);
                        	gen_dir->second["Eta_Lep"].fill(Eta_Lep);
			}
			if( ttbar.lep_b() ){
				BLep = ttbar.lep_b();
				Pt_BLep = BLep->Pt();
                       		gen_dir->second["Pt_BLep_vs_Mtt"].fill(ttbar.M(), Pt_BLep);
                        	gen_dir->second["Pt_BLep"].fill(Pt_BLep);

                        	Eta_BLep = BLep->Eta();
                        	gen_dir->second["Eta_BLep_vs_Mtt"].fill(ttbar.M(), Eta_BLep);
                        	gen_dir->second["Eta_BLep"].fill(Eta_BLep);
			}
			if( ttbar.had_b() ){
				BHad = ttbar.had_b();
				Pt_BHad = BHad->Pt();
                        	gen_dir->second["Pt_BHad_vs_Mtt"].fill(ttbar.M(), Pt_BHad);
                        	gen_dir->second["Pt_BHad"].fill(Pt_BHad);

                        	Eta_BHad = BHad->Eta();
                        	gen_dir->second["Eta_BHad_vs_Mtt"].fill(ttbar.M(), Eta_BHad);
                        	gen_dir->second["Eta_BHad"].fill(Eta_BHad);
//				cout << "BHad Energy, |p|: " << BHad->E() << ", " << sqrt(pow(BHad->Px(), 2.0)+pow(BHad->Py(),2.0)+pow(BHad->Pz(),2.0)) << endl;
			}
			if( ttbar.had_W() ){
				if( ttbar.had_W()->first && !ttbar.had_W()->second ){
					WJa = ttbar.had_W()->first;
					Pt_WJa = WJa->Pt();
                        		gen_dir->second["Pt_WJa_vs_Mtt"].fill(ttbar.M(), Pt_WJa);
                        		gen_dir->second["Pt_WJa"].fill(Pt_WJa);

                        		Eta_WJa = WJa->Eta();
                        		gen_dir->second["Eta_WJa_vs_Mtt"].fill(ttbar.M(), Eta_WJa);
                        		gen_dir->second["Eta_WJa"].fill(Eta_WJa);
				}
				if( ttbar.had_W()->second && !ttbar.had_W()->first ){
					WJb = ttbar.had_W()->second;
					Pt_WJb = WJb->Pt();
                        		gen_dir->second["Pt_WJb_vs_Mtt"].fill(ttbar.M(), Pt_WJb);
                        		gen_dir->second["Pt_WJb"].fill(Pt_WJb);

                        		Eta_WJb = WJb->Eta();
                        		gen_dir->second["Eta_WJb_vs_Mtt"].fill(ttbar.M(), Eta_WJb);
                        		gen_dir->second["Eta_WJb"].fill(Eta_WJb);
				}
				if( ttbar.had_W()->first && ttbar.had_W()->second ){
					WJa = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
					WJb = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;

					Pt_WJa = WJa->Pt();
                        		gen_dir->second["Pt_WJa_vs_Mtt"].fill(ttbar.M(), Pt_WJa);
                        		gen_dir->second["Pt_WJa"].fill(Pt_WJa);

                        		Eta_WJa = WJa->Eta();
                        		gen_dir->second["Eta_WJa_vs_Mtt"].fill(ttbar.M(), Eta_WJa);
                        		gen_dir->second["Eta_WJa"].fill(Eta_WJa);

					Pt_WJb = WJb->Pt();
                        		gen_dir->second["Pt_WJb_vs_Mtt"].fill(ttbar.M(), Pt_WJb);
                        		gen_dir->second["Pt_WJb"].fill(Pt_WJb);

                        		Eta_WJb = WJb->Eta();
                        		gen_dir->second["Eta_WJb_vs_Mtt"].fill(ttbar.M(), Eta_WJb);
                        		gen_dir->second["Eta_WJb"].fill(Eta_WJb);
				}
			}

				//DR
			if( !(lepton == 0) && !(BHad == 0) ){
				DR_LepBHad = lepton->DeltaR(*BHad);
				gen_dir->second["DR_LepBHad_vs_Mtt"].fill(ttbar.M(), DR_LepBHad);
				gen_dir->second["DR_LepBHad"].fill(DR_LepBHad);

				for( int i = 0; i < 4; i++ ){
					if( DR_LepBHad < DR_[i] ){
					        gen_dir->second[std::string("DR_LepBHad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
					}
					if( DR_LepBHad >= DR_[i] ){
					        gen_dir->second[std::string("DR_LepBHad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
					}
				}
			}
                        if( !(lepton == 0) && !(BLep == 0) ){ //check if lepton and lep b exist
			        DR_LepBLep = lepton->DeltaR(*BLep);
			        gen_dir->second["DR_LepBLep_vs_Mtt"].fill(ttbar.M(), DR_LepBLep);
			        gen_dir->second["DR_LepBLep"].fill(DR_LepBLep);
			        gen_dir->second["DRmin_tlep_vs_mttbar"].fill(ttbar.M(), DR_LepBLep);
			        gen_dir->second["DRmin_tlep_vs_pttlep"].fill(tlep->Pt(), DR_LepBLep);
			        gen_dir->second["DRmin_tlep"].fill(DR_LepBLep);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_LepBLep < DR_[i] ){
			        	        gen_dir->second[std::string("DRmin_tlep_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	        gen_dir->second[std::string("DRmin_tlep_l")+DRnames_[i]+"_vs_pttlep"].fill(tlep->Pt());
			        	        gen_dir->second[std::string("DR_LepBLep_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
			        	}
			        	if( DR_LepBLep >= DR_[i] ){
			        	        gen_dir->second[std::string("DRmin_tlep_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
			        	        gen_dir->second[std::string("DRmin_tlep_g")+DRnames_[i]+"_vs_pttlep"].fill(tlep->Pt());
			        	        gen_dir->second[std::string("DR_LepBLep_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
			        	}
				}
			}
			if( !(lepton == 0) && !(WJa == 0) ){ //check if lepton and WJa exist
			        DR_LepWJa = lepton->DeltaR(*WJa);
			        gen_dir->second["DR_LepWJa_vs_Mtt"].fill(ttbar.M(), DR_LepWJa);
			        gen_dir->second["DR_LepWJa"].fill(DR_LepWJa);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_LepWJa < DR_[i] ){
			        	        gen_dir->second[std::string("DR_LepWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
			        	if( DR_LepWJa >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_LepWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
				}
			}
			  if( !(lepton == 0) && !(WJb == 0) ){
			        DR_LepWJb = lepton->DeltaR(*WJb);
			        gen_dir->second["DR_LepWJb_vs_Mtt"].fill(ttbar.M(), DR_LepWJb);
			        gen_dir->second["DR_LepWJb"].fill(DR_LepWJb);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_LepWJb < DR_[i] ){
			        	        gen_dir->second[std::string("DR_LepWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
			        	if( DR_LepWJb >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_LepWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
				}
			}
                        if( !(BHad == 0) && !(BLep == 0) ){
			        DR_BHadBLep = BHad->DeltaR(*BLep);
			        gen_dir->second["DR_BHadBLep_vs_Mtt"].fill(ttbar.M(), DR_BHadBLep);
			        gen_dir->second["DR_BHadBLep"].fill(DR_BHadBLep);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_BHadBLep < DR_[i] ){
			        	        gen_dir->second[std::string("DR_BHadBLep_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
			        	if( DR_BHadBLep >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_BHadBLep_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
				}
			}
                        if( !(BHad == 0) && !(WJa == 0) ){
				DR_BHadWJa = BHad->DeltaR(*WJa);
				gen_dir->second["DR_BHadWJa_vs_Mtt"].fill(ttbar.M(), DR_BHadWJa);
				gen_dir->second["DR_BHadWJa"].fill(DR_BHadWJa);

				for( int i = 0; i < 4; i++ ){
					if( DR_BHadWJa < DR_[i] ){
					        gen_dir->second[std::string("DR_BHadWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
                        		}
					if( DR_BHadWJa >= DR_[i] ){
					        gen_dir->second[std::string("DR_BHadWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
					}
				}
			}
                        if( !(BHad == 0) && !(WJb == 0) ){
				DR_BHadWJb = BHad->DeltaR(*WJb);
				gen_dir->second["DR_BHadWJb_vs_Mtt"].fill(ttbar.M(), DR_BHadWJb);
				gen_dir->second["DR_BHadWJb"].fill(DR_BHadWJb);

				for( int i = 0; i < 4; i++ ){
					if( DR_BHadWJb < DR_[i] ){
					        gen_dir->second[std::string("DR_BHadWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
                        		}
					if( DR_BHadWJb >= DR_[i] ){
					        gen_dir->second[std::string("DR_BHadWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
					}
				}
			}
                        if( !(BLep == 0) && !(WJa == 0) ){
			        DR_BLepWJa = BLep->DeltaR(*WJa);
			        gen_dir->second["DR_BLepWJa_vs_Mtt"].fill(ttbar.M(), DR_BLepWJa);
			        gen_dir->second["DR_BLepWJa"].fill(DR_BLepWJa);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_BLepWJa < DR_[i] ){
			        	        gen_dir->second[std::string("DR_BLepWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
			        	if( DR_BLepWJa >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_BLepWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
				}
			}
			if( !(BLep == 0) && !(WJb == 0) ){
			        DR_BLepWJb = BLep->DeltaR(*WJb);
			        gen_dir->second["DR_BLepWJb_vs_Mtt"].fill(ttbar.M(), DR_BLepWJb);
			        gen_dir->second["DR_BLepWJb"].fill(DR_BLepWJb);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_BLepWJb < DR_[i] ){
			        	        gen_dir->second[std::string("DR_BLepWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
			        	if( DR_BLepWJb >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_BLepWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
			        	}
				}
			}
                        if( !(WJa == 0) && !(WJb == 0) ){
			        DR_WJaWJb = WJa->DeltaR(*WJb);
			        gen_dir->second["DR_WJaWJb_vs_Mtt"].fill(ttbar.M(), DR_WJaWJb);
			        gen_dir->second["DR_WJaWJb"].fill(DR_WJaWJb);

				for( int i = 0; i < 4; i++ ){
			        	if( DR_WJaWJb < DR_[i] ){
			        	        gen_dir->second[std::string("DR_WJaWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());   
			        	}
			        	if( DR_WJaWJb >= DR_[i] ){
			        	        gen_dir->second[std::string("DR_WJaWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
			        	}
				}
			}

			//find min value of DR for thad objects
			list<double> thad_DR_list;
			list<double>::iterator thad_DR_it;
			if( !(BHad == 0) && !(WJa == 0) ) thad_DR_list.push_back(DR_BHadWJa);
			if( !(BHad == 0) && !(WJb == 0) ) thad_DR_list.push_back(DR_BHadWJb);
			if( !(WJa == 0) && !(WJb == 0) ) thad_DR_list.push_back(DR_WJaWJb);
			double thad_DRmin = 1e2;
			double thad_DRmax = 0;
			for( thad_DR_it = thad_DR_list.begin(); thad_DR_it != thad_DR_list.end(); ++thad_DR_it ){
				if( *thad_DR_it < thad_DRmin ){
					thad_DRmin = *thad_DR_it;
				}
				if( *thad_DR_it > thad_DRmax ){
					thad_DRmax = *thad_DR_it;
				}
			}
			if( thad_DRmin < 0.4 ){
				for( thad_DR_it = thad_DR_list.begin(); thad_DR_it != thad_DR_list.end(); ++thad_DR_it ){
					if( *thad_DR_it > thad_DRmax ){
						thad_DRmax = *thad_DR_it;
					}
				}

			}
			if( !(thad_DRmax == 0) ){
				gen_dir->second["DRmax_thad_vs_mttbar"].fill(ttbar.M(), thad_DRmax);
				gen_dir->second["DRmax_thad_vs_ptthad"].fill(thad->Pt(), thad_DRmax);
				gen_dir->second["DRmax_thad"].fill(thad_DRmax);

				for( int i = 0; i < 4; i++ ){
					if( thad_DRmax < DR_[i] ){
						gen_dir->second[std::string("DRmax_thad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
						gen_dir->second[std::string("DRmax_thad_l")+DRnames_[i]+"_vs_ptthad"].fill(thad->Pt());
					}
					if( thad_DRmax >= DR_[i] ){
						gen_dir->second[std::string("DRmax_thad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
						gen_dir->second[std::string("DRmax_thad_g")+DRnames_[i]+"_vs_ptthad"].fill(thad->Pt());
					}
				}
			}
			if( !(thad_DRmin == 1e2) ){
				gen_dir->second["DRmin_thad_vs_mttbar"].fill(ttbar.M(), thad_DRmin);
				gen_dir->second["DRmin_thad_vs_ptthad"].fill(thad->Pt(), thad_DRmin);
				gen_dir->second["DRmin_thad"].fill(thad_DRmin);

				for( int i = 0; i < 4; i++ ){
					if( thad_DRmin < DR_[i] ){
						gen_dir->second[std::string("DRmin_thad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
						gen_dir->second[std::string("DRmin_thad_l")+DRnames_[i]+"_vs_ptthad"].fill(thad->Pt());
					}
					if( thad_DRmin >= DR_[i] ){
						gen_dir->second[std::string("DRmin_thad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
						gen_dir->second[std::string("DRmin_thad_g")+DRnames_[i]+"_vs_ptthad"].fill(thad->Pt());
					}
				}
			}


			//Perm Matching
			for( int i = 0; i < 4; i++ ){

				auto reco_dir = histos_.find(DRnames_[i]);
				//initialize permutation objects
					//DR=0.4
				const IDJet* best_BHad = 0;
				const IDJet* best_BLep = 0;
				const IDJet* best_WJa = 0;
				const IDJet* best_WJb = 0;
//
//					//DR=0.5
//				const IDJet* best_BHad_DRP5 = 0;
//				const IDJet* best_BLep_DRP5 = 0;
//				const IDJet* best_WJa_DRP5 = 0;
//				const IDJet* best_WJb_DRP5 = 0;
//
//					//DR=0.6
//				const IDJet* best_BHad_DRP6 = 0;
//				const IDJet* best_BLep_DRP6 = 0;
//				const IDJet* best_WJa_DRP6 = 0;
//				const IDJet* best_WJb_DRP6 = 0;
//
//					//DR=0.8
//				const IDJet* best_BHad_DRP8 = 0;
//				const IDJet* best_BLep_DRP8 = 0;
//				const IDJet* best_WJa_DRP8 = 0;
//				const IDJet* best_WJb_DRP8 = 0;
//
					//DR=0.4
				list<IDJet*> bhad_list;
				float bhad_DR = 1e10;// itialize dr to high number
				if( object_selector_.clean_jets().size() < 3 ) continue;
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

				//Reco Obj. Kinematic Vars
				if( !(best_BHad == 0) ){//reco BHad exists
					reco_dir->second["Matched_perm_BHad_pt"].fill(best_BHad->Pt());
					reco_dir->second["Matched_perm_BHad_eta"].fill(best_BHad->Eta());
//							++nmatched_objects;
					if( best_BHad == best_BLep || best_BHad == best_WJa || best_BHad == best_WJb ){//reco BHad merged with something
						if( best_BHad->BTagId(cut_loose_b_) ) reco_dir->second["BTag_less_BHad_loose_pass"].fill(ttbar.M());
						if( best_BHad->BTagId(cut_medium_b_) ) reco_dir->second["BTag_less_BHad_medium_pass"].fill(ttbar.M());	
						if( best_BHad->BTagId(cut_tight_b_) ) reco_dir->second["BTag_less_BHad_tight_pass"].fill(ttbar.M());
						if( !best_BHad->BTagId(cut_loose_b_) ) reco_dir->second["BTag_less_BHad_loose_fail"].fill(ttbar.M());
						if( !best_BHad->BTagId(cut_medium_b_) ) reco_dir->second["BTag_less_BHad_medium_fail"].fill(ttbar.M());	
						if( !best_BHad->BTagId(cut_tight_b_) ) reco_dir->second["BTag_less_BHad_tight_fail"].fill(ttbar.M());
					}
				}
				if( !(best_BLep == 0) ){//reco BLep exists
					reco_dir->second["Matched_perm_BLep_pt"].fill(best_BLep->Pt());
					reco_dir->second["Matched_perm_BLep_eta"].fill(best_BLep->Eta());
//					++nmatched_objects; 
					if( best_BLep == best_BHad || best_BLep == best_WJa || best_BLep == best_WJb ){//reco BLep merged with something
						if( best_BLep->BTagId(cut_loose_b_) ) reco_dir->second["BTag_less_BLep_loose_pass"].fill(ttbar.M());
						if( best_BLep->BTagId(cut_medium_b_) ) reco_dir->second["BTag_less_BLep_medium_pass"].fill(ttbar.M());	
						if( best_BLep->BTagId(cut_tight_b_) ) reco_dir->second["BTag_less_BLep_tight_pass"].fill(ttbar.M());
						if( !best_BLep->BTagId(cut_loose_b_) ) reco_dir->second["BTag_less_BLep_loose_fail"].fill(ttbar.M());
						if( !best_BLep->BTagId(cut_medium_b_) ) reco_dir->second["BTag_less_BLep_medium_fail"].fill(ttbar.M());	
						if( !best_BLep->BTagId(cut_tight_b_) ) reco_dir->second["BTag_less_BLep_tight_fail"].fill(ttbar.M());
					}
				}
				if( !(best_WJa == 0) ){//reco WJa exists
//							++nmatched_objects; 
					reco_dir->second["Matched_perm_WJa_pt"].fill(best_WJa->Pt());
					reco_dir->second["Matched_perm_WJa_eta"].fill(best_WJa->Eta());
                        		if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                        		        reco_dir->second["Matched_perm_WJa_ttbarM700_pt"].fill(best_WJa->Pt());
                        		        reco_dir->second["Matched_perm_WJa_ttbarM700_frac_p"].fill(best_WJa->P()/(ttbar.M()/2));
                        		}
                        		if( ttbar.M() >= 1000 ){
                        		        reco_dir->second["Matched_perm_WJa_ttbarM1000_pt"].fill(best_WJa->Pt());
                        		        reco_dir->second["Matched_perm_WJa_ttbarM1000_frac_p"].fill(best_WJa->P()/(ttbar.M()/2));
                        		}
				}
				if( !(best_WJb == 0) ){//reco WJb exists
//							++nmatched_objects; 
					reco_dir->second["Matched_perm_WJb_pt"].fill(best_WJb->Pt());
					reco_dir->second["Matched_perm_WJb_eta"].fill(best_WJb->Eta());
					if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
					        reco_dir->second["Matched_perm_WJb_ttbarM700_pt"].fill(best_WJb->Pt());
					        reco_dir->second["Matched_perm_WJb_ttbarM700_frac_p"].fill(best_WJb->P()/(ttbar.M()/2));
					}
					if( ttbar.M() >= 1000 ){
					        reco_dir->second["Matched_perm_WJb_ttbarM1000_pt"].fill(best_WJb->Pt());
					        reco_dir->second["Matched_perm_WJb_ttbarM1000_frac_p"].fill(best_WJb->P()/(ttbar.M()/2));
					}
				}
				if( !(best_BLep == 0) && !(best_BLep == best_BHad) && !(best_BLep == best_WJa) && !(best_BLep == best_WJb) ){//reco BLep not merged with anything
					if( best_BLep->BTagId(cut_loose_b_) ) reco_dir->second["BTag_great_BLep_loose_pass"].fill(ttbar.M());
					if( best_BLep->BTagId(cut_medium_b_) ) reco_dir->second["BTag_great_BLep_medium_pass"].fill(ttbar.M());	
					if( best_BLep->BTagId(cut_tight_b_) ) reco_dir->second["BTag_great_BLep_tight_pass"].fill(ttbar.M());
					if( !best_BLep->BTagId(cut_loose_b_) ) reco_dir->second["BTag_great_BLep_loose_fail"].fill(ttbar.M());
					if( !best_BLep->BTagId(cut_medium_b_) ) reco_dir->second["BTag_great_BLep_medium_fail"].fill(ttbar.M());	
					if( !best_BLep->BTagId(cut_tight_b_) ) reco_dir->second["BTag_great_BLep_tight_fail"].fill(ttbar.M());

					reco_dir->second["Unmerged_BLep_mass"].fill(best_BLep->M());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Unmerged_BLep_DRLepBLep"].fill(DR_LepBLep);
	
					reco_dir->second["Unmerged_BLep_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
					if( !(best_WJa == 0) ) reco_dir->second["Unmerged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJa).M());
					if( !(best_WJb == 0) ) reco_dir->second["Unmerged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJb).M());
					if( !(best_WJa == 0) && !(best_WJb == 0) ){
						if( best_WJa->M() > best_WJb->M() ) reco_dir->second["Unmerged_BLepWJet_highest_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJa).M());
						if( best_WJa->M() < best_WJb->M() ) reco_dir->second["Unmerged_BLepWJet_highest_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJb).M());
					}


				}
				if( !(best_BHad == 0) && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) ){//reco BHad not merged with anything
					if( best_BHad->BTagId(cut_loose_b_) ) reco_dir->second["BTag_great_BHad_loose_pass"].fill(ttbar.M());
					if( best_BHad->BTagId(cut_medium_b_) ) reco_dir->second["BTag_great_BHad_medium_pass"].fill(ttbar.M());	
					if( best_BHad->BTagId(cut_tight_b_) ) reco_dir->second["BTag_great_BHad_tight_pass"].fill(ttbar.M());
					if( !best_BHad->BTagId(cut_loose_b_) ) reco_dir->second["BTag_great_BHad_loose_fail"].fill(ttbar.M());
					if( !best_BHad->BTagId(cut_medium_b_) ) reco_dir->second["BTag_great_BHad_medium_fail"].fill(ttbar.M());	
					if( !best_BHad->BTagId(cut_tight_b_) ) reco_dir->second["BTag_great_BHad_tight_fail"].fill(ttbar.M());

					reco_dir->second["Unmerged_BHad_mass"].fill(best_BHad->M());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Unmerged_BHad_DRLepBLep"].fill(DR_LepBLep);
	
					reco_dir->second["Unmerged_BHad_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
					if( !(best_WJa == 0) ) reco_dir->second["Unmerged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJa).M());
					if( !(best_WJb == 0) ) reco_dir->second["Unmerged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJb).M());
					if( !(best_WJa == 0) && !(best_WJb == 0) ){
						if( best_WJa->M() > best_WJb->M() ) reco_dir->second["Unmerged_BHadWJet_highest_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJa).M());
						if( best_WJa->M() < best_WJb->M() ) reco_dir->second["Unmerged_BHadWJet_highest_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJb).M());
					}
				}
				if( !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) ){//reco WJa not merged with anything
					reco_dir->second["Unmerged_WJa_massDivpt"].fill(best_WJa->M()/best_WJa->Pt());
					reco_dir->second["Unmerged_WJa_mass"].fill(best_WJa->M());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Unmerged_WJa_DRLepBLep"].fill(DR_LepBLep);
				}
				if( !(best_WJb == 0) && !(best_WJb == best_BHad) && !(best_WJb == best_WJa) && !(best_WJb == best_BLep) ){//reco WJb not merged with anything
					reco_dir->second["Unmerged_WJb_massDivpt"].fill(best_WJb->M()/best_WJb->Pt());
					reco_dir->second["Unmerged_WJb_mass"].fill(best_WJb->M());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Unmerged_WJb_DRLepBLep"].fill(DR_LepBLep);
				}
				if( !(best_BHad == 0) && !(best_WJa == 0) ){ //reco BHad and WJa exist
					if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["All_Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["All_Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["All_Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
				}
				if( !(best_BHad == 0) && !(best_WJb == 0) ){ //reco BHad and WJb exist
					if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["All_Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["All_Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["All_Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
				}
				if( !(best_WJa == 0) && !(best_WJb == 0) ){ //reco WJa and WJb exist
					if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["All_Matched_WJaWJb_ptthad_3J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["All_Matched_WJaWJb_ptthad_4J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["All_Matched_WJaWJb_ptthad_5PJ"].fill(thad->Pt());
				}
				if( best_BHad == best_WJa && !(best_BHad == 0) && !(best_BHad == best_BLep) && !(best_BHad == best_WJb) ){//only reco BHad and WJa merged
					reco_dir->second["Merged_BHadWJa_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Merged_BHadWJa_perm_DRLepBLep"].fill(DR_LepBLep);
					reco_dir->second["Merged_BHadWJa_perm_mass"].fill(best_BHad->M());//merged b and wjet mass
					if( !(best_WJb == 0) && !(best_WJb == best_BHad) && !(best_WJb == best_BLep) ){//reco WJb exists and not merged
						reco_dir->second["Merged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJb).M());
						reco_dir->second["Merged_BHadWJa_perm_and_WJb_mass"].fill((*best_BHad+*best_WJb).M());//comb invariant mass of merged jet and other wjet
						reco_dir->second["Merged_WJb_mass"].fill(best_WJb->M());
					}
					if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
				}
				if( best_BHad == best_WJb && !(best_BHad == 0) && !(best_BHad == best_BLep) && !(best_BHad == best_WJa) ){//only reco BHad and WJb merged
					reco_dir->second["Merged_BHadWJb_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Merged_BHadWJb_perm_DRLepBLep"].fill(DR_LepBLep);
					reco_dir->second["Merged_BHadWJb_perm_mass"].fill(best_BHad->M());//merged b and wjet mass
					if( !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){//reco WJa exists and not merged
						reco_dir->second["Merged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJa).M());
						reco_dir->second["Merged_BHadWJb_perm_and_WJa_mass"].fill((*best_BHad+*best_WJa).M());//comb mass of merged jet and other wjet
						reco_dir->second["Merged_WJa_mass"].fill(best_WJa->M());
					}
					if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
					if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
				}
				if( best_WJa == best_WJb && !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){//only reco WJa and WJb merged
							if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["Matched_WJaWJb_ptthad_3J"].fill(thad->Pt());
							if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["Matched_WJaWJb_ptthad_4J"].fill(thad->Pt());
							if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["Matched_WJaWJb_ptthad_5PJ"].fill(thad->Pt());
				}

				if( best_BLep == best_WJa && !(best_BLep == 0) && !(best_BHad == best_BLep) && !(best_BLep == best_WJb) ){//only reco BLep and WJa merged
					reco_dir->second["Merged_BLepWJa_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Merged_BLepWJa_perm_DRLepBLep"].fill(DR_LepBLep);
					reco_dir->second["Merged_BLepWJa_perm_mass"].fill(best_BLep->M());//merged b and wjet mass
					if( !(best_WJb == 0) && !(best_WJb == best_BLep) && !(best_WJb == best_BHad) ){//reco WJb exists and not merged
						reco_dir->second["Merged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(), (*best_BLep+*best_WJb).M());
						reco_dir->second["Merged_BLepWJa_perm_and_WJb_mass"].fill((*best_BLep+*best_WJb).M());//comb invariant mass of merged jet and other wjet
						//reco_dir->second["Merged_WJb_mass"].fill(best_WJb_DRP4->M());
					}
					//if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
					//if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
					//if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
				}
				if( best_BLep == best_WJb && !(best_BLep == 0) && !(best_BHad == best_BLep) && !(best_BLep == best_WJa) ){//only reco BLep and WJb merged
					reco_dir->second["Merged_BLepWJb_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
					if( !(DR_LepBLep == -1) ) reco_dir->second["Merged_BLepWJb_perm_DRLepBLep"].fill(DR_LepBLep);
					reco_dir->second["Merged_BLepWJb_perm_mass"].fill(best_BLep->M());//merged b and wjet mass
					if( !(best_WJa == 0) && !(best_WJa == best_BLep) && !(best_WJa == best_BHad) ){//reco WJa exists and not merged
						reco_dir->second["Merged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(), (*best_BLep+*best_WJa).M());
						reco_dir->second["Merged_BLepWJb_perm_and_WJa_mass"].fill((*best_BLep+*best_WJa).M());//comb mass of merged jet and other wjet
					//	reco_dir->second["Merged_WJa_mass"].fill(best_WJa_DRP4->M());
					}
					//if( object_selector_.clean_jets().size() == 3 ) reco_dir->second["Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
					//if( object_selector_.clean_jets().size() == 4 ) reco_dir->second["Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
					//if( object_selector_.clean_jets().size() > 4 ) reco_dir->second["Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
				}
			}

//			if( !(nmatched_objects == 0) ){
//				if( object_selector_.clean_jets().size() == 3 ) nMatched_objects_3J_hist->Fill(nmatched_objects);
//				if( object_selector_.clean_jets().size() == 4 ) nMatched_objects_4J_hist->Fill(nmatched_objects);
//				if( object_selector_.clean_jets().size() > 4 ) nMatched_objects_5PJ_hist->Fill(nmatched_objects);
//			}

		}
		Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
	}

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
	URDriver<jet_effs> test;
	int thing = test.run();
	return thing;
}
