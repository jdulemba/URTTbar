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
//#include "Analyses/URTTbar/interface/TTBarPlots.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"

#include "TROOT.h"
#include "TTree.h"

using namespace TMath;

class Test_Analyzer : public AnalyzerBase
{
    private:
	//Counters
	unsigned long evt_idx_ = 0; //event index
	unsigned long gen_match_ = 0; //number of gen matched events

	//  double unmatched_muons_[22] =  {0.0};
	//  double matched_muons_[22] = {0.0};
	//  double Mu_Effs_[22] = {0.0};
	//  int bin_index_ = 0;
 
	//histograms 
	CutFlowTracker tracker_; //tracker_ tracks how many events pass a certain cut

	//gen var plots
		//DR
	TH2D *DR_LepBHad_vs_Mtt_hist = 0;
	TH2D *DR_LepBLep_vs_Mtt_hist = 0;
	TH2D *DR_LepWJa_vs_Mtt_hist = 0;
	TH2D *DR_LepWJb_vs_Mtt_hist = 0;
	TH2D *DR_BHadBLep_vs_Mtt_hist = 0;
	TH2D *DR_BHadWJa_vs_Mtt_hist = 0;
	TH2D *DR_BHadWJb_vs_Mtt_hist = 0;
	TH2D *DR_BLepWJa_vs_Mtt_hist = 0;
	TH2D *DR_BLepWJb_vs_Mtt_hist = 0;
	TH2D *DR_WJaWJb_vs_Mtt_hist = 0;

	TH1F *DR_LepBHad_hist = 0;
	TH1F *DR_LepBLep_hist = 0;
	TH1F *DR_LepWJa_hist = 0;
	TH1F *DR_LepWJb_hist = 0;
	TH1F *DR_BHadBLep_hist = 0;
	TH1F *DR_BHadWJa_hist = 0;
	TH1F *DR_BHadWJb_hist = 0;
	TH1F *DR_BLepWJa_hist = 0;
	TH1F *DR_BLepWJb_hist = 0;
	TH1F *DR_WJaWJb_hist = 0;
		//Pt
	TH2D *Pt_Lep_vs_Mtt_hist = 0;
	TH2D *Pt_BLep_vs_Mtt_hist = 0;
	TH2D *Pt_BHad_vs_Mtt_hist = 0;
	TH2D *Pt_WJa_vs_Mtt_hist = 0;
	TH2D *Pt_WJb_vs_Mtt_hist = 0;

	TH1F *Pt_Lep_hist = 0;
	TH1F *Pt_BLep_hist = 0;
	TH1F *Pt_BHad_hist = 0;
	TH1F *Pt_WJa_hist = 0;
	TH1F *Pt_WJb_hist = 0;
		//Eta
	TH2D *Eta_Lep_vs_Mtt_hist = 0;
	TH2D *Eta_BLep_vs_Mtt_hist = 0;
	TH2D *Eta_BHad_vs_Mtt_hist = 0;
	TH2D *Eta_WJa_vs_Mtt_hist = 0;
	TH2D *Eta_WJb_vs_Mtt_hist = 0;

	TH1F *Eta_Lep_hist = 0;
	TH1F *Eta_BLep_hist = 0;
	TH1F *Eta_BHad_hist = 0;
	TH1F *Eta_WJa_hist = 0;
	TH1F *Eta_WJb_hist = 0;

//	//Muon Hists
//	TH1F *MuISO_hist = 0;
//	TH2D *LepBLep_DR_vs_mtt_hist = 0;
//	TH2D *MuISO_vs_mtt_hist = 0;
//	TH2D *MuDR_vs_MuRelDPT_hist = 0;
//
//	TH1F *delta_pt_lep_hist = 0;
//	TH1F *delta_phi_lep_hist = 0;
//	TH1F *delta_y_lep_hist = 0;
//
//		//muon efficiency hists
//	TH1D *Unmatched_Muons_hist = 0;
//	TH1D *Matched_Muons_hist = 0;
//	TH1D *Non_ISO_Mu_hist = 0;
//	TH1D *ISO_Mu_hist = 0;
//
//	//BJet Hists
//	TH1F *delta_pt_BHad_hist = 0;
//	TH1F *delta_pt_BLep_hist = 0;
//	TH1F *delta_phi_BHad_hist = 0;
//	TH1F *delta_phi_BLep_hist = 0;
//	TH1F *delta_y_BHad_hist = 0;
//	TH1F *delta_y_BLep_hist = 0;
//	
//	TH1F *Unmatched_BLep_pt_hist = 0;
//	TH1F *Unmatched_BHad_pt_hist = 0;
//	TH1F *Unmatched_BLep_y_hist = 0;
//	TH1F *Unmatched_BHad_y_hist = 0;
//		//BJet efficiency hists
//	TH1D *Unmatched_BHad_hist = 0;
//	TH1D *Unmatched_BLep_hist = 0;
//	TH1D *Matched_BHad_hist = 0;
//	TH1D *Matched_BLep_hist = 0;
//
//	TH1D *Unmatched_BHad_zoom_hist = 0;
//	TH1D *Unmatched_BLep_zoom_hist = 0;
//	TH1D *Matched_BHad_zoom_hist = 0;
//	TH1D *Matched_BLep_zoom_hist = 0;
//
//	//WJet Hists
//	TH1F *delta_pt_WJa_hist = 0;
//	TH1F *delta_pt_WJb_hist = 0;
//	TH1F *delta_phi_WJa_hist = 0;
//	TH1F *delta_phi_WJb_hist = 0;
//	TH1F *delta_y_WJa_hist = 0;
//	TH1F *delta_y_WJb_hist = 0;
//
//	TH1F *Unmatched_WJa_pt_hist = 0;
//	TH1F *Unmatched_WJb_pt_hist = 0;
//	TH1F *Unmatched_WJa_y_hist = 0;
//	TH1F *Unmatched_WJb_y_hist = 0;
//		//WJet efficiency hist
//	TH1D *Unmatched_WJa_hist = 0;
//	TH1D *Unmatched_WJb_hist = 0;
//	TH1D *Matched_WJa_hist = 0;
//	TH1D *Matched_WJb_hist = 0;
//	
//
//	TH1D *Unmatched_WJa_zoom_hist = 0;
//	TH1D *Unmatched_WJb_zoom_hist = 0;
//	TH1D *Matched_WJa_zoom_hist = 0;
//	TH1D *Matched_WJb_zoom_hist = 0;
//
//	TH1D *Unmatched_2Wjets_hist = 0;
//	TH1D *Matched_2Wjets_hist = 0;
//	TH1D *Unmatched_2Wjets_zoom_hist = 0;
//	TH1D *Matched_2Wjets_zoom_hist = 0;

	//System Mass Plots
//	TH1F *Mu_Mttbar_hist = 0;
	TH1F *Jets_Mttbar_hist = 0;
	TH1F *nJets_hist = 0;

	//switches
	bool isData_, isTTbar_;
	
	//selectors and helpers
	TTObjectSelector object_selector_; //selects ttbar objects
	TTGenParticleSelector genp_selector_; //selects generator level objects
	TTGenMatcher matcher_; //matches particles on generator level
	TTPermutator permutator_;
	TTBarSolver solver_;
	
	float evt_weight_;
	MCWeightProducer mc_weights_;
	BTagSFProducer btag_sf_;

	//Gen Object Vars
	double DR_LepBHad_ = 0;
	double DR_LepBLep_ = 0;
	double DR_LepWJa_ = 0;
	double DR_LepWJb_ = 0;
	double DR_BHadBLep_ = 0;
	double DR_BHadWJa_ = 0;
	double DR_BHadWJb_ = 0;
	double DR_BLepWJa_ = 0;
	double DR_BLepWJb_ = 0;
	double DR_WJaWJb_ = 0;

	double Pt_Lep_ = 0;
	double Pt_BHad_ = 0;
	double Pt_BLep_ = 0;
	double Pt_WJa_ = 0;
	double Pt_WJb_ = 0;
	
	double Eta_Lep_ = 0;
	double Eta_BHad_ = 0;
	double Eta_BLep_ = 0;
	double Eta_WJa_ = 0;
	double Eta_WJb_ = 0;

	//  double ElISO_;
	double MuISO_ = 0;
	
	double delta_pt_lep_ = 0;
	double delta_pt_BHad_ = 0;
	double delta_pt_BLep_ = 0;
	double delta_pt_WJa_ = 0;
	double delta_pt_WJb_ = 0;
	
	double delta_phi_lep_ = 0;
	double delta_phi_BHad_ = 0;
	double delta_phi_BLep_ = 0;
	double delta_phi_WJa_ = 0;
	double delta_phi_WJb_ = 0;
	
	double delta_y_lep_ = 0;
	double delta_y_BHad_ = 0;
	double delta_y_BLep_ = 0;
	double delta_y_WJa_ = 0;
	double delta_y_WJb_ = 0;
	
	//  TTree *t1 = 0;
	
	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

    public:
	Test_Analyzer(const std::string output_filename):
		AnalyzerBase("Test_Analyzer", output_filename),
		tracker_(),
		object_selector_(),
		permutator_(),
		matcher_(),
		solver_(),
		genp_selector_(TTGenParticleSelector::SelMode::LHE),
		// adding weighting
		mc_weights_(),
		evt_weight_(1.),
		electron_sf_("electron_sf", false),
		muon_sf_("muon_sf")
		{
		//set tracker
			tracker_.use_weight(&evt_weight_);
			//find out which sample are we running on
//			opts::variables_map &values = parser.values();
			opts::variables_map &values = URParser::instance().values();

			string output_file = values["output"].as<std::string>();
			string sample = systematics::get_sample(output_file);
			bool isSignal = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
			isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal;
	
			isData_  = boost::starts_with(sample, "data");
//			isTTbar_ = boost::starts_with(sample, "ttJets");
			Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
//			systematics_ = get_systematics(output_file);
//			int nosys = values["nosys"].as<int>();
//			int nopdf = values["nopdf"].as<int>();
//			if(nosys == 1) {
//			        systematics_ = {systematics::SysShifts::NOSYS};
//			        Logger::log().info() << "DISABLING SYSTEMATICS!" << endl;
//			}
			if(isData_) {
				if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
				else object_selector_.lepton_type(1);
			}
//			pdfs_ = isTTbar_ && (systematics_.size() > 1) && (nopdf == 0);
 
			if(!isData_) mc_weights_.init(sample);

//			mc_weights_.init(sample);
//

			//Init solver
//			string filename = "prob_ttJets.root";
			string filename = "prob_ttJets.permProbComputer.test.root";
			Logger::log().debug() << "solver file: " << filename << endl;
			TFile probfile(DataFile(filename).path().c_str());
			TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
			// if(!td) td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());      
			solver_.Init(td, false, true, true, true, true); //probfile, btag, nusolv,massdis,anghad,anglep
		};
  
	//This method is called once per job at the beginning of the analysis
	//book here your histograms/tree and run every initialization needed
	virtual void begin()
	{
		Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
		
		outFile_.cd();

//		t1 = new TTree("tree1", "Tree with a few variables");

//		t1->Branch("ElISO", &ElISO_, "ElISO_/D");
//		t1->Branch("MuISO", &MuISO_, "MuISO_/D");

  		//histograms
  		int nbins = 22;
  		double mass_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};

		//Gen Obj. Hists
			//DR
		DR_LepBHad_vs_Mtt_hist = new TH2D("DR_LepBHad_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (l, b_{h})", 500, 250., 2000., 100, 0., 5.);
		DR_LepBLep_vs_Mtt_hist = new TH2D("DR_LepBLep_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (l, b_{l})", 500, 250., 2000., 100, 0., 5.);
		DR_LepWJa_vs_Mtt_hist = new TH2D("DR_LepWJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (l, WJa)", 500, 250., 2000., 100, 0., 5.);
		DR_LepWJb_vs_Mtt_hist = new TH2D("DR_LepWJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (l, WJb)", 500, 250., 2000., 100, 0., 5.);
		DR_BHadBLep_vs_Mtt_hist = new TH2D("DR_BHadBLep_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (b_{h}, b_{l})", 500, 250., 2000., 100, 0., 5.);
		DR_BHadWJa_vs_Mtt_hist = new TH2D("DR_BHadWJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (b_{h}, WJa)", 500, 250., 2000., 100, 0., 5.);
		DR_BHadWJb_vs_Mtt_hist = new TH2D("DR_BHadWJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (b_{h}, WJb)", 500, 250., 2000., 100, 0., 5.);
		DR_BLepWJa_vs_Mtt_hist = new TH2D("DR_BLepWJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (b_{l}, WJa)", 500, 250., 2000., 100, 0., 5.);
		DR_BLepWJb_vs_Mtt_hist = new TH2D("DR_BLepWJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (b_{l}, WJb)", 500, 250., 2000., 100, 0., 5.);
		DR_WJaWJb_vs_Mtt_hist = new TH2D("DR_WJaWJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; #Delta R (WJa, WJb)", 500, 250., 2000., 100, 0., 5.);

		DR_LepBHad_hist = new TH1F("DR_LepBHad_hist", "", 100, 0., 5.);
		DR_LepBLep_hist = new TH1F("DR_LepBLep_hist", "", 100, 0., 5.) ;
		DR_LepWJa_hist = new TH1F("DR_LepWJa_hist", "", 100, 0., 5.); 
		DR_LepWJb_hist = new TH1F("DR_LepWJb_hist", "", 100, 0., 5.);
		DR_BHadBLep_hist = new TH1F("DR_BHadBLep_hist", "", 100, 0., 5.);
		DR_BHadWJa_hist = new TH1F("DR_BHadWJa_hist", "", 100, 0., 5.); 
		DR_BHadWJb_hist = new TH1F("DR_BHadWJb_hist", "", 100, 0., 5.); 
		DR_BLepWJa_hist = new TH1F("DR_BLepWJa_hist", "", 100, 0., 5.); 
		DR_BLepWJb_hist = new TH1F("DR_BLepWJb_hist", "", 100, 0., 5.); 
		DR_WJaWJb_hist = new TH1F("DR_WJaWJb_hist", "", 100, 0., 5.); 

			//Pt
		Pt_Lep_vs_Mtt_hist =  new TH2D("Pt_Lep_vs_Mtt_hist", "m_{t#bar t} [GeV]; l p_{t} ", 500, 250., 2000., 100, 0., 300.);
		Pt_BLep_vs_Mtt_hist = new TH2D("Pt_BLep_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{l} p_{t} ", 500, 250., 2000., 100, 0., 300.); 
		Pt_BHad_vs_Mtt_hist = new TH2D("Pt_BHad_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{h} p_{t} ", 500, 250., 2000., 100, 0., 300.); 
		Pt_WJa_vs_Mtt_hist = new TH2D("Pt_WJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJa p_{t} ", 500, 250., 2000., 100, 0., 300.); 
		Pt_WJb_vs_Mtt_hist = new TH2D("Pt_WJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJb p_{t} ", 500, 250., 2000., 100, 0., 300.); 
	
		Pt_Lep_hist = new TH1F("Pt_Lep_hist", "", 100, 0., 300.);
		Pt_BLep_hist = new TH1F("Pt_BLep_hist", "", 100, 0., 300.); 
		Pt_BHad_hist = new TH1F("Pt_BHad_hist", "", 100, 0., 300.); 
		Pt_WJa_hist = new TH1F("Pt_WJa_hist", "", 100, 0., 300.); 
		Pt_WJb_hist = new TH1F("Pt_WJb_hist", "", 100, 0., 300.); 
		  //Eta
		Eta_Lep_vs_Mtt_hist = new TH2D("Eta_Lep_vs_Mtt_hist", "m_{t#bar t} [GeV]; l #eta ", 500, 250., 2000., 100, -3., 3.); 
		Eta_BLep_vs_Mtt_hist = new TH2D("Eta_BLep_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{l} #eta ", 500, 250., 2000., 100, -3., 3.);  
		Eta_BHad_vs_Mtt_hist = new TH2D("Eta_BHad_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{h} #eta ", 500, 250., 2000., 100, -3., 3.);  
		Eta_WJa_vs_Mtt_hist = new TH2D("Eta_WJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJa #eta ", 500, 250., 2000., 100, -3., 3.);  
		Eta_WJb_vs_Mtt_hist = new TH2D("Eta_WJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJb #eta ", 500, 250., 2000., 100, -3., 3.);  
	
		Eta_Lep_hist = new TH1F("Eta_Lep_hist", "", 100, -3, 3.);
		Eta_BLep_hist = new TH1F("Eta_BLep_hist", "", 100, -3, 3.); 
		Eta_BHad_hist = new TH1F("Eta_BHad_hist", "", 100, -3, 3.); 
		Eta_WJa_hist = new TH1F("Eta_WJa_hist", "", 100, -3, 3.); 
		Eta_WJb_hist = new TH1F("Eta_WJb_hist", "", 100, -3, 3.); 

//			//Muon Hists
//  		MuISO_hist = new TH1F("MuISO_hist", "", 500, 0., 10.);
//  		LepBLep_DR_vs_mtt_hist = new TH2D("LepBLep_DR_vs_mtt_hist", "m_{t#bar t} [GeV]; #Delta R (l,b_{l})", 500, 250., 2000., 300, 0., 4.);
//  		MuISO_vs_mtt_hist= new TH2D("MuISO_vs_mtt_hist", "m_{t#bar t} [GeV]; #mu ISO", 500, 250., 2000., 300, 0., 3.);
//  		MuDR_vs_MuRelDPT_hist = new TH2D("MuDR_vs_MuRelDPT_hist", "#mu #Delta R (gen,reco); #mu Rel #Delta P_{t}", 500, 0., 0.4, 500, 0., 40);
//
//  		delta_pt_lep_hist = new TH1F("delta_pt_lep_hist", "", 100, -300., 300.);
//  		delta_phi_lep_hist = new TH1F("delta_phi_lep_hist", "", 100, -1., 1.);
//  		delta_y_lep_hist = new TH1F("delta_y_lep_hist", "", 100, -1., 1.);
//				//Muon Eff Hists
//  		Unmatched_Muons_hist = new TH1D("Unmatched_Muons_hist","", nbins, mass_bins);
//		Unmatched_Muons_hist->Sumw2();
//  		Matched_Muons_hist = new TH1D("Matched_Muons_hist","", nbins, mass_bins);
//		Matched_Muons_hist->Sumw2();
//  		Non_ISO_Mu_hist = new TH1D("Non_ISO_Mu_hist", "", nbins, mass_bins);
//		Non_ISO_Mu_hist->Sumw2();
//  		ISO_Mu_hist = new TH1D("ISO_Mu_hist", "", nbins, mass_bins);
//		ISO_Mu_hist->Sumw2();
//
//			//BJet Hists
//  		delta_pt_BHad_hist = new TH1F("delta_pt_BHad_hist", "", 100, -300., 300.);
//  		delta_pt_BLep_hist = new TH1F("delta_pt_BLep_hist", "", 100, -300., 300.);
//  		delta_phi_BHad_hist = new TH1F("delta_phi_BHad_hist", "", 100, -1., 1.);
//  		delta_phi_BLep_hist = new TH1F("delta_phi_BLep_hist", "", 100, -1., 1.);
//  		delta_y_BHad_hist = new TH1F("delta_y_BHad_hist", "", 100, -1., 1.);
//  		delta_y_BLep_hist = new TH1F("delta_y_BLep_hist", "", 100, -1., 1.);
//
//		Unmatched_BHad_pt_hist = new TH1F("Unmatched_BHad_pt_hist", "", 100, 0., 400.);
//		Unmatched_BLep_pt_hist = new TH1F("Unmatched_BLep_pt_hist", "", 100, 0., 400.);
//		Unmatched_BHad_y_hist = new TH1F("Unmatched_BHad_y_hist", "", 100, -10., 10.);
//		Unmatched_BLep_y_hist = new TH1F("Unmatched_BLep_y_hist", "", 100, -10., 10.);
//
//				//BJet Eff Hists
//  		Unmatched_BHad_hist = new TH1D("Unmatched_BHad_hist","", nbins, mass_bins);
//		Unmatched_BHad_hist->Sumw2();
//  		Unmatched_BLep_hist = new TH1D("Unmatched_BLep_hist","", nbins, mass_bins);
//		Unmatched_BLep_hist->Sumw2();
//  		Matched_BHad_hist = new TH1D("Matched_BHad_hist","", nbins, mass_bins);
//		Matched_BHad_hist->Sumw2();
//  		Matched_BLep_hist = new TH1D("Matched_BLep_hist","", nbins, mass_bins);
//		Matched_BLep_hist->Sumw2();
//
//
//		Unmatched_BHad_zoom_hist = new TH1D("Unmatched_BHad_zoom_hist", "", 10, 250., 450.);
//		Unmatched_BHad_zoom_hist->Sumw2();
//		Unmatched_BLep_zoom_hist = new TH1D("Unmatched_BLep_zoom_hist", "", 10, 250., 450.);
//		Unmatched_BLep_zoom_hist->Sumw2();
//		Matched_BHad_zoom_hist = new TH1D("Matched_BHad_zoom_hist", "", 10, 250., 450.);
//		Matched_BHad_zoom_hist->Sumw2();
//		Matched_BLep_zoom_hist = new TH1D("Matched_BLep_zoom_hist", "", 10, 250., 450.);
//		Matched_BLep_zoom_hist->Sumw2();
//			//WJet Hists
//  		delta_pt_WJa_hist = new TH1F("delta_pt_WJa_hist", "", 100, -300., 300.);
//  		delta_pt_WJb_hist = new TH1F("delta_pt_WJb_hist", "", 100, -300., 300.);
//  		delta_phi_WJa_hist = new TH1F("delta_phi_WJa_hist", "", 100, -1., 1.);
//  		delta_phi_WJb_hist = new TH1F("delta_phi_WJb_hist", "", 100, -1., 1.);
//  		delta_y_WJa_hist = new TH1F("delta_y_WJa_hist", "", 100, -1., 1.);
//  		delta_y_WJb_hist = new TH1F("delta_y_WJb_hist", "", 100, -1., 1.);
//
//		Unmatched_WJa_pt_hist = new TH1F("Unmatched_WJa_pt_hist", "", 100, 0., 400.);
//		Unmatched_WJb_pt_hist = new TH1F("Unmatched_WJb_pt_hist", "", 100, 0., 400.);
//		Unmatched_WJa_y_hist = new TH1F("Unmatched_WJa_y_hist", "", 100, -10., 10.);
//		Unmatched_WJb_y_hist = new TH1F("Unmatched_WJb_y_hist", "", 100, -10., 10.);
//				//WJet Eff Hists
//  		Unmatched_WJa_hist = new TH1D("Unmatched_WJa_hist","", nbins, mass_bins);
//		Unmatched_WJa_hist->Sumw2();
//  		Unmatched_WJb_hist = new TH1D("Unmatched_WJb_hist","", nbins, mass_bins);
//		Unmatched_WJb_hist->Sumw2();
//  		Matched_WJa_hist = new TH1D("Matched_WJa_hist","", nbins, mass_bins);
//		Matched_WJa_hist->Sumw2();
//  		Matched_WJb_hist = new TH1D("Matched_WJb_hist","", nbins, mass_bins);
//		Matched_WJb_hist->Sumw2();
//  	
//		Unmatched_WJa_zoom_hist = new TH1D("Unmatched_WJa_zoom_hist", "", 10, 250., 450.);
//		Unmatched_WJa_zoom_hist->Sumw2();
//		Unmatched_WJb_zoom_hist = new TH1D("Unmatched_WJb_zoom_hist", "", 10, 250., 450.);
//		Unmatched_WJb_zoom_hist->Sumw2();
//		Matched_WJa_zoom_hist = new TH1D("Matched_WJa_zoom_hist", "", 10, 250., 450.);
//		Matched_WJa_zoom_hist->Sumw2();
//		Matched_WJb_zoom_hist = new TH1D("Matched_WJb_zoom_hist", "", 10, 250., 450.);
//		Matched_WJb_zoom_hist->Sumw2();
//
//		Unmatched_2Wjets_hist = new TH1D("Unmatched_2Wjets_hist", "", nbins, mass_bins);
//		Unmatched_2Wjets_hist->Sumw2();
//		Matched_2Wjets_hist = new TH1D("Matched_2Wjets_hist", "", nbins, mass_bins);
//		Matched_2Wjets_hist->Sumw2();
//
//		Unmatched_2Wjets_zoom_hist = new TH1D("Unmatched_2Wjets_zoom_hist", "", 10, 250., 450.);
//		Unmatched_2Wjets_zoom_hist->Sumw2();
//		Matched_2Wjets_zoom_hist = new TH1D("Matched_2Wjets_zoom_hist", "", 10, 250., 450.);
//		Matched_2Wjets_zoom_hist->Sumw2();
//
//		//System Masses
//		Mu_Mttbar_hist = new TH1F("Mu_Mttbar_hist", "", 500, 250., 2000.);
		Jets_Mttbar_hist = new TH1F("Jets_Mttbar_hist", "", 500, 250., 2000.);
		nJets_hist = new TH1F("nJets_hist", "", 10, 0., 10.);

//		
  		Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
	}

	//This method is called once every file, contains the event loop
	//run your proper analysis here
	virtual void analyze()
	{
		Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;
		
		opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();
		int report = values["report"].as<int>();
		if(evt_idx_ >= limit) return;
		
		URStreamer event(tree_);

		while(event.next() /*&& evt_idx_ < 80000*/)
		{
			if(limit > 0 && evt_idx_ > limit) return;
			evt_idx_++;
			if(skip > 0 && evt_idx_ < skip) continue;
			if(evt_idx_ % 10000 == 1) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

			Logger::log().debug() << "njets: " << object_selector_.clean_jets().size() << endl;
			
			//generator selection
			bool selection = genp_selector_.select(event);
			tracker_.track("gen selection");
			if( !selection ){
				Logger::log().debug() << "event has no selection " << endl;
				continue;
			}
			
			//full ttbar event selection trial (reco objects)
			if( !object_selector_.select(event) ) continue;
			tracker_.track("obj selection");
			
			if( !permutator_.preselection(				
				object_selector_.clean_jets(), 
				object_selector_.lepton(),
				object_selector_.met(),
				object_selector_.lepton_charge() ) ) continue;
			
			//find mc weights
			if(!isData_) {
				if(object_selector_.tight_muons().size() == 1)
					evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
				if(object_selector_.medium_electrons().size() == 1)
					evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
			}

			//Generator level matching (Gen matching)
//			Permutation matched_perm;
			GenTTBar &ttbar = genp_selector_.ttbar_system();
//			if( ttbar.type == GenTTBar::DecayType::SEMILEP){ //gets correct kind of events
//				matched_perm = matcher_.match( //needs 4 arguments
//					genp_selector_.ttbar_final_system(),
//					object_selector_.clean_jets(),
//					object_selector_.loose_electrons(),
//					object_selector_.loose_muons()
//					);
//			}
//			matched_perm.SetMET(object_selector_.met());
////			if( !matched_perm.IsComplete() ) continue;	
////			if( !matched_perm.WJa()->match() || !matched_perm.WJb()->match() ){
////				Logger::log().debug() << "matched_perm has no match " << endl;
////				continue;
////			}
////			matched_perm.Solve(solver_);
////			tracker_.track("gen matching");
////			++gen_match_;
//
//			//Find best permutation	
//			bool go_on = true;
//			Permutation best_permutation;
//			size_t ncycles_=0;
//			while(go_on) {
//				ncycles_++; 
//				Permutation test_perm = permutator_.next(go_on);
//				if(go_on) {
//					test_perm.Solve(solver_);
//					//if(test_perm.MassDiscr() < best_permutation.MassDiscr()){
//					if(test_perm.Prob() < best_permutation.Prob()){
//						best_permutation = test_perm;
//					}
//				}
//			}
//			if(!best_permutation.IsComplete() || best_permutation.Prob() > 1E9) continue; //FIXME, is right??? best_permutation.Prob() > 1E9
//			tracker_.track("best perm");

			// Reco-Gen vars/plots
				//require at least 3 jets in final state
			if( object_selector_.clean_jets().size() < 3 ) continue;
//			if( object_selector_.clean_jets().size() == 3 ){ cout << "njets = 3" << endl;}

			//Generator level matching for ID muons
			const GenObject* lepton = ttbar.lepton();
			const GenObject* lep_b = ttbar.lep_b();

//			if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
//
//				//match gen and reco muons
//				const vector<Muon>& muons = event.muons();
//				if( muons.size() == 0 ) continue;
//		
//				list<IDMuon> mu_list;
//				IDMuon *best_mu = 0;
//				float dr = 1e10; // initialize delta r to a really high number
//		
//				for( vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon ){
//		
//					IDMuon mu(*muon, event.rho().value());
//					if( mu.DeltaR(*lepton) > 0.4 ) continue; //cut on DR cone of 0.4
//					if( Abs(lepton->Pt() - mu.Pt())/mu.Pt() > 0.2*Abs(lepton->Pt() - mu.Pt()) ) continue;//Rel DPt < 20% of DPt
//					if( mu.DeltaR(*lepton) < dr ){
//						dr = mu.DeltaR(*lepton);
//						mu_list.push_back(mu);
//						best_mu = &(mu_list.back());
//					}
//				}
//		
//		//		double ttmass = genp_selector_.ttbar_final_system().M();
//				double ttmass = genp_selector_.ttbar_system().M();
//				if( best_mu == 0 ){
//					Unmatched_Muons_hist->Fill(ttmass);
//					continue;
//				}
//				best_mu->addMatch(lepton);
//				MuISO_ = best_mu->RelPFIsoDb();
//		
//				if( MuISO_ > 0.15 ){
//					Non_ISO_Mu_hist->Fill(ttmass);
//				}
//				if( MuISO_ <= 0.15 ){
//					ISO_Mu_hist->Fill(ttmass);
//				}
//		
//				Matched_Muons_hist->Fill(ttmass);
//
//				MuISO_hist->Fill(MuISO_);
//				MuDR_vs_MuRelDPT_hist->Fill(best_mu->DeltaR(*lepton), Abs(lepton->Pt() - best_mu->Pt())/best_mu->Pt());
//				MuISO_vs_mtt_hist->Fill(ttmass, MuISO_);
//				LepBLep_DR_vs_mtt_hist->Fill(ttmass, lepton->DeltaR(*lep_b));
//		
//				delta_pt_lep_ = lepton->Pt() - best_mu->Pt();
//				delta_phi_lep_ = lepton->Phi() - best_mu->Phi();
//				delta_y_lep_ = lepton->Rapidity() - best_mu->Rapidity();
//		
//				if( delta_phi_lep_ > M_PI ){
//					delta_phi_lep_ = delta_phi_lep_ - 2*M_PI;}
//				if( delta_phi_lep_ < -M_PI ){
//					delta_phi_lep_ = delta_phi_lep_ + 2*M_PI;}
//
//				Mu_Mttbar_hist->Fill(ttmass);
//
//				delta_pt_lep_hist->Fill(delta_pt_lep_);//reco-gen
//				delta_phi_lep_hist->Fill(delta_phi_lep_);//reco-gen
//				delta_y_lep_hist->Fill(delta_y_lep_);//reco-gen
//			}

			// Jet matching things
			const GenObject* had_b = ttbar.had_b();

			if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
//				if( !ttbar.had_W()->first || !ttbar.had_W()->second ){
//					cout << "Hi" << endl;
//				}
//					//BJets
//				if( !matched_perm.BHad() && !matched_perm.BLep() ){ //both Bs unmatched
//					Unmatched_BHad_hist->Fill(ttbar.M());
//					Unmatched_BLep_hist->Fill(ttbar.M());
//					if( ttbar.M() < 450 ){
//						Unmatched_BHad_zoom_hist->Fill(ttbar.M());
//						Unmatched_BLep_zoom_hist->Fill(ttbar.M());
//	
//						if( had_b ) {
//							Unmatched_BHad_pt_hist->Fill(had_b->Pt());
//							Unmatched_BHad_y_hist->Fill(had_b->Rapidity());
//						}
//						if( lep_b ){
//							Unmatched_BLep_pt_hist->Fill(lep_b->Pt());
//							Unmatched_BLep_y_hist->Fill(lep_b->Rapidity());
//						}
//					}
//				}
//				if( !matched_perm.BLep() && matched_perm.BHad() ){ //only BHad matched
//					Unmatched_BLep_hist->Fill(ttbar.M());
//	
//					delta_pt_BHad_ = matched_perm.BHad()->match()->Pt()-matched_perm.BHad()->Pt();
//					delta_phi_BHad_ = matched_perm.BHad()->match()->Phi()-matched_perm.BHad()->Phi();
//					delta_y_BHad_ = matched_perm.BHad()->match()->Rapidity()-matched_perm.BHad()->Rapidity();
//					Matched_BHad_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450 ){
//						Unmatched_BLep_zoom_hist->Fill(ttbar.M());
//						Matched_BHad_zoom_hist->Fill(ttbar.M());
//	
//						if( lep_b ){
//							Unmatched_BLep_pt_hist->Fill(lep_b->Pt());
//							Unmatched_BLep_y_hist->Fill(lep_b->Rapidity());
//						}
//					}
//				}
//				if( !matched_perm.BHad() && matched_perm.BLep() ){ //only BLep matched
//					Unmatched_BHad_hist->Fill(ttbar.M());
//	
//					delta_pt_BLep_ = matched_perm.BLep()->match()->Pt()-matched_perm.BLep()->Pt();
//					delta_phi_BLep_ = matched_perm.BLep()->match()->Phi()-matched_perm.BLep()->Phi();
//					delta_y_BLep_ = matched_perm.BLep()->match()->Rapidity()-matched_perm.BLep()->Rapidity();
//					Matched_BLep_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450 ){
//						Unmatched_BHad_zoom_hist->Fill(ttbar.M());
//						Matched_BLep_zoom_hist->Fill(ttbar.M());
//	
//						if( had_b ){
//							Unmatched_BHad_pt_hist->Fill(had_b->Pt());
//							Unmatched_BHad_y_hist->Fill(had_b->Rapidity());
//						}
//					}
//				}
//				if( matched_perm.BHad() && matched_perm.BLep() ){ //both Bs matched
//					if( (matched_perm.BHad()->E() == matched_perm.BLep()->E()) ) continue; //check reco objects aren't the same
//					delta_pt_BLep_ = matched_perm.BLep()->match()->Pt()-matched_perm.BLep()->Pt();
//					delta_phi_BLep_ = matched_perm.BLep()->match()->Phi()-matched_perm.BLep()->Phi();
//					delta_y_BLep_ = matched_perm.BLep()->match()->Rapidity()-matched_perm.BLep()->Rapidity();
//					Matched_BLep_hist->Fill(ttbar.M());
//					
//					delta_pt_BHad_ = matched_perm.BHad()->match()->Pt()-matched_perm.BHad()->Pt();
//					delta_phi_BHad_ = matched_perm.BHad()->match()->Phi()-matched_perm.BHad()->Phi();
//					delta_y_BHad_ = matched_perm.BHad()->match()->Rapidity()-matched_perm.BHad()->Rapidity();
//					Matched_BHad_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450){
//						Matched_BLep_zoom_hist->Fill(ttbar.M());
//						Matched_BHad_zoom_hist->Fill(ttbar.M());
//					}
//				}
					//WJets
				const GenObject* WJetA = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
				const GenObject* WJetB = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;

//				if( !matched_perm.WJa() && !matched_perm.WJb() ){//both Wjets unmatched
//					Unmatched_WJa_hist->Fill(ttbar.M());
//					Unmatched_WJb_hist->Fill(ttbar.M());
//					Unmatched_2Wjets_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450 ){
//						Unmatched_WJa_zoom_hist->Fill(ttbar.M());
//						Unmatched_WJb_zoom_hist->Fill(ttbar.M());
//						Unmatched_2Wjets_zoom_hist->Fill(ttbar.M());
//		
//						if( ttbar.had_W()->first && ttbar.had_W()->second ){//check both gen objects exist
//							const GenObject* WjetA = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
//							const GenObject* WjetB = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
//							Unmatched_WJa_pt_hist->Fill(WjetA->Pt());
//							Unmatched_WJa_y_hist->Fill(WjetA->Rapidity());
//	
//							Unmatched_WJb_pt_hist->Fill(WjetB->Pt());
//							Unmatched_WJb_y_hist->Fill(WjetB->Rapidity());
//						}
//						if( ttbar.had_W()->first && !ttbar.had_W()->second ){//if only first gen W exists
//							Unmatched_WJa_pt_hist->Fill(ttbar.had_W()->first->Pt());
//							Unmatched_WJa_y_hist->Fill(ttbar.had_W()->first->Rapidity());
//						}
//						if( !ttbar.had_W()->first && ttbar.had_W()->second ){//if only second gen W exists
//							Unmatched_WJb_pt_hist->Fill(ttbar.had_W()->second->Pt());
//							Unmatched_WJb_y_hist->Fill(ttbar.had_W()->second->Rapidity());
//						}
//					}
//				}
//				if( !matched_perm.WJa() && matched_perm.WJb() ){//only WJb matched
//					Unmatched_WJa_hist->Fill(ttbar.M());
//					Unmatched_2Wjets_hist->Fill(ttbar.M());
//	
//					delta_pt_WJb_ = matched_perm.WJb()->match()->Pt()-matched_perm.WJb()->Pt();
//					delta_phi_WJb_ = matched_perm.WJb()->match()->Phi()-matched_perm.WJb()->Phi();
//					delta_y_WJb_ = matched_perm.WJb()->match()->Rapidity()-matched_perm.WJb()->Rapidity();
//					Matched_WJb_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450 ){
//						Unmatched_WJa_zoom_hist->Fill(ttbar.M());
//						Unmatched_2Wjets_zoom_hist->Fill(ttbar.M());
//						Matched_WJb_zoom_hist->Fill(ttbar.M());
//		
//						if( ttbar.had_W()->first && ttbar.had_W()->second ){//check both gen objects exist
//							const GenObject* WjetA = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
//							const GenObject* WjetB = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
//							if( matched_perm.WJb()->match()->E() == WjetA->E() ){//if WJb matched to WJetA
//								Unmatched_WJa_pt_hist->Fill(WjetB->Pt());
//								Unmatched_WJa_y_hist->Fill(WjetB->Rapidity());
//							}
//							if( matched_perm.WJb()->match()->E() == WjetB->E() ){//if WJb matched to WjetB
//								Unmatched_WJa_pt_hist->Fill(WjetA->Pt());
//								Unmatched_WJa_y_hist->Fill(WjetA->Rapidity());
//							}
//						}
//						if( ttbar.had_W()->first && !ttbar.had_W()->second ){//if only first gen W exists
//							if( !(matched_perm.WJb()->match()->E() == ttbar.had_W()->first->E()) ){//WJb not matched to W first
//								Unmatched_WJa_pt_hist->Fill(ttbar.had_W()->first->Pt());
//								Unmatched_WJa_y_hist->Fill(ttbar.had_W()->first->Rapidity());
//							}
//						}
//						if( !ttbar.had_W()->first && ttbar.had_W()->second ){//if only second gen W exists
//							if( !(matched_perm.WJb()->match()->E() == ttbar.had_W()->second->E()) ){//WJb not matched to W second
//								Unmatched_WJa_pt_hist->Fill(ttbar.had_W()->second->Pt());
//								Unmatched_WJa_y_hist->Fill(ttbar.had_W()->second->Rapidity());
//							}
//						}
//					}
//				}
//				if( matched_perm.WJa() && !matched_perm.WJb() ){// only WJa matched
//					Unmatched_WJb_hist->Fill(ttbar.M());
//					Unmatched_2Wjets_hist->Fill(ttbar.M());
//	
//					delta_pt_WJa_ = matched_perm.WJa()->match()->Pt()-matched_perm.WJa()->Pt();
//					delta_phi_WJa_ = matched_perm.WJa()->match()->Phi()-matched_perm.WJa()->Phi();
//					delta_y_WJa_ = matched_perm.WJa()->match()->Rapidity()-matched_perm.WJa()->Rapidity();
//					Matched_WJa_hist->Fill(ttbar.M());
//	
//					if( ttbar.M() < 450 ){
//						Unmatched_WJb_zoom_hist->Fill(ttbar.M());
//						Unmatched_2Wjets_zoom_hist->Fill(ttbar.M());
//						Matched_WJa_zoom_hist->Fill(ttbar.M());
//		
//						if( ttbar.had_W()->first && ttbar.had_W()->second ){//check both gen objects exist
//							const GenObject* WjetA = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
//							const GenObject* WjetB = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
//							if( matched_perm.WJa()->match()->E() == WjetA->E() ){//if WJa matched to WJetA
//								Unmatched_WJb_pt_hist->Fill(WjetB->Pt());
//								Unmatched_WJb_y_hist->Fill(WjetB->Rapidity());
//							}
//							if( matched_perm.WJa()->match()->E() == WjetB->E() ){//if WJa matched to WjetB
//								Unmatched_WJb_pt_hist->Fill(WjetA->Pt());
//								Unmatched_WJb_y_hist->Fill(WjetA->Rapidity());
//							}
//						}
//						if( ttbar.had_W()->first && !ttbar.had_W()->second ){//if only first gen W exists
//							if( !(matched_perm.WJa()->match()->E() == ttbar.had_W()->first->E()) ){//WJa not matched to W first
//								Unmatched_WJb_pt_hist->Fill(ttbar.had_W()->first->Pt());
//								Unmatched_WJb_y_hist->Fill(ttbar.had_W()->first->Rapidity());
//							}
//						}
//						if( !ttbar.had_W()->first && ttbar.had_W()->second ){//if only second gen W exists
//							if( !(matched_perm.WJa()->match()->E() == ttbar.had_W()->second->E()) ){//WJa not matched to W second
//								Unmatched_WJb_pt_hist->Fill(ttbar.had_W()->second->Pt());
//								Unmatched_WJb_y_hist->Fill(ttbar.had_W()->second->Rapidity());
//							}
//						}
//					}
//				}
//				if( matched_perm.WJa() && matched_perm.WJb() ){//both Wjets matched
//					if( matched_perm.WJa()->E() == matched_perm.WJb()->E() ) continue;//check reco objects aren't same
//					Matched_2Wjets_hist->Fill(ttbar.M());
//	
//					IDJet* WJa = (matched_perm.WJa()->E() > matched_perm.WJb()->E()) ? matched_perm.WJa() : matched_perm.WJb();
//					IDJet* WJb = (matched_perm.WJa()->E() > matched_perm.WJb()->E()) ? matched_perm.WJb() : matched_perm.WJa();
//	
//					delta_pt_WJa_ = WJa->match()->Pt()-WJa->Pt();
//					delta_phi_WJa_ = WJa->match()->Phi()-WJa->Phi();
//					delta_y_WJa_ = WJa->match()->Rapidity()-WJa->Rapidity();
//					Matched_WJa_hist->Fill(ttbar.M());
//	
//					delta_pt_WJb_ = WJb->match()->Pt()-WJb->Pt();
//					delta_phi_WJb_ = WJb->match()->Phi()-WJb->Phi();
//					delta_y_WJb_ = WJb->match()->Rapidity()-WJb->Rapidity();
//					Matched_WJb_hist->Fill(ttbar.M());
//					
//					if( ttbar.M() < 450 ){
//						Matched_WJa_zoom_hist->Fill(ttbar.M());
//						Matched_WJb_zoom_hist->Fill(ttbar.M());
//						Matched_2Wjets_zoom_hist->Fill(ttbar.M());
//					}
//				}
//
//					//restrict delta phi between -Pi to Pi
//				if( delta_phi_BHad_ > M_PI ){
//					delta_phi_BHad_ = delta_phi_BHad_ - 2*M_PI;}
//				if( delta_phi_BHad_ < -M_PI ){
//					delta_phi_BHad_ = delta_phi_BHad_ + 2*M_PI;}
//				if( delta_phi_BLep_ > M_PI ){
//					delta_phi_BLep_ = delta_phi_BLep_ - 2*M_PI;}
//				if( delta_phi_BLep_ < -M_PI ){
//					delta_phi_BLep_ = delta_phi_BLep_ + 2*M_PI;}
//				if( delta_phi_WJa_ > M_PI ){
//					delta_phi_WJa_ = delta_phi_WJa_ - 2*M_PI;}
//				if( delta_phi_WJa_ < -M_PI ){
//					delta_phi_WJa_ = delta_phi_WJa_ + 2*M_PI;}
//				if( delta_phi_WJb_ > M_PI ){
//					delta_phi_WJb_ = delta_phi_WJb_ - 2*M_PI;}
//				if( delta_phi_WJb_ < -M_PI ){
//					delta_phi_WJb_ = delta_phi_WJb_ + 2*M_PI;}
//			
//					//check range of rapidity
//				if( Abs(delta_y_BHad_) > 2 ){
//					cout << "Delta BHad Rapidity: " << delta_y_BHad_ << endl;
//				}
//				if( Abs(delta_y_BLep_) > 2 ){
//					cout << "Delta BLep Rapidity: " << delta_y_BLep_ << endl;
//				}
//				if( Abs(delta_y_WJa_) > 2 ){
//					cout << "Delta WJa Rapidity: " << delta_y_WJa_ << endl;
//				}
//				if( Abs(delta_y_WJb_) > 2 ){
//					cout << "Delta WJb Rapidity: " << delta_y_WJb_ << endl;
//				}
//			
//				// Fill Histograms
//				delta_pt_BHad_hist->Fill(delta_pt_BHad_); //reco-gen
//				delta_pt_BLep_hist->Fill(delta_pt_BLep_);
//				delta_pt_WJa_hist->Fill(delta_pt_WJa_);
//				delta_pt_WJb_hist->Fill(delta_pt_WJb_);
//			
//				delta_phi_BHad_hist->Fill(delta_phi_BHad_);
//				delta_phi_BLep_hist->Fill(delta_phi_BLep_);
//				delta_phi_WJa_hist->Fill(delta_phi_WJa_);
//				delta_phi_WJb_hist->Fill(delta_phi_WJb_);
//			
//				delta_y_BHad_hist->Fill(delta_y_BHad_);
//				delta_y_BLep_hist->Fill(delta_y_BLep_);
//				delta_y_WJa_hist->Fill(delta_y_WJa_);
//				delta_y_WJb_hist->Fill(delta_y_WJb_);

				//Gen Obj Vars.
					//DR
				DR_LepBHad_ = lepton->DeltaR(*had_b);
				DR_LepBLep_ = lepton->DeltaR(*lep_b);
				DR_LepWJa_ = lepton->DeltaR(*WJetA);
				DR_LepWJb_ = lepton->DeltaR(*WJetB);
				DR_BHadBLep_ = had_b->DeltaR(*lep_b);
				DR_BHadWJa_ = had_b->DeltaR(*WJetA);
				DR_BHadWJb_ = had_b->DeltaR(*WJetB);
				DR_BLepWJa_ = lep_b->DeltaR(*WJetA);
				DR_BLepWJb_ = lep_b->DeltaR(*WJetB);
				DR_WJaWJb_ = WJetA->DeltaR(*WJetB);
					//Pt
				Pt_Lep_ = lepton->Pt();
				Pt_BHad_ = had_b->Pt();
				Pt_BLep_ = lep_b->Pt();
				Pt_WJa_ = WJetA->Pt();
				Pt_WJb_ = WJetB->Pt();
					//Eta
				Eta_Lep_ = lepton->Eta();
				Eta_BHad_ = had_b->Eta();
				Eta_BLep_ = lep_b->Eta();
				Eta_WJa_ = WJetA->Eta();
				Eta_WJb_ = WJetB->Eta();
				
				//Fill Hists
				DR_LepBHad_vs_Mtt_hist->Fill(ttbar.M(), DR_LepBHad_);
				DR_LepBLep_vs_Mtt_hist->Fill(ttbar.M(), DR_LepBLep_);
				DR_LepWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_LepWJa_);
				DR_LepWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_LepWJb_);
				DR_BHadBLep_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadBLep_);
				DR_BHadWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadWJa_);
				DR_BHadWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadWJb_);
				DR_BLepWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_BLepWJa_);
				DR_BLepWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_BLepWJb_);
				DR_WJaWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_WJaWJb_);
		
				DR_LepBHad_hist->Fill(DR_LepBHad_);
				DR_LepBLep_hist->Fill(DR_LepBLep_);
				DR_LepWJa_hist->Fill(DR_LepWJa_);
				DR_LepWJb_hist->Fill(DR_LepWJb_);
				DR_BHadBLep_hist->Fill(DR_BHadBLep_);
				DR_BHadWJa_hist->Fill(DR_BHadWJa_);
				DR_BHadWJb_hist->Fill(DR_BHadWJb_);
				DR_BLepWJa_hist->Fill(DR_BLepWJa_);
				DR_BLepWJb_hist->Fill(DR_BLepWJb_);
				DR_WJaWJb_hist->Fill(DR_WJaWJb_);
		
				Pt_Lep_vs_Mtt_hist->Fill(ttbar.M(), Pt_Lep_);
				Pt_BLep_vs_Mtt_hist->Fill(ttbar.M(), Pt_BLep_);
				Pt_BHad_vs_Mtt_hist->Fill(ttbar.M(), Pt_BHad_);
				Pt_WJa_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJa_);
				Pt_WJb_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJb_);
			
				Pt_Lep_hist->Fill(Pt_Lep_);
				Pt_BLep_hist->Fill(Pt_BLep_);
				Pt_BHad_hist->Fill(Pt_BHad_);
				Pt_WJa_hist->Fill(Pt_WJa_);
				Pt_WJb_hist->Fill(Pt_WJb_);
		
				Eta_Lep_vs_Mtt_hist->Fill(ttbar.M(), Eta_Lep_);
				Eta_BLep_vs_Mtt_hist->Fill(ttbar.M(), Eta_BLep_);
				Eta_BHad_vs_Mtt_hist->Fill(ttbar.M(), Eta_BHad_);
				Eta_WJa_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJa_);
				Eta_WJb_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJb_);
			
				Eta_Lep_hist->Fill(Eta_Lep_);
				Eta_BLep_hist->Fill(Eta_BLep_);
				Eta_BHad_hist->Fill(Eta_BHad_);
				Eta_WJa_hist->Fill(Eta_WJa_);
				Eta_WJb_hist->Fill(Eta_WJb_);
				
				Jets_Mttbar_hist->Fill(ttbar.M());
				nJets_hist->Fill(object_selector_.clean_jets().size());
			}
		//
		////	Mttbar_hist->Fill(ttang.M());
//			if( !matched_perm.IsComplete() ) continue;	
//			if( !matched_perm.WJa()->match() || !matched_perm.WJb()->match() ){
//				Logger::log().debug() << "matched_perm has no match " << endl;
//				continue;
//			}
//			matched_perm.Solve(solver_);
//			tracker_.track("gen matching");
//			++gen_match_;

//			//make sure one event is even and the other is odd
//			if( (Abs(matched_perm.WJa()->match()->pdgId()) + Abs(matched_perm.WJb()->match()->pdgId()) % 2 == 0) ){
//				Logger::log().debug() << "one jet isn't even and the other isn't odd " << endl;
//				continue;
//			}
//		
//			//select up and down type jets
		//	const IDJet* up_Wjet = ( Abs(matched_perm.WJa()->match()->pdgId()) % 2 == 0 ) ? matched_perm.WJa() : matched_perm.WJb();// up type
		//	tracker_.track("Up-type WJets");
		//
		//	const IDJet* down_Wjet = ( Abs(matched_perm.WJa()->match()->pdgId()) % 2 == 1 ) ? matched_perm.WJa() : matched_perm.WJb();// down type
		//	tracker_.track("Down-type WJets");
		
//			if( !(Abs(matched_perm.WJa()->match()->pdgId()) % 2 == 0) ) cout << "WJa not utype!!!!!!" << endl;
//		//	cout << "WJa pdgId: " << matched_perm.WJa()->match()->pdgId() << endl;
//			if( !(Abs(matched_perm.WJb()->match()->pdgId()) % 2 == 1) ) cout << "WJb not dtype!!!!!!" << endl;
		
		
		
		//		//require DeltaR >= 0.4 between all matched jets
		//	if(matched_perm.WJa()->match()->DeltaR(*(matched_perm.WJb()->match())) < 0.4 ) continue;
		////	if(matched_perm.WJa()->match()->DeltaR(*(matched_perm.WJb()->match())) < 0.4 ){
		////		cout << "DeltaR WJa->WJb: " << matched_perm.WJa()->match()->DeltaR(*(matched_perm.WJb()->match())) << endl;}
		//	if( matched_perm.WJa()->match()->DeltaR(*(matched_perm.BHad()->match())) < 0.4 ) continue;
		////	if( matched_perm.WJa()->match()->DeltaR(*(matched_perm.BHad()->match())) < 0.4 ){	
		////		cout << "DeltaR WJa->BHad: " << matched_perm.WJa()->match()->DeltaR(*(matched_perm.BHad()->match())) << endl;}
		//	if( matched_perm.WJa()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ) continue;
		////	if( matched_perm.WJa()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ){	
		////		cout << "DeltaR WJa->BLep: " << matched_perm.WJa()->match()->DeltaR(*(matched_perm.BLep()->match())) << endl;}
		//	if( matched_perm.WJb()->match()->DeltaR(*(matched_perm.BHad()->match())) < 0.4 ) continue;
		////	if( matched_perm.WJb()->match()->DeltaR(*(matched_perm.BHad()->match())) < 0.4 ){
		////		cout << "DeltaR WJb->BHad: " << matched_perm.WJb()->match()->DeltaR(*(matched_perm.BHad()->match())) << endl;}
		//	if( matched_perm.WJb()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ) continue;
		////	if( matched_perm.WJb()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ){
		////		cout << "DeltaR WJb->BLep: " << matched_perm.WJb()->match()->DeltaR(*(matched_perm.BLep()->match())) << endl;}
		//	if( matched_perm.BHad()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ) continue;
		////	if( matched_perm.BHad()->match()->DeltaR(*(matched_perm.BLep()->match())) < 0.4 ){
		////		cout << "DeltaR BHad->BLep: " << matched_perm.BHad()->match()->DeltaR(*(matched_perm.BLep()->match())) << endl;}
		//
		//	delta_pt_BHad_ = matched_perm.BHad()->match()->Pt()-matched_perm.BHad()->Pt();
		//	delta_pt_BLep_ = matched_perm.BLep()->match()->Pt()-matched_perm.BLep()->Pt();
		//	delta_pt_WJa_ = matched_perm.WJa()->match()->Pt()-matched_perm.WJa()->Pt();
		//	delta_pt_WJb_ = matched_perm.WJb()->match()->Pt()-matched_perm.WJb()->Pt();
		//
		//	delta_phi_BHad_ = matched_perm.BHad()->match()->Phi()-matched_perm.BHad()->Phi();
		//	delta_phi_BLep_ = matched_perm.BLep()->match()->Phi()-matched_perm.BLep()->Phi();
		//	delta_phi_WJa_ = matched_perm.WJa()->match()->Phi()-matched_perm.WJa()->Phi();
		//	delta_phi_WJb_ = matched_perm.WJb()->match()->Phi()-matched_perm.WJb()->Phi();
		//
		//	delta_y_BHad_ = matched_perm.BHad()->match()->Rapidity()-matched_perm.BHad()->Rapidity();
		//	delta_y_BLep_ = matched_perm.BLep()->match()->Rapidity()-matched_perm.BLep()->Rapidity();
		//	delta_y_WJa_ = matched_perm.WJa()->match()->Rapidity()-matched_perm.WJa()->Rapidity();
		//	delta_y_WJb_ = matched_perm.WJb()->match()->Rapidity()-matched_perm.WJb()->Rapidity();
		//
		}
		Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
	}

	//this method is called at the end of the job, by default saves
	//every histogram/tree produced, override it if you need something more
	virtual void end()
	{
	//	t1->Write();
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
  URDriver<Test_Analyzer> test;
  int thing = test.run();
  // auto files = gROOT->GetListOfFiles() ;
  // for(int i=0; i<files->GetSize(); i++) {
  //   TNamed *obj = (TNamed*) files->At(i);
  //   Logger::log().debug() << "file " << i << " " << obj->GetName() << std::endl;
  //   TFile *ff = (TFile*) obj;
  //   ff->Close();
  // }
  // gROOT->CloseFiles();
  // files = gROOT->GetListOfFiles();
  // Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl;
  return thing;
}

