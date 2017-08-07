#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include <map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include <algorithm>
//#include "NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include <list>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include <set>
#include "Analyses/URTTbar/interface/IDMet.h"
//#include "Analyses/URTTbar/interface/JetScaler.h"
#include "TUUID.h"   
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"

#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"
#include <sstream>
#include <unordered_map>
#include "URAnalysis/AnalysisFW/interface/EventList.h"
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "TROOT.h"
//#include <fstream>
//#include "JetMETCorrections/Modules/interface/JetResolution.h"
//#include "Analyses/URTTbar/interface/JERFile.h"

using namespace TMath;
using namespace systematics;
typedef SysShifts Sys;

class htt_simple : public AnalyzerBase
{
public:
	//enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
	struct SyncInfo {
		int Run;
		long int LumiSection;
		long int Event; 
		bool RecoSuccess; 
		float MassTT;
		bool hasMuon;
	};

private:
	//histograms and helpers
	unordered_map<string, unordered_map<string, RObject> > histos_;
	CutFlowTracker tracker_;

	//switches
	bool isTTbar_, isSignal_, isData_, runsys_, has_pdfs_, optim_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
	TTGenMatcher matcher_;
	TTPermutator permutator_;
  TTObjectSelector object_selector_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  BTagSFProducer btag_sf_;
  TTBarSolver solver_;
  PDFuncertainty pdf_uncs_;

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

//    JME::JetResolution pt_jet_res_; // initialize jet energy pt resolution
//    JME::JetResolutionScaleFactor jet_res_sf_; // initialize jet energy resolution scale factors


  unsigned long evt_idx_ = 0;
	vector<systematics::SysShifts> systematics_;

  //cuts
  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;

	//sync
	bool sync_;
	TTree *sync_tree_;
	SyncInfo sync_info_;

public:
	inline double MT(TLorentzVector *l, TLorentzVector *met) {
		return sqrt(pow(l->Pt() + met->Pt(), 2) - pow(l->Px() + met->Px(), 2) - pow(l->Py() + met->Py(), 2));
	}

  htt_simple(const std::string output_filename):
    AnalyzerBase("htt_simple", output_filename), 
		tracker_(),
		genp_selector_(),
		matcher_(),
		permutator_(),
    object_selector_(),
    mc_weights_(),
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_(),    
		btag_sf_("permutations.tightb", "permutations.looseb"),
    solver_(true),
		pdf_uncs_(254),
		systematics_(),
		sync_(false),
		sync_tree_(0),
		sync_info_()
	{
		object_selector_.allow_loose();//allow loose but not tight lepton events
		Logger::log().debug() << "htt_simple ctor" << endl;
    cut_tight_b_ = btag_sf_.tight_cut();
    cut_loose_b_ = btag_sf_.loose_cut();
    //cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << 
    //  " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;
    
    URParser &parser = URParser::instance();

//    parser.addCfgParameter<string>("JERC", "JER_SF","");
//    parser.addCfgParameter<string>("JERC", "PT_JER","");
//    parser.parseArguments();
//
//    DataFile jer_sf_fname_(parser.getCfgPar<string>("JERC","JER_SF"));
//    DataFile pt_jer_fname_(parser.getCfgPar<string>("JERC","PT_JER"));
//
//    jet_res_sf_ = JME::JetResolutionScaleFactor(jer_sf_fname_.path());
//    pt_jet_res_ = JME::JetResolution(pt_jer_fname_.path());

		TUUID id;  
		randomizer_.SetSeed(id.Hash());   

    //find out which sample are we running on
    opts::variables_map &values = parser.values();
		string output_file = values["output"].as<std::string>();
		string sample = systematics::get_sample(output_file);
    isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
		isTTbar_ = boost::starts_with(sample, "ttJets");
		isData_  = boost::starts_with(sample, "data");

    //set tracker
		if(!values.count("noweights")) tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
		permutator_.set_tracker(&tracker_);

		if(isData_) {
			if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
			else object_selector_.lepton_type(1);
		}
		sync_ = values.count("sync");
		optim_ = values.count("optimization");
		has_pdfs_ = !(
			boost::starts_with(sample, "QCD") || 
			boost::starts_with(sample, "singletbar_tW") || boost::starts_with(sample, "singlet_tW") ||
			boost::starts_with(sample, "WZ") || boost::starts_with(sample, "WW") || boost::starts_with(sample, "ZZ"));
		bool isTTShift = boost::starts_with(sample, "ttJets_");
		if(sync_) tracker_.verbose(true);
    if(!isData_) mc_weights_.init(sample, has_pdfs_);

		runsys_ = !(values.count("nosys") || isData_ || isTTShift || sync_);
		if(!runsys_)
			systematics_ = {systematics::SysShifts::NOSYS};
		else {
			systematics_ = {
				Sys::NOSYS, 
				Sys::JES_UP,  Sys::JES_DW, 
				Sys::JER_UP,  Sys::JER_DW, 
				Sys::MET_UP,  Sys::MET_DW,
				Sys::PU_UP,   Sys::PU_DW,
				Sys::BEFF_UP, Sys::BEFF_DW, 
				Sys::BFAKE_UP, Sys::BFAKE_DW, 
				Sys::LEPEFF_UP, Sys::LEPEFF_DW,
			};
			if(isTTbar_) {
				systematics_.push_back(Sys::HDAMP_UP);
				systematics_.push_back(Sys::HDAMP_DW);
				systematics_.push_back(Sys::RENORM_UP);	
				systematics_.push_back(Sys::RENORM_DW); 
				systematics_.push_back(Sys::FACTOR_UP);	
				systematics_.push_back(Sys::FACTOR_DW);
				systematics_.push_back(Sys::RENFACTOR_UP);
				systematics_.push_back(Sys::RENFACTOR_DW);
			}
		}
	}
  
	TDirectory* getDir(std::string path){
		TDirectory* out = (TDirectory*) outFile_.Get(path.c_str());
		if(out) return out;
		outFile_.mkdir(path.c_str());
		return (TDirectory*) outFile_.Get(path.c_str());
	}

	template <class H, typename ... Args>
	void book(std::string folder, std::string name, Args ... args)
	{
		getDir(folder)->cd();
		histos_[folder][name] = RObject::book<H>(name.c_str(), args ...);
	}

	void book_combo_plots(string folder){
		if(optim_) return;
		book<TH1F>(folder, "mass_discriminant", "", 20,   0., 20.);
		book<TH1F>(folder, "nu_discriminant", "", 20,   0., 20.);
		book<TH1F>(folder, "full_discriminant", "", 40,   0., 40.);
    
		book<TH1F>(folder, "tmasshad", "", 100, 0., 500);			
		book<TH1F>(folder, "Wmasshad", "", 100, 0., 500);			
		book<TH1F>(folder, "lbratio" , "", 100, 0., 10.);			
		book<TH1F>(folder, "j2bratio", "", 100, 0., 10.);			
		book<TH1F>(folder, "bjets_pt", "", 500, 0., 500.);			
		book<TH1F>(folder, "wjets_pt", "", 500, 0., 500.);			
	}

	void fill_combo_plots(string folder, const Permutation &hyp){
		if(optim_) return;
		auto dir = histos_.find(folder);
		dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
		dir->second["nu_discriminant"  ].fill(hyp.NuDiscr(), evt_weight_);
		dir->second["full_discriminant"].fill(hyp.Prob(), evt_weight_);
		
		double whad_mass = hyp.WHad().M();
		double thad_mass = hyp.THad().M();
		dir->second["Wmasshad"].fill(whad_mass, evt_weight_);
		dir->second["tmasshad"].fill(thad_mass, evt_weight_);
		dir->second["lbratio" ].fill(hyp.L()->Pt()/hyp.BLep()->Pt(), evt_weight_);
		double minpt = ((hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJb() : hyp.WJa())->Pt();
		dir->second["j2bratio"].fill(minpt/hyp.BHad()->Pt(), evt_weight_);
		dir->second["bjets_pt"].fill(hyp.BHad()->Pt(), evt_weight_);
		dir->second["bjets_pt"].fill(hyp.BLep()->Pt(), evt_weight_);
		dir->second["wjets_pt"].fill(hyp.WJa()->Pt(), evt_weight_);		
		dir->second["wjets_pt"].fill(hyp.WJb()->Pt(), evt_weight_);		
	}

  void book_presel_plots(string folder) {
		if(optim_) return;
    book<TH1F>(folder, "nvtx_noweight", "", 100, 0., 100.);
    book<TH1F>(folder, "nvtx", "", 100, 0., 100.);
    book<TH1F>(folder, "rho", "", 100, 0., 100.);
    book<TH1F>(folder, "weight", "", 100, 0., 5.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(l) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "lep_eta"  , ";#eta(l) (GeV)", 300, -3, 3);
		book<TH1F>(folder, "jets_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "jets_eta" , ";#eta(j) (GeV)",  300, -3, 3);
		book<TH1F>(folder, "lead_jet_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "lead_jet_eta" , ";#eta(j) (GeV)",  300, -3, 3);
		book<TH1F>(folder, "jets_CSV" , ";#eta(j) (GeV)",  200, -1, 1);
		book<TH1F>(folder, "max_jets_CSV" , ";#eta(j) (GeV)",  200, -1, 1);
		book<TH1F>(folder, "lep_iso", ";#eta(j) (GeV)",  100, 0, 10);
		book<TH1F>(folder, "lep_wp" , ";#eta(j) (GeV)",  4, 0, 4);
		book<TH1F>(folder, "MT" , ";#eta(j) (GeV)",  500, 0, 500);
		book<TH1F>(folder, "MET" , ";#eta(j) (GeV)",  500, 0, 1000);
		book<TH1F>(folder, "METPhi" , ";#eta(j) (GeV)",  314, -1*Pi(), Pi());

		book<TH1F>(folder, "cMVA" , "",  100, -1, 1.1);
		book<TH1F>(folder, "cMVA_p11" , "",  100, -1, 1.1);
		
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);

		book<TH2F>(folder, "MT_iso" , ";#eta(j) (GeV)"  ,  10, 0, 100, 10, 0, 1);
		book<TH2F>(folder, "MT_btag" , ";#eta(j) (GeV)" ,  10, 0, 100, 10, 0, 1);
		book<TH2F>(folder, "iso_btag" , ";#eta(j) (GeV)",  10, 0, 1,   10, 0, 1);
		book<TH2D>(folder, "jets_cMVA_WP", "", 4, 0., 4., 4, 0., 4.);
  }

	int btag_idval(const IDJet* jet) {
		if(jet->BTagId(IDJet::BTag::MVATIGHT)) return 3;
		else if(jet->BTagId(IDJet::BTag::MVAMEDIUM)) return 2;
		else if(jet->BTagId(IDJet::BTag::MVALOOSE) ) return 1;
		return 0;
	}

  void fill_presel_plots(string folder, URStreamer &event) {
		if(optim_) return;
    auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "could not find: " << folder << endl;
			throw 42;
		}
		dir->second["nvtx"].fill(event.vertexs().size(), evt_weight_);
		dir->second["nvtx_noweight"].fill(event.vertexs().size());
		dir->second["rho"].fill(event.rho().value(), evt_weight_);
        dir->second["weight"].fill(evt_weight_);
		dir->second["lep_pt"].fill(object_selector_.lepton()->Pt(), evt_weight_);
		dir->second["lep_eta"].fill(object_selector_.lepton()->Eta(), evt_weight_);
		dir->second["njets" ].fill(object_selector_.clean_jets().size(), evt_weight_);
		double max_csv = -10000;
		double max_pt = -1;
		double max_eta = 0;
    for(IDJet* jet : object_selector_.clean_jets()) {
            dir->second["jets_pt"].fill(jet->Pt(), evt_weight_);
            dir->second["jets_eta"].fill(jet->Eta(), evt_weight_);
			dir->second["jets_CSV"].fill(jet->csvIncl(), evt_weight_);
			if(jet->csvIncl() > max_csv) max_csv = jet->csvIncl();
			dir->second["cMVA"    ].fill(jet->CombinedMVA(), evt_weight_);
			dir->second["cMVA_p11"].fill(pow(jet->CombinedMVA(), 11), evt_weight_);
			if(jet->Pt() > max_pt) {
				max_pt = jet->Pt();
				max_eta = jet->Eta();
			}

//        float pt_res = pt_jet_res_.getResolution({{JME::Binning::JetPt, jet->Pt()}, {JME::Binning::JetEta, jet->Eta()}, {JME::Binning::Rho, event.rho().value()}});
//
//        float sf = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, jet->Eta()}});
//        float sf_down = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, jet->Eta()}}, Variation::DOWN);
//        float sf_up = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, jet->Eta()}}, Variation::UP);

    }
		dir->second["lead_jet_pt" ].fill(max_pt , evt_weight_);
		dir->second["lead_jet_eta"].fill(max_eta, evt_weight_);

		dir->second["MET"   ].fill(object_selector_.met()->Et() , evt_weight_);
		dir->second["METPhi"].fill(object_selector_.met()->Phi(), evt_weight_);
		double mt = MT(object_selector_.lepton(), object_selector_.met());
		dir->second["MT"].fill(mt, evt_weight_);
		dir->second["max_jets_CSV"].fill(max_csv, evt_weight_);

		auto &clean_jets = object_selector_.clean_jets();
		sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->CombinedMVA() > B->CombinedMVA());});		
		dir->second["jets_cMVA_WP"].fill(btag_idval(clean_jets[0]), btag_idval(clean_jets[1]), evt_weight_);

		double iso = -1.;
		if(object_selector_.lepton_type() == 1) {
			iso = object_selector_.muon()->RelPFIsoDb();
		}
		else {
			auto electron = object_selector_.electron();			
			iso = electron->PFIsolationRho2015();
			int elwp = -1;
			if(electron->TightID25ns()) elwp = 3;
			else if(electron->MediumID25ns()) elwp =2;
			else if(electron->LooseID25ns()) elwp =1;
			else if(electron->VetoID25ns()) elwp =0;
			dir->second["lep_wp" ].fill(elwp, evt_weight_);
		}
		dir->second["lep_iso"].fill(iso, evt_weight_);
		dir->second["MT_iso"  ].fill(mt, iso, evt_weight_);
		dir->second["MT_btag" ].fill(mt, max_csv, evt_weight_);
		dir->second["iso_btag"].fill(iso, max_csv, evt_weight_);
  }

	void book_selection_plots(string folder, bool pdf) {
		book_presel_plots(folder);
		book_combo_plots(folder);
		const vector<double> mbinning = {
			250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 
			480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 610.0, 640.0, 
			680.0, 730.0, 800.0, 920.0, 1200.0};

		//book<TH1F>(folder, "m_tt", "", 180, 200., 2000);			
		book<TH1F>(folder, "m_tt", "", 19, &mbinning[0]);			

		book<TH1F>(folder, "pt_thad" , "", 200, 0., 2000);			
		book<TH1F>(folder, "eta_thad", "", 200, -10., 10);			
		book<TH1F>(folder, "pt_tlep" , "", 200, 0., 2000);			
		book<TH1F>(folder, "eta_tlep", "", 200, -10., 10);			
		book<TH1F>(folder, "pt_tt" , "", 200, 0., 2000);			
		book<TH1F>(folder, "eta_tt", "", 200, -10., 10);			

		book<TH1F>(folder, "full_discriminant_4j", "" , 640, 0., 80.);
		book<TH1F>(folder, "full_discriminant_5j", "" , 640, 0., 80.);
		book<TH1F>(folder, "full_discriminant_6Pj", "", 640, 0., 80.);

		//Angular variables
		//top angles
    book<TH1F>(folder, "tlep_ctstar", "", 200, -1.0001, 1.0001);
    book<TH1F>(folder, "thad_ctstar", "", 200, -1.0001, 1.0001);

		//delta
    book<TH1F>(folder, "cdelta_ld", "", 200, -1.0001, 1.0001);

		//helicity frame
    book<TH1F>(folder, "hframe_ctheta_d", "", 200, -1.0001, 1.0001);

		//2D plots
		book<TH2F>(folder, "mtt_tlep_ctstar", "", 19, &mbinning[0], 20, 0., 1.0001);
		//book<TH2F>(folder, "mtt_tlep_ctstar", "", 19, &mbinning[0], 5, 0., 1.0001);

		//PDF uncertainties
    //book<TH1D>(folder, "PDFYields", "", 249, 0, 249);		
		if(pdf) {
			auto dir = histos_.find(folder);
			if(dir == histos_.end()) {
				throw std::runtime_error("Could not access folder");
			}
			pdf_uncs_.book_replicas(folder, "m_tt", dir->second["m_tt"]);
			pdf_uncs_.book_replicas(folder, "mtt_tlep_ctstar", dir->second["mtt_tlep_ctstar"]);
		}
	}

	void fill_selection_plots(string folder, URStreamer &event, Permutation &hyp, bool pdf) {
		fill_presel_plots(folder, event);
		fill_combo_plots(folder, hyp);
		auto dir = histos_.find(folder);

		//PDF uncertainties
		//if(has_pdfs_) {
		//	const vector<Mcweight>& ws =  event.MCWeights();
		//	for(size_t h = 0 ; h < ws.size() ; ++h) {
		//		dir->second["PDFYields"].fill(h, evt_weight_*ws[h].weights()/ws[0].weights());
		//	}
		//}

		dir->second["pt_thad" ].fill(hyp.THad().Pt() , evt_weight_);
		dir->second["eta_thad"].fill(hyp.THad().Eta(), evt_weight_);
		dir->second["pt_tlep" ].fill(hyp.TLep().Pt() , evt_weight_);
		dir->second["eta_tlep"].fill(hyp.TLep().Eta(), evt_weight_);
		dir->second["pt_tt"   ].fill(hyp.LVect().Pt() , evt_weight_);
		dir->second["eta_tt"  ].fill(hyp.LVect().Eta(), evt_weight_);
		dir->second["m_tt"].fill(hyp.LVect().M(), evt_weight_);
		if(pdf) pdf_uncs_.fill_replicas(folder, "m_tt", hyp.LVect().M(), evt_weight_, event);
		
		if(object_selector_.clean_jets().size() == 4)
			dir->second["full_discriminant_4j" ].fill(hyp.Prob(), evt_weight_);
		else if(object_selector_.clean_jets().size() == 5)
			dir->second["full_discriminant_5j" ].fill(hyp.Prob(), evt_weight_);
		else
			dir->second["full_discriminant_6Pj"].fill(hyp.Prob(), evt_weight_);

		//choose up/down
		TLorentzVector whad(hyp.WHad());
		TLorentzVector j1(*hyp.WJa());
		TLorentzVector j2(*hyp.WJb());
		j1.Boost(whad.BoostVector()*-1);
		j2.Boost(whad.BoostVector()*-1);
		if(j2.E() > j1.E()) hyp.SwapWJets();

		//Angular variables
		hyp::TTbar ttang(hyp);
		auto ttcm = ttang.to_CM(); 
		auto tlepcm = ttang.tlep().to_CM();
		auto thadcm = ttang.thad().to_CM();
		
		//top angles
		double tlep_ctstar = min(fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())), 0.99999);
        double thad_ctstar = min(fabs(ttang.unit3D().Dot(ttcm.thad().unit3D())), 0.99999);
        dir->second["tlep_ctstar"].fill(tlep_ctstar, evt_weight_);
        dir->second["thad_ctstar"].fill(thad_ctstar, evt_weight_);
        //dir->second["tlep_ctstar"].fill(ttang.unit3D().Dot(ttcm.tlep().unit3D()), evt_weight_);
        //dir->second["thad_ctstar"].fill(ttang.unit3D().Dot(ttcm.thad().unit3D()), evt_weight_);
		dir->second["mtt_tlep_ctstar"].fill(
			hyp.LVect().M(), tlep_ctstar_evt_weight_ 
			//fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())), 
			//evt_weight_
			);
		if(pdf) pdf_uncs_.fill_replicas2D(
			folder, "mtt_tlep_ctstar", 
			hyp.LVect().M(), tlep_ctstar,
			//hyp.LVect().M(), fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())),
			evt_weight_, event
			);

    dir->second["cdelta_ld"].fill(tlepcm.W().l().unit3D().Dot(thadcm.W().down().unit3D()), evt_weight_);

		//helicity frame
    dir->second["hframe_ctheta_d"].fill(
			thadcm.W().down().unit3D().Dot(ttcm.thad().unit3D()), 
			evt_weight_);
	}

	void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
	}

  virtual void begin()
  {
    Logger::log().debug() << "htt_simple::begin" << endl;
    outFile_.cd();
		if(sync_) {
			sync_tree_ = new TTree("sync", "sync");
			sync_tree_->Branch("Run", &sync_info_.Run, "Run/i");
			sync_tree_->Branch("LumiSection", &sync_info_.LumiSection, "LumiSection/l");
			sync_tree_->Branch("Event", &sync_info_.Event, "Event/l");
			sync_tree_->Branch("RecoSuccess", &sync_info_.RecoSuccess, "RecoSuccess/O");
			sync_tree_->Branch("MassTT", &sync_info_.MassTT, "MassTT/f");
			sync_tree_->Branch("hasMuon", &sync_info_.hasMuon, "hasMuon/O");
			return;
		}
		

		vector<string> leptons = {"electrons", "muons"};		
		vector<string> subs    ;
		if(isSignal_) subs = {"/positive", "/negative"};
		else subs = {"/right", "/matchable", "/unmatchable", "/noslep"};		
		vector<string> lepIDs  = {"looseNOTTight", "tight"};		
		vector<string> MTs     = {"MTHigh"};		
		//vector<string> tagging = {"ctagged", "notag"};
		for(auto& lepton : leptons) {
			for(auto& sys : systematics_) {
				string sys_name = systematics::shift_to_name.at(sys);
				string dname = lepton+"/"+sys_name+"/preselection";
				//Logger::log().debug() << "Booking histos in: " << dname << endl;
				book_presel_plots(dname);				
				book_combo_plots(dname);
				for(auto& subsample : subs) {
					string sub = (isTTbar_ || isSignal_) ? subsample : "";
					for(auto& lepid : lepIDs) {
						for(auto& mt : MTs) {
							if(sys != Sys::NOSYS && (lepid != "tight" || mt != "MTHigh")) {
								continue;
							}
							stringstream dstream;
							dstream << lepton << sub << "/";
							dstream << sys_name << "/" << lepid << "/" << mt;
							
							//Logger::log().debug() << "Booking histos in: " << dstream.str() << endl;
							bool runpdf = (runsys_ && isTTbar_ && sys == systematics::SysShifts::NOSYS && lepid == "tight" && mt == "MTHigh");
							book_selection_plots(dstream.str(), runpdf);
						}
					}
					if(!(isTTbar_ || isSignal_)) break;
				}
			}
		}
	}

	void process_evt(systematics::SysShifts shift, URStreamer &event){
		tracker_.track("start", "electrons");
		tracker_.track("start", "muons");
		//float weight = 1.;

    //select reco objects
    if( !object_selector_.select(event, shift, sync_) ) return;
		int njets = object_selector_.clean_jets().size();		
		string leptype = (object_selector_.lepton_type() == -1) ? "electrons" : "muons";
		bool lep_is_tight = (object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTMU || 
												 object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTEL);

		tracker_.track("before tight lep", leptype);
		if(lep_is_tight) {
			tracker_.group(leptype);
			tracker_.track("object selection", leptype);
		}
		
		// cout << "Muon: " << *object_selector_.lepton() << endl;
    //MC Weight for lepton selection
		float lep_weight=1;
    if(!isData_) { 
      if(object_selector_.tight_muons().size() == 1) {
        lep_weight = muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta());
				// if(-1.3 < object_selector_.muon()->Eta() && object_selector_.muon()->Eta() < -1)
				// 	cout << "Mu " << *object_selector_.muon() << " weight: " << lep_weight << " prev: " << evt_weight_ << endl;
			}
      if(object_selector_.tight_electrons().size() == 1)
        lep_weight = electron_sf_.get_sf(object_selector_.electron()->Pt(), object_selector_.electron()->etaSC());
		}
		evt_weight_ *= lep_weight;
    string sys_name = systematics::shift_to_name.at(shift);
    stringstream presel_dir;
		presel_dir << leptype << "/";
		presel_dir << sys_name << "/preselection";
    if(!sync_ && lep_is_tight){
        fill_presel_plots(presel_dir.str(), event);
        tracker_.track("not sync and tight lep");
    }

		//cut on btag
		auto &clean_jets = object_selector_.clean_jets();
		sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->CombinedMVA() > B->CombinedMVA());});
		if(!clean_jets[0]->BTagId(cut_tight_b_)){
            tracker_.track("first b not pass");
            return;
        }
		if(!clean_jets[1]->BTagId(cut_loose_b_)){
            tracker_.track("second b not pass");
            return;
        }
		if(lep_is_tight) tracker_.track("b cuts", leptype);

    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met(), object_selector_.lepton_charge(),
          lep_is_tight) ) return;
          //event.rho().value(), lep_is_tight) ) return;
    if(lep_is_tight) tracker_.track("perm preselection");
				
    //find mc weight for btag
		double bweight = 1;
    if(!isData_ && !sync_){
        bweight *= btag_sf_.scale_factor(clean_jets, shift);
        tracker_.track("not data and not sync");
    }
		evt_weight_ *= bweight;
    
		// cout << object_selector_.clean_jets().size() << " --> " << permutator_.capped_jets().size() << endl;
		// for(auto i : object_selector_.clean_jets()) {
		// 	cout << *i << endl;
		// }
		// cout << endl << endl;
    //Find best permutation
    Permutation best_permutation;
		int idx = 0;
		for(auto& test_perm : permutator_.permutations()) {
			solver_.Solve(test_perm);

			// cout << "Permutation #" << idx << "  " << test_perm.Prob() << " -- " << best_permutation.Prob() << endl;
			// idx++;
			// if(test_perm.IsComplete()) {
			// 	cout << test_perm << endl;
			// }
			// else {
			// 	cout << "Permutation incomplete!" << endl
			// 			 << "l: " << test_perm.L() << ", b_l: " << test_perm.BLep() << ", b_h: " 
			// 			 << test_perm.BHad() << ", j1: " << test_perm.WJa() << ", j2: " 
			// 			 << test_perm.WJb() << endl;
			// }
			if(!sync_ && lep_is_tight) fill_combo_plots(presel_dir.str(), test_perm);
			if(test_perm.Prob()  < best_permutation.Prob()) {
				best_permutation = test_perm;
				// cout << "This is the best so far!" << endl; 
			}
			// cout << endl;
    }

//        float wja_sf = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}});
//        float wja_sf_down = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}}, Variation::DOWN);
//        float wja_sf_up = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}}, Variation::UP);

		bool reco_success = (best_permutation.IsComplete() && best_permutation.Prob() <= 1E9);
    if(!sync_ && !reco_success){
        tracker_.track("not sync and not reco success");
        return;
    }
    if(lep_is_tight) tracker_.track("best perm");

//		cout << "BEST PERM" << endl;
//		if(best_permutation.IsComplete()){
//			cout << best_permutation << endl;
//		} 
//		else cout << "Permutation incomplete!" <<endl;
//        cout << "wja pt res: " << wja_pt_res << ", sf: " << wja_sf << ", sf down: " << wja_sf_down << ", sf up: " << wja_sf_up << endl;


		if(sync_) {
			if(lep_is_tight) {
				sync_info_.Run  = event.run;
				sync_info_.LumiSection = event.lumi;
				sync_info_.Event= event.evt;
				sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
				sync_info_.RecoSuccess = reco_success;
				if(reco_success)
					sync_info_.MassTT = best_permutation.LVect().M();
				else
					sync_info_.MassTT = -1;
				sync_tree_->Fill();
				tracker_.track("SYNC END", leptype);
			}
			return;
		}

		//create event dir (contains the dir path to be filled)
		stringstream evtdir;		
		evtdir << leptype;

    //Gen matching (TT events only)
    Permutation matched_perm;
    if(isTTbar_) {
			if(genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) {
				matched_perm = matcher_.match(
					genp_selector_.ttbar_final_system(),
					object_selector_.clean_jets(), 
					object_selector_.veto_electrons(),
					object_selector_.veto_muons()
					);
				matched_perm.SetMET(object_selector_.met());
				
				if(best_permutation.IsCorrect(matched_perm)) evtdir << "/right";
				else if(matched_perm.IsComplete()) evtdir << "/matchable";
				else evtdir << "/unmatchable";
			}
			else evtdir << "/noslep";
    } 
		else if (isSignal_) {
			if(evt_weight_ > 0.) evtdir << "/positive";
			else evtdir << "/negative";
		}
    if(lep_is_tight) tracker_.track("matched perm");

		//check isolation type
		bool runpdf = (runsys_ && isTTbar_ && shift == systematics::SysShifts::NOSYS);
		evtdir << "/" << sys_name << "/";
		bool tight;
		switch(object_selector_.event_type()) {
		case TTObjectSelector::EvtType::TIGHTMU: 
		case TTObjectSelector::EvtType::TIGHTEL: evtdir << "tight"; runpdf &= true; tight=true; break;
		case TTObjectSelector::EvtType::LOOSEMU: 
		case TTObjectSelector::EvtType::LOOSEEL: evtdir << "looseNOTTight"; runpdf &= false; tight=false; break;
		default: throw 42; break;
		}
		if(!tight && shift != Sys::NOSYS) return;

		//MT category (fixed, now)
		evtdir << "/MTHigh";

		//fill right category
		fill_selection_plots(evtdir.str(), event, best_permutation, runpdf);
		if(lep_is_tight) tracker_.track("END", leptype);
	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int report = values["report"].as<int>();
		int skip  = values["skip"].as<int>();
		string pick = values["pick"].as<string>();
		if(values.count("bystep")) tracker_.verbose(true);
		EventList picker;
		if(pick.size() != 0) {
			EventList nn(pick);
			picker = nn;
		}

    if(evt_idx_ >= limit) return;
    Logger::log().debug() << "htt_simple::analyze" << endl;
		Logger::log().debug() << tree_ << " " << tree_->GetEntries() << endl;
    URStreamer event(tree_);


    while(event.next()) {
			// if(evt_idx_ % 1000 == 0) Logger::log().warning() << "Beginning event: " <<
      //                           evt_idx_ << endl;

			if(picker.active()) {
				if(picker.contains(event.run, event.lumi, event.evt)) {
					Logger::log().fatal() << "Picking event " << " run: " << event.run << " lumisection: " << 
						event.lumi << " eventnumber: " << event.evt << endl;
				}
				else continue;
			}
			if(limit > 0 && evt_idx_ > limit) {
				return;
			}
			evt_idx_++;
			if(skip > 0 && evt_idx_ < skip) {
				continue;
			}
			evt_idx_--;
			if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << " run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
			evt_idx_++;

			if(isTTbar_){
				genp_selector_.select(event);			
			}

			for(auto shift : systematics_){
        evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);

				if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(shift, event);
				tracker_.group("");
				if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
			}
            //cout << "evt_idx: " << evt_idx_ << endl;
		} //while(event.next())
   }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  virtual void end(){
		outFile_.cd();
		tracker_.writeTo(outFile_);
		if(sync_) {
			sync_tree_->Write();
			return;
		}
		outFile_.Write();
		for(string ltype : {"electrons", "muons"}) {
			TDirectory* tdir = (TDirectory*) outFile_.Get(ltype.c_str());
			tracker_.writeTo(*tdir, ltype);
		}
	}

  //do you need command-line or cfg options? If so implement this 
  //method to book the options you need. CLI parsing is provided
  //by AnalysisFW/interface/URParser.h and uses boost::program_options
  //look here for a quickstart tutorial: 
  //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html
  static void setOptions() {
		URParser &parser = URParser::instance();
		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("nopdfs", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report,r", opts::value<int>()->default_value(10000), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("optimization", "run in optimization mode: reduced number of histograms")
      ("sync", "dump sync ntuple")
      ("bystep", "print every step")
      ("nosys", "do not run systematics")
      ("noweights", "do not run systematics")
      ("pick", opts::value<string>()->default_value(""), "pick from evtlist");

    parser.addCfgParameter<std::string>("general", "csv_sffile", "");
    parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");

	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<htt_simple> test;
	int excode = test.run();
	//Logger::log().debug() << "RUNNING DONE " << std::endl;
	auto files = gROOT->GetListOfFiles(); //make ROOT aware that some files do not exist, because REASONS
	Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl; //need to print out this otherwise ROOT loses its shit in 7.4.X (such I/O, much features)
  return excode;
}
