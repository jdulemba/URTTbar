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
#include "TUUID.h"   
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "Analyses/URTTbar/interface/TTBarPlots.h"
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

using namespace TMath;

class htt_simple : public AnalyzerBase
{
public:
	//enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
private:
	//histograms and helpers
	unordered_map<string, unordered_map<string, RObject> > histos_;
	CutFlowTracker tracker_;

	//switches
	bool isTTbar_, isSignal_, isData_;

  //selectors and helpers
  TTObjectSelector object_selector_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  BTagSFProducer btag_sf_;

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

  unsigned long evt_idx_ = 0;
	vector<systematics::SysShifts> systematics_;

  //cuts
  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;
	double cut_MT_=0;

	//sync
	bool sync_;
	TTree *sync_tree_;
	map<string, float> sync_info_;
	long long sync_evt_; //the only one cannot work with floats

public:
  htt_simple(const std::string output_filename):
    AnalyzerBase("htt_simple", output_filename), 
		tracker_(),
    object_selector_(),
    mc_weights_(),
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_(),    
		btag_sf_("bjets.tightb", "bjets.looseb"),
		systematics_(),
		sync_(false),
		sync_tree_(0),
		sync_info_()
	{
    cut_tight_b_ = btag_sf_.tight_cut();
    cut_loose_b_ = btag_sf_.loose_cut();
    //cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << 
    //  " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;
    
    URParser &parser = URParser::instance();
    parser.addCfgParameter<float>("event", "MT", "");
		parser.parseArguments();
		cut_MT_ = parser.getCfgPar<float>("event", "MT");

		TUUID id;  
		randomizer_.SetSeed(id.Hash());   

    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);

    //find out which sample are we running on
    opts::variables_map &values = parser.values();
		string output_file = values["output"].as<std::string>();
		string sample = systematics::get_sample(output_file);
    isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
		isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal_;
		isData_  = boost::starts_with(sample, "data");

		sync_ = values.count("sync");
    if(!isData_) mc_weights_.init(sample);

		systematics_ = {systematics::SysShifts::NOSYS};//systematics::get_systematics(output_file);
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

  void book_presel_plots(string folder) {
    book<TH1F>(folder, "nvtx_noweight", "", 100, 0., 100.);
    book<TH1F>(folder, "nvtx", "", 100, 0., 100.);
    book<TH1F>(folder, "rho", "", 100, 0., 100.);
    book<TH1F>(folder, "weight", "", 100, 0., 20.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(l) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "lep_eta"  , ";#eta(l) (GeV)", 300, -3, 3);
		book<TH1F>(folder, "jets_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "jets_eta" , ";#eta(j) (GeV)",  300, -3, 3);
		book<TH1F>(folder, "jets_CSV" , ";#eta(j) (GeV)",  200, -1, 1);
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
  }

  void fill_presel_plots(string folder, URStreamer &event) {
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
    for(IDJet* jet : object_selector_.clean_jets()) {
      dir->second["jets_pt"].fill(jet->Pt(), evt_weight_);
      dir->second["jets_eta"].fill(jet->Eta(), evt_weight_);
			dir->second["jets_CSV"].fill(jet->csvIncl(), evt_weight_);
    }
  }

	void book_selection_plots(string folder) {
		book_presel_plots(folder);
	}

	void fill_selection_plots(string folder, URStreamer &event) {
		fill_presel_plots(folder, event);
	}

	void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
	}

  virtual void begin()
  {
    Logger::log().debug() << "htt_simple::begin" << endl;
    outFile_.cd();
		if(sync_) {
			sync_tree_ = new TTree("sync", "sync");
			sync_info_["hasMuon"] = -1;
			sync_info_["finalweight"] = -1;
			sync_info_["leptonweight"] = -1;
			sync_info_["btagweight"] = -1;
			sync_info_["puweight"] = -1;
			sync_info_["mcweight"] = -1;
			sync_info_["lumi"] = -1;

			sync_tree_->Branch("evt", &sync_evt_);
			for(auto& entry :  sync_info_) {
				sync_tree_->Branch(entry.first.c_str(), &entry.second);
			}
			return;
		}
		

		vector<string> leptons = {"electrons", "muons"};		
		vector<string> njets = {/*"3jets",*/ "4jets", "5Pjets"};
		//vector<string> tagging = {"ctagged", "notag"};
		for(auto& lepton : leptons) {
			for(auto& sys : systematics_) {
				string sys_name = systematics::shift_to_name.at(sys);
				string dname = lepton+"/"+sys_name+"/preselection";
				Logger::log().debug() << "Booking histos in: " << dname << endl;
				book_presel_plots(dname);
				for(auto& njet : njets) {
					stringstream dstream;
					dstream << lepton << "/";
					dstream << sys_name<<"/"<<njet;
					
					Logger::log().debug() << "Booking histos in: " << dstream.str() << endl;
					book_selection_plots(dstream.str());
				}
			}			
		}
	}

	void process_evt(systematics::SysShifts shift, URStreamer &event){
		tracker_.track("start", "electrons");
		tracker_.track("start", "muons");
		//float weight = 1.;

    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
		int njets = object_selector_.clean_jets().size();		
		bool totrack = true;//(njets == 3);
		if(!totrack) tracker_.deactivate();
		string leptype = (object_selector_.lepton_type() == -1) ? "electrons" : "muons";
		tracker_.group(leptype);
		tracker_.track("object selection", leptype);

    //MC Weight for lepton selection
		float lep_weight=1;
    if(!isData_) { 
      if(object_selector_.tight_muons().size() == 1)
        lep_weight = muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
      if(object_selector_.medium_electrons().size() == 1)
        lep_weight = electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
		}
		evt_weight_ *= lep_weight;

    string sys_name = systematics::shift_to_name.at(shift);
    stringstream presel_dir;
		presel_dir << leptype << "/";
		presel_dir << sys_name << "/preselection";
    if(!sync_) fill_presel_plots(presel_dir.str(), event);

		//cut on btag
		auto &clean_jets = object_selector_.clean_jets();
		sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->csvIncl() > B->csvIncl());});
		if(!clean_jets[0]->BTagId(cut_tight_b_)) return;
		tracker_.track("tight b cut", leptype);
		if(!clean_jets[1]->BTagId(cut_loose_b_)) return;
		tracker_.track("loose b cut", leptype);

    //find mc weight for btag
		double bweight = 1;
    if(!isData_) bweight *= btag_sf_.scale_factor(clean_jets, shift);
		evt_weight_ *= bweight;

		//cut on MT
		double mt = (*object_selector_.lepton()+*object_selector_.met()).Mt() ;
		if(mt < cut_MT_) return;
		tracker_.track("MT cut", leptype);
    
		stringstream evtdir;		
		evtdir << leptype << "/" << sys_name << "/";
		switch(njets) {
		case 3: evtdir << "3jets"; break;
		case 4: evtdir << "4jets"; break;
		default:evtdir << "5Pjets"; break;
		}
		if(!sync_){
			fill_selection_plots(evtdir.str(), event);
		} else {
			cout << "passing event " << " / " << event.run << ":" << event.lumi << ":" << event.evt << 
				" / Muon: "<< (object_selector_.lepton_type() == -1) << endl;
			sync_info_["hasMuon"] = (object_selector_.lepton_type() == -1) ? 0 : 1;
			sync_info_["btagweight"] = bweight;
			sync_info_["leptonweight"] = lep_weight;
			sync_info_["finalweight"] = evt_weight_;
			sync_tree_->Fill();
		} 
	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int report = values["report"].as<int>();
		int skip  = values["skip"].as<int>();

    if(evt_idx_ >= limit) return;
    Logger::log().debug() << "htt_simple::analyze" << endl;
    URStreamer event(tree_);

    while(event.next()) {
			// if(evt_idx_ % 1000 == 0) Logger::log().warning() << "Beginning event: " <<
      //                           evt_idx_ << endl;
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
			if(sync_) {
				sync_info_["lumi"] = event.lumi;
				sync_evt_ = event.evt;
				sync_info_["puweight"] = (isData_) ? 1. : mc_weights_.pu_weight( event);
				sync_info_["mcweight"] = (isData_) ? 1. : mc_weights_.gen_weight(event);
			}
			for(auto shift : systematics_){
        evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);

				if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(shift, event);
				tracker_.group("");
				if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
			}
		} //while(event.next())
   }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  virtual void end(){
		outFile_.cd();
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
      ("nosys", opts::value<int>()->default_value(0), "do not run systematics")
      ("nopdfs", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report,r", opts::value<int>()->default_value(10000), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("sync", "dump sync ntuple");

    parser.addCfgParameter<std::string>("general", "csv_sffile", "");
    parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");

	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<htt_simple> test;
  return test.run();
}
