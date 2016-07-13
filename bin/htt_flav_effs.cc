#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/TTBarResponse.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "Analyses/URTTbar/interface/TTBarPlots.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "TRandom3.h"
#include "Analyses/URTTbar/interface/helper.h"
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include <unordered_map>

using namespace std;

class htt_flav_effs : public AnalyzerBase
{
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  //histos
  unordered_map<string, unordered_map<string, RObject> > histos_;

  //selectors and helpers
  TTObjectSelector object_selector_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;

	vector<systematics::SysShifts> systematics_;

	//Scale factors
  LeptonSF electron_sf_, muon_sf_;
  string loose_blabel_, tight_blabel_;

  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE;
  IDJet::BTag cut_loose_b_=IDJet::BTag::NONE;
	float cut_MT_;
	bool sameCut_;

public:
  htt_flav_effs(const std::string output_filename):
    AnalyzerBase("htt_flav_effs", output_filename),
		tracker_(),
    loose_blabel_(),
    tight_blabel_(),
    object_selector_(),
    evt_weight_(1.),
    mc_weights_(),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"){
    
    //tracker_.verbose(true);    
    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);

    opts::variables_map &values = URParser::instance().values();
		URParser& parser = URParser::instance();

    //find out which sample are we running on
		string output_file = values["output"].as<std::string>();
    //DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
		string sample = systematics::get_sample(output_file);
		bool isTTbar = boost::starts_with(sample, "ttJets");

    //choose systematics to run based on sample
    systematics_ = {systematics::SysShifts::NOSYS};
    if(!isTTbar) {
      Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
      throw 49;
    }
    
    mc_weights_.init(sample);
    
    parser.addCfgParameter<float>("event", "MT", "");
		parser.addCfgParameter<string>("bjets", "tightb", "");
		parser.addCfgParameter<string>("bjets", "looseb", "");
		parser.parseArguments();
		cut_MT_ = parser.getCfgPar<float>("event", "MT");

    cut_tight_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("bjets", "tightb"));
    cut_loose_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("bjets", "looseb"));
		sameCut_=false;
		if(cut_tight_b_ == cut_loose_b_) sameCut_=true;

    //get cut names
    tight_blabel_ = IDJet::tag2string(cut_tight_b_);
    loose_blabel_ = IDJet::tag2string(cut_loose_b_);
  };
  
  ~htt_flav_effs() {
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

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
		vector<string> dirs = {"alljets", "mtpass", "mtfail"};
		vector<string> ltypes = {"electrons", "muons"};
		for(auto ltype : ltypes) {
			for(auto shift : systematics_) {
				string sys_name = systematics::shift_to_name.at(shift);
				for(auto dir : dirs) {
					string dirname = ltype+"/"+sys_name+"/"+dir;
					TDirectory* tdir = outFile_.mkdir(dirname.c_str());
					tdir->cd();
					Logger::log().debug() << "booking: " << dirname << std::endl;					
					for(string cut : {loose_blabel_, tight_blabel_}) {
						for(string flav : {"_bjet_", "_cjet_", "_ljet_"}) {
							for(string cat : {"all", "pass"}) {
								string name = "btag_" + cut + flav + cat;
								Logger::log().debug() << "booking: " << name << std::endl;
								book<TH2D>(dirname, name, "btag SF input histograms;p_{T};#eta", 200, 0, 1000, 60, -3, 3);
							}
						}
						if(sameCut_) break;
					}
				}//for(auto dir : dirs) {
			}
		}
  }

	void fill(string folder) {
    auto plots = histos_.find(folder);
		if(plots == histos_.end()) {
			Logger::log().error() << "could not find: " << folder << endl;
			throw 42;
		}
    for(auto jet : object_selector_.clean_jets()) {
      int jet_flav = Abs(jet->hadronFlavour());
      string flav;
      if(jet_flav == ura::PDGID::b) flav = "_bjet_";
      else if(jet_flav == ura::PDGID::c) flav = "_cjet_";
      else flav = "_ljet_";
      
      if(jet->BTagId(cut_loose_b_)) 
        plots->second["btag_"+loose_blabel_+flav+"pass"].fill(jet->Pt(), jet->Eta(), evt_weight_);
      plots->second["btag_"+loose_blabel_+flav+"all"].fill(jet->Pt(), jet->Eta(), evt_weight_);
			if(!sameCut_) {
				if(jet->BTagId(cut_tight_b_))         
					plots->second["btag_"+tight_blabel_+flav+"pass"].fill(jet->Pt(), jet->Eta(), evt_weight_);
				plots->second["btag_"+tight_blabel_+flav+"all"].fill(jet->Pt(), jet->Eta(), evt_weight_);
			}
    } //for(auto jet : leading_jets)
	}
  
  void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) { 
    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("obj selection");
		string leptype = (object_selector_.lepton_type() == -1) ? "electrons" : "muons";
		string sys_name = systematics::shift_to_name.at(shift);
    //find mc weight
    if(object_selector_.tight_muons().size() == 1)
      evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    if(object_selector_.medium_electrons().size() == 1)
      evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    tracker_.track("MC weights");

		fill(leptype+"/"+sys_name+"/"+"alljets");

		double mt = (*object_selector_.lepton()+*object_selector_.met()).Mt() ;
		if(mt < cut_MT_) {
			fill(leptype+"/"+sys_name+"/"+"mtfail");
		} else {
			fill(leptype+"/"+sys_name+"/"+"mtpass");
		}    
    tracker_.track("end");
  }

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
		unsigned long evt_idx = 0;
    URStreamer event(tree_);

    Logger::log().debug() << "--retrieving running conditions--" << endl;
    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();
    Logger::log().debug() << "--DONE--" << endl;

    while(event.next()) {
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			if(evt_idx % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx << endl;
      tracker_.track("start");

      tracker_.deactivate();    
			for(auto shift : systematics_){
        evt_weight_ = mc_weights_.evt_weight(event, shift);
				if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(event, shift);
				if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
			}

    }  //while(event.next())
  }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  virtual void end(){
		outFile_.Write();
		tracker_.writeTo(outFile_);
	}

  //do you need command-line or cfg options? If so implement this 
  //method to book the options you need. CLI parsing is provided
  //by AnalysisFW/interface/URParser.h and uses boost::program_options
  //look here for a quickstart tutorial: 
  //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html
  static void setOptions() {
    URParser &parser = URParser::instance();
    parser.addCfgParameter<int>("general", "skew_ttpt_distribution", "Should I skew the pt distribution? (0/1)");
    parser.addCfgParameter<int>("general", "pseudotop", "should I use pseudo-top? (0/1)");

		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");
    parser.addCfgParameter<std::string>("best_permutation", "tightb", "");
    parser.addCfgParameter<std::string>("best_permutation", "looseb", "");
  }
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<htt_flav_effs> test;
  return test.run();
}
