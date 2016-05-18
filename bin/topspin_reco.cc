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

using namespace TMath;

class topspin_reco : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
private:
	//histograms and helpers
	map<string, map<string, RObject> > histos_;
	map<TTNaming, string> naming_;
	CutFlowTracker tracker_;

	//switches
	bool isTTbar_, isSignal_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTPermutator permutator_;
  TTGenMatcher matcher_;
  TTBarSolver solver_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  //BTagSFProducer btag_sf_;

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

  unsigned long evt_idx_ = 0;

public:
  topspin_reco(const std::string output_filename):
    AnalyzerBase("topspin_reco", output_filename), 
		tracker_(),
		working_points_(),
    object_selector_(),
    permutator_(),
    matcher_(),
    solver_(),
    mc_weights_(),
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_()    //btag_sf_("permutations.tightb", "permutations.looseb")
	{
    //cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << 
    //  " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;
    
    URParser &parser = URParser::instance();

		TUUID id;  
		randomizer_.SetSeed(id.Hash());   

    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

    //find out which sample are we running on
    opts::variables_map &values = parser.values();
		string output_file = values["output"].as<std::string>();
		string sample = systematics::get_sample(output_file);
    isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
		isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal_;
    if(!isTTbar_) {
      Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
      throw 49;
    }

    if(isSignal_) genp_selector_.setmode(TTGenParticleSelector::SelMode::MADGRAPHLHE);
    else genp_selector_.setmode(TTGenParticleSelector::SelMode::LHE);

    mc_weights_.init(sample);

    //Init solver
    string filename = "prob_ttJets.root";
    Logger::log().debug() << "solver file: " << filename << endl;
    TFile probfile(DataFile(filename).path().c_str());
    TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
    solver_.Init(td, false, true, true);

		naming_[TTNaming::RIGHT ] = "semilep_visible_right";
		naming_[TTNaming::RIGHT_THAD ] = "semilep_right_thad";
		naming_[TTNaming::RIGHT_WHAD ] = "semilep_right_whad";
		naming_[TTNaming::RIGHT_TLEP ] = "semilep_right_tlep";
		naming_[TTNaming::WRONG ] = "semilep_wrong" 			 ;
		naming_[TTNaming::OTHER ] = 	"other"              ;

    Logger::log().debug() << "Constructor completed" << std::endl;
	}
  
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

  virtual void begin()
  {
    Logger::log().debug() << "topspin_reco::begin" << endl;
    outFile_.cd();
		vector<string> folders;
		if(isTTbar_) {
      for(auto &entry : naming_) folders.push_back(entry.second);
    }
		else folders = {""};
  }

	void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
	}

	TTNaming get_ttdir_name(Permutation &gen, Permutation &reco) {
    if(reco.IsCorrect(gen)) return TTNaming::RIGHT;
    else if(reco.IsTHadCorrect(gen)) return TTNaming::RIGHT_THAD; 
    else if(gen.IsWHadComplete() && reco.IsWHadCorrect(gen)) return TTNaming::RIGHT_WHAD;
    else if(reco.IsBLepCorrect(gen)) return TTNaming::RIGHT_TLEP;
    else if(genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) return TTNaming::WRONG;
    else return TTNaming::OTHER;
	}
		

	void process_evt(systematics::SysShifts shift, URStreamer &event){
		tracker_.track("start");
		//float weight = 1.;

    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("object selection");

    string sys_name = systematics::shift_to_name.at(shift);
    string presel_dir = sys_name;
    if(isTTbar_){
      presel_dir = naming_.at(TTNaming::RIGHT) + "/" + sys_name;
    }

    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;
    tracker_.track("perm preselection");

				
    //Find best permutation
    bool go_on = true;
    Permutation best_permutation;
    size_t ncycles_=0;
    while(go_on) {
      ncycles_++;
      Permutation test_perm = permutator_.next(go_on);
      if(go_on) {
        test_perm.Solve(solver_);
        double bjet_lpt = Max(test_perm.BHad()->Pt(), test_perm.BLep()->Pt());
        fill_combo_plots(presel_dir+"/permutations", test_perm);
        //if(bjet_lpt < test_perm.WJa()->Pt()) continue;
        if(ordering_fcn_(test_perm, best_permutation)){
          best_permutation = test_perm;
        }
      }
    }
    if(!best_permutation.IsComplete() || best_permutation.Prob() > 1E9) return; //FIXME, is right??? best_permutation.Prob() > 1E9
    tracker_.track("best perm");
    
    size_t capped_jets_size = permutator_.capped_jets().size();
    best_permutation.permutating_jets(capped_jets_size);

    //Gen matching
    Permutation matched_perm;
    // if(isTTbar_) Logger::log().debug() << " --  GEN  --" << genp_selector_.ttbar_final_system() << std::endl;
    if(isTTbar_ && genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) {
      matched_perm = matcher_.match(
        genp_selector_.ttbar_final_system(),
        object_selector_.clean_jets(), 
        object_selector_.loose_electrons(),
        object_selector_.loose_muons()
        );
      matched_perm.SetMET(object_selector_.met());
    }
    tracker_.track("matched perm");

		string sys_dir = sys_name;
		//Logger::log().debug() << "\n\nShift: " << shift << " name: " << root_dir << endl;
    double weight = evt_weight_;
    string ttsubdir = "";

    TTNaming dir_id = get_ttdir_name(matched_perm, best_permutation);
    ttsubdir = naming_.at(dir_id) + "/";

	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();

    if(evt_idx_ >= limit) return;
    Logger::log().debug() << "topspin_reco::analyze" << endl;
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
			if(evt_idx_ % 10000 == 1) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

			//long and time consuming
      genp_selector_.select(event);			

      evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);
      process_evt(shift, event);
		} //while(event.next())
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
		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("nosys", opts::value<int>()->default_value(0), "do not run systematics")
      ("nopdfs", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");

    parser.addCfgParameter<std::string>("general", "csv_sffile", "");
    parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");

	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<topspin_reco> test;
  return test.run();
}
