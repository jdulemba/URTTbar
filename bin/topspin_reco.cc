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

class topspin_reco : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
private:
	//histograms and helpers
	unordered_map<string, unordered_map<string, RObject> > histos_;
	map<TTNaming, string> naming_;
	CutFlowTracker tracker_;

	//switches
	bool isTTbar_, isSignal_, isData_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTPermutator permutator_;
  TTGenMatcher matcher_;
  TTBarSolver solver_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  BTagSFProducer btag_sf_;

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

  unsigned long evt_idx_ = 0;
	vector<systematics::SysShifts> systematics_;

public:
  topspin_reco(const std::string output_filename):
    AnalyzerBase("topspin_reco", output_filename), 
		tracker_(),
    object_selector_(),
    permutator_(),
    matcher_(),
    solver_(true),
    mc_weights_(),
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_(),    
		btag_sf_("permutations.tightb", "permutations.looseb"),
		systematics_()
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
		isData_  = boost::starts_with(sample, "data");

    if(isSignal_) genp_selector_.setmode(TTGenParticleSelector::SelMode::MADGRAPHLHE);
    else genp_selector_.setmode(TTGenParticleSelector::SelMode::LHE);

    if(!isData_) mc_weights_.init(sample);

		systematics_ = {systematics::SysShifts::NOSYS};//systematics::get_systematics(output_file);

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

	void book_combo_plots(string folder){
		book<TH1F>(folder, "mass_discriminant", "", 20,   0., 20.);
		book<TH1F>(folder, "comb_discriminant", "", 60,   0., 30.);
		book<TH1F>(folder, "nu_chi2", "",  50,  0., 10.);
		book<TH1F>(folder, "tmasshad", "", 100, 0., 500);			
		book<TH1F>(folder, "Wmasshad", "", 100, 0., 500);			
	}

	void fill_combo_plots(string folder, const Permutation &hyp){
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "could not find: " << folder << endl;
			throw 42;
		}
		dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
		dir->second["comb_discriminant"].fill(hyp.Prob(), evt_weight_);
		dir->second["nu_chi2"].fill(hyp.NuDiscr(), evt_weight_);
		dir->second["Wmasshad"].fill(hyp.WHad().M(), evt_weight_);
		dir->second["tmasshad"].fill(hyp.THad().M(), evt_weight_);
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
    }
  }

	void book_selection_plots(string folder) {
		book_combo_plots(folder);
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		book<TH1F>(folder, "lep_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "had_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "wjets_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(#ell) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "Whad_DR"  , ";m_{W}(had) (GeV)", 100, 0., 10.);
		book<TH1F>(folder, "Whad_pt"  , ";m_{W}(had) (GeV)", 100, 0., 500.);
	}

	void fill_selection_plots(string folder, Permutation &hyp) {
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "could not find: " << folder << endl;
			throw 42;
		}
		fill_combo_plots(folder, hyp);
		dir->second["njets"    ].fill(object_selector_.clean_jets().size(), evt_weight_);
		dir->second["lep_b_pt" ].fill(hyp.BLep()->Pt(), evt_weight_);
		dir->second["had_b_pt" ].fill(hyp.BHad()->Pt(), evt_weight_);
		dir->second["wjets_pt" ].fill(hyp.WJa()->Pt(), evt_weight_);
		dir->second["lep_pt"   ].fill(hyp.L()->Pt(), evt_weight_);
		if(hyp.WJb()) {
			dir->second["wjets_pt" ].fill(hyp.WJb()->Pt(), evt_weight_);
			dir->second["Whad_DR"  ].fill(hyp.WJa()->DeltaR(*hyp.WJb()), evt_weight_);
			dir->second["Whad_pt"  ].fill(hyp.WHad().Pt(), evt_weight_);
		}
	}

	void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
	}

  virtual void begin()
  {
    Logger::log().debug() << "topspin_reco::begin" << endl;
    outFile_.cd();
		vector<string> leptons = {"electrons", "muons"};
		vector<string> types;
		if(isTTbar_) {
      for(auto &entry : naming_) types.push_back(entry.second);
    }
		else types = {""};
		
		vector<string> njets = {"3jets", "4jets", "5Pjets"};
		//vector<string> tagging = {"ctagged", "notag"};
		for(auto& lepton : leptons) {
			for(auto& type : types) {
				string base = (type.size() == 0 ) ? lepton : lepton+"/"+type;
				for(auto& sys : systematics_) {
					string sys_name = systematics::shift_to_name.at(sys);
					if(type == "semilep_visible_right" || type.size() == 0) {
						string dname = base+"/"+sys_name+"/preselection";
						Logger::log().debug() << "Booking histos in: " << dname << endl;
						book_presel_plots(dname);
						book_combo_plots(dname);
					}
					for(auto& njet : njets) {
						stringstream dstream;
						dstream << lepton << "/";
						if(type.size() != 0) {
							dstream << type << "/";
						} 
						dstream << sys_name<<"/"<<njet;
						
						Logger::log().debug() << "Booking histos in: " << dstream.str() << endl;
						book_selection_plots(dstream.str());
					}
				}
			}
		}
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
    if(!isData_) { 
      if(object_selector_.tight_muons().size() == 1)
        evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
      if(object_selector_.tight_electrons().size() == 1)
        evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
		}

    string sys_name = systematics::shift_to_name.at(shift);
    stringstream presel_dir;
		presel_dir << leptype << "/";
    if(isTTbar_){
      presel_dir << naming_.at(TTNaming::RIGHT) << "/"; 
    }
		presel_dir << sys_name << "/preselection";
    fill_presel_plots(presel_dir.str(), event);
		
    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;
    //find mc weight for btag
    if(!isData_) evt_weight_ *= btag_sf_.scale_factor(permutator_.capped_jets(), shift);
    tracker_.track("perm preselection", leptype);

    //Find best permutation
    bool go_on = true;
    Permutation best_permutation;
    size_t ncycles=0;
		bool lazy_solving = (njets == 3);
		auto ordering = (njets > 3) ? [](const Permutation &one, const Permutation &two) {return one.Prob() < two.Prob();} : \
		    [](const Permutation &one, const Permutation &two) {return one.NuDiscr() < two.NuDiscr();};
		for(auto test_perm : permutator_.permutations()) {
			ncycles++;
			solver_.Solve(test_perm, lazy_solving);
			fill_combo_plots(presel_dir.str(), test_perm);
			if(ordering(test_perm, best_permutation)){
				best_permutation = test_perm;
			}
    }
    if((njets != 3 && best_permutation.Prob() > 1E9) || ncycles == 0) return; //FIXME, is right??? best_permutation.Prob() > 1E9
    tracker_.track("best perm", leptype);
    
    size_t capped_jets_size = permutator_.capped_jets().size();
    best_permutation.permutating_jets(capped_jets_size);

    //Gen matching
    Permutation matched_perm;
    // if(isTTbar_) Logger::log().debug() << " --  GEN  --" << genp_selector_.ttbar_final_system() << std::endl;
    if(isTTbar_ && genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) {
      matched_perm = matcher_.match(
        genp_selector_.ttbar_final_system(),
        object_selector_.clean_jets(), 
        object_selector_.veto_electrons(),
        object_selector_.veto_muons()
        );
      matched_perm.SetMET(object_selector_.met());
    }
    tracker_.track("matched perm", leptype);

		stringstream evtdir;		
		evtdir << leptype << "/";
		if(isTTbar_) {
			TTNaming dir_id = get_ttdir_name(matched_perm, best_permutation);
			evtdir << naming_.at(dir_id) << "/";			
		}
		evtdir << sys_name << "/";
		switch(njets) {
		case 3: evtdir << "3jets"; break;
		case 4: evtdir << "4jets"; break;
		default:evtdir << "5Pjets"; break;
		}
		fill_selection_plots(evtdir.str(), best_permutation);
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
			evt_idx_--;
			if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
			evt_idx_++;

			//long and time consuming
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
		} //while(event.next())
   }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  virtual void end(){
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
