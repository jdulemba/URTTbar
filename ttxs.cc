#include <iostream>
#include "AnalyzerBase.h"
#include "URStreamer.h"
#include "URDriver.h"
#include "CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
#include "URParser.h"
#include "DataFile.h"
#include "TTBarResponse.h"
#include "TTObjectSelector.h"
#include "TTBarPlots.h"
#include "TTBarSolver.h"
#include "TTGenParticleSelector.h"
#include "TTPermutator.h"
#include "TTGenMatcher.h"
#include "TRandom3.h"
#include "helper.h"
//#include <map>

using namespace std;

class ttxs : public AnalyzerBase
{
public:
	enum TTNaming {NOTSET, RIGHT, RIGHT_THAD, RIGHT_TLEP, WRONG, OTHER};
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  vector<string> dir_names_ = {"RECO", "semilep_visible_right", "semilep_right_thad", "semilep_right_tlep", "semilep_wrong", "other_tt_decay"};
  map<TTNaming, map<TTObjectSelector::SysShifts, TTBarPlots> > ttplots_; //not optimal for write one read many, but still (better unordered map)
  map<TTObjectSelector::SysShifts, TTBarResponse> responses_;
  //binning vectors
  vector<double> topptbins_;
  vector<double> topybins_;
  vector<double> ttmbins_;
  vector<double> ttybins_;
  vector<double> ttptbins_;
  vector<double> metbins_;
  vector<double> jetbins_;
  vector<double> nobins_;

	//switches
	bool isData_, isTTbar_, skew_tt_distro_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTPermutator permutator_;
  TTGenMatcher matcher_;
  TTBarSolver solver_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);

	vector<TTObjectSelector::SysShifts> systematics_;

	//Scale factors
	TH1D *electron_sf_, *muon_sf_;
	TH1D *pu_sf_;
  
public:
  ttxs(const std::string output_filename):
    AnalyzerBase("ttxs", output_filename),
		tracker_(),
    object_selector_(),
    permutator_(),
    matcher_(),
    solver_(),
    evt_weight_(1.) {

    //set tracker
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

    opts::variables_map &values = URParser::instance().values();

    //switches
    skew_tt_distro_ = values["general.skew_ttpt_distribution"].as<int>();

    //find out which sample are we running on
		string output_file = values["output"].as<std::string>();
    DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
		size_t slash = output_file.rfind("/") + 1;
		string basename(output_file, slash, output_file.size() - slash);
		isData_  = boost::starts_with(basename, "data");
		isTTbar_ = boost::starts_with(basename, "ttJets");
    Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;

    //choose systematics to run based on sample
    vector<TTObjectSelector::SysShifts> nosys = {TTObjectSelector::SysShifts::NOSYS};
    vector<TTObjectSelector::SysShifts> allsys = {
          TTObjectSelector::SysShifts::NOSYS, TTObjectSelector::SysShifts::JES_UP, TTObjectSelector::SysShifts::JES_DW, 
          TTObjectSelector::SysShifts::JER_UP, TTObjectSelector::SysShifts::MET_UP, TTObjectSelector::SysShifts::MET_DW}; //FIXME JER_DW still to implement
    if(isData_) systematics_ = nosys;
    else {
      if(isTTbar_) { //not all tt samples need the ful sys, as they are ALREADY a systematic!
        if(basename.find("mtopdown") != string::npos ||  //FIXME: add the samples as they come
           basename.find("mtopup") != string::npos || 
           basename.find("scaledown") != string::npos || 
           basename.find("scaleup") != string::npos ) systematics_ = nosys;
        else systematics_ = allsys;
      }
      else systematics_ = allsys;
    }
    
    solver_.Init(solver_input.path(), false, true, true);
    if(values["general.pseudotop"].as<int>() == 0) genp_selector_ = TTGenParticleSelector(); //FIXME allow for herwig setting
    else genp_selector_ = TTGenParticleSelector(TTGenParticleSelector::SelMode::PSEUDOTOP);

    DataFile sf_filename(values["general.lepton_sf"].as<std::string>());
    Logger::log().debug() << "lepton sf file: " << sf_filename.path() << endl;
    TFile* sf_file = TFile::Open(sf_filename.path().c_str());
    TH1::AddDirectory(false);
    electron_sf_ = (TH1D*) ((TH1D*)sf_file->Get("Scale_ElTOT_Pt"))->Clone("electron_sf");
    muon_sf_ = (TH1D*) ((TH1D*)sf_file->Get("Scale_MuTOT_Pt"))->Clone("muon_sf");

    DataFile pu_filename("PUweight.root"); //FIXME, use better recipe
    Logger::log().debug() << "PU sf file: " << pu_filename.path() << endl;
    TFile* pu_file = TFile::Open(pu_filename.path().c_str()); 
    pu_sf_ = (TH1D*) ((TH1D*)pu_file->Get("PUweight"))->Clone("pu_weight");
    TH1::AddDirectory(true);
    
    //Set binning
    Logger::log().debug() << "-- Setting binning --" << endl;
    setbinning(topptbins_, 0., 800., 5.);
    setbinning(topybins_, 0., 2.5, 0.1);
    setbinning(ttmbins_, 250., 2000., 5.);
    setbinning(ttptbins_, 0., 500., 5.);
    setbinning(ttybins_, 0., 2.5, 0.1);
    setbinning(metbins_, 0., 600, 10.);
    jetbins_ = {-0.5, 0.5, 1.5, 2.5, 10.};
    nobins_ = {0., 13000.};
  };
  
  ~ttxs() {
    delete electron_sf_;
    delete muon_sf_;
    delete pu_sf_;
  }

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
    vector<TTNaming> evt_types;
    if(isTTbar_) evt_types = {RIGHT, RIGHT_THAD, RIGHT_TLEP, WRONG, OTHER};
    else evt_types = {NOTSET};

    for(auto evt_type : evt_types){
      TDirectory* dir_type = outFile_.mkdir(dir_names_[evt_type].c_str());
      dir_type->cd();
      for(auto shift : systematics_) {
        TDirectory* dir_sys = dir_type->mkdir(TTObjectSelector::shift_to_name.at(shift).c_str());
        dir_sys->cd();
        ttplots_[evt_type][shift];
        ttplots_[evt_type][shift].Init(topptbins_, topybins_, ttmbins_,
                                       ttybins_  , ttptbins_, metbins_,
                                       jetbins_  , nobins_);
        
        //create response matrixes
        if(isTTbar_ && evt_type == RIGHT){
          responses_[shift] = TTBarResponse("response");
          responses_[shift].AddMatrix("thadpt", topptbins_, topptbins_, "p_{T}(t_{h}) [GeV]");
          responses_[shift].AddMatrix("thady", topybins_, topybins_, "|y(t_{h})|");
          responses_[shift].AddMatrix("tleppt", topptbins_, topptbins_, "p_{T}(t_{l}) [GeV]");
          responses_[shift].AddMatrix("tlepy", topybins_, topybins_, "|y(t_{l})|");
          responses_[shift].AddMatrix("ttm", ttmbins_, ttmbins_, "m(t#bar{t}) [GeV]");
          responses_[shift].AddMatrix("ttpt", ttptbins_, ttptbins_, "p_{T}(t#bar{t}) [GeV]");
          responses_[shift].AddMatrix("tty", ttybins_, ttybins_, "|y(t#bar{t})|");
          responses_[shift].AddMatrix("njet", jetbins_, jetbins_, "n-jets");
          responses_[shift].AddMatrix("nobin", nobins_, nobins_, "total");
        }
      }
    }
  }
  
  void process_evt(URStreamer &event, TTObjectSelector::SysShifts shift=TTObjectSelector::SysShifts::NOSYS) {
    tracker_.track("start");
    //Fill truth of resp matrix
    if(isTTbar_ && genp_selector_.is_in_acceptance()) {
      GenTTBar &ttbar = genp_selector_.ttbar_system();

      if(skew_tt_distro_) evt_weight_ *= 1.+0.05*(ttbar.top.Pt()-200.)/1000.;
      size_t added_jets = genp_selector_.additional_jets().size();
      responses_[shift].FillTruth("thadpt", ttbar.had_top()->Pt(), added_jets, evt_weight_);
      responses_[shift].FillTruth("nobin", ttbar.had_top()->Pt(), added_jets, evt_weight_);
      responses_[shift].FillTruth("thady", Abs(ttbar.had_top()->Rapidity()), added_jets, evt_weight_);
      responses_[shift].FillTruth("tleppt", ttbar.lep_top()->Pt(), added_jets, evt_weight_);
      responses_[shift].FillTruth("tlepy", Abs(ttbar.lep_top()->Rapidity()), added_jets, evt_weight_);
      responses_[shift].FillTruth("ttm", ttbar.M(), added_jets, evt_weight_);
      responses_[shift].FillTruth("ttpt", ttbar.Pt(), added_jets, evt_weight_);
      responses_[shift].FillTruth("tty", Abs(ttbar.Rapidity()), added_jets, evt_weight_);
      responses_[shift].FillTruth("njet", genp_selector_.additional_jets().size(), added_jets, evt_weight_);
    }

    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("object selection");
    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;
    tracker_.track("perm preselection");

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
      
    //Find best permutation
    bool go_on = true;
    Permutation best_perm;
    size_t ncycles_=0;
    while(go_on) {
      ncycles_++;
      Permutation test_perm = permutator_.next(go_on);
      if(go_on) {
        test_perm.Solve(solver_);
        if(test_perm.MassDiscr() < best_perm.MassDiscr())	{
          best_perm = test_perm;
        }
      }
    }
    size_t capped_jets_size = permutator_.capped_jets().size();
    best_perm.permutating_jets(capped_jets_size);

    //Logger::log().debug() << "ncycles: " << ncycles_ << ", best_perm.L(): " << best_perm.L() << endl;
    //find mc weight
    if(!isData_) { //FIXME, use real PU!
      double npu = event.vertexs().size();
      if(npu > 0)	evt_weight_ *= pu_sf_->GetBinContent(pu_sf_->FindFixBin(npu));		
      if(object_selector_.tight_muons().size() == 1)
        evt_weight_ *= muon_sf_->GetBinContent(muon_sf_->FindFixBin(Min(best_perm.L()->Pt(), 95.)));
      if(object_selector_.medium_electrons().size() == 1)
        evt_weight_ *= electron_sf_->GetBinContent(electron_sf_->FindFixBin(Min(best_perm.L()->Pt(), 95.)));
		}

    //Fill appropriate plots
    if(!isTTbar_) {
      ttplots_[TTNaming::NOTSET][shift].Fill(best_perm, object_selector_, evt_weight_);
    } //if(!isTTbar_)
    else {
      // if(matched_perm.IsComplete())
      //   Logger::log().debug() << "-- Matched -- " << matched_perm << endl <<
      //     "-- Best --" << best_perm << endl;
      if(best_perm.IsCorrect(matched_perm)){
        ttplots_[TTNaming::RIGHT][shift].Fill(best_perm, object_selector_, evt_weight_);
        responses_[shift].FillReco("thadpt", best_perm.THad().Pt(), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("nobin", best_perm.THad().Pt(), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("thady", Abs(best_perm.THad().Rapidity()), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("tleppt", best_perm.TLep().Pt(), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("tlepy", Abs(best_perm.TLep().Rapidity()), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("ttm", (best_perm.THad() + best_perm.TLep()).M(), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("ttpt", (best_perm.THad() + best_perm.TLep()).Pt(), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("tty", Abs((best_perm.THad() + best_perm.TLep()).Rapidity()), capped_jets_size, evt_weight_);
        responses_[shift].FillReco("njet", object_selector_.clean_jets().size() - 4, capped_jets_size, evt_weight_);
      }
      else if(best_perm.IsTHadCorrect(matched_perm)) ttplots_[TTNaming::RIGHT_THAD][shift].Fill(best_perm, object_selector_, evt_weight_);
      else if(best_perm.IsBLepCorrect(matched_perm)) ttplots_[TTNaming::RIGHT_TLEP][shift].Fill(best_perm, object_selector_, evt_weight_);
      else if(genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) ttplots_[TTNaming::WRONG][shift].Fill(best_perm, object_selector_, evt_weight_);
      else ttplots_[TTNaming::OTHER][shift].Fill(best_perm, object_selector_, evt_weight_);
    } //else -- if(!isTTbar_)
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

		tracker_.deactivate();
    while(event.next()) {
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			if(evt_idx % 1000 == 0) Logger::log().debug() << "Beginning event " << evt_idx << endl;
      // Logger::log().debug() << " *************************************************** " << endl;
      // Logger::log().debug() << "Beginning event " << evt_idx << endl;

			//long and time consuming
			if(isTTbar_){
				genp_selector_.select(event);			
			}

			for(auto shift : systematics_){
        evt_weight_ = 1;
				if(shift == TTObjectSelector::SysShifts::NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(event, shift);
        responses_[shift].Flush();
				if(shift == TTObjectSelector::SysShifts::NOSYS) tracker_.deactivate();
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
    parser.addCfgParameter<string>("general", "ttsolver_input", "ttsolver input file");
    parser.addCfgParameter<string>("general", "lepton_sf", "lepton SF input file");
    parser.addCfgParameter<int>("general", "pseudotop", "should I use pseudo-top? (0/1)");

		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");
  }
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<ttxs> test;
  return test.run();
}
