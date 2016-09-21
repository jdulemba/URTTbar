#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
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
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "TROOT.h"
#include <algorithm>
//#include <map>

using namespace std;

class permProbComputer : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_TLEP, WRONG, OTHER, NOTSET};
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  vector<string> dir_names_ = {"semilep_visible_right", "semilep_right_thad", "semilep_right_tlep", "semilep_wrong", "other_tt_decay"};
  //histos
  map<systematics::SysShifts, map<TTNaming, map<string, RObject> > > histos_;


	//switches
	bool isTTbar_, skew_tt_distro_=false;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTPermutator permutator_;
  TTGenMatcher matcher_;
  TTBarSolver solver_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;

	vector<systematics::SysShifts> systematics_;

	//Scale factors
  LeptonSF electron_sf_, muon_sf_;
  
public:
  permProbComputer(const std::string output_filename):
    AnalyzerBase("permProbComputer", output_filename),
		tracker_(),
    object_selector_(),
    permutator_(),
    matcher_(),
    solver_(),
    evt_weight_(1.),
    mc_weights_(),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"){
    
    //tracker_.verbose(true);    
    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

    opts::variables_map &values = URParser::instance().values();

    //switches
    skew_tt_distro_ = values["general.skew_ttpt_distribution"].as<int>();

    //find out which sample are we running on
		string output_file = values["output"].as<std::string>();
    //DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
		string sample = systematics::get_sample(output_file);
    bool isSignal = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
		isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal;

    //choose systematics to run based on sample
    systematics_ = {systematics::SysShifts::NOSYS};
    if(!isTTbar_) {
      Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
      throw 49;
    }
    
    //Do not init the solver, as we do not have the files! 
    if(values["general.pseudotop"].as<int>() == 0) genp_selector_ = TTGenParticleSelector(); //TTGenParticleSelector::SelMode::LHE); //FIXME allow for herwig setting
    else genp_selector_ = TTGenParticleSelector(TTGenParticleSelector::SelMode::PSEUDOTOP);

    solver_.Init("", false, false, false);

    mc_weights_.init(sample);
  };
  
  ~permProbComputer() {
  }

	template <class T> 
	pair<const T,const T> Minmax(const T& a, const T& b) {
		return (b<a) ? std::make_pair(b, min(a, .99999)) : std::make_pair(a, min(b, .99999));
	}

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
    vector<TTNaming> evt_types = {RIGHT, RIGHT_THAD, RIGHT_TLEP, WRONG, OTHER};
    for(auto evt_type : evt_types){
      TDirectory* dir_type = outFile_.mkdir(dir_names_[evt_type].c_str());
      dir_type->cd();
      for(auto shift : systematics_) {
        TDirectory* dir_sys = dir_type->mkdir(systematics::shift_to_name.at(shift).c_str());
        dir_sys->cd();
        //TODO add plots
        histos_[shift][evt_type]["mWhad_vs_mtophad"] = RObject::book<TH2D>("mWhad_vs_mtophad", ";M(W_{had}) [GeV];M(t_{had}) [GeV]", 500, 0., 500., 500, 0., 500);
        histos_[shift][evt_type]["nusolver_chi2"] = RObject::book<TH1D>("nusolver_chi2", "#chi^{2};# Events", 75, 0., 150.);
        histos_[shift][evt_type]["wjets_cMVA"] = RObject::book<TH2D>("wjets_cMVA", "", 100, -1., 1, 100, -1., 1);
        histos_[shift][evt_type]["bjets_cMVA"] = RObject::book<TH2D>("bjets_cMVA", "", 100, -1., 1, 100, -1., 1);
        histos_[shift][evt_type]["wjets_qgt"] = RObject::book<TH2D>("wjets_qgt", "", 50, 0., 1, 50, 0., 1);
        histos_[shift][evt_type]["bjets_qgt"] = RObject::book<TH2D>("bjets_qgt", "", 50, 0., 1, 50, 0., 1);
        histos_[shift][evt_type]["wjets_bqgt"] = RObject::book<TH1D>("wjets_bqgt", "", 50, 0., 1);
        histos_[shift][evt_type]["wjets_wqgt"] = RObject::book<TH1D>("wjets_wqgt", "", 50, 0., 1);
        histos_[shift][evt_type]["wjets_cMVA_WP"] = RObject::book<TH2D>("wjets_cMVA_WP", ";M(W_{had}) [GeV];M(t_{had}) [GeV]", 4, 0., 4., 4, 0., 4.);
        histos_[shift][evt_type]["bjets_cMVA_WP"] = RObject::book<TH2D>("bjets_cMVA_WP", ";M(W_{had}) [GeV];M(t_{had}) [GeV]", 4, 0., 4., 4, 0., 4.);
        histos_[shift][evt_type]["lb_ratio"] = RObject::book<TH1D>("lb_ratio", "", 100, 0., 5.);
        histos_[shift][evt_type]["w1b_ratio"] = RObject::book<TH1D>("w1b_ratio", "", 100, 0., 5.);
        histos_[shift][evt_type]["w2b_ratio"] = RObject::book<TH1D>("w2b_ratio", "", 100, 0., 5.);
        histos_[shift][evt_type]["lb_w2b_ratio"] = RObject::book<TH2D>("lb_w2b_ratio", "", 100, 0., 5., 100, 0., 5.);
      }
    }
  }

	int btag_idval(const IDJet* jet) {
		if(jet->BTagId(IDJet::BTag::MVATIGHT)) return 3;
		else if(jet->BTagId(IDJet::BTag::MVAMEDIUM)) return 2;
		else if(jet->BTagId(IDJet::BTag::MVALOOSE) ) return 1;
		return 0;
	}
  
  void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) {
    //Fill truth of resp matrix
    GenTTBar &ttbar = genp_selector_.ttbar_system();
    if(skew_tt_distro_) evt_weight_ *= 1.+0.05*(ttbar.top.Pt()-200.)/1000.;
 
    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("obj selection");

    //find mc weight
    if(object_selector_.tight_muons().size() == 1)
      evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    if(object_selector_.tight_electrons().size() == 1)
      evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    tracker_.track("MC weights");

    bool preselection_pass = permutator_.preselection(
      object_selector_.clean_jets(), object_selector_.lepton(), object_selector_.met()
      );
    tracker_.track("permutation pre-selection done (not applied)");

    //get needed histo map
    auto plots = histos_.find(shift)->second;    

    if( !preselection_pass ) return;
    tracker_.track("perm preselection");

    //Gen matching
    Permutation matched_perm;
    if(ttbar.type == GenTTBar::DecayType::SEMILEP) {
      matched_perm = matcher_.match(
        genp_selector_.ttbar_final_system(),
        object_selector_.clean_jets(), 
        object_selector_.veto_electrons(),
        object_selector_.veto_muons()
        );
    }
    matched_perm.SetMET(object_selector_.met());
    tracker_.track("gen matching");
    
    //Find best permutation
    bool go_on = true;
    while(go_on) {
      Permutation test_perm = permutator_.next(go_on);
      if(go_on) {
				test_perm.LepCharge(object_selector_.lepton_charge());
        test_perm.Solve(solver_);
        TTNaming perm_status;
        if(test_perm.IsCorrect(matched_perm)) perm_status = TTNaming::RIGHT;
        else if(test_perm.IsTHadCorrect(matched_perm)) perm_status = TTNaming::RIGHT_THAD;
        else if(test_perm.IsBLepCorrect(matched_perm)) perm_status = TTNaming::RIGHT_TLEP;
        else if(ttbar.type == GenTTBar::DecayType::SEMILEP) perm_status = TTNaming::WRONG;
        else perm_status = TTNaming::OTHER;
        
        plots[perm_status]["mWhad_vs_mtophad"].fill(test_perm.WHad().M(), test_perm.THad().M(), evt_weight_);
        plots[perm_status]["nusolver_chi2"].fill(test_perm.NuChisq(), evt_weight_);

				auto b_mM = Minmax(test_perm.BHad()->CombinedMVA(), test_perm.BLep()->CombinedMVA());
				auto w_mM = Minmax(test_perm.WJa()->CombinedMVA(), test_perm.WJb()->CombinedMVA());
        plots[perm_status]["wjets_cMVA"].fill(pow(w_mM.first, 11), pow(w_mM.second, 11), evt_weight_);
        plots[perm_status]["bjets_cMVA"].fill(pow(b_mM.first, 11), pow(b_mM.second, 11), evt_weight_);

				auto b_mM_qg = Minmax(test_perm.BHad()->qgTag(), test_perm.BLep()->qgTag());
				auto w_mM_qg = Minmax(test_perm.WJa() ->qgTag(), test_perm.WJb() ->qgTag());
        plots[perm_status]["wjets_qgt"].fill(w_mM_qg.first, w_mM_qg.second, evt_weight_);
        plots[perm_status]["bjets_qgt"].fill(b_mM_qg.first, b_mM_qg.second, evt_weight_);
        plots[perm_status]["wjets_bqgt"].fill(w_mM_qg.second, evt_weight_);
        plots[perm_status]["wjets_wqgt"].fill(pow(w_mM_qg.first, 8),  evt_weight_);
				
				auto bwp_mM = Minmax(btag_idval(test_perm.BHad()), btag_idval(test_perm.BLep()));
				auto wwp_mM = Minmax(btag_idval(test_perm.WJa() ), btag_idval(test_perm.WJb() ));
        plots[perm_status]["wjets_cMVA_WP"].fill(wwp_mM.first, wwp_mM.second, evt_weight_);
        plots[perm_status]["bjets_cMVA_WP"].fill(bwp_mM.first, bwp_mM.second, evt_weight_);

				auto wpt_mM = Minmax(test_perm.WJa()->Pt(), test_perm.WJb()->Pt());
        plots[perm_status]["lb_ratio" ].fill(test_perm.L()->Pt()/test_perm.BLep()->Pt(), evt_weight_);
        plots[perm_status]["w1b_ratio"].fill(wpt_mM.first/test_perm.BHad()->Pt(), evt_weight_);
        plots[perm_status]["w2b_ratio"].fill(wpt_mM.second/test_perm.BHad()->Pt(), evt_weight_);
        plots[perm_status]["lb_w2b_ratio"].fill(test_perm.L()->Pt()/test_perm.BLep()->Pt(), wpt_mM.second/test_perm.BHad()->Pt(), evt_weight_);
			}
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
    int report = values["report"].as<int>();
    Logger::log().debug() << "-- DONE -- reporting every -- " << report << endl;
    while(event.next()) {
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			if(evt_idx % report == 0) Logger::log().debug() << "Beginning event " << evt_idx << " run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
      tracker_.track("start");		 

			//long and time consuming
			if(isTTbar_){
        bool selection = 	genp_selector_.select(event);			
        tracker_.track("gen selection");        
        if(!selection) {
          Logger::log().error() << "Error: TTGenParticleSelector was not able to find all the generated top decay products in event " << evt_idx << endl <<
            "run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
          continue;
        }
			}
	
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
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report", opts::value<int>()->default_value(10000), "report every in debug mode");
  }
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<permProbComputer> test;
	int excode = test.run();
	//Logger::log().debug() << "RUNNING DONE " << std::endl;
	auto files = gROOT->GetListOfFiles(); //make ROOT aware that some files do not exist, because REASONS
	Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl; //need to print out this otherwise ROOT loses its shit in 7.4.X (such I/O, much features)
  return excode;
}
