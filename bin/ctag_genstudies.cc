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

class ctag_genstudies : public AnalyzerBase
{
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  //histos
  unordered_map<string, RObject> histos_;

	//switches
	bool isTTbar_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTGenMatcher matcher_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;

	//Scale factors
  LeptonSF electron_sf_, muon_sf_;
  
public:
  ctag_genstudies(const std::string output_filename):
    AnalyzerBase("ctag_genstudies", output_filename),
    histos_(),
		tracker_(),
    object_selector_(),
    matcher_(),
    evt_weight_(1.),
    mc_weights_(),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"){
    
    //tracker_.verbose(true);    
    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);

    opts::variables_map &values = URParser::instance().values();

    //find out which sample are we running on
		string output_file = values["output"].as<std::string>();
    //DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
		string sample = systematics::get_sample(output_file);
		isTTbar_ = boost::starts_with(sample, "ttJets");

    if(!isTTbar_) {
      Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
      throw 49;
    }
    
    //Do not init the solver, as we do not have the files! 
    if(values["general.pseudotop"].as<int>() == 0) genp_selector_ = TTGenParticleSelector(); //FIXME allow for herwig setting
    else genp_selector_ = TTGenParticleSelector(TTGenParticleSelector::SelMode::PSEUDOTOP);

    mc_weights_.init(sample);
  };
  
  ~ctag_genstudies() {
  }

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
		histos_["Whad_jet_pts"   ] = RObject::book<TH2F>("Whad_jet_pts"    , ";lead pT; sublead pT", 100, 0., 500., 100, 0., 500.);
		histos_["BJet_jet_pts"   ] = RObject::book<TH2F>("BJet_jet_pts"    , ";lead pT; sublead pT", 100, 0., 500., 100, 0., 500.);
		histos_["leadB_leadW_pts"] = RObject::book<TH2F>("leadB_leadW_pts" , ";lead pT; sublead pT", 100, 0., 500., 100, 0., 500.);
		histos_["subB_subW_pts"  ] = RObject::book<TH2F>("subB_subW_pts"   , ";lead pT; sublead pT", 100, 0., 500., 100, 0., 500.);
		histos_["subB_leadW_pts" ] = RObject::book<TH2F>("subB_leadW_pts"  , ";lead pT; sublead pT", 100, 0., 500., 100, 0., 500.);
		histos_["bjets_BTag" ]  = RObject::book<TH2F>("bjets_BTag"  , ";best BTag; worst BTag", 4, -0.5, 3.5, 4, -0.5, 3.5);
    histos_["pt_first_idx"] = RObject::book<TH1D>("pt_first_idx", ";idx;# Events", 41, -0.5, 40.5);
    histos_["pt_last_idx" ] = RObject::book<TH1D>("pt_last_idx" , ";idx;# Events", 41, -0.5, 40.5);
  }

  int get_btag(const IDJet *jet) {
    if(jet->BTagId(IDJet::BTag::CSVTIGHT)) return 3;
    else if(jet->BTagId(IDJet::BTag::CSVMEDIUM)) return 2;
    else if(jet->BTagId(IDJet::BTag::CSVLOOSE)) return 1;
    else return 0;
  }
  
  void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) {
    //select reco objects
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("obj selection");

    //find mc weight
    if(object_selector_.tight_muons().size() == 1)
      evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    if(object_selector_.tight_electrons().size() == 1)
      evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    tracker_.track("MC weights");

    //Gen matching
    GenTTBar &ttbar = genp_selector_.ttbar_system();
    Permutation matched;
    if(ttbar.type == GenTTBar::DecayType::SEMILEP) {
      matched = matcher_.match(
        genp_selector_.ttbar_final_system(),
        object_selector_.clean_jets(), 
        object_selector_.veto_electrons(),
        object_selector_.veto_muons()
        );
    }
    matched.SetMET(object_selector_.met());
    tracker_.track("gen matching");
    if(!matched.IsComplete()) return;

    const IDJet *lead_wj = (matched.WJa()->Pt() > matched.WJb()->Pt()) ? matched.WJa() : matched.WJb();
    const IDJet *subl_wj = (matched.WJa()->Pt() > matched.WJb()->Pt()) ? matched.WJb() : matched.WJa();
    const IDJet *lead_bj = (matched.BLep()->Pt() > matched.BHad()->Pt()) ? matched.BLep() : matched.BHad();
    const IDJet *subl_bj = (matched.BLep()->Pt() > matched.BHad()->Pt()) ? matched.BHad() : matched.BLep();

		histos_["Whad_jet_pts"   ].fill( lead_wj->Pt(), subl_wj->Pt());
		histos_["BJet_jet_pts"   ].fill( lead_bj->Pt(), subl_bj->Pt());
		histos_["leadB_leadW_pts"].fill( lead_bj->Pt(), lead_wj->Pt());
		histos_["subB_subW_pts"  ].fill( subl_bj->Pt(), subl_wj->Pt());
		histos_["subB_leadW_pts" ].fill( subl_bj->Pt(), lead_wj->Pt());

    auto jets = object_selector_.clean_jets();
    sort(jets.begin(), jets.end(), [](IDJet* A, IDJet* B){return(A->Pt() > B->Pt());});
    int strt=9999;
    int stop=-1;
    vector<const IDJet*> cmb_jets = {lead_wj, subl_wj, lead_bj, subl_bj};
    for(int idx=0; idx<jets.size(); idx++) {
      if(any_of(cmb_jets.begin(), cmb_jets.end(), [&] (const IDJet*i) {return jets[idx] == i;}) ) {
        strt = Min(strt, idx); 
        stop = Max(stop, idx); 
      }
    }
    histos_["pt_first_idx"].fill(strt);
    histos_["pt_last_idx" ].fill(stop);
    
    int btagl = get_btag(lead_bj);
    int btags = get_btag(subl_bj);
    
		histos_["bjets_BTag"].fill(Max(btagl, btags), Min(btagl, btags));

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

      evt_weight_ = mc_weights_.evt_weight(event, systematics::SysShifts::NOSYS);
      process_evt(event);
    }//while(event.next())
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
  URDriver<ctag_genstudies> test;
  return test.run();
}
