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
#include "RObject.h"
#include "systematics.h"
#include "MCWeightProducer.h"
#include "LeptonSF.h"
#include <unordered_map>
#include "TwoBodyDecay.h"
#include "Hypotheses.h"
using namespace std;

class topspin_gen : public AnalyzerBase
{
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  //histos
  unordered_map<string, unordered_map<string, RObject> > histos_;

	//switches
	bool isTTbar_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTGenParticleSelector genp_selector_pst_;
  TTObjectSelector object_selector_;
  TTGenMatcher matcher_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  TTBarSolver solver_;
	//Scale factors
  LeptonSF electron_sf_, muon_sf_;
  TTPermutator permutator_;
  
public:
  topspin_gen(const std::string output_filename):
    AnalyzerBase("topspin_gen", output_filename),
    histos_(),
		tracker_(),
    object_selector_(),
    matcher_(),
    evt_weight_(1.),
    mc_weights_(),
    solver_(),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
    genp_selector_(),
    genp_selector_pst_(TTGenParticleSelector::SelMode::PSEUDOTOP),
    permutator_()
  {
    
    //tracker_.verbose(true);    
    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

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
    
    mc_weights_.init(sample);

    //Init solver
    string filename = "prob_ttJets.root";
    Logger::log().debug() << "solver file: " << filename << endl;
    TFile probfile(DataFile(filename).path().c_str());
    TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
    solver_.Init(td, false, true, true);

  };
  
  ~topspin_gen() {
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

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
    vector<double> mass_binning = {200,400,600,800,1000,2000}; 
    for(string pdir : {"parton", "part_acceptance", "matched", "wrong"}) {
      // for(string charge : {"lplus", "lminus"}) {
      string dir = pdir;//+"/"+charge;
      //cos theta (tt | tt)
      book<TH2F>(dir, "costheta_tops", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      //cos theta (jj | W)
      //book<TH2F>(dir, "costheta_jj", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      //book<TH2F>(dir, "costheta_ln", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);

      book<TH2F>(dir, "costheta_wl", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "costheta_wh", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);

      book<TH2F>(dir, "costheta_wl_wh_M200_400", "", 200, -1, 1, 200, -1, 1);
      book<TH2F>(dir, "costheta_wl_wh_M400_600", "", 200, -1, 1, 200, -1, 1);
      book<TH2F>(dir, "costheta_wl_wh_M600_800", "", 200, -1, 1, 200, -1, 1);
      book<TH2F>(dir, "costheta_wl_wh_M800Plus", "", 200, -1, 1, 200, -1, 1);

      //cos theta (bW | t)
      //book<TH2F>(dir, "costheta_bWl", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      //book<TH2F>(dir, "costheta_bWh", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);

      book<TH2F>(dir, "cosgamma_bWtlep", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "cosgamma_bWthad", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      //cos alpha bW -- bWbar planes
      book<TH2F>(dir, "cosalpha_tbW_tbarbW", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);

      //helicity and lab frames from TOP 14 023
      book<TH2F>(dir, "helframe_costheta_lep", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_costheta_nu" , "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_costheta_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_costheta_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "labframe_cosdeltaphi_lep_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "labframe_deltaphi_lep_dtype", "", 200, 0, Pi(), mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_cosdelta_lep_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_prodcosth_lep_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);      
      book<TH2F>(dir, "labframe_cosdeltaphi_lep_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "labframe_deltaphi_lep_utype", "", 200, 0, Pi(), mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_cosdelta_lep_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_prodcosth_lep_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);      

      book<TH2F>(dir, "labframe_cosdeltaphi_nu_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "labframe_deltaphi_nu_dtype", "", 200, 0, Pi(), mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_cosdelta_nu_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_prodcosth_nu_dtype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);      
      book<TH2F>(dir, "labframe_deltaphi_nu_utype", "", 200, 0, Pi(), mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "labframe_cosdeltaphi_nu_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_cosdelta_nu_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);
      book<TH2F>(dir, "helframe_prodcosth_nu_utype", "", 200, -1, 1, mass_binning.size()-1, &mass_binning[0]);      
      // }
    }
  }

  void fill(hyp::TTbar& ttbar, string dir, bool first_leptonic=true) {
    //string lcharge = (ttbar.lep_charge() > 0) ? "lplus" : "lminus";
    auto hdir = histos_.find(dir); //+"/"+lcharge);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }    
    
    auto leptopcm = ttbar.tlep().to_CM();
    auto hadtopcm = ttbar.thad().to_CM();
    auto ttcm = ttbar.to_CM();
    double mass = ttbar.M();

    //
    //Top pair CM-based variables
    //
    hdir->second["costheta_tops"].fill(ttbar.decay_opening_cm(), mass);    
    hdir->second["cosalpha_tbW_tbarbW"].fill(ttcm.top().decay_plane().Dot(ttcm.tbar().decay_plane()), mass);
    

    //
    //top-based variables
    //
    hdir->second["cosgamma_bWtlep"].fill(leptopcm.W().decay_plane().Dot(leptopcm.b().unit3D()), mass);
    hdir->second["cosgamma_bWtbar"].fill(hadtopcm.W().decay_plane().Dot(hadtopcm.b().unit3D()), mass);

    //
    //W-based
    //
    //W plus
    double ctheta_l = ttbar.tlep().W().decay_opening_cm();
    double ctheta_h = ttbar.thad().W().decay_opening_cm();
    hdir->second["costheta_wl"].fill(ctheta_l, mass);
    hdir->second["costheta_wh"].fill(ctheta_h, mass);
    
    if(mass < 400) hdir->second["costheta_wl_wh_M200_400"].fill(ctheta_l, ctheta_h);
    if(mass < 600) hdir->second["costheta_wl_wh_M400_600"].fill(ctheta_l, ctheta_h);
    if(mass < 800) hdir->second["costheta_wl_wh_M600_800"].fill(ctheta_l, ctheta_h);
    else hdir->second["costheta_wl_wh_M800Plus"].fill(ctheta_l, ctheta_h);
    
    //helicity and lab frames from TOP 14 023
    double cth_l = leptopcm.W().l().unit3D().Dot(ttcm.tlep().unit3D());
    double cth_n = leptopcm.W().nu().unit3D().Dot(ttcm.tlep().unit3D());
    double cth_d = hadtopcm.W().down().unit3D().Dot(ttcm.thad().unit3D());
    double cth_u = hadtopcm.W().up().unit3D().Dot(ttcm.thad().unit3D());

    if(cth_l != cth_l || cth_n != cth_n || cth_d != cth_d || cth_u != cth_u) {
      cout <<dir << " NANERROR! "<< cth_l << " " << cth_n << " " << cth_d << " " << cth_u << endl;
    //   cout << "input top:" << endl <<
    //     "TTbar: " << ttbar << endl <<
    //     "  TLep: " << ttbar.tlep() << endl <<
    //     "    b: " << ttbar.tlep().b() << endl <<
    //     "    W: " << ttbar.tlep().W() << endl <<
    //     "      l : " << ttbar.tlep().W().l() << endl << 
    //     "      nu: " << ttbar.tlep().W().nu() << endl << 
    //     "  THad: " << ttbar.thad() << endl <<
    //     "    b: " << ttbar.thad().b() << endl <<
    //     "    W: " << ttbar.thad().W() << endl <<
    //     "      up: " << ttbar.thad().W().up() << endl << 
    //     "      down: " << ttbar.thad().W().down() << endl;

    //   cout << "components" << endl <<
    //     "ttcm.tlep() " << ttcm.tlep() << endl <<
    //     "ttcm.thad() " << ttcm.thad() << endl <<
    //     "leptopcm.W().l() " << leptopcm.W().l() << endl <<
    //     "leptopcm.W().nu() " << leptopcm.W().nu() << endl <<
    //     "hadtopcm.W().down() " << hadtopcm.W().down() << endl <<
    //     "hadtopcm.W().up() " << hadtopcm.W().up() << endl;
    //   throw 40;
    }
    hdir->second["helframe_costheta_lep"  ].fill(cth_l, mass);
    hdir->second["helframe_costheta_nu"   ].fill(cth_n, mass);
    hdir->second["helframe_costheta_dtype"].fill(cth_d, mass);
    hdir->second["helframe_costheta_utype"].fill(cth_u, mass);
    hdir->second["helframe_cosdelta_lep_dtype"].fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()), mass);
    hdir->second["helframe_cosdelta_lep_utype"].fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D()), mass);
    hdir->second["helframe_prodcosth_lep_dtype"].fill(cth_l*cth_d, mass);      
    hdir->second["helframe_prodcosth_lep_utype"].fill(cth_l*cth_u, mass);      

    hdir->second["helframe_cosdelta_nu_dtype"].fill(leptopcm.W().nu().unit3D().Dot(hadtopcm.W().down().unit3D()), mass);
    hdir->second["helframe_cosdelta_nu_utype"].fill(leptopcm.W().nu().unit3D().Dot(hadtopcm.W().up().unit3D()), mass);
    hdir->second["helframe_prodcosth_nu_dtype"].fill(cth_n*cth_d, mass);      
    hdir->second["helframe_prodcosth_nu_utype"].fill(cth_n*cth_u, mass);      

    auto lep = ttbar.tlep().W().l();
    auto nu = ttbar.tlep().W().nu();
    auto up = ttbar.thad().W().up();
    auto dw = ttbar.thad().W().down();    
    hdir->second["labframe_cosdeltaphi_lep_dtype"].fill(Cos(lep.DeltaPhi(up)), mass);
    hdir->second["labframe_cosdeltaphi_lep_utype"].fill(Cos(lep.DeltaPhi(dw)), mass);
    hdir->second["labframe_deltaphi_lep_dtype"].fill(lep.DeltaPhi(up), mass);
    hdir->second["labframe_deltaphi_lep_utype"].fill(lep.DeltaPhi(dw), mass);

    hdir->second["labframe_cosdeltaphi_nu_dtype"].fill(Cos(nu.DeltaPhi(up)), mass);
    hdir->second["labframe_cosdeltaphi_nu_utype"].fill(Cos(nu.DeltaPhi(dw)), mass);
    hdir->second["labframe_deltaphi_nu_dtype"].fill(nu.DeltaPhi(up), mass);
    hdir->second["labframe_deltaphi_nu_utype"].fill(nu.DeltaPhi(dw), mass);
  }

  void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) {
    //select reco objects
    GenTTBar &ttbar = genp_selector_.ttbar_system();
    if(!(ttbar.type == GenTTBar::DecayType::SEMILEP)) return;

    hyp::TTbar ttgen(ttbar);
    fill(ttgen, "parton");

    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("obj selection");

    //find mc weight
    // if(object_selector_.tight_muons().size() == 1)
    //   evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    // if(object_selector_.medium_electrons().size() == 1)
    //   evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
    // tracker_.track("MC weights");

    //Gen matching
    Permutation matched;
    if(ttbar.type == GenTTBar::DecayType::SEMILEP) {
      matched = matcher_.match(
        genp_selector_.ttbar_final_system(),
        object_selector_.clean_jets(), 
        object_selector_.loose_electrons(),
        object_selector_.loose_muons()
        );
    }
    matched.SetMET(object_selector_.met());

    if(!matched.IsComplete()) return;
    matched.Solve(solver_);
    tracker_.track("gen matching");

    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;

    //fill gen
    hyp::TTbar ttgenfinal(genp_selector_.ttbar_final_system());
    fill(ttgenfinal, "part_acceptance");

    // //fill pst
    // //TwoBodyDecay ttpst(genp_selector_pst_.ttbar_final_system());
    // //fill(ttgen, "pseudotop");

    //fil reco
    hyp::TTbar ttreco(matched);
    fill(ttreco, "matched");

    //Find best permutation
    bool go_on = true;
    Permutation best_permutation;
    while(go_on) {
      Permutation test_perm = permutator_.next(go_on);
      test_perm.LepCharge(object_selector_.lepton_charge());
      if(go_on) {
        test_perm.Solve(solver_);
        if(test_perm.IsCorrect(matched)) continue;
        TLorentzVector jsum = *test_perm.WJa()+*test_perm.WJb();        
        if(Abs(jsum.M() - 75) > 25) continue;
        jsum += *test_perm.BHad();
        if(Abs(jsum.M() - 150) > 50) continue;
        hyp::TTbar ttwrong(test_perm);
        fill(ttwrong, "wrong");
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
        //genp_selector_pst_.select(event);
        tracker_.track("gen selection");        
        if(!selection) {
          Logger::log().error() << "Error: TTGenParticleSelector was not able to find all the generated top decay products in event " << evt_idx << endl <<
            "run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
          continue;
        }
			}

      // evt_weight_ = mc_weights_.evt_weight(event, systematics::SysShifts::NOSYS);
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
  URDriver<topspin_gen> test;
  return test.run();
}
