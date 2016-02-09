#include <iostream>
#include "AnalyzerBase.h"
#include "URStreamer.h"
#include "URDriver.h"
#include <map>
#include "RObject.h"
#include "IDMuon.h"
#include "IDElectron.h"
#include "IDJet.h"
#include <algorithm>
//#include "NeutrinoSolver.h"
#include "TTBarSolver.h"
#include <list>
#include "Logger.h"
#include "URParser.h"
#include "CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
#include "PDGID.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Permutation.h"
#include <set>
#include "IDMet.h"
#include "TUUID.h"   
#include "systematics.h"
#include "TTObjectSelector.h"
#include "TTBarPlots.h"
#include "TTBarSolver.h"
#include "TTGenParticleSelector.h"
#include "TTPermutator.h"
#include "TTGenMatcher.h"
#include "MCWeightProducer.h"
#include "BTagSFProducer.h"
#include "LeptonSF.h"
#include "DataFile.h"

using namespace TMath;

static map<string, IDJet::BTag> available_bjet_id = {
	{"none"     , IDJet::BTag::NONE},
	{"csvMedium", IDJet::BTag::CSVMEDIUM},
	{"csvLoose" , IDJet::BTag::CSVLOOSE},
	{"csvTight" , IDJet::BTag::CSVTIGHT},
};	

static map<string, std::function<bool(const Permutation &, const Permutation &)> > available_ordering = {
	{"full_discriminant", [](const Permutation &one, const Permutation &two) {return  one.Prob()  < two.Prob() ;}         },
	{"nu_chisq"         , [](const Permutation &one, const Permutation &two) {return  one.NuChisq() 	 < two.NuChisq() ;}	},
	{"nu_discriminant"	 , [](const Permutation &one, const Permutation &two) {return  one.NuDiscr() 	 < two.NuDiscr() ;}	},
	{"btag_discriminant", [](const Permutation &one, const Permutation &two) {return  one.BDiscr()  	 < two.BDiscr()  ;}	},
	{"mass_discriminant", [](const Permutation &one, const Permutation &two) {return  one.MassDiscr() < two.MassDiscr();} },
	};


class ctag_eff : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
private:
	//histograms and helpers
	map<string, map<string, RObject> > histos_;
	map<TTNaming, string> naming_;
	CutFlowTracker tracker_;

	//switches
	bool isData_, isTTbar_, run_jes_;

  //selectors and helpers
  TTGenParticleSelector genp_selector_;
  TTObjectSelector object_selector_;
  TTPermutator permutator_;
  TTGenMatcher matcher_;
  map<systematics::SysShifts, TTBarSolver> solvers_;
  float evt_weight_;
	TRandom3 randomizer_;// = TRandom3(98765);
  MCWeightProducer mc_weights_;
  BTagSFProducer btag_sf_;

	vector<systematics::SysShifts> systematics_;
	string cut_ordering_;
	std::function<bool(const Permutation &, const Permutation &)> ordering_fcn_;
	map<string, std::function<bool(const IDJet*)> > working_points_;

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

public:
  std::vector<systematics::SysShifts> get_systematics(std::string outname) {
    std::string sample = systematics::get_sample(outname);
    std::vector<systematics::SysShifts> full_sys = {
      systematics::SysShifts::NOSYS, systematics::SysShifts::JES_UP, systematics::SysShifts::JES_DW,
      systematics::SysShifts::BTAG_B_UP, systematics::SysShifts::BTAG_B_DW,
      systematics::SysShifts::BTAG_L_UP, systematics::SysShifts::BTAG_L_DW,
      systematics::SysShifts::BTAG_C_UP, systematics::SysShifts::BTAG_C_DW,
    };
    std::vector<systematics::SysShifts> nosys = {systematics::SysShifts::NOSYS};

		if(boost::starts_with(sample, "ttJets")) {
      if(sample.find("_") != std::string::npos) {
        //sys shifted sample!
        return nosys;
      }
      else return full_sys;
    }
    else if(!boost::starts_with(sample, "data")) {
      return full_sys;
    }
    else return nosys;
  }

  ctag_eff(const std::string output_filename):
    AnalyzerBase("ctag_eff", output_filename), 
		tracker_(),
		working_points_(),
    object_selector_(),
    permutator_(),
    matcher_(),
    solvers_(),
    mc_weights_(),
    btag_sf_(permutator_, 0.5, -1, -1), //float CTag +/- 50% the rest as SF table
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_() 
	{
		TUUID id;  
		randomizer_.SetSeed(id.Hash());   

    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

    //find out which sample are we running on
    opts::variables_map &values = URParser::instance().values();
		string output_file = values["output"].as<std::string>();
		string sample = systematics::get_sample(output_file);
		isData_  = boost::starts_with(sample, "data");
		isTTbar_ = boost::starts_with(sample, "ttJets");
    Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
    systematics_ = get_systematics(output_file);
    if(!isData_) mc_weights_.init(sample);

    //Init solver
    string filename = "prob_ttJets.root";
    Logger::log().debug() << "solver file: " << filename << endl;
    TFile probfile(DataFile(filename).path().c_str());
    // for(auto shift : systematics_) {
    TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
    //   if(!td) td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());      
    solvers_[systematics::SysShifts::NOSYS];
    solvers_[systematics::SysShifts::NOSYS].Init(td, false, true, true);
    // }

		//SET CUTS FROM CFG
		cut_ordering_ = URParser::instance().getCfgPar<string>("permutations.ordering");
		ordering_fcn_ = available_ordering[cut_ordering_];

		working_points_["notag"]     = [](const IDJet* jet) {return false;};
		working_points_["csvLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVLOOSE);};
		working_points_["csvTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVTIGHT);};
		working_points_["csvMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVMEDIUM);};
		// working_points_["rndm10"]    = [](const IDJet* jet) {return jet->rndm() < 0.1;};
		// working_points_["ssvHiPur"] = [](const Jet* jet) {return (bool) jet->ssvHiPur()};
		// working_points_["ssvHiEff"] = [](const Jet* jet) {return (bool) jet->ssvHiEff()};
		// working_points_["trkHiPur"] = [](const Jet* jet) {return (bool) jet->trkHiPur()};
		// working_points_["trkHiEff"] = [](const Jet* jet) {return (bool) jet->trkHiEff()};
		// // working_points_[] = [](const Jet* jet) {};

		naming_[TTNaming::RIGHT ] = "semilep_visible_right";
		naming_[TTNaming::RIGHT_THAD ] = "semilep_right_thad" 	 ;
		naming_[TTNaming::RIGHT_WHAD ] = "semilep_right_whad" 	 ;
		naming_[TTNaming::RIGHT_TLEP ] = "semilep_right_tlep" 	 ;
		naming_[TTNaming::WRONG ] = "semilep_wrong" 			 ;
		naming_[TTNaming::OTHER ] = 	"other"              ;

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
	void book_combo_plots(string folder){
		// book<TH1F>(folder, "nu_chisq"         , "", 520,  -2., 50.);
		// book<TH1F>(folder, "nu_discriminant"	, "", 110,  -1., 10.);
		// book<TH1F>(folder, "btag_discriminant", "", 300, -30.,  0.);
		book<TH1F>(folder, "mass_discriminant", "", 200,   0., 20.);
		// book<TH1F>(folder, "full_discriminant", "", 200, -10., 10.);
		
		// book<TH1F>(folder, "btag_value", "", 100, 0., 1.);
		book<TH2F>(folder, "Wmasshad_tmasshad", "", 500, 0., 500., 500, 0., 500);			
		// book<TH2F>(folder, "WtMass_special", "", 500, 0., 500., 500, 0., 500);			
		// book<TH2F>(folder, "WtMass_rphi", "", 500, 0., 500., 500,  0, 5);
	}

	void book_hyp_plots(string folder){
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		//book<TH1F>(folder, "nbjets"   , "", 50, 0., 50.);
		book<TH1F>(folder, "lep_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "had_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(#ell) (GeV)", 500, 0., 500.);
		// book<TH1F>(folder, "Wlep_mass", ";m_{W}(lep) (GeV)", 140, 0., 140.);
		// book<TH1F>(folder, "Wlep_char", ";charge_{W}(lep) (GeV)", 2, -1, 1);
		book<TH1F>(folder, "Whad_mass", ";m_{W}(had) (GeV)", 140, 0., 140.);
		book<TH1F>(folder, "Whad_DR"  , ";m_{W}(had) (GeV)", 100, 0., 10.);
		book<TH1F>(folder, "Whad_pt"  , ";m_{W}(had) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "jetCSV"   , "", 40, -20., 20.);
		// book<TH1F>(folder, "Whad_leading_DR", "", 100, 0., 10.);
		// book<TH1F>(folder, "Whad_sublead_DR", "", 100, 0., 10.);
		// book<TH2F>(folder, "Whad_lead_sub_DR", ";#DeltaR(leading jet, WHad); #DeltaR(subleading jet, WHad)", 100, 0., 10., 100, 0., 10.);
		// book<TH1F>(folder, "tlep_mass", ";m_{t}(lep) (GeV)", 300, 100., 400.);
		book<TH1F>(folder, "thad_mass", ";m_{t}(had) (GeV)", 300, 100., 400.);
		// book<TH1F>(folder, "hjet_pt"  , ";m_{t}(lep) (GeV)", 100, 0., 500.);
		// book<TH2F>(folder, "hjet_pts" , ";m_{t}(had) (GeV)", 100, 0., 500., 100, 0., 500.);
		// book<TH2F>(folder, "Whad_pt_lead_pt" , "", 100, 0., 500., 100, 0., 500.);
		// book<TH2F>(folder, "hjet_es"  , ";m_{t}(had) (GeV)", 100, 0., 500., 100, 0., 500.);
	}

	void book_jet_plots(string folder){
	}

  virtual void begin()
  {
    Logger::log().debug() << "ctag_eff::begin" << endl;
    outFile_.cd();
		vector<string> folders;

		//FIXME use folders as root dir, makes much more sense!
		//Would also be nice to have nothing instead of "all" for the others
		//semilep_visible_right becomes a SubdirectoryView in the plotter
		if(isTTbar_) folders = {"semilep_visible_right", "semilep_right_thad",
														"semilep_right_tlep", "semilep_right_whad", "semilep_wrong", "other"};
		else folders = {""};
		// string folders[] = {"all", "semilep_visible_right", "semilep_right_thad", 
		// 										"semilep_right_tlep", "semilep_right_whad", "semilep_wrong", "other"};
		string wjet_folders[] = {"leading", "subleading"};
		string tagging[] = {"lead_tagged", "sublead_tagged", "both_tagged", "both_untagged"};

    //TH1::AddDirectory(true);
		//if(isTTbar_) book_hyp_plots("gen");
		for(auto& genCategory : folders){			
			for(auto& sys : systematics_){
				string gcategory;
				if(!genCategory.empty()) gcategory  = genCategory +"/";
				string sys_name = systematics::shift_to_name.at(sys);
        Logger::log().debug() << "booking histos in " << gcategory+sys_name << endl;
				book<TH1F>(gcategory+sys_name, "nvtx_noweight", "", 100, 0., 100.);
				book<TH1F>(gcategory+sys_name, "nvtx", "", 100, 0., 100.);
				book<TH1F>(gcategory+sys_name, "weight", "", 100, 0., 20.);
				string criterion = cut_ordering_;
				for(auto& wp_item : working_points_) {
					string working_point = wp_item.first;
          if(working_point == "notag" && sys != systematics::SysShifts::NOSYS) continue;
					string base;
					if(!genCategory.empty()) base  = genCategory +"/";
					base += sys_name + "/" + criterion + "/" + working_point;
					book<TH1F>(base, "dummy"  , "", 1, 0., 500.);
					// book<TH2F>(base, "passing_jpt_massDiscriminant"  , "", 100, 0., 500., 200, 0., 20.);
					// book<TH2F>(base, "failing_jpt_massDiscriminant"  , "", 100, 0., 500., 200, 0., 20.);
					// book<TH2F>(base, "csvL_csvS"  , ";m_{t}(had) (GeV)", 100, 0., 1., 100, 0., 1.);
					for(auto& tag : tagging){
						if(working_point == "notag" && tag != "both_untagged") continue;
						string folder = base + "/" + tag;
						Logger::log().debug() << "creating plots for folder: " << folder << std::endl;
						
						book_combo_plots(folder);
						book_hyp_plots(folder);

						for(auto& w_folder : wjet_folders){
							string current = folder + "/" + w_folder;
							book<TH1F>(current, "eta"	,"eta"	, 100, -5, 5);
							book<TH1F>(current, "pt" 	,"pt" 	, 100, 0 , 500);
							book<TH1F>(current, "phi"	,"phi"	, 100, -4, 4);
							book<TH1F>(current, "pflav_smart","pflav", 55, -27.5, 27.5);
							book<TH1F>(current, "abs_pflav_smart","pflav", 28, -0.5, 27.5);
							book<TH1F>(current, "hadronflav","hflav", 28, -0.5, 27.5);
							book<TH1F>(current, "energy", ";E_{jet} (GeV)", 100, 0., 500.);
							book<TH1F>(current, "ncharged", "", 50, 0., 50.);						
							book<TH1F>(current, "nneutral", "", 50, 0., 50.);						
							book<TH1F>(current, "ntotal"  , "", 50, 0., 50.);						
						}
						Logger::log().debug() << "...Done" << endl;
						//if(genCategory == "semilep_visible_right"){
						// book<TH1F>(folder, "nu_DR"    , "#DeltaR between gen and reco #nu;#DeltaR;counts", 140, 0., 7.);
						// book<TH1F>(folder, "nu_DE"    , "#DeltaE between gen and reco #nu;#DeltaE (GeV);counts", 250, -250, 250.);
						//}
					}//for(auto& tag : tagging)
				}//for(auto& wp_item : working_points_)
			}//for(auto& sys : systematics){
		}//for(auto& genCategory : folders)

		/*string folder = "preselection";
		Logger::log().debug() << "creating plots for folder: " << folder << std::endl;
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		book<TH1F>(folder, "nleadjets", "", 10, 0., 10.);
		//book<TH1F>(folder, "nbjets"   , "", 50, 0., 50.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(#ell) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "lep_char" , ";charge_{W}(lep) (GeV)", 2, -1, 1);
		book<TH1F>(folder, "tlep_char", "Tight leptons charge;charge_{#ell}", 2, -1, 1);
		book<TH1F>(folder, "mu_char"  , "Tight muons charge;charge_{#mu}", 2, -1, 1);
		book<TH1F>(folder, "el_char"  , "Tight electrons charge;charge_{e}", 2, -1, 1);
		book<TH1F>(folder, "cjets_pt" , "Charm jet pt;p_{T} (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "cbjets_pt", "B-Tagged charm jet pt;p_{T} (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "nperms", "number of possible permutations", 40, -0.5, 39.5);
		book<TH1F>(folder, "matchable", "do we have the right permutation among the selected ones", 5, -0.5, 4.5);
		book<TH1F>(folder, "matching", "do we match the gen?", 5, -0.5, 4.5);
		book<TH1F>(folder, "decay", "which top decay is?", 5, -0.5, 4.5);*/
  }

	// void fill_gen_info(string folder, Permutation &reco, GenTTBar &gen, float weight){
	// 	auto dir = histos_.find(folder);
  //   if(dir == histos_.end()) {
	// 		Logger::log().error() << "histogram folder: " << folder <<
  //       " does not exists!" << endl;
  //     return;
  //   }
	// 	// dir->second["nu_DR"].fill(reco_wlep.second->DeltaR(*gen_wlep.second)  , weight);
	// 	// dir->second["nu_DE"].fill(reco_wlep.second->E() - gen_wlep.second->E(), weight);
	// }

	void fill_jet_info(string folder, const IDJet* jet, float weight){//, const Genparticle* genp=0){
		auto dir = histos_.find(folder);
		dir->second["eta"	 ].fill(jet->Eta(), weight);
		dir->second["pt" 	 ].fill(jet->Pt() , weight);
		dir->second["phi"	 ].fill(jet->Phi(), weight);
		//dir->second["pflav"].fill(jet->partonFlavour());
		dir->second["hflav"].fill(jet->hadronFlavour(), weight);
		dir->second["energy"].fill(jet->E(), weight);

		dir->second["ncharged"].fill(jet->numChargedHadrons(), weight);
		dir->second["nneutral"].fill(jet->numNeutralHadrons(), weight);
		dir->second["ntotal"  ].fill(jet->numChargedHadrons()+jet->numNeutralHadrons(), weight);

		//const IDJet* idjet = (const IDJet*) jet;
		dir->second["pflav_smart"].fill(jet->flavor(), weight);
		dir->second["abs_pflav_smart"].fill(fabs(jet->flavor()), weight);

    dir->second["hadronflav"].fill(fabs(jet->hadronFlavour()), weight);
	}

	void fill_discriminator_info(string folder, const Permutation &hyp, float weight){
		//Logger::log().warning() << "filling " << hyp.btag_discriminant << endl;
		auto dir = histos_.find(folder);
		// dir->second["full_discriminant"].fill(hyp.Prob(), weight );
		// dir->second["nu_chisq"         ].fill(hyp.NuChisq() );
		// dir->second["nu_discriminant"	 ].fill(hyp.NuDiscr()	);
		// dir->second["btag_discriminant"].fill(hyp.BDiscr() );
		dir->second["mass_discriminant"].fill(hyp.MassDiscr(), weight);

		double whad_mass = hyp.WHad().M();
		double thad_mass = hyp.THad().M();
		dir->second["Wmasshad_tmasshad"].fill(whad_mass, thad_mass, weight);
		// dir->second["WtMass_special"   ].fill(whad_mass + thad_mass, thad_mass - whad_mass);
		double r = Sqrt(pow(whad_mass,2) + pow(thad_mass,2));
		double phi = ATan(thad_mass/whad_mass);
		// dir->second["WtMass_rphi"      ].fill(r, phi, weight);
	}

	void fill(string folder, Permutation &hyp, size_t njets, float weight){//, TTbarHypothesis *genHyp=0) {
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "fill: histogram folder: " << folder <<
				" does not exists!" << endl;
			throw 40;
		}

		fill_discriminator_info(folder, hyp, weight);
		
		//dir->second["Wlep_char"].fill(charge);
		dir->second["njets"    ].fill(njets , weight);
		//dir->second["nbjets"   ].fill(nbjets, weight);
		dir->second["lep_b_pt" ].fill(hyp.BLep()->Pt(), weight);
		dir->second["had_b_pt" ].fill(hyp.BHad()->Pt(), weight);
		dir->second["lep_pt"   ].fill(hyp.L()->Pt(), weight);
		// dir->second["Wlep_mass"].fill(hyp.WLep().M());
		dir->second["Whad_mass"].fill(hyp.WHad().M(), weight);
		dir->second["Whad_DR"  ].fill(hyp.WJa()->DeltaR(*hyp.WJb()), weight);
		// dir->second["tlep_mass"].fill(hyp.TLep().M());
		dir->second["thad_mass"].fill(hyp.THad().M() , weight);
		// dir->second["hjet_pt"  ].fill(hyp.WJa()->Pt(), weight);
		// dir->second["hjet_pt"  ].fill(hyp.WJb()->Pt(), weight);
		dir->second["jetCSV"   ].fill(hyp.WJa()->csvIncl(), weight);
		dir->second["jetCSV"   ].fill(hyp.WJb()->csvIncl(), weight);
		float lead = (hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJa()->Pt() : hyp.WJb()->Pt();
		float sub  = (hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJb()->Pt() : hyp.WJa()->Pt();
		// dir->second["hjet_pts" ].fill(lead, sub, weight);//*/

		//if(folder == "gen") return;
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();

		// dir->second["Whad_leading_DR" ].fill(hyp.WHad().DeltaR(*leading));
		// dir->second["Whad_sublead_DR" ].fill(hyp.WHad().DeltaR(*subleading));
		// dir->second["Whad_lead_sub_DR"].fill(hyp.WHad().DeltaR(*leading), hyp.WHad().DeltaR(*subleading));
		// dir->second["Whad_pt_lead_pt" ].fill(hyp.WHad().Pt(), leading->Pt());		
		// dir->second["hjet_es"  ].fill(leading->E(), subleading->E(), weight);

		//FILL JET INFORMATION
		fill_jet_info(folder + "/leading", leading, weight);
		fill_jet_info(folder + "/subleading", subleading, weight);
	}

	void fill_other_jet_plots(string folder, Permutation &hyp, std::function<bool(const IDJet*)>& fcn, float weight) {
    vector<IDJet*>& jets = object_selector_.clean_jets();
		//folder += sys_name + "/" + criterion + "/" + working_point;
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "fill_other_jet_plots: histogram folder: " << folder <<
				" does not exists!" << endl;
			throw 40;
		}
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();
		// dir->second["csvL_csvS" ].fill(leading->csvIncl(), subleading->csvIncl(), weight);

		set<IDJet*> hypjets;
		hypjets.insert(hyp.WJa());
		hypjets.insert(hyp.WJb());
		hypjets.insert(hyp.BLep());
		hypjets.insert(hyp.BHad());
		// for(IDJet* jet : jets) {
		// 	if(hypjets.find(jet) == hypjets.end()){
		// 		if(fcn(jet))
		// 			dir->second["passing_jpt_massDiscriminant"].fill(jet->Pt(), hyp.MassDiscr(), weight);
		// 		else
		// 			dir->second["failing_jpt_massDiscriminant"].fill(jet->Pt(), hyp.MassDiscr(), weight);
		// 	}
		// }
	}

	string get_wjet_category(Permutation &hyp, std::function<bool(const IDJet*)>& fcn) {
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();
		bool lead_tag = fcn(leading   );
		bool sub_tag  = fcn(subleading);
		if(lead_tag && sub_tag) return "/both_tagged";
		else if(lead_tag) return "/lead_tagged";
		else if(sub_tag)  return "/sublead_tagged";
		return "/both_untagged";
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
    //LeptonType lep_type = (object_selector_.tight_muons().size() > 0) ? LeptonType::MUON : LeptonType::ELECTRON;

		// if(shift == systematics::SysShifts::NOSYS){
		// 	histos_["preselection"]["nleadjets"].fill(leading_jets.size(), weight);
		// 	histos_["preselection"]["njets"   ].fill(selected_jets.size(), weight);
		// 	histos_["preselection"]["nbjets"  ].fill(bjets.size(), weight);
		// 	histos_["preselection"]["lep_pt"  ].fill(lepton->Pt(), weight);
		// 	histos_["preselection"]["lep_char"].fill(lep_charge*0.5, weight);
		// }

    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;
    if(!isData_) evt_weight_ *= btag_sf_.scale_factor(permutator_, shift);
    tracker_.track("perm preselection");

    string sys_name = systematics::shift_to_name.at(shift);
		histos_[sys_name]["nvtx"].fill(event.vertexs().size(), evt_weight_);
		histos_[sys_name]["nvtx_noweight"].fill(event.vertexs().size());
				
    //Find best permutation
    bool go_on = true;
    Permutation best_permutation;
    size_t ncycles_=0;
    while(go_on) {
      ncycles_++;
      Permutation test_perm = permutator_.next(go_on);
      if(go_on) {
        test_perm.Solve(solvers_[systematics::SysShifts::NOSYS]);
        if(ordering_fcn_(test_perm, best_permutation)){
          best_permutation = test_perm;
        }
      }
    }
    if(!best_permutation.IsComplete() || best_permutation.Prob() > 1E9) return; //FIXME, is right??? best_permutation.Prob() > 1E9
    tracker_.track("best perm");

    size_t capped_jets_size = permutator_.capped_jets().size();
    best_permutation.permutating_jets(capped_jets_size);

    //find mc weight
    if(!isData_) { 
      if(object_selector_.tight_muons().size() == 1)
        evt_weight_ *= muon_sf_.get_sf(best_permutation.L()->Pt(), best_permutation.L()->Eta());
      if(object_selector_.medium_electrons().size() == 1)
        evt_weight_ *= electron_sf_.get_sf(best_permutation.L()->Pt(), best_permutation.L()->Eta());
		}

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
		histos_[sys_dir]["weight"].fill(evt_weight_);
		//Logger::log().debug() << "\n\nShift: " << shift << " name: " << root_dir << endl;
		if(!isTTbar_){
			for(auto& wpoint : working_points_){
        if(wpoint.first == "notag" && shift != systematics::SysShifts::NOSYS) continue;
				string jet_category = get_wjet_category(best_permutation, wpoint.second);
				string folder = sys_dir+"/"+cut_ordering_+"/"+wpoint.first;
				fill_other_jet_plots(folder, best_permutation, wpoint.second, evt_weight_);
				folder += jet_category;
				//Logger::log().debug() << "filling: " << folder << endl;
				fill(folder, best_permutation, object_selector_.clean_jets().size(), evt_weight_);
			}
		} else {
			//define which subdir we fall in
			TTNaming dir_id = get_ttdir_name(matched_perm, best_permutation);
			string ttsubdir = naming_.at(dir_id);

			//Logger::log().debug() << ttsubdir << " " << dir_id << " " << TTNaming::RIGHT_WHAD << endl;
			//if(dir_id <= TTNaming::RIGHT_WHAD) gen = &gen_hyp;

			for(auto& wpoint : working_points_){
        if(wpoint.first == "notag" && shift != systematics::SysShifts::NOSYS) continue;
				string jet_category = get_wjet_category(best_permutation, wpoint.second);
				string folder = ttsubdir+"/"+sys_dir+"/"+cut_ordering_+"/"+wpoint.first;
        
				fill_other_jet_plots(folder, best_permutation, wpoint.second, evt_weight_);
				folder += jet_category;

				//Logger::log().debug() << "filling: " << folder << endl;
				fill(folder, best_permutation, object_selector_.clean_jets().size(), evt_weight_);
			}
		} //if(!isTTbar_)
	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    Logger::log().debug() << "ctag_eff::analyze" << endl;
		unsigned long evt_idx = 0;
    URStreamer event(tree_);

    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();

		tracker_.deactivate();
    while(event.next())
    {
			// if(evt_idx % 1000 == 0) Logger::log().warning() << "Beginning event: " <<
      //                           evt_idx << endl;
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			if(evt_idx % 10000 == 1) Logger::log().debug() << "Beginning event " << evt_idx << endl;

			//long and time consuming
			if(isTTbar_){
				genp_selector_.select(event);			
			}

			for(auto shift : systematics_){
        evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);
				if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(shift, event);
				if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
			}

		} //while(event.next())
    cout << "Processed " << evt_idx << " events" << endl;
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

    parser.addCfgParameter<std::string>("permutations", "ordering", "ID to be applied");
    //("permutations.ordering", opts::value<string>()->default_value("mass_discriminant"));
	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<ctag_eff> test;
  return test.run();
}
