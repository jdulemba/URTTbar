#include <iostream>
#include "AnalyzerBase.h"
#include "URStreamer.h"
#include "URDriver.h"
#include "genSelection.h"
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


class btag_efficiency : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
	enum SysShift {NOSYS, JES_UP, JES_DOWN};
private:
	//histograms and helpers
	map<string, map<string, RObject> > histos_;
	map<TTNaming, string> naming_;
	CutFlowTracker tracker_;
	TTBarSolver solver_;
	TRandom3 randomizer_;// = TRandom3(98765);

	//switches
	bool isData_, isTTbar_, run_jes_;

	//CUTS
	float cut_lmuon_pt_, cut_lmuon_eta_, cut_tmuon_pt_, cut_tmuon_eta_, cut_whad_mass_, cut_thad_mass_;
	float cut_lelectron_pt_, cut_lelectron_eta_, cut_telectron_pt_, cut_telectron_eta_;
	float cut_jet_pt_, cut_jet_eta_, cut_leadB_pt_, cut_subB_pt_, cut_leadJet_pt_, cut_subJet_pt_;
	size_t cut_njets_permutated_;
	IDJet::BTag cut_leadB_ID_, cut_subB_ID_;
	string cut_ordering_;

	float DR_match_to_gen_ = 0.3;

	std::function<bool(const Permutation &, const Permutation &)> ordering_fcn_;
	map<string, std::function<bool(const IDJet*)> > working_points_;
	map<SysShift, string> sys_names_ = {
		{NOSYS, "nosys"}, 
		{JES_UP  , "jes_up"}, 
		{JES_DOWN, "jes_down"}
	};
	vector<SysShift> systematics_;

	//Scale factors
	TH1D *electron_sf_, *muon_sf_;
	TH1F *pu_sf_;
  // Nothing by default
  // Add your private variables/methods here

public:
  btag_efficiency(const std::string output_filename):
    AnalyzerBase("btag_efficiency", output_filename), 
		tracker_(),
		working_points_(),
		solver_(),
		randomizer_() 
	{
		TUUID id;  
		randomizer_.SetSeed(id.Hash());   

		solver_.Init("Prob_parton.root");
		opts::variables_map &values = URParser::instance().values();
		string output_file = values["output"].as<std::string>();
		size_t slash = output_file.rfind("/") + 1;
		string basename(output_file, slash, output_file.size() - slash);
		isData_  = boost::starts_with(basename, "data");
		isTTbar_ = boost::starts_with(basename, "ttJets");
		run_jes_ = !isData_;
		
		//A sort of FileInPath would be MUCH better!
		if(isData_){
			electron_sf_=0;
			muon_sf_=0;
			pu_sf_=0;			
		} else {
			TH1::AddDirectory(false);
			TFile lsf("Lepton_SF.root");
			electron_sf_ = (TH1D*) lsf.Get("Scale_ElTOT_Pt")->Clone();
			muon_sf_ = (TH1D*) lsf.Get("Scale_MuTOT_Pt")->Clone();
			TFile pusf("vtxs.root");
			pu_sf_ = (TH1F*) pusf.Get("vtx_reweigt")->Clone();
			TH1::AddDirectory(true);
		}

    // SYSTEMATCIS
		systematics_ = {NOSYS};
		if(run_jes_) {
			systematics_.push_back(JES_UP  );
			systematics_.push_back(JES_DOWN);
		}

		//SET CUTS FROM CFG
		cut_ordering_ = URParser::instance().getCfgPar<string>("permutations.ordering");
		cut_leadB_ID_ = available_bjet_id[
			URParser::instance().getCfgPar<string>("permutations.leadB_ID")];
		cut_subB_ID_  = available_bjet_id[
			URParser::instance().getCfgPar<string>("permutations.subB_ID")];
		ordering_fcn_ = available_ordering[cut_ordering_];

		cut_lmuon_pt_  		 		= URParser::instance().getCfgPar<float>("looseMuons.pt");
		cut_lmuon_eta_ 		 		= URParser::instance().getCfgPar<float>("looseMuons.abseta");
		cut_tmuon_pt_  		 		= URParser::instance().getCfgPar<float>("tightMuons.pt");
		cut_tmuon_eta_ 		 		= URParser::instance().getCfgPar<float>("tightMuons.abseta");
		cut_lelectron_pt_  		= URParser::instance().getCfgPar<float>("looseElectrons.pt");
		cut_lelectron_eta_ 		= URParser::instance().getCfgPar<float>("looseElectrons.abseta");
		cut_telectron_pt_  		= URParser::instance().getCfgPar<float>("tightElectrons.pt");
		cut_telectron_eta_ 		= URParser::instance().getCfgPar<float>("tightElectrons.abseta");
		cut_jet_pt_   		 		= URParser::instance().getCfgPar<float>("jets.pt");
		cut_jet_eta_  		 		= URParser::instance().getCfgPar<float>("jets.abseta");
		cut_leadB_pt_ 		 		= URParser::instance().getCfgPar<float>("permutations.leadB_pt");
		cut_subB_pt_  		 		= URParser::instance().getCfgPar<float>("permutations.subB_pt");
		cut_leadJet_pt_ 	 		= URParser::instance().getCfgPar<float>("permutations.leadJet_pt");
		cut_subJet_pt_  	 		= URParser::instance().getCfgPar<float>("permutations.subJet_pt");
		cut_njets_permutated_ = URParser::instance().getCfgPar<size_t>("permutations.njets_permutated");
		cut_whad_mass_        = URParser::instance().getCfgPar<float>("permutations.whad_mass");
		cut_thad_mass_        = URParser::instance().getCfgPar<float>("permutations.thad_mass");

		naming_[TTNaming::RIGHT ] = "semilep_visible_right";
		naming_[TTNaming::RIGHT_THAD ] = "semilep_right_thad" 	 ;
		naming_[TTNaming::RIGHT_WHAD ] = "semilep_right_whad" 	 ;
		naming_[TTNaming::RIGHT_TLEP ] = "semilep_right_tlep" 	 ;
		naming_[TTNaming::WRONG ] = "semilep_wrong" 			 ;
		naming_[TTNaming::OTHER ] = 	"other"              ;

		working_points_["notag"]     = [](const IDJet* jet) {return false;};
		working_points_["csvLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVLOOSE);};
		//working_points_["csvTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVTIGHT);};
		//working_points_["csvMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVMEDIUM);};
		// working_points_["rndm10"]    = [](const IDJet* jet) {return jet->rndm() < 0.1;};
		// working_points_["ssvHiPur"] = [](const Jet* jet) {return (bool) jet->ssvHiPur()};
		// working_points_["ssvHiEff"] = [](const Jet* jet) {return (bool) jet->ssvHiEff()};
		// working_points_["trkHiPur"] = [](const Jet* jet) {return (bool) jet->trkHiPur()};
		// working_points_["trkHiEff"] = [](const Jet* jet) {return (bool) jet->trkHiEff()};
		// // working_points_[] = [](const Jet* jet) {};
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
		book<TH1F>(folder, "full_discriminant", "", 200, -10., 10.);
		
		// book<TH1F>(folder, "btag_value", "", 100, 0., 1.);
		book<TH2F>(folder, "Wmasshad_tmasshad", "", 500, 0., 500., 500, 0., 500);			
		// book<TH2F>(folder, "WtMass_special", "", 500, 0., 500., 500, 0., 500);			
		book<TH2F>(folder, "WtMass_rphi", "", 500, 0., 500., 500,  0, 5);
	}

	void book_hyp_plots(string folder){
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		book<TH1F>(folder, "nbjets"   , "", 50, 0., 50.);
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
		book<TH1F>(folder, "hjet_pt"  , ";m_{t}(lep) (GeV)", 100, 0., 500.);
		book<TH2F>(folder, "hjet_pts" , ";m_{t}(had) (GeV)", 100, 0., 500., 100, 0., 500.);
		// book<TH2F>(folder, "Whad_pt_lead_pt" , "", 100, 0., 500., 100, 0., 500.);
		book<TH2F>(folder, "hjet_es"  , ";m_{t}(had) (GeV)", 100, 0., 500., 100, 0., 500.);
	}

	void book_jet_plots(string folder){
	}

  virtual void begin()
  {
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

		//if(isTTbar_) book_hyp_plots("gen");
		for(auto& genCategory : folders){			
			for(auto& sys : systematics_){
				string gcategory;
				if(!genCategory.empty()) gcategory  = genCategory +"/";
				string sys_name = sys_names_[sys];
				book<TH1F>(gcategory+sys_name, "nvtx_noweight", "", 100, 0., 100.);
				book<TH1F>(gcategory+sys_name, "nvtx", "", 100, 0., 100.);
				book<TH1F>(gcategory+sys_name, "weight", "", 100, 0., 20.);
				string criterion = cut_ordering_;
				for(auto& wp_item : working_points_) {
					string working_point = wp_item.first;
					string base;
					if(!genCategory.empty()) base  = genCategory +"/";
					base += sys_name + "/" + criterion + "/" + working_point;
					book<TH2F>(base, "passing_jpt_massDiscriminant"  , "", 100, 0., 500., 200, 0., 20.);
					book<TH2F>(base, "failing_jpt_massDiscriminant"  , "", 100, 0., 500., 200, 0., 20.);
					book<TH2F>(base, "csvL_csvS"  , ";m_{t}(had) (GeV)", 100, 0., 1., 100, 0., 1.);
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
							book<TH1F>(current, "hflav","hflav", 55, -27.5, 27.5);
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

		string folder = "preselection";
		Logger::log().debug() << "creating plots for folder: " << folder << std::endl;
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		book<TH1F>(folder, "nleadjets", "", 10, 0., 10.);
		book<TH1F>(folder, "nbjets"   , "", 50, 0., 50.);
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
		book<TH1F>(folder, "decay", "which top decay is?", 5, -0.5, 4.5);
  }

	void fill_gen_info(string folder, const TTbarHypothesis &reco, const TTbarHypothesis &gen, float weight){
		auto dir = histos_.find(folder);
    if(dir == histos_.end()) {
			Logger::log().error() << "histogram folder: " << folder <<
        " does not exists!" << endl;
      return;
    }
		const WHypothesis &reco_wlep = reco.wplus.isLeptonic ? reco.wplus : reco.wminus;
    const WHypothesis &gen_wlep  = reco.wplus.isLeptonic ? gen.wplus : gen.wminus;
		// dir->second["nu_DR"].fill(reco_wlep.second->DeltaR(*gen_wlep.second)  , weight);
		// dir->second["nu_DE"].fill(reco_wlep.second->E() - gen_wlep.second->E(), weight);
	}

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
	}

	void fill_discriminator_info(string folder, const Permutation &hyp, float weight){
		//Logger::log().warning() << "filling " << hyp.btag_discriminant << endl;
		auto dir = histos_.find(folder);
		dir->second["full_discriminant"].fill(hyp.Prob(), weight );
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
		dir->second["WtMass_rphi"      ].fill(r, phi, weight);
	}

	void fill(string folder, Permutation &hyp, size_t njets, size_t nbjets, float weight){//, TTbarHypothesis *genHyp=0) {
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "histogram folder: " << folder <<
				" does not exists!" << endl;
			return;
		}

		fill_discriminator_info(folder, hyp, weight);
		
		//dir->second["Wlep_char"].fill(charge);
		dir->second["njets"    ].fill(njets , weight);
		dir->second["nbjets"   ].fill(nbjets, weight);
		dir->second["lep_b_pt" ].fill(hyp.BLep()->Pt(), weight);
		dir->second["had_b_pt" ].fill(hyp.BHad()->Pt(), weight);
		dir->second["lep_pt"   ].fill(hyp.L()->Pt(), weight);
		// dir->second["Wlep_mass"].fill(hyp.WLep().M());
		dir->second["Whad_mass"].fill(hyp.WHad().M(), weight);
		dir->second["Whad_DR"  ].fill(hyp.WJa()->DeltaR(*hyp.WJb()), weight);
		// dir->second["tlep_mass"].fill(hyp.TLep().M());
		dir->second["thad_mass"].fill(hyp.THad().M() , weight);
		dir->second["hjet_pt"  ].fill(hyp.WJa()->Pt(), weight);
		dir->second["hjet_pt"  ].fill(hyp.WJb()->Pt(), weight);
		dir->second["jetCSV"   ].fill(hyp.WJa()->csvIncl(), weight);
		dir->second["jetCSV"   ].fill(hyp.WJb()->csvIncl(), weight);
		float lead = (hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJa()->Pt() : hyp.WJb()->Pt();
		float sub  = (hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJb()->Pt() : hyp.WJa()->Pt();
		dir->second["hjet_pts" ].fill(lead, sub, weight);//*/

		//if(folder == "gen") return;
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();

		// dir->second["Whad_leading_DR" ].fill(hyp.WHad().DeltaR(*leading));
		// dir->second["Whad_sublead_DR" ].fill(hyp.WHad().DeltaR(*subleading));
		// dir->second["Whad_lead_sub_DR"].fill(hyp.WHad().DeltaR(*leading), hyp.WHad().DeltaR(*subleading));
		// dir->second["Whad_pt_lead_pt" ].fill(hyp.WHad().Pt(), leading->Pt());		
		dir->second["hjet_es"  ].fill(leading->E(), subleading->E(), weight);

		//FILL JET INFORMATION
		fill_jet_info(folder + "/leading", leading, weight);
		fill_jet_info(folder + "/subleading", subleading, weight);
	}

	void fill_other_jet_plots(string folder, Permutation &hyp, vector<IDJet*>& jets, std::function<bool(const IDJet*)>& fcn, float weight) {
		//folder += sys_name + "/" + criterion + "/" + working_point;
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "histogram folder: " << folder <<
				" does not exists!" << endl;
			return;
		}
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();
		dir->second["csvL_csvS" ].fill(leading->csvIncl(), subleading->csvIncl(), weight);

		set<IDJet*> hypjets;
		hypjets.insert(hyp.WJa());
		hypjets.insert(hyp.WJb());
		hypjets.insert(hyp.BLep());
		hypjets.insert(hyp.BHad());
		for(IDJet* jet : jets) {
			if(hypjets.find(jet) == hypjets.end()){
				if(fcn(jet))
					dir->second["passing_jpt_massDiscriminant"].fill(jet->Pt(), hyp.MassDiscr(), weight);
				else
					dir->second["failing_jpt_massDiscriminant"].fill(jet->Pt(), hyp.MassDiscr(), weight);
			}
		}
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

	TTNaming get_ttdir_name(Permutation &gen, Permutation &reco, DecayType decay_type){
		// Logger::log().warning() << "fcn called: " << gen.decay << " (SEMI = "<< DecayType::SEMILEP <<")" << endl;
		if(decay_type == DecayType::SEMILEP){
			//if(!gen.hasMissingProng()){
			// Logger::log().debug() << "--------------------------------------" << endl;
			// Logger::log().debug() << "RECO: wh " << reco.whad()->first << " " << reco.whad()->second << " wl " <<
			// 	reco.wlep()->first << " b " << reco.b << " bb " << reco.bbar << endl;
			// Logger::log().debug() << "GEN: wh " << gen.whad()->first << " " << gen.whad()->second << " wl " <<
			// 	gen.wlep()->first << " b " << gen.b << " bb " << gen.bbar << endl;
				//}
			TTNaming name;
			if(gen.IsComplete() && reco.IsCorrect(gen)) name = TTNaming::RIGHT;//"semilep_visible_right";
			else if(gen.IsTHadComplete() && reco.IsTHadCorrect(gen)) name = TTNaming::RIGHT_THAD; //"semilep_right_thad";
			else if(gen.IsWHadComplete() && reco.IsWHadCorrect(gen)) name = TTNaming::RIGHT_WHAD; //"semilep_right_whad";
			else if(gen.IsTLepComplete() && reco.IsTLepCorrect(gen)) name = TTNaming::RIGHT_TLEP; //"semilep_right_tlep";
			else name = TTNaming::WRONG; //"semilep_wrong";

			/*bool whad_matches = reco.IsWHadCorrect(gen);
			bool b_matches    = reco.IsBHadCorrect(gen);
			bool lepb_matches = reco.IsBLepCorrect(gen);
			bool lep_matches  = reco.L() == gen.L();

			TTNaming xcheck;
			if(lepb_matches && lep_matches && whad_matches && b_matches) xcheck = TTNaming::RIGHT;//"semilep_visible_right";
			else if(whad_matches && b_matches) xcheck =  TTNaming::RIGHT_THAD; //"semilep_right_thad";
			else if(whad_matches) xcheck =  TTNaming::RIGHT_WHAD; //"semilep_right_whad";
			else if(lep_matches && lepb_matches) xcheck =  TTNaming::RIGHT_TLEP; //"semilep_right_tlep";
			else xcheck = TTNaming::WRONG; //"semilep_wrong";

			if(name == TTNaming::WRONG){
				Logger::log().error() << naming_[name] << " " << naming_[xcheck] << " differ." << endl;
				Logger::log().error() << "Gen completeness: " << endl <<
					"   L: " << gen.L() << " Bl: " << gen.BLep() << " Bh: " << gen.BHad() <<
					" Wj: " << gen.WJa() << ", " << gen.WJb() << endl << 
					"		gen.IsTHadComplete(): " << gen.IsTHadComplete() << endl <<
					"		gen.IsWHadComplete(): " << gen.IsWHadComplete() << endl <<
					"		gen.IsTLepComplete(): " << gen.IsTLepComplete() << endl <<
					"Gen matching: " << endl <<
					"   L: " << reco.L() << " Bl: " << reco.BLep() << " Bh: " << reco.BHad() <<
          " Wj: " << reco.WJa() << ", " << reco.WJb() << endl <<
					"   reco.IsTHadCorrect(gen): " << reco.IsTHadCorrect(gen)	<< endl <<
					"   reco.IsWHadCorrect(gen): " << reco.IsWHadCorrect(gen)	<< endl <<
					"   reco.IsBHadCorrect(gen): " << reco.IsBHadCorrect(gen)	<< endl <<
					"   reco.IsTLepCorrect(gen): " << reco.IsTLepCorrect(gen)	<< endl <<
					"   reco.IsBLepCorrect(gen): " << reco.IsBLepCorrect(gen)	<< endl <<
					"   reco.L() == gen.L()    : " << (reco.L() == gen.L())  	<< endl;
				char c;
				cin >> c;
				
			}//*/
			return name;
		}
		else {
			return TTNaming::OTHER; //"other";
		}
	}
		

	void process_evt(SysShift shift, URStreamer &event, TTbarHypothesis &gen_hyp){
		tracker_.track("start");
		float weight = 1.;

		bool mu_trg, el_trg;
		if(isData_){
			mu_trg = (event.trigger().HLT_IsoMu24_eta2p1() == 1);
			el_trg = (event.trigger().HLT_Ele27_eta2p1_WPLoose_Gsf() == 1);
		}	else {
			mu_trg = (event.trigger().HLT_IsoMu27() == 1);
			el_trg = (event.trigger().HLT_Ele27_WP85_Gsf() == 1);
		}
		if(!(mu_trg || el_trg)) return;
		tracker_.track("trigger");

		string sys_name = sys_names_[shift];
		if(!isData_) weight *= pu_sf_->GetBinContent(
			pu_sf_->FindFixBin(event.vertexs().size())
			);
		histos_[sys_name]["nvtx"].fill(event.vertexs().size(), weight);
		histos_[sys_name]["nvtx_noweight"].fill(event.vertexs().size());
		
		const double rho = event.rho().value();
		const vector<Muon>& muons = event.muons();
		list<  IDMuon> keep_muons;
		vector<IDMuon*> tight_muons; tight_muons.reserve(muons.size());
		vector<IDMuon*> loose_muons; loose_muons.reserve(muons.size());
		for(vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon){
			IDMuon mu(*muon, rho);
			if(mu.ID(IDMuon::LOOSE_12) && mu.Pt() > cut_lmuon_pt_ && abs(mu.Eta()) < cut_lmuon_eta_){
				keep_muons.push_back(mu);
				loose_muons.push_back(&keep_muons.back());
				if(mu.ID(IDMuon::TIGHT_12) && mu.Pt() > cut_tmuon_pt_ && abs(mu.Eta()) < cut_tmuon_eta_){
					if(!isData_) weight *= muon_sf_->GetBinContent(
						muon_sf_->FindFixBin(mu.Pt())
						);
					histos_["preselection"]["tlep_char"].fill(muon->charge()*0.5, weight);
					histos_["preselection"]["mu_char"].fill(muon->charge()*0.5  , weight);
					tight_muons.push_back(&keep_muons.back());
				}
			}
		}
		//Apply vetoes
		tracker_.track("selections");
		if(loose_muons.size() > 1) return;
		tracker_.track("muon veto");
			
		const vector<Electron>& electrons = event.electrons();
		list<  IDElectron> keep_electrons;
		vector<IDElectron*> tight_electrons; tight_electrons.reserve(electrons.size());
		vector<IDElectron*> loose_electrons; loose_electrons.reserve(electrons.size());
		for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron){
			IDElectron el(*electron, rho);
			if(el.ID(IDElectron::LOOSE_12) && el.Pt() > cut_lelectron_pt_ && abs(el.Eta()) < cut_lelectron_eta_){
				keep_electrons.push_back(el);
				loose_electrons.push_back(&keep_electrons.back());
				if(el.ID(IDElectron::MEDIUM_12) && el.Pt() > cut_telectron_pt_ && abs(el.Eta()) < cut_telectron_eta_){
					if(!isData_) weight *= electron_sf_->GetBinContent(
						electron_sf_->FindFixBin(el.Pt())
						);
					histos_["preselection"]["tlep_char"].fill(electron->charge()*0.5, weight);
					histos_["preselection"]["el_char"].fill(electron->charge()*0.5  , weight);
					tight_electrons.push_back(&keep_electrons.back());
				}
			}
		}
		if(loose_electrons.size() > 1) return;
		tracker_.track("electron veto");

		if(tight_muons.size() + tight_electrons.size() != 1) return;
		//remove mismatch in trgger and selected object 
		if(mu_trg && tight_electrons.size() == 1) return;
		if(el_trg && tight_muons.size() == 1) return;
		tracker_.track("One lepton");
		TLorentzVector* lepton = tight_muons.size() ? (TLorentzVector*) tight_muons[0] : (TLorentzVector*) tight_electrons[0];
		if(weight == 0){
			Logger::log().debug() << "PU*lepton weight: " << weight << endl;
			Logger::log().debug() << "# vtx: "<< event.vertexs().size() << ", muons "<< tight_muons.size() <<", electrons: " << tight_electrons.size()  << ", lepton.pt: " << lepton->Pt() << endl;
		}

		const vector<Jet>& jets = event.jets();
		list<IDJet> keep_jets;
		vector<IDJet*> selected_jets; selected_jets.reserve(jets.size());
		for(vector<Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet){
			IDJet idjet(*jet, randomizer_.Rndm());
			if(shift == JES_UP)	idjet.SetPxPyPzE(jet->Px()*1.02, jet->Py()*1.02, jet->Pz()*1.02, jet->E()*1.02);
			else if(shift == JES_DOWN) idjet.SetPxPyPzE(jet->Px()*0.98, jet->Py()*0.98, jet->Pz()*0.98, jet->E()*0.98);

			if(!(idjet.Pt() > cut_jet_pt_ && Abs(idjet.Eta()) < cut_jet_eta_)) continue;

			bool overlaps = false;
			for(auto mu : loose_muons){
				if(mu->DeltaR(idjet) < 0.4){
					overlaps = true;
					break;
				}
			}

			for(auto el : loose_electrons){
				if(el->DeltaR(idjet) < 0.4){
					overlaps = true;
					break;
				}
			}

			if(overlaps) continue;
				
			keep_jets.push_back(idjet);
			selected_jets.push_back(&keep_jets.back());
		}
		if(selected_jets.size() < 4) return;
		tracker_.track("4 jets");

		vector<Jet*> bjets; bjets.reserve(selected_jets.size());
		for(auto jet : selected_jets){
			if(jet->BTagId(cut_leadB_ID_) || jet->BTagId(cut_subB_ID_)) bjets.push_back(jet);
		}
			
		IDMet met(event.METs()[0]);

		if(bjets.size() < 2) return;
		tracker_.track("2 bjets");

		int lep_charge = tight_muons.size() ? tight_muons[0]->charge() : tight_electrons[0]->charge();
		//sort jets by pt
		sort(selected_jets.begin(), selected_jets.end(), [](const Jet* one, const Jet* two) {return  one->Pt() > two->Pt();});
		size_t pick = selected_jets.size() < cut_njets_permutated_  ? selected_jets.size() : cut_njets_permutated_; //avoid vector out of bounds
		vector<IDJet*> leading_jets(selected_jets.begin(), selected_jets.begin()+pick); //pick first njets_permutated
		if(leading_jets.size() > 20){
			Logger::log().debug() << "#Jets: " << leading_jets.size() << endl;
		}

		//check that the pt theshold are satisfied
		float pt_thresholds[4] = {cut_leadB_pt_, cut_subB_pt_, cut_leadJet_pt_, cut_subJet_pt_};
		sort(pt_thresholds, pt_thresholds+4);
		for(int i=0; i<4; i++){
			if(leading_jets[i]->Pt() < pt_thresholds[i]) return;
		}
		tracker_.track("jets pt cut");

		if(shift == NOSYS){
			histos_["preselection"]["nleadjets"].fill(leading_jets.size(), weight);
			histos_["preselection"]["njets"   ].fill(selected_jets.size(), weight);
			histos_["preselection"]["nbjets"  ].fill(bjets.size(), weight);
			histos_["preselection"]["lep_pt"  ].fill(lepton->Pt(), weight);
			histos_["preselection"]["lep_char"].fill(lep_charge*0.5, weight);
		}
		
		//sort leading jets to make next_permutation to work properly
		sort(leading_jets.begin(), leading_jets.end());

		Permutation best_permutation;
		double figure_of_merit = std::numeric_limits<double>::max();
		size_t ncombos=0;
		size_t n_passing_combos=0;
		for(auto bjet1 : leading_jets){
			for(auto bjet2 : leading_jets){
				if(bjet1 == bjet2) continue;
				for(auto wjet1 : leading_jets){
					if(wjet1 == bjet1) continue;
					if(wjet1 == bjet2) continue;
					for(auto wjet2 : leading_jets){
						if(wjet2 == bjet1) continue;
						if(wjet2 == bjet2) continue;
						if(wjet2 == wjet1) continue;
						if(wjet1->Pt() < wjet2->Pt()) continue;
						ncombos++;
						tracker_.track("Combination start");
						auto lead_b = (bjet1->Pt() > bjet2->Pt()) ? bjet1 : bjet2;
						auto sub_b  = (bjet1->Pt() > bjet2->Pt()) ? bjet2 : bjet1;
						
						if(!lead_b->BTagId(cut_leadB_ID_)) continue;
						if(!sub_b->BTagId(cut_subB_ID_)) continue;
						tracker_.track("Bjet Tag cuts");
						
						if(lead_b->Pt() < cut_leadB_pt_) continue;
						if(sub_b->Pt() < cut_subB_pt_) continue;
						tracker_.track("Bjet PT cuts");
						
						//remove W had ambiguity			
						if(wjet1->Pt() < cut_leadJet_pt_) continue;
						if(wjet2->Pt() < cut_subJet_pt_) continue;
						tracker_.track("Ambiguity cut");
						
						Permutation hyp(
							wjet1, 
							wjet2, 
							bjet1, 
							bjet2, 
							lepton, 
							&met
							);
						hyp.Solve(solver_);
						
						//Hypothesis selection (rough cuts so far)
						if(fabs(hyp.WHad().M() - 80) > cut_whad_mass_) continue;
						tracker_.track("W mass cut");
						if(fabs(hyp.THad().M() - 173) > cut_thad_mass_) continue;
						tracker_.track("Top mass cut");
						if(hyp.NuChisq() < 0) continue;
						tracker_.track("Neutrino chisquare cut");
						//if(hyp.mass_discriminant > 4) continue;
						
						n_passing_combos++;
						if(ordering_fcn_(hyp, best_permutation)){
							best_permutation = hyp;
						}
					} //for(auto wjet2 : leading_jets){
				} //for(auto wjet1 : leading_jets){
			} //for(auto bjet2 : leading_jets){
		} //for(auto bjet1 : leading_jets){


		histos_["preselection"]["nperms"].fill(n_passing_combos, weight);
		if(best_permutation.WJa() == 0) return;
		tracker_.track("one ttbar hypothesis");

		/*
			FETCH TTBAR GEN HYPOTHESIS (IF IT IS TT MC)
		*/
		Permutation matched_hyp;
		if(isTTbar_){
			//fill("gen", gen_hyp, selected_jets.size(), bjets.size());
			//WARNING! The matched object MAY change
			matched_hyp = match_to_gen(
				gen_hyp, 
				//leading_jets,
				selected_jets, 
				tight_electrons,
				tight_muons,
				DR_match_to_gen_
				);
			matched_hyp.SetMET(&met);

			size_t code=0;
			if(gen_hyp.decay == DecayType::SEMILEP){
				if(matched_hyp.IsComplete()) code=4;
				else if(matched_hyp.IsTHadComplete()) code=3;
				else if(matched_hyp.IsWHadComplete()) code=2;
				else if(matched_hyp.IsTLepComplete()) code=1;
				else code=0;
				histos_["preselection"]["matching"].fill(code, weight);
			}
			histos_["preselection"]["decay"].fill(gen_hyp.decay, weight);
		}

		string sys_dir = sys_names_[shift];
		histos_[sys_dir]["weight"].fill(weight);
		//Logger::log().debug() << "\n\nShift: " << shift << " name: " << root_dir << endl;
		if(!isTTbar_){
			for(auto& wpoint : working_points_){
				string jet_category = get_wjet_category(best_permutation, wpoint.second);
				string folder = sys_dir+"/"+cut_ordering_+"/"+wpoint.first;
				fill_other_jet_plots(folder, best_permutation, selected_jets, wpoint.second, weight);
				folder += jet_category;
				//Logger::log().debug() << "filling: " << folder << endl;
				fill(folder, best_permutation, selected_jets.size(), bjets.size(), weight);
			}
		} else {
			//define which subdir we fall in
			TTNaming dir_id = get_ttdir_name(matched_hyp, best_permutation, gen_hyp.decay);
			string ttsubdir = naming_[dir_id];
			TTbarHypothesis *gen = 0;

			//Logger::log().debug() << ttsubdir << " " << dir_id << " " << TTNaming::RIGHT_WHAD << endl;
			if(dir_id <= TTNaming::RIGHT_WHAD) gen = &gen_hyp;

			for(auto& wpoint : working_points_){
				string jet_category = get_wjet_category(best_permutation, wpoint.second);
				string folder = ttsubdir+"/"+sys_dir+"/"+cut_ordering_+"/"+wpoint.first;

				fill_other_jet_plots(folder, best_permutation, selected_jets, wpoint.second, weight);
				folder += jet_category;

				//Logger::log().debug() << "filling: " << folder << endl;
				fill(folder, best_permutation, selected_jets.size(), bjets.size(), weight);
			}
		} //if(!isTTbar_)
	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
		unsigned long evt_idx = 0;
    URStreamer event(tree_);

    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();

		tracker_.deactivate();
    while(event.next())
    {
			// if(evt_idx % 100 == 0) Logger::log().warning() << "Beginning event: " <<
			// 												 evt_idx << endl;
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			if(evt_idx % 100 == 0 || (evt_idx > 6900 && evt_idx < 7000)) Logger::log().debug() << "Beginning event " << evt_idx << endl;

			//long and time consuming
			TTbarHypothesis gen_hyp;
			if(isTTbar_){
				gen_hyp = SelectGenParticles(event);			
			}

			for(auto shift : systematics_){
				if(shift == NOSYS) tracker_.activate();
				//Logger::log().debug() << "processing: " << shift << endl;
				process_evt(shift, event, gen_hyp);
				if(shift == NOSYS) tracker_.deactivate();
			}

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
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");

		opts::options_description &loose_mu_cfg = parser.optionGroup("looseMuons", "Loose muons selections", URParser::CFG);
		loose_mu_cfg.add_options()
			("looseMuons.pt", opts::value<float>()->default_value(15.))
			("looseMuons.abseta", opts::value<float>()->default_value(8.));
			
		opts::options_description &tight_mu_cfg = parser.optionGroup("tightMuons", "Tight muons selections", URParser::CFG);
		tight_mu_cfg.add_options()
      ("tightMuons.pt", opts::value<float>()->default_value(30.))
      ("tightMuons.abseta", opts::value<float>()->default_value(8.));

		opts::options_description &loose_el_cfg = parser.optionGroup("looseElectrons", "Loose electrons selections", URParser::CFG);
		loose_el_cfg.add_options()
      ("looseElectrons.pt", opts::value<float>()->default_value(15.))
      ("looseElectrons.abseta", opts::value<float>()->default_value(8.));

		opts::options_description &tight_el_cfg = parser.optionGroup("tightElectrons", "Tight electrons selections", URParser::CFG);
		tight_el_cfg.add_options()
      ("tightElectrons.pt", opts::value<float>()->default_value(30.))
      ("tightElectrons.abseta", opts::value<float>()->default_value(8.));

		opts::options_description &jets_cfg = parser.optionGroup("jets", "Jets selections", URParser::CFG);
		tight_el_cfg.add_options()
      ("jets.pt", opts::value<float>()->default_value(20.))
      ("jets.abseta", opts::value<float>()->default_value(2.4));

		opts::options_description &permutations_cfg = parser.optionGroup("permutations", "permutations selections", URParser::CFG);
		permutations_cfg.add_options()
			("permutations.njets_permutated", opts::value<size_t>()->default_value(400)) //All of them
			("permutations.leadB_pt", opts::value<float>()->default_value(30.))
			("permutations.leadB_ID", opts::value<string>()->default_value("csvMedium"))
			("permutations.subB_pt", opts::value<float>()->default_value(30.))
			("permutations.subB_ID", opts::value<string>()->default_value("csvMedium"))
			("permutations.leadJet_pt", opts::value<float>()->default_value(20.))
			("permutations.subJet_pt", opts::value<float>()->default_value(20.))
			("permutations.whad_mass", opts::value<float>()->default_value(40.))
			("permutations.thad_mass", opts::value<float>()->default_value(50.))
			("permutations.ordering", opts::value<string>()->default_value("mass_discriminant"));
	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<btag_efficiency> test;
  return test.run();
}
