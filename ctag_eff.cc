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
#include "TTGenParticleSelector.h"
#include "TTPermutator.h"
#include "TTGenMatcher.h"
#include "MCWeightProducer.h"
#include "BTagSFProducer.h"
#include "LeptonSF.h"
#include "DataFile.h"
#include "PDFuncertainty.h"

using namespace TMath;

static map<string, IDJet::BTag> available_bjet_id = {
	{"none"     , IDJet::BTag::NONE},
	{"csvMedium", IDJet::BTag::CSVMEDIUM},
	{"csvLoose" , IDJet::BTag::CSVLOOSE},
	{"csvTight" , IDJet::BTag::CSVTIGHT},

	{"ctagMedium", IDJet::BTag::CTAGMEDIUM},
	{"ctagLoose" , IDJet::BTag::CTAGLOOSE},
	{"ctagTight" , IDJet::BTag::CTAGTIGHT},
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
	bool isData_, isTTbar_, pdfs_;

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
  PDFuncertainty pdf_uncs_;

	vector<systematics::SysShifts> systematics_;
	string cut_ordering_;
	std::function<bool(const Permutation &, const Permutation &)> ordering_fcn_;
	map<string, std::function<bool(const IDJet*)> > working_points_;
	map<string, BTagSFProducer> wp_SFs_;  

	//Scale factors
	LeptonSF electron_sf_, muon_sf_;

  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE;
  IDJet::BTag cut_loose_b_=IDJet::BTag::NONE;

  unsigned long evt_idx_ = 0;

public:
  std::vector<systematics::SysShifts> get_systematics(std::string outname) {
    std::string sample = systematics::get_sample(outname);
    std::vector<systematics::SysShifts> full_sys = {
      systematics::SysShifts::NOSYS, systematics::SysShifts::JES_UP, systematics::SysShifts::JES_DW,
      systematics::SysShifts::JER_UP, systematics::SysShifts::JER_DW,
      systematics::SysShifts::PU_UP, systematics::SysShifts::PU_DW,
      systematics::SysShifts::BTAG_UP, systematics::SysShifts::BTAG_DW, //Systematic on B-Jet selection
      systematics::SysShifts::BTAG_B_UP, systematics::SysShifts::BTAG_B_DW, //systematics on B/C-Tagging
      systematics::SysShifts::BTAG_L_UP, systematics::SysShifts::BTAG_L_DW, //systematics on B/C-Tagging
      systematics::SysShifts::BTAG_C_UP, systematics::SysShifts::BTAG_C_DW, //systematics on B/C-Tagging
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
    solver_(),
    mc_weights_(),
    evt_weight_(1.),
    electron_sf_("electron_sf", false),
    muon_sf_("muon_sf"),
		randomizer_(),
    btag_sf_("best_permutation.tightb", "best_permutation.looseb", 0.5, -1, -1),
    pdf_uncs_(249)
	{
    btag_sf_.ignore_partial_shifts();
    cut_tight_b_ = btag_sf_.tight_cut();
    cut_loose_b_ = btag_sf_.loose_cut();
    // cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << 
    //   " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;
    
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
		isData_  = boost::starts_with(sample, "data");
		isTTbar_ = boost::starts_with(sample, "ttJets");
    Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
    systematics_ = get_systematics(output_file);    
    int nosys = values["nosys"].as<int>();
    int nopdf = values["nopdf"].as<int>();
    if(nosys == 1) {
      systematics_ = {systematics::SysShifts::NOSYS};
      Logger::log().info() << "DISABLING SYSTEMATICS!" << endl;
    }
    if(isData_) {
      if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
      else object_selector_.lepton_type(1);
    }
    pdfs_ = isTTbar_ && (systematics_.size() > 1) && (nopdf == 0);

    if(!isData_) mc_weights_.init(sample);

    //Init solver
    string filename = "prob_ttJets.root";
    Logger::log().debug() << "solver file: " << filename << endl;
    TFile probfile(DataFile(filename).path().c_str());
    // for(auto shift : systematics_) {
    TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
    //   if(!td) td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());      
    solver_.Init(td, false, true, true);
    // }

		//SET CUTS FROM CFG
		cut_ordering_ = parser.getCfgPar<string>("permutations.ordering");
		ordering_fcn_ = available_ordering[cut_ordering_];

		working_points_["notag"]     = [](const IDJet* jet) {return false;};
		working_points_["csvLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVLOOSE);};
		working_points_["csvTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVTIGHT);};
		working_points_["csvMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVMEDIUM);};
		working_points_["ctagLoose"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGLOOSE);};
		working_points_["ctagTight"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGTIGHT);};
		working_points_["ctagMedium"] = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGMEDIUM);};

    //Get appropriate SFs for the probe working points
    DataFile csv_sfs(parser.getCfgPar<string>("general.csv_sffile"));
    DataFile wjet_efficiency(parser.getCfgPar<string>("general.wjets_efficiencies"));
    DataFile ctag_sfs(parser.getCfgPar<string>("general.ctag_sffile"));
    wp_SFs_["csvLoose" ] = BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVLOOSE, IDJet::BTag::NONE, 0.5, -1, -1); 
    wp_SFs_["csvTight" ] = BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVTIGHT, IDJet::BTag::NONE, 0.5, -1, -1); 
    wp_SFs_["csvMedium"] = BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVMEDIUM, IDJet::BTag::NONE, 0.5, -1, -1); 
    wp_SFs_["ctagLoose" ] = BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGLOOSE , IDJet::BTag::NONE, 0.5, -1, 0.5); 
    wp_SFs_["ctagTight" ] = BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGTIGHT , IDJet::BTag::NONE, 0.5, -1, 0.5); 
    wp_SFs_["ctagMedium"] = BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGMEDIUM, IDJet::BTag::NONE, 0.5, -1, 0.5); 

    wp_SFs_["csvLoose" ].ignore_general_shifts();
    wp_SFs_["csvTight" ].ignore_general_shifts();
    wp_SFs_["csvMedium"].ignore_general_shifts();
    wp_SFs_["ctagLoose" ].ignore_general_shifts(); 
    wp_SFs_["ctagTight" ].ignore_general_shifts(); 
    wp_SFs_["ctagMedium"].ignore_general_shifts(); 
    // // working_points_[] = [](const Jet* jet) {};

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

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
	void book_combo_plots(string folder){
		book<TH1F>(folder, "mass_discriminant", "", 20,   0., 20.);
    
		book<TH2F>(folder, "flav_map", ";leading;subleading", 3, 0., 3., 3, 0., 3.);

		book<TH1F>(folder, "tmasshad", "", 100, 0., 500);			
		book<TH1F>(folder, "Wmasshad", "", 100, 0., 500);			
	}

	void fill_combo_plots(string folder, const Permutation &hyp){
		auto dir = histos_.find(folder);
		dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
    float lflav = -1;
    switch(hyp.WJa()->hadronFlavour()) {
    case 5: lflav = 2.5; break;
    case 4: lflav = 1.5; break;
    default: lflav = 0.5; break;      
    }
    float sflav = -1;
    switch(hyp.WJb()->hadronFlavour()) {
    case 5:  sflav = 2.5; break;
    case 4:  sflav = 1.5; break;
    default: sflav = 0.5; break;      
    }
		dir->second["flav_map"].fill(lflav, sflav, evt_weight_);

		double whad_mass = hyp.WHad().M();
		double thad_mass = hyp.THad().M();
		dir->second["Wmasshad"].fill(whad_mass, evt_weight_);
		dir->second["tmasshad"].fill(thad_mass, evt_weight_);
	}

  void book_pdf_plots(string folder) {
    pdf_uncs_.book_replicas(folder, "mass_discriminant", 20,   0., 20.);
  }
  void fill_pdf_plots(string folder, const Permutation &hyp, URStreamer& streamer) {
    pdf_uncs_.fill_replicas(folder, "mass_discriminant", hyp.MassDiscr(), evt_weight_, streamer);
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
		book<TH1F>(folder, "jets_CvsL" , "", 55, -1., 1.1);
		book<TH1F>(folder, "jets_CvsB" , "", 55, -1., 1.1);
		book<TH1F>(folder, "jets_CSV"  , "", 55,  0., 1.1);

		book<TH1F>(folder, "jets_hflav_CvsL_B" , "", 55, -1., 1.1);
		book<TH1F>(folder, "jets_hflav_CvsL_C" , "", 55, -1., 1.1);
		book<TH1F>(folder, "jets_hflav_CvsL_S" , "", 55, -1., 1.1);
		book<TH1F>(folder, "jets_hflav_CvsL_L" , "", 55, -1., 1.1);
  }

  void fill_presel_plots(string folder, URStreamer &event) {
    auto dir = histos_.find(folder);
		dir->second["nvtx"].fill(event.vertexs().size(), evt_weight_);
		dir->second["nvtx_noweight"].fill(event.vertexs().size());
		dir->second["rho"].fill(event.rho().value(), evt_weight_);
    dir->second["weight"].fill(evt_weight_);
		dir->second["lep_pt"].fill(object_selector_.lepton()->Pt(), evt_weight_);
		dir->second["lep_eta"].fill(object_selector_.lepton()->Eta(), evt_weight_);
    for(IDJet* jet : object_selector_.clean_jets()) {
      dir->second["jets_pt"].fill(jet->Pt(), evt_weight_);
      dir->second["jets_eta"].fill(jet->Eta(), evt_weight_);
      dir->second["jets_CvsL"].fill(jet->CvsLtag(), evt_weight_);
      dir->second["jets_CvsB"].fill(jet->CvsBtag(), evt_weight_);
      dir->second["jets_CSV" ].fill(jet->csvIncl(), evt_weight_);

      int hflav = fabs(jet->hadronFlavour());
      int pflav = fabs(jet->partonFlavour());
      string hstr;
      if(hflav == 5) hstr="B";
      else if(hflav == 4) hstr="C";
      else if(pflav == 3) hstr="S";
      else hstr="L"; 
      dir->second["jets_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);
    }
		dir->second["njets" ].fill(object_selector_.clean_jets().size(), evt_weight_);
  }

	void book_notag_plots(string folder){
		book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
		book<TH1F>(folder, "lep_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "had_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "lep_pt"   , ";p_{T}(#ell) (GeV)", 500, 0., 500.);
		book<TH1F>(folder, "Whad_mass", ";m_{W}(had) (GeV)", 28, 0., 140.);
		book<TH1F>(folder, "Whad_DR"  , ";m_{W}(had) (GeV)", 100, 0., 10.);
		book<TH1F>(folder, "Whad_pt"  , ";m_{W}(had) (GeV)", 100, 0., 500.);
		book<TH1F>(folder, "thad_mass", ";m_{t}(had) (GeV)", 60, 100., 400.);
		book<TH2F>(folder, "Whad_jet_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
		book<TH2F>(folder, "BJet_jet_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
		book<TH2F>(folder, "leadB_leadW_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
		book<TH2F>(folder, "subB_subW_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);

    //info by flavor    
		book<TH1F>(folder, "Wjets_hflav_CvsL_B" , "", 42, -1., 1.1);
		book<TH1F>(folder, "Wjets_hflav_CvsL_C" , "", 42, -1., 1.1);
		book<TH1F>(folder, "Wjets_hflav_CvsL_L" , "", 42, -1., 1.1);

		book<TH1F>(folder, "Wjets_hflav_CvsB_B" , "", 42, -1., 1.1);
		book<TH1F>(folder, "Wjets_hflav_CvsB_C" , "", 42, -1., 1.1);
		book<TH1F>(folder, "Wjets_hflav_CvsB_L" , "", 42, -1., 1.1);

		book<TH1F>(folder, "Wjets_hflav_jpt_C" , "", 100, 0., 500.);    
		book<TH1F>(folder, "Wjets_hflav_jpt_L" , "", 100, 0., 500.);    
		book<TH1F>(folder, "Wjets_hflavJP_jpt_C" , "", 20, 25., 100.);    
		book<TH1F>(folder, "Wjets_hflavJP_jpt_L" , "", 20, 25., 100.);    
		book<TH2F>(folder, "Wjets_hflav_jpt_jeta_C" , "", 20, 25., 100., 16, -2.5, 2.5);    
		book<TH2F>(folder, "Wjets_hflav_jpt_jeta_L" , "", 20, 25., 100., 16, -2.5, 2.5);    
    

		book<TH1F>(folder, "Wja_hflav_CvsL_B" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wja_hflav_CvsL_C" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wja_hflav_CvsL_L" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wjb_hflav_CvsL_B" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wjb_hflav_CvsL_C" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wjb_hflav_CvsL_L" , "", 55, -1., 1.1);

		book<TH1F>(folder, "Wja_CvsL" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wjb_CvsL" , "", 55, -1., 1.1);

		book<TH1F>(folder, "Wjets_CvsL" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Wjets_CvsB" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Bjets_CvsL" , "", 55, -1., 1.1);
		book<TH1F>(folder, "Bjets_CvsB" , "", 55, -1., 1.1);
		book<TH1F>(folder, "WjetCSV"   , "", 40, -20., 20.);
	}

  void fill_notag_plots(string folder, Permutation &hyp){
    auto dir = histos_.find(folder);
		dir->second["njets"    ].fill(object_selector_.clean_jets().size(), evt_weight_);
		dir->second["lep_b_pt" ].fill(hyp.BLep()->Pt(), evt_weight_);
		dir->second["had_b_pt" ].fill(hyp.BHad()->Pt(), evt_weight_);
		dir->second["lep_pt"   ].fill(hyp.L()->Pt(), evt_weight_);
		dir->second["Whad_mass"].fill(hyp.WHad().M(), evt_weight_);
		dir->second["Whad_DR"  ].fill(hyp.WJa()->DeltaR(*hyp.WJb()), evt_weight_);
		dir->second["thad_mass"].fill(hyp.THad().M() , evt_weight_);
		dir->second["WjetCSV"  ].fill(hyp.WJa()->csvIncl(), evt_weight_);
		dir->second["WjetCSV"  ].fill(hyp.WJb()->csvIncl(), evt_weight_);

		const IDJet *lb = (hyp.BHad()->Pt() > hyp.BLep()->Pt()) ? hyp.BHad() : hyp.BLep();
		const IDJet *sb = (hyp.BHad()->Pt() > hyp.BLep()->Pt()) ? hyp.BLep() : hyp.BHad();
		const IDJet *lj = (hyp.WJa()->Pt()  > hyp.WJb()->Pt()) ? hyp.WJa() : hyp.WJb();
		const IDJet *sj = (hyp.WJa()->Pt()  > hyp.WJb()->Pt()) ? hyp.WJb() : hyp.WJa();

		dir->second["Whad_jet_pts" ].fill(lj->Pt(), sj->Pt()  , evt_weight_);		
		dir->second["BJet_jet_pts" ].fill(lb->Pt(), sb->Pt()  , evt_weight_);		
		dir->second["leadB_leadW_pts"].fill(lb->Pt(), lj->Pt(), evt_weight_);
		dir->second["subB_subW_pts"  ].fill(sb->Pt(), sj->Pt(), evt_weight_);

    string wj="Wja";
    for(IDJet* jet : {hyp.WJa(), hyp.WJb()}){
      int hflav = fabs(jet->hadronFlavour());
      int pflav = fabs(jet->partonFlavour());
      string hstr;
      if(hflav == 5) hstr="B";
      else if(hflav == 4) {
        hstr="C";
        dir->second["Wjets_hflav_jpt_C"].fill(jet->Pt(), evt_weight_);
        dir->second["Wjets_hflavJP_jpt_C"].fill(jet->Pt(), evt_weight_);
        dir->second["Wjets_hflav_jpt_jeta_C"].fill(jet->Pt(), jet->Eta(), evt_weight_);
      }
      // else if(pflav == 3) {
      //   hstr="S";
      //   dir->second["Wjets_hflav_jpt_LS"].fill(jet->Pt(), evt_weight_);
      // }
      else {
        hstr="L"; 
        dir->second["Wjets_hflav_jpt_L"].fill(jet->Pt(), evt_weight_);
        dir->second["Wjets_hflavJP_jpt_L"].fill(jet->Pt(), evt_weight_);
        dir->second["Wjets_hflav_jpt_jeta_L"].fill(jet->Pt(), jet->Eta(), evt_weight_);
      }

      dir->second["Wjets_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);
      dir->second["Wjets_hflav_CvsB_"+hstr].fill(jet->CvsBtag(), evt_weight_);
      dir->second[wj+"_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);
    }

		dir->second["Wja_CvsL"].fill(hyp.WJa()->CvsLtag() , evt_weight_);
		dir->second["Wjb_CvsL"].fill(hyp.WJb()->CvsLtag() , evt_weight_);

		dir->second["Wjets_CvsL"].fill(hyp.WJa()->CvsLtag() , evt_weight_);
		dir->second["Wjets_CvsL"].fill(hyp.WJb()->CvsLtag() , evt_weight_);
		dir->second["Wjets_CvsB"].fill(hyp.WJa()->CvsBtag() , evt_weight_);
		dir->second["Wjets_CvsB"].fill(hyp.WJb()->CvsBtag() , evt_weight_);
		dir->second["Bjets_CvsL"].fill(hyp.BHad()->CvsLtag(), evt_weight_);
		dir->second["Bjets_CvsL"].fill(hyp.BLep()->CvsLtag(), evt_weight_);
		dir->second["Bjets_CvsB"].fill(hyp.BHad()->CvsBtag(), evt_weight_);
		dir->second["Bjets_CvsB"].fill(hyp.BLep()->CvsBtag(), evt_weight_);
  }

	void book_jet_plots(string folder){
    book<TH1F>(folder, "eta"	,"eta"	, 100, -5, 5);
    book<TH1F>(folder, "pt" 	,"pt" 	, 100, 0 , 500);
    book<TH1F>(folder, "phi"	,"phi"	, 100, -4, 4);
    book<TH1F>(folder, "pflav_smart","pflav", 55, -27.5, 27.5);
    book<TH1F>(folder, "abs_pflav_smart","pflav", 28, -0.5, 27.5);
    book<TH1F>(folder, "hadronflav","hflav", 28, -0.5, 27.5);
    book<TH1F>(folder, "energy", ";E_{jet} (GeV)", 100, 0., 500.);
    book<TH1F>(folder, "ncharged", "", 50, 0., 50.);						
    book<TH1F>(folder, "nneutral", "", 50, 0., 50.);						
    book<TH1F>(folder, "ntotal"  , "", 50, 0., 50.);						
	}

	void fill_jet_plots(string folder, const IDJet* jet){//, const Genparticle* genp=0){
		auto dir = histos_.find(folder);
		dir->second["eta"	 ].fill(jet->Eta(), evt_weight_);
		dir->second["pt" 	 ].fill(jet->Pt() , evt_weight_);
		dir->second["phi"	 ].fill(jet->Phi(), evt_weight_);
		dir->second["energy"].fill(jet->E(), evt_weight_);

		dir->second["ncharged"].fill(jet->numChargedHadrons(), evt_weight_);
		dir->second["nneutral"].fill(jet->numNeutralHadrons(), evt_weight_);
		dir->second["ntotal"  ].fill(jet->numChargedHadrons()+jet->numNeutralHadrons(), evt_weight_);

		dir->second["pflav_smart"].fill(jet->flavor(), evt_weight_);
		dir->second["abs_pflav_smart"].fill(fabs(jet->flavor()), evt_weight_);
    dir->second["hadronflav"].fill(fabs(jet->hadronFlavour()), evt_weight_);
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
        if(genCategory.empty() || genCategory == "semilep_visible_right") {
          book_presel_plots(gcategory+sys_name+"/preselection");
          book_combo_plots(gcategory+sys_name+"/permutations");
        }
				string criterion = cut_ordering_;
				for(auto& wp_item : working_points_) {
					string working_point = wp_item.first;
					string base;
					if(!genCategory.empty()) base  = genCategory +"/";
					base += sys_name + "/" + criterion + "/" + working_point;
					book<TH1F>(base, "dummy"  , "", 1, 0., 500.);
					for(auto& tag : tagging){
						if(working_point == "notag" && tag != "both_untagged") continue;
						string folder = base + "/" + tag;
						
						book_combo_plots(folder);
            if(pdfs_ && sys == systematics::SysShifts::NOSYS) book_pdf_plots(folder);
						if(working_point == "notag") book_notag_plots(folder);

						for(auto& w_folder : wjet_folders){
							string current = folder + "/" + w_folder;
              book_jet_plots(current);
						}
					}//for(auto& tag : tagging)
				}//for(auto& wp_item : working_points_)
			}//for(auto& sys : systematics){
		}//for(auto& genCategory : folders)
  }

	void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
		auto dir = histos_.find(folder);
		if(dir == histos_.end()) {
			Logger::log().error() << "fill: histogram folder: " << folder <<
				" does not exists!" << endl;
			throw 40;
		}

		fill_combo_plots(folder, hyp);
		
		const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
		const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();

		//FILL JET INFORMATION
		fill_jet_plots(folder + "/leading",       leading);
		fill_jet_plots(folder + "/subleading", subleading);
	}

	/*void fill_other_jet_plots(string folder, Permutation &hyp, std::function<bool(const IDJet*)>& fcn, float weight) {
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
    }*/

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

    //MC Weight for lepton selection
    if(!isData_) { 
      if(object_selector_.tight_muons().size() == 1)
        evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
      if(object_selector_.medium_electrons().size() == 1)
        evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
		}

    string sys_name = systematics::shift_to_name.at(shift);
    string presel_dir = sys_name;
    if(isTTbar_){
      presel_dir = naming_.at(TTNaming::RIGHT) + "/" + sys_name;
    }
    fill_presel_plots(presel_dir+"/preselection", event);

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
    
    //NOW, cut on btagging
    if(!(best_permutation.BHad()->BTagId(cut_tight_b_) || best_permutation.BLep()->BTagId(cut_tight_b_))) {return;}
    if(!(best_permutation.BHad()->BTagId(cut_loose_b_) && best_permutation.BLep()->BTagId(cut_loose_b_))) {return;}
    tracker_.track("perm btag");

    size_t capped_jets_size = permutator_.capped_jets().size();
    best_permutation.permutating_jets(capped_jets_size);

    //find mc weight for btag
    if(!isData_) evt_weight_ *= btag_sf_.scale_factor({best_permutation.BHad(), best_permutation.BLep()}, shift);

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
    if(isTTbar_) { 			
      TTNaming dir_id = get_ttdir_name(matched_perm, best_permutation);
      ttsubdir = naming_.at(dir_id) + "/";
    }

    for(auto& wpoint : working_points_){        
      auto wp_sf = wp_SFs_.find(wpoint.first);
      if(!isData_ && wp_sf != wp_SFs_.end()) evt_weight_ = weight*wp_sf->second.scale_factor({best_permutation.WJa(), best_permutation.WJb()}, shift);
      else evt_weight_ = weight;
      string jet_category = get_wjet_category(best_permutation, wpoint.second);
      string folder = ttsubdir+sys_dir+"/"+cut_ordering_+"/"+wpoint.first;
      //fill_other_jet_plots(folder, best_permutation, wpoint.second, evt_weight_);
      folder += jet_category;
      //Logger::log().debug() << "filling: " << folder << endl;
      if(wpoint.first == "notag") fill_notag_plots(folder, best_permutation);
      fill(folder, best_permutation);
      if(pdfs_ && shift == systematics::SysShifts::NOSYS) fill_pdf_plots(folder, best_permutation, event);
    }
	}

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    opts::variables_map &values = URParser::instance().values();
		int limit = values["limit"].as<int>();
		int skip  = values["skip"].as<int>();

    if(evt_idx_ >= limit) return;
    Logger::log().debug() << "ctag_eff::analyze" << endl;
    URStreamer event(tree_);

		tracker_.deactivate();
    while(event.next())
    {
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
      ("nopdf", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");

    parser.addCfgParameter<std::string>("general", "csv_sffile", "");
    parser.addCfgParameter<std::string>("general", "ctag_sffile", "");
    parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");

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
