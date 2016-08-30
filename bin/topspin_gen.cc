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
#include <unordered_map>
#include "Analyses/URTTbar/interface/TwoBodyDecay.h"
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "Analyses/URTTbar/interface/LHEParticle.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "Analyses/URTTbar/interface/TTKinFitter.h"

using namespace std;
typedef TLorentzVector tlv;
typedef TVector3 tv3;

class topspin_gen : public AnalyzerBase
{
private:
	//histograms and helpers
	CutFlowTracker tracker_;
  //histos
  unordered_map<string, unordered_map<string, RObject> > histos_;

	//switches
	bool isTTbar_, isSignal_;

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
  TTKinFitter kinfitter_;

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
    permutator_(),
    kinfitter_()
  {
    
    //tracker_.verbose(true);    
    //set tracker
    tracker_.use_weight(&evt_weight_);
    object_selector_.set_tracker(&tracker_);
    permutator_.set_tracker(&tracker_);

    URParser &parser = URParser::instance();
    opts::variables_map &values = parser.values();

    //find out which sample are we running on
		string output_file = values["output"].as<std::string>();
    //DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
		string sample = systematics::get_sample(output_file);
    isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
		isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal_;

    if(!isTTbar_) {
      Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
      throw 49;
    }

    int lhe = parser.getCfgPar<int>("general", "uselhe");
    if(lhe > 0) {
      Logger::log().info() << "Using LHE information" << endl;
      if(isSignal_) genp_selector_.setmode(TTGenParticleSelector::SelMode::MADGRAPHLHE);
      else genp_selector_.setmode(TTGenParticleSelector::SelMode::LHE);
    } else {
      Logger::log().info() << "Using normal matching" << endl;
      if(isSignal_) genp_selector_.setmode(TTGenParticleSelector::SelMode::FULLDEP);
      else genp_selector_.setmode(TTGenParticleSelector::SelMode::NORMAL);
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

  void book_kin_plots(string dir) {
    book<TH1F>(dir,  "npass", "", 1, 0., 1);
    book<TH1F>(dir,  "lpt", "", 100, 0., 300);
    book<TH1F>(dir, "blpt", "", 100, 0., 300);
    book<TH1F>(dir, "bhpt", "", 100, 0., 300);
    book<TH1F>(dir, "wupt", "", 100, 0., 300);
    book<TH1F>(dir, "wdpt", "", 100, 0., 300);
    
    book<TH1F>(dir,  "leta", "", 100, 0., 3.);
    book<TH1F>(dir, "bleta", "", 100, 0., 3.);
    book<TH1F>(dir, "bheta", "", 100, 0., 3.);
    book<TH1F>(dir, "wueta", "", 100, 0., 3.);
    book<TH1F>(dir, "wdeta", "", 100, 0., 3.);
    
    book<TH1F>(dir, "ttm",   "", 100, 300., 1000);
    book<TH1F>(dir, "thm",   "",  80, 130., 210);
    book<TH1F>(dir, "tlm",   "",  80, 130., 210);
    book<TH1F>(dir, "Whm",   "",  80,  40., 120);
    book<TH1F>(dir, "ttpt",  "", 100, 0., 800);
    book<TH1F>(dir, "tteta", "", 100, 0., 3.);
  }

  void fill_kin_plots(string dir, hyp::TTbar& ttbar) {
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }    
    hdir->second["npass"].fill(0.5);

    hdir->second["thm"].fill(ttbar.thad().M(), evt_weight_);
    hdir->second["tlm"].fill(ttbar.tlep().M(), evt_weight_);
    hdir->second["Whm"].fill(ttbar.thad().W().M(), evt_weight_);

    hdir->second[ "lpt"].fill(ttbar.tlep().W().l().Pt(), evt_weight_);
    hdir->second["blpt"].fill(ttbar.tlep().b().Pt(), evt_weight_);
    hdir->second["bhpt"].fill(ttbar.thad().b().Pt(), evt_weight_);
    hdir->second["wupt"].fill(ttbar.thad().W().up().Pt(), evt_weight_);
    hdir->second["wdpt"].fill(ttbar.thad().W().down().Pt(), evt_weight_);
    
    hdir->second[ "leta"].fill(fabs(ttbar.tlep().W().l().Eta()) , evt_weight_);
    hdir->second["bleta"].fill(fabs(ttbar.tlep().b().Eta())     , evt_weight_);
    hdir->second["bheta"].fill(fabs(ttbar.thad().b().Eta())     , evt_weight_);
    hdir->second["wueta"].fill(fabs(ttbar.thad().W().up().Eta()), evt_weight_);
    hdir->second["wdeta"].fill(fabs(ttbar.thad().W().down().Eta()), evt_weight_);

    hdir->second["ttm"  ].fill(ttbar.M(), evt_weight_);
    hdir->second["ttpt" ].fill(ttbar.Pt(), evt_weight_);
    hdir->second["tteta"].fill(ttbar.Eta(), evt_weight_);
  }

  void book_angular_plots(string dir) {
    // book<TH1F>(dir, "costheta_wl", "", 200, -1, 1);
    // book<TH1F>(dir, "costheta_wh", "", 200, -1, 1);

    book<TH1F>(dir, "costheta_top_BA", "", 200, -1, 1);
    book<TH1F>(dir, "costheta_top_CM", "", 200, -1, 1);

    book<TH1F>(dir, "cosgamma_bWtlep", "", 200, -1, 1);
    book<TH1F>(dir, "cosgamma_bWthad", "", 200, -1, 1);
    //cos alpha bW -- bWbar planes
    book<TH1F>(dir, "cosalpha_tbW_tbarbW", "", 200, -1, 1);

    //helicity and lab frames from TOP 14 023
    book<TH1F>(dir, "helframe_costheta_lep", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_costheta_nu" , "", 200, -1, 1);
    book<TH1F>(dir, "helframe_costheta_dtype", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_costheta_utype", "", 200, -1, 1);
    book<TH1F>(dir, "labframe_cosdeltaphi_lep_dtype", "", 200, -1, 1);
    book<TH1F>(dir, "labframe_deltaphi_lep_dtype", "", 200, 0, Pi());
    book<TH1F>(dir, "helframe_cosdelta_lep_dtype", "", 200, -1, 1);

    book<TH1F>(dir, "helframe_prodcosth_lep_dtype", "", 200, -1, 1);      
    book<TH1F>(dir, "labframe_cosdeltaphi_lep_utype", "", 200, -1, 1);
    book<TH1F>(dir, "labframe_deltaphi_lep_utype", "", 200, 0, Pi());
    book<TH1F>(dir, "helframe_cosdelta_lep_utype", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_prodcosth_lep_utype", "", 200, -1, 1);      

    book<TH1F>(dir, "labframe_cosdeltaphi_nu_dtype", "", 200, -1, 1);
    book<TH1F>(dir, "labframe_deltaphi_nu_dtype", "", 200, 0, Pi());
    book<TH1F>(dir, "helframe_cosdelta_nu_dtype", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_prodcosth_nu_dtype", "", 200, -1, 1);      
    book<TH1F>(dir, "labframe_deltaphi_nu_utype", "", 200, 0, Pi());
    book<TH1F>(dir, "labframe_cosdeltaphi_nu_utype", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_cosdelta_nu_utype", "", 200, -1, 1);
    book<TH1F>(dir, "helframe_prodcosth_nu_utype", "", 200, -1, 1);      

    book<TH1F>(dir, "wframe_costheta_lep", "", 200, -1, 1);
    book<TH1F>(dir, "wframe_costheta_nu" , "", 200, -1, 1);
    book<TH1F>(dir, "wframe_costheta_dtype", "", 200, -1, 1);
    book<TH1F>(dir, "wframe_costheta_utype", "", 200, -1, 1);

		//energies
    book<TH1F>(dir, "wframe_energy_dtype", "", 200, 0, 100);
    book<TH1F>(dir, "wframe_energy_utype", "", 200, 0, 100);
    book<TH1F>(dir, "tframe_energy_dtype", "", 200, 0, 100);
    book<TH1F>(dir, "tframe_energy_utype", "", 200, 0, 100);
    book<TH1F>(dir, "labframe_energy_dtype", "", 200, 0, 300);
    book<TH1F>(dir, "labframe_energy_utype", "", 200, 0, 300);

		//W angle
    book<TH1F>(dir, "tframe_wlep_costheta", "", 200, -1, 1);
    book<TH1F>(dir, "tframe_whad_costheta", "", 200, -1, 1);

		//b angles
    book<TH1F>(dir, "tframe_blep_costheta", "", 200, -1, 1);
    book<TH1F>(dir, "tframe_bhad_costheta", "", 200, -1, 1);

		//top angles
    book<TH1F>(dir, "ttframe_top_costheta", "", 200, -1, 1);
    book<TH1F>(dir, "ttframe_tba_costheta", "", 200, -1, 1);
  }

  void fill_angular_plots(string dir, hyp::TTbar& ttbar) {
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }    
    
    auto leptopcm = ttbar.tlep().to_CM();
    auto hadtopcm = ttbar.thad().to_CM();
    auto ttcm = ttbar.to_CM();
    double mass = ttbar.M();
    TVector3 beam(0,0,1);
    
    hdir->second["costheta_top_BA"].fill(beam.Dot(ttcm.top().unit3D()), evt_weight_);
    hdir->second["costheta_top_CM"].fill(ttbar.unit3D().Dot(ttcm.top().unit3D()), evt_weight_);

    hdir->second["cosalpha_tbW_tbarbW"].fill(ttcm.top().decay_plane().Dot(ttcm.tbar().decay_plane()), evt_weight_);
    hdir->second["cosgamma_bWtlep"].fill(leptopcm.W().decay_plane().Dot(leptopcm.b().unit3D()), evt_weight_);
    hdir->second["cosgamma_bWthad"].fill(hadtopcm.W().decay_plane().Dot(hadtopcm.b().unit3D()), evt_weight_);
    double ctheta_l = ttbar.tlep().W().decay_opening_cm();
    double ctheta_h = ttbar.thad().W().decay_opening_cm();
    //hdir->second["costheta_wl"].fill(ctheta_l);
    //hdir->second["costheta_wh"].fill(ctheta_h);
        
    //helicity and lab frames from TOP 14 023
    double cth_l = leptopcm.W().l().unit3D().Dot(ttcm.tlep().unit3D());
    double cth_n = leptopcm.W().nu().unit3D().Dot(ttcm.tlep().unit3D());
    double cth_d = hadtopcm.W().down().unit3D().Dot(ttcm.thad().unit3D());
    double cth_u = hadtopcm.W().up().unit3D().Dot(ttcm.thad().unit3D());

		auto lepwcm = ttbar.tlep().W().to_CM();
		auto hadwcm = ttbar.thad().W().to_CM();		
    hdir->second["wframe_costheta_lep"  ].fill(leptopcm.W().unit3D().Dot(lepwcm.l().unit3D()));
    hdir->second["wframe_costheta_nu"   ].fill(leptopcm.W().unit3D().Dot(lepwcm.nu().unit3D()));
    hdir->second["wframe_costheta_dtype"].fill(hadtopcm.W().unit3D().Dot(hadwcm.down().unit3D()));
    hdir->second["wframe_costheta_utype"].fill(hadtopcm.W().unit3D().Dot(hadwcm.up().unit3D()));

    if(cth_l != cth_l || cth_n != cth_n || cth_d != cth_d || cth_u != cth_u) {
      Logger::log().error() <<dir << " NANERROR! "<< cth_l << " " << cth_n << " " << cth_d << " " << cth_u << endl;
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
      throw 40;
    }
    hdir->second["helframe_costheta_lep"  ].fill(cth_l, evt_weight_);
    hdir->second["helframe_costheta_nu"   ].fill(cth_n, evt_weight_);
    hdir->second["helframe_costheta_dtype"].fill(cth_d, evt_weight_);
    hdir->second["helframe_costheta_utype"].fill(cth_u, evt_weight_);
    hdir->second["helframe_cosdelta_lep_dtype"].fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()), evt_weight_);
    hdir->second["helframe_cosdelta_lep_utype"].fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D()), evt_weight_);
    hdir->second["helframe_prodcosth_lep_dtype"].fill(cth_l*cth_d, evt_weight_);      
    hdir->second["helframe_prodcosth_lep_utype"].fill(cth_l*cth_u, evt_weight_);      

    hdir->second["helframe_cosdelta_nu_dtype"].fill(leptopcm.W().nu().unit3D().Dot(hadtopcm.W().down().unit3D()), evt_weight_);
    hdir->second["helframe_cosdelta_nu_utype"].fill(leptopcm.W().nu().unit3D().Dot(hadtopcm.W().up().unit3D())  , evt_weight_);
    hdir->second["helframe_prodcosth_nu_dtype"].fill(cth_n*cth_d, evt_weight_);      
    hdir->second["helframe_prodcosth_nu_utype"].fill(cth_n*cth_u, evt_weight_);      

    auto lep = ttbar.tlep().W().l();
    auto nu = ttbar.tlep().W().nu();
    auto up = ttbar.thad().W().up();
    auto dw = ttbar.thad().W().down();    
    hdir->second["labframe_cosdeltaphi_lep_dtype"].fill(Cos(lep.DeltaPhi(up)), evt_weight_);
    hdir->second["labframe_cosdeltaphi_lep_utype"].fill(Cos(lep.DeltaPhi(dw)), evt_weight_);
    hdir->second["labframe_deltaphi_lep_dtype"].fill(lep.DeltaPhi(up), evt_weight_);
    hdir->second["labframe_deltaphi_lep_utype"].fill(lep.DeltaPhi(dw), evt_weight_);

    hdir->second["labframe_cosdeltaphi_nu_dtype"].fill(Cos(nu.DeltaPhi(up)), evt_weight_);
    hdir->second["labframe_cosdeltaphi_nu_utype"].fill(Cos(nu.DeltaPhi(dw)), evt_weight_);
    hdir->second["labframe_deltaphi_nu_dtype"].fill(nu.DeltaPhi(up), evt_weight_);
    hdir->second["labframe_deltaphi_nu_utype"].fill(nu.DeltaPhi(dw), evt_weight_);


		//energies
    auto hadWcm = ttbar.thad().W().to_CM();
    hdir->second["wframe_energy_dtype"].fill(hadWcm.down().E(), evt_weight_);
    hdir->second["wframe_energy_utype"].fill(hadWcm.up().E(), evt_weight_);
    hdir->second["tframe_energy_dtype"].fill(hadtopcm.W().down().E(), evt_weight_);
    hdir->second["tframe_energy_utype"].fill(hadtopcm.W().up().E(), evt_weight_);
    hdir->second["labframe_energy_dtype"].fill(ttbar.thad().W().down().E(), evt_weight_);
    hdir->second["labframe_energy_utype"].fill(ttbar.thad().W().up().E(), evt_weight_);

		//W angle
    hdir->second["tframe_wlep_costheta"].fill(leptopcm.W().unit3D().Dot(ttcm.tlep().unit3D()), evt_weight_);
    hdir->second["tframe_whad_costheta"].fill(hadtopcm.W().unit3D().Dot(ttcm.thad().unit3D()), evt_weight_);

		//b angles
    hdir->second["tframe_blep_costheta"].fill(leptopcm.b().unit3D().Dot(ttcm.tlep().unit3D()), evt_weight_);
    hdir->second["tframe_bhad_costheta"].fill(hadtopcm.b().unit3D().Dot(ttcm.thad().unit3D()), evt_weight_);

		//top angles
    hdir->second["ttframe_top_costheta"].fill(ttbar.unit3D().Dot(ttcm.top().unit3D()), evt_weight_);
    hdir->second["ttframe_tba_costheta"].fill(ttbar.unit3D().Dot(ttcm.tbar().unit3D()), evt_weight_);
  }

  void book_eff_plots(string dir) {
    book<TH1F>(dir, "eff_mtt" , "", 56, 200, 3000);
  }

  void fill_eff_plots(hyp::TTbar& gen, string dir) {
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }
    hdir->second["eff_mtt"].fill(gen.M(), evt_weight_);
  }

  void book_res_plots(string dir) {
    //matrices
    book<TH2F>(dir, "mat_mtt" , ";gen;reco", 16, 200, 1000, 16, 200, 1000);    
    book<TH2F>(dir, "mat_pttt", ";gen;reco", 10, 0, 500, 10, 0, 500);

    //resolution (wrt mass) 
    for(string name : {
        "res_mtt", "res_pttt", "res_ptl", "res_pth", "res_mth", 
          "res_mt", "res_mtl", "res_deltatt", "res_mwh"}) {
      book<TH2F>(dir, name, ";rel_delta;mass", 120, -2, 2, 16, 200, 1000);
    }
    book<TH2F>(dir, "res_jete", ";rel_delta;mass", 120, -2, 2, 100, 0, 1000);
    book<TH2F>(dir, "res_deltaw_ptw", ";rel_delta;mass", 120, -2, 2, 100, 0, 1000);
    book<TH2F>(dir, "res_mw_ptw", ";rel_delta;mass", 120, -2, 2, 100, 0, 1000);
    book<TH2F>(dir, "res_mth_ptth", ";rel_delta;mass", 120, -2, 2, 100, 0, 1000);
  }

  inline double res(double reco, double gen) {return (reco-gen)/gen;}

  void fill_res_plots(string dir, hyp::TTbar& ttbar, hyp::TTbar* gen) {
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }    

    //matrices
    hdir->second["mat_mtt" ].fill(gen->M(), ttbar.M(), evt_weight_);
    hdir->second["mat_pttt"].fill(gen->Pt(), ttbar.Pt(), evt_weight_);

    //resolution    
    hdir->second["res_mtt" ].fill(res(ttbar.M(), gen->M()), gen->M(), evt_weight_);
    hdir->second["res_pttt"].fill(res(ttbar.Pt(), gen->Pt()), gen->M(), evt_weight_);
    hdir->second["res_ptl" ].fill(res(ttbar.tlep().Vect().Mag(), gen->tlep().Vect().Mag()), gen->M(), evt_weight_);
    hdir->second["res_pth" ].fill(res(ttbar.thad().Vect().Mag(), gen->thad().Vect().Mag()), gen->M(), evt_weight_);
    hdir->second["res_mth" ].fill(res(ttbar.thad().M(), gen->thad().M()), gen->M(), evt_weight_);
    hdir->second["res_mt"  ].fill(res(ttbar.thad().M(), 173.21), gen->M(), evt_weight_);
    hdir->second["res_mwh" ].fill(res(ttbar.thad().W().M(), gen->thad().W().M()), gen->M(), evt_weight_);
    hdir->second["res_mtl" ].fill(res(ttbar.tlep().M(), gen->tlep().M()), gen->M(), evt_weight_);
    hdir->second["res_mth_ptth"].fill(res(ttbar.thad().M(), gen->thad().M()), gen->thad().M(), evt_weight_);
    hdir->second["res_deltatt"].fill( 
      res(
        ttbar.top().unit3D().Dot(ttbar.tbar().unit3D()), 
        gen->top().unit3D().Dot(gen->tbar().unit3D())
        ), 
      gen->M(), evt_weight_
      );

    // if(res(ttbar.thad().W().up().E(), gen->thad().W().up().E()) < -0.8) {
    //   cout << dir << endl
    //        << "reco up: " << ttbar.thad().W().up() << endl
    //        << "reco up: " << ttbar.thad().W().down() << endl
    //        << "gen : " << gen->thad().W().up() << endl
    //        << "gen : " << gen->thad().W().down() << endl;
    //   throw 40;
    // }

    hdir->second["res_mw_ptw"].fill(res(ttbar.thad().W().M(), gen->thad().W().M()), gen->thad().W().Pt(), evt_weight_);
    auto gup = gen->thad().W().up();
    auto gdown = gen->thad().W().down();
    auto rup = ttbar.thad().W().up();
    auto rdown = ttbar.thad().W().down();
    hdir->second["res_deltaw_ptw"].fill(
      res(
        rup.unit3D().Dot(rdown.unit3D()),
        gup.unit3D().Dot(gdown.unit3D())
        ),
      gen->thad().W().Pt(), evt_weight_
      );
    if(gup.DeltaR(rup) < gup.DeltaR(rdown)) {
      hdir->second["res_jete"].fill( res(rup.E(), gup.E()), gup.E(), evt_weight_);
      hdir->second["res_jete"].fill( res(rdown.E(), gdown.E()), gdown.E(), evt_weight_);
    } else {
      hdir->second["res_jete"].fill( res(rup.E(), gdown.E()), gdown.E(), evt_weight_);
      hdir->second["res_jete"].fill( res(rdown.E(), gup.E()), gup.E(), evt_weight_);
    }
  }

  void book_fit_plots(string dir) {
    book<TH1F>(dir, "niter" , "", 100, 0, 100);    
    book<TH1F>(dir, "chi2" , "", 100, 0, 50);    
    book<TH1F>(dir, "status" , "", 3, -0.5, 2.5);    
  }

  void fill_fit_plots(string dir, TTKinFitter::KFResult& result) {
    if(result.status == TTKinFitter::FAILED) return;
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }
    hdir->second["niter"].fill(result.niter);
    hdir->second["chi2"].fill(result.chi2);
    hdir->second["status"].fill(result.status);
  }

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
    outFile_.cd();
    vector<double> mass_binning = {200,400,600,800,1000,2000}; 
    for(string lepdir : {"semilep", "fulllep"}) {
      for(string pdir : {"parton", "matched", "selected", "wrong"}) {
        if(lepdir == "fulllep" && pdir != "parton") continue;
        // for(string charge : {"lplus", "lminus"}) {
        string dir = lepdir+"/"+pdir;//+"/"+charge;
        book_kin_plots(dir);        
        book_angular_plots(dir);
      }
    }
    
    string lepdir = "semilep";
    for(string pdir : {"parton", "matched", "selected", "wrong", "obj_selection", "mass_selection", "perm_selection"}) {
      string dir = lepdir+"/"+pdir;//+"/"+charge;
      book_eff_plots(dir);
    }

    for(string pdir : {"matched", "selected", "wrong"}) {
      string dir = lepdir+"/"+pdir;
      book_res_plots(dir);
    }

    for(string pdir : {"matchfit", "selfit", "wrongfit"}) {
      string dir = lepdir+"/"+pdir;
      book_res_plots(dir);
      book_fit_plots(dir);
      book_kin_plots(dir);
    }
  }

  void fill(hyp::TTbar& ttbar, string dir, hyp::TTbar* gen=0) {
    auto hdir = histos_.find(dir);
    if(hdir == histos_.end()) {
      Logger::log().fatal() << dir << "does not exists!" << endl;
      throw 42;
    }  
    fill_kin_plots(dir, ttbar);
    fill_angular_plots(dir, ttbar);    
    if(gen) {
      fill_res_plots(dir, ttbar, gen);
    }
  }

  void print_decay(URStreamer &event) {
    const vector<Genparticle>& gps = event.genParticles();
    const Genparticle *top = 0;
    const Genparticle *tbar = 0;
    for(auto &gp : gps) {
      if(!top && gp.pdgId() == ura::PDGID::t) {
        top = &gp;
      }
      else if(!tbar && gp.pdgId() == ura::PDGID::tbar) {
        tbar = &gp;
      }
    }
    //cout << *top << " " << *tbar << endl;
    cout << *top;
    while(top) {
      cout << " --> ";
      const Genparticle *next=0;
      for(auto &gp : gps) {
        for(auto idx : gp.momIdx()) {
          if(idx == top->idx()) {
            cout << gp << ", ";
            next = &gp;
          }
        }
      }
      top = next;
    }
    cout <<endl;
    throw 43;
  }

  void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) {
    //const vector<Genparticle>& gps = event.genParticles();
    //cout << "genp size: " << gps.size() << endl;
    //print_decay(event);
    tracker_.track("evt processing");
    //select reco objects    
    GenTTBar &ttbar = genp_selector_.ttbar_system();

    string decay = "";
    if(ttbar.type == GenTTBar::DecayType::SEMILEP) decay = "semilep/"; 
    else if(ttbar.type ==  GenTTBar::DecayType::FULLLEP) decay = "fulllep/";
    else return;
    
    tracker_.track("decay mode");
    
    hyp::TTbar ttgen(ttbar);
		// if(ttbar.type == GenTTBar::DecayType::SEMILEP) {
		// 	cout << "-------------------------------------------------------" << endl;
		// 	cout << *ttbar.top.W.first << " " << *ttbar.top.W.second << " " << *ttbar.tbar.W.first << " " << *ttbar.tbar.W.second << " " << endl;
		// 	cout << ttgen.thad().W().up() << " " << ttgen.thad().W().down() << " " << ttgen.tlep().W().l() << " " << ttgen.tlep().W().nu() << " " << endl;
		// }

    fill(ttgen, decay+"parton");
    if(ttbar.type ==  GenTTBar::DecayType::FULLLEP) return;
		//if( fabs(ttgen.M() - 400) > 50 ) return;
    fill_eff_plots(ttgen, decay+"parton");

    hyp::TTbar* gen_ptr = (ttbar.type == GenTTBar::DecayType::SEMILEP) ? &ttgen : 0;
    if( !object_selector_.select(event, shift) ) return;
    tracker_.track("obj selection");

    if( !permutator_.preselection(
          object_selector_.clean_jets(), object_selector_.lepton(), 
          object_selector_.met() ) ) return;

    fill_eff_plots(ttgen, decay+"obj_selection");
    //Gen matching
    Permutation matched;
    matched = matcher_.match(
      genp_selector_.ttbar_system(),
      permutator_.capped_jets(), 
      object_selector_.loose_electrons(),
      object_selector_.loose_muons()
      );
    matched.SetMET(object_selector_.met());

    tracker_.track("gen matching");
    if(matched.IsComplete()) {
      matched.Solve(solver_);
      //fil reco
      hyp::TTbar ttreco(matched);
			tracker_.track("matchfill");
      fill(ttreco, decay+"matched", gen_ptr);
      fill_eff_plots(ttgen, decay+"matched");

      auto result = kinfitter_.fit(matched);
      hyp::TTbar ttfit(matched);
      fill_res_plots(decay+"matchfit", ttfit , gen_ptr);
      fill_fit_plots(decay+"matchfit", result);
      fill_kin_plots(decay+"matchfit", ttfit);
      matched.WJa()->resetp4();
      matched.WJb()->resetp4();
      matched.BHad()->resetp4();
    }    
		tracker_.track("aftermfl");
  
    //Find best permutation
    bool go_on = true;
    Permutation best_permutation;
    size_t nperms = 0;
    while(go_on) {
      Permutation test_perm = permutator_.next(go_on);
      test_perm.LepCharge(object_selector_.lepton_charge());
      if(go_on) {
        nperms++;
        test_perm.Solve(solver_);
        if(test_perm.Prob() < best_permutation.Prob()){
          best_permutation = test_perm;
        }
      }
    }

    if(!best_permutation.IsComplete() || fabs(best_permutation.THad().M() - 173) > 50 || fabs(best_permutation.WHad().M() - 80) > 20) return;
    fill_eff_plots(ttgen, decay+"mass_selection");
    if(best_permutation.Prob() > 1E9) {
      // cout << best_permutation << endl
      //      << *best_permutation.L() << " " << *best_permutation.BLep() << endl;
      // throw 30;
      return;
    }
    fill_eff_plots(ttgen, decay+"perm_selection");

    hyp::TTbar best(best_permutation);
    auto result = kinfitter_.fit(best_permutation);
    hyp::TTbar ttfit(best_permutation);
    if(best_permutation.IsCorrect(matched)) {
      // if(res(ttreco.thad().W().up().E(), gen_ptr->thad().W().up().E()) < -0.6) {
      //   cout << "reco up: " << ttreco.thad().W().up() << endl
      //        << "reco up: " << ttreco.thad().W().down() << endl
      //        << "gen : " << gen_ptr->thad().W().up() << endl
      //        << "gen : " << gen_ptr->thad().W().down() << endl;
      //   throw 40;
      // }
      fill(best, decay+"selected", gen_ptr);
      fill_eff_plots(ttgen, decay+"selected");
      fill_res_plots( decay+"selfit", ttfit ,gen_ptr);
      fill_fit_plots( decay+"selfit", result);
      fill_kin_plots( decay+"selfit", ttfit);
    }
    else {
      fill(best, decay+"wrong", gen_ptr);
      fill_res_plots(decay+"wrongfit", ttfit , gen_ptr);
      fill_fit_plots(decay+"wrongfit", result);
      fill_kin_plots(decay+"wrongfit", ttfit);
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
    Logger::log().debug() << "--DONE--" << endl;

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

			//Long and time consuming
      bool selection = 	genp_selector_.select(event);			
      if(!selection) {
        //Logger::log().error() << "Error: TTGenParticleSelector was not able to find all the generated top decay products in event " << evt_idx << endl <<
        //  "run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
        continue;
			}
      tracker_.track("gen selection");        

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

		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report", opts::value<int>()->default_value(10000), "report every in debug mode");
    parser.addCfgParameter<int>("general", "uselhe", "", 0);
  }
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<topspin_gen> test;
  return test.run();
}
