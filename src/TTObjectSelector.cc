#include "TTObjectSelector.h"
#include "Logger.h"
#include "URParser.h"
#include "TMath.h"
#include "TRandom3.h"
#include "DataFile.h"

using namespace TMath;

const std::map<std::string, TTObjectSelector::SysShifts> TTObjectSelector::name_to_shift = {
  {"nosys", TTObjectSelector::SysShifts::NOSYS}, 
  {"jes_up", TTObjectSelector::SysShifts::JES_UP}, 
  {"jes_down", TTObjectSelector::SysShifts::JES_DW}, 
  {"jer_up", TTObjectSelector::SysShifts::JER_UP},
  {"jer_down", TTObjectSelector::SysShifts::JER_DW}, 
  {"met_up", TTObjectSelector::SysShifts::MET_UP},
  {"met_down", TTObjectSelector::SysShifts::MET_DW}
};
const std::map<TTObjectSelector::SysShifts, std::string> TTObjectSelector::shift_to_name = {
  {TTObjectSelector::SysShifts::NOSYS , "nosys"   }, 
  {TTObjectSelector::SysShifts::JES_UP, "jes_up"  }, 
  {TTObjectSelector::SysShifts::JES_DW, "jes_down"}, 
  {TTObjectSelector::SysShifts::JER_UP, "jer_up"  },
  {TTObjectSelector::SysShifts::JER_DW, "jer_down"}, 
  {TTObjectSelector::SysShifts::MET_UP, "met_up"  },
  {TTObjectSelector::SysShifts::MET_DW, "met_down"}
};

TTObjectSelector::TTObjectSelector():
  sel_jets_(),
  clean_jets_(),
  sel_muons_(),
  loose_muons_(),
  tight_muons_(),
  sel_electrons_(),
  loose_electrons_(),
  medium_electrons_(),
  jet_scaler_() {
    URParser &parser = URParser::instance();
    parser.addCfgParameter<std::string>("general", "jet_uncertainties", "source file for jet uncertainties");

    //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
    parser.addCfgParameter<std::string>("loose_muons", "id", "ID to be applied");
    parser.addCfgParameter<float>      ("loose_muons", "ptmin", "minimum pt");
    parser.addCfgParameter<float>      ("loose_muons", "etamax", "maximum eta");

    parser.addCfgParameter<std::string>("tight_muons", "id", "ID to be applied");
    parser.addCfgParameter<float>      ("tight_muons", "ptmin", "minimum pt");
    parser.addCfgParameter<float>      ("tight_muons", "etamax", "maximum eta");

    parser.addCfgParameter<std::string>("loose_electrons", "id", "ID to be applied");
    parser.addCfgParameter<float>      ("loose_electrons", "ptmin", "minimum pt");
    parser.addCfgParameter<float>      ("loose_electrons", "etamax", "maximum eta");
      
    parser.addCfgParameter<std::string>("tight_electrons", "id", "ID to be applied");
    parser.addCfgParameter<float>      ("tight_electrons", "ptmin", "minimum pt");
    parser.addCfgParameter<float>      ("tight_electrons", "etamax", "maximum eta");
      
    parser.addCfgParameter<int>  ("jets", "n_min", "minimum number of jets");
    parser.addCfgParameter<float>("jets", "ptmin", "minimum pt");
    parser.addCfgParameter<float>("jets", "etamax", "maximum eta");
    
    parser.parseArguments();
    parser.parseArguments();

    jet_scaler_.Init(DataFile(parser.getCfgPar<std::string>("general", "jet_uncertainties")).path());

    cut_loosemu_id_ = IDMuon::id(
      parser.getCfgPar<std::string>("loose_muons", "id"    ));
    cut_loosemu_ptmin_ = parser.getCfgPar<float>("loose_muons", "ptmin" );
    cut_loosemu_etamax_ = parser.getCfgPar<float>("loose_muons", "etamax");

    cut_tightmu_id_ = IDMuon::id(
      parser.getCfgPar<std::string>("tight_muons", "id"    ));
    cut_tightmu_ptmin_ = parser.getCfgPar<float>("tight_muons", "ptmin" );
    cut_tightmu_etamax_ = parser.getCfgPar<float>("tight_muons", "etamax");

    cut_looseel_id_ = IDElectron::id(
      parser.getCfgPar<std::string>("loose_electrons", "id"    ));
    cut_looseel_ptmin_ = parser.getCfgPar<float>("loose_electrons", "ptmin" );
    cut_looseel_etamax_ = parser.getCfgPar<float>("loose_electrons", "etamax");
    
    cut_tightel_id_ = IDElectron::id(
      parser.getCfgPar<std::string>("tight_electrons", "id"    ));
    cut_tightel_ptmin_ = parser.getCfgPar<float>("tight_electrons", "ptmin" );
    cut_tightel_etamax_ = parser.getCfgPar<float>("tight_electrons", "etamax");
    
    cut_nminjets_ = parser.getCfgPar<int>  ("jets", "n_min" );
    cut_jet_ptmin_ = parser.getCfgPar<float>("jets", "ptmin" );
    cut_jet_etamax_ = parser.getCfgPar<float>("jets", "etamax");  
}

void TTObjectSelector::reset() {
  sel_jets_.clear();
  clean_jets_.clear();
  sel_muons_.clear();
  loose_muons_.clear();
  tight_muons_.clear();
  sel_electrons_.clear();
  loose_electrons_.clear();
  medium_electrons_.clear();
}

bool TTObjectSelector::select_muons(URStreamer &event, SysShifts shift) {
	const vector<Muon>& muons = event.muons();
	for(vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon)	{
		IDMuon mu(*muon, event.rho().value());
		if(mu.ID(cut_loosemu_id_) && mu.Pt() > cut_loosemu_ptmin_ && Abs(mu.Eta()) < cut_loosemu_etamax_) {
			sel_muons_.push_back(mu);
			loose_muons_.push_back(&(sel_muons_.back()));
			if(mu.ID(cut_tightmu_id_) && mu.Pt() > cut_tightmu_ptmin_ && Abs(mu.Eta()) < cut_tightmu_etamax_)	{
				tight_muons_.push_back(&(sel_muons_.back()));
			}
		}
	}
  return loose_muons_.size() > 0;
}

bool TTObjectSelector::select_electrons(URStreamer &event, SysShifts shift) {
	const vector<Electron>& electrons = event.electrons();
	for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron)	{
		IDElectron el(*electron, event.rho().value());
		if(el.ID(cut_looseel_id_) && el.Pt() > cut_looseel_ptmin_ && Abs(el.Eta()) < cut_looseel_etamax_) {
			sel_electrons_.push_back(el);
			loose_electrons_.push_back(&(sel_electrons_.back()));
			if(el.ID(cut_tightel_id_) && el.Pt() > cut_tightel_ptmin_ && Abs(el.Eta()) < cut_tightel_etamax_) {
				medium_electrons_.push_back(&(sel_electrons_.back()));
			}
		}
	}
  return loose_electrons_.size() > 0;
}

bool TTObjectSelector::select_jetmet(URStreamer &event, SysShifts shift) {
	double  metcorrx = 0.;
	double  metcorry = 0.;
	const vector<Jet>& jets = event.jets();
	for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
	{
		IDJet jet(*jetit);
		//double sf = gRandom->Gaus(1., 0.05);
		double sf = 1; //csigmajet; ??
		if(shift == SysShifts::JES_DW){sf *= -1*Abs(jet_scaler_.GetUncM(jet));}
		if(shift == SysShifts::JES_UP){sf *= Abs(jet_scaler_.GetUncP(jet));}
		if(shift == SysShifts::JER_UP){sf *= gRandom->Gaus(sf, 0.1);} //FIXME, update as soon as available
		//if(shift == SysShifts::JER_DW){sf *= gRandom->Gaus(sf, cjetres);}
		sf += 1;

    //set corrections for MET
		metcorrx -= (sf-1)*jet.Px(); 
		metcorry -= (sf-1)*jet.Py(); 

    //Set shifted P
		jet.SetPxPyPzE(jet.Px()*sf, jet.Py()*sf, jet.Pz()*sf, jet.E()*sf);
		if(jet.Pt() < cut_jet_ptmin_ || Abs(jet.Eta()) > cut_jet_etamax_) {continue;}
		if(!jet.ID() || !jet.Clean(loose_muons_, loose_electrons_)) {continue;}

		sel_jets_.push_back(jet);
		clean_jets_.push_back(&(sel_jets_.back()));
	}
  if( clean_jets_.size() ==0) return false;

	//const vector<Nohfmet>& mets = event.NoHFMETs();
	const vector<Met>& mets = event.METs();
	if(mets.size() == 1) {
    float shift = 0;
    if(shift == SysShifts::MET_DW) shift = 11.0;
    if(shift == SysShifts::MET_UP) shift = 1.0;
		met_ = mets[0];
		met_.SetPx(met_.Px() + metcorrx + shift*met_.pxunc());
		met_.SetPy(met_.Py() + metcorry + shift*met_.pyunc());
		met_.SetE(Sqrt(met_.Px()*met_.Px() + met_.Py()*met_.Py()));
	}
  else {
    Logger::log().error() << "met collection has more than one entry! which I should pick?" << std::endl;
  }
  return true;
}


bool TTObjectSelector::select(URStreamer &event, SysShifts shift) {
  reset();

  //Trigger!
  /*bool isMC = (event.run == 1);
  bool electron_trig = (isMC) ? (event.trigger().HLT_Ele27_WP85_Gsf() == 1) : (event.trigger().HLT_Ele27_eta2p1_WPLoose_Gsf() == 1);
  bool muon_trig = (isMC) ? (event.trigger().HLT_IsoMu27() == 1) : (event.trigger().HLT_IsoMu24_eta2p1() == 1);
  if(!(electron_trig || muon_trig)) return false;
  if(tracker_) tracker_->track("trigger");*/

  bool has_muons = select_muons(event, shift);
  bool has_electrons = select_electrons(event, shift);
  if(!(has_muons || has_electrons)) return false;
  if(tracker_) tracker_->track("has leptons");

  if(tight_muons_.size()+medium_electrons_.size() != 1) return false;
  //1 tight lepton and no loose ones
  if(loose_electrons_.size() + loose_muons_.size() > 1) return false;
  if(tracker_) tracker_->track("lepton veto");
  
  //right trigger
  /*if(has_muons && !muon_trig) return false;
    if(has_electrons && !electron_trig) return false;
  if(tracker_) tracker_->track("right trigger");*/
  
  select_jetmet(event, shift);
  return (clean_jets_.size() >= cut_nminjets_);
}

TLorentzVector* TTObjectSelector::lepton() {
  return (tight_muons_.size() == 1) ? dynamic_cast<TLorentzVector*>(tight_muons_[0]) : dynamic_cast<TLorentzVector*>(medium_electrons_[0]);
}

int TTObjectSelector::lepton_charge() {
  return (tight_muons_.size() == 1) ? tight_muons_[0]->charge() : medium_electrons_[0]->charge();
}
