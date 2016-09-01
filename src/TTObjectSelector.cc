#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "TMath.h"
#include "TRandom3.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"

using namespace TMath;

TTObjectSelector::TTObjectSelector(int objsel):
  objsel_(objsel),
  sel_jets_(),
  clean_jets_(),
  sel_muons_(),
  veto_muons_(),
	vmu_sel_("veto_muons"),
  loose_muons_(),
	lmu_sel_("loose_muons"),
  tight_muons_(),
	tmu_sel_("tight_muons"),
  sel_electrons_(),
  veto_electrons_(),
	vel_sel_("veto_electrons"),
  loose_electrons_(),
	lel_sel_("loose_electrons"),
  tight_electrons_(),
	tel_sel_("tight_electrons") {
    URParser &parser = URParser::instance();

    //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
    parser.addCfgParameter<int>("event", "use_trig"   , "", 1);
    parser.addCfgParameter<int>("event", "use_filters", "", 1);
    parser.addCfgParameter<int>("event", "smear_met", "", 1);
		parser.addCfgParameter<int>("event", "trig_config", "", 2015);
		parser.addCfgParameter<int>("event", "allow_loose", "", 0);

    parser.addCfgParameter<int>  ("jets", "n_min", "minimum number of jets");
    parser.addCfgParameter<int>  ("jets", "applyJER", "");
    parser.addCfgParameter<float>("jets", "ptmin", "minimum pt");
    parser.addCfgParameter<float>("jets", "etamax", "maximum eta");
    
    parser.parseArguments();
		use_trg_     = (parser.getCfgPar<int>("event", "use_trig"   ) == 1);
		trg_cfg_     = parser.getCfgPar<int>("event", "trig_config");
		use_filters_ = (parser.getCfgPar<int>("event", "use_filters") == 1);
    smear_met_   = (parser.getCfgPar<int>("event", "smear_met") == 1);
		allow_loose_ = (parser.getCfgPar<int>("event", "allow_loose") == 1);
    
    cut_nminjets_ = parser.getCfgPar<int>  ("jets", "n_min" );
    apply_jer_ = (parser.getCfgPar<int>("jets", "applyJER" ) > 0.5);
    Logger::log().debug() << "Running JER: " << apply_jer_ << std::endl;
    cut_jet_ptmin_ = parser.getCfgPar<float>("jets", "ptmin" );
    cut_jet_etamax_ = parser.getCfgPar<float>("jets", "etamax");  
}

void TTObjectSelector::reset() {
	pass_lepton_=0;
	el_trg_ = false;
	mu_trg_ = false;
	evt_type_ = EvtType::NOTSET;
  sel_jets_.clear();
  clean_jets_.clear();
  sel_muons_.clear();
  veto_muons_.clear();
  loose_muons_.clear();
  tight_muons_.clear();
  sel_electrons_.clear();
  veto_electrons_.clear();
  loose_electrons_.clear();
  tight_electrons_.clear();
}

bool TTObjectSelector::select_muons(URStreamer &event, systematics::SysShifts shift) {
	const vector<Muon>& muons = event.muons();
	for(vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon)	{
		IDMuon mu(*muon, event.rho().value());
		bool isveto  = vmu_sel_.pass(mu);
		bool isloose = lmu_sel_.pass(mu);
		bool istight = tmu_sel_.pass(mu);
		if(isveto || isloose || istight) {
			sel_muons_.push_back(mu);
			if(isveto ) veto_muons_.push_back(&(sel_muons_.back()));
			if(isloose) loose_muons_.push_back(&(sel_muons_.back()));
			if(istight) tight_muons_.push_back(&(sel_muons_.back()));
		}
	}
  return (loose_muons_.size() > 0 || tight_muons_.size() > 0);
}

bool TTObjectSelector::select_electrons(URStreamer &event, systematics::SysShifts shift) {
	const vector<Electron>& electrons = event.electrons();
	for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron)	{
		if(event.rho().value() < 0.) cout << "RHO: " << event.rho().value() << endl;
		IDElectron el(*electron, event.rho().value());
		bool isveto  = vel_sel_.pass(el);
		bool isloose = lel_sel_.pass(el);
		bool istight = tel_sel_.pass(el);
		if(isveto || isloose || istight) {
			sel_electrons_.push_back(el);
			if(isveto ) veto_electrons_.push_back(&(sel_electrons_.back()));
			if(isloose) loose_electrons_.push_back(&(sel_electrons_.back()));
			if(istight) tight_electrons_.push_back(&(sel_electrons_.back()));
		}
	}
  return (loose_electrons_.size() > 0 || tight_electrons_.size() > 0);
}

bool TTObjectSelector::select_jetmet(URStreamer &event, systematics::SysShifts shift) {
	//SET JETS
	const vector<Jet>& jets = event.jets();
	for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
	{
		IDJet jet(*jetit);
		//set shifts and systematics
		if(shift == systematics::SysShifts::JER_UP) jet.scale(jet.JERUp());
		else if(shift == systematics::SysShifts::JER_DW) jet.scale(jet.JERDown());
		else if(event.run == 1 && apply_jer_) jet.scale(jet.JER());

		if(shift == systematics::SysShifts::JES_DW)      jet.scale(1-jet.JESUnc());
		else if(shift == systematics::SysShifts::JES_UP) jet.scale(1+jet.JESUnc());

		if(jet.Pt() < cut_jet_ptmin_ || Abs(jet.Eta()) > cut_jet_etamax_) {continue;}
		if(!jet.ID() || !jet.Clean(veto_muons_, veto_electrons_)) {continue;}

		sel_jets_.push_back(jet);
		clean_jets_.push_back(&(sel_jets_.back()));
	}
  if( clean_jets_.size() ==0) return false;

	//SET MET
	const vector<Met>& mets = event.METs();
	if(mets.size() == 1) {
		met_ = mets[0];
		if(event.run == 1 && apply_jer_ && smear_met_)   met_.setvect(met_.pxsmear(), met_.pysmear());
		if(shift == systematics::SysShifts::JER_UP)      met_.shiftvect(   met_.pxuncJER(),    met_.pyuncJER());
		else if(shift == systematics::SysShifts::JER_DW) met_.shiftvect(-1*met_.pxuncJER(), -1*met_.pyuncJER());
		else if(shift == systematics::SysShifts::JES_DW) met_.shiftvect(   met_.pxuncJES(),    met_.pyuncJES());
		else if(shift == systematics::SysShifts::JES_UP) met_.shiftvect(-1*met_.pxuncJES(), -1*met_.pyuncJES());
    else if(shift == systematics::SysShifts::MET_DW) met_.shiftvect(met_.pxunc(), met_.pyunc());
    else if(shift == systematics::SysShifts::MET_UP) met_.shiftvect(-1*met_.pxunc(), -1*met_.pyunc());
	}
  else {
    Logger::log().error() << "met collection has more than one entry! which I should pick?" << std::endl;
  }
  return true;
}

bool TTObjectSelector::pass_through(URStreamer &event, systematics::SysShifts shift) {
  reset();
  select_muons(event, shift);
  select_electrons(event, shift);
  select_jetmet(event, shift);
  return true;
}

bool TTObjectSelector::pass_trig(URStreamer &event, systematics::SysShifts shift) {
  //Trigger!
  bool isMC = (event.run == 1);
	switch(trg_cfg_) {
	case 2015:
		el_trg_ = (event.trigger().HLT_Ele23_WPLoose_Gsf() == 1);
		mu_trg_ = (event.trigger().HLT_IsoMu20() == 1 || event.trigger().HLT_IsoTkMu20() == 1);
		//std::cout << "TRG: " << event.trigger().HLT_IsoMu20() << " -- " << event.trigger().HLT_IsoTkMu20() << std::endl;
		break;
	case 2016:
		el_trg_ = (event.trigger().HLT_Ele27_WPLoose_Gsf() == 1);
		mu_trg_ = (event.trigger().HLT_IsoMu22() == 1 || event.trigger().HLT_IsoTkMu22() == 1);
		break;
	default:
		Logger::log().fatal() << "Trigger configuration can only be 2015 or 2016" << std::endl;
		throw 42;
	}

	return !((use_trg_ || !isMC) && !(el_trg_ || mu_trg_));
}

bool TTObjectSelector::pass_filter(URStreamer &event, systematics::SysShifts shift) {
	//Filters
	auto filters = event.filter();
	bool filter_answer = (filters.Flag_HBHENoiseFilter() == 1); 
	filter_answer &= (filters.Flag_HBHENoiseIsoFilter() == 1); 
	filter_answer &= (filters.Flag_CSCTightHalo2015Filter() == 1);
	filter_answer &= (filters.Flag_EcalDeadCellTriggerPrimitiveFilter() == 1);
	filter_answer &= (filters.Flag_eeBadScFilter() == 1);
	// std::cout << filters.Flag_HBHENoiseFilter() << " " << filters.Flag_HBHENoiseIsoFilter()  << " "
	// 					<< filters.Flag_CSCTightHaloFilter() << " " << filters.Flag_EcalDeadCellTriggerPrimitiveFilter() << " "
	// 					<< filters.Flag_eeBadScFilter() << std::endl;
	return filter_answer;
}

bool TTObjectSelector::select(URStreamer &event, systematics::SysShifts shift) {
  reset();

	if(!pass_trig(event, shift)) return false;
  if(tracker_) tracker_->track("trigger");

	bool filter_answer = pass_filter(event, shift);
	if(use_filters_ && !filter_answer) return false;
  if(tracker_) tracker_->track("Evt filters");

  bool has_muons = select_muons(event, shift);
  bool has_electrons = select_electrons(event, shift);
  
  if(!(has_muons || has_electrons)) return false;
  if(tracker_) tracker_->track("has leptons");

  //1 lepton and no veto ones
	if(tight_muons_.size()) {
		if(veto_electrons_.size()) return false;
		if(nvetos(veto_muons_, tight_muons_[0])) return false;
		evt_type_ = TIGHTMU;
		pass_lepton_ = 1;
	}
	else if(tight_electrons_.size()) {
		if(veto_muons_.size()) return false;
		if(nvetos(veto_electrons_, tight_electrons_[0])) return false;
		evt_type_ = TIGHTEL;
		pass_lepton_ = -1;
	}
	else if(loose_muons_.size()) {
		if(veto_electrons_.size()) return false;
		if(loose_electrons_.size()) return false;
		if(nvetos(veto_muons_, loose_muons_[0])) return false;
		vector<IDJet*> new_jets;
		for(auto jet : clean_jets_) {
			if(jet->Clean(loose_muons_, loose_electrons_)) new_jets.push_back(jet);
		}
		clean_jets_ = new_jets;
		evt_type_ = LOOSEMU;
		pass_lepton_ = 1;
	}
	else if(loose_electrons_.size()) {
		if(veto_muons_.size()) return false;
		if(loose_muons_.size()) return false;
		if(nvetos(veto_electrons_, loose_electrons_[0])) return false;
		vector<IDJet*> new_jets;
		for(auto jet : clean_jets_) {
			if(jet->Clean(loose_muons_, loose_electrons_)) new_jets.push_back(jet);
		}
		clean_jets_ = new_jets;
		evt_type_ = LOOSEEL;
		pass_lepton_ = -1;
	}
	if(!allow_loose_ && (evt_type_ == LOOSEMU || evt_type_ == LOOSEEL)) return false;

	//right lepton
  if(objsel_ == -1 && (evt_type_ == TIGHTMU || evt_type_ == LOOSEMU)) return false;
  if(objsel_ == 1  && (evt_type_ == TIGHTEL || evt_type_ == LOOSEEL)) return false;
  
  //right trigger 
  if(has_muons && !mu_trg_) return false;
  if(has_electrons && !el_trg_) return false;
  if(tracker_) tracker_->track("right trigger");
  
  select_jetmet(event, shift);
  if(clean_jets_.size() < cut_nminjets_) return false;
 
  if(tracker_) tracker_->track("has N jets");
  return true;
}

TLorentzVector* TTObjectSelector::lepton() {
  return (pass_lepton_ == 1) ? dynamic_cast<TLorentzVector*>(muon()) : dynamic_cast<TLorentzVector*>(electron());
}

int TTObjectSelector::lepton_charge() {
  return (pass_lepton_ == 1) ? muon()->charge() : electron()->charge();
}
