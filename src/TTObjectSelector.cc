#include "TTObjectSelector.h"
#include "Logger.h"
#include "URParser.h"
#include "TMath.h"
#include "TRandom3.h"
#include "DataFile.h"

using namespace TMath;

TTObjectSelector::TTObjectSelector(int objsel):
  objsel_(objsel),
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
    parser.addCfgParameter<std::string>("general", "JER", "source file for jet uncertainties");

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

    DataFile jes(parser.getCfgPar<std::string>("general", "jet_uncertainties"));
    DataFile jer(parser.getCfgPar<std::string>("general", "JER"));    
    jet_scaler_.Init(
      jes.path(), jer.path()
      );

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

bool TTObjectSelector::select_muons(URStreamer &event, systematics::SysShifts shift) {
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
  return (loose_muons_.size() > 0);
}

bool TTObjectSelector::select_electrons(URStreamer &event, systematics::SysShifts shift) {
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
  return (loose_electrons_.size() > 0);
}

bool TTObjectSelector::select_jetmet(URStreamer &event, systematics::SysShifts shift) {
	double  metcorrx = 0.;
	double  metcorry = 0.;
	const vector<Jet>& jets = event.jets();
	for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
	{
		IDJet jet(*jetit);
    double scale_factor = 0;
		if(shift == systematics::SysShifts::JES_DW)      scale_factor = jet_scaler_.GetUncP(jet);
		else if(shift == systematics::SysShifts::JES_UP) scale_factor = -1*jet_scaler_.GetUncM(jet);

    if(event.run == 1) { //is MC: apply JER
      double jer_sigma = 0;
      if(shift == systematics::SysShifts::JER_UP) jer_sigma = jet_scaler_.JERSigmaUp(jet); 
      else if(shift == systematics::SysShifts::JER_DW) jer_sigma = jet_scaler_.JERSigmaDw(jet); 
      else jer_sigma = jet_scaler_.JERSigma(jet);
      double jer_sf = gRandom->Gaus(0., jer_sigma);
      scale_factor += jer_sf;
    }

    //set corrections for MET
		metcorrx -= scale_factor*jet.Px(); 
		metcorry -= scale_factor*jet.Py(); 

    //Set shifted P
    scale_factor += 1;
		jet.SetPxPyPzE(
      jet.Px() * scale_factor, 
      jet.Py() * scale_factor, 
      jet.Pz() * scale_factor, 
      jet.E()  * scale_factor
      );
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
    if(shift == systematics::SysShifts::MET_DW) shift = -1.0;
    if(shift == systematics::SysShifts::MET_UP) shift = 1.0;
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

bool TTObjectSelector::pass_through(URStreamer &event, systematics::SysShifts shift) {
  reset();
  select_muons(event, shift);
  select_electrons(event, shift);
  select_jetmet(event, shift);
  return true;
}

bool TTObjectSelector::select(URStreamer &event, systematics::SysShifts shift) {
  reset();

  //Trigger!
  bool isMC = (event.run == 1);
  bool electron_trig = (isMC) ? (event.trigger().HLT_Ele27_WPLoose_Gsf() == 1) : (event.trigger().HLT_Ele27_eta2p1_WPLoose_Gsf() == 1);
  bool muon_trig = (isMC) ? (event.trigger().HLT_IsoMu27() == 1) : (event.trigger().HLT_IsoMu20() == 1 || event.trigger().HLT_IsoTkMu20());
  if(!(electron_trig || muon_trig)) return false;
  if(tracker_) tracker_->track("trigger");

  bool has_muons = select_muons(event, shift);
  bool has_electrons = select_electrons(event, shift);
  
  if(!(has_muons || has_electrons)) return false;
  if(tracker_) tracker_->track("has leptons");

  if(objsel_ == 0  && (tight_muons_.size()+medium_electrons_.size() != 1)) return false;
  if(objsel_ == -1 && (medium_electrons_.size() != 1)) return false;
  if(objsel_ == 1  && (tight_muons_.size() != 1)) return false;
  //1 tight lepton and no loose ones
  if(loose_electrons_.size() + loose_muons_.size() != 1) return false;
  if(tracker_) tracker_->track("lepton veto");
  
  //right trigger
  if(has_muons && !muon_trig) return false;
  if(has_electrons && !electron_trig) return false;
  if(tracker_) tracker_->track("right trigger");
  
  select_jetmet(event, shift);
  if(clean_jets_.size() < cut_nminjets_) return false;
 
  if(tracker_) tracker_->track("has N jets");
  return true;
}

TLorentzVector* TTObjectSelector::lepton() {
  return (tight_muons_.size() == 1) ? dynamic_cast<TLorentzVector*>(tight_muons_[0]) : dynamic_cast<TLorentzVector*>(medium_electrons_[0]);
}

int TTObjectSelector::lepton_charge() {
  return (tight_muons_.size() == 1) ? tight_muons_[0]->charge() : medium_electrons_[0]->charge();
}
