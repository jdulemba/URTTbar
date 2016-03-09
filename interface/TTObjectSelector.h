#ifndef TTObjectSelector_h
#define TTObjectSelector_h

#include "URStreamer.h"
#include <list>
#include <vector>
#include "IDMuon.h"
#include "IDElectron.h"
#include "IDJet.h"
#include "IDMet.h"
#include "JetScale.h"
#include "JetScaler.h"
#include <map>
#include <string>
#include "CutFlowTracker.h"
#include "systematics.h"

using namespace std;
class CutFlowTracker;

class TTObjectSelector {
public:
  TTObjectSelector(int objsel=0); //-1 e only, 1 mu only, 0 any
  
  void lepton_type(int ltype) {objsel_=ltype;} //-1 e only, 1 mu only, 0 any
  void reset();
  bool select(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);

  //getters for selected objects
  vector<IDMuon*>& loose_muons() {return loose_muons_;}
  vector<IDMuon*>& tight_muons() {return tight_muons_;}

  vector<IDElectron*>& loose_electrons() {return loose_electrons_;}
  vector<IDElectron*>& medium_electrons() {return medium_electrons_;}

  vector<IDJet*>& clean_jets() {return clean_jets_;}
  IDMet * met() {return &met_;}

  TLorentzVector* lepton();
  int lepton_charge();

  void set_tracker(CutFlowTracker *t) {tracker_ = t;}

private:

  bool select_muons(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
  bool select_electrons(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
  bool select_jetmet(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);

  list<IDJet> sel_jets_;
  vector<IDJet*> clean_jets_;
  list<IDMuon> sel_muons_;
  vector<IDMuon*> loose_muons_;
  vector<IDMuon*> tight_muons_;
  list<IDElectron> sel_electrons_;
  vector<IDElectron*> loose_electrons_;
  vector<IDElectron*> medium_electrons_;
  IDMet met_;
  JetScaler jet_scaler_;  

  //tracker
  CutFlowTracker *tracker_=0;

  //cuts
  bool is_configured_ = false;
  int objsel_=0;
  IDMuon::IDS cut_loosemu_id_, cut_tightmu_id_;
  float cut_loosemu_ptmin_, cut_loosemu_etamax_, cut_tightmu_ptmin_, cut_tightmu_etamax_;
  
  IDElectron::IDS cut_looseel_id_, cut_tightel_id_;
  float cut_looseel_ptmin_, cut_looseel_etamax_, cut_tightel_ptmin_, cut_tightel_etamax_;

  float cut_jet_ptmin_, cut_jet_etamax_;
  size_t cut_nminjets_;
};

#endif
