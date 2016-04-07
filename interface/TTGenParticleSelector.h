#ifndef TTGenParticleSelector_h
#define TTGenParticleSelector_h

#include "GenObject.h"
#include <vector>
#include <list>
#include "LHEParticle.h"

using namespace std;

class TTGenParticleSelector { 
public:
  enum SelMode {NORMAL, PSEUDOTOP, HERWIGPP, FULLDEP, MADGRAPH, LHE};
  TTGenParticleSelector(SelMode mode=NORMAL);
  bool select(URStreamer& event);

  bool is_in_acceptance(GenTTBar::DecayType decay_mode = GenTTBar::DecayType::SEMILEP);
  GenTTBar & ttbar_system() {return ttbar_;}
  GenTTBar & ttbar_final_system() {return ttbar_;}
  vector<Genjet*>& additional_jets() {return added_jets_;}
  void setmode(SelMode mode) {
    mode_ = mode;
    if(mode_ == MADGRAPH) w_decay_momid_=6;
  }

private:
  //GEN PARTICLE SELECTION
  void select_normal(URStreamer& event);
  void select_pstop(URStreamer& event);
  void select_herwig(URStreamer& event);

  int w_decay_momid_=24;
  SelMode mode_;

  //counters storage and stuff
  
  list<GenObject> selected_;
  vector<GenObject*> wpartons_;
  vector<GenObject*> charged_leps_;
  vector<GenObject*> neutral_leps_;
  vector<GenObject*> final_charged_leps_;
  list<Genjet> jets_;
  vector<Genjet*> added_jets_;
  GenObject* top_;
  GenObject* tbar_;
  GenObject* b_;
  GenObject* bbar_;
  GenTTBar ttbar_;
  GenTTBar ttbar_final_; //after ration
  //GenObject* wplus_, wminus_;
  //NEEDED?
  // GenObject* genbl_;
  // GenObject* genbh_;
  // TLorentzVector gentoplep;
  // TLorentzVector gentophad;
  int topcounter_;
  int lepdecays_;
  float jetptmin_, jetetamax_;

  //ACCEPTANCE
  int is_in_acceptance_=-1;
  float cut_lep_ptmin_=0., cut_lep_etamax_=999.;
  float cut_wjet_pthard_=0., cut_wjet_ptsoft_=0., cut_wjet_etamax_=999.;
  float cut_bjet_ptsoft_=0., cut_bjet_pthard_=0., cut_bjet_etamax_=999.;
  float cut_parton_separation_=0.;

};

#endif
