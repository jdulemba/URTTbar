#ifndef TTPermutator_h
#define TTPermutator_h

#include "Permutation.h"
#include <vector>
#include "URSelector.h"
#include <array>
#include "CutFlowTracker.h"

class TTPermutator: public URSelector {
public:
  TTPermutator();
  bool preselection(vector<IDJet*> jets, TLorentzVector* lepton, IDMet* met, int lc=0);
  Permutation next(bool &keep_going);
  void reset() {
    jets_.clear();
    capped_jets_.clear();
    lepton_ = 0;
    met_ = 0;
    jet_pos_ = {0, 0, 0, 0};
  }

  virtual void configure() override;
  std::vector<IDJet*>& capped_jets() {return capped_jets_;}

  void set_tracker(CutFlowTracker *t) {tracker_ = t;}

  IDJet::BTag tight_bID_cut() {return cut_tight_b_;}
  IDJet::BTag loose_bID_cut() {return cut_loose_b_;}

  size_t njets_max() {return cut_max_jets_;}

private:
  std::vector<IDJet*> jets_;
  std::vector<IDJet*> capped_jets_;
  TLorentzVector* lepton_ = 0;
  IDMet* met_ = 0;
  bool is_configured_ = false;
  int lcharge_=0;

  //tracker
  CutFlowTracker *tracker_=0;

  //permutation status
  std::array<size_t, 4> jet_pos_;

  //cuts
  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;
  size_t bjet_idx_limit_=0, cut_max_jets_=99999;
  float cut_bjetpt_hard_=0., cut_bjetpt_soft_=0., cut_wjetpt_hard_=0., cut_wjetpt_soft_=0.;
};

#endif
