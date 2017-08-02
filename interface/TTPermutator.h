#ifndef TTPermutator_h
#define TTPermutator_h

#include "Analyses/URTTbar/interface/Permutation.h"
#include <vector>
#include "URAnalysis/AnalysisFW/interface/URSelector.h"
#include <array>
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include <string>
#include <list>
#include <iostream>

class TTPermutator: public URSelector {
public:
  TTPermutator(std::string cfgname="permutations");
  bool preselection(vector<IDJet*> jets, TLorentzVector* lepton, IDMet* met, int lc=0, double rho=-1, bool track=true);
  list<Permutation>& permutations() {
  	if(!has_run_) permutate();
  	return permutations_;
  }

/// Joseph added for perm disc
  list<Permutation>& permutations_3J(IDJet* wj1, IDJet* wj2, IDJet* bj1, IDJet* bj2, TLorentzVector* lep, IDMet* met, int lc) {
  	if(!has_run_) permutate_3J(wj1, wj2, bj1, bj2, lep, met, lc);
  	return permutations_3J_;
  }

  void reset_3J() {
		has_run_ = false;
		permutations_3J_.clear();
    jets_.clear();
    capped_jets_.clear();
    lepton_ = 0;
    met_ = 0;
  }

///
  void reset() {
		has_run_ = false;
		permutations_.clear();
    jets_.clear();
    capped_jets_.clear();
    lepton_ = 0;
    met_ = 0;
  }

  virtual void configure() override;
  std::vector<IDJet*>& capped_jets() {return capped_jets_;}

  void set_tracker(CutFlowTracker *t) {tracker_ = t;}

  IDJet::BTag tight_bID_cut() {return cut_tight_b_;}
  IDJet::BTag loose_bID_cut() {return cut_loose_b_;}

  size_t njets_max() {return cut_max_jets_;}

private:
	void permutate();
	list<Permutation> permutations_;

/// Joseph added for perm disc
  void permutate_3J(IDJet* wj1, IDJet* wj2, IDJet* bj1, IDJet* bj2, TLorentzVector* lep, IDMet* met, int lc);
  list<Permutation> permutations_3J_;

//  void permutate_3J(IDJet* wja, IDJet* bjh, IDJet* bjl, TLorentzVector* lep);
//  list<Permutation> permutations_3J_;

///


  std::vector<IDJet*> jets_;
  std::vector<IDJet*> capped_jets_;
  TLorentzVector* lepton_ = 0;
  IDMet* met_ = 0;
  bool is_configured_ = false;
  int lcharge_=0;
  double rho_=-1;
	bool has_run_ = false;

  //tracker
  CutFlowTracker *tracker_=0;

  //cuts
  IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;
  size_t bjet_idx_limit_=0, cut_max_jets_=99999;
  float cut_bjetpt_hard_=0., cut_bjetpt_soft_=0., cut_wjetpt_hard_=0., cut_wjetpt_soft_=0.;
};

#endif
