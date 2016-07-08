#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "TMath.h"
#include <algorithm>
#include "URAnalysis/AnalysisFW/interface/Logger.h"

using namespace TMath;
using namespace std;

TTPermutator::TTPermutator():
  URSelector("permutations"),
  jets_(),
  capped_jets_() {
  jet_pos_ = {{0, 0, 0, 0}};
  configure();
}

void TTPermutator::configure() {
  URSelector::configure();
  if(!is_configured_) {    
    addCfgParameter<size_t>("max_jets", "maximum number of jets to be considered.", 999999);
    addCfgParameter<size_t>("max_bjets", "maximum number of b jets to be considered. (0=disabled)", 0);
    addCfgParameter<float>("hardb_pt", "min pt of the hardest b-jet", 0.);
    addCfgParameter<float>("softb_pt", "min pt of the softer b-jet", 0.);
    addCfgParameter<float>("hardw_pt", "min pt of the hardest w-jet", 0.);
    addCfgParameter<float>("softw_pt", "min pt of the softer w-jet", 0.);
    addCfgParameter<std::string>("tightb", "tightest BTag working point to be applied");
    addCfgParameter<std::string>("looseb", "loosest BTag working point to be applied");

    registerConfiguration();

    cut_max_jets_ = getCfgParameter<size_t>("max_jets" );
    bjet_idx_limit_ = getCfgParameter<size_t>("max_bjets");
    cut_bjetpt_hard_ = getCfgParameter<float>( "hardb_pt" );
    cut_bjetpt_soft_ = getCfgParameter<float>( "softb_pt" );
    cut_wjetpt_hard_ = getCfgParameter<float>( "hardw_pt" );
    cut_wjetpt_soft_ = getCfgParameter<float>( "softw_pt" );
    cut_tight_b_ = IDJet::tag(getCfgParameter<std::string>("tightb"));
    cut_loose_b_ = IDJet::tag(getCfgParameter<std::string>("looseb"));
    is_configured_ = true;
  }
}

bool TTPermutator::preselection(vector<IDJet*> jets, TLorentzVector* lepton, IDMet* met, int lc) {
	reset(); //clear everything
	jets_ = jets;
	lepton_ = lepton;
	met_ = met;
	lcharge_ = lc;
  //keeping only the n leading jets. 
	sort(jets_.begin(), jets_.end(), [](IDJet* A, IDJet* B){return(A->Pt() > B->Pt());});
	int reducedsize = Min(jets_.size(), cut_max_jets_);
	capped_jets_.resize(reducedsize);
	copy(jets_.begin(), jets_.begin()+reducedsize, capped_jets_.begin());
  //check b-tagging conditions
	sort(capped_jets_.begin(), capped_jets_.end(), [](IDJet* A, IDJet* B){return(A->csvIncl() > B->csvIncl());});
	if(!capped_jets_[0]->BTagId(cut_tight_b_)) return false;
	if(tracker_) tracker_->track("tight b cut");
	if(!capped_jets_[1]->BTagId(cut_loose_b_)) return false;
	if(tracker_) tracker_->track("loose b cut");

	return true;
}

Permutation TTPermutator::next(bool &keep_going) {
  //Logger::log().debug() << "<--" << jet_pos_[0] << " " << jet_pos_[1] << " " << jet_pos_[2] << " " << jet_pos_[3] << std::endl;
  keep_going = true;
	if(capped_jets_.size() == 3) capped_jets_.push_back(NULL);
  size_t bjet_cap = (bjet_idx_limit_ > 0) ? bjet_idx_limit_  : capped_jets_.size();
  for( ; jet_pos_[0] < bjet_cap ; ++jet_pos_[0]) {
    IDJet* bjet1 = capped_jets_[jet_pos_[0]];
		if(!bjet1) continue; //in 3 jet case it could be NULL, skip

    for( ; jet_pos_[1] < bjet_cap ; ++jet_pos_[1]) {
      if(jet_pos_[0] == jet_pos_[1]) continue; //same jet removal

      IDJet* bjet2 = capped_jets_[jet_pos_[1]];
			if(!bjet2) continue; //in 3 jet case it could be NULL, skip
			
      if(!(bjet1->BTagId(cut_loose_b_) && bjet2->BTagId(cut_loose_b_))) continue;
      if(!(bjet1->BTagId(cut_tight_b_) || bjet2->BTagId(cut_tight_b_))) continue;
      if(bjet1->Pt() < cut_bjetpt_hard_ && bjet2->Pt() < cut_bjetpt_hard_) continue;
      if(bjet1->Pt() < cut_bjetpt_soft_ || bjet2->Pt() < cut_bjetpt_soft_) continue;

      for( ; jet_pos_[2] < capped_jets_.size() ; ++jet_pos_[2]) {
        if(jet_pos_[0] == jet_pos_[2] || jet_pos_[1] == jet_pos_[2]) continue; //same jet removal
        IDJet* wjet1 = capped_jets_[jet_pos_[2]]; 
				if(!wjet1) continue; //in 3 jet case it could be NULL, skip

        for( ; jet_pos_[3] < capped_jets_.size() ; ++jet_pos_[3]) {
          if(jet_pos_[0] == jet_pos_[3] || jet_pos_[1] == jet_pos_[3] || jet_pos_[2] == jet_pos_[3]) continue;//same jet removal
          IDJet* wjet2 = capped_jets_[jet_pos_[3]]; 
          
					//in case of three jets it SHOULD be NULL
					if(wjet2) {
						//disambiguation W jets
						if(wjet2->Pt() > wjet1->Pt()) continue;
						if(wjet1->Pt() < cut_wjetpt_hard_ || wjet2->Pt() < cut_wjetpt_soft_) continue;
					}
					
          Permutation perm(
            wjet1, wjet2,
            bjet1, bjet2,
            lepton_, met_,
						lcharge_ //added this to pass lepton charge
            );                    
          //Logger::log().debug() <<"-->"<< jet_pos_[0] << " " << jet_pos_[1] << " " << jet_pos_[2] << " " << jet_pos_[3] << std::endl;
          //step forward the state before returning, otherwise gets stuck
          //I LOVE python's generators
          jet_pos_[3]++;
          if(jet_pos_[3] < capped_jets_.size()){
            jet_pos_[3] = 0;
            jet_pos_[2]++;
            if(jet_pos_[2] < capped_jets_.size()){
              jet_pos_[2] = 0;
              jet_pos_[1]++;
              if(jet_pos_[1] < bjet_cap){
                jet_pos_[1] = 0;
                jet_pos_[0]++; //NO RESET, if it gets beyond the limit we are done
              }
            }
          }
          return perm;
        }
        jet_pos_[3] = 0; //reset state
      }
      jet_pos_[2] = 0;
    }
    jet_pos_[1] = 0;
  }
  keep_going = false;
	if(!capped_jets_[3]) capped_jets_.pop_back();
  return Permutation();
}
