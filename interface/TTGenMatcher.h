#ifndef TTGenMatcher_h
#define TTGenMatcher_h

#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include <string>

class GenTTBar;
class GenObject;

class TTGenMatcher {
public:
  enum MatchMode {DR_DR, DR_DPT, DR_PTMAX};
  static const map<std::string, MatchMode> name_to_mode;

  TTGenMatcher(); //configuration
  Permutation match(GenTTBar& gen_hyp, std::vector<IDJet*> &jets, std::vector<IDElectron*> &electrons, std::vector<IDMuon*> &muons);

  template <class T>
  T* best_match(const GenObject* gen, std::vector<T*> &recos){
    T* best = NULL;
    float dpt = 1e100;
    for(auto reco : recos){
      if(reco->DeltaR(*gen) > max_dr_)	continue;
 
      switch(mode_) {
      // CLOSEST IN DPT
      case TTGenMatcher::DR_DPT :
        if(fabs(reco->Pt() - gen->Pt()) < dpt){
          dpt = fabs(reco->Pt() - gen->Pt());
          best = reco;
        }
        break;
        
      //HIGHEST PT
      case TTGenMatcher::DR_PTMAX :
        if(!best || reco->Pt() > best->Pt()){
          best = reco;
        }
        break;
      
      // CLOSEST IN DR
      case TTGenMatcher::DR_DR :
        if(reco->DeltaR(*gen) < dpt){
          dpt  = reco->DeltaR(*gen);
          best = reco;
        }
        break;
      }
    }
		return best;
	}	

  template <class T>
  T* gen_match(const GenObject* gen, std::vector<T*> &recos){
		T* best = best_match(gen, recos);
    if(best) best->addMatch(gen);
    return best;
  }

private:
  float max_dr_;
  MatchMode mode_;
};

#endif
