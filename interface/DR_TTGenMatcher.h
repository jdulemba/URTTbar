#ifndef DR_TTGenMatcher_h
#define DR_TTGenMatcher_h

#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include <string>

class GenTTBar;
class GenObject;

class DR_TTGenMatcher {
private:
    float max_dr_;

public:

    DR_TTGenMatcher(); //configuration
    Permutation dr_match(GenTTBar& gen_hyp, std::vector<IDJet*> &jets, TLorentzVector *lep, IDMet *met, int lepcharge);

    template <class T>
    T* best_DR_match( const GenObject* gen, std::vector<T*> &recos){
        T* best = NULL;
        float gen_reco_DR = 1e10; // initialize DR to high number
        for( auto reco : recos ){
            if( reco->DeltaR(*gen) > max_dr_ ) continue;
            if( reco->DeltaR(*gen) < gen_reco_DR ){
                gen_reco_DR = reco->DeltaR(*gen);
                best = reco;
            } 
        }
        return best;
    }

    template <class T>
    T* gen_DR_match(const GenObject* gen, std::vector<T*> &recos){
        T* best = best_DR_match(gen, recos);
        if( best ) best->addMatch(gen);
        return best;
    }

};
#endif
