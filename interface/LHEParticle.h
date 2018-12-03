#ifndef LHEParticle_h
#define LHEParticle_h

/*
Only adds the TLorentzVector interface to the Lhepaticle class. Not provided by default by autocode as it is provided in px py py e instead of pt eta phi mass
 */

#include "TLorentzVector.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include <vector>
#include <iostream>

class LHEParticle: public Genparts, public TLorentzVector {
private:
  int pdgid_;
  int status_;
  int first_mom_;
  int last_mom_;
public:
  LHEParticle(const Genparts &lhe):
		Genparts(lhe),
        //{}
    TLorentzVector(px(),py(),pz(),e()) {}

  //std::pair<int, int> mothers_range() const {return std::make_pair(fmother()-1, lmother()-1);}
  //static std::vector<LHEParticle> LHEParticles(URStreamer &event);
	int pdgId() const {return pdgId();} //stupid patch
  friend std::ostream & operator<<(std::ostream &os, const LHEParticle &obj);
};
#endif
