#ifndef LHEParticle_h
#define LHEParticle_h

#include "TLorentzVector.h"
#include "URStreamer.h"
#include <vector>
#include <iostream>

class LHEParticle: public Lhepaticle, public TLorentzVector {
private:
  int pdgid_;
  int status_;
  int first_mom_;
  int last_mom_;
public:
  LHEParticle(const Lhepaticle &lhe):
		Lhepaticle(lhe),
    TLorentzVector(px(),py(),pz(),e()) {}

  std::pair<int, int> mothers_range() const {return std::make_pair(fmother()-1, lmother()-1);}
  //static std::vector<LHEParticle> LHEParticles(URStreamer &event);
	int pdgId() const {return pdgid();} //stupid patch
  friend std::ostream & operator<<(std::ostream &os, const LHEParticle &obj);
};
#endif
