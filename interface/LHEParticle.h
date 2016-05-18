#ifndef LHEParticle_h
#define LHEParticle_h

#include "TLorentzVector.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include <vector>
#include <iostream>

class LHEParticle: public TLorentzVector {
private:
  int pdgid_;
  int status_;
  int first_mom_;
  int last_mom_;
public:
  LHEParticle(float x, float y, float z, float e, int id, int status, int fm, int lm):
    TLorentzVector(x,y,z,e),
    pdgid_(id),
    status_(status),
    first_mom_(fm),
    last_mom_(lm) {}

  int status() const {return status_;}
  int pdgId() const {return pdgid_;}
  std::pair<int, int> mothers_range() const {return std::make_pair(first_mom_, last_mom_);}
  static std::vector<LHEParticle> LHEParticles(URStreamer &event);
  friend std::ostream & operator<<(std::ostream &os, const LHEParticle &obj);
};
#endif
