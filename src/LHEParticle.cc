#include "LHEParticle.h"

std::vector<LHEParticle> LHEParticle::LHEParticles(URStreamer &event) {
  auto &xs  = event.PXLHEs();
  auto &ys  = event.PYLHEs();
  auto &zs  = event.PZLHEs();
  auto &es  = event.ELHEs() ;
  auto &ids = event.PDGIDLHEs();
  auto &sts = event.STATUSLHEs();
  auto &fms = event.FMOTHLHEs();
  auto &lms = event.LMOTHLHEs();
  size_t nparticles = xs.size();
  std::vector<LHEParticle> ret;
  ret.reserve(nparticles);
  for(size_t idx=0; idx<nparticles; ++idx) {
    ret.emplace_back(
      xs[idx].px(), ys[idx].py(), zs[idx].pz(), es[idx].e(), 
      ids[idx].pdgid(), sts[idx].status(), fms[idx].fmother(), lms[idx].lmother()
      );
  }
  return ret;
}
