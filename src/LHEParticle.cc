#include "Analyses/URTTbar/interface/LHEParticle.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"

/*std::vector<LHEParticle> LHEParticle::LHEParticles(URStreamer &event) {
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
      ids[idx].pdgid(), sts[idx].status(), fms[idx].fmother()-1, lms[idx].lmother()-1
      );
  }
  return ret;
	}*/

std::ostream & operator<<(std::ostream &os, const LHEParticle &obj) {
  os << "LHEParticle(" << ura::pdg_names.at(obj.pdgId()) << ")";// ", " 
     //<< obj.status() << ", " << obj.mothers_range().first 
     //<< " --> " << obj.mothers_range().second << ") ";
  return os;
}
