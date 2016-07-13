#ifndef MCWeightProducer_h
#define MCWeightProducer_h

#include "Analyses/URTTbar/interface/systematics.h"
#include "URAnalysis/AnalysisFW/interface/PUReweighter.h"
#include <string>

class URStreamer;

class MCWeightProducer {
public:
  MCWeightProducer():
    pu_sf_(),
    pu_sf_up_(),
    pu_sf_dw_() {}

  void init(std::string sample) {
    std::string sfile = sample+".meta.pu.root";
    pu_sf_.init(sfile);
    pu_sf_up_.init(sfile, "data.meta.pu_up.root");
    pu_sf_dw_.init(sfile, "data.meta.pu_down.root");
  }  
	float gen_weight(URStreamer &evt, systematics::SysShifts shift = systematics::SysShifts::NOSYS);
	float pu_weight(URStreamer &evt, systematics::SysShifts shift = systematics::SysShifts::NOSYS);
  float evt_weight(URStreamer &evt, systematics::SysShifts shift = systematics::SysShifts::NOSYS) {
		return gen_weight(evt, shift)*pu_weight(evt, shift);
	}

private:
  PUReweighter pu_sf_;
  PUReweighter pu_sf_up_;
  PUReweighter pu_sf_dw_;
};

#endif
