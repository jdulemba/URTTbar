#ifndef MCWeightProducer_h
#define MCWeightProducer_h

#include "systematics.h"
#include "PUReweighter.h"
#include <string>

class URStreamer;

class MCWeightProducer {
public:
  MCWeightProducer():
    pu_sf_() {}

  void init(std::string sample) {pu_sf_.init(sample);}  
  float evt_weight(URStreamer &evt, systematics::SysShifts shift = systematics::SysShifts::NOSYS);

private:
  PUReweighter pu_sf_;
};

#endif
