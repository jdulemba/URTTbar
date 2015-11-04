#ifndef MCWeightProducer_h
#define MCWeightProducer_h

#include "TH1D.h"
#include "systematics.h"

class DataFile;
class URStreamer;

class MCWeightProducer {
public:
  MCWeightProducer():
    pu_sf_(0) {}

  ~MCWeightProducer() { 
    if(pu_sf_) delete pu_sf_;
  }

  void init(DataFile &fname);  
  float evt_weight(URStreamer &evt, systematics::SysShifts shift = systematics::SysShifts::NOSYS);

private:
  TH1D *pu_sf_;
};

#endif
