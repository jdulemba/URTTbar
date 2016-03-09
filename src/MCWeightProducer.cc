#include "MCWeightProducer.h"
#include "DataFile.h"
#include "URStreamer.h"
#include "TFile.h"
#include "TMath.h"

using namespace TMath;

float MCWeightProducer::evt_weight(URStreamer &evt, systematics::SysShifts shift) {
  const Geninfo& info = evt.genInfo();
  float ret = info.weight()/Abs(info.weight());
  
  const vector<Mcweight>& ws =  evt.MCWeights();
  //TODO: check if correct!
  switch(shift) {
  case systematics::SysShifts::FACTOR_DW: ret *= ret*ws[2].weights()/ws[0].weights(); break;
  case systematics::SysShifts::FACTOR_UP: ret *= ret*ws[1].weights()/ws[0].weights(); break;
  case systematics::SysShifts::RENORM_DW: ret *= ret*ws[6].weights()/ws[0].weights(); break;
  case systematics::SysShifts::RENORM_UP: ret *= ret*ws[3].weights()/ws[0].weights(); break;
  default: break;
  }

  //FIXME: add shifts up/down
  switch(shift) {
  case systematics::SysShifts::PU_UP: ret *= pu_sf_up_.weight(evt.PUInfos()[0].nInteractions()); break;
  case systematics::SysShifts::PU_DW: ret *= pu_sf_dw_.weight(evt.PUInfos()[0].nInteractions()); break;
  default: ret *= pu_sf_.weight(evt.PUInfos()[0].nInteractions()); break;
  }

  return ret;
}
