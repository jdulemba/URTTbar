#include "MCWeightProducer.h"
#include "DataFile.h"
#include "URStreamer.h"
#include "TFile.h"
#include "TMath.h"

using namespace TMath;

void MCWeightProducer::init(DataFile &fname) {
  TH1::AddDirectory(false);
  TFile* pu_file = TFile::Open(fname.path().c_str()); 
  pu_sf_ = (TH1D*) ((TH1D*)pu_file->Get("PUweight"))->Clone("pu_weight");
  TH1::AddDirectory(true);
}

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

  double npu = evt.vertexs().size(); //FIXME, use real PU!
  if(npu > 0)	ret *= pu_sf_->GetBinContent(pu_sf_->FindFixBin(npu));		
  
  return ret;
}
