#ifndef LeptonSF_h
#define LeptonSF_h

#include "TH2.h"
#include <string>

class LeptonSF {
public:
  LeptonSF(std::string parname, bool ptx=true);
  ~LeptonSF() {
    if(id_) delete id_;
    if(iso_) delete iso_;
    if(trig_) delete trig_;
  }
  double get_sf(double pt, double eta) const;
private:
  double get_weight(TH2 *h, double pt, double eta) const;
  TH2 *id_=0, *iso_=0, *trig_=0;
  bool pt_as_x_=true;
};

#endif
