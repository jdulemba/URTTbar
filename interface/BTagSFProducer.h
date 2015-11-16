#ifndef BTagSFProducer_h
#define BTagSFProducer_h

#include "BTagCalibrationStandalone.h"
#include "TH2D.h"
#include "systematics.h"
#include "TFile.h"
#include <string>

class TTPermutator;

class BTagSFProducer {
public:
  BTagSFProducer(TTPermutator &permutator);
  ~BTagSFProducer();    
  double scale_factor(TTPermutator &permutator, systematics::SysShifts shift);

private:
  template <class T>
  T* get_from(TFile &file, std::string path, std::string newname) {
    return (T*) (dynamic_cast<T*>( file.Get( path.c_str() ) )->Clone(newname.c_str()));
  }

  bool no_loose_cut_=false;
  
  BTagCalibration calibration_;
  BTagCalibrationReader readers_tight_[3]; //[down, central, up]
  BTagCalibrationReader readers_loose_[3]; //[down, central, up]

  TH2D *eff_light_loose, *eff_light_tight;
  TH2D *eff_charm_loose, *eff_charm_tight;
  TH2D *eff_bottom_loose, *eff_bottom_tight;
};

#endif
