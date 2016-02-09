#ifndef BTagSFProducer_h
#define BTagSFProducer_h

#include "BTagCalibrationStandalone.h"
#include "TH2D.h"
#include "systematics.h"
#include "TFile.h"
#include <string>
#include "Logger.h"

class TTPermutator;

class BTagSFProducer {
public:
  BTagSFProducer(TTPermutator &permutator, float float_c=-1, float float_l=-1, float float_b=-1);
  ~BTagSFProducer();    
  double scale_factor(TTPermutator &permutator, systematics::SysShifts shift);

private:
  template <class T>
  T* get_from(TFile &file, std::string path, std::string newname) {
    T* ptr = (T*) ( (T*) ( file.Get( path.c_str() ) )->Clone(newname.c_str()));
    if(!ptr) {
      Logger::log().fatal() << "Could not get " << path << " from the file!" << std::endl;
      throw 42;
    }
    return ptr;
  }

  bool no_loose_cut_=false;
  //Systematics use
  //How much to float the single SF components
  //-1: how much the SF tell you to do. 0: Do not float at all (disable systematic). 0< <1 always return such value (Up > 1/Down < 1). >1 Invalid!
  float float_b_=-1;
  float float_c_=-1;
  float float_l_=-1;

  BTagCalibration calibration_;
  BTagCalibrationReader readers_tight_[3]; //[down, central, up]
  BTagCalibrationReader readers_loose_[3]; //[down, central, up]

  TH2D *eff_light_loose, *eff_light_tight;
  TH2D *eff_charm_loose, *eff_charm_tight;
  TH2D *eff_bottom_loose, *eff_bottom_tight;
};

#endif
