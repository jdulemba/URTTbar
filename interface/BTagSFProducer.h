#ifndef BTagSFProducer_h
#define BTagSFProducer_h

#include "BTagCalibrationStandalone.h"
#include "TH2D.h"
#include "systematics.h"
#include "TFile.h"
#include <string>
#include "Logger.h"
#include "IDJet.h"
#include "DataFile.h"
#include "TFile.h"

class TTPermutator;

class BTagSFProducer {
public:
  BTagSFProducer():
    calibration_(),
    readers_tight_(),
    readers_loose_() {}
  BTagSFProducer(TTPermutator &permutator, float float_c=-1, float float_l=-1, float float_b=-1);
  BTagSFProducer(std::string tight, std::string loose, float float_c=-1, float float_l=-1, float float_b=-1);
  BTagSFProducer(const DataFile &sf_file, const DataFile &eff_file, IDJet::BTag tighttag, IDJet::BTag loosetag=IDJet::BTag::NONE, float float_c=-1, float float_l=-1, float float_b=-1);
  ~BTagSFProducer();    
  double scale_factor(const std::vector<IDJet*> &jets, systematics::SysShifts shift);
  IDJet::BTag tight_cut() {return tight_;}
  IDJet::BTag loose_cut() {return loose_;}

private:
  void configure(const DataFile &sf_file, const DataFile &eff_file, IDJet::BTag tighttag, IDJet::BTag loosetag=IDJet::BTag::NONE, float float_c=-1, float float_l=-1, float float_b=-1);

  template <class T>
  T* get_from(TFile &file, std::string path, std::string newname) {
    T* ptr = (T*) ((T*) ( file.Get( path.c_str() ) ) )->Clone(newname.c_str());
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
  float float_c_=-1;
  float float_l_=-1;
  float float_b_=-1;

  IDJet::BTag tight_ = IDJet::BTag::NONE;
  IDJet::BTag loose_ = IDJet::BTag::NONE;

  BTagCalibration calibration_;
  BTagCalibrationReader readers_tight_[3]; //[down, central, up]
  BTagCalibrationReader readers_loose_[3]; //[down, central, up]

  TH2D *eff_light_loose=0; 
  TH2D *eff_light_tight=0;
  TH2D *eff_charm_loose=0; 
  TH2D *eff_charm_tight=0;
  TH2D *eff_bottom_loose=0; 
  TH2D *eff_bottom_tight=0;
};

#endif
