#ifndef BTagSFProducer_h
#define BTagSFProducer_h

#include "URAnalysis/AnalysisFW/interface/BTagCalibrationStandalone.h"
#include "TH2D.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "TFile.h"
#include <string>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "TFile.h"
#include <memory>

class TTPermutator;

class BTagSFProducer {
public:
  BTagSFProducer():
    calibration_(),
    readers_tight_(),
    readers_loose_() {}
  BTagSFProducer(TTPermutator &permutator, float float_c=-1, float float_l=-1, float float_b=-1);
  BTagSFProducer(std::string tight, std::string loose="", float float_c=-1, float float_l=-1, float float_b=-1);
  BTagSFProducer(const DataFile &sf_file, const DataFile &eff_file, IDJet::BTag tighttag, IDJet::BTag loosetag=IDJet::BTag::NONE, float float_c=-1, float float_l=-1, float float_b=-1);
  ~BTagSFProducer();    
  double scale_factor(const std::vector<IDJet*> &jets, systematics::SysShifts shift);
  IDJet::BTag tight_cut() {return tight_;}
  IDJet::BTag loose_cut() {return loose_;}
  void ignore_partial_shifts() {ignore_partial_shifts_=true;}
  void ignore_general_shifts() {ignore_general_shifts_=true;}

private:
  void configure(const DataFile &sf_file, const DataFile &eff_file, IDJet::BTag tighttag, IDJet::BTag loosetag=IDJet::BTag::NONE, float float_c=-1, float float_l=-1, float float_b=-1);

  template <class T>
  std::shared_ptr<T> get_from(TFile &file, std::string path, std::string newname) {
    T* original = (T*) file.Get( path.c_str() );
    if(!original) {
      Logger::log().fatal() << "Could not get " << path << " from the file " << file.GetName() << "!" << std::endl;
      throw 42;
    }
    std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
    return ptr;
  }

  bool ignore_partial_shifts_=false;
  bool ignore_general_shifts_=false;
  bool no_loose_cut_=false, no_tight_cut_=false;
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

  std::shared_ptr<TH2D> eff_light_loose; 
  std::shared_ptr<TH2D> eff_light_tight;
  std::shared_ptr<TH2D> eff_charm_loose; 
  std::shared_ptr<TH2D> eff_charm_tight;
  std::shared_ptr<TH2D> eff_bottom_loose; 
  std::shared_ptr<TH2D> eff_bottom_tight;
};

#endif
