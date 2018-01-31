#ifndef LeptonSF_h
#define LeptonSF_h

#include "TH2.h"
#include <string>
#include <memory>
#include "TFile.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "Analyses/URTTbar/interface/systematics.h"

class LeptonSF {
public:
	typedef systematics::SysShifts Sys;
  LeptonSF(std::string parname, bool ptx=true);
  ~LeptonSF() {
		Logger::log().debug() << "LeptonSF Destructor" << std::endl;
	}
  double get_sf(double pt, double eta, Sys shift = Sys::NOSYS) const;
private:

  template <class T>
  std::shared_ptr<T> get_from(TFile &file, std::string path, std::string newname) {
    T* original = (T*) file.Get( path.c_str() );
    if(!original) {
      Logger::log().warning() << "Could not get " << path << " from the file " << file.GetName() << "!" << std::endl;
			std::shared_ptr<T> ptr;
			return ptr;
    }
    std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
    return ptr;
  }

  double get_2d_weight(std::shared_ptr<TH2> h, double pt, double eta, int shift) const;
  double get_1d_weight(std::shared_ptr<TH1> h, double pt, double eta, int shift) const;
	std::shared_ptr<TH1> trk_; 
	std::shared_ptr<TH2> id_; 
	std::shared_ptr<TH2> iso_; 
	std::shared_ptr<TH2> trig_;
  bool pt_as_x_=true;
  bool abs_etas_[4] = {true, true, true, true};
};

#endif
