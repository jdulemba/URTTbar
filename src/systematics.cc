#include "systematics.h"
#include <boost/algorithm/string/predicate.hpp>

namespace systematics {
  std::string get_sample(std::string outname) {
    size_t slash = outname.rfind("/") + 1;
    std::string basename(outname, slash, outname.size() - slash);
    size_t dot = basename.find(".");
    std::string sample(basename, 0, dot);
    if(sample.find("_out_") != std::string::npos){
      size_t out = sample.find("_out_");
      std::string real(sample, 0, out);
      sample = real;
    }
    return sample;
  }

  std::vector<SysShifts> get_systematics(std::string outname) {
    std::string sample = get_sample(outname);
    
		if(boost::starts_with(sample, "ttJets")) {
      if(sample.find("_") != std::string::npos) {
        //sys shifted sample!
        return {SysShifts::NOSYS};
      }
      else return {SysShifts::NOSYS, SysShifts::JES_UP, SysShifts::JES_DW,
          SysShifts::JER_UP, SysShifts::MET_UP, SysShifts::MET_DW};
    }
    else return {SysShifts::NOSYS};
  }
}
