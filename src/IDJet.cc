#include "IDJet.h"
#include "Logger.h"

const std::unordered_map<std::string, IDJet::BTag> IDJet::tag_names = {
  {"NONE",      IDJet::BTag::NONE},
  {"CSVLOOSE",  IDJet::BTag::CSVLOOSE},
  {"CSVMEDIUM", IDJet::BTag::CSVMEDIUM},
  {"CSVTIGHT",  IDJet::BTag::CSVTIGHT},
  {"CTAGLOOSE" , IDJet::BTag::CTAGLOOSE }, 
  {"CTAGMEDIUM", IDJet::BTag::CTAGMEDIUM}, 
  {"CTAGTIGHT" , IDJet::BTag::CTAGTIGHT }
};

IDJet::BTag IDJet::tag(const std::string label) {
  try {
    return IDJet::tag_names.at(label);
  }
  catch (std::out_of_range e){
    Logger::log().error() << "Requested IDJet ID label: " << label << " does not exist!" << std::endl;
    throw 49;
  }
}

std::string IDJet::tag2string(BTag id) {
  for(auto &entry : tag_names) {
    if(entry.second == id) return entry.first;
  }
  return "";
}

IDJet::IDType IDJet::id_type(BTag id) {
  switch(id) {
  case CSVLOOSE: 
  case CSVMEDIUM: 
  case CSVTIGHT: return IDType::CSV;
  case CTAGLOOSE: 
  case CTAGMEDIUM: 
  case CTAGTIGHT: return IDType::CTAG;
  default: return IDType::NOTSET;
  }
  return IDType::NOTSET;
}

std::string IDJet::id_string(BTag id) {
  switch(id) {
  case CSVLOOSE: 
  case CSVMEDIUM: 
  case CSVTIGHT: return "csvv2";
  case CTAGLOOSE: 
  case CTAGMEDIUM: 
  case CTAGTIGHT: return "ctag";
  default: return "";
  }
  return "";
}

BTagEntry::OperatingPoint IDJet::tag_tightness(BTag id) {
  BTagEntry::OperatingPoint val = BTagEntry::OperatingPoint::OP_NOTSET;
  switch(id) {
  case BTag::CSVLOOSE:   val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::CTAGLOOSE:  val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::CTAGMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::CSVMEDIUM:  val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::CSVTIGHT:   val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::CTAGTIGHT:  val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  default:  val = BTagEntry::OperatingPoint::OP_NOTSET; break;
  }
  return val;
}


bool IDJet::BTagId(BTag wp) const {
  double threshold = -1.;
  if(wp == BTag::NONE) return true;
  else if(wp == BTag::CSVLOOSE)  threshold = 0.460;
  else if(wp == BTag::CSVMEDIUM) threshold = 0.800;
  else if(wp == BTag::CSVTIGHT)  threshold = 0.935;
  else {
    Logger::log().fatal() << wp << "Is not a valid b-tagging working point!"<< std::endl;
    throw 42;
  }
  return csvIncl() > threshold;
}

bool IDJet::CTagId(BTag wp) const	{
  double cvsl_thr = -1.;
  double cvsb_thr = -1.;
  if(wp == BTag::NONE) return true;
  else if(wp == BTag::CTAGLOOSE)  {cvsl_thr = -0.337; cvsb_thr = -0.356;}
  else if(wp == BTag::CTAGMEDIUM) {cvsl_thr = -0.073; cvsb_thr = -0.302;}
  else if(wp == BTag::CTAGTIGHT)  {cvsl_thr = 0.294; cvsb_thr = -0.682;}
  else {
    Logger::log().fatal() << wp << "Is not a valid C-tagging working point!"<< std::endl;
    throw 42;
  }
  return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);
}

bool IDJet::TagId(BTag wp) const {
  switch(id_type(wp)) {
  case CSV: return BTagId(wp);
  case CTAG: return CTagId(wp);
  default: return false;
  }
}
